//! The partial-order alignment graph and its heaviest-path consensus.
//!
//! Every graph has a virtual source and sink node, and every segment threaded into the graph
//! runs from the source to the sink. A node's weight is the number of segments through it.

use std::collections::VecDeque;

/// The id of the virtual source node, which precedes every segment.
pub(crate) const SOURCE: usize = 0;
/// The id of the virtual sink node, which follows every segment.
pub(crate) const SINK: usize = 1;

/// The base stored on the source and sink; it never matches a segment base.
const SENTINEL_BASE: u8 = 0;

/// The byte every non-`ACGT` input base is normalised to.
pub(crate) const UNKNOWN_BASE: u8 = b'N';

/// Normalises an input byte: `ACGT` in either case become upper case, anything else `N`.
pub(crate) const fn normalize_base(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'A',
        b'C' | b'c' => b'C',
        b'G' | b'g' => b'G',
        b'T' | b't' => b'T',
        _ => UNKNOWN_BASE,
    }
}

/// Returns whether two normalised bases match; `N` and the sentinel match nothing.
pub(crate) fn bases_match(a: u8, b: u8) -> bool {
    a == b && matches!(a, b'A' | b'C' | b'G' | b'T')
}

/// One step of an alignment of a segment to the graph.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum AlignedPair {
    /// The segment base at `pos` is aligned (matched or mismatched) to `node`.
    Aligned { node: usize, pos: usize },
    /// The node is missing from the segment.
    Deletion { node: usize },
    /// The segment base at `pos` is not in the graph.
    Insertion { pos: usize },
}

/// A weighted edge, stored on its target node.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Edge {
    /// The node the edge leaves.
    pub(crate) from: usize,
    /// The number of segments that traverse the edge.
    pub(crate) weight: u32,
}

#[derive(Clone, Debug)]
struct Node {
    base: u8,
    preds: Vec<Edge>,
    succs: Vec<usize>,
    /// Nodes in the same alignment column with a different base.
    aligned: Vec<usize>,
}

impl Node {
    const fn new(base: u8) -> Self {
        Self {
            base,
            preds: Vec::new(),
            succs: Vec::new(),
            aligned: Vec::new(),
        }
    }
}

/// A partial-order alignment graph (a DAG of bases) with a cached topological order.
#[derive(Debug, Default)]
pub(crate) struct PoaGraph {
    nodes: Vec<Node>,
    /// Topological order of every node except the sink, starting with the source.
    rows: Vec<usize>,
    /// Position of each node in `rows` (`usize::MAX` for the sink).
    rank: Vec<usize>,
    /// The number of segments threaded into the graph.
    segments: u32,
    /// The length of the first segment, whose nodes are ids `2..2 + backbone_len` in order.
    backbone_len: usize,
    indegree: Vec<usize>,
    queue: VecDeque<usize>,
}

impl PoaGraph {
    /// Empties the graph and threads `seq` through it as the first segment.
    pub(crate) fn reset(&mut self, seq: &[u8]) {
        self.nodes.clear();
        self.nodes.push(Node::new(SENTINEL_BASE));
        self.nodes.push(Node::new(SENTINEL_BASE));
        self.segments = 1;
        self.backbone_len = seq.len();
        let mut prev = SOURCE;
        for &base in seq {
            let node = self.add_node(base);
            self.add_edge(prev, node);
            prev = node;
        }
        self.add_edge(prev, SINK);
    }

    /// Returns the base of a node.
    pub(crate) fn base(&self, node: usize) -> u8 {
        self.nodes[node].base
    }

    /// Returns the in-edges of a node.
    pub(crate) fn preds(&self, node: usize) -> &[Edge] {
        &self.nodes[node].preds
    }

    /// Returns every node but the sink in topological order, starting with the source.
    ///
    /// Valid after a successful [`PoaGraph::topo_sort`] with no modification since.
    pub(crate) fn rows(&self) -> &[usize] {
        &self.rows
    }

    /// Returns the position of `node` in [`PoaGraph::rows`].
    pub(crate) fn rank(&self, node: usize) -> usize {
        self.rank[node]
    }

    /// Computes the topological order with Kahn's algorithm.
    ///
    /// Returns `false` if the graph has a cycle, which alignments should never create.
    pub(crate) fn topo_sort(&mut self) -> bool {
        let n = self.nodes.len();
        self.indegree.clear();
        self.indegree
            .extend(self.nodes.iter().map(|node| node.preds.len()));
        self.rank.clear();
        self.rank.resize(n, usize::MAX);
        self.rows.clear();
        self.queue.clear();
        self.queue.push_back(SOURCE);
        let mut visited = 0;
        while let Some(node) = self.queue.pop_front() {
            visited += 1;
            if node != SINK {
                self.rank[node] = self.rows.len();
                self.rows.push(node);
            }
            for &succ in &self.nodes[node].succs {
                self.indegree[succ] -= 1;
                if self.indegree[succ] == 0 {
                    self.queue.push_back(succ);
                }
            }
        }
        visited == n
    }

    /// Threads a segment into the graph along an alignment from the aligner.
    ///
    /// Aligned bases reuse the node they align to, or a node of the same base in its column;
    /// otherwise a new node joins the column. Insertions become new nodes.
    pub(crate) fn add_alignment(&mut self, seq: &[u8], alignment: &[AlignedPair]) {
        let mut prev = SOURCE;
        for &pair in alignment {
            let (target, pos) = match pair {
                AlignedPair::Aligned { node, pos } => (Some(node), pos),
                AlignedPair::Insertion { pos } => (None, pos),
                AlignedPair::Deletion { .. } => continue,
            };
            let base = seq[pos];
            let node = match target {
                Some(target) if bases_match(self.nodes[target].base, base) => target,
                Some(target) => {
                    let same_base = self.nodes[target]
                        .aligned
                        .iter()
                        .copied()
                        .find(|&other| bases_match(self.nodes[other].base, base));
                    match same_base {
                        Some(other) => other,
                        None => self.add_aligned_node(base, target),
                    }
                }
                None => self.add_node(base),
            };
            self.add_edge(prev, node);
            prev = node;
        }
        self.add_edge(prev, SINK);
        self.segments = self.segments.saturating_add(1);
    }

    /// Returns the consensus: the source-to-sink path with the highest total node score.
    ///
    /// A node of weight `w` in an alignment column of total weight `c` (its own weight plus that
    /// of the nodes aligned to it), in a graph of `n` segments, scores `w - (n - c)`: the
    /// segments with its base minus the segments that skip its column. So a base enters the
    /// consensus when more segments have it than skip the column, and of the parallel bases of
    /// a column the most common wins. Lone insertions and minority deletions lose.
    ///
    /// Ties go to the predecessor with the heavier edge, then to the edge of the backbone (the
    /// first segment), then to the lower node id. So where the votes split evenly, at a
    /// substitution or an indel alike, the consensus follows the backbone. Valid after a
    /// successful [`PoaGraph::topo_sort`].
    pub(crate) fn consensus(&self) -> Vec<u8> {
        let n_nodes = self.nodes.len();
        let weight: Vec<i64> = self
            .nodes
            .iter()
            .map(|node| node.preds.iter().map(|edge| i64::from(edge.weight)).sum())
            .collect();
        let segments = i64::from(self.segments);
        let mut best = vec![i64::MIN; n_nodes];
        let mut best_pred = vec![usize::MAX; n_nodes];
        best[SOURCE] = 0;
        let order = self.rows.iter().skip(1).copied().chain([SINK]);
        for node in order {
            let pred = self.nodes[node]
                .preds
                .iter()
                .filter(|edge| best[edge.from] != i64::MIN)
                .max_by(|a, b| {
                    best[a.from]
                        .cmp(&best[b.from])
                        .then(a.weight.cmp(&b.weight))
                        .then(
                            self.is_backbone_edge(a.from, node)
                                .cmp(&self.is_backbone_edge(b.from, node)),
                        )
                        .then(b.from.cmp(&a.from))
                })
                .map(|edge| edge.from);
            let Some(pred) = pred else { continue };
            let score = if node == SINK {
                0
            } else {
                let column: i64 = self.nodes[node]
                    .aligned
                    .iter()
                    .map(|&other| weight[other])
                    .sum::<i64>()
                    + weight[node];
                weight[node] - (segments - column)
            };
            best[node] = best[pred] + score;
            best_pred[node] = pred;
        }
        let mut out = Vec::new();
        let mut node = best_pred[SINK];
        // Each step moves to a node earlier in topological order, so this terminates.
        while node != SOURCE && node != usize::MAX {
            out.push(self.nodes[node].base);
            node = best_pred[node];
        }
        out.reverse();
        out
    }

    /// Returns whether `from -> to` belongs to the backbone chain: the source, nodes
    /// `2..2 + backbone_len` in order, then the sink.
    fn is_backbone_edge(&self, from: usize, to: usize) -> bool {
        let chain_pos = |node: usize| match node {
            SOURCE => Some(0),
            SINK => Some(self.backbone_len + 1),
            node if node < self.backbone_len + 2 => Some(node - 1),
            _ => None,
        };
        matches!((chain_pos(from), chain_pos(to)), (Some(a), Some(b)) if b == a + 1)
    }

    fn add_node(&mut self, base: u8) -> usize {
        self.nodes.push(Node::new(base));
        self.nodes.len() - 1
    }

    /// Adds a node with `base` to the alignment column of `target`.
    fn add_aligned_node(&mut self, base: u8, target: usize) -> usize {
        let node = self.add_node(base);
        for i in 0..self.nodes[target].aligned.len() {
            let other = self.nodes[target].aligned[i];
            self.nodes[other].aligned.push(node);
            self.nodes[node].aligned.push(other);
        }
        self.nodes[target].aligned.push(node);
        self.nodes[node].aligned.push(target);
        node
    }

    fn add_edge(&mut self, from: usize, to: usize) {
        let preds = &mut self.nodes[to].preds;
        if let Some(edge) = preds.iter_mut().find(|edge| edge.from == from) {
            edge.weight = edge.weight.saturating_add(1);
        } else {
            preds.push(Edge { from, weight: 1 });
            self.nodes[from].succs.push(to);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pair(node: Option<usize>, pos: Option<usize>) -> AlignedPair {
        match (node, pos) {
            (Some(node), Some(pos)) => AlignedPair::Aligned { node, pos },
            (Some(node), None) => AlignedPair::Deletion { node },
            (None, Some(pos)) => AlignedPair::Insertion { pos },
            (None, None) => panic!("a test alignment step needs a node or a position"),
        }
    }

    #[test]
    fn normalizes_bases() {
        assert_eq!(normalize_base(b'a'), b'A');
        assert_eq!(normalize_base(b'T'), b'T');
        assert_eq!(normalize_base(b'n'), UNKNOWN_BASE);
        assert_eq!(normalize_base(b'-'), UNKNOWN_BASE);
        assert!(bases_match(b'G', b'G'));
        assert!(!bases_match(b'N', b'N'));
        assert!(!bases_match(SENTINEL_BASE, SENTINEL_BASE));
    }

    #[test]
    fn majority_wins_at_a_mismatch_column() {
        // Nodes 2, 3, 4 hold A, C, G.
        let mut graph = PoaGraph::default();
        graph.reset(b"ACG");
        let alt = [
            pair(Some(2), Some(0)),
            pair(Some(3), Some(1)),
            pair(Some(4), Some(2)),
        ];
        graph.add_alignment(b"ATG", &alt);
        graph.add_alignment(b"ATG", &alt);
        assert!(graph.topo_sort());
        assert_eq!(graph.consensus(), b"ATG");
        // The second ATG reused the T node rather than adding another.
        assert_eq!(graph.nodes.len(), 6);
    }

    #[test]
    fn lone_insertions_and_deletions_lose() {
        let mut graph = PoaGraph::default();
        graph.reset(b"ACGT");
        let same = [
            pair(Some(2), Some(0)),
            pair(Some(3), Some(1)),
            pair(Some(4), Some(2)),
            pair(Some(5), Some(3)),
        ];
        graph.add_alignment(b"ACGT", &same);
        // One segment with a trailing insertion, one missing its first base.
        let insertion = [
            pair(Some(2), Some(0)),
            pair(Some(3), Some(1)),
            pair(Some(4), Some(2)),
            pair(Some(5), Some(3)),
            pair(None, Some(4)),
        ];
        graph.add_alignment(b"ACGTT", &insertion);
        let deletion = [
            pair(Some(2), None),
            pair(Some(3), Some(0)),
            pair(Some(4), Some(1)),
            pair(Some(5), Some(2)),
        ];
        graph.add_alignment(b"CGT", &deletion);
        assert!(graph.topo_sort());
        assert_eq!(graph.rows()[0], SOURCE);
        assert_eq!(graph.consensus(), b"ACGT");
    }

    #[test]
    fn even_indel_splits_follow_the_backbone() {
        // Nodes 2..=5 hold A C G T. The other segment lacks the G: the tie keeps it.
        let mut graph = PoaGraph::default();
        graph.reset(b"ACGT");
        let deletion = [
            pair(Some(2), Some(0)),
            pair(Some(3), Some(1)),
            pair(Some(4), None),
            pair(Some(5), Some(2)),
        ];
        graph.add_alignment(b"ACT", &deletion);
        assert!(graph.topo_sort());
        assert_eq!(graph.consensus(), b"ACGT");

        // Nodes 2..=4 hold A C T. The other segment has an extra G: the tie drops it.
        graph.reset(b"ACT");
        let insertion = [
            pair(Some(2), Some(0)),
            pair(Some(3), Some(1)),
            pair(None, Some(2)),
            pair(Some(4), Some(3)),
        ];
        graph.add_alignment(b"ACGT", &insertion);
        assert!(graph.topo_sort());
        assert_eq!(graph.consensus(), b"ACT");
    }

    #[test]
    fn even_substitution_splits_follow_the_backbone() {
        let mut graph = PoaGraph::default();
        graph.reset(b"ACGT");
        let alt = [
            pair(Some(2), Some(0)),
            pair(Some(3), Some(1)),
            pair(Some(4), Some(2)),
            pair(Some(5), Some(3)),
        ];
        graph.add_alignment(b"ACCT", &alt);
        assert!(graph.topo_sort());
        assert_eq!(graph.consensus(), b"ACGT");
    }

    #[test]
    fn detects_a_cycle() {
        let mut graph = PoaGraph::default();
        graph.reset(b"AC");
        assert!(graph.topo_sort());
        // Nodes 2 and 3 hold A and C; an edge back from C to A closes a cycle.
        graph.add_edge(3, 2);
        assert!(!graph.topo_sort());
    }

    #[test]
    fn split_votes_do_not_favour_insertions() {
        // Nodes 2..=8 hold T C G A A G T. Four of six segments have a single A, but they reach
        // the second A from two different predecessors, so no in-edge carries a majority.
        let mut graph = PoaGraph::default();
        graph.reset(b"TCGAAGT");
        let one_a = [
            pair(Some(2), Some(0)),
            pair(Some(3), Some(1)),
            pair(Some(4), Some(2)),
            pair(Some(5), None),
            pair(Some(6), Some(3)),
            pair(Some(7), Some(4)),
            pair(Some(8), Some(5)),
        ];
        let no_g = [
            pair(Some(2), Some(0)),
            pair(Some(3), Some(1)),
            pair(Some(4), None),
            pair(Some(5), None),
            pair(Some(6), Some(2)),
            pair(Some(7), Some(3)),
            pair(Some(8), Some(4)),
        ];
        let same: Vec<AlignedPair> = (0..7).map(|i| pair(Some(i + 2), Some(i))).collect();
        graph.add_alignment(b"TCGAGT", &one_a);
        graph.add_alignment(b"TCGAGT", &one_a);
        graph.add_alignment(b"TCAGT", &no_g);
        graph.add_alignment(b"TCAGT", &no_g);
        graph.add_alignment(b"TCGAAGT", &same);
        assert!(graph.topo_sort());
        assert_eq!(graph.consensus(), b"TCGAGT");
    }
}
