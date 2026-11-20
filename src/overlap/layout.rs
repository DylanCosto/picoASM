//! The string graph of corrected reads, and its unitigs as contigs (Myers 2005, as hifiasm
//! lays out its reads).
//!
//! Each read is two vertices, one per orientation; an overlap that runs off one read's end
//! onto the other's start (a dovetail) is an edge, with its complement between the reverse
//! orientations. The caller says how each overlap is used ([`Usage`]): after correction, those
//! joining reads of one copy with almost no difference left are the graph's edges, and those
//! between reads that differ at a site (two haplotypes, or two copies) are kept for rescue.
//! Then, in order:
//! 0. a chimeric read is dropped: one with a point inside it that no overlap (differing at no
//!    more than [`LayoutParams::max_divergence`]) runs across with
//!    [`LayoutParams::chimera_flank`] bases either side, while overlaps cover it on both sides
//!    of that point (its two parts come from two places in the genome);
//! 1. a read held whole inside another is dropped (of two reads that hold each other, the
//!    later one; a chimeric read holds none);
//! 2. an edge `v → x` is dropped where edges `v → w → x` reach `x` at the same place, give or
//!    take [`LayoutParams::fuzz`] bases (transitive reduction);
//! 3. where a read has two edges out or more, those whose overlap is shorter than
//!    [`LayoutParams::min_overlap_ratio`] of the longest are dropped (miniasm's rule: the short
//!    one is a read the longer's own overlaps missed, or another copy);
//! 4. a dead-end path of at most [`LayoutParams::max_tip_reads`] reads that runs into a vertex
//!    with another way in, one that is not such a short dead end itself, is dropped (a tip),
//!    [`LayoutParams::tip_rounds`] times; then a path that ends with no way on is joined to one
//!    that starts with no way in by the longest [`Usage::Rescue`] overlap between their end
//!    reads, longest first (a primary assembly runs on across a haplotype's sites, as hifiasm's
//!    primary contigs do); then, where the paths out of a read all meet again within
//!    [`LayoutParams::max_bubble_len`] bases with no other way in or out (a bubble: two
//!    haplotypes, two close copies, or an error the correction left), the path with the most
//!    reads is kept and the reads only on the others are dropped; then tips again;
//! 5. each maximal path through vertices with one way in and one way out is a contig, its
//!    bases each read's up to where the next starts, and the last read's in full; contigs of
//!    fewer than [`LayoutParams::min_contig_reads`] reads are left out.
//!
//! Reads and edges are visited in index order, so the contigs are the same on every run.

use crate::seq::{Strand, reverse_complement_into};

use crate::overlap::overlap::Overlap;
use crate::overlap::reads::{ReadIndex, Reads};

/// Settings for [`lay_out`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayoutParams {
    /// Two paths to a read agree on where it starts within this many bases.
    pub fuzz: u32,
    /// A dead-end path of at most this many reads is a tip.
    pub max_tip_reads: usize,
    /// Tips are clipped this many times, each round from what the last left.
    pub tip_rounds: usize,
    /// Of a read's edges out, those overlapping it by less than this share of the longest
    /// overlap are dropped.
    pub min_overlap_ratio: f64,
    /// Bubbles whose paths meet again within this many bases are popped.
    pub max_bubble_len: u64,
    /// An overlap runs across a point of a read when it covers this many bases either side.
    pub chimera_flank: u32,
    /// A contig of fewer reads is left out: a read no other continues is no assembly, and on
    /// the fly most such long reads were chimeric in ways their overlaps could not show (a
    /// hairpin, read on into its own reverse complement, has no overlap at all).
    pub min_contig_reads: usize,
    /// The overlaps that show a read is not chimeric differ at no more than this share of
    /// their columns: [`MAX_DIVERGENCE`] for corrected reads, more for raw reads whose
    /// overlaps all differ at a few percent.
    pub max_divergence: f64,
}

impl Default for LayoutParams {
    fn default() -> Self {
        Self {
            fuzz: 32,
            max_tip_reads: 3,
            tip_rounds: 3,
            min_overlap_ratio: 0.7,
            max_bubble_len: 50_000,
            chimera_flank: 1_000,
            min_contig_reads: 2,
            max_divergence: MAX_DIVERGENCE,
        }
    }
}

/// How [`lay_out`] uses an overlap.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Usage {
    /// An edge of the graph.
    Edge,
    /// An edge only between a path's dead end and another's dead start, once the graph is
    /// cleaned.
    Rescue,
    /// Left out.
    Ignore,
}

/// Overlaps differing at more than this share of their columns are not used: after correction,
/// reads of one copy differ at almost no base (two-thirds of the fly's final overlaps at none),
/// and on the fly 0.05% spanned fewer gaps (470 of 739) than 0.1% (491) or 0.2% (490).
pub const MAX_DIVERGENCE: f64 = 0.001;

/// Overlaps that are not edges and differ at no more than this share of their columns may
/// rescue a dead end. Past the fly's dead ends most reads that run on differ from the last read
/// at 0.1-1% (another haplotype or copy), but crossing to them adds misjoins about as fast as
/// gaps: at 0.1%, 0.2%, 0.3%, 0.5% and 1%, 604, 610, 615, 620 and 625 of the 739 gaps spanned,
/// with 77, 91, 104, 111 and 139 misjoins against the reference (hifiasm: 675 and 105).
pub const MAX_RESCUE_DIVERGENCE: f64 = 0.001;

/// Returns an overlap's differences as a share of its columns.
fn divergence(o: &Overlap) -> f64 {
    #[expect(
        clippy::cast_precision_loss,
        reason = "a share of an overlap's columns"
    )]
    let d = o.differences() as f64 / o.columns().max(1) as f64;
    d
}

/// Returns how [`lay_out`] uses each of `overlaps`: an edge where its reads agree at every
/// site (`consistent`, from [`crate::overlap::correct_in_rounds`]) and it differs at no more than
/// [`MAX_DIVERGENCE`] of its columns; for rescue where it differs at no more than
/// `max_rescue_divergence` (as [`MAX_RESCUE_DIVERGENCE`]); otherwise not at all.
///
/// # Panics
///
/// Panics if `consistent` is not as long as `overlaps`.
#[must_use]
pub fn usage(overlaps: &[Overlap], consistent: &[bool], max_rescue_divergence: f64) -> Vec<Usage> {
    assert_eq!(consistent.len(), overlaps.len(), "one flag per overlap");
    overlaps
        .iter()
        .zip(consistent)
        .map(|(o, &consistent)| {
            let d = divergence(o);
            if consistent && d <= MAX_DIVERGENCE {
                Usage::Edge
            } else if d <= max_rescue_divergence {
                Usage::Rescue
            } else {
                Usage::Ignore
            }
        })
        .collect()
}

/// A read in a contig.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Placed {
    /// The read.
    pub read: ReadIndex,
    /// [`Strand::Reverse`] if the contig holds it reverse-complemented.
    pub strand: Strand,
    /// Where its first base (in the contig's orientation) lies in the contig.
    pub start: u64,
}

/// A contig: a path through the string graph.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Contig {
    /// Its reads, in order.
    pub reads: Vec<Placed>,
    /// Its bases.
    pub seq: Vec<u8>,
    /// The path closes on itself: its last read runs on into its first.
    pub circular: bool,
    /// How its start and its end stop (both [`Stop::Circular`] for a circular contig).
    pub stops: [Stop; 2],
}

/// How a contig end stops.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Stop {
    /// No edge out: no read reaches on past the end.
    Dead,
    /// Two edges out or more.
    Fork,
    /// One edge out, into a read with two ways in or more.
    Join,
    /// The contig is circular.
    Circular,
}

/// What became of a read in [`lay_out`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Fate {
    /// In a contig.
    Contig,
    /// Dropped as chimeric.
    Chimeric,
    /// Held whole inside another read.
    Contained,
    /// Dropped in a tip.
    Tip,
    /// Dropped off a bubble's kept path.
    Bubble,
    /// In a contig of too few reads, left out.
    Alone,
}

/// What [`lay_out`] built, and how much each step removed.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Layout {
    /// The contigs, longest first (then by first read).
    pub contigs: Vec<Contig>,
    /// Chimeric reads dropped.
    pub chimeric: u64,
    /// Reads held whole inside another.
    pub contained: u64,
    /// Edges (counted once per pair of complements) before transitive reduction.
    pub edges: u64,
    /// Edges dropped as transitive.
    pub transitive: u64,
    /// Edges dropped as short beside another (counted once per pair of complements).
    pub short: u64,
    /// Reads dropped in tips.
    pub tip_reads: u64,
    /// Reads dropped off the kept path of a bubble.
    pub bubble_reads: u64,
    /// Dead ends joined to dead starts by [`Usage::Rescue`] overlaps.
    pub rescued: u64,
    /// Reads in contigs of fewer than [`LayoutParams::min_contig_reads`], left out.
    pub lone_reads: u64,
    /// What became of each read, by read index.
    pub fates: Vec<Fate>,
    /// How the contigs' ends (two per contig that is not circular) stop.
    pub ends: Ends,
}

/// How contig ends stop.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Ends {
    /// No edge out: no read reaches on past the end.
    pub dead: u64,
    /// Two edges out or more.
    pub forks: u64,
    /// One edge out, into a read with two ways in or more.
    pub joins: u64,
}

impl Stop {
    /// Returns how a path ending at `last` stops.
    fn at(graph: &Graph, last: Vertex) -> Self {
        match graph.out[last as usize][..] {
            [] => Self::Dead,
            [_] => Self::Join,
            _ => Self::Fork,
        }
    }
}

impl Ends {
    fn count(&mut self, stop: Stop) {
        match stop {
            Stop::Dead => self.dead += 1,
            Stop::Fork => self.forks += 1,
            Stop::Join => self.joins += 1,
            Stop::Circular => {}
        }
    }
}

/// One orientation of a read: `2 * read`, plus 1 for the reverse complement.
type Vertex = u32;

const fn vertex(read: ReadIndex, strand: Strand) -> Vertex {
    2 * read.0
        + match strand {
            Strand::Forward => 0,
            Strand::Reverse => 1,
        }
}

const fn complement(v: Vertex) -> Vertex {
    v ^ 1
}

const fn read_of(v: Vertex) -> ReadIndex {
    ReadIndex(v / 2)
}

const fn strand_of(v: Vertex) -> Strand {
    if v & 1 == 1 {
        Strand::Reverse
    } else {
        Strand::Forward
    }
}

/// An edge `v → to`: `to` starts `offset` bases into `v`, and runs on past `v`'s end.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Edge {
    offset: u32,
    to: Vertex,
}

/// The graph: each vertex's edges out, sorted by offset.
struct Graph {
    out: Vec<Vec<Edge>>,
    /// Reads still in the graph.
    live: Vec<bool>,
}

impl Graph {
    /// Returns the number of edges into `v`: those out of its complement.
    fn in_degree(&self, v: Vertex) -> usize {
        self.out[complement(v) as usize].len()
    }

    fn out_degree(&self, v: Vertex) -> usize {
        self.out[v as usize].len()
    }

    /// Returns the one vertex with an edge into `v`, if there is exactly one.
    fn only_predecessor(&self, v: Vertex) -> Option<Vertex> {
        match self.out[complement(v) as usize][..] {
            [e] => Some(complement(e.to)),
            _ => None,
        }
    }

    /// Returns the one edge out of `v`, if there is exactly one.
    fn only_edge(&self, v: Vertex) -> Option<Edge> {
        match self.out[v as usize][..] {
            [e] => Some(e),
            _ => None,
        }
    }

    /// Returns the first vertex of the path through `v`: back while the one way in comes from a
    /// vertex with one way out (around a circle, back to the vertex after `v`).
    fn path_start(&self, v: Vertex) -> Vertex {
        let mut start = v;
        while let Some(u) = self.only_predecessor(start) {
            if self.out_degree(u) != 1 || u == v {
                break;
            }
            start = u;
        }
        start
    }

    /// Returns `true` if walking back from `v` through vertices with one way in, each from a
    /// vertex with one way out, reaches a vertex with no way in within `limit` reads.
    fn short_dead_end(&self, v: Vertex, limit: usize) -> bool {
        let mut at = v;
        for _ in 0..limit {
            if self.in_degree(at) == 0 {
                return true;
            }
            match self.only_predecessor(at) {
                Some(u) if self.out_degree(u) == 1 && u != v => at = u,
                _ => return false,
            }
        }
        false
    }

    /// Drops a read and every edge into or out of either of its orientations.
    fn drop_read(&mut self, read: ReadIndex) {
        self.live[read.index()] = false;
        for v in [vertex(read, Strand::Forward), vertex(read, Strand::Reverse)] {
            let edges = std::mem::take(&mut self.out[v as usize]);
            for e in edges {
                // The complement of v → e.to is e.to' → v'.
                self.out[complement(e.to) as usize].retain(|f| f.to != complement(v));
            }
        }
    }
}

/// Where an overlap puts its two reads, both oriented as the target is read forward.
struct Geometry {
    /// The query's bases covered, in the orientation it is aligned in.
    query_start: u32,
    query_end: u32,
    query_len: u32,
    target_start: u32,
    target_end: u32,
    target_len: u32,
}

impl Geometry {
    fn new(reads: &Reads, o: &Overlap) -> Self {
        let query_len = reads.read_len(o.query);
        let (query_start, query_end) = if o.strand.is_reverse() {
            (query_len - o.query_end, query_len - o.query_start)
        } else {
            (o.query_start, o.query_end)
        };
        Self {
            query_start,
            query_end,
            query_len,
            target_start: o.target_start,
            target_end: o.target_end,
            target_len: reads.read_len(o.target),
        }
    }

    fn query_held(&self) -> bool {
        self.query_start == 0 && self.query_end == self.query_len
    }

    fn target_held(&self) -> bool {
        self.target_start == 0 && self.target_end == self.target_len
    }
}

/// Returns `true` if the stretches of a read of `len` bases that its overlaps cover (`spans`,
/// sorted here) leave a hole inside it: a point that none runs across with `flank` bases either
/// side, with such points covered both before and after it.
fn is_chimeric(spans: &mut [(u32, u32)], len: u32, flank: u32) -> bool {
    // The points each overlap runs across: its span less the flank at both ends.
    spans.sort_unstable();
    let mut reach: Option<u32> = None;
    for &(start, end) in spans.iter() {
        let (from, to) = (start.saturating_add(flank), end.saturating_sub(flank));
        if from >= to || to <= flank || from >= len.saturating_sub(flank) {
            continue;
        }
        match reach {
            Some(r) if from > r => return true,
            Some(r) => reach = Some(r.max(to)),
            None => reach = Some(to),
        }
    }
    false
}

/// Returns the edge an overlap makes, `from → to` at `offset`, and its complement's offset, if
/// it is a dovetail.
fn dovetail(reads: &Reads, o: &Overlap) -> Option<(Vertex, Vertex, u32, u32)> {
    let g = Geometry::new(reads, o);
    let (t, q) = (vertex(o.target, Strand::Forward), vertex(o.query, o.strand));
    if g.target_start > 0
        && g.target_end == g.target_len
        && g.query_start == 0
        && g.query_end < g.query_len
    {
        // The target, then the query.
        Some((t, q, g.target_start, g.query_len - g.query_end))
    } else if g.query_start > 0
        && g.query_end == g.query_len
        && g.target_start == 0
        && g.target_end < g.target_len
    {
        // The query, then the target.
        Some((q, t, g.query_start, g.target_len - g.target_end))
    } else {
        None
    }
}

/// Lays out `reads` (corrected) along the `overlaps` found between them, each used as `usage`
/// says, and returns the contigs.
///
/// # Panics
///
/// Panics if `usage` is not as long as `overlaps`, or an overlap names a read not in `reads`.
#[must_use]
pub fn lay_out(
    reads: &Reads,
    overlaps: &[Overlap],
    usage: &[Usage],
    params: &LayoutParams,
) -> Layout {
    assert_eq!(usage.len(), overlaps.len(), "one usage per overlap");
    let used = |how: Usage| {
        overlaps
            .iter()
            .zip(usage)
            .filter(move |(_, u)| **u == how)
            .map(|(o, _)| o)
    };
    let kept = || used(Usage::Edge);
    let mut layout = Layout::default();

    // 0. Chimeric reads, from every overlap used that differs at almost no base (reads of
    // another copy cover a chimera's join of two copies).
    let mut covered: Vec<Vec<(u32, u32)>> = vec![Vec::new(); reads.len()];
    for o in overlaps
        .iter()
        .zip(usage)
        .filter(|(o, u)| **u != Usage::Ignore && divergence(o) <= params.max_divergence)
        .map(|(o, _)| o)
    {
        covered[o.query.index()].push((o.query_start, o.query_end));
        covered[o.target.index()].push((o.target_start, o.target_end));
    }
    let mut live: Vec<bool> = covered
        .iter_mut()
        .zip(reads.indices())
        .map(|(spans, read)| !is_chimeric(spans, reads.read_len(read), params.chimera_flank))
        .collect();
    drop(covered);
    layout.chimeric = live.iter().filter(|&&l| !l).count() as u64;
    let chimeric: Vec<bool> = live.iter().map(|&l| !l).collect();

    // 1. Contained reads.
    for o in kept() {
        if chimeric[o.query.index()] || chimeric[o.target.index()] {
            continue;
        }
        let g = Geometry::new(reads, o);
        // Of two reads that hold each other, the target (the later read) goes.
        if g.target_held() {
            live[o.target.index()] = false;
        } else if g.query_held() {
            live[o.query.index()] = false;
        }
    }
    layout.contained = live.iter().filter(|&&l| !l).count() as u64 - layout.chimeric;
    // Reads still in the graph are settled below; the rest were contained or chimeric.
    layout.fates = chimeric
        .iter()
        .map(|&c| if c { Fate::Chimeric } else { Fate::Contained })
        .collect();

    let mut graph = Graph {
        out: vec![Vec::new(); 2 * reads.len()],
        live,
    };
    for o in kept() {
        if !graph.live[o.query.index()] || !graph.live[o.target.index()] {
            continue;
        }
        let Some((from, to, offset, back_offset)) = dovetail(reads, o) else {
            continue;
        };
        graph.out[from as usize].push(Edge { offset, to });
        graph.out[complement(to) as usize].push(Edge {
            offset: back_offset,
            to: complement(from),
        });
        layout.edges += 1;
    }
    for edges in &mut graph.out {
        edges.sort_unstable();
    }

    // 2. Transitive reduction, marked on the whole graph before anything is dropped.
    let mut transitive: Vec<(Vertex, Vertex)> = Vec::new();
    for (v, edges) in graph.out.iter().enumerate() {
        let v = Vertex::try_from(v).expect("vertices fit in u32");
        for &x in edges {
            let reached = edges.iter().any(|&w| {
                w.to != x.to
                    && w.offset < x.offset
                    && graph.out[w.to as usize].iter().any(|&wx| {
                        wx.to == x.to && (w.offset + wx.offset).abs_diff(x.offset) <= params.fuzz
                    })
            });
            if reached {
                transitive.push((v, x.to));
            }
        }
    }
    for &(v, x) in &transitive {
        graph.out[v as usize].retain(|e| e.to != x);
        graph.out[complement(x) as usize].retain(|e| e.to != complement(v));
    }
    // Each pair of complements was marked twice, once from each side.
    layout.transitive = layout.edges - graph.out.iter().map(Vec::len).sum::<usize>() as u64 / 2;

    // 3. Short edges beside long ones, marked on the whole graph before anything is dropped.
    let mut short: Vec<(Vertex, Vertex)> = Vec::new();
    for (v, edges) in graph.out.iter().enumerate() {
        if edges.len() < 2 {
            continue;
        }
        let v = Vertex::try_from(v).expect("vertices fit in u32");
        let len = reads.read_len(read_of(v));
        // Edges are sorted by offset: the first overlaps v the most.
        let longest = f64::from(len - edges[0].offset);
        for e in &edges[1..] {
            if f64::from(len - e.offset) < params.min_overlap_ratio * longest {
                short.push((v, e.to));
            }
        }
    }
    let before = graph.out.iter().map(Vec::len).sum::<usize>();
    for &(v, x) in &short {
        graph.out[v as usize].retain(|e| e.to != x);
        graph.out[complement(x) as usize].retain(|e| e.to != complement(v));
    }
    layout.short = (before - graph.out.iter().map(Vec::len).sum::<usize>()) as u64 / 2;

    // 4. Tips, bubbles, and the tips the bubbles leave.
    let before = graph.live.clone();
    layout.tip_reads = clip_tips(&mut graph, params);
    let after_tips = graph.live.clone();
    // 5. Rescue, before bubbles (on the fly 604 gaps spanned against 600 after them).
    let rescue: Vec<&Overlap> = used(Usage::Rescue).collect();
    layout.rescued = rescue_dead_ends(&mut graph, reads, &rescue);
    layout.bubble_reads = pop_bubbles(&mut graph, params.max_bubble_len);
    let after_bubbles = graph.live.clone();
    layout.tip_reads += clip_tips(&mut graph, params);
    for (i, fate) in layout.fates.iter_mut().enumerate() {
        if before[i] && !graph.live[i] {
            *fate = if after_tips[i] && !after_bubbles[i] {
                Fate::Bubble
            } else {
                Fate::Tip
            };
        }
    }

    // 6. Unitigs.
    let mut seen = vec![false; graph.out.len()];
    let mut bases = Vec::new();
    for v in 0..Vertex::try_from(graph.out.len()).expect("vertices fit in u32") {
        if seen[v as usize] || !graph.live[read_of(v).index()] {
            continue;
        }
        let start = graph.path_start(v);
        let mut contig = Contig {
            reads: Vec::new(),
            seq: Vec::new(),
            circular: false,
            stops: [Stop::Circular; 2],
        };
        let mut at = start;
        loop {
            seen[at as usize] = true;
            seen[complement(at) as usize] = true;
            let read = read_of(at);
            let strand = strand_of(at);
            let bases_of = |bases: &mut Vec<u8>| match strand {
                Strand::Forward => {
                    bases.clear();
                    bases.extend_from_slice(reads.get(read));
                }
                Strand::Reverse => reverse_complement_into(reads.get(read), bases),
            };
            bases_of(&mut bases);
            contig.reads.push(Placed {
                read,
                strand,
                start: contig.seq.len() as u64,
            });
            let next = graph
                .only_edge(at)
                .filter(|e| graph.in_degree(e.to) == 1 && !seen[e.to as usize]);
            match next {
                Some(e) => {
                    contig.seq.extend_from_slice(&bases[..e.offset as usize]);
                    at = e.to;
                }
                None => {
                    match graph.only_edge(at) {
                        // The last read runs on into the first: its bases up to there close
                        // the circle.
                        Some(e) if e.to == start && graph.in_degree(start) == 1 => {
                            contig.circular = true;
                            contig.seq.extend_from_slice(&bases[..e.offset as usize]);
                        }
                        _ => contig.seq.extend_from_slice(&bases),
                    }
                    break;
                }
            }
        }
        if !contig.circular {
            contig.stops = [Stop::at(&graph, complement(start)), Stop::at(&graph, at)];
        }
        let fate = if contig.reads.len() < params.min_contig_reads {
            Fate::Alone
        } else {
            Fate::Contig
        };
        for placed in &contig.reads {
            layout.fates[placed.read.index()] = fate;
        }
        if fate == Fate::Alone {
            layout.lone_reads += contig.reads.len() as u64;
            continue;
        }
        for stop in contig.stops {
            layout.ends.count(stop);
        }
        layout.contigs.push(contig);
    }
    layout.contigs.sort_by(|a, b| {
        b.seq
            .len()
            .cmp(&a.seq.len())
            .then(a.reads[0].read.cmp(&b.reads[0].read))
    });
    layout
}

/// Joins dead ends to dead starts by `rescue` overlaps, the longest first, and returns how many
/// were joined; never onto the start of the dead end's own path (a circle of a tandem array's
/// reads). Joining a dead end into any read instead, a bubble for popping to settle, spanned no
/// more gaps on the fly (603 against 604), and at 1% divergence fewer (589).
fn rescue_dead_ends(graph: &mut Graph, reads: &Reads, rescue: &[&Overlap]) -> u64 {
    let mut joins: Vec<(u32, Vertex, Vertex, u32, u32)> = rescue
        .iter()
        .filter(|o| graph.live[o.query.index()] && graph.live[o.target.index()])
        .filter_map(|o| dovetail(reads, o))
        .map(|(from, to, offset, back_offset)| {
            let overlap = reads.read_len(read_of(from)) - offset;
            (overlap, from, to, offset, back_offset)
        })
        .collect();
    joins.sort_unstable_by(|a, b| b.0.cmp(&a.0).then((a.1, a.2).cmp(&(b.1, b.2))));
    let mut joined = 0;
    for (_, from, to, offset, back_offset) in joins {
        if read_of(from) != read_of(to)
            && graph.out_degree(from) == 0
            && graph.in_degree(to) == 0
            && graph.path_start(from) != to
        {
            graph.out[from as usize].push(Edge { offset, to });
            let back = &mut graph.out[complement(to) as usize];
            back.push(Edge {
                offset: back_offset,
                to: complement(from),
            });
            back.sort_unstable();
            joined += 1;
        }
    }
    joined
}

/// Clips tips (see the module docs) and returns the reads dropped.
fn clip_tips(graph: &mut Graph, params: &LayoutParams) -> u64 {
    let mut dropped = 0;
    for _ in 0..params.tip_rounds {
        let mut tips: Vec<ReadIndex> = Vec::new();
        for v in 0..Vertex::try_from(graph.out.len()).expect("vertices fit in u32") {
            if !graph.live[read_of(v).index()] || graph.in_degree(v) != 0 {
                continue;
            }
            // Walk out of the dead end while the path does not branch.
            let mut path = vec![read_of(v)];
            let mut at = v;
            let is_tip = loop {
                let Some(e) = graph.only_edge(at) else {
                    break false;
                };
                if graph.in_degree(e.to) > 1 {
                    // Another way into the vertex the path runs into, not a short dead end.
                    let others = &graph.out[complement(e.to) as usize];
                    break others.iter().any(|o| {
                        let p = complement(o.to);
                        p != at && !graph.short_dead_end(p, params.max_tip_reads)
                    });
                }
                if path.len() >= params.max_tip_reads {
                    break false;
                }
                path.push(read_of(e.to));
                at = e.to;
            };
            if is_tip {
                tips.extend(path);
            }
        }
        tips.sort_unstable();
        tips.dedup();
        if tips.is_empty() {
            break;
        }
        dropped += tips.len() as u64;
        for read in tips {
            graph.drop_read(read);
        }
    }
    dropped
}

/// Pops bubbles (see the module docs) and returns the reads dropped.
fn pop_bubbles(graph: &mut Graph, max_len: u64) -> u64 {
    let mut dropped = 0;
    let mut search = BubbleSearch::default();
    for v in 0..Vertex::try_from(graph.out.len()).expect("vertices fit in u32") {
        if !graph.live[read_of(v).index()] || graph.out_degree(v) < 2 {
            continue;
        }
        let Some(drop) = search.find(graph, v, max_len) else {
            continue;
        };
        dropped += drop.len() as u64;
        for read in drop {
            graph.drop_read(read);
        }
    }
    dropped
}

/// Buffers for finding one bubble, reused from source to source.
#[derive(Default)]
struct BubbleSearch {
    /// For each vertex reached: its distance from the source in bases (the longest way), its
    /// in-edges not yet walked, the reads on its best path from the source and the vertex that
    /// path comes from.
    reached: std::collections::HashMap<Vertex, Reached>,
    ready: Vec<Vertex>,
}

#[derive(Clone, Copy)]
struct Reached {
    distance: u64,
    waiting: usize,
    reads: usize,
    from: Vertex,
}

impl BubbleSearch {
    /// Returns the reads to drop if `source` opens a bubble that closes within `max_len`
    /// bases: every path out of it meets again at one vertex (the sink), with no way into the
    /// paths from outside and no way out but through the sink. The path with the most reads is
    /// kept; the reads only on the others are returned.
    fn find(&mut self, graph: &Graph, source: Vertex, max_len: u64) -> Option<Vec<ReadIndex>> {
        self.reached.clear();
        self.ready.clear();
        self.reached.insert(
            source,
            Reached {
                distance: 0,
                waiting: 0,
                reads: 1,
                from: source,
            },
        );
        self.ready.push(source);
        // Vertices reached whose in-edges are not all walked yet.
        let mut pending = 0_usize;
        let sink = loop {
            let at = self.ready.pop()?;
            if at != source && graph.out_degree(at) == 0 {
                // A dead end inside: a tip, not a bubble.
                return None;
            }
            let here = self.reached[&at];
            for e in &graph.out[at as usize] {
                if e.to == source || e.to == complement(source) {
                    return None;
                }
                let distance = here.distance + u64::from(e.offset);
                if distance > max_len {
                    return None;
                }
                let next = self.reached.entry(e.to).or_insert_with(|| {
                    pending += 1;
                    Reached {
                        distance,
                        waiting: graph.in_degree(e.to),
                        reads: 0,
                        from: at,
                    }
                });
                next.distance = next.distance.max(distance);
                if here.reads + 1 > next.reads {
                    next.reads = here.reads + 1;
                    next.from = at;
                }
                next.waiting -= 1;
                if next.waiting == 0 {
                    pending -= 1;
                    self.ready.push(e.to);
                }
            }
            if self.ready.len() == 1 && pending == 0 {
                let last = self.ready[0];
                // The sink: every path met here.
                if last == source {
                    return None;
                }
                break last;
            }
            if self.ready.is_empty() {
                return None;
            }
        };
        // The kept path, back from the sink.
        let mut kept = std::collections::HashSet::new();
        let mut at = sink;
        while at != source {
            kept.insert(read_of(at));
            at = self.reached[&at].from;
        }
        kept.insert(read_of(source));
        let mut drop: Vec<ReadIndex> = self
            .reached
            .keys()
            .map(|&v| read_of(v))
            .filter(|r| !kept.contains(r))
            .collect();
        drop.sort_unstable();
        drop.dedup();
        // A read in the bubble both ways round (an inversion) is left alone.
        let both = self
            .reached
            .keys()
            .any(|&v| self.reached.contains_key(&complement(v)));
        if both || drop.is_empty() {
            return None;
        }
        Some(drop)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::overlap::align::{Edit, EditOp};

    /// The default settings, with flanks short enough for the tests' sparse tilings (reads 2 kb
    /// apart, 5 kb long).
    fn params() -> LayoutParams {
        LayoutParams {
            chimera_flank: 200,
            min_contig_reads: 1,
            ..LayoutParams::default()
        }
    }

    /// A random genome from a fixed linear congruential sequence.
    fn genome(len: usize, seed: u64) -> Vec<u8> {
        let mut state = seed;
        (0..len)
            .map(|_| {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                b"ACGT"[usize::try_from(state >> 62).expect("two bits")]
            })
            .collect()
    }

    /// The exact overlap of reads `q` and `t` (q < t), placed at `q_at` and `t_at` in the
    /// genome, the query reverse-complemented when `reverse`.
    fn overlap(
        reads: &Reads,
        (q, q_at): (u32, u32),
        (t, t_at): (u32, u32),
        strand: Strand,
    ) -> Overlap {
        let (ql, tl) = (reads.read_len(ReadIndex(q)), reads.read_len(ReadIndex(t)));
        let from = q_at.max(t_at);
        let to = (q_at + ql).min(t_at + tl);
        let (mut qs, mut qe) = (from - q_at, to - q_at);
        if strand.is_reverse() {
            (qs, qe) = (ql - qe, ql - qs);
        }
        Overlap {
            query: ReadIndex(q),
            target: ReadIndex(t),
            strand,
            query_start: qs,
            query_end: qe,
            target_start: from - t_at,
            target_end: to - t_at,
            edits: vec![Edit::new(EditOp::Match, to - from)],
        }
    }

    /// Reads tiling `g` at `starts`, each `len` long, every other one reverse-complemented,
    /// and the overlaps between every pair that shares 500 bases or more.
    fn tiling(g: &[u8], starts: &[u32], len: u32) -> (Reads, Vec<Overlap>) {
        let strand = |i: usize| {
            if i % 2 == 1 {
                Strand::Reverse
            } else {
                Strand::Forward
            }
        };
        let seqs: Vec<Vec<u8>> = starts
            .iter()
            .enumerate()
            .map(|(i, &s)| {
                let fwd = &g[s as usize..(s + len) as usize];
                if strand(i).is_reverse() {
                    crate::seq::reverse_complement(fwd)
                } else {
                    fwd.to_vec()
                }
            })
            .collect();
        let reads: Reads = seqs.iter().map(Vec::as_slice).collect();
        let mut overlaps = Vec::new();
        for i in 0..starts.len() {
            for j in i + 1..starts.len() {
                let (a, b) = (starts[i], starts[j]);
                if a.max(b) + 500 <= (a + len).min(b + len) {
                    // Aligned to the target as it is read: the query is reversed relative to it
                    // when their strands differ.
                    let rel = if strand(i) == strand(j) {
                        Strand::Forward
                    } else {
                        Strand::Reverse
                    };
                    let mut o = overlap(
                        &reads,
                        (u32::try_from(i).expect("small"), a),
                        (u32::try_from(j).expect("small"), b),
                        Strand::Forward,
                    );
                    // Express the overlap in each read's own orientation.
                    let flip = |s: u32, e: u32, l: u32, reverse: bool| {
                        if reverse { (l - e, l - s) } else { (s, e) }
                    };
                    (o.query_start, o.query_end) =
                        flip(o.query_start, o.query_end, len, strand(i).is_reverse());
                    (o.target_start, o.target_end) =
                        flip(o.target_start, o.target_end, len, strand(j).is_reverse());
                    o.strand = rel;
                    overlaps.push(o);
                }
            }
        }
        (reads, overlaps)
    }

    #[test]
    fn a_tiling_lays_out_as_one_contig_spelling_the_genome() {
        let g = genome(12_000, 5);
        let starts = [0, 1_500, 3_000, 4_000, 5_500, 7_000];
        let (reads, overlaps) = tiling(&g, &starts, 5_000);
        let keep = vec![Usage::Edge; overlaps.len()];
        let layout = lay_out(&reads, &overlaps, &keep, &params());
        assert_eq!(layout.contained, 0);
        assert!(layout.transitive > 0);
        assert_eq!(layout.contigs.len(), 1);
        let contig = &layout.contigs[0];
        assert_eq!(contig.reads.len(), starts.len());
        let spelled = &contig.seq;
        let want = &g[..12_000];
        assert!(
            spelled == want || *spelled == crate::seq::reverse_complement(want),
            "the contig spells the genome or its reverse complement"
        );
    }

    #[test]
    fn contained_reads_and_dropped_overlaps_leave_the_layout() {
        let g = genome(9_000, 6);
        // Read 1 lies inside read 0, which overlaps read 2.
        let (reads, overlaps) = {
            let (r, mut o) = tiling(&g, &[0, 4_000], 5_000);
            let inner = &g[1_000..4_500];
            let mut seqs: Vec<&[u8]> = vec![r.get(ReadIndex(0)), inner, r.get(ReadIndex(1))];
            let reads: Reads = seqs.drain(..).collect();
            // Renumber: the old read 1 is now read 2.
            for x in &mut o {
                x.target = ReadIndex(2);
            }
            let mut held = overlap(&reads, (0, 0), (1, 1_000), Strand::Forward);
            held.strand = Strand::Forward;
            o.insert(0, held);
            (reads, o)
        };
        let keep = vec![Usage::Edge; overlaps.len()];
        let layout = lay_out(&reads, &overlaps, &keep, &params());
        assert_eq!(layout.contained, 1);
        assert_eq!(layout.contigs.len(), 1);
        assert_eq!(layout.contigs[0].seq.len(), 9_000);

        // Without its one overlap, the pair falls apart into two contigs.
        let keep = [Usage::Edge, Usage::Ignore];
        let layout = lay_out(&reads, &overlaps, &keep, &params());
        assert_eq!(layout.contigs.len(), 2);
        // Which are a read each, left out where a contig needs two.
        let two = LayoutParams {
            min_contig_reads: 2,
            ..params()
        };
        let layout = lay_out(&reads, &overlaps, &keep, &two);
        assert_eq!((layout.contigs.len(), layout.lone_reads), (0, 2));
    }

    #[test]
    fn a_bubble_keeps_one_path() {
        let g = genome(9_000, 9);
        // Reads 1 and 2 hold the same stretch, a base apart at 4,500, and do not overlap each
        // other: two paths from read 0 to read 3.
        let (mut reads, overlaps) = tiling(&g, &[0, 2_000, 2_100, 4_000], 5_000);
        let mut two = reads.get(ReadIndex(2)).to_vec();
        two[2_400] = if two[2_400] == b'A' { b'C' } else { b'A' };
        reads = (0..4)
            .map(|i| {
                if i == 2 {
                    two.clone()
                } else {
                    reads.get(ReadIndex(i)).to_vec()
                }
            })
            .collect::<Vec<_>>()
            .iter()
            .map(Vec::as_slice)
            .collect();
        let keep: Vec<Usage> = overlaps
            .iter()
            .map(|o| {
                if (o.query.0, o.target.0) == (1, 2) {
                    Usage::Ignore
                } else {
                    Usage::Edge
                }
            })
            .collect();
        let layout = lay_out(&reads, &overlaps, &keep, &params());
        assert_eq!(layout.bubble_reads, 1);
        assert_eq!(layout.contigs.len(), 1);
        assert_eq!(layout.contigs[0].reads.len(), 3);
        assert_eq!(layout.contigs[0].seq.len(), 9_000);

        let unpopped = LayoutParams {
            max_bubble_len: 1_000,
            ..params()
        };
        let layout = lay_out(&reads, &overlaps, &keep, &unpopped);
        assert_eq!(layout.bubble_reads, 0);
        assert_eq!(layout.contigs.len(), 4);
    }

    #[test]
    fn a_rescue_overlap_joins_a_dead_end_to_a_dead_start() {
        let g = genome(11_000, 10);
        let (reads, overlaps) = tiling(&g, &[0, 2_000, 4_000, 6_000], 5_000);
        // Reads 0 and 1, and 2 and 3, are joined by edges; the overlaps between the pairs differ
        // at a site.
        let between = |o: &Overlap| (o.query.0 < 2) != (o.target.0 < 2);
        let with = |how: Usage| -> Vec<Usage> {
            overlaps
                .iter()
                .map(|o| if between(o) { how } else { Usage::Edge })
                .collect()
        };
        let layout = lay_out(&reads, &overlaps, &with(Usage::Ignore), &params());
        assert_eq!((layout.rescued, layout.contigs.len()), (0, 2));
        let layout = lay_out(&reads, &overlaps, &with(Usage::Rescue), &params());
        assert_eq!((layout.rescued, layout.contigs.len()), (1, 1));
        assert_eq!(layout.contigs[0].seq.len(), 11_000);
    }

    #[test]
    fn a_read_whose_overlaps_leave_a_hole_inside_it_is_chimeric() {
        // A 10 kb read: overlaps over its first 6 kb and its last 6 kb run across 1-5 kb and
        // 5-9 kb, and meet.
        assert!(!is_chimeric(
            &mut [(0, 6_000), (4_000, 10_000)],
            10_000,
            1_000
        ));
        // Over its first 5 kb and its last 5 kb: nothing runs across 4-6 kb.
        assert!(is_chimeric(
            &mut [(5_000, 10_000), (0, 5_000)],
            10_000,
            1_000
        ));
        // Only one side covered: a read running out of the gathered reads, not a chimera.
        assert!(!is_chimeric(&mut [(0, 5_000), (0, 4_000)], 10_000, 1_000));
        // Overlaps too short to run across anything count for nothing.
        assert!(!is_chimeric(
            &mut [(0, 5_000), (6_000, 7_500)],
            10_000,
            1_000
        ));
        assert!(!is_chimeric(&mut [], 10_000, 1_000));
    }

    #[test]
    fn a_chimeric_read_is_dropped_and_holds_nothing() {
        let g = genome(9_000, 12);
        let other = genome(9_000, 13);
        let (reads, mut overlaps) = tiling(&g, &[0, 4_000], 5_000);
        // Read 2: the first read's 2-5 kb, then 3 kb from elsewhere; read 3 is those 3 kb with
        // more of elsewhere after, and read 4 lies inside read 2's first part.
        let mut chimera = g[2_000..5_000].to_vec();
        chimera.extend_from_slice(&other[0..3_000]);
        let seqs: Vec<&[u8]> = vec![
            reads.get(ReadIndex(0)),
            reads.get(ReadIndex(1)),
            &chimera,
            &other[0..5_000],
            &g[2_000..4_000],
        ];
        let reads: Reads = seqs.into_iter().collect();
        let exact = |query, target, (qs, qe), (ts, te)| Overlap {
            query: ReadIndex(query),
            target: ReadIndex(target),
            strand: Strand::Forward,
            query_start: qs,
            query_end: qe,
            target_start: ts,
            target_end: te,
            edits: vec![Edit::new(EditOp::Match, qe - qs)],
        };
        overlaps.push(exact(0, 2, (2_000, 5_000), (0, 3_000)));
        overlaps.push(exact(2, 3, (3_000, 6_000), (0, 3_000)));
        overlaps.push(exact(2, 4, (0, 2_000), (0, 2_000)));
        overlaps.push(exact(0, 4, (2_000, 4_000), (0, 2_000)));
        let keep = vec![Usage::Edge; overlaps.len()];
        let params = LayoutParams {
            chimera_flank: 500,
            ..params()
        };
        let layout = lay_out(&reads, &overlaps, &keep, &params);
        assert_eq!(layout.chimeric, 1);
        // Read 4 is held by read 0, not only by the chimera.
        assert_eq!(layout.contained, 1);
        let lengths: Vec<usize> = layout.contigs.iter().map(|c| c.seq.len()).collect();
        assert_eq!(lengths, [9_000, 5_000]);
    }

    #[test]
    fn a_circle_is_spelled_once_round() {
        // Three reads around a 6 kb circle, each starting 2 kb after the last.
        let g = genome(6_000, 11);
        let around: Vec<u8> = g.iter().chain(&g).copied().collect();
        let (reads, overlaps) = tiling(&around, &[0, 2_000, 4_000], 5_000);
        // The third read's last 3 kb are the first's first 3 kb: an overlap past the end.
        let mut overlaps = overlaps;
        overlaps.push(Overlap {
            query: ReadIndex(0),
            target: ReadIndex(2),
            strand: Strand::Forward,
            query_start: 0,
            query_end: 3_000,
            target_start: 2_000,
            target_end: 5_000,
            edits: vec![Edit::new(EditOp::Match, 3_000)],
        });
        let keep = vec![Usage::Edge; overlaps.len()];
        let layout = lay_out(&reads, &overlaps, &keep, &params());
        assert_eq!(layout.contigs.len(), 1);
        let contig = &layout.contigs[0];
        assert!(contig.circular);
        assert_eq!(contig.seq.len(), 6_000);
        assert_eq!(contig.stops, [Stop::Circular; 2]);
    }

    #[test]
    fn rescue_does_not_close_a_path_onto_itself() {
        // The same circle, its closing overlap differing at a site: it stays a line.
        let g = genome(6_000, 11);
        let around: Vec<u8> = g.iter().chain(&g).copied().collect();
        let (reads, mut overlaps) = tiling(&around, &[0, 2_000, 4_000], 5_000);
        overlaps.push(Overlap {
            query: ReadIndex(0),
            target: ReadIndex(2),
            strand: Strand::Forward,
            query_start: 0,
            query_end: 3_000,
            target_start: 2_000,
            target_end: 5_000,
            edits: vec![Edit::new(EditOp::Match, 3_000)],
        });
        let mut usage = vec![Usage::Edge; overlaps.len()];
        *usage.last_mut().expect("one overlap") = Usage::Rescue;
        let layout = lay_out(&reads, &overlaps, &usage, &params());
        assert_eq!((layout.rescued, layout.contigs.len()), (0, 1));
        assert!(!layout.contigs[0].circular);
        assert_eq!(layout.contigs[0].seq.len(), 9_000);
    }

    #[test]
    fn a_short_dead_end_off_the_path_is_clipped() {
        let g = genome(20_000, 7);
        let other = genome(3_000, 8);
        // Reads 0-7 tile the genome; read 8 shares read 1's end and then runs elsewhere.
        let starts: Vec<u32> = (0..8).map(|i| i * 2_000).collect();
        let (r, mut overlaps) = tiling(&g, &starts, 5_000);
        let mut tip = g[5_000..7_000].to_vec();
        tip.extend_from_slice(&other);
        let mut seqs: Vec<&[u8]> = (0..8).map(|i| r.get(ReadIndex(i))).collect();
        seqs.push(&tip);
        let reads: Reads = seqs.into_iter().collect();
        // Read 1 (reverse-complemented in the tiling, at 2,000-7,000) holds the tip's first
        // 2 kb at its own 3,000-5,000 forward, which reverse-complemented is 0-2,000.
        overlaps.push(Overlap {
            query: ReadIndex(1),
            target: ReadIndex(8),
            strand: Strand::Reverse,
            query_start: 0,
            query_end: 2_000,
            target_start: 0,
            target_end: 2_000,
            edits: vec![Edit::new(EditOp::Match, 2_000)],
        });
        let keep = vec![Usage::Edge; overlaps.len()];
        let no_ratio = LayoutParams {
            min_overlap_ratio: 0.0,
            ..params()
        };
        let layout = lay_out(&reads, &overlaps, &keep, &no_ratio);
        assert_eq!((layout.short, layout.tip_reads), (0, 1));
        assert_eq!(layout.contigs.len(), 1);
        assert_eq!(layout.contigs[0].seq.len(), 19_000);

        // The tip overlaps read 1 by 2 kb beside read 2's 3 kb: under 70%, its edge goes first,
        // and it stands alone.
        let layout = lay_out(&reads, &overlaps, &keep, &params());
        assert_eq!((layout.short, layout.tip_reads), (1, 0));
        assert_eq!(layout.contigs.len(), 2);
        assert_eq!(layout.contigs[0].seq.len(), 19_000);
    }
}
