//! Pairwise alignment of a query to a target, with the partial-order aligner.
//!
//! A target is a partial-order graph of one segment, so aligning a query to it is the same
//! banded, affine-gap global alignment the consensus uses (see the [crate documentation](crate)),
//! with the same band widening, cell cap and tie-breaking.

use crate::align::Aligner;
use crate::graph::{AlignedPair, PoaGraph, bases_match, normalize_base};
use crate::{InvalidPoaParamsError, PoaParams};

/// The first node id of a target's bases in its graph: ids 0 and 1 are the source and sink.
const FIRST_BASE: usize = 2;

/// One column of a pairwise alignment, in target order.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum PairOp {
    /// A target base and a query base that match.
    Match,
    /// A target base and a query base that differ (or either is not `ACGT`).
    Mismatch,
    /// A query base the target lacks.
    Insertion,
    /// A target base the query lacks.
    Deletion,
}

/// Reusable buffers for global pairwise alignment.
///
/// # Examples
///
/// ```
/// use picoasm_consensus::{PairAligner, PairOp, PoaParams};
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let mut aligner = PairAligner::new(PoaParams::default())?;
/// let mut ops = Vec::new();
/// assert!(aligner.align(b"ACGTTGCA", b"ACGATGCA", &mut ops));
/// assert_eq!(ops[3], PairOp::Mismatch);
/// assert_eq!(ops.len(), 8);
/// # Ok(())
/// # }
/// ```
#[derive(Debug)]
pub struct PairAligner {
    params: PoaParams,
    graph: PoaGraph,
    aligner: Aligner,
    target: Vec<u8>,
    query: Vec<u8>,
}

impl PairAligner {
    /// Creates an aligner with validated parameters.
    ///
    /// # Errors
    ///
    /// Returns the error from [`PoaParams::validate`] if the parameters are unusable.
    pub fn new(params: PoaParams) -> Result<Self, InvalidPoaParamsError> {
        params.validate()?;
        Ok(Self {
            params,
            graph: PoaGraph::default(),
            aligner: Aligner::default(),
            target: Vec::new(),
            query: Vec::new(),
        })
    }

    /// Aligns all of `query` to all of `target` and writes the columns, in target order, to
    /// `ops` (cleared first). Bases are compared as the consensus compares them: `ACGT` in
    /// either case, anything else mismatching everything.
    ///
    /// Returns `false`, leaving `ops` empty, if the alignment does not fit the aligner's cell
    /// cap (48 MiB of matrices) even with its narrowest band.
    pub fn align(&mut self, target: &[u8], query: &[u8], ops: &mut Vec<PairOp>) -> bool {
        ops.clear();
        self.target.clear();
        self.target
            .extend(target.iter().map(|&b| normalize_base(b)));
        self.query.clear();
        self.query.extend(query.iter().map(|&b| normalize_base(b)));
        self.graph.reset(&self.target);
        if !self.graph.topo_sort() {
            return false;
        }
        let Some(path) = self.aligner.align(&self.graph, &self.query, &self.params) else {
            return false;
        };
        ops.extend(path.iter().map(|&pair| match pair {
            AlignedPair::Aligned { node, pos } => {
                if bases_match(self.target[node - FIRST_BASE], self.query[pos]) {
                    PairOp::Match
                } else {
                    PairOp::Mismatch
                }
            }
            AlignedPair::Deletion { .. } => PairOp::Deletion,
            AlignedPair::Insertion { .. } => PairOp::Insertion,
        }));
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn align(target: &[u8], query: &[u8]) -> Vec<PairOp> {
        let mut aligner = PairAligner::new(PoaParams::default()).expect("default parameters");
        let mut ops = Vec::new();
        assert!(aligner.align(target, query, &mut ops));
        ops
    }

    /// Counts each kind of column: matches, mismatches, insertions, deletions.
    fn counts(ops: &[PairOp]) -> [usize; 4] {
        let mut n = [0; 4];
        for op in ops {
            n[*op as usize] += 1;
        }
        n
    }

    #[test]
    fn identical_sequences_align_base_to_base() {
        let seq = b"ACGTTGCATGCAAGTCAGT";
        assert_eq!(align(seq, seq), vec![PairOp::Match; seq.len()]);
        // Lower case compares as upper case.
        assert_eq!(
            align(seq, b"acgttgcatgcaagtcagt"),
            vec![PairOp::Match; seq.len()]
        );
    }

    #[test]
    fn a_substitution_is_one_mismatch() {
        let ops = align(b"ACGTTGCATGCAAGTC", b"ACGTTGCTTGCAAGTC");
        assert_eq!(counts(&ops), [15, 1, 0, 0]);
        assert_eq!(ops[7], PairOp::Mismatch);
    }

    #[test]
    fn short_indels_are_insertions_and_deletions() {
        // Two bases more in the query, then one base fewer.
        let target = b"ACGTTGCATGCAAGTCAGTCCATG";
        let query = b"ACGTTGCAGGTGCAAGTCAGCCATG";
        let ops = align(target, query);
        let [matched, mismatched, inserted, deleted] = counts(&ops);
        assert_eq!((mismatched, inserted, deleted), (0, 2, 1));
        assert_eq!(matched + deleted, target.len());
        assert_eq!(matched + inserted, query.len());
    }

    #[test]
    fn a_homopolymer_length_error_is_one_indel() {
        let ops = align(b"ACGTTAAAAAGCATGC", b"ACGTTAAAAGCATGC");
        assert_eq!(counts(&ops), [15, 0, 0, 1]);
    }

    #[test]
    fn non_acgt_bases_mismatch() {
        let ops = align(b"ACGTNACGT", b"ACGTNACGT");
        assert_eq!(counts(&ops), [8, 1, 0, 0]);
    }

    #[test]
    fn repeated_calls_give_the_same_columns() {
        let mut aligner = PairAligner::new(PoaParams::default()).expect("default parameters");
        let (target, query) = (b"ACGTTGCATGCATGCAAGTC", b"ACGTGCATGCATTGCAAGTC");
        let (mut first, mut again) = (Vec::new(), Vec::new());
        assert!(aligner.align(target, query, &mut first));
        assert!(aligner.align(b"TTTT", b"TTAT", &mut again));
        assert!(aligner.align(target, query, &mut again));
        assert_eq!(first, again);
    }

    #[test]
    fn empty_sequences_are_all_gaps() {
        assert_eq!(align(b"ACG", b""), vec![PairOp::Deletion; 3]);
        assert_eq!(align(b"", b"AC"), vec![PairOp::Insertion; 2]);
    }
}
