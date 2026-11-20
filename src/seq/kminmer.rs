//! k-min-mers: windows of consecutive anchors, the nodes of the minimizer-space graph.

use crate::seq::{Anchor, Strand, combine};

/// A window of consecutive anchors of one read.
///
/// The derived order compares `hash`, then `strand`, then `first_anchor`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct KMinMer {
    /// Canonical hash: the smaller of the hashes of the anchor hashes read forward and backward.
    pub hash: u64,
    /// [`Strand::Forward`] if the forward reading produced [`KMinMer::hash`].
    pub strand: Strand,
    /// Index, in the read's anchor list, of the window's first anchor.
    pub first_anchor: u32,
}

/// Replaces the contents of `out` with the k-min-mers of `anchors`, `window` anchors each.
///
/// A k-min-mer covers anchors `first_anchor..first_anchor + window`, and `out` lists them in
/// increasing [`KMinMer::first_anchor`] order. Its hash is the smaller of [`combine`] over the
/// anchor hashes read forward and read backward, and its strand says which reading won.
///
/// The k-min-mer at position `i` of `out` need not start at anchor `i`: windows whose two
/// readings hash equally (a palindromic hash tuple such as `h1 h2 h1`) have no orientation and
/// are skipped. A window of one anchor is the exception: its two readings always agree, so it
/// takes the anchor's own strand and is never skipped.
///
/// Only anchor hashes enter the hash, not anchor strands, as in mdBG. Windows with the same hash
/// tuple but different anchor strands, which inverted repeats at anchor scale can produce,
/// therefore share a hash, and a palindromic hash tuple is skipped even when its strands would
/// orient it.
///
/// # Panics
///
/// Panics if `window` is zero or `anchors` has more than `u32::MAX` entries. Neither happens
/// with [`SketchParams::window`](crate::seq::SketchParams::window) and the output of
/// [`Sketcher::sketch`](crate::seq::Sketcher::sketch), which never exceeds the read length.
pub fn kminmers(anchors: &[Anchor], window: usize, out: &mut Vec<KMinMer>) {
    assert!(window > 0, "k-min-mer window must be at least 1");
    out.clear();
    for (first, chunk) in anchors.windows(window).enumerate() {
        let (hash, strand) = if let [anchor] = chunk {
            (combine([anchor.hash]), anchor.strand)
        } else {
            let forward = combine(chunk.iter().map(|a| a.hash));
            let backward = combine(chunk.iter().rev().map(|a| a.hash));
            match forward.cmp(&backward) {
                std::cmp::Ordering::Less => (forward, Strand::Forward),
                std::cmp::Ordering::Greater => (backward, Strand::Reverse),
                std::cmp::Ordering::Equal => continue,
            }
        };
        let first_anchor = u32::try_from(first)
            .unwrap_or_else(|_| panic!("too many anchors in one read: {}", anchors.len()));
        out.push(KMinMer {
            hash,
            strand,
            first_anchor,
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn anchors(hashes: &[u64]) -> Vec<Anchor> {
        hashes
            .iter()
            .zip(0..)
            .map(|(&hash, i)| Anchor {
                hash,
                strand: Strand::Forward,
                raw_start: i * 10,
                raw_end: i * 10 + 5,
            })
            .collect()
    }

    #[test]
    fn reverse_read_gives_same_kminmers_flipped() {
        let hashes = [5, 9, 2, 7, 11, 3];
        let forward_anchors = anchors(&hashes);
        let mut reversed = hashes;
        reversed.reverse();
        let reverse_anchors = anchors(&reversed);

        let mut forward = Vec::new();
        kminmers(&forward_anchors, 3, &mut forward);
        let mut reverse = Vec::new();
        kminmers(&reverse_anchors, 3, &mut reverse);

        assert_eq!(forward.len(), 4);
        let n = u32::try_from(hashes.len() - 3).unwrap();
        let mirrored: Vec<KMinMer> = reverse
            .iter()
            .rev()
            .map(|m| KMinMer {
                hash: m.hash,
                strand: !m.strand,
                first_anchor: n - m.first_anchor,
            })
            .collect();
        assert_eq!(forward, mirrored);
    }

    #[test]
    fn skips_palindromic_windows_and_short_reads() {
        let mut out = Vec::new();
        kminmers(&anchors(&[1, 2, 1]), 3, &mut out);
        assert!(out.is_empty());
        kminmers(&anchors(&[1, 2]), 3, &mut out);
        assert!(out.is_empty());
    }

    #[test]
    fn window_of_one_keeps_every_anchor_with_its_strand() {
        let mut input = anchors(&[4, 8, 15]);
        input[1].strand = Strand::Reverse;
        let mut out = Vec::new();
        kminmers(&input, 1, &mut out);
        let expected: Vec<KMinMer> = input
            .iter()
            .zip(0..)
            .map(|(a, i)| KMinMer {
                hash: combine([a.hash]),
                strand: a.strand,
                first_anchor: i,
            })
            .collect();
        assert_eq!(out, expected);
    }
}
