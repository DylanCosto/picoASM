//! The read segments of one consensus window.

/// The read segments covering one consensus window, all oriented to the same strand.
///
/// Segments are plain ASCII bases. `A`, `C`, `G` and `T` are matched in either case; any other
/// byte (such as `N`) mismatches everything, including itself. Empty segments are allowed and
/// ignored by the backends.
///
/// # Segment order
///
/// The consensus depends on the order of the segments, not only on the set: the backbone is
/// chosen by length with ties broken by position, and even vote splits follow the backbone.
/// For output that is reproducible across runs and thread counts, add the segments in a
/// canonical order, for example sorted by read id or by the hash used to select them, never in
/// the order a parallel pass happened to collect them.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct Window {
    segments: Vec<Vec<u8>>,
}

impl Window {
    /// Creates a window from its segments.
    #[must_use]
    pub const fn new(segments: Vec<Vec<u8>>) -> Self {
        Self { segments }
    }

    /// Returns the segments, in the order they were added.
    #[must_use]
    pub fn segments(&self) -> &[Vec<u8>] {
        &self.segments
    }

    /// Appends a segment.
    pub fn push(&mut self, segment: Vec<u8>) {
        self.segments.push(segment);
    }

    /// Returns the number of segments, counting empty ones.
    #[must_use]
    pub fn len(&self) -> usize {
        self.segments.len()
    }

    /// Returns `true` if the window has no segments.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.segments.is_empty()
    }

    /// Consumes the window and returns its segments.
    #[must_use]
    pub fn into_segments(self) -> Vec<Vec<u8>> {
        self.segments
    }
}

impl From<Vec<Vec<u8>>> for Window {
    fn from(segments: Vec<Vec<u8>>) -> Self {
        Self::new(segments)
    }
}

impl FromIterator<Vec<u8>> for Window {
    fn from_iter<I: IntoIterator<Item = Vec<u8>>>(iter: I) -> Self {
        Self::new(iter.into_iter().collect())
    }
}

impl Extend<Vec<u8>> for Window {
    fn extend<I: IntoIterator<Item = Vec<u8>>>(&mut self, iter: I) {
        self.segments.extend(iter);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn collects_and_extends() {
        let mut window: Window = [b"ACGT".to_vec(), b"ACG".to_vec()].into_iter().collect();
        assert_eq!(window.len(), 2);
        window.extend([Vec::new()]);
        window.push(b"T".to_vec());
        assert_eq!(window.len(), 4);
        assert!(!window.is_empty());
        assert_eq!(window.segments()[3], b"T");
        assert_eq!(Window::from(window.clone().into_segments()), window);
        assert!(Window::default().is_empty());
    }
}
