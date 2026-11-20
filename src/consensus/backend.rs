//! The batched consensus backend interface.

use crate::consensus::Window;

/// Computes base-level consensus sequences for batches of windows.
///
/// Backends work on batches so that a backend with a high per-call cost, such as a GPU, can
/// amortise it over many windows. Implementations must be deterministic: the output for a window
/// may depend only on that window and the backend's parameters, never on the rest of the batch,
/// the thread count or scheduling.
///
/// The output for a window may depend on the order of its segments (backends pick a backbone
/// and break ties by it), so a caller that needs reproducible output must fill each [`Window`]
/// in a canonical order; see [`Window`].
pub trait ConsensusBackend: Send + Sync {
    /// Returns one consensus sequence per window: `output[i]` belongs to `windows[i]`.
    ///
    /// The output has exactly `windows.len()` entries. A window without any non-empty segment
    /// yields an empty sequence.
    #[must_use]
    fn consensus(&self, windows: &[Window]) -> Vec<Vec<u8>>;
}
