#![doc = "Generated code for the projective2 algebra."]
pub mod conversions;
pub mod traits;
pub mod types;

// Domain-specific extensions (hand-written)
mod extensions;

// Re-export types at module level
pub use types::*;
