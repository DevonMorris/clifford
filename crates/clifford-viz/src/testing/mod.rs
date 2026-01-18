//! Testing utilities for visualization correctness.
//!
//! This module provides tools for verifying that visualizations correctly
//! represent the underlying geometric algebra operations. Rather than using
//! golden image tests (which only detect change, not correctness), we focus on:
//!
//! - **Coordinate assertions**: Verify algebraic coordinates map correctly
//! - **Visual invariants**: Property-based tests for geometric properties
//! - **Scene graph assertions**: Test that primitives exist with correct properties
//!
//! # Sample Limiting for Graphics Tests
//!
//! Graphics tests can be expensive, so we limit proptest samples. Use the
//! [`VIZ_PROPTEST_CASES`] constant or create a config with [`viz_proptest_config`].
//!
//! ```ignore
//! use clifford_viz::testing::{viz_proptest_config, VIZ_PROPTEST_CASES};
//!
//! proptest! {
//!     #![proptest_config(viz_proptest_config())]
//!     #[test]
//!     fn my_test(x in -10.0f32..10.0) {
//!         // Test with limited samples
//!     }
//! }
//! ```

mod inspector;
mod invariants;
mod scene;

pub use inspector::*;
pub use invariants::*;
pub use scene::*;

use proptest::test_runner::Config as ProptestConfig;

/// Default number of proptest cases for visualization tests.
///
/// We use a lower number than the default (256) because graphics tests
/// can involve more computation. This balances coverage with test speed.
pub const VIZ_PROPTEST_CASES: u32 = 32;

/// Create a proptest config suitable for visualization tests.
///
/// Uses [`VIZ_PROPTEST_CASES`] as the case count.
#[must_use]
pub fn viz_proptest_config() -> ProptestConfig {
    ProptestConfig::with_cases(VIZ_PROPTEST_CASES)
}

/// Create a proptest config with a custom case count.
#[must_use]
pub fn viz_proptest_config_with_cases(cases: u32) -> ProptestConfig {
    ProptestConfig::with_cases(cases)
}
