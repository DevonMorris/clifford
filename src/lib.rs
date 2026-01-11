#![doc = include_str!("../README.md")]
#![doc(
    html_logo_url = "https://raw.githubusercontent.com/DevonMorris/clifford/main/assets/clifford.png",
    html_favicon_url = "https://raw.githubusercontent.com/DevonMorris/clifford/main/assets/clifford.png"
)]

// Ensure nalgebra feature flags are mutually exclusive
#[cfg(any(
    all(feature = "nalgebra-0_32", feature = "nalgebra-0_33"),
    all(feature = "nalgebra-0_32", feature = "nalgebra-0_34"),
    all(feature = "nalgebra-0_33", feature = "nalgebra-0_34"),
))]
compile_error!(
    "Features `nalgebra-0_32`, `nalgebra-0_33`, and `nalgebra-0_34` are mutually exclusive. Enable only one."
);

pub mod algebra;
pub mod basis;
pub mod prelude;
pub mod scalar;
pub mod signature;
pub mod specialized;

/// Test utilities available only during testing.
#[cfg(test)]
pub(crate) mod test_utils {
    /// Standard epsilon for absolute difference comparisons in tests.
    ///
    /// Use this constant instead of magic numbers like `1e-10` or `1e-9`.
    /// This value is chosen to be strict enough for most operations while
    /// allowing for reasonable floating-point accumulation in compound operations.
    pub const ABS_DIFF_EQ_EPS: f64 = 1e-10;
}
