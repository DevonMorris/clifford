//! Demo implementations for clifford-viz.
//!
//! Each demo demonstrates a specific geometric algebra concept.
//! These can be run natively or compiled to WASM for web deployment.

mod complex_domain;
mod complex_fractal;
mod conformal2_circles;
mod conformal2_intersection;
mod conformal2_inversion;
mod conformal2_mobius;
mod dual_autodiff;
mod euclidean2;
#[cfg(feature = "three-d")]
mod euclidean3;
mod minkowski2_diagram;
mod minkowski2_dilation;
mod projective2;
mod projective2_robot;
#[cfg(feature = "three-d")]
mod projective3_lines;
#[cfg(feature = "three-d")]
mod projective3_motor;
#[cfg(feature = "three-d")]
mod projective3_robot;

pub use complex_domain::ComplexDomainDemo;
pub use complex_fractal::ComplexFractalDemo;
pub use conformal2_circles::Conformal2CirclesDemo;
pub use conformal2_intersection::Conformal2IntersectionDemo;
pub use conformal2_inversion::Conformal2InversionDemo;
pub use conformal2_mobius::Conformal2MobiusDemo;
pub use dual_autodiff::DualAutodiffDemo;
pub use euclidean2::Euclidean2Demo;
#[cfg(feature = "three-d")]
pub use euclidean3::Euclidean3Demo;
pub use minkowski2_diagram::Minkowski2DiagramDemo;
pub use minkowski2_dilation::Minkowski2DilationDemo;
pub use projective2::Projective2Demo;
pub use projective2_robot::RobotArmDemo;
#[cfg(feature = "three-d")]
pub use projective3_lines::Projective3LinesDemo;
#[cfg(feature = "three-d")]
pub use projective3_motor::Projective3MotorDemo;
#[cfg(feature = "three-d")]
pub use projective3_robot::Projective3RobotDemo;
