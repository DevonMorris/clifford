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
mod euclidean3;
#[cfg(feature = "native")]
mod menu;
mod minkowski2_diagram;
mod minkowski2_dilation;
mod projective2;
mod projective2_robot;
mod test_3d;

pub use complex_domain::ComplexDomainDemo;
pub use complex_fractal::ComplexFractalDemo;
pub use conformal2_circles::Conformal2CirclesDemo;
pub use conformal2_intersection::Conformal2IntersectionDemo;
pub use conformal2_inversion::Conformal2InversionDemo;
pub use conformal2_mobius::Conformal2MobiusDemo;
pub use dual_autodiff::DualAutodiffDemo;
pub use euclidean2::Euclidean2Demo;
pub use euclidean3::Euclidean3Demo;
#[cfg(feature = "native")]
pub use menu::DemoMenu;
pub use minkowski2_diagram::Minkowski2DiagramDemo;
pub use minkowski2_dilation::Minkowski2DilationDemo;
pub use projective2::Projective2Demo;
pub use projective2_robot::RobotArmDemo;
pub use test_3d::Test3DDemo;
