//! Demo implementations for clifford-viz.
//!
//! Each demo demonstrates a specific geometric algebra concept.
//! These can be run natively or compiled to WASM for web deployment.

mod conformal2_circles;
mod conformal2_inversion;
mod euclidean2;
mod menu;
mod projective2;
mod projective2_robot;

pub use conformal2_circles::Conformal2CirclesDemo;
pub use conformal2_inversion::Conformal2InversionDemo;
pub use euclidean2::Euclidean2Demo;
pub use menu::DemoMenu;
pub use projective2::Projective2Demo;
pub use projective2_robot::RobotArmDemo;
