//! Demo implementations for clifford-viz.
//!
//! Each demo demonstrates a specific geometric algebra concept.
//! These can be run natively or compiled to WASM for web deployment.

mod euclidean2;
mod menu;
mod projective2;
mod projective2_robot;

pub use euclidean2::Euclidean2Demo;
pub use menu::DemoMenu;
pub use projective2::Projective2Demo;
pub use projective2_robot::RobotArmDemo;
