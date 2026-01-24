//! 3D Point-Line-Plane operations visualization using native three-d rendering.
//!
//! Run with:
//! ```bash
//! cargo run -p clifford-viz --example projective3_geometry --features three-d --release
//! ```

use clifford_viz::common::app_three_d::run_three_d_app_3d;
use clifford_viz::demos::Projective3GeometryDemo;

fn main() {
    run_three_d_app_3d::<Projective3GeometryDemo>("3D PGA Point-Line-Plane");
}
