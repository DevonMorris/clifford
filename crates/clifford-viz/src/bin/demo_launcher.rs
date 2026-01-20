//! Web demo launcher - single WASM entry point for all demos.
//!
//! Loads the appropriate demo based on URL query parameter.
//! Usage: `https://example.com/demos/?demo=euclidean2`
//!
//! This binary is only meant to be compiled for WASM. Running it natively
//! will print an error message.

#[cfg(not(target_arch = "wasm32"))]
fn main() {
    eprintln!("This binary is only meant to be run in a web browser via WASM.");
    eprintln!("For native demos, use:");
    eprintln!("  cargo run -p clifford-viz --example euclidean2 --release");
    eprintln!("  cargo run -p clifford-viz --example projective2 --release");
    eprintln!("  cargo run -p clifford-viz --example projective2_robot --release");
    eprintln!("  cargo run -p clifford-viz --example conformal2_circles --release");
    eprintln!("  cargo run -p clifford-viz --example conformal2_inversion --release");
    eprintln!("  cargo run -p clifford-viz --example conformal2_mobius --release");
}

#[cfg(target_arch = "wasm32")]
use clifford_viz::common::prelude::AppWrapper;
#[cfg(target_arch = "wasm32")]
use clifford_viz::demos::{
    ComplexDomainDemo, ComplexFractalDemo, Conformal2CirclesDemo, Conformal2IntersectionDemo,
    Conformal2InversionDemo, Conformal2MobiusDemo, DemoMenu, DualAutodiffDemo, Euclidean2Demo,
    Euclidean3Demo, Minkowski2DiagramDemo, Minkowski2DilationDemo, Projective2Demo, RobotArmDemo,
    Test3DDemo,
};
#[cfg(target_arch = "wasm32")]
use wasm_bindgen::JsCast;

#[cfg(target_arch = "wasm32")]
fn main() {
    // Configure logging for browser console
    eframe::WebLogger::init(log::LevelFilter::Debug).ok();

    let demo_name = get_demo_from_url().unwrap_or_else(|| "menu".to_string());

    let web_options = eframe::WebOptions::default();

    wasm_bindgen_futures::spawn_local(async move {
        // Get the canvas element from the DOM
        let document = web_sys::window()
            .expect("No window")
            .document()
            .expect("No document");
        let canvas = document
            .get_element_by_id("the_canvas_id")
            .expect("No canvas element with id 'the_canvas_id'");
        let canvas: web_sys::HtmlCanvasElement = canvas
            .dyn_into()
            .expect("Element is not an HtmlCanvasElement");

        let result = eframe::WebRunner::new()
            .start(
                canvas,
                web_options,
                Box::new(move |cc| Ok(create_app(&demo_name, cc))),
            )
            .await;

        // Remove loading indicator
        if let Some(loading) = web_sys::window()
            .and_then(|w| w.document())
            .and_then(|d| d.get_element_by_id("loading"))
        {
            loading.remove();
        }

        if let Err(e) = result {
            log::error!("Failed to start eframe: {:?}", e);
        }
    });
}

/// Parse demo name from URL query string.
#[cfg(target_arch = "wasm32")]
fn get_demo_from_url() -> Option<String> {
    let window = web_sys::window()?;
    let search = window.location().search().ok()?;

    // Parse "?demo=name" or "?demo=name&other=value"
    search.strip_prefix('?').and_then(|s| {
        s.split('&').find_map(|pair| {
            let mut parts = pair.splitn(2, '=');
            match (parts.next(), parts.next()) {
                (Some("demo"), Some(value)) => Some(value.to_string()),
                _ => None,
            }
        })
    })
}

/// Create the appropriate demo based on name.
#[cfg(target_arch = "wasm32")]
fn create_app(name: &str, cc: &eframe::CreationContext<'_>) -> Box<dyn eframe::App> {
    match name {
        // Euclidean demos
        "euclidean2" => Box::new(AppWrapper::<Euclidean2Demo>::default()),
        "euclidean3" => Box::new(AppWrapper::<Euclidean3Demo>::default()),
        // Projective demos
        "projective2" => Box::new(AppWrapper::<Projective2Demo>::default()),
        "projective2_robot" => Box::new(AppWrapper::<RobotArmDemo>::default()),
        // Conformal demos
        "conformal2_circles" => Box::new(AppWrapper::<Conformal2CirclesDemo>::default()),
        "conformal2_inversion" => Box::new(AppWrapper::<Conformal2InversionDemo>::default()),
        "conformal2_mobius" => Box::new(AppWrapper::<Conformal2MobiusDemo>::default()),
        "conformal2_intersection" => Box::new(AppWrapper::<Conformal2IntersectionDemo>::default()),
        // Complex and dual number demos
        "complex_domain" => Box::new(AppWrapper::<ComplexDomainDemo>::default()),
        "complex_fractal" => Box::new(AppWrapper::<ComplexFractalDemo>::default()),
        "dual_autodiff" => Box::new(AppWrapper::<DualAutodiffDemo>::default()),
        // Minkowski demos
        "minkowski2_diagram" => Box::new(AppWrapper::<Minkowski2DiagramDemo>::default()),
        "minkowski2_dilation" => Box::new(AppWrapper::<Minkowski2DilationDemo>::default()),
        // 3D test demo
        "test_3d" => Box::new(AppWrapper::<Test3DDemo>::default()),
        // Default: show menu
        _ => Box::new(DemoMenu::new(cc)),
    }
}
