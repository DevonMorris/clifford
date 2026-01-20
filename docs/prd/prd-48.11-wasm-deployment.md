# PRD-48.11: WASM Web Deployment

**Status**: Complete
**Parent**: PRD-48
**Depends on**: PRD-48.1, PRD-48.10
**Goal**: Deploy interactive demos to the web via WebAssembly

## Status

WASM deployment is complete for all 2D demos. The demos are accessible at the GitHub Pages URL.

**Note**: When 3D demos are added (using three-d), WASM compatibility will need to be validated and potentially updated.

## Motivation

Web deployment maximizes accessibility:
- No installation required
- Shareable links to specific demos
- Works on any device with a browser
- Great for documentation and marketing

## Build Tool: Trunk

[Trunk](https://trunkrs.dev/) is the standard tool for building Rust WASM apps:

```bash
# Install
cargo install trunk

# Build
trunk build --release

# Serve locally
trunk serve
```

## Crate Configuration

### Cargo.toml Updates

```toml
# crates/clifford-viz/Cargo.toml

[package]
name = "clifford-viz"
# ...

[lib]
crate-type = ["cdylib", "rlib"]  # cdylib for WASM

[dependencies]
clifford = { path = "../clifford" }
eframe = { version = "0.30", default-features = false, features = [
    "default_fonts",
    "glow",          # OpenGL backend (works in WASM)
    "persistence",   # Save state to localStorage
] }
egui = "0.30"
egui_plot = "0.30"

[target.'cfg(target_arch = "wasm32")'.dependencies]
wasm-bindgen-futures = "0.4"
web-sys = { version = "0.3", features = ["Window", "Document"] }

[features]
default = ["native"]
native = ["eframe/default"]
web = ["eframe/wasm"]
```

### index.html Template

```html
<!-- crates/clifford-viz/index.html -->
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Clifford Visualization</title>
    <link data-trunk rel="css" href="style.css">
    <link data-trunk rel="icon" href="favicon.ico">
</head>
<body>
    <div id="loading">Loading...</div>
    <link data-trunk rel="rust" data-wasm-opt="z" data-bin="demo_launcher"/>
</body>
</html>
```

### Demo Launcher

Single entry point that can launch any demo:

```rust
// crates/clifford-viz/src/bin/demo_launcher.rs

use eframe::wasm_bindgen::JsCast;
use clifford_viz::demos::*;

fn main() {
    // Get demo name from URL query param: ?demo=euclidean2
    let demo_name = get_query_param("demo").unwrap_or("menu".to_string());

    let options = eframe::WebOptions::default();

    wasm_bindgen_futures::spawn_local(async move {
        eframe::WebRunner::new()
            .start(
                "canvas",
                options,
                Box::new(|cc| Ok(create_demo(&demo_name, cc))),
            )
            .await
            .expect("Failed to start eframe");
    });
}

fn create_demo(name: &str, cc: &eframe::CreationContext) -> Box<dyn eframe::App> {
    match name {
        "euclidean2" => Box::new(Euclidean2Demo::new(cc)),
        "euclidean3" => Box::new(Euclidean3Demo::new(cc)),
        "projective2" => Box::new(Projective2Demo::new(cc)),
        "projective3" => Box::new(Projective3Demo::new(cc)),
        "conformal3" => Box::new(Conformal3Demo::new(cc)),
        "quaternion" => Box::new(QuaternionDemo::new(cc)),
        "complex" => Box::new(ComplexDemo::new(cc)),
        // ... etc
        _ => Box::new(DemoMenu::new(cc)),  // Show menu if unknown
    }
}

fn get_query_param(name: &str) -> Option<String> {
    let window = web_sys::window()?;
    let search = window.location().search().ok()?;
    // Parse query string...
    todo!()
}
```

### Demo Menu

Landing page with links to all demos:

```rust
// crates/clifford-viz/src/demos/menu.rs

pub struct DemoMenu;

impl eframe::App for DemoMenu {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.heading("Clifford Algebra Visualizations");
            ui.separator();

            ui.heading("Euclidean Geometry");
            if ui.link("2D Rotors").clicked() {
                navigate_to("?demo=euclidean2");
            }
            if ui.link("3D Rotors & Gimbal Lock").clicked() {
                navigate_to("?demo=euclidean3");
            }

            ui.heading("Projective Geometry (PGA)");
            if ui.link("2D Point-Line Geometry").clicked() {
                navigate_to("?demo=projective2");
            }
            // ... etc
        });
    }
}
```

## GitHub Pages Deployment

### Directory Structure

```
docs/                      # GitHub Pages root
  demos/
    index.html             # Demo launcher
    demo_launcher.js       # Generated WASM bindings
    demo_launcher_bg.wasm  # Compiled WASM
    style.css
```

### GitHub Actions Workflow

```yaml
# .github/workflows/deploy-demos.yml
name: Deploy Web Demos

on:
  push:
    branches: [main]
    paths:
      - 'crates/clifford-viz/**'
  workflow_dispatch:  # Manual trigger

permissions:
  contents: read
  pages: write
  id-token: write

jobs:
  build:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - name: Install Rust
        uses: dtolnay/rust-toolchain@stable
        with:
          targets: wasm32-unknown-unknown

      - name: Install Trunk
        run: cargo install trunk

      - name: Build WASM
        run: |
          cd crates/clifford-viz
          trunk build --release --public-url /clifford/demos/

      - name: Upload artifact
        uses: actions/upload-pages-artifact@v3
        with:
          path: crates/clifford-viz/dist

  deploy:
    needs: build
    runs-on: ubuntu-latest
    environment:
      name: github-pages
      url: ${{ steps.deployment.outputs.page_url }}
    steps:
      - name: Deploy to GitHub Pages
        id: deployment
        uses: actions/deploy-pages@v4
```

### Repository Settings

Enable GitHub Pages:
1. Settings → Pages
2. Source: GitHub Actions

Demo URL: `https://devonmorris.github.io/clifford/demos/?demo=euclidean2`

## Trunk Configuration

```toml
# crates/clifford-viz/Trunk.toml
[build]
target = "index.html"
dist = "dist"

[watch]
ignore = ["dist", "target"]

[[hooks]]
stage = "pre_build"
command = "sh"
command_arguments = ["-c", "echo 'Building WASM...'"]

[tools]
wasm_opt = "version_118"
```

## Performance Considerations

### WASM Size Optimization

```toml
# crates/clifford-viz/Cargo.toml

[profile.release]
opt-level = "z"      # Optimize for size
lto = true           # Link-time optimization
codegen-units = 1    # Better optimization
strip = true         # Strip symbols

[profile.release.package."*"]
opt-level = "z"
```

### Lazy Loading (Future)

For many demos, consider code splitting:
- Load only the selected demo
- Use dynamic imports

## Local Development

```bash
cd crates/clifford-viz

# Install trunk
cargo install trunk

# Serve with hot reload
trunk serve --open

# Build for production
trunk build --release
```

## Implementation Tasks

1. [ ] Add `web` feature to clifford-viz Cargo.toml
2. [ ] Create `index.html` template
3. [ ] Create `demo_launcher.rs` binary
4. [ ] Create `DemoMenu` landing page
5. [ ] Create `Trunk.toml` configuration
6. [ ] Add deploy-demos.yml workflow
7. [ ] Enable GitHub Pages in repo settings
8. [ ] Add WASM size optimizations to Cargo.toml
9. [ ] Test all demos in browser (Chrome, Firefox, Safari)
10. [ ] Add link to web demos in main README

## Verification

```bash
# Local test
cd crates/clifford-viz
trunk serve --open
# Navigate to http://localhost:8080/?demo=euclidean2

# Build size check
trunk build --release
ls -lh dist/*.wasm  # Should be < 5MB ideally
```

### Browser Compatibility

Test in:
- [ ] Chrome (latest)
- [ ] Firefox (latest)
- [ ] Safari (latest)
- [ ] Edge (latest)
- [ ] Mobile Chrome/Safari

## Demo URLs

After deployment:

| Demo | URL |
|------|-----|
| Menu | `https://devonmorris.github.io/clifford/demos/` |
| Euclidean 2D | `.../demos/?demo=euclidean2` |
| Euclidean 3D | `.../demos/?demo=euclidean3` |
| Projective 2D | `.../demos/?demo=projective2` |
| Projective 3D | `.../demos/?demo=projective3` |
| Conformal 3D | `.../demos/?demo=conformal3` |
| Quaternion | `.../demos/?demo=quaternion` |
| Complex | `.../demos/?demo=complex` |
| ... | ... |
