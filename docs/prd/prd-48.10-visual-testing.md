# PRD-48.10: Visual Regression Testing

**Status**: Draft
**Parent**: PRD-48
**Depends on**: PRD-48.1
**Goal**: Automated visual testing to catch rendering regressions

## Overview

Visual output is the primary deliverable of the visualization demos. Unit tests can verify calculations, but only visual tests can catch:
- Rendering artifacts
- Color/style regressions
- Layout issues
- Animation glitches at specific frames

## Approach: Golden Image Comparison

1. **Generate golden images**: Run demos in deterministic state, capture screenshots
2. **CI comparison**: Re-run demos, capture new screenshots, compare pixel-by-pixel
3. **Threshold tolerance**: Allow small differences (anti-aliasing, font rendering)
4. **Failure artifacts**: Save diff images for human review

## Tools

### Screenshot Capture

**Linux (CI and local)**:
```bash
# scrot - simple and reliable
scrot --focused --delay 1 screenshot.png

# maim - more features
maim --window $(xdotool getactivewindow) screenshot.png

# For headless CI: Xvfb (virtual framebuffer)
Xvfb :99 -screen 0 1024x768x24 &
export DISPLAY=:99
cargo run --example euclidean2 &
sleep 2
scrot screenshot.png
```

**egui built-in** (preferred for determinism):
```rust
// egui can render to an image buffer directly
use egui::ColorImage;

fn capture_frame(ctx: &egui::Context) -> ColorImage {
    // Use egui's frame capture API
    ctx.tessellate(shapes)
    // Render to image...
}
```

### Image Comparison

**Rust crates**:
```toml
[dev-dependencies]
image = "0.25"
img_diff = "0.2"  # or pixelmatch-rs
```

**Comparison implementation**:
```rust
use image::{GenericImageView, Rgba};

/// Compare two images, return percentage of differing pixels
fn compare_images(golden: &str, actual: &str, threshold: u8) -> f64 {
    let golden_img = image::open(golden).unwrap();
    let actual_img = image::open(actual).unwrap();

    assert_eq!(golden_img.dimensions(), actual_img.dimensions(),
        "Image dimensions must match");

    let (width, height) = golden_img.dimensions();
    let total_pixels = (width * height) as f64;
    let mut diff_pixels = 0u64;

    for y in 0..height {
        for x in 0..width {
            let g = golden_img.get_pixel(x, y);
            let a = actual_img.get_pixel(x, y);
            if !pixels_similar(g, a, threshold) {
                diff_pixels += 1;
            }
        }
    }

    (diff_pixels as f64 / total_pixels) * 100.0
}

fn pixels_similar(a: Rgba<u8>, b: Rgba<u8>, threshold: u8) -> bool {
    (a[0] as i16 - b[0] as i16).abs() <= threshold as i16 &&
    (a[1] as i16 - b[1] as i16).abs() <= threshold as i16 &&
    (a[2] as i16 - b[2] as i16).abs() <= threshold as i16
}
```

## Test Infrastructure

### Directory Structure

```
tests/
  visual/
    mod.rs                    # Visual test harness
    golden/                   # Reference images (committed to git)
      euclidean2/
        rotor_0deg.png
        rotor_45deg.png
        rotor_90deg.png
      projective2/
        point_line_join.png
        motor_animation_frame0.png
      ...
    actual/                   # Generated during test (gitignored)
    diff/                     # Diff images on failure (gitignored)
```

### Test Harness

```rust
// tests/visual/mod.rs

use std::path::PathBuf;

pub struct VisualTest {
    name: String,
    golden_dir: PathBuf,
    actual_dir: PathBuf,
    diff_dir: PathBuf,
    tolerance_percent: f64,  // Max allowed difference
    pixel_threshold: u8,     // Per-pixel color threshold
}

impl VisualTest {
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            golden_dir: PathBuf::from("tests/visual/golden"),
            actual_dir: PathBuf::from("tests/visual/actual"),
            diff_dir: PathBuf::from("tests/visual/diff"),
            tolerance_percent: 0.1,  // 0.1% pixel difference allowed
            pixel_threshold: 2,      // RGB values can differ by 2
        }
    }

    pub fn compare(&self, test_name: &str) -> Result<(), VisualTestError> {
        let golden = self.golden_dir.join(&self.name).join(format!("{}.png", test_name));
        let actual = self.actual_dir.join(&self.name).join(format!("{}.png", test_name));

        if !golden.exists() {
            return Err(VisualTestError::MissingGolden(golden));
        }

        let diff_percent = compare_images(
            golden.to_str().unwrap(),
            actual.to_str().unwrap(),
            self.pixel_threshold
        );

        if diff_percent > self.tolerance_percent {
            // Generate diff image for debugging
            let diff_path = self.diff_dir.join(&self.name).join(format!("{}_diff.png", test_name));
            generate_diff_image(&golden, &actual, &diff_path);

            return Err(VisualTestError::DifferenceExceeded {
                test: test_name.to_string(),
                expected: self.tolerance_percent,
                actual: diff_percent,
                diff_image: diff_path,
            });
        }

        Ok(())
    }
}

#[derive(Debug)]
pub enum VisualTestError {
    MissingGolden(PathBuf),
    DifferenceExceeded {
        test: String,
        expected: f64,
        actual: f64,
        diff_image: PathBuf,
    },
}
```

### Deterministic Rendering

For reproducible screenshots, demos must support deterministic mode:

```rust
// examples/visualization/common/app.rs

pub struct DeterministicConfig {
    pub window_size: (u32, u32),
    pub animation_time: f32,      // Fixed time for animation state
    pub random_seed: Option<u64>, // If any randomness
}

pub trait VisualizationApp {
    // ... existing methods ...

    /// Set up deterministic state for visual testing
    fn set_deterministic(&mut self, config: DeterministicConfig) {
        // Default: no-op, override in demos that need it
    }

    /// Capture current frame as image (for visual testing)
    fn capture_frame(&self, ctx: &egui::Context) -> egui::ColorImage;
}
```

### Example Visual Test

```rust
// tests/visual_euclidean2.rs

use clifford_visualization::euclidean2::Euclidean2Demo;
use visual_test::VisualTest;

#[test]
fn test_euclidean2_rotor_angles() {
    let test = VisualTest::new("euclidean2");

    // Test rotor at 0°
    let mut demo = Euclidean2Demo::default();
    demo.set_deterministic(DeterministicConfig {
        window_size: (800, 600),
        animation_time: 0.0,
        random_seed: None,
    });
    demo.angle = 0.0;
    demo.capture_to("tests/visual/actual/euclidean2/rotor_0deg.png");
    test.compare("rotor_0deg").unwrap();

    // Test rotor at 45°
    demo.angle = std::f32::consts::FRAC_PI_4;
    demo.capture_to("tests/visual/actual/euclidean2/rotor_45deg.png");
    test.compare("rotor_45deg").unwrap();

    // Test rotor at 90°
    demo.angle = std::f32::consts::FRAC_PI_2;
    demo.capture_to("tests/visual/actual/euclidean2/rotor_90deg.png");
    test.compare("rotor_90deg").unwrap();
}
```

## CI Integration

### GitHub Actions Workflow

```yaml
# .github/workflows/visual-tests.yml
name: Visual Tests

on:
  push:
    paths:
      - 'examples/visualization/**'
      - 'tests/visual/**'
  pull_request:
    paths:
      - 'examples/visualization/**'
      - 'tests/visual/**'

jobs:
  visual-tests:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - name: Install dependencies
        run: |
          sudo apt-get update
          sudo apt-get install -y xvfb libxkbcommon-x11-0 libgl1-mesa-dri

      - name: Run visual tests
        run: |
          xvfb-run --auto-servernum cargo test --test visual_tests

      - name: Upload diff images on failure
        if: failure()
        uses: actions/upload-artifact@v4
        with:
          name: visual-test-diffs
          path: tests/visual/diff/
          retention-days: 7
```

### Local Golden Image Update

```bash
# Script to update golden images after intentional changes
# scripts/bless-visual-tests.sh

#!/bin/bash
set -e

echo "Blessing visual tests (updating baseline images)..."

# Run each demo in deterministic mode, capture screenshots
for demo in euclidean2 euclidean3 projective2 projective3; do
    echo "Capturing $demo..."
    cargo run --example $demo -- --bless
done

echo "Baseline images updated in tests/visual/golden/"
echo "Review changes with: git diff --stat tests/visual/golden/"
```

## Demo Integration

Each demo needs a `--bless` flag for automated screenshot capture:

```rust
// examples/visualization/euclidean2.rs

fn main() {
    let args: Vec<String> = std::env::args().collect();

    if args.contains(&"--bless".to_string()) {
        capture_golden_images();
        return;
    }

    // Normal interactive mode
    run_app::<Euclidean2Demo>()
}

fn capture_golden_images() {
    let test_cases = vec![
        ("rotor_0deg", 0.0_f32),
        ("rotor_45deg", std::f32::consts::FRAC_PI_4),
        ("rotor_90deg", std::f32::consts::FRAC_PI_2),
        ("rotor_180deg", std::f32::consts::PI),
    ];

    for (name, angle) in test_cases {
        let mut demo = Euclidean2Demo::default();
        demo.angle = angle;
        demo.set_deterministic(DeterministicConfig {
            window_size: (800, 600),
            animation_time: 0.0,
            random_seed: None,
        });

        let path = format!("tests/visual/golden/euclidean2/{}.png", name);
        demo.render_to_file(&path);
        println!("Captured: {}", path);
    }
}
```

## Test Cases Per Demo

### euclidean2
- `rotor_0deg.png` - No rotation
- `rotor_45deg.png` - 45° rotation
- `rotor_90deg.png` - 90° rotation
- `rotor_180deg.png` - 180° rotation
- `bivector_visible.png` - With bivector arc shown
- `animation_frame_25.png` - Animation at 25% progress

### projective2
- `empty.png` - Initial empty state
- `two_points.png` - Two points placed
- `point_line_join.png` - Line joining two points
- `line_intersection.png` - Two lines meeting at a point
- `motor_identity.png` - Motor at identity
- `motor_rotated.png` - Motor with rotation applied

### conformal3
- `sphere_only.png` - Inversion sphere
- `sphere_inversion.png` - Sphere being inverted
- `line_to_circle.png` - Line transforming to circle
- `mobius_animation_50.png` - Möbius transform at 50%

### (similar for all other demos)

## Implementation Tasks

1. [ ] Add `image` and comparison crate to dev-dependencies
2. [ ] Create `tests/visual/mod.rs` harness
3. [ ] Create golden image directory structure
4. [ ] Add `DeterministicConfig` to `VisualizationApp` trait
5. [ ] Add `--bless` flag to each demo
6. [ ] Implement `render_to_file()` for headless capture
7. [ ] Create `scripts/bless-visual-tests.sh`
8. [ ] Add GitHub Actions workflow
9. [ ] Generate initial golden images for each demo
10. [ ] Document golden image update process in CONTRIBUTING.md

## Verification

```bash
# Run visual tests locally
cargo test --test visual_tests

# Update golden images after intentional changes
./scripts/bless-visual-tests.sh
git diff tests/visual/golden/  # Review changes
git add tests/visual/golden/
git commit -m "Update golden images for XYZ change"
```

## Benefits

1. **Catch regressions**: Any visual change triggers test failure
2. **Review process**: Diff images make it easy to review visual changes
3. **Documentation**: Golden images serve as visual documentation
4. **Confidence**: Refactor rendering code with confidence
5. **Cross-platform**: Can detect platform-specific rendering issues

## Considerations

### Git LFS for Golden Images

Golden images can grow large. Consider Git LFS:

```bash
# .gitattributes
tests/visual/golden/**/*.png filter=lfs diff=lfs merge=lfs -text
```

### Platform Differences

Font rendering and anti-aliasing differ between platforms. Options:
1. Only run visual tests on Linux CI
2. Use platform-specific golden images
3. Increase tolerance threshold
4. Use custom fonts embedded in the demo

### Flaky Tests

Visual tests can be flaky due to timing. Mitigations:
1. Deterministic animation state (fixed time, not real-time)
2. Disable animations during capture
3. Multiple retries with exponential backoff
