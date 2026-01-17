# PRD-48.10: Visual Testing via Invariants and Coordinate Assertions

**Status**: Draft
**Parent**: PRD-48
**Depends on**: PRD-48.1
**Goal**: Automated testing that proves correctness, not just detects change

## Philosophy: Correctness Over Change Detection

Golden image ("blessed") tests are fundamentally change detectors, not correctness proofs. They can tell you *something changed* but not *whether the change is correct*. For a mathematical visualization library, we need tests that verify:

1. **Coordinate correctness**: A point at (3, 4) renders at the correct screen position
2. **Visual invariants**: Geometric properties hold visually (orthogonal lines look orthogonal)
3. **Scene graph correctness**: The right primitives exist with the right properties
4. **Mathematical relationships**: Visual output matches algebraic computation

## Testing Strategies

### 1. Coordinate Mapping Tests

Verify that algebraic coordinates map to correct screen coordinates:

```rust
#[test]
fn test_point_renders_at_correct_position() {
    let demo = Euclidean2Demo::new();
    let point = Vector::new(3.0, 4.0);

    // Get the screen position where this point would render
    let screen_pos = demo.world_to_screen(point);

    // Verify the mapping is correct given the viewport
    let expected = demo.viewport.world_to_screen(3.0, 4.0);
    assert_eq!(screen_pos, expected);
}

#[test]
fn test_rotor_rotates_vector_correctly() {
    let demo = Euclidean2Demo::new();
    let rotor = Rotor::from_angle(std::f32::consts::FRAC_PI_2);
    let input = Vector::new(1.0, 0.0);

    // Apply rotor algebraically
    let rotated = rotor.transform(&input);

    // Verify the visual representation matches
    let input_screen = demo.world_to_screen(input);
    let rotated_screen = demo.world_to_screen(rotated);

    // 90 degree rotation: (1, 0) -> (0, 1)
    relative_eq!(rotated.x(), 0.0, epsilon = 1e-6);
    relative_eq!(rotated.y(), 1.0, epsilon = 1e-6);
}
```

### 2. Visual Invariant Property Tests

Use proptest to verify geometric invariants hold across random inputs:

```rust
proptest! {
    /// Rotation preserves distance from origin
    #[test]
    fn rotation_preserves_distance(
        angle in -PI..PI,
        x in -10.0f32..10.0,
        y in -10.0f32..10.0,
    ) {
        let demo = Euclidean2Demo::new();
        let rotor = Rotor::from_angle(angle);
        let point = Vector::new(x, y);

        let original_dist = (x * x + y * y).sqrt();
        let rotated = rotor.transform(&point);
        let rotated_dist = (rotated.x().powi(2) + rotated.y().powi(2)).sqrt();

        prop_assert!(relative_eq!(original_dist, rotated_dist, epsilon = 1e-5));
    }

    /// Two points joined by a line should both be incident with that line
    #[test]
    fn join_creates_incident_line(
        x1 in -10.0f32..10.0,
        y1 in -10.0f32..10.0,
        x2 in -10.0f32..10.0,
        y2 in -10.0f32..10.0,
    ) {
        let p1 = Point::from_cartesian(x1, y1);
        let p2 = Point::from_cartesian(x2, y2);
        let line = p1.wedge(&p2);

        // Both points should be incident with the line
        // (regressive product should be near zero)
        let incident1 = p1.antiwedge(&line);
        let incident2 = p2.antiwedge(&line);

        prop_assert!(incident1.is_near_zero(1e-5));
        prop_assert!(incident2.is_near_zero(1e-5));
    }

    /// Parallel lines should not intersect (meet at ideal point)
    #[test]
    fn parallel_lines_meet_at_infinity(
        x in -10.0f32..10.0,
        y in -10.0f32..10.0,
        dx in -1.0f32..1.0,
        dy in -1.0f32..1.0,
        offset in 0.1f32..5.0,
    ) {
        let p1 = Point::from_cartesian(x, y);
        let p2 = Point::from_cartesian(x + dx, y + dy);
        let line1 = p1.wedge(&p2);

        // Offset line perpendicular to direction
        let perp_x = -dy * offset;
        let perp_y = dx * offset;
        let p3 = Point::from_cartesian(x + perp_x, y + perp_y);
        let p4 = Point::from_cartesian(x + dx + perp_x, y + dy + perp_y);
        let line2 = p3.wedge(&p4);

        let intersection = line1.antiwedge(&line2);

        // Intersection should be an ideal point (weight = 0)
        prop_assert!(intersection.is_ideal(1e-5));
    }
}
```

### 3. Scene Graph Assertions

Test the structure of what gets rendered, not the pixels:

```rust
#[test]
fn test_euclidean2_scene_structure() {
    let mut demo = Euclidean2Demo::default();
    demo.angle = std::f32::consts::FRAC_PI_4;

    let scene = demo.build_scene();

    // Verify expected primitives exist
    assert!(scene.has_vector("original"));
    assert!(scene.has_vector("rotated"));
    assert!(scene.has_arc("rotation_arc"));

    // Verify vector properties
    let original = scene.get_vector("original");
    assert_eq!(original.start, Point2::ORIGIN);
    assert_eq!(original.color, Color::BLUE);

    let rotated = scene.get_vector("rotated");
    assert_eq!(rotated.start, Point2::ORIGIN);
    assert_eq!(rotated.color, Color::GREEN);

    // Verify the rotated vector is actually rotated 45 degrees
    let angle_diff = original.angle_to(&rotated);
    relative_eq!(angle_diff, std::f32::consts::FRAC_PI_4, epsilon = 1e-5);
}

#[test]
fn test_projective2_join_scene() {
    let mut demo = Projective2Demo::default();
    demo.add_point(1.0, 2.0);
    demo.add_point(4.0, 5.0);

    let scene = demo.build_scene();

    // Two points should create a line
    assert_eq!(scene.point_count(), 2);
    assert_eq!(scene.line_count(), 1);

    // The line should pass through both points
    let line = scene.get_line(0);
    assert!(line.contains_point(1.0, 2.0, 1e-5));
    assert!(line.contains_point(4.0, 5.0, 1e-5));
}
```

### 4. Animation Frame Tests

Test that animations produce correct intermediate states:

```rust
#[test]
fn test_motor_interpolation_midpoint() {
    let demo = Projective3Demo::default();

    let start = Motor::identity();
    let end = Motor::from_translation(10.0, 0.0, 0.0);

    // At t=0.5, should be halfway
    let mid = start.slerp(&end, 0.5);
    let mid_translation = mid.extract_translation();

    relative_eq!(mid_translation.x(), 5.0, epsilon = 1e-5);
    relative_eq!(mid_translation.y(), 0.0, epsilon = 1e-5);
    relative_eq!(mid_translation.z(), 0.0, epsilon = 1e-5);
}

#[test]
fn test_animation_frame_sequence() {
    let mut demo = Euclidean2Demo::default();
    demo.set_animation_duration(1.0);

    // Frame at t=0
    demo.set_time(0.0);
    let scene_0 = demo.build_scene();
    let angle_0 = scene_0.get_rotor_angle();

    // Frame at t=0.5
    demo.set_time(0.5);
    let scene_50 = demo.build_scene();
    let angle_50 = scene_50.get_rotor_angle();

    // Frame at t=1.0
    demo.set_time(1.0);
    let scene_100 = demo.build_scene();
    let angle_100 = scene_100.get_rotor_angle();

    // Verify smooth interpolation
    assert!(angle_0 < angle_50);
    assert!(angle_50 < angle_100);
}
```

## Automated Visual Inspection

For issues that are hard to express as coordinate assertions, we provide automated first-pass inspection:

### Screenshot Analysis Tool

```rust
/// Automated visual analysis that flags potential issues
pub struct VisualInspector {
    tolerance: f32,
}

impl VisualInspector {
    /// Analyze a screenshot for common rendering issues
    pub fn inspect(&self, image: &ColorImage) -> Vec<VisualIssue> {
        let mut issues = Vec::new();

        // Check for blank/all-black frames
        if self.is_blank(image) {
            issues.push(VisualIssue::BlankFrame);
        }

        // Check for NaN artifacts (often appear as black or white pixels)
        if self.has_nan_artifacts(image) {
            issues.push(VisualIssue::PossibleNanArtifacts);
        }

        // Check for clipping (content outside viewport)
        if self.has_edge_clipping(image) {
            issues.push(VisualIssue::EdgeClipping);
        }

        // Check for z-fighting (flickering depth)
        if self.has_z_fighting_patterns(image) {
            issues.push(VisualIssue::PossibleZFighting);
        }

        issues
    }

    fn is_blank(&self, image: &ColorImage) -> bool {
        let pixels = &image.pixels;
        let first = pixels[0];
        pixels.iter().all(|p| *p == first)
    }

    fn has_nan_artifacts(&self, image: &ColorImage) -> bool {
        // NaN often renders as solid black or white clusters
        // Look for suspicious uniform rectangles
        // ...
        false
    }
}

#[derive(Debug)]
pub enum VisualIssue {
    BlankFrame,
    PossibleNanArtifacts,
    EdgeClipping,
    PossibleZFighting,
    UnexpectedColor { expected: Color, found: Color },
}
```

### CI Integration

```yaml
# .github/workflows/visual-tests.yml
name: Visual Tests

on:
  push:
    paths:
      - 'crates/clifford-viz/**'
  pull_request:
    paths:
      - 'crates/clifford-viz/**'

jobs:
  visual-tests:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - name: Install dependencies
        run: |
          sudo apt-get update
          sudo apt-get install -y xvfb libxkbcommon-x11-0 libgl1-mesa-dri

      - name: Run coordinate and invariant tests
        run: |
          cargo test -p clifford-viz --lib

      - name: Run visual inspection
        run: |
          xvfb-run --auto-servernum cargo test -p clifford-viz --test visual_inspection
        continue-on-error: true  # Flag for review, don't block

      - name: Upload inspection results
        if: always()
        uses: actions/upload-artifact@v4
        with:
          name: visual-inspection-results
          path: target/visual-inspection/
          retention-days: 7
```

## Test Infrastructure

### Directory Structure

```
crates/clifford-viz/
  src/
    testing/
      mod.rs              # Test utilities
      scene.rs            # Scene graph for assertions
      inspector.rs        # Automated visual inspection
      viewport.rs         # Coordinate mapping tests
  tests/
    coordinate_tests.rs   # World-to-screen mapping
    invariant_tests.rs    # Property-based visual invariants
    scene_tests.rs        # Scene graph assertions
    visual_inspection.rs  # Automated screenshot analysis
```

### Scene Graph for Testing

```rust
/// Testable scene representation
pub struct TestScene {
    vectors: HashMap<String, TestVector>,
    points: Vec<TestPoint>,
    lines: Vec<TestLine>,
    arcs: HashMap<String, TestArc>,
}

impl TestScene {
    pub fn has_vector(&self, name: &str) -> bool {
        self.vectors.contains_key(name)
    }

    pub fn get_vector(&self, name: &str) -> &TestVector {
        &self.vectors[name]
    }

    pub fn point_count(&self) -> usize {
        self.points.len()
    }

    pub fn line_count(&self) -> usize {
        self.lines.len()
    }
}

pub struct TestVector {
    pub start: Point2,
    pub end: Point2,
    pub color: Color,
}

impl TestVector {
    pub fn angle_to(&self, other: &TestVector) -> f32 {
        let self_dir = (self.end - self.start).normalize();
        let other_dir = (other.end - other.start).normalize();
        self_dir.angle_between(other_dir)
    }
}
```

## Implementation Tasks

1. [ ] Create `crates/clifford-viz/src/testing/` module
2. [ ] Implement `TestScene` and test primitives
3. [ ] Add `build_scene()` method to each demo
4. [ ] Write coordinate mapping tests for each algebra
5. [ ] Write property-based visual invariant tests
6. [ ] Implement `VisualInspector` for automated analysis
7. [ ] Add GitHub Actions workflow
8. [ ] Document test patterns in CONTRIBUTING.md

## Verification

```bash
# Run all visualization tests
cargo test -p clifford-viz

# Run specific test categories
cargo test -p clifford-viz coordinate
cargo test -p clifford-viz invariant
cargo test -p clifford-viz scene

# Run visual inspection (requires display)
xvfb-run cargo test -p clifford-viz --test visual_inspection
```

## Benefits Over Golden Images

| Aspect | Golden Images | Invariant Tests |
|--------|---------------|-----------------|
| Detects bugs | Sometimes | Yes |
| Proves correctness | No | Yes |
| Handles rendering differences | Poorly | N/A |
| Maintenance burden | High | Low |
| Documents intent | No | Yes |
| Works headless | Awkwardly | Mostly |
| Test failures are actionable | "Something changed" | "This invariant violated" |

## When Human Review is Needed

The automated inspection flags issues but doesn't block CI. When issues are flagged:

1. CI uploads screenshots to artifacts
2. Reviewer checks flagged issues
3. If issue is real: fix the bug
4. If false positive: improve inspector heuristics

This gives automated first-pass analysis while keeping humans in the loop for subjective visual quality.
