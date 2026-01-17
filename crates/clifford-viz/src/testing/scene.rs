//! Scene graph representation for testing.
//!
//! Provides a testable representation of what a visualization renders,
//! allowing assertions on the structure and properties of rendered primitives
//! without needing to inspect actual pixels.

use std::collections::HashMap;

use egui::Color32;

/// A 2D point for testing purposes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point2 {
    /// X coordinate.
    pub x: f64,
    /// Y coordinate.
    pub y: f64,
}

impl Point2 {
    /// Create a new point.
    #[must_use]
    pub const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    /// The origin point (0, 0).
    pub const ORIGIN: Self = Self { x: 0.0, y: 0.0 };

    /// Distance from another point.
    #[must_use]
    pub fn distance_to(&self, other: &Self) -> f64 {
        ((self.x - other.x).powi(2) + (self.y - other.y).powi(2)).sqrt()
    }

    /// Distance from origin.
    #[must_use]
    pub fn magnitude(&self) -> f64 {
        (self.x * self.x + self.y * self.y).sqrt()
    }

    /// Angle from positive x-axis.
    #[must_use]
    pub fn angle(&self) -> f64 {
        self.y.atan2(self.x)
    }
}

/// A vector in the test scene (with start and end points).
#[derive(Debug, Clone)]
pub struct TestVector {
    /// Starting point of the vector.
    pub start: Point2,
    /// Ending point of the vector.
    pub end: Point2,
    /// Color of the vector.
    pub color: Color32,
    /// Optional name/label.
    pub name: Option<String>,
}

impl TestVector {
    /// Create a new test vector.
    #[must_use]
    pub fn new(start: Point2, end: Point2, color: Color32) -> Self {
        Self {
            start,
            end,
            color,
            name: None,
        }
    }

    /// Create a test vector with a name.
    #[must_use]
    pub fn with_name(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    /// Get the direction vector (end - start).
    #[must_use]
    pub fn direction(&self) -> Point2 {
        Point2::new(self.end.x - self.start.x, self.end.y - self.start.y)
    }

    /// Get the length of the vector.
    #[must_use]
    pub fn length(&self) -> f64 {
        self.direction().magnitude()
    }

    /// Get the angle of the vector direction.
    #[must_use]
    pub fn angle(&self) -> f64 {
        self.direction().angle()
    }

    /// Compute the angle to another vector.
    #[must_use]
    pub fn angle_to(&self, other: &Self) -> f64 {
        let d1 = self.direction();
        let d2 = other.direction();
        let dot = d1.x * d2.x + d1.y * d2.y;
        let mag1 = d1.magnitude();
        let mag2 = d2.magnitude();
        if mag1 < 1e-10 || mag2 < 1e-10 {
            return 0.0;
        }
        (dot / (mag1 * mag2)).clamp(-1.0, 1.0).acos()
    }
}

/// A point marker in the test scene.
#[derive(Debug, Clone)]
pub struct TestPoint {
    /// Position of the point.
    pub position: Point2,
    /// Color of the point.
    pub color: Color32,
    /// Radius of the marker.
    pub radius: f32,
    /// Optional name/label.
    pub name: Option<String>,
}

impl TestPoint {
    /// Create a new test point.
    #[must_use]
    pub fn new(x: f64, y: f64, color: Color32, radius: f32) -> Self {
        Self {
            position: Point2::new(x, y),
            color,
            radius,
            name: None,
        }
    }

    /// Create a test point with a name.
    #[must_use]
    pub fn with_name(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }
}

/// An arc in the test scene.
#[derive(Debug, Clone)]
pub struct TestArc {
    /// Center of the arc.
    pub center: Point2,
    /// Radius of the arc.
    pub radius: f64,
    /// Starting angle in radians.
    pub start_angle: f64,
    /// Ending angle in radians.
    pub end_angle: f64,
    /// Color of the arc.
    pub color: Color32,
    /// Optional name/label.
    pub name: Option<String>,
}

impl TestArc {
    /// Create a new test arc.
    #[must_use]
    pub fn new(
        center: Point2,
        radius: f64,
        start_angle: f64,
        end_angle: f64,
        color: Color32,
    ) -> Self {
        Self {
            center,
            radius,
            start_angle,
            end_angle,
            color,
            name: None,
        }
    }

    /// Get the arc angle (sweep).
    #[must_use]
    pub fn sweep(&self) -> f64 {
        self.end_angle - self.start_angle
    }
}

/// A line in the test scene.
#[derive(Debug, Clone)]
pub struct TestLine {
    /// A point on the line.
    pub point: Point2,
    /// Direction of the line.
    pub direction: Point2,
    /// Color of the line.
    pub color: Color32,
    /// Optional name/label.
    pub name: Option<String>,
}

impl TestLine {
    /// Check if a point lies on this line (within tolerance).
    #[must_use]
    pub fn contains_point(&self, x: f64, y: f64, epsilon: f64) -> bool {
        // Distance from point to line
        let dx = x - self.point.x;
        let dy = y - self.point.y;
        let dir_mag = self.direction.magnitude();
        if dir_mag < epsilon {
            // Degenerate line
            return (dx * dx + dy * dy).sqrt() < epsilon;
        }
        // Cross product magnitude gives distance * direction magnitude
        let cross = dx * self.direction.y - dy * self.direction.x;
        (cross / dir_mag).abs() < epsilon
    }
}

/// A testable scene representation.
///
/// Collects the primitives that would be rendered, allowing assertions
/// on scene structure without rendering to pixels.
#[derive(Debug, Default)]
pub struct TestScene {
    /// Named vectors in the scene.
    pub vectors: HashMap<String, TestVector>,
    /// All vectors (including unnamed).
    pub all_vectors: Vec<TestVector>,
    /// Points in the scene.
    pub points: Vec<TestPoint>,
    /// Named arcs in the scene.
    pub arcs: HashMap<String, TestArc>,
    /// All arcs (including unnamed).
    pub all_arcs: Vec<TestArc>,
    /// Lines in the scene.
    pub lines: Vec<TestLine>,
}

impl TestScene {
    /// Create an empty test scene.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a vector to the scene.
    pub fn add_vector(&mut self, vector: TestVector) {
        if let Some(ref name) = vector.name {
            self.vectors.insert(name.clone(), vector.clone());
        }
        self.all_vectors.push(vector);
    }

    /// Add a point to the scene.
    pub fn add_point(&mut self, point: TestPoint) {
        self.points.push(point);
    }

    /// Add an arc to the scene.
    pub fn add_arc(&mut self, arc: TestArc) {
        if let Some(ref name) = arc.name {
            self.arcs.insert(name.clone(), arc.clone());
        }
        self.all_arcs.push(arc);
    }

    /// Add a line to the scene.
    pub fn add_line(&mut self, line: TestLine) {
        self.lines.push(line);
    }

    /// Check if a named vector exists.
    #[must_use]
    pub fn has_vector(&self, name: &str) -> bool {
        self.vectors.contains_key(name)
    }

    /// Get a named vector.
    #[must_use]
    pub fn get_vector(&self, name: &str) -> Option<&TestVector> {
        self.vectors.get(name)
    }

    /// Check if a named arc exists.
    #[must_use]
    pub fn has_arc(&self, name: &str) -> bool {
        self.arcs.contains_key(name)
    }

    /// Get a named arc.
    #[must_use]
    pub fn get_arc(&self, name: &str) -> Option<&TestArc> {
        self.arcs.get(name)
    }

    /// Get the count of points.
    #[must_use]
    pub fn point_count(&self) -> usize {
        self.points.len()
    }

    /// Get the count of vectors.
    #[must_use]
    pub fn vector_count(&self) -> usize {
        self.all_vectors.len()
    }

    /// Get the count of lines.
    #[must_use]
    pub fn line_count(&self) -> usize {
        self.lines.len()
    }

    /// Get the count of arcs.
    #[must_use]
    pub fn arc_count(&self) -> usize {
        self.all_arcs.len()
    }
}
