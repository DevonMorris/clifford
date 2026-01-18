//! Automated visual inspection utilities.
//!
//! This module provides tools for automated first-pass analysis of rendered frames.
//! It flags potential issues without blocking CI, allowing humans to review
//! flagged issues while catching obvious problems automatically.
//!
//! # Philosophy
//!
//! The inspector doesn't determine correctness (that's what invariant tests do).
//! Instead, it flags potential rendering issues that warrant human review:
//!
//! - **Blank frames**: Nothing rendered (possible rendering failure)
//! - **NaN artifacts**: Suspicious pixel patterns from invalid math
//! - **Clipping**: Content cut off at viewport edges
//! - **Uniform regions**: Large areas of identical color (possible fill errors)
//!
//! # Usage
//!
//! ```ignore
//! use clifford_viz::testing::VisualInspector;
//!
//! let inspector = VisualInspector::default();
//! let issues = inspector.inspect(&image);
//!
//! if !issues.is_empty() {
//!     println!("Visual issues found: {:?}", issues);
//! }
//! ```

use egui::ColorImage;

/// Region bounds: (x, y, width, height).
pub type RegionBounds = (usize, usize, usize, usize);

/// A potential visual issue detected by automated inspection.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VisualIssue {
    /// The frame appears to be blank (all pixels the same color).
    BlankFrame,

    /// Possible NaN artifacts detected (suspicious uniform rectangular regions).
    PossibleNanArtifacts {
        /// Location of the suspicious region (x, y, width, height).
        region: RegionBounds,
    },

    /// Content appears to be clipped at the viewport edge.
    EdgeClipping {
        /// Which edge has the clipping (top, bottom, left, right).
        edge: Edge,
    },

    /// Large uniform color region that may indicate a fill error.
    UnexpectedUniformRegion {
        /// Location and size of the region.
        region: RegionBounds,
        /// The uniform color.
        color: [u8; 4],
    },

    /// Very low color diversity (image may not be rendering properly).
    LowColorDiversity {
        /// Number of unique colors found.
        unique_colors: usize,
    },
}

/// Viewport edge for clipping detection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Edge {
    /// Top edge.
    Top,
    /// Bottom edge.
    Bottom,
    /// Left edge.
    Left,
    /// Right edge.
    Right,
}

/// Automated visual inspection tool.
///
/// Analyzes rendered frames for common issues that may indicate rendering problems.
/// This is a first-pass filter to flag potential issues for human review.
#[derive(Debug, Clone)]
pub struct VisualInspector {
    /// Minimum number of unique colors expected in a valid frame.
    pub min_color_diversity: usize,

    /// Threshold for considering a region "uniform" (pixels must match exactly).
    pub uniform_region_min_size: usize,

    /// Percentage of edge pixels that must be non-background to flag clipping.
    pub edge_clipping_threshold: f32,

    /// Background color to ignore when checking for clipping.
    pub background_color: [u8; 4],
}

impl Default for VisualInspector {
    fn default() -> Self {
        Self {
            min_color_diversity: 5,
            uniform_region_min_size: 100,
            edge_clipping_threshold: 0.1,
            background_color: [255, 255, 255, 255], // White background
        }
    }
}

impl VisualInspector {
    /// Create a new visual inspector with default settings.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the minimum color diversity threshold.
    #[must_use]
    pub fn with_min_color_diversity(mut self, min: usize) -> Self {
        self.min_color_diversity = min;
        self
    }

    /// Set the background color for clipping detection.
    #[must_use]
    pub fn with_background_color(mut self, color: [u8; 4]) -> Self {
        self.background_color = color;
        self
    }

    /// Analyze an image for potential visual issues.
    ///
    /// Returns a list of detected issues. An empty list means no obvious problems
    /// were found (but doesn't guarantee correctness).
    #[must_use]
    pub fn inspect(&self, image: &ColorImage) -> Vec<VisualIssue> {
        let mut issues = Vec::new();

        // Check for blank frame
        if self.is_blank(image) {
            issues.push(VisualIssue::BlankFrame);
            return issues; // No point checking further if blank
        }

        // Check color diversity
        let unique_colors = self.count_unique_colors(image);
        if unique_colors < self.min_color_diversity {
            issues.push(VisualIssue::LowColorDiversity { unique_colors });
        }

        // Check for edge clipping
        for edge in [Edge::Top, Edge::Bottom, Edge::Left, Edge::Right] {
            if self.has_edge_clipping(image, edge) {
                issues.push(VisualIssue::EdgeClipping { edge });
            }
        }

        // Check for suspicious uniform regions (potential NaN artifacts)
        if let Some((region, color)) = self.find_suspicious_uniform_region(image) {
            // Solid black or white rectangular regions may be NaN artifacts
            if color == [0, 0, 0, 255] || color == [255, 255, 255, 255] {
                issues.push(VisualIssue::PossibleNanArtifacts { region });
            } else {
                issues.push(VisualIssue::UnexpectedUniformRegion { region, color });
            }
        }

        issues
    }

    /// Check if the image is completely blank (all pixels identical).
    #[must_use]
    pub fn is_blank(&self, image: &ColorImage) -> bool {
        if image.pixels.is_empty() {
            return true;
        }

        let first = image.pixels[0];
        image.pixels.iter().all(|p| *p == first)
    }

    /// Count the number of unique colors in the image.
    ///
    /// For performance, this samples the image rather than checking every pixel
    /// on large images.
    #[must_use]
    pub fn count_unique_colors(&self, image: &ColorImage) -> usize {
        use std::collections::HashSet;

        let mut colors = HashSet::new();
        let total_pixels = image.pixels.len();

        // Sample at most 10000 pixels for performance
        let sample_step = (total_pixels / 10000).max(1);

        for (i, pixel) in image.pixels.iter().enumerate() {
            if i % sample_step == 0 {
                colors.insert([pixel.r(), pixel.g(), pixel.b(), pixel.a()]);
            }

            // Early exit if we've found enough diversity
            if colors.len() >= self.min_color_diversity * 2 {
                break;
            }
        }

        colors.len()
    }

    /// Check if content appears to be clipped at the given edge.
    ///
    /// This heuristic looks for non-background pixels at the very edge of the viewport,
    /// which may indicate content is being cut off.
    #[must_use]
    pub fn has_edge_clipping(&self, image: &ColorImage, edge: Edge) -> bool {
        let [width, height] = image.size;
        if width == 0 || height == 0 {
            return false;
        }

        let edge_pixels: Vec<_> = match edge {
            Edge::Top => (0..width).map(|x| image.pixels[x]).collect(),
            Edge::Bottom => {
                let start = (height - 1) * width;
                (0..width).map(|x| image.pixels[start + x]).collect()
            }
            Edge::Left => (0..height).map(|y| image.pixels[y * width]).collect(),
            Edge::Right => (0..height)
                .map(|y| image.pixels[y * width + width - 1])
                .collect(),
        };

        let bg = egui::Color32::from_rgba_unmultiplied(
            self.background_color[0],
            self.background_color[1],
            self.background_color[2],
            self.background_color[3],
        );

        let non_bg_count = edge_pixels.iter().filter(|&&p| p != bg).count();
        let non_bg_ratio = non_bg_count as f32 / edge_pixels.len() as f32;

        non_bg_ratio > self.edge_clipping_threshold
    }

    /// Find a suspicious uniform rectangular region in the image.
    ///
    /// Returns the region bounds and color if found.
    fn find_suspicious_uniform_region(
        &self,
        image: &ColorImage,
    ) -> Option<(RegionBounds, [u8; 4])> {
        let [width, height] = image.size;
        if width < 10 || height < 10 {
            return None;
        }

        // Simple heuristic: scan for runs of identical pixels
        // This is a simplified version - a production implementation would use
        // more sophisticated region detection

        let min_run_length = self.uniform_region_min_size.min(width / 4);

        for y in (0..height).step_by(height / 10 + 1) {
            let mut run_start = 0;
            let mut run_color = image.pixels[y * width];

            for x in 1..width {
                let pixel = image.pixels[y * width + x];
                if pixel == run_color {
                    continue;
                }

                // End of run
                let run_length = x - run_start;
                if run_length >= min_run_length {
                    let color = [run_color.r(), run_color.g(), run_color.b(), run_color.a()];
                    // Skip if it's the background color
                    if color != self.background_color {
                        return Some(((run_start, y, run_length, 1), color));
                    }
                }

                run_start = x;
                run_color = pixel;
            }

            // Check final run
            let run_length = width - run_start;
            if run_length >= min_run_length {
                let color = [run_color.r(), run_color.g(), run_color.b(), run_color.a()];
                if color != self.background_color {
                    return Some(((run_start, y, run_length, 1), color));
                }
            }
        }

        None
    }
}

/// Result of a visual inspection run.
#[derive(Debug, Clone)]
pub struct InspectionResult {
    /// The name/identifier of the inspected frame.
    pub frame_name: String,

    /// List of issues found.
    pub issues: Vec<VisualIssue>,

    /// Whether this result should be flagged for human review.
    pub needs_review: bool,
}

impl InspectionResult {
    /// Create a new inspection result.
    #[must_use]
    pub fn new(frame_name: impl Into<String>, issues: Vec<VisualIssue>) -> Self {
        let needs_review = !issues.is_empty();
        Self {
            frame_name: frame_name.into(),
            issues,
            needs_review,
        }
    }

    /// Check if any issues were found.
    #[must_use]
    pub fn has_issues(&self) -> bool {
        !self.issues.is_empty()
    }

    /// Get a human-readable summary of the issues.
    #[must_use]
    pub fn summary(&self) -> String {
        if self.issues.is_empty() {
            return format!("{}: OK", self.frame_name);
        }

        let issue_strs: Vec<_> = self
            .issues
            .iter()
            .map(|issue| match issue {
                VisualIssue::BlankFrame => "blank frame".to_string(),
                VisualIssue::PossibleNanArtifacts { region } => {
                    format!("possible NaN artifacts at {:?}", region)
                }
                VisualIssue::EdgeClipping { edge } => format!("clipping at {:?} edge", edge),
                VisualIssue::UnexpectedUniformRegion { region, .. } => {
                    format!("uniform region at {:?}", region)
                }
                VisualIssue::LowColorDiversity { unique_colors } => {
                    format!("low color diversity ({} colors)", unique_colors)
                }
            })
            .collect();

        format!("{}: {}", self.frame_name, issue_strs.join(", "))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use egui::Color32;

    fn create_test_image(width: usize, height: usize, color: Color32) -> ColorImage {
        ColorImage {
            size: [width, height],
            pixels: vec![color; width * height],
        }
    }

    #[test]
    fn detects_blank_frame() {
        let inspector = VisualInspector::new();
        let image = create_test_image(100, 100, Color32::WHITE);

        let issues = inspector.inspect(&image);

        assert!(issues.contains(&VisualIssue::BlankFrame));
    }

    #[test]
    fn non_blank_frame_not_flagged_as_blank() {
        let inspector = VisualInspector::new();
        let mut image = create_test_image(100, 100, Color32::WHITE);

        // Add some non-white pixels
        for i in 0..50 {
            image.pixels[i] = Color32::RED;
        }
        for i in 50..100 {
            image.pixels[i] = Color32::BLUE;
        }
        for i in 100..150 {
            image.pixels[i] = Color32::GREEN;
        }
        for i in 150..200 {
            image.pixels[i] = Color32::BLACK;
        }
        for i in 200..250 {
            image.pixels[i] = Color32::YELLOW;
        }

        let issues = inspector.inspect(&image);

        assert!(!issues.contains(&VisualIssue::BlankFrame));
    }

    #[test]
    fn detects_low_color_diversity() {
        let inspector = VisualInspector::new().with_min_color_diversity(5);
        let mut image = create_test_image(100, 100, Color32::WHITE);

        // Add only 3 colors (below threshold of 5)
        for i in 0..3000 {
            image.pixels[i] = Color32::RED;
        }
        for i in 3000..6000 {
            image.pixels[i] = Color32::BLUE;
        }

        let issues = inspector.inspect(&image);

        assert!(
            issues
                .iter()
                .any(|i| matches!(i, VisualIssue::LowColorDiversity { .. })),
            "Should detect low color diversity"
        );
    }

    #[test]
    fn count_unique_colors_works() {
        let inspector = VisualInspector::new();
        let mut image = create_test_image(100, 100, Color32::WHITE);

        // Add 5 distinct colors
        for i in 0..2000 {
            image.pixels[i] = Color32::RED;
        }
        for i in 2000..4000 {
            image.pixels[i] = Color32::GREEN;
        }
        for i in 4000..6000 {
            image.pixels[i] = Color32::BLUE;
        }
        for i in 6000..8000 {
            image.pixels[i] = Color32::BLACK;
        }
        // Remaining pixels are white (5th color)

        let count = inspector.count_unique_colors(&image);
        assert!(
            count >= 5,
            "Should find at least 5 unique colors, got {}",
            count
        );
    }

    #[test]
    fn inspection_result_summary() {
        let result = InspectionResult::new(
            "test_frame",
            vec![
                VisualIssue::BlankFrame,
                VisualIssue::LowColorDiversity { unique_colors: 2 },
            ],
        );

        let summary = result.summary();
        assert!(summary.contains("test_frame"));
        assert!(summary.contains("blank frame"));
        assert!(summary.contains("low color diversity"));
    }

    #[test]
    fn empty_inspection_result() {
        let result = InspectionResult::new("good_frame", vec![]);

        assert!(!result.has_issues());
        assert!(!result.needs_review);
        assert!(result.summary().contains("OK"));
    }
}
