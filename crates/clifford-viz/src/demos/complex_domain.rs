//! Complex Domain Coloring Visualization
//!
//! This demo visualizes complex functions using domain coloring, where:
//! - Hue encodes the argument (phase angle) of f(z)
//! - Brightness encodes the magnitude |f(z)|
//!
//! This visualization technique reveals the structure of complex functions:
//! zeros appear as points where all colors meet, and poles show infinite
//! color cycling.

use crate::common::prelude::*;
use clifford::specialized::complex::Complex;
use egui::Color32;
use egui_plot::{Plot, PlotImage, PlotPoint};

/// Fixed viewport bounds for stable viewing experience.
const DEFAULT_VIEW_SCALE: f64 = 3.0;

/// Image resolution for domain coloring (pixels per side).
const IMAGE_RESOLUTION: usize = 256;

/// Available complex functions for visualization.
#[derive(Clone, Copy, PartialEq, Eq, Default)]
enum ComplexFunction {
    /// f(z) = z²
    #[default]
    Square,
    /// f(z) = z³
    Cube,
    /// f(z) = 1/z
    Reciprocal,
    /// f(z) = eᶻ
    Exp,
    /// f(z) = sin(z)
    Sin,
    /// f(z) = cos(z)
    Cos,
    /// f(z) = log(z)
    Log,
    /// f(z) = √z
    Sqrt,
    /// f(z) = z (identity)
    Identity,
}

impl ComplexFunction {
    /// Returns the display name of the function.
    fn name(&self) -> &'static str {
        match self {
            Self::Identity => "z",
            Self::Square => "z\u{00b2}",
            Self::Cube => "z\u{00b3}",
            Self::Reciprocal => "1/z",
            Self::Exp => "e\u{1d63}",
            Self::Sin => "sin(z)",
            Self::Cos => "cos(z)",
            Self::Log => "log(z)",
            Self::Sqrt => "\u{221a}z",
        }
    }

    /// Evaluates the function at z.
    fn evaluate(&self, z: Complex<f64>) -> Complex<f64> {
        match self {
            Self::Identity => z,
            Self::Square => z * z,
            Self::Cube => z * z * z,
            Self::Reciprocal => Complex::one() / z,
            Self::Exp => z.exp(),
            Self::Sin => z.sin(),
            Self::Cos => z.cos(),
            Self::Log => z.ln(),
            Self::Sqrt => z.sqrt(),
        }
    }

    /// Returns all available functions.
    fn all() -> &'static [ComplexFunction] {
        &[
            Self::Identity,
            Self::Square,
            Self::Cube,
            Self::Reciprocal,
            Self::Exp,
            Self::Sin,
            Self::Cos,
            Self::Log,
            Self::Sqrt,
        ]
    }
}

/// Domain coloring visualization for complex functions.
pub struct ComplexDomainDemo {
    /// Current function selection.
    function: ComplexFunction,
    /// View center in complex plane.
    view_center: (f64, f64),
    /// View scale (half-width of visible region).
    view_scale: f64,
    /// Whether to show contour lines.
    show_contours: bool,
    /// Whether to show coordinate grid.
    show_grid: bool,
    /// Cached image texture.
    cached_image: Option<egui::TextureHandle>,
    /// Whether cache needs refresh.
    needs_refresh: bool,
    /// Hover position in complex plane.
    hover_z: Option<Complex<f64>>,
    /// Brightness scaling factor.
    brightness_scale: f32,
}

impl Default for ComplexDomainDemo {
    fn default() -> Self {
        Self {
            function: ComplexFunction::default(),
            view_center: (0.0, 0.0),
            view_scale: DEFAULT_VIEW_SCALE,
            show_contours: true,
            show_grid: true,
            cached_image: None,
            needs_refresh: true,
            hover_z: None,
            brightness_scale: 0.3,
        }
    }
}

impl ComplexDomainDemo {
    /// Computes the domain coloring for the current function and view.
    fn compute_domain_coloring(&self, ctx: &egui::Context) -> egui::TextureHandle {
        let mut pixels = vec![Color32::BLACK; IMAGE_RESOLUTION * IMAGE_RESOLUTION];

        let center_x = self.view_center.0;
        let center_y = self.view_center.1;
        let scale = self.view_scale;

        for py in 0..IMAGE_RESOLUTION {
            for px in 0..IMAGE_RESOLUTION {
                // Map pixel to complex plane
                let x = center_x + scale * (2.0 * (px as f64) / (IMAGE_RESOLUTION as f64) - 1.0);
                let y = center_y + scale * (1.0 - 2.0 * (py as f64) / (IMAGE_RESOLUTION as f64));

                let z = Complex::new(x, y);
                let fz = self.function.evaluate(z);

                let color = self.domain_color(fz);
                pixels[py * IMAGE_RESOLUTION + px] = color;
            }
        }

        let image = egui::ColorImage {
            size: [IMAGE_RESOLUTION, IMAGE_RESOLUTION],
            pixels,
        };

        ctx.load_texture("domain_coloring", image, egui::TextureOptions::LINEAR)
    }

    /// Computes the color for a complex value using domain coloring.
    fn domain_color(&self, z: Complex<f64>) -> Color32 {
        let arg = z.arg();
        let mag = z.norm();

        // Hue from argument: 0 = red, pi/2 = yellow, pi = cyan, 3pi/2 = blue
        let hue = (arg / std::f64::consts::TAU + 1.0) % 1.0;

        // Lightness from magnitude with logarithmic scaling
        let log_mag = mag.ln();
        let brightness_scale = f64::from(self.brightness_scale);
        let lightness = if self.show_contours {
            // Add contour lines
            let contour = (log_mag * 2.0).fract().abs();
            let base_lightness =
                0.5 + 0.4 * (1.0 / (1.0 + (-log_mag * brightness_scale).exp()) - 0.5);
            let contour_factor = if !(0.05..=0.95).contains(&contour) {
                0.7
            } else {
                1.0
            };
            base_lightness * contour_factor
        } else {
            0.5 + 0.4 * (1.0 / (1.0 + (-log_mag * brightness_scale).exp()) - 0.5)
        };

        // Handle special cases
        if !mag.is_finite() || mag < 1e-10 {
            return Color32::BLACK;
        }
        if mag > 1e10 {
            return Color32::WHITE;
        }

        hsl_to_color(hue as f32, 1.0, lightness.clamp(0.1, 0.9) as f32)
    }
}

impl VisualizationApp for ComplexDomainDemo {
    fn name(&self) -> &'static str {
        "Complex Domain Coloring"
    }

    fn update(&mut self, _dt: f32) {
        // No animation in this demo
    }

    fn render(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();

        // Refresh texture if needed
        if self.needs_refresh || self.cached_image.is_none() {
            self.cached_image = Some(self.compute_domain_coloring(&ctx));
            self.needs_refresh = false;
        }

        let scale = self.view_scale;
        let center_x = self.view_center.0;
        let center_y = self.view_center.1;

        Plot::new("complex_domain_plot")
            .data_aspect(1.0)
            .show_axes(false)
            .show_grid(false)
            .allow_zoom(true)
            .allow_drag(true)
            .auto_bounds(egui::Vec2b::new(false, false))
            .include_x(center_x - scale)
            .include_x(center_x + scale)
            .include_y(center_y - scale)
            .include_y(center_y + scale)
            .show(ui, |plot_ui| {
                // Draw the domain coloring image
                if let Some(texture) = &self.cached_image {
                    let image = PlotImage::new(
                        texture,
                        PlotPoint::new(center_x, center_y),
                        [(2.0 * scale) as f32, (2.0 * scale) as f32],
                    );
                    plot_ui.image(image);
                }

                // Draw coordinate grid on top of the image
                if self.show_grid {
                    for l in grid_2d(&ctx, scale, 1.0) {
                        plot_ui.line(l);
                    }
                    for axis in axes_2d(&ctx, scale) {
                        plot_ui.line(axis);
                    }
                }

                // Track hover position
                if let Some(pos) = plot_ui.pointer_coordinate() {
                    self.hover_z = Some(Complex::new(pos.x, pos.y));
                } else {
                    self.hover_z = None;
                }
            });
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        group_header(ui, "Function");

        // Function selection grid
        egui::Grid::new("function_grid")
            .num_columns(3)
            .spacing([8.0, 4.0])
            .show(ui, |ui| {
                for (i, func) in ComplexFunction::all().iter().enumerate() {
                    if ui
                        .selectable_label(self.function == *func, func.name())
                        .clicked()
                    {
                        self.function = *func;
                        self.needs_refresh = true;
                    }
                    if (i + 1) % 3 == 0 {
                        ui.end_row();
                    }
                }
            });

        section_separator(ui, Some("View"));

        // View controls
        ui.horizontal(|ui| {
            ui.label("Scale:");
            let mut log_scale = self.view_scale.ln();
            if ui
                .add(egui::Slider::new(&mut log_scale, (-2.0)..=3.0).logarithmic(false))
                .changed()
            {
                self.view_scale = log_scale.exp();
                self.needs_refresh = true;
            }
        });

        ui.horizontal(|ui| {
            ui.label("Brightness:");
            if ui
                .add(egui::Slider::new(&mut self.brightness_scale, 0.1..=1.0))
                .changed()
            {
                self.needs_refresh = true;
            }
        });

        // Reset view button
        if ui.button("Reset View").clicked() {
            self.view_center = (0.0, 0.0);
            self.view_scale = DEFAULT_VIEW_SCALE;
            self.needs_refresh = true;
        }

        section_separator(ui, Some("Display"));

        ui.horizontal(|ui| {
            if ui.checkbox(&mut self.show_grid, "Grid").changed() {
                // Grid doesn't need texture refresh
            }
            if ui.checkbox(&mut self.show_contours, "Contours").changed() {
                self.needs_refresh = true;
            }
        });

        // Hover info
        section_separator(ui, Some("Hover Info"));

        if let Some(z) = self.hover_z {
            let fz = self.function.evaluate(z);
            value_display(ui, "z", z.real() as f32, 3);
            value_display(ui, "+ i", z.imag() as f32, 3);
            ui.add_space(spacing::XS);
            value_display(ui, "f(z)", fz.real() as f32, 3);
            value_display(ui, "+ i", fz.imag() as f32, 3);
            ui.add_space(spacing::XS);
            value_display(ui, "|f(z)|", fz.norm() as f32, 3);
            value_display(ui, "arg(f(z))", fz.arg() as f32, 3);
        } else {
            ui.label("Hover over plot to see values");
        }

        // Color legend
        section_separator(ui, Some("Color Legend"));
        info_box(
            ui,
            "Hue = argument (phase)\n\
             0\u{00b0} = red (positive real)\n\
             90\u{00b0} = yellow (positive imag)\n\
             180\u{00b0} = cyan (negative real)\n\
             270\u{00b0} = blue (negative imag)\n\n\
             Brightness = magnitude\n\
             Dark = small |f(z)|\n\
             Bright = large |f(z)|",
        );
    }

    fn info(&self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.label(format!("f(z) = {}", self.function.name()));
            ui.separator();
            ui.label(format!(
                "View: [{:.1}, {:.1}]",
                self.view_center.0 - self.view_scale,
                self.view_center.0 + self.view_scale
            ));
            ui.separator();
            if let Some(z) = self.hover_z {
                ui.label(format!("z = {:.2} + {:.2}i", z.real(), z.imag()));
            }
        });
    }

    fn educational_content(&self) -> Option<EducationalContent> {
        Some(COMPLEX_DOMAIN_EDUCATION)
    }
}

/// Educational content for the complex domain coloring demo.
const COMPLEX_DOMAIN_EDUCATION: EducationalContent = EducationalContent {
    title: "Domain Coloring for Complex Functions",

    overview: "\
Domain coloring is a technique for visualizing complex functions f: C \u{2192} C.
Since both input and output are 2-dimensional, we can't use a simple 2D graph.
Instead, we color each point z in the domain based on the value f(z).

The COLOR (hue) shows the ARGUMENT (phase angle) of f(z), while the
BRIGHTNESS shows the MAGNITUDE |f(z)|. This reveals beautiful patterns
that expose the structure of complex functions.",

    math_background: "\
Every complex number can be written in polar form:
    z = r\u{00b7}e^(i\u{03b8}) = r(cos\u{03b8} + i\u{00b7}sin\u{03b8})

where r = |z| is the magnitude and \u{03b8} = arg(z) is the argument.

Domain coloring maps:
  \u{2022} arg(f(z)) \u{2192} hue (color wheel)
  \u{2022} |f(z)| \u{2192} brightness (log scale)

Key features to recognize:
  \u{2022} ZEROS: All colors meet at a point (f(z) = 0)
  \u{2022} POLES: Colors cycle infinitely (f(z) \u{2192} \u{221e})
  \u{2022} BRANCH CUTS: Color discontinuities (e.g., log, sqrt)
  \u{2022} Essential singularities: Wild color chaos

For z\u{00b2}, colors cycle TWICE around the origin (degree 2).
For z\u{00b3}, colors cycle THREE times (degree 3).",

    how_to_use: "\
\u{2022} Select a function from the grid
\u{2022} Hover to see z and f(z) values
\u{2022} Adjust scale to zoom in/out
\u{2022} Toggle contours to see magnitude levels
\u{2022} Look for zeros (all colors meet) and poles (color cycles)",

    key_concepts: "\
\u{2022} Complex functions map C \u{2192} C (4D problem!)
\u{2022} Hue = argument, Brightness = magnitude
\u{2022} Zeros appear as points where all colors meet
\u{2022} Poles show infinite color cycling
\u{2022} The number of color cycles = degree of zero/pole
\u{2022} Branch cuts create color discontinuities",

    resources: &[
        (
            "Domain Coloring - Wikipedia",
            "https://en.wikipedia.org/wiki/Domain_coloring",
        ),
        (
            "Visual Complex Analysis",
            "https://www.visual-complex-analysis.com/",
        ),
    ],
};
