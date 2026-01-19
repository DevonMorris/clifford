//! Mandelbrot and Julia Set Visualization
//!
//! This demo explores the Mandelbrot set and its corresponding Julia sets,
//! demonstrating the relationship between these famous fractals through
//! complex number iteration.

use crate::common::prelude::*;
use clifford::specialized::complex::Complex;
use egui::Color32;
use egui_plot::{Plot, PlotImage, PlotPoint, Points};

/// Image resolution for fractal rendering.
const IMAGE_RESOLUTION: usize = 256;

/// Default maximum iterations.
const DEFAULT_MAX_ITER: usize = 100;

/// Fractal visualization modes.
#[derive(Clone, Copy, PartialEq, Eq, Default)]
enum FractalMode {
    /// Mandelbrot set view.
    #[default]
    Mandelbrot,
    /// Julia set view.
    Julia,
    /// Side-by-side comparison.
    SideBySide,
}

impl FractalMode {
    /// Returns the display name of the mode.
    fn name(&self) -> &'static str {
        match self {
            Self::Mandelbrot => "Mandelbrot",
            Self::Julia => "Julia",
            Self::SideBySide => "Side by Side",
        }
    }
}

/// Mandelbrot and Julia set visualization.
pub struct ComplexFractalDemo {
    /// Current visualization mode.
    mode: FractalMode,
    /// Julia set parameter c.
    julia_c: (f64, f64),
    /// Maximum iteration count.
    max_iterations: usize,
    /// Mandelbrot view center.
    mandelbrot_center: (f64, f64),
    /// Mandelbrot view scale.
    mandelbrot_scale: f64,
    /// Julia view center.
    julia_center: (f64, f64),
    /// Julia view scale.
    julia_scale: f64,
    /// Cached Mandelbrot texture.
    mandelbrot_texture: Option<egui::TextureHandle>,
    /// Cached Julia texture.
    julia_texture: Option<egui::TextureHandle>,
    /// Whether Mandelbrot needs refresh.
    mandelbrot_dirty: bool,
    /// Whether Julia needs refresh.
    julia_dirty: bool,
    /// Animation state for c parameter.
    animation: Animation,
    /// Animate c along a path.
    animate_c: bool,
    /// Hover position in Mandelbrot set.
    hover_c: Option<(f64, f64)>,
}

impl Default for ComplexFractalDemo {
    fn default() -> Self {
        Self {
            mode: FractalMode::default(),
            julia_c: (-0.7, 0.27),
            max_iterations: DEFAULT_MAX_ITER,
            mandelbrot_center: (-0.5, 0.0),
            mandelbrot_scale: 1.5,
            julia_center: (0.0, 0.0),
            julia_scale: 1.5,
            mandelbrot_texture: None,
            julia_texture: None,
            mandelbrot_dirty: true,
            julia_dirty: true,
            animation: Animation::with_duration(10.0),
            animate_c: false,
            hover_c: None,
        }
    }
}

impl ComplexFractalDemo {
    /// Computes escape time for Mandelbrot iteration.
    fn mandelbrot_escape(&self, c: Complex<f64>) -> Option<usize> {
        let mut z = Complex::zero();
        for i in 0..self.max_iterations {
            z = z * z + c;
            if z.norm_squared() > 4.0 {
                return Some(i);
            }
        }
        None
    }

    /// Computes escape time for Julia iteration.
    fn julia_escape(&self, z0: Complex<f64>, c: Complex<f64>) -> Option<usize> {
        let mut z = z0;
        for i in 0..self.max_iterations {
            z = z * z + c;
            if z.norm_squared() > 4.0 {
                return Some(i);
            }
        }
        None
    }

    /// Maps iteration count to color.
    fn escape_color(&self, escape: Option<usize>) -> Color32 {
        match escape {
            None => Color32::BLACK,
            Some(i) => {
                let t = (i as f32) / (self.max_iterations as f32);
                // Smooth coloring using HSL
                let hue = 0.6 + 0.4 * t;
                let saturation = 0.8;
                let lightness = 0.1 + 0.7 * (1.0 - t);
                hsl_to_color(hue, saturation, lightness)
            }
        }
    }

    /// Renders Mandelbrot set to texture.
    fn render_mandelbrot(&self, ctx: &egui::Context) -> egui::TextureHandle {
        let mut pixels = vec![Color32::BLACK; IMAGE_RESOLUTION * IMAGE_RESOLUTION];

        let cx = self.mandelbrot_center.0;
        let cy = self.mandelbrot_center.1;
        let scale = self.mandelbrot_scale;

        for py in 0..IMAGE_RESOLUTION {
            for px in 0..IMAGE_RESOLUTION {
                let x = cx + scale * (2.0 * (px as f64) / (IMAGE_RESOLUTION as f64) - 1.0);
                let y = cy + scale * (1.0 - 2.0 * (py as f64) / (IMAGE_RESOLUTION as f64));

                let c = Complex::new(x, y);
                let escape = self.mandelbrot_escape(c);
                pixels[py * IMAGE_RESOLUTION + px] = self.escape_color(escape);
            }
        }

        let image = egui::ColorImage {
            size: [IMAGE_RESOLUTION, IMAGE_RESOLUTION],
            pixels,
        };

        ctx.load_texture("mandelbrot", image, egui::TextureOptions::LINEAR)
    }

    /// Renders Julia set to texture.
    fn render_julia(&self, ctx: &egui::Context) -> egui::TextureHandle {
        let mut pixels = vec![Color32::BLACK; IMAGE_RESOLUTION * IMAGE_RESOLUTION];

        let cx = self.julia_center.0;
        let cy = self.julia_center.1;
        let scale = self.julia_scale;
        let c = Complex::new(self.julia_c.0, self.julia_c.1);

        for py in 0..IMAGE_RESOLUTION {
            for px in 0..IMAGE_RESOLUTION {
                let x = cx + scale * (2.0 * (px as f64) / (IMAGE_RESOLUTION as f64) - 1.0);
                let y = cy + scale * (1.0 - 2.0 * (py as f64) / (IMAGE_RESOLUTION as f64));

                let z0 = Complex::new(x, y);
                let escape = self.julia_escape(z0, c);
                pixels[py * IMAGE_RESOLUTION + px] = self.escape_color(escape);
            }
        }

        let image = egui::ColorImage {
            size: [IMAGE_RESOLUTION, IMAGE_RESOLUTION],
            pixels,
        };

        ctx.load_texture("julia", image, egui::TextureOptions::LINEAR)
    }

    /// Renders a single fractal plot.
    fn render_fractal_plot(
        &mut self,
        ui: &mut egui::Ui,
        is_mandelbrot: bool,
        _ctx: &egui::Context,
    ) {
        let (center, scale, texture) = if is_mandelbrot {
            (
                self.mandelbrot_center,
                self.mandelbrot_scale,
                &self.mandelbrot_texture,
            )
        } else {
            (self.julia_center, self.julia_scale, &self.julia_texture)
        };

        let plot_id = if is_mandelbrot {
            "mandelbrot_plot"
        } else {
            "julia_plot"
        };

        Plot::new(plot_id)
            .data_aspect(1.0)
            .show_axes(false)
            .show_grid(false)
            .allow_zoom(true)
            .allow_drag(true)
            .auto_bounds(egui::Vec2b::new(false, false))
            .include_x(center.0 - scale)
            .include_x(center.0 + scale)
            .include_y(center.1 - scale)
            .include_y(center.1 + scale)
            .show(ui, |plot_ui| {
                if let Some(texture) = texture {
                    let image = PlotImage::new(
                        texture,
                        PlotPoint::new(center.0, center.1),
                        [(2.0 * scale) as f32, (2.0 * scale) as f32],
                    );
                    plot_ui.image(image);
                }

                // Show c parameter marker on Mandelbrot
                if is_mandelbrot {
                    plot_ui.points(
                        Points::new(vec![[self.julia_c.0, self.julia_c.1]])
                            .color(Color32::WHITE)
                            .radius(5.0)
                            .filled(true)
                            .name("c parameter"),
                    );

                    // Track hover for c selection
                    if let Some(pos) = plot_ui.pointer_coordinate() {
                        self.hover_c = Some((pos.x, pos.y));
                    } else {
                        self.hover_c = None;
                    }
                }
            });
    }
}

impl VisualizationApp for ComplexFractalDemo {
    fn name(&self) -> &'static str {
        "Mandelbrot & Julia Sets"
    }

    fn update(&mut self, dt: f32) {
        self.animation.update(dt);

        if self.animate_c && self.animation.playing {
            // Animate c along a cardioid path near the Mandelbrot boundary
            let t = self.animation.angle();
            let r = 0.26;
            self.julia_c.0 = -0.75 + r * (2.0 * t as f64).cos();
            self.julia_c.1 = r * (2.0 * t as f64).sin();
            self.julia_dirty = true;
        }
    }

    fn render(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();

        // Refresh textures if needed
        if self.mandelbrot_dirty || self.mandelbrot_texture.is_none() {
            self.mandelbrot_texture = Some(self.render_mandelbrot(&ctx));
            self.mandelbrot_dirty = false;
        }

        if self.julia_dirty || self.julia_texture.is_none() {
            self.julia_texture = Some(self.render_julia(&ctx));
            self.julia_dirty = false;
        }

        match self.mode {
            FractalMode::Mandelbrot => {
                self.render_fractal_plot(ui, true, &ctx);
            }
            FractalMode::Julia => {
                self.render_fractal_plot(ui, false, &ctx);
            }
            FractalMode::SideBySide => {
                let available = ui.available_size();
                let is_wide = available.x > available.y;

                if is_wide {
                    // Side by side (columns) for wide layouts
                    ui.columns(2, |columns| {
                        columns[0].vertical_centered(|ui| {
                            ui.label("Mandelbrot Set");
                        });
                        self.render_fractal_plot(&mut columns[0], true, &ctx);

                        columns[1].vertical_centered(|ui| {
                            ui.label("Julia Set");
                        });
                        self.render_fractal_plot(&mut columns[1], false, &ctx);
                    });
                } else {
                    // Stacked (rows) for tall layouts
                    ui.vertical_centered(|ui| {
                        ui.label("Mandelbrot Set");
                    });
                    self.render_fractal_plot(ui, true, &ctx);

                    ui.add_space(8.0);

                    ui.vertical_centered(|ui| {
                        ui.label("Julia Set");
                    });
                    self.render_fractal_plot(ui, false, &ctx);
                }
            }
        }
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        group_header(ui, "Mode");

        ui.horizontal(|ui| {
            for mode in [
                FractalMode::Mandelbrot,
                FractalMode::Julia,
                FractalMode::SideBySide,
            ] {
                if ui
                    .selectable_label(self.mode == mode, mode.name())
                    .clicked()
                {
                    self.mode = mode;
                }
            }
        });

        section_separator(ui, Some("Julia Parameter c"));

        // c parameter sliders
        ui.horizontal(|ui| {
            ui.label("Re(c):");
            if ui
                .add(egui::Slider::new(&mut self.julia_c.0, -2.0..=2.0).step_by(0.01))
                .changed()
            {
                self.julia_dirty = true;
            }
        });

        ui.horizontal(|ui| {
            ui.label("Im(c):");
            if ui
                .add(egui::Slider::new(&mut self.julia_c.1, -2.0..=2.0).step_by(0.01))
                .changed()
            {
                self.julia_dirty = true;
            }
        });

        // Click to set c
        if self.mode != FractalMode::Julia {
            if let Some((cx, cy)) = self.hover_c {
                if ui.input(|i| i.pointer.button_clicked(egui::PointerButton::Primary)) {
                    self.julia_c = (cx, cy);
                    self.julia_dirty = true;
                }
            }
            ui.label("Click on Mandelbrot to set c");
        }

        // Preset c values
        ui.horizontal(|ui| {
            if ui.button("Dendrite").clicked() {
                self.julia_c = (0.0, 1.0);
                self.julia_dirty = true;
            }
            if ui.button("Rabbit").clicked() {
                self.julia_c = (-0.123, 0.745);
                self.julia_dirty = true;
            }
            if ui.button("San Marco").clicked() {
                self.julia_c = (-0.75, 0.0);
                self.julia_dirty = true;
            }
        });

        section_separator(ui, Some("Animation"));

        ui.checkbox(&mut self.animate_c, "Animate c");
        if self.animate_c {
            animation_controls(ui, &mut self.animation);
        }

        section_separator(ui, Some("Rendering"));

        ui.horizontal(|ui| {
            ui.label("Max iterations:");
            if ui
                .add(egui::Slider::new(&mut self.max_iterations, 20..=500).logarithmic(true))
                .changed()
            {
                self.mandelbrot_dirty = true;
                self.julia_dirty = true;
            }
        });

        // Reset buttons
        ui.horizontal(|ui| {
            if ui.button("Reset Mandelbrot").clicked() {
                self.mandelbrot_center = (-0.5, 0.0);
                self.mandelbrot_scale = 1.5;
                self.mandelbrot_dirty = true;
            }
            if ui.button("Reset Julia").clicked() {
                self.julia_center = (0.0, 0.0);
                self.julia_scale = 1.5;
                self.julia_dirty = true;
            }
        });

        // Current c display
        section_separator(ui, Some("Current c"));
        info_box(
            ui,
            &format!(
                "c = {:.4} + {:.4}i\n\n\
                 The Julia set J_c is connected\n\
                 if and only if c is IN the\n\
                 Mandelbrot set.",
                self.julia_c.0, self.julia_c.1
            ),
        );
    }

    fn info(&self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.label(format!(
                "c = {:.3} + {:.3}i",
                self.julia_c.0, self.julia_c.1
            ));
            ui.separator();
            ui.label(format!("Max iter: {}", self.max_iterations));
            ui.separator();
            ui.label(self.mode.name());
        });
    }

    fn educational_content(&self) -> Option<EducationalContent> {
        Some(COMPLEX_FRACTAL_EDUCATION)
    }
}

/// Educational content for the Mandelbrot/Julia fractal demo.
const COMPLEX_FRACTAL_EDUCATION: EducationalContent = EducationalContent {
    title: "Mandelbrot and Julia Sets",

    overview: "\
The Mandelbrot set and Julia sets are iconic fractals defined by the simple
iteration z \u{2192} z\u{00b2} + c. Despite this simplicity, they exhibit
infinite complexity and self-similarity at all scales.

The Mandelbrot set is the 'master catalog' of Julia sets: each point c in the
complex plane corresponds to a unique Julia set J_c. Points INSIDE the
Mandelbrot set have CONNECTED Julia sets; points OUTSIDE have 'dust' Julia sets.",

    math_background: "\
Both fractals use the iteration:
    z_{n+1} = z_n\u{00b2} + c

MANDELBROT SET: Fix z\u{2080} = 0, vary c.
    M = {c \u{2208} C : |z_n| stays bounded}

JULIA SET: Fix c, vary z\u{2080}.
    J_c = {z\u{2080} \u{2208} C : |z_n| stays bounded}

The ESCAPE TIME (how many iterations before |z| > 2) creates
the colorful images. Points that never escape are in the set.

KEY THEOREM: J_c is connected \u{21d4} c \u{2208} M

This beautiful connection means you can 'see' the Julia set's
structure just by looking at where c lies in the Mandelbrot set.",

    how_to_use: "\
\u{2022} Click on Mandelbrot to select c for Julia set
\u{2022} Use presets to see famous Julia sets
\u{2022} Animate c to watch Julia sets transform
\u{2022} Increase iterations for more detail
\u{2022} Zoom in to see infinite detail",

    key_concepts: "\
\u{2022} z \u{2192} z\u{00b2} + c: simple rule, infinite complexity
\u{2022} Mandelbrot: c varies, z\u{2080} = 0
\u{2022} Julia: z\u{2080} varies, c fixed
\u{2022} Connected Julia \u{21d4} c in Mandelbrot
\u{2022} Self-similar at all zoom levels
\u{2022} Boundary has Hausdorff dimension 2",

    resources: &[
        (
            "Mandelbrot Set - Wikipedia",
            "https://en.wikipedia.org/wiki/Mandelbrot_set",
        ),
        (
            "Julia Set - Wikipedia",
            "https://en.wikipedia.org/wiki/Julia_set",
        ),
    ],
};
