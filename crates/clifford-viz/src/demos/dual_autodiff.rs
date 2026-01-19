//! Dual Number Automatic Differentiation Visualization
//!
//! This demo visualizes how dual numbers compute exact derivatives through
//! simple arithmetic. The key insight is that for any smooth function f:
//!
//!     f(x + ε) = f(x) + f'(x)·ε
//!
//! where ε² = 0. This gives us the exact derivative with no truncation error!

use crate::common::prelude::*;
use clifford::specialized::dual::Dual;
use egui_plot::{Line, Plot, PlotPoints, Points};

/// Viewport X bounds for the function plot.
const VIEWPORT_X: f64 = 5.0;
/// Viewport Y bounds for the function plot.
const VIEWPORT_Y: f64 = 3.0;

/// Available functions for differentiation.
#[derive(Clone, Copy, PartialEq, Eq, Default)]
enum DiffFunction {
    /// f(x) = x²
    #[default]
    Square,
    /// f(x) = x³
    Cube,
    /// f(x) = sin(x)
    Sin,
    /// f(x) = cos(x)
    Cos,
    /// f(x) = eˣ
    Exp,
    /// f(x) = ln(x)
    Log,
    /// f(x) = √x
    Sqrt,
    /// f(x) = x³ - 2x + 1
    Custom,
}

impl DiffFunction {
    /// Returns the display name of the function.
    fn name(&self) -> &'static str {
        match self {
            Self::Square => "x\u{00b2}",
            Self::Cube => "x\u{00b3}",
            Self::Sin => "sin(x)",
            Self::Cos => "cos(x)",
            Self::Exp => "e\u{02e3}",
            Self::Log => "ln(x)",
            Self::Sqrt => "\u{221a}x",
            Self::Custom => "x\u{00b3} - 2x + 1",
        }
    }

    /// Returns the display name of the derivative formula.
    fn derivative_name(&self) -> &'static str {
        match self {
            Self::Square => "2x",
            Self::Cube => "3x\u{00b2}",
            Self::Sin => "cos(x)",
            Self::Cos => "-sin(x)",
            Self::Exp => "e\u{02e3}",
            Self::Log => "1/x",
            Self::Sqrt => "1/(2\u{221a}x)",
            Self::Custom => "3x\u{00b2} - 2",
        }
    }

    /// Evaluates the function at x using standard arithmetic.
    fn evaluate(&self, x: f64) -> f64 {
        match self {
            Self::Square => x * x,
            Self::Cube => x * x * x,
            Self::Sin => x.sin(),
            Self::Cos => x.cos(),
            Self::Exp => x.exp(),
            Self::Log => x.ln(),
            Self::Sqrt => x.sqrt(),
            Self::Custom => x * x * x - 2.0 * x + 1.0,
        }
    }

    /// Evaluates the function using dual number arithmetic.
    fn evaluate_dual(&self, d: Dual<f64>) -> Dual<f64> {
        match self {
            Self::Square => d * d,
            Self::Cube => d * d * d,
            Self::Sin => d.sin(),
            Self::Cos => d.cos(),
            Self::Exp => d.exp(),
            Self::Log => d.ln(),
            Self::Sqrt => d.sqrt(),
            Self::Custom => d * d * d - Dual::from_real(2.0) * d + Dual::from_real(1.0),
        }
    }

    /// Returns the minimum x value for this function's domain.
    fn min_x(&self) -> f64 {
        match self {
            Self::Log | Self::Sqrt => 0.1,
            _ => -VIEWPORT_X,
        }
    }

    /// Returns all available functions.
    fn all() -> &'static [DiffFunction] {
        &[
            Self::Square,
            Self::Cube,
            Self::Sin,
            Self::Cos,
            Self::Exp,
            Self::Log,
            Self::Sqrt,
            Self::Custom,
        ]
    }
}

/// Dual number automatic differentiation visualization.
pub struct DualAutodiffDemo {
    /// Selected function.
    function: DiffFunction,
    /// Evaluation point x.
    x: f32,
    /// Whether to show the tangent line.
    show_tangent: bool,
    /// Whether to show the secant line (numerical approximation).
    show_secant: bool,
    /// Step size for numerical differentiation.
    secant_h: f32,
    /// Animation state.
    animation: Animation,
    /// Whether x is being animated.
    animate_x: bool,
}

impl Default for DualAutodiffDemo {
    fn default() -> Self {
        Self {
            function: DiffFunction::default(),
            x: 1.0,
            show_tangent: true,
            show_secant: true,
            secant_h: 0.5,
            animation: Animation::with_duration(6.0),
            animate_x: false,
        }
    }
}

impl DualAutodiffDemo {
    /// Computes f(x) and f'(x) using dual numbers.
    fn compute_derivative(&self, x: f64) -> (f64, f64) {
        let dual_x = Dual::variable(x);
        let result = self.function.evaluate_dual(dual_x);
        (result.real(), result.dual())
    }

    /// Computes numerical derivative using forward difference.
    fn numerical_derivative(&self, x: f64, h: f64) -> f64 {
        let f_x = self.function.evaluate(x);
        let f_xh = self.function.evaluate(x + h);
        (f_xh - f_x) / h
    }

    /// Generates points for the function curve.
    fn function_curve(&self) -> Vec<[f64; 2]> {
        let min_x = self.function.min_x();
        let n = 200;
        (0..n)
            .filter_map(|i| {
                let x = min_x + (VIEWPORT_X - min_x) * 2.0 * (i as f64) / (n as f64);
                let y = self.function.evaluate(x);
                if y.is_finite() && y.abs() < VIEWPORT_Y * 2.0 {
                    Some([x, y])
                } else {
                    None
                }
            })
            .collect()
    }

    /// Generates points for the tangent line at x.
    fn tangent_line(&self, x: f64, f_x: f64, f_prime_x: f64) -> Vec<[f64; 2]> {
        let x_range = VIEWPORT_X;
        vec![
            [x - x_range, f_x - f_prime_x * x_range],
            [x + x_range, f_x + f_prime_x * x_range],
        ]
    }

    /// Generates points for the secant line.
    fn secant_line(&self, x: f64, h: f64) -> Vec<[f64; 2]> {
        let f_x = self.function.evaluate(x);
        let f_xh = self.function.evaluate(x + h);
        let slope = (f_xh - f_x) / h;
        let x_range = VIEWPORT_X;
        vec![
            [x - x_range, f_x - slope * x_range],
            [x + x_range, f_x + slope * x_range],
        ]
    }
}

impl VisualizationApp for DualAutodiffDemo {
    fn name(&self) -> &'static str {
        "Dual Numbers - Automatic Differentiation"
    }

    fn update(&mut self, dt: f32) {
        self.animation.update(dt);

        if self.animate_x && self.animation.playing {
            let min_x = self.function.min_x() as f32 + 0.1;
            let max_x = VIEWPORT_X as f32 - 0.1;
            let t = self.animation.progress();
            // Oscillate back and forth
            let t_osc = (t * 2.0 * std::f32::consts::PI).sin() * 0.5 + 0.5;
            self.x = min_x + (max_x - min_x) * t_osc;
        }
    }

    fn render(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        let x = f64::from(self.x);
        let (f_x, f_prime_x) = self.compute_derivative(x);

        Plot::new("dual_autodiff_plot")
            .data_aspect(1.0)
            .show_axes(true)
            .show_grid(false)
            .allow_zoom(true)
            .allow_drag(true)
            .auto_bounds(egui::Vec2b::new(false, false))
            .include_x(-VIEWPORT_X)
            .include_x(VIEWPORT_X)
            .include_y(-VIEWPORT_Y)
            .include_y(VIEWPORT_Y)
            .show(ui, |plot_ui| {
                // Draw coordinate grid
                for l in grid_2d(&ctx, VIEWPORT_X, 1.0) {
                    plot_ui.line(l);
                }

                // Draw function curve
                let curve_points = self.function_curve();
                plot_ui.line(
                    Line::new(PlotPoints::new(curve_points))
                        .color(point(&ctx))
                        .width(2.0)
                        .name(format!("f(x) = {}", self.function.name())),
                );

                // Draw tangent line (from dual numbers)
                if self.show_tangent && f_x.is_finite() && f_prime_x.is_finite() {
                    let tangent_points = self.tangent_line(x, f_x, f_prime_x);
                    plot_ui.line(
                        Line::new(PlotPoints::new(tangent_points))
                            .color(plane(&ctx))
                            .width(1.5)
                            .name("Tangent (exact)"),
                    );
                }

                // Draw secant line (numerical approximation)
                if self.show_secant && f_x.is_finite() {
                    let h = f64::from(self.secant_h);
                    let secant_points = self.secant_line(x, h);
                    plot_ui.line(
                        Line::new(PlotPoints::new(secant_points))
                            .color(with_alpha(line(&ctx), 150))
                            .width(1.5)
                            .name(format!("Secant (h = {:.2})", h)),
                    );

                    // Mark the two points used for secant
                    let f_xh = self.function.evaluate(x + h);
                    if f_xh.is_finite() {
                        plot_ui.points(
                            Points::new(vec![[x + h, f_xh]])
                                .color(line(&ctx))
                                .radius(4.0)
                                .filled(true)
                                .name("(x+h, f(x+h))"),
                        );
                    }
                }

                // Mark the evaluation point
                if f_x.is_finite() {
                    plot_ui.points(
                        Points::new(vec![[x, f_x]])
                            .color(motor(&ctx))
                            .radius(6.0)
                            .filled(true)
                            .name("(x, f(x))"),
                    );
                }
            });
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();

        group_header(ui, "Function");

        // Function selection grid
        egui::Grid::new("function_grid")
            .num_columns(4)
            .spacing([8.0, 4.0])
            .show(ui, |ui| {
                for (i, func) in DiffFunction::all().iter().enumerate() {
                    if ui
                        .selectable_label(self.function == *func, func.name())
                        .clicked()
                    {
                        self.function = *func;
                        // Adjust x if it's out of range for new function
                        let min_x = func.min_x() as f32;
                        if self.x < min_x {
                            self.x = min_x + 0.5;
                        }
                    }
                    if (i + 1) % 4 == 0 {
                        ui.end_row();
                    }
                }
            });

        section_separator(ui, Some("Evaluation Point"));

        // x slider
        let min_x = self.function.min_x() as f32;
        ui.horizontal(|ui| {
            ui.label("x =");
            ui.add(egui::Slider::new(&mut self.x, min_x..=4.5).step_by(0.01));
        });

        // Animation controls
        ui.checkbox(&mut self.animate_x, "Animate x");
        if self.animate_x {
            animation_controls(ui, &mut self.animation);
        }

        section_separator(ui, Some("Dual Number Computation"));

        // Show the dual computation
        let x = f64::from(self.x);
        let (f_x, f_prime_x) = self.compute_derivative(x);

        info_box(
            ui,
            &format!(
                "Input: x + \u{03b5} = {:.4} + \u{03b5}\n\n\
                 f(x + \u{03b5}) = {}\n\n\
                 Result: {:.4} + {:.4}\u{03b5}\n\n\
                 Therefore:\n\
                   f(x)  = {:.4}\n\
                   f'(x) = {:.4}",
                x,
                self.function.name(),
                f_x,
                f_prime_x,
                f_x,
                f_prime_x
            ),
        );

        section_separator(ui, Some("Display Options"));

        ui.horizontal(|ui| {
            ui.checkbox(&mut self.show_tangent, "Tangent");
            ui.colored_label(plane(&ctx), "(exact)");
        });

        ui.horizontal(|ui| {
            ui.checkbox(&mut self.show_secant, "Secant");
            ui.colored_label(line(&ctx), "(numerical)");
        });

        if self.show_secant {
            ui.horizontal(|ui| {
                ui.label("h =");
                ui.add(egui::Slider::new(&mut self.secant_h, 0.01..=1.0).logarithmic(true));
            });

            // Compare derivatives
            let numerical = self.numerical_derivative(x, f64::from(self.secant_h));
            let error = (numerical - f_prime_x).abs();

            section_separator(ui, Some("Comparison"));
            info_box(
                ui,
                &format!(
                    "Dual (exact):     {:.6}\n\
                     Numerical (h={:.2}): {:.6}\n\
                     Error:            {:.2e}\n\n\
                     Dual numbers give EXACT\n\
                     derivatives with NO\n\
                     truncation error!",
                    f_prime_x, self.secant_h, numerical, error
                ),
            );
        }

        // Derivative formula
        section_separator(ui, Some("Derivative"));
        ui.horizontal(|ui| {
            ui.label("f(x) =");
            ui.label(self.function.name());
        });
        ui.horizontal(|ui| {
            ui.label("f'(x) =");
            ui.label(self.function.derivative_name());
        });
    }

    fn info(&self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        let x = f64::from(self.x);
        let (f_x, f_prime_x) = self.compute_derivative(x);

        ui.horizontal(|ui| {
            ui.label(format!("f(x) = {}", self.function.name()));
            ui.separator();
            ui.label(format!("x = {:.2}", x));
            ui.separator();
            ui.label(format!("f(x) = {:.3}", f_x));
            ui.separator();
            ui.colored_label(plane(&ctx), format!("f'(x) = {:.3}", f_prime_x));
        });
    }

    fn educational_content(&self) -> Option<EducationalContent> {
        Some(DUAL_AUTODIFF_EDUCATION)
    }
}

/// Educational content for the dual autodiff demo.
const DUAL_AUTODIFF_EDUCATION: EducationalContent = EducationalContent {
    title: "Automatic Differentiation with Dual Numbers",

    overview: "\
Dual numbers provide a beautiful way to compute EXACT derivatives without
limits, finite differences, or symbolic manipulation. By extending the
real numbers with an infinitesimal \u{03b5} where \u{03b5}\u{00b2} = 0, function
evaluation automatically tracks the derivative!

This is 'forward-mode automatic differentiation' - the foundation for
computing gradients in machine learning, physics simulations, and more.",

    math_background: "\
A DUAL NUMBER has the form a + b\u{03b5} where \u{03b5}\u{00b2} = 0.

Multiplication follows naturally:
    (a + b\u{03b5})(c + d\u{03b5}) = ac + (ad + bc)\u{03b5}

The magic: for ANY smooth function f, Taylor expansion gives:
    f(x + \u{03b5}) = f(x) + f'(x)\u{03b5} + f''(x)\u{03b5}\u{00b2}/2 + ...
                = f(x) + f'(x)\u{03b5}  (since \u{03b5}\u{00b2} = 0!)

So just evaluate f(x + \u{03b5}) and extract:
  \u{2022} Real part = f(x)
  \u{2022} Dual part = f'(x)

Example: f(x) = x\u{00b2}
    f(x + \u{03b5}) = (x + \u{03b5})\u{00b2}
                = x\u{00b2} + 2x\u{03b5} + \u{03b5}\u{00b2}
                = x\u{00b2} + 2x\u{03b5}  (\u{03b5}\u{00b2} = 0)
    \u{2192} f'(x) = 2x  \u{2714}",

    how_to_use: "\
\u{2022} Select a function to differentiate
\u{2022} Adjust x to see the tangent line move
\u{2022} Enable 'Secant' to compare with numerical differentiation
\u{2022} Adjust h to see numerical error increase/decrease
\u{2022} Animate x to watch the tangent track the curve",

    key_concepts: "\
\u{2022} \u{03b5}\u{00b2} = 0: the nilpotent infinitesimal
\u{2022} f(x + \u{03b5}) = f(x) + f'(x)\u{03b5}: automatic derivatives
\u{2022} NO truncation error (unlike finite differences)
\u{2022} NO symbolic manipulation needed
\u{2022} Chain rule works automatically through composition
\u{2022} Foundation of forward-mode autodiff in ML",

    resources: &[
        (
            "Dual Numbers - Wikipedia",
            "https://en.wikipedia.org/wiki/Dual_number",
        ),
        (
            "Automatic Differentiation - Wikipedia",
            "https://en.wikipedia.org/wiki/Automatic_differentiation",
        ),
    ],
};
