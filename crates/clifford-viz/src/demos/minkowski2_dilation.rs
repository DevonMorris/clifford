//! Time Dilation and Twin Paradox Visualization
//!
//! This demo visualizes the famous twin paradox from special relativity:
//!
//! - A "stay-at-home" twin follows a vertical worldline (stationary in the lab frame)
//! - A "traveling" twin takes a round trip at relativistic speed
//! - The traveling twin ages LESS due to time dilation!
//!
//! ## Mathematical Background
//!
//! Proper time along a worldline is the integral:
//!
//! ```text
//! τ = ∫ √(dt² - dx²) = ∫ √(1 - v²) dt
//! ```
//!
//! For the stay-at-home twin (dx = 0): τ = T (coordinate time)
//! For the traveler moving at speed v: τ = T/γ where γ = 1/√(1-v²)
//!
//! The KEY INSIGHT: proper time is a PATH LENGTH in spacetime!
//! The straight path (stationary) has the LONGEST proper time.
//! This is opposite to Euclidean geometry where straight lines are shortest.

use crate::common::prelude::*;
use clifford::specialized::minkowski::dim2::Vector;
use egui::Color32;
use egui_plot::{Line, Plot, PlotPoints, Points, Text};

/// Fixed viewport bounds.
const VIEWPORT_BOUNDS: f64 = 12.0;

/// A point on a worldline with accumulated proper time.
///
/// This represents a point in 1+1D spacetime along with the accumulated
/// proper time (τ) experienced by an observer following this worldline.
/// The spacetime position is stored using the clifford library's [`Vector`] type.
#[derive(Clone, Debug)]
struct WorldlinePoint {
    /// Position in spacetime as a Minkowski vector.
    /// Uses clifford's Vector<f64> with spacelike `x` (e₁² = -1) and timelike `t` (e₂² = +1).
    /// This gives the physics (timelike-positive) convention: v² = t² - x².
    position: Vector<f64>,
    /// Accumulated proper time along worldline.
    tau: f64,
}

impl WorldlinePoint {
    /// Creates a new worldline point at coordinates (x, t) with accumulated proper time τ.
    fn new(x: f64, t: f64, tau: f64) -> Self {
        Self {
            position: Vector::new(x, t),
            tau,
        }
    }

    /// Returns the spatial coordinate.
    fn x(&self) -> f64 {
        self.position.x()
    }

    /// Returns the time coordinate.
    fn t(&self) -> f64 {
        self.position.t()
    }
}

/// Demo state for time dilation visualization.
pub struct Minkowski2DilationDemo {
    /// Travel velocity (as fraction of c).
    travel_velocity: f32,
    /// Total coordinate time for the journey.
    total_time: f32,
    /// Current animation progress (0 to 1).
    journey_progress: f32,
    /// Animation state.
    animation: Animation,
    /// Whether to show light cones.
    show_light_cones: bool,
    /// Whether to show coordinate grid.
    show_grid: bool,
    /// Whether to show simultaneity lines from traveler.
    show_simultaneity: bool,
    /// Whether to show proper time ticks.
    show_proper_time_ticks: bool,
    /// Number of proper time ticks to show.
    tick_count: i32,
}

impl Default for Minkowski2DilationDemo {
    fn default() -> Self {
        Self {
            travel_velocity: 0.8,
            total_time: 10.0,
            journey_progress: 1.0, // Show full journey by default
            animation: Animation::with_duration(6.0),
            show_light_cones: false,
            show_grid: true,
            show_simultaneity: true,
            show_proper_time_ticks: true,
            tick_count: 10,
        }
    }
}

impl Minkowski2DilationDemo {
    /// Compute the stay-at-home twin's worldline (vertical).
    fn stationary_worldline(&self, progress: f32) -> Vec<WorldlinePoint> {
        let t_max = f64::from(self.total_time) * f64::from(progress);
        let n = 100;

        (0..=n)
            .map(|i| {
                let t = t_max * (i as f64) / (n as f64);
                // For stationary observer, proper time = coordinate time
                WorldlinePoint::new(0.0, t, t)
            })
            .collect()
    }

    /// Compute the traveling twin's worldline (outbound + return).
    fn traveler_worldline(&self, progress: f32) -> Vec<WorldlinePoint> {
        let v = f64::from(self.travel_velocity);
        let t_total = f64::from(self.total_time) * f64::from(progress);
        let t_turn = f64::from(self.total_time) / 2.0; // Turnaround at midpoint
        let gamma = 1.0 / (1.0 - v * v).sqrt();
        let n = 100;

        let mut points = Vec::with_capacity(n + 1);

        for i in 0..=n {
            let t = t_total * (i as f64) / (n as f64);

            // Outbound phase
            let (x, tau) = if t <= t_turn.min(t_total) {
                let x = v * t;
                // Proper time: τ = t/γ for constant velocity
                let tau = t / gamma;
                (x, tau)
            } else {
                // Return phase
                let t_out = t_turn;
                let tau_out = t_out / gamma;
                let x_turn = v * t_turn;

                // Return leg (traveling at -v)
                let dt_return = t - t_turn;
                let x = x_turn - v * dt_return;
                let tau = tau_out + dt_return / gamma;
                (x, tau)
            };

            points.push(WorldlinePoint::new(x, t, tau));
        }

        points
    }

    /// Compute proper time for the full journey.
    fn compute_proper_times(&self) -> (f64, f64) {
        let v = f64::from(self.travel_velocity);
        let t_total = f64::from(self.total_time);
        let gamma = 1.0 / (1.0 - v * v).sqrt();

        let tau_stationary = t_total;
        let tau_traveler = t_total / gamma; // Same γ for both legs

        (tau_stationary, tau_traveler)
    }

    /// Draw light cone from origin.
    fn draw_light_cone(&self, plot_ui: &mut egui_plot::PlotUi, ctx: &egui::Context) {
        let light_color = with_alpha(circle(ctx), 80);
        let extent = VIEWPORT_BOUNDS;

        // Future light cone
        let future_right: PlotPoints = (0..=50)
            .map(|i| {
                let t = extent * (i as f64) / 50.0;
                [t, t]
            })
            .collect();
        let future_left: PlotPoints = (0..=50)
            .map(|i| {
                let t = extent * (i as f64) / 50.0;
                [-t, t]
            })
            .collect();

        plot_ui.line(Line::new(future_right).color(light_color).width(1.5));
        plot_ui.line(Line::new(future_left).color(light_color).width(1.5));
    }

    /// Draw lines of simultaneity from the traveler's perspective.
    fn draw_simultaneity_lines(
        &self,
        plot_ui: &mut egui_plot::PlotUi,
        ctx: &egui::Context,
        traveler_points: &[WorldlinePoint],
    ) {
        let v = f64::from(self.travel_velocity);
        let t_turn = f64::from(self.total_time) / 2.0;
        let sim_color = with_alpha(motor(ctx), 60);

        // Draw a few simultaneity lines along the journey
        let num_lines = 5;
        let t_max = traveler_points.last().map_or(0.0, |p| p.t());

        for i in 1..num_lines {
            let fraction = (i as f64) / (num_lines as f64);
            let t = t_max * fraction;

            // Find the traveler's position at this time
            let x_traveler = if t <= t_turn {
                v * t
            } else {
                let x_turn = v * t_turn;
                x_turn - v * (t - t_turn)
            };

            // Simultaneity line slope depends on velocity direction
            let current_v = if t <= t_turn { v } else { -v };

            // Line of simultaneity has slope v in the traveler's frame
            // t_sim = t + v*(x - x_traveler) in lab frame
            let extent = 8.0;
            let points: PlotPoints = vec![
                [x_traveler - extent, t - current_v * extent],
                [x_traveler + extent, t + current_v * extent],
            ]
            .into();

            plot_ui.line(
                Line::new(points)
                    .color(sim_color)
                    .width(1.0)
                    .style(egui_plot::LineStyle::Dashed { length: 4.0 }),
            );
        }
    }

    /// Draw proper time tick marks along a worldline.
    fn draw_proper_time_ticks(
        &self,
        plot_ui: &mut egui_plot::PlotUi,
        points: &[WorldlinePoint],
        color: Color32,
        offset_x: f64,
    ) {
        if points.len() < 2 {
            return;
        }

        let tau_max = points.last().unwrap().tau;
        if tau_max < 0.01 {
            return;
        }

        let tick_interval = tau_max / (self.tick_count as f64);
        let mut current_tick = tick_interval;

        // Interpolate tick positions along the worldline
        for window in points.windows(2) {
            let p1 = &window[0];
            let p2 = &window[1];

            while current_tick <= p2.tau && current_tick > p1.tau {
                // Interpolate position
                let frac = (current_tick - p1.tau) / (p2.tau - p1.tau);
                let x = p1.x() + frac * (p2.x() - p1.x());
                let t = p1.t() + frac * (p2.t() - p1.t());

                // Draw tick mark perpendicular to worldline
                let dx = p2.x() - p1.x();
                let dt = p2.t() - p1.t();
                let len = (dx * dx + dt * dt).sqrt();
                let tick_size = 0.15;
                // Perpendicular direction (rotate 90°)
                let nx = -dt / len * tick_size;
                let nt = dx / len * tick_size;

                let tick_points: PlotPoints =
                    vec![[x + nx + offset_x, t + nt], [x - nx + offset_x, t - nt]].into();

                plot_ui.line(Line::new(tick_points).color(color).width(2.0));

                current_tick += tick_interval;
            }
        }
    }
}

impl VisualizationApp for Minkowski2DilationDemo {
    fn name(&self) -> &'static str {
        "Time Dilation (Twin Paradox)"
    }

    fn update(&mut self, dt: f32) {
        self.animation.update(dt);
        if self.animation.playing {
            self.journey_progress = self.animation.progress();
        }
    }

    fn render(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();

        let stationary_points = self.stationary_worldline(self.journey_progress);
        let traveler_points = self.traveler_worldline(self.journey_progress);

        Plot::new("minkowski2_dilation")
            .data_aspect(1.0)
            .show_axes(false)
            .show_grid(false)
            .allow_zoom(true)
            .allow_drag(true)
            .x_axis_label("x (space)")
            .y_axis_label("t (time)")
            .show(ui, |plot_ui| {
                // Coordinate grid
                if self.show_grid {
                    for l in grid_2d(&ctx, VIEWPORT_BOUNDS, 1.0) {
                        plot_ui.line(l);
                    }
                    // Axes (renamed to avoid shadowing color functions)
                    let x_axis_line: PlotPoints =
                        vec![[-VIEWPORT_BOUNDS, 0.0], [VIEWPORT_BOUNDS, 0.0]].into();
                    plot_ui.line(
                        Line::new(x_axis_line)
                            .color(x_axis(&ctx))
                            .width(2.0)
                            .name("x (space)"),
                    );
                    let t_axis_line: PlotPoints = vec![[0.0, -1.0], [0.0, VIEWPORT_BOUNDS]].into();
                    plot_ui.line(
                        Line::new(t_axis_line)
                            .color(t_axis(&ctx))
                            .width(2.0)
                            .name("t (time)"),
                    );
                }

                // Light cones
                if self.show_light_cones {
                    self.draw_light_cone(plot_ui, &ctx);
                }

                // Simultaneity lines
                if self.show_simultaneity && !traveler_points.is_empty() {
                    self.draw_simultaneity_lines(plot_ui, &ctx, &traveler_points);
                }

                // Stay-at-home twin's worldline (vertical, blue)
                let stationary_plot: PlotPoints =
                    stationary_points.iter().map(|p| [p.x(), p.t()]).collect();
                plot_ui.line(
                    Line::new(stationary_plot)
                        .color(point(&ctx))
                        .width(3.0)
                        .name("Stay-at-home"),
                );

                // Traveler's worldline (V-shape, green)
                let traveler_plot: PlotPoints =
                    traveler_points.iter().map(|p| [p.x(), p.t()]).collect();
                plot_ui.line(
                    Line::new(traveler_plot)
                        .color(plane(&ctx))
                        .width(3.0)
                        .name("Traveler"),
                );

                // Proper time ticks
                if self.show_proper_time_ticks {
                    self.draw_proper_time_ticks(plot_ui, &stationary_points, point(&ctx), -0.2);
                    self.draw_proper_time_ticks(plot_ui, &traveler_points, plane(&ctx), 0.0);
                }

                // Event markers
                // Start event
                plot_ui.points(
                    Points::new(vec![[0.0, 0.0]])
                        .color(active(&ctx))
                        .radius(8.0)
                        .filled(true)
                        .name("Departure"),
                );
                plot_ui.text(
                    Text::new(egui_plot::PlotPoint::new(-0.5, -0.5), "Start")
                        .color(text_primary(&ctx)),
                );

                // Turnaround event (if journey has progressed past midpoint)
                let t_turn = f64::from(self.total_time) / 2.0;
                let t_current = f64::from(self.total_time) * f64::from(self.journey_progress);
                if t_current > t_turn {
                    let x_turn = f64::from(self.travel_velocity) * t_turn;
                    plot_ui.points(
                        Points::new(vec![[x_turn, t_turn]])
                            .color(circle(&ctx))
                            .radius(8.0)
                            .filled(true)
                            .name("Turnaround"),
                    );
                    plot_ui.text(
                        Text::new(
                            egui_plot::PlotPoint::new(x_turn + 0.3, t_turn + 0.3),
                            "Turn",
                        )
                        .color(text_primary(&ctx)),
                    );
                }

                // Reunion event (if journey complete)
                if self.journey_progress > 0.99 {
                    let t_total = f64::from(self.total_time);
                    plot_ui.points(
                        Points::new(vec![[0.0, t_total]])
                            .color(active(&ctx))
                            .radius(8.0)
                            .filled(true)
                            .name("Reunion"),
                    );
                    plot_ui.text(
                        Text::new(egui_plot::PlotPoint::new(0.5, t_total + 0.3), "Reunion")
                            .color(text_primary(&ctx)),
                    );

                    // Show ages at reunion
                    let (tau_home, tau_travel) = self.compute_proper_times();
                    plot_ui.text(
                        Text::new(
                            egui_plot::PlotPoint::new(-2.0, t_total / 2.0),
                            format!("tau = {:.1}", tau_home),
                        )
                        .color(point(&ctx)),
                    );
                    plot_ui.text(
                        Text::new(
                            egui_plot::PlotPoint::new(3.5, t_total / 2.0),
                            format!("tau = {:.1}", tau_travel),
                        )
                        .color(plane(&ctx)),
                    );
                }
            });
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        // === Journey Parameters ===
        group_header(ui, "Journey");

        ui.horizontal(|ui| {
            ui.label("Travel speed:");
            ui.add(
                egui::Slider::new(&mut self.travel_velocity, 0.1..=0.99)
                    .fixed_decimals(2)
                    .suffix("c"),
            );
        });

        let gamma = 1.0 / (1.0 - self.travel_velocity * self.travel_velocity).sqrt();
        value_display(ui, "gamma", gamma, 3);

        ui.horizontal(|ui| {
            ui.label("Total time:");
            ui.add(
                egui::Slider::new(&mut self.total_time, 4.0..=20.0)
                    .fixed_decimals(1)
                    .suffix(" units"),
            );
        });

        // === Animation ===
        section_separator(ui, Some("Animation"));

        ui.horizontal(|ui| {
            ui.label("Progress:");
            ui.add(egui::Slider::new(&mut self.journey_progress, 0.0..=1.0).fixed_decimals(2));
        });

        animation_controls(ui, &mut self.animation);
        progress_slider(ui, &mut self.animation);

        // === Proper Time Results ===
        section_separator(ui, Some("Proper Time"));

        let (tau_home, tau_travel) = self.compute_proper_times();

        info_box(
            ui,
            &format!(
                "Stay-at-home: tau = {:.2}\nTraveler: tau = {:.2}\n\nDifference: {:.2} ({:.0}% younger!)",
                tau_home,
                tau_travel,
                tau_home - tau_travel,
                100.0 * (1.0 - tau_travel / tau_home)
            ),
        );

        ui.add_space(spacing::XS);
        ui.label("The traveler ages LESS because:");
        info_box(
            ui,
            "- Proper time = path length in spacetime\n\
             - Straight worldline = LONGEST proper time\n\
             - This is opposite to Euclidean geometry!\n\
             - tau_travel = T/gamma < T = tau_home",
        );

        // === Display Options ===
        section_separator(ui, Some("Display"));

        ui.horizontal(|ui| {
            ui.checkbox(&mut self.show_grid, "Grid");
            ui.checkbox(&mut self.show_light_cones, "Light cones");
        });
        ui.horizontal(|ui| {
            ui.checkbox(&mut self.show_simultaneity, "Simultaneity lines");
            ui.checkbox(&mut self.show_proper_time_ticks, "Proper time ticks");
        });

        if self.show_proper_time_ticks {
            ui.horizontal(|ui| {
                ui.label("Tick count:");
                ui.add(egui::Slider::new(&mut self.tick_count, 5..=20));
            });
        }
    }

    fn info(&self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        let (tau_home, tau_travel) = self.compute_proper_times();
        let gamma = 1.0 / (1.0 - self.travel_velocity * self.travel_velocity).sqrt();

        ui.horizontal(|ui| {
            ui.label("Twin Paradox");
            ui.separator();
            ui.label(format!(
                "v = {:.2}c, gamma = {:.2}",
                self.travel_velocity, gamma
            ));
            ui.separator();
            ui.colored_label(point(&ctx), format!("Home: tau={:.1}", tau_home));
            ui.separator();
            ui.colored_label(plane(&ctx), format!("Travel: tau={:.1}", tau_travel));
        });
    }

    fn educational_content(&self) -> Option<EducationalContent> {
        Some(MINKOWSKI2_DILATION_EDUCATION)
    }
}

/// Educational content for the time dilation visualization.
const MINKOWSKI2_DILATION_EDUCATION: EducationalContent = EducationalContent {
    title: "Time Dilation and the Twin Paradox",

    overview: "\
This visualization demonstrates the famous TWIN PARADOX from special relativity.

Two twins start at the same point in spacetime. One stays home (vertical worldline), \
while the other travels away at high speed and returns. When they reunite, the \
traveling twin has aged LESS!

This is not a paradox - it's a consequence of the geometry of spacetime. The key \
insight is that PROPER TIME (experienced time) is a PATH LENGTH in spacetime, and \
the straight path has the LONGEST proper time.",

    math_background: "\
Proper time along a worldline is:

    tau = integral sqrt(dt^2 - dx^2) = integral sqrt(1 - v^2) dt

For the STAY-AT-HOME twin (v = 0):
    tau_home = T (total coordinate time)

For the TRAVELER (speed v, round trip):
    tau_travel = T/gamma = Tsqrt(1 - v^2)

Since gamma > 1 for any v > 0, we have:
    tau_travel < tau_home

The traveler ages less!

Example with v = 0.8c:
    gamma = 1/sqrt(1 - 0.64) = 1/sqrt0.36 = 5/3 ~= 1.67
    If home twin ages 10 years, traveler ages only 6 years!",

    how_to_use: "\
- Adjust travel speed to see how time dilation changes
- The tick marks show equal intervals of PROPER TIME
- Count ticks: traveler has FEWER ticks!
- Dashed lines show traveler's 'now' (simultaneity)
- Watch simultaneity lines jump at turnaround
- Click Play to animate the journey",

    key_concepts: "\
- PROPER TIME: Experienced time along a worldline
- Proper time = path length in spacetime
- STRAIGHT path has LONGEST proper time (unlike Euclidean!)
- Time dilation factor: tau = T/gamma
- NOT a paradox: acceleration breaks the symmetry
- The traveler feels acceleration; the stay-at-home doesn't
- Simultaneity lines show the asymmetry clearly",

    resources: &[
        (
            "Twin Paradox - Wikipedia",
            "https://en.wikipedia.org/wiki/Twin_paradox",
        ),
        (
            "Time Dilation",
            "https://en.wikipedia.org/wiki/Time_dilation",
        ),
        ("Proper Time", "https://en.wikipedia.org/wiki/Proper_time"),
    ],
};
