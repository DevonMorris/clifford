//! 1+1D Spacetime Diagram Visualization
//!
//! This demo visualizes Minkowski spacetime in 1+1 dimensions (one space + one time),
//! demonstrating the causal structure of special relativity:
//!
//! - Light cones at 45 degrees (c = 1)
//! - Timelike, spacelike, and lightlike intervals
//! - Lorentz boost transformations
//! - Spacetime events and worldlines
//!
//! ## Mathematical Background
//!
//! In Minkowski spacetime, the interval between two events is:
//!
//! ```text
//! Δs² = Δt² - Δx²  (with c = 1)
//! ```
//!
//! This interval is invariant under Lorentz transformations:
//! - Δs² > 0: timelike (events can be causally connected)
//! - Δs² < 0: spacelike (events cannot be causally connected)
//! - Δs² = 0: lightlike (connected by light rays)
//!
//! A Lorentz boost with velocity v transforms coordinates:
//!
//! ```text
//! t' = γ(t - vx)
//! x' = γ(x - vt)
//! where γ = 1/√(1 - v²)
//! ```

use crate::common::prelude::*;
use clifford::specialized::minkowski::dim2::Vector;
use egui::Color32;
use egui_plot::{Line, Plot, PlotPoints, Points};

/// Fixed viewport bounds for stable viewing experience.
const VIEWPORT_BOUNDS: f64 = 6.0;

/// A spacetime event (a point in spacetime).
///
/// This stores the spacetime position using the clifford library's [`Vector`] type
/// from the Minkowski Cl(1,1) algebra, along with a label for visualization.
#[derive(Clone, Debug)]
struct SpacetimeEvent {
    /// Position in spacetime as a Minkowski vector.
    /// Uses clifford's Vector<f64> with spacelike `x` (e₁² = -1) and timelike `t` (e₂² = +1).
    /// This gives the physics (timelike-positive) convention: v² = t² - x².
    position: Vector<f64>,
    /// Label for the event.
    label: String,
}

impl SpacetimeEvent {
    /// Creates a new spacetime event at coordinates (x, t).
    fn new(x: f64, t: f64, label: impl Into<String>) -> Self {
        Self {
            position: Vector::new(x, t),
            label: label.into(),
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

    /// Sets the spatial coordinate.
    fn set_x(&mut self, x: f64) {
        self.position = Vector::new(x, self.position.t());
    }

    /// Sets the time coordinate.
    fn set_t(&mut self, t: f64) {
        self.position = Vector::new(self.position.x(), t);
    }

    /// Apply a Lorentz boost with velocity v (in units of c).
    ///
    /// The Lorentz transformation preserves the spacetime interval:
    /// - t' = γ(t - vx)
    /// - x' = γ(x - vt)
    ///
    /// where γ = 1/√(1 - v²) is the Lorentz factor.
    fn boost(&self, v: f64) -> Self {
        let gamma = 1.0 / (1.0 - v * v).sqrt();
        let x = self.position.x();
        let t = self.position.t();
        Self {
            position: Vector::new(gamma * (x - v * t), gamma * (t - v * x)),
            label: self.label.clone(),
        }
    }

    /// Compute the spacetime interval squared to another event.
    ///
    /// Uses the clifford library's native norm_squared() which computes Δs² = Δt² - Δx²
    /// (physics/timelike-positive convention).
    fn interval_squared(&self, other: &SpacetimeEvent) -> f64 {
        let delta = Vector::new(
            other.position.x() - self.position.x(),
            other.position.t() - self.position.t(),
        );
        delta.norm_squared()
    }

    /// Compute the proper time between events (if timelike).
    ///
    /// Returns `Some(τ)` where τ = √(Δs²) for timelike intervals (Δs² > 0).
    /// Returns `None` for spacelike intervals (no proper time exists).
    fn proper_time_to(&self, other: &SpacetimeEvent) -> Option<f64> {
        let s2 = self.interval_squared(other);
        if s2 > 0.0 {
            Some(s2.sqrt())
        } else {
            None // Spacelike interval has no proper time
        }
    }
}

/// Classification of spacetime intervals.
#[derive(Clone, Copy, Debug, PartialEq)]
enum IntervalType {
    /// Δs² > 0: events can be causally connected.
    Timelike,
    /// Δs² < 0: events cannot be causally connected.
    Spacelike,
    /// Δs² = 0: events on the light cone.
    Lightlike,
}

impl IntervalType {
    /// Determines interval type from squared interval.
    fn from_interval_squared(s2: f64) -> Self {
        const EPSILON: f64 = 1e-6;
        if s2 > EPSILON {
            Self::Timelike
        } else if s2 < -EPSILON {
            Self::Spacelike
        } else {
            Self::Lightlike
        }
    }

    /// Human-readable description.
    fn description(&self) -> &'static str {
        match self {
            Self::Timelike => "timelike (causally connected)",
            Self::Spacelike => "spacelike (causally disconnected)",
            Self::Lightlike => "lightlike (on light cone)",
        }
    }
}

/// A worldline (trajectory through spacetime).
#[derive(Clone, Debug)]
struct Worldline {
    /// Starting event.
    start: SpacetimeEvent,
    /// Ending event.
    end: SpacetimeEvent,
    /// Color for display.
    color: Color32,
}

impl Worldline {
    /// Apply a Lorentz boost to both endpoints.
    fn boost(&self, v: f64) -> Self {
        Self {
            start: self.start.boost(v),
            end: self.end.boost(v),
            color: self.color,
        }
    }
}

/// Demo state for 1+1D spacetime diagram visualization.
pub struct Minkowski2DiagramDemo {
    /// List of events in the diagram.
    events: Vec<SpacetimeEvent>,
    /// List of worldlines.
    worldlines: Vec<Worldline>,
    /// Boost velocity (as fraction of c).
    boost_velocity: f32,
    /// Whether to show light cones from origin.
    show_light_cones: bool,
    /// Whether to show coordinate grid.
    show_grid: bool,
    /// Whether to show boosted frame axes.
    show_boosted_axes: bool,
    /// Selected event index (for interval calculation).
    selected_event: Option<usize>,
    /// Second selected event (for interval display).
    second_event: Option<usize>,
    /// Animation state for boost.
    animation: Animation,
}

impl Default for Minkowski2DiagramDemo {
    fn default() -> Self {
        // Default events demonstrating different interval types
        let events = vec![
            SpacetimeEvent::new(0.0, 0.0, "O"),
            SpacetimeEvent::new(0.0, 2.0, "A"),
            SpacetimeEvent::new(1.5, 3.0, "B"),
            SpacetimeEvent::new(3.0, 1.0, "C"),
            SpacetimeEvent::new(2.0, 2.0, "D"),
        ];

        Self {
            events,
            worldlines: Vec::new(),
            boost_velocity: 0.0,
            show_light_cones: true,
            show_grid: true,
            show_boosted_axes: true,
            selected_event: Some(0),
            second_event: Some(2),
            animation: Animation::with_duration(8.0),
        }
    }
}

impl Minkowski2DiagramDemo {
    /// Get events, applying boost if non-zero.
    fn boosted_events(&self) -> Vec<SpacetimeEvent> {
        if self.boost_velocity.abs() < 0.001 {
            self.events.clone()
        } else {
            self.events
                .iter()
                .map(|e| e.boost(f64::from(self.boost_velocity)))
                .collect()
        }
    }

    /// Get worldlines, applying boost if non-zero.
    fn boosted_worldlines(&self) -> Vec<Worldline> {
        if self.boost_velocity.abs() < 0.001 {
            self.worldlines.clone()
        } else {
            self.worldlines
                .iter()
                .map(|w| w.boost(f64::from(self.boost_velocity)))
                .collect()
        }
    }

    /// Draw light cone from origin.
    fn draw_light_cone(&self, plot_ui: &mut egui_plot::PlotUi, ctx: &egui::Context, extent: f64) {
        // Use circle color (gold) for light cones
        let light_color = with_alpha(circle(ctx), 150);

        // Future light cone (45° lines going up)
        let future_right: PlotPoints = (0..=50)
            .map(|i| {
                let t = extent * (i as f64) / 50.0;
                [t, t] // x = t (rightward light ray)
            })
            .collect();
        let future_left: PlotPoints = (0..=50)
            .map(|i| {
                let t = extent * (i as f64) / 50.0;
                [-t, t] // x = -t (leftward light ray)
            })
            .collect();

        // Past light cone (45° lines going down)
        let past_right: PlotPoints = (0..=50)
            .map(|i| {
                let t = -extent * (i as f64) / 50.0;
                [-t, t] // x = -t, t < 0
            })
            .collect();
        let past_left: PlotPoints = (0..=50)
            .map(|i| {
                let t = -extent * (i as f64) / 50.0;
                [t, t] // x = t, t < 0
            })
            .collect();

        plot_ui.line(
            Line::new(future_right)
                .color(light_color)
                .width(2.0)
                .name("Future light cone"),
        );
        plot_ui.line(
            Line::new(future_left)
                .color(light_color)
                .width(2.0)
                .name("Future light cone"),
        );
        plot_ui.line(
            Line::new(past_right)
                .color(light_color)
                .width(2.0)
                .name("Past light cone"),
        );
        plot_ui.line(
            Line::new(past_left)
                .color(light_color)
                .width(2.0)
                .name("Past light cone"),
        );

        // Shade the light cone regions with very light fill
        // Future timelike region
        let future_fill: PlotPoints =
            vec![[0.0, 0.0], [-extent, extent], [extent, extent], [0.0, 0.0]]
                .into_iter()
                .collect();
        plot_ui.line(
            Line::new(future_fill)
                .color(with_alpha(light_color, 20))
                .fill(0.0),
        );
    }

    /// Draw boosted frame axes (t' and x').
    fn draw_boosted_axes(
        &self,
        plot_ui: &mut egui_plot::PlotUi,
        ctx: &egui::Context,
        v: f64,
        extent: f64,
    ) {
        if v.abs() < 0.001 {
            return;
        }

        // The t' axis has slope 1/v (events at same x' have x = vt)
        // The x' axis has slope v (events at same t' have t = vx)

        let axis_color = with_alpha(motor(ctx), 150);

        // t' axis (worldline of boosted origin)
        // Parametrize by t: x = vt
        let t_prime_axis: PlotPoints = (-50..=50)
            .map(|i| {
                let t = extent * (i as f64) / 50.0;
                let x = v * t;
                [x, t]
            })
            .collect();

        // x' axis (line of simultaneity in boosted frame)
        // Parametrize by x: t = vx
        let x_prime_axis: PlotPoints = (-50..=50)
            .map(|i| {
                let x = extent * (i as f64) / 50.0;
                let t = v * x;
                [x, t]
            })
            .collect();

        plot_ui.line(
            Line::new(t_prime_axis)
                .color(axis_color)
                .width(1.5)
                .style(egui_plot::LineStyle::Dashed { length: 8.0 })
                .name("t' axis"),
        );
        plot_ui.line(
            Line::new(x_prime_axis)
                .color(axis_color)
                .width(1.5)
                .style(egui_plot::LineStyle::Dashed { length: 8.0 })
                .name("x' axis"),
        );
    }
}

impl VisualizationApp for Minkowski2DiagramDemo {
    fn name(&self) -> &'static str {
        "Minkowski 1+1D Spacetime"
    }

    fn update(&mut self, dt: f32) {
        self.animation.update(dt);
        if self.animation.playing {
            // Oscillate boost velocity between -0.8 and 0.8
            let t = self.animation.progress();
            self.boost_velocity = 0.7 * (t * std::f32::consts::TAU).sin();
        }
    }

    fn render(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        let boosted_events = self.boosted_events();
        let boosted_worldlines = self.boosted_worldlines();

        Plot::new("minkowski2_diagram")
            .data_aspect(1.0)
            .show_axes(false)
            .show_grid(false)
            .allow_zoom(true)
            .allow_drag(true)
            .x_axis_label("x (space)")
            .y_axis_label("t (time)")
            .show(ui, |plot_ui| {
                // Draw coordinate grid
                if self.show_grid {
                    for l in grid_2d(&ctx, VIEWPORT_BOUNDS, 1.0) {
                        plot_ui.line(l);
                    }
                    // Custom axes with labels (renamed to avoid shadowing color functions)
                    // x-axis (space)
                    let x_axis_line: PlotPoints =
                        vec![[-VIEWPORT_BOUNDS, 0.0], [VIEWPORT_BOUNDS, 0.0]].into();
                    plot_ui.line(
                        Line::new(x_axis_line)
                            .color(x_axis(&ctx))
                            .width(2.0)
                            .name("x (space)"),
                    );
                    // t-axis (time)
                    let t_axis_line: PlotPoints =
                        vec![[0.0, -VIEWPORT_BOUNDS], [0.0, VIEWPORT_BOUNDS]].into();
                    plot_ui.line(
                        Line::new(t_axis_line)
                            .color(t_axis(&ctx))
                            .width(2.0)
                            .name("t (time)"),
                    );
                }

                // Draw light cones from origin
                if self.show_light_cones {
                    self.draw_light_cone(plot_ui, &ctx, VIEWPORT_BOUNDS);
                }

                // Draw boosted frame axes
                if self.show_boosted_axes {
                    self.draw_boosted_axes(
                        plot_ui,
                        &ctx,
                        f64::from(self.boost_velocity),
                        VIEWPORT_BOUNDS,
                    );
                }

                // Draw worldlines
                for worldline in &boosted_worldlines {
                    let points: PlotPoints = vec![
                        [worldline.start.x(), worldline.start.t()],
                        [worldline.end.x(), worldline.end.t()],
                    ]
                    .into();
                    plot_ui.line(Line::new(points).color(worldline.color).width(2.0));
                }

                // Draw interval between selected events
                if let (Some(i), Some(j)) = (self.selected_event, self.second_event) {
                    if i < boosted_events.len() && j < boosted_events.len() {
                        let e1 = &boosted_events[i];
                        let e2 = &boosted_events[j];
                        let s2 = e1.interval_squared(e2);
                        let interval_type = IntervalType::from_interval_squared(s2);

                        let interval_color = match interval_type {
                            IntervalType::Timelike => plane(&ctx),
                            IntervalType::Spacelike => line_secondary(&ctx),
                            IntervalType::Lightlike => circle(&ctx),
                        };

                        let points: PlotPoints = vec![[e1.x(), e1.t()], [e2.x(), e2.t()]].into();
                        plot_ui.line(
                            Line::new(points)
                                .color(interval_color)
                                .width(2.5)
                                .style(egui_plot::LineStyle::Dashed { length: 6.0 })
                                .name("Interval"),
                        );
                    }
                }

                // Draw events
                for (i, event) in boosted_events.iter().enumerate() {
                    let is_selected =
                        self.selected_event == Some(i) || self.second_event == Some(i);
                    let event_color = if is_selected {
                        active(&ctx)
                    } else {
                        point(&ctx)
                    };
                    let radius = if is_selected { 8.0 } else { 6.0 };

                    plot_ui.points(
                        Points::new(vec![[event.x(), event.t()]])
                            .color(event_color)
                            .radius(radius)
                            .filled(true)
                            .name(&event.label),
                    );

                    // Label the event
                    plot_ui.text(
                        egui_plot::Text::new(
                            egui_plot::PlotPoint::new(event.x() + 0.15, event.t() + 0.15),
                            event.label.clone(),
                        )
                        .color(text_primary(&ctx)),
                    );
                }
            });
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        // === Lorentz Boost ===
        group_header(ui, "Lorentz Boost");

        ui.horizontal(|ui| {
            ui.label("Velocity v/c:");
            ui.add(
                egui::Slider::new(&mut self.boost_velocity, -0.9..=0.9)
                    .fixed_decimals(2)
                    .suffix("c"),
            );
        });

        let gamma = 1.0 / (1.0 - self.boost_velocity * self.boost_velocity).sqrt();
        value_display(ui, "\u{03b3}", gamma, 3);

        ui.add_space(spacing::XS);
        animation_controls(ui, &mut self.animation);
        progress_slider(ui, &mut self.animation);

        // === Interval Calculation ===
        section_separator(ui, Some("Interval"));

        // Event selection dropdowns
        ui.horizontal(|ui| {
            ui.label("From:");
            egui::ComboBox::from_id_salt("event1")
                .selected_text(
                    self.selected_event
                        .and_then(|i| self.events.get(i))
                        .map_or("None", |e| &e.label),
                )
                .show_ui(ui, |ui| {
                    for (i, event) in self.events.iter().enumerate() {
                        ui.selectable_value(&mut self.selected_event, Some(i), &event.label);
                    }
                });
        });

        ui.horizontal(|ui| {
            ui.label("To:");
            egui::ComboBox::from_id_salt("event2")
                .selected_text(
                    self.second_event
                        .and_then(|i| self.events.get(i))
                        .map_or("None", |e| &e.label),
                )
                .show_ui(ui, |ui| {
                    for (i, event) in self.events.iter().enumerate() {
                        ui.selectable_value(&mut self.second_event, Some(i), &event.label);
                    }
                });
        });

        // Show interval information
        if let (Some(i), Some(j)) = (self.selected_event, self.second_event) {
            if i < self.events.len() && j < self.events.len() && i != j {
                // Use original (unboosted) events for invariant calculation
                let e1 = &self.events[i];
                let e2 = &self.events[j];
                let s2 = e1.interval_squared(e2);
                let interval_type = IntervalType::from_interval_squared(s2);

                ui.add_space(spacing::XS);
                value_display(ui, "\u{0394}s\u{00b2}", s2 as f32, 3);
                ui.label(format!("Type: {}", interval_type.description()));

                if let Some(tau) = e1.proper_time_to(e2) {
                    value_display(ui, "\u{03c4}", tau as f32, 3);
                }

                info_box(
                    ui,
                    "\u{0394}s\u{00b2} is invariant under boosts!\nTry changing v and watch it stay constant.",
                );
            }
        }

        // === Events Editor ===
        section_separator(ui, Some("Events"));

        let mut events_to_remove = Vec::new();
        let mut event_updates: Vec<(usize, f64, f64)> = Vec::new();

        for (i, event) in self.events.iter().enumerate() {
            let mut x = event.x();
            let mut t = event.t();
            let mut changed = false;

            ui.horizontal(|ui| {
                ui.label(&event.label);
                if ui
                    .add(egui::DragValue::new(&mut x).prefix("x:").speed(0.1))
                    .changed()
                {
                    changed = true;
                }
                if ui
                    .add(egui::DragValue::new(&mut t).prefix("t:").speed(0.1))
                    .changed()
                {
                    changed = true;
                }
                if ui.small_button("\u{2715}").clicked() {
                    events_to_remove.push(i);
                }
            });

            if changed {
                event_updates.push((i, x, t));
            }
        }

        // Apply coordinate updates
        for (i, x, t) in event_updates {
            self.events[i].set_x(x);
            self.events[i].set_t(t);
        }

        // Remove events (in reverse order to preserve indices)
        for i in events_to_remove.into_iter().rev() {
            self.events.remove(i);
            // Adjust selection indices
            if self.selected_event == Some(i) {
                self.selected_event = None;
            } else if let Some(sel) = self.selected_event {
                if sel > i {
                    self.selected_event = Some(sel - 1);
                }
            }
            if self.second_event == Some(i) {
                self.second_event = None;
            } else if let Some(sel) = self.second_event {
                if sel > i {
                    self.second_event = Some(sel - 1);
                }
            }
        }

        if ui.button("Add Event").clicked() {
            let label = format!("E{}", self.events.len());
            self.events.push(SpacetimeEvent::new(1.0, 1.0, label));
        }

        // === Display Options ===
        section_separator(ui, Some("Display"));

        ui.horizontal(|ui| {
            ui.checkbox(&mut self.show_grid, "Grid");
            ui.checkbox(&mut self.show_light_cones, "Light cones");
        });
        ui.checkbox(&mut self.show_boosted_axes, "Boosted axes (t', x')");
    }

    fn info(&self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        ui.horizontal(|ui| {
            ui.label("1+1D Spacetime (c = 1)");
            ui.separator();
            ui.label(format!(
                "v = {:.2}c, \u{03b3} = {:.2}",
                self.boost_velocity,
                1.0 / (1.0 - self.boost_velocity * self.boost_velocity).sqrt()
            ));
            ui.separator();
            ui.colored_label(circle(&ctx), "Light cone");
            ui.separator();
            ui.colored_label(plane(&ctx), "Timelike");
            ui.separator();
            ui.colored_label(line_secondary(&ctx), "Spacelike");
        });
    }

    fn educational_content(&self) -> Option<EducationalContent> {
        Some(MINKOWSKI2_DIAGRAM_EDUCATION)
    }
}

/// Educational content for the spacetime diagram visualization.
const MINKOWSKI2_DIAGRAM_EDUCATION: EducationalContent = EducationalContent {
    title: "Minkowski Spacetime in 1+1D",

    overview: "\
This visualization demonstrates the structure of MINKOWSKI SPACETIME in 1+1 dimensions \
(one spatial dimension + time). This is the geometry of special relativity.

The key features are:
\u{2022} Light travels at 45\u{00b0} (we use units where c = 1)
\u{2022} The LIGHT CONE divides spacetime into causally connected and disconnected regions
\u{2022} LORENTZ BOOSTS transform between reference frames moving at different velocities
\u{2022} The SPACETIME INTERVAL \u{0394}s\u{00b2} is INVARIANT under boosts",

    math_background: "\
The spacetime interval between two events is:

    \u{0394}s\u{00b2} = \u{0394}t\u{00b2} - \u{0394}x\u{00b2}   (with c = 1)

This is INVARIANT under Lorentz transformations!

Interval types:
  \u{2022} \u{0394}s\u{00b2} > 0: TIMELIKE (inside light cone)
      - Events CAN be causally connected
      - \u{03c4} = \u{221a}(\u{0394}s\u{00b2}) is the proper time
  \u{2022} \u{0394}s\u{00b2} < 0: SPACELIKE (outside light cone)
      - Events CANNOT be causally connected
  \u{2022} \u{0394}s\u{00b2} = 0: LIGHTLIKE (on the light cone)
      - Events connected by light rays

Lorentz boost with velocity v:
    t' = \u{03b3}(t - vx)
    x' = \u{03b3}(x - vt)
    where \u{03b3} = 1/\u{221a}(1 - v\u{00b2})",

    how_to_use: "\
\u{2022} Use the velocity slider to apply a Lorentz boost
\u{2022} Watch events move but the INTERVAL stays constant!
\u{2022} The light cone (45\u{00b0} lines) is the same in all frames
\u{2022} Select two events to see their spacetime interval
\u{2022} Dashed lines show the boosted coordinate axes (t', x')
\u{2022} Click 'Play' to animate continuous boosts",

    key_concepts: "\
\u{2022} LIGHT CONE: 45\u{00b0} lines from origin, boundary of causality
\u{2022} TIMELIKE: Inside cone, accessible, has proper time
\u{2022} SPACELIKE: Outside cone, inaccessible, no causal connection
\u{2022} INVARIANT INTERVAL: \u{0394}s\u{00b2} = \u{0394}t\u{00b2} - \u{0394}x\u{00b2} same in all frames
\u{2022} LORENTZ BOOST: Transforms coordinates, preserves interval
\u{2022} The factor \u{03b3} = 1/\u{221a}(1-v\u{00b2}) appears in time dilation",

    resources: &[
        (
            "Minkowski Space - Wikipedia",
            "https://en.wikipedia.org/wiki/Minkowski_space",
        ),
        (
            "Spacetime Diagrams",
            "https://en.wikipedia.org/wiki/Spacetime_diagram",
        ),
        (
            "Lorentz Transformation",
            "https://en.wikipedia.org/wiki/Lorentz_transformation",
        ),
    ],
};
