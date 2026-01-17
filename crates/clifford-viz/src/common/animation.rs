//! Animation utilities for continuous playback and timing.
//!
//! Provides a simple [`Animation`] state machine and UI controls for play/pause,
//! speed adjustment, and reset functionality.

use std::f32::consts::TAU;

/// Animation state for continuous playback.
///
/// Tracks time, speed, and loop duration for animations that cycle continuously.
///
/// # Example
/// ```ignore
/// let mut anim = Animation::default();
/// anim.playing = true;
///
/// // In update loop:
/// anim.update(dt);
/// let angle = anim.angle(); // 0 to 2π over loop_duration
/// ```
#[derive(Debug, Clone)]
pub struct Animation {
    /// Whether the animation is currently playing.
    pub playing: bool,
    /// Current time in seconds within the loop.
    pub time: f32,
    /// Playback speed multiplier (1.0 = normal speed).
    pub speed: f32,
    /// Duration of one complete loop in seconds.
    pub loop_duration: f32,
}

impl Default for Animation {
    fn default() -> Self {
        Self {
            playing: false,
            time: 0.0,
            speed: 1.0,
            loop_duration: 4.0,
        }
    }
}

impl Animation {
    /// Create a new animation with specified loop duration.
    #[must_use]
    pub fn with_duration(loop_duration: f32) -> Self {
        Self {
            loop_duration,
            ..Default::default()
        }
    }

    /// Update the animation state by the given delta time.
    ///
    /// If playing, advances time and wraps around at `loop_duration`.
    pub fn update(&mut self, dt: f32) {
        if self.playing {
            self.time += dt * self.speed;
            if self.time > self.loop_duration {
                self.time -= self.loop_duration;
            }
            if self.time < 0.0 {
                self.time += self.loop_duration;
            }
        }
    }

    /// Get normalized progress (0.0 to 1.0) through the animation loop.
    #[must_use]
    pub fn progress(&self) -> f32 {
        self.time / self.loop_duration
    }

    /// Get progress as an angle (0 to 2π) for rotation animations.
    #[must_use]
    pub fn angle(&self) -> f32 {
        self.progress() * TAU
    }

    /// Toggle between playing and paused states.
    pub fn toggle(&mut self) {
        self.playing = !self.playing;
    }

    /// Reset animation to the beginning.
    pub fn reset(&mut self) {
        self.time = 0.0;
    }

    /// Set time directly (useful for manual scrubbing).
    pub fn set_progress(&mut self, progress: f32) {
        self.time = progress.clamp(0.0, 1.0) * self.loop_duration;
    }
}

/// Render animation controls (play/pause, reset, speed slider).
///
/// # Example
/// ```ignore
/// animation_controls(ui, &mut my_animation);
/// ```
pub fn animation_controls(ui: &mut egui::Ui, anim: &mut Animation) {
    ui.horizontal(|ui| {
        // Use separate Play and Pause buttons to avoid toggle issues
        if anim.playing {
            if ui.button("\u{23f8} Pause").clicked() {
                anim.playing = false;
            }
        } else if ui.button("\u{25b6} Play").clicked() {
            anim.playing = true;
        }
        if ui.button("\u{23ee} Reset").clicked() {
            anim.time = 0.0;
        }
        ui.add(
            egui::Slider::new(&mut anim.speed, 0.1..=3.0)
                .text("Speed")
                .fixed_decimals(1),
        );
    });
}

/// Render a progress slider for manual animation scrubbing.
///
/// Returns true if the user is actively dragging the slider.
pub fn progress_slider(ui: &mut egui::Ui, anim: &mut Animation) -> bool {
    let mut progress = anim.progress();
    let response = ui.add(
        egui::Slider::new(&mut progress, 0.0..=1.0)
            .text("Progress")
            .fixed_decimals(2),
    );
    // Only update when user is actively dragging, not on any change
    // (fixed_decimals can cause spurious "changed" events due to rounding)
    if response.dragged() {
        anim.set_progress(progress);
    }
    response.dragged()
}

/// Easing functions for smooth animations.
pub mod easing {
    /// Linear interpolation (no easing).
    #[must_use]
    pub fn linear(t: f32) -> f32 {
        t
    }

    /// Ease in (slow start).
    #[must_use]
    pub fn ease_in_quad(t: f32) -> f32 {
        t * t
    }

    /// Ease out (slow end).
    #[must_use]
    pub fn ease_out_quad(t: f32) -> f32 {
        1.0 - (1.0 - t) * (1.0 - t)
    }

    /// Ease in and out (slow start and end).
    #[must_use]
    pub fn ease_in_out_quad(t: f32) -> f32 {
        if t < 0.5 {
            2.0 * t * t
        } else {
            1.0 - (-2.0 * t + 2.0).powi(2) / 2.0
        }
    }

    /// Smooth step (Hermite interpolation).
    #[must_use]
    pub fn smooth_step(t: f32) -> f32 {
        t * t * (3.0 - 2.0 * t)
    }

    /// Smoother step (Ken Perlin's improved version).
    #[must_use]
    pub fn smoother_step(t: f32) -> f32 {
        t * t * t * (t * (t * 6.0 - 15.0) + 10.0)
    }
}
