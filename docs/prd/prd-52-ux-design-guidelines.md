# PRD-52: UX & Graphic Design Guidelines for Visualizations

**Status**: Draft
**Goal**: Transform visualization demos from "engineer brain" to "designer brain" quality
**Depends on**: PRD-48.1 (visualization framework)

## Problem Statement

The current visualization demos are functional but suffer from several UX and graphic design issues:

1. **High-contrast color palette** - Pure reds against pure blacks, saturated primaries that clash
2. **Dynamic viewport sizing** - Bounds calculated from entities causes disorienting viewport jumps
3. **Visual hierarchy issues** - Everything draws with equal visual weight
4. **Lack of visual polish** - Missing subtle cues like shadows, gradients, anti-aliased rendering
5. **Control panel UX** - Controls are functional but don't guide the user through a learning journey
6. **Missing visual feedback** - No hover states, no transition animations, no confirmation of actions

These issues reflect "engineer brain" thinking: the visualization technically works, so it's done. A professional designer would ask: does it *feel* good to use? Does it invite exploration? Is the visual hierarchy clear?

## Design Philosophy

### 1. Visual Hierarchy Through Contrast

Not all elements are equally important. Create hierarchy through:

- **Size** - Important elements are larger
- **Saturation** - Primary focus is saturated, secondary is desaturated
- **Opacity** - Supporting elements fade into background
- **Line weight** - Primary objects have thicker strokes

**Current Problem:**
```
Grid lines:       ████ (66, 66, 66)     - Dark gray, visible
Axes:             ████ (244, 67, 54)    - Pure saturated red
Points:           ████ (66, 133, 244)   - Saturated blue
```

Everything is highly saturated, nothing "pops" because everything tries to.

**Solution:** Desaturate supporting elements, reserve saturation for focus.

### 2. Color Palette Redesign

#### Current Palette Problems

| Element | Current Color | Issue |
|---------|---------------|-------|
| Background | `rgb(30, 30, 30)` | Too dark, harsh contrast with saturated colors |
| Grid | `rgb(66, 66, 66)` | Visible but fights with other elements |
| X Axis | `rgb(244, 67, 54)` | Pure Material red, too aggressive |
| Y Axis | `rgb(76, 175, 80)` | Pure Material green, clashes with red |
| Point | `rgb(66, 133, 244)` | Google blue, okay but saturated |
| Line | `rgb(234, 67, 53)` | Same as X axis, confusing |
| Rotor | `rgb(171, 71, 188)` | Neon purple, distracting |

#### Proposed Palette: "Soft Technical"

**Design Principles:**
- Warm neutral background (not pure gray)
- Desaturated versions of semantic colors
- Primary colors only for interactive elements and focus
- Use opacity rather than different grays for hierarchy

**Base Colors:**

| Role | Hex | RGB | Rationale |
|------|-----|-----|-----------|
| Background | `#1e1e22` | `(30, 30, 34)` | Slightly warm, easier on eyes |
| Surface | `#252529` | `(37, 37, 41)` | For cards, panels |
| Grid Major | `#3a3a40` | `(58, 58, 64)` | Subtle, doesn't compete |
| Grid Minor | `#2e2e33` | `(46, 46, 51)` | Very subtle |
| Text Primary | `#e0e0e0` | `(224, 224, 224)` | Off-white, less harsh |
| Text Secondary | `#9e9e9e` | `(158, 158, 158)` | For labels |

**Semantic Colors (Desaturated):**

| Role | Current | Proposed | Rationale |
|------|---------|----------|-----------|
| X Axis | `#f44336` | `#b85450` | Muted brick red |
| Y Axis | `#4caf50` | `#5a9a5e` | Muted sage green |
| Z Axis | `#2196f3` | `#5085b0` | Muted steel blue |
| Point (primary) | `#4285f4` | `#6a9bd4` | Softer blue |
| Point (secondary) | `#64b5f6` | `#8ab4d0` | Even softer |
| Line | `#ea4335` | `#c47068` | Muted coral |
| Plane | `#34a853` | `#68a878` | Muted mint |
| Rotor/Transform | `#ab47bc` | `#9878a8` | Muted lavender |

**Accent Colors (Full saturation, use sparingly):**

| Role | Hex | When to Use |
|------|-----|-------------|
| Selection | `#ffd54f` | Currently selected element |
| Hover | `#81d4fa` | Mouse-over highlight |
| Active | `#aed581` | Currently animating |

### 3. Fixed Viewport with Optional Auto-Fit

**Problem:** Dynamic bounds cause disorienting jumps when entities move/scale.

**Solution:** Fixed viewport by default, explicit user control for zoom/fit.

```rust
pub struct ViewportConfig {
    /// Fixed bounds (-bounds to +bounds on each axis)
    pub bounds: f64,
    /// Whether to allow user zoom/pan
    pub interactive: bool,
    /// Whether to show a "fit to content" button (not auto-fit!)
    pub show_fit_button: bool,
}

impl Default for ViewportConfig {
    fn default() -> Self {
        Self {
            bounds: 5.0,          // Fixed at ±5 by default
            interactive: true,    // User can zoom/pan
            show_fit_button: true, // Button to fit (not automatic)
        }
    }
}
```

**UI Pattern:**
- Viewport starts at fixed bounds (±5.0)
- User can zoom/pan manually
- "Fit to Content" button (icon: ⊡) available in corner
- Never auto-fit - user maintains control

### 4. Animation & Transitions

**Problem:** State changes are instant, feels jarring.

**Solution:** Smooth transitions for all state changes.

```rust
pub struct Transition<T> {
    current: T,
    target: T,
    progress: f32,    // 0.0 to 1.0
    duration: f32,    // seconds
    easing: Easing,
}

impl<T: Lerp> Transition<T> {
    pub fn set_target(&mut self, target: T) {
        self.current = self.value(); // Start from current position
        self.target = target;
        self.progress = 0.0;
    }

    pub fn update(&mut self, dt: f32) {
        self.progress = (self.progress + dt / self.duration).min(1.0);
    }

    pub fn value(&self) -> T {
        self.current.lerp(&self.target, self.easing.apply(self.progress))
    }
}

pub enum Easing {
    Linear,
    EaseOutQuad,     // Decelerate at end
    EaseInOutCubic,  // Smooth both ends
}
```

**Apply transitions to:**
- Slider value changes (angle, dilation)
- Toggling display options
- Viewport pan/zoom
- Color changes

### 5. Visual Feedback & Affordances

#### Hover States

When cursor hovers over interactive elements:
- Geometric objects: subtle glow/outline
- Control values: highlight background
- Buttons: lighten by 10%

#### Cursor Changes

| Context | Cursor |
|---------|--------|
| Over draggable object | `grab` |
| Dragging object | `grabbing` |
| Over zoomable area | `zoom-in` / `zoom-out` |
| Over pan area | `move` |

#### Tooltips

Rich tooltips for geometric objects:
```
┌─────────────────────────────┐
│ Point v                     │
│ ─────────────────────────── │
│ Position: (2.000, 0.500)    │
│ Magnitude: 2.062            │
│ Angle: 14.0°                │
│ ─────────────────────────── │
│ Drag to move                │
└─────────────────────────────┘
```

### 6. Control Panel UX

#### Current Problems
- All controls dumped in one column
- No visual grouping
- No indication of what controls affect what
- Dense text, hard to scan

#### Proposed Structure

```
┌─────────────────────────────────┐
│ ROTATION                    ▼   │  ← Collapsible section
│ ┌─────────────────────────────┐ │
│ │ Angle θ      [-180°───180°] │ │  ← Primary control
│ │                    45.0°    │ │
│ └─────────────────────────────┘ │
│                                 │
│ ┌─────────────────────────────┐ │
│ │ ▶ Play   ⏮ Reset   1.0x ─┤ │ │  ← Animation controls
│ │ ├────────────●────────────┤ │ │  ← Progress scrubber
│ └─────────────────────────────┘ │
├─────────────────────────────────┤
│ VECTOR                      ▼   │
│ Input v   x: 2.00   y: 0.50     │
│                                 │
│ Output    (1.414, 1.414)        │  ← Read-only display
│ Magnitude 2.062 → 2.062         │
├─────────────────────────────────┤
│ DISPLAY                     ▼   │
│ ☑ Grid      ☑ Axes              │  ← Toggle group
│ ☑ Bivector  ☐ Labels            │
└─────────────────────────────────┘
```

#### Design Principles

1. **Group related controls** with clear section headers
2. **Primary controls first** - angle slider is the main interaction
3. **Compact toggle groups** - checkboxes in rows, not columns
4. **Clear visual hierarchy** - sections separated, labels aligned
5. **Live feedback** - output values update as input changes

### 7. Grid & Axes Redesign

#### Problems
- Grid is uniform, no hierarchy between major/minor
- Axes use same colors as geometric objects (confusing)
- No subtle markers at integer positions

#### Solution

**Three-level grid hierarchy:**
```rust
pub struct GridStyle {
    pub minor_color: Color32,    // Very subtle, every 0.5 units
    pub major_color: Color32,    // Visible, every 1.0 unit
    pub axis_color: Color32,     // Distinct, at x=0 and y=0
    pub minor_width: f32,        // 0.5px
    pub major_width: f32,        // 1.0px
    pub axis_width: f32,         // 1.5px
}
```

**Axis styling:**
- Subtle arrows at positive ends
- Small tick marks at integers (not numbers, too cluttered)
- Origin marker (small circle or crosshair)

### 8. Typography

#### Current Issues
- System font, inconsistent sizing
- Mathematical symbols look out of place
- Dense text blocks

#### Proposed

| Element | Font | Size | Weight |
|---------|------|------|--------|
| Section headers | System | 14px | Bold |
| Control labels | System | 12px | Normal |
| Values/numbers | Monospace | 12px | Normal |
| Math formulas | Monospace | 11px | Normal |
| Tooltips | System | 11px | Normal |

**Unicode math symbols** (consistent use):
- θ for angles (not "theta")
- × for multiplication
- ° for degrees
- ₁₂ for subscripts (e₁₂)
- → for "maps to"

### 9. Mobile/Small Screen Considerations

While primarily a desktop app, handle small windows gracefully:

- **< 800px width**: Stack controls below visualization
- **< 600px width**: Hide non-essential controls, show "More" button
- **Minimum**: 640×480 (current), but ensure usability

### 10. Loading & Performance

- Show loading indicator if first render takes > 100ms
- Debounce rapid slider changes (don't re-render every pixel)
- Precompute static elements (grid lines don't need recalculation)

## Implementation Plan

### Phase 1: Color Palette Update
**Files:** `colors.rs`

1. Add new color constants for the "Soft Technical" palette
2. Keep old colors as `legacy::` for comparison
3. Add helper functions:
   - `desaturate(color, amount)` - reduce saturation
   - `with_opacity(color, opacity)` - apply opacity
   - `highlight(color)` - lighten for hover
4. Update all demos to use new palette

### Phase 2: Fixed Viewport
**Files:** `app.rs`, all examples

1. Add `ViewportConfig` struct
2. Modify `Plot` setup to use fixed bounds by default
3. Add "Fit to Content" button (not auto-fit)
4. Update euclidean2 demo to use fixed bounds

### Phase 3: Grid Redesign
**Files:** `grid.rs`

1. Add `GridStyle` struct
2. Implement three-level hierarchy (minor/major/axis)
3. Add subtle origin marker
4. Add axis arrows and tick marks

### Phase 4: Control Panel Improvements
**Files:** `widgets.rs`, `app.rs`, all examples

1. Add `collapsible_section` with styling
2. Create `toggle_group` for compact checkboxes
3. Improve section separators
4. Reorganize euclidean2 controls into groups

### Phase 5: Transitions & Feedback
**Files:** `animation.rs` (new: `transition.rs`)

1. Add `Transition<T>` struct
2. Add easing functions
3. Apply to slider values in demos
4. Add hover highlighting

### Phase 6: Polish
**Files:** Various

1. Consistent tooltip implementation
2. Cursor changes for interactive areas
3. Performance optimizations (debouncing, caching)
4. Final color tuning based on testing

## Success Criteria

### Before/After Comparison

**Before (Current):**
- Saturated clashing colors
- Viewport jumps when values change
- Flat, uniform visual weight
- Dense, unorganized controls

**After (Goal):**
- Harmonious, desaturated palette
- Stable viewport under user control
- Clear visual hierarchy
- Organized, scannable controls

### User Experience Tests

1. **First impression** - Does it look professional?
2. **Learnability** - Can a new user figure out the controls?
3. **Comfort** - Can you use it for 30+ minutes without eye strain?
4. **Feedback** - Do you know what's happening at all times?

## Appendix A: Color Accessibility

All color combinations must meet WCAG AA contrast requirements:
- Normal text: 4.5:1 contrast ratio
- Large text/graphics: 3:1 contrast ratio

**Verification:**
```
Text Primary (#e0e0e0) on Background (#1e1e22): 11.7:1 ✓
Text Secondary (#9e9e9e) on Background (#1e1e22): 6.1:1 ✓
Point Color (#6a9bd4) on Background (#1e1e22): 5.8:1 ✓
```

## Appendix B: Reference Implementations

Study these for design inspiration:

1. **Desmos** (desmos.com) - Excellent math visualization UX
2. **GeoGebra** - Dynamic geometry with good controls
3. **3Blue1Brown's Manim** - Beautiful mathematical animations
4. **Observable** - Clean, modern data visualization

## Appendix C: egui Styling Reference

```rust
// Custom egui visuals for clifford-viz
fn configure_visuals(ctx: &egui::Context) {
    let mut visuals = egui::Visuals::dark();

    // Background colors
    visuals.window_fill = Color32::from_rgb(37, 37, 41);
    visuals.panel_fill = Color32::from_rgb(30, 30, 34);

    // Widget colors
    visuals.widgets.noninteractive.bg_fill = Color32::from_rgb(46, 46, 51);
    visuals.widgets.inactive.bg_fill = Color32::from_rgb(58, 58, 64);
    visuals.widgets.hovered.bg_fill = Color32::from_rgb(70, 70, 76);
    visuals.widgets.active.bg_fill = Color32::from_rgb(82, 82, 88);

    // Selection
    visuals.selection.bg_fill = Color32::from_rgb(70, 100, 140);
    visuals.selection.stroke = Stroke::new(1.0, Color32::from_rgb(100, 140, 180));

    ctx.set_visuals(visuals);
}
```

## Open Questions

1. Should we support a light theme, or is dark-only acceptable?
2. Should geometric objects have configurable colors, or enforce the palette?
3. How much animation is too much? (Risk of being distracting)
4. Should we show rulers/measurements along axes?
