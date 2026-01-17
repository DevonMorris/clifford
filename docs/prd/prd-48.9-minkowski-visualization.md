# PRD-48.9: Minkowski Spacetime Visualization

**Status**: Draft
**Parent**: PRD-48
**Depends on**: PRD-48.1
**Goal**: Interactive demos for Minkowski spacetime (special relativity)

## Overview

Visualize spacetime diagrams, light cones, Lorentz transformations, and relativistic effects using the Minkowski algebras (1+1D and 3+1D).

---

## Demo 1: 1+1D Spacetime Diagram

**File**: `examples/visualization/minkowski2_diagram.rs`

### Features

1. **Spacetime Axes**
   - Vertical: time (t)
   - Horizontal: space (x)
   - Units: c = 1 (light travels at 45°)

2. **Light Cones**
   - 45° lines from origin
   - Future light cone (upward)
   - Past light cone (downward)
   - Causal structure visualization

3. **Worldlines**
   - Draw particle trajectories
   - Stationary objects: vertical lines
   - Moving objects: tilted lines
   - Light: 45° lines

4. **Events**
   - Click to place events
   - Label events
   - Connect with worldlines

5. **Lorentz Boost**
   - Slider for velocity v
   - Watch spacetime diagram transform
   - Axes tilt, events move

### Light Cone Structure

```
           t (time)
           ↑
           │      ╱│╲
           │     ╱ │ ╲ future
           │    ╱  │  ╲  light
           │   ╱   │   ╲  cone
           │  ╱    │    ╲
           │ ╱     │     ╲
  ─────────●───────────────→ x (space)
           │╲     │     ╱
           │ ╲    │    ╱
           │  ╲   │   ╱
           │   ╲  │  ╱ past
           │    ╲ │ ╱  light
           │     ╲│╱   cone
           │

Inside light cone: timelike (accessible)
Outside light cone: spacelike (inaccessible)
On light cone: lightlike (light rays)
```

### Implementation

```rust
use clifford::specialized::minkowski::dim2::{Vector, Bivector, Eventor};

pub struct Minkowski2Demo {
    events: Vec<SpacetimeEvent>,
    worldlines: Vec<Worldline>,
    boost_velocity: f32,
    show_light_cones: bool,
    show_grid: bool,
}

struct SpacetimeEvent {
    x: f32,
    t: f32,
    label: String,
}

impl SpacetimeEvent {
    fn to_vector(&self) -> Vector<f32> {
        Vector::new(self.x, self.t)
    }

    fn apply_boost(&self, v: f32) -> SpacetimeEvent {
        // Lorentz boost
        let gamma = 1.0 / (1.0 - v * v).sqrt();
        SpacetimeEvent {
            x: gamma * (self.x - v * self.t),
            t: gamma * (self.t - v * self.x),
            label: self.label.clone(),
        }
    }

    fn proper_time_to(&self, other: &SpacetimeEvent) -> f32 {
        let dt = other.t - self.t;
        let dx = other.x - self.x;
        // τ² = Δt² - Δx² (with c=1)
        let tau_sq = dt * dt - dx * dx;
        if tau_sq > 0.0 {
            tau_sq.sqrt() // Timelike interval
        } else {
            f32::NAN // Spacelike interval (no proper time)
        }
    }
}
```

### Visualization Layout

```
┌─────────────────────────────────────────────────────────────┐
│  1+1D Spacetime Diagram                                [×]  │
├─────────────────────────────────────┬───────────────────────┤
│        t                            │ Controls              │
│        ↑      ╱                     │                       │
│        │     ╱                      │ Boost velocity:       │
│        │    ╱ light                 │ v = [====●====] 0.5c  │
│        │   ╱                        │                       │
│     B  ●  ╱                         │ [Apply Boost]         │
│        │ ╱                          │                       │
│        │╱                           │ Events:               │
│ ───────●────────→ x                │ ● A: (0, 0)           │
│       ╱│                            │ ● B: (0, 2)           │
│      ╱ │                            │ ● C: (1.5, 3)         │
│     ╱  ●─────● C                    │                       │
│    ╱   A                            │ Interval A→C:         │
│   ╱                                 │ Δs² = 6.75 (timelike) │
│                                     │ τ = 2.60              │
│                                     │                       │
│                                     │ [✓] Show light cones  │
│                                     │ [✓] Show grid         │
│                                     │ [Add Event]           │
└─────────────────────────────────────┴───────────────────────┘
```

---

## Demo 2: Time Dilation

**File**: `examples/visualization/minkowski2_dilation.rs`

### Features

1. **Twin Paradox Setup**
   - Stationary twin (vertical worldline)
   - Traveling twin (diagonal, then return)

2. **Clock Comparison**
   - Animate proper time along each worldline
   - Show traveling twin ages less!

3. **Simultaneity**
   - Lines of simultaneity tilt under boost
   - Events simultaneous in one frame not in another

### Twin Paradox Visualization

```
        t
        ↑
        │         stay-at-home
    10  ●─────────● reunion
        │        ╱
        │       ╱
        │      ╱ traveler
        │     ╱  returns
        │    ╱
     5  │   ●  turnaround
        │    ╲
        │     ╲
        │      ╲ traveler
        │       ╲ departs
    0   ●────────────→ x
      start

Stay-at-home: ages 10 years
Traveler: ages 6 years (!)

Proper time is path length in spacetime.
Straight path (stationary) is LONGEST!
```

---

## Demo 3: 3+1D Spacetime

**File**: `examples/visualization/minkowski3_spacetime.rs`

### Features

1. **3D Space + Time**
   - Time as vertical axis or color
   - 3D spatial slices

2. **Light Cone (3D)**
   - Cone surface in 4D projected to 3D
   - Future/past cones

3. **Lorentz Transformations**
   - Boost in any direction
   - Spatial rotations
   - Combined transformations

4. **Worldlines in 3D**
   - Particle trajectories
   - Circular motion (accelerated)

### Implementation

```rust
use clifford::specialized::minkowski::dim3::{Vector, Bivector, Eventor};

pub struct Minkowski3Demo {
    events: Vec<Spacetime4Event>,
    worldlines: Vec<Worldline3D>,
    boost_direction: [f32; 3],
    boost_speed: f32,
    camera: Camera3D,
}

struct Spacetime4Event {
    x: f32,
    y: f32,
    z: f32,
    t: f32,
}

impl Spacetime4Event {
    fn to_vector(&self) -> Vector<f32> {
        Vector::new(self.x, self.y, self.z, self.t)
    }

    fn interval_squared(&self, other: &Self) -> f32 {
        let dt = other.t - self.t;
        let dx = other.x - self.x;
        let dy = other.y - self.y;
        let dz = other.z - self.z;
        dt * dt - dx * dx - dy * dy - dz * dz
    }
}
```

### 3D Light Cone

```
         t
         ↑
         │    ╱│╲
         │   ╱ │ ╲
         │  ╱  │  ╲
         │ ╱   │   ╲
         │╱    │    ╲
    ─────●─────────────── y
        ╱│╲
       ╱ │ ╲
      ╱  │  ╲
     ╱   │   ╲
    x

Light cone is a true cone in 3D!
(4D cone projected to 3D spatial slice)
```

---

## Demo 4: Relativistic Effects Gallery

**File**: `examples/visualization/minkowski_effects.rs`

### Features

1. **Length Contraction**
   - Moving rod appears shorter
   - Show in spacetime diagram
   - Slider for velocity

2. **Time Dilation**
   - Moving clock runs slow
   - Compare tick intervals
   - Animate clock faces

3. **Relativity of Simultaneity**
   - Train/platform thought experiment
   - Lightning strikes at ends
   - Different observers disagree!

4. **Velocity Addition**
   - Two boosts don't add linearly
   - v₁ ⊕ v₂ ≠ v₁ + v₂
   - Show Minkowski sum formula

### Length Contraction

```
Rest frame:                  Moving frame (v = 0.8c):
┌────────────────────┐       ┌────────────┐
│   rod (L₀ = 10m)   │       │ rod (L = 6m)│  L = L₀√(1-v²/c²)
└────────────────────┘       └────────────┘

Spacetime view:
        t                           t'
        ↑                           ↑
        │  ║                       ╱│
        │  ║ rod                  ╱ │
        │  ║ (at rest)           ╱  │ rod
        │  ║                    ╱   │ (moving)
    ────║══════║────→ x     ───╱════╱───→ x'
     end1    end2              │   │
                            appears
                            shorter!
```

### Relativity of Simultaneity

```
Platform frame:                Train frame:
     ⚡ simultaneous ⚡              ⚡ NOT simultaneous ⚡
     │               │              │                 │
 ════●═══════════════●════      ════●═════════════════●════
    rear    TRAIN   front          rear    TRAIN    front
                                   (sees front flash first!)

In spacetime:
        t                           t'
        │   ⚡       ⚡               │  ⚡
        │   │       │              ╱│      ⚡
        │   │       │             ╱ │     ╱
    ────●═══●═══════●───→ x   ───●══●════●───→ x'
        simultaneous             not simultaneous!
```

---

## Implementation Tasks

### minkowski2_diagram.rs
1. [ ] Spacetime axes rendering (t vs x)
2. [ ] Light cone rendering (45° lines)
3. [ ] Event placement
4. [ ] Worldline drawing
5. [ ] Lorentz boost transformation
6. [ ] Interval calculation (timelike/spacelike)
7. [ ] Proper time display

### minkowski2_dilation.rs
1. [ ] Twin paradox worldlines
2. [ ] Proper time integration along path
3. [ ] Clock animation
4. [ ] Simultaneity lines
5. [ ] Age comparison display

### minkowski3_spacetime.rs
1. [ ] 3D+t visualization (time as height or color)
2. [ ] 3D light cone rendering
3. [ ] Boost in arbitrary direction
4. [ ] 3D worldlines
5. [ ] Camera controls

### minkowski_effects.rs
1. [ ] Length contraction demo
2. [ ] Time dilation demo
3. [ ] Simultaneity demo
4. [ ] Velocity addition demo
5. [ ] Interactive velocity slider

## Verification

```bash
cargo run --example minkowski2_diagram --release
cargo run --example minkowski2_dilation --release
cargo run --example minkowski3_spacetime --release
cargo run --example minkowski_effects --release
```

### Checklist
- [ ] Light rays travel at 45° (c = 1)
- [ ] Lorentz boost transforms coordinates correctly
- [ ] Interval Δs² is invariant under boosts
- [ ] Time dilation formula: Δt' = γΔτ
- [ ] Length contraction formula: L = L₀/γ
- [ ] Velocity addition: u ⊕ v = (u+v)/(1+uv/c²)
- [ ] Twin paradox shows correct age difference
- [ ] Simultaneity lines tilt under boost
