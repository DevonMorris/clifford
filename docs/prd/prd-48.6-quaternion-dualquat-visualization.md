# PRD-48.6: Quaternion and Dual Quaternion Visualization

**Status**: Draft
**Parent**: PRD-48
**Depends on**: PRD-48.1
**Goal**: Interactive demos for quaternion rotations and dual quaternion rigid transforms

## Overview

Demonstrate quaternion rotation interpolation and dual quaternion screw motion. Focus on practical benefits over matrices/Euler angles.

---

## Demo 1: Quaternion SLERP

**File**: `examples/visualization/quaternion_slerp.rs`

### Features

1. **Start/End Rotation Setup**
   - Define two orientations
   - Show as coordinate frames or objects

2. **Interpolation Comparison**
   - SLERP (spherical linear interpolation) - correct
   - LERP + normalize - incorrect (variable speed)
   - Euler angle interpolation - gimbal lock issues

3. **Interpolation Parameter**
   - Slider for t ∈ [0, 1]
   - Animation mode

4. **Path Visualization**
   - Show interpolation path on unit sphere
   - SLERP = great circle arc
   - LERP = chord through sphere

### SLERP vs LERP

```
Unit quaternion sphere (S³ projected to S²):

SLERP (correct):                LERP (wrong):
        ╭───╮                         ╭───╮
       ╱     ╲                       ╱     ╲
      │   ╭───╲                     │   ╭┈┈┈╲
      │  ╱  arc ╲                   │  ╱chord ╲
   q₀ ●─╯        ● q₁            q₀ ●┈┈┈┈┈┈┈┈● q₁
      │           │                 │           │
       ╲         ╱                   ╲         ╱
        ╰───────╯                     ╰───────╯

SLERP: constant angular velocity
LERP: variable speed (slows in middle)
```

### Implementation

```rust
use clifford::specialized::quaternion::Quaternion;

pub struct QuaternionSlerpDemo {
    q_start: Quaternion<f32>,
    q_end: Quaternion<f32>,
    t: f32,
    method: InterpolationMethod,
    animation: Animation,
}

#[derive(PartialEq)]
enum InterpolationMethod {
    Slerp,      // Spherical linear interpolation
    NLerp,      // Normalized linear interpolation
    Euler,      // Euler angle interpolation
}

impl QuaternionSlerpDemo {
    fn interpolate(&self) -> Quaternion<f32> {
        match self.method {
            InterpolationMethod::Slerp => {
                self.q_start.slerp(&self.q_end, self.t)
            }
            InterpolationMethod::NLerp => {
                let q = self.q_start.lerp(&self.q_end, self.t);
                q.normalize()
            }
            InterpolationMethod::Euler => {
                // Extract Euler angles, interpolate, reconstruct
                // (shows gimbal lock problems)
                todo!()
            }
        }
    }
}
```

### Visualization Layout

```
┌─────────────────────────────────────────────────────────────┐
│  Quaternion SLERP                                      [×]  │
├─────────────────────────────────────┬───────────────────────┤
│                                     │ Setup                 │
│         ╭───────╮                   │                       │
│        ╱   ↗    ╲  ← q(t)          │ Start rotation:       │
│       │  ╱       │                  │   Axis: [x][y][z]     │
│   q₀ →●╱─────────● ← q₁            │   Angle: [====●====]  │
│       │╲         │                  │                       │
│        ╲ ╲      ╱                   │ End rotation:         │
│         ╰─╲────╯                    │   Axis: [x][y][z]     │
│            ╲ path                   │   Angle: [====●====]  │
│                                     │                       │
│   ┌─────┐      ┌─────┐             │ Interpolation         │
│   │start│      │ end │             │ t: [====●====] 0.50   │
│   └─────┘      └─────┘             │                       │
│       ↓                            │ Method:               │
│   ┌─────┐ interpolated             │ ● SLERP (correct)     │
│   │     │                          │ ○ NLERP (approx)      │
│   └─────┘                          │ ○ Euler (broken)      │
│                                     │                       │
│                                     │ [▶ Animate]           │
└─────────────────────────────────────┴───────────────────────┘
```

---

## Demo 2: Quaternion Double Cover

**File**: `examples/visualization/quaternion_double.rs`

### Features

1. **360° vs 720° Rotation**
   - Animate object rotating
   - At 360°: object back to start, but q flipped to -q
   - At 720°: both object AND quaternion back to start

2. **Belt Trick Visualization**
   - Classic "plate on palm" demonstration
   - 360° rotation tangles imaginary belt
   - 720° rotation untangles it

3. **Quaternion Path**
   - Show path in quaternion space
   - 360° = path from q to -q
   - 720° = closed loop back to q

### The Double Cover Phenomenon

```
Physical rotation:     Quaternion:

0° → 360°:            q → -q (antipodal!)
  Object: back to     Different quaternions,
  original pose       same rotation!

360° → 720°:          -q → q (back to start)
  Object: still       Now quaternion also
  original pose       back to start!

Quaternions are a DOUBLE COVER of rotations.
q and -q represent the SAME rotation.
```

---

## Demo 3: Dual Quaternion Screw Motion

**File**: `examples/visualization/dualquat_screw.rs`

### Features

1. **Screw Axis Definition**
   - Line in 3D space (point + direction)
   - Pitch: translation per rotation

2. **Screw Motion Animation**
   - Object follows helical path
   - Rotation and translation coupled

3. **Parameter Control**
   - Rotation angle
   - Pitch (translation/rotation ratio)
   - Screw axis position and direction

4. **Special Cases**
   - Pitch = 0: pure rotation
   - Pitch = ∞: pure translation
   - Pitch > 0: right-handed screw
   - Pitch < 0: left-handed screw

### Screw Motion Visualization

```
Pitch = 0 (pure rotation):     Pitch > 0 (screw):
        │                             │
    ┌───┼───┐                     ┌───┼───┐
    │   │   │  rotates            │   │   │
    │   ●───│──around             │   ●   │  spiral
    │   │   │  axis               │   │↘  │  path
    └───┼───┘                     └───┼─↘─┘
        │                             │  ↘
                                      ↓

Pitch = ∞ (pure translation):
    ┌───────┐     ┌───────┐
    │       │ ──► │       │
    │   ●   │     │   ●   │
    │       │     │       │
    └───────┘     └───────┘
```

### Implementation

```rust
use clifford::specialized::dualquat::DualQuaternion;

pub struct DualQuatScrewDemo {
    // Screw parameters
    axis_point: [f32; 3],
    axis_direction: [f32; 3],
    angle: f32,
    pitch: f32,

    animation: Animation,
}

impl DualQuatScrewDemo {
    fn compute_dual_quaternion(&self) -> DualQuaternion<f32> {
        // Build dual quaternion from screw parameters
        DualQuaternion::from_screw(
            self.axis_point,
            self.axis_direction,
            self.angle,
            self.pitch,
        )
    }
}
```

### Visualization Layout

```
┌─────────────────────────────────────────────────────────────┐
│  Dual Quaternion Screw Motion                          [×]  │
├─────────────────────────────────────┬───────────────────────┤
│                                     │ Screw Parameters      │
│                                     │                       │
│           axis                      │ Axis Point:           │
│            ║                        │   [0.0][0.0][0.0]     │
│            ║                        │                       │
│   ┌────────╫────────┐               │ Axis Direction:       │
│   │        ║        │               │   [0.0][0.0][1.0]     │
│   │    ┌───╫───┐    │               │                       │
│   │    │   ║   │    │               │ Rotation: [===●===]   │
│   │    │   ●   │    │               │           90.0°       │
│   │    │  ╱║   │    │               │                       │
│   │    └─╱─╫───┘    │               │ Pitch: [====●====]    │
│   │     ╱  ║        │               │        0.5 units/rad  │
│   │    ╱   ║        │               │                       │
│   │   ↙    ║        │ trail         │ ─────────────────     │
│   └────────╫────────┘               │ Motion Type:          │
│            ║                        │   ● Screw (general)   │
│            ║                        │   ○ Pure rotation     │
│                                     │   ○ Pure translation  │
│                                     │                       │
│                                     │ [▶ Animate]           │
└─────────────────────────────────────┴───────────────────────┘
```

---

## Demo 4: Dual Quaternion Blending

**File**: `examples/visualization/dualquat_blend.rs`

### Features

1. **Multiple Poses**
   - Define 2-4 poses as dual quaternions
   - Each pose has a weight

2. **DLB (Dual Quaternion Linear Blending)**
   - Blend poses with weights
   - Compare with matrix blending (shows artifacts)

3. **Character Skinning Preview**
   - Simple mesh deformed by blended transforms
   - Show advantage over linear blend skinning

### ScLERP (Screw Linear Interpolation)

```
Dual quaternion interpolation preserves:
• Rigidity (no shearing/scaling)
• Volume (no collapse)
• Smooth motion (no discontinuities)

Matrix interpolation problems:
• Can shrink ("candy wrapper" effect)
• Can shear
• Discontinuous at boundaries
```

---

## Implementation Tasks

### quaternion_slerp.rs
1. [ ] Start/end quaternion definition UI
2. [ ] SLERP implementation
3. [ ] NLERP implementation
4. [ ] Euler interpolation (for comparison)
5. [ ] Interpolation path visualization
6. [ ] Rotating object display
7. [ ] Animation support

### quaternion_double.rs
1. [ ] Continuous rotation animation
2. [ ] Quaternion value display
3. [ ] 360° vs 720° visualization
4. [ ] Quaternion space path drawing
5. [ ] Belt trick visualization (optional)

### dualquat_screw.rs
1. [ ] Screw axis definition UI
2. [ ] Dual quaternion from screw params
3. [ ] Screw motion animation
4. [ ] Motion trail rendering
5. [ ] Pitch variation (pure rotation to pure translation)
6. [ ] Object transformation display

### dualquat_blend.rs
1. [ ] Multiple pose definition
2. [ ] Weight sliders
3. [ ] DLB implementation
4. [ ] Matrix blend comparison
5. [ ] Simple mesh deformation (optional)

## Verification

```bash
cargo run --example quaternion_slerp --release
cargo run --example quaternion_double --release
cargo run --example dualquat_screw --release
cargo run --example dualquat_blend --release
```

### Checklist
- [ ] SLERP produces constant angular velocity
- [ ] NLERP shows variable speed (compare to SLERP)
- [ ] Euler interpolation shows gimbal lock
- [ ] Double cover: q at 0°, -q at 360°, q at 720°
- [ ] Screw motion follows correct helical path
- [ ] Pitch=0 gives pure rotation
- [ ] DLB produces rigid transforms
