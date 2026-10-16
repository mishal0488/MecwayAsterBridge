# Beam offsets: linear static and natural frequencies

The bridge reads Mecway beam **Offset** definitions and applies them to
**linear static 3D** and ordinary **modal 3D** analyses. Select `POU_D_E`
(Euler-Bernoulli) or `POU_D_T` (Timoshenko), two-node beams, and the **LIML**
results route. For natural frequencies, turn **modal prestress off** and
use real modes (not complex modes).

Create the offset while editing a static Mecway model, then retain its
`<offset>` definition and element selection in the modal model. Mecway 34
was checked to preserve those definitions when opening/saving a model
whose analysis type was changed in its LIML file from S30 to V30 and back.
This does not add a modal-offset editor to Mecway. Inspect the bridge log
for `Beam offsets:` to confirm that offsets were applied.

## Meaning and implementation

The offset components U, V and W are constant lengths along each beam's
local axes, including its orientation. Global X, Y and Z offsets are also
read directly in the global coordinate system. Do not mix the two systems
in one offset definition. They move the physical beam
centreline away from the original reference nodes. Original nodes still
carry nodal loads, supports and connections. The bridge creates internal
physical nodes and connects them by linear rigid-arm equations:

```
u_physical = u_reference + rotation_reference cross offset
rotation_physical = rotation_reference
```

Both stiffness and mass are evaluated on the physical centreline. The
constraints therefore include eccentric inertia in natural-frequency
calculations. They do not introduce artificial stiffness or mass.
Separate beams may have different offsets at a common reference node.

Results use the original Mecway model and node numbers. Displacements
and rotations belong to the reference nodes; beam forces, shear, torsion
and bending moments belong to the physical beam section in its local
axes. Support reactions are transferred to reference nodes, including
the moment-arm contribution. Internal generated nodes are excluded from
Mecway results. The solver mesh can be inspected separately to see the
physical centreline. Mecway's own modal drawing rules still apply; this
change does not implement a new offset renderer in Mecway.

Rotational beam-end releases can be combined with offsets. The release is
applied at the original reference junction; the rigid arm connects the
released side to the offset physical beam. Both intermediate and physical
nodes are mapped back to the original node for result import.

Unsupported offset settings display a dialog explaining why Code_Aster
was not started, as well as recording the reason in the bridge log.

## Current limits

- Linear static and ordinary real modal analyses only. Nonlinear,
  prestressed modal, buckling and response analyses with offsets stop
  with an explanatory error. The separate geometric-nonlinear beam
  feature remains available for beams **without offsets**.
- Use node selections for loads and supports on offset beams. Gravity is
  supported. Distributed/element-selected beam loads are rejected.
- Constant, finite offsets on line2 POU_D_E/POU_D_T elements only. No
  overlapping offset definitions on an element.
- No contact, bonded/node-surface coupling, cyclic symmetry, XFEM,
  adaptive refinement or automatic mesh-order upgrades in offset models.
  This includes COQUE_3D shell conversion; Mecway's supplied quadratic
  CompositeBeam sample therefore currently stops at this check.
- Models without offset definitions follow their existing translation.

## Numerical verification

Run from a Python environment with NumPy and the installed Code_Aster:

```
python TestFiles/Verificatie/check_beam_offsets.py PATH_TO_RUN_ASTER NEW_OUTPUT_DIRECTORY
```

The output directory must not already exist. The script runs seven jobs:
static eccentric loading for both beam formulations, offset/zero-offset
cantilever modal comparisons for both, and a constrained Euler-Bernoulli
model compared with independently assembled stiffness and consistent mass
matrices. It checks displacement, force/moment equilibrium and frequencies.

With Code_Aster 17.4, a 1 m steel cantilever, 20 by 30 mm section, 100 mm
local-W offset and 100 N axial force at the reference tip gives
`ux = 0.000106613756614 m`, agreeing with
`F L / EA + F e^2 L / EI`. The root axial reaction is -100 N and the net
reference-node reaction moment is zero. The six constrained-model
frequencies agree with the independent calculation to better than 0.001%.
These checks validate the tested formulations and boundary conditions;
convergence and model-specific checks remain necessary for new models.

The offset/release combination was also checked with two 1 m rectangular
beam spans, fixed outer ends, a V/W release at their junction, a 180 mm
global-Z offset and a 100 N global-Y junction force. Code_Aster 17.4 and
Mecway 34 agreed on all original-node displacements, rotations and support
reactions (largest displacement difference below 1e-15 m). The combined
ordinary-modal case ran with eight positive modes; this is a run check,
not a native Mecway modal-offset equivalence check.
