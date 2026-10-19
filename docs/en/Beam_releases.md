# Flexible joints / beam-end releases

Mecway's **Flexible joint on beam** is translated for linear static and
real modal 3D analyses, including modal prestress when no offsets are present. Use two-node `POU_D_E`
or `POU_D_T` beams and the LIML results route. Apply the joint to beam end
faces, as in Mecway; no additional bridge setting is needed.

The bridge reads `<rotu/>`, `<rotv/>` and `<rotw/>` from the `flexiblejoint`
card. The selected rotations are free in the beam's local axes, including
selection-based orientation cards. All three translations remain connected.
Unreleased rotations remain connected too. Releasing V and W therefore
transmits torque but no bending moments; releasing U removes torque transfer.
Suppressed joints are ignored.

Each released beam end receives an independent, coincident solver node.
Exact linear equations connect translations and retained local rotations.
No penalty stiffness, artificial springs or artificial mass are added.
Internal end nodes are mapped back to the original elements for section
forces and moments. Nodal results refer to the original junction nodes;
reactions include the contributions from the connected beam ends. The
log reports `Beam end releases: ... independent beam ends`.

Offsets and releases can be combined for linear static and ordinary real
modal analyses. Turn modal prestress off when offsets are present. The
release acts at the reference junction, before the rigid offset arm.

## Current limits

- Each released junction must also connect to at least one unreleased beam
  end. Isolated released ends and junctions where every incident beam end
  is released stop with an explanatory error; their unused reference
  rotations require additional elimination.
- Use node selections for loads/supports at released ends. Face-selected
  loads/supports directly on those ends are rejected.
- Only free rotational releases are supported, not finite joint stiffness
  or translational releases.
- Selective nonlinear releases U, V, W, UV and UW, transient, buckling
  and complex modal releases are not yet supported. Neither are combinations with mesh-order upgrades,
  contact/bonded/node-surface couplings, cyclic symmetry, adaptive
  refinement or XFEM. These combinations stop before solving.

## Verification

The numerical check uses two 1 m beam spans, fixed at the outer ends and
joined at the centre. Releasing the appropriate bending rotation at one
side of the junction gives a centre displacement of `F L^3 / (6 E I)`;
the released end moment is zero. A U-only release transmits no torque to
the released span. The modal check compares the first symmetric bending
mode with the analytical cantilever frequency.

```
python TestFiles/Verificatie/check_beam_releases.py PATH_TO_RUN_ASTER NEW_OUTPUT_DIRECTORY
```

Seven Code_Aster checks cover U, V, W, V+W and U+V+W releases, static
equilibrium, section resultants, real modes and Timoshenko shear compliance.
Translation tests additionally cover rotated local axes, suppression and
rejection of unsupported cases. A separate 30-mode frame comparison with
Mecway 34 verified both the frequencies and displacement/rotation mode
shapes; the largest frequency difference was below 0.00002%.

## Geometric nonlinear static: exact spherical joints

With `POU_D_T_GD`, `GROT_GDEP`, elastic materials and no offsets, full
**U+V+W** releases are supported. Each released end has a coincident
independent node. Only its three translations are tied to the junction;
all rotations are independent. These coincidence equations are exact at
finite rotations. No stiff springs, short flexible beams or penalty
stiffness are introduced. Existing restrictions on unsupported mesh/load
features and isolated/all-released junctions still apply.

This does **not** enable U, V, W, UV or UW releases in general 3D
large-rotation analysis. They are rejected with a specific message, without
changing them to spherical joints. See [the investigation](Nonlinear_release_investigation.md).

Run `python TestFiles/Verificatie/check_nonlinear_releases.py RUN_ASTER NEW_OUTPUT_DIRECTORY`.
Three Code_Aster 17.4 cases verify large-deflection bending against an
independent elastica solution, invariance under a 3D rotation of the model,
and a 90-degree relative twist with no transmitted joint torque. The
benchmark checks coincident end displacement and all three nodal joint
moments directly in the solver result. It does not infer hinge equilibrium
from the extrapolated `SIEF_ELNO` section moment at a finite-length element.

## Geometric nonlinear static: custom exact V+W joints

Combined **V+W** releases (U retained) now use a custom nonlinear constraint
solve around Code_Aster `STAT_NON_LINE`. Select `POU_D_T_GD`, `GROT_GDEP`,
elastic materials and the LIML route. Apply the normal Mecway flexible-joint
card with V and W checked. **Offsets are excluded.** No extra script setting
or spring stiffness is required.

At finite rotations, a joint needs an explicit relative-orientation convention.
This implementation uses **zero relative swing-twist torsion**: two bending
(swing) freedoms with no relative twist about the beam's initial U axis.
For unit rotation quaternions `qa`, `qb` and initial unit axis `u`, it solves
`C = 2 u · vector(conjugate(qa) qb) = 0`. This is not a claim that every
physical universal/Cardan joint has the same finite-rotation kinematics.
The chart excludes a relative swing of 180 degrees.

Translations coincide exactly. The rotational constraint uses native Lagrange
multipliers and the spatial derivative of C, giving equal and opposite joint
couples. The direction is updated in an outer iteration; it is not a fixed
global rotation-component tie. At each accepted load step, the runtime checks
the quaternion constraint and the released-end moment against that derivative.
At large relative swing, the allowed reaction direction follows this joint
convention; do not interpret a nonzero extrapolated local section moment as
proof that the release is locked.

No penalty springs, short connection elements, artificial mass or finite
connection stiffness are used. Only converged joint states enter the result
history. A native or outer convergence failure stops the run instead of
exporting trial states. The additional nonlinear iterations take extra time.
The inner equilibrium check uses relative tolerance 1e-8 or absolute 1e-5;
joint rotation-update and constraint tolerances are 1e-8 and 1e-9 respectively.
"Exact" refers to the connection formulation, within these numerical tolerances.

The initial integration supports beam-only elastic static models, with the
restrictions above. Stability extraction, automatic imperfections, thermal
histories and custom analysis code are excluded. U+V+W releases can coexist
with V+W releases. Single V or single W releases remain unsupported here.
Retained junction rotations must still belong to an unreleased incident end.

Refine both the beam mesh and load increments for your model. Native internal
substeps are trial states; the custom constraint is converged at the requested
output increments. The combined benchmark's largest change in the joint
rotation components fell from 0.000860 rad (40 to 80 increments) to 0.000222 rad
(80 to 160 increments). These are validation observations, not a universal
accuracy bound or a substitute for checking a particular structure.

### Checks performed with Code_Aster 17.4

- Analytical torsional sharing: two cantilevers both rotate pi/4, with equal
  torque sharing; a second test retains pi/4 with V/W releases at both ends
  of an intermediate beam.
- Combined bending and torsion, including a beam released at both ends:
  converged quaternion constraints and released-end moment equilibrium.
- Rigid rotation of the whole combined-load model: transformed displacement
  and rotation differences below 1e-10; reaction differences below 2e-8.
- Small-load comparison with the existing linear V/W formulation: the
  transverse-displacement difference decreases from approximately 0.60%
  with six elements per span to 0.038% with 24 elements per span.
- Load-step refinement, complete accepted-step result export, and Mecway
  IronPython command generation, mesh-fit checking and LIML result merging.
- The supplied 360-beam frame, tested as a separate nonlinear copy with its
  existing gravity load and V/W releases, completed with no solver errors.
  This is a compatibility check, not independent certification of that frame.

Reproduce the bridge/solver/export checks with:

```
python TestFiles/Verificatie/check_nonlinear_vw.py RUN_ASTER NEW_OUTPUT_DIRECTORY
```

Pass individual case names to run only selected benchmarks, for example
`both_ends_compound`. Quaternion objectivity and derivative checks also run in
the standard `tests/test_aster_joint_kinematics.py` suite.
