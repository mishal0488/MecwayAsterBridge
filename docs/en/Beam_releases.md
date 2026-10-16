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
- Nonlinear, transient, buckling and complex modal releases are not yet
  supported. Neither are combinations with mesh-order upgrades,
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
