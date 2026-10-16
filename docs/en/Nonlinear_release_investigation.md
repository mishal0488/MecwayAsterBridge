# Exact nonlinear beam releases, without offsets

Investigation: Code_Aster 17.4, 16 October 2026. Scope: static elastic
POU_D_T_GD beams in general 3D finite rotations; ideal connections only.

| Released local rotations | Bridge result |
| --- | --- |
| U, V, W (each individually) | No validated exact native implementation found; rejected |
| UV, UW, VW | No validated exact native implementation found; rejected |
| UVW | Exact coincident-translation constraint implemented and benchmarked |

## Findings

The installed `AFFE_CHAR_MECA` catalogue provides `LIAISON_DDL`,
`LIAISON_GROUP`, `LIAISON_UNIF` and related linear constraints. The
[official constraint documentation](https://codeaster.gitlab.io/doc/docaster/manuals/man_u/u4/u4.44.01/G_n_ralit_s.html)
states that their coefficients do not depend on deformation and remain
constant during the calculation. The installed beam-property catalogue
has no end-release keyword for POU_D_T_GD. Discrete spring elements and
their constitutive laws supply finite stiffness; that does not implement
an ideal retained rotational constraint and is outside the requested scope.
No suitable native finite-rotation selective-joint operator was established.
This is a bounded finding about the investigated native operators, not a
proof that a custom Code_Aster solver extension is impossible.

For selective releases, taking components of the difference of two global
rotation vectors is not equivalent to constraining relative orientation.
For example, a 0.6 rad permitted local rotation following a common 1.1 rad
rotation produces a spurious retained component of about 0.32958 rad in
the old linear equations. Equivalent counterexamples exist for all six
selective combinations. Purely planar examples or infinitesimal rotations
can conceal this error. Passing such examples does not validate a general
3D joint implementation. The old equations therefore remain disabled for
selective geometrically nonlinear releases.

For UVW, only x_end = x_junction is required. Initially coincident nodes
can enforce this with three homogeneous translation equations at any
rotation. No rotational constraints are necessary. This gives an ideal
spherical joint within the beam formulation and solver tolerances.

The separate CHARGES_36 rigid-offset probe is not evidence that selective
end releases are impossible. That result applies to follower
LIAISON_SOLIDE on nodes with rotational degrees of freedom. Offsets are
excluded from this implementation and from the spherical-joint benchmarks.

## Verification

`TestFiles/Verificatie/check_nonlinear_releases.py` builds two independent
coincident 1 m cantilevers with a shared tip and a UVW release on one tip.
The 40-element-per-beam bending model receives 700 N total transverse force.
It gives tip displacement (-0.16060934, 0.49345016, 0) m, versus the
independent inextensible elastica (-0.16064172, 0.49345748, 0) m. The maximum
component difference is about 0.033 mm; the reference neglects the small
shear and extension terms present in the FE model. An obliquely rotated
version gives the same displacement after transformation back. A separate
pure-torsion case gives pi/2 rotation on the loaded side and zero twist on
the released side. Raw nodal moments at the released end vanish; nodal
translations coincide and support reactions balance applied forces.

Translation regression tests cover rejection of every selective combination,
acceptance of UVW, command generation and rejection of nonlinear offsets.
The existing linear and modal release cases remain supported.

## Remaining work for selective releases

A dedicated nonlinear constraint integration must define the relative
orientation convention for each joint, consistent residuals and Jacobians,
constraint tangent contributions, and reaction recovery with the solver's
finite-rotation updates. It must then pass moving-frame and compound-rotation
tests, as well as force and moment equilibrium and independent frame
benchmarks. This work is not implemented by the UVW addition.
