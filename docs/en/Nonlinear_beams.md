# Elastic beams with geometric nonlinearity

Select **Elements > Beams (line2) > POU_D_T_GD**, choose **nonlinear
static**, and set **DEFORMATION = GROT_GDEP**. The existing POU_D_E default
is unchanged. Use an elastic material and two-node beam elements.

The bridge generates Code_Aster's existing large-displacement,
large-rotation Timoshenko beam formulation:

```python
MODELISATION='POU_D_T_GD'
RELATION='ELAS_POUTRE_GR'
DEFORMATION='GROT_GDEP'
```

This enables geometric nonlinearity and second-order effects. It does not
add yielding or plasticity. Shear deformation is included, so small-load
results need not exactly match Euler-Bernoulli POU_D_E for short beams.
Section properties and local orientations follow the existing beam
translation. General sections receive the existing AY/AZ coefficients
(estimated for I/T/C/L profiles; 1.0 with a warning for general A/I/J input).

Use incremental loading, preferably an explicit Mecway load history, and
check convergence and mesh refinement. Gravity and line loads use the
existing bridge load preparation. This change does not automatically turn
local loads into follower loads: check whether a load is fixed in space or
must rotate with the structure. Rigid/interpolation couplings retain their
existing constraint equations and need particular care for large motion.

## Results

Displacements and reactions use the existing Mecway result routes. Enable
**MED export** to inspect local beam section resultants `N`, `VY`, `VZ`,
`MT`, `MFY`, `MFZ` in `SIEF_ELNO`. These are forces/moments, not stress
tensor components. The Code_Aster 17.4 large-rotation beam does not implement
`SIPM_ELNO`; extreme beam stresses are therefore not exported to Mecway for
this formulation. Existing beam stress output remains unchanged for
POU_D_E/POU_D_T.

## Initial scope

The bridge rejects this option for linear, modal, buckling and dynamic
analyses, plastic beam materials, non-SEG2 beams, and centrifugal loading.
These restrictions apply to this initial integration, not to all possible
Code_Aster workflows.

Automatic buckling-mode imperfections are not available with this option:
the current pre-solve requires `RIGI_MECA` and `RIGI_GEOM`, which are absent
from this element's 17.4 catalogue. An initially curved input mesh can be
used directly for second-order analysis. Large prescribed rotations and
rigid connections require compatible finite-rotation boundary conditions;
changing the beam formulation alone does not change those constraints.

## Verification

Offline command-generation tests:

```text
python tests/test_bridge_nonlinear_beams.py
```

Numerical acceptance cases are in
`TestFiles/Verificatie/tests_nonlinear_beams.py`, for the existing verification
runner: GB-01 (small-load Timoshenko cantilever), GB-02 (90-degree pure
bending), GB-03 (initially curved column amplification).

Validated on 15 October 2026 with Code_Aster 17.4 on Windows: all three
jobs reported `OK` and reached the full load at t=1, with MED export enabled.

| Case | Quantity | Computed | Reference | Relative error |
| --- | --- | ---: | ---: | ---: |
| GB-01 | Tip DY (m) | 0.00095230554 | 0.00095245524 | -0.0157% |
| GB-02 | Tip DX (m) | -0.36333933 | -0.36338023 | 0.0113% |
| GB-02 | Tip DY (m) | 0.63666068 | 0.63661977 | 0.00643% |
| GB-03 | Total offset / initial offset | 1.9995730 | 2.0 | -0.0213% |

GB-01 uses the existing relative-or-absolute convergence criterion with an
absolute residual of 1e-6 N: excessively strict relative-only tolerances
can stagnate at round-off for small loads. These benchmarks validate the
tested rectangular cantilevers; they do not establish numerical accuracy
for all section shapes, distributed/follower loads, or connection types.
All seven offline test scripts also passed, including ten new beam tests.

References checked against Code_Aster 17.4:

- [Large-rotation beam element catalogue](https://gitlab.com/codeaster/src/-/blob/17.4.0/catalo/cataelem/Elements/meca_pou_d_t_gd.py)
- [SSNL103A large-rotation beam test](https://gitlab.com/codeaster/src/-/blob/17.4.0/astest/ssnl103a.comm)
- [Nonlinear behaviour documentation](https://codeaster.gitlab.io/doc/docaster/manuals/man_u/u4/u4.51.11/Mot_cl__COMPORTEMENT.html)
