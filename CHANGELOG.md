# Changelog

## V2026.1 – first public release

First public version of the MecwayAsterBridge: a Mecway 33 model (`.liml`) is solved with
Code-Aster 17.4 for Windows and the results come back into Mecway. Standalone: only the
Python of Mecway and of Code-Aster, no `pip` packages.

- **Start screen** in English and Dutch, with settings per analysis, tooltips, a
  convergence and messages window, and the documentation built in.
- **Analyses**: linear and nonlinear static (contact, plasticity, large displacements,
  element failure), implicit and explicit nonlinear dynamics, mode superposition and
  harmonic response, thermal (steady state and transient) and thermo-mechanical, modal
  (also prestressed), linear and nonlinear buckling (CRIT_STAB, imperfection).
- **Acoustics**: cavity modes and harmonic pressure field (sources, absorbing boundaries,
  microphones); vibroacoustics, coupled modes and harmonic response.
- **Welding and metallurgy**: moving heat source along a weld bead, residual stress and
  distortion, filler material that is born, phases and hardness.
- **Fracture mechanics**: K and G along a meshed crack, and XFEM cracks without a mesh.
- **Fatigue**: damage and life per node (Wöhler curve, rainflow, Miner, critical-plane
  criteria, Manson-Coffin).
- **Cyclic symmetry** for modal analysis of a sector.
- **Response spectrum** (earthquake, shock) on a frequency model.
- **Error estimator** per element (linear static) and **adaptive refinement** of the
  mesh in Mecway.
- **`#ASTER` lines** in Mecway's *Custom model definition* for Code-Aster input that
  Mecway has no field for.
- **Verification set**: 272 test cases with 1629 checks against analytical solutions
  and benchmarks, plus 19 official Code-Aster validation tests
  ([`docs/en/Verification.md`](docs/en/Verification.md)).
