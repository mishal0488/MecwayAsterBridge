# -*- coding: utf-8 -*-
"""Elastic geometric-nonlinear SEG2 beam benchmarks.

Run through run_verificatie.py. Validated with Code_Aster 17.4 on Windows
(2026-10-15; see docs/en/Nonlinear_beams.md). No plasticity or automatic buckling
imperfection is used. Length and section refinement should be checked when
extending this implementation to new section/load combinations.
"""
import math
from verif_lib import Case, Check, Model, linemesh

GROUP = "Geometric nonlinear beams"
E, NU = 210e9, 0.3
L, HY, HZ = 1.0, 0.01, 0.02
IZ = HZ * HY ** 3 / 12.0
AREA = HY * HZ
F_SMALL = 1.0
GENERAL = {"elem_beam": "POU_D_T_GD", "med_export": True}
SETTINGS = {"deformation": "GROT_GDEP", "substeps_per_interval": 40,
            "newton_max_iter": 40, "convergence_resi_rela": 1e-7}


def model(imperfection=0.0):
    m = Model("L30")
    m.elset("Beam", "Steel")
    _, nodes = linemesh(m, "Beam", 40, (0, 0, 0), (L, 0, 0))
    if imperfection:
        # First cantilever buckling mode, zero displacement/slope at root.
        m.nodes = [(x, imperfection * (1 - math.cos(math.pi * x / (2 * L))), z)
                   for x, y, z in m.nodes]
    m.mat("Steel", '<geometric type="RectBar" recta="%g m" rectb="%g m" />'
          % (HY, HZ), '<mechanical type="Isotropic" youngsmodulus="%g Pa" '
          'poissonratio="%g" />' % (E, NU))
    m.node_selection("Root", [nodes[0]])
    m.node_selection("Tip", [nodes[-1]])
    m.add('<fixed selection="Root" />')
    return m


def moment_model(theta):
    m = model()
    moment = theta * E * IZ / L
    m.add('<moment selection="Tip"><x>0 N.m</x><y>0 N.m</y>'
          '<z time="0 s">0 N.m</z><z time="1 s">%.12g N.m</z></moment>'
          % moment)
    return m


def small_force():
    m = model()
    m.add('<force selection="Tip"><x>0 N</x><y time="0 s">0 N</y>'
          '<y time="1 s">%g N</y><z>0 N</z></force>' % F_SMALL)
    return m


E0 = 0.0001 * L
ALPHA = 0.5
PCR = math.pi ** 2 * E * IZ / (4 * L ** 2)


def imperfect_column():
    m = model(E0)
    m.add('<force selection="Tip"><x time="0 s">0 N</x>'
          '<x time="1 s">%.12g N</x><y>0 N</y><z>0 N</z></force>'
          % (-ALPHA * PCR))
    return m


def tip(result, field):
    return result.node(field, 41)


THETA = math.pi / 2
TESTS = [
    Case("GB-01", GROUP, "Small-load Timoshenko cantilever", "nonlin_static_3d",
         small_force,
         [Check("tip DY", lambda r: tip(r, "disply"),
                ref=F_SMALL * (L ** 3 / (3 * E * IZ) +
                            1.2 * L / ((E / (2 * (1 + NU))) * AREA)),
                tol=0.01, unit="m", source="FL^3/(3EI) + F L/(kGA), k=5/6")],
         settings=dict(SETTINGS, convergence_criterion="RELA_OF_MAXI",
                       convergence_resi_maxi=1e-6), general=GENERAL,
         note="Absolute residual 1e-6 N handles round-off at small loads."),
    Case("GB-02", GROUP, "Cantilever bent through 90 degrees by an end moment",
         "nonlin_static_3d", lambda: moment_model(THETA),
         [Check("tip DX", lambda r: tip(r, "displx"),
                ref=L * (math.sin(THETA) / THETA - 1), tol=0.01, unit="m",
                source="Circular arc under pure bending: x=L sin(theta)/theta"),
          Check("tip DY", lambda r: tip(r, "disply"),
                ref=L * (1 - math.cos(THETA)) / THETA, tol=0.01, unit="m",
                source="Circular arc under pure bending: y=L(1-cos(theta))/theta")],
         settings=SETTINGS, general=GENERAL),
    Case("GB-03", GROUP, "Initially curved column: second-order amplification",
         "nonlin_static_3d", imperfect_column,
         [Check("total tip offset / e0", lambda r: 1 + tip(r, "disply") / E0,
                ref=1 / (1 - ALPHA), tol=0.02, unit="-",
                source="Elastic small-slope amplification 1/(1-P/Pcr)")],
         settings=SETTINGS, general=GENERAL,
         note="Initial curvature is in the input mesh; no eigenmode pre-solve."),
]
