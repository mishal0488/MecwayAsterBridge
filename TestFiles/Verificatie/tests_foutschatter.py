# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Foutschatter (lineair statisch, CALC_ERREUR ERME; bridge-instelling
error_estimator, docs/Foutschatter.md). Toelichting: docs/Verificatie.md
hoofdstuk Foutschatter (FS).

FS-01  Ingeklemde balk hex8 (het model van LS-02-h8, n = 1..16): de
       globale schatting daalt met orde ~1 (lineaire elementen, energienorm
       O(h)), de effectiviteitsindex (geschatte / werkelijke relatieve fout
       in de energienorm) is O(1) en constant, de grootste fout zit bij de
       inklemming; plus de doorbuiging tegen de balktheorie (Roark).
FS-02  Dikwandige cilinder hex20 (het model van LS-04, Lamé): orde ~2
       (kwadratische elementen), grootste fout aan de binnenwand.
FS-03  Code-Aster-testcase sslv153a (primair probleem) via de bridge
       nagebouwd: zelfde mesh (hexa/penta/pyram/tetra, LINE_QUAD), lasten
       en vrije randen; ERREST van ERME_NOEU en DZ in knoop MESURE gelijk
       aan de TEST_RESU-waarden (NON_REGRESSION).

Werkelijke fout FS-01: met een gelijkmatige traction t = P/A op het
eindvlak is de uitwendige arbeid P * w (w = gemiddelde doorbuiging van het
eindvlak), dus ||u||_E^2 = P w en (Galerkin-orthogonaliteit)
||u - u_h||_E^2 = P (w_inf - w_h): relatieve fout sqrt((w_inf - w_h) /
w_inf). w_inf = limiet van het 3D-model (Richardson op LS-02-h20 n = 1, 2,
4: 0,0019057 m; LS-02-h8 n = 4, 8, 16: 0,0019067 m), hier 0,001906 m.
"""
import glob
import math
import os
import sys

import verif_lib as V
from verif_lib import Case, Check, Model
import tests_lineair_statisch as LS

if V.SCRIPTS not in sys.path:
    sys.path.insert(0, V.SCRIPTS)

GROUP = "Foutschatter"
W_INF = 0.001906        # m, 3D-limiet van de LS-02-balk (zie de kop)


def _est(r):
    return V.error_summary(r)


def _true_rel(r):
    """Werkelijke relatieve fout in de energienorm (FS-01, zie de kop)."""
    w = LS._mean(r, "displz", "TipKnopen")
    return math.sqrt(max(0.0, (W_INF - (-w)) / W_INF))


def _prev_est(r, case_id):
    """Globale schatting van een andere case uit dezelfde werkmap (de
    grovere mesh van de reeks; die draait eerder)."""
    import bridge_results_liml
    root = os.path.dirname(os.path.dirname(r.mess))
    mess = glob.glob(os.path.join(root, case_id, "*.mess"))
    found = bridge_results_liml.read_error_summary(mess[0]) if mess else []
    if not found:
        raise ValueError("geen foutschatting van %s (eerst die case draaien)"
                         % case_id)
    return found[-1]


def _order(case_prev, ratio=2.0):
    """Orde p van de absolute globale schatting eta tussen de vorige mesh en
    deze (h / 2): p = ln(eta_vorig / eta) / ln 2."""
    def value(r):
        return math.log(_prev_est(r, case_prev)["abs"] / _est(r)["abs"]) / \
            math.log(ratio)
    return value


# ------------------------------------------------------------------
# FS-01 Ingeklemde balk hex8, mesh-reeks
# ------------------------------------------------------------------
def _fs01_cases():
    out = []
    ns = (1, 2, 4, 8, 16)
    for k, n in enumerate(ns):
        h = LS.H2 / n
        checks = [
            Check("doorbuiging tip (gem.)",
                  lambda r: -LS._mean(r, "displz", "TipKnopen"),
                  ref=LS.W2, tol=LS._tol(0.02, h, LS.H2 / 16, 2), unit="m",
                  source="Timoshenko-balk PL^3/3EI + PL/kGA (Roark tabel 8.1 "
                         "geval 1a + afschuiving), zoals LS-02",
                  in_series=False),
            Check("effectiviteitsindex",
                  lambda r: (_est(r)["rel"] / 100.0) / _true_rel(r),
                  ref=2.0, tol=0.5, unit="",
                  source="residuschatter: index O(1) en ~constant "
                         "(Ainsworth & Oden 2000, hfst. 2); werkelijke fout "
                         "sqrt((w_inf - w) / w_inf)", in_series=False),
            Check("grootste fout bij inklemming (x)",
                  lambda r: r.argmax("errorestimate")[0], ref=0.0,
                  tol_abs=LS.L2 / (10 * n) + 1e-9, unit="m",
                  source="singuliere randen van de inklemming",
                  in_series=False),
        ]
        if n >= 4:
            checks.append(Check(
                "orde globale schatting", _order("FS-01-h8-n%d" % ns[k - 1]),
                ref=1.0, tol=0.2, unit="",
                source="lineaire elementen: energienorm O(h) (Verfuerth "
                       "1996; R4.10.02)", in_series=False))
        out.append(Case(
            "FS-01-h8-n%d" % n, GROUP,
            "Foutschatter: ingeklemde balk hex8 %dx%dx%d" % (10 * n, n, n),
            "lin_static_3d", LS.build_ls02(n, 1), checks,
            general={"elem_hex8": "3D"}, h=h,
            note="Model van LS-02-h8; globale schatting uit de .mess "
                 "(CALC_ERREUR)."))
    return out


TESTS = _fs01_cases()


# ------------------------------------------------------------------
# FS-02 Dikwandige cilinder hex20 (Lamé), mesh-reeks
# ------------------------------------------------------------------
def _fs02_cases():
    out = []
    ns = (1, 2, 4, 8)
    for k, n in enumerate(ns):
        h = (LS.RB4 - LS.RA4) / n
        checks = [
            Check("u_r binnenwand", lambda r: r.at("displx", (LS.RA4, 0, 0)),
                  ref=LS.UR4, tol=LS._tol(0.001, h, (LS.RB4 - LS.RA4) / 8, 3),
                  unit="m", source="Lamé, vlakke rek, zoals LS-04",
                  in_series=False),
            Check("grootste fout aan de binnenwand (r)",
                  lambda r: math.hypot(*r.argmax("errorestimate")[:2]),
                  ref=LS.RA4, tol_abs=h + 1e-9, unit="m",
                  source="grootste spanningsgradiënt bij r = a (Lamé)",
                  in_series=False),
        ]
        if n >= 2:
            checks.append(Check(
                "orde globale schatting", _order("FS-02-h20-n%d" % ns[k - 1]),
                ref=2.0, tol=0.3, unit="",
                source="kwadratische elementen, gladde oplossing: "
                       "energienorm O(h^2)", in_series=False))
        out.append(Case(
            "FS-02-h20-n%d" % n, GROUP,
            "Foutschatter: dikwandige cilinder hex20 %dx%dx1" % (n, 2 * n),
            "lin_static_3d", LS.build_ls04(n, 2 * n), checks, h=h,
            note="Model van LS-04 (symmetry en frictionless op de "
                 "randvlakken, druk binnen, buitenwand vrij)."))
    return out


TESTS += _fs02_cases()


# ------------------------------------------------------------------
# FS-03 Code-Aster-testcase sslv153a via de bridge
# ------------------------------------------------------------------
# Mesh uit sslv153a.mmed (share/aster/tests, Code-Aster 17.4), uitgelezen
# met mesh.getCoordinates()/getConnectivity()/getCells(): knopen 1..11,
# cellen in Code-Aster-volgorde (= Mecway voor 3D, bridge_mesh._ASTER_TYPES).
_S153_XYZ = [(0, 100, 100), (0, 0, 100), (100, 0, 100), (100, 100, 100),
             (0, 100, 0), (0, 0, 0), (100, 0, 0), (100, 100, 0),
             (-50, 100, 50), (-50, 0, 50), (50, -50, 50)]
_S153_VOLUMES = [("tet4", [2, 6, 11, 10]), ("wedge6", [5, 1, 9, 6, 2, 10]),
                 ("pyr5", [2, 6, 7, 3, 11]),
                 ("hex8", [1, 4, 3, 2, 5, 8, 7, 6])]
_S153_FACES = {"CHARPYRA": [[2, 11, 3]], "CHARTETR": [[10, 11, 2]],
               "LIBRE": [[11, 7, 3], [11, 6, 7], [11, 10, 6], [5, 8, 7, 6],
                         [9, 5, 6, 10]],
               "MURPENT": [[1, 5, 9]], "MURHEXA": [[1, 4, 8, 5]],
               "CHARPENT": [[1, 9, 10, 2]], "CHARHEXA": [[1, 2, 3, 4]],
               "LIBRE2": [[3, 7, 8, 4]]}
_S153_PRES = (("CHARTETR", 13.0), ("CHARPYRA", 14.0), ("CHARPENT", 15.0),
              ("CHARHEXA", 16.0))
# TEST_RESU van sslv153a (B.4, primair probleem, NON_REGRESSION)
S153_DZ = -0.018667672103377813
S153_ERREST = 9339.361119203888
# Globale fout van hetzelfde probleem in de .mess van sslv153a (gedraaid
# met Code-Aster 17.4 bij de bouw van de foutschatter)
S153_REL = 51.518


def build_fs03():
    import bridge_mesh
    m = Model("S30")
    m.elset("VOLUME", "Mat")
    for p in _S153_XYZ:
        m.node(*p)
    eids = []
    for shape, nids in _S153_VOLUMES:
        eids.append((m.elem("VOLUME", shape, nids), shape, nids))
    atypes = {"tet4": "TETRA", "wedge6": "PENTA", "pyr5": "PYRAM",
              "hex8": "HEXA"}
    for name, faces in _S153_FACES.items():
        sel = []
        for face in faces:
            for eid, shape, nids in eids:
                table = bridge_mesh._FACE_CORNERS[atypes[shape]]
                fid = [f for f, c in table.items()
                       if set(nids[i] for i in c) == set(face)]
                if fid:
                    sel.append((eid, fid[0]))
                    break
        assert len(sel) == len(faces), name
        m.face_selection(name, sel)
    m.node_selection("MESURE", [6])
    m.mat("Mat", '<mechanical type="Isotropic" youngsmodulus="400000 Pa" '
                 'poissonratio="0.3" />')
    m.add('<fixed selection="MURHEXA" />')
    m.add('<fixed selection="MURPENT" />')
    for name, p in _S153_PRES:
        m.add('<surfacepressure selection="%s"><normal>%g Pa</normal>'
              '</surfacepressure>' % (name, p))
    return m


TESTS.append(Case(
    "FS-03", GROUP, "Foutschatter: Code-Aster-testcase sslv153a",
    "lin_static_3d", build_fs03,
    [Check("DZ knoop MESURE", lambda r: r.node("displz", 6), ref=S153_DZ,
           tol=1e-6, unit="", source="sslv153a TEST_RESU (NON_REGRESSION)"),
     Check("ERREST ERME_NOEU knoop MESURE",
           lambda r: r.node("errorestimate", 6), ref=S153_ERREST, tol=1e-6,
           unit="", source="sslv153a TEST_RESU (NON_REGRESSION)"),
     Check("globale relatieve fout", lambda r: _est(r)["rel"], ref=S153_REL,
           tol=1e-4, unit="%", source="sslv153a .mess (CALC_ERREUR)")],
    general={"elem_hex8": "LINE_QUAD"},
    note="Mesh van sslv153a (11 knopen, hexa/penta/pyram/tetra, door "
         "LINE_QUAD kwadratisch zoals de test); de groepen LIBRE/LIBRE2 "
         "(druk 1e-14 in de test) worden in de bridge vanzelf vrije vlakken."))
