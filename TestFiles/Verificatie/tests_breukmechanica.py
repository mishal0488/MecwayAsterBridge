# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Breukmechanica met een gemeshte scheur (lineair statisch, SCHEUR_-
selecties -> DEFI_FOND_FISS, CALC_G, POST_K1_K2_K3; bridge_fracture.py,
docs/Breukmechanica.md). Toelichting: docs/Verificatie.md hoofdstuk
Breukmechanica (BR).

BR-01  Plaat met een doorgaande middenscheur onder trek, vlakke rek (DZ = 0
       op beide zijvlakken), hex20, beide lippen (de scheur is open zoals
       na Mecway's Disconnect elements: dubbele knopen op de lippen).
       K_I tegen Tada: K = sigma sqrt(pi a) F(a/b), F = sqrt(sec(pi a /
       2b)) (1 - 0,025 (a/b)^2 + 0,06 (a/b)^4) (Tada, Paris & Irwin, The
       Stress Analysis of Cracks Handbook, 3e druk 2000, p. 2.2; nauwkeurig
       0,1 %); G = K^2 (1 - nu^2) / E.
BR-02  Hetzelfde als een half model: alleen de bovenhelft, het ligament
       met een symmetrievoorwaarde, alleen de bovenlip (SYME='OUI').
BR-03  Dikke plaat (dikte B = 4a, half model in de dikte) met vrije
       zijvlakken: in het midden vlakke rek (K_I volgens Tada), naar het
       vrije oppervlak toe daalt K_I (randeffect).
BR-04  Code-Aster-testcase forma07a (schijfvormige scheur in een groot
       blok, kwartmodel, hex20/penta15 met kwartpuntknopen) via de bridge:
       zelfde mesh, lasten en voorwaarden; K1 MAX van POST_K1_K2_K3 gelijk
       aan de TEST_TABLE-waarde.
BR-05  BR-01 met lineaire elementen (hex8, 8 over de dikte).
BR-06  Middenscheur over de volle breedte: twee fronten als twee scheuren
       met een achtervoegsel (_L, _R).
"""
import io
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

import verif_lib as V
from verif_lib import Case, Check, Model, hexmesh

if V.SCRIPTS not in sys.path:
    sys.path.insert(0, V.SCRIPTS)

GROUP = "Breukmechanica"
E, NU = 210e9, 0.3
E_PRIME = E / (1.0 - NU ** 2)          # vlakke rek
STEEL = ('<mechanical type="Isotropic" youngsmodulus="210000 MPa" '
         'poissonratio="0.3" />')


# ------------------------------------------------------------------
# Tabel van de run (<model>_verif_scheur.csv naast de .mess)
# ------------------------------------------------------------------
def _rows(r):
    import bridge_fracture
    if not r.mess:
        raise ValueError("geen .mess")
    rows = bridge_fracture.read_csv(r.mess[:-len(".mess")] +
                                    bridge_fracture.CSV_SUFFIX)
    if not rows:
        raise ValueError("geen tabel %s" % bridge_fracture.CSV_SUFFIX)
    return sorted(rows, key=lambda row: row["s"])


def _with(r, key):
    """Rijen met een waarde voor key (POST_K1_K2_K3 geeft bij een
    kwadratisch front alleen waarden in de hoekknopen)."""
    rows = [q for q in _rows(r) if isinstance(q.get(key), float)]
    if not rows:
        raise ValueError("geen waarden voor %s" % key)
    return rows


def _at(key, frac=0.5):
    """Waarde van key in het frontpunt het dichtst bij een fractie van de
    frontlengte (0 = begin, 0,5 = midden)."""
    def value(r):
        s_end = _rows(r)[-1]["s"]
        row = min(_with(r, key), key=lambda q: abs(q["s"] - frac * s_end))
        return row[key]
    return value


def _at_crack(key, crack):
    """Waarde van key in het middelste punt van scheur crack (kolom
    "crack" van de tabel, bijv. "_L")."""
    def value(r):
        rows = [q for q in _with(r, key) if str(q.get("crack")) == crack]
        if not rows:
            raise ValueError("geen scheur %s in de tabel" % crack)
        return rows[len(rows) // 2][key]
    return value


def _at_z(key, z):
    """Waarde van key in het frontpunt het dichtst bij hoogte z."""
    def value(r):
        return min(_with(r, key), key=lambda q: abs(q["z"] - z))[key]
    return value


def _spread(key):
    """(max - min) / gemiddelde van key langs het front."""
    def value(r):
        vals = [q[key] for q in _with(r, key)]
        return (max(vals) - min(vals)) / (sum(vals) / len(vals))
    return value


def _max(key):
    def value(r):
        return max(q[key] for q in _with(r, key))
    return value


# ------------------------------------------------------------------
# Mesh-hulpmiddelen
# ------------------------------------------------------------------
def _graded(x0, x1, n, h0, at_start=True):
    """n intervallen van x0 naar x1, geometrisch oplopend vanaf h0 aan de
    kant x0 (at_start) of x1. Groeifactor q uit h0 (bisectie)."""
    length = x1 - x0
    lo, hi = 1.0, 3.0
    for _i in range(200):
        q = 0.5 * (lo + hi)
        total = h0 * (q ** n - 1.0) / (q - 1.0) if q != 1.0 else h0 * n
        if total > length:
            hi = q
        else:
            lo = q
    steps = [h0 * q ** i for i in range(n)]
    scale = length / sum(steps)
    steps = [s * scale for s in steps]
    if not at_start:
        steps.reverse()
    xs = [x0]
    for s in steps:
        xs.append(xs[-1] + s)
    xs[-1] = x1
    return xs


def _interp(xs):
    n = len(xs) - 1

    def f(u):
        t = u * n
        k = min(max(int(math.floor(t + 1e-9)), 0), n - 1)
        return xs[k] + (t - k) * (xs[k + 1] - xs[k])
    return f


def _grid(xs, ys, zs):
    fx, fy, fz = _interp(xs), _interp(ys), _interp(zs)

    def f(u, v, w):
        return (fx(u), fy(v), fz(w))
    return f


def _split_node(m, nid):
    """Nieuwe knoop op dezelfde plek (zoals Disconnect elements in Mecway):
    buiten de samenvoeg-index van Model.node."""
    m.nodes.append(m.nodes[nid - 1])
    return len(m.nodes)


def disconnect(m, elset, eids, keep):
    """
    Geef de elementen eids van elset eigen knopen voor alle knopen
    waarvoor keep(x, y, z) False is en die ook door andere elementen
    gebruikt worden: de scheur gaat open (Mecway: de elementen aan een kant
    selecteren, Disconnect elements, daarna het ligament en het front weer
    samenvoegen met Merge nearby nodes). Returns {oud: nieuw}.
    """
    eidset = set(eids)
    used_elsewhere = set()
    for _name, _mat, elems in m.elsets:
        for eid, _shape, nids, _a in elems:
            if eid not in eidset:
                used_elsewhere.update(nids)
    new = {}
    for _name, _mat, elems in m.elsets:
        for k, (eid, shape, nids, attrs) in enumerate(elems):
            if eid not in eidset:
                continue
            out = []
            for n in nids:
                if n in used_elsewhere and not keep(*m.nodes[n - 1]):
                    if n not in new:
                        new[n] = _split_node(m, n)
                    out.append(new[n])
                else:
                    out.append(n)
            elems[k] = (eid, shape, out, attrs)
    return new


# ------------------------------------------------------------------
# BR-01..03 Plaat met doorgaande middenscheur (half model in x)
# ------------------------------------------------------------------
A, B_HALF, H = 0.05, 0.10, 0.30     # halve scheurlengte a, halve breedte b,
#                                     halve hoogte (H/b = 3)
SIGMA = 100e6
ALPHA = A / B_HALF


def tada_f(alpha):
    """Tada (2000) p. 2.2, middenscheur in een strook (Feddersen-vorm met
    correctie); nauwkeurig 0,1 % voor elke a/b."""
    return (math.sqrt(1.0 / math.cos(math.pi * alpha / 2.0)) *
            (1.0 - 0.025 * alpha ** 2 + 0.06 * alpha ** 4))


K_TADA = SIGMA * math.sqrt(math.pi * A) * tada_f(ALPHA)
G_TADA = K_TADA ** 2 / E_PRIME
SRC_K = ("Tada, Paris & Irwin (2000) p. 2.2: K = sigma sqrt(pi a) F(a/b), "
         "a/b = 0,5, F = 1,1867")
SRC_G = "G = K^2 (1 - nu^2) / E (vlakke rek; Irwin)"


def build_plate(lower=True, thickness=None, nz=2, plane_strain=True,
                h_tip=A / 40.0, order=2, open_crack=True):
    """
    Plaat x in [0, b] (symmetrie x = 0, door het scheurmidden), y in [-H, H]
    (lower: beide helften, de lippen op y = 0 voor x < a; anders alleen
    y >= 0 met symmetrie op het ligament), z in [0, t]. plane_strain: DZ =
    0 op z = 0 en z = t; anders symmetrie op z = 0 (middenvlak) en z = t
    vrij (t = halve dikte).
    """
    t = thickness if thickness is not None else A
    xs = _graded(0.0, A, 12, h_tip, at_start=False)[:-1] + \
        _graded(A, B_HALF, 12, h_tip)
    yu = _graded(0.0, H, 16, h_tip)
    if plane_strain:
        zs = [t * k / float(nz) for k in range(nz + 1)]
    else:
        zs = _graded(0.0, t, nz, t / (2.0 * nz) if nz > 1 else t,
                     at_start=False) if nz > 1 else [0.0, t]
    m = Model("S30")
    m.elset("Plaat", "Staal")
    up = hexmesh(m, "Plaat", len(xs) - 1, len(yu) - 1, len(zs) - 1,
                 _grid(xs, yu, zs), order=order)
    tol = 1e-9
    i_crack = [i for i in range(len(xs) - 1) if xs[i + 1] <= A + tol]
    lip_up = [(up["elems"][(i, 0, k)], V.HEX_FACE["v0"]) for i in i_crack
              for k in range(len(zs) - 1)]
    lig_up = [(up["elems"][(i, 0, k)], V.HEX_FACE["v0"])
              for i in range(len(xs) - 1) if i not in i_crack
              for k in range(len(zs) - 1)]
    faces = {"SymX": list(up["faces"]["u0"]), "Boven": up["faces"]["v1"],
             "Z0": list(up["faces"]["w0"]), "Z1": list(up["faces"]["w1"])}
    if lower:
        yl = [-y for y in reversed(yu)]
        lo = hexmesh(m, "Plaat", len(xs) - 1, len(yl) - 1, len(zs) - 1,
                     _grid(xs, yl, zs), order=order)
        lip_lo = [(lo["elems"][(i, len(yl) - 2, k)], V.HEX_FACE["v1"])
                  for i in i_crack for k in range(len(zs) - 1)]
        # Scheur open: de onderste helft krijgt eigen knopen op y = 0 voor
        # x < a (front x = a en ligament blijven gedeeld)
        if open_crack:
            disconnect(m, "Plaat", list(lo["elems"].values()),
                       lambda x, y, z: abs(y) > tol or x >= A - tol)
        faces["SymX"] += lo["faces"]["u0"]
        faces["Onder"] = lo["faces"]["v0"]
        faces["Z0"] += lo["faces"]["w0"]
        faces["Z1"] += lo["faces"]["w1"]
    for name, fl in sorted(faces.items()):
        m.face_selection(name, fl)
    front = [i + 1 for i, p in enumerate(m.nodes)
             if abs(p[0] - A) < tol and abs(p[1]) < tol]
    m.node_selection("SCHEUR_FRONT", front)
    m.face_selection("SCHEUR_LIP_BOVEN", lip_up)
    if lower:
        m.face_selection("SCHEUR_LIP_ONDER", lip_lo)
    else:
        m.face_selection("Ligament", lig_up)
    m.mat("Staal", STEEL)
    m.add('<symmetry selection="SymX"><x>1</x><y>0</y><z>0</z></symmetry>')
    m.add('<traction selection="Boven"><x>0 Pa</x><y>%r Pa</y><z>0 Pa</z>'
          '</traction>' % SIGMA)
    if lower:
        m.add('<traction selection="Onder"><x>0 Pa</x><y>%r Pa</y>'
              '<z>0 Pa</z></traction>' % -SIGMA)
        # Starre verplaatsing in y: een knoop op het ligament
        anchor = [i + 1 for i, p in enumerate(m.nodes)
                  if abs(p[0] - B_HALF) < tol and abs(p[1]) < tol and
                  abs(p[2]) < tol]
        m.node_selection("Anker", anchor)
        m.add('<displacement selection="Anker"><value>0 m</value><x>0</x>'
              '<y>1</y><z>0</z></displacement>')
    else:
        m.add('<symmetry selection="Ligament"><x>0</x><y>1</y><z>0</z>'
              '</symmetry>')
    if plane_strain:
        m.add('<symmetry selection="Z0"><x>0</x><y>0</y><z>1</z></symmetry>')
        m.add('<symmetry selection="Z1"><x>0</x><y>0</y><z>1</z></symmetry>')
    else:
        m.add('<symmetry selection="Z0"><x>0</x><y>0</y><z>1</z></symmetry>')
    return m


def build_two_fronts(h_tip=A / 40.0):
    """
    Plaat over de volle breedte x in [-b, b] (geen symmetrie in x), vlakke
    rek: de middenscheur heeft twee fronten (x = -a en x = a), ingevoerd als
    twee scheuren met achtervoegsel _L en _R, elk met de halve lippen.
    Starre beweging: DX en DY op de lijn x = b, y = 0 en DY op x = -b, y = 0.
    """
    tol = 1e-9
    xr = _graded(0.0, A, 12, h_tip, at_start=False)[:-1] + \
        _graded(A, B_HALF, 12, h_tip)
    xs = [-x for x in reversed(xr)][:-1] + xr
    yu = _graded(0.0, H, 16, h_tip)
    zs = [0.0, A / 2.0, A]
    m = Model("S30")
    m.elset("Plaat", "Staal")
    nu, nv = len(xs) - 1, len(yu) - 1
    up = hexmesh(m, "Plaat", nu, nv, 2, _grid(xs, yu, zs), order=2)
    yl = [-y for y in reversed(yu)]
    lo = hexmesh(m, "Plaat", nu, nv, 2, _grid(xs, yl, zs), order=2)
    disconnect(m, "Plaat", list(lo["elems"].values()),
               lambda x, y, z: abs(y) > tol or abs(x) >= A - tol)
    for side, sgn in (("_L", -1.0), ("_R", 1.0)):
        ic = [i for i in range(nu) if min(xs[i] * sgn, xs[i + 1] * sgn) >=
              -tol and max(abs(xs[i]), abs(xs[i + 1])) <= A + tol]
        m.face_selection("SCHEUR_LIP_BOVEN" + side,
                         [(up["elems"][(i, 0, k)], V.HEX_FACE["v0"])
                          for i in ic for k in range(2)])
        m.face_selection("SCHEUR_LIP_ONDER" + side,
                         [(lo["elems"][(i, nv - 1, k)], V.HEX_FACE["v1"])
                          for i in ic for k in range(2)])
        m.node_selection("SCHEUR_FRONT" + side,
                         [i + 1 for i, p in enumerate(m.nodes)
                          if abs(p[0] - sgn * A) < tol and abs(p[1]) < tol])
    m.face_selection("Boven", up["faces"]["v1"])
    m.face_selection("Onder", lo["faces"]["v0"])
    m.face_selection("Z0", up["faces"]["w0"] + lo["faces"]["w0"])
    m.face_selection("Z1", up["faces"]["w1"] + lo["faces"]["w1"])
    m.mat("Staal", STEEL)
    m.add('<traction selection="Boven"><x>0 Pa</x><y>%r Pa</y><z>0 Pa</z>'
          '</traction>' % SIGMA)
    m.add('<traction selection="Onder"><x>0 Pa</x><y>%r Pa</y><z>0 Pa</z>'
          '</traction>' % -SIGMA)
    for sel in ("Z0", "Z1"):
        m.add('<symmetry selection="%s"><x>0</x><y>0</y><z>1</z></symmetry>'
              % sel)
    for name, xa, dirs in (("AnkerR", B_HALF, ((1, 0), (0, 1))),
                           ("AnkerL", -B_HALF, ((0, 1),))):
        m.node_selection(name, [i + 1 for i, p in enumerate(m.nodes)
                                if abs(p[0] - xa) < tol and abs(p[1]) < tol])
        for dx, dy in dirs:
            # x/y/z van <displacement> is een richting (LIAISON_OBLIQUE),
            # dus per richting een eigen voorwaarde
            m.add('<displacement selection="%s"><value>0 m</value><x>%d</x>'
                  '<y>%d</y><z>0</z></displacement>' % (name, dx, dy))
    return m


def _plate_checks(tol_k=0.03, tol_g=0.02, z=None, label="midden"):
    """K_I (CALC_G) en G tegen Tada (tol_k), G tegen K_I^2/E' (tol_g; de
    theta-methode en de interactie-integraal van CALC_G zijn verschillende
    berekeningen) en POST_K1_K2_K3 als controle (5 %: zonder kwartpunt-
    elementen is de extrapolatie van de lipverplaatsing minder nauwkeurig,
    docs/Breukmechanica.md)."""
    at = _at if z is None else (lambda key: _at_z(key, z))
    return [
        Check("K_I CALC_G %s" % label, at("K1"), ref=K_TADA, tol=tol_k,
              unit="Pa.m^0.5", source=SRC_K),
        Check("G CALC_G %s" % label, at("G"), ref=G_TADA, tol=tol_k,
              unit="J/m2", source=SRC_K + "; " + SRC_G),
        Check("G / (K_I^2/E') CALC_G %s" % label,
              (lambda r: at("G")(r) / (at("K1")(r) ** 2 / E_PRIME)),
              ref=1.0, tol=tol_g, unit="",
              source=SRC_G + "; G theta-methode, K interactie-integraal"),
        Check("K_I POST_K1_K2_K3 %s" % label, at("K1_POST_K"), ref=K_TADA,
              tol=0.05, unit="Pa.m^0.5",
              source=SRC_K + " (controlemethode, hex20 zonder kwartpunten)"),
    ]


TESTS = [
    Case("BR-01", GROUP, "Middenscheur in een plaat, vlakke rek, beide "
         "lippen (hex20)", "lin_static_3d",
         lambda: build_plate(lower=True),
         _plate_checks() + [
             Check("spreiding K_I langs het front", _spread("K1"), ref=0.0,
                   tol_abs=0.01, unit="",
                   source="vlakke rek: K_I constant langs het front"),
             Check("K_II midden", _at("K2"), ref=0.0,
                   tol_abs=0.01 * K_TADA, unit="Pa.m^0.5",
                   source="symmetrische belasting: modus I")],
         note="Half model in x (symmetrie door het scheurmidden), a/b = "
              "0,5, H/b = 3; de onderhelft heeft eigen knopen op de lippen "
              "(dubbele knopen, zoals Disconnect elements)."),
    Case("BR-02", GROUP, "Middenscheur, vlakke rek, kwartmodel met "
         "symmetrie op het ligament (SYME='OUI')", "lin_static_3d",
         lambda: build_plate(lower=False), _plate_checks(),
         note="Alleen y >= 0; ligament symmetry, alleen SCHEUR_LIP_BOVEN."),
    Case("BR-03", GROUP, "Middenscheur in een dikke plaat (B = 4a) met vrije "
         "zijvlakken", "lin_static_3d",
         lambda: build_plate(lower=True, thickness=2 * A, nz=10,
                             plane_strain=False),
         _plate_checks(z=0.0, label="middenvlak") + [
             Check("K_I CALC_G op a van het oppervlak", _at_z("K1", A),
                   ref=K_TADA, tol=0.03, unit="Pa.m^0.5", source=SRC_K)],
         note="Half model in de dikte (symmetrie op het middenvlak z = 0, "
              "vrij oppervlak z = 2a)."),
    Case("BR-06", GROUP, "Middenscheur over de volle breedte: twee fronten "
         "als twee scheuren (_L, _R)", "lin_static_3d",
         lambda: build_two_fronts(),
         [Check("K_I CALC_G midden %s" % side,
                (lambda r, s=side: _at_crack("K1", s)(r)), ref=K_TADA,
                tol=0.03, unit="Pa.m^0.5", source=SRC_K)
          for side in ("_L", "_R")] +
         [Check("G CALC_G midden %s" % side,
                (lambda r, s=side: _at_crack("G", s)(r)), ref=G_TADA,
                tol=0.03, unit="J/m2", source=SRC_K + "; " + SRC_G)
          for side in ("_L", "_R")],
         note="Plaat x in [-b, b] zonder symmetrie in x; per front een eigen "
              "achtervoegsel (SCHEUR_FRONT_L met SCHEUR_LIP_BOVEN_L/_ONDER_L "
              "over de halve lippen, idem _R)."),
    Case("BR-05", GROUP, "Middenscheur, vlakke rek, lineaire elementen "
         "(hex8, 8 over de dikte)", "lin_static_3d",
         lambda: build_plate(lower=True, nz=8, order=1),
         _plate_checks()[:3],
         general={"elem_hex8": "3D"},
         note="Zelfde plaat als BR-01 met hex8 (volledig geintegreerd). "
              "POST_K1_K2_K3 geeft met hex8 ~10 % te laag (geen check; "
              "docs/Breukmechanica.md)."),
]


# ------------------------------------------------------------------
# BR-04 Code-Aster-testcase forma07a via de bridge
# ------------------------------------------------------------------
# forma07a (share/aster/tests, Code-Aster 17.4): schijfvormige scheur
# (straal 2 m) in een blok, kwartmodel; FACE_AV/FACE_LAT DNOR = 0, knoop D
# DZ = 0, FORCE_FACE FZ = +-1e6 Pa op FACE_SUP/FACE_INF, E = 2e11, nu = 0,3,
# kwartpuntknopen (MODI_MAILLAGE NOEUD_QUART op het front LFF).
F07_K1_CALC = 1607011.55098     # TEST_TABLE K1 MAX, VALE_CALC (POST_K1_K2_K3)
F07_K1_REFE = 1.5957e6          # VALE_REFE (analytisch: 2 sigma sqrt(a/pi))
F07_G_REF = 11.58648            # G_REF in forma07a.comm (oneindig medium)
F07_R_INF, F07_R_SUP = 0.2, 0.5  # zoals CALC_G in forma07a


def _aster_tests_dir():
    import bridge_utils
    exe = bridge_utils.exe_aster_path()
    return os.path.join(os.path.dirname(os.path.dirname(exe)), "share",
                        "aster", "tests")


_CONVERT = """DEBUT()
MA = LIRE_MAILLAGE(FORMAT='MED', UNITE=20)
MA = MODI_MAILLAGE(reuse=MA, MAILLAGE=MA,
                   MODI_MAILLE=_F(OPTION='NOEUD_QUART', GROUP_MA_FOND='LFF'))
IMPR_RESU(FORMAT='ASTER', UNITE=21, RESU=_F(MAILLAGE=MA))
FIN()
"""


def _forma07a_mail():
    """Mesh van forma07a met de kwartpuntknopen in ASTER-formaat: met
    Code-Aster omgezet uit forma07a.mmed (in een tijdelijke map, bewaard)."""
    import bridge_utils
    work = os.path.join(tempfile.gettempdir(), "mab_forma07a")
    mail = os.path.join(work, "forma07a_q.mail")
    if os.path.isfile(mail) and os.path.getsize(mail) > 0:
        return mail
    if os.path.isdir(work):
        shutil.rmtree(work)
    os.makedirs(work)
    shutil.copy(os.path.join(_aster_tests_dir(), "forma07a.mmed"),
                os.path.join(work, "forma07a.mmed"))
    io.open(os.path.join(work, "conv.comm"), "w").write(_CONVERT)
    io.open(os.path.join(work, "conv.export"), "w").write(
        "P time_limit 120\nP memory_limit 1024\nP ncpus 1\nP mpi_nbcpu 1\n"
        "F comm %s D 1\nF mmed %s D 20\nF mail %s R 21\nF mess %s R 6\n" % (
            os.path.join(work, "conv.comm"),
            os.path.join(work, "forma07a.mmed"), mail,
            os.path.join(work, "conv.mess")))
    with open(os.path.join(work, "console.txt"), "w") as f:
        subprocess.call([bridge_utils.exe_aster_path(),
                         os.path.join(work, "conv.export")], stdout=f,
                        stderr=subprocess.STDOUT, cwd=work)
    if not os.path.isfile(mail):
        raise RuntimeError("forma07a-mesh niet omgezet (zie %s)" % work)
    return mail


_ASTER_N = {"SEG3": 3, "TRIA6": 6, "QUAD8": 8, "PENTA15": 15, "HEXA20": 20}


def read_aster_mail(path):
    """(knopen {naam: (x, y, z)} in volgorde, cellen {naam: (type,
    [knoopnamen])} in volgorde, groepen {('GROUP_MA'|'GROUP_NO', naam):
    [leden]}) uit een mesh in ASTER-formaat (IMPR_RESU FORMAT='ASTER')."""
    from collections import OrderedDict
    nodes, cells, groups = OrderedDict(), OrderedDict(), OrderedDict()
    section, tokens = None, []

    def flush():
        if section == "COOR_3D":
            for k in range(0, len(tokens), 4):
                nodes[tokens[k]] = tuple(float(v) for v in
                                         tokens[k + 1:k + 4])
        elif section in _ASTER_N:
            n = _ASTER_N[section] + 1
            for k in range(0, len(tokens), n):
                cells[tokens[k]] = (section, tokens[k + 1:k + n])
        elif section in ("GROUP_MA", "GROUP_NO"):
            groups[(section, tokens[0])] = tokens[1:]
    for raw in io.open(path, encoding="latin-1"):
        line = raw.strip()
        if not line or line.startswith("%"):
            continue
        if section is None:
            section, tokens = line.split()[0], []
            continue
        if line == "FINSF":
            flush()
            section = None
            continue
        if section != "TITRE":
            tokens += line.split()
    return nodes, cells, groups


def build_forma07a():
    import bridge_mesh
    nodes, cells, groups = read_aster_mail(_forma07a_mail())
    m = Model("S30")
    index = {}
    for name, xyz in nodes.items():
        m.nodes.append(xyz)
        index[name] = len(m.nodes)
    m.elset("Blok", "Staal")
    shapes = {"HEXA20": ("hex20", "HEXA"), "PENTA15": ("wedge15", "PENTA")}
    face_of = {}
    for name, (atype, cn) in cells.items():
        if atype not in shapes:
            continue
        shape, base = shapes[atype]
        nids = [index[n] for n in cn]
        eid = m.elem("Blok", shape, nids)
        for fid, corners in bridge_mesh._FACE_CORNERS[base].items():
            face_of[frozenset(nids[i] for i in corners)] = (eid, fid)

    def faces(group):
        out = []
        for cname in groups[("GROUP_MA", group)]:
            atype, cn = cells[cname]
            nc = 4 if atype.startswith("QUAD") else 3
            out.append(face_of[frozenset(index[n] for n in cn[:nc])])
        return out
    for group, sel in (("FACE_AV", "FACE_AV"), ("FACE_LAT", "FACE_LAT"),
                       ("FACE_SUP", "FACE_SUP"), ("FACE_INF", "FACE_INF"),
                       ("LEV_SUP", "SCHEUR_LIP_BOVEN"),
                       ("LEV_INF", "SCHEUR_LIP_ONDER")):
        m.face_selection(sel, faces(group))
    front = []
    for cname in groups[("GROUP_MA", "LFF")]:
        front += [index[n] for n in cells[cname][1]]
    m.node_selection("SCHEUR_FRONT", front)
    m.node_selection("D", [index[n] for n in groups[("GROUP_NO", "D")]])
    m.mat("Staal", '<mechanical type="Isotropic" youngsmodulus="200000 MPa" '
                   'poissonratio="0.3" />')
    for sel in ("FACE_AV", "FACE_LAT"):
        m.add('<frictionless selection="%s"><displacement>0 m'
              '</displacement></frictionless>' % sel)
    m.add('<displacement selection="D"><value>0 m</value><x>0</x><y>0</y>'
          '<z>1</z></displacement>')
    m.add('<traction selection="FACE_SUP"><x>0 Pa</x><y>0 Pa</y>'
          '<z>1000000 Pa</z></traction>')
    m.add('<traction selection="FACE_INF"><x>0 Pa</x><y>0 Pa</y>'
          '<z>-1000000 Pa</z></traction>')
    # R_INF/R_SUP van forma07a per model (docs/ASTER_regels_voorbeelden.md
    # voorbeeld 16)
    m.add("<ccxmodeldefinition>** ASTER-OPTIE fracture_r_inf = %r\n"
          "** ASTER-OPTIE fracture_r_sup = %r</ccxmodeldefinition>"
          % (F07_R_INF, F07_R_SUP))
    return m


TESTS.append(Case(
    "BR-04", GROUP, "Code-Aster-testcase forma07a (schijfvormige scheur)",
    "lin_static_3d", build_forma07a,
    [Check("K1 MAX POST_K1_K2_K3", _max("K1_POST_K"), ref=F07_K1_CALC,
           tol=0.01, unit="Pa.m^0.5",
           source="forma07a TEST_TABLE VALE_CALC (NON_REGRESSION)"),
     Check("K1 MAX POST_K1_K2_K3 (analytisch)", _max("K1_POST_K"),
           ref=F07_K1_REFE, tol=0.008, unit="Pa.m^0.5",
           source="forma07a TEST_TABLE VALE_REFE (ANALYTIQUE, 0,8 %): "
                  "K = 2 sigma sqrt(a/pi)"),
     Check("K_I CALC_G midden (FORM_FISS CERCLE)", _at("K1"),
           ref=F07_K1_REFE, tol=0.01, unit="Pa.m^0.5",
           source="analytisch 2 sigma sqrt(a/pi) (forma07a VALE_REFE)"),
     Check("G CALC_G midden", _at("G"), ref=F07_G_REF, tol=0.015,
           unit="J/m2", source="forma07a G_REF (K^2 (1 - nu^2) / E)")],
    note="Mesh uit forma07a.mmed met de kwartpuntknopen (NOEUD_QUART) al in "
         "de .liml; LEV_SUP/LEV_INF/LFF als SCHEUR_-selecties, FACE_AV/"
         "FACE_LAT frictionless, R_INF/R_SUP = 0,2/0,5 m zoals de test (met "
         "#ASTER-OPTIE)."))

TESTS.sort(key=lambda c: c.id)
