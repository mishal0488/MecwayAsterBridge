# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Akoestische holtemodes (Mecway A30 -> PHENOMENE='ACOUSTIQUE',
CALC_MODES; bridge-analyse acoustic_modes_3d, bridge_acoustic.py).
Toelichting: docs/Verificatie.md par. 7 (AK) en docs/Akoestiek.md.

Referenties:
  - buis met lengte L en harde wanden (vlakke golven, 1D Helmholtz
    p'' + (w/c)^2 p = 0): gesloten-gesloten (dp/dx = 0 aan beide einden)
    f_n = n c / (2L), n = 0, 1, 2, ...; open-gesloten (p = 0 aan een
    einde) f_n = (2n - 1) c / (4L); open-open f_n = n c / (2L), n >= 1
    (afgeleid uit de randvoorwaarden; standaardresultaat, o.a. Kinsler,
    Frey, Coppens & Sanders, Fundamentals of Acoustics, niet lokaal
    ingezien);
  - rechthoekige ruimte Lx x Ly x Lz met harde wanden:
    f_lmn = c/2 sqrt((l/Lx)^2 + (m/Ly)^2 + (n/Lz)^2) (scheiding van
    variabelen met cos-modes; zelfde bron);
  - buis van twee media (c1 over L1, c2 over L2, gelijke dichtheid),
    harde einden: p = A cos(k1 x) en B cos(k2 (L - x)), druk en dp/dx
    continu op het grensvlak: k1 tan(k1 L1) + k2 tan(k2 L2) = 0 (zelf
    afgeleid; numeriek opgelost in _two_media_roots);
  - Mecway 33 headless (TestFiles/Tools/mw_solve.py, interne solver) op
    dezelfde .liml (modellen/AK-*.liml); waarden hieronder in MECWAY.
"""
import math
from verif_lib import Case, Check, Model, hexmesh, box

GROUP = "Akoestiek"
C_AIR = 340.0
AIR = ['<mechanical type="Isotropic" speedofsound="340 m/s" />',
       u'<density type="Simple" density="1.2 kg/m³" />']

# Mecway 33 headless (mw_solve.py) op modellen/AK-*.liml, 7 okt 2026:
# frequenties [Hz] per mode en de grootste |p| van mode 1 en 2 (AK-05: met
# twee materialen geeft Mecway geen oplossing).
MECWAY = {
    "AK-01": {"f": (8.68850761007296e-05, 170.0000717803, 340.002288125291,
                    510.017264978249, 680.072115195093),
              "pmax": {1: 3400.000000000438, 2: 4808.330166125849}},
    "AK-02": {"f": (85.0000022453305, 255.000544204582, 425.006962736163,
                    595.037163275001, 765.129266801552),
              "pmax": {1: 4808.326365995432, 2: 4808.346561842661}},
    "AK-03": {"f": (170.02796525741, 340.223754974219, 510.75535633919,
                    681.79107331688),
              "pmax": {1: 4809.9082481448995, 2: 4805.1571505164975}},
    "AK-04": {"f": (0.00364451821145843, 170.001144062646, 242.86391388667,
                    296.451172986836, 340.036057636623, 417.861076849524,
                    425.108815966069, 457.840781140722),
              "pmax": {1: 642.539604372501, 2: 908.7003750294504}},
    "AK-06": {"f": (3.88561872782948e-05, 170.040644855932, 340.320027206749,
                    511.053488072442, 682.418851692191),
              "pmax": {1: 3400.0000000003733, 2: 4820.042226939056}},
    "AK-07": {"f": (0.0, 170.03832332782, 340.297526952729, 510.958428953003,
                    682.142726725895),
              "pmax": {1: 3399.9999999998663, 2: 4826.2850659891665}},
}


def _mecway_checks(case_id, n, first=1):
    ref = MECWAY.get(case_id)
    if not ref:
        return []
    out = [Check("f%d t.o.v. Mecway" % k, lambda r, k=k: r.mode(k),
                 ref=ref["f"][k - 1], tol=0.01, unit="Hz", in_series=False,
                 source="Mecway 33 headless (mw_solve.py), zelfde .liml")
           for k in range(first, n + 1)]
    for k, pmax in sorted(ref.get("pmax", {}).items()):
        out.append(Check("max|p| mode %d t.o.v. Mecway" % k,
                         lambda r, k=k: _pmax(r, k), ref=pmax, tol=0.001,
                         unit="Pa", in_series=False,
                         source="Mecway 33 headless: normering p^T M p = 1, "
                         "M = int N N / c^2 dV"))
    return out


def _pmax(r, k):
    """Grootste |p| van mode k (de stappen van een modaal resultaat hebben
    t = mode-id)."""
    return max(abs(v) for v in r.step(float(k))["pressure"].values())


def _f_checks(refs, tol, source, first=1):
    out = []
    for k, f in enumerate(refs, start=first):
        if f == 0.0:
            # constante druk: numeriek 0 (|f| van orde 1e-5 Hz)
            out.append(Check("f%d (constante druk)" % k,
                             lambda r, k=k: r.mode(k), ref=0.0,
                             tol_abs=0.01, unit="Hz", source=source))
        else:
            out.append(Check("f%d" % k, lambda r, k=k: r.mode(k), ref=f,
                             tol=tol, unit="Hz", source=source))
    return out


# ------------------------------------------------------------------
# Buis L = 1 m langs X, doorsnede 0,1 x 0,1 m (eerste dwarsmode bij
# c / (2 x 0,1) = 1700 Hz, ver boven de gevraagde modes).
# ------------------------------------------------------------------
LT, WT = 1.0, 0.1


def build_tube(order, n, ends, maxmode):
    """ends: "gg" gesloten-gesloten, "og" p = 0 op x = L (knoopselectie),
    "oo" p = 0 op beide eindvlakken (vlakselecties)."""
    def build():
        m = Model("A30", attrs={"maxmode": maxmode})
        m.elset("Lucht", "Lucht")
        h = hexmesh(m, "Lucht", n, 1, 1, box(0, LT, 0, WT, 0, WT),
                    order=order)
        m.mat("Lucht", *AIR)
        if ends == "og":
            m.node_selection("Open", sorted(h["nodes"]["u1"]))
            m.add('<fluidpressure selection="Open"><value>0 Pa</value>'
                  '</fluidpressure>')
        elif ends == "oo":
            for side, name in (("u0", "Open0"), ("u1", "OpenL")):
                m.face_selection(name, h["faces"][side])
                m.add('<fluidpressure selection="%s"><value>0 Pa</value>'
                      '</fluidpressure>' % name)
        return m
    return build


TUBE = "buis L = 1 m, c = 340 m/s, vlakke golven (1D Helmholtz)"
TESTS = [
    Case("AK-01", GROUP, "Buis gesloten-gesloten, hex20 20x1x1",
         "acoustic_modes_3d", build_tube(2, 20, "gg", 5),
         _f_checks([0.0] + [n * C_AIR / (2 * LT) for n in (1, 2, 3, 4)],
                   0.005, "f_n = n c / (2L); " + TUBE) +
         [Check("max|p| mode 1 = c / sqrt(V)", lambda r: _pmax(r, 1),
                ref=C_AIR / math.sqrt(LT * WT * WT), tol=1e-6, unit="Pa",
                in_series=False,
                source="normering p^T M p = 1 met M = int N N / c^2 dV "
                       "(zoals Mecway): constante mode p = c / sqrt(V)")] +
         _mecway_checks("AK-01", 5, first=2),
         note="Harde wanden, geen drukknoop: mode 1 = constante druk "
              "(0 Hz); de bridge zoekt rond een negatieve frequentie."),
    Case("AK-02", GROUP, "Buis open-gesloten (p = 0 op x = L), hex20 20x1x1",
         "acoustic_modes_3d", build_tube(2, 20, "og", 5),
         _f_checks([(2 * n - 1) * C_AIR / (4 * LT) for n in range(1, 6)],
                   0.005, "f_n = (2n - 1) c / (4L); " + TUBE) +
         _mecway_checks("AK-02", 5),
         note="fluid pressure 0 op een knoopselectie."),
    Case("AK-03", GROUP, "Buis open-open (p = 0 op beide eindvlakken), hex8 "
         "50x1x1", "acoustic_modes_3d", build_tube(1, 50, "oo", 4),
         _f_checks([n * C_AIR / (2 * LT) for n in range(1, 5)], 0.005,
                   "f_n = n c / (2L), n >= 1; " + TUBE) +
         _mecway_checks("AK-03", 4),
         note="hex8 (lineair, geen incompatibele modes nodig: de druk is "
              "een scalair veld); fluid pressure 0 op vlakselecties."),
]


# ------------------------------------------------------------------
# AK-04 Rechthoekige ruimte 1,0 x 0,7 x 0,4 m, harde wanden, hex20
# 10 x 7 x 4: ook schuine en dwarse modes.
# ------------------------------------------------------------------
LX, LY, LZ = 1.0, 0.7, 0.4


def room_modes(count):
    out = []
    for l in range(4):
        for mm in range(4):
            for n in range(3):
                f = C_AIR / 2.0 * math.sqrt((l / LX) ** 2 + (mm / LY) ** 2 +
                                            (n / LZ) ** 2)
                out.append((f, (l, mm, n)))
    return sorted(out)[:count]


def build_room():
    m = Model("A30", attrs={"maxmode": 8})
    m.elset("Lucht", "Lucht")
    hexmesh(m, "Lucht", 10, 7, 4, box(0, LX, 0, LY, 0, LZ), order=2)
    m.mat("Lucht", *AIR)
    return m


TESTS.append(Case(
    "AK-04", GROUP, "Rechthoekige ruimte 1,0 x 0,7 x 0,4 m, hex20 10x7x4",
    "acoustic_modes_3d", build_room,
    [Check("f%d (constante druk)" % 1, lambda r: r.mode(1), ref=0.0,
           tol_abs=0.01, unit="Hz", source="harde wanden")] +
    [Check("f%d mode %s" % (k + 1, "".join(str(i) for i in lmn)),
           lambda r, k=k: r.mode(k + 1), ref=f, tol=0.005, unit="Hz",
           source="f_lmn = c/2 sqrt((l/Lx)^2 + (m/Ly)^2 + (n/Lz)^2) "
                  "(harde wanden, cos-modes)")
     for k, (f, lmn) in enumerate(room_modes(8)) if k > 0] +
    _mecway_checks("AK-04", 8, first=2),
    note="Harde wanden; de modes (l, m, n) volgen uit de afmetingen."))


# ------------------------------------------------------------------
# AK-05 Buis van twee media (c 340 en 680 m/s, elk 0,5 m, geen
# dichtheid): Code-Aster koppelt ze op het grensvlak. Mecway zelf rekent
# akoestiek alleen met een materiaal (headless: geen oplossing).
# ------------------------------------------------------------------
C2 = 680.0


def _two_media_roots(count, c1=C_AIR, c2=C2, l1=0.5, l2=0.5):
    """Nulpunten f > 0 van k1 sin(k1 l1) cos(k2 l2) + k2 sin(k2 l2)
    cos(k1 l1) (vorm zonder polen van k1 tan + k2 tan = 0)."""
    def g(f):
        w = 2 * math.pi * f
        k1, k2 = w / c1, w / c2
        return (k1 * math.sin(k1 * l1) * math.cos(k2 * l2) +
                k2 * math.sin(k2 * l2) * math.cos(k1 * l1))
    out, f, prev = [], 1.0, g(1.0)
    while len(out) < count:
        f2 = f + 0.1
        v = g(f2)
        if prev * v < 0.0:
            a, b = f, f2
            for _i in range(80):
                mid = 0.5 * (a + b)
                if g(a) * g(mid) <= 0.0:
                    b = mid
                else:
                    a = mid
            out.append(0.5 * (a + b))
        f, prev = f2, v
    return out


def build_two_media():
    m = Model("A30", attrs={"maxmode": 5})
    m.elset("Lucht1", "Medium1")
    m.elset("Lucht2", "Medium2")
    hexmesh(m, "Lucht1", 10, 1, 1, box(0, 0.5, 0, WT, 0, WT), order=2)
    hexmesh(m, "Lucht2", 10, 1, 1, box(0.5, 1.0, 0, WT, 0, WT), order=2)
    m.mat("Medium1", '<mechanical type="Isotropic" speedofsound="340 m/s" />')
    m.mat("Medium2", '<mechanical type="Isotropic" speedofsound="680 m/s" />')
    return m


TESTS.append(Case(
    "AK-05", GROUP, "Buis van twee media (c = 340 / 680 m/s), hex20 20x1x1",
    "acoustic_modes_3d", build_two_media,
    _f_checks([0.0] + _two_media_roots(4), 0.005,
              "k1 tan(k1 L1) + k2 tan(k2 L2) = 0 (druk en dp/dx continu, "
              "gelijke dichtheid)"),
    note="Twee materialen zonder dichtheid; Mecway zelf weigert meer dan "
         "een akoestisch materiaal."))


# ------------------------------------------------------------------
# AK-06 / AK-07 De buis van AK-01 met wedge6 en tet4: elk hex8-blok van
# 40 x 1 x 1 gesplitst in 2 wedges of 6 tetraeders (Kuhn). Lineaire
# elementen: discretisatiefout als hex8 (AK-03).
# ------------------------------------------------------------------
_CUBE = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1),
         (1, 1, 1), (0, 1, 1)]
_WEDGES = [(0, 1, 2, 4, 5, 6), (0, 2, 3, 4, 6, 7)]
_TETS = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6), (0, 7, 4, 6), (0, 4, 5, 6),
         (0, 5, 1, 6)]


def _tet_volume(p):
    a, b, c = [[p[i][k] - p[0][k] for k in range(3)] for i in (1, 2, 3)]
    return (a[0] * (b[1] * c[2] - b[2] * c[1]) -
            a[1] * (b[0] * c[2] - b[2] * c[0]) +
            a[2] * (b[0] * c[1] - b[1] * c[0]))


def build_split_tube(shape, n=40):
    def build():
        m = Model("A30", attrs={"maxmode": 5})
        m.elset("Lucht", "Lucht")
        for i in range(n):
            pts = [(LT * (i + u) / n, WT * v, WT * w) for u, v, w in _CUBE]
            nids = [m.node(*p) for p in pts]
            for cell in (_WEDGES if shape == "wedge6" else _TETS):
                cell = list(cell)
                if shape == "tet4" and _tet_volume([pts[k] for k in cell]) < 0:
                    cell[1], cell[2] = cell[2], cell[1]
                m.elem("Lucht", shape, [nids[k] for k in cell])
        m.mat("Lucht", *AIR)
        return m
    return build


for _cid, _shape in (("AK-06", "wedge6"), ("AK-07", "tet4")):
    TESTS.append(Case(
        _cid, GROUP, "Buis gesloten-gesloten, %s (40 hex8-blokken "
        "gesplitst)" % _shape, "acoustic_modes_3d", build_split_tube(_shape),
        _f_checks([0.0] + [n * C_AIR / (2 * LT) for n in (1, 2, 3, 4)],
                  0.005, "f_n = n c / (2L); " + TUBE) +
        _mecway_checks(_cid, 5, first=2),
        note="Lineaire %s; de fout groeit als (kh)^2, zoals bij hex8." % _shape))
