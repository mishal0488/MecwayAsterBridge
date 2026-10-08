# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Modaal met voorspanning (Mecway V30 met belastingen -> MECA_STATIQUE,
SIEF_ELGA, RIGI_GEOM (+ RIGI_ROTA), CALC_MODES met K + Kg; bridge-analyse
modal_3d, instelling modal_prestress). Toelichting: docs/Verificatie.md
par. 7 (MP) en docs/Modaal_en_knik.md.

Referenties:
  - scharnierend opgelegde Euler-balk onder axiale kracht P (trek > 0):
    omega_n^2 = (n pi/L)^4 EI/(rho A) + (n pi/L)^2 P/(rho A), dus
    f_n = f_n0 sqrt(1 + P/Pcr_n) met Pcr_n = n^2 pi^2 EI/L^2 (eigenvorm
    sin(n pi x/L); dezelfde vergelijking geeft de VALE_REFE 'ANALYTIQUE'
    van Code-Aster-testcase sdll101a "corde a piano": nagerekend 4,9713
    en 5,0728 Hz voor F = 0 en 10 N);
  - vrij opgelegde vierkante Kirchhoff-plaat onder gelijkmatige biaxiale
    membraankracht N (druk > 0): rho t omega_mn^2 = pi^4 D (m^2+n^2)^2/a^4
    - N pi^2 (m^2+n^2)/a^2 (Navier-vorm sin sin; grensgevallen N = 0:
    Leissa, Vibration of Plates, NASA SP-160, 1969, par. 4.1; omega = 0
    bij N = Ncr = 2 pi^2 D/a^2: Timoshenko & Gere 1961, par. 9.2);
  - Mecway 33 headless (TestFiles/Tools/mw_solve.py, interne solver) op
    dezelfde .liml (modellen/MP-*.liml); waarden hieronder in MECWAY.
"""
import math
from verif_lib import Case, Check, Model, linemesh, hexmesh, quadmesh, box, \
    rect

GROUP = "Modaal (voorspanning)"
E, NU, RHO, ALPHA_T = 210e9, 0.3, 7850.0, 1.2e-5
STEEL = ('<mechanical type="Isotropic" youngsmodulus="210000 MPa" '
         'poissonratio="0.3" />')
STEEL_T = ('<mechanical type="Isotropic" youngsmodulus="210000 MPa" '
           'poissonratio="0.3" thermalexpansion="%g 1/K" />' % ALPHA_T)
DENS = '<density type="Simple" density="7850 kg/m³" />'

# Mecway 33 headless (mw_solve.py) op modellen/MP-*.liml, 6 okt 2026.
MECWAY = {
    "MP-01": (66.3359813450715, 84.5621227004543, 209.773880217595,
              296.665266081589),
    "MP-02": (33.168011686957, 62.0516492413176, 175.50982514665,
              273.512153297547),
    "MP-03": (66.3049688805252, 84.4730803761122),
    "MP-04": (20.9030140819913, 42.6347744242016, 128.573053550737,
              176.786288528252),
    "MP-05": (86.9719003542034, 275.43032880853, 275.43032880868),
    "MP-07": (20.9751008350118, 42.651647148566, 128.551043402919,
              176.317385505606),
}


def _mode_check(k, ref, tol, source, name=None):
    return Check(name or "f%d" % k, lambda r, k=k: r.mode(k), ref=ref,
                 tol=tol, unit="Hz", source=source)


def _mecway_checks(case_id, n, tol=0.01):
    ref = MECWAY.get(case_id)
    if not ref:
        return []
    return [Check("f%d t.o.v. Mecway" % (k + 1), lambda r, k=k: r.mode(k + 1),
                  ref=ref[k], tol=tol, unit="Hz", in_series=False,
                  source="Mecway 33 headless (mw_solve.py), zelfde .liml")
            for k in range(n)]


# ------------------------------------------------------------------
# Balk L = 1 m langs X, doorsnede a x b = 0,02 x 0,03 m (a in Y = recta),
# scharnierend: knik/buiging in Y (zwakke as, I_z) en in Z (I_y = 2,25 I_z).
# ------------------------------------------------------------------
L1, A1, B1 = 1.0, 0.02, 0.03
AREA1 = A1 * B1
IZ1 = B1 * A1 ** 3 / 12.0
IY1 = A1 * B1 ** 3 / 12.0
PCR1 = math.pi ** 2 * E * IZ1 / L1 ** 2     # 41452 N, mode Y1


def f_pp(n, inertia, p):
    """Scharnierende Euler-balk met axiale kracht p (trek > 0)."""
    k = n * math.pi / L1
    w2 = (k ** 4 * E * inertia + k ** 2 * p) / (RHO * AREA1)
    return math.sqrt(w2) / (2 * math.pi)


def modes_pp(p, count=4):
    """De laagste modes (frequentie, naam) van de scharnierende balk."""
    return sorted([(f_pp(n, IZ1, p), "Y%d" % n) for n in (1, 2, 3)] +
                  [(f_pp(n, IY1, p), "Z%d" % n) for n in (1, 2, 3)])[:count]


def build_beam(p):
    def build():
        m = Model("V30", attrs={"maxmode": 4})
        m.elset("Balk", "Staal")
        _e, nids = linemesh(m, "Balk", 20, (0, 0, 0), (L1, 0, 0))
        m.mat("Staal", '<geometric type="RectBar" recta="%g m" rectb="%g m" '
              '/>' % (A1, B1), STEEL, DENS)
        m.node_selection("Voet", [nids[0]])
        m.node_selection("Kop", [nids[-1]])
        # Scharnierend zoals BU-02: voet X, Y, Z en torsie vast; kop Y, Z.
        m.add('<displacement selection="Voet"><value>0 m</value><x>1</x>'
              '<y>0</y><z>0</z></displacement>')
        for sel in ("Voet", "Kop"):
            for d in ("<x>0</x><y>1</y><z>0</z>", "<x>0</x><y>0</y><z>1</z>"):
                m.add('<displacement selection="%s"><value>0 m</value>%s'
                      '</displacement>' % (sel, d))
        m.add('<noderotation selection="Voet"><value>0 rad</value><x>1</x>'
              '<y>0</y><z>0</z></noderotation>')
        m.add('<force selection="Kop"><x>%r N</x><y>0 N</y><z>0 N</z>'
              '</force>' % p)
        return m
    return build


def _beam_case(cid, title, p, note, settings=None, ref_p=None):
    ref_p = p if ref_p is None else ref_p
    checks = [_mode_check(k + 1, f, 0.002,
                          "Euler-balk scharnierend, P = %+.0f N: f_n0 "
                          "sqrt(1 + P/Pcr_n), mode %s (sdll101a-vergelijking)"
                          % (ref_p, plane))
              for k, (f, plane) in enumerate(modes_pp(ref_p))]
    if not settings:
        checks += _mecway_checks(cid, 4)
    return Case(cid, GROUP, title, "modal_3d", build_beam(p), checks,
                settings=settings, note=note)


TESTS = [
    _beam_case("MP-01", "Scharnierende balk onder trek P = Pcr1 "
               "(20 balkelementen)", PCR1,
               "POU_D_E; trek = 1 x de eerste kniklast: f1 = sqrt(2) f10."),
    _beam_case("MP-02", "Scharnierende balk onder druk P = 0,5 Pcr1 "
               "(20 balkelementen)", -0.5 * PCR1,
               "POU_D_E; druk = 0,5 x de eerste kniklast: f1 = f10/sqrt(2)."),
    # Instelling uit: de belasting telt niet mee (de onbelaste frequenties).
    _beam_case("MP-06", "Scharnierende balk onder trek, voorspanning uit",
               PCR1, "modal_prestress = False: frequenties zonder last.",
               settings={"modal_prestress": False}, ref_p=0.0),
]


# ------------------------------------------------------------------
# MP-03 Dezelfde balk als volume (hex20 50 x 2 x 4), scharnierend op de
# lijn z = 0 van de eindvlakken; trek als traction op beide eindvlakken
# (zelfevenwicht; DX = 0 alleen in het midden van x = 0).
# ------------------------------------------------------------------
def build_mp03():
    m = Model("V30", attrs={"maxmode": 2})
    m.elset("Balk", "Staal")
    mesh = hexmesh(m, "Balk", 50, 2, 4,
                   box(0, L1, -A1 / 2, A1 / 2, -B1 / 2, B1 / 2), order=2)
    m.mat("Staal", STEEL, DENS)

    def nodes_at(x, line):
        out = []
        for i, p in enumerate(m.nodes):
            if abs(p[0] - x) < 1e-9 and abs(p[2]) < 1e-9 and (
                    line or abs(p[1]) < 1e-9):
                out.append(i + 1)
        return out
    m.node_selection("Lijn0", nodes_at(0.0, True))
    m.node_selection("LijnL", nodes_at(L1, True))
    m.node_selection("Mid0", nodes_at(0.0, False))
    m.node_selection("MidL", nodes_at(L1, False))
    for sel in ("Lijn0", "LijnL"):
        m.add('<displacement selection="%s"><value>0 m</value><x>0</x>'
              '<y>0</y><z>1</z></displacement>' % sel)
    for sel in ("Mid0", "MidL"):
        m.add('<displacement selection="%s"><value>0 m</value><x>0</x>'
              '<y>1</y><z>0</z></displacement>' % sel)
    m.add('<displacement selection="Mid0"><value>0 m</value><x>1</x>'
          '<y>0</y><z>0</z></displacement>')
    m.face_selection("Kop", mesh["faces"]["u1"])
    m.face_selection("Voet", mesh["faces"]["u0"])
    s = PCR1 / AREA1
    m.add('<traction selection="Kop"><x>%r Pa</x><y>0 Pa</y><z>0 Pa</z>'
          '</traction>' % s)
    m.add('<traction selection="Voet"><x>%r Pa</x><y>0 Pa</y><z>0 Pa</z>'
          '</traction>' % -s)
    return m


TESTS.append(Case(
    "MP-03", GROUP, "Scharnierende balk onder trek, hex20 50x2x4",
    "modal_3d", build_mp03,
    # Modelfout 3D t.o.v. Euler bij L/h = 33..50 (afschuiving,
    # rotatietraagheid, scharnier op een lijn): ~0,1-0,3%; tolerantie 1%.
    [_mode_check(k + 1, f, 0.01, "Euler-balk scharnierend, P = Pcr1: "
                 "f_n0 sqrt(1 + P/Pcr_n), mode %s" % plane)
     for k, (f, plane) in enumerate(modes_pp(PCR1, 2))] +
    _mecway_checks("MP-03", 2),
    note="Volume-elementen; traction op beide eindvlakken."))


# ------------------------------------------------------------------
# MP-04 Ingeklemde balk (MO-01) die om Z draait met 200 rad/s
# (centrifugal): spanningsverstijving (trek naar de inklemming toe) en
# spin softening (RIGI_ROTA) in het rotatievlak XY. Geen gesloten
# oplossing: referentie Mecway; een eigen FE-berekening (Hermite-balk,
# consistente matrices, 20 en 200 elementen) gaf 20,918 / 20,909 Hz
# (in het vlak) en 42,635 / 42,629 Hz (uit het vlak).
# ------------------------------------------------------------------
W4 = 200.0


def build_mp04():
    m = Model("V30", attrs={"maxmode": 4})
    m.elset("Balk", "Staal")
    _e, nids = linemesh(m, "Balk", 20, (0, 0, 0), (L1, 0, 0))
    m.mat("Staal", '<geometric type="RectBar" recta="%g m" rectb="%g m" />'
          % (A1, B1), STEEL, DENS)
    m.node_selection("Inklemming", [nids[0]])
    m.add('<fixed selection="Inklemming" />')
    m.add('<centrifugal><x>0 rad/s</x><y>0 rad/s</y><z>%g rad/s</z>'
          '</centrifugal>' % W4)
    return m


TESTS.append(Case(
    "MP-04", GROUP, "Draaiende ingeklemde balk, 200 rad/s om Z (centrifugal)",
    "modal_3d", build_mp04, _mecway_checks("MP-04", 4),
    note="POU_D_E; RIGI_GEOM + RIGI_ROTA (spin softening)."))


# MP-07 Dezelfde draaiende balk als volume (hex20 50 x 2 x 3, zoals
# MO-02-h20-n2), ingeklemd op x = 0.
def build_mp07():
    m = Model("V30", attrs={"maxmode": 4})
    m.elset("Balk", "Staal")
    mesh = hexmesh(m, "Balk", 50, 2, 3,
                   box(0, L1, -A1 / 2, A1 / 2, -B1 / 2, B1 / 2), order=2)
    m.mat("Staal", STEEL, DENS)
    m.face_selection("Inklemming", mesh["faces"]["u0"])
    m.add('<fixed selection="Inklemming" />')
    m.add('<centrifugal><x>0 rad/s</x><y>0 rad/s</y><z>%g rad/s</z>'
          '</centrifugal>' % W4)
    return m


TESTS.append(Case(
    "MP-07", GROUP, "Draaiende ingeklemde balk, hex20 50x2x3, 200 rad/s om Z",
    "modal_3d", build_mp07, _mecway_checks("MP-07", 4),
    note="Volume-elementen; RIGI_GEOM + RIGI_ROTA."))


# ------------------------------------------------------------------
# MP-05 Vrij opgelegde vierkante plaat a = 0,4 m, t = 4 mm (DKT 40 x 40),
# randen in het vlak vast (DX = DY = 0) en uniform verwarmd: biaxiale druk
# N = E alpha dT t / (1 - nu) = 0,5 Ncr (thermische voorspanning, zoals
# Mecway's VibratingMembrane).
# ------------------------------------------------------------------
A5, T5 = 0.4, 0.004
D5 = E * T5 ** 3 / (12 * (1 - NU ** 2))
NCR5 = 2 * math.pi ** 2 * D5 / A5 ** 2
N5 = 0.5 * NCR5
DT5 = N5 * (1 - NU) / (E * ALPHA_T * T5)
TREF5 = 293.15


def f_plate(mm, nn, n_comp):
    s = (mm ** 2 + nn ** 2) / A5 ** 2
    w2 = (math.pi ** 4 * D5 * s ** 2 - n_comp * math.pi ** 2 * s) / (RHO * T5)
    return math.sqrt(w2) / (2 * math.pi)


MODES5 = sorted((f_plate(mm, nn, N5), "(%d,%d)" % (mm, nn))
                for mm in (1, 2, 3) for nn in (1, 2, 3))[:3]


def build_mp05():
    m = Model("V30", attrs={"maxmode": 3})
    m.elset("Plaat", "Plaat4")
    # 40 x 40: Mecway's quad4 wijkt bij 20 x 20 al zonder last +0,25% /
    # +0,93% af (mode 1 / 2; met de voorspanning +0,39% / +1,08%), DKT
    # -0,07% / -0,20%; bij 40 x 40 beide < 0,25%.
    mesh = quadmesh(m, "Plaat", 40, 40, rect(0, A5, 0, A5))
    m.mat("Plaat4", '<geometric type="Plate" thickness="%g m" />' % T5,
          STEEL_T, DENS)
    rand = set()
    for s in ("u0", "u1", "v0", "v1"):
        rand |= mesh["nodes"][s]
    m.node_selection("Rand", sorted(rand))
    for d in ("<x>1</x><y>0</y><z>0</z>", "<x>0</x><y>1</y><z>0</z>",
              "<x>0</x><y>0</y><z>1</z>"):
        m.add('<displacement selection="Rand"><value>0 m</value>%s'
              '</displacement>' % d)
    m.node_selection("RandX", sorted(mesh["nodes"]["u0"] |
                                     mesh["nodes"]["u1"]))
    m.node_selection("RandY", sorted(mesh["nodes"]["v0"] |
                                     mesh["nodes"]["v1"]))
    m.add('<noderotation selection="RandX"><value>0 rad</value>'
          '<x>1</x><y>0</y><z>0</z></noderotation>')
    m.add('<noderotation selection="RandY"><value>0 rad</value>'
          '<x>0</x><y>1</y><z>0</z></noderotation>')
    m.element_selection("Alles", sorted(mesh["elems"].values()))
    m.add('<thermalstress><referencetemperature>%g K</referencetemperature>'
          '</thermalstress>' % TREF5)
    m.add('<temperature2 selection="Alles"><value>%.12g K</value>'
          '<starttime>0 s</starttime><endtime>0 s</endtime></temperature2>'
          % (TREF5 + DT5))
    return m


TESTS.append(Case(
    "MP-05", GROUP, "Vrij opgelegde plaat, thermische biaxiale druk 0,5 Ncr "
    "(DKT 40x40)", "modal_3d", build_mp05,
    # DKT 40 x 40: discretisatie < 0,1% (MO-03); tolerantie 1%.
    [_mode_check(k + 1, f, 0.01, "Kirchhoff-plaat met biaxiale druk N = "
                 "0,5 Ncr, mode %s (Navier; Leissa 4.1)" % mn)
     for k, (f, mn) in enumerate(MODES5)] + _mecway_checks("MP-05", 3),
    note="Thermische rek (thermal stress) als voorspanning."))
