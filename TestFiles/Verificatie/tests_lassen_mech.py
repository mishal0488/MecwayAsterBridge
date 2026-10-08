# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Lassen (stap S4b): keten thermisch -> mechanisch (Mecway T31 met een
Thermal stress -> ther_3d: THER_LINEAIRE/THER_NON_LINE, daarna
STAT_NON_LINE met AFFE_VARC TEMP; bridge_thermomech.py).

Staaf 100 x 10 x 10 mm (hex8 4 x 1 x 1), staal E = 210 GPa, nu = 0,3,
alpha = 1,2e-5 1/K. De temperatuur komt uit de thermische solve: een
opgelegde temperatuur op alle knopen met een tijdstabel (uniform,
T0 -> T0 + dT -> T0). Opleggingen: symmetrievlakken x = 0, y = 0, z = 0
(vrij uitzetten, WE-07) en daarbij x = L vast in x (eenassig verhinderd,
WE-08/WE-09: dwarsrichting vrij).

Referenties (handberekening, eenassig):
  vrij:        eps = alpha dT (Timoshenko & Goodier 1970, par. 150),
               sigma = 0;
  verhinderd:  eps_mech = -eps_th; elastisch sigma = -E eps_th tot
               |sigma| = SY; daarboven bilineair (plastische modulus H =
               Mecway's tangent modulus): p = (eps_th - SY/E) / (1 + H/E),
               |sigma| = SY + H p; terugkoelen elastisch (geen omgekeerd
               vloeien als |sigma_res| < SY + H p): sigma_res = E p.
  E(T), SY(T), alpha(T) (WE-09, #ASTER-MATERIAAL, H = 0): elastisch
  sigma = -E(T) alpha(T) (T - T0) (alpha = secans t.o.v. T0); plastisch
  |sigma| = SY(T), p = alpha(T) dT - SY(T)/E(T); afkoelen elastisch:
  sigma(T) = E(T) (p - alpha(T) (T - T0)), restspanning E(T0) p.
"""
import math
import os
from verif_lib import Case, Check, Model, hexmesh, box

GROUP = "Lassen"

LB, BB = 0.1, 0.01
E, NU, ALPHA = 210e9, 0.3, 1.2e-5
SY, HM = 250e6, 20e9               # vloeigrens, plastische modulus (Mecway)
T0 = 293.15
K, RHO, CP = 50.0, 7800.0, 500.0


def _bar(m, nx=4, mat_lines=()):
    m.elset("Staaf", "Staal")
    g = hexmesh(m, "Staaf", nx, 1, 1, box(0, LB, 0, BB, 0, BB))
    m.mat("Staal",
          '<mechanical type="Isotropic" youngsmodulus="%g Pa" '
          'poissonratio="%g" thermalexpansion="%g 1/K" />' % (E, NU, ALPHA),
          '<density type="Simple" density="%g kg/m³" />' % RHO,
          '<thermal type="Isotropic" thermalconductivity="%g W/m/K" '
          'specificheat="%g J/kg/K" />' % (K, CP), *mat_lines)
    return g


def _supports(m, g, clamped):
    m.face_selection("X0", g["faces"]["u0"])
    m.face_selection("Y0", g["faces"]["v0"])
    m.face_selection("Z0", g["faces"]["w0"])
    m.node_selection("XL", sorted(g["nodes"]["u1"]))
    faces = [("X0", 1, 0, 0), ("Y0", 0, 1, 0), ("Z0", 0, 0, 1)]
    if clamped:
        m.face_selection("XLv", g["faces"]["u1"])
        faces.append(("XLv", 1, 0, 0))
    for sel, x, y, z in faces:
        m.add('<displacement selection="%s"><value>0 m</value><x>%d</x>'
              '<y>%d</y><z>%d</z></displacement>' % (sel, x, y, z))


def build_bar(history, clamped, dt=1.0, mat_lines=(), custom=None):
    """history: [(t, T in K)] van de opgelegde temperatuur op alle
    knopen."""
    def build():
        m = Model("T31", timesteps={"timeperiod": "%g s" % history[-1][0],
                                    "stepsize": "%g s" % dt})
        g = _bar(m, mat_lines=mat_lines)
        _supports(m, g, clamped)
        m.node_selection("Alle", list(range(1, len(m.nodes) + 1)))
        m.add('<temperature2 selection="Alle">%s<starttime>0 s</starttime>'
              '<endtime>0 s</endtime></temperature2>' % "".join(
                  '<value time="%g s">%.12g K</value>' % tv
                  for tv in history))
        m.add(u'<initialconditions>\n    <temperature2 selection="Alle">'
              u'<value>%.12g K</value><starttime>0 s</starttime><endtime>'
              u'0 s</endtime></temperature2>\n  </initialconditions>' % T0)
        m.add('<thermalstress><referencetemperature>%g K'
              '</referencetemperature></thermalstress>' % T0)
        if custom:
            m.add(u"<ccxmodeldefinition>%s</ccxmodeldefinition>" % "\n".join(
                u"** " + line for line in custom))
        return m
    return build


# ------------------------------------------------------------------
# WE-07: vrij uitzetten, uniforme opwarming en terugkoelen
# ------------------------------------------------------------------
DT7 = 100.0
HIST7 = [(0.0, T0), (5.0, T0 + DT7), (10.0, T0)]


def _max_abs(var):
    return lambda r, t=None: max(abs(v) for v in r.step(t)[var].values())


WE07 = Case(
    "WE-07", GROUP, "Keten thermisch-mechanisch: vrij uitzetten",
    "ther_3d", build_bar(HIST7, clamped=False),
    [Check("u_x eind t=5 s", lambda r: r.at("displx", (LB, 0, 0), time=5.0),
           ref=ALPHA * DT7 * LB, tol=1e-3, unit="m",
           source="alpha dT L (Timoshenko & Goodier 1970, par. 150)"),
     Check("u_y zijkant t=5 s",
           lambda r: r.at("disply", (LB, BB, 0), time=5.0),
           ref=ALPHA * DT7 * BB, tol=1e-3, unit="m", source="alpha dT b"),
     Check("u_x eind t=2 s", lambda r: r.at("displx", (LB, 0, 0), time=2.0),
           ref=ALPHA * DT7 * 0.4 * LB, tol=1e-3, unit="m",
           source="alpha dT(t) L, dT lineair in t (tijdstabel)"),
     Check("max |sigma_xx| t=5 s",
           lambda r: max(abs(v) for v in r.step(5.0)["stressxx"].values()),
           ref=0.0, tol_abs=1e-4 * E * ALPHA * DT7, unit="Pa",
           source="vrij uitzetten: spanningsloos (tolerantie 1e-4 E alpha "
                  "dT)"),
     Check("u_x eind na terugkoelen",
           lambda r: r.at("displx", (LB, 0, 0), time=10.0),
           ref=0.0, tol_abs=1e-4 * ALPHA * DT7 * LB, unit="m",
           source="terug op T0: geen rek"),
     Check("temperatuur t=5 s", lambda r: r.at("temperature", (LB / 2, 0, 0),
                                               time=5.0),
           ref=T0 + DT7, tol=1e-6, unit="K",
           source="opgelegde temperatuur (thermisch resultaat res_ther)")],
    note="T31 met Thermal stress (referentie T0); opgelegde temperatuur "
         "op alle knopen T0 -> T0 + 100 K (5 s) -> T0 (10 s), dt = 1 s; "
         "drie symmetrievlakken (statisch bepaald)")


# ------------------------------------------------------------------
# WE-08: eenassig verhinderd, bilineair isotroop (Mecway
# BilinearIsotropic, tangent modulus = plastische modulus H)
# ------------------------------------------------------------------
DT8 = 150.0
HIST8 = [(0.0, T0), (15.0, T0 + DT8), (30.0, T0)]
PLASTIC8 = ('<plastic type="BilinearIsotropic" yieldstrength="%g Pa" '
            'tangentmodulus="%g Pa" />' % (SY, HM))


def restrained(dt_k):
    """(sigma_xx, p) eenassig verhinderd bij opwarming dt_k (monotoon,
    bilineair isotroop met plastische modulus HM)."""
    eps = ALPHA * dt_k
    if E * eps <= SY:
        return -E * eps, 0.0
    p = (eps - SY / E) / (1.0 + HM / E)
    return -(SY + HM * p), p


SIG8, P8 = restrained(DT8)
SIGRES8 = E * P8                       # terugkoelen elastisch
assert SIGRES8 < SY + HM * P8          # geen omgekeerd vloeien
# Dwarsrek: eps_yy = alpha dT - nu sigma_xx / E + p / 2 (plastisch
# volumebehoudend: eps_p,yy = -eps_p,xx / 2)
EPSY8 = ALPHA * DT8 - NU * SIG8 / E + P8 / 2.0
EPSY8_RES = -NU * SIGRES8 / E + P8 / 2.0
MID = (LB / 2, BB / 2, BB)

WE08 = Case(
    "WE-08", GROUP, "Keten: verhinderde staaf, bilineair plastisch, "
    "restspanning", "ther_3d",
    build_bar(HIST8, clamped=True, mat_lines=(PLASTIC8,)),
    [Check("sigma_xx elastisch dT=50 K",
           lambda r: r.at("stressxx", MID, time=5.0),
           ref=-E * ALPHA * 50.0, tol=0.01, unit="Pa",
           source="-E alpha dT (verhinderde uitzetting, dwars vrij)"),
     Check("sigma_xx plastisch dT=150 K",
           lambda r: r.at("stressxx", MID, time=15.0),
           ref=SIG8, tol=0.01, unit="Pa",
           source="bilineair: p = (alpha dT - SY/E)/(1 + H/E), |sigma| = "
                  "SY + H p (H = Mecway tangent modulus)"),
     Check("eq. plastische rek dT=150 K",
           lambda r: r.at("equivalentplasticstrain", MID, time=15.0),
           ref=P8, tol=0.01, unit="-", source="p, idem"),
     Check("u_y zijkant dT=150 K",
           lambda r: r.at("disply", (LB / 2, BB, 0), time=15.0),
           ref=EPSY8 * BB, tol=0.01, unit="m",
           source="eps_yy = alpha dT - nu sigma/E + p/2"),
     Check("restspanning sigma_xx na terugkoelen",
           lambda r: r.at("stressxx", MID, time=30.0),
           ref=SIGRES8, tol=0.01, unit="Pa",
           source="elastisch ontlasten: sigma_res = E p (geen omgekeerd "
                  "vloeien: E p < SY + H p)"),
     Check("plastische rek xx na terugkoelen",
           lambda r: r.at("plasticstrainxx", MID, time=30.0),
           ref=-P8, tol=0.01, unit="-",
           source="eps_p,xx = -p (blijft na ontlasten)"),
     Check("u_y zijkant na terugkoelen",
           lambda r: r.at("disply", (LB / 2, BB, 0), time=30.0),
           ref=EPSY8_RES * BB, tol=0.01, unit="m",
           source="eps_yy = -nu sigma_res/E + p/2"),
     Check("|reactie x| op x = L na terugkoelen",
           lambda r: abs(r.total("reactionforcex", "XL", time=30.0)),
           ref=SIGRES8 * BB * BB, tol=0.01, unit="N",
           source="sigma_res A")],
    note="x = 0 en x = L vast in x, y = 0 / z = 0 symmetrie; T0 -> T0 + "
         "150 K (15 s) -> T0 (30 s), dt = 1 s; SY = 250 MPa, H = 20 GPa; "
         "vloeien vanaf dT = 99,2 K")


# ------------------------------------------------------------------
# WE-09: E(T), SY(T), alpha(T) via #ASTER-MATERIAAL (ideaal plastisch)
# ------------------------------------------------------------------
TC9 = (20.0, 520.0)
E9 = (210e9, 150e9)
SY9 = (300e6, 100e6)
AL9 = (1.2e-5, 1.4e-5)            # secans t.o.v. 20 graden C (= T0)
DT9 = 400.0
HIST9 = [(0.0, T0), (40.0, T0 + DT9), (80.0, T0)]
MAT9 = ["ASTER-MATERIAAL Staal: TC=(%g, %g), E=(%g, %g)," % (TC9 + E9),
        "ASTER-MATERIAAL     SY=(%g, %g), H=0, ALPHA=(%g, %g)"
        % (SY9 + AL9)]


def _lin(pair, dt_k):
    th = dt_k / (TC9[1] - TC9[0])
    return pair[0] + (pair[1] - pair[0]) * th


def tdep_path(dt_peak=DT9, n=4000):
    """Handberekening WE-09, eenassig verhinderd, ideaal plastisch met
    E(T), SY(T), alpha(T) lineair in T (secans t.o.v. T0): opwarmen tot
    dt_peak, dan afkoelen tot T0. Returns dict met sigma op de
    controlepunten."""
    out = {}
    # opwarmen: elastisch tot E alpha dT = SY, daarna sigma = -SY,
    # p = alpha dT - SY/E (monotoon)
    def heat(dt_k):
        e, a, sy = _lin(E9, dt_k), _lin(AL9, dt_k), _lin(SY9, dt_k)
        if e * a * dt_k <= sy:
            return -e * a * dt_k, 0.0
        return -sy, a * dt_k - sy / e
    out["heat50"] = heat(50.0)[0]
    out["peak"], p = heat(dt_peak)
    # afkoelen in kleine stappen: elastisch sigma = E (p - eps_th), tot
    # sigma = +SY (omgekeerd vloeien), daarna sigma = SY(T)
    for k in range(1, n + 1):
        dt_k = dt_peak * (1.0 - k / float(n))
        e, a, sy = _lin(E9, dt_k), _lin(AL9, dt_k), _lin(SY9, dt_k)
        sig = e * (p - a * dt_k)
        if sig > sy:
            p = a * dt_k + sy / e           # p daalt (omgekeerd vloeien)
            sig = sy
        if abs(dt_k - 300.0) < 0.5 * dt_peak / n:
            out["cool300"] = sig
        out["res"] = sig
    out["p_res"] = p
    return out


REF9 = tdep_path()

WE09 = Case(
    "WE-09", GROUP, "Keten: E(T), SY(T), alpha(T) (#ASTER-MATERIAAL)",
    "ther_3d",
    build_bar(HIST9, clamped=True, custom=MAT9),
    [Check("sigma_xx elastisch dT=50 K",
           lambda r: r.at("stressxx", MID, time=5.0),
           ref=REF9["heat50"], tol=0.01, unit="Pa",
           source="-E(T) alpha(T) dT, alpha secans t.o.v. T0 (U4.43.01 "
                  "ELAS_FO TEMP_DEF_ALPHA)"),
     Check("sigma_xx plastisch dT=400 K",
           lambda r: r.at("stressxx", MID, time=40.0),
           ref=REF9["peak"], tol=0.01, unit="Pa", source="-SY(T)"),
     Check("sigma_xx afkoelen dT=300 K (elastisch)",
           lambda r: r.at("stressxx", MID, time=50.0),
           ref=REF9["cool300"], tol=0.01, unit="Pa",
           source="E(T) (p - alpha(T) dT), p van de piek"),
     Check("restspanning sigma_xx na terugkoelen",
           lambda r: r.at("stressxx", MID, time=80.0),
           ref=REF9["res"], tol=0.01, unit="Pa",
           source="omgekeerd vloeien tijdens afkoelen: sigma_res = SY(T0)"),
     Check("plastische rek xx na terugkoelen",
           lambda r: r.at("plasticstrainxx", MID, time=80.0),
           ref=-REF9["p_res"], tol=0.01, unit="-",
           source="eps_p = -(sigma_res / E(T0)) (totale mechanische rek 0)")],
    note="E 210 -> 150 GPa, SY 300 -> 100 MPa, alpha 1,2 -> 1,4e-5 1/K "
         "lineair tussen 20 en 520 graden C; H = 0; T0 -> T0 + 400 K (40 "
         "s) -> T0 (80 s), dt = 1 s")

TESTS = [WE07, WE08, WE09]


# ------------------------------------------------------------------
# WE-10: lasrups op een plaat (#ASTER-LASBRON + E(T), SY(T), alpha(T)):
# restspanningspatroon en vervorming, kwalitatief (Radaj 1992, Masubuchi
# 1980) plus evenwicht van de restspanning (vrije plaat)
# ------------------------------------------------------------------
from tests_lassen import _axis, _grid

LX10, BY10, TZ10 = 0.040, 0.015, 0.005       # halve plaat, y = 0 symmetrie
X10S = _axis([(0.0, LX10, 40)])
Y10S = _axis([(0.0, 0.006, 6), (0.006, BY10, 4)])
Z10S = _axis([(-TZ10, 0.0, 5)])
V10, Q10 = 0.005, 1000.0                     # m/s, W (eta Q, hele bron)
SRC10 = ["ASTER-LASBRON Q=%g, eta=1, v=%g, start=(0.005, 0, 0), "
         "richting=(1, 0, 0), lengte=0.03" % (Q10, V10),
         "ASTER-LASBRON af=0.002, ar=0.004, b=0.002, c=0.002, symmetrie=1"]
# Staal, temperatuurafhankelijk (illustratief, verloop zoals in de
# literatuur voor constructiestaal, bijv. S355; geen meetwaarden)
TC10 = (20, 200, 400, 600, 800, 1000, 1200, 1500)
E10 = (210e9, 200e9, 180e9, 120e9, 60e9, 30e9, 15e9, 10e9)
SY10 = (355e6, 300e6, 250e6, 150e6, 50e6, 25e6, 15e6, 10e6)
AL10 = (1.2e-5, 1.25e-5, 1.3e-5, 1.35e-5, 1.4e-5, 1.45e-5, 1.5e-5, 1.5e-5)
MAT10 = ["ASTER-MATERIAAL Staal: TC=(%s)," % ", ".join("%g" % t for t in TC10),
         "ASTER-MATERIAAL     E=(%s)," % ", ".join("%g" % v for v in E10),
         "ASTER-MATERIAAL     SY=(%s)," % ", ".join("%g" % v for v in SY10),
         "ASTER-MATERIAAL     ALPHA=(%s), H=1e9, TC_CUT=1500"
         % ", ".join("%g" % v for v in AL10)]
TEND10, DT10 = 60.0, 0.3


def build_we10():
    m = Model("T31", timesteps={"timeperiod": "%g s" % TEND10,
                                "stepsize": "%g s" % DT10})
    m.elset("Plaat", "Staal")
    g = hexmesh(m, "Plaat", len(X10S) - 1, len(Y10S) - 1, len(Z10S) - 1,
                _grid(X10S, Y10S, Z10S))
    m.mat("Staal",
          '<mechanical type="Isotropic" youngsmodulus="%g Pa" '
          'poissonratio="%g" thermalexpansion="%g 1/K" />' % (E, NU, ALPHA),
          '<density type="Simple" density="%g kg/m³" />' % 7800.0,
          '<thermal type="Isotropic" thermalconductivity="%g W/m/K" '
          'specificheat="%g J/kg/K" />' % (30.0, 600.0))
    m.node_selection("Alle", list(range(1, len(m.nodes) + 1)))
    m.add(u'<initialconditions>\n    <temperature2 selection="Alle">'
          u'<value>%.12g K</value><starttime>0 s</starttime><endtime>'
          u'0 s</endtime></temperature2>\n  </initialconditions>' % T0)
    # Warmteafvoer: de rand y = 15 mm op T0 (de rest van een bredere,
    # koude plaat); de andere vlakken geisoleerd.
    m.face_selection("Rand", g["faces"]["v1"])
    m.add('<temperature2 selection="Rand"><value>%.12g K</value><starttime>'
          '0 s</starttime><endtime>0 s</endtime></temperature2>' % T0)
    # Mechanisch: symmetrievlak y = 0, verder statisch bepaald (3-2-1)
    m.face_selection("Sym", g["faces"]["v0"])
    m.add('<displacement selection="Sym"><value>0 m</value><x>0</x><y>1</y>'
          '<z>0</z></displacement>')
    a = m.node(0.0, 0.0, -TZ10)
    b = m.node(LX10, 0.0, -TZ10)
    m.node_selection("A", [a])
    m.node_selection("B", [b])
    # Mecway's displacement is een verplaatsing langs een richting (x, y,
    # z): per vastgehouden component een eigen voorwaarde.
    m.add('<displacement selection="A"><value>0 m</value><x>1</x><y>0</y>'
          '<z>0</z></displacement>')
    m.add('<displacement selection="A"><value>0 m</value><x>0</x><y>0</y>'
          '<z>1</z></displacement>')
    m.add('<displacement selection="B"><value>0 m</value><x>0</x><y>0</y>'
          '<z>1</z></displacement>')
    m.add('<thermalstress><referencetemperature>%g K</referencetemperature>'
          '</thermalstress>' % T0)
    m.add(u"<ccxmodeldefinition>%s</ccxmodeldefinition>" % "\n".join(
        u"** " + line for line in SRC10 + MAT10))
    return m


def _section(r, xc=LX10 / 2, time=None, moment=False):
    """int sigma_xx dA (of int sigma_xx (z - zc) dA, zc = -t/2) over de
    doorsnede x = xc van de halve plaat (trapeziumregel op de knopen)."""
    d = r.step(time)["stressxx"]
    zc = -TZ10 / 2

    def w(arr, i):
        lo = arr[i] - arr[max(i - 1, 0)]
        hi = arr[min(i + 1, len(arr) - 1)] - arr[i]
        return 0.5 * (lo + hi)
    total = 0.0
    for j in range(len(Y10S)):
        for k in range(len(Z10S)):
            s = d[r.nearest((xc, Y10S[j], Z10S[k]))]
            total += s * w(Y10S, j) * w(Z10S, k) * (
                (Z10S[k] - zc) if moment else 1.0)
    return total


def _flag(cond):
    return 1.0 if cond else 0.0


MID10 = LX10 / 2
HM10 = 1e9
# Bovengrens dwarskrimp (hele plaat): de warmte per lengte q/v als uniforme
# opwarming van een strook over de dikte t die volledig krimpt:
# alpha q / (v rho c t) (afgeleid: int rho c t dT dy = q/v).
SHRINK10 = ALPHA * Q10 / (V10 * 7800.0 * 600.0 * TZ10)


def _we10_checks(lq):
    out = [
        Check("evenwicht: int sigma_xx dA (x = 20 mm) / (SY A)",
              lambda r: _section(r) / (SY10[0] * BY10 * TZ10),
              ref=0.0, tol_abs=0.03, unit="-", in_series=False,
              source="vrije plaat: de restspanning is in evenwicht (netto "
                     "langskracht 0); tolerantie 3% van SY A"),
        Check("evenwicht: int sigma_xx (z - zc) dA / (SY A t)",
              lambda r: _section(r, moment=True) /
              (SY10[0] * BY10 * TZ10 * TZ10),
              ref=0.0, tol_abs=0.03, unit="-", in_series=False,
              source="vrije plaat: netto buigend moment 0"),
        Check("von Mises lasas boven / (SY(T0) + H p): vloeit",
              lambda r: r.at("vonmises", (MID10, 0.0, 0.0)) /
              (SY10[0] + HM10 * r.at("equivalentplasticstrain",
                                     (MID10, 0.0, 0.0))),
              ref=1.0, tol=0.05, unit="-", in_series=False,
              source="in de las vloeit het materiaal bij het afkoelen: von "
                     "Mises = actuele vloeigrens SY + H p (knoopwaarde; 5%)"),
        Check("sigma_xx lasas boven (x = 20 mm) > 0 (trek)",
              lambda r: _flag(r.at("stressxx", (MID10, 0.0, 0.0)) > 0.0),
              ref=1.0, tol_abs=0.0, unit="-", in_series=False,
              source="Radaj (1992), Masubuchi (1980): langsspanning in de las "
                     "trek"),
        Check("sigma_xx op y >= 8 mm (x = 20 mm, midden dikte) < 0 (druk)",
              lambda r: _flag(all(
                  r.at("stressxx", (MID10, y, -TZ10 / 2)) < 0.0
                  for y in Y10S if y >= 0.008)),
              ref=1.0, tol_abs=0.0, unit="-", in_series=False,
              source="naast de las druk (evenwicht met de trek in de las)"),
        Check("hoekverdraaiing: dz(y = 15) - dz(y = 0) onder > 0",
              lambda r: _flag(r.at("displz", (MID10, BY10, -TZ10)) -
                              r.at("displz", (MID10, 0.0, -TZ10)) > 0.0),
              ref=1.0, tol_abs=0.0, unit="-", in_series=False,
              source="Masubuchi (1980): bij een rups op een vlak "
                     "(gedeeltelijke inbranding) draaien de randen naar de "
                     "laskant"),
        Check("langskrimp: dx(x = 40) - dx(x = 0) op de lasas < 0",
              lambda r: _flag(r.at("displx", (LX10, 0.0, 0.0)) -
                              r.at("displx", (0.0, 0.0, 0.0)) < 0.0),
              ref=1.0, tol_abs=0.0, unit="-", in_series=False,
              source="krimp in de lasrichting (Radaj 1992)"),
        Check("dwarskrimp hele plaat / (alpha q/(v rho c t))",
              lambda r: -2.0 * r.at("disply", (MID10, BY10, -TZ10 / 2)) /
              SHRINK10,
              ref=0.6, tol_abs=0.4, unit="-", in_series=False,
              source="krimp dwars op de las (Radaj 1992), tussen 0,2 en 1 x "
                     "de bovengrens alpha q/(v rho c t) = %.3g mm (volledig "
                     "krimpende strook; de rups is niet door de dikte)"
                     % (SHRINK10 * 1e3)),
    ]
    if lq:
        out.append(Check(
            "sigma_xx lasas boven (x = 20 mm) / SY(T0)",
            lambda r: r.at("stressxx", (MID10, 0.0, 0.0)) / SY10[0],
            ref=1.0, tol=0.35, unit="-", in_series=False,
            source="Radaj (1992), Masubuchi (1980): langsspanning in de las "
                   "ongeveer de vloeigrens (kwalitatief, +-35%)"))
    return out


WE10 = Case(
    "WE-10", GROUP, "Keten: lasrups op een plaat, restspanning en vervorming",
    "ther_3d", build_we10, _we10_checks(False),
    timeout=3600,
    note="halve plaat 40 x 15 x 5 mm (hex8 1 mm langs de las), Goldak "
         "2/4/2/2 mm, 1000 W, 5 mm/s over 30 mm (6 s), afkoelen tot 60 s "
         "(rand y = 15 mm op T0), dt = 0,3 s; E/SY/alpha(T) via "
         "#ASTER-MATERIAAL, H = 1 GPa, TC_CUT = 1500; symmetrie + statisch "
         "bepaald; hex8 3D (locking, docs/Lassen.md par. 9.5)")

TESTS.append(WE10)
# WE-10b: dezelfde plaat mechanisch kwadratisch (LINE_QUAD); ~1 uur, daarom
# alleen met de omgevingsvariabele VERIF_ZWAAR=1 (docs/Verificatie.md par. 8).
WE10B = Case(
    "WE-10b", GROUP, WE10.title + " (LINE_QUAD)", "ther_3d", build_we10,
    _we10_checks(True), general={"elem_hex8": "LINE_QUAD"},
    settings={"time_limit": 18000}, timeout=18000,
    note=WE10.note.replace("hex8 3D (locking, docs/Lassen.md par. 9.5)",
                           "mechanisch kwadratisch (LINE_QUAD, temperatuur "
                           "via PROJ_CHAMP)"))
if os.environ.get("VERIF_ZWAAR"):
    TESTS.append(WE10B)
