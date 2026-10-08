# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Lassen (stap S4c): toevoegmateriaal dat geboren wordt (#ASTER-LASBRON
geboorte=..., "quiet element"-methode; bridge_welding._births,
bridge_welding.birth_materials / ther_loop, bridge_thermomech._birth_solve).

Regel van de bridge (docs/Lassen.md par. 10): een element van het
toevoegmateriaal wordt geboren op tb = max(t0, t0 + (s_min - voorloop)/v),
s_min = de kleinste booglengte van zijn hoekknopen langs de baan (rechte
baan: x - x_start, voor de start 0, voorbij het einde doorgetrokken); het
leeft vanaf de eerste staptijd na tb (tb = t0: vanaf t0) en wordt nooit
geboren als tb na het einde van de bron ligt. Voor de geboorte: rho*cp x 1e-4, lambda x 1e-8, E x 1e-4; de bron
werkt er niet. Mechanisch spanningsvrij in de toestand van het begin van de
geboortestap. De tests rekenen de geboortetijden zelf (onafhankelijk van de
bridge) uit deze regel.

Referenties:
  energiebalans (geisoleerd): int_levend rho (H(T) - H(T0)) dV = eta Q t;
  eenassige staaf (Timoshenko & Goodier 1970, par. 150): vrij uitzetten
  u = alpha dT x; verhinderd na de geboorte op T1 en afkoelen naar T0:
  sigma = E alpha (T1 - T0) (alpha dT t.o.v. de geboortetemperatuur),
  bilineair: p = (alpha dT - SY/E)/(1 + H/E), sigma = SY + H p;
  restspanningspatroon (kwalitatief): Radaj (1992), Masubuchi (1980).
"""
import math
from verif_lib import Case, Check, Model, hexmesh, box
from tests_lassen import _axis, _grid, _enthalpy_fn
import tests_lassen_mech as tm

GROUP = "Lassen"

T0 = 293.15
RHO, CP, K = 7800.0, 600.0, 30.0
QW = 500.0                      # eta Q van de hele bron [W]; halve plaat: 250 W
V = 0.005                       # m/s
AF, AR, BB_, CC = 0.0015, 0.003, 0.0015, 0.0015
LEAD = 2.0 * AF                 # standaard voorloop (2 x af)
THETA, DT = 0.57, 0.1


def _esc(line):
    return line.replace("<", "&lt;").replace(">", "&gt;")


def _custom(m, lines):
    m.add(u"<ccxmodeldefinition>%s</ccxmodeldefinition>" % "\n".join(
        u"** " + _esc(line) for line in lines))


def _initial(m):
    m.node_selection("Alle", list(range(1, len(m.nodes) + 1)))
    m.add(u'<initialconditions>\n    <temperature2 selection="Alle"><value>'
          u'%.12g K</value><starttime>0 s</starttime><endtime>0 s</endtime>'
          u'</temperature2>\n  </initialconditions>' % T0)


def _thermal_mat(tempdep=None):
    lines = ['<density type="Simple" density="%g kg/m³" />' % RHO]
    if tempdep:
        rows = "".join('<thermalconductivity temp="%g K" value="%g W/m/K" />'
                       % tv for tv in tempdep[0])
        rows += "".join('<specificheat temp="%g K" value="%g J/kg/K" />'
                        % tv for tv in tempdep[1])
        lines.append('<thermal type="TempDependent" thermalconductivity='
                     '"%g W/m/K" specificheat="%g J/kg/K">%s</thermal>'
                     % (K, CP, rows))
    else:
        lines.append('<thermal type="Isotropic" thermalconductivity='
                     '"%g W/m/K" specificheat="%g J/kg/K" />' % (K, CP))
    return lines


# ------------------------------------------------------------------
# Geboorte en energie (onafhankelijk van de bridge)
# ------------------------------------------------------------------
class Pass(object):
    """Rechte baan langs x: start x0 (y = 0), lengte L, starttijd t0."""

    def __init__(self, elset, x0, length, t0, v=V, lead=LEAD):
        self.elset, self.x0, self.length, self.t0 = elset, x0, length, t0
        self.v, self.lead = v, lead
        self.t_end = t0 + length / v

    def t_born(self, pts):
        # voor de start begrensd op 0, voorbij het einde doorgetrokken (daar
        # komt de bron niet: tb > t_eind = nooit geboren)
        s = min(max(p[0] - self.x0, 0.0) for p in pts)
        return max(self.t0, self.t0 + (s - self.lead) / self.v)


def _alive(passes, name, pts, t):
    """Levend op staptijd t: t > tb, of (geboren bij de start, tb = t0)
    t >= t0."""
    for p in passes:
        if p.elset == name:
            tb = p.t_born(pts)
            if tb > p.t_end:
                return False
            return t > tb + 1e-9 or (tb <= p.t0 and t >= p.t0 - 1e-9)
    return True


def _q_theta(passes, t, q, tend, dt=DT, theta=THETA):
    """Ingebrachte energie tot staptijd t volgens het theta-schema: sum dt
    (theta F(t_k) + (1 - theta) F(t_k-1)), F = q als de bron op t_k aan
    staat (bridge_welding GOLDAK_FUNC: aan op [t0, t_eind], ook op t_eind
    zelf). Een bron die eerder stopt dan de analyse krijgt in de stap erna
    nog (1 - theta) F(t_eind): na het einde q (t_aan + (1 - theta) dt)."""
    n_all = int(round(tend / dt))
    total = 0.0
    for k in range(1, n_all + 1):
        if tend * k / float(n_all) > t + 1e-9:
            break
        for p in passes:
            length = min(p.length, p.v * (p.t_end - p.t0))

            def on(tk):
                s = p.v * (tk - p.t0)
                return 0.0 if (s < 0.0 or s > length * (1.0 + 1e-9)) \
                    else 1.0
            total += dt * q * (theta * on(tend * k / float(n_all)) +
                               (1.0 - theta) * on(tend * (k - 1) /
                                                  float(n_all)))
    return total


_GP = (0.5 - 0.5 / math.sqrt(3.0), 0.5 + 0.5 / math.sqrt(3.0))


def _energy(r, t, passes, alive=True, cp_table=None):
    """int rho (H(T) - H(T0)) dV over de levende (alive=True) of de nog niet
    geboren (False) elementen, met de ECHTE rho*cp (ook voor de ongeboren:
    'zoveel zou erin zitten als het echt materiaal was'). Hex8-blokken:
    constante cp exact met het gemiddelde van de hoekknopen, c(T) met 2x2x2
    Gauss-punten (zoals tests_lassen._energy)."""
    temps = r.step(t)["temperature"]
    nodes = r.model.nodes
    h = _enthalpy_fn(cp_table)
    total = 0.0
    for name, _mat, elems in r.model.elsets:
        for _eid, _shape, nids, _a in elems:
            pts = [nodes[n - 1] for n in nids]
            if _alive(passes, name, pts, t) != alive:
                continue
            vol = abs((pts[1][0] - pts[0][0]) * (pts[3][1] - pts[0][1]) *
                      (pts[4][2] - pts[0][2]))
            tn = [temps[n] for n in nids]
            if cp_table is None:
                total += vol * CP * sum(x - T0 for x in tn) / 8.0
                continue
            for u in _GP:
                for v in _GP:
                    for w in _GP:
                        wts = [(1 - u) * (1 - v) * (1 - w),
                               u * (1 - v) * (1 - w), u * v * (1 - w),
                               (1 - u) * v * (1 - w), (1 - u) * (1 - v) * w,
                               u * (1 - v) * w, u * v * w, (1 - u) * v * w]
                        total += vol / 8.0 * h(sum(a * b for a, b in
                                                   zip(wts, tn)))
    return RHO * total


def _top_nodes(r, elset, z, x_min=-1.0):
    """Knopen van elset op hoogte z (alleen van dat deel: niet op het
    grensvlak x = x_min met een ander deel)."""
    nodes = r.model.nodes
    out = set()
    for name, _m, elems in r.model.elsets:
        if name != elset:
            continue
        for _e, _s, nids, _a in elems:
            out.update(n for n in nids if abs(nodes[n - 1][2] - z) < 1e-9 and
                       nodes[n - 1][0] > x_min + 1e-9)
    return sorted(out)


def _power(index):
    def value(r):
        import sys
        from verif_lib import SCRIPTS
        if SCRIPTS not in sys.path:
            sys.path.insert(0, SCRIPTS)
        import bridge_welding
        rows = bridge_welding.read_energy(r.mess)[index]
        return min(p for _t, p, _q in rows)
    return value


# ------------------------------------------------------------------
# WE-11: rups in twee stukken op een plaat (thermisch)
# ------------------------------------------------------------------
# Halve plaat (y >= 0) 30 x 8 x 4 mm; rups 30 x 2,25 x 1,5 mm erop in twee
# stukken (LASNAAD_1: x 0..15 mm, LASNAAD_2: x 15..30 mm), elk met een eigen
# bron: stuk 1 van t = 0 tot 1,8 s (x 3..12 mm), stuk 2 van t = 4 tot 5,8 s
# (x 18..27 mm); afkoelen tot 9 s. Alle vlakken geisoleerd.
XS = _axis([(0.0, 0.030, 40)])
YP = _axis([(0.0, 0.00225, 3), (0.00225, 0.008, 5)])
ZP = _axis([(-0.004, -0.0015, 2), (-0.0015, 0.0, 2)])
YR = _axis([(0.0, 0.00225, 3)])
ZR = _axis([(0.0, 0.0015, 2)])
PASSES11 = [Pass("LASNAAD_1", 0.003, 0.009, 0.0),
            Pass("LASNAAD_2", 0.018, 0.009, 4.0)]
TEND11 = 9.0
SRC11 = [
    "ASTER-LASBRON Q=%g, v=%g, start=(0.003, 0, 0.0015), richting=(1, 0, 0),"
    " lengte=0.009" % (QW, V),
    "ASTER-LASBRON af=%g, ar=%g, b=%g, c=%g, symmetrie=1, "
    "geboorte=<<LASNAAD_1>>" % (AF, AR, BB_, CC),
    "ASTER-LASBRON Q=%g, v=%g, start=(0.018, 0, 0.0015), richting=(1, 0, 0),"
    " lengte=0.009, t0=4" % (QW, V),
    "ASTER-LASBRON af=%g, ar=%g, b=%g, c=%g, symmetrie=1, "
    "geboorte=<<LASNAAD_2>>" % (AF, AR, BB_, CC)]
# Temperatuurafhankelijk (WE-11b): lambda 30 -> 25 W/m/K, c 500 -> 800
# J/kg/K tussen 300 en 1500 K (illustratief)
TEMPDEP11 = ([(300.0, 30.0), (1500.0, 25.0)], [(300.0, 500.0), (1500.0, 800.0)])


def build_we11(tempdep=None):
    def build():
        m = Model("T31", timesteps={"timeperiod": "%g s" % TEND11,
                                    "stepsize": "%g s" % DT})
        m.elset("Plaat", "Staal")
        m.elset("LASNAAD_1", "Staal")
        m.elset("LASNAAD_2", "Staal")
        hexmesh(m, "Plaat", len(XS) - 1, len(YP) - 1, len(ZP) - 1,
                _grid(XS, YP, ZP))
        half = len(XS) // 2
        hexmesh(m, "LASNAAD_1", half, len(YR) - 1, len(ZR) - 1,
                _grid(XS[:half + 1], YR, ZR))
        hexmesh(m, "LASNAAD_2", half, len(YR) - 1, len(ZR) - 1,
                _grid(XS[half:], YR, ZR))
        m.mat("Staal", *_thermal_mat(tempdep))
        _initial(m)
        _custom(m, SRC11)
        return m
    return build


def _we11_checks(cp_table=None):
    q = 0.5 * QW                                   # halve plaat
    src = ("energiebalans, geisoleerd: int_levend rho (H(T) - H(T0)) dV = "
           "ingebrachte energie (halve plaat Q/2), theta-schema met de bron "
           "aan op [t0, t_eind] (_q_theta)")
    return [
        Check("energie levend t=3.5 s (stuk 2 nog niet geboren)",
              lambda r: _energy(r, 3.5, PASSES11, True, cp_table),
              ref=_q_theta(PASSES11, 3.5, q, TEND11), tol=0.02, unit="J",
              source=src + "; bron 1 klaar (1,8 s), stuk 2 ongeboren"),
        Check("energie levend t=5 s (stuk 2 wordt geboren)",
              lambda r: _energy(r, 5.0, PASSES11, True, cp_table),
              ref=_q_theta(PASSES11, 5.0, q, TEND11), tol=0.02, unit="J",
              source=src + "; tijdens bron 2"),
        Check("warmte in ongeboren stuk 2 als echt materiaal t=3.5 s / (Q t)",
              lambda r: _energy(r, 3.5, PASSES11, False, cp_table) /
              _q_theta(PASSES11, 3.5, q, TEND11),
              ref=0.0, tol_abs=0.01, unit="-", in_series=False,
              source="voor de geboorte geen warmte in het toevoegmateriaal "
                     "(werkelijk opgeslagen: x 1e-4); zelfs met de echte "
                     "rho*cp < 1% van de ingebrachte energie: alleen de "
                     "knopen op het grensvlak met stuk 1 en de plaat zijn "
                     "warm"),
        Check("max T bovenop ongeboren stuk 2 (x > 15 mm) t=3.5 s",
              lambda r: max(r.step(3.5)["temperature"][n] for n in
                            _top_nodes(r, "LASNAAD_2", ZR[-1], 0.015)),
              ref=T0, tol_abs=1.0, unit="K", in_series=False,
              source="ongeboren: geen geleiding en een diagonale "
                     "warmtecapaciteit (3D_DIAG), dus de begintemperatuur, ook "
                     "naast het warme stuk 1 (tolerantie 1 K)"),
        Check("energie levend t=9 s (beide stukken geboren)",
              lambda r: _energy(r, TEND11, PASSES11, True, cp_table),
              ref=_q_theta(PASSES11, TEND11, q, TEND11), tol=0.02, unit="J",
              source=src + "; de geboorte brengt de warmte van de knopen op "
                     "het grensvlak mee (docs/Lassen.md par. 10.4)"),
        Check("opgenomen vermogen bron 1 (job)", _power(1), ref=q, tol=0.02,
              unit="W", in_series=False,
              source="BRIDGE_LASBRON in de .mess (CALC_VECT_ELEM CHAR_THER, "
                     "met de bron alleen in geboren toevoegmateriaal)"),
        Check("opgenomen vermogen bron 2 (job)", _power(2), ref=q, tol=0.02,
              unit="W", in_series=False,
              source="idem, bron 2"),
    ]


TESTS = [
    Case("WE-11", GROUP, "Toevoegmateriaal: rups in twee stukken, "
         "energiebalans (thermisch)", "ther_3d", build_we11(),
         _we11_checks(),
         note="halve plaat 30 x 8 x 4 mm + rups 30 x 2,25 x 1,5 mm in twee "
              "stukken (LASNAAD_1/LASNAAD_2, geboorte=<<..>>), Goldak "
              "1,5/3/1,5/1,5 mm, 500 W, 5 mm/s, stuk 1 0..2,4 s, stuk 2 "
              "4..6,4 s, dt = 0,1 s; THER_LINEAIRE (THER_FO van INST en "
              "NEUT1)"),
    Case("WE-11b", GROUP, "Toevoegmateriaal: twee stukken met lambda(T), "
         "c(T) (THER_NON_LINE)", "ther_3d", build_we11(TEMPDEP11),
         _we11_checks(TEMPDEP11[1]),
         note="als WE-11 met lambda 30 -> 25 W/m/K en c 500 -> 800 J/kg/K "
              "(300..1500 K); THER_NON_LINE per geboortestap met een stil "
              "THER_NL-materiaal (bridge_welding.ther_loop)"),
]


# ------------------------------------------------------------------
# WE-12: staaf die halverwege geboren wordt (mechanisch)
# ------------------------------------------------------------------
# Staaf 100 x 10 x 10 mm (hex8 4 x 1 x 1): deel A (x 0..50 mm) bestaat vanaf
# het begin, deel B (LASNAAD_B, x 50..100 mm) wordt op t = 6 s geboren
# (geboorte=ja, voorloop groot: in een keer). Opgelegde temperatuur op alle
# knopen: T0 -> T0 + 100 K (5 s), vast tot 6 s, dan terug naar T0 (11 s).
# x = 0 en x = L vast in x, y = 0 en z = 0 symmetrie.
DT12 = 100.0
HIST12 = [(0.0, T0), (5.0, T0 + DT12), (6.0, T0 + DT12), (11.0, T0)]
LB, BB = tm.LB, tm.BB
SIG12 = tm.E * tm.ALPHA * DT12           # 252 MPa
# Bron alleen voor de geboortetijd (de temperatuur is overal opgelegd)
SRC12 = ["ASTER-LASBRON Q=1, v=0.01, start=(0.05, 0, 0.01), "
         "richting=(1, 0, 0), lengte=0.05, r=0.04, t0=6, geboorte=ja, "
         "voorloop=1"]


def build_we12(plastic=False):
    def build():
        m = Model("T31", timesteps={"timeperiod": "%g s" % HIST12[-1][0],
                                    "stepsize": "1 s"})
        m.elset("A", "Staal")
        m.elset("LASNAAD_B", "Staal")
        ga = hexmesh(m, "A", 2, 1, 1, box(0, LB / 2, 0, BB, 0, BB))
        gb = hexmesh(m, "LASNAAD_B", 2, 1, 1, box(LB / 2, LB, 0, BB, 0, BB))
        m.mat("Staal",
              '<mechanical type="Isotropic" youngsmodulus="%g Pa" '
              'poissonratio="%g" thermalexpansion="%g 1/K" />'
              % (tm.E, tm.NU, tm.ALPHA),
              '<density type="Simple" density="%g kg/m³" />' % tm.RHO,
              '<thermal type="Isotropic" thermalconductivity="%g W/m/K" '
              'specificheat="%g J/kg/K" />' % (tm.K, tm.CP),
              *((PLASTIC12[plastic],) if plastic else ()))
        m.face_selection("X0", ga["faces"]["u0"])
        m.face_selection("XL", gb["faces"]["u1"])
        m.face_selection("Y0", ga["faces"]["v0"] + gb["faces"]["v0"])
        m.face_selection("Z0", ga["faces"]["w0"] + gb["faces"]["w0"])
        for sel, x, y, z in (("X0", 1, 0, 0), ("XL", 1, 0, 0),
                             ("Y0", 0, 1, 0), ("Z0", 0, 0, 1)):
            m.add('<displacement selection="%s"><value>0 m</value><x>%d</x>'
                  '<y>%d</y><z>%d</z></displacement>' % (sel, x, y, z))
        _initial(m)
        m.add('<temperature2 selection="Alle">%s<starttime>0 s</starttime>'
              '<endtime>0 s</endtime></temperature2>' % "".join(
                  '<value time="%g s">%.12g K</value>' % tv for tv in HIST12))
        m.add('<thermalstress><referencetemperature>%g K'
              '</referencetemperature></thermalstress>' % T0)
        _custom(m, SRC12)
        return m
    return build


# Plasticiteit: bilineair (Mecway BilinearIsotropic, H = tangent modulus) of
# een hardingscurve (Mecway Isotropic: rek = plastische rek, de standaard
# 'auto' van plastic_strain_definition bij een curve vanaf rek 0)
CURVE12 = (250e6, 270e6, 0.01)      # SY, spanning bij plastische rek 0,01
PLASTIC12 = {
    "bilineair": tm.PLASTIC8,
    "curve": '<plastic type="Isotropic"><hardeningcurve strain="0" '
             'stress="%g Pa" /><hardeningcurve strain="%g" stress="%g Pa" />'
             '</plastic>' % (CURVE12[0], CURVE12[2], CURVE12[1]),
}
MID_A = (LB / 4, BB / 2, BB)        # knoop midden op deel A
MID_B = (0.75 * LB, BB / 2, BB)     # knoop midden op deel B (alleen in B)


def _we12_checks(plastic):
    out = [
        Check("u_x(L/2) t=5 s (B ongeboren: A zet vrij uit)",
              lambda r: r.at("displx", (LB / 2, 0, 0), time=5.0),
              ref=tm.ALPHA * DT12 * LB / 2, tol=0.01, unit="m",
              source="vrij uitzetten van deel A: alpha dT L/2 (Timoshenko & "
                     "Goodier 1970, par. 150); B (E x 1e-4) houdt niet tegen"),
        Check("|sigma_xx| in ongeboren B t=5 s / (E alpha dT)",
              lambda r: abs(r.at("stressxx", MID_B, time=5.0)) / SIG12,
              ref=0.0, tol_abs=1e-3, unit="-", in_series=False,
              source="voor de geboorte geen spanning (tolerantie 1e-3 x E "
                     "alpha dT)"),
        Check("sigma_xx B t=7 s (1 s na de geboorte, 20 K afgekoeld)",
              lambda r: r.at("stressxx", MID_B, time=7.0),
              ref=tm.E * tm.ALPHA * 20.0, tol=0.01, unit="Pa",
              source="spanningsvrij geboren op T0 + 100 K: sigma = E alpha "
                     "(T_geboorte - T), hier 20 K"),
    ]
    if plastic:
        eps = tm.ALPHA * DT12
        if plastic == "curve":
            sy, hm = CURVE12[0], (CURVE12[1] - CURVE12[0]) / CURVE12[2]
        else:
            sy, hm = tm.SY, tm.HM
        p = (eps - sy / tm.E) / (1.0 + hm / tm.E)
        out.append(Check(
            "restspanning sigma_xx B t=11 s (%s)" % plastic,
            lambda r: r.at("stressxx", MID_B, time=11.0),
            ref=sy + hm * p, tol=0.01, unit="Pa",
            source="krimp alpha dT t.o.v. de geboortetemperatuur, lineaire "
                   "versteviging H: p = (alpha dT - SY/E)/(1 + H/E), sigma "
                   "= SY + H p"))
        out.append(Check(
            "eq. plastische rek B t=11 s",
            lambda r: r.at("equivalentplasticstrain", MID_B, time=11.0),
            ref=p, tol=0.01, unit="-", in_series=False, source="p, idem"))
    else:
        out += [
            Check("restspanning sigma_xx B t=11 s",
                  lambda r: r.at("stressxx", MID_B, time=11.0),
                  ref=SIG12, tol=0.01, unit="Pa",
                  source="krimpspanning t.o.v. de geboortetemperatuur: E "
                         "alpha (T_geboorte - T0) (zonder geboorte zou het 0 "
                         "zijn)"),
            Check("restspanning sigma_xx A t=11 s",
                  lambda r: r.at("stressxx", MID_A, time=11.0),
                  ref=SIG12, tol=0.01, unit="Pa",
                  source="evenwicht: dezelfde spanning in A"),
        ]
    return out


TESTS += [
    Case("WE-12", GROUP, "Toevoegmateriaal: staaf halverwege geboren, "
         "krimpspanning (elastisch)", "ther_3d", build_we12(False),
         _we12_checks(False),
         note="staaf 100 x 10 x 10 mm, hex8 4 x 1 x 1; B (x > 50 mm) geboren "
              "op t = 6 s bij T0 + 100 K, afkoelen tot T0 (11 s); x = 0 en "
              "x = L vast in x; temperatuur opgelegd (bron alleen voor de "
              "geboortetijd)"),
    Case("WE-12b", GROUP, "Toevoegmateriaal: staaf halverwege geboren, "
         "bilineair plastisch", "ther_3d", build_we12("bilineair"),
         _we12_checks("bilineair"),
         note="als WE-12 met SY = 250 MPa, H = 20 GPa (vloeit bij 99,2 K)"),
    Case("WE-12c", GROUP, "Toevoegmateriaal: staaf halverwege geboren, "
         "hardingscurve (TRACTION)", "ther_3d", build_we12("curve"),
         _we12_checks("curve"),
         note="als WE-12 met een Mecway-hardingscurve 250 -> 270 MPa bij "
              "plastische rek 0,01 (VMIS_ISOT_TRAC; stille kopie met de rek / "
              "1e-4)"),
]


# ------------------------------------------------------------------
# WE-13: twee lagen op een plaat (thermisch + mechanisch, kwalitatief)
# ------------------------------------------------------------------
# Halve plaat 30 x 8 x 4 mm (hex8 1 mm langs de las), laag 1 30 x 2,25 x 1,5
# mm (LASNAAD_1) en laag 2 erop 30 x 1,5 x 1,5 mm (LASNAAD_2). Laag 1 van t
# = 0 tot 4 s (x 5..25 mm), laag 2 van 6 tot 10 s; afkoelen tot 30 s (rand
# y = 8 mm op T0). Materiaal als WE-10 (E/SY/alpha(T), H = 1 GPa, TC_CUT).
X13 = _axis([(0.0, 0.030, 30)])
Y13 = _axis([(0.0, 0.00225, 3), (0.00225, 0.008, 4)])
Z13 = _axis([(-0.004, 0.0, 4)])
Y13L1 = _axis([(0.0, 0.00225, 3)])
Y13L2 = _axis([(0.0, 0.0015, 2)])
Z13L1 = _axis([(0.0, 0.0015, 2)])
Z13L2 = _axis([(0.0015, 0.003, 2)])
PASSES13 = [Pass("LASNAAD_1", 0.005, 0.020, 0.0),
            Pass("LASNAAD_2", 0.005, 0.020, 6.0)]
TEND13, DT13 = 30.0, 0.2
SRC13 = [
    "ASTER-LASBRON Q=%g, v=%g, start=(0.005, 0, 0.0015), richting=(1, 0, 0),"
    " lengte=0.02" % (QW, V),
    "ASTER-LASBRON af=%g, ar=%g, b=%g, c=%g, symmetrie=1, "
    "geboorte=<<LASNAAD_1>>" % (AF, AR, BB_, CC),
    "ASTER-LASBRON Q=%g, v=%g, start=(0.005, 0, 0.003), richting=(1, 0, 0),"
    " lengte=0.02, t0=6" % (QW, V),
    "ASTER-LASBRON af=%g, ar=%g, b=%g, c=%g, symmetrie=1, "
    "geboorte=<<LASNAAD_2>>" % (AF, AR, BB_, CC)]


def build_we13():
    m = Model("T31", timesteps={"timeperiod": "%g s" % TEND13,
                                "stepsize": "%g s" % DT13})
    m.elset("Plaat", "Staal")
    m.elset("LASNAAD_1", "Staal")
    m.elset("LASNAAD_2", "Staal")
    g = hexmesh(m, "Plaat", len(X13) - 1, len(Y13) - 1, len(Z13) - 1,
                _grid(X13, Y13, Z13))
    g1 = hexmesh(m, "LASNAAD_1", len(X13) - 1, len(Y13L1) - 1,
                 len(Z13L1) - 1, _grid(X13, Y13L1, Z13L1))
    g2 = hexmesh(m, "LASNAAD_2", len(X13) - 1, len(Y13L2) - 1,
                 len(Z13L2) - 1, _grid(X13, Y13L2, Z13L2))
    m.mat("Staal",
          '<mechanical type="Isotropic" youngsmodulus="%g Pa" '
          'poissonratio="%g" thermalexpansion="%g 1/K" />'
          % (tm.E, tm.NU, tm.ALPHA), *_thermal_mat())
    _initial(m)
    m.face_selection("Rand", g["faces"]["v1"])
    m.add('<temperature2 selection="Rand"><value>%.12g K</value><starttime>'
          '0 s</starttime><endtime>0 s</endtime></temperature2>' % T0)
    m.face_selection("Sym", g["faces"]["v0"] + g1["faces"]["v0"] +
                     g2["faces"]["v0"])
    m.add('<displacement selection="Sym"><value>0 m</value><x>0</x><y>1</y>'
          '<z>0</z></displacement>')
    a = m.node(0.0, 0.0, -0.004)
    b = m.node(0.030, 0.0, -0.004)
    m.node_selection("A", [a])
    m.node_selection("B", [b])
    for sel, x, z in (("A", 1, 0), ("A", 0, 1), ("B", 0, 1)):
        m.add('<displacement selection="%s"><value>0 m</value><x>%d</x>'
              '<y>0</y><z>%d</z></displacement>' % (sel, x, z))
    m.add('<thermalstress><referencetemperature>%g K</referencetemperature>'
          '</thermalstress>' % T0)
    _custom(m, SRC13 + tm.MAT10)
    return m


def _flag(cond):
    return 1.0 if cond else 0.0


P_L1 = (0.015, 0.0, 0.0015)       # bovenkant laag 1, midden van de las
P_L2 = (0.015, 0.0, 0.003)        # bovenkant laag 2


def _history(r, var, p):
    return [(t, r.at(var, p, time=t)) for t in r.times()]


def _layer2_top_nodes(r):
    return _top_nodes(r, "LASNAAD_2", 0.003)


def _section_force(r, xc=0.015, time=None):
    """int sigma_xx dA over de doorsnede x = xc (plaat + lagen; per element
    het gemiddelde van de knopen van zijn vlak x = xc)."""
    st = r.step(time)["stressxx"]
    nodes = r.model.nodes
    total = 0.0
    for _name, _m, elems in r.model.elsets:
        for _e, _s, nids, _a in elems:
            pts = [nodes[n - 1] for n in nids]
            face = [n for n, p in zip(nids, pts) if abs(p[0] - xc) < 1e-9]
            xs = [p[0] for p in pts]
            if len(face) != 4 or min(xs) < xc - 1e-9:
                continue          # alleen het element rechts van xc
            ys = [nodes[n - 1][1] for n in face]
            zs = [nodes[n - 1][2] for n in face]
            area = (max(ys) - min(ys)) * (max(zs) - min(zs))
            total += area * sum(st[n] for n in face) / 4.0
    return total


SY13 = tm.SY10[0]
A13 = 0.008 * 0.004 + 0.00225 * 0.0015 + 0.0015 * 0.0015


def _we13_checks():
    def peak(r, t_lo, t_hi):
        return max(v for t, v in _history(r, "temperature", P_L1)
                   if t_lo <= t <= t_hi)

    def t_at(r, t):
        return r.at("temperature", P_L1, time=t)
    return [
        Check("laag 1 smelt bij de eerste laag: T_max(x=15, boven laag 1) / "
              "1773 K", lambda r: peak(r, 0.0, 6.0) / 1773.0,
              ref=1.5, tol_abs=0.5, unit="-", in_series=False,
              source="smeltbad: de top van laag 1 komt boven ~1500 graden C "
                     "(kwalitatief: 1..2 x 1773 K)"),
        Check("laag 1 opnieuw opgewarmd door laag 2: T_max(6..12 s) > T(6 s)",
              lambda r: _flag(peak(r, 6.0, 12.0) > t_at(r, 6.0) + 200.0),
              ref=1.0, tol_abs=0.0, unit="-", in_series=False,
              source="meerlaags lassen: de vorige laag wordt opnieuw "
                     "opgewarmd (Radaj 1992)"),
        Check("laag 1 afgekoeld tussen de lagen: T(6 s) < T_max(0..6 s) / 2",
              lambda r: _flag(t_at(r, 6.0) < 0.5 * peak(r, 0.0, 6.0)),
              ref=1.0, tol_abs=0.0, unit="-", in_series=False,
              source="afkoelen tussen de lagen (2 s pauze)"),
        Check("max von Mises bovenop ongeboren laag 2 t=5.8 s / SY",
              lambda r: max(r.step(5.8)["vonmises"][n]
                            for n in _layer2_top_nodes(r)) / SY13,
              ref=0.0, tol_abs=1e-3, unit="-", in_series=False,
              source="voor de geboorte geen spanning in laag 2 (knopen die "
                     "alleen bij laag 2 horen; tolerantie 1e-3 SY)"),
        Check("T bovenop ongeboren laag 2 t=5.8 s",
              lambda r: max(r.step(5.8)["temperature"][n]
                            for n in _layer2_top_nodes(r)),
              ref=T0, tol_abs=1.0, unit="K", in_series=False,
              source="ongeboren: geen geleiding, begintemperatuur"),
        Check("evenwicht: int sigma_xx dA (x = 15 mm) t=30 s / (SY A)",
              lambda r: _section_force(r, time=TEND13) / (SY13 * A13),
              ref=0.0, tol_abs=0.03, unit="-", in_series=False,
              source="vrije plaat: netto langskracht 0 (3% van SY A)"),
        Check("sigma_xx bovenop laag 2 (x = 15 mm) t=30 s > 0 (trek)",
              lambda r: _flag(r.at("stressxx", P_L2, time=TEND13) > 0.0),
              ref=1.0, tol_abs=0.0, unit="-", in_series=False,
              source="Radaj (1992), Masubuchi (1980): langsspanning in de las "
                     "trek"),
        Check("sigma_xx rand y = 8 mm midden dikte t=30 s < 0 (druk)",
              lambda r: _flag(r.at("stressxx", (0.015, 0.008, -0.002),
                                   time=TEND13) < 0.0),
              ref=1.0, tol_abs=0.0, unit="-", in_series=False,
              source="naast de las druk (evenwicht)"),
    ]


TESTS.append(Case(
    "WE-13", GROUP, "Toevoegmateriaal: twee lagen op een plaat, "
    "temperatuurgeschiedenis en restspanning", "ther_3d", build_we13,
    _we13_checks(), timeout=3600,
    note="halve plaat 30 x 8 x 4 mm + laag 1 (30 x 2,25 x 1,5 mm) en laag 2 "
         "(30 x 1,5 x 1,5 mm), elk geboren met de eigen bron (500 W, 5 mm/s, "
         "x 5..25 mm; laag 1 0..4 s, laag 2 6..10 s), afkoelen tot 30 s, dt "
         "= 0,2 s; E/SY/alpha(T) als WE-10; hex8 3D (locking, kwalitatief)"))
