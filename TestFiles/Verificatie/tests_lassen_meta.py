# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Lassen (stap S5): metallurgie van staal (#ASTER-METALLURGIE ->
CALC_META met het TRC-diagram en META_ACIER van een standaard-dataset;
bridge_metallurgy.py), en als optie de metallurgische wet in de keten
thermisch -> mechanisch (META_P_IL met AFFE_VARC M_ACIER).

Staaf 100 x 10 x 10 mm (hex8 4 x 1 x 1) met een opgelegde temperatuur op
alle knopen (tijdstabel), zoals de testcases hsnv101a / mtlp100b (de
temperatuur is overal gelijk, dus het fasenverloop is dat van een punt).

  MT-01  hsnv101a in 3D: afkoelen 900 -> 0 graden C in 180 s, austeniet ->
         bainiet; dataset HSNV101A (TRC, META_ACIER en de mechanische
         gegevens van de testcase), mechanica=ja, trekspanning 0 -> 360 MPa
         (60 s). Referentie: de fasen van hsnv101a zelf (AXIS, META_NOEU,
         gedraaid met Code-Aster 17.4: de testcase heeft geen TEST_RESU op
         de fasen) en de TEST_RESU-waarden van hsnv101a (SIYY, EPYY).
  MT-02  snelle afkoeling (200 C/s, sneller dan de snelste TRC-geschiedenis
         met omzetting) van 16MND5: alleen martensiet, Koistinen-Marburger
         Z_M = 1 - exp(ALPHA (Ms - T)) met de hand; hardheid als lineaire
         mengregel van DURT_META.
  MT-03  austenitiseren bij opwarmen (16MND5, ferriet, 10 C/s tot boven
         Ac3): Z_a uit dZ/dt = (Z_eq(T) - Z)/tau met de hand.
  MT-04  transformatierek: vrije staaf, 16MND5 mechanica=ja, opwarmen tot
         austeniet en afschrikken tot martensiet; rek met de hand uit de
         formule van ELAS_META_FO (PHASE_REFE FROID).
  MT-05  idem, ingeklemd in x en elastisch: sigma_xx = -E eps_th.
"""
import math
from verif_lib import Case, Check, Model, hexmesh, box

GROUP = "Metallurgie"

LB, BB = 0.1, 0.01
E, NU = 200e9, 0.3
K, RHO, CP = 50.0, 7800.0, 500.0
KELVIN = 273.15
MID = (LB / 2, BB / 2, BB)


def _bar(m, alpha=1.5e-5, plastic=""):
    m.elset("Staaf", "Staal")
    g = hexmesh(m, "Staaf", 4, 1, 1, box(0, LB, 0, BB, 0, BB))
    lines = ['<mechanical type="Isotropic" youngsmodulus="%g Pa" '
             'poissonratio="%g" thermalexpansion="%g 1/K" />' % (E, NU, alpha),
             '<density type="Simple" density="%g kg/m³" />' % RHO,
             '<thermal type="Isotropic" thermalconductivity="%g W/m/K" '
             'specificheat="%g J/kg/K" />' % (K, CP)]
    if plastic:
        lines.append(plastic)
    m.mat("Staal", *lines)
    return g


def build_bar(history_c, dt, custom, mech=None, clamped=False, load=None,
              plastic=""):
    """history_c: [(t, T in graden C)] op alle knopen. mech: None (alleen
    thermisch + metallurgie) of de referentietemperatuur [C] van de
    Thermal stress. load: [(t, trekspanning Pa)] op x = L."""
    def build():
        m = Model("T31", timesteps={"timeperiod": "%g s" % history_c[-1][0],
                                    "stepsize": "%g s" % dt})
        g = _bar(m, plastic=plastic)
        m.node_selection("Alle", list(range(1, len(m.nodes) + 1)))
        m.add('<temperature2 selection="Alle">%s<starttime>0 s</starttime>'
              '<endtime>0 s</endtime></temperature2>' % "".join(
                  '<value time="%g s">%.12g K</value>' % (t, v + KELVIN)
                  for t, v in history_c))
        m.add(u'<initialconditions>\n    <temperature2 selection="Alle">'
              u'<value>%.12g K</value><starttime>0 s</starttime><endtime>'
              u'0 s</endtime></temperature2>\n  </initialconditions>' %
              (history_c[0][1] + KELVIN))
        if mech is not None:
            m.face_selection("X0", g["faces"]["u0"])
            m.face_selection("Y0", g["faces"]["v0"])
            m.face_selection("Z0", g["faces"]["w0"])
            m.face_selection("XL", g["faces"]["u1"])
            faces = [("X0", 1, 0, 0), ("Y0", 0, 1, 0), ("Z0", 0, 0, 1)]
            if clamped:
                faces.append(("XL", 1, 0, 0))
            for sel, x, y, z in faces:
                m.add('<displacement selection="%s"><value>0 m</value>'
                      '<x>%d</x><y>%d</y><z>%d</z></displacement>' %
                      (sel, x, y, z))
            if load:
                m.add('<traction selection="XL">%s<y>0 Pa</y><z>0 Pa</z>'
                      '</traction>' % "".join(
                          '<x time="%g s">%.12g Pa</x>' % tv for tv in load))
            m.add('<thermalstress><referencetemperature>%.12g K'
                  '</referencetemperature></thermalstress>' % (mech + KELVIN))
        m.add(u"<ccxmodeldefinition>%s</ccxmodeldefinition>" % "\n".join(
            u"** " + line for line in custom))
        return m
    return build


def _strain_x(r, t):
    """eps_xx = u_x(L) / L (uniaxiaal, homogeen)."""
    return r.at("displx", (LB, BB, BB), time=t) / LB


# ------------------------------------------------------------------
# MT-01: hsnv101a in 3D
# ------------------------------------------------------------------
# Fasen van hsnv101a (AXIS, Code-Aster 17.4, META_NOEU in knoop N13; de
# testcase zelf gedraaid met een extra CREA_TABLE): bainiet bij t (T = 900
# - 5 t graden C).
HSNV101A_BAINITE = [(64.0, 0.0784006), (80.0, 0.392003), (90.0, 0.588005),
                    (100.0, 0.784006), (110.0, 0.980008), (112.0, 1.0)]
# TEST_RESU van hsnv101a (VALE_CALC, knoop N13 van M1): SIYY en EPYY
HSNV101A_SIG = [(47.0, 2.8199986897636e08), (64.0, 3.5999984904968e08),
                (114.0, 3.5999984462008e08)]
HSNV101A_EPS = [(47.0, -4.1125003852496e-03), (48.0, -9.3464484572225e-04),
                (64.0, 0.034690852679456), (114.0, 0.036841623953657)]

MT01 = Case(
    "MT-01", GROUP, "Metallurgie: hsnv101a in 3D (bainiet, META_P_IL)",
    "ther_3d",
    build_bar([(0.0, 900.0), (176.0, 900.0 - 5.0 * 176.0)], 1.0,
              ["ASTER-METALLURGIE Staal: staal=HSNV101A, begin=austeniet, "
               "mechanica=ja"],
              mech=900.0, load=[(0.0, 0.0), (60.0, 360e6)]),
    [Check("bainiet t=%g s" % t, (lambda t: lambda r: r.at(
        "bainite", MID, time=t))(t), ref=z, tol=0.01, unit="-",
        source="hsnv101a (AXIS, Code-Aster 17.4): META_NOEU V3 in N13")
     for t, z in HSNV101A_BAINITE] +
    [Check("austeniet t=100 s", lambda r: r.at("austenite", MID, time=100.0),
           ref=1.0 - 0.784006, tol=0.01, unit="-",
           source="hsnv101a: META_NOEU V5 = 1 - V3")] +
    [Check("sigma_xx t=%g s" % t, (lambda t: lambda r: r.at(
        "stressxx", MID, time=t))(t), ref=s, tol=0.01, unit="Pa",
        source="hsnv101a TEST_RESU SIYY (N13)") for t, s in HSNV101A_SIG] +
    [Check("eps_xx t=%g s" % t, (lambda t: lambda r: _strain_x(r, t))(t),
           ref=e, tol=0.01, unit="-",
           source="hsnv101a TEST_RESU EPYY (N13)") for t, e in HSNV101A_EPS],
    timeout=1800,
    note="staaf 100 x 10 x 10 mm, 4 hex8; T = 900 - 5 t graden C op alle "
         "knopen, dt = 1 s; #ASTER-METALLURGIE staal=HSNV101A, "
         "begin=austeniet, mechanica=ja; Thermal stress T_ref = 900 C; "
         "traction 0 -> 360 MPa (0..60 s) op x = L, symmetrievlakken")


# ------------------------------------------------------------------
# MT-02: martensiet (Koistinen-Marburger), 16MND5
# ------------------------------------------------------------------
MS0_16, KM_16 = 415.0, -0.0249        # DATASETS["16MND5"] (mtlp103a)
HV_16 = (184.0, 184.0, 309.0, 450.0, 100.0)
RATE2 = 1000.0                        # C/s


def _km(t_c):
    """Z_M = 1 - exp(ALPHA (Ms - T)) onder Ms (ALPHA < 0)."""
    return 1.0 - math.exp(KM_16 * (MS0_16 - t_c)) if t_c < MS0_16 else 0.0


def _t2(time):
    return 900.0 - RATE2 * time


MT02_TIMES = (0.5, 0.54, 0.58, 0.62, 0.66, 0.88)  # T = 400 .. 240, 20 C

MT02 = Case(
    "MT-02", GROUP, "Metallurgie: snel afkoelen, martensiet "
    "(Koistinen-Marburger)", "ther_3d",
    build_bar([(0.0, 900.0), (0.88, 20.0)], 0.02,
              ["ASTER-METALLURGIE Staal: staal=16MND5, begin=austeniet"]),
    [Check("martensiet T=%g C" % _t2(t), (lambda t: lambda r: r.at(
        "martensite", MID, time=t))(t), ref=_km(_t2(t)), tol=0.01, unit="-",
        source="Koistinen-Marburger 1 - exp(ALPHA (MS0 - T)), ALPHA = "
               "-0.0249, MS0 = 415 C (dataset 16MND5, mtlp103a)")
     for t in MT02_TIMES] +
    [Check("bainiet + ferriet T=20 C", lambda r: r.at(
        "bainite", MID, time=0.88) + r.at("ferrite", MID, time=0.88),
        ref=0.0, tol_abs=1e-6, unit="-",
        source="1000 C/s: sneller dan de TRC-geschiedenissen met omzetting "
               "(bij 200 C/s vormt Code-Aster nog 3% bainiet, gemeten)"),
     Check("hardheid T=280 C", lambda r: r.at("hardness", MID, time=0.62),
           ref=_km(280.0) * HV_16[3] + (1.0 - _km(280.0)) * HV_16[4],
           tol=0.01, unit="HV",
           source="mengregel Z_M HV_M + Z_A HV_A (DURT_META 450 / 100 HV)")],
    note="staaf 4 hex8; T = 900 - 1000 t graden C op alle knopen (0..0,88 "
         "s), dt = 0,02 s; staal=16MND5, begin=austeniet")


# ------------------------------------------------------------------
# MT-03: austenitiseren bij opwarmen, 16MND5
# ------------------------------------------------------------------
AC1_16, AC3_16, TAU_16 = 724.0, 846.0, 0.34
RATE3 = 10.0                              # C/s
T03 = [(0.0, 600.0), (40.0, 1000.0)]


def _aust_lin(time, t_start, rate):
    """Z_a(t) bij T = t_start + rate t vanuit ferriet: dZ/dt = (Z_eq(T) -
    Z)/tau met Z_eq lineair tussen Ac1 en Ac3 en tau constant (TAUX_1 =
    TAUX_3); analytisch Z = Z_eq(t) - a tau (1 - exp(-(t - t1)/tau))
    tussen Ac1 en Ac3 (a = dZ_eq/dt), daarboven exponentieel naar 1."""
    a = rate / (AC3_16 - AC1_16)
    t1 = (AC1_16 - t_start) / rate
    t3 = (AC3_16 - t_start) / rate
    if time <= t1:
        return 0.0
    if time <= t3:
        return (a * (time - t1) -
                a * TAU_16 * (1.0 - math.exp(-(time - t1) / TAU_16)))
    z3 = _aust_lin(t3, t_start, rate)
    return 1.0 - (1.0 - z3) * math.exp(-(time - t3) / TAU_16)


def _aust(time):
    return _aust_lin(time, 600.0, RATE3)


MT03_TIMES = (14.0, 18.0, 22.0, 24.0, 25.0, 26.0)

MT03 = Case(
    "MT-03", GROUP, "Metallurgie: austenitiseren bij opwarmen", "ther_3d",
    build_bar(T03, 0.1,
              ["ASTER-METALLURGIE Staal: staal=16MND5, begin=ferriet"]),
    [Check("austeniet T=%g C" % (600.0 + RATE3 * t),
           (lambda t: lambda r: r.at("austenite", MID, time=t))(t),
           ref=_aust(t), tol=0.01, unit="-",
           source="dZ/dt = (Z_eq(T) - Z)/tau, Z_eq lineair Ac1..Ac3 = "
                  "724..846 C, tau = 0,34 s (dataset 16MND5)")
     for t in MT03_TIMES] +
    [Check("ferriet T=1000 C", lambda r: r.at("ferrite", MID, time=40.0),
           ref=0.0, tol_abs=1e-3, unit="-",
           source="boven Ac3 volledig austeniet")],
    note="staaf 4 hex8; T = 600 + 10 t graden C (0..40 s), dt = 0,1 s; "
         "staal=16MND5, begin=ferriet")


# ------------------------------------------------------------------
# MT-04 / MT-05: transformatierek (16MND5 + eigen mechanische gegevens)
# ------------------------------------------------------------------
ALPHA_F4, ALPHA_A4, DEPS4 = 1.4e-5, 2.2e-5, 1.0e-2
TREF4 = 20.0
# opwarmen 20 -> 1000 C (10 s), 2 s vasthouden, afschrikken 1000 C/s
T04 = [(0.0, 20.0), (10.0, 1000.0), (12.0, 1000.0), (12.98, 20.0)]


def _t4(time):
    for (t0, v0), (t1, v1) in zip(T04, T04[1:]):
        if time <= t1:
            return v0 + (v1 - v0) * (time - t0) / (t1 - t0)
    return T04[-1][1]


def _eps4(t_c, z_a):
    """ELAS_META_FO, PHASE_REFE FROID (begin ferriet bij T_ref):
    eps = (1 - Z_a) alpha_f (T - T_ref) + Z_a (alpha_a (T - T_ref) - DEPS)."""
    d = t_c - TREF4
    return (1.0 - z_a) * ALPHA_F4 * d + z_a * (ALPHA_A4 * d - DEPS4)


Z_HOT4 = _aust_lin(12.0 - 2.0, 20.0, 98.0)      # na het opwarmen (10 s)
Z_HOT4 = 1.0 - (1.0 - Z_HOT4) * math.exp(-2.0 / TAU_16)   # 2 s vasthouden
MT04_POINTS = [  # (tijd, Z_a met de hand)
    (12.0, Z_HOT4),                          # 1000 C, austeniet
    (12.4, Z_HOT4),                          # 600 C, boven Ms
    (12.6, Z_HOT4 * (1.0 - _km(_t4(12.6)))),   # 400 C
    (12.66, Z_HOT4 * (1.0 - _km(_t4(12.66)))),  # 340 C
    (12.7, Z_HOT4 * (1.0 - _km(_t4(12.7))))]    # 300 C
META4 = ["ASTER-METALLURGIE Staal: staal=16MND5, begin=ferriet, "
         "mechanica=ja,",
         "ASTER-METALLURGIE     ALPHA_F=%g, ALPHA_A=%g, DEPS=%g, SY=1e10" %
         (ALPHA_F4, ALPHA_A4, DEPS4)]

MT04 = Case(
    "MT-04", GROUP, "Metallurgie: transformatierek, vrije staaf",
    "ther_3d", build_bar(T04, 0.02, META4, mech=TREF4),
    [Check("eps_xx t=%g s (T=%g C)" % (t, _t4(t)),
           (lambda t: lambda r: _strain_x(r, t))(t),
           ref=_eps4(_t4(t), z), tol=0.01, unit="-",
           source="(1 - Z_a) alpha_f dT + Z_a (alpha_a dT - DEPS), Z_a "
                  "uit Koistinen-Marburger (ELAS_META_FO FROID)")
     for t, z in MT04_POINTS if abs(_eps4(_t4(t), z)) > 1e-5] +
    [Check("eps_xx na afschrikken", lambda r: _strain_x(r, 12.98),
           ref=_eps4(20.0, 1.0 - _km(20.0)), tol_abs=0.01 * DEPS4, unit="-",
           source="T = T_ref, martensiet: rek ~ 0 (koude fasen hebben "
                  "dezelfde rek; tolerantie 1% van DEPS: Code-Aster laat "
                  "~0,7% austeniet over, MT-02)"),
     Check("max |sigma_xx| t=12,6 s", lambda r: max(
         abs(v) for v in r.step(12.6)["stressxx"].values()),
         ref=0.0, tol_abs=1e5, unit="Pa", source="vrije staaf: geen "
         "spanning (tolerantie 0,1 MPa)")],
    note="staaf 4 hex8; T: 20 -> 1000 C (10 s), 2 s vast, -> 20 C in "
         "0,98 s (1000 C/s); dt = 0,02 s; 16MND5 mechanica=ja, ALPHA_F = "
         "1,4e-5, ALPHA_A = 2,2e-5, DEPS = 1e-2, SY = 1e10 (elastisch); "
         "Thermal stress T_ref = 20 C; symmetrievlakken")

MT05 = Case(
    "MT-05", GROUP, "Metallurgie: transformatierek, ingeklemde staaf "
    "(elastisch)", "ther_3d",
    build_bar(T04, 0.02, META4, mech=TREF4, clamped=True),
    [Check("sigma_xx t=%g s (T=%g C)" % (t, _t4(t)),
           (lambda t: lambda r: r.at("stressxx", MID, time=t))(t),
           ref=-E * _eps4(_t4(t), z), tol=0.01, unit="Pa",
           source="-E eps_th (eenassig verhinderd, dwars vrij), eps_th "
                  "zoals MT-04")
     for t, z in MT04_POINTS if abs(_eps4(_t4(t), z)) > 1e-5],
    note="als MT-04, x = L vast in x")


# ------------------------------------------------------------------
# MT-06: lasrups op een plaat (WE-10, alleen thermisch) met 16MND5:
# warmte-beinvloede zone en hardheid, consistent met Ac1/Ac3
# ------------------------------------------------------------------
import tests_lassen_mech as WE

META6 = ["ASTER-METALLURGIE Staal: staal=16MND5, begin=ferriet, "
         "TAUX_1=0.001, TAUX_3=0.001"]


def build_mt06():
    m = WE.build_we10()
    # alleen thermisch + metallurgie: Thermal stress, opleggingen en
    # #ASTER-MATERIAAL eruit, de metallurgie erbij
    m.items = [x for x in m.items if "thermalstress" not in x and
               "ccxmodeldefinition" not in x and "<displacement" not in x]
    m.add(u"<ccxmodeldefinition>%s</ccxmodeldefinition>" % "\n".join(
        u"** " + line for line in WE.SRC10 + META6))
    return m


def _line_nodes(r):
    """Knopen op de lijn x = 20 mm, z = 0 (bovenvlak), dwars op de las."""
    return [r.nearest((WE.MID10, y, 0.0)) for y in WE.Y10S]


def _peak(r, nid, var):
    return max(d[var][nid] for _t, d in r.steps if var in d)


def _haz_violations(r):
    """Knopen waar de fasen niet kloppen met de piektemperatuur: onder Ac1
    moet het ferriet blijven (1), boven Ac3 + 50 C volledig austeniet
    worden (max >= 0,99)."""
    bad = 0
    for nid in _line_nodes(r):
        tp = _peak(r, nid, "temperature") - KELVIN
        if tp < AC1_16 - 1.0 and abs(r.node("ferrite", nid) - 1.0) > 1e-9:
            bad += 1
        if tp > AC3_16 + 50.0 and _peak(r, nid, "austenite") < 0.99:
            bad += 1
    return bad


def _haz_width(r):
    """Breedte [m] van de zone die (deels) austeniet werd (max > 1%)."""
    ys = [WE.Y10S[i] for i, nid in enumerate(_line_nodes(r))
          if _peak(r, nid, "austenite") > 0.01]
    return max(ys) if ys else 0.0


WELD6 = (WE.MID10, 0.0, 0.0)
BASE6 = (WE.MID10, WE.BY10, -WE.TZ10)

MT06 = Case(
    "MT-06", GROUP, "Metallurgie: lasrups op een plaat, warmte-beinvloede "
    "zone (16MND5)", "ther_3d", build_mt06,
    [Check("knopen in strijd met Ac1/Ac3 (lijn x = 20 mm, boven)",
           _haz_violations, ref=0.0, tol_abs=0.0, unit="-", in_series=False,
           source="piektemperatuur < Ac1 (724 C): ferriet blijft; > Ac3 + "
                  "50 C (846 C): volledig austeniet (dataset 16MND5)"),
     Check("ferriet basis (y = 15 mm, onder)", lambda r: r.at(
         "ferrite", BASE6), ref=1.0, tol_abs=1e-9, unit="-",
           in_series=False, source="koud gebleven: geen omzetting"),
     Check("som fasen lasas (eind)", lambda r: sum(
         r.at(p, WELD6) for p in ("ferrite", "pearlite", "bainite",
                                  "martensite", "austenite")),
           ref=1.0, tol_abs=1e-6, unit="-", in_series=False,
           source="fracties tellen op tot 1"),
     Check("austeniet lasas (eind, afgekoeld)", lambda r: r.at(
         "austenite", WELD6), ref=0.0, tol_abs=0.01, unit="-",
           in_series=False, source="afgekoeld tot ~T0: omgezet"),
     Check("hardheid lasas > basis (184 HV)", lambda r: 1.0 if r.at(
         "hardness", WELD6) > HV_16[0] + 1.0 else 0.0, ref=1.0, tol_abs=0.0,
           unit="-", in_series=False,
           source="harde fasen (bainiet 309, martensiet 450 HV) in de las "
                  "(kwalitatief; Radaj 1992)"),
     Check("breedte austenietzone boven (x = 20 mm)", _haz_width,
           ref=0.003, tol_abs=0.002, unit="m", in_series=False,
           source="kwalitatief: van de orde van de bronbreedte b = 2 mm "
                  "(tussen 1 en 5 mm)")],
    timeout=1800,
    note="WE-10 (halve plaat 40 x 15 x 5 mm, Goldak 1000 W, 5 mm/s) alleen "
         "thermisch, met #ASTER-METALLURGIE staal=16MND5, begin=ferriet")


# ------------------------------------------------------------------
# MT-07: lasrups WE-10 met de metallurgische wet (restspanning met
# transformatierek en vloeigrens per fase), kwalitatief t.o.v. WE-10
# ------------------------------------------------------------------
# Mechanische gegevens per fase: illustratief (geen meetwaarden). ALPHA_A
# en DEPS volgen uit hsnv101a (alpha_f 15e-6, alpha_a 23,5e-6, DEPS 2,52e-3
# bij 900 C) omgerekend naar T_ref = 20 C: DEPS(20 C) = 2,52e-3 + (15e-6 -
# 23,5e-6) (20 - 900) = 1,0e-2. Vloeigrens koude fasen: de WE-10-tabel
# (#ASTER-MATERIAAL) voor ferriet/perliet, hoger voor bainiet en
# martensiet, lager voor austeniet.
SYB7 = (550e6, 480e6, 400e6, 220e6, 75e6, 35e6, 15e6, 10e6)
SYM7 = (900e6, 800e6, 650e6, 350e6, 120e6, 50e6, 20e6, 10e6)
SYA7 = (200e6, 170e6, 140e6, 90e6, 40e6, 20e6, 12e6, 10e6)
META7 = ["ASTER-METALLURGIE Staal: staal=16MND5, begin=ferriet, "
         "TAUX_1=0.001, TAUX_3=0.001, mechanica=ja,",
         "ASTER-METALLURGIE     ALPHA_A=2.35e-5, DEPS=1.0e-2, TC=(%s),"
         % ", ".join("%g" % t for t in WE.TC10),
         "ASTER-METALLURGIE     SY_B=(%s)," % ", ".join("%g" % v
                                                     for v in SYB7),
         "ASTER-METALLURGIE     SY_M=(%s)," % ", ".join("%g" % v
                                                     for v in SYM7),
         "ASTER-METALLURGIE     SY_A=(%s)" % ", ".join("%g" % v
                                                    for v in SYA7)]
# WE-10 (hex8 3D, zelfde plaat zonder metallurgie): sigma_xx lasas boven
# x = 20 mm, t = 60 s (docs/Lassen.md par. 9.5)
SXX_WE10 = 466e6


def build_mt07():
    m = WE.build_we10()
    m.items = [x for x in m.items if "ccxmodeldefinition" not in x]
    m.add(u"<ccxmodeldefinition>%s</ccxmodeldefinition>" % "\n".join(
        u"** " + line for line in WE.SRC10 + WE.MAT10 + META7))
    return m


MT07 = Case(
    "MT-07", GROUP, "Metallurgie: lasrups met metallurgische wet "
    "(restspanning)", "ther_3d", build_mt07,
    [Check("evenwicht: int sigma_xx dA (x = 20 mm) / (SY A)",
           lambda r: WE._section(r) / (WE.SY10[0] * WE.BY10 * WE.TZ10),
           ref=0.0, tol_abs=0.03, unit="-", in_series=False,
           source="vrije plaat: netto langskracht 0 (3% van SY A)"),
     Check("evenwicht: int sigma_xx (z - zc) dA / (SY A t)",
           lambda r: WE._section(r, moment=True) /
           (WE.SY10[0] * WE.BY10 * WE.TZ10 * WE.TZ10),
           ref=0.0, tol_abs=0.03, unit="-", in_series=False,
           source="vrije plaat: netto buigend moment 0"),
     Check("martensiet + bainiet lasas boven", lambda r: r.at(
         "martensite", WELD6) + r.at("bainite", WELD6), ref=1.0,
           tol_abs=0.02, unit="-", in_series=False,
           source="zoals MT-06: de las is volledig omgezet"),
     Check("sigma_xx lasas boven < WE-10 (466 MPa)", lambda r: 1.0 if r.at(
         "stressxx", WELD6) < SXX_WE10 else 0.0, ref=1.0, tol_abs=0.0,
           unit="-", in_series=False,
           source="de uitzetting bij de martensietvorming (onder Ms = "
                  "415 C) verlaagt de langsspanning in de las (Radaj "
                  "1992; kwalitatief)")],
    timeout=3600,
    note="WE-10 met #ASTER-METALLURGIE staal=16MND5, begin=ferriet, TAUX "
         "0,001 s, mechanica=ja (ALPHA_A 2,35e-5, DEPS 1e-2, SY per fase "
         "illustratief); hex8 3D")


# ------------------------------------------------------------------
# MT-08: eigen TRC-diagram (trc=... uit #ASTER-code) met eigen
# META_ACIER-waarden: Koistinen-Marburger met MS0 = 300 C, ALPHA = -0,011
# ------------------------------------------------------------------
MS8, KM8 = 300.0, -0.011
TRC8 = ["ASTER TRC_EIGEN = DEFI_TRC(",
        "ASTER     HIST_EXP=(",
        "ASTER         _F(VALE=(-100.0, 11.0, 0., 0., 0., 0., 0., 0.,",
        "ASTER                  0., 0., 0., 830.0)),",
        "ASTER         _F(VALE=(-400.0, 11.0, 0., 0., 0., 0., 0., 0.,",
        "ASTER                  0., 0., 0., 830.0)),",
        "ASTER     ),",
        "ASTER     TEMP_MS=_F(SEUIL=0.38, AKM=-30.791, BKM=11.698, "
        "TPLM=-1.0),",
        "ASTER )",
        "ASTER-METALLURGIE Staal: trc=TRC_EIGEN, AR3=830, MS0=%g, "
        "ALPHA=%g," % (MS8, KM8),
        "ASTER-METALLURGIE     AC1=724, AC3=846, TAUX_1=0.34, TAUX_3=0.34, "
        "begin=austeniet"]


def _km8(t_c):
    return 1.0 - math.exp(KM8 * (MS8 - t_c)) if t_c < MS8 else 0.0


MT08_TIMES = (6.5, 7.5, 8.5, 8.8)          # T = 250, 150, 50, 20 C

MT08 = Case(
    "MT-08", GROUP, "Metallurgie: eigen TRC-diagram (#ASTER-code)",
    "ther_3d", build_bar([(0.0, 900.0), (8.8, 20.0)], 0.1, TRC8),
    [Check("martensiet T=%g C" % (900.0 - 100.0 * t),
           (lambda t: lambda r: r.at("martensite", MID, time=t))(t),
           ref=_km8(900.0 - 100.0 * t), tol=0.01, unit="-",
           source="Koistinen-Marburger 1 - exp(ALPHA (MS0 - T)), ALPHA = "
                  "-0,011, MS0 = 300 C (eigen META_ACIER)")
     for t in MT08_TIMES] +
    [Check("bainiet + ferriet T=20 C", lambda r: r.at(
        "bainite", MID, time=8.8) + r.at("ferrite", MID, time=8.8),
        ref=0.0, tol_abs=1e-6, unit="-",
        source="TRC zonder omzetting (alleen de geschiedenissen -100 en "
               "-400 C/s van 16MND5)")],
    note="staaf 4 hex8; T = 900 - 100 t graden C (0..8,8 s), dt = 0,1 s; "
         "DEFI_TRC in #ASTER-code, #ASTER-METALLURGIE trc=TRC_EIGEN met "
         "eigen META_ACIER")


TESTS = [MT01, MT02, MT03, MT04, MT05, MT06, MT07, MT08]
