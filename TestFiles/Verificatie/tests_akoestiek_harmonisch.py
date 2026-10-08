# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Akoestiek harmonisch: drukveld per frequentie (Mecway A30-model ->
PHENOMENE='ACOUSTIQUE', DYNA_VIBRA HARM; bridge-analyse acoustic_harm_3d,
bridge_acoustic.build_acoustic_harm). Toelichting: docs/Verificatie.md
par. 7 (AH) en docs/Akoestiek.md.

Referenties (vlakke golven in een buis, lengte L langs x, p(t) =
Re(P e^(i w t)), k = w / c, bron: normaalsnelheid v0 de fluid in op x = 0;
1D Helmholtz P'' + k^2 P = 0 met i w rho u = -dP/dx, zelf afgeleid,
standaardresultaat o.a. Kinsler, Frey, Coppens & Sanders, Fundamentals of
Acoustics, niet lokaal ingezien):
  - absorberend einde (Z = rho c): P = rho c v0 e^(-i k x), |P| constant,
    fase -k x;
  - harde wand op x = L: P = -i rho c v0 cos(k (L - x)) / sin(k L);
  - impedantie Z op x = L: P = A (e^(-ikx) + R e^(-ik(2L - x))),
    R = (Z - rho c) / (Z + rho c), A = rho c v0 / (1 - R e^(-2ikL));
  - opgelegde druk p0 op x = 0 en absorberend einde: P = p0 e^(-ikx);
  - Code-Aster-testcase ahlv100a (golfgeleider 1 x 0,1 x 0,2 m, hex20,
    rho 1,3, c 343, VNOR 0,014 m/s, 500 Hz): TEST_RESU-waarden (VALE_REFE
    analytisch, VALE_CALC van Code-Aster zelf), nagebouwd op hetzelfde
    3D-mesh; ahlv100a legt VNOR naar buiten op (EXCIT COEF_MULT=-1), de
    bridge de fluid in: P_bridge = -P_ahlv100a.
Geluidsdrukniveau: SPL = 20 log10(|P| / (sqrt(2) 20 uPa)); Code-Aster's DB
(PRAC_ELNO) is 20 log10(|P| / 20 uPa), 3,0103 dB hoger.
"""
import io
import csv
import math
import cmath
from verif_lib import Case, Check, Model, hexmesh, box

GROUP = "Akoestiek harmonisch"
RHO, C = 1.2, 340.0
AIR = ['<mechanical type="Isotropic" speedofsound="340 m/s" />',
       u'<density type="Simple" density="1.2 kg/m³" />']
V0 = 1.0e-3              # source_velocity (standaard)
P_REF = 2.0e-5
LT, WT = 1.0, 0.1        # buis 1 x 0,1 x 0,1 m (dwarsmode vanaf 1700 Hz)


def _spl(amp):
    return 20.0 * math.log10(amp / math.sqrt(2.0) / P_REF)


# ------------------------------------------------------------------
# Resultaten: <mode id=k frequency=f> met pressure (Re) en pressurei (Im)
# ------------------------------------------------------------------
def _k_of(r, f):
    """mode-id van frequentie f (de stappen hebben t = mode-id)."""
    return min(r.modes, key=lambda k: abs(r.modes[k] - f))


def _p(r, f, nid):
    st = r.step(float(_k_of(r, f)))
    return complex(st["pressure"][nid], st["pressurei"][nid])


def _axis(r):
    """[(x, nid)] van de knoopselectie As (y = z = 0), oplopend in x."""
    return sorted((r.model.nodes[n - 1][0], n) for n in r.model.nodesel["As"])


def _amp_range(r, f, fn):
    return fn(abs(_p(r, f, n)) for _x, n in _axis(r))


def _phase_lag(r, f):
    """Opgetelde faseval van x = 0 tot x = L langs de as [rad]."""
    pts = _axis(r)
    lag = 0.0
    for (_x0, a), (_x1, b) in zip(pts, pts[1:]):
        lag -= cmath.phase(_p(r, f, b) / _p(r, f, a))
    return lag


def _max_dev(r, f, exact):
    """max |P_fem - P_exact| / max |P_exact| langs de as."""
    pts = _axis(r)
    ref = max(abs(exact(x)) for x, _n in pts)
    return max(abs(_p(r, f, n) - exact(x)) for x, n in pts) / ref


def _mic(r, name, f, col):
    """Waarde uit de microfoon-CSV (aster_acoustic_export) van deze run."""
    path = r.mess[:-len(".mess")] + "_microfoons.csv"
    rows = list(csv.reader(io.open(path, encoding="utf-8")))
    head = rows[0]
    for row in rows[1:]:
        if row[0] == name and abs(float(row[1]) - f) < 1e-6 * f:
            return float(row[head.index(col)])
    raise KeyError((name, f))


# ------------------------------------------------------------------
# Modellen
# ------------------------------------------------------------------
def build_tube(order, n, end, src="BRON_Zuiger", extra=(), size=None,
               mat=AIR, mics=False):
    """Buis langs x. end: "abs" (ABSORBEREND_Einde), "hard" (geen
    voorwaarde), ("Z", waarde) (Mecway impedance op x = L). src: naam van
    de vlakselectie op x = 0; extra: ruwe XML (lasten, #ASTER-tekst)."""
    lx, wy, wz = size or (LT, WT, WT)
    ny, nz = (1, 1) if size is None else (2, 2)

    def build():
        m = Model("A30", attrs={"maxmode": 9})
        m.elset("Lucht", "Lucht")
        h = hexmesh(m, "Lucht", n, ny, nz, box(0, lx, 0, wy, 0, wz),
                    order=order)
        m.mat("Lucht", *mat)
        m.face_selection(src, h["faces"]["u0"])
        if end == "abs":
            m.face_selection("ABSORBEREND_Einde", h["faces"]["u1"])
        elif isinstance(end, tuple):
            m.face_selection("Einde", h["faces"]["u1"])
            m.add('<impedance selection="Einde"><specificimpedancer>%r'
                  '</specificimpedancer><specificimpedancei>0'
                  '</specificimpedancei></impedance>' % end[1])
        step = lx / (n * order)
        m.node_selection("As", [m.node(i * step, 0, 0)
                                for i in range(n * order + 1)])
        if mics:
            m.node_selection("MICROFOON_Eind", [m.node(lx, 0, 0)])
            m.node_selection("MICROFOON_Midden", sorted(
                nid for nid, p in enumerate(m.nodes, 1)
                if abs(p[0] - 0.5 * lx) < 1e-9))
        for xml in extra:
            m.add(xml)
        return m
    return build


def _settings(freqs, **kw):
    s = {"harm_freq_list": ", ".join("%r" % f for f in freqs)}
    s.update(kw)
    return s


def _plane_checks(freqs, amp, tol, phase0=0.0, src="rho c v0"):
    """|P| constant = amp langs de as en faseval k L bij elke f."""
    out = []
    for f in freqs:
        k = 2.0 * math.pi * f / C
        lam = C / f
        out += [Check("f=%g Hz: max|p| (lambda = %s m)" % (f, "%g" % lam),
                      lambda r, f=f: _amp_range(r, f, max), ref=amp, tol=tol,
                      unit="Pa", source="vlakke golf, absorberend einde: "
                      "|P| = %s" % src),
                Check("f=%g Hz: min|p|" % f,
                      lambda r, f=f: _amp_range(r, f, min), ref=amp, tol=tol,
                      unit="Pa", source="vlakke golf: |P| constant"),
                Check("f=%g Hz: faseval over L" % f,
                      lambda r, f=f: _phase_lag(r, f), ref=k * LT, tol=tol,
                      unit="rad", source="P ~ e^(-ikx): fase -k x")]
    return out


ABS_F = (100.0, 340.0, 680.0)
TESTS = [
    Case("AH-01", GROUP, "Vlakke golf, absorberend einde, hex20 20x1x1",
         "acoustic_harm_3d", build_tube(2, 20, "abs"),
         _plane_checks(ABS_F, RHO * C * V0, 0.01) +
         [Check("f=340 Hz: spl op x = 0", lambda r: r.at(
             "spl", (0, 0, 0), float(_k_of(r, 340.0))),
             ref=_spl(RHO * C * V0), tol_abs=0.01, unit="dB",
             source="SPL = 20 log10(|P| / (sqrt(2) 20 uPa))"),
          Check("f=340 Hz: pressuremag = |pressure + i pressurei|",
                lambda r: r.at("pressuremag", (LT, 0, 0),
                               float(_k_of(r, 340.0))),
                ref=RHO * C * V0, tol=0.01, unit="Pa",
                source="vlakke golf")],
         settings=_settings(ABS_F),
         note="BRON_-vlak (source_velocity 1 mm/s), ABSORBEREND_-vlak; "
              "680 Hz = 10 hex20 per golflengte."),
    Case("AH-01b", GROUP, "Vlakke golf, absorberend einde, hex20 10x1x1 "
         "(5 per golflengte bij 680 Hz)", "acoustic_harm_3d",
         build_tube(2, 10, "abs"),
         _plane_checks(ABS_F, RHO * C * V0, 0.01),
         settings=_settings(ABS_F),
         note="grens van de WARNING voor kwadratische elementen."),
    Case("AH-02a", GROUP, "Vlakke golf, absorberend einde, hex8 20x1x1 "
         "(10 per golflengte bij 680 Hz)", "acoustic_harm_3d",
         build_tube(1, 20, "abs"),
         _plane_checks(ABS_F, RHO * C * V0, 0.02),
         settings=_settings(ABS_F),
         note="hex8: fasefout ~ (kh)^2/24 per golflengte; 2% grens."),
    Case("AH-02b", GROUP, "Vlakke golf, absorberend einde, hex8 40x1x1 "
         "(20 per golflengte bij 680 Hz)", "acoustic_harm_3d",
         build_tube(1, 40, "abs"),
         _plane_checks(ABS_F, RHO * C * V0, 0.01),
         settings=_settings(ABS_F)),
]


# ------------------------------------------------------------------
# AH-03 harde wand op x = L (resonanties bij n c / (2L) = 170, 340, 510)
# ------------------------------------------------------------------
HARD_F = (100.0, 250.0, 420.0)


def _hard(f):
    k = 2.0 * math.pi * f / C
    return lambda x: -1j * RHO * C * V0 * math.cos(k * (LT - x)) / \
        math.sin(k * LT)


def _hard_checks():
    out = []
    for f in HARD_F:
        k = 2.0 * math.pi * f / C
        out += [Check("f=%g Hz: |p(0)|" % f,
                      lambda r, f=f: abs(_p(r, f, _axis(r)[0][1])),
                      ref=abs(_hard(f)(0.0)), tol=0.01, unit="Pa",
                      source="|P| = rho c v0 |cos(k(L-x)) / sin(kL)|"),
                Check("f=%g Hz: |p(L)|" % f,
                      lambda r, f=f: abs(_p(r, f, _axis(r)[-1][1])),
                      ref=RHO * C * V0 / abs(math.sin(k * LT)), tol=0.01,
                      unit="Pa", source="|P(L)| = rho c v0 / |sin kL|"),
                Check("f=%g Hz: max|P - P_exact| / max|P_exact|" % f,
                      lambda r, f=f: _max_dev(r, f, _hard(f)), ref=0.0,
                      tol_abs=0.01, unit="-",
                      source="P = -i rho c v0 cos(k(L-x)) / sin(kL)")]
    return out


TESTS.append(Case(
    "AH-03", GROUP, "Buis met harde wand op x = L, hex20 20x1x1",
    "acoustic_harm_3d", build_tube(2, 20, "hard"), _hard_checks(),
    settings=_settings(HARD_F),
    note="Geen demping (WARNING); frequenties tussen de resonanties."))


# ------------------------------------------------------------------
# AH-04 / AH-05 ahlv100a nagebouwd: 1 x 0,1 x 0,2 m, hex20 15x2x2
# ------------------------------------------------------------------
AHLV_AIR = ['<mechanical type="Isotropic" speedofsound="343 m/s" />',
            u'<density type="Simple" density="1.3 kg/m³" />']
AHLV_SIZE = (1.0, 0.1, 0.2)
# (punt, VALE_REFE, VALE_CALC) van ahlv100a.comm (DYNA_VIBRA, PRES), met
# het teken van de bridge (-P, zie de kop)
AHLV_PTS = [
    ("NO1", (0.0, 0.0, 0.2), complex(-6.2426, 0.0),
     complex(-6.2425145509513, -3.1779206466043e-04)),
    ("NO780", (0.0, 0.05, 0.1), complex(-6.2426, 0.0),
     complex(-6.2425145732145, -3.1772413652398e-04)),
    ("NO751", (1.0, 0.0, 0.2), complex(6.02368, 1.6387),
     complex(6.0221726627638, 1.6440674349052)),
    ("NO763", (1.0, 0.05, 0.1), complex(6.02368, 1.6387),
     complex(6.0221725935735, 1.6440673833277)),
]
AHLV_DB = 109.88659136479       # PRAC_ELNO DB in NO1 (VALE_CALC)


def _ahlv_checks():
    out = []
    for name, pt, refe, calc in AHLV_PTS:
        def val(r, pt=pt, part="re"):
            p = -complex(r.at("pressure", pt), r.at("pressurei", pt))
            return p.real if part == "re" else p.imag
        out.append(Check("%s Re p (ahlv100a VALE_REFE)" % name,
                         lambda r, v=val: v(r), ref=refe.real, tol=1e-3,
                         unit="Pa", source="ahlv100a TEST_RESU, analytisch "
                         "(PRECISION 1e-3)"))
        out.append(Check("%s Re p (ahlv100a VALE_CALC)" % name,
                         lambda r, v=val: v(r), ref=calc.real, tol=1e-6,
                         unit="Pa", in_series=False,
                         source="ahlv100a TEST_RESU, Code-Aster zelf"))
        if refe.imag == 0.0:
            out.append(Check("%s Im p (ahlv100a VALE_REFE = 0)" % name,
                             lambda r, v=val: v(r, part="im"), ref=0.0,
                             tol_abs=1e-3 * abs(refe), unit="Pa",
                             source="ahlv100a TEST_RESU (CRITERE ABSOLU)"))
        else:
            out.append(Check("%s Im p (ahlv100a VALE_REFE)" % name,
                             lambda r, v=val: v(r, part="im"),
                             ref=refe.imag, tol=3.3e-3, unit="Pa",
                             source="ahlv100a TEST_RESU, analytisch "
                             "(PRECISION 3,3e-3)"))
        # absoluut: Im p in NO1/NO780 is ~3e-4 Pa (numerieke fasefout);
        # 1e-6 van |p| is daar 6e-6 Pa
        out.append(Check("%s Im p (ahlv100a VALE_CALC)" % name,
                         lambda r, v=val: v(r, part="im"), ref=calc.imag,
                         tol_abs=1e-6 * abs(calc), unit="Pa",
                         in_series=False,
                         source="ahlv100a TEST_RESU, Code-Aster zelf"))
    out.append(Check("NO1 spl + 3,0103 dB = ahlv100a DB",
                     lambda r: r.at("spl", (0.0, 0.0, 0.2)) +
                     20.0 * math.log10(math.sqrt(2.0)), ref=AHLV_DB,
                     tol_abs=0.01, unit="dB",
                     source="ahlv100a PRAC_ELNO DB = 20 log10(|P|/20 uPa)"))
    return out


# ahlv100a CALC_MODES (harde wanden, BANDE 1..1000 Hz): VALE_CALC en
# VALE_REFE van mode 1..8; bij de bridge is mode 1 de constante druk.
AHLV_MODES = [(171.50023005253, 171.5), (343.00726352824, 343.0),
              (514.55456202514, 514.5), (686.22640369374, 686.0),
              (858.17772650315, 857.5), (860.7191717656, 857.5),
              (877.64006231519, 874.482), (926.56793276098, 923.556)]

TESTS += [
    Case("AH-04", GROUP, "ahlv100a harmonisch nagebouwd (500 Hz), hex20 "
         "15x2x2", "acoustic_harm_3d",
         build_tube(2, 15, "abs", src="BRON_Entree", size=AHLV_SIZE,
                    mat=AHLV_AIR),
         _ahlv_checks(),
         settings=_settings((500.0,), source_velocity=0.014),
         note="Zelfde mesh, materiaal en bron als ahlv100a; teken "
              "omgekeerd (VNOR de fluid in)."),
    Case("AH-05", GROUP, "ahlv100a holtemodes nagebouwd, hex20 15x2x2",
         "acoustic_modes_3d",
         build_tube(2, 15, "hard", src="Entree", size=AHLV_SIZE,
                    mat=AHLV_AIR),
         [Check("f%d (ahlv100a mode %d, VALE_CALC)" % (k + 2, k + 1),
                lambda r, k=k: r.mode(k + 2), ref=calc, tol=1e-5, unit="Hz",
                in_series=False, source="ahlv100a TEST_RESU, Code-Aster "
                "zelf (zelfde mesh)")
          for k, (calc, _refe) in enumerate(AHLV_MODES)] +
         [Check("f%d (ahlv100a mode %d, VALE_REFE)" % (k + 2, k + 1),
                lambda r, k=k: r.mode(k + 2), ref=refe, tol=0.01, unit="Hz",
                source="ahlv100a TEST_RESU, analytisch (c/(2L) n, "
                "dwarsmodes)")
          for k, (_calc, refe) in enumerate(AHLV_MODES)],
         note="De CALC_MODES-stap van ahlv100a (harde wanden); de "
              "bridge-mode 1 is de constante druk (0 Hz)."),
]


# ------------------------------------------------------------------
# AH-06 impedantie Z = 2 rho c op x = L (Mecway impedance, reeel)
# ------------------------------------------------------------------
Z_F = (100.0, 250.0)
Z = 2.0 * RHO * C


def _imp(f):
    k = 2.0 * math.pi * f / C
    R = (Z - RHO * C) / (Z + RHO * C)
    A = RHO * C * V0 / (1.0 - R * cmath.exp(-2j * k * LT))
    return lambda x: A * (cmath.exp(-1j * k * x) +
                          R * cmath.exp(-1j * k * (2 * LT - x)))


def _imp_checks():
    out = []
    for f in Z_F:
        out += [Check("f=%g Hz: |p(0)|" % f,
                      lambda r, f=f: abs(_p(r, f, _axis(r)[0][1])),
                      ref=abs(_imp(f)(0.0)), tol=0.01, unit="Pa",
                      source="reflectie R = (Z - rho c)/(Z + rho c) = 1/3"),
                Check("f=%g Hz: |p(L)|" % f,
                      lambda r, f=f: abs(_p(r, f, _axis(r)[-1][1])),
                      ref=abs(_imp(f)(LT)), tol=0.01, unit="Pa",
                      source="idem"),
                Check("f=%g Hz: max|P - P_exact| / max|P_exact|" % f,
                      lambda r, f=f: _max_dev(r, f, _imp(f)), ref=0.0,
                      tol_abs=0.01, unit="-",
                      source="P = A (e^(-ikx) + R e^(-ik(2L-x)))")]
    out.append(Check("max|p| / min|p| langs de as (250 Hz, staande golf)",
                     lambda r: _amp_range(r, 250.0, max) /
                     _amp_range(r, 250.0, min), ref=2.0, tol=0.01, unit="-",
                     source="(1 + |R|) / (1 - |R|) = 2 (de as beslaat meer "
                            "dan een halve golflengte, 0,68 m)"))
    return out


TESTS.append(Case(
    "AH-06", GROUP, "Buis met impedantie Z = 2 rho c op x = L (Mecway "
    "impedance), hex20 20x1x1", "acoustic_harm_3d",
    build_tube(2, 20, ("Z", Z)), _imp_checks(), settings=_settings(Z_F),
    note="Z = 816 Pa.s/m reeel -> 3D_ABSO met c = Z / rho."))


# ------------------------------------------------------------------
# AH-07 opgelegde druk (fluid pressure 2 Pa op x = 0), absorberend einde
# AH-08 #ASTER-BRON met fase + microfoons
# ------------------------------------------------------------------
P0_F = (170.0, 500.0)
B_F = (200.0, 500.0)
V_B = 2.0e-3
TESTS += [
    Case("AH-07", GROUP, "Opgelegde druk 2 Pa op x = 0 (fluid pressure), "
         "absorberend einde, hex20 20x1x1", "acoustic_harm_3d",
         build_tube(2, 20, "abs", src="Speaker", extra=(
             '<fluidpressure selection="Speaker"><value>2 Pa</value>'
             '</fluidpressure>',)),
         _plane_checks(P0_F, 2.0, 0.01, src="p0"),
         settings=_settings(P0_F),
         note="fluid pressure != 0 is in de harmonische analyse een "
              "opgelegde druk (bron)."),
    Case("AH-08", GROUP, "#ASTER-BRON 2 mm/s met fase 90 graden, "
         "microfoons, hex20 20x1x1", "acoustic_harm_3d",
         build_tube(2, 20, "abs", src="Luidspreker", mics=True, extra=(
             "<ccxmodeldefinition>** ASTER-BRON Luidspreker = 2 mm/s, 90\n"
             "** ASTER-OPTIE harm_freq_list = 200, 500"
             "</ccxmodeldefinition>",)),
         [Check("f=%g Hz: Im p(0) = rho c v" % f,
                lambda r, f=f: _p(r, f, _axis(r)[0][1]).imag,
                ref=RHO * C * V_B, tol=0.01, unit="Pa",
                source="v = 2 mm/s e^(i 90 graden): P(0) = i rho c v")
          for f in B_F] +
         [Check("f=%g Hz: Re p(0) = 0" % f,
                lambda r, f=f: _p(r, f, _axis(r)[0][1]).real, ref=0.0,
                tol_abs=0.01 * RHO * C * V_B, unit="Pa", source="idem")
          for f in B_F] +
         [Check("f=%g Hz: microfoon Eind SPL" % f,
                lambda r, f=f: _mic(r, "MICROFOON_Eind", f, "spl_dB"),
                ref=_spl(RHO * C * V_B), tol_abs=0.05, unit="dB",
                source="vlakke golf |P| = rho c v; SPL op de effectieve "
                       "waarde")
          for f in B_F] +
         [Check("f=%g Hz: microfoon Midden SPL (8 knopen)" % f,
                lambda r, f=f: _mic(r, "MICROFOON_Midden", f, "spl_dB"),
                ref=_spl(RHO * C * V_B), tol_abs=0.05, unit="dB",
                source="idem (kwadratisch gemiddeld over de doorsnede)")
          for f in B_F] +
         [Check("f=%g Hz: microfoon Eind fase" % f,
                lambda r, f=f: (_mic(r, "MICROFOON_Eind", f, "fase_graden")
                                - (90.0 - math.degrees(2 * math.pi * f / C *
                                                       LT)) + 180.0) % 360.0
                - 180.0, ref=0.0, tol_abs=1.0, unit="graden",
                source="fase = 90 - k L (graden)")
          for f in B_F],
         note="#ASTER-BRON op een vlakselectie zonder voorvoegsel, "
              "frequenties met #ASTER-OPTIE; MICROFOON_-knoopselecties -> "
              "CSV."),
]
