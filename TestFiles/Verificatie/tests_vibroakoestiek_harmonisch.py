# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Vibro-akoestiek harmonisch: respons per frequentie van een constructie en
een luchtkolom (Mecway Frequency-model V30 met een fluidum-onderdeel ->
3D_FLUIDE + FLUI_STRU + 3D_FLUI_ABSO, FORMULATION U_P_PHI, DYNA_VIBRA HARM
PHYS; bridge-analyse vibroacoustic_harm_3d, bridge_vibroacoustic_harm.py).
Toelichting: docs/Verificatie.md par. 7 (VH) en docs/Akoestiek.md par. 7.

Referenties (vlakke golven in een kolom langs x met harde zijwanden,
p(t) = Re(P e^(i w t)), k = w / c; 1D Helmholtz P'' + k^2 P = 0 met
i w rho v = -dP/dx; zelf afgeleid, standaardresultaat o.a. Kinsler, Frey,
Coppens & Sanders, Fundamentals of Acoustics, niet lokaal ingezien):
  - stijve zuiger (massa m, veer k, kracht F) op een kolom met een
    rho-c-einde: de kolom drukt terug met p = rho c v (een golf die
    wegloopt), dus een demper rho c A:
        U = F / (k - m w^2 + i w rho c A),  P(x) = i w rho c U e^(-ikx);
  - idem met een impedantie Z op x = L: de kolom heeft aan de zuiger de
    ingangsimpedantie Z_in = rho c (1 + R e^(-2ikL)) / (1 - R e^(-2ikL)),
    R = (Z - rho c) / (Z + rho c):  U = F / (k - m w^2 + i w A Z_in);
  - met Rayleigh-demping C = alpha M op de zuiger (alleen de constructie):
    U = F / (k - m w^2 + i w (alpha m + rho c A));
  - invallende golf met amplitude p_in op x = 0 (INVALLEND_, Code-Aster
    ONDE_FLUI met PRES = 2 p_in) en een harde wand op x = L: P(x) = p_in
    (e^(-ikx) + e^(-ik(2L - x))), |P(L)| = 2 p_in;
  - Code-Aster-testcases ahlv101a (golfgeleider 1 x 0,1 x 0,2 m, hex20
    15x2x2, rho 1,3, c 343, ONDE_FLUI PRES 1 op de inlaat, 3D_FLUI_ABSO op
    de uitlaat, 500 Hz: P = 0,5 e^(-ikx)) en ahlv100b (dezelfde buis,
    VITE_FACE VNOR 0,014 m/s: P = -rho c v, VNOR naar buiten; de bridge
    legt de snelheid de fluid in op: P_bridge = -P_ahlv100b), TEST_RESU
    (VALE_REFE analytisch, VALE_CALC van Code-Aster zelf), nagebouwd op
    hetzelfde 3D-mesh als AH-04 (tests_akoestiek_harmonisch.py).
Geluidsdrukniveau: SPL = 20 log10(|P| / (sqrt(2) 20 uPa)) (zoals S3a).
"""
import io
import csv
import math
import cmath
from verif_lib import Case, Check, Model, hexmesh, quadmesh, box

GROUP = "Vibro-akoestiek harmonisch"
P_REF = 2.0e-5


def _spl(amp):
    return 20.0 * math.log10(amp / math.sqrt(2.0) / P_REF)


def _k_of(r, f):
    """mode-id van frequentie f (de stappen hebben t = mode-id)."""
    return min(r.modes, key=lambda k: abs(r.modes[k] - f))


def _p(r, f, point):
    st = r.step(float(_k_of(r, f)))
    nid = r.nearest(point)
    return complex(st["pressure"][nid], st["pressurei"][nid])


def _u(r, f, point):
    st = r.step(float(_k_of(r, f)))
    nid = r.nearest(point)
    return complex(st["displx"][nid], st["displix"][nid])


def _csv(r, suffix, name, f, col):
    """Waarde uit een CSV van deze run (microfoons of meetpunten)."""
    path = r.mess[:-len(".mess")] + suffix
    rows = list(csv.reader(io.open(path, encoding="utf-8")))
    head = rows[0]
    for row in rows[1:]:
        if row[0] == name and abs(float(row[1]) - f) < 1e-6 * f:
            return float(row[head.index(col)])
    raise KeyError((name, f))


def _rel(a, b):
    return abs(a - b) / abs(b)


def _phase_dev(a, b):
    """Faseverschil a t.o.v. b [graden], -180..180."""
    return (math.degrees(cmath.phase(a / b)) + 180.0) % 360.0 - 180.0


def _settings(freqs, **kw):
    s = {"harm_freq_list": ", ".join("%r" % f for f in freqs)}
    s.update(kw)
    return s


# ------------------------------------------------------------------
# VH-01 / VH-01b: ahlv101a en ahlv100b nagebouwd (alleen fluidum)
# ------------------------------------------------------------------
AHLV_AIR = ['<mechanical type="Isotropic" speedofsound="343 m/s" />',
            u'<density type="Simple" density="1.3 kg/m³" />']
AHLV_SIZE = (1.0, 0.1, 0.2)
# ahlv101a: (knoop, punt, VALE_REFE, VALE_CALC) van DYNA_VIBRA DEPL PRES
AHLV101_PTS = [
    ("NO1", (0.0, 0.0, 0.2), complex(0.5, 0.0),
     complex(0.49999657918807, 1.2726248366786e-05)),
    ("NO780", (0.0, 0.05, 0.1), complex(0.5, 0.0),
     complex(0.49999657781829, 1.2723087612133e-05)),
    ("NO751", (1.0, 0.0, 0.2), complex(-0.4825, -0.1313),
     complex(-0.48235155981834, -0.13166993000291)),
    ("NO763", (1.0, 0.05, 0.1), complex(-0.4825, -0.1313),
     complex(-0.48235155427641, -0.13166992587189)),
]
# ahlv100b (VNOR 0,014 naar buiten): met het teken van de bridge (-P)
AHLV100B_PTS = [
    ("NO1", (0.0, 0.0, 0.2), complex(6.2426, 0.0),
     -complex(-6.2425145509512, -3.177920646438e-04)),
    ("NO780", (0.0, 0.05, 0.1), complex(6.2426, 0.0),
     -complex(-6.2425145732145, -3.1772413645936e-04)),
    ("NO751", (1.0, 0.0, 0.2), -cmath.rect(6.2426, math.radians(15.2186)),
     -complex(6.0221726627638, 1.6440674349053)),
    ("NO763", (1.0, 0.05, 0.1), -cmath.rect(6.2426, math.radians(15.2186)),
     -complex(6.0221725935735, 1.6440673833277)),
]


def build_ahlv(inlet, extra=()):
    """Golfgeleider 1 x 0,1 x 0,2 m, hex20 15x2x2 (mesh van ahlv100a),
    alleen lucht; inlet = naam van de vlakselectie op x = 0, uitlaat
    ABSORBEREND_Uitlaat op x = 1; extra: ruwe XML (#ASTER-tekst)."""
    def build():
        m = Model("V30", attrs={"maxmode": 4})
        m.elset("Lucht", "Lucht")
        lx, wy, wz = AHLV_SIZE
        h = hexmesh(m, "Lucht", 15, 2, 2, box(0, lx, 0, wy, 0, wz), order=2)
        m.mat("Lucht", *AHLV_AIR)
        m.face_selection(inlet, h["faces"]["u0"])
        m.face_selection("ABSORBEREND_Uitlaat", h["faces"]["u1"])
        m.node_selection("MICROFOON_Uitlaat", [m.node(lx, 0.0, wz)])
        for xml in extra:
            m.add(xml)
        return m
    return build


def _ahlv_checks(points, prec_refe, what):
    out = []
    for name, pt, refe, calc in points:
        out.append(Check("%s |p - p_refe| / |p_refe| (%s VALE_REFE)" % (
            name, what), lambda r, pt=pt, refe=refe: _rel(_p(r, 500.0, pt),
                                                         refe),
            ref=0.0, tol_abs=prec_refe, unit="-",
            source="%s TEST_RESU, analytisch (PRECISION %g)" % (
                what, prec_refe)))
        out.append(Check("%s |p - p_calc| / |p_calc| (%s VALE_CALC)" % (
            name, what), lambda r, pt=pt, calc=calc: _rel(_p(r, 500.0, pt),
                                                         calc),
            ref=0.0, tol_abs=1e-6, unit="-", in_series=False,
            source="%s TEST_RESU, Code-Aster zelf" % what))
    return out


TESTS = [
    Case("VH-01", GROUP, "ahlv101a nagebouwd: invallende golf (INVALLEND_, "
         "ONDE_FLUI) en absorberende uitlaat, 500 Hz, hex20 15x2x2",
         "vibroacoustic_harm_3d", build_ahlv("INVALLEND_Inlaat"),
         _ahlv_checks(AHLV101_PTS, 1e-3, "ahlv101a") +
         [Check("microfoon uitlaat: SPL", lambda r: _csv(
             r, "_microfoons.csv", "MICROFOON_Uitlaat", 500.0, "spl_dB"),
             ref=_spl(0.5), tol_abs=0.01, unit="dB",
             source="|P| = 0,5 Pa: SPL = 20 log10(0,5 / (sqrt(2) 20 uPa))")],
         settings=_settings([500.0], incident_pressure=0.5),
         timeout=600,
         note="incident_pressure 0,5 Pa = ONDE_FLUI PRES 1 (zoals ahlv101a)."),
    Case("VH-01b", GROUP, "ahlv100b nagebouwd: snelheidsbron (BRON_, "
         "VITE_FACE) en absorberende uitlaat, 500 Hz, hex20 15x2x2",
         "vibroacoustic_harm_3d", build_ahlv("BRON_Inlaat"),
         _ahlv_checks(AHLV100B_PTS, 1e-3, "ahlv100b"),
         settings=_settings([500.0], source_velocity=0.014),
         timeout=600,
         note="source_velocity 0,014 m/s de fluid in; ahlv100b VNOR naar "
              "buiten (teken omgekeerd)."),
]


def _rotated(points, deg):
    """VALE_CALC-punten met een extra fase (bron met fase)."""
    rot = cmath.rect(1.0, math.radians(deg))
    return [(n, pt, refe * rot, calc * rot) for n, pt, refe, calc in points]


TESTS += [
    Case("VH-07", GROUP, "#ASTER-BRON 14 mm/s met fase 90 graden op een "
         "vlakselectie zonder voorvoegsel (ahlv100b x i)",
         "vibroacoustic_harm_3d", build_ahlv("Inlaat", extra=(
             "<ccxmodeldefinition>** ASTER-BRON Inlaat = 14 mm/s, 90\n"
             "** ASTER-OPTIE harm_freq_list = 500</ccxmodeldefinition>",)),
         [c for c in _ahlv_checks(_rotated(AHLV100B_PTS, 90.0), 1e-3,
                                  "ahlv100b") if "VALE_CALC" in c.name] +
         [Check("microfoon uitlaat: fase (CSV)", lambda r: _csv(
             r, "_microfoons.csv", "MICROFOON_Uitlaat", 500.0,
             "fase_graden"), ref=math.degrees(cmath.phase(
                 AHLV100B_PTS[2][3] * 1j)), tol_abs=0.01, unit="graden",
             source="ahlv100b VALE_CALC in NO751, 90 graden gedraaid")],
         settings={"harm_freq_list": "100"}, timeout=600,
         note="De frequentie komt uit #ASTER-OPTIE (instelling 100 Hz)."),
    Case("VH-08", GROUP, "#ASTER-BRON 0,5 Pa met fase -30 graden: "
         "invallende golf (ahlv101a x e^(-i 30))", "vibroacoustic_harm_3d",
         build_ahlv("Inlaat", extra=(
             "<ccxmodeldefinition>** ASTER-BRON Inlaat = 0.5 Pa, -30"
             "</ccxmodeldefinition>",)),
         [c for c in _ahlv_checks(_rotated(AHLV101_PTS, -30.0), 1e-3,
                                  "ahlv101a") if "VALE_CALC" in c.name],
         settings=_settings([500.0]), timeout=600,
         note="Een druk-eenheid in #ASTER-BRON maakt een invallende golf."),
]


# ------------------------------------------------------------------
# VH-02..04 Stijve zuiger (massa m, veer k) op een luchtkolom 1 x 0,1 x
# 0,1 m (hex20 10x1x1), kracht F0 in x op de zuiger.
# ------------------------------------------------------------------
RHO, C = 1.2, 340.0
AIR = ['<mechanical type="Isotropic" speedofsound="340 m/s" />',
       u'<density type="Simple" density="1.2 kg/m³" />']
LC, WC = 1.0, 0.1                   # kolom
AREA = WC * WC
TP = 0.01                           # zuigerdikte
RHO_P = 500.0                       # m = 0,05 kg
MP = RHO_P * AREA * TP
KP = 3000.0                         # N/m: f0 = 38,98 Hz
F0 = 1.0
PISTON_NODE = (-TP, 0.0, 0.0)
MIC_X = 0.5
PISTON_F = (5.0, 20.0, 35.0, 39.0, 45.0, 60.0, 100.0, 200.0)


def _zin(f, z_end=None):
    """Ingangsimpedantie van de kolom aan de zuiger (Pa.s/m)."""
    if z_end is None:
        return RHO * C
    k = 2.0 * math.pi * f / C
    refl = (z_end - RHO * C) / (z_end + RHO * C)
    e = cmath.exp(-2j * k * LC)
    return RHO * C * (1.0 + refl * e) / (1.0 - refl * e)


def piston_u(f, z_end=None, c_s=0.0):
    """U van de zuiger; c_s = demping van de constructie [N.s/m]."""
    w = 2.0 * math.pi * f
    return F0 / (KP - MP * w * w + 1j * w * (c_s + AREA * _zin(f, z_end)))


def piston_p(f, x, z_end=None, c_s=0.0):
    """P(x) in de kolom (rho-c-einde: wegloper; met Z: heen + terug)."""
    w = 2.0 * math.pi * f
    k = w / C
    v = 1j * w * piston_u(f, z_end, c_s)
    if z_end is None:
        return RHO * C * v * cmath.exp(-1j * k * x)
    refl = (z_end - RHO * C) / (z_end + RHO * C)
    e = cmath.exp(-2j * k * LC)
    a = RHO * C * v / (1.0 - refl * e)
    return a * (cmath.exp(-1j * k * x) + refl * cmath.exp(-1j * k *
                                                          (2 * LC - x)))


def build_piston(end="abs", alpha=0.0, beta=0.0, damper=0.0, force=None,
                 extra=()):
    def build():
        m = Model("V30", attrs={"maxmode": 4})
        m.elset("Lucht", "Lucht")
        m.elset("Zuiger", "Zuiger")
        m.elset("Veer", "Veer")
        h = hexmesh(m, "Lucht", 10, 1, 1, box(0, LC, 0, WC, 0, WC), order=2)
        hexmesh(m, "Zuiger", 1, 1, 1, box(-TP, 0, 0, WC, 0, WC), order=2)
        m.mat("Lucht", *AIR)
        m.mat("Zuiger", '<mechanical type="Isotropic" youngsmodulus='
              '"200000 MPa" poissonratio="0" />',
              u'<density type="Simple" density="%r kg/m³" />' % RHO_P)
        m.mat("Veer", '<mechanical type="Spring" springconstant="%r N/m" />'
              % KP)
        g = m.node(-TP - 0.5, 0.0, 0.0)
        n = m.node(*PISTON_NODE)
        m.elem("Veer", "line2", [g, n])
        m.node_selection("Grond", [g])
        m.add('<fixed selection="Grond" />')
        piston = [i + 1 for i, p in enumerate(m.nodes)
                  if -TP - 1e-9 <= p[0] <= 1e-9]
        m.node_selection("Zuigerknopen", piston)
        for d in ((0, 1, 0), (0, 0, 1)):
            m.add('<displacement selection="Zuigerknopen"><value>0 m'
                  '</value><x>%d</x><y>%d</y><z>%d</z></displacement>' % d)
        m.node_selection("Last", [n])
        m.add('<force selection="Last">%s<y>0 N</y><z>0 N</z></force>' % (
            force or "<x>%r N</x>" % F0))
        for xml in extra:
            m.add(xml)
        m.node_selection("MEETPUNT_Zuiger", [n])
        m.node_selection("MICROFOON_Midden", sorted(
            nid for nid, p in enumerate(m.nodes, 1)
            if abs(p[0] - MIC_X) < 1e-9))
        if end == "abs":
            m.face_selection("ABSORBEREND_Einde", h["faces"]["u1"])
        elif end is not None:
            m.face_selection("Einde", h["faces"]["u1"])
            m.add('<impedance selection="Einde"><specificimpedancer>%r'
                  '</specificimpedancer><specificimpedancei>0'
                  '</specificimpedancei></impedance>' % end)
        if alpha or beta:
            m.add('<rayleighdamping selection=""><alpha>%r rad/s</alpha>'
                  '<beta>%r s</beta></rayleighdamping>' % (alpha, beta))
        if damper:
            # demper parallel aan de veer (zelfde knopen, langs x)
            m.elset("Demper", "Demper")
            m.mat("Demper", '<mechanical type="Damper" dampingcoefficient='
                  '"%r N.s/m" />' % damper)
            m.elem("Demper", "line2", [g, n])
        return m
    return build


def _piston_checks(z_end=None, c_s=0.0, what="rho-c-einde",
                   freqs=PISTON_F, mics=True):
    src = ("U = F / (k - m w^2 + i w %s), stijve zuiger, %s" % (
        "(c_s + A Z_in)" if c_s else "A Z_in", what))
    out = []
    for f in freqs:
        ue = piston_u(f, z_end, c_s)
        out += [
            Check("f=%g Hz: |U - U_exact| / |U_exact| zuiger" % f,
                  lambda r, f=f, ue=ue: _rel(_u(r, f, PISTON_NODE), ue),
                  ref=0.0, tol_abs=0.01, unit="-", source=src),
            Check("f=%g Hz: fase U zuiger" % f,
                  lambda r, f=f, ue=ue: _phase_dev(_u(r, f, PISTON_NODE), ue),
                  ref=0.0, tol_abs=0.5, unit="graden", source=src)]
        pe = piston_p(f, MIC_X, z_end, c_s)
        out.append(Check(
            "f=%g Hz: |P - P_exact| / |P_exact| op x = %g" % (f, MIC_X),
            lambda r, f=f, pe=pe: _rel(_p(r, f, (MIC_X, 0.0, 0.0)), pe),
            ref=0.0, tol_abs=0.01, unit="-",
            source="P(x) uit de zuigersnelheid i w U (%s)" % what))
    if mics:
        f = 39.0
        pe = piston_p(f, MIC_X, z_end, c_s)
        out += [
            Check("microfoon midden 39 Hz: SPL (CSV)", lambda r: _csv(
                r, "_microfoons.csv", "MICROFOON_Midden", 39.0, "spl_dB"),
                ref=_spl(abs(pe)), tol_abs=0.09, unit="dB",
                source="SPL van |P_exact| (0,09 dB = 1% in |P|)"),
            Check("meetpunt zuiger 39 Hz: |U| (CSV)", lambda r: _csv(
                r, "_meetpunten.csv", "MEETPUNT_Zuiger", 39.0, "u_amp_m"),
                ref=abs(piston_u(39.0, z_end, c_s)), tol=0.01, unit="m",
                source=src)]
    return out


ZETA_RAD = RHO * C * AREA / (2.0 * math.sqrt(KP * MP))
TAB_F = [(0.0, 1.0), (200.0, 3.0)]          # F(f) = 1 + f / 100 [N]


def _tab_f(f):
    return TAB_F[0][1] + (TAB_F[1][1] - TAB_F[0][1]) * f / TAB_F[1][0]


TESTS += [
    Case("VH-02", GROUP, "Zuiger (massa-veer) op een luchtkolom met "
         "rho-c-einde: stralingsdemping", "vibroacoustic_harm_3d",
         build_piston("abs"), _piston_checks(),
         settings=_settings(PISTON_F, fields=["displacements", "stresses",
                                              "vonmises"]),
         timeout=600,
         note="m = %.3g kg, k = %g N/m (f0 = %.2f Hz), rho c A = %.3g N.s/m "
              "(zeta = %.3f), F0 = 1 N." % (
                  MP, KP, math.sqrt(KP / MP) / (2 * math.pi),
                  RHO * C * AREA, ZETA_RAD)),
    Case("VH-03", GROUP, "Zuiger op een luchtkolom met impedantie Z = 2 rho "
         "c op het einde (Mecway impedance)", "vibroacoustic_harm_3d",
         build_piston(2.0 * RHO * C),
         _piston_checks(z_end=2.0 * RHO * C, what="Z = 2 rho c op x = L",
                        mics=False),
         settings=_settings(PISTON_F), timeout=600,
         note="Z_in = rho c (1 + R e^(-2ikL)) / (1 - R e^(-2ikL)), R = 1/3."),
    Case("VH-04", GROUP, "Zuiger op een luchtkolom met rho-c-einde en "
         "Rayleigh-demping alpha op de constructie", "vibroacoustic_harm_3d",
         build_piston("abs", alpha=200.0),
         _piston_checks(c_s=200.0 * MP, what="rho-c-einde, C = alpha M op "
                        "de zuiger", mics=False),
         settings=_settings(PISTON_F), timeout=600,
         note="alpha = 200 1/s: alpha m = 10 N.s/m naast rho c A = 4,08 "
              "N.s/m; de demping werkt alleen op de constructie."),
    Case("VH-02t", GROUP, "Zuiger op een luchtkolom met rho-c-einde: kracht "
         "uit een Mecway-tabel (amplitude tegen de frequentie) met fase 90 "
         "graden (#ASTER-FASE)", "vibroacoustic_harm_3d",
         build_piston("abs", force=''.join('<x time="%r s">%r N</x>' % p
                                           for p in TAB_F),
                      extra=("<ccxmodeldefinition>** ASTER-FASE Last = 90"
                             "</ccxmodeldefinition>",)),
         [Check("f=%g Hz: |U - U_exact| / |U_exact| zuiger" % f,
                lambda r, f=f: _rel(_u(r, f, PISTON_NODE),
                                    _tab_f(f) * 1j * piston_u(f)),
                ref=0.0, tol_abs=0.01, unit="-",
                source="U = F(f) e^(i 90) / (k - m w^2 + i w rho c A), "
                       "F(f) = 1 + f / 100 N (tabel)")
          for f in (5.0, 39.0, 100.0, 200.0)],
         settings=_settings((5.0, 39.0, 100.0, 200.0)), timeout=600,
         note="Tabel: kolom 'time' = frequentie [Hz] (DEFI_FONCTION "
              "NOM_PARA='FREQ' als FONC_MULT), zoals harmonic_3d (HR-04)."),
    Case("VH-04b", GROUP, "Zuiger op een luchtkolom met rho-c-einde en "
         "Rayleigh-demping beta (C = beta K: de veer)",
         "vibroacoustic_harm_3d", build_piston("abs", beta=2.0e-3),
         _piston_checks(c_s=2.0e-3 * KP, what="rho-c-einde, C = beta K (de "
                        "starre zuiger vervormt niet: c = beta k)",
                        mics=False),
         settings=_settings(PISTON_F), timeout=600,
         note="beta = 2e-3 s: beta k = 6 N.s/m."),
    Case("VH-04c", GROUP, "Zuiger op een luchtkolom met rho-c-einde en een "
         "Mecway-demper (Damper) parallel aan de veer",
         "vibroacoustic_harm_3d", build_piston("abs", damper=10.0),
         _piston_checks(c_s=10.0, what="rho-c-einde, demper c = 10 N.s/m",
                        mics=False),
         settings=_settings(PISTON_F), timeout=600,
         note="AMOR_MECA van de demper (DIS_T); alarm MODELISA8_13: veer en "
              "demper op dezelfde twee knopen (onschuldig)."),
]


# ------------------------------------------------------------------
# VH-05 invallende golf (INVALLEND_, p_in = 1 Pa) en een harde wand op
# x = L (lucht 1 x 0,1 x 0,1 m, hex20 20x1x1)
# ------------------------------------------------------------------
WAVE_F = (100.0, 250.0, 420.0)


def build_wave(open_end=False):
    def build():
        m = Model("V30", attrs={"maxmode": 4})
        m.elset("Lucht", "Lucht")
        h = hexmesh(m, "Lucht", 20, 1, 1, box(0, LC, 0, WC, 0, WC),
                    order=2)
        m.mat("Lucht", *AIR)
        m.face_selection("INVALLEND_Inlaat", h["faces"]["u0"])
        if open_end:
            m.face_selection("OPEN_Einde", h["faces"]["u1"])
        return m
    return build


def _wave(f, x, refl=1.0):
    k = 2.0 * math.pi * f / C
    return cmath.exp(-1j * k * x) + refl * cmath.exp(-1j * k * (2 * LC - x))


TESTS.append(Case(
    "VH-05", GROUP, "Invallende golf p_in = 1 Pa en een harde wand op x = L",
    "vibroacoustic_harm_3d", build_wave(),
    [Check("f=%g Hz: |P - P_exact| / max|P_exact| op x = %g" % (f, x),
           lambda r, f=f, x=x: abs(_p(r, f, (x, 0.0, 0.0)) - _wave(f, x)) /
           2.0, ref=0.0, tol_abs=0.01, unit="-",
           source="P = p_in (e^(-ikx) + e^(-ik(2L-x))): de inlaat laat de "
                  "teruggekaatste golf door")
     for f in WAVE_F for x in (0.0, 0.5, LC)],
    settings=_settings(WAVE_F), timeout=600,
    note="incident_pressure 1 Pa (standaard); geen andere demping dan de "
         "inlaat (ONDE_FLUI is ook absorberend)."))
TESTS.append(Case(
    "VH-05b", GROUP, "Invallende golf p_in = 1 Pa en p = 0 op x = L "
    "(OPEN_-vlakselectie)", "vibroacoustic_harm_3d", build_wave(True),
    [Check("f=%g Hz: |P - P_exact| / max|P_exact| op x = %g" % (f, x),
           lambda r, f=f, x=x: abs(_p(r, f, (x, 0.0, 0.0)) -
                                   _wave(f, x, -1.0)) / 2.0,
           ref=0.0, tol_abs=0.01, unit="-",
           source="P = p_in (e^(-ikx) - e^(-ik(2L-x))): drukknoop op x = L")
     for f in WAVE_F for x in (0.0, 0.5)],
    settings=_settings(WAVE_F), timeout=600,
    note="OPEN_ = DDL_IMPO PRES = PHI = 0, zoals in de modes (S3b)."))


# ------------------------------------------------------------------
# VH-06 Schaalzuiger (DKT quad4 2x2, t = 10 mm, rho 500: m = 0,05 kg) op
# x = 0 van dezelfde kolom met rho-c-einde, lucht hex8 20x2x2; veer k en
# kracht F0 in het midden van de plaat. De plaat is geleid (DY = DZ = 0,
# rotaties vast, zoals VA-03): een stijve zuiger, dus dezelfde formule als
# VH-02.
# ------------------------------------------------------------------
SHELL_NODE = (0.0, 0.5 * WC, 0.5 * WC)
SHELL_F = (5.0, 39.0, 100.0, 200.0)


def build_shell_piston():
    m = Model("V30", attrs={"maxmode": 4})
    m.elset("Plaat", "Plaat")
    m.elset("Lucht", "Lucht")
    m.elset("Veer", "Veer")
    quadmesh(m, "Plaat", 2, 2, lambda u, v: (0.0, WC * u, WC * v))
    h = hexmesh(m, "Lucht", 20, 2, 2, box(0, LC, 0, WC, 0, WC), order=1)
    m.mat("Plaat", '<geometric type="Plate" thickness="%r mm" />' % (
        TP * 1000.0),
        '<mechanical type="Isotropic" youngsmodulus="210000 MPa" '
        'poissonratio="0.3" />',
        u'<density type="Simple" density="%r kg/m³" />' % RHO_P)
    m.mat("Lucht", *AIR)
    m.mat("Veer", '<mechanical type="Spring" springconstant="%r N/m" />'
          % KP)
    g = m.node(-0.5, 0.5 * WC, 0.5 * WC)
    n = m.node(*SHELL_NODE)
    m.elem("Veer", "line2", [g, n])
    m.node_selection("Grond", [g])
    m.add('<fixed selection="Grond" />')
    plate = sorted(i + 1 for i, p in enumerate(m.nodes) if abs(p[0]) < 1e-9)
    m.node_selection("Plaatknopen", plate)
    for d in ((0, 1, 0), (0, 0, 1)):
        m.add('<displacement selection="Plaatknopen"><value>0 m</value>'
              '<x>%d</x><y>%d</y><z>%d</z></displacement>' % d)
    for d in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
        m.add('<noderotation selection="Plaatknopen"><value>0 rad</value>'
              '<x>%d</x><y>%d</y><z>%d</z></noderotation>' % d)
    m.node_selection("Last", [n])
    m.add('<force selection="Last"><x>%r N</x><y>0 N</y><z>0 N</z>'
          '</force>' % F0)
    m.face_selection("ABSORBEREND_Einde", h["faces"]["u1"])
    return m


def _shell_checks():
    src = "U = F / (k - m w^2 + i w rho c A), stijve schaalzuiger (DKT)"
    out = []
    for f in SHELL_F:
        ue = piston_u(f)
        out += [
            Check("f=%g Hz: |U - U_exact| / |U_exact| plaat" % f,
                  lambda r, f=f, ue=ue: _rel(_u(r, f, SHELL_NODE), ue),
                  ref=0.0, tol_abs=0.01, unit="-", source=src),
            Check("f=%g Hz: |P - P_exact| / |P_exact| op x = %g" % (f, MIC_X),
                  lambda r, f=f: _rel(_p(r, f, (MIC_X, 0.0, 0.0)),
                                      piston_p(f, MIC_X)),
                  ref=0.0, tol_abs=0.01, unit="-",
                  source="P(x) = i w rho c U e^(-ikx)")]
    return out


TESTS.append(Case(
    "VH-06", GROUP, "Schaalzuiger (DKT) met veer op een luchtkolom met "
    "rho-c-einde, hex8 20x2x2", "vibroacoustic_harm_3d", build_shell_piston,
    _shell_checks(), settings=_settings(SHELL_F), timeout=600,
    note="Schaal op een vlak van het fluidum: automatisch koppelvlak; "
         "spanningen alleen op volumes (hier geen)."))
