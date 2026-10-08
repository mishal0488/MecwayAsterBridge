# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Vibro-akoestiek: gekoppelde modes van een constructie en een lucht- of
vloeistofkolom (Mecway Frequency-model V30 met een fluidum-onderdeel ->
3D_FLUIDE + FLUI_STRU, FORMULATION U_P_PHI; bridge-analyse
vibroacoustic_modes_3d, bridge_vibroacoustic.py). Toelichting:
docs/Verificatie.md par. 7 (VA) en docs/Akoestiek.md par. 6.

Referenties (vlakke golven in een kolom met harde zijwanden, per
oppervlakte-eenheid; p(x) en de verplaatsing u volgen uit rho u'' = -dp/dx,
dus -rho w^2 u = -dp/dx; k = w / c; zelf afgeleid uit de 1D-golfvergelijking
met de randvoorwaarden, zoals Blevins, Formulas for Natural Frequency and
Mode Shape, en Kinsler e.a., Fundamentals of Acoustics; niet lokaal
ingezien):
  - stijve zuiger (massa m per oppervlak) tussen een waterkolom L1 met
    vrij oppervlak (p = 0) en een gesloten luchtkolom L2 (adlv100a):
      m + rho1 tan(k1 L1) / k1 - rho2 cot(k2 L2) / k2 = 0;
    benadering voor lage frequentie (massa m + rho1 L1 op de luchtveer
    rho2 c2^2 / L2): f = 1/(2 pi) sqrt(rho2 c2^2 / (L2 (m + rho1 L1))) =
    13,8287 Hz; adlv100a geeft VALE_REFE 13,8285 Hz (analytisch) en zelf
    13,827698 Hz; de exacte wortel is 13,82768 Hz. Druk/verplaatsing op
    het zuigervlak: lucht p/u = rho2 w^2 cot(k2 L2) / k2, water
    p/u = rho1 w^2 tan(k1 L1) / k1;
  - elastische staaf (lengte Ls, Zs = rho_s c_s, ingeklemd op x = 0,
    zijdelings geleid, nu = 0) met een waterkolom Lf die gesloten is:
      Zs cot(ks Ls) + Zf cot(kf Lf) = 0;
  - schaalzuiger (massa m = rho t per oppervlak) op een gesloten
    luchtkolom L: m k tan(k L) = rho;
  - buis van twee media (alleen fluidum, harde einden), druk en
    verplaatsing continu op het grensvlak:
      tan(k1 L1) / Z1 + tan(k2 L2) / Z2 = 0   (Z = rho c).
De wortels rekent _roots (bisectie op een vorm zonder polen).
"""
import math
from verif_lib import Case, Check, Model, hexmesh, quadmesh, box

GROUP = "Vibro-akoestiek"


def _roots(fun, fmax, df=0.05):
    """Nulpunten van fun(f) op (0, fmax] (tekenwissel + bisectie)."""
    out = []
    f, prev = df, fun(df)
    while f < fmax:
        f2 = f + df
        cur = fun(f2)
        if prev == 0.0 or prev * cur < 0.0:
            a, b = f, f2
            for _ in range(100):
                c = 0.5 * (a + b)
                if fun(a) * fun(c) <= 0.0:
                    b = c
                else:
                    a = c
            out.append(0.5 * (a + b))
        prev, f = cur, f2
    return out


def _mode_checks(refs, tol, source, first=1):
    return [Check("f%d" % k, lambda r, k=k: r.mode(k), ref=f, tol=tol,
                  unit="Hz", source=source)
            for k, f in enumerate(refs, start=first)]


def _ratio(r, mode, point):
    """pressure / displx in de knoop het dichtst bij point (mode)."""
    st = r.step(float(mode))
    nid = r.nearest(point)
    return st["pressure"][nid] / st["displx"][nid]


# ------------------------------------------------------------------
# VA-01 adlv100a nagebouwd: water x 0..0,05 m (vrij oppervlak op x = 0),
# stalen zuiger 0,05..0,075 m (zijdelings geleid), lucht 0,075..0,15 m
# (gesloten), doorsnede 0,05 x 0,05 m; hex20.
# ------------------------------------------------------------------
RHO1, C1, L1 = 1000.0, 1400.0, 0.05
RHO2, C2, L2 = 1.2, 340.0, 0.075
RHOS, TP = 7800.0, 0.025
WP = 0.05
MP = RHOS * TP                      # 195 kg/m2
WATER = ['<mechanical type="Isotropic" speedofsound="1400 m/s" />',
         u'<density type="Simple" density="1000 kg/m³" />']
AIR = ['<mechanical type="Isotropic" speedofsound="340 m/s" />',
       u'<density type="Simple" density="1.2 kg/m³" />']
STEEL = ['<mechanical type="Isotropic" youngsmodulus="200000 MPa" '
         'poissonratio="0.3" />',
         u'<density type="Simple" density="7800 kg/m³" />']


def _piston_eq(f):
    w = 2.0 * math.pi * f
    k1, k2 = w / C1, w / C2
    # (m + rho1 tan(k1L1)/k1 - rho2 cot(k2L2)/k2) x k1 k2 cos(k1L1) sin(k2L2)
    return (MP * k1 * k2 * math.cos(k1 * L1) * math.sin(k2 * L2) +
            RHO1 * k2 * math.sin(k1 * L1) * math.sin(k2 * L2) -
            RHO2 * k1 * math.cos(k1 * L1) * math.cos(k2 * L2))


PISTON = _roots(_piston_eq, 3000.0)
F_LUMPED = math.sqrt(RHO2 * C2 ** 2 / (L2 * (MP + RHO1 * L1))) / (2 * math.pi)
F_ADLV = 13.8285                    # adlv100a VALE_REFE (analytisch)


def build_va01(coupling=False):
    def build():
        m = Model("V30", attrs={"maxmode": 2})
        m.elset("Water", "Water")
        m.elset("Zuiger", "Staal")
        m.elset("Lucht", "Lucht")
        hw = hexmesh(m, "Water", 4, 1, 1, box(0, L1, 0, WP, 0, WP), order=2)
        hexmesh(m, "Zuiger", 1, 1, 1, box(L1, L1 + TP, 0, WP, 0, WP),
                order=2)
        ha = hexmesh(m, "Lucht", 6, 1, 1, box(L1 + TP, L1 + TP + L2, 0, WP,
                                              0, WP), order=2)
        m.mat("Water", *WATER)
        m.mat("Lucht", *AIR)
        m.mat("Staal", *STEEL)
        m.face_selection("OPEN_Vrij", hw["faces"]["u0"])
        piston = [i + 1 for i, p in enumerate(m.nodes)
                  if L1 - 1e-9 <= p[0] <= L1 + TP + 1e-9]
        m.node_selection("Zuigerknopen", piston)
        for d in ((0, 1, 0), (0, 0, 1)):
            m.add('<displacement selection="Zuigerknopen"><value>0 m</value>'
                  '<x>%d</x><y>%d</y><z>%d</z></displacement>' % d)
        if coupling:
            # handmatige koppelvlakken: de vlakken van het fluidum aan de
            # zuiger (water x = L1, lucht x = L1 + TP)
            m.face_selection("KOPPEL_Water", hw["faces"]["u1"])
            m.face_selection("KOPPEL_Lucht", ha["faces"]["u0"])
        return m
    return build


SRC_PISTON = ("zuiger tussen water (vrij oppervlak) en gesloten lucht: m + "
              "rho1 tan(k1L1)/k1 - rho2 cot(k2L2)/k2 = 0 (1D, exact)")
W1 = 2 * math.pi * PISTON[0]
TESTS = [
    Case("VA-01", GROUP, "adlv100a nagebouwd: zuiger + water + lucht, hex20",
         "vibroacoustic_modes_3d", build_va01(),
         [Check("f1 t.o.v. adlv100a", lambda r: r.mode(1), ref=F_ADLV,
                tol=0.005, unit="Hz",
                source="adlv100a VALE_REFE 13,8285 Hz (analytisch, "
                       "PRECISION 0,5%)"),
          Check("f1 exact", lambda r: r.mode(1), ref=PISTON[0], tol=0.001,
                unit="Hz", source=SRC_PISTON),
          Check("f2 exact", lambda r: r.mode(2), ref=PISTON[1], tol=0.005,
                unit="Hz", source=SRC_PISTON),
          Check("p/u lucht op de zuiger (mode 1)",
                lambda r: _ratio(r, 1, (L1 + TP, 0, 0)),
                ref=RHO2 * W1 ** 2 / math.tan(W1 / C2 * L2) / (W1 / C2),
                tol=0.005, unit="Pa/m", in_series=False,
                source="p = rho2 w^2 U cot(k2 L2) / k2 (gesloten kolom)"),
          Check("p/u water op de zuiger (mode 1)",
                lambda r: _ratio(r, 1, (L1, 0, 0)),
                ref=RHO1 * W1 ** 2 * math.tan(W1 / C1 * L1) / (W1 / C1),
                tol=0.005, unit="Pa/m", in_series=False,
                source="p = rho1 w^2 U tan(k1 L1) / k1 (vrij oppervlak)")],
         note="Lage-frequentiebenadering %.4f Hz; exacte wortel %.5f Hz." % (
             F_LUMPED, PISTON[0])),
    Case("VA-01k", GROUP, "adlv100a nagebouwd met KOPPEL_-vlakselecties",
         "vibroacoustic_modes_3d", build_va01(coupling=True),
         [Check("f1 exact", lambda r: r.mode(1), ref=PISTON[0], tol=0.001,
                unit="Hz", source=SRC_PISTON),
          Check("f2 exact", lambda r: r.mode(2), ref=PISTON[1], tol=0.005,
                unit="Hz", source=SRC_PISTON)],
         note="Zelfde model als VA-01, koppelvlakken handmatig (KOPPEL_)."),
]


# ------------------------------------------------------------------
# VA-02 Elastische staaf (PVC-achtig, E 3 GPa, rho 1400, nu 0) 0..0,5 m,
# ingeklemd op x = 0 en zijdelings geleid, met een gesloten waterkolom
# 0,5..1,5 m; doorsnede 0,1 x 0,1 m; hex20. Zs / Zf = 1,38: sterke
# koppeling.
# ------------------------------------------------------------------
ES, RHOS2, LS = 3.0e9, 1400.0, 0.5
RHOF, CF, LF = 1000.0, 1480.0, 1.0
CS = math.sqrt(ES / RHOS2)
ZS, ZF = RHOS2 * CS, RHOF * CF


def _rod_eq(f):
    w = 2.0 * math.pi * f
    a, b = w * LS / CS, w * LF / CF
    # Zs cot(a) + Zf cot(b) = 0, vermenigvuldigd met sin(a) sin(b)
    return ZS * math.cos(a) * math.sin(b) + ZF * math.sin(a) * math.cos(b)


ROD = _roots(_rod_eq, 2500.0)


def build_va02():
    m = Model("V30", attrs={"maxmode": 4})
    m.elset("Staaf", "PVC")
    m.elset("Water", "Water")
    hs = hexmesh(m, "Staaf", 10, 1, 1, box(0, LS, 0, 0.1, 0, 0.1), order=2)
    hexmesh(m, "Water", 20, 1, 1, box(LS, LS + LF, 0, 0.1, 0, 0.1), order=2)
    m.mat("PVC", '<mechanical type="Isotropic" youngsmodulus="3000 MPa" '
          'poissonratio="0" />',
          u'<density type="Simple" density="1400 kg/m³" />')
    m.mat("Water", '<mechanical type="Isotropic" speedofsound="1480 m/s" />',
          u'<density type="Simple" density="1000 kg/m³" />')
    m.node_selection("Inklemming", sorted(hs["nodes"]["u0"]))
    m.add('<fixed selection="Inklemming" />')
    rod = [i + 1 for i, p in enumerate(m.nodes) if p[0] <= LS + 1e-9]
    m.node_selection("Staafknopen", rod)
    for d in ((0, 1, 0), (0, 0, 1)):
        m.add('<displacement selection="Staafknopen"><value>0 m</value>'
              '<x>%d</x><y>%d</y><z>%d</z></displacement>' % d)
    return m


TESTS.append(Case(
    "VA-02", GROUP, "Elastische staaf + gesloten waterkolom, hex20 10+20",
    "vibroacoustic_modes_3d", build_va02,
    _mode_checks(ROD[:4], 0.01, "Zs cot(ks Ls) + Zf cot(kf Lf) = 0 "
                 "(staaf ingeklemd, kolom gesloten; 1D, exact)"),
    note="Staaf zijdelings geleid (DY = DZ = 0), nu = 0: zuiver 1D."))


# ------------------------------------------------------------------
# VA-03 Schaalzuiger (DKT quad4 2x2, t = 0,01 m, staal) op x = 0 van een
# gesloten luchtkolom 0..0,5 m, doorsnede 0,1 x 0,1 m; lucht hex8 50x2x2.
# De plaat is geleid (DY = DZ = 0) en zijn rotaties liggen vast: hij kan
# alleen als zuiger bewegen. Koppeling via een schaal op een fluidumvlak.
# ------------------------------------------------------------------
TS, LA = 0.01, 0.5
MS = 7800.0 * TS                    # 78 kg/m2


def _shell_eq(f):
    k = 2.0 * math.pi * f / C2
    # m k tan(kL) - rho = 0, vermenigvuldigd met cos(kL)
    return MS * k * math.sin(k * LA) - RHO2 * math.cos(k * LA)


SHELL = _roots(_shell_eq, 800.0)


def build_va03():
    m = Model("V30", attrs={"maxmode": 3})
    m.elset("Plaat", "Plaatstaal")
    m.elset("Lucht", "Lucht")
    q = quadmesh(m, "Plaat", 2, 2, lambda u, v: (0.0, 0.1 * u, 0.1 * v))
    hexmesh(m, "Lucht", 50, 2, 2, box(0, LA, 0, 0.1, 0, 0.1), order=1)
    m.mat("Plaatstaal", '<geometric type="Plate" thickness="10 mm" />',
          '<mechanical type="Isotropic" youngsmodulus="210000 MPa" '
          'poissonratio="0.3" />',
          u'<density type="Simple" density="7800 kg/m³" />')
    m.mat("Lucht", *AIR)
    nodes = sorted(set(n for side in q["nodes"].values() for n in side) |
                   set(i + 1 for i, p in enumerate(m.nodes)
                       if abs(p[0]) < 1e-9))
    m.node_selection("Plaatknopen", nodes)
    for d in ((0, 1, 0), (0, 0, 1)):
        m.add('<displacement selection="Plaatknopen"><value>0 m</value>'
              '<x>%d</x><y>%d</y><z>%d</z></displacement>' % d)
    # noderotation met een richting legt de rotatie om die as vast
    # (bridge: LIAISON_OBLIQUE DRX); daarom drie keer
    for d in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
        m.add('<noderotation selection="Plaatknopen"><value>0 rad</value>'
              '<x>%d</x><y>%d</y><z>%d</z></noderotation>' % d)
    return m


TESTS.append(Case(
    "VA-03", GROUP, "Schaalzuiger (DKT) op een gesloten luchtkolom, hex8 "
    "50x2x2", "vibroacoustic_modes_3d", build_va03,
    _mode_checks(SHELL[:3], 0.01, "m k tan(kL) = rho (zuiger op een "
                 "gesloten kolom; 1D, exact)"),
    note="Schaal (quad4, DKT) op een vlak van het fluidum: automatisch "
         "koppelvlak."))


# ------------------------------------------------------------------
# VA-04 Buis van twee media zonder constructie: lucht (rho 1,2, c 340)
# 0..0,5 m en een tweede gas (rho 1,2, c 680) 0,5..1,0 m, harde einden,
# doorsnede 0,1 x 0,1 m; hex20 20x1x1. Gemeten: met een dichtheidssprong
# (rho 12 in het tweede gas) gaf dezelfde buis 188,23 Hz i.p.v. 306,47 Hz:
# bij gedeelde knopen is PHI continu, wat met p = rho w^2 PHI alleen bij
# gelijke dichtheid klopt; de bridge geeft dan een ERROR.
# ------------------------------------------------------------------
RA, CA_, RB, CB = 1.2, 340.0, 1.2, 680.0


def _two_eq(f):
    w = 2.0 * math.pi * f
    a, b = w / CA_ * 0.5, w / CB * 0.5
    # tan(a)/Za + tan(b)/Zb = 0, vermenigvuldigd met cos(a) cos(b)
    return (math.sin(a) * math.cos(b) / (RA * CA_) +
            math.cos(a) * math.sin(b) / (RB * CB))


TWO = _roots(_two_eq, 1500.0)


def build_va04():
    m = Model("V30", attrs={"maxmode": 3})
    m.elset("Lucht", "Lucht")
    m.elset("Gas", "Gas")
    hexmesh(m, "Lucht", 10, 1, 1, box(0, 0.5, 0, 0.1, 0, 0.1), order=2)
    hexmesh(m, "Gas", 10, 1, 1, box(0.5, 1.0, 0, 0.1, 0, 0.1), order=2)
    m.mat("Lucht", *AIR)
    m.mat("Gas", '<mechanical type="Isotropic" speedofsound="680 m/s" />',
          u'<density type="Simple" density="1.2 kg/m³" />')
    return m


TESTS.append(Case(
    "VA-04", GROUP, "Buis van twee media (c 340 / 680 m/s, gelijke "
    "dichtheid, alleen fluidum), hex20 20x1x1", "vibroacoustic_modes_3d",
    build_va04,
    _mode_checks(TWO[:3], 0.005, "tan(k1L1)/Z1 + tan(k2L2)/Z2 = 0 (druk "
                 "en verplaatsing continu; 1D, exact)"),
    note="Geen constructie: 3D_FLUIDE met twee materialen; de 0-Hz-mode "
         "(constante potentiaal) valt onder vibro_freq_min."))


# ------------------------------------------------------------------
# VA-05 Paneel op een luchtholte: doos 0,5 x 0,4 x 0,3 m (harde wanden),
# bovenwand z = 0,3 m een vrij opgelegde stalen plaat (DKT quad4 30x24,
# t = 2 mm); lucht hex8 30x24x12. VA-05: lucht met een 1000 x kleinere
# dichtheid (rho 0,0012): de koppeling verdwijnt, dus de modes zijn de
# plaatmodes in vacuum (Kirchhoff, Leissa par. 4.1) en de holtemodes met
# harde wanden (f_lmn). VA-05b: gewone lucht (rho 1,2): de pompmode (1,1)
# wordt stijver door de luchtveer.
# ------------------------------------------------------------------
AX, BY, HZ, TPL = 0.5, 0.4, 0.3, 0.002
EP, NUP, RHOP = 210e9, 0.3, 7800.0
DP = EP * TPL ** 3 / (12 * (1 - NUP ** 2))


def _plate(mm, nn):
    return math.pi / 2 * (mm ** 2 / AX ** 2 + nn ** 2 / BY ** 2) * math.sqrt(
        DP / (RHOP * TPL))


def _cavity(l, mm, n):
    return C2 / 2 * math.sqrt((l / AX) ** 2 + (mm / BY) ** 2 + (n / HZ) ** 2)


PANEL_MODES = sorted(
    [(_plate(a, b), "plaat (%d,%d)" % (a, b)) for a in range(1, 6)
     for b in range(1, 6)] +
    [(_cavity(l, a, n), "holte (%d,%d,%d)" % (l, a, n)) for l in range(3)
     for a in range(3) for n in range(3) if l + a + n > 0])[:8]
# Pompmode (1,1) met de luchtveer (een mode, uniforme druk): w^2 = w0^2 +
# rho c^2 (int phi dA)^2 / (V M), phi = sin sin: int phi = 4ab/pi^2,
# M = rho_p t ab/4 (benadering; zelf afgeleid, zoals de klassieke
# paneel-holte-koppeling, o.a. Pretlove 1965; niet lokaal ingezien).
_W0 = 2 * math.pi * _plate(1, 1)
F_PUMP = math.sqrt(_W0 ** 2 + RHO2 * C2 ** 2 * 64.0 * AX * BY / (
    math.pi ** 4 * AX * BY * HZ * RHOP * TPL)) / (2 * math.pi)


def build_va05(rho_air, maxmode=8, aster=None):
    def build():
        m = Model("V30", attrs={"maxmode": maxmode})
        if aster:
            m.add("<ccxmodeldefinition>%s</ccxmodeldefinition>" % aster)
        m.elset("Paneel", "Paneelstaal")
        m.elset("Lucht", "Lucht")
        q = quadmesh(m, "Paneel", 30, 24, lambda u, v: (AX * u, BY * v, HZ))
        hexmesh(m, "Lucht", 30, 24, 12, box(0, AX, 0, BY, 0, HZ), order=1)
        m.mat("Paneelstaal", '<geometric type="Plate" thickness="2 mm" />',
              '<mechanical type="Isotropic" youngsmodulus="210000 MPa" '
              'poissonratio="0.3" />',
              u'<density type="Simple" density="7800 kg/m³" />')
        m.mat("Lucht", '<mechanical type="Isotropic" speedofsound="340 m/s" '
              '/>', u'<density type="Simple" density="%g kg/m³" />'
              % rho_air)
        rand = set()
        for s in ("u0", "u1", "v0", "v1"):
            rand |= q["nodes"][s]
        m.node_selection("Rand", sorted(rand))
        for d in ("<x>1</x><y>0</y><z>0</z>", "<x>0</x><y>1</y><z>0</z>",
                  "<x>0</x><y>0</y><z>1</z>"):
            m.add('<displacement selection="Rand"><value>0 m</value>%s'
                  '</displacement>' % d)
        m.node_selection("RandX", sorted(q["nodes"]["u0"] |
                                         q["nodes"]["u1"]))
        m.node_selection("RandY", sorted(q["nodes"]["v0"] |
                                         q["nodes"]["v1"]))
        m.add('<noderotation selection="RandX"><value>0 rad</value>'
              '<x>1</x><y>0</y><z>0</z></noderotation>')
        m.add('<noderotation selection="RandY"><value>0 rad</value>'
              '<x>0</x><y>1</y><z>0</z></noderotation>')
        return m
    return build


TESTS.append(Case(
    "VA-05", GROUP, "Paneel op een luchtholte, lucht 1000x lichter "
    "(ongekoppeld), DKT 30x24 + hex8 30x24x12", "vibroacoustic_modes_3d",
    build_va05(0.0012),
    [Check("f%d %s" % (k + 1, what), lambda r, k=k: r.mode(k + 1), ref=f,
           tol=0.01, unit="Hz",
           source="plaat vrij opgelegd in vacuum: f_mn = pi/2 (m^2/a^2 + "
                  "n^2/b^2) sqrt(D/(rho t)) (Leissa par. 4.1); holte met "
                  "harde wanden: f_lmn = c/2 sqrt((l/a)^2 + (m/b)^2 + "
                  "(n/h)^2)")
     for k, (f, what) in enumerate(PANEL_MODES)],
    note="Bij een verwaarloosbare koppeling zijn de modes de losse plaat- en "
         "holtemodes, in volgorde."))
TESTS.append(Case(
    "VA-05b", GROUP, "Paneel op een luchtholte, gewone lucht, DKT 30x24 + "
    "hex8 30x24x12", "vibroacoustic_modes_3d", build_va05(1.2),
    [Check("f1 pompmode (1,1) met luchtveer", lambda r: r.mode(1),
           ref=F_PUMP, tol=0.03, unit="Hz",
           source="w^2 = w0^2 + rho c^2 (int phi dA)^2 / (V M), een mode "
                  "(benadering, zelf afgeleid)")],
    note="Pompmode in vacuum %.2f Hz; met de luchtveer (benadering) %.2f Hz."
         % (_plate(1, 1), F_PUMP)))
# VA-05c: als VA-05b, de modes boven 100 Hz met #ASTER-OPTIE (voorbeeld 14
# van docs/ASTER_regels_voorbeelden.md): plaatmodes (2,1), (1,2), (2,2); de
# lucht verlaagt ze 0,2-0,5% (toegevoegde massa), de schaaldiscretisatie
# 0,2-0,4% (VA-05), dus 1% t.o.v. de plaat in vacuum.
TESTS.append(Case(
    "VA-05c", GROUP, "Als VA-05b, modes boven 100 Hz (#ASTER-OPTIE "
    "vibro_freq_min)", "vibroacoustic_modes_3d",
    build_va05(1.2, maxmode=3,
               aster="** ASTER-OPTIE vibro_freq_min = 100"),
    [Check("f%d plaat %s" % (k + 1, mn), lambda r, k=k: r.mode(k + 1),
           ref=_plate(*mn), tol=0.01, unit="Hz",
           source="plaat vrij opgelegd in vacuum (Leissa par. 4.1); lucht "
                  "en discretisatie samen < 1%")
     for k, mn in enumerate(((2, 1), (1, 2), (2, 2)))],
    note="De eerste mode (pompmode 54,96 Hz) valt onder de ondergrens."))
