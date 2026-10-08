# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Vermoeiing (STAPPENPLAN S8; bridge_fatigue.py + aster_fatigue.py,
docs/Vermoeiing.md). Schade D en levensduur per knoop met een Wohler-curve
(#ASTER-WOHLER) en Miner, uit blokken of een signaal op een statische run
(#ASTER-BLOK / #ASTER-SIGNAAL) of uit de tijdstappen van een niet-lineaire
of dynamische run.

Model: een staaf (hex8) met symmetrievlakken x = 0, y = 0, z = 0 en een
traction op x = L: eenassige, homogene spanning sigma = traction, dus de
schade is in elke knoop gelijk en met de hand na te rekenen.

Referenties: Miner (D = som n/N) met de opgegeven curve, met de hand;
ASTM E1049-85 (2017) par. 5.4.4 en fig. 6 (rainflow-voorbeeld: bereiken 3,
4, 6, 8 en 9 met 0,5 / 1,5 / 0,5 / 1 / 0,5 cycli); Code-Aster-testcases
szlz100a (POST_FATIGUE, Wohler-tabel, GOODMAN/GERBER) en szlz105a
(CALC_FATIGUE CUMUL_DOMMAGE, Downing-Socie-signaal); Code-Aster's eigen
POST_FATIGUE/CALC_FATIGUE in dezelfde job (#ASTER-POST, regel VERIF_CA_D
in de .mess).
"""
import io
import math
from verif_lib import Case, Check, Model, hexmesh, box, quadmesh, rect

GROUP = "Vermoeiing"
E = 210e9
STEEL = ('<mechanical type="Isotropic" youngsmodulus="210000 MPa" '
         'poissonratio="0.3" />')
DENS = u'<density type="Simple" density="7850 kg/m³" />'
L, A = 0.1, 0.01
MPA = 1.0e6

# Wohler-tabel T1 (amplitude [MPa], N), log-log
T1_S = [100.0, 150.0, 250.0, 400.0]
T1_N = [1.0e7, 2.0e6, 2.0e5, 1.0e4]
T1 = "S = [100, 150, 250, 400] MPa, N = [1e7, 2e6, 2e5, 1e4]"
# szlz100a / sslv135a: WHOL (amplitude [MPa], N), INTERPOL LOG, PROL LINEAIRE
WHOL = [138.0, 1000000.0, 152.0, 500000.0, 165.0, 200000.0, 180.0, 100000.0,
        200.0, 50000.0, 250.0, 20000.0, 295.0, 12000.0, 305.0, 10000.0,
        340.0, 5000.0, 430.0, 2000.0, 540.0, 1000.0, 690.0, 500.0, 930.0,
        200.0, 1210.0, 100.0, 1590.0, 50.0, 2210.0, 20.0, 2900.0, 10.0]
WHOL_S, WHOL_N = WHOL[0::2], WHOL[1::2]


def _loglog(s, xs, ys, left="EXCLU"):
    """N(S) met log-log-interpolatie (met de hand, los van de bridge);
    buiten de tabel: eerste/laatste stuk doorgetrokken, onder de tabel
    EXCLU = geen schade."""
    if s < xs[0] and left == "EXCLU":
        return float("inf")
    k = 0
    while k < len(xs) - 2 and s > xs[k + 1]:
        k += 1
    t = (math.log(s) - math.log(xs[k])) / (math.log(xs[k + 1]) -
                                            math.log(xs[k]))
    return math.exp(math.log(ys[k]) + t * (math.log(ys[k + 1]) -
                                           math.log(ys[k])))


def _n_t1(s_mpa):
    return _loglog(s_mpa, T1_S, T1_N)


def _basquin(s, sf, b):
    """S = Sf (2N)^b -> N."""
    return 0.5 * (s / sf) ** (1.0 / b)


# ------------------------------------------------------------------
# Model
# ------------------------------------------------------------------
def _bar(analysis, lines, tx=None, ty=None, nx=2, density=False,
         timesteps=None, size=(L, A, A), steel=STEEL, plastic=None):
    """Staaf met symmetrievlakken en traction tx (x = L) en ty (y = A);
    een waarde is een getal [Pa] of een lijst (t, waarde). lines: de
    #ASTER-regels (zonder '** ')."""
    m = Model(analysis, timesteps=timesteps)
    m.elset("Staaf", "Staal")
    g = hexmesh(m, "Staaf", nx, 1, 1, box(0, size[0], 0, size[1], 0,
                                          size[2]))
    m.mat("Staal", steel, *(((DENS,) if density else ()) +
                            ((plastic,) if plastic else ())))
    for name, side, ax in (("SymX", "u0", "x"), ("SymY", "v0", "y"),
                           ("SymZ", "w0", "z")):
        m.node_selection(name, sorted(g["nodes"][side]))
        d = dict((a, 1 if a == ax else 0) for a in "xyz")
        m.add('<displacement selection="{0}"><value>0 mm</value><x>{x}</x>'
              '<y>{y}</y><z>{z}</z></displacement>'.format(name, **d))

    def value(v):
        if isinstance(v, (list, tuple)):
            return "".join('<x time="%g s">%.12g Pa</x>' % (t, s)
                           for t, s in v)
        return "<x>%.12g Pa</x>" % v
    if tx is not None:
        m.face_selection("EindX", g["faces"]["u1"])
        m.add('<traction selection="EindX">{0}<y>0 Pa</y><z>0 Pa</z>'
              '</traction>'.format(value(tx)))
    if ty is not None:
        m.face_selection("EindY", g["faces"]["v1"])
        yv = value(ty).replace("<x", "<y").replace("</x>", "</y>")
        m.add('<traction selection="EindY"><x>0 Pa</x>%s<z>0 Pa</z>'
              '</traction>' % yv)
    m.add(u"<ccxmodeldefinition>" + u"\n".join(
        u"** ASTER" + l for l in lines) + u"</ccxmodeldefinition>")
    return m


def _all_equal(r, var):
    """Waarde van var als die in alle knopen gelijk is (homogene staaf);
    anders een fout."""
    vals = list(r.step()[var].values())
    lo, hi = min(vals), max(vals)
    if hi - lo > 1e-9 * max(abs(hi), 1e-300):
        raise ValueError("%s niet homogeen: %g .. %g" % (var, lo, hi))
    return hi


def _csv_rows(r):
    """Cyclustabel van de run (<model>_<run>_vermoeiing.csv naast de
    .mess)."""
    path = r.mess[:-len(".mess")] + "_vermoeiing.csv"
    lines = [l.strip() for l in io.open(path, encoding="utf-8") if l.strip()]
    head = lines[0].split(",")
    out = []
    for l in lines[1:]:
        row = dict(zip(head, l.split(",")))
        out.append(dict((k, float(v)) for k, v in row.items()
                        if k not in ("material", "layer")))
    return out


def _count(r, rng, scale):
    """Aantal cycli met bereik rng x scale in de cyclustabel."""
    return sum(row["count"] for row in _csv_rows(r)
               if abs(row["s_max"] - row["s_min"] - rng * scale) <
               1e-6 * scale)


def _ca_value(r, key="VERIF_CA_D"):
    """Waarde die #ASTER-POST-code in de .mess print (Code-Aster's eigen
    POST_FATIGUE/CALC_FATIGUE in dezelfde job)."""
    for line in io.open(r.mess, encoding="utf-8", errors="replace"):
        if line.startswith(key + " "):
            return float(line.split()[1])
    raise ValueError("geen %s in de .mess" % key)


def _miner(cycles, n_of_s):
    """sum n / N(S_a) over (laag, hoog, n) in Pa."""
    d = 0.0
    for lo, hi, n in cycles:
        N = n_of_s(0.5 * (hi - lo) / MPA)
        d += n / N
    return d


# ------------------------------------------------------------------
# VM-01 Constante amplitude (blok) op een statische staaf
# ------------------------------------------------------------------
S01 = 100.0 * MPA                  # traction = sigma_xx
N01 = _n_t1(200.0)                 # BLOK 2.0 -> S_a = 200 MPa
D01 = 1.0e4 / N01


def build_vm01():
    return _bar("S30", ["-WOHLER Staal: " + T1, "-BLOK 2.0, 0.0, 1e4"],
                tx=S01)


S01B = 120.0 * MPA                 # BLOK 1.25 -> S_a = 150 MPa
N01B = _basquin(150.0, 900.0, -0.1)
D01B = 1000.0 / N01B


def build_vm01b():
    return _bar("S30", ["-WOHLER Staal: Sf = 900 MPa, b = -0.1",
                        "-BLOK 1.25, 0.0, 1000"], tx=S01B)


# Machtsvorm (detailcategorie-achtig): S_ref = 40 MPa bij 2e6, m = 3; een
# blok boven en een onder Se = 25 MPa
N01C = 2.0e6 * (40.0 / 60.0) ** 3


def build_vm01c():
    return _bar("S30", ["-WOHLER Staal: S = 40 MPa, N = 2e6, m = 3, "
                        "Se = 25 MPa",
                        "-BLOK 0.6, 0.0, 5e4",
                        "-BLOK 0.2, 0.0, 1e9"], tx=S01)


# ------------------------------------------------------------------
# VM-02 Twee-niveau-spectrum
# ------------------------------------------------------------------
# (a) tabel T1, zonder gemiddelde: 1e5 cycli op 120 MPa + 5e3 op 300 MPa
D02A = 1.0e5 / _n_t1(120.0) + 5.0e3 / _n_t1(300.0)


def build_vm02a():
    return _bar("S30", ["-WOHLER Staal: " + T1, "-BLOK 1.2, 0.0, 1e5",
                        "-BLOK 3.0, 0.0, 5e3"], tx=S01)


# (b) Basquin + Goodman (Su = 600 MPa), gemiddelde > 0, 5 herhalingen
SU02 = 600.0
SEQ1 = 150.0 / (1.0 - 50.0 / SU02)
SEQ2 = 80.0 / (1.0 - 100.0 / SU02)
DREP02 = (1.0e4 / _basquin(SEQ1, 900.0, -0.1) +
          2.0e5 / _basquin(SEQ2, 900.0, -0.1))


def build_vm02b():
    return _bar("S30", ["-WOHLER Staal: Sf = 900 MPa, b = -0.1, "
                        "Su = 600 MPa, correctie = GOODMAN",
                        "-VERMOEIING herhalingen = 5",
                        "-BLOK 1.5, 0.5, 1e4", "-BLOK 0.8, 1.0, 2e5"],
                tx=S01)


# ------------------------------------------------------------------
# VM-03 Rainflow op het ASTM E1049-signaal (statische staaf x signaal)
# ------------------------------------------------------------------
ASTM = [-2.0, 1.0, -3.0, 5.0, -1.0, 3.0, -4.0, 4.0, -2.0]
S03 = 55.0 * MPA                   # geen amplitude op een tabelpunt
# ASTM E1049 fig. 6 (een doorgang): (bereik, aantal); bereik 3 is -2..1
ASTM_CYCLES = [(-2, 1, 0.5), (-3, 1, 0.5), (-1, 3, 1.0), (-3, 5, 0.5),
               (-4, 5, 0.5), (-4, 4, 0.5), (-2, 4, 0.5)]
ASTM_TABLE = {3: 0.5, 4: 1.5, 6: 0.5, 8: 1.0, 9: 0.5}
# Periodiek (POST_FATIGUE RAINFLOW, gemeten in Code-Aster 17.4)
PER_CYCLES = [(-1, 3, 1.0), (-2, 1, 1.0), (-3, 4, 1.0), (-4, 5, 1.0)]


def _scaled(cycles, s):
    return [(lo * s, hi * s, n) for lo, hi, n in cycles]


def _astm_table_check():
    """Controle van ASTM_CYCLES tegen de tabel van de norm."""
    t = {}
    for lo, hi, n in ASTM_CYCLES:
        t[hi - lo] = t.get(hi - lo, 0.0) + n
    assert t == ASTM_TABLE, t


_astm_table_check()
D03A = _miner(_scaled(ASTM_CYCLES, S03), _n_t1)
D03B = _miner(_scaled(PER_CYCLES, S03), _n_t1)

# Code-Aster's POST_FATIGUE op hetzelfde signaal x S03 (in Pa)
_POST_FATIGUE = [
    " WH_ = DEFI_FONCTION(NOM_PARA='SIGM', INTERPOL='LOG', VALE=(%s))" %
    ", ".join("%r, %r" % (s * MPA, n) for s, n in zip(T1_S, T1_N)),
    " SG_ = DEFI_FONCTION(NOM_PARA='INST', VALE=(%s))" %
    ", ".join("%r, %r" % (float(i), v * S03) for i, v in enumerate(ASTM)),
    " MT_ = DEFI_MATERIAU(FATIGUE=_F(WOHLER=WH_))",
    " TF_ = POST_FATIGUE(CHARGEMENT='UNIAXIAL', HISTOIRE=_F(SIGM=SG_), "
    "COMPTAGE='RAINFLOW', DOMMAGE='WOHLER', MATER=MT_, CUMUL='LINEAIRE')",
    " print('VERIF_CA_D', [v for v in TF_.EXTR_TABLE().values()"
    "['DOMM_CUMU'] if v is not None][0])",
]


def build_vm03a():
    return _bar("S30", ["-WOHLER Staal: " + T1,
                        "-VERMOEIING telling = ASTM",
                        "-SIGNAAL " + ", ".join("%g" % v for v in ASTM)],
                tx=S03)


def build_vm03b():
    return _bar("S30", ["-WOHLER Staal: " + T1,
                        "-SIGNAAL " + ", ".join("%g" % v for v in ASTM[:5]),
                        "-SIGNAAL " + ", ".join("%g" % v for v in ASTM[5:])]
                + ["-POST" + l for l in _POST_FATIGUE], tx=S03)


# ------------------------------------------------------------------
# VM-04 szlz105a in 3D: Downing-Socie-signaal, niet-lineair statisch
# ------------------------------------------------------------------
# szlz105a: sigma(t) = 3, 1, 4, 2, 4, 1, 2, -3 (t = 1..8; hier vanaf 0 op
# t = 0, zoals de bridge rekent), Wohler N = 1000 - 100 S (LIN), CALC_FATIGUE
# CUMUL_DOMMAGE VMIS_SG: D = 4,8133e-3 (VALE_REFE). Hier S in MPa.
SZ105 = [(0, 0.0), (1, 3.0), (2, 1.0), (3, 4.0), (4, 2.0), (5, 4.0), (6, 1.0),
         (7, 2.0), (8, -3.0)]
D105 = 4.8133e-3
_CALC_FATIGUE = [
    " WH_ = DEFI_FONCTION(NOM_PARA='SIGM', PROL_DROITE='LINEAIRE', "
    "PROL_GAUCHE='LINEAIRE', VALE=(0.0, 1000.0, 10.0e6, 0.0))",
    " MT_ = DEFI_MATERIAU(FATIGUE=_F(WOHLER=WH_))",
    " DM_ = CALC_FATIGUE(TYPE_CALCUL='CUMUL_DOMMAGE', OPTION='DOMA_ELNO_SIGM',"
    " HISTOIRE=_F(RESULTAT=res, EQUI_GD='VMIS_SG'), DOMMAGE='WOHLER', "
    "MATER=MT_)",
    " print('VERIF_CA_D', max(DM_.getValues()))",
]


def build_vm04():
    return _bar("L30", ["-WOHLER Staal: S = [0, 10] MPa, N = [1000, 0], "
                        "INTERPOL = LIN"] +
                ["-POST" + l for l in _CALC_FATIGUE],
                tx=[(t, s * MPA) for t, s in SZ105],
                timesteps={"timeperiod": "8 s", "stepsize": "1 s"})


# ------------------------------------------------------------------
# VM-05 szlz100a in 3D: cyclus 600/-500 MPa, Wohler WHOL, Goodman/Gerber
# ------------------------------------------------------------------
SZ100 = [50.0, 600.0, 50.0, -500.0, 50.0]
D100 = {"GEEN": 1.0532567e-3, "GERBER": 1.063631e-3, "GOODMAN": 1.250219e-3}


def build_vm05(corr):
    def build():
        extra = (", Su = 850 MPa, correctie = %s" % corr
                 if corr != "GEEN" else "")
        return _bar("S30", [
            "-WOHLER Staal: S = [%s] MPa, N = [%s], PROL_GAUCHE = LINEAIRE"
            % (", ".join("%g" % s for s in WHOL_S),
               ", ".join("%g" % n for n in WHOL_N)) + extra,
            "-SIGNAAL " + ", ".join("%g" % v for v in SZ100)], tx=MPA)
    return build


# ------------------------------------------------------------------
# VM-06 Dynamisch (DYNA_NON_LINE): schade = Code-Aster's CALC_FATIGUE
# ------------------------------------------------------------------
# Een hex8 (alle knopen in een element: SIGM_NOEU = SIGM_ELNO), een
# traction met keerpunten; de dynamica geeft er trillingen bij. Referentie:
# CALC_FATIGUE CUMUL_DOMMAGE (VMIS_SG, DOMA_ELNO_SIGM) op hetzelfde
# resultaat met dezelfde curve (log-log, onder de tabel geen schade: de
# standaard PROL EXCLU van beide; Code-Aster's PROL 'LINEAIRE' trekt in de
# vermoeiingsroutines lineair in (S, N) door, niet log-log: gemeten, zie
# docs/Vermoeiing.md).
TX06 = [(0.0, 0.0), (0.002, 300e6), (0.004, -150e6), (0.006, 250e6),
        (0.008, -50e6), (0.010, 0.0)]
_CALC_FATIGUE_06 = [
    " WH_ = DEFI_FONCTION(NOM_PARA='SIGM', INTERPOL='LOG', VALE=(%s))" %
    ", ".join("%r, %r" % (s * MPA, n) for s, n in zip(T1_S, T1_N)),
    " MT_ = DEFI_MATERIAU(FATIGUE=_F(WOHLER=WH_))",
    " DM_ = CALC_FATIGUE(TYPE_CALCUL='CUMUL_DOMMAGE', OPTION='DOMA_ELNO_SIGM',"
    " HISTOIRE=_F(RESULTAT=res, EQUI_GD='VMIS_SG'), DOMMAGE='WOHLER', "
    "MATER=MT_)",
    " print('VERIF_CA_D', max(DM_.getValues()))",
]


def build_vm06():
    return _bar("D31", ["-WOHLER Staal: " + T1]
                + ["-POST" + l for l in _CALC_FATIGUE_06],
                tx=TX06, nx=1, density=True, size=(0.01, 0.01, 0.01),
                timesteps={"timeperiod": "0.01 s",
                           "automatictimestepping": "false",
                           "outputtimes": "regular",
                           "timepointsstart": "0 s",
                           "timepointsincrement": "0.0001 s"})


# ------------------------------------------------------------------
# VM-07 Equivalente spanning: VMIS_SG en PRIN (tweeassig)
# ------------------------------------------------------------------
SX07, SY07 = 100.0, -150.0
VM07 = math.sqrt(SX07 ** 2 + SY07 ** 2 - SX07 * SY07)   # spoor < 0: teken -
D07 = {"VMIS_SG": 1000.0 / _basquin(VM07, 900.0, -0.1),
       "PRIN": 1000.0 / _basquin(abs(SY07), 900.0, -0.1)}


def build_vm07(kind):
    def build():
        return _bar("S30", ["-WOHLER Staal: Sf = 900 MPa, b = -0.1",
                            "-VERMOEIING spanning = %s" % kind,
                            "-BLOK 1.0, 0.0, 1000"],
                    tx=SX07 * MPA, ty=SY07 * MPA)
    return build


SRC_MINER = "Miner D = n/N, N met de hand uit de opgegeven curve"
TESTS = [
    Case("VM-01a", GROUP, "Constante amplitude (blok), Wohler-tabel",
         "lin_static_3d", build_vm01,
         [Check("schade D", lambda r: _all_equal(r, "damage"), ref=D01,
                tol=1e-3, source=SRC_MINER + " (log-log tussen 150 en "
                "250 MPa)"),
          Check("levensduur", lambda r: _all_equal(r, "life"),
                ref=1.0 / D01, tol=1e-3, unit="herhalingen",
                source="1 / D")],
         note="sigma = 100 MPa, #ASTER-BLOK 2.0, 0, 1e4: S_a = 200 MPa."),
    Case("VM-01b", GROUP, "Constante amplitude, Basquin", "lin_static_3d",
         build_vm01b,
         [Check("schade D", lambda r: _all_equal(r, "damage"), ref=D01B,
                tol=1e-3, source="S = Sf (2N)^b, Miner")]),
    Case("VM-01c", GROUP, "Machtsvorm met Se: blok onder Se zonder schade",
         "lin_static_3d", build_vm01c,
         [Check("schade D", lambda r: _all_equal(r, "damage"),
                ref=5.0e4 / N01C, tol=1e-3,
                source="N = 2e6 (40/S)^3; S = 20 MPa < Se telt niet")]),
    Case("VM-02a", GROUP, "Twee-niveau-spectrum, Wohler-tabel",
         "lin_static_3d", build_vm02a,
         [Check("schade D", lambda r: _all_equal(r, "damage"), ref=D02A,
                tol=1e-3, source=SRC_MINER + " (twee blokken)")]),
    Case("VM-02b", GROUP, "Twee-niveau-spectrum, Basquin + Goodman, 5 "
         "herhalingen", "lin_static_3d", build_vm02b,
         [Check("schade D", lambda r: _all_equal(r, "damage"),
                ref=5.0 * DREP02, tol=1e-3,
                source="Goodman S_a/(1 - S_m/Su), Basquin, Miner x 5"),
          Check("levensduur", lambda r: _all_equal(r, "life"),
                ref=1.0 / DREP02, tol=1e-3, unit="herhalingen",
                source="1 / D per herhaling")]),
    Case("VM-03a", GROUP, "Rainflow ASTM E1049 (een doorgang)",
         "lin_static_3d", build_vm03a,
         [Check("cycli bereik %d" % rng,
                (lambda g: lambda r: _count(r, g, S03))(rng), ref=n,
                tol_abs=1e-12, unit="cycli",
                source="ASTM E1049-85 fig. 6")
          for rng, n in sorted(ASTM_TABLE.items())] +
         [Check("schade D", lambda r: _all_equal(r, "damage"), ref=D03A,
                tol=1e-3, source=SRC_MINER + " op de ASTM-cycli")]),
    Case("VM-03b", GROUP, "Rainflow periodiek op het ASTM-signaal",
         "lin_static_3d", build_vm03b,
         [Check("aantal cycli", lambda r: sum(row["count"]
                                             for row in _csv_rows(r)),
                ref=4.0, tol_abs=1e-12, unit="cycli",
                source="POST_FATIGUE RAINFLOW (4 gesloten cycli)"),
          Check("schade D (hand)", lambda r: _all_equal(r, "damage"),
                ref=D03B, tol=1e-3, source=SRC_MINER),
          Check("schade D / POST_FATIGUE", lambda r: _all_equal(
              r, "damage") / _ca_value(r), ref=1.0, tol=1e-6,
                source="Code-Aster POST_FATIGUE RAINFLOW in dezelfde job")]),
    Case("VM-04", GROUP, "szlz105a in 3D (niet-lineair statisch)",
         "nonlin_static_3d", build_vm04,
         [Check("schade D", lambda r: _all_equal(r, "damage"), ref=D105,
                tol=1e-2, source="szlz105a VALE_REFE (Downing-Socie)"),
          Check("schade D / CALC_FATIGUE", lambda r: _all_equal(
              r, "damage") / _ca_value(r), ref=1.0, tol=1e-6,
                source="Code-Aster CALC_FATIGUE CUMUL_DOMMAGE op res")],
         settings={"deformation": "PETIT"}),
    Case("VM-05a", GROUP, "szlz100a in 3D: zonder correctie",
         "lin_static_3d", build_vm05("GEEN"),
         [Check("schade D", lambda r: _all_equal(r, "damage"),
                ref=D100["GEEN"], tol=1e-2,
                source="szlz100a VALE_REFE (TAB_1, cyclus 1)")]),
    Case("VM-05b", GROUP, "szlz100a in 3D: Gerber", "lin_static_3d",
         build_vm05("GERBER"),
         [Check("schade D", lambda r: _all_equal(r, "damage"),
                ref=D100["GERBER"], tol=1e-2,
                source="szlz100a VALE_REFE (TAB_3)")]),
    Case("VM-05c", GROUP, "szlz100a in 3D: Goodman", "lin_static_3d",
         build_vm05("GOODMAN"),
         [Check("schade D", lambda r: _all_equal(r, "damage"),
                ref=D100["GOODMAN"], tol=1e-2,
                source="szlz100a VALE_REFE (TAB_4)")]),
    Case("VM-06", GROUP, "Dynamisch: schade = CALC_FATIGUE",
         "dyna_impl_3d", build_vm06,
         [Check("schade D / CALC_FATIGUE", lambda r: r.extreme("damage") /
                _ca_value(r), ref=1.0, tol=1e-6,
                source="Code-Aster CALC_FATIGUE CUMUL_DOMMAGE op res")]),
    Case("VM-07a", GROUP, "Tweeassig, VMIS_SG", "lin_static_3d",
         build_vm07("VMIS_SG"),
         [Check("schade D", lambda r: _all_equal(r, "damage"),
                ref=D07["VMIS_SG"], tol=1e-3,
                source="von Mises sqrt(sx^2 + sy^2 - sx sy), Basquin")]),
    Case("VM-07b", GROUP, "Tweeassig, PRIN", "lin_static_3d",
         build_vm07("PRIN"),
         [Check("schade D", lambda r: _all_equal(r, "damage"),
                ref=D07["PRIN"], tol=1e-3,
                source="grootste |hoofdspanning| = 150 MPa, Basquin")]),
]


# ==================================================================
# S8b: multiaxiaal (CALC_FATIGUE FATIGUE_MULTI), schalen, Manson-Coffin
# ==================================================================
STEEL200 = ('<mechanical type="Isotropic" youngsmodulus="200000 MPa" '
            'poissonratio="0.3" />')
WHOL_LINE = ("-WOHLER Staal: S = [%s] MPa, N = [%s], PROL_GAUCHE = LINEAIRE"
             % (", ".join("%g" % s for s in WHOL_S),
                ", ".join("%g" % n for n in WHOL_N)))


def _n_whol(s_mpa):
    return _loglog(s_mpa, WHOL_S, WHOL_N, left="LINEAIRE")


# ------------------------------------------------------------------
# VM-08 sslv135a in 3D: kubus, tweeassige wisselende trek (sigma_xx =
# +-100, sigma_yy = -+200 MPa), Matake en Dang Van (periodiek), Dang Van
# niet-periodiek. Referentie: TEST_RESU van sslv135a (knoop N1, DOMA_NOEUD):
# Matake SIGEQ1 300, NBRUP1 10946, ENDO1 9,1356e-5; Dang Van 275, 14903,
# 6,7100e-5 (VALE_REFE, a = 1, c = 1,5).
# ------------------------------------------------------------------
TX08 = [(0.0, 0.0), (1.0, 100e6), (2.0, 0.0), (3.0, -100e6), (4.0, 0.0)]
TY08 = [(0.0, 0.0), (1.0, -200e6), (2.0, 0.0), (3.0, 200e6), (4.0, 0.0)]
SSLV135A = {"MATAKE": (300.0, 10946.0, 9.1356470832402e-05),
            "DANG_VAN": (275.0, 14903.0, 6.7099587679274e-05)}


def build_vm08(crit, loading="PERIODIEK", ramp=False):
    def build():
        tx, ty, ts = TX08, TY08, {"timeperiod": "4 s", "stepsize": "0.5 s"}
        if ramp:
            # alleen 0 -> belasting op t = 1: de cyclus is de begintoestand
            # (t = 0) en t = 1
            tx, ty = TX08[:2], TY08[:2]
            ts = {"timeperiod": "1 s", "stepsize": "0.5 s"}
        return _bar("L30", [WHOL_LINE,
                            "-VERMOEIING criterium = %s, a = 1, c = 1.5, "
                            "belasting = %s" % (crit, loading)],
                    tx=tx, ty=ty, nx=1, size=(0.01, 0.01, 0.01),
                    steel=STEEL200, timesteps=ts)
    return build


def _vm08_checks(crit):
    eq, nb, endo = SSLV135A[crit]
    src = "sslv135a VALE_REFE (CALC_FATIGUE %s_MODI_AC, knoop N1)" % crit
    return [Check("SIGEQ", lambda r: _all_equal(r, "fatigue_eq") / MPA,
                  ref=eq, tol=1e-2, unit="MPa", source=src + " SIGEQ1"),
            Check("levensduur", lambda r: _all_equal(r, "life"), ref=nb,
                  tol=1e-2, unit="cycli", source=src + " NBRUP1 (afgerond)"),
            Check("schade D", lambda r: _all_equal(r, "damage"), ref=endo,
                  tol=1e-2, source=src + " ENDO1"),
            Check("schade D (hand)", lambda r: _all_equal(r, "damage"),
                  ref=1.0 / _n_whol(eq), tol=1e-6,
                  source="1/N(SIGEQ), N log-log uit WHOL")]


# ------------------------------------------------------------------
# VM-09 Trek + afschuiving (proportioneel, zoals trek-torsie in een
# wandpunt): homogene sigma_xx = 120, tau_xy = 80 MPa op een kubus met
# 3-2-1-opleggingen, lineair statisch met een signaal of blokken.
# Hand: R = sqrt((s/2)^2 + t^2) = 100 MPa (grootste schuifspanning, Mohr);
# bij lambda in [lmin, lmax]: dtau/2 = R (lmax - lmin)/2, P_max =
# lmax s/3, N_max (vlak van de grootste schuifspanning) = lmax s/2;
# Dang Van SIGEQ = c (dtau/2 + a max(P_max, 0)), Matake c (dtau/2 + a
# N_max), Crossland c (sqrt(J2_a) + a P_max) met sqrt(J2) = sqrt(s^2/3 +
# t^2) (lmax - lmin)/2.
# ------------------------------------------------------------------
SX09, TX09 = 120.0, 80.0
A09, C09, B09 = 0.3, 1.5, 110.0
R09 = math.sqrt((SX09 / 2) ** 2 + TX09 ** 2)
J09 = math.sqrt(SX09 ** 2 / 3.0 + TX09 ** 2)
SIG09 = [0.0, 1.0, 0.0, -1.0, 0.0]


def _eq09(crit, lmin, lmax):
    amp = 0.5 * (lmax - lmin)
    if crit == "DANG_VAN":
        return C09 * (R09 * amp + A09 * max(lmax * SX09 / 3.0, 0.0))
    if crit == "MATAKE":
        return C09 * (R09 * amp + A09 * max(lmax * SX09 / 2.0, 0.0))
    return C09 * (J09 * amp + A09 * max(lmax * SX09 / 3.0, 0.0))


def build_vm09(crit, blocks=None, curve=T1):
    def build():
        m = Model("S30")
        m.elset("Kubus", "Staal")
        g = hexmesh(m, "Kubus", 2, 2, 2, box(0, A, 0, A, 0, A))
        m.mat("Staal", STEEL)
        for name, pt, d in (("P1", (0, 0, 0), (1, 1, 1)),
                            ("P2", (A, 0, 0), (0, 1, 1)),
                            ("P3", (0, A, 0), (0, 0, 1))):
            m.node_selection(name, [m.node(*pt)])
            for ax, on in zip("xyz", d):
                if on:
                    dd = dict((a, 1 if a == ax else 0) for a in "xyz")
                    m.add('<displacement selection="{0}"><value>0 mm</value>'
                          '<x>{x}</x><y>{y}</y><z>{z}</z></displacement>'
                          .format(name, **dd))
        for side, t in (("u1", (SX09, TX09)), ("u0", (-SX09, -TX09)),
                        ("v1", (TX09, 0.0)), ("v0", (-TX09, 0.0))):
            m.face_selection("F" + side, g["faces"][side])
            m.add('<traction selection="F%s"><x>%.12g Pa</x><y>%.12g Pa</y>'
                  '<z>0 Pa</z></traction>' % (side, t[0] * MPA, t[1] * MPA))
        lines = ["-WOHLER Staal: " + curve,
                 "-VERMOEIING criterium = %s, a = %g, b = %g MPa, c = %g"
                 % (crit, A09, B09, C09)]
        if blocks:
            lines += ["-BLOK %g, %g, %g" % b for b in blocks]
        else:
            lines.append("-SIGNAAL " + ", ".join("%g" % v for v in SIG09))
        m.add(u"<ccxmodeldefinition>" + u"\n".join(
            u"** ASTER" + l for l in lines) + u"</ccxmodeldefinition>")
        return m
    return build


def _vm09_checks(crit):
    eq = _eq09(crit, -1.0, 1.0)
    # Matake: N_max op het gevonden vlak; Code-Aster zoekt de normaal op
    # een discreet raster (gemeten: SIGEQ +0,06 %, Dang Van -2e-6), dus de
    # tolerantie van het klaar-criterium (1 %)
    tol = 1e-2 if crit == "MATAKE" else 1e-3
    return [Check("SIGEQ", lambda r: _all_equal(r, "fatigue_eq") / MPA,
                  ref=eq, tol=tol, unit="MPa",
                  source="hand: %s, R = 100 MPa (Mohr), P_max = 40 MPa"
                  % crit),
            Check("criterium (tau_a + a P)/b", lambda r: _all_equal(
                r, "fatigue_ratio"), ref=eq / C09 / B09, tol=tol,
                  source="hand: SIGEQ / (c b)"),
            Check("schade D", lambda r: _all_equal(r, "damage"),
                  ref=1.0 / _n_t1(eq), tol=tol,
                  source="Miner 1/N(SIGEQ), tabel T1 log-log")]


B09B = [(0.5, 0.5, 1.0e4), (1.0, 0.0, 100.0)]
EQ09B = [_eq09("DANG_VAN", m_ - a_, m_ + a_) for a_, m_, _n in B09B]
D09B = sum(n_ / _basquin(e, 900.0, -0.1) for e, (_a, _m, n_) in
           zip(EQ09B, B09B))


# ------------------------------------------------------------------
# VM-10 Schaalplaat onder zuivere buiging: boven- en onderzijde
# ------------------------------------------------------------------
# Vrije plaat L x b x t (DKT quad4) met gelijke en tegengestelde
# randmomenten m [N m/m] om y op x = 0 en x = L (consistente knoopmomenten:
# m h/2 op de hoekknopen, m h op de middenknoop) en 3-2-1-opleggingen:
# zuivere buiging, M_x = m, M_y = 0 (anticlastisch vrij), sigma = 6 m / t^2
# boven (trek bij +M_y op x = L) en onder (druk). Blok 0..1 x de last met
# Goodman: boven S_m = +S_a, onder S_m = -S_a.
L10, B10, T10 = 0.1, 0.02, 0.002
S10 = 150.0                         # MPa
M10 = S10 * MPA * T10 ** 2 / 6.0    # N m/m
SU10 = 600.0
D10 = {"top": 1.0e5 / _basquin(75.0 / (1.0 - 75.0 / SU10), 900.0, -0.1),
       "bottom": 1.0e5 / _basquin(75.0 / (1.0 + 75.0 / SU10), 900.0,
                                  -0.1)}


def build_vm10():
    m = Model("S30")
    m.elset("Plaat", "Plaat2")
    nu, nv = 10, 2
    quadmesh(m, "Plaat", nu, nv, rect(0, L10, 0, B10))
    m.mat("Plaat2", '<geometric type="Plate" thickness="%g m" />' % T10,
          STEEL)
    h = B10 / nv
    for x, sgn in ((0.0, -1.0), (L10, 1.0)):
        for j in range(nv + 1):
            nid = m.node(x, j * h, 0.0)
            name = "M%d_%d" % (int(x > 0), j)
            m.node_selection(name, [nid])
            w = h / 2.0 if j in (0, nv) else h
            m.add('<moment selection="%s"><x>0 N.m</x><y>%.12g N.m</y>'
                  '<z>0 N.m</z></moment>' % (name, sgn * M10 * w))
    for name, pt, d in (("P1", (0, 0, 0), (1, 1, 1)),
                        ("P2", (L10, 0, 0), (0, 1, 1)),
                        ("P3", (0, B10, 0), (0, 0, 1))):
        m.node_selection(name, [m.node(*pt)])
        for ax, on in zip("xyz", d):
            if on:
                dd = dict((a, 1 if a == ax else 0) for a in "xyz")
                m.add('<displacement selection="{0}"><value>0 mm</value>'
                      '<x>{x}</x><y>{y}</y><z>{z}</z></displacement>'
                      .format(name, **dd))
    lines = ["-WOHLER Plaat2: Sf = 900 MPa, b = -0.1, Su = %g MPa, "
             "correctie = GOODMAN" % SU10, "-BLOK 0.5, 0.5, 1e5"]
    m.add(u"<ccxmodeldefinition>" + u"\n".join(
        u"** ASTER" + l for l in lines) + u"</ccxmodeldefinition>")
    return m


# ------------------------------------------------------------------
# VM-11 Rek-levensduur (Manson-Coffin-Basquin)
# ------------------------------------------------------------------
# (a) staaf met opgelegde verlenging 0 -> +u -> -u -> 0 (een gesloten
# cyclus), bilineair plastisch (vloeigrens 250 MPa), niet-lineair
# statisch met DEFORMATION = PETIT: eps_xx = u/L exact (verplaatsings-
# gestuurd, homogeen); PRIN = eps_xx. (b) elastische staaf, lineair
# statisch, signaal 0, 2, -2, 0 x sigma = 100 MPa, rek = VMIS_SG:
# eps_eq = 2/3 (1 + nu) eps_xx. Hand: eps_a = Sf/E (2N)^b + ef (2N)^c,
# N opgelost met Newton in ln(2N).
SF11, B11, EF11, C11 = 900.0e6, -0.09, 0.3, -0.55
EA11 = 0.006
U11 = EA11 * L


def _manson(ea, E=E):
    x = 10.0
    for _k in range(100):
        f = SF11 / E * math.exp(B11 * x) + EF11 * math.exp(C11 * x) - ea
        df = (SF11 / E * B11 * math.exp(B11 * x) +
              EF11 * C11 * math.exp(C11 * x))
        x -= f / df
    return 0.5 * math.exp(x)


MC11 = ("-WOHLER Staal: Sf = 900 MPa, b = -0.09, ef = 0.3, c = -0.55")
PLAST11 = ('<plastic type="BilinearIsotropic" yieldstrength="250 MPa" '
           'tangentmodulus="10000 MPa" />')


def build_vm11a():
    m = _bar("L30", [MC11], nx=2, plastic=PLAST11,
             timesteps={"timeperiod": "3 s", "stepsize": "0.25 s"})
    g_nodes = [i + 1 for i, p in enumerate(m.nodes)
               if abs(p[0] - L) < 1e-12]
    m.node_selection("Trek", sorted(g_nodes))
    m.add('<displacement selection="Trek">%s<x>1</x><y>0</y><z>0</z>'
          '</displacement>' % "".join(
              '<value time="%g s">%.12g mm</value>' % (t, u * 1e3)
              for t, u in ((0.0, 0.0), (1.0, U11), (2.0, -U11), (3.0, 0.0))))
    return m


EB11 = 2.0 * 100.0e6 / E                  # eps_xx bij factor 2
EQ11B = 2.0 / 3.0 * 1.3 * EB11


def build_vm11b():
    return _bar("S30", [MC11, "-VERMOEIING rek = VMIS_SG",
                        "-SIGNAAL 0, 2, -2, 0"], tx=100.0 * MPA)


def _csv_amp(r):
    rows = _csv_rows(r)
    return max(row["amplitude"] for row in rows)


TESTS += [
    Case("VM-08a", GROUP, "sslv135a in 3D: Matake (periodiek)",
         "nonlin_static_3d", build_vm08("MATAKE"), _vm08_checks("MATAKE"),
         settings={"deformation": "PETIT"},
         note="CALC_FATIGUE FATIGUE_MULTI DOMA_NOEUD, MATAKE_MODI_AC."),
    Case("VM-08b", GROUP, "sslv135a in 3D: Dang Van (periodiek)",
         "nonlin_static_3d", build_vm08("DANG_VAN"),
         _vm08_checks("DANG_VAN"), settings={"deformation": "PETIT"},
         note="CALC_FATIGUE FATIGUE_MULTI DOMA_NOEUD, DANG_VAN_MODI_AC."),
    Case("VM-08c", GROUP, "sslv135a in 3D: Dang Van niet-periodiek",
         "nonlin_static_3d", build_vm08("DANG_VAN", "NIET_PERIODIEK"),
         [Check("schade D", lambda r: _all_equal(r, "damage"),
                ref=1.0 / _n_whol(275.0), tol=1e-2,
                source="een volledige cyclus: 1/N(275 MPa), WHOL log-log "
                "(Code-Aster ENDO1, NBRUP afgerond)")],
         settings={"deformation": "PETIT"},
         note="DANG_VAN_MODI_AV, PROJECTION UN_AXE."),
    Case("VM-08d", GROUP, "sslv135a-kubus, ramp 0 -> 1: begintoestand telt",
         "nonlin_static_3d", build_vm08("MATAKE", ramp=True),
         [Check("SIGEQ", lambda r: _all_equal(r, "fatigue_eq") / MPA,
                ref=1.5 * 75.0, tol=1e-6, unit="MPa",
                source="hand: dtau/2 = 150/2 (0..1), N_max = -50 < 0: "
                "1,5 x 75")],
         settings={"deformation": "PETIT"},
         note="Zonder kopie (CREA_RESU) slaat CALC_FATIGUE t = 0 over: "
         "56,25 MPa (gemeten)."),
    Case("VM-09a", GROUP, "Trek + afschuiving: Dang Van met de hand",
         "lin_static_3d", build_vm09("DANG_VAN"), _vm09_checks("DANG_VAN"),
         note="sigma = 120, tau = 80 MPa; signaal 0, 1, 0, -1, 0 (CREA_RESU)."),
    Case("VM-09b", GROUP, "Trek + afschuiving: Dang Van, twee blokken",
         "lin_static_3d", build_vm09("DANG_VAN", B09B,
                                     "Sf = 900 MPa, b = -0.1"),
         [Check("SIGEQ (grootste blok)", lambda r: _all_equal(
             r, "fatigue_eq") / MPA, ref=max(EQ09B), tol=1e-3, unit="MPa",
                source="hand: Dang Van, blok -1..1"),
          Check("schade D", lambda r: _all_equal(r, "damage"), ref=D09B,
                tol=1e-3, source="Miner over de blokken, Basquin, SIGEQ "
                "per blok met de hand (93 en 168 MPa)")]),
    Case("VM-09c", GROUP, "Trek + afschuiving: Matake met de hand",
         "lin_static_3d", build_vm09("MATAKE"), _vm09_checks("MATAKE")),
    Case("VM-09d", GROUP, "Trek + afschuiving: Crossland met de hand",
         "lin_static_3d", build_vm09("CROSSLAND"),
         _vm09_checks("CROSSLAND")),
    Case("VM-10", GROUP, "Schaalplaat, zuivere buiging: boven- en onderzijde",
         "lin_static_3d", build_vm10,
         [Check("sigma boven", lambda r: _all_equal(r, "vonmisesu1") / MPA,
                ref=S10, tol=1e-3, unit="MPa", source="sigma = 6 m / t^2"),
          Check("schade boven", lambda r: _all_equal(r, "damage_top"),
                ref=D10["top"], tol=1e-3,
                source="Goodman met S_m = +S_a, Basquin, Miner"),
          Check("schade onder", lambda r: _all_equal(r, "damage_bottom"),
                ref=D10["bottom"], tol=1e-3,
                source="Goodman met S_m = -S_a, Basquin, Miner"),
          Check("schade (grootste)", lambda r: _all_equal(r, "damage"),
                ref=D10["top"], tol=1e-3, source="max(boven, onder)")],
         note="DKT quad4 10 x 2, t = 2 mm, randmomenten."),
    Case("VM-11a", GROUP, "Manson-Coffin: opgelegde rek, plastisch",
         "nonlin_static_3d", build_vm11a,
         [Check("rekamplitude", _csv_amp, ref=EA11, tol=1e-6,
                source="u / L (verplaatsingsgestuurd)"),
          Check("schade D", lambda r: _all_equal(r, "damage"),
                ref=1.0 / _manson(EA11), tol=1e-3,
                source="Manson-Coffin-Basquin, N met Newton opgelost")],
         settings={"deformation": "PETIT"}),
    Case("VM-11b", GROUP, "Manson-Coffin: elastisch, rek = VMIS_SG",
         "lin_static_3d", build_vm11b,
         [Check("schade D", lambda r: _all_equal(r, "damage"),
                ref=1.0 / _manson(EQ11B), tol=1e-3,
                source="eps_eq = 2/3 (1 + nu) sigma/E x 2, Manson-Coffin")]),
]
