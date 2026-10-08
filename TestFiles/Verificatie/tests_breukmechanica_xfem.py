# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Breukmechanica met XFEM (lineair statisch, #ASTER-XFEM -> DEFI_FISS_XFEM,
MODI_MODELE_XFEM, CALC_G_XFEM, POST_K1_K2_K3; bridge_xfem.py,
docs/Breukmechanica.md par. 7). Toelichting: docs/Verificatie.md hoofdstuk
Breukmechanica (BR-07..BR-12). Het mesh heeft geen scheur: de scheur loopt
als levelset door de elementen.

BR-07  Plaat met doorgaande middenscheur onder trek (zelfde plaat als
       BR-01: a/b = 0,5, vlakke rek), half model in x, hex8 (h = a/20),
       scheur als HALFVLAK; het scheurvlak en het front lopen midden door
       de elementen. K_I tegen Tada (2000) p. 2.2.
BR-08  BR-07 met hex20.
BR-09  BR-07 met het scheurvlak door een knooplaag en het front op
       elementgrenzen (het minst gunstige geval voor de levelsets).
BR-10  Cirkelvormige scheur (penny crack, straal a) in een blok (5a) onder
       trek, kwartmodel, hex8 (h = a/20), vorm CIRKEL: K_I = 2 sigma
       sqrt(a/pi) (Tada p. 24.2), en de opening van de scheur in het midden
       uit het verplaatsingsveld dat Mecway krijgt: delta = 8 (1 - nu^2)
       sigma a / (pi E) (Sneddon; Tada p. 24.2).
BR-11  Code-Aster-testcase forma07b (dezelfde schijfvormige scheur als
       forma07a, XFEM op een met HOMARD verfijnd mesh van hex8, pyram5 en
       tetra4, forma07b.41) via de bridge: K1 MAX van CALC_G_XFEM (LAGRANGE,
       NB_POINT_FOND = 20) en POST_K1_K2_K3 (NB_POINT_FOND = 20) tegen de
       TEST_TABLE-waarden van forma07b.com2.
BR-12  Middenscheur over de volle breedte als RECHTHOEK die de dikte
       overspant: een scheur met twee fronten (NUME_FOND 1 en 2).
"""
import io
import math
import os
import shutil
import subprocess
import sys
import tempfile

import verif_lib as V
from verif_lib import Case, Check, Model, hexmesh
import tests_breukmechanica as BR

if V.SCRIPTS not in sys.path:
    sys.path.insert(0, V.SCRIPTS)

GROUP = "Breukmechanica"
E, NU = BR.E, BR.NU
E_PRIME = BR.E_PRIME
STEEL = BR.STEEL
A, B_HALF, H, SIGMA = BR.A, BR.B_HALF, BR.H, BR.SIGMA
K_TADA, G_TADA = BR.K_TADA, BR.G_TADA
SRC_K, SRC_G = BR.SRC_K, BR.SRC_G


# ------------------------------------------------------------------
# Tabel van de run (<model>_verif_xfem.csv naast de .mess)
# ------------------------------------------------------------------
def _rows(r, front=None):
    import bridge_xfem
    if not r.mess:
        raise ValueError("geen .mess")
    import bridge_fracture
    rows = bridge_fracture.read_csv(r.mess[:-len(".mess")] +
                                    bridge_xfem.CSV_SUFFIX)
    if not rows:
        raise ValueError("geen tabel %s" % bridge_xfem.CSV_SUFFIX)
    if front is not None:
        rows = [q for q in rows if int(q["front"]) == front]
    return sorted(rows, key=lambda q: q["s"])


def _with(r, key, front=None):
    rows = [q for q in _rows(r, front) if isinstance(q.get(key), float)]
    if not rows:
        raise ValueError("geen waarden voor %s" % key)
    return rows


def _mid(key, front=None):
    """Waarde van key in het frontpunt het dichtst bij het midden van de
    frontlengte."""
    def value(r):
        rows = _with(r, key, front)
        s_mid = 0.5 * (rows[0]["s"] + rows[-1]["s"])
        return min(rows, key=lambda q: abs(q["s"] - s_mid))[key]
    return value


def _max(key):
    def value(r):
        return max(q[key] for q in _with(r, key))
    return value


def _spread(key):
    """(max - min) / gemiddelde van key langs het front."""
    def value(r):
        vals = [q[key] for q in _with(r, key)]
        return (max(vals) - min(vals)) / (sum(vals) / len(vals))
    return value


def _k_from_g(front=None):
    """K uit G in het midden van het front: sqrt(G E') (vlakke rek)."""
    def value(r):
        return math.sqrt(_mid("G", front)(r) * E_PRIME)
    return value


def _xfem_line(*lines):
    return ("<ccxmodeldefinition>" + "\n".join("** ASTER-XFEM " + l
                                              for l in lines) +
            "</ccxmodeldefinition>")


# ------------------------------------------------------------------
# BR-07..09, BR-12: plaat zonder scheur in het mesh
# ------------------------------------------------------------------
def _xgrid(h, xoff, x_end=B_HALF, mirror=False):
    """x-rooster: uniform h rond x = a (verschoven over xoff), daarbuiten
    oplopend naar 0 en x_end. mirror: ook de negatieve kant (volle
    breedte, symmetrisch rond 0)."""
    x_lo, x_hi = A - 8 * h + xoff, A + 8 * h + xoff
    xs = BR._graded(0.0, x_lo, 8, h, at_start=False)[:-1] + \
        [x_lo + k * h for k in range(17)][:-1] + \
        BR._graded(x_hi, x_end, 8, h)
    if mirror:
        xs = [-x for x in reversed(xs)][:-1] + xs
    return xs


def _ygrid(h, through_nodes):
    """y-rooster symmetrisch rond 0: een element van hoogte h midden over
    y = 0 (het scheurvlak loopt door de elementen), of een knooplaag op
    y = 0 (through_nodes)."""
    if through_nodes:
        yu = [k * h for k in range(4)]
        yu = yu[:-1] + BR._graded(yu[-1], H, 12, h)
        return [-y for y in reversed(yu)][:-1] + yu
    yu = [h / 2.0 + k * h for k in range(4)]
    yu = yu[:-1] + BR._graded(yu[-1], H, 12, h)
    return [-y for y in reversed(yu)] + yu


def build_plate(order=1, h=A / 20.0, xoff=None, through_nodes=False,
                nz=4, full=False, crack=None):
    """
    Plaat x in [0, b] (symmetrie x = 0) of [-b, b] (full), y in [-H, H],
    z in [0, a] met DZ = 0 op beide zijvlakken (vlakke rek); trek sigma op
    y = +-H. Geen scheur in het mesh: #ASTER-XFEM (crack = de regels; anders
    een HALFVLAK met front x = a). xoff: verschuiving van het rooster in x
    (standaard h/2: het front loopt midden door de elementen).
    """
    if xoff is None:
        xoff = 0.0 if through_nodes else h / 2.0
    xs = _xgrid(h, xoff, mirror=full)
    ys = _ygrid(h, through_nodes)
    zs = [A * k / float(nz) for k in range(nz + 1)]
    m = Model("S30")
    m.elset("Plaat", "Staal")
    b = hexmesh(m, "Plaat", len(xs) - 1, len(ys) - 1, len(zs) - 1,
                BR._grid(xs, ys, zs), order=order)
    m.face_selection("Boven", b["faces"]["v1"])
    m.face_selection("Onder", b["faces"]["v0"])
    m.face_selection("Z0", b["faces"]["w0"])
    m.face_selection("Z1", b["faces"]["w1"])
    m.mat("Staal", STEEL)
    m.add('<traction selection="Boven"><x>0 Pa</x><y>%r Pa</y><z>0 Pa</z>'
          '</traction>' % SIGMA)
    m.add('<traction selection="Onder"><x>0 Pa</x><y>%r Pa</y><z>0 Pa</z>'
          '</traction>' % -SIGMA)
    m.add('<symmetry selection="Z0"><x>0</x><y>0</y><z>1</z></symmetry>')
    m.add('<symmetry selection="Z1"><x>0</x><y>0</y><z>1</z></symmetry>')
    tol = 1e-9

    def anchor(name, x, dirs):
        cand = [(abs(p[1]), i + 1) for i, p in enumerate(m.nodes)
                if abs(p[0] - x) < tol and abs(p[2]) < tol]
        m.node_selection(name, [min(cand)[1]])
        for dx, dy in dirs:
            m.add('<displacement selection="%s"><value>0 m</value><x>%d</x>'
                  '<y>%d</y><z>0</z></displacement>' % (name, dx, dy))
    if full:
        anchor("AnkerR", B_HALF, ((1, 0), (0, 1)))
        anchor("AnkerL", -B_HALF, ((0, 1),))
    else:
        m.face_selection("SymX", b["faces"]["u0"])
        m.add('<symmetry selection="SymX"><x>1</x><y>0</y><z>0</z>'
              '</symmetry>')
        anchor("Anker", B_HALF, ((0, 1),))
    if crack == "selectie":
        # Scheur uit een vlakselectie XFEM_midden: de vlakken op y = 0 voor
        # x < a (knooplaag, through_nodes), afwisselend van het element
        # eronder (vlak v1) en erboven (v0), zodat de bridge de normalen
        # gelijk moet richten
        j0 = [j for j, y in enumerate(ys) if abs(y) < 1e-12][0]
        faces = []
        for (i, j, k), eid in sorted(b["elems"].items()):
            if xs[i + 1] > A + 1e-12:
                continue
            if j == j0 - 1 and (i + k) % 2 == 0:
                faces.append((eid, V.HEX_FACE["v1"]))
            elif j == j0 and (i + k) % 2 == 1:
                faces.append((eid, V.HEX_FACE["v0"]))
        m.face_selection("XFEM_midden", faces)
        return m
    if crack is None:
        crack = ["vorm=HALFVLAK, front=(%r, 0, 0), normaal=(0, 1, 0), "
                 "richting=(1, 0, 0)" % A]
    m.add(_xfem_line(*crack))
    return m


def _plate_checks(front=None, label="midden"):
    return [
        Check("K_I CALC_G_XFEM %s" % label, _mid("K1", front), ref=K_TADA,
              tol=0.05, unit="Pa.m^0.5", source=SRC_K),
        Check("G CALC_G_XFEM %s" % label, _mid("G", front), ref=G_TADA,
              tol=0.05, unit="J/m2", source=SRC_K + "; " + SRC_G),
        Check("K_I POST_K1_K2_K3 %s" % label, _mid("K1_POST_K", front),
              ref=K_TADA, tol=0.05, unit="Pa.m^0.5", source=SRC_K),
        Check("G / (K_I^2/E') CALC_G_XFEM %s" % label,
              (lambda r: _mid("G", front)(r) /
               (_mid("K1", front)(r) ** 2 / E_PRIME)),
              ref=1.0, tol=0.02, unit="",
              source=SRC_G + "; G theta-methode, K interactie-integraal"),
    ]


def _lips_symmetric(r, y=A / 40.0):
    """Opening van de boven- en onderlip in het midden van de scheur (x =
    0, y = +-h/2), elk t.o.v. dezelfde hoogte aan de rand x = b (dan valt de
    starre verschuiving weg): -1 als Mecway de fysische verplaatsing van
    elke lip krijgt (de XFEM-vrijheidsgraden zelf geven ~0). y: hoogte van
    de knooplaag (bij een scheurvlak op een knooplaag de laag erboven, y =
    h: de knopen op y = 0 krijgen de waarde van een kant)."""
    z = A / 2.0
    up = r.at("disply", (0.0, y, z)) - r.at("disply", (B_HALF, y, z))
    lo = r.at("disply", (0.0, -y, z)) - r.at("disply", (B_HALF, -y, z))
    return up / lo


TESTS = [
    Case("BR-07", GROUP, "XFEM: middenscheur in een plaat (HALFVLAK), vlakke "
         "rek, hex8", "lin_static_3d", lambda: build_plate(order=1),
         _plate_checks() + [
             Check("spreiding K_I langs het front", _spread("K1"), ref=0.0,
                   tol_abs=0.01, unit="",
                   source="vlakke rek: K_I constant langs het front"),
             Check("opening boven- / onderlip (x = 0, t.o.v. x = b)",
                   _lips_symmetric, ref=-1.0, tol=0.01, unit="",
                   source="symmetrie om y = 0: de lippen gaan even ver open "
                          "(fysische verplaatsing per lip, PROJ_CHAMP)")],
         note="Zelfde plaat als BR-01 (a/b = 0,5, H/b = 3) zonder scheur in "
              "het mesh; #ASTER-XFEM vorm=HALFVLAK, front=(a, 0, 0). Het "
              "scheurvlak y = 0 en het front x = a lopen midden door de "
              "elementen (h = a/20, 4 elementen over de dikte)."),
    Case("BR-08", GROUP, "XFEM: middenscheur in een plaat, vlakke rek, hex20",
         "lin_static_3d", lambda: build_plate(order=2),
         _plate_checks(),
         note="Als BR-07 met hex20: Code-Aster 17.4 kent kwadratische "
              "XFEM-elementen (alarm XFEM_63, ajustement van de levelset)."),
    Case("BR-09", GROUP, "XFEM: middenscheur, scheurvlak door een knooplaag, "
         "front op elementgrenzen (hex8)", "lin_static_3d",
         lambda: build_plate(order=1, through_nodes=True),
         _plate_checks(),
         note="Als BR-07, maar met knopen op y = 0 en x = a: de levelset is "
              "daar 0."),
    Case("BR-12", GROUP, "XFEM: middenscheur over de volle breedte als "
         "RECHTHOEK (twee fronten), hex8", "lin_static_3d",
         lambda: build_plate(order=1, full=True, crack=[
             "vorm=RECHTHOEK, centrum=(0, 0, %r), a=%r, b=%r, "
             "vect_x=(1, 0, 0), vect_y=(0, 0, 1)" % (A / 2.0, A, A)]),
         _plate_checks(front=1, label="midden front 1") +
         _plate_checks(front=2, label="midden front 2"),
         note="Plaat x in [-b, b]; de rechthoek (halve assen a in x en a in "
              "z) steekt in z buiten de plaat (dikte a), dus er zijn twee "
              "rechte fronten x = -a en x = a (NUME_FOND 1 en 2)."),
    Case("BR-13", GROUP, "XFEM: middenscheur uit een vlakselectie (XFEM_), "
         "scheurvlak op een knooplaag, hex8", "lin_static_3d",
         lambda: build_plate(order=1, through_nodes=True, crack="selectie"),
         _plate_checks() + [
             Check("opening boven- / onderlip (x = 0, y = +-h, t.o.v. "
                   "x = b)", lambda r: _lips_symmetric(r, A / 20.0),
                   ref=-1.0, tol=0.01, unit="",
                   source="symmetrie om y = 0: de lippen gaan even ver open "
                          "(knooplagen naast het scheurvlak)")],
         note="Mesh van BR-09 (knopen op y = 0 en x = a) zonder "
              "#ASTER-XFEM: de scheur is de vlakselectie XFEM_midden (de "
              "vlakken op y = 0 voor x < a, afwisselend van de elementen "
              "onder en boven het vlak gekozen). DEFI_FISS GROUP_MA_FISS/"
              "GROUP_MA_FOND; vergelijk met BR-09 (dezelfde scheur als "
              "HALFVLAK)."),
    Case("BR-14", GROUP, "XFEM: MED-bestand voor ParaView met het mesh dat "
         "de scheur volgt (med_export)", "lin_static_3d",
         lambda: build_plate(order=1),
         [Check("velden DEPL, SIGM_NOEU, SIEQ_NOEU in _xfem.rmed",
                lambda r: float(len(_med(r)["wanted"])), ref=3.0,
                tol_abs=0.0, unit="",
                source="IMPR_RESU RESULTAT=res_xvisu (bridge_xfem."
                       "result_block); teruggelezen met MEDCoupling"),
          Check("groep Plaat in het MED-mesh", lambda r: float(
              "Plaat" in _med(r)["groups"]), ref=1.0, tol_abs=0.0, unit="",
                source="POST_MAIL_XFEM houdt de groepen van het mesh"),
          Check("opening in het MED-bestand (dubbele knopen, x = 0) / "
                "opening Mecway (y = +-h/2)", lambda r: _med_opening_ratio(r),
                ref=1.0,
                tol=0.03, unit="",
                source="de lippen hebben in het MED-mesh eigen knopen; de "
                       "knopen van Mecway liggen h/2 van het scheurvlak "
                       "(verschil sigma h / E, ~1 %)")],
         general={"med_export": True},
         note="Plaat van BR-07 (HALFVLAK, scheurvlak midden door de "
              "elementen) met de instelling med_export: het mesh van "
              "POST_MAIL_XFEM met DEPL, SIGM_NOEU en SIEQ_NOEU als "
              "<model>_<run_id>_xfem.rmed, teruggelezen met de MEDCoupling "
              "van de Code-Aster-installatie (Python 3.11 van Code-Aster)."),
]


# ------------------------------------------------------------------
# BR-14: MED-bestand teruglezen (MEDCoupling van Code-Aster)
# ------------------------------------------------------------------
_MED_READER = r'''
import json, os, sys
root, path, out = sys.argv[1:4]
for d in ("medcoupling/lib", "lib/aster", "tools", "Python311"):
    d = os.path.join(root, d)
    if os.path.isdir(d):
        os.add_dll_directory(d)
sys.path.insert(0, os.path.join(root, "medcoupling", "lib", "python3.11",
                                "site-packages"))
import medcoupling as mc
name = mc.GetMeshNames(path)[0]
mm = mc.MEDFileMesh.New(path, name)
res = {"mesh": name, "groups": list(mm.getGroupsNames()),
       "coords": mm.getCoords().getValues(), "fields": {}}
for fname in mc.GetAllFieldNames(path):
    f1 = mc.MEDFileField1TS.New(path, fname)
    arr = f1.field(mm).getArray()
    res["fields"][fname] = {"ncomp": arr.getNumberOfComponents(),
                            "comps": list(arr.getInfoOnComponents()),
                            "values": arr.getValues()}
with open(out, "w") as f:
    json.dump(res, f)
'''


def _med(r):
    """Inhoud van <model>_verif_xfem.rmed (MEDCoupling in de Python van
    Code-Aster; gecached op r)."""
    cached = getattr(r, "_xfem_med", None)
    if cached is not None:
        return cached
    import json
    import bridge_utils
    import bridge_xfem
    if not r.mess:
        raise ValueError("geen .mess")
    path = r.mess[:-len(".mess")] + bridge_xfem.MED_SUFFIX
    if not os.path.isfile(path):
        raise ValueError("geen %s" % os.path.basename(path))
    root = os.path.dirname(os.path.dirname(
        bridge_utils.exe_aster_path()))
    work = os.path.dirname(path)
    script = os.path.join(work, "med_lezen.py")
    out = os.path.join(work, "med_lezen.json")
    io.open(script, "w", encoding="utf-8").write(_MED_READER)
    subprocess.call([os.path.join(root, "Python311", "python.exe"), script,
                     root, path, out], cwd=work)
    if not os.path.isfile(out):
        raise ValueError("MED niet gelezen (MEDCoupling)")
    with io.open(out, encoding="utf-8") as f:
        data = json.load(f)
    data["wanted"] = [c for c in ("DEPL", "SIGM_NOEU", "SIEQ_NOEU")
                      if any(n.endswith(c) for n in data["fields"])]
    r._xfem_med = data
    return data


def _med_opening_ratio(r):
    """Opening van de scheur in het midden (x = 0, z = a/2) uit het MED-
    mesh: het verschil in DY tussen de knopen van de twee lippen op y = 0,
    gedeeld door de opening die Mecway krijgt (knopen op y = +-h/2)."""
    data = _med(r)
    xyz = data["coords"]
    pts = [tuple(xyz[3 * i:3 * i + 3]) for i in range(len(xyz) // 3)]
    depl = [v for n, v in data["fields"].items() if n.endswith("DEPL")][0]
    iy = [c.strip().upper() for c in depl["comps"]].index("DY")
    nc = depl["ncomp"]
    z = A / 2.0
    tol = 1e-9 * B_HALF
    lips = [i for i, p in enumerate(pts) if abs(p[0]) < tol and
            abs(p[1]) < tol and abs(p[2] - z) < tol]
    if len(lips) < 2:
        raise ValueError("%d knopen op (0, 0, a/2) in het MED-mesh, "
                         "verwacht minstens 2 (lippen)" % len(lips))
    # (op een elementgrens z = a/2 heeft elke lip een knoop per element)
    dy = [depl["values"][nc * i + iy] for i in lips]
    y = A / 40.0
    mecway = r.at("disply", (0.0, y, z)) - r.at("disply", (0.0, -y, z))
    return (max(dy) - min(dy)) / mecway


# ------------------------------------------------------------------
# BR-10 Cirkelvormige scheur in een blok
# ------------------------------------------------------------------
AP = 0.02          # straal a
WP = 5 * AP        # halve breedte en halve hoogte van het blok
SIG_P = 1e6
K_PENNY = 2.0 * SIG_P * math.sqrt(AP / math.pi)
COD_PENNY = 8.0 * (1.0 - NU ** 2) * SIG_P * AP / (math.pi * E)
SRC_PENNY = ("Tada, Paris & Irwin (2000) p. 24.2: K = 2 sigma sqrt(a/pi) "
             "(oneindig medium; blok 5a)")
SRC_COD = ("Sneddon (1946), Tada p. 24.2: opening in het midden "
           "8 (1 - nu^2) sigma a / (pi E); knopen op +-h/2 (+ sigma h / E)")


def build_penny(nfine=30, order=1):
    """Kwartmodel x, y in [0, 5a] (symmetrie x = 0 en y = 0), z in [-5a,
    5a], trek op z = +-5a; uniform h = 1,5a/nfine tot r = 1,5a en |z| <
    5,5h, het scheurvlak z = 0 midden door de elementen."""
    h = 1.5 * AP / nfine
    xs = [k * h for k in range(nfine + 1)]
    xs = xs[:-1] + BR._graded(xs[-1], WP, 8, h)
    zu = [h / 2.0 + k * h for k in range(6)]
    zu = zu[:-1] + BR._graded(zu[-1], WP, 10, h)
    zs = [-z for z in reversed(zu)] + zu
    m = Model("S30")
    m.elset("Blok", "Staal")
    b = hexmesh(m, "Blok", len(xs) - 1, len(xs) - 1, len(zs) - 1,
                BR._grid(xs, xs, zs), order=order)
    m.face_selection("SymX", b["faces"]["u0"])
    m.face_selection("SymY", b["faces"]["v0"])
    m.face_selection("Boven", b["faces"]["w1"])
    m.face_selection("Onder", b["faces"]["w0"])
    m.mat("Staal", STEEL)
    m.add('<symmetry selection="SymX"><x>1</x><y>0</y><z>0</z></symmetry>')
    m.add('<symmetry selection="SymY"><x>0</x><y>1</y><z>0</z></symmetry>')
    m.add('<traction selection="Boven"><x>0 Pa</x><y>0 Pa</y><z>%r Pa</z>'
          '</traction>' % SIG_P)
    m.add('<traction selection="Onder"><x>0 Pa</x><y>0 Pa</y><z>%r Pa</z>'
          '</traction>' % -SIG_P)
    tol = 1e-9
    cand = [(abs(p[2]), i + 1) for i, p in enumerate(m.nodes)
            if abs(p[0] - WP) < tol and abs(p[1] - WP) < tol]
    m.node_selection("Anker", [min(cand)[1]])
    m.add('<displacement selection="Anker"><value>0 m</value><x>0</x>'
          '<y>0</y><z>1</z></displacement>')
    m.add(_xfem_line("vorm=CIRKEL, centrum=(0, 0, 0), r=20 mm, "
                     "normaal=(0, 0, 1), naam=penny"))
    return m


def _cod(r):
    h = 1.5 * AP / 30
    return (r.at("displz", (0.0, 0.0, h / 2.0)) -
            r.at("displz", (0.0, 0.0, -h / 2.0)))


TESTS.append(Case(
    "BR-10", GROUP, "XFEM: cirkelvormige scheur (penny crack) in een blok "
    "onder trek, hex8", "lin_static_3d", build_penny,
    [Check("K_I POST_K1_K2_K3 midden", _mid("K1_POST_K"), ref=K_PENNY,
           tol=0.05, unit="Pa.m^0.5", source=SRC_PENNY),
     Check("K uit G (sqrt(G E')) midden", _k_from_g(), ref=K_PENNY, tol=0.05,
           unit="Pa.m^0.5", source=SRC_PENNY + "; G = K^2 (1 - nu^2)/E"),
     Check("K_I CALC_G_XFEM midden (geen krommingscorrectie)", _mid("K1"),
           ref=K_PENNY, tol=0.05, unit="Pa.m^0.5", known="KNOWN_BUGS 83",
           source=SRC_PENNY + "; CALC_G_XFEM zonder krommingscorrectie: "
                  "bekend te hoog (+5,3% bij h = a/20; "
                  "docs/Breukmechanica.md par. 7.4)"),
     Check("opening van de scheur in het midden (displz)", _cod,
           ref=COD_PENNY, tol=0.05, unit="m", source=SRC_COD)],
    timeout=3600,
    note="Kwartmodel, blok 5a x 5a x 10a, a = 20 mm, h = a/20 aan het "
         "front (hex8); #ASTER-XFEM vorm=CIRKEL, r=20 mm (lengte met "
         "eenheid). Het front snijdt de elementen schuin."))


# ------------------------------------------------------------------
# BR-11 Code-Aster-testcase forma07b via de bridge
# ------------------------------------------------------------------
# forma07b (share/aster/tests, Code-Aster 17.4): schijfvormige scheur
# (straal 2 m) in een kwartblok 10 x 10 x 20 m als XFEM-levelset
# (DEFI_FISS_XFEM ELLIPSE), op het met HOMARD verfijnde mesh forma07b.41
# (com1); DY = 0 op FACE_AV, DX = 0 op FACE_LAT, DZ = 0 in knoop D,
# PRES_REP -1e6 Pa op FACE_SUP/FACE_INF, E = 2e11, nu = 0,3. com2:
# CALC_G_XFEM en POST_K1_K2_K3 met R_INF = 2h, R_SUP = 5h (h = 0,078).
F07B_H = 0.078
F07B_R_INF, F07B_R_SUP = 2.0 * F07B_H, 5.0 * F07B_H
F07B_K1_LAG20 = 1671075.8789433436   # TEST_TABLE G_LAG_20 K1 MAX (VALE_CALC)
F07B_K1_PK20 = 1607458.569834848     # TEST_TABLE PK_20 K1 MAX (VALE_CALC)
F07B_K1_REFE = 1.5957e6              # VALE_REFE (analytisch, 1 % voor PK_20)

_CONVERT = """DEBUT()
MA = LIRE_MAILLAGE(FORMAT='MED', UNITE=20)
IMPR_RESU(FORMAT='ASTER', UNITE=21, RESU=_F(MAILLAGE=MA))
FIN()
"""

_ASTER_N = {"HEXA8": 8, "PYRAM5": 5, "TETRA4": 4, "QUAD4": 4, "TRIA3": 3,
            "SEG2": 2, "POI1": 1}


def _forma07b_mail():
    """Mesh forma07b.41 in ASTER-formaat (met Code-Aster omgezet, in een
    tijdelijke map bewaard)."""
    import bridge_utils
    work = os.path.join(tempfile.gettempdir(), "mab_forma07b")
    mail = os.path.join(work, "forma07b.mail")
    if os.path.isfile(mail) and os.path.getsize(mail) > 0:
        return mail
    if os.path.isdir(work):
        shutil.rmtree(work)
    os.makedirs(work)
    shutil.copy(os.path.join(BR._aster_tests_dir(), "forma07b.41"),
                os.path.join(work, "forma07b.41"))
    io.open(os.path.join(work, "conv.comm"), "w").write(_CONVERT)
    io.open(os.path.join(work, "conv.export"), "w").write(
        "P time_limit 120\nP memory_limit 1024\nP ncpus 1\nP mpi_nbcpu 1\n"
        "F comm %s D 1\nF mmed %s D 20\nF mail %s R 21\nF mess %s R 6\n" % (
            os.path.join(work, "conv.comm"),
            os.path.join(work, "forma07b.41"), mail,
            os.path.join(work, "conv.mess")))
    with open(os.path.join(work, "console.txt"), "w") as f:
        subprocess.call([bridge_utils.exe_aster_path(),
                         os.path.join(work, "conv.export")], stdout=f,
                        stderr=subprocess.STDOUT, cwd=work)
    if not os.path.isfile(mail):
        raise RuntimeError("forma07b-mesh niet omgezet (zie %s)" % work)
    return mail


def build_forma07b():
    import bridge_mesh
    BR._ASTER_N.update(_ASTER_N)
    nodes, cells, groups = BR.read_aster_mail(_forma07b_mail())
    m = Model("S30")
    index = {}
    for name, xyz in nodes.items():
        m.nodes.append(xyz)
        index[name] = len(m.nodes)
    m.elset("Blok", "Staal")
    shapes = {"HEXA8": ("hex8", "HEXA"), "PYRAM5": ("pyr5", "PYRAM"),
              "TETRA4": ("tet4", "TETRA")}
    face_of = {}
    for name, (atype, cn) in cells.items():
        if atype not in shapes:
            continue
        shape, base = shapes[atype]
        nids = [index[n] for n in cn]
        eid = m.elem("Blok", shape, nids)
        for fid, corners in bridge_mesh._FACE_CORNERS[base].items():
            key = frozenset(nids[i] for i in corners)
            if len(key) == len(corners):
                face_of.setdefault(key, (eid, fid))

    def faces(group):
        return [face_of[frozenset(index[n] for n in cells[c][1])]
                for c in groups[("GROUP_MA", group)]]
    for group in ("FACE_AV", "FACE_LAT", "FACE_SUP", "FACE_INF"):
        m.face_selection(group, faces(group))
    m.node_selection("D", [index[n] for n in groups[("GROUP_NO", "D")]])
    m.mat("Staal", '<mechanical type="Isotropic" youngsmodulus="200000 MPa" '
                   'poissonratio="0.3" />')
    m.add('<symmetry selection="FACE_AV"><x>0</x><y>1</y><z>0</z>'
          '</symmetry>')
    m.add('<symmetry selection="FACE_LAT"><x>1</x><y>0</y><z>0</z>'
          '</symmetry>')
    m.add('<displacement selection="D"><value>0 m</value><x>0</x><y>0</y>'
          '<z>1</z></displacement>')
    # PRES_REP = -1e6 Pa: trek, als Mecway-druk van -1 MPa
    for sel in ("FACE_SUP", "FACE_INF"):
        m.add('<surfacepressure selection="%s"><normal>-1000000 Pa</normal>'
              '</surfacepressure>' % sel)
    m.add("<ccxmodeldefinition>"
          "** ASTER-XFEM vorm=ELLIPSE, centrum=(0, 0, 0), a=2, b=2, "
          "vect_x=(1, 0, 0), vect_y=(0, 1, 0), glad=LAGRANGE, punten=20\n"
          "** ASTER-OPTIE fracture_r_inf = %r\n"
          "** ASTER-OPTIE fracture_r_sup = %r</ccxmodeldefinition>"
          % (F07B_R_INF, F07B_R_SUP))
    return m


TESTS.append(Case(
    "BR-11", GROUP, "XFEM: Code-Aster-testcase forma07b (schijfvormige "
    "scheur, levelset)", "lin_static_3d", build_forma07b,
    [Check("K1 MAX CALC_G_XFEM (LAGRANGE, 20 punten)", _max("K1"),
           ref=F07B_K1_LAG20, tol=0.01, unit="Pa.m^0.5",
           source="forma07b TEST_TABLE G_LAG_20 VALE_CALC (NON_REGRESSION)"),
     Check("K1 MAX POST_K1_K2_K3 (20 punten)", _max("K1_POST_K"),
           ref=F07B_K1_PK20, tol=0.01, unit="Pa.m^0.5",
           source="forma07b TEST_TABLE PK_20 VALE_CALC (NON_REGRESSION)"),
     Check("K1 MAX POST_K1_K2_K3 (analytisch)", _max("K1_POST_K"),
           ref=F07B_K1_REFE, tol=0.01, unit="Pa.m^0.5",
           source="forma07b TEST_TABLE PK_20 VALE_REFE (ANALYTIQUE, 1 %): "
                  "K = 2 sigma sqrt(a/pi)")],
    timeout=3600,
    note="Mesh forma07b.41 (HOMARD-verfijnd: hex8, pyr5, tet4) als .liml; "
         "FACE_AV/FACE_LAT symmetry, D DZ = 0, druk -1 MPa op FACE_SUP/"
         "FACE_INF; #ASTER-XFEM vorm=ELLIPSE (a = b = 2 m), glad=LAGRANGE, "
         "punten=20; R_INF/R_SUP = 2h/5h met #ASTER-OPTIE."))

TESTS.sort(key=lambda c: c.id)
