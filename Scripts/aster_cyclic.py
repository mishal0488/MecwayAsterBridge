# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# aster_cyclic.py
#
# Verantwoordelijkheid: de stappen van de modale analyse met cyclische
# symmetrie (bridge_cyclic.py) die BINNEN de Code-Aster-job draaien
# (CPython 3.11, zoals aster_results_export.py; NIET in IronPython):
#   collect()  de modes van MODE_ITER_CYCL (via REST_SOUS_STRUC SECTEUR=1)
#              in Mecway's volgorde als een MODE_MECA-resultaat op de
#              sector, of (instelling cyclic_full_wheel) op het hele wiel
#              (REST_SOUS_STRUC SQUELETTE), terug in het globale
#              assenstelsel van het model, met een regel BRIDGE_MODE <nr>
#              FREQ=.. ND=.. per mode.
#   wheel_model(), node_rows()  het hele wiel in Mecway's nummering
#              (ook gebruikt door aster_results_export).
#
# HET HELE WIEL ZOALS MECWAY (headless gemeten met Mecway 33 op
# TestFiles/Verificatie/modellen/CY-01, CY-04, CY-05, Samples/Impeller en
# varianten met twee elsets, een losse knoop en verwisselde master/slave):
# de <solution> bevat N kopieen van het model. Sector j (0 .. N-1) is de
# sector gedraaid over +j x 360/N graden (rechtsom-positief om de as door
# de oorsprong, onafhankelijk van welke kant master is); knoop i van
# sector j heet i + j x n (n = aantal knopen van het model), element e
# heet e + j x ne (ne = hoogste elementnummer). Elke knoop wordt
# gekopieerd, ook die op de snijvlakken en op de as (dubbel, niet
# samengevoegd). Per elset staan de elementen per sector achter elkaar
# (sector 0, 1, ...). Code-Aster's skelet (DEFI_SQUELETTE CYCLIQUE,
# SECTEUR TOUT='OUI') heeft (gezien met CY-01) dezelfde opbouw: N x n
# knopen, blok j = de sector gedraaid over +j x 360/N om Z, en op blok j
# het veld van REST_SOUS_STRUC SECTEUR=j+1 gedraaid naar dat blok;
# node_rows() koppelt de knopen toch op coordinaat.
#
# Volgorde (gemeten met Code-Aster 17.4, testcase sdls04a en de ring van
# docs/Verificatie.md, CY): REST_SOUS_STRUC SECTEUR=1 geeft de modes per
# aantal nodal diameters (ND) in de volgorde van CALCUL NB_DIAM/TOUT_DIAM
# (0, 1, ..., N/2), per ND NMAX_FREQ modes; voor 0 < ND < N/2 staat elke
# mode twee keer (het paar staande golven cos/sin, zelfde FREQ), voor ND
# 0 en ND N/2 een keer. Binnen een ND zijn ze niet altijd oplopend
# (gezien: 1420,8 voor 1404,7 Hz). Mecway (manual v33 par. 9.41, headless
# nagelopen) geeft per ND het gevraagde aantal modes, elk een keer,
# oplopend in ND en binnen een ND in frequentie: collect() doet hetzelfde
# (van een paar het eerste lid).
#
# Gebruik geen f-strings: zelfde stijl als de rest van de bridge.

# Rekenframe: MODE_ITER_CYCL eist de Z-as als symmetrieas (U4.52.05 par.
# 1). Bij een X- of Y-as draait de .comm het mesh vooraf met een
# cyclische verwisseling van de assen (MODI_MAILLAGE ROTATION 120 graden
# om (1, 1, 1); bridge_cyclic.ROTATION):
#   as Y: p' = (z, x, y), dus u = (u'_y, u'_z, u'_x);
#   as X: p' = (y, z, x), dus u = (u'_z, u'_x, u'_y).
# _BACK: de componenten van het rekenframe die DX, DY, DZ van het model
# worden.
from aster_results_export import tr   # taal van de meldingen

_BACK = {"Z": ("DX", "DY", "DZ"), "Y": ("DY", "DZ", "DX"),
         "X": ("DZ", "DX", "DY")}
_AXIS_INDEX = {"X": 0, "Y": 1, "Z": 2}

# Koppeling van de knopen op coordinaat (node_rows): binnen _MATCH x de
# modelgrootte (de bridge zet coordinaten < 1e-12 x de modelgrootte op 0,
# bridge_cyclic._rewrite_mesh; draaien geeft ~1e-16).
_MATCH = 1.0e-6


# ==========================================
# HET HELE WIEL (Mecway's nummering, zie de kop)
# ==========================================
def wheel_coords(coords, axis, n_sectors):
    """(N x n, 3)-array: rij j x n + i = knoop i+1 van sector j, gedraaid
    over +j x 360/N graden om de as (index 0/1/2 = X/Y/Z) door de
    oorsprong (zoals bridge_cyclic._rotate)."""
    import numpy as np
    coords = np.asarray(coords, dtype=float)
    i, k = (axis + 1) % 3, (axis + 2) % 3
    out = []
    for j in range(n_sectors):
        a = 2.0 * np.pi * j / n_sectors
        c, s = np.cos(a), np.sin(a)
        q = coords.copy()
        q[:, i] = c * coords[:, i] - s * coords[:, k]
        q[:, k] = s * coords[:, i] + c * coords[:, k]
        out.append(q)
    return np.concatenate(out)


def wheel_model(model, axis, n_sectors):
    """
    Het model van aster_results_export.read_liml_model ("coords",
    "elsets", "elements") uitgebreid tot het hele wiel in Mecway's
    nummering: knoop i + j x n, element e + j x ne (ne = hoogste
    elementnummer), per elset de sectoren achter elkaar.
    """
    from collections import OrderedDict
    n = model["coords"].shape[0]
    ne = max(model["elements"]) if model["elements"] else 0
    elements = {}
    elsets = OrderedDict()
    for name, eids in model["elsets"].items():
        out = elsets[name] = []
        for j in range(n_sectors):
            for e in eids:
                shape, nids = model["elements"][e]
                elements[e + j * ne] = (shape, [m + j * n for m in nids])
                out.append(e + j * ne)
    return {"coords": wheel_coords(model["coords"], axis, n_sectors),
            "elsets": elsets, "elements": elements}


def node_rows(mesh_xyz, coords):
    """
    (rijen, fout): rijen[i] = de knoop (rij) van het mesh op de plaats van
    knoop i+1 van coords (binnen _MATCH x de modelgrootte), of (None,
    melding). Eerst de identieke volgorde (zo bouwt Code-Aster het skelet),
    anders op coordinaat (raster met de buurcellen).
    """
    import numpy as np
    mesh_xyz = np.asarray(mesh_xyz, dtype=float).reshape(-1, 3)
    coords = np.asarray(coords, dtype=float)
    size = max(1.0e-30, float(np.abs(coords).max()) if coords.size else 0.0)
    tol = _MATCH * size
    if mesh_xyz.shape == coords.shape and \
            float(np.abs(mesh_xyz - coords).max()) <= tol:
        return np.arange(coords.shape[0]), None
    cells = {}
    keys = np.floor(mesh_xyz / tol).astype(np.int64)
    for r, key in enumerate(map(tuple, keys)):
        cells.setdefault(key, []).append(r)
    rows = np.empty(coords.shape[0], dtype=np.int64)
    for i, p in enumerate(coords):
        k0 = np.floor(p / tol).astype(np.int64)
        best, dist = None, None
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for r in cells.get((k0[0] + dx, k0[1] + dy, k0[2] + dz),
                                       ()):
                        d = float(np.abs(mesh_xyz[r] - p).max())
                        if dist is None or d < dist:
                            best, dist = r, d
        if best is None or dist > tol:
            return None, tr("knoop {0} van het wiel ({1}) niet in het mesh "
                            "van het resultaat",
                            "node {0} of the wheel ({1}) not in the mesh of "
                            "the result").format(i + 1, tuple(p))
        rows[i] = best
    return rows, None


def _multiplicity(nd, n_sectors):
    """1 voor ND 0 en ND N/2 (een staande golf), anders 2 (paar)."""
    return 1 if (nd == 0 or 2 * nd == n_sectors) else 2


def mode_table(freqs, n_sectors, nmax_freq, diameters):
    """
    [(frequentie, positie in REST_SOUS_STRUC (0-based), ND)] in Mecway's
    volgorde, of None als het aantal modes niet klopt met de volgorde
    hierboven (dan is de ND niet te bepalen).
    """
    expected = sum(nmax_freq * _multiplicity(d, n_sectors)
                   for d in diameters)
    if len(freqs) != expected:
        return None
    out = []
    k = 0
    for d in diameters:
        group = []
        mult = _multiplicity(d, n_sectors)
        for _j in range(nmax_freq):
            group.append((float(freqs[k]), k, d))
            k += mult
        group.sort()
        out.extend(group)
    return out


def _attach(modsec, mesh, numbering):
    """
    Het resultaat van REST_SOUS_STRUC SECTEUR=1 heeft in de Python-laag van
    Code-Aster 17.4 geen mesh en geen velden: getFieldsNames() was leeg,
    getField() gaf KeyError en CREA_CHAMP EXTR <EXCEPTION> JEVEUX_26 (geen
    mesh, gezien). Met setMesh, setDOFNumbering (de nummering van de
    sectorbasis, waarop REST_SOUS_STRUC de velden schrijft) en build() zijn
    de velden DEPL leesbaar (gezien; build() meldt dan "NUME_EQUA already
    exists", onschuldig).
    """
    modsec.setMesh(mesh)
    modsec.setDOFNumbering(numbering)
    try:
        modsec.build()
    except Exception as e:
        if "already exists" not in str(e):
            raise


def _rotate_back(target, axis, angle):
    """Mesh of skelet terug naar het assenstelsel van het model
    (MODI_MAILLAGE ROTATION om (1, 1, 1) met -angle graden)."""
    from code_aster.Commands import MODI_MAILLAGE
    if axis != "Z":
        MODI_MAILLAGE(reuse=target, MAILLAGE=target,
                      ROTATION={"POIN_1": (0.0, 0.0, 0.0),
                                "DIR": (1.0, 1.0, 1.0),
                                "ANGLE": -float(angle)})


def _modes_on(source, target, table, axis):
    """MODE_MECA op mesh target met de velden DEPL van resultaat source in
    de volgorde van table (mode_table), componenten terug (_BACK)."""
    from code_aster.Commands import CREA_CHAMP, CREA_RESU
    orders = list(source.getAccessParameters()["NUME_ORDRE"])
    affe = []
    for k, (f, pos, _nd) in enumerate(table):
        fld = source.getField("DEPL", orders[pos])
        # Een eigen veld op het mesh (ook bij de Z-as), zonder de Lagrange-
        # vrijheidsgraden van de sectorbasis, met de componenten terug in
        # het assenstelsel van het model.
        new = CREA_CHAMP(OPERATION="ASSE", TYPE_CHAM="NOEU_DEPL_R",
                         MAILLAGE=target,
                         ASSE={"CHAM_GD": fld, "TOUT": "OUI",
                               "NOM_CMP": _BACK[axis],
                               "NOM_CMP_RESU": ("DX", "DY", "DZ")})
        affe.append({"NOM_CHAM": "DEPL", "CHAM_GD": new,
                     "NUME_MODE": k + 1, "FREQ": f})
    return CREA_RESU(OPERATION="AFFE", TYPE_RESU="MODE_MECA", AFFE=affe)


def _wheel(wheel, skeleton, table, axis, angle, n_sectors, liml):
    """
    De modes op het hele wiel (REST_SOUS_STRUC SQUELETTE op het skelet van
    DEFI_SQUELETTE CYCLIQUE), of (None, reden). Controle vooraf: zelfde
    frequenties als de sector, en elke knoop van Mecway's wiel (wheel_model
    op de .liml) heeft een knoop van het skelet op dezelfde plaats.
    """
    import numpy as np
    import aster_results_export
    acc = wheel.getAccessParameters()
    freqs = [float(f) for f in acc["FREQ"]]
    if len(freqs) <= max(pos for _f, pos, _nd in table):
        return None, tr("{0} modes op het skelet",
                        "{0} modes on the skeleton").format(len(freqs))
    worst = max(abs(freqs[pos] - f) / max(abs(f), 1.0)
                for f, pos, _nd in table)
    if worst > 1.0e-9:
        return None, tr("andere frequenties op het skelet ({0:.3g})",
                        "different frequencies on the skeleton "
                        "({0:.3g})").format(worst)
    _rotate_back(skeleton, axis, angle)
    model = aster_results_export.read_liml_model(liml)
    coords = wheel_coords(model["coords"], _AXIS_INDEX[axis], n_sectors)
    xyz = np.asarray(skeleton.getCoordinates().getValues(), dtype=float)
    rows, err = node_rows(xyz, coords)
    if err:
        return None, err
    order = "zelfde volgorde" if (len(rows) == len(xyz) // 3 and np.all(
        rows == np.arange(len(rows)))) else "op coordinaat gekoppeld"
    print("aster_cyclic: hele wiel: {0} sectoren, {1} knopen ({2}).".format(
        n_sectors, len(rows), order))
    return _modes_on(wheel, skeleton, table, axis), None


def collect(modsec, mesh, numbering, n_sectors, nmax_freq, diameters,
            axis="Z", angle=0.0, wheel=None, skeleton=None, liml=None):
    """
    MODE_MECA met de modes in Mecway's volgorde (mode_table), in het
    assenstelsel van het model: op het hele wiel als wheel (REST_SOUS_STRUC
    SQUELETTE), skeleton (het skelet) en liml (de .liml van het model)
    gegeven zijn en het wiel klopt (_wheel), anders op de sector (mesh).
    numbering: de NUME_DDL van de sectorbasis (_attach). axis: "X", "Y" of
    "Z" (de symmetrieas in Mecway); bij X/Y draaien mesh en skelet terug
    (MODI_MAILLAGE ROTATION om (1, 1, 1) met -angle graden) en worden de
    componenten van DEPL teruggezet (_BACK). Drukt per mode
    BRIDGE_MODE <nr> FREQ=<Hz> ND=<k> (bridge_results_liml.read_mode_summary)
    en BRIDGE_CYCL N=.. NMAX_FREQ=.. MODES=.. WIEL=0|1 af.
    """
    _attach(modsec, mesh, numbering)
    acc = modsec.getAccessParameters()
    freqs = list(acc["FREQ"])
    table = mode_table(freqs, n_sectors, nmax_freq, diameters)
    if table is None:
        print(tr("aster_cyclic: WAARSCHUWING: {0} modes uit "
                 "MODE_ITER_CYCL, verwacht {1}; de nodal diameters zijn "
                 "onbekend en alle modes gaan op volgorde naar Mecway.",
                 "aster_cyclic: WARNING: {0} modes from MODE_ITER_CYCL, "
                 "expected {1}; the nodal diameters are unknown and all "
                 "modes go to Mecway in order.").format(
                  len(freqs), sum(nmax_freq * _multiplicity(d, n_sectors)
                                  for d in diameters)))
        table = [(float(f), k, -1) for k, f in enumerate(freqs)]
    res = None
    if wheel is not None and skeleton is not None and liml:
        try:
            res, why = _wheel(wheel, skeleton, table, axis, angle,
                              n_sectors, liml)
        except Exception as e:
            res, why = None, repr(e)
        if res is None:
            print(tr("aster_cyclic: WAARSCHUWING: het hele wiel is niet "
                     "opgebouwd ({0}); naar Mecway gaat alleen de sector.",
                     "aster_cyclic: WARNING: the whole wheel was not built "
                     "({0}); only the sector goes to Mecway.").format(why))
    full = res is not None
    if res is None:
        _rotate_back(mesh, axis, angle)
        res = _modes_on(modsec, mesh, table, axis)
    for k, (f, _pos, nd) in enumerate(table):
        print("BRIDGE_MODE {0} FREQ={1!r} ND={2}".format(k + 1, f, nd))
    print("BRIDGE_CYCL N={0} NMAX_FREQ={1} MODES={2} WIEL={3}".format(
        n_sectors, nmax_freq, len(table), 1 if full else 0))
    return res
