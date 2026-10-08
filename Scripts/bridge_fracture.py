# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# bridge_fracture.py
# Breukmechanica met een gemeshte scheur (STAPPENPLAN S6a): K_I/K_II/K_III
# en G langs het scheurfront als nabewerking van een lineair statische
# analyse (lin_static_3d, MECA_STATIQUE). Gebruikersdoc:
# docs/Breukmechanica.md (EN: docs/en/Fracture_mechanics.md).
# Compatibel met IronPython 2.7.9 (Mecway); de .comm-code die
# fracture_block() schrijft draait in Code-Aster (Python 3).
#
# Invoer (Mecway heeft geen breukmechanica, manual v33; wel *Disconnect
# elements*, par. 4.2.10, en *Merge nearby nodes* met een knoopselectie,
# par. 4.2.11, waarmee je een scheur opent):
#   knoopselectie  SCHEUR_FRONT     (of CRACK_FRONT)      knopen van het front
#   vlakselectie   SCHEUR_LIP_BOVEN (of CRACK_LIP_UPPER)  bovenlip
#   vlakselectie   SCHEUR_LIP_ONDER (of CRACK_LIP_LOWER)  onderlip; zonder
#                  onderlip een half model met symmetrie (SYME='OUI')
# Meer scheuren: hetzelfde achtervoegsel achter alle drie de namen
# (SCHEUR_FRONT_A, SCHEUR_LIP_BOVEN_A, ...). Herkend op de oorspronkelijke
# Mecway-naam, hoofdletterongevoelig.
#
# Vertaling (catalogus Code-Aster 17.4, Cata/Commands):
#   - bridge_mesh.write_mesh schrijft de ribben van de volume-elementen
#     langs het front als SEG2/SEG3-cellen in groep BRSCHEUR<n>
#     (front_cells); DEFI_FOND_FISS FOND_FISS kent alleen GROUP_MA of een
#     geordende GROUP_NO (defi_fond_fiss.py), en Mecway kent geen ribben.
#   - na MECA_STATIQUE en de resultaatexport (fracture_block):
#     DEFI_FOND_FISS (CONFIG_INIT='COLLEE'), CALC_G OPTION=('G', 'K')
#     (calc_g.py; THETA R_INF/R_SUP uit fracture_r_inf/_sup of automatisch)
#     en POST_K1_K2_K3 (post_k1_k2_k3.py), tabellen in de .mess en als
#     <model>_<run_id>_scheur.csv; report() zet een samenvatting in de log en
#     het statusvak.

import io
import math
import os
from collections import OrderedDict

from bridge_utils import (log, tr, run_name, CONFIG_BRIDGE,
                          get_active_analysis, work_dir_path)

# Namen van de selecties (prefix, zonder achtervoegsel)
FRONT_NAMES = ("SCHEUR_FRONT", "CRACK_FRONT")
UPPER_NAMES = ("SCHEUR_LIP_BOVEN", "CRACK_LIP_UPPER")
LOWER_NAMES = ("SCHEUR_LIP_ONDER", "CRACK_LIP_LOWER")
# Groep met de SEG-cellen van het front (+ volgnummer); gereserveerd in
# bridge_mesh._aster_names.
FRONT_GROUP = "BRSCHEUR"
# Achtervoegsel van de tabel in de werkmap
CSV_SUFFIX = "_scheur.csv"

# Ribben per volumetype (hoekindices in Code-Aster-volgorde; de k-de rib
# heeft bij een kwadratisch element middenknoop n_hoeken + k), zoals
# bridge_mesh._MID_EDGES.
_EDGES = {
    "TETRA": [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)],
    "PYRAM": [(0, 1), (1, 2), (2, 3), (0, 3), (0, 4), (1, 4), (2, 4),
              (3, 4)],
    "PENTA": [(0, 1), (1, 2), (0, 2), (0, 3), (1, 4), (2, 5), (3, 4),
              (4, 5), (3, 5)],
    "HEXA":  [(0, 1), (1, 2), (2, 3), (0, 3), (0, 4), (1, 5), (2, 6),
              (3, 7), (4, 5), (5, 6), (6, 7), (4, 7)],
}
_CORNERS = {"TETRA": 4, "PYRAM": 5, "PENTA": 6, "HEXA": 8}

# Kolommen van de tabel (CSV): (kop, tabel, Code-Aster-parameter).
# "g" = CALC_G OPTION=('G', 'K') (G: theta-methode; K1..K3:
# interactie-integraal; G_IRWIN = (K1^2 + K2^2)(1 - nu^2)/E + K3^2 (1 + nu)
# /E), "pk" = POST_K1_K2_K3 (extrapolatie van de sprong in de
# verplaatsing over de lippen; zijn G is ook uit K berekend).
CSV_COLUMNS = (
    ("crack", None, None), ("point", "g", "NUM_PT"), ("node", "g", "NOEUD"),
    ("x", "g", "COOR_X"), ("y", "g", "COOR_Y"), ("z", "g", "COOR_Z"),
    ("s", "g", "ABSC_CURV"), ("G", "g", "G"), ("K1", "g", "K1"),
    ("K2", "g", "K2"), ("K3", "g", "K3"), ("G_IRWIN", "g", "G_IRWIN"),
    ("K1_POST_K", "pk", "K1"), ("K2_POST_K", "pk", "K2"),
    ("K3_POST_K", "pk", "K3"), ("G_POST_K", "pk", "G"),
)


# ==========================================
# HERKENNEN
# ==========================================
def _match(name, prefixes):
    """Achtervoegsel (hoofdletters) als name met een van de prefixen
    begint, anders None."""
    up = (name or u"").upper()
    for p in prefixes:
        if up.startswith(p):
            return up[len(p):]
    return None


def find_cracks(model_data):
    """
    Scheuren uit de Mecway-selecties: lijst van dicts {"key": achtervoegsel,
    "front": knoopselectie, "upper": vlakselectie, "lower": vlakselectie of
    None} (oorspronkelijke Mecway-namen), gesorteerd op achtervoegsel.
    Onvolledige invoer: ERROR en None. Geen scheur: [].
    """
    nodes = model_data.get("node_groups", {})
    faces = model_data.get("face_selection_faces", {})
    elems = model_data.get("element_selections", {})
    found = OrderedDict()
    ok = True
    for kind, prefixes in (("front", FRONT_NAMES), ("upper", UPPER_NAMES),
                           ("lower", LOWER_NAMES)):
        # Een knoopselectie als lip of een vlakselectie als front is een
        # vergissing; melden met de juiste soort.
        for name in list(nodes) + list(faces) + list(elems):
            key = _match(name, prefixes)
            if key is None:
                continue
            want = nodes if kind == "front" else faces
            if name not in want:
                log(tr("ERROR: Breukmechanica: '{0}' moet een {1} zijn "
                       "(in Mecway: {2}).",
                       "ERROR: Fracture mechanics: '{0}' must be a {1} "
                       "(in Mecway: {2}).").format(
                    name, tr("knoopselectie", "node selection")
                    if kind == "front" else
                    tr("vlakselectie", "face selection"),
                    tr("selecteer de knopen van het front en kies Named "
                       "selection", "select the front nodes and choose "
                       "Named selection") if kind == "front" else
                    tr("selecteer de vlakken van de lip (Select faces) en "
                       "kies Named selection", "select the faces of the lip "
                       "(Select faces) and choose Named selection")))
                ok = False
                continue
            crack = found.setdefault(key, {"key": key, "front": None,
                                           "upper": None, "lower": None})
            if crack[kind] is not None and crack[kind] != name:
                log(tr("ERROR: Breukmechanica: twee selecties voor dezelfde "
                       "scheur: '{0}' en '{1}'. Geef elke scheur een eigen "
                       "achtervoegsel (bijv. SCHEUR_FRONT_A).",
                       "ERROR: Fracture mechanics: two selections for the "
                       "same crack: '{0}' and '{1}'. Give each crack its own "
                       "suffix (e.g. CRACK_FRONT_A).").format(
                           crack[kind], name))
                ok = False
            crack[kind] = name
    out = []
    for key in sorted(found):
        c = found[key]
        if c["front"] is None or c["upper"] is None:
            log(tr("ERROR: Breukmechanica: scheur '{0}' is onvolledig "
                   "(front: {1}, bovenlip: {2}, onderlip: {3}). Nodig: een "
                   "knoopselectie SCHEUR_FRONT{0} en een vlakselectie "
                   "SCHEUR_LIP_BOVEN{0}; SCHEUR_LIP_ONDER{0} voor een "
                   "scheur met twee lippen.",
                   "ERROR: Fracture mechanics: crack '{0}' is incomplete "
                   "(front: {1}, upper lip: {2}, lower lip: {3}). Needed: a "
                   "node selection CRACK_FRONT{0} and a face selection "
                   "CRACK_LIP_UPPER{0}; CRACK_LIP_LOWER{0} for a crack with "
                   "two lips.").format(c["key"], c["front"] or "-",
                                       c["upper"] or "-", c["lower"] or "-"))
            ok = False
            continue
        out.append(c)
    return out if ok else None


def is_requested(model_data):
    """True als het model scheur-selecties heeft (ook onvolledige)."""
    names = (list(model_data.get("node_groups", {})) +
             list(model_data.get("face_selection_faces", {})))
    return any(_match(n, FRONT_NAMES + UPPER_NAMES + LOWER_NAMES)
               is not None for n in names)


def active_here():
    """True als de actieve analyse lineair statisch (MECA_STATIQUE) is:
    alleen daar bouwt de bridge de breukmechanica."""
    try:
        return get_active_analysis().get("solver_command") == "MECA_STATIQUE"
    except KeyError:
        return False


# ==========================================
# MESH: FRONT ALS SEG-CELLEN
# ==========================================
def _nid(value):
    return int(str(value).lstrip("N"))


def _chain(edges):
    """
    Orden de ribben [(a, b, mid)] tot een keten. Returns (geordende
    ribben, gesloten) of (None, reden) bij een vertakking of meerdere
    stukken.
    """
    adj = {}
    for k, (a, b, _m) in enumerate(edges):
        adj.setdefault(a, []).append(k)
        adj.setdefault(b, []).append(k)
    if any(len(v) > 2 for v in adj.values()):
        return None, "branch"
    ends = [n for n, v in adj.items() if len(v) == 1]
    closed = not ends
    start = min(adj) if closed else min(ends)
    used = set()
    out = []
    node = start
    while True:
        nxt = [k for k in adj[node] if k not in used]
        if not nxt:
            break
        k = nxt[0]
        used.add(k)
        a, b, m = edges[k]
        if a != node:
            a, b = b, a
        out.append((a, b, m))
        node = b
        if closed and node == start:
            break
    if len(out) != len(edges):
        return None, "pieces"
    return out, closed


def front_cells(model_data, conn):
    """
    SEG2/SEG3-cellen langs het scheurfront per scheur, voor
    bridge_mesh.write_mesh. Een rib van een volume-element hoort bij het
    front als zijn twee hoekknopen in SCHEUR_FRONT zitten (bij een
    kwadratisch element komt de middenknoop mee, ook als die niet
    geselecteerd is). Vult model_data["cracks"] (met "front_group",
    "closed", "front_xyz", "quadratic").

    Returns OrderedDict {groep: [(type, [knopen])]} of None bij een fout
    (ERROR in de log). Geen scheur of geen lineair statische analyse: {}.
    """
    model_data["cracks"] = []
    if not is_requested(model_data):
        return OrderedDict()
    if not active_here():
        # generate_AsterSolver_files meldt het (WARNING) voor de analyse
        return OrderedDict()
    cracks = find_cracks(model_data)
    if cracks is None:
        return None
    nodes = model_data.get("mesh_nodes", {})
    out = OrderedDict()
    for k, crack in enumerate(cracks):
        front = set(_nid(n) for n in
                    model_data["node_groups"].get(crack["front"], []))
        edges = OrderedDict()
        front_eids = set()
        for eid in sorted(conn):
            atype, cn = conn[eid]
            base = None
            for prefix in _EDGES:
                if atype.startswith(prefix):
                    base = prefix
            if base is None:
                continue
            nc = _CORNERS[base]
            quad = len(cn) > nc
            for j, (p, q) in enumerate(_EDGES[base]):
                a, b = cn[p], cn[q]
                if a in front and b in front:
                    front_eids.add(eid)
                    key = frozenset((a, b))
                    if key not in edges:
                        edges[key] = (a, b, cn[nc + j] if quad else None)
        if not edges:
            log(tr("ERROR: Breukmechanica: geen rib van een volume-element "
                   "met beide eindknopen in '{0}'. Selecteer alle knopen "
                   "langs het front (minstens twee naast elkaar op een "
                   "elementrand).",
                   "ERROR: Fracture mechanics: no edge of a solid element "
                   "with both end nodes in '{0}'. Select all nodes along "
                   "the front (at least two neighbours on an element "
                   "edge).").format(crack["front"]))
            return None
        chain, closed = _chain(list(edges.values()))
        if chain is None:
            log(tr("ERROR: Breukmechanica: de knopen van '{0}' vormen geen "
                   "enkele lijn ({1}). Een front is een open of gesloten "
                   "lijn van elementranden; gebruik voor een tweede front "
                   "een eigen achtervoegsel.",
                   "ERROR: Fracture mechanics: the nodes of '{0}' do not "
                   "form a single line ({1}). A front is an open or closed "
                   "line of element edges; use its own suffix for a second "
                   "front.").format(
                       crack["front"],
                       tr("vertakking", "branch") if closed == "branch"
                       else tr("meerdere stukken", "several pieces")))
            return None
        on_line = set()
        for a, b, m in chain:
            on_line.update(n for n in (a, b, m) if n is not None)
        loose = sorted(front - on_line)
        if loose:
            log(tr("WARNING: Breukmechanica: {0} knoop/knopen van '{1}' "
                   "liggen niet op een rib langs het front en tellen niet "
                   "mee (bijv. N{2}).",
                   "WARNING: Fracture mechanics: {0} node(s) of '{1}' are "
                   "not on an edge along the front and are ignored (e.g. "
                   "N{2}).").format(len(loose), crack["front"], loose[0]))
        quadratic = all(m is not None for _a, _b, m in chain)
        group = "{0}{1}".format(FRONT_GROUP, k + 1)
        out[group] = [("SEG3", [a, b, m]) if m is not None else
                      ("SEG2", [a, b]) for a, b, m in chain]
        pts = [nodes[chain[0][0]]] + [nodes[b] for _a, b, _m in chain]
        crack.update({"front_group": group, "closed": closed,
                      "front_xyz": pts, "quadratic": quadratic,
                      "n_edges": len(chain), "front_eids": front_eids})
        if not _check_lips(model_data, crack, conn, on_line):
            return None
        log("INFO: Breukmechanica: scheur {0}: front {1} ({2} ribben, "
            "{3}, {4}), lippen {5}{6}.".format(
                k + 1, crack["front"], len(chain),
                "gesloten" if closed else "open",
                "kwadratisch" if quadratic else "lineair", crack["upper"],
                " + " + crack["lower"] if crack["lower"] else
                " (half model, SYME='OUI')"))
        model_data["cracks"].append(crack)
    return out


def _face_nodes(model_data, conn, fname):
    """Alle knopen van de vlakken van vlakselectie fname (hoek- en
    middenknopen, uit de tabel van bridge_mesh), en de soorten elementen."""
    import bridge_mesh
    shapes = model_data.get("elem_shapes", {})
    out = set()
    kinds = set()
    for eid, fid in model_data.get("face_selection_faces", {}).get(fname, []):
        if eid not in conn:
            continue
        kinds.add("solid" if shapes.get(eid) in bridge_mesh._SOLID_SHAPES
                  else "other")
        atype, cn = conn[eid]
        corners = bridge_mesh._face_from_table(atype, cn, fid)
        if corners:
            stype, snodes = bridge_mesh._skin_element(atype, cn, corners)
            out.update(snodes)
    return out, kinds


def _check_lips(model_data, crack, conn, front_nodes):
    """
    Lippen: vlakken van volume-elementen; bij twee lippen geen gedeelde
    knopen behalve op het front (anders is de scheur dicht: Disconnect
    elements vergeten). Fout: ERROR en False.
    """
    upper, kinds_u = _face_nodes(model_data, conn, crack["upper"])
    lower, kinds_l = (_face_nodes(model_data, conn, crack["lower"])
                      if crack["lower"] else (set(), set()))
    if "other" in (kinds_u | kinds_l) or not upper:
        log(tr("ERROR: Breukmechanica: de lippen van scheur '{0}' moeten "
               "vlakken van volume-elementen zijn (geen schalen of balken).",
               "ERROR: Fracture mechanics: the lips of crack '{0}' must be "
               "faces of solid elements (no shells or beams).").format(
                   crack["front"]))
        return False
    if not (upper & front_nodes):
        log(tr("WARNING: Breukmechanica: bovenlip '{0}' raakt het front "
               "'{1}' niet; controleer de selecties.",
               "WARNING: Fracture mechanics: upper lip '{0}' does not touch "
               "the front '{1}'; check the selections.").format(
                   crack["upper"], crack["front"]))
    if crack["lower"]:
        shared = (upper & lower) - front_nodes
        if shared:
            log(tr("ERROR: Breukmechanica: boven- en onderlip van scheur "
                   "'{0}' delen {1} knopen buiten het front (bijv. N{2}): "
                   "de scheur is niet open. Maak hem open in Mecway: "
                   "selecteer de elementen aan een kant van het scheurvlak, "
                   "Mesh tools > Disconnect elements, en voeg daarna alleen "
                   "de knopen buiten de scheur (ligament en front) weer "
                   "samen met Merge nearby nodes (docs/Breukmechanica.md).",
                   "ERROR: Fracture mechanics: upper and lower lip of crack "
                   "'{0}' share {1} nodes outside the front (e.g. N{2}): "
                   "the crack is not open. Open it in Mecway: select the "
                   "elements on one side of the crack plane, Mesh tools > "
                   "Disconnect elements, then merge only the nodes outside "
                   "the crack (ligament and front) again with Merge nearby "
                   "nodes (docs/en/Fracture_mechanics.md).").format(
                       crack["front"], len(shared), sorted(shared)[0]))
            return False
        if not (upper & lower):
            log(tr("WARNING: Breukmechanica: de lippen van scheur '{0}' "
                   "delen geen frontknopen; is het front ook losgemaakt? "
                   "Het front moet aan beide lippen vastzitten.",
                   "WARNING: Fracture mechanics: the lips of crack '{0}' "
                   "share no front nodes; was the front disconnected too? "
                   "The front must be attached to both lips.").format(
                       crack["front"]))
    return True


# ==========================================
# GEOMETRIE VAN HET FRONT
# ==========================================
def _sub(a, b):
    return [a[i] - b[i] for i in range(3)]


def _dot(a, b):
    return sum(a[i] * b[i] for i in range(3))


def _cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]]


def _norm(a):
    return math.sqrt(_dot(a, a))


def front_shape(pts, tol=0.01):
    """
    Vorm van het front: ("recht", None), ("cirkel", (straal, middelpunt))
    of ("krom", None). Recht: alle punten binnen tol x lengte van de lijn
    door de eindpunten (bij een gesloten front nooit). Cirkel: in een vlak
    en op een straal binnen tol (middelpunt van de omgeschreven cirkel van
    het eerste, middelste en laatste punt; bij een gesloten front het
    eerste, een derde en twee derde).
    """
    if len(pts) < 3:
        return "recht", None
    size = max(_norm(_sub(p, pts[0])) for p in pts)
    if size <= 0.0:
        return "krom", None
    a, b = pts[0], pts[-1]
    ab = _sub(b, a)
    lab = _norm(ab)
    if lab > tol * size:
        u = [c / lab for c in ab]
        dev = max(_norm(_cross(_sub(p, a), u)) for p in pts)
        if dev <= tol * size:
            return "recht", None
    n = len(pts)
    if lab <= tol * size:     # gesloten
        p1, p2, p3 = pts[0], pts[n // 3], pts[(2 * n) // 3]
    else:
        p1, p2, p3 = pts[0], pts[n // 2], pts[-1]
    # Omgeschreven cirkel van p1, p2, p3
    u, v = _sub(p2, p1), _sub(p3, p1)
    w = _cross(u, v)
    ww = _dot(w, w)
    if ww <= (tol * size) ** 4:
        return "krom", None
    t1 = [c * _dot(v, v) for c in _cross(w, u)]
    t2 = [c * _dot(u, u) for c in _cross(v, w)]
    centre = [p1[i] + (t1[i] + t2[i]) / (2.0 * ww) for i in range(3)]
    radius = _norm(_sub(p1, centre))
    wn = [c / math.sqrt(ww) for c in w]
    for p in pts:
        if abs(_norm(_sub(p, centre)) - radius) > tol * radius or \
                abs(_dot(_sub(p, centre), wn)) > tol * radius:
            return "krom", None
    return "cirkel", (radius, centre)


# ==========================================
# .COMM
# ==========================================
# Minimum aantal knopen langs het front voor DISCRETISATION='LINEAIRE':
# met 3 knopen (plaat met 2 elementen over de dikte) oscilleerde G per knoop
# (11739 / 4696 / 25826 J/m2 waar 9391 verwacht; hex8 en hex20 met
# NB_POINT_FOND=3), met 5 knopen en meer niet (gemeten bij de bouw,
# docs/Breukmechanica.md par. 3).
_MIN_POINTS_LINEAIRE = 5


def _n_points(crack, line_quad=False):
    """Aantal knopen langs het front (hoek- en middenknopen; met LINE_QUAD
    krijgt elke rib een middenknoop)."""
    quad = crack.get("quadratic") or line_quad
    n = crack.get("n_edges", 0) * (2 if quad else 1)
    return n if crack.get("closed") else n + 1


def _smoothing(analysis, crack, line_quad=False):
    """THETA-tekst voor DISCRETISATION ("" = LINEAIRE, de standaard van
    CALC_G, catalogus calc_g.py), met een WARNING bij te weinig punten of
    een te hoge graad. None bij een onbekende keuze (ERROR)."""
    choice = str(analysis.get("fracture_smoothing", "LINEAIRE")).upper()
    npts = _n_points(crack, line_quad)
    if choice == "LINEAIRE":
        if npts < _MIN_POINTS_LINEAIRE:
            log(tr("WARNING: Breukmechanica: het front van scheur '{0}' heeft "
                   "maar {1} knopen; met fracture_smoothing LINEAIRE "
                   "oscilleren G en K dan van knoop tot knoop. Verfijn het "
                   "mesh langs het front (minstens 4 elementen) of kies "
                   "LEGENDRE met graad 0 (een constante waarde).",
                   "WARNING: Fracture mechanics: the front of crack '{0}' "
                   "has only {1} nodes; with fracture_smoothing LINEAIRE, G "
                   "and K then oscillate from node to node. Refine the mesh "
                   "along the front (at least 4 elements) or choose LEGENDRE "
                   "with degree 0 (a constant value).").format(
                       crack["front"], npts))
        return ""
    if choice == "LEGENDRE":
        degree = int(analysis.get("fracture_legendre_degree", 5))
        if degree < 0 or degree > 7:
            log(tr("ERROR: Breukmechanica: fracture_legendre_degree {0}: "
                   "toegestaan 0 .. 7 (CALC_G DEGRE).",
                   "ERROR: Fracture mechanics: fracture_legendre_degree {0}: "
                   "allowed 0 .. 7 (CALC_G DEGRE).").format(degree))
            return None
        if degree >= npts:
            log(tr("WARNING: Breukmechanica: Legendre-graad {0} bij {1} "
                   "knopen langs het front van scheur '{2}': de waarden "
                   "zijn dan onbetrouwbaar (gemeten: graad 5 bij 3 knopen "
                   "-8 tot -22 %). Kies een graad kleiner dan het aantal "
                   "knopen.",
                   "WARNING: Fracture mechanics: Legendre degree {0} with "
                   "{1} nodes along the front of crack '{2}': the values "
                   "are then unreliable (measured: degree 5 with 3 nodes "
                   "-8 to -22 %). Choose a degree below the number of "
                   "nodes.").format(degree, npts, crack["front"]))
        return ", DISCRETISATION='LEGENDRE', DEGRE={0}".format(degree)
    log(tr("ERROR: Breukmechanica: fracture_smoothing '{0}' onbekend; kies "
           "LINEAIRE of LEGENDRE.",
           "ERROR: Fracture mechanics: fracture_smoothing '{0}' unknown; "
           "choose LINEAIRE or LEGENDRE.").format(choice))
    return None


def _r_bounds(analysis):
    """(R_INF, R_SUP) in m uit de instellingen, of None (automatisch).
    Onzinnige combinatie: ERROR en False."""
    r_inf = float(analysis.get("fracture_r_inf", 0.0) or 0.0)
    r_sup = float(analysis.get("fracture_r_sup", 0.0) or 0.0)
    if r_inf == 0.0 and r_sup == 0.0:
        return None
    if r_inf < 0.0 or r_sup <= r_inf:
        log(tr("ERROR: Breukmechanica: fracture_r_inf ({0} m) en "
               "fracture_r_sup ({1} m): nodig 0 <= R_INF < R_SUP, of beide 0 "
               "(automatisch).",
               "ERROR: Fracture mechanics: fracture_r_inf ({0} m) and "
               "fracture_r_sup ({1} m): need 0 <= R_INF < R_SUP, or both 0 "
               "(automatic).").format(r_inf, r_sup))
        return False
    return r_inf, r_sup


def _py(text):
    return "'{0}'".format(text.replace("\\", "\\\\").replace("'", "\\'"))


# CSV-schrijver in de job (Python 3 van Code-Aster). EXTR_TABLE().values()
# geeft per parameter een lijst (None waar leeg); de rijen van CALC_G en
# POST_K1_K2_K3 worden op knoopnaam (NOEUD) gekoppeld.
_CSV_WRITER = '''
def _bridge_scheur_csv(path, cracks, columns):
    import io
    lines = [",".join(c[0] for c in columns)]
    for key, t_g, t_pk in cracks:
        vg = t_g.EXTR_TABLE().values()
        vpk = {}
        if t_pk is not None:
            vp = t_pk.EXTR_TABLE().values()
            for i, nd in enumerate(vp.get("NOEUD", [])):
                vpk[str(nd).strip()] = dict(
                    (k, vp[k][i]) for k in vp if i < len(vp[k]))
        for i, nd in enumerate(vg.get("NOEUD", [])):
            row = []
            for head, tab, para in columns:
                if tab is None:
                    v = key
                elif tab == "g":
                    v = vg.get(para, [None] * (i + 1))[i]
                else:
                    v = vpk.get(str(nd).strip(), {}).get(para)
                if v is None:
                    row.append("")
                elif isinstance(v, float):
                    row.append("%.10g" % v)
                else:
                    row.append(str(v).strip())
            lines.append(",".join(row))
    with io.open(path, "w", encoding="utf-8") as f:
        f.write("\\n".join(lines) + "\\n")
    print("BRIDGE_SCHEUR_CSV " + path)
'''


def _material_block(model_data, concepts, excit_clause, cara_ref):
    """
    Eigen materiaalveld en resultaat voor CALC_G/POST_K1_K2_K3. Waarom:
    CALC_G rekent ook op de randcellen (lippen, belaste vlakken) in het
    theta-gebied en vraagt daar een materiaal (<EXCEPTION> CALCUL_20 "Le
    matériau est nécessaire sur la maille", option CALC_K_G, MECA_FACE8;
    gezien bij de bouw), terwijl de bridge het materiaal per elset geeft.
    Hier krijgt eerst alles (TOUT='OUI') het materiaal van de elementen aan
    het front, daarna elke elset zijn eigen materiaal (later wint); de
    verplaatsingen en lasten zijn die van res (CREA_RESU AFFE met EXCIT,
    catalogus crea_resu.py b_evol_elas). Returns (tekst, materiaalnamen aan
    het front).
    """
    elem_elset = model_data.get("_elem_elset", {})
    front_mats = []
    for c in model_data.get("cracks") or []:
        for eid in sorted(c.get("front_eids", ())):
            mat = model_data.get("elsets", {}).get(elem_elset.get(eid))
            if mat in concepts and mat not in front_mats:
                front_mats.append(mat)
    entries, varc = model_data.get("_fieldmat_parts") or ([], "")
    if not front_mats:
        return None, front_mats
    first = "        _F(TOUT='OUI', MATER={0}),\n".format(
        concepts[front_mats[0]])
    text = (
        "fieldmat_fiss = AFFE_MATERIAU(\n"
        "    MAILLAGE=mesh,\n"
        "    MODELE=model,\n"
        "    AFFE=(\n" + first + "".join(entries) + "    ),\n" + varc + ")\n"
        "u_fiss = CREA_CHAMP(TYPE_CHAM='NOEU_DEPL_R', OPERATION='EXTR',\n"
        "                    RESULTAT=res, NOM_CHAM='DEPL', NUME_ORDRE=1)\n"
        "res_fiss = CREA_RESU(\n"
        "    OPERATION='AFFE',\n"
        "    TYPE_RESU='EVOL_ELAS',\n"
        "    AFFE=_F(NOM_CHAM='DEPL', CHAM_GD=u_fiss, MODELE=model,\n"
        "            CHAM_MATER=fieldmat_fiss,{0} INST=0.0),\n"
        "{1}"
        ")\n"
        "# CALC_G OPTION 'G' leest SIEF_ELGA (<EXCEPTION> RUPTURE0_94 zonder)\n"
        "res_fiss = CALC_CHAMP(reuse=res_fiss, RESULTAT=res_fiss,\n"
        "                      CONTRAINTE=('SIEF_ELGA',))\n").format(
            " CARA_ELEM=cara," if cara_ref else "", excit_clause)
    return text, front_mats


def fracture_block(model_data, analysis, work_dir, concepts=None,
                   excit_clause="", cara_ref=""):
    """
    .comm-code na MECA_STATIQUE en de resultaatexport: per scheur
    DEFI_FOND_FISS, CALC_G en POST_K1_K2_K3 op een eigen resultaat
    (_material_block), de tabellen in de .mess en de CSV (csv_path). Na de
    export, zodat een fout hier de resultaten in Mecway niet tegenhoudt.
    concepts: materiaalnaam -> DEFI_MATERIAU-concept (MAT_n); excit_clause
    en cara_ref zoals in MECA_STATIQUE. "" zonder scheur; None bij een
    fout.
    """
    cracks = model_data.get("cracks") or []
    if not cracks:
        return ""
    bounds = _r_bounds(analysis)
    if bounds is False:
        return None
    mat_text, front_mats = _material_block(model_data, concepts or {},
                                           excit_clause, cara_ref)
    if mat_text is None:
        log(tr("ERROR: Breukmechanica: de elementen aan het scheurfront "
               "hebben geen (mechanisch) materiaal.",
               "ERROR: Fracture mechanics: the elements at the crack front "
               "have no (mechanical) material."))
        return None
    if len(front_mats) > 1:
        log(tr("WARNING: Breukmechanica: het front raakt elementen met "
               "verschillende materialen ({0}); een scheur op een grensvlak "
               "is NIET geverifieerd en POST_K1_K2_K3 vraagt een homogeen "
               "materiaal. De randcellen krijgen materiaal '{1}'.",
               "WARNING: Fracture mechanics: the front touches elements "
               "with different materials ({0}); a crack on an interface is "
               "NOT verified and POST_K1_K2_K3 needs a homogeneous "
               "material. The boundary cells get material '{1}'.").format(
                   ", ".join(front_mats), front_mats[0]))
    rn = model_data.get("selection_renames", {})
    theta_r = ("" if bounds is None else
               ", R_INF={0!r}, R_SUP={1!r}".format(bounds[0], bounds[1]))
    lines = ["", "# " + "=" * 53,
             "# BREUKMECHANICA (bridge: SCHEUR_-selecties, bridge_fracture.py,",
             "# docs/Breukmechanica.md): K en G langs het scheurfront",
             "# " + "=" * 53, mat_text]
    entries = []
    for k, c in enumerate(cracks):
        i = k + 1
        smooth = _smoothing(analysis, c, bool(model_data.get("line_quad")))
        if smooth is None:
            return None
        upper = rn.get(c["upper"], c["upper"])
        lower = rn.get(c["lower"], c["lower"]) if c["lower"] else None
        shape, info = front_shape(c["front_xyz"])
        form = ""
        if shape == "cirkel" and _disc_crack(model_data, c, info):
            # Krommingscorrectie van de interactie-integraal (calc_g.py
            # b_courbure; forma07a zonder: K1 +5,6 %, met: +0,1 %, gemeten
            # bij de bouw)
            form = ", FORM_FISS='CERCLE', RAYON={0!r}".format(info[0])
        c["front_shape"] = shape
        if shape == "krom" or c.get("front_shape_note") == "cirkel_buiten":
            log(tr("WARNING: Breukmechanica: het front van scheur '{0}' is "
                   "geen rechte lijn en geen schijfvormige cirkel: K uit "
                   "CALC_G (interactie-integraal) heeft dan geen "
                   "krommingscorrectie en kan enkele procenten afwijken. "
                   "Gebruik G (CALC_G) of K uit POST_K1_K2_K3.",
                   "WARNING: Fracture mechanics: the front of crack '{0}' "
                   "is not a straight line nor a disc-shaped circle: K from "
                   "CALC_G (interaction integral) then has no curvature "
                   "correction and can be a few percent off. Use G (CALC_G) "
                   "or K from POST_K1_K2_K3.").format(c["front"]))
        lines += [
            "# Scheur {0}: front {1} ({2}), lippen {3}{4}".format(
                i, c["front"], shape, c["upper"],
                " / " + c["lower"] if lower else " (symmetrie)"),
            "fond_{0} = DEFI_FOND_FISS(".format(i),
            "    MAILLAGE=mesh,",
            "    FOND_FISS=_F(GROUP_MA='{0}'{1}),".format(
                c["front_group"], ", TYPE_FOND='FERME'" if c["closed"] else ""),
            "    SYME='{0}',".format("NON" if lower else "OUI"),
            "    LEVRE_SUP=_F(GROUP_MA='{0}'),".format(upper),
        ]
        if lower:
            lines.append("    LEVRE_INF=_F(GROUP_MA='{0}'),".format(lower))
        lines += [
            "    CONFIG_INIT='COLLEE',",
            ")",
            "scheur_g_{0} = CALC_G(".format(i),
            "    RESULTAT=res_fiss,",
            "    OPTION=('G', 'K'),",
            "    THETA=_F(FISSURE=fond_{0}{1}{2}){3},".format(
                i, theta_r, smooth, form),
            ")",
            "IMPR_TABLE(TABLE=scheur_g_{0}, UNITE=6)".format(i),
            "# POST_K1_K2_K3 (extrapolatie van de lipverplaatsingen) als",
            "# controle; een fout hier stopt de tabel van CALC_G niet",
            "try:",
            "    scheur_pk_{0} = POST_K1_K2_K3(RESULTAT=res_fiss, "
            "FOND_FISS=fond_{0})".format(i),
            "    IMPR_TABLE(TABLE=scheur_pk_{0}, UNITE=6)".format(i),
            "except Exception as _e:",
            "    scheur_pk_{0} = None".format(i),
            "    print('BRIDGE_SCHEUR_POSTK_FOUT {0} ' + str(_e)"
            ".replace(chr(10), ' ')[:300])".format(i),
        ]
        entries.append("({0}, scheur_g_{1}, scheur_pk_{1})".format(
            _py(c["key"] or str(i)), i))
    lines.append(_CSV_WRITER)
    lines.append("_bridge_scheur_csv({0}, [{1}], {2!r})".format(
        _py(csv_path(work_dir)), ", ".join(entries),
        [list(c) for c in CSV_COLUMNS]))
    return "\n".join(lines) + "\n"


def _disc_crack(model_data, crack, info):
    """True als de lip binnen de cirkel van het front ligt (een
    schijfvormige scheur: FORM_FISS='CERCLE' van CALC_G). Ligt de lip
    erbuiten (bijv. een omtrekscheur van buiten), dan geen correctie
    (WARNING in report via front_shape 'cirkel_buiten')."""
    radius, centre = info
    nodes = model_data.get("mesh_nodes", {})
    conn = model_data.get("_conn") or {}
    lip, _k = _face_nodes(model_data, conn, crack["upper"])
    if not lip:
        return False
    inside = sum(1 for n in lip
                 if _norm(_sub(nodes[n], centre)) < radius * 0.999)
    outside = sum(1 for n in lip
                  if _norm(_sub(nodes[n], centre)) > radius * 1.001)
    if inside > outside:
        return True
    crack["front_shape_note"] = "cirkel_buiten"
    return False


def csv_path(work_dir=None):
    """Pad van de tabel van deze run (<model>_<run_id>_scheur.csv)."""
    return os.path.join(work_dir or work_dir_path(),
                        run_name(CSV_SUFFIX))


# ==========================================
# NA DE RUN: TABEL LEZEN EN MELDEN
# ==========================================
def read_csv(path):
    """Rijen van de tabel als dicts (getallen als float; leeg = None), of
    [] als het bestand ontbreekt."""
    if not path or not os.path.isfile(path):
        return []
    with io.open(path, encoding="utf-8") as f:
        lines = [l.strip() for l in f if l.strip()]
    if not lines:
        return []
    head = lines[0].split(",")
    rows = []
    for line in lines[1:]:
        vals = line.split(",")
        row = {}
        for h, v in zip(head, vals):
            if v == "":
                row[h] = None
                continue
            try:
                row[h] = float(v)
            except ValueError:
                row[h] = v
        rows.append(row)
    return rows


def _range(rows, key, scale=1.0):
    vals = [r[key] / scale for r in rows if isinstance(r.get(key), float)]
    if not vals:
        return None
    return min(vals), max(vals)


def _fmt_range(rg, digits=4):
    if rg is None:
        return "-"
    fmt = "%." + str(digits) + "g"
    # Gelijk binnen de getoonde cijfers: een waarde
    if fmt % rg[0] == fmt % rg[1]:
        return fmt % rg[0]
    return (fmt % rg[0]) + " .. " + (fmt % rg[1])


def report(work_dir=None, mess_path=None):
    """
    Samenvatting van de tabel van deze run voor de log en het statusvak:
    lijst teksten (leeg zonder scheur). WARNING's bij negatieve K1,
    ontbrekende POST_K1_K2_K3 en een groot verschil tussen K uit CALC_G en
    uit POST_K1_K2_K3 of tussen G en G_IRWIN. Waarden: Pa.sqrt(m) en J/m2
    (SI); in de tekst MPa.sqrt(m).
    """
    path = csv_path(work_dir)
    rows = read_csv(path)
    if not rows:
        return []
    out = []
    postk_errors = {}
    if mess_path and os.path.isfile(mess_path):
        with io.open(mess_path, encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("BRIDGE_SCHEUR_POSTK_FOUT"):
                    parts = line.split(None, 2)
                    postk_errors[parts[1] if len(parts) > 1 else "?"] = (
                        parts[2].strip() if len(parts) > 2 else "")
    keys = []
    for r in rows:
        k = r.get("crack")
        k = "%g" % k if isinstance(k, float) else (k or "")
        r["crack"] = k
        if k not in keys:
            keys.append(k)
    for n, key in enumerate(keys):
        sub = [r for r in rows if r["crack"] == key]
        k1 = _range(sub, "K1", 1.0e6)
        k1p = _range(sub, "K1_POST_K", 1.0e6)
        g = _range(sub, "G")
        line = tr(u"Breukmechanica, scheur {0}: {1} punten langs het front. "
                  u"K_I {2} MPa.m^0,5 (CALC_G), {3} MPa.m^0,5 "
                  u"(POST_K1_K2_K3); G {4} J/m2. K_II tot {5}, K_III tot {6} "
                  u"MPa.m^0,5. Tabel: {7}",
                  u"Fracture mechanics, crack {0}: {1} points along the "
                  u"front. K_I {2} MPa.m^0.5 (CALC_G), {3} MPa.m^0.5 "
                  u"(POST_K1_K2_K3); G {4} J/m2. K_II up to {5}, K_III up "
                  u"to {6} MPa.m^0.5. Table: {7}").format(
                      key or (n + 1), len(sub), _fmt_range(k1),
                      _fmt_range(k1p), _fmt_range(g),
                      _fmt_abs(sub, "K2"), _fmt_abs(sub, "K3"),
                      os.path.basename(path))
        log(u"INFO: " + line)
        out.append(line)
        for text in _warnings(sub, key or str(n + 1), postk_errors):
            log(text)
            out.append(text[len("WARNING: "):])
    return out


def _fmt_abs(rows, key):
    vals = [abs(r[key]) / 1.0e6 for r in rows if isinstance(r.get(key), float)]
    return "%.3g" % max(vals) if vals else "-"


# Grenzen voor de meldingen (keuze van de bridge; verificatie
# docs/Verificatie.md BR: op een goed mesh < 2 %)
_DIFF_WARN = 0.05


def _warnings(rows, key, postk_errors):
    out = []
    k1 = [r["K1"] for r in rows if isinstance(r.get("K1"), float)]
    if k1 and min(k1) < 0.0:
        out.append(tr(u"WARNING: Breukmechanica, scheur {0}: K_I is negatief "
                      u"op {1} van {2} punten: de lippen gaan daar over "
                      u"elkaar (druk op de scheur; de bridge rekent zonder "
                      u"contact op de lippen, de waarde is dan fysisch "
                      u"zinloos). Controleer de belasting; boven- en "
                      u"onderlip verwisselen verandert het teken niet.",
                      u"WARNING: Fracture mechanics, crack {0}: K_I is "
                      u"negative at {1} of {2} points: the lips overlap "
                      u"there (compression on the crack; the bridge has no "
                      u"contact on the lips, the value is then physically "
                      u"meaningless). Check the load; swapping upper and "
                      u"lower lip does not change the sign.").format(
                          key, sum(1 for v in k1 if v < 0.0), len(k1)))
    if key in postk_errors or all(r.get("K1_POST_K") is None for r in rows):
        out.append(tr(u"WARNING: Breukmechanica, scheur {0}: POST_K1_K2_K3 "
                      u"gaf geen waarden ({1}); alleen CALC_G.",
                      u"WARNING: Fracture mechanics, crack {0}: "
                      u"POST_K1_K2_K3 gave no values ({1}); CALC_G only.")
                   .format(key, postk_errors.get(key, "-")))
    pairs = [(r["K1"], r["K1_POST_K"]) for r in rows
             if isinstance(r.get("K1"), float) and
             isinstance(r.get("K1_POST_K"), float)]
    ref = max([abs(a) for a, _b in pairs] + [0.0])
    if pairs and ref > 0.0:
        # Alleen binnenpunten: aan een vrij oppervlak verschillen de
        # methoden altijd (randeffect, docs/Breukmechanica.md)
        inner = pairs[1:-1] if len(pairs) > 2 else pairs
        diff = max(abs(a - b) for a, b in inner) / ref
        if diff > _DIFF_WARN:
            out.append(tr(u"WARNING: Breukmechanica, scheur {0}: K_I uit "
                          u"CALC_G en uit POST_K1_K2_K3 verschillen tot {1}% "
                          u"(binnen het front). Verfijn het mesh rond het "
                          u"front of controleer R_INF/R_SUP; bij een krom "
                          u"front (geen rechte lijn of cirkel) is K uit "
                          u"CALC_G zonder krommingscorrectie.",
                          u"WARNING: Fracture mechanics, crack {0}: K_I from "
                          u"CALC_G and from POST_K1_K2_K3 differ by up to "
                          u"{1}% (inside the front). Refine the mesh around "
                          u"the front or check R_INF/R_SUP; for a curved "
                          u"front (not a straight line or circle) K from "
                          u"CALC_G has no curvature correction.").format(
                              key, "%.1f" % (100.0 * diff)))
    return out
