# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# bridge_adapt.py
# Verantwoordelijkheid: adaptieve verfijning (stap S7) van een lineair
# statisch model (lin_static_3d) met Mecway's eigen verfijning: na een run
# met de foutschatter (CALC_ERREUR ERME, bridge_to_aster.
# _error_estimator_block) de elementen met de grootste geschatte fout
# kiezen, ze in Mecway verfijnen (mw.refine_x2) en opnieuw rekenen, tot een
# doel (globale relatieve fout) of een grens (aantal lussen, elementen).
# Gebruikersdoc: docs/Foutschatter.md par. 6.
#
#   adapt_settings(analysis, options) - instellingen (+ #ASTER-OPTIE), of
#                                       None als de verfijning uit staat
#   read_element_errors(fragment, ..) - eta_e per element uit het
#                                       resultaatbestand (ERME_ELNO)
#   select_elements(errors, ..)       - te verfijnen elementen (Dorfler,
#                                       fractie of factor x gemiddelde)
#   node_bc_conflicts(model_data, ..) - knoopselecties met een voorwaarde
#                                       op een rand die verfijnd wordt
#   after_run(history, ...)           - beslissing na een run (lus van het
#                                       hoofdscript, MecwayAsterBridgeMain)
#   refine_in_mecway(mw, eids)        - mw.refine_x2 + controle
#
# Wat Mecway's refine_x2 doet (headless gemeten, mecway.exe <model> script,
# Mecway 33; docs/Foutschatter.md par. 6.4):
# - elk gekozen element wordt 2 x 2 x 2 elementen van dezelfde vorm; de
#   buren krijgen overgangselementen (hex8 -> hex8 + pyr5 + tet4, hex20 ->
#   hex20 + pyr13 + tet10), het mesh blijft conform (geen hangende knopen;
#   manual par. 4.4.1 "correctly connected"). Middenknopen van pyr13/tet10
#   liggen op het midden van hun ribbe, in de volgorde van bridge_mesh.
#   _MID_EDGES (nagelopen op alle elementen van een verfijnd hex20-blok);
# - vlakselecties (en dus vlaklasten en -voorwaarden), elementselecties en
#   onderdelen (elset, materiaal) gaan mee naar de nieuwe elementen;
#   knoopselecties NIET: de nieuwe knopen op een rand of vlak van een
#   knoopselectie horen er niet bij (ook manual par. 9: "loads and
#   constraints applied to nodes remain at the same point locations");
# - nieuwe knopen liggen op de (kwadratische) elementgeometrie: bij hex20
#   blijft een cirkelrand een cirkel (afwijking 5e-4 r na een verfijning
#   van een 45-gradenboog), bij hex8 op de koorde (de rand blijft de
#   veelhoek van het grove mesh);
# - elementen en knopen worden hernummerd (de verfijnde en de
#   overgangselementen krijgen nieuwe nummers);
# - na een paar lussen kunnen de overgangen slechte elementen geven: bij
#   dunne elementen (dikte 1/40 van de breedte) na 2 lussen bijna platte
#   hex20 (geschaalde Jacobiaan 0,015) en een stijgende foutschatting; bij
#   tet10 al na de eerste lus (overgangen met wedge15/pyr13); bij hex20 met
#   een verhouding ~1 (verificatie AD-02) daalde de fout in 5 lussen
#   monotoon (docs/Verificatie.md, AD).
#
# Compatibel met IronPython 2.7.9 en CPython 3 (verificatie; alleen de
# standaardbibliotheek).

import base64
import math
import re
import struct

from bridge_utils import log, tr, parse_yes_no

# Standaardinstellingen (ook in CONFIG_BRIDGE["analyses"]["lin_static_3d"];
# hier als terugval voor een oudere bridge_utils.py).
DEFAULTS = {
    "adapt_refine":       False,
    "adapt_target_error": 2.0,
    "adapt_max_loops":    5,
    "adapt_select":       "dorfler",
    "adapt_fraction":     0.5,
    "adapt_neighbours":   False,
    "adapt_max_elements": 200000,
}
SELECT_METHODS = ("dorfler", "fraction", "factor")

# Foutveld per element in het resultaatbestand (aster_results_export
# FIELD_TABLE: ERME_ELNO ERREST -> elementnodevariable "errorestimate").
ERROR_VARIABLE = "errorestimate"

# Groei van het aantal elementen per verfijnd element (8 kinderen plus de
# overgangselementen in de buren), voor de grens adapt_max_elements
# vooraf. Gemeten (headless Mecway 33): 2 van 8 hex8 -> 46 elementen (19
# per verfijnd element); plaat met gat hex20 38 tot 70 per element. Na de
# eerste lus telt de gemeten groei van de vorige lus (_growth).
_GROWTH = 40.0

# Hoekparen van de ribben per volumevorm (Mecway-knoopvolgorde = Code-
# Aster-volgorde voor volumes, bridge_mesh._ASTER_TYPES; ribben als
# bridge_mesh._MID_EDGES).
_EDGES = {
    "hex": [(0, 1), (1, 2), (2, 3), (0, 3), (0, 4), (1, 5), (2, 6), (3, 7),
            (4, 5), (5, 6), (6, 7), (4, 7)],
    "tet": [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)],
    "pyr": [(0, 1), (1, 2), (2, 3), (0, 3), (0, 4), (1, 4), (2, 4), (3, 4)],
    "wedge": [(0, 1), (1, 2), (0, 2), (0, 3), (1, 4), (2, 5), (3, 4), (4, 5),
              (3, 5)],
}


# ==========================================
# INSTELLINGEN
# ==========================================
def _option_value(old, raw):
    """#ASTER-OPTIE-tekst naar het type van de standaardwaarde (zoals
    bridge_custom.apply_options). Geen geldig ja/nee-woord of getal:
    ValueError (de standaardwaarde blijft; apply_options geeft de
    WARNING)."""
    if isinstance(old, bool):
        value = parse_yes_no(raw)
        if value is None:
            raise ValueError(raw)
        return value
    if isinstance(old, int):
        return int(raw)
    if isinstance(old, float):
        return float(raw)
    return raw.strip("'\"")


def adapt_settings(analysis, options=None):
    """
    Instellingen van de adaptieve verfijning voor de analyse (dict uit
    CONFIG_BRIDGE["analyses"]) met de #ASTER-OPTIE-waarden van het model
    (model_data["aster_custom"]["options"]; die winnen, ongeldige waarden
    meldt bridge_custom.apply_options al). Returns de dict met de adapt_*-
    sleutels, of None als de verfijning uit staat of niet kan (geen
    MECA_STATIQUE of geen foutschatter).
    """
    options = options or {}
    out = {}
    for key, default in DEFAULTS.items():
        value = analysis.get(key, default)
        raw = options.get(key)
        if raw is not None:
            try:
                value = _option_value(default, "{0}".format(raw))
            except ValueError:
                pass
        out[key] = value
    if not out["adapt_refine"]:
        return None
    est = analysis.get("error_estimator", False)
    raw = parse_yes_no(options.get("error_estimator"))
    if raw is not None:
        est = raw
    if analysis.get("solver_command") != "MECA_STATIQUE" or not est:
        log(tr("WARNING: Adaptieve verfijning (adapt_refine) vraagt een "
               "lineair statische analyse met de foutschatter "
               "(error_estimator aan); niet verfijnd.",
               "WARNING: Adaptive refinement (adapt_refine) needs a linear "
               "static analysis with the error estimator (error_estimator "
               "on); not refined."))
        return None
    if out["adapt_select"] not in SELECT_METHODS:
        log(tr("WARNING: adapt_select '{0}' onbekend (bekend: {1}); "
               "dorfler gebruikt.",
               "WARNING: adapt_select '{0}' unknown (known: {1}); dorfler "
               "used.").format(out["adapt_select"],
                               ", ".join(SELECT_METHODS)))
        out["adapt_select"] = "dorfler"
    return out


# ==========================================
# FOUT PER ELEMENT UIT HET RESULTAATBESTAND
# ==========================================
_TAG = re.compile(r'<(results|elementnodeids|elementnodevariable|'
                  r'nodevariable|valuessingle|valuesdouble)\b([^>]*)>')
_NAME = re.compile(r'name="([^"]*)"')


def _text(line):
    if not isinstance(line, str):           # CPython 3: bytes
        line = line.decode("utf-8", "replace")
    return line


def _unpack(chunks, code):
    raw = base64.b64decode("".join(chunks))
    size = struct.calcsize(code)
    return struct.unpack("<%d%s" % (len(raw) // size, code), raw)


def read_element_errors(fragment_path, elements, variable=ERROR_VARIABLE):
    """
    Geschatte fout eta_e per element (ERME_ELNO, component ERREST; per
    element constant) uit een resultaatbestand van de bridge
    (aster_results_export, ook een .liml met <results>).

    elements : [(eid, vorm, [nid, ...], elset)] van het model
               (model_data["mesh_elements"]): het aantal waarden per element
               is zijn aantal knopen (elementnodevariable, knopen in de
               volgorde van <elem nodes>).
    Returns {eid: eta_e}; {} als het veld ontbreekt (foutschatter uit, of
    route 'tabel').
    """
    n_nodes = dict((e[0], len(e[2])) for e in elements)
    ids_chunks, val_chunks = None, None
    code = "f"
    target = None
    in_var = False
    f = open(fragment_path, "rb")
    try:
        for line in f:
            line = _text(line)
            if target is not None:
                part = line.split("<", 1)[0].strip()
                if part:
                    target.append(part)
                if "<" in line:
                    target = None
                continue
            m = _TAG.search(line)
            if m is None:
                continue
            kind, attrs = m.group(1), m.group(2)
            name = _NAME.search(attrs)
            name = name.group(1) if name else None
            sink = None
            if kind == "results":
                code = "d" if 'precision="double"' in attrs else "f"
            elif kind == "elementnodeids" and name == variable:
                ids_chunks = sink = []
            elif kind == "elementnodevariable":
                in_var = name == variable and val_chunks is None
            elif kind == "nodevariable":
                in_var = False
            elif kind in ("valuessingle", "valuesdouble") and in_var:
                code = "d" if kind == "valuesdouble" else "f"
                val_chunks = sink = []
                in_var = False
            if sink is not None:
                rest = line[m.end():]
                part = rest.split("<", 1)[0].strip()
                if part:
                    sink.append(part)
                if "<" not in rest:
                    target = sink
    finally:
        f.close()
    if not ids_chunks or not val_chunks:
        return {}
    eids = _unpack(ids_chunks, "i")
    vals = _unpack(val_chunks, code)
    out = {}
    pos = 0
    for eid in eids:
        n = n_nodes.get(eid)
        if n is None or pos + n > len(vals):
            log(tr("WARNING: Adaptieve verfijning: element {0} van het "
                   "resultaat past niet bij het model; foutwaarden niet "
                   "gelezen.",
                   "WARNING: Adaptive refinement: element {0} of the result "
                   "does not fit the model; error values not read.").format(
                       eid))
            return {}
        # float(): IronPython geeft uit struct een System.Single/Double.
        out[eid] = float(vals[pos])
        pos += n
    return out


# ==========================================
# KEUZE VAN DE ELEMENTEN
# ==========================================
def select_elements(errors, method="dorfler", value=0.5):
    """
    Te verfijnen elementen uit {eid: eta_e}:
      dorfler  - de kleinste groep elementen met de grootste fout die samen
                 een fractie value (0 < value <= 1) van de totale fout in
                 het kwadraat (sum eta_e^2) dragen ("bulk criterion",
                 Dorfler 1996); standaard 0,5;
      fraction - de fractie value van de elementen met de grootste fout
                 (0,1 = 10%);
      factor   - de elementen met eta_e > value x het gemiddelde van eta_e.
    Gelijke waarden: lager elementnummer eerst. Returns gesorteerde lijst.
    """
    items = sorted(errors.items(), key=lambda kv: (-kv[1], kv[0]))
    if not items:
        return []
    if method == "fraction":
        n = int(math.ceil(max(0.0, min(1.0, value)) * len(items) - 1e-9))
        chosen = items[:max(1, n)]
    elif method == "factor":
        mean = sum(v for _e, v in items) / float(len(items))
        chosen = [kv for kv in items if kv[1] > value * mean]
    else:
        total = sum(v * v for _e, v in items)
        goal = max(0.0, min(1.0, value)) * total
        chosen, acc = [], 0.0
        for kv in items:
            if acc >= goal and chosen:
                break
            chosen.append(kv)
            acc += kv[1] * kv[1]
    return sorted(e for e, _v in chosen)


def _shape_family(shape):
    for fam in ("hex", "tet", "pyr", "wedge"):
        if shape.startswith(fam):
            return fam
    return None


def add_neighbours(elements, eids):
    """
    eids plus alle volume-elementen die een knoop met een van eids delen
    (een bufferlaag; instelling adapt_neighbours). Gemeten (headless
    Mecway 33, dunne plaat met gat, hex20 van 1/40 van de breedte dik):
    refine_x2 op alleen de gekozen elementen gaf in de tweede lus, waar de
    keuze aan de overgangselementen van de eerste lus grensde, bijna
    platte hex20 (geschaalde Jacobiaan 0,015; globale schatting van 3,7%
    naar 19,6%); met de bufferlaag bleef de slechtste hex tot en met lus 3
    0,68 en kwam het pas in lus 4 (0,009). Bij de dikke plaat (AD-02)
    kost de buffer alleen elementen: Kt binnen 2% bij 22956 i.p.v. 15810
    elementen. elements: model_data["mesh_elements"].
    """
    wanted = set(eids)
    by_node = {}
    solid = []
    for eid, shape, nids, _elset in elements:
        if _shape_family(shape) is None:
            continue
        solid.append((eid, nids))
        for n in nids:
            by_node.setdefault(n, []).append(eid)
    out = set(wanted)
    for eid, nids in solid:
        if eid in wanted:
            for n in nids:
                out.update(by_node[n])
    return sorted(out)


# Hoekknopen van elk hoekpunt en zijn drie ribben (Mecway-volgorde):
# hex via de afgeleiden van de trilineaire afbeelding in de hoeken; tet,
# wedge, pyr (basis) met de drie ribben uit dat hoekpunt.
_HEX_RST = [(-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
            (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)]
_CORNER_EDGES = {
    "tet": [(0, 1, 2, 3), (1, 2, 0, 3), (2, 0, 1, 3), (3, 0, 1, 2)],
    "wedge": [(0, 1, 2, 3), (1, 2, 0, 4), (2, 0, 1, 5), (3, 5, 4, 0),
              (4, 3, 5, 1), (5, 4, 3, 2)],
    "pyr": [(0, 1, 3, 4), (1, 2, 0, 4), (2, 3, 1, 4), (3, 0, 2, 4)],
}


def _det(a, b, c):
    return (a[0] * (b[1] * c[2] - b[2] * c[1]) -
            a[1] * (b[0] * c[2] - b[2] * c[0]) +
            a[2] * (b[0] * c[1] - b[1] * c[0]))


def _scaled(a, b, c):
    la = math.sqrt(a[0] ** 2 + a[1] ** 2 + a[2] ** 2)
    lb = math.sqrt(b[0] ** 2 + b[1] ** 2 + b[2] ** 2)
    lc = math.sqrt(c[0] ** 2 + c[1] ** 2 + c[2] ** 2)
    if la * lb * lc <= 0.0:
        return 0.0
    return abs(_det(a, b, c)) / (la * lb * lc)


def scaled_jacobian(shape, pts):
    """
    Kleinste geschaalde Jacobiaan in de hoeken (|det| van de drie ribben /
    hun lengtes; 1 = rechte hoeken, 0 = plat), alleen de hoekknopen
    (middenknopen niet). pts: knoopcoordinaten in Mecway-volgorde.
    """
    fam = _shape_family(shape)
    if fam == "hex":
        worst = 1.0
        for r, s, t in _HEX_RST:
            d = [[0.0, 0.0, 0.0] for _k in range(3)]
            for (ri, si, ti), p in zip(_HEX_RST, pts[:8]):
                w = (ri * (1 + si * s) * (1 + ti * t),
                     si * (1 + ri * r) * (1 + ti * t),
                     ti * (1 + ri * r) * (1 + si * s))
                for k in range(3):
                    for j in range(3):
                        d[k][j] += w[k] * p[j]
            worst = min(worst, _scaled(*d))
        return worst
    worst = 1.0
    for i, a, b, c in _CORNER_EDGES.get(fam, []):
        o = pts[i]
        worst = min(worst, _scaled(*[[pts[x][j] - o[j] for j in range(3)]
                                     for x in (a, b, c)]))
    return worst


# Grens voor een "slecht" element (geschaalde Jacobiaan van een hex):
# keuze van de bridge; het grove plaatmodel heeft 0,71, de slechte
# overgangen 0,002 tot 0,02. Alleen hex: bij tet/pyr/wedge hangt de maat
# ook van de slankheid af (een gave tet in een dunne laag geeft al 0,02).
BAD_ELEMENT = 0.1


def bad_elements(model_data):
    """Hex-elementen met een geschaalde Jacobiaan < BAD_ELEMENT:
    [(eid, waarde)], slechtste eerst."""
    nodes = model_data.get("mesh_nodes", {})
    out = []
    for eid, shape, nids, _elset in model_data.get("mesh_elements", []):
        if _shape_family(shape) != "hex":
            continue
        try:
            q = scaled_jacobian(shape, [nodes[n] for n in nids])
        except KeyError:
            continue
        if q < BAD_ELEMENT:
            out.append((eid, q))
    return sorted(out, key=lambda x: x[1])


def node_bc_conflicts(model_data, eids):
    """
    Knoopselecties met een belasting of voorwaarde (aster_target GROUP_NO)
    waarvan twee hoekknopen een ribbe van een te verfijnen element vormen:
    Mecway zet de nieuwe knoop op die ribbe niet in de selectie, dus na het
    verfijnen zou de voorwaarde daar ontbreken. Returns [(selectie, soort)].
    """
    groups = model_data.get("node_groups", {})
    node_bcs = []
    for bc in model_data.get("bcs", []):
        if bc.get("aster_target") != "GROUP_NO":
            continue
        # node_groups: Code-Aster-namen "N<nid>" (bridge_from_liml STAP 1)
        nodes = set()
        for n in groups.get(bc.get("selection"), []):
            try:
                nodes.add(int("{0}".format(n).lstrip("N")))
            except ValueError:
                pass
        if len(nodes) > 1:
            node_bcs.append((bc.get("selection"), bc.get("type"), nodes))
    if not node_bcs:
        return []
    wanted = set(eids)
    edges = set()
    for eid, shape, nids, _elset in model_data.get("mesh_elements", []):
        if eid not in wanted:
            continue
        fam = _shape_family(shape)
        for a, b in _EDGES.get(fam, []):
            if a < len(nids) and b < len(nids):
                edges.add((nids[a], nids[b]))
    out = []
    for sel, kind, nodes in node_bcs:
        if any(a in nodes and b in nodes for a, b in edges):
            if (sel, kind) not in out:
                out.append((sel, kind))
    return out


# ==========================================
# BESLISSING NA EEN RUN
# ==========================================
def _fmt(v, digits=3):
    # %-opmaak: "{0:.3g}".format() volgt in IronPython de Windows-
    # landinstelling (MecwayAsterBridgeMain._report_modes).
    return ("%." + str(digits) + "g") % v


def summary_lines(history):
    """Tabel van de lussen voor log en statusvak."""
    lines = [tr(u"Adaptieve verfijning: lus, elementen (volume), globale "
                u"relatieve fout, verfijnd:",
                u"Adaptive refinement: loop, elements (solid), global "
                u"relative error, refined:")]
    for h in history:
        lines.append(u"  {0}: {1}, {2}%, {3}".format(
            h["loop"], h["n_elem"], _fmt(h["rel"]),
            h.get("n_refined", "-")))
    return lines


def _growth(history):
    """Nieuwe elementen per verfijnd element in de vorige lus (history[-2]
    -> history[-1]), anders _GROWTH. Telt alleen volume-elementen."""
    if len(history) >= 2 and history[-2].get("n_refined"):
        grow = (history[-1]["n_elem"] - history[-2]["n_elem"]) / float(
            history[-2]["n_refined"])
        if grow > 0:
            return grow
    return _GROWTH


def after_run(history, settings, model_data, fragment_path, error_summary):
    """
    Beslissing na een run van de lus (MecwayAsterBridgeMain._adapt_step).

    history       : lijst van eerdere lussen (wordt hier aangevuld met een
                    dict {"loop", "n_elem", "rel", "abs", "n_refined"})
    settings      : adapt_settings()
    model_data    : van deze run (mesh_elements, bcs, node_groups)
    fragment_path : resultaatbestand van deze run
    error_summary : laatste dict van bridge_results_liml.read_error_summary
                    (globale fout uit de .mess) of None

    Returns {"refine": bool, "elements": [eid], "notes": [tekst voor het
    statusvak]}; alle meldingen staan ook in de log.
    """
    notes = []

    def stop(text, level="INFO"):
        log(u"{0}: {1}".format(level, text))
        notes.append(text)
        notes.extend(summary_lines(history))
        for line in summary_lines(history):
            log(u"INFO: " + line)
        return {"refine": False, "elements": [], "notes": notes}

    loop = len(history)
    if not error_summary:
        history.append({"loop": loop, "n_elem": 0, "rel": float("nan"),
                        "abs": float("nan")})
        return stop(tr(u"Adaptieve verfijning gestopt: geen globale "
                       u"foutschatting in de .mess (foutschatter niet "
                       u"gedraaid of afgebroken).",
                       u"Adaptive refinement stopped: no global error "
                       u"estimate in the .mess (error estimator did not run "
                       u"or was aborted)."), "WARNING")
    errors = read_element_errors(fragment_path,
                                 model_data.get("mesh_elements", []))
    rec = {"loop": loop, "n_elem": len(errors), "rel": error_summary["rel"],
           "abs": error_summary.get("abs")}
    history.append(rec)
    log(u"INFO: Adaptieve verfijning lus {0}: {1} volume-elementen, "
        u"globale relatieve fout {2}%.".format(loop, len(errors),
                                               _fmt(rec["rel"])))
    bad = bad_elements(model_data)
    if bad:
        text = tr(u"WARNING: Adaptieve verfijning: {0} element(en) met een "
                  u"slechte vorm (geschaalde Jacobiaan < {1}; slechtste: "
                  u"element {2}, {3}). Daar zijn spanning en foutschatting "
                  u"onbetrouwbaar; bekijk het mesh in Mecway.",
                  u"WARNING: Adaptive refinement: {0} element(s) with a poor "
                  u"shape (scaled Jacobian < {1}; worst: element {2}, {3}). "
                  u"Stress and error estimate are unreliable there; check "
                  u"the mesh in Mecway.").format(
                      len(bad), BAD_ELEMENT, bad[0][0], _fmt(bad[0][1], 2))
        log(text)
        notes.append(text[len("WARNING: "):])
    target = float(settings["adapt_target_error"])
    if rec["rel"] <= target:
        return stop(tr(u"Adaptieve verfijning: doel bereikt (globale "
                       u"relatieve fout {0}% <= {1}%).",
                       u"Adaptive refinement: target reached (global "
                       u"relative error {0}% <= {1}%).").format(
                           _fmt(rec["rel"]), _fmt(target)))
    if loop >= 1 and rec["rel"] >= history[-2]["rel"]:
        return stop(tr(u"WARNING: Adaptieve verfijning gestopt: de globale "
                       u"fout daalt niet meer ({0}% -> {1}%). Mogelijk een "
                       u"singulariteit (scherpe hoek, puntlast, rand van een "
                       u"inklemming) of een last die de schatter niet kent "
                       u"(docs/Foutschatter.md par. 2.3).",
                       u"WARNING: Adaptive refinement stopped: the global "
                       u"error no longer decreases ({0}% -> {1}%). Possibly "
                       u"a singularity (sharp corner, point load, edge of a "
                       u"clamp) or a load the estimator does not know "
                       u"(docs/en/Error_estimator.md sec. 2.3).").format(
                           _fmt(history[-2]["rel"]), _fmt(rec["rel"]))
                    [len("WARNING: "):], "WARNING")
    if loop >= int(settings["adapt_max_loops"]):
        return stop(tr(u"Adaptieve verfijning: maximaal aantal lussen ({0}) "
                       u"bereikt; globale relatieve fout {1}% (doel {2}%).",
                       u"Adaptive refinement: maximum number of loops ({0}) "
                       u"reached; global relative error {1}% (target "
                       u"{2}%).").format(settings["adapt_max_loops"],
                                         _fmt(rec["rel"]), _fmt(target)))
    if not errors:
        return stop(tr(u"Adaptieve verfijning gestopt: geen fout per element "
                       u"in het resultaatbestand (veld errorestimate; "
                       u"resultaatroute 'liml' nodig).",
                       u"Adaptive refinement stopped: no error per element "
                       u"in the result file (field errorestimate; result "
                       u"route 'liml' needed)."), "WARNING")
    chosen = select_elements(errors, settings["adapt_select"],
                             float(settings["adapt_fraction"]))
    # Optioneel een bufferlaag (add_neighbours; dunne elementen).
    eids = chosen
    if settings.get("adapt_neighbours"):
        eids = add_neighbours(model_data.get("mesh_elements", []), chosen)
    n_all = len(model_data.get("mesh_elements", []))
    estimate = n_all + _growth(history) * len(eids)
    if estimate > int(settings["adapt_max_elements"]):
        return stop(tr(u"Adaptieve verfijning gestopt: na verfijnen van {0} "
                       u"elementen zou het model ongeveer {1} elementen "
                       u"hebben (grens adapt_max_elements {2}).",
                       u"Adaptive refinement stopped: after refining {0} "
                       u"elements the model would have about {1} elements "
                       u"(limit adapt_max_elements {2}).").format(
                           len(eids), int(estimate),
                           settings["adapt_max_elements"]))
    conflicts = node_bc_conflicts(model_data, eids)
    if conflicts:
        return stop(tr(u"ERROR: Adaptieve verfijning gestopt: de "
                       u"knoopselectie(s) {0} met een belasting of "
                       u"voorwaarde liggen op ribben die verfijnd worden. "
                       u"Mecway zet de nieuwe knopen daar niet in de "
                       u"selectie, dus de voorwaarde zou na het verfijnen "
                       u"gaten krijgen. Leg de voorwaarde op vlakken "
                       u"(vlakselectie; die verfijnt Mecway mee) en start "
                       u"opnieuw.",
                       u"ERROR: Adaptive refinement stopped: the node "
                       u"selection(s) {0} with a load or constraint lie on "
                       u"edges that would be refined. Mecway does not add "
                       u"the new nodes there to the selection, so the "
                       u"constraint would get gaps after refining. Apply "
                       u"it to faces (face selection; Mecway refines those "
                       u"along) and start again.").format(
                           u", ".join(u"'{0}' ({1})".format(s, k)
                                      for s, k in conflicts))
                    [len("ERROR: "):], "ERROR")
    rec["n_refined"] = len(eids)
    rec["n_chosen"] = len(chosen)
    text = tr(u"Adaptieve verfijning lus {0}: globale relatieve fout {1}% > "
              u"doel {2}%; {3} van {4} volume-elementen worden verfijnd "
              u"(mw.refine_x2, keuze {5} {6}{7}).",
              u"Adaptive refinement loop {0}: global relative error {1}% > "
              u"target {2}%; {3} of {4} solid elements are refined "
              u"(mw.refine_x2, selection {5} {6}{7}).").format(
                  loop, _fmt(rec["rel"]), _fmt(target), len(eids),
                  len(errors), settings["adapt_select"],
                  _fmt(float(settings["adapt_fraction"])),
                  tr(u"; {0} gekozen plus de buren",
                     u"; {0} chosen plus their neighbours").format(
                         len(chosen)) if len(eids) != len(chosen) else u"")
    log(u"INFO: " + text)
    notes.append(text)
    return {"refine": True, "elements": eids, "notes": notes}


# ==========================================
# VERFIJNEN IN MECWAY
# ==========================================
def refine_in_mecway(mw, eids):
    """
    mw.refine_x2 op de elementen eids (API: refine_x2(node_ids, face_ids,
    element_ids), precies een van de drie). Returns het aantal elementen
    erna, of None bij een fout (ERROR in de log).
    """
    before = len(mw.all_elements())
    try:
        mw.refine_x2(None, None, list(eids))
    except Exception as e:
        log(tr(u"ERROR: Mecway-verfijning (refine_x2) mislukt: {0}",
               u"ERROR: Mecway refinement (refine_x2) failed: {0}").format(e))
        return None
    after = len(mw.all_elements())
    log(u"INFO: Adaptieve verfijning: refine_x2 op {0} elementen: {1} -> "
        u"{2} elementen, {3} knopen.".format(len(eids), before, after,
                                              len(mw.all_nodes())))
    if after <= before:
        log(tr(u"ERROR: Mecway-verfijning (refine_x2) gaf geen nieuwe "
               u"elementen.",
               u"ERROR: Mecway refinement (refine_x2) gave no new "
               u"elements."))
        return None
    return after
