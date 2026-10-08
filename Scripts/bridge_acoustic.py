# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# bridge_acoustic.py
# Verantwoordelijkheid: Mecway Acoustic Resonance 3D (A30, holtemodes) naar
# Code-Aster PHENOMENE='ACOUSTIQUE' + CALC_MODES, en hetzelfde model als
# harmonisch drukveld (bron, absorberende randen, microfoons) met
# DYNA_VIBRA HARM (Mecway heeft daar geen analysetype voor).
#
#   note_tag(item, model_data)        - aangeroepen door bridge_from_liml.py
#                                       voor elke belasting/voorwaarde in een
#                                       akoestische run (melden of bewaren).
#   parse_acoustic(root, model_data)  - na de lus in bridge_from_liml.py:
#                                       materialen (geluidssnelheid,
#                                       dichtheid) en drukknopen.
#   build_acoustic_modes(...)         - builder; bridge_to_aster._build_modal
#                                       stuurt een akoestische analyse
#                                       hierheen (CALC_MODES is de sleutel in
#                                       SOLVER_BUILDERS_ASTER).
#   build_acoustic_harm(...)          - builder harmonisch drukveld
#                                       (acoustic_harm_3d, DYNA_VIBRA HARM;
#                                       via bridge_to_aster._build_dyna_modal);
#                                       zie het blok HARMONISCH DRUKVELD.
#   resolve_auto(), report_harmonic() - MecwayAsterBridgeMain: 'auto' met
#                                       een bron, microfoons in de log.
#
# Mecway-XML (Mecway 33, vastgesteld met headless 'mecway <model> solve',
# TestFiles/Tools/mw_solve.py; tagnamen ook in mecway.exe):
#   <analysis type="A30" maxmode="5" />            (A20 = 2D, niet vertaald)
#   <mat><mechanical type="Isotropic" speedofsound="340 m/s" />
#        <density type="Simple" density="1.2 kg/m³" /></mat>
#     (dichtheid optioneel: Mecway rekent zonder; alleen een impedantie
#     vraagt hem. Eén materiaal: met twee materialen geeft Mecway geen
#     oplossing, gezien met headless solve.)
#   <fluidpressure selection="knopen of vlakken"><value>0 Pa</value>
#     (manual v33 par. 9.34: drukknoop p = 0; een andere waarde weigert
#     Mecway, "Pressure must be zero.")
#   <impedance selection="vlakken"><specificimpedancer>408</...r>
#     <specificimpedancei>0</...i></impedance>    (manual par. 9.35; getal
#     zonder eenheid = Pa.s/m; Z = r + i i; met impedantie complexe modes:
#     velden pressure + pressurei, <mode .. frequency=".." decayrate="..">;
#     vertaald in _complex_modes_comm + aster_acoustic_complex.py, S2b.)
#   Oplossing: per mode <mode id="k" frequency="f Hz"> met de knoopvariabele
#   "pressure"; Mecway normeert p^T M p = 1 met M = int N^T N / c^2 dV
#   (nagerekend: constante mode van een buis met volume V: p = c / sqrt(V);
#   manual v33 par. 6.9.1: amplitudes zijn willekeurig).
#
# Vertaling (Code-Aster 17.4; syntax als testcases ahlv100a en adlv101f):
#   - AFFE_MODELE PHENOMENE='ACOUSTIQUE', MODELISATION='3D' op de volumes
#     (elementen ACOU_HEXA8/20, ACOU_PENTA6/15, ACOU_TETRA4/10 in elem.1;
#     geen piramides).
#   - DEFI_MATERIAU FLUIDE=_F(RHO, CELE_R) (catalogus defi_materiau.py;
#     RHO is verplicht). Zonder dichtheid RHO = 1. RIGI_ACOU en MASS_ACOU
#     bevatten rho niet (gemeten, NIET in de broncode nagelopen: RHO 1,2
#     gaf dezelfde frequenties en MASS_GENE als RHO 1, en een buis van twee
#     media met rho 1,2 / 12 gaf de frequenties van gelijke dichtheid,
#     248,98 Hz i.p.v. analytisch 306,47 Hz). Net als Mecway (K = int
#     grad N^T grad N, manual v33 par. 6.9.1) is het grensvlak tussen twee
#     media dus dp/dn continu: alleen juist bij gelijke dichtheid. Media
#     met verschillende dichtheid: ERROR (_materials).
#   - Drukknoop: AFFE_CHAR_ACOU PRES_IMPO PRES=0j (adlv101f).
#   - ASSEMBLAGE RIGI_ACOU / MASS_ACOU (complexe matrices) en
#     COMB_MATR_ASSE PARTIE='REEL' (ahlv100a), dan CALC_MODES (resultaat
#     MODE_ACOU, veld PRES, component PRES; gezien in een run).
#   - Harde wanden zonder drukknoop: de stijfheid is singulier (constante
#     druk, 0 Hz). PLUS_PETITE stopt dan met <EXCEPTION> ALGELINE3_81
#     (gezien); de bridge zoekt rond een kleine negatieve frequentie
#     (OPTION 'CENTRE', _free_shift), zoals vrij-vrij bij modaal.
#   - NORM_MODE werkt niet voor MODE_ACOU (MASS_GENE als norm gaf toch
#     max |p| = 1 en een MASS_GENE van 8,65e-8, gezien): aster_results_
#     export schaalt de vormen zelf met 1/sqrt(MASS_GENE), dus p^T M p = 1
#     zoals Mecway.
#
# Compatibel met IronPython 2.7.9 (alleen standaardbibliotheek).

import math
import re
import os
from collections import OrderedDict

from bridge_utils import (log, parse_and_convert, CONFIG_BRIDGE, tr,
                          run_name, _normalize_unit_string)


# Mecway-belastingen die alleen in Acoustic Resonance bestaan (manual v33
# par. 9.34, 9.35); in een mechanische of thermische analyse spelen ze
# geen rol (INFO in bridge_from_liml).
ACOUSTIC_TAGS = frozenset(["fluidpressure", "impedance"])

# Voorwaarden die Mecway ook in de akoestiek gebruikt (manual v33 par. 9.9:
# bonded contact koppelt de druk; testmodellen in mecway.exe: "A30
# node-surface coupling", "A30 Constraint Equations"), maar die de bridge
# daar niet vertaalt: melden en overslaan.
_NOT_IN_ACOUSTICS = frozenset(["contact2", "nodesurfacecoupling",
                               "constraintequation"])

# Mecway-vormen met een akoestisch Code-Aster-element (elem.1 van 17.4:
# ACOU_HEXA8/20/27, ACOU_PENTA6/15, ACOU_TETRA4/10; geen piramide).
_ACOUSTIC_SHAPES = frozenset(["hex8", "hex20", "wedge6", "wedge15",
                              "tet4", "tet10"])
_LINEAR_SHAPES = frozenset(["hex8", "wedge6", "tet4"])

# Harmonische analyse (S3a): analysenaam en de voorvoegsels van speciale
# Mecway-selecties (hoofdletters maken niet uit; NL en EN).
HARM_ANALYSIS = "acoustic_harm_3d"
_PREFIXES = (("source", ("BRON_", "SOURCE_")),
             ("absorbing", ("ABSORBEREND_", "ABSORBING_")),
             ("microphone", ("MICROFOON_", "MICROPHONE_")))

# #ASTER-BRON (source_value; ook bridge_vibroacoustic_harm): getal met
# eenheid, de hele tekst. Snelheidseenheden die parse_and_convert(..,
# "velocity") omrekent (leeg = m/s) en druk-eenheden die parse_and_convert(..,
# "pressure") omrekent (kleine letters; mPa/MPa beslist parse_and_convert
# op de hoofdletter). Hier geeft een druk een duidelijke ERROR; in
# vibroacoustic_harm_3d is een druk een invallende golf.
_NUMBER_UNIT_RE = re.compile(
    r"^[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?\s*([A-Za-z0-9/]*)\s*$")
_VELOCITY_UNITS = ("", "m/s", "mm/s", "mm/ms", "in/s", "ft/s")
_PRESSURE_UNITS = ("pa", "n/m2", "kpa", "kn/m2", "mpa", "n/mm2", "gpa",
                   "bar", "mbar", "psi")

# Elementen per golflengte bij de hoogste frequentie waaronder de bridge
# waarschuwt (lineair / kwadratisch). Gemeten (verificatie AH-01/01b/02,
# docs/Akoestiek.md par. 5.4), vlakke golf over 2 golflengten: hex8 10 per
# golflengte fase -1,6% (~ (kh)^2/24) en |p| +1,7%, 20 per golflengte
# -0,4%; hex20 10 per golflengte <= 0,02%, 5 per golflengte |p| -0,25% en
# fase -0,16%.
_MIN_PER_WAVELENGTH = {True: 10.0, False: 5.0}


def is_acoustic(analysis_cfg):
    """True als de analyse-config akoestisch is (sleutel "phenomene")."""
    return str(analysis_cfg.get("phenomene", "")).upper() == "ACOUSTIQUE"


def _r(x):
    """Getal voor de .comm (15 cijfers, decimale punt; bridge_mesh._real)."""
    import bridge_mesh
    return bridge_mesh._real(float(x))


def _skip(model_data, text):
    model_data["skipped_loads"].append(text)


# ==========================================
# PARSEN
# ==========================================
def note_tag(item, model_data):
    """
    Een belasting of voorwaarde in een akoestische run (bridge_from_liml
    slaat ze daar over). Akoestische tags leest parse_acoustic; wat Mecway
    in de akoestiek gebruikt maar de bridge niet vertaalt: ERROR en
    skipped_loads; mechanische en thermische belastingen: INFO (ze doen in
    Acoustic Resonance niet mee, manual v33 hoofdstuk 9: analysetypen per
    belasting).
    """
    tag = item.tag
    # contact2: selectionmaster/-slave, nodesurfacecoupling:
    # selectionsurface/-node (bridge_from_liml)
    sel = item.get("selection") or " / ".join(
        item.get(k) for k in ("selectionmaster", "selectionslave",
                              "selectionsurface", "selectionnode")
        if item.get(k)) or ""
    if tag in ACOUSTIC_TAGS or item.get("suppressed", "") == "true":
        return
    if tag in _NOT_IN_ACOUSTICS:
        _skip(model_data, tr("{0} op '{1}' (akoestiek)",
                             "{0} on '{1}' (acoustics)").format(tag, sel))
        log(tr("ERROR: '{0}' op '{1}' wordt in een akoestische analyse niet "
               "vertaald en is NIET meegenomen (gebruik een doorlopend mesh "
               "zonder contact of koppelingen).",
               "ERROR: '{0}' on '{1}' is not translated in an acoustic "
               "analysis and is NOT included (use a continuous mesh without "
               "contact or couplings).").format(tag, sel))
        return
    log("INFO: Belasting '{0}' op '{1}' telt niet mee in een akoestische "
        "analyse.".format(tag, sel))


def _selection_kind(model_data, name):
    if name in model_data.get("node_groups", {}):
        return "node"
    if name in model_data.get("face_selection_faces", {}):
        return "face"
    if name in model_data.get("element_selections", {}):
        return "element"
    return None


def _prefix_role(name):
    """'source'/'absorbing'/'microphone' als de Mecway-naam met een
    voorvoegsel uit _PREFIXES begint (hoofdletters maken niet uit), anders
    None."""
    upper = (name or "").upper()
    for role, prefixes in _PREFIXES:
        for p in prefixes:
            if upper.startswith(p):
                return role
    return None


def _named_groups(model_data):
    """{rol: [(Mecway-naam, soort)]} van de selecties met een speciaal
    voorvoegsel (op de oorspronkelijke Mecway-naam, STAPPENPLAN B.2)."""
    out = {"source": [], "absorbing": [], "microphone": []}
    seen = set()
    for kind, names in (("face", model_data.get("face_selection_faces", {})),
                        ("node", model_data.get("node_groups", {})),
                        ("element",
                         model_data.get("element_selections", {}))):
        for name in names:
            role = _prefix_role(name)
            if role and (name, kind) not in seen:
                seen.add((name, kind))
                out[role].append((name, kind))
    return out


def source_value(text):
    """
    Waarde van een #ASTER-BRON: 'getal [eenheid]' als hele tekst ->
    (soort, waarde SI, eenheid zoals getypt). soort 'velocity' (eenheid uit
    _VELOCITY_UNITS, zonder eenheid m/s; waarde in m/s), 'pressure'
    (_PRESSURE_UNITS; waarde in Pa, een keer omgerekend door
    parse_and_convert) of None (geen getal, of een andere eenheid; waarde
    None). parse_and_convert zou een onbekende eenheid ongewijzigd en tekst
    zonder getal als 0.0 teruggeven, daarom eerst deze controle.
    """
    text = _normalize_unit_string(u"{0}".format(text).strip())
    m = _NUMBER_UNIT_RE.match(text)
    if not m:
        return None, None, ""
    unit = m.group(1)
    if unit.lower() in _VELOCITY_UNITS:
        kind = "velocity"
    elif unit.lower() in _PRESSURE_UNITS:
        kind = "pressure"
    else:
        return None, None, unit
    return kind, float(parse_and_convert(text, kind)), unit


def _parse_aster_sources(model_data, ac):
    """#ASTER-BRON selectie = snelheid [eenheid] [, fase] (bridge_custom)
    -> ac["aster_sources"] = [(selectie, v [m/s], fase [graden])].
    Alleen een snelheid: zonder eenheid m/s, verder de snelheidseenheden
    van parse_and_convert (_VELOCITY_UNITS). Een druk-eenheid (in
    vibroacoustic_harm_3d een invallende golf) of een andere eenheid en
    tekst zonder getal geven een ERROR (parse_and_convert zou de waarde
    ongewijzigd als m/s of als 0.0 teruggeven)."""
    for body in model_data.get("aster_custom", {}).get("sources", []):
        # 'Sel = 2 mm/s' of 'Sel = 2 mm/s, 90' of 'Sel = 2 mm/s, fase = 90'
        name, sep, rest = body.partition("=")
        parts = [p.strip() for p in rest.split(",")]
        v = None
        phase = 0.0
        if sep and name.strip() and parts[0] and len(parts) <= 2:
            kind, value, unit = source_value(parts[0])
            if kind == "pressure":
                ac["bad"].append(tr(
                    "'#ASTER-BRON {0}': een druk ({1}) kan hier niet als "
                    "bron; #ASTER-BRON is in deze analyse een "
                    "normaalsnelheid [m/s of mm/s]. Een opgelegde druk: "
                    "Mecway's fluid pressure != 0 op de selectie",
                    "'#ASTER-BRON {0}': a pressure ({1}) cannot be a "
                    "source here; in this analysis #ASTER-BRON is a "
                    "normal velocity [m/s or mm/s]. An imposed pressure: "
                    "Mecway's fluid pressure != 0 on the "
                    "selection").format(body.strip(), unit))
                continue
            if kind == "velocity":
                try:
                    v = value
                    if len(parts) == 2:
                        phase = float(parts[1].split("=")[-1].strip())
                except (TypeError, ValueError):
                    v = None
        if v is None:
            ac["bad"].append(tr(
                "'#ASTER-BRON {0}': verwacht 'selectie = snelheid [m/s of "
                "mm/s] [, fase in graden]', bijv. '#ASTER-BRON Luidspreker "
                "= 2 mm/s, 90'",
                "'#ASTER-BRON {0}': expected 'selection = velocity [m/s or "
                "mm/s] [, phase in degrees]', e.g. '#ASTER-BRON Speaker = "
                "2 mm/s, 90'").format(body.strip()))
            continue
        ac["aster_sources"].append((name.strip(), v, phase))


def parse_acoustic(root, model_data):
    """
    Materialen (geluidssnelheid, dichtheid) en akoestische voorwaarden uit
    de .liml naar model_data["acoustic_data"]:
      "materials":  {mat: {"c": m/s of None, "rho": kg/m3 of None}}
      "pressure0":  [(selectie, "node"/"face")]   drukknopen p = 0
      "pressure_src": [(selectie, soort, p [Pa], tekst)]  fluid pressure
                    != 0: in de harmonische analyse een opgelegde druk, in
                    de holtemodes een fout
      "impedances": [selectie]
      "impedance_values": [(selectie, Re Z, Im Z)]  (Pa.s/m)
      "named":      {rol: [(selectie, soort)]}  BRON_/ABSORBEREND_/
                    MICROFOON_ (harmonisch, _named_groups)
      "aster_sources": [(selectie, v, fase)]   #ASTER-BRON
      "bad":        [tekst]                       fatale invoerfouten
    """
    ac = {"materials": OrderedDict(), "pressure0": [], "pressure_src": [],
          "impedances": [], "impedance_values": [],
          "named": _named_groups(model_data), "aster_sources": [],
          "bad": []}
    _parse_aster_sources(model_data, ac)
    for mat in root.findall("mat"):
        name = mat.get("name", "UnknownMat")
        mech = mat.find("mechanical")
        dens = mat.find("density")
        c = mech.get("speedofsound") if mech is not None else None
        rho = dens.get("density") if dens is not None else None
        ac["materials"][name] = {
            "c": parse_and_convert(c, "velocity") if c else None,
            "rho": parse_and_convert(rho, "density") if rho else None}
        log("INFO: Akoestisch materiaal '{0}': c = {1} m/s, rho = {2} "
            "kg/m3.".format(name, ac["materials"][name]["c"],
                            ac["materials"][name]["rho"]))
    for item in root:
        if item.tag not in ACOUSTIC_TAGS:
            continue
        sel = item.get("selection", "")
        if item.get("suppressed", "") == "true":
            log("INFO: {0} op '{1}' overgeslagen (suppressed).".format(
                item.tag, sel))
            continue
        kind = _selection_kind(model_data, sel)
        if item.tag == "impedance":
            ac["impedances"].append(sel)
            ac["impedance_values"].append(
                (sel, _number(item.findtext("specificimpedancer")),
                 _number(item.findtext("specificimpedancei"))))
            continue
        text = item.findtext("value")
        value = parse_and_convert(text, "pressure") if text else 0.0
        try:
            value = float(value)
        except (TypeError, ValueError):
            ac["bad"].append(tr(
                "fluid pressure op '{0}': waarde '{1}' is geen druk",
                "fluid pressure on '{0}': value '{1}' is not a "
                "pressure").format(sel, text))
            continue
        if kind not in ("node", "face"):
            ac["bad"].append(tr(
                "fluid pressure op '{0}': kies een knoop- of vlakselectie",
                "fluid pressure on '{0}': choose a node or face "
                "selection").format(sel))
        elif abs(value) > 0.0:
            # Holtemodes: Mecway weigert dit ook ("Pressure must be
            # zero."), build_acoustic_modes geeft een ERROR; harmonisch:
            # een opgelegde druk (bron).
            ac["pressure_src"].append((sel, kind, value, text))
        else:
            ac["pressure0"].append((sel, kind))
    model_data["acoustic_data"] = ac


def _number(text):
    """Getal uit een Mecway-veld zonder eenheid (impedantie), 0 als leeg;
    None als het geen getal is."""
    if text is None or not text.strip():
        return 0.0
    try:
        return float(parse_and_convert(text, "dimensionless"))
    except (TypeError, ValueError):
        return None


def has_source(model_data):
    """True als het model een harmonische bron heeft (vlakselectie BRON_/
    SOURCE_, #ASTER-BRON of fluid pressure != 0): active_analysis 'auto'
    kiest dan bij A30 de harmonische analyse (MecwayAsterBridgeMain,
    resolve_auto)."""
    ac = model_data.get("acoustic_data") or {}
    return bool(ac.get("named", {}).get("source") or
                ac.get("aster_sources") or ac.get("pressure_src"))


def resolve_auto(resolved, model_data):
    """Analyse voor 'auto': A30 met een bron -> acoustic_harm_3d (als die
    bestaat), anders ongewijzigd."""
    if (resolved == "acoustic_modes_3d" and has_source(model_data) and
            HARM_ANALYSIS in CONFIG_BRIDGE.get("analyses", {})):
        log("INFO: Akoestiek: het model heeft een bron (BRON_-selectie, "
            "#ASTER-BRON of fluid pressure != 0); 'auto' kiest de "
            "harmonische analyse '{0}'. Voor "
            "de holtemodes: kies 'acoustic_modes_3d' als actieve "
            "analyse.".format(HARM_ANALYSIS))
        return HARM_ANALYSIS
    return resolved


# ==========================================
# BUILDER
# ==========================================
def _fmt_groups(groups):
    return "({0},)".format(", ".join("'{0}'".format(g) for g in groups))


def _elsets(model_data):
    """[(Code-Aster-elset, Mecway-materiaal)] van de volumes, of None bij
    een fout (ERROR gelogd). Na write_mesh zijn de sleutels van "elsets"
    de Code-Aster-namen (bridge_mesh._apply_to_model_data)."""
    shapes_of = OrderedDict()
    elem_elset = model_data.get("_elem_elset", {})
    for eid, shape, _nodes, elset in model_data.get("mesh_elements", []):
        name = elem_elset.get(eid, elset)
        shapes_of.setdefault(name, set()).add(shape)
    out, bad = [], []
    for elset, shapes in shapes_of.items():
        wrong = sorted(shapes - _ACOUSTIC_SHAPES)
        if wrong:
            bad.append("'{0}' ({1})".format(elset, ", ".join(wrong)))
            continue
        mat = model_data.get("elsets", {}).get(elset)
        if not mat:
            bad.append(tr("'{0}' (geen materiaal)",
                          "'{0}' (no material)").format(elset))
            continue
        out.append((elset, mat))
    if bad:
        log(tr("ERROR: Akoestiek: alleen volume-elementen (hex, wedge, tet; "
               "geen piramides, schalen of balken) met een materiaal: {0}. "
               "Mesh in Mecway alleen het lucht- of vloeistofvolume.",
               "ERROR: Acoustics: only solid elements (hex, wedge, tet; no "
               "pyramids, shells or beams) with a material: {0}. In Mecway, "
               "mesh only the air or fluid volume.").format(", ".join(bad)))
        return None
    if not out:
        log(tr("ERROR: Akoestiek: geen volume-elementen gevonden.",
               "ERROR: Acoustics: no solid elements found."))
        return None
    return out


def _materials(elsets, ac):
    """{Mecway-materiaal: (c, rho)} of None (ERROR). Zonder dichtheid rho
    = 1. Meer media: alleen met gelijke (of geen) dichtheid, want de
    matrices bevatten rho niet (zie de kop)."""
    used = []
    for _es, mat in elsets:
        if mat not in used:
            used.append(mat)
    props = ac["materials"]
    no_c = [m for m in used if not (props.get(m, {}).get("c") or 0.0) > 0.0]
    if no_c:
        log(tr("ERROR: Akoestiek: materiaal {0} heeft geen geluidssnelheid. "
               "Vul in Mecway bij het materiaal (Isotropic) 'Speed of sound' "
               "in.",
               "ERROR: Acoustics: material {0} has no speed of sound. Fill "
               "in 'Speed of sound' for the material (Isotropic) in "
               "Mecway.").format(", ".join("'%s'" % m for m in no_c)))
        return None
    rhos = dict((m, props[m].get("rho") or 0.0) for m in used)
    given = sorted(set(v for v in rhos.values() if v > 0.0))
    if len(given) > 1 and given[-1] > given[0] * (1.0 + 1.0e-6):
        log(tr("ERROR: Akoestiek: de media hebben een verschillende "
               "dichtheid ({0}). Code-Aster's akoestische drukformulering "
               "(en Mecway) neemt een dichtheidssprong op het grensvlak niet "
               "mee; de frequenties zouden fout zijn. Geef alle materialen "
               "dezelfde dichtheid (of geen) als het verschil klein is, of "
               "reken de media apart.",
               "ERROR: Acoustics: the media have different densities ({0}). "
               "Code-Aster's acoustic pressure formulation (and Mecway) does "
               "not include a density jump at the interface; the "
               "frequencies would be wrong. Give all materials the same "
               "density (or none) if the difference is small, or compute "
               "the media separately.").format(", ".join(
                   "%s %s kg/m3" % (m, "%.6g" % rhos[m]) for m in used)))
        return None
    missing = [m for m in used if not rhos[m] > 0.0]
    if given and missing:
        # Eén medium met en een zonder dichtheid: de matrices bevatten rho
        # niet (zie de kop), dus het grensvlak (p en dp/dn continu) rekent
        # alsof de dichtheden gelijk zijn. Het medium zonder dichtheid
        # krijgt RHO = 1 in DEFI_MATERIAU; dat telt alleen voor een bron
        # of impedantie op dat medium, en die eisen zelf een dichtheid.
        # Doorgerekend (07-10-2026): AK-05 met 1,2 kg/m3 / geen geeft de
        # frequenties van AK-05 (gelijke dichtheid), RHO 1.2 en 1.0.
        log(tr("WARNING: Akoestiek: medium {0} heeft geen dichtheid, "
               "medium {1} wel ({2}). Op het grensvlak rekent de bridge "
               "alsof de dichtheden gelijk zijn (druk en dp/dn continu): "
               "dat klopt alleen als {0} dezelfde dichtheid heeft; bij een "
               "andere dichtheid zijn de resultaten fout (docs/Akoestiek.md "
               "par. 3.2). Code-Aster krijgt voor {0} RHO = 1 kg/m3. Vul bij "
               "elk medium de dichtheid in. De berekening gaat door.",
               "WARNING: Acoustics: medium {0} has no density, medium {1} "
               "has ({2}). At the interface the bridge computes as if the "
               "densities are equal (pressure and dp/dn continuous): that "
               "is only correct if {0} has the same density; with a "
               "different density the results are wrong (docs/en/"
               "Acoustics.md section 3.2). Code-Aster gets RHO = 1 kg/m3 "
               "for {0}. Fill in the density for every medium. The "
               "computation continues.").format(
                   ", ".join("'%s'" % m for m in missing),
                   ", ".join("'%s'" % m for m in used if m not in missing),
                   ", ".join("%s kg/m3" % ("%.6g" % v) for v in given)))
    out = OrderedDict()
    for m in used:
        out[m] = (props[m]["c"], rhos[m] if rhos[m] > 0.0 else 1.0)
    if len(used) > 1:
        log("INFO: Akoestiek: {0} media ({1}) met gelijke (of zonder) "
            "dichtheid; op het "
            "grensvlak is de druk en dp/dn continu. Mecway zelf rekent "
            "akoestiek met maar een materiaal.".format(
                len(used), ", ".join(used)))
    return out


def _free_shift(model_data, elsets, mats):
    """
    Zoekfrequentie [Hz] < 0 voor harde wanden zonder drukknoop: -0,1 x
    c_min / (2 D), D = diagonaal van het omhullende blok van de volumes.
    Voor een convex volume ligt de eerste mode met druk-variatie boven
    c / (2 D) (Payne-Weinberger, eigenwaarde >= pi^2 / D^2); de shift ligt
    dus ruim onder die mode en dicht bij de 0-Hz-mode. Gemeten: een buis
    van 1 m met CENTRE -10 Hz gaf dezelfde modes als Mecway (verificatie
    AK-01).
    """
    keep = set(es for es, _m in elsets)
    elem_elset = model_data.get("_elem_elset", {})
    nodes = model_data.get("mesh_nodes", {})
    lo, hi = [1.0e300] * 3, [-1.0e300] * 3
    for eid, _shape, enodes, elset in model_data.get("mesh_elements", []):
        if elem_elset.get(eid, elset) not in keep:
            continue
        for n in enodes:
            p = nodes.get(n)
            if p is None:
                continue
            for k in range(3):
                lo[k] = min(lo[k], p[k])
                hi[k] = max(hi[k], p[k])
    diag = math.sqrt(sum((hi[k] - lo[k]) ** 2 for k in range(3)))
    c_min = min(c for c, _rho in mats.values())
    if diag <= 0.0:
        return -1.0
    return -0.1 * c_min / (2.0 * diag)


def _calc_modes_block(model_data, analysis, shift_free):
    """res = CALC_MODES(...) zoals bridge_to_aster._calc_modes_block
    (modaal), zonder NORM_MODE (werkt niet voor MODE_ACOU, zie de kop)."""
    import bridge_to_aster as bta
    n = bta._mode_count(model_data, analysis)
    shift = float(analysis.get("modal_shift_hz", 0.0) or 0.0)
    lines = ["res = CALC_MODES(",
             "    TYPE_RESU='DYNAMIQUE',",
             "    MATR_RIGI=rigidite,",
             "    MATR_MASS=masse,"]
    if shift == 0.0 and shift_free is not None:
        shift = shift_free
        log("INFO: Akoestiek: harde wanden zonder drukknoop; de eerste "
            "mode is de constante druk (0 Hz). Zoekfrequentie {0} Hz "
            "(OPTION 'CENTRE').".format("%.6g" % shift))
    if shift != 0.0:
        lines += ["    OPTION='CENTRE',",
                  "    CALC_FREQ=_F(FREQ={0}, NMAX_FREQ={1}),".format(
                      _r(shift), n)]
    else:
        lines += ["    OPTION='PLUS_PETITE',",
                  "    CALC_FREQ=_F(NMAX_FREQ={0}),".format(n)]
    lines += ["    VERI_MODE=_F(SEUIL={0}, STOP_ERREUR='NON'),".format(
                  "%.10g" % float(analysis.get("modal_tolerance", 1.0e-6))),
              "    SOLVEUR=_F(METHODE='MUMPS', NPREC={0}),".format(
                  int(analysis.get("solver_nprec", 8))),
              ")"]
    log("INFO: CALC_MODES akoestisch: {0} modes.".format(n))
    return "\n".join(lines) + "\n"


def build_acoustic_modes(work_dir, model_data, msh_path, resu_name, analysis):
    """
    Genereer study.comm en study.export voor de holtemodes (Mecway A30).
    Returns het pad naar study.export, of None bij een fout.
    """
    import bridge_to_aster as bta
    log("--- Generating Aster Files (CALC_MODES akoestisch) ---")
    if CONFIG_BRIDGE.get("result_route", "liml") != "liml":
        log(tr("ERROR: De resultaatroute 'tabel' werkt niet voor modale "
               "analyses. Kies 'liml' (pagina Algemeen).",
               "ERROR: The result route 'tabel' does not work for modal "
               "analyses. Choose 'liml' (General page)."))
        return None
    ac = model_data.get("acoustic_data")
    if ac is None:
        log(tr("ERROR: Geen akoestische gegevens gelezen; is het Mecway-"
               "analysetype A30 (Acoustic Resonance 3D) of de actieve "
               "analyse 'acoustic_modes_3d'?",
               "ERROR: No acoustic data read; is the Mecway analysis type "
               "A30 (Acoustic Resonance 3D) or the active analysis "
               "'acoustic_modes_3d'?"))
        return None
    bad = list(ac["bad"])
    # Mecway zelf weigert dit ook ("Pressure must be zero.").
    bad += [tr("fluid pressure op '{0}' is {1}; alleen 0 kan (drukknoop; "
               "als bron: de harmonische analyse '{2}')",
               "fluid pressure on '{0}' is {1}; only 0 is possible "
               "(pressure node; as a source: the harmonic analysis "
               "'{2}')").format(sel, text, HARM_ANALYSIS)
            for sel, _k, _v, text in ac.get("pressure_src", [])]
    if bad:
        log(tr("ERROR: Akoestiek: {0}.", "ERROR: Acoustics: {0}.").format(
            "; ".join(bad)))
        return None
    # (geen list comprehensions met losse namen: in Python 2.7 lekken die
    # naar de functie, bijv. _r over de module-functie _r)
    named = list(item[0] for role in ("source", "absorbing")
                 for item in ac.get("named", {}).get(role, []))
    named += list(item[0] for item in ac.get("aster_sources", []))
    if named:
        log("INFO: Akoestiek: bron- en absorberende selecties ({0}) tellen "
            "alleen in de harmonische analyse '{1}'; bij de holtemodes zijn "
            "het harde wanden.".format(", ".join(named), HARM_ANALYSIS))
    if model_data.get("analysis_type") == "A20":
        log(tr("ERROR: Acoustic Resonance 2D (A20) wordt niet vertaald; de "
               "bridge rekent alleen 3D (A30).",
               "ERROR: Acoustic Resonance 2D (A20) is not translated; the "
               "bridge only computes 3D (A30)."))
        return None
    elsets = _elsets(model_data)
    if elsets is None:
        return None
    mats = _materials(elsets, ac)
    if mats is None:
        return None
    if ac["impedances"]:
        # Impedantie: complexe modes (S2b), eigen route (zie
        # _complex_modes_comm en aster_acoustic_complex.py).
        return _complex_modes_comm(work_dir, model_data, msh_path, resu_name,
                                   analysis, ac, elsets, mats)

    # --- model, materialen ---
    model_block = ("model = AFFE_MODELE(\n    MAILLAGE=mesh,\n    AFFE=_F(\n"
                   "        GROUP_MA={0},\n        PHENOMENE='ACOUSTIQUE',\n"
                   "        MODELISATION='3D',\n    ),\n)\n").format(
                       _fmt_groups([es for es, _m in elsets]))
    names = OrderedDict()
    mat_lines = []
    for k, (mat, (c, rho)) in enumerate(mats.items()):
        names[mat] = "MAC{0}".format(k + 1)
        mat_lines.append("{0}{1} = DEFI_MATERIAU(FLUIDE=_F(RHO={2}, "
                         "CELE_R={3}))".format(bta._mat_comment(mat),
                                               names[mat], _r(rho), _r(c)))
    fieldmat_block = ("fieldmat = AFFE_MATERIAU(\n    MAILLAGE=mesh,\n"
                      "    MODELE=model,\n    AFFE=(\n{0}    ),\n)\n").format(
                          "".join("        _F(GROUP_MA='{0}', MATER={1}),\n"
                                  .format(es, names[m]) for es, m in elsets))

    # --- drukknopen ---
    renames = model_data.get("selection_renames",
                             model_data.get("group_renames", {}))
    known = model_data.get("mesh_group_names", set())
    entries = []
    for sel, kind in ac["pressure0"]:
        name = renames.get(sel, sel)
        if name not in known:
            log(tr("WARNING: Akoestiek: selectie '{0}' van fluid pressure is "
                   "leeg; overgeslagen.",
                   "WARNING: Acoustics: selection '{0}' of fluid pressure is "
                   "empty; skipped.").format(sel))
            continue
        entries.append("_F({0}='{1}', PRES=0j)".format(
            "GROUP_NO" if kind == "node" else "GROUP_MA", name))
    charges = []
    loads_block = ""
    if entries:
        loads_block = ("# Drukknopen (Mecway fluid pressure = 0)\n"
                       "bc_pres = AFFE_CHAR_ACOU(\n    MODELE=model,\n"
                       "    PRES_IMPO=(\n        {0},\n    ),\n)\n").format(
                           ",\n        ".join(entries))
        charges.append("bc_pres")
    extra = model_data.get("aster_custom", {}).get("excit", [])
    charges += [n for n in extra if n not in charges]
    loads_block += bta._custom_code(model_data, "code")
    if not loads_block:
        loads_block = "#-- geen drukknopen (harde wanden) --\n"
    charge_ref = ("    CHARGE=({0},),\n".format(", ".join(charges))
                  if charges else "")
    assemblage_block = (
        "ASSEMBLAGE(\n"
        "    MODELE=model,\n"
        "    CHAM_MATER=fieldmat,\n"
        "{0}"
        "    NUME_DDL=CO('numedd'),\n"
        "    MATR_ASSE=(\n"
        "        _F(MATRICE=CO('rigi_c'), OPTION='RIGI_ACOU'),\n"
        "        _F(MATRICE=CO('mass_c'), OPTION='MASS_ACOU'),\n"
        "    ),\n"
        ")\n"
        "# Reele delen (de akoestische matrices zijn complex; ahlv100a)\n"
        "rigidite = COMB_MATR_ASSE(COMB_R=_F(MATR_ASSE=rigi_c, "
        "PARTIE='REEL', COEF_R=1.0))\n"
        "masse = COMB_MATR_ASSE(COMB_R=_F(MATR_ASSE=mass_c, "
        "PARTIE='REEL', COEF_R=1.0))\n").format(charge_ref)
    shift_free = None
    if not charges:
        shift_free = _free_shift(model_data, elsets, mats)
    calc_modes_block = _calc_modes_block(model_data, analysis, shift_free)
    export_block = (bta._custom_code(model_data, "post") +
                    _export_block(work_dir, model_data))
    if CONFIG_BRIDGE.get("med_export", False):
        impr_resu_block = ("IMPR_RESU(FORMAT='MED', UNITE=80, "
                           "RESU=_F(RESULTAT=res, NOM_CHAM=('PRES',)))\n")
    else:
        impr_resu_block = "#-- geen MED-export (instelling med_export uit) --\n"

    comm_path = os.path.join(work_dir, run_name(".comm"))
    try:
        template = bta._read_template(analysis["template"])
        bta._write_comm(comm_path, template.format(
            model_block=model_block,
            materials_block="\n".join(mat_lines) + "\n",
            fieldmat_block=fieldmat_block,
            loads_block=loads_block,
            assemblage_block=assemblage_block,
            calc_modes_block=calc_modes_block,
            export_block=export_block,
            impr_resu_block=impr_resu_block,
        ))
        log("SUCCESS: .comm geschreven naar '{0}'.".format(comm_path))
        export_path = os.path.join(work_dir, run_name(".export"))
        with open(export_path, "w") as f:
            f.write("\n".join([
                "P actions make_etude",
                "P time_limit {0}".format(analysis.get("time_limit", 900)),
                "P memory_limit {0}".format(analysis.get("memory_limit", 2048)),
                "P mpi_nbcpu 1",
                "P ncpus {0}".format(CONFIG_BRIDGE["ncpus"]),
                "P no-mpi true",
                "F comm {0} D 1".format(comm_path),
                "F mail {0} D 20".format(msh_path),
                "F mess {0} R 6".format(os.path.join(work_dir,
                                                     run_name(".mess"))),
            ]) + bta._export_result_lines(work_dir, resu_name))
        log("SUCCESS: .export geschreven naar '{0}'.".format(export_path))
        return export_path
    except KeyError as e:
        log(tr("ERROR: Template placeholder ontbreekt: {0}.",
               "ERROR: Template placeholder missing: {0}.").format(str(e)))
    except Exception as e:
        log(tr("ERROR: Genereren akoestische solver-bestanden mislukt: {0}",
               "ERROR: Generating the acoustic solver files failed: {0}")
            .format(str(e)))
    return None


def _export_block(work_dir, model_data):
    """aster_results_export: modes als <mode> met de knoopvariabele
    'pressure' (Mecway's naam), genormeerd p^T M p = 1 (export_results
    doet dat voor MODE_ACOU met MASS_GENE), en de BRIDGE_MODE-regels."""
    import bridge_to_aster as bta
    scripts = os.path.dirname(os.path.abspath(__file__))
    return ("# Resultaten in Mecway's eigen formaat (result_route 'liml')\n"
            "import sys\n"
            "sys.path.insert(0, {0})\n"
            "import aster_results_export\n"
            "aster_results_export.print_modes(res)\n"
            "aster_results_export.export_results(\n"
            "    res, {1}, {2},\n"
            "    fields=['pressure'], precision='double')\n").format(
                bta._py_str(scripts), bta._py_str(bta.results_file(work_dir)),
                bta._py_str(model_data.get("liml_path", "")))


# ==========================================
# COMPLEXE HOLTEMODES MET IMPEDANTIE (S2b)
#
# Mecway A30 met <impedance>: kwadratisch eigenwaardeprobleem (K + lambda C
# + lambda^2 M) P = 0 met C = int rho/Z N N^T dA (manual v33 par. 6.9.2).
# CALC_MODES kan dat niet (MATR_AMOR alleen matr_asse_depl_r/gene_r,
# catalogus calc_modes.py v17.4; AFFE_CHAR_ACOU kent geen impedantie;
# een complexe Z geeft een complexe C). Route:
#   - AFFE_MODELE: volumes '3D' en de impedantievlakken '3D_ABSO';
#   - ASSEMBLAGE RIGI_ACOU / MASS_ACOU zonder CHARGE (de drukknopen
#     elimineert aster_acoustic_complex zelf);
#   - per impedantie CALC_MATR_ELEM AMOR_ACOU met een eigen CHAM_MATER:
#     CELE_R = 1 op dat vlak (de 3D_ABSO-matrix is (1/c) int N N^T dA, zonder
#     rho, gemeten in S3a; met c = 1 het oppervlaktemassa-achtige int N N^T
#     dA, som = oppervlak, gemeten), _CELE_OTHER op de andere
#     impedantievlakken (bijdrage 1e-30: verwaarloosbaar);
#   - aster_acoustic_complex.solve_modes (scipy, in de job): C = sum
#     rho/Z_g A_g, linearisering + shift-invert, export als Mecway's
#     complexe modes (frequency + decayrate, pressure + pressurei).
# Gemeten (buis hard / Z, hex20; verificatie AI): gelijk aan Mecway 33
# headless tot 1e-9 en aan de analytische wortels.
# ==========================================
_CELE_OTHER = 1.0e30


def _impedance_inputs(model_data, ac, mats):
    """
    Impedantieranden voor de complexe modes: [(Code-Aster-groep, rho / Z
    (complex), omschrijving, Mecway-selectie)], of None (ERROR gelogd).
    Mecway's Z = specificimpedancer + i specificimpedancei (Pa.s/m, getal
    zonder eenheid); headless gemeten: Z = 2000 - 1000i geeft dezelfde
    eigenwaarden als C = rho/Z int N N dA in een eigen berekening.
    """
    renames = model_data.get("selection_renames",
                             model_data.get("group_renames", {}))
    known = model_data.get("mesh_group_names", set())
    faces = model_data.get("face_selection_faces", {})
    out, bad = [], []
    for sel, zr, zi in ac.get("impedance_values", []):
        if sel not in faces:
            bad.append(tr("impedantie op '{0}': kies een vlakselectie",
                          "impedance on '{0}': choose a face "
                          "selection").format(sel))
            continue
        if zr is None or zi is None:
            bad.append(tr("impedantie op '{0}': Re Z en Im Z moeten getallen "
                          "zijn (Pa.s/m, zonder eenheid)",
                          "impedance on '{0}': Re Z and Im Z must be numbers "
                          "(Pa.s/m, without unit)").format(sel))
            continue
        z = complex(zr, zi)
        if abs(z) == 0.0:
            bad.append(tr("impedantie op '{0}' is 0; gebruik fluid pressure "
                          "= 0 voor een drukknoop (open einde)",
                          "impedance on '{0}' is 0; use fluid pressure = 0 "
                          "for a pressure node (open end)").format(sel))
            continue
        media = _face_media(model_data, sel, mats)
        if len(media) != 1:
            bad.append(tr("impedantie '{0}' ligt aan {1} media ({2}); maak "
                          "per medium een eigen selectie",
                          "impedance '{0}' borders {1} media ({2}); make a "
                          "separate selection per medium").format(
                              sel, len(media), ", ".join(media) or "-"))
            continue
        mat = media[0]
        rho = ac["materials"].get(mat, {}).get("rho") or 0.0
        if rho <= 0.0:
            # Mecway rekent zonder dichtheid met rho = 0: de impedantie
            # werkt dan niet (headless gezien: Z = 2000 zonder dichtheid gaf
            # de modes van een harde wand, decay rate ~1e-10).
            bad.append(tr("impedantie op '{0}': materiaal '{1}' heeft geen "
                          "dichtheid; de randterm is rho / Z (Mecway rekent "
                          "zonder dichtheid met rho = 0: een harde wand)",
                          "impedance on '{0}': material '{1}' has no "
                          "density; the boundary term is rho / Z (without a "
                          "density Mecway uses rho = 0: a hard wall)").format(
                              sel, mat))
            continue
        name = renames.get(sel, sel)
        if name not in known:
            log(tr("WARNING: Akoestiek: impedantieselectie '{0}' is leeg; "
                   "overgeslagen.",
                   "WARNING: Acoustics: impedance selection '{0}' is empty; "
                   "skipped.").format(sel))
            continue
        c = mats[mat][0]
        ratio = z / (rho * c)
        refl = (z - rho * c) / (z + rho * c)
        log("INFO: Akoestiek: impedantie Z = {0} {1} {2}i Pa.s/m op '{3}' "
            "(Z / (rho c) = {4} {5} {6}i; reflectiefactor bij loodrechte "
            "inval |R| = {7}).".format(
                "%.6g" % z.real, "-" if z.imag < 0 else "+",
                "%.6g" % abs(z.imag), sel, "%.4g" % ratio.real,
                "-" if ratio.imag < 0 else "+", "%.4g" % abs(ratio.imag),
                "%.4g" % abs(refl)))
        if z.real < 0.0:
            log(tr("WARNING: Akoestiek: impedantie op '{0}' heeft Re Z < 0: "
                   "de rand voert energie toe; modes groeien aan (negatieve "
                   "decay rate).",
                   "WARNING: Acoustics: impedance on '{0}' has Re Z < 0: "
                   "the boundary adds energy; modes grow (negative decay "
                   "rate).").format(sel))
        out.append((name, rho / z, "{0}: Z = {1}".format(
            sel, ("%.6g" % z.real) + ("%+.6gi" % z.imag if z.imag else "")),
            sel, mat))
    if bad:
        log(tr("ERROR: Akoestiek: {0}.", "ERROR: Acoustics: {0}.").format(
            "; ".join(bad)))
        return None
    if not out:
        log(tr("ERROR: Akoestiek: alle impedantieselecties zijn leeg.",
               "ERROR: Acoustics: all impedance selections are empty."))
        return None
    return out


def _cplx(z):
    """complex(..) als Python-literal voor de .comm."""
    return "complex({0}, {1})".format(_r(z.real), _r(z.imag))


def _complex_modes_comm(work_dir, model_data, msh_path, resu_name, analysis,
                        ac, elsets, mats):
    """
    .comm en .export voor de complexe holtemodes (impedantie, zie het blok
    hierboven). Zelfde template als de reele holtemodes
    (AsterTemplate_Akoestiek.comm); de 'CALC_MODES'-sectie roept
    aster_acoustic_complex aan. Returns het pad naar de .export of None.
    """
    import bridge_to_aster as bta
    imps = _impedance_inputs(model_data, ac, mats)
    if imps is None:
        return None
    vol = [es for es, _m in elsets]
    faces = []
    for g, _coef, _label, _sel, _mat in imps:
        if g not in faces:
            faces.append(g)
    model_block = ("model = AFFE_MODELE(\n    MAILLAGE=mesh,\n    AFFE=(\n"
                   "        _F(GROUP_MA={0}, PHENOMENE='ACOUSTIQUE', "
                   "MODELISATION='3D'),\n"
                   "        _F(GROUP_MA={1}, PHENOMENE='ACOUSTIQUE', "
                   "MODELISATION='3D_ABSO'),\n    ),\n)\n").format(
                       _fmt_groups(vol), _fmt_groups(faces))
    names = OrderedDict()
    mat_lines = []
    for k, (mat, (c, rho)) in enumerate(mats.items()):
        names[mat] = "MAC{0}".format(k + 1)
        mat_lines.append("{0}{1} = DEFI_MATERIAU(FLUIDE=_F(RHO={2}, "
                         "CELE_R={3}))".format(bta._mat_comment(mat),
                                               names[mat], _r(rho), _r(c)))
    mat_lines.append(
        "# Impedantievlakken: CELE_R = 1 geeft AMOR_ACOU = int N N dA "
        "(bridge_acoustic,\n# complexe modes); de andere impedantievlakken "
        "in dezelfde matrix ~0.\n"
        "MAZ1 = DEFI_MATERIAU(FLUIDE=_F(RHO=1.0, CELE_R=1.0))\n"
        "MAZ0 = DEFI_MATERIAU(FLUIDE=_F(RHO=1.0, CELE_R={0}))".format(
            _r(_CELE_OTHER)))
    face_mat = dict((g, names[mat]) for g, _c, _l, _s, mat in imps)
    vol_assign = "".join("        _F(GROUP_MA='{0}', MATER={1}),\n".format(
        es, names[m]) for es, m in elsets)
    fieldmat_block = ("fieldmat = AFFE_MATERIAU(\n    MAILLAGE=mesh,\n"
                      "    MODELE=model,\n    AFFE=(\n{0}{1}    ),\n)\n").format(
                          vol_assign, "".join(
                              "        _F(GROUP_MA='{0}', MATER={1}),\n"
                              .format(g, face_mat[g]) for g in faces))

    # --- drukknopen: in aster_acoustic_complex geelimineerd ---
    dirichlet = []
    renames = model_data.get("selection_renames",
                             model_data.get("group_renames", {}))
    known = model_data.get("mesh_group_names", set())
    for sel, kind in ac["pressure0"]:
        name = renames.get(sel, sel)
        if name not in known:
            log(tr("WARNING: Akoestiek: selectie '{0}' van fluid pressure is "
                   "leeg; overgeslagen.",
                   "WARNING: Acoustics: selection '{0}' of fluid pressure is "
                   "empty; skipped.").format(sel))
            continue
        dirichlet.append(("GROUP_NO" if kind == "node" else "GROUP_MA", name))
    extra = model_data.get("aster_custom", {}).get("excit", [])
    if extra:
        log(tr("WARNING: Akoestiek: #ASTER-EXCIT ({0}) wordt bij complexe "
               "modes (impedantie) niet gebruikt: de matrices worden zonder "
               "lasten geassembleerd.",
               "WARNING: Acoustics: #ASTER-EXCIT ({0}) is not used for "
               "complex modes (impedance): the matrices are assembled "
               "without loads.").format(", ".join(extra)))
    loads_block = ("#-- drukknopen: in aster_acoustic_complex geelimineerd "
                   "({0}) --\n".format(", ".join(n for _k, n in dirichlet)
                                       or "geen") +
                   bta._custom_code(model_data, "code"))

    lines = ["ASSEMBLAGE(",
             "    MODELE=model,",
             "    CHAM_MATER=fieldmat,",
             "    NUME_DDL=CO('numedd'),",
             "    MATR_ASSE=(",
             "        _F(MATRICE=CO('rigi_c'), OPTION='RIGI_ACOU'),",
             "        _F(MATRICE=CO('mass_c'), OPTION='MASS_ACOU'),",
             "    ),",
             ")"]
    damping = []
    for k, (g, coef, label, _sel, _mat) in enumerate(imps):
        idx = k + 1
        lines += ["# Impedantie {0}".format(label),
                  "chz{0} = AFFE_MATERIAU(".format(idx),
                  "    MAILLAGE=mesh,",
                  "    MODELE=model,",
                  "    AFFE=(",
                  vol_assign.rstrip("\n")]
        lines += ["        _F(GROUP_MA='{0}', MATER={1}),".format(
            f, "MAZ1" if f == g else "MAZ0") for f in faces]
        lines += ["    ),",
                  ")",
                  "mez{0} = CALC_MATR_ELEM(OPTION='AMOR_ACOU', MODELE=model, "
                  "CHAM_MATER=chz{0})".format(idx),
                  "amz{0} = ASSE_MATRICE(MATR_ELEM=mez{0}, NUME_DDL=numedd)"
                  .format(idx)]
        damping.append("(amz{0}, {1}, {2})".format(idx, _cplx(coef),
                                                   bta._py_str(label)))
    assemblage_block = "\n".join(lines) + "\n"

    n = bta._mode_count(model_data, analysis)
    shift = float(analysis.get("modal_shift_hz", 0.0) or 0.0)
    free = _free_shift(model_data, elsets, mats)
    c_min = min(c for c, _rho in mats.values())
    size = -0.1 * c_min / (2.0 * free) if free < 0.0 else 1.0
    scripts = os.path.dirname(os.path.abspath(__file__))
    calc_modes_block = (
        "# Complexe modes (impedantie): kwadratisch eigenwaardeprobleem\n"
        "# (K + l C + l^2 M) P = 0, C = sum rho/Z int N N dA (Mecway manual\n"
        "# par. 6.9.2), opgelost met scipy in aster_acoustic_complex\n"
        "# (CALC_MODES neemt MATR_AMOR niet voor drukmatrices).\n"
        "import sys\n"
        "sys.path.insert(0, {0})\n"
        "import aster_acoustic_complex\n"
        "aster_acoustic_complex.solve_modes(\n"
        "    rigi_c, mass_c, [\n{1}    ],\n"
        "    numedd, mesh, {2}, {3},\n"
        "    n_modes={4}, dirichlet={5},\n"
        "    shift_hz={6}, free_shift_hz={7}, c_ref={8}, size_ref={9},\n"
        "    tolerance={10}, precision='double')\n").format(
            bta._py_str(scripts),
            "".join("        {0},\n".format(d) for d in damping),
            bta._py_str(bta.results_file(work_dir)),
            bta._py_str(model_data.get("liml_path", "")), n,
            repr(dirichlet), _r(shift), _r(free), _r(c_min), _r(size),
            "%.10g" % float(analysis.get("modal_tolerance", 1.0e-6)))
    log("INFO: Akoestiek: {0} impedantierand(en): complexe modes "
        "(frequentie en decay rate), {1} modes{2}.".format(
            len(imps), n, (" rond {0} Hz".format("%.6g" % shift)
                           if shift > 0.0 else " met de kleinste |lambda| "
                           "(zoals Mecway)")))
    if shift < 0.0:
        log(tr("WARNING: Akoestiek: een negatieve zoekfrequentie telt bij "
               "complexe modes niet; de modes met de kleinste |lambda| "
               "(zoals Mecway).",
               "WARNING: Acoustics: a negative search frequency does not "
               "apply to complex modes; the modes with the smallest "
               "|lambda| (like Mecway)."))
    export_block = bta._custom_code(model_data, "post")
    if CONFIG_BRIDGE.get("med_export", False):
        log(tr("WARNING: Akoestiek: geen MED-export bij complexe modes "
               "(impedantie): er is geen Code-Aster-resultaat; de modes "
               "staan alleen in Mecway.",
               "WARNING: Acoustics: no MED export for complex modes "
               "(impedance): there is no Code-Aster result; the modes are "
               "only in Mecway."))
    blocks = dict(model_block=model_block,
                  materials_block="\n".join(mat_lines) + "\n",
                  fieldmat_block=fieldmat_block, loads_block=loads_block,
                  assemblage_block=assemblage_block,
                  calc_modes_block=calc_modes_block,
                  export_block=export_block or "#-- geen #ASTER-POST --\n",
                  impr_resu_block="#-- geen MED-export (complexe modes: geen "
                                  "Code-Aster-resultaat) --\n")
    return _write_modes_files(work_dir, analysis, msh_path, resu_name, blocks,
                              med=False)


def _write_modes_files(work_dir, analysis, msh_path, resu_name, blocks,
                       med=True):
    """De .comm (template analysis["template"], placeholders blocks) en de
    .export; returns het pad naar de .export of None (ERROR gelogd). med
    False: geen .rmed in de .export (er wordt niets naar UNITE 80
    geschreven)."""
    import bridge_to_aster as bta
    result_lines = bta._export_result_lines(work_dir, resu_name)
    if not med:
        result_lines = "".join("\n" + line for line in
                               result_lines.split("\n")
                               if line and not line.startswith("F rmed"))
    comm_path = os.path.join(work_dir, run_name(".comm"))
    try:
        template = bta._read_template(analysis["template"])
        bta._write_comm(comm_path, template.format(**blocks))
        log("SUCCESS: .comm geschreven naar '{0}'.".format(comm_path))
        export_path = os.path.join(work_dir, run_name(".export"))
        with open(export_path, "w") as f:
            f.write("\n".join([
                "P actions make_etude",
                "P time_limit {0}".format(analysis.get("time_limit", 900)),
                "P memory_limit {0}".format(analysis.get("memory_limit", 2048)),
                "P mpi_nbcpu 1",
                "P ncpus {0}".format(CONFIG_BRIDGE["ncpus"]),
                "P no-mpi true",
                "F comm {0} D 1".format(comm_path),
                "F mail {0} D 20".format(msh_path),
                "F mess {0} R 6".format(os.path.join(work_dir,
                                                     run_name(".mess"))),
            ]) + result_lines)
        log("SUCCESS: .export geschreven naar '{0}'.".format(export_path))
        return export_path
    except KeyError as e:
        log(tr("ERROR: Template placeholder ontbreekt: {0}.",
               "ERROR: Template placeholder missing: {0}.").format(str(e)))
    except Exception as e:
        log(tr("ERROR: Genereren akoestische solver-bestanden mislukt: {0}",
               "ERROR: Generating the acoustic solver files failed: {0}")
            .format(str(e)))
    return None


# ==========================================
# HARMONISCH DRUKVELD (S3a): DYNA_VIBRA TYPE_CALCUL='HARM'
#
# Mecway heeft geen harmonische analyse (manual v33 hoofdstuk 6); de bridge
# rekent het A30-model (materialen, fluid pressure, impedance) met:
#   bron        vlakselectie BRON_... / SOURCE_... (normaalsnelheid
#               source_velocity, de fluid in) of #ASTER-BRON selectie =
#               v [, fase]  -> AFFE_CHAR_ACOU VITE_FACE VNOR
#               fluid pressure != 0  -> opgelegde druk (PRES_IMPO)
#   absorberend vlakselectie ABSORBEREND_... / ABSORBING_... -> 3D_ABSO
#               (rho c van het medium); Mecway impedance met reele Z
#               -> 3D_ABSO met een eigen materiaal CELE_R = Z / rho
#   microfoon   knoopselectie MICROFOON_... / MICROPHONE_... -> CSV en
#               BRIDGE_MIC-regels (aster_acoustic_export)
# Code-Aster-route: testcase ahlv100a (ASSEMBLAGE RIGI_ACOU, MASS_ACOU,
# AMOR_ACOU en VECT_ASSE CHAR_ACOU; DYNA_VIBRA HARM PHYS met EXCIT
# VECT_ASSE, PUIS_PULS=1, PHAS_DEG=90).
# Gemeten in runs (Code-Aster 17.4, buis 1 x 0,1 x 0,1 m, hex20):
#   - de matrices bevatten rho niet (zie de kop); de bron wel: VITE_FACE
#     VNOR = v met EXCIT COEF_MULT=+1, PUIS_PULS=1, PHAS_DEG=90 (= i w v)
#     gaf p = +rho c v op het bronvlak bij een absorberend einde; met
#     COEF_MULT=-1 (ahlv100a) -rho c v. +1 betekent dus: v de fluid in
#     (een zuiger die duwt), p(t) = Re(P e^(i w t)).
#   - VITE_FACE vraagt de vlakken in het model (MODELISATION '3D',
#     elementen ACOUNA_FACE*) met een materiaal (anders <EXCEPTION>
#     CALCUL_20, gezien).
#   - 3D_ABSO: AMOR_ACOU = (1/c) int N N dA met c van het materiaal van
#     het vlak, zonder rho (RHO 2,4 i.p.v. 1,2 op het vlak: geen verschil;
#     CELE_R 680 i.p.v. 340: reflectie 1/3, max/min van |p| = 2,000, zoals
#     een impedantie Z = 2 rho c). Een impedantie Z wordt dus CELE_R =
#     Z / rho van het medium (Mecway: C = int rho/Z N N dA, manual v33
#     par. 6.9.2). CELE_I (complexe c) gaf NIET de verwachte reflectie van
#     een complexe Z (gemeten 2,19 i.p.v. 2,62): complexe Z -> ERROR.
#   - PRES_IMPO != 0 in VECT_ASSE CHAR_ACOU (met de matrices geassembleerd
#     met CHARGE) geeft de opgelegde druk; als aparte EXCIT zonder
#     PUIS_PULS (p = 2 Pa: |p| = 2,000000 langs de buis, gezien).
# ==========================================
_EDGES = {
    "hex": ((0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7)),
    "wedge": ((0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3),
              (1, 4), (2, 5)),
    "tet": ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)),
}


def _frequencies(analysis, label=None):
    """Oplopende lijst frequenties [Hz] uit harm_freq_list (tekst, "0" =
    niet gebruikt) of harm_freq_start/_end/_count/_spacing; None bij een
    fout (ERROR gelogd). label: (nl, en) voor de melding (standaard
    akoestiek harmonisch; ook bridge_vibroacoustic_harm)."""
    nl, en = label or ("Akoestiek harmonisch", "Acoustics harmonic")
    text = "{0}".format(analysis.get("harm_freq_list", "0") or "0").strip()
    vals = []
    if text not in ("", "0", "0.0"):
        try:
            vals = [float(t) for t in re.split(r"[;,\s]+", text) if t]
        except ValueError:
            vals = None
        if not vals or min(vals) <= 0.0:
            log(tr("ERROR: {1}: frequentielijst "
                   "(harm_freq_list) '{0}' is ongeldig; verwacht getallen > 0 "
                   "in Hz gescheiden door komma's, bijv. '100, 250, 500' "
                   "(punt als decimaalteken), of 0 voor het bereik.",
                   "ERROR: {1}: frequency list "
                   "(harm_freq_list) '{0}' is invalid; expected numbers > 0 "
                   "in Hz separated by commas, e.g. '100, 250, 500' (point "
                   "as decimal separator), or 0 for the range.").format(
                       text, tr(nl, en)))
            return None
    else:
        try:
            f0 = float(analysis.get("harm_freq_start", 0.0))
            f1 = float(analysis.get("harm_freq_end", 0.0))
            n = int(analysis.get("harm_freq_count", 1))
        except (TypeError, ValueError):
            f0, f1, n = 0.0, 0.0, 0
        spacing = "{0}".format(analysis.get("harm_freq_spacing", "lin"))
        if f0 <= 0.0 or f1 < f0 or n < 1:
            log(tr("ERROR: {3}: frequentiebereik {0} .. {1} "
                   "Hz met {2} stappen is ongeldig (start > 0, eind >= "
                   "start, aantal >= 1).",
                   "ERROR: {3}: frequency range {0} .. {1} Hz "
                   "with {2} steps is invalid (start > 0, end >= start, "
                   "count >= 1).").format(f0, f1, n, tr(nl, en)))
            return None
        if n == 1 or f1 == f0:
            vals = [f0]
        elif spacing == "log":
            vals = [f0 * (f1 / f0) ** (k / float(n - 1)) for k in range(n)]
        else:
            vals = [f0 + (f1 - f0) * k / float(n - 1) for k in range(n)]
    out = []
    for v in sorted(vals):
        if not out or v > out[-1] * (1.0 + 1.0e-9):
            out.append(v)
    return out


def _face_media(model_data, sel, mats):
    """Mecway-materialen van de volume-elementen waar de vlakken van de
    selectie (Mecway-naam) op liggen."""
    elem_elset = model_data.get("_elem_elset", {})
    elsets = model_data.get("elsets", {})
    out = []
    for eid, _fid in model_data.get("face_selection_faces", {}).get(sel, []):
        mat = elsets.get(elem_elset.get(eid))
        if mat in mats and mat not in out:
            out.append(mat)
    return out


def _wavelength_check(model_data, elsets, mats, f_max):
    """WARNING als de grootste elementrand meer is dan lambda / N bij de
    hoogste frequentie (N = _MIN_PER_WAVELENGTH, lineair/kwadratisch);
    INFO met het kleinste aantal elementen per golflengte."""
    keep = dict(elsets)
    nodes = model_data.get("mesh_nodes", {})
    elem_elset = model_data.get("_elem_elset", {})
    worst = None
    n_bad = 0
    for eid, shape, enodes, elset in model_data.get("mesh_elements", []):
        name = elem_elset.get(eid, elset)
        mat = keep.get(name)
        if mat is None:
            continue
        edges = _EDGES.get(shape.rstrip("0123456789"))
        if not edges:
            continue
        h = 0.0
        for a, b in edges:
            pa, pb = nodes.get(enodes[a]), nodes.get(enodes[b])
            if pa is None or pb is None:
                continue
            h = max(h, math.sqrt(sum((pa[k] - pb[k]) ** 2
                                     for k in range(3))))
        if h <= 0.0:
            continue
        per = mats[mat][0] / f_max / h
        linear = shape in _LINEAR_SHAPES
        if per < _MIN_PER_WAVELENGTH[linear]:
            n_bad += 1
        if worst is None or per < worst[0]:
            worst = (per, linear)
    if worst is None:
        return
    log("INFO: Akoestiek harmonisch: kleinste aantal elementen per "
        "golflengte bij {0} Hz: {1} (grootste elementrand).".format(
            "%.6g" % f_max, "%.3g" % worst[0]))
    if n_bad:
        log(tr("WARNING: Akoestiek harmonisch: {0} element(en) hebben bij de "
               "hoogste frequentie ({1} Hz) minder dan {2} (lineair) / {3} "
               "(kwadratisch) elementen per golflengte (kleinste: {4}). De "
               "druk en vooral de fase worden dan onnauwkeurig; verfijn het "
               "mesh of verlaag de hoogste frequentie (docs/Akoestiek.md).",
               "WARNING: Acoustics harmonic: {0} element(s) have fewer than "
               "{2} (linear) / {3} (quadratic) elements per wavelength at "
               "the highest frequency ({1} Hz) (smallest: {4}). The pressure "
               "and especially the phase become inaccurate; refine the mesh "
               "or lower the highest frequency (docs/en/Acoustics.md).")
            .format(n_bad, "%.6g" % f_max,
                    "%g" % _MIN_PER_WAVELENGTH[True],
                    "%g" % _MIN_PER_WAVELENGTH[False], "%.3g" % worst[0]))


def _ri(z):
    """Complex getal als ('RI', re, im) voor de .comm (ahlv100a)."""
    return "('RI', {0}, {1})".format(_r(z.real), _r(z.imag))


def _harm_inputs(model_data, analysis, ac, mats):
    """
    Bronnen, absorberende vlakken en microfoons naar Code-Aster-namen.
    Returns dict of None (ERROR gelogd):
      "vite":  [(groep, complexe v [m/s], omschrijving)]
      "pres":  [(GROUP_NO/GROUP_MA, groep, complexe p [Pa])]  (ook p = 0)
      "abso":  [(groep, c [m/s], rho, omschrijving)]
      "mics":  OrderedDict {Mecway-naam: [knoopnummers]}
      "source_media": {groep: Mecway-materiaal} (materiaal op het bronvlak)
    """
    renames = model_data.get("selection_renames",
                             model_data.get("group_renames", {}))
    known = model_data.get("mesh_group_names", set())
    faces = model_data.get("face_selection_faces", {})
    named = ac.get("named", {})
    out = {"vite": [], "pres": [], "abso": [], "mics": OrderedDict(),
           "source_media": OrderedDict(), "abso_media": OrderedDict()}
    bad = []
    roles = {}

    def aster(sel):
        name = renames.get(sel, sel)
        return name if name in known else None

    def claim(sel, role):
        if sel in roles and roles[sel] != role:
            bad.append(tr("selectie '{0}' is zowel {1} als {2}",
                          "selection '{0}' is both {1} and {2}").format(
                              sel, roles[sel], role))
        roles[sel] = role

    def one_medium(sel, what):
        media = _face_media(model_data, sel, mats)
        if len(media) != 1:
            bad.append(tr(
                "{0} '{1}' ligt aan {2} media ({3}); maak per medium een "
                "eigen selectie",
                "{0} '{1}' borders {2} media ({3}); make a separate "
                "selection per medium").format(
                    what, sel, len(media), ", ".join(media) or "-"))
            return None
        return media[0]

    # --- bronnen: snelheid ---
    v0 = float(analysis.get("source_velocity", 1.0e-3))
    how = tr("BRON_-selectie", "SOURCE_ selection")
    vel = [(sel, complex(v0, 0.0), how)
           for sel, kind in named.get("source", []) if kind == "face"]
    for sel, kind in named.get("source", []):
        if kind != "face":
            bad.append(tr("bron '{0}' is een {1}selectie; een bron moet "
                          "een vlakselectie zijn",
                          "source '{0}' is a {1} selection; a source must "
                          "be a face selection").format(sel, kind))
    for sel, v, phase in ac.get("aster_sources", []):
        if sel not in faces:
            bad.append(tr("#ASTER-BRON: '{0}' is geen vlakselectie",
                          "#ASTER-BRON: '{0}' is not a face "
                          "selection").format(sel))
            continue
        rad = math.radians(phase)
        # cos(90 graden) = 6e-17: afronden, zodat 90 graden zuiver imaginair is
        re, im = [0.0 if abs(t) < 1.0e-12 else t
                  for t in (math.cos(rad), math.sin(rad))]
        vel = [x for x in vel if x[0] != sel]
        vel.append((sel, complex(v * re, v * im), "#ASTER-BRON"))
    for sel, z, how in vel:
        claim(sel, tr("bron", "source"))
        name = aster(sel)
        if name is None:
            log(tr("WARNING: Akoestiek: bronselectie '{0}' is leeg; "
                   "overgeslagen.",
                   "WARNING: Acoustics: source selection '{0}' is empty; "
                   "skipped.").format(sel))
            continue
        mat = one_medium(sel, tr("bron", "source"))
        if mat is None:
            continue
        out["vite"].append((name, z, "{0} {1}".format(how, sel)))
        out["source_media"][name] = mat
        log("INFO: Akoestiek: bron '{0}' ({1}): normaalsnelheid {2} m/s, "
            "fase {3} graden (de fluid in).".format(
                sel, how, "%.6g" % abs(z),
                "%.6g" % math.degrees(math.atan2(z.imag, z.real))))

    # --- opgelegde druk (fluid pressure; 0 = drukknoop) ---
    for sel, kind in ac.get("pressure0", []):
        name = aster(sel)
        if name is None:
            log(tr("WARNING: Akoestiek: selectie '{0}' van fluid pressure is "
                   "leeg; overgeslagen.",
                   "WARNING: Acoustics: selection '{0}' of fluid pressure is "
                   "empty; skipped.").format(sel))
            continue
        out["pres"].append(("GROUP_NO" if kind == "node" else "GROUP_MA",
                            name, 0j))
    for sel, kind, value, _text in ac.get("pressure_src", []):
        claim(sel, tr("opgelegde druk", "imposed pressure"))
        name = aster(sel)
        if name is None:
            continue
        out["pres"].append(("GROUP_NO" if kind == "node" else "GROUP_MA",
                            name, complex(value, 0.0)))
        log("INFO: Akoestiek: opgelegde druk {0} Pa op '{1}' (fluid "
            "pressure; bron).".format("%.6g" % value, sel))

    # --- absorberende randen ---
    for sel, kind in named.get("absorbing", []):
        if kind != "face":
            bad.append(tr("absorberende rand '{0}' is een {1}selectie; "
                          "kies een vlakselectie",
                          "absorbing boundary '{0}' is a {1} selection; "
                          "choose a face selection").format(sel, kind))
            continue
        claim(sel, tr("absorberend", "absorbing"))
        name = aster(sel)
        mat = one_medium(sel, tr("absorberende rand", "absorbing boundary"))
        if name is None or mat is None:
            continue
        c, rho = mats[mat]
        out["abso"].append((name, c, rho, "rho c ({0})".format(mat)))
        out["abso_media"][name] = mat
    for sel, zr, zi in ac.get("impedance_values", []):
        claim(sel, tr("impedantie", "impedance"))
        name = aster(sel)
        if sel not in faces:
            bad.append(tr("impedantie op '{0}': kies een vlakselectie",
                          "impedance on '{0}': choose a face "
                          "selection").format(sel))
            continue
        if zr is None or zi is None or zr <= 0.0:
            bad.append(tr("impedantie op '{0}': Re Z moet een getal > 0 zijn "
                          "(Pa.s/m)",
                          "impedance on '{0}': Re Z must be a number > 0 "
                          "(Pa.s/m)").format(sel))
            continue
        if zi != 0.0:
            bad.append(tr(
                "impedantie op '{0}' is complex (Im Z = {1}); de bridge "
                "vertaalt in de harmonische analyse alleen een reele Z "
                "(KNOWN_BUGS 57)",
                "impedance on '{0}' is complex (Im Z = {1}); in the "
                "harmonic analysis the bridge only translates a real Z "
                "(KNOWN_BUGS 57)").format(sel, "%g" % zi))
            continue
        mat = one_medium(sel, tr("impedantie", "impedance"))
        if name is None or mat is None:
            continue
        rho_given = ac["materials"].get(mat, {}).get("rho") or 0.0
        if rho_given <= 0.0:
            bad.append(tr(
                "impedantie op '{0}': materiaal '{1}' heeft geen dichtheid "
                "(nodig voor Z / rho, zoals in Mecway)",
                "impedance on '{0}': material '{1}' has no density (needed "
                "for Z / rho, as in Mecway)").format(sel, mat))
            continue
        out["abso"].append((name, zr / rho_given, rho_given,
                            "Z = {0} Pa.s/m ({1})".format("%g" % zr, mat)))
        out["abso_media"][name] = mat
        log("INFO: Akoestiek: impedantie Z = {0} Pa.s/m op '{1}' "
            "(Z / (rho c) = {2}): absorberende rand 3D_ABSO met c = Z / rho "
            "= {3} m/s.".format("%g" % zr, sel,
                                "%.4g" % (zr / (rho_given * mats[mat][0])),
                                "%.6g" % (zr / rho_given)))

    # --- microfoons ---
    for sel, kind in named.get("microphone", []):
        if kind != "node":
            log(tr("WARNING: Akoestiek: microfoon '{0}' is een {1}selectie; "
                   "alleen knoopselecties tellen; overgeslagen.",
                   "WARNING: Acoustics: microphone '{0}' is a {1} "
                   "selection; only node selections count; "
                   "skipped.").format(sel, kind))
            continue
        nids = [int(str(n).lstrip("N")) for n in
                model_data.get("node_groups", {}).get(sel, [])]
        if nids:
            out["mics"][sel] = nids
    if bad:
        log(tr("ERROR: Akoestiek harmonisch: {0}.",
               "ERROR: Acoustics harmonic: {0}.").format("; ".join(bad)))
        return None
    if not out["vite"] and not any(abs(p) > 0.0 for _g, _n, p
                                   in out["pres"]):
        log(tr("ERROR: Akoestiek harmonisch: geen bron. Geef een "
               "vlakselectie de naam BRON_... (normaalsnelheid "
               "source_velocity), gebruik '#ASTER-BRON selectie = "
               "snelheid', of een fluid pressure != 0 (opgelegde druk).",
               "ERROR: Acoustics harmonic: no source. Name a face "
               "selection SOURCE_... (normal velocity source_velocity), "
               "use '#ASTER-BRON selection = velocity', or a fluid "
               "pressure != 0 (imposed pressure)."))
        return None
    if not out["abso"]:
        log(tr("WARNING: Akoestiek harmonisch: geen absorberende rand "
               "(ABSORBEREND_-selectie of impedance): zonder demping is "
               "de respons bij een holtemode oneindig groot. Kies "
               "frequenties tussen de resonanties of voeg een "
               "absorberende rand toe.",
               "WARNING: Acoustics harmonic: no absorbing boundary "
               "(ABSORBING_ selection or impedance): without damping the "
               "response at a cavity mode is infinite. Choose frequencies "
               "between the resonances or add an absorbing boundary."))
    return out


def build_acoustic_harm(work_dir, model_data, msh_path, resu_name, analysis):
    """
    Genereer de .comm en .export voor het harmonische drukveld
    (acoustic_harm_3d; bridge_to_aster._build_dyna_modal stuurt een
    akoestische DYNA_VIBRA-analyse hierheen). Returns het pad naar de
    .export, of None bij een fout.
    """
    import bridge_to_aster as bta
    log("--- Generating Aster Files (DYNA_VIBRA HARM akoestisch) ---")
    if CONFIG_BRIDGE.get("result_route", "liml") != "liml":
        log(tr("ERROR: De resultaatroute 'tabel' werkt niet voor de "
               "harmonische akoestiek. Kies 'liml' (pagina Algemeen).",
               "ERROR: The result route 'tabel' does not work for harmonic "
               "acoustics. Choose 'liml' (General page)."))
        return None
    ac = model_data.get("acoustic_data")
    if ac is None:
        log(tr("ERROR: Geen akoestische gegevens gelezen; is het Mecway-"
               "analysetype A30 (Acoustic Resonance 3D)?",
               "ERROR: No acoustic data read; is the Mecway analysis type "
               "A30 (Acoustic Resonance 3D)?"))
        return None
    if ac["bad"]:
        log(tr("ERROR: Akoestiek: {0}.", "ERROR: Acoustics: {0}.").format(
            "; ".join(ac["bad"])))
        return None
    if model_data.get("analysis_type") == "A20":
        log(tr("ERROR: Acoustic Resonance 2D (A20) wordt niet vertaald; de "
               "bridge rekent alleen 3D (A30).",
               "ERROR: Acoustic Resonance 2D (A20) is not translated; the "
               "bridge only computes 3D (A30)."))
        return None
    elsets = _elsets(model_data)
    if elsets is None:
        return None
    mats = _materials(elsets, ac)
    if mats is None:
        return None
    freqs = _frequencies(analysis)
    if freqs is None:
        return None
    inp = _harm_inputs(model_data, analysis, ac, mats)
    if inp is None:
        return None
    no_rho = [m for m in mats if not ac["materials"].get(m, {}).get("rho")]
    if inp["vite"] and no_rho:
        # De bron is i w rho v (VITE_FACE): zonder dichtheid zou p = c v
        # i.p.v. rho c v zijn (RHO = 1).
        log(tr("ERROR: Akoestiek harmonisch: materiaal {0} heeft geen "
               "dichtheid; een snelheidsbron geeft p = rho c v. Vul in "
               "Mecway bij het materiaal een dichtheid in (bijv. 1.2 kg/m3 "
               "voor lucht).",
               "ERROR: Acoustics harmonic: material {0} has no density; a "
               "velocity source gives p = rho c v. Fill in a density for "
               "the material in Mecway (e.g. 1.2 kg/m3 for air).").format(
                   ", ".join("'%s'" % m for m in no_rho)))
        return None
    _wavelength_check(model_data, elsets, mats, freqs[-1])
    log("INFO: Akoestiek harmonisch: {0} frequentie(s) van {1} tot {2} "
        "Hz.".format(len(freqs), "%.6g" % freqs[0], "%.6g" % freqs[-1]))

    # --- model: volumes en bronvlakken 3D, absorberende vlakken 3D_ABSO ---
    vol = [es for es, _m in elsets]
    src = list(item[0] for item in inp["vite"] if item[0] not in vol)
    src = [g for k, g in enumerate(src) if g not in src[:k]]
    abso = [a[0] for a in inp["abso"]]
    model_block = ("model = AFFE_MODELE(\n    MAILLAGE=mesh,\n    AFFE=(\n"
                   "        _F(GROUP_MA={0}, PHENOMENE='ACOUSTIQUE', "
                   "MODELISATION='3D'),\n").format(_fmt_groups(vol + src))
    if abso:
        model_block += ("        _F(GROUP_MA={0}, PHENOMENE='ACOUSTIQUE', "
                        "MODELISATION='3D_ABSO'),\n").format(
                            _fmt_groups(abso))
    model_block += "    ),\n)\n"

    # --- materialen ---
    names = OrderedDict()
    mat_lines = []
    for k, (mat, (c, rho)) in enumerate(mats.items()):
        names[mat] = "MAC{0}".format(k + 1)
        mat_lines.append("{0}{1} = DEFI_MATERIAU(FLUIDE=_F(RHO={2}, "
                         "CELE_R={3}))".format(bta._mat_comment(mat),
                                               names[mat], _r(rho), _r(c)))
    assign = ["        _F(GROUP_MA='{0}', MATER={1}),\n".format(es, names[m])
              for es, m in elsets]
    assign += ["        _F(GROUP_MA='{0}', MATER={1}),\n".format(
        g, names[inp["source_media"][g]]) for g in src]
    for k, (g, c, rho, desc) in enumerate(inp["abso"]):
        mat = inp["abso_media"][g]
        if abs(c - mats[mat][0]) <= 1.0e-12 * mats[mat][0]:
            mname = names[mat]
        else:
            # impedantie: eigen materiaal met c = Z / rho (zie hierboven)
            mname = "MAZ{0}".format(k + 1)
            mat_lines.append("# absorberende rand '{0}': {1}\n{2} = "
                             "DEFI_MATERIAU(FLUIDE=_F(RHO={3}, CELE_R={4}))"
                             .format(g, desc, mname, _r(rho), _r(c)))
        assign.append("        _F(GROUP_MA='{0}', MATER={1}),\n".format(
            g, mname))
    fieldmat_block = ("fieldmat = AFFE_MATERIAU(\n    MAILLAGE=mesh,\n"
                      "    MODELE=model,\n    AFFE=(\n{0}    ),\n)\n").format(
                          "".join(assign))

    # --- lasten ---
    loads = []
    if inp["vite"]:
        loads.append("# Bron: normaalsnelheid de fluid in (VITE_FACE VNOR)\n"
                     "ch_vite = AFFE_CHAR_ACOU(\n    MODELE=model,\n"
                     "    VITE_FACE=(\n{0}    ),\n)\n".format("".join(
                         "        _F(GROUP_MA='{0}', VNOR={1}),  # {2}\n"
                         .format(g, _ri(v), d) for g, v, d in inp["vite"])))
    if inp["pres"]:
        loads.append("# Opgelegde druk (fluid pressure; 0 = drukknoop)\n"
                     "bc_pres = AFFE_CHAR_ACOU(\n    MODELE=model,\n"
                     "    PRES_IMPO=(\n{0}    ),\n)\n".format("".join(
                         "        _F({0}='{1}', PRES={2}),\n".format(
                             kw, g, _ri(p)) for kw, g, p in inp["pres"])))
    loads_block = "\n".join(loads) + bta._custom_code(model_data, "code")
    if model_data.get("aster_custom", {}).get("excit"):
        # DYNA_VIBRA HARM neemt hier alleen VECT_ASSE (catalogus
        # dyna_vibra.py: CHARGE alleen char_meca); eigen akoestische lasten
        # kan de bridge niet in EXCIT zetten.
        log(tr("WARNING: Akoestiek harmonisch: #ASTER-EXCIT wordt niet "
               "gebruikt ({0}); gebruik #ASTER-BRON of een eigen DYNA_VIBRA "
               "in #ASTER-POST.",
               "WARNING: Acoustics harmonic: #ASTER-EXCIT is not used "
               "({0}); use #ASTER-BRON or your own DYNA_VIBRA in "
               "#ASTER-POST.").format(", ".join(
                   model_data["aster_custom"]["excit"])))
    p_source = any(abs(p) > 0.0 for _kw, _g, p in inp["pres"])

    # --- matrices en vectoren (ahlv100a) ---
    vect = []
    if inp["vite"]:
        vect.append("        _F(VECTEUR=CO('vec_v'), OPTION='CHAR_ACOU', "
                    "CHARGE=ch_vite),\n")
    if p_source:
        vect.append("        _F(VECTEUR=CO('vec_p'), OPTION='CHAR_ACOU', "
                    "CHARGE=bc_pres),\n")
    assemblage_block = (
        "ASSEMBLAGE(\n"
        "    MODELE=model,\n"
        "    CHAM_MATER=fieldmat,\n"
        "{0}"
        "    NUME_DDL=CO('numedd'),\n"
        "    MATR_ASSE=(\n"
        "        _F(MATRICE=CO('rigi_c'), OPTION='RIGI_ACOU'),\n"
        "        _F(MATRICE=CO('mass_c'), OPTION='MASS_ACOU'),\n"
        "{1}"
        "    ),\n"
        "    VECT_ASSE=(\n{2}    ),\n"
        ")\n").format(
            "    CHARGE=(bc_pres,),\n" if inp["pres"] else "",
            "        _F(MATRICE=CO('amor_c'), OPTION='AMOR_ACOU'),\n"
            if abso else "", "".join(vect))

    # --- harmonische berekening ---
    excit = []
    if inp["vite"]:
        # i w x vec_v (zie de kop van dit blok: v de fluid in)
        excit.append("        _F(VECT_ASSE=vec_v, COEF_MULT=1.0, "
                     "PUIS_PULS=1, PHAS_DEG=90.0),\n")
    if p_source:
        excit.append("        _F(VECT_ASSE=vec_p, COEF_MULT=1.0),\n")
    harm_block = (
        "freqs = DEFI_LIST_REEL(VALE=({0},))\n"
        "res = DYNA_VIBRA(\n"
        "    TYPE_CALCUL='HARM',\n"
        "    BASE_CALCUL='PHYS',\n"
        "    CHAM_MATER=fieldmat,\n"
        "    MATR_MASS=mass_c,\n"
        "    MATR_RIGI=rigi_c,\n"
        "{1}"
        "    LIST_FREQ=freqs,\n"
        "    EXCIT=(\n{2}    ),\n"
        "    SOLVEUR=_F(METHODE='MUMPS', NPREC={3}),\n"
        ")\n").format(", ".join(_r(f) for f in freqs),
                      "    MATR_AMOR=amor_c,\n" if abso else "",
                      "".join(excit), int(analysis.get("solver_nprec", 8)))

    export_block = (bta._custom_code(model_data, "post") +
                    _harm_export_block(work_dir, model_data, inp["mics"]))
    if CONFIG_BRIDGE.get("med_export", False):
        impr_resu_block = (
            "IMPR_RESU(FORMAT='MED', UNITE=80, RESU=(\n"
            "    _F(RESULTAT=res, NOM_CHAM='PRES', PARTIE='REEL', "
            "NOM_CHAM_MED='PRES_REEL'),\n"
            "    _F(RESULTAT=res, NOM_CHAM='PRES', PARTIE='IMAG', "
            "NOM_CHAM_MED='PRES_IMAG'),\n))\n")
    else:
        impr_resu_block = "#-- geen MED-export (instelling med_export uit) --\n"

    comm_path = os.path.join(work_dir, run_name(".comm"))
    try:
        template = bta._read_template(analysis["template"])
        bta._write_comm(comm_path, template.format(
            model_block=model_block,
            materials_block="\n".join(mat_lines) + "\n",
            fieldmat_block=fieldmat_block,
            loads_block=loads_block,
            assemblage_block=assemblage_block,
            harm_block=harm_block,
            export_block=export_block,
            impr_resu_block=impr_resu_block,
        ))
        log("SUCCESS: .comm geschreven naar '{0}'.".format(comm_path))
        export_path = os.path.join(work_dir, run_name(".export"))
        with open(export_path, "w") as f:
            f.write("\n".join([
                "P actions make_etude",
                "P time_limit {0}".format(analysis.get("time_limit", 900)),
                "P memory_limit {0}".format(analysis.get("memory_limit", 2048)),
                "P mpi_nbcpu 1",
                "P ncpus {0}".format(CONFIG_BRIDGE["ncpus"]),
                "P no-mpi true",
                "F comm {0} D 1".format(comm_path),
                "F mail {0} D 20".format(msh_path),
                "F mess {0} R 6".format(os.path.join(work_dir,
                                                     run_name(".mess"))),
            ]) + bta._export_result_lines(work_dir, resu_name))
        log("SUCCESS: .export geschreven naar '{0}'.".format(export_path))
        return export_path
    except KeyError as e:
        log(tr("ERROR: Template placeholder ontbreekt: {0}.",
               "ERROR: Template placeholder missing: {0}.").format(str(e)))
    except Exception as e:
        log(tr("ERROR: Genereren akoestische solver-bestanden mislukt: {0}",
               "ERROR: Generating the acoustic solver files failed: {0}")
            .format(str(e)))
    return None


MIC_SUFFIX = "_microfoons.csv"


def _harm_export_block(work_dir, model_data, mics):
    """aster_acoustic_export: een <mode> per frequentie (pressure,
    pressurei, pressuremag, pressurephase, spl) en de microfoons."""
    import bridge_to_aster as bta
    scripts = os.path.dirname(os.path.abspath(__file__))
    mic_text = "{" + ", ".join("{0}: [{1}]".format(
        bta._py_str(name), ", ".join(str(n) for n in nids))
        for name, nids in mics.items()) + "}"
    return ("# Resultaten in Mecway's eigen formaat (result_route 'liml')\n"
            "import sys\n"
            "sys.path.insert(0, {0})\n"
            "import aster_acoustic_export\n"
            "aster_acoustic_export.export_harmonic(\n"
            "    res, {1}, {2},\n"
            "    microphones={3},\n"
            "    csv_path={4},\n"
            "    language={5})\n").format(
                bta._py_str(scripts), bta._py_str(bta.results_file(work_dir)),
                bta._py_str(model_data.get("liml_path", "")), mic_text,
                bta._py_str(os.path.join(work_dir, run_name(MIC_SUFFIX))),
                bta._py_str(CONFIG_BRIDGE.get("language", "en")))


# ==========================================
# NA DE RUN: MICROFOONS IN DE LOG (IronPython)
# ==========================================
def read_microphones(csv_path):
    """[(naam, f, amplitude, fase, spl)] uit de CSV van
    aster_acoustic_export; [] als hij ontbreekt."""
    import csv
    out = []
    if not csv_path or not os.path.isfile(csv_path):
        return out
    f = open(csv_path, "r")
    try:
        rows = list(csv.reader(f))
    finally:
        f.close()
    for row in rows[1:]:
        try:
            out.append((row[0], float(row[1]), float(row[4]),
                        float(row[5]), float(row[6])))
        except (IndexError, ValueError):
            continue
    return out


def report_harmonic(work_dir, frequencies):
    """
    Na een harmonische run: aantal frequenties en per microfoon het
    hoogste geluidsdrukniveau in de log; de regels voor het statusvak.
    frequencies: [(id, "frequency", f)] uit read_modes().
    """
    path = os.path.join(work_dir, run_name(MIC_SUFFIX))
    lines = [tr(u"Akoestiek harmonisch: {0} frequentie(s), {1} .. {2} Hz. "
                u"Velden in Mecway per frequentie: pressure (Re), pressurei "
                u"(Im), pressuremag (|p|), pressurephase (graden), spl (dB "
                u"re 20 uPa).",
                u"Acoustics harmonic: {0} frequency(ies), {1} .. {2} Hz. "
                u"Fields in Mecway per frequency: pressure (Re), pressurei "
                u"(Im), pressuremag (|p|), pressurephase (degrees), spl (dB "
                u"re 20 uPa).").format(
                    len(frequencies),
                    "%.6g" % min(f for _i, _k, f in frequencies)
                    if frequencies else "-",
                    "%.6g" % max(f for _i, _k, f in frequencies)
                    if frequencies else "-")]
    mics = read_microphones(path)
    by_name = OrderedDict()
    for name, f, amp, phase, level in mics:
        by_name.setdefault(name, []).append((f, amp, phase, level))
    for name, rows in by_name.items():
        f, amp, _ph, level = max(rows, key=lambda r: r[3])
        lines.append(tr(u"Microfoon {0}: hoogste niveau {1} dB bij {2} Hz "
                        u"(|p| = {3} Pa).",
                        u"Microphone {0}: highest level {1} dB at {2} Hz "
                        u"(|p| = {3} Pa).").format(
                            name, "%.1f" % level, "%.6g" % f, "%.4g" % amp))
        for f, amp, phase, level in rows:
            log(u"INFO: Microfoon {0}: f = {1} Hz, |p| = {2} Pa, fase = {3} "
                u"graden, SPL = {4} dB.".format(
                    name, "%.6g" % f, "%.6g" % amp, "%.2f" % phase,
                    "%.2f" % level))
    if by_name:
        lines.append(tr(u"Frequentierespons per microfoon: {0}",
                        u"Frequency response per microphone: {0}").format(
                            path))
    for line in lines:
        log(u"INFO: " + line)
    return lines
