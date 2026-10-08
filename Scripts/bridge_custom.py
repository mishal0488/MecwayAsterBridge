# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# bridge_custom.py
# Verantwoordelijkheid: #ASTER-regels lezen die de gebruiker in Mecway in een
# vrij tekstveld zet, en ze klaarzetten voor bridge_to_aster.py. Zo kan een
# model Code-Aster-functies gebruiken waar Mecway geen invoer voor heeft
# (bijv. LIAISON_ELEM voor balk-volume-verbindingen) zonder de bridge aan te
# passen.
#
# Waar: de vrije tekstvelden van Mecway (XML-tags, gevonden in mecway.exe):
#   <ccxmodeldefinition>      CCX > Custom model definition
#   <ccxstep>                 CCX > Custom step contents
#   <radiossmodeldefinition>  OpenRadioss > Custom model definition
#   <radiossengineinput>      OpenRadioss > Custom engine input
#   <nastranbulkdata>         Mystran > Custom bulk data
# Een regel begint met '#' (commentaar voor OpenRadioss), '**' (CCX) of '$'
# (Mystran/Nastran), daarna ASTER; de Mecway-solver negeert de regel dus.
# Praktisch: dynamica/contact in een OpenRadioss-analyse, statica in een
# Mystran-analyse; de bridge leest de velden ongeacht de gekozen solver.
#
# Regels (voorvoegsel '#', '**' of '$', spaties ervoor mogen):
#   #ASTER <code>            Code-Aster-code, letterlijk vóór de solve
#                            (na de lasten van de bridge). Opeenvolgende
#                            regels vormen een blok; de inspringing na
#                            '#ASTER ' blijft behouden.
#   #ASTER-EXCIT naam, ...   last(en) uit #ASTER-code toevoegen aan EXCIT
#   #ASTER-CONTACT K=V, ...  extra sleutelwoorden voor DEFI_CONTACT (globaal)
#   #ASTER-ZONE K=V, ...     extra sleutelwoorden voor elke contactzone
#   #ASTER-OPTIE sleutel = waarde
#                            instelling van de actieve analyse voor dit
#                            model (bijv. newton_max_iter = 50)
#   #ASTER-POST <code>       Code-Aster-code na de solve (vóór IMPR_RESU)
#   #ASTER-BEZWIJKEN materiaal = rek
#                            elementen van dit Mecway-materiaal bezwijken
#                            bij deze vergelijkingsplastische rek ('*' =
#                            alle plastische materialen); ook
#                            #ASTER-FAILURE. Zie bridge_failure.py.
#   #ASTER-SPECTRE X, Y = naam [, K=V, ...]
#                            responsspectrum (analyse spectrum_3d) voor
#                            deze richtingen: een DEFI_NAPPE/DEFI_FONCTION/
#                            FORMULE uit #ASTER-code (SPEC_OSCI), met
#                            eventueel extra sleutelwoorden voor SPECTRE
#                            (ECHELLE, NATURE, CORR_FREQ); ook
#                            #ASTER-SPECTRUM. Zie bridge_spectrum.py.
#   #ASTER-BRON selectie = snelheid [, fase]
#                            akoestiek harmonisch (acoustic_harm_3d): een
#                            vlakselectie als bron met normaalsnelheid de
#                            fluid in (m/s of mm/s) en fase in graden; ook
#                            #ASTER-SOURCE. Zie bridge_acoustic.py. In de
#                            vibro-akoestiek harmonisch (vibroacoustic_
#                            harm_3d) ook met een druk (Pa): een invallende
#                            golf. Zie bridge_vibroacoustic_harm.py.
#   #ASTER-LASBRON K=V, ...  bewegende lasbron (transient thermisch; ook
#                            #ASTER-WELD). Zie bridge_welding.py.
#   #ASTER-FASE selectie = graden
#                            harmonische respons (harmonic_3d): de
#                            belasting op deze selectie krijgt deze fase
#                            (EXCIT PHAS_DEG); ook #ASTER-PHASE. Zie
#                            bridge_harmonic.py.
#   #ASTER-MATERIAAL naam: K=V, ...
#                            temperatuurafhankelijke E, nu, alpha, SY, ET
#                            van een Mecway-materiaal (keten thermisch ->
#                            mechanisch; ook #ASTER-MATERIAL). Zie
#                            bridge_thermomech.py.
#   #ASTER-WOHLER materiaal: K=V, ...
#                            Wohler-curve voor vermoeiing; ook
#                            #ASTER-VERMOEIING/-FATIGUE K=V (opties),
#                            #ASTER-BLOK/-BLOCK a, m, n en #ASTER-SIGNAAL/
#                            -SIGNAL f1, f2, ... Zie bridge_fatigue.py.
#   #ASTER-METALLURGIE naam: K=V, ...
#                            fase-omzettingen in staal (TRC, CALC_META) voor
#                            een Mecway-materiaal (transient thermisch; ook
#                            #ASTER-METALLURGY). Zie bridge_metallurgy.py.
#   #ASTER-XFEM K=V, ...     XFEM-scheur zonder hem te meshen (lineair
#                            statisch: K en G langs het front). Zie
#                            bridge_xfem.py.
#   <<Mecway-selectie>>      in alle regels: de Code-Aster-groepsnaam van
#                            die selectie (bridge_mesh hernoemt namen)
#
# Voorbeeld (balkknoop aan een volumevlak, LIAISON_ELEM OPTION='3D_POU'):
#   #ASTER verb = AFFE_CHAR_MECA(MODELE=model,
#   #ASTER     LIAISON_ELEM=_F(OPTION='3D_POU', GROUP_MA_1='<<Vlak>>',
#   #ASTER                     GROUP_NO_2='<<Knoop>>'))
#   #ASTER-EXCIT verb
#
# De bridge controleert de Code-Aster-syntax van #ASTER-code niet; fouten
# komen uit Code-Aster (study.mess).
#
# Compatibel met IronPython 2.7.9.

import re

from bridge_utils import log, tr, parse_yes_no

CUSTOM_TAGS = ("ccxmodeldefinition", "ccxstep", "radiossmodeldefinition",
               "radiossengineinput", "nastranbulkdata")

# Soort ook in kleine letters herkennen, zodat '#ASTER-Excit' niet stil
# wegvalt maar als soort (of als onbekende soort) wordt gemeld.
_LINE_RE = re.compile(r"^\s*(?:#|\*\*|\$)\s*ASTER(-[A-Za-z]+)?(?:\s(.*))?$")
# Lijkt op een #ASTER-regel (voorvoegsel direct gevolgd door ASTER) maar
# past niet in _LINE_RE, bijv. '#ASTER_EXCIT' of '#ASTER:': melden i.p.v.
# stil overslaan. Zonder spatie na het voorvoegsel, zodat gewone CCX-
# commentaarregels ('** ...') niet geraakt worden.
_NEAR_RE = re.compile(r"^\s*(?:#|\*\*|\$)ASTER")
_PLACEHOLDER_RE = re.compile(r"<<([^<>]+)>>")
_KINDS = ("", "-EXCIT", "-CONTACT", "-ZONE", "-OPTIE", "-POST",
          "-BEZWIJKEN", "-FAILURE", "-SPECTRE", "-SPECTRUM",
          "-BRON", "-SOURCE", "-LASBRON", "-WELD",
          "-FASE", "-PHASE", "-MATERIAAL", "-MATERIAL",
          "-WOHLER", "-VERMOEIING", "-FATIGUE", "-BLOK", "-BLOCK",
          "-SIGNAAL", "-SIGNAL", "-METALLURGIE", "-METALLURGY",
          "-XFEM")
# Vermoeiing (bridge_fatigue.py): soort -> sleutel in out["fatigue"]
_FATIGUE_KINDS = {"-WOHLER": "wohler", "-VERMOEIING": "options",
                  "-FATIGUE": "options", "-BLOK": "block", "-BLOCK": "block",
                  "-SIGNAAL": "signal", "-SIGNAL": "signal"}


def empty():
    """Lege structuur zoals parse() hem teruggeeft."""
    return {"code": [], "excit": [], "contact": [], "zone": [],
            "options": {}, "post": [], "failure": {}, "spectrum": {},
            "sources": [], "welding": [],
            "phases": [], "materials": [],
            "fatigue": [], "metallurgy": [],
            "xfem": []}


def parse(texts):
    """
    #ASTER-regels uit de vrije tekstvelden.

    texts : lijst van strings (inhoud van de CUSTOM_TAGS)
    Returns dict met "code"/"post" (regels), "excit" (namen),
    "contact"/"zone" (stukken 'K=V'), "options" ({sleutel: tekst}).
    """
    out = empty()
    for text in texts:
        for raw in (text or "").splitlines():
            m = _LINE_RE.match(raw.rstrip("\r"))
            if not m:
                if _NEAR_RE.match(raw):
                    log(tr("ERROR: Regel lijkt op een #ASTER-regel maar wordt "
                           "niet herkend en wordt overgeslagen (verwacht "
                           "'#ASTER <code>' of '#ASTER-SOORT ...'): {0}",
                           "ERROR: Line looks like an #ASTER line but is "
                           "not recognised and is skipped (expected "
                           "'#ASTER <code>' or '#ASTER-KIND ...'): {0}")
                        .format(raw.strip()))
                continue
            kind = (m.group(1) or "").upper()
            body = m.group(2) or ""
            if kind not in _KINDS:
                log(tr("ERROR: Onbekende regel '#ASTER{0}' wordt "
                       "overgeslagen: {1}",
                       "ERROR: Unknown line '#ASTER{0}' is skipped: "
                       "{1}").format(kind, raw.strip()))
                continue
            if kind == "":
                out["code"].append(body)
            elif kind == "-POST":
                out["post"].append(body)
            elif kind == "-EXCIT":
                out["excit"] += [n.strip() for n in body.split(",")
                                 if n.strip()]
            elif kind == "-CONTACT":
                if body.strip():
                    out["contact"].append(body.strip().rstrip(","))
            elif kind == "-ZONE":
                if body.strip():
                    out["zone"].append(body.strip().rstrip(","))
            elif kind in ("-BEZWIJKEN", "-FAILURE"):
                _parse_failure(body, out["failure"])
            elif kind in ("-SPECTRE", "-SPECTRUM"):
                _parse_spectrum(body, out["spectrum"])
            elif kind in ("-BRON", "-SOURCE"):
                # uitgewerkt in bridge_acoustic._parse_aster_sources en
                # bridge_vibroacoustic_harm._aster_sources
                if body.strip():
                    out["sources"].append(body.strip())
            elif kind in ("-LASBRON", "-WELD"):
                if body.strip():
                    out["welding"].append(body.strip())
            elif kind in _FATIGUE_KINDS:
                # uitgewerkt in bridge_fatigue.parse
                if body.strip():
                    out["fatigue"].append((_FATIGUE_KINDS[kind],
                                           body.strip()))
            elif kind in ("-FASE", "-PHASE"):
                # uitgewerkt in bridge_harmonic._phase_bodies
                if body.strip():
                    out["phases"].append(body.strip())
            elif kind == "-XFEM":
                # uitgewerkt in bridge_xfem.parse_specs
                if body.strip():
                    out["xfem"].append(body.strip())
            elif kind in ("-MATERIAAL", "-MATERIAL"):
                # uitgewerkt in bridge_thermomech.parse_materials
                if body.strip():
                    out["materials"].append(body.strip())
            elif kind in ("-METALLURGIE", "-METALLURGY"):
                # uitgewerkt in bridge_metallurgy.parse_specs
                if body.strip():
                    out["metallurgy"].append(body.strip())
            else:
                key, sep, value = body.partition("=")
                if not sep or not key.strip():
                    log(tr("ERROR: '#ASTER-OPTIE {0}': verwacht 'sleutel = "
                           "waarde'.",
                           "ERROR: '#ASTER-OPTIE {0}': expected 'key = "
                           "value'.").format(body.strip()))
                    continue
                out["options"][key.strip()] = value.strip()
    n = sum(len(v) for v in out.values())
    if n:
        log("INFO: #ASTER-regels: {0} code, {1} EXCIT, {2} CONTACT, {3} ZONE, "
            "{4} OPTIE, {5} POST, {6} BEZWIJKEN, {7} SPECTRE, {8} BRON, "
            "{9} LASBRON, {10} FASE, {11} MATERIAAL, {12} VERMOEIING, "
            "{13} METALLURGIE, {14} XFEM.".format(
                len(out["code"]), len(out["excit"]), len(out["contact"]),
                len(out["zone"]), len(out["options"]), len(out["post"]),
                len(out["failure"]), len(out["spectrum"]),
                len(out["sources"]), len(out["welding"]),
                len(out["phases"]), len(out["materials"]),
                len(out["fatigue"]), len(out["metallurgy"]),
                len(out["xfem"])))
    return out


_SPECTRUM_RE = re.compile(r"^\s*([XYZxyz](?:\s*,\s*[XYZxyz])*)\s*=\s*"
                          r"([A-Za-z_]\w*)\s*(?:,(.*))?$")


def _parse_spectrum(body, spectrum):
    """'#ASTER-SPECTRE X, Y = naam, ECHELLE=9.81' -> spectrum[richting] =
    (naam, ['ECHELLE=9.81']) (bridge_spectrum.py: SPECTRE van
    COMB_SISM_MODAL)."""
    m = _SPECTRUM_RE.match(body)
    extra = []
    if m and m.group(3):
        extra = [p.strip() for p in m.group(3).split(",") if p.strip()]
    if not m or any("=" not in p for p in extra):
        log(tr("ERROR: '#ASTER-SPECTRE {0}': verwacht 'X, Y = naam' met "
               "eventueel ', K=V' (bijv. '#ASTER-SPECTRE X = spec, "
               "ECHELLE=9.81'); genegeerd.",
               "ERROR: '#ASTER-SPECTRE {0}': expected 'X, Y = name' with "
               "optionally ', K=V' (e.g. '#ASTER-SPECTRE X = spec, "
               "ECHELLE=9.81'); ignored.").format(body.strip()))
        return
    for axis in re.findall(r"[XYZxyz]", m.group(1)):
        spectrum[axis.upper()] = (m.group(2), extra)


def _parse_failure(body, failure):
    """'#ASTER-BEZWIJKEN materiaal = rek' -> failure[materiaal] = rek."""
    name, sep, value = body.rpartition("=")
    try:
        limit = float(value.strip().replace(",", "."))
    except ValueError:
        limit = None
    if not sep or not name.strip() or limit is None or limit <= 0.0:
        log(tr("ERROR: '#ASTER-BEZWIJKEN {0}': verwacht 'materiaal = rek' "
               "met een rek > 0 (bijv. '#ASTER-BEZWIJKEN Staal = 0.2'); "
               "genegeerd.",
               "ERROR: '#ASTER-BEZWIJKEN {0}': expected 'material = "
               "strain' with a strain > 0 (e.g. '#ASTER-BEZWIJKEN Steel = "
               "0.2'); ignored.").format(body.strip()))
        return
    failure[name.strip()] = limit


def substitute(text, renames, known):
    """
    Vervang <<Mecway-selectie>> door de Code-Aster-groepsnaam.

    renames : {Mecway-naam: Code-Aster-naam} (bridge_mesh)
    known   : verzameling Code-Aster-groepsnamen in het mesh
    Een onbekende selectie wordt gelogd als ERROR en blijft staan, zodat
    Code-Aster er een duidelijke fout op geeft.
    """
    def repl(m):
        name = m.group(1).strip()
        aster = renames.get(name, name)
        if aster not in known:
            log(tr("ERROR: #ASTER: selectie '{0}' bestaat niet in het "
                   "mesh.",
                   "ERROR: #ASTER: selection '{0}' does not exist in the "
                   "mesh.").format(name))
        return aster
    return _PLACEHOLDER_RE.sub(repl, text)


def apply_options(analysis, options):
    """
    Kopie van de analyse-instellingen met de #ASTER-OPTIE-waarden. Alleen
    sleutels die de analyse al heeft; het type volgt de bestaande waarde.
    Een ja/nee-instelling kent de woorden van bridge_utils.parse_yes_no
    (ja/nee, aan/uit, yes/no, on/off, true/false, 1/0, oui/non); een andere
    tekst geeft een WARNING en de bestaande waarde blijft.
    """
    result = dict(analysis)
    for key, raw in options.items():
        if key not in analysis:
            log(tr("ERROR: #ASTER-OPTIE '{0}' bestaat niet voor deze "
                   "analyse; genegeerd. Bekend: {1}",
                   "ERROR: #ASTER-OPTIE '{0}' does not exist for this "
                   "analysis; ignored. Known: {1}").format(
                    key, ", ".join(sorted(analysis))))
            continue
        old = analysis[key]
        try:
            if isinstance(old, bool):
                value = parse_yes_no(raw)
                if value is None:
                    log(tr("WARNING: #ASTER-OPTIE {0} = {1}: geen ja/nee-"
                           "waarde (ja/nee, aan/uit, yes/no, on/off, "
                           "true/false, 1/0, oui/non); de instelling blijft "
                           "{2}.",
                           "WARNING: #ASTER-OPTIE {0} = {1}: not a yes/no "
                           "value (yes/no, on/off, true/false, 1/0, ja/nee, "
                           "aan/uit, oui/non); the setting stays "
                           "{2}.").format(key, raw, old))
                    continue
            elif isinstance(old, int):
                value = int(raw)
            elif isinstance(old, float):
                value = float(raw)
            else:
                value = raw.strip("'\"")
        except ValueError:
            log(tr("ERROR: #ASTER-OPTIE {0} = {1}: geen geldige waarde; "
                   "genegeerd.",
                   "ERROR: #ASTER-OPTIE {0} = {1}: not a valid value; "
                   "ignored.").format(key, raw))
            continue
        result[key] = value
        log("INFO: #ASTER-OPTIE {0} = {1} (was {2}).".format(key, value, old))
    return result
