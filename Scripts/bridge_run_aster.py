# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# bridge_run_aster.py
# Verantwoordelijkheid: start Code-Aster als extern proces via de
# run_aster.bat launcher, en zoekt na de run fouten in study.mess
# (aster_mess_errors). Geisoleerd van bestandsgeneratie (bridge_to_aster.py)
# en resultaat-parsing (bridge_from_aster.py) omdat falen hier een andere
# oorzaak heeft (launcher niet gevonden, permissie-probleem) dan falen
# in die andere modules.
#
# Synchronisatie: Code-Aster draait in een apart cmd-venster. Deze module
# wacht NIET tot de solver klaar is; in plaats daarvan toont mw.message()
# een dialoog die de gebruiker handmatig moet bevestigen. Dat is een
# bewuste keuze omdat run_aster.bat geen betrouwbaar exit-signaal geeft
# en de solver-output (.mess) asynchroon wordt weggeschreven.
#
# Console-uitvoer: de uitvoer van run_aster (ook stderr) gaat via "tee" uit
# de Code-Aster-installatie (tools\tee.exe) naar het venster EN naar
# study_console.txt. Een Python-fout in study.comm staat alleen daar (niet
# in study.mess; KNOWN_BUGS 13 opgelost).
#
# Convergentiescherm: bij elke run (instelling convergence_monitor) start
# start_convergence_monitor() aster_convergence_monitor.py als apart proces
# (de CPython van Code-Aster), zie de kop van dat bestand; achteraf via de
# knop Convergentiescherm in het hoofdvenster (bridge_gui).
#
# Compatibel met IronPython 2.7.9.

import os
import re
import subprocess
import time

from bridge_utils import (log, CONFIG_BRIDGE, tr, run_name, run_names,
                          exe_aster_path)
from bridge_convergence import python_errors


# Meldingen in study.mess staan als "<F> <CODE>" in een kader; <A> is een
# alarm. Een <EXCEPTION> stopt de commando's erna, maar de launcher meldt
# dan soms alleen "DIAGNOSTIC JOB : <A>_ALARM" (gezien met Code-Aster 17.4
# op Windows). Daarom wordt het .mess-bestand zelf doorzocht.
# Ook <ConvergenceError> en andere <...Error>: STAT_NON_LINE stopte met
# <ConvergenceError> <MECANONLINE9_7> terwijl de DIAGNOSTIC-regel "OK" gaf
# (gezien in Mecway).
_MESS_ERROR_RE = re.compile(
    r"<(EXCEPTION|F|S|E|[A-Za-z]+Error)>\s*<([A-Za-z0-9_]+)>")
# Zonder meldingskader, na een Python-exceptie in een commando.
_MESS_UNCAUGHT = "after an uncaught exception"
_MESS_ALARM_RE = re.compile(r"<A>\s*<([A-Za-z0-9_]+)>")
# Uitleg bij alarmen van CALC_MODES (bridge_to_aster._calc_modes_block:
# VERI_MODE STOP_ERREUR='NON', dus alarmen i.p.v. fouten). Tekst van de
# meldingen uit Code-Aster 17.4 (gezien met Mecway's BucklingBeam).
_ALARM_HINTS = {
    "ALGELINE5_23": (
        u"de Sturm-controle telt in het interval van de berekende modes "
        u"meer eigenwaarden dan berekend. Meestal valt een meervoudige mode "
        u"(symmetrie) op de grens van het aantal modes; anders is er een "
        u"mode gemist. Vraag een of twee modes meer en vergelijk.",
        u"the Sturm check counts more eigenvalues in the interval of the "
        u"computed modes than were computed. Usually a multiple mode "
        u"(symmetry) falls on the limit of the number of modes; otherwise a "
        u"mode was missed. Ask for one or two more modes and compare."),
    # Gezien met een vrij-vrij model en een negatieve zoekfrequentie
    # (docs-review): de ondergrens van de Sturm-band valt op een starre mode.
    "ALGELINE3_82": (
        u"de ondergrens van het zoekinterval valt (bijna) op een eigenwaarde, "
        u"meestal een starre-lichaamsmode bij ~0 Hz (akoestiek: de constante "
        u"druk bij harde wanden). Meestal onschuldig: de "
        u"Sturm-telling kan daardoor afwijken. Controleer het aantal modes "
        u"of verschuif de zoekfrequentie iets.",
        u"the lower bound of the search interval (nearly) coincides with an "
        u"eigenvalue, usually a rigid body mode at ~0 Hz (acoustics: the "
        u"constant pressure with hard walls). Usually harmless: "
        u"the Sturm count may be off because of it. Check the number of "
        u"modes or shift the search frequency slightly."),
    # Dubbele cellen (zelfde knopen). Vibro-akoestiek: een schaal op een
    # fluidumvlak heeft een eigen koppelvlak-cel (FLUI_STRU) met dezelfde
    # knopen (bridge_vibroacoustic.interface_faces); gezien in verificatie
    # VA-03/VA-05, frequenties exact.
    "MODELISA4_8": (
        u"dubbele cellen (zelfde knopen). Bij vibro-akoestiek normaal voor "
        u"een schaal op een lucht- of vloeistofvlak (schaal en koppelvlak "
        u"zijn aparte cellen); anders: controleer het mesh op dubbele "
        u"elementen.",
        u"duplicate cells (same nodes). Normal in vibroacoustics for a shell "
        u"on an air or fluid face (shell and coupling face are separate "
        u"cells); otherwise: check the mesh for duplicate elements."),
    "MODELISA8_13": (
        u"zie MODELISA4_8 (dubbele cellen).",
        u"see MODELISA4_8 (duplicate cells)."),
    "ALGELINE2_74": (
        u"de controle van de modes (VERI_MODE: fout per mode en Sturm) is "
        u"niet gehaald; zie de tabel in study.mess.",
        u"the check of the modes (VERI_MODE: error per mode and Sturm) "
        u"failed; see the table in study.mess."),
    # Complexe modes (bridge_modal_damping): Code-Aster geeft alleen de
    # oscillerende modes (gezien: beta = 1e-3 s op een hex20-console gaf 3
    # van de 6 modes; vrij-vrij met alpha: de 6 starre modes vielen weg).
    "ALGELINE5_79": (
        u"complexe modes: minder modes dan gevraagd. Overkritisch gedempte "
        u"modes (zeta >= 1, vaak hoge modes met Rayleigh-beta of lage met "
        u"alpha) en starre modes trillen niet en worden niet gegeven. "
        u"Verlaag de demping als je die modes wilt zien; zeta per mode staat "
        u"in de log.",
        u"complex modes: fewer modes than requested. Overdamped modes (zeta "
        u">= 1, often high modes with Rayleigh beta or low ones with alpha) "
        u"and rigid body modes do not oscillate and are not given. Lower the "
        u"damping to see those modes; zeta per mode is in the log."),
    # Gezien met complexe modes in een vrij-vrij model (de stijfheid is
    # singulier); de elastische modes waren gelijk aan Mecway (MC-05).
    "FACTOR_11": (
        u"een (bijna) singuliere matrix bij het ontbinden. Bij complexe "
        u"modes in een model zonder opleggingen normaal (starre "
        u"bewegingen); anders: controleer de opleggingen.",
        u"a (nearly) singular matrix in the factorisation. Normal for "
        u"complex modes in a model without supports (rigid body motion); "
        u"otherwise: check the supports."),
}
_MESS_BOX_END = u"╚"        # onderrand van het meldingskader
_MESS_BOX_SIDE = u"║ "      # zijkant van het kader en spatie
# Laatste regel die run_aster in study.mess zet, ook na een <EXCEPTION>.
_MESS_DONE = "--- DIAGNOSTIC JOB :"
# Uitleg bij meldingen die meestal op het model wijzen (tekst uit
# Messages/factor.py van Code-Aster 17.4). Een vrij mechanisme geeft met
# NPREC=-1 (geen singulariteitsdetectie) FACTOR_57 i.p.v. FACTOR_10/11/42.
# Teksten als paar (nl, en); tr(*paar) kiest bij gebruik de taal.
_SINGULAR_HINT = (u"De stijfheidsmatrix is (bijna) singulier: waarschijnlijk "
                  u"kan een deel van het model vrij bewegen (te weinig "
                  u"opleggingen; een pinned support houdt alleen de "
                  u"gemiddelde verplaatsing vast) of zijn opleggingen dubbel. "
                  u"Zie docs/Foutmeldingen_CodeAster.md.",
                  u"The stiffness matrix is (nearly) singular: probably a "
                  u"part of the model can move freely (too few supports; a "
                  u"pinned support only fixes the average displacement) or "
                  u"supports are duplicated. See "
                  u"docs/en/CodeAster_error_messages.md.")
# Newton convergeert niet (STAT_NON_LINE/DYNA_NON_LINE). De bridge laat de
# job doorgaan met de stappen tot het laatste geconvergeerde tijdstip
# (BRIDGE_ONVOLLEDIG); advies volgens U2.04.04 par. 4.1 en 4.5.
_NEWTON_HINT = (u"Newton convergeert niet; de resultaten tot het laatst "
                u"bereikte tijdstip worden ingeladen. Controleer eerst het "
                u"model: "
                u"contactvlakken (master = het stijve gereedschap, slave = "
                u"het fijnere of vervormende vlak), beginspleten en "
                u"-penetraties, en of delen alleen via contact op hun plaats "
                u"blijven. Reken stap voor stap: eerst contact met lineair "
                u"materiaal, dan plasticiteit. Doorsnijden of bezwijken: geef "
                u"het materiaal een bezwijkgrens (docs/Bezwijken.md); na een "
                u"bezwijkmoment helpen kleinere stappen. Blijft RESI_GLOB_RELA "
                u"in de eerste stappen hoog terwijl RESI_GLOB_MAXI heel klein "
                u"is (een beweging zonder weerstand), kies dan "
                u"Convergentiecriterium = MAXI. Stagneert RESI_GLOB_RELA net "
                u"boven het criterium (bijv. 2e-6 bij 1e-6, gezien bij "
                u"ontlasten met GDEF_LOG), gebruik dan de standaard 1e-5, of "
                u"Convergentiecriterium = RELA_OF_MAXI met RESI_GLOB_MAXI "
                u"boven de absolute rest waarop het stagneert (kolom "
                u"RESI_GLOB_MAXI in de .mess). Bij contact kan "
                u"Contactmethode = PENALISATION (Mecway's stijfheid) helpen. "
                u"Zie "
                u"docs/Instellingen_contact_en_nietlineair.md en "
                u"docs/Foutmeldingen_CodeAster.md.",
                u"Newton does not converge; the results up to the last "
                u"time reached are loaded. Check the model first: "
                u"contact surfaces (master = the stiff tool, slave = the "
                u"finer or deforming surface), initial gaps and "
                u"penetrations, and whether parts are held in place only by "
                u"contact. Work step by step: first contact with linear "
                u"material, then plasticity. Cutting through or collapse: "
                u"give the material a failure limit (docs/en/Failure.md); "
                u"after a failure event smaller steps help. If "
                u"RESI_GLOB_RELA stays high in the first steps while "
                u"RESI_GLOB_MAXI is very small (a motion without "
                u"resistance), choose Convergence criterion = MAXI. If "
                u"RESI_GLOB_RELA stalls just above the criterion (e.g. 2e-6 "
                u"at 1e-6, seen when unloading with GDEF_LOG), use the "
                u"default 1e-5, or Convergence criterion = RELA_OF_MAXI with "
                u"RESI_GLOB_MAXI above the absolute residual where it stalls "
                u"(column RESI_GLOB_MAXI in the .mess). With "
                u"contact, Contact method = PENALISATION (Mecway's "
                u"stiffness) can help. See "
                u"docs/en/Settings_contact_and_nonlinear.md and "
                u"docs/en/CodeAster_error_messages.md.")
# MECANONLINE9_4 met contact: pivot op LAGS_C (contact-Lagrange) bij
# lijncontact (Tutorials/AssemblyWithContact Fixed_support, headless
# vergeleken); met PENALISATION rekent het door. NPREC -1 rekent ook door
# maar is onnauwkeurig (docs-review: -1,6% / +2,8% op een stijfheidscontrast),
# dus geen advies.
_NEWTON_SINGULAR_HINT = (
    _SINGULAR_HINT[0] + u" Met contact en 'LAGS_C' in de melding: een "
    u"contactvrijheidsgraad zonder steun (bijv. lijncontact); kies "
    u"Contactmethode = PENALISATION (Mecway's stijfheid).",
    _SINGULAR_HINT[1] + u" With contact and 'LAGS_C' in the message: a "
    u"contact degree of freedom without support (e.g. line contact); "
    u"choose Contact method = PENALISATION (Mecway's stiffness).")
# CALC_MODES (bridge_to_aster._calc_modes_block): de stijfheid is bij de
# zoekfrequentie singulier (ALGELINE3_81, gezien met een balk die alleen
# scharnierend was vastgezet), of Sorensen vond te weinig modes
# (ALGELINE5_49 / ALGELINE3_78, gezien met een zoekfrequentie ver van de
# modes).
_MODES_SINGULAR_HINT = (
    u"De stijfheidsmatrix is singulier: het model heeft een mechanisme of "
    u"starre-lichaamsbeweging (te weinig opleggingen, een los onderdeel, "
    u"een scharnier). Leg het model beter vast, of kies op de pagina "
    u"Modaal een negatieve zoekfrequentie (bijv. -10 x de verwachte eerste "
    u"frequentie); de starre modes komen dan met ~0 Hz mee. Dat helpt alleen "
    u"als het losse deel massa heeft (dichtheid of puntmassa). Een "
    u"verbinding uit #ASTER-code telt alleen mee via #ASTER-EXCIT.",
    u"The stiffness matrix is singular: the model has a mechanism or rigid "
    u"body motion (too few supports, a loose part, a hinge). Support the "
    u"model better, or choose a negative search frequency on the Modal page "
    u"(e.g. -10 x the expected first frequency); the rigid modes then "
    u"appear with ~0 Hz. That only helps if the loose part has mass "
    u"(density or point mass). A connection from #ASTER code only counts "
    u"via #ASTER-EXCIT.")
_MODES_FEW_HINT = (
    u"Code-Aster vond (te) weinig modes rond de zoekfrequentie. Kies een "
    u"zoekfrequentie dichter bij de verwachte frequenties (0 = de "
    u"laagste), of minder modes. Met zoekfrequentie 0 kunnen veel valse "
    u"modes bij ~0 Hz de oorzaak zijn (bijv. een membraan met een zeer "
    u"kleine E en voorspanning, KNOWN_BUGS 49): kies dan een zoekfrequentie "
    u"rond de verwachte frequentie.",
    u"Code-Aster found (too) few modes around the search frequency. Choose "
    u"a search frequency closer to the expected frequencies (0 = the "
    u"lowest), or fewer modes. With search frequency 0, many spurious "
    u"modes at ~0 Hz can be the cause (e.g. a membrane with a very small E "
    u"and prestress, KNOWN_BUGS 49): then choose a search frequency around "
    u"the expected frequency.")
# MODE_ITER_CYCL (bridge_cyclic.py): de controle van de cyclische ligging
# (gezien met Code-Aster 17.4: ALGORITH16_52 bij ruis rond 0 in de axiale
# coordinaat, ALGORITH16_58 bij verwisselde snijvlakken, ALGORITH14_95 bij 8
# segmenten met knopen op de as).
_CYCLIC_HINT = (
    u"Cyclische symmetrie (MODE_ITER_CYCL): de snijvlakken passen niet als "
    u"herhaling om de as. Controleer in Mecway de as, het aantal segmenten "
    u"en dat de knopen van master en slave op elkaar passen; bij 8 "
    u"segmenten met knopen op de as: maak een klein gat op de as "
    u"(KNOWN_BUGS 61).",
    u"Cyclic symmetry (MODE_ITER_CYCL): the cut faces do not repeat about "
    u"the axis. Check in Mecway the axis, the number of segments and that "
    u"the master and slave nodes match; with 8 segments and nodes on the "
    u"axis: make a small hole on the axis (KNOWN_BUGS 61).")
_MESS_HINTS = {"FACTOR_10": _SINGULAR_HINT, "FACTOR_11": _SINGULAR_HINT,
               "ALGORITH16_52": _CYCLIC_HINT, "ALGORITH16_58": _CYCLIC_HINT,
               "ALGORITH14_94": _CYCLIC_HINT,
               "ALGELINE3_81": _MODES_SINGULAR_HINT,
               "ALGELINE5_49": _MODES_FEW_HINT,
               "ALGELINE3_78": _MODES_FEW_HINT,
               "FACTOR_42": _SINGULAR_HINT, "FACTOR_57": _SINGULAR_HINT,
               "MECANONLINE9_7": _NEWTON_HINT,
               # "Arret pour cause de matrice non inversible" (STAT_NON_LINE
               # met contact, gezien in study.mess); oorzaken als FACTOR_*,
               # of een pivot op een contact-Lagrange-vrijheidsgraad.
               "MECANONLINE9_4": _NEWTON_SINGULAR_HINT}


def mess_hints():
    """Uitleg per Code-Aster-meldingscode in de gekozen taal (dict)."""
    return dict((code, tr(*pair)) for code, pair in _MESS_HINTS.items())
# Regel die de .comm print na een opgevangen ConvergenceError
# (bridge_to_aster._build_nonlin_static).
_MESS_PARTIAL_RE = re.compile(r"BRIDGE_ONVOLLEDIG\s+(\S+)")
# Het FIN-commando schrijft "<I> <FIN> ARRET NORMAL DANS "FIN" ..." (in elke
# afgeronde study.mess van de testmodellen, ook na een opgevangen
# ConvergenceError). Ontbreekt het, dan is study.comm niet tot het eind
# uitgevoerd (bijv. SyntaxError: study.mess heeft dan alleen
# "DIAGNOSTIC JOB : OK"; gezien met Code-Aster 17.4).
_MESS_FIN = "<FIN>"
# Regel die de .comm van de respons op modale basis print voordat hij stopt
# (bridge_dyna_modal._dyna_block): een mode zonder frequentie of een
# dempingsverhouding >= 1.
_MESS_MODAL_RE = re.compile(r"BRIDGE_DYNA_MODAL_FOUT\s+(\w+)\s*(.*)")
_MODAL_HINTS = {
    "FREQ": (u"Respons op modale basis: een mode met frequentie <= 0 ({0}): "
             u"een starre beweging of mechanisme. Mecway staat dat ook niet "
             u"toe (manual par. 6.4). Leg het model volledig vast.",
             u"Mode superposition: a mode with frequency <= 0 ({0}): a rigid "
             u"body motion or mechanism. Mecway does not allow it either "
             u"(manual section 6.4). Constrain the model fully."),
    "DEMPING": (u"Respons op modale basis: dempingsverhouding >= 1 ({0}), "
                u"overkritisch. Mecway vraagt xi < 1 voor elke mode (manual "
                u"par. 6.4). Verlaag de Rayleigh-demping (alpha dempt de "
                u"laagste, beta de hoogste modes; xi per mode staat in "
                u"study.mess, BRIDGE_DEMPING), de instelling modal_damping, "
                u"of het aantal modes.",
                u"Mode superposition: damping ratio >= 1 ({0}), overdamped. "
                u"Mecway requires xi < 1 for every mode (manual section "
                u"6.4). Lower the Rayleigh damping (alpha damps the lowest, "
                u"beta the highest modes; xi per mode is in study.mess, "
                u"BRIDGE_DEMPING), the setting modal_damping, or the number "
                u"of modes."),
}
# Idem voor de responsspectrumanalyse (bridge_spectrum._spectrum_block).
_MESS_SPECTRUM_RE = re.compile(r"BRIDGE_SPECTRUM_FOUT\s+(\w+)\s*(.*)")
_SPECTRUM_HINTS = {
    "FREQ": (u"Spectrumanalyse: een mode met frequentie <= 0 ({0}): een "
             u"starre beweging of mechanisme. Leg het model volledig vast.",
             u"Response spectrum analysis: a mode with frequency <= 0 ({0}): "
             u"a rigid body motion or mechanism. Constrain the model "
             u"fully."),
    "DEMPING": (u"Spectrumanalyse: dempingsverhouding >= 1 ({0}). Verlaag "
                u"spectrum_damping of Mecway's Rayleigh-demping (beta dempt "
                u"de hoogste modes), of het aantal modes.",
                u"Response spectrum analysis: damping ratio >= 1 ({0}). "
                u"Lower spectrum_damping or Mecway's Rayleigh damping (beta "
                u"damps the highest modes), or the number of modes."),
}
# Console-uitvoer van run_aster (stdout + stderr), geschreven door tee.
# Console-uitvoer van de run: <model>_<run_id>_console.txt (run_name).
CONSOLE_SUFFIX = "_console.txt"
# Een fout in aster_results_export.py (in de job) is een fout van de bridge,
# niet van het model; die module print dan "aster_results_export: FOUT"
# (taal nl) of "aster_results_export: ERROR" (en; de taal van de jobscripts
# volgt CONFIG_BRIDGE["language"] via study.comm, bridge_to_aster).
_EXPORT_FAIL = (u"aster_results_export: FOUT", u"aster_results_export: ERROR")
_EXPORT_HINT = (u"Het wegschrijven van de resultaten naar Mecway "
                u"(aster_results_export.py) mislukte: een fout in de bridge, "
                u"niet in het model. Zie study_console.txt in de werkmap.",
                u"Writing the results to Mecway (aster_results_export.py) "
                u"failed: an error in the bridge, not in the model. See "
                u"study_console.txt in the work folder.")
_PYTHON_HINT = (u"study.comm is niet (volledig) uitgevoerd. Een fout in "
                u"#ASTER-regels? Zie het Code-Aster-venster of "
                u"study_console.txt in de werkmap.",
                u"study.comm was not (fully) executed. An error in "
                u"#ASTER lines? See the Code-Aster window or "
                u"study_console.txt in the work folder.")


def _read_mess(mess_path):
    """Tekst van study.mess (UTF-8, onleesbare tekens vervangen)."""
    f = open(mess_path, "rb")
    try:
        return f.read().decode("utf-8", "replace")
    finally:
        f.close()


def _mess_diagnostic(text):
    """Tekst na 'DIAGNOSTIC JOB :' (bijv. '<A>_ALARM'), of ''."""
    if _MESS_DONE not in text:
        return u""
    rest = text[text.rindex(_MESS_DONE) + len(_MESS_DONE):].strip()
    return rest.splitlines()[0].strip() if rest else u""


def _diag_failed(diag):
    """True als de DIAGNOSTIC-regel een fout aangeeft (<F>, <S> of <E>)."""
    return diag[:3] in (u"<F>", u"<S>", u"<E>")


def _mess_finished(work_dir, since):
    """True als study.mess van deze run is en met DIAGNOSTIC JOB eindigt."""
    mess_path = os.path.join(work_dir, run_name(".mess"))
    try:
        # run_aster kopieert fort.6 pas aan het eind naar study.mess
        # ("Copying results" in het venster); 1 s marge voor de klok.
        if os.path.getmtime(mess_path) < since - 1.0:
            return False
        return _MESS_DONE in _read_mess(mess_path)
    except Exception:
        return False


def _mess_box_text(lines, start, max_chars=300):
    """Inhoud van het meldingskader na regel 'start', op een regel."""
    body = []
    for line in lines[start + 1:start + 40]:
        if _MESS_BOX_END in line:
            break
        part = line.strip().strip(_MESS_BOX_SIDE).strip()
        if part:
            body.append(part)
    text = u" ".join(body)
    if len(text) > max_chars:
        text = text[:max_chars] + u" ..."
    return text


# ==========================================
# CONSOLE-UITVOER EN CONVERGENTIESCHERM
# ==========================================
def aster_root():
    """Installatiemap van Code-Aster (exe_aster = <map>\\bin\\run_aster.bat)."""
    return os.path.dirname(os.path.dirname(
        os.path.abspath(exe_aster_path())))


def _launch_command(export_path, console_path):
    """
    cmd-regel voor run_aster in een eigen venster. 'cmd /k' houdt het
    venster open zodat de gebruiker fouten kan lezen; geneste quotes zijn
    nodig omdat de paden spaties kunnen bevatten (cmd haalt het buitenste
    paar weg). Met tools\\tee.exe uit de installatie gaat de uitvoer
    (stdout en stderr, 2>&1) ook naar console_path.

    Returns (commando, True als de uitvoer wordt vastgelegd).
    """
    run = '"{0}" "{1}"'.format(exe_aster_path(), export_path)
    tee = os.path.join(aster_root(), "tools", "tee.exe")
    if console_path and os.path.isfile(tee):
        return ('cmd /k "{0} 2>&1 | "{1}" "{2}""'.format(
            run, tee, console_path), True)
    return 'cmd /k "{0}"'.format(run), False


def _console_text(work_dir, since=None):
    """Tekst van de console-uitvoer van deze run, of ""."""
    path = os.path.join(work_dir, run_name(CONSOLE_SUFFIX))
    try:
        if not os.path.exists(path) or (
                since is not None and os.path.getmtime(path) < since - 1.0):
            return u""
        return _read_mess(path)
    except Exception:
        return u""


def console_errors(work_dir, since=None):
    """
    Python-fouten in study_console.txt van deze run (bijv. een SyntaxError
    in study.comm), via bridge_convergence.python_errors().

    Returns
    -------
    list of unicode - leeg als er geen fouten zijn of het bestand ontbreekt
                      of van een eerdere run is.
    """
    path = os.path.join(work_dir, run_name(CONSOLE_SUFFIX))
    try:
        if not os.path.exists(path):
            return []
        if since is not None and os.path.getmtime(path) < since - 1.0:
            return []
        return python_errors(_read_mess(path))
    except Exception as e:
        log(tr("WARNING: {0} niet leesbaar: {1}",
               "WARNING: {0} not readable: {1}").format(
                   os.path.basename(path), str(e)))
        return []


def _monitor_python():
    """pythonw.exe (zonder consolevenster) van de Code-Aster-installatie,
    anders python.exe; None als geen van beide bestaat."""
    root = aster_root()
    try:
        subdirs = sorted(d for d in os.listdir(root)
                         if d.lower().startswith("python"))
    except OSError:
        return None
    for exe in ("pythonw.exe", "python.exe"):
        for d in subdirs:
            path = os.path.join(root, d, exe)
            if os.path.isfile(path):
                return path
    return None


def latest_mess(work_dir, model=None):
    """Pad van het nieuwste .mess-bestand in work_dir (de laatste run), of
    None. model: naam van het Mecway-bestand (bridge_utils.model_stem); dan
    bij voorkeur de laatste run van dat model."""
    def mtime(p):
        try:
            return os.path.getmtime(p)
        except OSError:
            return 0.0
    try:
        found = [os.path.join(work_dir, n) for n in os.listdir(work_dir)
                 if n.lower().endswith(".mess")]
    except OSError:
        return None
    if model:
        own = [p for p in found
               if os.path.basename(p).startswith(model + "_")]
        found = own or found
    if not found:
        return None
    return max(found, key=mtime)


def _monitor_log(mess_path):
    """Logbestand van het convergentiescherm naast het .mess-bestand:
    <model>_<run_id>.mess -> <model>_<run_id>_monitor_log.txt."""
    stem = os.path.splitext(os.path.abspath(mess_path))[0]
    return stem + "_monitor_log.txt"


def start_convergence_monitor(mess_path, since=0.0, console_path=None):
    """
    Start het convergentiescherm (aster_convergence_monitor.py) als apart
    proces met de CPython van Code-Aster. Het proces is los van Mecway:
    het blijft werken als het Mecway-script klaar is, en sluiten ervan
    heeft geen invloed op de run.

    Parameters
    ----------
    mess_path    : str   - .mess-bestand uit de .export (live: de fort.6
                           van de run die daarheen schrijft).
    since        : float - time.time() bij de start van de run; 0 = achteraf
                           (het bestaande .mess-bestand tonen).
    console_path : str   - study_console.txt (Python-fouten), of None.

    Returns
    -------
    bool - True als het proces gestart is (fouten worden gelogd, nooit
           doorgegeven: de run gaat altijd door).
    """
    try:
        python = _monitor_python()
        script = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "aster_convergence_monitor.py")
        if python is None or not os.path.isfile(script):
            log(tr("WARNING: Convergentiescherm niet gestart: Python van "
                   "Code-Aster ({0}) of {1} niet gevonden.",
                   "WARNING: Convergence window not started: Python of "
                   "Code-Aster ({0}) or {1} not found.").format(
                       aster_root(), script))
            return False
        work_dir = os.path.dirname(os.path.abspath(mess_path))
        args = [python, script, "--mess", mess_path,
                # %-opmaak: altijd een punt (landinstelling, zie bridge_gui)
                "--since", "%.3f" % since,
                "--log", _monitor_log(mess_path),
                "--lang", tr("nl", "en")]
        if console_path:
            args += ["--console", console_path]
        # PYTHONHOME/PYTHONPATH van een andere Python mogen de Python van
        # Code-Aster niet in de war brengen.
        env = dict((k, v) for k, v in os.environ.items()
                   if k.upper() not in ("PYTHONHOME", "PYTHONPATH"))
        proc = subprocess.Popen(args, cwd=work_dir, env=env,
                                creationflags=0x08000000)  # CREATE_NO_WINDOW
        log("INFO: Convergentiescherm gestart (pid {0}).".format(proc.pid))
        return True
    except Exception as e:
        log(tr("WARNING: Convergentiescherm niet gestart: {0}",
               "WARNING: Convergence window not started: {0}").format(
                   str(e)))
        return False


# ==========================================
# START CODE-ASTER
# ==========================================
def run_AsterSolver(export_path, work_dir, resu_name, mw):
    """
    Start Code-Aster via de run_aster.bat launcher in een apart cmd-venster.

    Parameters
    ----------
    export_path : str  - Pad naar het door bridge_to_aster.py geschreven
                         study.export-bestand.
    work_dir    : str  - Werkmap (cwd voor het subprocess).
    resu_name   : str  - Bestandsnaam van het .resu-bestand (uniek per run
                         door de tijdstempel; alleen gebruikt voor de log).
    mw          : obj  - Mecway API-object voor mw.message() (synchrone
                         dialoog terwijl Code-Aster draait).

    Returns
    -------
    bool - True als het subprocess succesvol gestart is, False bij fout.
           Let op: deze return-waarde zegt NIETS over of de solver
           succesvol heeft gerekend.
    """
    log("--- Running Code_Aster ---")
    exe = exe_aster_path()
    if not os.path.isfile(exe):
        # 'cmd /k' start altijd; zonder deze controle zou het venster met
        # "niet herkend" verschijnen en daarna gesloten worden.
        log(tr("ERROR: Code-Aster-launcher niet gevonden: {0} (instelling "
               "op de pagina Algemeen).",
               "ERROR: Code-Aster launcher not found: {0} (setting on the "
               "General page).").format(exe))
        return False

    # Het .resu-bestand wordt NIET vooraf aangemaakt: ontbreekt het na de
    # run, dan meldt AsterToMecway_results() dat Code-Aster niets schreef.
    # Alleen de tabelroute schrijft een .resu (bridge_to_aster.
    # _export_result_lines); de liml-route een resultaatfragment.
    if CONFIG_BRIDGE.get("result_route", "liml") == "tabel":
        log("INFO: Verwacht resultaatbestand: {0}".format(resu_name))

    try:
        console_path = os.path.join(work_dir, run_name(CONSOLE_SUFFIX))
        try:
            if os.path.exists(console_path):
                os.remove(console_path)   # niet verwarren met een oude run
        except OSError:
            pass
        cmd_command, captured = _launch_command(export_path, console_path)
        if not captured:
            log(tr("WARNING: tee.exe niet gevonden in de "
                   "Code-Aster-installatie; de console-uitvoer "
                   "(Python-fouten) wordt niet vastgelegd.",
                   "WARNING: tee.exe not found in the Code-Aster "
                   "installation; the console output (Python errors) is "
                   "not captured."))

        start_time = time.time()
        process = subprocess.Popen(
            cmd_command,
            cwd=work_dir,
            creationflags=0x00000010,  # CREATE_NEW_CONSOLE
        )

        monitor = False
        # Ook bij lineaire analyses: dan de statistiek van de lineaire
        # oplosser (MECA_STATIQUE met INFO=2), de meldingen en de diagnose.
        if CONFIG_BRIDGE.get("convergence_monitor", True):
            monitor = start_convergence_monitor(
                os.path.join(work_dir, run_name(".mess")), start_time,
                console_path if captured else None)

        mw.message(
            tr("Code_Aster is gestart in een apart venster.\n",
               "Code_Aster has been started in a separate window.\n") +
            (tr("Het convergentiescherm toont de voortgang.\n",
                "The convergence window shows the progress.\n") if monitor
             else "") +
            tr("Klik op OK zodra de berekening klaar is.",
               "Click OK once the calculation has finished.")
        )
        # Te vroeg op OK geklikt? Dan nog een keer vragen, want
        # process.terminate() hieronder sluit het venster (gezien in Mecway).
        if not _mess_finished(work_dir, start_time):
            log(tr("WARNING: study.mess is na de eerste OK nog niet "
                   "afgesloten.",
                   "WARNING: study.mess was not yet closed after the first "
                   "OK."))
            mw.message(run_names(tr(
                "Code-Aster is nog niet klaar (study.mess is niet "
                "afgesloten).\nKlik op OK zodra in het venster "
                "'DIAGNOSTIC JOB' staat.",
                "Code-Aster has not finished yet (study.mess is not "
                "closed).\nClick OK once the window shows "
                "'DIAGNOSTIC JOB'."
            )))

        # Bij een abort of een Python-fout (geen <FIN>) blijft het venster
        # open: de traceback staat daar (en in study_console.txt).
        mess_path = os.path.join(work_dir, run_name(".mess"))
        try:
            if (not os.path.exists(mess_path) or
                    os.path.getmtime(mess_path) < start_time - 1.0):
                # Geen .mess van deze run: Code-Aster is niet (goed)
                # gestart; de oorzaak staat alleen in het venster.
                aborted = True
            else:
                text = _read_mess(mess_path)
                aborted = (_diag_failed(_mess_diagnostic(text)) or
                           (_MESS_DONE in text and _MESS_FIN not in text))
        except Exception:
            aborted = False
        if not aborted and console_errors(work_dir, start_time):
            aborted = True
        if aborted:
            log(tr("WARNING: Code-Aster afgebroken; het venster blijft open.",
                   "WARNING: Code-Aster aborted; the window stays open."))
            return True

        # process.terminate() sluit het cmd-venster. Dat mag: fouten uit
        # study.mess meldt main() via aster_mess_errors(). De Popen-handle
        # wijst alleen naar het cmd-proces van het aparte venster, niet naar
        # de Code-Aster-processen die daaruit gestart zijn. Een fout (bijv.
        # venster al gesloten) wordt genegeerd.
        try:
            process.terminate()
        except Exception:
            pass
        return True

    except Exception as e:
        log(tr("ERROR: Kon Code-Aster niet starten: {0}",
               "ERROR: Could not start Code-Aster: {0}").format(str(e)))
        return False


# ==========================================
# FOUTEN IN study.mess
# ==========================================
def aster_mess_errors(work_dir, since=None):
    """
    Zoek fouten (<EXCEPTION>, <F>, <S>, <E>) in study.mess.

    Parameters
    ----------
    work_dir : str          - Werkmap met study.mess.
    since    : float of None - Starttijd van de run (time.time()). Is
                               study.mess ouder, dan hoort het bij een eerdere
                               run en wordt het niet beoordeeld.

    Returns
    -------
    list of str - Unieke meldingen in volgorde van voorkomen, bijv.
                  "<EXCEPTION> <CALCULEL_58>"; leeg als er geen fouten zijn
                  of het bestand ontbreekt/verouderd is (dat wordt gelogd).
    """
    mess_path = os.path.join(work_dir, run_name(".mess"))
    if not os.path.exists(mess_path) or (
            since is not None and os.path.getmtime(mess_path) < since):
        # Geen .mess van deze run: Code-Aster startte niet of stopte
        # direct; de oorzaak staat in het Code-Aster-venster (dat blijft
        # open, run_AsterSolver).
        msg = run_names(tr(u"Code-Aster schreef geen study.mess voor deze "
                           u"run; de oorzaak staat in het Code-Aster-venster "
                           u"(bijv. een fout in de installatie of de "
                           u"werkmap).",
                           u"Code-Aster wrote no study.mess for this run; "
                           u"the cause is shown in the Code-Aster window "
                           u"(e.g. an error in the installation or the work "
                           u"folder)."))
        log(u"ERROR: " + msg)
        return [msg]
    hints = mess_hints()
    python_hint = tr(*_PYTHON_HINT)
    errors = []
    seen = set()
    alarms = []
    try:
        text = _read_mess(mess_path)
        lines = text.splitlines()
        for i, line in enumerate(lines):
            m = _MESS_ERROR_RE.search(line)
            if m:
                code = u"<{0}> <{1}>".format(m.group(1), m.group(2))
                if code not in seen:
                    seen.add(code)
                    msg = u"{0}: {1}".format(code, _mess_box_text(lines, i))
                    if m.group(2) in hints:
                        msg += u" -> " + hints[m.group(2)]
                    partial = _MESS_PARTIAL_RE.search(text)
                    if partial and m.group(2).startswith("MECANONLINE"):
                        msg += tr(u" Laatste tijdstip met resultaten: {0} s.",
                                  u" Last time with results: {0} s.").format(
                            partial.group(1))
                    errors.append(msg)
                continue
            m = _MESS_ALARM_RE.search(line)
            if m and m.group(1) not in alarms:
                alarms.append(m.group(1))
        modal = _MESS_MODAL_RE.search(text)
        hints = _MODAL_HINTS
        if not modal:
            modal = _MESS_SPECTRUM_RE.search(text)
            hints = _SPECTRUM_HINTS
        if modal:
            pair = hints.get(modal.group(1))
            if pair:
                errors.append(tr(*pair).format(modal.group(2).strip()))
        # Python-fout in study.comm: staat alleen in de console-uitvoer.
        console = _console_text(work_dir, since)
        export_failed = any(key in console for key in _EXPORT_FAIL)
        for err in ([] if modal else console_errors(work_dir, since)):
            if export_failed:
                errors.append(u"{0} -> {1}".format(err, tr(*_EXPORT_HINT)))
            else:
                errors.append(tr(u"Python-fout in study.comm: {0} -> {1}",
                                 u"Python error in study.comm: {0} -> "
                                 u"{1}").format(err, python_hint))
        if not errors and _MESS_UNCAUGHT in text:
            errors.append(tr(u"Code-Aster stopte op een exceptie (\"Trying "
                             u"to close the database after an uncaught "
                             u"exception\"); zie study.mess.",
                             u"Code-Aster stopped on an exception (\"Trying "
                             u"to close the database after an uncaught "
                             u"exception\"); see study.mess."))
        if _MESS_DONE not in text:
            errors.append(tr(u"study.mess is niet afgesloten: Code-Aster is "
                             u"niet klaar of afgebroken.",
                             u"study.mess is not closed: Code-Aster has not "
                             u"finished or was aborted."))
        elif not errors and _diag_failed(_mess_diagnostic(text)):
            # Een crash in Python of C (bijv. "Windows fatal exception:
            # access violation", gezien in Mecway) geeft geen meldingskader,
            # alleen "DIAGNOSTIC JOB : <F>_ABNORMAL_ABORT". De traceback
            # staat dan alleen in het Code-Aster-venster.
            errors.append(tr(u"DIAGNOSTIC JOB : {0}. Er staat geen melding "
                             u"in study.mess; de oorzaak staat in het "
                             u"Code-Aster-venster.",
                             u"DIAGNOSTIC JOB : {0}. There is no message in "
                             u"study.mess; the cause is shown in the "
                             u"Code-Aster window.").format(
                                 _mess_diagnostic(text)))
        elif not errors and _MESS_FIN not in text:
            # Zonder vastgelegde console-uitvoer (geen tee.exe) is dit het
            # enige teken van een Python-fout (SyntaxError: DIAGNOSTIC OK).
            errors.append(tr(u"Code-Aster kwam niet tot FIN (geen '<FIN>' "
                             u"in study.mess, DIAGNOSTIC JOB : {0}): {1}",
                             u"Code-Aster did not reach FIN (no '<FIN>' in "
                             u"study.mess, DIAGNOSTIC JOB : {0}): {1}").format(
                                 _mess_diagnostic(text), python_hint))
    except Exception as e:
        log(tr("WARNING: study.mess niet leesbaar: {0}",
               "WARNING: study.mess not readable: {0}").format(str(e)))
        return []
    for alarm in alarms:
        if alarm in _ALARM_HINTS:
            log(tr(u"WARNING: Code-Aster-alarm {0} in study.mess: {1}",
                   u"WARNING: Code-Aster alarm {0} in study.mess: {1}").format(
                       alarm, tr(*_ALARM_HINTS[alarm])))
            continue
        log(tr(u"WARNING: Code-Aster-alarm {0} in study.mess.",
               u"WARNING: Code-Aster alarm {0} in study.mess.").format(alarm))
    # Meldingen met de bestandsnamen van deze run (getoond in mw.message)
    errors = [run_names(msg) for msg in errors]
    for msg in errors:
        log(tr(u"ERROR: Code-Aster meldt {0}",
               u"ERROR: Code-Aster reports {0}").format(msg))
    return errors
