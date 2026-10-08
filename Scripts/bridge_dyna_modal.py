# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# bridge_dyna_modal.py
# Verantwoordelijkheid: respons op modale basis, Mecway Dynamic Response
# (D30) met de optie Mode superposition (manual v33 par. 6.4), naar
# Code-Aster:
#   ASSEMBLAGE (RIGI_MECA, MASS_MECA)  -> CALC_MODES (de laagste maxmode
#   modes) -> belastingsvectoren (CALC_VECT_ELEM + ASSE_VECTEUR) ->
#   PROJ_BASE (gegeneraliseerde massa, stijfheid en krachten) ->
#   DYNA_VIBRA TYPE_CALCUL='TRAN', BASE_CALCUL='GENE' (modale demping,
#   NEWMARK) -> REST_GENE_PHYS (verplaatsingen en snelheden op de knopen)
#   -> CALC_CHAMP (spanningen) -> aster_results_export (tijdstappen in de
#   .liml, zoals Mecway's eigen oplossing).
# Syntax: Code-Aster-testcase forma12d (17.4, share/aster/tests) en de
# catalogi dyna_vibra.py, proj_base.py, rest_gene_phys.py, calc_vect_elem.py
# van code_aster 17.4.
#
#   build_dyna_modal(...) - builder voor SOLVER_BUILDERS_ASTER["DYNA_VIBRA"]
#                           (bridge_to_aster._build_dyna_modal).
#
# Gedrag van Mecway (headless gemeten met Mecway 33, mecway.exe <model>
# solve; verificatie MD in docs/Verificatie.md):
#   - De modale vergelijkingen worden exact geintegreerd: een staplast op
#     een massa aan een balk gaf met stepsize 0,02 s (omega dt = 0,71) en
#     0,001 s dezelfde respons, binnen 1e-4 van u_st (1 - cos omega t).
#     Code-Aster integreert met NEWMARK (gemiddelde versnelling); de bridge
#     kiest de rekenstap zo klein dat de periodefout verwaarloosbaar is
#     (steps_per_period per periode van de hoogste mode in de basis; de
#     periodefout van NEWMARK is (omega dt)^2 / 12, Hughes 2000, par. 9.3).
#   - Geen statische correctie voor de weggelaten modes: met 1 mode bleef
#     de eindverplaatsing van een uitkragende balk (puntlast op de tip,
#     gedempt) 0,9707 x de statische (de bijdrage van mode 1), met 6 modes
#     0,9986. De instelling static_correction (standaard uit) voegt die
#     correctie toe (MACRO_ELAS_MULT als MODE_CORR, CORR_STAT, testcase
#     forma12d).
#   - Rayleigh-demping (alleen op het hele model, manual par. 9.26) als
#     modale dempingsverhouding xi_i = (alpha / omega_i + beta omega_i) / 2
#     (manual par. 9.26); xi_i < 1 voor elke mode (manual par. 6.4).
#   - Uitvoer: t = 0 en elke keepeverynthtimestep-de stap van stepsize
#     (bridge_to_aster._output_times); displacement, rotation, velocity en
#     von Mises; geen versnellingen (manual par. 6.4).
#
# Compatibel met IronPython 2.7.9 (alleen standaardbibliotheek).

import math
import os
import re
from collections import OrderedDict

import bridge_to_aster as b2a
from bridge_utils import log, tr, CONFIG_BRIDGE, run_name


# Velden van de fysieke respons (REST_GENE_PHYS NOM_CHAM) per groep uit
# "fields"; reacties en plastische rek horen niet bij deze analyse.
_STRESS_GROUPS = ("stresses", "vonmises")
_FIELD_GROUPS = ("displacements", "velocities") + _STRESS_GROUPS


def _value_names(text):
    """Namen van concepten die als waarde in een _F(...) staan (KW=naam);
    een groepsnaam staat tussen aanhalingstekens en telt niet mee."""
    return set(re.findall(r"=\s*([A-Za-z_]\w*)", text))


def _substitute(kws, mapping):
    """Kopie van een excit-dict (_build_charges) met de functienamen uit
    mapping vervangen (alleen als waarde: na '=', niet in een groepsnaam)."""
    out = OrderedDict()
    for kw, entries in kws.items():
        new = []
        for e in entries:
            for old, rep in mapping.items():
                e = re.sub(r"(=\s*){0}\b".format(re.escape(old)),
                           r"\g<1>{0}".format(rep), e)
            new.append(e)
        out[kw] = new
    return out


def _excitations(charges, timed):
    """
    Splits de lasten van _build_charges in vectoren met elk een
    tijdsfunctie: F(t) = sum_j V_j f_j(t) (lineair in de functiewaarden).
    timed: namen van de DEFI_FONCTION's in de tijd. Returns (opleggingen,
    excitaties, extra lasten, fouten):
      opleggingen : namen van de lasten zonder belasting (load_k)
      excitaties  : [(last, functie of None)] (None = constant vanaf t = 0)
      extra       : [(naam, commando, excit-dict)] nieuwe lasten (gload_*):
                    een AFFE_CHAR_MECA_F met een functie op 1 en de andere
                    tijdsfuncties op 0
      fouten      : lasten met een tijdsafhankelijke oplegging
    """
    constraints, excitations, extra, errors = [], [], [], []
    for k, (cmd, kws, fonc) in enumerate(charges):
        name = "load_{0}".format(k + 1)
        if not b2a._is_load_charge((cmd, kws, fonc)):
            if cmd == "AFFE_CHAR_MECA_F":
                errors.append(name)
            constraints.append(name)
            continue
        if cmd == "AFFE_CHAR_MECA":
            excitations.append((name, fonc))
            continue
        text = " ".join(e for es in kws.values() for e in es)
        names = _value_names(text)
        used = [f for f in timed if f in names]
        others = sorted(n for n in names if n not in timed)
        if others:
            # Constanten (DEFI_CONSTANTE van balklasten) en plaatsfuncties
            # (hydrostatische druk): vanaf t = 0, de tijdsfuncties op 0.
            gname = "gload_{0}_0".format(k + 1)
            extra.append((gname, cmd, _substitute(
                kws, dict((f, "bridge_f0") for f in used))))
            excitations.append((gname, None))
        for i, f in enumerate(used):
            gname = "gload_{0}_{1}".format(k + 1, i + 1)
            mapping = dict((g, "bridge_f0") for g in used + others)
            mapping[f] = "bridge_f1"
            extra.append((gname, cmd, _substitute(kws, mapping)))
            excitations.append((gname, f))
    return constraints, excitations, extra, errors


def _check_model(model_data, analysis):
    """Wat de modale respons niet kan: ERROR met de oplossing. Returns
    True als het model te rekenen is."""
    ok = True
    if not b2a._check_linear(model_data, tr("een respons op modale basis",
                                            "a mode superposition analysis")):
        ok = False
    dampers = [e for e, s in model_data.get("discrete_sections", {}).items()
               if s.get("kind") == "damper"]
    if dampers:
        # Mecway: "Damper elements cannot be used" (manual v33 par. 6.4).
        log(tr("ERROR: Dempers ({0}) kunnen niet in een respons op modale "
               "basis (ook niet in Mecway, manual par. 6.4). Gebruik "
               "Rayleigh-demping, of de directe analyse dyna_impl_3d "
               "(Mecway: Mode superposition uit).",
               "ERROR: Dampers ({0}) cannot be used in a mode superposition "
               "analysis (not in Mecway either, manual section 6.4). Use "
               "Rayleigh damping, or the direct analysis dyna_impl_3d "
               "(Mecway: Mode superposition off).").format(
                   ", ".join(dampers)))
        ok = False
    if model_data.get("initial_velocities"):
        log(tr("ERROR: Een beginsnelheid wordt in de respons op modale basis "
               "niet vertaald. Gebruik de directe analyse dyna_impl_3d.",
               "ERROR: An initial velocity is not translated in the mode "
               "superposition analysis. Use the direct analysis "
               "dyna_impl_3d."))
        ok = False
    if b2a._thermal_active(model_data):
        log(tr("ERROR: Thermische rek (Thermal stress) wordt in de respons "
               "op modale basis niet vertaald. Gebruik dyna_impl_3d.",
               "ERROR: Thermal strain (Thermal stress) is not translated in "
               "the mode superposition analysis. Use dyna_impl_3d."))
        ok = False
    if CONFIG_BRIDGE.get("result_route", "liml") != "liml":
        log(tr("ERROR: De resultaatroute 'tabel' werkt niet voor de respons "
               "op modale basis. Kies 'liml' (pagina Algemeen).",
               "ERROR: The result route 'tabel' does not work for the mode "
               "superposition analysis. Choose 'liml' (General page)."))
        ok = False
    if model_data.get("analysis_type") == "D30" and not model_data.get(
            "mode_superposition"):
        log(tr("WARNING: Het Mecway-model rekent Dynamic Response zonder Mode "
               "superposition (Newmark op het hele model); de bridge rekent "
               "op modale basis met {0} modes. Voor de directe berekening: "
               "dyna_impl_3d.",
               "WARNING: The Mecway model computes Dynamic Response without "
               "Mode superposition (Newmark on the whole model); the bridge "
               "computes on a modal basis with {0} modes. For the direct "
               "computation: dyna_impl_3d.").format(
                   b2a._mode_count(model_data, analysis)))
    return ok


def _regular_outputs(out):
    """(interval, aantal) als de uitvoertijden out = k * interval (k = 1..n)
    zijn, anders None. DYNA_VIBRA GENE rekent met een vaste stap (PAS); de
    uitvoertijden moeten daarop vallen."""
    if not out:
        return None
    inc = out[0]
    for k, t in enumerate(out):
        if abs(t - (k + 1) * inc) > 1.0e-6 * inc:
            return None
    return inc, len(out)


def _damping_text(model_data, analysis):
    """(alpha [1/s], beta [s], extra modale dempingsverhouding)."""
    ray = model_data.get("rayleigh") or {}
    zeta = float(analysis.get("modal_damping", 0.0) or 0.0)
    return float(ray.get("alpha", 0.0) or 0.0), \
        float(ray.get("beta", 0.0) or 0.0), zeta


def _dyna_block(analysis, excitations, inc, n_out, alpha, beta, zeta,
                static_corr, fields):
    """
    Demping, tijdstap, DYNA_VIBRA en REST_GENE_PHYS (Python in de .comm:
    de frequenties zijn pas na CALC_MODES bekend).
    """
    nsteps = max(4, int(analysis.get("steps_per_period", 50)))
    lines = [
        "# Modale demping: Mecway's Rayleigh-demping C = alpha M + beta K is",
        "# per mode xi_i = (alpha / w_i + beta w_i) / 2 (Mecway manual par.",
        "# 9.26), plus de instelling modal_damping.",
        "import math as _math",
        "bridge_freqs = [float(_f) for _f in "
        "modes.getAccessParameters()['FREQ']]",
        "bridge_alpha, bridge_beta, bridge_zeta = {0!r}, {1!r}, {2!r}".format(
            alpha, beta, zeta),
        "bridge_xi = []",
        "for _f in bridge_freqs:",
        "    _w = 2.0 * _math.pi * _f",
        "    if _w <= 0.0:",
        "        # Starre beweging of mechanisme: Mecway staat dat ook niet toe",
        "        # (manual par. 6.4); bridge_run_aster meldt de oorzaak.",
        "        print('BRIDGE_DYNA_MODAL_FOUT FREQ mode=%d f=%g' % "
        "(len(bridge_xi) + 1, _f))",
        "        raise ValueError('mode superposition: mode met f <= 0')",
        "    bridge_xi.append(0.5 * (bridge_alpha / _w + bridge_beta * _w) + "
        "bridge_zeta)",
        "for _i, (_f, _x) in enumerate(zip(bridge_freqs, bridge_xi)):",
        "    print('BRIDGE_DEMPING mode %d: %.6g Hz, xi = %.6g' % "
        "(_i + 1, _f, _x))",
        "if max(bridge_xi) >= 1.0:",
        "    # Mecway: xi < 1 voor elke mode (manual par. 6.4, 9.26).",
        "    print('BRIDGE_DYNA_MODAL_FOUT DEMPING mode=%d xi=%g' % "
        "(bridge_xi.index(max(bridge_xi)) + 1, max(bridge_xi)))",
        "    raise ValueError('mode superposition: dempingsverhouding >= 1')",
        "",
        "# Rekenstap: {0} stappen per periode van de hoogste mode (NEWMARK,".format(
            nsteps),
        "# periodefout (w dt)^2 / 12), een geheel aantal per uitvoerinterval.",
        "bridge_inc, bridge_nout = {0!r}, {1}".format(inc, n_out),
        "bridge_nsub = max(1, int(_math.ceil(bridge_inc * max(bridge_freqs) * {0})))"
        .format(nsteps),
        "bridge_pas = bridge_inc / bridge_nsub",
        "print('BRIDGE_DYNA_MODAL rekenstap %g s (%d per uitvoerinterval), "
        "%d stappen' % (bridge_pas, bridge_nsub, bridge_nsub * bridge_nout))",
        "",
    ]
    excit = []
    for k, (charge, fonc) in enumerate(excitations):
        f = fonc or "bridge_f1"
        entry = "_F(VECT_ASSE_GENE=fgen_{0}, FONC_MULT={1}".format(k + 1, f)
        if static_corr:
            entry += ", CORR_STAT='OUI', D_FONC_DT={0}, D_FONC_DT2={1}".format(
                "dfdt_{0}".format(k + 1) if fonc else "bridge_f0",
                "d2fdt_{0}".format(k + 1) if fonc else "bridge_f0")
        excit.append(entry + ")")
    if static_corr:
        # Statische correctie (testcase forma12d): de statische respons op
        # elke belastingsvector (MACRO_ELAS_MULT) als MODE_CORR, met de
        # eerste en tweede afgeleide van de tijdsfunctie (CALC_FONCTION
        # DERIVE; voor een constante last 0).
        for k, (charge, fonc) in enumerate(excitations):
            if fonc:
                lines += [
                    "dfdt_{0} = CALC_FONCTION(DERIVE=_F(FONCTION={1}))".format(
                        k + 1, fonc),
                    "d2fdt_{0} = CALC_FONCTION(DERIVE=_F(FONCTION=dfdt_{0}))"
                    .format(k + 1),
                ]
    lines += [
        "tran = DYNA_VIBRA(",
        "    TYPE_CALCUL='TRAN',",
        "    BASE_CALCUL='GENE',",
        "    SCHEMA_TEMPS=_F(SCHEMA='NEWMARK'),",
        "    MATR_MASS=mgen,",
        "    MATR_RIGI=kgen,",
        "    AMOR_MODAL=_F(AMOR_REDUIT=bridge_xi),",
        "    INCREMENT=_F(PAS=bridge_pas, INST_FIN=bridge_inc * bridge_nout),",
        "    ARCHIVAGE=_F(PAS_ARCH=bridge_nsub),",
        "    EXCIT=(",
    ] + ["        {0},".format(e) for e in excit] + ["    ),"]
    if static_corr:
        lines.append("    MODE_CORR=corr,")
    lines += [")", ""]
    cham = ["DEPL"]
    if "velocities" in fields:
        cham.append("VITE")
    if any(f in fields for f in _STRESS_GROUPS):
        # SIPM_ELNO (balken) op een DYNA_TRANS zonder versnellingen gaf
        # alarm CALCULEL3_1 "Il manque les accélérations" (Code-Aster 17.4,
        # gezien); de versnellingen gaan niet naar Mecway (zoals Mecway).
        cham.append("ACCE")
    lines += [
        "res = REST_GENE_PHYS(",
        "    RESU_GENE=tran,",
        "    TOUT_INST='OUI',",
        "    NOM_CHAM=({0},),".format(", ".join("'{0}'".format(c)
                                               for c in cham)),
    ]
    if static_corr:
        lines.append("    CORR_STAT='OUI',")
    lines += [")", ""]
    return "\n".join(lines)


def build_dyna_modal(work_dir, model_data, msh_path, resu_name, analysis):
    """
    Genereer de .comm en .export voor de respons op modale basis
    (analyse dyna_modal_3d, solver_command DYNA_VIBRA). Returns het pad van
    de .export, of None (ERROR gelogd).
    """
    log("--- Generating Aster Files (DYNA_VIBRA GENE) ---")
    if not _check_model(model_data, analysis):
        return None
    comm_path = os.path.join(work_dir, run_name(".comm"))

    load_funcs = b2a._collect_load_funcs(model_data)
    func_names = b2a._load_func_names(model_data)
    timed = [name for name, _h in load_funcs]
    charges = b2a._build_charges(model_data, func_names)
    constraints, excitations, extra, errors = _excitations(charges, timed)
    if errors:
        log(tr("ERROR: Een tijdsafhankelijke opgelegde verplaatsing kan niet "
               "in de respons op modale basis (de modes hebben die "
               "vrijheidsgraden vast). Gebruik dyna_impl_3d.",
               "ERROR: A time-dependent imposed displacement is not possible "
               "in the mode superposition analysis (the modes have those "
               "degrees of freedom fixed). Use dyna_impl_3d."))
        return None
    nonzero = b2a._nonzero_imposed(model_data)
    if nonzero:
        log(tr("ERROR: Opgelegde verplaatsing of rotatie <> 0 op {0}: dat "
               "kan niet in de respons op modale basis (de modes hebben die "
               "vrijheidsgraden vast). Zet de waarde op 0, of gebruik "
               "dyna_impl_3d.",
               "ERROR: Imposed displacement or rotation <> 0 on {0}: not "
               "possible in the mode superposition analysis (the modes have "
               "those degrees of freedom fixed). Set the value to 0, or use "
               "dyna_impl_3d.").format(
                   ", ".join("'{0}'".format(s) for s in nonzero)))
        return None
    custom_excit = model_data.get("aster_custom", {}).get("excit", [])
    if not constraints and not custom_excit:
        # Mecway: "Rigid body motion is not allowed so the model must be
        # fully constrained" (manual v33 par. 6.4).
        log(tr("ERROR: Respons op modale basis zonder opleggingen: starre "
               "bewegingen kunnen niet (ook niet in Mecway, manual par. "
               "6.4). Leg het model in Mecway vast.",
               "ERROR: Mode superposition analysis without supports: rigid "
               "body motion is not possible (not in Mecway either, manual "
               "section 6.4). Add supports in Mecway."))
        return None
    if not excitations:
        log(tr("WARNING: Geen belastingen: de respons is overal 0.",
               "WARNING: No loads: the response is zero everywhere."))

    out = b2a._output_times(model_data, analysis)
    reg = _regular_outputs(out)
    if reg is None:
        log(tr("ERROR: De uitvoertijden ({0} ... {1} s) liggen niet op een "
               "regelmatig raster vanaf 0; DYNA_VIBRA rekent met een vaste "
               "stap. Kies in Mecway een eindtijd die een veelvoud is van "
               "de tijdstap.",
               "ERROR: The output times ({0} ... {1} s) are not on a regular "
               "grid from 0; DYNA_VIBRA uses a constant step. Choose an end "
               "time in Mecway that is a multiple of the time step.").format(
                   out[0] if out else "?", out[-1] if out else "?"))
        return None
    inc, n_out = reg
    log("INFO: Respons op modale basis: {0} uitvoertijden, interval {1} s, "
        "eindtijd {2} s.".format(n_out, inc, out[-1]))

    fields = [f for f in analysis.get("fields", ["displacements"])
              if f in _FIELD_GROUPS]
    ignored = [f for f in analysis.get("fields", []) if f not in
               _FIELD_GROUPS]
    if "displacements" not in fields:
        fields.insert(0, "displacements")
    if ignored:
        log("INFO: Respons op modale basis: velden {0} worden niet "
            "berekend (geen reacties of plastische rek).".format(
                ", ".join(ignored)))
    static_corr = bool(analysis.get("static_correction", False))
    alpha, beta, zeta = _damping_text(model_data, analysis)

    defigroup_block = b2a._build_defigroup_block(model_data)
    model_block = b2a._build_model_block(model_data)
    cara_elem_block = b2a._build_cara_elem_block(model_data)
    cara_ref = "    CARA_ELEM=cara,\n" if b2a._has_cara(model_data) else ""
    materials_block = b2a._build_materials_block(model_data, analysis)
    fieldmat_block = b2a._build_fieldmat_block(model_data)
    load_funcs_block = b2a._build_load_funcs_block(load_funcs)
    loads_block, _excit = b2a._render_charges(charges, model_data)
    consts = ("bridge_f0 = DEFI_CONSTANTE(VALE=0.0)\n"
              "bridge_f1 = DEFI_CONSTANTE(VALE=1.0)\n")
    loads_block = consts + loads_block + "\n" + "".join(
        b2a._render_char_meca(n, c, k) for n, c, k in extra)

    charge_names = list(constraints) + [n for n in custom_excit
                                        if n not in constraints]
    charge_ref = ""
    if charge_names:
        charge_ref = "    CHARGE=({0},),\n".format(", ".join(charge_names))
    assemblage_block = (
        "ASSEMBLAGE(\n"
        "    MODELE=model,\n"
        "    CHAM_MATER=fieldmat,\n"
        "{cara}{charge}"
        "    NUME_DDL=CO('numedd'),\n"
        "    MATR_ASSE=(\n"
        "        _F(MATRICE=CO('rigidite'), OPTION='RIGI_MECA'),\n"
        "        _F(MATRICE=CO('masse'), OPTION='MASS_MECA'),\n"
        "    ),\n"
        ")\n").format(cara=cara_ref, charge=charge_ref)

    scripts = os.path.dirname(os.path.abspath(__file__))
    calc_modes_block = b2a._calc_modes_block(
        model_data, analysis, False, cara_ref).replace(
            "res = CALC_MODES(", "modes = CALC_MODES(", 1)
    calc_modes_block += (
        "# Modeoverzicht (frequentie, effectieve massa) in de .mess\n"
        "import sys\n"
        "sys.path.insert(0, {0})\n"
        "import aster_results_export\n"
        "aster_results_export.print_modes(modes, mass=masse)\n").format(
            b2a._py_str(scripts))

    proj = []
    for k, (charge, _fonc) in enumerate(excitations):
        proj += [
            "vel_{0} = CALC_VECT_ELEM(OPTION='CHAR_MECA', CHARGE={1}, "
            "CHAM_MATER=fieldmat{2})".format(
                k + 1, charge, ", CARA_ELEM=cara" if cara_ref else ""),
            "vec_{0} = ASSE_VECTEUR(VECT_ELEM=vel_{0}, NUME_DDL=numedd)"
            .format(k + 1),
        ]
    proj += [
        "PROJ_BASE(",
        "    BASE=modes,",
        "    STOCKAGE='DIAG',",
        "    MATR_ASSE_GENE=(",
        "        _F(MATRICE=CO('mgen'), MATR_ASSE=masse),",
        "        _F(MATRICE=CO('kgen'), MATR_ASSE=rigidite),",
        "    ),",
    ]
    if excitations:
        proj.append("    VECT_ASSE_GENE=(")
        proj += ["        _F(VECTEUR=CO('fgen_{0}'), VECT_ASSE=vec_{0}, "
                 "TYPE_VECT='FORC'),".format(k + 1)
                 for k in range(len(excitations))]
        proj.append("    ),")
    proj += [")", ""]
    if static_corr and excitations:
        # Statische respons per belastingsvector (testcase forma12d:
        # MACRO_ELAS_MULT met de opleggingen als CHAR_MECA_GLOBAL).
        proj += [
            "corr = MACRO_ELAS_MULT(",
            "    MODELE=model,",
            "    CHAM_MATER=fieldmat,",
        ] + (["    CARA_ELEM=cara,"] if cara_ref else []) + [
            "    NUME_DDL=numedd,",
            "    CHAR_MECA_GLOBAL=({0},),".format(", ".join(charge_names)),
            "    CAS_CHARGE=(",
        ] + ["        _F(NOM_CAS='c{0}', CHAR_MECA={1}),".format(k + 1, c)
             for k, (c, _f) in enumerate(excitations)] + [
            "    ),",
            ")",
            "",
        ]
        log("INFO: Statische correctie voor de weggelaten modes "
            "(MACRO_ELAS_MULT, CORR_STAT).")
    elif static_corr:
        static_corr = False
    projection_block = "\n".join(proj)
    if not excitations:
        # Zonder belasting een nulvector, zodat DYNA_VIBRA een EXCIT heeft.
        excitations = []
    dyna_block = _dyna_block(analysis, excitations, inc, n_out, alpha, beta,
                             zeta, static_corr, fields)
    if not excitations:
        dyna_block = dyna_block.replace("    EXCIT=(\n    ),\n", "")

    calc_champ_block = "#-- geen spanningen gevraagd (fields) --\n"
    if any(f in fields for f in _STRESS_GROUPS):
        calc_champ_block = (
            "# Spanningen op de fysieke respons (lineair in de "
            "verplaatsingen)\n" + _with_model(b2a._build_calc_champ_block(
                model_data, reactions=False), bool(cara_ref)))
    export_block = b2a._custom_code(model_data, "post") + (
        "# Resultaten in Mecway's eigen formaat (bridge: "
        "result_route='liml')@N"
        "aster_results_export.export_results(@N"
        "    res, {fragment}, {liml},@N"
        "    extra_results=dict((k, globals()[k]) for k in@N"
        "                       ('res_sup', 'res_moy', 'res_inf')@N"
        "                       if k in globals()),@N"
        "    precision='single', static=False, fields={fields!r})@N"
    ).replace("@N", chr(10)).format(
        fragment=b2a._py_str(b2a.results_file(work_dir)),
        liml=b2a._py_str(model_data.get("liml_path", "")),
        fields=[str(f) for f in fields])
    if CONFIG_BRIDGE.get("med_export", False):
        impr_resu_block = ("IMPR_RESU(FORMAT='MED', UNITE=80, "
                           "RESU=_F(RESULTAT=res, NOM_CHAM=('DEPL',)))\n")
    else:
        impr_resu_block = "#-- geen MED-export (instelling med_export uit) --\n"

    log("INFO: Respons op modale basis: {0} modes, {1} belastingsvector(en), "
        "Rayleigh alpha {2} 1/s, beta {3} s, extra modale demping {4}{5}."
        .format(b2a._mode_count(model_data, analysis), len(excitations),
                alpha, beta, zeta, ", statische correctie" if static_corr
                else ""))
    try:
        template = b2a._read_template(analysis["template"])
        b2a._write_comm(comm_path, template.format(
            defigroup_block=defigroup_block,
            model_block=model_block,
            cara_elem_block=cara_elem_block,
            materials_block=materials_block,
            fieldmat_block=fieldmat_block,
            load_funcs_block=load_funcs_block,
            loads_block=loads_block,
            assemblage_block=assemblage_block,
            calc_modes_block=calc_modes_block,
            projection_block=projection_block,
            dyna_block=dyna_block,
            calc_champ_block=calc_champ_block,
            export_block=export_block,
            impr_resu_block=impr_resu_block,
        ))
        log("SUCCESS: .comm geschreven naar '{0}'.".format(comm_path))
        export_path = os.path.join(work_dir, run_name(".export"))
        export_content = (
            "P actions make_etude\n"
            "P time_limit {time_limit}\n"
            "P memory_limit {memory_limit}\n"
            "P mpi_nbcpu 1\n"
            "P ncpus {ncpus}\n"
            "P no-mpi true\n"
            "F comm {comm} D 1\n"
            "F mail {msh} D 20\n"
            "F mess {mess} R 6"
        ).format(
            time_limit=analysis.get("time_limit", 1800),
            memory_limit=analysis.get("memory_limit", 2048),
            ncpus=CONFIG_BRIDGE["ncpus"],
            comm=comm_path, msh=msh_path,
            mess=os.path.join(work_dir, run_name(".mess")),
        ) + b2a._export_result_lines(work_dir, resu_name)
        with open(export_path, "w") as f:
            f.write(export_content)
        log("SUCCESS: .export geschreven naar '{0}'.".format(export_path))
        return export_path
    except KeyError as e:
        log(tr("ERROR: Template placeholder ontbreekt: {0}.",
               "ERROR: Template placeholder missing: {0}.").format(str(e)))
    except Exception as e:
        log(tr("ERROR: Genereren solver-bestanden mislukt: {0}",
               "ERROR: Generating solver files failed: {0}").format(str(e)))
    return None


def _with_model(block, cara):
    """CALC_CHAMP op het resultaat van REST_GENE_PHYS: MODELE, CHAM_MATER
    (en CARA_ELEM) expliciet, achter elke 'RESULTAT=res,'."""
    head = "res = CALC_CHAMP(\n    reuse=res,\n    RESULTAT=res,\n"
    add = "    MODELE=model,\n    CHAM_MATER=fieldmat,\n"
    if cara:
        add += "    CARA_ELEM=cara,\n"
    return block.replace(head, head + add)
