# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# bridge_modal_damping.py
# Verantwoordelijkheid: modale analyse met demping (complexe modes), de
# instelling modal_complex van modal_3d (Mecway V30). Hulpfuncties voor
# bridge_to_aster._build_modal:
#   plan(model_data, analysis, buckling) - None (reele modes, zoals Mecway)
#                                 of {"alpha", "beta", "dampers"}
#   assemblage_entry(plan)       - AMOR_MECA van de dempers in ASSEMBLAGE
#   damping_block(plan)          - de dempingsmatrix 'amor'
#   calc_modes_lines(...)        - MATR_AMOR en de zoekopties van CALC_MODES
#   zero_damping_entries(...)    - A_*-entries 0 voor de andere discreten
#
# Gedrag van Mecway (headless gemeten met Mecway 33, mw_solve.py, 7 okt
# 2026; docs/Modaal_en_knik.md par. 6): een Frequency-analyse (V30) negeert
# Rayleigh-demping (massa op een balk met alpha = 3,55 1/s: dezelfde
# frequenties 5,6494 / 8,4741 / 564,94 Hz als zonder, geen decay rate) en
# rekent met een Damper-element (line2) geen oplossing; Rayleigh hoort
# volgens de manual alleen bij Dynamic Response (manual v33 par. 9.26
# "Analysis types"). Complexe modes kent Mecway alleen in de akoestiek met
# impedantie (par. 6.9.2): <mode frequency=".. Hz" decayrate=".. 1/s"> met
# de frequentie = Im(s)/2pi en decay rate = -Re(s) van de eigenwaarde s
# (headless gezien: buis met Z = 2000 Pa.s/m, decayrate 70,35 1/s =
# (c/L) artanh(rho c/Z)). Dat formaat gebruikt de bridge ook hier.
#
# Code-Aster (catalogus calc_modes.py v17.4; testcases sdld21c, sdld23a,
# sdld27a): CALC_MODES TYPE_RESU='DYNAMIQUE' met MATR_AMOR (matr_asse_depl_r)
# lost het kwadratische eigenwaardeprobleem (K + s C + s^2 M) phi = 0 op;
# het resultaat is een MODE_MECA_C met complexe DEPL. Gemeten (massa 10 kg
# op een balk, k = 12 600 N/m, demper c = 213 N.s/m, zeta = 0,3):
#   FREQ = de GEDEMPTE frequentie Im(s)/2pi (5,38919 = 5,64944 sqrt(1 -
#   0,09) Hz) en AMOR_REDUIT = zeta = -Re(s)/|s| (0,29999), in de tabel
#   (RECU_TABLE NOM_PARA AMOR_REDUIT; niet in getAccessParameters).
#   Alleen geldig voor ondergedempte modes; overkritisch gedempte modes
#   (reele eigenwaarden) komen niet mee: alarm ALGELINE5_79 (minder modes
#   dan gevraagd; gezien met beta = 0,017 s op de axiale mode, zeta = 30).
#   SORENSEN met OPTION='PLUS_PETITE' kan niet (<EXCEPTION> ALGELINE2_71:
#   "la borne inferieure ... est nulle"); TRI_DIAG en QZ wel (gelijk tot
#   1e-9), SORENSEN met 'CENTRE' ook.
#
# Compatibel met IronPython 2.7.9 (alleen standaardbibliotheek).

import re

from bridge_utils import log, tr


def _dampers(model_data):
    """Elsets met Mecway Damper-materiaal (line2, DIS_T A_T_D_L)."""
    return sorted(e for e, s in model_data.get("discrete_sections",
                                                {}).items()
                  if s.get("kind") == "damper")


def _rayleigh(model_data):
    """(alpha [1/s], beta [s]) van Mecway's Rayleigh-demping (het hele
    model, bridge_from_liml), of (0, 0)."""
    ray = model_data.get("rayleigh") or {}
    return (float(ray.get("alpha", 0.0) or 0.0),
            float(ray.get("beta", 0.0) or 0.0))


def plan(model_data, analysis, buckling):
    """
    None: reele modes zonder demping (zoals Mecway, standaard). Anders een
    dict {"alpha", "beta", "dampers"} voor complexe modes: instelling
    modal_complex aan en demping in het model (Rayleigh en/of dempers).
    Zet dan model_data["damping_matrix"] (bridge_to_aster._build_cara_elem
    _block: dempingskarakteristiek 0 voor de andere discreten, anders
    alarm DISCRETS_28 per knoopmassa).
    """
    if buckling:
        return None
    alpha, beta = _rayleigh(model_data)
    dampers = _dampers(model_data)
    what = []
    if alpha or beta:
        what.append(tr("Rayleigh-demping (alpha {0} 1/s, beta {1} s)",
                       "Rayleigh damping (alpha {0} 1/s, beta {1} s)").format(
                           "%.6g" % alpha, "%.6g" % beta))
    if dampers:
        what.append(tr("dempers ({0})", "dampers ({0})").format(
            ", ".join(dampers)))
    if not analysis.get("modal_complex", False):
        if what:
            log(tr("INFO: Modaal: {0} telt niet mee: reele modes zonder "
                   "demping, zoals Mecway's Frequency-analyse. Voor de "
                   "gedempte frequentie en de dempingsverhouding per mode: "
                   "zet 'Complexe modes met demping' aan (pagina Modaal).",
                   "INFO: Modal: {0} is not included: real modes without "
                   "damping, like Mecway's Frequency analysis. For the "
                   "damped frequency and the damping ratio per mode: switch "
                   "on 'Complex modes with damping' (Modal page).").format(
                       tr(" en ", " and ").join(what)))
        return None
    if not what:
        log(tr("INFO: Complexe modes gevraagd, maar het model heeft geen "
               "demping (Rayleigh of dempers): reele modes.",
               "INFO: Complex modes requested, but the model has no "
               "damping (Rayleigh or dampers): real modes."))
        return None
    if alpha < 0.0 or beta < 0.0:
        log(tr("WARNING: Negatieve Rayleigh-coefficient (alpha {0}, beta "
               "{1}): negatieve demping, modes kunnen aangroeien.",
               "WARNING: Negative Rayleigh coefficient (alpha {0}, beta "
               "{1}): negative damping, modes may grow.").format(
                   "%.6g" % alpha, "%.6g" % beta))
    model_data["damping_matrix"] = True
    log(tr("INFO: Complexe modes (kwadratisch eigenwaardeprobleem, "
           "CALC_MODES met MATR_AMOR) met {0}: C = {1}.",
           "INFO: Complex modes (quadratic eigenvalue problem, CALC_MODES "
           "with MATR_AMOR) with {0}: C = {1}.").format(
               tr(" en ", " and ").join(what), _c_text(alpha, beta, dampers)))
    return {"alpha": alpha, "beta": beta, "dampers": dampers}


def _c_text(alpha, beta, dampers):
    terms = []
    if alpha:
        terms.append("alpha M")
    if beta:
        terms.append("beta K")
    if dampers:
        terms.append("C_dempers")
    return " + ".join(terms)


def assemblage_entry(plan_):
    """ASSEMBLAGE-entry voor de dempers: AMOR_MECA geeft de A_T_D_L van de
    DIS_T-elementen (testcase sdld21c). Zonder Rayleigh in DEFI_MATERIAU
    (bridge_to_aster._build_materials_block zet AMOR_ALPHA/BETA alleen bij
    DYNA_NON_LINE) dragen de andere elementen niets bij (gemeten: balk +
    demper, zeta van de mode zonder demper 1e-12)."""
    if not plan_ or not plan_["dampers"]:
        return ""
    return "        _F(MATRICE=CO('amor_dis'), OPTION='AMOR_MECA'),\n"


def damping_block(plan_):
    """
    De dempingsmatrix 'amor' na ASSEMBLAGE (en na de voorgespannen
    stijfheid): C = alpha M + beta K (Mecway manual v33 par. 9.26, het hele
    model, dus ook knoopmassa's en veren) + C_dempers. K is de stijfheid
    van het eigenwaardeprobleem ('rigidite'; met voorspanning K + Kg), zodat
    voor Rayleigh zeta_i = alpha / (2 w_i) + beta w_i / 2 met w_i de
    ongedempte eigenfrequentie, zoals bij de respons op modale basis
    (bridge_dyna_modal._dyna_block). COMB_MATR_ASSE met de
    Lagrange-vrijheidsgraden van de opleggingen werkt (gemeten: alpha M en
    beta K op een ingeklemde balk, zeta 0,30000 / 0,20000 / 0,0030000 =
    alpha / 2w per mode).
    """
    if not plan_:
        return ""
    terms = []
    if plan_["alpha"]:
        terms.append(("masse", plan_["alpha"]))
    if plan_["beta"]:
        terms.append(("rigidite", plan_["beta"]))
    if plan_["dampers"]:
        terms.append(("amor_dis", 1.0))
    head = ("# Dempingsmatrix voor de complexe modes (instelling "
            "modal_complex):\n# C = {0} (bridge_modal_damping)\n".format(
                _c_text(plan_["alpha"], plan_["beta"], plan_["dampers"])))
    if terms == [("amor_dis", 1.0)]:
        return head + "amor = amor_dis\n"
    return head + "amor = COMB_MATR_ASSE(COMB_R=(\n{0}))\n".format("".join(
        "    _F(MATR_ASSE={0}, COEF_R={1!r}),\n".format(m, float(c))
        for m, c in terms))


def calc_modes_lines(n, shift, free):
    """
    Regels voor CALC_MODES met demping (na MATR_MASS): MATR_AMOR, OPTION,
    CALC_FREQ en SOLVEUR_MODAL. shift: modal_shift_hz; free: geen
    opleggingen (de bridge zoekt dan rond -bridge_fshift,
    bridge_to_aster._FREE_SHIFT_BLOCK).
      shift 0 en opgelegd: PLUS_PETITE met TRI_DIAG (Lanczos; SORENSEN kan
        niet met PLUS_PETITE bij een kwadratisch probleem, ALGELINE2_71);
      anders: CENTRE met SORENSEN (testcase sdld21c).
    """
    lines = ["    MATR_AMOR=amor,"]
    if shift != 0.0:
        lines += ["    OPTION='CENTRE',",
                  "    CALC_FREQ=_F(FREQ={0!r}, NMAX_FREQ={1}),".format(
                      shift, n)]
    elif free:
        lines += ["    OPTION='CENTRE',",
                  "    CALC_FREQ=_F(FREQ=-bridge_fshift, NMAX_FREQ={0}),"
                  .format(n)]
    else:
        lines += ["    OPTION='PLUS_PETITE',",
                  "    CALC_FREQ=_F(NMAX_FREQ={0}),".format(n),
                  "    SOLVEUR_MODAL=_F(METHODE='TRI_DIAG'),"]
    return lines


# DISCRET-entry in bridge_to_aster._build_cara_elem_block (drie vormen,
# met of zonder REPERE).
_ENTRY_RE = re.compile(r"GROUP_MA='([^']+)',\s*CARA='([A-Z_]+)',\s*"
                       r"(?:REPERE='\w+',\s*)?VALE=\(([^)]*)\)")


def zero_damping_entries(entries):
    """
    Extra DISCRET-entries met demping 0 (A_* naast elke K_*) voor groepen
    zonder eigen A_*-entry. Zonder geeft AMOR_MECA per discreet element
    alarm DISCRETS_28 ("vous n'avez pas affecte les caracteristiques
    d'amortissements. Par defaut l'amortissement est nul", gezien met een
    knoopmassa). Catalogus affe_cara_elem.py v17.4: bij elke K_T_D_N,
    K_TR_D_N, K_T_N, K_T_D_L, ... hoort A_ met hetzelfde aantal waarden.
    """
    found = []
    damped = set()
    for text in entries:
        for group, cara, vals in _ENTRY_RE.findall(text):
            if cara.startswith("A_"):
                damped.add(group)
            elif cara.startswith("K_"):
                n = len([v for v in vals.split(",") if v.strip()])
                found.append((group, "A_" + cara[2:], n))
    out = []
    for group, cara, n in found:
        if group in damped:
            continue
        damped.add(group)
        out.append("        _F(GROUP_MA='{0}', CARA='{1}', VALE=({2},)),"
                   .format(group, cara, ", ".join(["0.0"] * n)))
    return out
