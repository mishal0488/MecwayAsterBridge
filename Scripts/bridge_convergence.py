# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# bridge_convergence.py
# Verantwoordelijkheid: lees de convergentie-informatie van STAT_NON_LINE /
# DYNA_NON_LINE uit de uitvoer van Code-Aster (fort.6 tijdens de run,
# study.mess erna). Gebruikt door het convergentiescherm
# (aster_convergence_monitor.py); zelf zonder GUI en zonder bestanden te
# schrijven.
#
# Puur Python, alleen de standaardbibliotheek; draait in Python 2.7
# (IronPython 2.7.9, Mecway) en Python 3 (de CPython 3.11 van Code-Aster).
#
# Wat Code-Aster 17.4 per rekenstap schrijft (gezien in study.mess van de
# testmodellen, o.a. statisch en impliciet dynamisch met contact):
#
#   Instant de calcul:  1.250000000000e-03 - Niveau de découpe: 2
#   |     CONTACT    |     NEWTON     |     RESIDU     | ...   (kop, 3 regels)
#   |     1        X |     0        X | 3.37162E-10    | ...   (een regel per
#                                                               Newton-iteratie)
#   Critère(s) de convergence atteint(s)        of  <Erreur> ... / <Action> ...
#   Le résidu de type <RESI_GLOB_RELA> vaut ... au noeud et degré de liberté <298 DY>
#   Temps CPU consommé dans ce pas de temps    : 2.000 s
#   * Nombre d'itérations de Newton            : 32
#   * Nombre de liaisons de contact            : 8
#   [  0%] Instant calculé : 1.25000e-03, ...
#
# De kolommen verschillen per model (zonder contact: 4 kolommen; discreet
# contact, frottement, RESI_REFE_RELA ...). De kop wordt daarom per tabel
# gelezen en via _COLUMN_KEYS op vaste sleutels gezet; onbekende kolommen
# komen onder hun kopnaam in it["extra"]. Een "X" achter een waarde betekent
# dat dat criterium nog niet gehaald is, "E" dat de iteratie in fout eindigt
# (NIET geverifieerd in de U-doc; afgeleid uit de uitvoer: "E" staat op de
# laatste iteratie voor "<Erreur> Le nombre maximum d'itérations ...").

import io
import os
import re

try:
    from bridge_utils import tr
except ImportError:      # pragma: no cover - bridge_utils naast dit bestand
    def tr(nl, en):
        return nl


# ==========================================
# HERKENNING VAN REGELS
# ==========================================
_INSTANT_RE = re.compile(
    r"Instant de calcul:\s*([-+0-9.eEdD]+)\s*-\s*Niveau de d.coupe:\s*(\d+)")
_PROGRESS_RE = re.compile(
    r"^\[\s*(\d+)%\]\s*Instant calcul.\s*:\s*([-+0-9.eE]+)")
_CPU_RE = re.compile(r"Temps CPU consomm. dans ce pas de temps\s*:\s*([0-9.]+)")
_SUMMARY_RE = re.compile(r"^\*\s*(.+?)\s*:\s*([-+0-9.eE]+)")
_CONTACT_INIT_RE = re.compile(
    r"Il y a (\d+) points initialement en contact et (\d+) points exclus")
_RESI_LOC_RE = re.compile(
    r"Le r.sidu de type <(\w+)> vaut\s*([-+0-9.eE]+) au noeud et "
    r"degr. de libert. <\s*(\S+)\s+(\S+?)\s*>")
_PIVOT_RE = re.compile(
    r"associ\S* au noeud (\S+) et \S+ la composante (\S+?)\.?\s*$")
_BOX_CODE_RE = re.compile(
    r"<(EXCEPTION|F|S|E|A|[A-Za-z]+Error)>\s*<([A-Za-z0-9_]+)>")
_COMMAND_RE = re.compile(r"\b(STAT_NON_LINE|DYNA_NON_LINE|THER_NON_LINE)\s*\(")
_RESI_CRIT_RE = re.compile(r"\bRESI_GLOB_RELA\s*=\s*([-+0-9.eE]+)")
_ITER_MAX_RE = re.compile(r"\bITER_GLOB_MAXI\s*=\s*(\d+)")
_PARTIAL_RE = re.compile(r"BRIDGE_ONVOLLEDIG\s+(\S+)")
# Bezwijken van elementen (bridge_failure.py, print in de .comm)
_FAILURE_RE = re.compile(r"BRIDGE_BEZWIJKEN t=([-+0-9.eE]+) elset=(\S+) "
                         r"nieuw=(\d+) totaal=(\d+)")
# Eindtijd van de bezwijklus (bridge_failure.py): die rekent met een
# solve-aanroep per tijdstap, en '[ n%]' geldt dan per aanroep.
_TEND_RE = re.compile(r"^BRIDGE_TEND\s+([-+0-9.eE]+)\s*$")
# Expliciet schema in de commando-echo (geen Newton-iteraties).
_EXPLICIT_RE = re.compile(r"SCHEMA\s*=\s*'(DIFF_CENT|TCHAMWA)'")
_DIAG_MARK = "--- DIAGNOSTIC JOB :"
_BOX_TOP = u"╔"     # bovenrand meldingskader
_BOX_BOTTOM = u"╚"  # onderrand meldingskader
_BOX_SIDE = u"║"    # zijkant meldingskader
_CONVERGED = u"Critère(s) de convergence atteint(s)"
_NUM_CELL_RE = re.compile(r"^\s*([-+]?\d[\d.]*(?:[EeDd][-+]?\d+)?)\s*([A-Z]?)\s*$")

# Kolomkop (samengevoegd over de drie kopregels) -> sleutel. Volgorde telt:
# de eerste passende regel wint.
_COLUMN_KEYS = [
    ("BCL. GEOM.", "geom"),
    ("BCL. CONT.", "cont_loop"),
    ("BCL. FROT.", "frot_loop"),
    ("NEWTON ITERATION", "newton"),
    ("RESI_GLOB_RELA", "resi_rela"),
    ("RESI_GLOB_MAXI", "resi_maxi"),
    ("RESI_REFE_RELA", "resi_refe"),
    ("RESI_GEOM", "resi_geom"),
    ("OPTION", "option"),
    ("VARI. CONT.", "cont_vari"),
    ("CRIT. GEOM.", "cont_crit_geom"),
    ("PRESSURE", "cont_press"),
    ("CRITERE VALEUR", "cont_crit"),
    ("PENETRATION", "cont_pen"),
    ("DISCRET", "cont_discrete_iter"),
]

# Samenvatting na elke stap ("* Nombre de ... : n") -> sleutel.
_SUMMARY_KEYS = [
    (u"Nombre d'it", "newton_count"),
    (u"Nombre de liaisons de contact", "contact_links"),
    (u"contact/pas contact", "cycles_contact"),
    (u"glissement/adh", "cycles_slip_stick"),
    (u"glissement avant/arri", "cycles_slip_dir"),
    (u"point fixe contact", "cycles_fixed_point"),
]

# Gebeurtenissen in de convergentietabel (Frans in study.mess) met uitleg
# (nl, en). Teksten gezien in de study.mess-bestanden van de testmodellen.
EVENT_TEXTS = [
    (u"nombre maximum d'it", u"Maximaal aantal Newton-iteraties "
     u"(ITER_GLOB_MAXI) bereikt zonder convergentie.",
     u"Maximum number of Newton iterations (ITER_GLOB_MAXI) reached "
     u"without convergence."),
    (u"matrice du syst", u"De stijfheidsmatrix is singulier (pivot nul); "
     u"zie de pivot-melding (knoop en component).",
     u"The stiffness matrix is singular (zero pivot); see the pivot "
     u"message (node and component)."),
    (u"essaie de découper", u"Code-Aster probeert de tijdstap op te "
     u"delen.", u"Code-Aster tries to subdivide the time step."),
    (u"On découpe le pas", u"Tijdstap opgedeeld.",
     u"Time step subdivided."),
    (u"chec dans la tentative de découper", u"Opdelen lukt niet meer "
     u"(maximaal aantal opdelingsniveaus bereikt).",
     u"Subdividing no longer possible (maximum number of subdivision "
     u"levels reached)."),
    (u"On arrête le calcul", u"Code-Aster stopt de berekening.",
     u"Code-Aster stops the calculation."),
]

# Uitleg bij een pivot nul, per soort vrijheidsgraad. De raad volgt de
# tekst die Code-Aster zelf bij de pivot-melding geeft (study.mess:
# "Si vous faites du contact, il ne faut pas que la structure ne 'tienne'
# que par le contact" en "exclure certains noeuds des conditions de contact
# (mots clés SANS_NOEUD et SANS_GROUP_NO)").
_PIVOT_HINT_CONTACT = (u"Contact-Lagrange (LAGS_C/LAGS_F): de contactvoorwaarde "
                       u"is overbepaald (knoop zit ook in een oplegging of in "
                       u"twee contactparen) of het deel hangt alleen aan "
                       u"contact. Sluit zulke knopen uit (SANS_GROUP_NO) of "
                       u"geef het deel een oplegging.",
                       u"Contact Lagrange (LAGS_C/LAGS_F): the contact "
                       u"condition is overdetermined (the node is also in a "
                       u"support or in two contact pairs) or the part hangs "
                       u"on contact only. Exclude such nodes (SANS_GROUP_NO) "
                       u"or give the part a support.")
_PIVOT_HINT_PHYS = (u"Fysieke vrijheidsgraad: waarschijnlijk een niet "
                    u"vastgehouden starre-lichaamsbeweging (te weinig "
                    u"opleggingen, of een deel dat alleen via contact op "
                    u"zijn plaats blijft).",
                    u"Physical degree of freedom: probably an unrestrained "
                    u"rigid-body motion (too few supports, or a part held "
                    u"in place only by contact).")
_PIVOT_HINT_LAGR = (u"Lagrange-vrijheidsgraad: waarschijnlijk een dubbele of "
                    u"tegenstrijdige randvoorwaarde.",
                    u"Lagrange degree of freedom: probably a duplicated or "
                    u"conflicting constraint.")


def _to_float(text):
    try:
        return float(text.replace("D", "E").replace("d", "e"))
    except (ValueError, AttributeError):
        return None


# ==========================================
# LINEAIRE OPLOSSER (MECA_STATIQUE met INFO=2)
#
# Wat Code-Aster 17.4 schrijft binnen "res = MECA_STATIQUE(...)" tot
# "# Fin commande" (gezien in de .mess van een klein 3D-testmodel en een
# blok van 37 389 vergelijkingen; MUMPS 5.6.2 met ICNTL(4)=2):
#
#   Le système linéaire à résoudre a 252 degrés de liberté:
#      - 228 sont des degrés de liberté physiques (... 76 noeuds ...)
#      - 24 sont les couples de paramètres de Lagrange ... 12 relations
#   La matrice est de taille 252 équations.
#     Elle contient 4182 termes non nuls si elle est symétrique ...
#     Soit un taux de remplissage de  12.774 %.
#   TEMPS CPU/SYS/ELAPSED ASSEMBLAGE M        :  0.1000D-02 0.0000D+00 ...
#   Entering DMUMPS 5.6.2 with JOB, N, NNZ =   1 ...   (1 analyse, 2
#   ... Elapsed time in analysis driver=  0.0010      factorisatie,
#   ... INFOG (12) Number of negative pivots = 24    3 oplossen)
#   ERROR ANALYSIS BEFORE ITERATIVE REFINEMENT / STATISTICS AFTER ...
#      RINFOG(6):SCALED RESIDUAL, (7)/(8) backward errors W1/W2,
#      (9) Upper bound ERROR, (10)/(11) CONDITION NUMBER (1)/(2)
#   Champ stocké <DEPL> ...                            (oplossing klaar)
#   Statistics:  - _computeMatrix called 1 times in 0.002843 s  ...
#   # Fin commande #0011   user+syst: 0.06s (syst: 0.00s, elaps: 0.00s)
#
# GCPC (gcpc.F90): "* GCPC NORME DU RESIDU = ...", per iteratie
# "* ITERATION n NORME DU RESIDU EN ABS/RELA = a r" (bij INFO=2 alleen als
# het residu met een factor 10 gedaald is: paraaf = 0.1d0), daarna
# "* NORME DU RESIDU INITIAL/FINAL/RELATIF=" en "* CONVERGENCE EN n
# ITERATIONS".
#
# Betekenis (U2.08.03 par. 7.2.1): de foutgrens RINFOG(9) is ongeveer
# conditiegetal x 'backward error' en wordt door Code-Aster met
# SOLVEUR/RESI_RELA vergeleken (U4.50.01 par. 3.5.3).
# ==========================================
_LIN_CMD_RE = re.compile(r"^\s*(?:[A-Za-z_]\w*\s*=\s*)?MECA_STATIQUE\s*\(")
# Modaal/knik (bridge_to_aster._build_modal): de echo van CALC_MODES en de
# regels die aster_results_export.print_modes na de berekening print:
#   BRIDGE_MODE 1 FREQ=7.0737 MASS_EFFE_UN_DX=.. MASS_EFFE_UN_DY=.. ..
#   BRIDGE_MODE 1 CHAR_CRIT=41124.7
_MODES_CMD_RE = re.compile(r"^\s*(?:[A-Za-z_]\w*\s*=\s*)?CALC_MODES\s*\(")
_MODE_LINE_RE = re.compile(r"^BRIDGE_MODE\s+(\d+)\s+(.*)$")
_LIN_END = u"# Fin commande"
_LIN_TXT_KW_RE = re.compile(r"\b(METHODE|PRE_COND|POSTTRAITEMENTS)\s*=\s*"
                            r"'(\w+)'")
_LIN_NUM_KW_RE = re.compile(r"\b(RESI_RELA|NPREC|NMAX_ITER|NIVE_REMPLISSAGE)"
                            r"\s*=\s*([-+0-9.eE]+)")
_LIN_DOF_RE = re.compile(u"système linéaire à résoudre a "
                         r"(\d+) degr")
_LIN_PHYS_RE = re.compile(r"^-\s*(\d+) sont des degr")
_LIN_NODES_RE = re.compile(r"port\S* par (\d+) noeuds")
_LIN_LAGR_RE = re.compile(r"^-\s*(\d+) sont les couples")
_LIN_REL_RE = re.compile(r"aux (\d+) relations lin")
_LIN_EQ_RE = re.compile(r"La matrice est de taille (\d+)")
_LIN_NNZ_RE = re.compile(r"contient (\d+) termes non nuls")
_LIN_FILL_RE = re.compile(r"taux de remplissage de\s*([0-9.]+)\s*%")
_LIN_ASSEMB_RE = re.compile(r"TEMPS CPU/SYS/ELAPSED ASSEMBLAGE ([MV])\s*:"
                            r"\s*(\S+)\s+(\S+)\s+(\S+)")
_LIN_JOB_RE = re.compile(r"Entering [SDCZ]MUMPS (\S+) with JOB\S*\s*"
                         r"(?:,\s*N,\s*NNZ\s*)?=\s*(-?\d+)")
_LIN_DRIVER_RE = re.compile(r"Elapsed time in (analysis|factorization|solve)"
                            r" driver\s*=\s*([0-9.eEdD+-]+)")
_LIN_ORDER_RE = re.compile(r"Ordering based on (\S+)")
_LIN_FLOPS_RE = re.compile(r"RINFOG\(1\)\)?\s*=\s*([-+0-9.eEdD]+)")
_LIN_INFOG_RE = re.compile(r"(?<![RA-Z])INFOG\s*\(\s*(\d+)\s*\)[^=]*="
                           r"\s*([-+0-9.eEdD]+)")
_LIN_MEM_RE = re.compile(r"\(INFOG\((17|19|22)\)\):\s*(\d+)")
# Schattingen uit de analysefase (al bekend tijdens de factorisatie)
_LIN_FLOPS_EST_RE = re.compile(r"RINFOG\(1\) Operations during elimination "
                               r"\(estim\)\s*=\s*(\S+)")
_LIN_ENTRIES_EST_RE = re.compile(r"\(20\) Number of entries in factors "
                                 r"\(estim\.\)\s*=\s*(\d+)")
_LIN_REFINE_STEPS_RE = re.compile(r"NUMBER OF STEPS OF ITERATIVE REFINEMENTS"
                                  r" PERFORMED\s+(\d+)")
_LIN_RESID_RE = re.compile(r"^RESIDUAL IS \.+ \(MAX-NORM\)\s*=\s*(\S+)")
_LIN_RINFOG_RE = re.compile(r"^(?:RINFOG|-+)\(\s*(\d+)\)\s*:.*=\s*(\S+)\s*$")
_LIN_GCPC_INIT_RE = re.compile(r"^\*\s*GCPC\s+NORME DU RESIDU =\s*(\S+)")
_LIN_GCPC_TARGET_RE = re.compile(r"NORME DU RESIDU A ATTEINDRE EN ABS/RELA="
                                 r"\s*(\S+)\s+(\S+)")
_LIN_GCPC_ITER_RE = re.compile(r"^\*\s*ITERATION\s+(\d+) NORME DU RESIDU EN "
                               r"ABS/RELA =\s*(\S+)\s+(\S+)")
_LIN_GCPC_FINAL_RE = re.compile(r"NORME DU RESIDU INITIAL/FINAL/RELATIF=\s*"
                                r"(\S+)\s+(\S+)\s+(\S+)")
_LIN_GCPC_CONV_RE = re.compile(r"^\*\s*CONVERGENCE EN\s+(\d+) ITERATIONS")
_LIN_STAT_GROUP_RE = re.compile(r"^- (\w+)\.\s*$")
_LIN_STAT_RE = re.compile(r"^(\s*)- (\w+) called (\d+) times in "
                          r"([-+0-9.eE]+) s")
_LIN_FIN_RE = re.compile(r"user\+syst:\s*([0-9.]+)s.*elaps:\s*([0-9.]+)s")

# Fasen van een lineaire oplossing, in volgorde (sleutel, nl, en).
LIN_PHASES = [
    ("matrix", u"Stijfheidsmatrix (elementen + assemblage)",
     u"Stiffness matrix (elements + assembly)"),
    ("analysis", u"Analyse (ordening, symbolisch)",
     u"Analysis (ordering, symbolic)"),
    ("factor", u"Factorisatie", u"Factorization"),
    ("rhs", u"Belastingsvector", u"Load vector"),
    ("solve", u"Oplossen", u"Solve"),
    ("post", u"Spanningen + opslag", u"Stresses + storage"),
]
_LIN_PHASE_KEYS = [p[0] for p in LIN_PHASES]


def lin_phase_title(key, method="MUMPS"):
    """Titel van een fase in de gekozen taal; bij GCPC is 'factor' de
    voorconditionering (onvolledige factorisatie)."""
    if key == "factor" and method == "GCPC":
        return tr(u"Voorconditionering (LDLT_INC)",
                  u"Preconditioner (LDLT_INC)")
    if key == "solve" and method == "GCPC":
        return tr(u"Oplossen (GCPC-iteraties)", u"Solve (GCPC iterations)")
    for k, nl, en in LIN_PHASES:
        if k == key:
            return tr(nl, en)
    return key


class LinearSolve(object):
    """Een MECA_STATIQUE-aanroep: matrix, oplosser, fasen en kwaliteit."""

    def __init__(self, seen_at=None):
        self.command = "MECA_STATIQUE"
        self.started_at = seen_at     # wandkloktijd (live), anders None
        self.status = "bezig"         # bezig / klaar / mislukt
        self.solver = {}              # METHODE, PRE_COND, RESI_RELA, NPREC ..
        self.dofs = None              # vrijheidsgraden totaal
        self.dofs_phys = None
        self.nodes = None
        self.lagr_pairs = None        # Lagrange-parameters ("couples de
        #                               paramètres de Lagrange", 2 per relatie)
        self.relations = None         # gedualiseerde lineaire relaties
        self.equations = None
        self.nnz = None               # niet-nul-termen (symmetrisch)
        self.fill_pct = None          # vullingsgraad [%]
        self.mumps = {}               # version, ordering, flops, entries,
        #                               neg_pivots, delayed, front, mem_*
        self.refine_before = {}       # residu en kwaliteit voor verfijning
        self.refine_after = {}        # ... na verfijning
        self.refine_steps = None
        self.gcpc_init = None         # beginresidu (absoluut)
        self.gcpc_target = None       # (absoluut, relatief)
        self.gcpc_iters = []          # (iteratie, absoluut, relatief)
        self.gcpc_final = None        # (begin, eind, relatief)
        self.gcpc_converged = None    # aantal iteraties bij convergentie
        self.phases = {}              # sleutel -> {start, end, seconds}
        self.current = None           # lopende fase
        self.stats = {}               # "MumpsSolver.factorize" -> s
        self.cpu_total = None
        self.elapsed_total = None
        self.detailed = False         # INFO=2-uitvoer gezien (ASSEMBLAGE ..)

    # ---------- fasen ----------
    def method(self):
        return self.solver.get("METHODE", u"MUMPS")

    def phase_keys(self):
        """De fasen die bij deze oplosser horen (volgorde van
        LIN_PHASES)."""
        keys = list(_LIN_PHASE_KEYS)
        if self.method() != "MUMPS":
            keys.remove("analysis")
        return keys

    def phase_seconds(self, key, now=None):
        """Duur van een fase [s]: de waarde uit de uitvoer, anders de
        wandklok (start..eind, of start..now als de fase loopt); None als
        onbekend."""
        ph = self.phases.get(key)
        if not ph:
            return None
        if ph.get("seconds") is not None:
            return ph["seconds"]
        if ph.get("start") is not None:
            end = ph.get("end")
            if end is None and key == self.current and now is not None:
                end = now
            if end is not None:
                return max(0.0, end - ph["start"])
        return None

    def finished(self):
        return self.status != "bezig"

    def forward_error(self):
        """Geschatte relatieve fout op de oplossing (MUMPS RINFOG(9) na
        verfijning), of None."""
        return self.refine_after.get(9)

    def condition(self):
        return self.refine_after.get(10)


def _e(value, digits):
    """Getal in e-notatie met een punt (%-opmaak: in IronPython geeft
    str.format de decimale komma van de landinstelling); '-' bij None."""
    if value is None:
        return u"-"
    return u"%.*e" % (digits, value)


def assess_linear(lin):
    """
    Beoordeling van een lineaire oplossing: lijst (niveau, tekst) met
    niveau "ok", "info", "let_op" of "fout", in de gekozen taal.

    Bronnen: U4.50.01 par. 3.5.3 (RESI_RELA: foutgrens van MUMPS; groter
    = Code-Aster stopt), U2.08.03 par. 7.2.1 (conditiegetal 1e4 = goed
    geconditioneerd, tot 1e20 = zeer slecht; 'backward error' bij een
    goede oplossing in de orde van de machineprecisie 1e-15) en par. 7.1
    (conditiegetal > 1e12 verraadt een singulariteit, gemeten zonder
    voorbehandeling). De grens 1e-12 voor de 'backward error' is een
    keuze van de bridge (NIET geverifieerd in de U-doc).
    """
    out = []
    resi = lin.solver.get("RESI_RELA")
    method = lin.method()
    nprec = lin.solver.get("NPREC")
    if method == "MUMPS":
        if nprec is not None and nprec < 0:
            if resi is not None and resi < 0:
                out.append(("fout", tr(
                    u"NPREC < 0 en RESI_RELA < 0: geen enkele controle op de "
                    u"oplossing (singuliere matrix wordt niet gemeld).",
                    u"NPREC < 0 and RESI_RELA < 0: no check of the solution "
                    u"at all (a singular matrix is not reported).")))
            else:
                out.append(("info", tr(
                    u"NPREC = {0}: detectie van een singuliere matrix staat "
                    u"uit; de foutgrens (RESI_RELA) is de controle.",
                    u"NPREC = {0}: detection of a singular matrix is off; "
                    u"the error bound (RESI_RELA) is the check.").format(
                        int(nprec))))
        err = lin.forward_error()
        if err is not None and resi is not None and resi > 0:
            if err > resi:
                out.append(("fout", tr(
                    u"Geschatte fout op de oplossing {0} > RESI_RELA "
                    u"{1}: de oplossing is onnauwkeurig (Code-Aster "
                    u"stopt). Controleer opleggingen en materiaal.",
                    u"Estimated error of the solution {0} > RESI_RELA "
                    u"{1}: the solution is inaccurate (Code-Aster "
                    u"stops). Check supports and material.").format(
                        _e(err, 2), _e(resi, 0))))
            else:
                out.append(("ok", tr(
                    u"Geschatte fout op de oplossing {0} <= RESI_RELA "
                    u"{1}.",
                    u"Estimated error of the solution {0} <= RESI_RELA "
                    u"{1}.").format(_e(err, 2), _e(resi, 0))))
        cond = lin.condition()
        if cond is not None:
            if cond > 1e12:
                out.append(("let_op", tr(
                    u"Conditiegetal {0} > 1e12: wijst op een (bijna) "
                    u"singuliere matrix, bijv. een deel dat niet of "
                    u"nauwelijks is vastgehouden, of zeer grote "
                    u"stijfheidsverschillen.",
                    u"Condition number {0} > 1e12: indicates a "
                    u"(nearly) singular matrix, e.g. a part that is not or "
                    u"barely held, or very large stiffness "
                    u"differences.").format(_e(cond, 1))))
            elif cond > 1e4:
                out.append(("info", tr(
                    u"Conditiegetal {0}: gebruikelijk voor een "
                    u"FE-model (1e4 = goed, 1e12 en hoger = verdacht).",
                    u"Condition number {0}: usual for an FE model (1e4 "
                    u"= good, 1e12 and above = suspicious).").format(
                        _e(cond, 1))))
            else:
                out.append(("ok", tr(
                    u"Conditiegetal {0}: goed geconditioneerd.",
                    u"Condition number {0}: well conditioned.").format(
                        _e(cond, 1))))
        bw = lin.refine_after.get(7)
        if bw is not None and bw > 1e-12:
            out.append(("let_op", tr(
                u"'Backward error' {0} is groot (normaal ~1e-15): de "
                u"factorisatie had numerieke moeite.",
                u"Backward error {0} is large (normally ~1e-15): the "
                u"factorization had numerical difficulty.").format(
                    _e(bw, 1))))
        # Negatieve pivots: een stelsel met Lagrange-multiplicatoren is
        # niet positief definiet; in de testruns was het aantal negatieve
        # pivots gelijk aan het aantal Lagrange-parameters (252 en 268 569
        # vergelijkingen); met een mechanisme (alleen DX vastgehouden,
        # FACTOR_57) 10 bij 8. Meer negatieve pivots wijst op een niet positief
        # definiete stijfheid (traagheidsstelling van Sylvester; NIET
        # geverifieerd in de U-doc).
        neg = lin.mumps.get("neg_pivots")
        if (neg is not None and lin.lagr_pairs is not None and
                neg > lin.lagr_pairs):
            out.append(("let_op", tr(
                u"{0} negatieve pivots, meer dan de {1} Lagrange-parameters: "
                u"de stijfheidsmatrix is niet positief definiet (negatieve "
                u"stijfheid, of een mechanisme dat alleen numeriek "
                u"vastzit).",
                u"{0} negative pivots, more than the {1} Lagrange "
                u"parameters: the stiffness matrix is not positive definite "
                u"(negative stiffness, or a mechanism held only "
                u"numerically).").format(int(neg), lin.lagr_pairs)))
        delayed = lin.mumps.get("delayed")
        if delayed:
            out.append(("info", tr(
                u"{0} uitgestelde pivots (numeriek pivoteren nodig).",
                u"{0} delayed pivots (numerical pivoting needed).").format(
                    int(delayed))))
    elif method == "GCPC":
        rel = lin.gcpc_final[2] if lin.gcpc_final else None
        if lin.gcpc_converged is not None:
            out.append(("ok", tr(
                u"GCPC geconvergeerd in {0} iteraties (residu {1}, "
                u"criterium {2}).",
                u"GCPC converged in {0} iterations (residual {1}, "
                u"criterion {2}).").format(
                    lin.gcpc_converged, _e(rel, 2), _e(resi, 0))))
            # Grens 20 % van de vergelijkingen: keuze van de bridge (NIET
            # geverifieerd in de U-doc; NMAX_ITER 0 = 50 %).
            if lin.equations and lin.gcpc_converged > 0.2 * lin.equations:
                out.append(("let_op", tr(
                    u"Veel iteraties ({0} bij {1} vergelijkingen): trage "
                    u"convergentie (slecht geconditioneerd); MUMPS is dan "
                    u"meestal sneller.",
                    u"Many iterations ({0} for {1} equations): slow "
                    u"convergence (poorly conditioned); MUMPS is then "
                    u"usually faster.").format(lin.gcpc_converged,
                                               lin.equations)))
        elif lin.status == "mislukt":
            out.append(("fout", tr(
                u"GCPC niet geconvergeerd (maximaal aantal iteraties, "
                u"NMAX_ITER = 0: vergelijkingen / 2). Gebruik MUMPS.",
                u"GCPC did not converge (maximum number of iterations, "
                u"NMAX_ITER = 0: equations / 2). Use MUMPS.")))
    if lin.status == "mislukt" and not [o for o in out if o[0] == "fout"]:
        out.append(("fout", tr(
            u"MECA_STATIQUE is mislukt; zie de meldingen hieronder.",
            u"MECA_STATIQUE failed; see the messages below.")))
    return out


def _column_key(name):
    for pattern, key in _COLUMN_KEYS:
        if pattern in name:
            return key
    return None


def event_hint(text):
    """Uitleg (gekozen taal) bij een <Erreur>/<Action>-tekst, of u''."""
    for pattern, nl, en in EVENT_TEXTS:
        if pattern in text:
            return tr(nl, en)
    return u""


def pivot_hint(component):
    """Uitleg bij een pivot nul op deze component (bijv. 'LAGS_C', 'DZ')."""
    comp = component.upper()
    if comp.startswith("LAGS_") or comp.startswith("LAG_C"):
        return tr(*_PIVOT_HINT_CONTACT)
    if comp.startswith("LAGR") or comp.startswith("LAG"):
        return tr(*_PIVOT_HINT_LAGR)
    return tr(*_PIVOT_HINT_PHYS)


# ==========================================
# DATASTRUCTUREN
# ==========================================
class Step(object):
    """Een rekenpoging ("Instant de calcul"); een opgedeelde stap geeft
    meerdere pogingen met een hoger niveau."""

    def __init__(self, index, instant, level, seen_at=None):
        self.index = index            # volgnummer van de poging (0, 1, ...)
        self.instant = instant        # doeltijdstip van de poging
        self.level = level            # Niveau de découpe (1 = niet opgedeeld)
        self.seen_at = seen_at        # wandkloktijd (live), anders None
        self.iterations = []          # dicts per Newton-iteratie
        self.status = "bezig"         # bezig / geconvergeerd / mislukt
        self.error = u""              # <Erreur>-tekst (Frans)
        self.actions = []             # <Action>-teksten (Frans)
        self.cpu = None               # CPU-tijd van de poging [s]
        self.summary = {}             # newton_count, contact_links, cycles_*
        self.resi_location = u""      # bijv. "knoop 298 DY"
        self.contact_initial = None   # (in contact, uitgesloten)
        self.energy = {}              # TOTAL-regel van BILAN D'ENERGIE

    def newton_count(self):
        """Aantal Newton-iteraties (samenvatting, anders tabelregels)."""
        n = self.summary.get("newton_count")
        if n is not None:
            return int(n)
        return len(self.iterations)

    def geom_loops(self):
        """Aantal geometrische contactlussen (BCL. GEOM.), of None."""
        vals = [it["geom"] for it in self.iterations
                if it.get("geom") is not None]
        return int(max(vals)) if vals else None

    def status_changes(self):
        """Som van de contactstatuswisselingen (VARI. CONT.), of None."""
        vals = [it["cont_vari"] for it in self.iterations
                if it.get("cont_vari") is not None]
        return int(sum(vals)) if vals else None


class ConvergenceData(object):
    """Alles wat tot nu toe uit de uitvoer gelezen is."""

    def __init__(self):
        self.steps = []
        self.messages = []        # dicts: kind, code, text, hint, step
        self.progress = []        # (procent, instant)
        self.command = None       # STAT_NON_LINE / DYNA_NON_LINE
        self.resi_crit = None     # RESI_GLOB_RELA uit het commando
        self.iter_max = None      # ITER_GLOB_MAXI uit het commando
        self.diagnostic = u""     # DIAGNOSTIC JOB, leeg zolang bezig
        self.partial_until = None  # BRIDGE_ONVOLLEDIG-tijdstip
        self.t_end = None          # BRIDGE_TEND (bezwijklus), anders None
        self.explicit = False      # expliciet schema gezien
        self.fin = False           # "<I> <FIN>" gezien (study.comm tot FIN)
        self.linear = []           # LinearSolve per MECA_STATIQUE
        self.modes = []            # (nr, {FREQ/CHAR_CRIT/MASS_EFFE_UN_..})
        self.lines = 0

    def finished(self):
        return bool(self.diagnostic)

    def last_linear(self):
        return self.linear[-1] if self.linear else None

    def overall_progress(self):
        """(procent, tijdstip) van de hele run, of None. Met BRIDGE_TEND
        (bezwijklus, een solve per tijdstap) uit tijdstip / eindtijd; anders
        de laatste '[ n%]' van Code-Aster."""
        if not self.progress:
            return None
        pct, inst = self.progress[-1]
        if self.t_end and inst is not None and self.t_end > 0.0:
            pct = int(min(100.0, max(0.0, 100.0 * inst / self.t_end)))
        return pct, inst

    def last_converged(self):
        for step in reversed(self.steps):
            if step.status == "geconvergeerd":
                return step
        return None


# ==========================================
# PARSER (incrementeel)
# ==========================================
class ConvergenceParser(object):
    """
    Leest de uitvoer regel voor regel; feed() mag met stukken tekst worden
    aangeroepen die midden in een regel eindigen (live fort.6).

    hints : dict code -> uitleg in de gekozen taal (bijv.
            bridge_run_aster.mess_hints()); optioneel.
    clock : functie die de wandkloktijd geeft (live), of None.
    """

    def __init__(self, hints=None, clock=None):
        self.data = ConvergenceData()
        self._hints = hints or {}
        self._clock = clock
        self._partial = u""
        self._columns = None       # sleutels van de huidige tabel
        self._names = None         # kopnamen van de huidige tabel
        self._header = []          # kopregels in opbouw
        self._energy_cols = None
        self._box = None           # [kind, code, regels] in opbouw
        self._pending = u""        # regel die op de volgende doorloopt
        self._msg_seen = {}
        self._lin = None           # LinearSolve binnen MECA_STATIQUE
        self._lin_refine = None    # "before" / "after" (MUMPS-verfijning)
        self._lin_group = u""      # groep in het Statistics-blok
        self._lin_prev = u""       # vorige regel binnen MECA_STATIQUE

    # ---------- invoer ----------
    def feed(self, text):
        """Verwerk een stuk tekst (unicode); een onvolledige laatste regel
        wordt bewaard tot het volgende stuk."""
        text = self._partial + text
        lines = text.split(u"\n")
        self._partial = lines.pop()
        for line in lines:
            self.feed_line(line.rstrip(u"\r"))

    def close(self):
        """Verwerk een eventuele onvolledige laatste regel."""
        if self._partial:
            self.feed_line(self._partial.rstrip(u"\r"))
            self._partial = u""
        if self._pending:
            line, self._pending = self._pending, u""
            self._handle(line)

    def feed_line(self, line):
        self.data.lines += 1
        if self._pending:
            # Doorloopregel (residu-plaats of "au numéro d'ordre :").
            line = self._pending.rstrip() + u" " + line.strip()
            self._pending = u""
        elif (u"au noeud et degr" in line and u">" not in
              line.split(u"libert", 1)[-1]):
            self._pending = line
            return
        self._handle(line)

    # ---------- verwerking ----------
    def _current(self):
        return self.data.steps[-1] if self.data.steps else None

    def _add_message(self, kind, code, text, hint=u""):
        """Melding toevoegen; een herhaalde melding telt op (count) en
        krijgt de laatste poging/tijdstip en een nieuw volgnummer (seq)."""
        step = self._current()
        self._seq += 1
        key = (kind, code, text[:80])
        msg = self._msg_seen.get(key)
        if msg is None:
            msg = {"kind": kind, "code": code, "text": text, "hint": hint,
                   "count": 0}
            self._msg_seen[key] = msg
            self.data.messages.append(msg)
        msg["count"] += 1
        msg["seq"] = self._seq
        msg["step"] = step.index if step else None
        msg["instant"] = step.instant if step else None

    _seq = 0

    def _handle(self, line):
        data = self.data
        stripped = line.strip()
        if stripped.startswith(u"<I> <FIN>"):
            data.fin = True
            return

        # --- meldingskader ---
        if self._box is not None:
            if _BOX_BOTTOM in line:
                kind, code, body = self._box
                self._box = None
                text = u" ".join(body)
                if len(text) > 400:
                    text = text[:400] + u" ..."
                if kind is None:
                    return   # <I>-kader (bijv. lijst van alarmen)
                self._add_message(kind, code, text, self._hints.get(code, u""))
                if kind == "fout":
                    step = self._current()
                    if step is not None and step.status == "bezig":
                        step.status = "mislukt"
                    if self._lin is not None:
                        self._lin.status = "mislukt"
                return
            part = stripped.strip(_BOX_SIDE).strip()
            m = _BOX_CODE_RE.search(part)
            if m and not self._box[1]:
                self._box[0] = self._kind(m.group(1))
                self._box[1] = m.group(2)
            elif part:
                self._box[2].append(part)
            return
        if _BOX_TOP in line:
            self._box = [None, u"", []]
            return
        m = _BOX_CODE_RE.search(line)
        if m and _BOX_SIDE in line:
            # Kader zonder herkende bovenrand (afgekapt bestand).
            self._box = [self._kind(m.group(1)), m.group(2), []]
            return

        # --- lineaire oplosser (MECA_STATIQUE) ---
        if self._lin is not None:
            if self._lin_line(line, stripped):
                return
        elif _LIN_CMD_RE.match(line):
            self._lin_start()
            self._lin_line(line, stripped)
            return

        # --- modaal / knik ---
        if _MODES_CMD_RE.match(line):
            data.command = u"CALC_MODES"
            return
        m = _MODE_LINE_RE.match(stripped)
        if m:
            values = {}
            for part in m.group(2).split():
                key, _sep, val = part.partition(u"=")
                if _sep:
                    values[key] = _to_float(val)
            data.modes.append((int(m.group(1)), values))
            return

        # --- convergentietabel ---
        if stripped.startswith(u"|"):
            self._table_line(stripped)
            return
        if stripped.startswith(u"---"):
            if self._header:
                self._finish_header()
            if _DIAG_MARK in stripped:
                data.diagnostic = stripped.split(u":", 1)[1].strip()
            return

        m = _INSTANT_RE.search(line)
        if m:
            step = Step(len(data.steps), _to_float(m.group(1)),
                        int(m.group(2)),
                        self._clock() if self._clock else None)
            if self._contact_init is not None:
                step.contact_initial = self._contact_init
                self._contact_init = None
            data.steps.append(step)
            self._columns = None
            return

        step = self._current()
        if stripped.startswith(u"<Erreur>"):
            text = stripped[len(u"<Erreur>"):].strip()
            if step is not None:
                step.error = text
                step.status = "mislukt"
            # Een mislukte integratie van de materiaalwet lost Code-Aster op
            # door de stap op te delen (o.a. na elk bezwijkmoment); pas als
            # de run stopt is het een fout, en die komt apart als kader.
            kind = "alarm" if u"intégration" in text else "fout"
            self._add_message(kind, u"", text, event_hint(text))
            return
        if stripped.startswith(u"<Action>"):
            text = stripped[len(u"<Action>"):].strip()
            text = text.replace(u"<Échec>", u"").strip()
            if step is not None:
                step.actions.append(text)
            kind = "fout" if (u"chec" in text or u"arrête" in text) \
                else "actie"
            self._add_message(kind, u"", text, event_hint(text))
            return
        if stripped.startswith(_CONVERGED):
            if step is not None:
                step.status = "geconvergeerd"
            return
        if stripped.startswith(u"Découpe uniforme") and step is not None:
            step.actions.append(stripped)
            return

        m = _RESI_LOC_RE.search(line)
        if m and step is not None:
            step.resi_location = tr(u"knoop {0} {1}",
                                    u"node {0} {1}").format(m.group(3),
                                                            m.group(4))
            return
        m = _PIVOT_RE.search(stripped)
        if m:
            text = tr(u"Pivot nul bij knoop {0}, component {1}.",
                      u"Zero pivot at node {0}, component {1}.").format(
                m.group(1), m.group(2))
            self._add_message("fout", u"PIVOT", text, pivot_hint(m.group(2)))
            return
        m = _CPU_RE.search(line)
        if m and step is not None:
            step.cpu = _to_float(m.group(1))
            return
        if stripped.startswith(u"*"):
            m = _SUMMARY_RE.match(stripped)
            if m and step is not None:
                for pattern, key in _SUMMARY_KEYS:
                    if pattern in m.group(1):
                        step.summary[key] = _to_float(m.group(2))
                        break
            return
        m = _PROGRESS_RE.match(stripped)
        if m:
            data.progress.append((int(m.group(1)), _to_float(m.group(2))))
            return
        m = _TEND_RE.match(stripped)
        if m:
            data.t_end = _to_float(m.group(1))
            return
        if not data.explicit and _EXPLICIT_RE.search(line):
            data.explicit = True
        m = _CONTACT_INIT_RE.search(line)
        if m:
            self._contact_init = (int(m.group(1)), int(m.group(2)))
            return
        if data.command is None or not data.steps:
            m = _COMMAND_RE.search(line)
            if m and data.command in (None, "MECA_STATIQUE"):
                data.command = m.group(1)
            if data.command is not None and not data.steps:
                m = _RESI_CRIT_RE.search(line)
                if m:
                    data.resi_crit = _to_float(m.group(1))
                m = _ITER_MAX_RE.search(line)
                if m:
                    data.iter_max = int(m.group(1))
        m = _PARTIAL_RE.search(line)
        if m:
            data.partial_until = _to_float(m.group(1))
            self._add_message("fout", u"BRIDGE_ONVOLLEDIG",
                              tr(u"Niet geconvergeerd; resultaten tot "
                                 u"t = {0} s.",
                                 u"Not converged; results up to "
                                 u"t = {0} s.").format(m.group(1)))
            return
        m = _FAILURE_RE.search(line)
        if m:
            self._add_message("alarm", u"BRIDGE_BEZWIJKEN", tr(
                u"Bezwijken bij t = {0} s: {1} element(en) van elset {2} "
                u"(totaal {3}).",
                u"Failure at t = {0} s: {1} element(s) of elset {2} "
                u"(total {3}).").format(m.group(1), m.group(3), m.group(2),
                                        m.group(4)), tr(
                u"Hun stijfheid is verlaagd en hun spanning is nul; de "
                u"volgende stap begint met een elastische voorspelling en "
                u"wordt zo nodig opgedeeld.",
                u"Their stiffness has been reduced and their stress is "
                u"zero; the next step starts with an elastic prediction "
                u"and is subdivided if needed."))
            return
        if _DIAG_MARK in line:
            data.diagnostic = line.split(u":", 1)[1].strip()

    _contact_init = None

    @staticmethod
    def _kind(tag):
        return "alarm" if tag == "A" else "fout"

    # ---------- lineaire oplosser ----------
    def _now(self):
        return self._clock() if self._clock else None

    def _lin_start(self):
        lin = LinearSolve(self._now())
        self.data.linear.append(lin)
        if self.data.command is None:
            self.data.command = "MECA_STATIQUE"
        self._lin = lin
        self._lin_refine = None
        self._lin_group = u""
        self._lin_phase("matrix")

    def _lin_phase(self, key):
        """Fase 'key' begint (de lopende fase eindigt)."""
        lin = self._lin
        now = self._now()
        if lin.current is not None and lin.current != key:
            ph = lin.phases.get(lin.current)
            if ph is not None and ph.get("end") is None:
                ph["end"] = now
        ph = lin.phases.setdefault(key, {"start": now, "end": None,
                                         "seconds": None})
        if ph.get("start") is None:
            ph["start"] = now
        ph["end"] = None
        lin.current = key

    def _lin_phase_end(self, key, seconds=None):
        lin = self._lin
        ph = lin.phases.setdefault(key, {"start": None, "end": None,
                                         "seconds": None})
        ph["end"] = self._now()
        if seconds is not None:
            ph["seconds"] = seconds
        if lin.current == key:
            lin.current = None

    def _lin_finish(self, stripped):
        lin = self._lin
        m = _LIN_FIN_RE.search(stripped)
        if m:
            lin.cpu_total = _to_float(m.group(1))
            lin.elapsed_total = _to_float(m.group(2))
        if lin.current is not None:
            self._lin_phase_end(lin.current)
        # Tijden uit het Statistics-blok (Python-timers van Code-Aster,
        # wandklok) gaan voor de klok van het scherm; de MUMPS-fasen houden
        # de 'Elapsed time in ... driver' van MUMPS zelf.
        st = lin.stats
        solver = u"GcpcSolver" if lin.method() == "GCPC" else u"MumpsSolver"

        def put(key, value):
            if value is not None:
                ph = lin.phases.setdefault(key, {"start": None, "end": None,
                                                 "seconds": None})
                ph["seconds"] = value
        put("matrix", st.get(u"_computeMatrix"))
        put("rhs", st.get(u"_computeRhs"))
        post = [st[k] for k in (u"_computeStress",
                                u"StorageManager.storeState") if k in st]
        if post:
            put("post", sum(post))
        if lin.method() == "GCPC":
            put("factor", st.get(solver + u".factorize"))
            put("solve", st.get(solver + u".solve"))
        else:
            fac = lin.phases.get("factor", {})
            if fac.get("seconds") is None and solver + u".factorize" in st:
                ana = lin.phases.get("analysis", {}).get("seconds") or 0.0
                put("factor", max(0.0, st[solver + u".factorize"] - ana))
            if lin.phases.get("solve", {}).get("seconds") is None:
                put("solve", st.get(solver + u".solve"))
        if lin.status == "bezig":
            lin.status = "klaar"
        self._lin = None
        self._lin_refine = None

    def _lin_line(self, line, stripped):
        """Regel binnen MECA_STATIQUE; True als hij verwerkt is."""
        lin = self._lin
        # Code-Aster breekt lange zinnen af ("taux de remplissage" /
        # "de   0.028 %."): de vorige regel gaat mee voor _LIN_FILL_RE.
        joined = self._lin_prev + u" " + stripped
        self._lin_prev = stripped
        if stripped.startswith(_LIN_END):
            self._lin_finish(stripped)
            return True
        if lin.dofs is None:
            # Commando-echo: de waarden van SOLVEUR (incl. de defaults).
            for m in _LIN_TXT_KW_RE.finditer(line):
                lin.solver.setdefault(m.group(1), m.group(2))
            for m in _LIN_NUM_KW_RE.finditer(line):
                value = _to_float(m.group(2))
                if value is not None:
                    if m.group(1) in ("NPREC", "NMAX_ITER",
                                      "NIVE_REMPLISSAGE"):
                        value = int(value)
                    lin.solver.setdefault(m.group(1), value)
            m = _LIN_DOF_RE.search(line)
            if m:
                lin.dofs = int(m.group(1))
                return True
            return False
        # --- matrix (eerste keer; GCPC herhaalt dit na RCMK) ---
        if stripped.startswith(u"-"):
            m = _LIN_PHYS_RE.match(stripped)
            if m and lin.dofs_phys is None:
                lin.dofs_phys = int(m.group(1))
                return True
            m = _LIN_LAGR_RE.match(stripped)
            if m and lin.lagr_pairs is None:
                lin.lagr_pairs = int(m.group(1))
                return True
        for rx, attr, conv in ((_LIN_NODES_RE, "nodes", int),
                               (_LIN_REL_RE, "relations", int),
                               (_LIN_EQ_RE, "equations", int),
                               (_LIN_NNZ_RE, "nnz", int),
                               (_LIN_FILL_RE, "fill_pct", float)):
            m = rx.search(joined if attr == "fill_pct" else line)
            if m:
                if getattr(lin, attr) is None:
                    setattr(lin, attr, conv(m.group(1)))
                return True
        m = _LIN_ASSEMB_RE.search(line)
        if m:
            lin.detailed = True
            if m.group(1) == "M":
                self._lin_phase_end("matrix")
            elif lin.current != "rhs":
                self._lin_phase("rhs")
            return True
        if u"RENUMEROTATION DES NOEUDS" in line:
            self._lin_phase("factor")      # GCPC: RCMK + LDLT_INC
            return True
        m = _LIN_JOB_RE.search(line)
        if m:
            lin.detailed = True
            lin.mumps.setdefault("version", m.group(1))
            job = int(m.group(2))
            key = {1: "analysis", 4: "analysis", 6: "analysis",
                   2: "factor", 5: "factor", 3: "solve"}.get(job)
            if key is not None:
                self._lin_phase(key)
            return True
        m = _LIN_DRIVER_RE.search(line)
        if m:
            key = {"analysis": "analysis", "factorization": "factor",
                   "solve": "solve"}[m.group(1)]
            self._lin_phase_end(key, _to_float(m.group(2)))
            return True
        m = _LIN_ORDER_RE.search(line)
        if m:
            lin.mumps["ordering"] = m.group(1)
            return True
        if u"Start factorization" in line:
            m = _LIN_FLOPS_RE.search(line)
            if m:
                lin.mumps["flops"] = _to_float(m.group(1))
            return True
        m = _LIN_MEM_RE.search(line)
        if m:
            key = {"17": "mem_est", "19": "mem_alloc",
                   "22": "mem_used"}[m.group(1)]
            lin.mumps[key] = int(m.group(2))
            return True
        m = _LIN_FLOPS_EST_RE.search(line)
        if m:
            lin.mumps.setdefault("flops", _to_float(m.group(1)))
            return True
        m = _LIN_ENTRIES_EST_RE.search(line)
        if m:
            lin.mumps["entries_est"] = int(m.group(1))
            return True
        m = _LIN_INFOG_RE.search(line)
        if m:
            key = {"1": "infog1", "11": "front", "12": "neg_pivots",
                   "13": "delayed", "29": "entries"}.get(m.group(1))
            if key is not None:
                lin.mumps[key] = _to_float(m.group(2))
            return True
        if u"ERROR ANALYSIS BEFORE ITERATIVE REFINEMENT" in line:
            self._lin_refine = "before"
            return True
        if u"STATISTICS AFTER ITERATIVE REFINEMENT" in line:
            self._lin_refine = "after"
            return True
        if u"END   ITERATIVE REFINEMENT" in line or \
                u"END ITERATIVE REFINEMENT" in line:
            self._lin_refine = None
            return True
        m = _LIN_REFINE_STEPS_RE.search(line)
        if m:
            lin.refine_steps = int(m.group(1))
            return True
        if self._lin_refine is not None:
            target = (lin.refine_before if self._lin_refine == "before"
                      else lin.refine_after)
            m = _LIN_RESID_RE.match(stripped)
            if m:
                target["resid"] = _to_float(m.group(1))
                return True
            m = _LIN_RINFOG_RE.match(stripped)
            if m:
                target[int(m.group(1))] = _to_float(m.group(2))
                return True
        # --- GCPC ---
        m = _LIN_GCPC_INIT_RE.match(stripped)
        if m:
            lin.detailed = True
            lin.gcpc_init = _to_float(m.group(1))
            self._lin_phase("solve")
            return True
        m = _LIN_GCPC_TARGET_RE.search(stripped)
        if m:
            lin.gcpc_target = (_to_float(m.group(1)), _to_float(m.group(2)))
            return True
        m = _LIN_GCPC_ITER_RE.match(stripped)
        if m:
            lin.gcpc_iters.append((int(m.group(1)), _to_float(m.group(2)),
                                   _to_float(m.group(3))))
            return True
        m = _LIN_GCPC_FINAL_RE.search(stripped)
        if m:
            lin.gcpc_final = tuple(_to_float(m.group(i)) for i in (1, 2, 3))
            return True
        m = _LIN_GCPC_CONV_RE.match(stripped)
        if m:
            lin.gcpc_converged = int(m.group(1))
            self._lin_phase_end("solve")
            return True
        # --- na de oplossing ---
        if stripped.startswith(u"Champ stock"):
            self._lin_phase("post")
            return True
        m = _LIN_STAT_GROUP_RE.match(line)
        if m:
            self._lin_group = m.group(1)
            return True
        m = _LIN_STAT_RE.match(line)
        if m:
            name = m.group(2)
            if m.group(1):
                name = u"{0}.{1}".format(self._lin_group, name)
            else:
                self._lin_group = u""
            lin.stats[name] = _to_float(m.group(4))
            return True
        return False

    # ---------- tabellen ----------
    def _cells(self, line):
        return line.strip().strip(u"|").split(u"|")

    def _table_line(self, line):
        cells = self._cells(line)
        first = cells[0].strip() if cells else u""
        if first.startswith(u"BILAN D'ENERGIE"):
            self._energy_cols = [c.strip() for c in cells[1:]]
            return
        if first in (u"PAS COURANT", u"TOTAL"):
            step = self._current()
            if first == u"TOTAL" and step is not None and self._energy_cols:
                for name, cell in zip(self._energy_cols, cells[1:]):
                    step.energy[name] = _to_float(cell.strip())
            return
        if not _NUM_CELL_RE.match(cells[0]):
            if self._header and len(cells) != len(self._header[0]):
                self._header = []
            self._header.append(cells)
            return
        if self._header:
            self._finish_header()
        if not self._columns or len(cells) != len(self._columns):
            return
        step = self._current()
        if step is None:
            return
        it = {"extra": {}, "flags": {}}
        for key, name, cell in zip(self._columns, self._names, cells):
            m = _NUM_CELL_RE.match(cell)
            if key == "option":
                it["option"] = cell.strip()
                continue
            value = None
            if m:
                value = _to_float(m.group(1))
                if m.group(2):
                    it["flags"][key or name] = m.group(2)
            if key is None:
                if value is not None:
                    it["extra"][name] = value
            else:
                it[key] = value
        if it.get("newton") is None:
            return
        step.iterations.append(it)

    def _finish_header(self):
        rows = self._header
        self._header = []
        if not rows:
            return
        n = len(rows[0])
        names = []
        for i in range(n):
            parts = [r[i].strip() for r in rows if i < len(r) and r[i].strip()]
            names.append(u" ".join(parts))
        self._names = names
        self._columns = [_column_key(nm) for nm in names]


# ==========================================
# HULPFUNCTIES
# ==========================================
def read_text(path):
    """Tekst van een .mess/fort.6 (UTF-8, onleesbare tekens vervangen)."""
    f = io.open(path, "rb")
    try:
        return f.read().decode("utf-8", "replace")
    finally:
        f.close()


def parse_file(path, hints=None):
    """Lees een volledig .mess-bestand. Returns ConvergenceData."""
    parser = ConvergenceParser(hints=hints)
    parser.feed(read_text(path))
    parser.close()
    return parser.data


def modes_text(data, limit=6):
    """Korte tekst met de eerste modes ('1: 7.074 Hz, 2: ...' of
    kniklastfactoren), of '' zonder modes."""
    out = []
    for nr, values in data.modes[:limit]:
        if values.get("FREQ") is not None:
            # %-opmaak: format() volgt in IronPython de landinstelling
            out.append(u"{0}: {1} Hz".format(nr, "%.4g" % values["FREQ"]))
        elif values.get("CHAR_CRIT") is not None:
            out.append(u"{0}: {1}".format(nr, "%.4g" % values["CHAR_CRIT"]))
    if len(data.modes) > limit:
        out.append(u"...")
    return u", ".join(out)


def summarize(data):
    """Korte tekstsamenvatting (log, tests)."""
    steps = data.steps
    conv = [s for s in steps if s.status == "geconvergeerd"]
    fail = [s for s in steps if s.status == "mislukt"]
    its = sum(len(s.iterations) for s in steps)
    last = data.last_converged()
    lines = [
        u"commando={0} RESI_GLOB_RELA={1} ITER_GLOB_MAXI={2}".format(
            data.command, data.resi_crit, data.iter_max),
        u"pogingen={0} geconvergeerd={1} mislukt={2} iteraties={3} "
        u"max_niveau={4}".format(
            len(steps), len(conv), len(fail), its,
            max([s.level for s in steps] or [0])),
        u"laatste geconvergeerd t={0} voortgang={1} diagnose={2} "
        u"fin={3}".format(
            last.instant if last else None,
            data.overall_progress(),
            data.diagnostic or u"(bezig)", data.fin),
    ]
    for lin in data.linear:
        lines.append(
            u"lineair {0} {1}: vergelijkingen={2} nnz={3} status={4} "
            u"fout={5} cond={6} gcpc_iter={7} fasen={8}".format(
                lin.command, lin.method(), lin.equations, lin.nnz,
                lin.status, lin.forward_error(), lin.condition(),
                lin.gcpc_converged,
                u",".join(u"%s:%.3g" % (k, lin.phase_seconds(k))
                          for k in lin.phase_keys()
                          if lin.phase_seconds(k) is not None)))
    for nr, values in data.modes:
        lines.append(u"mode {0}: {1}".format(nr, u" ".join(
            u"{0}={1}".format(k, values[k]) for k in sorted(values))))
    codes = []
    for msg in data.messages:
        c = msg["code"] or msg["text"][:40]
        if c not in codes:
            codes.append(c)
    lines.append(u"meldingen: " + u"; ".join(codes))
    return u"\n".join(lines)


# ==========================================
# PYTHON-FOUTEN IN DE CONSOLE-UITVOER
#
# Een Python-fout in study.comm (bijv. een tikfout in #ASTER-regels) staat
# niet in fort.6/study.mess: run_aster stuurt alleen stdout van de job via
# "tee" naar fort.6, de traceback gaat naar stderr (het Code-Aster-venster).
# Gezien met Code-Aster 17.4: een SyntaxError geeft dan "DIAGNOSTIC JOB : OK"
# zonder <FIN>, een NameError "<F>_ABNORMAL_ABORT". De bridge legt de
# console-uitvoer vast in study_console.txt (bridge_run_aster.py).
# ==========================================
_PY_ERR_RE = re.compile(
    r"^([A-Za-z_][\w.]*(?:Error|Exception|Interrupt|Exit)):\s*(.*)$")
_PY_FILE_RE = re.compile(r'^\s*File "([^"]+)", line (\d+)')


def python_errors(text, max_errors=5):
    """
    Python-fouten (traceback of SyntaxError) uit console-uitvoer.

    Alleen een foutregel ("NameError: ...") aan het begin van een regel,
    na een traceback of na een 'File "...", line n'-regel telt; zo geven
    Code-Aster-meldingen zoals <ConvergenceError> geen valse treffers.

    Returns
    -------
    list of unicode - bijv. "SyntaxError: keyword argument repeated:
                      INCREMENT (regel: 'INCREMENT=_F(...)')", uniek, in
                      volgorde van voorkomen.
    """
    lines = text.splitlines()
    found = []
    for i, line in enumerate(lines):
        m = _PY_ERR_RE.match(line.rstrip())
        if not m or m.group(1).startswith(u"libaster."):
            # libaster.AsterError e.d.: Code-Aster-fout met een eigen
            # meldingskader in study.mess (gezien bij FACTOR_11, DVP_2).
            continue
        frame = None
        in_traceback = False
        for j in range(i - 1, max(-1, i - 60), -1):
            prev = lines[j]
            if prev.startswith(u"Traceback (most recent call last)"):
                in_traceback = True
                if j > 0 and lines[j - 1].strip() == u"--- Logging error ---":
                    # logging.Handler.handleError: een logregel kon niet
                    # naar de console (geen console, bijv. een job zonder
                    # venster: "De ingang is ongeldig"); de job rekent
                    # gewoon door (DIAGNOSTIC JOB : OK, gezien 17.4)
                    in_traceback = None
                break
            fm = _PY_FILE_RE.match(prev)
            if fm and (frame is None or u".comm" in fm.group(1)):
                frame = (fm.group(1), j)    # frame in study.comm wint
            if not fm and not prev.startswith(u" ") and prev.strip():
                break
        if in_traceback is None or (frame is None and not in_traceback):
            continue
        msg = u"{0}: {1}".format(m.group(1), m.group(2).strip())
        if frame is not None and frame[1] + 1 < i:
            code = lines[frame[1] + 1].strip()
            if (code and not code.startswith(u"^") and
                    not code.startswith(u"File ")):
                if len(code) > 120:
                    code = code[:120] + u"..."
                msg += tr(u" (regel: '{0}')", u" (line: '{0}')").format(code)
        if msg not in found:
            found.append(msg)
        if len(found) >= max_errors:
            break
    return found


# ==========================================
# LIVE fort.6 VAN EEN LOPENDE RUN
#
# run_aster (Code-Aster 17.4, lib/aster/run_aster/run.py) draait de job in
# een tijdelijke map (tempfile.mkdtemp(prefix="run_aster_"), in de map
# "tmpdir" uit share/aster/config.json, anders %TEMP%) en schrijft de
# uitvoer via "tee -a fort.6" live in die map. Pas aan het eind kopieert
# hij fort.6 naar het .mess-bestand uit de .export en verwijdert hij de
# map. In de map staat een kopie van de .export ("<pid>.export") met de
# regel "F mess <pad> R 6"; daaraan herkennen we de run.
# ==========================================
def _norm(path):
    return os.path.normcase(os.path.normpath(path.replace("/", os.sep)))


def run_dir_mess(run_dir):
    """Het .mess-pad uit de .export in een run_aster-map, of None."""
    try:
        names = os.listdir(run_dir)
    except OSError:
        return None
    for name in names:
        if not name.endswith(".export"):
            continue
        try:
            text = read_text(os.path.join(run_dir, name))
        except (IOError, OSError):
            continue
        for line in text.splitlines():
            parts = line.split()
            if len(parts) >= 4 and parts[0] == "F" and parts[1] == "mess":
                return u" ".join(parts[2:-2])
    return None


def find_live_fort6(tmp_dirs, mess_path, since=0.0):
    """
    Zoek de fort.6 van de lopende run die naar mess_path schrijft.

    tmp_dirs : lijst mappen waarin run_aster zijn run_aster_*-map maakt.
    since    : alleen mappen die na dit tijdstip (time.time()) zijn gemaakt
               (2 s marge); 0 = alle.

    Returns
    -------
    (pad naar fort.6, aanmaaktijd van de map) van de nieuwste passende
    map, of (None, None).
    """
    target = _norm(mess_path)
    best = (None, None)
    for tmp in tmp_dirs:
        try:
            names = os.listdir(tmp)
        except OSError:
            continue
        for name in names:
            if not name.startswith("run_aster_"):
                continue
            run_dir = os.path.join(tmp, name)
            fort6 = os.path.join(run_dir, "fort.6")
            try:
                created = os.path.getctime(run_dir)
            except OSError:
                continue
            if created < since - 2.0 or not os.path.isfile(fort6):
                continue
            mess = run_dir_mess(run_dir)
            if mess is None or _norm(mess) != target:
                continue
            if best[1] is None or created > best[1]:
                best = (fort6, created)
    return best
