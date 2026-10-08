# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# bridge_radioss.py
# Verantwoordelijkheid: OpenRadioss-invoer uit de .liml lezen die Mecway
# voor de OpenRadioss-solver bewaart, voor zover de bridge die naar
# Code-Aster kan vertalen:
#   - <radiossmodeldefinition>: /FUNCT/<id> (tabellen) en /BEGIN (alleen de
#     eenheden, voor een melding); andere kaarten worden alleen opgesomd
#     (melding "niet vertaald");
#   - <mat><openradioss type="Custom" law="PLAS_TAB">: /MAT/LAW36
#     (elastisch-plastisch met een tabel vloeispanning - plastische rek)
#     -> E, NU, RHO en een hardingscurve; een /FAIL/JOHNSON met een vaste
#     breukrek in hetzelfde veld -> eps_p_max (bridge_failure.py);
#   - law="VOID": materiaal zonder stijfheid (in Radioss o.a. voor starre
#     lichamen met /RBODY); niet te vertalen.
#
# Kaartformaat (OpenRadioss/Altair Radioss Reference Guide, /MAT/LAW36,
# /FUNCT, /FAIL/JOHNSON, /BEGIN): vaste kolommen van 10 of 20 tekens;
# regels die met '#' beginnen zijn commentaar, een regel die met '/' begint
# start een nieuwe kaart. De waarden staan in het eenhedenstelsel van het
# Radioss-model. Mecway rekent in SI (de templates: RHO 7800, E 206E9,
# vloeispanning 228.03E6) en zet bij de export zelf een /BEGIN ('kg' staat
# bij de /BEGIN-tekst in mecway.exe); de bridge rekent niets om en meldt
# waarden die op een ander stelsel wijzen (_check_units).
#
# Compatibel met IronPython 2.7.9.

from bridge_utils import log, tr

# Radioss-standaard voor "niet actief" bij breukrekken (Eps_p_max, Eps_t,
# Eps_m, Eps_f: 1e20 tot 3e20 volgens de Reference Guide /MAT/LAW36); een
# waarde vanaf deze grens telt als niet ingevuld.
_RADIOSS_INFINITE = 1.0e19

# SI-invoereenheden in /BEGIN (massa, lengte, tijd)
_SI_UNITS = ("kg", "m", "s")


def _data_lines(text):
    """Regels zonder commentaar (#), tot de eerste kaartregel (/...)."""
    out = []
    for line in (text or "").splitlines():
        if line.startswith("#"):
            continue
        if line.startswith("/"):
            break
        out.append(line.rstrip("\r"))
    return out


def _field(line, start, width):
    """Veld uit een vaste kolom als float, of None als het leeg is."""
    txt = line[start:start + width].strip() if line else ""
    if not txt:
        return None
    try:
        return float(txt.replace("d", "e").replace("D", "E"))
    except ValueError:
        return None


def _active(value):
    """True voor een ingevulde breukrek (> 0 en onder de Radioss-standaard
    'oneindig'); 0 of leeg betekent in Radioss de standaardwaarde."""
    return value is not None and 0.0 < value < _RADIOSS_INFINITE


def parse_cards(text):
    """
    [(kaartnaam, [dataregels])] uit een Radioss-tekst, bijv.
    ("/FUNCT/99", ["TRUE_Stress_strain", "  0.0  228.03E6", ...]).
    Commentaarregels vallen weg.
    """
    cards = []
    for line in (text or "").splitlines():
        line = line.rstrip("\r")
        if line.startswith("#"):
            continue
        if line.startswith("/"):
            cards.append((line.split()[0], []))
        elif cards:
            cards[-1][1].append(line)
    return cards


def _point(line):
    """
    (x, y) van een /FUNCT-regel, of None. Eerst de vaste kolommen van de
    Reference Guide (X in kolom 1-20, Y in 21-40): twee volle velden staan
    dan zonder spatie tegen elkaar ('-1.0000000000E+00-2.0E+00'), wat
    split() niet kan scheiden. Lukt dat niet (vrij formaat, bijv.
    '  0.0  228.03E6', dat over de grens van kolom 20 loopt), dan split().
    Vrij formaat dat in de vaste kolommen wel leesbaar is, geeft dezelfde
    waarden: elk veld moet dan precies een getal bevatten.
    """
    x, y = _field(line, 0, 20), _field(line, 20, 20)
    if x is not None and y is not None:
        return (x, y)
    vals = line.split()
    if len(vals) < 2:
        return None
    try:
        return (float(vals[0].replace("d", "e").replace("D", "E")),
                float(vals[1].replace("d", "e").replace("D", "E")))
    except ValueError:
        return None


def parse_functions(text):
    """{id: [(x, y), ...]} van alle /FUNCT/<id>-kaarten."""
    funcs = {}
    for name, lines in parse_cards(text):
        parts = name.split("/")
        if len(parts) < 3 or parts[1] != "FUNCT":
            continue
        try:
            fid = int(parts[2])
        except ValueError:
            continue
        points = []
        for line in lines[1:]:          # eerste regel = titel
            point = _point(line)
            if point is None:
                if line.strip():
                    log(tr("WARNING: /FUNCT/{0}: regel '{1}' is geen punt "
                           "(X, Y) en wordt overgeslagen.",
                           "WARNING: /FUNCT/{0}: line '{1}' is not a point "
                           "(X, Y) and is skipped.").format(
                               fid, line.strip()))
                continue
            points.append(point)
        funcs[fid] = points
    return funcs


def other_cards(text):
    """Namen van kaarten die de bridge niet vertaalt (alles behalve
    /FUNCT)."""
    return [name for name, _ in parse_cards(text)
            if not name.startswith("/FUNCT/")]


def begin_units(text):
    """
    Invoereenheden (massa, lengte, tijd) uit een /BEGIN-kaart, of None.
    /BEGIN (Reference Guide): titel, Invers/Irun, invoereenheden (3 velden
    van 20 tekens), werkeenheden.
    """
    for name, lines in parse_cards(text):
        if name != "/BEGIN" or len(lines) < 3:
            continue
        units = tuple(lines[2][20 * k:20 * k + 20].strip()
                      for k in range(3))
        if all(units):
            return units
    return None


def _check_units(mat_name, e_mod, rho, units, notes):
    """
    Meld een waarschijnlijk niet-SI-eenhedenstelsel (de waarden worden
    ongewijzigd gebruikt). Grenzen: E van vaste stoffen ligt in SI tussen
    ~1e6 Pa (rubber) en ~1e12 Pa (diamant), dichtheden tussen ~1 en
    ~25000 kg/m3; N-mm-t geeft E ~2e5 en RHO ~8e-9.
    """
    reasons = []
    if units is not None and tuple(u.lower() for u in units) != _SI_UNITS:
        reasons.append("/BEGIN {0} {1} {2}".format(*units))
    if e_mod < 1.0e6 or e_mod > 1.0e13:
        reasons.append("E=%.4g" % e_mod)
    if rho and (rho < 1.0 or rho > 1.0e5):
        reasons.append("RHO=%.4g" % rho)
    if not reasons:
        return
    log(tr("WARNING: PLAS_TAB van materiaal '{0}': de waarden lijken niet "
           "in SI (kg, m, s) te staan ({1}). De bridge rekent niet om; "
           "Mecway rekent in SI. Zet E, RHO en de /FUNCT-spanningen in Pa "
           "en kg/m3.",
           "WARNING: PLAS_TAB of material '{0}': the values do not seem to "
           "be in SI (kg, m, s) ({1}). The bridge does not convert; Mecway "
           "works in SI. Give E, RHO and the /FUNCT stresses in Pa and "
           "kg/m3.").format(mat_name, ", ".join(reasons)))
    notes.append(tr("PLAS_TAB van '{0}': waarschijnlijk geen SI-eenheden "
                    "({1}), niet omgerekend",
                    "PLAS_TAB of '{0}': probably no SI units ({1}), not "
                    "converted").format(mat_name, ", ".join(reasons)))


def _ignored_fields(mat_name, lines, n_funct, notes):
    """
    Meld velden van /MAT/LAW36 die de bridge niet vertaalt (Reference Guide
    /MAT/LAW36; standaarden tussen haakjes):
      Eps_t, Eps_m (1e20, 2e20): spanning neemt af vanaf hoofdrek Eps_t en
          is nul bij Eps_m;
      Eps_f (3e20): element weg bij trekrek Eps_f;
      C_hard (0 = isotroop, 1 = kinematisch Prager-Ziegler, ertussen
          gemengd);
      fct_IDp (0): vloeigrens afhankelijk van de druk;
      Fct_IDE (0), EInf, CE: E neemt af met de plastische rek;
      F_smooth, F_cut, VP: alleen voor de reksnelheid.
    Wat de berekening verandert gaat ook in notes (skipped_loads).
    """
    eps_t = _field(lines[1], 60, 20)
    eps_m = _field(lines[1], 80, 20)
    f_smooth = _field(lines[2], 10, 10)
    c_hard = _field(lines[2], 20, 20)
    f_cut = _field(lines[2], 40, 20)
    eps_f = _field(lines[2], 60, 20)
    vp = _field(lines[2], 80, 20)
    fct_idp = _field(lines[3], 0, 10)
    fct_ide = _field(lines[3], 30, 10)
    e_inf = _field(lines[3], 40, 20)
    ce = _field(lines[3], 60, 20)

    changed = []          # (veld, Nederlandse uitleg, Engelse uitleg)
    if _active(eps_t) or _active(eps_m):
        changed.append(("Eps_t/Eps_m",
                        "afname van de spanning bij een grote hoofdrek",
                        "stress reduction at a large principal strain"))
    if _active(eps_f):
        changed.append(("Eps_f", "element weg bij trekrek Eps_f",
                        "element deletion at tensile strain Eps_f"))
    if c_hard:
        changed.append(("C_hard=%g" % c_hard,
                        "kinematische/gemengde versteviging; isotroop "
                        "gerekend (gelijk bij monotone belasting, anders bij "
                        "ontlasten of wisselbelasting)",
                        "kinematic/mixed hardening; computed as isotropic "
                        "(equal under monotonic loading, different on "
                        "unloading or load reversal)"))
    if fct_idp:
        changed.append(("fct_IDp", "drukafhankelijke vloeigrens",
                        "pressure-dependent yield stress"))
    if fct_ide or (ce and e_inf):
        changed.append(("Fct_IDE/EInf/CE",
                        "E afhankelijk van de plastische rek",
                        "E depending on the plastic strain"))
    for field, nl, en in changed:
        log(tr("WARNING: PLAS_TAB van materiaal '{0}': {1} ({2}) wordt niet "
               "vertaald.",
               "WARNING: PLAS_TAB of material '{0}': {1} ({2}) is not "
               "translated.").format(mat_name, field, tr(nl, en)))
        notes.append(tr("PLAS_TAB van '{0}': {1} ({2}) niet vertaald",
                        "PLAS_TAB of '{0}': {1} ({2}) not translated").format(
                            mat_name, field, tr(nl, en)))
    rate = [n for n, v in (("F_smooth", f_smooth), ("F_cut", f_cut),
                           ("VP", vp)) if v]
    if rate:
        log("INFO: PLAS_TAB '{0}': {1} (reksnelheid) niet gebruikt; "
            "Code-Aster rekent met een curve ({2} functie(s)).".format(
                mat_name, ", ".join(rate), n_funct))


def _johnson_strain(mat_name, lines):
    """
    Vaste breukrek uit /FAIL/JOHNSON, of None. Reference Guide
    /FAIL/JOHNSON: regel 1 D1..D5 (20 tekens elk);
    eps_f = [D1 + D2 exp(D3 s*)] [1 + D4 ln(eps_dot*)] [1 + D5 T*] en
    schade D = som(d eps_p / eps_f), bezwijken bij D >= 1. Met D4 = D5 = 0
    en D2 = 0 of D3 = 0 is eps_f constant (D1 + D2) en is D >= 1 hetzelfde
    als eps_p >= eps_f: precies het criterium van bridge_failure.py
    (equivalente plastische rek). Anders hangt eps_f af van de
    triaxialiteit, reksnelheid of temperatuur: niet vertaald.
    """
    if not lines:
        return None
    d = [_field(lines[0], 20 * k, 20) or 0.0 for k in range(5)]
    if d[3] != 0.0 or d[4] != 0.0 or (d[1] != 0.0 and d[2] != 0.0):
        return None
    eps_f = d[0] + (d[1] if d[2] == 0.0 else 0.0)
    if eps_f <= 0.0:
        return None
    log("INFO: PLAS_TAB '{0}': /FAIL/JOHNSON met vaste breukrek {1} "
        "(D1={2}, D2={3}, D3=D4=D5=0) -> bezwijkgrens plastische rek."
        .format(mat_name, eps_f, d[0], d[1]))
    return eps_f


def _failure_cards(mat_name, fieldcontents, notes):
    """
    Kaarten na de LAW36-gegevens in hetzelfde <fieldcontents> (in het
    Mecway-sjabloon staat daar een uitgecommentarieerde /FAIL/BIQUAD).
    Returns de breukrek van een vertaalbare /FAIL/JOHNSON of None; de rest
    wordt gemeld (ERROR + notes).
    """
    eps = None
    for name, lines in parse_cards(fieldcontents):
        value = None
        if name.startswith("/FAIL/JOHNSON"):
            value = _johnson_strain(mat_name, lines)
        if value is not None:
            eps = value if eps is None else min(eps, value)
            continue
        if name.startswith("/FAIL/JOHNSON"):
            why = tr("breukrek hangt af van triaxialiteit, reksnelheid of "
                     "temperatuur (D2*D3, D4 of D5 niet 0)",
                     "failure strain depends on triaxiality, strain rate "
                     "or temperature (D2*D3, D4 or D5 not 0)")
        elif name.startswith("/FAIL/"):
            why = tr("alleen /FAIL/JOHNSON met een vaste breukrek wordt "
                     "vertaald", "only /FAIL/JOHNSON with a constant "
                     "failure strain is translated")
        else:
            why = tr("kaart in het materiaalveld", "card in the material "
                     "field")
        log(tr("ERROR: PLAS_TAB van materiaal '{0}': {1} wordt niet "
               "vertaald ({2}).",
               "ERROR: PLAS_TAB of material '{0}': {1} is not translated "
               "({2}).").format(mat_name, name, why))
        notes.append(tr("PLAS_TAB van '{0}': {1} niet vertaald ({2})",
                        "PLAS_TAB of '{0}': {1} not translated ({2})"
                        ).format(mat_name, name, why))
    return eps


def parse_plas_tab(fieldcontents, functions, mat_name, units=None):
    """
    /MAT/LAW36 (PLAS_TAB) uit <fieldcontents>.

    Kolommen (OpenRadioss Reference Guide, /MAT/LAW36):
      1: RHO_I(20)
      2: E(20) Nu(20) Eps_p_max(20) Eps_t(20) Eps_m(20)
      3: N_funct(10) F_smooth(10) C_hard(20) F_cut(20) Eps_f(20) VP(20)
      4: fct_IDp(10) Fscale(20) Fct_IDE(10) EInf(20) CE(20)
      5: func_ID1..5 (10 elk)
      6: Fscale_1..5 (20 elk)
      7: Eps_dot_1..5 (20 elk)
    De functies geven de vloeispanning als functie van de plastische rek,
    per reksnelheid Eps_dot_i. Code-Aster (VMIS_ISOT_TRAC) is
    reksnelheidsonafhankelijk: de functie met de laagste reksnelheid wordt
    gebruikt. Niet vertaalde velden: _ignored_fields; /FAIL-kaarten in
    hetzelfde veld: _failure_cards. units: invoereenheden uit /BEGIN
    (begin_units) of None.

    Returns
    -------
    (props, notes)
        props : {"E", "NU", "rho", "plastic_curve", "curve_definition"} of
                None bij een fout
        notes : [str] - wat niet vertaald is (voor skipped_loads)
    """
    notes = []
    lines = _data_lines(fieldcontents)
    if len(lines) < 7:
        log(tr("ERROR: PLAS_TAB van materiaal '{0}': {1} dataregels, "
               "verwacht minstens 7.",
               "ERROR: PLAS_TAB of material '{0}': {1} data lines, "
               "expected at least 7.").format(mat_name, len(lines)))
        return None, [tr("PLAS_TAB van '{0}' onleesbaar",
                         "PLAS_TAB of '{0}' unreadable").format(mat_name)]
    rho = _field(lines[0], 0, 20)
    e_mod = _field(lines[1], 0, 20)
    nu = _field(lines[1], 20, 20)
    eps_p_max = _field(lines[1], 40, 20)
    n_funct = int(_field(lines[2], 0, 10) or 0)
    if e_mod is None or nu is None or n_funct < 1:
        log(tr("ERROR: PLAS_TAB van materiaal '{0}': E, Nu of N_funct "
               "ontbreekt.",
               "ERROR: PLAS_TAB of material '{0}': E, Nu or N_funct "
               "is missing.").format(mat_name))
        return None, [tr("PLAS_TAB van '{0}' onvolledig",
                         "PLAS_TAB of '{0}' incomplete").format(mat_name)]
    _check_units(mat_name, e_mod, rho, units, notes)
    _ignored_fields(mat_name, lines, n_funct, notes)
    if n_funct > 5:
        notes.append(tr("PLAS_TAB van '{0}': alleen de eerste 5 van {1} "
                        "functies gelezen",
                        "PLAS_TAB of '{0}': only the first 5 of {1} "
                        "functions read").format(mat_name, n_funct))
        n_funct = 5
    ids = [_field(lines[4], 10 * k, 10) for k in range(n_funct)]
    scales = [_field(lines[5], 20 * k, 20) for k in range(n_funct)]
    rates = [_field(lines[6], 20 * k, 20) for k in range(n_funct)]
    choice = sorted(range(n_funct), key=lambda k: rates[k] or 0.0)[0]
    fid = int(ids[choice] or 0)
    scale = scales[choice] if scales[choice] else 1.0
    if fid not in functions:
        log(tr("ERROR: PLAS_TAB van materiaal '{0}': /FUNCT/{1} niet "
               "gevonden in <radiossmodeldefinition>.",
               "ERROR: PLAS_TAB of material '{0}': /FUNCT/{1} not found "
               "in <radiossmodeldefinition>.").format(mat_name, fid))
        return None, [tr("PLAS_TAB van '{0}': /FUNCT/{1} ontbreekt",
                         "PLAS_TAB of '{0}': /FUNCT/{1} is missing").format(
            mat_name, fid)]
    if n_funct > 1:
        notes.append(tr("reksnelheidsafhankelijkheid van '{0}' (alleen de "
                        "curve bij Eps_dot={1} gebruikt)",
                        "strain-rate dependence of '{0}' (only the curve "
                        "at Eps_dot={1} used)").format(
                         mat_name, rates[choice] or 0.0))
    eps_fail = _failure_cards(mat_name, fieldcontents, notes)
    if eps_fail is not None and (not eps_p_max or eps_fail < eps_p_max):
        eps_p_max = eps_fail
    curve = [(x, y * scale) for x, y in functions[fid]]
    props = {
        "E": e_mod,
        "NU": nu,
        "Texp": 0.0,
        "rho": rho or 0.0,
        "plastic_type": "Isotropic",
        "plastic_curve": sorted(curve),
        # Radioss LAW36: x = plastische rek (Reference Guide)
        "curve_definition": "plastic",
        # Bezwijkgrens (plastische rek); bridge_failure.py maakt er
        # bezwijken van elementen van. 0 = geen bezwijken.
        "eps_p_max": eps_p_max or 0.0,
    }
    log("INFO: PLAS_TAB '{0}': E={1} Pa, NU={2}, rho={3}, /FUNCT/{4} met "
        "{5} punten (schaal {6}).".format(mat_name, e_mod, nu, rho, fid,
                                          len(curve), scale))
    return props, notes
