# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# aster_fatigue.py
#
# Vermoeiing (STAPPENPLAN S8) BINNEN de Code-Aster-job (CPython 3.11 met
# numpy, meegeleverd met Code-Aster 17.4), NIET in IronPython. De invoer
# leest bridge_fatigue.py (IronPython) uit de #ASTER-regels; die schrijft
# in de .comm na CALC_CHAMP:
#     import aster_fatigue
#     _bridge_vermoeiing = aster_fatigue.compute(res, liml, csv, curves,
#                                                options, blocks, signal)
# en geeft _bridge_vermoeiing als node_fields aan
# aster_results_export.export_results (velden damage en life in Mecway).
# Gebruikersdoc: docs/Vermoeiing.md (EN: docs/en/Fatigue.md).
#
# Werkwijze per knoop van een volume-element met een Wohler-curve:
#   1. equivalente spanning per tijdstap uit SIGM_NOEU: VMIS_SG (von Mises
#      met het teken van het spoor, zoals CALC_FATIGUE CUMUL_DOMMAGE,
#      catalogus calc_fatigue.py EQUI_GD='VMIS_SG') of PRIN (de
#      hoofdspanning met de grootste absolute waarde, met teken);
#   2. cycli: rainflow (RAINFLOW = het signaal als periodiek, beginnend bij
#      de grootste absolute waarde: alleen gesloten cycli, zoals
#      POST_FATIGUE COMPTAGE='RAINFLOW' en CALC_FATIGUE; gemeten: het
#      ASTM-signaal geeft in POST_FATIGUE dezelfde 4 cycli, szlz105a
#      dezelfde schade) of ASTM (ASTM E1049-85 par. 5.4.4, een doorgang,
#      restcycli als halve cyclus); bij een statische run een signaal of
#      blokken als factor op de statische spanning;
#   3. per cyclus N uit de Wohler-curve (amplitude S_a = (max - min)/2,
#      zoals DEFI_MATERIAU FATIGUE WOHLER: szlz100a cyclus 600/-500 ->
#      N(550)), eventueel met de gemiddelde-spanningscorrectie van
#      POST_FATIGUE (GOODMAN S_a/(1 - S_m/S_u), GERBER S_a/(1 - (S_m/
#      S_u)^2); gemeten: Code-Aster neemt het teken van S_m mee, ook bij
#      druk);
#   4. Miner: D = herhalingen x sum(n_i / N_i); life = 1 / sum(n_i / N_i)
#      herhalingen van het signaal of blokkenspectrum (begrensd op
#      LIFE_MAX).
# De knoop met de grootste schade krijgt een cyclustabel
# (<model>_<run_id>_vermoeiing.csv); regels BRIDGE_VERMOEIING_* in de
# .mess voor bridge_fatigue.report.
#
# Gebruik geen f-strings: zelfde stijl als de rest van de bridge.

import io
import math
import xml.etree.ElementTree as ET

import numpy as np

from aster_results_export import tr   # taal van de meldingen

# Levensduur zonder schade (D = 0) of daarboven: deze waarde (herhalingen);
# zo blijft de contour in Mecway leesbaar (keuze van de bridge,
# docs/Vermoeiing.md).
LIFE_MAX = 1.0e12

SOLID_SHAPES = frozenset(["hex8", "hex20", "tet4", "tet10", "pyr5", "pyr13",
                          "wedge6", "wedge15"])

_STRESS = ("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ")


# ==========================================
# CYCLI TELLEN
# ==========================================
def reversals(x):
    """Keerpunten van de reeks x (gelijke en tussenliggende punten op een
    monotoon stuk vallen weg). Returns lijst floats."""
    out = []
    for v in x:
        v = float(v)
        if out and v == out[-1]:
            continue
        if len(out) >= 2 and (out[-1] - out[-2]) * (v - out[-1]) > 0.0:
            out[-1] = v
        else:
            out.append(v)
    return out


def rainflow_astm(x):
    """
    Rainflow volgens ASTM E1049-85 (2017) par. 5.4.4, een doorgang: een
    bereik met het beginpunt telt als halve cyclus, de rest aan het eind
    ook. Returns [(laag, hoog, aantal)] met aantal 0,5 of 1.
    """
    st = []
    out = []
    for p in reversals(x):
        st.append(p)
        while len(st) >= 3:
            xr = abs(st[-1] - st[-2])
            yr = abs(st[-2] - st[-3])
            if xr < yr:
                break
            a, b = st[-3], st[-2]
            if len(st) == 3:
                out.append((min(a, b), max(a, b), 0.5))
                st.pop(0)
            else:
                out.append((min(a, b), max(a, b), 1.0))
                last = st.pop()
                st.pop()
                st.pop()
                st.append(last)
    for a, b in zip(st[:-1], st[1:]):
        out.append((min(a, b), max(a, b), 0.5))
    return out


def rainflow_periodic(x):
    """
    Rainflow voor een periodiek herhaald signaal (ASTM E1049 par. 5.4.4
    voor herhalende reeksen; POST_FATIGUE COMPTAGE='RAINFLOW'): begin bij
    de grootste absolute waarde, sluit de reeks met datzelfde punt; dan
    zijn alle cycli gesloten. Returns [(laag, hoog, 1.0)].
    """
    r = reversals(x)
    if len(r) < 2:
        return []
    k = max(range(len(r)), key=lambda i: abs(r[i]))
    r = reversals(r[k:] + r[:k] + [r[k]])
    st = []
    out = []
    for p in r:
        st.append(p)
        while len(st) >= 3:
            xr = abs(st[-1] - st[-2])
            yr = abs(st[-2] - st[-3])
            if xr < yr:
                break
            a, b = st[-3], st[-2]
            out.append((min(a, b), max(a, b), 1.0))
            last = st.pop()
            st.pop()
            st.pop()
            st.append(last)
    return out


def count_cycles(x, method):
    """Cycli van reeks x: method 'RAINFLOW' (periodiek) of 'ASTM'."""
    if method == "ASTM":
        return rainflow_astm(x)
    return rainflow_periodic(x)


# ==========================================
# WOHLER-CURVE
# ==========================================
class Curve(object):
    """
    Wohler-curve uit de dict van bridge_fatigue (SI: spanning in Pa,
    amplitude). Vormen:
      "table"   S = [..], N = [..] (S oplopend), interpolatie LOG (log-log,
                DEFI_FONCTION INTERPOL='LOG') of LIN; PROL_GAUCHE (onder de
                tabel) EXCLU (geen schade, zoals Code-Aster's standaard:
                gemeten D = 0), CONSTANT (N van het eerste punt) of
                LINEAIRE (eerste stuk doorgetrokken); PROL_DROITE (boven de
                tabel) LINEAIRE of CONSTANT.
      "power"   N = N_ref (S_ref / S)^m, onder Se geen schade.
      "basquin" S = Sf (2 N)^b, dus N = 0,5 (S / Sf)^(1/b), onder Se geen
                schade.
      "strain"  rek-levensduur (Manson-Coffin-Basquin, S8b): de amplitude
                is een rekamplitude ea = Sf/E (2N)^b + ef (2N)^c; N
                opgelost per cyclus (bisectie in log(2N)).
    correction: GEEN, GOODMAN of GERBER met Su.
    """

    def __init__(self, d):
        self.d = d
        self.form = d["form"]
        self.material = d["material"]
        self.su = d.get("Su")
        self.correction = d.get("correction", "GEEN")
        if self.form == "table":
            self.s = np.asarray(d["S"], dtype=float)
            self.n = np.asarray(d["N"], dtype=float)
            self.log = d.get("interpol", "LOG") == "LOG"
            self.left = d.get("prol_gauche", "EXCLU")
            self.right = d.get("prol_droite", "LINEAIRE")
        self.se = d.get("Se")
        # grootheid van de amplitude: spanning [Pa] of rek [-]
        self.quantity = "strain" if self.form == "strain" else "stress"

    def lowest_damaging(self):
        """Kleinste amplitude met schade (zonder correctie), of 0."""
        if self.form == "table":
            return float(self.s[0]) if self.left == "EXCLU" else 0.0
        return float(self.se or 0.0)

    def cycles(self, s):
        """N(S) per amplitude (np.ndarray); inf = geen schade. Returns
        (N, aantal boven de tabel)."""
        s = np.asarray(s, dtype=float)
        n = np.full(s.shape, np.inf)
        pos = s > 0.0
        above = 0
        if self.form == "table":
            sx, ny = self.s, self.n
            if self.log:
                fx, fy = np.log(sx), np.log(ny)
                xs = np.log(np.where(pos, s, 1.0))
            else:
                fx, fy = sx, ny
                xs = s
            y = np.interp(xs, fx, fy)
            lo = pos & (s < sx[0])
            hi = pos & (s > sx[-1])
            above = int(hi.sum())
            if len(sx) >= 2:
                k0 = (fy[1] - fy[0]) / (fx[1] - fx[0])
                k1 = (fy[-1] - fy[-2]) / (fx[-1] - fx[-2])
            else:
                k0 = k1 = 0.0
            if self.left == "LINEAIRE":
                y = np.where(lo, fy[0] + k0 * (xs - fx[0]), y)
            if self.right == "LINEAIRE":
                y = np.where(hi, fy[-1] + k1 * (xs - fx[-1]), y)
            vals = np.exp(y) if self.log else y
            n = np.where(pos, vals, np.inf)
            if self.left == "EXCLU":
                n = np.where(lo, np.inf, n)
        elif self.form == "power":
            n = np.where(pos, self.d["N"] * (self.d["S"] /
                                             np.where(pos, s, 1.0)) **
                         self.d["m"], np.inf)
        elif self.form == "basquin":
            n = np.where(pos, 0.5 * (np.where(pos, s, 1.0) /
                                     self.d["Sf"]) ** (1.0 / self.d["b"]),
                         np.inf)
        elif self.form == "strain":
            n = np.where(pos, self.strain_life(np.where(pos, s, 1.0)),
                         np.inf)
        if self.se is not None and self.form != "table":
            n = np.where(s < self.se, np.inf, n)
        return n, above

    def strain_life(self, ea):
        """N uit de rekamplitude ea (> 0): ea = Sf/E (2N)^b + ef (2N)^c.
        Rechterlid daalt strikt in x = ln(2N) (b, c < 0): bisectie op
        [ln 1, ln 2e30] tot een interval < 1e-13 (relatief in N). Boven
        de curve bij N = 0,5 (2N = 1) komt N = 0,5 terug (daarna N = 1 in
        damage()); onder 2N = 2e30: N = inf."""
        d = self.d
        ea = np.asarray(ea, dtype=float)
        se, b, ef, c = d["Sf"] / d["E"], d["b"], d["ef"], d["c"]

        def g(x):
            return se * np.exp(b * x) + ef * np.exp(c * x)
        lo = np.zeros(ea.shape)
        hi = np.full(ea.shape, math.log(2.0e30))
        for _k in range(200):
            mid = 0.5 * (lo + hi)
            big = g(mid) > ea
            lo = np.where(big, mid, lo)
            hi = np.where(big, hi, mid)
            if float(np.max(hi - lo)) < 1e-13:
                break
        n = 0.5 * np.exp(0.5 * (lo + hi))
        n = np.where(ea >= g(0.0), 0.5, n)
        return np.where(ea <= g(math.log(2.0e30)), np.inf, n)

    def equivalent_amplitude(self, sa, sm):
        """Amplitude na de gemiddelde-spanningscorrectie (POST_FATIGUE
        CORR_SIGM_MOYE, gemeten in Code-Aster 17.4: GOODMAN S_a/(1 -
        S_m/S_u), GERBER S_a/(1 - (S_m/S_u)^2)). Returns (S_eq, aantal
        cycli met |S_m| >= S_u: daar N = 1)."""
        sa = np.asarray(sa, dtype=float)
        sm = np.asarray(sm, dtype=float)
        if self.correction == "GOODMAN":
            den = 1.0 - sm / self.su
        elif self.correction == "GERBER":
            den = 1.0 - (sm / self.su) ** 2
        else:
            return sa, 0
        bad = den <= 0.0
        return np.where(bad, np.inf, sa / np.where(bad, 1.0, den)), \
            int(bad.sum())

    def damage(self, sa, sm, count):
        """(som n/N, N per cyclus, S_eq per cyclus, boven de tabel,
        gemiddelde >= Su, N < 1)."""
        seq, bad = self.equivalent_amplitude(sa, sm)
        finite = np.isfinite(seq)
        n, above = self.cycles(np.where(finite, seq, 0.0))
        n = np.where(finite, n, 1.0)
        low = int((n < 1.0).sum())
        n = np.maximum(n, 1.0)
        d = float(np.sum(np.asarray(count) / n))
        return d, n, seq, above, bad, low


# ==========================================
# SPANNINGEN EN REKKEN
# ==========================================
def _principal(t):
    """Eigenwaarden per rij van t (n, 6) in volgorde xx, yy, zz, xy, xz,
    yz (tensorcomponenten). Returns (n, 3), oplopend."""
    mats = np.empty((t.shape[0], 3, 3))
    mats[:, 0, 0], mats[:, 1, 1], mats[:, 2, 2] = t[:, 0], t[:, 1], t[:, 2]
    mats[:, 0, 1] = mats[:, 1, 0] = t[:, 3]
    mats[:, 0, 2] = mats[:, 2, 0] = t[:, 4]
    mats[:, 1, 2] = mats[:, 2, 1] = t[:, 5]
    return np.linalg.eigvalsh(mats)


def _largest_abs(ev):
    pick = np.argmax(np.abs(ev), axis=1)
    return ev[np.arange(ev.shape[0]), pick]


def equivalent(sig, kind):
    """Equivalente spanning per rij van sig (n, 6) in volgorde _STRESS:
    VMIS_SG (von Mises x teken van het spoor) of PRIN (hoofdspanning met de
    grootste absolute waarde)."""
    sxx, syy, szz, sxy, sxz, syz = [sig[:, k] for k in range(6)]
    if kind == "PRIN":
        return _largest_abs(_principal(sig))
    vm = np.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 +
                        (szz - sxx) ** 2) +
                 3.0 * (sxy ** 2 + sxz ** 2 + syz ** 2))
    return np.where(sxx + syy + szz >= 0.0, vm, -vm)


def equivalent_strain(eps, kind):
    """Equivalente rek per rij van eps (n, 6) in volgorde _STRAIN
    (tensorcomponenten: EPXY = eps_xy, de helft van de hoekverandering;
    Code-Aster R4.20.01, NIET in de U-doc nagelezen, wel gemeten in de
    verificatie VM-11): PRIN = de hoofdrek met de grootste absolute waarde
    (bij een staaf de rek in de lengterichting); VMIS_SG = sqrt(2/3 e:e)
    van de deviator e, met het teken van het spoor."""
    if kind == "PRIN":
        return _largest_abs(_principal(eps))
    tr = eps[:, 0] + eps[:, 1] + eps[:, 2]
    dev = eps[:, :3] - tr[:, None] / 3.0
    ee = np.sum(dev ** 2, axis=1) + 2.0 * np.sum(eps[:, 3:] ** 2, axis=1)
    eq = np.sqrt(2.0 / 3.0 * ee)
    return np.where(tr >= 0.0, eq, -eq)


_STRAIN = ("EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ")


def _rows(res, field, comps_wanted, index, n):
    """Waarden van een knoopveld in tijdstap index: (waarden (n, k), masker
    (n,)) of None (component ontbreekt). Rij i = Mecway-knoop i+1
    (bridge_mesh-volgorde)."""
    fld = res.getField(field, index)
    sf = fld.toSimpleFieldOnNodes()
    vals, mask = sf.getValues()
    comps = list(sf.getComponents())
    try:
        cols = [comps.index(c) for c in comps_wanted]
    except ValueError:
        return None
    vals = np.asarray(vals)[:n][:, cols]
    mask = np.asarray(mask)[:n][:, cols].all(axis=1)
    return vals, mask


def _stress_rows(res, index, n):
    """SIGM_NOEU van tijdstap index (zie _rows)."""
    return _rows(res, "SIGM_NOEU", _STRESS, index, n)


def _history(res, indexes, n, source, kind):
    """Equivalente grootheid per knoop en tijdstap: (hist (n, stappen),
    masker (n,)). source: "stress" (SIGM_NOEU), "strain" (EPSI_NOEU) of
    "shell" (SIEQ_NOEU van res_sup/res_inf: VMIS_SG of de PRIN_i met de
    grootste absolute waarde; invarianten, dus onafhankelijk van het lokale
    assenstelsel van de schaal)."""
    hist = []
    have = np.ones(n, dtype=bool)
    for idx in indexes:
        if source == "stress":
            got = _rows(res, "SIGM_NOEU", _STRESS, idx, n)
        elif source == "strain":
            got = _rows(res, "EPSI_NOEU", _STRAIN, idx, n)
        else:
            got = _rows(res, "SIEQ_NOEU", (("VMIS_SG",) if kind == "VMIS_SG"
                                           else ("PRIN_1", "PRIN_2",
                                                 "PRIN_3")), idx, n)
        if got is None:
            raise ValueError(tr("%s zonder de benodigde componenten",
                                "%s without the required components") % {
                "stress": "SIGM_NOEU", "strain": "EPSI_NOEU",
                "shell": "SIEQ_NOEU"}[source])
        vals, mask = got
        have &= mask
        vals = np.where(mask[:, None], vals, 0.0)
        if source == "stress":
            hist.append(equivalent(vals, kind))
        elif source == "strain":
            hist.append(equivalent_strain(vals, kind))
        elif kind == "VMIS_SG":
            hist.append(vals[:, 0])
        else:
            hist.append(_largest_abs(vals))
    return np.array(hist).T, have


# ==========================================
# MODEL
# ==========================================
SHELL_SHAPES = frozenset(["tri3", "tri6", "quad4", "quad8", "quad9"])


def read_model(liml_path):
    """Knopen, elsets (naam -> materiaal) en elementen uit de .liml
    (top-level). Returns (coords (n, 3), {elset: materiaal}, {elset:
    [(vorm, [nid])]})."""
    nodes = {}
    mats = {}
    elems = {}
    depth = 0
    with open(liml_path, "rb") as f:
        for event, el in ET.iterparse(f, events=("start", "end")):
            if event == "start":
                depth += 1
                if depth == 2 and el.tag == "solution":
                    break
                continue
            depth -= 1
            if depth != 1:
                continue
            if el.tag == "node":
                nodes[int(el.get("nid"))] = (float(el.get("x")),
                                             float(el.get("y")),
                                             float(el.get("z")))
            elif el.tag == "elset":
                name = el.get("name", "")
                mats[name] = el.get("material", "")
                lst = elems.setdefault(name, [])
                for e in el.findall("elem"):
                    lst.append((e.get("shape"),
                                [int(v) for v in e.get("nodes").split()]))
            el.clear()
    n = max(nodes) if nodes else 0
    coords = np.array([nodes.get(i, (0.0, 0.0, 0.0))
                       for i in range(1, n + 1)], dtype=float)
    return coords, mats, elems


def _node_curves(n, mats, elems, curves, shapes=SOLID_SHAPES):
    """Per curve-index een booleaans masker van de knopen (elementen van
    deze vormen en van de materialen van die curve; '*' = alle materialen
    zonder eigen curve)."""
    named = set(c.material for c in curves if c.material != "*")
    out = []
    for c in curves:
        mask = np.zeros(n, dtype=bool)
        for elset, mat in mats.items():
            if c.material == "*":
                if mat in named:
                    continue
            elif mat != c.material:
                continue
            for shape, nids in elems.get(elset, []):
                if shape in shapes:
                    mask[np.asarray(nids, dtype=int) - 1] = True
        out.append(mask)
    return out


# ==========================================
# HOOFDFUNCTIE (uniaxiaal: spanning of rek, volumes en schalen)
# ==========================================
def compute(res, liml_path, csv_path, curves, options, blocks=None,
            signal=None, shells=None):
    """
    Schade en levensduur per knoop (zie de kop).

    res      : Code-Aster-resultaat met SIGM_NOEU (alle tijdstappen); met
               een rekcurve ook EPSI_NOEU.
    curves   : lijst dicts (bridge_fatigue._parse_curve; zie Curve).
    options  : {"stress": VMIS_SG|PRIN, "strain": PRIN|VMIS_SG,
                "counting": RAINFLOW|ASTM, "repeats": float}
    blocks   : [(amplitudefactor, gemiddelde factor, aantal)] of None
               (statische run: factoren op de statische spanning).
    signal   : [factor, ...] of None (statische run: rainflow op het
               signaal x de statische spanning).
    shells   : {"res_sup": .., "res_inf": ..} (POST_CHAMP EXTR_COQUE met
               SIEQ_NOEU) of None: schalen krijgen schade op de boven- en
               onderzijde (velden damage_top, damage_bottom; damage en
               life = de grootste schade van de twee).
    Zonder blocks en signal: de tijdstappen van res zijn het signaal.

    Returns OrderedDict {"damage": (nids, D), "life": (nids, L), ...} voor
    export_results(node_fields=...), of None bij een fout (regel
    BRIDGE_VERMOEIING_FOUT in de .mess).
    """
    try:
        return _compute(res, liml_path, csv_path, curves, options, blocks,
                        signal, shells)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print("BRIDGE_VERMOEIING_FOUT " + str(e).replace("\n", " ")[:300])
        return None


def _node_damage(h, curve_ids, cv, route, blocks, base_cycles, counting):
    """Schade per herhaling van een knoop met geschiedenis h (equivalente
    grootheid per tijdstap). Returns None (geen cyclus met schade) of
    (d, k, lo, hi, sa, sm, cnt, seq, nn, above, bad, low)."""
    if route == "blocks":
        s0 = h[0]
        sa = np.array([abs(a * s0) for a, _m, _c in blocks])
        sm = np.array([m_ * s0 for _a, m_, _c in blocks])
        cnt = np.array([c for _a, _m, c in blocks])
        lo, hi = sm - sa, sm + sa
    else:
        if route == "signal":
            s0 = h[0]
            cyc = [(lo_ * s0, hi_ * s0, c) if s0 >= 0.0 else
                   (hi_ * s0, lo_ * s0, c)
                   for lo_, hi_, c in base_cycles]
        else:
            # Snel overslaan: geen enkele cyclus kan schade geven
            span = 0.5 * (float(h.max()) - float(h.min()))
            if all(cv[k].correction == "GEEN" and
                   span < cv[k].lowest_damaging() for k in curve_ids):
                return None
            cyc = count_cycles(h, counting)
        if not cyc:
            return None
        lo = np.array([c[0] for c in cyc])
        hi = np.array([c[1] for c in cyc])
        cnt = np.array([c[2] for c in cyc])
        sa = 0.5 * (hi - lo)
        sm = 0.5 * (hi + lo)
    worst = None
    for k in curve_ids:
        d, nn, seq, ab, bd, lw = cv[k].damage(sa, sm, cnt)
        if worst is None or d > worst[0]:
            worst = (d, k, lo, hi, sa, sm, cnt, seq, nn, ab, bd, lw)
    return worst


def _compute(res, liml_path, csv_path, curves, options, blocks, signal,
             shells=None):
    from collections import OrderedDict
    coords, mats, elems = read_model(liml_path)
    n = coords.shape[0]
    cv = [Curve(d) for d in curves]
    masks = _node_curves(n, mats, elems, cv)
    stress_kind = options.get("stress", "VMIS_SG")
    strain_kind = options.get("strain", "PRIN")
    counting = options.get("counting", "RAINFLOW")
    repeats = float(options.get("repeats", 1.0))
    indexes = list(res.getIndexes())
    if blocks or signal:
        indexes = indexes[:1]
    route = "blocks" if blocks else ("signal" if signal else "history")
    base_cycles = count_cycles(signal, counting) if signal else None
    # Bronnen: (laag, resultaat, soort veld, equivalent, curves, maskers)
    sources = []
    for qty, kind in (("stress", stress_kind), ("strain", strain_kind)):
        ids = [k for k, c in enumerate(cv) if c.quantity == qty]
        if ids and any(masks[k].any() for k in ids):
            sources.append(("volume", res, qty, kind, ids, masks))
    shell_masks = None
    if shells:
        shell_masks = _node_curves(n, mats, elems, cv, SHELL_SHAPES)
        ids = [k for k, c in enumerate(cv) if c.quantity == "stress" and
               shell_masks[k].any()]
        for layer, key in (("top", "res_sup"), ("bottom", "res_inf")):
            if ids and shells.get(key) is not None:
                sources.append((layer, shells[key], "shell", stress_kind,
                                ids, shell_masks))
    dmg = np.zeros(n)
    sel = np.zeros(n, dtype=bool)
    layer_dmg = {}
    above = bad = low = 0
    best = (-1.0, None, None, "")
    for layer, r, field, kind, ids, mks in sources:
        hist, have = _history(r, indexes, n, field, kind)
        here = np.zeros(n, dtype=bool)
        for k in ids:
            here |= mks[k]
        here &= have
        sel |= here
        d_layer = np.zeros(n)
        for i in np.nonzero(here)[0]:
            got = _node_damage(hist[i], [k for k in ids if mks[k][i]], cv,
                               route, blocks, base_cycles, counting)
            if got is None:
                continue
            d, k = got[0], got[1]
            above += 1 if got[9] else 0
            bad += 1 if got[10] else 0
            low += 1 if got[11] else 0
            d_layer[i] = d
            if d > best[0]:
                best = (d, i, (k,) + tuple(got[2:9]), layer)
        dmg = np.maximum(dmg, d_layer)
        if layer in ("top", "bottom"):
            layer_dmg[layer] = (here, d_layer)
    nodes = np.nonzero(sel)[0]
    if not len(nodes):
        raise ValueError(tr("geen knopen van volume- of schaalelementen met "
                            "een Wohler-curve en een spanning",
                            "no nodes of volume or shell elements with a "
                            "Wohler curve and a stress"))
    d_rep = dmg[nodes]
    damage = repeats * d_rep
    life = _life(d_rep)
    nids = (nodes + 1).astype("<i4")
    _write_csv(csv_path, best, cv, repeats)
    i = best[1]
    if i is None:
        i = int(nodes[0])
    k = best[2][0] if best[2] else 0
    _print_max(best, repeats, i, coords, cv[k].material,
               ("%g" % float(np.sum(best[2][5]))) if best[2] else 0)
    print("BRIDGE_VERMOEIING_INFO NODES={0} FAILED={1} ROUTE={2} STEPS={3} "
          "COUNTING={4} STRESS={5} REPEATS={6!r} LAYER={7}".format(
              len(nodes), int((damage >= 1.0).sum()), route, len(indexes),
              counting, ("%s/%s" % (stress_kind, strain_kind)
                         if any(c.quantity == "strain" for c in cv)
                         else stress_kind), repeats, best[3] or "volume"))
    if above or bad or low:
        print("BRIDGE_VERMOEIING_WARN ABOVE={0} MEAN_SU={1} N_BELOW_1={2}"
              .format(above, bad, low))
    print("BRIDGE_VERMOEIING_CSV " + csv_path)
    out = OrderedDict([("damage", (nids, damage)), ("life", (nids, life))])
    for layer in ("top", "bottom"):
        if layer in layer_dmg:
            here, d_layer = layer_dmg[layer]
            ln = np.nonzero(here)[0]
            out["damage_" + layer] = ((ln + 1).astype("<i4"),
                                      repeats * d_layer[ln])
    return out


def _life(d_rep):
    """Levensduur in herhalingen: 1 / D per herhaling, begrensd op
    LIFE_MAX (ook bij D = 0)."""
    with np.errstate(divide="ignore"):
        life = np.where(d_rep > 0.0, 1.0 / np.where(d_rep > 0.0, d_rep, 1.0),
                        LIFE_MAX)
    return np.minimum(life, LIFE_MAX)


def _print_max(best, repeats, i, coords, material, cycles):
    print("BRIDGE_VERMOEIING_MAX D={0!r} LIFE={1!r} NODE={2} X={3!r} "
          "Y={4!r} Z={5!r} CYCLES={6} MATERIAL={7}".format(
              float(repeats * max(best[0], 0.0)),
              float(min(1.0 / best[0], LIFE_MAX) if best[0] > 0.0
                    else LIFE_MAX),
              i + 1, float(coords[i, 0]), float(coords[i, 1]),
              float(coords[i, 2]), cycles, material))


def _write_csv(path, best, curves, repeats):
    """Cyclustabel van de knoop met de grootste schade (laag: volume, top
    of bottom; bij een rekcurve zijn s_* rekken)."""
    head = ("node,material,layer,cycle,s_min,s_max,amplitude,mean,count,"
            "s_eq,n_failure,damage")
    lines = [head]
    if best[2] is not None:
        k, lo, hi, sa, sm, cnt, seq, nn = best[2]
        mat = curves[k].material
        for j in range(len(lo)):
            lines.append(",".join([
                str(best[1] + 1), mat.replace(",", ";"), best[3] or "volume",
                str(j + 1)] +
                ["%.10g" % float(v) for v in (lo[j], hi[j], sa[j], sm[j],
                                              cnt[j], seq[j], nn[j])] +
                ["%.10g" % float(repeats * cnt[j] / nn[j])]))
    with io.open(path, "w", encoding="utf-8") as f:
        f.write(u"\n".join(lines) + u"\n")


# ==========================================
# MULTIAXIAAL (CALC_FATIGUE FATIGUE_MULTI, S8b)
# ==========================================
def _cham_values(cham, n, comps_wanted):
    """{component: (n,) array} en masker (n,) van een cham_no van
    CALC_FATIGUE (DOMA_NOEUD)."""
    sf = cham.toSimpleFieldOnNodes()
    vals, mask = sf.getValues()
    comps = list(sf.getComponents())
    vals = np.asarray(vals)[:n]
    mask = np.asarray(mask)[:n]
    out = {}
    ok = np.ones(vals.shape[0], dtype=bool)
    for c in comps_wanted:
        j = comps.index(c)
        out[c] = vals[:, j]
        ok &= mask[:, j]
    return out, ok


def compute_multi(chams, liml_path, csv_path, curves, options):
    """
    Schade en levensduur per knoop uit CALC_FATIGUE TYPE_CALCUL=
    'FATIGUE_MULTI' OPTION='DOMA_NOEUD' (kritiek vlak in Code-Aster).

    chams   : [(cham_no, aantal)]: een veld per resultaat (de tijdstappen
              van de run of het signaal: aantal 1; een blok: aantal n).
    curves  : dicts van bridge_fatigue (Wohler: tabel, macht, Basquin).
    options : #ASTER-VERMOEIING (criterion, loading, a, b, c, repeats).

    PERIODIQUE: de equivalente spanning SIGEQ = max(SIGEQ1, SIGEQ2) (de
    twee vlakken met de grootste schuifamplitude; gemeten sslv135a: SIGEQ =
    c (dtau/2 + a max(P_max, 0)) bij DANG_VAN_MODI_AC, met N_max bij
    MATAKE_MODI_AC; D_VAN_B/MATAKE_B tellen niet mee) is de amplitude op de
    Wohler-curve van de bridge (Curve.cycles, dezelfde als uniaxiaal): N per
    resultaat, D = herhalingen x sum(aantal / N). Code-Aster's eigen NBRUP1
    is afgerond op een geheel getal (gemeten 10946,0 bij N = 10946,1).
    NON_PERIODIQUE: D = herhalingen x ENDO1 van Code-Aster (rainflow op het
    kritieke vlak, met de Wohler-functie van bridge_fatigue.wohler_points).
    Extra velden: fatigue_eq (SIGEQ, periodiek) en fatigue_ratio = SIGEQ /
    (c b) (DANG_VAN, MATAKE, CROSSLAND) of SIGEQ / b (FORMULE) als b
    gegeven is.
    """
    try:
        return _compute_multi(chams, liml_path, csv_path, curves, options)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print("BRIDGE_VERMOEIING_FOUT " + str(e).replace("\n", " ")[:300])
        return None


def _compute_multi(chams, liml_path, csv_path, curves, options):
    from collections import OrderedDict
    coords, mats, elems = read_model(liml_path)
    n = coords.shape[0]
    cv = [Curve(d) for d in curves]
    masks = _node_curves(n, mats, elems, cv)
    repeats = float(options.get("repeats", 1.0))
    periodic = options.get("loading", "PERIODIQUE") == "PERIODIQUE"
    crit = options.get("criterion")
    sel = np.zeros(n, dtype=bool)
    for m in masks:
        sel |= m
    eqs, endo, counts = [], [], []
    for cham, cnt in chams:
        if periodic:
            v, ok = _cham_values(cham, n, ("SIGEQ1", "SIGEQ2"))
            eqs.append(np.maximum(v["SIGEQ1"], v["SIGEQ2"]))
        else:
            v, ok = _cham_values(cham, n, ("ENDO1",))
            # onder de Wohler-tabel (EXCLU): NBRUP = 0, ENDO ~ 1e-308
            # (gemeten) = geen schade
            endo.append(np.where(v["ENDO1"] > 1e-300, v["ENDO1"], 0.0))
        counts.append(float(cnt))
        sel &= ok
    nodes = np.nonzero(sel)[0]
    if not len(nodes):
        raise ValueError(tr("geen knopen van volume-elementen met een "
                            "Wohler-curve in het CALC_FATIGUE-veld",
                            "no nodes of volume elements with a Wohler curve "
                            "in the CALC_FATIGUE field"))
    counts = np.array(counts)
    above = low = 0
    if periodic:
        eq = np.array(eqs)                 # (resultaten, n)
        dk_all = np.full((len(cv), n), -1.0)
        nn_all = []
        for k, c in enumerate(cv):
            nn, _ab = c.cycles(np.maximum(eq, 0.0).ravel())
            nn = nn.reshape(eq.shape)
            low_k = nn < 1.0
            nn = np.maximum(nn, 1.0)
            nn_all.append(nn)
            dk_all[k] = np.where(masks[k], np.sum(counts[:, None] / nn,
                                                  axis=0), -1.0)
            if c.form == "table":
                above += int(np.sum(masks[k] & (eq > c.s[-1]).any(axis=0)))
            low += int(np.sum(masks[k] & low_k.any(axis=0)))
        used = np.argmax(dk_all, axis=0)
        dmg = np.maximum(dk_all[used, np.arange(n)], 0.0)
    else:
        used = np.argmax(np.array(masks), axis=0)
        dmg = np.sum(counts[:, None] * np.array(endo), axis=0)
    d_rep = dmg[nodes]
    damage = repeats * d_rep
    life = _life(d_rep)
    nids = (nodes + 1).astype("<i4")
    i = int(nodes[np.argmax(d_rep)])
    best = float(dmg[i])
    material = cv[int(used[i])].material
    _print_max((best,), repeats, i, coords, material,
               "%g" % float(np.sum(counts)))
    print("BRIDGE_VERMOEIING_INFO NODES={0} FAILED={1} ROUTE={2} STEPS={3} "
          "COUNTING={4} STRESS={5} REPEATS={6!r} LAYER=volume".format(
              len(nodes), int((damage >= 1.0).sum()),
              options.get("route", "history"), len(chams),
              options.get("loading"), crit, repeats))
    out = OrderedDict([("damage", (nids, damage)), ("life", (nids, life))])
    lines = ["node,material,layer,cycle,s_min,s_max,amplitude,mean,count,"
             "s_eq,n_failure,damage"]
    if periodic:
        eq_max = np.max(eq, axis=0)
        out["fatigue_eq"] = (nids, eq_max[nodes])
        line = "BRIDGE_VERMOEIING_MULTI CRIT={0} EQ={1!r}".format(
            crit, float(eq_max[i]))
        b = options.get("b")
        if b:
            scale = b * (1.0 if crit == "FORMULE" else
                         float(options.get("c", 1.0)))
            ratio = eq_max / scale
            out["fatigue_ratio"] = (nids, ratio[nodes])
            line += " RATIO={0!r} RATIO_MAX={1!r}".format(
                float(ratio[i]), float(np.max(ratio[nodes])))
        print(line)
        per_n = nn_all[int(used[i])][:, i]
        for j in range(eq.shape[0]):
            lines.append(",".join([str(i + 1), material.replace(",", ";"),
                                   "volume", str(j + 1), "nan", "nan",
                                   "%.10g" % float(eq[j, i]), "nan",
                                   "%.10g" % counts[j],
                                   "%.10g" % float(eq[j, i]),
                                   "%.10g" % float(per_n[j]),
                                   "%.10g" % float(repeats * counts[j] /
                                                   per_n[j])]))
    else:
        print("BRIDGE_VERMOEIING_MULTI CRIT={0}".format(crit))
        en = np.array(endo)
        for j in range(en.shape[0]):
            dj = float(en[j, i])
            lines.append(",".join([str(i + 1), material.replace(",", ";"),
                                   "volume", str(j + 1), "nan", "nan", "nan",
                                   "nan", "%.10g" % counts[j], "nan",
                                   "%.10g" % (1.0 / dj if dj > 0.0 else
                                              float("inf")),
                                   "%.10g" % (repeats * counts[j] * dj)]))
    with io.open(csv_path, "w", encoding="utf-8") as f:
        f.write(u"\n".join(lines) + u"\n")
    if above or low:
        print("BRIDGE_VERMOEIING_WARN ABOVE={0} MEAN_SU=0 N_BELOW_1={1}"
              .format(above, low))
    print("BRIDGE_VERMOEIING_CSV " + csv_path)
    return out
