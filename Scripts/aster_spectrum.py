# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# aster_spectrum.py
#
# Verantwoordelijkheid: de stappen van de responsspectrumanalyse
# (bridge_spectrum.py) die BINNEN de Code-Aster-job draaien (CPython 3.11
# met numpy, zoals aster_results_export.py; NIET in IronPython):
#   table()     spectrumwaarde uit een Mecway-tabel (lineair, constant
#               buiten de tabel), gebruikt in de FORMULE per richting;
#   report()    regels BRIDGE_SPECTRUM_* in de .mess: spectrale waarde per
#               mode, effectieve massa per richting en de gecombineerde
#               totale reactiekracht (base shear);
#   maxima()    de maxima van COMB_SISM_MODAL als een resultaat met een
#               index (CREA_RESU DYNA_TRANS: DEPL, ACCE, REAC_NODA,
#               SIGM_ELNO), eventueel met de richtingscombinatie 100-30-30;
#   von_mises() von Mises per knoop, per mode berekend en gecombineerd
#               (zie daar), als extra veld SIEQ_NOEU voor
#               aster_results_export.
#
# De combinatieregels zijn die van Code-Aster zelf: CombModalResponse en
# BaseRunner._comb_directions uit code_aster/MacroCommands/
# comb_sism_modal_ops.py (17.4), zodat de totale reactie en von Mises met
# dezelfde formules gecombineerd worden als de velden van COMB_SISM_MODAL.
#
# Gebruik geen f-strings: zelfde stijl als de rest van de bridge.

import math

import numpy as np

from aster_results_export import tr   # taal van de meldingen

# Componenten van SIGM_NOEU (volgorde voor von Mises).
_SIG = ("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ")
# COMB_SISM_MODAL OPTION -> NOM_CHAM in het resultaat voor Mecway.
_FIELDS = (("DEPL", "DEPL"), ("ACCE_ABSOLU", "ACCE"),
           ("REAC_NODA", "REAC_NODA"), ("SIGM_ELNO", "SIGM_ELNO"))


def table(x, xs, ys):
    """Waarde van een Mecway-tabel in x: lineair tussen de punten, constant
    buiten de tabel (Mecway manual v33 hoofdstuk 9, tabellen)."""
    return float(np.interp(x, xs, ys))


# ==========================================
# MODALE GROOTHEDEN
# ==========================================
def _modal(modes, amors):
    """(frequenties, {richting: participatiefactoren}, {richting: fractie
    effectieve massa}) in de volgorde van NUME_ORDRE (CALC_MODES
    PLUS_PETITE: oplopend, zoals COMB_SISM_MODAL ze sorteert)."""
    para = modes.LIST_PARA()
    freqs = np.array(para["FREQ"], dtype=float)
    gam, frac = {}, {}
    for d in "XYZ":
        gam[d] = np.array(para["FACT_PARTICI_D" + d], dtype=float)
        key = "MASS_EFFE_UN_D" + d
        frac[d] = np.array(para[key], dtype=float) if key in para else None
    return freqs, gam, frac


def _spectral(spec, amors, freqs):
    """(S, w_r) per mode zoals BaseRunner.compute (MONO_APPUI): de
    spectrumwaarde maal ECHELLE, met CORR_FREQ en NATURE omgerekend naar
    een versnelling, en de (eventueel gedempte) hoekfrequentie."""
    s = np.array([spec["spec"](a, f) * spec["scale"]
                  for a, f in zip(amors, freqs)], dtype=float)
    w = 2.0 * np.pi * freqs
    if spec.get("corr_freq") == "OUI":
        c = np.sqrt(1.0 - amors ** 2)
        w = w * c
        s = s * c
    if spec.get("nature") == "DEPL":
        s = (2.0 * np.pi * freqs) ** 2 * s
    elif spec.get("nature") == "VITE":
        s = 2.0 * np.pi * freqs * s
    return s, w


def _s_cut(spec, amor, f_cut):
    """Spectrumwaarde op de afsnijfrequentie (BaseRunner._s_r_freq_cut)."""
    s = spec["spec"](amor, f_cut) * spec["scale"]
    if spec.get("corr_freq") == "OUI":
        s = s * math.sqrt(1.0 - amor ** 2)
    if spec.get("nature") == "DEPL":
        s = (2.0 * math.pi * f_cut) ** 2 * s
    elif spec.get("nature") == "VITE":
        s = 2.0 * math.pi * f_cut * s
    return s


def _comb_modes(rule, amors, freqs):
    from code_aster.MacroCommands.comb_sism_modal_ops import (
        CombModalResponse)
    return CombModalResponse({"TYPE": rule}, "MONO_APPUI", np.asarray(amors),
                             np.asarray(freqs))


def comb_directions(rule, arrays):
    """Richtingscombinatie van de maxima per richting (arrays >= 0):
    'QUAD' (SRSS), 'ABS', 'NEWMARK' (100-40-40, Code-Aster) zoals
    BaseRunner._comb_directions, of '100-30-30' (EN 1998-1 par.
    4.3.3.5.2 (4): E_a + 0,30 E_b + 0,30 E_c, de grootste over de
    richtingen a)."""
    arrays = [np.asarray(a, dtype=float) for a in arrays]
    if len(arrays) == 1:
        return arrays[0]
    if rule == "100-30-30":
        total = sum(arrays)
        return np.max([a + 0.3 * (total - a) for a in arrays], axis=0)
    from code_aster.MacroCommands.comb_sism_modal_ops import BaseRunner
    return BaseRunner._comb_directions(rule, arrays)[0]


def _node_values(field, comps):
    """(waarden n x k, masker n) van een knoopveld voor de componenten
    comps (ontbrekende componenten 0)."""
    sf = field.toSimpleFieldOnNodes()
    vals, mask = sf.getValues()
    names = list(sf.getComponents())
    out = np.zeros((vals.shape[0], len(comps)))
    have = np.zeros(vals.shape[0], dtype=bool)
    for k, c in enumerate(comps):
        if c in names:
            j = names.index(c)
            out[:, k] = np.where(mask[:, j], vals[:, j], 0.0)
            have |= mask[:, j]
    return out, have


def _pseudo_index(psmod, d):
    """NUME_ORDRE van de pseudo-mode van richting d (NOEUD_CMP 'ACCE    X',
    zoals BaseRunner._corr_pseudo_mode)."""
    acc = psmod.getAccessParameters()
    for nr, nc in zip(acc["NUME_MODE"], acc["NOEUD_CMP"]):
        if str(nc).split() == ["ACCE", d]:
            return nr
    return None


# ==========================================
# OVERZICHT IN DE .MESS
# ==========================================
def report(modes, specs, amors, comb_mode, comb_dir, psmod=None,
           f_cut=None):
    """
    BRIDGE_SPECTRUM-regels (gelezen door bridge_spectrum.report_summary):
      BRIDGE_SPECTRUM_MODE <nr> FREQ=.. XI=.. SA_X=.. [SA_Y=..]
      BRIDGE_SPECTRUM_MASSA X=<som fractie> ..
      BRIDGE_SPECTRUM_CUT FREQ=.. SA_X=..      (missing mass)
      BRIDGE_SPECTRUM_REACTIE FX=.. FY=.. FZ=..
    De totale reactie is de som van REAC_NODA over alle knopen, per mode
    (en voor de pseudo-mode) opgeteld en daarna gecombineerd zoals
    COMB_SISM_MODAL een knoopveld combineert; de som van de gecombineerde
    knoopreacties zelf is groter (maxima zonder teken).
    """
    try:
        amors = np.asarray(amors, dtype=float)
        freqs, gam, frac = _modal(modes, amors)
        sa = dict((s["dir"], _spectral(s, amors, freqs)) for s in specs)
        for i, f in enumerate(freqs):
            parts = ["BRIDGE_SPECTRUM_MODE {0} FREQ={1:.6g} XI={2:.6g}".format(
                i + 1, f, amors[i])]
            parts += ["SA_{0}={1:.6g}".format(d, sa[d][0][i]) for d in sa]
            print(" ".join(parts))
        print("BRIDGE_SPECTRUM_MASSA " + " ".join(
            "{0}={1:.6g}".format(d, float(np.sum(frac[d])))
            for d in sa if frac[d] is not None))
        if psmod is not None and f_cut is not None:
            print("BRIDGE_SPECTRUM_CUT FREQ={0:.6g} ".format(f_cut) + " ".join(
                "SA_{0}={1:.6g}".format(s["dir"], _s_cut(s, amors[-1], f_cut))
                for s in specs))
        # Totale reactie per mode: som van REAC_NODA over de knopen
        idx = list(modes.getIndexes())
        r_modes = np.array([_node_values(modes.getField("REAC_NODA", k),
                                         ("DX", "DY", "DZ"))[0].sum(axis=0)
                            for k in idx])
        comb = _comb_modes(comb_mode, amors, freqs)
        per_dir = []
        for s in specs:
            d = s["dir"]
            S, w = sa[d]
            r_mi = (S * gam[d] / w ** 2)[:, None] * r_modes
            r_m2, _qs = comb.get(r_mi)
            r_c = np.zeros(3)
            if psmod is not None and f_cut is not None:
                k = _pseudo_index(psmod, d)
                if k is not None:
                    r_ps = _node_values(psmod.getField("REAC_NODA", k),
                                        ("DX", "DY", "DZ"))[0].sum(axis=0)
                    pr = np.sum((gam[d] / (2.0 * np.pi * freqs) ** 2)[:, None]
                                * r_modes, axis=0)
                    r_c = (r_ps - pr) * _s_cut(s, amors[-1], f_cut)
            per_dir.append(np.sqrt(r_m2 + r_c ** 2))
        total = comb_directions(comb_dir, per_dir)
        print("BRIDGE_SPECTRUM_REACTIE FX={0:.6g} FY={1:.6g} FZ={2:.6g}"
              .format(*[float(v) for v in total]))
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(tr("aster_spectrum: WAARSCHUWING: overzicht niet volledig: "
                 "{0}",
                 "aster_spectrum: WARNING: summary not complete: {0}").format(
            e))


# ==========================================
# MAXIMA ALS RESULTAAT
# ==========================================
def maxima(sism, options, comb_dir, dirs, **affe):
    """
    CREA_RESU (TYPE_RESU 'DYNA_TRANS', een index) met de gecombineerde
    maxima van COMB_SISM_MODAL (NOM_CAS 'TOTA'); bij comb_dir '100-30-30'
    combineert deze functie de richtingen zelf uit de NOM_CAS 'DIRE_X' ..
    (TYPE_RESU VALE_DIRE). ACCE_ABSOLU heet in het resultaat ACCE. affe:
    MODELE, CHAM_MATER, CARA_ELEM voor AFFE.
    """
    from code_aster.Commands import CREA_RESU
    acc = sism.getAccessParameters()
    cases = dict(zip(acc["NOM_CAS"], acc["NUME_ORDRE"]))
    res = None
    for opt, cham in _FIELDS:
        if opt not in options:
            continue
        fld = sism.getField(opt, cases["TOTA"]).copy()
        if comb_dir == "100-30-30" and len(dirs) > 1:
            parts = [np.asarray(sism.getField(opt, cases["DIRE_" + d])
                                .getValues()) for d in dirs]
            fld.setValues(comb_directions("100-30-30", parts))
        kw = dict(OPERATION="AFFE", TYPE_RESU="DYNA_TRANS", NOM_CHAM=cham,
                  AFFE=dict(CHAM_GD=fld, INST=0.0, **affe))
        if res is None:
            res = CREA_RESU(**kw)
        else:
            res = CREA_RESU(reuse=res, RESULTAT=res, **kw)
    return res


class _NodeValues(object):
    """Minimaal knoopveld voor aster_results_export (_Fields.get):
    toSimpleFieldOnNodes().getValues() / getComponents()."""

    def __init__(self, comps, vals, mask):
        self.comps, self.vals, self.mask = list(comps), vals, mask

    def toSimpleFieldOnNodes(self):
        return self

    def getValues(self):
        return self.vals, self.mask

    def getComponents(self):
        return self.comps


class Result(object):
    """Resultaat voor aster_results_export: res (CREA_RESU) plus eigen
    velden (bijv. von Mises) als extra NOM_CHAM."""

    def __init__(self, res, extra):
        self._res, self._extra = res, dict(extra)

    def getField(self, name, index):
        if name in self._extra:
            return self._extra[name]
        return self._res.getField(name, index)

    def getFieldsNames(self):
        return list(self._res.getFieldsNames()) + list(self._extra)

    def __getattr__(self, name):
        return getattr(self._res, name)


def _vm(sig):
    """von Mises uit SIXX SIYY SIZZ SIXY SIXZ SIYZ (rijen)."""
    sxx, syy, szz, sxy, sxz, syz = [sig[:, k] for k in range(6)]
    return np.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 +
                          (szz - sxx) ** 2) +
                   3.0 * (sxy ** 2 + sxz ** 2 + syz ** 2))


def von_mises(modes, specs, amors, comb_mode, comb_dir, psmod=None,
              f_cut=None):
    """
    von Mises per knoop (veld SIEQ_NOEU, component VMIS) uit SIGM_NOEU van
    de modes (en de pseudo-modes), of None. Per richting: de modale
    bijdragen R_i = |S_i Gamma_i / w_i^2| VM(sigma_i) met de modale regel
    gecombineerd, plus de missing mass VM((sigma_ps - sum Gamma_i/w_i^2
    sigma_i) S_cut) kwadratisch (zoals COMB_SISM_MODAL R_x = sqrt(R_m^2 +
    R_c^2)); daarna de richtingsregel. Met SRSS is dat de verwachtingswaarde
    van sigma^T P sigma (von Mises is een kwadratische vorm, modes
    ongecorreleerd); met CQC een bovengrens (|sigma_i^T P sigma_j| <=
    VM_i VM_j, rho_ij >= 0). Von Mises uit de gecombineerde componenten
    (maxima zonder teken, elk op een ander moment) is dat niet.
    """
    try:
        amors = np.asarray(amors, dtype=float)
        freqs, gam, _frac = _modal(modes, amors)
        idx = list(modes.getIndexes())
        sig, have = [], None
        for k in idx:
            v, h = _node_values(modes.getField("SIGM_NOEU", k), _SIG)
            sig.append(v)
            have = h if have is None else (have | h)
        vm_modes = np.array([_vm(s) for s in sig])
        comb = _comb_modes(comb_mode, amors, freqs)
        per_dir = []
        for s in specs:
            d = s["dir"]
            S, w = _spectral(s, amors, freqs)
            r_mi = np.abs(S * gam[d] / w ** 2)[:, None] * vm_modes
            r_m2, _qs = comb.get(r_mi)
            r_c2 = 0.0
            if psmod is not None and f_cut is not None:
                k = _pseudo_index(psmod, d)
                if k is not None:
                    resid = _node_values(psmod.getField("SIGM_NOEU", k),
                                         _SIG)[0]
                    for i, si in enumerate(sig):
                        resid = resid - gam[d][i] / (
                            2.0 * np.pi * freqs[i]) ** 2 * si
                    r_c2 = (_vm(resid) * _s_cut(s, amors[-1], f_cut)) ** 2
            per_dir.append(np.sqrt(r_m2 + r_c2))
        total = comb_directions(comb_dir, per_dir)
        return _NodeValues(["VMIS"], total[:, None], have[:, None])
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(tr("aster_spectrum: WAARSCHUWING: von Mises niet berekend: "
                 "{0}",
                 "aster_spectrum: WARNING: von Mises not computed: {0}")
              .format(e))
        return None
