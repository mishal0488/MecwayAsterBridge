# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# aster_acoustic_export.py
#
# Verantwoordelijkheid: het harmonische akoestische resultaat (DYNA_VIBRA
# TYPE_CALCUL='HARM' op PHENOMENE='ACOUSTIQUE', bridge_acoustic.
# build_acoustic_harm) naar Mecway's resultaatformaat, plus de
# frequentierespons op de microfoonknopen (CSV en BRIDGE_MIC-regels in de
# .mess). Draait BINNEN de Code-Aster-job (CPython 3.11 met numpy, zoals
# aster_results_export.py), NIET in IronPython.
#
# Code-Aster 17.4 (gezien in een run): DYNA_VIBRA geeft een resultaat van
# het type ACOU_HARMO (FullHarmonicAcousticResult) met per frequentie
# (NUME_ORDRE, toegangsparameter FREQ) het complexe veld PRES
# (FieldOnNodesComplex, component PRES). Houd de tussenobjecten vast:
# res.getField(..).toSimpleFieldOnNodes().getValues() in een keten gaf een
# keer een <EXCEPTION> DVP_2 (floating point exception), los niet.
#
# Uitvoer per frequentie een blok
#     <mode id="k" frequency="f Hz">
# (aster_results_export._write_fragment, zoals de modes: Mecway toont ze als
# lijst met de frequentie) met de knoopvariabelen
#     pressure       Re(P)   [Pa]  (Mecway's naam voor de druk; met
#     pressurei      Im(P)   [Pa]   pressurei zoals Mecway's complexe modes)
#     pressuremag    |P|     [Pa]  amplitude
#     pressurephase  arg(P)  [graden], -180..180
#     spl            20 log10(|P| / (sqrt(2) p_ref)) [dB], p_ref = 20 uPa
# met p(t) = Re(P e^(i w t)) (Code-Aster; zelfde conventie als Mecway's
# complexe modes, manual v33 par. 6.9.2). Het geluidsdrukniveau is op de
# effectieve waarde |P|/sqrt(2) gebaseerd; Code-Aster's PRAC_ELNO DB
# (ahlv100a) gebruikt |P| zelf en ligt dus 20 log10(sqrt(2)) = 3,0103 dB
# hoger (nagerekend: |P| 6,2426 Pa -> DB 109,886).
#
# Gebruik geen f-strings: zelfde stijl als de rest van de bridge.

import csv
import math
import os
import time
import traceback
from collections import OrderedDict

import numpy as np

import aster_results_export as are

P_REF = 2.0e-5           # Pa, referentiedruk lucht (ISO 1683)
_SQRT2 = math.sqrt(2.0)
# Ondergrens voor de log in het geluidsdrukniveau (een drukknoop heeft
# p = 0): -200 dB i.p.v. -inf, zodat er geen inf/NaN in de .liml komt.
_SPL_FLOOR = P_REF * 1.0e-10

_CSV_HEAD = {
    "nl": ("microfoon", "frequentie_Hz", "p_re_Pa", "p_im_Pa", "p_amp_Pa",
           "fase_graden", "spl_dB", "knopen"),
    "en": ("microphone", "frequency_Hz", "p_re_Pa", "p_im_Pa", "p_amp_Pa",
           "phase_deg", "spl_dB", "nodes"),
}


def spl(amplitude):
    """Geluidsdrukniveau [dB re 20 uPa] van een amplitude |P| (numpy of
    float): 20 log10(|P| / (sqrt(2) p_ref))."""
    a = np.maximum(np.asarray(amplitude, dtype=float) / _SQRT2, _SPL_FLOOR)
    return 20.0 * np.log10(a / P_REF)


class _Harm(object):
    """Complexe druk per frequentie (rij i = Mecway-knoop i+1); dezelfde
    rol als aster_results_export._Fields voor _write_fragment."""

    def __init__(self, pressures):
        self.pressures = pressures
        self.index = None

    def set_index(self, index):
        self.index = index

    def current(self):
        return self.pressures[self.index]


def _getter(rows, fn):
    def getter(data):
        return fn(data.current()[rows])
    return getter


_VARIABLES = (
    ("pressure", lambda p: p.real),
    ("pressurei", lambda p: p.imag),
    ("pressuremag", lambda p: np.abs(p)),
    ("pressurephase", lambda p: np.degrees(np.angle(p))),
    ("spl", lambda p: spl(np.abs(p))),
)


def _pressures(res, n):
    """({index: complexe druk (n,)}, rijen met een waarde, {index: f})."""
    indexes = list(res.getIndexes())
    acc = res.getAccessParameters()
    freqs = dict(zip(acc["NUME_ORDRE"], acc["FREQ"]))
    out = OrderedDict()
    rows = None
    for idx in indexes:
        fld = res.getField("PRES", idx)
        sf = fld.toSimpleFieldOnNodes()
        vals, mask = sf.getValues()
        col = list(sf.getComponents()).index("PRES")
        p = np.zeros(n, dtype=complex)
        m = np.asarray(mask)[:n, col]
        p[m] = np.asarray(vals)[:n, col][m]
        out[idx] = p
        if rows is None:
            rows = np.nonzero(m)[0]
    return out, rows, freqs


def microphone_rows(pressures, freqs, microphones):
    """
    [(naam, f, P, aantal knopen)] per microfoon en frequentie. P is de
    complexe druk van een knoop; bij meer knopen het gemiddelde van P met
    als amplitude de kwadratisch gemiddelde |P| (sqrt(mean |P|^2), de
    energie over de knopen; zo geeft een selectie op een plaats met
    tegengestelde fasen geen schijnbare stilte).
    """
    out = []
    for name, nids in (microphones or {}).items():
        rows = [int(k) - 1 for k in nids]
        for idx, p in pressures.items():
            vals = p[rows]
            mean = complex(vals.mean())
            amp = float(np.sqrt(np.mean(np.abs(vals) ** 2)))
            if abs(mean) > 0.0:
                value = mean / abs(mean) * amp
            else:
                value = complex(amp, 0.0)
            out.append((name, float(freqs.get(idx, idx)), value, len(rows)))
    return out


def export_harmonic(res, out_path, liml_path, microphones=None,
                    csv_path=None, language="en", precision="double"):
    """
    Schrijf het harmonische drukveld als Mecway-<results>-fragment (een
    <mode> per frequentie) en de microfoonrespons.

    microphones : {naam: [Mecway-knoopnummers]} (bridge_acoustic,
                  selecties MICROFOON_...), of None.
    csv_path    : CSV met de microfoonrespons (komma's, punt als
                  decimaalteken, zoals _vergelijking.csv), of None.
    Returns dict met statistiek, of None bij een fout (melding in de
    .mess).
    """
    t0 = time.time()
    try:
        return _export(res, out_path, liml_path, microphones, csv_path,
                       language, precision, t0)
    except Exception as e:
        traceback.print_exc()
        print(are.tr("aster_acoustic_export: FOUT: {0}",
                     "aster_acoustic_export: ERROR: {0}").format(e))
        return None


def _export(res, out_path, liml_path, microphones, csv_path, language,
            precision, t0):
    model = are.read_liml_model(liml_path)
    n = model["coords"].shape[0]
    err = are._check_node_order(res.getMesh(), model["coords"])
    if err:
        print(are.tr("aster_acoustic_export: FOUT: {0}",
                     "aster_acoustic_export: ERROR: {0}").format(err))
        return None
    pressures, rows, freqs = _pressures(res, n)
    if not pressures or rows is None or not len(rows):
        print(are.tr("aster_acoustic_export: FOUT: geen drukveld in het "
                     "resultaat.",
                     "aster_acoustic_export: ERROR: no pressure field in the "
                     "result."))
        return None
    ids = (rows + 1).astype("<i4")
    variables = [are._Var(name, False, ids, _getter(rows, fn))
                 for name, fn in _VARIABLES]
    elset_nodes, free_nodes = are._group_nodes(model, n)
    for var in variables:
        are._assign_groups(var, model, elset_nodes, free_nodes)
    data = _Harm(pressures)
    indexes = list(pressures.keys())
    tmp = out_path + ".tmp"
    with open(tmp, "wb") as f:
        are._write_fragment(f, variables, data, indexes, None, precision,
                            model, ("frequency", freqs))
        size = f.tell()
    os.replace(tmp, out_path)
    mic = microphone_rows(pressures, freqs, microphones)
    _print_microphones(mic)
    if csv_path and mic:
        _write_csv(csv_path, mic, language)
    print("BRIDGE_HARM FREQS={0} FMIN={1} FMAX={2} MICROFOONS={3}".format(
        len(indexes), are._num(min(freqs.values())),
        are._num(max(freqs.values())), len(microphones or {})))
    print("aster_acoustic_export: {0} geschreven: {1} frequentie(s), "
          "{2:.1f} MB in {3:.2f} s.".format(out_path, len(indexes),
                                            size / 1.0e6, time.time() - t0))
    return {"stappen": len(indexes), "bytes": size}


def _print_microphones(mic):
    """BRIDGE_MIC-regels in de .mess (naam als laatste: mag spaties
    hebben)."""
    for name, f, p, count in mic:
        print("BRIDGE_MIC F={0} RE={1} IM={2} AMP={3} FASE={4} SPL={5} "
              "N={6} NAAM={7}".format(
                  are._num(f), are._num(p.real), are._num(p.imag),
                  are._num(abs(p)), "%.6f" % math.degrees(math.atan2(
                      p.imag, p.real)), "%.6f" % float(spl(abs(p))), count,
                  name))


def _write_csv(path, mic, language):
    head = _CSV_HEAD.get(language, _CSV_HEAD["en"])
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(head)
        for name, freq, p, count in mic:
            w.writerow([name, are._num(freq), are._num(p.real),
                        are._num(p.imag), are._num(abs(p)),
                        "%.6f" % math.degrees(math.atan2(p.imag, p.real)),
                        "%.6f" % float(spl(abs(p))), count])
