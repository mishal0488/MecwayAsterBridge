# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Offline test van de #ASTER-invoer (Python 2.7 en 3, zonder Mecway en
# Code-Aster):
#   - #ASTER-BRON-waarden (bridge_acoustic.source_value,
#     bridge_vibroacoustic_harm._aster_sources, bridge_acoustic.
#     _parse_aster_sources): getal + bekende eenheid, druk een keer naar Pa,
#     anders een fout;
#   - ja/nee-woorden van #ASTER-OPTIE (bridge_utils.parse_yes_no,
#     bridge_custom.apply_options, bridge_adapt): bekende woorden in NL/EN/FR,
#     anders een WARNING en de standaardwaarde.
#
#   python tests/test_bridge_invoer.py

from __future__ import print_function

import os
import sys

_HERE = os.path.dirname(os.path.abspath(globals().get("__file__",
                                                     sys.argv[0])))
sys.path.insert(0, os.path.join(_HERE, "..", "Scripts"))
import bridge_utils  # noqa: E402
bridge_utils.CONFIG_BRIDGE["language"] = "nl"

# Logregels opvangen (alleen in het geheugen; geen bridge_log.txt).
LOG = []


def _log(msg):
    LOG.append(msg)


bridge_utils.log = _log
import bridge_custom  # noqa: E402
import bridge_adapt  # noqa: E402
import bridge_acoustic as bac  # noqa: E402
import bridge_vibroacoustic_harm as bvh  # noqa: E402
for _mod in (bridge_custom, bridge_adapt, bac, bvh):
    _mod.log = _log

FAILS = []


def check(cond, text):
    if not cond:
        FAILS.append(text)


def close(a, b):
    return a is not None and abs(a - b) <= 1e-9 * max(1.0, abs(b))


def test_source_value():
    psi = 4.4482216152605 / 0.0254 ** 2
    good = [("2", "velocity", 2.0), ("2 m/s", "velocity", 2.0),
            ("2 mm/s", "velocity", 2.0e-3), ("2mm/s", "velocity", 2.0e-3),
            ("1 in/s", "velocity", 0.0254), ("1 ft/s", "velocity", 0.3048),
            ("3 Pa", "pressure", 3.0), ("2 N/m2", "pressure", 2.0),
            (u"2 N/m²", "pressure", 2.0),
            ("2 kPa", "pressure", 2.0e3), ("2 kN/m2", "pressure", 2.0e3),
            ("2 MPa", "pressure", 2.0e6), ("2 mPa", "pressure", 2.0e-3),
            ("2 N/mm2", "pressure", 2.0e6), ("2 GPa", "pressure", 2.0e9),
            ("1 bar", "pressure", 1.0e5), ("1 mbar", "pressure", 1.0e2),
            ("2 psi", "pressure", 2.0 * psi), ("-1.5e-3 Pa", "pressure",
                                                -1.5e-3)]
    for text, kind, value in good:
        k, v, _u = bac.source_value(text)
        check(k == kind and close(v, value),
              "source_value({0!r}) = ({1}, {2}), verwacht ({3}, {4})".format(
                  text, k, v, kind, value))
    for text in ("2 km/h", "2 m", "2 N", "abc", "", "2 bar x", "2 Pa/s"):
        k, v, _u = bac.source_value(text)
        check(k is None and v is None,
              "source_value({0!r}) = ({1}, {2}), verwacht een fout".format(
                  text, k, v))


def test_vibro_sources():
    md = {"aster_custom": {"sources": [
        "A = 1 bar", "B = 2 mm/s, 90", "C = 2 N/m2, fase = 45",
        "D = 2 GPa", "E = 2 psi", "F = 1 mbar",
        "X = 2 km/h", "Y = 2 m", "Z = 2 Pa, x", "W = geen getal"]}}
    del LOG[:]
    out, bad = bvh._aster_sources(md)
    got = dict((s, (k, v, p)) for s, k, v, p in out)
    psi = 4.4482216152605 / 0.0254 ** 2
    expect = {"A": ("pressure", 1.0e5, 0.0), "B": ("velocity", 2.0e-3, 90.0),
              "C": ("pressure", 2.0, 45.0), "D": ("pressure", 2.0e9, 0.0),
              "E": ("pressure", 2.0 * psi, 0.0),
              "F": ("pressure", 100.0, 0.0)}
    for sel, (kind, value, phase) in expect.items():
        k, v, p = got.get(sel, (None, None, None))
        check(k == kind and close(v, value) and p == phase,
              "#ASTER-BRON {0}: ({1}, {2}, {3}), verwacht ({4}, {5}, "
              "{6})".format(sel, k, v, p, kind, value, phase))
    check(sorted(got) == sorted(expect),
          "#ASTER-BRON: onverwachte bronnen {0}".format(sorted(got)))
    check(len(bad) == 4, "#ASTER-BRON: {0} fouten, verwacht 4: {1}".format(
        len(bad), bad))
    check(any("X = 2 km/h" in b and "'km/h' niet herkend" in b for b in bad),
          "#ASTER-BRON X = 2 km/h: geen melding van de eenheid: {0}".format(
              bad))
    check(not any("niet herkend" in b for b in bad if "Z = " in b),
          "#ASTER-BRON Z (foute fase): onterecht 'eenheid niet herkend'")
    check(not any("Eenheid niet herkend" in m for m in LOG),
          "#ASTER-BRON: onterechte WARNING 'Eenheid niet herkend': "
          "{0}".format(LOG))


def test_acoustic_sources():
    md = {"aster_custom": {"sources": [
        "A = 2 mm/s, 90", "B = 1 bar", "C = 2 km/h", "D = 3"]}}
    ac = {"aster_sources": [], "bad": []}
    bac._parse_aster_sources(md, ac)
    check([(s, round(v, 12), p) for s, v, p in ac["aster_sources"]] ==
          [("A", 0.002, 90.0), ("D", 3.0, 0.0)],
          "akoestiek #ASTER-BRON: {0}".format(ac["aster_sources"]))
    check(len(ac["bad"]) == 2 and "druk" in ac["bad"][0],
          "akoestiek #ASTER-BRON fouten: {0}".format(ac["bad"]))


def test_yes_no():
    for word in ("1", "true", "True", "yes", "YES", "on", "ja", "Aan",
                 "waar", "oui", "vrai", " ja ", "'ja'", '"on"'):
        check(bridge_utils.parse_yes_no(word) is True,
              "parse_yes_no({0!r}) is niet True".format(word))
    for word in ("0", "false", "no", "off", "Nee", "uit", "onwaar", "non",
                 "faux", "OFF"):
        check(bridge_utils.parse_yes_no(word) is False,
              "parse_yes_no({0!r}) is niet False".format(word))
    for word in (None, "", "j", "jaa", "2", "misschien", "aan uit"):
        check(bridge_utils.parse_yes_no(word) is None,
              "parse_yes_no({0!r}) is niet None".format(word))


def test_apply_options():
    analysis = {"a": False, "b": True, "c": False, "d": True, "n": 3}
    del LOG[:]
    res = bridge_custom.apply_options(
        analysis, {"a": "aan", "b": "uit", "c": "jaa", "d": "misschien",
                   "n": "5"})
    check(res["a"] is True and res["b"] is False,
          "apply_options aan/uit: {0}".format(res))
    check(res["c"] is False and res["d"] is True,
          "apply_options ongeldig: standaardwaarde niet behouden: "
          "{0}".format(res))
    check(res["n"] == 5, "apply_options int: {0}".format(res))
    warns = [m for m in LOG if m.startswith("WARNING:")]
    check(len(warns) == 2 and any("c = jaa" in w for w in warns) and
          any("d = misschien" in w for w in warns),
          "apply_options: verwacht 2 WARNINGs (c, d): {0}".format(LOG))


def test_adapt():
    analysis = dict(bridge_adapt.DEFAULTS)
    analysis.update({"solver_command": "MECA_STATIQUE",
                     "error_estimator": False})
    del LOG[:]
    s = bridge_adapt.adapt_settings(
        analysis, {"adapt_refine": "on", "error_estimator": "aan",
                   "adapt_neighbours": "misschien"})
    check(s is not None and s["adapt_refine"] is True,
          "adapt_settings adapt_refine = on: {0}".format(s))
    check(s is not None and s["adapt_neighbours"] is False,
          "adapt_settings ongeldige adapt_neighbours: {0}".format(s))
    s = bridge_adapt.adapt_settings(analysis, {"adapt_refine": "jaa"})
    check(s is None, "adapt_settings adapt_refine = jaa: niet uit: "
          "{0}".format(s))


def main():
    for test in (test_source_value, test_vibro_sources, test_acoustic_sources,
                 test_yes_no, test_apply_options, test_adapt):
        test()
    if FAILS:
        for f in FAILS:
            print(u"FOUT: {0}".format(f).encode("ascii", "replace").decode(
                "ascii"))
        return 1
    print("OK: #ASTER-invoer (#ASTER-BRON-eenheden akoestiek en "
          "vibro-akoestiek, ja/nee-woorden #ASTER-OPTIE en adaptief)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
