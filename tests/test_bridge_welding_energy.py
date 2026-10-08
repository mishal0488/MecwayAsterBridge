# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Energiecontrole van de lasbron (bridge_welding.report_lines): neemt het
mesh de helft of een kwart van het vermogen op, dan noemt de WARNING
symmetrie= (docs/Lassen.md par. 8: de baan op een symmetrievlak zonder
symmetrie=). Gevonden bij de docs-controle S11 (WE-02 zonder symmetrie=1:
50,0%, de WARNING noemde alleen een grof mesh).

Gebruik: python tests/test_bridge_welding_energy.py (CPython 3 of 2.7).
"""
import io
import os
import shutil
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "Scripts"))

import bridge_utils  # noqa: E402
import bridge_welding  # noqa: E402

_LOG = []
bridge_welding.log = _LOG.append


def _run(fraction, lang):
    """WARNING-regels voor een .mess waarin het mesh 'fraction' van 1000 W
    opneemt."""
    bridge_utils.CONFIG_BRIDGE["language"] = lang
    d = tempfile.mkdtemp()
    try:
        mess = os.path.join(d, "x.mess")
        f = io.open(mess, "w", encoding="utf-8")
        for t in (0.6, 1.8, 3.0):
            f.write(u"BRIDGE_LASBRON 1 {0} {1} 1000.0\n".format(
                t, 1000.0 * fraction))
        f.close()
        del _LOG[:]
        bridge_welding.report_lines(mess)
    finally:
        shutil.rmtree(d)
    return [l for l in _LOG if l.startswith("WARNING")]


def main():
    fails = []
    for lang, key in (("nl", u"symmetrie=1"), ("en", u"symmetry=1")):
        w = _run(0.5, lang)
        if len(w) != 1 or key not in w[0]:
            fails.append((lang, 0.5, w))
        w = _run(0.25, lang)
        if len(w) != 1 or key.replace("1", "2") not in w[0]:
            fails.append((lang, 0.25, w))
        w = _run(0.8, lang)
        if len(w) != 1 or u"symmetr" in w[0]:
            fails.append((lang, 0.8, w))
        w = _run(1.005, lang)
        if w:
            fails.append((lang, 1.005, w))
    for f in fails:
        print("FOUT: {0}".format(f))
    if not fails:
        print("OK: energiecontrole lasbron (symmetrie-hint nl/en)")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
