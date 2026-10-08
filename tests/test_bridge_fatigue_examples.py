# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
De voorbeelden van #ASTER-WOHLER in de documentatie (docs/Vermoeiing.md,
docs/en/Fatigue.md) en in de meldingen van Scripts/bridge_fatigue.py moeten
zelf geldig zijn: een gebruiker neemt ze letterlijk over. Gevonden bij de
docs-controle S11: het eerste voorbeeld (S aflopend) en de voorbeelden in
twee meldingen gaven de ERROR 'S moet oplopen'.

Gebruik: python tests/test_bridge_fatigue_examples.py (CPython 3 of 2.7).
"""
import io
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "Scripts"))

import bridge_utils  # noqa: E402
import bridge_fatigue  # noqa: E402

# Logregels niet naar een logbestand van de werkmap schrijven
_LOG = []
bridge_fatigue.log = _LOG.append

_LINE = re.compile(u"#ASTER-WOHLER ([^`'\n]+)")


def _examples(path):
    text = io.open(path, encoding="utf-8").read()
    out = []
    for m in _LINE.finditer(text):
        body = m.group(1).strip()
        if ".." in body or ":" not in body:
            continue        # schrijfwijze met [..] of een foutvoorbeeld
        out.append(body)
    return out


def _message_examples():
    """Voorbeelden uit de meldingen: de fout uitlokken en het voorbeeld uit
    de logregel halen."""
    out = []
    triggers = [u"Staal S = [200, 400] MPa, N = [1e6, 1e4]",
                u"Staal: S = [400, 200] MPa, N = [1e4, 1e6]",
                u"Staal: S = [200, 400], N = [1e6, 1e4]"]
    for body in triggers:
        del _LOG[:]
        assert bridge_fatigue._parse_curve(body, {"Staal": {}}) is None
        msg = u" ".join(_LOG)
        m = re.search(u"'#ASTER-WOHLER ([^']+)'", msg)
        if m:
            out.append(m.group(1))
            continue
        # alleen het voorbeeld achter "bijv."/"Voorbeeld:" (of "e.g."/
        # "Example:"), niet de regel van de gebruiker vooraan
        tail = re.split(r"bijv\.|Voorbeeld:|e\.g\.|Example:", msg)[-1]
        m = re.search(r"(S = \[[^\]]+\] MPa)(, N = \[[^\]]+\])?", tail)
        if not m:
            out.append(u"Staal: (geen voorbeeld in: {0})".format(msg))
            continue
        s = m.group(1)
        n = m.group(2)
        if not n:
            k = len(s.split(","))
            n = u", N = [{0}]".format(u", ".join(
                u"1e{0}".format(6 - i) for i in range(k)))
        out.append(u"Staal: " + s + n)
    return out


def main():
    files = [os.path.join(ROOT, "docs", "Vermoeiing.md"),
             os.path.join(ROOT, "docs", "en", "Fatigue.md")]
    failed = 0
    count = 0
    for lang in ("nl", "en"):
        bridge_utils.CONFIG_BRIDGE["language"] = lang
        cases = [(os.path.basename(p), b) for p in files
                 for b in _examples(p)]
        cases += [("bridge_fatigue.py", b) for b in _message_examples()]
        for where, body in cases:
            name = body.partition(":")[0].strip()
            materials = {name: {}} if name != "*" else {"Staal": {}}
            del _LOG[:]
            curve = bridge_fatigue._parse_curve(body, materials)
            count += 1
            if curve is None:
                failed += 1
                print(u"FOUT {0} ({1}): '{2}' -> {3}".format(
                    where, lang, body, u" | ".join(_LOG)))
    if count < 10:
        print("FOUT: maar {0} voorbeelden gevonden".format(count))
        return 1
    print("{0} voorbeelden, {1} fout".format(count, failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
