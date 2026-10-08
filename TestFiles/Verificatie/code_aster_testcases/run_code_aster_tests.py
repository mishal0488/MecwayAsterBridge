# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Officiele Code-Aster-validatietests draaien die de onderdelen raken die de
MecwayAsterBridge gebruikt (CPython 3, ontwikkelgereedschap; draait niet in
Mecway).

    python run_code_aster_tests.py               alle tests uit TESTS
    python run_code_aster_tests.py ssnp sdls115a  alleen namen die zo beginnen
    python run_code_aster_tests.py -j 4           4 tests tegelijk
    python run_code_aster_tests.py --keep         werkmap niet opruimen
    python run_code_aster_tests.py --report       alleen het rapport opnieuw

Elke test draait zoals ctest het doet (share/aster/run_aster_for_ctest.bat):
    run_aster.bat --ctest <tests>/<naam>.export  > <naam>.mess
in een eigen tijdelijke map (%TEMP%/ca_testcases/<naam>; niets in de repo).
Met --ctest gelden de bestanden van de .export relatief aan de testmap van de
installatie (run_aster/export.py, Export._abspath) en schrijft run_aster ook
<naam>.code (eenheid 15: alle gebruikte sleutelwoorden, ook standaardwaarden).

Uit de .mess: de diagnose (--- DIAGNOSTIC JOB : OK/NOOK_TEST_RESU/...) en elke
TEST_RESU/TEST_TABLE-regel " OK   <REFERENTIE> ..." of " NOOK <REFERENTIE> ..."
(geteld per referentietype: ANALYTIQUE, SOURCE_EXTERNE, AUTRE_ASTER,
NON_DEFINI, NON_REGRESSION). Uit de .code: welke MODELISATION, RELATION,
DEFORMATION, SCHEMA, FORMULATION en solver (SOLVEUR METHODE) echt gebruikt
zijn. Resultaat: resultaten.md en resultaten.csv naast dit script, en
resultaten.json (ruwe gegevens, per test bijgewerkt).

Code-Aster: CONFIG_BRIDGE["exe_aster"] (bridge_utils.exe_aster_path()) uit Scripts/bridge_utils.py, of de
omgevingsvariabele ASTER_RUN (pad naar run_aster.bat).
"""
import sys
import os
import io
import re
import csv
import json
import time
import shutil
import tempfile
import datetime
import subprocess
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
RAW = os.path.join(HERE, "resultaten.json")

# ------------------------------------------------------------------
# Gekozen tests. "bridge": welk onderdeel van de bridge dit dekt.
# V-documentatie: naam van de test = naam van de V-doc op code-aster.org
# (bijv. SSLL102 = V3.01.102; V2 = lineaire dynamica, V3 = lineaire
# statica, V4 = thermisch, V6 = niet-lineaire statica).
# ------------------------------------------------------------------
TESTS = [
    dict(name="sslv04b", bridge="MECA_STATIQUE 3D (hex20/penta15)",
         onderwerp="Dikwandige holle cilinder onder inwendige druk 60 MPa",
         wat="radiale verplaatsing en spanningen (Lame), 3D kwadratisch, PRES_REP",
         ref="analytisch (Lame)"),
    dict(name="ssll102a", bridge="MECA_STATIQUE POU_D_E",
         onderwerp="Ingeklemde balk met eenheidskrachten en -momenten op het vrije eind",
         wat="verplaatsingen, rotaties en snedekrachten van Euler-Bernoulli-balken",
         ref="analytisch (balktheorie)"),
    dict(name="ssls101b", bridge="MECA_STATIQUE DKT",
         onderwerp="Ronde plaat (kwart model) onder druk, F3 en eigen gewicht",
         wat="doorbuiging en momenten in DKT-schaalelementen; FORCE_COQUE, PESANTEUR",
         ref="analytisch (plaattheorie)"),
    dict(name="ssls105a", bridge="MECA_STATIQUE DKT (dubbel gekromd)",
         onderwerp="NAFEMS LE3: halve bol, dubbel geknepen (Morley)",
         wat="radiale verplaatsing in het belastingspunt; grote rotaties in het vlak",
         ref="NAFEMS LE3 (SOURCE_EXTERNE)"),
    dict(name="ssls114a", bridge="MECA_STATIQUE COQUE_3D",
         onderwerp="Cilinderschaal onder inwendige druk (PRES_REP en FORCE_COQUE)",
         wat="verplaatsingen, membraankrachten en rekken in COQUE_3D",
         ref="analytisch (membraantheorie)"),
    dict(name="ssnp15a", bridge="STAT_NON_LINE VMIS_ISOT_LINE",
         onderwerp="Vierkante plaat, niet-proportionele trek + afschuiving, 3D",
         wat="spanningen en plastische rek bij isotrope versteviging (Von Mises)",
         ref="analytisch (integratie van de constitutieve wet)"),
    dict(name="ssnp14b", bridge="STAT_NON_LINE VMIS_CINE_LINE",
         onderwerp="NAFEMS: vierkante plaat, trek + afschuiving, kinematische versteviging",
         wat="spanningen en plastische rek bij kinematische versteviging, penta6",
         ref="NAFEMS (in de .comm: NON_DEFINI)"),
    dict(name="hsnv100g", bridge="STAT_NON_LINE VMIS_ISOT_TRAC + thermische rek",
         onderwerp="Thermoplasticiteit: element met opgelegde temperatuur en trekkromme",
         wat="spanning en plastische rek; VMIS_ISOT_TRAC, ALPHA, AFFE_VARC TEMP",
         ref="analytisch"),
    dict(name="ssnp159a", bridge="STAT_NON_LINE GDEF_LOG",
         onderwerp="Blok in trek tot grote rek, VMIS_ISOT_LINE, PETIT / GDEF_LOG / SIMO_MIEHE",
         wat="spanning en rek bij grote vervorming (logaritmische rek)",
         ref="analytisch"),
    dict(name="ssnv122a", bridge="STAT_NON_LINE GREEN_LAGRANGE, meebewegende druk",
         onderwerp="Starre rotatie van een blok + druk in de vervormde configuratie",
         wat="objectiviteit van GREEN_LAGRANGE; PRES_REP met TYPE_CHARGE='SUIV'",
         ref="analytisch"),
    dict(name="ssnv167d", bridge="contact DEFI_CONTACT CONTINUE",
         onderwerp="Twee blokken in contact, opgelegde verplaatsing, 3D",
         wat="contactdruk (LAGS_C) en verplaatsingen; continue formulering, NEWTON",
         ref="analytisch"),
    dict(name="sdls115a", bridge="DYNA_NON_LINE DIFF_CENT en HHT",
         onderwerp="Plaatstrook (DKT) met sinusvormige oplossing u = q0 sin(kx) sin(wt)",
         wat="verplaatsing in de tijd met expliciet (DIFF_CENT) en HHT compleet",
         ref="analytisch"),
    dict(name="sdlv121d", bridge="DYNA_NON_LINE HHT / NEWMARK (continuum)",
         onderwerp="Vlakke golf in een elastische kolom: voortplanting, reflectie, absorptie",
         wat="verplaatsing/snelheid van de golf, D_PLAN, HHT en NEWMARK",
         ref="analytisch (golfvoortplanting)"),
    dict(name="tpll01e", bridge="THER_LINEAIRE stationair (3D)",
         onderwerp="Wand met opgelegde temperatuur, flux en convectie (1D-oplossing in 3D)",
         wat="temperatuur; TEMP_IMPO, FLUX_REP, ECHANGE",
         ref="analytisch"),
    dict(name="ttlv01a", bridge="THER_LINEAIRE transient (3D)",
         onderwerp="Afkoelend/opwarmend lichaam met convectie (penta6/tetra4)",
         wat="temperatuur in de tijd; ECHANGE, THETA-schema, ETAT_INIT",
         ref="analytisch (reeksoplossing)"),
    dict(name="tpna01b", bridge="THER_NON_LINE straling (3D)",
         onderwerp="Stationaire geleiding met straling naar de omgeving",
         wat="temperatuur en flux; RAYONNEMENT (T+273,15 intern), THER_NL",
         ref="analytisch"),
    dict(name="tplv101c", bridge="THER_NON_LINE PLAN (thermische schaal in de bridge)",
         onderwerp="Vlak model met bronnen, flux, opgelegde temperatuur en wandcontact",
         wat="temperatuur; MODELISATION='PLAN', SOURCE, FLUX_REP, TEMP_IMPO",
         ref="analytisch"),
    dict(name="sslv153a", bridge="CALC_ERREUR ERME (foutschatter)",
         onderwerp="Blok van hexa/penta/pyram/tetra (LINE_QUAD), druk, vrije randen met druk ~0",
         wat="ERREST van ERME_NOEU en QIRE_NOEU; ook via de bridge nagebouwd (FS-03)",
         ref="NON_REGRESSION"),
    dict(name="sslv111f", bridge="CALC_ERREUR ERME en ERZ1 (3D)",
         onderwerp="Plaat met gat (Kirsch), hexa8/pyram5/tetra4, FORCE_FACE",
         wat="ERME_ELEM/ERME_ELNO en ZZ1-spanningen in 3D",
         ref="analytisch (Kirsch) / NON_REGRESSION"),
]

REF_TYPES = ["ANALYTIQUE", "SOURCE_EXTERNE", "AUTRE_ASTER", "NON_DEFINI",
             "NON_REGRESSION"]
_LINE = re.compile(r"^ (OK|NOOK|SKIP)\s+([A-Z_]+)\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)"
                   r"\s+(\S+)")


def aster_exe():
    if os.environ.get("ASTER_RUN"):
        return os.environ["ASTER_RUN"]
    sys.path.insert(0, os.path.join(REPO, "Scripts"))
    import bridge_utils
    return bridge_utils.exe_aster_path()


def tests_dir(exe):
    # <root>/bin/run_aster.bat -> <root>/share/aster/tests
    return os.path.join(os.path.dirname(os.path.dirname(exe)), "share",
                        "aster", "tests")


def _export_info(path):
    txt = io.open(path, encoding="utf-8", errors="replace").read()
    tl = re.search(r"P time_limit (\S+)", txt)
    return float(tl.group(1)) if tl else 60.0


def parse_mess(text):
    """Diagnose en TEST_RESU-regels uit een .mess."""
    m = re.search(r"--- DIAGNOSTIC JOB : (\S+)", text)
    out = {"diagnostic": m.group(1) if m else "?", "lines": []}
    for line in text.splitlines():
        mm = _LINE.match(line)
        if mm:
            out["lines"].append({"status": mm.group(1), "ref": mm.group(2),
                                 "refe": mm.group(3), "calc": mm.group(4),
                                 "err": mm.group(5), "tole": mm.group(6)})
    out["errors"] = sorted(set(re.findall(
        r"<(?:F|EXCEPTION|S)>\s*<([A-Za-z0-9_]+)>", text)))
    out["alarms"] = sorted(set(re.findall(r"<A>\s*<([A-Za-z0-9_]+)>", text)))
    return out


def parse_code(text):
    """Gebruikte sleutelwoordwaarden uit de .code (TEST <CMD> <FACT> <KW> 'v')."""
    keys = ("MODELISATION", "RELATION", "DEFORMATION", "SCHEMA",
            "FORMULATION", "METHODE")
    found = dict((k, set()) for k in keys)
    for line in text.splitlines():
        p = line.split()
        if len(p) >= 5 and p[0] == "TEST" and p[3] in keys:
            v = p[4].strip("'")
            if p[3] == "METHODE" and p[2] != "SOLVEUR":
                continue
            if p[3] == "DEFORMATION" and p[2] != "COMPORTEMENT":
                continue        # CALC_CHAMP DEFORMATION=('EPSI_ELNO', ..)
            if re.match(r"^[A-Z0-9][A-Z0-9_]*$", v):
                found[p[3]].add(v)
    return dict((k, sorted(v)) for k, v in found.items())


def run_one(test, exe, root):
    name = test["name"]
    export = os.path.join(tests_dir(exe), name + ".export")
    work = os.path.join(root, name)
    if os.path.isdir(work):
        shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    row = dict(test)
    row.update({"diagnostic": "", "seconds": None, "counts": {},
                "max_err": {}, "features": {}, "errors": [], "alarms": []})
    if not os.path.isfile(export):
        row["diagnostic"] = "GEEN_EXPORT"
        return row
    tl = _export_info(export)
    mess = os.path.join(work, name + ".mess")
    t0 = time.time()
    with open(mess, "wb") as f:
        try:
            subprocess.call([exe, "--ctest", export], stdout=f,
                            stderr=subprocess.STDOUT, cwd=work,
                            timeout=max(600.0, 10.0 * tl))
        except subprocess.TimeoutExpired:
            row["errors"].append("TIMEOUT")
    row["seconds"] = round(time.time() - t0, 1)
    text = io.open(mess, encoding="utf-8", errors="replace").read()
    res = parse_mess(text)
    row["diagnostic"] = res["diagnostic"]
    row["errors"] += res["errors"]
    row["alarms"] = res["alarms"]
    counts = {}
    max_err = {}
    for ln in res["lines"]:
        c = counts.setdefault(ln["ref"], {"OK": 0, "NOOK": 0, "SKIP": 0})
        c[ln["status"]] += 1
        # ERREUR en TOLE: relatief in % (CRITERE='RELATIF') of absoluut;
        # bewaard: grootste relatieve fout met zijn tolerantie.
        e, tol = ln["err"], ln["tole"]
        if ln["ref"] != "NON_REGRESSION" and e.endswith("%"):
            try:
                val = abs(float(e.rstrip("%")))
                if val >= max_err.get(ln["ref"], (-1.0, 0.0))[0]:
                    max_err[ln["ref"]] = (val, float(tol.rstrip("%")))
            except ValueError:
                pass
    row["counts"] = counts
    row["max_err"] = max_err
    code = os.path.join(work, name + ".code")
    if os.path.isfile(code):
        row["features"] = parse_code(io.open(code, encoding="utf-8",
                                             errors="replace").read())
    return row


def _verdict(r):
    # <A>_ALARM: geslaagd met alarmen (run_aster geeft dan returncode 0,
    # ctest telt de test als geslaagd; run_aster/status.py).
    nook = sum(c.get("NOOK", 0) for c in r["counts"].values())
    ok = sum(c.get("OK", 0) for c in r["counts"].values())
    if r["diagnostic"] in ("OK", "<A>_ALARM") and nook == 0 and ok > 0:
        return "OK"
    if nook:
        return "NOOK"
    return r["diagnostic"] or "FOUT"


def _counts_txt(r, sep=", "):
    parts = []
    for ref in REF_TYPES:
        c = r["counts"].get(ref)
        if c:
            s = "%s %d" % (ref, c["OK"] + c["NOOK"] + c["SKIP"])
            if c["NOOK"]:
                s += " (%d NOOK)" % c["NOOK"]
            parts.append(s)
    return sep.join(parts)


def _err_txt(r):
    return ", ".join("%s %.2g%% (tol. %.2g%%)" % (k, v[0], v[1])
                     for k, v in sorted(r["max_err"].items()))


def _features_txt(r):
    f = r.get("features") or {}
    parts = []
    for k, lab in (("MODELISATION", "mod"), ("RELATION", "rel"),
                   ("DEFORMATION", "def"), ("SCHEMA", "schema"),
                   ("FORMULATION", "form"), ("METHODE", "solver")):
        if f.get(k):
            parts.append("%s %s" % (lab, "/".join(f[k])))
    return "; ".join(parts)


def _anon_path(path):
    """Pad zonder gebruikers- of machinegegevens voor het rapport: het
    begin %LOCALAPPDATA% / %USERPROFILE% of <REPO> in plaats van de map."""
    norm = os.path.normcase(os.path.normpath(path))
    for var in ("LOCALAPPDATA", "USERPROFILE"):
        root = os.environ.get(var)
        if root and norm.startswith(os.path.normcase(
                os.path.normpath(root)) + os.sep):
            return "%" + var + "%" + os.path.normpath(path)[len(
                os.path.normpath(root)):]
    if norm.startswith(os.path.normcase(REPO) + os.sep):
        return "<REPO>" + os.path.normpath(path)[len(REPO):]
    return path


def write_report(rows, exe):
    order = [t["name"] for t in TESTS]
    rows = [rows[n] for n in order if n in rows]
    n_ok = sum(1 for r in rows if _verdict(r) == "OK")
    n_lines = sum(sum(c["OK"] + c["NOOK"] for c in r["counts"].values())
                  for r in rows)
    n_nook = sum(sum(c["NOOK"] for c in r["counts"].values()) for r in rows)
    out = [u"# Code-Aster-validatietests (installatie)", u"",
           u"Gemaakt door `run_code_aster_tests.py` op %s met `%s`. Elke test "
           u"draait als `run_aster.bat --ctest <test>.export` (zoals ctest) in "
           u"een tijdelijke map. Toelichting: `docs/Verificatie.md` (hoofdstuk Code-Aster-tests)." % (
               datetime.date.today().isoformat(), _anon_path(exe)), u"",
           u"**%d van %d tests OK**; %d TEST_RESU-controles, waarvan %d NOOK."
           % (n_ok, len(rows), n_lines, n_nook), u"",
           u"| Test | Bridge-onderdeel | Onderwerp | Wat geverifieerd wordt | "
           u"Referentie | Controles (type: aantal) | Max. fout t.o.v. referentie | "
           u"Uitkomst | Tijd [s] |",
           u"|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        out.append(u"| %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
            r["name"], r["bridge"], r["onderwerp"], r["wat"], r["ref"],
            _counts_txt(r), _err_txt(r), _verdict(r),
            "" if r["seconds"] is None else "%.0f" % r["seconds"]))
    out += [u"", u"## Gebruikte Code-Aster-opties (uit de `.code`)", u"",
            u"| Test | Opties |", u"|---|---|"]
    for r in rows:
        out.append(u"| %s | %s |" % (r["name"], _features_txt(r)))
    al = [r for r in rows if r.get("alarms")]
    if al:
        out += [u"", u"## Alarmen (diagnose <A>_ALARM, test geslaagd)", u"",
                u"| Test | Alarmen |", u"|---|---|"]
        for r in al:
            out.append(u"| %s | %s |" % (r["name"], ", ".join(r["alarms"])))
    bad = [r for r in rows if _verdict(r) != "OK" or r["errors"]]
    if bad:
        out += [u"", u"## Afwijkingen", u""]
        for r in bad:
            out.append(u"- %s: diagnose %s, %s; meldingen %s" % (
                r["name"], r["diagnostic"], _counts_txt(r),
                ", ".join(r["errors"]) or "-"))
    io.open(os.path.join(HERE, "resultaten.md"), "w", encoding="utf-8",
            newline="\n").write(u"\n".join(out) + u"\n")
    with io.open(os.path.join(HERE, "resultaten.csv"), "w", encoding="utf-8",
                 newline="") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["test", "bridge_onderdeel", "onderwerp", "wat",
                    "referentie", "diagnose", "uitkomst"] +
                   ["n_" + t for t in REF_TYPES] + ["n_nook", "max_fout",
                                                     "opties", "alarmen", "seconden"])
        for r in rows:
            w.writerow([r["name"], r["bridge"], r["onderwerp"], r["wat"],
                        r["ref"], r["diagnostic"], _verdict(r)] +
                       [sum(r["counts"].get(t, {}).values()) for t in REF_TYPES] +
                       [sum(c["NOOK"] for c in r["counts"].values()),
                        _err_txt(r), _features_txt(r),
                        " ".join(r.get("alarms", [])), r["seconds"]])


def main(argv):
    exe = aster_exe()
    rows = {}
    if os.path.isfile(RAW):
        rows = json.load(io.open(RAW, encoding="utf-8"))
    if "--report" in argv:
        write_report(rows, exe)
        return
    jobs = 1
    if "-j" in argv:
        jobs = int(argv[argv.index("-j") + 1])
        argv = argv[:argv.index("-j")] + argv[argv.index("-j") + 2:]
    sel = [a for a in argv if not a.startswith("--")]
    tests = [t for t in TESTS if not sel or
             any(t["name"].startswith(s) for s in sel)]
    root = os.environ.get("CA_TESTS_WORK") or os.path.join(
        tempfile.gettempdir(), "ca_testcases")
    if not os.path.isdir(root):
        os.makedirs(root)

    def job(t):
        r = run_one(t, exe, root)
        print("%-10s %-8s %6ss  %s" % (t["name"], _verdict(r), r["seconds"],
                                       _counts_txt(r)))
        sys.stdout.flush()
        return r

    with ThreadPoolExecutor(max_workers=jobs) as ex:
        for r in ex.map(job, tests):
            rows[r["name"]] = r
    json.dump(rows, io.open(RAW, "w", encoding="utf-8"), indent=1)
    write_report(rows, exe)
    if "--keep" not in argv:
        shutil.rmtree(root, ignore_errors=True)
    else:
        print("werkmap:", root)


if __name__ == "__main__":
    main(sys.argv[1:])
