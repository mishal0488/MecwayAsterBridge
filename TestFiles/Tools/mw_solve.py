# Mecway headless laten rekenen (referentie voor verificatie, geen GUI nodig).
# gebruik (CPython): python mw_solve.py <in.liml> <uit_basisnaam>  -> <uit>_sol.liml
# Mecway rekent met de solver uit <analysis> (solver="ccx" voor CCX). Het
# resultaat lezen: aster_results_export-formaat, zie TestFiles/Verificatie/verif_lib.py.
# Omgeving: MWT = tijdslimiet in s (standaard 300).
import sys, os, io, re, subprocess, time
src, base = sys.argv[1], os.path.abspath(sys.argv[2])
txt = io.open(src, encoding="utf-8").read()
txt = re.sub(r"\s*<solution>.*?</solution>", "", txt, flags=re.S)
txt = re.sub(r"\s*<solutionfile[^>]*/>", "", txt, flags=re.S)
sol = base + "_sol.liml"
if os.path.exists(sol): os.remove(sol)
txt = txt.replace("</liml>", '  <solutionfile filename="{0}" />\n</liml>'.format(sol.replace("/", "\\")))
inp = base + "_in.liml"
io.open(inp, "w", encoding="utf-8").write(txt)
t0 = time.time()
try:
    r = subprocess.run([r"C:\Program Files\Mecway\mecway.exe", inp.replace("/", "\\"), "solve"], timeout=float(os.environ.get("MWT", "300")))
    print("exit", r.returncode, "%.1f s" % (time.time() - t0))
except subprocess.TimeoutExpired:
    print("TIMEOUT")
print("solution:", os.path.exists(sol), os.path.getsize(sol) if os.path.exists(sol) else 0)
