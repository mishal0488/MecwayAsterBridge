# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Modaal met cyclische symmetrie (Mecway V30 met <cyclicsymmetry> ->
MODE_ITER_CYCL, Scripts/bridge_cyclic.py; bridge-analyse modal_3d).
Toelichting: docs/Verificatie.md par. 7 (CY) en docs/Modaal_en_knik.md
par. 7.

Mode-volgorde (zoals Mecway): per aantal nodal diameters ND = 0 .. N/2
het aantal modes maxmode, oplopend; mode k van ND d is dus mode
d * maxmode + k.

Referenties:
  - dunne ring in het vlak, vrij: f_n = 1/(2 pi) sqrt(EI/(rho A R^4))
    n (n^2 - 1)/sqrt(n^2 + 1), I = h t^3/12 voor buiging in het vlak
    (Blevins, Formulas for Natural Frequency and Mode Shape, 1979, tabel
    9-2); n = ND;
  - ingeklemde cirkelvormige plaat: f = lambda^2/(2 pi a^2) sqrt(D/(rho
    h)), D = E h^3/(12 (1 - nu^2)), lambda^2 van (n, s) = (0, 1) 10,2158,
    (1, 1) 21,260, (2, 1) 34,877, (0, 2) 39,771, (3, 1) 51,030, (4, 1)
    69,666 (Leissa, Vibration of Plates, NASA SP-160, 1969, tabel 2.1);
    n = ND;
  - het volledige model (de sector N keer, dezelfde elementen) met de
    gewone modale route van de bridge (CY-02, CY-06): elke cyclische
    frequentie moet in het volledige spectrum zitten;
  - Mecway 33 headless (TestFiles/Tools/mw_solve.py, interne solver,
    Mecway's eigen cyclische symmetrie) op dezelfde .liml.
"""
import io
import json
import math
import os
import tempfile
from verif_lib import Case, Check, Model, Results, hexmesh

GROUP = "Cyclische symmetrie"
E, NU, RHO = 210e9, 0.3, 7800.0
STEEL = ('<mechanical type="Isotropic" youngsmodulus="210000 MPa" '
         'poissonratio="0.3" />')
DENS = '<density type="Simple" density="7800 kg/m³" />'

# Ring: middenstraal R, radiale dikte T, hoogte H (as-richting).
R, T, H = 0.5, 0.01, 0.03
# Plaat: straal A, dikte HP; schijf: naaf R0 .. A.
A, HP, R0 = 0.5, 0.01, 0.1
LEISSA = {(0, 1): 10.2158, (1, 1): 21.260, (2, 1): 34.877, (0, 2): 39.771,
          (3, 1): 51.030, (4, 1): 69.666}


def _axis_map(axis):
    """Afbeelding van (x, y, z) met de symmetrieas langs z naar het model
    met de as langs axis (cyclische verwisseling, zelfde handigheid)."""
    if axis == "z":
        return lambda p: p
    if axis == "y":
        return lambda p: (p[1], p[2], p[0])
    return lambda p: (p[2], p[0], p[1])


def _cyclic(m, nseg, axis, nodes0, nodes1):
    m.node_selection("Cyclic_master", sorted(nodes0))
    m.node_selection("Cyclic_slave", sorted(nodes1))
    m.add('<cyclicsymmetry axis="%d" segments="%d" selectionmaster='
          '"Cyclic_master" selectionslave="Cyclic_slave" />'
          % ({"x": 1, "y": 2, "z": 3}[axis], nseg))


# ==========================================
# RING (vrij, buiging in het vlak)
# ==========================================
def ring(nseg=12, nel=6, sectors=1, axis="z", maxmode=3):
    """Ring van hex20 (1 radiaal, 1 axiaal, nel per segment); sectors =
    nseg: het volledige model, anders een sector met cyclische
    symmetrie."""
    def build():
        m = Model("V30", attrs={"maxmode": maxmode})
        m.elset("Ring", "Staal")
        m.mat("Staal", STEEL, DENS)
        ang = 2.0 * math.pi * sectors / nseg
        f = _axis_map(axis)

        def mp(u, v, w):
            r = R - T / 2.0 + u * T
            th = v * ang
            return f((r * math.cos(th), r * math.sin(th), -H / 2.0 + w * H))
        g = hexmesh(m, "Ring", 1, nel * sectors, 1, mp, order=2)
        if sectors < nseg:
            _cyclic(m, nseg, axis, g["nodes"]["v0"], g["nodes"]["v1"])
        return m
    return build


def blevins_ring(n):
    i = H * T ** 3 / 12.0
    a = T * H
    return (math.sqrt(E * i / (RHO * a * R ** 4)) * n * (n * n - 1) /
            math.sqrt(n * n + 1.0) / (2.0 * math.pi))


# ==========================================
# PLAAT EN SCHIJF (hex20, in het midden wedge15)
# ==========================================
def _wedges(m, elset, nt, ang, r1, f, nz=1, order=2):
    """Middelste ring van wedge15 (Code-Aster PENTA15-volgorde, zoals
    bridge_mesh._ASTER_TYPES): onder c, p1, p2; boven c', p1', p2'."""
    for j in range(nt):
        t0, t1 = ang * j / nt, ang * (j + 1) / nt
        tm = 0.5 * (t0 + t1)
        for k in range(nz):
            z0, z1 = HP * k / nz, HP * (k + 1) / nz
            zm = 0.5 * (z0 + z1)

            def P(r, t, z):
                return m.node(*f((r * math.cos(t), r * math.sin(t), z)))
            c = [P(0, 0, z0), P(r1, t0, z0), P(r1, t1, z0),
                 P(0, 0, z1), P(r1, t0, z1), P(r1, t1, z1)]
            if order == 1:
                m.elem(elset, "wedge6", c)
                continue
            mids = [P(r1 / 2, t0, z0), P(r1, tm, z0), P(r1 / 2, t1, z0),
                    P(0, 0, zm), P(r1, t0, zm), P(r1, t1, zm),
                    P(r1 / 2, t0, z1), P(r1, tm, z1), P(r1 / 2, t1, z1)]
            m.elem(elset, "wedge15", c + mids)


def disc(nseg=8, nr=8, nt=2, sectors=1, axis="z", maxmode=2, hub=0.0,
         rpm=0.0, clamp="outer", order=2):
    """
    Cirkelvormige plaat (hub = 0: midden van wedge15) of schijf (hub > 0)
    met straal A en dikte HP; nr elementen radiaal, nt per segment. clamp:
    "outer" (rand ingeklemd) of "inner" (naaf ingeklemd). rpm > 0:
    centrifugaalkracht om de symmetrieas (voorspanning).
    """
    def build():
        m = Model("V30", attrs={"maxmode": maxmode})
        m.elset("Plaat", "Staal")
        m.mat("Staal", STEEL, DENS)
        ang = 2.0 * math.pi * sectors / nseg
        f = _axis_map(axis)
        n_t = nt * sectors
        if hub > 0.0:
            r_in, n_hex = hub, nr
        else:
            r_in, n_hex = A / nr, nr - 1
            _wedges(m, "Plaat", n_t, ang, r_in, f, order=order)

        def mp(u, v, w):
            r = r_in + u * (A - r_in)
            th = v * ang
            return f((r * math.cos(th), r * math.sin(th), w * HP))
        g = hexmesh(m, "Plaat", n_hex, n_t, 1, mp, order=order)
        if clamp:
            side = "u1" if clamp == "outer" else "u0"
            m.face_selection("Inklemming", g["faces"][side])
            m.add('<fixed selection="Inklemming" />')
        if rpm:
            vec = f((0.0, 0.0, 1.0))
            m.add('<centrifugal><x>%r r/min</x><y>%r r/min</y><z>%r r/min'
                  '</z></centrifugal>' % tuple(rpm * c for c in vec))
        if sectors < nseg:
            # snijvlakken op coordinaat: hoek 0 en de sectorhoek
            back = {"z": lambda p: p, "y": lambda p: (p[2], p[0], p[1]),
                    "x": lambda p: (p[1], p[2], p[0])}[axis]
            n0, n1 = [], []
            for nid, p in enumerate(m.nodes, 1):
                x, y, _z = back(p)
                r = math.hypot(x, y)
                if r < 1e-9 or abs(y) < 1e-9 * A and x > 0:
                    n0.append(nid)
                if r < 1e-9 or abs(math.atan2(y, x) - ang) < 1e-9:
                    n1.append(nid)
            _cyclic(m, nseg, axis, n0, n1)
        return m
    return build


def leissa(n, s):
    d = E * HP ** 3 / (12.0 * (1.0 - NU ** 2))
    return LEISSA[(n, s)] / (2.0 * math.pi * A ** 2) * math.sqrt(
        d / (RHO * HP))


# ==========================================
# CONTROLES
# ==========================================
def _mode(maxmode, nd, k):
    """Modenummer (vanaf 1) van de k-de mode (vanaf 1) van ND nd."""
    return nd * maxmode + k


def _check(name, mode, ref, tol, source, tol_abs=None):
    return Check(name, lambda r, mode=mode: r.mode(mode), ref=ref,
                 tol=tol, tol_abs=tol_abs, unit="Hz", source=source)


def _nearest(target):
    """Frequentie van het volledige model die het dichtst bij target ligt."""
    def value(r):
        return min(r.modes.values(), key=lambda f: abs(f - target))
    return value


def _mecway_checks(values, tol=0.01, zero=0.05):
    out = []
    for k, f in enumerate(values):
        if abs(f) < zero:
            out.append(Check("Mecway f%d" % (k + 1),
                             lambda r, k=k: r.mode(k + 1), ref=0.0,
                             tol_abs=zero, unit="Hz",
                             source="Mecway 33 headless (starre mode)"))
        else:
            out.append(_check("Mecway f%d" % (k + 1), k + 1, f, tol,
                              "Mecway 33 headless (mw_solve.py), zelfde "
                              ".liml"))
    return out


# Mecway 33 headless (mw_solve.py) op modellen/CY-*.liml, 7 okt 2026.
MECWAY = {
    # ringsector om Z, 12 segmenten, maxmode 2: per ND 0 .. 6 twee modes
    "CY-01": (0.00096290910183046, 0.00309681618712625, 0.0, 0.0,
              25.7384754968607, 52.8790985334051, 72.8659040204307,
              173.475377794062, 139.891464017042, 358.343483306869,
              226.604217983236, 459.542969643256, 333.089053558176,
              333.089053571122),
    # dezelfde ring om X
    "CY-03": (0.0, 0.00397827377439815, 0.0, 0.00278032243522335,
              25.738475655284, 52.8790985087909, 72.8659039074771,
              173.475377777203, 139.89146393205, 358.343483302234,
              226.604217937623, 459.542969654955, 333.089053545519,
              333.0890535642),
    # ingeklemde plaat, 6 segmenten, knopen op de as
    "CY-04": (102.768587808574, 400.008035037742, 213.801648843324,
              611.468050286028, 350.652482022328, 702.582953314467,
              513.433362037578, 513.43336203981),
    # draaiende schijf om Y (3000 r/min), naaf ingeklemd, 12 segmenten
    "CY-05": (71.9145913242386, 365.191485188694, 75.4345914705494,
              387.439411607865, 100.636767436413, 460.551903038563,
              162.594098760447, 594.379095799437, 258.626226704608,
              788.964251329896, 387.416328902474, 761.075646768225,
              552.639912163602, 552.639912181229),
    # volledige draaiende schijf (CY-06), de eerste 30 modes
    "CY-06": (71.9145912364453, 75.4345914895843, 75.4345914961295,
              100.636767416218, 100.636767435022, 162.594098743408,
              162.594098748921, 258.626226689256, 258.62622669762,
              365.191485193235, 387.41632891446, 387.416328916547,
              387.439411607947, 387.439411608633, 460.551903051982,
              460.551903058373, 552.639912164824, 552.639912175433,
              594.379095782603, 594.379095784395, 610.630843864692,
              761.075646765921, 761.075646767283, 788.964251330983,
              788.964251334768, 1020.98545202429, 1020.98545204173,
              1039.3696363479, 1039.3696363488, 1056.72846482094),
}

# Het volledige model met de gewone modale route van de bridge (Code-Aster
# CALC_MODES, CY-02 en CY-06, 7 okt 2026): per ND de frequenties die bij de
# cyclische modes horen (paren: beide gelijk).
FULL_RING = {2: (25.7384758244533, 52.879098477051),
             3: (72.8659040148293, 173.475377776945),
             4: (139.891464012766, 358.343483310295),
             5: (226.604217967619, 459.542969649303),
             6: (333.089053559821, 333.089053566949)}
FULL_DISC = {0: (71.9145912589741, 365.191485194442),
             1: (75.4345915090768, 387.439411608979),
             2: (100.636767428944, 460.55190305468),
             3: (162.594098752413, 594.379095783622),
             4: (258.626226692684, 788.964251330937),
             5: (387.416328919116, 761.075646769153),
             6: (552.639912171029, 552.639912176722)}


def _ring_checks(maxmode=2):
    """Buiging in het vlak (eerste mode van ND n >= 2) tegen Blevins, en
    beide modes per ND tegen het volledige model (CY-02)."""
    out = []
    for n in range(2, 7):
        out.append(_check("ND%d in het vlak (Blevins)" % n,
                          _mode(maxmode, n, 1), blevins_ring(n), 0.02,
                          "Blevins 1979, tabel 9-2: dunne ring, buiging in "
                          "het vlak, n = ND"))
    for n, refs in sorted(FULL_RING.items()):
        for k, ref in enumerate(refs):
            out.append(_check("ND%d mode %d = volledig model" % (n, k + 1),
                              _mode(maxmode, n, k + 1), ref, 0.01,
                              "volledige ring met modal_3d (CY-02)"))
    return out


def _plate_checks(maxmode=2):
    keys = [((0, 1), 0, 1), ((0, 2), 0, 2), ((1, 1), 1, 1), ((2, 1), 2, 1),
            ((3, 1), 3, 1)]
    return [_check("ND%d mode %d (Leissa %d,%d)" % (nd, k, n, s),
                   _mode(maxmode, nd, k), leissa(n, s), 0.02,
                   "Leissa 1969, tabel 2.1: ingeklemde cirkelvormige plaat")
            for (n, s), nd, k in keys]


def _disc_full_checks(maxmode=2):
    out = []
    for nd, refs in sorted(FULL_DISC.items()):
        for k, ref in enumerate(refs):
            out.append(_check("ND%d mode %d = volledig model" % (nd, k + 1),
                              _mode(maxmode, nd, k + 1), ref, 0.01,
                              "volledige draaiende schijf met modal_3d "
                              "(CY-06)"))
    return out


# ==========================================
# HET HELE WIEL (S10i, instelling cyclic_full_wheel, standaard aan)
# ==========================================
# De bridge zet elke mode op alle N sectoren in Mecway's nummering: knoop
# i van sector j = i + j * n, gedraaid over +j * 360/N om de as (headless
# gemeten aan Mecway's eigen oplossing, Scripts/aster_cyclic.py).
_AXIS = {"x": 0, "y": 1, "z": 2}


def _wheel_xyz(model, axis, nseg):
    """[(x, y, z)] van de N * n knopen van het wiel (index = nid - 1)."""
    i, k = (_AXIS[axis] + 1) % 3, (_AXIS[axis] + 2) % 3
    out = []
    for j in range(nseg):
        a = 2.0 * math.pi * j / nseg
        c, s = math.cos(a), math.sin(a)
        for p in model.nodes:
            q = list(p)
            q[i] = c * p[i] - s * p[k]
            q[k] = s * p[i] + c * p[k]
            out.append(tuple(q))
    return out


class _Grid(object):
    """Knopen op coordinaat zoeken (binnen tol, buurcellen)."""

    def __init__(self, pts, tol):
        self.tol, self.pts, self.cells = tol, pts, {}
        for n, p in enumerate(pts):
            self.cells.setdefault(self._key(p), []).append(n)

    def _key(self, p):
        return tuple(int(math.floor(v / self.tol)) for v in p)

    def near(self, p):
        """Alle punten binnen tol van p (indexen)."""
        k = self._key(p)
        out = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for n in self.cells.get((k[0] + dx, k[1] + dy,
                                             k[2] + dz), ()):
                        if max(abs(a - b) for a, b in
                               zip(self.pts[n], p)) <= self.tol:
                            out.append(n)
        return out


def _vec(step, nids):
    """Verplaatsingsvector (displx, disply, displz per knoop) van een mode."""
    out = []
    for n in nids:
        out += [step["displx"][n], step["disply"][n], step["displz"][n]]
    return out


def _wheel_nodes(r):
    return len(r.step(1)["displx"])


def _wheel_jump(r, axis, on_axis):
    """Grootste verschil in verplaatsing tussen knopen op dezelfde plaats,
    over alle modes, relatief t.o.v. de grootste verplaatsing van de mode:
    on_axis False = de snijvlakken tussen de sectoren, True = de knopen op
    de as (N kopieen van dezelfde knoop)."""
    nids = sorted(r.step(1)["displx"])
    pts = r.wheel_xyz
    size = max(abs(v) for p in pts for v in p)
    grid = _Grid([pts[n - 1] for n in nids], 1e-7 * size)
    i, k = (_AXIS[axis] + 1) % 3, (_AXIS[axis] + 2) % 3
    pairs = set()
    for a in range(len(nids)):
        p = grid.pts[a]
        if (math.hypot(p[i], p[k]) <= 1e-7 * size) != on_axis:
            continue
        for b in grid.near(p):
            if b > a:
                pairs.add((nids[a], nids[b]))
    if not pairs:
        raise ValueError("geen samenvallende knopen in het wiel")
    worst = 0.0
    for _t, step in r.steps:
        top = max(abs(v) for name in ("displx", "disply", "displz")
                  for v in step[name].values())
        for a, b in pairs:
            d = max(abs(step[name][a] - step[name][b])
                    for name in ("displx", "disply", "displz"))
            worst = max(worst, d / top)
    return worst


def _wheel_checks(build, axis, nseg, axis_nodes=False):
    """Controles op het hele wiel van een cyclische case (zelfstandig).
    axis_nodes: het model heeft knopen op de as (extra controle)."""
    model = build()
    n = len(model.nodes)

    def jump(on_axis):
        def value(r):
            r.wheel_xyz = _wheel_xyz(model, axis, nseg)
            return _wheel_jump(r, axis, on_axis)
        return value
    out = [
        Check("hele wiel: knopen met verplaatsing", _wheel_nodes,
              ref=nseg * n, tol_abs=0.5,
              source="Mecway's eigen cyclische oplossing (headless, %d x %d "
                     "knopen)" % (nseg, n), in_series=False),
        Check("hele wiel: sprong op de snijvlakken", jump(False), ref=0.0,
              tol_abs=1e-6, source="doorlopende vorm: samenvallende knopen "
              "van naburige sectoren gelijk (relatief)", in_series=False)]
    if axis_nodes:
        # Mecway's eigen oplossing van CY-04 (headless, 7 okt 2026) heeft
        # op de as in ND 1 dezelfde orde: 1,9e-4 en 9,4e-4 (modes 3, 4);
        # de bridge 1,2e-4 en 6,2e-4. In de andere ND's < 1e-14.
        out.append(Check("hele wiel: sprong op de as", jump(True), ref=0.0,
                         tol_abs=2e-3, source="Mecway's eigen cyclische "
                         "oplossing: tot 9,4e-4 (ND 1, relatief)",
                         in_series=False))
    return out


_WHEEL_CACHE = {}


def _work_root():
    return os.environ.get("VERIF_WORK") or os.path.join(
        tempfile.gettempdir(), "mab_verificatie")


def _wheel_result(case_id, build):
    """Results van het hele wiel van een eerder gedraaide cyclische case
    (werkmap van run_verificatie), met de knoopcoordinaten van het wiel."""
    if case_id not in _WHEEL_CACHE:
        run = os.path.join(_work_root(), case_id, "run.json")
        path = json.load(io.open(run, encoding="utf-8")).get("results")
        if not path:
            raise ValueError("%s heeft geen resultaten" % case_id)
        _WHEEL_CACHE[case_id] = Results(path, build())
    return _WHEEL_CACHE[case_id]


def _subspace_mac(u, basis):
    """||P u||^2 / ||u||^2 met P de projectie op de kolommen van basis
    (Gram-Schmidt): 1 = u ligt in de ruimte van de modes van het volledige
    model (bij een paar of een meervoudige frequentie meer dan een mode)."""
    q = []
    for b in basis:
        v = list(b)
        for e in q:
            d = sum(x * y for x, y in zip(v, e))
            v = [x - d * y for x, y in zip(v, e)]
        norm = math.sqrt(sum(x * x for x in v))
        if norm > 1e-12 * math.sqrt(sum(x * x for x in b)):
            q.append([x / norm for x in v])
    uu = sum(x * x for x in u)
    return sum(sum(x * y for x, y in zip(u, e)) ** 2 for e in q) / uu


def _wheel_mac(cyc_id, build, axis, nseg, mode, rel=1e-4, zero=0.05):
    """
    Waarde voor een controle in de case van het volledige model: 1 - de
    subspace-MAC van mode `mode` van het hele wiel van case cyc_id met de
    modes van het volledige model (r) binnen de frequentieband van die
    mode (relatief rel, of < zero Hz voor starre modes).
    """
    def value(r):
        w = _wheel_result(cyc_id, build)
        if not hasattr(w, "wheel_xyz"):
            w.wheel_xyz = _wheel_xyz(w.model, axis, nseg)
        f = w.modes[mode]
        full = [k for k, g in r.modes.items()
                if (abs(g) < zero if abs(f) < zero
                    else abs(g - f) <= rel * abs(f))]
        if not full:
            raise ValueError("geen mode van het volledige model bij %.6g Hz"
                             % f)
        # knoop van het volledige model bij elke knoop van het wiel
        size = max(abs(v) for p in r.model.nodes for v in p)
        grid = _Grid(r.model.nodes, 1e-7 * size)
        step = w.step(mode)
        nids = sorted(step["displx"])
        match = []
        for nid in nids:
            hit = grid.near(w.wheel_xyz[nid - 1])
            if not hit:
                raise ValueError("knoop %d van het wiel niet in het volledige"
                                 " model" % nid)
            match.append(hit[0] + 1)
        u = _vec(step, nids)
        basis = [_vec(r.step(k), match) for k in sorted(full)]
        return 1.0 - _subspace_mac(u, basis)
    return value


def _mac_checks(cyc_id, build, axis, nseg, freqs_nd):
    """Per mode van het wiel van cyc_id: MAC >= 0,999 met het volledige
    model (deze case), als 1 - MAC <= 1e-3."""
    return [Check("%s mode %d (ND %d): 1 - MAC hele wiel" % (cyc_id, k + 1,
                                                             nd),
                  _wheel_mac(cyc_id, build, axis, nseg, k + 1), ref=0.0,
                  tol_abs=1e-3, source="hele wiel van %s (S10i) tegen de "
                  "modes van dit volledige model, subspace-MAC binnen de "
                  "frequentieband" % cyc_id, in_series=False)
            for k, nd in enumerate(freqs_nd)]


def _full_blevins():
    """Volledige ring: de frequentie het dichtst bij Blevins, n = 2 .. 6."""
    return [Check("n=%d in het vlak (Blevins)" % n, _nearest(blevins_ring(n)),
                  ref=blevins_ring(n), tol=0.02, unit="Hz",
                  source="Blevins 1979, tabel 9-2") for n in range(2, 7)]


# Modellen met het hele wiel (S10i): dezelfde builders in de case en in de
# controles van het volledige model.
RING_Z = ring(axis="z", maxmode=2)
RING_X = ring(axis="x", maxmode=2)
PLATE = disc(nseg=6, nr=16, nt=6, maxmode=2)
DISC_Y = disc(nseg=12, nr=6, nt=2, axis="y", maxmode=2, hub=R0, rpm=3000.0,
              clamp="inner")
ND_12 = [k // 2 for k in range(14)]     # ND per mode: 12 segmenten, 2 per ND

TESTS = [
    Case("CY-01", GROUP, "Vrije ring, sector 30 graden (12 segmenten) om Z",
         "modal_3d", RING_Z,
         _ring_checks() + _mecway_checks(MECWAY["CY-01"]) +
         _wheel_checks(RING_Z, "z", 12),
         note="hex20, 6 per segment; starre modes in ND 0 en 1; het hele "
              "wiel naar Mecway"),
    Case("CY-02", GROUP, "Vrije ring, volledig model (referentie CY-01)",
         "modal_3d", ring(sectors=12, maxmode=40),
         _full_blevins() + _mac_checks("CY-01", RING_Z, "z", 12, ND_12),
         note="zelfde mesh als 12 x CY-01, gewone modale route; MAC met het "
              "hele wiel van CY-01 (draai CY-01 eerst)"),
    Case("CY-03", GROUP, "Vrije ring, sector om X (as gedraaid)",
         "modal_3d", RING_X,
         _ring_checks() + _mecway_checks(MECWAY["CY-03"]) +
         _wheel_checks(RING_X, "x", 12),
         note="bridge draait het model naar de Z-as en terug"),
    Case("CY-04", GROUP, "Ingeklemde cirkelvormige plaat, 6 segmenten, "
         "knopen op de as", "modal_3d", PLATE,
         _plate_checks() + _mecway_checks(MECWAY["CY-04"]) +
         _wheel_checks(PLATE, "z", 6, axis_nodes=True),
         note="a = 0,5 m, h = 0,01 m; hex20 + wedge15 in het midden; "
              "interface AXE"),
    Case("CY-05", GROUP, "Draaiende schijf om Y (3000 r/min), naaf "
         "ingeklemd, 12 segmenten", "modal_3d", DISC_Y,
         _disc_full_checks() + _mecway_checks(MECWAY["CY-05"]) +
         _wheel_checks(DISC_Y, "y", 12),
         note="voorspanning (stress stiffening + spin softening), paren met "
              "een vaste knoop op de naaf"),
    Case("CY-06", GROUP, "Draaiende schijf, volledig model (referentie "
         "CY-05)", "modal_3d",
         disc(nseg=12, nr=6, nt=2, sectors=12, axis="y", maxmode=30,
              hub=R0, rpm=3000.0, clamp="inner"),
         _mecway_checks(MECWAY["CY-06"]) +
         _mac_checks("CY-05", DISC_Y, "y", 12, ND_12),
         note="zelfde mesh als 12 x CY-05, gewone modale route; MAC met het "
              "hele wiel van CY-05 (draai CY-05 eerst)"),
]
