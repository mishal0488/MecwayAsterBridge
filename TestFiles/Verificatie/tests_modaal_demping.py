# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Modaal met demping: complexe modes (Mecway V30 + instelling modal_complex
-> CALC_MODES met MATR_AMOR, Scripts/bridge_modal_damping.py). Toelichting:
docs/Verificatie.md par. 7 (MC).

Grootheden per mode (aster_results_export._print_complex_modes, regels
BRIDGE_MODE in de .mess): de gedempte frequentie f_d (ook in <mode
frequency>), de dempingsverhouding zeta (AMOR_REDUIT) en f_n = |s| / 2 pi.

Referenties:
  - 1-vrijheidsgraadsysteem k, m, c: w_n = sqrt(k/m), zeta = c / (2
    sqrt(k m)), w_d = w_n sqrt(1 - zeta^2) (Clough & Penzien, Dynamics of
    Structures, 2e druk 1993, hoofdstuk 2); Rayleigh C = alpha M + beta K:
    zeta = alpha / (2 w) + beta w / 2 (Mecway manual v33 par. 9.26;
    Clough & Penzien par. 12-3);
  - 2 vrijheidsgraden met een demper op een massa (niet-proportioneel): de
    wortels s van det(s^2 M + s C + K) = 0 (hier met Durand-Kerner), f_d =
    Im(s) / 2 pi, zeta = -Re(s) / |s|;
  - hex20-console met Rayleigh: f_n = Mecway's eigen (ongedempte)
    frequenties (Mecway 33 headless, mw_solve.py: Mecway negeert de
    Rayleigh-demping in een Frequency-analyse) en zeta_n uit de formule
    hierboven met die f_n.
"""
import math
from verif_lib import Case, Check, Model, linemesh, hexmesh, box

GROUP = "Modaal met demping"
SETTINGS = {"modal_complex": True}

# 1-vrijheidsgraad: grond A (0,0,0) - veer k - massa (1,0,0) - demper c -
# grond B (2,0,0); de massaknoop alleen vrij in X.
K1, M1 = 1000.0, 2.0
W1 = math.sqrt(K1 / M1)                 # 22,36 rad/s
F1 = W1 / (2.0 * math.pi)               # 3,559 Hz
MASS_X = (1.0, 0.0, 0.0)


def _spring(m, mat, k):
    m.mat(mat, '<mechanical type="Spring" springconstant="%r N/m" />' % k)


def _damper(m, mat, c):
    m.mat(mat, '<mechanical type="Damper" dampingcoefficient="%r N.s/m" />'
          % c)


def _rayleigh(m, alpha, beta):
    m.add('<rayleighdamping selection=""><alpha>%r rad/s</alpha>'
          '<beta>%r s</beta></rayleighdamping>' % (alpha, beta))


def _only_x(m, name, nids):
    """Knopen alleen vrij in X (Mecway: displacement 0 in Y en Z)."""
    m.node_selection(name, nids)
    for d in ((0, 1, 0), (0, 0, 1)):
        m.add('<displacement selection="%s"><value>0 m</value>'
              '<x>%d</x><y>%d</y><z>%d</z></displacement>' % ((name,) + d))


def build_sdof(c=None, alpha=0.0, beta=0.0):
    def build():
        m = Model("V30", attrs={"maxmode": 1})
        m.elset("Veer", "Veer")
        a, n = m.node(0, 0, 0), m.node(1, 0, 0)
        m.elem("Veer", "line2", [a, n])
        _spring(m, "Veer", K1)
        ground = [a]
        if c:
            m.elset("Demper", "Demper")
            b = m.node(2, 0, 0)
            m.elem("Demper", "line2", [n, b])
            _damper(m, "Demper", c)
            ground.append(b)
        m.node_selection("Grond", ground)
        m.add('<fixed selection="Grond" />')
        _only_x(m, "Massa", [n])
        m.add('<mass2 selection="Massa"><value>%r kg</value></mass2>' % M1)
        if alpha or beta:
            _rayleigh(m, alpha, beta)
        return m
    return build


def _zeta(r, k):
    return r.mode_param(k, "AMOR_REDUIT")


def _fn(r, k):
    return r.mode_param(k, "FREQ_NAT")


def _sdof_checks(zeta, source):
    return [
        Check("zeta", lambda r: _zeta(r, 1), ref=zeta, tol=0.01, unit="-",
              source=source),
        Check("f_d (gedempt, <mode frequency>)", lambda r: r.mode(1),
              ref=F1 * math.sqrt(1.0 - zeta ** 2), tol=0.01, unit="Hz",
              source="f_d = f_n sqrt(1 - zeta^2), f_n = sqrt(k/m) / 2 pi "
                     "(Clough & Penzien hfst. 2)"),
        Check("f_n = |s| / 2 pi", lambda r: _fn(r, 1), ref=F1, tol=0.01,
              unit="Hz", source="f_n = sqrt(k/m) / 2 pi"),
        Check("vorm reeel (|Im| / |Re|)", lambda r: r.mode_param(1, "IMAG"),
              ref=0.0, tol_abs=1e-6, unit="-", in_series=False,
              source="1 vrijheidsgraad: de vorm is reeel"),
    ]


TESTS = []
for _id, _z in (("MC-01a", 0.1), ("MC-01b", 0.5)):
    _c = 2.0 * _z * math.sqrt(K1 * M1)
    TESTS.append(Case(
        _id, GROUP, "Massa-veer-demper, zeta = %g (Mecway Damper)" % _z,
        "modal_3d", build_sdof(c=_c),
        _sdof_checks(_z, "zeta = c / (2 sqrt(k m)), c = %.6g N.s/m "
                     "(Clough & Penzien hfst. 2)" % _c),
        settings=SETTINGS, timeout=300,
        note="Veer %g N/m, mass2 %g kg, Damper %.4g N.s/m (line2)."
             % (K1, M1, _c)))

# Rayleigh alpha (C = alpha M) en beta (C = beta K) op hetzelfde systeem
# zonder demper: zeta = alpha / (2 w) en zeta = beta w / 2.
_ZR = 0.2
TESTS.append(Case(
    "MC-02a", GROUP, "Massa-veer, Rayleigh alpha (zeta = 0,2)", "modal_3d",
    build_sdof(alpha=2.0 * _ZR * W1),
    _sdof_checks(_ZR, "zeta = alpha / (2 w_n), alpha = 2 zeta w_n "
                 "(Mecway manual par. 9.26)"),
    settings=SETTINGS, timeout=300,
    note="C = alpha M op de knoopmassa (mass2)."))
TESTS.append(Case(
    "MC-02b", GROUP, "Massa-veer, Rayleigh beta (zeta = 0,2)", "modal_3d",
    build_sdof(beta=2.0 * _ZR / W1),
    _sdof_checks(_ZR, "zeta = beta w_n / 2, beta = 2 zeta / w_n "
                 "(Mecway manual par. 9.26)"),
    settings=SETTINGS, timeout=300,
    note="C = beta K op de veer (Spring)."))


# ------------------------------------------------------------------
# MC-03 Twee vrijheidsgraden, niet-proportionele demping: grond - k1 - m1 -
# k2 - m2 - demper c - grond, alles langs X. Complexe modes met een
# echt complexe vorm (IMAG > 0).
# ------------------------------------------------------------------
K31, K32, M31, M32, C3 = 2000.0, 1000.0, 2.0, 1.0, 30.0


def _poly_mul(p, q):
    out = [0.0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return out


def _roots(coef):
    """Wortels van een polynoom (coef: hoogste macht eerst), Durand-Kerner."""
    a = [c / coef[0] for c in coef]
    n = len(a) - 1
    z = [complex(0.4, 0.9) ** k for k in range(n)]
    for _it in range(500):
        new = []
        for i, zi in enumerate(z):
            num = sum(c * zi ** (n - k) for k, c in enumerate(a))
            den = 1.0 + 0j
            for j, zj in enumerate(z):
                if j != i:
                    den *= (zi - zj)
            new.append(zi - num / den)
        if max(abs(x - y) for x, y in zip(new, z)) < 1e-14:
            z = new
            break
        z = new
    return z


def _mc03_modes():
    """[(f_d, zeta)] oplopend in f_d: det(s^2 M + s C + K) = (m1 s^2 + k1 +
    k2)(m2 s^2 + c s + k2) - k2^2 = 0."""
    p = _poly_mul([M31, 0.0, K31 + K32], [M32, C3, K32])
    p[-1] -= K32 ** 2
    out = []
    for s in _roots(p):
        if s.imag > 0:
            out.append((s.imag / (2 * math.pi), -s.real / abs(s)))
    return sorted(out)


MC03 = _mc03_modes()


def build_mc03():
    m = Model("V30", attrs={"maxmode": 2})
    m.elset("Veer1", "Veer1")
    m.elset("Veer2", "Veer2")
    m.elset("Demper", "Demper")
    g0, n1, n2, g3 = [m.node(x, 0, 0) for x in (0, 1, 2, 3)]
    m.elem("Veer1", "line2", [g0, n1])
    m.elem("Veer2", "line2", [n1, n2])
    m.elem("Demper", "line2", [n2, g3])
    _spring(m, "Veer1", K31)
    _spring(m, "Veer2", K32)
    _damper(m, "Demper", C3)
    m.node_selection("Grond", [g0, g3])
    m.add('<fixed selection="Grond" />')
    _only_x(m, "Massas", [n1, n2])
    for n, mass in ((n1, M31), (n2, M32)):
        m.node_selection("M%d" % n, [n])
        m.add('<mass2 selection="M%d"><value>%r kg</value></mass2>'
              % (n, mass))
    return m


_SRC3 = ("wortels van det(s^2 M + s C + K) = 0, f_d = Im(s) / 2 pi, zeta = "
         "-Re(s) / |s| (k1 = %g, k2 = %g N/m, m1 = %g, m2 = %g kg, c = %g "
         "N.s/m op m2)" % (K31, K32, M31, M32, C3))
TESTS.append(Case(
    "MC-03", GROUP, "Twee massa's, demper op een massa (niet-proportioneel)",
    "modal_3d", build_mc03,
    [Check("f_d mode %d" % (k + 1), lambda r, k=k: r.mode(k + 1),
           ref=MC03[k][0], tol=1e-3, unit="Hz", source=_SRC3)
     for k in range(2)] +
    [Check("zeta mode %d" % (k + 1), lambda r, k=k: _zeta(r, k + 1),
           ref=MC03[k][1], tol=1e-3, unit="-", source=_SRC3)
     for k in range(2)] +
    [Check("vorm complex (|Im| / |Re| > 0,01)",
           lambda r: min(1.0, r.mode_param(1, "IMAG") / 0.01), ref=1.0,
           tol=1e-9, unit="-", in_series=False,
           source="niet-proportionele demping: de vorm is niet reeel")],
    settings=SETTINGS, timeout=300,
    note="Veren en demper (line2), mass2; de demper zit alleen aan m2."))


# ------------------------------------------------------------------
# MC-04 Hex20-console 1 x 0,1 x 0,05 m (10 x 2 x 1, als MD-05), 6 modes,
# Rayleigh alpha = 20 1/s, beta = 1e-4 s: zeta_n = alpha / (2 w_n) + beta
# w_n / 2 per mode, met f_n van Mecway (headless, Rayleigh genegeerd).
# ------------------------------------------------------------------
ALPHA4, BETA4 = 20.0, 1.0e-4
# Mecway 33 headless (mw_solve.py) op modellen/MC-04.liml, 7 okt 2026.
MECWAY4 = (42.2805978662499, 83.6447249461647, 262.851054349367,
           502.193289423994, 619.279712262026, 727.614987544409)


def build_mc04():
    m = Model("V30", attrs={"maxmode": 6})
    m.elset("Balk", "Staal")
    msh = hexmesh(m, "Balk", 10, 2, 1, box(0, 1.0, 0, 0.1, 0, 0.05), order=2)
    m.mat("Staal", '<mechanical type="Isotropic" youngsmodulus="210000 MPa" '
          'poissonratio="0.3" />',
          u'<density type="Simple" density="7850 kg/m³" />')
    m.face_selection("Inklemming", msh["faces"]["u0"])
    m.add('<fixed selection="Inklemming" />')
    _rayleigh(m, ALPHA4, BETA4)
    return m


def _zeta_ray(f):
    w = 2.0 * math.pi * f
    return ALPHA4 / (2.0 * w) + BETA4 * w / 2.0


def _mc04_checks():
    out = []
    src = ("zeta_n = alpha / (2 w_n) + beta w_n / 2 (Mecway manual par. "
           "9.26), w_n = 2 pi f_n(Mecway 33 headless, ongedempt)")
    for k, f in enumerate(MECWAY4):
        z = _zeta_ray(f)
        out += [
            Check("zeta mode %d" % (k + 1), lambda r, k=k: _zeta(r, k + 1),
                  ref=z, tol=0.01, unit="-", source=src),
            Check("f_d mode %d" % (k + 1), lambda r, k=k: r.mode(k + 1),
                  ref=f * math.sqrt(1.0 - z ** 2), tol=0.01, unit="Hz",
                  source="f_d = f_n sqrt(1 - zeta^2), f_n Mecway"),
            Check("f_n mode %d t.o.v. Mecway" % (k + 1),
                  lambda r, k=k: _fn(r, k + 1), ref=f, tol=0.001, unit="Hz",
                  in_series=False, source="Mecway 33 headless (ongedempt)"),
        ]
    return out


TESTS.append(Case(
    "MC-04", GROUP, "Hex20-console, Rayleigh alpha en beta, 6 modes",
    "modal_3d", build_mc04, _mc04_checks(), settings=SETTINGS, timeout=900,
    note="10 x 2 x 1 hex20, ingeklemd; alpha = 20 1/s, beta = 1e-4 s."))


# ------------------------------------------------------------------
# MC-05 Als MC-04 zonder oplegging (vrij-vrij), 10 modes gevraagd. De 6
# starre modes zijn met demping niet oscillerend (s = 0 en s = -alpha,
# reeel) en komen niet mee (Code-Aster alarm ALGELINE5_79): mode 1 is de
# eerste elastische mode. Mecway (headless, ongedempt): modes 7, 8, ...
# ------------------------------------------------------------------
# Mecway 33 headless (mw_solve.py) op modellen/MC-05.liml, 7 okt 2026:
# modes 1-6 0 Hz, dan de elastische modes 7-10.
MECWAY5 = (264.040844874676, 514.076767391355, 720.525093226985,
           1222.29931334231)


def build_mc05():
    m = build_mc04()
    m.items = [i for i in m.items if not i.startswith("  <fixed")]
    m.attrs["maxmode"] = 10
    return m


def _mc05_checks():
    out = []
    for k, f in enumerate(MECWAY5):
        z = _zeta_ray(f)
        out += [
            Check("zeta elastische mode %d" % (k + 1),
                  lambda r, k=k: _zeta(r, k + 1), ref=z, tol=0.01, unit="-",
                  source="zeta_n = alpha / (2 w_n) + beta w_n / 2, f_n = "
                         "Mecway mode %d (headless, ongedempt)" % (k + 7)),
            Check("f_n mode %d t.o.v. Mecway mode %d" % (k + 1, k + 7),
                  lambda r, k=k: _fn(r, k + 1), ref=f, tol=0.001,
                  unit="Hz", in_series=False,
                  source="Mecway 33 headless (ongedempt, starre modes 1-6)"),
        ]
    return out


TESTS.append(Case(
    "MC-05", GROUP, "Hex20-blok vrij-vrij, Rayleigh: starre modes vallen "
    "weg", "modal_3d", build_mc05, _mc05_checks(), settings=SETTINGS,
    timeout=900, note="Als MC-04 zonder inklemming, maxmode 10."))
