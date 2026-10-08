# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# aster_acoustic_complex.py
#
# Verantwoordelijkheid: complexe akoestische holtemodes met impedantieranden
# (Mecway A30 met <impedance>, bridge_acoustic.build_acoustic_modes). Draait
# BINNEN de Code-Aster-job (CPython 3.11 met numpy en scipy, beide
# meegeleverd met Code-Aster 17.4 in ...\v2025\Python311), NIET in
# IronPython.
#
# Waarom niet CALC_MODES: het kwadratische eigenwaardeprobleem met een
# dempingsmatrix kan CALC_MODES alleen met MATR_AMOR van het type
# matr_asse_depl_r of matr_asse_gene_r (catalogus calc_modes.py v17.4); de
# akoestische matrices zijn matr_asse_pres_c, AFFE_CHAR_ACOU kent geen
# impedantie (alleen PRES_IMPO, VITE_FACE, LIAISON_UNIF) en een complexe Z
# geeft een complexe dempingsmatrix (KNOWN_BUGS 55). Code-Aster maakt
# daarom alleen de matrices (ASSEMBLAGE RIGI_ACOU/MASS_ACOU en per
# impedantie AMOR_ACOU met CELE_R = 1 op dat vlak = int N N^T dA); dit
# script lost het eigenwaardeprobleem met scipy op.
#
# Mecway's probleem (manual v33 par. 6.9.2):
#     (K + lambda C + lambda^2 M) P = 0,   lambda = i w - sigma
#     K = int grad N grad N^T dV,  M = int N N^T / c^2 dV,
#     C = int rho / Z N N^T dA  (Z = specificimpedancer + i
#                                specificimpedancei, Pa.s/m)
#     p(t) = Re(P e^(lambda t)): frequentie w / (2 pi), decay rate sigma.
# Code-Aster's RIGI_ACOU en MASS_ACOU zijn K en M (gemeten: som van M = V /
# c^2 tot 1e-15; zonder rho, docs/Akoestiek.md par. 3.2); AMOR_ACOU met
# CELE_R = 1 is int N N^T dA (gemeten: som = oppervlak 0,01 m2). Gemeten
# (buis 1 x 0,1 x 0,1 m hex20, hard / Z): dezelfde eigenwaarden als Mecway
# 33 headless tot 1e-9 (Z = 2000, 200, 2000 - 1000i) en dezelfde vormen op
# een complexe factor na (spreiding 1e-9).
#
# Oplossing: lineariseren (eerste companion-vorm) A x = lambda B x met
#     A = [[0, I], [-K, -C]],  B = [[I, 0], [0, M]],  x = [P; lambda P]
# en shift-invert rond sigma0: (A - sigma0 B)^-1 B y vraagt per stap een
# oplossing met Q(sigma0) = K + sigma0 C + sigma0^2 M (scipy splu); de
# eigenwaarden met de kleinste |lambda - doel| (scipy.sparse.linalg.eigs,
# ARPACK). Kleine modellen (<= _DENSE_MAX vrijheidsgraden): dicht met
# scipy.linalg.eig.
#
# Keuze van de modes (zoals Mecway, headless nagegaan met maxmode 3/6/8/10
# op buizen en een ruimte): Mecway neemt de 2 x maxmode eigenwaarden met de
# kleinste |lambda| en toont die met frequentie >= 0, oplopend (dus niet
# altijd precies maxmode: buis Z = 200 met maxmode 6 gaf 5 modes, de
# ruimte AI-10 met maxmode 8 gaf 10); modes met 0 Hz (de constante druk
# lambda = 0 zonder drukknoop, niet-oscillerende modes) vallen door
# afrondingsruis in Im(lambda) wisselend wel of niet weg (manual:
# "randomly omits or duplicates 0 Hz rigid body modes"). De bridge doet
# hetzelfde, maar telt een reele eigenwaarde (|Im| ruis) altijd een keer
# mee, als 0 Hz, en een conjugatiepaar op de grens helemaal. Een
# eigenwaarde met negatieve frequentie hoort bij een complexe Z bij de
# toegevoegde impedantie conj(Z) (Re(P e^(lambda t)) = Re(conj(P)
# e^(conj(lambda) t))) en valt weg, zoals bij Mecway. Met een zoekfrequentie
# (shift_hz > 0, eigen keuze van de bridge): de n_modes modes met f >= 0 het
# dichtst bij i 2 pi shift_hz.
#
# Normering (Mecway schaalt willekeurig: Euclidische norm 1 met een
# willekeurige fase, gemeten): fase zo dat P^T M P reeel en positief is,
# dan P^H M P = 1 en de grootste component van Re(P) positief. Bij
# Z -> oneindig is dat precies de reele holtemode van bridge_acoustic
# (p^T M p = 1, zoals Mecway's reele modes).
#
# Uitvoer: een <mode id="k" frequency="f Hz" decayrate="sigma 1/s"> per mode
# met de knoopvariabelen pressure = Re(P) en pressurei = Im(P) (Mecway's
# namen en formaat, headless gezien), via aster_results_export.
# _write_fragment, en per mode een regel
#     BRIDGE_MODE k FREQ=f DECAY=sigma ZETA=sigma/|lambda| Q=|lambda|/(2 sigma)
#                 T60=3 ln(10)/sigma RESIDU=..
# in de .mess (bridge_results_liml.read_mode_summary; T60 en Q alleen bij
# sigma > 0).
#
# Gebruik geen f-strings: zelfde stijl als de rest van de bridge.

import math
import os
import time
import traceback

import numpy as np

import aster_results_export as are

_DENSE_MAX = 300         # vrijheidsgraden waaronder dicht opgelost wordt
_T60 = 3.0 * math.log(10.0)     # 60 dB: 20 log10(e^(-sigma T)) = -60


def _sparse(matrix):
    """Geassembleerde Code-Aster-matrix -> scipy csc (reeel deel)."""
    import scipy.sparse as sp
    vals, rows, cols, dim = matrix.getValuesWithDescription()[:4]
    vals = np.asarray(vals)
    if np.iscomplexobj(vals):
        vals = vals.real
    return sp.csc_matrix((vals, (np.asarray(rows), np.asarray(cols))),
                         shape=(dim, dim))


def _dof_nodes(numbering):
    """[mesh-knoopindex (0-based)] per vrijheidsgraad (alleen PRES)."""
    (nodes, _cmps), dofs = numbering.getEquationNumbering() \
        .getDOFsWithDescription(["PRES"])
    n = numbering.getNumberOfDOFs()
    out = np.full(n, -1, dtype=int)
    out[np.asarray(dofs, dtype=int)] = np.asarray(nodes, dtype=int)
    return out


def _group_nodes(mesh, kind, name):
    """Knoopindices (0-based) van een GROUP_NO of van de knopen van een
    GROUP_MA."""
    if kind == "GROUP_NO":
        return list(mesh.getNodes(name))
    return list(mesh.getNodesFromCells(name))


def _eig_dense(K, C, M, sigma0, k):
    """Alle 2n eigenwaarden van het gelineariseerde probleem (dicht)."""
    import scipy.linalg as sl
    n = K.shape[0]
    Kd, Cd, Md = [np.asarray(x.todense()) for x in (K, C, M)]
    A = np.zeros((2 * n, 2 * n), dtype=complex)
    B = np.zeros((2 * n, 2 * n), dtype=complex)
    A[:n, n:] = np.eye(n)
    A[n:, :n] = -Kd
    A[n:, n:] = -Cd
    B[:n, :n] = np.eye(n)
    B[n:, n:] = Md
    lam, X = sl.eig(A, B)
    keep = np.isfinite(lam)
    return lam[keep], X[:n, keep]


def _eig_sparse(K, C, M, sigma0, k):
    """k eigenwaarden dicht bij sigma0 (shift-invert, ARPACK)."""
    import scipy.sparse.linalg as sla
    n = K.shape[0]
    Kc, Cc, Mc = [x.astype(complex).tocsc() for x in (K, C, M)]
    lu = None
    for factor in (1.0, 1.37, 1.91):
        s = sigma0 * factor
        try:
            lu = sla.splu((Kc + s * Cc + s * s * Mc).tocsc())
            sigma0 = s
            break
        except RuntimeError:
            lu = None       # singulier (shift precies op een eigenwaarde)
    if lu is None:
        raise RuntimeError(are.tr("Q(sigma0) singulier; kies een andere "
                                  "zoekfrequentie.",
                                  "Q(sigma0) singular; choose another "
                                  "search frequency."))
    CsM = (Cc + sigma0 * Mc).tocsr()
    Mr = Mc.tocsr()

    def op(x):
        u, w = x[:n], x[n:]
        uu = lu.solve(-(Mr @ w) - CsM @ u)
        return np.concatenate([uu, u + sigma0 * uu])

    L = sla.LinearOperator((2 * n, 2 * n), matvec=op, dtype=complex)
    ncv = min(2 * n - 1, max(2 * k + 1, 40))
    try:
        mu, X = sla.eigs(L, k=k, which="LM", ncv=ncv)
    except sla.ArpackNoConvergence:
        # tweede poging met een grotere Krylov-ruimte (scipy-documentatie:
        # ncv groter maakt convergentie waarschijnlijker)
        print(are.tr("aster_acoustic_complex: ARPACK niet geconvergeerd; "
                     "tweede poging met ncv = {0}.",
                     "aster_acoustic_complex: ARPACK did not converge; "
                     "second attempt with ncv = {0}.").format(
            min(2 * n - 1, 2 * ncv)))
        mu, X = sla.eigs(L, k=k, which="LM", ncv=min(2 * n - 1, 2 * ncv),
                         maxiter=max(1000, 20 * n))
    lam = sigma0 + 1.0 / mu
    return lam, X[:n, :]


def _normalize(p, M):
    """Fase zo dat p^T M p reeel en positief is, p^H M p = 1, grootste
    component van Re(p) positief (zie de kop)."""
    mp = M @ p
    t = complex(p.dot(mp))
    if abs(t) > 1.0e-12 * float(np.real(np.vdot(p, mp))):
        p = p * np.exp(-0.5j * np.angle(t))
    else:
        k = int(np.argmax(np.abs(p)))
        p = p * np.exp(-1j * np.angle(p[k]))
    h = float(np.real(np.vdot(p, M @ p)))
    if h > 0.0:
        p = p / math.sqrt(h)
    k = int(np.argmax(np.abs(p.real)))
    if p.real[k] < 0.0:
        p = -p
    return p


def _select(lam, X, n_modes, target, w_ref, like_mecway):
    """
    (lambda, kolom)-paren, oplopend in frequentie (zie de kop).
    like_mecway (zoekfrequentie 0): de 2 n_modes eigenwaarden met de
    kleinste |lambda| (een conjugatiepaar dat op de grens valt telt
    helemaal mee) en daarvan die met frequentie >= 0; het aantal is dus
    ongeveer n_modes (bij reele Z n_modes + de helft van het aantal reele
    eigenwaarden). Anders: de n_modes met frequentie >= 0 het dichtst bij
    target. Een reele eigenwaarde (|Im| ruis) telt een keer, als 0 Hz.
    Returns (lijst, aantal kandidaten met f >= 0).
    """
    scale = max(w_ref, float(np.max(np.abs(lam))) if len(lam) else 1.0)
    tol = 1.0e-7 * scale
    order = sorted(range(len(lam)), key=lambda j: abs(lam[j] - target))
    if like_mecway and len(order) > 2 * n_modes:
        cut = abs(lam[order[2 * n_modes - 1]] - target)
        order = [j for j in order
                 if abs(lam[j] - target) <= cut * (1.0 + 1.0e-8) + tol]
    out = []
    n_pos = 0
    for j in order:
        l = complex(lam[j])
        if abs(l.imag) <= tol:
            l = complex(l.real, 0.0)        # reeel: 0 Hz, een keer
        elif l.imag < 0.0:
            continue                        # negatieve frequentie
        n_pos += 1
        out.append((l, j))
        if not like_mecway and len(out) >= n_modes:
            break
    out.sort(key=lambda t: (t[0].imag, -t[0].real))
    return out, n_pos


def solve_modes(rigi, mass, damping, numbering, mesh, out_path, liml_path,
                n_modes, dirichlet=(), shift_hz=0.0, free_shift_hz=-1.0,
                c_ref=340.0, size_ref=1.0, tolerance=1.0e-6,
                precision="double"):
    """
    Los de complexe holtemodes op en schrijf ze als Mecway-<results>-
    fragment (zie de kop).

    rigi, mass    : RIGI_ACOU / MASS_ACOU (ASSEMBLAGE, zonder CHARGE).
    damping       : [(AMOR_ACOU met CELE_R = 1 op het vlak, rho / Z
                    (complex), omschrijving)] per impedantierand.
    numbering     : de NUME_DDL van de matrices.
    dirichlet     : [("GROUP_NO" | "GROUP_MA", naam)] drukknopen (p = 0).
    shift_hz      : 0 = de modes met de kleinste |lambda| (zoals Mecway);
                    > 0 = de modes het dichtst bij i 2 pi shift_hz.
    free_shift_hz : negatieve hulpfrequentie [Hz] voor de factorisatie bij
                    shift_hz = 0 (bridge_acoustic._free_shift).
    c_ref, size_ref: geluidssnelheid en grootte van het model (toleranties).
    Returns dict met statistiek. Bij een fout: melding in de .mess en de
    uitzondering opnieuw (de job stopt; er komt geen resultaatfragment, de
    bridge meldt dat).
    """
    t0 = time.time()
    try:
        out = _solve(rigi, mass, damping, numbering, mesh, out_path,
                     liml_path, int(n_modes), dirichlet, float(shift_hz),
                     float(free_shift_hz), float(c_ref), float(size_ref),
                     float(tolerance), precision, t0)
    except Exception as e:
        traceback.print_exc()
        print(are.tr("aster_acoustic_complex: FOUT: {0}",
                     "aster_acoustic_complex: ERROR: {0}").format(e))
        raise
    if out is None:
        raise RuntimeError(are.tr(
            "aster_acoustic_complex: geen complexe modes (zie de FOUT-regel "
            "hierboven).",
            "aster_acoustic_complex: no complex modes (see the ERROR line "
            "above)."))
    return out


def _solve(rigi, mass, damping, numbering, mesh, out_path, liml_path,
           n_modes, dirichlet, shift_hz, free_shift_hz, c_ref, size_ref,
           tolerance, precision, t0):
    import scipy.sparse as sp
    model = are.read_liml_model(liml_path)
    n_liml = model["coords"].shape[0]
    err = are._check_node_order(mesh, model["coords"])
    if err:
        print(are.tr("aster_acoustic_complex: FOUT: {0}",
                     "aster_acoustic_complex: ERROR: {0}").format(err))
        return None
    K = _sparse(rigi)
    M = _sparse(mass)
    C = sp.csc_matrix(K.shape, dtype=complex)
    for amor, coef, label in damping:
        A = _sparse(amor)
        C = C + complex(coef) * A
        print("aster_acoustic_complex: impedantie {0}: rho/Z = {1}, "
              "oppervlak {2} m2".format(label, complex(coef),
                                        are._num(A.sum())))
    dof_node = _dof_nodes(numbering)
    fixed = set()
    for kind, name in dirichlet:
        fixed.update(_group_nodes(mesh, kind, name))
    free = np.array([d for d in range(len(dof_node))
                     if dof_node[d] >= 0 and dof_node[d] not in fixed],
                    dtype=int)
    if not len(free):
        print(are.tr("aster_acoustic_complex: FOUT: geen vrije drukknopen.",
                     "aster_acoustic_complex: ERROR: no free pressure "
                     "nodes."))
        return None
    K, C, M = [x[free][:, free].tocsc() for x in (K, C, M)]
    n = K.shape[0]
    w_ref = c_ref / max(size_ref, 1.0e-30)
    if shift_hz > 0.0:
        target = complex(0.0, 2.0 * math.pi * shift_hz)
        sigma0 = target
    else:
        target = 0j
        sigma0 = complex(2.0 * math.pi * min(free_shift_hz, -1.0e-6), 0.0)
    like_mecway = shift_hz <= 0.0
    # Ruim meer dan 2 n_modes: de grens van de selectie (rond 0, de shift
    # ligt op sigma0) en een gesplitst conjugatiepaar moeten erin vallen.
    want = 2 * n_modes + 8
    found = []
    tries = 0
    while True:
        tries += 1
        if 2 * n <= _DENSE_MAX:
            lam, X = _eig_dense(K, C, M, sigma0, want)
        else:
            k = min(want, 2 * n - 2)
            lam, X = _eig_sparse(K, C, M, sigma0, k)
        found, n_pos = _select(lam, X, n_modes, target, w_ref, like_mecway)
        if (like_mecway or len(found) >= n_modes or 2 * n <= _DENSE_MAX or
                want >= 2 * n - 2 or tries >= 3):
            break
        want *= 2
    if not found:
        print(are.tr("aster_acoustic_complex: FOUT: geen modes gevonden.",
                     "aster_acoustic_complex: ERROR: no modes found."))
        return None

    nodes = dof_node[free]
    rows = np.sort(np.unique(dof_node[dof_node >= 0]))
    rows = rows[rows < n_liml]
    pos = np.searchsorted(rows, nodes)
    modes = []
    # Relatief residu ||Q(l) p|| / ((|K| + |l| |C| + |l|^2 |M|) ||p||) met
    # 1-normen (backward error van het kwadratische probleem; ook zinvol
    # voor l = 0, waar K p alleen afrondingsruis is).
    import scipy.sparse.linalg as sla
    nk, nc, nm = [float(sla.norm(x, 1)) for x in (K, C, M)]
    for l, j in found:
        p = _normalize(np.asarray(X[:, j], dtype=complex), M)
        r = (K @ p + l * (C @ p) + l * l * (M @ p))
        den = (nk + abs(l) * nc + abs(l) ** 2 * nm) * np.linalg.norm(p)
        resid = float(np.linalg.norm(r) / den) if den > 0.0 else 0.0
        full = np.zeros(len(rows), dtype=complex)
        ok = nodes < n_liml
        full[pos[ok]] = p[ok]
        modes.append((l, full, resid))

    indexes = list(range(1, len(modes) + 1))
    freqs = dict((i, m[0].imag / (2.0 * math.pi))
                 for i, m in zip(indexes, modes))
    decays = dict((i, _clean(-m[0].real, w_ref))
                  for i, m in zip(indexes, modes))
    ids = (rows + 1).astype("<i4")
    data = _Modes([m[1] for m in modes])
    variables = [are._Var("pressure", False, ids, lambda d: d.current().real),
                 are._Var("pressurei", False, ids,
                          lambda d: d.current().imag)]
    elset_nodes, free_nodes = are._group_nodes(model, n_liml)
    for var in variables:
        are._assign_groups(var, model, elset_nodes, free_nodes)
    tmp = out_path + ".tmp"
    with open(tmp, "wb") as f:
        are._write_fragment(f, variables, data, list(range(len(modes))),
                            None, precision, model,
                            ("frequency", dict((i - 1, freqs[i])
                                               for i in indexes),
                             dict((i - 1, decays[i]) for i in indexes)))
        size = f.tell()
    os.replace(tmp, out_path)
    worst = 0.0
    for i, (l, _p, resid) in zip(indexes, modes):
        worst = max(worst, resid)
        sigma = decays[i]
        parts = ["BRIDGE_MODE {0} FREQ={1} DECAY={2}".format(
            i, are._num(freqs[i]), are._num(sigma))]
        if abs(l) > 0.0:
            parts.append("ZETA={0}".format(are._num(sigma / abs(l))))
        if sigma > 0.0:
            parts.append("Q={0} T60={1}".format(
                are._num(abs(l) / (2.0 * sigma)), are._num(_T60 / sigma)))
        parts.append("RESIDU={0:.3g}".format(resid))
        print(" ".join(parts))
    if worst > tolerance:
        print("BRIDGE_COMPLEX_RESIDU max={0:.3g} grens={1:.3g}".format(
            worst, tolerance))
    print("BRIDGE_COMPLEX_ACOU MODES={0} DOF={1} DRUKKNOPEN={2} "
          "IMPEDANTIES={3} KANDIDATEN={4}".format(
              len(modes), n, len(dof_node) - n, len(damping), n_pos))
    print("aster_acoustic_complex: {0} geschreven: {1} complexe mode(s), "
          "{2:.1f} MB in {3:.2f} s.".format(out_path, len(modes),
                                            size / 1.0e6, time.time() - t0))
    return {"stappen": len(modes), "bytes": size}


def _clean(x, w_ref):
    """Numerieke nul (|x| < 1e-9 w_ref) als 0."""
    return 0.0 if abs(x) < 1.0e-9 * w_ref else float(x)


class _Modes(object):
    """Complexe druk per mode (rij i = Mecway-knoop ids[i]); dezelfde rol
    als aster_results_export._Fields voor _write_fragment."""

    def __init__(self, pressures):
        self.pressures = pressures
        self.index = None

    def set_index(self, index):
        self.index = index

    def current(self):
        return self.pressures[self.index]
