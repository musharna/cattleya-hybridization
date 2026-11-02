import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "analysis"))
import prereg_sims as ps  # noqa: E402

RNG = np.random.default_rng(0)
T = np.repeat(np.arange(1900, 2021), 20).astype(float)


def test_hinge_F_matches_direct_least_squares():
    y = RNG.normal(size=(len(T), 2)) + 0.01 * (T - 1900)[:, None]
    F, b2 = ps.hinge_F(y, ps.hinge_basis(T))
    b = 1947
    X = np.column_stack([np.ones_like(T), T, np.clip(T - b, 0, None)])
    L = np.column_stack([np.ones_like(T), T])
    sse = lambda A: ((y[:, 0] - A @ np.linalg.lstsq(A, y[:, 0], rcond=None)[0]) ** 2).sum()
    direct = (sse(L) - sse(X)) / (sse(X) / (len(T) - 3))
    assert abs(F[list(ps.BREAKS).index(b), 0] - direct) < 1e-6 * max(1, direct)
    X60 = np.column_stack([np.ones_like(T), T, np.clip(T - 1960, 0, None)])
    assert abs(b2[0] - np.linalg.lstsq(X60, y[:, 0], rcond=None)[0][2]) < 1e-9


def test_planted_1960_hinge_passes_and_linear_twin_rarely_does():
    base = (10 + 0.1 * (T - 1900))[:, None]
    noise = RNG.normal(0, 1, (len(T), 200))
    basis = ps.hinge_basis(T)
    planted = ps.b_pass(*ps.hinge_F(base + 0.3 * np.clip(T - 1960, 0, None)[:, None] + noise, basis))
    twin = ps.b_pass(*ps.hinge_F(base + noise, basis))
    assert planted.mean() > 0.95 and twin.mean() < 0.15


def test_placebo_break_does_not_pass():
    y = 10 + 0.3 * np.clip(T - 1980, 0, None) + RNG.normal(0, 1, len(T))
    F, b2 = ps.hinge_F(y[:, None], ps.hinge_basis(T))
    assert not ps.b_pass(F, b2)[0]


def test_irls_recovers_planted_logistic():
    n = 20000
    X = np.column_stack([np.ones(n), RNG.integers(0, 2, n), RNG.normal(size=n)])
    beta = np.array([-0.4, np.log(0.5), 0.8])
    y = (RNG.random(n) < 1 / (1 + np.exp(-X @ beta))).astype(float)
    b, se = ps.irls(X, y)
    assert np.all(np.abs(b - beta) < 4 * se) and se[1] < 0.05


def test_negative_1960_slope_change_does_not_pass():
    y = 10 + 0.1 * (T - 1900) - 0.3 * np.clip(T - 1960, 0, None) + RNG.normal(0, 1, len(T))
    F, b2 = ps.hinge_F(y[:, None], ps.hinge_basis(T))
    i60 = list(ps.BREAKS).index(1960)
    assert F[i60, 0] > F[ps.PLACEBO, 0].max() and not ps.b_pass(F, b2)[0]
