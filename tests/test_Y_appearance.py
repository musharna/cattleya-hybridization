import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "analysis"))
import Y_appearance as Y  # noqa: E402


def _unit(x):
    return x / np.linalg.norm(x, axis=-1, keepdims=True)


def _keys(n):
    """Every offspring has its own two parent photos."""
    return {"ka": [("a", i) for i in range(n)], "kb": [("b", i) for i in range(n)]}


def test_split_recovers_planted_selection_with_zero_transmission():
    rng = np.random.default_rng(0)
    d = 16
    axis = np.zeros(d)
    axis[0] = 1
    pool = rng.normal(size=(300, d))
    pool[:, 0] -= pool[:, 0].mean()  # pool mean exactly 0 on the axis
    pa = rng.normal(size=(50, d))
    pb = rng.normal(size=(50, d))
    pa[:, 0] += 0.8  # breeders chose parents 0.8 above the pool
    pb[:, 0] += 0.8
    off = (
        pa + pb
    ) / 2  # transmission-free: offspring sit exactly on their parents' midpoint
    r = Y.split(
        {1960: {"pool": pool, "off": off, "pa": pa, "pb": pb, **_keys(50)}}, axis, rng, n_boot=200
    )[1960]
    assert abs(r["T"]) < 1e-12 and r["T_ci"][0] > -1e-12 and r["T_ci"][1] < 1e-12
    assert abs(r["S"] - ((pa + pb) / 2)[:, 0].mean()) < 1e-12
    assert r["S_ci"][0] > 0.4  # planted 0.8 (plus noise) is clearly above 0
    assert Y.label(r) == "selection-led"


def test_split_recovers_planted_transmission():
    rng = np.random.default_rng(1)
    d = 8
    axis = np.eye(d)[0]
    pool = rng.normal(size=(400, d))
    pa, pb = (
        pool[:60].copy(),
        pool[60:120].copy(),
    )  # parents drawn from the pool: no selection
    off = (pa + pb) / 2
    off[:, 0] += 1.0  # offspring shifted 1.0 off their midpoint
    r = Y.split(
        {1990: {"pool": pool, "off": off, "pa": pa, "pb": pb, **_keys(60)}}, axis, rng, n_boot=200
    )[1990]
    assert abs(r["T"] - 1.0) < 1e-9
    assert Y.label(r) == "transmission-led", r


def test_label_small_decade_is_unlabelled():
    assert Y.label(
        {"n": 19, "S": 1, "S_ci": [0.5, 1.5], "T": 0, "T_ci": [-0.1, 0.1]}
    ).startswith("unlabelled")


def test_gate_auc_high_for_true_parents_and_chance_for_random():
    rng = np.random.default_rng(2)
    sp = _unit(rng.normal(size=(60, 32)))
    a, b = rng.integers(0, 60, 300), rng.integers(0, 60, 300)
    off = _unit((sp[a] + sp[b]) / 2 + 0.05 * rng.normal(size=(300, 32)))
    hi = Y.gate_auc(off, sp[a], sp[b], sp, rng, reps=50)[0]
    rnd = _unit(rng.normal(size=(300, 32)))
    lo = Y.gate_auc(rnd, sp[a], sp[b], sp, rng, reps=50)[0]
    assert hi > 0.95 and 0.4 < lo < 0.6, (hi, lo)


def test_shuffle_null_flags_chosen_parents_not_random_ones():
    rng = np.random.default_rng(3)
    decade_of = {g: 1950 for g in range(200)}
    proj = {("g", g): float(rng.normal()) for g in range(200)}
    top = sorted(range(200), key=lambda g: -proj[("g", g)])[
        :40
    ]  # breeders chose the 40 highest
    chosen = {
        "pool_mean": 0.0,
        "pa": [("g", g) for g in top[:20]],
        "pb": [("g", g) for g in top[20:]],
    }
    rand = {
        "pool_mean": 0.0,
        "pa": [("g", g) for g in range(20)],
        "pb": [("g", g) for g in range(20, 40)],
    }
    null = Y.shuffle_null({1960: chosen, 1970: rand}, proj, decade_of, rng, reps=200)
    s_chosen = np.mean([proj[("g", g)] for g in top])
    s_rand = np.mean([proj[("g", g)] for g in range(40)])
    assert s_chosen > np.percentile(null[1960], 99)
    assert np.percentile(null[1970], 5) < s_rand < np.percentile(null[1970], 95)


def test_trajectory_axis_points_from_first_to_last_decade():
    means = {D: np.array([D - 1880.0, 0.0, 0.0]) for D in range(1880, 2030, 10)}
    v = Y.trajectory_axis(means)
    assert v[0] > 0.99


def test_gate_random_comparison_is_a_two_species_midpoint():
    # offspring = species centroid, own parents = a random distinct pair: own and random midpoints are
    # exchangeable, so AUC ~ 0.5. A random "pair" that is one species twice sits farther from the centroid
    # than any two-species midpoint and would push AUC far above 0.5.
    rng = np.random.default_rng(4)
    sp = _unit(rng.normal(size=(80, 64)))
    n = 400
    a = rng.integers(0, 80, n)
    b = (a + rng.integers(1, 80, n)) % 80
    off = np.repeat(_unit(sp.mean(0, keepdims=True)), n, 0)
    auc = Y.gate_auc(off, sp[a], sp[b], sp, rng, reps=50)[0]
    assert 0.4 < auc < 0.6, auc


def test_parent_unit_bootstrap_covers_truth_when_parents_are_reused():
    # 200 offspring from pairs of 12 parents drawn from a pool with mean 0 on the axis: E[S] = 0. Parent reuse
    # makes offspring dependent; an offspring-level bootstrap (unique keys per offspring) is far too narrow.
    d, cover_unit, cover_off = 4, 0, 0
    axis = np.eye(d)[0]
    for s in range(60):
        rng = np.random.default_rng(100 + s)
        pool = rng.normal(size=(2000, d))
        pool[:, 0] -= pool[:, 0].mean()
        par = rng.normal(size=(12, d))
        a, b = rng.integers(0, 12, 200), rng.integers(0, 12, 200)
        off = (par[a] + par[b]) / 2 + 0.1 * rng.normal(size=(200, d))
        base = {"pool": pool, "off": off, "pa": par[a], "pb": par[b]}
        unit = Y.split({1: {**base, "ka": list(a), "kb": list(b)}}, axis, rng, n_boot=300)[1]
        offb = Y.split({1: {**base, **_keys(200)}}, axis, rng, n_boot=300)[1]
        cover_unit += unit["S_ci"][0] <= 0 <= unit["S_ci"][1]
        cover_off += offb["S_ci"][0] <= 0 <= offb["S_ci"][1]
    assert cover_unit / 60 >= 0.8 and cover_off / 60 <= 0.5, (cover_unit, cover_off)
