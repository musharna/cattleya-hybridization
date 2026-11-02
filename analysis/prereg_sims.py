"""Planted-truth power simulations for B (generation-interval break) and D (fertility ladder).

What: runs the pre-registered statistic on outcomes simulated over the real register's structure
(child years, dated grex parents per child, primaries' pairs, rank classes, parent popularity).
When: before opening the real B or D statistic (prereg/B_D_sim_first.md).
Output: per-cell pass rates; the open/closed decision follows the frozen bars. The real outcome
is never read here.
Test: tests/test_prereg_sims.py (hinge statistic and IRLS against known-answer cases).
Usage: python analysis/prereg_sims.py HYBRIDS.csv SPECIES.csv OUT_JSON
"""

import itertools
import json
import os
import signal
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import founders as fo  # noqa: E402
import register_qc as rq  # noqa: E402

BREAKS = np.arange(1930, 1991)
PLACEBO = np.abs(BREAKS - 1960) > 5


def hinge_basis(t):
    """Per break year b: orthonormal basis of [1, t, (t-b)+] and the β2 row of its pseudo-inverse."""
    t = np.asarray(t, float)
    Q, beta2 = [], {}
    for b in BREAKS:
        X = np.column_stack([np.ones_like(t), t - 1960, np.clip(t - b, 0, None)])
        q, _ = np.linalg.qr(X)
        Q.append(q)
        if b == 1960:
            beta2 = np.linalg.pinv(X)[2]
    lin, _ = np.linalg.qr(np.column_stack([np.ones_like(t), t - 1960]))
    return np.stack(Q), lin, beta2


def hinge_F(Y, basis):
    """Y: n × reps. Returns F (len(BREAKS) × reps) and β2 at 1960 (reps)."""
    Q, lin, beta2 = basis
    n = Y.shape[0]
    yy = (Y * Y).sum(0)
    sse_lin = yy - ((lin.T @ Y) ** 2).sum(0)
    sse_h = yy[None, :] - (np.einsum("bnk,nr->bkr", Q, Y) ** 2).sum(1)
    return (sse_lin[None, :] - sse_h) / (sse_h / (n - 3)), beta2 @ Y


def b_pass(F, b2):
    i60 = int(np.where(BREAKS == 1960)[0][0])
    return (F[i60] > F[PLACEBO].max(0)) & (b2 > 0)


def sim_b(t, k, a, b1, c, delta, reps, rng):
    t = np.asarray(t, float)
    mu = np.clip(
        a + b1 * (t - 1900) + c * (t - 1900) ** 2 + delta * np.clip(t - 1960, 0, None),
        1,
        None,
    )
    out = np.zeros((len(t), reps))
    for m in np.unique(k):
        rows = np.where(k == m)[0]
        ages = rng.gamma(
            2.0, (mu[rows] / 2.0)[:, None, None], size=(len(rows), reps, m)
        )
        out[rows] = ages.mean(2)
    return out


def irls(X, y, iters=50):
    """Logistic regression by IRLS: returns (beta, se)."""
    beta = np.zeros(X.shape[1])
    for _ in range(iters):
        eta = np.clip(X @ beta, -30, 30)
        p = 1 / (1 + np.exp(-eta))
        W = p * (1 - p) + 1e-12
        H = X.T @ (X * W[:, None])
        step = np.linalg.solve(H, X.T @ (y - p))
        beta = beta + step
        if np.abs(step).max() < 1e-10:
            break
    return beta, np.sqrt(np.diag(np.linalg.inv(H)))


def d_design(diff, pop1, pop2, decade):
    decs = sorted(set(decade))[1:]
    cols = [np.ones(len(diff)), np.asarray(diff, float), np.log1p(pop1), np.log1p(pop2)]
    cols += [np.asarray([d == x for d in decade], float) for x in decs]
    return np.column_stack(cols)


def sim_d(X, gamma, kappa, reps, rng, base=0.4):
    diff, lp = X[:, 1], X[:, 2] + X[:, 3]
    lp = lp - lp.mean()
    alpha = np.log(base / (1 - base)) - gamma * diff.mean()
    p = 1 / (1 + np.exp(-(alpha + gamma * diff + kappa * lp)))
    hits = 0
    for _ in range(reps):
        y = (rng.random(len(p)) < p).astype(float)
        beta, se = irls(X, y)
        z = beta[1] / se[1]
        sig = abs(z) > 1.959964
        hits += sig and (gamma == 0 or np.sign(beta[1]) == np.sign(gamma))
    return hits / reps


def structure(hybrids_csv, species_csv, here):
    records, species = rq.load(hybrids_csv, species_csv)
    ped = fo.build(records, species)
    years = {i: r["year"] for i, r in enumerate(records) if r["year"] is not None}
    tb, kb = [], []
    for i, (y, ps) in ped.items():
        m = sum(
            1 for k, t in ps if k == "g" and t in years
        )  # dated grex parents: structure only
        if y is not None and m:
            tb.append(y)
            kb.append(m)
    traits = {
        r["Binomial"]: r
        for r in __import__("csv").DictReader(
            open(os.path.join(here, "..", "data", "species_traits.csv"))
        )
    }
    prim = []
    for i, (y, ps) in ped.items():
        if y is None or [k for k, _ in ps] != ["sp", "sp"]:
            continue
        a, b = ps[0][1], ps[1][1]
        if a in traits and b in traits and a != b:
            prim.append(
                (
                    int(traits[a]["LastSubgenus"] != traits[b]["LastSubgenus"]),
                    int(traits[a]["primary_T1"]),
                    int(traits[b]["primary_T1"]),
                    y // 10 * 10,
                )
            )
    return np.array(tb), np.array(kb), prim


def main(hybrids_csv, species_csv, out_json, reps_b=200, reps_d=500):
    rng = np.random.default_rng(20261025)
    t, k, prim = structure(hybrids_csv, species_csv, os.path.dirname(__file__))
    basis = hinge_basis(t)
    res = {
        "B": {"n_children": int(len(t)), "cells": []},
        "D": {"n_primaries": len(prim), "cells": []},
    }
    for a, b1, c, delta in itertools.product(
        (10, 20), (0, 0.15), (0, 0.001, -0.001), (0, 0.1, 0.2)
    ):
        F, b2 = hinge_F(sim_b(t, k, a, b1, c, delta, reps_b, rng), basis)
        res["B"]["cells"].append(
            {
                "a": a,
                "b1": b1,
                "c": c,
                "delta": delta,
                "pass_rate": float(b_pass(F, b2).mean()),
            }
        )
        print("B", res["B"]["cells"][-1], flush=True)
    cells = res["B"]["cells"]
    res["B"]["power_min_delta0.2"] = min(
        x["pass_rate"] for x in cells if x["delta"] == 0.2
    )
    res["B"]["fpr_max_delta0"] = max(x["pass_rate"] for x in cells if x["delta"] == 0)
    res["B"]["open_real"] = (
        res["B"]["power_min_delta0.2"] >= 0.8 and res["B"]["fpr_max_delta0"] <= 0.10
    )
    diff, p1, p2, dec = zip(*prim)
    X = d_design(diff, p1, p2, dec)
    res["D"]["n_different_rank"] = int(sum(diff))
    for kappa, gamma in itertools.product((0, 0.5), (0, np.log(0.7), np.log(0.5))):
        r = sim_d(X, gamma, kappa, reps_d, rng)
        res["D"]["cells"].append({"kappa": kappa, "gamma": gamma, "rate": r})
        print("D", res["D"]["cells"][-1], flush=True)
    dc = res["D"]["cells"]
    res["D"]["power_min_gamma_log0.7"] = min(
        x["rate"] for x in dc if abs(x["gamma"] - np.log(0.7)) < 1e-9
    )
    res["D"]["fpr_max"] = max(x["rate"] for x in dc if x["gamma"] == 0)
    res["D"]["open_real"] = (
        res["D"]["power_min_gamma_log0.7"] >= 0.8 and res["D"]["fpr_max"] <= 0.07
    )
    with open(out_json, "w") as f:
        json.dump(res, f, indent=1, default=lambda o: o.item())  # numpy scalars
    print(
        json.dumps(
            {s: {k: v for k, v in res[s].items() if k != "cells"} for s in res},
            default=lambda o: o.item(),
        )
    )


if __name__ == "__main__":
    signal.signal(
        signal.SIGALRM,
        lambda *_: (sys.stderr.write("aborting: walltime guard\n"), sys.exit(2)),
    )
    signal.alarm(5400)
    main(*sys.argv[1:4])
