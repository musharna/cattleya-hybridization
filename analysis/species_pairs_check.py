"""Second route for prereg W: refit the real-data distance coefficient with a numpy IRLS.

What: Poisson GLM (numpy IRLS; statsmodels is broken against the installed scipy) on data/species_pairs_real_design.csv (kept pairs, written by species_pairs_real.R):
n ~ d + species effects (two-member indicators) + offset(logE). Shares only the design table with R
(deepest shared layer: the pair counts, cophenetic distances and exposures written there); the fit is
independent. Writes data/species_pairs_check_py.json; the readout requires |beta_py - beta_R| < 0.01.
Usage: python -I analysis/species_pairs_check.py
"""

import json
import os
import signal
import sys

import numpy as np
import pandas as pd

signal.signal(
    signal.SIGALRM,
    lambda *_: (sys.stderr.write("aborting: walltime guard\n"), sys.exit(2)),
)
signal.alarm(600)

here = os.path.join(os.path.dirname(__file__), "..", "data")
sfx = sys.argv[1] if len(sys.argv) > 1 else ""  # "_DRY" for the plumbing test
df = pd.read_csv(os.path.join(here, f"species_pairs_real_design{sfx}.csv"))
sp = sorted(set(df.sp1) | set(df.sp2))
ix = {s: i for i, s in enumerate(sp)}
X = np.zeros((len(df), len(sp)))
X[np.arange(len(df)), df.sp1.map(ix)] = 1
X[np.arange(len(df)), df.sp2.map(ix)] = 1
design = np.column_stack(
    [df.d.values, X]
)  # no separate intercept: the indicator sums span it
# IRLS for a Poisson log-link GLM with offset; start from the species-only rates
y, off = df.n.values.astype(float), df.logE.values
eta = np.log((y + 0.5) / np.exp(off)) + off
for it in range(200):
    mu = np.exp(eta)
    z = eta - off + (y - mu) / mu
    w = np.sqrt(mu)
    params = np.linalg.lstsq(design * w[:, None], z * w, rcond=None)[0]
    eta_new = design @ params + off
    if np.max(np.abs(eta_new - eta)) < 1e-12:
        break
    eta = eta_new
converged = it < 199
out = {
    "beta": float(params[0]),
    "rr": float(np.exp(params[0])),
    "n_pairs": len(df),
    "n_species": len(sp),
    "converged": bool(converged), "iterations": it,
}
assert out["converged"], out
json.dump(out, open(os.path.join(here, f"species_pairs_check_py{sfx}.json"), "w"), indent=1)
print(out)
