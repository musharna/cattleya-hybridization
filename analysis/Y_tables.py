"""Derived tables for the Y figures (pre-registration Y + amendments 1-2).

What: data/Y_split.csv (engine, axis, decade, quantity S|T, estimate, 95% parent-unit CI, n, label) from
data/Y_appearance.json; data/Y_decade_pcs.csv (engine, decade, PC1, PC2, n) = first two principal components of the
decade means of all analysed grex photos (descriptive, prereg figure 1; no thumbnails, photos are not redistributed).
Usage: python -I analysis/Y_tables.py HYBRIDS.csv SPECIES.csv PHOTO_DIR
"""

import csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import Y_appearance as Y  # noqa: E402

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def main(hybrids_csv, species_csv, photo_dir):
    res = json.load(open(os.path.join(ROOT, "data", "Y_appearance.json")))
    with open(os.path.join(ROOT, "data", "Y_split.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "engine",
                "axis",
                "decade",
                "quantity",
                "estimate",
                "lo",
                "hi",
                "n",
                "label",
            ]
        )
        for e in ("v8", "dino"):
            for ax in ("trajectory", "minicatt"):
                for D, d in res[e][ax].items():
                    for q in ("S", "T"):
                        lab = d.get("label", "") if ax == "trajectory" else ""
                        w.writerow(
                            [
                                e,
                                ax,
                                D,
                                q,
                                d[q],
                                d[q + "_ci"][0],
                                d[q + "_ci"][1],
                                d["n"],
                                lab,
                            ]
                        )
    records, ped, pids, sp_year, axis_pids = Y.load_all(
        hybrids_csv, species_csv, photo_dir
    )
    with open(os.path.join(ROOT, "data", "Y_decade_pcs.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["engine", "decade", "PC1", "PC2", "n"])
        for e in ("v8", "dino"):
            emb = Y.load_emb(photo_dir, e)
            _, _, gz, _ = Y.build_units(records, ped, pids, sp_year, emb, axis_pids)
            Ds = [
                D
                for D in Y.DECADES
                if any(D <= records[i]["year"] < D + 10 for i in gz)
            ]
            X = np.stack(
                [
                    np.mean([gz[i] for i in gz if D <= records[i]["year"] < D + 10], 0)
                    for D in Ds
                ]
            )
            Xc = X - X.mean(0)
            U, s, Vt = np.linalg.svd(Xc, full_matrices=False)
            pcs = Xc @ Vt[:2].T
            if pcs[-1, 0] < pcs[0, 0]:  # display only: 2020s to the right
                pcs[:, 0] *= -1
            for D, (a, b) in zip(Ds, pcs):
                w.writerow(
                    [e, D, a, b, sum(D <= records[i]["year"] < D + 10 for i in gz)]
                )
            print(
                e,
                "PC1/PC2 share of decade-mean variance: %.2f / %.2f"
                % tuple((s[:2] ** 2) / (s**2).sum()),
            )


if __name__ == "__main__":
    main(*sys.argv[1:4])
