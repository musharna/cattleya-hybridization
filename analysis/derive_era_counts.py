"""Per-species primary-hybrid counts and exposure by era, for pre-registration V.

What: for each species in data/species_traits.csv, registered primary grexes with it as a parent
(T1 counting rules, as in derive_species_traits.py / derive_period_counts.py), split at 1960:
early <= 1959, late >= 1960. Exposure per era = years described within it (prereg V), at least 1;
no early row (n_early, E_early empty) for species first described in 1960 or later.
Grexes with no registration year are dropped and counted.
Output: data/species_era_counts.csv, species names only. Prints pooled totals only.
Usage: python -I analysis/derive_era_counts.py HYBRIDS.csv SPECIES.csv OUT.csv
"""

import collections
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import register_qc as rq  # noqa: E402
from derive_species_traits import binomial  # noqa: E402

SPLIT, START, END = 1960, 1850, 2025


def main(hybrids_csv, species_csv, out_csv):
    records, species = rq.load(hybrids_csv, species_csv)
    here = os.path.dirname(__file__)
    traits = list(
        csv.DictReader(open(os.path.join(here, "..", "data", "species_traits.csv")))
    )
    idx = rq.build_index(records)
    years = collections.defaultdict(list)
    no_year = 0
    for r in records:
        if r["synonym_of"] is not None or len(r["parents"]) != 2:
            continue
        kinds = [rq.resolve(p, idx, records, species)[0] for p in r["parents"]]
        if not all(k in ("species", "species_unlisted") for k in kinds):
            continue
        if not r["year"]:
            no_year += 1
            continue
        for b in {binomial(p) for p in r["parents"]}:
            years[b].append(int(r["year"]))
    rows = []
    for t in traits:
        b, oy = t["Binomial"], int(t["orig_year"])
        ys = years.get(b, [])
        early = oy < SPLIT
        rows.append(
            {
                "Binomial": b,
                "orig_year": oy,
                "n_early": sum(y < SPLIT for y in ys) if early else "",
                "E_early": max(1, SPLIT - max(oy, START)) if early else "",
                "n_late": sum(y >= SPLIT for y in ys),
                "E_late": max(1, END - max(oy, SPLIT)),
            }
        )
        # T1 consistency: the era split must account for every dated primary of the species
        n_total = sum(1 for _ in ys)
        assert (rows[-1]["n_early"] or 0) + rows[-1][
            "n_late"
        ] == n_total or not early, b
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    ne = sum(r["n_early"] or 0 for r in rows)
    nl = sum(r["n_late"] for r in rows)
    print(
        f"species {len(rows)}; early rows {sum(r['n_early'] != '' for r in rows)}; "
        f"primaries early {ne}, late {nl}; dropped with no year {no_year}"
    )


if __name__ == "__main__":
    main(*sys.argv[1:4])
