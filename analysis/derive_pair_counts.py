"""Registered primary-hybrid counts per species pair, for pre-registration W.

What: for every unordered pair of species in data/species_traits.csv, the number of registered
primary grexes with those two species as parents (T1 counting rules, as in derive_era_counts.py:
synonyms dropped, both parents resolve to species; intraspecific crosses excluded), plus the
split at 1960 (early <= 1959, late >= 1960; undated grexes count in n only).
Output: data/species_pair_counts.csv (sp1, sp2, n, n_early, n_late), species names only; pairs with
n = 0 are not written (the model script fills them in). Prints totals only.
Usage: python -I analysis/derive_pair_counts.py HYBRIDS.csv SPECIES.csv OUT.csv
"""

import collections
import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(__file__))
import register_qc as rq  # noqa: E402
from derive_species_traits import binomial  # noqa: E402

SPLIT = 1960


def main(hybrids_csv, species_csv, out_csv):
    records, species = rq.load(hybrids_csv, species_csv)
    here = os.path.dirname(__file__)
    traits = {
        t["Binomial"]: int(t["primary_T1"])
        for t in csv.DictReader(
            open(os.path.join(here, "..", "data", "species_traits.csv"))
        )
    }
    orig = {
        t["Binomial"]: int(t["orig_year"])
        for t in csv.DictReader(open(os.path.join(here, "..", "data", "species_traits.csv")))
    }
    tips = set(
        s.replace("_", " ")
        for s in re.findall(
            r"Cattleya_[a-z]+",
            open(os.path.join(here, "..", "data", "molecular_combined.nwk")).read(),
        )
    )
    assert len(tips) == 91 and tips <= set(traits), (len(tips), tips - set(traits))
    idx = rq.build_index(records)
    n, early, late = collections.Counter(), collections.Counter(), collections.Counter()
    self_cross = predated = 0
    for r in records:
        if r["synonym_of"] is not None or len(r["parents"]) != 2:
            continue
        kinds = [rq.resolve(p, idx, records, species)[0] for p in r["parents"]]
        if not all(k in ("species", "species_unlisted") for k in kinds):
            continue
        a, b = (binomial(p) for p in r["parents"])
        if a == b:
            self_cross += 1
            continue
        if a not in traits or b not in traits:
            continue
        key = tuple(sorted((a, b)))
        n[key] += 1
        if r["year"]:
            # prereg W Amendment 2: a grex dated before a parent species was first described has a wrong date
            # (e.g. a reused grex name); it stays in n but is left out of the era split
            if int(r["year"]) < max(orig[a], orig[b]):
                predated += 1
                continue
            (early if int(r["year"]) < SPLIT else late)[key] += 1
    # T1 consistency: no species can have more pair-counted primaries than its T1 total
    per_sp = collections.Counter()
    for (a, b), c in n.items():
        per_sp[a] += c
        per_sp[b] += c
    over = {s: (per_sp[s], traits[s]) for s in per_sp if per_sp[s] > traits[s]}
    assert not over, over
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sp1", "sp2", "n", "n_early", "n_late"])
        for a, b in sorted(n):
            w.writerow([a, b, n[(a, b)], early[(a, b)], late[(a, b)]])
    on_tree = sum(c for (a, b), c in n.items() if a in tips and b in tips)
    print(
        f"pairs with >=1 primary: {len(n)} ({sum(n.values())} primaries); "
        f"both on the 91-species tree: {sum(1 for a, b in n if a in tips and b in tips)} pairs, {on_tree} primaries; "
        f"intraspecific excluded {self_cross}; dated before a parent's description (era split only) {predated}; "
        f"pair-counted share of T1 per species min "
        f"{min(per_sp[s] / traits[s] for s in per_sp if traits[s]):.2f}"
    )


if __name__ == "__main__":
    main(*sys.argv[1:4])
