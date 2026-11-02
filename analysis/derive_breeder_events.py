"""Breeder-years, breeder network and adoption years for pre-registration X.

What: breeders = normalized Originator (or Registrant with --registrant). For each breeder: grexes per
registration year; adoption year = first year a grex of theirs has a mini-catt founder species
(C. coccinea, walkeriana, aclandiae, intermedia) as a direct parent (also walkeriana alone); network edges
b -> c with the first year b registered a grex using one of c's grexes as a parent.
Output (local only, names people): OUT_DIR/breeder_years.csv, OUT_DIR/edges.csv, OUT_DIR/adoption.csv.
Prints aggregates only.
Usage: python -I analysis/derive_breeder_events.py HYBRIDS.csv SPECIES.csv OUT_DIR [--registrant]
"""

import collections
import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(__file__))
import register_qc as rq  # noqa: E402
from derive_species_traits import binomial  # noqa: E402

M = {
    "Cattleya coccinea",
    "Cattleya walkeriana",
    "Cattleya aclandiae",
    "Cattleya intermedia",
}


def norm(s):
    s = re.sub(r"[^a-z]", "", (s or "").lower())
    return "" if s in ("", "ou") else s


def main(hybrids_csv, species_csv, out_dir, field="Originator"):
    records, species = rq.load(hybrids_csv, species_csv)
    who = collections.defaultdict(set)
    for r in csv.DictReader(open(hybrids_csv, encoding="utf-8")):
        who[r["Grex Name"]].add(norm(r[field]))
    dup = sum(1 for v in who.values() if len(v) > 1)
    breeder = [
        next(iter(who[r["name"]])) if len(who.get(r["name"], ())) == 1 else ""
        for r in records
    ]
    idx = rq.build_index(records)
    n_year = collections.Counter()
    adopt, adopt_w = {}, {}
    adopt_sp = collections.defaultdict(dict)  # per founder species, for first-use curves
    n_mini = collections.Counter()  # grexes from 1950 with a founder species as a direct parent
    edge = {}
    for i, r in enumerate(records):
        b = breeder[i]
        if r["synonym_of"] is not None or not r["year"] or not b:
            continue
        y = int(r["year"])
        n_year[(b, y)] += 1
        res = [rq.resolve(p, idx, records, species) for p in r["parents"]]
        sp = {binomial(t) for k, t in res if k in ("species", "species_unlisted")}
        if sp & M:
            adopt[b] = min(adopt.get(b, 9999), y)
            if y >= 1950:
                n_mini[b] += 1
        for m in sp & M:
            adopt_sp[m][b] = min(adopt_sp[m].get(b, 9999), y)
        if "Cattleya walkeriana" in sp:
            adopt_w[b] = min(adopt_w.get(b, 9999), y)
        for k, t in res:
            if k == "grex":
                c = breeder[t]
                if c and c != b:
                    edge[(b, c)] = min(edge.get((b, c), 9999), y)
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "breeder_years.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["breeder", "year", "n"])
        for (b, y), n in sorted(n_year.items()):
            w.writerow([b, y, n])
    with open(os.path.join(out_dir, "edges.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["breeder", "contact", "first_year"])
        for (b, c), y in sorted(edge.items()):
            w.writerow([b, c, y])
    # contacts seen only on undated or synonym grexes are breeders too (never at risk, never adopt)
    breeders = sorted({b for b, _ in n_year} | {c for _, c in edge})
    with open(os.path.join(out_dir, "adoption.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["breeder", "adopt_year", "adopt_year_walkeriana", "n_mini_from1950"] + [f"adopt_{m.split()[1]}" for m in sorted(M)])
        for b in breeders:
            w.writerow([b, adopt.get(b, ""), adopt_w.get(b, ""), n_mini[b]] + [adopt_sp[m].get(b, "") for m in sorted(M)])
    print(
        f"field {field}: breeders {len(breeders)}; grex names with conflicting breeder {dup}; "
        f"edges {len(edge)}; adopters {len(adopt)} ({sum(y >= 1950 for y in adopt.values())} from 1950); "
        f"walkeriana adopters {len(adopt_w)}"
    )


if __name__ == "__main__":
    a = [x for x in sys.argv[1:] if not x.startswith("--")]
    main(*a[:3], field="Registrant" if "--registrant" in sys.argv else "Originator")
