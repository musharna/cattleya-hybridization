"""Per-species primary-hybrid counts before and after the 2008 genus moves (pre-registration T).

What: for each v1 species, registered primary grexes with it as a parent, split by registration
year at the species' cutoff (WCVP v15 year of its Cattleya combination if it was moved in 2006-2014,
else 2008), at 2014 (sensitivity), and at 1985 within <= 2008 (placebo). Counting rules as in
derive_species_traits.py (T1). Grexes with no registration year are dropped and counted.
Output: data/species_period_counts.csv, species names only. Prints pooled totals only, so the
moved-vs-rest split is not seen before the planted-truth check passes.
Usage: python analysis/derive_period_counts.py HYBRIDS.csv SPECIES.csv WCVP_NAMES.csv OUT.csv
"""

import collections
import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(__file__))
import register_qc as rq  # noqa: E402
from derive_species_traits import binomial  # noqa: E402

MOVED_FROM, MOVED_TO = 2006, 2014  # combination years that count as the modern moves


def combination_years(wcvp_csv, names):
    out = {}
    with open(wcvp_csv, encoding="utf-8", errors="replace") as f:
        for r in csv.DictReader(f, delimiter="|"):
            if (
                r["genus"] == "Cattleya"
                and r["taxon_rank"] == "Species"
                and r["taxon_status"] == "Accepted"
                and r["taxon_name"] in names
            ):
                m = re.search(r"\((\d{4})\)", r["first_published"] or "")
                out[r["taxon_name"]] = int(m.group(1)) if m else None
    return out


def main(hybrids_csv, species_csv, wcvp_csv, out_csv):
    records, species = rq.load(hybrids_csv, species_csv)
    here = os.path.dirname(__file__)
    traits = list(
        csv.DictReader(open(os.path.join(here, "..", "data", "species_traits.csv")))
    )
    names = {t["Binomial"] for t in traits}
    comb = combination_years(wcvp_csv, names)
    idx = rq.build_index(records)
    years = collections.defaultdict(
        list
    )  # species -> registration years of its primaries
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
        b = t["Binomial"]
        c = comb.get(b)
        moved = int(
            t["transferred"] == "1" and c is not None and MOVED_FROM <= c <= MOVED_TO
        )
        cut = c if moved else 2008
        ys = years.get(b, [])
        rows.append(
            {
                "Binomial": b,
                "moved": moved,
                "cutoff": cut,
                "descr_year": t["descr_year"],
                "n_pre": sum(y <= cut for y in ys),
                "n_post": sum(y > cut for y in ys),
                "n_pre14": sum(y <= 2014 for y in ys),
                "n_post14": sum(y > 2014 for y in ys),
                "n_p1985": sum(y <= 1985 for y in ys),
                "n_p2008": sum(1985 < y <= 2008 for y in ys),
            }
        )
    total = sum(len(v) for v in years.values())
    assert sum(r["n_pre"] + r["n_post"] for r in rows) <= total
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(
        f"species {len(rows)}; moved {sum(r['moved'] for r in rows)}; primaries without a year {no_year}"
    )
    print(
        f"pooled species-primary counts: pre {sum(r['n_pre'] for r in rows)}, "
        f"post {sum(r['n_post'] for r in rows)}"
    )


if __name__ == "__main__":
    main(*sys.argv[1:5])
