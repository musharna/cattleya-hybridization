"""Per-species traits for the corrected signal test (pre-registration C).

What: T1 = registered primary grexes with the species as a parent, T2 = grexes of any type with it
as an immediate parent (parents resolved by register_qc; var./subsp./f. forms map to the binomial),
description year, WCVP basionym genus and a `transferred` flag (basionym genus is not Cattleya).
Output: data/species_traits.csv, species names only (no grex records), one row per v1 species.
Usage: python analysis/derive_species_traits.py HYBRIDS.csv SPECIES.csv WCVP_NAMES.csv OUT.csv
"""

import collections
import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(__file__))
import register_qc as rq  # noqa: E402


def pub_year(r):
    m = re.search(r"\((\d{4})\)", r["first_published"] or "")
    return int(m.group(1)) if m else None


def binomial(ref):
    return " ".join(ref.split()[:2])


def main(hybrids_csv, species_csv, wcvp_csv, out_csv):
    records, species = rq.load(hybrids_csv, species_csv)
    with open(species_csv, encoding="utf-8", errors="replace") as f:
        year = {rq.norm(r["Binomial"]): r["Year"] for r in csv.DictReader(f)}
    v1 = list(
        csv.DictReader(
            open(
                os.path.join(
                    os.path.dirname(__file__), "..", "data", "species_hybrid_counts.csv"
                )
            )
        )
    )
    idx = rq.build_index(records)
    t1, t2 = collections.Counter(), collections.Counter()
    n_primary = 0
    for r in records:
        if r["synonym_of"] is not None or len(r["parents"]) != 2:
            continue
        kinds = [rq.resolve(p, idx, records, species)[0] for p in r["parents"]]
        sp = {
            binomial(p)
            for p, k in zip(r["parents"], kinds)
            if k in ("species", "species_unlisted")
        }
        for b in sp:  # T2: each grex once per distinct parent species
            t2[b] += 1
        if all(k in ("species", "species_unlisted") for k in kinds):
            n_primary += 1
            for b in sp:
                t1[b] += 1
    wc = {}
    by_id = {}
    with open(wcvp_csv, encoding="utf-8", errors="replace") as f:
        for r in csv.DictReader(f, delimiter="|"):
            by_id[r["plant_name_id"]] = r
            if (
                r["genus"] == "Cattleya"
                and r["taxon_rank"] == "Species"
                and not r["species_hybrid"]
            ):
                wc.setdefault(r["taxon_name"], []).append(r)
    rows, missing = [], 0
    for v in v1:
        b = v["Binomial"]
        acc = [r for r in wc.get(b, []) if r["taxon_status"] == "Accepted"]
        bg, oy = "", None
        if acc:
            bid = acc[0]["basionym_plant_name_id"]
            bg = by_id[bid]["genus"] if bid in by_id else "Cattleya"
            oy = pub_year(by_id[bid] if bid in by_id else acc[0])
        else:
            missing += 1
        rows.append(
            {
                "Binomial": b,
                "LastSubgenus": v["LastSubgenus"],
                "Descendants": v["Descendants"],
                "primary_T1": t1[b],
                "direct_T2": t2[b],
                "descr_year": year.get(b, ""),
                # prereg C Amendment 5: the species list's Year is the year of the current name
                # (2006-2014 for the 2008 moves); exposure uses the original description instead
                "orig_year": oy or year.get(b, ""),
                "basionym_genus": bg,
                "transferred": int(bool(bg) and bg != "Cattleya"),
            }
        )
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(
        f"species {len(rows)}; primary grexes {n_primary}; not accepted in WCVP {missing}; "
        f"transferred {sum(r['transferred'] for r in rows)}; zero-T1 {sum(r['primary_T1'] == 0 for r in rows)}; "
        f"missing year {sum(not r['descr_year'] for r in rows)}; no WCVP original year {sum(r['orig_year'] == year.get(r['Binomial'], '') and r['transferred'] for r in rows)} transferred"
    )


if __name__ == "__main__":
    main(*sys.argv[1:5])
