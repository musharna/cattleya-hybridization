"""Species-level predictors from WCVP for prereg U (what predicts breeders' use of a species).

What: for each species in data/species_traits.csv, WCVP's accepted Cattleya species record gives
`brazil` (native to any Brazilian level-3 area), `wet` (climate "wet tropical"), `litho` (life form
mentions lithophyte) and `n_areas` (native, non-extinct level-3 areas). Hybrid counts are not read.
Usage: python -I analysis/derive_species_covariates.py WCVP_NAMES.csv WCVP_DISTRIBUTION.csv OUT.csv
Fails loudly if a species has no accepted WCVP record or no native area.
"""

import csv
import sys


def main(names_csv, dist_csv, out_csv):
    species = [r["Binomial"] for r in csv.DictReader(open("data/species_traits.csv"))]
    want = set(species)
    rec = {}
    with open(names_csv, encoding="utf-8", errors="replace") as f:
        for r in csv.DictReader(f, delimiter="|"):
            n = "Cattleya " + r["species"]
            if (
                r["genus"] == "Cattleya"
                and r["taxon_rank"] == "Species"
                and r["taxon_status"] == "Accepted"
                and n in want
            ):
                assert n not in rec, f"two accepted records for {n}"
                rec[n] = r
    missing = want - set(rec)
    areas = {n: set() for n in rec}
    by_id = {r["plant_name_id"]: n for n, r in rec.items()}
    with open(dist_csv, encoding="utf-8", errors="replace") as f:
        for r in csv.DictReader(f, delimiter="|"):
            n = by_id.get(r["plant_name_id"])
            if n and r["introduced"] == "0" and r["extinct"] == "0":
                areas[n].add(r["area_code_l3"])
    no_area = sorted(n for n, a in areas.items() if not a)
    rows = []
    for n in species:
        if n not in rec:
            continue
        r = rec[n]
        assert r["climate_description"] in (
            "wet tropical",
            "seasonally dry tropical",
        ), (n, r["climate_description"])
        rows.append(
            {
                "Binomial": n,
                "brazil": int(any(a.startswith("BZ") for a in areas[n])),
                "wet": int(r["climate_description"] == "wet tropical"),
                "litho": int("lithophyte" in r["lifeform_description"]),
                "n_areas": len(areas[n]),
            }
        )
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(
        f"species {len(species)}; written {len(rows)}; no accepted WCVP record {len(missing)}: {sorted(missing)}; "
        f"no native area {len(no_area)}: {no_area}"
    )


if __name__ == "__main__":
    main(*sys.argv[1:4])
