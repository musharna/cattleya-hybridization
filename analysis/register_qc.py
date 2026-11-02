"""Register QC for the orchidroots Cattleya scrape (pre-registration: analysis/prereg/G_register_qc.md).

What: parses grex records, resolves each parent reference by exact name including its "(YYYY)"
homonym tag, and flags date-order violations, cycles, self-parents, ambiguous and unresolved
references, missing years/parentage and synonym chains.
When: before any pedigree analysis (A, B, D) and as the test-case report for the Research Data List.
Output: flags per record; a flag is a curation candidate, not a verdict. A wrong-but-plausible
parent (no date or cycle violation) is invisible to these checks.
Test: tests/test_register_qc.py; `plant` re-runs the checks on a copy with known errors.

Usage: python analysis/register_qc.py HYBRIDS.csv SPECIES.csv OUT_DIR SUMMARY.csv
"""

import collections
import csv
import json
import os
import random
import re
import signal
import sys

SEP = re.compile(r"\s+×\s+")
TAG = re.compile(r"\s*\(\d{4}\)$")


def norm(s):
    return " ".join(s.split())


def parse_record(row):
    name = norm(row["Grex Name"])
    year = int(row["Year"]) if row.get("Year", "").strip().isdigit() else None
    par = norm(row.get("Parentage") or "")
    syn = par[len("Synonym of ") :] if par.startswith("Synonym of ") else None
    parents = (
        [] if syn is not None or "×" not in par else [norm(p) for p in SEP.split(par)]
    )
    return {
        "name": name,
        "year": year,
        "parentage": par,
        "synonym_of": syn,
        "parents": parents,
    }


def build_index(records):
    idx = collections.defaultdict(list)
    for i, r in enumerate(records):
        idx[r["name"]].append(i)
    return idx


def _is_species(ref):
    w = ref.split()
    return len(w) > 1 and w[1][:1].islower()


def resolve(ref, idx, records, species):
    ref = norm(ref)
    hits = idx.get(ref, [])
    live = [i for i in hits if records[i]["synonym_of"] is None]
    if len(live) == 1:
        return "grex", live[0]
    if len(live) > 1:
        return "ambiguous", tuple(live)
    if hits:  # only synonym records carry this name
        return "grex", hits[0]
    if _is_species(ref):
        # amendment 1 (prereg G): a species-form name absent from the species list is a gap in
        # our reference list, not a register error
        return ("species", ref) if ref in species else ("species_unlisted", ref)
    return "unresolved", ref


def _edges(records, idx, species):
    """child -> list of (ref, kind, target)."""
    out = []
    for r in records:
        out.append([(p, *resolve(p, idx, records, species)) for p in r["parents"]])
    return out


def _cycle_nodes(n, parents_of):
    """Records on a directed cycle (incl. self-loops), iterative Tarjan over parent links."""
    index, low, on, stack, res, counter = {}, {}, set(), [], set(), [0]
    for root in range(n):
        if root in index:
            continue
        work = [(root, iter(parents_of[root]))]
        index[root] = low[root] = counter[0]
        counter[0] += 1
        stack.append(root)
        on.add(root)
        while work:
            v, it = work[-1]
            nxt = next(it, None)
            if nxt is not None:
                if nxt not in index:
                    index[nxt] = low[nxt] = counter[0]
                    counter[0] += 1
                    stack.append(nxt)
                    on.add(nxt)
                    work.append((nxt, iter(parents_of[nxt])))
                elif nxt in on:
                    low[v] = min(low[v], index[nxt])
                continue
            work.pop()
            if work:
                low[work[-1][0]] = min(low[work[-1][0]], low[v])
            if low[v] == index[v]:
                comp = []
                while True:
                    w = stack.pop()
                    on.discard(w)
                    comp.append(w)
                    if w == v:
                        break
                if len(comp) > 1 or v in parents_of[v]:
                    res.update(comp)
    return res


def qc(records, species):
    idx = build_index(records)
    edges = _edges(records, idx, species)
    flags = []
    parents_of = [[t for _, k, t in e if k == "grex"] for e in edges]
    for i, (r, e) in enumerate(zip(records, edges)):
        if r["synonym_of"] is not None:
            k, t = resolve(r["synonym_of"], idx, records, species)
            if k != "grex" or records[t]["synonym_of"] is not None:
                flags.append((i, "synonym_chain", r["synonym_of"]))
            continue
        if r["year"] is None:
            flags.append((i, "missing_year", ""))
        if not r["parents"]:
            flags.append((i, "missing_parentage", r["parentage"]))
        for ref, k, t in e:
            if k == "ambiguous":
                flags.append((i, "ambiguous_ref", ref))
            elif k == "unresolved":
                flags.append((i, "unresolved_ref", ref))
            elif k == "species_unlisted":
                flags.append((i, "species_unlisted", ref))
            elif k == "grex":
                if t == i:
                    flags.append((i, "self_parent", ref))
                py = records[t]["year"]
                if py is not None and r["year"] is not None and py > r["year"]:
                    flags.append((i, "date_order", f"{ref} {py} > {r['year']}"))
    for i in sorted(_cycle_nodes(len(records), parents_of)):
        flags.append((i, "cycle", ""))
    return flags


def plant(
    records, rng, n_date=50, n_cycle=20, n_typo=100, n_untag=20, species=frozenset()
):
    """Copy with known errors; returns (copy, {class: set(record idx)})."""
    rs = [dict(r, parents=list(r["parents"])) for r in records]
    idx = build_index(rs)
    edges = _edges(rs, idx, species)
    clean = {i for i, _, _ in qc(rs, species)}
    gx = [
        (i, j, t)
        for i, e in enumerate(edges)
        for j, (_, k, t) in enumerate(e)
        if k == "grex" and i not in clean and t not in clean and t != i
    ]
    truth = {c: set() for c in ("date_order", "cycle", "unresolved_ref", "untag")}
    used = set()
    for i, j, t in rng.sample(gx, len(gx)):
        if len(truth["date_order"]) >= n_date:
            break
        if i in used or rs[t]["year"] is None:
            continue
        rs[i]["year"] = rs[t]["year"] - 1
        used.add(i)
        truth["date_order"].add(i)
    for i, j, t in rng.sample(gx, len(gx)):
        if len(truth["cycle"]) >= n_cycle:
            break
        if i in used or t in used or not rs[t]["parents"]:
            continue
        rs[t]["parents"][0] = rs[i]["name"]
        used.update((i, t))
        truth["cycle"].add(t)
    for i, j, t in rng.sample(gx, len(gx)):
        if len(truth["unresolved_ref"]) >= n_typo:
            break
        if i in used:
            continue
        ref = rs[i]["parents"][j]
        w = ref.split()
        if len(w[-1]) < 3:
            continue
        pos = rng.randrange(1, len(w[-1]) - 1)
        w[-1] = w[-1][:pos] + ("q" if w[-1][pos] != "q" else "z") + w[-1][pos + 1 :]
        typo = " ".join(w)
        if typo in idx:
            continue
        rs[i]["parents"][j] = typo
        used.add(i)
        truth["unresolved_ref"].add(i)
    for i, j, t in rng.sample(gx, len(gx)):
        if len(truth["untag"]) >= n_untag:
            break
        ref = rs[i]["parents"][j]
        if i in used or not TAG.search(ref) or len(idx.get(TAG.sub("", ref), [])) == 0:
            continue
        rs[i]["parents"][j] = TAG.sub("", ref)
        used.add(i)
        truth["untag"].add(i)
    return rs, truth


def load(hybrids_csv, species_csv):
    with open(hybrids_csv, encoding="utf-8", errors="replace") as f:
        records = [parse_record(r) for r in csv.DictReader(f)]
    with open(species_csv, encoding="utf-8", errors="replace") as f:
        species = {norm(r["Binomial"]) for r in csv.DictReader(f)}
    return records, species


def main(hybrids_csv, species_csv, out_dir, summary_csv):
    records, species = load(hybrids_csv, species_csv)
    flags = qc(records, species)
    assert flags == qc(records, species), "QC is not deterministic"
    planted, truth = plant(records, random.Random(0), species=species)
    pflags = qc(planted, species)
    real = {(i, c) for i, c, _ in flags}
    new = {(i, c) for i, c, _ in pflags} - real
    by = collections.defaultdict(set)
    for i, c in new:
        by[c].add(i)
    idx = build_index(planted)
    nbr = set().union(*truth.values())
    for i in list(
        nbr
    ):  # direct parents/children of planted records are mechanically implied
        for p in planted[i]["parents"]:
            for k in idx.get(p, []):
                nbr.add(k)
    for k, r in enumerate(planted):
        if any(
            t in set().union(*truth.values())
            for p in r["parents"]
            for t in idx.get(p, [])
        ):
            nbr.add(k)
    recall = {
        c: len(truth[c] & by[c]) / len(truth[c])
        for c in ("date_order", "cycle", "unresolved_ref")
    }
    untag_hit = len(truth["untag"] & (by["ambiguous_ref"] | by["date_order"])) / max(
        len(truth["untag"]), 1
    )
    unexplained = sorted({i for i, _ in new if i not in nbr})
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "qc_flags.tsv"), "w") as f:
        f.write("grex\tyear\tclass\tdetail\tparentage\n")
        for i, c, d in flags:
            r = records[i]
            f.write(f"{r['name']}\t{r['year']}\t{c}\t{d}\t{r['parentage']}\n")
    counts = collections.Counter(c for _, c, _ in flags)
    with open(summary_csv, "w") as f:
        f.write("class,records_flagged\n")
        for c in (
            "date_order",
            "cycle",
            "self_parent",
            "ambiguous_ref",
            "unresolved_ref",
            "species_unlisted",
            "missing_year",
            "missing_parentage",
            "synonym_chain",
        ):
            f.write(f"{c},{len({i for i, cc, _ in flags if cc == c})}\n")
    control = {
        "n_records": len(records),
        "n_synonym": sum(r["synonym_of"] is not None for r in records),
        "flags_by_class": counts,
        "planted": {c: len(v) for c, v in truth.items()},
        "recall": recall,
        "untag_caught_as_ambiguous_or_date": untag_hit,
        "unexplained_new_flags": len(unexplained),
        "PASS": all(v >= 0.95 for v in recall.values()) and not unexplained,
    }
    with open(os.path.join(out_dir, "qc_control.json"), "w") as f:
        json.dump(control, f, indent=1, default=dict)
    print(json.dumps(control, default=dict))


if __name__ == "__main__":
    signal.signal(
        signal.SIGALRM,
        lambda *_: (sys.stderr.write("aborting: walltime guard\n"), sys.exit(2)),
    )
    signal.alarm(900)
    main(*sys.argv[1:5])
