"""Post-hoc triage of register_qc flags: which are artefacts, which have a benign reading, which remain.

What: re-reads each G flag against the raw list and WCVP and gives it a subclass, e.g. a species-style
  (natural-hybrid) parent whose Year is a publication date, an older grex of the same cross sunk as a
  synonym into a later-registered parent, a parent registered late by the breeder who used it.
When: after register_qc.py, before sending flags to curators. NOT pre-registered (G prereg, post-hoc note).
Output: TSV with grex, year, class, detail, subclass, note (names grexes: keep it out of the repo).
  "probable_register_error" is what is left after every benign reading was tried.
Test: none of its own. The frozen curator spot check tests it: rows it calls artefacts should come back not_error.
Usage: python -I analysis/register_qc_resolve.py HYBRIDS.csv SPECIES.csv CRAWL.csv WCVP_NAMES OUT.tsv
"""

import csv
import collections
import os
import re
import sys
import signal

signal.signal(
    signal.SIGALRM,
    lambda *_: (sys.stderr.write("aborting: walltime guard\n"), sys.exit(2)),
)
signal.alarm(300)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import register_qc as rq

csv.field_size_limit(10**8)
HYB, SPP, CRAWL, WCVP, OUT = sys.argv[1:6]
raw = list(csv.DictReader(open(HYB, encoding="utf-8", errors="replace")))
records, species = rq.load(HYB, SPP)
assert len(raw) == len(records)
idx = rq.build_index(records)
base = collections.defaultdict(list)
for i, r in enumerate(records):
    base[rq.TAG.sub("", r["name"])].append(i)
crawl = {rq.norm(d["grex_name"]) for d in csv.DictReader(open(CRAWL, encoding="utf-8"))}
W = [
    d
    for d in csv.DictReader(open(WCVP, encoding="utf-8"), delimiter="|")
    if d["family"] == "Orchidaceae"
]
wname = collections.defaultdict(list)
for d in W:
    wname[d["taxon_name"]].append(d)
syn_into = collections.defaultdict(list)
for r in records:
    if r["synonym_of"] is not None:
        syn_into[r["synonym_of"]].append((r["name"], r["year"]))
tagbad = [
    (r["name"], r["year"])
    for r in records
    if re.search(r"\((\d{4})\)$", r["name"])
    and r["year"] is not None
    and int(re.search(r"\((\d{4})\)$", r["name"]).group(1)) != r["year"]
]
print(
    "records whose (YYYY) tag differs from Year:",
    len(tagbad),
    tagbad[:15],
    file=sys.stderr,
)
sp_style = lambda n: len(n.split()) > 1 and n.split()[1][:1].islower()
flags = rq.qc(records, species)
out = []  # (grex, year, class, detail, subclass, note)
for i, cls, det in flags:
    r = records[i]
    sub, note = "", ""
    if cls == "date_order":
        ref = det.rsplit(" ", 3)[0]
        k, t = rq.resolve(ref, idx, records, species)
        p = records[t]
        gap = p["year"] - r["year"]
        older = [
            j
            for j in base[rq.TAG.sub("", p["name"])]
            if j != t
            and records[j]["year"] is not None
            and records[j]["year"] <= r["year"]
        ]
        tagm = re.search(r"\((\d{4})\)$", p["name"])
        sunk = [
            (n, y)
            for n, y in syn_into.get(p["name"], [])
            if y is not None and y <= r["year"]
        ]
        if sp_style(p["name"]):
            sub, note = (
                "artefact_parent_natural_hybrid",
                f"parent is a species-style entry; registrant '{raw[t]['Registrant']}', originator '{raw[t]['Originator']}'",
            )
        elif sp_style(r["name"]):
            sub, note = (
                "artefact_child_natural_hybrid",
                f"child is a species-style entry; registrant '{raw[i]['Registrant']}'",
            )
        elif (
            tagm and int(tagm.group(1)) != p["year"] and int(tagm.group(1)) <= r["year"]
        ):
            sub, note = (
                "parent_year_contradicts_own_tag",
                f"parent named with ({tagm.group(1)}) but Year field {p['year']}",
            )
        elif sunk:
            sub, note = (
                "synonym_consolidation",
                "older grex(es) of the same cross sunk into the parent: "
                + "; ".join(f"{n} {y}" for n, y in sunk),
            )
        elif older:
            sub, note = (
                "probable_wrong_homonym",
                "older same-name record(s): "
                + "; ".join(
                    f"{records[j]['name']} {records[j]['year']}" for j in older
                ),
            )
        elif gap <= 2:
            sub, note = "registration_lag_le2y", f"gap {gap} y"
        else:
            sub, note = "probable_register_error", f"gap {gap} y"
    elif cls == "ambiguous_ref":
        live = [j for j in idx[rq.norm(det)] if records[j]["synonym_of"] is None]
        ok = [
            j
            for j in live
            if records[j]["year"] is not None
            and r["year"] is not None
            and records[j]["year"] <= r["year"]
        ]
        yrs = ", ".join(str(records[j]["year"]) for j in live)
        if len(ok) == 1:
            sub, note = (
                "resolved_by_year",
                f"only the {records[ok[0]]['year']} record predates the child (candidates {yrs})",
            )
        elif not ok:
            sub, note = "no_candidate_predates", f"candidates {yrs}"
        else:
            sub, note = (
                "still_ambiguous",
                f"{len(ok)} candidates predate (candidates {yrs})",
            )
    elif cls == "unresolved_ref":
        g = det.split()[0]
        b = rq.TAG.sub("", det)
        if b != det and b in base:
            tag = det[len(b) :].strip(" ()")
            sub, note = (
                "tag_mismatch",
                f"no record '{det}'; same-name records years "
                + ", ".join(str(records[j]["year"]) for j in base[b])
                + f" (tag {tag})",
            )
        elif det in crawl:
            sub, note = (
                "other_genus_in_crawl",
                "parent exists in the OrchidRoots crawl under another genus",
            )
        elif g in ("Calanthe", "Paphiopedilum"):
            sub, note = (
                "probable_register_error",
                f"{g} cannot be a Cattleya-alliance parent",
            )
        else:
            sub, note = (
                "other_genus_out_of_scope",
                f"{g} parent; not in the Cattleya file",
            )
    elif cls == "synonym_chain":
        g = det.split()[0]
        if g != "Cattleya":
            sub, note = (
                "target_other_genus",
                f"synonym target is in {g} (genus transfer); outside the Cattleya file"
                + ("; found in crawl" if det in crawl else ""),
            )
        else:
            hits = idx.get(rq.norm(det), [])
            if not hits:
                sub, note = "target_missing", "no record of that name"
            elif all(records[j]["synonym_of"] is not None for j in hits):
                sub, note = (
                    "synonym_of_synonym",
                    "target is itself a synonym of " + records[hits[0]]["synonym_of"],
                )
            else:
                sub, note = "target_ambiguous", f"{len(hits)} records of that name"
    elif cls == "species_unlisted":
        ws = wname.get(det, [])
        if not ws:
            # try without infraspecific part
            bin_ = " ".join(det.split()[:2])
            ws2 = wname.get(bin_, [])
            sub, note = (
                ("wcvp_binomial_only" if ws2 else "not_in_wcvp"),
                (
                    f"binomial {bin_}: "
                    + "; ".join(f"{d['taxon_status']}" for d in ws2)
                    if ws2
                    else ""
                ),
            )
        else:
            acc = []
            for d in ws:
                a = (
                    d["plant_name_id"]
                    if d["taxon_status"] == "Accepted"
                    else d["accepted_plant_name_id"]
                )
                acc += (
                    [x["taxon_name"] for x in W if x["plant_name_id"] == a] if a else []
                )
            sub, note = (
                "wcvp_" + ws[0]["taxon_status"].lower().replace(" ", "_"),
                "accepted as " + ", ".join(sorted(set(acc))),
            )
        if raw[i]["Parentage"].count(det.split()[0]) == 2 and det.split()[0] not in (
            "Cattleya",
        ):
            sub += "+both_parents_other_genus"
        if det.split()[0] not in (
            "Cattleya",
            "Laelia",
            "Sophronitis",
            "Guarianthe",
            "Brassavola",
            "Encyclia",
            "Rhyncholaelia",
            "Epidendrum",
            "Broughtonia",
            "Prosthechea",
            "Hoffmannseggella",
        ):
            sub += "+alliance_doubt"
    out.append((r["name"], r["year"], cls, det, sub, note))
# late registration: a parent whose own registrant also bred one of its too-early children was in
# use before it was registered; every too-early child of that parent gets the benign reading
breeder = lambda s: re.sub(
    r"[^a-z]", "", s.lower().split("[")[0].split("(")[0].split("/")[0]
)[:6]
rawn = {rq.norm(x["Grex Name"]): x for x in raw}
late = set()
for g, y, cls, det, sub, note in out:
    if sub == "probable_register_error" and cls == "date_order":
        p = det.rsplit(" ", 3)[0]
        if breeder(rawn[g]["Registrant"]) and breeder(rawn[g]["Registrant"]) == breeder(
            rawn[p]["Registrant"]
        ):
            late.add(p)
out = [
    (
        g,
        y,
        cls,
        det,
        "late_registered_parent",
        note + "; parent's registrant bred an earlier child",
    )
    if sub == "probable_register_error"
    and cls == "date_order"
    and det.rsplit(" ", 3)[0] in late
    else (g, y, cls, det, sub, note)
    for g, y, cls, det, sub, note in out
]


# duplicate synonym names (added after the Louis Chaton case): one name, several synonym records.
# Same registrant and year but different target crosses = a conflict (one link is wrong); otherwise
# the records are usually different grexes that shared a name, each sunk into its own cross.
def target_cross(t):
    live = [j for j in idx.get(rq.norm(t), []) if records[j]["synonym_of"] is None]
    return frozenset(records[live[0]]["parents"]) if len(live) == 1 else None


syn_recs = collections.defaultdict(list)
for i, r in enumerate(records):
    if r["synonym_of"] is not None:
        syn_recs[r["name"]].append(i)
for n, ii in sorted(syn_recs.items()):
    if len(ii) < 2:
        continue
    targets = {records[i]["synonym_of"] for i in ii}
    crosses = [target_cross(records[i]["synonym_of"]) for i in ii]
    who = {(breeder(raw[i]["Registrant"]), records[i]["year"]) for i in ii}
    if len(targets) == 1:
        sub = "same_target_repeated"
    elif None in crosses:
        sub = "target_unresolvable"
    elif len(set(crosses)) == 1:
        sub = "targets_same_cross"
    elif len(who) == 1:
        sub = "conflict_same_registrant_year"
    else:
        sub = "homonym_different_grexes"
    note = "; ".join(
        f"{records[i]['year']} {raw[i]['Registrant']} -> {records[i]['synonym_of']}"
        for i in ii
    )
    out.append(
        (
            n,
            records[ii[0]]["year"],
            "duplicate_synonym",
            ", ".join(sorted(targets)),
            sub,
            note,
        )
    )
# duplicate crosses: one species × species cross carried by several live names. Grex names follow
# parentage, so only one should be live, unless one entry is the natural hybrid (a species-style name).
by_cross = collections.defaultdict(list)
for i, r in enumerate(records):
    if (
        r["synonym_of"] is None
        and len(r["parents"]) == 2
        and all(sp_style(p) for p in r["parents"])
    ):
        by_cross[frozenset(r["parents"])].append(i)
for k, ii in sorted(by_cross.items(), key=lambda kv: sorted(kv[0])):
    if len(ii) < 2:
        continue
    nat = [i for i in ii if sp_style(records[i]["name"])]
    gx = [i for i in ii if not sp_style(records[i]["name"])]
    ys = sorted(records[i]["year"] or 9999 for i in gx)
    if len(k) == 1:
        sub = "self"
    elif len({records[i]["name"] for i in gx}) < len(gx):
        sub = "duplicate_record"
    elif nat and len(gx) <= 1:
        sub = "natural_hybrid_and_grex" if len(nat) == 1 else "several_natural_entries"
    elif nat:
        sub = "natural_hybrid_and_several_grexes"
    elif ys[0] < 1960 <= ys[-1]:
        sub = "old_grex_registered_again_from_1960"
    else:
        sub = "grexes_same_era"
    ii = sorted(ii, key=lambda i: records[i]["year"] or 9999)
    note = "; ".join(
        f"{records[i]['name']} {records[i]['year']} {raw[i]['Registrant']}" for i in ii
    )
    out.append(
        (
            records[ii[0]]["name"],
            records[ii[0]]["year"],
            "duplicate_cross",
            " × ".join(sorted(k)),
            sub,
            note,
        )
    )
w = csv.writer(open(OUT, "w", newline=""), delimiter="\t")
w.writerow(["grex", "year", "class", "detail", "subclass", "note"])
w.writerows(out)
c = collections.Counter((o[2], o[4]) for o in out)
for k in sorted(c):
    print(k[0], k[1], c[k], sep="\t")
