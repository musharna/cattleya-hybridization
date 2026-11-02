"""Fetch GenBank ITS / matK / trnL-F for the C species and outgroups (prereg C, Amendment 2).

What: searches nuccore for each marker under Cattleya and its former genera plus Laeliinae outgroups,
  maps each record's organism to a WCVP v15 accepted name, drops hybrids and short or plastome-length
  records, and keeps the longest record per species per marker.
When: once, before building the molecular tree. Sequences are public; the accession list is the record.
Output: ACCESSIONS_CSV (marker, species, accession, length, organism, title) and OUT_DIR/<marker>.fasta
  with headers ">Genus_species". A species missing from a marker is simply absent from that file.
Test: the coverage counts printed at the end must match the 2026-10-25 probe (ITS 92/100) within a few.
Usage: python -I analysis/molecular_fetch.py WCVP_NAMES SPECIES_TRAITS.csv OUT_DIR ACCESSIONS_CSV
"""

import csv
import json
import signal
import sys
import time
import urllib.parse
import urllib.request

E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
UA = {"User-Agent": "orchid-research-probe"}
INGROUP = [
    "Cattleya",
    "Laelia",
    "Sophronitis",
    "Hadrolaelia",
    "Hoffmannseggella",
    "Sophronitella",
    "Dungsia",
    "Microlaelia",
    "Brasilaelia",
    "Chironiella",
    "Schluckebieria",
    "Cattleyella",
]
OUTGROUP = ["Guarianthe", "Brassavola", "Rhyncholaelia"]
MARKERS = {
    "ITS": (
        '(ITS1[All Fields] OR ITS2[All Fields] OR "internal transcribed spacer"[All Fields])',
        400,
    ),
    "matK": ("matK[All Fields]", 500),
    "trnLF": ("(trnL[All Fields] OR trnL-F[All Fields])", 500),
}
# the field search also matches feature text (a trnK/matK record hits "trnL"), so the title must name the marker
# trnL-F: the trnL-trnF spacer or the trnL intron; rpl32-trnL(UAG) records also say "trnL" but are another region
TITLE = {
    "ITS": lambda t: any(k in t for k in ("internal transcribed spacer", "its1", "its2", "5.8s")),
    "matK": lambda t: "matk" in t,
    "trnLF": lambda t: ("trnl-trnf" in t or "trnl-f" in t or ("trnl" in t and "intron" in t)) and "rpl32" not in t,
}
MAX_LEN = 3000  # longer records are plastomes or multi-region assemblies; not this marker alone


def call(path, **p):
    p["tool"] = "orchid-research-probe"
    data = urllib.parse.urlencode(p).encode()
    req = urllib.request.Request(E + path, data=data, headers=UA)
    time.sleep(0.5)  # NCBI: max 3 requests/s without a key
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read().decode()


def accepted_map(wcvp):
    rows = [
        d
        for d in csv.DictReader(open(wcvp, encoding="utf-8"), delimiter="|")
        if d["family"] == "Orchidaceae" and d["taxon_rank"] == "Species"
    ]
    byid = {d["plant_name_id"]: d for d in rows}
    acc = {}
    for d in rows:
        a = (
            d["plant_name_id"]
            if d["taxon_status"] == "Accepted"
            else d["accepted_plant_name_id"]
        )
        if a in byid:
            acc.setdefault(d["taxon_name"], set()).add(byid[a]["taxon_name"])
    return acc


def parse_fasta(text):
    seqs, name = {}, None
    for line in text.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def main(wcvp, traits, out_dir, acc_csv):
    acc = accepted_map(wcvp)
    ours = {r["Binomial"] for r in csv.DictReader(open(traits))}
    org = " OR ".join(f"{g}[Organism]" for g in INGROUP + OUTGROUP)
    rows = []
    for m, (q, min_len) in MARKERS.items():
        ids = json.loads(
            call(
                "esearch.fcgi",
                db="nuccore",
                term=f"({org}) AND {q}",
                retmax=10000,
                retmode="json",
            )
        )["esearchresult"]["idlist"]
        best = {}
        for i in range(0, len(ids), 200):
            chunk = ids[i : i + 200]
            summ = json.loads(
                call("esummary.fcgi", db="nuccore", id=",".join(chunk), retmode="json")
            )["result"]
            for u in summ.get("uids", []):
                s = summ[u]
                organism, title, n = (
                    s.get("organism", ""),
                    s.get("title", ""),
                    int(s.get("slen", 0)),
                )
                if "×" in organism or " x " in organism or "hybrid" in title.lower():
                    continue
                if not (min_len <= n <= MAX_LEN):
                    continue
                if not TITLE[m](title.lower()):
                    continue
                b = " ".join(organism.split()[:2])
                genus = b.split()[0]
                targets = acc.get(b, set())
                if genus in OUTGROUP:
                    sp = b if targets == {b} or not targets else sorted(targets)[0]
                    if sp.split()[0] not in OUTGROUP:
                        continue
                else:
                    hit = sorted(t for t in targets if t in ours)
                    if len(hit) != 1:
                        continue
                    sp = hit[0]
                if sp not in best or n > best[sp][1]:
                    best[sp] = (s["accessionversion"], n, organism, title)
        seqs = parse_fasta(
            call(
                "efetch.fcgi",
                db="nuccore",
                id=",".join(v[0] for v in best.values()),
                rettype="fasta",
                retmode="text",
            )
        )
        with open(f"{out_dir}/{m}.fasta", "w") as f:
            for sp, (a, n, organism, title) in sorted(best.items()):
                seq = seqs[a]
                assert len(seq) == n, (a, len(seq), n)
                f.write(f">{sp.replace(' ', '_')}\n{seq}\n")
                rows.append((m, sp, a, n, organism, title))
        n_ours = sum(sp in ours for sp in best)
        n_out = sum(sp.split()[0] in OUTGROUP for sp in best)
        print(
            f"{m}: {len(ids)} records searched; kept {n_ours}/{len(ours)} study species + {n_out} outgroup species"
        )
    with open(acc_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["marker", "species", "accession", "length", "organism", "title"])
        w.writerows(rows)


if __name__ == "__main__":
    signal.signal(
        signal.SIGALRM,
        lambda *_: (sys.stderr.write("aborting: walltime guard\n"), sys.exit(2)),
    )
    signal.alarm(1200)
    main(*sys.argv[1:5])
