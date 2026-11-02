"""Concatenate per-marker alignments into one supermatrix with a NEXUS partition file.

What: reads aligned FASTA files (mafft output; a "_R_" prefix from --adjustdirection is stripped),
  pads taxa missing from a marker with gaps, and writes the supermatrix plus charsets.
When: between alignment and IQ-TREE in analysis/molecular_tree.sh (prereg C, Amendment 2).
Output: OUT.fasta and OUT.nex; prints taxa and sites per marker. Duplicate taxon names fail loudly.
Usage: python -I analysis/molecular_concat.py OUT_PREFIX MARKER=ALN.fasta [MARKER=ALN.fasta ...]
"""

import sys


def read_aln(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            if name.startswith("_R_"):
                name = name[3:]
            assert name not in seqs, f"duplicate taxon {name} in {path}"
            seqs[name] = []
        elif name:
            seqs[name].append(line)
    seqs = {k: "".join(v).upper() for k, v in seqs.items()}
    lens = {len(v) for v in seqs.values()}
    assert len(lens) == 1, f"{path}: unequal aligned lengths {sorted(lens)[:5]}"
    return seqs, lens.pop()


def main(prefix, specs):
    parts = []
    for spec in specs:
        m, path = spec.split("=", 1)
        seqs, n = read_aln(path)
        parts.append((m, seqs, n))
        print(f"{m}: {len(seqs)} taxa, {n} aligned sites")
    taxa = sorted(set().union(*(p[1] for p in parts)))
    with open(prefix + ".fasta", "w") as f:
        for t in taxa:
            f.write(f">{t}\n" + "".join(s.get(t, "-" * n) for _, s, n in parts) + "\n")
    with open(prefix + ".nex", "w") as f:
        f.write("#nexus\nbegin sets;\n")
        start = 1
        for m, _, n in parts:
            f.write(f"  charset {m} = {start}-{start + n - 1};\n")
            start += n
        f.write("end;\n")
    print(f"supermatrix: {len(taxa)} taxa, {start - 1} sites")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
