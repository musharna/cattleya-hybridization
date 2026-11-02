#!/usr/bin/env bash
# Molecular Cattleya trees for prereg C, Amendment 2: align each marker, build combined, ITS-only and
# plastid-only ML trees with ultrafast bootstraps, rooted on the Laeliinae outgroups.
# Usage: analysis/molecular_tree.sh SEQ_DIR OUT_DIR   (SEQ_DIR holds ITS.fasta, matK.fasta, trnLF.fasta
#   from molecular_fetch.py). Output: OUT_DIR/{combined,ITS,plastid}.treefile and .ufboot, plus logs.
set -euo pipefail
SEQ=$1; OUT=$2
BIN=$HOME/miniconda3/envs/orthofinder/bin
PY=$HOME/miniconda3/envs/plantcad2/bin/python
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT"
for m in ITS matK trnLF; do
  "$BIN/mafft" --auto --adjustdirection --thread 4 "$SEQ/$m.fasta" > "$OUT/$m.aln.fasta" 2> "$OUT/$m.mafft.log"
done
"$PY" -I "$HERE/molecular_concat.py" "$OUT/combined" ITS="$OUT/ITS.aln.fasta" matK="$OUT/matK.aln.fasta" trnLF="$OUT/trnLF.aln.fasta"
"$PY" -I "$HERE/molecular_concat.py" "$OUT/ITS" ITS="$OUT/ITS.aln.fasta"
"$PY" -I "$HERE/molecular_concat.py" "$OUT/plastid" matK="$OUT/matK.aln.fasta" trnLF="$OUT/trnLF.aln.fasta"
for t in combined ITS plastid; do
  og=$(grep '^>' "$OUT/$t.fasta" | sed 's/^>//' | grep -E '^(Guarianthe|Brassavola|Rhyncholaelia)_' | paste -sd, -)
  [ -n "$og" ] || { echo "no outgroup taxa in $t" >&2; exit 1; }
  # no -o: IQ-TREE 3.1.3 aborts in setRootNode with a multi-taxon outgroup (job 7417); -o only sets the
  # drawn root, and molecular_chronos.R roots on these outgroups anyway
  "$BIN/iqtree" -s "$OUT/$t.fasta" -p "$OUT/$t.nex" -m MFP -B 1000 --wbtl --seed 1 -T 4 \
    --prefix "$OUT/$t" -redo > "$OUT/$t.iqtree.stdout"
done
echo done
