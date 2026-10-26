#!/usr/bin/env bash
# Build the RIsearch2 suffix-array index that target prediction needs.
#
# The repository ships the reference FASTA (GRCh38.CDS.fa, GRCh38.3UTR.fa) but
# not the .suf index: the CDS index alone is ~500 MB, too large to keep in git.
# Build it once into a mounted directory and reuse it on every later run.
#
# Inside the container:
#   docker run --rm -v "$PWD/data:/data" heracles:1.0 \
#     bash /opt/heracles/docker/build_index.sh /data/RIsearch2_index
#
# Without a container:
#   bash docker/build_index.sh ./RIsearch2_index
#
# Existing .suf files are skipped, so the script is safe to re-run.

set -euo pipefail

OUT=${1:-/data/RIsearch2_index}
HERE=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$OUT"

echo "index output directory: $OUT"
echo

for ref in CDS 3UTR; do
    fa="$HERE/RIsearch2_index/GRCh38.${ref}.fa"
    suf="$OUT/GRCh38.${ref}.suf"

    if [ ! -s "$fa" ]; then
        echo "  skip ${ref}: reference FASTA not found at $fa"
        continue
    fi
    if [ -s "$suf" ]; then
        echo "  skip ${ref}: $suf already exists"
        continue
    fi

    echo "  building ${ref} index (this takes several minutes for CDS) ..."
    RIsearch2 -c "$fa" -o "$suf"
    # RIsearch2 exits 0 even when it writes nothing; check the result.
    if [ ! -s "$suf" ]; then
        echo "  FAILED ${ref}: no index written" >&2
        exit 1
    fi
    echo "  done ${ref}: $(du -h "$suf" | cut -f1)"
done

echo
echo "Point --predict-index at:"
for ref in CDS 3UTR; do
    [ -s "$OUT/GRCh38.${ref}.suf" ] && echo "  $OUT/GRCh38.${ref}.suf"
done
