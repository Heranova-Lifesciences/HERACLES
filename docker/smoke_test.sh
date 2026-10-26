#!/usr/bin/env bash
# Smoke test for a HERACLES installation (container or conda env).
# Usage:  bash smoke_test.sh            # checks tools only
#         bash smoke_test.sh /data      # also checks that reference data exist
set -u

fail=0
check_out() {  # check_out <label> <expected-substring> <command...>
    local label=$1 want=$2; shift 2
    if "$@" 2>&1 | grep -q "$want"; then
        printf '  ok    %s\n' "$label"
    else
        printf '  FAIL  %s\n' "$label"
        fail=$((fail+1))
    fi
}

check() {  # check <label> <command...>
    local label=$1; shift
    if "$@" >/dev/null 2>&1; then
        printf '  ok    %s\n' "$label"
    else
        printf '  FAIL  %s\n' "$label"
        fail=$((fail+1))
    fi
}

echo "External tools"
check "bowtie"       bowtie --version
check "cutadapt"     cutadapt --version
check "trim_galore"  trim_galore --version
check_out "RIsearch2" "RIsearch v" RIsearch2

echo "Python packages"
check "pandas"    python -c "import pandas"
check "numpy"     python -c "import numpy"
check "biopython" python -c "import Bio"
check "networkx"  python -c "import networkx"
check "matplotlib" python -c "import matplotlib"
check "seaborn"   python -c "import seaborn"
check "pydeseq2"  python -c "import pydeseq2"
check "gseapy"    python -c "import gseapy"

echo "HERACLES entry point"
HC=${HERACLES_HOME:-$(cd "$(dirname "$0")/.." && pwd)}
check "main.py --help" python "$HC/main.py" --help

if [ $# -ge 1 ]; then
    D=$1
    echo "Reference data under $D"
    check "tRNA fasta"      test -s "$D/tRNA_index/hg38-tRNA.fa"
    check "bowtie index"    test -s "$D/tRNA_index/mature.1.ebwt"
    check "CDS index (.suf)" test -s "$D/RIsearch2_index/GRCh38.CDS.suf"
fi

echo
if [ "$fail" -eq 0 ]; then
    echo "all checks passed"
    exit 0
else
    echo "$fail check(s) failed"
    exit 1
fi
