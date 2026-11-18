#!/usr/bin/env bash
# Step 3: phase and impute the 100 test people with Selphi 2, then score the result against the truth.
set -euo pipefail
cd "$(dirname "$0")/.."
T=${THREADS:-8}
SELPHI=${SELPHI:-selphi}
MAP=maps/chr_in_chrom_field/plink.chrchr22.GRCh38.map
mkdir -p results

# 3a. Convert the panel once into Selphi's own format (.srp2). Reuse it for every run after this.
"$SELPHI" --prepare-reference-from work/panel.chr22.vcf.gz --out work/panel.chr22.srp2 --threads "$T"

# 3b. Phase + impute. Writes results/selphi2.chr22.vcf.gz (GT, DS, AP1, AP2; INFO has AF and the imputation R2).
"$SELPHI" --refpanel work/panel.chr22.srp2 --input work/array.chr22.vcf.gz --map "$MAP" \
  --out results/selphi2.chr22 --threads "$T"

# 3c. Score against the 30x truth, imputed sites only (the array sites are left out).
"$SELPHI" --evaluate results/selphi2.chr22.vcf.gz --truth work/truth.chr22.vcf.gz \
  --exclude-sites work/array.chr22.vcf.gz --out results/selphi2.eval
