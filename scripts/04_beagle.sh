#!/usr/bin/env bash
# Step 4 (optional): the same job with Beagle 5.5, scored the same way, to compare.
set -euo pipefail
cd "$(dirname "$0")/.."
T=${THREADS:-8}
MAP=maps/chr_in_chrom_field/plink.chrchr22.GRCh38.map
BEAGLE=tools/beagle.27Feb25.75f.jar
BREF3=tools/bref3.27Feb25.75f.jar
mkdir -p tools results
wget -nc -P tools "https://faculty.washington.edu/browning/beagle/$(basename $BEAGLE)" \
                  "https://faculty.washington.edu/browning/beagle/$(basename $BREF3)"

# 4a. Panel in Beagle's binary format.
java -jar "$BREF3" work/panel.chr22.vcf.gz > work/panel.chr22.bref3

# 4b. Phase + impute.
java -Xmx16g -jar "$BEAGLE" gt=work/array.chr22.vcf.gz ref=work/panel.chr22.bref3 map="$MAP" \
  out=results/beagle.chr22 nthreads="$T"
bcftools index -t results/beagle.chr22.vcf.gz

# 4c. Same scoring as step 3.
"${SELPHI:-selphi}" --evaluate results/beagle.chr22.vcf.gz --truth work/truth.chr22.vcf.gz \
  --exclude-sites work/array.chr22.vcf.gz --out results/beagle.eval
