#!/usr/bin/env bash
# Step 5 (optional): impute straight from low-coverage reads, no variant calling.
# 20 of the test people (lowcov_samples.txt) were sequenced at 30x; we keep a random 1/30 of their chr22 reads (~1x) and impute from those.
set -euo pipefail
cd "$(dirname "$0")/.."
T=${THREADS:-8}
SELPHI=${SELPHI:-selphi}
COV=${COV:-1}                       # target coverage: 0.5, 1, 2 ...
MAP=maps/chr_in_chrom_field/plink.chrchr22.GRCh38.map
KG=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage
REF=ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
mkdir -p ref reads results

# The genome the 1000 Genomes CRAMs were aligned to (3.1 GB, once).
REFURL=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/$(basename $REF)
wget -nc -P ref "$REFURL" "$REFURL.fai"

# Where each person's 30x CRAM lives.
curl -fsS "$KG/1000G_2504_high_coverage.sequence.index" | grep -v '^#' |
  awk -F'\t' '{u = $1; sub(/^ftp:\/\//, "https://", u); print $10 "\t" u}' > reads/cram_urls.tsv
[[ -s reads/cram_urls.tsv ]] || { echo "Could not download the CRAM index (see the README for CURL_CA_BUNDLE)"; exit 1; }

# Each person's chr22 reads at 30x, streamed from the CRAM (about 1 GB each, once),
# then a random 1/30*COV of the read pairs kept (samtools -s SEED.FRACTION).
FRAC=$(awk -v c="$COV" 'BEGIN{printf "%.4f", c / 30}')
while read -r s; do
  full=reads/$s.chr22.30x.cram
  out=reads/$s.chr22.${COV}x.cram
  if [[ ! -s $full ]]; then
    url=$(awk -v s="$s" '$1 == s {print $2}' reads/cram_urls.tsv)
    [[ -n $url ]] || { echo "No CRAM listed for $s" >&2; exit 1; }
    samtools view -T "$REF" -C -o "$full.part" "$url" chr22 && mv "$full.part" "$full"
  fi
  samtools view -T "$REF" -s "7${FRAC#0}" -C -o "$out" "$full"
  samtools index "$out"
  echo "$out"
done < lowcov_samples.txt > reads/crams_${COV}x.txt

# Impute the 20 from their reads in one run, then score.
"$SELPHI" --bam-list reads/crams_${COV}x.txt --reference "$REF" --refpanel work/panel.chr22.srp2 \
  --map "$MAP" --out results/lowcov_${COV}x.chr22 --threads "$T"
"$SELPHI" --evaluate results/lowcov_${COV}x.chr22.vcf.gz --truth work/truth.chr22.vcf.gz \
  --out results/lowcov_${COV}x.eval
