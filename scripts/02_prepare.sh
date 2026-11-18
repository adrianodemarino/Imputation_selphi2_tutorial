#!/usr/bin/env bash
# Step 2: split the 1000 Genomes people into a reference panel and 100 test people,
# then make an "array" for the test people: what a genotyping chip would have read.
set -euo pipefail
cd "$(dirname "$0")/.."
T=${THREADS:-8}
IN=data/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
mkdir -p work

# Reference panel: everyone except the 100 test people (none of them has a parent or child in 1000 Genomes),
# biallelic SNPs and indels, sites still polymorphic once the test people are out.
bcftools view --threads "$T" -S ^test_samples.txt -v snps,indels -m2 -M2 "$IN" -Ou |
  bcftools view --threads "$T" -c1 -Oz -o work/panel.chr22.vcf.gz && bcftools index -t work/panel.chr22.vcf.gz

# Truth: the 100 test people at every site, as sequenced at 30x.
bcftools view --threads "$T" -S test_samples.txt -v snps,indels -m2 -M2 "$IN" -Oz -o work/truth.chr22.vcf.gz && bcftools index -t work/truth.chr22.vcf.gz

# Array sites: common SNPs (panel MAF >= 1%), at most one every 5 kb, like a genotyping chip (~7,000 on chr22).
bcftools view -v snps -q 0.01:minor work/panel.chr22.vcf.gz -Ou |
  bcftools query -f '%CHROM\t%POS\n' |
  awk 'BEGIN{OFS="\t"} $2 >= last + 5000 {print; last = $2}' > work/array_sites.chr22.tsv

# The test people's genotypes at those sites, unphased, as a chip reports them.
bcftools view -T work/array_sites.chr22.tsv work/truth.chr22.vcf.gz -Ou |
  bcftools +setGT -Ou - -- -t a -n u |
  bcftools annotate -x INFO -Oz -o work/array.chr22.vcf.gz && bcftools index -t work/array.chr22.vcf.gz

for f in panel truth array; do
  printf '%-6s %6s people  %9s sites\n' "$f" \
    "$(bcftools query -l work/$f.chr22.vcf.gz | wc -l)" \
    "$(bcftools index -n work/$f.chr22.vcf.gz)"
done
