#!/usr/bin/env bash
# Step 1: public inputs. 1000 Genomes high-coverage phased panel (chr22), sample sheet, genetic map.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data maps

KG=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage
PANEL=$KG/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz

wget -nc -P data "$PANEL" "$PANEL.tbi"                                   # 425 MB, 3,202 people
wget -nc -P data "$KG/20130606_g1k_3202_samples_ped_population.txt"      # sex, population, family

wget -nc -P maps https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip
unzip -o -q maps/plink.GRCh38.map.zip -d maps   # maps/chr_in_chrom_field/plink.chrchr22.GRCh38.map ("chr22", as in the panel)

echo "Downloaded: $(ls data maps | wc -l) files"
