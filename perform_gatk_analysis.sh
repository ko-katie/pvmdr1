# !/usr/bin/env bash

#set path to genomicsdb_workspace
genomicsdb_workspace_path="/path/to/genomicsdb_workspace/"

#run GenomicsDBImport to combine individual VCF files listed in sample_map.txt, mask multigene families based on unmasked_regions.intervals
python3 gatk GenomicsDBImport --sample-name-map sample_map.txt --genomicsdb-workspace-path ${genomicsdb_workspace_path} --intervals unmasked_regions.intervals --batch-size 50

#finish combining files into all_samples_variants.output
python3 gatk GenotypeGVCFs -R PlasmoDB-67_PvivaxP01_Genome.fasta -V /gendb:///${genomicsdb_workspace_path} -O all_samples_variants.vcf

#run gatk Variant Filtration
python3 gatk VariantFiltration -R PlasmoDB-67_PvivaxP01_Genome.fasta -V all_samples_variants.vcf -window 35 -cluster 3 --filter-name FS -filter "FS > 60.0" --filter-name QD -filter "QD < 2.0" -O all_samples_variants_filtered.vcf

#run vcftools on filtered variants
vcftools --vcf all_samples_variants_filtered.vcf --remove-indels --max-alleles 2 --max-missing 0.2 --recode --recode-INFO-all --out all_samples_variants_filtered2.vcf
