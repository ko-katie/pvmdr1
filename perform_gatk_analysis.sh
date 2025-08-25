# !/usr/bin/env bash

#run GenomicsDBImport to combine individual VCF files listed in sample_map.txt, mask multigene families based on unmasked_regions.intervals
python3 /path/to/packages/gatk-4.2.2.0/gatk GenomicsDBImport --sample-name-map /path/to/sample_map.txt --genomicsdb-workspace-path /path/to/fws_working_dir/genomicsdb_workspace/ --intervals /path/to/fws_working_dir/unmasked_regions.intervals --batch-size 50

#finish combining files into all_samples_variants.output
python3 /usr/local/packages/gatk-4.2.2.0/gatk GenotypeGVCFs -R /path/to/fws_working_dir/PlasmoDB-67_PvivaxP01_Genome.fasta -V /gendb:///path/togenomicsdb_database -O /path/to/fws_working_dir/all_samples_variants.output

#run gatk Variant Filtration
python3 /usr/local/packages/gatk-4.2.2.0/gatk VariantFiltration -R /path/to/PlasmoDB-67_PvivaxP01_Genome.fasta -V /path/to/fws_working_dir/all_samples_variants.vcf -window 35 -cluster 3 --filter-name FS -filter "FS > 60.0" --filter-name QD -filter "QD < 2.0" -O /path/to/fws_working_dir/all_samples_variants_filtered.vcf

#run vcftools on filtered variants
/usr/local/packages/vcftools/bin/vcftools --vcf /path/to/fws_working_dir/all_samples_variants_filtered.vcf --remove-indels --max-alleles 2 --max-missing 0.2 --recode --recode-INFO-all --out /path/to/fws_working_dir/all_samples_variants_filtered2.vcf

gzip /path/to/all_samples_variants_filtered2.vcf
