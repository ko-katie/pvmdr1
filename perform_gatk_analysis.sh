# !/usr/bin/env bash

#run GenomicsDBImport to combine individual VCF files listed in sample_map.txt, mask multigene families based on unmasked_regions.intervals
python3 /usr/local/packages/gatk-4.2.2.0/gatk GenomicsDBImport --sample-name-map /local/fws_working_dir/sample_map.txt --genomicsdb-workspace-path /local/fws_working_dir/genomicsdb_workspace/ --intervals /local/fws_working_dir/unmasked_regions.intervals --batch-size 50

#finish combining files into all_samples_variants.output
python3 /usr/local/packages/gatk-4.2.2.0/gatk GenotypeGVCFs -R /local/fws_working_dir/PlasmoDB-67_PvivaxP01_Genome.fasta -V /gendb:///local/fws_working_dir/genomicsdb_database -O /local/fws_working_dir/all_samples_variants.output

#run gatk Variant Filtration
python3 /usr/local/packages/gatk-4.2.2.0/gatk VariantFiltration -R /local/fws_working_dir/PlasmoDB-67_PvivaxP01_Genome.fasta -V /local/fws_working_dir/all_samples_variants.vcf -window 35 -cluster 3 --filter-name FS -filter "FS > 60.0" --filter-name QD -filter "QD < 2.0" -O /local/projects-t3/SerreDLab-3/kko/mdr_del_search_dna_whole_genome/gatk_polyclonality/all_samples_variants_filtered.vcf

#run vcftools on filtered variants
/usr/local/packages/vcftools/bin/vcftools --vcf /local/fws_working_dir/all_samples_variants_filtered.vcf --remove-indels --max-alleles 2 --max-missing 0.2 --recode --recode-INFO-all --out /local/fws_working_dir/all_samples_variants_filtered2.vcf

gzip /local/fws_working_dir/all_samples_variants_filtered2.vcf
