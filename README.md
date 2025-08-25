# pvmdr1
Code used for characterizing pvmdr1 in manuscript: A common DNA deletion altering the 3’UTR of _mdr1_ is associated with reduced mefloquine susceptibility in _Plasmodium vivax_ parasites from Cambodian patients

**Screening for deletions and tandem duplications in Cambodian patient samples**\
_Screen WGS data from_ P. vivax _patient samples for deletions and tandem duplications of at least 1kb_
- local_wgs_rearrangement_search.py: Gets paths to paired end sequencing data, maps samples to P01 reference genome using Hisat2, and runs local_wgs_flag_search_by_window_final_local.py on resulting sam file to identify possible deletions and duplications
- local_wgs_flag_search_by_window.py: Iterates through sam file from local_wgs_rearrangement_search.py and identifies reads with appropriate flag and insert size that are indicative of deletions or tandem duplications

**Screening for _mdr1_ rearrangements from MalariaGEN PV4 Dataset**\
_Screen 60kb region surrounding_ mdr1 _in samples from MalariaGEN PV4 dataset_
- sra_wgs_rearrangement_search.py: Uses MalariaGEN data information with SRA accessions to map samples to P01 reference genome using Hisat2, and runs sra_wgs_flag_search_by_window_final_local.py on resulting sam file to identify possible deletions and duplications
- sra_wgs_flag_search_by_window_final.py: Iterates through sam file from sra_wgs_rearrangement_search.py and identifies reads with appropriate flag and insert size that are indicative of deletions or tandem duplications

**GATK Analysis**\
_Calculate fws to determine clonality of samples based on WGS_
- fix_bam_format_slurm.slurm: Used to run fix_bam_format_slurm.sh as array job in SLURM workload manager. Allows fix_bam_format_slurm.sh to run on multiple samples in parallel
- fix_bam_format_slurm.sh: Reformats BAM files in preparation for GATK and VCFTools analysis, then performs HaplotypeCaller to generate VCF file for sample
- perform_gatk_analysis.sh: Bash script containing commands to combine outputted VCF files from fix_bam_format_slurm.sh and generate VCF file representing variants between all samples
- moimix_analysis: Takes final VCF file from perform_gatk_analysis.sh and calculates FWS of samples

**Generate mdr1 coding sequence consensus**\
_Get consensus sequence of mdr1 coding sequence from WGS of samples_
- run_mpileup_mdr1_slurm.slurm: Used to run mpileup_mdr1_slurm.sh as array job in SLURM workload manager. Allows mpileup_mdr1_slurm.sh to run on multiple samples in parallel
- mpileup_mdr1_slurm.sh: Runs mpileup on region of BAM file pertaining to mdr1 coding sequence, then pipes to mpileup_to_fasta.py
- mpileup_to_fasta.py: Python script to process VCF file outputted by mpileup and rebuild mdr1 coding sequence consensus
