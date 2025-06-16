# pvmdr1
Code used for characterizing pvmdr1 in manuscript: A common DNA deletion altering the 3’UTR of mdr1 is associated with reduced mefloquine susceptibility in Plasmodium vivax parasites from Cambodian patients

**

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
