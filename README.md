# pvmdr1
Code used for characterizing pvmdr1 in manuscript: A common DNA deletion altering the 3’UTR of _mdr1_ is associated with reduced mefloquine susceptibility in _Plasmodium vivax_ parasites from Cambodian patients

**Screening for deletions and tandem duplications in Cambodian patient samples**\
_Screen WGS data from_ P. vivax _patient samples for deletions and tandem duplications of at least 1kb_
- local_wgs_rearrangement_search.py: Gets paths to paired end sequencing data, maps samples to P01 reference genome using Hisat2, and runs local_wgs_flag_search_by_window_final_local.py on resulting sam file to identify possible deletions and duplications
- local_wgs_flag_search_by_window.py: Iterates through sam file from local_wgs_rearrangement_search.py and identifies reads with appropriate flag and insert size that are indicative of deletions or tandem duplications

  Peform analysis by running: python3 local_wgs_rearrangement_search.py /path/to/working_directory/ /path/to/list_of_r1_fastqs.txt /path/to/list_of_r2_fastqs.txt\
  **/path/to/working_directory/** is the path to the directory where files should be outputted\
  **/path/to/list_of_r1_fastqs.txt** and **/path/to/list_of_r2_fastqs.txt** should be text files providing the path to fastq.gz files for R1 and R2, respectively, for each sample. Each file should be on a new line, and the samples should be in the same order in both files.\

**Screening for _mdr1_ rearrangements from MalariaGEN PV4 Dataset**\
_Screen 60kb region surrounding_ mdr1 _in samples from MalariaGEN PV4 dataset_
- sra_wgs_rearrangement_search.py: Uses MalariaGEN data information with SRA accessions to map samples to P01 reference genome using Hisat2, and runs sra_wgs_flag_search_by_window_final_local.py on resulting sam file to identify possible deletions and duplications
- sra_wgs_flag_search_by_window_final.py: Iterates through sam file from sra_wgs_rearrangement_search.py and identifies reads with appropriate flag and insert size that are indicative of deletions or tandem duplications

  Peform analysis by running: python3 sra_wgs_rearrangement_search.py /path/to/working_directory/ /path/to/list_of_r1_fastqs.txt /path/to/list_of_r2_fastqs.txt\
  **/path/to/working_directory/** is the path to the directory where files should be outputted\
  **/path/to/list_of_r1_fastqs.txt** and **/path/to/list_of_r2_fastqs.txt** should be text files providing the path to fastq.gz files for R1 and R2, respectively, for each sample. Each file should be on a new line, and the samples should be in the same order in both files.\

**GATK Analysis**\
_Calculate fws to determine clonality of samples based on WGS_
- fix_bam_format_slurm.slurm: Used to run fix_bam_format_slurm.sh as array job in SLURM workload manager. Allows fix_bam_format_slurm.sh to run on multiple samples in parallel
- fix_bam_format_slurm.sh: Reformats BAM files in preparation for GATK and VCFTools analysis, then performs HaplotypeCaller to generate VCF file for sample
- perform_gatk_analysis.sh: Bash script containing commands to combine outputted VCF files from fix_bam_format_slurm.sh and generate VCF file representing variants between all samples
- moimix_analysis: Takes final VCF file from perform_gatk_analysis.sh and calculates FWS of samples

    Perform analysis by:\
    Change fix_bam_format_slurm.slurm lines 3, 4, and 10 to desired paths to SLURM outputs and with downloaded location of fix_bam_format_slurm.sh script\
    Change fix_bam_format_slurm.sh lines 4, 5, 6, and 7 to desired paths to SLURM outputs and line 8 to path to a list of bam files to be run\
      Bam file should contain one path on each line\
    Change /path/to/sample_map.txt on line 4 of perform_gatk_analysis.sh to sample map as outlined [here](https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport), then run: bash perform_gatk_analysis.sh\
    Change paths in moimix_analysis.sh to reflect those outputted by perform_gatk_analysis.sh and run to obtain fws values\

**Generate mdr1 coding sequence consensus**\
_Get consensus sequence of mdr1 coding sequence from WGS of samples_
- run_mpileup_mdr1_slurm.slurm: Used to run mpileup_mdr1_slurm.sh as array job in SLURM workload manager. Allows mpileup_mdr1_slurm.sh to run on multiple samples in parallel
- mpileup_mdr1_slurm.sh: Runs mpileup on region of BAM file pertaining to mdr1 coding sequence, then pipes to mpileup_to_fasta.py
- mpileup_to_fasta.py: Python script to process VCF file outputted by mpileup and rebuild mdr1 coding sequence consensus

    Perform analysis by:\
    Change run_mpileup_mdr1_slurm.slurm lines 3 and 4 to desired paths to SLURM output and line 10 to downloaded location of mpileup_mdr1_slurm.sh\
    Change mpileup_mdr1_slurm.sh line 4 to include path to working directory and line 11 to include path to input file\
    Input file should be a tab delimited file with each line containing one sample, sample name in first column and path to sample bam file in second column\
