# pvmdr1
Code used for characterizing pvmdr1 in manuscript: A common DNA deletion altering the 3’UTR of _mdr1_ is associated with reduced mefloquine susceptibility in _Plasmodium vivax_ parasites from Cambodian patients

Code developed using: [Hisat2 version 2.2.1](https://daehwankimlab.github.io/hisat2/download/), [gatk version 4.2.2.0](https://github.com/broadinstitute/gatk/releases), [samtools version 1.9](https://www.htslib.org/download/), [R version 4.4.2](https://cran.rstudio.com/), [picard version 2.9.4](https://broadinstitute.github.io/picard/), [vcftools version 0.1.15](https://sourceforge.net/projects/vcftools/), [slurm version 23.11.6](https://slurm.schedmd.com/) 

**Prior to running any screening for deletions and tandem duplications**
- Make sure to generate hisat2 index for reference genome using command:
- hisat2-build PlasmoDB-67_PvivaxP01_Genome.fasta PvivaxP01_hisat_index

**Screening for deletions and tandem duplications in Cambodian patient samples**\
_Screen WGS data from_ P. vivax _patient samples for deletions and tandem duplications of at least 1kb_
- local_wgs_rearrangement_search.py: Gets paths to paired end sequencing data, maps samples to P01 reference genome using Hisat2, and runs local_wgs_flag_search_by_window_final_local.py on resulting sam file to identify possible deletions and duplications
- local_wgs_flag_search_by_window.py: Iterates through sam file from local_wgs_rearrangement_search.py and identifies reads with appropriate flag and insert size that are indicative of deletions or tandem duplications

Perform analysis:
- Ensure list_of_r1_fastqs.txt and list_of_r2_fastqs.txt are in current directory
   - /path/to/list_of_r1_fastqs.txt and /path/to/list_of_r2_fastqs.txt should be text files providing the path to fastq.gz files for R1 and R2, respectively, for each sample. Each file should be on a new line, and the samples should be in the same order in both files

```
python3 sra_wgs_rearrangement_search.py /path/to/working_directory/ list_of_r1_fastqs.txt list_of_r2_fastqs.txt
```
- /path/to/working_directory/ is the path to the directory where files should be outputted


**Screening for _mdr1_ rearrangements from MalariaGEN PV4 Dataset**\
_Screen 60kb region surrounding_ mdr1 _in samples from MalariaGEN PV4 dataset_
- sra_wgs_rearrangement_search.py: Uses MalariaGEN data information with SRA accessions to map samples to P01 reference genome using Hisat2, and runs sra_wgs_flag_search_by_window_final_local.py on resulting sam file to identify possible deletions and duplications
- sra_wgs_flag_search_by_window_final.py: Iterates through sam file from sra_wgs_rearrangement_search.py and identifies reads with appropriate flag and insert size that are indicative of deletions or tandem duplications

Perform analysis:

- Download Pv4_samples.txt from [MalariaGen PV4 dataset page](https://www.malariagen.net/data_package/open-dataset-plasmodium-vivax-v4-0/) (click "Download sample provenance and sequencing metadata") and ensure in current directory
  
```
python3 sra_wgs_rearrangement_search.py /path/to/working_directory/ Pv4_samples.txt
```
- /path/to/working_directory/ is the path to the directory where files should be outputted

**GATK Analysis**\
_Calculate fws to determine clonality of samples based on WGS_
- fix_bam_format_slurm.slurm: Used to run fix_bam_format_slurm.sh as array job in SLURM workload manager. Allows fix_bam_format_slurm.sh to run on multiple samples in parallel
- fix_bam_format_slurm.sh: Reformats BAM files in preparation for GATK and VCFTools analysis, then performs HaplotypeCaller to generate VCF file for sample
- perform_gatk_analysis.sh: Bash script containing commands to combine outputted VCF files from fix_bam_format_slurm.sh and generate VCF file representing variants between all samples
- unmasked_intervals.intervals: .intervals file for use with perform_gatk_analysis.sh, contains regions that were not masked in analysis (ie does not belong to multigene family)
- moimix_analysis: Takes final VCF file from perform_gatk_analysis.sh and calculates FWS of samples

Perform analysis:

- Set line 4 of fix_bam_format_slurm.sh and line 3 of fix_bam_format_slurm.slurm to path for desired working directory
- Ensure PlasmoDB-67_PvivaxP01_Genome.fasta, unmasked_regions.intervals, all_bam_paths.txt, and sample_map.txt are in current directory
  - all_bam_paths.txt should contain paths to bam files to be processed, each on a new line
  - sample_map.txt should be formatted as outlined [here](https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport)
- Set X in line 11 of fix_bam_format_slurm.slurm to number of samples in all_bam_paths.txt
- Set line 3 of perform_gatk_analysis.sh to desired genomicsdb_workspace for running gatk, should be empty or nonexistent

```
sbatch --mem=44G fix_bam_format_slurm.slurm
bash perform_gatk_analysis.sh
```

- Run moimix_analysis.r using Rstudio
  
**Generate mdr1 coding sequence consensus**\
_Get consensus sequence of mdr1 coding sequence from WGS of samples_
- run_mpileup_mdr1_slurm.slurm: Used to run mpileup_mdr1_slurm.sh as array job in SLURM workload manager. Allows mpileup_mdr1_slurm.sh to run on multiple samples in parallel
- mpileup_mdr1_slurm.sh: Runs mpileup on region of BAM file pertaining to mdr1 coding sequence, then pipes to mpileup_to_fasta.py
- mpileup_to_fasta.py: Python script to process VCF file outputted by mpileup and rebuild mdr1 coding sequence consensus

Perform analysis:

- Set run_mpileup_mdr1_slurm.slurm line 3 to desired path for working directory and X in line 10 to number of samples to be run
- Set mpileup_mdr1_slurm.sh line 4 to include path to desired working directory
- Ensure samples_to_run.txt and PlasmoDB-67_PvivaxP01_Genome.fasta files are in current directory
   - samples_to_run.txt should be a tab-delimited file with each line containing one sample, sample name in first column and path to sample bam file in second column
