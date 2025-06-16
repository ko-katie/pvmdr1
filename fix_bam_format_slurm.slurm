#!/bin/bash
#SBATCH --job-name=fix_bam_format                     # Job name
#SBATCH --output=/local/projects-t3/SerreDLab-3/kko/mdr_del_search_dna_whole_genome/gatk_polyclonality/fix_bam_files_dir/fix_bam.out   # Standard output log
#SBATCH --error=/local/projects-t3/SerreDLab-3/kko/mdr_del_search_dna_whole_genome/gatk_polyclonality/fix_bam_files_dir/fix_bam.err     # Standard error log
#SBATCH --mail-type=BEGIN,END --mail-user=katie.ko@som.umaryland.edu
#SBATCH --cpus-per-task=24				# Number of CPUs per task
#SBATCH --mem=40G                                        # Memory per node
#SBATCH --array=1-232

/local/projects-t3/SerreDLab-3/kko/mdr_del_search_dna_whole_genome/gatk_polyclonality/fix_bam_format_slurm.sh

echo "All bams fixed"