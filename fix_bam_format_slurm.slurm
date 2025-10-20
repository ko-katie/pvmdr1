#!/bin/bash
#SBATCH --job-name=fix_bam_format                     # Job name
#SBATCH --output=/path/to/fix_bam_files_dir/fix_bam.out   # Standard output log
#SBATCH --error=/path/to/fix_bam_files_dir/fix_bam.err     # Standard error log
#SBATCH --mail-type=BEGIN,END --mail-user=your@email.com
#SBATCH --cpus-per-task=24				# Number of CPUs per task
#SBATCH --mem=40G                                        # Memory per node
#SBATCH --array=1-232

/path/to/fix_bam_format_slurm.sh

echo "All bams fixed"
