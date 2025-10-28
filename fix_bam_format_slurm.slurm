#!/bin/bash

WORKING_DIR_PATH="/path/to/working_dir"      # Set to path to working_dir, must match working_dir path used in fix_bam_format_slurm.sh

#SBATCH --job-name=fix_bam_format                     # Job name
#SBATCH --output=${WORKING_DIR_PATH}/fix_bam.out   # Standard output log
#SBATCH --error=${WORKING_DIR_PATH}/fix_bam.err     # Standard error log
#SBATCH --mail-type=BEGIN,END --mail-user=your@email.com
#SBATCH --cpus-per-task=24				# Number of CPUs per task
#SBATCH --mem=40G                                        # Memory per node
#SBATCH --array=1-X                # Set X to number of samples to be processed

fix_bam_format_slurm.sh

echo "All bams fixed"
