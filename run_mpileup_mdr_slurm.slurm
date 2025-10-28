#!/bin/bash

WORKING_DIR_PATH="/path/to/working_dir"

#SBATCH --job-name=mp_mdr1                     # Job name
#SBATCH --output=${WORKING_DIR_PATH}/mpileup.out   # Standard output log
#SBATCH --error=${WORKING_DIR_PATH}/mpileup.err    # Standard error log
#SBATCH --cpus-per-task=6				# Number of CPUs per task
#SBATCH --mem=12G                                       # Memory per node
#SBATCH --array=1-X                       # Set X to number of samples to be processed


mpileup_mdr_slurm.sh

echo "All samples processed"
