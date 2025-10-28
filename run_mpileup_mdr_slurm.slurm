#!/bin/bash
#SBATCH --job-name=mp_mdr1                     # Job name
#SBATCH --output=/path/to/mpileup_dir/mpileup.out   # Standard output log
#SBATCH --error=/path/to/mpileup_dir/mpileup.err    # Standard error log
#SBATCH --cpus-per-task=6				# Number of CPUs per task
#SBATCH --mem=12G                                       # Memory per node
#SBATCH --array=1-112                       #Number of samples


mpileup_mdr_slurm.sh

echo "All samples processed"
