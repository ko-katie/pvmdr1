#!/bin/bash
#SBATCH --job-name=mp_mdr1                     # Job name
#SBATCH --output=/local/mpileup_dir/mpileup.out   # Standard output log
#SBATCH --error=/local/mpileup_dir/mpileup.err    # Standard error log
#SBATCH --cpus-per-task=6				# Number of CPUs per task
#SBATCH --mem=12G                                       # Memory per node
#SBATCH --array=1-112                       #Number of samples


/local/mpileup_dir/mpileup_mdr_slurm.sh

echo "All samples processed"