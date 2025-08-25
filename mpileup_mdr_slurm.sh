#!/bin/bash

# Paths and settings
WORKING_DIR="/path/to/mpileup_dir/"  # Directory for running script

# Create directory
mkdir -p "$WORKING_DIR" 

# Get the sample line for the slurm array task ID
#samples_to_run.txt is formatted as a tab-delimited file with the sample name in the first column and the path to the bam file in the second column
SAMPLE_LINE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /path/to/mpileup_dir/samples_to_run.txt)

echo "Processing sample: $SAMPLE_LINE"

# Get strings in tab-separated line
LINE_ARRAY=( $( echo "$SAMPLE_LINE" | awk -F'\t' '{print $1, $2}' ) )
SAMPLE_NAME=${LINE_ARRAY[0]}
SAMPLE_BAM=${LINE_ARRAY[1]}

# Check if bam file exists
if [ ! -f "$SAMPLE_BAM" ]; then
    echo "Skipping $SAMPLE_NAME: Missing BAM file"
    exit 1
fi

echo "SAMPLE_NAME: $SAMPLE_NAME"
echo "SAMPLE_BAM: $SAMPLE_BAM"

SAMPLE_VCF_DIR="${WORKING_DIR}${SAMPLE_NAME}_VCF_dir/"

mkdir -p "$SAMPLE_VCF_DIR"

CONSENSUS_FASTA="${WORKING_DIR}${SAMPLE_NAME}_consensus.fasta"

#run mpileup on bam file for mdr1 cds region and pipe to script which will make consensus sequence
bcftools mpileup -Ov -r PvP01_10_v2:478539-483333 -f /local/Pv_P01_Index/PlasmoDB-67_PvivaxP01_Genome.fasta $SAMPLE_BAM | python3 /local/mpileup_dir/mpileup_to_fasta.py $CONSENSUS_FASTA $SAMPLE_NAME

rm -r $SAMPLE_VCF_DIR

echo "Processing completed for $SAMPLE_NAME."

exit 0
