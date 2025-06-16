#!/bin/bash

# Paths and settings
WORKING_DIR="/local/projects-t3/SerreDLab-3/kko/mdr_del_search_dna_whole_genome/gatk_polyclonality/fix_bam_files_dir/"  # Directory for fixing BAM files
SUM_DIR="/local/projects-t3/SerreDLab-3/kko/mdr_del_search_dna_whole_genome/gatk_polyclonality/fix_bam_files_dir/sum_dir"
LOG_DIR="/local/projects-t3/SerreDLab-3/kko/mdr_del_search_dna_whole_genome/gatk_polyclonality/fix_bam_files_dir/log_dir"
TMP_DIR="/local/projects-t3/SerreDLab-3/kko/mdr_del_search_dna_whole_genome/gatk_polyclonality/fix_bam_files_dir/tmp"
BAM_LIST_FILE="/local/projects-t3/SerreDLab-3/kko/mdr_del_search_dna_whole_genome/all_bam_paths.txt"

# Create necessary directories
mkdir -p "$WORKING_DIR"
mkdir -p "$SUM_DIR"
mkdir -p "$LOG_DIR"
mkdir -p "$TMP_DIR"

# Get the subfolder name from the corresponding line in BAM_LIST_FILE
BAM_FILE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $BAM_LIST_FILE)

if [ -n "$BAM_FILE" ]; then

    echo "Bam file name: $BAM_FILE"

    #get sample name from bam file path
    if [[ "$BAM_FILE" =~ dir/(.+).sorted.bam ]]; then
        SAMPLE_NAME=${BASH_REMATCH[1]}
        echo "Sample name: $SAMPLE_NAME"
    else
        echo "Could not extract sample name"
    fi

    #make file paths for fixed bam files
    MATE_BAM=${WORKING_DIR}${SAMPLE_NAME}"_mate.sorted.bam"
    MATE_RG_BAM=${WORKING_DIR}${SAMPLE_NAME}"_mate_rg.sorted.bam"

    #run picard functions to fix bam files
    java -jar /usr/local/packages/picard/picard.jar FixMateInformation I=$BAM_FILE O=$MATE_BAM tmp_dir=$TMP_DIR
    java -jar /usr/local/packages/picard/picard.jar AddOrReplaceReadGroups I=$MATE_BAM O=$MATE_RG_BAM RGLB=laneX RGPL=illumina RGPU=NONE RGSM=$SAMPLE_NAME

    # Remove intermediate BAM
    rm $MATE_BAM

    #make bam index file
    INDEX_FILE=${WORKING_DIR}${SAMPLE_NAME}_mate_rg.sorted.bai
    samtools index $MATE_RG_BAM $INDEX_FILE

    #run HaplotypeCaller
    VCF_PATH=${WORKING_DIR}${SAMPLE_NAME}_output.g.vcf.gz
    python3 /usr/local/packages/gatk-4.2.2.0/gatk --java-options "-Xmx4g" HaplotypeCaller -R /local/projects-t3/SerreDLab-3/kko/Pv_P01_Index/PlasmoDB-67_PvivaxP01_Genome.fasta -I $MATE_RG_BAM -L /local/projects-t3/SerreDLab-3/kko/mdr_del_search_dna_whole_genome/gatk_polyclonality/unmasked_regions.intervals -ERC GVCF --dont-use-soft-clipped-bases true -stand-call-conf 20.0 -O $VCF_PATH

    echo "Processing completed for $SAMPLE_NAME."

else
    echo "File $BAM_FILE does not exist."
    exit 1

fi

echo "All processes completed. Aligned files saved in $WORKING_DIR."
