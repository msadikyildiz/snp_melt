#!/bin/bash
#SBATCH --partition=512GB
#SBATCH --time=14-00:00
#SBATCH --array=0-2
#SBATCH --nodes=1
#SBATCH --job-name 33-snp-melt
#SBATCH -o /project/greencenter/Toprak_lab/shared/TEM1_Combinatorial_Mutagenesis/log/20240625/33-snp-melt-%j-%a.out
#SBATCH -e /project/greencenter/Toprak_lab/shared/TEM1_Combinatorial_Mutagenesis/log/20240625/33-snp-melt-%j-%a.err

# Set up the environment
EXP_ID="20240625"
PROJECT_PATH="/project/greencenter/Toprak_lab/shared/TEM1_Combinatorial_Mutagenesis"
DATA_PATH="$PROJECT_PATH/data/$EXP_ID/alignment/end-to-end"
# FASTA_PATH="$PROJECT_PATH/data/$EXP_ID/reference/CTXM/pBAD33-CTX-M-15_CML.fa"
FASTA_PATH="$PROJECT_PATH/data/$EXP_ID/reference/TEM1/pBR322_TEM1-TwistBio+barcode.primer.fa"
OUTPUT_DIR="$PROJECT_PATH/data/$EXP_ID/snp_melt"
mkdir -p $OUTPUT_DIR
# Specify the snp_melt
SNP_MELT_PATH="$PROJECT_PATH/tools/snp_melt"

source /home2/s414024/.bashrc

BAMS=($DATA_PATH/T*.bam)
# BioHPC allows for 1000 job array maximum, submit sbatch jobs in 1000 increments
# by manually updating the offset to +0, +1000, +2000, etc.
let "TASK_ID = $SLURM_ARRAY_TASK_ID + 0"
SAMPLE_PATH=${BAMS[$TASK_ID]}
SAMPLE_ID=${SAMPLE_PATH##*/}
SAMPLE_ID=${SAMPLE_ID%.bam}

echo "TASKID: $TASK_ID \t SAMPLE_ID: $SAMPLE_ID"

$SNP_MELT_PATH/snp_melt                 \
    -b $SAMPLE_PATH                     \
    -r $FASTA_PATH                      \
    -t 72                               \
    -m 50                              \
    | gzip > $OUTPUT_DIR/$SAMPLE_ID.csv.gz
