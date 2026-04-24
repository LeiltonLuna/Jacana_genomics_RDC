#!/bin/bash
#SBATCH --account=account_name
#SBATCH --job-name=merge_fastq
#SBATCH --output=merge_fastq_%A_%a.out
#SBATCH --error=merge_fastq_%A_%a.err
#SBATCH --array=0-4 #number of samples to be processed in parallel
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# Define samples
samples=(SAMPLE_1 SAMPLE_2 SAMPLE_3 SAMPLE_4 SAMPLE_N)

sample=${samples[$SLURM_ARRAY_TASK_ID]}

echo "Processing sample: $sample"

# Merge R1 and R2 across lanes
cat ${sample}_*_L00*_R1_001.fastq.gz > ${sample}_R1.fastq.gz
cat ${sample}_*_L00*_R2_001.fastq.gz > ${sample}_R2.fastq.gz

echo "Finished sample: $sample"
