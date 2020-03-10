#!/bin/bash -l

#SBATCH --array=1-146
#SBATCH -o slurm.%A.%a.out
#SBATCH -e slurm.%A.%a.err
#SBATCH --mem=64G
#SBATCH --time=0-02:00:00
#SBATCH -p short

cd merged_reads
FILE=$(ls *fastq | grep "merged" | awk "NR==$SLURM_ARRAY_TASK_ID")
NEWFILE=$(ls *fastq | grep "merged" | sed 's/.fq.assembled.fastq/.subsmp.fa/' | awk "NR==$SLURM_ARRAY_TASK_ID")

paste - - - - < $FILE | cut -f 1,2 | sed 's/^@/>/' | shuf -n 200000 | tr "\t" "\n" > $NEWFILE


