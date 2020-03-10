#!/bin/bash -l

#SBATCH --array=1-146
#SBATCH -o slurm.%A.%a.out
#SBATCH -e slurm.%A.%a.err
#SBATCH --time=0-02:00:00
#SBATCH -p short

cd aligned_reads
FILE=$(ls *sam | awk "NR==$SLURM_ARRAY_TASK_ID")
NEWNAME=$(ls *sam | sed 's/sam/bam/' | awk "NR==$SLURM_ARRAY_TASK_ID")

module load samtools

samtools view -b $FILE > $NEWNAME

