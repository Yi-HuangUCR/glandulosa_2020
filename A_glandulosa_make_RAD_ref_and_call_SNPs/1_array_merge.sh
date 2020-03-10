#!/bin/bash -l

#SBATCH --array=1-146
#SBATCH -o slurm.%A.%a.out
#SBATCH -e slurm.%A.%a.err
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH -p intel

module load pear

FORWARD=$(ls forward_reads/ | grep ".fq" | grep -v "fai" | awk "NR==$SLURM_ARRAY_TASK_ID")
REVERSE=$(ls reverse_reads/ | grep ".fq" | grep -v "fai" | awk "NR==$SLURM_ARRAY_TASK_ID")
OUTPUT=$(ls forward_reads/ | grep ".fq" | grep -v "fai" | sed 's@\.1\.@.merged.@g' | awk "NR==$SLURM_ARRAY_TASK_ID")

pear -f forward_reads/$FORWARD -r reverse_reads/$REVERSE -o merged_reads/$OUTPUT
