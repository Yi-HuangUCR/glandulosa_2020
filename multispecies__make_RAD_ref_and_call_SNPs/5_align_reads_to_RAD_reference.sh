#!/bin/bash -l

#SBATCH --array=1-17
#SBATCH -o slurm.%A.%a.out
#SBATCH -e slurm.%A.%a.err
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH -p intel

module load samtools

samtools faidx A_gland_RAD_ref.fa

# Align reads to the RAD reference file.
module load bwa
cd forward_reads
FORWARD=$(ls *fq | awk "NR==$SLURM_ARRAY_TASK_ID")
cd ..
cd reverse_reads
REVERSE=$(ls *fq | awk "NR==$SLURM_ARRAY_TASK_ID")
cd ..
cd forward_reads
OUTPUT=$(ls *fq | sed 's@\.1.fq\.@.aln_w_rg.sam@g' | awk "NR==$SLURM_ARRAY_TASK_ID")
cd ..

sample=$(echo $FORWARD | cut -d. -f1)

bwa index forward_reads/$FORWARD
bwa index reverse_reads/$REVERSE
bwa mem -M -t 16 multispecies_RAD_ref.fa forward_reads/$FORWARD reverse_reads/$REVERSE -R '@RG\tID:'$sample'\tSM:'$sample'\tLB:'$sample'\' > aligned_reads/$OUTPUT
