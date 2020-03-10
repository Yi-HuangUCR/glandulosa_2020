#!/bin/bash -l

#SBATCH --time=1-00:00:00
#SBATCH --mem=256G
#SBATCH -p highmem

module load samtools

samtools merge A_glandulosa_all_files_merged.bam *bam

samtools sort A_glandulosa_all_files_merged.bam > A_glandulosa_all_sorted.bam

