#!/bin/bash -l

#SBATCH --time=1-00:00:00
#SBATCH --mem=256G
#SBATCH -p highmem

module load samtools

samtools merge multispecies_all_files_merged.bam *bam

samtools sort multispecies__all_files_merged.bam > multispecies_all_sorted.bam

