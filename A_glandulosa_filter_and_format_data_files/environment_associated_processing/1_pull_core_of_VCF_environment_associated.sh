#!/bin/bash -l

#SBATCH -p short

grep "#" A_glandulosa_environment_associated.vcf | tail --l 1 | cut -c 2- | tr '-' '_' > header_delete.txt

grep -v "#" A_glandulosa_environment_associated.vcf  | grep -v "mnp" | grep -v "del" | grep -v "ins" | grep -v "complex" > snps_delete.txt

cat header_delete.txt snps_delete.txt > output_from_step_one_environment_associated.txt


