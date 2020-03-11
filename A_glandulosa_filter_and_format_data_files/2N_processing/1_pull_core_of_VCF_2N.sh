#!/bin/bash -l

#SBATCH -p short

grep "#" A_glandulosa_2N.vcf | tail --l 1 | cut -c 2- | tr '-' '_' > header_delete.txt

grep "TYPE=snp" A_glandulosa_2N.vcf | grep -v "mnp" | grep -v "del" | grep -v "ins" | grep -v "complex" > snps_delete.txt

cat header_delete.txt snps_delete.txt > output_from_step_one_2N.txt

rm header_delete.txt
rm snps_delete.txt


