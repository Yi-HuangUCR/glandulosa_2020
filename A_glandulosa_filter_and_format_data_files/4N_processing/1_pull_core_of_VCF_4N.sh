#!/bin/bash -l

#SBATCH -p short

grep "#" A_glandulosa_4N.vcf | tail --l 1 | cut -c 2- | tr '-' '_' > header_delete.txt

grep "TYPE=snp" A_glandulosa_4N.vcf | grep -v "mnp" | grep -v "del" | grep -v "ins" | grep -v "complex" > snps_delete.txt

cat header_delete.txt snps_delete.txt > output_from_step_one_4N.txt

rm header_delete.txt
rm snps_delete.txt


