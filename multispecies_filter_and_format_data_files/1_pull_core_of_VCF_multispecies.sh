#!/bin/bash -l

#SBATCH -p short

grep "#" multispecies.vcf | tail --l 1 | cut -c 2- | tr '-' '_' > header_delete.txt

grep "TYPE=snp" multispecies.vcf | grep -v "mnp" | grep -v "del" | grep -v "ins" | grep -v "complex" > snps_delete.txt

cat header_delete.txt snps_delete.txt > output_from_step_one_multispecies.delete.txt

cat output_from_step_one_2N.delete.txt | sed 's@forward_reads/@@g' > output_from_step_one_multispecies.txt

rm output_from_step_one_multispecies.delete.txt
rm header_delete.txt
rm snps_delete.txt


