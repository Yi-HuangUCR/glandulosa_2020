
#!/bin/bash -l

#SBATCH --mem=64
#SBATCH --time=0-02:00:00
#SBATCH --output=catsubsmptoone.stdout
#SBATCH -p short

cd merged_reads
cat *.subsmp.fa > multispecies_singlefile.fa
mv multispecies_singlefile.fa .. 
