#!/bin/bash -l

#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=clustreads.stdout
#SBATCH -p intel

module load cd-hit

cd-hit-est -i A_glandulosa_singlefile.fa -o A_glandulosa_RAD_ref.fa -c 0.95 -n 10 -M 0 -T 16 -B 1

