#!/bin/bash -l

#SBATCH --time=0-02:00:00
#SBATCH -p short

Rscript code_2_filter_percent_missing_2N.R
