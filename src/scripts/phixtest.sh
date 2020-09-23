#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=3
#SBATCH --time=48:00:00
#SBATCH --mem=20GB
#SBATCH --job-name=LIHC
#SBATCH --mail-type=END
#SBATCH --mail-user=jjb509@nyu.edu
#SBATCH --output=slurm_%j.out

ulimit -s 1300000000

export GOMP_STACKSIZE=2000000
# rename the phixer script
# fit the 
echo "compile phixer"
gcc -Wall pphi_bs.c -fopenmp -o phixer.out
echo "running phixer"
./phixer.out test_100_585.txt
# Rscript grant_support_stage_PFI.R lung_cancers bic 0.01 200 25 LUAD
