#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=3
#SBATCH --time=48:00:00
#SBATCH --mem=20GB
#SBATCH --job-name=LIHC
#SBATCH --mail-type=END
#SBATCH --mail-user=jjb509@nyu.edu
#SBATCH --output=slurm_%j.out


# rename the phixer script
# fit the 
python edit_phixer_script.py LIHC
./phixer_LIHC.out ../data/input_data/tcga/LIHC/normal_expression.txt
# Rscript grant_support_stage_PFI.R lung_cancers bic 0.01 200 25 LUAD
