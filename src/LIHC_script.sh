#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=3
#SBATCH --time=48:00:00
#SBATCH --mem=20GB
#SBATCH --job-name=LIHCN
#SBATCH --mail-type=END
#SBATCH --mail-user=jjb509@nyu.edu
#SBATCH --output=slurmouts/slurm_%j.out

ulimit -s 1300000000

export GOMP_STACKSIZE=2000000
mkdir ../data/results_data/tcga/LIHC 
echo "python edit"
python edit_phixer_script.py LIHC normal
echo "running phixer"
./phixer_LIHC.out ../data/input_data/tcga/LIHC/normal_expression.txt
mv pruned*normal_expression.txt ../data/results_data/tcga/LIHC/
python edit_phixer_script.py LIHC tumor 
./phixer_LIHC.out ../data/input_data/tcga/LIHC/tumor_expression.txt
mv pruned*tumor_expression.txt ../data/results_data/tcga/LIHC/
