#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=3
#SBATCH --time=48:00:00
#SBATCH --mem=20GB
#SBATCH --job-name=LUSCT
#SBATCH --mail-type=END
#SBATCH --mail-user=jjb509@nyu.edu
#SBATCH --output=slurmouts/slurm_%j.out

ulimit -s 1300000000
export GOMP_STACKSIZE=2000000
echo "python edit"
python edit_phixer_script.py LUSC tumor
echo "running phixer"
./phixer_LUSC_t.out ../data/input_data/tcga/LUSC/tumor_expression.txt
mv pruned*tumor_expression.txt ../data/results_data/tcga/LUSC/
