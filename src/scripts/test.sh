#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00
#SBATCH --mem=20GB
#SBATCH --job-name=LIHC
#SBATCH --mail-type=END
#SBATCH --mail-user=jjb509@nyu.edu
#SBATCH --output=slurm_%j.out	

ulimit -s 1300000000
export GOMP_STACKSIZE=2000000
echo "hello"
python u.py
