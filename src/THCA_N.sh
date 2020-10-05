#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=128:00:00
#SBATCH --mem=20GB
#SBATCH --job-name=THCAN
#SBATCH --mail-type=END
#SBATCH --mail-user=jjb509@nyu.edu
#SBATCH --output=slurmouts/slurm_%j.out

ulimit -s 1300000000
export GOMP_STACKSIZE=2000000
python agony_pipeline.py THCA normal