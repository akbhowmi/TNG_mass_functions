#!/bin/bash
#SBATCH --export=NONE
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00
#SBATCH --job-name=test
#SBATCH --output=bh_MF.stdout
#SBATCH --error=bh_MF.stderr
#SBATCH --mem-per-cpu=1000mb
#SBATCH --time=48:00:00
#SBATCH --mail-user=akbhowmi@andrew.cmu.edu
#SBATCH --mail-type=END,FAIL

python generate_SM_HM_relation.py
