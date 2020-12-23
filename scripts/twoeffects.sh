#!/bin/bash

#SBATCH -A mlynch11
#SBATCH -p cmecpu1
#SBATCH -q cmeqos
#SBATCH -n 1
#SBATCH -t 0-4:00

module load intel/2020.2

echo "Running the script"
./twoeff $SLURM_ARRAY_TASK_ID
