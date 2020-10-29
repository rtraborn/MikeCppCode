#! /bin/bash

#SBATCH -A rraborn
#SBATCH -p cmecpu1
#SBATCH -q cmeqos
#SBATCH -n 19
#SBATCH -t 0-4:00

module load intel/2019.4

echo "Running the script in parallel"
./twosites $SLURM_ARRAY_TASK_ID
