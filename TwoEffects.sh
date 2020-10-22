#! /bin/bash

#SBATCH -A rraborn
#SBATCH -n 8
#SBATCH -t 0-4:00
#SBATCH --mail-user=rtraborn@asu.edu

module load intel/2019.4
icc -o twoeff TwoEffects.cpp
echo "Running the script"
./twoeff

