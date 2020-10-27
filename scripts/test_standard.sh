#!/bin/bash
#SBATCH -c 1
#SBATCH -t 00-04:00
#SBATCH -o %j.OUT
#SBATCH -e %j.ERROR

module load intel/2019.4
./a.out 
