### New repo for a sample of Mike Lynch's C++ PopGen simulation code
#### Code written by Michael Lynch, Repo initialized and hosted by R. Taylor Raborn. Gil Speyer and Rebecca Belshe at ASU Research Computing have contributed updates and improvements to the code.

All of the code is found in the `scripts/` directory.
To naviate there, type the following:
`cd scripts`

You have the choice to run the code in serial or in parallel.

To run `TwoSites`, in serial, adjust the necessary SBATCH lines (i.e. `#SBATCH -t 0-4:00` for time) in the file `twosites.sh` and type:
`sbatch twosites.sh`

To run `TwoSites`, in parallel, adjust the necessary SBATCH lines (i.e. `#SBATCH -t 0-4:00` for time) in the file `twosites.sh` and type:
`sbatch --array=1-19 twosites_par.sh`

You can follow the progress of the job using the command
`squeue <username>` where <username> refers to your userID on the machine

This will submit the job on Agave.
Of course, if you're not using Agave you'll need to made the adjustments depending on which job scheduler you use, etc.
