### New repo for a sample of Mike Lynch's C++ PopGen simulation code
#### Code written by Michael Lynch, Repo initialized by R. Taylor Raborn

To run this code, adjust the necessary SBATCH lines (including the first and fourth one, replacing my userID and email address with your own, for this use nano or your preferred CLI text editor) in the file `TwoEffects.sh` and type:
`sbatch TwoEffects.sh`

You can follow the progress of the job using the command
`squeue <username>` where <username> refers to your userID on the machine

This will submit the job on Agave.     Of course, if you're not using Agave you'll need to made the adjustments depending on which job scheduler you use, etc.
