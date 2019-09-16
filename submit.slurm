#! /bin/bash -l
#SBATCH -q premium
#SBATCH -N 1
#SBATCH -t 01:00:00
#SBATCH -J srun_test_62_files
#SBATCH -L project,SCRATCH
#SBATCH -C haswell
#SBATCH --mail-user=skohn@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH --array=0-0

srun -n 1 -c 64 job.sh
