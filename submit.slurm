#! /bin/bash -l
#SBATCH -q premium
#SBATCH -N 21 -c 64
#SBATCH -t 01:00:00
#SBATCH -J taskfarmer_first_2500
#SBATCH -L project,SCRATCH
#SBATCH -C haswell
#SBATCH --mail-user=skohn@lbl.gov
#SBATCH --mail-type=ALL

export PATH=$PATH:/usr/common/tig/taskfarmer/1.5/bin:$(pwd)
export THREADS=1280
runcommands.sh tasks.txt
