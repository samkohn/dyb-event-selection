#! /bin/bash -l
#SBATCH -q premium
#SBATCH -N 6 -c 64
#SBATCH -t 01:00:00
#SBATCH -J medium_test_run
#SBATCH -L project,SCRATCH
#SBATCH -C haswell
#SBATCH --image=docker:samkohn/pyroot:2
#SBATCH --mail-user=skohn@lbl.gov
#SBATCH --mail-type=ALL

export PATH=$PATH:/usr/common/tig/taskfarmer/1.5/bin:$(pwd)
export THREADS=320
runcommands.sh tasks_1_100000.txt
