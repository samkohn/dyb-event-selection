#! /bin/bash -l
#SBATCH -q debug
#SBATCH -N 3 -c 64
#SBATCH -t 0:30:00
#SBATCH -J basic_test_run
#SBATCH -L project,SCRATCH
#SBATCH -C haswell
#SBATCH --image=docker:samkohn/pyroot:2
#SBATCH --mail-user=skohn@lbl.gov
#SBATCH --mail-type=ALL

export PATH=$PATH:/usr/common/tig/taskfarmer/1.5/bin:$(pwd)
export THREADS=128
runcommands.sh tasks_nh_thu_1000.txt
