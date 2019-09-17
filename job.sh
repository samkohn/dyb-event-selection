#! /bin/bash -l

module load python
module load root

cd $SCRATCH/dyb11
python "$@"
