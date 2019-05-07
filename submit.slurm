#! /bin/bash -l
#SBATCH -q regular
#SBATCH -N 1
#SBATCH -t 01:00:00
#SBATCH -J skohn-dyb-test
#SBATCH -L project,SCRATCH
#SBATCH -C haswell
#SBATCH --mail-user=skohn@lbl.gov
#SBATCH --mail-type=ALL
####### IGNORE SBATCH --array=0-14

export DIR="/global/homes/s/skohn/dyb-event-selection/"
export RUNDIR=$SCRATCH/dyb4
module load root

NEVENTS=-1
for i in `seq 1 64`
do
    $DIR/job.sh $DIR $RUNDIR $i $NEVENTS &
done

until [ $(jobs | wc -l) -eq 0 ]
do
    echo "BASE: waiting for `jobs | wc -l` background jobs to finish:"
    jobs
    sleep 60
done
