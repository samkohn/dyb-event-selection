#! /bin/bash -l
#SBATCH -q regular
#SBATCH -N 1
#SBATCH -t 01:00:00
#SBATCH -J skohn-dyb-test
#SBATCH -L project,SCRATCH
#SBATCH -C haswell
#SBATCH --mail-user=skohn@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH --array=1-13

export DIR="/global/homes/s/skohn/dyb-event-selection-production/"
export RUNDIR=$SCRATCH/dyb6
SUBLIST=$SLURM_ARRAY_TASK_ID
module load root

#for i in `seq 1 64`
while IFS= read -r line
do
    #NEXT_JOB=`python $DIR/next_job.py --sublist $SUBLIST --srcdir $DIR --outdir $RUNDIR`
    NEXT_JOB=`python $DIR/next_job.py --lineno $line --srcdir $DIR --outdir $RUNDIR`
    echo $NEXT_JOB
    $NEXT_JOB &
    sleep 3
done < $DIR/progress/notdone-${SLURM_ARRAY_TASK_ID}.txt

until [ $(jobs | wc -l) -eq 0 ]
do
    echo "BASE: waiting for `jobs | wc -l` background jobs to finish:"
    jobs
    sleep 60
done
