#! /bin/bash -l

module load python
module load root

export SRCDIR="/global/homes/s/skohn/dyb-event-selection-production/"
export OUTDIR=$SCRATCH/dyb11
export NWORKERS=60

python $SRCDIR/job-producer.py $SRCDIR/run_list.txt \
    --task $SLURM_ARRAY_TASK_ID \
    --npertask 60 \
    --jobscript $SRCDIR/fullstack.py \
    --nworkers $NWORKERS &

for i in `seq 1 $NWORKERS`
do
    python $SRCDIR/worker.py --outdir $OUTDIR &
done
wait
