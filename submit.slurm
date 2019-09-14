#! /bin/bash -l
#SBATCH -q premium
#SBATCH -N 1
#SBATCH -t 01:00:00
#SBATCH -J skohn-dyb-test
#SBATCH -L project,SCRATCH
#SBATCH -C haswell
#SBATCH --mail-user=skohn@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH --array=0-0

module load python
module load root

export SRCDIR="/global/homes/s/skohn/dyb-event-selection-production/"
export OUTDIR=$SCRATCH/dyb11
export NWORKERS=62

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
