#! /bin/bash -l

SRC_DIR=$1
OUT_DIR=$2
RUNNO=$3
FILENO=$4
NEVENTS=$5
SUBLIST=$6
SITE=$7
RUNFILE_LINENO=$8

function log {
    echo "${RUNFILE_LINENO}.${RUNNO}.${FILENO}: $@"
}

log $SRC_DIR $OUT_DIR $FILENO $NEVENTS

log "creating progress file $SRC_DIR/progress/__in_progress_${RUNFILE_LINENO}__"
touch $SRC_DIR/progress/__in_progress_${RUNFILE_LINENO}__

log "executing cd $OUT_DIR"
cd $OUT_DIR

FIND_FILE=~mkramer/projscratch/p17b/code/p17b_find/p17b_find
log "executing time python $SRC_DIR/fullstack.py -i `$FIND_FILE $RUNNO $FILENO` -n $NEVENTS --site $SITE --lineno $RUNFILE_LINENO"
time python $SRC_DIR/fullstack.py -i `$FIND_FILE $RUNNO $FILENO` -n $NEVENTS --site $SITE --lineno $RUNFILE_LINENO

log "deleting progress file $SRC_DIR/progress/__in_progress_${RUNFILE_LINENO}__"
rm $SRC_DIR/progress/__in_progress_${RUNFILE_LINENO}__
log "deleted progress file $SRC_DIR/progress/__in_progress_${RUNFILE_LINENO}__"
log "done"
NEXT_JOB=`python $SRC_DIR/next_job.py --sublist $SUBLIST --srcdir $SRC_DIR --outdir $OUT_DIR`
log $NEXT_JOB
####$NEXT_JOB
exit 0
