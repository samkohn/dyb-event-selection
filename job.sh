#! /bin/bash -l

SRC_DIR=$1
OUT_DIR=$2
RUNNO=$3
FILENO=$4
NEVENTS=$5
SUBLIST=$6
SITE=$7
LINENO=$8

function log {
  echo "$FILENO: $@"
}

log $SRC_DIR $OUT_DIR $FILENO $NEVENTS

log "creating progress file $SRC_DIR/progress/__in_progress_${LINENO}__"
touch $SRC_DIR/progress/__in_progress_${LINENO}__

log "executing cd $OUT_DIR"
cd $OUT_DIR

FIND_FILE=~mkramer/projscratch/p17b/code/p17b_find/p17b_find
log "executing time python $SRC_DIR/fullstack.py -i `$FIND_FILE $RUNNO $FILENO` -n $NEVENTS"
time python $SRC_DIR/fullstack.py -i `$FIND_FILE $RUNNO $FILENO` -n $NEVENTS --site $SITE

log "deleting progress file $SRC_DIR/progress/__in_progress_${LINENO}__"
rm $SRC_DIR/progress/__in_progress_${LINENO}__
log "deleted progress file $SRC_DIR/progress/__in_progress_${LINENO}__"
log "done"
NEXT_JOB=`python $SRC_DIR/next_job.py --sublist $SUBLIST --srcdir $SRC_DIR --outdir $OUT_DIR`
log $NEXT_JOB
$NEXT_JOB
exit 0
