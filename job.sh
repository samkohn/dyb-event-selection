#! /bin/bash -l

SRC_DIR=$1
OUT_DIR=$2
RUNNO=$3
FILENO=$4
NEVENTS=$5
SUBLIST=$6

function log {
  echo "$FILENO: $@"
}

log $SRC_DIR $OUT_DIR $FILENO $NEVENTS

log "executing cd $OUT_DIR"
cd $OUT_DIR

FIND_FILE=~mkramer/projscratch/p17b/code/p17b_find/p17b_find
log "executing time python $SRC_DIR/fullstack.py -i `$FIND_FILE $RUNNO $FILENO` -n $NEVENTS"
time python $SRC_DIR/fullstack.py -i `$FIND_FILE $RUNNO $FILENO` -n $NEVENTS

log "done"
NEXT_JOB=`python $SRC_DIR/next_job.py --sublist $SUBLIST --srcdir $SRC_DIR --outdir $OUT_DIR`
log $NEXT_JOB
$NEXT_JOB
exit 0
