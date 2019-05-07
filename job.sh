#! /bin/bash -l

DIR=$1
RUNDIR=$2
ARRAY_NUM=$3
NEVENTS=$4

function log {
  echo "$ARRAY_NUM: $@"
}

log $DIR $RUNDIR $ARRAY_NUM $NEVENTS

log "executing cd $RUNDIR"
cd $RUNDIR

FIND_FILE=~mkramer/projscratch/p17b/code/p17b_find/p17b_find
log "executing time python $DIR/fullstack.py -i `$FIND_FILE 72442 $ARRAY_NUM` -n $NEVENTS"
time python $DIR/fullstack.py -i `$FIND_FILE 72442 $ARRAY_NUM` -n $NEVENTS

log "done"
exit 0
