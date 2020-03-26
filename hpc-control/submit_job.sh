#! /bin/bash
ARRAY_SIZE=$(nwsetup test | tr -d -c 0-9)
qsub -t 1-$ARRAY_SIZE array_script.sh $1
