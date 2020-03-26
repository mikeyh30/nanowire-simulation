#! /bin/bash
ARRAY_SIZE=$(nwsetup $1 | tr -d -c 0-9)
qsub -t 1-$ARRAY_SIZE /home/ucapmhy/nanowire-simulations/nanowire-simulation/hpc-control/array_script.sh $1
