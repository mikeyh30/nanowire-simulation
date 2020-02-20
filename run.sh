#!/bin/bash
DATE=$1

python nwsetup $DATE
python simulate /scratch/mhynes/nanowire-simulation/$DATE/$DATE.csv $DATE
