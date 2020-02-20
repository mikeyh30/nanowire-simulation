#!/bin/bash
DATE=$1

nwsetup $DATE
simulate /scratch/mhynes/nanowire-simulation/$DATE/$DATE.csv $DATE
