#!/bin/bash
DATE=$1

nwsetup $DATE
simulate /home/michael/Documents/UCL/hydrogen/scratch/$DATE/$DATE.csv $DATE
