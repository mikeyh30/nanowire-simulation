#!/bin/bash
DATE=$1

python src/setup-sim.py $DATE
python ./src/simulation_all.py ./Scratch/$DATE/$DATE.csv $DATE