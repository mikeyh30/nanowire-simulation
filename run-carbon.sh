#!/bin/bash
DATE=$1

nwsetup $DATE
simulate /home/michael/OneDrive/quantum-devices/simulations/nanowire-simulation/data/carbon/$DATE/$DATE.csv $DATE
