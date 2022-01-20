#!/bin/bash
DATE=$1

nwsetup $DATE
simulate /home/michael/OneDrive/quantum-devices/simulations/nanowire-simulation/data/helium/$DATE/$DATE.csv $DATE
