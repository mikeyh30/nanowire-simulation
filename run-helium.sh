#!/bin/bash
DATE=$1

nwsetup $DATE
simulate /home/michael/Documents/nanowire-simulation/data/helium/$DATE/$DATE.csv $DATE
