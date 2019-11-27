#!/bin/bash -l
# Batch script to run a serial array job with the upgraded
# software stack under SGE.
# 1. Force bash
#$ -S /bin/bash

# 2. Request wallclock time (format hours:minutes:seconds).
#$ -l h_rt=12:00:00

# 3. Request 1 gigabyte of RAM (must be an integer)
#$ -l mem=1G

# 4. Request 2 gigabyte of TMPDIR space (default is 10 GB)

# 5. Set up the job array.  In this instance we have requested 10000 tasks
# numbered 1 to 10000.
#$ -t 1-7

# 6. Set the name of the job.
#$ -N 'params'

# 7. Set the working directory to somewhere in your scratch space.  This is
# a necessary step with the upgraded software stack as compute nodes cannot
# write to $HOME.
# Replace "<your_UCL_id>" with your UCL user ID :)
#$ -wd /home/ucapmhy/Scratch/output/logs

# 8. Setup email notifications
#$ -M michael.hynes.18@ucl.ac.uk

# 9. Setup virtual environment
#$ -V
# maybe need this? /home/ucapmhy/nanowire-simulations/virtualenv3/bin/activate

# 8. Run the application. Commented version is to determine the required resources.
# /usr/bin/time --verbose python /home/ucapmhy/nanowire-simulations/nanowire-simulation/hpc_simulation.py /home/ucapmhy/nanowire-simulations/nanowire-simulation/2019-11-15.csv $((SGE_TASK_ID-1))
python /home/ucapmhy/nanowire-simulations/nanowire-simulation/hpc_simulation.py /home/ucapmhy/nanowire-simulations/nanowire-simulation/2019-11-15.csv $((SGE_TASK_ID-1))
