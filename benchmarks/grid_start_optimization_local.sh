#!/bin/bash

#$ -N local_optimization_duration
#$ -q c_main.q
#$ -e /scr/minos/jango/Software/RNAdesign-toolbox/benchmarks/optimization_duration_output
#$ -o /scr/minos/jango/Software/RNAdesign-toolbox/benchmarks/optimization_duration_output
#$ -t 1-500

#################################### DONT FORGET TO CHANGE ALL FOLDERS AND SETTINGS!!

###$SGE_TASK_ID
DPATH=/scr/minos/jango/Software/RNAdesign-toolbox/benchmarks/optimization_duration_output/error_local.log
SCRIPT=/scr/minos/jango/Software/RNAdesign-toolbox/benchmarks/optimization_duration.py
START=$(($SGE_TASK_ID*1000-1000))
END=$(($SGE_TASK_ID*1000-1))

python $SCRIPT -x $START -y $END --local > $DPATH/local_optimization_duration_$START.out;
