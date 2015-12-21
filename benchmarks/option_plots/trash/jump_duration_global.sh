#!/bin/bash

#$ -N jump_duration_global
#$ -q c_main.q
#$ -e /scr/minos/jango/Software/RNAdesign-toolbox/benchmarks/hoehner-multistable-example/results/jump_duration/error_global.log
#$ -o /scr/minos/jango/Software/RNAdesign-toolbox/benchmarks/hoehner-multistable-example/results/jump_duration
#$ -t 1-51

#################################### DONT FORGET TO CHANGE ALL FOLDERS AND SETTINGS!!

###$SGE_TASK_ID
DPATH=/scr/minos/jango/Software/RNAdesign-toolbox/benchmarks/hoehner-multistable-example/results/jump_duration
SCRIPT=/scr/minos/jango/Software/RNAdesign-toolbox/benchmarks/hoehner-multistable-example/jump_duration.py
START=$(($SGE_TASK_ID*60-60))
END=$(($SGE_TASK_ID*60-1))

python $SCRIPT -x $START -y $END > $DPATH/jump_duration_global_$START.out;
