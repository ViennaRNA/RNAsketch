#!/bin/bash

#$ -N random_optimization_global
#$ -q c_main.q
#$ -e /scr/minos/jango/Software/RNAdesign-toolbox/benchmarks/hoehner-multistable-example/results/random_optimization/error_global.log
#$ -o /scr/minos/jango/Software/RNAdesign-toolbox/benchmarks/hoehner-multistable-example/results/random_optimization
#$ -t 1-51

#################################### DONT FORGET TO CHANGE ALL FOLDERS AND SETTINGS!!

###$SGE_TASK_ID
DPATH=/scr/minos/jango/Software/RNAdesign-toolbox/benchmarks/hoehner-multistable-example/results/random_optimization
SCRIPT=/scr/minos/jango/Software/RNAdesign-toolbox/benchmarks/hoehner-multistable-example/random_optimization.py
START=$(($SGE_TASK_ID*500-500))
END=$(($SGE_TASK_ID*500-1))

python $SCRIPT -m mutate_global -x $START -y $END > $DPATH/random_optimization_global_$START.out;
