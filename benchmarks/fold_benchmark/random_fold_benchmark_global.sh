#!/bin/bash

#$ -N random_fold_global
#$ -q c_main.q
#$ -l h_vmem=2G
#$ -e /scr/minos/jango/Software/RNAdesign-toolbox/benchmarks/fold_benchmark/results/error_global.log
#$ -o /scr/minos/jango/Software/RNAdesign-toolbox/benchmarks/fold_benchmark/results/
#$ -t 2-10

#################################### DONT FORGET TO CHANGE ALL FOLDERS AND SETTINGS!!

###$SGE_TASK_ID
DPATH=/scr/minos/jango/Software/RNAdesign-toolbox/benchmarks/fold_benchmark/results
SCRIPT=/scr/minos/jango/Software/RNAdesign-toolbox/benchmarks/fold_benchmark/random_fold_benchmark.py
LENGTH=$(($SGE_TASK_ID*10))
STRUCT=$1

if [ ! -f $DPATH/random_fold_global_${LENGTH}_${STRUCT}.out ]
then
    python $SCRIPT -l $LENGTH -s $STRUCT -m sample_global > $DPATH/random_fold_global_${LENGTH}_${STRUCT}.out;
fi
