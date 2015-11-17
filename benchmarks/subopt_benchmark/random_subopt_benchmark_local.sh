#!/bin/bash

#$ -N random_subopt_local
#$ -q c_main.q
#$ -e /scr/minos/jango/Software/RNAdesign-toolbox/benchmarks/subopt_benchmark/results/error_global.log
#$ -o /scr/minos/jango/Software/RNAdesign-toolbox/benchmarks/subopt_benchmark/results
#$ -t 3-30

#################################### DONT FORGET TO CHANGE ALL FOLDERS AND SETTINGS!!

###$SGE_TASK_ID
DPATH=/scr/minos/jango/Software/RNAdesign-toolbox/benchmarks/subopt_benchmark/results
SCRIPT=/scr/minos/jango/Software/RNAdesign-toolbox/benchmarks/subopt_benchmark/random_subopt_benchmark.py
LENGTH=$(($SGE_TASK_ID*10))

for STRUCT in {2..10}; do
    python $SCRIPT -l $LENGTH -s $STRUCT -m mutate_local > $DPATH/random_subopt_local_$LENGTH_$STRUCT.out;
done
