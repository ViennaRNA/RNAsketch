#!/bin/bash

#$ -N modena_comp_str4
#$ -q c_highmem.q
#$ -l h=archer
#$ -l h_vmem=4G
#$ -e /scr/minos/jango/Software/RNAdesign-toolbox/benchmarks/MODENA/results_4str/error_global.log
#$ -o /scr/minos/jango/Software/RNAdesign-toolbox/benchmarks/MODENA/results_4str

#################################### DONT FORGET TO CHANGE ALL FOLDERS AND SETTINGS!!

DPATH=/scr/minos/jango/Software/RNAdesign-toolbox/benchmarks/MODENA/results_4str
IPATH=/scr/minos/jango/Software/RNAdesign-toolbox/benchmarks/MODENA/inputData/RNAdesign/4str
SCRIPT=/scr/minos/jango/Software/RNAdesign-toolbox/benchmarks/MODENA/modena_comp_hoehner.py
FILE=$1

python $SCRIPT -f ${IPATH}/${FILE} > $DPATH/modena_comp_${FILE}.out

