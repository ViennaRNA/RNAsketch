#!/bin/bash

#$ -N modena_comp_str4
#$ -q c_highmem.q
#$ -l h=archer
#$ -l h_vmem=4G
#$ -e /scr/minos/bioinf/results_4str/error_global.log
#$ -o /scr/minos/bioinf/results_4str

#################################### DONT FORGET TO CHANGE ALL FOLDERS AND SETTINGS!!

DPATH=/scr/minos/bioinf/results_4str
IPATH=/home/bioinf/Birgit/rnadesign-toolbox/benchmarks/MODENA/inputData/RNAdesign/4str
SCRIPT=/home/bioinf/Birgit/rnadesign-toolbox/benchmarks/MODENA/modena_comp_hoehner.py
FILE=$1

python $SCRIPT -f ${IPATH}/${FILE} > $DPATH/modena_comp_${FILE}.out

