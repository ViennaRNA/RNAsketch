#!/bin/bash

#$ -N modena_comp_tab2
#$ -q c_highmem.q
#$ -l h=archer
#$ -l h_vmem=4G
#$ -e /scr/minos/bioinf/results_tab2/error_global.log
#$ -o /scr/minos/bioinf/results_tab2

#################################### DONT FORGET TO CHANGE ALL FOLDERS AND SETTINGS!!

DPATH=/scr/minos/bioinf/results_tab2
IPATH=/home/bioinf/Birgit/rnadesign-toolbox/benchmarks/MODENA/inputData/RNAtabupath
SCRIPT=/home/bioinf/Birgit/rnadesign-toolbox/benchmarks/MODENA/modena_comp_hoehner.py
FILE=$1

python $SCRIPT -f ${IPATH}/${FILE} > $DPATH/modena_comp_${FILE}.out

