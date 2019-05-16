#!/bin/bash
NUM_THREADS=$1
NUM_THREADS_LAST=$2
PROC_ID_COARSE=$3
if [ $PMI_RANK -eq $PROC_ID_COARSE ]
then
   export MKL_NUM_THREADS=$NUM_THREADS_LAST
   export OMP_NUM_THREADS=$NUM_THREADS_LAST
   export KMP_AFFINITY=compact,verbose
   eval $4
else
   export MKL_NUM_THREADS=$NUM_THREADS
   export OMP_NUM_THREADS=$NUM_THREADS
   export KMP_AFFINITY=compact
   eval $4
fi
