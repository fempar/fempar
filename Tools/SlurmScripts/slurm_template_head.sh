#!/bin/bash
#
#SBATCH --qos=class_a
#SBATCH --time=HH:MM:00
#SBATCH --nodes=NNO
#SBATCH --tasks-per-node=TPN
#SBATCH --job-name=NAME
#SBATCH --error=%j_job-stderr.txt 
#SBATCH --output=%j_job-stdout.txt
#SBATCH -x
#
# Template script controlling execution environment. The following variables must be replaced
# (etiher manually or by a wraping script):
#
# ODD: an extra node with only one task is reserved.
# NNO: number of nodes 
# TPN: Tasks per node 
# TNT: Total number of tasks
# CPT: cores per task (= 48/TPN in MN4 if all the cores are to be used)
# NRP: number of repetitions of the execution
#
# Configure IMPI
export I_MPI_PIN=on
export I_MPI_DEBUG=0
#export LD_PRELOAD=/apps/INTEL/2018.0.128/itac/2018.0.015/slib/libvtunwind.so
#export LD_PRELOAD=/apps/INTEL/2018.0.128/itac/2018.0.015/slib/libdwarf.so
#export LD_LIBRARY_PATH=/apps/INTEL/2018.0.128/itac/2018.0.015/slib:$LD_LIBRARY_PATH
#export LD_PRELOAD=/apps/INTEL/2018.0.128/itac/2018.0.015/slib/libVT.so
#export LD_PRELOAD=/apps/INTEL/2017.4/itac/2017.3.030/intel64/slib/libVTmc.so
#
#
############################
