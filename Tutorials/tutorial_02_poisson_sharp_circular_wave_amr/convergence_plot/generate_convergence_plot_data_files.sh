#!/bin/bash
#Absolute PATH to the root directory of FEMPAR's Git repository
FEMPAR_SOURCES_ROOT=~/git-repos/fempar

#Absolute PATH to the compiled tutorial_02 program (for performance reasons, highly recommended that it is compiled in RELEASE mode)
TUTORIAL_02_PATH=~/git-repos/build-gnu9-release/TUTORIAL_02/bin/tutorial_02_poisson_sharp_circular_wave_amr

#Load file from FEMPAR's Git repository with all the names of FEMPAR Command-Line-Arguments (CLAs)
. $FEMPAR_SOURCES_ROOT/Sources/Tests/Scripts/fempar_cla_names

#Set the number of OpenMP threads to be forked by Intel MKL PARDISO direct solver
#(this should be adjusted to the particular number of cores in the multi-core CPU computer architecture this script is running on)
export MKL_NUM_THREADS=4

declare -A num_amr_steps
#Number of AMR steps required to achieve a ~1.0e+06 DOFs problem or the error is below 1.0e-06, assuming 0.1 and 0.05 as refinement and coarsening fractions, resp.,
#for fixed fraction refinement strategy
num_amr_steps[1]=30
num_amr_steps[2]=25
num_amr_steps[4]=19
num_amr_steps[8]=10

prefix_data_files="amr_error_order_"
for order in 1 2 4 8
do
  data_file_name="$prefix_data_files""$order"
  rm -f $data_file_name    
  command="$TUTORIAL_02_PATH $p4est_triang_num_dims_flag 2 --NUM_UNIFORM_REFINEMENT_STEPS 4 --NUM_AMR_STEPS ${num_amr_steps[$order]} \
           $fes_ref_fe_orders_flag $order --FE_FORMULATION CG --ALPHA 200.0 --CIRCLE_RADIUS 0.7 --CIRCLE_CENTER -0.05 -0.05 --WRITE_POSTPROCESS_DATA FALSE"
  echo $command
  eval $command | tee /tmp/data
  NDOFS=$(cat /tmp/data | grep NUM_DOFS | cut -f2 -d: |  sed s/" "" "*/""/g)
  ERRORS=$(cat /tmp/data | grep "GLOBAL" | cut -f2 -d: | sed s/" "" "*/""/g)
  j=1
  for i in $NDOFS
  do
    A=$i
    ERROR=$(echo $ERRORS|cut -f$j -d" ")
    ACTUAL_NDOFS_SQRT=$(awk "BEGIN {printf \"%f\n\", sqrt($A)}")
    echo $ACTUAL_NDOFS_SQRT $ERROR >> $data_file_name
    let j=j+1
  done
done
rm -f /tmp/data
