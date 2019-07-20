#!/bin/bash
#Absolute PATH to the root directory of FEMPAR's Git repository
FEMPAR_SOURCES_ROOT=~/git-repos/fempar

#Absolute PATH to the compiled tutorial_03 program (for performance reasons, highly recommended that it is compiled in RELEASE mode)
TUTORIAL_03_PATH=~/git-repos/build-gnu9-release/TUTORIAL_03/bin/tutorial_03_poisson_sharp_circular_wave

#Load file from FEMPAR's Git repository with all the names of FEMPAR Command-Line-Arguments (CLAs)
. $FEMPAR_SOURCES_ROOT/Sources/Tests/Scripts/fempar_cla_names

#Set the number of OpenMP threads to be forked by Intel MKL PARDISO direct solver
#(this should be adjusted to the particular number of cores in the multi-core CPU computer architecture this script is running on)
export MKL_NUM_THREADS=4

prefix_data_files="uniform_mesh_error_order_"
for order in 1 2 4 8
do
  data_file_name="$prefix_data_files""$order"
  rm -f $data_file_name    
  for i in 10 100 1000
  do       	
    for j in 1 2 3 4 5 6 7 8 9
    do
      let NX=$i*$j
      let NY=$i*$j
      command="$TUTORIAL_03_PATH $struct_hex_mesh_generator_num_dims_flag 2 \
       	       $struct_hex_mesh_generator_num_cells_x_dim_flag $NX $NY \
               $fes_ref_fe_orders_flag $order \
               --FE_FORMULATION CG --ALPHA 200.0 --CIRCLE_RADIUS 0.7 --CIRCLE_CENTER -0.05 -0.05 --WRITE_POSTPROCESS_DATA FALSE"
      echo $command
      eval $command | tee /tmp/data
      ACTUAL_NDOFS=$(cat /tmp/data | grep "NUM_DOFS" | cut -f2 -d: | sed s/" "" "*/""/g)
      ACTUAL_NDOFS_SQRT=$(awk "BEGIN {printf \"%f\n\", sqrt($ACTUAL_NDOFS)}")
      ERROR=$(cat /tmp/data | grep "GLOBAL" | cut -f2 -d: | sed s/" "" "*/""/g)
      echo $ACTUAL_NDOFS_SQRT $ERROR >> $data_file_name
      if [ 1 -eq "$(echo "${ACTUAL_NDOFS_SQRT} > 1000" | bc)" ]
      then  
        continue 3
      fi
    done
  done 
done
rm -f /tmp/data
