#!/bin/bash

# Default parameters
num_meshes=4
num_dim=2

# Parse inputs
i=0
for arg in "$@"
do
  (( i++ ))
  argpos=$((i+1))
  if [ "$arg" == "--num_meshes" ]
  then
    num_meshes=${!argpos}
  elif [ "$arg" == "--num_dim" ]
  then
    num_dim=${!argpos}
  fi
done

echo "Series of $num_meshes meshes of dimension $num_dim"

driver="${HOME}/Workspace/code/fempar/build/FEMPAR/bin/test_poisson_unfitted_serial"
for (( i=1, nx=10; i<=$num_meshes; i++, nx=nx*2 ))
do
  out_file="out_mesh${i}.txt"
  echo "Running mesh id $i, with $nx elements per direction"
  $driver -tt structured -dim $num_dim -nx $nx -ny $nx -nz $nx -in_space .false. -check .false. > $out_file
done

