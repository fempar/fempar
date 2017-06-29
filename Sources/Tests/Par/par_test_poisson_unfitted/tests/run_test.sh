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

if [ "$num_dim" == "2" ]
then
  np=5
elif [ "$num_dim" == "3" ]
then
  np=9
fi

echo "Series of $num_meshes meshes of dimension $num_dim"

driver="${HOME}/Workspace/code/fempar/build/FEMPAR/bin/par_test_poisson_unfitted"
for (( i=1, nx=10; i<=$num_meshes; i++, nx=nx*2 ))
do
  out_file="out_mesh${i}.txt"
  echo "Running mesh id $i, with $nx elements per direction"
  mpirun -n $np $driver -tt 1 -l 2 -dm $num_dim -n $nx $nx $nx -np 2 2 2 1 1 1 -in_space .false. -check .false. > $out_file
done

