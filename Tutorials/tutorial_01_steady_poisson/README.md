This tutorial aims to introduce the user on **`FEMPAR`** style driver programs. 
Following the steps of [The Commented Code](#the-commented-code) would help the user to understand 
the steps to solve a simple Finite Element problem, as the basis for understanding more complex **`FEMPAR`**  data structures.

### Compilation and execution
In order to compile this tutorial, first compile **`FEMPAR`**, see the [README](../../README.html#compilation).

Then, move to the  **`FEMPAR`** build directory, create the tutorials build subfolder, and compile the tutorial.
```
cd build
mkdir TUTORIALS
cd TUTORIALS
cmake ../../fempar/Tutorials -DFEMPAR_DIR=../FEMPAR -DFEMPAR_TUTORIAL=tutorial_01_steady_poisson
make
```

Once the tutorial is compiled, you can run it and play around with it:
```
bin/tutorial_01_steady_poisson
```


### What this program does

This tutorial solves the Poisson Equation using the Finite Element tools implememented in **`FEMPAR`**,

![Poisson equation](https://latex.codecogs.com/svg.latex?%5Cbegin%7Baligned%7D%20-%20%5CDelta%20u%20%26%3D%20f%20%5Cqquad%20%26in%20%5Cquad%20%26%5COmega%20%5C%5C%20u%20%26%3D%20u_D%20%5Cqquad%20%26on%20%5Cquad%20%26%5Cpartial%20%5COmega%20%5Cend%7Baligned%7D)

where the source term, ![Source term](https://latex.codecogs.com/svg.latex?%5Cinline%20%28f%20%3D%200%29), 
and the Dirichlet boundary conditions, ![Dirichlet conditions](https://latex.codecogs.com/svg.latex?%5Cinline%20%28u_D%20%3D%20x%20+%20y%29), 
are chosen such that the the solution is ![Solution](https://latex.codecogs.com/svg.latex?%5Cinline%20%28u%3D%20x%20+%20y%29).

The Poisson equation is solved within the default `FEMPAR` settings and the following FE setup, which can be found within the source code description.

+ **Domain  ![Omega](https://latex.codecogs.com/svg.latex?%5Cinline%20%28%5COmega%29)**: a hexahedron structured mesh of `10x10` in 2D. 
+ **Reference FE**: Lagrangian FE of 1st order
+ **Boundary Conditions**: Dirichled  BC with value of ![Dirichlet value](https://latex.codecogs.com/svg.latex?%5Cinline%20%28u_D%29)
+ **Solver**: Iterative linear solver




