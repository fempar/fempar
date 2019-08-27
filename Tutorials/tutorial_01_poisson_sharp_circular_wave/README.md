This tutorial aims to introduce the user on **`FEMPAR`** style driver programs. 
Following the steps of [The Commented Code](#the-commented-code) would help the user to understand 
the steps to solve a simple Finite Element problem, as the basis for understanding more complex **`FEMPAR`**  data structures.

### Compilation and execution
In order to compile this tutorial, first compile **`FEMPAR`**, see the [README](../../01_Readme.html#compilation).

Then, move to the  **`FEMPAR`** build directory, create the tutorials build subfolder, and compile the tutorial.
```
cd build
mkdir TUTORIALS
cd TUTORIALS
cmake ../../fempar/Tutorials -DFEMPAR_DIR=../FEMPAR -DFEMPAR_TUTORIAL=tutorial_01_poisson_sharp_circular_wave
make
```

Once the tutorial is compiled, you can run it and play around with it:
```
bin/tutorial_01_poisson_sharp_circular_wave
```


### What this program does

This tutorial tackles the Poisson problem. In strong form this problem reads: find ![u](https://latex.codecogs.com/svg.latex?%5Cinline%20u) such that

 ![strongform](https://latex.codecogs.com/svg.latex?%5Cinline%20-%5CDelta%20u%20%3D%20f%20%5Cqquad%20%5Chbox%7Bin%20%7D%20%5C%2C%20%5COmega%2C) 

where ![sourceterm](https://latex.codecogs.com/svg.latex?%5Cinline%20f%20%3A%20%5COmega%20%5Crightarrow%20%5Cmathbb%7BR%7D) is a given source term, and ![domain](https://latex.codecogs.com/svg.latex?%5Cinline%20%24%5COmega%3A%3D%5B0%2C1%5D%5Ed%24) is the unit box domain, with ![dimensions](https://latex.codecogs.com/svg.latex?%5Cinline%20d%3A%3D2%2C3) being the number of space dimensions.
Poisson equation problem is supplied with inhomogeneous Dirichlet BCs ![dirichlet](https://latex.codecogs.com/svg.latex?%5Cinline%20u%3Dg) on ![boundary](https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cpartial%20%5COmega), with ![function](https://latex.codecogs.com/svg.latex?%5Cinline%20g%20%3A%20%5Cpartial%20%5COmega%20%5Crightarrow%20%5Cmathbb%7BR%7D) a given function defined on the domain boundary.

@note
Other BCs, e.g., Neumann or Robin (mixed) conditions can also be considered for the Poisson problem. While these sort of BCs are supported by **FEMPAR** as well, we do not consider them in this tutorial for simplicity.
@endnote

The source term ![f](https://latex.codecogs.com/svg.latex?%5Cinline%20f) and Dirichlet function ![g](https://latex.codecogs.com/svg.latex?%5Cinline%20g) are chosen such that the exact (manufactured) solution of poisson equation is:

![solution](https://latex.codecogs.com/svg.latex?%5Cinline%20u%28x%29%20%3A%3D%20%5Cmathrm%7Barctan%7D%28%5Calpha%28%5Csqrt%7B%28x-x_c%29%5Ccdot%28x-x_c%29%7D-r%29%29)

This solution has a sharp circular/spherical wave front of radius ![radius](https://latex.codecogs.com/svg.latex?%5Cinline%20r) centered at ![x_c](https://latex.codecogs.com/svg.latex?%5Cinline%20x_c).

Figure 1 and Figure 2 illustrate the solution for ![dimensions](https://latex.codecogs.com/svg.latex?%5Cinline%20d%3A%3D2%2C3), resp., and parameter values
![alpha](https://latex.codecogs.com/svg.latex?%5Cinline%20%5Calpha%3D200), ![r](https://latex.codecogs.com/svg.latex?%5Cinline%20r%3D0.7), and ![x_c](https://latex.codecogs.com/svg.latex?%5Cinline%20x_c%20%3D%20%28-0.05%2C%20-0.05%29), ![x_c](https://latex.codecogs.com/svg.latex?%5Cinline%20x_c%20%3D%20%28-0.05%2C%20-0.05%2C%20-0.05%29) for ![dimensions](https://latex.codecogs.com/svg.latex?%5Cinline%20d%3A%3D2%2C3), resp.

<img src="media/circular_sharp_wave_2d.png" alt="2D benchmark problem." width=45%/>
<img src="media/circular_sharp_wave_3d.png" alt="3D benchmark problem." width=45%/>




