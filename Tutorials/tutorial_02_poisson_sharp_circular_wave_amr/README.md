This tutorial aims to introduce the user on **`FEMPAR`** style driver programs. 
Following the steps of [The Commented Code](#the-commented-code) would help the user to understand 
the steps to solve a Finite Element problem, as the basis for understanding more complex **`FEMPAR`**  data structures.

### Compilation and execution

In order to compile this tutorial you need to have access to **`FEMPAR`** library.
The [Quick start](../../01_Readme.html#quick-start) on the main [README](../../01_Readme.html) is the easiest way to get **`FEMPAR`**.
You can see a more custom and complex installation process in the [Native compilation](../../01_Readme.html#native-compilation) section.

Then, the very first point to start with are the [tutorial programs](https://gitlab.com/fempar/fempar/tree/experimental/Tutorials) available at the official FEMPAR repository. 

After completing the previous compilation steps, you can compile and run FEMPAR tutorials with the following steps:

```bash
...
$ FEMPAR_TUTORIALS_DIR=$SOURCES_DIR/Tutorials
$ mkdir -p $FEMPAR_TUTORIALS_DIR
$ cd $FEMPAR_TUTORIALS_DIR
$ cmake -DFEMPAR_DIR=$FEMPAR_DIR -DFEMPAR_TUTORIAL=tutorial_02_poisson_sharp_circular_wave_amr $SOURCES_DIR/Tutorials
$ make -j 4
$ bin/tutorial_02_poisson_sharp_circular_wave_amr --help                              # get informative message on screen
$ bin/tutorial_02_poisson_sharp_circular_wave_amr [optional command line arguments]   # execute the tutorial
```

### What this program does

This tutorial tackles the Poisson problem. 


#### Model problem

In strong form this problem reads: 

> find ![u](https://latex.codecogs.com/svg.latex?%5Cinline%20u) such that ![strongform](https://latex.codecogs.com/svg.latex?%5Cinline%20-%5CDelta%20u%20%3D%20f%20%5Cqquad%20%5Chbox%7Bin%20%7D%20%5C%2C%20%5COmega%2C) where ![sourceterm](https://latex.codecogs.com/svg.latex?%5Cinline%20f%20%3A%20%5COmega%20%5Crightarrow%20%5Cmathbb%7BR%7D) is a given source term, and ![domain](https://latex.codecogs.com/svg.latex?%5Cinline%20%24%5COmega%3A%3D%5B0%2C1%5D%5Ed%24) is the unit box domain, with ![dimensions](https://latex.codecogs.com/svg.latex?%5Cinline%20d%3A%3D2%2C3) being the number of space dimensions.

Poisson equation problem is supplied with inhomogeneous Dirichlet BCs ![dirichlet](https://latex.codecogs.com/svg.latex?%5Cinline%20u%3Dg) on ![boundary](https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cpartial%20%5COmega), with ![function](https://latex.codecogs.com/svg.latex?%5Cinline%20g%20%3A%20%5Cpartial%20%5COmega%20%5Crightarrow%20%5Cmathbb%7BR%7D) a given function defined on the domain boundary.

@note
Other BCs, e.g., Neumann or Robin (mixed) conditions can also be considered for the Poisson problem. While these sort of BCs are supported by **FEMPAR** as well, we do not consider them in this tutorial for simplicity.
@endnote

The source term ![f](https://latex.codecogs.com/svg.latex?%5Cinline%20f) and Dirichlet function ![g](https://latex.codecogs.com/svg.latex?%5Cinline%20g) are chosen such that the exact (manufactured) solution of poisson equation is:

![solution](https://latex.codecogs.com/svg.latex?%5Cinline%20u%28x%29%20%3A%3D%20%5Cmathrm%7Barctan%7D%28%5Calpha%28%5Csqrt%7B%28x-x_c%29%5Ccdot%28x-x_c%29%7D-r%29%29)

This solution has a sharp circular/spherical wave front of radius ![radius](https://latex.codecogs.com/svg.latex?%5Cinline%20r) centered at ![x_c](https://latex.codecogs.com/svg.latex?%5Cinline%20x_c).

Figures below illustrate the solution for ![dimensions](https://latex.codecogs.com/svg.latex?%5Cinline%20d%3A%3D2%2C3), resp., and parameter values
![alpha](https://latex.codecogs.com/svg.latex?%5Cinline%20%5Calpha%3D200), ![r](https://latex.codecogs.com/svg.latex?%5Cinline%20r%3D0.7), and ![x_c](https://latex.codecogs.com/svg.latex?%5Cinline%20x_c%20%3D%20%28-0.05%2C%20-0.05%29), ![x_c](https://latex.codecogs.com/svg.latex?%5Cinline%20x_c%20%3D%20%28-0.05%2C%20-0.05%2C%20-0.05%29) for ![dimensions](https://latex.codecogs.com/svg.latex?%5Cinline%20d%3A%3D2%2C3), resp.

<img src="media/circular_sharp_wave_2d.png" alt="2D benchmark problem." width=45%/>
<img src="media/circular_sharp_wave_3d.png" alt="3D benchmark problem." width=45%/>


#### FE discretization

While [tutorial 01](../tutorial_01_poisson_sharp_circular_wave/index.html) uses a uniform (thus conforming) mesh ![triang](https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cmathcal%7BT%7D_h) of ![Omega](https://latex.codecogs.com/svg.latex?%5Cinline%20%5COmega), this tutorial combines the [two FE formulations presented](../tutorial_01_poisson_sharp_circular_wave/index.html#fe-discretization) with a more clever/efficient domain discretization approach.

Given that the solution of Poisson equation exhibits highly localized features, in particular an internal sharp circular/spherical wave front layer, this tutorial exploits a **`FEMPAR`** triangulation data structure that efficiently supports dynamic ![h](https://latex.codecogs.com/svg.latex?%5Cinline%20h)-adaptivity techniques (a.k.a. Adaptive Mesh Refinement and coarsening (AMR)), i.e., the ability of the mesh to be refined _in the course of the simulation_ in those regions of the domain that present a complex behaviour (e.g., the internal layer in the case of Poisson problem}), and to be coarsened  in those areas where essentially nothing relevant happens (e.g., those areas away from the internal layer).

This tutorial restricts itself to ![h](https://latex.codecogs.com/svg.latex?%5Cinline%20h)-adaptivity techniques with a fixed polynomial order. This is in contrast to ![hp](https://latex.codecogs.com/svg.latex?%5Cinline%20hp)-adaptivity techniques, in which the local FE space polynomial order ![h](https://latex.codecogs.com/svg.latex?%5Cinline%20p) also varies among cells.
In its first public release, the support of ![hp](https://latex.codecogs.com/svg.latex?%5Cinline%20hp)-adaptivity techniques in **`FEMPAR`** is restricted to non-conforming FE formulations.



