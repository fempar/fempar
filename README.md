# Fempar

**Finite Element Multiphysics PARallel solvers**

[![build status](https://gitlab.com/fempar/fempar/badges/experimental/build.svg)](https://gitlab.com/fempar/fempar/commits/experimental)
[![coverage report](https://gitlab.com/fempar/fempar/badges/experimental/coverage.svg)](https://gitlab.com/fempar/fempar/commits/experimental)

**FEMPAR** is a scientific software for the simulation of problems governed by partial differential equations (PDEs). It provides a set of state-of-the-art numerical discretizations, including finite element methods, discontinuous Galerkin methods, XFEM, and spline-based functional spaces. The library was originally designed to efficiently exploit distributed-memory supercomputers and easily handle multiphysics problems. It also provides a set of highly scalable numerical linear algebra solvers based on multilevel domain decomposition for the systems of equations that arise from PDE discretizations. Some applications of FEMPAR include the simulation of metal additive manufacturing processes, superconductor devices, breeding blankets in fusion reactors, and nuclear waste repositories.

## Software design

For those who are interested on the design and rationale behind the software abstractions in FEMPAR, a very through presentation is available at:

Santiago Badia, Alberto F. Martín and Javier Principe. 
FEMPAR: An object-oriented parallel finite element framework. Archives of Computational Methods in Engineering 25, 2 (2018), 195–271. 
[ArXiv link](arxiv.org/abs/1708.01773) [DOI](https://link.springer.com/article/10.1007%2Fs11831-017-9244-1)

## Links

- [Web page](http://www.fempar.org/)
- [Official repository](https://www.gitlab.com/fempar/fempar)
- [Mirror repository](https://www.github.com/fempar/fempar)
- [Testing dashboard](https://cdash.cimne.upc.edu/user.php)

## Quick Start

The quickest and easiest way to start with FEMPAR is using [Docker](https://www.docker.com). [Docker](https://opensource.com/resources/what-docker) is a tool designed to make it easier to create, deploy, and run applications by using containers.

FEMPAR provides a Docker [container with the required environment](https://hub.docker.com/u/fempar) to compile the project source code and to run tutorials and drivers. 

Please, follow the steps below to compile FEMPAR using the docker container:

```
$ docker run -ti fempar/fempar-env:gnu 
$ WORKDIR=/data
$ SOURCES_DIR=$WORKDIR/sources
$ FEMPAR_DIR=$WORKDIR/FEMPAR
$ git clone --recursive https://github.com/fempar/fempar $SOURCES_DIR
$ cd $WORKDIR
$ cmake $SOURCES_DIR/SuperBuild
$ make -j 4
```

If you are new to FEMPAR, the very first point to start with are the [tutorial programs](https://gitlab.com/fempar/fempar/tree/experimental/Tutorials) available at the official FEMPAR repository. After completing the previous compilation steps, you can compile and run FEMPAR tutorials with the following steps:

```
...
$ FEMPAR_TUTORIALS_DIR=$SOURCES_DIR/Tutorials
$ mkdir -p $FEMPAR_TUTORIALS_DIR
$ cd $FEMPAR_TUTORIALS_DIR
$ cmake -DFEMPAR_DIR=$FEMPAR_DIR -DFEMPAR_TUTORIAL=tutorial_01_steady_poisson $SOURCES_DIR/Tutorials
$ make -j 4
```

## Native compilation

**FEMPAR** uses [CMake](https://cmake.org/) as a portable compilation system. 

Native compilation is only recommended for experienced users. It requires to set up (configure, compile, install, etc) in your own all the mandatory (and optional) software dependencies required to deploy FEMPAR in your Desktop/Laptopp or HPC cluster computing environment. This approach is not fully documented yet.

We strongly recommend to use the `configure` script in `$FEMPARDIR/Tools`, where we assume hereafter that `FEMPARDIR` is an environment variable pointing to the path of the root directory of FEMPAR's git repository. Information of the script can be obtained by typing `configure --help`.

We also strongly recommend to use the `module` functionality, in order to easily switch between different compilers and compiler versions, and to automatically have a well-defined environment.
Instructions to set up this functionality can be found [here](https://gitlab.com/fempar/fempar/wikis/how%20to%20setup%20modules%20environment).
For instance, the compilation of **FEMPAR** using the previous two functionalities using a ++GNU++ compiler would read:

```
$ mkdir build
$ cd build
$ module load gcc/8.2.0
$ $FEMPARDIR/Tools/configure -c GNU -s $FEMPARDIR/SuperBuild
$ make
```

Instead, for INTEL compilers, the last part would read:

```
$ module load ifc/19.0.0
$ $FEMPARDIR/Tools/configure -s $FEMPARDIR/SuperBuild
$ make
```
The compiler by default is the **INTEL** compiler. The tests are compiled unless otherwise stated with `--without-tests`.
One driver can be included in the compilation by adding `-d drivername`.

For stubborn users that want not to take profit from the previous tools, a manual way to compile **FEMPAR** under Linux is (with compilation of tests included):
```
$ mkdir build
$ cd build
$ cmake ../fempar/SuperBuild -DFEMPAR_ENABLE_TESTS=ON
$ make
```
assuming that the right environment configuration is in place through, e.g., `.bashrc`. But this approach is prone to error when switching environments. 
Another approach is to directly include in the `cmake` call the right compilation flags and values (see one of the notes below).

In order to configure the compilation of a driver (right after the previous steps):

```
$ cd build
$ cmake . -DFEMPAR_DRIVER=driver_folder_name
$ make
```

In order to compile FEMPAR library and tests if the first block of commands has been executed:

```
$ cd build/FEMPAR
$ make -jP
```
with ```P``` being the number of parallel processes involved in the compilation


**FEMPAR** compiles with GNU Fortran compiler 5.3.0 (and newer versions) and Intel Fortran compiler 16.0.0 (and newer versions).

## Run tests

In order to run the tests, we need the right environment. If we are relying on `module` functionalities, we must be sure that the required modules are load in the terminal in which we want to run the tests.
It seems that the only module that has to be loaded is `mkl`, whereas the path to other dynamic libraries is hard-coded, e.g., `openmpi`. In any case, one can use the initialization suggested in the `module` manual and pre-load some modules in `modulerc`.
Assuming the right environment is in place, to run all tests in fast mode, we just do:

```
$ cd build/FEMPAR
$ ctest -R fast -VV
```

or 

```
$ cd build/FEMPAR
$ ctest -R test_name
```

to run a particular test.

## Run drivers

Given a driver ```driver_name```, to run it (assuming it has been compiled, see above), we do:

```
$ cd build/DRIVERS/driver_name/bin
$ mpirun -np P ./driver_name [options]
```

where ```P``` is the number of MPI processes to be used. Clearly, ```mpirun -np P``` must be eliminated to run serial drivers.

To see the different options and default values we can do

```
$ ./driver_name --help
```

## Testing dashboard (CDash)

This project offers to its users/developers a testing dashboard service. This service is powered by [CDash](https://www.cdash.org/) on a 
server hosted by CIMNE. Click [here](https://cdash.cimne.upc.edu/user.php) in order to access to the web interface of the service. 
The CDash server gathers and displays rich information regarding the execution of tests which are performed each
time you push into a branch of the fempar repository, and thus, lets you know, e.g., which tests failed,
with which compiler, amount of code, and code covered by the tests, memory defects (e.g., leaks), etc.
If you want to access to this (highly recommended) service, then you have to follow the instructions available
[here](https://cdash.cimne.upc.edu/user.php). Once you are provided with a new user account, then [e-mail us](mailto:amartin@cimne.upc.edu),
so that we can add you to the fempar project at the CDash server. 

## Known issues

**NOTE**: we have detected that some tests (e.g., `test_poisson_unffited`) do **NOT** pass with GNU Fortran compiler 5.5.0, 6.3.0, & 7.3.0 for `experimental` 
commit f7b4199e  due to what it seems to be a compiler BUG. See issue #259 for more details.
Please also note that we do not actually know since 
which commit in `experimental` this is happening, but only that it happens at this one. Thus, avoid using these GNU Fortran compiler versions.
We neither know whether this also happens for GNU compiler version different from the ones above. It does NOT happen with 5.4.0.
**UPDATE**: As pointed out by @principe, the issue disappears with `gfortran` v8.1.0. This version can easily be installed in Ubuntu as:
```
$ sudo apt-get install gfortran-8
$ sudo apt-get install g++-8
```
If you do not currently use or do not plan to install the modules environment (see https://gitlab.com/fempar/fempar/wikis/how%20to%20setup%20modules%20environment for details), 
then one can use this compiler version when configuring FEMPAR as:

```
cmake -DCMAKE_Fortran_COMPILER=gfortran-8 -DCMAKE_C_COMPILER=gcc-8 \
-DCMAKE_CXX_COMPILER=g++-8 -DFEMPAR_ENABLE_TESTS=ON -DFEMPAR_ENABLE_OPENMP=OFF \
-DFORTRAN_EXTRA_FLAGS= -DC_EXTRA_FLAGS= -DCMAKE_BUILD_TYPE=Debug -DMPIEXEC_PREFLAGS= \
-DFEMPAR_ENABLE_BLAS=ON -DFEMPAR_ENABLE_LAPACK=ON -DFEMPAR_ENABLE_MKL=ON \
-DFEMPAR_ENABLE_P4EST=ON -DFEMPAR_ENABLE_UMFPACK=ON -DFEMPAR_ENABLE_METIS=ON \
-DFEMPAR_ENABLE_GIDPOST=OFF ../fempar/SuperBuild/
```

**NOTE**: there is also an open issue with gfortran 6.4.1 and gfortran 7.3.1 (https://gitlab.com/fempar/XH5For/issues/7). An internal
compiler error raises when compiling FoX, a third party library of XH5For.  **UPDATE**: This issue has been already by-passed in FEMPAR's experimental
branch from commit 47b947f2. See https://gitlab.com/fempar/XH5For/issues/7 for more details

**NOTE**: we also detected a BUG with Intel Fortran compiler 18.0.0 related to missing initialization of member variables to default values in the case of 
polymorphic allocatable variables. For example, `test_transient_poisson` do not pass with Intel Fortran compiler 18.0.0 for commit 5176d2976659c64f45e35022bfea5dcb1e72045e.
due to a compiler BUG (see issue #250). Thus, avoid using this Intel Fortran compiler version. With Intel compiler 18.0.1 this issue disappears

**NOTE**: if you plan to use `Intel Parallel Studio XE 2019` in order to compile FEMPAR with the Intel compilers in your machine (this is indeed the only version currently supported 
by `Ubuntu 18.04`), please note the following. The most annoying issue is related to the compilation of `SISL`. If you use the `icc` compiler, the `icpc` C++ compiler must be used as well.
This is already achieved by fempar's `configure` script in `Tools` whenever you specify `-c Intel`, but it won't be if you call `cmake` directly.
In the latter case, you must specify `-DCMAKE_CXX_COMPILER=icpc` explicitly when invoking cmake. 


