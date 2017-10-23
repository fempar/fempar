# Fempar

Finite Element Multiphysics PARallel solvers

[![build status](https://gitlab.com/fempar/fempar/badges/experimental/build.svg)](https://gitlab.com/fempar/fempar/commits/experimental)
[![coverage report](https://gitlab.com/fempar/fempar/badges/experimental/coverage.svg)](https://gitlab.com/fempar/fempar/commits/experimental)

## Links

- [Web page](http://www.fempar.org/)
- [Wiki](https://gitlab.com/fempar/fempar/wikis/home)
- [Source code documentation](http://fempar.org/documentation/)
- [Issue tracker](https://gitlab.com/fempar/fempar/issues)
- [Continuous integration dashboard](https://gitlab.com/fempar/fempar/builds)
- [Testing dashboard](http://my.cdash.org/index.php?project=Fempar)

## Compilation

**FEMPAR** compiles with GNU Fortran compiler 5.3.0 (and newer versions) and Intel Fortran compiler 16.0.0 (and newer versions).

**FEMPAR** uses [CMake](https://cmake.org/) as a portable compilation system. 

The easiest way to compile **FEMPAR** under Linux is (with compilation of tests included):

```
$ mkdir build
$ cd build
$ cmake ../fempar/SuperBuild -DFEMPAR_ENABLE_TESTS=ON
$ make
```

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

In order to run the all tests in fast mode:

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












