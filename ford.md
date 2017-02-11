Title: Fempar
Project: Fempar
Author: Large Scale Scientific Computing
Autor_description: The Large Scale Scientific Computing group develops novel finite element formulations for solid mechanics and fluid dynamics (turbulent incompressible and compressible flows). It is particularly focused on the scalability of the whole simulation process on the largest supercomputers today. In this sense, it develops novel domain decomposition preconditioners and implementations that are scalable at extreme scales.
Summary: Finite Element Multiphysics PARallel solvers
Date:    January 10, 2017
base_url: http://fempar.gitlab.io/
github: https://gitlab.com/fempar/fempar
website: http://fempar.gitlab.io/
blank-value: 
docmark: <
display: public
search: false
preprocess: true
source: true
graph: true
print_creation_date: false
fpp_extensions: f90
                i90
src_dir: ./Sources/Lib
output_dir: ./docs
exclude: sort.f90
exclude_dir: ./Sources/Lib/Generic
include: ./Sources/Include
         ./Sources/Lib/Generic
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
            FLAP:http://szaghi.github.io/FLAP
            PENF:http://szaghi.github.io/PENF
            XH5For:https://gitlab.com/fempar/XH5For
            FPL:https://gitlab.com/fempar/FPL
md_extensions: markdown.extensions.toc(anchorlink=False)
               markdown.extensions.extra
               markdown.extensions.footnotes


Finite Element Multiphysics PARallel solvers

[![build status](https://gitlab.com/fempar/fempar/badges/experimental/build.svg)](https://gitlab.com/fempar/fempar/commits/experimental)
[![coverage report](https://gitlab.com/fempar/fempar/badges/experimental/coverage.svg)](https://gitlab.com/fempar/fempar/commits/experimental)

## Links

- [Web page](http://fempar.gitlab.io)
- [Wiki](https://gitlab.com/fempar/fempar/wikis/home)
- [Source code documentation](http://fempar.gitlab.io/documentation)
- [Issue tracker](https://gitlab.com/fempar/fempar/issues)
- [Continuous integration dashboard](https://gitlab.com/fempar/fempar/builds)
- [Testing dashboard](http://my.cdash.org/index.php?project=Fempar)
