
## Downloading and installing FEMPAR and its tutorial programs

The quickest and easiest way to start with **`FEMPAR`** is using Docker. Docker is a tool designed to easily create, deploy, and run applications by using containers. **`FEMPAR`** provides a set of Docker containers with the required environment (serial or parallel, debug or release) to compile the project source code and to run tutorials and tests. A detailed and very simple installation guide can be found in  [FEMPAR README](https://gitlab.com/fempar/fempar/blob/experimental/README.md), together with instructions for the compilation of the tutorial programs explained below.

## Common structure and usage instructions of FEMPAR tutorials

A **`FEMPAR`** tutorial is a FE application program that uses the tools (i.e., Fortran200X derived data types and their TBPs) provided by **`FEMPAR`** in order to approximate the solution of a PDE (or, more generally, a system of such equations). We strive to follow a common structure for all **`FEMPAR`** tutorials in the seek of uniformity and ease of presentation of the different features of the library. Such structure is sketched in code below. Most of the code of tutorial programs is encompassed within a single program unit. Such unit is in turn composed of four parts:

1. import of external module symbols,
2. declaration of tutorial parameter constants and variables.
3. the main executable code of the tutorial, and 
4. implementation of helper procedures within the `contains` section of the program unit.

```fortran
program tutorial_#_...
  use fempar_names
  use tutorial_#_support_ module1_names
  use tutorial_#_support_ module2_names
  ... ! Usage of the rest of support modules
  implicit none
  ... ! Declaration of tutorial_#_... parameter constants
  ... ! Declaration of data type instances required by tutorial_#_ ...
  ... ! Declaration of variables storing CLA values particular to tutorial_#_ ...
  call fempar_init() !(i.e., construct system-wide objects)
  call setup_parameter_handler()
  call get_tutorial_cla_values()
  ... ! Calls to the rest of helper procedures within the contains section
      ! in order to drive all the necessary steps in the FE simulation
  call fempar_finalize() ! (i.e., destroy system-wide objects)
contains
  ... ! Implementation of helper procedures
end program tutorial_ # _ ...
```

### Modules

In part 1), the tutorial uses the `fempar_names` module to import all **`FEMPAR`** library symbols (i.e., derived types, parameter constants, system-wide variables, etc.), and a set of tutorial-specific module units, which are not part of the **`FEMPAR`** library, but developed specifically for the problem at hand. Each of these modules defines a tutorial-specific data type and its TBPs. Although not necessarily, these are typically type extensions (i.e., subclasses) of parent classes defined within **`FEMPAR`**. These data type extensions let the user define problem-specific ingredients such as, e.g., the source term of the PDE, the function to be imposed on the Dirichlet and/or Neumann boundaries, or the definition of a discrete weak form suitable for the problem at hand, while leveraging (re-using) the code within **`FEMPAR`** by means of Fortran200X native support of run-time polymorphism.

### Constants and variables.

In part 2), the tutorial declares a set of parameter constants, typically the tutorial name, authors, problem description, etc. (to be output on screen on demand by the user), the tutorial data type instances in charge of the FE simulation, such as the triangulation (mesh) of the computational domain or the FE space from which the approximate solution of the PDE is sought, and a set of variables to hold the values of the CLAs (_Command Line Arguments_) which are specific to the tutorial. As covered in the sequel in more detail, tutorial users are provided with a CLI (_Command Line Interface_). Such interface constitutes the main communication mechanism to provide the input required by tutorial programs apart from, e.g., mesh data files generated from the GiD unstructured mesh generator (if the application problem requires such kind of meshes).

### Main code 

Part 3) contains the main tutorial executable code, which is in charge of driving all the necessary FE simulation steps. This code in turn relies on part 4), i.e., a set of helper procedures implemented within the `contains` section of the program unit. The main tasks of a FE program (and thus, a **`FEMPAR`** tutorial), even for transient, non-linear PDE problems, typically encompass:

1. to set up a mesh and a FE space;
2. to assemble a discrete linear or linearized algebraic system of equations; 
3. to solve the system built in 2); 
4. to numerically post-process and/or visualize the solution. There is an almost one-to-one correspondence among these tasks and the helper procedures within the `contains` section.

The main executable code of the prototypical **`FEMPAR`** tutorial is (and _must be_) encompassed within calls to `fempar_init()` and `fempar_finalize()`. The former constructs/initializes all system-wide objects, while the latter performs the reverse operation. For example, in the call to `fempar_init()`, a system-wide dictionary of creational methods for iterative linear solver instances is set up. Such dictionary lays at the kernel of a _Creational OO design pattern_  that lets **`FEMPAR`** users to add new iterative linear solver implementations _without the need to recompile the library at all_. Apart from these two calls, the tutorial main executable code also calls the `setup_parameter_handler()` and `get_tutorial_cla_values()` helper procedures in, which are related to CLI processing.

### Parameter handler

The code of the `setup_parameter_handler()` helper procedure is shown in the code snippet below. It sets up the so-called `parameter_handler` system-wide object, which is  directly connected with the tutorial CLI. The `process_parameters` TBP registers/defines a set of CLAs to be parsed, parses the CLAs provided by the tutorial user through the CLI, and internally stores their values into a parameter dictionary of _<key,value>_ pairs; a pointer to such dictionary can be obtained by calling `parameter_handler%get_values()` later on.

```fortran
subroutine setup_parameter_handler()
  call parameter_handler%process_parameters(&
       define_user_parameters_procedure=define_tutorial_clas ,&
       progname = tutorial_name ,&
       version = tutorial_version ,&
       description = tutorial_description ,&
       authors = tutorial_authors)
end subroutine setup_parameter_handler

```

@note
This parameter dictionary, with type name `parameterlist_t`, is provided by a stand-alone external software library called [FPL](https://gitlab.com/fempar/FPL)


There are essentially two kind of CLAs registered by `process_parameters`. _On the one hand_, **`FEMPAR`** itself registers a large bunch of CLAs. Each of these CLA corresponds one-to-one to a particular **`FEMPAR`** derived type. The data type a CLA is linked with can be easily inferred from the convention followed for **`FEMPAR`** CLA names, which prefixes the name of the data type (or an abbreviation of it) to the CLA name. Many of the **`FEMPAR`** data types require a set of parameter values in order to customize their behaviour and/or the way they are set up. These data types are designed such that these parameter values may be provided by an instance of the aforementioned parameter dictionary. Thus, by extracting the parameter dictionary stored within `parameter\_handler`, and passing it to the **`FEMPAR`** data type instances, one directly connects the CLI with the instances managed by the FE program. This is indeed the mechanism followed by all tutorial programs. In any case, **`FEMPAR`** users are not forced to use this mechanism in their FE application programs. They can always build and pass an ad-hoc parameter dictionary to the corresponding instance, thus by-passing the parameter values provided to the CLI.

_On the other hand_, the tutorial program itself (or, in general, any FE application program) may optionally register tutorial-specific CLAs. This is achieved by providing a user-declared procedure to the optional `define_user_parameters_procedure` dummy argument of `process_parameters`. In `setup_parameter_handler` subroutine example provided, the particular procedure passed is called `define_tutorial_clas`. In the code below, the reader may observe that registering a CLA involves defining a parameter dictionary _key_ (`FE\_FORMULATION`), a CLA name (`--FE\_FORMULATION`), a default value for the CLA in case it is not passed (`CG`), a help message, and (optionally) a set of admissible choices for the CLA.

```fortran
subroutine define_tutorial_clas()
  call parameter_handler%add( "FE_FORMULATION" , "-- FE_FORMULATION" , "CG" , &
           help = "Select Finite Element formulation for the problem at hand ; &
                   either Continuous (CG) or Discontinuous Galerkin (DG)" , &
           choices = "CG,DG")
  ! ... Register the rest of tutorial-specific CLAs
end subroutine define_tutorial_clas
```

The parameter dictionary key passed when registering a CLA can be used later on in order to get the value of the corresponding CLA or to override it with a fixed value, thus ignoring the value provided to the CLA. This is achieved by means of the `get...()` and `update()` TBPs of `parameter_handler. The subroutine in the code below uses the `getasstring()` TBP of `parameter_handler` in order to obtain the _string_ passed by the tutorial user to the `--FE\_FORMULATION` tutorial specific CLA. 

```fortran
subroutine get_tutorial_cla_values()
  call parameter_handler%getasstring("FE_ FORMULATION", fe_formulation)
  call parameter_handler%get("ALPHA", alpha)
  ... ! Obtain the rest of tutorial-specific CLA values
end subroutine get_tutorial_cla_values
```

The full set of tutorial CLAs, along with rich help messages, can be output on screen by calling the tutorial program with the `--help` CLA, while the full list of parameter dictionary of _<key,value>_ pairs _after parsing_, with the `--PARAMETER_HANDLER_PRINT_VALUES` one. This latter CLA may be useful to confirm that the tutorial program invocation from command-line produces the desired effect on the values actually handled by the program.
