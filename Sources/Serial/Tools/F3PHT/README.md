# F3PHT
**F**ortran **P**olymorphic **P**rocedure **P**ointer **H**ash **T**able

## License

[LGPL v3](http://www.gnu.org/licenses/lgpl-3.0.txt)

## What is FPHT?

**F3PHT** is pure fortran 2003 library that implements a polymorphic procedure pointer hash table.

**F3PHT** is a container (dictionary) of ```<Key, Proc>``` pairs, where the *Key* is a character string and *Proc* is a polymorphic procedure pointer

## How to get F3PHT

```git clone --recursive http://servercomfus/LSSC/F3PHT.git ```

## Compilation

**F3PHT** compile with GNU Fortran compiler 5.3 (and newer versions), Intel Fortran compiler 16.0.0 and IBM XLF 14.0.1 (and newer versions).

**F3PHT** uses [CMake](https://cmake.org/) as a portable compilation system. 

The easiest way to compile **F3PHT** under Linux is:

```
$ cd F3PHT
$ mkdir build
$ cd build
$ cmake ../
$ make
```

*To compile FPHT under Windows use de equivalent commands*

## Getting started with F3PHT


### Using FPHT in your program

```fortran
USE F3PHT

type(ProcedurePointer_t) :: My_Dictionary

call My_Dictionary%Init()

... [Program boddy]

call My_Dictionary%Free()
```

### Setting tuples

```fortran
...
subroutine one_int_arg_subroutine(a)
    integer, intent(in) :: a
    write(*,*) 'Calling a subroutine with 1 integer argument:', a
end subroutine

subroutine one_real_arg_subroutine(a)
    real, intent(in) :: a
    write(*,*) 'Calling a subroutine with 1 real argument:', a
end subroutine
...
procedure(), pointer :: my_procedure_pointer => NULL()
...
my_procedure_pointer => one_int_arg_subroutine
call My_Dictionary%Set(Key='Procedure1', Proc=my_procedure_pointer)
my_procedure_pointer => one_reañ_arg_subroutine
call My_Dictionary%Set(Key='Procedure2', Proc=my_procedure_pointer)
...
```

### Getting values

```fortran
procedure(), pointer :: my_procedure_pointer => NULL()

call My_Dictionary%Get(Key='Procedure1', Value=my_procedure_pointer)
```

### Deleting tuples

```fortran
call My_Dictionary%Del(Key='Procedure1')
```

### Checking if a key is present

```fortran
logical :: procedure1_is_present

procedure1_is_present = My_Dictionary%isPresent(Key='procedure1')
```

### Print the content of the dictionary

```fortran
call My_Dictionary%Print()
```

