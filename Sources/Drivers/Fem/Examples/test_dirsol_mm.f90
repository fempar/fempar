! Copyright (C) 2014 Santiago Badia, Alberto F. Martín and Javier Principe
!
! This file is part of FEMPAR (Finite Element Multiphysics PARallel library)
!
! FEMPAR is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! FEMPAR is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with FEMPAR. If not, see <http://www.gnu.org/licenses/>.
!
! Additional permission under GNU GPL version 3 section 7
!
! If you modify this Program, or any covered work, by linking or combining it 
! with the Intel Math Kernel Library and/or the Watson Sparse Matrix Package 
! and/or the HSL Mathematical Software Library (or a modified version of them), 
! containing parts covered by the terms of their respective licenses, the
! licensors of this Program grant you additional permission to convey the 
! resulting work. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_dirsol_mm
  !-----------------------------------------------------------------------
  !
  !  Test program for direct solvers wrapped in FEMPAR
  !  using a matrix read from matrix market
  !
  !-----------------------------------------------------------------------
  use fem
  implicit none
  ! Files
  integer(ip)              :: lunio

  type(fem_matrix)         :: mmmat
  type(fem_graph)          :: mmgraph
  type(fem_vector)         :: fevec
  type(fem_vector)         :: feunk

  type(fem_precond)        :: feprec
  type(fem_precond_params) :: ppars
  type(solver_control)     :: sctrl

  integer(ip)              :: solver 
  integer(ip)              :: driver 
  integer(ip)              :: i, j 

  character(len=256)       :: dir_path
  character(len=256)       :: prefix
  character(len=256)       :: name 

  integer(ip)      :: smoother
  logical          :: one_pass_coarsen
  real(rp)         :: st_parameter
  real(8)          :: t1, t2

  call meminit
  
  call read_pars_cl_test_dirsol_mm ( solver, driver, dir_path, prefix, smoother, one_pass_coarsen, st_parameter)

  ! Read a symmetric matrix
  call fem_matrix_compose_name_matrix_market ( prefix, name ) 
  lunio = io_open(trim(dir_path) // '/' // trim(name),status='old')
  ! call fem_matrix_read_matrix_market (lunio, mmmat, mmgraph, symm_true, positive_definite)
  call fem_matrix_read_matrix_market (lunio, mmmat, mmgraph, symm_false)

  call io_close(lunio)

!!$  prefix = trim(prefix) // '.out'
!!$  write(*,*) prefix
!!$  call fem_matrix_compose_name_matrix_market ( prefix , name ) 
!!$  lunio = io_open(trim(dir_path) // '/' // trim(name))
!!$  call fem_matrix_print_matrix_market (lunio, mmmat)
!!$  call io_close(lunio)

  ! Alloc vectors
  call fem_vector_alloc (scal,1,mmmat%gr%nv,fevec)
  call fem_vector_alloc (scal,1,mmmat%gr%nv,feunk)

  ! Solve using the higher level interface
  select case(solver)
  case(0)
     ppars%type = no_prec
  case(1)
     ppars%type = pardiso_mkl_prec
  case(2)
     ppars%type = wsmp_prec
  case(3)
     ppars%type             = hsl_mi20_prec
     ppars%pre_smoothing    = 1
     ppars%post_smoothing   = 1 
     ppars%smoother         = smoother  
     ppars%one_pass_coarsen = one_pass_coarsen
     ppars%st_parameter     = st_parameter
     ppars%verbosity        = 0 
     ppars%c_fail           = 2
  case(4)
     ppars%type = hsl_ma87_prec 
  case(5)
     ppars%type = umfpack_prec
  case default
     write(*,*) 'solver must be 0 (none), 1 (pardiso), 2 (wsmp), 3(hsl_mi20), 4(hsl_ma87), 5(umfpack)'
     stop
  end select

  sctrl%method=driver
  sctrl%trace=1
  sctrl%itmax=200
  sctrl%dkrymax=200
  ! sctrl%stopc=res_rhs
  sctrl%orto=icgs

  do i=1,1

     call fem_precond_create  (mmmat, feprec, ppars)
     t1 = wtime()
     call fem_precond_symbolic(mmmat, feprec)
     call fem_precond_numeric (mmmat, feprec)
     t2 = wtime() 

     call fem_precond_log_info(feprec)

     write(*,*) 'Set-up preconditioner time (secs.):', t2-t1
     
     t1 = wtime()
     feunk%b=1.0_rp
     fevec%b=1.0_rp
     call solve(mmmat,feprec,fevec,feunk,sctrl)
     t2 = wtime() 
     write(*,'(a,e15.7)') 'Generic Iterative solution time (secs.):', t2-t1 

     ! call solver_control_log_conv_his(sctrl)
     call solver_control_free_conv_his(sctrl)

     t1 = wtime()
     feunk%b=1.0_rp
     fevec%b=1.0_rp
     call abstract_solve(mmmat,feprec,fevec,feunk,sctrl)
     t2 = wtime() 
     write(*,'(a,e15.7)') 'Abstract Iterative solution time (secs.):', t2-t1 

     ! call solver_control_log_conv_his(sctrl)
     call solver_control_free_conv_his(sctrl)

     call fem_precond_free ( precond_free_values, feprec)
     call fem_precond_free ( precond_free_struct, feprec)
     call fem_precond_free ( precond_free_clean, feprec)

  end do


  call fem_graph_free ( mmgraph )
  call fem_matrix_free ( mmmat ) 


contains

  ! *****************************************************************************!
  ! Read params from command-line options.                                       ! 
  ! Command-line options processing for f90 is discussed, e.g.,                  !
  ! on the following URL: http://people.sc.fsu.edu/~jburkardt/f_src/args/args.f90!
  ! Still to confirm whether this support is standard in f90 or depends          !
  ! on the compiler (i.e., INTEL, GNU, etc.)                                     !
  ! *****************************************************************************!
  subroutine read_pars_cl_test_dirsol_mm (solver, driver, dir_path, prefix, & 
                                          smoother, one_pass_coarsen, st_parameter)
    implicit none
    character*(*), intent(out)   :: dir_path, prefix
    integer(ip)  , intent(out)   :: solver, driver
    character(len=256)           :: program_name
    character(len=256)           :: argument 
    integer                      :: numargs,iargc
    integer(ip), intent(out)     :: smoother
    logical    , intent(out)     :: one_pass_coarsen
    real(rp)   , intent(out)     :: st_parameter

    numargs = iargc()
    call getarg(0, program_name)
    if ((numargs < 4) ) then 
       call print_usage(program_name)
       stop
    end if

    call getarg(1, argument)
    read (argument,*) solver     
    
    call getarg(2, argument)
    read (argument,*) driver     
    
    call getarg(3, argument)
    dir_path = trim(argument)
    
    call getarg(4, argument)
    prefix = trim(argument)
    
    if (solver==3) then
       if (.not. (numargs == 7) ) then
          call print_usage(program_name)
       end if
       
       call getarg(5,argument)
       read (argument,*) smoother
       
       call getarg(6, argument)
       if ( trim(argument) .eq. 'T' ) then
          one_pass_coarsen = .true.
       else if ( trim(argument) .eq. 'F' ) then
          one_pass_coarsen = .false.
       end if
       
       call getarg(7, argument)
       read (argument,*) st_parameter
    else
       if (.not. (numargs == 4) ) then
          call print_usage(program_name)
          stop
       end if
    end if
  end subroutine read_pars_cl_test_dirsol_mm

  subroutine print_usage(program_name)
    implicit none
    character*(*), intent(in) :: program_name
    write (6,'(a)') 'Usage: ', trim(program_name), ' solver driver dir_path prefix [smoother one_pass_coarsen st_parameter]'
    write (6,'(a)') '          where solver=1 is pardiso, solver=2 is wsmp, solver=3 is hsl_mi20, solver=4 is hsl_ma87, solver=4 is umfpack'
    write (6,'(a)') '          and driver is cg=1, lgmres=2, rgmres=3, fgmres=4, richard=5, direct=6, icg=7, lfom=8, minres=9'
  end subroutine print_usage


  function wtime ( )

  !*****************************************************************************80
  !
  !! WTIME returns a reading of the wall clock time.
  !
  !  Discussion:
  !
  !    To get the elapsed wall clock time, call WTIME before and after a given
  !    operation, and subtract the first reading from the second.
  !
  !    This function is meant to suggest the similar routines:
  !
  !      "omp_get_wtime ( )" in OpenMP,
  !      "MPI_Wtime ( )" in MPI,
  !      and "tic" and "toc" in MATLAB.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    27 April 2009
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Output, real ( kind = 8 ) WTIME, the wall clock reading, in seconds.
  !
    implicit none
  
    integer ( kind = 4 ) clock_max
    integer ( kind = 4 ) clock_rate
    integer ( kind = 4 ) clock_reading
    real ( kind = 8 ) wtime

    call system_clock ( clock_reading, clock_rate, clock_max )

    wtime = real ( clock_reading, kind = 8 ) &
          / real ( clock_rate, kind = 8 )

    return
  end function wtime

end program test_dirsol_mm
