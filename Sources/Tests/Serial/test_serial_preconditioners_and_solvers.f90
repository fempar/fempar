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
program test_serial_preconditioners_and_solvers
  !-----------------------------------------------------------------------
  !
  !  Test program for solvers wrapped in FEMPAR
  !  using a matrix read from matrix market
  !
  !-----------------------------------------------------------------------
  use serial_names
# include "debug.i90"

  implicit none
  ! Files
  type conv_t
    integer(ip), allocatable :: list(:)
  end type conv_t

  type(conv_t), allocatable :: methodstopc(:)
  integer(ip) , allocatable :: list_solvers_to_be_tested(:)
  integer(ip)               :: lunio

  type(serial_scalar_matrix_t), target :: mmmat
  type(serial_scalar_array_t), target :: b
  type(serial_scalar_array_t), target :: x
  type(serial_scalar_array_t), target :: exact_solution 

  class(vector_t) , pointer :: x_base, b_base, exact_solution_base
  class(operator_t), pointer :: A

  type(preconditioner_t)        :: feprec
  type(preconditioner_params_t) :: ppars
  type(solver_control_t)     :: sctrl
  type(serial_environment_t) :: senv

  integer(ip)              :: solver 
  integer(ip)              :: driver 
  integer(ip)              :: i, j, k ,l

  character(len=256)       :: dir_path
  character(len=256)       :: prefix
  character(len=256)       :: name 

  integer(ip)      :: smoother = 2               ! Default value (Symmetric Gauss-Seidel)
  logical          :: one_pass_coarsen = .false. ! Default value
  real(rp)         :: st_parameter=0.25_rp       ! Default value
  real(8)          :: t1, t2, gerror, aerror

  logical          :: from_file
  
  call meminit
  
  call read_pars_cl_test_dirsol_mm ( solver, driver, dir_path, prefix, smoother, one_pass_coarsen, st_parameter, from_file)
  name = trim(prefix) // '.mtx'
  lunio = io_open(trim(dir_path) // '/' // trim(name),status='old')
  call mmmat%read_matrix_market (lunio, symmetric_storage=.false., is_symmetric=.false., sign=indefinite)
  call io_close(lunio)

  ! Alloc vectors
  call b%create_and_allocate (mmmat%graph%nv)
  call x%create_and_allocate (mmmat%graph%nv)
  call exact_solution%create_and_allocate(mmmat%graph%nv)
  call exact_solution%init(1.0_rp)

  A      => mmmat
  x_base => exact_solution 
  b_base => b
  b_base = A*x_base
  x_base => x
  exact_solution_base => exact_solution

  ! Solve using the higher level interface
  sctrl%method=driver
  sctrl%trace=1
  sctrl%itmax=200
  sctrl%dkrymax=200
  ! sctrl%stopc=res_res
  sctrl%orto=icgs

  if(.not.from_file) then

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

     call preconditioner_create  (mmmat, feprec, ppars)
     t1 = wtime()
     call preconditioner_symbolic(mmmat, feprec)
     !call preconditioner_numeric (mmmat, feprec)
     call preconditioner_numeric (feprec)
     t2 = wtime() 

     call preconditioner_log_info(feprec)

     write(*,*) 'Set-up preconditioner time (secs.):', t2-t1
     
     t1 = wtime()
     x%b=0.0_rp
     call abstract_solve(mmmat,feprec,b,x,sctrl,senv)
     t2 = wtime() 
     write(*,'(a,e15.7)') 'Abstract Iterative solution time (secs.):', t2-t1 

     call preconditioner_free_in_stages ( feprec, preconditioner_free_values )
     call preconditioner_free_in_stages ( feprec, preconditioner_free_struct )
     call preconditioner_free_in_stages ( feprec, preconditioner_free_clean )

  else

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! List of Krylov subspace methods available
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! cg               = 1
    ! lgmres           = 2
    ! rgmres           = 3
    ! fgmres           = 4
    ! icg              = 7
    ! lfom             = 8  ! Left preconditioned Full Orthogonalization Method
    ! minres           = 9  ! Preconditioned MINimal RESidual method
    ! Not actually Krylov methods
    ! richard          = 5  ! Richardson (fixed-point iteration)
    ! direct           = 6  ! Apply preconditioner directly

    allocate( methodstopc(9) )

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! List of convergence criteria available for iterative solvers 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! res_nrmgiven_rhs_nrmgiven  = 1  ! ||  r(i) ||g <= rtol*||  b    ||g + atol 
    ! res_nrmgiven_res_nrmgiven  = 2  ! ||  r(i) ||g <= rtol*||  r(0) ||g + atol   
    ! delta_rhs                  = 3  ! || dx(i) ||  <= rtol*||  b  || + atol
    ! delta_delta                = 4  ! || dx(i) ||  <= rtol*||dx(1)|| + atol
    ! res_res                    = 5  ! ||  r(i) ||  <= rtol*|| r(0)|| + atol
    ! res_rhs                    = 6  ! ||  r(i) ||  <= rtol*||  b  || + atol
    ! delta_rhs_and_res_res      = 7  ! delta_rhs    AND res_res
    ! delta_rhs_and_res_rhs      = 8  ! delta_rhs    AND res_rhs
    ! delta_delta_and_res_res    = 9  ! delta_delta  AND res_res
    ! delta_delta_and_res_rhs    = 10 ! delta_delta  AND res_rhs 

    allocate( methodstopc( cg )%list(10) );     methodstopc( cg )%list =      (/1,2,3,4,5,6,7,8,9,10/)
    allocate( methodstopc( lgmres )%list(4) );  methodstopc( lgmres )%list =  (/1,2,5,6/)
    allocate( methodstopc( rgmres )%list(2) );  methodstopc( rgmres )%list =  (/1,2/)
    allocate( methodstopc( fgmres )%list(2) );  methodstopc( fgmres )%list =  (/1,2/)
    allocate( methodstopc( richard )%list(2) ); methodstopc( richard )%list = (/5,6/)
    allocate( methodstopc( direct )%list(10) ); methodstopc( direct )%list =  (/1,2,3,4,5,6,7,8,9,10/)
    allocate( methodstopc( icg )%list(2) );    methodstopc( icg )%list =     (/5,6/)
    allocate( methodstopc( lfom )%list(2) );    methodstopc( lfom )%list =    (/5,6/)
    allocate( methodstopc( minres )%list(1) );  methodstopc( minres )%list =  (/5/)

    allocate(list_solvers_to_be_tested(7))
    list_solvers_to_be_tested = (/cg, lgmres, rgmres, fgmres, icg, lfom, minres/)

    ! Loop in krylov methods
    do i=1,size(list_solvers_to_be_tested)
        sctrl%method=list_solvers_to_be_tested(i) 

        ! Loop in allowed stop conditions for each krylov methods
        do j=1,size(methodstopc(list_solvers_to_be_tested(i))%list,1)
            sctrl%stopc=methodstopc(list_solvers_to_be_tested(i))%list(j)

            ! Loop in allowed solvers. When a solver is not allowed: cycle
            do k=0,5
                select case(k)
                    case(0)
                        ppars%type = no_prec
                    case(1)
#ifdef ENABLE_MKL
                        ppars%type = pardiso_mkl_prec
#else
                        cycle
#endif
                    case(2)
#ifdef ENABLE_WSMP
                        ppars%type = wsmp_prec
#else
                        cycle
#endif
                    case(3)
#ifdef ENABLE_HSL_MI20
                        ppars%type             = hsl_mi20_prec
                        ppars%pre_smoothing    = 1
                        ppars%post_smoothing   = 1 
                        ppars%smoother         = smoother  
                        ppars%one_pass_coarsen = one_pass_coarsen
                        ppars%st_parameter     = st_parameter
                        ppars%verbosity        = 0 
                        ppars%c_fail           = 2
#else
                        cycle
#endif
                    case(4)
#ifdef ENABLE_HSL_MA87
!						 Needs a symmetric positive definite matrix
!                        ppars%type = hsl_ma87_prec 
                        cycle
#else
                        cycle
#endif
                    case(5)
#ifdef ENABLE_UMFPACK
                        ppars%type = umfpack_prec
#else
                        cycle
#endif
                    case default
                        write(*,*) 'solver must be 0 (none), 1 (pardiso), 2 (wsmp), 3(hsl_mi20), 4(hsl_ma87), 5(umfpack)'
                        stop
                end select

                ! Loop in ortogonalization methods
                do l=1,2
                    sctrl%orto=l

                    call preconditioner_create  (mmmat, feprec, ppars)
                    t1 = wtime()
                    call preconditioner_symbolic(mmmat, feprec)
                    !call preconditioner_numeric (mmmat, feprec)
                    call preconditioner_numeric (feprec)
                    t2 = wtime() 
        
                    call preconditioner_log_info(feprec)
        
                    write(*,*) 'Set-up preconditioner time (secs.):', t2-t1
 
                    t1 = wtime()
                    call x%init(0.0_rp)
                    call abstract_solve(mmmat,feprec,b,x,sctrl,senv)
                    t2 = wtime() 
                    write(*,'(a,e15.7)') 'Abstract Iterative solution time (secs.):', t2-t1 
        
                    call preconditioner_free_in_stages ( feprec, preconditioner_free_values )
                    call preconditioner_free_in_stages ( feprec, preconditioner_free_struct )
                    call preconditioner_free_in_stages ( feprec, preconditioner_free_clean )
        
                    x_base = x_base-exact_solution_base
                    if (x%nrm2()/exact_solution%nrm2() > 1.e-06 ) then
                        check(.false.)
                    endif

                enddo

            enddo

        enddo
        deallocate(methodstopc(i)%list)
    enddo

    deallocate( methodstopc )
    deallocate( list_solvers_to_be_tested )
  endif


  call mmmat%free()
  call b%free()
  call x%free()
  call exact_solution%free()
  call memstatus
contains

  ! *****************************************************************************!
  ! Read params from command-line options.                                       ! 
  ! Command-line options processing for f90 is discussed, e.g.,                  !
  ! on the following URL: http://people.sc.fsu.edu/~jburkardt/f_src/args/args.f90!
  ! Still to confirm whether this support is standard in f90 or depends          !
  ! on the compiler (i.e., INTEL, GNU, etc.)                                     !
  ! *****************************************************************************!
  subroutine read_pars_cl_test_dirsol_mm (solver, driver, dir_path, prefix, & 
                                          smoother, one_pass_coarsen, st_parameter,from_file)
    implicit none
    character*(*), intent(out)   :: dir_path, prefix
    integer(ip)  , intent(out)   :: solver, driver
    character(len=256)           :: program_name
    character(len=256)           :: argument 
    integer                      :: numargs,iargc
    integer(ip), intent(out)     :: smoother
    logical    , intent(out)     :: one_pass_coarsen
    real(rp)   , intent(out)     :: st_parameter
    logical, intent(out)         :: from_file
    logical                      :: file_exists


    from_file=.false.
    numargs = iargc()
    call getarg(0, program_name)
    if((numargs == 1)) then
        ! Check if file exists
        call getarg(1, argument)
        call read_pars_cl_test_dirsol_mm_from_file (argument, dir_path, prefix)
        from_file=.true.
    elseif ((numargs < 4) ) then 
        call print_usage(program_name)
        check(.false.)
    else

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
              check(.false.)
           end if
        end if

    end if

  end subroutine read_pars_cl_test_dirsol_mm

  subroutine print_usage(program_name)
    implicit none
    character*(*), intent(in) :: program_name
    write (6,'(a)') 'Usage: ', trim(program_name)//' solver driver dir_path prefix [smoother one_pass_coarsen st_parameter]'
    write (6,'(a)') '          where solver=1 is pardiso, solver=2 is wsmp, solver=3 is hsl_mi20, solver=4 is hsl_ma87, solver=5 is umfpack'
    write (6,'(a)') '          and driver is cg=1, lgmres=2, rgmres=3, fgmres=4, richard=5, direct=6, icg=7, lfom=8, minres=9'
    write (6,'(a)') 'Or: ', trim(program_name)//' filename.txt'
    write (6,'(a)') '          where filename.txt is a ASCII text file containing the command line parameters, one per line, '
    write (6,'(a)') '          in the same order as in the command line.'
  end subroutine print_usage


  subroutine read_pars_cl_test_dirsol_mm_from_file (filename, dir_path, prefix)
    implicit none
    character(len=256), intent(in)           :: filename
    character*(*), intent(out)   :: dir_path, prefix
    integer                      :: iu, ios
    logical                      :: file_exists


        inquire(file=trim(filename), exist=file_exists)
        if(.not. file_exists) then
            write (6,'(a)') ' ERROR: Input file '//trim(filename)//' not found!'
            check(.false.)
        else

            iu = io_open(trim(filename),status='old')
!            open(unit=10, file=trim(filename), form='FORMATTED', access='SEQUENTIAL', action = 'READ', iostat = ios)
!            if(ios /= 0) call io_check(ios,filename)

            read(unit=iu, fmt='(a)', iostat = ios) dir_path
            call io_check(ios,filename)
            read(unit=iu, fmt=*, iostat = ios) prefix
            call io_check(ios,filename)
            
            if (solver==3) then
                read(unit=iu, fmt=*, iostat = ios) smoother
                call io_check(ios,filename)
                read(unit=iu, fmt=*, iostat = ios) one_pass_coarsen
                call io_check(ios,filename)
                read(unit=iu, fmt=*, iostat = ios) st_parameter
                call io_check(ios,filename)
            end if
            close(iu)
        endif

  end subroutine read_pars_cl_test_dirsol_mm_from_file


  subroutine io_check(ios, filename)
    implicit none
    integer            :: ios
    character(len=256), intent(in) :: filename

        if(ios /= 0) then 
            write (6,'(a,i3,a)') ' ERROR: Reading text file '//trim(filename)//' (',ios,')'
            check(.false.)
        endif

  end subroutine io_check


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

end program test_serial_preconditioners_and_solvers
