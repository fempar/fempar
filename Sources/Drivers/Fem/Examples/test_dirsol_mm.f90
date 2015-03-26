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
# include "debug.i90"

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
  real(8)          :: t1, t2, gerror, aerror

  logical          :: from_file


  call meminit
  
  call read_pars_cl_test_dirsol_mm ( solver, driver, dir_path, prefix, smoother, one_pass_coarsen, st_parameter, from_file)

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
  ! sctrl%stopc=res_res
  sctrl%orto=icgs

  if(.not.from_file) then

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

  else
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Iterate on this methods
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! List of Krylov subspace methods available
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
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i=1,9
        sctrl%method=i !driver
        sctrl%stopc=res_res
        if(i==3 .or. i==4) sctrl%stopc=res_nrmgiven_rhs_nrmgiven

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
        gerror = sctrl%err1
        ! call solver_control_log_conv_his(sctrl)
        call solver_control_free_conv_his(sctrl)

        t1 = wtime()
        feunk%b=1.0_rp
        fevec%b=1.0_rp
        call abstract_solve(mmmat,feprec,fevec,feunk,sctrl)
        t2 = wtime() 
        write(*,'(a,e15.7)') 'Abstract Iterative solution time (secs.):', t2-t1 
        aerror = sctrl%err1
        ! call solver_control_log_conv_his(sctrl)
        call solver_control_free_conv_his(sctrl)

        call fem_precond_free ( precond_free_values, feprec)
        call fem_precond_free ( precond_free_struct, feprec)
        call fem_precond_free ( precond_free_clean, feprec)

        if((gerror-aerror)>1000*epsilon(gerror)) then
            ! check generic-abstract
            check(.false.)
        endif
    enddo

  endif


  call fem_graph_free ( mmgraph )
  call fem_matrix_free ( mmmat ) 
  call fem_vector_free (fevec)
  call fem_vector_free (feunk)

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
        call read_pars_cl_test_dirsol_mm_from_file (argument, solver, driver, dir_path, prefix, & 
                                          smoother, one_pass_coarsen, st_parameter)
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
    write (6,'(a)') '          where solver=1 is pardiso, solver=2 is wsmp, solver=3 is hsl_mi20, solver=4 is hsl_ma87, solver=4 is umfpack'
    write (6,'(a)') '          and driver is cg=1, lgmres=2, rgmres=3, fgmres=4, richard=5, direct=6, icg=7, lfom=8, minres=9'
    write (6,'(a)') 'Or: ', trim(program_name)//' filename.txt'
    write (6,'(a)') '          where filename.txt is a ASCII text file containing the command line parameters, one per line, '
    write (6,'(a)') '          in the same order as in the command line.'
  end subroutine print_usage


  subroutine read_pars_cl_test_dirsol_mm_from_file (filename, solver, driver, dir_path, prefix, & 
                                          smoother, one_pass_coarsen, st_parameter)
    implicit none
    character(len=256), intent(in)           :: filename
    character*(*), intent(out)   :: dir_path, prefix
    integer(ip)  , intent(out)   :: solver, driver
    character(len=256)           :: program_name
    character(len=256)           :: argument 
    integer                      :: numargs,iargc, iu, ios
    integer(ip), intent(out)     :: smoother
    logical    , intent(out)     :: one_pass_coarsen
    real(rp)   , intent(out)     :: st_parameter
    logical                      :: file_exists


        inquire(file=trim(filename), exist=file_exists)
        if(.not. file_exists) then
            write (6,'(a)') ' ERROR: Input file '//trim(filename)//' not found!'
            check(.false.)
        else

            iu = io_open(trim(filename),status='old')
!            open(unit=10, file=trim(filename), form='FORMATTED', access='SEQUENTIAL', action = 'READ', iostat = ios)
!            if(ios /= 0) call io_check(ios,filename)


            read(unit=iu, fmt=*, iostat = ios) solver
            call io_check(ios,filename)
            read(unit=iu, fmt=*, iostat = ios) driver
            call io_check(ios,filename)
            read(unit=iu, fmt=*, iostat = ios) dir_path
            call io_check(ios,filename)
            ! dirpath is a relative path from filename folder
            dir_path=trim(filename(1:index( filename, '/', back=.true.)))//trim(dir_path)
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

end program test_dirsol_mm
