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
module umfpack_class
  ! This module is a wrapper in which we define functions 
  ! to interact with UMFPACK using our data structures.
  ! Error control on calling parameters must (ideally)
  ! be performed here.

  ! Serial modules
  use iso_c_binding
  use umfpack_interface
  use types
  use memor
  use fem_matrix_class
  use fem_vector_class 
  use fem_graph_class
  use renum_class

# include "debug.i90"
  
  implicit none
  private

  ! Possible states of umfpack_context
  integer(ip), parameter :: not_created   = 0
  integer(ip), parameter :: created       = 1 
  integer(ip), parameter :: symb_computed = 2 ! Symbolic data already computed
  integer(ip), parameter :: num_computed  = 3 ! Numerical data already computed 


  type umfpack_context
     ! Our components
     integer(ip)    :: state = not_created
     type(c_ptr)    :: Symbolic
     type(c_ptr)    :: Numeric
#ifdef ENABLE_UMFPACK
     real(c_double) :: Control(0:UMFPACK_CONTROL-1)            
     real(c_double) :: Info(0:UMFPACK_INFO-1)            
#endif
  end type umfpack_context

  ! Types
  public :: umfpack_context

  ! Possible actions that can be perfomed by solve_umfpack
  integer(ip), parameter :: umfpack_init              = 1  ! Construct solve_umfpack_state
  integer(ip), parameter :: umfpack_finalize          = 2  ! Destruct  solve_umfpack_state     
  integer(ip), parameter :: umfpack_compute_symb      = 3  ! Compute symb. fact.  
  integer(ip), parameter :: umfpack_compute_num       = 4  ! Compute numerical fact. 
  integer(ip), parameter :: umfpack_solve             = 6  ! Fwd./Bck. substitution 
  integer(ip), parameter :: umfpack_free_values       = 7
  integer(ip), parameter :: umfpack_free_struct       = 8
  integer(ip), parameter :: umfpack_free_clean        = 9

  interface umfpack
     module procedure umfpack_vector, umfpack_r2, umfpack_r1
  end interface umfpack

  ! Constants (actions)
  public :: umfpack_init,  &
       &    umfpack_finalize, &
       &    umfpack_compute_symb, & 
       &    umfpack_compute_num, &
       &    umfpack_solve, &
       &    umfpack_free_values,&
       &    umfpack_free_struct, &
       &    umfpack_free_clean

  ! Functions
  public :: umfpack

contains
  !=============================================================================
  subroutine umfpack_vector ( action, context, A, b, x )
    implicit none
    ! Mandatory Parameters
    type(umfpack_context), intent(inout) :: context   ! Information required between calls
    integer(ip)           , intent(in)    :: action    ! Action to be performed 
                                                   ! (see public constants above)
    type(fem_matrix)  , intent(in)    :: A         ! Linear system coefficient matrix
    type(fem_vector)  , intent(in)    :: b         ! RHS (Right-hand-side)
    type(fem_vector)  , intent(inout) :: x         ! LHS (Left-hand-side)

    select case(action)

    case ( umfpack_init )
       ! Check pre-conditions
       assert ( context%state == not_created )
       call umfpack_ini ( context, A )
       ! State transition 
       context%state = created
    case ( umfpack_finalize )

       ! Check pre-conditions
       assert ( context%state /= not_created )
       select case (context%state)
       case (created)
          call  umfpack_free ( umfpack_free_clean , context )
       case (symb_computed)
          call  umfpack_free ( umfpack_free_struct, context )
          call  umfpack_free ( umfpack_free_clean , context )
       case (num_computed)
          call  umfpack_free ( umfpack_free_values, context )
          call  umfpack_free ( umfpack_free_struct, context )
          call  umfpack_free ( umfpack_free_clean , context )
       end select
       ! State transition 
       context%state = not_created

    case ( umfpack_compute_symb )

       ! Check pre-conditions
       assert (context%state==created)
       ! Create sparsity pattern of LU sparse direct factorization
       call umfpack_analysis ( context, A )
       ! State transition 
       context%state = symb_computed

    case ( umfpack_compute_num )

       ! Check pre-conditions
       assert (context%state==symb_computed.or.context%state==num_computed)
       ! Compute LU sparse direct factorization
       call umfpack_factorization ( context, A )
       ! State transition 
       context%state = num_computed

    case ( umfpack_solve )
       ! Check pre-conditions
       assert ( context%state ==  num_computed )
       call umfpack_solution ( context, A, b, x )

    case ( umfpack_free_values )

       call  umfpack_free ( umfpack_free_values, context )
       context%state=symb_computed

    case ( umfpack_free_struct )

       call  umfpack_free ( umfpack_free_struct, context )
       context%state=created

    case ( umfpack_free_clean )

       call  umfpack_free ( umfpack_free_clean, context )
       context%state=not_created

    case default

       ! Write an error message and stop ?      

    end select

  end subroutine umfpack_vector
  !=============================================================================
  subroutine umfpack_r2 ( action, context, A, nrhs, b, ldb, x, ldx)
    implicit none
    ! Mandatory Parameters
    type(umfpack_context), intent(inout) :: context   ! Information required between calls
    integer(ip)       , intent(in)    :: action    ! Action to be performed 
                                                   ! (see public constants above)
    type(fem_matrix)  , intent(in)    :: A         ! Linear system coefficient matrix
!    type(fem_vector)  , intent(in)    :: b         ! RHS (Right-hand-side)
!    type(fem_vector)  , intent(inout) :: x         ! LHS (Left-hand-side)
    integer(ip)       , intent(in)    :: nrhs, ldb, ldx
    real(rp)          , intent(in)    :: b (ldb, nrhs)
    real(rp)          , intent(inout) :: x (ldx, nrhs)


    select case(action)

    case ( umfpack_solve )
       
       ! Check pre-conditions
       assert ( context%state ==  num_computed )
       call umfpack_solution_several_rhs ( context, A, nrhs, b, ldb, x, ldx  )
       
    case default
       
       ! Write an error message and stop ?      

    end select

  end subroutine umfpack_r2

  !=============================================================================
  subroutine umfpack_r1 ( action, context, A, b, x )
    implicit none
    ! Mandatory Parameters
    type(umfpack_context), intent(inout) :: context   ! Information required between calls
    integer(ip)       , intent(in)    :: action    ! Action to be performed 
                                                   ! (see public constants above)
    type(fem_matrix)  , intent(in)    :: A         ! Linear system coefficient matrix
    real(rp)          , intent(in)    :: b (A%gr%nv)
    real(rp)          , intent(inout) :: x (A%gr%nv)


    select case(action)

    case ( umfpack_solve )

       ! Check pre-conditions
       assert ( context%state ==  num_computed )
       call umfpack_solution_real ( context, A, b, x  )

    case default

       ! Write an error message and stop ?      

    end select
    
  end subroutine umfpack_r1

  !=============================================================================
#ifndef ENABLE_UMFPACK
  subroutine enable_umfpack_error_message
      implicit none
      write (0,*) 'Error: Fem was not compiled with -DENABLE_UMFPACK.'
      write (0,*) "Error: You must activate this cpp macro in order to use UMFPACK"
      check(1==0)
  end subroutine
#endif

  !=============================================================================
  subroutine umfpack_ini ( context, matrix )
    implicit none
    ! Parameters
    type(umfpack_context) , intent(inout) :: context
    type(fem_matrix)      , intent(in)    :: matrix

#ifdef ENABLE_UMFPACK
    !
    !  Set the default control parameters.
    !
    call umfpack_di_defaults (context%Control)
#else
    call enable_umfpack_error_message
#endif
  end subroutine umfpack_ini

  !=============================================================================
  subroutine umfpack_free ( mode, context)
    implicit none
    ! Parameters
    integer(ip)          , intent(in)    :: mode
    type(umfpack_context), intent(inout) :: context

#ifdef ENABLE_UMFPACK
    ! Free umfpack_context structures
    if ( mode == umfpack_free_clean ) then
      
    else if ( mode == umfpack_free_struct  ) then
      !
      !  Free the memory associated with the symbolic factorization.
      !
      call umfpack_di_free_symbolic ( context%Symbolic )
    else if ( mode == umfpack_free_values ) then
      !
      !  Free the memory associated with the numeric factorization.
      !
      call umfpack_di_free_numeric ( context%Numeric )
    end if
#else
    call enable_umfpack_error_message
#endif

  end subroutine umfpack_free
 
  !=============================================================================
  subroutine increment_array ( v )
    implicit none
    integer(ip) :: v(:)
    v = v + 1
  end subroutine increment_array

  !=============================================================================
  subroutine decrement_array ( v )
    implicit none
    integer(ip) :: v(:)
    v = v - 1
  end subroutine decrement_array

  !=============================================================================
  subroutine umfpack_analysis ( context, matrix )
    use fem_mesh_partition_base
    use graph_renum
    implicit none

    ! Parameters 
    type(umfpack_context)  , intent(inout) :: context
    type(fem_matrix)       , intent(in)    :: matrix

    ! Locals
    integer(ip) :: status
    
    assert ( matrix%gr%type == csr )

#ifdef ENABLE_UMFPACK
    ! Fortran to C numbering 
    call decrement_array( matrix%gr%ia )
    call decrement_array( matrix%gr%ja )

    !
    !  From the matrix data, create the symbolic factorization information.
    !
    status = umfpack_di_symbolic ( matrix%gr%nv, matrix%gr%nv2, matrix%gr%ia, matrix%gr%ja, C_NULL_PTR, context%Symbolic, context%Control, context%Info)
    if ( status < 0 ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'UMFPACK - Fatal error!'
      write ( *, '(a,i10)' ) '  umfpack_di_symbolic returns STATUS = ', status
      check ( status == UMFPACK_OK )  
    end if
    
     ! C to Fortran numbering 
    call increment_array( matrix%gr%ia )
    call increment_array( matrix%gr%ja )
#else
    call enable_umfpack_error_message
#endif


  end subroutine umfpack_analysis

  !=============================================================================
  subroutine umfpack_factorization ( context, matrix )
    implicit none
    ! Parameters 
    type(umfpack_context) , intent(inout) :: context
    type(fem_matrix)       , target, intent(in) :: matrix

    ! Locals
    real(rp), pointer :: a_(:)
    integer(ip)       :: status

#ifdef ENABLE_UMFPACK
    assert ( matrix%type == csr_mat )
    assert ( matrix%symm == symm_false )

    ! Fortran to C numbering 
    call decrement_array( matrix%gr%ia )
    call decrement_array( matrix%gr%ja )

    a_ => matrix%a(1,1,:)

    !
    !  From the symbolic factorization information, carry out the numeric factorization.
    !
    status = umfpack_di_numeric ( matrix%gr%ia, matrix%gr%ja, a_, context%Symbolic, context%Numeric, context%Control, context%Info )
  
    if ( status < 0 ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'UMFPACK - Fatal error!'
      write ( *, '(a,i10)' ) '  umfpack_di_numeric returns status = ', status
      check ( status == UMFPACK_OK ) 
    end if

     ! C to Fortran numbering 
    call increment_array( matrix%gr%ia )
    call increment_array( matrix%gr%ja )
#else
    call enable_umfpack_error_message
#endif

  end subroutine umfpack_factorization

  !=============================================================================
  subroutine umfpack_solution ( context, matrix, x, y )
    ! Computes y <- A^-1 * x, using previously computed LU factorization
    implicit none

    ! Parameters 
    type(umfpack_context), intent(inout)         :: context
    type(fem_matrix)     , intent(in), target    :: matrix
    type(fem_vector)     , intent(in), target    :: x
    type(fem_vector)     , intent(inout), target :: y

    ! Locals
    real(rp), pointer :: x_(:)
    real(rp), pointer :: y_(:)
    real(rp), pointer :: a_(:)
    integer(ip)       :: status

#ifdef ENABLE_UMFPACK
    ! Fortran to C numbering 
    call decrement_array( matrix%gr%ia )
    call decrement_array( matrix%gr%ja )

    x_ => x%b(1,:)
    y_ => y%b(1,:)
    a_ => matrix%a(1,1,:)

    !
    ! Solve the linear system.
    !
    status = umfpack_di_solve ( UMFPACK_At, matrix%gr%ia, matrix%gr%ja, a_, y_, x_, context%Numeric, context%Control, context%Info )

    if ( status < 0 ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'UMFPACK - Fatal error!'
      write ( *, '(a,i10)' ) '  umfpack_di_solve returns status = ', status
      check ( status == UMFPACK_OK ) 
    end if

    ! C to Fortran numbering 
    call increment_array( matrix%gr%ia )
    call increment_array( matrix%gr%ja )
#else
    call enable_umfpack_error_message
#endif

  end subroutine umfpack_solution

  !=============================================================================
  subroutine umfpack_solution_real ( context, matrix, rhs, sol )
    implicit none
    ! Parameters 
    type(umfpack_context), intent(inout)   :: context
    type(fem_matrix)  , intent(in), target :: matrix
    real(rp)          , intent(in)         :: rhs (matrix%gr%nv)
    real(rp)          , intent(inout)      :: sol (matrix%gr%nv)
    
    ! Locals 
    real(rp), pointer :: a_(:)
    integer(ip)       :: status

#ifdef ENABLE_UMFPACK
    ! Fortran to C numbering 
    call decrement_array( matrix%gr%ia )
    call decrement_array( matrix%gr%ja )

    a_ => matrix%a(1,1,:)

    !
    ! Solve the linear system.
    !
    status = umfpack_di_solve ( UMFPACK_At, matrix%gr%ia, matrix%gr%ja, a_, sol, rhs, context%Numeric, context%Control, context%Info )

    if ( status < 0 ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'UMFPACK - Fatal error!'
      write ( *, '(a,i10)' ) '  umfpack_di_solve returns status = ', status
      check ( status == UMFPACK_OK ) 
    end if

    ! C to Fortran numbering 
    call increment_array( matrix%gr%ia )
    call increment_array( matrix%gr%ja )
#else
    call enable_umfpack_error_message
#endif

  end subroutine umfpack_solution_real

  !=============================================================================
  subroutine umfpack_solution_several_rhs ( context, matrix, nrhs, rhs, ldrhs, sol, ldsol )
    implicit none
    ! Parameters 
    type(umfpack_context), intent(inout), target :: context
    type(fem_matrix)  , intent(in)   , target :: matrix
    integer(ip)       , intent(in)            :: nrhs, ldrhs, ldsol
    real(rp)          , intent(in)   , target :: rhs (ldrhs, nrhs)
    real(rp)          , intent(inout), target :: sol (ldsol, nrhs)

    ! Locals
    integer(ip) :: i 
    real(rp), pointer :: a_(:)
    integer(ip)       :: status


#ifdef ENABLE_UMFPACK
    ! Fortran to C numbering 
    call decrement_array( matrix%gr%ia )
    call decrement_array( matrix%gr%ja )

    a_ => matrix%a(1,1,:)

    do i=1,nrhs
      !
      ! Solve the linear system.
      !
      status = umfpack_di_solve ( UMFPACK_At, matrix%gr%ia, matrix%gr%ja, a_, sol(1,i), rhs(1,i), context%Numeric, context%Control, context%Info )
      if ( status < 0 ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'UMFPACK - Fatal error!'
        write ( *, '(a,i10)' ) '  umfpack_di_solve returns status = ', status
        check ( status == UMFPACK_OK ) 
      end if
    end do

    ! C to Fortran numbering 
    call increment_array( matrix%gr%ia )
    call increment_array( matrix%gr%ja )
#else
    call enable_umfpack_error_message
#endif

  end subroutine umfpack_solution_several_rhs

end module umfpack_class


