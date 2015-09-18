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
module hsl_mi20_names
  ! This module is a wrapper in which we define functions 
  ! to interact with HSL_MI20 using our data structures.
  ! Error control on calling parameters must (ideally)
  ! be performed here.

  ! Serial modules
  use types_names
  use memor_names
  use serial_scalar_matrix_names
  use serial_scalar_array_names 
  use graph_names

#ifdef ENABLE_HSL_MI20
  use hsl_mi20_double
#endif

# include "debug.i90"
  
  implicit none
  private

  ! Possible states of hsl_mi20_context
  integer(ip), parameter :: not_created   = 0
  integer(ip), parameter :: created       = 1 
  integer(ip), parameter :: symb_computed = 2 ! Symbolic data already computed
  integer(ip), parameter :: num_computed  = 3 ! Numerical data already computed 


  type hsl_mi20_context_t
     ! Our components
     integer(ip)     :: state = not_created
#ifdef ENABLE_HSL_MI20
     ! HSL_MI20 components
     type(zd11_type) :: zd11_mat
     type(mi20_keep) :: keep
#endif
  end type hsl_mi20_context_t

  type hsl_mi20_control_t
#ifdef ENABLE_HSL_MI20
     type(MI20_control) :: control
#endif
  end type hsl_mi20_control_t

  type hsl_mi20_info_t
#ifdef ENABLE_HSL_MI20
     type(MI20_info) :: info
#endif
  end type hsl_mi20_info_t

  type hsl_mi20_data_t
#ifdef ENABLE_HSL_MI20
     type(MI20_data), allocatable :: coarse_data(:) 
#endif     
  end type hsl_mi20_data_t

  ! Types
  public :: hsl_mi20_context_t, hsl_mi20_control_t, hsl_mi20_info_t, hsl_mi20_data_t

  ! Possible actions that can be perfomed by solve_hsl_mi20
  integer(ip), parameter :: hsl_mi20_init              = 1  ! Construct solve_hsl_mi20_state
  integer(ip), parameter :: hsl_mi20_finalize          = 2  ! Destruct  solve_hsl_mi20_state     
  integer(ip), parameter :: hsl_mi20_compute_symb      = 3  ! Compute symb. fact.  
  integer(ip), parameter :: hsl_mi20_compute_num       = 4  ! Compute numerical fact. 
  integer(ip), parameter :: hsl_mi20_solve             = 6  ! Fwd./Bck. substitution 
  integer(ip), parameter :: hsl_mi20_free_values       = 7
  integer(ip), parameter :: hsl_mi20_free_struct       = 8
  integer(ip), parameter :: hsl_mi20_free_clean        = 9

  interface hsl_mi20
     module procedure hsl_mi20_vector, hsl_mi20_r2, hsl_mi20_r1
  end interface hsl_mi20

  ! Constants (actions)
  public :: hsl_mi20_init,  &
       &    hsl_mi20_finalize, &
       &    hsl_mi20_compute_symb, & 
       &    hsl_mi20_compute_num, &
       &    hsl_mi20_solve, &
       &    hsl_mi20_free_values,&
       &    hsl_mi20_free_struct, &
       &    hsl_mi20_free_clean

  ! Functions
  public :: hsl_mi20

contains
  ! HSL_MI20 calling routines and their arguments
  !
  ! call wgsmp (n, ia, ja, avals, b, ldb, nrhs, rmisc, iparm, dparm)      
  !
  ! n, ia, ja, avals define the matrix
  ! ldb nrhs b(ldb,nrhs) define the vector 
  ! rmisc(n,nrhs) is the backward error, only used if iparm(25)=1 in
  !               the triangular solution and iterative refinement.
  !
  ! call wssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
  !             aux, naux, mrp, iparm, dparm)
  !
  ! n, ia, ja, avals define the matrix
  ! ldb nrhs b(ldb,nrhs) define the vector 
  ! perm iperm
  ! aux obsolete not used
  ! naux = 0 not used
  ! mrp(n) pivot information, only used if iparm(11)=2 in the factorization
  !
  !=============================================================================
  subroutine hsl_mi20_vector ( action, context, A, b, x, data, ctrl, info )
    implicit none
    ! Mandatory Parameters
    type(hsl_mi20_context_t), intent(inout) :: context   ! Information required between calls
    integer(ip)           , intent(in)    :: action    ! Action to be performed 
                                                   ! (see public constants above)
    type(serial_scalar_matrix_t)  , intent(in)    :: A         ! Linear system coefficient matrix
    type(serial_scalar_array_t)  , intent(in)    :: b         ! RHS (Right-hand-side)
    type(serial_scalar_array_t)  , intent(inout) :: x         ! LHS (Left-hand-side)

    type(hsl_mi20_data_t)    , intent(inout) :: data
    type(hsl_mi20_control_t) , intent(in)    :: ctrl
    type(hsl_mi20_info_t)    , intent(inout) :: info

    select case(action)

    case ( hsl_mi20_init )

       ! Check pre-conditions
       assert ( context%state == not_created )
       call hsl_mi20_ini ( context, A )
       ! State transition 
       context%state = created
    case ( hsl_mi20_finalize )

       ! Check pre-conditions
       assert ( context%state /= not_created )
       select case (context%state)
       case (created)
          call  hsl_mi20_free ( hsl_mi20_free_clean , context, data, ctrl, info )
       case (symb_computed)
          call  hsl_mi20_free ( hsl_mi20_free_struct, context, data, ctrl, info )
          call  hsl_mi20_free ( hsl_mi20_free_clean , context, data, ctrl, info )
       case (num_computed)
          call  hsl_mi20_free ( hsl_mi20_free_values, context, data, ctrl, info )
          call  hsl_mi20_free ( hsl_mi20_free_struct, context, data, ctrl, info )
          call  hsl_mi20_free ( hsl_mi20_free_clean , context, data, ctrl, info )
       end select
       ! State transition 
       context%state = not_created

    case ( hsl_mi20_compute_symb )

       ! Check pre-conditions
       assert (context%state==created)
       ! Create sparsity pattern of LU sparse direct factorization
       call hsl_mi20_analysis ( context, A )
       ! State transition 
       context%state = symb_computed

    case ( hsl_mi20_compute_num )

       ! Check pre-conditions
       assert (context%state==symb_computed.or.context%state==num_computed)
       ! Compute LU sparse direct factorization
       call hsl_mi20_factorization ( context, A, data, ctrl, info )
       ! State transition 
       context%state = num_computed

    case ( hsl_mi20_solve )
       ! Check pre-conditions
       assert ( context%state ==  num_computed )
       call hsl_mi20_solution ( context, A, b, x, data, ctrl, info )

    case ( hsl_mi20_free_values )

       call  hsl_mi20_free ( hsl_mi20_free_values, context, data, ctrl, info )
       context%state=symb_computed

    case ( hsl_mi20_free_struct )

       call  hsl_mi20_free ( hsl_mi20_free_struct, context, data, ctrl, info )
       context%state=created

    case ( hsl_mi20_free_clean )

       call  hsl_mi20_free ( hsl_mi20_free_clean, context, data, ctrl, info )
       context%state=not_created

    case default

       ! Write an error message and stop ?      

    end select

  end subroutine hsl_mi20_vector
  !=============================================================================
  subroutine hsl_mi20_r2 ( action, context, A, nrhs, b, ldb, x, ldx, data, ctrl, info)
    implicit none
    ! Mandatory Parameters
    type(hsl_mi20_context_t), intent(inout) :: context   ! Information required between calls
    integer(ip)       , intent(in)    :: action    ! Action to be performed 
                                                   ! (see public constants above)
    type(serial_scalar_matrix_t)  , intent(in)    :: A         ! Linear system coefficient matrix
!    type(vector_t)  , intent(in)    :: b         ! RHS (Right-hand-side)
!    type(vector_t)  , intent(inout) :: x         ! LHS (Left-hand-side)
    integer(ip)       , intent(in)    :: nrhs, ldb, ldx
    real(rp)          , intent(in)    :: b (ldb, nrhs)
    real(rp)          , intent(inout) :: x (ldx, nrhs)

    type(hsl_mi20_data_t)    , intent(in)          :: data
    type(hsl_mi20_control_t) , intent(in)          :: ctrl
    type(hsl_mi20_info_t)    , intent(out)         :: info

    select case(action)

    case ( hsl_mi20_solve )

       ! Check pre-conditions
       assert ( context%state ==  num_computed )
       call hsl_mi20_solution_several_rhs ( context, A, nrhs, b, ldb, x, ldx, data, ctrl, info  )

    case default

       ! Write an error message and stop ?      

    end select

  end subroutine hsl_mi20_r2

  !=============================================================================
  subroutine hsl_mi20_r1 ( action, context, A, b, x, data, ctrl, info )
    implicit none
    ! Mandatory Parameters
    type(hsl_mi20_context_t), intent(inout) :: context   ! Information required between calls
    integer(ip)       , intent(in)    :: action    ! Action to be performed 
                                                   ! (see public constants above)
    type(serial_scalar_matrix_t)  , intent(in)    :: A         ! Linear system coefficient matrix
    real(rp)          , intent(in)    :: b (A%graph%nv)
    real(rp)          , intent(inout) :: x (A%graph%nv)

    type(hsl_mi20_data_t)    , intent(in)          :: data
    type(hsl_mi20_control_t) , intent(in)          :: ctrl
    type(hsl_mi20_info_t)    , intent(out)         :: info

    select case(action)

    case ( hsl_mi20_solve )

       ! Check pre-conditions
       assert ( context%state ==  num_computed )
       call hsl_mi20_solution_real ( context, A, b, x, data, ctrl, info  )

    case default

       ! Write an error message and stop ?      

    end select
    
  end subroutine hsl_mi20_r1

  !=============================================================================
#ifndef ENABLE_HSL_MI20
  subroutine enable_hsl_mi20_error_message
      implicit none
      write (0,*) 'Error: Fem was not compiled with -DENABLE_HSL_MI20.'
      write (0,*) "Error: You must activate this cpp macro in order to use HSL_MI20"
      stop
  end subroutine
#endif

  !=============================================================================
  subroutine hsl_mi20_ini ( context, matrix )
    implicit none
    ! Parameters
    type(hsl_mi20_context_t), intent(inout), target :: context
    type(serial_scalar_matrix_t)      , intent(in)            :: matrix

#ifdef ENABLE_HSL_MI20
#else
    call enable_hsl_mi20_error_message
#endif
  end subroutine hsl_mi20_ini

  !=============================================================================
  subroutine hsl_mi20_free ( mode, context, data, ctrl, info )
    implicit none
    ! Parameters
    integer(ip)            , intent(in)    :: mode
    type(hsl_mi20_context_t) , intent(inout) :: context
    type(hsl_mi20_data_t)    , intent(inout) :: data
    type(hsl_mi20_control_t) , intent(in)    :: ctrl
    type(hsl_mi20_info_t)    , intent(out)   :: info

#ifdef ENABLE_HSL_MI20
    ! Free hsl_mi20_context structures
    if ( mode == hsl_mi20_free_clean ) then

    else if ( mode == hsl_mi20_free_struct  ) then

    else if ( mode == hsl_mi20_free_values ) then
       call memfree( context%zd11_mat%ptr,__FILE__,__LINE__)

       call memfree( context%zd11_mat%col,__FILE__,__LINE__)
       
       call memfree( context%zd11_mat%val,__FILE__,__LINE__)

       call mi20_finalize( data%coarse_data, & 
                           context%keep, & 
                           ctrl%control, & 
                           info%info)
    end if
#else
    call enable_hsl_mi20_error_message
#endif

  end subroutine hsl_mi20_free

  !=============================================================================
  subroutine hsl_mi20_analysis ( context, matrix )
    implicit none
    ! Parameters 
    type(hsl_mi20_context_t), intent(inout) :: context
    type(serial_scalar_matrix_t)      , intent(in)    :: matrix

#ifdef ENABLE_HSL_MI20
#else
    call enable_hsl_mi20_error_message
#endif

  end subroutine hsl_mi20_analysis

  !=============================================================================
  subroutine hsl_mi20_factorization ( context, matrix, data, ctrl, info )
    implicit none
    ! Parameters 
    type(hsl_mi20_context_t) , intent(inout) :: context
    type(serial_scalar_matrix_t)       , intent(in)    :: matrix
    type(hsl_mi20_data_t)    , intent(out)   :: data
    type(hsl_mi20_control_t) , intent(in) :: ctrl
    type(hsl_mi20_info_t)    , intent(out)   :: info

#ifdef ENABLE_HSL_MI20
    assert ( .not. matrix%graph%symmetric_storage )
    
    ! Copy matrix_t to type(ZD11_type)
    context%zd11_mat%m = matrix%graph%nv
    context%zd11_mat%n = matrix%graph%nv2

    call memalloc( context%zd11_mat%m+1,context%zd11_mat%ptr,__FILE__,__LINE__)

    context%zd11_mat%ptr = matrix%graph%ia

    call memalloc( context%zd11_mat%ptr(context%zd11_mat%m+1)-1,context%zd11_mat%col, __FILE__,__LINE__)

    context%zd11_mat%col = matrix%graph%ja

    call memalloc( context%zd11_mat%ptr(context%zd11_mat%m+1)-1, context%zd11_mat%val,__FILE__,__LINE__)

    context%zd11_mat%val = matrix%a(:)

    ! Setup AMG preconditioner
    call mi20_setup( context%zd11_mat, & 
                     data%coarse_data, & 
                     context%keep, &   
                     ctrl%control, &
                     info%info )

    if ( info%info%flag < 0 ) then
       write (0,*) 'Error, HSL_MI20: the following ERROR was triggered: ', & 
                   info%info%flag, 'during mi20_setup_stage'
       stop
    else if (info%info%flag > 0 ) then
       write (0,*) 'Warning, HSL_MI20: the following WARNING was triggered: ', & 
                   info%info%flag, 'during mi20_setup_stage'
    end if
#else
    call enable_hsl_mi20_error_message
#endif

  end subroutine hsl_mi20_factorization

  !=============================================================================
  subroutine hsl_mi20_solution ( context, matrix, x, y, data, ctrl, info )
    ! Computes y <- A^-1 * x, using previously computed LU factorization
    implicit none
    ! Parameters 
    type(hsl_mi20_context_t), intent(inout)        :: context
    type(serial_scalar_matrix_t)  , intent(in)               :: matrix
    type(serial_scalar_array_t)  , intent(in), target       :: x
    type(serial_scalar_array_t)  , intent(inout), target    :: y
    type(hsl_mi20_data_t)    , intent(in)          :: data
    type(hsl_mi20_control_t) , intent(in)          :: ctrl
    type(hsl_mi20_info_t)    , intent(out)         :: info

    ! Locals
    real(rp), pointer :: x_(:)
    real(rp), pointer :: y_(:)

#ifdef ENABLE_HSL_MI20
    x_ => x%b(:)
    y_ => y%b(:)
    
    ! Solve with AMG preconditioner
    call mi20_precondition( context%zd11_mat, & 
                            data%coarse_data, & 
                            x_, &   
                            y_, &
                            context%keep, &
                            ctrl%control, &
                            info%info ) 

    if ( info%info%flag < 0 ) then
       write (0,*) 'Error, HSL_MI20: the following ERROR was triggered: ', & 
                   info%info%flag, 'during mi20_setup_stage'
       stop
    else if (info%info%flag > 0 ) then
       write (0,*) 'Warning, HSL_MI20: the following WARNING was triggered: ', & 
                   info%info%flag, 'during mi20_setup_stage'
    end if
#else
    call enable_hsl_mi20_error_message
#endif

  end subroutine hsl_mi20_solution

  !=============================================================================
  subroutine hsl_mi20_solution_real ( context, matrix, rhs, sol, data, ctrl, info )
    implicit none
    ! Parameters 
    type(hsl_mi20_context_t), intent(inout) :: context
    type(serial_scalar_matrix_t)  , intent(in)   :: matrix
    real(rp)          , intent(in)    :: rhs (matrix%graph%nv)
    real(rp)          , intent(inout) :: sol (matrix%graph%nv)
    type(hsl_mi20_data_t)    , intent(in) :: data
    type(hsl_mi20_control_t) , intent(in) :: ctrl
    type(hsl_mi20_info_t)    , intent(out):: info

#ifdef ENABLE_HSL_MI20
    ! Solve with AMG preconditioner
    call mi20_precondition( context%zd11_mat, & 
                            data%coarse_data, & 
                            rhs, &   
                            sol, &
                            context%keep, &
                            ctrl%control, &
                            info%info ) 
#else
    call enable_hsl_mi20_error_message
#endif

  end subroutine hsl_mi20_solution_real

  !=============================================================================
  subroutine hsl_mi20_solution_several_rhs ( context, matrix, nrhs, rhs, ldrhs, sol, ldsol, data, ctrl, info )
    implicit none
    ! Parameters 
    type(hsl_mi20_context_t), intent(inout), target :: context
    type(serial_scalar_matrix_t)  , intent(in)   , target :: matrix
    integer(ip)       , intent(in)            :: nrhs, ldrhs, ldsol
    real(rp)          , intent(in)   , target :: rhs (ldrhs, nrhs)
    real(rp)          , intent(inout), target :: sol (ldsol, nrhs)

    type(hsl_mi20_data_t)    , intent(in)          :: data
    type(hsl_mi20_control_t) , intent(in)          :: ctrl
    type(hsl_mi20_info_t)    , intent(out)         :: info

    ! Locals
    integer(ip)       :: irhs
    real(rp), pointer :: rhs_(:)
    real(rp), pointer :: sol_(:)
    
#ifdef ENABLE_HSL_MI20
    do irhs=1, nrhs
       rhs_ => rhs(1:matrix%graph%nv,irhs)
       sol_ => sol(1:matrix%graph%nv,irhs)
       call hsl_mi20_solution_real ( context, matrix, rhs_, sol_, data, ctrl, info  )
    end do
#else
    call enable_hsl_mi20_error_message
#endif

  end subroutine hsl_mi20_solution_several_rhs

end module hsl_mi20_names


