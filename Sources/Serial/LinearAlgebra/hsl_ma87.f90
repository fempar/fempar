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
module hsl_ma87_names
  ! This module is a wrapper in which we define functions 
  ! to interact with HSL_MA87 using our data structures.
  ! Error control on calling parameters must (ideally)
  ! be performed here.

  ! Serial modules
  use types_names
  use memor_names
  use serial_scalar_matrix_names
  use serial_scalar_array_names 
  use graph_names
  use renumbering_names


#ifdef ENABLE_HSL_MA87
use hsl_ma87_double
#endif

# include "debug.i90"
  
  implicit none
  private

  ! Possible states of hsl_ma87_context
  integer(ip), parameter :: not_created   = 0
  integer(ip), parameter :: created       = 1 
  integer(ip), parameter :: symb_computed = 2 ! Symbolic data already computed
  integer(ip), parameter :: num_computed  = 3 ! Numerical data already computed 


  type hsl_ma87_context_t
     ! Our components
     integer(ip)              :: state = not_created
     type(renumbering_t)              :: renumbering
#ifdef ENABLE_HSL_MA87
     ! HSL_MA87 components
     type(MA87_keep) :: keep
#endif
  end type hsl_ma87_context_t

  type hsl_ma87_control_t
#ifdef ENABLE_HSL_MA87
     type(MA87_control) :: control
#endif
  end type hsl_ma87_control_t

  type hsl_ma87_info_t
#ifdef ENABLE_HSL_MA87
     type(MA87_info) :: info
#endif
  end type hsl_ma87_info_t

  ! Types
  public :: hsl_ma87_context_t, hsl_ma87_control_t, hsl_ma87_info_t

  ! Possible actions that can be perfomed by solve_hsl_ma87
  integer(ip), parameter :: hsl_ma87_init              = 1  ! Construct solve_hsl_ma87_state
  integer(ip), parameter :: hsl_ma87_finalize          = 2  ! Destruct  solve_hsl_ma87_state     
  integer(ip), parameter :: hsl_ma87_compute_symb      = 3  ! Compute symb. fact.  
  integer(ip), parameter :: hsl_ma87_compute_num       = 4  ! Compute numerical fact. 
  integer(ip), parameter :: hsl_ma87_solve             = 6  ! Fwd./Bck. substitution 
  integer(ip), parameter :: hsl_ma87_free_values       = 7
  integer(ip), parameter :: hsl_ma87_free_struct       = 8
  integer(ip), parameter :: hsl_ma87_free_clean        = 9

  interface hsl_ma87
     module procedure hsl_ma87_vector, hsl_ma87_r2, hsl_ma87_r1
  end interface hsl_ma87

  ! Constants (actions)
  public :: hsl_ma87_init,  &
       &    hsl_ma87_finalize, &
       &    hsl_ma87_compute_symb, & 
       &    hsl_ma87_compute_num, &
       &    hsl_ma87_solve, &
       &    hsl_ma87_free_values,&
       &    hsl_ma87_free_struct, &
       &    hsl_ma87_free_clean

  ! Functions
  public :: hsl_ma87

contains
  ! HSL_MA87 calling routines and their arguments
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
  subroutine hsl_ma87_vector ( action, context, A, b, x, ctrl, info )
    implicit none
    ! Mandatory Parameters
    type(hsl_ma87_context_t), intent(inout) :: context   ! Information required between calls
    integer(ip)           , intent(in)    :: action    ! Action to be performed 
                                                   ! (see public constants above)
    type(serial_scalar_matrix_t)  , intent(in)    :: A         ! Linear system coefficient matrix
    type(serial_scalar_array_t)  , intent(in)    :: b         ! RHS (Right-hand-side)
    type(serial_scalar_array_t)  , intent(inout) :: x         ! LHS (Left-hand-side)

    type(hsl_ma87_control_t) , intent(in)    :: ctrl
    type(hsl_ma87_info_t)    , intent(inout) :: info

    select case(action)

    case ( hsl_ma87_init )
       ! Check pre-conditions
       assert ( context%state == not_created )
       call hsl_ma87_ini ( context, A )
       ! State transition 
       context%state = created
    case ( hsl_ma87_finalize )

       ! Check pre-conditions
       assert ( context%state /= not_created )
       select case (context%state)
       case (created)
          call  hsl_ma87_free ( hsl_ma87_free_clean , context, ctrl )
       case (symb_computed)
          call  hsl_ma87_free ( hsl_ma87_free_struct, context, ctrl )
          call  hsl_ma87_free ( hsl_ma87_free_clean , context, ctrl )
       case (num_computed)
          call  hsl_ma87_free ( hsl_ma87_free_values, context, ctrl )
          call  hsl_ma87_free ( hsl_ma87_free_struct, context, ctrl )
          call  hsl_ma87_free ( hsl_ma87_free_clean , context, ctrl )
       end select
       ! State transition 
       context%state = not_created

    case ( hsl_ma87_compute_symb )

       ! Check pre-conditions
       assert (context%state==created)
       ! Create sparsity pattern of LU sparse direct factorization
       call hsl_ma87_analysis ( context, A, ctrl, info )
       ! State transition 
       context%state = symb_computed

    case ( hsl_ma87_compute_num )

       ! Check pre-conditions
       assert (context%state==symb_computed.or.context%state==num_computed)
       ! Compute LU sparse direct factorization
       call hsl_ma87_factorization ( context, A, ctrl, info )
       ! State transition 
       context%state = num_computed

    case ( hsl_ma87_solve )
       ! Check pre-conditions
       assert ( context%state ==  num_computed )
       call hsl_ma87_solution ( context, A, b, x, ctrl, info )

    case ( hsl_ma87_free_values )

       call  hsl_ma87_free ( hsl_ma87_free_values, context, ctrl )
       context%state=symb_computed

    case ( hsl_ma87_free_struct )

       call  hsl_ma87_free ( hsl_ma87_free_struct, context, ctrl )
       context%state=created

    case ( hsl_ma87_free_clean )

       call  hsl_ma87_free ( hsl_ma87_free_clean, context, ctrl )
       context%state=not_created

    case default

       ! Write an error message and stop ?      

    end select

  end subroutine hsl_ma87_vector
  !=============================================================================
  subroutine hsl_ma87_r2 ( action, context, A, nrhs, b, ldb, x, ldx, ctrl, info)
    implicit none
    ! Mandatory Parameters
    type(hsl_ma87_context_t), intent(inout) :: context   ! Information required between calls
    integer(ip)       , intent(in)    :: action    ! Action to be performed 
                                                   ! (see public constants above)
    type(serial_scalar_matrix_t)  , intent(in)    :: A         ! Linear system coefficient matrix
!    type(vector_t)  , intent(in)    :: b         ! RHS (Right-hand-side)
!    type(vector_t)  , intent(inout) :: x         ! LHS (Left-hand-side)
    integer(ip)       , intent(in)    :: nrhs, ldb, ldx
    real(rp)          , intent(in)    :: b (ldb, nrhs)
    real(rp)          , intent(inout) :: x (ldx, nrhs)

    type(hsl_ma87_control_t) , intent(in)          :: ctrl
    type(hsl_ma87_info_t)    , intent(inout)         :: info

    select case(action)

    case ( hsl_ma87_solve )
       
       ! Check pre-conditions
       assert ( context%state ==  num_computed )
       call hsl_ma87_solution_several_rhs ( context, A, nrhs, b, ldb, x, ldx, ctrl, info  )
       
    case default
       
       ! Write an error message and stop ?      

    end select

  end subroutine hsl_ma87_r2

  !=============================================================================
  subroutine hsl_ma87_r1 ( action, context, A, b, x, ctrl, info )
    implicit none
    ! Mandatory Parameters
    type(hsl_ma87_context_t), intent(inout) :: context   ! Information required between calls
    integer(ip)       , intent(in)    :: action    ! Action to be performed 
                                                   ! (see public constants above)
    type(serial_scalar_matrix_t)  , intent(in)    :: A         ! Linear system coefficient matrix
    real(rp)          , intent(in)    :: b (A%graph%nv)
    real(rp)          , intent(inout) :: x (A%graph%nv)

    type(hsl_ma87_control_t) , intent(in)          :: ctrl
    type(hsl_ma87_info_t)    , intent(inout)       :: info

    select case(action)

    case ( hsl_ma87_solve )

       ! Check pre-conditions
       assert ( context%state ==  num_computed )
       call hsl_ma87_solution_real ( context, A, b, x, ctrl, info  )

    case default

       ! Write an error message and stop ?      

    end select
    
  end subroutine hsl_ma87_r1

  !=============================================================================
#ifndef ENABLE_HSL_MA87
  subroutine enable_hsl_ma87_error_message
      implicit none
      write (0,*) 'Error: Fem was not compiled with -DENABLE_HSL_MA87.'
      write (0,*) "Error: You must activate this cpp macro in order to use HSL_MA87"
      check(1==0)
  end subroutine
#endif

  !=============================================================================
  subroutine hsl_ma87_ini ( context, matrix )
    implicit none
    ! Parameters
    type(hsl_ma87_context_t), intent(inout), target :: context
    type(serial_scalar_matrix_t)      , intent(in)            :: matrix

#ifdef ENABLE_HSL_MA87

#else
    call enable_hsl_ma87_error_message
#endif
  end subroutine hsl_ma87_ini

  !=============================================================================
  subroutine hsl_ma87_free ( mode, context, ctrl )
    implicit none
    ! Parameters
    integer(ip)            , intent(in)    :: mode
    type(hsl_ma87_context_t) , intent(inout) :: context
    type(hsl_ma87_control_t) , intent(in)    :: ctrl

#ifdef ENABLE_HSL_MA87
    ! Free hsl_ma87_context structures
    if ( mode == hsl_ma87_free_clean ) then

    else if ( mode == hsl_ma87_free_struct  ) then

    else if ( mode == hsl_ma87_free_values ) then
       ! ** IMPORTANT NOTE: for non-linear iterations/transient problems
       !    it may have more sense to call
       !    ma87_finalize on the free_clean
       !    section instead of here. TO THINK.
       call renumbering_free( context%renumbering )
       call ma87_finalise( context%keep, ctrl%control )
    end if
#else
    call enable_hsl_ma87_error_message
#endif

  end subroutine hsl_ma87_free

  !=============================================================================
  subroutine hsl_ma87_analysis ( context, matrix, ctrl, info )
use partitioning_params_names
use graph_renumbering_names
    implicit none

    ! Parameters 
    type(hsl_ma87_context_t) , intent(inout) :: context
    type(serial_scalar_matrix_t)       , intent(in)    :: matrix
    type(hsl_ma87_control_t) , intent(in)    :: ctrl
    type(hsl_ma87_info_t)    , intent(out)   :: info

    ! Locals (required for the call to graph_nd_renumbering)
    type (graph_t)                  :: aux_graph
    type(partitioning_params_t)             :: prt_parts
    integer(ip)                   :: i

    assert ( matrix%graph%symmetric_storage )

#ifdef ENABLE_HSL_MA87
    call renumbering_alloc( matrix%graph%nv, context%renumbering )

    ! Call to graph_nd_renumbering
    ! Set-up graph
    aux_graph%symmetric_storage = .false.
    aux_graph%nv   = matrix%graph%nv
    aux_graph%nv2  = matrix%graph%nv2
    call memalloc (aux_graph%nv+1, aux_graph%ia, __FILE__,__LINE__)
    call memalloc ( (matrix%graph%ia(aux_graph%nv+1)-1-aux_graph%nv)*2, aux_graph%ja, __FILE__,__LINE__)

    call half_to_full ( aux_graph%nv, matrix%graph%ia, matrix%graph%ja, aux_graph%ia, aux_graph%ja )

    ! call graph_print (6, matrix%graph)
    ! call graph_print (6, aux_graph)
    
    call graph_nd_renumbering(prt_parts,aux_graph,context%renumbering) 

    ! De-allocate aux_graph
    call aux_graph%free()

    call ma87_analyse(matrix%graph%nv, matrix%graph%ia, matrix%graph%ja, &
                      context%renumbering%lperm, context%keep, ctrl%control, info%info)

    if ( info%info%flag < 0 ) then
       write (0,*) 'Error, HSL_MA87: the following ERROR was triggered: ', & 
                   info%info%flag, 'during ma87_analyse'
       check ( info%info%flag >= 0 )
    else if (info%info%flag > 0 ) then
       write (0,*) 'Warning, HSL_MA87: the following WARNING was triggered: ', & 
                   info%info%flag, 'during ma87_analyse'
    end if
#else
    call enable_hsl_ma87_error_message
#endif


  end subroutine hsl_ma87_analysis


  ! Convert a matrix in half storage to one in full storage.
  ! Drops any diagonal entries.
  subroutine half_to_full(n, ptr, row, ptr2, row2)
    implicit none
    integer(ip), intent(in) :: n
    integer(ip), dimension(n+1), intent(in) :: ptr
    integer(ip), dimension(ptr(n+1)-1), intent(in) :: row
    integer(ip), dimension(*), intent(out) :: ptr2
    integer(ip), dimension(*), intent(out) :: row2

    integer(ip) :: i, j, k

    ! Set ptr2(j) to hold no. nonzeros in column j
    ptr2(1:n+1) = 0
    DO j = 1, n
       DO k = ptr(j), ptr(j+1) - 1
          i = row(k)
          IF (j/=i) THEN
             ptr2(i) = ptr2(i) + 1
             ptr2(j) = ptr2(j) + 1
          END IF
       END DO
    END DO

    ! Set ptr2(j) to point to where row indices will end in row2
    DO j = 2, n
       ptr2(j) = ptr2(j-1) + ptr2(j)
    END DO
    ptr2(n+1) = ptr2(n) + 1

    ! Fill ptr2 and row2
    DO j = 1, n
       DO k = ptr(j), ptr(j+1) - 1
          i = row(k)
          IF (j/=i) THEN
             row2(ptr2(i)) = j
             row2(ptr2(j)) = i
             ptr2(i) = ptr2(i) - 1
             ptr2(j) = ptr2(j) - 1
          END IF
       END DO
    END DO
    DO j = 1, n
       ptr2(j) = ptr2(j) + 1
    END DO

  end subroutine half_to_full

  !=============================================================================
  subroutine hsl_ma87_factorization ( context, matrix, ctrl, info )
    implicit none
    ! Parameters 
    type(hsl_ma87_context_t) , intent(inout) :: context
    type(serial_scalar_matrix_t)       , target, intent(in) :: matrix
    type(hsl_ma87_control_t) , intent(in)    :: ctrl
    type(hsl_ma87_info_t)    , intent(inout) :: info

    ! Locals
    real(rp), pointer :: a_(:)

    assert ( matrix%is_symmetric .and. matrix%sign == positive_definite .and. matrix%graph%symmetric_storage )

#ifdef ENABLE_HSL_MA87
    a_ => matrix%a(:)
   
    ! Factor
    call ma87_factor(matrix%graph%nv, matrix%graph%ia, matrix%graph%ja, &
                     a_, context%renumbering%lperm, context%keep, ctrl%control, info%info)


    if ( info%info%flag < 0 ) then
       write (0,*) 'Error, HSL_MA87: the following ERROR was triggered: ', & 
                   info%info%flag, 'during ma87_factor'
       check( info%info%flag >= 0 )
    else if (info%info%flag > 0 ) then
       write (0,*) 'Warning, HSL_MA87: the following WARNING was triggered: ', & 
                   info%info%flag, 'during ma87_factor'
    end if
#else
    call enable_hsl_ma87_error_message
#endif

  end subroutine hsl_ma87_factorization

  !=============================================================================
  subroutine hsl_ma87_solution ( context, matrix, x, y, ctrl, info )
    ! Computes y <- A^-1 * x, using previously computed LU factorization
    implicit none
    ! Parameters 
    type(hsl_ma87_context_t), intent(inout)        :: context
    type(serial_scalar_matrix_t)  , intent(in)               :: matrix
    type(serial_scalar_array_t)  , intent(in), target       :: x
    type(serial_scalar_array_t)  , intent(inout), target    :: y
    type(hsl_ma87_control_t) , intent(in)          :: ctrl
    type(hsl_ma87_info_t)    , intent(inout)       :: info

    ! Locals
    real(rp), pointer :: x_(:)
    real(rp), pointer :: y_(:)

#ifdef ENABLE_HSL_MA87
    x_ => x%b(:)
    y_ => y%b(:)
    
    y_ = x_ 
    
    ! Solve
    call ma87_solve(y_, context%renumbering%lperm, context%keep, ctrl%control, info%info)

    if ( info%info%flag < 0 ) then
       write (0,*) 'Error, HSL_MA87: the following ERROR was triggered: ', & 
                   info%info%flag, 'during ma87_solve'
       check ( info%info%flag >= 0 )
    else if (info%info%flag > 0 ) then
       write (0,*) 'Warning, HSL_MA87: the following WARNING was triggered: ', & 
                   info%info%flag, 'during ma87_solve'
    end if
#else
    call enable_hsl_ma87_error_message
#endif

  end subroutine hsl_ma87_solution

  !=============================================================================
  subroutine hsl_ma87_solution_real ( context, matrix, rhs, sol, ctrl, info )
    implicit none
    ! Parameters 
    type(hsl_ma87_context_t), intent(inout) :: context
    type(serial_scalar_matrix_t)  , intent(in)   :: matrix
    real(rp)          , intent(in)    :: rhs (matrix%graph%nv)
    real(rp)          , intent(inout) :: sol (matrix%graph%nv)
    type(hsl_ma87_control_t) , intent(in) :: ctrl
    type(hsl_ma87_info_t)    , intent(inout):: info

#ifdef ENABLE_HSL_MA87
    sol = rhs
    ! Solve
    call ma87_solve(sol, context%renumbering%lperm, context%keep, ctrl%control, info%info)

    if ( info%info%flag < 0 ) then
       write (0,*) 'Error, HSL_MA87: the following ERROR was triggered: ', & 
                   info%info%flag, 'during ma87_solve'
       check ( info%info%flag >= 0 )
    else if (info%info%flag > 0 ) then
       write (0,*) 'Warning, HSL_MA87: the following WARNING was triggered: ', & 
                   info%info%flag, 'during ma87_solve'
    end if
#else
    call enable_hsl_ma87_error_message
#endif

  end subroutine hsl_ma87_solution_real

  !=============================================================================
  subroutine hsl_ma87_solution_several_rhs ( context, matrix, nrhs, rhs, ldrhs, sol, ldsol, ctrl, info )
    implicit none
    ! Parameters 
    type(hsl_ma87_context_t), intent(inout), target :: context
    type(serial_scalar_matrix_t)  , intent(in)   , target :: matrix
    integer(ip)       , intent(in)            :: nrhs, ldrhs, ldsol
    real(rp)          , intent(in)   , target :: rhs (ldrhs, nrhs)
    real(rp)          , intent(inout), target :: sol (ldsol, nrhs)

    type(hsl_ma87_control_t) , intent(in)          :: ctrl
    type(hsl_ma87_info_t)    , intent(out)         :: info

    ! Locals
    integer(ip) :: i 

#ifdef ENABLE_HSL_MA87
    do i=1,nrhs
       sol(1:matrix%graph%nv,i) = rhs(1:matrix%graph%nv,i)
    end do

    ! Solve
    call ma87_solve(nrhs, ldsol, sol, context%renumbering%lperm, context%keep, ctrl%control, info%info)

    if ( info%info%flag < 0 ) then
       write (0,*) 'Error, HSL_MA87: the following ERROR was triggered: ', & 
                   info%info%flag, 'during ma87_solve'
       check ( info%info%flag >= 0 )
    else if (info%info%flag > 0 ) then
       write (0,*) 'Warning, HSL_MA87: the following WARNING was triggered: ', & 
                   info%info%flag, 'during ma87_solve'
    end if
#else
    call enable_hsl_ma87_error_message
#endif

  end subroutine hsl_ma87_solution_several_rhs

end module hsl_ma87_names


