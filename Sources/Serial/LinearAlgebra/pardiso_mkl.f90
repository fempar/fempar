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
#ifdef ENABLE_MKL
  include 'mkl_pardiso.f90'
#endif
module pardiso_mkl_names

  ! This module is a wrapper in which we define functions 
  ! to interact with Pardiso using our data structures.
  ! Error control on calling parameters must (ideally)
  ! be performed here.

  ! Serial modules
  use types_names
  use memor_names
  use serial_scalar_matrix_names
  use serial_scalar_array_names 

  ! F90 interface to Intel MKL PARDISO
#ifdef ENABLE_MKL
  use mkl_pardiso
#endif
  ! Temporary solution for using more than one wrapper
  use mkl_names
  
  implicit none
# include "debug.i90"

  private

  ! Possible states of pardiso_mkl_instance
  integer(ip), parameter :: not_created   = 0
  integer(ip), parameter :: created       = 1 
  integer(ip), parameter :: symb_computed = 2 ! Symbolic data already computed
  integer(ip), parameter :: num_computed  = 3 ! Numerical data already computed  

  ! Derived data type which stores information required
  ! between calls. Initialization of state to not_created 
  ! is a MUST for this module to work correctly.
  type pardiso_mkl_context_t
     ! Our components
     integer(ip) :: state = not_created
     ! Pardiso components (send as arguments)
     integer     :: mtype   
     integer     :: n       
#ifdef ENABLE_MKL
     type ( mkl_pardiso_handle ), allocatable  :: pt(:)
#endif
  end type pardiso_mkl_context_t

  ! Types
  public :: pardiso_mkl_context_t

  ! *** Begin Pardiso MKL interface constants *** !
  ! The following constant is always declared independently
  ! of the macro ENABLE_MKL as it is required by
  ! the definition of the headers of some subroutines
  integer , parameter :: dp      = 8 ! kind (1.0D0)

  ! The following ones are always declared independently of the
  ! macro ENABLE_MKL as they are required by other modules
  ! that use the present module
  integer , parameter    :: pardiso_mkl_spd =  2  ! Real Symmetric positive definite 
  integer , parameter    :: pardiso_mkl_sin = -2  ! Real Symmetric indefinite
  integer , parameter    :: pardiso_mkl_uss = 1   ! Real Unsymmetric, structurally symmetric
  integer , parameter    :: pardiso_mkl_uns = 11  ! Real Unsymmetric, structurally unsymmetric

  ! Possible actions that can be perfomed by pardiso_mkl_solver
  integer(ip), parameter :: pardiso_mkl_initialize        = 1  ! Construct pardiso_mkl_solver
  integer(ip), parameter :: pardiso_mkl_finalize          = 2  ! Destruct  pardiso_mkl_solver     
  integer(ip), parameter :: pardiso_mkl_compute_symb      = 3  ! Compute symb. fact.  
  integer(ip), parameter :: pardiso_mkl_compute_num       = 4  ! Compute numerical fact. 
  integer(ip), parameter :: pardiso_mkl_compute_symb_num  = 5  ! Compute symb.+num. fact.
  integer(ip), parameter :: pardiso_mkl_solve             = 6  ! Fwd./Bck. substitution 

  interface pardiso_mkl
     module procedure pardiso_mkl_vector, pardiso_mkl_r2, pardiso_mkl_r1
  end interface pardiso_mkl

  ! Constants (matrix types)
  public :: pardiso_mkl_spd, &
       &    pardiso_mkl_sin, &
       &    pardiso_mkl_uss, &
       &    pardiso_mkl_uns
  ! Constants (actions)
  public :: pardiso_mkl_initialize, &
       &    pardiso_mkl_finalize, &
       &    pardiso_mkl_compute_symb, & 
       &    pardiso_mkl_compute_num, &
       &    pardiso_mkl_compute_symb_num, &
       &    pardiso_mkl_solve

  ! Functions
  public :: pardiso_mkl ! , pardiso_mkl_solution_several_rhs

#ifdef ENABLE_MKL
  integer , parameter :: maxfct  = 1
  integer , parameter :: mnum    = 1
  integer , parameter :: nrhsp   = 1
#endif
  ! *** End Pardiso MKL interface constants *** !

contains

  !=============================================================================
  subroutine pardiso_mkl_vector ( action, context, A, b, x, iparm, msglvl, perm )
    implicit none

    ! Mandatory Parameters
    integer(ip)              , intent(in)    :: action  ! Action to be performed
                                                        ! (see public constants above)
    type(pardiso_mkl_context_t), intent(inout) :: context ! Information required between calls
    type(serial_scalar_matrix_t), intent(in)             :: A       ! Linear system coefficient matrix
    type(serial_scalar_array_t), intent(in)             :: b       ! RHS (Right-hand-side)
    type(serial_scalar_array_t), intent(inout)          :: x       ! LHS (Left-hand-side)
    integer         , intent(inout), target, optional :: iparm(64)
    integer         , optional                     :: msglvl
    integer         , intent(in), target, optional :: perm(*)
!    integer         , intent(in), target, optional :: perm(A%graph%nv)

    select case(action)
    case ( pardiso_mkl_initialize )
       ! Check pre-conditions
       assert ( context%state == not_created )
       ! Create a context and fill default parameters in iparm (if present)
       call pardiso_mkl_init ( context, A, iparm)
       ! State transition 
       context%state = created
    case ( pardiso_mkl_finalize )
       ! Check pre-conditions
       assert ( context%state /= not_created )
       select case (context%state)
       case (created)
          call  pardiso_mkl_free ( free_clean , context, iparm, msglvl )
       case (symb_computed)
          call  pardiso_mkl_free ( free_symbolic_setup, context, iparm, msglvl )
          call  pardiso_mkl_free ( free_clean , context, iparm, msglvl )
       case (num_computed)
          call  pardiso_mkl_free ( free_numerical_setup, context, iparm, msglvl )
          call  pardiso_mkl_free ( free_symbolic_setup, context, iparm, msglvl )
          call  pardiso_mkl_free ( free_clean , context, iparm, msglvl )
       end select
       ! State transition 
       context%state = not_created
    case ( pardiso_mkl_compute_symb )
       ! Check pre-conditions
       assert (context%state==created.or.context%state==symb_computed.or.context%state==num_computed)
       ! Create sparsity pattern of LU sparse direct factorization
       select case ( context%state )
       case (created)       ! Compute symbolic info
          call pardiso_mkl_analysis ( context, A, iparm, msglvl, perm)
       case (symb_computed) ! Recompute symbolic info
          call pardiso_mkl_free     ( free_symbolic_setup, context, iparm, msglvl )
          call pardiso_mkl_analysis ( context, A, iparm, msglvl, perm)
       case (num_computed)  ! Recompute symbolic info
          call pardiso_mkl_free     ( free_numerical_setup, context, iparm, msglvl )
          call pardiso_mkl_free     ( free_symbolic_setup, context, iparm, msglvl )
          call pardiso_mkl_analysis ( context, A, iparm, msglvl, perm)
       end select
       ! State transition 
       context%state = symb_computed
    case ( pardiso_mkl_compute_num )
       ! Check pre-conditions
       assert (context%state==symb_computed.or.context%state==num_computed)
       ! Compute LU sparse direct factorization
       select case (context%state)
       case (symb_computed) ! Compute numerical info
          call pardiso_mkl_factorization ( context, A, iparm, msglvl)             
       case (num_computed)  ! Recompute numerical info
          call pardiso_mkl_free          ( free_numerical_setup, context, iparm, msglvl )
          call pardiso_mkl_factorization ( context, A, iparm, msglvl)
       end select
       ! State transition 
       context%state = num_computed
    case ( pardiso_mkl_compute_symb_num )
       ! Check pre-conditions
       assert (context%state==created.or.context%state==symb_computed.or.context%state==num_computed)
       select case (context%state)
       case (created)       ! Compute symbolic info
          call pardiso_mkl_analysis      ( context, A, iparm, msglvl, perm)
          call pardiso_mkl_factorization ( context, A, iparm, msglvl )
       case (symb_computed) ! Recompute symbolic info
          call pardiso_mkl_free          ( free_symbolic_setup, context, iparm, msglvl)
          call pardiso_mkl_analysis      ( context, A, iparm, msglvl, perm)
          call pardiso_mkl_factorization ( context, A, iparm, msglvl )
       case (num_computed)  ! Recompute symbolic info
          call pardiso_mkl_free          ( free_numerical_setup, context, iparm, msglvl )
          call pardiso_mkl_free          ( free_symbolic_setup, context, iparm, msglvl )
          call pardiso_mkl_analysis      ( context, A, iparm, msglvl, perm )
          call pardiso_mkl_factorization ( context, A, iparm, msglvl )
       end select
       ! State transition 
       context%state = num_computed
    case ( pardiso_mkl_solve )
       ! Check pre-conditions
       assert ( context%state ==  num_computed )
       call pardiso_mkl_solution ( context, A, b, x, iparm, msglvl )
    case ( free_numerical_setup )
       call  pardiso_mkl_free ( free_numerical_setup, context, iparm, msglvl )
       context%state=symb_computed
    case ( free_symbolic_setup )
       call  pardiso_mkl_free ( free_symbolic_setup, context, iparm, msglvl )
       context%state=created
    case ( free_clean )
       call  pardiso_mkl_free ( free_clean, context, iparm, msglvl )
       context%state=not_created
    case default
       ! Write an error message and stop ?
    end select
  end subroutine pardiso_mkl_vector
  !=============================================================================
  subroutine pardiso_mkl_r2 ( action, context, A, nrhs, b, ldb, x, ldx, iparm, msglvl, perm )
    implicit none

    ! Mandatory Parameters
    integer(ip)              , intent(in)    :: action
    type(pardiso_mkl_context_t), intent(inout) :: context
    type(serial_scalar_matrix_t), intent(in)             :: A 
    integer(ip)     , intent(in)             :: nrhs, ldb, ldx
    real(rp)        , intent(in)             :: b (ldb, nrhs)
    real(rp)        , intent(inout)          :: x (ldx, nrhs)
    integer         , intent(inout), target, optional :: iparm(64)
    integer                                , optional :: msglvl
    integer         , intent(in)   , target, optional :: perm(A%graph%nv)

    select case(action)

    case ( pardiso_mkl_solve )
       ! Check pre-conditions
       assert ( context%state ==  num_computed )
       call pardiso_mkl_solution_several_rhs ( context, A, nrhs, b, ldb, x, ldx, iparm, msglvl )
    case default
       ! Write an error message and stop ?
    end select

  end subroutine pardiso_mkl_r2
  !=============================================================================
  subroutine pardiso_mkl_r1 ( action, context, A, b, x, iparm, msglvl, perm )
    implicit none

    ! Mandatory Parameters
    integer(ip)              , intent(in)    :: action
    type(pardiso_mkl_context_t), intent(inout) :: context
    type(serial_scalar_matrix_t), intent(in)             :: A 
    real(rp)        , intent(in)             :: b (A%graph%nv)
    real(rp)        , intent(inout)          :: x (A%graph%nv)
    integer         , intent(inout), target, optional :: iparm(64)
    integer                                , optional :: msglvl
    integer         , intent(in)   , target, optional :: perm(A%graph%nv)

    select case(action)

    case ( pardiso_mkl_solve )
       ! Check pre-conditions
       assert ( context%state ==  num_computed )
       call pardiso_mkl_solution_real ( context, A, b, x, iparm, msglvl )
    case default
       ! Write an error message and stop ?
    end select

  end subroutine pardiso_mkl_r1

  !=============================================================================
#ifndef ENABLE_MKL
  subroutine enable_mkl_error_message
      implicit none
      write (0,*) 'Error: Fem was not compiled with -DENABLE_MKL.'
      write (0,*) "Error: You must activate this cpp macro in order to use Intel MKL's interface to PARDISO"
      check(1==0)
  end subroutine
#endif

  !=============================================================================
  subroutine pardiso_mkl_init ( context, matrix, iparm)
    implicit none
    ! Parameters
    type(serial_scalar_matrix_t)                         , intent(in)   :: matrix
    type(pardiso_mkl_context_t)            , intent(out)  :: context
    integer                     , optional , intent(out)  :: iparm(64)
    ! Locals
    integer :: mtype, i

#ifdef ENABLE_MKL
    ! Initiliaze the internal solver memory pointer. This is only
    ! necessary before FIRST call of PARDISO.
    allocate  ( context%pt(64) )
    do i = 1, 64
       context%pt(i)%dummy =  0 
    end do

    ! Choose pardiso_mkl matrix according to matrix
    if(matrix%graph%symmetric_storage.and.matrix%is_symmetric.and.matrix%sign == positive_definite) then
       mtype = pardiso_mkl_spd
    else if(matrix%graph%symmetric_storage.and.matrix%is_symmetric.and.matrix%sign /= positive_definite) then
       mtype = pardiso_mkl_sin
    else ! if(.not. matrix%graph%symmetric_storage) then
       mtype = pardiso_mkl_uss
       !!!!!!!!!! PROVISIONAL for uns mtype !!!!!!!!!!
       ! mtype = pardiso_mkl_uns
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end if

    context%mtype = mtype

    ! Set PARDISO default parameters
    if ( present(iparm) ) then
       do i = 1, 64
          iparm(i) = 0
       end do
       if ( mtype == pardiso_mkl_uns ) then
!!$       !!!!!!!!!! PROVISIONAL for uns mtype !!!!!!!!!!
          iparm(1) = 1 ! no solver default
          iparm(2) = 2 ! fill-in reordering from METIS
          iparm(3) = 1 ! numbers of processors
          iparm(4) = 0 ! no iterative-direct algorithm
          iparm(5) = 0 ! no user fill-in reducing permutation
          iparm(6) = 0 ! =0 solution on the first n compoments of x
          iparm(7) = 0 ! not in use
          iparm(8) = 9 ! numbers of iterative refinement steps
          iparm(9) = 0 ! not in use
          iparm(10) = 13 ! perturbe the pivot elements with 1E-13
          iparm(11) = 0 ! not use nonsymmetric permutation and scaling MPS
          iparm(12) = 0 ! not in use
          iparm(13) = 0 ! maximum weighted matching algorithm is switched-on (default for non-symmetric)
          iparm(14) = 0 ! Output: number of perturbed pivots
          iparm(15) = 0 ! not in use
          iparm(16) = 0 ! not in use
          iparm(17) = 0 ! not in use
          iparm(18) = -1 ! Output: number of nonzeros in the factor LU
          iparm(19) = -1 ! Output: Mflops for LU factorization
          iparm(20) = 0 ! Output: Numbers of CG Iterations
!!$       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       end if
    end if
#else
    call enable_mkl_error_message
#endif
  end subroutine pardiso_mkl_init

  !=============================================================================
  subroutine pardiso_mkl_free ( mode, context, iparm, msglvl )
    implicit none

    ! Parameters
    integer(ip)              , intent(in)                   :: mode
    type(pardiso_mkl_context_t), intent(inout)                :: context
    integer                  , intent(in), target, optional :: iparm(64)
    integer                  , optional                     :: msglvl

    ! Locals
    integer          :: i
    integer          :: phase
    integer          :: error  
    logical          :: alloc_iparm
    integer, pointer :: iparm_(:)
    integer          :: msglvl_
    integer          :: idum(1)
    real(dp)         :: ddum(1)

#ifdef ENABLE_MKL
    ! Free PARDISO internal state allocatable

    if ( mode == free_clean ) then
         deallocate ( context%pt     )
         context%n      = -1
         context%mtype  = -1500
         return  
    end if 

    if ( present(msglvl) ) then
       msglvl_ = msglvl
    else
       msglvl_ = 0
    end if

    alloc_iparm = .false. 
    if ( present(iparm) ) then
       iparm_ => iparm
    else
       call memallocp (64, iparm_, __FILE__,__LINE__ )
       alloc_iparm = .true.
       ! Set PARDISO defaults
       do i = 1, 64
          iparm_(i) = 0
       end do
    end if

    if ( mode == free_symbolic_setup  ) then
       ! Termination and release of memory
       phase = -1 ! release internal memory
       if (context%n>1) then
          call pardiso ( context%pt, maxfct, mnum, context%mtype,phase,    & 
               &         context%n, ddum, idum, idum, idum, nrhsp, iparm_, &  
               &         msglvl_, ddum, ddum, error)
          if (error /= 0) then
             write (0,*) 'Error, PARDISO_MKL: the following ERROR was detected: ', & 
                  error, 'during stage', phase
             stop
          end if
       end if
    else if ( mode == free_numerical_setup ) then
       ! Release internal memory only for L and U factors
       phase = 0
       if (context%n>1) then
          call pardiso ( context%pt, maxfct, mnum, context%mtype, phase, & 
               &         context%n, ddum, idum, idum, idum, nrhsp, iparm_, &  
               &         msglvl_, ddum, ddum, error)
          if (error /= 0) then
             write (0,*) 'Error, PARDISO_MKL: the following ERROR was detected: ', & 
                  error, 'during stage', phase
             stop
          end if
       end if
    end if

    if ( alloc_iparm ) then
       call memfreep ( iparm_,__FILE__,__LINE__)
    end if

#else
    call enable_mkl_error_message
#endif
  
  end subroutine pardiso_mkl_free

  !=============================================================================
  subroutine pardiso_mkl_analysis ( context, matrix, iparm, msglvl, perm )
    implicit none

    ! Parameters 
    type(pardiso_mkl_context_t), intent(inout)                 :: context
    type(serial_scalar_matrix_t)         , intent(in)                    :: matrix
    integer                  , intent(in), target, optional  :: perm(matrix%graph%nv)
    integer                  , intent(in), target, optional  :: iparm(64)
    integer                  , intent(in), optional          :: msglvl

    ! Locals (PARDISO-RELATED)
    integer, pointer  :: perm_ (:)
    integer, pointer  :: iparm_(:)
    integer           :: msglvl_
    integer           :: error
    integer           :: phase

    ! Locals
    integer         :: i
    logical         :: alloc_iparm
    integer, target :: idum(1)
    real(dp)        :: ddum(1)   

#ifdef ENABLE_MKL

    ! Process optional parameters
    if ( present(msglvl) ) then
       msglvl_ = msglvl
    else
       msglvl_ = 0
    end if

    if ( present (perm) ) then
       perm_ => perm
    else
       perm_ => idum
    end if

    alloc_iparm = .false. 
    if ( present(iparm) ) then
       iparm_ => iparm
    else
       call memallocp (64, iparm_, __FILE__,__LINE__ )
       alloc_iparm = .true.
       ! Set PARDISO defaults
       do i = 1, 64
          iparm_(i) = 0
       end do
    end if

    context%n = matrix%graph%nv

    ! Reordering and symbolic factorization, this step also allocates 
    ! all memory that is necessary for the factorization
    phase = 11 ! only reordering and symbolic factorization
    call pardiso ( context%pt, maxfct, mnum, context%mtype, phase, & 
         &         matrix%graph%nv, ddum, matrix%graph%ia, matrix%graph%ja, &
         &         perm_, nrhsp,  iparm_, msglvl_, ddum, ddum, error )
    if (error /= 0) then
       write (0,*) 'Error, PARDISO_MKL: the following ERROR was detected: ', & 
            error, 'during stage', phase
       stop
    end if

    if ( alloc_iparm ) then
       call memfreep ( iparm_,__FILE__,__LINE__)
    end if
#else
    call enable_mkl_error_message
#endif

  end subroutine pardiso_mkl_analysis

  !=============================================================================
  subroutine pardiso_mkl_factorization ( context, matrix, iparm, msglvl )
    implicit none
    ! Parameters
    type(pardiso_mkl_context_t), intent(inout)                :: context
    type(serial_scalar_matrix_t)         , intent(in), target           :: matrix
    integer                  , intent(in), target, optional :: iparm(64)
    integer                  , intent(in), optional         :: msglvl

    ! Locals (PARDISO-RELATED)
    integer, pointer  :: iparm_(:)
    integer           :: msglvl_
    integer           :: error
    integer           :: idum(1)
    real(dp)          :: ddum(1) 
    integer           :: phase

    ! Locals
    integer           :: i
    logical           :: alloc_iparm 
    real(dp), pointer :: a_(:)


#ifdef ENABLE_MKL

    ! Process optional parameters
    if ( present(msglvl) ) then
       msglvl_ = msglvl
    else
       msglvl_ = 0
    end if

    alloc_iparm = .false. 
    if ( present(iparm) ) then
       iparm_ => iparm
    else
       call memallocp (64, iparm_, __FILE__,__LINE__ )
       alloc_iparm = .true.
       ! Set PARDISO defaults
       do i = 1, 64
          iparm_(i) = 0
       end do
    end if

    ! Factorization.
    phase = 22 ! only factorization
    a_ => matrix%a(:)

    call pardiso (  context%pt, maxfct, mnum, context%mtype, phase,  & 
         &          matrix%graph%nv, a_, matrix%graph%ia, matrix%graph%ja,  &
         &          idum, nrhsp, iparm_, msglvl_, ddum, ddum, error)
    if (error /= 0) then
       write (0,*) 'Error, PARDISO_MKL: the following ERROR was detected: ', & 
            error, 'during stage', phase, 'in row', iparm_(30)
       stop
    end if

    if ( alloc_iparm ) then
       call memfreep ( iparm_,__FILE__,__LINE__)
    end if

#else
    call enable_mkl_error_message
#endif
  end subroutine pardiso_mkl_factorization

  !=============================================================================
  subroutine pardiso_mkl_solution ( context, matrix, x, y, iparm, msglvl )
    ! Computes y <- A^-1 * x, using previously computed LU factorization
    implicit none
    ! Parameters 
    type(pardiso_mkl_context_t), intent(inout)                   :: context
    type(serial_scalar_matrix_t)         , intent(in)   , target           :: matrix
    type(serial_scalar_array_t)         , intent(in)   , target           :: x
    type(serial_scalar_array_t)         , intent(inout), target           :: y
    integer                  , intent(in)   , target, optional :: iparm(64)
    integer                  , intent(in)   , optional         :: msglvl

    ! Locals (PARDISO-RELATED)
    integer, pointer  :: iparm_(:)
    integer           :: msglvl_
    integer           :: error
    integer           :: idum(1)
    real(dp)          :: ddum(1) 
    integer           :: phase
    real(dp), pointer :: a_(:)
    real(dp), pointer :: x_(:)
    real(dp), pointer :: y_(:)

    ! Locals
    integer           :: i
    logical           :: alloc_iparm

#ifdef ENABLE_MKL
    ! Process optional parameters
    if ( present(msglvl) ) then
       msglvl_ = msglvl
    else
       msglvl_ = 0
    end if

    alloc_iparm = .false. 
    if ( present(iparm) ) then
       iparm_ => iparm
    else
       call memallocp (64, iparm_, __FILE__,__LINE__ )
       alloc_iparm = .true.
       ! Set PARDISO defaults
       do i = 1, 64
          iparm_(i) = 0
       end do
    end if

    ! (c) y  <- A^-1 * x
    phase = 33 ! only Fwd/Bck substitution
    a_ => matrix%a(:)
    x_ => x%b(:)
    y_ => y%b(:)
    call pardiso ( context%pt, maxfct, mnum, context%mtype, phase,   & 
         &         matrix%graph%nv, a_, matrix%graph%ia, matrix%graph%ja,  & 
         &         idum, nrhsp, iparm_, msglvl_, x_, y_, error)
    if (error /= 0) then
       write (0,*) 'Error, PARDISO_MKL: the following ERROR was detected: ', & 
            error, 'during stage', phase
       stop
    end if

    if ( alloc_iparm ) then
       call memfreep ( iparm_,__FILE__,__LINE__)
    end if

#else
    call enable_mkl_error_message
#endif

  end subroutine pardiso_mkl_solution

  !=============================================================================
  subroutine pardiso_mkl_solution_real ( context, matrix, rhs, sol, iparm, msglvl )
    ! Computes y <- A^-1 * x, using previously computed LU factorization
    implicit none
    ! Parameters 
    type(pardiso_mkl_context_t), intent(inout)                   :: context
    type(serial_scalar_matrix_t)         , intent(in)   , target           :: matrix
    real(rp)                 , intent(in)   , target           :: rhs (matrix%graph%nv)
    real(rp)                 , intent(inout), target           :: sol (matrix%graph%nv)
    integer                  , intent(in)   , target, optional :: iparm(64)
    integer                  , intent(in)   , optional         :: msglvl

    ! Locals (PARDISO-RELATED)
    integer, pointer  :: iparm_(:)
    integer           :: msglvl_
    integer           :: error
    integer           :: idum(1)
    real(dp)          :: ddum(1) 
    integer           :: nrhs
    integer           :: phase
    real(dp), pointer :: a_(:)
    real(dp), pointer :: x_(:)
    real(dp), pointer :: y_(:)

    ! Locals
    integer           :: i
    logical           :: alloc_iparm

#ifdef ENABLE_MKL

    ! Process optional parameters
    if ( present(msglvl) ) then
       msglvl_ = msglvl
    else
       msglvl_ = 0
    end if

    alloc_iparm = .false. 
    if ( present(iparm) ) then
       iparm_ => iparm
    else
       call memallocp (64, iparm_, __FILE__,__LINE__ )
       alloc_iparm = .true.
       ! Set PARDISO defaults
       do i = 1, 64
          iparm_(i) = 0
       end do
    end if

    ! (c) y  <- A^-1 * x
    phase = 33 ! only Fwd/Bck substitution
    a_ => matrix%a(:)
    x_ => rhs
    y_ => sol
    nrhs = 1

    call pardiso ( context%pt, & 
         maxfct, &
         mnum, &
         context%mtype, & 
         phase, & 
         matrix%graph%nv, &
         a_, &
         matrix%graph%ia, & 
         matrix%graph%ja, &
         idum, &
         nrhs, & 
         iparm_, &
         msglvl_, &
         x_, &       
         y_, & 
         error)

    if (error /= 0) then
       write (0,*) 'Error, PARDISO_MKL: the following ERROR was detected: ', & 
            error, 'during stage', phase
       stop
    end if

    if ( alloc_iparm ) then
       call memfreep ( iparm_,__FILE__,__LINE__)
    end if

#else
    call enable_mkl_error_message
#endif

  end subroutine pardiso_mkl_solution_real

  !=============================================================================
  subroutine pardiso_mkl_solution_several_rhs ( context, matrix, nrhs, rhs, ldrhs, sol, ldsol, iparm, msglvl )
    ! Computes y <- A^-1 * x, using previously computed LU factorization
    implicit none
    ! Parameters 
    type(pardiso_mkl_context_t), intent(inout)                   :: context
    type(serial_scalar_matrix_t)         , intent(in)   , target           :: matrix
    integer(ip)              , intent(in)                      :: nrhs, ldrhs, ldsol
    real(rp)                 , intent(in)   , target           :: rhs (ldrhs, nrhs)
    real(rp)                 , intent(inout), target           :: sol (ldsol, nrhs)
    integer                  , intent(in)   , target, optional :: iparm(64)
    integer                  , intent(in)   , optional         :: msglvl

    ! Locals (PARDISO-RELATED)
    integer, pointer      :: iparm_(:)
    integer               :: msglvl_
    integer               :: error
    integer               :: idum(1)
    real(dp)              :: ddum(1) 
    integer               :: phase
    real(dp), pointer     :: a_(:)
    real(dp), pointer     :: x_(:)
    real(dp), pointer     :: y_(:)
    real(dp), allocatable, target :: workrhs(:,:)
    real(dp), allocatable, target :: worksol(:,:)

    ! Locals
    integer           :: i
    logical           :: alloc_iparm

#ifdef ENABLE_MKL

    ! Process optional parameters
    if ( present(msglvl) ) then
       msglvl_ = msglvl
    else
       msglvl_ = 0
    end if

    alloc_iparm = .false. 
    if ( present(iparm) ) then
       iparm_ => iparm
    else
       call memallocp (64, iparm_, __FILE__,__LINE__ )
       alloc_iparm = .true.
       ! Set PARDISO defaults
       do i = 1, 64
          iparm_(i) = 0
       end do
    end if

    ! (c) y  <- A^-1 * x
    phase = 33 ! only Fwd/Bck substitution
    a_ => matrix%a(:)

    if ( ldrhs == matrix%graph%nv ) then
       x_ => rhs(1:1,1)
    else
       call memalloc (matrix%graph%nv, nrhs, workrhs, __FILE__,__LINE__ )
       do i=1,nrhs
          workrhs(:,i) = rhs(1:matrix%graph%nv,i)
       end do
       x_ => workrhs(1:1,1)
    end if

    if ( ldsol == matrix%graph%nv ) then
       y_ => sol(1:1,1)
    else
       call memalloc (matrix%graph%nv, nrhs, worksol, __FILE__,__LINE__ )
       y_ => worksol(1:1,1)
    end if

    call pardiso ( context%pt, & 
         maxfct, &
         mnum, &
         context%mtype, & 
         phase, & 
         matrix%graph%nv, &
         a_, &
         matrix%graph%ia, & 
         matrix%graph%ja, &
         idum, &
         nrhs, & 
         iparm_, &
         msglvl_, &
         x_, &       
         y_, & 
         error)

    if (error /= 0) then
       write (0,*) 'Error, PARDISO_MKL: the following ERROR was detected: ', & 
            error, 'during stage', phase
       stop
    end if

    if ( ldrhs /= matrix%graph%nv ) then
       call memfree (workrhs, __FILE__,__LINE__ )
    end if

    if ( ldsol /= matrix%graph%nv ) then
       do i=1,nrhs
          sol(1:matrix%graph%nv,i) = worksol(:,i)
       end do
       call memfree (worksol, __FILE__,__LINE__ )
    end if

    if ( alloc_iparm ) then
       call memfreep ( iparm_,__FILE__,__LINE__)
    end if

#else
    call enable_mkl_error_message
#endif

  end subroutine pardiso_mkl_solution_several_rhs

end module pardiso_mkl_names
