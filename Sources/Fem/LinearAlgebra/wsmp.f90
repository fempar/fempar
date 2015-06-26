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
module wsmp_names
  ! This module is a wrapper in which we define functions 
  ! to interact with WSMP using our data structures.
  ! Error control on calling parameters must (ideally)
  ! be performed here.

  ! Serial modules
use types_names
use memor_names
  use fem_matrix_names
  use fem_vector_names 
!  use fem_graph_names
!  use fem_graph_partition
# include "debug.i90"
  
  implicit none
  private

  ! Type definition for wsmp args
  integer , parameter :: dp      = kind (1.0D0)

  ! Total number of context initialized
  integer(ip), save      :: num_contexts_symm   = 0
  integer(ip), save      :: num_contexts_unsy   = 0
  integer(ip), save      :: context_status (64) = 0
  integer(ip), parameter :: max_contexts_symm   = 63
  integer(ip), parameter :: max_contexts_unsy   = 1

  ! Posible matrix types managed in wsmp_context_t
  integer , parameter    :: wsmp_symm = 2  ! Symmetric
  integer , parameter    :: wsmp_unsy = 1  ! Unsymmetric

  ! Posible matrix signs managed in wsmp_context
  integer , parameter    :: wsmp_undefined = 0
  integer , parameter    :: wsmp_positive  = 1

  ! Possible states of wsmp_context
  integer(ip), parameter :: not_created   = 0
  integer(ip), parameter :: created       = 1 
  integer(ip), parameter :: symb_computed = 2 ! Symbolic data already computed
  integer(ip), parameter :: num_computed  = 3 ! Numerical data already computed 

  ! This interface does not compile with IBM XLF compilers, so that
  ! I decided to comment it out. Any reason ?
!!$xlf90_r -c -q64 -qrealsize=8 -qsuffix=f=f90:cpp=f90 -WF,-DENABLE_BLAS  -WF,-DENABLE_WSMP -I../../Sources/Generic -IXLF/Objects_O -qmoddir=XLF/Objects_O -O3 -qstrict -qtune=ppc970 -qarch=ppc970 -qcache=auto -o XLF/Objects_O/wsmp.o ../../Sources/Fem/wsmp.f90
!!$"../../Sources/Fem/wsmp.f90", line 377.41: 1513-061 (S) Actual argument attributes do not match those specified by an accessible explicit interface.
!!$"../../Sources/Fem/wsmp.f90", line 408.44: 1513-061 (S) Actual argument attributes do not match those specified by an accessible explicit interface.
!!$"../../Sources/Fem/wsmp.f90", line 574.63: 1513-061 (S) Actual argument attributes do not match those specified by an accessible explicit interface.
!!$"../../Sources/Fem/wsmp.f90", line 690.63: 1513-061 (S) Actual argument attributes do not match those specified by an accessible explicit interface.
!!$"../../Sources/Fem/wsmp.f90", line 812.63: 1513-061 (S) Actual argument attributes do not match those specified by an accessible explicit interface.
!!$"../../Sources/Fem/wsmp.f90", line 936.63: 1513-061 (S) Actual argument attributes do not match those specified by an accessible explicit interface.
!!$"../../Sources/Fem/wsmp.f90", line 1064.63: 1513-061 (S) Actual argument attributes do not match those specified by an accessible explicit interface.
!!$** wsmp_names   === End of Compilation 1 ===
!!$1501-511  Compilation failed for file wsmp.f90.
!!$make[2]: *** [XLF/Objects_O/wsmp.o] Error 1
!!$make[1]: *** [release_serial] Error 2
!!$make: *** [release] Error 2



!!$  interface 
!!$     subroutine wssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs, &
!!$          &                    aux, naux, mrp, iparm, dparm)
!!$       implicit none
!!$       integer n, ia(*), ja(*)
!!$       double precision avals(*), diag(*)
!!$       integer ldb, nrhs
!!$       double precision b(ldb,*)
!!$       integer perm(*), invp(*)
!!$       double precision aux
!!$       integer naux
!!$       integer mrp(*)
!!$       integer iparm(*)
!!$       double precision dparm(*)
!!$     end subroutine wssmp
!!$  end interface

  type wsmp_context_t
     ! Our components
     integer(ip) :: mtype = wsmp_symm      ! Matrix symmetry (check consistency between calls)
     integer(ip) :: sign  = wsmp_undefined ! Matrix symmetry (check consistency between calls)
     ! WSMP components
     integer     :: id    = 0              ! wsmp internal context id
     integer     :: n     = 0              ! Number of rows/cols of the matrix    
     integer     :: state = not_created
     ! Permutation and parameters arrays are used in the symmetric
     ! case only and must be consistent between calls (the output of
     ! ordering phase is to be send to symbolic phase where they are
     ! modified again and need to be send to numerical phase).
     integer, pointer :: perm(:)  => NULL()
     integer, pointer :: iperm(:) => NULL()
     integer, pointer :: iparm(:) => NULL()
     real   , pointer :: dparm(:) => NULL()
  end type wsmp_context_t

  ! Types
  public :: wsmp_context_t

  ! Constants (matrix types)
  public :: wsmp_symm, &
       &    wsmp_unsy


  ! Possible actions that can be perfomed by solve_wsmp
  integer(ip), parameter :: wsmp_init              = 1  ! Construct solve_wsmp_state
  integer(ip), parameter :: wsmp_finalize          = 2  ! Destruct  solve_wsmp_state     
  integer(ip), parameter :: wsmp_compute_symb      = 3  ! Compute symb. fact.  
  integer(ip), parameter :: wsmp_compute_num       = 4  ! Compute numerical fact. 
  integer(ip), parameter :: wsmp_solve             = 6  ! Fwd./Bck. substitution 
  integer(ip), parameter :: wsmp_free_values       = 7
  integer(ip), parameter :: wsmp_free_struct       = 8
  integer(ip), parameter :: wsmp_free_clean        = 9

  interface wsmp
     module procedure wsmp_vector, wsmp_r2, wsmp_r1
  end interface wsmp

  ! Constants (actions)
  public :: wsmp_init,  &
       &    wsmp_finalize, &
       &    wsmp_compute_symb, & 
       &    wsmp_compute_num, &
       &    wsmp_solve, &
       &    wsmp_free_values,&
       &    wsmp_free_struct, &
       &    wsmp_free_clean

  ! Functions
  public :: wsmp ! , wsmp_solution_several_rhs

contains
  ! WSMP calling routines and their arguments
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
  subroutine wsmp_vector ( action, context, A, b, x, iparm, dparm, perm, iperm, mrp, rmisc)
    implicit none
    ! Mandatory Parameters
    type(wsmp_context_t), intent(inout) :: context   ! Information required between calls
    integer(ip)       , intent(in)    :: action    ! Action to be performed 
                                                   ! (see public constants above)
    type(fem_matrix_t)  , intent(in)    :: A         ! Linear system coefficient matrix
    type(fem_vector_t)  , intent(in)    :: b         ! RHS (Right-hand-side)
    type(fem_vector_t)  , intent(inout) :: x         ! LHS (Left-hand-side)
    ! Optional parameters
!!$    integer , intent(inout), target, optional :: perm(A%gr%nv)
!!$    integer , intent(inout), target, optional :: iperm(A%gr%nv)
!!$    integer , intent(inout), target, optional :: mrp(A%gr%nv)
!!$    real(8) , intent(inout), target, optional :: rmisc(A%gr%nv)
    integer , intent(inout), target, optional :: perm(*)
    integer , intent(inout), target, optional :: iperm(*)
    integer , intent(inout), target, optional :: mrp(*)
    real(8) , intent(inout), target, optional :: rmisc(*)
    integer , intent(inout), target, optional :: iparm(64)
    real(8) , intent(inout), target, optional :: dparm(64)

    select case(action)

    case ( wsmp_init )

       ! Check pre-conditions
       assert ( context%state == not_created )
       call wsmp_ini ( context, A, iparm, dparm)
       ! State transition 
       context%state = created

    case ( wsmp_finalize )

       ! Check pre-conditions
       assert ( context%state /= not_created )
       select case (context%state)
       case (created)
          call  wsmp_free ( wsmp_free_clean , context )
       case (symb_computed)
          call  wsmp_free ( wsmp_free_struct, context )
          call  wsmp_free ( wsmp_free_clean , context )
       case (num_computed)
          call  wsmp_free ( wsmp_free_values, context )
          call  wsmp_free ( wsmp_free_struct, context )
          call  wsmp_free ( wsmp_free_clean , context )
       end select
       ! State transition 
       context%state = not_created

    case ( wsmp_compute_symb )

       ! Check pre-conditions
       assert (context%state==created)
       ! Create sparsity pattern of LU sparse direct factorization
       call wsmp_analysis ( context, A, iparm, dparm, perm, iperm )
       ! State transition 
       context%state = symb_computed

    case ( wsmp_compute_num )

       ! Check pre-conditions
       assert (context%state==symb_computed.or.context%state==num_computed)
       ! Compute LU sparse direct factorization
          call wsmp_factorization ( context, A, iparm, dparm, perm, iperm, mrp )             
       ! State transition 
       context%state = num_computed

    case ( wsmp_solve )

       ! Check pre-conditions
       assert ( context%state ==  num_computed )
       call wsmp_solution ( context, A, b, x, iparm, dparm, perm, iperm, rmisc )

    case ( wsmp_free_values )

       call  wsmp_free ( wsmp_free_values, context )
       context%state=symb_computed

    case ( wsmp_free_struct )

       call  wsmp_free ( wsmp_free_struct, context )
       context%state=created

    case ( wsmp_free_clean )

       call  wsmp_free ( wsmp_free_clean, context )
       context%state=not_created

    case default

       ! Write an error message and stop ?      

    end select

  end subroutine wsmp_vector
  !=============================================================================
  subroutine wsmp_r2 ( action, context, A, nrhs, b, x, iparm, dparm, perm, iperm, mrp, rmisc)
    implicit none
    ! Mandatory Parameters
    type(wsmp_context_t), intent(inout) :: context   ! Information required between calls
    integer(ip)       , intent(in)    :: action    ! Action to be performed 
                                                   ! (see public constants above)
    type(fem_matrix_t)  , intent(in)    :: A         ! Linear system coefficient matrix
!    type(fem_vector_t)  , intent(in)    :: b         ! RHS (Right-hand-side)
!    type(fem_vector_t)  , intent(inout) :: x         ! LHS (Left-hand-side)
    integer(ip)       , intent(in)            :: nrhs
    real(rp)          , intent(in)    :: b (A%gr%nv, nrhs)
    real(rp)          , intent(inout) :: x (A%gr%nv, nrhs)

    ! Optional parameters
    integer , intent(inout), target, optional :: perm(A%gr%nv)
    integer , intent(inout), target, optional :: iperm(A%gr%nv)
    integer , intent(inout), target, optional :: mrp(A%gr%nv)
    real(8) , intent(inout), target, optional :: rmisc(A%gr%nv)
    integer , intent(inout), target, optional :: iparm(64)
    real(8) , intent(inout), target, optional :: dparm(64)

    select case(action)

    case ( wsmp_solve )

       ! Check pre-conditions
       assert ( context%state ==  num_computed )
       call wsmp_solution_several_rhs ( context, A, nrhs, b, x, iparm, dparm, perm, iperm, rmisc )

    case default

       ! Write an error message and stop ?      

    end select

  end subroutine wsmp_r2

  !=============================================================================
  subroutine wsmp_r1 ( action, context, A, b, x, iparm, dparm, perm, iperm, mrp, rmisc)
    implicit none
    ! Mandatory Parameters
    type(wsmp_context_t), intent(inout) :: context   ! Information required between calls
    integer(ip)       , intent(in)    :: action    ! Action to be performed 
                                                   ! (see public constants above)
    type(fem_matrix_t)  , intent(in)    :: A         ! Linear system coefficient matrix
    real(rp)          , intent(in)    :: b (A%gr%nv)
    real(rp)          , intent(inout) :: x (A%gr%nv)

    ! Optional parameters
    integer , intent(inout), target, optional :: perm(A%gr%nv)
    integer , intent(inout), target, optional :: iperm(A%gr%nv)
    integer , intent(inout), target, optional :: mrp(A%gr%nv)
    real(8) , intent(inout), target, optional :: rmisc(A%gr%nv)
    integer , intent(inout), target, optional :: iparm(64)
    real(8) , intent(inout), target, optional :: dparm(64)

    select case(action)

    case ( wsmp_solve )

       ! Check pre-conditions
       assert ( context%state ==  num_computed )
       call wsmp_solution_real ( context, A, b, x, iparm, dparm, perm, iperm, rmisc )

    case default

       ! Write an error message and stop ?      

    end select
    
  end subroutine wsmp_r1

  !=============================================================================
#ifndef ENABLE_WSMP
  subroutine enable_wsmp_error_message
      implicit none
      write (0,*) 'Error: Fem was not compiled with -DENABLE_WSMP.'
      write (0,*) "Error: You must activate this cpp macro in order to use WSMP"
      stop
  end subroutine
#endif

  !=============================================================================
  subroutine wsmp_ini ( context, matrix, iparm, dparm)
    implicit none
    ! Parameters
    type(wsmp_context_t), intent(inout), target :: context
    type(fem_matrix_t)  , intent(in)            :: matrix
    integer, intent(out), target, optional    :: iparm(64)
    real   , intent(out), target, optional    :: dparm(64)
    ! Locals
    integer , pointer :: iparm_(:)
    real(8) , pointer :: dparm_(:)
    integer  :: i1dum(1),i0dum, nrhs, naux, idc, info
    real(8)  :: d1dum(1),d0dum, aux 
    naux = 0
    nrhs = 1

#ifdef ENABLE_WSMP

    ! Choose wsmp context according to matrix
    if(matrix%symm == symm_true) then
       if ( num_contexts_symm == 0 ) then 
          call wsmp_initialize
!          write(*,*) 'wsmp_initialize'
       end if
       do idc=0, max_contexts_symm-1
          if ( context_status(idc+1) == 0 ) then
             context_status(idc+1) = 1
             context%id = idc
             context%mtype = wsmp_symm
             num_contexts_symm = num_contexts_symm +1
             exit
          end if
       end do
!!$       if(num_contexts_symm < max_contexts_symm) then
!!$          num_contexts_symm = num_contexts_symm +1
!!$          context%id = num_contexts_symm
!!$          context%mtype = wsmp_symm
!!$       else
!!$          write (0,*) 'Error, WSMP: reached maximum number of symmetric contexts'
!!$          stop
!!$       end if
    else if(matrix%symm == symm_false) then
       if(num_contexts_unsy < max_contexts_unsy) then
          num_contexts_unsy = num_contexts_unsy +1
          context%id = num_contexts_unsy
          context%mtype = wsmp_unsy
       else
          write (0,*) 'Error, WSMP: reached maximum number of unsymmetric contexts'
          stop
       end if
    end if
    if(matrix%sign == positive_definite) then
       context%sign = wsmp_positive
    else if(matrix%sign /= positive_definite) then
       context%sign = wsmp_undefined
    end if

    ! Allocate and fill default parameters just in case
    ! they are not given by the user (in subsequent calls)
    call memallocp (64, context%iparm, __FILE__,__LINE__ )
    call memallocp (64, context%dparm, __FILE__,__LINE__ )
    context%iparm(1) = 0
    context%iparm(2) = 0
    context%iparm(3) = 0
    if(context%mtype==wsmp_symm) then
       call wssmp (i0dum, i1dum, i1dum, d1dum, d1dum, i1dum, i1dum, d1dum, &
            &      i0dum, nrhs,aux, naux, i1dum, context%iparm, context%dparm)
       if(context%sign == wsmp_undefined) then
          context%iparm(31)=2 ! LDLt factorization with pivoting (default is 0, cholesky)
       end if
    else if(context%mtype==wsmp_unsy) then
       call wgsmp (i0dum, i1dum, i1dum, d1dum, d1dum, i0dum, nrhs, d1dum, &
            &      context%iparm, context%dparm)      
    end if
    if (context%iparm(64) /= 0) then
       write (0,*) 'Error, WSMP: the following ERROR was detected: ', & 
            &       context%iparm(64), 'during initialization'
       stop
    end if

    ! Now fill user parameters, if present
    if ( present(iparm).or.present(dparm) ) then
       if ( present(iparm) ) then
          iparm_ => iparm
       else
          iparm_ => context%iparm
       end if
       if ( present(dparm) ) then
          dparm_ => dparm
       else
          dparm_ => context%dparm
       end if
       iparm_(1) = 0
       iparm_(2) = 0
       iparm_(3) = 0
       if(context%mtype==wsmp_symm) then
          call wssmp (i0dum, i1dum, i1dum, d1dum, d1dum, i1dum, i1dum, d1dum, &
               &      i0dum, nrhs,aux, naux, i1dum, iparm_, dparm_)
          if(context%sign == wsmp_undefined) then
             iparm_(31)=2 ! LDLt factorization with pivoting (default is 0, cholesky)
          end if
       else if(context%mtype==wsmp_unsy) then
          call wgsmp (i0dum, i1dum, i1dum, d1dum, d1dum, i0dum, nrhs, d1dum, &
               &      iparm_, dparm_)      
       end if
       if (iparm_(64) /= 0) then
          write (0,*) 'Error, WSMP: the following ERROR was detected: ', & 
               &       iparm_(64), 'during initialization'
          stop
       end if

    end if

#else
    call enable_wsmp_error_message
#endif
  end subroutine wsmp_ini

  !=============================================================================
  subroutine wsmp_free ( mode, context )
    implicit none

    ! Parameters
    integer(ip)       , intent(in)    :: mode
    type(wsmp_context_t), intent(inout) :: context
    integer                           :: i0dum, nrhs, naux, i, info


#ifdef ENABLE_WSMP
    ! Free wsmp_context structures
    if ( mode == wsmp_free_clean ) then
       context_status(context%id+1) = 0
       num_contexts_symm = num_contexts_symm - 1
       ! Recall context
       call wrecallmat (context%id, info)
       ! write (*,*) 'wsmp_free'
       if(num_contexts_symm == 0 .and. num_contexts_unsy == 0) then
          call wsmp_clear
          ! write (*,*) 'wsmp_clear'
       end if
       call memfreep ( context%iparm,__FILE__,__LINE__)
       call memfreep ( context%dparm,__FILE__,__LINE__)
       if(associated(context%perm)) then
          call memfreep ( context%perm ,__FILE__,__LINE__)
          call memfreep ( context%iperm,__FILE__,__LINE__)
       end if
       return
    end if

    if ( mode == wsmp_free_struct  ) then
       if(context%mtype==wsmp_symm) then
          ! Recall context
          call wrecallmat (context%id, info)
          if (info/=0) then
             write (0,*) 'Error, WSMP: error recalling context',context%id 
             stop
          end if
          !if(num_contexts_symm==0) call wssfree()
          call wssfree()
          call wstoremat (context%id, info)
          if(info/=0) then
             write (0,*) 'WSMP: error storing context',context%id , info
             stop
          end if
       else if(context%mtype==wsmp_unsy) then
          !num_contexts_unsy = num_contexts_unsy - 1
          !if(num_contexts_unsy==0) call wgsfree()
          call wgsfree()
       end if
    else if ( mode == wsmp_free_values ) then
       ! Release internal memory only for L and U factors of current context
       if(context%mtype==wsmp_symm) then
          ! Recall context
          call wrecallmat (context%id, info)
          if (info/=0) then
             write (0,*) 'Error, WSMP: error recalling context',context%id 
             stop
          end if
          call wsffree()
          call wstoremat (context%id , info)
          if(info/=0) then
             write (0,*) 'WSMP: error storing context',context%id , info
             stop
          end if
       else if(context%mtype==wsmp_unsy) then
          call wgffree()
       end if
    end if

#else
    call enable_wsmp_error_message
#endif

  end subroutine wsmp_free

  !=============================================================================
  subroutine wsmp_analysis ( context, matrix, iparm, dparm , perm, iperm)
    implicit none
    ! Parameters 
    type(wsmp_context_t), intent(inout), target :: context
    type(fem_matrix_t)  , intent(in), target    :: matrix
    integer, intent(inout), target, optional  :: perm(matrix%gr%nv)
    integer, intent(inout), target, optional  :: iperm(matrix%gr%nv)
    integer, intent(inout), target, optional  :: iparm(64)
    real   , intent(inout), target, optional  :: dparm(64)
    ! Locals
    integer , pointer :: iparm_(:)
    real(8) , pointer :: dparm_(:)
    integer  :: i1dum(1),i0dum, nrhs, naux, i, info
    real(8)  :: d1dum(1),d0dum, aux 
    real(8) , target  :: adumm(1)
    real(8) , pointer :: a_(:)
    naux = 0
    nrhs = 1

#ifdef ENABLE_WSMP

    ! Check a correct matrix type
    !assert (matrix%type == csr_mat)

    if ( present(perm) ) then
       context%perm => perm
    else
       call memallocp (matrix%gr%nv, context%perm, __FILE__,__LINE__ )
    end if

    if ( present(iperm) ) then
       context%iperm => iperm
    else
       call memallocp (matrix%gr%nv, context%iperm, __FILE__,__LINE__ )
    end if

    if ( present(iparm) ) then
       iparm_ => iparm
    else
       iparm_ => context%iparm
    end if
    if ( present(dparm) ) then
       dparm_ => dparm
    else
       dparm_ => context%dparm
    end if
    
    ! Allocated does not work here!!!
    !if( allocated( matrix%a ) ) then
    !   a_ => matrix%a(1,1,:)
    !else
    !   a_ => adumm
    !end if

    iparm_(2) = 1
    iparm_(3) = 2

    if(context%mtype==wsmp_symm) then

       a_ => adumm
       iparm_(10) = 2

       !write(*,*) 'Performing analysis', iparm_(31)

       call wssmp (matrix%gr%nv, matrix%gr%ia, matrix%gr%ja , a_, d1dum, &
            &      context%perm, context%iperm, d1dum, i0dum, nrhs, aux, &
            &      naux, i1dum, iparm_, dparm_)



       ! call wssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
       !             aux, naux, mrp, context%iparm, context%dparm)
       ! n, ia, ja, avals in matrix
       ! ldb=1 nrhs=1 b in vector 
       ! perm iperm
       ! aux obsolete not used
       ! naux = 0 not used
       ! mrp(n) pivot information

       call wstoremat (context%id , info)
       if(info/=0) then
          write (0,*) 'WSMP analysis: error storing context',context%id , info
          stop
       end if

    else if(context%mtype==wsmp_unsy) then

       a_ => matrix%a(:)

       call wgsmp (matrix%gr%nv, matrix%gr%ia, matrix%gr%ja , a_, d1dum, &
            &      i0dum, nrhs, d1dum, iparm_, dparm_)
       ! call wgsmp (n, ia, ja, avals, b, ldb, nrhs, rmisc, iparm, dparm)
       ! n, ia, ja, avals in matrix
       ! ldb=1 nrhs=1 b in vector 
       ! rmisc(n,nrhs) is the backward error if iparm(25)=1

    end if

    if (iparm_(64) /= 0) then
       write (0,*) 'Error, WSMP: the following ERROR was detected: ', & 
            iparm_(64), 'during analysis'
       stop
    end if

#else
    call enable_wsmp_error_message
#endif

  end subroutine wsmp_analysis

  !=============================================================================
  subroutine wsmp_factorization ( context, matrix, iparm, dparm, perm, iperm, mrp )
    implicit none
    ! Parameters 
    type(wsmp_context_t), intent(inout),target :: context
    type(fem_matrix_t)  , intent(in), target   :: matrix
    integer, intent(inout), target, optional :: perm(matrix%gr%nv)
    integer, intent(inout), target, optional :: iperm(matrix%gr%nv)
    integer, intent(inout), target, optional :: mrp(matrix%gr%nv)
    integer, intent(inout), target, optional :: iparm(64)
    real   , intent(inout), target, optional :: dparm(64)
    ! Locals
    integer  :: i0dum, nrhs, naux, i, info
    real(8)  :: d0dum, aux 
    integer , pointer :: mrp_(:)
    integer , pointer :: iparm_(:)
    real(8) , pointer :: dparm_(:)
    integer , target  :: i1dum(1)
    real(8) , target  :: d1dum(1)
    real(8) , pointer :: a_(:)
    naux = 0
    nrhs = 1

#ifdef ENABLE_WSMP

    ! Check a correct matrix type
    assert (matrix%type == csr_mat)

    if ( present(perm).and. (.not.associated(context%perm,perm)) ) then
       write (0,*) 'Error, WSMP: array perm cannot be changed'
       stop
    end if

    if ( present(iperm).and. (.not.associated(context%iperm,iperm)) ) then
       write (0,*) 'Error, WSMP: array iperm cannot be changed'
       stop
    end if

    if ( present(iparm) ) then
       iparm_ => iparm
    else
       iparm_ => context%iparm
    end if

    if ( present(dparm) ) then
       dparm_ => dparm
    else
       dparm_ => context%dparm
    end if

    ! Point to matrix
    a_ => matrix%a(:)
    iparm_(2) = 3
    iparm_(3) = 3

    if(context%mtype==wsmp_symm) then

       ! Recall context
       call wrecallmat (context%id, info)
       if (info/=0) then
          write (0,*) 'Error, WSMP: error recalling context',context%id 
          stop
       end if

       if ( present(mrp) ) then
          mrp_ => mrp
       else
          mrp_ => i1dum
       end if

       !write(*,*) 'Performing factorization', iparm_(31)
       !write (*,*) 'iparm_', iparm_

       call wssmp (matrix%gr%nv, matrix%gr%ia, matrix%gr%ja , a_, d1dum, &
            &      context%perm, context%iperm, d1dum, i0dum, nrhs, aux, &
            &      naux, mrp, iparm_, dparm_)


       ! call wssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
       !             aux, naux, mrp, context%iparm, context%dparm)
       ! n, ia, ja, avals in matrix
       ! ldb=1 nrhs=1 b in vector 
       ! perm iperm
       ! aux obsolete not used
       ! naux = 0 not used
       ! mrp(n) pivot information

       call wstoremat (context%id , info)
       if(info/=0) then
          write (0,*) 'WSMP factorization: error storing context',context%id , info
          stop
       end if

    else if(context%mtype==wsmp_unsy) then

       call wgsmp (matrix%gr%nv, matrix%gr%ia, matrix%gr%ja , a_, d1dum, &
            &      i0dum, nrhs, d1dum, iparm_, dparm_)
       ! call wgsmp (n, ia, ja, avals, b, ldb, nrhs, rmisc, iparm, dparm)
       ! n, ia, ja, avals in matrix
       ! ldb=1 nrhs=1 b in vector 
       ! rmisc(n,nrhs) is the backward error if iparm(25)=1

    end if

    if (iparm_(64) /= 0) then
       write (0,*) 'Error, WSMP: the following ERROR was detected: ', & 
            &       iparm_(64), 'during factorization'
       stop
    end if

#else
    call enable_wsmp_error_message
#endif

  end subroutine wsmp_factorization

  !=============================================================================
  subroutine wsmp_solution ( context, matrix, x, y, iparm, dparm, perm, iperm, rmisc )
    ! Computes y <- A^-1 * x, using previously computed LU factorization
    implicit none
    ! Parameters 
    type(wsmp_context_t), intent(inout), target :: context
    type(fem_matrix_t)  , intent(in)   , target :: matrix
!    type(fem_vector_t)  , intent(in)   , target :: x
    type(fem_vector_t)  , intent(in)            :: x
    type(fem_vector_t)  , intent(inout), target :: y
    integer, intent(inout), target, optional :: perm(matrix%gr%nv)
    integer, intent(inout), target, optional :: iperm(matrix%gr%nv)
    integer, intent(inout), target, optional :: iparm(64)
    real   , intent(inout), target, optional :: dparm(64)
    real   , intent(inout), target, optional :: rmisc(matrix%gr%nv)
    ! Locals
    integer  :: i0dum, nrhs, naux, i, info
    real(8)  :: d0dum, aux 
    integer , pointer :: iparm_(:)
    real(8) , pointer :: dparm_(:)
    real(8) , pointer :: rmisc_(:)
    integer , target  :: i1dum(1)
    real(8) , target  :: d1dum(1)
    real(8) , pointer :: a_(:)
    real(8) , pointer :: y_(:)
    real(rp) :: xnorm,ynorm
    naux = 0
    nrhs = 1
    aux = 0.0_rp
    d1dum(1) = 0.0_rp

#ifdef ENABLE_WSMP

     write(*,*) 'wsmp solution 0'

    ! Check a correct matrix type
    assert (matrix%type == csr_mat)

    if ( present(perm).and. (.not.associated(context%perm,perm)) ) then
       write (0,*) 'Error, WSMP: array perm cannot be changed'
       stop
    end if

    if ( present(iperm).and. (.not.associated(context%iperm,iperm)) ) then
       write (0,*) 'Error, WSMP: array iperm cannot be changed'
       stop
    end if

    if ( present(iparm) ) then
       iparm_ => iparm
    else
       iparm_ => context%iparm
    end if

    if ( present(dparm) ) then
       dparm_ => dparm
    else
       dparm_ => context%dparm
    end if

     write(*,*) 'wsmp solution 1'

    ! Point to matrix
    a_ => matrix%a(:)
     write(*,*) 'wsmp solution 2'
     write(*,*) 'y',size(y%b)
     write(*,*) 'x',size(x%b)
    ! Solution will overwrite rhs in wsmp
    y%b(:) = x%b(:)
    write(*,*) 'wsmp solution 3'
    y_ => y%b(:)
    write(*,*) 'wsmp solution 4'

    iparm_(2) = 4
    iparm_(3) = 4 ! Set this value to 5 to perform iterative refinement...

    if(context%mtype==wsmp_symm) then

       ! Recall context
       call wrecallmat (context%id, info)
       if (info/=0) then
          write (0,*) 'Error, WSMP: error recalling context',context%id 
          stop
       end if

       !call fem_vector_nrm2(x,xnorm)
       !call fem_vector_nrm2(y,ynorm)

       write(*,*) 'Performing solution', iparm_(10), iparm_(31)

       call wssmp (matrix%gr%nv, matrix%gr%ia, matrix%gr%ja , a_, d1dum, &
            &      context%perm, context%iperm, y_, matrix%gr%nv, nrhs,  &
            &      aux, naux, i1dum, iparm_, dparm_)
       ! call wssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
       !             aux, naux, mrp, context%iparm, context%dparm)
       ! n, ia, ja, avals in matrix
       ! ldb=1 nrhs=1 b in vector 
       ! perm iperm
       ! aux obsolete not used
       ! naux = 0 not used
       ! mrp(n) pivot information

       call wstoremat (context%id , info)
       if(info/=0) then
          write (0,*) 'WSMP solution: error storing context',context%id , info
          stop
       end if

    else if(context%mtype==wsmp_unsy) then

       if ( present(rmisc) ) then
          rmisc_ => rmisc
       else
          rmisc_ => d1dum
       end if

       call wgsmp (matrix%gr%nv, matrix%gr%ia, matrix%gr%ja , a_, y_, &
            &      matrix%gr%nv, nrhs, rmisc_, iparm_, dparm_)
       ! call wgsmp (n, ia, ja, avals, b, ldb, nrhs, rmisc, iparm, dparm)
       ! n, ia, ja, avals in matrix
       ! ldb=1 nrhs=1 b in vector 
       ! rmisc(n,nrhs) is the backward error if iparm(25)=1

    end if

    if (iparm_(64) /= 0) then
       write (0,*) 'Error, WSMP: the following ERROR was detected: ', & 
            &       iparm_(64), 'during solution'
       stop
    end if

#else
    call enable_wsmp_error_message
#endif

  end subroutine wsmp_solution

  !=============================================================================
  subroutine wsmp_solution_real ( context, matrix, rhs, sol, &
       &                          iparm, dparm, perm, iperm ,rmisc )
    implicit none
    ! Parameters 
    type(wsmp_context_t), intent(inout), target :: context
    type(fem_matrix_t)  , intent(in)   , target :: matrix
    real(rp)          , intent(in)   , target :: rhs (matrix%gr%nv)
    real(rp)          , intent(inout), target :: sol (matrix%gr%nv)
    integer, intent(inout), target, optional :: perm(matrix%gr%nv)
    integer, intent(inout), target, optional :: iperm(matrix%gr%nv)
    integer, intent(inout), target, optional :: iparm(64)
    real   , intent(inout), target, optional :: dparm(64)
    real   , intent(inout), target, optional :: rmisc(matrix%gr%nv)
    ! Locals
    integer  :: i0dum, naux, i, info
    real(8)  :: d0dum, aux 
    integer  :: nrhs
    integer , pointer :: iparm_(:)
    real(8) , pointer :: dparm_(:)
    integer , target  :: i1dum(1)
    real(8) , target  :: d1dum(1)
    real(8) , pointer :: a_(:)
    real(8) , pointer :: y_(:)
    real(8) , pointer :: rmisc_(:)
    naux = 0
    nrhs = 1

#ifdef ENABLE_WSMP

    ! Check a correct matrix type
    assert (matrix%type == csr_mat)

    if ( present(perm).and. (.not.associated(context%perm,perm)) ) then
       write (0,*) 'Error, WSMP: array perm cannot be changed'
       stop
    end if

    if ( present(iperm).and. (.not.associated(context%iperm,iperm)) ) then
       write (0,*) 'Error, WSMP: array iperm cannot be changed'
       stop
    end if

    if ( present(iparm) ) then
       iparm_ => iparm
    else
       iparm_ => context%iparm
    end if

    if ( present(dparm) ) then
       dparm_ => dparm
    else
       dparm_ => context%dparm
    end if
    
    ! Point to matrix
    a_ => matrix%a(:)
    ! Solution will overwrite rhs in wsmp
    sol = rhs
    y_ => sol(1:1)

    iparm_(2) = 4
    iparm_(3) = 4 ! Set this value to 5 to perform iterative refinement...

    if(context%mtype==wsmp_symm) then

       ! Recall context
       call wrecallmat (context%id, info)
       if (info/=0) then
          write (0,*) 'Error, WSMP: error recalling context',context%id 
          stop
       end if

       ! Perform both ordering and symbolic factorizations
       ! They are different tasks in the symmetric version
       ! but the same task in the unsymmetric one.
       call wssmp (matrix%gr%nv, matrix%gr%ia, matrix%gr%ja , a_, d1dum, &
            &      context%perm, context%iperm, y_, matrix%gr%nv, nrhs,  &
            &      aux, naux, i1dum, iparm_, dparm_)
       ! call wssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
       !             aux, naux, mrp, context%iparm, context%dparm)
       ! n, ia, ja, avals in matrix
       ! ldb=1 nrhs=1 b in vector 
       ! perm iperm
       ! aux obsolete not used
       ! naux = 0 not used
       ! mrp(n) pivot information

       call wstoremat (context%id , info)
       if(info/=0) then
          write (0,*) 'WSMP solution: error storing context',context%id , info
          stop
       end if

    else if(context%mtype==wsmp_unsy) then

       if ( present(rmisc) ) then
          rmisc_ => rmisc
       else
          rmisc_ => d1dum
       end if

       call wgsmp (matrix%gr%nv, matrix%gr%ia, matrix%gr%ja , a_, y_, &
            &      matrix%gr%nv, nrhs, rmisc_ , iparm_, dparm_)
       ! call wgsmp (n, ia, ja, avals, b, ldb, nrhs, rmisc, iparm, dparm)
       ! n, ia, ja, avals in matrix
       ! ldb=1 nrhs=1 b in vector 
       ! rmisc(n,nrhs) is the backward error if iparm(25)=1

    end if

    if (iparm_(64) /= 0) then
       write (0,*) 'Error, WSMP: the following ERROR was detected: ', & 
            &       iparm_(64), 'during solution with several rhs'
       stop
    end if

#else
    call enable_wsmp_error_message
#endif

  end subroutine wsmp_solution_real

  !=============================================================================
  subroutine wsmp_solution_several_rhs ( context, matrix, nrhs, rhs, sol, &
       &                                 iparm, dparm, perm, iperm ,rmisc )
    implicit none
    ! Parameters 
    type(wsmp_context_t), intent(inout), target :: context
    type(fem_matrix_t)  , intent(in)   , target :: matrix
    integer(ip)       , intent(in)            :: nrhs
    real(rp)          , intent(in)   , target :: rhs (matrix%gr%nv, nrhs)
    real(rp)          , intent(inout), target :: sol (matrix%gr%nv, nrhs)
    integer, intent(inout), target, optional :: perm(matrix%gr%nv)
    integer, intent(inout), target, optional :: iperm(matrix%gr%nv)
    integer, intent(inout), target, optional :: iparm(64)
    real   , intent(inout), target, optional :: dparm(64)
    real   , intent(inout), target, optional :: rmisc(matrix%gr%nv)
    ! Locals
    integer  :: i0dum, naux, i, info
    real(8)  :: d0dum, aux 
    integer , pointer :: iparm_(:)
    real(8) , pointer :: dparm_(:)
    integer , target  :: i1dum(1)
    real(8) , target  :: d1dum(1)
    real(8) , pointer :: a_(:)
    real(8) , pointer :: y_(:)
    real(8) , pointer :: rmisc_(:)

    integer :: irhs
!!$    do irhs = 1, nrhs
!!$      call wsmp_solution_real ( context, matrix, rhs(1,irhs), sol(1,irhs) )
!!$    end do

    naux = 0
#ifdef ENABLE_WSMP

    ! Check a correct matrix type
    assert (matrix%type == csr_mat)

    if ( present(perm).and. (.not.associated(context%perm,perm)) ) then
       write (0,*) 'Error, WSMP: array perm cannot be changed'
       stop
    end if

    if ( present(iperm).and. (.not.associated(context%iperm,iperm)) ) then
       write (0,*) 'Error, WSMP: array iperm cannot be changed'
       stop
    end if

    if ( present(iparm) ) then
       iparm_ => iparm
    else
       iparm_ => context%iparm
    end if

    if ( present(dparm) ) then
       dparm_ => dparm
    else
       dparm_ => context%dparm
    end if
    
    ! Point to matrix
    a_ => matrix%a(:)
    ! Solution will overwrite rhs in wsmp
    sol = rhs
    y_ => sol(1:1,1)

    iparm_(2) = 4
    iparm_(3) = 4 ! Set this value to 5 to perform iterative refinement...

    if(context%mtype==wsmp_symm) then

       ! Recall context
       call wrecallmat (context%id, info)
       if (info/=0) then
          write (0,*) 'Error, WSMP: error recalling context',context%id 
          stop
       end if

       ! Perform both ordering and symbolic factorizations
       ! They are different tasks in the symmetric version
       ! but the same task in the unsymmetric one.
       call wssmp (matrix%gr%nv, matrix%gr%ia, matrix%gr%ja , a_, d1dum, &
            &      context%perm, context%iperm, y_, matrix%gr%nv, nrhs,  &
            &      aux, naux, i1dum, iparm_, dparm_)
       ! call wssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
       !             aux, naux, mrp, context%iparm, context%dparm)
       ! n, ia, ja, avals in matrix
       ! ldb=1 nrhs=1 b in vector 
       ! perm iperm
       ! aux obsolete not used
       ! naux = 0 not used
       ! mrp(n) pivot information

       call wstoremat (context%id , info)
       if(info/=0) then
          write (0,*) 'WSMP solution: error storing context',context%id , info
          stop
       end if

    else if(context%mtype==wsmp_unsy) then

       if ( present(rmisc) ) then
          rmisc_ => rmisc
       else
          rmisc_ => d1dum
       end if

       call wgsmp (matrix%gr%nv, matrix%gr%ia, matrix%gr%ja , a_, y_, &
            &      matrix%gr%nv, nrhs, rmisc_ , iparm_, dparm_)
       ! call wgsmp (n, ia, ja, avals, b, ldb, nrhs, rmisc, iparm, dparm)
       ! n, ia, ja, avals in matrix
       ! ldb=1 nrhs=1 b in vector 
       ! rmisc(n,nrhs) is the backward error if iparm(25)=1

    end if

    if (iparm_(64) /= 0) then
       write (0,*) 'Error, WSMP: the following ERROR was detected: ', & 
            &       iparm_(64), 'during solution with several rhs'
       stop
    end if

#else
    call enable_wsmp_error_message
#endif

  end subroutine wsmp_solution_several_rhs

end module wsmp_names


