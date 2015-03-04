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
module fem_precond_names
  ! Serial modules
  use types
  use memor
  use fem_graph_names
  use fem_matrix_names
  use fem_vector_names
  use pardiso_mkl_names
  use wsmp_names
  use hsl_mi20_names
  use hsl_ma87_names
  use umfpack_interface
  use umfpack_names

# include "debug.i90"

  implicit none
  private

  ! Preconditioner type
  integer(ip), parameter :: no_prec = 0
  integer(ip), parameter :: diag_prec = 1
  integer(ip), parameter :: pardiso_mkl_prec = 2
  integer(ip), parameter :: wsmp_prec = 3
  integer(ip), parameter :: hsl_mi20_prec = 4
  integer(ip), parameter :: hsl_ma87_prec = 5
  integer(ip), parameter :: umfpack_prec = 6 

!!$  ! WARNING Set a default WARNING
!!$#ifdef ENABLE_PARDISO_MKL
!!$  integer(ip), parameter :: default_prec=pardiso_mkl_prec
!!$#else
!!$#ifdef ENABLE_WSMP
!!$  integer(ip), parameter :: default_prec=wsmp_prec
!!$#else
!!$  integer(ip), parameter :: default_prec=no_prec
!!$#endif
!!$#endif
  integer(ip), parameter :: default_prec=pardiso_mkl_prec

  ! Verbosity level
  integer(ip), parameter :: no_verbose = 0
  integer(ip), parameter :: verbose = 1

  ! Release level
  integer (ip), parameter  :: precond_free_values = 7
  integer (ip), parameter  :: precond_free_struct = 8
  integer (ip), parameter  :: precond_free_clean  = 9

  type fem_precond
     ! Preconditioner type (none, diagonal, ILU, etc.)
     integer(ip)          :: type = -1 ! Undefined

     integer (ip)          :: nd1 = 0   ! Number of degrees of freedom, ndof1
     integer (ip)          :: nd2 = 0   ! Number of degrees of freedom, ndof2
     integer (ip)          :: storage   ! Storage layout (blk: block; scal: scalar)
     real(rp), allocatable :: d(:,:)    ! Inverse of main diagonal (ndof,neq) if sto=blk
                                        !                          (1   ,neq) if sto=scal

     ! Info direct solvers
     integer(ip)          :: mem_peak_symb
     integer(ip)          :: mem_perm_symb
     integer(ip)          :: nz_factors   
     integer(ip)          :: mem_peak_num 
     real(rp)             :: Mflops

     ! Info AMG preconditioners
     real(rp)    :: cs  ! Average stencil size
     real(rp)    :: cg  ! Grid complexity
     real(rp)    :: ca  ! Operator complexity
     integer(ip) :: lev ! Number of levels in the AMG hierarchy

     ! If prec_type == pardiso_mkl_prec store pardiso_mkl state
     type (pardiso_mkl_context) :: pardiso_mkl_ctxt
     integer                    :: pardiso_mkl_iparm(64)

     ! If prec_type == wsmp_prec store wsmp context
     type (wsmp_context) :: wsmp_ctxt
     integer             :: wsmp_iparm(64)
     real                :: wsmp_rparm(64)

     ! If prec_type == hsl_mi20_prec store hsl_mi20 context
!!$     integer(ip)              :: deflation = 1 
     integer(ip), allocatable :: dirichlet_nodes(:)
     type(hsl_mi20_context)   :: hsl_mi20_ctxt
     type(hsl_mi20_control)   :: hsl_mi20_ctrl ! hsl_mi20 params
     type(hsl_mi20_info)      :: hsl_mi20_info ! hsl_mi20_info
     type(hsl_mi20_data)      :: hsl_mi20_data ! hsl_mi20_data
!!$     type(fem_vector)       :: a             ! K w 
!!$     real(rp)               :: B_inv         ! (w^t K w)^{-1}  

     ! If prec_type == hsl_ma87_prec store hsl_ma87 context
     type(hsl_ma87_context) :: hsl_ma87_ctxt
     type(hsl_ma87_control) :: hsl_ma87_ctrl ! hsl_ma87 params
     type(hsl_ma87_info)    :: hsl_ma87_info ! hsl_ma87_info 

     ! If prec_type == umfpack_prec umfpack_context
     type(umfpack_context)  :: umfpack_ctxt 
     
  end type fem_precond

  type fem_precond_params
     ! The objective of this type is to define a generic 
     ! set of parameters and implement their "translation"
     ! to the specific parameter list of external libraries.
     ! It is a big TO DO.
     integer(ip) :: type      = default_prec
     integer(ip) :: verbosity = 0

!!$     ! Deflation, see Dohrmann's paper (An approximate BDDC preconditioner)
!!$     integer(ip) :: deflation = 1

     ! HSL-MI20-params (defaults consistent with HSL-MI20 defaults, except for c_fail == 2)
     real    :: st_parameter      = 0.25 
     logical :: one_pass_coarsen  = .false.
     integer :: smoother          = 2         ! 1. DJ. 2. GS.
     integer :: pre_smoothing     = 2
     integer :: post_smoothing    = 2
     integer :: v_iterations      = 1         ! number of V-cycles
     integer :: c_fail            = 1 
  end type fem_precond_params

  interface fem_precond_apply
     module procedure fem_precond_apply_vector, fem_precond_apply_r2, fem_precond_apply_r1
  end interface fem_precond_apply

  ! Constants
  public :: no_prec, diag_prec, pardiso_mkl_prec, wsmp_prec, hsl_mi20_prec, hsl_ma87_prec, umfpack_prec
  public :: precond_free_values
  public :: precond_free_struct
  public :: precond_free_clean

  ! Types
  public :: fem_precond, fem_precond_params

  ! Functions
  public :: fem_precond_create, fem_precond_free, fem_precond_symbolic, &
       &    fem_precond_numeric, fem_precond_apply, fem_precond_log_info, &
       &    fem_precond_bcast, fem_precond_fine_task,  extract_diagonal, &
            invert_diagonal, extract_diagonal_scal, apply_diagonal

contains

  !=============================================================================
  ! Dummy method required to specialize Krylov subspace methods
  subroutine fem_precond_bcast(prec,conv)
    implicit none
    type(fem_precond) , intent(in)      :: prec
    logical           , intent( inout ) :: conv
  end subroutine fem_precond_bcast

  ! Dummy method required to specialize Krylov subspace methods
  ! Needs to be filled with the abs operator machinery.
  function fem_precond_fine_task(prec)
    implicit none
    type(fem_precond) , intent(in) :: prec
    logical                        :: fem_precond_fine_task
    fem_precond_fine_task = .true. 
  end function fem_precond_fine_task

  !=============================================================================
  subroutine  fem_precond_log_info (prec)
    implicit none
    ! Parameters
    type(fem_precond)       , intent(in)          :: prec

    if(prec%type==pardiso_mkl_prec  .or. prec%type==wsmp_prec .or. prec%type==umfpack_prec ) then
       write (*,'(a,i10)') 'Peak mem.      in KBytes (symb fact) = ', prec%mem_peak_symb
       write (*,'(a,i10)') 'Permanent mem. in KBytes (symb fact) = ', prec%mem_perm_symb
       write (*,'(a,i10)') 'Peak mem.      in KBytes (num fact)  = ', prec%mem_peak_num 
       write (*,'(a,i10)') 'Size of factors (thousands)          = ', prec%nz_factors   
       write (*,'(a,f10.2)') 'MFlops for factorization             = ', prec%Mflops        
    else if (prec%type==hsl_mi20_prec) then
       write (*,'(a,f10.3)') 'Average stencil size (cS) = ', prec%cs
       write (*,'(a,f10.3)') 'Grid complexity      (cG) = ', prec%cg
       write (*,'(a,f10.3)') 'Operator complexity  (cA) = ', prec%ca
       write (*,'(a,i10)')   'Number of AMG levels      = ', prec%lev  
    else if (prec%type==hsl_ma87_prec) then
#ifdef ENABLE_HSL_MA87 
       write (*,'(a,f10.2)')   'Size of factors (thousands)          = ', & 
            real(prec%hsl_ma87_info%info%num_factor,rp)/1.0e+03   
       write (*,'(a,f10.2)') 'Number of Flops (millions)         = '  , &
            real(prec%hsl_ma87_info%info%num_flops,rp)/1.0e+06 
#endif
    end if
    
  end subroutine fem_precond_log_info
  !=============================================================================
  subroutine  fem_precond_create (mat, prec, pars)
    implicit none
    ! Parameters
    type(fem_matrix)        , intent(in)           :: mat
    type(fem_precond)       , intent(inout)        :: prec
    type(fem_precond_params), intent(in), optional :: pars

    ! Locals
    type (fem_vector) :: dum

    prec%nd1      = mat%nd1
    prec%nd2      = mat%nd2
    prec%storage  = mat%storage

    ! Save type
    if(present(pars)) then
       prec%type = pars%type
    else
       prec%type = default_prec 
    end if

    if(prec%type==pardiso_mkl_prec) then
       call pardiso_mkl ( pardiso_mkl_initialize, prec%pardiso_mkl_ctxt, &
            &             mat, dum, dum, prec%pardiso_mkl_iparm)
       prec%pardiso_mkl_iparm(18) = -1
       prec%pardiso_mkl_iparm(19) = -1

    else if (prec%type == wsmp_prec) then
       call wsmp ( wsmp_init, prec%wsmp_ctxt, mat, dum, dum, &
            &      prec%wsmp_iparm, prec%wsmp_rparm)

    else if (prec%type == hsl_mi20_prec) then
       
       call hsl_mi20 ( hsl_mi20_init, prec%hsl_mi20_ctxt, mat, dum, dum, &
            &          prec%hsl_mi20_data, prec%hsl_mi20_ctrl, prec%hsl_mi20_info )

#ifdef ENABLE_HSL_MI20 
       ! This set of instructions could may be go in
       ! a hsl_mi20_set_parameters subroutine within
       ! a module hsl_mi20_names. TODO ???
       if ( present(pars) ) then
         prec%hsl_mi20_ctrl%control%st_parameter     = pars%st_parameter
         prec%hsl_mi20_ctrl%control%one_pass_coarsen = pars%one_pass_coarsen
         prec%hsl_mi20_ctrl%control%smoother         = pars%smoother
         prec%hsl_mi20_ctrl%control%pre_smoothing    = pars%pre_smoothing
         prec%hsl_mi20_ctrl%control%post_smoothing   = pars%post_smoothing
         prec%hsl_mi20_ctrl%control%v_iterations     = pars%v_iterations
         prec%hsl_mi20_ctrl%control%c_fail           = pars%c_fail
         if ( pars%verbosity == 1 ) then
            prec%hsl_mi20_ctrl%control%print_level = 2
         end if
!!$         prec%deflation = pars%deflation
         prec%hsl_mi20_ctrl%control%error = -1
         prec%hsl_mi20_ctrl%control%print = -1
      end if
#endif
    else if (prec%type == hsl_ma87_prec) then
       call hsl_ma87 ( hsl_ma87_init, prec%hsl_ma87_ctxt, mat, dum, dum, &
            &          prec%hsl_ma87_ctrl, prec%hsl_ma87_info )
    else if (prec%type == umfpack_prec) then
       call umfpack ( umfpack_init, prec%umfpack_ctxt, mat, dum, dum)
    else if(prec%type/=no_prec .and. prec%type /= diag_prec) then
       write (0,*) 'Error: preconditioner type not supported'
       stop
    end if

  end subroutine fem_precond_create

  !=============================================================================
  subroutine fem_precond_free ( action, prec )
    implicit none

    ! Parameters
    type(fem_precond), intent(inout) :: prec
    integer(ip)      , intent(in)    :: action

    ! Locals
    type (fem_matrix) :: adum 
    type (fem_vector) :: vdum 

    if(prec%type==pardiso_mkl_prec) then

       if ( action == precond_free_clean ) then
          call pardiso_mkl ( pardiso_mkl_free_clean, prec%pardiso_mkl_ctxt, &
               &                   adum, vdum, vdum, prec%pardiso_mkl_iparm)
          return  
       end if
       if ( action == precond_free_struct  ) then
          call pardiso_mkl ( pardiso_mkl_free_struct, prec%pardiso_mkl_ctxt, &
               &                   adum, vdum, vdum, prec%pardiso_mkl_iparm)

       else if ( action == precond_free_values ) then
          call pardiso_mkl ( pardiso_mkl_free_values, prec%pardiso_mkl_ctxt, &
               &                   adum, vdum, vdum, prec%pardiso_mkl_iparm)
       end if

    else if(prec%type==wsmp_prec) then

       if ( action == precond_free_clean ) then
          call wsmp ( wsmp_free_clean, prec%wsmp_ctxt, adum, vdum, &
               &            vdum, prec%wsmp_iparm, prec%wsmp_rparm)
          return  
       end if
       if ( action == precond_free_struct  ) then
          call wsmp ( wsmp_free_struct, prec%wsmp_ctxt, adum, vdum, &
               &            vdum, prec%wsmp_iparm, prec%wsmp_rparm)

       else if ( action == precond_free_values ) then
          call wsmp ( wsmp_free_values, prec%wsmp_ctxt, adum, vdum, &
               &            vdum, prec%wsmp_iparm, prec%wsmp_rparm)
       end if

    else if(prec%type==hsl_mi20_prec) then
       if ( action == precond_free_clean ) then
          call hsl_mi20 ( hsl_mi20_free_clean, prec%hsl_mi20_ctxt, adum, vdum, vdum, &
               &          prec%hsl_mi20_data, prec%hsl_mi20_ctrl, prec%hsl_mi20_info )
          return  
       end if
       if ( action == precond_free_struct  ) then
          call hsl_mi20 ( hsl_mi20_free_struct, prec%hsl_mi20_ctxt, adum, vdum, vdum, &
               &          prec%hsl_mi20_data, prec%hsl_mi20_ctrl, prec%hsl_mi20_info )
       else if ( action == precond_free_values ) then
          call hsl_mi20 ( hsl_mi20_free_values, prec%hsl_mi20_ctxt, adum, vdum, vdum, &
               &          prec%hsl_mi20_data, prec%hsl_mi20_ctrl, prec%hsl_mi20_info )
!!$          if ( prec%deflation == 1 ) then
!!$             call fem_vector_free ( prec%a )
!!$             call memfree ( prec%dirichlet_nodes, __FILE__, __LINE__)
!!$          end if
       end if
    else if(prec%type==hsl_ma87_prec) then
       if ( action == precond_free_clean ) then
          call hsl_ma87 ( hsl_ma87_free_clean, prec%hsl_ma87_ctxt, adum, vdum, vdum, &
               &          prec%hsl_ma87_ctrl, prec%hsl_ma87_info )
          return  
       end if
       if ( action == precond_free_struct  ) then
          call hsl_ma87 ( hsl_ma87_free_struct, prec%hsl_ma87_ctxt, adum, vdum, vdum, &
               &          prec%hsl_ma87_ctrl, prec%hsl_ma87_info )
       else if ( action == precond_free_values ) then
          call hsl_ma87 ( hsl_ma87_free_values, prec%hsl_ma87_ctxt, adum, vdum, vdum, &
               &          prec%hsl_ma87_ctrl, prec%hsl_ma87_info )
       end if
    else if(prec%type==umfpack_prec) then
       if ( action == precond_free_clean ) then
          call umfpack ( umfpack_free_clean, prec%umfpack_ctxt, adum, vdum, vdum )
          return  
       end if
       if ( action == precond_free_struct  ) then
          call umfpack ( umfpack_free_struct, prec%umfpack_ctxt, adum, vdum, vdum )
       else if ( action == precond_free_values ) then
          call umfpack ( umfpack_free_values, prec%umfpack_ctxt, adum, vdum, vdum )
       end if
    else if ( prec%type == diag_prec ) then
       if ( action == precond_free_values ) then
          call memfree ( prec%d,__FILE__,__LINE__)
       end if
    else if(prec%type/=no_prec) then
       write (0,*) 'Error: preconditioner type not supported'
       stop
    end if

  end subroutine fem_precond_free

  !=============================================================================

  subroutine fem_precond_symbolic(mat, prec)
    implicit none
    ! Parameters
    type(fem_matrix)      , intent(in)    :: mat
    type(fem_precond)     , intent(inout) :: prec
    ! Locals
    type (fem_vector) :: vdum 

    if(prec%type==pardiso_mkl_prec) then
       call pardiso_mkl ( pardiso_mkl_compute_symb, prec%pardiso_mkl_ctxt, &
            &             mat, vdum, vdum, prec%pardiso_mkl_iparm )
       prec%mem_peak_symb = prec%pardiso_mkl_iparm(15)
       prec%mem_perm_symb = prec%pardiso_mkl_iparm(16)
       prec%nz_factors    = prec%pardiso_mkl_iparm(18)/1e3
    else if(prec%type==wsmp_prec) then
       call wsmp ( wsmp_compute_symb, prec%wsmp_ctxt, mat, vdum, &
            &      vdum, prec%wsmp_iparm, prec%wsmp_rparm )
       prec%mem_peak_symb = 8*prec%wsmp_iparm(23)
       prec%nz_factors    = prec%wsmp_iparm(24)
       !write(*,*) prec%wsmp_iparm(23)
       !write(*,*) prec%wsmp_iparm(24)
    else if (prec%type==hsl_mi20_prec) then
       call hsl_mi20 ( hsl_mi20_compute_symb, prec%hsl_mi20_ctxt, mat, vdum, vdum, &
            &          prec%hsl_mi20_data, prec%hsl_mi20_ctrl, prec%hsl_mi20_info )
    else if (prec%type==hsl_ma87_prec) then
       call hsl_ma87 ( hsl_ma87_compute_symb, prec%hsl_ma87_ctxt, mat, vdum, vdum, &
            &          prec%hsl_ma87_ctrl, prec%hsl_ma87_info )   
    else if (prec%type==umfpack_prec) then
       call umfpack ( umfpack_compute_symb, prec%umfpack_ctxt, mat, vdum, vdum)
#ifdef ENABLE_UMFPACK 
       prec%mem_peak_symb = (prec%umfpack_ctxt%Info(UMFPACK_SYMBOLIC_PEAK_MEMORY)*prec%umfpack_ctxt%Info(UMFPACK_SIZE_OF_UNIT))/1024.0_rp
       prec%mem_perm_symb = (prec%umfpack_ctxt%Info(UMFPACK_SYMBOLIC_SIZE)*prec%umfpack_ctxt%Info(UMFPACK_SIZE_OF_UNIT))/1024.0_rp
#endif
    else if(prec%type/=no_prec .and. prec%type /= diag_prec) then
       write (0,*) 'Error: preconditioner type not supported'
       stop
    end if

  end subroutine fem_precond_symbolic

  !=============================================================================
  subroutine fem_precond_numeric(mat, prec)
    implicit none
    ! Parameters
    type(fem_matrix)      , intent(in)    :: mat
    type(fem_precond)     , intent(inout) :: prec
    ! Locals
    type (fem_vector) :: vdum 
    integer(ip)       :: ilev, n, nnz
    integer(ip)       :: i, j
    real(rp)          :: diag
    
    if(prec%type==pardiso_mkl_prec) then
       call pardiso_mkl ( pardiso_mkl_compute_num, prec%pardiso_mkl_ctxt, &
            mat, vdum, vdum, prec%pardiso_mkl_iparm )
       prec%mem_peak_num = prec%pardiso_mkl_iparm(16)+prec%pardiso_mkl_iparm(17)
       prec%Mflops       = real(prec%pardiso_mkl_iparm(19))/1.0e3_rp
    else if(prec%type==wsmp_prec) then
       call wsmp ( wsmp_compute_num, prec%wsmp_ctxt, mat, vdum, &
            &      vdum, prec%wsmp_iparm, prec%wsmp_rparm )
       prec%mem_peak_num = 8*prec%wsmp_iparm(23)
       prec%Mflops       = real(prec%wsmp_rparm(23))/1.0e9_rp
       !write(*,*) prec%wsmp_iparm(23)
       !write(*,*) prec%wsmp_rparm(23)
    else if (prec%type==hsl_mi20_prec) then
       call hsl_mi20 ( hsl_mi20_compute_num, prec%hsl_mi20_ctxt, mat, vdum, vdum, &
            &          prec%hsl_mi20_data, prec%hsl_mi20_ctrl, prec%hsl_mi20_info )

#ifdef ENABLE_HSL_MI20 
       ! This set of instructions could may be go in
       ! a hsl_mi20_get_info subroutine within
       ! a module hsl_mi20_names. TODO ???
       prec%lev = prec%hsl_mi20_info%info%clevels
       prec%cs  = 0.0 
       prec%cg  = 0.0        
       prec%ca  = 0.0
       do ilev=1, prec%lev
          n   = prec%hsl_mi20_data%coarse_data(ilev)%A_mat%m
          nnz = prec%hsl_mi20_data%coarse_data(ilev)%A_mat%ptr(n+1)-1
          ! write (*,*) 'XXX', ilev, n, nnz, mat%gr%nv, mat%gr%ia(mat%gr%nv+1)-1
          prec%cs = prec%cs + dble(nnz)/dble(n)
          prec%cg = prec%cg + dble(n)
          prec%ca = prec%ca + dble(nnz)
       end do
       prec%cs = prec%cs/dble(prec%lev)
       prec%cg = prec%cg/dble(mat%gr%nv)
       prec%ca = prec%ca/dble(mat%gr%ia(mat%gr%nv+1)-1)
#endif

!!$       if (prec%deflation == 1) then
!!$          if (prec%storage == blk) then
!!$             ! To be implemented
!!$             assert (1 == 0)
!!$          else if (prec%storage == scal) then
!!$             ! The next portion of code only works in the case of
!!$             ! unsymmetric storage. This is precisely the case of
!!$             ! HSL_MI20
!!$             assert ( mat%gr%type == csr )
!!$             call fem_vector_alloc ( scal, mat%nd1, mat%gr%nv, prec%a )
!!$             call memalloc (  mat%gr%nv, prec%dirichlet_nodes, __FILE__, __LINE__)
!!$             prec%a%b   = 0.0_rp
!!$             prec%B_inv = 0.0_rp 
!!$             if ( associated (mat%bcs) ) then
!!$                prec%dirichlet_nodes = 0
!!$                do j=1, mat%gr%nv
!!$                   if ( mat%bcs%code(1,j) /= 1 ) then
!!$                      do i=mat%gr%ia(j), mat%gr%ia(j+1)-1
!!$                         prec%a%b(1,j) = prec%a%b(1,j) + mat%a(1,1,i)
!!$                      end do
!!$                   end if
!!$                   prec%B_inv = prec%B_inv + prec%a%b(1,j)
!!$                end do
!!$             else
!!$             prec%dirichlet_nodes = 0
!!$             do j=1, mat%gr%nv
!!$                do i=mat%gr%ia(j), mat%gr%ia(j+1)-1
!!$                   prec%a%b(1,j) = prec%a%b(1,j) + mat%a(1,1,i)
!!$                   if (j==mat%gr%ja(i)) then
!!$                      diag = mat%a(1,1,i)
!!$                   end if
!!$                end do
!!$                if ( abs(prec%a%b(1,j)-diag) < 1.0e-08 ) then
!!$                   prec%dirichlet_nodes(j) = 1
!!$                   prec%a%b(1,j) = 0.0_rp
!!$                end if
!!$                prec%B_inv = prec%B_inv + prec%a%b(1,j)
!!$             end do
!!$             end if
!!$             prec%B_inv = 1.0_rp/prec%B_inv
!!$          end if
!!$       end if

    else if (prec%type==hsl_ma87_prec) then
       call hsl_ma87 ( hsl_ma87_compute_num, prec%hsl_ma87_ctxt, mat, vdum, vdum, &
            &          prec%hsl_ma87_ctrl, prec%hsl_ma87_info )
    else if (prec%type==umfpack_prec) then
       call umfpack ( umfpack_compute_num, prec%umfpack_ctxt, mat, vdum, vdum)
#ifdef ENABLE_UMFPACK
       prec%mem_peak_num = (prec%umfpack_ctxt%Info(UMFPACK_PEAK_MEMORY)*prec%umfpack_ctxt%Info(UMFPACK_SIZE_OF_UNIT))/1024.0_rp
       prec%Mflops       = prec%umfpack_ctxt%Info(UMFPACK_FLOPS)/(1.0e+06_rp * prec%umfpack_ctxt%Info(UMFPACK_NUMERIC_TIME))
       prec%nz_factors   = (prec%umfpack_ctxt%Info(UMFPACK_UNZ)+prec%umfpack_ctxt%Info(UMFPACK_LNZ))/1.0e+03_rp
#endif
    else if (prec%type==diag_prec) then
       ! Allocate + extract
       if (prec%storage == blk) then
          call memalloc ( prec%nd1, mat%gr%nv/mat%nd1, prec%d, __FILE__,__LINE__)
          call extract_diagonal  ( mat, prec%nd1, mat%gr%nv/mat%nd1, prec%d )
       else if (prec%storage == scal) then
          call memalloc ( 1,  mat%gr%nv, prec%d, __FILE__,__LINE__)
          call extract_diagonal_scal ( mat, mat%nd1, mat%gr%nv, prec%d )
       end if

       ! Invert diagonal
       call invert_diagonal  ( mat%gr%nv*mat%nd1, prec%d )
    else if(prec%type/=no_prec) then
       write (0,*) 'Error: preconditioner type not supported'
       stop
    end if
    
  end subroutine fem_precond_numeric

  !=============================================================================
  subroutine fem_precond_apply_vector (mat, prec, x, y)
    implicit none
    ! Parameters
    type(fem_matrix)      , intent(in)    :: mat
    type(fem_precond)     , intent(inout) :: prec
    type(fem_vector)      , intent(in)    :: x
    type(fem_vector)      , intent(inout) :: y
    ! Locals
    type (fem_vector) :: vdum 
    type (fem_vector) :: E_r
    real (rp)         :: alpha, beta
    integer(ip)       :: j 

    ! write(*,*) 'Applying precond'

    if(prec%type==pardiso_mkl_prec) then
       call pardiso_mkl ( pardiso_mkl_solve, prec%pardiso_mkl_ctxt,  &
            &             mat, x, y, prec%pardiso_mkl_iparm )
    else if(prec%type==wsmp_prec) then
       call wsmp ( wsmp_solve, prec%wsmp_ctxt, mat, x, y, &
            &      prec%wsmp_iparm, prec%wsmp_rparm )
    else if(prec%type==no_prec) then
       call fem_vector_copy (x,y)
    else if ( prec%type==diag_prec ) then
       call apply_diagonal  ( mat%gr%nv*mat%nd1, prec%d, x%b, y%b )
    else if (prec%type==hsl_mi20_prec) then
!!$       if ( prec%deflation == 1 ) then
!!$          if (prec%storage == blk) then
!!$             ! To be implemented
!!$             assert (1 == 0)
!!$          else if (prec%storage == scal) then
!!$             ! B^{-1} E r 
!!$             call fem_vector_alloc ( scal, mat%nd1, mat%gr%nv, E_r )
!!$             alpha = 0.0_rp
!!$             do j=1, mat%gr%nv
!!$                alpha = alpha + x%b(1,j) 
!!$             end do
!!$             alpha = prec%B_inv*alpha
!!$
!!$             do j=1, mat%gr%nv
!!$                 E_r%b(1,j) = x%b(1,j)-prec%a%b(1,j)*alpha  
!!$             end do
!!$          end if
!!$          
!!$          call hsl_mi20 ( hsl_mi20_solve, prec%hsl_mi20_ctxt, mat, E_r, y, &
!!$               &       prec%hsl_mi20_data, prec%hsl_mi20_ctrl, prec%hsl_mi20_info )
!!$
!!$          ! beta = a^{T} * v
!!$          call fem_vector_dot ( prec%a, y, beta )
!!$
!!$          ! beta = B^{-1} * a^{T} * v
!!$          beta = prec%B_inv * beta
!!$
!!$          if (prec%storage == blk) then
!!$             ! To be implemented
!!$             assert (1 == 0)
!!$          else if (prec%storage == scal) then             
!!$             if ( associated (mat%bcs) ) then
!!$                do j=1, mat%gr%nv
!!$                   if (mat%bcs%code(1,j) /= 1) then
!!$                      y%b(1,j) = alpha +  y%b(1,j) - beta
!!$                   end if
!!$                end do
!!$             else
!!$             do j=1, mat%gr%nv
!!$                if ( prec%dirichlet_nodes(j) == 0 ) then
!!$                   y%b(1,j) = alpha +  y%b(1,j) - beta
!!$                end if
!!$             end do
!!$          end if
!!$          end if
!!$
!!$          call fem_vector_free ( E_r )

!!$       else ! deflation deactivated
          call hsl_mi20 ( hsl_mi20_solve, prec%hsl_mi20_ctxt, mat, x, y, &
               &       prec%hsl_mi20_data, prec%hsl_mi20_ctrl, prec%hsl_mi20_info )
!!$       end if
    else if (prec%type==hsl_ma87_prec) then
      call hsl_ma87 ( hsl_ma87_solve, prec%hsl_ma87_ctxt, mat, x, y, &
            &          prec%hsl_ma87_ctrl, prec%hsl_ma87_info )
    else if (prec%type==umfpack_prec) then
      call umfpack ( umfpack_solve, prec%umfpack_ctxt, mat, x, y )
    else
       write (0,*) 'Error: precondtioner type not supported'
       stop
    end if
    
  end subroutine fem_precond_apply_vector

  !=============================================================================
  subroutine fem_precond_apply_r2 (mat, prec, nrhs, x, ldx, y, ldy)
    implicit none
    ! Parameters
    type(fem_matrix)      , intent(in)    :: mat
    type(fem_precond)     , intent(inout) :: prec
    integer(ip)       , intent(in)        :: nrhs, ldx, ldy
    real(rp)          , intent(in)        :: x (ldx, nrhs)
    real(rp)          , intent(inout)     :: y (ldy, nrhs)

    ! Locals
    type (fem_vector)      :: vdum
    integer(ip)            :: i, j 
    real(rp) , allocatable :: E_r(:,:) 
    real (rp), allocatable :: alpha(:), beta(:)

    if(prec%type==pardiso_mkl_prec) then
       call pardiso_mkl ( pardiso_mkl_solve, prec%pardiso_mkl_ctxt,  &
            &             mat, nrhs, x, ldX, y, ldY, prec%pardiso_mkl_iparm )
    else if(prec%type==wsmp_prec) then
       ! AFM : I did not modify wsmp interface in such
       ! a way that it is able to handle non-contiguous
       ! 2D arrays. PENDING!!!
       assert ( mat%gr%nv == ldx )
       assert ( mat%gr%nv == ldy )
       call wsmp ( wsmp_solve, prec%wsmp_ctxt, mat, nrhs, x, y, &
            &      prec%wsmp_iparm, prec%wsmp_rparm )
    else if(prec%type==no_prec) then
       do i=1, nrhs
          y(1:mat%gr%nv,i) = x(1:mat%gr%nv,i)
       end do
    else if(prec%type==diag_prec) then
       do i=1, nrhs
          call apply_diagonal ( mat%gr%nv*mat%nd1, prec%d, x(1,i), y(1,i) )
       end do
    else if (prec%type==hsl_mi20_prec) then
!!$       if ( prec%deflation == 1 ) then
!!$          if (prec%storage == blk) then
!!$             ! To be implemented
!!$             assert (1 == 0)
!!$          else if (prec%storage == scal) then
!!$             call memalloc ( mat%gr%nv, nrhs, E_r, __FILE__, __LINE__)
!!$             call memalloc ( nrhs, alpha, __FILE__, __LINE__)
!!$             call memalloc ( nrhs, beta , __FILE__, __LINE__)
!!$
!!$             alpha = 0.0_rp
!!$             do i=1, nrhs
!!$                do j=1, mat%gr%nv
!!$                   alpha(i) = alpha(i) + x(j,i) 
!!$                end do
!!$                alpha(i) = prec%B_inv*alpha(i)
!!$             end do
!!$
!!$             do j=1, nrhs
!!$                E_r(:,j) = x(:,j)-prec%a%b(1,:)*alpha(j)  
!!$             end do
!!$
!!$          end if
!!$
!!$          call hsl_mi20 ( hsl_mi20_solve, prec%hsl_mi20_ctxt, mat, nrhs, E_r, ldx, y, ldy, &
!!$               &          prec%hsl_mi20_data, prec%hsl_mi20_ctrl, prec%hsl_mi20_info )
!!$
!!$          ! beta = a^{T} * v
!!$          beta = 0.0_rp
!!$          do i=1, nrhs
!!$             do j=1, mat%gr%nv
!!$                beta(i) = beta(i) + prec%a%b(1,j)*y(j,i)
!!$             end do
!!$          end do
!!$          
!!$          ! beta = B^{-1} * a^{T} * v
!!$          beta = prec%B_inv * beta
!!$
!!$          if (prec%storage == blk) then
!!$             ! To be implemented
!!$             assert (1 == 0)
!!$          else if (prec%storage == scal) then
!!$             do i=1,nrhs
!!$                do j=1, mat%gr%nv
!!$                   if ( prec%dirichlet_nodes(j) == 0 ) then
!!$                      y(j,i) = alpha(i) +  y(j,i) - beta(i)
!!$                   end if
!!$                end do
!!$             end do
!!$          end if
!!$
!!$
!!$          call memfree ( E_r )
!!$          call memfree ( alpha )
!!$          call memfree ( beta )
!!$
!!$       else
          call hsl_mi20 ( hsl_mi20_solve, prec%hsl_mi20_ctxt, mat, nrhs, x, ldx, y, ldy, &
               &          prec%hsl_mi20_data, prec%hsl_mi20_ctrl, prec%hsl_mi20_info )
          
!!$       end if

    else if (prec%type==hsl_ma87_prec) then
       call hsl_ma87 ( hsl_ma87_solve, prec%hsl_ma87_ctxt, mat, nrhs, x, ldx, y, ldy, &
            &          prec%hsl_ma87_ctrl, prec%hsl_ma87_info )
    else if (prec%type==umfpack_prec) then
       call umfpack ( umfpack_solve, prec%umfpack_ctxt, mat, nrhs, x, ldx, y, ldy)
    else
       write (0,*) 'Error: precondtioner type not supported'
       stop
    end if

  end subroutine fem_precond_apply_r2

  !=============================================================================
  subroutine fem_precond_apply_r1 (mat, prec, x, y)
    implicit none
    ! Parameters
    type(fem_matrix) , intent(in)    :: mat
    type(fem_precond), intent(inout) :: prec
    real(rp)         , intent(in)    :: x (mat%gr%nv)
    real(rp)         , intent(inout) :: y (mat%gr%nv)

    ! Locals
    type (fem_vector)     :: vdum 
    real(rp), allocatable :: E_r(:) 
    real (rp)             :: alpha, beta
    integer(ip)           :: j 

    if(prec%type==pardiso_mkl_prec) then
       call pardiso_mkl ( pardiso_mkl_solve, prec%pardiso_mkl_ctxt,  &
            &             mat,  x, y, prec%pardiso_mkl_iparm )
    else if(prec%type==wsmp_prec) then
       call wsmp ( wsmp_solve, prec%wsmp_ctxt, mat, x, y, &
            &      prec%wsmp_iparm, prec%wsmp_rparm )
    else if(prec%type==no_prec) then
       y=x
    else if(prec%type==diag_prec) then
       call apply_diagonal ( mat%gr%nv*mat%nd1, prec%d, x, y )
    else if (prec%type==hsl_mi20_prec) then
!!$       if ( prec%deflation == 1 ) then
!!$          if (prec%storage == blk) then
!!$             ! To be implemented
!!$             assert (1 == 0)
!!$          else if (prec%storage == scal) then
!!$             call memalloc ( mat%gr%nv, E_r, __FILE__, __LINE__)
!!$
!!$             alpha = 0.0_rp
!!$             do j=1, mat%gr%nv
!!$                alpha = alpha + x(j) 
!!$             end do
!!$             alpha = prec%B_inv*alpha
!!$             
!!$             do j=1, mat%gr%nv
!!$                E_r(j) = x(j)-prec%a%b(1,j)*alpha  
!!$             end do
!!$
!!$          end if
!!$
!!$          call hsl_mi20 ( hsl_mi20_solve, prec%hsl_mi20_ctxt, mat, E_r, y, &
!!$               &       prec%hsl_mi20_data, prec%hsl_mi20_ctrl, prec%hsl_mi20_info )
!!$
!!$          ! beta = a^{T} * v
!!$          beta = 0.0_rp
!!$          do j=1, mat%gr%nv
!!$             beta = beta + prec%a%b(1,j)*y(j)
!!$          end do
!!$          
!!$          ! beta = B^{-1} * a^{T} * v
!!$          beta = prec%B_inv * beta
!!$
!!$          if (prec%storage == blk) then
!!$             ! To be implemented
!!$             assert (1 == 0)
!!$          else if (prec%storage == scal) then
!!$             do j=1, mat%gr%nv
!!$                if ( prec%dirichlet_nodes(j) == 0 ) then
!!$                   y(j) = alpha +  y(j) - beta   
!!$                end if
!!$             end do
!!$          end if
!!$           
!!$          call memfree ( E_r ) 
!!$       else
          call hsl_mi20 ( hsl_mi20_solve, prec%hsl_mi20_ctxt, mat,  x, y, &
               &          prec%hsl_mi20_data, prec%hsl_mi20_ctrl, prec%hsl_mi20_info )
!!$       end if
    else if (prec%type==hsl_ma87_prec) then
       call hsl_ma87 ( hsl_ma87_solve, prec%hsl_ma87_ctxt, mat,  x, y, &
            &          prec%hsl_ma87_ctrl, prec%hsl_ma87_info )
    else if (prec%type==umfpack_prec) then
       call umfpack ( umfpack_solve, prec%umfpack_ctxt, mat,  x, y )
    else
       write (0,*) 'Error: precondtioner type not supported'
       stop
    end if

  end subroutine fem_precond_apply_r1

    ! Auxiliary routines
  subroutine invert_diagonal (n, d)
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: n  
    real(rp)   , intent(inout) :: d(n)

    ! Locals
    integer(ip) :: i
    do i=1,n 
       d(i) = 1/d(i)
    end do
  end subroutine invert_diagonal

  subroutine apply_diagonal (n,d,x,y)
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: n
    real(rp)   , intent(in)    :: d(n)
    real(rp)   , intent(in)    :: x(n)
    real(rp)   , intent(inout) :: y(n)

    ! Locals
    integer(ip)                :: i
    do i = 1 ,n 
       y(i) = x(i) * d(i)
    end do

  end subroutine apply_diagonal

  subroutine extract_diagonal (f_mat, d_nd, d_nv, d)
    implicit none
    ! Parameters
    type(fem_matrix), intent(in)  :: f_mat
    integer(ip)     , intent(in)  :: d_nd, d_nv
    real(rp)        , intent(out) :: d(d_nd, d_nv)

    if( f_mat%type == css_mat ) then
      call extract_diagonal_css (f_mat%symm,f_mat%nd1,f_mat%nd2,f_mat%gr%nv,f_mat%d,d_nd,d_nv,d)
    else if( f_mat%type == csr_mat ) then
      call extract_diagonal_csr (f_mat%symm,f_mat%nd1,f_mat%nd2,f_mat%gr%nv,f_mat%gr%nv2,f_mat%gr%ia,f_mat%gr%ja,f_mat%a,d_nd,d_nv,d)
    else if( f_mat%type == csc_mat ) then
      call extract_diagonal_csc (f_mat%symm,f_mat%nd1,f_mat%nd2,f_mat%gr%nv,f_mat%gr%nv2,f_mat%gr%ia,f_mat%gr%ja,f_mat%a,d_nd,d_nv,d)
    end if

  end subroutine extract_diagonal

  subroutine extract_diagonal_scal (f_mat, d_nd, d_nv, d)
    implicit none
    ! Parameters
    type(fem_matrix), intent(in)  :: f_mat
    integer(ip)     , intent(in)  :: d_nd, d_nv
    real(rp)        , intent(out) :: d(1, d_nv)

    if ( f_mat%type == css_mat ) then
      call extract_diagonal_css_scal (f_mat%symm,f_mat%nd1,f_mat%nd2,f_mat%gr%nv,f_mat%d,d_nd,d_nv,d)
    else if ( f_mat%type == csr_mat ) then
      call extract_diagonal_csr_scal (f_mat%symm,f_mat%nd1,f_mat%nd2,f_mat%gr%nv,f_mat%gr%nv2,f_mat%gr%ia,f_mat%gr%ja,f_mat%a,d_nd,d_nv,d)
    else if ( f_mat%type == csc_mat ) then
      call extract_diagonal_csc_scal (f_mat%symm,f_mat%nd1,f_mat%nd2,f_mat%gr%nv,f_mat%gr%nv2,f_mat%gr%ia,f_mat%gr%ja,f_mat%a,d_nd,d_nv,d)
    end if

  end subroutine extract_diagonal_scal

  subroutine extract_diagonal_css (ks,nd1,nd2,nv,da,d_nd,d_nv,d)
    implicit none
    ! Parameters
    integer(ip), intent(in)  :: ks,nd1,nd2,nv,d_nd,d_nv
    real(rp)   , intent(in)  :: da(nd1,nd2,nv)
    real(rp)   , intent(out) :: d(d_nd,d_nv)

    write (0,*) 'Error: the body of extract_diagonal_css in par_precond.f90 still to be written'
    write (0,*) 'Error: volunteers are welcome !!!'
    stop
  end subroutine extract_diagonal_css

  subroutine extract_diagonal_css_scal (ks,nd1,nd2,nv,da,d_nd,d_nv,d)
    implicit none
    ! Parameters
    integer(ip), intent(in)  :: ks,nd1,nd2,nv,d_nd,d_nv
    real(rp)   , intent(in)  :: da(1,1,nv)
    real(rp)   , intent(out) :: d (1,d_nv)

    ! Locals
    integer(ip) :: iv
    do iv = 1, d_nv
       d(1,iv) =  da(1,1,iv)
    end do
  end subroutine extract_diagonal_css_scal

  subroutine extract_diagonal_csr (ks,nd1,nd2,nv,nv2,ia,ja,a,d_nd,d_nv,d)
    implicit none
    ! Parameters
    integer(ip), intent(in)  :: ks,nd1,nd2,nv,nv2,d_nd,d_nv
    integer(ip), intent(in)  :: ia(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: a(nd2,nd1,ia(nv+1)-1)
    real(rp)   , intent(out) :: d(d_nd,d_nv)

    write (0,*) 'Error: the body of extract_diagonal_csr in par_precond.f90 still to be written'
    write (0,*) 'Error: volunteers are welcome !!!'
    stop

  end subroutine extract_diagonal_csr

  subroutine extract_diagonal_csr_scal (ks,nd1,nd2,nv,nv2,ia,ja,a,d_nd,d_nv,d)
    implicit none
    ! Parameters
    integer(ip), intent(in)  :: ks,nd1,nd2,nv,nv2,d_nd,d_nv
    integer(ip), intent(in)  :: ia(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: a(1,1,ia(nv+1)-1)
    real(rp)   , intent(out) :: d(1,d_nv)

    ! Locals
    integer(ip)              :: iv, iz, of, izc, ivc


    if(ks==1) then                     ! Unsymmetric 
       do iv = 1, nv, nd1
          iz   = ia(iv)
          of   = 0
          do while( ja(iz) /= iv )
             iz = iz + 1
             of = of + 1
          end do

          do ivc = iv, iv + nd1 - 1
             izc      = ia(ivc) + of
             d(1,ivc) = a(1,1,izc) 
             of       = of + 1
          end do ! ivc

       end do ! iv
    else if (ks==0) then                  ! Symmetric
       do iv = 1, nv
          izc     = ia(iv)
          assert(ja(izc)==iv)
          d(1,iv) = a(1,1,izc)
       end do ! iv
    end if

  end subroutine extract_diagonal_csr_scal

  subroutine extract_diagonal_csc (ks,nd1,nd2,nv,nv2,ia,ja,a,d_nd,d_nv,d)
    implicit none
    ! Parameters
    integer(ip), intent(in)  :: ks,nd1,nd2,nv,nv2,d_nd,d_nv
    integer(ip), intent(in)  :: ia(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: a(nd2,nd1,ia(nv+1)-1)
    real(rp)   , intent(out) :: d(d_nd,d_nv)

    write (0,*) 'Error: the body of extract_diagonal_csc in par_precond.f90 still to be written'
    write (0,*) 'Error: volunteers are welcome !!!'
    stop

  end subroutine extract_diagonal_csc

  subroutine extract_diagonal_csc_scal (ks,nd1,nd2,nv,nv2,ia,ja,a,d_nd,d_nv,d)
    implicit none
    ! Parameters
    integer(ip), intent(in)  :: ks,nd1,nd2,nv,nv2,d_nd,d_nv
    integer(ip), intent(in)  :: ia(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: a(1,1,ia(nv+1)-1)
    real(rp)   , intent(out) :: d(1,d_nv)

    write (0,*) 'Error: the body of extract_diagonal_csc_scal in par_precond.f90 still to be written'
    write (0,*) 'Error: volunteers are welcome !!!'
    stop

  end subroutine extract_diagonal_csc_scal

end module fem_precond_names
