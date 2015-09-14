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
module preconditioner_names
  ! Serial modules
  use types_names
  use memor_names
  use graph_names
  use matrix_names
  use serial_scalar_array_names
  use pardiso_mkl_names
  use wsmp_names
  use hsl_mi20_names
  use hsl_ma87_names
  use umfpack_interface_names
  use umfpack_names

  ! Abstract modules
  use abstract_vector_names
  use abstract_operator_names

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
  integer (ip), parameter  :: preconditioner_free_values = 7
  integer (ip), parameter  :: preconditioner_free_struct = 8
  integer (ip), parameter  :: preconditioner_free_clean  = 9

  type, extends(abstract_operator_t) :: preconditioner_t
     ! Preconditioner type (none, diagonal, ILU, etc.)
     integer(ip)          :: type = -1 ! Undefined

     real(rp), allocatable :: d(:)    ! Inverse of main diagonal (neq)

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
     type (pardiso_mkl_context_t), pointer :: pardiso_mkl_ctxt
     integer                   , pointer :: pardiso_mkl_iparm(:)

     ! If prec_type == wsmp_prec store wsmp context
     type (wsmp_context_t), pointer :: wsmp_ctxt
     integer            , pointer :: wsmp_iparm(:)
     real               , pointer :: wsmp_rparm(:)

     ! If prec_type == hsl_mi20_prec store hsl_mi20 context
     type(hsl_mi20_context_t), pointer   :: hsl_mi20_ctxt
     type(hsl_mi20_control_t), pointer   :: hsl_mi20_ctrl ! hsl_mi20 params
     type(hsl_mi20_info_t)   , pointer   :: hsl_mi20_info ! hsl_mi20_info
     type(hsl_mi20_data_t)   , pointer   :: hsl_mi20_data ! hsl_mi20_datat

     ! If prec_type == hsl_ma87_prec store hsl_ma87 context
     type(hsl_ma87_context_t), pointer :: hsl_ma87_ctxt
     type(hsl_ma87_control_t), pointer :: hsl_ma87_ctrl ! hsl_ma87 params
     type(hsl_ma87_info_t)   , pointer :: hsl_ma87_info ! hsl_ma87_info

     ! If prec_type == umfpack_prec umfpack_context_t
     type(umfpack_context_t), pointer  :: umfpack_ctxt 

     ! AFM: I had to add a pointer to the linear system coefficient matrix within the
     !      preconditioner. The linear coefficient matrix is no longer passed to 
     !      precond%apply(r,z) in the abstract implementation of Krylov subspace methods, but
     !      it is still required.
     type(matrix_t), pointer :: mat

   contains
     procedure :: apply => preconditioner_apply_tbp
     procedure :: apply_fun => preconditioner_apply_fun_tbp
     procedure :: free => preconditioner_free_tbp
  end type preconditioner_t

  type preconditioner_params_t
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
     real    :: damping           = 0.8       ! Damping factor for smoother==DJ

  end type preconditioner_params_t

  interface preconditioner_apply
     module procedure preconditioner_apply_vector, preconditioner_apply_r2, preconditioner_apply_r1
  end interface preconditioner_apply

  ! Constants
  public :: no_prec, diag_prec, pardiso_mkl_prec, wsmp_prec, hsl_mi20_prec, hsl_ma87_prec, umfpack_prec
  public :: preconditioner_free_values
  public :: preconditioner_free_struct
  public :: preconditioner_free_clean

  ! Types
  public :: preconditioner_t, preconditioner_params_t

  ! Functions
  public :: preconditioner_create, preconditioner_free, preconditioner_symbolic, &
       &    preconditioner_numeric, preconditioner_apply, preconditioner_log_info, &
       &    preconditioner_bcast, preconditioner_fine_task,  extract_diagonal, &
            invert_diagonal, apply_diagonal

contains

  !=============================================================================
  ! Dummy method required to specialize Krylov subspace methods
  subroutine preconditioner_bcast(prec,conv)
    implicit none
    type(preconditioner_t) , intent(in)      :: prec
    logical           , intent( inout ) :: conv
  end subroutine preconditioner_bcast

  ! Dummy method required to specialize Krylov subspace methods
  ! Needs to be filled with the abs operator machinery.
  function preconditioner_fine_task(prec)
    implicit none
    type(preconditioner_t) , intent(in) :: prec
    logical                        :: preconditioner_fine_task
    preconditioner_fine_task = .true. 
  end function preconditioner_fine_task

  !=============================================================================
  subroutine  preconditioner_log_info (prec)
    implicit none
    ! Parameters
    type(preconditioner_t)       , intent(in)          :: prec

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
    
  end subroutine preconditioner_log_info
  !=============================================================================
  subroutine  preconditioner_create (mat, prec, pars)
    implicit none
    ! Parameters
    type(matrix_t)        , target, intent(in)           :: mat
    type(preconditioner_t)       , intent(inout)        :: prec
    type(preconditioner_params_t), intent(in), optional :: pars

    ! Locals
    type (serial_scalar_array_t) :: dum

    prec%mat => mat

    ! Save type
    if(present(pars)) then
       prec%type = pars%type
    else
       prec%type = default_prec 
    end if

    if(prec%type==pardiso_mkl_prec) then
       allocate(prec%pardiso_mkl_ctxt)
       allocate(prec%pardiso_mkl_iparm(64))
       call pardiso_mkl ( pardiso_mkl_initialize, prec%pardiso_mkl_ctxt, &
            &             mat, dum, dum, prec%pardiso_mkl_iparm)
       prec%pardiso_mkl_iparm(18) = -1
       prec%pardiso_mkl_iparm(19) = -1

    else if (prec%type == wsmp_prec) then
       allocate(prec%wsmp_ctxt)
       allocate(prec%wsmp_iparm(64))
       allocate(prec%wsmp_rparm(64))
       call wsmp ( wsmp_init, prec%wsmp_ctxt, mat, dum, dum, &
            &      prec%wsmp_iparm, prec%wsmp_rparm)

    else if (prec%type == hsl_mi20_prec) then
       allocate(prec%hsl_mi20_ctxt)
       allocate(prec%hsl_mi20_data)
       allocate(prec%hsl_mi20_ctrl)
       allocate(prec%hsl_mi20_info)
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
         prec%hsl_mi20_ctrl%control%damping          = pars%damping

         if ( pars%verbosity == 1 ) then
            prec%hsl_mi20_ctrl%control%print_level = 2
         end if
         prec%hsl_mi20_ctrl%control%error = -1
         prec%hsl_mi20_ctrl%control%print = -1
      end if
#endif
    else if (prec%type == hsl_ma87_prec) then
       allocate(prec%hsl_ma87_ctxt)
       allocate(prec%hsl_ma87_ctrl)
       allocate(prec%hsl_ma87_info)
       call hsl_ma87 ( hsl_ma87_init, prec%hsl_ma87_ctxt, mat, dum, dum, &
            &          prec%hsl_ma87_ctrl, prec%hsl_ma87_info )
    else if (prec%type == umfpack_prec) then
       allocate(prec%umfpack_ctxt)
       call umfpack ( umfpack_init, prec%umfpack_ctxt, mat, dum, dum)
    else if(prec%type/=no_prec .and. prec%type /= diag_prec) then
       write (0,*) 'Error: preconditioner type not supported'
       check(1==0)
    end if
  end subroutine preconditioner_create

  !=============================================================================
  subroutine preconditioner_free ( action, prec )
    implicit none

    ! Parameters
    type(preconditioner_t), intent(inout) :: prec
    integer(ip)      , intent(in)    :: action

    ! Locals
    type (matrix_t) :: adum 
    type (serial_scalar_array_t) :: vdum 

    if ( action == preconditioner_free_clean ) then
       nullify(prec%mat)
    end if

    if(prec%type==pardiso_mkl_prec) then
       if ( action == preconditioner_free_clean ) then
          call pardiso_mkl ( pardiso_mkl_free_clean, prec%pardiso_mkl_ctxt, &
               &                   adum, vdum, vdum, prec%pardiso_mkl_iparm)
          deallocate(prec%pardiso_mkl_ctxt)
          deallocate(prec%pardiso_mkl_iparm)
          return  
       end if
       if ( action == preconditioner_free_struct  ) then
          call pardiso_mkl ( pardiso_mkl_free_struct, prec%pardiso_mkl_ctxt, &
               &                   adum, vdum, vdum, prec%pardiso_mkl_iparm)
       else if ( action == preconditioner_free_values ) then
          call pardiso_mkl ( pardiso_mkl_free_values, prec%pardiso_mkl_ctxt, &
               &                   adum, vdum, vdum, prec%pardiso_mkl_iparm)
       end if

    else if(prec%type==wsmp_prec) then
       if ( action == preconditioner_free_clean ) then
          call wsmp ( wsmp_free_clean, prec%wsmp_ctxt, adum, vdum, &
               &            vdum, prec%wsmp_iparm, prec%wsmp_rparm)
          deallocate(prec%wsmp_ctxt)
          deallocate(prec%wsmp_iparm)
          deallocate(prec%wsmp_rparm)
          return  
       end if
       if ( action == preconditioner_free_struct  ) then
          call wsmp ( wsmp_free_struct, prec%wsmp_ctxt, adum, vdum, &
               &            vdum, prec%wsmp_iparm, prec%wsmp_rparm)

       else if ( action == preconditioner_free_values ) then
          call wsmp ( wsmp_free_values, prec%wsmp_ctxt, adum, vdum, &
               &            vdum, prec%wsmp_iparm, prec%wsmp_rparm)
       end if

    else if(prec%type==hsl_mi20_prec) then
       if ( action == preconditioner_free_clean ) then
          call hsl_mi20 ( hsl_mi20_free_clean, prec%hsl_mi20_ctxt, adum, vdum, vdum, &
               &          prec%hsl_mi20_data, prec%hsl_mi20_ctrl, prec%hsl_mi20_info )
          deallocate(prec%hsl_mi20_ctxt)
          deallocate(prec%hsl_mi20_data)
          deallocate(prec%hsl_mi20_ctrl)
          deallocate(prec%hsl_mi20_info)
          return  
       end if
       if ( action == preconditioner_free_struct  ) then
          call hsl_mi20 ( hsl_mi20_free_struct, prec%hsl_mi20_ctxt, adum, vdum, vdum, &
               &          prec%hsl_mi20_data, prec%hsl_mi20_ctrl, prec%hsl_mi20_info )
       else if ( action == preconditioner_free_values ) then
          call hsl_mi20 ( hsl_mi20_free_values, prec%hsl_mi20_ctxt, adum, vdum, vdum, &
               &          prec%hsl_mi20_data, prec%hsl_mi20_ctrl, prec%hsl_mi20_info )
       end if
    else if(prec%type==hsl_ma87_prec) then
       if ( action == preconditioner_free_clean ) then
          call hsl_ma87 ( hsl_ma87_free_clean, prec%hsl_ma87_ctxt, adum, vdum, vdum, &
               &          prec%hsl_ma87_ctrl, prec%hsl_ma87_info )
          deallocate(prec%hsl_ma87_ctxt)
          deallocate(prec%hsl_ma87_ctrl)
          deallocate(prec%hsl_ma87_info)
          return  
       end if
       if ( action == preconditioner_free_struct  ) then
          call hsl_ma87 ( hsl_ma87_free_struct, prec%hsl_ma87_ctxt, adum, vdum, vdum, &
               &          prec%hsl_ma87_ctrl, prec%hsl_ma87_info )
       else if ( action == preconditioner_free_values ) then
          call hsl_ma87 ( hsl_ma87_free_values, prec%hsl_ma87_ctxt, adum, vdum, vdum, &
               &          prec%hsl_ma87_ctrl, prec%hsl_ma87_info )
       end if
    else if(prec%type==umfpack_prec) then
       if ( action == preconditioner_free_clean ) then
          call umfpack ( umfpack_free_clean, prec%umfpack_ctxt, adum, vdum, vdum )
          deallocate(prec%umfpack_ctxt)
          return  
       end if
       if ( action == preconditioner_free_struct  ) then
          call umfpack ( umfpack_free_struct, prec%umfpack_ctxt, adum, vdum, vdum )
       else if ( action == preconditioner_free_values ) then
          call umfpack ( umfpack_free_values, prec%umfpack_ctxt, adum, vdum, vdum )
       end if
    else if ( prec%type == diag_prec ) then
       if ( action == preconditioner_free_values ) then
          call memfree ( prec%d,__FILE__,__LINE__)
       end if
    else if(prec%type/=no_prec) then
       write (0,*) 'Error: preconditioner type not supported'
       check(1==0)
    end if

  end subroutine preconditioner_free

  !=============================================================================
  subroutine preconditioner_symbolic(mat, prec)
    implicit none
    ! Parameters
    type(matrix_t)      , intent(in), target    :: mat
    type(preconditioner_t)     , intent(inout) :: prec
    ! Locals
    type (serial_scalar_array_t) :: vdum 

    prec%mat => mat

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
       check(1==0)
    end if

  end subroutine preconditioner_symbolic

  !=============================================================================
  subroutine preconditioner_numeric( prec)
    implicit none
    ! Parameters
    type(preconditioner_t), target, intent(inout) :: prec
    ! Locals
    type(matrix_t), pointer :: mat
    type (serial_scalar_array_t) :: vdum 
    integer(ip)       :: ilev, n, nnz
    integer(ip)       :: i, j
    real(rp)          :: diag
    
    !prec%mat => mat
    mat => prec%mat

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
       call memalloc ( mat%gr%nv, prec%d, __FILE__,__LINE__)
       call extract_diagonal(mat%gr%symmetric_storage,mat%gr%nv,mat%gr%ia,mat%gr%ja,mat%a,mat%gr%nv,prec%d)
       ! Invert diagonal
       call invert_diagonal  ( mat%gr%nv, prec%d )
    else if(prec%type/=no_prec) then
       write (0,*) 'Error: preconditioner type not supported'
       check(1==0)
    end if
    
  end subroutine preconditioner_numeric

  !=============================================================================
  subroutine preconditioner_apply_vector (mat, prec, x, y)
    implicit none
    ! Parameters
    type(matrix_t)      , intent(in)    :: mat
    type(preconditioner_t)     , intent(inout) :: prec
    type(serial_scalar_array_t)      , intent(in)    :: x
    type(serial_scalar_array_t)      , intent(inout) :: y
    ! Locals
    type (serial_scalar_array_t) :: vdum 
    type (serial_scalar_array_t) :: E_r
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
       call serial_scalar_array_copy (x,y)
    else if ( prec%type==diag_prec ) then
       call apply_diagonal  ( mat%gr%nv, prec%d, x%b, y%b )
    else if (prec%type==hsl_mi20_prec) then
          call hsl_mi20 ( hsl_mi20_solve, prec%hsl_mi20_ctxt, mat, x, y, &
               &       prec%hsl_mi20_data, prec%hsl_mi20_ctrl, prec%hsl_mi20_info )
    else if (prec%type==hsl_ma87_prec) then
      call hsl_ma87 ( hsl_ma87_solve, prec%hsl_ma87_ctxt, mat, x, y, &
            &          prec%hsl_ma87_ctrl, prec%hsl_ma87_info )
    else if (prec%type==umfpack_prec) then
      call umfpack ( umfpack_solve, prec%umfpack_ctxt, mat, x, y )
    else
       write (0,*) 'Error: preconditioner type not supported'
       check(1==0)
    end if
    
  end subroutine preconditioner_apply_vector

  !=============================================================================
  subroutine preconditioner_apply_r2 (mat, prec, nrhs, x, ldx, y, ldy)
    implicit none
    ! Parameters
    type(matrix_t)      , intent(in)    :: mat
    type(preconditioner_t)     , intent(inout) :: prec
    integer(ip)       , intent(in)        :: nrhs, ldx, ldy
    real(rp)          , intent(in)        :: x (ldx, nrhs)
    real(rp)          , intent(inout)     :: y (ldy, nrhs)

    ! Locals
    type (serial_scalar_array_t)      :: vdum
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
          call apply_diagonal ( mat%gr%nv, prec%d, x(1,i), y(1,i) )
       end do
    else if (prec%type==hsl_mi20_prec) then
          call hsl_mi20 ( hsl_mi20_solve, prec%hsl_mi20_ctxt, mat, nrhs, x, ldx, y, ldy, &
               &          prec%hsl_mi20_data, prec%hsl_mi20_ctrl, prec%hsl_mi20_info )
    else if (prec%type==hsl_ma87_prec) then
       call hsl_ma87 ( hsl_ma87_solve, prec%hsl_ma87_ctxt, mat, nrhs, x, ldx, y, ldy, &
            &          prec%hsl_ma87_ctrl, prec%hsl_ma87_info )
    else if (prec%type==umfpack_prec) then
       call umfpack ( umfpack_solve, prec%umfpack_ctxt, mat, nrhs, x, ldx, y, ldy)
    else
       write (0,*) 'Error: precondtioner type not supported'
       check(1==0)
    end if

  end subroutine preconditioner_apply_r2

  !=============================================================================
  subroutine preconditioner_apply_r1 (mat, prec, x, y)
    implicit none
    ! Parameters
    type(matrix_t) , intent(in)    :: mat
    type(preconditioner_t), intent(inout) :: prec
    real(rp)         , intent(in)    :: x (mat%gr%nv)
    real(rp)         , intent(inout) :: y (mat%gr%nv)

    ! Locals
    type (serial_scalar_array_t)     :: vdum 
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
       call apply_diagonal ( mat%gr%nv, prec%d, x, y )
    else if (prec%type==hsl_mi20_prec) then
          call hsl_mi20 ( hsl_mi20_solve, prec%hsl_mi20_ctxt, mat,  x, y, &
               &          prec%hsl_mi20_data, prec%hsl_mi20_ctrl, prec%hsl_mi20_info )
    else if (prec%type==hsl_ma87_prec) then
       call hsl_ma87 ( hsl_ma87_solve, prec%hsl_ma87_ctxt, mat,  x, y, &
            &          prec%hsl_ma87_ctrl, prec%hsl_ma87_info )
    else if (prec%type==umfpack_prec) then
       call umfpack ( umfpack_solve, prec%umfpack_ctxt, mat,  x, y )
    else
       write (0,*) 'Error: precondtioner type not supported'
       check(1==0)
    end if

  end subroutine preconditioner_apply_r1

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

  subroutine extract_diagonal (graph_symmetric_storage,nv,ia,ja,a,d_nv,d)
    implicit none
    ! Parameters
	logical    , intent(in)  :: graph_symmetric_storage
    integer(ip), intent(in)  :: nv,d_nv
    integer(ip), intent(in)  :: ia(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: a(ia(nv+1)-1)
    real(rp)   , intent(out) :: d(d_nv)

    ! Locals
    integer(ip)              :: iv, iz, of, izc, ivc

    if(.not. graph_symmetric_storage) then
       do iv = 1, nv
          iz   = ia(iv)
          of   = 0
          do while( ja(iz) /= iv )
             iz = iz + 1
             of = of + 1
          end do

          do ivc = iv, iv 
             izc      = ia(ivc) + of
             d(ivc) = a(izc) 
             of       = of + 1
          end do ! ivc
       end do ! iv
    else
       do iv = 1, nv
          izc     = ia(iv)
          assert(ja(izc)==iv)
          d(iv) = a(izc)
       end do ! iv
    end if

  end subroutine extract_diagonal

  !=============================================================================
  subroutine preconditioner_apply_tbp (op, x, y)
    implicit none
    ! Parameters
    class(preconditioner_t)    , intent(in)    :: op
    class(abstract_vector_t)   , intent(in)    :: x
    class(abstract_vector_t)   , intent(inout) :: y
    
    assert (associated(op%mat))

    call x%GuardTemp()

    select type(x)
    class is (serial_scalar_array_t)
       select type(y)
       class is(serial_scalar_array_t)
          if(op%type==pardiso_mkl_prec) then
             call pardiso_mkl ( pardiso_mkl_solve, op%pardiso_mkl_ctxt,  &
                  &             op%mat, x, y, op%pardiso_mkl_iparm ) 
          else if(op%type==wsmp_prec) then
             call wsmp ( wsmp_solve, op%wsmp_ctxt, op%mat, x, y, &
                  &      op%wsmp_iparm, op%wsmp_rparm )
          else if(op%type==no_prec) then
             call y%copy(x)
          else if ( op%type==diag_prec) then
             call apply_diagonal  ( op%mat%gr%nv, op%d, x%b, y%b )
          else if (op%type==hsl_mi20_prec) then
             call hsl_mi20 ( hsl_mi20_solve, op%hsl_mi20_ctxt, op%mat, x, y, &
                  &       op%hsl_mi20_data, op%hsl_mi20_ctrl, op%hsl_mi20_info )
          else if (op%type==hsl_ma87_prec) then
             call hsl_ma87 ( hsl_ma87_solve, op%hsl_ma87_ctxt, op%mat, x, y, &
                  &          op%hsl_ma87_ctrl, op%hsl_ma87_info )
          else if (op%type==umfpack_prec) then
             call umfpack ( umfpack_solve, op%umfpack_ctxt, op%mat, x, y )
          else
             write (0,*) 'Error: preconditioner type not supported'
             check(1==0)
          end if
       class default
          write(0,'(a)') 'matrix_t%apply: unsupported y class'
          check(1==0)
       end select
    class default
       write(0,'(a)') 'preconditioner_t%apply: unsupported x class'
       check(1==0)
    end select

    call x%CleanTemp()
  end subroutine preconditioner_apply_tbp

  !=============================================================================
  function preconditioner_apply_fun_tbp (op, x) result(y)
    implicit none
    ! Parameters
    class(preconditioner_t), intent(in)   :: op
    class(abstract_vector_t), intent(in)  :: x
    class(abstract_vector_t), allocatable :: y
    type(serial_scalar_array_t), allocatable :: local_y

    
    assert (associated(op%mat))

    call x%GuardTemp()

    select type(x)
    class is (serial_scalar_array_t)
       allocate(local_y)
       call serial_scalar_array_alloc ( op%mat%gr%nv, local_y)
       call op%apply(x, local_y)
       call move_alloc(local_y, y)
       call y%SetTemp()
    class default
       write(0,'(a)') 'preconditioner_t%apply_fun: unsupported x class'
       check(1==0)
    end select

    call x%CleanTemp()
  end function preconditioner_apply_fun_tbp

  subroutine preconditioner_free_tbp(this)
    implicit none
    class(preconditioner_t), intent(inout) :: this
  end subroutine preconditioner_free_tbp

end module preconditioner_names
