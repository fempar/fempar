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
!=============================================================================
!
! Base solver class containing parameters, control type and control functions.
!
!=============================================================================
module solver_base_names
use types_names
use memor_names

  ! List of convergence criteria available for iterative solvers 
  integer(ip), parameter :: res_nrmgiven_rhs_nrmgiven  = 1  ! ||  r(i) ||g <= rtol*||  b    ||g + atol 
  integer(ip), parameter :: res_nrmgiven_res_nrmgiven  = 2  ! ||  r(i) ||g <= rtol*||  r(0) ||g + atol   
  integer(ip), parameter :: delta_rhs                  = 3  ! || dx(i) ||  <= rtol*||  b  || + atol
  integer(ip), parameter :: delta_delta                = 4  ! || dx(i) ||  <= rtol*||dx(1)|| + atol
  integer(ip), parameter :: res_res                    = 5  ! ||  r(i) ||  <= rtol*|| r(0)|| + atol
  integer(ip), parameter :: res_rhs                    = 6  ! ||  r(i) ||  <= rtol*||  b  || + atol
  integer(ip), parameter :: delta_rhs_and_res_res      = 7  ! delta_rhs    AND res_res
  integer(ip), parameter :: delta_rhs_and_res_rhs      = 8  ! delta_rhs    AND res_rhs
  integer(ip), parameter :: delta_delta_and_res_res    = 9  ! delta_delta  AND res_res
  integer(ip), parameter :: delta_delta_and_res_rhs    = 10 ! delta_delta  AND res_rhs 
                                                            ! ||.|| is the 2-norm, dx(i) = x(i) - x(i-1),
                                                            ! r(i) is the residual at the i-th iteration


  integer (ip), parameter :: mgs  = 1 ! mgs : Modified Gram-Schmidt 
                                      !       (appropriate for serial GMRES)
  integer (ip), parameter :: icgs = 2 ! icgs: Iterative Classical Gram-Schmidt 
                                      !       (appropriate for distributed GMRES)


  logical     , parameter :: default_track_conv_his  = .false.
  real    (rp), parameter :: default_relax           = 1.0_rp
  real    (rp), parameter :: default_atol            = 0.0_rp
  real    (rp), parameter :: default_rtol            = 1.0e-08_rp
  integer (ip), parameter :: default_luout           = 6
  integer (ip), parameter :: default_stopc           = res_nrmgiven_res_nrmgiven
  integer (ip), parameter :: default_trace           = 0
  integer (ip), parameter :: default_itmax           = 1000
  integer (ip), parameter :: default_dkrymax         = 30
  integer (ip), parameter :: default_orto            = icgs
  integer (ip), parameter :: default_nrhs            = 1

  ! List of Krylov subspace methods available
  integer (ip), parameter :: cg               = 1
  integer (ip), parameter :: lgmres           = 2
  integer (ip), parameter :: rgmres           = 3
  integer (ip), parameter :: fgmres           = 4
  integer (ip), parameter :: icg              = 7
  integer (ip), parameter :: lfom             = 8  ! Left preconditioned Full Orthogonalization Method
  integer (ip), parameter :: minres           = 9  ! Preconditioned MINimal RESidual method

  ! Not actually Krylov methods
  integer (ip), parameter :: richard          = 5  ! Richardson (fixed-point iteration)
  integer (ip), parameter :: direct           = 6  ! Apply preconditioner directly


  integer (ip), parameter :: default_kry_meth = direct

  character(len=*), parameter  :: methdname(9)=(/ &
       & 'PCG       ', &
       & 'LGMRES    ', &
       & 'RGMRES    ', &
       & 'FGMRES    ', &
       & 'RICHARDSON', &
       & 'DIRECT    ', &
       & 'IPCG      ', &
       & 'LFOM      ', &
       & 'MINRES    ' /)

  ! Control data (always initialized to default values)
  type solver_control_t
     ! Is the solver machinery going to track the convergence history?
     logical :: track_conv_his = default_track_conv_his

     ! Input parameters
     real(rp)      :: relax    = default_relax              ! Relaxation parameter
     real(rp)      :: rtol     = default_rtol               ! Relative tolerance
     real(rp)      :: atol     = default_atol               ! Absolute tolerance
     integer(ip)   :: dkrymax  = default_dkrymax            ! Maximum dimension of the Krylov subspace
     integer(ip)   :: orto     = default_orto               ! Ortonormalization strategy (mgs,icgs)
     integer(ip)   :: luout    = default_luout              ! Logical unit to output
     integer(ip)   :: stopc    = default_stopc              ! Stopping criteria
     integer(ip)   :: trace    = default_trace              ! Message every trace iterations
     integer(ip)   :: itmax    = default_itmax              ! Max. # of iterations
     integer(ip)   :: method   = default_kry_meth           ! Krylov subspace method
     integer(ip)   :: nrhs     = default_nrhs               ! Number of simultaneous RHS
     ! Output parameters
     logical               :: converged = .false.           ! Converged?
     integer(ip)           :: it        = 0                 ! # of iterations to converge
     real(rp)              :: bn2       = 0.0_rp            ! RHS 2-norm
     real(rp)              :: rn2       = 0.0_rp            ! Residual 2-norm
     real(rp)              :: dxn2      = 0.0_rp            ! dx 2-norm
     real(rp)              :: tol1      = default_atol      ! Tolerance 1
     real(rp)              :: tol2      = default_atol      ! Tolerance 2
     real(rp)              :: err1      = 0.0_rp            ! Current error estimate 1
     real(rp)              :: err2      = 0.0_rp            ! Current error estimate 2
     real(rp), allocatable :: err1h(:)                      ! Error estimates 1 history
     real(rp), allocatable :: err2h(:)                      ! Error estimates 2 history
  end type solver_control_t

contains

  subroutine solver_control_allocate_conv_his( ctrl )
    implicit none
    type(solver_control_t), intent(inout) :: ctrl

    if (ctrl%track_conv_his) then
       call memalloc(ctrl%itmax,ctrl%err1h,__FILE__,__LINE__)
       call memalloc(ctrl%itmax,ctrl%err2h,__FILE__,__LINE__)
       ctrl%err1h=0.0_rp
       ctrl%err2h=0.0_rp
    end if
  end subroutine solver_control_allocate_conv_his

  subroutine solver_control_free_conv_his( ctrl )
    implicit none
    type(solver_control_t), intent(inout) :: ctrl

    if (ctrl%track_conv_his) then
       call memfree(ctrl%err1h,__FILE__,__LINE__)
       call memfree(ctrl%err2h,__FILE__,__LINE__)
    end if
  end subroutine solver_control_free_conv_his

  subroutine solver_control_log_header( ctrl )
    implicit none 

    ! Parameters
    type(solver_control_t), intent(in) :: ctrl

    ! Local variables
    character(len=*), parameter    :: fmt1='(a18,1x,a4,3(2x,a15))'
    character(len=*), parameter    :: fmt2='(a18,1x,a4,3(2x,a15),3(2x,a15))'
    integer, parameter             :: outlen=18 
    character(len=len(methdname))  :: mname
    character(len=outlen)          :: outname

    mname = adjustl(trim(methdname(ctrl%method)))
    write(outname,'(a)') mname(1:min(len_trim(mname),outlen-1))//':'

    select case(ctrl%stopc)

    case ( delta_rhs, delta_delta, res_res, res_rhs, res_nrmgiven_rhs_nrmgiven, &
         & res_nrmgiven_res_nrmgiven )

       write(ctrl%luout,fmt1) adjustl(outname),'Iteration','Error Estimate','Tolerance'

    case ( delta_rhs_and_res_res, delta_rhs_and_res_rhs,  & 
         delta_delta_and_res_res, delta_delta_and_res_rhs )

       write(ctrl%luout,fmt2) adjustl(outname), 'Iteration', 'Error Estimate', 'Tolerance', &
            & 'Error Estimate', 'Tolerance' 

    case default
       ! Write an error message and stop ?      
    end select
  end subroutine solver_control_log_header

  subroutine solver_control_log_conv ( ctrl )
    implicit none 

    ! Parameters
    type(solver_control_t), intent(in) :: ctrl

    ! Local variables
    character(len=*), parameter   :: fmt1='(a18,1x,i4,3(2x,es16.9))'
    character(len=*), parameter   :: fmt2='(a18,1x,i4,3(2x,es16.9),3(2x,es16.9))'
    integer, parameter            :: outlen=18 
    character(len=len(methdname)) :: mname
    character(len=outlen)         :: outname

    if ( (mod(ctrl%it,ctrl%trace) == 0).or.ctrl%converged.or.(ctrl%it>=ctrl%itmax)) then 
       mname = adjustl(trim(methdname(ctrl%method)))
       write(outname,'(a)') mname(1:min(len_trim(mname),outlen-1))//':'
       select case(ctrl%stopc)
       case ( delta_rhs, delta_delta, res_res, res_rhs, & 
            & res_nrmgiven_rhs_nrmgiven, res_nrmgiven_res_nrmgiven )
          write(ctrl%luout,fmt1) adjustl(outname), ctrl%it, ctrl%err1, ctrl%tol1
       case ( delta_rhs_and_res_res  , delta_rhs_and_res_rhs, & 
            & delta_delta_and_res_res, delta_delta_and_res_rhs )
          write(ctrl%luout,fmt2) adjustl(outname), ctrl%it, ctrl%err1, ctrl%tol1, &
               &  ctrl%err2, ctrl%tol2 
       case default
          ! Write an error message and stop ?
       end select
    endif

  end subroutine solver_control_log_conv

  subroutine solver_control_log_end ( ctrl )
    implicit none
    ! Parameters
    type(solver_control_t), intent(in) :: ctrl

    character(len=*), parameter  :: fmt11='(a,2x,es16.9,1x,a,1x,i4,1x,a)'
    character(len=*), parameter  :: fmt12='(a,3(2x,es16.9))'

    character(len=*), parameter  :: fmt21='(a,2x,es16.9,1x,es16.9,1x,a,1x,i4,1x,a)'
    character(len=*), parameter  :: fmt22='(a,3(2x,es16.9),3(2x,es16.9))'


    select case( ctrl%stopc )

    case ( delta_rhs,delta_delta,res_res,res_rhs,&
         & res_nrmgiven_rhs_nrmgiven, res_nrmgiven_res_nrmgiven)

       if ( ctrl%converged ) then
          write(ctrl%luout,fmt11) trim(methdname(ctrl%method))//' converged to ', &
               & ctrl%tol1,' in ',ctrl%it,' iterations. '
          write(ctrl%luout,fmt12) 'Last iteration error estimate: ', ctrl%err1
       else
          write(ctrl%luout,fmt11) trim(methdname(ctrl%method))//' failed to converge to ', &
               & ctrl%tol1,' in ',ctrl%it,' iterations. '
          write(ctrl%luout,fmt12) 'Last iteration error estimate: ', ctrl%err1 
       end if

    case ( delta_rhs_and_res_res  , delta_rhs_and_res_rhs, & 
         & delta_delta_and_res_res, delta_delta_and_res_rhs )

       if ( ctrl%converged ) then
          write(ctrl%luout,fmt21) trim(methdname(ctrl%method))//' converged to ', &
               & ctrl%tol1, ctrl%tol2, ' in ', ctrl%it ,' iterations. '
          write(ctrl%luout,fmt22) 'Last iteration error estimates: ', ctrl%err1, ctrl%err2
       else
             write(ctrl%luout,fmt21) trim(methdname(ctrl%method))//' failed to converge to ', &
                  & ctrl%tol1, ctrl%tol2, ' in ', ctrl%it ,' iterations. '
             write(ctrl%luout,fmt22) 'Last iteration error estimates: ', ctrl%err1, ctrl%err2
       end if

    case default
       ! Write an error message and stop ?      

    end select

  end subroutine solver_control_log_end

  subroutine solver_control_log_conv_his ( ctrl )
    implicit none 

    ! Parameters
    type(solver_control_t), intent(in) :: ctrl

    ! Local variables
    character(len=*), parameter   :: fmt1='(a18,1x,i4,3(2x,es16.9))'
    character(len=*), parameter   :: fmt2='(a18,1x,i4,3(2x,es16.9),3(2x,es16.9))'
    integer, parameter            :: outlen=18 
    character(len=len(methdname)) :: mname
    character(len=outlen)         :: outname
    integer(ip)                   :: i

    if (ctrl%track_conv_his) then
       call solver_control_log_header(ctrl)

       mname = adjustl(trim(methdname(ctrl%method)))
       write(outname,'(a)') mname(1:min(len_trim(mname),outlen-1))//':'
       select case(ctrl%stopc)
       case ( delta_rhs, delta_delta, res_res, res_rhs, & 
            & res_nrmgiven_rhs_nrmgiven, res_nrmgiven_res_nrmgiven )
          do i=1,ctrl%it
             write(ctrl%luout,fmt1) adjustl(outname), i, ctrl%err1h(i), ctrl%tol1
          end do
       case ( delta_rhs_and_res_res  , delta_rhs_and_res_rhs, & 
            & delta_delta_and_res_res, delta_delta_and_res_rhs )
          do i=1,ctrl%it
             write(ctrl%luout,fmt2) adjustl(outname), i, ctrl%err1h(i), ctrl%tol1, &
                  &                                      ctrl%err2h(i), ctrl%tol2 
          end do
       case default
          ! Write an error message and stop ?
       end select

       call solver_control_log_end(ctrl)
    end if

  end subroutine solver_control_log_conv_his

  !Taken from SPARSKIT
  subroutine givens(x,y,c,s)
    real(rp) x,y,c,s

    !     Given x and y, this subroutine generates a Givens' rotation c, s.
    !     And apply the rotation on (x,y) ==> (sqrt(x**2 + y**2), 0).
    !     (See P 202 of "matrix computations" by Golub and van Loan.)

    real(rp) t,one,rzero
    parameter (rzero=0.0_rp,one=1.0_rp)

    if (x.eq.rzero .and. y.eq.rzero) then
       c = one
       s = zero
    else if (abs(y).gt.abs(x)) then
       t = x / y
       x = sqrt(one+t*t)
       s = sign(one / x, y)
       c = t*s
    else if (abs(y).le.abs(x)) then
       t = y / x
       y = sqrt(one+t*t)
       c = sign(one / y, x)
       s = t*c
    else
       ! X or Y must be an invalid floating-point number, set both to zero
       x = rzero
       y = rzero
       c = one
       s = rzero
    endif
    x = abs(x*y)
    y = rzero
  end subroutine givens

end module solver_base_names

