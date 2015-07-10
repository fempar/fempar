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
module rungekutta_names
  use types_names
  use memor_names
  use time_integration_names
  implicit none
# include "debug.i90"
  private

  type rungekutta_table_t
     integer(ip)           :: &
          schem,              & ! Type of Runge-Kutta method (1: Explicit, 2: Implicit)
          stage,              & ! Number of Runge-Kutta stages
          order,              & ! Order of accuracy of the Runge-Kutta method
          flag,               & ! Flag to distinguish between same scheme, order and stage methods
          adapt                 ! Flag to set adaptive time-step coefficients
     real(rp), allocatable :: &
          A(:,:),             & ! a_{ij} coefficients of the Butcher tableau
          b(:)  ,             & ! b_i coefficients of the Butcher tableau
          c(:)  ,             & ! c_i coefficients of the Butcher tableau
          d(:)                  ! d_i coefficients of the Butcher tableau
   contains
     procedure :: create => rk_table_create
     procedure :: free   => rk_table_free
     procedure :: fill   => rk_table_fill
  end type rungekutta_table_t

  type rungekutta_table_pointer_t
     type(rungekutta_table_t), pointer :: p => NULL()
  end type rungekutta_table_pointer_t

  type rungekutta_terms_t
     integer(ip)          :: &
          ltype,             & ! Flag (0: Linear,  1:Nonlinear)
          hsite                ! Flag (0: Explicit, 1:Implicit)
   contains
     procedure :: create => rk_terms_create
  end type rungekutta_terms_t

  type, extends(time_integration_t) :: rungekutta_integrator_t
     type(rungekutta_table_pointer_t), allocatable :: rktable(:) ! Butcher tableau
     type(rungekutta_terms_t)        , allocatable :: rkterms(:) ! Terms info
     type(rungekutta_table_t)                      :: rktable_explicit
     type(rungekutta_table_t)                      :: rktable_implicit
     integer(ip)                                   :: istge = 0
   contains
     procedure :: create => rk_integ_create
     procedure :: free   => rk_integ_free
  end type rungekutta_integrator_t

  ! Type
  public :: rungekutta_integrator_t

contains

  !==================================================================================================
  subroutine rk_table_create(rktable,settb,adapt)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine creates a Butcher tableau.                                                  !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(rungekutta_table_t), intent(out) :: rktable
    integer(ip)              , intent(in)  :: settb(4)
    integer(ip), optional    , intent(in)  :: adapt
    ! Locals
    integer(ip) :: adapt_

    if(present(adapt)) then
       adapt_ = adapt
    else
       adapt_ = 0
    end if
        
    ! Fill rktable
    rktable%schem = settb(1)
    rktable%stage = settb(2)
    rktable%order = settb(3)
    rktable%flag  = settb(4)
    rktable%adapt = adapt_
    call memalloc(rktable%stage,rktable%stage,rktable%A,__FILE__,__LINE__)
    call memalloc(rktable%stage,rktable%b,__FILE__,__LINE__)
    call memalloc(rktable%stage,rktable%c,__FILE__,__LINE__)
    call memalloc(rktable%stage,rktable%d,__FILE__,__LINE__)
    
    ! Initialize A, b, c and d
    rktable%A = 0.0_rp
    rktable%b = 0.0_rp
    rktable%c = 0.0_rp
    rktable%d = 0.0_rp

  end subroutine rk_table_create

  !==================================================================================================
  subroutine rk_table_free(rktable)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine frees a Butcher tableau.                                                    !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(rungekutta_table_t), intent(inout) :: rktable

    rktable%schem = 0
    rktable%stage = 0
    rktable%order = 0
    rktable%flag  = 0
    rktable%adapt = 0
    call memfree(rktable%A,__FILE__,__LINE__)
    call memfree(rktable%b,__FILE__,__LINE__)
    call memfree(rktable%c,__FILE__,__LINE__)
    call memfree(rktable%d,__FILE__,__LINE__)

  end subroutine rk_table_free
    
  !==================================================================================================
  subroutine rk_table_fill(rktable)
    !-----------------------------------------------------------------------------------------------!
    !   This routine fill the Butcher tableau for a given RK method.                                !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(rungekutta_table_t), intent(inout) :: rktable

    ! Local variables
    real(rp) :: gamma,delta

    ! Fill Butcher tableau
    if(rktable%schem==1.and.rktable%stage==2.and.rktable%order==1.and.rktable%flag==0) then 
       ! Forward Euler (1-stage)
       rktable%A(2,1) = 1.0_rp
       rktable%b(2)   = 1.0_rp
       rktable%c(2)   = 1.0_rp
    elseif(rktable%schem==2.and.rktable%stage==2.and.rktable%order==1.and.rktable%flag==0) then 
       ! Backward Euler (1-stage)
       rktable%A(2,2) = 1.0_rp
       rktable%b(2)   = 1.0_rp
       rktable%c(2)   = 1.0_rp
    elseif(rktable%schem==1.and.rktable%stage==2.and.rktable%order==1.and.rktable%flag==1) then 
       ! Velocity Correction (1-stage)
       rktable%A(2,1) = 1.0_rp
       rktable%b(1)   = 1.0_rp
       rktable%c(2)   = 1.0_rp
    elseif(rktable%schem==2.and.rktable%stage==2.and.rktable%order==1.and.rktable%flag==1) then 
       ! Velocity Correction (1-stage)
       rktable%A(2,2) = 1.0_rp
       rktable%b(1)   = 1.0_rp
       rktable%c(2)   = 1.0_rp
    elseif(rktable%schem==1.and.rktable%stage==2.and.rktable%order==2) then 
       ! Expl. Midpoint (1-stage)
       rktable%A(2,1) = 0.5_rp
       rktable%b(2)   = 1.0_rp
       rktable%c(2)   = 0.5_rp
    elseif(rktable%schem==2.and.rktable%stage==2.and.rktable%order==2) then 
       ! Impl. Midpoint (1-stage)
       rktable%A(2,2) = 0.5_rp
       rktable%b(2)   = 1.0_rp
       rktable%c(2)   = 0.5_rp
    elseif(rktable%schem==1.and.rktable%stage==3.and.rktable%order==3.and.rktable%flag==0) then 
       ! Expl. 3rd-order (2-stage)
       gamma = (3.0_rp+sqrt(3.0_rp))/6.0_rp
       rktable%A(2,1) = gamma
       rktable%A(3,1) = gamma - 1.0_rp
       rktable%A(3,2) = 2.0_rp*(1.0_rp - gamma)
       rktable%b(2)   = 0.5_rp
       rktable%b(3)   = 0.5_rp
       rktable%c(2)   = gamma
       rktable%c(3)   = 1.0_rp - gamma
    elseif(rktable%schem==1.and.rktable%stage==3.and.rktable%order==3.and.rktable%flag==1) then 
       ! Expl. 3rd-order (2-stage)
       rktable%A(2,1) = 1.0_rp/3.0_rp
       rktable%A(3,1) = - 1.0_rp
       rktable%A(3,2) = 2.0_rp
       rktable%b(2)   = 3.0_rp/4.0_rp
       rktable%b(3)   = 1.0_rp/4.0_rp
       rktable%c(2)   = 1.0_rp/3.0_rp
       rktable%c(3)   = 1.0_rp
    elseif(rktable%schem==2.and.rktable%stage==3.and.rktable%order==3) then 
       ! Impl. 3rd-order (2-stage)
       gamma = (3.0_rp+sqrt(3.0_rp))/6.0_rp
       rktable%A(2,2) = gamma
       rktable%A(3,2) = 1.0_rp - 2.0_rp*gamma
       rktable%A(3,3) = gamma
       rktable%b(2)   = 0.5_rp
       rktable%b(3)   = 0.5_rp
       rktable%c(2)   = gamma
       rktable%c(3)   = 1.0_rp - gamma
    elseif(rktable%schem==1.and.rktable%stage==3.and.rktable%order==2.and.rktable%flag==0) then 
       ! Expl. 2nd-order, L-stab (2-stage) (flag=0)
       gamma = (2.0_rp-sqrt(2.0_rp))/2.0_rp
       delta = -2.0_rp*sqrt(2.0_rp)/3.0_rp
       rktable%A(2,1) = gamma
       rktable%A(3,1) = delta
       rktable%A(3,2) = 1.0_rp - delta
       rktable%b(2)   = 1.0_rp - gamma
       rktable%b(3)   = gamma
       rktable%c(2)   = gamma
       rktable%c(3)   = 1.0_rp
    elseif(rktable%schem==1.and.rktable%stage==3.and.rktable%order==2.and.rktable%flag==1) then 
       ! Expl. 2nd-order, L-stab (2-stage) (flag=1)
       gamma = (2.0_rp-sqrt(2.0_rp))/2.0_rp
       delta = 1.0_rp-1.0_rp/(2.0_rp*gamma)
       rktable%A(2,1) = gamma
       rktable%A(3,1) = delta
       rktable%A(3,2) = 1.0_rp - delta
       rktable%b(1)   = delta
       rktable%b(2)   = 1.0_rp - delta
       rktable%c(2)   = gamma
       rktable%c(3)   = 1.0_rp
    elseif(rktable%schem==2.and.rktable%stage==3.and.rktable%order==2.and.rktable%flag<=1) then 
       ! Impl. 2nd-order, L-stab (2-stage)
       gamma = (2.0_rp-sqrt(2.0_rp))/2.0_rp
       rktable%A(2,2) = gamma
       rktable%A(3,2) = 1.0_rp - gamma
       rktable%A(3,3) = gamma
       rktable%b(2)   = 1.0_rp - gamma
       rktable%b(3)   = gamma
       rktable%c(2)   = gamma
       rktable%c(3)   = 1.0_rp
    elseif(rktable%schem==2.and.rktable%stage==3.and.rktable%order==2.and.rktable%flag==2) then 
       ! Impl. 2nd-order (2-stage)
       rktable%A(2,2) = 1.0_rp/4.0_rp
       rktable%A(3,2) = 5.0_rp/12.0_rp
       rktable%A(3,3) = 1.0_rp/3.0_rp
       rktable%b(2)   = 0.5_rp
       rktable%b(3)   = 0.5_rp
       rktable%c(2)   = 1.0_rp/4.0_rp
       rktable%c(3)   = 3.0_rp/4.0_rp
    elseif(rktable%schem==2.and.rktable%stage==4.and.rktable%order==3) then 
       ! Impl. 3rd-order, L-stab (3-stage)
       rktable%A(2,2) = 0.4358665215_rp
       rktable%A(3,2) = 0.2820667392_rp
       rktable%A(3,3) = 0.4358665215_rp
       rktable%A(4,2) = 1.208496649_rp
       rktable%A(4,3) = -0.644363171_rp
       rktable%A(4,4) = 0.4358665215_rp
       rktable%b(2)   = 1.208496649_rp
       rktable%b(3)   = -0.644363171_rp
       rktable%b(4)   = 0.4358665215_rp
       rktable%c(2)   = 0.4358665215_rp
       rktable%c(3)   = 0.7179332608_rp
       rktable%c(4)   = 1.0_rp
       rktable%d(2)   = 0.4358665215_rp/(1.0_rp-0.4358665215_rp)
       rktable%d(3)   = (2.0_rp*0.4358665215_rp-1)/(0.4358665215_rp-1.0_rp)
    elseif(rktable%schem==1.and.rktable%stage==4.and.rktable%order==3.and.rktable%flag==0) then 
       ! Expl. 3rd-order, L-stab (3-stage)
       rktable%A(2,1) = 0.4358665215_rp
       rktable%A(3,1) = 0.3212788860_rp
       rktable%A(3,2) = 0.3966543747_rp
       rktable%A(4,1) = -0.105858296_rp
       rktable%A(4,2) = 0.5529291479_rp
       rktable%A(4,3) = 0.5529291479_rp
       rktable%b(2)   = 1.208496649_rp
       rktable%b(3)   = -0.644363171_rp
       rktable%b(4)   = 0.4358665215_rp
       rktable%c(2)   = 0.4358665215_rp
       rktable%c(3)   = 0.7179332608_rp
       rktable%c(4)   = 1.0_rp
       rktable%d(2)   = 0.4358665215_rp/(1.0_rp-0.4358665215_rp)
       rktable%d(3)   = (2.0_rp*0.4358665215_rp-1)/(0.4358665215_rp-1.0_rp)
    elseif(rktable%schem==1.and.rktable%stage==4.and.rktable%order==3.and.rktable%flag==1) then 
       ! Expl. 3rd-order, L-stab (3-stage)
       rktable%A(2,1) = 0.5_rp
       rktable%A(3,2) = 0.5_rp
       rktable%A(4,3) = 1.0_rp
       rktable%b(1)   = 1.0_rp/6.0_rp
       rktable%b(2)   = 1.0_rp/3.0_rp
       rktable%b(3)   = 1.0_rp/3.0_rp
       rktable%b(4)   = 1.0_rp/6.0_rp
       rktable%c(2)   = 0.5_rp
       rktable%c(3)   = 0.5_rp
       rktable%c(4)   = 1.0_rp
    elseif(rktable%schem==1.and.rktable%stage==4.and.rktable%order==4) then 
       ! Expl. 4rd-order, L-stab (4-stage)
       rktable%A(2,1) = 1.0_rp
       rktable%A(3,1) = 3.0_rp/8.0_rp
       rktable%A(3,2) = 1.0_rp/8.0_rp
       rktable%A(4,1) = -1.0_rp/8.0_rp
       rktable%A(4,2) = -3.0_rp/8.0_rp
       rktable%A(4,3) = 3.0_rp/2.0_rp
       rktable%b(1)   = 1.0_rp/6.0_rp
       rktable%b(2)   = 1.0_rp/18.0_rp
       rktable%b(3)   = 2.0_rp/3.0_rp
       rktable%b(4)   = 2.0_rp/9.0_rp
       rktable%c(2)   = 1.0_rp
       rktable%c(3)   = 0.5_rp
       rktable%c(4)   = 1.0_rp
    elseif(rktable%schem==2.and.rktable%stage==5.and.rktable%order==3.and.rktable%flag==0) then 
       ! Impl. 3rd-order, L-stab (4-stage)
       rktable%A(2,2) = 0.5_rp
       rktable%A(3,2) = 1.0_rp/6.0_rp
       rktable%A(3,3) = 0.5_rp
       rktable%A(4,2) = -0.5_rp
       rktable%A(4,3) = 0.5_rp
       rktable%A(4,4) = 0.5_rp
       rktable%A(5,2) = 1.5_rp
       rktable%A(5,3) = -1.5_rp
       rktable%A(5,4) = 0.5_rp
       rktable%A(5,5) = 0.5_rp
       rktable%b(2)   = 1.5_rp
       rktable%b(3)   = -1.5_rp
       rktable%b(4)   = 0.5_rp
       rktable%b(5)   = 0.5_rp
       rktable%c(2)   = 0.5_rp
       rktable%c(3)   = 2.0_rp/3.0_rp
       rktable%c(4)   = 0.5_rp
       rktable%c(5)   = 1.0_rp
    elseif(rktable%schem==1.and.rktable%stage==5.and.rktable%order==3.and.rktable%flag==0) then 
       ! Expl. 3rd-order, L-stab (4-stage)
       rktable%A(2,1) = 0.5_rp
       rktable%A(3,1) = 11.0_rp/18.0_rp
       rktable%A(3,2) = 1.0_rp/18.0_rp
       rktable%A(4,1) = 5.0_rp/6.0_rp
       rktable%A(4,2) = -5.0_rp/6.0_rp
       rktable%A(4,3) = 0.5_rp
       rktable%A(5,1) = 1.0_rp/4.0_rp
       rktable%A(5,2) = 7.0_rp/4.0_rp
       rktable%A(5,3) = 3.0_rp/4.0_rp
       rktable%A(5,4) = -7.0_rp/4.0_rp
       rktable%b(1)   = 1.0_rp/4.0_rp
       rktable%b(2)   = 7.0_rp/4.0_rp
       rktable%b(3)   = 3.0_rp/4.0_rp
       rktable%b(4)   = -7.0_rp/4.0_rp
       rktable%c(2)   = 0.5_rp
       rktable%c(3)   = 2.0_rp/3.0_rp
       rktable%c(4)   = 0.5_rp
       rktable%c(5)   = 1.0_rp
    elseif(rktable%schem==2.and.rktable%stage==5.and.rktable%order==3.and.rktable%flag==1) then 
       ! Impl. 3rd-order, BHR (5-stage)
       gamma = 424782.0_rp/974569.0_rp
       rktable%A(2,1) = gamma
       rktable%A(2,2) = gamma
       rktable%A(3,1) = gamma
       rktable%A(3,2) = -31733082319927313.0_rp/455705377221960889379854647102.0_rp
       rktable%A(3,3) = gamma
       rktable%A(4,1) = -3012378541084922027361996761794919360516301377809610.0_rp/45123394056585269977907753045030512597955897345819349.0_rp
       rktable%A(4,2) = -62865589297807153294268.0_rp/102559673441610672305587327019095047.0_rp
       rktable%A(4,3) = 418769796920855299603146267001414900945214277000.0_rp/212454360385257708555954598099874818603217167139.0_rp
       rktable%A(4,4) = gamma
       rktable%b(1)   = 487698502336740678603511.0_rp/1181159636928185920260208.0_rp
       rktable%b(2)   = 0.0_rp
       rktable%b(3)   = 302987763081184622639300143137943089.0_rp/1535359944203293318639180129368156500.0_rp
       rktable%b(4)   = -105235928335100616072938218863.0_rp/2282554452064661756575727198000.0_rp
       rktable%b(5)   = gamma
       rktable%A(5,1) = rktable%b(1)!487698502336740678603511.0_rp/1181159636928185920260208.0_rp
       rktable%A(5,2) = rktable%b(2)!0.0_rp
       rktable%A(5,3) = rktable%b(3)!0.0_rp
       rktable%A(5,4) = rktable%b(4)!302987763081184622639300143137943089.0_rp/1535359944203293318639180129368156500.0_rp
       rktable%A(5,5) = rktable%b(5)!gamma
       rktable%c(2)   = 2.0_rp*gamma
       rktable%c(3)   = 902905985686.0_rp/1035759735069.0_rp
       rktable%c(4)   = 2684624.0_rp/1147171.0_rp
       rktable%c(5)   = 1.0_rp
    elseif(rktable%schem==1.and.rktable%stage==5.and.rktable%order==3.and.rktable%flag==1) then 
       ! Expl. 3rd-order, BHR (5-stage)
       gamma = 424782.0_rp/974569.0_rp
       rktable%A(2,1) = 2.0_rp*gamma
       rktable%A(3,1) = gamma
       rktable%A(3,2) = gamma
       rktable%A(4,1) = -475883375220285986033264.0_rp/594112726933437845704163.0_rp
       rktable%A(4,2) = 0.0_rp
       rktable%A(4,3) = 1866233449822026827708736.0_rp/594112726933437845704163.0_rp
       rktable%A(5,1) = 62828845818073169585635881686091391737610308247.0_rp/176112910684412105319781630311686343715753056000.0_rp
       rktable%A(5,2) = -302987763081184622639300143137943089.0_rp/1535359944203293318639180129368156500.0_rp
       rktable%A(5,3) = 262315887293043739337088563996093207.0_rp/297427554730376353252081786906492000.0_rp
       rktable%A(5,4) = -987618231894176581438124717087.0/23877337660202969319526901856000.0_rp
       rktable%b(1)   = 487698502336740678603511.0_rp/1181159636928185920260208.0_rp
       rktable%b(2)   = 0.0_rp
       rktable%b(3)   = 302987763081184622639300143137943089.0_rp/1535359944203293318639180129368156500.0_rp
       rktable%b(4)   = -105235928335100616072938218863.0_rp/2282554452064661756575727198000.0_rp
       rktable%b(5)   = gamma
       rktable%c(2)   = 2.0_rp*gamma
       rktable%c(3)   = 902905985686.0_rp/1035759735069.0_rp
       rktable%c(4)   = 2684624.0_rp/1147171.0_rp
       rktable%c(5)   = 1.0_rp
    else
       write(*,*) 'There is no matching Runge-Kutta method for this parameters:'
       write(*,*) '   Scheme: ', rktable%schem
       write(*,*) '   Stages: ', rktable%stage
       write(*,*) '   Order : ', rktable%order
       check(.false.)
    end if
    
  end subroutine rk_table_fill

  !==================================================================================================
  subroutine rk_terms_create(rkterms,setterms)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine generates a rungekutta_terms_t type.                                        !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(rungekutta_terms_t), intent(out) :: rkterms
    integer(ip)              , intent(in)  :: setterms(2)
    
    ! Fill rkterms
    rkterms%ltype = setterms(1)  ! Linearity type (0: linear, 1: Nonlinear)
    rkterms%hsite = setterms(2)  ! Hand site of the operator (0: RHS, 1:LHS)

  end subroutine rk_terms_create

  !==================================================================================================
  subroutine rk_integ_create(rkinteg,setterms,settable,adapt)
    !-----------------------------------------------------------------------------------------------!
    !   This routine generates a rungekutta_integration_t type.                                     !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(rungekutta_integrator_t),target, intent(out) :: rkinteg
    integer(ip)                          , intent(in)  :: setterms(:,:),settable(3)
    integer(ip), optional                , intent(in)  :: adapt

    ! Local variables
    integer(ip) :: iterm,nterms
    integer(ip) :: auxterms(2),auxtb(4)

    ! Create Explicit and Implicit Runge-Kutta butcher tableaus
    auxtb(2:4) = settable(:)
    auxtb(1)   = 1
    call rkinteg%rktable_explicit%create(auxtb,adapt)
    call rkinteg%rktable_explicit%fill()
    auxtb(1)   = 2
    call rkinteg%rktable_implicit%create(auxtb,adapt)
    call rkinteg%rktable_implicit%fill()
    
    ! Allocate variables
    nterms = size(setterms,1)
    assert(size(setterms,2)==2)
    allocate(rkinteg%rkterms(nterms))
    allocate(rkinteg%rktable(nterms))
    
    ! Fill Runge-Kutta time integration type
    do iterm=1,nterms
       ! Term info
       auxterms = setterms(iterm,:)
       call rkinteg%rkterms(iterm)%create(auxterms)
       ! Point Butcher tableaus
       if(auxterms(2) == 0) then     ! Explicit term
          rkinteg%rktable(iterm)%p => rkinteg%rktable_explicit
       elseif(auxterms(2) == 1) then ! Implicit term
          rkinteg%rktable(iterm)%p => rkinteg%rktable_implicit
       end if
    end do
    
  end subroutine rk_integ_create

  !==================================================================================================
  subroutine rk_integ_free(rkinteg)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine deallocates a rungekutta_integrator_t type.                                 !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(rungekutta_integrator_t), intent(inout) :: rkinteg

    ! Local variables
    integer(ip) :: iterm

    ! Deallocate Runge-Kutta time integration subtypes
    do iterm=1,size(rkinteg%rkterms,1)
       ! Terms info
       rkinteg%rkterms(iterm)%ltype=0
       rkinteg%rkterms(iterm)%hsite=0 
       ! Butcher Tables
       rkinteg%rktable(iterm)%p => null()
    end do
    call rkinteg%rktable_explicit%free()
    call rkinteg%rktable_implicit%free()
    deallocate(rkinteg%rkterms)
    deallocate(rkinteg%rktable)
    
  end subroutine rk_integ_free

end module rungekutta_names
