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
# include "debug.i90"
module nsi_names
 use types
 use memor
 use problem_names
 implicit none
 private 

 type, extends(physical_problem) :: nsi_problem
    integer(ip) ::   & 
         ksnsi,      & ! Symmetry flag (+-1) (1: no simetrica definida, -1: simetrica no definida)
         kfl_conv,   & ! Flag for enabling advection (Stokes=0; NavSto=1)
         kfl_symg,   & ! Flag for symmetric grad-grad term
         kfl_tder,   & ! Flag for time derivative computation
         kfl_skew,   & ! Flag for enabling skewsymmetric convective terms (Off=0; type1=1, type2=2)
         kfl_vort,   & ! Flag for vorticity computation
         case_veloc, & ! Exact velocity
         case_press, & ! Exact pressure
         case_tempo, & ! Exact temporal solution
         case_t_pre    ! Exact temporal solution for pressure
    real(rp) ::         &
         react,         & ! Reaction
         diffu,         & ! Diffusion
         par_veloc(30), & ! Exact velocity field components (vector)
         par_press(10), & ! Exact pressure field components (scalar)
         par_tempo(3),  & ! Exact temporal field components (scalar)
         gravi(3)         ! Gravity field
    contains
      procedure :: create => nsi_create
 end type nsi_problem

 public :: nsi_problem

contains

  !=================================================================================================
  subroutine nsi_create(prob,ndime)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine contains definitions of the Navier-Stokes problem                          !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_problem), intent(out) :: prob
    integer(ip)       , intent(in)  :: ndime
    integer(ip) :: i,istat

    ! Fill default problem data
    prob%ndime = ndime
    prob%ntens = prob%ndime*(prob%ndime+1)/2 ! Number of tensor components 
    prob%nunks = 2                           ! Velocity and pressure
    prob%nvars = prob%ndime+1                ! Number of degrees of freedom

    call memalloc(prob%nunks,prob%vars_of_unk,__FILE__,__LINE__)
    prob%vars_of_unk(1) = ndime
    prob%vars_of_unk(2) = 1

    !prob%vars_of_unk(2) = prob%vars_of_unk(1) + 1
    !prob%vars_of_unk(3) = prob%vars_of_unk(2) + 1

    ! Should be overwritten by the driver in a multiphyiscs context
    call memalloc(prob%nvars,prob%l2g_var,__FILE__,__LINE__)
    do i = 1,prob%nvars
       prob%l2g_var(i) = i
    end do

    ! allocate(prob%unkno_names(nunks),stat=istat)
    ! check(istat==0)
    ! prob%unkno_names(1)='velocity'
    ! prob%unkno_names(2)='pressure'

    ! Flags
    prob%ksnsi = 1    ! Symmetry flag
    prob%kfl_conv = 1 ! Enabling advection
    prob%kfl_vort = 0 ! Vorticity not computed
    prob%kfl_symg = 0 ! Symmetric grad-grad term (On=1, Off=0)
    prob%kfl_tder = 0 ! Time derivative not computed 
    prob%kfl_skew = 0 ! Enabling skewsymmetric convective terms: Off

    ! Problem variables
    prob%react  = 0.0_rp  ! Reaction
    prob%diffu  = 1.0_rp  ! Diffusion
    prob%gravi  = 0.0_rp  ! Gravity field

    ! Analytical field variables
    prob%case_veloc   = 0      ! Velocity field (see exact.f90)
    prob%case_press   = 0      ! Pressure field
    prob%case_tempo   = 0      ! Temporal field
    prob%case_t_pre   = 0      ! Temporal field for pressure
    prob%par_veloc(:) = 0.0_rp ! Exact velocity field
    prob%par_press(:) = 0.0_rp ! Exact pressure field
    prob%par_tempo(:) = 0.0_rp ! Exact temporal field
    
  end subroutine nsi_create

end module nsi_names
