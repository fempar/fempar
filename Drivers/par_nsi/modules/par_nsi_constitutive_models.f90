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
module par_nsi_constitutive_models_names
  use fempar_names
  implicit none
# include "debug.i90"
  private

  real(rp), parameter         :: E  = 1.0_rp
  real(rp), parameter         :: nu = 0.2_rp
  real(rp), parameter, public :: lambda = (nu*E)/((1+nu)*(1-2*nu))
  real(rp), parameter, public :: mu     = E/(2*(1+nu))
  real(rp), parameter, public :: inv_K = 1.0_rp/(lambda + 2*mu/3)
  real(rp), parameter, public :: one_third = 1.0_rp/3.0_rp

end module par_nsi_constitutive_models_names
