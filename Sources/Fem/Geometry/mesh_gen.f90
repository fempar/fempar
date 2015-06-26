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
module fem_mesh_gen_names
  use types_names
  use memor_names
!  use fem_conditions_names
# include "debug.i90"
  implicit none
  private

  integer(ip), parameter :: periodic = -1

  type mesh_size_t
     integer(ip)                :: &
        lverb=0,                   &         ! Verbosity level
        ntdix=0,                   &         ! Type of discretization in x (0=uniform, 1=cubic, 2=tanh, 3=imh+unif, 4:imh+tanh)
        ntdiy=0,                   &         ! Type of discretization in y (0=uniform, 1=cubic, 2=tanh, 3=imh+unif, 4:imh+tanh)
        ntdiz=0,                   &         ! Type of discretization in z (0=uniform, 1=cubic, 2=tanh, 3=imh+unif, 4:imh+tanh)
        pdegr=1,                   &         ! Interpolation order p (1=linear, 2=quadratic, 3=cubic) 
        mater=0,                   &         ! Material case
        ierrc=0,                   &         ! Error code
        neblx=0,                   &         ! Number of elements in the x boundary layer
        nebly=0,                   &         ! Number of elements in the y boundary layer
        neblz=0                              ! Number of elements in the z boundary layer

     real(rp)                   :: &
        xleng   = 1.0_rp,          &         ! Size of the domain in x
        yleng   = 1.0_rp,          &         ! Size of the domain in y
        zleng   = 1.0_rp,          &         ! Size of the domain in z
        zx1     = 0.1_rp,          &         ! size of the elements at x=0   (left)  
        zx2     = 0.1_rp,          &         ! size of the elements at x=a/2 (center)
        zy1     = 0.1_rp,          &         ! size of the elements at y=0   (bottom)
        zy2     = 0.1_rp,          &         ! size of the elements at y=b/2 (center)
        x0      = 0.0_rp,          &         ! Origin x-coordinate
        y0      = 0.0_rp,          &         ! Origin y-coordinate
        z0      = 0.0_rp,          &         ! Origin z-coordinate
        atol    = 1.0e-8_rp,       &         ! Absolute tolerance
        stret   = 0.0_rp,          &         ! Stretching parameter (old)
        xstret  = 2.75_rp,         &         ! Stretching parameter
        ystret  = 2.75_rp,         &         ! Stretching parameter
        zstret  = 2.75_rp,         &         ! Stretching parameter
        xlengbl = 0.0_rp,          &         ! Size of the boundary layer in x
        ylengbl = 0.0_rp,          &         ! Size of the boundary layer in y
        zlengbl = 0.0_rp                     ! Size of the boundary layer in z
  end type mesh_size_t

  ! Types
  public :: mesh_size_t

end module fem_mesh_gen_names

