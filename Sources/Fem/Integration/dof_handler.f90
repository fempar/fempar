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
module dof_handler_names
  use types
  use memor

  implicit none
# include "debug.i90"
  private
  ! Data needed to go from geometrical info (import, partition, mesh) to 
  ! a dof info (dof_import, dof_partition, dof_mesh)

 !  type dof_block

!       integer(ip) :: nvars                    ! Total number of variables per block

!      integer(ip), allocatable ::     &
! !          ndofs,                     &       ! Total number of dofs (nprob)
!           nvars_prob(:)                       ! Number of vars per problem vars_prob(nprob)

!      type(array_ip2) ::              &
!           vars_prob(:)                        ! Variables per problem

!   end type dof_block

  type physical_problem
     integer(ip)        ::           &
          nvars                              ! Number of different problems
     integer(ip), allocatable ::     &
          l2g_var(:)                    ! Order chosen for variables (size nvars)
  end type physical_problem

  type dof_handler

     integer(ip) ::     &
          nblocks,                   &       ! Number of blocks
          nprobs,                     &       ! Number of problems
          nvars_global                      ! Total number of different physical variables
          !nvars_prob                         ! Number of physical variables per problem

     integer(ip), allocatable ::     &
!          l2g_vars                   &       ! Local (in block) to global Id of variable
!          phys_prob(:),              &       ! List of physical problems (nprob) 
          dof_coupl(:,:),            &        ! Dof_coupling(nvar,nvar) for avoiding allocation & assembly of zero blocks
          vars_block(:)                      ! Parameter per unknown (size nvars)

!     type(dof_handler), allocatable :: dof_blocks(:)

     type(physical_problem), allocatable :: problems(:)
!     type(fem_blocks)  ::           &
!          blocks                             ! blocks ordering

  end type dof_handler

  ! Types
  public :: dof_handler

  ! Functions
  !public :: dof_handler_print, dof_handler_create, dof_handler_fill, dof_handler_free

end module dof_handler_names
