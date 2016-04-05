
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

module vtk_field

USE types_names
USE memor_names

implicit none

#include "debug.i90"

private

    ! Type for storing field descriptors
    type vtk_field_t
    private
        character(len=:), allocatable :: var_location              ! 'Node' or 'Cell' field
        character(len=:), allocatable :: var_name                  ! Name of the field
        character(len=:), allocatable :: field_type                ! Field data type 'Float32', 'Float64', 'Int32', etc.
        integer(ip)                   :: number_components = 0     ! Number of components
        logical                       :: filled = .false.          ! Field data is already filled
    contains
    private
        procedure, non_overridable, public :: free => vtk_field_free
    end type vtk_field_t

public :: vtk_field_t

contains

    subroutine vtk_field_free(this)
    !-----------------------------------------------------------------
    !< Free the vtk_field_t derived type
    !-----------------------------------------------------------------
        class(vtk_field_t), intent(inout) :: this
        integer(ip)                       :: error
    !-----------------------------------------------------------------
        if(allocated(this%var_location)) deallocate(this%var_location, stat=error)
        assert(error==0)
        if(allocated(this%var_name))     deallocate(this%var_name, stat=error)
        assert(error==0)
        if(allocated(this%field_type))   deallocate(this%field_type, stat=error)
        assert(error==0)
        this%number_components = 0
        this%filled            = .false.
    end subroutine

end module vtk_field
