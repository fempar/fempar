
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

module field_metadata

USE types_names
USE memor_names

implicit none

#include "debug.i90"

private

    ! Type for storing field metadata
    type field_metadata_t
    private
        character(len=:), allocatable :: data_type
        character(len=:), allocatable :: name
        integer(ip)                   :: number_components
        logical                       :: filled =.false.
    contains
        procedure, non_overridable    :: set                   => field_metadata_set
        procedure, non_overridable    :: is_filled             => field_metadata_is_filled
        procedure, non_overridable    :: get_data_type         => field_metadata_get_data_type
        procedure, non_overridable    :: get_name              => field_metadata_get_name
        procedure, non_overridable    :: get_number_components => field_metadata_get_number_components
        procedure, non_overridable    :: free                  => field_metadata_free
    end type

public :: field_metadata_t

contains

    subroutine field_metadata_set(this, name, data_type, number_components)
    !-----------------------------------------------------------------
    !< Set the data of the field
    !-----------------------------------------------------------------
        class(field_metadata_t), intent(INOUT) :: this
        character(len=*) ,       intent(IN)    :: name
        character(len=*) ,       intent(IN)    :: data_type
        integer(ip) ,            intent(IN)    :: number_components
    !-----------------------------------------------------------------
        this%name              = name
        this%data_type         = data_type
        this%number_components = number_components
        this%filled            = .true.
    end subroutine field_metadata_set


    function field_metadata_is_filled(this) result(filled)
    !-----------------------------------------------------------------
    !< Check if field was filled
    !-----------------------------------------------------------------
        class(field_metadata_t), intent(IN) :: this
        logical                        :: filled
    !-----------------------------------------------------------------
        filled = this%filled
    end function field_metadata_is_filled


    function field_metadata_get_data_type(this) result(data_type)
    !-----------------------------------------------------------------
    !< Return the field data type
    !-----------------------------------------------------------------
        class(field_metadata_t),         intent(IN) :: this
        character(len=:) ,  allocatable             :: data_type
    !-----------------------------------------------------------------
        data_type = this%data_type
    end function field_metadata_get_data_type


    function field_metadata_get_name(this) result(name)
    !-----------------------------------------------------------------
    !< Return the name of the field
    !-----------------------------------------------------------------
        class(field_metadata_t),       intent(IN) :: this
        character(len=:), allocatable             :: name
    !-----------------------------------------------------------------
        name = this%name
    end function field_metadata_get_name


    function field_metadata_get_number_components(this) result(number_components)
    !-----------------------------------------------------------------
    !< Return the number of components of the field
    !-----------------------------------------------------------------
        class(field_metadata_t), intent(IN) :: this
        integer(ip)                         :: number_components
    !-----------------------------------------------------------------
        number_components = this%number_components
    end function field_metadata_get_number_components


    subroutine field_metadata_free(this)
    !-----------------------------------------------------------------
    !< Free the field metada derived type
    !-----------------------------------------------------------------
        class(field_metadata_t), intent(INOUT) :: this
    !-----------------------------------------------------------------
        if(allocated(this%name))      deallocate(this%name)
        if(allocated(this%data_type)) deallocate(this%data_type)
        this%number_components = 0
        this%filled = .false.
    end subroutine field_metadata_free

end module field_metadata
