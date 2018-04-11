!   Copyright (C) 2014 Santiago Badia, Alberto F. Martín and Javier Principe
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
module output_handler_field_generator_names
use types_names
use output_handler_patch_names
use output_handler_fe_cell_function_names
use field_names
USE reference_fe_names,          only: field_type_scalar, field_type_vector, field_type_tensor, field_type_symmetric_tensor
implicit none
#include "debug.i90"
private

    type, abstract :: output_handler_field_generator_t
    private
    contains
      procedure(get_field_type_interface)      , deferred :: get_field_type
      procedure(generate_patch_field_interface), deferred :: generate_patch_field
    end type
    
    abstract interface
       function get_field_type_interface(this) result(field_type)
         import :: output_handler_field_generator_t
         implicit none
         class(output_handler_field_generator_t), intent(in)    :: this
         character(len=:), allocatable :: field_type
       end function get_field_type_interface
       
       subroutine generate_patch_field_interface(this, oh_cell_fe_function, visualization_points_coordinates, patch_field)
         import :: output_handler_field_generator_t, & 
                   output_handler_fe_cell_function_t, & 
                   output_handler_patch_field_t, &
                   point_t
         implicit none          
         class(output_handler_field_generator_t), intent(in)    :: this
         type(output_handler_fe_cell_function_t), intent(in)    :: oh_cell_fe_function
         type(point_t)                          , intent(in)    :: visualization_points_coordinates(:)
         type(output_handler_patch_field_t)     , intent(inout) :: patch_field
       end subroutine generate_patch_field_interface
    end interface
    
    type :: output_handler_field_generator_info_t
    private
        character(len=:)                       ,  allocatable :: name
        class(output_handler_field_generator_t), pointer      :: field_generator
    contains
    private
        procedure, non_overridable         ::                          oh_field_generator_info_assign
        procedure, non_overridable, public :: create                => oh_field_generator_info_create
        procedure, non_overridable, public :: get_name              => oh_field_generator_info_get_name
        procedure, non_overridable, public :: get_field_type        => oh_field_generator_info_get_field_type
        procedure, non_overridable, public :: get_field_generator   => oh_field_generator_info_get_field_generator
        procedure, non_overridable, public :: get_num_components    => oh_field_generator_info_get_num_components
        procedure, non_overridable, public :: free                  => oh_field_generator_info_free
        generic,                    public :: assignment(=)         => oh_field_generator_info_assign
    end type
    
    public :: output_handler_field_generator_t, output_handler_field_generator_info_t
    
contains

!---------------------------------------------------------------------
!< output_handler_field_generator_info_t PROCEDURES
!---------------------------------------------------------------------

    subroutine oh_field_generator_info_free(this)
    !-----------------------------------------------------------------
    !< Free output_handler_field_generator_info_t derived type
    !-----------------------------------------------------------------
        class(output_handler_field_generator_info_t), intent(inout) :: this
    !-----------------------------------------------------------------
        if(allocated(this%name))          deallocate(this%name)
        nullify(this%field_generator)
    end subroutine oh_field_generator_info_free


    subroutine oh_field_generator_info_assign(this, oh_field_generator_info)
    !----------------------------------------------------------------- 
    !< Assign operator overloading for output_handler_field_generator_info_t
    !----------------------------------------------------------------- 
        class(output_handler_field_generator_info_t), intent(inout) :: this
        type(output_handler_field_generator_info_t),  intent(in)    :: oh_field_generator_info
        class(output_handler_field_generator_t), pointer            :: field_generator
        character(len=:), allocatable                   :: name
    !----------------------------------------------------------------- 
        call this%free()
        field_generator => oh_field_generator_info%get_field_generator() 
        name        =  oh_field_generator_info%get_name()
        call this%create(name, field_generator)
    end subroutine oh_field_generator_info_assign


    subroutine oh_field_generator_info_create(this, name, field_generator)
    !-----------------------------------------------------------------
    !< Associate a fe_function with a field name
    !-----------------------------------------------------------------
        class(output_handler_field_generator_info_t)    ,  intent(inout) :: this
        character(len=*)                                ,  intent(in)    :: name
        class(output_handler_field_generator_t), target,   intent(in)    :: field_generator
        
    !-----------------------------------------------------------------
        call this%free()
        this%field_generator => field_generator
        this%name          = name
    end subroutine oh_field_generator_info_create


    function oh_field_generator_info_get_name(this) result(name)
    !-----------------------------------------------------------------
    !< Return the name of a field associated with a fe_function
    !-----------------------------------------------------------------
        class(output_handler_field_generator_info_t), intent(in) :: this
        character(len=:), allocatable                :: name
    !-----------------------------------------------------------------
        assert(allocated(this%name))
        name = this%name
    end function oh_field_generator_info_get_name

    
    function oh_field_generator_info_get_field_type(this) result(field_type)
    !-----------------------------------------------------------------
    !< Return the field_type of a field associated with a fe_function
    !-----------------------------------------------------------------
        class(output_handler_field_generator_info_t), intent(in) :: this
        character(len=:), allocatable                :: field_type
    !-----------------------------------------------------------------
        field_type = this%field_generator%get_field_type()
    end function oh_field_generator_info_get_field_type


    function oh_field_generator_info_get_field_generator(this) result(field_generator)
    !-----------------------------------------------------------------
    !< Return a fe_function pointer
    !-----------------------------------------------------------------
        class(output_handler_field_generator_info_t), target, intent(in) :: this
        class(output_handler_field_generator_t), pointer                 :: field_generator
    !-----------------------------------------------------------------
        assert(associated(this%field_generator))
        field_generator => this%field_generator
    end function oh_field_generator_info_get_field_generator


    function oh_field_generator_info_get_num_components(this) result(num_components)
    !-----------------------------------------------------------------
    !< Return the number of components of the field
    !-----------------------------------------------------------------
        class(output_handler_field_generator_info_t), intent(in) :: this
        integer(ip)                                  :: num_components
    !-----------------------------------------------------------------
 
        select case(this%get_field_type())
            case ( field_type_scalar )
                num_components = 1
            case ( field_type_vector )
                num_components = SPACE_DIM
            case ( field_type_tensor )
                num_components = SPACE_DIM*SPACE_DIM
        end select
        
    end function oh_field_generator_info_get_num_components

end module output_handler_field_generator_names
