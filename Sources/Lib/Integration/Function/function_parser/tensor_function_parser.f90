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
module tensor_function_parser_names
    use types_names
    use field_names
    use function_names
    use scalar_function_parser_names

    implicit none
# include "debug.i90"

    private

    type, extends(tensor_function_t) :: tensor_function_parser_t
    private
        type(p_scalar_function_parser_t)  :: components(SPACE_DIM,SPACE_DIM)
        integer                           :: ncomponents = -1
    contains
        procedure :: create_2D            => tensor_function_parser_create_2D
        procedure :: create_3D            => tensor_function_parser_create_3D
        procedure :: create_isotropic     => tensor_function_parser_create_isotropic
        procedure :: get_value_space      => tensor_function_parser_get_value_space
        procedure :: get_value_space_time => tensor_function_parser_get_value_space_time
        procedure :: free                 => tensor_function_parser_free
        generic   :: create               => create_2D, create_3D, create_isotropic
    end type tensor_function_parser_t

    public :: tensor_function_parser_t
        
contains

    subroutine tensor_function_parser_create_2D( this, component_11, component_12, component_21, component_22)
    !-----------------------------------------------------------------
    !< Initialize a time independant tensor analytical function
    !-----------------------------------------------------------------
        class(tensor_function_parser_t),        intent(inout) :: this
        type(scalar_function_parser_t), target, intent(in)    :: component_11
        type(scalar_function_parser_t), target, intent(in)    :: component_12
        type(scalar_function_parser_t), target, intent(in)    :: component_21
        type(scalar_function_parser_t), target, intent(in)    :: component_22
    !----------------------------------------------------------------- 
        assert( component_11%get_num_dims()==2 .and. component_12%get_num_dims()==2 .and. component_21%get_num_dims()==2 .and. component_22%get_num_dims()==2)
        call this%free()
        this%ncomponents = 4
        call this%set_num_dims(2)
        this%components(1,1)%function => component_11
        this%components(1,2)%function => component_12
        this%components(2,1)%function => component_21
        this%components(2,2)%function => component_22
    end subroutine tensor_function_parser_create_2D


    subroutine tensor_function_parser_create_3D( this, component_11, component_12, component_13, &
                                                       component_21, component_22, component_23, &
                                                       component_31, component_32, component_33)
    !-----------------------------------------------------------------
    !< Initialize a time independant tensor analytical function
    !-----------------------------------------------------------------
        class(tensor_function_parser_t),        intent(inout) :: this
        type(scalar_function_parser_t), target, intent(in)    :: component_11
        type(scalar_function_parser_t), target, intent(in)    :: component_12
        type(scalar_function_parser_t), target, intent(in)    :: component_13
        type(scalar_function_parser_t), target, intent(in)    :: component_21
        type(scalar_function_parser_t), target, intent(in)    :: component_22
        type(scalar_function_parser_t), target, intent(in)    :: component_23
        type(scalar_function_parser_t), target, intent(in)    :: component_31
        type(scalar_function_parser_t), target, intent(in)    :: component_32
        type(scalar_function_parser_t), target, intent(in)    :: component_33
    !----------------------------------------------------------------- 
        assert( component_11%get_num_dims()==3 .and. component_12%get_num_dims()==3 .and. component_13%get_num_dims()==3 )
        assert( component_21%get_num_dims()==3 .and. component_22%get_num_dims()==3 .and. component_23%get_num_dims()==3 )
        assert( component_31%get_num_dims()==3 .and. component_32%get_num_dims()==3 .and. component_33%get_num_dims()==3 )
        call this%free()
        this%ncomponents = 9
        call this%set_num_dims(3)
        this%components(1,1)%function => component_11
        this%components(1,2)%function => component_12
        this%components(1,3)%function => component_13
        this%components(2,1)%function => component_21
        this%components(2,2)%function => component_22
        this%components(2,3)%function => component_23
        this%components(3,1)%function => component_31
        this%components(3,2)%function => component_32
        this%components(3,3)%function => component_33
    end subroutine tensor_function_parser_create_3D


    subroutine tensor_function_parser_create_isotropic( this, component)
    !-----------------------------------------------------------------
    !< Initialize a time independant tensor analytical function
    !-----------------------------------------------------------------
        class(tensor_function_parser_t),        intent(inout) :: this
        type(scalar_function_parser_t), target, intent(in)    :: component
        integer                                               :: num_dims
    !-----------------------------------------------------------------
        num_dims = component%get_num_dims()
        assert(num_dims == 2 .or. num_dims == 3)
        if(num_dims == 2) then
            call this%create(component_11=component, component_12=component, &
                             component_21=component, component_22=component)
        elseif(num_dims == 3) then
            call this%create(component_11=component, component_12=component, component_13=component, &
                             component_21=component, component_22=component, component_23=component, &
                             component_31=component, component_32=component, component_33=component)
        endif 
    end subroutine tensor_function_parser_create_isotropic


    subroutine tensor_function_parser_get_value_space( this, point, result )
    !-----------------------------------------------------------------
    !< Evaluate the time independant tensor analytical function in a given point
    !-----------------------------------------------------------------
        class(tensor_function_parser_t),     intent(in)    :: this
        type(point_t),                       intent(in)    :: point
        type(tensor_field_t),                intent(inout) :: result
        real(rp)                                           :: tmp_result
        integer                                            :: i, j, dims
    !-----------------------------------------------------------------
        dims = this%get_num_dims()
        do j=1, dims
            do i=1, dims
                call this%components(i,j)%function%get_value_space(point, tmp_result)
                call result%set(i, j, tmp_result)
            enddo
        enddo
    end subroutine tensor_function_parser_get_value_space


    subroutine tensor_function_parser_get_value_space_time( this, point, time, result )
    !-----------------------------------------------------------------
    !< Evaluate the time dependant tensor analytical function in a given point and time step
    !-----------------------------------------------------------------
        class(tensor_function_parser_t), intent(in)    :: this
        type(point_t),                   intent(in)    :: point
        real(rp),                        intent(in)    :: time
        type(tensor_field_t),            intent(inout) :: result
        real(rp)                                       :: tmp_result
        integer                                        :: i, j, dims
    !-----------------------------------------------------------------
        dims = this%get_num_dims()
        do j=1, dims
            do i=1, dims
                call this%components(i,j)%function%get_value_space_time(point, time,  tmp_result)
                call result%set(i, j, tmp_result)
            enddo
        enddo
    end subroutine tensor_function_parser_get_value_space_time


    subroutine tensor_function_parser_free( this)
    !-----------------------------------------------------------------
    !< Free tensor analytical function
    !-----------------------------------------------------------------
        class(tensor_function_parser_t), intent(inout) :: this
        integer                                        :: i, j, dims
    !-----------------------------------------------------------------
        dims = this%get_num_dims() 
        do i=1, dims
            do j=1, dims
                nullify(this%components(i,j)%function)
            enddo
        enddo
        this%ncomponents = -1
        call this%set_num_dims(-1) 
    end subroutine tensor_function_parser_free

end module tensor_function_parser_names
