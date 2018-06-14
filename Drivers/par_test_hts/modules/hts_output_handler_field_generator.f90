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
module hts_output_handler_field_generator_names
  use fempar_names

# include "debug.i90"
  implicit none
  private

  type, extends(output_handler_field_generator_t) :: resistivity_field_generator_t
     private
     type(fe_function_t),       pointer :: magnetic_field 
     real(rp),                  pointer :: nodal_values(:)
     integer(ip)                        :: field_id = 1
     ! Parameter values 
     real(rp)                           :: nonlinear_exponent 
     real(rp)                           :: critical_current
     real(rp)                           :: critical_electric_field 
   contains
     procedure :: set_parameter_values      =>  resistivity_field_generator_set_parameter_values
     procedure :: set_magnetic_field        =>  resistivity_field_generator_set_magnetic_field 
     procedure :: set_nodal_values          =>  resistivity_field_generator_set_nodal_values
     procedure :: get_field_type            =>  resistivity_field_generator_get_field_type    
     procedure :: allocate_function_values  =>  resistivity_field_generator_allocate_function_values
     procedure :: compute_function_values   =>  resistivity_field_generator_compute_function_values
     procedure :: generate_patch_field      =>  resistivity_field_generator_generate_patch_field
  end type resistivity_field_generator_t

public :: resistivity_field_generator_t 

contains

  subroutine resistivity_field_generator_set_parameter_values(this,nonlinear_exponent, critical_current, critical_electric_field)
    implicit none
    class(resistivity_field_generator_t) ,  intent(inout) :: this
    real(rp)                             ,  intent(in)    :: nonlinear_exponent 
    real(rp)                             ,  intent(in)    :: critical_current 
    real(rp)                             ,  intent(in)    :: critical_electric_field 
    this%nonlinear_exponent      = nonlinear_exponent 
    this%critical_current        = critical_current
    this%critical_electric_field = critical_electric_field
  end subroutine resistivity_field_generator_set_parameter_values
  
  subroutine resistivity_field_generator_set_magnetic_field(this,magnetic_field)
    implicit none
    class(resistivity_field_generator_t) ,  intent(inout) :: this
    type(fe_function_t),    target,         intent(in)    :: magnetic_field
    this%magnetic_field => magnetic_field
  end subroutine resistivity_field_generator_set_magnetic_field
  
  subroutine resistivity_field_generator_set_nodal_values(this,oh_cell_fe_function,patch_field)
    implicit none
    class(resistivity_field_generator_t)         ,  intent(inout) :: this
    type(output_handler_fe_cell_function_t)      ,  intent(in)    :: oh_cell_fe_function
    type(output_handler_patch_field_t)           ,  intent(inout) :: patch_field
    class(fe_cell_iterator_t)                    ,  pointer       :: current_fe
    class(reference_fe_t)                        ,  pointer       :: reference_fe
    type(std_vector_real_rp_t)                   ,  pointer       :: patch_field_nodal_values

    ! Get reference_Fe
    current_fe => oh_cell_fe_function%get_fe()
    reference_fe => current_fe%get_reference_fe(this%field_id)
    assert(reference_fe%get_field_type() == field_type_vector)

    ! Gather DoFs of current cell + field_id on nodal_values 
    patch_field_nodal_values => patch_field%get_nodal_values()
    call patch_field_nodal_values%resize(reference_fe%get_num_shape_functions())
    this%nodal_values => patch_field_nodal_values%get_pointer()
    call this%magnetic_field%gather_nodal_values(current_fe, this%field_id, this%nodal_values)

  end subroutine resistivity_field_generator_set_nodal_values
     
  function resistivity_field_generator_get_field_type(this) result(field_type)
    implicit none
    class(resistivity_field_generator_t),   intent(in)  :: this
    character(len=:),                       allocatable :: field_type
    field_type = field_type_scalar
  end function resistivity_field_generator_get_field_type

  subroutine resistivity_field_generator_generate_patch_field(this, oh_cell_fe_function, visualization_points_coordinates, patch_field)
    implicit none          
    class(resistivity_field_generator_t)   , intent(inout)  :: this
    type(output_handler_fe_cell_function_t), intent(in)     :: oh_cell_fe_function
    type(point_t)                          , intent(in)     :: visualization_points_coordinates(:)
    type(output_handler_patch_field_t)     , intent(inout)  :: patch_field

    type(allocatable_array_rp1_t),                 pointer :: patch_field_scalar_function_values
    real(rp),                                  allocatable :: scalar_function_values(:)
    integer(ip)                                            :: npoints

    npoints = size(visualization_points_coordinates)
    call this%set_nodal_values(oh_cell_fe_function,patch_field)
    call patch_field%set_field_type(this%get_field_type())
    
    patch_field_scalar_function_values => patch_field%get_scalar_function_values()
    call this%allocate_function_values(npoints,patch_field_scalar_function_values,scalar_function_values)

    call this%compute_function_values(oh_cell_fe_function,visualization_points_coordinates,scalar_function_values)

    call patch_field_scalar_function_values%move_alloc_in(scalar_function_values) 

  end subroutine resistivity_field_generator_generate_patch_field
  

  subroutine resistivity_field_generator_allocate_function_values(this,npoints,scalar_function_values_aa,scalar_function_values)
    implicit none
    class(resistivity_field_generator_t),               intent(inout) :: this
    integer(ip),                                           intent(in) :: npoints
    type(allocatable_array_rp1_t),                   intent(inout)    :: scalar_function_values_aa
    real(rp),                              allocatable, intent(inout) :: scalar_function_values(:)
  
    call scalar_function_values_aa%move_alloc_out(scalar_function_values) 

    if (.NOT. allocated(scalar_function_values)) then
       call memalloc(npoints, scalar_function_values, __FILE__, __LINE__ )
    end if

  end subroutine resistivity_field_generator_allocate_function_values
 
  subroutine resistivity_field_generator_compute_function_values(this,oh_cell_fe_function,visualization_points_coordinates,scalar_function_values)
    implicit none
    class(resistivity_field_generator_t),                           intent(in)     :: this
    type(output_handler_fe_cell_function_t),                        intent(in)     :: oh_cell_fe_function
    type(point_t)                                                 , intent(in)     :: visualization_points_coordinates(:)  
    real(rp),                                 allocatable,          intent(inout)  :: scalar_function_values(:)
    
    type(vector_field_t)                                                           :: curl_values 
    type(tensor_field_t),                     allocatable                          :: grad_function_values(:)
    type(cell_integrator_t),                                              pointer  :: cell_integrator
    integer(ip)                                                                    :: npoints     
    integer(ip)                                                                    :: ipoint

    cell_integrator => oh_cell_fe_function%get_cell_integrator(this%field_id) 
    call cell_integrator%evaluate_gradient_fe_function(this%nodal_values,grad_function_values)
    npoints = size(visualization_points_coordinates)
    do ipoint = 1,npoints
       call curl_values%set(1, grad_function_values(ipoint)%get(2,3) - grad_function_values(ipoint)%get(3,2) ) 
       call curl_values%set(2, grad_function_values(ipoint)%get(3,1) - grad_function_values(ipoint)%get(1,3) ) 
       call curl_values%set(3, grad_function_values(ipoint)%get(1,2) - grad_function_values(ipoint)%get(2,1) ) 
       scalar_function_values(ipoint) = this%critical_electric_field/this%critical_current*           & 
                                       (curl_values%nrm2()/this%critical_current)**this%nonlinear_exponent   
    end do

  end subroutine resistivity_field_generator_compute_function_values

end module hts_output_handler_field_generator_names
