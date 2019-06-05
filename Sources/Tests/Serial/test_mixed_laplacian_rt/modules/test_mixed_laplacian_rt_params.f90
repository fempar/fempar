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
module mixed_laplacian_rt_params_names
  use fempar_names
# include "debug.i90"

  implicit none
  private

  character(len=*), parameter :: reference_fe_geo_order_key    = 'reference_fe_geo_order'
  character(len=*), parameter :: reference_fe_order_key        = 'reference_fe_order' 

  type :: mixed_laplacian_rt_params_t      
   contains
     procedure, non_overridable             :: process_parameters
     procedure, non_overridable             :: get_parameter_list
     procedure, non_overridable             :: get_dir_path
     procedure, non_overridable             :: get_prefix
     procedure, non_overridable             :: get_dir_path_out
     procedure, non_overridable             :: get_reference_fe_geo_order
     procedure, non_overridable             :: get_reference_fe_order
  end type mixed_laplacian_rt_params_t

  ! Types
  public :: mixed_laplacian_rt_params_t

contains
  
  !==================================================================================================
  subroutine mixed_laplacian_rt_define_parameters()
    implicit none

    ! IO parameters
    call parameter_handler%add(reference_fe_geo_order_key, '--reference-fe-geo-order', 1, 'Order of the triangulation reference fe', switch_ab='-gorder')
    call parameter_handler%add(reference_fe_order_key, '--reference-fe-order', 1, 'Order of the fe space reference fe', switch_ab='-order')     
  end subroutine mixed_laplacian_rt_define_parameters

  !==================================================================================================

  subroutine process_parameters(this)
    implicit none
    class(mixed_laplacian_rt_params_t) , intent(in)  :: this
    call parameter_handler%process_parameters(mixed_laplacian_rt_define_parameters)
  end subroutine process_parameters

  !==================================================================================================

  function get_parameter_list(this)
    implicit none
    class(mixed_laplacian_rt_params_t) , intent(in) :: this
    type(ParameterList_t), pointer                  :: get_parameter_list
    get_parameter_list  => parameter_handler%get_values()
  end function get_parameter_list
  
  ! GETTERS *****************************************************************************************
  !==================================================================================================
  function get_dir_path(this)
    implicit none
    class(mixed_laplacian_rt_params_t) , intent(in) :: this
    character(len=:), allocatable                   :: get_dir_path
    get_dir_path = parameter_handler%get_dir_path()
  end function get_dir_path

  !==================================================================================================
  function get_prefix(this)
    implicit none
    class(mixed_laplacian_rt_params_t) , intent(in) :: this
    character(len=:), allocatable                   :: get_prefix
    get_prefix = parameter_handler%get_prefix()
  end function get_prefix

  !==================================================================================================
  function get_dir_path_out(this)
    implicit none
    class(mixed_laplacian_rt_params_t) , intent(in) :: this
    character(len=:), allocatable                   :: get_dir_path_out
    get_dir_path_out = parameter_handler%get_dir_path_out()
  end function get_dir_path_out

  !==================================================================================================
  function get_reference_fe_geo_order(this)
    implicit none
    class(mixed_laplacian_rt_params_t) , intent(in) :: this
    integer(ip)                                     :: get_reference_fe_geo_order
    call parameter_handler%Get(key = reference_fe_geo_order_key, Value = get_reference_fe_geo_order)
  end function get_reference_fe_geo_order
  
  !==================================================================================================
  function get_reference_fe_order(this)
    implicit none
    class(mixed_laplacian_rt_params_t) , intent(in) :: this
    integer(ip)                                     :: get_reference_fe_order
    call parameter_handler%Get(key = reference_fe_order_key, Value = get_reference_fe_order)
  end function get_reference_fe_order
  
end module mixed_laplacian_rt_params_names
