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
module maxwell_nedelec_params_names
  use fempar_names
# include "debug.i90"

  implicit none
  private
  
  character(len=*), parameter :: reference_fe_geo_order_key    = 'reference_fe_geo_order'
  character(len=*), parameter :: reference_fe_order_key        = 'reference_fe_order'
  character(len=*), parameter :: write_solution_key            = 'write_solution'
  character(len=*), parameter :: analytical_function_case_key  = 'analytical_function_case'
     
  type, extends(parameter_handler_t) :: maxwell_nedelec_params_t 
     private 
     contains
       procedure :: define_parameters  => maxwell_nedelec_params_define_parameters
       procedure, non_overridable             :: get_dir_path
       procedure, non_overridable             :: get_prefix
       procedure, non_overridable             :: get_dir_path_out 
       procedure, non_overridable             :: get_triangulation_type 
       procedure, non_overridable             :: get_reference_fe_geo_order
       procedure, non_overridable             :: get_reference_fe_order
       procedure, non_overridable             :: get_write_solution
       procedure, non_overridable             :: get_analytical_function_case
  end type maxwell_nedelec_params_t

  ! Types
  public :: maxwell_nedelec_params_t

contains

 !==================================================================================================
  subroutine maxwell_nedelec_params_define_parameters(this)
    implicit none
    class(maxwell_nedelec_params_t), intent(inout) :: this
    type(ParameterList_t), pointer :: list, switches, switches_ab, helpers, required
    integer(ip)    :: error
    character(len=:), allocatable            :: msg

    list        => this%get_values()
    switches    => this%get_switches()
    switches_ab => this%get_switches_ab()
    helpers     => this%get_helpers()
    required    => this%get_required()

    error = list%set(key = dir_path_key                      , value = '.') ; check(error==0)
    error = list%set(key = prefix_key                        , value = 'square') ; check(error==0)
    error = list%set(key = dir_path_out_key                  , value = '.') ; check(error==0)
    error = list%set(key = struct_hex_triang_num_dims_key    , value =  2)                   ; check(error==0)
    error = list%set(key = struct_hex_triang_num_cells_dir   , value =  [2,2,2])             ; check(error==0)
    error = list%set(key = struct_hex_triang_is_dir_periodic_key, value =  [0,0,0])             ; check(error==0)
    error = list%set(key = reference_fe_geo_order_key        , value =  1)                   ; check(error==0)
    error = list%set(key = reference_fe_order_key            , value =  1)                   ; check(error==0)
    error = list%set(key = write_solution_key                , value =  .false.)             ; check(error==0)
    error = list%set(key = triang_generate_key        , value =  triangulation_generate_structured ) ; check(error==0)
    error = list%set(key = analytical_function_case_key      , value = 'in_fe_space' ) ; check(error==0)

    ! Only some of them are controlled from cli
    error = switches%set(key = dir_path_key                  , value = '--dir-path')                 ; check(error==0)
    error = switches%set(key = prefix_key                    , value = '--prefix')                   ; check(error==0)
    error = switches%set(key = dir_path_out_key              , value = '--dir-path-out')             ; check(error==0)
    error = switches%set(key = struct_hex_triang_num_dims_key , value = '--dim')                      ; check(error==0)
    error = switches%set(key = struct_hex_triang_num_cells_dir, value = '--num_cells')          ; check(error==0)
    error = switches%set(key = struct_hex_triang_is_dir_periodic_key, value =  '--dir-periodic-key') ; check(error==0)
    error = switches%set(key = reference_fe_geo_order_key    , value = '--reference-fe-geo-order')   ; check(error==0)
    error = switches%set(key = reference_fe_order_key        , value = '--reference-fe-order'    )   ; check(error==0)
    error = switches%set(key = write_solution_key            , value = '--write-solution'        )   ; check(error==0)
    error = switches%set(key = triang_generate_key    , value = '--triangulation-type'    )   ; check(error==0)
    error = switches%set(key = analytical_function_case_key  , value = '--analytical_function_case' ); check(error==0)

    error = switches_ab%set(key = dir_path_key               , value = '-d')        ; check(error==0) 
    error = switches_ab%set(key = prefix_key                 , value = '-p')        ; check(error==0) 
    error = switches_ab%set(key = dir_path_out_key           , value = '-o')        ; check(error==0) 
    error = switches_ab%set(key = struct_hex_triang_num_dims_key , value = '-dim')       ; check(error==0)
    error = switches_ab%set(key = struct_hex_triang_num_cells_dir, value = '-n')         ; check(error==0) 
    error = switches_Ab%set(key = struct_hex_triang_is_dir_periodic_key, value =  '-dpk'); check(error==0)
    error = switches_ab%set(key = reference_fe_geo_order_key  , value = '-gorder')   ; check(error==0)
    error = switches_ab%set(key = reference_fe_order_key      , value = '-order')    ; check(error==0)
    error = switches_ab%set(key = write_solution_key          , value = '-wsolution'); check(error==0)
    error = switches_ab%set(key = triang_generate_key  , value = '-tt')       ; check(error==0)
    error = switches_ab%set(key = analytical_function_case_key, value = '-function-case' ); check(error==0)

    error = helpers%set(key = dir_path_key                   , value = 'Directory of the source files')            ; check(error==0)
    error = helpers%set(key = prefix_key                     , value = 'Name of the GiD files')                    ; check(error==0)
    error = helpers%set(key = dir_path_out_key               , value = 'Output Directory')                         ; check(error==0)
    error = helpers%set(key = struct_hex_triang_num_dims_key , value = 'Number of space dimensions')               ; check(error==0)
    error = helpers%set(key = struct_hex_triang_num_cells_dir, value = 'Number of cells per dir')                  ; check(error==0)
    error = helpers%set(key = struct_hex_triang_is_dir_periodic_key, value =  'Is dir periodic key')              ; check(error==0)
    error = helpers%set(key = reference_fe_geo_order_key     , value = 'Order of the triangulation reference fe')  ; check(error==0)
    error = helpers%set(key = reference_fe_order_key         , value = 'Order of the fe space reference fe')       ; check(error==0)
    error = helpers%set(key = write_solution_key             , value = 'Write solution in VTK format')             ; check(error==0)
    error = helpers%set(key = analytical_function_case_key   , value  = 'Select analytical solution case. Possible values: in_fe_space, fichera_2D, fichera_3D' ); check(error==0)
    
    msg = 'structured (1) or unstructured (0) triangulation?'
    error = helpers%set(key = triang_generate_key     , value = msg)  ; check(error==0)
      
    error = required%set(key = dir_path_key                  , value = .false.) ; check(error==0)
    error = required%set(key = prefix_key                    , value = .false.) ; check(error==0)
    error = required%set(key = dir_path_out_key              , value = .false.) ; check(error==0)
    error = required%set(key = struct_hex_triang_num_dims_key, value = .false.) ; check(error==0)
    error = required%set(key = struct_hex_triang_num_cells_dir, value = .false.) ; check(error==0)
    error = required%set(key = struct_hex_triang_is_dir_periodic_key, value =  .false.); check(error==0)
    error = required%set(key = reference_fe_geo_order_key    , value = .false.) ; check(error==0)
    error = required%set(key = reference_fe_order_key        , value = .false.) ; check(error==0)
    error = required%set(key = write_solution_key            , value = .false.) ; check(error==0)
    error = required%set(key = triang_generate_key    , value = .false.) ; check(error==0)
    error = required%set(key = analytical_function_case_key  , value = .false.) ; check(error==0)

  end subroutine maxwell_nedelec_params_define_parameters

  ! GETTERS *****************************************************************************************
  function get_dir_path(this)
    implicit none
    class(maxwell_nedelec_params_t) , intent(in) :: this
    character(len=:),      allocatable            :: get_dir_path
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(dir_path_key, 'string'))
    error = list%GetAsString(key = dir_path_key, string = get_dir_path)
    assert(error==0)
  end function get_dir_path
  
  !==================================================================================================
  function get_dir_path_out(this)
    implicit none
    class(maxwell_nedelec_params_t) , intent(in) :: this
    character(len=:),      allocatable            :: get_dir_path_out
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(dir_path_out_key, 'string'))
    error = list%GetAsString(key = dir_path_out_key, string = get_dir_path_out)
    assert(error==0)
  end function get_dir_path_out

  !==================================================================================================
  function get_prefix(this)
    implicit none
    class(maxwell_nedelec_params_t) , intent(in) :: this
    character(len=:),      allocatable            :: get_prefix
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(prefix_key, 'string'))
    error = list%GetAsString(key = prefix_key, string = get_prefix)
    assert(error==0)
  end function get_prefix

   !==================================================================================================
  function get_triangulation_type(this)
    implicit none
    class(maxwell_nedelec_params_t)  , intent(in) :: this
    integer(ip)                                   :: get_triangulation_type
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(triang_generate_key, get_triangulation_type))
    error = list%Get(key = triang_generate_key, Value = get_triangulation_type)
    assert(error==0)
  end function get_triangulation_type 
  
    !==================================================================================================
  function get_reference_fe_geo_order(this)
    implicit none
    class(maxwell_nedelec_params_t) , intent(in) :: this
    integer(ip)                                   :: get_reference_fe_geo_order
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(reference_fe_geo_order_key, get_reference_fe_geo_order))
    error = list%Get(key = reference_fe_geo_order_key, Value = get_reference_fe_geo_order)
    assert(error==0)
  end function get_reference_fe_geo_order
  
  !==================================================================================================
  function get_reference_fe_order(this)
    implicit none
    class(maxwell_nedelec_params_t) , intent(in) :: this
    integer(ip)                                   :: get_reference_fe_order
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(reference_fe_order_key, get_reference_fe_order))
    error = list%Get(key = reference_fe_order_key, Value = get_reference_fe_order)
    assert(error==0)
  end function get_reference_fe_order
  
  !==================================================================================================
  function get_write_solution(this)
    implicit none
    class(maxwell_nedelec_params_t) , intent(in) :: this
    logical                                       :: get_write_solution
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    logical                                       :: is_present
    logical                                       :: same_data_type
    integer(ip), allocatable                      :: shape(:)
    list  => this%get_values()
    assert(list%isAssignable(write_solution_key, get_write_solution))
    error = list%Get(key = write_solution_key, Value = get_write_solution)
    assert(error==0)
  end function get_write_solution
  
  !==================================================================================================
  function get_analytical_function_case(this)
    implicit none
    class(maxwell_nedelec_params_t) , intent(in) :: this
    character(len=:), allocatable                            :: get_analytical_function_case
    type(ParameterList_t), pointer                           :: list
    integer(ip)                                              :: error
    character(1) :: dummy_string
    list  => this%get_values()
    assert(list%isAssignable(analytical_function_case_key, dummy_string))
    error = list%GetAsString(key = analytical_function_case_key, string = get_analytical_function_case)
    assert(error==0)
  end function get_analytical_function_case
  
  
end module maxwell_nedelec_params_names
