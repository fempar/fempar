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

module par_pb_bddc_maxwell_analytical_functions_names
  use fempar_names
  use par_pb_bddc_maxwell_params_names
  implicit none
# include "debug.i90"
  private
    
  type, extends(scalar_function_t) :: curl_curl_coeff_function_t
   private 
    character(len=:), allocatable :: coefficient_case 
    real(rp)                      :: default_value 
    integer(ip)                   :: num_peaks
   contains
     procedure :: set_value          => curl_curl_coeff_function_set_default_value
     procedure :: set_num_peaks      => curl_curl_coeff_function_set_num_peaks 
     procedure :: set_coefficient_case => curl_curl_coeff_function_set_coefficient_case 
  end type curl_curl_coeff_function_t 
  
  type :: curl_curl_coeff_holder_t 
     class(curl_curl_coeff_function_t), pointer :: p => NULL()  
  end type curl_curl_coeff_holder_t 
  
  type, extends(curl_curl_coeff_function_t) :: curl_curl_coeff_function_white_t
    private 
   contains
     procedure :: get_value_space    => curl_curl_coeff_function_white_get_value_space
  end type curl_curl_coeff_function_white_t 
  
  type, extends(curl_curl_coeff_function_t) :: curl_curl_coeff_function_black_t
    private 
   contains
     procedure :: get_value_space    => curl_curl_coeff_function_black_get_value_space
  end type curl_curl_coeff_function_black_t 
  
  type, extends(scalar_function_t) :: mass_coeff_function_t
   private 
    character(len=:), allocatable :: coefficient_case
    real(rp)                      :: default_value 
    integer(ip)                   :: num_peaks 
   contains
     procedure :: set_value          => mass_coeff_function_set_default_value
     procedure :: set_num_peaks      => mass_coeff_function_set_num_peaks 
     procedure :: set_coefficient_case => mass_coeff_function_set_coefficient_case 
  end type mass_coeff_function_t 
  
  type :: mass_coeff_holder_t 
     class(mass_coeff_function_t), pointer :: p => NULL()  
  end type mass_coeff_holder_t 
  
  type, extends(mass_coeff_function_t) :: mass_coeff_function_white_t
    private 
   contains
     procedure :: get_value_space    => mass_coeff_function_white_get_value_space
  end type mass_coeff_function_white_t 
  
  type, extends(mass_coeff_function_t) :: mass_coeff_function_black_t
    private 
   contains
     procedure :: get_value_space    => mass_coeff_function_black_get_value_space
  end type mass_coeff_function_black_t 
  
  ! Boundary scalar functions definition 
    type, extends(scalar_function_t) :: boundary_function_Hx_t
    private 
   contains
     procedure :: get_value_space    => boundary_function_Hx_get_value_space
  end type boundary_function_Hx_t 
  
    type, extends(scalar_function_t) :: boundary_function_Hy_t
    private 
   contains
     procedure :: get_value_space    => boundary_function_Hy_get_value_space
  end type boundary_function_Hy_t 
  
   type, extends(scalar_function_t) :: boundary_function_Hz_t
    private 
   contains
     procedure :: get_value_space    => boundary_function_Hz_get_value_space
  end type boundary_function_Hz_t 
  
  ! Vector functions 
    type, extends(vector_function_t) :: source_term_t
     real(rp) :: mass_coeff 
     real(rp) :: curl_curl_coeff 
   contains
     procedure :: get_value_space => source_term_get_value_space
  end type source_term_t

   type, extends(vector_function_t) :: solution_t
   contains
     procedure :: get_value_space    => solution_get_value_space
     procedure :: get_gradient_space => solution_get_gradient_space
  end type solution_t
  
  type par_pb_bddc_maxwell_analytical_functions_t
     private
     type(curl_curl_coeff_holder_t), pointer :: curl_curl_coeff(:)
     type(mass_coeff_holder_t), pointer      :: mass_coeff(:) 
     type(boundary_function_Hx_t)            :: boundary_function_Hx
     type(boundary_function_Hy_t)            :: boundary_function_Hy
     type(boundary_function_Hz_t)            :: boundary_function_Hz
     type(source_term_t)                     :: source_term
     type(solution_t)                        :: solution
   contains
     procedure :: set_num_dims                     => mn_set_num_dims
     procedure :: set_curl_curl_coeff              => mn_set_curl_curl_coeff_holder
     procedure :: set_mass_coeff                   => mn_set_mass_coeff_holder 
     ! Getters 
     procedure :: get_curl_curl_coeff              => mn_get_curl_curl_coeff
     procedure :: get_mass_coeff                   => mn_get_mass_coeff 
     procedure :: get_boundary_function_Hx         => mn_get_boundary_function_Hx
     procedure :: get_boundary_function_Hy         => mn_get_boundary_function_Hy
     procedure :: get_boundary_function_Hz         => mn_get_boundary_function_Hz
     procedure :: get_source_term                  => mn_get_source_term
     procedure :: get_solution_function            => mn_get_solution_function
  end type par_pb_bddc_maxwell_analytical_functions_t

  public :: par_pb_bddc_maxwell_analytical_functions_t
  public :: curl_curl_coeff_holder_t, curl_curl_coeff_function_white_t, curl_curl_coeff_function_black_t 
  public :: mass_coeff_holder_t, mass_coeff_function_white_t, mass_coeff_function_black_t 

contains  
  !===============================================================================================
  subroutine curl_curl_coeff_function_set_coefficient_case ( this, coefficient_case )
    implicit none
    class(curl_curl_coeff_function_t), intent(inout)    :: this
    character(len=*), intent(in) ::  coefficient_case
    this%coefficient_case = coefficient_case
  end subroutine curl_curl_coeff_function_set_coefficient_case 
  
    !===============================================================================================
  subroutine mass_coeff_function_set_coefficient_case ( this, coefficient_case )
    implicit none
    class(mass_coeff_function_t), intent(inout)    :: this
    character(len=*), intent(in) ::  coefficient_case
    this%coefficient_case = coefficient_case
  end subroutine mass_coeff_function_set_coefficient_case
  
    !===============================================================================================
  subroutine curl_curl_coeff_function_set_num_peaks ( this, num_peaks )
    implicit none
    class(curl_curl_coeff_function_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_peaks
    this%num_peaks = num_peaks
  end subroutine curl_curl_coeff_function_set_num_peaks 
  
      !===============================================================================================
  subroutine mass_coeff_function_set_num_peaks ( this, num_peaks )
    implicit none
    class(mass_coeff_function_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_peaks
    this%num_peaks = num_peaks
  end subroutine mass_coeff_function_set_num_peaks 
    
  !===============================================================================================
  subroutine curl_curl_coeff_function_set_default_value ( this, default_value )
    implicit none
    class(curl_curl_coeff_function_t), intent(inout)    :: this
    real(rp), intent(in) ::  default_value
    this%default_value = default_value
  end subroutine curl_curl_coeff_function_set_default_value
  
    !===============================================================================================
  subroutine mass_coeff_function_set_default_value ( this, default_value )
    implicit none
    class(mass_coeff_function_t), intent(inout)    :: this
    real(rp), intent(in) ::  default_value
    this%default_value = default_value
  end subroutine mass_coeff_function_set_default_value
    
   !===============================================================================================
  subroutine curl_curl_coeff_function_white_get_value_space( this, point, result )
    implicit none 
    class(curl_curl_coeff_function_white_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(inout) :: result 
  
    select case ( this%coefficient_case )
    case ( constant ) 
      result = this%default_value 
    case ( sinusoidal ) 
      result = 10**( this%default_value*sin(this%num_peaks*pi*(point%get(1))))
    case DEFAULT 
      massert( .false. , 'Invalid curl_curl_coeff coefficient case' ) 
    end select 

  end subroutine curl_curl_coeff_function_white_get_value_space
  
    !===============================================================================================
  subroutine curl_curl_coeff_function_black_get_value_space( this, point, result )
    implicit none 
    class(curl_curl_coeff_function_black_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(inout) :: result 
    result = this%default_value 
    
    select case ( this%coefficient_case ) 
    case ( constant ) 
      result = this%default_value
    case DEFAULT 
      massert( .false., 'Black subdomains only allocate constant functions' ) 
    end select 
  end subroutine curl_curl_coeff_function_black_get_value_space
  
     !===============================================================================================
  subroutine mass_coeff_function_white_get_value_space( this, point, result )
    implicit none 
    class(mass_coeff_function_white_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(inout) :: result 

    select case ( this%coefficient_case )
    case ( constant ) 
      result = this%default_value 
    case ( sinusoidal ) 
      result = 10**( this%default_value*sin(this%num_peaks*pi*(point%get(2))))
    case DEFAULT 
      massert( .false. , 'Invalid curl_curl_coeff coefficient case' ) 
    end select 

  end subroutine mass_coeff_function_white_get_value_space
  
    !===============================================================================================
  subroutine mass_coeff_function_black_get_value_space( this, point, result )
    implicit none 
    class(mass_coeff_function_black_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(inout) :: result 
 
    select case ( this%coefficient_case ) 
    case ( constant ) 
      result = this%default_value
    case DEFAULT 
      massert( .false., 'Black subdomains only allocate constant functions' ) 
    end select 
    
  end subroutine mass_coeff_function_black_get_value_space
  
      !===============================================================================================
  subroutine boundary_function_Hx_get_value_space( this, point, result )
    implicit none 
    class(boundary_function_Hx_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(inout) :: result 
    real(rp) :: x,y,z 
    x = point%get(1); y=point%get(2); z=point%get(3)
    result = -y 
  end subroutine boundary_function_Hx_get_value_space
  
  !===============================================================================================
  subroutine boundary_function_Hy_get_value_space( this, point, result )
    implicit none 
    class(boundary_function_Hy_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(inout) :: result 
    real(rp) :: x,y,z 
    x = point%get(1); y=point%get(2); z=point%get(3)
    result = x
  end subroutine boundary_function_Hy_get_value_space

  !===============================================================================================
  subroutine boundary_function_Hz_get_value_space( this, point, result )
    implicit none 
    class(boundary_function_Hz_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(inout) :: result 
    real(rp) :: x,y,z 
    x = point%get(1); y=point%get(2); z=point%get(3)
    result = 0.0_rp
  end subroutine boundary_function_Hz_get_value_space

 !===============================================================================================
  subroutine source_term_get_value_space ( this, point, result )
    implicit none
    class(source_term_t)    , intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    type(vector_field_t)    , intent(inout) :: result
    real(rp) :: x,y,z 
    x = point%get(1); y=point%get(2); z=point%get(3)
    call result%init(0.0_rp)
    call result%set(1, -y)
    call result%set(2,  x) 
  end subroutine source_term_get_value_space

  !===============================================================================================
  subroutine solution_get_value_space ( this, point, result )
    implicit none
    class(solution_t)       , intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    type(vector_field_t)    , intent(inout) :: result
 
    real(rp) :: x,y,z 
    x = point%get(1); y=point%get(2); z=point%get(3) 
    ! call result%init(0.0_rp) 
    call result%set(1, -y ) 
    call result%set(2,  x ) 
    call result%set(3,  0.0_rp )
  
  end subroutine solution_get_value_space

  !===============================================================================================
  subroutine solution_get_gradient_space ( this, point, result )
    implicit none
    class(solution_t)   , intent(in)    :: this
    type(point_t)       , intent(in)    :: point
    type(tensor_field_t), intent(inout) :: result
 
    real(rp) :: x,y,z 
    x = point%get(1); y=point%get(2); z=point%get(3)
    call result%init(0.0_rp)  
    call result%set(2,1, -1.0_rp)
    call result%set(1,2,  1.0_rp)

  end subroutine solution_get_gradient_space
  
  !===============================================================================================
  subroutine mn_set_num_dims ( this, num_dims )
    implicit none
    class(par_pb_bddc_maxwell_analytical_functions_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dims
    call this%source_term%set_num_dims(num_dims)
    call this%solution%set_num_dims(num_dims)
  end subroutine mn_set_num_dims
  
    !===============================================================================================
  subroutine mn_set_curl_curl_coeff_holder ( this, curl_curl_coeff_holder )
    implicit none
    class(par_pb_bddc_maxwell_analytical_functions_t), intent(inout)    :: this
    type(curl_curl_coeff_holder_t), target   , intent(in)       :: curl_curl_coeff_holder(:)
    this%curl_curl_coeff => curl_curl_coeff_holder 
  end subroutine mn_set_curl_curl_coeff_holder
  
      !===============================================================================================
  subroutine mn_set_mass_coeff_holder ( this, mass_coeff_holder )
    implicit none
    class(par_pb_bddc_maxwell_analytical_functions_t), intent(inout)    :: this
    type(mass_coeff_holder_t), target   , intent(in)      :: mass_coeff_holder(:)
    this%mass_coeff => mass_coeff_holder 
  end subroutine mn_set_mass_coeff_holder
  
      !===============================================================================================
  function mn_get_curl_curl_coeff ( this )
    implicit none
    class(par_pb_bddc_maxwell_analytical_functions_t), target, intent(in)    :: this
    type(curl_curl_coeff_holder_t), pointer     :: mn_get_curl_curl_coeff(:)
    mn_get_curl_curl_coeff => this%curl_curl_coeff 
  end function mn_get_curl_curl_coeff
  
        !===============================================================================================
  function mn_get_mass_coeff ( this )
    implicit none
    class(par_pb_bddc_maxwell_analytical_functions_t), target, intent(in)    :: this
    type(mass_coeff_holder_t), pointer     :: mn_get_mass_coeff(:)
    mn_get_mass_coeff => this%mass_coeff 
  end function mn_get_mass_coeff
  
  !===============================================================================================
  function mn_get_boundary_function_Hx ( this )
    implicit none
    class(par_pb_bddc_maxwell_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: mn_get_boundary_function_Hx
    mn_get_boundary_function_Hx => this%boundary_function_Hx 
  end function mn_get_boundary_function_Hx

  !===============================================================================================
  function mn_get_boundary_function_Hy ( this )
    implicit none
    class(par_pb_bddc_maxwell_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: mn_get_boundary_function_Hy
    mn_get_boundary_function_Hy => this%boundary_function_Hy 
  end function mn_get_boundary_function_Hy

  !===============================================================================================
  function mn_get_boundary_function_Hz ( this )
    implicit none
    class(par_pb_bddc_maxwell_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: mn_get_boundary_function_Hz
    mn_get_boundary_function_Hz => this%boundary_function_Hz 
  end function mn_get_boundary_function_Hz
  
  !===============================================================================================
  function mn_get_solution_function ( this )
    implicit none
    class(par_pb_bddc_maxwell_analytical_functions_t), target, intent(in)    :: this
    class(vector_function_t), pointer :: mn_get_solution_function
    mn_get_solution_function => this%solution
  end function mn_get_solution_function

  !===============================================================================================
  function mn_get_source_term ( this )
    implicit none
    class(par_pb_bddc_maxwell_analytical_functions_t), target, intent(in)    :: this
    class(vector_function_t), pointer :: mn_get_source_term
    mn_get_source_term => this%source_term
  end function mn_get_source_term

end module par_pb_bddc_maxwell_analytical_functions_names
!***************************************************************************************************
