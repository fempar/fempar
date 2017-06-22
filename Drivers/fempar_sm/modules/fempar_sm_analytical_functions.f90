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
module fempar_sm_analytical_functions_names
  use fempar_names
  use fempar_sm_constitutive_models_names
  implicit none
# include "debug.i90"
  private

  type, extends(scalar_function_t) :: base_scalar_function_t
    private
    integer(ip) :: num_dimensions = -1  
  contains
    procedure :: set_num_dimensions    => base_scalar_function_set_num_dimensions
  end type base_scalar_function_t

  type, extends(vector_function_t) :: base_vector_function_t
    integer(ip) :: num_dimensions = -1  
  contains
    procedure :: set_num_dimensions    => base_vector_function_set_num_dimensions
  end type base_vector_function_t

  !===============================================================================================
  type, extends(scalar_function_t) :: zero_scalar_function_t
  contains
     procedure :: get_value_space  => zero_scalar_function_get_value_space
  end type zero_scalar_function_t

  type, extends(vector_function_t) :: zero_vector_function_t
  contains
     procedure :: get_value_space  => zero_vector_function_get_value_space
  end type zero_vector_function_t

  !===============================================================================================
  type, extends(base_scalar_function_t) :: source_term_p_t
    private 
   contains
     procedure :: get_value_space  => source_term_p_get_value_space
  end type source_term_p_t
  
  type, extends(base_vector_function_t) :: source_term_u_t
   contains
     procedure :: get_value_space => source_term_u_get_value_space
  end type source_term_u_t

  !===============================================================================================
  type, extends(base_scalar_function_t) :: solution_function_p_t
    private 
   contains
     procedure :: get_value_space    => solution_function_p_get_value_space
     procedure :: get_gradient_space => solution_function_p_get_gradient_space
  end type solution_function_p_t

  type, extends(base_vector_function_t) :: solution_function_u_t
   contains
     procedure :: get_value_space    => solution_function_u_get_value_space
     procedure :: get_gradient_space => solution_function_u_get_gradient_space
  end type solution_function_u_t
  !===============================================================================================

  type fempar_sm_analytical_functions_t
     private
     type(source_term_p_t)        :: source_term_p
     type(source_term_u_t)        :: source_term_u
     type(solution_function_p_t)  :: solution_function_p
     type(solution_function_u_t)  :: solution_function_u
     type(zero_vector_function_t) :: zero_u
     type(zero_scalar_function_t) :: zero_p
   contains
     procedure :: set_num_dimensions        => fempar_sm_analytical_functions_set_num_dimensions
     procedure :: get_source_term_p         => fempar_sm_analytical_functions_get_source_term_p
     procedure :: get_source_term_u         => fempar_sm_analytical_functions_get_source_term_u
     procedure :: get_solution_function_p   => fempar_sm_analytical_functions_get_solution_function_p
     procedure :: get_solution_function_u   => fempar_sm_analytical_functions_get_solution_function_u
     procedure :: get_zero_function_u       => fempar_sm_analytical_functions_get_zero_function_u
     procedure :: get_zero_function_p       => fempar_sm_analytical_functions_get_zero_function_p
  end type fempar_sm_analytical_functions_t

  public :: fempar_sm_analytical_functions_t

contains  

  subroutine base_scalar_function_set_num_dimensions ( this, num_dimensions )
    implicit none
    class(base_scalar_function_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dimensions
    this%num_dimensions = num_dimensions
  end subroutine base_scalar_function_set_num_dimensions

  subroutine base_vector_function_set_num_dimensions ( this, num_dimensions )
    implicit none
    class(base_vector_function_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dimensions
    this%num_dimensions = num_dimensions
  end subroutine base_vector_function_set_num_dimensions

 !===============================================================================================
  subroutine zero_scalar_function_get_value_space ( this, point, result )
    implicit none
    class(zero_scalar_function_t), intent(in)    :: this
    type(point_t), intent(in)    :: point
    real(rp)     , intent(inout) :: result
    result = 0.0
  end subroutine zero_scalar_function_get_value_space

  subroutine zero_vector_function_get_value_space( this, point, result )
    implicit none
    class(zero_vector_function_t), intent(in) :: this
    type(point_t)             , intent(in)    :: point
    type(vector_field_t)      , intent(inout) :: result
    call result%init(0.0)
  end subroutine zero_vector_function_get_value_space

  !===============================================================================================
  subroutine source_term_p_get_value_space ( this, point, result )
    implicit none
    class(source_term_p_t), intent(in)    :: this
    type(point_t)       , intent(in)    :: point
    real(rp)            , intent(inout) :: result
    assert ( this%num_dimensions == 2 .or. this%num_dimensions == 3 )
    result = 0.0_rp 
  end subroutine source_term_p_get_value_space

  subroutine source_term_u_get_value_space ( this, point, result )
    implicit none
    class(source_term_u_t), intent(in)    :: this
    type(point_t)       , intent(in)    :: point
    type(vector_field_t), intent(inout) :: result
    if ( this%num_dimensions == 2 ) then
      call result%set(1,0.0_rp)
      call result%set(2,0.0_rp) !2 * ( pi**2 ) * sin ( pi * point%get(1) ) * sin ( pi * point%get(2) )
    else
      call result%set(1,0.0_rp)
      call result%set(2,0.0_rp) !2 * ( pi**2 ) * sin ( pi * point%get(1) ) * sin ( pi * point%get(2) )
      call result%set(3,0.0_rp) !2 * ( pi**2 ) * sin ( pi * point%get(1) ) * sin ( pi * point%get(2) )
    end if  
  end subroutine source_term_u_get_value_space

  !===============================================================================================
  subroutine solution_function_p_get_value_space ( this, point, result )
    implicit none
    class(solution_function_p_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    real(rp)                  , intent(inout) :: result
    assert ( this%num_dimensions == 2 .or. this%num_dimensions == 3 )
    result = (lambda+2*one_third*mu)*2.0_rp ! div u
  end subroutine solution_function_p_get_value_space

  subroutine solution_function_u_get_value_space ( this, point, result )
    implicit none
    class(solution_function_u_t), intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    type(vector_field_t)    , intent(inout) :: result
    if ( this%num_dimensions == 2 ) then
      call result%set(1, point%get(1)+point%get(2) ) 
      call result%set(2, point%get(1)+point%get(2) ) 
    else
      call result%set(1, point%get(1)+point%get(2)+point%get(3) ) 
      call result%set(2, point%get(1)+point%get(2)+point%get(3) ) 
      call result%set(3, point%get(1)+point%get(2)+point%get(3) )
    end if
  end subroutine solution_function_u_get_value_space

  !===============================================================================================
  subroutine solution_function_p_get_gradient_space ( this, point, result )
    implicit none
    class(solution_function_p_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    type(vector_field_t)      , intent(inout) :: result
    assert ( this%num_dimensions == 2 .or. this%num_dimensions == 3 )
    result = 0.0_rp
  end subroutine solution_function_p_get_gradient_space

  subroutine solution_function_u_get_gradient_space ( this, point, result )
    implicit none
    class(solution_function_u_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    type(tensor_field_t)      , intent(inout) :: result
    if ( this%num_dimensions == 2 ) then
      call result%set( 1, 1, 1.0_rp ) 
      call result%set( 2, 1, 1.0_rp )
      call result%set( 1, 2, 1.0_rp ) 
      call result%set( 2, 2, 1.0_rp )
    else
      call result%init(1.0_rp)
    end if
    !if ( this%num_dimensions == 2 ) then
    !  call result%set( 1, 1, 1.0_rp ) 
    !  call result%set( 2, 1, 0.0_rp )
    !  call result%set( 1, 2, 1.0_rp ) 
    !  call result%set( 2, 2, 0.0_rp )
    !else
    !  call result%set( 1, 1, 1.0_rp ) 
    !  call result%set( 2, 1, 0.0_rp )
    !  call result%set( 3, 1, 0.0_rp )
    !  call result%set( 1, 2, 1.0_rp ) 
    !  call result%set( 2, 2, 0.0_rp )
    !  call result%set( 3, 2, 0.0_rp )
    !  call result%set( 1, 3, 1.0_rp ) 
    !  call result%set( 2, 3, 0.0_rp )
    !  call result%set( 3, 3, 0.0_rp )
    !end if
  end subroutine solution_function_u_get_gradient_space
  
  !===============================================================================================
  subroutine fempar_sm_analytical_functions_set_num_dimensions ( this, num_dimensions )
    implicit none
    class(fempar_sm_analytical_functions_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dimensions
    call this%source_term_p%set_num_dimensions(num_dimensions)
    call this%source_term_u%set_num_dimensions(num_dimensions)
    call this%solution_function_p%set_num_dimensions(num_dimensions)
    call this%solution_function_u%set_num_dimensions(num_dimensions)
  end subroutine fempar_sm_analytical_functions_set_num_dimensions 
  
  !===============================================================================================
  function fempar_sm_analytical_functions_get_source_term_p ( this )
    implicit none
    class(fempar_sm_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: fempar_sm_analytical_functions_get_source_term_p
    fempar_sm_analytical_functions_get_source_term_p => this%source_term_p
  end function fempar_sm_analytical_functions_get_source_term_p

  function fempar_sm_analytical_functions_get_source_term_u ( this )
    implicit none
    class(fempar_sm_analytical_functions_t), target, intent(in)    :: this
    class(vector_function_t), pointer :: fempar_sm_analytical_functions_get_source_term_u
    fempar_sm_analytical_functions_get_source_term_u => this%source_term_u
  end function fempar_sm_analytical_functions_get_source_term_u
  
  !===============================================================================================
  function fempar_sm_analytical_functions_get_solution_function_p ( this )
    implicit none
    class(fempar_sm_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: fempar_sm_analytical_functions_get_solution_function_p
    fempar_sm_analytical_functions_get_solution_function_p => this%solution_function_p
  end function fempar_sm_analytical_functions_get_solution_function_p

  function fempar_sm_analytical_functions_get_solution_function_u ( this )
    implicit none
    class(fempar_sm_analytical_functions_t), target, intent(in)    :: this
    class(vector_function_t), pointer :: fempar_sm_analytical_functions_get_solution_function_u
    fempar_sm_analytical_functions_get_solution_function_u => this%solution_function_u
  end function fempar_sm_analytical_functions_get_solution_function_u

  !===============================================================================================

  function fempar_sm_analytical_functions_get_zero_function_u ( this )
    implicit none
    class(fempar_sm_analytical_functions_t), target, intent(in)    :: this
    class(vector_function_t), pointer :: fempar_sm_analytical_functions_get_zero_function_u
    fempar_sm_analytical_functions_get_zero_function_u => this%zero_u
  end function fempar_sm_analytical_functions_get_zero_function_u

  function fempar_sm_analytical_functions_get_zero_function_p ( this )
    implicit none
    class(fempar_sm_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: fempar_sm_analytical_functions_get_zero_function_p
    fempar_sm_analytical_functions_get_zero_function_p => this%zero_p
  end function fempar_sm_analytical_functions_get_zero_function_p

end module fempar_sm_analytical_functions_names
!***************************************************************************************************
