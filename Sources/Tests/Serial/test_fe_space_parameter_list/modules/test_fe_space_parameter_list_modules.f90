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

!****************************************************************************************************

module poisson_analytical_functions_names
  use fempar_names
  implicit none
# include "debug.i90"
  private
  type, extends(scalar_function_t) :: source_term_t
    private 
   contains
     procedure :: get_value_space    => source_term_get_value_space
  end type source_term_t

  type, extends(scalar_function_t) :: boundary_function_t
    private
   contains
     procedure :: get_value_space => boundary_function_get_value_space
  end type boundary_function_t

  type, extends(scalar_function_t) :: solution_function_t
    private 
   contains
     procedure :: get_value_space    => solution_function_get_value_space
     procedure :: get_gradient_space => solution_function_get_gradient_space
  end type solution_function_t

  public :: source_term_t, boundary_function_t, solution_function_t
contains  

  !===============================================================================================
  subroutine source_term_get_value_space ( this, point, result )
    implicit none
    class(source_term_t), intent(in)    :: this
    type(point_t)       , intent(in)    :: point
    real(rp)            , intent(inout) :: result
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    result = 0.0_rp 
  end subroutine source_term_get_value_space

  !===============================================================================================
  subroutine boundary_function_get_value_space ( this, point, result )
    implicit none
    class(boundary_function_t), intent(in)  :: this
    type(point_t)           , intent(in)    :: point
    real(rp)                , intent(inout) :: result
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    if ( this%get_num_dims() == 2 ) then
      result = point%get(1)+ point%get(2) ! x+y
    else if ( this%get_num_dims() == 3 ) then
      result = point%get(1)+ point%get(2) + point%get(3) ! x+y+z
    end if  
  end subroutine boundary_function_get_value_space 

  !===============================================================================================
  subroutine solution_function_get_value_space ( this, point, result )
    implicit none
    class(solution_function_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    real(rp)                  , intent(inout) :: result
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    if ( this%get_num_dims() == 2 ) then
      result = point%get(1)+ point%get(2) ! x+y 
    else if ( this%get_num_dims() == 3 ) then
      result = point%get(1)+ point%get(2) + point%get(3) ! x+y+z
    end if  
      
  end subroutine solution_function_get_value_space
  
  !===============================================================================================
  subroutine solution_function_get_gradient_space ( this, point, result )
    implicit none
    class(solution_function_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    type(vector_field_t)      , intent(inout) :: result
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    if ( this%get_num_dims() == 2 ) then
      call result%set( 1, 1.0_rp ) 
      call result%set( 2, 1.0_rp )
    else if ( this%get_num_dims() == 3 ) then
      call result%set( 1, 1.0_rp ) 
      call result%set( 2, 1.0_rp )
      call result%set( 3, 1.0_rp ) 
    end if
  end subroutine solution_function_get_gradient_space

end module poisson_analytical_functions_names
!***************************************************************************************************

module poisson_discrete_integration_names
  use fempar_names
  use poisson_analytical_functions_names
  
  implicit none
# include "debug.i90"
  private
  type, extends(discrete_integration_t) :: poisson_discrete_integration_t
     type(source_term_t), pointer :: source_term 
     type(fe_function_t)                 , pointer :: fe_function          => NULL()
   contains
     procedure :: set_source_term
     procedure :: set_fe_function
     procedure :: integrate_galerkin
  end type poisson_discrete_integration_t
  
  public :: poisson_discrete_integration_t
  
contains
  subroutine set_source_term ( this, source_term )
     implicit none
     class(poisson_discrete_integration_t)        , intent(inout) :: this
     type(source_term_t)    , target, intent(in)    :: source_term
     this%source_term => source_term
  end subroutine set_source_term

  subroutine set_fe_function (this, fe_function)
     implicit none
     class(poisson_discrete_integration_t), intent(inout) :: this
     type(fe_function_t)             , target, intent(in)    :: fe_function
     this%fe_function => fe_function
  end subroutine set_fe_function

  subroutine integrate_galerkin ( this, fe_space, assembler )
    implicit none
    class(poisson_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)         , intent(inout) :: fe_space
    class(assembler_t)      , intent(inout) :: assembler

    ! FE space traversal-related data types
    class(fe_cell_iterator_t), allocatable :: fe

    ! FE integration-related data types
    type(quadrature_t)       , pointer :: quad
    type(point_t)            , pointer :: quad_coords(:)
    type(vector_field_t), allocatable  :: shape_gradients(:,:)
    real(rp)            , allocatable  :: shape_values(:,:)

    ! FE matrix and vector i.e., A_K + f_K
    real(rp), allocatable              :: elmat(:,:), elvec(:)

    integer(ip)  :: istat
    integer(ip)  :: qpoint, num_quad_points
    integer(ip)  :: idof, jdof, num_dofs, max_num_dofs
    real(rp)     :: factor
    real(rp)     :: source_term_value
    
    call fe_space%create_fe_cell_iterator(fe)
    max_num_dofs = fe_space%get_max_num_dofs_on_a_cell()
    call memalloc ( max_num_dofs, max_num_dofs, elmat, __FILE__, __LINE__ )
    call memalloc ( max_num_dofs, elvec, __FILE__, __LINE__ )

    do while ( .not. fe%has_finished() )

       if ( fe%is_local() ) then
          ! Update FE-integration related data structures
          call fe%update_integration()

          ! Very important: this has to be inside the loop, as different FEs can be present!
          quad            => fe%get_quadrature()
          num_quad_points = quad%get_num_quadrature_points()
          num_dofs = fe%get_num_dofs()
          
          ! Get quadrature coordinates to evaluate source_term
          quad_coords => fe%get_quadrature_points_coordinates()
                    
          ! Compute element matrix and vector
          elmat = 0.0_rp
          elvec = 0.0_rp
          call fe%get_gradients(shape_gradients)
          call fe%get_values(shape_values)
          do qpoint = 1, num_quad_points
             factor = fe%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
             do idof = 1, num_dofs
                do jdof = 1, num_dofs
                   ! A_K(i,j) = (grad(phi_i),grad(phi_j))
                   elmat(idof,jdof) = elmat(idof,jdof) + factor * shape_gradients(jdof,qpoint) * shape_gradients(idof,qpoint)
                end do
             end do
             
             ! Source term
             call this%source_term%get_value_space(quad_coords(qpoint),source_term_value)
             do idof = 1, num_dofs
                elvec(idof) = elvec(idof) + factor * source_term_value * shape_values(idof,qpoint) 
             end do
          end do
          
          call fe%assembly( this%fe_function, elmat, elvec, assembler )
       end if
       call fe%next()
    end do
    call fe_space%free_fe_cell_iterator(fe)
    call memfree(shape_values, __FILE__, __LINE__)
    deallocate (shape_gradients, stat=istat); check(istat==0);
    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )
  end subroutine integrate_galerkin
  
end module poisson_discrete_integration_names
