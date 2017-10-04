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
module hts_discrete_integration_names
  use fempar_names
  use hts_analytical_functions_names
		use hts_theta_method_names 
  
  implicit none
# include "debug.i90"
  private
		
		integer(ip), parameter :: hts = 1
  integer(ip), parameter :: air = 2
		
		
  type, extends(discrete_integration_t) :: hts_cG_discrete_integration_t
  type(hts_analytical_functions_t), pointer :: analytical_functions => NULL()
	 type(fe_function_t)             , pointer :: H_current          => NULL()
		type(fe_function_t)             , pointer :: H_previous         => NULL()
		type(theta_method_t)            , pointer :: theta_method       => NULL() 
   contains
     procedure :: set_analytical_functions
	    procedure :: set_fe_functions 
					procedure :: set_theta_method 
     procedure :: integrate_galerkin 
  end type hts_cG_discrete_integration_t
  
  public :: hts_cG_discrete_integration_t
  
contains
   
  subroutine set_analytical_functions ( this, analytical_functions )
     implicit none
     class(hts_cG_discrete_integration_t)    ,intent(inout)  :: this
     type(hts_analytical_functions_t), target, intent(in)    :: analytical_functions
     this%analytical_functions => analytical_functions
  end subroutine set_analytical_functions
  
  subroutine set_fe_functions (this, H_previous, H_current)
    implicit none
    class(hts_cG_discrete_integration_t), intent(inout)       :: this
    type(fe_function_t)  , target, intent(in)  :: H_previous
				type(fe_function_t)  , target, intent(in)  :: H_current 
				
    this%H_previous => H_previous 
				this%H_current  => H_current 
  end subroutine set_fe_functions

		  subroutine set_theta_method (this, theta_method)
    implicit none
    class(hts_cG_discrete_integration_t), intent(inout)    :: this
				type(theta_method_t) , target       , intent(in)       :: theta_method
    this%theta_method => theta_method 
  end subroutine set_theta_method
		
  subroutine integrate_galerkin ( this, fe_space, assembler )
    implicit none
    class(hts_cG_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)                , intent(inout) :: fe_space
    class(assembler_t)                      , intent(inout) :: assembler

	! FE space traversal-related data types
    class(fe_cell_iterator_t), allocatable :: fe

    ! FE integration-related data types
    type(quadrature_t)       , pointer :: quad
    type(point_t)            , pointer :: quad_coords(:)
    type(vector_field_t), allocatable  :: shape_curls(:,:)
    type(vector_field_t), allocatable  :: shape_values(:,:)
				type(vector_field_t)               :: H_value_previous
				type(fe_cell_function_vector_t)    :: fe_cell_function_previous

    ! FE matrix and vector i.e., A_K + f_K
    real(rp), allocatable              :: elmat(:,:), elvec(:)

    integer(ip)  :: istat
    integer(ip)  :: qpoint, num_quad_points
    integer(ip)  :: idof, jdof, num_dofs, max_num_dofs 
    real(rp)     :: factor, time_factor 
    type(vector_field_t), allocatable :: source_term_values(:,:)
				real(rp)     :: current_time(1)

    integer(ip)  :: number_fields
    class(vector_function_t), pointer :: source_term
    
    assert (associated(this%analytical_functions)) 
	   assert (associated(this%H_previous))
				assert (associated(this%H_current))

    source_term => this%analytical_functions%get_source_term()
	
	   call fe_space%set_up_cell_integration()
    call fe_space%create_fe_cell_iterator(fe)
				call fe_cell_function_previous%create(fe_space, 1)
    
    max_num_dofs = fe_space%get_max_num_dofs_on_a_cell()
	   call memalloc ( max_num_dofs, max_num_dofs, elmat, __FILE__, __LINE__ )
    call memalloc ( max_num_dofs, elvec, __FILE__, __LINE__ )
				current_time(1) = this%theta_method%get_current_time()
	
	quad            => fe%get_quadrature()
	num_quad_points = quad%get_num_quadrature_points()
	allocate (source_term_values(num_quad_points,1), stat=istat); check(istat==0)
	time_factor = this%theta_method%get_theta() * this%theta_method%get_time_step() 

    do while ( .not. fe%has_finished())
       if ( fe%is_local() ) then
	   
          ! Update FE-integration related data structures
          call fe%update_integration()
          call fe_cell_function_previous%update(fe, this%H_previous)

          quad            => fe%get_quadrature()
          num_quad_points =  quad%get_num_quadrature_points()
          num_dofs        =  fe%get_num_dofs()
		          
          ! Get quadrature coordinates to evaluate source_term
										quad_coords => fe%get_quadrature_points_coordinates()
							
										! Evaluate pressure source term at quadrature points
          call source_term%get_values_set( quad_coords, current_time, source_term_values)
										                    
          ! Compute element matrix and vector
          elmat = 0.0_rp
          elvec = 0.0_rp
          call fe%get_curls(shape_curls)
          call fe%get_values(shape_values)
          do qpoint = 1, num_quad_points
             factor = fe%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
													
										! Previous solution to integrate RHS contribution 
          call fe_cell_function_previous%get_value( qpoint, H_value_previous )
										
             do idof = 1, num_dofs
                do jdof = 1, num_dofs
                   ! A_K(i,j) =  (curl(phi_i),curl(phi_j)) + (phi_i,phi_j)
                   elmat(idof,jdof) = elmat(idof,jdof) + factor * ( 1.0_rp/time_factor * shape_values(jdof,qpoint)*shape_values(idof,qpoint) + & 
																																																																			                      shape_curls(jdof,qpoint)*shape_curls(idof,qpoint))
                end do
             end do
             
             do idof = 1, num_dofs
			       ! F_K(i) = (f(i),phi_i)
                elvec(idof) = elvec(idof) + factor * (source_term_values(qpoint,1) + 1.0_rp/time_factor*H_value_previous) * shape_values(idof,qpoint)
             end do
          end do
          
          ! Apply boundary conditions
		  call fe%assembly( this%H_current, elmat, elvec, assembler )
       end if
       call fe%next()
    end do
	call fe_space%free_fe_cell_iterator(fe)

    deallocate (shape_values, stat=istat);       check(istat==0)
    deallocate (shape_curls, stat=istat);        check(istat==0)
	   deallocate (source_term_values, stat=istat); check(istat==0)
    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )
				call fe_cell_function_previous%free() 
  end subroutine integrate_galerkin
  
end module hts_discrete_integration_names
