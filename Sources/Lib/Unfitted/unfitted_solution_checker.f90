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

!***************************************************************************************************
module unfitted_solution_checker_names
  use fempar_names
  use unfitted_fe_spaces_names
  implicit none
# include "debug.i90"
  private

  ! TODO this is a temporal class to compute error norms
  type :: unfitted_solution_checker_t
    private
    class(serial_fe_space_t), pointer :: fe_space       => null()
    class(fe_function_t),     pointer :: fe_function    => null()
    class(scalar_function_t), pointer :: exact_solution => null()
  contains
    procedure, non_overridable :: create              => unfitted_solution_checker_create
    procedure, non_overridable :: free                => unfitted_solution_checker_free
    procedure, non_overridable :: compute_error_norms => unfitted_solution_checker_compute_error_norms
  end type unfitted_solution_checker_t

  public :: unfitted_solution_checker_t

contains

  !=================================================================================================
  subroutine unfitted_solution_checker_create(this,fe_space,fe_function,exact_solution)
    implicit none
    class(unfitted_solution_checker_t),        intent(inout) :: this
    class(serial_fe_space_t), target, intent(in) :: fe_space
    class(fe_function_t),     target, intent(in) :: fe_function
    class(scalar_function_t), target, intent(in) :: exact_solution
    this%fe_space       => fe_space
    this%fe_function    => fe_function
    this%exact_solution => exact_solution
  end subroutine unfitted_solution_checker_create

  !=================================================================================================
  subroutine unfitted_solution_checker_free(this)
    implicit none
    class(unfitted_solution_checker_t),         intent(inout) :: this
    this%fe_space       => null()
    this%fe_function    => null()
    this%exact_solution => null()
  end subroutine unfitted_solution_checker_free

  !=================================================================================================
  subroutine unfitted_solution_checker_compute_error_norms(this, error_h1_semi_norm, error_l2_norm, h1_semi_norm, l2_norm)
    implicit none
    class(unfitted_solution_checker_t), intent(in) :: this
    real(rp), intent(inout) :: error_h1_semi_norm
    real(rp), intent(inout) :: error_l2_norm
    real(rp), intent(inout) :: h1_semi_norm
    real(rp), intent(inout) :: l2_norm

    type(unfitted_fe_iterator_t) :: fe_iterator
    type(unfitted_fe_accessor_t), target :: fe
    class(fe_accessor_t), pointer :: fe_ptr
    real(rp), allocatable :: nodal_vals(:)
    real(rp), allocatable :: element_vals(:)
    type(vector_field_t), allocatable :: grad_element_vals(:)
    type(fe_map_t)           , pointer :: fe_map
    type(quadrature_t)       , pointer :: quad
    type(point_t)            , pointer :: quad_coords(:)
    type(volume_integrator_t), pointer :: vol_int
    integer(ip) :: num_elem_nodes
    integer(ip) :: num_quad_points
    integer(ip) :: num_dime
    integer(ip) :: istat, igp, idime
    real(rp) :: dV
    real(rp) :: u_exact_gp
    type(vector_field_t) :: grad_u_exact_gp
    real(rp) :: error_l2_gp, error_h1sn_gp
    real(rp) :: l2_gp, h1sn_gp
    class(environment_t), pointer :: environment
    class(serial_fe_space_t), pointer :: fe_space
    class(par_unfitted_fe_space_t),    pointer :: par_unf_fe_space
    class(serial_unfitted_fe_space_t), pointer :: unf_fe_space

    ! Skip if it is not l1 task
    environment => this%fe_space%get_environment()
    if ( .not. environment%am_i_l1_task() ) return

    ! Initialize
    error_h1_semi_norm = 0.0
    error_l2_norm = 0.0
    h1_semi_norm = 0.0
    l2_norm = 0.0

    num_dime = this%fe_space%get_num_dimensions()
    num_elem_nodes =  this%fe_space%get_max_number_shape_functions()
    call memalloc ( num_elem_nodes, nodal_vals, __FILE__, __LINE__ )

    fe_space => this%fe_space
    select type(fe_space)
      class is(serial_unfitted_fe_space_t)
        fe_iterator = fe_space%create_unfitted_fe_iterator()
      class is(par_unfitted_fe_space_t)
        fe_iterator = fe_space%create_unfitted_fe_iterator()
      class default
        check(.false.)
    end select

    do while ( .not. fe_iterator%has_finished() )

       ! Get current FE
       call fe_iterator%current(fe)
       
       ! Skip ghost cells
       if (fe%is_ghost()) then
         call fe_iterator%next(); cycle
       end if

       ! Update FE-integration related data structures
       call fe%update_integration()

       !This cannot be outside the loop
       quad            => fe%get_quadrature()
       num_quad_points = quad%get_number_quadrature_points()
       fe_map          => fe%get_fe_map()
       vol_int         => fe%get_volume_integrator(1)

       ! Recover nodal values
       !TODO we assume a single field
       fe_ptr => fe
       call this%fe_function%gather_nodal_values(fe_ptr,1,nodal_vals)

       !TODO move outside
       call memalloc(num_quad_points,element_vals,__FILE__, __LINE__)
       allocate(grad_element_vals(1:num_quad_points), stat=istat); check(istat==0)

       ! Recover values of the FE solution at integration points
       call vol_int%evaluate_fe_function(nodal_vals,element_vals)
       call vol_int%evaluate_gradient_fe_function(nodal_vals,grad_element_vals)

       ! Physical coordinates of the quadrature points
       quad_coords => fe_map%get_quadrature_points_coordinates()

       ! Loop in quadrature points
       do igp = 1, num_quad_points

         ! Evaluate exact solution at quadrature points
         call this%exact_solution%get_value(quad_coords(igp),u_exact_gp)
         call this%exact_solution%get_gradient(quad_coords(igp),grad_u_exact_gp)

         ! Integrand at gp
         l2_gp = u_exact_gp**2
         error_l2_gp   = (u_exact_gp - element_vals(igp))**2
         h1sn_gp = 0.0
         error_h1sn_gp = 0.0
         ! TODO a way of avoiding this loop? Maybe using TBP of vector_field_t??
         do idime = 1, num_dime
           h1sn_gp = h1sn_gp + ( grad_u_exact_gp%get(idime) )**2
           error_h1sn_gp = error_h1sn_gp + ( grad_u_exact_gp%get(idime) - grad_element_vals(igp)%get(idime) )**2
         end do

         ! Integrate
         dV = fe_map%get_det_jacobian(igp) * quad%get_weight(igp)
         error_l2_norm      = error_l2_norm      + error_l2_gp  *dV
         error_h1_semi_norm = error_h1_semi_norm + error_h1sn_gp*dV
         h1_semi_norm       = h1_semi_norm       + l2_gp        *dV
         l2_norm            = l2_norm            + h1sn_gp      *dV

       end do

       !TODO move outside
       call memfree(element_vals,__FILE__, __LINE__)
       deallocate(grad_element_vals,stat=istat); check(istat==0)

       call fe_iterator%next()
    end do

    call memfree ( nodal_vals, __FILE__, __LINE__ )

    ! All reduce with a sum between l1 tasks
    call environment%l1_sum(error_l2_norm     )
    call environment%l1_sum(error_h1_semi_norm)
    call environment%l1_sum(h1_semi_norm      )
    call environment%l1_sum(l2_norm           )

    ! Do not forget to do the square root
    error_l2_norm      = sqrt(error_l2_norm)
    error_h1_semi_norm = sqrt(error_h1_semi_norm)
    l2_norm            = sqrt(l2_norm)
    h1_semi_norm       = sqrt(h1_semi_norm)

  end subroutine unfitted_solution_checker_compute_error_norms

end module unfitted_solution_checker_names
!***************************************************************************************************
