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
  use types_names
  use memor_names
  use field_names
  
  use reference_fe_names
  use fe_space_names
  use function_names
  use environment_names
  
  use unfitted_fe_spaces_names
  use piecewise_cell_map_names
  implicit none
# include "debug.i90"
  private

  ! TODO this is a temporal class to compute error norms
  type :: unfitted_solution_checker_t
    private
    class(serial_fe_space_t), pointer :: fe_space       => null()
    class(fe_function_t),     pointer :: fe_function    => null()
    class(scalar_function_t), pointer :: exact_solution => null()
    integer(ip) :: field_id
  contains
    procedure, non_overridable :: create              => unfitted_solution_checker_create
    procedure, non_overridable :: free                => unfitted_solution_checker_free
    procedure, non_overridable :: compute_error_norms => unfitted_solution_checker_compute_error_norms
  end type unfitted_solution_checker_t

  public :: unfitted_solution_checker_t

contains

  !=================================================================================================
  subroutine unfitted_solution_checker_create(this,fe_space,fe_function,exact_solution,field_id)
    implicit none
    class(unfitted_solution_checker_t),        intent(inout) :: this
    class(serial_fe_space_t), target, intent(in) :: fe_space
    class(fe_function_t),     target, intent(in) :: fe_function
    class(scalar_function_t), target, intent(in) :: exact_solution
    integer(ip)          , optional , intent(in) :: field_id
    this%fe_space       => fe_space
    this%fe_function    => fe_function
    this%exact_solution => exact_solution
    if (present(field_id)) then
      this%field_id = field_id
    else
      this%field_id = 1
    end if
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
  subroutine unfitted_solution_checker_compute_error_norms(this,&
      error_h1_semi_norm, error_l2_norm, h1_semi_norm, l2_norm,&
      error_h1_semi_norm_boundary, error_l2_norm_boundary, h1_semi_norm_boundary, l2_norm_boundary)
    implicit none
    class(unfitted_solution_checker_t), intent(in) :: this
    real(rp), intent(inout) :: error_h1_semi_norm
    real(rp), intent(inout) :: error_l2_norm
    real(rp), intent(inout) :: h1_semi_norm
    real(rp), intent(inout) :: l2_norm

    real(rp), optional , intent(inout) :: error_h1_semi_norm_boundary
    real(rp), optional , intent(inout) :: error_l2_norm_boundary
    real(rp), optional , intent(inout) :: h1_semi_norm_boundary
    real(rp), optional , intent(inout) :: l2_norm_boundary

    class(fe_cell_iterator_t), allocatable :: fe
    real(rp), allocatable :: nodal_vals(:)
    real(rp), allocatable :: element_vals(:)
    type(vector_field_t), allocatable :: grad_element_vals(:)
    type(cell_map_t)           , pointer :: cell_map
    type(piecewise_cell_map_t) , pointer :: pw_cell_map
    
    type(quadrature_t)       , pointer :: quad
    type(point_t)            , pointer :: quad_coords(:)
    type(cell_integrator_t), pointer :: cell_int
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

    ! Skip if it is not l1 task
    environment => this%fe_space%get_environment()
    if ( .not. environment%am_i_l1_task() ) return

    ! Initialize
    error_h1_semi_norm = 0.0
    error_l2_norm = 0.0
    h1_semi_norm = 0.0
    l2_norm = 0.0

    if (present(error_l2_norm_boundary)) then
      error_h1_semi_norm_boundary= 0.0
      error_l2_norm_boundary     = 0.0
      h1_semi_norm_boundary      = 0.0
      l2_norm_boundary           = 0.0
    end if
    
    num_dime = this%fe_space%get_num_dims()
    num_elem_nodes =  this%fe_space%get_max_num_shape_functions()
    call memalloc ( num_elem_nodes, nodal_vals, __FILE__, __LINE__ )


    call this%fe_space%create_fe_cell_iterator(fe)
    do while ( .not. fe%has_finished() )

       ! Skip ghost cells
       if (fe%is_ghost()) then
         call fe%next(); cycle
       end if

       ! Update FE-integration related data structures
       call fe%update_integration()

       !This cannot be outside the loop
       quad            => fe%get_quadrature()
       num_quad_points = quad%get_num_quadrature_points()
       cell_int         => fe%get_cell_integrator(this%field_id)

       ! Recover nodal values
       !TODO we assume a single field
       call this%fe_function%gather_nodal_values(fe,this%field_id,nodal_vals)

       !TODO move outside
       call memalloc(num_quad_points,element_vals,__FILE__, __LINE__)
       allocate(grad_element_vals(1:num_quad_points), stat=istat); check(istat==0)

       ! Recover values of the FE solution at integration points
       call cell_int%evaluate_fe_function(nodal_vals,element_vals)
       call cell_int%evaluate_gradient_fe_function(nodal_vals,grad_element_vals)

       ! Physical coordinates of the quadrature points
       quad_coords => fe%get_quadrature_points_coordinates()

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
         dV = fe%get_det_jacobian(igp) * quad%get_weight(igp)
         error_l2_norm      = error_l2_norm      + error_l2_gp  *dV
         error_h1_semi_norm = error_h1_semi_norm + error_h1sn_gp*dV
         h1_semi_norm       = h1_semi_norm       + h1sn_gp*dV
         l2_norm            = l2_norm            + l2_gp  *dV

       end do

       !TODO move outside
       call memfree(element_vals,__FILE__, __LINE__)
       deallocate(grad_element_vals,stat=istat); check(istat==0)


       if (present(error_l2_norm_boundary)) then
         if (fe%is_cut()) then

           call fe%update_boundary_integration()

           ! Get info on the unfitted boundary for integrating BCs
           quad            => fe%get_boundary_quadrature()
           num_quad_points = quad%get_num_quadrature_points()
           pw_cell_map       => fe%get_boundary_piecewise_cell_map()
           quad_coords     => pw_cell_map%get_quadrature_points_coordinates()
           cell_int         => fe%get_boundary_cell_integrator(this%field_id)

           !TODO move outside
           call memalloc(num_quad_points,element_vals,__FILE__, __LINE__)
           allocate(grad_element_vals(1:num_quad_points), stat=istat); check(istat==0)

           ! Recover values of the FE solution at integration points
           call cell_int%evaluate_fe_function(nodal_vals,element_vals)
           call cell_int%evaluate_gradient_fe_function(nodal_vals,grad_element_vals)

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
             dV = pw_cell_map%get_det_jacobian(igp) * quad%get_weight(igp)
             error_l2_norm_boundary     = error_l2_norm_boundary      + error_l2_gp  *dV
             error_h1_semi_norm_boundary= error_h1_semi_norm_boundary + error_h1sn_gp*dV
             h1_semi_norm_boundary      = h1_semi_norm_boundary       + h1sn_gp*dV
             l2_norm_boundary           = l2_norm_boundary            + l2_gp  *dV

           end do

           !TODO move outside
           call memfree(element_vals,__FILE__, __LINE__)
           deallocate(grad_element_vals,stat=istat); check(istat==0)

         end if
       end if

       call fe%next()
    end do
    call this%fe_space%free_fe_cell_iterator(fe)

    call memfree ( nodal_vals, __FILE__, __LINE__ )

    ! All reduce with a sum between l1 tasks
    call environment%l1_sum(error_l2_norm     )
    call environment%l1_sum(error_h1_semi_norm)
    call environment%l1_sum(h1_semi_norm      )
    call environment%l1_sum(l2_norm           )

    if (present(error_l2_norm_boundary)) then
      call environment%l1_sum(error_l2_norm_boundary     )
      call environment%l1_sum(error_h1_semi_norm_boundary)
      call environment%l1_sum(h1_semi_norm_boundary      )
      call environment%l1_sum(l2_norm_boundary           )
    end if

    ! Do not forget to do the square root
    error_l2_norm      = sqrt(error_l2_norm)
    error_h1_semi_norm = sqrt(error_h1_semi_norm)
    l2_norm            = sqrt(l2_norm)
    h1_semi_norm       = sqrt(h1_semi_norm)

    if (present(error_l2_norm_boundary)) then
      error_l2_norm_boundary     = sqrt(error_l2_norm_boundary     )
      error_h1_semi_norm_boundary= sqrt(error_h1_semi_norm_boundary)
      h1_semi_norm_boundary      = sqrt(h1_semi_norm_boundary      )
      l2_norm_boundary           = sqrt(l2_norm_boundary           )
    end if

  end subroutine unfitted_solution_checker_compute_error_norms

end module unfitted_solution_checker_names
!***************************************************************************************************
