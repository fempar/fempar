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
module poisson_unfitted_discrete_integration_names
  use fempar_names
  use unfitted_temporary_names
  use poisson_unfitted_analytical_functions_names
  use par_test_poisson_unfitted_params_names
  use unfitted_triangulations_names
  use unfitted_fe_spaces_names
  use piecewise_fe_map_names
  use blas77_interfaces_names
  use gen_eigenvalue_solver_names
  
  implicit none
# include "debug.i90"
  private
  type, extends(discrete_integration_t) :: poisson_unfitted_cG_discrete_integration_t
     type(poisson_unfitted_analytical_functions_t), pointer :: analytical_functions => NULL()
     type(fe_function_t)                          , pointer :: fe_function          => NULL()
     type(par_test_poisson_unfitted_params_t)     , pointer :: test_params => null()
   contains
     procedure :: set_analytical_functions
     procedure :: set_fe_function
     procedure :: set_test_params
     procedure :: integrate_galerkin
  end type poisson_unfitted_cG_discrete_integration_t
  
  public :: poisson_unfitted_cG_discrete_integration_t
  
contains
   
  subroutine set_analytical_functions ( this, analytical_functions )
     implicit none
     class(poisson_unfitted_cG_discrete_integration_t)    ,intent(inout)  :: this
     type(poisson_unfitted_analytical_functions_t), target, intent(in)    :: analytical_functions
     this%analytical_functions => analytical_functions
  end subroutine set_analytical_functions

  subroutine set_test_params ( this, test_params )
     implicit none
     class(poisson_unfitted_cG_discrete_integration_t)    ,intent(inout)  :: this
     type(par_test_poisson_unfitted_params_t), target, intent(in)    :: test_params
     this%test_params => test_params
  end subroutine set_test_params
  
  subroutine set_fe_function (this, fe_function)
     implicit none
     class(poisson_unfitted_cG_discrete_integration_t)       , intent(inout) :: this
     type(fe_function_t)                             , target, intent(in)    :: fe_function
     this%fe_function => fe_function
  end subroutine set_fe_function
  
  subroutine integrate_galerkin ( this, fe_space, matrix_array_assembler )
    implicit none
    class(poisson_unfitted_cG_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)         , intent(inout) :: fe_space
    class(matrix_array_assembler_t)      , intent(inout) :: matrix_array_assembler

    ! FE space traversal-related data types
    ! TODO We need this because the accesors and iterators are not polymorphic
    class(fe_iterator_t), allocatable  :: fe

    ! FE integration-related data types
    type(fe_map_t)           , pointer :: fe_map
    type(piecewise_fe_map_t) , pointer :: pw_fe_map
    type(quadrature_t)       , pointer :: quad
    type(point_t)            , pointer :: quad_coords(:)
    type(cell_integrator_t), pointer :: cell_int
    type(vector_field_t), allocatable  :: shape_gradients(:,:)
    real(rp)            , allocatable  :: shape_values(:,:)
    real(rp)            , allocatable  :: boundary_shape_values(:,:)
    type(vector_field_t), allocatable  :: boundary_shape_gradients(:,:)
    type(vector_field_t)               :: exact_gradient_gp
    type(vector_field_t)               :: normal_vec
    real(rp)                           :: normal_d

    ! FE matrix and vector i.e., A_K + f_K
    real(rp), allocatable              :: elmat(:,:), elvec(:)

    integer(ip)  :: istat
    integer(ip)  :: qpoint, num_quad_points
    integer(ip)  :: idof, jdof, num_dofs
    real(rp)     :: dV, dS
    real(rp)     :: volume, surface
    real(rp)     :: source_term_value

    class(scalar_function_t), pointer :: source_term
    class(scalar_function_t), pointer :: exact_sol

    ! For Nitsche
    class(reference_fe_t), pointer :: ref_fe
    class(quadrature_t), pointer :: nodal_quad
    real(rp), allocatable :: elmatB(:,:), elmatV(:,:), elmatB_pre(:,:)
    real(rp), allocatable, target :: shape2mono(:,:)
    real(rp), pointer :: shape2mono_fixed(:,:)
    real(rp) :: beta_coef
    real(rp) :: beta
    real(rp) :: exact_sol_gp
    real(rp), pointer :: lambdas(:,:)
    type(gen_eigenvalue_solver_t) :: eigs

    assert (associated(this%analytical_functions))
    assert (associated(this%test_params))
    assert (associated(this%fe_function))
    beta_coef = this%test_params%get_nitsche_beta_factor()

    ! TODO We will delete this once implemented the fake methods in the father class
    call fe_space%create_fe_iterator(fe)

    source_term => this%analytical_functions%get_source_term()
    exact_sol   => this%analytical_functions%get_solution_function()

    ! Find the first non-void FE
    call fe%first_local_non_void(field_id = 1)
    if (fe%has_finished()) then 
      call fe_space%free_fe_iterator(fe)
      return
    end if
    quad            => fe%get_quadrature()
    num_quad_points = quad%get_number_quadrature_points()

    ! TODO We assume that all non-void FEs are the same...
    num_dofs = fe%get_number_dofs()
    call memalloc ( num_dofs, num_dofs, elmat, __FILE__, __LINE__ )
    call memalloc ( num_dofs, elvec, __FILE__, __LINE__ )

    !This is for the Nitsche's BCs
    ! TODO  We assume same ref element for all cells, and for all fields
    ref_fe => fe%get_reference_fe(1)
    nodal_quad => ref_fe%get_nodal_quadrature()
    call memalloc ( num_dofs, num_dofs  , shape2mono, __FILE__, __LINE__ )
    call evaluate_monomials(nodal_quad,shape2mono,degree=this%analytical_functions%get_degree())
    ! TODO  We assume that the constant monomial is the first
    shape2mono_fixed => shape2mono(:,2:)
    ! Allocate the eigenvalue solver
    call eigs%create(num_dofs - 1)

    call memalloc ( num_dofs, num_dofs, elmatB_pre, __FILE__, __LINE__ )
    call memalloc ( num_dofs-1, num_dofs-1, elmatB, __FILE__, __LINE__ )
    call memalloc ( num_dofs-1, num_dofs-1, elmatV, __FILE__, __LINE__ )

    call fe%first()
    volume = 0.0
    surface = 0.0
    do while ( .not. fe%has_finished() )

       if ( fe%is_local() ) then

         ! Update FE-integration related data structures
         call fe%update_integration()

         !WARNING This has to be inside the loop
         quad            => fe%get_quadrature()
         num_quad_points = quad%get_number_quadrature_points()
         fe_map          => fe%get_fe_map()
         cell_int         => fe%get_cell_integrator(1)
         num_dofs = fe%get_number_dofs()

         ! Get quadrature coordinates to evaluate source_term
         quad_coords => fe_map%get_quadrature_coordinates()

         ! Compute element matrix and vector
         elmat = 0.0_rp
         elvec = 0.0_rp
         call cell_int%get_gradients(shape_gradients)
         call cell_int%get_values(shape_values)
         do qpoint = 1, num_quad_points
            dV = fe_map%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
            volume = volume + dV
            do idof = 1, num_dofs
               do jdof = 1, num_dofs
                  ! A_K(i,j) = (grad(phi_i),grad(phi_j))
                  elmat(idof,jdof) = elmat(idof,jdof) + dV * shape_gradients(jdof,qpoint) * shape_gradients(idof,qpoint)
               end do
            end do

            ! Source term
            call source_term%get_value(quad_coords(qpoint),source_term_value)
            do idof = 1, num_dofs
               elvec(idof) = elvec(idof) + dV * source_term_value * shape_values(idof,qpoint)
            end do
         end do

         ! Update FE boundary integration related data structures
         ! Only for cut elements
         ! TODO @fverdugo FEMPAR PRIORITY LOW EFFORT HIGH
         ! Create iterator for cut and for full elements? Then we can remove this if
         if (fe%is_cut()) then

           call fe%update_boundary_integration()

           ! Get info on the unfitted boundary for integrating BCs
           quad            => fe%get_boundary_quadrature()
           num_quad_points = quad%get_number_quadrature_points()
           pw_fe_map       => fe%get_boundary_piecewise_fe_map()
           quad_coords     => pw_fe_map%get_quadrature_points_coordinates()
           cell_int         => fe%get_boundary_cell_integrator(1)
           call cell_int%get_values(boundary_shape_values)
           call cell_int%get_gradients(boundary_shape_gradients)

           select case (trim(this%test_params%get_unfitted_boundary_type()))
           case ('neumann')

             ! Neumann BCs unfitted boundary
             do qpoint = 1, num_quad_points

               ! Surface measure
               dS = pw_fe_map%get_det_jacobian(qpoint) * quad%get_weight(qpoint)

               ! Value of the gradient of the solution at the boundary
               call exact_sol%get_gradient(quad_coords(qpoint),exact_gradient_gp)

               ! Get the boundary normals
               call pw_fe_map%get_normal(qpoint,normal_vec)

               ! Normal derivative
               ! It is save to do so in 2d only if the 3rd component is set to 0
               ! in at least one of the 2 vectors
               normal_d = normal_vec*exact_gradient_gp

               ! Integration
                do idof = 1, num_dofs
                   elvec(idof) = elvec(idof) + normal_d * boundary_shape_values(idof,qpoint) * dS
                end do

             end do

           case ('dirichlet') ! Nitsche on the unfitted boundary

             ! Nitsche beta

             ! Integrate the matrix associated with the normal derivatives
             elmatB_pre(:,:)=0.0_rp
             do qpoint = 1, num_quad_points
               dS = pw_fe_map%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
               call pw_fe_map%get_normal(qpoint,normal_vec)
                do idof = 1, num_dofs
                   do jdof = 1, num_dofs
                      ! B_K(i,j) = (n*grad(phi_i),n*grad(phi_j))_{\partial\Omega}
                      elmatB_pre(idof,jdof) = elmatB_pre(idof,jdof) + &
                        dS *( (normal_vec*boundary_shape_gradients(jdof,qpoint)) * (normal_vec*boundary_shape_gradients(idof,qpoint)) )
                   end do
                end do
             end do

             ! Compute the matrices without the kernel
             call At_times_B_times_A(shape2mono_fixed,elmat,elmatV)
             call At_times_B_times_A(shape2mono_fixed,elmatB_pre,elmatB)

             ! Solve the eigenvalue problem
             lambdas => eigs%solve(elmatB,elmatV,istat)
             if (istat .ne. 0) then
               write(*,*) 'istat = ', istat
               write(*,*) 'lid   = ', fe%get_lid()
               !write(*,*) 'elmatB = '
               !do idof = 1,size(elmatB,1)
               !  write(*,*) elmatB(idof,:)
               !end do
               !write(*,*) 'elmatV = '
               !do idof = 1,size(elmatV,1)
               !  write(*,*) elmatV(idof,:)
               !end do
             end if
             mcheck(istat == 0,'Failed to solve the generalized eigenvalue problem')

             ! The eigenvalue should be real. Thus, it is save to take only the real part.
             beta = beta_coef*maxval(lambdas(:,1))
             assert(beta>=0)

             ! Once we have the beta, we can compute Nitsche's terms
             do qpoint = 1, num_quad_points

               ! Get info at quadrature point
               dS = pw_fe_map%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
               surface = surface + dS
               call pw_fe_map%get_normal(qpoint,normal_vec)
               call exact_sol%get_value(quad_coords(qpoint),exact_sol_gp)

               ! Elem matrix
                do idof = 1, num_dofs
                   do jdof = 1, num_dofs
                      ! A_K(i,j)=(beta*phi_i,phi_j)_{\partial\Omega} - (phi_i,n*grad(phi_j))_{\partial\Omega}  - (phi_j,n*grad(phi_i))_{\partial\Omega}
                      elmat(idof,jdof) = elmat(idof,jdof) &
                        + dS*beta*(boundary_shape_values(idof,qpoint)*boundary_shape_values(jdof,qpoint)) &
                        - dS*(boundary_shape_values(idof,qpoint)*(normal_vec*boundary_shape_gradients(jdof,qpoint))) &
                        - dS*(boundary_shape_values(jdof,qpoint)*(normal_vec*boundary_shape_gradients(idof,qpoint)))
                   end do
                end do

                ! Elem vector
                do idof = 1, num_dofs
                   ! f_k(i) = (beta*ufun,phi_i)_{\partial\Omega} - (ufun,n*grad(phi_i))_{\partial\Omega}
                   elvec(idof) = elvec(idof) &
                     + dS*beta*exact_sol_gp*boundary_shape_values(idof,qpoint) &
                     - dS*exact_sol_gp*(normal_vec*boundary_shape_gradients(idof,qpoint))
                end do

             end do

           case default
             mcheck(.false.,'Unknown unfitted boundary type `'//trim(this%test_params%get_unfitted_boundary_type())//'`')
           end select

         end if ! Only for cut elems

         call fe%assemble( this%fe_function, elmat, elvec, matrix_array_assembler )

       end if
       call fe%next()

    end do

    !write(*,*) "Domain volume    = ", volume
    !write(*,*) "Boundary surface = ", surface

    ! These are allocated by the vol integrators
    ! We have to check if they are actually allocated before deallocate them,
    ! since we cannot assure that they have been allocated
    if (allocated(shape_values            )) call memfree(shape_values            , __FILE__, __LINE__)
    if (allocated(boundary_shape_values   )) call memfree(boundary_shape_values   , __FILE__, __LINE__)
    if (allocated(shape_gradients         )) deallocate  (shape_gradients         , stat=istat); check(istat==0);
    if (allocated(boundary_shape_gradients)) deallocate  (boundary_shape_gradients, stat=istat); check(istat==0);

    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )
    call memfree ( elmatB_pre, __FILE__, __LINE__ )
    call memfree ( elmatB, __FILE__, __LINE__ )
    call memfree ( elmatV, __FILE__, __LINE__ )
    call memfree ( shape2mono, __FILE__, __LINE__ )
    call eigs%free()
    call fe_space%free_fe_iterator(fe)

  end subroutine integrate_galerkin
  
end module poisson_unfitted_discrete_integration_names
