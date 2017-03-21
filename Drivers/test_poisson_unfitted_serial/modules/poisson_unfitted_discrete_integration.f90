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
module poisson_unfitted_cG_discrete_integration_names
  use fempar_names
  use poisson_unfitted_analytical_functions_names
  use serial_unfitted_fe_space_names
  use piecewise_fe_map_names
  use blas77_interfaces_names
  use gen_eigenvalue_solver_names

  implicit none
# include "debug.i90"
  private
  type, extends(discrete_integration_t) :: poisson_unfitted_cG_discrete_integration_t
     type(poisson_unfitted_analytical_functions_t), pointer :: analytical_functions => NULL()
   contains
     procedure :: set_analytical_functions
     procedure :: integrate
  end type poisson_unfitted_cG_discrete_integration_t

  public :: poisson_unfitted_cG_discrete_integration_t

contains

!========================================================================================
! TODO @fverdugo FEMPAR PRIORITY MEDIUM EFFORT HIGH
! A better way to do this?
! This subroutine assumes ref elem with hex topology and linear order
! Where to put this info? At the reference element?
subroutine evaluate_monomials(points,monomials)
  implicit none
  type(quadrature_t), intent(in) :: points
  real(rp), allocatable, intent(inout) :: monomials(:,:)
  integer(ip) :: q_point
  real(rp), pointer :: quad_coords(:,:)
  assert(allocated(monomials))
  assert(size(monomials,1)==points%get_number_quadrature_points())
  quad_coords => points%get_coordinates()
  select case(points%get_number_dimensions())
    case(1)
      assert(size(monomials,2)==2)
      do q_point = 1, points%get_number_quadrature_points()
        monomials(q_point,1) = 1 ! 1
        monomials(q_point,2) = quad_coords(1,q_point) ! x
      end do
    case(2)
      assert(size(monomials,2)==4)
      do q_point = 1, points%get_number_quadrature_points()
        monomials(q_point,1) = 1 ! 1
        monomials(q_point,2) = quad_coords(1,q_point) ! x
        monomials(q_point,3) = quad_coords(2,q_point) ! y
        monomials(q_point,4) = quad_coords(1,q_point)*quad_coords(2,q_point) ! xy
      end do
    case(3)
      assert(size(monomials,2)==8)
      do q_point = 1, points%get_number_quadrature_points()
        monomials(q_point,1) = 1 ! 1
        monomials(q_point,2) = quad_coords(1,q_point) ! x
        monomials(q_point,3) = quad_coords(2,q_point) ! y
        monomials(q_point,5) = quad_coords(1,q_point)*quad_coords(2,q_point) ! xy
        monomials(q_point,4) = quad_coords(3,q_point) ! z
        monomials(q_point,6) = quad_coords(1,q_point)*quad_coords(3,q_point) ! xz
        monomials(q_point,7) = quad_coords(2,q_point)*quad_coords(3,q_point) ! yz
        monomials(q_point,8) = quad_coords(1,q_point)*quad_coords(2,q_point)*quad_coords(3,q_point) ! xyz
      end do
    case default
      check(.false.)
  end select
end subroutine evaluate_monomials

!========================================================================================
  subroutine At_times_B_times_A(A,B,C)
  !TODO @fverdugo FEMPAR PRIORITY HIGH EFFORT LOW
  ! Is there a way in FEMPAR to multiply two dense matrices?
  ! This can be optimized
    implicit none
    real(rp), intent(in)    :: A(:,:)
    real(rp), intent(in)    :: B(:,:)
    real(rp), intent(inout) :: C(:,:)
    integer(ip) :: i,j,m,n
    assert( lbound(C,1)==lbound(A,2) .and. ubound(C,1)==ubound(A,2) )
    assert( lbound(A,1)==lbound(B,1) .and. ubound(A,1)==ubound(B,1) )
    assert( lbound(B,2)==lbound(A,1) .and. ubound(B,2)==ubound(A,1) )
    assert( lbound(C,2)==lbound(A,2) .and. ubound(C,2)==ubound(A,2) )
    do i = lbound(C,1), ubound(C,1)
      do j = lbound(C,2), ubound(C,2)
        C(i,j) = 0
        do m = lbound(A,1), ubound(A,1)
          do n = lbound(A,1), ubound(A,1)
            C(i,j) = C(i,j) + A(m,i)*B(m,n)*A(n,j)
          end do
        end do
      end do
    end do
  end subroutine At_times_B_times_A

!========================================================================================
  subroutine set_analytical_functions ( this, analytical_functions )
     implicit none
     class(poisson_unfitted_cG_discrete_integration_t)    , intent(inout) :: this
     type(poisson_unfitted_analytical_functions_t), target, intent(in)    :: analytical_functions
     this%analytical_functions => analytical_functions
  end subroutine set_analytical_functions

!========================================================================================
  subroutine integrate ( this, fe_space, matrix_array_assembler )
    ! TODO @fverdugo
    ! Clean up of this routine required
    implicit none
    class(poisson_unfitted_cG_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)         , intent(inout) :: fe_space
    class(matrix_array_assembler_t)      , intent(inout) :: matrix_array_assembler

    ! FE space traversal-related data types
    ! TODO We need this because the accesors and iterators are not polymorphic
    type(unfitted_fe_iterator_t) :: fe_iterator
    type(unfitted_fe_accessor_t) :: fe

    ! FE integration-related data types
    type(fe_map_t)           , pointer :: fe_map
    type(piecewise_fe_map_t) , pointer :: pw_fe_map
    type(quadrature_t)       , pointer :: quad
    type(point_t)            , pointer :: quad_coords(:)
    type(volume_integrator_t), pointer :: vol_int
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
    real(rp)     :: source_term_value

    integer(ip)  :: number_fields

    integer(ip), pointer :: field_blocks(:)
    logical    , pointer :: field_coupling(:,:)

    type(i1p_t), allocatable :: elem2dof(:)
    integer(ip), allocatable :: num_dofs_per_field(:)
    class(scalar_function_t), pointer :: source_term
    class(scalar_function_t), pointer :: exact_sol

    ! For Nitsche
    class(reference_fe_t), pointer :: ref_fe
    class(quadrature_t), pointer :: nodal_quad
    real(rp), allocatable :: elmatB(:,:), elmatV(:,:), elmatB_aux(:,:)
    real(rp), allocatable :: elmatB_work(:,:), elmatV_work(:,:)! TODO remove this ones, only for debug
    real(rp), allocatable, target :: shape2mono(:,:)
    real(rp), pointer :: shape2mono_fixed(:,:)
    real(rp), allocatable :: eigenvals(:)
    integer(ip) :: lwork
    real(rp), allocatable :: work(:)
    real(rp), parameter::beta_coef=2.0_rp
    real(rp) :: beta
    real(rp) :: exact_sol_gp
    real(rp), pointer :: lambdas(:,:)
    real(rp) :: volele
    type(gen_eigenvalue_solver_t) :: eigs

    assert (associated(this%analytical_functions))


    ! TODO @fverdugo FEMPAR PRIORITY HIGH EFFORT HIGH
    ! we need this because iterator is not polymorpfic
    select type(fe_space)
      class is (serial_unfitted_fe_space_t)
        fe_iterator = fe_space%create_unfitted_fe_iterator()
      class default
        check(.false.)
    end select

    source_term => this%analytical_functions%get_source_term()
    exact_sol   => this%analytical_functions%get_solution_function()

    number_fields = fe_space%get_number_fields()
    allocate( elem2dof(number_fields), stat=istat); check(istat==0);
    field_blocks => fe_space%get_field_blocks()
    field_coupling => fe_space%get_field_coupling()

    call fe_iterator%current(fe) ! TODO assumes constant element ....
    num_dofs = fe%get_number_dofs()
    call memalloc ( num_dofs, num_dofs, elmat, __FILE__, __LINE__ )
    call memalloc ( num_dofs, elvec, __FILE__, __LINE__ )
    call memalloc ( number_fields, num_dofs_per_field, __FILE__, __LINE__ )
    call fe%get_number_dofs_per_field(num_dofs_per_field)

    !This is for the Nitsche's BCs
    ! TODO @fverdugo DRIVER PRIORITY LOW EFFORT HIGH
    ! We assume same ref element for all cells, and for all fields
    ref_fe => fe%get_reference_fe(1)
    nodal_quad => ref_fe%get_nodal_quadrature()
    call memalloc ( num_dofs, num_dofs  , shape2mono, __FILE__, __LINE__ )
    call evaluate_monomials(nodal_quad,shape2mono)
    ! TODO @fverdugo DRIVER PRIORITY LOW EFFORT HIGH
    ! We assume that the constant monomial is the first
    shape2mono_fixed => shape2mono(:,2:)
    lwork = 3*(num_dofs-1)-1
    ! Allocate the eigenvalue solver
    call eigs%create(num_dofs - 1)

    call memalloc ( num_dofs, num_dofs, elmatB_aux, __FILE__, __LINE__ )
    call memalloc ( num_dofs-1, num_dofs-1, elmatB, __FILE__, __LINE__ )
    call memalloc ( num_dofs-1, num_dofs-1, elmatV, __FILE__, __LINE__ )
    call memalloc ( num_dofs-1, eigenvals, __FILE__, __LINE__ )
    call memalloc ( lwork, work, __FILE__, __LINE__ )
    call memalloc ( num_dofs-1, num_dofs-1, elmatB_work, __FILE__, __LINE__ )
    call memalloc ( num_dofs-1, num_dofs-1, elmatV_work, __FILE__, __LINE__ )

    do while ( .not. fe_iterator%has_finished() )
       ! Get current FE
       call fe_iterator%current(fe)

       ! Assemble only for active elements
       ! TODO a better way to do that? with a custom iterator?
       if (.not. fe%is_active()) then
         call fe_iterator%next()
         cycle
       end if

       ! Update FE-integration related data structures
       call fe%update_integration()

       !WARNING This has to be inside the loop
       quad            => fe%get_quadrature()
       num_quad_points = quad%get_number_quadrature_points()
       fe_map          => fe%get_fe_map()
       vol_int         => fe%get_volume_integrator(1)

       ! Get DoF numbering within current FE
       call fe%get_elem2dof(elem2dof)

       ! Get quadrature coordinates to evaluate source_term
       quad_coords => fe_map%get_quadrature_coordinates()

       ! Compute element matrix and vector
       elmat = 0.0_rp
       elvec = 0.0_rp
       volele = 0.0_rp
       call vol_int%get_gradients(shape_gradients)
       call vol_int%get_values(shape_values)
       do qpoint = 1, num_quad_points
          dV = fe_map%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
          volele = volele + dV
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
       ! TODO to this part only for cut elements
       call fe%update_boundary_integration()

       ! Get info on the unfitted boundary for integrating BCs
       quad            => fe%get_boundary_quadrature()
       num_quad_points = quad%get_number_quadrature_points()
       pw_fe_map       => fe%get_boundary_piecewise_fe_map()
       quad_coords     => pw_fe_map%get_quadrature_points_coordinates()
       vol_int         => fe%get_boundary_volume_integrator(1)
       call vol_int%get_values(boundary_shape_values)
       call vol_int%get_gradients(boundary_shape_gradients)

       ! TODO @fverdugo DRIVER PRIORITY HIGH EFFORT MEDIUM
       ! We assume that the unfitted boundary is Nitche
       if (.false.) then
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
       end if

       ! Nitsche beta

       ! Integrate the matrix associated with the normal derivatives
       elmatB_aux(:,:)=0.0_rp
       do qpoint = 1, num_quad_points
         dS = pw_fe_map%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
         call pw_fe_map%get_normal(qpoint,normal_vec)
          do idof = 1, num_dofs
             do jdof = 1, num_dofs
                ! B_K(i,j) = (n*grad(phi_i),n*grad(phi_j))_{\partial\Omega}
                elmatB_aux(idof,jdof) = elmatB_aux(idof,jdof) + &
                  dS *( (normal_vec*boundary_shape_gradients(jdof,qpoint)) * (normal_vec*boundary_shape_gradients(idof,qpoint)) )
             end do
          end do
       end do

       ! Compute the matrices without the kernel
       call At_times_B_times_A(shape2mono_fixed,elmat,elmatV)
       call At_times_B_times_A(shape2mono_fixed,elmatB_aux,elmatB)

       ! Solve the eigenvalue problem
       lambdas => eigs%solve(elmatB,elmatV,istat)
       if( (istat .ne. 0) ) then
         write(*,*) 'Eig Solver returned code ', istat
         write(*,*) 'shape2mono_fixed ='
         do idof = 1, (num_dofs)
           write(*,*) shape2mono_fixed(idof,:)
         end do
         write(*,*) 'B matrix_pre ='
         do idof = 1, (num_dofs)
           write(*,*) elmatB_aux(idof,:)
         end do
         write(*,*) 'V matrix_pre ='
         do idof = 1, (num_dofs)
           write(*,*) elmat(idof,:)
         end do
         write(*,*) 'B matrix ='
         do idof = 1, (num_dofs-1)
           write(*,*) elmatB(idof,:)
         end do
         write(*,*) 'V matrix ='
         do idof = 1, (num_dofs-1)
           write(*,*) elmatV(idof,:)
         end do
         write(*,*) 'Normalized volume trimmed cell = ', volele/(4.0)
         write(*,*) 'V matrix ='
         do idof = 1, (num_dofs-1)
           write(*,*) lambdas(idof,:)
         end do
       end if
       check(istat == 0)

       !! Solve the eigenvalue problem (old method)
       !! TODO @fverdugo FEMPAT PRIORITY HIGH EFFORT LOW
       !! How to avoid implicit interface warning?
       !! How to know if we have to call the double or single precision version?
       !! Wrap the Lapack call in a subroutine?
       !elmatB_work(:,:) = elmatB(:,:)
       !elmatV_work(:,:) = elmatV(:,:)
       !call DSYGV( 1, 'N', 'L', num_dofs-1, elmatB_work, num_dofs-1, elmatV_work, num_dofs-1, eigenvals, work, lwork, istat )
       !if(istat /= 0) then
       !  write(*,*) 'DSYGV returned code ', istat
       !  write(*,*) 'B matrix ='
       !  do idof = 1, (num_dofs-1)
       !    write(*,*) elmatB(idof,:)
       !  end do
       !  write(*,*) 'V matrix ='
       !  do idof = 1, (num_dofs-1)
       !    write(*,*) elmatV(idof,:)
       !  end do
       !end if
       !check(istat==0)

       !if (abs( maxval(eigenvals) - maxval(lambdas(:,1)))>=1e-10_rp) then
       !  write(*,*) maxval(eigenvals), maxval(lambdas(:,1))
       !end if
       !check(abs( maxval(eigenvals) - maxval(lambdas(:,1)))<1e-10_rp)

       ! Set beta (the maximum eigenvalue is stored in the last position)
       !beta = beta_coef*eigenvals(num_dofs-1)

       ! The eigenvalue should be real. Thus, it is save to take only the real part.
       beta = beta_coef*maxval(lambdas(:,1))
       assert(beta>=0)

       ! Once we have the beta, we can compute Nitsche's terms
       do qpoint = 1, num_quad_points

         ! Get info at quadrature point
         dS = pw_fe_map%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
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

       ! Apply boundary conditions
       call fe%impose_strong_dirichlet_bcs( elmat, elvec )
       call matrix_array_assembler%assembly( number_fields, num_dofs_per_field, elem2dof, field_blocks, field_coupling, elmat, elvec )
       call fe_iterator%next()
    end do

    ! TODO Why these are not allocated??
    call memfree(shape_values, __FILE__, __LINE__)
    call memfree(boundary_shape_values, __FILE__, __LINE__)
    deallocate (shape_gradients, stat=istat); check(istat==0);
    deallocate (boundary_shape_gradients, stat=istat); check(istat==0);

    deallocate (elem2dof, stat=istat); check(istat==0);
    call memfree ( num_dofs_per_field, __FILE__, __LINE__ )
    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )
    call memfree ( elmatB_aux, __FILE__, __LINE__ )
    call memfree ( elmatB, __FILE__, __LINE__ )
    call memfree ( elmatV, __FILE__, __LINE__ )
    call memfree ( elmatB_work, __FILE__, __LINE__ )
    call memfree ( elmatV_work, __FILE__, __LINE__ )
    call memfree ( shape2mono, __FILE__, __LINE__ )
    call memfree ( eigenvals, __FILE__, __LINE__ )
    call memfree ( work, __FILE__, __LINE__ )
    call eigs%free()
  end subroutine integrate

end module poisson_unfitted_cG_discrete_integration_names
