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
module stokes_cG_discrete_integration_names
  use fempar_names
  use unfitted_temporary_names
  use stokes_analytical_functions_names
  use unfitted_triangulations_names
  use unfitted_fe_spaces_names
  use piecewise_cell_map_names
  use blas77_interfaces_names
  use gen_eigenvalue_solver_names

  implicit none
# include "debug.i90"
  private

  integer(ip), parameter, public :: U_FIELD_ID = 1
  integer(ip), parameter, public :: P_FIELD_ID = 2
  
  type, extends(discrete_integration_t) :: stokes_cG_discrete_integration_t
     type(stokes_analytical_functions_t), pointer :: analytical_functions => NULL()
     type(fe_function_t)                          , pointer :: fe_function => NULL()    
     logical :: unfitted_boundary_is_dirichlet = .true.
     logical :: is_constant_nitches_beta       = .false.
     logical :: use_face_stabilization = .false.
     logical, pointer :: is_in_aggregate(:) => null()
     integer(ip), pointer :: aggregate_root_ids(:) => null()
     real(rp) :: face_stab_constant = 1.0_rp
     real(rp) :: cell_stab_constant = 1.0_rp
   contains
     procedure :: set_analytical_functions
     procedure :: set_stabilization_constants
     procedure :: set_fe_function
     procedure :: set_unfitted_boundary_is_dirichlet
     procedure :: set_is_constant_nitches_beta
     procedure :: set_is_in_aggregate
     procedure :: set_aggregate_root_ids
     procedure :: set_use_face_stabilization
     procedure :: integrate_galerkin
  end type stokes_cG_discrete_integration_t

  public :: stokes_cG_discrete_integration_t

contains

!========================================================================================
  subroutine set_analytical_functions ( this, analytical_functions )
     implicit none
     class(stokes_cG_discrete_integration_t)    , intent(inout) :: this
     type(stokes_analytical_functions_t), target, intent(in)    :: analytical_functions
     this%analytical_functions => analytical_functions
  end subroutine set_analytical_functions

!========================================================================================
  subroutine set_stabilization_constants( this, face_stab_constant)
     implicit none
     class(stokes_cG_discrete_integration_t)    , intent(inout) :: this
     real(rp), intent(in) :: face_stab_constant
     this%face_stab_constant = face_stab_constant
  end subroutine set_stabilization_constants

!========================================================================================
  subroutine set_unfitted_boundary_is_dirichlet ( this, is_dirichlet )
     implicit none
     class(stokes_cG_discrete_integration_t)    , intent(inout) :: this
     logical, intent(in) :: is_dirichlet
     this%unfitted_boundary_is_dirichlet = is_dirichlet
  end subroutine set_unfitted_boundary_is_dirichlet

!========================================================================================
  subroutine set_is_constant_nitches_beta ( this, is_constant )
     implicit none
     class(stokes_cG_discrete_integration_t)    , intent(inout) :: this
     logical, intent(in) :: is_constant
     this%is_constant_nitches_beta = is_constant
  end subroutine set_is_constant_nitches_beta

!========================================================================================
  subroutine set_is_in_aggregate ( this, is_in_aggregate )
     implicit none
     class(stokes_cG_discrete_integration_t)    , intent(inout) :: this
     logical, target, intent(in) :: is_in_aggregate(:)
     this%is_in_aggregate => is_in_aggregate
  end subroutine set_is_in_aggregate

!========================================================================================
  subroutine set_aggregate_root_ids ( this, aggregate_root_ids )
     implicit none
     class(stokes_cG_discrete_integration_t)    , intent(inout) :: this
     integer(ip), target, intent(in) :: aggregate_root_ids(:)
     this%aggregate_root_ids => aggregate_root_ids
  end subroutine set_aggregate_root_ids

!========================================================================================
  subroutine set_use_face_stabilization ( this, use_face_stabilization )
     implicit none
     class(stokes_cG_discrete_integration_t)    , intent(inout) :: this
     logical, intent(in) :: use_face_stabilization
     this%use_face_stabilization = use_face_stabilization
  end subroutine set_use_face_stabilization

!========================================================================================
  subroutine set_fe_function (this, fe_function)
     implicit none
     class(stokes_cG_discrete_integration_t)       , intent(inout) :: this
     type(fe_function_t)                             , target, intent(in)    :: fe_function
     this%fe_function => fe_function
  end subroutine set_fe_function

!========================================================================================
  subroutine integrate_galerkin ( this, fe_space, assembler )
    implicit none
    class(stokes_cG_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)               , intent(inout) :: fe_space
    class(assembler_t)                     , intent(inout) :: assembler

    class(fe_cell_iterator_t), allocatable  :: fe
    type(cell_map_t)         , pointer      :: cell_map
    type(quadrature_t)       , pointer      :: quad
    type(point_t)            , pointer      :: quad_coords(:)
    type(cell_integrator_t)  , pointer      :: cell_int_u
    type(cell_integrator_t)  , pointer      :: cell_int_p
    type(vector_field_t)     , allocatable  :: shape_values_u(:,:)
    type(tensor_field_t)     , allocatable  :: shape_gradients_u(:,:)
    real(rp)                 , allocatable  :: shape_values_p(:,:)
    real(rp)                 , allocatable  :: elmat(:,:), elvec(:)
    class(vector_function_t) , pointer      :: source_term_u
    class(scalar_function_t) , pointer      :: source_term_p
    type(vector_field_t)                    :: source_term_value_u
    real(rp)                                :: source_term_value_p
    type(tensor_field_t)                    :: epsi_u
    type(tensor_field_t)                    :: epsi_v
    real(rp)                                :: div_u
    real(rp)                                :: div_v

    integer(ip)  :: istat
    integer(ip)  :: qpoint, num_quad_points
    integer(ip)  :: idof, jdof, num_dofs
    integer(ip)  :: idof_u, jdof_u
    integer(ip)  :: idof_p, jdof_p
    real(rp)     :: dV

    real(rp), parameter :: viscosity = 1.0

    !! For Nitsche
    real(rp)    :: dS
    real(rp)    :: tau
    class(vector_function_t)   , pointer      :: exact_sol_u
    type(vector_field_t)                      :: exact_sol_value_u
    type(piecewise_cell_map_t) , pointer      :: pw_cell_map
    type(vector_field_t)       , allocatable  :: boundary_shape_values_u(:,:)
    type(tensor_field_t)       , allocatable  :: boundary_shape_gradients_u(:,:)
    real(rp)                   , allocatable  :: boundary_shape_values_p(:,:)
    type(vector_field_t)                      :: normal_vec

    ! For Neumann facet integration
    class(fe_facet_iterator_t), allocatable :: fe_facet
    type(tensor_field_t) :: exact_sol_gradient_u
    real(rp) :: exact_sol_value_p
    class(scalar_function_t) , pointer      :: exact_sol_p
    type(vector_field_t)              :: normals(2)
    
    integer(ip) :: icell, ncells, ipo
    class(triangulation_t), pointer :: triangulation
    integer(ip), parameter :: Npo = 20
    integer(ip) :: print_points(Npo)

    real(rp), allocatable :: facemat(:,:,:,:), facevec(:,:)
    real(rp), allocatable, target :: shape_values_first(:,:), shape_values_second(:,:)
    real(rp), pointer :: shape_values_ineigh(:,:),shape_values_jneigh(:,:)
    real(rp) :: h_length
    integer(ip) :: ineigh
    integer(ip) :: jneigh
    integer(ip) :: n_interior_facets
    integer(ip) :: n_stabilized_facets
    logical :: stabilize_facet
    class(vef_iterator_t), allocatable :: vef
    class(cell_iterator_t), allocatable :: cell_root
    class(cell_iterator_t), allocatable :: cell
    integer(ip) :: ivef
    integer(ip) :: root2_id

    !class(scalar_function_t) , pointer      :: exact_sol
    !type(piecewise_cell_map_t) , pointer :: pw_cell_map
    !real(rp)            , allocatable  :: boundary_shape_values(:,:)
    !type(vector_field_t), allocatable  :: boundary_shape_gradients(:,:)
    !type(vector_field_t)               :: exact_gradient_gp
    !type(vector_field_t)               :: normal_vec
    !real(rp)                           :: normal_d
    !class(reference_fe_t), pointer :: ref_fe
    !class(quadrature_t), pointer :: nodal_quad
    !real(rp), allocatable :: elmatB(:,:), elmatV(:,:), elmatB_pre(:,:)
    !real(rp), allocatable, target :: shape2mono(:,:)
    !real(rp), pointer :: shape2mono_fixed(:,:)
    !real(rp), parameter::beta_coef=2.0_rp
    !real(rp) :: beta
    !real(rp) :: exact_sol_gp
    !real(rp), pointer :: lambdas(:,:)
    !type(gen_eigenvalue_solver_t) :: eigs

    assert (associated(this%analytical_functions))
    assert (associated(this%fe_function))

    call fe_space%create_fe_cell_iterator(fe)

    source_term_u => this%analytical_functions%get_source_term_u()
    source_term_p => this%analytical_functions%get_source_term_p()
    exact_sol_u   => this%analytical_functions%get_solution_function_u()
    exact_sol_p   => this%analytical_functions%get_solution_function_p()

    num_dofs = fe_space%get_max_num_dofs_on_a_cell()
    call memalloc ( num_dofs, num_dofs, elmat, __FILE__, __LINE__ )
    call memalloc ( num_dofs, elvec, __FILE__, __LINE__ )

    call memalloc ( num_dofs, num_dofs, 2, 2, facemat, __FILE__, __LINE__ )
    call memalloc ( num_dofs,              2, facevec, __FILE__, __LINE__ )


    triangulation => fe_space%get_triangulation()
    ncells = triangulation%get_num_cells()
    do ipo = 1, Npo
      icell = (ipo*ncells)/Npo
      print_points(ipo) = icell
    end do
    write(*,*) 'Discrete integration ...'; flush(stdout)
    write(*,*) 'num. cells: ', ncells, ' (this can take a while)'; flush(stdout)
    icell = 0
    ipo = 0
    write(*,'(F8.1x1a1)') 100.0*(real(ipo)/real(Npo)), "%"; flush(stdout)

    call fe%first()
    do while ( .not. fe%has_finished() )


       call fe%update_integration()

       quad            => fe%get_quadrature()
       num_quad_points =  quad%get_num_quadrature_points()
       cell_map        => fe%get_cell_map()
       cell_int_u      => fe%get_cell_integrator(U_FIELD_ID)
       cell_int_p      => fe%get_cell_integrator(P_FIELD_ID)
       quad_coords     => cell_map%get_quadrature_points_coordinates()

       call cell_int_u%get_values(shape_values_u)
       call cell_int_u%get_gradients(shape_gradients_u)
       call cell_int_p%get_values(shape_values_p)

       elmat(:,:) = 0.0_rp
       elvec(:)   = 0.0_rp
       do qpoint = 1, num_quad_points

         dV = cell_map%get_det_jacobian(qpoint) * quad%get_weight(qpoint)

         do idof_u = 1, fe%get_num_dofs_field(U_FIELD_ID)
           idof = idof_u
           epsi_v = (shape_gradients_u(idof_u,qpoint))
           div_v  = trace(epsi_v)
           do jdof_u = 1, fe%get_num_dofs_field(U_FIELD_ID)
             jdof = jdof_u
             epsi_u = (shape_gradients_u(jdof_u,qpoint))
             elmat(idof,jdof) = elmat(idof,jdof) + dV * double_contract(epsi_v,viscosity*epsi_u)
           end do
           do jdof_p = 1, fe%get_num_dofs_field(P_FIELD_ID)
             jdof = jdof_p + fe%get_num_dofs_field(U_FIELD_ID)
             elmat(idof,jdof) = elmat(idof,jdof) + dV * div_v * shape_values_p(jdof_p,qpoint)
           end do
         end do

         do idof_p = 1, fe%get_num_dofs_field(P_FIELD_ID)
           idof = idof_p + fe%get_num_dofs_field(U_FIELD_ID)
           do jdof_u = 1, fe%get_num_dofs_field(U_FIELD_ID)
             jdof = jdof_u
             epsi_u = (shape_gradients_u(jdof_u,qpoint))
             div_u  = trace(epsi_u)
             elmat(idof,jdof) = elmat(idof,jdof) + dV * shape_values_p(idof_p,qpoint) * div_u
           end do
         end do

         call source_term_u%get_value(quad_coords(qpoint),source_term_value_u)
         do idof_u = 1, fe%get_num_dofs_field(U_FIELD_ID)
           idof = idof_u
           elvec(idof) = elvec(idof) + dV * shape_values_u(idof_u,qpoint) * source_term_value_u
         end do

         call source_term_p%get_value(quad_coords(qpoint),source_term_value_p)
         do idof_p = 1, fe%get_num_dofs_field(P_FIELD_ID)
           idof = idof_p + fe%get_num_dofs_field(U_FIELD_ID)
           elvec(idof) = elvec(idof) + dV * shape_values_p(idof_p,qpoint) * source_term_value_p
         end do

       end do

       if (fe%is_cut()) then

         call fe%update_boundary_integration()

         quad            => fe%get_boundary_quadrature()
         num_quad_points = quad%get_num_quadrature_points()
         pw_cell_map     => fe%get_boundary_piecewise_cell_map()
         quad_coords     => pw_cell_map%get_quadrature_points_coordinates()
         cell_int_u      => fe%get_boundary_cell_integrator(U_FIELD_ID)
         cell_int_p      => fe%get_boundary_cell_integrator(P_FIELD_ID)

         call cell_int_u%get_values(boundary_shape_values_u)
         call cell_int_u%get_gradients(boundary_shape_gradients_u)
         call cell_int_p%get_values(boundary_shape_values_p)

         tau = 100.0/cell_map%compute_h(1)

         do qpoint = 1, num_quad_points

           dS = pw_cell_map%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
           call pw_cell_map%get_normal(qpoint,normal_vec)
           call exact_sol_u%get_value(quad_coords(qpoint),exact_sol_value_u)

           do idof_u = 1, fe%get_num_dofs_field(U_FIELD_ID)
             idof = idof_u
             do jdof_u = 1, fe%get_num_dofs_field(U_FIELD_ID)
               jdof = jdof_u
               elmat(idof,jdof) = elmat(idof,jdof) + dS * tau * boundary_shape_values_u(idof_u,qpoint) * boundary_shape_values_u(jdof_u,qpoint)
               elmat(idof,jdof) = elmat(idof,jdof) - dS * ( viscosity * boundary_shape_gradients_u(jdof_u,qpoint) * boundary_shape_values_u(idof_u,qpoint) ) * normal_vec
               elmat(idof,jdof) = elmat(idof,jdof) - dS * ( viscosity * boundary_shape_gradients_u(idof_u,qpoint) * boundary_shape_values_u(jdof_u,qpoint) ) * normal_vec
             end do
             do jdof_p = 1, fe%get_num_dofs_field(P_FIELD_ID)
               jdof = jdof_p + fe%get_num_dofs_field(U_FIELD_ID)
               elmat(idof,jdof) = elmat(idof,jdof) - dS * ( boundary_shape_values_p(jdof_p,qpoint) * boundary_shape_values_u(idof_u,qpoint) ) * normal_vec
             end do
           end do

           do idof_p = 1, fe%get_num_dofs_field(P_FIELD_ID)
             idof = idof_p + fe%get_num_dofs_field(U_FIELD_ID)
             do jdof_u = 1, fe%get_num_dofs_field(U_FIELD_ID)
               jdof = jdof_u
               elmat(idof,jdof) = elmat(idof,jdof) - dS * ( boundary_shape_values_p(idof_p,qpoint) * boundary_shape_values_u(jdof_u,qpoint) ) * normal_vec
             end do
           end do

           do idof_u = 1, fe%get_num_dofs_field(U_FIELD_ID)
             idof = idof_u
             elvec(idof) = elvec(idof) + dS * tau * boundary_shape_values_u(idof_u,qpoint) * exact_sol_value_u
             elvec(idof) = elvec(idof) - dS * ( viscosity * boundary_shape_gradients_u(idof_u,qpoint) * exact_sol_value_u ) * normal_vec
           end do

           do idof_p = 1, fe%get_num_dofs_field(P_FIELD_ID)
             idof = idof_p + fe%get_num_dofs_field(U_FIELD_ID)
             elvec(idof) = elvec(idof) - dS * ( boundary_shape_values_p(idof_p,qpoint) * exact_sol_value_u ) * normal_vec
           end do

         end do

       end if

       call fe%assembly( this%fe_function, elmat, elvec, assembler )
       call fe%next()

       icell = icell + 1
       if ( any(icell == print_points(:))  ) then
         ipo = ipo + 1
         write(*,'(F8.1x1a1)') 100.0*(real(ipo)/real(Npo)), "%"; flush(stdout)
       end if

    end do

    call fe_space%create_fe_facet_iterator(fe_facet)

    if (this%use_face_stabilization) then
      ! Integrate face stabilization terms
      write(*,*) 'Computing face stabilization ...'
      assert (associated(this%is_in_aggregate))
      assert (associated(this%aggregate_root_ids))

      call triangulation%create_vef_iterator(vef)
      call triangulation%create_cell_iterator(cell_root)
      call triangulation%create_cell_iterator(cell)

      n_interior_facets = 0_ip
      n_stabilized_facets = 0_ip
      call fe_facet%first()
      do while ( .not. fe_facet%has_finished() ) 

        if ( fe_facet%is_at_field_interior(1) ) then

          n_interior_facets = n_interior_facets + 1

          assert( fe_facet%get_num_cells_around() == 2 )

          stabilize_facet = .false.

          call fe_facet%get_cell_around(1,cell)
          call cell_root%set_gid( this%aggregate_root_ids(cell%get_gid()) )
          stabilize_facet = stabilize_facet .or. this%is_in_aggregate(cell%get_gid())

          call fe_facet%get_cell_around(2,cell)
          root2_id = this%aggregate_root_ids(cell%get_gid())
          stabilize_facet = stabilize_facet .or. this%is_in_aggregate(cell%get_gid())

          if (cell_root%get_gid() == root2_id) then
            stabilize_facet = .false.
          end if

          if (triangulation%get_num_dims() == 3) then
            if (stabilize_facet) then
              loop1: do ivef = 1,cell_root%get_num_vefs()
                call cell_root%get_vef(ivef,vef)
                if ( vef%is_facet() ) then
                  if ( vef%get_num_cells_around() == 2 ) then
                    do jneigh =1,vef%get_num_cells_around()
                      call vef%get_cell_around(jneigh,cell)
                      if ( cell_root%get_gid() /= cell%get_gid() ) then
                        if (cell%get_gid() == root2_id) then
                          stabilize_facet = .false.
                          exit loop1
                        end if
                      end if
                    end do
                  end if
                end if
              end do loop1
            end if
          end if

          if (stabilize_facet) then

            n_stabilized_facets = n_stabilized_facets + 1

            call fe_facet%update_integration()    
            quad            => fe_facet%get_quadrature()
            num_quad_points = quad%get_num_quadrature_points()

            call fe_facet%get_values(1,shape_values_first , P_FIELD_ID)
            call fe_facet%get_values(2,shape_values_second, P_FIELD_ID)

            facemat(:,:,:,:) = 0.0_rp
            facevec(:,:) = 0.0_rp
            do qpoint = 1, num_quad_points

              call fe_facet%get_normals(qpoint,normals)
              h_length = fe_facet%compute_characteristic_length(qpoint)
              dS = fe_facet%get_det_jacobian(qpoint) * quad%get_weight(qpoint)

              do ineigh = 1, fe_facet%get_num_cells_around()
                if (ineigh==1) then
                  shape_values_ineigh    => shape_values_first
                else if (ineigh==2) then
                  shape_values_ineigh    => shape_values_second
                end if

                do jneigh = 1, fe_facet%get_num_cells_around()

                  if (jneigh==1) then
                    shape_values_jneigh    => shape_values_first
                  else if (jneigh==2) then
                    shape_values_jneigh    => shape_values_second
                  end if

                  do idof_p = 1, fe_facet%get_num_dofs_field(ineigh,P_FIELD_ID)
                    idof = idof_p + fe_facet%get_num_dofs_field(ineigh,U_FIELD_ID)
                    do jdof_p = 1, fe_facet%get_num_dofs_field(jneigh,P_FIELD_ID)
                      jdof = jdof_p + fe_facet%get_num_dofs_field(jneigh,U_FIELD_ID) 
                      facemat(idof,jdof,ineigh,jneigh) = facemat(idof,jdof,ineigh,jneigh) +     &
                        & dS * this%face_stab_constant * h_length &
                        & * shape_values_jneigh(jdof_p,qpoint) &
                        & * shape_values_ineigh(idof_p,qpoint) &
                        & * normals(ineigh)*normals(jneigh) 
                    end do
                  end do

                end do
              end do

            end do

            call fe_facet%assembly( facemat, facevec, assembler )

          end if
        end if

        call fe_facet%next()
      end do
      write(*,*) 'Computing face stabilization ... OK'
      write(*,*) 'n_interior_facets: ', n_interior_facets
      write(*,*) 'n_stabilized_facets: ', n_stabilized_facets
      call triangulation%free_vef_iterator(vef)
      call triangulation%free_cell_iterator(cell_root)
      call triangulation%free_cell_iterator(cell)
    end if

    ! Integrate Neumann boundary conditions
    ! Loop in faces
    call fe_facet%first()
    do while ( .not. fe_facet%has_finished() )

      ! Skip faces that are not in the Neumann boundary
      if ( fe_facet%get_set_id() /= -1 ) then
        call fe_facet%next(); cycle
      end if

      ! Update FE-integration related data structures
      call fe_facet%update_integration()

      quad            => fe_facet%get_quadrature()
      num_quad_points = quad%get_num_quadrature_points()

      ! Get quadrature coordinates to evaluate boundary value
      quad_coords => fe_facet%get_quadrature_points_coordinates()

      ! Get shape functions at quadrature points
      call fe_facet%get_values(1,shape_values_u,U_FIELD_ID)

      ! Compute element vector
      facemat(:,:,:,:) = 0.0_rp
      facevec(:,:) = 0.0_rp
      do qpoint = 1, num_quad_points

        dS = fe_facet%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
        call fe_facet%get_normals(qpoint,normals)
        call exact_sol_u%get_gradient(quad_coords(qpoint),exact_sol_gradient_u)
        call exact_sol_p%get_value(quad_coords(qpoint),exact_sol_value_p)

        do idof_u = 1, fe_facet%get_num_dofs_field(1,U_FIELD_ID)
           idof = idof_u
           facevec(idof,1) = facevec(idof,1) + dS * ( viscosity * exact_sol_gradient_u * shape_values_u(idof_u,qpoint) )*normals(1)
           facevec(idof,1) = facevec(idof,1) + dS * ( exact_sol_value_p * shape_values_u(idof_u,qpoint) )*normals(1)
        end do

      end do

      call fe_facet%assembly( facemat, facevec, assembler )

      call fe_facet%next()
    end do

    call fe_space%free_fe_facet_iterator(fe_facet)


    if (allocated(shape_values_u    )) deallocate  (shape_values_u          , stat=istat); check(istat==0)
    if (allocated(shape_gradients_u )) deallocate  (shape_gradients_u       , stat=istat); check(istat==0)
    if (allocated(shape_values_p    )) call memfree(shape_values_p          , __FILE__, __LINE__)

    if (allocated(boundary_shape_values_u    )) deallocate  (boundary_shape_values_u          , stat=istat); check(istat==0)
    if (allocated(boundary_shape_gradients_u )) deallocate  (boundary_shape_gradients_u       , stat=istat); check(istat==0)
    if (allocated(boundary_shape_values_p    )) call memfree(boundary_shape_values_p          , __FILE__, __LINE__)

    !if (allocated(boundary_shape_values   )) call memfree(boundary_shape_values   , __FILE__, __LINE__)
    !if (allocated(boundary_shape_gradients)) deallocate  (boundary_shape_gradients, stat=istat); check(istat==0);
    !call memfree ( elmatB_pre, __FILE__, __LINE__ )
    !call memfree ( elmatB, __FILE__, __LINE__ )
    !call memfree ( elmatV, __FILE__, __LINE__ )
    !call memfree ( shape2mono, __FILE__, __LINE__ )
    !call eigs%free()

    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )
    call memfree ( facemat, __FILE__, __LINE__ )
    call memfree ( facevec, __FILE__, __LINE__ )

    if (allocated(shape_values_first)) then
      call memfree(shape_values_first, __FILE__, __LINE__) 
    end if

    if (allocated(shape_values_second)) then
      call memfree(shape_values_second, __FILE__, __LINE__) 
    end if

    call fe_space%free_fe_cell_iterator(fe)


    write(*,*) 'Discrete integration ... OK'; flush(stdout)
  end subroutine integrate_galerkin

end module stokes_cG_discrete_integration_names
