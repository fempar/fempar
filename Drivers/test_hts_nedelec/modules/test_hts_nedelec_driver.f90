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
module test_hts_nedelec_driver_names
  use fempar_names
  use hts_nedelec_params_names
  use hts_nedelec_analytical_functions_names
  use hts_nedelec_discrete_integration_names
  use hts_nedelec_conditions_names
  use hts_nonlinear_solver_names
  use hts_theta_method_names 
  use base_sparse_matrix_names 
# include "debug.i90"

  implicit none
  private

  integer(ip), parameter :: hts = 1
  integer(ip), parameter :: air = 2
  
  type test_hts_nedelec_driver_t 
     private 

     ! Place-holder for parameter-value set provided through command-line interface
     type(hts_nedelec_params_t)               :: test_params

     ! Cells and lower dimension objects container
     type(p4est_serial_triangulation_t)       :: triangulation
					type(ParameterList_t), pointer           :: parameter_list

     ! Analytical functions of the problem
     type(hts_nedelec_analytical_functions_t) :: problem_functions 

     ! Discrete weak problem integration-related data type instances 
     type(serial_fe_space_t)                  :: fe_space 
     type(p_reference_fe_t) , allocatable     :: reference_fes(:) 
     type(hts_nedelec_discrete_integration_t) :: hts_nedelec_integration
     type(hts_nedelec_conditions_t)           :: hts_nedelec_conditions

     ! Place-holder for the coefficient matrix and RHS of the linear system
     type(fe_affine_operator_t)                  :: fe_affine_operator
     type(coo_sparse_matrix_t)                   :: constraint_matrix
     type(serial_scalar_array_t)                 :: constraint_vector
     
     ! Temporal and Nonlinear solver data type 
     type(hts_nonlinear_solver_t)              :: nonlinear_solver 
     type(theta_method_t)                      :: theta_method 
  
     ! HTS problem solution FE function
     type(fe_function_t)            :: H_current 
     type(fe_function_t)            :: H_previous 
     
     ! Writing solution 
     type(output_handler_t)         :: oh 

   contains
     procedure                  :: run_simulation
     procedure        , private :: parse_command_line_parameters
     procedure        , private :: setup_triangulation
					procedure        , private :: set_cells_for_refinement 
     procedure        , private :: setup_reference_fes
     procedure        , private :: setup_fe_space
     procedure        , private :: setup_system
     procedure        , private :: setup_constraint_matrix 
     procedure        , private :: setup_nonlinear_solver
     procedure        , private :: setup_theta_method 
     procedure        , private :: assemble_system
     procedure        , private :: solve_system
     procedure        , private :: solve_nonlinear_system 
     procedure        , private :: compute_hysteresis_data 
     procedure        , private :: check_solution
     procedure        , private :: initialize_output  
     procedure        , private :: write_time_step_solution
     procedure        , private :: finalize_output 
     procedure        , private :: free
  end type test_hts_nedelec_driver_t

  ! Types
  public :: test_hts_nedelec_driver_t

contains

  subroutine parse_command_line_parameters(this)
    implicit none
    class(test_hts_nedelec_driver_t ), intent(inout) :: this
    call this%test_params%create()
				this%parameter_list => this%test_params%get_values()
  end subroutine parse_command_line_parameters

  ! -----------------------------------------------------------------------------------------------
  subroutine setup_triangulation(this)
    implicit none
    class(test_hts_nedelec_driver_t), intent(inout) :: this
    ! Locals 
    integer(ip)                 , allocatable :: cells_set(:) 
    class(cell_iterator_t)      , allocatable :: cell
    type(point_t), allocatable                :: cell_coordinates(:)
    integer(ip)                               :: inode
    integer(ip)       :: icell, icoord 
    real(rp)          :: cx, cy, cz 
    integer(ip)       :: i, istat, idime  
   	real(rp)          :: R, h, x0, y0, z0
				real(rp)          :: domain(6)
				real(rp)          :: domain_length(0:SPACE_DIM-1) 
				real(rp)          :: hts_domain_length(0:SPACE_DIM-1)
				real(rp)          :: hts_lx, hts_ly, hts_lz 

    ! Create a structured mesh with a custom domain 
				domain_length     = this%test_params%get_domain_length() 
    domain = [ 0.0_rp, domain_length(0), 0.0_rp, domain_length(1), 0.0_rp, domain_length(2) ]  
    istat = this%parameter_list%set(key = hex_mesh_domain_limits_key , value = domain); check(istat==0)
    call this%triangulation%create(this%parameter_list)

				do i = 1, this%test_params%get_num_refinements() 
      call this%set_cells_for_refinement()
      call this%triangulation%refine_and_coarsen()
      call this%triangulation%clear_refinement_and_coarsening_flags()
    end do
				
	if ( this%test_params%get_triangulation_type() == triangulation_generate_structured ) then
					hts_domain_length = this%test_params%get_hts_domain_length() 
					hts_lx = hts_domain_length(0) 
					hts_ly = hts_domain_length(1)
					hts_lz = hts_domain_length(2)
					
    ! Assign subset_id to different cells for the created structured mesh 
    allocate(cells_set(this%triangulation%get_num_cells() ), stat=istat); check(istat==0)
    call this%triangulation%create_cell_iterator(cell)
    allocate(cell_coordinates( cell%get_num_nodes() ) , stat=istat); check(istat==0) 
    
    do while ( .not. cell%has_finished() )
       call cell%get_nodes_coordinates(cell_coordinates)
       ! Compute center of the element coordinates 
       cx = 0.0_rp
       cy = 0.0_rp 
	      cz = 0.0_rp 
       do inode=1,cell%get_num_nodes()  
          cx = cx + cell_coordinates(inode)%get(1)
          cy = cy + cell_coordinates(inode)%get(2)
		        cz = cz + cell_coordinates(inode)%get(3)
       end do
       cx = cx/real(cell%get_num_nodes(),rp)
       cy = cy/real(cell%get_num_nodes(),rp)
	      cz = cz/real(cell%get_num_nodes(),rp)

       ! Select material case: HTS TAPE in the center 
       if ( ( (domain_length(0) - hts_lx)/2.0_rp < cx .and. cx < (domain_length(0) + hts_lx)/2.0_rp ) & 
											.and. ( (domain_length(1) - hts_ly)/2.0_rp < cy .and. cy < (domain_length(1) + hts_ly)/2.0_rp) ) then
          cells_set( cell%get_gid() ) = hts 
       else 
          cells_set( cell%get_gid() ) = air
       end if
    
							call cell%next() 
    end do

    call this%triangulation%fill_cells_set(cells_set)  
    deallocate(cells_set, stat=istat); check(istat==0) 
    deallocate(cell_coordinates, stat=istat); check(istat==0) 
    call this%triangulation%free_cell_iterator(cell)
	
	end if 

  end subroutine setup_triangulation
		
		 subroutine set_cells_for_refinement(this)
    implicit none
    class(test_hts_nedelec_driver_t), intent(inout) :: this
    class(cell_iterator_t)      , allocatable :: cell
    type(point_t), allocatable :: coords(:)
    integer(ip) :: istat, k
    real(rp) ::  x,y
    integer(ip), parameter :: max_num_cell_nodes = 4
    integer(ip), parameter :: max_level = 4
			 
				real(rp)          :: x0, xL, y0, yL, x_eps, y_eps 
				real(rp)          :: domain_length(0:SPACE_DIM-1) 
				real(rp)          :: hts_domain_length(0:SPACE_DIM-1)
				real(rp)          :: epsilon_length(0:SPACE_DIM-1) 

				! Set refinemenent if the cell contains centered hts device 
				domain_length        = this%test_params%get_domain_length()
				hts_domain_length    = this%test_params%get_hts_domain_length()
				epsilon_length       = this%test_params%get_eps_hts_domain_length()
	
				x0 = (domain_length(0)-hts_domain_length(0))/2.0_rp - epsilon_length(0) 
				xL = (domain_length(0)+hts_domain_length(0))/2.0_rp + epsilon_length(0)
				y0 = (domain_length(1)-hts_domain_length(1))/2.0_rp - epsilon_length(1)
				yL = (domain_length(1)+hts_domain_length(1))/2.0_rp + epsilon_length(1) 
				
    call this%triangulation%create_cell_iterator(cell)
    if (this%triangulation%get_num_dims() == 2) then
      allocate(coords(max_num_cell_nodes),stat=istat); check(istat==0)

      do while ( .not. cell%has_finished() )

        call cell%get_nodes_coordinates(coords)
								if ( cell%get_level() < this%test_params%get_num_min_refinements() ) then 
										call cell%set_for_refinement()
								else 
        do k=1,max_num_cell_nodes
								! If cell contains device, refine! 
         if ( ( x0 < coords(k)%get(1) .and. coords(k)%get(1) < xL ) .and. ( y0 < coords(k)%get(2) .and. coords(k)%get(2) < yL ) ) then 
										call cell%set_for_refinement(); exit 
									end if 
        end do
								end if 
																
        call cell%next()
      end do
      deallocate(coords,stat=istat); check(istat==0)

    else if (this%triangulation%get_num_dims() == 3) then    
      ! To define a 3D algorithm 
    else
      mcheck(.false.,'Only for 2D and 3D')

    end if

    call this%triangulation%free_cell_iterator(cell)

  end subroutine set_cells_for_refinement

  ! -----------------------------------------------------------------------------------------------
  subroutine setup_reference_fes(this)
    implicit none
    class(test_hts_nedelec_driver_t), intent(inout) :: this
    integer(ip) :: istat, ivef
				class(vef_iterator_t), allocatable :: vef

    allocate(this%reference_fes(1), stat=istat)
    check(istat==0)
    
    this%reference_fes(1) =  make_reference_fe ( topology = topology_hex,                                          &
                                                 fe_type = fe_type_nedelec,                                        &
                                                 num_dims = this%triangulation%get_num_dims(),      &
                                                 order = this%test_params%get_magnetic_field_reference_fe_order(), &
                                                 field_type = field_type_vector,                                   &
                                                 conformity = .true. ) 
   
    if ( this%test_params%get_triangulation_type() == triangulation_generate_structured ) then
       call this%triangulation%create_vef_iterator(vef)
       do while ( .not. vef%has_finished() )
          ! In the 3D case, vefs asociated to faces 21,22 are Neumann boundary (2D case set_id <= 9)
       !  if ( vef%is_at_boundary() .and. ( vef%get_set_id() .ne. 21 .and. vef%get_set_id() .ne. 22) ) then 
          if ( vef%is_at_boundary() ) then 
		        call vef%set_set_id(1)
          else
             call vef%set_set_id(0)
          end if
          call vef%next()
       end do
       call this%triangulation%free_vef_iterator(vef)
    end if    
 
  end subroutine setup_reference_fes

  ! -----------------------------------------------------------------------------------------------
  subroutine setup_fe_space(this)
    implicit none
    class(test_hts_nedelec_driver_t), intent(inout) :: this

    call this%hts_nedelec_conditions%set_num_dims( this%triangulation%get_num_dims() )
    call this%problem_functions%initialize( H  = this%test_params%get_external_magnetic_field_amplitude(),  &
                                            wH = this%test_params%get_external_magnetic_field_frequency(),  &
                                            J  = this%test_params%get_external_current_amplitude(),         &
                                            wJ = this%test_params%get_external_current_frequency()  )  
    
    call this%fe_space%create( triangulation = this%triangulation, &
                               reference_fes = this%reference_fes, &
                               conditions    = this%hts_nedelec_conditions )

    call this%fe_space%set_up_cell_integration()
    call this%fe_space%set_up_facet_integration() 
       
  end subroutine setup_fe_space 
    
  ! -----------------------------------------------------------------------------------------------
  subroutine setup_system (this)
    implicit none
    class(test_hts_nedelec_driver_t), intent(inout) :: this 
    
    call this%problem_functions%set_num_dims(this%triangulation%get_num_dims())
    call this%hts_nedelec_integration%create( this%theta_method, this%H_current, this%H_previous, &
                                              this%test_params, this%problem_functions%get_source_term() )
    call this%fe_affine_operator%create ( sparse_matrix_storage_format      = csr_format, &
                                          diagonal_blocks_symmetric_storage = [ .false.  ], &
                                          diagonal_blocks_symmetric         = [ .false. ], &
                                          diagonal_blocks_sign              = [ SPARSE_MATRIX_SIGN_UNKNOWN ], &
                                          fe_space                          = this%fe_space,           &
                                          discrete_integration              = this%hts_nedelec_integration )
    
    call this%hts_nedelec_conditions%set_boundary_function_Hx(this%problem_functions%get_boundary_function_Hx())
    call this%hts_nedelec_conditions%set_boundary_function_Hy(this%problem_functions%get_boundary_function_Hy())
    if ( this%triangulation%get_num_dims() == 3) then 
       call this%hts_nedelec_conditions%set_boundary_function_Hz(this%problem_functions%get_boundary_function_Hz())
    end if
    ! Create H_previous with initial time (t0) boundary conditions 
    call this%H_previous%create(this%fe_space)
				call this%fe_space%interpolate_vector_function( 1, this%problem_functions%get_solution(), &
                                                       this%H_previous,                       &  
                                                       time=this%theta_method%get_initial_time() ) 

    call this%fe_space%interpolate_dirichlet_values(this%H_previous, time=this%theta_method%get_initial_time() ) 
    ! Update fe_space to the current time (t1) boundary conditions, create H_current
    call this%H_current%create(this%fe_space)
				call this%fe_space%interpolate_vector_function( 1, this%problem_functions%get_solution(), &
                                                       this%H_current,                        &  
                                                       time=this%theta_method%get_initial_time() ) 
    call this%fe_space%interpolate_dirichlet_values(this%H_current, time=this%theta_method%get_current_time()  )
    
    ! Setup constraint matrix if the problem is defined constrained 
    if (this%test_params%get_apply_current_density_constraint() ) then 
    call this%setup_constraint_matrix() 
    end if
    
  end subroutine setup_system
  
  ! -----------------------------------------------------------------------------------------------
  subroutine setup_constraint_matrix(this)
    implicit none
    class(test_hts_nedelec_driver_t), intent(inout) :: this

    class(matrix_t)      , pointer :: matrix
    type(sparse_matrix_t), pointer :: coefficient_matrix

    ! Integration loop 
    class(fe_cell_iterator_t)     , allocatable :: fe
    class(fe_facet_iterator_t), allocatable :: fe_face 
    integer(ip) :: ielem 
    type(quadrature_t)       , pointer     :: quad
    type(vector_field_t)     , allocatable :: curl_values(:,:)
    integer(ip)                            :: qpoin, num_qpoints, idof 
    type(i1p_t)              , pointer     :: fe_dofs(:)
    integer(ip)                            :: i, inode, vector_size
    integer(ip)                            :: num_dofs, num_fields 
    real(rp)                 , allocatable :: elvec(:), facevec(:) 
    real(rp)                               :: factor  
    integer(ip)  :: istat 


    matrix=> this%fe_affine_operator%get_matrix()
    select type(matrix)
    type is ( sparse_matrix_t) 
       coefficient_matrix => matrix
       class default
       check(.false.)
    end select

    ! Free any dynamic memory that constraint_matrix may have inside
    call this%constraint_matrix%free()
    call this%constraint_vector%free() 

    ! Create constraint matrix (transposed)
    call this%constraint_matrix%create ( num_rows=coefficient_matrix%get_num_rows(), num_cols=1 )
    call this%constraint_vector%create_and_allocate( coefficient_matrix%get_num_rows() )
    call this%constraint_vector%init(0.0_rp) 

        ! Initialize
    call this%fe_space%set_up_cell_integration()
    call this%fe_space%create_fe_cell_iterator(fe)
    
    num_fields         =  this%fe_space%get_num_fields()
    num_dofs              =  fe%get_num_dofs()
    call memalloc ( num_dofs, elvec, __FILE__, __LINE__ )
    allocate( fe_dofs(num_fields), stat=istat); check(istat==0);
    
        ! ================================  2D CASE, integrate over entire HTS section ================
    if ( this%triangulation%get_num_dims() == 2) then  
    
    quad           => fe%get_quadrature()
    num_qpoints =  quad%get_num_quadrature_points()
    
    ! Loop over elements
    do while ( .not. fe%has_finished())

       if ( fe%get_set_id() == hts ) then  
          ! Update finite structures
          call fe%update_integration()		               
          call fe%get_fe_dofs(fe_dofs) 
          call fe%get_curls(curl_values)

          elvec      = 0.0_rp 
          ! Integrate J over the hts subdomain 
          do qpoin=1, num_qpoints
             factor = fe%get_det_jacobian(qpoin) * quad%get_weight(qpoin) 						
             do inode = 1, fe%get_num_dofs_field(1)  
                elvec(inode) = elvec(inode) + factor * curl_values(inode,qpoin)%get(3) 
             end do
          end do

          ! Add element contribution to matrix and vector 
          do i = 1, fe%get_num_dofs_field(1) 
             idof = fe_dofs(1)%p(i) 
             if ( idof > 0 ) then 
                 call this%constraint_matrix%insert( idof, 1, elvec(i) )
                 call this%constraint_vector%add(idof, elvec(i))
             end if
          end do

       end if
       call fe%next()
    end do    
    ! ================================   3D CASE, only integrate over z-normal faces ===================
    elseif ( this%triangulation%get_num_dims() == 3) then 
    
       call this%fe_space%set_up_facet_integration()

       ! Search for the first boundary face
       call this%fe_space%create_fe_facet_iterator(fe_face)
       do while ( .not. fe_face%is_at_boundary() ) 
          call fe_face%next()
       end do

       num_dofs              =  fe%get_num_dofs() 
       call memalloc ( num_dofs, facevec, __FILE__, __LINE__ )
       quad            => fe_face%get_quadrature()
       num_qpoints  =  quad%get_num_quadrature_points()
       call fe_face%get_curls(1,curl_values) 

       do while ( .not. fe_face%has_finished() )
          facevec = 0.0_rp
          if ( fe_face%is_at_boundary() .and. fe_face%get_set_id() == 0 ) then
             
             call fe_face%get_cell_around(1, fe)
             if ( fe%get_set_id() == hts ) then 

                call fe_face%update_integration()    
                do qpoin = 1, num_qpoints
                   factor = fe_face%get_det_jacobian(qpoin) * quad%get_weight(qpoin)
                   do idof = 1, fe%get_num_dofs_field(1) 
                      facevec(idof) = facevec(idof) + factor * curl_values(idof,qpoin)%get(3) 
                   end do
                end do

                call fe_face%get_fe_dofs(1, fe_dofs)

                ! Add element contribution to vector 
                do i = 1, fe%get_num_dofs_field(1) 
                   idof = fe_dofs(1)%p(i) 
                   if ( idof > 0 ) then 
                      call this%constraint_matrix%insert( idof, 1, facevec(i) )
                      call this%constraint_vector%add(idof, facevec(i))
                   end if
                end do

             end if
          end if
          call fe_face%next()
       end do
       call memfree ( facevec, __FILE__, __LINE__ )
    end if 
    call this%fe_space%free_fe_cell_iterator(fe)
    call this%fe_space%free_fe_facet_iterator(fe_face)
    ! Sum duplicates, re-order by rows, and leave the matrix in a final state
    call this%constraint_matrix%sort_and_compress()
    ! call this%constraint_vector%print(6) 
    ! =============================================================================================
    call memfree ( elvec, __FILE__, __LINE__ )
    deallocate (fe_dofs, stat=istat); check(istat==0)
    deallocate (curl_values, stat=istat); check(istat==0)

  end subroutine setup_constraint_matrix
  
  ! -----------------------------------------------------------------------------------------------
  subroutine setup_nonlinear_solver (this)
    implicit none
    class(test_hts_nedelec_driver_t), intent(inout) :: this

    call this%nonlinear_solver%create( convergence_criteria = this%test_params%get_nonlinear_convergence_criteria() , &
                                       abs_tol = this%test_params%get_absolute_nonlinear_tolerance(),                 &
                                       rel_tol = this%test_params%get_relative_nonlinear_tolerance(),                 &
                                       max_iters = this%test_params%get_max_nonlinear_iterations(),                   &
                                       ideal_iters = this%test_params%get_stepping_parameter(),                       &
                                       fe_affine_operator = this%fe_affine_operator,                                  &
                                       current_dof_values = this%H_current%get_free_dof_values(),                          &
                                       apply_constraint = this%test_params%get_apply_current_density_constraint(),    & 
                                       constraint_vector = this%constraint_vector                                 )
    
    
  end subroutine setup_nonlinear_solver
  
  ! -----------------------------------------------------------------------------------------------
  subroutine setup_theta_method(this) 
  implicit none 
  class(test_hts_nedelec_driver_t), intent(inout) :: this
  
    call this%theta_method%create( this%test_params%get_theta_value(),          &               
                                   this%test_params%get_initial_time(),         &
                                   this%test_params%get_final_time(),           & 
                                   this%test_params%get_num_time_steps(),    &
                                   this%test_params%get_max_time_step(),        & 
                                   this%test_params%get_min_time_step(),        &
                                   this%test_params%get_save_solution_n_steps() )
  
  end subroutine setup_theta_method 

  ! -----------------------------------------------------------------------------------------------
  subroutine assemble_system (this)
    implicit none
    class(test_hts_nedelec_driver_t), intent(inout) :: this
    class(matrix_t), pointer       :: matrix
    class(vector_t), pointer       :: rhs
    call this%fe_affine_operator%numerical_setup()
    !rhs    => this%fe_affine_operator%get_translation()
    !matrix => this%fe_affine_operator%get_matrix()
    !select type(matrix)
    !class is (sparse_matrix_t)  
    !   call matrix%print_matrix_market(6) 
    !class DEFAULT
    !   assert(.false.) 
    !end select
  end subroutine assemble_system
  
  ! -----------------------------------------------------------------------------------------------
  subroutine solve_system(this) 
  implicit none 
  class(test_hts_nedelec_driver_t), intent(inout) :: this
    
  temporal: do while ( .not. this%theta_method%finished() ) 
     call this%theta_method%print(6) 
    
     call this%solve_nonlinear_system()

     if (this%nonlinear_solver%converged() ) then  ! Theta method goes forward 
        call this%compute_hysteresis_data() 
        call this%theta_method%update_solutions(this%H_current, this%H_previous)
								call this%fe_space%update_hanging_dof_values(this%H_current)
        call this%write_time_step_solution() 
        call this%theta_method%move_time_forward( this%nonlinear_solver%get_current_iteration(), &
                                                  this%nonlinear_solver%get_ideal_num_iterations() ) 
     elseif (.not. this%nonlinear_solver%converged()) then ! Theta method goes backwards and restarts   
        call this%theta_method%move_time_backwards(this%H_current, this%H_previous)
     end if

     if (.not. this%theta_method%finished() ) then 
					   call this%fe_space%interpolate_dirichlet_values(this%H_current, time=this%theta_method%get_current_time()  )
        call this%assemble_system() 
     end if

  end do temporal
 
  end subroutine solve_system 
 
  ! -----------------------------------------------------------------------------------------------
   subroutine solve_nonlinear_system(this)
    implicit none
    class(test_hts_nedelec_driver_t), intent(inout) :: this
    type(constraint_value_t), pointer :: constraint_value_function 
    real(rp)                          :: constraint_value 
    
    ! Get constraint value for the system 
    if (this%nonlinear_solver%get_apply_current_constraint() ) then 
    constraint_value_function => this%problem_functions%get_constraint_value() 
    call constraint_value_function%get_constraint_value(this%theta_method%get_current_time(), constraint_value)
    call this%nonlinear_solver%initialize(constraint_value) 
    else 
    call this%nonlinear_solver%initialize() 
    end if 
    

    nonlinear: do while ( .not. this%nonlinear_solver%finished() )
    ! 0 - Update initial residual 
    call this%nonlinear_solver%start_new_iteration() 
    ! 1 - Integrate Jacobian
    call this%nonlinear_solver%compute_jacobian(this%hts_nedelec_integration)
    ! 2 - Solve tangent system 
    if (this%nonlinear_solver%get_apply_current_constraint() ) then 
    call this%nonlinear_solver%solve_constrained_tangent_system( this%constraint_matrix ) 
    else 
    call this%nonlinear_solver%solve_tangent_system() 
    end if 
    ! 3 - Determine step length to update solution 
    call this%nonlinear_solver%line_search%cubic_backtracking( this%nonlinear_solver )
    ! 4 - Update solution 
    call this%nonlinear_solver%update_solution() 
    ! 5 - New picard iterate with updated solution 
    call this%assemble_system()  
    ! 6 - Evaluate new residual 
    call this%nonlinear_solver%compute_residual()
    ! 7 - Print current output 
    call this%nonlinear_solver%print_current_iteration_output() 
   
    end do nonlinear 
    
    call this%nonlinear_solver%print_final_output() 

    
  end subroutine solve_nonlinear_system
  
  ! -----------------------------------------------------------------------------------------------
  subroutine compute_hysteresis_data(this)
    implicit none 
    class(test_hts_nedelec_driver_t)   , intent(inout) :: this
    class(vector_t),      pointer                      :: dof_values_current_solution     
    class(fe_cell_iterator_t), allocatable :: fe
    ! Integration loop 
    type(quadrature_t)       , pointer     :: quad
    type(fe_cell_function_vector_t)        :: fe_cell_function_current
    integer(ip)                            :: qpoin, num_quad_points, idof 
    type(point_t)            , pointer     :: quad_coords(:)
    type(point_t)         , allocatable    :: aux_quad_coords(:)
    integer(ip)                            :: inode, num_nodes 
    real(rp)                               :: factor 
    type(vector_field_t)                   :: H_value, H_curl 
    ! Hysteresis variables for final computations 
				real(rp)                               :: resistivity 
    real(rp)                               :: Hy_average, xJ_average
    real(rp)                               :: hts_volume
    real(rp)                               :: Happ, AC_loss 
    real(rp)                               :: hts_domain_length(0:SPACE_DIM-1)
    class(scalar_function_t) , pointer     :: boundary_function_Hy
    type(constraint_value_t) , pointer     :: constraint_value_function 
    real(rp) :: constraint_value 
    
    integer(ip) :: istat 

    ! Integrate structures needed 
    call fe_cell_function_current%create(this%fe_space,  1)
    call this%fe_space%set_up_cell_integration()
    call this%fe_space%create_fe_cell_iterator(fe)
				call fe%update_integration() 
    quad             => fe%get_quadrature()
    num_quad_points  = quad%get_num_quadrature_points()
    quad_coords      => fe%get_quadrature_points_coordinates()
    aux_quad_coords  = quad_coords

    ! Loop over elements
    Hy_average  = 0
    xJ_average  = 0
    hts_volume  = 0.0_rp 
				AC_loss     = 0.0_rp 
    do while ( .not. fe%has_finished())

       if ( fe%get_set_id() == hts ) then  ! Integrate only in HTS device DOMAIN 
          ! Update FE-integration related data structures
          call fe%update_integration()
          call fe_cell_function_current%update(fe, this%H_current)

          ! Get quadrature coordinates to evaluate boundary value
          quad_coords => fe%get_quadrature_points_coordinates()

          ! Integrate cell contribution to H_y, x·J_z average 
          do qpoin=1, num_quad_points
             factor = fe%get_det_jacobian(qpoin) * quad%get_weight(qpoin) 		
             call fe_cell_function_current%get_value(qpoin, H_value)
             call fe_cell_function_current%compute_curl(qpoin, H_curl)
													resistivity  = this%hts_nedelec_integration%compute_resistivity( H_curl, HTS )
													
             Hy_average  = Hy_average  + factor*H_value%get(2)          
             xJ_average  = xJ_average  + factor*quad_coords(qpoin)%get(1)*H_curl%get(3)   
													AC_loss     = AC_loss     + factor*resistivity*H_curl%get(3)*H_curl%get(3) 
          end do
           hts_volume = hts_volume + fe%compute_volume()
       end if
       call fe%next()
    end do
    call this%fe_space%free_fe_cell_iterator(fe)

    ! Coordinates of quadrature does influence the constant value Happ(t) 
    boundary_function_Hy => this%problem_functions%get_boundary_function_Hy()
    call boundary_function_Hy%get_value_space_time( aux_quad_coords(1), this%theta_method%get_current_time() , Happ )
    constraint_value_function => this%problem_functions%get_constraint_value() 
    call constraint_value_function%get_constraint_value(this%theta_method%get_current_time(), constraint_value)
    write(*,*) 'Hysteresis Data -----------------------------------------'
    write(*,*) 'mu0·(Hy-Happ)', this%test_params%get_air_permeability()*(Hy_average/hts_volume-Happ), 'Happ', Happ 
    write(*,*) 'xJ', xJ_average/hts_volume, 'AC_loss', AC_loss  
    write(*,*) ' --------------------------------------------------------' 

  end subroutine compute_hysteresis_data
  
  ! -----------------------------------------------------------------------------------------------
  subroutine check_solution(this)
    implicit none
    class(test_hts_nedelec_driver_t), intent(inout) :: this
    class(vector_function_t), pointer :: H_exact_function
    
    type(error_norms_vector_t)  :: H_error_norm
    type(error_norms_scalar_t)  :: p_error_norm 
    real(rp) :: l2, hcurl, l2p
    
    H_exact_function => this%problem_functions%get_solution()
    
    call H_error_norm%create(this%fe_space,1)

    l2 = H_error_norm%compute(H_exact_function, this%H_current, l2_norm, time=this%theta_method%get_current_time() - this%theta_method%get_time_step() )   
    hcurl = H_error_norm%compute(H_exact_function, this%H_current, hcurl_seminorm, time=this%theta_method%get_current_time() - this%theta_method%get_time_step() )    
    
    write(*,*) 'H ERROR NORMS **********************' 
    write(*,'(a20,f20.16)') 'l2_norm(H):', l2;        
    write(*,'(a20,f20.16)') 'hcurl_norm(H):', hcurl;      
    
    call H_error_norm%free()
  end subroutine check_solution 
  
    ! -----------------------------------------------------------------------------------------------
  subroutine initialize_output(this)
    implicit none
    class(test_hts_nedelec_driver_t), intent(inout) :: this
    integer(ip)                                      :: err
    if(this%test_params%get_write_solution()) then
       call  this%oh%create(VTK) 
       call  this%oh%attach_fe_space(this%fe_space)
       call  this%oh%add_fe_function(this%H_current, 1, 'H')
       call  this%oh%add_fe_function(this%H_current, 1, 'J',       curl_diff_operator)
       call  this%oh%add_fe_function(this%H_current, 1, 'div(H)',  div_diff_operator)
       call  this%oh%open(this%test_params%get_dir_path_out(), this%test_params%get_prefix())
    endif
  end subroutine initialize_output
  
  ! -----------------------------------------------------------------------------------------------
  subroutine write_time_step_solution(this)
    implicit none
    class(test_hts_nedelec_driver_t), intent(inout) :: this
    integer(ip)                                      :: err
    
    ! Define customized output by the user 
    if( this%test_params%get_write_solution() .and. this%theta_method%print_this_step() ) then
        call this%oh%append_time_step(this%theta_method%get_current_time())
        call this%oh%write()
        call this%theta_method%update_time_to_be_printed() 
    endif
   
  end subroutine write_time_step_solution
  
      ! -----------------------------------------------------------------------------------------------
  subroutine finalize_output(this)
    implicit none
    class(test_hts_nedelec_driver_t), intent(inout) :: this
    integer(ip)                                      :: err
    if(this%test_params%get_write_solution()) then
    call this%oh%close()
    call this%oh%free()
    endif
  end subroutine finalize_output

  ! -----------------------------------------------------------------------------------------------
  subroutine run_simulation(this) 
    implicit none
    class(test_hts_nedelec_driver_t), intent(inout) :: this
    call this%free()
    call this%parse_command_line_parameters()
    call this%setup_triangulation()
    call this%setup_reference_fes()
    call this%setup_theta_method() 
    call this%setup_fe_space()
    call this%setup_system()
    call this%assemble_system()   
    call this%setup_nonlinear_solver()
    call this%initialize_output() 
    call this%solve_system()
    call this%finalize_output()
   ! call this%check_solution() 
    call this%free()
  end subroutine run_simulation

  ! -----------------------------------------------------------------------------------------------
  subroutine free(this)
    implicit none
    class(test_hts_nedelec_driver_t), intent(inout) :: this
    integer(ip) :: i, istat
    call this%H_previous%free()
    call this%H_current%free() 
    call this%fe_affine_operator%free()
    call this%fe_space%free()
    if ( allocated(this%reference_fes) ) then
       do i=1, size(this%reference_fes)
          call this%reference_fes(i)%p%free()
       end do
       deallocate(this%reference_fes, stat=istat)
       check(istat==0)
    end if
    call this%triangulation%free()
    call this%test_params%free()
    if ( this%nonlinear_solver%get_apply_current_constraint() ) call this%constraint_matrix%free()
    if ( this%nonlinear_solver%get_apply_current_constraint() ) call this%constraint_vector%free()
    call this%nonlinear_solver%free()
  end subroutine free

end module test_hts_nedelec_driver_names
