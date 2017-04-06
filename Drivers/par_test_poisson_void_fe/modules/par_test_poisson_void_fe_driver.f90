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
module par_test_poisson_void_fe_driver_names
  use fempar_names
  use list_types_names
  use par_poisson_static_parameters_names
  use par_test_poisson_void_fe_params_names
  use poisson_void_fe_discrete_integration_names
  use poisson_void_fe_conditions_names
  use poisson_void_fe_analytical_functions_names
# include "debug.i90"

  implicit none
  private

  type par_test_poisson_void_fe_fe_driver_t 
     private 
     
     ! Place-holder for parameter-value set provided through command-line interface
     type(par_test_poisson_void_fe_params_t)      :: test_params
     type(ParameterList_t), pointer       :: parameter_list
     
     ! Cells and lower dimension objects container
     type(par_triangulation_t)             :: triangulation
     integer(ip), allocatable                  :: cell_set_ids(:)
     
     ! Discrete weak problem integration-related data type instances 
     type(par_fe_space_t)                      :: fe_space 
     type(p_reference_fe_t), allocatable       :: reference_fes(:) 
     type(standard_l1_coarse_fe_handler_t)     :: l1_coarse_fe_handler
     type(poisson_void_fe_CG_discrete_integration_t)   :: poisson_void_fe_integration
     type(poisson_void_fe_conditions_t)                :: poisson_void_fe_conditions
     type(poisson_void_fe_analytical_functions_t)      :: poisson_void_fe_analytical_functions

     
     ! Place-holder for the coefficient matrix and RHS of the linear system
     type(fe_affine_operator_t)            :: fe_affine_operator
     
!#ifdef ENABLE_MKL     
     ! MLBDDC preconditioner
     type(mlbddc_t)                            :: mlbddc
!#endif  
    
     ! Iterative linear solvers data type
     type(iterative_linear_solver_t)           :: iterative_linear_solver
 
     ! poisson_void_fe problem solution FE function
     type(fe_function_t)                   :: solution
     
     ! Environment required for fe_affine_operator + vtk_handler
     !type(par_context_t)                       :: w_context
     type(environment_t), pointer           :: par_environment
   contains
     procedure                  :: run_simulation
     procedure        , private :: parse_command_line_parameters
     procedure        , private :: setup_triangulation
     procedure        , private :: setup_reference_fes
     procedure        , private :: setup_fe_space
     procedure        , private :: setup_system
     procedure        , private :: setup_solver
     procedure        , private :: assemble_system
     procedure        , private :: solve_system
     procedure        , private :: check_solution
     procedure        , private :: write_solution
     procedure        , private :: free
     procedure, nopass, private :: ls_fun => par_test_poisson_void_fe_driver_ls_fun
  end type par_test_poisson_void_fe_fe_driver_t

  ! Types
  public :: par_test_poisson_void_fe_fe_driver_t

contains

  subroutine parse_command_line_parameters(this)
    implicit none
    class(par_test_poisson_void_fe_fe_driver_t), intent(inout) :: this
    call this%test_params%create()
    !call this%test_params%parse(this%parameter_list)
    this%parameter_list => this%test_params%get_values()
  end subroutine parse_command_line_parameters
   
  subroutine setup_triangulation(this)
    implicit none
    class(par_test_poisson_void_fe_fe_driver_t), intent(inout) :: this

    type(cell_iterator_t) :: cell_iter
    type(cell_accessor_t) :: cell
    type(point_t), allocatable :: cell_coords(:)
    integer(ip) :: istat
    integer(ip) :: set_id
    real(rp) :: x, y
    integer(ip) :: num_void_neigs

    integer(ip)           :: ivef
    type(vef_iterator_t)  :: vef_iterator
    type(vef_accessor_t)  :: vef, vef_of_vef
    type(list_t), pointer :: vefs_of_vef
    type(list_t), pointer :: vertices_of_line
    type(list_iterator_t) :: vefs_of_vef_iterator
    type(list_iterator_t) :: vertices_of_line_iterator
    class(lagrangian_reference_fe_t), pointer :: reference_fe_geo
    integer(ip) :: ivef_pos_in_cell, vef_of_vef_pos_in_cell
    integer(ip) :: vertex_pos_in_cell, icell_arround
    integer(ip) :: inode, num

    call this%triangulation%create(this%parameter_list)
    this%par_environment => this%triangulation%get_par_environment()

    ! Set the cell ids
    if ( this%par_environment%am_i_l1_task() ) then
      call memalloc(this%triangulation%get_num_local_cells(),this%cell_set_ids)
      cell_iter = this%triangulation%create_cell_iterator()
      call cell_iter%current(cell)
      allocate(cell_coords(1:cell%get_num_nodes()),stat=istat); check(istat == 0)
      do while( .not. cell_iter%has_finished() )
        call cell_iter%current(cell)
        if (cell%is_local()) then
          set_id = PAR_POISSON_SET_ID_VOID
          call cell%get_coordinates(cell_coords)
          do inode = 1,cell%get_num_nodes()
            if ( this%ls_fun(cell_coords(inode),&
              this%triangulation%get_num_dimensions()) < 0.0 ) then
              set_id = PAR_POISSON_SET_ID_FULL
              exit
            end if
          end do
          this%cell_set_ids(cell%get_lid()) = set_id
        end if
        call cell_iter%next()
      end do
      deallocate(cell_coords, stat = istat); check(istat == 0)
      call this%triangulation%fill_cells_set(this%cell_set_ids)
    end if

    ! Initialize all the vefs set ids to SET_ID_FREE
    vef_iterator = this%triangulation%create_vef_iterator()
    do while ( .not. vef_iterator%has_finished() )
      call vef_iterator%current(vef)
      call vef%set_set_id(PAR_POISSON_SET_ID_FREE)
      call vef_iterator%next()
    end do

    ! Fix all the vefs at the boundary
    if ( this%test_params%get_triangulation_type() == triangulation_generate_structured ) then
      vef_iterator = this%triangulation%create_vef_iterator()
      do while ( .not. vef_iterator%has_finished() )
        call vef_iterator%current(vef)
        if(vef%is_at_boundary()) call vef%set_set_id(PAR_POISSON_SET_ID_DIRI)
        call vef_iterator%next()
      end do
    end if

    ! Set all the vefs on the interface between full/void
    vef_iterator = this%triangulation%create_vef_iterator()
    do while ( .not. vef_iterator%has_finished() )
       call vef_iterator%current(vef)
       ! If it is an INTERIOR face
       if( vef%get_dimension() == this%triangulation%get_num_dimensions()-1 .and. vef%get_num_cells_around()==2 ) then

         ! Compute number of void neighbors
         num_void_neigs = 0
         do icell_arround = 1,vef%get_num_cells_around()
           call vef%get_cell_around(icell_arround,cell)
           if (cell%get_set_id() == PAR_POISSON_SET_ID_VOID) num_void_neigs = num_void_neigs + 1
         end do

         if(num_void_neigs==1) then ! If vef (face) is between a full and a void cell

             ! Set this face as Dirichlet boundary
             call vef%set_set_id(PAR_POISSON_SET_ID_DIRI)

             ! Do a loop on all edges in 3D (vertex in 2D) of the face
             ivef = vef%get_lid()
             call vef%get_cell_around(1,cell) ! There is always one cell around
             reference_fe_geo => cell%get_reference_fe_geo()
             ivef_pos_in_cell = cell%find_lpos_vef_lid(ivef)
             vefs_of_vef => reference_fe_geo%get_n_faces_n_face()
             vefs_of_vef_iterator = vefs_of_vef%create_iterator(ivef_pos_in_cell)
             do while( .not. vefs_of_vef_iterator%is_upper_bound() )

                ! Set edge (resp. vertex) as Dirichlet
                vef_of_vef_pos_in_cell = vefs_of_vef_iterator%get_current()
                call cell%get_vef(vef_of_vef_pos_in_cell, vef_of_vef)
                call vef_of_vef%set_set_id(PAR_POISSON_SET_ID_DIRI)

                ! If 3D, traverse vertices of current line
                if ( this%triangulation%get_num_dimensions() == 3 ) then
                  vertices_of_line          => reference_fe_geo%get_vertices_n_face()
                  vertices_of_line_iterator = vertices_of_line%create_iterator(vef_of_vef_pos_in_cell)
                  do while( .not. vertices_of_line_iterator%is_upper_bound() )

                    ! Set vertex as Dirichlet
                    vertex_pos_in_cell = vertices_of_line_iterator%get_current()
                    call cell%get_vef(vertex_pos_in_cell, vef_of_vef)
                    call vef_of_vef%set_set_id(PAR_POISSON_SET_ID_DIRI)

                    call vertices_of_line_iterator%next()
                  end do ! Loop in vertices in 3D only
                end if

                call vefs_of_vef_iterator%next()
             end do ! Loop in edges (resp. vertices)

         end if ! If face on void/full boundary
       end if ! If vef is an interior face

       call vef_iterator%next()
    end do ! Loop in vefs

    call this%triangulation%setup_coarse_triangulation()
    
  end subroutine setup_triangulation
  
  subroutine setup_reference_fes(this)
    implicit none
    class(par_test_poisson_void_fe_fe_driver_t), intent(inout) :: this
    integer(ip) :: istat
    type(cell_iterator_t)                     :: cell_iterator
    type(cell_accessor_t)                     :: cell
    class(lagrangian_reference_fe_t), pointer :: reference_fe_geo
    
    allocate(this%reference_fes(2), stat=istat)
    check(istat==0)
    
    if ( this%par_environment%am_i_l1_task() ) then
      cell_iterator = this%triangulation%create_cell_iterator()
      call cell_iterator%current(cell)
      reference_fe_geo => cell%get_reference_fe_geo()
      this%reference_fes(PAR_POISSON_SET_ID_FULL) =  make_reference_fe ( topology = reference_fe_geo%get_topology(), &
                                                   fe_type = fe_type_lagrangian, &
                                                   number_dimensions = this%triangulation%get_num_dimensions(), &
                                                   order = this%test_params%get_reference_fe_order(), &
                                                   field_type = field_type_scalar, &
                                                   continuity = .true. )

      this%reference_fes(PAR_POISSON_SET_ID_VOID) =  make_reference_fe ( topology = reference_fe_geo%get_topology(), &
                                                   fe_type = fe_type_void, &
                                                   number_dimensions = this%triangulation%get_num_dimensions(), &
                                                   order = -1, &
                                                   field_type = field_type_scalar, &
                                                   continuity = .true. )
    end if  
  end subroutine setup_reference_fes

  subroutine setup_fe_space(this)
    implicit none
    class(par_test_poisson_void_fe_fe_driver_t), intent(inout) :: this

    integer(ip) :: set_ids_to_reference_fes(1,2)

    set_ids_to_reference_fes(1,PAR_POISSON_SET_ID_FULL) = PAR_POISSON_SET_ID_FULL
    set_ids_to_reference_fes(1,PAR_POISSON_SET_ID_VOID) = PAR_POISSON_SET_ID_VOID

    
    call this%fe_space%create( triangulation            = this%triangulation, &
                               conditions               = this%poisson_void_fe_conditions, &
                               reference_fes            = this%reference_fes, &
                               set_ids_to_reference_fes = set_ids_to_reference_fes, &
                               coarse_fe_handler        = this%l1_coarse_fe_handler)
    
    call this%fe_space%fill_dof_info() 
    call this%fe_space%setup_coarse_fe_space(this%parameter_list)
    call this%fe_space%initialize_fe_integration()
    
    call this%poisson_void_fe_analytical_functions%set_num_dimensions(this%triangulation%get_num_dimensions())
    call this%poisson_void_fe_conditions%set_boundary_function(this%poisson_void_fe_analytical_functions%get_boundary_function())
    call this%fe_space%interpolate_dirichlet_values(this%poisson_void_fe_conditions)    
    !call this%fe_space%print()
  end subroutine setup_fe_space
  
  subroutine setup_system (this)
    implicit none
    class(par_test_poisson_void_fe_fe_driver_t), intent(inout) :: this
    
    call this%poisson_void_fe_integration%set_analytical_functions(this%poisson_void_fe_analytical_functions)
    
    ! if (test_single_scalar_valued_reference_fe) then
    call this%fe_affine_operator%create ( sparse_matrix_storage_format      = csr_format, &
                                          diagonal_blocks_symmetric_storage = [ .true. ], &
                                          diagonal_blocks_symmetric         = [ .true. ], &
                                          diagonal_blocks_sign              = [ SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE ], &
                                          fe_space                          = this%fe_space, &
                                          discrete_integration              = this%poisson_void_fe_integration )
  end subroutine setup_system
  
  subroutine setup_solver (this)
    implicit none
    class(par_test_poisson_void_fe_fe_driver_t), intent(inout) :: this
    type(parameterlist_t) :: parameter_list
    integer(ip) :: FPLError

!#ifdef ENABLE_MKL   
    ! Set-up MLBDDC preconditioner
    call this%mlbddc%create(this%fe_affine_operator, this%parameter_list)
    call this%mlbddc%symbolic_setup()
    call this%mlbddc%numerical_setup()
!#endif    
   
    call this%iterative_linear_solver%create(this%fe_space%get_environment())
    call this%iterative_linear_solver%set_type_from_string(cg_name)

!#ifdef ENABLE_MKL
    call this%iterative_linear_solver%set_operators(this%fe_affine_operator, this%mlbddc) 
!#else
!    call parameter_list%init()
!    FPLError = parameter_list%set(key = ils_rtol, value = 1.0e-12_rp)
!    FPLError = parameter_list%set(key = ils_max_num_iterations, value = 5000)
!    assert(FPLError == 0)
!    call this%iterative_linear_solver%set_parameters_from_pl(parameter_list)
!    call this%iterative_linear_solver%set_operators(this%fe_affine_operator, .identity. this%fe_affine_operator) 
!    call parameter_list%free()
!#endif   
    
  end subroutine setup_solver
  
  
  subroutine assemble_system (this)
    implicit none
    class(par_test_poisson_void_fe_fe_driver_t), intent(inout) :: this
    class(matrix_t)                  , pointer       :: matrix
    class(vector_t)                  , pointer       :: rhs
    call this%fe_affine_operator%numerical_setup()
    rhs                => this%fe_affine_operator%get_translation()
    matrix             => this%fe_affine_operator%get_matrix()
    
    !select type(matrix)
    !class is (sparse_matrix_t)  
    !   call matrix%print_matrix_market(6) 
    !class DEFAULT
    !   assert(.false.) 
    !end select
    
    !select type(rhs)
    !class is (serial_scalar_array_t)  
    !   call rhs%print(6) 
    !class DEFAULT
    !   assert(.false.) 
    !end select
  end subroutine assemble_system
  
  
  subroutine solve_system(this)
    implicit none
    class(par_test_poisson_void_fe_fe_driver_t), intent(inout) :: this
    class(matrix_t)                         , pointer       :: matrix
    class(vector_t)                         , pointer       :: rhs
    class(vector_t)                         , pointer       :: dof_values

    matrix     => this%fe_affine_operator%get_matrix()
    rhs        => this%fe_affine_operator%get_translation()
    dof_values => this%solution%get_dof_values()
    call this%iterative_linear_solver%solve(this%fe_affine_operator%get_translation(), &
                                            dof_values)
    
    !select type (dof_values)
    !class is (serial_scalar_array_t)  
    !   call dof_values%print(6)
    !class DEFAULT
    !   assert(.false.) 
    !end select
    
    !select type (matrix)
    !class is (sparse_matrix_t)  
    !   call this%direct_solver%update_matrix(matrix, same_nonzero_pattern=.true.)
    !   call this%direct_solver%solve(rhs , dof_values )
    !class DEFAULT
    !   assert(.false.) 
    !end select
  end subroutine solve_system
   
  subroutine check_solution(this)
    implicit none
    class(par_test_poisson_void_fe_fe_driver_t), intent(inout) :: this
    type(error_norms_scalar_t) :: error_norm 
    real(rp) :: mean, l1, l2, lp, linfty, h1, h1_s, w1p_s, w1p, w1infty_s, w1infty
    
    call error_norm%create(this%fe_space,1)    
    mean = error_norm%compute(this%poisson_void_fe_analytical_functions%get_solution_function(), this%solution, mean_norm)   
    l1 = error_norm%compute(this%poisson_void_fe_analytical_functions%get_solution_function(), this%solution, l1_norm)   
    l2 = error_norm%compute(this%poisson_void_fe_analytical_functions%get_solution_function(), this%solution, l2_norm)   
    lp = error_norm%compute(this%poisson_void_fe_analytical_functions%get_solution_function(), this%solution, lp_norm)   
    linfty = error_norm%compute(this%poisson_void_fe_analytical_functions%get_solution_function(), this%solution, linfty_norm)   
    h1_s = error_norm%compute(this%poisson_void_fe_analytical_functions%get_solution_function(), this%solution, h1_seminorm) 
    h1 = error_norm%compute(this%poisson_void_fe_analytical_functions%get_solution_function(), this%solution, h1_norm) 
    w1p_s = error_norm%compute(this%poisson_void_fe_analytical_functions%get_solution_function(), this%solution, w1p_seminorm)   
    w1p = error_norm%compute(this%poisson_void_fe_analytical_functions%get_solution_function(), this%solution, w1p_norm)   
    w1infty_s = error_norm%compute(this%poisson_void_fe_analytical_functions%get_solution_function(), this%solution, w1infty_seminorm) 
    w1infty = error_norm%compute(this%poisson_void_fe_analytical_functions%get_solution_function(), this%solution, w1infty_norm)  
    if ( this%par_environment%am_i_l1_root() ) then
      write(*,'(a20,e32.25)') 'mean_norm:', mean; check ( abs(mean) < 1.0e-04 )
      write(*,'(a20,e32.25)') 'l1_norm:', l1; check ( l1 < 1.0e-04 )
      write(*,'(a20,e32.25)') 'l2_norm:', l2; check ( l2 < 1.0e-04 )
      write(*,'(a20,e32.25)') 'lp_norm:', lp; check ( lp < 1.0e-04 )
      write(*,'(a20,e32.25)') 'linfnty_norm:', linfty; check ( linfty < 1.0e-04 )
      write(*,'(a20,e32.25)') 'h1_seminorm:', h1_s; check ( h1_s < 1.0e-04 )
      write(*,'(a20,e32.25)') 'h1_norm:', h1; check ( h1 < 1.0e-04 )
      write(*,'(a20,e32.25)') 'w1p_seminorm:', w1p_s; check ( w1p_s < 1.0e-04 )
      write(*,'(a20,e32.25)') 'w1p_norm:', w1p; check ( w1p < 1.0e-04 )
      write(*,'(a20,e32.25)') 'w1infty_seminorm:', w1infty_s; check ( w1infty_s < 1.0e-04 )
      write(*,'(a20,e32.25)') 'w1infty_norm:', w1infty; check ( w1infty < 1.0e-04 )
    end if  
    call error_norm%free()
  end subroutine check_solution
  
  subroutine write_solution(this)
    implicit none
    class(par_test_poisson_void_fe_fe_driver_t), intent(in) :: this

    type(output_handler_t)                          :: oh
    real(rp),allocatable :: cell_vector(:)
    real(rp),allocatable :: mypart_vector(:)
    integer(ip) :: istat

    if(this%test_params%get_write_solution() .and. &
       this%par_environment%am_i_l1_task()) then

        allocate(cell_vector(1:size(this%cell_set_ids)),stat=istat ); check(istat == 0)
        cell_vector(:) = this%cell_set_ids(:)

        allocate(mypart_vector(1:size(this%cell_set_ids)),stat=istat ); check(istat == 0)
        mypart_vector(:) = this%par_environment%get_l1_rank()

        call oh%create()
        call oh%attach_fe_space(this%fe_space)
        call oh%add_fe_function(this%solution, 1, 'solution')
        call oh%add_cell_vector(cell_vector,'cell_set_ids')
        call oh%add_cell_vector(mypart_vector,'l1_rank')
        call oh%open(this%test_params%get_dir_path(), this%test_params%get_prefix())
        call oh%write()
        call oh%close()
        call oh%free()

        deallocate(cell_vector,stat=istat); check(istat == 0)
        deallocate(mypart_vector,stat=istat); check(istat == 0)

    endif

  end subroutine write_solution
  
  subroutine run_simulation(this) 
    implicit none
    class(par_test_poisson_void_fe_fe_driver_t), intent(inout) :: this
    !call this%free()
    call this%parse_command_line_parameters()
    !call this%setup_context()
    !call this%setup_par_environment()
    call this%setup_triangulation()
    call this%setup_reference_fes()
    call this%setup_fe_space()
    call this%setup_system()
    call this%assemble_system()
    call this%setup_solver()
    call this%solution%create(this%fe_space) 
    call this%solve_system()
    call this%check_solution()
    call this%write_solution()
    call this%free()
  end subroutine run_simulation
  
  subroutine free(this)
    implicit none
    class(par_test_poisson_void_fe_fe_driver_t), intent(inout) :: this
    integer(ip) :: i, istat
    
    call this%solution%free()
!#ifdef ENABLE_MKL    
    call this%mlbddc%free()
!#endif    
    call this%iterative_linear_solver%free()
    call this%fe_affine_operator%free()
    call this%fe_space%free()
    if ( allocated(this%reference_fes) ) then
      do i=1, size(this%reference_fes)
        call this%reference_fes(i)%free()
      end do
      deallocate(this%reference_fes, stat=istat)
      check(istat==0)
    end if
    call this%triangulation%free()
    call this%test_params%free()
    if (allocated(this%cell_set_ids)) call memfree(this%cell_set_ids,__FILE__,__LINE__)
    call this%par_environment%free() 
    !call this%w_context%free(.true.)
  end subroutine free  

  function par_test_poisson_void_fe_driver_ls_fun(point,num_dim) result (val)
    implicit none
    type(point_t), intent(in) :: point
    integer(ip),   intent(in) :: num_dim
    real(rp) :: val
    type(point_t) :: p
    real(rp) :: x, y, z
    real(rp) :: xk, yk, zk
    real(rp) :: r0, sg, A
    integer(ip) :: k

    p = point
    if (num_dim < 3) call p%set(3,0.62)
    !x = 1.85*( 2.0*p%get(1) - 1.0 )
    !y = 1.85*( 2.0*p%get(2) - 1.0 )
    !z = 1.85*( 2.0*p%get(3) - 1.0 )
    !val = (x**2+y**2-4)**2 + (z**2-1.2)**2 + (y**2+z**2-4)**2 +&
    !      (x**2-1.2)**2 + (z**2+x**2-4)**2 + (y**2-1.2)**2 - 12
    x = ( 2.0*p%get(1) - 1.0 )
    y = ( 2.0*p%get(2) - 1.0 )
    z = ( 2.0*p%get(3) - 1.0 )
    r0 = 0.6
    sg = 0.2
    A  = 2.0
    val = sqrt(x**2 + y**2 + z**2) - r0
    do k = 0,11
        if (0 <= k .and. k <= 4) then
            xk = (r0/sqrt(5.0))*2.0*cos(2.0*k*pi/5.0)
            yk = (r0/sqrt(5.0))*2.0*sin(2.0*k*pi/5.0)
            zk = (r0/sqrt(5.0))
        else if (5 <= k .and. k <= 9) then
            xk = (r0/sqrt(5.0))*2.0*cos((2.0*(k-5)-1.0)*pi/5.0)
            yk = (r0/sqrt(5.0))*2.0*sin((2.0*(k-5)-1.0)*pi/5.0)
            zk =-(r0/sqrt(5.0))
        else if (k == 10) then
            xk = 0
            yk = 0
            zk = r0
        else
            xk = 0
            yk = 0
            zk = -r0
        end if
        val = val - A*exp( -( (x - xk)**2  + (y - yk)**2 + (z - zk)**2 )/(sg**2) )
    end do
  end function par_test_poisson_void_fe_driver_ls_fun
  
end module par_test_poisson_void_fe_driver_names
