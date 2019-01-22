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
module test_h_adaptive_poisson_driver_names
  use fempar_names
  use test_poisson_params_names
  use poisson_cG_discrete_integration_names
  use poisson_dG_discrete_integration_names
  use poisson_conditions_names
  use poisson_analytical_functions_names
  use vector_poisson_discrete_integration_names
  use vector_poisson_conditions_names
  use vector_poisson_analytical_functions_names
  use IR_Precision ! VTK_IO
  use Lib_VTK_IO ! VTK_IO
    
# include "debug.i90"

  implicit none
  private
  
  integer(ip), parameter :: TEST_POISSON_FULL = 1
  integer(ip), parameter :: TEST_POISSON_VOID = 2
  
  type test_h_adaptive_poisson_driver_t 
     private 
     
     ! Place-holder for parameter-value set provided through command-line interface
     type(test_poisson_params_t)   :: test_params
     type(ParameterList_t)         :: parameter_list
     
     ! Cells and lower dim objects container
     type(p4est_serial_triangulation_t)           :: triangulation
     
     ! Discrete weak problem integration-related data type instances 
     type(serial_fe_space_t)                      :: fe_space 
     type(p_reference_fe_t), allocatable          :: reference_fes(:) 
     
     type(poisson_cG_discrete_integration_t)      :: poisson_cG_integration
     type(poisson_dG_discrete_integration_t)      :: poisson_dG_integration
     type(poisson_conditions_t)                   :: poisson_conditions
     type(poisson_analytical_functions_t)         :: poisson_analytical_functions
     
     type(vector_poisson_discrete_integration_t)  :: vector_poisson_integration
     type(vector_poisson_analytical_functions_t)  :: vector_poisson_analytical_functions
     type(vector_poisson_conditions_t)            :: vector_poisson_conditions
     
     ! Place-holder for the coefficient matrix and RHS of the linear system
     type(fe_affine_operator_t)                   :: fe_affine_operator
     
     ! Direct and Iterative linear solvers data type
     type(environment_t)                          :: serial_environment
#ifdef ENABLE_MKL     
     type(direct_solver_t)                        :: direct_solver
#else     
     type(iterative_linear_solver_t)              :: iterative_linear_solver
#endif     
 
     ! Poisson problem solution FE function
     type(fe_function_t)                          :: solution
   contains
     procedure                  :: run_simulation
     procedure                  :: parse_command_line_parameters
     procedure                  :: setup_environment
     procedure                  :: free_environment
     procedure        , private :: setup_triangulation
     procedure        , private :: set_cells_for_refinement
     procedure        , private :: set_cells_for_coarsening
     procedure        , private :: fill_cells_set
     procedure        , private :: setup_reference_fes
     procedure        , private :: setup_fe_space
     procedure        , private :: refine_and_coarsen
     procedure        , private :: setup_system
     procedure        , private :: setup_solver
     procedure        , private :: assemble_system
     procedure        , private :: solve_system     
     procedure        , private :: check_solution
     procedure        , private :: check_solution_vector
     procedure        , private :: write_solution
     procedure        , private :: write_filling_curve
     procedure        , private :: free
  end type test_h_adaptive_poisson_driver_t

  ! Types
  public :: test_h_adaptive_poisson_driver_t

contains

  subroutine parse_command_line_parameters(this)
    implicit none
    class(test_h_adaptive_poisson_driver_t ), intent(inout) :: this
    call this%test_params%create()
    call this%test_params%parse(this%parameter_list)
  end subroutine parse_command_line_parameters
  
  subroutine setup_environment(this, world_context)
    implicit none
    class(test_h_adaptive_poisson_driver_t ), intent(inout) :: this
    class(execution_context_t)  , intent(in)    :: world_context
    integer(ip) :: ierr
    call this%serial_environment%create(world_context, this%parameter_list)
  end subroutine setup_environment
  
  subroutine free_environment(this)
    implicit none
    class(test_h_adaptive_poisson_driver_t ), intent(inout) :: this
    call this%serial_environment%free()
  end subroutine free_environment
  
  subroutine setup_triangulation(this)
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(inout) :: this

    class(cell_iterator_t), allocatable :: cell
    integer(ip), allocatable :: cell_set_ids(:)
    type(point_t), allocatable :: cell_coords(:)
    integer(ip) :: istat
    integer(ip) :: set_id
    real(rp) :: x, y
    integer(ip) :: num_void_neigs
    
    integer(ip)           :: ivef
    class(vef_iterator_t),allocatable :: vef, vef_of_vef
    type(list_t), pointer :: vefs_of_vef
    type(list_t), pointer :: vertices_of_line
    type(list_iterator_t) :: vefs_of_vef_iterator
    type(list_iterator_t) :: vertices_of_line_iterator
    class(reference_fe_t), pointer :: reference_fe_geo
    integer(ip) :: ivef_pos_in_cell, vef_of_vef_pos_in_cell
    integer(ip) :: vertex_pos_in_cell, icell_arround
    
    integer(ip), parameter :: num_nodes_x_cell = 4
    integer(ip) :: i
    
    call this%triangulation%create(this%serial_environment,this%parameter_list)
    
    if ( .not. this%test_params%get_use_void_fes() ) then
      call this%triangulation%create_vef_iterator(vef)
      do while ( .not. vef%has_finished() )
        if(vef%is_at_boundary()) then
          call vef%set_set_id(1)
        else
          call vef%set_set_id(0)
        end if
        call vef%next()
      end do
      call this%triangulation%free_vef_iterator(vef)
    end if 
   
    do i = 1,2
      call this%set_cells_for_refinement()
      call this%triangulation%refine_and_coarsen()
      call this%triangulation%clear_refinement_and_coarsening_flags()
    end do
    
    if (this%test_params%get_use_void_fes()) then
        call memalloc(this%triangulation%get_num_cells(),cell_set_ids)
        call this%triangulation%create_cell_iterator(cell)
        allocate(cell_coords(1:num_nodes_x_cell),stat=istat); check(istat == 0)
        do while( .not. cell%has_finished() )
          if (cell%is_local()) then
            set_id = TEST_POISSON_VOID
            call cell%get_nodes_coordinates(cell_coords)
            select case (this%test_params%get_use_void_fes_case())
            case ('half')
              y = cell_coords(1)%get(2)
              if (y>=0.5) set_id = TEST_POISSON_FULL
            case ('quarter')
              x = cell_coords(1)%get(1)
              y = cell_coords(1)%get(2)
              if (x>=0.5 .and. y>=0.5) set_id = TEST_POISSON_FULL
            case default
              check(.false.)
            end select
            cell_set_ids(cell%get_gid()) = set_id
          end if
          call cell%next()
        end do
        call this%triangulation%free_cell_iterator(cell)
        deallocate(cell_coords, stat = istat); check(istat == 0)
        call this%triangulation%fill_cells_set(cell_set_ids)
    end if
    
    if ( this%test_params%get_triangulation_type() == 'structured') then

      ! Set all the vefs on the interface between full/void if there are void fes
      if ( this%test_params%get_use_void_fes()) then
        ! WARNING: Consider updating this piece of code for h-adaptivity, because 
        ! the vef set_IDs are not correctly assigned on improper vefs and on proper 
        ! vefs containing them. At the moment, all examples are such that only proper 
        ! vefs are on the interface void/full.
        call this%triangulation%create_vef_iterator(vef)
        call this%triangulation%create_vef_iterator(vef_of_vef)
        call this%triangulation%create_cell_iterator(cell)
        do while ( .not. vef%has_finished() )
                                         
           ! If it is an INTERIOR face
           if( vef%is_facet() .and. vef%get_num_cells_around()==2 ) then

             ! Compute number of void neighbors
             num_void_neigs = 0
             do icell_arround = 1,vef%get_num_cells_around()
               call vef%get_cell_around(icell_arround,cell)
               if (cell%get_set_id() == TEST_POISSON_VOID) num_void_neigs = num_void_neigs + 1
             end do

             if(num_void_neigs==1) then ! If vef (face) is between a full and a void cell

                 ! Set this face as Dirichlet boundary
                 call vef%set_set_id(1)

                 ! Do a loop on all edges in 3D (vertex in 2D) of the face
                 ivef = vef%get_gid()
                 call vef%get_cell_around(1,cell) ! There is always one cell around
                 reference_fe_geo => cell%get_reference_fe()
                 ivef_pos_in_cell = cell%get_vef_lid_from_gid(ivef)
                 vefs_of_vef => reference_fe_geo%get_facets_n_face()
                 vefs_of_vef_iterator = vefs_of_vef%create_iterator(ivef_pos_in_cell)
                 do while( .not. vefs_of_vef_iterator%is_upper_bound() )

                    ! Set edge (resp. vertex) as Dirichlet
                    vef_of_vef_pos_in_cell = vefs_of_vef_iterator%get_current()
                    call cell%get_vef(vef_of_vef_pos_in_cell, vef_of_vef)
                    call vef_of_vef%set_set_id(1)

                    ! If 3D, traverse vertices of current line
                    if ( this%triangulation%get_num_dims() == 3 ) then
                      vertices_of_line          => reference_fe_geo%get_vertices_n_face()
                      vertices_of_line_iterator = vertices_of_line%create_iterator(vef_of_vef_pos_in_cell)
                      do while( .not. vertices_of_line_iterator%is_upper_bound() )

                        ! Set vertex as Dirichlet
                        vertex_pos_in_cell = vertices_of_line_iterator%get_current()
                        call cell%get_vef(vertex_pos_in_cell, vef_of_vef)
                        call vef_of_vef%set_set_id(1)

                        call vertices_of_line_iterator%next()
                      end do ! Loop in vertices in 3D only
                    end if

                    call vefs_of_vef_iterator%next()
                 end do ! Loop in edges (resp. vertices)

             end if ! If face on void/full boundary
           end if ! If vef is an interior face
           
           if(vef%is_at_boundary()) then
              call vef%get_cell_around(1,cell)
              if (cell%get_set_id() == TEST_POISSON_FULL) then
                 call vef%set_set_id(1)
              end if
           end if
           
           call vef%next()
        end do ! Loop in vefs
        call this%triangulation%free_vef_iterator(vef)
        call this%triangulation%free_vef_iterator(vef_of_vef)
        call this%triangulation%free_cell_iterator(cell)
        
      end if
    
    end if
    
    if (allocated(cell_set_ids)) call memfree(cell_set_ids,__FILE__,__LINE__)
    
  end subroutine setup_triangulation
  
  subroutine set_cells_for_refinement(this)
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(inout) :: this
    class(cell_iterator_t)      , allocatable :: cell
    type(point_t), allocatable :: coords(:)
    integer(ip) :: istat, k
    real(rp) ::  x,y
    real(rp), parameter :: Re = 0.46875
    real(rp), parameter :: Ri = 0.15625
    real(rp) :: R
    integer(ip), parameter :: max_num_cell_nodes = 4
    integer(ip), parameter :: max_level = 4

    call this%triangulation%create_cell_iterator(cell)
    if (this%triangulation%get_num_dims() == 2) then
      allocate(coords(max_num_cell_nodes),stat=istat); check(istat==0)

      do while ( .not. cell%has_finished() )
        !if ( mod(cell%get_gid()-1,2) == 0 ) then
        !  call cell%set_for_refinement()
        !end if

        call cell%get_nodes_coordinates(coords)
        x = 0.0
        y = 0.0
        do k=1,max_num_cell_nodes
         x = x + (1.0/max_num_cell_nodes)*coords(k)%get(1)
         y = y + (1.0/max_num_cell_nodes)*coords(k)%get(2)
        end do
        R = sqrt( (x-0.5)**2 + (y-0.5)**2 )
       
        if ( ((R - Re) < 0.0) .and. ((R - Ri) > 0.0) .and. (cell%get_level()<= max_level) .or. (cell%get_level() == 0) )then
          call cell%set_for_refinement()
        end if


        call cell%next()
      end do
      deallocate(coords,stat=istat); check(istat==0)

    else if (this%triangulation%get_num_dims() == 3) then    
      do while ( .not. cell%has_finished() )
          if ( (cell%get_ggid()==1) .or. (cell%get_ggid()==4) .or. (cell%get_ggid()==5) .or. (cell%get_ggid()==8) )then
          call cell%set_for_refinement()
        end if
        call cell%next()
      end do

    else
      mcheck(.false.,'Only for 2D and 3D')

    end if

    call this%triangulation%free_cell_iterator(cell)

  end subroutine set_cells_for_refinement
  
  subroutine set_cells_for_coarsening(this)
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(inout) :: this
    class(cell_iterator_t)      , allocatable :: cell
    call this%triangulation%create_cell_iterator(cell)
    do while ( .not. cell%has_finished() )
      call cell%set_for_coarsening()
      call cell%next()
    end do
    call this%triangulation%free_cell_iterator(cell)
  end subroutine set_cells_for_coarsening
  
  subroutine fill_cells_set(this)
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(inout) :: this
    integer(ip), allocatable :: cell_set_ids(:)
    class(cell_iterator_t), allocatable :: cell
    
    call memalloc(this%triangulation%get_num_cells(),cell_set_ids)
    call this%triangulation%create_cell_iterator(cell)
    do while( .not. cell%has_finished() )
      if (cell%is_local()) then
         cell_set_ids(cell%get_gid()) = cell%get_gid()
      end if
      call cell%next()
    end do
    call this%triangulation%free_cell_iterator(cell)
    call this%triangulation%fill_cells_set(cell_set_ids)
    call memfree(cell_set_ids)
    
  end subroutine fill_cells_set
  
  subroutine setup_reference_fes(this)
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(inout) :: this
    integer(ip)                                   :: istat
    logical                                       :: conformity
    class(cell_iterator_t)          , allocatable :: cell
    class(reference_fe_t), pointer     :: reference_fe_geo
    character(:)                    , allocatable :: field_type
    
    type(interpolation_t), pointer :: h_refinement_interpolation
    integer(ip), pointer :: h_refinement_subfacet_permutation(:,:,:)
    integer(ip), pointer :: h_refinement_subedge_permutation(:,:,:)
    
    if (this%test_params%get_use_void_fes()) then
      allocate(this%reference_fes(2), stat=istat)
    else
      allocate(this%reference_fes(1), stat=istat)
    end if
    check(istat==0)
    
    conformity = .true.
    if ( this%test_params%get_fe_formulation() == 'dG' ) then
      conformity = .false.
    end if
    
    field_type = field_type_scalar
    if ( this%test_params%get_laplacian_type() == 'vector' ) then
      field_type = field_type_vector
    end if
    
    call this%triangulation%create_cell_iterator(cell)
    reference_fe_geo => cell%get_reference_fe()
    this%reference_fes(TEST_POISSON_FULL) =  make_reference_fe ( topology = reference_fe_geo%get_topology(),                  &
                                                                 fe_type = fe_type_lagrangian,                                &
                                                                 num_dims = this%triangulation%get_num_dims(), &
                                                                 order = this%test_params%get_reference_fe_order(),           &
                                                                 field_type = field_type,                                     &
                                                                 conformity = conformity )
    
    if (this%test_params%get_use_void_fes()) then
      this%reference_fes(TEST_POISSON_VOID) =  make_reference_fe ( topology = reference_fe_geo%get_topology(),                  &
                                                                   fe_type = fe_type_void,                                      &
                                                                   num_dims = this%triangulation%get_num_dims(), &
                                                                   order = -1,                                                  &
                                                                   field_type = field_type,                                     &
                                                                   conformity = conformity )
    end if
    call this%triangulation%free_cell_iterator(cell)
    
    ! Take h_refinement arrays to local scope to debug with DDT
    select type( reference_fe => this%reference_fes(1)%p )
    type is (hex_lagrangian_reference_fe_t)
       h_refinement_interpolation       => reference_fe%get_h_refinement_interpolation()
       h_refinement_subfacet_permutation => reference_fe%get_h_refinement_subfacet_permutation()
       h_refinement_subedge_permutation => reference_fe%get_h_refinement_subedge_permutation()
    class default
      assert(.false.)
    end select
    
  end subroutine setup_reference_fes

  subroutine setup_fe_space(this)
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(inout) :: this
    integer(ip) :: set_ids_to_reference_fes(1,2)

    if ( this%test_params%get_laplacian_type() == 'scalar' ) then
      call this%poisson_analytical_functions%set_num_dims(this%triangulation%get_num_dims())
      call this%poisson_conditions%set_boundary_function(this%poisson_analytical_functions%get_boundary_function())
      if ( this%test_params%get_fe_formulation() == 'cG' ) then
        if ( this%test_params%get_use_void_fes() ) then
          set_ids_to_reference_fes(1,TEST_POISSON_FULL) = TEST_POISSON_FULL
          set_ids_to_reference_fes(1,TEST_POISSON_VOID) = TEST_POISSON_VOID
          call this%fe_space%create( triangulation            = this%triangulation,       &
                                     reference_fes            = this%reference_fes,       &
                                     set_ids_to_reference_fes = set_ids_to_reference_fes, &
                                     conditions               = this%poisson_conditions )
        else 
          call this%fe_space%create( triangulation       = this%triangulation, &
                                     reference_fes       = this%reference_fes, &
                                     conditions          = this%poisson_conditions )
        end if
      end if
    else
      call this%vector_poisson_analytical_functions%set_num_dims(this%triangulation%get_num_dims())
      call this%vector_poisson_conditions%set_boundary_function(this%vector_poisson_analytical_functions%get_boundary_function()) 
      if ( this%test_params%get_fe_formulation() == 'cG' ) then
        if ( this%test_params%get_use_void_fes() ) then
          set_ids_to_reference_fes(1,TEST_POISSON_FULL) = TEST_POISSON_FULL
          set_ids_to_reference_fes(1,TEST_POISSON_VOID) = TEST_POISSON_VOID
          call this%fe_space%create( triangulation       = this%triangulation,            &
                                     reference_fes       = this%reference_fes,            &
                                     set_ids_to_reference_fes = set_ids_to_reference_fes, &
                                     conditions          = this%vector_poisson_conditions )
        else
          call this%fe_space%create( triangulation       = this%triangulation, &
                                     reference_fes       = this%reference_fes, &
                                     conditions          = this%vector_poisson_conditions )
        end if
      end if
    end if
    
    if ( this%test_params%get_fe_formulation() == 'dG' ) then
      if ( this%test_params%get_use_void_fes() ) then
        set_ids_to_reference_fes(1,TEST_POISSON_FULL) = TEST_POISSON_FULL
        set_ids_to_reference_fes(1,TEST_POISSON_VOID) = TEST_POISSON_VOID
        call this%fe_space%create( triangulation            = this%triangulation,       &
                                   reference_fes            = this%reference_fes,       &
                                   set_ids_to_reference_fes = set_ids_to_reference_fes )
      else 
        call this%fe_space%create( triangulation       = this%triangulation, &
                                   reference_fes       = this%reference_fes )
      end if
    end if
    
    call this%fe_space%set_up_cell_integration()
    if ( this%test_params%get_fe_formulation() == 'dG' ) then
      call this%fe_space%set_up_facet_integration()
    end if
  end subroutine setup_fe_space
  
  subroutine refine_and_coarsen(this)
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(inout) :: this
    integer(ip) :: i
    
    do i=1,3
       
       if ( i == 2 ) then
         call this%set_cells_for_coarsening()
       else
         call this%set_cells_for_refinement()
       end if

       call this%triangulation%refine_and_coarsen()
       
       call this%fe_space%refine_and_coarsen( this%solution ) 
       
       call this%fe_space%set_up_cell_integration()
       
       if ( this%test_params%get_laplacian_type() == 'scalar' ) then
         call this%check_solution()
       else
         call this%check_solution_vector()
       end if
       
    end do  
    
  end subroutine refine_and_coarsen
  
  subroutine setup_system (this)
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(inout) :: this

    if ( this%test_params%get_laplacian_type() == 'scalar' ) then    
      if ( this%test_params%get_fe_formulation() == 'cG' ) then
        call this%poisson_cG_integration%set_analytical_functions(this%poisson_analytical_functions)
        call this%fe_affine_operator%create ( sparse_matrix_storage_format      = csr_format,                               &
                                              diagonal_blocks_symmetric_storage = [ .true. ],                               &
                                              diagonal_blocks_symmetric         = [ .true. ],                               &
                                              diagonal_blocks_sign              = [ SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE ], &
                                              fe_space                          = this%fe_space,                            &
                                              discrete_integration              = this%poisson_cG_integration )
      else
        call this%poisson_dG_integration%set_analytical_functions(this%poisson_analytical_functions)
        call this%poisson_dG_integration%set_poisson_conditions(this%poisson_conditions)
        call this%fe_affine_operator%create ( sparse_matrix_storage_format      = csr_format, &
                                              diagonal_blocks_symmetric_storage = [ .true. ], &
                                              diagonal_blocks_symmetric         = [ .true. ], &
                                              diagonal_blocks_sign              = [ SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE ], &
                                              fe_space                          = this%fe_space, &
                                              discrete_integration              = this%poisson_dG_integration )
      end if
    else
      if ( this%test_params%get_fe_formulation() == 'cG' ) then
        call this%vector_poisson_integration%set_source_term(this%vector_poisson_analytical_functions%get_source_term())
        call this%fe_affine_operator%create ( sparse_matrix_storage_format      = csr_format,                               &
                                              diagonal_blocks_symmetric_storage = [ .true. ],                               &
                                              diagonal_blocks_symmetric         = [ .true. ],                               &
                                              diagonal_blocks_sign              = [ SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE ], &
                                              fe_space                          = this%fe_space,                            &
                                              discrete_integration              = this%vector_poisson_integration )
      else
        mcheck(.false.,'Vector poisson dG integration is not implemented')
      end if
    end if
    call this%solution%create(this%fe_space) 
    if ( this%test_params%get_fe_formulation() == 'cG' ) then
      call this%fe_space%interpolate_dirichlet_values(this%solution)
      if ( this%test_params%get_laplacian_type() == 'scalar' ) then
        call this%poisson_cG_integration%set_fe_function(this%solution)
      else
        call this%vector_poisson_integration%set_fe_function(this%solution)
      end if
    end if
  end subroutine setup_system
  
  subroutine setup_solver (this)
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(inout) :: this
    integer :: FPLError
    type(parameterlist_t) :: parameter_list
    integer :: iparm(64)
    class(matrix_t), pointer       :: matrix

    call parameter_list%init()
#ifdef ENABLE_MKL
    FPLError = parameter_list%set(key = dls_type_key,        value = pardiso_mkl)
    FPLError = FPLError + parameter_list%set(key = pardiso_mkl_matrix_type,   value = pardiso_mkl_spd)
    FPLError = FPLError + parameter_list%set(key = pardiso_mkl_message_level, value = 0)
    iparm = 0
    FPLError = FPLError + parameter_list%set(key = pardiso_mkl_iparm,         value = iparm)
    assert(FPLError == 0)
    
    call this%direct_solver%set_type_from_pl(parameter_list)
    call this%direct_solver%set_parameters_from_pl(parameter_list)
    
    matrix => this%fe_affine_operator%get_matrix()
    select type(matrix)
    class is (sparse_matrix_t)  
       call this%direct_solver%set_matrix(matrix)
    class DEFAULT
       assert(.false.) 
    end select
#else    
    FPLError = parameter_list%set(key = ils_rtol_key, value = 1.0e-12_rp)
    !FPLError = FPLError + parameter_list%set(key = ils_output_frequency, value = 30)
    FPLError = parameter_list%set(key = ils_max_num_iterations_key, value = 5000)
    assert(FPLError == 0)
    call this%iterative_linear_solver%create(this%fe_space%get_environment())
    call this%iterative_linear_solver%set_type_from_string(cg_name)
    call this%iterative_linear_solver%set_parameters_from_pl(parameter_list)
    call this%iterative_linear_solver%set_operators(this%fe_affine_operator%get_tangent(), .identity. this%fe_affine_operator) 
#endif
    call parameter_list%free()
  end subroutine setup_solver
  
  
  subroutine assemble_system (this)
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(inout) :: this
    class(matrix_t)                  , pointer       :: matrix
    class(vector_t)                  , pointer       :: rhs
    call this%fe_affine_operator%compute()
    rhs                => this%fe_affine_operator%get_translation()
    matrix             => this%fe_affine_operator%get_matrix()
    
    select type(matrix)
    class is (sparse_matrix_t)  
       !call matrix%print_matrix_market(6) 
    class DEFAULT
       assert(.false.) 
    end select
    
    select type(rhs)
    class is (serial_scalar_array_t)  
    !call rhs%print(6) 
    class DEFAULT
       assert(.false.) 
    end select
  end subroutine assemble_system
  
  
  subroutine solve_system(this)
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(inout) :: this
    class(matrix_t)                         , pointer       :: matrix
    class(vector_t)                         , pointer       :: rhs
    class(vector_t)                         , pointer       :: dof_values

    matrix     => this%fe_affine_operator%get_matrix()
    rhs        => this%fe_affine_operator%get_translation()
    dof_values => this%solution%get_free_dof_values()
    
#ifdef ENABLE_MKL    
    call this%direct_solver%solve(this%fe_affine_operator%get_translation(), dof_values)
#else
    call this%iterative_linear_solver%solve(this%fe_affine_operator%get_translation(), &
                                            dof_values)
#endif    
    
    call this%fe_space%update_hanging_dof_values(this%solution)
    !call this%solution%update_fixed_dof_values(this%fe_space)
    
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
    class(test_h_adaptive_poisson_driver_t), intent(inout) :: this
    type(error_norms_scalar_t) :: error_norm
    real(rp) :: mean, l1, l2, lp, linfty, h1, h1_s, w1p_s, w1p, w1infty_s, w1infty
    real(rp) :: error_tolerance
    
    call error_norm%create(this%fe_space,1)
    mean = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, mean_norm)   
    l1 = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, l1_norm)   
    l2 = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, l2_norm)   
    lp = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, lp_norm)   
    linfty = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, linfty_norm)   
    h1_s = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, h1_seminorm) 
    h1 = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, h1_norm) 
    w1p_s = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, w1p_seminorm)   
    w1p = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, w1p_norm)   
    w1infty_s = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, w1infty_seminorm) 
    w1infty = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, w1infty_norm)

#ifdef ENABLE_MKL    
    error_tolerance = 1.0e-08
#else
    error_tolerance = 1.0e-06
#endif    
    
    write(*,'(a20,e32.25)') 'mean_norm:', mean; check ( abs(mean) < error_tolerance )
    write(*,'(a20,e32.25)') 'l1_norm:', l1; check ( l1 < error_tolerance )
    write(*,'(a20,e32.25)') 'l2_norm:', l2; check ( l2 < error_tolerance )
    write(*,'(a20,e32.25)') 'lp_norm:', lp; check ( lp < error_tolerance )
    write(*,'(a20,e32.25)') 'linfnty_norm:', linfty; check ( linfty < error_tolerance )
    write(*,'(a20,e32.25)') 'h1_seminorm:', h1_s; check ( h1_s < error_tolerance )
    write(*,'(a20,e32.25)') 'h1_norm:', h1; check ( h1 < error_tolerance )
    write(*,'(a20,e32.25)') 'w1p_seminorm:', w1p_s; check ( w1p_s < error_tolerance )
    write(*,'(a20,e32.25)') 'w1p_norm:', w1p; check ( w1p < error_tolerance )
    write(*,'(a20,e32.25)') 'w1infty_seminorm:', w1infty_s; check ( w1infty_s < error_tolerance )
    write(*,'(a20,e32.25)') 'w1infty_norm:', w1infty; check ( w1infty < error_tolerance )
    call error_norm%free()
  end subroutine check_solution
  
  subroutine check_solution_vector(this)
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(in) :: this
    type(error_norms_vector_t) :: error_norm
    real(rp) :: mean, l1, l2, lp, linfty, h1, h1_s, w1p_s, w1p, w1infty_s, w1infty
    real(rp) :: error_tolerance
    
    call error_norm%create(this%fe_space,1)
    mean = error_norm%compute(this%vector_poisson_analytical_functions%get_solution_function(), this%solution, mean_norm)   
    l1 = error_norm%compute(this%vector_poisson_analytical_functions%get_solution_function(), this%solution, l1_norm)   
    l2 = error_norm%compute(this%vector_poisson_analytical_functions%get_solution_function(), this%solution, l2_norm)   
    lp = error_norm%compute(this%vector_poisson_analytical_functions%get_solution_function(), this%solution, lp_norm)   
    linfty = error_norm%compute(this%vector_poisson_analytical_functions%get_solution_function(), this%solution, linfty_norm)   
    h1_s = error_norm%compute(this%vector_poisson_analytical_functions%get_solution_function(), this%solution, h1_seminorm) 
    h1 = error_norm%compute(this%vector_poisson_analytical_functions%get_solution_function(), this%solution, h1_norm) 
    w1p_s = error_norm%compute(this%vector_poisson_analytical_functions%get_solution_function(), this%solution, w1p_seminorm)   
    w1p = error_norm%compute(this%vector_poisson_analytical_functions%get_solution_function(), this%solution, w1p_norm)   
    w1infty_s = error_norm%compute(this%vector_poisson_analytical_functions%get_solution_function(), this%solution, w1infty_seminorm) 
    w1infty = error_norm%compute(this%vector_poisson_analytical_functions%get_solution_function(), this%solution, w1infty_norm)

#ifdef ENABLE_MKL    
    error_tolerance = 1.0e-08
#else
    error_tolerance = 1.0e-06
#endif    
    
    write(*,'(a20,e32.25)') 'mean_norm:', mean; check ( abs(mean) < error_tolerance )
    write(*,'(a20,e32.25)') 'l1_norm:', l1; check ( l1 < error_tolerance )
    write(*,'(a20,e32.25)') 'l2_norm:', l2; check ( l2 < error_tolerance )
    write(*,'(a20,e32.25)') 'lp_norm:', lp; check ( lp < error_tolerance )
    write(*,'(a20,e32.25)') 'linfnty_norm:', linfty; check ( linfty < error_tolerance )
    write(*,'(a20,e32.25)') 'h1_seminorm:', h1_s; check ( h1_s < error_tolerance )
    write(*,'(a20,e32.25)') 'h1_norm:', h1; check ( h1 < error_tolerance )
    write(*,'(a20,e32.25)') 'w1p_seminorm:', w1p_s; check ( w1p_s < error_tolerance )
    write(*,'(a20,e32.25)') 'w1p_norm:', w1p; check ( w1p < error_tolerance )
    write(*,'(a20,e32.25)') 'w1infty_seminorm:', w1infty_s; check ( w1infty_s < error_tolerance )
    write(*,'(a20,e32.25)') 'w1infty_norm:', w1infty; check ( w1infty < error_tolerance )
    call error_norm%free()
  end subroutine check_solution_vector
  
  subroutine write_solution(this)
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(in) :: this
    type(output_handler_t)                   :: oh
    character(len=:), allocatable            :: path
    character(len=:), allocatable            :: prefix
    real(rp),allocatable :: cell_vector(:)
    integer(ip) :: N, P, pid, i
    class(cell_iterator_t), allocatable :: cell
    if(this%test_params%get_write_solution()) then
        path = this%test_params%get_dir_path_out()
        prefix = this%test_params%get_prefix()
        call oh%create()
        call oh%attach_fe_space(this%fe_space)
        call oh%add_fe_function(this%solution, 1, 'solution')
        call oh%add_fe_function(this%solution, 1, 'grad_solution', grad_diff_operator)
        call memalloc(this%triangulation%get_num_cells(),cell_vector,__FILE__,__LINE__)
        
        call this%triangulation%create_cell_iterator(cell)
        if (this%test_params%get_use_void_fes()) then
          do while( .not. cell%has_finished() )
            cell_vector(cell%get_gid()) = cell%get_set_id()
            call cell%next()
          end do
          call oh%add_cell_vector(cell_vector,'cell_set_ids')
        else
          N=this%triangulation%get_num_cells()
          P=6
          do pid=0, P-1
              i=0
              do while ( i < (N*(pid+1))/P - (N*pid)/P ) 
                cell_vector(cell%get_gid()) = pid 
                call cell%next()
                i=i+1
              end do
          end do
          call oh%add_cell_vector(cell_vector,'cell_set_ids')
        end if
        call this%triangulation%free_cell_iterator(cell)
        
        call oh%open(path, prefix)
        call oh%write()
        call oh%close()
        call oh%free()
        call memfree(cell_vector,__FILE__,__LINE__)
    endif
  end subroutine write_solution

  subroutine write_filling_curve(this)
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(in) :: this

    integer(ip) :: Nn, Ne
    real(rp), allocatable :: x(:), y(:), z(:)
    integer(ip), allocatable :: cell_type(:), offset(:), connect(:)
    class(cell_iterator_t)      , allocatable :: cell
    type(point_t), allocatable :: coords(:)
    integer(ip) :: istat, k
    real(rp) ::  xc,yc
    integer(ip), parameter :: max_num_cell_nodes = 4
    integer(ip), parameter :: vtk_1d_elem_id = 3
    character(len=*), parameter :: filename_out = 'output/filling_curve.vtu'
    integer(ip) :: E_IO
    if(this%test_params%get_write_solution()) then

      Nn = this%triangulation%get_num_cells()
      Ne = Nn - 1

      call memalloc ( Nn, x, __FILE__, __LINE__ )
      call memalloc ( Nn, y, __FILE__, __LINE__ )
      call memalloc ( Nn, z, __FILE__, __LINE__ )
      call memalloc ( Ne, cell_type, __FILE__, __LINE__ )
      call memalloc ( Ne, offset   , __FILE__, __LINE__ )
      call memalloc ( 2*Ne, connect  , __FILE__, __LINE__ )

      call this%triangulation%create_cell_iterator(cell)
      allocate(coords(max_num_cell_nodes),stat=istat); check(istat==0)

      do while ( .not. cell%has_finished() )

        call cell%get_nodes_coordinates(coords)
        xc = 0.0
        yc = 0.0
        do k=1,max_num_cell_nodes
          xc = xc + (1.0/max_num_cell_nodes)*coords(k)%get(1)
          yc = yc + (1.0/max_num_cell_nodes)*coords(k)%get(2)
        end do

        x(cell%get_gid()) = xc;
        y(cell%get_gid()) = yc;
        z(cell%get_gid()) = 0.0;

        if (cell%get_gid()>1) then
          connect(  2*(cell%get_gid()-1)-1  ) = cell%get_gid()-2
          connect(  2*(cell%get_gid()-1)    ) = cell%get_gid()-1
          offset( cell%get_gid()-1 ) = 2*(cell%get_gid()-1)
          cell_type( cell%get_gid()-1 ) = vtk_1d_elem_id
        end if

        call cell%next()
      end do

      deallocate(coords,stat=istat); check(istat==0)
      call this%triangulation%free_cell_iterator(cell)

      E_IO = VTK_INI_XML(output_format = 'ascii', filename = filename_out, mesh_topology = 'UnstructuredGrid')
      E_IO = VTK_GEO_XML(NN = Nn, NC = Ne, X = x, Y = y, Z = z)
      E_IO = VTK_CON_XML(NC = Ne, connect = connect, offset = offset, cell_type = int(cell_type,I1P) )
      E_IO = VTK_GEO_XML()
      E_IO = VTK_END_XML()

      call memfree ( x, __FILE__, __LINE__ )
      call memfree ( y, __FILE__, __LINE__ )
      call memfree ( z, __FILE__, __LINE__ )
      call memfree ( cell_type, __FILE__, __LINE__ )
      call memfree ( offset   , __FILE__, __LINE__ )
      call memfree ( connect  , __FILE__, __LINE__ )

    endif
  end subroutine write_filling_curve
  
  subroutine run_simulation(this) 
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(inout) :: this    
    call this%free()
    call this%setup_triangulation()
    call this%setup_reference_fes()
    call this%setup_fe_space()
    call this%setup_system()
    call this%assemble_system()
    call this%setup_solver()
    call this%solve_system()
    if ( this%test_params%get_laplacian_type() == 'scalar' ) then
      call this%check_solution()
    else
      call this%check_solution_vector()
    end if
    call this%refine_and_coarsen()
    call this%write_solution()
    !call this%write_filling_curve()
    call this%free()
  end subroutine run_simulation
  
  subroutine free(this)
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(inout) :: this
    integer(ip) :: i, istat
    
    call this%solution%free()
    
#ifdef ENABLE_MKL        
    call this%direct_solver%free()
#else
    call this%iterative_linear_solver%free()
#endif
    
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
  end subroutine free  
  

  
end module test_h_adaptive_poisson_driver_names
