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

module stokes_driver_names
  use fempar_names
  use unfitted_triangulations_names
  use unfitted_fe_spaces_names
  use level_set_functions_gallery_names
  use unfitted_vtk_writer_names
  use unfitted_solution_checker_names
  use unfitted_solution_checker_vector_names
  use level_set_functions_gallery_names
  use unfitted_vtk_writer_names
  use stokes_params_names
  use stokes_cG_discrete_integration_names
  use stokes_conditions_names
  use stokes_analytical_functions_names
  use IR_Precision ! VTK_IO
  use Lib_VTK_IO ! VTK_IO
  use plot_aggregates_utils_names
    
# include "debug.i90"

  implicit none
  private

  integer(ip), parameter :: SET_ID_FULL = 1
  integer(ip), parameter :: SET_ID_VOID = 2

  integer(ip), parameter :: POS_FULL_U  = 1
  integer(ip), parameter :: POS_VOID_U  = 2
  integer(ip), parameter :: POS_FULL_P  = 3
  integer(ip), parameter :: POS_VOID_P  = 4


  type stokes_driver_t 
     private 
     
     type(stokes_params_t)                        :: test_params
     type(ParameterList_t)                        :: parameter_list
     type(unfitted_p4est_serial_triangulation_t)  :: triangulation
     class(level_set_function_t), allocatable     :: level_set_function
     type(serial_unfitted_hp_adaptive_fe_space_t) :: fe_space 
     type(p_reference_fe_t), allocatable          :: reference_fes(:) 
     type(stokes_cG_discrete_integration_t)       :: cG_integration
     type(stokes_conditions_t)                    :: conditions
     type(stokes_analytical_functions_t)          :: analytical_functions
     type(fe_affine_operator_t)                   :: fe_affine_operator
#ifdef ENABLE_MKL     
     type(direct_solver_t)                        :: direct_solver
     type(iterative_linear_solver_t)              :: iterative_linear_solver
#else     
     type(iterative_linear_solver_t)              :: iterative_linear_solver
#endif     
     type(fe_function_t)                          :: solution

   contains
     procedure                  :: run_simulation
     procedure        , private :: parse_command_line_parameters
     procedure        , private :: setup_levelset
     procedure        , private :: setup_triangulation
     procedure        , private :: fill_cells_set
     procedure        , private :: setup_reference_fes
     procedure        , private :: setup_fe_space
     procedure        , private :: setup_system
     procedure        , private :: setup_solver
     procedure        , private :: assemble_system
     procedure        , private :: solve_system     
     procedure        , private :: fix_pressure
     procedure        , private :: check_solution
     procedure        , private :: write_solution
     procedure        , private :: free
  end type stokes_driver_t

  ! Types
  public :: stokes_driver_t

contains

  subroutine parse_command_line_parameters(this)
    implicit none
    class(stokes_driver_t ), intent(inout) :: this
    call this%test_params%create()
    call this%test_params%parse(this%parameter_list)
  end subroutine parse_command_line_parameters

  subroutine setup_levelset(this)
    implicit none
    class(stokes_driver_t ), target, intent(inout) :: this

    integer(ip) :: num_dime
    integer(ip) :: istat
    class(level_set_function_t), pointer :: levset
    type(level_set_function_factory_t) :: level_set_factory
    real(rp) :: dom1d(2)
    real(rp) :: dom3d(6)

    ! Get number of dimensions form input
    massert( this%parameter_list%isPresent   (key = num_dims_key), 'Use -tt structured' )
    assert( this%parameter_list%isAssignable (key = num_dims_key, value=num_dime) )
    istat = this%parameter_list%get          (key = num_dims_key, value=num_dime); check(istat==0)

    ! Create the desired type of level set function
    call level_set_factory%create(this%test_params%get_levelset_function_type(), this%level_set_function)

    ! Set options of the base class
    call this%level_set_function%set_num_dims(num_dime)
    call this%level_set_function%set_tolerance(this%test_params%get_levelset_tolerance())
    dom1d = this%test_params%get_domain_limits()
    mcheck(dom1d(2)>dom1d(1),'Upper limit has to be bigger than lower limit')
    dom3d(1) = dom1d(1)
    dom3d(2) = dom1d(2)
    dom3d(3) = dom1d(1)
    dom3d(4) = dom1d(2)
    dom3d(5) = dom1d(1)
    dom3d(6) = dom1d(2)
    call this%level_set_function%set_domain(dom3d)

    ! Set options of the derived classes
    ! TODO a parameter list would be better to define the level set function together with its parameters
    levset => this%level_set_function
    select type ( levset )
      class is (level_set_sphere_t)
        call levset%set_radius(0.9)
        call levset%set_center([0.0,0.0,0.0])
    end select

  end subroutine setup_levelset

  subroutine setup_triangulation(this)
    implicit none
    class(stokes_driver_t), intent(inout) :: this
    class(vef_iterator_t),allocatable :: vef

    class(cell_iterator_t), allocatable :: cell
    integer(ip) :: ilev
    integer(ip) :: max_levels
    integer(ip) :: diri_set_id_u
    integer(ip) :: diri_set_id_u_and_p
    logical :: first_interior_vertex
    integer(ip) :: ivef
    real(rp) :: target_size

    ! Create the triangulation, with the levelset function
    call this%triangulation%create(this%parameter_list,this%level_set_function)


    ! Create initial refined mesh
    select case ( trim(this%test_params%get_refinement_pattern()) )
      case ('uniform')

        max_levels = this%test_params%get_max_level()
        do ilev = 1, max_levels
          call this%triangulation%create_cell_iterator(cell)
          do while (.not. cell%has_finished())
            call cell%set_for_refinement()
            call cell%next()
          end do
          call this%triangulation%refine_and_coarsen()
          call this%triangulation%clear_refinement_and_coarsening_flags()
          call this%triangulation%free_cell_iterator(cell)
        end do
        call this%triangulation%update_cut_cells(this%level_set_function)

      case ('adaptive-1')

        max_levels = this%test_params%get_max_level()
        do ilev = 1, max_levels
          call this%triangulation%create_cell_iterator(cell)
          do while (.not. cell%has_finished())
            if (ilev <= 2) then
              call cell%set_for_refinement()
            else if (ilev == max_levels) then
              if (cell%is_interior()) then
                call cell%set_for_do_nothing()
              else if (cell%is_cut()) then
                call cell%set_for_refinement()
              else
                call cell%set_for_coarsening()
              end if
            else
              if (cell%is_interior()) then
                call cell%set_for_refinement()
              else if (cell%is_cut()) then
                call cell%set_for_refinement()
              else
                call cell%set_for_coarsening()
              end if
            end if
            call cell%next()
          end do
          call this%triangulation%refine_and_coarsen()
          call this%triangulation%clear_refinement_and_coarsening_flags()
          call this%triangulation%update_cut_cells(this%level_set_function)
          call this%triangulation%free_cell_iterator(cell)
        end do

      case ('adaptive-2')

        max_levels = this%test_params%get_max_level()
        do ilev = 1, max_levels
          call this%triangulation%create_cell_iterator(cell)
          do while (.not. cell%has_finished())
            if (ilev <= 2) then
              call cell%set_for_refinement()
            else
              if (cell%is_interior()) then
                call cell%set_for_refinement()
              else if (cell%is_cut()) then
                call cell%set_for_refinement()
              else
                call cell%set_for_coarsening()
              end if
            end if
            call cell%next()
          end do
          call this%triangulation%refine_and_coarsen()
          call this%triangulation%clear_refinement_and_coarsening_flags()
          call this%triangulation%update_cut_cells(this%level_set_function)
          call this%triangulation%free_cell_iterator(cell)
        end do

      case ('adaptive-3')

        max_levels = this%test_params%get_max_level()
        do ilev = 1, max_levels
          call this%triangulation%create_cell_iterator(cell)
          do while (.not. cell%has_finished())
            if (ilev <= 2) then
              call cell%set_for_refinement()
            else
              if (cell%is_interior()) then
                call cell%set_for_refinement()
              else if (cell%is_cut()) then
                call cell%set_for_refinement()
              else
                call cell%set_for_coarsening()
              end if
            end if
            call cell%next()
          end do
          call this%triangulation%refine_and_coarsen()
          call this%triangulation%clear_refinement_and_coarsening_flags()
          call this%triangulation%update_cut_cells(this%level_set_function)
          call this%triangulation%free_cell_iterator(cell)
        end do

        target_size = 1.0/(2.0**this%test_params%get_max_level())
        call this%fe_space%refine_mesh_for_small_aggregates(this%triangulation,target_size,this%level_set_function)

      case ('debug-1')

        call this%triangulation%create_cell_iterator(cell)
        call cell%set_for_refinement()
        call this%triangulation%refine_and_coarsen()
        call this%triangulation%clear_refinement_and_coarsening_flags()
        call this%triangulation%update_cut_cells(this%level_set_function)
        call cell%set_gid(2)
        call cell%set_for_refinement()
        call this%triangulation%refine_and_coarsen()
        call this%triangulation%clear_refinement_and_coarsening_flags()
        call this%triangulation%update_cut_cells(this%level_set_function)
        call this%triangulation%free_cell_iterator(cell)

      case default
            mcheck(.false.,'Refinement pattern `'//trim(this%test_params%get_refinement_pattern())//'` not known')
    end select

    ! Impose Dirichlet
    if (this%test_params%is_strong_dirichlet_on_fitted_boundary()) then
      diri_set_id_u = 1
      diri_set_id_u_and_p = 2
    else
      diri_set_id_u = 0
      diri_set_id_u_and_p = 0
    end if
    if ( trim(this%test_params%get_triangulation_type()) == 'structured' ) then

       call this%triangulation%create_vef_iterator(vef)
       call this%triangulation%create_cell_iterator(cell)

       ! For velocities
       do while ( .not. vef%has_finished() )
         if(vef%is_at_boundary()) then
            call vef%set_set_id(diri_set_id_u)
         else
            call vef%set_set_id(0)
         end if
         call vef%next()
       end do

       ! For pressures
       call cell%first()
       first_interior_vertex = .true.
       do while ( .not. cell%has_finished() )
         if (cell%is_interior()) then
           do ivef = 1, cell%get_num_vefs()
             call cell%get_vef(ivef,vef)
             if (vef%is_proper() .and. vef%get_dim() == 0) then
               if (first_interior_vertex) then
                 call vef%set_set_id(diri_set_id_u_and_p)
                 first_interior_vertex = .false.
                 exit
               end if
             end if
           end do
           if (.not. first_interior_vertex) then
             exit
           end if
         end if
         call cell%next()
       end do

       massert(.not. first_interior_vertex,'No interior vertex in interior cell found: refine your mesh!')

       call this%triangulation%free_vef_iterator(vef)
       call this%triangulation%free_cell_iterator(cell)

    end if

    
    !call this%triangulation%print()

  end subroutine setup_triangulation
  
  subroutine fill_cells_set(this)
    implicit none
    class(stokes_driver_t), intent(inout) :: this

    integer(ip), allocatable :: cell_set_ids(:)
    class(cell_iterator_t), allocatable :: cell
    integer(ip) :: set_id
    
    call memalloc(this%triangulation%get_num_cells(),cell_set_ids)
    call this%triangulation%create_cell_iterator(cell)
    do while( .not. cell%has_finished() )
      if (cell%is_exterior()) then
        set_id = SET_ID_VOID
      else
        set_id = SET_ID_FULL
      end if
      cell_set_ids(cell%get_gid()) = set_id
      call cell%next()
    end do
    call this%triangulation%fill_cells_set(cell_set_ids)
    call this%triangulation%free_cell_iterator(cell)
    call memfree(cell_set_ids)
    
  end subroutine fill_cells_set
  
  subroutine setup_reference_fes(this)
    implicit none
    class(stokes_driver_t), intent(inout) :: this
    integer(ip) :: istat    
    class(cell_iterator_t), allocatable :: cell
    class(reference_fe_t),  pointer     :: reference_fe
    character(:),           allocatable :: field_type
    integer(ip) :: order_u
    integer(ip) :: order_p

    order_p = this%test_params%get_reference_fe_order()
    order_u = order_p + 1
    
    allocate(this%reference_fes(4), stat=istat)
    check(istat==0)
    
    call this%triangulation%create_cell_iterator(cell)
    reference_fe => cell%get_reference_fe()
    this%reference_fes(POS_FULL_U) =  make_reference_fe ( topology = reference_fe%get_topology(), &
                                                 fe_type = fe_type_lagrangian, &
                                                 num_dims = this%triangulation%get_num_dims(), &
                                                 order = order_u, &
                                                 field_type = field_type_vector, &
                                                 conformity = .true. )
    this%reference_fes(POS_VOID_U) =  make_reference_fe ( topology = reference_fe%get_topology(), &
                                                 fe_type = fe_type_void, &
                                                 num_dims = this%triangulation%get_num_dims(), &
                                                 order = order_u, &
                                                 field_type = field_type_vector, &
                                                 conformity = .true. )
    this%reference_fes(POS_FULL_P) =  make_reference_fe ( topology = reference_fe%get_topology(), &
                                                 fe_type = fe_type_lagrangian, &
                                                 num_dims = this%triangulation%get_num_dims(), &
                                                 order = order_p, &
                                                 field_type = field_type_scalar, &
                                                 conformity = .true. )
    this%reference_fes(POS_VOID_P) =  make_reference_fe ( topology = reference_fe%get_topology(), &
                                                 fe_type = fe_type_void, &
                                                 num_dims = this%triangulation%get_num_dims(), &
                                                 order = order_p, &
                                                 field_type = field_type_scalar, &
                                                 conformity = .true. )
    call this%triangulation%free_cell_iterator(cell)
    
  end subroutine setup_reference_fes

  subroutine setup_fe_space(this)
    implicit none
    class(stokes_driver_t), intent(inout) :: this

    integer(ip) :: set_ids_to_reference_fes(2,2) ! num_fields x void/non_void
    class(vector_function_t) , pointer :: fun_u
    class(scalar_function_t) , pointer :: fun_p

    set_ids_to_reference_fes(U_FIELD_ID,SET_ID_FULL) = POS_FULL_U
    set_ids_to_reference_fes(U_FIELD_ID,SET_ID_VOID) = POS_VOID_U
    set_ids_to_reference_fes(P_FIELD_ID,SET_ID_FULL) = POS_FULL_P
    set_ids_to_reference_fes(P_FIELD_ID,SET_ID_VOID) = POS_VOID_P
    
    call this%analytical_functions%set_num_dims(this%triangulation%get_num_dims())
    call this%analytical_functions%set_case_id(this%test_params%get_case_id())
    call this%analytical_functions%set_degree(this%test_params%get_reference_fe_order())

    fun_u => this%analytical_functions%get_solution_function_u()
    fun_p => this%analytical_functions%get_solution_function_p()
    call this%conditions%set_num_dims(this%triangulation%get_num_dims())
    call this%conditions%set_boundary_function(fun_u,fun_p)

    call this%fe_space%set_use_constraints(this%test_params%get_use_constraints())
    call this%fe_space%create( triangulation            = this%triangulation,      &
                               conditions               = this%conditions, &
                               reference_fes            = this%reference_fes,&
                               set_ids_to_reference_fes = set_ids_to_reference_fes)
    call this%fe_space%set_up_cell_integration()    
    
  end subroutine setup_fe_space
  
  subroutine setup_system (this)
    implicit none
    class(stokes_driver_t), intent(inout) :: this

    integer(ip) :: iounit

    call this%cG_integration%set_analytical_functions(this%analytical_functions)
    call this%cG_integration%set_unfitted_boundary_is_dirichlet(this%test_params%get_unfitted_boundary_is_dirichlet())
    call this%cG_integration%set_is_constant_nitches_beta(this%test_params%get_is_constant_nitches_beta())

    call this%fe_affine_operator%create (   sparse_matrix_storage_format      = csr_format,                                  &
                                            diagonal_blocks_symmetric_storage = [ .true. ],                               &
                                            diagonal_blocks_symmetric         = [ .true. ],                               &
                                            diagonal_blocks_sign              = [ SPARSE_MATRIX_SIGN_INDEFINITE ],         &
                                            fe_space                          = this%fe_space,                             &
                                            discrete_integration              = this%cG_integration )
    call this%solution%create(this%fe_space)
    call this%fe_space%interpolate_dirichlet_values(this%solution)
    call this%cG_integration%set_fe_function(this%solution)

    ! Write some info
    if (this%test_params%get_write_aggr_info()) then
      iounit = io_open(file=this%test_params%get_dir_path_out()//this%test_params%get_prefix()//'_aggr_info.csv',action='write')
      check(iounit>0)
      call this%fe_space%print_debug_info(iounit)
      call io_close(iounit)
    end if
    
  end subroutine setup_system
  
  subroutine setup_solver (this)
    implicit none
    class(stokes_driver_t), intent(inout) :: this
    integer :: FPLError
    type(parameterlist_t) :: parameter_list
    integer :: iparm(64)
    class(matrix_t), pointer       :: matrix

    call parameter_list%init()
#ifdef ENABLE_MKL

    if (this%test_params%get_lin_solver_type()=='pardiso') then
      FPLError = parameter_list%set(key = direct_solver_type,        value = pardiso_mkl)
      FPLError = FPLError + parameter_list%set(key = pardiso_mkl_matrix_type,   value = pardiso_mkl_sin)
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

    else if (this%test_params%get_lin_solver_type()=='minres') then
      FPLError = parameter_list%set(key = ils_rtol, value = 1.0e-10_rp)
      FPLError = FPLError + parameter_list%set(key = ils_output_frequency, value = 50)
      FPLError = FPLError + parameter_list%set(key = ils_max_num_iterations, value = 5000000)
      assert(FPLError == 0)
      call this%iterative_linear_solver%create(this%fe_space%get_environment())
      call this%iterative_linear_solver%set_type_from_string(minres_name)
      call this%iterative_linear_solver%set_parameters_from_pl(parameter_list)
      call this%iterative_linear_solver%set_operators(this%fe_affine_operator, .identity. this%fe_affine_operator) 
    else
      mcheck(.false.,'Unknown linear solver type: '//this%test_params%get_lin_solver_type())
    end if
    
#else    
    FPLError = parameter_list%set(key = ils_rtol, value = 1.0e-12_rp)
    FPLError = FPLError + parameter_list%set(key = ils_output_frequency, value = 5)
    FPLError = FPLError + parameter_list%set(key = ils_max_num_iterations, value = 50000)
    assert(FPLError == 0)
    call this%iterative_linear_solver%create(this%fe_space%get_environment())
    call this%iterative_linear_solver%set_type_from_string(minres_name)
    call this%iterative_linear_solver%set_parameters_from_pl(parameter_list)
    call this%iterative_linear_solver%set_operators(this%fe_affine_operator, .identity. this%fe_affine_operator) 
#endif
    call parameter_list%free()

  end subroutine setup_solver
  
  
  subroutine assemble_system (this)
    implicit none
    class(stokes_driver_t), intent(inout) :: this
    class(matrix_t)                  , pointer       :: matrix
    class(vector_t)                  , pointer       :: rhs
    !class(vector_t)                  , pointer       :: sol
    !type(serial_scalar_array_t)     :: r

    integer(ip) :: iounit


    call this%fe_affine_operator%numerical_setup()
    rhs                => this%fe_affine_operator%get_translation()
    matrix             => this%fe_affine_operator%get_matrix()

    select type(matrix)
    class is (sparse_matrix_t)  
       if (this%test_params%get_write_matrix()) then
       iounit = io_open(file=this%test_params%get_dir_path_out()//this%test_params%get_prefix()//'_matrix.mm',action='write')
       check(iounit>0)
       call matrix%print_matrix_market(iounit) 
       call io_close(iounit)
       end if
    class DEFAULT
       assert(.false.) 
    end select
    
    select type(rhs)
    class is (serial_scalar_array_t)  
       !if (this%test_params%get_write_matrix()) then
       !iounit = io_open(file=this%test_params%get_dir_path_out()//this%test_params%get_prefix()//'_vector.mm',action='write')
       !check(iounit>0)
       !call rhs%print(iounit) 
       !call io_close(iounit)
       !end if
    class DEFAULT
       assert(.false.) 
    end select

    !! Interpolate solution
    !call this%fe_space%interpolate(U_FIELD_ID,this%analytical_functions%get_solution_function_u(),this%solution)
    !call this%fe_space%interpolate(P_FIELD_ID,this%analytical_functions%get_solution_function_p(),this%solution)
    !sol => this%solution%get_free_dof_values()
    !select type(matrix)
    !class is (sparse_matrix_t)  
    !  select type(rhs)
    !  class is (serial_scalar_array_t)  
    !  select type(sol)
    !    class is (serial_scalar_array_t)  
    !      call r%clone(rhs)
    !      call r%scal(-1.0,rhs)
    !      call matrix%apply_add(sol,r)
    !      call r%print(6)
    !    class DEFAULT
    !       assert(.false.) 
    !    end select
    !  class DEFAULT
    !     assert(.false.) 
    !  end select
    !class DEFAULT
    !   assert(.false.) 
    !end select
    !call r%free()

  end subroutine assemble_system
  
  
  subroutine solve_system(this)
    implicit none
    class(stokes_driver_t), intent(inout) :: this
    class(matrix_t)                         , pointer       :: matrix
    class(vector_t)                         , pointer       :: rhs
    class(vector_t)                         , pointer       :: dof_values

    matrix     => this%fe_affine_operator%get_matrix()
    rhs        => this%fe_affine_operator%get_translation()
    dof_values => this%solution%get_free_dof_values()
    
#ifdef ENABLE_MKL    
    if (this%test_params%get_lin_solver_type()=='pardiso') then
      call this%direct_solver%solve(this%fe_affine_operator%get_translation(), dof_values)
    else if (this%test_params%get_lin_solver_type()=='minres') then
      call this%iterative_linear_solver%solve(this%fe_affine_operator%get_translation(), &
                                            dof_values)
    else
      mcheck(.false.,'Unknown linear solver type: '//this%test_params%get_lin_solver_type())
    end if
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

  subroutine fix_pressure(this)
    implicit none
    class(stokes_driver_t), intent(inout) :: this

    class(fe_cell_iterator_t), allocatable  :: fe
    type(cell_map_t)         , pointer      :: cell_map
    type(quadrature_t)       , pointer      :: quad
    integer(ip)              , pointer      :: fe_dofs(:)
    type(point_t)            , pointer      :: quad_coords(:)
    class(scalar_function_t) , pointer      :: exact_pressure
    logical                  , allocatable  :: is_pressure_dof(:)
    class(vector_t)          , pointer      :: free_dof_values
    real(rp)                 , pointer      :: free_dof_entries(:)
    type(fe_cell_function_scalar_t)         :: fe_pressure

    integer(ip)  :: qpoint, num_quad_points
    real(rp)     :: dV, V
    real(rp) :: p_exact_gp
    real(rp) :: p_fe_gp
    real(rp) :: correction
    integer(ip) :: idof_p
    integer(ip) :: gdof

    call memalloc(this%fe_space%get_total_num_dofs(),is_pressure_dof,__FILE__,__LINE__)
    is_pressure_dof(:) = .false.

    exact_pressure => this%analytical_functions%get_solution_function_p()
    call fe_pressure%create(this%fe_space,P_FIELD_ID)
    call this%fe_space%create_fe_cell_iterator(fe)
    correction = 0.0
    V = 0.0
    do while ( .not. fe%has_finished() )

       call fe%update_integration()
       quad            => fe%get_quadrature()
       num_quad_points =  quad%get_num_quadrature_points()
       cell_map        => fe%get_cell_map()
       quad_coords     => fe%get_quadrature_points_coordinates()

       call fe_pressure%update(fe,this%solution)
       call fe%get_field_fe_dofs(P_FIELD_ID,fe_dofs)

       do idof_p = 1, fe%get_num_dofs_field(P_FIELD_ID)
         gdof = fe_dofs(idof_p)
         if (gdof>0) then
           is_pressure_dof(gdof) = .true.
         end if
       end do

       do qpoint = 1, num_quad_points
         dV = cell_map%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
         call exact_pressure%get_value(quad_coords(qpoint),p_exact_gp)
         call fe_pressure%get_value(qpoint,p_fe_gp)
         correction = correction + ( p_fe_gp  - p_exact_gp )*dV
         V = V + dV
       end do

       call fe%next()
    end do
    call this%fe_space%free_fe_cell_iterator(fe)
    call fe_pressure%free()

    correction = -correction / V

    free_dof_values => this%solution%get_free_dof_values()
    select type (free_dof_values)
      class is (serial_scalar_array_t)
        free_dof_entries => free_dof_values%get_entries()
      class DEFAULT
        check(.false.)
    end select

    do gdof = 1, this%fe_space%get_total_num_dofs()
      if (is_pressure_dof(gdof)) then
        free_dof_entries(gdof) = free_dof_entries(gdof) + correction
      end if
    end do

    call memfree(is_pressure_dof,__FILE__,__LINE__)

  end subroutine fix_pressure

  subroutine check_solution(this)
    implicit none
    class(stokes_driver_t), intent(inout) :: this
    !TODO do it for stokes

    type(unfitted_solution_checker_vector_t) :: solution_checker_u
    type(unfitted_solution_checker_t) :: solution_checker_p

    real(rp) :: error_h1_semi_norm
    real(rp) :: error_l2_norm
    real(rp) :: h1_semi_norm
    real(rp) :: l2_norm

    real(rp) :: l2_norm_boundary           
    real(rp) :: h1_semi_norm_boundary      
    real(rp) :: error_l2_norm_boundary     
    real(rp) :: error_h1_semi_norm_boundary

    real(rp) :: error_tolerance, tol
    integer(ip) :: iounit

#ifdef ENABLE_MKL
    if (this%test_params%get_lin_solver_type()=='pardiso') then
      error_tolerance = 1.0e-08
    else if (this%test_params%get_lin_solver_type()=='minres') then
      error_tolerance = 1.0e-04
    else
      mcheck(.false.,'Unknown linear solver type: '//this%test_params%get_lin_solver_type())
    end if
#else
    error_tolerance = 1.0e-04
#endif

    call solution_checker_u%create(this%fe_space,this%solution,this%analytical_functions%get_solution_function_u(),U_FIELD_ID)
    call solution_checker_u%compute_error_norms(error_h1_semi_norm,error_l2_norm,h1_semi_norm,l2_norm,&
           error_h1_semi_norm_boundary, error_l2_norm_boundary, h1_semi_norm_boundary, l2_norm_boundary)
    call solution_checker_u%free()

    write(*,'(a,e32.25)') 'u: l2_norm:               ', l2_norm
    write(*,'(a,e32.25)') 'u: h1_semi_norm:          ', h1_semi_norm
    write(*,'(a,e32.25)') 'u: error_l2_norm:         ', error_l2_norm
    write(*,'(a,e32.25)') 'u: error_h1_semi_norm:    ', error_h1_semi_norm
    write(*,'(a,e32.25)') 'u: rel_error_l2_norm:     ', error_l2_norm/l2_norm
    write(*,'(a,e32.25)') 'u: rel_error_h1_semi_norm:', error_h1_semi_norm/h1_semi_norm

    write(*,'(a,e32.25)') 'u: l2_norm_boundary:               ', l2_norm_boundary               
    write(*,'(a,e32.25)') 'u: h1_semi_norm_boundary:          ', h1_semi_norm_boundary          
    write(*,'(a,e32.25)') 'u: error_l2_norm_boundary:         ', error_l2_norm_boundary         
    write(*,'(a,e32.25)') 'u: error_h1_semi_norm_boundary:    ', error_h1_semi_norm_boundary    
    write(*,'(a,e32.25)') 'u: rel_error_l2_norm_boundary:     ', error_l2_norm_boundary      /l2_norm_boundary
    write(*,'(a,e32.25)') 'u: rel_error_h1_semi_norm_boundary:', error_h1_semi_norm_boundary /h1_semi_norm_boundary

    if (this%test_params%get_write_error_norms()) then
      iounit = io_open(file=this%test_params%get_dir_path_out()//this%test_params%get_prefix()//'_error_norms_u.csv',action='write')
      check(iounit>0)
      write(iounit,'(a,e32.25)') 'u: l2_norm                ;', l2_norm
      write(iounit,'(a,e32.25)') 'u: h1_semi_norm           ;', h1_semi_norm
      write(iounit,'(a,e32.25)') 'u: error_l2_norm          ;', error_l2_norm
      write(iounit,'(a,e32.25)') 'u: error_h1_semi_norm     ;', error_h1_semi_norm
      write(iounit,'(a,e32.25)') 'u: rel_error_l2_norm      ;', error_l2_norm/l2_norm
      write(iounit,'(a,e32.25)') 'u: rel_error_h1_semi_norm ;', error_h1_semi_norm/h1_semi_norm
      write(iounit,'(a,e32.25)') 'u: l2_norm_boundary               ;', l2_norm_boundary               
      write(iounit,'(a,e32.25)') 'u: h1_semi_norm_boundary          ;', h1_semi_norm_boundary          
      write(iounit,'(a,e32.25)') 'u: error_l2_norm_boundary         ;', error_l2_norm_boundary         
      write(iounit,'(a,e32.25)') 'u: error_h1_semi_norm_boundary    ;', error_h1_semi_norm_boundary    
      write(iounit,'(a,e32.25)') 'u: rel_error_l2_norm_boundary     ;', error_l2_norm_boundary      /l2_norm_boundary
      write(iounit,'(a,e32.25)') 'u: rel_error_h1_semi_norm_boundary;', error_h1_semi_norm_boundary /h1_semi_norm_boundary
      call io_close(iounit)
    end if

    if ( this%test_params%are_checks_active() ) then
      if (l2_norm == 0.0) then
        check( error_l2_norm < error_tolerance )
      else
        tol = error_tolerance*l2_norm
        check( error_l2_norm < tol )
      end if
      if (h1_semi_norm == 0.0) then
        check( error_h1_semi_norm < error_tolerance )
      else
        tol = error_tolerance*h1_semi_norm
        check( error_h1_semi_norm < tol )
      end if
    end if

    call solution_checker_p%create(this%fe_space,this%solution,this%analytical_functions%get_solution_function_p(),P_FIELD_ID)
    call solution_checker_p%compute_error_norms(error_h1_semi_norm,error_l2_norm,h1_semi_norm,l2_norm,&
           error_h1_semi_norm_boundary, error_l2_norm_boundary, h1_semi_norm_boundary, l2_norm_boundary)
    call solution_checker_p%free()

    write(*,'(a,e32.25)') 'p: l2_norm:               ', l2_norm
    write(*,'(a,e32.25)') 'p: h1_semi_norm:          ', h1_semi_norm
    write(*,'(a,e32.25)') 'p: error_l2_norm:         ', error_l2_norm
    write(*,'(a,e32.25)') 'p: error_h1_semi_norm:    ', error_h1_semi_norm
    write(*,'(a,e32.25)') 'p: rel_error_l2_norm:     ', error_l2_norm/l2_norm
    write(*,'(a,e32.25)') 'p: rel_error_h1_semi_norm:', error_h1_semi_norm/h1_semi_norm

    write(*,'(a,e32.25)') 'p: l2_norm_boundary:               ', l2_norm_boundary               
    write(*,'(a,e32.25)') 'p: h1_semi_norm_boundary:          ', h1_semi_norm_boundary          
    write(*,'(a,e32.25)') 'p: error_l2_norm_boundary:         ', error_l2_norm_boundary         
    write(*,'(a,e32.25)') 'p: error_h1_semi_norm_boundary:    ', error_h1_semi_norm_boundary    
    write(*,'(a,e32.25)') 'p: rel_error_l2_norm_boundary:     ', error_l2_norm_boundary      /l2_norm_boundary
    write(*,'(a,e32.25)') 'p: rel_error_h1_semi_norm_boundary:', error_h1_semi_norm_boundary /h1_semi_norm_boundary

    if (this%test_params%get_write_error_norms()) then
      iounit = io_open(file=this%test_params%get_dir_path_out()//this%test_params%get_prefix()//'_error_norms_p.csv',action='write')
      check(iounit>0)
      write(iounit,'(a,e32.25)') 'p: l2_norm                ;', l2_norm
      write(iounit,'(a,e32.25)') 'p: h1_semi_norm           ;', h1_semi_norm
      write(iounit,'(a,e32.25)') 'p: error_l2_norm          ;', error_l2_norm
      write(iounit,'(a,e32.25)') 'p: error_h1_semi_norm     ;', error_h1_semi_norm
      write(iounit,'(a,e32.25)') 'p: rel_error_l2_norm      ;', error_l2_norm/l2_norm
      write(iounit,'(a,e32.25)') 'p: rel_error_h1_semi_norm ;', error_h1_semi_norm/h1_semi_norm
      write(iounit,'(a,e32.25)') 'p: l2_norm_boundary               ;', l2_norm_boundary               
      write(iounit,'(a,e32.25)') 'p: h1_semi_norm_boundary          ;', h1_semi_norm_boundary          
      write(iounit,'(a,e32.25)') 'p: error_l2_norm_boundary         ;', error_l2_norm_boundary         
      write(iounit,'(a,e32.25)') 'p: error_h1_semi_norm_boundary    ;', error_h1_semi_norm_boundary    
      write(iounit,'(a,e32.25)') 'p: rel_error_l2_norm_boundary     ;', error_l2_norm_boundary      /l2_norm_boundary
      write(iounit,'(a,e32.25)') 'p: rel_error_h1_semi_norm_boundary;', error_h1_semi_norm_boundary /h1_semi_norm_boundary
      call io_close(iounit)
    end if

    if ( this%test_params%are_checks_active() ) then
      if (l2_norm == 0.0) then
        check( error_l2_norm < error_tolerance )
      else
        tol = error_tolerance*l2_norm
        check( error_l2_norm < tol )
      end if
      if (h1_semi_norm == 0.0) then
        check( error_h1_semi_norm < error_tolerance )
      else
        tol = error_tolerance*h1_semi_norm
        check( error_h1_semi_norm < tol )
      end if
    end if

  end subroutine check_solution
  
  subroutine write_solution(this)
    implicit none
    class(stokes_driver_t), intent(in) :: this
    type(output_handler_t)                   :: oh
    character(len=:), allocatable            :: path
    character(len=:), allocatable            :: prefix
    real(rp),allocatable :: cell_vector(:)
    real(rp),allocatable :: cell_vector_set_ids(:)
    real(rp), allocatable :: cell_rel_pos(:)
    real(rp), allocatable :: cell_in_aggregate(:)
    integer(ip) :: N, P, pid, i
    class(cell_iterator_t), allocatable :: cell
    
    real(rp),allocatable :: aggrs_ids(:)
    real(rp),allocatable :: aggrs_ids_color(:)
    integer(ip), pointer :: aggregate_ids(:)
    integer(ip), allocatable :: aggregate_ids_color(:)
    
    type(unfitted_vtk_writer_t) :: vtk_writer
    real(rp), pointer :: aggregate_size_ptr(:)
    real(rp), allocatable :: aggregate_size(:)

    if(this%test_params%get_write_solution()) then
        path = this%test_params%get_dir_path_out()
        prefix = this%test_params%get_prefix()
        call oh%create()
        call oh%attach_fe_space(this%fe_space)
        call oh%add_fe_function(this%solution, U_FIELD_ID, 'u')
        call oh%add_fe_function(this%solution, U_FIELD_ID, 'grad_u', grad_diff_operator)
        call oh%add_fe_function(this%solution, P_FIELD_ID, 'p')
        call memalloc(this%triangulation%get_num_cells(),cell_vector,__FILE__,__LINE__)
        call memalloc(this%triangulation%get_num_cells(),cell_vector_set_ids,__FILE__,__LINE__)
        call memalloc(this%triangulation%get_num_cells(),cell_rel_pos,__FILE__,__LINE__)
        call memalloc(this%triangulation%get_num_cells(),cell_in_aggregate,__FILE__,__LINE__)
        call memalloc(this%triangulation%get_num_cells(),aggrs_ids,__FILE__,__LINE__)
        call memalloc(this%triangulation%get_num_cells(),aggrs_ids_color,__FILE__,__LINE__)
        call memalloc(this%triangulation%get_num_cells(),aggregate_ids_color,__FILE__,__LINE__)
        call memalloc(this%triangulation%get_num_cells(),aggregate_size,__FILE__,__LINE__)
        
        if (this%test_params%get_use_constraints()) then
          aggregate_ids => this%fe_space%get_aggregate_ids()
          aggrs_ids(:) = real(aggregate_ids,kind=rp)
          aggregate_ids_color(:) = aggregate_ids
          call colorize_aggregate_ids(this%triangulation,aggregate_ids_color)
          aggrs_ids_color(:) = real(aggregate_ids_color,kind=rp)
        end if
        
        N=this%triangulation%get_num_cells()
        P=6
        call this%triangulation%create_cell_iterator(cell)
        do pid=0, P-1
            i=0
            do while ( i < (N*(pid+1))/P - (N*pid)/P ) 
              cell_vector(cell%get_gid()) = pid 
              call cell%next()
              i=i+1
            end do
        end do
        call this%triangulation%free_cell_iterator(cell)

        cell_rel_pos(:) = 0.0_rp
        call this%triangulation%create_cell_iterator(cell)
        do while (.not. cell%has_finished())
          cell_vector_set_ids(cell%get_gid()) = cell%get_set_id()
          if (cell%is_cut()) then
            cell_rel_pos(cell%get_gid()) = 0.0_rp
          else if (cell%is_interior()) then
            cell_rel_pos(cell%get_gid()) = -1.0_rp
          else if (cell%is_exterior()) then
            cell_rel_pos(cell%get_gid()) = 1.0_rp
          else
            mcheck(.false.,'Cell can only be either interior, exterior or cut')
          end if
          call cell%next()
        end do
        
        if (this%test_params%get_use_constraints()) then
          cell_in_aggregate(:) = 0.0_rp
          call cell%first()
          do while (.not. cell%has_finished())
            if (cell%is_cut()) then
              cell_in_aggregate(cell%get_gid()) = 1.0_rp
              cell_in_aggregate(aggregate_ids(cell%get_gid())) = 1.0_rp
            end if
            call cell%next()
          end do
        end if

        call this%triangulation%free_cell_iterator(cell)

        call oh%add_cell_vector(cell_vector,'cell_ids')
        call oh%add_cell_vector(cell_vector_set_ids,'cell_set_ids')
        call oh%add_cell_vector(cell_rel_pos,'cell_rel_pos')
        
        if (this%test_params%get_use_constraints()) then
          call oh%add_cell_vector(cell_in_aggregate,'cell_in_aggregate')
          aggregate_size_ptr => this%fe_space%get_aggregate_size()
          aggregate_size(:) = aggregate_size_ptr(:)
          call oh%add_cell_vector(aggregate_size,'aggregate_size')
        
          call oh%add_cell_vector(aggrs_ids,'aggregate_ids')
          call oh%add_cell_vector(aggrs_ids_color,'aggregate_ids_color')
        end if

        call oh%open(path, prefix)
        call oh%write()
        call oh%close()
        call oh%free()
        call memfree(cell_vector,__FILE__,__LINE__)
        call memfree(cell_vector_set_ids,__FILE__,__LINE__)
        call memfree(cell_rel_pos,__FILE__,__LINE__)
        call memfree(cell_in_aggregate,__FILE__,__LINE__)
        call memfree(aggrs_ids,__FILE__,__LINE__)
        call memfree(aggrs_ids_color,__FILE__,__LINE__)
        call memfree(aggregate_ids_color,__FILE__,__LINE__)
        call memfree(aggregate_size,__FILE__,__LINE__)

        ! Write the unfitted mesh
        call vtk_writer%attach_triangulation(this%triangulation)
        call vtk_writer%write_to_vtk_file(this%test_params%get_dir_path_out()//this%test_params%get_prefix()//'_mesh.vtu')
        call vtk_writer%free()

        ! TODO make it work when no cut cells
        ! Write the unfitted mesh
        call vtk_writer%attach_boundary_faces(this%triangulation)
        call vtk_writer%write_to_vtk_file(this%test_params%get_dir_path_out()//this%test_params%get_prefix()//'_boundary_faces.vtu')
        call vtk_writer%free()
        
        ! TODO do it for stokes
        !! Write the solution
        !call vtk_writer%attach_fe_function(this%solution,this%fe_space)
        !call vtk_writer%write_to_vtk_file(this%test_params%get_dir_path_out()//this%test_params%get_prefix()//'_mesh_solution.vtu')
        !call vtk_writer%free()
        

    endif
  end subroutine write_solution
  
  subroutine run_simulation(this) 
    implicit none
    class(stokes_driver_t), intent(inout) :: this
    call this%free()
    call this%parse_command_line_parameters()
    call this%setup_levelset()
    call this%setup_triangulation()
    call this%fill_cells_set()
    call this%setup_reference_fes()
    call this%setup_fe_space()
    
    
    call this%setup_system()

    if ( .not. this%test_params%get_only_setup() ) then
      call this%assemble_system()
      call this%setup_solver()
      call this%solve_system()
      call this%fix_pressure()
      call this%check_solution()
    end if

    call this%write_solution()
    call this%free()
  end subroutine run_simulation
  
  subroutine free(this)
    implicit none
    class(stokes_driver_t), intent(inout) :: this
    integer(ip) :: i, istat
    
    call this%solution%free()
    
#ifdef ENABLE_MKL        
      call this%direct_solver%free()
      call this%iterative_linear_solver%free()
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
  

  
end module stokes_driver_names
