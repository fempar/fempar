module par_test_poisson_unfitted_params_names
  use fempar_names
  use level_set_functions_gallery_names

  implicit none
#include "debug.i90" 
  private

  character(len=*), parameter :: reference_fe_geo_order_key = 'reference_fe_geo_order'
  character(len=*), parameter :: reference_fe_order_key     = 'reference_fe_order'    
  character(len=*), parameter :: write_solution_key         = 'write_solution'        
  character(len=*), parameter :: triangulation_type_key     = 'triangulation_type'    
  character(len=*), parameter :: coarse_fe_handler_type_key = 'coarse_fe_handler_type'    
  character(len=*), parameter :: level_set_function_type_key= 'level_set_function_type'    
  character(len=*), parameter :: use_preconditioner_key     = 'use_preconditioner'    
  character(len=*), parameter :: unfitted_boundary_type_key = 'unfitted_boundary_type'    
  character(len=*), parameter :: nitsche_beta_factor_key    = 'nitsche_beta_factor'    
  character(len=*), parameter :: levelset_tolerance_key     = 'levelset_tolerance'    
  character(len=*), parameter :: num_runs_key            = 'num_runs'    
  character(len=*), parameter :: is_in_fe_space_key         = 'is_in_fe_space'    
  character(len=*), parameter :: are_checks_active_key      = 'are_checks_active'    

  character(len=*), public, parameter :: unfitted_coarse_fe_handler_value = 'unfitted'    
  character(len=*), public, parameter :: standard_coarse_fe_handler_value = 'standard'    
  character(len=*), public, parameter :: stiffness_coarse_fe_handler_value = 'stiffness'

  type, extends(parameter_handler_t) :: par_test_poisson_unfitted_params_t
     private
     contains
       procedure :: define_parameters  => par_test_poisson_unfitted_params_define_parameters
       procedure, non_overridable             :: get_dir_path
       procedure, non_overridable             :: get_prefix
       procedure, non_overridable             :: get_reference_fe_geo_order
       procedure, non_overridable             :: get_reference_fe_order
       procedure, non_overridable             :: get_write_solution
       procedure, non_overridable             :: get_triangulation_type
       procedure, non_overridable             :: get_coarse_fe_handler_type
       procedure, non_overridable             :: get_level_set_function_type
       procedure, non_overridable             :: get_use_preconditioner
       procedure, non_overridable             :: get_unfitted_boundary_type
       procedure, non_overridable             :: get_nitsche_beta_factor
       procedure, non_overridable             :: get_levelset_tolerance
       procedure, non_overridable             :: get_num_runs
       procedure, non_overridable             :: print
       procedure, non_overridable, private    :: print_character_switch
       procedure, non_overridable, private    :: print_integer_switch
       procedure, non_overridable, private    :: print_real_switch
       procedure, non_overridable, private    :: print_logical_switch
       !procedure, non_overridable             :: get_num_dims
       procedure, non_overridable             :: is_in_fe_space
       procedure, non_overridable             :: are_checks_active
  end type par_test_poisson_unfitted_params_t

  ! Types
  public :: par_test_poisson_unfitted_params_t

contains

  !==================================================================================================
  subroutine par_test_poisson_unfitted_params_define_parameters(this)
    implicit none
    class(par_test_poisson_unfitted_params_t), intent(inout) :: this
    character(len=:), allocatable                            :: msg

    msg = 'structured (*) or unstructured (*) triangulation?'
    write(msg(13:13),'(i1)') triangulation_generate_structured
    write(msg(33:33),'(i1)') triangulation_generate_from_mesh

    ! common
    call this%add(dir_path_key,'--dir-path', '.', 'Directory of the source files', switch_ab='-d')
    call this%add(prefix_key, '--prefix', 'square', 'Name of the GiD files', switch_ab='-p')
    call this%add(dir_path_out_key, '--dir-path-out', '.', 'Output Directory', switch_ab='-o')
    call this%add(struct_hex_triang_num_dims_key, '--dim', 2, 'Number of space dimensions', switch_ab='-dm')
    call this%add(struct_hex_triang_num_cells_dir, '--num_cells', [12,12,12], 'Number of cells per dir', switch_ab='-n')
    call this%add(struct_hex_triang_is_dir_periodic_key, '--STRUCT_HEX_TRIANG_IS_DIR_PERIODIC', [0,0,0], 'Is the mesh periodic for every dimension')           
    call this%add(struct_hex_triang_num_levels_key, '--num_levels', 3, 'Number of levels', switch_ab='-l')
    call this%add(struct_hex_triang_num_parts_x_dir_key, '--num_parts_x_dir', [4,4,0,2,2,0,1,1,0], 'Number of parts per dir and per level', switch_ab='-np')
    call this%add(reference_fe_geo_order_key, '--reference-fe-geo-order', 1, 'Order of the triangulation reference fe', switch_ab='-gorder')
    call this%add(reference_fe_order_key, '--reference-fe-order', 1, 'Order of the fe space reference fe', switch_ab='-order')
    call this%add(write_solution_key, '--write-solution', .false., 'Write solution in VTK format', switch_ab='-wsolution')
    call this%add(triang_generate_key, '--trinagulation-type', triangulation_generate_from_mesh, msg, switch_ab='-tt')
    call this%add(coarse_space_use_vertices_key, '--coarse-space-use-vertices', .true., 'Include vertex coarse DoFs in coarse FE space', switch_ab='-use-vertices')
    call this%add(coarse_space_use_edges_key, '--coarse-space-use-edges', .true., 'Include edge coarse DoFs in coarse FE space', switch_ab='-use-edges')
    call this%add(coarse_space_use_faces_key, '--coarse-space-use-faces', .true., 'Include face coarse DoFs in coarse FE space', switch_ab='-use-faces')

    ! specific
    call this%add(coarse_fe_handler_type_key, '--coarse-fe-handler', unfitted_coarse_fe_handler_value, 'Type of coarse fe handler. `'//&
      standard_coarse_fe_handler_value//'`, `'//&
      unfitted_coarse_fe_handler_value//'` or `'//&
      stiffness_coarse_fe_handler_value//'` ?', switch_ab='-chandler')
    call this%add(level_set_function_type_key, '--level-set-function', level_set_sphere_str, 'Type of levelset to be used.'//&
      ' The possible values are the public character constants defined in the `level_set_functions_gallery_names` module', switch_ab='-levelset')
    call this%add(use_preconditioner_key, '--use-preconditioner', .true., 'Use (T) or not (F) a preconditioner', switch_ab='-precond')
    call this%add(unfitted_boundary_type_key, '--unfitted-boundary', 'dirichlet', 'Use (dirichlet) or not (neumann) boundary conditions on the unfitted boundary', switch_ab='-uboundary')
    call this%add(nitsche_beta_factor_key, '--nitsche-beta', 2.0_rp, 'Set the value of the factor to compute nitches beta', switch_ab='-beta')
    call this%add(levelset_tolerance_key, '--level-set-tol',  1.0e-2_rp, 'Set the tolerance of the levelset', switch_ab='-levelsettol')
    call this%add(num_runs_key, '--number-runs', 1, 'Number of times the simulation is repeated (useful when measuring times)', switch_ab='-nruns')
    call this%add(is_in_fe_space_key, '--solution_in_fe_space', .true., 'Is the solution in fe space', switch_ab='-in_space')
    call this%add(are_checks_active_key, '--check_solution', .true., 'Check or not the solution', '-check')


    
  end subroutine par_test_poisson_unfitted_params_define_parameters

  ! GETTERS *****************************************************************************************
  function get_dir_path(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    character(len=:),      allocatable            :: get_dir_path
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(dir_path_key, 'string'))
    error = list%GetAsString(key = dir_path_key, string = get_dir_path)
    assert(error==0)
  end function get_dir_path

  !==================================================================================================
  function get_prefix(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    character(len=:),      allocatable            :: get_prefix
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(prefix_key, 'string'))
    error = list%GetAsString(key = prefix_key, string = get_prefix)
    assert(error==0)
  end function get_prefix

    !==================================================================================================
  function get_reference_fe_geo_order(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    integer(ip)                                   :: get_reference_fe_geo_order
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(reference_fe_geo_order_key, get_reference_fe_geo_order))
    error = list%Get(key = reference_fe_geo_order_key, Value = get_reference_fe_geo_order)
    assert(error==0)
  end function get_reference_fe_geo_order
  
  !==================================================================================================
  function get_reference_fe_order(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    integer(ip)                                   :: get_reference_fe_order
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(reference_fe_order_key, get_reference_fe_order))
    error = list%Get(key = reference_fe_order_key, Value = get_reference_fe_order)
    assert(error==0)
  end function get_reference_fe_order
  
  !==================================================================================================
  function get_write_solution(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    logical                                       :: get_write_solution
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    logical                                       :: is_present
    logical                                       :: same_data_type
    integer(ip), allocatable                      :: shape(:)
    list  => this%get_values()
    assert(list%isAssignable(write_solution_key, get_write_solution))
    error = list%Get(key = write_solution_key, Value = get_write_solution)
    assert(error==0)
  end function get_write_solution

  !==================================================================================================
  function get_triangulation_type(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    integer(ip)                                   :: get_triangulation_type
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(triang_generate_key, get_triangulation_type))
    error = list%Get(key = triang_generate_key, Value = get_triangulation_type)
    assert(error==0)
  end function get_triangulation_type 

  !==================================================================================================
  function get_coarse_fe_handler_type(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    character(len=:),      allocatable            :: get_coarse_fe_handler_type
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(coarse_fe_handler_type_key, 'string'))
    error = list%GetAsString(key = coarse_fe_handler_type_key, string = get_coarse_fe_handler_type)
    assert(error==0)
  end function get_coarse_fe_handler_type

  !==================================================================================================
  function get_level_set_function_type(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    character(len=:),      allocatable            :: get_level_set_function_type
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(level_set_function_type_key, 'string'))
    error = list%GetAsString(key = level_set_function_type_key, string = get_level_set_function_type)
    assert(error==0)
  end function get_level_set_function_type

  !==================================================================================================
  function get_use_preconditioner(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    logical                                       :: get_use_preconditioner
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(use_preconditioner_key, get_use_preconditioner))
    error = list%Get(key = use_preconditioner_key, Value = get_use_preconditioner)
    assert(error==0)
  end function get_use_preconditioner

  !==================================================================================================
  function get_unfitted_boundary_type(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    character(len=:),      allocatable            :: get_unfitted_boundary_type
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(level_set_function_type_key, 'string'))
    error = list%GetAsString(key = unfitted_boundary_type_key, string = get_unfitted_boundary_type)
    assert(error==0)
  end function get_unfitted_boundary_type

  !==================================================================================================
  function get_nitsche_beta_factor(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    real(rp)                                      :: get_nitsche_beta_factor
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(nitsche_beta_factor_key, get_nitsche_beta_factor))
    error = list%Get(key = nitsche_beta_factor_key, Value = get_nitsche_beta_factor)
    assert(error==0)
  end function get_nitsche_beta_factor

  !==================================================================================================
  function get_levelset_tolerance(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    real(rp)                                      :: get_levelset_tolerance
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(levelset_tolerance_key, get_levelset_tolerance))
    error = list%Get(key = levelset_tolerance_key, Value = get_levelset_tolerance)
    assert(error==0)
  end function get_levelset_tolerance

  !==================================================================================================
  function get_num_runs(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    integer(ip)                                   :: get_num_runs
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(num_runs_key, get_num_runs))
    error = list%Get(key = num_runs_key, Value = get_num_runs)
    assert(error==0)
  end function get_num_runs

  !==================================================================================================
  subroutine print(this, environment)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    class(environment_t) :: environment

    if (environment%get_l1_rank()==0) then
      call this%print_character_switch(coarse_fe_handler_type_key )
      call this%print_character_switch(level_set_function_type_key)
      call this%print_logical_switch  (use_preconditioner_key     )
      call this%print_character_switch(unfitted_boundary_type_key )
      call this%print_real_switch     (nitsche_beta_factor_key    )
      call this%print_real_switch     (levelset_tolerance_key     )
      call this%print_integer_switch  (num_runs_key            )
    end if

  end subroutine print

  !==================================================================================================
  subroutine print_character_switch(this,switch_key)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    character(len=*),                           intent(in) :: switch_key
    type(ParameterList_t), pointer :: values, switches, switches_ab
    integer(ip)                    :: error
    character(len=:), allocatable  :: switch, switch_ab
    character(len=:), allocatable  :: val
    character(len=30) :: charaux, charaux_ab
    values   => this%get_values()
    switches => this%get_switches()
    switches_ab => this%get_switches_ab()
    assert(switches%isAssignable(key = switch_key, value  = 'string'))
    error= switches%GetAsString (key = switch_key, string = switch  )
    assert(switches_ab%isAssignable(key = switch_key, value  = 'string'))
    error= switches_ab%GetAsString (key = switch_key, string = switch_ab)
    assert(values%isAssignable  (key = switch_key, value  = 'string'))
    error= values%GetAsString   (key = switch_key, string = val     )
    write (charaux, '(a30)') switch
    write (charaux_ab, '(a30)') switch_ab
    write(*,'(2a30,a20)') adjustl(charaux),adjustl(charaux_ab), val
  end subroutine print_character_switch

  !==================================================================================================
  subroutine print_integer_switch(this,switch_key)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    character(len=*),                           intent(in) :: switch_key
    type(ParameterList_t), pointer :: values, switches, switches_ab
    integer(ip)                    :: error
    character(len=:), allocatable  :: switch, switch_ab
    integer(ip)  :: val
    character(len=30) :: charaux, charaux_ab
    values   => this%get_values()
    switches => this%get_switches()
    switches_ab => this%get_switches_ab()
    assert(switches%isAssignable(key = switch_key, value  = 'string'))
    error= switches%GetAsString (key = switch_key, string = switch  )
    assert(switches_ab%isAssignable(key = switch_key, value  = 'string'))
    error= switches_ab%GetAsString (key = switch_key, string = switch_ab)
    assert(values%isAssignable  (key = switch_key, value  = val     ))
    error= values%Get           (key = switch_key, value  = val     )
    write (charaux, '(a30)') switch
    write (charaux_ab, '(a30)') switch_ab
    write(*,'(2a30,i20)') adjustl(charaux),adjustl(charaux_ab), val
  end subroutine print_integer_switch

  !==================================================================================================
  subroutine print_real_switch(this,switch_key)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    character(len=*),                           intent(in) :: switch_key
    type(ParameterList_t), pointer :: values, switches, switches_ab
    integer(ip)                    :: error
    character(len=:), allocatable  :: switch, switch_ab
    real(rp)  :: val
    character(len=30) :: charaux, charaux_ab
    values   => this%get_values()
    switches => this%get_switches()
    switches_ab => this%get_switches_ab()
    assert(switches%isAssignable(key = switch_key, value  = 'string'))
    error= switches%GetAsString (key = switch_key, string = switch  )
    assert(switches_ab%isAssignable(key = switch_key, value  = 'string'))
    error= switches_ab%GetAsString (key = switch_key, string = switch_ab)
    assert(values%isAssignable  (key = switch_key, value  = val     ))
    error= values%Get           (key = switch_key, value  = val     )
    write (charaux, '(a30)') switch
    write (charaux_ab, '(a30)') switch_ab
    write(*,'(2a30,e20.5)') adjustl(charaux),adjustl(charaux_ab), val
  end subroutine print_real_switch

  !==================================================================================================
  subroutine print_logical_switch(this,switch_key)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    character(len=*),                           intent(in) :: switch_key
    type(ParameterList_t), pointer :: values, switches, switches_ab
    integer(ip)                    :: error
    character(len=:), allocatable  :: switch, switch_ab
    logical  :: val
    character(len=30) :: charaux, charaux_ab
    values   => this%get_values()
    switches => this%get_switches()
    switches_ab => this%get_switches_ab()
    assert(switches%isAssignable(key = switch_key, value  = 'string'))
    error= switches%GetAsString (key = switch_key, string = switch  )
    assert(switches_ab%isAssignable(key = switch_key, value  = 'string'))
    error= switches_ab%GetAsString (key = switch_key, string = switch_ab)
    assert(values%isAssignable  (key = switch_key, value  = val     ))
    error= values%Get           (key = switch_key, value  = val     )
    write (charaux, '(a30)') switch
    write (charaux_ab, '(a30)') switch_ab
    write(*,'(2a30,l20)') adjustl(charaux),adjustl(charaux_ab), val
  end subroutine print_logical_switch

  !==================================================================================================
  function is_in_fe_space(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    logical :: is_in_fe_space
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(is_in_fe_space_key, is_in_fe_space))
    error = list%Get(key = is_in_fe_space_key, Value = is_in_fe_space)
    assert(error==0)
  end function is_in_fe_space

  !==================================================================================================
  function are_checks_active(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    logical :: are_checks_active
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(are_checks_active_key, are_checks_active))
    error = list%Get(key = are_checks_active_key, Value = are_checks_active)
    assert(error==0)
  end function are_checks_active

end module par_test_poisson_unfitted_params_names
