module test_geometry_params
    use fempar_names
# include "debug.i90"

implicit none
private

    type,extends(parameter_generator_t) :: test_geometry_params_t
        character(len=:),   allocatable :: default_dir_path 
        character(len=:),   allocatable :: default_prefix
        character(len=:),   allocatable :: default_dir_path_out
        character(len=str_cla_len)      :: dir_path
        character(len=str_cla_len)      :: prefix
        character(len=str_cla_len)      :: dir_path_out
    contains
        procedure :: set_default      => test_geometry_params_set_default
        procedure :: get_dir_path     => test_geometry_params_get_dir_path
        procedure :: get_prefix       => test_geometry_params_get_prefix
        procedure :: get_dir_path_out => test_geometry_params_get_dir_path_out
    end type test_geometry_params_t

  ! Types
  public :: test_geometry_params_t

contains


  !==================================================================================================
  subroutine test_geometry_params_set_default(this)
    implicit none
    class(test_geometry_params_t), intent(inout) :: this
    type(ParameterList_t), pointer :: list, switches, switches_ab, helpers, required
    integer(ip)    :: error
    character(len=:), allocatable            :: msg

    list        => this%get_parameters()
    switches    => this%get_switches()
    switches_ab => this%get_switches_ab()
    helpers     => this%get_helpers()
    required    => this%get_required()

    error = list%set(key = dir_path_key,     value = '.') ;      check(error==0)
    error = list%set(key = prefix_key,       value = 'square') ; check(error==0)
    error = list%set(key = dir_path_out_key, value = '.') ;      check(error==0)

    ! Only some of them are controlled from cli
    error = switches%set(key = dir_path_key,     value = '--dir-path') ;     check(error==0)
    error = switches%set(key = prefix_key,       value = '--prefix') ;       check(error==0)
    error = switches%set(key = dir_path_out_key, value = '--dir-path-out') ; check(error==0)
                                                             
    error = switches_ab%set(key = dir_path_key,     value = '-d') ; check(error==0) 
    error = switches_ab%set(key = prefix_key,       value = '-p') ; check(error==0) 
    error = switches_ab%set(key = dir_path_out_key, value = '-o') ; check(error==0) 

    error = helpers%set(key = dir_path_key,     value = 'Directory of the source files') ; check(error==0)
    error = helpers%set(key = prefix_key,       value = 'Name of the GiD files') ;         check(error==0)
    error = helpers%set(key = dir_path_out_key, value = 'Output Directory') ;              check(error==0)

    
    error = required%set(key = dir_path_key,     value = .false.) ; check(error==0)
    error = required%set(key = prefix_key,       value = .false.) ; check(error==0)
    error = required%set(key = dir_path_out_key, value = .false.) ; check(error==0)


  end subroutine test_geometry_params_set_default


  function test_geometry_params_get_dir_path(this)
    implicit none
    class(test_geometry_params_t) , intent(in) :: this
    character(len=:),      allocatable         :: test_geometry_params_get_dir_path
    type(ParameterList_t), pointer             :: list
    integer(ip)                                :: error
    list  => this%get_parameters()
    assert(list%isAssignable(dir_path_key, 'string'))
    error = list%GetAsString(key = dir_path_key, string = test_geometry_params_get_dir_path)
    assert(error==0)
  end function test_geometry_params_get_dir_path


  function test_geometry_params_get_prefix(this)
    implicit none
    class(test_geometry_params_t) , intent(in) :: this
    character(len=:),      allocatable         :: test_geometry_params_get_prefix
    type(ParameterList_t), pointer             :: list
    integer(ip)                                :: error
    list  => this%get_parameters()
    assert(list%isAssignable(prefix_key, 'string'))
    error = list%GetAsString(key = prefix_key, string = test_geometry_params_get_prefix)
    assert(error==0)
  end function test_geometry_params_get_prefix


  function test_geometry_params_get_dir_path_out(this)
    implicit none
    class(test_geometry_params_t) , intent(in) :: this
    character(len=:),      allocatable         :: test_geometry_params_get_dir_path_out
    type(ParameterList_t), pointer             :: list
    integer(ip)                                :: error
    list  => this%get_parameters()
    assert(list%isAssignable(dir_path_out_key, 'string'))
    error = list%GetAsString(key = dir_path_out_key, string = test_geometry_params_get_dir_path_out)
    assert(error==0)
  end function test_geometry_params_get_dir_path_out

end module test_geometry_params

