module par_test_reference_fe_params_names
  use serial_names
  implicit none
#include "debug.i90" 
  private
  
  ! Types
  type par_test_reference_fe_params_t
     private
  
     ! IO parameters
     character(len=:), allocatable :: default_dir_path
     character(len=:), allocatable :: default_prefix
     character(len=:), allocatable :: default_nparts

     ! IO parameters
     character(len=256)            :: dir_path
     character(len=256)            :: prefix
     integer(ip)                   :: nparts
 
     type(Command_Line_Interface)  :: cli
  contains
     procedure, non_overridable             :: create       => par_test_reference_fe_params_create
     procedure, non_overridable, private    :: set_default  => par_test_reference_fe_params_set_default
     procedure, non_overridable, private    :: add_to_cli   => par_test_reference_fe_params_add_to_cli
     procedure, non_overridable             :: parse        => par_test_reference_fe_params_parse 
     procedure, non_overridable             :: free         => par_test_reference_fe_params_free
     procedure, non_overridable             :: get_dir_path
     procedure, non_overridable             :: get_prefix
     procedure, non_overridable             :: get_nparts
  end type par_test_reference_fe_params_t

  ! Types
  public :: par_test_reference_fe_params_t

contains
  subroutine par_test_reference_fe_params_create ( this )
    implicit none
    class(par_test_reference_fe_params_t), intent(inout) :: this
     
    call this%free()
    
    ! Initialize Command Line Interface
    call this%cli%init(progname    = 'par_test_reference_fe', &
         &             version     = '', &
         &             authors     = '', &
         &             license     = '', &
         &             description = "FEMPAR parallel test driver", &
         &             examples    = ['par_test_reference_fe -h'] )
    call this%set_default()
    call this%add_to_cli()
  end subroutine par_test_reference_fe_params_create 

  !==================================================================================================
  subroutine par_test_reference_fe_params_set_default(this)
    implicit none
    class(par_test_reference_fe_params_t), intent(inout) :: this
    this%default_dir_path     = 'PARTS_4/'
    this%default_prefix       = 'square'
    this%default_nparts       = '4'
  end subroutine par_test_reference_fe_params_set_default

  !==================================================================================================
  subroutine par_test_reference_fe_params_add_to_cli(this)
    implicit none
    class(par_test_reference_fe_params_t), intent(inout) :: this
    integer(ip) :: error

    ! Set Command Line Arguments
    call this%cli%add(switch='--dir-path',switch_ab='-dir-path',help='Absolute or relative PATH to the partitioned&
                       & problem. Must end with /',required=.false., act='store', def=trim(this%default_dir_path), error=error)
    check(error==0)

    call this%cli%add(switch='--prefix',switch_ab='-prefix',help='Prefix for all input files (mesh, partition, etc.).& 
                       & E.g., if these files were generated from square.gid GiD project, then --prefix square.',& 
                       & required=.false., act='store', def=trim(this%default_prefix), error=error)
    check(error==0)

    call this%cli%add(switch='--nparts',switch_ab='-nparts',help='Number of parts in which the problem was split.',& 
                       & required=.false., act='store', def=trim(this%default_nparts), error=error)
    check(error==0)
  end subroutine par_test_reference_fe_params_add_to_cli

  !==================================================================================================
  subroutine par_test_reference_fe_params_parse(this)
    implicit none
    class(par_test_reference_fe_params_t), intent(inout) :: this
    integer(ip) :: error
    call this%cli%parse(error=error)
    check(error==0)
    call this%cli%get(switch='-dir-path',val=this%dir_path,error=error); check(error==0)
    call this%cli%get(switch='-prefix',val=this%prefix,error=error); check(error==0)
    call this%cli%get(switch='-nparts',val=this%nparts,error=error); check(error==0)
  end subroutine par_test_reference_fe_params_parse
  
  subroutine par_test_reference_fe_params_free(this)
    implicit none
    class(par_test_reference_fe_params_t), intent(inout) :: this
    if(allocated(this%default_dir_path)) deallocate(this%default_dir_path)              
    if(allocated(this%default_prefix)) deallocate(this%default_prefix)                    
    if(allocated(this%default_nparts)) deallocate(this%default_nparts)
    call this%cli%free()
  end subroutine par_test_reference_fe_params_free
  
  ! GETTERS *****************************************************************************************
  function get_dir_path(this)
    implicit none
    class(par_test_reference_fe_params_t) , intent(in) :: this
    character(len=256) :: get_dir_path
    get_dir_path = this%dir_path
  end function get_dir_path

  !==================================================================================================
  function get_prefix(this)
    implicit none
    class(par_test_reference_fe_params_t) , intent(in) :: this
    character(len=256) :: get_prefix
    get_prefix = this%prefix
  end function get_prefix

  !==================================================================================================
  function get_nparts(this)
    implicit none
    class(par_test_reference_fe_params_t) , intent(in) :: this
    integer(ip) :: get_nparts
    get_nparts = this%nparts
  end function get_nparts
  
end module par_test_reference_fe_params_names
