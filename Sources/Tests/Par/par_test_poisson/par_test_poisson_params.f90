module par_test_poisson_params_names
  use serial_names
  implicit none
#include "debug.i90" 
  private
  
  ! Types
  type par_test_poisson_params_t
     private
  
     ! IO parameters
     character(len=:), allocatable :: default_dir_path
     character(len=:), allocatable :: default_prefix
     character(len=:), allocatable :: default_nparts
     character(len=:), allocatable :: default_reference_fe_geo_order
     character(len=:), allocatable :: default_reference_fe_order

     ! IO parameters
     character(len=256)            :: dir_path
     character(len=256)            :: prefix
     integer(ip)                   :: nparts
     integer(ip)                   :: reference_fe_geo_order
     integer(ip)                   :: reference_fe_order
     
 
     type(Command_Line_Interface)  :: cli
  contains
     procedure, non_overridable             :: create       => par_test_poisson_params_create
     procedure, non_overridable, private    :: set_default  => par_test_poisson_params_set_default
     procedure, non_overridable, private    :: add_to_cli   => par_test_poisson_params_add_to_cli
     procedure, non_overridable             :: parse        => par_test_poisson_params_parse 
     procedure, non_overridable             :: free         => par_test_poisson_params_free
     procedure, non_overridable             :: get_dir_path
     procedure, non_overridable             :: get_prefix
     procedure, non_overridable             :: get_nparts
     procedure, non_overridable             :: get_reference_fe_geo_order
     procedure, non_overridable             :: get_reference_fe_order
  end type par_test_poisson_params_t

  ! Types
  public :: par_test_poisson_params_t

contains
  subroutine par_test_poisson_params_create ( this )
    implicit none
    class(par_test_poisson_params_t), intent(inout) :: this
     
    call this%free()
    
    ! Initialize Command Line Interface
    call this%cli%init(progname    = 'par_test_poisson', &
                       version     = '', &
                       authors     = '', &
                       license     = '', &
                       description =  'FEMPAR parallel test to solve the 2D Poisson PDE with known analytical solution. &
                                       Boundary set ID 1 MUST BE ASSIGNED to the whole boundary.', & 
                       examples    = ['par_test_poisson -h'] )
    call this%set_default()
    call this%add_to_cli()
  end subroutine par_test_poisson_params_create 

  !==================================================================================================
  subroutine par_test_poisson_params_set_default(this)
    implicit none
    class(par_test_poisson_params_t), intent(inout) :: this
    this%default_dir_path     = 'PARTS_4/'
    this%default_prefix       = 'square'
    this%default_nparts       = '4'
    this%default_reference_fe_geo_order = '1'
    this%default_reference_fe_order = '1'
  end subroutine par_test_poisson_params_set_default

  !==================================================================================================
  subroutine par_test_poisson_params_add_to_cli(this)
    implicit none
    class(par_test_poisson_params_t), intent(inout) :: this
    integer(ip) :: error

    ! Set Command Line Arguments
    call this%cli%add(switch='--dir-path',switch_ab='-d',help='Absolute or relative PATH to the partitioned&
                       & problem. Must end with /',required=.false., act='store', def=trim(this%default_dir_path), error=error)
    check(error==0)

    call this%cli%add(switch='--prefix',switch_ab='-p',help='Prefix for all input files (mesh, partition, etc.).& 
                       & E.g., if these files were generated from square.gid GiD project, then --prefix square.',& 
                       & required=.false., act='store', def=trim(this%default_prefix), error=error)
    check(error==0)

    call this%cli%add(switch='--nparts',switch_ab='-n',help='Number of parts in which the problem was split.',& 
                       & required=.false., act='store', def=trim(this%default_nparts), error=error)
    check(error==0)
    call this%cli%add(switch='--reference-fe-geo-order',switch_ab='-gorder',help='Order of the triangulation reference fe',&
         &            required=.false.,act='store',def=trim(this%default_reference_fe_geo_order),error=error)
    check(error==0)  
    call this%cli%add(switch='--reference-fe-order',switch_ab='-order',help='Order of the fe space reference fe',&
         &            required=.false.,act='store',def=trim(this%default_reference_fe_order),error=error) 
    check(error==0) 
    
  end subroutine par_test_poisson_params_add_to_cli

  !==================================================================================================
  subroutine par_test_poisson_params_parse(this)
    implicit none
    class(par_test_poisson_params_t), intent(inout) :: this
    integer(ip) :: error
    call this%cli%parse(error=error)
    check(error==0)
    call this%cli%get(switch='-d',val=this%dir_path,error=error); check(error==0)
    call this%cli%get(switch='-p',val=this%prefix,error=error); check(error==0)
    call this%cli%get(switch='-n',val=this%nparts,error=error); check(error==0)
    call this%cli%get(switch='-gorder',val=this%reference_fe_geo_order,error=error); check(error==0)
    call this%cli%get(switch='-order',val=this%reference_fe_order,error=error); check(error==0)
  end subroutine par_test_poisson_params_parse
  
  subroutine par_test_poisson_params_free(this)
    implicit none
    class(par_test_poisson_params_t), intent(inout) :: this
    if(allocated(this%default_dir_path)) deallocate(this%default_dir_path)              
    if(allocated(this%default_prefix)) deallocate(this%default_prefix)                    
    if(allocated(this%default_nparts)) deallocate(this%default_nparts)
    if(allocated(this%default_reference_fe_geo_order)) deallocate(this%default_reference_fe_geo_order)
    if(allocated(this%default_reference_fe_order)) deallocate(this%default_reference_fe_order)
    call this%cli%free()
  end subroutine par_test_poisson_params_free
  
  ! GETTERS *****************************************************************************************
  function get_dir_path(this)
    implicit none
    class(par_test_poisson_params_t) , intent(in) :: this
    character(len=256) :: get_dir_path
    get_dir_path = this%dir_path
  end function get_dir_path

  !==================================================================================================
  function get_prefix(this)
    implicit none
    class(par_test_poisson_params_t) , intent(in) :: this
    character(len=256) :: get_prefix
    get_prefix = this%prefix
  end function get_prefix

  !==================================================================================================
  function get_nparts(this)
    implicit none
    class(par_test_poisson_params_t) , intent(in) :: this
    integer(ip) :: get_nparts
    get_nparts = this%nparts
  end function get_nparts
  
    !==================================================================================================
  function get_reference_fe_geo_order(this)
    implicit none
    class(par_test_poisson_params_t) , intent(in) :: this
    integer(ip) :: get_reference_fe_geo_order
    get_reference_fe_geo_order = this%reference_fe_geo_order
  end function get_reference_fe_geo_order
  
  !==================================================================================================
  function get_reference_fe_order(this)
    implicit none
    class(par_test_poisson_params_t) , intent(in) :: this
    integer(ip) :: get_reference_fe_order
    get_reference_fe_order = this%reference_fe_order
  end function get_reference_fe_order
  
  
end module par_test_poisson_params_names
