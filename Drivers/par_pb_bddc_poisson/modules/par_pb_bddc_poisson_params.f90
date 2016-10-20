module par_pb_bddc_poisson_params_names
  use fempar_names
  implicit none
#include "debug.i90" 
  private
  
  ! Types
  type par_pb_bddc_poisson_params_t
     private
  
     ! IO parameters
     character(len=:), allocatable :: default_dir_path
     character(len=:), allocatable :: default_prefix
     character(len=:), allocatable :: default_nparts
     character(len=:), allocatable :: default_reference_fe_geo_order
     character(len=:), allocatable :: default_reference_fe_order
     character(len=:), allocatable :: default_write_solution

     character(len=:), allocatable :: default_triangulation_type
     character(len=:), allocatable :: default_num_dimensions
     character(len=:), allocatable :: default_nx
     character(len=:), allocatable :: default_ny
     character(len=:), allocatable :: default_nz
     character(len=:), allocatable :: default_npx
     character(len=:), allocatable :: default_npy
     character(len=:), allocatable :: default_npz
     character(len=:), allocatable :: default_is_periodic_in_x
     character(len=:), allocatable :: default_is_periodic_in_y
     character(len=:), allocatable :: default_is_periodic_in_z

     ! IO parameters
     character(len=256)            :: dir_path
     character(len=256)            :: prefix
     integer(ip)                   :: nparts
     integer(ip)                   :: reference_fe_geo_order
     integer(ip)                   :: reference_fe_order
     logical                       :: write_solution
     
     character(len=256)            :: triangulation_type
     integer(ip) :: num_dimensions     
     integer(ip) :: number_of_cells_per_dir(0:SPACE_DIM-1)
     integer(ip) :: number_of_parts_per_dir(0:SPACE_DIM-1)
     integer(ip) :: is_dir_periodic(0:SPACE_DIM-1) 

     type(Command_Line_Interface)  :: cli
  contains
     procedure, non_overridable             :: create       => par_pb_bddc_poisson_params_create
     procedure, non_overridable, private    :: set_default  => par_pb_bddc_poisson_params_set_default
     procedure, non_overridable, private    :: add_to_cli   => par_pb_bddc_poisson_params_add_to_cli
     procedure, non_overridable             :: parse        => par_pb_bddc_poisson_params_parse 
     procedure, non_overridable             :: free         => par_pb_bddc_poisson_params_free
     procedure, non_overridable             :: get_dir_path
     procedure, non_overridable             :: get_prefix
     procedure, non_overridable             :: get_nparts
     procedure, non_overridable             :: get_reference_fe_geo_order
     procedure, non_overridable             :: get_reference_fe_order
     procedure, non_overridable             :: get_write_solution
     procedure, non_overridable             :: get_triangulation_type
     procedure, non_overridable             :: get_num_dimensions
  end type par_pb_bddc_poisson_params_t

  ! Types
  public :: par_pb_bddc_poisson_params_t

contains
  subroutine par_pb_bddc_poisson_params_create ( this )
    implicit none
    class(par_pb_bddc_poisson_params_t), intent(inout) :: this
     
    call this%free()
    
    ! Initialize Command Line Interface
    call this%cli%init(progname    = 'par_pb_bddc_poisson', &
                       version     = '', &
                       authors     = '', &
                       license     = '', &
                       description =  'FEMPAR parallel test to solve the 2D Poisson PDE with known analytical solution. &
                                       Boundary set ID 1 MUST BE ASSIGNED to the whole boundary.', & 
                       examples    = ['par_pb_bddc_poisson -h'] )
    call this%set_default()
    call this%add_to_cli()
  end subroutine par_pb_bddc_poisson_params_create 

  !==================================================================================================
  subroutine par_pb_bddc_poisson_params_set_default(this)
    implicit none
    class(par_pb_bddc_poisson_params_t), intent(inout) :: this
    this%default_dir_path     = 'PARTS_4/'
    this%default_prefix       = 'square'
    this%default_nparts       = '4'
    this%default_reference_fe_geo_order = '1'
    this%default_reference_fe_order = '1'
    this%default_write_solution = '.false.'

    this%default_triangulation_type = 'unstructured'
    this%default_num_dimensions = '2'
    this%default_nx = '1'
    this%default_ny = '1'
    this%default_nz = '1'
    this%default_npx = '1'
    this%default_npy = '1'
    this%default_npz = '1'
    this%default_is_periodic_in_x = '0'
    this%default_is_periodic_in_y = '0'
    this%default_is_periodic_in_z = '0'

  end subroutine par_pb_bddc_poisson_params_set_default

  !==================================================================================================
  subroutine par_pb_bddc_poisson_params_add_to_cli(this)
    implicit none
    class(par_pb_bddc_poisson_params_t), intent(inout) :: this
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
    call this%cli%add(switch='--write-solution',switch_ab='-wsolution',help='Write solution in VTK format',&
         &            required=.false.,act='store',def=trim(this%default_write_solution),error=error) 
    check(error==0) 
    
        call this%cli%add(switch='--trinagulation-type',switch_ab='-tt',help='Structured or unstructured (GiD) triangulation?',&
         &            required=.false.,act='store',def=trim(this%default_triangulation_type),choices='structured,unstructured',error=error) 
    check(error==0) 
        call this%cli%add(switch='--number_of_dimensions',switch_ab='-dim',help='Number of space dimensions',&
         &            required=.false.,act='store',def=trim(this%default_num_dimensions),error=error) 
    check(error==0) 
        call this%cli%add(switch='--number_of_cells_in_x',switch_ab='-nx',help='Number of cells in x',&
         &            required=.false.,act='store',def=trim(this%default_nx),error=error) 
    check(error==0) 
        call this%cli%add(switch='--number_of_cells_in_y',switch_ab='-ny',help='Number of cells in y',&
         &            required=.false.,act='store',def=trim(this%default_ny),error=error) 
    check(error==0) 
        call this%cli%add(switch='--number_of_cells_in_z',switch_ab='-nz',help='Number of cells in z',&
         &            required=.false.,act='store',def=trim(this%default_nz),error=error) 
    check(error==0) 

        call this%cli%add(switch='--number_of_parts_in_x',switch_ab='-npx',help='Number of parts in x',&
         &            required=.false.,act='store',def=trim(this%default_npx),error=error) 
    check(error==0) 
        call this%cli%add(switch='--number_of_parts_in_y',switch_ab='-npy',help='Number of parts in y',&
         &            required=.false.,act='store',def=trim(this%default_npy),error=error) 
    check(error==0) 
        call this%cli%add(switch='--number_of_parts_in_z',switch_ab='-npz',help='Number of parts in z',&
         &            required=.false.,act='store',def=trim(this%default_npz),error=error) 
    check(error==0) 
    
    call this%cli%add(switch='--periodic_in_x',switch_ab='-px',help='Is the mesh periodic in x',&
         &            required=.false.,act='store',def=trim(this%default_is_periodic_in_x),error=error) 
    check(error==0) 
        call this%cli%add(switch='--periodic_in_y',switch_ab='-py',help='Is the mesh periodic in y',&
         &            required=.false.,act='store',def=trim(this%default_is_periodic_in_y),error=error) 
    check(error==0) 
        call this%cli%add(switch='--periodic_in_z',switch_ab='-pz',help='Is the mesh periodic in z',&
         &            required=.false.,act='store',def=trim(this%default_is_periodic_in_z),error=error) 
    check(error==0) 

  end subroutine par_pb_bddc_poisson_params_add_to_cli

  !==================================================================================================
  subroutine par_pb_bddc_poisson_params_parse(this,parameter_list)
    implicit none
    class(par_pb_bddc_poisson_params_t), intent(inout) :: this
    type(ParameterList_t)       , intent(inout) :: parameter_list

    integer(ip) :: error, istat

    call this%cli%parse(error=error)
    check(error==0)
    call this%cli%get(switch='-d',val=this%dir_path,error=error); check(error==0)
    call this%cli%get(switch='-p',val=this%prefix,error=error); check(error==0)
    call this%cli%get(switch='-n',val=this%nparts,error=error); check(error==0)
    call this%cli%get(switch='-gorder',val=this%reference_fe_geo_order,error=error); check(error==0)
    call this%cli%get(switch='-order',val=this%reference_fe_order,error=error); check(error==0)
    call this%cli%get(switch='-wsolution',val=this%write_solution,error=error); check(error==0)

    call this%cli%get(switch='-tt',val=this%triangulation_type,error=istat); check(istat==0)
    call this%cli%get(switch='-dim',val=this%num_dimensions,error=istat); check(istat==0)
    call this%cli%get(switch='-nx',val=this%number_of_cells_per_dir(0),error=istat); check(istat==0)
    call this%cli%get(switch='-ny',val=this%number_of_cells_per_dir(1),error=istat); check(istat==0)
    call this%cli%get(switch='-nz',val=this%number_of_cells_per_dir(2),error=istat); check(istat==0)
    call this%cli%get(switch='-npx',val=this%number_of_parts_per_dir(0),error=istat); check(istat==0)
    call this%cli%get(switch='-npy',val=this%number_of_parts_per_dir(1),error=istat); check(istat==0)
    call this%cli%get(switch='-npz',val=this%number_of_parts_per_dir(2),error=istat); check(istat==0)
    call this%cli%get(switch='-px',val=this%is_dir_periodic(0),error=istat); check(istat==0)
    call this%cli%get(switch='-py',val=this%is_dir_periodic(1),error=istat); check(istat==0)
    call this%cli%get(switch='-pz',val=this%is_dir_periodic(2),error=istat); check(istat==0)

    call parameter_list%init()
    istat = 0
    istat = istat + parameter_list%set(key = dir_path_key, value = this%dir_path)
    istat = istat + parameter_list%set(key = prefix_key  , value = this%prefix)
    istat = istat + parameter_list%set(key = geometry_interpolation_order_key  , value = this%reference_fe_geo_order)
    check(istat==0)

    if(trim(this%triangulation_type)=='unstructured') then
       istat = parameter_list%set(key = triangulation_generate_key, value = triangulation_generate_from_mesh)
    else if(trim(this%triangulation_type)=='structured') then
       istat = parameter_list%set(key = triangulation_generate_key         , value = triangulation_generate_structured)
       istat = istat + parameter_list%set(key = number_of_dimensions_key   , value = this%num_dimensions)
       istat = istat + parameter_list%set(key = number_of_cells_per_dir_key, value = this%number_of_cells_per_dir)
       istat = istat + parameter_list%set(key = number_of_parts_per_dir_key, value = this%number_of_parts_per_dir)
       istat = istat + parameter_list%set(key = is_dir_periodic_key        , value = this%is_dir_periodic)
    end if
    check(istat==0)
    
  end subroutine par_pb_bddc_poisson_params_parse
  
  subroutine par_pb_bddc_poisson_params_free(this)
    implicit none
    class(par_pb_bddc_poisson_params_t), intent(inout) :: this
    if(allocated(this%default_dir_path)) deallocate(this%default_dir_path)              
    if(allocated(this%default_prefix)) deallocate(this%default_prefix)                    
    if(allocated(this%default_nparts)) deallocate(this%default_nparts)
    if(allocated(this%default_reference_fe_geo_order)) deallocate(this%default_reference_fe_geo_order)
    if(allocated(this%default_reference_fe_order)) deallocate(this%default_reference_fe_order)
    if(allocated(this%default_write_solution)) deallocate(this%default_write_solution)
    call this%cli%free()
  end subroutine par_pb_bddc_poisson_params_free
  
  ! GETTERS *****************************************************************************************
  function get_dir_path(this)
    implicit none
    class(par_pb_bddc_poisson_params_t) , intent(in) :: this
    character(len=256) :: get_dir_path
    get_dir_path = this%dir_path
  end function get_dir_path

  !==================================================================================================
  function get_prefix(this)
    implicit none
    class(par_pb_bddc_poisson_params_t) , intent(in) :: this
    character(len=256) :: get_prefix
    get_prefix = this%prefix
  end function get_prefix

  !==================================================================================================
  function get_nparts(this)
    implicit none
    class(par_pb_bddc_poisson_params_t) , intent(in) :: this
    integer(ip) :: get_nparts
    get_nparts = this%nparts
  end function get_nparts
  
    !==================================================================================================
  function get_reference_fe_geo_order(this)
    implicit none
    class(par_pb_bddc_poisson_params_t) , intent(in) :: this
    integer(ip) :: get_reference_fe_geo_order
    get_reference_fe_geo_order = this%reference_fe_geo_order
  end function get_reference_fe_geo_order
  
  !==================================================================================================
  function get_reference_fe_order(this)
    implicit none
    class(par_pb_bddc_poisson_params_t) , intent(in) :: this
    integer(ip) :: get_reference_fe_order
    get_reference_fe_order = this%reference_fe_order
  end function get_reference_fe_order
  
  !==================================================================================================
  function get_write_solution(this)
    implicit none
    class(par_pb_bddc_poisson_params_t) , intent(in) :: this
    logical :: get_write_solution
    get_write_solution = this%write_solution
  end function get_write_solution

  !==================================================================================================
  function get_triangulation_type(this)
    implicit none
    class(par_pb_bddc_poisson_params_t) , intent(in) :: this
    character(len=256) :: get_triangulation_type
    get_triangulation_type = this%triangulation_type
  end function get_triangulation_type 

  !==================================================================================================
  function get_num_dimensions(this)
    implicit none
    class(par_pb_bddc_poisson_params_t) , intent(in) :: this
    integer(ip) :: get_num_dimensions
    get_num_dimensions = this%num_dimensions
  end function get_num_dimensions
  
end module par_pb_bddc_poisson_params_names
