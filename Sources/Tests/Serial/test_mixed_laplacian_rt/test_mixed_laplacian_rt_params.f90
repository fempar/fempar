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
module mixed_laplacian_rt_params_names
  use fempar_names
# include "debug.i90"

  implicit none
  private

  type mixed_laplacian_rt_params_t 
     private 
     ! IO parameters
     character(len=:), allocatable :: default_dir_path
     character(len=:), allocatable :: default_prefix
     character(len=:), allocatable :: default_dir_path_out
     character(len=:), allocatable :: default_reference_fe_geo_order
     character(len=:), allocatable :: default_reference_fe_order
     character(len=:), allocatable :: default_triangulation_type
     character(len=:), allocatable :: default_num_dimensions
     character(len=:), allocatable :: default_nx
     character(len=:), allocatable :: default_ny
     character(len=:), allocatable :: default_nz
     character(len=:), allocatable :: default_is_periodic_in_x
     character(len=:), allocatable :: default_is_periodic_in_y
     character(len=:), allocatable :: default_is_periodic_in_z

     
     type(Command_Line_Interface):: cli 

     ! IO parameters
     character(len=256)            :: dir_path
     character(len=256)            :: prefix
     character(len=256)            :: dir_path_out
     integer(ip)                   :: reference_fe_geo_order
     integer(ip)                   :: reference_fe_order
     character(len=256)            :: triangulation_type
     integer(ip)                   :: num_dimensions     
     integer(ip)                   :: number_of_cells_per_dir(0:SPACE_DIM-1)
     integer(ip)                   :: is_dir_periodic(0:SPACE_DIM-1)
     
   contains
     procedure, non_overridable             :: create       => mixed_laplacian_rt_create
     procedure, non_overridable, private    :: set_default  => mixed_laplacian_rt_set_default
     procedure, non_overridable, private    :: add_to_cli   => mixed_laplacian_rt_add_to_cli
     procedure, non_overridable             :: parse        => mixed_laplacian_rt_parse 
     procedure, non_overridable             :: free         => mixed_laplacian_rt_free
     procedure, non_overridable             :: get_dir_path
     procedure, non_overridable             :: get_prefix
     procedure, non_overridable             :: get_dir_path_out
     procedure, non_overridable             :: get_reference_fe_geo_order
     procedure, non_overridable             :: get_reference_fe_order
  end type mixed_laplacian_rt_params_t

  ! Types
  public :: mixed_laplacian_rt_params_t

contains

  subroutine mixed_laplacian_rt_create(this)
    implicit none
    class(mixed_laplacian_rt_params_t), intent(inout) :: this
    
    call this%free()
    
     ! Initialize Command Line Interface
    call this%cli%init(progname    = 'mixed_laplacian_rt',                                                     &
         &        version     = '',                                                                 &
         &        authors     = '',                                                                 &
         &        license     = '',                                                                 &
         &        description =  'FEMPAR test to solve the 2D Mixed Laplacian PDE with known analytical solution. &
                                  Boundary set ID 1 MUST BE ASSIGNED to the whole boundary.', &
         &        examples    = ['mixed_laplacian_rt -h  ', 'mixed_laplacian_rt -h  ' ])
    
    call this%set_default()
    call this%add_to_cli()
  end subroutine mixed_laplacian_rt_create
  
  subroutine mixed_laplacian_rt_set_default(this)
    implicit none
    class(mixed_laplacian_rt_params_t), intent(inout) :: this
    ! IO parameters
    this%default_dir_path       = 'data/'
    this%default_prefix         = 'square'
    this%default_dir_path_out   = 'output/'
    this%default_reference_fe_geo_order = '1'
    this%default_reference_fe_order = '1'
    this%default_triangulation_type = 'unstructured'
    this%default_num_dimensions = '2'
    this%default_nx = '1'
    this%default_ny = '1'
    this%default_nz = '1'
    this%default_is_periodic_in_x = '0'
    this%default_is_periodic_in_y = '0'
    this%default_is_periodic_in_z = '0'
  end subroutine mixed_laplacian_rt_set_default
  
  !==================================================================================================
  subroutine mixed_laplacian_rt_add_to_cli(this)
    implicit none
    class(mixed_laplacian_rt_params_t) , intent(inout) :: this

    ! Locals
    integer(ip) :: error

    ! IO parameters
    call this%cli%add(switch='--dir-path',switch_ab='-d',                              &
         &            help='Directory of the source files',required=.false., act='store',                &
         &            def=trim(this%default_dir_path),error=error)
    check(error==0)
    call this%cli%add(switch='--prefix',switch_ab='-p',help='Name of the GiD files',  &
         &            required=.false.,act='store',def=trim(this%default_prefix),error=error) 
    check(error==0)
    call this%cli%add(switch='--dir-path-out',switch_ab='-o',help='Output Directory',&
         &            required=.false.,act='store',def=trim(this%default_dir_path_out),error=error)
    check(error==0)
    call this%cli%add(switch='--reference-fe-geo-order',switch_ab='-gorder',help='Order of the triangulation reference fe',&
         &            required=.false.,act='store',def=trim(this%default_reference_fe_geo_order),error=error)
    check(error==0)  
    call this%cli%add(switch='--reference-fe-order',switch_ab='-order',help='Order of the fe space reference fe',&
         &            required=.false.,act='store',def=trim(this%default_reference_fe_order),error=error) 
    check(error==0) 
    call this%cli%add(switch='--triangulation-type',switch_ab='-tt',help='Structured or unstructured (GiD) triangulation?',&
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
    call this%cli%add(switch='--periodic_in_x',switch_ab='-px',help='Is the mesh periodic in x',&
         &            required=.false.,act='store',def=trim(this%default_is_periodic_in_x),error=error) 
    check(error==0) 
    call this%cli%add(switch='--periodic_in_y',switch_ab='-py',help='Is the mesh periodic in y',&
         &            required=.false.,act='store',def=trim(this%default_is_periodic_in_y),error=error) 
    check(error==0) 
    call this%cli%add(switch='--periodic_in_z',switch_ab='-pz',help='Is the mesh periodic in z',&
         &            required=.false.,act='store',def=trim(this%default_is_periodic_in_z),error=error) 
    check(error==0) 
    
  end subroutine mixed_laplacian_rt_add_to_cli
  
  subroutine mixed_laplacian_rt_parse(this,parameter_list)
    implicit none
    class(mixed_laplacian_rt_params_t), intent(inout) :: this
    type(ParameterList_t)       , intent(inout) :: parameter_list
    integer(ip) :: istat
    
    call this%cli%parse(error=istat); check(istat==0)
    
    ! IO parameters
    call this%cli%get(switch='-d',val=this%dir_path    ,error=istat); check(istat==0)
    call this%cli%get(switch='-p',val=this%prefix      ,error=istat); check(istat==0)
    call this%cli%get(switch='-o',val=this%dir_path_out,error=istat); check(istat==0)
    call this%cli%get(switch='-gorder',val=this%reference_fe_geo_order,error=istat); check(istat==0)
    call this%cli%get(switch='-order',val=this%reference_fe_order,error=istat); check(istat==0)
    call this%cli%get(switch='-tt',val=this%triangulation_type,error=istat); check(istat==0)
    call this%cli%get(switch='-dim',val=this%num_dimensions,error=istat); check(istat==0)
    call this%cli%get(switch='-nx',val=this%number_of_cells_per_dir(0),error=istat); check(istat==0)
    call this%cli%get(switch='-ny',val=this%number_of_cells_per_dir(1),error=istat); check(istat==0)
    call this%cli%get(switch='-nz',val=this%number_of_cells_per_dir(2),error=istat); check(istat==0)
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
       istat = istat + parameter_list%set(key = is_dir_periodic_key        , value = this%is_dir_periodic)
    end if
    check(istat==0)
    
  end subroutine mixed_laplacian_rt_parse  

  subroutine mixed_laplacian_rt_free(this)
    implicit none
    class(mixed_laplacian_rt_params_t), intent(inout) :: this
    if(allocated(this%default_dir_path)) deallocate(this%default_dir_path)              
    if(allocated(this%default_prefix)) deallocate(this%default_prefix)                    
    if(allocated(this%default_dir_path_out)) deallocate(this%default_dir_path_out)
    if(allocated(this%default_reference_fe_geo_order)) deallocate(this%default_reference_fe_geo_order)
    if(allocated(this%default_reference_fe_order)) deallocate(this%default_reference_fe_order)
    call this%cli%free()
  end subroutine mixed_laplacian_rt_free

  ! GETTERS *****************************************************************************************
  function get_dir_path(this)
    implicit none
    class(mixed_laplacian_rt_params_t) , intent(in) :: this
    character(len=256) :: get_dir_path
    get_dir_path = this%dir_path
  end function get_dir_path

  !==================================================================================================
  function get_prefix(this)
    implicit none
    class(mixed_laplacian_rt_params_t) , intent(in) :: this
    character(len=256) :: get_prefix
    get_prefix = this%prefix
  end function get_prefix

  !==================================================================================================
  function get_dir_path_out(this)
    implicit none
    class(mixed_laplacian_rt_params_t) , intent(in) :: this
    character(len=256) :: get_dir_path_out
    get_dir_path_out = this%dir_path_out
  end function get_dir_path_out
  
  !==================================================================================================
  function get_reference_fe_geo_order(this)
    implicit none
    class(mixed_laplacian_rt_params_t) , intent(in) :: this
    integer(ip) :: get_reference_fe_geo_order
    get_reference_fe_geo_order = this%reference_fe_geo_order
  end function get_reference_fe_geo_order
  
  !==================================================================================================
  function get_reference_fe_order(this)
    implicit none
    class(mixed_laplacian_rt_params_t) , intent(in) :: this
    integer(ip) :: get_reference_fe_order
    get_reference_fe_order = this%reference_fe_order
  end function get_reference_fe_order
  
end module mixed_laplacian_rt_params_names
