module geometry_names
  use iso_c_binding
  use types_names
  use stdio_names
  use memor_names
  use hash_table_names
  use sisl_names
  use field_names
  use FPL
  use parameters_consistency_names
  implicit none
  private
#include "debug.i90"

  type geometric_point_t
     private
     integer(ip) :: id
     real(rp)    :: coord(SPACE_DIM)
   contains
     procedure, non_overridable :: read => point_read
  end type geometric_point_t 

  type line_t
     private
     integer(ip) :: id
     integer(ip) :: point(2)
     ! Nurbs description not available for stlines
     ! which is detected checking when n==0
     integer(ip) :: n=0
     integer(ip) :: p=0                           ! p as in sisl=order+1
     real(rp), allocatable :: control_points(:)   ! (number_space_dimensions+1)*n, including weights
     real(rp), allocatable :: knots(:)            ! nu+pu+1
     type(c_ptr)               :: sisl_ptr = c_null_ptr
     type(geometry_t), pointer :: geometry
   contains
     procedure, non_overridable :: read => line_read
     procedure, non_overridable :: set_geometry_pointer_to => line_set_geometry_pointer_to
     procedure, non_overridable :: init => line_init
     procedure, non_overridable :: get_parameter => line_get_parameter
     procedure, non_overridable :: evaluate      => line_evaluate
  end type line_t

  type surface_t
     private
     integer(ip) :: id
     ! Composition (orientation?)
     integer(ip) :: num_lines
     integer(ip), allocatable :: lines_ids(:)
     ! Nurbs description
     integer(ip) :: nu=0
     integer(ip) :: nv=0
     integer(ip) :: pu=0                          ! p as in sisl=order+1
     integer(ip) :: pv=0
     real(rp), allocatable :: control_points(:)   ! (ndime+1)*nu*nv, including weights
     real(rp), allocatable :: u_knots(:)          ! nu+pu+1
     real(rp), allocatable :: v_knots(:)          ! nv+pv+1
     type(c_ptr) :: sisl_ptr = c_null_ptr
   contains
     procedure, non_overridable :: read => surface_read
  end type surface_t

  type volume_t
     private
     integer(ip) :: id
     ! Composition (orientation?)
     integer(ip) :: num_surfaces
     integer(ip), allocatable :: surfaces_ids(:)
   contains
     procedure, non_overridable :: read => volume_read
  end type volume_t

  type geometry_t
     private
     integer(ip)     :: num_points=0
     integer(ip)     :: num_lines=0
     integer(ip)     :: num_surfaces=0
     integer(ip)     :: num_volumes=0
     type(hash_table_ip_ip_t) :: point_index
     type(hash_table_ip_ip_t) :: line_index
     type(hash_table_ip_ip_t) :: surface_index
     type(hash_table_ip_ip_t) :: volume_index
     type(geometric_point_t) , allocatable :: points(:)
     type(line_t)      , allocatable :: lines(:)
     type(surface_t)   , allocatable :: surfaces(:)
     type(volume_t)    , allocatable :: volumes(:)
   contains
     procedure, non_overridable :: read_from_file  =>  geometry_read_from_file
     procedure, non_overridable :: read_from_unit  =>  geometry_read_from_unit
     generic :: read => read_from_file, read_from_unit
     procedure, non_overridable :: get_point   => geometry_get_point
     procedure, non_overridable :: get_line    => geometry_get_line
     procedure, non_overridable :: get_surface => geometry_get_surface
     procedure, non_overridable :: get_volume  => geometry_get_volume
  end type geometry_t

  public :: geometry_t, geometric_point_t, line_t, surface_t, volume_t

contains

  !=============================================================================
  !
  ! Point TBP
  !
  !=============================================================================
  subroutine point_read(point,unit)
    implicit none
    class(geometric_point_t), intent(inout) :: point
    integer(ip)             , intent(in)    :: unit
    integer(ip)     :: i 
    character(256)  :: tel
    read(unit,'(a)') tel
    do while(tel(1:5)/='Coord')
       read(unit,'(a)') tel
    end do
    read(tel(7:),*) (point%coord(i),i=1,SPACE_DIM)
   
    end subroutine point_read

  !=============================================================================
  !
  ! Line TBP
  !
  !=============================================================================
  subroutine line_read(line,unit)
    implicit none
    class(line_t), intent(inout) :: line
    integer(ip)  , intent(in)  :: unit
    character(256)  :: tel
    character(7)    :: dum
    integer(ip)     :: i,j

    ! Line end points
    read(unit,'(a)') tel
    do while(tel(1:6)/='Points')
       read(unit,'(a)') tel
    end do
    read(tel(8:),*) line%point(1),line%point(2)

    ! Control points
    read(unit,'(a)') tel
    do while(tel(1:25)/='Number of Control Points=')
       if(tel(1:3)=='END') return
       read(unit,'(a)') tel
    end do
    read(tel(26:),*) line%n,dum,line%p
    allocate(line%control_points(line%n*(SPACE_DIM+1)),stat=i)
    !if(i/=0) check(.false.)

    do i=1,line%n
       read(unit,'(a)') tel
       read(tel(16:),*) (line%control_points((SPACE_DIM+1)*(i-1)+j),j=1,SPACE_DIM)
       !read(tel(16:),*) line%control_points((number_space_dimensions+1)*(i-1)+1), &
       !     &           line%control_points((number_space_dimensions+1)*(i-1)+2), &
       !     &           line%control_points((number_space_dimensions+1)*(i-1)+3)
    end do

    ! knots
    read(unit,'(a)') tel
    read(tel(17:),*) i ! Number of knots
    assert(i==line%n+line%p+1)
    allocate(line%knots(line%n+line%p+1),stat=i)
    do i=1,line%n+line%p+1
       read(unit,'(a)') tel
       if(i<10) then
          read(tel(14:),*) line%knots(i)
       else
          read(tel(15:),*) line%knots(i)
       end if
    end do

    ! Weights
    read(unit,'(a)') tel
    do while(tel(1:17)/='Rational weights:')
       if(tel(1:3)=='END') return
       read(unit,'(a)') tel
    end do
    do i=1,line%n
       read(unit,'(a)') tel
       read(tel,*) line%control_points(4*i)
    end do

    ! Define 4D control poitns (i.e. multiply by weights)
    do i=1,line%n
       do j=1,SPACE_DIM
          line%control_points((SPACE_DIM+1)*(i-1)+j) = line%control_points((SPACE_DIM+1)*(i-1)+j) * line%control_points((SPACE_DIM+1)*i)
       end do
    end do

  end subroutine line_read

  !=============================================================================
  subroutine line_set_geometry_pointer_to(line,geometry)
    implicit none
    class(line_t), intent(inout) :: line
    type(geometry_t), target, intent(in) :: geometry
    line%geometry => geometry
  end subroutine line_set_geometry_pointer_to

  !=============================================================================
  subroutine line_init(line)
    implicit none
    class(line_t), intent(inout) :: line
    type(geometric_point_t), pointer :: point 
    integer(ip) :: i

    ! The magic constants: 2 means nurbs, always in 3D and 0 means point (not copy).
    if(line%n>0) then ! It is a nurbs
       line%sisl_ptr = new_curve(line%n,line%p+1,line%knots,line%control_points,2,SPACE_DIM,0)
    else              ! It is a STLINE, create a linear spline using extremes
       line%n = 2
       line%p = 1
       allocate(line%control_points(line%n*SPACE_DIM),stat=i)
       allocate(line%knots(line%n+line%p+1),stat=i)
       line%knots = (/0.0,0.0,1.0,1.0/)
       point => line%geometry%get_point(line%point(1))
       line%control_points(1:SPACE_DIM) =  point%coord
       point => line%geometry%get_point(line%point(2))
       line%control_points(SPACE_DIM+1:2*SPACE_DIM) = point%coord
       line%sisl_ptr = new_curve(line%n,line%p+1,line%knots,line%control_points,1,SPACE_DIM,0)
    end if
  end subroutine line_init

  !=============================================================================
  !function line_get_parameter(line,point_coords,tol)
  function line_get_parameter(line,point,tol)
    implicit none
    class(line_t), intent(in) :: line
    type(point_t), intent(in) :: point
    real(rp)     , intent(in) :: tol
    real(rp)     line_get_parameter

    !real(rp)          :: point_coords(number_space_dimensions)
    integer(ip)       :: p_shape(1)
    real(rp), pointer :: param(:)
    type(c_ptr)       :: p_param
    integer(ip)       :: istat,num_int,num_curves
    type(c_ptr)       :: wcurve

    !point_coords = point%get_value()
    !call point_intersection(line%sisl_ptr, point_coords, number_space_dimensions, tol, num_int, p_param, num_curves, wcurve, istat)
    call point_intersection(line%sisl_ptr, point%get_value(), SPACE_DIM, tol, num_int, p_param, num_curves, wcurve, istat)
    !write(*,*) num_int, istat
    assert(istat==0)
    assert(num_int==1)
    assert(num_curves==0)
    p_shape(1) = num_int
    call c_f_pointer(p_param,param,p_shape)

    line_get_parameter = param(1)
  end function line_get_parameter

  !=============================================================================
  subroutine line_evaluate(line,param,point,tangent)
    implicit none
    class(line_t)     , intent(in)  :: line
    real(rp)          , intent(in)  :: param
    type(point_t)     , intent(out) :: point
    type(vector_field_t), optional, intent(out) :: tangent
    !real(rp)          , intent(out) :: point(3)
    !real(rp), optional, intent(out) :: tangent(3)
    real(rp)    :: values(6)
    integer(ip) :: leftknot = 0
    integer(ip) :: stat
    integer(ip) :: number_of_derivatives_to_compute

    number_of_derivatives_to_compute = 0
    if(present(tangent)) number_of_derivatives_to_compute=1
    
    call curve_left_evaluation(line%sisl_ptr, number_of_derivatives_to_compute, param, leftknot, values, stat) 

    assert(stat==0)
    call point%init(values(1:SPACE_DIM))
    if(present(tangent)) then
       call tangent%init(values(SPACE_DIM+1:2*SPACE_DIM))
    end if
    
    !point = values(1:3)
    !if(present(tangent)) tangent = values(4:6)
  end subroutine line_evaluate

  !=============================================================================
  ! 
  ! Surface TBP
  ! 
  !=============================================================================
  subroutine surface_read(surface,unit)
    implicit none
    class(surface_t), intent(inout) :: surface
    integer(ip)     , intent(in)  :: unit
  end subroutine surface_read

  !=============================================================================
  ! 
  ! Volume TBP
  ! 
  !=============================================================================
  subroutine volume_read(volume,unit)
    implicit none
    class(volume_t), intent(inout) :: volume
    integer(ip)    , intent(in)  :: unit
  end subroutine volume_read

  !=============================================================================
  ! 
  ! Geometry TBP
  ! 
  !=============================================================================
  function geometry_get_point(geometry,id)
    implicit none
    class(geometry_t), target, intent(in) :: geometry
    integer(ip)              , intent(in) :: id
    type(geometric_point_t) , pointer :: geometry_get_point
    integer(ip) :: index,istat
    call geometry%point_index%get(key=id,val=index,stat=istat)
    assert(istat==key_found)
    geometry_get_point => geometry%points(index)
  end function geometry_get_point

  !=============================================================================
  function geometry_get_line(geometry,id)
    implicit none
    class(geometry_t), target, intent(in) :: geometry
    integer(ip)              , intent(in) :: id
    type(line_t) , pointer :: geometry_get_line
    integer(ip) :: index,istat
    call geometry%line_index%get(key=id,val=index,stat=istat)
    assert(istat==key_found)
    geometry_get_line => geometry%lines(index)
  end function geometry_get_line

  !=============================================================================
  function geometry_get_surface(geometry,id)
    implicit none
    class(geometry_t), target, intent(in) :: geometry
    integer(ip)              , intent(in) :: id
    type(surface_t) , pointer :: geometry_get_surface
    integer(ip) :: index,istat
    call geometry%surface_index%get(key=id,val=index,stat=istat)
    assert(istat==key_found)
    geometry_get_surface => geometry%surfaces(index)
  end function geometry_get_surface

  !=============================================================================
  function geometry_get_volume(geometry,id)
    implicit none
    class(geometry_t), target, intent(in) :: geometry
    integer(ip)              , intent(in) :: id
    type(volume_t) , pointer :: geometry_get_volume
    integer(ip) :: index,istat
    call geometry%volume_index%get(key=id,val=index,stat=istat)
    assert(istat==was_stored)
    geometry_get_volume => geometry%volumes(index)
  end function geometry_get_volume
  !=============================================================================
  subroutine geometry_compose_name ( prefix, name ) 
    implicit none
    character(len=*)             , intent(in)    :: prefix 
    character(len=:), allocatable, intent(inout) :: name
    name = trim(prefix) // '.txt'
  end subroutine geometry_compose_name
  !=============================================================================
  subroutine geometry_read_from_file (geometry, parameter_list) ! dir_path, prefix )
     implicit none 
     ! Parameters
     !character (*)    , intent(in)  :: dir_path
     !character (*)    , intent(in)  :: prefix
     type(ParameterList_t), intent(in)    :: parameter_list
     class(geometry_t), intent(out) :: geometry
     ! Locals
     !integer(ip)                    :: lunio
     !character(len=:), allocatable  :: name

     ! Locals
     integer(ip)          :: istat
     character(len=256)   :: dir_path
     character(len=256)   :: prefix
     character(len=:), allocatable   :: name
     integer(ip)                    :: lunio

     ! Mandatory parameters
     assert(parameter_consistency(parameter_list, dir_path_key, dir_path))
     istat = istat + parameter_list%get(key = dir_path_key, value = dir_path)
     check(istat == 0)
     
     assert(parameter_consistency(parameter_list, prefix_key, prefix))
     istat = istat + parameter_list%get(key = prefix_key  , value = prefix)
     check(istat==0)
     
     ! Read geometry
     call geometry_compose_name ( prefix, name )
     lunio = io_open( trim(dir_path)//'/'//trim(name), 'read', status='old' )
     call geometry%read(lunio)
     call io_close(lunio)
   end subroutine geometry_read_from_file
  !=============================================================================
  subroutine geometry_read_from_unit(geometry,unit)
    !------------------------------------------------------------------------
    !
    ! This routine reads a geometry writen by GiD using data report.
    !
    !------------------------------------------------------------------------
    implicit none
    integer(ip)      , intent(in)  :: unit
    class(geometry_t), intent(out) :: geometry
    character(4)    :: dum1
    character(256)  :: tel
    integer(ip)     :: i,istat

    !write(*,*) 'Read geometry:'

    ! Count points, lines, surfaces and volumes
    read(unit,'(a)') tel
    do while(tel(1:12)/='END ENTITIES')
       read(unit,'(a)') tel
       if(tel(1:5)=='POINT') then
          geometry%num_points   = geometry%num_points+1
       else if(tel(1:6)=='STLINE'.or.tel(1:8)=='NURBLINE') then
          geometry%num_lines    = geometry%num_lines+1
       else if(tel(1:11)=='NURBSURFACE') then
          geometry%num_surfaces = geometry%num_surfaces+1
       else if(tel(1:6)=='VOLUME')       then
          geometry%num_volumes  = geometry%num_volumes+1
       end if
    end do
    call io_rewind(unit)

    allocate(geometry%points  (geometry%num_points  ) ,stat=istat)
    allocate(geometry%lines   (geometry%num_lines   ) ,stat=istat)
    allocate(geometry%surfaces(geometry%num_surfaces) ,stat=istat)
    allocate(geometry%volumes (geometry%num_volumes ) ,stat=istat)

    call geometry%point_index%init(geometry%num_points)
    call geometry%line_index%init(geometry%num_lines)
    call geometry%surface_index%init(geometry%num_surfaces)
    call geometry%volume_index%init(geometry%num_volumes)

    geometry%num_points  =0
    geometry%num_lines   =0
    geometry%num_surfaces=0
    geometry%num_volumes =0

    read(unit,'(a)') tel
    do while(tel(1:12)/='END ENTITIES')
       read(unit,'(a)') tel
       if(tel(1:5)=='POINT') then
          geometry%num_points   = geometry%num_points+1
          read(unit,'(a)') tel
          read(tel,*) dum1,geometry%points(geometry%num_points)%id
          !write(*,*) 'Read point', geometry%points(geometry%num_points)%id, geometry%num_points
          call geometry%point_index%put(key=geometry%points(geometry%num_points)%id, &
               &                        val=geometry%num_points,stat=istat)
          call geometry%points(geometry%num_points)%read(unit)
       else if(tel(1:6)=='STLINE'.or.tel(1:8)=='NURBLINE') then
          geometry%num_lines    = geometry%num_lines+1
          read(unit,'(a)') tel
          read(tel,*) dum1,geometry%lines(geometry%num_lines)%id
          !write(*,*) 'Read line', geometry%lines(geometry%num_lines)%id, geometry%num_lines
          call geometry%line_index%put(key=geometry%lines(geometry%num_lines)%id, &
               &                       val=geometry%num_lines,stat=istat)
          call geometry%lines(geometry%num_lines)%set_geometry_pointer_to(geometry)
          call geometry%lines(geometry%num_lines)%read(unit)
       else if(tel(1:11)=='NURBSURFACE') then
          geometry%num_surfaces = geometry%num_surfaces+1
          read(unit,'(a)') tel
          read(tel,*) dum1,geometry%surfaces(geometry%num_surfaces)%id
          call geometry%surface_index%put(key=geometry%surfaces(geometry%num_surfaces)%id, &
               &                          val=geometry%num_surfaces,stat=istat)
          call geometry%surfaces(geometry%num_surfaces)%read(unit)
       else if(tel(1:6)=='VOLUME')       then
          geometry%num_volumes  = geometry%num_volumes+1
          read(unit,'(a)') tel
          read(tel,*) dum1,geometry%volumes(geometry%num_volumes)%id
          call geometry%volume_index%put(key=geometry%volumes(geometry%num_volumes)%id, &
               &                         val=geometry%num_volumes,stat=istat)
          call geometry%volumes(geometry%num_volumes)%read(unit)
       end if
    end do

    do i=1,geometry%num_lines
       call geometry%lines(i)%init()
    end do
    !do i=1,geometry%num_surfaces
    !end do

    !if(geometry%num_points>0)   write(*,*) 'Finally read points:'  , geometry%num_points
    !if(geometry%num_lines>0)    write(*,*) 'Finally read lines:'   , geometry%num_lines   
    !if(geometry%num_surfaces>0) write(*,*) 'Finally read surfaces:', geometry%num_surfaces
    !if(geometry%num_volumes>0)  write(*,*) 'Finally read volumes:' , geometry%num_volumes 

  end subroutine geometry_read_from_unit

end module geometry_names
