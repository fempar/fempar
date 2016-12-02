module geometry_names
  use iso_c_binding
  use types_names
  use stdio_names
  use memor_names
  use hash_table_names
  use sisl_names
  use field_names
  use FPL
  use PENF
  implicit none
  private
#include "debug.i90"

  type geometric_point_t
     private
     integer(ip) :: id
     real(rp)    :: coord(SPACE_DIM)
   contains
     procedure, non_overridable :: create => point_create
     procedure, non_overridable :: read   => point_read
     procedure, non_overridable :: print  => point_print
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
     procedure, non_overridable :: create                  => line_create
     procedure, non_overridable :: read                    => line_read
     procedure, non_overridable :: print                   => line_print
     procedure, non_overridable :: set_geometry_pointer_to => line_set_geometry_pointer_to
     procedure, non_overridable :: init                    => line_init
     procedure, non_overridable :: get_parameter           => line_get_parameter
     procedure, non_overridable :: evaluate                => line_evaluate
  end type line_t

  type surface_t
     private
     integer(ip) :: id
     ! Composition (orientation?)
     integer(ip) :: num_lines
     integer(ip), allocatable :: lines_ids(:)
     integer(ip), allocatable :: lines_orientation(:)
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
     procedure, non_overridable :: read  => surface_read
     procedure, non_overridable :: print => surface_print
  end type surface_t

  type volume_t
     private
     integer(ip) :: id
     ! Composition (orientation?)
     integer(ip) :: num_surfaces
     integer(ip), allocatable :: surfaces_ids(:)
     integer(ip), allocatable :: surfaces_orientation(:)
   contains
     procedure, non_overridable :: read  => volume_read
     procedure, non_overridable :: print => volume_print
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

    subroutine move_forward_to_find_string(unit, string, stopstring, line, position)
    !-----------------------------------------------------------------
    !< Move forward line by line in order to find the given string
    !< and return the position of the string in the line
    !-----------------------------------------------------------------
        integer(ip),    intent(in)    :: unit
        character(*),   intent(in)    :: string
        character(*),   intent(in)    :: stopstring
        character(256), intent(inout) :: line
        integer(ip),    intent(out)   :: position
        integer(ip)                   :: error
    !-----------------------------------------------------------------
        position = index(string=line, substring=string, back=.false., kind=ip)
        do while(position==0 .and. index(string=line, substring=stopstring, back=.false., kind=ip)==0)
            read(unit=unit,fmt='(a)',iostat=error) line
            if(IS_IOSTAT_END(error) .or. IS_IOSTAT_EOR(error)) exit
            position = index(string=line, substring=string, back=.false., kind=ip)
        end do
    end subroutine

  !=============================================================================
  ! Point TBP's
  !=============================================================================

    subroutine point_create(this,id, coords)
    !-----------------------------------------------------------------
    !< Read a point
    !-----------------------------------------------------------------
        class(geometric_point_t), intent(inout) :: this
        integer(ip)             , intent(in)    :: id
        real(rp)                , intent(in)    :: coords(SPACE_DIM)
    !-----------------------------------------------------------------
        this%id = id
        this%coord = coords
    end subroutine point_create


    subroutine point_read(this,unit)
    !-----------------------------------------------------------------
    !< Read a point
    !-----------------------------------------------------------------
        class(geometric_point_t), intent(inout) :: this
        integer(ip)             , intent(in)    :: unit
        integer(ip)                             :: i, pos
        character(256)                          :: string
    !-----------------------------------------------------------------
        read(unit,'(a)') string
        ! Read point ID
        call move_forward_to_find_string(unit,'Num:', 'END', string, pos); assert(pos/=0)
        read(string(pos+4:),*) this%id
        ! Read point coordinates
        call move_forward_to_find_string(unit, 'Coord:', 'END', string, pos); assert(pos/=0)
        read(string(pos+6:),*) (this%coord(i),i=1,SPACE_DIM)
    end subroutine point_read


    subroutine point_print(this,unit)
    !-----------------------------------------------------------------
    !< Print a point
    !-----------------------------------------------------------------
        class(geometric_point_t), intent(inout) :: this
        integer(ip), optional,    intent(in)    :: unit
        integer(ip)                             :: i , the_unit
    !-----------------------------------------------------------------
        the_unit = stdout
        if(present(unit)) the_unit = unit
        write(the_unit,fmt='(A)') 'Point ID: '//trim(adjustl(str(this%Id)))//&
                                  ' Coords: '//trim(adjustl(str(this%coord)))
    end subroutine point_print

    !=============================================================================
    ! Line TBP's
    !=============================================================================

    subroutine line_create(this, id, first_coord, last_coord)
    !-----------------------------------------------------------------
    !< Create a line
    !-----------------------------------------------------------------
        class(line_t), intent(inout) :: this
        integer(ip),   intent(in)    :: id
        real(rp),      intent(in)    :: first_coord(SPACE_DIM)
        real(rp),      intent(in)    :: last_coord(SPACE_DIM)
    !-----------------------------------------------------------------
        this%id = id
        this%point(1) = 2*(id-1)+1
        this%point(2) = 2*(id-1)+2
        this%n = 0
        this%p = 0
    end subroutine line_create


    subroutine line_read(this,unit)
    !-----------------------------------------------------------------
    !< Read a line
    !-----------------------------------------------------------------
        class(line_t), intent(inout) :: this
        integer(ip)  , intent(in)    :: unit
        character(256)               :: string
        integer(ip)                  :: i, j, pos, n, error
    !-----------------------------------------------------------------
        read(unit,'(a)') string
        ! Read line ID
        call move_forward_to_find_string(unit, 'Num:', 'END', string, pos); assert(pos/=0)
        read(string(pos+4:),*) this%id
        ! Read line points
        call move_forward_to_find_string(unit, 'Points:', 'END', string, pos); assert(pos/=0)
        read(string(pos+7:),*) this%point(1),this%point(2)
        ! Read number of control points
        call move_forward_to_find_string(unit, 'Number of Control Points=', 'END', string, pos)
        if(pos /= 0) then ! if has control points is a NURBS
            read(string(pos+25:),*) this%n
            allocate(this%control_points(this%n*(SPACE_DIM+1)), stat=error); assert(error==0)

            ! Read NURBS degree
            call move_forward_to_find_string(unit, 'Degree=', 'END', string, pos); assert(pos/=0)
            read(string(pos+7:),*) this%p

            ! Read NURBS control points coords
            do i=1,this%n
                read(unit,'(a)') string
                call move_forward_to_find_string(unit, 'Point', 'END', string, pos); assert(pos/=0)
                read(string(pos+5:),*) n
                call move_forward_to_find_string(unit, 'coords:', 'END', string, pos); assert(pos/=0)
                read(string(pos+7:),*) (this%control_points((SPACE_DIM+1)*(n-1)+j),j=1,SPACE_DIM)
            end do

            ! Read NURBS number of knots
            call move_forward_to_find_string(unit, 'Number of knots=', 'END', string, pos); assert(pos/=0)
            read(string(pos+16:),*) i
            assert(i==this%n+this%p+1)

            ! Read NURBS knots
            allocate(this%knots(this%n+this%p+1), stat=error); assert(error==0)
            do i=1,this%n+this%p+1
                read(unit,'(a)') string
                call move_forward_to_find_string(unit, 'knot', 'END', string, pos); assert(pos/=0)
                read(string(pos+4:),*) n
                call move_forward_to_find_string(unit, 'value=', 'END', string, pos); assert(pos/=0)
                read(string(pos+6:),*) this%knots(n)
            end do

            ! Read NURBS Weights
            call move_forward_to_find_string(unit, 'Rational weights:', 'END', string, pos); assert(pos/=0)
            do i=1,this%n
                read(unit,'(a)') string
                read(string,*) this%control_points(4*i)
            end do

            ! Define 4D control points (i.e. multiply by weights)
            do i=1,this%n
                do j=1,SPACE_DIM
                    this%control_points((SPACE_DIM+1)*(i-1)+j) = this%control_points((SPACE_DIM+1)*(i-1)+j) * this%control_points((SPACE_DIM+1)*i)
                end do
            end do
        endif

    end subroutine line_read


    subroutine line_print(this,unit)
    !-----------------------------------------------------------------
    !< Print a line
    !-----------------------------------------------------------------
        class(line_t),          intent(inout) :: this
        integer(ip), optional,  intent(in)    :: unit
        integer(ip)                           :: i , the_unit
    !-----------------------------------------------------------------
        the_unit = stdout
        if(present(unit)) the_unit = unit

        write(the_unit,fmt='(A)') 'Line ID: '//trim(adjustl(str(this%Id)))//&
                                   ' Points: '//trim(adjustl(str(this%point)))//&
                                   ' Degree: '//trim(adjustl(str(this%p)))
        if(this%p>0) then
            write(the_unit,fmt='(A)') '    Control points: '
            do i=1, this%n
                write(the_unit,fmt='(A)') '        '//trim(adjustl(str(n=i,no_sign=.true.)))//' '//trim(adjustl(str(this%control_points((i-1)*SPACE_DIM+1:i*SPACE_DIM))))
            enddo
            write(the_unit,fmt='(A)') '    Knots: '
            do i=1, this%n+this%p+1
                write(the_unit,fmt='(A)') '        '//trim(adjustl(str(n=i,no_sign=.true.)))//' '//trim(adjustl(str(this%knots(i))))
            enddo
        endif
  end subroutine line_print

  !=============================================================================
  subroutine line_set_geometry_pointer_to(this,geometry)
    implicit none
    class(line_t),            intent(inout) :: this
    type(geometry_t), target, intent(in)    :: geometry
    this%geometry => geometry
  end subroutine line_set_geometry_pointer_to

  !=============================================================================
  subroutine line_init(this)
    implicit none
    class(line_t), intent(inout) :: this
    type(geometric_point_t), pointer :: point 
    integer(ip) :: i

    ! The magic constants: 2 means nurbs, always in 3D and 0 means point (not copy).
    if(this%n>0) then ! It is a nurbs
       this%sisl_ptr = new_curve(this%n,this%p+1,this%knots,this%control_points,2,SPACE_DIM,0)
    else              ! It is a STLINE, create a linear spline using extremes
       this%n = 2
       this%p = 1
       allocate(this%control_points(this%n*SPACE_DIM),stat=i)
       allocate(this%knots(this%n+this%p+1),stat=i)
       this%knots = (/0.0,0.0,1.0,1.0/)
       point => this%geometry%get_point(this%point(1))
       this%control_points(1:SPACE_DIM) =  point%coord
       point => this%geometry%get_point(this%point(2))
       this%control_points(SPACE_DIM+1:2*SPACE_DIM) = point%coord
       this%sisl_ptr = new_curve(this%n,this%p+1,this%knots,this%control_points,1,SPACE_DIM,0)
    end if
  end subroutine line_init

  !=============================================================================
  function line_get_parameter(this,point,tol)
    implicit none
    class(line_t), intent(in) :: this
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
    call point_intersection(this%sisl_ptr, point%get_value(), SPACE_DIM, tol, num_int, p_param, num_curves, wcurve, istat)
    !write(*,*) num_int, istat
    assert(istat==0)
    assert(num_int==1)
    assert(num_curves==0)
    p_shape(1) = num_int
    call c_f_pointer(p_param,param,p_shape)

    line_get_parameter = param(1)
  end function line_get_parameter

  !=============================================================================
  subroutine line_evaluate(this,param,point,tangent)
    implicit none
    class(line_t)     , intent(in)  :: this
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
    
    call curve_left_evaluation(this%sisl_ptr, number_of_derivatives_to_compute, param, leftknot, values, stat) 

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

    subroutine surface_read(this, unit)
    !-----------------------------------------------------------------
    !< Read a surface
    !-----------------------------------------------------------------
        class(surface_t), intent(inout) :: this
        integer(ip)  ,    intent(in)    :: unit
        character(256)                  :: string
        character(256)                  :: tmpstring
        integer(ip)                     :: i, j, k, pos, number_knots, error, nu, nv, counter
    !-----------------------------------------------------------------
        read(unit,'(a)') string
        ! Read surface ID
        call move_forward_to_find_string(unit, 'Num:', 'END', string, pos); assert(pos/=0)
        read(string(pos+4:),*) this%id
        ! Read surface number of lines
        call move_forward_to_find_string(unit, 'NumLines:', 'END', string, pos); assert(pos/=0)
        read(string(pos+9:),*) this%num_lines
        allocate(this%lines_ids(this%num_lines),stat=error); assert(error==0)
        allocate(this%lines_orientation(this%num_lines),stat=error); assert(error==0)

        ! Read surface lines ID's
        do i=1, this%num_lines
            read(unit, '(a)') string
            call move_forward_to_find_string(unit, 'Line:', 'END', string, pos); assert(pos/=0)
            read(string(pos+5:),*) this%lines_ids(i)
            call move_forward_to_find_string(unit, 'Orientation:', 'END', string, pos); assert(pos/=0)
            read(string(pos+12:),*) tmpstring
            if(trim(adjustl(tmpString)) == 'SAME1ST') then
                this%lines_orientation(i) = 1
            elseif(trim(adjustl(tmpString)) == 'DIFF1ST') then
                this%lines_orientation(i) = -1
            else
                assert(.false.)
            endif
        enddo

        ! Read number of control points
        call move_forward_to_find_string(unit, 'Number of Control Points=', 'END', string, pos)
        if(pos /= 0) then ! if it has control points is a NURBS
            read(string(pos+25:),*) this%nu, this%nv
            allocate(this%control_points(this%nu*this%nv*(SPACE_DIM+1)),stat=error); assert(error==0)
            call move_forward_to_find_string(unit, 'Number of Control Points=', 'END', string, pos)
            ! Read NURBS degree
            call move_forward_to_find_string(unit, 'Degree=', 'END', string, pos); assert(pos/=0)
            read(string(pos+7:),*) this%pu, this%pv

            counter = 0
            do i=1, this%nu
                do j=1, this%nv
                    read(unit, '(a)') string
                    call move_forward_to_find_string(unit, 'Point', 'END', string, pos)
                    read(string(pos+5:),*) nu, nv
                    call move_forward_to_find_string(unit, 'coords:', 'END', string, pos); assert(pos/=0)
                    read(string(pos+7:),*) (this%control_points(counter*SPACE_DIM+k),k=1,SPACE_DIM)
                    counter = counter + 1
                enddo
            enddo

            call move_forward_to_find_string(unit, 'Number of knots in U=', 'END', string, pos)
            read(string(pos+21:),*) number_knots
            assert(number_knots == this%nu+this%pu+1)

            ! Read NURBS knots
            allocate(this%u_knots(this%nu+this%pu+1), stat=error); assert(error==0)
            do i=1, this%nu+this%pu+1
                read(unit,'(a)') string
                call move_forward_to_find_string(unit, 'knot', 'END', string, pos); assert(pos/=0)
                read(string(pos+4:),*) k
                call move_forward_to_find_string(unit, 'value=', 'END', string, pos); assert(pos/=0)
                read(string(pos+6:),*) this%u_knots(k)
            end do

            call move_forward_to_find_string(unit, 'Number of knots in V=', 'END', string, pos)
            read(string(pos+21:),*) number_knots
            assert(number_knots == this%nv+this%pv+1)

            ! Read NURBS knots
            allocate(this%v_knots(this%nv+this%pv+1), stat=error); assert(error==0)
            do i=1, this%nv+this%pv+1
                read(unit,'(a)') string
                call move_forward_to_find_string(unit, 'knot', 'END', string, pos); assert(pos/=0)
                read(string(pos+4:),*) k
                call move_forward_to_find_string(unit, 'value=', 'END', string, pos); assert(pos/=0)
                read(string(pos+6:),*) this%v_knots(k)
            end do

        endif

    end subroutine surface_read


    subroutine surface_print(this,unit)
    !-----------------------------------------------------------------
    !< Print a surface
    !-----------------------------------------------------------------
        class(surface_t),       intent(inout) :: this
        integer(ip), optional,  intent(in)    :: unit
        integer(ip)                           :: i, j, counter, the_unit
    !-----------------------------------------------------------------
        the_unit = stdout
        if(present(unit)) the_unit = unit

        write(the_unit,fmt='(A)') 'Surface ID: '//trim(adjustl(str(this%Id)))//&
                                   ' Lines: '//trim(adjustl(str(this%lines_ids)))//&
                                   ' Orientation: '//trim(adjustl(str(this%lines_orientation)))//&
                                   ' Degree: '//trim(adjustl(str([this%pu,this%pv])))
        if(this%pu>0 .or. this%pv>0) then
            write(the_unit,fmt='(A)') '    Control points: '
            counter = 0
            do i=1, this%nu
                do j=1, this%nv
                    write(the_unit,fmt='(A)') '        '//trim(adjustl(str(n=[i,j],no_sign=.true.)))//' '//&
                                    trim(adjustl(str(this%control_points(counter*SPACE_DIM+1:(counter+1)*SPACE_DIM))))
                    counter = counter + 1
                enddo
            enddo
            write(the_unit,fmt='(A)') '    U_Knots: '
            do i=1, this%nu+this%pu+1
                write(the_unit,fmt='(A)') '        '//trim(adjustl(str(n=i,no_sign=.true.)))//' '//trim(adjustl(str(this%u_knots(i))))
            enddo
            write(the_unit,fmt='(A)') '    V_Knots: '
            do i=1, this%nu+this%pu+1
                write(the_unit,fmt='(A)') '        '//trim(adjustl(str(n=i,no_sign=.true.)))//' '//trim(adjustl(str(this%u_knots(i))))
            enddo
        endif
  end subroutine surface_print


  !=============================================================================
  ! Volume TBP
  !=============================================================================

    subroutine volume_read(this, unit)
    !-----------------------------------------------------------------
    !< Read a volume
    !-----------------------------------------------------------------
        class(volume_t), intent(inout) :: this
        integer(ip)  ,    intent(in)    :: unit
        character(256)                  :: string
        character(256)                  :: tmpstring
        integer(ip)                     :: i, j, k, pos, number_knots, error, nu, nv, counter
    !-----------------------------------------------------------------
        read(unit,'(a)') string
        ! Read volume ID
        call move_forward_to_find_string(unit, 'Num:', 'END', string, pos); assert(pos/=0)
        read(string(pos+4:),*) this%id
        ! Read volume number of surfaces
        call move_forward_to_find_string(unit, 'NumSurfaces:', 'END', string, pos); assert(pos/=0)
        read(string(pos+12:),*) this%num_surfaces
        allocate(this%surfaces_ids(this%num_surfaces),stat=error); assert(error==0)
        allocate(this%surfaces_orientation(this%num_surfaces),stat=error); assert(error==0)

        ! Read volume surfaces ID's
        do i=1, this%num_surfaces
            read(unit, '(a)') string
            call move_forward_to_find_string(unit, 'Surface:', 'END', string, pos); assert(pos/=0)
            read(string(pos+8:),*) this%surfaces_ids(i)
            call move_forward_to_find_string(unit, 'Orientation:', 'END', string, pos); assert(pos/=0)
            read(string(pos+12:),*) tmpstring
            if(trim(adjustl(tmpString)) == 'SAME1ST') then
                this%surfaces_orientation(i) = 1
            elseif(trim(adjustl(tmpString)) == 'DIFF1ST') then
                this%surfaces_orientation(i) = -1
            else
                assert(.false.)
            endif
        enddo
    end subroutine volume_read


    subroutine volume_print(this, unit)
    !-----------------------------------------------------------------
    !< Print a volume
    !-----------------------------------------------------------------
        class(volume_t),        intent(inout) :: this
        integer(ip), optional,  intent(in)    :: unit
        integer(ip)                           :: i, j, counter, the_unit
    !-----------------------------------------------------------------
        the_unit = stdout
        if(present(unit)) the_unit = unit

        write(the_unit,fmt='(A)') 'Volume ID: '//trim(adjustl(str(this%Id)))//&
                                  ' Surfaces: '//trim(adjustl(str(this%surfaces_ids)))//&
                                  ' Orientation: '//trim(adjustl(str(this%surfaces_orientation)))
  end subroutine volume_print

  !=============================================================================
  ! Geometry TBP
  !=============================================================================

  function geometry_get_point(this,id)
    implicit none
    class(geometry_t), target, intent(in) :: this
    integer(ip)              , intent(in) :: id
    type(geometric_point_t) , pointer :: geometry_get_point
    integer(ip) :: index,istat
    call this%point_index%get(key=id,val=index,stat=istat)
    assert(istat==key_found)
    geometry_get_point => this%points(index)
  end function geometry_get_point

  !=============================================================================
  function geometry_get_line(this,id)
    implicit none
    class(geometry_t), target, intent(in) :: this
    integer(ip)              , intent(in) :: id
    type(line_t) , pointer :: geometry_get_line
    integer(ip) :: index,istat
    call this%line_index%get(key=id,val=index,stat=istat)
    assert(istat==key_found)
    geometry_get_line => this%lines(index)
  end function geometry_get_line

  !=============================================================================
  function geometry_get_surface(this, id)
    implicit none
    class(geometry_t), target, intent(in) :: this
    integer(ip)              , intent(in) :: id
    type(surface_t) , pointer :: geometry_get_surface
    integer(ip) :: index,istat
    call this%surface_index%get(key=id,val=index,stat=istat)
    assert(istat==key_found)
    geometry_get_surface => this%surfaces(index)
  end function geometry_get_surface

  !=============================================================================
  function geometry_get_volume(this, id)
    implicit none
    class(geometry_t), target, intent(in) :: this
    integer(ip)              , intent(in) :: id
    type(volume_t) , pointer :: geometry_get_volume
    integer(ip) :: index,istat
    call this%volume_index%get(key=id,val=index,stat=istat)
    assert(istat==was_stored)
    geometry_get_volume => this%volumes(index)
  end function geometry_get_volume

  !=============================================================================
  subroutine geometry_compose_name ( prefix, name ) 
    implicit none
    character(len=*)             , intent(in)    :: prefix 
    character(len=:), allocatable, intent(inout) :: name
    name = trim(prefix) // '.txt'
  end subroutine geometry_compose_name

  !=============================================================================
  subroutine geometry_read_from_file (this, parameter_list) ! dir_path, prefix )
     implicit none 
     ! Parameters
     !character (*)    , intent(in)  :: dir_path
     !character (*)    , intent(in)  :: prefix
     type(ParameterList_t), intent(in)  :: parameter_list
     class(geometry_t),     intent(out) :: this
     ! Locals
     !integer(ip)                    :: lunio
     !character(len=:), allocatable  :: name

     ! Locals
     integer(ip)          :: istat
     character(len=:), allocatable   :: dir_path
     character(len=:), allocatable   :: prefix
     character(len=:), allocatable   :: name
     integer(ip)                     :: lunio

     ! Mandatory parameters
     assert(parameter_list%isAssignable(dir_path_key, 'string'))
     istat = parameter_list%getAsString(key = dir_path_key, string = dir_path)
     assert(istat == 0)
     
     assert(parameter_list%isAssignable(prefix_key, 'string'))
     istat = parameter_list%getAsString(key = prefix_key  , string = prefix)
     assert(istat==0)
     
     ! Read geometry
     call geometry_compose_name ( prefix, name )
     lunio = io_open( trim(dir_path)//'/'//trim(name), 'read', status='old' )
     call this%read(lunio)
     call io_close(lunio)
   end subroutine geometry_read_from_file

  !=============================================================================
  subroutine geometry_read_from_unit(this,unit)
    !------------------------------------------------------------------------
    !< This routine reads a geometry writen by GiD using data report.
    !------------------------------------------------------------------------
    implicit none
    integer(ip)      , intent(in)  :: unit
    class(geometry_t), intent(out) :: this
    character(4)    :: dum1
    character(256)  :: tel
    integer(ip)     :: i,istat

    !write(*,*) 'Read geometry:'

    ! Count points, lines, surfaces and volumes
    read(unit,'(a)') tel
    do while(tel(1:12)/='END ENTITIES')
       read(unit,'(a)') tel
       if(tel(1:5)=='POINT') then
          this%num_points   = this%num_points+1
       else if(tel(1:6)=='STLINE'.or.tel(1:8)=='NURBLINE') then
          this%num_lines    = this%num_lines+1
       else if(tel(1:11)=='NURBSURFACE') then
          this%num_surfaces = this%num_surfaces+1
       else if(tel(1:6)=='VOLUME')       then
          this%num_volumes  = this%num_volumes+1
       end if
    end do
    call io_rewind(unit)

    allocate(this%points  (this%num_points  ) ,stat=istat)
    allocate(this%lines   (this%num_lines   ) ,stat=istat)
    allocate(this%surfaces(this%num_surfaces) ,stat=istat)
    allocate(this%volumes (this%num_volumes ) ,stat=istat)

    call this%point_index%init(this%num_points)
    call this%line_index%init(this%num_lines)
    call this%surface_index%init(this%num_surfaces)
    call this%volume_index%init(this%num_volumes)

    this%num_points  =0
    this%num_lines   =0
    this%num_surfaces=0
    this%num_volumes =0

    read(unit,'(a)') tel
    do while(tel(1:12)/='END ENTITIES')
       read(unit,'(a)') tel
       if(tel(1:5)=='POINT') then
          this%num_points   = this%num_points+1
          call this%points(this%num_points)%read(unit)
          call this%point_index%put(key=this%points(this%num_points)%id, &
               &                        val=this%num_points,stat=istat)
       else if(tel(1:6)=='STLINE'.or.tel(1:8)=='NURBLINE') then
          this%num_lines    = this%num_lines+1
          !write(*,*) 'Read line', geometry%lines(geometry%num_lines)%id, geometry%num_lines
          call this%lines(this%num_lines)%set_geometry_pointer_to(this)
          call this%lines(this%num_lines)%read(unit)
          call this%line_index%put(key=this%lines(this%num_lines)%id, &
               &                       val=this%num_lines,stat=istat)
       else if(tel(1:11)=='NURBSURFACE') then
          this%num_surfaces = this%num_surfaces+1
          call this%surfaces(this%num_surfaces)%read(unit)
          call this%surface_index%put(key=this%surfaces(this%num_surfaces)%id, &
               &                          val=this%num_surfaces,stat=istat)
       else if(tel(1:6)=='VOLUME')       then
          this%num_volumes  = this%num_volumes+1
          call this%volumes(this%num_volumes)%read(unit)
          call this%volume_index%put(key=this%volumes(this%num_volumes)%id, &
               &                         val=this%num_volumes,stat=istat)
       end if
    end do

    do i=1,this%num_lines
       call this%lines(i)%init()
    end do
    !do i=1,geometry%num_surfaces
    !end do

    !if(geometry%num_points>0)   write(*,*) 'Finally read points:'  , geometry%num_points
    !if(geometry%num_lines>0)    write(*,*) 'Finally read lines:'   , geometry%num_lines   
    !if(geometry%num_surfaces>0) write(*,*) 'Finally read surfaces:', geometry%num_surfaces
    !if(geometry%num_volumes>0)  write(*,*) 'Finally read volumes:' , geometry%num_volumes 

  end subroutine geometry_read_from_unit

end module geometry_names
