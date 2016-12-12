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

  integer(ip), parameter, private :: hexa_num_points = 8
  integer(ip), parameter, private :: hexa_num_lines  = 12
  integer(ip), parameter, private :: hexa_num_faces  = 6

  type geometric_point_t
     private
     integer(ip) :: id
     real(rp)    :: coord(SPACE_DIM)
   contains
     procedure, non_overridable :: create => point_create
     procedure, non_overridable :: print  => point_print
  end type geometric_point_t 

  type line_t
     private
     integer(ip)               :: id
     integer(ip)               :: point(2)
     ! Nurbs description not available for stlines
     ! which is detected checking when n==0
     integer(ip)               :: n = 0
     integer(ip)               :: p = 0               ! p as in sisl=order+1
     real(rp), allocatable     :: control_points(:)   ! (number_space_dimensions+1)*n, including weights
     real(rp), allocatable     :: knots(:)            ! nu+pu+1
     type(c_ptr)               :: sisl_ptr = c_null_ptr
     type(geometry_t), pointer :: geometry => null()
   contains
     procedure, non_overridable, private :: create_linear           => line_create_linear
     procedure, non_overridable, private :: create_nurbs            => line_create_nurbs
     procedure, non_overridable          :: print                   => line_print
     procedure, non_overridable          :: set_geometry_pointer_to => line_set_geometry_pointer_to
     procedure, non_overridable          :: init                    => line_init
     procedure, non_overridable          :: get_parameter           => line_get_parameter
     procedure, non_overridable          :: evaluate                => line_evaluate
     procedure, non_overridable          :: free                    => line_free
     generic                             :: create                  => create_linear, create_nurbs
  end type line_t

  type surface_t
     private
     integer(ip)               :: id
     ! Composition (orientation?)
     integer(ip)               :: num_lines
     integer(ip), allocatable  :: lines_ids(:)
     integer(ip), allocatable  :: lines_orientation(:)
     ! Nurbs description
     integer(ip)               :: nu=0
     integer(ip)               :: nv=0
     integer(ip)               :: pu=0                ! p as in sisl=order+1
     integer(ip)               :: pv=0
     real(rp), allocatable     :: control_points(:)   ! (ndime+1)*nu*nv, including weights
     real(rp), allocatable     :: u_knots(:)          ! nu+pu+1
     real(rp), allocatable     :: v_knots(:)          ! nv+pv+1
     type(geometry_t), pointer :: geometry => null()
     type(c_ptr)               :: sisl_ptr = c_null_ptr
   contains
     procedure, non_overridable, private :: create_linear           => surface_create_linear
     procedure, non_overridable, private :: create_nurbs            => surface_create_nurbs
     procedure, non_overridable          :: init                    => surface_init
     procedure, non_overridable          :: print                   => surface_print
     procedure, non_overridable          :: set_geometry_pointer_to => surface_set_geometry_pointer_to
     procedure, non_overridable          :: get_parameter           => surface_get_parameter
     procedure, non_overridable          :: evaluate                => surface_evaluate
     procedure, non_overridable          :: free                    => surface_free
     generic                             :: create                  => create_linear, create_nurbs
  end type surface_t

  type volume_t
     private
     integer(ip)                :: id
     ! Composition (orientation?)
     integer(ip)                :: num_surfaces
     integer(ip), allocatable   :: surfaces_ids(:)
     integer(ip), allocatable   :: surfaces_orientation(:)
   contains
     procedure, non_overridable, private :: create_linear      => volume_create_linear
     procedure, non_overridable          :: print              => volume_print
     procedure, non_overridable          :: free               => volume_free
     generic                             :: create             => create_linear
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
     procedure, non_overridable          :: init                           => geometry_init
     procedure, non_overridable          :: create                         => geometry_create
     procedure, non_overridable          :: add_point                      => geometry_add_point
     procedure, non_overridable, private :: add_linear_line_from_coords    => geometry_add_linear_line_from_coords
     procedure, non_overridable, private :: add_linear_line_from_point_ids => geometry_add_linear_line_from_point_ids
     procedure, non_overridable, private :: add_linear_surf_from_line_ids  => geometry_add_linear_surf_from_line_ids
     procedure, non_overridable, private :: add_linear_vol_from_surf_ids   => geometry_add_linear_vol_from_surf_ids
     procedure, non_overridable, private :: add_nurbs_line_from_coords     => geometry_add_nurbs_line_from_coords
     procedure, non_overridable, private :: add_nurbs_line_from_point_ids  => geometry_add_nurbs_line_from_point_ids
     procedure, non_overridable, private :: add_nurbs_surf_from_line_ids   => geometry_add_nurbs_surf_from_line_ids
     procedure, non_overridable          :: add_quad                       => geometry_add_quad
     procedure, non_overridable          :: add_hexa                       => geometry_add_hexa
     procedure, non_overridable          :: get_point                      => geometry_get_point
     procedure, non_overridable          :: get_line                       => geometry_get_line
     procedure, non_overridable          :: get_surface                    => geometry_get_surface
     procedure, non_overridable          :: get_volume                     => geometry_get_volume
     procedure, non_overridable          :: free                           => geometry_free
     generic                             :: add_line                       => add_linear_line_from_coords,    &
                                                                              add_linear_line_from_point_ids, &
                                                                              add_nurbs_line_from_coords,     &
                                                                              add_nurbs_line_from_point_ids
     generic                             :: add_surface                    => add_linear_surf_from_line_ids, &
                                                                              add_nurbs_surf_from_line_ids
     generic                             :: add_volume                     => add_linear_vol_from_surf_ids
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

    subroutine point_create(this,id, coord)
    !-----------------------------------------------------------------
    !< Read a point
    !-----------------------------------------------------------------
        class(geometric_point_t), intent(inout) :: this
        integer(ip)             , intent(in)    :: id
        real(rp)                , intent(in)    :: coord(SPACE_DIM)
    !-----------------------------------------------------------------
        this%id = id
        this%coord = coord
    end subroutine point_create


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

    subroutine line_free(this)
    !-----------------------------------------------------------------
    !< Free a line
    !-----------------------------------------------------------------
        class(line_t), intent(inout) :: this
    !-----------------------------------------------------------------
        if(allocated(this%control_points)) call memfree(this%control_points, __FILE__, __LINE__)
        if(allocated(this%knots))          call memfree(this%knots, __FILE__, __LINE__)
        this%sisl_ptr = c_null_ptr
        nullify(this%geometry)
        this%id = 0
        this%point = 0
        this%n=0
        this%p=0
    end subroutine line_free

    subroutine line_create_linear(this, line_id, point_ids)
    !-----------------------------------------------------------------
    !< Create a straight line
    !-----------------------------------------------------------------
        class(line_t), intent(inout) :: this
        integer(ip),   intent(in)    :: line_id
        integer(ip),   intent(in)    :: point_ids(2)
    !-----------------------------------------------------------------
        call this%free()
        this%id = line_id
        this%point = point_ids
    end subroutine line_create_linear


    subroutine line_create_nurbs(this, line_id, point_ids, n, p, control_points, knots)
    !-----------------------------------------------------------------
    !< Create a NURBS line
    !-----------------------------------------------------------------
        class(line_t), intent(inout) :: this
        integer(ip),   intent(in)    :: line_id
        integer(ip),   intent(in)    :: point_ids(2)
        integer(ip),   intent(in)    :: n
        integer(ip),   intent(in)    :: p
        real(rp),      intent(IN)    :: control_points(n*(SPACE_DIM+1))
        real(rp),      intent(IN)    :: knots(n+p+1)
    !-----------------------------------------------------------------
        call this%create(line_id, point_ids)
        this%n = n
        this%p = p
        call memalloc(this%n*(SPACE_DIM+1), this%control_points, __FILE__, __LINE__)
        call memalloc(this%n+this%p+1, this%knots, __FILE__, __LINE__)
        this%control_points = control_points
        this%knots = knots
    end subroutine line_create_nurbs


    subroutine line_print(this,unit)
    !-----------------------------------------------------------------
    !< Print a line
    !-----------------------------------------------------------------
        class(line_t),          intent(in) :: this
        integer(ip), optional,  intent(in) :: unit
        integer(ip)                        :: i , the_unit
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
            write(the_unit,fmt='(A)') '    Rational weights: '
            do i=this%n*SPACE_DIM+1, this%n*(SPACE_DIM+1)
                write(the_unit,fmt='(A)') '        '//trim(adjustl(str(this%control_points(i))))
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
       call memalloc(this%n*SPACE_DIM, this%control_points, __FILE__, __LINE__)
       call memalloc(this%n+this%p+1, this%knots, __FILE__, __LINE__)
       this%knots = [0.0,0.0,1.0,1.0]
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
    real(rp)                  :: line_get_parameter

    !real(rp)          :: point_coords(number_space_dimensions)
    integer(ip)       :: p_shape(1)
    real(rp), pointer :: param(:)
    type(c_ptr)       :: p_param
    integer(ip)       :: istat,num_int,num_curves
    type(c_ptr)       :: wcurve

    !point_coords = point%get_value()
    !call curve_point_intersection(line%sisl_ptr, point_coords, number_space_dimensions, tol, num_int, p_param, num_curves, wcurve, istat)
    call curve_point_intersection(this%sisl_ptr, point%get_value(), SPACE_DIM, tol, num_int, p_param, num_curves, wcurve, istat)
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
  ! Surface TBP
  !=============================================================================

    subroutine surface_free(this)
    !-----------------------------------------------------------------
    !< Create a line
    !-----------------------------------------------------------------
        class(surface_t), intent(inout) :: this
    !-----------------------------------------------------------------
        if(allocated(this%lines_ids))         call memfree(this%lines_ids, __FILE__, __LINE__)
        if(allocated(this%lines_orientation)) call memfree(this%lines_orientation, __FILE__, __LINE__)
        if(allocated(this%control_points))    call memfree(this%control_points, __FILE__, __LINE__)
        if(allocated(this%u_knots))           call memfree(this%u_knots, __FILE__, __LINE__)
        if(allocated(this%v_knots))           call memfree(this%v_knots, __FILE__, __LINE__)
        this%sisl_ptr  = c_null_ptr
        this%id        = 0
        this%num_lines = 0
        this%nu        = 0
        this%nv        = 0
        this%pu        = 0
        this%pv        = 0
    end subroutine surface_free


    subroutine surface_create_linear(this, surface_id, number_lines, lines_ids, lines_orientation)
    !-----------------------------------------------------------------
    !< Create a linear quad
    !-----------------------------------------------------------------
        class(surface_t), intent(inout) :: this
        integer(ip),      intent(in)    :: surface_id
        integer(ip),      intent(in)    :: number_lines
        integer(ip),      intent(in)    :: lines_ids(number_lines)
        integer(ip),      intent(in)    :: lines_orientation(number_lines)
    !-----------------------------------------------------------------
        call this%free()
        this%id        = surface_id
        this%num_lines = number_lines
        call memalloc(this%num_lines, this%lines_ids, __FILE__, __LINE__)
        call memalloc(this%num_lines, this%lines_orientation, __FILE__, __LINE__)
        this%lines_ids = lines_ids
        this%lines_orientation = lines_orientation
    end subroutine surface_create_linear



    subroutine surface_create_nurbs(this, surface_id, number_lines, lines_ids, lines_orientation, nu, nv, pu, pv, control_points, u_knots, v_knots)
    !-----------------------------------------------------------------
    !< Create a nurbs quad
    !-----------------------------------------------------------------
        class(surface_t), intent(inout) :: this
        integer(ip),      intent(in)    :: surface_id
        integer(ip),      intent(in)    :: number_lines
        integer(ip),      intent(in)    :: lines_ids(number_lines)
        integer(ip),      intent(in)    :: lines_orientation(number_lines)
        integer(ip),      intent(in)    :: nu
        integer(ip),      intent(in)    :: nv
        integer(ip),      intent(in)    :: pu
        integer(ip),      intent(in)    :: pv
        real(rp),         intent(in)    :: control_points(nu*nv*(SPACE_DIM+1))
        real(rp),         intent(in)    :: u_knots(nu+pu+1)
        real(rp),         intent(in)    :: v_knots(nv+pv+1)
    !-----------------------------------------------------------------
        call this%create(surface_id, number_lines, lines_ids, lines_orientation)
        this%nu = nu
        this%nv = nv
        this%pu = pu
        this%pv = pv
        call memalloc(this%nu*this%nv*(SPACE_DIM+1), this%control_points, __FILE__, __LINE__)
        call memalloc(this%nu+this%pu+1, this%u_knots, __FILE__, __LINE__)
        call memalloc(this%nv+this%pv+1, this%v_knots, __FILE__, __LINE__)
        this%control_points = control_points
        this%u_knots = u_knots
        this%v_knots = v_knots
    end subroutine surface_create_nurbs


    subroutine surface_init(this)
    !-----------------------------------------------------------------
    !< Initialize a surface
    !-----------------------------------------------------------------
        class(surface_t), intent(inout)  :: this
        type(line_t),            pointer :: line
        type(geometric_point_t), pointer :: point
        integer(ip), allocatable         :: unique_point_ids(:)
        integer(ip)                      :: i, num_unique_points
    !-----------------------------------------------------------------
        ! The magic constants: 2 means nurbs, always in 3D and 0 means point (not copy).
        if(this%nu>0 .or. this%nv>0) then ! It is a nurbs
            this%sisl_ptr = new_surface(this%nu, this%nv, this%pu+1, this%pv+1, this%u_knots, this%v_knots, this%control_points, 2, SPACE_DIM, 0)
        else              ! It is a STLINE, create a linear spline using extremes
            this%nu = 2
            this%nv = 2
            this%pu = 1
            this%pv = 1
            call memalloc(this%nu*this%nv*SPACE_DIM, this%control_points, __FILE__, __LINE__)
            call memalloc(this%nu+this%pu+1, this%u_knots, __FILE__, __LINE__)
            call memalloc(this%nv+this%pv+1, this%v_knots, __FILE__, __LINE__)
            this%u_knots = [0.0,0.0,1.0,1.0]
            this%v_knots = [0.0,0.0,1.0,1.0]
            call memalloc(this%num_lines, unique_point_ids, __FILE__, __LINE__)
            unique_point_ids = 0
            num_unique_points = 0
            do i=1, this%num_lines
                line => this%geometry%get_line(this%lines_ids(i))
                if(.not. any(unique_point_ids == line%point(1))) then
                    num_unique_points = num_unique_points+1
                    unique_point_ids(num_unique_points) = line%point(1)
                endif
                if(num_unique_points == this%num_lines) exit
                if(.not. any(unique_point_ids == line%point(2))) then
                    num_unique_points = num_unique_points+1
                    unique_point_ids(num_unique_points) = line%point(2)
                endif
                if(num_unique_points == this%num_lines) exit
            enddo
            do i=1, num_unique_points
                point => this%geometry%get_point(unique_point_ids(i))
                this%control_points((i-1)*SPACE_DIM+1:i*SPACE_DIM) = point%coord
            enddo
            this%sisl_ptr = new_surface(this%nu, this%nv, this%pu+1, this%pv+1, this%u_knots, this%v_knots, this%control_points, 1, SPACE_DIM, 0)
            call memfree(unique_point_ids, __FILE__, __LINE__)
        end if
    end subroutine surface_init


    function surface_get_parameter(this, point, tol)
    !-----------------------------------------------------------------
    !< Surface get parameter
    !-----------------------------------------------------------------
        class(surface_t), intent(in) :: this
        type(point_t),    intent(in) :: point
        real(rp),         intent(in) :: tol
        real(rp)                     :: surface_get_parameter
        integer(ip)                  :: p_shape(1)
        real(rp), pointer            :: param(:)
        type(c_ptr)                  :: p_param
        integer(ip)                  :: istat,num_int,num_curves
        type(c_ptr)                  :: wcurve
    !-----------------------------------------------------------------
        call surface_point_intersection(this%sisl_ptr, point%get_value(), SPACE_DIM, tol, num_int, p_param, num_curves, wcurve, istat)
        assert(istat==0)
        assert(num_int==1)
        assert(num_curves==0)
        p_shape(1) = num_int
        call c_f_pointer(p_param,param,p_shape)

        surface_get_parameter = param(1)
    end function surface_get_parameter


    subroutine surface_evaluate(this, param, point, tangent)
    !-----------------------------------------------------------------
    !< Surface evaluation
    !-----------------------------------------------------------------
        class(surface_t),               intent(in)    :: this
        real(rp),                       intent(in)    :: param(2)
        type(point_t),                  intent(inout) :: point
        type(vector_field_t), optional, intent(inout) :: tangent
        real(rp)                                      :: values(6)
        real(rp)                                      :: normal(3)
        integer(ip)                                   :: leftknot1
        integer(ip)                                   :: leftknot2
        integer(ip)                                   :: stat
        integer(ip)                                   :: number_of_derivatives_to_compute
    !-----------------------------------------------------------------
        leftknot1 = 0
        leftknot2 = 0
        number_of_derivatives_to_compute = 0
        if(present(tangent)) number_of_derivatives_to_compute=1

        call surface_left_evaluation(this%sisl_ptr, number_of_derivatives_to_compute, param, leftknot1, leftknot2, values, normal, stat) 

        assert(stat==0)
        call point%init(values(1:SPACE_DIM))
        if(present(tangent)) then
           call tangent%init(values(SPACE_DIM+1:2*SPACE_DIM))
        end if
    end subroutine surface_evaluate


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
                                    trim(adjustl(str(this%control_points(counter*(SPACE_DIM+1)+1:(counter+1)*(SPACE_DIM+1)-1))))
                    counter = counter + 1
                enddo
            enddo
            write(the_unit,fmt='(A)') '    Rational wights: '
            do i=1, this%nu*this%nv
                write(the_unit,fmt='(A)') '        '//trim(adjustl(str(this%control_points(i*SPACE_DIM+1))))
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


    subroutine surface_set_geometry_pointer_to(this,geometry)
    !-----------------------------------------------------------------
    !< Set a pointer to the geometry
    !-----------------------------------------------------------------
        class(surface_t),            intent(inout) :: this
        type(geometry_t), target, intent(in)    :: geometry
    !-----------------------------------------------------------------
        this%geometry => geometry
    end subroutine surface_set_geometry_pointer_to


  !=============================================================================
  ! Volume TBP
  !=============================================================================

    subroutine volume_free(this)
    !-----------------------------------------------------------------
    !< Free a volume
    !-----------------------------------------------------------------
        class(volume_t),  intent(inout) :: this
    !-----------------------------------------------------------------
        if(allocated(this%surfaces_ids)) call memfree(this%surfaces_ids, __FILE__, __LINE__)
        if(allocated(this%surfaces_orientation)) call memfree(this%surfaces_orientation, __FILE__, __LINE__)
        this%id           = 0
        this%num_surfaces = 0
    end subroutine volume_free


    subroutine volume_create_linear(this, volume_id, number_surfaces, surfaces_ids, surfaces_orientation)
    !-----------------------------------------------------------------
    !< Create a valume
    !-----------------------------------------------------------------
        class(volume_t),  intent(inout) :: this
        integer(ip),      intent(in)    :: volume_id
        integer(ip),      intent(in)    :: number_surfaces
        integer(ip),      intent(in)    :: surfaces_ids(number_surfaces)
        integer(ip),      intent(in)    :: surfaces_orientation(number_surfaces)
    !-----------------------------------------------------------------
        call this%free()
        this%id = volume_id
        this%num_surfaces = number_surfaces
        call memalloc(this%num_surfaces, this%surfaces_ids, __FILE__, __LINE__)
        call memalloc(this%num_surfaces, this%surfaces_orientation, __FILE__, __LINE__)
        this%surfaces_ids = surfaces_ids
        this%surfaces_orientation = surfaces_orientation
    end subroutine volume_create_linear


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

    subroutine geometry_free(this)
    !-----------------------------------------------------------------
    !< Free the geometry
    !-----------------------------------------------------------------
        class(geometry_t),        intent(inout) :: this
        integer(ip)                             :: i
    !-----------------------------------------------------------------

        if(allocated(this%points)) deallocate(this%points)

        if(allocated(this%lines)) then
            do i=1, size(this%lines)
                call this%lines(i)%free()
            enddo
            deallocate(this%lines)
        endif

        if(allocated(this%surfaces)) then
            do i=1, size(this%surfaces)
                call this%surfaces(i)%free()
            enddo
            deallocate(this%surfaces)
        endif

        if(allocated(this%volumes)) then
            do i=1, size(this%volumes)
                call this%volumes(i)%free()
            enddo
            deallocate(this%volumes)
        endif

        call this%point_index%free()
        call this%line_index%free()
        call this%surface_index%free()
        call this%volume_index%free()

        this%num_points   = 0
        this%num_lines    = 0
        this%num_surfaces = 0
        this%num_volumes  = 0
    end subroutine geometry_free


    subroutine geometry_create(this, num_points, num_lines, num_surfaces, num_volumes)
    !-----------------------------------------------------------------
    !< Create the geometry and allocate its members
    !-----------------------------------------------------------------
        class(geometry_t),        intent(inout) :: this
        integer(ip),              intent(in)    :: num_points
        integer(ip),              intent(in)    :: num_lines
        integer(ip),              intent(in)    :: num_surfaces
        integer(ip),              intent(in)    :: num_volumes
        integer(ip)                             :: error
    !-----------------------------------------------------------------
        call this%free()

        allocate(this%points  (num_points  ) ,stat=error)
        allocate(this%lines   (num_lines   ) ,stat=error)
        allocate(this%surfaces(num_surfaces) ,stat=error)
        allocate(this%volumes (num_volumes ) ,stat=error)

        call this%point_index%init(num_points)
        call this%line_index%init(num_lines)
        call this%surface_index%init(num_surfaces)
        call this%volume_index%init(num_volumes)
    end subroutine geometry_create

    subroutine geometry_init(this)
    !-----------------------------------------------------------------
    !< Initialize SISL nubs
    !-----------------------------------------------------------------
        class(geometry_t),        intent(inout) :: this
        integer(ip)                             :: i
    !-----------------------------------------------------------------
        do i=1,this%num_lines
            call this%lines(i)%init()
        end do
        do i=1,this%num_surfaces
            call this%surfaces(i)%init()
        end do
    end subroutine geometry_init


    subroutine geometry_add_point(this, point_id, coord)
    !-----------------------------------------------------------------
    !< Create a point
    !-----------------------------------------------------------------
        class(geometry_t),        intent(inout) :: this
        integer(ip),              intent(in)    :: point_id
        real(rp),                 intent(in)    :: coord(SPACE_DIM)
        integer(ip)                             :: i
        integer(ip)                             :: error
    !-----------------------------------------------------------------
        assert(allocated(this%points) .and. size(this%points)>= this%num_points+1)
        this%num_points = this%num_points+1
        call this%points(this%num_points)%create(id=point_id, coord=coord)
        call this%point_index%put(key=this%points(this%num_points)%id, val=this%num_points, stat=error)
    end subroutine geometry_add_point


    subroutine geometry_add_linear_line_from_coords(this, line_id, coord1, coord2)
    !-----------------------------------------------------------------
    !< Create a straight line given its coords
    !-----------------------------------------------------------------
        class(geometry_t),         intent(inout) :: this
        integer(ip),               intent(in)    :: line_id
        real(rp),                  intent(in)    :: coord1(SPACE_DIM)
        real(rp),                  intent(in)    :: coord2(SPACE_DIM)
        integer(ip)                              :: i, num_points
        integer(ip)                              :: pid1, pid2
        integer(ip)                              :: error
    !-----------------------------------------------------------------
        num_points = this%num_points
        pid1 = num_points+1; call this%add_point(pid1, coord1)
        pid2 = num_points+2; call this%add_point(pid2, coord2)
        call this%add_line(line_id, [pid1,pid2])
    end subroutine geometry_add_linear_line_from_coords


    subroutine geometry_add_nurbs_line_from_coords(this, line_id, coord1, coord2, n, p, control_points, knots)
    !-----------------------------------------------------------------
    !< Create a NURBS line given its coords
    !-----------------------------------------------------------------
        class(geometry_t),         intent(inout) :: this
        integer(ip),               intent(in)    :: line_id
        real(rp),                  intent(in)    :: coord1(SPACE_DIM)
        real(rp),                  intent(in)    :: coord2(SPACE_DIM)
        integer(ip),               intent(in)    :: n
        integer(ip),               intent(in)    :: p
        real(rp),                  intent(IN)    :: control_points(n*(SPACE_DIM+1))
        real(rp),                  intent(IN)    :: knots(n+p+1)
        integer(ip)                              :: i
        integer(ip)                              :: pid1, pid2
        integer(ip)                              :: error
    !-----------------------------------------------------------------
        call this%add_line(line_id, coord1, coord2)
        call this%add_line(line_id, [pid1,pid2], n, p, control_points, knots)
    end subroutine geometry_add_nurbs_line_from_coords


    subroutine geometry_add_linear_line_from_point_ids(this, line_id, point_ids)
    !-----------------------------------------------------------------
    !< Create a line
    !-----------------------------------------------------------------
        class(geometry_t),         intent(inout) :: this
        integer(ip),               intent(in)    :: line_id
        integer(ip),               intent(in)    :: point_ids(2)
        integer(ip)                              :: i
        type(line_t), allocatable                :: tmplines(:)
        integer(ip)                              :: error
    !-----------------------------------------------------------------
        assert(allocated(this%lines) .and. size(this%lines)>= this%num_lines+1)
        this%num_lines = this%num_lines+1
        call this%lines(this%num_lines)%create(line_id=line_id, point_ids=point_ids)
        call this%lines(this%num_lines)%set_geometry_pointer_to(this)
        call this%line_index%put(key=this%lines(this%num_lines)%id, val=this%num_lines,stat=error)
    end subroutine geometry_add_linear_line_from_point_ids


    subroutine geometry_add_nurbs_line_from_point_ids(this, line_id, point_ids, n, p, control_points, knots)
    !-----------------------------------------------------------------
    !< Create a NURBS line given its coords
    !-----------------------------------------------------------------
        class(geometry_t),         intent(inout) :: this
        integer(ip),               intent(in)    :: line_id
        integer(ip),               intent(in)    :: point_ids(2)
        integer(ip),               intent(in)    :: n
        integer(ip),               intent(in)    :: p
        real(rp),                  intent(IN)    :: control_points(n*(SPACE_DIM+1))
        real(rp),                  intent(IN)    :: knots(n+p+1)
        integer(ip)                              :: i
        integer(ip)                              :: pid1, pid2
        integer(ip)                              :: error
    !-----------------------------------------------------------------
        assert(allocated(this%lines) .and. size(this%lines)>= this%num_lines+1)
        this%num_lines = this%num_lines+1
        call this%lines(this%num_lines)%create(line_id, point_ids, n, p, control_points, knots)
        call this%lines(this%num_lines)%set_geometry_pointer_to(this)
        call this%line_index%put(key=this%lines(this%num_lines)%id, val=this%num_lines,stat=error)
    end subroutine geometry_add_nurbs_line_from_point_ids


    subroutine geometry_add_quad(this, coord1, coord2, coord3, coord4)
    !-----------------------------------------------------------------
    !< Create a quad
    !< 3-------4     .---2---.
    !< |       |     |       |
    !< |       |     3       4
    !< |       |     |       |
    !< 1-------2     .---1---.
    !-----------------------------------------------------------------
        class(geometry_t),         intent(inout) :: this
        real(rp),                  intent(in)    :: coord1(SPACE_DIM)
        real(rp),                  intent(in)    :: coord2(SPACE_DIM)
        real(rp),                  intent(in)    :: coord3(SPACE_DIM)
        real(rp),                  intent(in)    :: coord4(SPACE_DIM)
        integer(ip)                              :: pid1, pid2, pid3, pid4
        integer(ip)                              :: lid1, lid2, lid3, lid4
        integer(ip)                              :: i, num_points, num_lines
        integer(ip)                              :: error
    !-----------------------------------------------------------------
        ! Create points
        !< 3-------4 
        !< |       | 
        !< |       | 
        !< |       | 
        !< 1-------2 
        num_points = this%num_points
        pid1 = num_points+1; call this%add_point(pid1, coord1)
        pid2 = num_points+2; call this%add_point(pid2, coord2)
        pid3 = num_points+3; call this%add_point(pid3, coord3)
        pid4 = num_points+4; call this%add_point(pid4, coord4)
        ! Create lines
        !< .---2---.
        !< |       |
        !< 3       4
        !< |       |
        !< .---1---.
        num_lines = this%num_lines
        lid1 = num_lines+1; call this%add_line(lid1, point_ids=[pid1, pid2])
        lid2 = num_lines+2; call this%add_line(lid2, point_ids=[pid3, pid4])
        lid3 = num_lines+3; call this%add_line(lid3, point_ids=[pid1, pid3])
        lid4 = num_lines+4; call this%add_line(lid4, point_ids=[pid2, pid4])
        ! Create quad
        call this%add_surface(this%num_surfaces+1, 4, [lid1, lid2, lid3, lid4], lines_orientation=[1,1,1,1])
    end subroutine geometry_add_quad


    subroutine geometry_add_linear_surf_from_line_ids(this, surface_id, number_lines, lines_ids, lines_orientation)
    !-----------------------------------------------------------------
    !< Create a linear surface given the information about its lines
    !-----------------------------------------------------------------
        class(geometry_t),         intent(inout) :: this
        integer(ip),               intent(in)    :: surface_id
        integer(ip),               intent(in)    :: number_lines
        integer(ip),               intent(in)    :: lines_ids(number_lines)
        integer(ip),               intent(in)    :: lines_orientation(number_lines)
        integer(ip)                              :: i
        integer(ip)                              :: error
    !-----------------------------------------------------------------
        assert(allocated(this%surfaces) .and. size(this%surfaces)>= this%num_surfaces+1)
        this%num_surfaces = this%num_surfaces+1
        call this%surfaces(this%num_surfaces)%create(surface_id=surface_id,      &
                                                     number_lines=number_lines,  &
                                                     lines_ids=lines_ids,        &
                                                     lines_orientation=lines_orientation)
        call this%surfaces(this%num_surfaces)%set_geometry_pointer_to(this)
        call this%surface_index%put(key=this%surfaces(this%num_surfaces)%id, val=surface_id,stat=error)
    end subroutine geometry_add_linear_surf_from_line_ids


    subroutine geometry_add_nurbs_surf_from_line_ids(this, surface_id, number_lines, lines_ids, lines_orientation, &
                                                     nu, nv, pu, pv, control_points, u_knots, v_knots)
    !-----------------------------------------------------------------
    !< Create a nurbs surface given the information about its lines
    !-----------------------------------------------------------------
        class(geometry_t),         intent(inout) :: this
        integer(ip),               intent(in)    :: surface_id
        integer(ip),               intent(in)    :: number_lines
        integer(ip),               intent(in)    :: lines_ids(number_lines)
        integer(ip),               intent(in)    :: lines_orientation(number_lines)
        integer(ip),               intent(in)    :: nu
        integer(ip),               intent(in)    :: nv
        integer(ip),               intent(in)    :: pu
        integer(ip),               intent(in)    :: pv
        real(rp),                  intent(IN)    :: control_points(nu*nv*(SPACE_DIM+1))
        real(rp),                  intent(IN)    :: u_knots(nu+pu+1)
        real(rp),                  intent(IN)    :: v_knots(nv+pv+1)
        integer(ip)                              :: i
        integer(ip)                              :: error
    !-----------------------------------------------------------------
        assert(allocated(this%surfaces) .and. size(this%surfaces)>= this%num_surfaces+1)
        this%num_surfaces = this%num_surfaces+1
        call this%surfaces(this%num_surfaces)%create(surface_id=surface_id,               &
                                                     number_lines=number_lines,           &
                                                     lines_ids=lines_ids,                 &
                                                     lines_orientation=lines_orientation, &
                                                     nu=nu, nv=nv, pu=pu, pv=pv,          &
                                                     control_points=control_points,       &
                                                     u_knots=u_knots, v_knots=v_knots)
        call this%surfaces(this%num_surfaces)%set_geometry_pointer_to(this)
        call this%surface_index%put(key=this%surfaces(this%num_surfaces)%id, val=this%num_surfaces,stat=error)
    end subroutine geometry_add_nurbs_surf_from_line_ids


    subroutine geometry_add_hexa(this, coord1, coord2, coord3, coord4, coord5, coord6, coord7, coord8)
    !-----------------------------------------------------------------
    !< Create a hexa
    !<     7-------8       .---4---.         .-------.
    !<    /       /|      /|      /|        /|      /|
    !<   /       / |     7 11    8 12      / | 2   / |<6
    !<  /       /  |    /  |    /  |    3>/  |    /  |
    !< 5---3---6___4   .---3---.___.     .-------. 4 .
    !< |  /    |  /    |  /    |  /      |  /    |  /
    !  | /     | /     9 5    10 6       | / 5   | /
    !< |/      |/      |/      |/        |/      |/
    !< 1-------2       .---1---.         .-------.                    
    !<                                       ^ 
    !<                                       1
    !-----------------------------------------------------------------
        class(geometry_t),         intent(inout) :: this
        real(rp),                  intent(in)    :: coord1(SPACE_DIM)
        real(rp),                  intent(in)    :: coord2(SPACE_DIM)
        real(rp),                  intent(in)    :: coord3(SPACE_DIM)
        real(rp),                  intent(in)    :: coord4(SPACE_DIM)
        real(rp),                  intent(in)    :: coord5(SPACE_DIM)
        real(rp),                  intent(in)    :: coord6(SPACE_DIM)
        real(rp),                  intent(in)    :: coord7(SPACE_DIM)
        real(rp),                  intent(in)    :: coord8(SPACE_DIM)
        integer(ip)                              :: pid1, pid2, pid3, pid4, pid5, pid6, pid7, pid8
        integer(ip)                              :: lid1, lid2, lid3, lid4, lid5, lid6
        integer(ip)                              :: lid7, lid8, lid9, lid10, lid11, lid12
        integer(ip)                              :: sid1, sid2, sid3, sid4, sid5, sid6
        integer(ip)                              :: i, num_points, num_lines, num_surfaces
        integer(ip)                              :: error
    !-----------------------------------------------------------------
        assert(allocated(this%volumes) .and. size(this%volumes)>= this%num_volumes+1)
        ! Create 8xpoints
        !<     7-------8
        !<    /       /|
        !<   /       / |
        !<  /       /  |
        !< 5---3---6___4
        !< |  /    |  / 
        !< | /     | / 
        !< |/      |/ 
        !< 1-------2 
        num_points = this%num_points
        pid1 = num_points+1; call this%add_point(pid1, coord1)
        pid2 = num_points+2; call this%add_point(pid2, coord2)
        pid3 = num_points+3; call this%add_point(pid3, coord3)
        pid4 = num_points+4; call this%add_point(pid4, coord4)
        pid5 = num_points+5; call this%add_point(pid5, coord5)
        pid6 = num_points+6; call this%add_point(pid6, coord6)
        pid7 = num_points+7; call this%add_point(pid7, coord7)
        pid8 = num_points+8; call this%add_point(pid8, coord8)
        ! Create 12xline
        !<     .---4---.     
        !<    /|      /|     
        !<   7 11    8 12    
        !<  /  |    /  |     
        !< .---3---.___.     
        !< |  /    |  /      
        !< 9 5    10 6       
        !< |/      |/        
        !< .---1---.          
        num_lines = this%lines(this%num_lines)%id
        lid1  = num_lines+1 ; call this%add_line(lid1,  point_ids=[pid1, pid2])
        lid2  = num_lines+2 ; call this%add_line(lid2,  point_ids=[pid3, pid4])
        lid3  = num_lines+3 ; call this%add_line(lid3,  point_ids=[pid5, pid6])
        lid4  = num_lines+4 ; call this%add_line(lid4,  point_ids=[pid7, pid8])
        lid5  = num_lines+5 ; call this%add_line(lid5,  point_ids=[pid1, pid3])
        lid6  = num_lines+6 ; call this%add_line(lid6,  point_ids=[pid2, pid4])
        lid7  = num_lines+7 ; call this%add_line(lid7,  point_ids=[pid5, pid7])
        lid8  = num_lines+8 ; call this%add_line(lid8,  point_ids=[pid6, pid8])
        lid9  = num_lines+9 ; call this%add_line(lid9,  point_ids=[pid1, pid5])
        lid10 = num_lines+10; call this%add_line(lid10, point_ids=[pid2, pid6])
        lid11 = num_lines+11; call this%add_line(lid11, point_ids=[pid3, pid7])
        lid12 = num_lines+12; call this%add_line(lid12, point_ids=[pid4, pid8])
        ! Create 6xquads
        !<      .-------.
        !<     /|      /|
        !<    / | 2   / |<6
        !< 3>/  |    /  |
        !<  .-------. 4 .
        !<  |  /    |  /
        !   | / 5   | /
        !<  |/      |/
        !<  .-------.                    
        !<      ^ 
        !<      1
        num_surfaces =  this%num_surfaces
        sid1 = num_surfaces+1; call this%add_surface(sid1, 4, [lid1, lid2, lid5, lid6],  [1,1,1,1])
        sid2 = num_surfaces+2; call this%add_surface(sid2, 4, [lid3, lid4, lid7, lid8],  [1,1,1,1])
        sid3 = num_surfaces+3; call this%add_surface(sid3, 4, [lid5, lid6, lid9, lid11], [1,1,1,1])
        sid4 = num_surfaces+4; call this%add_surface(sid4, 4, [lid6, lid8, lid10,lid12], [1,1,1,1])
        sid5 = num_surfaces+5; call this%add_surface(sid5, 4, [lid1, lid3, lid9, lid10], [1,1,1,1])
        sid6 = num_surfaces+6; call this%add_surface(sid6, 4, [lid2, lid4, lid11,lid12], [1,1,1,1])
        ! Create hexa
        this%num_volumes = this%num_volumes+1
        call this%volumes(this%num_volumes)%create(volume_id=this%num_volumes,                   &
                                                   number_surfaces=6,                            &
                                                   surfaces_ids=[sid1,sid2,sid3,sid4,sid5,sid6], &
                                                   surfaces_orientation=[1,1,1,1,1,1])
        call this%volume_index%put(key=this%volumes(this%num_volumes)%id, val=this%num_volumes,stat=error)
    end subroutine geometry_add_hexa


    subroutine geometry_add_linear_vol_from_surf_ids(this, volume_id, number_faces, surfaces_ids, surfaces_orientation)
    !-----------------------------------------------------------------
    !< Create a linear volume given the information about its surfaces
    !-----------------------------------------------------------------
        class(geometry_t),         intent(inout) :: this
        integer(ip),               intent(in)    :: volume_id
        integer(ip),               intent(in)    :: number_faces
        integer(ip),               intent(in)    :: surfaces_ids(number_faces)
        integer(ip),               intent(in)    :: surfaces_orientation(number_faces)
        integer(ip)                              :: i
        integer(ip)                              :: error
    !-----------------------------------------------------------------
        assert(allocated(this%surfaces) .and. size(this%surfaces)>= this%num_surfaces+1)
        this%num_volumes = this%num_volumes+1
        call this%volumes(this%num_volumes)%create(volume_id=volume_id,          &
                                                   number_surfaces=number_faces, &
                                                   surfaces_ids=surfaces_ids,    &
                                                   surfaces_orientation=surfaces_orientation)
        call this%volume_index%put(key=this%volumes(this%num_volumes)%id, val=volume_id,stat=error)
    end subroutine geometry_add_linear_vol_from_surf_ids


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


end module geometry_names
