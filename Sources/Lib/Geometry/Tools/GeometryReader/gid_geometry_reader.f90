module gid_geometry_reader_names

    use types_names
    use stdio_names
    use memor_names
    use geometry_names
    use FPL

implicit none
private
#include "debug.i90"

    integer(ip), parameter :: gid_line_len = 256

    type :: gid_geometry_reader_t
    contains
    private
        procedure, non_overridable         :: read_and_fill            => gid_geometry_reader_read_and_fill
        procedure, non_overridable         :: read_point               => gid_geometry_reader_read_point
        procedure, non_overridable         :: read_line                => gid_geometry_reader_read_line
        procedure, non_overridable         :: read_surface             => gid_geometry_reader_read_surface
        procedure, non_overridable         :: read_volume              => gid_geometry_reader_read_volume
        procedure, non_overridable         :: read_point_conditions    => gid_geometry_reader_read_point_conditions
        procedure, non_overridable         :: read_line_conditions     => gid_geometry_reader_read_line_conditions
        procedure, non_overridable         :: read_surface_conditions  => gid_geometry_reader_read_surface_conditions
        procedure, non_overridable         :: read_volume_conditions   => gid_geometry_reader_read_volume_conditions
        procedure,                  public :: fill_geometry => gid_geometry_reader_fill_geometry
    end type

    type(gid_geometry_reader_t), save :: the_GiD_geometry_reader
    
public :: the_GiD_geometry_reader

contains

    subroutine geometry_compose_name ( prefix, name ) 
    !-----------------------------------------------------------------
    ! compose the filename of the GiD geometry given its prefix
    !-----------------------------------------------------------------
        character(len=*)             , intent(in)    :: prefix 
        character(len=:), allocatable, intent(inout) :: name
    !-----------------------------------------------------------------
        name = trim(prefix) // '.txt'
    end subroutine geometry_compose_name

    subroutine move_forward_to_find_string(unit, string, stopstring, line, position)
    !-----------------------------------------------------------------
    !< Move forward line by line in order to find the given string
    !< and return the position of the string in the line
    !-----------------------------------------------------------------
        integer(ip),             intent(in)    :: unit
        character(*),            intent(in)    :: string
        character(*),            intent(in)    :: stopstring
        character(gid_line_len), intent(inout) :: line
        integer(ip),             intent(out)   :: position
        integer(ip)                            :: error
    !-----------------------------------------------------------------
        position = index(string=line, substring=string, back=.false., kind=ip)
        do while(position==0 .and. index(string=line, substring=stopstring, back=.false., kind=ip)==0)
            read(unit=unit,fmt='(a)',iostat=error) line
            if(IS_IOSTAT_END(error) .or. IS_IOSTAT_EOR(error)) exit
            position = index(string=line, substring=string, back=.false., kind=ip)
        end do
    end subroutine


    subroutine gid_geometry_reader_read_point(this, unit, geometry)
    !-----------------------------------------------------------------
    !< Read a point
    !-----------------------------------------------------------------
        class(gid_geometry_reader_t), intent(in)    :: this
        integer(ip),                  intent(in)    :: unit
        type(geometry_t),             intent(inout) :: geometry
        integer(ip)                                 :: i, pos
        character(gid_line_len)                     :: string
        integer(ip)                                 :: id
        real(rp)                                    :: coord(SPACE_DIM)
    !-----------------------------------------------------------------
        read(unit,'(a)') string
        ! Read point ID
        call move_forward_to_find_string(unit,'Num:', 'END', string, pos); assert(pos/=0)
        read(string(pos+4:),*) id
        ! Read point coordinates
        call move_forward_to_find_string(unit, 'Coord:', 'END', string, pos); assert(pos/=0)
        read(string(pos+6:),*) (coord(i),i=1,SPACE_DIM)
        call geometry%add_point(id, coord)
    end subroutine gid_geometry_reader_read_point


    subroutine gid_geometry_reader_read_line(this, unit, geometry)
    !-----------------------------------------------------------------
    !< Read a line
    !-----------------------------------------------------------------
        class(gid_geometry_reader_t), intent(in)    :: this
        integer(ip),                  intent(in)    :: unit
        type(geometry_t),             intent(inout) :: geometry
        character(gid_line_len)                     :: string
        integer(ip)                                 :: id
        integer(ip)                                 :: point(2)
        integer(ip)                                 :: n
        integer(ip)                                 :: p
        real(rp), allocatable                       :: control_points(:)
        real(rp), allocatable                       :: knots(:)
        integer(ip)                                 :: i, j, pos, num, error
    !-----------------------------------------------------------------
        read(unit,'(a)') string
        ! Read line ID
        call move_forward_to_find_string(unit, 'Num:', 'END', string, pos); assert(pos/=0)
        read(string(pos+4:),*) id
        ! Read line points
        call move_forward_to_find_string(unit, 'Points:', 'END', string, pos); assert(pos/=0)
        read(string(pos+7:),*) point(1), point(2)
        ! Read number of control points
        call move_forward_to_find_string(unit, 'Number of Control Points=', 'END', string, pos)
        if(pos == 0) then ! if has not control points is a STLINE
            call geometry%add_line(id, point)
        else              ! if has not control points is a STLINE
            read(string(pos+25:),*) n
            call memalloc(n*(SPACE_DIM+1), control_points, __FILE__, __LINE__)

            ! Read NURBS degree
            call move_forward_to_find_string(unit, 'Degree=', 'END', string, pos); assert(pos/=0)
            read(string(pos+7:),*) p

            ! Read NURBS control points coords
            do i=1, n
                read(unit,'(a)') string
                call move_forward_to_find_string(unit, 'Point', 'END', string, pos); assert(pos/=0)
                read(string(pos+5:),*) num
                call move_forward_to_find_string(unit, 'coords:', 'END', string, pos); assert(pos/=0)
                read(string(pos+7:),*) (control_points((SPACE_DIM+1)*(num-1)+j), j=1, SPACE_DIM)
            end do

            ! Read NURBS number of knots
            call move_forward_to_find_string(unit, 'Number of knots=', 'END', string, pos); assert(pos/=0)
            read(string(pos+16:),*) i
            assert(i==n+p+1)

            ! Read NURBS knots
            call memalloc(n+p+1, knots, __FILE__, __LINE__)
            do i=1, n+p+1
                read(unit,'(a)') string
                call move_forward_to_find_string(unit, 'knot', 'END', string, pos); assert(pos/=0)
                read(string(pos+4:),*) num
                call move_forward_to_find_string(unit, 'value=', 'END', string, pos); assert(pos/=0)
                read(string(pos+6:),*) knots(num)
            end do

            ! Read NURBS Weights
            call move_forward_to_find_string(unit, 'Rational weights:', 'END', string, pos)
            if(pos == 0) then
                do i=1, n
                    control_points((SPACE_DIM+1)*i) = 1._rp
                end do
            else
                do i=1, n
                    read(unit,'(a)') string
                    read(string,*) control_points((SPACE_DIM+1)*i)
                end do
            endif

            ! Define 4D control points (i.e. multiply by weights)
            do i=1, n
                do j=1,SPACE_DIM
                    control_points((SPACE_DIM+1)*(i-1)+j) = control_points((SPACE_DIM+1)*(i-1)+j) * control_points((SPACE_DIM+1)*i)
                end do
            end do

            call geometry%add_line(id, point, n, p, control_points, knots)
            call memfree(control_points, __FILE__, __LINE__)
            call memfree(knots, __FILE__, __LINE__)
        endif

    end subroutine gid_geometry_reader_read_line


    subroutine gid_geometry_reader_read_surface(this, unit, geometry)
    !-----------------------------------------------------------------
    !< Read a surface
    !-----------------------------------------------------------------
        class(gid_geometry_reader_t), intent(in)    :: this
        integer(ip),                  intent(in)    :: unit
        type(geometry_t),             intent(inout) :: geometry
        integer(ip)                                 :: id
        integer(ip)                                 :: num_lines
        integer(ip), allocatable                    :: lines_ids(:)
        integer(ip), allocatable                    :: lines_orientation(:)
        integer(ip)                                 :: nu
        integer(ip)                                 :: nv
        integer(ip)                                 :: pu
        integer(ip)                                 :: pv
        real(rp),    allocatable                    :: control_points(:)
        real(rp),    allocatable                    :: u_knots(:)
        real(rp),    allocatable                    :: v_knots(:)
        character(gid_line_len)                     :: string
        character(gid_line_len)                     :: tmpstring
        integer(ip)                                 :: i, j, k, pos, num_knots, error, idu, idv, counter
    !-----------------------------------------------------------------
        read(unit,'(a)') string
        ! Read surface ID
        call move_forward_to_find_string(unit, 'Num:', 'END', string, pos); assert(pos/=0)
        read(string(pos+4:),*) id
        ! Read surface number of lines
        call move_forward_to_find_string(unit, 'NumLines:', 'END', string, pos); assert(pos/=0)
        read(string(pos+9:),*) num_lines
        call memalloc(num_lines, lines_ids, __FILE__, __LINE__)
        call memalloc(num_lines, lines_orientation, __FILE__, __LINE__)

        ! Read surface lines ID's
        do i=1, num_lines
            read(unit, '(a)') string
            call move_forward_to_find_string(unit, 'Line:', 'END', string, pos); assert(pos/=0)
            read(string(pos+5:),*) lines_ids(i)
            call move_forward_to_find_string(unit, 'Orientation:', 'END', string, pos); assert(pos/=0)
            read(string(pos+12:),*) tmpstring
            if(trim(adjustl(tmpString)) == 'SAME1ST') then
                lines_orientation(i) = 1
            elseif(trim(adjustl(tmpString)) == 'DIFF1ST') then
                lines_orientation(i) = -1
            else
                assert(.false.)
            endif
        enddo

        ! Read number of control points
        call move_forward_to_find_string(unit, 'Number of Control Points=', 'END', string, pos)
        if(pos == 0) then ! if it has NOT control points is a STSURFACE
            call geometry%add_surface(id, num_lines, lines_ids, lines_orientation)
        else              ! if it has control points is a NURBS SURFACE
            read(string(pos+25:),*) nu, nv
            call memalloc(nu*nv*(SPACE_DIM+1), control_points, __FILE__, __LINE__)
            call move_forward_to_find_string(unit, 'Number of Control Points=', 'END', string, pos)
            ! Read NURBS degree
            call move_forward_to_find_string(unit, 'Degree=', 'END', string, pos); assert(pos/=0)
            read(string(pos+7:),*) pu, pv

            counter = 0
            do i=1, nu
                do j=1, nv
                    read(unit, '(a)') string
                    call move_forward_to_find_string(unit, 'Point', 'END', string, pos)
                    read(string(pos+5:),*) idu, idv
                    call move_forward_to_find_string(unit, 'coords:', 'END', string, pos); assert(pos/=0)
                    read(string(pos+7:),*) (control_points(counter*(SPACE_DIM+1)+k),k=1,SPACE_DIM)
                    counter = counter + 1
                enddo
            enddo

            call move_forward_to_find_string(unit, 'Number of knots in U=', 'END', string, pos)
            read(string(pos+21:),*) num_knots
            assert(num_knots == nu+pu+1)

            ! Read NURBS knots
            call memalloc(nu+pu+1, u_knots, __FILE__, __LINE__)
            do i=1, nu+pu+1
                read(unit,'(a)') string
                call move_forward_to_find_string(unit, 'knot', 'END', string, pos); assert(pos/=0)
                read(string(pos+4:),*) k
                call move_forward_to_find_string(unit, 'value=', 'END', string, pos); assert(pos/=0)
                read(string(pos+6:),*) u_knots(k)
            end do

            call move_forward_to_find_string(unit, 'Number of knots in V=', 'END', string, pos)
            read(string(pos+21:),*) num_knots
            assert(num_knots == nv+pv+1)

            ! Read NURBS knots
            call memalloc(nv+pv+1, v_knots, __FILE__, __LINE__)
            do i=1, nv+pv+1
                read(unit,'(a)') string
                call move_forward_to_find_string(unit, 'knot', 'END', string, pos); assert(pos/=0)
                read(string(pos+4:),*) k
                call move_forward_to_find_string(unit, 'value=', 'END', string, pos); assert(pos/=0)
                read(string(pos+6:),*) v_knots(k)
            end do

            ! Read NURBS Weights
            call move_forward_to_find_string(unit, 'Rational weights:', 'END', string, pos)
            if(pos == 0) then ! B-Spline Rational weights by default
                do i=1, nu*nv
                    control_points((SPACE_DIM+1)*i) = 1.0_rp 
                enddo
            else
                do i=1, nu*nv
                    read(unit,'(a)') string
                    read(string,*) control_points((SPACE_DIM+1)*i)
                end do
            ENDIF
            call geometry%add_surface(id, num_lines, lines_ids, lines_orientation, nu, nv, pu, pv, control_points, u_knots, v_knots)
            call memfree(control_points, __FILE__, __LINE__)
            call memfree(u_knots, __FILE__, __LINE__)
            call memfree(v_knots, __FILE__, __LINE__)
        endif
        call memfree(lines_ids, __FILE__, __LINE__)
        call memfree(lines_orientation, __FILE__, __LINE__)

    end subroutine gid_geometry_reader_read_surface


    subroutine gid_geometry_reader_read_volume(this, unit, geometry)
    !-----------------------------------------------------------------
    !< Read a volume
    !-----------------------------------------------------------------
        class(gid_geometry_reader_t), intent(in)    :: this
        integer(ip)  ,                intent(in)    :: unit
        type(geometry_t),             intent(inout) :: geometry
        integer(ip)                                 :: id
        integer(ip)                                 :: num_surfaces
        integer(ip), allocatable                    :: surfaces_ids(:)
        integer(ip), allocatable                    :: surfaces_orientation(:)
        character(gid_line_len)                     :: string
        character(gid_line_len)                     :: tmpstring
        integer(ip)                                 :: i, j, k, pos, num_knots, error, nu, nv, counter
    !-----------------------------------------------------------------
        read(unit,'(a)') string
        ! Read volume ID
        call move_forward_to_find_string(unit, 'Num:', 'END', string, pos); assert(pos/=0)
        read(string(pos+4:),*) id
        ! Read volume number of surfaces
        call move_forward_to_find_string(unit, 'NumSurfaces:', 'END', string, pos); assert(pos/=0)
        read(string(pos+12:),*) num_surfaces
        call memalloc(num_surfaces, surfaces_ids, __FILE__, __LINE__)
        call memalloc(num_surfaces, surfaces_orientation, __FILE__, __LINE__)

        ! Read volume surfaces ID's
        do i=1, num_surfaces
            read(unit, '(a)') string
            call move_forward_to_find_string(unit, 'Surface:', 'END', string, pos); assert(pos/=0)
            read(string(pos+8:),*) surfaces_ids(i)
            call move_forward_to_find_string(unit, 'Orientation:', 'END', string, pos); assert(pos/=0)
            read(string(pos+12:),*) tmpstring
            if(trim(adjustl(tmpString)) == 'SAME1ST') then
                surfaces_orientation(i) = 1
            elseif(trim(adjustl(tmpString)) == 'DIFF1ST') then
                surfaces_orientation(i) = -1
            else
                assert(.false.)
            endif
        enddo

        call geometry%add_volume(id, num_surfaces, surfaces_ids, surfaces_orientation)
        call memfree(surfaces_ids, __FILE__, __LINE__)
        call memfree(surfaces_orientation, __FILE__, __LINE__)
    end subroutine gid_geometry_reader_read_volume


    subroutine gid_geometry_reader_read_point_conditions(this, unit, geometry)
    !-----------------------------------------------------------------
    !< Read point conditions
    !-----------------------------------------------------------------
        class(gid_geometry_reader_t), intent(in)    :: this
        integer(ip),                  intent(in)    :: unit
        type(geometry_t),             intent(inout) :: geometry
        type(geometric_point_t), pointer            :: point 
        integer(ip)                                 :: i, pos
        character(gid_line_len)                     :: string
        integer(ip)                                 :: id, set_id
    !-----------------------------------------------------------------
        call move_forward_to_find_string(unit,'END CONDITION', 'END INTERVAL DATA', string, pos); assert(pos/=0)
        call move_forward_to_find_string(unit,'Geometry Entities', 'END INTERVAL DATA', string, pos); assert(pos/=0)
        do while(.true.)
            read(unit, '(a)') string
            call move_forward_to_find_string(unit,'Conds:', 'End Geometry Entities', string, pos)
            if(pos == 0) exit
            read(string(:pos),*) id
            read(string(pos+6:),*) set_id
            point => geometry%get_point(id)
            call point%set_condition(set_id)
        enddo
    end subroutine gid_geometry_reader_read_point_conditions


    subroutine gid_geometry_reader_read_line_conditions(this, unit, geometry)
    !-----------------------------------------------------------------
    !< Read point conditions
    !-----------------------------------------------------------------
        class(gid_geometry_reader_t), intent(in)    :: this
        integer(ip),                  intent(in)    :: unit
        type(geometry_t),             intent(inout) :: geometry
        type(line_t), pointer                       :: line
        integer(ip)                                 :: i, pos
        character(gid_line_len)                     :: string
        integer(ip)                                 :: id, set_id
    !-----------------------------------------------------------------
        call move_forward_to_find_string(unit,'END CONDITION', 'END INTERVAL DATA', string, pos); assert(pos/=0)
        call move_forward_to_find_string(unit,'Geometry Entities', 'END INTERVAL DATA', string, pos); assert(pos/=0)
        do while(.true.)
            read(unit, '(a)') string
            call move_forward_to_find_string(unit,'Conds:', 'End Geometry Entities', string, pos)
            if(pos == 0) exit
            read(string(:pos),*) id
            read(string(pos+6:),*) set_id
            line => geometry%get_line(id)
            call line%set_condition(set_id)
        enddo
    end subroutine gid_geometry_reader_read_line_conditions


    subroutine gid_geometry_reader_read_surface_conditions(this, unit, geometry)
    !-----------------------------------------------------------------
    !< Read point conditions
    !-----------------------------------------------------------------
        class(gid_geometry_reader_t), intent(in)    :: this
        integer(ip),                  intent(in)    :: unit
        type(geometry_t),             intent(inout) :: geometry
        type(surface_t), pointer                    :: surface
        integer(ip)                                 :: i, pos
        character(gid_line_len)                     :: string
        integer(ip)                                 :: id, set_id
    !-----------------------------------------------------------------
        call move_forward_to_find_string(unit,'END CONDITION', 'END INTERVAL DATA', string, pos); assert(pos/=0)
        call move_forward_to_find_string(unit,'Geometry Entities', 'END INTERVAL DATA', string, pos); assert(pos/=0)
        do while(.true.)
            read(unit, '(a)') string
            call move_forward_to_find_string(unit,'Conds:', 'End Geometry Entities', string, pos)
            if(pos == 0) exit
            read(string(:pos),*) id
            read(string(pos+6:),*) set_id
            surface => geometry%get_surface(id)
            call surface%set_condition(set_id)
        enddo
    end subroutine gid_geometry_reader_read_surface_conditions


    subroutine gid_geometry_reader_read_volume_conditions(this, unit, geometry)
    !-----------------------------------------------------------------
    !< Read point conditions
    !-----------------------------------------------------------------
        class(gid_geometry_reader_t), intent(in)    :: this
        integer(ip),                  intent(in)    :: unit
        type(geometry_t),             intent(inout) :: geometry
        type(volume_t), pointer                     :: volume
        integer(ip)                                 :: i, pos
        character(gid_line_len)                     :: string
        integer(ip)                                 :: id, set_id
    !-----------------------------------------------------------------
        call move_forward_to_find_string(unit,'END CONDITION', 'END INTERVAL DATA', string, pos); assert(pos/=0)
        call move_forward_to_find_string(unit,'Geometry Entities', 'END INTERVAL DATA', string, pos); assert(pos/=0)
        do while(.true.)
            read(unit, '(a)') string
            call move_forward_to_find_string(unit,'Conds:', 'End Geometry Entities', string, pos)
            if(pos == 0) exit
            read(string(:pos),*) id
            read(string(pos+6:),*) set_id
            volume => geometry%get_volume(id)
            call volume%set_condition(set_id)
        enddo
    end subroutine gid_geometry_reader_read_volume_conditions


    subroutine gid_geometry_reader_fill_geometry(this, parameter_list, geometry)
    !------------------------------------------------------------------------
    !< Open the unit of the GiD file and fill the geometry
    !------------------------------------------------------------------------
        class(gid_geometry_reader_t), intent(in)    :: this
        type(ParameterList_t),        intent(in)    :: parameter_list
        type(geometry_t),             intent(inout) :: geometry
        integer(ip)                                 :: istat
        character(len=:), allocatable               :: dir_path
        character(len=:), allocatable               :: prefix
        character(len=:), allocatable               :: name
        integer(ip)                                 :: lunio
    !------------------------------------------------------------------------
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
        call this%read_and_fill(lunio, geometry)
        call io_close(lunio)
    end subroutine gid_geometry_reader_fill_geometry


    subroutine gid_geometry_reader_read_and_fill(this, unit, geometry)
    !------------------------------------------------------------------------
    !< This routine reads a geometry writen by GiD using data report.
    !------------------------------------------------------------------------
        class(gid_geometry_reader_t), intent(in)    :: this
        integer(ip),                  intent(in)    :: unit
        type(geometry_t),             intent(inout) :: geometry
        character(gid_line_len)                     :: string
        integer(ip)                                 :: num_points, num_lines, num_surfaces, num_volumes
        integer(ip)                                 :: i, pos, istat
    !------------------------------------------------------------------------
        call geometry%free()
        num_points   = 0
        num_lines    = 0
        num_surfaces = 0
        num_volumes  = 0
        ! Count points, lines, surfaces and volumes
        read(unit,'(a)') string
        do while(string(1:12)/='END ENTITIES')
           read(unit,'(a)') string
           if(string(1:5)=='POINT') then
              num_points   = num_points+1
           else if(string(1:6)=='STLINE'.or.string(1:8)=='NURBLINE') then
              num_lines    = num_lines+1
           else if(string(1:11)=='NURBSURFACE') then
              num_surfaces = num_surfaces+1
           else if(string(1:6)=='VOLUME')       then
              num_volumes  = num_volumes+1
           end if
        end do
        call io_rewind(unit)

        call geometry%create(num_points, num_lines, num_surfaces, num_volumes)

        read(unit,'(a)') string
        do while(string(1:12)/='END ENTITIES')
           read(unit,'(a)') string
           if(string(1:5)=='POINT') then
              call this%read_point(unit, geometry)
           else if(string(1:6)=='STLINE'.or.string(1:8)=='NURBLINE') then
              call this%read_line(unit, geometry)
           else if(string(1:11)=='NURBSURFACE') then
              call this%read_surface(unit, geometry)
           else if(string(1:6)=='VOLUME')       then
              call this%read_volume(unit, geometry)
           end if
        end do

        read(unit,'(a)') string
       call move_forward_to_find_string(unit,'CONDITION:', 'END INTERVAL DATA', string, pos)
        do while(pos/=0)
           read(unit,'(a)') string
           call move_forward_to_find_string(unit,'CONDTYPE:', 'END INTERVAL DATA', string, pos)
           if(trim(adjustl(string(pos+9:)))=='over points') then
               call this%read_point_conditions(unit, geometry)
           else if(trim(adjustl(string(pos+9:)))=='over lines') then
               call this%read_line_conditions(unit, geometry)
           else if(trim(adjustl(string(pos+9:)))=='over surfaces') then
               call this%read_surface_conditions(unit, geometry)
           else if(trim(adjustl(string(pos+9:)))=='over volumes')       then
               call this%read_volume_conditions(unit, geometry)
           end if
           call move_forward_to_find_string(unit,'CONDITION:', 'END INTERVAL DATA', string, pos)
        end do

        call geometry%init()
    end subroutine gid_geometry_reader_read_and_fill


end module gid_geometry_reader_names
