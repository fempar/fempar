
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

module vtk_handler_names

USE iso_c_binding
USE lib_vtk_io
USE types_names
USE memor_names
USE field_metadata
USE vtk_mesh_and_field_generator
USE iso_fortran_env,             only: error_unit
USE ir_precision,                only: str
USE environment_names,           only: environment_t
USE fe_space_names,              only: serial_fe_space_t
USE fe_function_names,           only: fe_function_t
USE vtk_parameters_names

implicit none
#include "debug.i90"

private
  
    interface
        function mkdir_recursive(path) bind(c,name="mkdir_recursive")
            use iso_c_binding
            integer(kind=c_int) :: mkdir_recursive
            character(kind=c_char,len=1), intent(IN) :: path(*)
        end function mkdir_recursive
    end interface  

    ! File extensions and time prefix
    character(len=5) :: time_prefix = 'time_'
    character(len=4) :: vtk_ext     = '.vtu'
    character(len=4) :: pvd_ext     = '.pvd'
    character(len=5) :: pvtu_ext    = '.pvtu'

    !-----------------------------------------------------------------
    ! State transition diagram for type(vtk_handler_t)
    !-----------------------------------------------------------------
    ! Input State         | Action                | Output State 
    !-----------------------------------------------------------------
    ! START               | Initialize            | INITIALIZED
    ! START               | Free                  | START

    ! INITIALIZED         | Open                  | OPEN 
    ! INITIALIZED         | Write_pvtu            | INITIALIZED
    ! INITIALIZED         | Write_pvd             | INITIALIZED
    ! INITIALIZED         | Free                  | START

    ! OPEN                | Write_mesh            | GEO_OPEN   
    ! OPEN                | Write_pvtu            | OPEN         
    ! OPEN                | Write_pvd             | OPEN       
    ! OPEN                | Close                 | CLOSE
    ! OPEN                | Free                  | START

    ! GEO_OPEN            | Write_node_field      | POINTDATA_OPEN
    ! GEO_OPEN            | Write_pvtu            | GEO_OPEN
    ! GEO_OPEN            | Write_pvd             | GEO_OPEN
    ! GEO_OPEN            | Close                 | CLOSE
    ! GEO_OPEN            | Free                  | START

    ! POINTDATA_OPEN      | Write_node_field      | POINTDATA_OPEN
    ! POINTDATA_OPEN      | Write_pvtu            | POINTDATA_OPEN
    ! POINTDATA_OPEN      | Write_pvd             | POINTDATA_OPEN
    ! POINTDATA_OPEN      | Close                 | CLOSE
    ! POINTDATA_OPEN      | Free                  | START

    ! CLOSE               | Open                  | OPEN
    ! CLOSE               | Write_pvtu            | CLOSE        
    ! CLOSE               | Write_pvd             | CLOSE
    ! CLOSE               | Free                  | START

    ! STATE PARAMETERS
    integer(ip), parameter :: vtk_handler_state_start                 = 0
    integer(ip), parameter :: vtk_handler_state_initialized           = 1
    integer(ip), parameter :: vtk_handler_state_write_open            = 2
    integer(ip), parameter :: vtk_handler_state_write_geo_open        = 3
    integer(ip), parameter :: vtk_handler_state_write_pointdata_open  = 4
    integer(ip), parameter :: vtk_handler_state_write_pointdata_close = 5
    integer(ip), parameter :: vtk_handler_state_write_geo_close       = 6
    integer(ip), parameter :: vtk_handler_state_write_close           = 7


    ! Type for storing several mesh data with its field descriptors
    ! It also contains information about the number of parts (PVTU) and time steps (PVD)
    ! It stores the directory path and the prefix where to write in disk
    type vtk_handler_t
    private
        class(environment_t),       pointer :: env      => NULL() ! Poins to fe_space_t
        character(len=:), allocatable       :: path               ! Output path
        character(len=:), allocatable       :: prefix             ! XVTX filename prefix
        integer(ip)                         :: file_id            ! VTU file descriptor
        type(vtk_mesh_and_field_generator_t):: mesh               ! mesh data 
        type(field_metadata_t), allocatable :: field(:)           ! field metadata
        real(rp),               allocatable :: steps(:)           ! Array of parameters (time, eigenvalues,etc.)
        integer(ip)                         :: steps_counter = 0  ! time steps counter
        integer(ip)                         :: num_meshes    = 0  ! Number of VTK meshes stored
        integer(ip)                         :: num_steps     = 0  ! Number of time steps
        integer(ip)                         :: num_parts     = 0  ! Number of parts
        integer(ip)                         :: root_task     = 0  ! Root processor
        integer(ip)                         :: state         = vtk_handler_state_start
    contains
    private
        procedure, public :: create                           => vtk_handler_create
        procedure, public :: open_vtu                         => vtk_handler_open_vtu
        procedure, public :: write_vtu_mesh                   => vtk_handler_write_vtu_mesh
        procedure, public :: write_vtu_node_field             => vtk_handler_write_vtu_node_field
        procedure, public :: close_vtu                        => vtk_handler_close_vtu
        procedure, public :: write_pvtu                       => vtk_handler_write_pvtu
        procedure, public :: write_pvd                        => vtk_handler_write_pvd
        procedure, public :: free                             => vtk_handler_free
        procedure         :: get_vtk_time_output_path         => vtk_handler_get_vtk_time_output_path
        procedure         :: get_pvd_time_output_path         => vtk_handler_get_pvd_time_output_path
        procedure         :: get_vtk_filename                 => vtk_handler_get_vtk_filename
        procedure         :: get_pvtu_filename                => vtk_handler_get_pvtu_filename
        procedure         :: set_root_task                    => vtk_handler_set_root_task
        procedure         :: set_num_steps                    => vtk_handler_set_num_steps
        procedure         :: set_num_parts                    => vtk_handler_set_num_parts
        procedure         :: append_step                      => vtk_handler_append_step
        procedure         :: create_directory                 => vtk_handler_create_dir_hierarchy_on_root_task
    end type vtk_handler_t

public :: vtk_handler_t

contains


    function vtk_handler_create_dir_hierarchy_on_root_task(this, path, issue_final_barrier) result(res)
    !-----------------------------------------------------------------
    !< The root process create a hierarchy of directories
    !-----------------------------------------------------------------
        class(vtk_handler_t), intent(INOUT) :: this
        character(len=*),     intent(IN)    :: path
        logical, optional,    intent(IN)    :: issue_final_barrier
        logical                             :: ifb
        integer(kind=c_int)                 :: res
    !-----------------------------------------------------------------
        ifb = .False.
        if(present(issue_final_barrier)) ifb = issue_final_barrier
        assert(associated(this%env))

        res=0

        if(this%env%am_i_l1_root()) then
            res = mkdir_recursive(path//C_NULL_CHAR)
            check ( res == 0 ) 
        end if

        if(ifb) call this%env%l1_barrier()
    end function vtk_handler_create_dir_hierarchy_on_root_task


    function vtk_handler_get_VTK_time_output_path(this, time_step) result(dp)
    !-----------------------------------------------------------------
    !< Build time output dir path for the vtk files in each timestep
    !-----------------------------------------------------------------
        class(vtk_handler_t),       intent(INOUT) :: this
        real(rp),         optional, intent(IN)    :: time_step
        character(len=:), allocatable             :: dp
        character(len=:), allocatable             :: fp
        real(rp)                                  :: ts
    !-----------------------------------------------------------------
        ts = this%num_steps
        if(present(time_step)) ts = time_step 
        dp = trim(adjustl(this%path))//'/'//time_prefix//trim(adjustl(str(no_sign=.true., n=ts)))//'/'
    end function vtk_handler_get_VTK_time_output_path


    function vtk_handler_get_PVD_time_output_path(this, path, time_step) result(dp)
    !-----------------------------------------------------------------
    !< Build output dir path for the PVD files
    !-----------------------------------------------------------------
        class(vtk_handler_t),       intent(INOUT) :: this
        character(len=*), optional, intent(IN)    :: path
        real(RP),         optional, intent(IN)    :: time_step
        character(len=:), allocatable             :: dp
        character(len=:), allocatable             :: fp
        real(rp)                                  :: ts
    !-----------------------------------------------------------------
        ts = this%num_steps
        if(present(time_step)) ts = time_step 
        dp = time_prefix//trim(adjustl(str(no_sign=.true., n=ts)))//'/'
    end function vtk_handler_get_PVD_time_output_path


    function vtk_handler_get_VTK_filename(this, part_number) result(fn)
    !-----------------------------------------------------------------
    !< Build VTK filename
    !-----------------------------------------------------------------
        class(vtk_handler_t),       intent(INOUT) :: this
        integer(ip), optional,      intent(IN)    :: part_number
        character(len=:), allocatable             :: fn
        character(len=:), allocatable             :: fp
        integer(ip)                               :: me !, np
    !-----------------------------------------------------------------
        !me = default_root_task
        !np = default_number_of_tasks
        !call this%env%info(me, np)
        me = this%env%get_l1_rank()

        if(present(part_number)) me = part_number
        fn = trim(adjustl(this%prefix))//'_'//trim(adjustl(str(no_sign=.true., n=me)))//vtk_ext
    end function vtk_handler_get_VTK_filename


      ! Build VTK filename
    function vtk_handler_get_pvtu_filename(this, time_step) result(fn)
    !-----------------------------------------------------------------
    !< Build pvtu filename
    !-----------------------------------------------------------------
        class(vtk_handler_t),       intent(INOUT) :: this
        real(rp),         optional, intent(IN)    :: time_step
        character(len=:), allocatable             :: fn
        real(rp)                                  :: ts
    !-----------------------------------------------------------------
        ts = default_step_value
        if(allocated(this%steps)) then
            if(this%steps_counter >0 .and. this%steps_counter <= size(this%steps,1)) &
                ts = this%steps(this%steps_counter)
        endif
        if(present(time_step)) ts = time_step
        fn = trim(adjustl(this%prefix))//'_'//trim(adjustl(str(no_sign=.true., n=ts)))//pvtu_ext
    end function vtk_handler_get_pvtu_filename


    subroutine vtk_handler_set_root_task(this, root)
    !-----------------------------------------------------------------
    !< Set the root processor
    !-----------------------------------------------------------------
        class(vtk_handler_t), intent(INOUT) :: this
        integer(ip),          intent(IN)    :: root
    !-----------------------------------------------------------------
        this%root_task = root
    end subroutine vtk_handler_set_root_task


    subroutine vtk_handler_set_num_steps(this, steps)
    !-----------------------------------------------------------------
    !< Set the number of time steps of the simulation to be writen in the PVD
    !-----------------------------------------------------------------
        class(vtk_handler_t), intent(INOUT) :: this
        integer(ip),          intent(IN)    :: steps
    !-----------------------------------------------------------------
        this%num_steps = steps
        if(.not.allocated(this%steps)) then
            call memalloc ( this%num_steps, this%steps, __FILE__,__LINE__)
        elseif(size(this%steps)<this%num_steps) then
            call memrealloc ( steps, this%steps, __FILE__,__LINE__)
        endif
    end subroutine vtk_handler_set_num_steps


    subroutine vtk_handler_append_step(this, current_step)
    !-----------------------------------------------------------------
    !< Append a new time stepstep
    !-----------------------------------------------------------------
        class(vtk_handler_t), intent(INOUT) :: this
        real(rp),             intent(IN)    :: current_step !Current time step
    !-----------------------------------------------------------------
        this%steps_counter = this%steps_counter + 1
        this%num_steps = max(this%num_steps, this%steps_counter)
        if(.not.allocated(this%steps)) then
            call memalloc ( this%num_steps, this%steps, __FILE__,__LINE__)
        elseif(size(this%steps)<this%num_steps) then
            call memrealloc (max(this%steps_counter, int(size(this%steps)*1.5)), this%steps, __FILE__,__LINE__)
        endif        
        this%steps(this%steps_counter) = current_step
    end subroutine vtk_handler_append_step


    subroutine vtk_handler_set_num_parts(this, number_of_parts)
    !-----------------------------------------------------------------
    !< Set the number of parts of the partitioned mesh to be writen in the PVTU
    !-----------------------------------------------------------------
        class(vtk_handler_t), intent(INOUT) :: this
        integer(ip),          intent(IN)    :: number_of_parts
    !-----------------------------------------------------------------
        this%num_parts = number_of_parts
    end subroutine vtk_handler_set_num_parts


    subroutine vtk_handler_create(this, fe_space, path, prefix, root_task, number_of_steps)
    !-----------------------------------------------------------------
    !< Initialize the vtk_handler_t derived type
    !-----------------------------------------------------------------
        class(vtk_handler_t),             intent(INOUT) :: this
        class(serial_fe_space_t), target, intent(IN)    :: fe_space
        character(len=*),                 intent(IN)    :: path
        character(len=*),                 intent(IN)    :: prefix  
        integer(ip),      optional,       intent(IN)    :: root_task
        integer(ip),      optional,       intent(IN)    :: number_of_steps
        class(environment_t),     pointer               :: environment
        integer(ip)                                     :: me, np, st, rp
    !-----------------------------------------------------------------
        assert(this%state == vtk_handler_state_start)

        ! Default values
        rp = default_root_task
        st = default_guess_number_of_steps

        ! Optional arguments
        if(present(root_task))       rp = root_task
        if(present(number_of_steps)) st = number_of_steps

        this%env => fe_space%get_environment()
        !call this%env%info(me,np) 
        call this%set_root_task(rp)

        if(this%env%am_i_l1_task() ) then
           
            call this%mesh%set_fe_space(fe_space)
            allocate(this%field(this%mesh%get_number_fields()))

            this%path = path
            this%prefix = prefix
            call this%set_num_parts(this%env%get_l1_size())
            call this%set_num_steps(st)
        endif
        this%state = vtk_handler_state_initialized
    end subroutine vtk_handler_create


    function vtk_handler_open_vtu(this, file_name, time_step, format) result(E_IO)
    !-----------------------------------------------------------------
    !< Start the writing of a single VTK file to disk ( VTK_INI_XML)
    !< (only if it's a fine MPI task)
    !-----------------------------------------------------------------
        class(vtk_handler_t),       intent(INOUT) :: this        !< vtk_handler_t derived type
        character(len=*), optional, intent(IN)    :: file_name   !< VTK File NAME
        real(rp),         optional, intent(IN)    :: time_step   !< Time STEP value
        character(len=*), optional, intent(IN)    :: format      !< Ouput ForMaT
        character(len=:), allocatable             :: fn          !< Real File Name
        character(len=:), allocatable             :: dp          !< Real Directory Path
        character(len=:), allocatable             :: of          !< Real Output Format
        real(rp)                                  :: ts          !< Real Time Step
        integer(ip)                               :: me          !< Task identifier
        !integer(ip)                               :: np          !< Number of processors
        integer(ip)                               :: E_IO        !< Error IO
    !-----------------------------------------------------------------
        assert(this%state == vtk_handler_state_initialized .or. this%state == vtk_handler_state_write_close)
        assert(associated(this%env))

        E_IO = 0

        if(this%env%am_i_l1_task() ) then

            assert(allocated(this%path))
            assert(allocated(this%prefix))

            ts = default_step_value
            of = default_vtk_format

            ! Optional arguments
            if(present(time_step)) ts = time_step 
            if(present(file_name)) fn = file_name
            if(present(format)) of = trim(adjustl(format))

            !call this%env%info(me,np) 
            call this%append_step(ts)

            ! Build path
            dp = this%get_VTK_time_output_path(time_step=ts)
            fn = this%get_VTK_filename(part_number=this%env%get_l1_rank()+1)
            fn = dp//fn

            if( this%create_directory(dp,issue_final_barrier=.True.) == 0) then    
                E_IO = VTK_INI_XML(output_format = trim(adjustl(of)),   &
                                   filename = trim(adjustl(fn)),        &
                                   mesh_topology = 'UnstructuredGrid',  &
                                   cf=this%file_id)
            endif
        endif
        this%state = vtk_handler_state_write_open
    end function vtk_handler_open_vtu


    function vtk_handler_write_vtu_mesh(this, fe_function) result(E_IO)
    !-----------------------------------------------------------------
    !< Write VTU mesh (VTK_GEO_XML, VTK_CON_XML )
    !-----------------------------------------------------------------
        class(vtk_handler_t),       intent(INOUT) :: this        !< vtk_handler_t derived type
        type(fe_function_t),        intent(IN)    :: fe_function
        integer(ip)                               :: E_IO        !< Error IO
      ! ----------------------------------------------------------------------------------
        assert(this%state == vtk_handler_state_write_open)
        assert(associated(this%env))

        E_IO = 0

        if(this%env%am_i_l1_task()) then
            call this%mesh%generate_mesh(fe_function)

            E_IO = VTK_GEO_XML(NN = this%mesh%get_number_nodes(),    &
                               NC = this%mesh%get_number_elements(), &
                               X  = this%mesh%get_X_coordinates(),   &
                               Y  = this%mesh%get_Y_coordinates(),   &
                               Z  = this%mesh%get_Z_coordinates(),   &
                               cf = this%file_id)

            E_IO = VTK_CON_XML(NC        = this%mesh%get_number_elements(), &
                               connect   = this%mesh%get_connectivities(),  &
                               offset    = this%mesh%get_offset(),          &
                               cell_type = this%mesh%get_cell_types(),      &
                               cf        = this%file_id)
        endif
        this%state = vtk_handler_state_write_geo_open
    end function vtk_handler_write_vtu_mesh


    function vtk_handler_write_vtu_node_field(this, fe_function, field_id, field_name) result(E_IO)
    !-----------------------------------------------------------------
    !< Write node field to file
    !-----------------------------------------------------------------
        class(vtk_handler_t),       intent(INOUT) :: this           !< vtk_handler_t derived type
        type(fe_function_t),        intent(IN)    :: fe_function    !< fe_function containing the field to be written
        integer(ip),                intent(IN)    :: field_id       !< Field id
        character(len=*),           intent(IN)    :: field_name     !< name of the field
        real(rp), allocatable                     :: field(:,:)     !< Raw data of the field
        integer(ip)                               :: nc             !< number_of_components
        logical                                   :: ft             !< Fine Task
        integer(ip)                               :: E_IO           !< IO Error
    !-----------------------------------------------------------------
        assert(this%state == vtk_handler_state_write_geo_open .or. this%state == vtk_handler_state_write_pointdata_open)
        assert(associated(this%env))

        E_IO = 0

        if(this%env%am_i_l1_task() ) then        
            assert(allocated(this%field))
            assert(field_id>0 .and. field_id<=size(this%field))
            E_IO = this%mesh%generate_field(fe_function, field_id, field_name, field, number_components=nc)
            if(this%state == vtk_handler_state_write_geo_open) then
                E_IO = VTK_DAT_XML(var_location='node',var_block_action='open', cf=this%file_id)
                this%state = vtk_handler_state_write_pointdata_open
            endif
            E_IO = VTK_VAR_XML(NC_NN=this%mesh%get_number_nodes(), N_COL=nc, varname=field_name, var=field, cf=this%file_id)
            call this%field(field_id)%set(name=field_name, data_type='Float64', number_components=nc)
            if(allocated(field)) call memfree(field, __FILE__, __LINE__)
        endif
        this%state = vtk_handler_state_write_pointdata_open
    end function vtk_handler_write_vtu_node_field


    function vtk_handler_close_vtu(this) result(E_IO)
    !-----------------------------------------------------------------
    !< Ends the writing of a single VTK file to disk (if I am fine MPI task)
    !< Closes geometry ( VTK_END_XML, VTK_GEO_XML )
    !-----------------------------------------------------------------
        class(vtk_handler_t),       intent(INOUT) :: this        !< vtk_handler_t derived type
        integer(ip)                               :: E_IO        !< IO Error
    !-----------------------------------------------------------------
        assert(this%state == vtk_handler_state_write_pointdata_open .or. this%state == vtk_handler_state_write_geo_open .or. this%state == vtk_handler_state_write_open)
        assert(associated(this%env))

        E_IO = 0

        if (this%env%am_i_l1_task() ) then
            if(this%state == vtk_handler_state_write_pointdata_open) then
                E_IO = VTK_DAT_XML(var_location='node',var_block_action='close', cf=this%file_id)
                ! this%state = vtk_handler_state_write_pointdata_close
                E_IO = VTK_GEO_XML(cf=this%file_id)
                ! this%state = vtk_handler_state_write_geo_close
                E_IO = VTK_END_XML(cf=this%file_id)
                ! this%state = vtk_handler_state_close
            elseif(this%state == vtk_handler_state_write_geo_open) then
                E_IO = VTK_GEO_XML(cf=this%file_id)
                ! this%state = vtk_handler_state_write_geo_close
                E_IO = VTK_END_XML(cf=this%file_id)
                ! this%state = vtk_handler_state_write_close
            elseif(this%state == vtk_handler_state_write_open) then
                E_IO = VTK_END_XML(cf=this%file_id)
                ! this%state = vtk_handler_state_write_close
            endif
        endif
        this%state = vtk_handler_state_write_close
    end function vtk_handler_close_vtu


    function vtk_handler_write_pvtu(this, file_name, time_step) result(E_IO)
    !-----------------------------------------------------------------
    !< Write the pvtu file containing the number of parts
    !< (only root processor)
    !-----------------------------------------------------------------
        class(vtk_handler_t),       intent(INOUT) :: this
        character(len=*), optional, intent(IN)    :: file_name
        real(rp),         optional, intent(IN)    :: time_step
        integer(ip)                               :: rf
        character(len=:),allocatable              :: fn ,dp
        real(rp)                                  :: ts
        integer(ip)                               :: i, fid, E_IO
        !integer(ip)                               :: me, np
        logical                                   :: isDir
    !-----------------------------------------------------------------
        assert(.not. this%state == vtk_handler_state_start)
        check(associated(this%env))

        E_IO = 0

        !call this%env%info(me,np) 
        !if( this%env%am_i_l1_task() .and. me == this%root_task) then
        if( this%env%am_i_l1_root()) then
            check(allocated(this%path))
            check(allocated(this%prefix))
            check(allocated(this%field))
            if(allocated(this%steps)) then
                if(this%steps_counter >0 .and. this%steps_counter <= size(this%steps,1)) &
                    ts = this%steps(this%steps_counter)
            endif

            ts = default_step_value
            if(present(time_step)) ts = time_step 

            ! Build path
            dp = this%get_VTK_time_output_path(time_step=ts)
            fn = this%get_pvtu_filename(time_step=ts)
            fn = dp//fn
            if(present(file_name)) fn = file_name

            if(this%create_directory(trim(dp), issue_final_barrier=.False.) == 0) then
                E_IO = PVTK_INI_XML(filename = trim(adjustl(fn)), mesh_topology = 'PUnstructuredGrid', tp='Float64', cf=rf)
                do i=0, this%num_parts-1
                    E_IO = PVTK_GEO_XML(source=trim(adjustl(this%get_VTK_filename(part_number=i))), cf=rf)
                enddo
                E_IO = PVTK_DAT_XML(var_location = 'Node', var_block_action = 'OPEN', cf=rf)
                ! Write point data fields
                do i=1, size(this%field)
                    if(this%field(i)%is_filled()) then
                        E_IO = PVTK_VAR_XML(varname = trim(adjustl(this%field(i)%get_name())),      &
                                            tp      = trim(adjustl(this%field(i)%get_data_type())), &
                                            Nc      = this%field(i)%get_number_components(),        &
                                            cf      = rf )
                    endif
                enddo
                E_IO = PVTK_DAT_XML(var_location = 'Node', var_block_action = 'CLOSE', cf=rf)
                E_IO = PVTK_END_XML(cf=rf)
            endif
        endif
    end function vtk_handler_write_pvtu


    function vtk_handler_write_pvd(this, file_name) result(E_IO)
    !-----------------------------------------------------------------
    !< Write the PVD file referencing several PVTU files in a timeline
    !< (only root processor)
    !-----------------------------------------------------------------
        class(vtk_handler_t),       intent(INOUT) :: this
        character(len=*), optional, intent(IN)    :: file_name
        integer(ip)                               :: rf
        character(len=:),allocatable              :: pvdfn, pvtufn ,dp
        integer(ip)                               :: ts, E_IO
        !integer(ip)                               :: me, np
        logical                                   :: isDir
    !-----------------------------------------------------------------
        assert(.not. this%state == vtk_handler_state_start)
        !me = default_root_task
        !np = default_number_of_tasks
        check(associated(this%env))
        !call this%env%info(me,np) 

        E_IO = 0

        !if(this%env%am_i_l1_task() .and. me == this%root_task) then
        if( this%env%am_i_l1_root()) then
            pvdfn = trim(adjustl(this%path))//'/'//trim(adjustl(this%prefix))//pvd_ext
            if(present(file_name)) pvdfn = file_name

            if(this%create_directory(trim(adjustl(this%path)), issue_final_barrier=.False.) == 0) then
                if(allocated(this%steps)) then
                    if(size(this%steps,1) >= min(this%num_steps,this%steps_counter)) then
                        E_IO = PVD_INI_XML(filename=trim(adjustl(pvdfn)),cf=rf)
                        do ts=1, min(this%num_steps,this%steps_counter)
                            dp = this%get_PVD_time_output_path(path=this%path, time_step=this%steps(ts))
                            pvtufn = this%get_pvtu_filename(time_step=this%steps(ts))
                            pvtufn = dp//pvtufn
                            E_IO = PVD_DAT_XML(filename=trim(adjustl(pvtufn)),timestep=this%steps(ts), cf=rf)
                        enddo
                        E_IO = PVD_END_XML(cf=rf)
                    endif
                endif
            endif
        endif
    end function vtk_handler_write_pvd


    subroutine VTK_handler_free (this)
    !-----------------------------------------------------------------
    !< Free the vtk_handler_t derived type
    !-----------------------------------------------------------------
        class(vtk_handler_t), intent(inout) :: this
        integer(ip)                         :: i, j 
    !-----------------------------------------------------------------
        call this%mesh%free()
        if(allocated(this%steps))  call memfree(this%steps, __FILE__,__LINE__)
        if(allocated(this%path))   deallocate(this%path)
        if(allocated(this%prefix)) deallocate(this%prefix)
        if(allocated(this%field)) then
            do i=1, size(this%field)
                call this%field(i)%free()
            enddo
            deallocate(this%field)
        endif
        nullify(this%env)
        this%num_meshes = 0
        this%num_steps = 0
        this%num_parts = 0
        this%root_task = 0
        this%state = vtk_handler_state_start
    end subroutine VTK_handler_free


end module vtk_handler_names
