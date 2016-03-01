
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

module lib_vtk_io_interface_names

USE iso_c_binding
USE types_names
USE memor_names
USE ir_precision,               only: str
USE environment_names,          only: environment_t
USE triangulation_names,        only: triangulation_t
USE serial_fe_space_names,      only: serial_fe_space_t
USE Lib_VTK_IO
USE reference_fe_names,         only: topology_quad, topology_tet, fe_type_lagrangian


implicit none
# include "debug.i90"

private

    ! STATE PARAMETERS
    integer(ip), parameter :: VTK_STATE_UNKNOWN          = 0
    integer(ip), parameter :: VTK_STATE_WRITE_STARTED    = 1
    integer(ip), parameter :: VTK_STATE_POINTDATA_OPENED = 2
    integer(ip), parameter :: VTK_STATE_POINTDATA_CLOSED = 3
    integer(ip), parameter :: VTK_STATE_ENDED            = 4
    integer(ip), parameter :: VTK_STATE_READ_STARTED     = 5
  
    interface
        function mkdir_recursive(path) bind(c,name="mkdir_recursive")
            use iso_c_binding
            integer(kind=c_int) :: mkdir_recursive
            character(kind=c_char,len=1), intent(IN) :: path(*)
        end function mkdir_recursive
    end interface

    ! Type for storing field descriptors
    type vtk_field_t
    private
        character(len=:), allocatable :: var_location   ! 'Node' or 'Cell' field
        character(len=:), allocatable :: var_name       ! Name of the field
        character(len=:), allocatable :: field_type     ! Field data type 'Float32', 'Float64', 'Int32', etc.
        integer(ip)                   :: num_comp = 0   ! Number of components
        logical                       :: filled = .False.
    end type vtk_field_t
    
    ! Type for storing mesh data
    type vtk_mesh_t
    private
        character(len=:), allocatable :: directory_path        ! Directory where the results are going to be stored
        character(len=:), allocatable :: name_prefix           ! Name prefix of the VTK files
        real(rp)   , allocatable      :: X(:),Y(:),Z(:)        ! Coordinates of the mesh
        integer(ip), allocatable      :: connectivities(:)     ! Connectivity matrix
        integer(ip), allocatable      :: offset(:)             ! VTK element offset
        integer(1) , allocatable      :: cell_types(:)         ! VTK element type
        integer(ip)                   :: number_of_nodes       ! Number of nodes
        integer(ip)                   :: dimensions            ! Dimensions of the mesh
        type(vtk_field_t), allocatable:: unknowns(:)           ! Array storing field_ts info
        type(vtk_field_t), allocatable:: postprocess_fields(:) ! Array storing postprocess field_ts info
        logical                       :: linear_order = .False.
        logical                       :: filled = .False.
        integer(ip)                   :: status = VTK_STATE_UNKNOWN      ! Status of the write process
    end type vtk_mesh_t

    ! Type for storing several mesh data with its field descriptors
    ! It also contains information about the number of parts (PVTK) and time steps (PVD)
    ! It stores the directory path and the prefix where to write in disk
    type vtk_t
        type(vtk_mesh_t), allocatable        :: mesh(:)            ! VTK mesh data and field_t descriptors
         type(serial_fe_space_t),    pointer :: fe_space => NULL() ! Poins to fe_space_t
         class(environment_t),       pointer :: env => NULL()      ! Poins to fe_space_t
         real(rp),         allocatable       :: steps(:)           ! Array of parameters (time, eigenvalues,etc.)
         integer(ip)                         :: steps_counter=0    ! time steps counter
         integer(ip)                         :: num_meshes = 0     ! Number of VTK meshes stored
         integer(ip)                         :: num_steps = 0      ! Number of time steps
         integer(ip)                         :: num_parts = 0      ! Number of parts
         integer(ip)                         :: root_proc = 0      ! Root processor
    contains
        procedure          :: initialize               => vtk_initialize
        procedure          :: free                     => VTK_free
        procedure          :: begin_write              => vtk_begin_write
        procedure          :: end_write                => vtk_end_write
        procedure          :: begin_read               => vtk_begin_read
        procedure          :: end_read                 => vtk_end_read
        procedure, private :: fill_mesh_linear_order   => vtk_fill_mesh_linear_order
        procedure, private :: get_VTK_time_output_path => vtk_get_VTK_time_output_path
        procedure, private :: get_PVD_time_output_path => vtk_get_PVD_time_output_path
        procedure, private :: get_VTK_filename         => vtk_get_VTK_filename
        procedure, private :: get_PVTK_filename        => vtk_get_PVTK_filename
        procedure, private :: set_root_proc            => vtk_set_root_proc
        procedure, private :: set_num_steps            => vtk_set_num_steps
        procedure, private :: set_num_parts            => vtk_set_num_parts
        procedure, private :: set_path                 => vtk_set_path
        procedure, private :: set_prefix               => vtk_set_prefix
        procedure, private :: append_step              => vtk_append_step
        procedure, private :: create_directory         => vtk_create_dir_hierarchy_on_root_process
    end type vtk_t

    character(len=5) :: time_prefix = 'time_'
    character(len=4) :: vtk_ext     = '.vtu'
    character(len=4) :: pvd_ext     = '.pvd'
    character(len=5) :: pvtk_ext    = '.pvtu'

public :: vtk_t

contains


    subroutine vtk_set_path(this, path, mesh_number)
    !-----------------------------------------------------------------
    !< Set the name of the output directory
    !-----------------------------------------------------------------
        class(vtk_t),          intent(INOUT) :: this
        character(len=*),      intent(IN)    :: path
        integer(ip), optional, intent(IN)    :: mesh_number
        integer(ip)                          :: nm = 1
    !-----------------------------------------------------------------
        if(present(mesh_number)) nm = mesh_number
        this%mesh(nm)%directory_path = path
    end subroutine vtk_set_path


    subroutine vtk_set_prefix(this, prefix, mesh_number)
    !-----------------------------------------------------------------
    !< Set the name of the output directory
    !-----------------------------------------------------------------
        class(vtk_t),          intent(INOUT) :: this
        character(len=*),      intent(IN)    :: prefix
        integer(ip), optional, intent(IN)    :: mesh_number
        integer(ip)                          :: nm = 1
    !-----------------------------------------------------------------
        if(present(mesh_number)) nm = mesh_number
        this%mesh(nm)%name_prefix = prefix
    end subroutine vtk_set_prefix


    function vtk_create_dir_hierarchy_on_root_process(this, path, issue_final_barrier) result(res)
    !-----------------------------------------------------------------
    !< The root process create a hierarchy of directories
    !-----------------------------------------------------------------
        class(vtk_t),      intent(INOUT) :: this
        character(len=*),  intent(IN)    :: path
        logical, optional, intent(IN)    :: issue_final_barrier
        logical                          :: ft, ifb
        integer(kind=c_int)              :: res
        integer(ip)                      :: me, np
    !-----------------------------------------------------------------
        me = 0; np = 1; ft = .False.; ifb = .False.
        if(present(issue_final_barrier)) ifb = issue_final_barrier
        check(associated(this%env))
        call this%env%info(me,np) 
        check(this%root_proc <= np-1)

        res=0
        if(me == this%root_proc) then
            res = mkdir_recursive(path//C_NULL_CHAR)
            check ( res == 0 ) 
        end if

        if(ifb) call this%env%first_level_barrier()
    end function vtk_create_dir_hierarchy_on_root_process


    function vtk_get_VTK_time_output_path(this, path, time_step, mesh_number) result(dp)
    !-----------------------------------------------------------------
    !< Build time output dir path for the vtk files in each timestep
    !-----------------------------------------------------------------
        implicit none
        class(vtk_t),               intent(INOUT) :: this
        character(len=*), optional, intent(IN)    :: path
        real(rp),         optional, intent(IN)    :: time_step
        integer(ip),      optional, intent(IN)    :: mesh_number
        character(len=:), allocatable             :: dp
        character(len=:), allocatable             :: fp
        integer(ip)                               :: nm
        real(rp)                                  :: ts
    !-----------------------------------------------------------------¡
        nm = this%num_meshes
        if(present(mesh_number)) nm = mesh_number

        ts = this%num_steps
        if(present(time_step)) ts = time_step 

        fp = this%mesh(nm)%directory_path
        if(present(path)) fp = path

        dp = trim(adjustl(fp))//'/'//time_prefix//trim(adjustl(str(no_sign=.true., n=ts)))//'/'
    end function vtk_get_VTK_time_output_path


    function vtk_get_PVD_time_output_path(this, path, time_step) result(dp)
    !-----------------------------------------------------------------
    !< Build output dir path for the PVD files
    !-----------------------------------------------------------------
        class(vtk_t),               intent(INOUT) :: this
        character(len=*), optional, intent(IN)    :: path
        real(RP),         optional, intent(IN)    :: time_step
        character(len=:), allocatable             :: dp
        character(len=:), allocatable             :: fp
        real(rp)                                  :: ts
    !-----------------------------------------------------------------
        ts = this%num_steps
        if(present(time_step)) ts = time_step 
        dp = time_prefix//trim(adjustl(str(no_sign=.true., n=ts)))//'/'
    end function vtk_get_PVD_time_output_path


    function vtk_get_VTK_filename(this, prefix, part_number, mesh_number) result(fn)
    !-----------------------------------------------------------------
    !< Build VTK filename
    !-----------------------------------------------------------------
        class(vtk_t),               intent(INOUT) :: this
        character(len=*), optional, intent(IN)    :: prefix
        integer(ip),      optional, intent(IN)    :: part_number
        integer(ip),      optional, intent(IN)    :: mesh_number
        character(len=:), allocatable             :: fn
        character(len=:), allocatable             :: fp
        integer(ip)                               :: nm, me, np
    !-----------------------------------------------------------------
        nm = this%num_meshes
        if(present(mesh_number)) nm = mesh_number

        me = 0; np = 1
        call this%env%info(me, np)
        np = me
        if(present(part_number)) np = part_number

        fp = this%mesh(nm)%name_prefix
        if(present(prefix)) fp = prefix

        fn = trim(adjustl(fp))//'_'//trim(adjustl(str(no_sign=.true., n=nm)))//'_'//trim(adjustl(str(no_sign=.true., n=np)))//vtk_ext
    end function vtk_get_VTK_filename


      ! Build VTK filename
    function vtk_get_PVTK_filename(this, prefix, mesh_number, time_step) result(fn)
    !-----------------------------------------------------------------
    !< Build PVTK filename
    !-----------------------------------------------------------------
        class(vtk_t),             intent(INOUT) :: this
        character(len=*), optional, intent(IN)    :: prefix
        integer(ip),      optional, intent(IN)    :: mesh_number
        real(rp),         optional, intent(IN)    :: time_step
        character(len=:), allocatable             :: fn
        character(len=:), allocatable             :: fp
        integer(ip)                               :: nm
        real(rp)                                  :: ts
    !-----------------------------------------------------------------
        nm = this%num_meshes
        if(present(mesh_number)) nm = mesh_number

        ts = 0._rp

        if(allocated(this%steps)) then
            if(this%steps_counter >0 .and. this%steps_counter <= size(this%steps,1)) &
                ts = this%steps(this%steps_counter)
        endif
        if(present(time_step)) ts = time_step

        fp = this%mesh(nm)%name_prefix
        if(present(prefix)) fp = prefix

        fn = trim(adjustl(fp))//'_'//trim(adjustl(str(no_sign=.true., n=nm)))//'_'//trim(adjustl(str(no_sign=.true., n=ts)))//pvtk_ext
    end function vtk_get_PVTK_filename


    subroutine vtk_set_root_proc(this, root)
    !-----------------------------------------------------------------
    !< Set the root processor
    !-----------------------------------------------------------------
        class(vtk_t), intent(INOUT) :: this
        integer(ip),  intent(IN)    :: root
    !-----------------------------------------------------------------
        this%root_proc = root
    end subroutine vtk_set_root_proc


    subroutine vtk_set_num_steps(this, steps)
    !-----------------------------------------------------------------
    !< Set the number of time steps of the simulation to be writen in the PVD
    !-----------------------------------------------------------------
        class(vtk_t), intent(INOUT) :: this
        integer(ip),  intent(IN)    :: steps
    !-----------------------------------------------------------------
        this%num_steps = steps
        if(.not.allocated(this%steps)) then
            call memalloc ( this%num_steps, this%steps, __FILE__,__LINE__)
        elseif(size(this%steps)<this%num_steps) then
            call memrealloc ( steps, this%steps, __FILE__,__LINE__)
        endif
    end subroutine vtk_set_num_steps


    subroutine vtk_append_step(this, current_step)
    !-----------------------------------------------------------------
    !< Append a new time stepstep
    !-----------------------------------------------------------------
        class(vtk_t), intent(INOUT) :: this
        real(rp),     intent(IN)    :: current_step !Current time step
    !-----------------------------------------------------------------
        this%steps_counter = this%steps_counter + 1
        this%num_steps = max(this%num_steps, this%steps_counter)
        if(.not.allocated(this%steps)) then
            call memalloc ( this%num_steps, this%steps, __FILE__,__LINE__)
        elseif(size(this%steps)<this%num_steps) then
            call memrealloc (this%num_steps, this%steps, __FILE__,__LINE__)
        endif        
        this%steps(this%steps_counter) = current_step
    end subroutine vtk_append_step


    subroutine vtk_set_num_parts(this, number_of_parts)
    !-----------------------------------------------------------------
    !< Set the number of parts of the partitioned mesh to be writen in the PVTK
    !-----------------------------------------------------------------
        class(vtk_t), intent(INOUT) :: this
        integer(ip),  intent(IN)    :: number_of_parts
    !-----------------------------------------------------------------
        this%num_parts = number_of_parts
    end subroutine vtk_set_num_parts


    subroutine vtk_initialize(this, triangulation, fe_space, env, path, prefix, root_proc, number_of_parts, number_of_steps, mesh_number, linear_order)
    !-----------------------------------------------------------------
    !< Initialize the vtk_t derived type
    !-----------------------------------------------------------------
        class(vtk_t),                    intent(INOUT) :: this
        type(triangulation_t),           intent(IN)    :: triangulation
        type(serial_fe_space_t), target, intent(INOUT) :: fe_space
        class(environment_t),    target, intent(IN)    :: env
        character(len=*),                intent(IN)    :: path
        character(len=*),                intent(IN)    :: prefix  
        integer(ip),      optional,      intent(IN)    :: root_proc
        integer(ip),      optional,      intent(IN)    :: number_of_parts
        integer(ip),      optional,      intent(IN)    :: number_of_steps
        logical,          optional,      intent(IN)    :: linear_order
        integer(ip),      optional,      intent(OUT)   :: mesh_number
        integer(ip)                                    :: nm = 1
        logical                                        :: lo = .False., ft = .False.
        integer(ip)                                    :: me, np, st, rp
      ! ----------------------------------------------------------------------------------
        if(present(linear_order)) lo = linear_order

        this%fe_space => fe_space
        this%env => env

        me = 0; np = 1; rp = 0
        if(associated(this%env)) then 
            call this%env%info(me,np) 
            ft =  this%env%am_i_fine_task() 
        endif
        if(present(root_proc)) rp = root_proc
        call this%set_root_proc(rp)

        if(ft) then
            if(lo) then 
                call this%fill_mesh_linear_order(triangulation=triangulation, mesh_number=nm)
            else
                check(.false.)
!                call this%initialize_superlinear_order(fe_space=fe_space, mesh_number=nm)
            endif
            if(present(mesh_number)) mesh_number = nm

            call this%set_path(path,nm)
            call this%set_prefix(prefix,nm)
            if(present(number_of_parts)) np = number_of_parts
            call this%set_num_parts(np)
            st = 1
            if(present(number_of_steps)) st = number_of_steps
            call this%set_num_steps(st)

        endif

    end subroutine vtk_initialize


    subroutine initialize_linear_order(this, triangulation, mesh_number)
    !-----------------------------------------------------------------
    !< vtk_t derived type linear order initialization
    !-----------------------------------------------------------------
        class(vtk_t),                    intent(INOUT) :: this
        type(triangulation_t),           intent(IN)    :: triangulation
        integer(ip), optional,           intent(OUT)   :: mesh_number
        integer(ip)                                    :: nm
    !-----------------------------------------------------------------
        call this%fill_mesh_linear_order(triangulation, nm)
        if(present(mesh_number)) mesh_number = nm
    end subroutine initialize_linear_order


    subroutine vtk_fill_mesh_linear_order(this, triangulation, mesh_number)
    !-----------------------------------------------------------------
    !< Store a linear_order mesh in a vtk_t derived type from a triangulation
    !-----------------------------------------------------------------
        class(vtk_t),          intent(INOUT) :: this
        type(triangulation_t), intent(IN)    :: triangulation
        integer(ip), optional, intent(OUT)   :: mesh_number
        type(vtk_mesh_t), allocatable        :: f_vtk_tmp(:)
        integer(ip)                          :: i
        integer(ip)                          :: j
        integer(ip)                          :: tnnod
        integer(ip)                          :: counter
        character(:), pointer                :: topology
    !-----------------------------------------------------------------
        ! Meshes allocation
        if(this%num_meshes == 0) then 
            this%num_meshes = 1
            if(allocated(this%mesh)) deallocate(this%mesh)
            allocate(this%mesh(this%num_meshes))
        else
            this%num_meshes = this%num_meshes + 1 
            call move_alloc(from=this%mesh, to=f_vtk_tmp)
            allocate(this%mesh(this%num_meshes))
            this%mesh(1:size(f_vtk_tmp,dim=1)) = f_vtk_tmp(:)
            deallocate(f_vtk_tmp)
        endif
        if(present(mesh_number)) mesh_number = this%num_meshes
        this%mesh(this%num_meshes)%linear_order = .True.

        call memalloc ( triangulation%num_elems, this%mesh(this%num_meshes)%cell_types, __FILE__,__LINE__)
        call memalloc ( triangulation%num_elems, this%mesh(this%num_meshes)%offset, __FILE__,__LINE__)
        this%mesh(this%num_meshes)%cell_types = 0; this%mesh(this%num_meshes)%offset = 0

        ! Fill VTK cell type and offset arrays and and count nodes
        tnnod = 0
        do i=1, triangulation%num_elems
            topology => triangulation%elems(i)%reference_fe_geo%get_topology()
            if(topology == topology_quad) then
                this%mesh(this%num_meshes)%cell_types(i) = 8 ! VTK_VOXEL
            else
                write(*,*) 'fill_mesh_from_triangulation: Topology not supported -> ', topology
                check(.false.)
            endif
            tnnod = tnnod + triangulation%elems(i)%reference_fe_geo%get_number_vertices()
            this%mesh(this%num_meshes)%offset(i:) = tnnod
        enddo
        this%mesh(this%num_meshes)%number_of_nodes = tnnod
        this%mesh(this%num_meshes)%dimensions = triangulation%num_dims

        call memalloc ( tnnod, this%mesh(this%num_meshes)%connectivities, __FILE__,__LINE__)
        call memalloc ( tnnod, this%mesh(this%num_meshes)%X, __FILE__,__LINE__)
        call memalloc ( tnnod, this%mesh(this%num_meshes)%Y, __FILE__,__LINE__)
        call memalloc ( tnnod, this%mesh(this%num_meshes)%Z, __FILE__,__LINE__)
        this%mesh(this%num_meshes)%X = 0._rp; this%mesh(this%num_meshes)%Y = 0._rp; this%mesh(this%num_meshes)%Z = 0._rp
        counter = 1
        tnnod = 0

        ! Fill VTK coordinate arrays
        do i=1, triangulation%num_elems
            do j=1, triangulation%elems(i)%reference_fe_geo%get_number_vertices()
                this%mesh(this%num_meshes)%connectivities(tnnod+j) = j + tnnod - 1
                if (triangulation%num_dims >=1) this%mesh(this%num_meshes)%X(counter) = triangulation%elems(i)%coordinates(1,j)
                if (triangulation%num_dims >=2) this%mesh(this%num_meshes)%Y(counter) = triangulation%elems(i)%coordinates(2,j)
                if (triangulation%num_dims >=3) this%mesh(this%num_meshes)%Z(counter) = triangulation%elems(i)%coordinates(3,j)
                counter = counter + 1
            enddo
            tnnod = tnnod + triangulation%elems(i)%reference_fe_geo%get_number_vertices()
        enddo
        this%mesh(this%num_meshes)%filled = .True.
    end subroutine vtk_fill_mesh_linear_order


    function vtk_begin_write(this, file_name, part_number, time_step, mesh_number, format, f_id) result(E_IO)
    !-----------------------------------------------------------------
    !< Start the writing of a single VTK file to disk (if I am fine MPI task)
    !< Writes connectivities and coordinates ( VTK_INI_XML, 
    !< VTK_GEO_XML, VTK_CON_XML )
    !-----------------------------------------------------------------
        class(vtk_t),               intent(INOUT) :: this        !< VTK_t derived type
        character(len=*), optional, intent(IN)    :: file_name   !< VTK File NAME
        integer(ip),      optional, intent(IN)    :: part_number !< Number of the PART
        real(rp),         optional, intent(IN)    :: time_step   !< Time STEP value
        integer(ip),      optional, intent(IN)    :: mesh_number !< Number of the MESH
        character(len=*), optional, intent(IN)    :: format      !< Ouput ForMaT
        integer(ip),      optional, intent(OUT)   :: f_id        !< File ID
        character(len=:), allocatable             :: fn          !< Real File Name
        character(len=:), allocatable             :: dp          !< Real Directory Path
        character(len=:), allocatable             :: of          !< Real Output Format
        real(rp)                                  :: ts          !< Real Time Step
        logical                                   :: ft          !< Fine Task
        integer(ip)                               :: nm          !< Real Number of the Mesh
        integer(ip)                               :: np          !< Real Number of the Part
        integer(ip)                               :: me          !< Task identifier
        integer(ip)                               :: fid         !< Real File ID
        integer(ip)                               :: nnods       !< Number of NODeS
        integer(ip)                               :: nels        !< Number of ELementS
        integer(ip)                               :: E_IO        !< Error IO
      ! ----------------------------------------------------------------------------------
        check(associated(this%env))
     
        ft =  this%env%am_i_fine_task() 

        E_IO = 0
        fid = -1
        
        if(ft) then
            me = 0; np = 1
            call this%env%info(me,np) 
            np = me
            if(present(part_number)) np = part_number

            nm = this%num_meshes
            if(present(mesh_number)) nm = mesh_number
        
            ts = 0._rp
            if(present(time_step)) ts = time_step 
            call this%append_step(ts)


            if(this%mesh(nm)%status == VTK_STATE_UNKNOWN .or. this%mesh(nm)%status == VTK_STATE_ENDED) then
        
                dp = this%get_VTK_time_output_path(path=this%mesh(nm)%directory_path, time_step=ts, mesh_number=nm)
                fn = this%get_VTK_filename(prefix=this%mesh(nm)%name_prefix, part_number=np, mesh_number=nm)
                fn = dp//fn
                if(present(file_name)) fn = file_name

                if( this%create_directory(dp,issue_final_barrier=.True.) == 0) then    
                    of = 'raw'
                    if(present(format)) of = trim(adjustl(format))
            
                    nnods = this%mesh(nm)%number_of_nodes
                    nels = size(this%mesh(nm)%cell_types, dim=1)
            
                    E_IO = VTK_INI_XML(output_format = trim(adjustl(of)),  &
                                       filename = trim(adjustl(fn)),       &
                                       mesh_topology = 'UnstructuredGrid', &
                                       cf=fid)
                    E_IO = VTK_GEO_XML(NN = nnods,           &
                                       NC = nels,            &
                                       X  = this%mesh(nm)%X, &
                                       Y  = this%mesh(nm)%Y, &
                                       Z  = this%mesh(nm)%Z, &
                                       cf = fid)
                    E_IO = VTK_CON_XML(NC = nels,                                &
                                       connect   = this%mesh(nm)%connectivities, &
                                       offset    = this%mesh(nm)%offset,         &
                                       cell_type = this%mesh(nm)%cell_types,     &
                                       cf        = fid)
                    this%mesh(nm)%status = VTK_STATE_WRITE_STARTED
                endif
            endif
            if(present(f_id)) f_id = fid

        endif
    end function vtk_begin_write


    function vtk_end_write(this, mesh_number, f_id) result(E_IO)
    !-----------------------------------------------------------------
    !< Ends the writing of a single VTK file to disk (if I am fine MPI task)
    !< Closes geometry ( VTK_END_XML, VTK_GEO_XML )
    !-----------------------------------------------------------------
        implicit none
        class(vtk_t),               intent(INOUT) :: this        !< VTK_t derived type
        integer(ip),      optional, intent(IN)    :: mesh_number !< Number of MESH
        integer(ip),      optional, intent(INOUT) :: f_id        !< File ID
        integer(ip)                               :: nm          !< Real Number of Mesh
        integer(ip)                               :: E_IO        !< IO Error
        logical                                   :: ft          !< Fine Task
      ! ----------------------------------------------------------------------------------
        check(associated(this%env))
        ft =  this%env%am_i_fine_task() 

        E_IO = 0

        if (ft) then
           
           nm = this%num_meshes
           if(present(mesh_number)) nm = mesh_number
           
           if ((this%mesh(nm)%status >= VTK_STATE_WRITE_STARTED) .and. (this%mesh(nm)%status /= VTK_STATE_ENDED)) then
              
              if(this%mesh(nm)%status == VTK_STATE_POINTDATA_OPENED) then
                 if(present(f_id)) then
                    E_IO = VTK_DAT_XML(var_location='node',var_block_action='close', cf=f_id)
                 else
                    E_IO = VTK_DAT_XML(var_location='node',var_block_action='close')
                 endif
                 this%mesh(nm)%status = VTK_STATE_POINTDATA_CLOSED
              endif
              
              if(present(f_id)) then       
                 E_IO = VTK_GEO_XML(cf=f_id)
                 E_IO = VTK_END_XML(cf=f_id)
              else
                 E_IO = VTK_GEO_XML()
                 E_IO = VTK_END_XML()
              endif
              
              this%mesh(nm)%status = VTK_STATE_ENDED
           endif
        endif
    end function vtk_end_write


    function vtk_begin_read(this, filename, part_number, time_step, mesh_number, format, f_id) result(E_IO)
    !-----------------------------------------------------------------
    !< Begin the reading of a single VTK file to disk (if I am fine MPI task)
    !-----------------------------------------------------------------
        class(vtk_t),             intent(INOUT)   :: this        !< VTK_t derived type
        character(len=*), optional, intent(IN)    :: filename    !< VTK File NAME
        integer(ip),      optional, intent(IN)    :: part_number !< Number of the PART
        real(rp),         optional, intent(IN)    :: time_step   !< Time STEP value
        integer(ip),      optional, intent(IN)    :: mesh_number !< Number of the MESH
        character(len=*), optional, intent(IN)    :: format      !< Ouput ForMaT
        integer(ip),      optional, intent(OUT)   :: f_id        !< File ID
        character(len=:), allocatable             :: fn          !< Real File Name
        character(len=:), allocatable             :: dp          !< Real Directory Path
        character(len=:), allocatable             :: of          !< Real Output Format
        real(rp)                                  :: ts          !< Real Time Step
        logical                                   :: ft          !< Fine Task
        integer(ip)                               :: nm          !< Real Number of the Mesh
        integer(ip)                               :: np          !< Real Number of the Part
        integer(ip)                               :: me          !< Task identifier
        integer(ip)                               :: fid         !< Real File ID
        integer(ip)                               :: nnods       !< Number of NODeS
        integer(ip)                               :: nels        !< Number of ELementS
        integer(ip)                               :: E_IO        !< Error IO
    !-----------------------------------------------------------------
        check(associated(this%env))
     
        ft =  this%env%am_i_fine_task() 

        E_IO = 0
        fid = -1
        
        if(ft) then
            me = 0; np = 1
            call this%env%info(me,np) 
            np = me
            if(present(part_number)) np = part_number

            nm = this%num_meshes
            if(present(mesh_number)) nm = mesh_number
        
            ts = 0._rp
            if(present(time_step)) ts = time_step
            call this%append_step(ts)

            if(this%mesh(nm)%status == VTK_STATE_UNKNOWN .or. this%mesh(nm)%status == VTK_STATE_ENDED) then
        
                dp = this%get_VTK_time_output_path(path=this%mesh(nm)%directory_path, time_step=ts, mesh_number=nm)
                fn = this%get_VTK_filename(prefix=this%mesh(nm)%name_prefix, part_number=np, mesh_number=nm)
                fn = dp//fn
                if(present(filename)) fn = filename

                if( this%create_directory(dp,issue_final_barrier=.True.) == 0) then    
                    of = 'raw'
                    if(present(format)) of = trim(adjustl(format))
            
                    nnods = this%mesh(nm)%number_of_nodes
                    nels = size(this%mesh(nm)%cell_types, dim=1)
            
                    E_IO = VTK_INI_XML_READ(input_format = trim(adjustl(of)), filename = trim(adjustl(fn)), mesh_topology = 'UnstructuredGrid', cf=fid)
                    this%mesh(nm)%status = VTK_STATE_READ_STARTED
                endif
            endif
            if(present(f_id)) f_id = fid

        endif
    end function  vtk_begin_read


    function vtk_end_read(this, mesh_number, f_id) result(E_IO)
    !-----------------------------------------------------------------
    !< End the reading of a single VTK file to disk (if I am fine MPI task)
    !-----------------------------------------------------------------
        class(vtk_t),               intent(INOUT) :: this        !< VTK_t derived type
        integer(ip),      optional, intent(IN)    :: mesh_number !< Number of MESH
        integer(ip),      optional, intent(INOUT) :: f_id        !< File ID
        integer(ip)                               :: nm          !< Real Number of Mesh
        integer(ip)                               :: E_IO        !< IO Error
        logical                                   :: ft          !< Fine Task
    !-----------------------------------------------------------------
        check(associated(this%env))
        ft =  this%env%am_i_fine_task() 

        E_IO = 0

        nm = this%num_meshes
        if(present(mesh_number)) nm = mesh_number
        
        if(ft .and. (this%mesh(nm)%status >= VTK_STATE_READ_STARTED) .and. (this%mesh(nm)%status /= VTK_STATE_ENDED)) then

            if(present(f_id)) then       
                E_IO = VTK_END_XML_READ(cf=f_id)
            else
                E_IO = VTK_END_XML_READ()
            endif

            this%mesh(nm)%status = VTK_STATE_ENDED
        endif
    end function vtk_end_read


    subroutine VTK_free (this)
    !-----------------------------------------------------------------
    !< Free the vtk_t derived type
    !-----------------------------------------------------------------
        class(vtk_t), intent(inout) :: this
        integer(ip)                 :: i
        integer(ip)                 :: j
        logical                     :: ft
    !-----------------------------------------------------------------
        check(associated(this%env))
        ft = this%env%am_i_fine_task() 
    
        if(ft) then
    
            if(allocated(this%mesh)) then
                do i=1, size(this%mesh,1)
                    call memfree(this%mesh(i)%X, __FILE__,__LINE__)
                    call memfree(this%mesh(i)%Y, __FILE__,__LINE__)
                    call memfree(this%mesh(i)%Z, __FILE__,__LINE__)
                    call memfree(this%mesh(i)%connectivities, __FILE__,__LINE__)
                    call memfree(this%mesh(i)%offset, __FILE__,__LINE__)
                    call memfree(this%mesh(i)%cell_types, __FILE__,__LINE__)
                    if(allocated(this%mesh(i)%unknowns)) then
                        do j=1, size(this%mesh(i)%unknowns)
                            if(allocated(this%mesh(i)%unknowns(j)%var_location)) deallocate(this%mesh(i)%unknowns(j)%var_location)
                            if(allocated(this%mesh(i)%unknowns(j)%var_name)) deallocate(this%mesh(i)%unknowns(j)%var_name)
                            if(allocated(this%mesh(i)%unknowns(j)%field_type)) deallocate(this%mesh(i)%unknowns(j)%field_type)
                            this%mesh(i)%unknowns(j)%filled = .False.
                        enddo
                        deallocate(this%mesh(i)%unknowns)
                    endif
                    this%mesh(i)%filled = .False.
                enddo
                deallocate(this%mesh)
            endif

            if (allocated(this%steps)) call memfree(this%steps, __FILE__,__LINE__)

        endif
        this%num_meshes = 0
        this%num_steps = 0
        this%num_parts = 0
        this%root_proc = 0
        this%fe_space => NULL()
        this%env => NULL()
    end subroutine VTK_free


end module lib_vtk_io_interface_names
