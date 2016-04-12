
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
USE types_names
USE memor_names
USE ir_precision,                only: str
USE vector_names,                only: vector_t
USE environment_names,           only: environment_t
USE triangulation_names,         only: triangulation_t
USE serial_fe_space_names,       only: serial_fe_space_t, finite_element_t, fe_function_t
USE serial_scalar_array_names,   only: serial_scalar_array_t
USE reference_fe_names,          only: reference_fe_t, quad_lagrangian_reference_fe_t, fe_map_t, quadrature_t, interpolation_t
USE field_names,                 only: point_t
USE allocatable_array_ip2_names, only: allocatable_array_ip2_t
USE reference_fe_names,          only: topology_quad, topology_tet, fe_type_lagrangian
USE vtk_mesh


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

    ! Type for storing several mesh data with its field descriptors
    ! It also contains information about the number of parts (PVTK) and time steps (PVD)
    ! It stores the directory path and the prefix where to write in disk
    type vtk_handler_t
    private
        type(serial_fe_space_t),    pointer :: fe_space => NULL() ! Poins to fe_space_t
        class(environment_t),       pointer :: env      => NULL() ! Poins to fe_space_t
        type(vtk_mesh_t), allocatable       :: mesh(:)            ! VTK mesh data and field_t descriptors
        real(rp),         allocatable       :: steps(:)           ! Array of parameters (time, eigenvalues,etc.)
        integer(ip)                         :: steps_counter = 0  ! time steps counter
        integer(ip)                         :: num_meshes    = 0  ! Number of VTK meshes stored
        integer(ip)                         :: num_steps     = 0  ! Number of time steps
        integer(ip)                         :: num_parts     = 0  ! Number of parts
        integer(ip)                         :: root_proc     = 0  ! Root processor
    contains
    private
        procedure, public :: initialize                   => vtk_initialize
        procedure, public :: begin_write                  => vtk_begin_write
        procedure, public :: write_node_field             => vtk_write_node_field
        procedure, public :: end_write                    => vtk_end_write
        procedure, public :: free                         => vtk_free
        procedure         :: reallocate_meshes            => vtk_reallocate_meshes
        procedure         :: fill_mesh_linear_order       => vtk_fill_mesh_linear_order
        procedure         :: fill_mesh_superlinear_order  => vtk_fill_mesh_superlinear_order
        procedure         :: write_node_field_superlinear => vtk_write_node_field_superlinear
        procedure         :: get_VTK_time_output_path     => vtk_get_VTK_time_output_path
        procedure         :: get_PVD_time_output_path     => vtk_get_PVD_time_output_path
        procedure         :: get_VTK_filename             => vtk_get_VTK_filename
        procedure         :: get_PVTK_filename            => vtk_get_PVTK_filename
        procedure         :: set_root_proc                => vtk_set_root_proc
        procedure         :: set_num_steps                => vtk_set_num_steps
        procedure         :: set_num_parts                => vtk_set_num_parts
        procedure         :: set_path                     => vtk_set_path
        procedure         :: set_prefix                   => vtk_set_prefix
        procedure         :: append_step                  => vtk_append_step
        procedure         :: create_directory             => vtk_create_dir_hierarchy_on_root_process
    end type vtk_handler_t

    character(len=5) :: time_prefix = 'time_'
    character(len=4) :: vtk_ext     = '.vtu'
    character(len=4) :: pvd_ext     = '.pvd'
    character(len=5) :: pvtk_ext    = '.pvtu'

public :: vtk_handler_t

contains


    subroutine vtk_set_path(this, path, mesh_number)
    !-----------------------------------------------------------------
    !< Set the name of the output directory
    !-----------------------------------------------------------------
        class(vtk_handler_t),  intent(INOUT) :: this
        character(len=*),      intent(IN)    :: path
        integer(ip), optional, intent(IN)    :: mesh_number
        integer(ip)                          :: nm = 1
    !-----------------------------------------------------------------
        if(present(mesh_number)) nm = mesh_number
        call this%mesh(nm)%set_path(path=path)
    end subroutine vtk_set_path


    subroutine vtk_set_prefix(this, prefix, mesh_number)
    !-----------------------------------------------------------------
    !< Set the name of the output directory
    !-----------------------------------------------------------------
        class(vtk_handler_t),  intent(INOUT) :: this
        character(len=*),      intent(IN)    :: prefix
        integer(ip), optional, intent(IN)    :: mesh_number
        integer(ip)                          :: nm = 1
    !-----------------------------------------------------------------
        if(present(mesh_number)) nm = mesh_number
        call this%mesh(nm)%set_prefix(prefix=prefix)
    end subroutine vtk_set_prefix


    function vtk_create_dir_hierarchy_on_root_process(this, path, issue_final_barrier) result(res)
    !-----------------------------------------------------------------
    !< The root process create a hierarchy of directories
    !-----------------------------------------------------------------
        class(vtk_handler_t), intent(INOUT) :: this
        character(len=*),     intent(IN)    :: path
        logical, optional,    intent(IN)    :: issue_final_barrier
        logical                             :: ft, ifb
        integer(kind=c_int)                 :: res
        integer(ip)                         :: me, np
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
        class(vtk_handler_t),       intent(INOUT) :: this
        character(len=*), optional, intent(IN)    :: path
        real(rp),         optional, intent(IN)    :: time_step
        integer(ip),      optional, intent(IN)    :: mesh_number
        character(len=:), allocatable             :: dp
        character(len=:), allocatable             :: fp
        integer(ip)                               :: nm
        real(rp)                                  :: ts
    !-----------------------------------------------------------------
        nm = this%num_meshes
        if(present(mesh_number)) nm = mesh_number

        ts = this%num_steps
        if(present(time_step)) ts = time_step 

        call this%mesh(nm)%get_path(fp)
        if(present(path)) fp = path

        dp = trim(adjustl(fp))//'/'//time_prefix//trim(adjustl(str(no_sign=.true., n=ts)))//'/'
    end function vtk_get_VTK_time_output_path


    function vtk_get_PVD_time_output_path(this, path, time_step) result(dp)
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
    end function vtk_get_PVD_time_output_path


    function vtk_get_VTK_filename(this, prefix, part_number, mesh_number) result(fn)
    !-----------------------------------------------------------------
    !< Build VTK filename
    !-----------------------------------------------------------------
        class(vtk_handler_t),       intent(INOUT) :: this
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

        call this%mesh(nm)%get_prefix(fp)
        if(present(prefix)) fp = prefix

        fn = trim(adjustl(fp))//'_'//trim(adjustl(str(no_sign=.true., n=nm)))//'_'//trim(adjustl(str(no_sign=.true., n=np)))//vtk_ext
    end function vtk_get_VTK_filename


      ! Build VTK filename
    function vtk_get_PVTK_filename(this, prefix, mesh_number, time_step) result(fn)
    !-----------------------------------------------------------------
    !< Build PVTK filename
    !-----------------------------------------------------------------
        class(vtk_handler_t),       intent(INOUT) :: this
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

        call this%mesh(nm)%get_prefix(fp)
        if(present(prefix)) fp = prefix

        fn = trim(adjustl(fp))//'_'//trim(adjustl(str(no_sign=.true., n=nm)))//'_'//trim(adjustl(str(no_sign=.true., n=ts)))//pvtk_ext
    end function vtk_get_PVTK_filename


    subroutine vtk_set_root_proc(this, root)
    !-----------------------------------------------------------------
    !< Set the root processor
    !-----------------------------------------------------------------
        class(vtk_handler_t), intent(INOUT) :: this
        integer(ip),          intent(IN)    :: root
    !-----------------------------------------------------------------
        this%root_proc = root
    end subroutine vtk_set_root_proc


    subroutine vtk_set_num_steps(this, steps)
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
    end subroutine vtk_set_num_steps


    subroutine vtk_append_step(this, current_step)
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
            call memrealloc (this%num_steps, this%steps, __FILE__,__LINE__)
        endif        
        this%steps(this%steps_counter) = current_step
    end subroutine vtk_append_step


    subroutine vtk_set_num_parts(this, number_of_parts)
    !-----------------------------------------------------------------
    !< Set the number of parts of the partitioned mesh to be writen in the PVTK
    !-----------------------------------------------------------------
        class(vtk_handler_t), intent(INOUT) :: this
        integer(ip),          intent(IN)    :: number_of_parts
    !-----------------------------------------------------------------
        this%num_parts = number_of_parts
    end subroutine vtk_set_num_parts


    subroutine vtk_reallocate_meshes(this)
    !-----------------------------------------------------------------
    !< Reallocate meshes array to fit a new mesh
    !-----------------------------------------------------------------
        class(vtk_handler_t), intent(INOUT) :: this
        type(vtk_mesh_t), allocatable       :: f_vtk_tmp(:)
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
    end subroutine vtk_reallocate_meshes


    subroutine vtk_initialize(this, triangulation, fe_space, env, path, prefix, root_proc, number_of_parts, number_of_steps, mesh_number, linear_order)
    !-----------------------------------------------------------------
    !< Initialize the vtk_handler_t derived type
    !-----------------------------------------------------------------
        class(vtk_handler_t),            intent(INOUT) :: this
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
    !-----------------------------------------------------------------
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
                call this%fill_mesh_superlinear_order(fe_space=fe_space, nmesh=nm)
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
    !< vtk_handler_t derived type linear order initialization
    !-----------------------------------------------------------------
        class(vtk_handler_t),  intent(INOUT) :: this
        type(triangulation_t), intent(IN)    :: triangulation
        integer(ip), optional, intent(OUT)   :: mesh_number
        integer(ip)                          :: nm
    !-----------------------------------------------------------------
        call this%fill_mesh_linear_order(triangulation, nm)
        if(present(mesh_number)) mesh_number = nm
    end subroutine initialize_linear_order


    subroutine vtk_fill_mesh_linear_order(this, triangulation, mesh_number)
    !-----------------------------------------------------------------
    !< Store a linear_order mesh in a vtk_handler_t derived type from a triangulation
    !-----------------------------------------------------------------
        class(vtk_handler_t),  intent(INOUT) :: this
        type(triangulation_t), intent(IN)    :: triangulation
        integer(ip), optional, intent(OUT)   :: mesh_number
        integer(ip)                          :: i
        integer(ip)                          :: j
        integer(ip)                          :: number_nodes
        integer(ip)                          :: counter
        character(:), pointer                :: topology
    !-----------------------------------------------------------------
        call this%reallocate_meshes()
        if(present(mesh_number)) mesh_number = this%num_meshes
        call this%mesh(this%num_meshes)%set_linear_order(.True.)

        call this%mesh(this%num_meshes)%set_dimensions(triangulation%num_dims)
        call this%mesh(this%num_meshes)%set_number_elements(triangulation%num_elems)
        call this%mesh(this%num_meshes)%allocate_elemental_arrays()

        ! Fill VTK cell type and offset arrays and and count nodes
        number_nodes = 0
        do i=1, this%mesh(this%num_meshes)%get_number_elements()
            topology => triangulation%elems(i)%reference_fe_geo%get_topology()
            if(topology == topology_quad) then
                call this%mesh(this%num_meshes)%set_cell_type(index=i, type=8) ! VTK_VOXEL
            else
                write(*,*) 'fill_mesh_from_triangulation: Topology not supported -> ', topology
                check(.false.)
            endif
            number_nodes = number_nodes + triangulation%elems(i)%reference_fe_geo%get_number_vertices()
            call this%mesh(this%num_meshes)%set_offset(i, number_nodes)
        enddo
        call this%mesh(this%num_meshes)%set_number_nodes(number_nodes)
        call this%mesh(this%num_meshes)%allocate_nodal_arrays()


        counter = 1
        number_nodes = 0

        ! Fill VTK coordinate arrays
        do i=1, this%mesh(this%num_meshes)%get_number_elements()
            do j=1, triangulation%elems(i)%reference_fe_geo%get_number_vertices()
                call this%mesh(this%num_meshes)%set_connectivity(number_nodes+j, j+number_nodes-1)
                if (triangulation%num_dims >=1) call this%mesh(this%num_meshes)%set_x_coordinate(counter, triangulation%elems(i)%coordinates(1,j))
                if (triangulation%num_dims >=2) call this%mesh(this%num_meshes)%set_y_coordinate(counter, triangulation%elems(i)%coordinates(2,j))
                if (triangulation%num_dims >=3) call this%mesh(this%num_meshes)%set_z_coordinate(counter, triangulation%elems(i)%coordinates(3,j))
                counter = counter + 1
            enddo
            number_nodes = number_nodes + triangulation%elems(i)%reference_fe_geo%get_number_vertices()
        enddo
        call this%mesh(this%num_meshes)%set_filled(.true.)
    end subroutine vtk_fill_mesh_linear_order


    subroutine vtk_fill_mesh_superlinear_order(this, fe_space, nmesh)
    !-----------------------------------------------------------------
    !< Store a superlinear_order mesh in a vtk_handler_t derived type from a fe_space
    !-----------------------------------------------------------------
        class(vtk_handler_t),           intent(inout) :: this
        type(serial_fe_space_t),        intent(inout) :: fe_space
        integer(ip), optional,          intent(out)   :: nmesh
        type(point_t),          pointer               :: vertex_coordinates(:) => null()
        type(point_t),          pointer               :: nodal_coordinates(:) => null()
        type(quadrature_t),     pointer               :: nodal_quadrature => null()
        class(reference_fe_t),  pointer               :: reference_fe_geo => null()
        type(finite_element_t), pointer               :: fe => null()
        integer(ip),            pointer               :: subelements_connectivity(:,:) => null()
        type(fe_map_t)                                :: fe_map
        integer(ip)                                   :: num_elements
        integer(ip)                                   :: num_nodes_per_element
        integer(ip)                                   :: num_vertices_per_element
        integer(ip)                                   :: num_subelements_per_element
        integer(ip)                                   :: elements_counter
        integer(ip)                                   :: vertex
        integer(ip)                                   :: max_order
        integer(ip)                                   :: dimensions
        integer(ip)                                   :: subelement_vertex
        integer(ip)                                   :: fe_space_index
        integer(ip)                                   :: element_index
        integer(ip)                                   :: subelement_index
        integer(ip)                                   :: max_order_fe_space_index
        integer(ip)                                   :: nodes_counter
    !-----------------------------------------------------------------
        call this%reallocate_meshes()
        if(present(nmesh)) nmesh = this%num_meshes
        call this%mesh(this%num_meshes)%set_linear_order(.false.)


        ! Get dimensions from the geometric reference finite element
        fe => fe_space%get_finite_element(1)
        reference_fe_geo => fe%get_reference_fe_geo()
        dimensions = reference_fe_geo%get_number_dimensions()
        call this%mesh(this%num_meshes)%set_dimensions(dimensions)

        ! Get max order and its fe_space_index
        max_order = 0
        do fe_space_index=1,fe_space%get_number_fe_spaces()
            if(max_order<fe_space%get_order(fe_space_index)) max_order_fe_space_index = fe_space_index
            max_order = max(max_order,fe_space%get_order(fe_space_index))
        end do

        ! Getters
        nodal_quadrature => fe_space%get_nodal_quadrature(max_order_fe_space_index)
     
        ! Calculate the number of subelems and points for the postprocess
        num_elements = fe_space%get_number_elements()
        !num_nodes_per_element = fe_space%get_number_nodes(max_order_fe_space_index)
        num_vertices_per_element = reference_fe_geo%get_number_vertices()
        num_subelements_per_element = fe_space%get_number_subelements(max_order_fe_space_index)
        call this%mesh(this%num_meshes)%set_number_elements(num_elements * num_subelements_per_element)
        call this%mesh(this%num_meshes)%set_number_nodes(num_vertices_per_element * this%mesh(this%num_meshes)%get_number_elements())

        ! Allocate VTK geometry and connectivity data
        call this%mesh(this%num_meshes)%allocate_nodal_arrays()
        call this%mesh(this%num_meshes)%allocate_elemental_arrays()
        

        ! Get the connectivity of the subelements
        call this%mesh(this%num_meshes)%allocate_subelements_connectivity(num_vertices_per_element, num_subelements_per_element)
        subelements_connectivity => this%mesh(this%num_meshes)%get_subelements_connectivity()
        call fe_space%get_subelements_connectivity(max_order_fe_space_index, subelements_connectivity)

        ! Create FE map
        call fe_map%create(nodal_quadrature, reference_fe_geo)

        ! Translate coordinates and connectivities to VTK format
        nodes_counter = 0
        elements_counter = 0
        do element_index = 1, num_elements
            ! Get Finite element
            fe => fe_space%get_finite_element(element_index)  

            ! Interpolate coordinates
            vertex_coordinates => fe_map%get_coordinates()
            call fe%get_cell_coordinates(vertex_coordinates)
            call fe_map%compute_quadrature_coordinates()
            nodal_coordinates => fe_map%get_quadrature_coordinates()

            ! Fill VTK mesh
            do subelement_index = 1, num_subelements_per_element
                elements_counter = elements_counter + 1
                do vertex = 1, num_vertices_per_element
                    subelement_vertex = subelements_connectivity(vertex, subelement_index)
                    nodes_counter = nodes_counter + 1
                    if (dimensions >= 1) call this%mesh(this%num_meshes)%set_x_coordinate(nodes_counter, nodal_coordinates(subelement_vertex)%get(1))
                    if (dimensions >= 2) call this%mesh(this%num_meshes)%set_y_coordinate(nodes_counter, nodal_coordinates(subelement_vertex)%get(2))
                    if (dimensions >= 3) call this%mesh(this%num_meshes)%set_z_coordinate(nodes_counter, nodal_coordinates(subelement_vertex)%get(3))
                    call this%mesh(this%num_meshes)%set_connectivity(nodes_counter, nodes_counter-1)
                end do

                ! Store the type of element
                call this%mesh(this%num_meshes)%set_cell_type(elements_counter, 8)  ! VTK_VOXEL

                ! Fill offset
                call this%mesh(this%num_meshes)%set_offset(elements_counter, nodes_counter)
            end do
        end do

        call fe_map%free()
    end subroutine vtk_fill_mesh_superlinear_order


    function vtk_begin_write(this, file_name, part_number, time_step, mesh_number, format, f_id) result(E_IO)
    !-----------------------------------------------------------------
    !< Start the writing of a single VTK file to disk (if I am fine MPI task)
    !< Writes connectivities and coordinates ( VTK_INI_XML, 
    !< VTK_GEO_XML, VTK_CON_XML )
    !-----------------------------------------------------------------
        class(vtk_handler_t),       intent(INOUT) :: this        !< vtk_handler_t derived type
        character(len=*), optional, intent(IN)    :: file_name   !< VTK File NAME
        integer(ip),      optional, intent(IN)    :: part_number !< Number of the PART
        real(rp),         optional, intent(IN)    :: time_step   !< Time STEP value
        integer(ip),      optional, intent(IN)    :: mesh_number !< Number of the MESH
        character(len=*), optional, intent(IN)    :: format      !< Ouput ForMaT
        integer(ip),      optional, intent(OUT)   :: f_id        !< File ID
        character(len=:), allocatable             :: path        !< Directory path
        character(len=:), allocatable             :: prefix      !< Filename prefix
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

            call this%mesh(nm)%get_path(path)
            call this%mesh(nm)%get_prefix(prefix)
            dp = this%get_VTK_time_output_path(path=path, time_step=ts, mesh_number=nm)
            fn = this%get_VTK_filename(prefix=prefix, part_number=np, mesh_number=nm)
            fn = dp//fn

            if(present(file_name)) fn = file_name

            if( this%create_directory(dp,issue_final_barrier=.True.) == 0) then    
                E_IO = this%mesh(nm)%begin_write(fn, np, ts, format, f_id)
            endif

            if(present(f_id)) f_id = fid
        endif
    end function vtk_begin_write


    function vtk_write_node_field(this, fe_function, field_name, mesh_number, f_id) result(E_IO)
    !-----------------------------------------------------------------
    !< Write node field to file
    !-----------------------------------------------------------------
        class(vtk_handler_t),       intent(INOUT) :: this           !< vtk_handler_t derived type
        type(fe_function_t),        intent(INOUT) :: fe_function    !< fe_function containing the field to be written
        character(len=*), optional, intent(IN)    :: field_name     !< name of the field
        integer(ip),      optional, intent(IN)    :: mesh_number    !< Real Number of Mesh
        integer(ip),      optional, intent(IN)    :: f_id           !< File ID
        integer(ip)                               :: nm             !< Aux number of mesh
        character(len=:), allocatable             :: fn             !< Aux name of the field
        logical                                   :: ft             !< Fine Task
        integer(ip)                               :: E_IO           !< IO Error
    !-----------------------------------------------------------------
        E_IO = 0
        ft =  this%env%am_i_fine_task() 
        if(ft) then        
           nm = this%num_meshes
           if(present(mesh_number)) nm = mesh_number

           fn = 'unknown'
           if(present(field_name)) fn = field_name

            if(this%mesh(nm)%is_linear_order()) then
!                E_IO = this%write_field_linear(this, fe_function, mesh_number, f_id)
                check(.false.)
            else
                E_IO = this%write_node_field_superlinear(fe_function, fn, nm, f_id)
            endif

        endif
    end function vtk_write_node_field


    function vtk_write_node_field_superlinear(this, fe_function, field_name, mesh_number, f_id) result(E_IO)
    !-----------------------------------------------------------------
    !< Write superlinear field to file
    !-----------------------------------------------------------------
        implicit none
        class(vtk_handler_t),       intent(INOUT) :: this                                         !< vtk_handler_t derived type
        type(fe_function_t),        intent(INOUT) :: fe_function                                  !< Postprocess field structure to be written
        character(len=*), optional, intent(IN)    :: field_name                                   !< name of the field
        integer(ip),                intent(IN)    :: mesh_number                                  !< Real Number of Mesh
        integer(ip),      optional, intent(IN)    :: f_id                                         !< File ID
        real(rp), allocatable                     :: field(:,:)                                   !< FIELD(ncomp,nnod)
        real(rp), allocatable                     :: nodal_values_origin(:)                       !< nodal values of the origin fe_space
        real(rp), allocatable                     :: nodal_values_target(:)                       !< nodal values for the interpolation
        type(quadrature_t),          pointer      :: nodal_quadrature_target            => null() !< Nodal quadrature
        type(serial_scalar_array_t), pointer      :: strong_dirichlet_values            => null() !< Strong dirichlet values
        type(finite_element_t),      pointer      :: fe                                 => null() !< finite element
        class(reference_fe_t),       pointer      :: reference_fe_phy_origin            => null() !< reference finite element
        class(vector_t),             pointer      :: fe_function_dof_values             => null() !< dof values of the fe_function
        type(i1p_t),                 pointer      :: elem2dof(:)                        => null() !< element 2 dof translator
        type(interpolation_t),       pointer      :: interpolation_geometry             => null() !< interpolator
        integer(ip),                 pointer      :: subelements_connectivity(:,:)      => null() !< connectivity of subelements
        integer(ip),                 pointer      :: field_blocks(:)                    => null() !< field blocks
        real(rp),                    pointer      :: strong_dirichlet_values_entries(:) => null() !< strong dirichlet values
        type(fe_map_t)                            :: fe_map                                       !< fe map
        integer(ip)                               :: number_elements                              !< number of elements
        integer(ip)                               :: number_subelements                           !< number of subelements per element
        integer(ip)                               :: number_vertices                              !< number of geo vertex
        integer(ip)                               :: number_components                            !< number of components
        integer(ip)                               :: number_nodes_scalar                          !< number of scalar nodes
        integer(ip)                               :: number_nodes_fe_space                        !< number of nodes per fe space
        integer(ip)                               :: fe_space_index                               !< fe_space index
        integer(ip)                               :: element_index                                !< element index
        integer(ip)                               :: component_index                              !< component index
        integer(ip)                               :: node_index                                   !< node index
        integer(ip)                               :: subelement_index                             !< subelement index
        integer(ip)                               :: subnode_index                                !< subelement node index
        integer(ip)                               :: max_order                                    !< max fe_space order
        integer(ip)                               :: max_order_fe_space_index                     !< index of the fe_space with max order
        integer(ip)                               :: idx                                          !< Indices
        integer(ip)                               :: E_IO                                         !< IO Error
        logical                                   :: ft                                           !< Fine Task
    !-----------------------------------------------------------------
        ! Get max order and its fe_space_index
        max_order = 0
        do fe_space_index=1,this%fe_space%get_number_fe_spaces()
            if(max_order<this%fe_space%get_order(fe_space_index)) max_order_fe_space_index = fe_space_index
            max_order = max(max_order,this%fe_space%get_order(fe_space_index))
        end do

        ! Getters
        nodal_quadrature_target => this%fe_space%get_nodal_quadrature(max_order_fe_space_index)

        ! Calculate number nodes of target field given a certain order of interpolation, the fe_type and the number of components
        number_components = this%fe_space%get_number_field_components(max_order_fe_space_index)
        number_elements = this%fe_space%get_number_elements()
        number_subelements = this%fe_space%get_number_subelements(max_order_fe_space_index)
        number_nodes_scalar = this%fe_space%get_number_nodes_scalar(max_order_fe_space_index)
        call memalloc(number_nodes_scalar*number_components, nodal_values_target, __FILE__, __LINE__)

        subelements_connectivity => this%mesh(mesh_number)%get_subelements_connectivity()

        do fe_space_index=1, this%fe_space%get_number_fe_spaces()
            number_nodes_fe_space = this%fe_space%get_number_nodes(fe_space_index)
            call memalloc(number_nodes_fe_space, nodal_values_origin, __FILE__, __LINE__)

            reference_fe_phy_origin => this%fe_space%get_reference_fe_phy(fe_space_index)
            number_vertices = reference_fe_phy_origin%get_number_vertices()
            call fe_map%create(nodal_quadrature_target, reference_fe_phy_origin)
            call memalloc(number_components, number_subelements*number_elements*number_vertices , field, __FILE__, __LINE__)

            field_blocks => this%fe_space%get_field_blocks()

            ! Extract nodal values associated to dirichlet bcs
            strong_dirichlet_values => fe_function%get_strong_dirichlet_values()
            strong_dirichlet_values_entries => strong_dirichlet_values%get_entries()

            do element_index=1, number_elements
                fe => this%fe_space%get_finite_element(element_index)  
                elem2dof => fe%get_elem2dof()

                 ! Extract nodal values associated to dofs
                fe_function_dof_values => fe_function%get_dof_values()
                call fe_function_dof_values%extract_subvector ( field_blocks(fe_space_index), &
                                                                & number_nodes_fe_space,      &
                                                                & elem2dof(fe_space_index)%p, &
                                                                & nodal_values_origin )

                ! Fill nodal values with strong dirichlet values
                do node_index = 1, number_nodes_fe_space
                    if ( elem2dof(fe_space_index)%p(node_index) < 0 ) then
                        nodal_values_origin(node_index) = strong_dirichlet_values_entries(abs(elem2dof(fe_space_index)%p(node_index)))
                    end if
                end do

                ! interpolate nodal values if needed
                if(this%fe_space%get_order(fe_space_index) /= max_order) then
                    interpolation_geometry => fe_map%get_interpolation_geometry()
                    call reference_fe_phy_origin%interpolate_nodal_values( &
                                                        interpolation_geometry, &
                                                        nodal_values_origin, &
                                                        nodal_values_target ) 
                else
                    assert(number_nodes_fe_space==number_nodes_scalar*number_components)
                    nodal_values_target=nodal_values_origin
                endif

                ! Build field in VTK-like format
                ! Loop over number of components
                do component_index=1, number_components
                    ! Loop over subelements
                    do subelement_index = 1, number_subelements
                        ! Loop over geometrical nodes in subelement
                        do node_index=1, number_vertices
                            idx = subelements_connectivity(node_index,subelement_index)
                            subnode_index=( (element_index-1)*number_subelements+(subelement_index-1) )*number_vertices+node_index
                            field(component_index , subnode_index) = nodal_values_target(idx+(component_index-1)*number_nodes_scalar)
                        enddo
                    end do
                end do
            enddo

            call fe_map%free()
            E_IO = this%mesh(mesh_number)%write_node_field(field, field_name, f_id)
            call memfree(nodal_values_origin, __FILE__, __LINE__)
            call memfree(field, __FILE__, __LINE__)
        enddo
        call memfree(nodal_values_target, __FILE__, __LINE__)
    end function vtk_write_node_field_superlinear


    function vtk_end_write(this, mesh_number, f_id) result(E_IO)
    !-----------------------------------------------------------------
    !< Ends the writing of a single VTK file to disk (if I am fine MPI task)
    !< Closes geometry ( VTK_END_XML, VTK_GEO_XML )
    !-----------------------------------------------------------------
        implicit none
        class(vtk_handler_t),       intent(INOUT) :: this        !< vtk_handler_t derived type
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
           
           E_IO = this%mesh(nm)%end_write()
        endif
    end function vtk_end_write


    subroutine VTK_free (this)
    !-----------------------------------------------------------------
    !< Free the vtk_handler_t derived type
    !-----------------------------------------------------------------
        class(vtk_handler_t), intent(inout) :: this
        integer(ip)                         :: i, j 
        logical                             :: ft
    !-----------------------------------------------------------------
        check(associated(this%env))
        ft = this%env%am_i_fine_task() 
    
        if(ft) then
            if(allocated(this%mesh)) then
                do i=1, size(this%mesh,1)
                    call this%mesh(i)%free()
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


end module vtk_handler_names
