
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
USE iso_fortran_env,             only: error_unit
USE types_names
USE list_types_names
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
USE lib_vtk_io


implicit none
#include "debug.i90"

private

    ! VTK cell type parameters
    integer(ip), parameter :: vtk_vertex               = 1
    integer(ip), parameter :: vtk_poly_vertex          = 2
    integer(ip), parameter :: vtk_line                 = 3
    integer(ip), parameter :: vtk_poly_line            = 4
    integer(ip), parameter :: vtk_triangle             = 5
    integer(ip), parameter :: vtk_triangle_strip       = 6
    integer(ip), parameter :: vtk_polygon              = 7
    integer(ip), parameter :: vtk_pixel                = 8
    integer(ip), parameter :: vtk_quad                 = 9
    integer(ip), parameter :: vtk_tetra                = 10
    integer(ip), parameter :: vtk_voxel                = 11
    integer(ip), parameter :: vtk_hexahedron           = 12
    integer(ip), parameter :: vtk_wedge                = 13
    integer(ip), parameter :: vtk_pyramid              = 14
    integer(ip), parameter :: vtk_quadratic_edge       = 21
    integer(ip), parameter :: vtk_quadratic_triangle   = 22
    integer(ip), parameter :: vtk_quadratic_quad       = 23
    integer(ip), parameter :: vtk_quadratic_tetra      = 24
    integer(ip), parameter :: vtk_quadratic_hexahedron = 25
  
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
        procedure, public :: write_pvtk                   => vtk_write_pvtk
        procedure, public :: write_pvd                    => vtk_write_pvd
        procedure, public :: free                         => vtk_free
        procedure         :: reallocate_meshes            => vtk_reallocate_meshes
        procedure         :: fill_mesh_linear_order       => vtk_fill_mesh_linear_order
        procedure         :: fill_mesh_superlinear_order  => vtk_fill_mesh_superlinear_order
        procedure         :: write_node_field_linear      => vtk_write_node_field_linear
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
        assert(associated(this%env))
        call this%env%info(me,np) 
        assert(this%root_proc <= np-1)

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


    subroutine vtk_reallocate_meshes(this, mesh_number)
    !-----------------------------------------------------------------
    !< Reallocate meshes array to fit a new mesh
    !-----------------------------------------------------------------
        class(vtk_handler_t), intent(INOUT) :: this
        integer(ip),          intent(OUT)   :: mesh_number
        type(vtk_mesh_t), allocatable       :: f_vtk_tmp(:)
        integer(ip)                         :: i
    !-----------------------------------------------------------------
        ! Meshes allocation
        if(this%num_meshes == 0) then 
            this%num_meshes = 1
            if(allocated(this%mesh)) deallocate(this%mesh)
            allocate(this%mesh(this%num_meshes))
        else
            allocate(f_vtk_tmp(this%num_meshes))
            do i=1, this%num_meshes
                call this%mesh(i)%move_to(f_vtk_tmp(i))
            enddo
            deallocate(this%mesh)

            allocate(this%mesh(this%num_meshes+1))
            do i=1, this%num_meshes
                call f_vtk_tmp(i)%move_to(this%mesh(i))
            enddo
            this%num_meshes = this%num_meshes + 1 
            deallocate(f_vtk_tmp)
        endif
        mesh_number = this%num_meshes
    end subroutine vtk_reallocate_meshes


    subroutine vtk_initialize(this, fe_space, env, path, prefix, root_proc, number_of_parts, number_of_steps, linear_order, mesh_number)
    !-----------------------------------------------------------------
    !< Initialize the vtk_handler_t derived type
    !-----------------------------------------------------------------
        class(vtk_handler_t),             intent(INOUT) :: this
        class(serial_fe_space_t), target, intent(INOUT) :: fe_space
        class(environment_t),     target, intent(IN)    :: env
        character(len=*),                 intent(IN)    :: path
        character(len=*),                 intent(IN)    :: prefix  
        integer(ip),      optional,       intent(IN)    :: root_proc
        integer(ip),      optional,       intent(IN)    :: number_of_parts
        integer(ip),      optional,       intent(IN)    :: number_of_steps
        logical,          optional,       intent(IN)    :: linear_order
        integer(ip),      optional,       intent(OUT)   :: mesh_number
        integer(ip)                                     :: nm
        logical                                         :: lo, ft
        integer(ip)                                     :: me, np, st, rp
    !-----------------------------------------------------------------
        lo = .False. !< Default linear order = .false.
        ft = .False. !< Default fine task = .false.

        if(present(linear_order)) lo = linear_order

        this%env => env

        me = 0; np = 1; rp = 0
        if(associated(this%env)) then 
            call this%env%info(me,np) 
            ft =  this%env%am_i_fine_task() 
        endif
        if(present(root_proc)) rp = root_proc
        call this%set_root_proc(rp)

        if(ft) then
            call this%reallocate_meshes(nm)
            call this%mesh(nm)%set_fe_space(fe_space)
            call this%mesh(nm)%set_linear_order(lo)

            if(this%mesh(nm)%is_linear_order()) then 
                call this%fill_mesh_linear_order(mesh_number=nm)
            else
                call this%fill_mesh_superlinear_order(mesh_number=nm)
            endif

            call this%set_path(path,nm)
            call this%set_prefix(prefix,nm)
            if(present(number_of_parts)) np = number_of_parts
            call this%set_num_parts(np)
            st = 1
            if(present(number_of_steps)) st = number_of_steps
            call this%set_num_steps(st)

            if(present(mesh_number)) mesh_number = nm
        endif
    end subroutine vtk_initialize


    subroutine vtk_fill_mesh_linear_order(this, mesh_number)
    !-----------------------------------------------------------------
    !< Store a linear_order mesh in a vtk_handler_t derived type from a triangulation
    !-----------------------------------------------------------------
        class(vtk_handler_t),             intent(INOUT) :: this
        integer(ip),                      intent(IN)    :: mesh_number
        type(triangulation_t), pointer                  :: triangulation
        integer(ip)                                     :: i
        integer(ip)                                     :: j
        integer(ip)                                     :: number_nodes
        integer(ip)                                     :: counter
    !-----------------------------------------------------------------
        triangulation => this%mesh(mesh_number)%get_triangulation()
        assert(associated(triangulation))

        call this%mesh(mesh_number)%set_dimensions(triangulation%num_dims)
        call this%mesh(mesh_number)%set_number_elements(triangulation%num_elems)
        call this%mesh(mesh_number)%allocate_elemental_arrays()

        ! Fill VTK cell type and offset arrays and and count nodes
        number_nodes = 0
        do i=1, this%mesh(mesh_number)%get_number_elements()
            if(triangulation%elems(i)%reference_fe_geo%get_topology() == topology_quad) then
                call this%mesh(mesh_number)%set_cell_type(index=i, type=vtk_pixel) 
            else
                write(error_unit,*) 'fill_mesh_from_triangulation: Topology not supported'
                check(.false.)
            endif
            number_nodes = number_nodes + triangulation%elems(i)%reference_fe_geo%get_number_vertices()
            call this%mesh(mesh_number)%set_offset(i, number_nodes)
        enddo
        call this%mesh(mesh_number)%set_number_nodes(number_nodes)
        call this%mesh(mesh_number)%allocate_nodal_arrays()
        call this%mesh(mesh_number)%initialize_coordinates()

        ! Fill VTK coordinate arrays
        counter = 1
        do i=1, this%mesh(mesh_number)%get_number_elements()
            do j=1, triangulation%elems(i)%reference_fe_geo%get_number_vertices()
                call this%mesh(mesh_number)%set_connectivity(counter, counter-1)
                if (triangulation%num_dims >=1) call this%mesh(mesh_number)%set_x_coordinate(counter, triangulation%elems(i)%coordinates(1,j))
                if (triangulation%num_dims >=2) call this%mesh(mesh_number)%set_y_coordinate(counter, triangulation%elems(i)%coordinates(2,j))
                if (triangulation%num_dims >=3) call this%mesh(mesh_number)%set_z_coordinate(counter, triangulation%elems(i)%coordinates(3,j))
                counter = counter + 1
            enddo
            number_nodes = number_nodes + triangulation%elems(i)%reference_fe_geo%get_number_vertices()
        enddo
        call this%mesh(mesh_number)%set_filled(.true.)
    end subroutine vtk_fill_mesh_linear_order


    subroutine vtk_fill_mesh_superlinear_order(this, mesh_number)
    !-----------------------------------------------------------------
    !< Store a superlinear_order mesh in a vtk_handler_t derived type from a fe_space
    !-----------------------------------------------------------------
        class(vtk_handler_t),             intent(INOUT) :: this
        integer(ip), optional,            intent(IN)    :: mesh_number
        class(serial_fe_space_t), pointer               :: fe_space
        type(point_t),            pointer               :: vertex_coordinates(:)
        type(point_t),            pointer               :: nodal_coordinates(:)
        type(quadrature_t),       pointer               :: nodal_quadrature
        class(reference_fe_t),    pointer               :: reference_fe_geo
        type(finite_element_t),   pointer               :: fe
        class(reference_fe_t),    pointer               :: max_order_reference_fe_phy
        integer(ip),              pointer               :: subelements_connectivity(:,:)
        type(fe_map_t)                                  :: fe_map
        integer(ip)                                     :: num_elements
        integer(ip)                                     :: num_nodes_per_element
        integer(ip)                                     :: num_vertices_per_element
        integer(ip)                                     :: num_subelements_per_element
        integer(ip)                                     :: elements_counter
        integer(ip)                                     :: vertex
        integer(ip)                                     :: max_order
        integer(ip)                                     :: dimensions
        integer(ip)                                     :: subelement_vertex
        integer(ip)                                     :: fe_space_index
        integer(ip)                                     :: element_index
        integer(ip)                                     :: subelement_index
        integer(ip)                                     :: max_order_fe_space_index
        integer(ip)                                     :: nodes_counter
    !-----------------------------------------------------------------
        fe_space => this%mesh(mesh_number)%get_fe_space()
        assert(associated(fe_space))

        nullify(vertex_coordinates)
        nullify(nodal_coordinates)
        nullify(nodal_quadrature)
        nullify(reference_fe_geo)
        nullify(fe)
        nullify(max_order_reference_fe_phy)
        nullify(subelements_connectivity)

        ! Get max order and its fe_space_index
        max_order                = fe_space%get_max_order()
        max_order_fe_space_index = fe_space%get_max_order_fe_space_component()
        
        ! Getters
        max_order_reference_fe_phy => fe_space%get_reference_fe_phy(max_order_fe_space_index)
        reference_fe_geo           => fe_space%get_reference_fe_geo()
        nodal_quadrature           => max_order_reference_fe_phy%get_nodal_quadrature()

        ! Get dimensions from the geometric reference finite element
        dimensions = reference_fe_geo%get_number_dimensions()
        call this%mesh(mesh_number)%set_dimensions(dimensions)
     
        ! Calculate the number of subelems and points for the postprocess
        num_elements                = fe_space%get_number_elements()
        num_vertices_per_element    = reference_fe_geo%get_number_vertices()
        num_subelements_per_element = max_order_reference_fe_phy%get_number_subelements()
        call this%mesh(mesh_number)%set_number_elements(num_elements * num_subelements_per_element)
        call this%mesh(mesh_number)%set_number_nodes(num_vertices_per_element * this%mesh(mesh_number)%get_number_elements())

        ! Allocate VTK geometry and connectivity data
        call this%mesh(mesh_number)%allocate_elemental_arrays()
        call this%mesh(mesh_number)%allocate_nodal_arrays()
        call this%mesh(mesh_number)%initialize_coordinates()
        
        ! Get the connectivity of the subelements
        call this%mesh(mesh_number)%allocate_subelements_connectivity(num_vertices_per_element, num_subelements_per_element)
        subelements_connectivity => this%mesh(mesh_number)%get_subelements_connectivity()
        call max_order_reference_fe_phy%get_subelements_connectivity(subelements_connectivity)

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
                    if(dimensions>=1) call this%mesh(mesh_number)%set_x_coordinate(nodes_counter, nodal_coordinates(subelement_vertex)%get(1))
                    if(dimensions>=2) call this%mesh(mesh_number)%set_y_coordinate(nodes_counter, nodal_coordinates(subelement_vertex)%get(2))
                    if(dimensions>=3) call this%mesh(mesh_number)%set_z_coordinate(nodes_counter, nodal_coordinates(subelement_vertex)%get(3))
                    call this%mesh(mesh_number)%set_connectivity(nodes_counter, nodes_counter-1)
                end do

                ! Store the type of element
                call this%mesh(mesh_number)%set_cell_type(elements_counter, type=vtk_pixel)

                ! Fill offset
                call this%mesh(mesh_number)%set_offset(elements_counter, nodes_counter)
            end do
        end do

        call fe_map%free()
    end subroutine vtk_fill_mesh_superlinear_order


    function vtk_begin_write(this, file_name, part_number, time_step, mesh_number, format) result(E_IO)
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
        assert(associated(this%env))
     
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
                E_IO = this%mesh(nm)%begin_write(fn, np, ts, format)
            endif
        endif
    end function vtk_begin_write


    function vtk_write_node_field(this, fe_function, fe_space_index, field_name, mesh_number) result(E_IO)
    !-----------------------------------------------------------------
    !< Write node field to file
    !-----------------------------------------------------------------
        class(vtk_handler_t),       intent(INOUT) :: this           !< vtk_handler_t derived type
        type(fe_function_t),        intent(INOUT) :: fe_function    !< fe_function containing the field to be written
        integer(ip),                intent(IN)    :: fe_space_index !< Fe space index
        character(len=*),           intent(IN)    :: field_name     !< name of the field
        integer(ip),      optional, intent(IN)    :: mesh_number    !< Real Number of Mesh
        integer(ip)                               :: nm             !< Aux number of mesh
        logical                                   :: ft             !< Fine Task
        integer(ip)                               :: E_IO           !< IO Error
    !-----------------------------------------------------------------
        assert(associated(this%env))
        E_IO = 0
        ft =  this%env%am_i_fine_task() 
        if(ft) then        
           nm = this%num_meshes
           if(present(mesh_number)) nm = mesh_number

            if(this%mesh(nm)%is_linear_order()) then
                E_IO = this%write_node_field_linear(fe_function, fe_space_index, field_name, nm)
            else
                E_IO = this%write_node_field_superlinear(fe_function, fe_space_index, field_name, nm)
            endif

        endif
    end function vtk_write_node_field


    function vtk_write_node_field_linear(this, fe_function, fe_space_index, field_name, mesh_number) result(E_IO)
    !-----------------------------------------------------------------
    !< Write linear field to file
    !-----------------------------------------------------------------
        implicit none
        class(vtk_handler_t),       intent(INOUT) :: this                                         !< vtk_handler_t derived type
        type(fe_function_t),        intent(INOUT) :: fe_function                                  !< Postprocess field structure to be written
        integer(ip),                intent(IN)    :: fe_space_index                               !< Fe space index
        character(len=*),           intent(IN)    :: field_name                                   !< name of the field
        integer(ip),      optional, intent(IN)    :: mesh_number                                  !< Real Number of Mesh
        real(rp), allocatable                     :: field(:,:)                                   !< FIELD(ncomp,nnod)
        real(rp), allocatable                     :: nodal_values(:)                              !< nodal values
        class(serial_fe_space_t),    pointer      :: fe_space                                     !< fe space
        type(serial_scalar_array_t), pointer      :: strong_dirichlet_values                      !< Strong dirichlet values
        type(finite_element_t),      pointer      :: fe                                           !< finite element
        class(reference_fe_t),       pointer      :: reference_fe_phy_origin                      !< reference finite element
        class(vector_t),             pointer      :: fe_function_dof_values                       !< dof values of the fe_function
        type(i1p_t),                 pointer      :: elem2dof(:)                                  !< element 2 dof translator
        integer(ip),                 pointer      :: field_blocks(:)                              !< field blocks
        real(rp),                    pointer      :: strong_dirichlet_values_entries(:)           !< strong dirichlet values
        type(list_t),                pointer      :: nodes_vef                                    !< list of reference_fe_phy nodes
        type(list_iterator_t)                     :: nodes_vertex_iterator                        !< iterator on vertex nodes of the reference_fe_phy
        integer(ip)                               :: number_elements                              !< number of elements
        integer(ip)                               :: number_vertices                              !< number of geo vertex
        integer(ip)                               :: number_components                            !< number of components
        integer(ip)                               :: number_nodes                                 !< number of nodes per fe space
        integer(ip)                               :: element_index                                !< element index
        integer(ip)                               :: component_index                              !< component index
        integer(ip)                               :: vertex_index                                 !< vertex index
        integer(ip)                               :: node_index                                   !< node index
        integer(ip)                               :: E_IO                                         !< IO Error
    !-----------------------------------------------------------------
        fe_space => this%mesh(mesh_number)%get_fe_space()
        assert(associated(fe_space))
        assert(this%mesh(mesh_number)%is_linear_order())

        nullify(strong_dirichlet_values)
        nullify(fe)
        nullify(reference_fe_phy_origin)
        nullify(fe_function_dof_values)
        nullify(elem2dof)
        nullify(field_blocks)
        nullify(strong_dirichlet_values_entries)

        ! Point to some fe_space content
        reference_fe_phy_origin => fe_space%get_reference_fe_phy(fe_space_index)
        field_blocks            => fe_space%get_field_blocks()
        nodes_vef               => reference_fe_phy_origin%get_nodes_vef()


        ! Extract nodal values associated to dirichlet bcs and dof values
        strong_dirichlet_values         => fe_function%get_strong_dirichlet_values()
        strong_dirichlet_values_entries => strong_dirichlet_values%get_entries()
        fe_function_dof_values          => fe_function%get_dof_values()

        ! Get number components, elements, nodes and vertices
        number_elements   = fe_space%get_number_elements()
        number_components = reference_fe_phy_origin%get_number_field_components()
        number_nodes      = reference_fe_phy_origin%get_number_nodes()
        number_vertices   = reference_fe_phy_origin%get_number_vertices()

        ! Allocate nodal values per finite element and VTK field
        call memalloc(number_nodes, nodal_values, __FILE__, __LINE__)
        call memalloc(number_components, number_elements*number_vertices , field, __FILE__, __LINE__)

        do element_index=1, number_elements
            fe => fe_space%get_finite_element(element_index)  
            elem2dof => fe%get_elem2dof()

            ! Extract nodal values associated to dofs
            call fe_function_dof_values%extract_subvector ( field_blocks(fe_space_index), number_nodes, elem2dof(fe_space_index)%p, nodal_values )

            ! Build field in VTK-like format
            ! Loop on geometrical nodes in subelement
            do vertex_index=1, number_vertices
                nodes_vertex_iterator = nodes_vef%get_iterator(vertex_index)
                assert(nodes_vertex_iterator%get_size() == number_components)
                ! Loop on field components
                do while(.not. nodes_vertex_iterator%is_upper_bound())
                    node_index = nodes_vertex_iterator%get_current()
                    if ( elem2dof(fe_space_index)%p(node_index) < 0 ) then
                        ! Fill field with strong dirichlet values
                        field(reference_fe_phy_origin%get_component_node(node_index), (element_index-1)*number_vertices+vertex_index) = &
                            strong_dirichlet_values_entries(-elem2dof(fe_space_index)%p(node_index))
                    else
                        ! Fill field with nodal values
                        field(reference_fe_phy_origin%get_component_node(node_index), (element_index-1)*number_vertices+vertex_index) = nodal_values(node_index)
                    endif
                    call nodes_vertex_iterator%next()
                end do
            end do
        enddo

        E_IO = this%mesh(mesh_number)%write_node_field(fe_space_index, field, field_name)
        call memfree(nodal_values, __FILE__, __LINE__)
        call memfree(field, __FILE__, __LINE__)
    end function vtk_write_node_field_linear


    function vtk_write_node_field_superlinear(this, fe_function, fe_space_index, field_name, mesh_number) result(E_IO)
    !-----------------------------------------------------------------
    !< Write superlinear field to file
    !-----------------------------------------------------------------
        implicit none
        class(vtk_handler_t),       intent(INOUT) :: this                                         !< vtk_handler_t derived type
        type(fe_function_t),        intent(INOUT) :: fe_function                                  !< Postprocess field structure to be written
        integer(ip),                intent(IN)    :: fe_space_index                               !< Fe space index
        character(len=*),           intent(IN)    :: field_name                                   !< name of the field
        integer(ip),      optional, intent(IN)    :: mesh_number                                  !< Real Number of 
        real(rp), allocatable                     :: field(:,:)                                   !< FIELD(ncomp,nnod)
        real(rp), allocatable                     :: nodal_values_origin(:)                       !< nodal values of the origin fe_space
        real(rp), allocatable                     :: nodal_values_target(:)                       !< nodal values for the interpolation
        class(serial_fe_space_t),    pointer      :: fe_space                                     !< fe space
        type(quadrature_t),          pointer      :: nodal_quadrature_target                      !< Nodal quadrature
        type(serial_scalar_array_t), pointer      :: strong_dirichlet_values                      !< Strong dirichlet values
        type(finite_element_t),      pointer      :: fe                                           !< finite element
        class(reference_fe_t),       pointer      :: reference_fe_phy_origin                      !< reference finite element
        class(reference_fe_t),       pointer      :: reference_fe_phy_target                      !< reference finite element
        class(vector_t),             pointer      :: fe_function_dof_values                       !< dof values of the fe_function
        type(i1p_t),                 pointer      :: elem2dof(:)                                  !< element 2 dof translator
        type(interpolation_t)                     :: interpolation                                !< interpolator
        integer(ip),                 pointer      :: subelements_connectivity(:,:)                !< connectivity of subelements
        integer(ip),                 pointer      :: field_blocks(:)                              !< field blocks
        real(rp),                    pointer      :: strong_dirichlet_values_entries(:)           !< strong dirichlet values
        integer(ip)                               :: number_elements                              !< number of elements
        integer(ip)                               :: number_subelements                           !< number of subelements per element
        integer(ip)                               :: number_vertices                              !< number of geo vertex
        integer(ip)                               :: number_components                            !< number of components
        integer(ip)                               :: number_nodes_scalar                          !< number of scalar nodes
        integer(ip)                               :: number_nodes                                 !< number of nodes per fe space
        integer(ip)                               :: element_index                                !< element index
        integer(ip)                               :: component_index                              !< component index
        integer(ip)                               :: node_index                                   !< node index
        integer(ip)                               :: subelement_index                             !< subelement index
        integer(ip)                               :: subnode_index                                !< subelement node index
        integer(ip)                               :: order                                        !< current fe_space order
        integer(ip)                               :: max_order                                    !< max fe_space order
        integer(ip)                               :: max_order_fe_space_index                     !< index of the fe_space with max order
        integer(ip)                               :: idx                                          !< Indices
        integer(ip)                               :: E_IO                                         !< IO Error
        logical                                   :: ft                                           !< Fine Task
    !-----------------------------------------------------------------
        fe_space => this%mesh(mesh_number)%get_fe_space()
        assert(associated(fe_space))
        assert(.not. this%mesh(mesh_number)%is_linear_order())

        nullify(nodal_quadrature_target)
        nullify(strong_dirichlet_values)
        nullify(fe)
        nullify(reference_fe_phy_origin)
        nullify(reference_fe_phy_target)
        nullify(fe_function_dof_values)
        nullify(elem2dof)
        nullify(subelements_connectivity)
        nullify(field_blocks)
        nullify(strong_dirichlet_values_entries)

        ! Get max order and its fe_space_index
        max_order = fe_space%get_max_order()
        max_order_fe_space_index = fe_space%get_max_order_fe_space_component()

        ! Point physical referece fe, field blocks, nodal quadrature and subelements connectivity
        reference_fe_phy_origin  => fe_space%get_reference_fe_phy(fe_space_index)
        reference_fe_phy_target  => fe_space%get_reference_fe_phy(max_order_fe_space_index)
        field_blocks             => fe_space%get_field_blocks()
        nodal_quadrature_target  => reference_fe_phy_target%get_nodal_quadrature()
        subelements_connectivity => this%mesh(mesh_number)%get_subelements_connectivity()

        ! Extract nodal values associated to dirichlet bcs and dof values
        strong_dirichlet_values         => fe_function%get_strong_dirichlet_values()
        strong_dirichlet_values_entries => strong_dirichlet_values%get_entries()
        fe_function_dof_values          => fe_function%get_dof_values()
        


        ! Calculate number nodes of target field given a certain order of interpolation, the fe_type and the number of components
        number_elements     = fe_space%get_number_elements()
        number_components   = reference_fe_phy_origin%get_number_field_components()
        number_subelements  = reference_fe_phy_target%get_number_subelements()
        number_nodes_scalar = reference_fe_phy_origin%get_number_nodes_scalar()
        number_nodes        = reference_fe_phy_origin%get_number_nodes()
        number_vertices     = reference_fe_phy_target%get_number_vertices()
        order               = reference_fe_phy_origin%get_order()

        ! Allocate VTK field, origin and target nodal values
        call memalloc(number_nodes, nodal_values_origin, __FILE__, __LINE__)
        call memalloc(number_components*nodal_quadrature_target%get_number_quadrature_points(), nodal_values_target, __FILE__, __LINE__)
        call memalloc(number_components, number_subelements*number_elements*number_vertices , field, __FILE__, __LINE__)

        ! Create interpolation from nodal quadrature
        call reference_fe_phy_origin%create_interpolation(nodal_quadrature_target, interpolation)
        
        ! Loop on elements
        subnode_index=1
        do element_index=1, number_elements
            fe => fe_space%get_finite_element(element_index)  
            elem2dof => fe%get_elem2dof()

            ! Extract nodal values associated to dofs
            call fe_function_dof_values%extract_subvector ( field_blocks(fe_space_index), number_nodes, elem2dof(fe_space_index)%p, nodal_values_origin )

            ! Fill nodal values with strong dirichlet values
            do node_index = 1, number_nodes
                if ( elem2dof(fe_space_index)%p(node_index) < 0 ) then
                    nodal_values_origin(node_index) = strong_dirichlet_values_entries(-elem2dof(fe_space_index)%p(node_index))
                end if
            end do

            ! interpolate nodal values if needed
            if(order /= max_order) then
                call reference_fe_phy_origin%interpolate_nodal_values( interpolation, nodal_values_origin, nodal_values_target ) 
            else
                nodal_values_target=nodal_values_origin
            endif

            ! Build field in VTK-like format
            ! Loop on subelements
            do subelement_index = 1, number_subelements
                ! Loop on geometrical nodes per subelement
                do node_index=1, number_vertices
                    ! Loop on components
                    do component_index=1, number_components
                        idx = subelements_connectivity(node_index,subelement_index) + (component_index-1)*number_nodes_scalar
                        field(reference_fe_phy_target%get_component_node(idx) , subnode_index) = nodal_values_target(idx)
                    enddo
                    subnode_index=subnode_index+1
                end do
            end do
        enddo

        E_IO = this%mesh(mesh_number)%write_node_field(fe_space_index, field, field_name)
        call interpolation%free()
        call memfree(nodal_values_origin, __FILE__, __LINE__)
        call memfree(field, __FILE__, __LINE__)
        call memfree(nodal_values_target, __FILE__, __LINE__)
    end function vtk_write_node_field_superlinear


    function vtk_end_write(this, mesh_number) result(E_IO)
    !-----------------------------------------------------------------
    !< Ends the writing of a single VTK file to disk (if I am fine MPI task)
    !< Closes geometry ( VTK_END_XML, VTK_GEO_XML )
    !-----------------------------------------------------------------
        class(vtk_handler_t),       intent(INOUT) :: this        !< vtk_handler_t derived type
        integer(ip),      optional, intent(IN)    :: mesh_number !< Number of MESH
        integer(ip)                               :: nm          !< Real Number of Mesh
        integer(ip)                               :: E_IO        !< IO Error
        logical                                   :: ft          !< Fine Task
      ! ----------------------------------------------------------------------------------
        assert(associated(this%env))
        ft =  this%env%am_i_fine_task() 

        E_IO = 0

        if (ft) then
           nm = this%num_meshes
           if(present(mesh_number)) nm = mesh_number
           
           E_IO = this%mesh(nm)%end_write()
        endif
    end function vtk_end_write


    function vtk_write_PVTK(this, file_name, mesh_number, time_step) result(E_IO)
    !-----------------------------------------------------------------
    !< Write the PVTK file containing the number of parts
    !-----------------------------------------------------------------
        class(vtk_handler_t),       intent(INOUT) :: this
        character(len=*), optional, intent(IN)    :: file_name
        integer(ip),      optional, intent(IN)    :: mesh_number
        real(rp),         optional, intent(IN)    :: time_step
        integer(ip)                               :: nm, rf
        character(len=:),allocatable              :: path
        character(len=:),allocatable              :: prefix
        character(len=:),allocatable              :: var_name
        character(len=:),allocatable              :: field_type
        character(len=:),allocatable              :: fn ,dp
        real(rp)                                  :: ts
        integer(ip)                               :: i, fid, nnods, nels, E_IO
        integer(ip)                               :: me, np
        logical                                   :: isDir
    !-----------------------------------------------------------------

        me = 0; np = 1
        check(associated(this%env))
        call this%env%info(me,np) 

        E_IO = 0

        if( this%env%am_i_fine_task() .and. me == this%root_proc) then

            nm = this%num_meshes
            if(present(mesh_number)) nm = mesh_number
            ts = 0_rp
            if(allocated(this%steps)) then
                if(this%steps_counter >0 .and. this%steps_counter <= size(this%steps,1)) &
                    ts = this%steps(this%steps_counter)
            endif
            if(present(time_step)) ts = time_step 

            call this%mesh(nm)%get_path(path)
            call this%mesh(nm)%get_prefix(prefix)

            dp = this%get_VTK_time_output_path(path=path, time_step=ts, mesh_number=nm)
            fn = this%get_PVTK_filename(prefix=prefix, mesh_number=nm, time_step=ts)
            fn = dp//fn
            if(present(file_name)) fn = file_name

            nnods = this%mesh(nm)%get_number_nodes()
            nels = this%mesh(nm)%get_number_elements()

    !        inquire( file=trim(dp)//'/.', exist=isDir ) 
            if(this%create_directory(trim(dp), issue_final_barrier=.False.) == 0) then
                ! pvtu
                E_IO = PVTK_INI_XML(filename = trim(adjustl(fn)), mesh_topology = 'PUnstructuredGrid', tp='Float64', cf=rf)
                do i=0, this%num_parts-1
                    E_IO = PVTK_GEO_XML(source=trim(adjustl(this%get_VTK_filename(prefix=prefix, part_number=i, mesh_number=nm))), cf=rf)
                enddo

                E_IO = PVTK_DAT_XML(var_location = 'Node', var_block_action = 'OPEN', cf=rf)
                ! Write fields point data
                do i=1, this%mesh(nm)%get_number_fields()
                    if(this%mesh(nm)%field_is_filled(i)) then
                        call this%mesh(nm)%get_field_name(i, var_name)
                        call this%mesh(nm)%get_field_type(i, field_type)
                        E_IO = PVTK_VAR_XML(varname = trim(adjustl(var_name)), tp=trim(adjustl(field_type)), Nc=this%mesh(nm)%get_field_number_components(i) , cf=rf)
                    endif
                enddo
                E_IO = PVTK_DAT_XML(var_location = 'Node', var_block_action = 'CLOSE', cf=rf)
                E_IO = PVTK_END_XML(cf=rf)
            endif

        endif
    end function vtk_write_PVTK


    function vtk_write_PVD(this, file_name, mesh_number) result(E_IO)
    !-----------------------------------------------------------------
    !< Write the PVD file referencing several PVTK files in a timeline
    !< (only root processor)
    !-----------------------------------------------------------------
        class(vtk_handler_t),       intent(INOUT) :: this
        character(len=*), optional, intent(IN)    :: file_name
        integer(ip),      optional, intent(IN)    :: mesh_number
        integer(ip)                               :: nm, rf
        character(len=:),allocatable              :: path
        character(len=:),allocatable              :: prefix
        character(len=:),allocatable              :: var_name
        character(len=:),allocatable              :: pvdfn, pvtkfn ,dp
        integer(ip)                               :: i, fid, nnods, nels, ts, E_IO
        integer(ip)                               :: me, np
        logical                                   :: isDir
    !-----------------------------------------------------------------
        me = 0; np = 1
        check(associated(this%env))
        call this%env%info(me,np) 

        E_IO = 0

        if(this%env%am_i_fine_task() .and. me == this%root_proc) then
            nm = this%num_meshes
            if(present(mesh_number)) nm = mesh_number

            call this%mesh(nm)%get_path(path)
            call this%mesh(nm)%get_prefix(prefix)
        
            pvdfn = trim(adjustl(path))//'/'//trim(adjustl(prefix))//pvd_ext
            if(present(file_name)) pvdfn = file_name

    !        inquire( file=trim(adjustl(this%mesh(nm)%dir_path))//'/.', exist=isDir )
            if(this%create_directory(trim(adjustl(path)), issue_final_barrier=.False.) == 0) then
                if(allocated(this%steps)) then
                    if(size(this%steps,1) >= min(this%num_steps,this%steps_counter)) then
                        E_IO = PVD_INI_XML(filename=trim(adjustl(pvdfn)),cf=rf)
                        do ts=1, min(this%num_steps,this%steps_counter)
                            dp = this%get_PVD_time_output_path(path=path, time_step=this%steps(ts))
                            pvtkfn = this%get_PVTK_filename(prefix=prefix, mesh_number=nm, time_step=this%steps(ts))
                            pvtkfn = dp//pvtkfn
                            E_IO = PVD_DAT_XML(filename=trim(adjustl(pvtkfn)),timestep=ts, cf=rf)
                        enddo
                        E_IO = PVD_END_XML(cf=rf)
                    endif
                endif
            endif
        endif
    end function vtk_write_PVD


    subroutine VTK_free (this)
    !-----------------------------------------------------------------
    !< Free the vtk_handler_t derived type
    !-----------------------------------------------------------------
        class(vtk_handler_t), intent(inout) :: this
        integer(ip)                         :: i, j 
        logical                             :: ft
    !-----------------------------------------------------------------
        assert(associated(this%env))
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
        this%env => NULL()
    end subroutine VTK_free


end module vtk_handler_names
