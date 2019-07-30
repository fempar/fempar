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
module mesh_partitioner_names
  use types_names
  use list_types_names
  use allocatable_array_names
  use hash_table_names
  use rcm_renumbering_names
  use memor_names
  use stdio_names
  use metis_names
  use mesh_names
  use mesh_distribution_names
  use mesh_partitioner_parameters_names
  use environment_names
  use lib_vtk_io
  use vtk_utils_names
  use vtk_parameters_names
  use std_vector_names
  use FPL
  implicit none
# include "debug.i90"
  private
  
  type mesh_partitioner_t
     private
     
     type(mesh_t), pointer                      :: mesh      ! Global mesh to be partitioned
     
     integer(ip)                                :: nparts
     integer(ip)                                :: num_levels
     integer(ip)                  , allocatable :: num_parts_x_level (:)
     type(allocatable_array_ip1_t), allocatable :: parts_mapping (:)
     
     type(mesh_t)                 , allocatable :: part_meshes(:)
     type(mesh_distribution_t)    , allocatable :: mesh_distributions(:)
     type(list_t)                 , allocatable :: graph_hierarchy(:)
     type(allocatable_array_ip1_t), allocatable :: cells_part(:)
     
     character(len=:), allocatable :: strat                   ! Partitioning algorithm (part_kway,part_recursive,part_strip,part_rcm_strip)
     integer(ip)                   :: metis_option_ufactor    ! Imbalance tol of metis_option_ufactor/1000 + 1
     integer(ip)                   :: metis_option_minconn 
     integer(ip)                   :: metis_option_contig  
     integer(ip)                   :: metis_option_ctype 
     integer(ip)                   :: metis_option_iptype
     integer(ip)                   :: metis_option_debug
     character(len=:), allocatable :: vtk_format
  contains 
     procedure, non_overridable                   :: create         => mesh_partitioner_create
     procedure, non_overridable                   :: partition_mesh => mesh_partitioner_partition_mesh
     procedure, non_overridable, private          :: set_mesh
     procedure, non_overridable, private          :: set_default_parameter_values
     procedure, non_overridable, private          :: set_parameters_from_pl
     procedure, non_overridable, private          :: allocate_graph_hierarchy
     procedure, non_overridable, private          :: allocate_cells_part
     procedure, non_overridable, private          :: allocate_parts_mapping
     procedure, non_overridable, private          :: setup_parts_mapping
     procedure, non_overridable, private          :: build_and_partition_graph_hierarchy
     procedure, non_overridable, private          :: partition_graph
     procedure, non_overridable, private, nopass  :: build_parts_graph
     procedure, non_overridable, private          :: setup_part_l2g_cells_and_vertices
     procedure, non_overridable, private          :: setup_part_interface
     procedure, non_overridable, private          :: setup_part_mesh
     procedure, non_overridable, private          :: allocate_part_meshes
     procedure, non_overridable, private          :: allocate_mesh_distributions
     procedure, non_overridable, private          :: write_mesh_parts_fempar_gid_problem_type_format_dir_path_prefix
     procedure, non_overridable, private          :: write_mesh_parts_fempar_gid_problem_type_format_pl
     procedure, non_overridable, private          :: write_environment_related_data_file_unit
     generic                                      :: write_mesh_parts_fempar_gid_problem_type_format & 
                                                       => write_mesh_parts_fempar_gid_problem_type_format_dir_path_prefix, &
                                                          write_mesh_parts_fempar_gid_problem_type_format_pl
                                                          
     procedure, non_overridable, private          :: write_mesh_parts_gid_postprocess_format_dir_path_prefix
     procedure, non_overridable, private          :: write_mesh_parts_gid_postprocess_format_pl
     generic                                      :: write_mesh_parts_gid_postprocess_format & 
                                                       => write_mesh_parts_gid_postprocess_format_dir_path_prefix, &
                                                          write_mesh_parts_gid_postprocess_format_pl
                                                          
     procedure, non_overridable, private          :: vtk_output_handler_write_vtu
     procedure, non_overridable, private          :: write_mesh_partition_vtk_format_dir_path_prefix
     procedure, non_overridable, private          :: write_mesh_partition_vtk_format_pl
     generic                                      :: write_mesh_partition_vtk_format &
                                                        => write_mesh_partition_vtk_format_dir_path_prefix, &
                                                           write_mesh_partition_vtk_format_pl
                                                           
     procedure, non_overridable          :: free                                  => mesh_partitioner_free
     procedure, non_overridable, nopass  :: get_dir_path_and_prefix_from_pl       => mesh_partitioner_get_dir_path_and_prefix_from_pl
     procedure, non_overridable, private :: graph_hierarchy_free
     procedure, non_overridable, private :: cells_part_free
     procedure, non_overridable, private :: parts_mapping_free
     procedure, non_overridable, private :: num_parts_x_level_free
     procedure, non_overridable, private :: part_meshes_free
     procedure, non_overridable, private :: mesh_distributions_free
  end type mesh_partitioner_t
  
  public :: mesh_partitioner_t
  
contains

  subroutine mesh_partitioner_create ( this, mesh, parameter_list )
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    type(mesh_t)             , intent(in)    :: mesh
    type(parameterlist_t)    , intent(in)    :: parameter_list
    call this%free()
    call this%set_mesh(mesh)
    call this%set_parameters_from_pl(parameter_list)
  end subroutine mesh_partitioner_create


  !=============================================================================
  subroutine mesh_partitioner_partition_mesh ( this ) 
     implicit none
     class(mesh_partitioner_t), intent(inout) :: this
     integer(ip)  :: ipart, istat
     integer(ip)  :: num_local_vertices, num_local_cells
     integer(igp) :: num_global_vertices, num_global_cells
     integer(igp), allocatable :: l2g_vertices(:)
     integer(igp), allocatable :: l2g_cells(:)
     integer(ip) , allocatable :: lebou(:), lnbou(:)
     integer(ip) , allocatable :: pextn(:), lextp(:)
     integer(igp), allocatable :: lextn(:)
     integer(ip) :: nebou, nnbou
     
     assert ( associated(this%mesh) ) 
     call this%build_and_partition_graph_hierarchy(this%mesh)
     call this%setup_parts_mapping()
    
     call this%allocate_part_meshes()
     call this%allocate_mesh_distributions()

     num_global_vertices = int(this%mesh%get_num_vertices(),igp)
     num_global_cells    = int(this%mesh%get_num_cells(),igp)
    do ipart=1,this%nparts
       call this%setup_part_l2g_cells_and_vertices(ipart, &
                                                   num_local_vertices,  &
                                                   num_local_cells,     &
                                                   l2g_vertices,        &
                                                   l2g_cells)                          
       
       call this%setup_part_mesh(ipart, &
                                 num_local_vertices, &
                                 num_local_cells, &
                                 l2g_vertices, &
                                 l2g_cells)
       
       call this%setup_part_interface (ipart, &
                                       l2g_vertices, &
                                       l2g_cells,    &
                                       nebou,        &
                                       nnbou,        &
                                       lebou,        &
                                       lnbou,        &
                                       pextn,        &
                                       lextn,        &
                                       lextp )
       
       call this%mesh_distributions(ipart)%create (ipart, &
                                       this%nparts, &
                                       num_local_cells, &
                                       num_global_cells, &
                                       nebou, &
                                       num_local_vertices, &
                                       num_global_vertices, &
                                       nnbou, &
                                       l2g_cells, &
                                       l2g_vertices, &
                                       lebou, &
                                       lnbou, &
                                       pextn, &
                                       lextn, &
                                       lextp)
       
    end do
    call memfree(l2g_cells, __FILE__, __LINE__ )
    call memfree(l2g_vertices, __FILE__, __LINE__ )
    call memfree(lebou, __FILE__, __LINE__ )
    call memfree(lnbou, __FILE__, __LINE__ )
    call memfree(pextn, __FILE__, __LINE__ )
    call memfree(lextn, __FILE__, __LINE__ )
    call memfree(lextp, __FILE__, __LINE__ )
  end subroutine mesh_partitioner_partition_mesh
  
  subroutine allocate_graph_hierarchy(this)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    integer(ip) :: istat
    call this%graph_hierarchy_free()
    allocate(this%graph_hierarchy(this%num_levels), stat=istat); check(istat==0);
  end subroutine allocate_graph_hierarchy
  
  subroutine allocate_cells_part( this )
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    integer(ip) :: istat
    call this%cells_part_free()
    allocate(this%cells_part(this%num_levels), stat=istat); check(istat==0);
  end subroutine allocate_cells_part
  
  subroutine allocate_parts_mapping(this)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    integer(ip) :: i, istat   
    call this%parts_mapping_free()
    allocate(this%parts_mapping(sum(this%num_parts_x_level)), stat=istat); check(istat==0)
    do i=1, sum(this%num_parts_x_level)
      call this%parts_mapping(i)%create(this%num_levels)
    end do
  end subroutine allocate_parts_mapping
  
  subroutine allocate_part_meshes(this)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    integer(ip) :: istat
    call this%part_meshes_free()
    allocate(this%part_meshes(this%nparts), stat=istat); check(istat==0);
  end subroutine allocate_part_meshes
  
  subroutine allocate_mesh_distributions(this)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    integer(ip) :: istat
    call this%mesh_distributions_free()
    allocate(this%mesh_distributions(this%nparts), stat=istat); check(istat==0);
  end subroutine allocate_mesh_distributions
    
  subroutine build_and_partition_graph_hierarchy(this, mesh)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    type(mesh_t)             , intent(in)    :: mesh
    
    integer(ip), allocatable :: pelpo(:)
    integer(ip), allocatable :: lelpo(:)
    integer(ip) :: i, ilevel, istat
    
    call this%allocate_graph_hierarchy()
    call this%allocate_cells_part()
    
    ! Generate dual mesh (i.e., list of elements around points)
    call mesh%build_dual_mesh(pelpo, lelpo)
    
    ! Create dual (i.e. list of elements around elements)
    call mesh%build_dual_graph(pelpo, lelpo, this%graph_hierarchy(1))
    
    ! Free dual mesh
    call memfree(pelpo, __FILE__, __LINE__ )
    call memfree(lelpo, __FILE__, __LINE__ )
    
    ! Partition dual graph to assign a domain to each element (in cell_parts(1))
    call this%cells_part(1)%create(mesh%get_num_cells()) 
    call this%partition_graph(this%num_parts_x_level(1), &
                                   this%graph_hierarchy(1), &
                                   this%cells_part(1)%a)
    do ilevel=1,this%num_levels-1
       call this%cells_part(ilevel+1)%create(this%num_parts_x_level(ilevel)) 
       ! Typically in the last level there is only one part
       if(this%num_parts_x_level(ilevel+1)>1) then
          call this%build_parts_graph (this%num_parts_x_level(ilevel), &
                                       this%cells_part(ilevel)%a, &
                                       this%graph_hierarchy(ilevel), &
                                       this%graph_hierarchy(ilevel+1))

          call this%partition_graph(this%num_parts_x_level(ilevel+1), &
                                         this%graph_hierarchy(ilevel+1), &
                                         this%cells_part(ilevel+1)%a)
       else
          this%cells_part(ilevel+1)%a(:) = 1
       end if
    end do
  end subroutine build_and_partition_graph_hierarchy
  
  subroutine setup_parts_mapping(this)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    integer(ip) :: itask, ilevel, jlevel, ipart
    call this%allocate_parts_mapping()
    itask = 0
    do ilevel=1,this%num_levels
       do ipart = 1, this%num_parts_x_level(ilevel)
          itask = itask+1
          do jlevel = 1 , ilevel - 1 
             this%parts_mapping(itask)%a(jlevel) = 0
          end do
          this%parts_mapping(itask)%a(ilevel) = ipart
          do jlevel = ilevel+1 , this%num_levels
             this%parts_mapping(itask)%a(jlevel) = this%cells_part(jlevel)%a( this%parts_mapping(itask)%a(jlevel-1) )
          end do
       end do
    end do
  end subroutine setup_parts_mapping
  
  !=============================================================================
  subroutine set_parameters_from_pl(this,parameter_list)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    type(ParameterList_t)    , intent(in)    :: parameter_list
    ! Locals
    integer(ip)              :: istat
    integer(ip), allocatable :: param_size(:), param(:)

    ! Mandatory parameters num_levels
    assert(parameter_list%isPresent(key = mesh_partitioner_num_levels_key))
    assert(parameter_list%isAssignable(mesh_partitioner_num_levels_key, this%num_levels))
    istat = parameter_list%get(key = mesh_partitioner_num_levels_key  , value = this%num_levels)
    assert(istat==0)
       
    assert(parameter_list%isPresent(key = mesh_partitioner_num_parts_x_level_key ))
    assert(parameter_list%GetDimensions(key = mesh_partitioner_num_parts_x_level_key) == 1)

    ! Get the array using the local variable
    istat =  parameter_list%GetShape(key = mesh_partitioner_num_parts_x_level_key, shape = param_size ); check(istat==0)
    call memalloc(param_size(1), param,__FILE__,__LINE__)
    assert(parameter_list%isAssignable(mesh_partitioner_num_parts_x_level_key, param))
    istat = parameter_list%get(key = mesh_partitioner_num_parts_x_level_key, value = param)
    assert(istat==0)

    call memalloc(this%num_levels, this%num_parts_x_level,__FILE__,__LINE__)
    this%num_parts_x_level = param(1:this%num_levels)
    call memfree(param,__FILE__,__LINE__)

    this%nparts = this%num_parts_x_level(1)

    if( parameter_list%isPresent(mesh_partitioner_strategy_key) ) then
       assert(parameter_list%isAssignable(mesh_partitioner_strategy_key, this%strat))
       istat = parameter_list%get(key = mesh_partitioner_strategy_key  , value = this%strat)
       assert(istat==0)
       assert(this%strat==metis_part_kway.or.this%strat==metis_part_recursive.or.this%strat==part_strip.or.this%strat==part_rcm_strip)
    end if

    if( parameter_list%isPresent(mesh_partitioner_metis_option_debug_key) ) then
       assert(parameter_list%isAssignable(mesh_partitioner_metis_option_debug_key, this%metis_option_debug))
       istat = parameter_list%get(key = mesh_partitioner_metis_option_debug_key  , value = this%metis_option_debug)
       check(istat==0)
    end if

    if( parameter_list%isPresent(mesh_partitioner_metis_option_ufactor_key) ) then
       assert(parameter_list%isAssignable(mesh_partitioner_metis_option_ufactor_key, this%metis_option_ufactor))
       istat = parameter_list%get(key = mesh_partitioner_metis_option_ufactor_key, value = this%metis_option_ufactor)
       assert(istat==0)
    end if

    if( parameter_list%isPresent(mesh_partitioner_metis_option_minconn_key) ) then
       assert(parameter_list%isAssignable(mesh_partitioner_metis_option_minconn_key, this%metis_option_minconn))
       istat = parameter_list%get(key = mesh_partitioner_metis_option_minconn_key, value = this%metis_option_minconn)
       check(istat==0)
    end if

    if( parameter_list%isPresent(mesh_partitioner_metis_option_contig_key) ) then
       assert(parameter_list%isAssignable(mesh_partitioner_metis_option_contig_key, this%metis_option_contig))
       istat = parameter_list%get(key = mesh_partitioner_metis_option_contig_key , value = this%metis_option_contig)
       assert(istat==0)
    end if

    if( parameter_list%isPresent(mesh_partitioner_metis_option_ctype_key) ) then
       assert(parameter_list%isAssignable(mesh_partitioner_metis_option_ctype_key, this%metis_option_ctype))
       istat = parameter_list%get(key = mesh_partitioner_metis_option_ctype_key  , value = this%metis_option_ctype)
       assert(istat==0)
    end if

    this%vtk_format = mesh_partitioner_default_vtk_format
    if(parameter_list%isPresent(mesh_partitioner_vtk_format_key)) then
        assert(parameter_list%isAssignable(mesh_partitioner_vtk_format_key, 'string'))
        istat = parameter_list%GetAsString(Key=mesh_partitioner_vtk_format_key, String=this%vtk_format)
        assert(istat == 0)
    endif
  end subroutine set_parameters_from_pl
  
  !=============================================================================
  subroutine set_mesh(this, mesh)
    implicit none
    class(mesh_partitioner_t), intent(inout)      :: this
    type(mesh_t)             , target, intent(in) :: mesh
    this%mesh => mesh
  end subroutine set_mesh
  
  !=============================================================================
  subroutine set_default_parameter_values(this)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    this%strat = mesh_partitioner_default_strat
    this%metis_option_ufactor = mesh_partitioner_default_metis_option_ufactor
    this%metis_option_minconn = mesh_partitioner_default_metis_option_minconn
    this%metis_option_contig = mesh_partitioner_default_metis_option_contig
    this%metis_option_ctype = mesh_partitioner_default_metis_option_ctype
    this%metis_option_iptype = mesh_partitioner_default_metis_option_iptype
    this%metis_option_debug = mesh_partitioner_default_metis_option_debug
  end subroutine set_default_parameter_values
  
  !=============================================================================
  subroutine write_mesh_parts_fempar_gid_problem_type_format_pl(this, parameter_list)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    type(ParameterList_t)    , intent(in)    :: parameter_list 
    character(len=:), allocatable :: dir_path
    character(len=:), allocatable :: prefix
    call this%get_dir_path_and_prefix_from_pl(parameter_list, dir_path, prefix)
    call this%write_mesh_parts_fempar_gid_problem_type_format(dir_path,prefix)
  end subroutine write_mesh_parts_fempar_gid_problem_type_format_pl
  
  !=============================================================================
  subroutine write_mesh_parts_fempar_gid_problem_type_format_dir_path_prefix(this, dir_path, prefix)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    character(len=*)         , intent(in)    :: dir_path
    character(len=*)         , intent(in)    :: prefix
   
    character(len=:), allocatable  :: name, rename
    integer(ip) :: lunio
    integer(ip) :: i

    ! Write part meshes
    call this%part_meshes(1)%mesh_fempar_gid_problem_type_format_compose_name ( prefix, name )
    do i=this%nparts, 1, -1  
       rename=name
       call numbered_filename_compose(i,this%nparts,rename)
       lunio = io_open( trim(dir_path) // '/' // trim(rename), 'write' ); check(lunio>0)
       call this%part_meshes(i)%write_fempar_gid_problem_type_format(lunio)
       call io_close(lunio)
    end do
    
    ! Write mesh_distributions
    call this%mesh_distributions(1)%mesh_distribution_compose_name ( prefix, name )
    do i=this%nparts, 1, -1  
       rename=name
       call numbered_filename_compose(i,this%nparts,rename)
       lunio = io_open( trim(dir_path) // '/' // trim(rename), 'write' ); check(lunio>0)
       call this%mesh_distributions(i)%write(lunio)
       call io_close(lunio)
    end do
    
    ! Write environment related data
    call environment_compose_name ( prefix, name )
    do i=sum(this%num_parts_x_level), 1, -1  
       rename=name
       call numbered_filename_compose(i,sum(this%num_parts_x_level),rename)
       lunio = io_open( trim(dir_path) // '/' // trim(rename), 'write' ); check(lunio>0)
       call this%write_environment_related_data_file_unit(i, lunio)
       call io_close(lunio)
    end do
    
  end subroutine write_mesh_parts_fempar_gid_problem_type_format_dir_path_prefix
  
  !=============================================================================
  subroutine write_mesh_parts_gid_postprocess_format_pl(this, parameter_list)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    type(ParameterList_t)    , intent(in)    :: parameter_list 
    character(len=:), allocatable :: dir_path
    character(len=:), allocatable :: prefix
    call this%get_dir_path_and_prefix_from_pl(parameter_list, dir_path, prefix)
    call this%write_mesh_parts_gid_postprocess_format(dir_path,prefix)
  end subroutine write_mesh_parts_gid_postprocess_format_pl
  
  !=============================================================================
  subroutine write_mesh_parts_gid_postprocess_format_dir_path_prefix(this, dir_path, prefix)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    character(len=*)         , intent(in)    :: dir_path
    character(len=*)         , intent(in)    :: prefix
   
    character(len=:), allocatable  :: name, rename
    integer(ip) :: lunio
    integer(ip) :: i

    ! Write part meshes
    call this%part_meshes(1)%mesh_gid_postprocess_format_compose_name ( prefix, name )
    do i=this%nparts, 1, -1  
       rename=name
       call numbered_filename_compose(i,this%nparts,rename)
       lunio = io_open( trim(dir_path) // '/' // trim(rename), 'write' ); check(lunio>0)
       call this%part_meshes(i)%write_gid_postprocess_format(lunio)
       call io_close(lunio)
    end do
  end subroutine write_mesh_parts_gid_postprocess_format_dir_path_prefix
  
  !=============================================================================
  subroutine write_environment_related_data_file_unit (this, level, lunio)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    integer(ip)              , intent(in)    :: level
    integer(ip)              , intent(in)    :: lunio
    assert( this%num_levels>0)
    write ( lunio, '(10i10)' ) this%num_levels
    write ( lunio, '(10i10)' ) this%num_parts_x_level
    write ( lunio, '(10i10)' ) this%parts_mapping(level)%a
  end subroutine   write_environment_related_data_file_unit 
  
  !=============================================================================
    subroutine vtk_output_handler_write_vtu(this, dir_path, prefix, task_id)
    !-----------------------------------------------------------------
    !< Write the vtu mesh partitions file 
    !-----------------------------------------------------------------
        class(mesh_partitioner_t),       intent(in) :: this
        character(len=*),                intent(in) :: dir_path
        character(len=*),                intent(in) :: prefix
        integer(ip),                     intent(in) :: task_id
        character(len=:), allocatable               :: filename
        integer(ip)                                 :: num_components
        integer(ip)                                 :: file_id
        real(rp), pointer                           :: FieldValue(:,:)
        real(rp), pointer                           :: CellValue(:)
        integer(ip)                                 :: E_IO, i

        integer(ip)                                 :: ndim
        integer(ip)                                 :: ielem, nelem
        integer(ip)                                 :: inode, nnode, nnodes
        real(rp),      pointer                      :: coord(:,:)
        integer(ip),   pointer                      :: pnods(:), lnods(:), conn(:)
        integer(ip), allocatable                    :: offset(:)
        integer(1),  allocatable                    :: celltypes(:)
        type(std_vector_integer_ip_t)               :: connectivities
    !-----------------------------------------------------------------
        ! Create dir_path directory if needed
        E_IO = 0
        E_IO = create_directory(dir_path, 0)
        assert(E_IO == 0)

        filename = dir_path//'/'//prefix//'.post.vtu'

        ndim =  this%mesh%get_num_dims()
        nnodes = this%mesh%get_num_vertices()
        nelem = this%mesh%get_num_cells()

        ! Write VTU
        E_IO = VTK_INI_XML(output_format = this%vtk_format,    &
                           filename = filename,                &
                           mesh_topology = 'UnstructuredGrid', &
                           cf=file_id)
        assert(E_IO == 0)

        ! Write coordinates
        coord => this%mesh%get_vertex_coordinates()
        E_IO = VTK_GEO_XML(NN = nnodes, &
                           NC = nelem, &
                           X  = coord(1,:),                  &
                           Y  = coord(2,:),                  &
                           Z  = coord(3,:),                  &
                           cf = file_id)
        assert(E_IO == 0)

        lnods => this%mesh%get_vertices_x_cell()
        pnods => this%mesh%get_vertices_x_cell_pointers()

        call memalloc(nelem+1,   offset,       __FILE__, __LINE__, valin=0, lb1=0)
        call memalloc(nelem,   celltypes,      __FILE__, __LINE__)

        ! Build connectivities from mesh
        do ielem=1, nelem
            nnode = pnods(ielem+1)-pnods(ielem)
            celltypes(ielem) = nnodes_to_vtk_celltype(nnode, ndim)
            offset(ielem) = offset(ielem-1) + nnode
            do inode=1, nnode
                call connectivities%push_back(lnods(pnods(ielem)+inode-1)-1)
            enddo
        end do

        call connectivities%shrink_to_fit()
        conn => connectivities%get_pointer()

        ! Write connectivities
        E_IO = VTK_CON_XML(NC        = nelem,      &
                           connect   = conn,       &
                           offset    = offset(1:), &
                           cell_type = cellTypes,  &
                           cf        = file_id)
        assert(E_IO == 0)


        call memfree(offset,         __FILE__, __LINE__)
        call memfree(celltypes,      __FILE__, __LINE__)
        call connectivities%free()


        ! Write colored partition cell field
        E_IO = VTK_DAT_XML(var_location='Cell',var_block_action='OPEN', cf=file_id)
        assert(E_IO == 0)
        E_IO = VTK_VAR_XML(NC_NN=nelem, varname='partitions', var=this%cells_part(1)%a, cf=file_id)
        assert(E_IO == 0)
        E_IO = VTK_DAT_XML(var_location='Cell', var_block_action='CLOSE', cf=file_id)
        assert(E_IO == 0)

        E_IO = VTK_GEO_XML(cf=file_id)
        assert(E_IO == 0)
        E_IO = VTK_END_XML(cf=file_id)
        assert(E_IO == 0)
    end subroutine vtk_output_handler_write_vtu

  !=============================================================================
  subroutine write_mesh_partition_vtk_format_pl(this, parameter_list)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    type(ParameterList_t)    , intent(in)    :: parameter_list 
    character(len=:), allocatable :: dir_path
    character(len=:), allocatable :: prefix
    call this%get_dir_path_and_prefix_from_pl(parameter_list, dir_path, prefix)
    call this%write_mesh_partition_vtk_format_dir_path_prefix(dir_path,prefix)
  end subroutine write_mesh_partition_vtk_format_pl
  
  !=============================================================================
  subroutine write_mesh_partition_vtk_format_dir_path_prefix(this, dir_path, prefix)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    character(len=*)         , intent(in)    :: dir_path
    character(len=*)         , intent(in)    :: prefix  
    ! Locals
    character(len=:), allocatable :: file_path
    assert ( associated(this%mesh) )
    call this%vtk_output_handler_write_vtu(dir_path, prefix,1)
  end subroutine write_mesh_partition_vtk_format_dir_path_prefix
  
  !=============================================================================
  subroutine mesh_partitioner_get_dir_path_and_prefix_from_pl( parameter_list, dir_path, prefix ) 
     implicit none
     type(ParameterList_t),         intent(in)    :: parameter_list
     character(len=:), allocatable, intent(inout) :: dir_path
     character(len=:), allocatable, intent(inout) :: prefix
     ! Locals
     integer(ip)                                  :: error

     ! Mandatory parameters
     assert(parameter_list%isAssignable(mesh_partitioner_dir_path_key, 'string'))
     error = parameter_list%GetAsString(key = mesh_partitioner_dir_path_key, string = dir_path)
     assert(error==0)

     assert(parameter_list%isAssignable(mesh_partitioner_prefix_key, 'string'))
     error = parameter_list%GetAsString(key = mesh_partitioner_prefix_key, string = prefix)
     assert(error==0)
  end subroutine mesh_partitioner_get_dir_path_and_prefix_from_pl
  
  
  !=============================================================================
  subroutine mesh_partitioner_free(this)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    if ( allocated(this%num_parts_x_level) ) call memfree(this%num_parts_x_level,__FILE__,__LINE__)
    call this%set_default_parameter_values()
    call this%graph_hierarchy_free()
    call this%cells_part_free()
    call this%parts_mapping_free()
    call this%num_parts_x_level_free()
    call this%part_meshes_free()
    call this%mesh_distributions_free()
    nullify(this%mesh)
  end subroutine mesh_partitioner_free
  
    subroutine graph_hierarchy_free ( this )
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    integer(ip) :: i, istat
    if ( allocated(this%graph_hierarchy) ) then
       do i = 1, size(this%graph_hierarchy)
         call this%graph_hierarchy(i)%free()
       end do
       deallocate(this%graph_hierarchy, stat=istat); check(istat==0);
    end if
  end subroutine graph_hierarchy_free 
  
  subroutine part_meshes_free(this)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    integer(ip) :: i, istat
    if ( allocated(this%part_meshes) ) then
       do i=1,size(this%part_meshes)
         call this%part_meshes(i)%free()
       end do
       deallocate(this%part_meshes, stat=istat); check(istat==0);
    end if 
  end subroutine part_meshes_free
  
  subroutine mesh_distributions_free(this)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    integer(ip) :: i, istat
    if ( allocated(this%mesh_distributions) ) then
       do i=1,size(this%mesh_distributions)
         call this%mesh_distributions(i)%free()
       end do
       deallocate(this%mesh_distributions, stat=istat); check(istat==0);
    end if
  end subroutine mesh_distributions_free
  
  subroutine cells_part_free( this )
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    integer(ip) :: i, istat
    if ( allocated(this%cells_part) ) then
       do i = 1, size(this%cells_part)
         call this%cells_part(i)%free()
       end do
       deallocate(this%cells_part, stat=istat); check(istat==0) 
    end if
  end subroutine cells_part_free
  
  subroutine parts_mapping_free(this)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    integer(ip) :: i, istat
    if ( allocated(this%parts_mapping) ) then
      do i=1, size(this%parts_mapping)
        call this%parts_mapping(i)%free()
      end do
      deallocate(this%parts_mapping, stat=istat); check(istat==0) 
    end if  
  end subroutine parts_mapping_free
 
  subroutine num_parts_x_level_free(this)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    if ( allocated(this%num_parts_x_level) ) & 
      call memfree(this%num_parts_x_level,__FILE__,__LINE__)
  end subroutine num_parts_x_level_free
  
  !================================================================================================
  subroutine setup_part_l2g_cells_and_vertices( this, &
                                                ipart, &
                                                num_local_vertices, &
                                                num_local_cells, &
                                                l2g_vertices, &
                                                l2g_cells)
    implicit none
    class(mesh_partitioner_t)  , intent(in)    :: this 
    integer(ip)                , intent(in)    :: ipart
    integer(ip)                , intent(inout) :: num_local_vertices
    integer(ip)                , intent(inout) :: num_local_cells
    integer(igp), allocatable  , intent(inout) :: l2g_vertices(:)
    integer(igp), allocatable  , intent(inout) :: l2g_cells(:)
    
    integer(igp), allocatable :: work1(:)
    integer(igp), allocatable :: work2(:)
    
    integer(ip), pointer :: pnods(:)
    integer(ip), pointer :: lnods(:)
    integer(ip) :: ielem, inode
    
    if (allocated(l2g_vertices)) & 
         call memfree(l2g_vertices, __FILE__, __LINE__ )
    if (allocated(l2g_cells))& 
         call memfree(l2g_cells, __FILE__, __LINE__ )  

    num_local_cells=0
    do ielem=1,this%mesh%get_num_cells()
       if ( this%cells_part(1)%a(ielem) == ipart ) & 
          num_local_cells = num_local_cells + 1 
    end do
    
    call memalloc(num_local_cells, l2g_cells, __FILE__, __LINE__)
    num_local_cells = 0
    do ielem=1,this%mesh%get_num_cells()
       if ( this%cells_part(1)%a(ielem) == ipart ) then 
          num_local_cells = num_local_cells + 1 
          l2g_cells(num_local_cells) = ielem
       end if   
    end do
   
    call memalloc ( this%mesh%get_num_vertices(), work1,__FILE__,__LINE__)
    call memalloc ( this%mesh%get_num_vertices(), work2,__FILE__,__LINE__)

    lnods => this%mesh%get_vertices_x_cell()
    pnods => this%mesh%get_vertices_x_cell_pointers()
    assert(associated(lnods) .and. associated(pnods))

    num_local_vertices=0
    work1 = 0
    work2 = 0
    do ielem=1,this%mesh%get_num_cells()
      if(this%cells_part(1)%a(ielem)==ipart) then
         do inode = pnods(ielem), pnods(ielem+1) - 1 
           if(work1(lnods(inode)) == 0 ) then
              num_local_vertices = num_local_vertices+1
              work1(lnods(inode)) = 1
              work2(num_local_vertices) = lnods(inode)
           end if
         end do
      end if
    end do
    call memalloc(num_local_vertices, l2g_vertices, __FILE__, __LINE__)
    l2g_vertices = work2(1:num_local_vertices)
    call memfree ( work1,__FILE__,__LINE__)
    call memfree ( work2,__FILE__,__LINE__)
  end subroutine setup_part_l2g_cells_and_vertices
  
  !================================================================================================
  subroutine setup_part_interface (this, &
                                   ipart, &
                                   l2g_vertices, &
                                   l2g_cells, &
                                   nebou, &
                                   nnbou, &
                                   lebou, &
                                   lnbou, &
                                   pextn, &
                                   lextn, &
                                   lextp)
    implicit none
    class(mesh_partitioner_t)  , intent(in)    :: this 
    integer(ip)                , intent(in)    :: ipart
    integer(igp)               , intent(in)    :: l2g_vertices(this%part_meshes(ipart)%get_num_vertices())
    integer(igp)               , intent(in)    :: l2g_cells(this%part_meshes(ipart)%get_num_cells())
    integer(ip)                , intent(out)   :: nebou
    integer(ip)                , intent(out)   :: nnbou
    integer(ip)   , allocatable, intent(inout) :: lebou(:)    ! List of boundary elements
    integer(ip)   , allocatable, intent(inout) :: lnbou(:)    ! List of boundary vertices
    integer(ip)   , allocatable, intent(inout) :: pextn(:)    ! Pointers to the lextn
    integer(igp)  , allocatable, intent(inout) :: lextn(:)    ! List of (GID of) external neighbors
    integer(ip)   , allocatable, intent(inout) :: lextp(:)    ! List of parts of external neighbors

    integer(ip) :: lelem, ielem, jelem, pelem, pnode, inode1, inode2, ipoin, lpoin, jpart, iebou, istat, touch
    integer(ip) :: nextn, nexte, nepos
    integer(ip), allocatable :: local_visited(:)
    integer(ip), pointer :: pnods(:), lnods(:)
    integer(ip), allocatable :: pelpo(:), lelpo(:)
    type(hash_table_ip_ip_t)   :: external_visited
    
    if (allocated(lebou)) call memfree(lebou, __FILE__, __LINE__)
    if (allocated(lnbou)) call memfree(lnbou, __FILE__, __LINE__)
    if (allocated(pextn)) call memfree(pextn, __FILE__, __LINE__)
    if (allocated(lextn)) call memfree(lextn, __FILE__, __LINE__)
    if (allocated(lextp)) call memfree(lextp, __FILE__, __LINE__)

    call this%mesh%build_dual_mesh(pelpo, lelpo)
    
    lnods => this%mesh%get_vertices_x_cell()
    pnods => this%mesh%get_vertices_x_cell_pointers()
    assert(associated(lnods) .and. associated(pnods))

    ! Count boundary vertices
    nnbou = 0 
    do lpoin=1, this%part_meshes(ipart)%get_num_vertices()
       ipoin = l2g_vertices(lpoin)
       do pelem = pelpo(ipoin), pelpo(ipoin+1) - 1
          ielem = lelpo(pelem)
          jpart = this%cells_part(1)%a(ielem)
          if ( jpart /= ipart ) then 
             nnbou = nnbou +1
             exit
          end if
       end do
    end do

    ! List boundary vertices
    call memalloc ( nnbou, lnbou, __FILE__, __LINE__ ) 
    nnbou = 0
    do lpoin=1, this%part_meshes(ipart)%get_num_vertices()
       ipoin = l2g_vertices(lpoin)
       do pelem = pelpo(ipoin), pelpo(ipoin+1) - 1
          ielem = lelpo(pelem)
          jpart = this%cells_part(1)%a(ielem)
          if ( jpart /= ipart ) then 
             lnbou(nnbou+1) = ipoin
             nnbou = nnbou +1
             exit
          end if
       end do
    end do

    ! As the dual mesh is given with global IDs we need a hash table to do the touch.
    call memalloc(this%part_meshes(ipart)%get_num_cells(), local_visited,__FILE__,__LINE__)
    local_visited = 0
    call external_visited%init(20)

    ! 1) Count boundary elements and external edges
    touch = 1
    nebou = 0 ! number of boundary elements
    nextn = 0 ! number of external edges
    do lelem = 1, this%part_meshes(ipart)%get_num_cells()
       nexte = 0   ! number of external neighbours of this element
       ielem = l2g_cells(lelem)
       inode1 = pnods(ielem)
       inode2 = pnods(ielem+1)-1
       do pnode = inode1, inode2
          ipoin = lnods(pnode)
          do pelem = pelpo(ipoin), pelpo(ipoin+1) - 1
             jelem = lelpo(pelem)
             if(jelem/=ielem) then
                jpart = this%cells_part(1)%a(jelem)
                if(jpart/=ipart) then                                   ! This is an external element
                   if(local_visited(lelem) == 0 ) nebou = nebou +1        ! Count it
                   !call external_visited%put(key=jelem,val=1, stat=istat) ! Touch jelem as external neighbor of lelem.
                   call external_visited%put(key=jelem,val=touch, stat=istat) ! Touch jelem as external neighbor of lelem.
                   if(istat==now_stored) nexte = nexte + 1                ! Count external neighbours of lelem
                   local_visited(lelem) = nexte                           ! Touch lelem also storing the number
                end if                                                    ! of external neighbours it has
             end if
          end do
       end do
       nextn = nextn + nexte
       ! Clean hash table
       if(local_visited(lelem) /= 0 ) then 
          do pnode = inode1, inode2
             ipoin = lnods(pnode)
             do pelem = pelpo(ipoin), pelpo(ipoin+1) - 1
                jelem = lelpo(pelem)
                if(jelem/=ielem) then
                   jpart = this%cells_part(1)%a(jelem)
                   if(jpart/=ipart) then
                      call external_visited%del(key=jelem, stat=istat)
                   end if
                end if
             end do
          end do
       end if
       call external_visited%print
    end do

    if(ipart==0) then
       write(*,*)  'Visited (boundary) elements:'
       do lelem=1,this%part_meshes(ipart)%get_num_cells()
          write(*,*)  local_visited(lelem)
       end do
    end if

    ! 2) Allocate arrays and store list and pointers to externals
    call memalloc(nebou  , lebou,__FILE__,__LINE__)
    call memalloc(nebou+1, pextn,__FILE__,__LINE__)
    call memalloc(nextn  , lextn,__FILE__,__LINE__)
    call memalloc(nextn  , lextp,__FILE__,__LINE__)

    iebou = 0
    pextn(1) = 1
    do lelem = 1, this%part_meshes(ipart)%get_num_cells()
       if(local_visited(lelem) /= 0 ) then
          iebou = iebou +1
          lebou(iebou) = lelem
          pextn(iebou+1) = local_visited(lelem) + pextn(iebou)
       end if
    end do

    if(ipart==0) then
       write(*,*)  'Boundary elements:'
       do iebou=1,nebou
          write(*,*)  lebou(iebou)
       end do
    end if

    ! 3) Store boundary elements and external edges
    do iebou = 1, nebou
       lelem = lebou(iebou)
       ielem = l2g_cells(lelem)
       nexte = 0   ! number of external neighbours of this element
       inode1 = pnods(ielem)
       inode2 = pnods(ielem+1)-1
       do pnode = inode1, inode2
          ipoin = lnods(pnode)
          do pelem = pelpo(ipoin), pelpo(ipoin+1) - 1
             jelem = lelpo(pelem)
             if(jelem/=ielem) then
                jpart = this%cells_part(1)%a(jelem)
                if(jpart/=ipart) then                                   ! This is an external element
                   call external_visited%put(key=jelem,val=touch, stat=istat) ! Touch jelem as external neighbor of lelem.
                   if(istat==now_stored) then
                      lextn(pextn(iebou)+nexte) = jelem
                      lextp(pextn(iebou)+nexte) = jpart
                      nexte = nexte + 1
                   end if
                end if
             end if
          end do
       end do
       ! Clean hash table
       do pnode = inode1, inode2
          ipoin = lnods(pnode)
          do pelem = pelpo(ipoin), pelpo(ipoin+1) - 1
             jelem = lelpo(pelem)
             if(jelem/=ielem) then
                jpart = this%cells_part(1)%a(jelem)
                if(jpart/=ipart) then                                   ! This is an external element
                   call external_visited%del(key=jelem, stat=istat)
                end if
             end if
          end do
       end do
    end do
    call external_visited%free()
    call memfree(local_visited,__FILE__,__LINE__)
    call memfree(pelpo)
    call memfree(lelpo)
  end subroutine setup_part_interface

  ! Inspired on http://en.wikipedia.org/wiki/Breadth-first_search.
  ! Given a mesh (m) and its dual graph (g), it computes the list 
  ! of vertices (lconn) of each connected component in m. Can be very
  ! useful as a tool to determine whether the mesh partitioning process
  ! leads to disconnected parts or not.
  subroutine mesh_graph_compute_connected_components (m, g, lconn)
    implicit none

    ! Parameters
    type(mesh_t) , intent(in)   :: m   
    type(list_t),  intent(in)   :: g
    type(list_t),  intent(out)  :: lconn

    ! Locals
    integer(ip), allocatable :: auxv(:), auxe(:), e(:)
    integer(ip), allocatable :: emarked(:), vmarked(:)
    integer(ip), allocatable :: q(:)
    integer(ip)              :: head, tail, i, esize, vsize, current, & 
         j, l, k, inods1d, inods2d, p_ipoin, ipoin, graph_num_rows, lconnn
    type(list_iterator_t)    :: graph_column_iterator
    type(list_iterator_t)    :: lconn_iterator
    integer(ip), pointer :: pnods(:), lnods(:)

    graph_num_rows = g%get_num_pointers()
    call memalloc ( graph_num_rows   , auxe     , __FILE__,__LINE__)
    call memalloc ( graph_num_rows   , auxv     , __FILE__,__LINE__)
    call memalloc ( graph_num_rows   , q        , __FILE__,__LINE__)
    call memalloc ( graph_num_rows   , emarked  , __FILE__,__LINE__)
    call memalloc ( m%get_num_vertices(), vmarked  , __FILE__,__LINE__)
    call memalloc ( graph_num_rows   ,  e       , __FILE__,__LINE__)

    lconnn   = 0
    emarked  = 0
    current  = 1
    lnods   => m%get_vertices_x_cell()
    pnods   => m%get_vertices_x_cell_pointers()
    assert(associated(lnods) .and. associated(pnods)) 

    do i=1, graph_num_rows
       if (emarked(i) == 0) then
          ! New connected component
          lconnn = lconnn +1
          esize   = 0
          vsize   = 0
          vmarked = 0 
!!$1  procedure BFS(G,v):
!!$2      create a queue Q
          head=1
          tail=1
!!$3      enqueue v onto Q
          q(tail)=i
          tail=tail+1
!!$4      mark v
          emarked(i)=1
          e(current)=i
          esize  = esize + 1
          current = current + 1  

!!$5      while Q is not empty:
          do while (head/=tail)
!!$6         t ← Q.dequeue()
             j=q(head)
             head = head + 1

             ! Traverse the vertices of the element number j
             inods1d = pnods(j)
             inods2d = pnods(j+1)-1

             do p_ipoin = inods1d, inods2d
                ipoin = lnods(p_ipoin)
                if (vmarked(ipoin)==0) then
                   vmarked(ipoin)=1
                   vsize = vsize+1
                end if
             end do

!!$9         for all edges e in G.adjacentEdges(t) do
             graph_column_iterator = g%create_iterator(j)
             do while(.not. graph_column_iterator%is_upper_bound())
!!$12           u ← G.adjacentVertex(t,e)
                l=graph_column_iterator%get_current()
!!$13           if u is not emarked:
                if (emarked(l)==0) then
!!$14              mark u
                   emarked(l)=1
                   e(current)=l
                   esize  = esize + 1
                   current = current + 1  

!!$15              enqueue u onto Q
                   q(tail)=l
                   tail=tail+1
                end if
                call graph_column_iterator%next()
             end do
          end do
          auxe(lconnn) = esize
          auxv(lconnn) = vsize
       end if
    end do

    call lconn%create(lconnn)

    do i=1, lconnn
       call lconn%sum_to_pointer_index(i, auxv(i))
    end do

    call memfree( auxv   ,__FILE__,__LINE__)
    call memfree( q      ,__FILE__,__LINE__)
    call memfree( emarked,__FILE__,__LINE__)

    call lconn%calculate_header()
    call lconn%allocate_list_from_pointer()

    current=1
    lconn_iterator = lconn%create_iterator()
    do i=1, lconn_iterator%get_size()
       vmarked = 0
       ! Traverse elements of current connected component  
       do current=current,current+auxe(i)-1
          j=e(current)

          ! Traverse the vertices of the element number j
          inods1d = pnods(j)
          inods2d = pnods(j+1)-1

          do p_ipoin = inods1d, inods2d
             ipoin = lnods(p_ipoin)
             if (vmarked(ipoin)==0) then
                vmarked(ipoin)=1
                call lconn_iterator%set_current(ipoin)
                call lconn_iterator%next()
             end if
          end do

       end do
    end do

    call memfree( auxe,__FILE__,__LINE__)
    call memfree( e,__FILE__,__LINE__)
    call memfree( vmarked,__FILE__,__LINE__)
  end subroutine mesh_graph_compute_connected_components

  !=================================================================================================
  subroutine partition_graph(this, &
                                  nparts, &
                                  graph, &
                                  cell_parts, &
                                  cell_weights)
    !-----------------------------------------------------------------------
    ! This routine computes a nparts-way-partitioning of the input graph gp
    !-----------------------------------------------------------------------
    implicit none
    class(mesh_partitioner_t)       , target, intent(in)    :: this
    integer(ip)                     , target, intent(in)    :: nparts
    type(list_t)                    , target, intent(in)    :: graph
    integer(ip)                     , target, intent(inout) :: cell_parts(graph%get_num_pointers())
    integer(ip),            optional, target, intent(in)    :: cell_weights(graph%get_num_pointers())

    ! Local variables 
    integer(ip), target      :: kedge
    integer(ip)              :: idumm,iv
    integer(ip), allocatable :: lwork(:)
    integer(ip)              :: i, j, m, k, ipart
    integer(ip), allocatable :: iperm(:)
#ifdef ENABLE_METIS
    integer(c_int),target :: options(0:METIS_NOPTIONS-1)
    integer(c_int),target :: ncon 
    integer(c_int)        :: ierr
#endif    
   
#ifdef ENABLE_METIS
    ierr = metis_setdefaultoptions(c_loc(options))
    assert(ierr == METIS_OK) 

!!$      From METIS 5.0 manual:
!!$
!!$      The following options are valid for METIS PartGraphRecursive:
!!$      
!!$      METIS_OPTION_CTYPE, METIS_OPTION_IPTYPE, METIS_OPTION_RTYPE,
!!$      METIS_OPTION_NO2HOP, METIS_OPTION_NCUTS, METIS_OPTION_NITER,
!!$      METIS_OPTION_SEED, METIS_OPTION_UFACTOR, METIS_OPTION_NUMBERING,
!!$      METIS_OPTION_DBGLVL
!!$     
!!$      The following options are valid for METIS PartGraphKway:
!!$ 
!!$      METIS_OPTION_OBJTYPE, METIS_OPTION_CTYPE, METIS_OPTION_IPTYPE,
!!$      METIS_OPTION_RTYPE, METIS_OPTION_NO2HOP, METIS_OPTION_NCUTS,
!!$      METIS_OPTION_NITER, METIS_OPTION_UFACTOR, METIS_OPTION_MINCONN,
!!$      METIS_OPTION_CONTIG, METIS_OPTION_SEED, METIS_OPTION_NUMBERING,
!!$      METIS_OPTION_DBGLVL

    if ( this%strat == metis_part_kway ) then
       options(METIS_OPTION_NUMBERING) = 1
       options(METIS_OPTION_DBGLVL)    = this%metis_option_debug
       
       ! Enforce contiguous partititions
       options(METIS_OPTION_CONTIG)    = this%metis_option_contig
       
       ! Explicitly minimize the maximum degree of the part graph
       options(METIS_OPTION_MINCONN)   = this%metis_option_minconn
       if ( this%metis_option_ufactor /= -1 ) options(METIS_OPTION_UFACTOR) = this%metis_option_ufactor

       ! Select random (default) or sorted heavy edge matching
       options(METIS_OPTION_CTYPE)     = this%metis_option_ctype
       options(METIS_OPTION_IPTYPE)    = this%metis_option_iptype
       
       ncon = 1
       if(present(cell_weights)) then
          options(METIS_OPTION_NITER) = 100
          ierr = metis_partgraphkway( graph%get_num_pointers_c_loc(), c_loc(ncon), graph%get_pointers_c_loc(), graph%get_list_c_loc() , & 
                                      C_NULL_PTR  , C_NULL_PTR , c_loc(cell_weights) , c_loc(nparts), &
                                      C_NULL_PTR  , C_NULL_PTR , c_loc(options), c_loc(kedge), c_loc(cell_parts) )
       else
          ierr = metis_partgraphkway( graph%get_num_pointers_c_loc(), c_loc(ncon), graph%get_pointers_c_loc(), graph%get_list_c_loc() , & 
                                      C_NULL_PTR  , C_NULL_PTR , C_NULL_PTR    , c_loc(nparts), &
                                      C_NULL_PTR  , C_NULL_PTR , c_loc(options), c_loc(kedge), c_loc(cell_parts) )
       end if

       assert(ierr == METIS_OK) 
       
    else if ( this%strat == metis_part_recursive ) then
       options(METIS_OPTION_NUMBERING) = 1
       options(METIS_OPTION_DBGLVL)    = this%metis_option_debug
       if ( this%metis_option_ufactor /= -1 ) options(METIS_OPTION_UFACTOR) = this%metis_option_ufactor

       ncon = 1 
       ierr = metis_partgraphrecursive( graph%get_num_pointers_c_loc(), c_loc(ncon), graph%get_pointers_c_loc(), graph%get_list_c_loc() , & 
                                        C_NULL_PTR  , C_NULL_PTR , C_NULL_PTR    , c_loc(nparts), &
                                        C_NULL_PTR  , C_NULL_PTR , c_loc(options), c_loc(kedge), c_loc(cell_parts) )
    end if    
#else
    call enable_metis_error_message
#endif

    if ( this%strat == part_strip ) then
       j = graph%get_num_pointers()
       m = 0
       do ipart=1,nparts
          k = j / (nparts-ipart+1)
          do i = 1, k
             cell_parts(m+i) = ipart
          end do
          m = m + k
          j = j - k
       end do
    else if ( this%strat == part_rcm_strip ) then
       call memalloc ( graph%get_num_pointers(), iperm, __FILE__,__LINE__ )
       call genrcm ( graph, iperm )
       j = graph%get_num_pointers()
       m = 0
       do ipart=1,nparts
          k = j / (nparts-ipart+1)
          do i = 1, k
             cell_parts(iperm(m+i)) = ipart
          end do
          m = m + k
          j = j - k
       end do
       call memfree ( iperm,__FILE__,__LINE__)
    end if
  end subroutine partition_graph


  !================================================================================================
  subroutine setup_part_mesh(this, &
                             ipart, &
                             num_local_vertices, &
                             num_local_cells, &
                             l2g_vertices, &
                             l2g_cells)
    implicit none
    class(mesh_partitioner_t),  intent(inout)  :: this
    integeR(ip)              ,  intent(in)     :: ipart
    integer(ip)              ,  intent(in)     :: num_local_vertices
    integer(ip)              ,  intent(in)     :: num_local_cells
    integer(igp)             ,  intent(in)     :: l2g_vertices(num_local_vertices)
    integer(igp)             ,  intent(in)     :: l2g_cells(num_local_cells)
    
    
    type(hash_table_igp_ip_t)      :: ws_inmap
    type(hash_table_igp_ip_t)      :: el_inmap
    integer(ip)    , allocatable   :: node_list(:)
    integer(ip)                    :: aux, ipoin,inode,knode,kvef_size,lvef_size,istat
    integer(ip)                    :: ielem_lmesh,ielem_gmesh,ivef_lmesh,ivef_gmesh
    integer(ip)                    :: p_ipoin_lmesh,p_ipoin_gmesh
    type(list_iterator_t)          :: given_vefs_iterator
    logical :: count_it
    integer(ip), pointer  :: gpnods(:), glnods(:), glegeo(:), gleset(:), glst_vefs_geo(:), glst_vefs_set(:)
    type(list_t), pointer :: ggiven_vefs
    real(rp), pointer     :: gcoord(:,:)
    
    ! Sizes
    integer(ip) :: order 
    integer(ip) :: nelty
    integer(ip) :: ndime
    integer(ip) :: npoin
    integer(ip) :: nelem
    
    ! Elements
    integer(ip), allocatable :: pnods(:)
    integer(ip), allocatable :: lnods(:)
    integer(ip), allocatable :: legeo(:)
    integer(ip), allocatable :: leset(:)
      
    ! Boundary vefs-related data
    type(list_t)             :: lgiven_vefs
    integer(ip), allocatable :: lst_vefs_geo(:)
    integer(ip), allocatable :: lst_vefs_set(:)    
      
    ! Vertex coordinates
    real(rp),  allocatable   :: coord(:,:)
    
    integer(ip) :: nnodb, nnode

    order=this%mesh%get_element_order()
    nelty=this%mesh%get_num_element_types()
    ndime=this%mesh%get_num_dims()
    npoin=num_local_vertices
    nelem=num_local_cells

    call ws_inmap%init(max(int(num_local_vertices*0.25,ip),10))
    do ipoin=1,num_local_vertices
       ! aux is used to avoid compiler warning related to val being an intent(inout) argument
       aux = ipoin
       call ws_inmap%put(key=l2g_vertices(ipoin),val=aux,stat=istat) 
    end do

    call el_inmap%init(max(int(num_local_cells*0.25,ip),10))
    do ipoin=1,num_local_cells
       ! aux is used to avoid compiler warning related to val being an intent(inout) argument
       aux = ipoin
       call el_inmap%put(key=l2g_cells(ipoin),val=aux,stat=istat) 
    end do

    ! Elements
    call memalloc(nelem+1, pnods, __FILE__,__LINE__)
    call memalloc(nelem  , legeo, __FILE__,__LINE__)
    call memalloc(nelem  , leset, __FILE__,__LINE__)
    nnode=0
    pnods=0
    pnods(1)=1

    glnods   => this%mesh%get_vertices_x_cell()
    gpnods   => this%mesh%get_vertices_x_cell_pointers()
    assert(associated(glnods) .and. associated(gpnods))
    glegeo  => this%mesh%get_cells_geometry_id()
    gleset => this%mesh%get_cells_set()

    do ielem_lmesh=1,nelem
       ielem_gmesh = l2g_cells(ielem_lmesh)
       knode = gpnods(ielem_gmesh+1)-gpnods(ielem_gmesh)
       pnods(ielem_lmesh+1)=pnods(ielem_lmesh)+knode
       nnode=max(nnode,knode)
       legeo(ielem_lmesh)=glegeo(ielem_gmesh)
       leset(ielem_lmesh)=gleset(ielem_gmesh)
    end do
    call memalloc (pnods(nelem+1)-1, lnods, __FILE__,__LINE__)
    do ielem_lmesh=1,nelem
       ielem_gmesh = l2g_cells(ielem_lmesh)
       p_ipoin_gmesh = gpnods(ielem_gmesh)-1
       p_ipoin_lmesh = pnods(ielem_lmesh)-1
       knode = gpnods(ielem_gmesh+1)-gpnods(ielem_gmesh)
       do inode=1,knode
          call ws_inmap%get(key=int(glnods(p_ipoin_gmesh+inode),igp),val=lnods(p_ipoin_lmesh+inode),stat=istat) 
       end do
    end do

    ! Boundary elements
    ivef_lmesh=0
    nnodb=0
    lvef_size=0
    ggiven_vefs => this%mesh%get_boundary_vefs()
    assert(associated(ggiven_vefs))
    do ivef_gmesh=1,ggiven_vefs%get_num_pointers()
       given_vefs_iterator = ggiven_vefs%create_iterator(ivef_gmesh)
       kvef_size = given_vefs_iterator%get_size()
       count_it=.true.
       do while(.not. given_vefs_iterator%is_upper_bound())
          call ws_inmap%get(key=int(given_vefs_iterator%get_current(),igp),val=knode,stat=istat)
          call given_vefs_iterator%next()
          if(istat==key_not_found) then
             count_it=.false.
             exit
          end if
       end do
       if(count_it) then
          lvef_size=lvef_size+kvef_size
          nnodb=max(nnodb,kvef_size)
          ivef_lmesh=ivef_lmesh+1
       end if
    end do

    if(ivef_lmesh>0) then

       call memalloc (  nnodb,   node_list, __FILE__,__LINE__)
       call memalloc(   ivef_lmesh, lst_vefs_geo, __FILE__,__LINE__)
       call memalloc(   ivef_lmesh, lst_vefs_set, __FILE__,__LINE__)

       glst_vefs_geo => this%mesh%get_boundary_vefs_geometry_id()
       glst_vefs_set => this%mesh%get_boundary_vefs_set()
       assert(associated(glst_vefs_geo) .and. associated(glst_vefs_set))

       call lgiven_vefs%create(ivef_lmesh)

       ivef_lmesh=1
       do ivef_gmesh=1,ggiven_vefs%get_num_pointers()
          given_vefs_iterator = ggiven_vefs%create_iterator(ivef_gmesh)
          kvef_size = given_vefs_iterator%get_size()
          count_it=.true.
          do inode=1,kvef_size
             call ws_inmap%get(key=int(given_vefs_iterator%get_current(),igp),val=node_list(inode),stat=istat)
             call given_vefs_iterator%next()
             if(istat==key_not_found) then
                count_it=.false.
                exit
             end if
          end do
          if(count_it) then
             call lgiven_vefs%sum_to_pointer_index(ivef_lmesh, kvef_size)
             lst_vefs_geo(ivef_lmesh)=glst_vefs_geo(ivef_gmesh)
             lst_vefs_set(ivef_lmesh)=glst_vefs_set(ivef_gmesh)
             ivef_lmesh=ivef_lmesh+1
          end if
       end do

       call lgiven_vefs%calculate_header()
       call lgiven_vefs%allocate_list_from_pointer()

       ivef_lmesh=1
       do ivef_gmesh=1,ggiven_vefs%get_num_pointers()
          given_vefs_iterator = ggiven_vefs%create_iterator(ivef_gmesh)
          kvef_size = given_vefs_iterator%get_size()
          count_it=.true.
          do inode=1,kvef_size
             call ws_inmap%get(key=int(given_vefs_iterator%get_current(),igp),val=node_list(inode),stat=istat)
             call given_vefs_iterator%next()
             if(istat==key_not_found) then
                count_it=.false.
                exit
             end if
          end do
          if(count_it) then
             given_vefs_iterator = lgiven_vefs%create_iterator(ivef_lmesh)
             do inode=1,kvef_size
                call given_vefs_iterator%set_current(node_list(inode))
                call given_vefs_iterator%next()
             enddo
             ivef_lmesh=ivef_lmesh+1
          end if
       end do
       call memfree (node_list, __FILE__,__LINE__)
    end if
    
    call ws_inmap%free()
    call el_inmap%free()

    gcoord => this%mesh%get_vertex_coordinates()
    assert(associated(gcoord))

    call memalloc(SPACE_DIM, npoin, coord, __FILE__,__LINE__)
    do ipoin=1,num_local_vertices
       coord(:,ipoin)=gcoord(:,l2g_vertices(ipoin))
    end do
    
    call this%part_meshes(ipart)%create(order, &
                           nelty, &
                           ndime, &
                           npoin, &
                           nelem, &
                           pnods, &
                           lnods, &
                           legeo, &
                           leset, &
                           lgiven_vefs, &
                           lst_vefs_geo, &
                           lst_vefs_set, &
                           coord)
    
    call memfree(pnods, __FILE__, __LINE__)
    call memfree(lnods, __FILE__, __LINE__)
    call memfree(legeo, __FILE__, __LINE__)
    call memfree(leset, __FILE__, __LINE__)
    call memfree(lst_vefs_geo, __FILE__, __LINE__)
    call memfree(lst_vefs_set, __FILE__, __LINE__)
    call memfree(coord, __FILE__, __LINE__)
    call lgiven_vefs%free()
  end subroutine setup_part_mesh

  subroutine build_parts_graph (nparts, cell_parts, graph, parts_graph)
    implicit none
    integer(ip)             , intent(in)    :: nparts
    type(list_t)            , intent(in)    :: graph
    integer(ip)             , intent(in)    :: cell_parts(:)
    type(list_t)            , intent(inout) :: parts_graph

    integer(ip)              :: istat,ielem,jelem,ipart,jpart
    integer(ip)              :: num_parts_around, touched
    type(list_iterator_t)                    :: graph_iterator
    type(list_iterator_t)                    :: parts_graph_iterator
    type(position_hash_table_t), allocatable :: visited_parts_touched(:)
    type(hash_table_ip_ip_t)   , allocatable :: visited_parts_numbers(:)

    call parts_graph%create(nparts)

    ! The maximum number of parts around a part can be estimated from the
    ! maximum number of elements connected to an element, that is, by the 
    ! maximum degree of elements graph. Note, however, that it can be bigger.
    num_parts_around=0
    do ielem=1,graph%get_num_pointers()
       graph_iterator = graph%create_iterator(ielem)
       num_parts_around = max(num_parts_around,graph_iterator%get_size())
    end do
    allocate(visited_parts_touched(nparts),stat=istat); assert(istat==0);
    allocate(visited_parts_numbers(nparts),stat=istat); assert(istat==0);
    do ipart=1,nparts
       call visited_parts_touched(ipart)%init(num_parts_around)
       call visited_parts_numbers(ipart)%init(num_parts_around)
    end do

    ! Now compute graph pointers and fill tables
    do ielem=1,graph%get_num_pointers()
       ipart = cell_parts(ielem)
       graph_iterator = graph%create_iterator(ielem)
       do while(.not.graph_iterator%is_upper_bound())
          jelem = graph_iterator%get_current()
          jpart = cell_parts(jelem)
          call visited_parts_touched(ipart)%get(key=jpart,val=num_parts_around,stat=istat) ! Touch it (jpart is around ipart)
          if(istat==new_index) then
             call visited_parts_numbers(ipart)%put(key=num_parts_around,val=jpart,stat=istat) ! Store it
             assert(istat==now_stored)
          end if
          call graph_iterator%next()
       end do
    end do
    do ipart=1,nparts
       call parts_graph%sum_to_pointer_index(ipart,visited_parts_touched(ipart)%last())
    end do

    ! Fill graph from tables
    call parts_graph%calculate_header()
    call parts_graph%allocate_list_from_pointer()
    do ipart=1,nparts
       num_parts_around = 0
       parts_graph_iterator = parts_graph%create_iterator(ipart)
       do while(.not.parts_graph_iterator%is_upper_bound())
          num_parts_around = num_parts_around + 1
          call visited_parts_numbers(ipart)%get(key=num_parts_around,val=jpart,stat=istat) 
          assert(istat==key_found)
          call parts_graph_iterator%set_current(jpart)
          call parts_graph_iterator%next()
       end do
       assert(num_parts_around==visited_parts_touched(ipart)%last())
       call visited_parts_touched(ipart)%free() ! This could be done before, as far as the assert is eliminated
       call visited_parts_numbers(ipart)%free()
    end do

    deallocate(visited_parts_touched,stat=istat)
    deallocate(visited_parts_numbers,stat=istat)
  end subroutine build_parts_graph
  
end module mesh_partitioner_names
