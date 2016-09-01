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
module uniform_hex_mesh_generator_names
  ! Serial modules
  use types_names
  use memor_names
  use hex_boundary_set_ids_descriptor_names
  use reference_fe_names
  use FPL
  implicit none
# include "debug.i90"
  private
  
  ! WE SHOULD TRY TO MINIMIZE THE MODULES/DATA TYPES ON WHICH THIS MODULE DEPENDS.
  ! IN PARTICULAR, IT MUST NOT DEPEND ON THE FEMPAR'S TRIANGULATION MODULES, BUT 
  ! PROVIDE THE DATA AS PLAIN ARRAYS THAT THE TRIANGULATION REQUIRES IN ORDER TO BE
  ! GENERATED

  ! FPL KEY names
  ! character(len=*), parameter :: number_elements_name          = 'number_elements'          
  ! character(len=*), parameter :: number_parts_name             = 'number_parts'             
  ! character(len=*), parameter :: number_sockets_name           = 'number_sockets'           
  ! character(len=*), parameter :: discretization_type_name      = 'discretization_type'      
  ! character(len=*), parameter :: periodic_boundaries_name      = 'periodic_boundaries'      
  ! character(len=*), parameter :: number_elements_boundary_name = 'number_elements_boundary' 
  ! character(len=*), parameter :: material_case_name            = 'material_case'            
  ! character(len=*), parameter :: domain_length_name            = 'domain_length'            
  ! character(len=*), parameter :: origin_name                   = 'origin'                   
  ! character(len=*), parameter :: stretching_parameter_name     = 'stretching_parameter'     
  ! character(len=*), parameter :: size_boundary_name            = 'size_boundary'  

  character(len=*), parameter :: number_of_dimensions_key    = 'number_of_dimensions'
  character(len=*), parameter :: number_of_cells_per_dir_key = 'number_of_cells_per_dir'
  character(len=*), parameter :: number_of_parts_per_dir_key = 'number_of_parts_per_dir'
  character(len=*), parameter :: is_dir_periodic_key         = 'is_dir_periodic'
  character(len=*), parameter :: interpolation_order_key     = 'interpolation_order'

  public :: number_of_dimensions_key
  public :: number_of_cells_per_dir_key 
  public :: number_of_parts_per_dir_key 
  public :: is_dir_periodic_key         
  public :: interpolation_order_key     

  ! Formely uniform_mesh_descriptor_t. Provided here for convenience (so that legacy 
  ! code can be re-used here AS IS). However, the user of FEMPAR should not be aware of this 
  ! data type. Instead, it would be better that he uses FPL to parametrize the creation of 
  ! the triangulation. I would like to avoid having multiple data types for algorithms parameters 
  ! (e.g., one per algorithm/module), using a single one (i.e., FPL) for all of them ...
  type uniform_hex_mesh_descriptor_t 
    private 
    integer(ip) :: number_of_dimensions
    integer(ip) :: interpolation_order
    integer(ip) :: number_of_cells_per_dir(0:SPACE_DIM-1)
    integer(ip) :: number_of_parts_per_dir(0:SPACE_DIM-1)
    integer(ip) :: is_dir_periodic(0:SPACE_DIM-1)

    ! integer(ip)            ::  &
    !       ntdix=0,             &         ! Type of discretization in x (0=uniform, 1=cubic, 2=tanh, 3=imh+unif, 4:imh+tanh)
    !       ntdiy=0,             &         ! Type of discretization in y (0=uniform, 1=cubic, 2=tanh, 3=imh+unif, 4:imh+tanh)
    !       ntdiz=0,             &         ! Type of discretization in z (0=uniform, 1=cubic, 2=tanh, 3=imh+unif, 4:imh+tanh)
    !       mater=0,             &         ! Material case
    !       neblx=0,             &         ! Number of elements in the x boundary layer
    !       nebly=0,             &         ! Number of elements in the y boundary layer
    !       neblz=0,             &         ! Number of elements in the z boundary layer
    !       ndime,               &         ! Number of dimensions
    !       nparts=1,            &         ! Number of partitions
    !       nedir(3),            &         ! Number of elements in each direction
    !       npdir(3)=[1,1,1],    &         ! Number of parts on each direction           
    !       nsckt(3)=[1,1,1],    &         ! Number of parts on each socket and direction
    !       isper(3)=[0,0,0],    &         ! Flag for periodic boundary conditions on each direction
    !       pdegr=1                        ! Order of interpolation
    ! real(rp)               :: &
    !       xleng   = 1.0_rp,    &         ! Size of the domain in x
    !       yleng   = 1.0_rp,    &         ! Size of the domain in y
    !       zleng   = 1.0_rp,    &         ! Size of the domain in z
    !       zx1     = 0.1_rp,    &         ! size of the elements at x=0   (left)  
    !       zx2     = 0.1_rp,    &         ! size of the elements at x=a/2 (center)
    !       zy1     = 0.1_rp,    &         ! size of the elements at y=0   (bottom)
    !       zy2     = 0.1_rp,    &         ! size of the elements at y=b/2 (center)
    !       x0      = 0.0_rp,    &         ! Origin x-coordinate
    !       y0      = 0.0_rp,    &         ! Origin y-coordinate
    !       z0      = 0.0_rp,    &         ! Origin z-coordinate
    !       xstret  = 2.75_rp,   &         ! Stretching parameter
    !       ystret  = 2.75_rp,   &         ! Stretching parameter
    !       zstret  = 2.75_rp,   &         ! Stretching parameter
    !       xlengbl = 0.0_rp,    &         ! Size of the boundary layer in x
    !       ylengbl = 0.0_rp,    &         ! Size of the boundary layer in y
    !       zlengbl = 0.0_rp               ! Size of the boundary layer in z
  contains   
    !procedure, non_overridable          :: describe_mesh                => uniform_hex_mesh_descriptor_describe_mesh
    procedure, non_overridable          :: get_data_from_parameter_list => uniform_hex_mesh_descriptor_get_data_from_parameter_list
    !procedure, non_overridable, private :: generate_geom_and_topo_size  => uniform_hex_mesh_descriptor_generate_geom_and_topo_size
  end type uniform_hex_mesh_descriptor_t
  
  ! ! Private data type (required by helper legacy subroutines)
  ! type geom_size_t 
  !    integer(ip)            :: &
  !         nnode,               &         ! Number of nodes on each element
  !         nedir(3),            &         ! Number of elements in each direction
  !         npdir(3),            &         ! Number of parts on each direction           
  !         nsckt(3),            &         ! Number of parts on each socket and direction
  !         npsoc(3),            &         ! Number of sockets on each direction
  !         nedom(3),            &         ! Number of elements on each part and direction
  !         npdom(3),            &         ! Number of points on each part and direction
  !         ncorn(3),            &         ! Number of partition corners on each direction
  !         nedge(3,3),          &         ! Number of partition edges on each direction
  !         nface(3,3),          &         ! Number of partition faces on each direction
  !         npdomt,              &         ! Number of points on each domain
  !         nedomt,              &         ! Number of elements on each domain
  !         ncornt,              &         ! Number of partition corners on each domain
  !         nedget(3),           &         ! Number of edges on each direction
  !         nfacet(3),           &         ! Number of faces on each direction
  !         nedgett,             &         ! Number of edges on each domain
  !         nfacett,             &         ! Number of faces on each domain
  !         neghost                        ! Maximum of elements on each domain counting ghosts
  ! end type geom_size_t

  ! ! Private data type (required by helper legacy subroutines)
  ! type topo_size_t
  !    integer(ip)            :: &
  !         notot,               &         ! Total amount of elemental objects of the partition
  !         nctot,               &         ! Total amount of elemental corners of the partition
  !         ncglb,               &         ! Total amount of elemental corners of the domain
  !         noglb,               &         ! Total amount of elemental corners of the domain
  !         neglb,               &         ! Total amount of elements of the domain
  !         ndtot,               &         ! Total amount of elemental edges of the partition
  !         nftot,               &         ! Total amount of elemental faces of the partition
  !         nddir(3),            &         ! Total amount of elemental edges of the partition for each direction
  !         nfdir(3),            &         ! Total amount of elemental faces of the partition for each direction
  !         nddom(3,3),          &         ! # edges for each direction given the edge direction (local)
  !         nfdom(3,3),          &         ! # faces for each direction given the face normal direction (local)
  !         ndglb(3,3),          &         ! # edges for each direction given the edge direction (global)
  !         nfglb(3,3),          &         ! # faces for each direction given the face normal direction (global)
  !         ndsum(3),            &         ! Total amount of elemental edges of the domain for each direction
  !         nfsum(3)                       ! Total amount of elemental faces of the domain for each direction
  ! end type topo_size_t
  
  interface ijk_to_spatial_numbering
     module procedure ijk_to_spatial_numbering_ip !, ijk_to_spatial_numbering_igp
  end interface
  
  public :: uniform_hex_mesh_generator_generate_connectivities
  
contains

  subroutine uniform_hex_mesh_descriptor_get_data_from_parameter_list(this,parameter_list)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine generates geometry data to construct a structured mesh                      !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(uniform_hex_mesh_descriptor_t), intent(inout) :: this
    type(ParameterList_t)               , intent(in)    :: parameter_list
    ! Locals
    integer(ip)          :: istat
    logical              :: is_present

    ! Mandatory parameters
    is_present = .true.
    is_present =  is_present.and. parameter_list%isPresent(key = number_of_dimensions_key )
    is_present =  is_present.and. parameter_list%isPresent(key = number_of_cells_per_dir_key )
    is_present =  is_present.and. parameter_list%isPresent(key = is_dir_periodic_key )
    is_present =  is_present.and. parameter_list%isPresent(key = interpolation_order_key )
    !is_present =  is_present.and. parameter_list%isPresent(key = 
    assert(is_present)

    istat = 0
    istat = istat + parameter_list%get(key = number_of_dimensions_key   , value = this%number_of_dimensions)
    istat = istat + parameter_list%get(key = number_of_cells_per_dir_key, value = this%number_of_cells_per_dir)
    istat = istat + parameter_list%get(key = is_dir_periodic_key        , value = this%is_dir_periodic)
    istat = istat + parameter_list%get(key = interpolation_order_key    , value = this%interpolation_order)
    check(istat==0)

    ! Optional parameters
    if( parameter_list%isPresent(key = number_of_parts_per_dir_key) ) then
       istat = parameter_list%get(key = number_of_parts_per_dir_key , value = this%number_of_parts_per_dir)
    else
       this%number_of_parts_per_dir = 1
    end if

  end subroutine uniform_hex_mesh_descriptor_get_data_from_parameter_list

  ! Main driver subroutine of this module
  subroutine uniform_hex_mesh_generator_generate_connectivities(parameter_list,        &
                                                                num_local_cells,       &
                                                                num_local_vefs,        &
                                                                num_vertices,          &
                                                                ptr_vefs_per_cell,     &
                                                                lst_vefs_lids,         &
                                                                boundary_id,           &
                                                                coordinates,           &
                                                                num_ghost_cells,       &
                                                                cells_gids,            &
                                                                vefs_gids,             &
                                                                num_itfc_cells,        &
                                                                lst_itfc_cells,        &
                                                                ptr_ext_neighs_per_itfc_cell, &
                                                                lst_ext_neighs_gids,          &
                                                                lst_ext_neighs_part_ids,      &
                                                                part_id)
    implicit none
    type(ParameterList_t)     , intent(in)    :: parameter_list
    integer(ip)               , intent(out)   :: num_local_cells
    integer(ip)               , intent(out)   :: num_local_vefs
    integer(ip)               , intent(out)   :: num_vertices
    integer(ip)  , allocatable, intent(inout) :: ptr_vefs_per_cell(:)            ! Size = num_local_cells + 1
    integer(ip)  , allocatable, intent(inout) :: lst_vefs_lids(:)                ! Size = ptr_vefs_per_cell(num_local_cells+1)-1
    integer(ip)  , allocatable, intent(inout) :: boundary_id(:)                  ! Size = num_local_vefs
    real(rp)     , allocatable, intent(inout) :: coordinates(:,:)

    integer(ip)               , optional, intent(out)   :: num_ghost_cells
    integer(ip)               , optional, intent(out)   :: num_itfc_cells
    integer(igp) , allocatable, optional, intent(inout) :: cells_gids(:)                   ! Size = num_local_cells 
    integer(igp) , allocatable, optional, intent(inout) :: vefs_gids(:)                    ! Size = num_local_vefs
    integer(ip)  , allocatable, optional, intent(inout) :: lst_itfc_cells(:)              
    integer(ip)  , allocatable, optional, intent(inout) :: ptr_ext_neighs_per_itfc_cell(:)
    integer(igp) , allocatable, optional, intent(inout) :: lst_ext_neighs_gids(:)         
    integer(ip)  , allocatable, optional, intent(inout) :: lst_ext_neighs_part_ids(:)
    integer(ip)               , optional, intent(in)    :: part_id

    type(uniform_hex_mesh_descriptor_t) :: input_data
    integer(ip), allocatable :: cell_is_ghost(:)

    integer(ip) :: part_ijk(0:SPACE_DIM-1)
    integer(ip) :: cell_ijk(0:SPACE_DIM-1)
    integer(ip) :: neighbor_ijk(0:SPACE_DIM-1)
    integer(ip) :: neighbor_part_ijk(0:SPACE_DIM-1)
    integer(ip) :: nface_ijk(0:SPACE_DIM-1)
    integer(ip) :: first_cell_ijk(0:SPACE_DIM-1)
    integer(ip) :: last_cell_ijk(0:SPACE_DIM-1)

    ! Here total=local+ghost (if any) refers to the things I have,
    ! whereas global to the whole distributed mesh.
    integer(ip) :: num_total_cells_per_dir(0:SPACE_DIM-1)
    integer(ip) :: num_local_cells_per_dir(0:SPACE_DIM-1)
    integer(ip) :: num_left_parts_per_dir(0:SPACE_DIM-1)
    integer(ip) :: num_right_parts_per_dir(0:SPACE_DIM-1)

    integer(igp), allocatable  :: num_global_n_faces(:)
    integer(ip) , allocatable  :: num_total_n_faces(:)
    integer(ip) , allocatable  :: num_global_nfaces_per_dir(:,:)
    integer(ip) , allocatable  :: num_total_nfaces_per_dir(:,:)

    integer(ip)               :: topology, num_nface_types, partial_count
    integer(ip)               :: ighost_cell, ilocal_cell
    integer(ip)               :: num_ghost_cells_ , has_left_ghost, has_right_ghost
    integer(ip)               :: idime, jdime, icell, iface, iface_of_itype, index, itype, itfc_cells

    type(polytope_tree_t)     :: polytope_tree
    type(node_array_t)        :: node_array
    integer(ip)               :: ones(SPACE_DIM)
    logical                   :: count_it

    if(present(num_ghost_cells)) then
       assert(present(num_itfc_cells))
       assert(present(cells_gids))
       assert(present(vefs_gids))
       assert(present(lst_itfc_cells))
       assert(present(ptr_ext_neighs_per_itfc_cell))
       assert(present(lst_ext_neighs_gids))  
       assert(present(lst_ext_neighs_part_ids))
    end if

    call input_data%get_data_from_parameter_list(parameter_list)

    ones = 1
    topology = 2**input_data%number_of_dimensions-1  ! Hexahedral
    call polytope_tree%create( input_data%number_of_dimensions, topology )  
    call node_array%create ( polytope_tree, ones*input_data%interpolation_order )

    ! PARTS
    ! =====
    ! Get my part coordinates (make it 0-based, assuming part_id is 1-based) and the number of parts I have around (if any)
    if(present(part_id)) then
       call spatial_to_ijk_numbering(input_data%number_of_dimensions, input_data%number_of_parts_per_dir, part_id, part_ijk)
    else
       part_ijk = 0
    end if
    num_left_parts_per_dir=1
    num_right_parts_per_dir=1
    do idime = 0, input_data%number_of_dimensions - 1 
       if(input_data%is_dir_periodic(idime)==0.or.input_data%number_of_parts_per_dir(idime)==1) then ! Not periodic
          if(part_ijk(idime)==0) num_left_parts_per_dir(idime)=0 
          if(part_ijk(idime)==input_data%number_of_parts_per_dir(idime)-1) num_right_parts_per_dir(idime)=0 
       end if
    end do

    ! CELLS
    ! =====
    ! Global and local number of cells (per direction and total; local, ghost and global)
    do idime = 0, input_data%number_of_dimensions - 1 
       num_local_cells_per_dir(idime) = input_data%number_of_cells_per_dir(idime) / input_data%number_of_parts_per_dir(idime) 
       first_cell_ijk(idime) =  part_ijk(idime) * num_local_cells_per_dir(idime)
       !last_cell_ijk(idime) =  first_cell_ijk(idime) + num_local_cells_per_dir(idime)
    end do
    num_local_cells = 1
    do idime = 0, input_data%number_of_dimensions - 1
       num_local_cells = num_local_cells * num_local_cells_per_dir(idime)
    end do
    num_total_cells_per_dir = num_local_cells_per_dir + num_left_parts_per_dir + num_right_parts_per_dir
    first_cell_ijk = first_cell_ijk - num_left_parts_per_dir

    if(present(num_ghost_cells)) then
       num_ghost_cells_ = 1
       do idime = 0, input_data%number_of_dimensions - 1
          num_ghost_cells_ = num_ghost_cells_ * num_total_cells_per_dir(idime)
       end do
       num_ghost_cells_ = num_ghost_cells_ - num_local_cells
       num_ghost_cells  = num_ghost_cells_
    else
       num_ghost_cells_ = 0
    end if

    ! The following paragraph is not needed but I keep it because the loop can be useful
    ! num_ghost_cells_ = 0
    ! do iface=1,polytope_tree%get_number_n_faces()
    !    if(polytope_tree%get_n_face_dimension(iface)<input_data%number_of_dimensions) then ! do not include the polytope itself
    !       partial_count=1
    !       count_it = .true.
    !       do idime = 0, input_data%number_of_dimensions - 1 
    !          if(polytope_tree%n_face_dir_is_fixed(iface,idime)==1) then
    !             partial_count=partial_count*num_local_cells_per_dir(idime)
    !          else
    !             if( (polytope_tree%n_face_dir_coordinate(iface,idime)==0.and.num_left_parts_per_dir(idime)==0) .or. &
    !                 (polytope_tree%n_face_dir_coordinate(iface,idime)==1.and.num_right_parts_per_dir(idime)==0) ) count_it = .false.
    !          end if
    !       end do
    !       if(count_it) num_ghost_cells_ = num_ghost_cells_ + partial_count
    !    end if
    ! end do

    ! N_FACES
    ! =======
    ! Global and local number of n_faces (per direction and total)
    num_nface_types = 0
    do iface=1,polytope_tree%get_number_n_faces()
       if(polytope_tree%get_n_face_dimension(iface)<input_data%number_of_dimensions.and. &
            & polytope_tree%n_face_coordinate(iface)==0) num_nface_types = num_nface_types + 1
    end do

    call memalloc( num_nface_types+1, num_global_n_faces, __FILE__,__LINE__,lb1=0)
    call memalloc( num_nface_types+1, num_total_n_faces, __FILE__,__LINE__,lb1=0)
    call memalloc( input_data%number_of_dimensions, num_nface_types, num_global_nfaces_per_dir, __FILE__,__LINE__,lb1=0,lb2=0)
    call memalloc( input_data%number_of_dimensions, num_nface_types, num_total_nfaces_per_dir, __FILE__,__LINE__,lb1=0,lb2=0)
    do iface=1,polytope_tree%get_number_n_faces()
       if(polytope_tree%get_n_face_dimension(iface)<input_data%number_of_dimensions.and. &
            & polytope_tree%n_face_coordinate(iface)==0) then
          itype = polytope_tree%n_face_type(iface)
          do idime = 0, input_data%number_of_dimensions - 1
             num_global_nfaces_per_dir(idime,itype) = &
                  & input_data%number_of_cells_per_dir(idime)  + &
                  & 1 - max(polytope_tree%n_face_dir_is_fixed(iface,idime),input_data%is_dir_periodic(idime))
             num_total_nfaces_per_dir(idime,itype) =  &
                  & num_total_cells_per_dir(idime) + &
                  & 1 - max(polytope_tree%n_face_dir_is_fixed(iface,idime),input_data%is_dir_periodic(idime)/input_data%number_of_parts_per_dir(idime))
          end do
          num_global_n_faces(itype+1) = 1
          num_total_n_faces(itype+1) = 1
          do idime = 0, input_data%number_of_dimensions - 1
             num_global_n_faces(itype+1) = num_global_n_faces(itype+1) * num_global_nfaces_per_dir(idime,itype)
             num_total_n_faces(itype+1) = num_total_n_faces(itype+1) * num_total_nfaces_per_dir(idime,itype)
          end do
       end if
    end do
    num_global_n_faces(0) = 1
    num_total_n_faces(0) = 1
    do iface=1,polytope_tree%get_number_n_faces()
       if(polytope_tree%get_n_face_dimension(iface)<input_data%number_of_dimensions.and. &
            & polytope_tree%n_face_coordinate(iface)==0) then
          itype = polytope_tree%n_face_type(iface)
          num_global_n_faces(itype+1) = num_global_n_faces(itype+1) + num_global_n_faces(itype)
          num_total_n_faces(itype+1) = num_total_n_faces(itype+1) + num_total_n_faces(itype)
       end if
    end do
    num_local_vefs = num_total_n_faces(num_nface_types) - 1
    num_vertices = num_total_n_faces(1) - 1

    ! FILL ARRAYS
    ! Construct local numbering
    call memalloc(num_local_cells+num_ghost_cells_+1, ptr_vefs_per_cell, __FILE__,__LINE__)
    ptr_vefs_per_cell = polytope_tree%get_number_n_faces() - 1 ! the cell itself does not count
    ptr_vefs_per_cell(1) = 1
    do icell = 1, num_local_cells+num_ghost_cells_
       ptr_vefs_per_cell(icell+1) = ptr_vefs_per_cell(icell+1) + ptr_vefs_per_cell(icell)
    end do
    call memalloc( ptr_vefs_per_cell(num_local_cells+num_ghost_cells_+1)-1, lst_vefs_lids, __FILE__,__LINE__)
    do icell = 1, num_local_cells+num_ghost_cells_
       call spatial_to_ijk_numbering(input_data%number_of_dimensions, num_total_cells_per_dir, icell, cell_ijk)
       do iface=1,polytope_tree%get_number_n_faces()
          if(polytope_tree%get_n_face_dimension(iface)<input_data%number_of_dimensions) then ! do not include the polytope itself
             do idime = 0, input_data%number_of_dimensions - 1
                nface_ijk(idime) = mod(cell_ijk(idime) + polytope_tree%n_face_dir_coordinate(iface,idime),num_total_nfaces_per_dir(idime,itype))
             end do
             itype = polytope_tree%n_face_type(iface)
             lst_vefs_lids(ptr_vefs_per_cell(icell)+iface-1) = num_total_n_faces(itype) + &
                  &  ijk_to_spatial_numbering(input_data%number_of_dimensions, num_total_nfaces_per_dir(:,itype), nface_ijk)
          end if
       end do
    end do

    if(present(num_ghost_cells)) then

       ! Cells global numbering
       call memalloc(num_local_cells+num_ghost_cells_,cells_gids,__FILE__,__LINE__)
       do icell = 1, num_local_cells+num_ghost_cells_
          call spatial_to_ijk_numbering(input_data%number_of_dimensions, num_total_cells_per_dir, icell, cell_ijk)
          cell_ijk = first_cell_ijk + cell_ijk
          cells_gids(icell) = 1 + ijk_to_spatial_numbering(input_data%number_of_dimensions, input_data%number_of_cells_per_dir, cell_ijk)
       end do

       ! Number of interface cells (=local-interior)
       num_itfc_cells = 1
       do idime = 0, input_data%number_of_dimensions - 1
          if(num_local_cells_per_dir(idime) > &
               & num_left_parts_per_dir(idime) + num_right_parts_per_dir(idime)) then
             num_itfc_cells = num_itfc_cells * &
                  & (num_local_cells_per_dir(idime) - num_left_parts_per_dir(idime) - num_right_parts_per_dir(idime))
          !if(num_local_cells_per_dir(idime)>1) then
          !   num_itfc_cells = num_itfc_cells * (num_local_cells_per_dir(idime)-2)
          else
             num_itfc_cells = 0
             exit
          end if
       end do
       num_itfc_cells  = num_local_cells - num_itfc_cells

       ! List ghost cells (a permutation array to reorder cells with ghost at the end)
       call memalloc(num_local_cells+num_ghost_cells_, cell_is_ghost, __FILE__,__LINE__)
       cell_is_ghost=0
       ighost_cell = 0
       ilocal_cell = 0
       itfc_cells = 0
       do icell = 1, num_local_cells+num_ghost_cells_
          call spatial_to_ijk_numbering(input_data%number_of_dimensions, num_total_cells_per_dir, icell, cell_ijk)
          do idime = 0, input_data%number_of_dimensions - 1
             if(input_data%is_dir_periodic(idime)==0) then
                if(    (num_left_parts_per_dir(idime)==1 .and.cell_ijk(idime)==0).or. &
                     & (num_right_parts_per_dir(idime)==1.and.cell_ijk(idime)==num_total_cells_per_dir(idime)-1) ) then ! cell is ghost
                   cell_is_ghost(icell) = num_local_cells + num_ghost_cells_- ighost_cell
                   ighost_cell = ighost_cell + 1
                   exit
                end if
             end if
          end do
          if(cell_is_ghost(icell)==0) then
             do idime = 0, input_data%number_of_dimensions - 1
                if(input_data%is_dir_periodic(idime)==0) then
                   if((num_left_parts_per_dir(idime)==1 .and.cell_ijk(idime)==1).or. &
                        &  (num_right_parts_per_dir(idime)==1.and.cell_ijk(idime)==num_total_cells_per_dir(idime)-2) ) then ! cell is interface
                      cell_is_ghost(icell) = num_local_cells - itfc_cells
                      itfc_cells = itfc_cells + 1
                      exit
                   end if
                end if
             end do
          end if
          if(cell_is_ghost(icell)==0) then
             ilocal_cell = ilocal_cell + 1
             cell_is_ghost(icell) = ilocal_cell
          end if
       end do
       assert(ilocal_cell == num_local_cells-num_itfc_cells)
       assert(itfc_cells == num_itfc_cells)
       assert(ighost_cell == num_ghost_cells_)

       ! List ghost cells and compute interface cells pointers
       call memalloc(num_itfc_cells,lst_itfc_cells, __FILE__,__LINE__)
       call memalloc(num_itfc_cells+1,ptr_ext_neighs_per_itfc_cell,__FILE__,__LINE__)
       itfc_cells = 0
       do icell = 1, num_local_cells+num_ghost_cells_
          if(cell_is_ghost(icell)>num_local_cells-num_itfc_cells.and.cell_is_ghost(icell)<=num_local_cells) then ! cell is interface
             call spatial_to_ijk_numbering(input_data%number_of_dimensions, num_total_cells_per_dir, icell, cell_ijk)
             index = 0
             do iface=1,polytope_tree%get_number_n_faces()
                if(polytope_tree%get_n_face_dimension(iface)<input_data%number_of_dimensions) then
                   count_it = .false.
                   do idime = 0, input_data%number_of_dimensions - 1
                   !   if( input_data%is_dir_periodic(idime)==0) then
                   !      neighbor_ijk(idime) = cell_ijk(idime) - 1 + &
                   !           & 2 * polytope_tree%n_face_dir_coordinate(iface,idime) + &
                   !           & polytope_tree%n_face_dir_is_fixed(iface,idime)
                   !   else if(input_data%is_dir_periodic(idime)==1) then 
                   !      neighbor_ijk(idime) = mod(cell_ijk(idime) - 1 + &
                   !           & 2 * polytope_tree%n_face_dir_coordinate(iface,idime) + &
                   !           & polytope_tree%n_face_dir_is_fixed(iface,idime), &
                   !           & num_total_nfaces_per_dir(idime,itype))
                   !   end if
                      neighbor_ijk(idime) = cell_ijk(idime) - 1 + &
                           & 2 * polytope_tree%n_face_dir_coordinate(iface,idime) + &
                           & polytope_tree%n_face_dir_is_fixed(iface,idime)
                      if(neighbor_ijk(idime)<0.or.neighbor_ijk(idime)>num_total_cells_per_dir(idime)-1) then
                         count_it = .false.  ! the neighbor is out of the domain.
                         exit
                      else if( (num_left_parts_per_dir(idime)==1.and.neighbor_ijk(idime)==0) .or. &
                           &   (num_right_parts_per_dir(idime)==1.and.neighbor_ijk(idime)==num_total_cells_per_dir(idime)-1)) then
                         count_it = .true.
                      end if
                   end do
                   if(count_it) index = index + 1
                end if
             end do
             if(index>0) then
                itfc_cells = itfc_cells + 1
                lst_itfc_cells(itfc_cells) = icell
                ptr_ext_neighs_per_itfc_cell(itfc_cells+1)=index
             end if
          end if
       end do
       assert(itfc_cells==num_itfc_cells)

       ! Point to head
       ptr_ext_neighs_per_itfc_cell(1)=1
       do itfc_cells=1,num_itfc_cells
          ptr_ext_neighs_per_itfc_cell(itfc_cells+1) = ptr_ext_neighs_per_itfc_cell(itfc_cells+1) + ptr_ext_neighs_per_itfc_cell(itfc_cells)
       end do

       ! List interface cells neighbors and neighbor parts
       call memalloc(ptr_ext_neighs_per_itfc_cell(num_itfc_cells+1)-1, lst_ext_neighs_gids ,__FILE__,__LINE__)
       call memalloc(ptr_ext_neighs_per_itfc_cell(num_itfc_cells+1)-1, lst_ext_neighs_part_ids ,__FILE__,__LINE__)
       itfc_cells = 1
       do icell = 1, num_local_cells+num_ghost_cells_
          if(cell_is_ghost(icell)>num_local_cells-num_itfc_cells.and.cell_is_ghost(icell)<=num_local_cells) then ! cell is interface
             call spatial_to_ijk_numbering(input_data%number_of_dimensions, num_total_cells_per_dir, icell, cell_ijk)
             index = 0
             do iface=1,polytope_tree%get_number_n_faces()
                if(polytope_tree%get_n_face_dimension(iface)<input_data%number_of_dimensions) then
                   count_it = .false.
                   do idime = 0, input_data%number_of_dimensions - 1
                      neighbor_ijk(idime) = cell_ijk(idime) - 1 + &
                           & 2 * polytope_tree%n_face_dir_coordinate(iface,idime) + &
                           & polytope_tree%n_face_dir_is_fixed(iface,idime)
                      if(neighbor_ijk(idime)<0.or.neighbor_ijk(idime)>num_total_cells_per_dir(idime)-1) then
                         count_it = .false.  ! the neighbor is out of the domain.
                         exit
                      else if( (num_left_parts_per_dir(idime)==1.and.neighbor_ijk(idime)==0)) then
                         neighbor_part_ijk(idime)=part_ijk(idime)-1
                         count_it = .true.
                      else if( (num_right_parts_per_dir(idime)==1.and.neighbor_ijk(idime)==num_total_cells_per_dir(idime)-1)) then
                         neighbor_part_ijk(idime)=part_ijk(idime)+1
                         count_it = .true.
                      else
                         neighbor_part_ijk(idime)=part_ijk(idime)
                      end if
                   end do
                   if(count_it) then
                      lst_ext_neighs_gids(ptr_ext_neighs_per_itfc_cell(itfc_cells)+index)= cells_gids(1 + &
                           &   ijk_to_spatial_numbering( input_data%number_of_dimensions, &
                           &                             num_total_cells_per_dir, neighbor_ijk))

                      ! This should work too (test it!):
                      ! neighbor_ijk = first_cell_ijk + neighbor_ijk
                      ! lst_ext_neighs_gids(ptr_ext_neighs_per_itfc_cell(itfc_cells)+index)= 1 + &
                      !      &   ijk_to_spatial_numbering( input_data%number_of_dimensions, &
                      !      &                             input_data%number_of_cells_per_dir, neighbor_ijk)

                      lst_ext_neighs_part_ids(ptr_ext_neighs_per_itfc_cell(itfc_cells)+index)= 1 + &
                           &   ijk_to_spatial_numbering( input_data%number_of_dimensions, &
                           &                             input_data%number_of_parts_per_dir, neighbor_part_ijk)
                      index = index + 1
                   end if
                end if
             end do
             if(index>0) itfc_cells = itfc_cells + 1
          end if
       end do
       assert(itfc_cells==num_itfc_cells+1)

    end if

    ! vef global numbering (if needed), coordinates and boundary ids
    if(present(num_ghost_cells)) call memalloc(num_local_vefs,vefs_gids,__FILE__,__LINE__)
    call memalloc(SPACE_DIM,num_vertices,coordinates,__FILE__,__LINE__)
    call memalloc(num_local_vefs,boundary_id,__FILE__,__LINE__)
    boundary_id=-1
    do iface=1,polytope_tree%get_number_n_faces()
       if(polytope_tree%get_n_face_dimension(iface)<input_data%number_of_dimensions.and. &
            & polytope_tree%n_face_coordinate(iface)==0) then
          itype = polytope_tree%n_face_type(iface)
          do iface_of_itype = num_total_n_faces(itype), num_total_n_faces(itype+1) - 1
             call spatial_to_ijk_numbering(input_data%number_of_dimensions, num_total_nfaces_per_dir(:,itype), &
                  &                        iface_of_itype + 1 - num_total_n_faces(itype), nface_ijk)
             nface_ijk = first_cell_ijk + nface_ijk
             if(present(num_ghost_cells)) &
             vefs_gids(iface_of_itype) = num_global_n_faces(itype) + &
                  &                      ijk_to_spatial_numbering( input_data%number_of_dimensions, &
                  &                                                num_global_nfaces_per_dir(:,itype), &
                  &                                                nface_ijk )
             if(itype==0) then
                do idime = 0, input_data%number_of_dimensions - 1 
                   coordinates(idime+1,iface_of_itype) = real(nface_ijk(idime),rp) / real(input_data%number_of_cells_per_dir(idime),rp)
                end do
             end if
             index = 0
             do idime = 0, input_data%number_of_dimensions - 1 
                if(input_data%is_dir_periodic(idime)==0) then ! Not periodic
                   if(  polytope_tree%n_face_dir_is_fixed(iface,idime)==0.and.nface_ijk(idime)==0) then 
                      ! idime bit is already 0
                   else if(polytope_tree%n_face_dir_is_fixed(iface,idime)==0.and.nface_ijk(idime)==num_global_nfaces_per_dir(idime,itype)-1) then 
                      index = ibset( index, idime )
                   else
                      index = ibset( index, input_data%number_of_dimensions + idime ) ! Fix this coordinate
                   end if
                end if
             end do
             boundary_id(iface_of_itype) = index
          end do
       end if
    end do

    call memfree( num_global_n_faces, __FILE__,__LINE__)
    call memfree( num_total_n_faces, __FILE__,__LINE__)
    call memfree( num_global_nfaces_per_dir, __FILE__,__LINE__)
    call memfree( num_total_nfaces_per_dir, __FILE__,__LINE__)

    call node_array%free()
    call polytope_tree%free()

  end subroutine uniform_hex_mesh_generator_generate_connectivities

  pure function ijk_to_spatial_numbering_ip(num_dimensions, num_per_dim, ijk)
    implicit none
    integer(ip)           , intent(in) :: num_dimensions
    integer(ip)           , intent(in) :: num_per_dim(0:SPACE_DIM-1) 
    integer(ip)           , intent(in) :: ijk(0:SPACE_DIM-1) 
    integer(ip) :: ijk_to_spatial_numbering_ip
    integer(ip) :: idime, jdime
    integer(ip) :: previous
    ijk_to_spatial_numbering_ip = 0
    do idime = 0, num_dimensions - 1
       previous = 1
       do jdime = 0, idime - 1 
          previous = previous * num_per_dim(jdime)
       end do
       ijk_to_spatial_numbering_ip = ijk_to_spatial_numbering_ip + previous*ijk(idime)
    end do
  end function ijk_to_spatial_numbering_ip

  pure function ijk_to_spatial_numbering_igp(num_dimensions, num_per_dim, ijk)
    implicit none
    integer(ip)           , intent(in) :: num_dimensions
    integer(igp)          , intent(in) :: num_per_dim(0:SPACE_DIM-1) 
    integer(ip)           , intent(in) :: ijk(0:SPACE_DIM-1) 
    integer(igp) :: ijk_to_spatial_numbering_igp
    integer(ip)  :: idime, jdime
    integer(igp) :: previous
    ijk_to_spatial_numbering_igp = 0
    do idime = 0, num_dimensions - 1
       previous = 1
       do jdime = 0, idime - 1 
          previous = previous * num_per_dim(jdime)
       end do
       ijk_to_spatial_numbering_igp = ijk_to_spatial_numbering_igp + previous*ijk(idime)
    end do
  end function ijk_to_spatial_numbering_igp

  pure subroutine spatial_to_ijk_numbering(num_dimensions, num_per_dim, spatial_numbering, ijk)
    implicit none
    integer(ip)           , intent(in)  :: num_dimensions
    integer(ip)           , intent(in)  :: num_per_dim(0:SPACE_DIM-1) 
    integer(ip)           , intent(in)  :: spatial_numbering
    integer(ip)           , intent(out) :: ijk(0:SPACE_DIM-1) 
    integer(ip) :: idime,j

    j = spatial_numbering - 1          ! To make it 0-based (assuming spatial_numbering is 1-based)
    do idime = 0, num_dimensions - 1
       ijk(idime) = mod(j,num_per_dim(idime))
       j = j / num_per_dim(idime)
    end do

  end subroutine spatial_to_ijk_numbering
  
  
  ! subroutine uniform_hex_mesh_generator_generate_fempar_triangulation_arrays(part_id,                         &
  !                                                                            uniform_hex_mesh_descriptor,     &
  !                                                                            hex_boundary_set_ids_descriptor, &
  !                                                                            num_local_cells,                 &
  !                                                                            num_local_vefs,                  &
  !                                                                            ptr_vefs_per_cell,               &
  !                                                                            lst_vefs_lids,                   &
  !                                                                            vef_lid_to_vertex_lid,           &
  !                                                                            cell_gids,                       &
  !                                                                            vefs_gids,                       &
  !                                                                            vefs_boundary_set_ids,           &
  !                                                                            vefs_geometry_ids,               &
  !                                                                            vertex_coordinates,              &  
  !                                                                            num_itfc_cells,                  &
  !                                                                            lst_itfc_cells,                  &
  !                                                                            ptr_ext_neighs_per_itfc_cell,    &
  !                                                                            lst_ext_neighs_gids,             &
  !                                                                            lst_ext_neighs_part_ids)
                                                                             
  !   implicit none
  !   integer(ip)                            , intent(in)    :: part_id
  !   type(uniform_hex_mesh_descriptor_t)    , intent(in)    :: uniform_hex_mesh_descriptor
  !   type(hex_boundary_set_ids_descriptor_t), intent(in)    :: hex_boundary_set_ids_descriptor
  !   integer(ip)                            , intent(out)   :: num_local_cells
  !   integer(ip)                            , intent(out)   :: num_local_vefs
  !   integer(ip)           , allocatable    , intent(inout) :: ptr_vefs_per_cell(:)            ! Size = num_local_cells + 1
  !   integer(ip)           , allocatable    , intent(inout) :: lst_vefs_lids(:)                ! Size = ptr_vefs_per_cell(num_local_cells+1)-1
  !   integer(ip)           , allocatable    , intent(inout) :: vef_lid_to_vertex_lid(:)        ! Size = num_local_vefs (-1 if vef_lid is not a vertex)
  !   integer(igp)          , allocatable    , intent(inout) :: cell_gids(:)                    ! Size = num_local_cells 
  !   integer(igp)          , allocatable    , intent(inout) :: vefs_gids(:)                    ! Size = num_local_vefs
  !   integer(ip)           , allocatable    , intent(inout) :: vefs_boundary_set_ids(:)        ! Size = num_local_vefs
  !   integer(ip)           , allocatable    , intent(inout) :: vefs_geometry_ids(:)            ! Size = num_local_vefs
  !   real(rp)              , allocatable    , intent(inout) :: vertex_coordinates(:,:)         ! Size = (number_dimensions, num_local_vertices)
  !   integer(ip) , optional                 , intent(out)   :: num_itfc_cells                  ! NONE or ALL OPTIONAL ARGUMENTS MUST BE PRESENT LOGIC
  !   integer(ip) , optional, allocatable    , intent(inout) :: lst_itfc_cells(:)               ! Size = num_itfc_cells 
  !   integer(ip) , optional, allocatable    , intent(inout) :: ptr_ext_neighs_per_itfc_cell(:) ! Size = num_itfc_cells + 1
  !   integer(ip) , optional, allocatable    , intent(inout) :: lst_ext_neighs_gids(:)          ! Size = ptr_ext_neighs_per_itfc_cell(num_itfc_cells + 1)-1
  !   integer(ip) , optional, allocatable    , intent(inout) :: lst_ext_neighs_part_ids(:)      ! Size = ptr_ext_neighs_per_itfc_cell(num_itfc_cells + 1)-1
  ! end subroutine uniform_hex_mesh_generator_generate_fempar_triangulation_arrays

  ! ********* class(uniform_hex_mesh_descriptor_t) TBPS ********
  ! ************************************************************
  ! subroutine uniform_hex_mesh_descriptor_describe_mesh(this,parameter_list)
  !   !-----------------------------------------------------------------------------------------------!
  !   !   This subroutine generates geometry data to construct a structured mesh                      !
  !   !-----------------------------------------------------------------------------------------------!
  !   implicit none
  !   class(uniform_hex_mesh_descriptor_t), intent(inout) :: this
  !   type(ParameterList_t)               , intent(in)    :: parameter_list
  !   ! Locals
  !   integer(ip)          :: istat,mc
  !   integer(ip), allocatable :: ne_size(:),np_size(:),ns_size(:),disc_size(:),peri_size(:),nb_size(:)
  !   integer(ip), allocatable :: dl_size(:),o_size(:),st_size(:),sb_size(:)
  !   integer(ip), allocatable :: ne(:),np(:),ns(:),disc(:),peri(:),nb(:)
  !   real(rp)   , allocatable :: dl(:),o(:),st(:),sb(:)

  !   ! Fill uniform_mesh_descriptor


  !   istat = 0
  !   istat = istat + parameter_list%getshape(key = number_elements_name, shape = ne_size)
  !   istat = istat + parameter_list%getshape(key = number_parts_name, shape = np_size)
  !   istat = istat + parameter_list%getshape(key = number_sockets_name, shape = ns_size)
  !   istat = istat + parameter_list%getshape(key = discretization_type_name, shape = disc_size)
  !   istat = istat + parameter_list%getshape(key = periodic_boundaries_name, shape = peri_size)
  !   istat = istat + parameter_list%getshape(key = number_elements_boundary_name, shape = nb_size)
  !   istat = istat + parameter_list%get(key = material_case_name, value = mc)
  !   istat = istat + parameter_list%getshape(key = domain_length_name, shape = dl_size)
  !   istat = istat + parameter_list%getshape(key = origin_name, shape = o_size)
  !   istat = istat + parameter_list%getshape(key = stretching_parameter_name, shape = st_size)
  !   istat = istat + parameter_list%getshape(key = size_boundary_name, shape = sb_size)
  !   check(istat==0)      
  !   call memalloc(ne_size(1),ne,__FILE__,__LINE__)
  !   call memalloc(np_size(1),np,__FILE__,__LINE__)
  !   call memalloc(ns_size(1),ns,__FILE__,__LINE__)
  !   call memalloc(disc_size(1),disc,__FILE__,__LINE__)
  !   call memalloc(peri_size(1),peri,__FILE__,__LINE__)
  !   call memalloc(nb_size(1),nb,__FILE__,__LINE__)
  !   call memalloc(dl_size(1),dl,__FILE__,__LINE__)
  !   call memalloc(o_size(1),o,__FILE__,__LINE__)
  !   call memalloc(st_size(1),st,__FILE__,__LINE__)
  !   call memalloc(sb_size(1),sb,__FILE__,__LINE__)
  !   istat = 0
  !   istat = istat + parameter_list%get(key = 'number_elements', value = ne)
  !   istat = istat + parameter_list%get(key = 'number_parts', value = np)
  !   istat = istat + parameter_list%get(key = 'number_sockets', value = ns)
  !   istat = istat + parameter_list%get(key = 'discretization_type', value = disc)
  !   istat = istat + parameter_list%get(key = 'periodic_boundaries', value = peri)
  !   istat = istat + parameter_list%get(key = 'number_elements_boundary', value = nb)
  !   istat = istat + parameter_list%get(key = 'domain_length', value = dl)
  !   istat = istat + parameter_list%get(key = 'origin', value = o)
  !   istat = istat + parameter_list%get(key = 'stretching_parameter', value = st)
  !   istat = istat + parameter_list%get(key = 'size_boundary', value = sb)
  !   check(istat==0)

  !   this%ntdix   = disc(1) ! Type of discretization in x (0=uniform, 1=cubic, 2=tanh, 3=imh+unif, 4:imh+tanh)
  !   this%ntdiy   = disc(2) ! Type of discretization in y (0=uniform, 1=cubic, 2=tanh, 3=imh+unif, 4:imh+tanh)
  !   this%mater   = mc      ! Material case
  !   this%neblx   = nb(1)   ! Number of elements in the x boundary layer
  !   this%nebly   = nb(2)   ! Number of elements in the y boundary layer
  !   this%xleng   = dl(1)   ! Size of the domain in x
  !   this%yleng   = dl(2)   ! Size of the domain in y
  !   this%x0      = o(1)    ! Origin x-coordinate
  !   this%y0      = o(2)    ! Origin y-coordinate
  !   this%xstret  = st(1)   ! Stretching parameter
  !   this%ystret  = st(2)   ! Stretching parameter
  !   this%xlengbl = sb(1)   ! Size of the boundary layer in x
  !   this%ylengbl = sb(2)   ! Size of the boundary layer in y
  !   if(size(disc,1)==3) this%ntdiz   = disc(3) ! Type of discretization in z (0=uniform, 1=cubic, 2=tanh, 3=imh+unif, 4:imh+tanh)
  !   if(size(nb,1)==3)   this%neblz   = nb(3)   ! Number of elements in the z boundary layer
  !   if(size(o,1)==3)    this%z0      = o(3)    ! Origin z-coordinate
  !   if(size(dl,1)==3)   this%zleng   = dl(3)   ! Size of the domain in z
  !   if(size(st,1)==3)   this%zstret  = st(3)   ! Stretching parameter
  !   if(size(sb,1)==3)   this%zlengbl = sb(3)   ! Size of the boundary layer in z
  
  !   ! Dimension
  !   if(size(ne,1)==2) then
  !      this%ndime=2     
  !   else
  !      if(ne(3)==0) then
  !         this%ndime=2
  !      else
  !         this%ndime=3
  !      end if
  !   end if 

  !   ! Geometry order of interpolations
  !   this%pdegr = 1

  !   ! Partition info
  !   this%nedir(1) = ne(1)
  !   this%nedir(2) = ne(2)
  !   if(size(ne,1)==3) this%nedir(3) = ne(3)
  !   this%npdir(1) = np(1)
  !   this%npdir(2) = np(2)
  !   if(size(np,1)==3) this%npdir(3) = np(3)
  !   check(ns(1)>0)
  !   this%nsckt(1) = this%npdir(1)/ns(1)
  !   check(ns(2)>0)
  !   this%nsckt(2) = this%npdir(2)/ns(2)
  !   this%nparts = this%npdir(1)*this%npdir(2)
  !   if(this%ndime==3) then
  !      check(size(ns,1)==3)
  !      check(ns(3)>0)
  !      this%nsckt(3) = this%npdir(3)/ns(3)
  !      this%nparts = this%nparts*this%npdir(3)
  !   else
  !      this%nsckt(3) = 1
  !      this%npdir(3) = 1
  !   end if
  !   this%isper(1) = peri(1)
  !   this%isper(2) = peri(2)
  !   if(size(peri,1)==3) this%isper(3) = peri(3)   

  !   ! Deallocate
  !   deallocate(ne_size)
  !   deallocate(np_size)
  !   deallocate(ns_size)
  !   deallocate(disc_size)
  !   deallocate(peri_size)
  !   deallocate(nb_size)
  !   deallocate(dl_size)
  !   deallocate(o_size)
  !   deallocate(st_size)
  !   deallocate(sb_size)
  !   call memfree(ne,__FILE__,__LINE__)
  !   call memfree(np,__FILE__,__LINE__)
  !   call memfree(ns,__FILE__,__LINE__)
  !   call memfree(disc,__FILE__,__LINE__)
  !   call memfree(peri,__FILE__,__LINE__)
  !   call memfree(nb,__FILE__,__LINE__)
  !   call memfree(dl,__FILE__,__LINE__)
  !   call memfree(o,__FILE__,__LINE__)
  !   call memfree(st,__FILE__,__LINE__)
  !   call memfree(sb,__FILE__,__LINE__)
  ! end subroutine uniform_hex_mesh_descriptor_describe_mesh

  
  ! ! ********* Helper stand-alone subroutines (legacy code) *****
  ! ! ************************************************************
  !  !==================================================================================================
  ! subroutine uniform_hex_mesh_descriptor_generate_geom_and_topo_size(this,geom_size,topo_size)
  !   !-----------------------------------------------------------------------------------------------!
  !   !   This subroutine generates geom_size_t type from uniform_mesh_descriptor_t and reference_element_t types            !
  !   !-----------------------------------------------------------------------------------------------!
  !   implicit none
  !   class(uniform_hex_mesh_descriptor_t), intent(in)  :: this
  !   type(geom_size_t)                   , intent(out) :: geom_size
  !   type(topo_size_t)                   , intent(out) :: topo_size


  !   ! Local variables
  !   integer(ip) :: pdegr,idime,jdime
    
  !   ! Local variables
  !   integer(ip) :: pdime,i,nfaux,ndaux
  !   integer(ip), allocatable :: auxv(:,:)

  !   pdegr = this%pdegr

  !   ! Directional sizes
  !   geom_size%npdir = this%npdir
  !   geom_size%nedir = this%nedir
  !   geom_size%nsckt = this%nsckt
  !   geom_size%nnode = (pdegr+1)**this%ndime  
  !   geom_size%npsoc = this%npdir/this%nsckt   
  !   geom_size%nedom = this%nedir/this%npdir   
  !   geom_size%npdom = pdegr*geom_size%nedom + 1  
  !   geom_size%ncorn = this%npdir + (1-this%isper) 
  !   do idime=1,this%ndime
  !      geom_size%nedge(idime,:) = this%npdir + 1 - this%isper
  !      geom_size%nface(idime,:) = this%npdir
  !      if(this%isper(idime)==0) then
  !         geom_size%nface(idime,idime) = geom_size%nface(idime,idime) + 1 
  !      end if
  !   end do
  !   if(this%ndime==2) geom_size%nface = 0

  !   ! Total sizes
  !   geom_size%npdomt  = 1
  !   geom_size%nedomt  = 1
  !   geom_size%ncornt  = 1
  !   geom_size%nedget  = 1
  !   geom_size%nfacet  = 1
  !   geom_size%nedgett = 0
  !   geom_size%nfacett = 0
  !   geom_size%neghost = 1
  !   do idime=1,this%ndime
  !      geom_size%npdomt = geom_size%npdomt*geom_size%npdom(idime)
  !      geom_size%nedomt = geom_size%nedomt*geom_size%nedom(idime)
  !      geom_size%ncornt = geom_size%ncornt*geom_size%ncorn(idime)
  !      geom_size%neghost = geom_size%neghost*(geom_size%nedom(idime)+2)
  !      do jdime=1,this%ndime
  !         geom_size%nedget(idime) = geom_size%nedget(idime)*geom_size%nedge(idime,jdime)
  !         geom_size%nfacet(idime) = geom_size%nfacet(idime)*geom_size%nface(idime,jdime) 
  !      end do
  !      geom_size%nedgett = geom_size%nedgett + geom_size%nedget(idime)
  !      geom_size%nfacett = geom_size%nfacett + geom_size%nfacet(idime)
  !   end do
    
  !   ! Auxiliar vector of components
  !   if(this%ndime==2) then
  !      call memalloc(2,1,auxv,__FILE__,__LINE__)
  !      auxv(1,1) = 2
  !      auxv(2,1) = 1
  !   else if(this%ndime==3) then
  !      call memalloc(3,2,auxv,__FILE__,__LINE__)
  !      auxv(1,:) = (/2,3/)
  !      auxv(2,:) = (/1,3/)
  !      auxv(3,:) = (/1,2/)
  !   end if

  !   topo_size%notot = 0
  !   topo_size%nctot = 1
  !   topo_size%ncglb = 1
  !   topo_size%neglb = 1
  !   topo_size%ndtot = 0
  !   topo_size%nftot = 0
  !   do pdime=this%ndime,1,-1
  !      nfaux = 1
  !      ndaux = 1
  !      do i=1,this%ndime-1
  !         nfaux = nfaux*geom_size%nedom(auxv(pdime,i))
  !         ndaux = ndaux*(geom_size%nedom(auxv(pdime,i))+1)
  !         topo_size%nfdom(auxv(pdime,i),pdime) = geom_size%nedom(auxv(pdime,i))
  !         topo_size%nfglb(auxv(pdime,i),pdime) = geom_size%nedir(auxv(pdime,i))          
  !         topo_size%nddom(auxv(pdime,i),pdime) = (geom_size%nedom(auxv(pdime,i))+1)
  !         topo_size%ndglb(auxv(pdime,i),pdime) = (geom_size%nedir(auxv(pdime,i))+1)        
  !      end do
  !      topo_size%nfdir(pdime) = (geom_size%nedom(pdime)+1)*nfaux
  !      topo_size%nddir(pdime) = geom_size%nedom(pdime)*ndaux
  !      topo_size%nfdom(pdime,pdime) = (geom_size%nedom(pdime)+1)
  !      topo_size%nfglb(pdime,pdime) = (geom_size%nedir(pdime)+1)
  !      topo_size%nddom(pdime,pdime) = geom_size%nedom(pdime)
  !      topo_size%ndglb(pdime,pdime) = geom_size%nedir(pdime)
  !      topo_size%ncglb = topo_size%ncglb*(geom_size%nedir(pdime)+1)
  !      topo_size%neglb = topo_size%neglb*geom_size%nedir(pdime)
  !      topo_size%nctot = topo_size%nctot*geom_size%npdom(pdime)
  !      topo_size%ndtot = topo_size%ndtot + topo_size%nddir(pdime)
  !      topo_size%nftot = topo_size%nftot + topo_size%nfdir(pdime)
  !   end do 
  !   do pdime=this%ndime,1,-1
  !      nfaux = 1
  !      ndaux = 1
  !      do i=1,this%ndime-1
  !         nfaux = nfaux*geom_size%nedir(auxv(pdime,i))
  !         ndaux = ndaux*(geom_size%nedir(auxv(pdime,i))+1)
  !      end do
  !      topo_size%ndsum(pdime) = geom_size%nedir(pdime)*ndaux
  !      topo_size%nfsum(pdime) = (geom_size%nedir(pdime)+1)*nfaux
  !   end do
  !   if(this%ndime==2) then
  !      topo_size%nftot = 0
  !      topo_size%nfdir = 0
  !      topo_size%nfdom = 0
  !      topo_size%nfglb = 0
  !   end if
  !   topo_size%notot = topo_size%nftot + topo_size%ndtot + topo_size%nctot
  !   topo_size%noglb = topo_size%notot*this%nparts

  !   ! Deallocate auxiliar vector of components
  !   call memfree(auxv,__FILE__,__LINE__)
  ! end subroutine uniform_hex_mesh_descriptor_generate_geom_and_topo_size


!  !==================================================================================================
!  subroutine generate_uniform_triangulation(lpart,gdata,bdata,trian,bcond,mdist)
!    !-----------------------------------------------------------------------------------------------!
!    !   This subroutine generates a triangulation, boundary conditions and mesh_distribution for a  !
!    !   structured mesh                                                                             !
!    !-----------------------------------------------------------------------------------------------!
!    implicit none
!    integer(ip)                          , intent(in)  :: lpart
!    type(uniform_mesh_descriptor_t)      , intent(in)  :: gdata
!    type(uniform_conditions_descriptor_t), intent(in)  :: bdata
!    type(triangulation_t)                , intent(out) :: trian
!    type(conditions_t)                   , intent(out) :: bcond
!    type(mesh_distribution_t), optional  , intent(out) :: mdist
!    
!    ! Locals
!    type(geom_size_t) :: gsize
!    type(topo_size_t) :: tsize
!    integer(ip)     :: ijkpart(3),ielem,inode

!    ! Allocatables
!    integer(ip) , allocatable   :: pextn(:),lextp(:),lexte(:),lextm(:)
!    integer(igp), allocatable   :: lextn(:)

!    ! Geometrical sizes
!    call structured_geom_size_create(gdata,gsize)

!    ! Topological sizes
!    call structured_topo_size_create(gdata,gsize,tsize)

!    ! Transform lpart to ijk numeration
!    call global_to_ijk(lpart,gdata%nsckt(1:gdata%ndime),gsize%npsoc(1:gdata%ndime),gdata%ndime,ijkpart)

!    ! Create triangulation
!    if(present(mdist)) then
!       call triangulation_create(gsize%neghost,trian)
!    else
!       call triangulation_create(gsize%nedomt,trian)
!    end if
!    trian%num_elems = gsize%nedomt
!    trian%num_dims  = gdata%ndime

!    ! Create geometric reference element
!    call uniform_mesh_create_reference_fe(gdata,trian%reference_fe_geo_list(1))

!    ! Triangulation, bcond and distribution generation
!    if(present(mdist)) then
!       call structured_mesh_gen(ijkpart,gdata,gsize,tsize,bdata%poin,bdata%line,bdata%surf,    &
!            &                   trian,bcond,nmap=mdist%nmap,emap=mdist%emap,pextn=mdist%pextn, &
!            &                   lextn=mdist%lextn,lextp=mdist%lextp)
!    else
!       call structured_mesh_gen(ijkpart,gdata,gsize,tsize,bdata%poin,bdata%line,bdata%surf, &
!            &                   trian,bcond)
!    end if

!    ! Dual triangulation
!    if(.not.present(mdist)) call triangulation_to_dual(trian)

!    ! Create mesh_distribution
!    if(present(mdist)) then
!       mdist%ipart  = lpart
!       mdist%nparts = gdata%nparts
!       mdist%nebou  = mdist%emap%nb
!       mdist%nnbou  = mdist%nmap%nb
!       call memalloc(mdist%nebou,mdist%lebou,__FILE__,__LINE__)
!       do ielem = mdist%emap%ni+1,mdist%emap%ni+mdist%emap%nb
!          mdist%lebou(ielem-mdist%emap%ni) = ielem
!       end do
!       call memalloc(mdist%nnbou,mdist%lnbou,__FILE__,__LINE__)
!       do inode = mdist%nmap%ni+1,mdist%nmap%ni+mdist%nmap%nb
!          mdist%lnbou(inode-mdist%nmap%ni) = inode
!       end do
!    end if

!  end subroutine generate_uniform_triangulation

!  !================================================================================================
!  subroutine structured_mesh_gen(ijkpart,gdata,gsize,tsize,poin,line,surf,trian,nodes, &
!       &                         nmap,emap,pextn,lextn,lextp)
!    !-----------------------------------------------------------------------
!    ! 
!    !-----------------------------------------------------------------------
!    implicit none
!    integer(ip)                        , intent(in)    :: ijkpart(3)
!    type(uniform_mesh_descriptor_t)    , intent(in)    :: gdata
!    type(geom_size_t)                  , intent(in)    :: gsize
!    type(topo_size_t)                  , intent(in)    :: tsize
!    type(triangulation_t)              , intent(inout) :: trian
!    type(conditions_t)                 , intent(in)    :: poin,line,surf
!    type(conditions_t)                 , intent(out)   :: nodes
!    type(map_igp_t)          , optional, intent(inout) :: nmap,emap
!    integer(ip) , allocatable, optional, intent(out)   :: pextn(:),lextp(:)
!    integer(igp), allocatable, optional, intent(out)   :: lextn(:)

!    ! Local variables
!    integer(ip)              :: nparts,ndime,isper(3)
!    integer(ip)              :: cnt(3),subgl(gdata%ndime),idime,lpart,pni,eni,pnb,enb
!    integer(ip)              :: pextn_cnt(2),nelbl(3)

!    ! Local allocatables 
!    integer(ip) , allocatable :: npnumg(:),npnumt(:),nenum(:)    
!    integer(igp), allocatable :: l2ge(:),l2gp(:)   
!    real(rp)    , allocatable :: coord(:,:)

!    ! Unpack variables
!    nparts = gdata%nparts
!    ndime  = gdata%ndime
!    isper  = gdata%isper

!    ! Allocate 
!    !call memalloc(gsize%npdomt,l2gp,__FILE__,__LINE__)
!    call memalloc(tsize%notot,l2gp,__FILE__,__LINE__)
!    call memalloc(gsize%nedomt,l2ge,__FILE__,__LINE__)
!    call memalloc(gsize%nedomt,nenum,__FILE__,__LINE__)
!    call memalloc(gsize%npdomt,npnumg,__FILE__,__LINE__)
!    call memalloc(tsize%notot,npnumt,__FILE__,__LINE__)
!    call memalloc(ndime,gsize%npdomt,coord,__FILE__,__LINE__)

!    ! Initialize counter
!    cnt(1) = 1; cnt(2) = 1; cnt(3) = 1

!    ! Part id
!    call globalid_2l(ijkpart,gsize%nsckt,gsize%npsoc,ndime,lpart)

!    ! Create bounary conditions
!    call conditions_create(poin%ncode,poin%nvalu,tsize%notot,nodes)
!   
!    ! Interior nodes/elements
!    call volu_loop(ijkpart,ndime,gsize,tsize,gdata,npnumg,npnumt,nenum,coord,cnt)
!    call face_loop(ijkpart,ndime,gsize,tsize,gdata,gdata%isper,surf,npnumg,npnumt,nenum,coord, &
!         &         cnt,inter,nodes)
!    call edge_loop(ijkpart,ndime,gsize,tsize,gdata,gdata%isper,line,surf,npnumg,npnumt,nenum,  &
!         &         coord,cnt,inter,nodes)
!    call corn_loop(ijkpart,ndime,gsize,tsize,gdata,gdata%isper,poin,line,surf,npnumg,npnumt,   &
!         &         nenum,coord,cnt,inter,nodes)
!    !pni = cnt(3)-1
!    pni = cnt(1)-1
!    eni = cnt(2)-1

!    ! Boundary nodes/elements
!    call face_loop(ijkpart,ndime,gsize,tsize,gdata,gdata%isper,surf,npnumg,npnumt,nenum,coord, &
!         &         cnt,bound)
!    call edge_loop(ijkpart,ndime,gsize,tsize,gdata,gdata%isper,line,surf,npnumg,npnumt,nenum,  &
!         &         coord,cnt,bound,nodes)
!    call corn_loop(ijkpart,ndime,gsize,tsize,gdata,gdata%isper,poin,line,surf,npnumg,npnumt,   &
!         &         nenum,coord,cnt,bound,nodes)
!    !pnb = cnt(3) - pni - 1
!    pnb = cnt(1) - pni - 1
!    enb = cnt(2) - eni - 1

!    ! Global part numeration
!    do idime=1,ndime
!       subgl(idime) = (ijkpart(idime)-1)*gsize%nedom(idime)
!    end do

!    ! Calculate lnods, coord and l2g vector for emap and nmap
!    call generic_l2g(subgl,npnumg,npnumt,nenum,ndime,gdata%pdegr,gsize,tsize,gdata,l2ge,l2gp,trian, &
!         &           coord)

!    ! Fill nmap
!    if(present(nmap)) then
!       !call map_alloc(gsize%npdomt,int(tsize%ncglb,igp),nmap)
!       call map_alloc(tsize%notot,int(tsize%noglb,igp),nmap)
!       nmap%ni  = pni
!       nmap%nb  = pnb
!       nmap%ne  = 0
!       nmap%l2g = l2gp
!    end if

!    ! Fill emap
!    if(present(emap)) then
!       call map_alloc(gsize%nedomt,int(tsize%neglb,igp),emap)
!       emap%ni  = eni
!       emap%nb  = enb
!       emap%ne  = 0
!       emap%l2g = l2ge
!    end if
!    
!    ! Compute adjacencies
!    if(present(pextn)) then
!       check(present(lextp))
!       check(present(lextn))
!       call memalloc(enb+1,pextn,__FILE__,__LINE__)
!       nelbl(1) = gdata%neblx
!       nelbl(2) = gdata%nebly
!       nelbl(3) = gdata%neblz
!       pextn_cnt = 0
!       pextn(1) = 1
!       call face_loop_adj(ijkpart,ndime,gsize,tsize,gdata%isper,emap%nb,pextn_cnt,pextn,0)
!       call edge_loop_adj(ijkpart,ndime,gsize,tsize,gdata%isper,emap%nb,pextn_cnt,pextn,0)
!       call corn_loop_adj(ijkpart,ndime,gsize,tsize,gdata%isper,emap%nb,pextn_cnt,pextn,0)
!       pextn_cnt = 0
!       call memalloc(pextn(enb+1)-1,lextn,__FILE__,__LINE__)
!       call memalloc(pextn(enb+1)-1,lextp,__FILE__,__LINE__)
!       lextn = 0; lextp = 0
!       call face_loop_adj(ijkpart,ndime,gsize,tsize,gdata%isper,emap%nb,pextn_cnt,pextn,1, &
!            &             lextn,lextp,gdata%mater,nelbl)
!       call edge_loop_adj(ijkpart,ndime,gsize,tsize,gdata%isper,emap%nb,pextn_cnt,pextn,1, &
!            &             lextn,lextp,gdata%mater,nelbl)
!       call corn_loop_adj(ijkpart,ndime,gsize,tsize,gdata%isper,emap%nb,pextn_cnt,pextn,1, &
!            &             lextn,lextp,gdata%mater,nelbl)
!    end if

!    ! Deallocate auxiliar vectors
!    call memfree( l2gp,__FILE__,__LINE__)
!    call memfree( l2ge,__FILE__,__LINE__)
!    call memfree(npnumt,__FILE__,__LINE__)
!    call memfree(npnumg,__FILE__,__LINE__)
!    call memfree(nenum,__FILE__,__LINE__)
!    call memfree(coord,__FILE__,__LINE__)

!  end subroutine structured_mesh_gen

!  !==================================================================================================
!  subroutine volu_loop(ijkpart,ndime,gsize,tsize,gdata,npnumg,npnumt,nenum,coord,cnt)
!    !-----------------------------------------------------------------------
!    ! 
!    !-----------------------------------------------------------------------
!    implicit none
!    integer(ip)                    , intent(in)    :: ijkpart(3),ndime
!    type(geom_size_t)              , intent(in)    :: gsize
!    type(topo_size_t)              , intent(in)    :: tsize
!    type(uniform_mesh_descriptor_t), intent(in)    :: gdata
!    integer(ip)                    , intent(inout) :: npnumg(:),npnumt(:),nenum(:),cnt(3)
!    real(rp)                       , intent(inout) :: coord(:,:)

!    integer(ip) :: i,j,k,ijkpoin(3),ijkelem(3),glnum,npdomk,nedomk,nddomk,aux_cnt
!    integer(ip) :: pdime,ijkface(3),ijkedge(3),jdime
!    integer(ip), allocatable :: auxv(:,:)

!    ! Auxiliar
!    if(ndime==2) then
!       npdomk=2
!       nedomk=2
!       call memalloc(2,1,auxv,__FILE__,__LINE__)
!       auxv(1,1) = 2
!       auxv(2,1) = 1
!    else if(ndime==3) then
!       npdomk = gsize%npdom(3) - 1
!       nedomk = gsize%nedom(3) - 1
!       call memalloc(3,2,auxv,__FILE__,__LINE__)
!       auxv(1,:) = (/2,3/)
!       auxv(2,:) = (/1,3/)
!       auxv(3,:) = (/1,2/)
!    end if
!    
!    ! Interior nodes loop
!    do i=2,gsize%npdom(1)-1
!       do j=2,gsize%npdom(2)-1
!          do k=2,npdomk

!             ! Point identifier (inside the part)
!             ijkpoin = (/i,j,k/)
!             call globalid(ijkpoin,gsize%npdom,ndime,glnum)

!             ! Global to local numeration (inside the part)
!             npnumg(glnum) = cnt(3)
!             npnumt(glnum) = cnt(1)

!             ! Coordinates
!             call coord_ijk(ijkpoin,ijkpart,gsize%npdir,gsize%nedom,gsize%nedir, &
!                  &         ndime,gdata,coord(:,cnt(3)))

!             ! Update counter
!             cnt(1) = cnt(1) + 1
!             cnt(3) = cnt(3) + 1

!          end do
!       end do
!    end do

!    ! Interior edges loop
!    aux_cnt = tsize%nctot
!    do pdime=ndime,1,-1
!       if(ndime==2) then
!          nddomk=2
!       elseif(ndime==3) then
!          nddomk=tsize%nddom(auxv(pdime,2),pdime)-1
!       end if
!       do i=1,tsize%nddom(pdime,pdime)
!          do j=2,tsize%nddom(auxv(pdime,1),pdime)-1
!             do k=2,nddomk

!                ! Face identifier (inside the part)
!                ijkedge(pdime) = i
!                if(ndime==2) then
!                   ijkedge(auxv(pdime,:)) =  j
!                else
!                   ijkedge(auxv(pdime,:)) =  (/j,k/)
!                end if
!                call globalid(ijkedge,tsize%nddom(:,pdime),ndime,glnum)

!                ! Global to local numeration (inside the part)
!                npnumt(aux_cnt+glnum) = cnt(1)

!                ! Update counter
!                cnt(1) = cnt(1) + 1

!             end do
!          end do
!       end do
!       aux_cnt = aux_cnt + tsize%nddir(pdime)
!    end do

!    ! Interior faces loop
!    aux_cnt = tsize%nctot + tsize%ndtot
!    if(ndime==3) then
!       do pdime=ndime,1,-1
!          do i=2,tsize%nfdom(pdime,pdime)-1
!             do j=1,tsize%nfdom(auxv(pdime,1),pdime)
!                do k=1,tsize%nfdom(auxv(pdime,2),pdime)

!                   ! Face identifier (inside the part)
!                   ijkface(pdime) = i
!                   ijkface(auxv(pdime,:)) = (/j,k/)
!                   call globalid(ijkface,tsize%nfdom(:,pdime),ndime,glnum)

!                   ! Global to local numeration (inside the part)
!                   npnumt(aux_cnt+glnum) = cnt(1)

!                   ! Update counter
!                   cnt(1) = cnt(1) + 1

!                end do
!             end do
!          end do
!          aux_cnt = aux_cnt + tsize%nfdir(pdime)
!       end do
!    end if
!    
!    ! Interior elements loop
!    do i=2,gsize%nedom(1)-1
!       do j=2,gsize%nedom(2)-1
!          do k=2,nedomk

!             ! Element identifier (inside the part)
!             ijkelem = (/i,j,k/)
!             call globalid(ijkelem,gsize%nedom,ndime,glnum)

!             ! Global to local numeration (inside the part)
!             nenum(glnum) = cnt(2)

!             ! Update counter
!             cnt(2) = cnt(2) + 1

!          end do
!       end do
!    end do

!    ! Deallocate auxiliar vector
!    call memfree(auxv,__FILE__,__LINE__)

!  end subroutine volu_loop

!  !==================================================================================================
!  subroutine face_loop(ijkpart,ndime,gsize,tsize,gdata,isper,surf,npnumg,npnumt,nenum,coord, &
!       &               cnt,case,nodes)
!    !-----------------------------------------------------------------------
!    ! 
!    !-----------------------------------------------------------------------
!    implicit none
!    integer(ip)                    , intent(in)    :: ijkpart(3),isper(3),ndime,case
!    type(geom_size_t)              , intent(in)    :: gsize
!    type(topo_size_t)              , intent(in)    :: tsize
!    type(uniform_mesh_descriptor_t), intent(in)    :: gdata
!    type(conditions_t)             , intent(in)    :: surf
!    integer(ip)                    , intent(inout) :: npnumg(:),npnumt(:),nenum(:),cnt(3)
!    real(rp)                       , intent(inout) :: coord(:,:)
!    type(conditions_t), optional   , intent(inout) :: nodes
!    
!    integer(ip) :: auxv(3,2),i,j,k,pdime,iface,ijkpoin(3),ijkelem(3),ijkface(3),glnum,flag,neigh(2)
!    integer(ip) :: lface(2),surf_cnt
!    integer(ip) :: idir,jdir,idime,jdime,aux_cnt_f,aux_cnt_d,ijkedge(3)

!    ! Auxiliar vector of components
!    auxv(1,:) = (/2,3/)
!    auxv(2,:) = (/1,3/)
!    auxv(3,:) = (/1,2/)

!    ! Face nodes loop
!    surf_cnt = 1
!    aux_cnt_f = tsize%nctot + tsize%ndtot
!    if(ndime==3) then
!       do pdime = ndime,1,-1
!          do iface = 0,1

!             ! Check if it is a boundary
!             call check_face_boundary(ijkpart,pdime,iface,gsize%npdir,gsize%nsckt,gsize%npsoc,ndime, &
!                  &                   isper,auxv,flag,neigh,lface)

!             ! Check case
!             if(flag==case) then

!                ! Point identifier (inside the part)
!                ijkpoin(pdime) = (gsize%npdom(pdime)-1)*iface+1

!                do i=2,gsize%npdom(auxv(pdime,1))-1
!                   do j=2,gsize%npdom(auxv(pdime,2))-1

!                      ! Point identifier (inside the part)
!                      ijkpoin(auxv(pdime,:)) = (/i,j/)
!                      call globalid(ijkpoin,gsize%npdom,ndime,glnum)

!                      ! Global to local numeration (inside the part)
!                      npnumg(glnum) = cnt(3)
!                      npnumt(glnum) = cnt(1)

!                      ! Coordinates
!                      call coord_ijk(ijkpoin,ijkpart,gsize%npdir,gsize%nedom,gsize%nedir, &
!                           &         ndime,gdata,coord(:,cnt(3)))

!                      ! Boundary conditions
!                      if(case==inter) then
!                         nodes%code(:,cnt(1)) = surf%code(:,surf_cnt)
!                         nodes%valu(:,cnt(1)) = surf%valu(:,surf_cnt)
!                      end if

!                      ! Update point counter
!                      cnt(1) = cnt(1) + 1
!                      cnt(3) = cnt(3) + 1

!                   end do
!                end do

!                ! Elemental edges loop
!                do idir = 1,2
!                   jdir = 3-idir
!                   idime = auxv(pdime,idir)
!                   i = (tsize%nddom(pdime,idime)-1)*iface+1
!                   aux_cnt_d = tsize%nctot
!                   do jdime=1,3-idime
!                      aux_cnt_d = aux_cnt_d + tsize%nddir(4-jdime)
!                   end do
!                   do j=1,tsize%nddom(idime,idime)
!                      do k=2,tsize%nddom(auxv(pdime,jdir),idime)-1

!                         ! Face identifier (inside the part)
!                         ijkedge(pdime) = i
!                         ijkedge(auxv(pdime,idir)) = j
!                         ijkedge(auxv(pdime,jdir)) = k
!                         call globalid(ijkedge,tsize%nddom(:,idime),ndime,glnum)

!                         ! Global to local numeration (inside the part)
!                         npnumt(aux_cnt_d+glnum) = cnt(1)

!                         ! Boundary conditions
!                         if(case==inter) then
!                            nodes%code(:,cnt(1)) = surf%code(:,surf_cnt)
!                            nodes%valu(:,cnt(1)) = surf%valu(:,surf_cnt)
!                         end if

!                         ! Update counter
!                         cnt(1) = cnt(1) + 1

!                      end do
!                   end do
!                end do

!                ! Elemental faces loop
!                i = (tsize%nfdom(pdime,pdime)-1)*iface+1
!                do j=1,tsize%nfdom(auxv(pdime,1),pdime)
!                   do k=1,tsize%nfdom(auxv(pdime,2),pdime)

!                      ! Face identifier (inside the part)
!                      ijkface(pdime) = i
!                      ijkface(auxv(pdime,:)) = (/j,k/)
!                      call globalid(ijkface,tsize%nfdom(:,pdime),ndime,glnum)

!                      ! Global to local numeration (inside the part)
!                      npnumt(aux_cnt_f+glnum) = cnt(1)

!                      ! Boundary conditions
!                      if(case==inter) then
!                         nodes%code(:,cnt(1)) = surf%code(:,surf_cnt)
!                         nodes%valu(:,cnt(1)) = surf%valu(:,surf_cnt)
!                      end if

!                      ! Update counter
!                      cnt(1) = cnt(1) + 1

!                   end do
!                end do

!                if(iface.gt.min(1,abs(gsize%nedom(pdime)-1))) cycle

!                ! Element identifier (inside the part)
!                ijkelem(pdime) = (gsize%nedom(pdime)-1)*iface+1

!                do i=2,gsize%nedom(auxv(pdime,1))-1
!                   do j=2,gsize%nedom(auxv(pdime,2))-1

!                      ! Point identifier (inside the part)
!                      ijkelem(auxv(pdime,:)) = (/i,j/)
!                      call globalid(ijkelem,gsize%nedom,ndime,glnum)

!                      ! Global to local numeration (inside the part)
!                      nenum(glnum) = cnt(2)

!                      ! Update element counter
!                      cnt(2) = cnt(2) + 1

!                   end do
!                end do
!             end if

!             ! Update surface counter
!             surf_cnt = surf_cnt + 1

!          end do
!          
!          ! Update elemental faces counter
!          aux_cnt_f = aux_cnt_f + tsize%nfdir(pdime)
!       end do
!    end if

!  end subroutine face_loop

! !===================================================================================================
!  subroutine face_loop_adj(ijkpart,ndime,gsize,tsize,isper,nb,pextn_cnt,pextn,case,lextn,lextp, &
!       &                   mcase,nelbl)
!    !-----------------------------------------------------------------------
!    ! 
!    !-----------------------------------------------------------------------
!    implicit none
!    integer(ip)           , intent(in)    :: ijkpart(3),isper(3),ndime,case,nb
!    type(geom_size_t)       , intent(in)    :: gsize
!    type(topo_size_t)       , intent(in)    :: tsize
!    integer(ip)           , intent(inout) :: pextn_cnt(2),pextn(nb+1)
!    integer(ip) , optional, intent(inout) :: lextp(:)
!    integer(igp), optional, intent(inout) :: lextn(:)
!    integer(ip) , optional, intent(in)    :: mcase,nelbl(3)
!    
!    integer(ip)  :: auxv(3,2),ijkelem(3),ijkneigh(3),subgl(3),subgl_ijk(3),ijkpart_neigh(3)
!    integer(ip)  :: i,j,pdime,iface,ielem,jelem,glnum,flag,gpart,lface(2),k,kedge,idime,neigh(2)
!    integer(igp) :: gneigh 

!    if(ndime==3) then

!       ! Auxiliar vector of components
!       auxv(1,:) = (/2,3/)
!       auxv(2,:) = (/1,3/)
!       auxv(3,:) = (/1,2/)

!       ! Face nodes loop
!       do pdime = ndime,1,-1
!          do iface = 0,1

!             ! Check if it is a boundary
!             call check_face_boundary(ijkpart,pdime,iface,gsize%npdir,gsize%nsckt,gsize%npsoc,ndime, &
!                  &                   isper,auxv,flag,neigh,lface)

!             ! Only boundary elements
!             if(flag==bound) then

!                if(case==do_count) then

!                   do i=2,gsize%nedom(auxv(pdime,1))-1
!                      do j=2,gsize%nedom(auxv(pdime,2))-1

!                         ! Update boundary element counter
!                         pextn_cnt(1) = pextn_cnt(1) + 1

!                         do ielem = -1,1
!                            do jelem = -1,1

!                               ! Update boundary neighbour counter
!                               pextn_cnt(2) = pextn_cnt(2) + 1

!                            end do
!                         end do

!                         ! Fill pextn
!                         pextn(pextn_cnt(1)+1) = pextn_cnt(2) + 1

!                      end do
!                   end do

!                elseif(case==do_list) then

!                   ! Element identifier (inside the part)
!                   ijkelem(pdime) = (gsize%nedom(pdime)-1)*iface+1

!                   ! Neighbour element identifier
!                   ijkneigh(pdime) = (gsize%nedom(pdime)-1)*(1-iface)+1

!                   ! Neighbour part
!                   ijkpart_neigh = ijkpart
!                   ijkpart_neigh(pdime) = ijkpart(pdime) + (2*iface-1)
!                   call check_part_boundary(ijkpart_neigh,gsize%npdir,ndime,isper)
!                   call globalid_2l(ijkpart_neigh,gsize%nsckt,gsize%npsoc,ndime,gpart)

!                   ! Global part numeration
!                   do i=1,ndime
!                      subgl(i) = (ijkpart_neigh(i)-1)*gsize%nedom(i)
!                   end do

!                   do i=2,gsize%nedom(auxv(pdime,1))-1
!                      do j=2,gsize%nedom(auxv(pdime,2))-1

!                         ! Update boundary element counter
!                         pextn_cnt(1) = pextn_cnt(1) + 1

!                         ! Point identifier (inside the part)
!                         ijkelem(auxv(pdime,:)) = (/i,j/)
!                         call globalid(ijkelem,gsize%nedom,ndime,glnum)

!                         do ielem = -1,1
!                            do jelem = -1,1

!                               ! Update boundary neighbour counter
!                               pextn_cnt(2) = pextn_cnt(2) + 1

!                               ! Neighbour element identifier (local)
!                               ijkneigh(auxv(pdime,:)) = (/i+ielem,j+jelem/)
!                               subgl_ijk = subgl + ijkneigh
!                               call globalid(subgl_ijk,gsize%nedir,ndime,gneigh)

!                               ! Fill adjacencies info
!                               lextn(pextn_cnt(2)) = gneigh
!                               lextp(pextn_cnt(2)) = gpart

!                            end do
!                         end do

!                      end do
!                   end do
!                end if
!             end if
!          end do
!       end do

!    end if

!  end subroutine face_loop_adj

!  !==================================================================================================
!  subroutine edge_loop(ijkpart,ndime,gsize,tsize,gdata,isper,line,surf,npnumg,npnumt,nenum, &
!       &               coord,cnt,case,nodes)
!    !-----------------------------------------------------------------------
!    ! 
!    !-----------------------------------------------------------------------
!    implicit none
!    integer(ip)                    , intent(in)    :: ijkpart(3),isper(3),ndime,case
!    type(geom_size_t)              , intent(in)    :: gsize
!    type(topo_size_t)              , intent(in)    :: tsize
!    type(uniform_mesh_descriptor_t), intent(in)    :: gdata
!    type(conditions_t)             , intent(in)    :: line,surf
!    integer(ip)                    , intent(inout) :: npnumg(:),npnumt(:),nenum(:),cnt(3)
!    real(rp)                       , intent(inout) :: coord(:,:)
!    type(conditions_t)             , intent(inout) :: nodes
! 
!    integer(ip) :: i,j,k,pdime,iedge,jedge,ijedge(2),ijkpoin(3),ijkelem(3),ijkvert(3),ijkedge(3)
!    integer(ip) :: inipo,ledge(2,2),line_cnt,neigh(2*(ndime-1))
!    integer(ip) :: glnum,flag,kron(ndime),gsurf,inode,aux_surf(4,2),gsurf_aux(4),aux_cnt
!    integer(ip), allocatable :: auxv(:,:)

!    ! Auxiliar vector of components 
!    if(ndime==2) then
!       call memalloc(2,1,auxv,__FILE__,__LINE__)
!       auxv(1,1) = 2
!       auxv(2,1) = 1
!    else if(ndime==3) then
!       call memalloc(3,2,auxv,__FILE__,__LINE__)
!       auxv(1,:) = (/2,3/)
!       auxv(2,:) = (/1,3/)
!       auxv(3,:) = (/1,2/)
!    end if

!    ! Edge nodes loop
!    line_cnt = 0
!    aux_cnt = tsize%nctot
!    do pdime = ndime,1,-1
!       do iedge = 0,1!min(1,nedom(auxv(pdime,1))-1)
!          do jedge =0,1!min(1,nedom(auxv(pdime,2))-1)

!             ! Check if it is a boundary
!             ijedge = (/iedge,jedge/)
!             call check_edge_boundary(ijkpart,pdime,ijedge,gsize%npdir,gsize%nsckt,gsize%npsoc,ndime, &
!                  &                   isper,auxv,flag,neigh,ledge)

!             ! Update edge counter
!             line_cnt = line_cnt + 1

!             ! Check case
!             if(flag==case) then

!                ! Initial object point
!                inipo = cnt(1)

!                ! Point identifier (inside the part)
!                ijkvert(auxv(pdime,:)) = ijedge(1:ndime-1)
!                ijkvert(pdime) = 0
!                do i=1,ndime
!                   ijkpoin(i) = ijkvert(i)*(gsize%npdom(i)-1) + 1
!                end do

!                do i=1,gsize%npdom(pdime)-2

!                   ! Point identifier (inside the part)
!                   call kronec(pdime,ndime,kron)
!                   ijkpoin(1:ndime) = ijkpoin(1:ndime) + kron
!                   call globalid(ijkpoin,gsize%npdom,ndime,glnum)

!                   ! Global to local numeration (inside the part)
!                   npnumg(glnum) = cnt(3)
!                   npnumt(glnum) = cnt(1)

!                   ! Coordinates
!                   call coord_ijk(ijkpoin,ijkpart,gsize%npdir,gsize%nedom,gsize%nedir, &
!                        &         ndime,gdata,coord(:,cnt(3)))

!                   ! Boundary conditions
!                   if(case==inter) then
!                      nodes%code(:,cnt(1)) = line%code(:,line_cnt)
!                      nodes%valu(:,cnt(1)) = line%valu(:,line_cnt)
!                   end if

!                   ! Update point counter
!                   cnt(1) = cnt(1) + 1
!                   cnt(3) = cnt(3) + 1
!                end do

!                ! Elemental edges loop
!                j = iedge*(tsize%nddom(auxv(pdime,1),pdime)-1) + 1
!                if(ndime==3) k = jedge*(tsize%nddom(auxv(pdime,2),pdime)-1) + 1
!                do i = 1,tsize%nddom(pdime,pdime)

!                   ! Face identifier (inside the part)
!                   ijkedge(pdime) = i
!                   ijkedge(auxv(pdime,1)) = j
!                   if(ndime==3) ijkedge(auxv(pdime,2)) = k
!                   call globalid(ijkedge,tsize%nddom(:,pdime),ndime,glnum)

!                   ! Global to local numeration (inside the part)
!                   npnumt(aux_cnt+glnum) = cnt(1)

!                   ! Boundary conditions
!                   if(case==inter) then
!                      nodes%code(:,cnt(1)) = line%code(:,line_cnt)
!                      nodes%valu(:,cnt(1)) = line%valu(:,line_cnt)
!                   end if

!                   ! Update counter
!                   cnt(1) = cnt(1) + 1

!                end do

!                if(case==bound) then

!                   ! Boundary conditions for non-global edges
!                   k = 0
!                   do i=1,2*(ndime-1)
!                      if(neigh(i)>0) then
!                         k=k+1
!                      end if
!                   end do
!                   if(ndime==3.and.k==2) then
!                      call aux_surf_num(pdime,aux_surf,gsurf_aux)
!                      do i=1,4
!                         if(neigh(aux_surf(i,1))>0.and.neigh(aux_surf(i,2))>0) then
!                            gsurf = gsurf_aux(i)
!                         end if
!                      end do
!                      do inode=inipo,cnt(1)-1
!                         nodes%code(:,inode) = surf%code(:,gsurf)
!                         nodes%valu(:,inode) = surf%valu(:,gsurf)
!                      end do
!                   end if

!                end if

!                if(iedge.gt.min(1,abs(gsize%nedom(auxv(pdime,1))-1))) cycle
!                if(ndime==3) then
!                   if(jedge.gt.min(1,abs(gsize%nedom(auxv(pdime,2))-1))) cycle
!                end if

!                ! Element identifier (inside the part)
!                do i=1,ndime
!                   ijkelem(i) = ijkvert(i)*(gsize%nedom(i)-1) + 1
!                end do

!                do i=2,gsize%nedom(pdime)-1

!                   ! Element identifier (inside the part)
!                   call kronec(pdime,ndime,kron)
!                   ijkelem(1:ndime) = ijkelem(1:ndime) + kron
!                   call globalid(ijkelem,gsize%nedom,ndime,glnum)

!                   ! Global to local numeration (inside the part)
!                   nenum(glnum) = cnt(2)

!                   ! Update element counter
!                   cnt(2) = cnt(2) + 1
!                end do

!             end if

!             if(ndime==2) exit
!          end do
!       end do

!       ! Update elemental edge counter
!       aux_cnt = aux_cnt + tsize%nddir(pdime)
!    end do

!    ! Deallocate auxiliar vector
!    call memfree(auxv,__FILE__,__LINE__)

!  end subroutine edge_loop

!  !===================================================================================================
!  subroutine edge_loop_adj(ijkpart,ndime,gsize,tsize,isper,nb,pextn_cnt,pextn,case,lextn,lextp, &
!       &                   mcase,nelbl)
!    !-----------------------------------------------------------------------
!    ! 
!    !-----------------------------------------------------------------------
!    implicit none
!    integer(ip)           , intent(in)    :: ijkpart(3),isper(3),ndime,case,nb
!    type(geom_size_t)       , intent(in)    :: gsize
!    type(topo_size_t)       , intent(in)    :: tsize
!    integer(ip)           , intent(inout) :: pextn_cnt(2),pextn(nb+1)
!    integer(ip) , optional, intent(inout) :: lextp(:)
!    integer(igp), optional, intent(inout) :: lextn(:)
!    integer(ip) , optional, intent(in)    :: mcase,nelbl(3)
!    
!    integer(ip) :: ijkelem(3),neigh(2*(ndime-1)),ijedge(2),ijkvert(3),ijkneigh(3)
!    integer(ip) :: subgl(3),subgl_ijk(3),ijkpart_neigh(3),ledge(2,2),kron(ndime),lpart
!    integer(ip) :: i,j,pdime,iedge,jedge,kedge,medge,ielem,glnum,flag,gpart,idime
!    integer(ip) :: iaux,jaux,iaux2,jaux2
!    integer(igp) :: gneigh
!    integer(ip), allocatable :: auxv(:,:)

!    ! Auxiliar vector of components
!    if(ndime==2) then
!       call memalloc(2,1,auxv,__FILE__,__LINE__)
!       auxv(1,1) = 2
!       auxv(2,1) = 1
!    else if(ndime==3) then
!       call memalloc(3,2,auxv,__FILE__,__LINE__)
!       auxv(1,:) = (/2,3/)
!       auxv(2,:) = (/1,3/)
!       auxv(3,:) = (/1,2/)
!    end if

!    ! Global numbering
!    call globalid_2l(ijkpart,gsize%nsckt,gsize%npsoc,ndime,lpart)

!    ! Edge nodes loop
!    do pdime = ndime,1,-1
!       do iedge = 0,1
!          iaux = 2*iedge-1
!          do jedge =0,1
!             jaux = 2*jedge-1

!             ! Check if it is a boundary
!             ijedge = (/iedge,jedge/)
!             call check_edge_boundary(ijkpart,pdime,ijedge,gsize%npdir,gsize%nsckt,gsize%npsoc, &
!                  &                   ndime,isper,auxv,flag,neigh,ledge)

!             ! Only boundary elements
!             if(flag==bound) then

!                if(case==do_count) then

!                   do i=2,gsize%nedom(pdime)-1

!                      ! Update boundary element counter
!                      pextn_cnt(1) = pextn_cnt(1) + 1

!                      do ielem = -1,1
!                         do kedge = -1,1
!                            do medge = -1,1
!                               if(kedge==iaux.or.medge==jaux) then

!                                  iaux2 = kedge
!                                  jaux2 = medge
!                                  if(kedge.ne.iaux) iaux2 = 0
!                                  if(medge.ne.jaux) jaux2 = 0

!                                  ! Neighbour part
!                                  ijkpart_neigh(pdime) = ijkpart(pdime)
!                                  if(ndime==2) then
!                                     ijkpart_neigh(auxv(pdime,:)) = ijkpart(auxv(pdime,:)) + iaux2
!                                  else
!                                     ijkpart_neigh(auxv(pdime,:)) = ijkpart(auxv(pdime,:)) + (/iaux2,jaux2/)
!                                  end if
!                                  call check_part_boundary(ijkpart_neigh,gsize%npdir,ndime,isper)
!                                  call globalid_2l(ijkpart_neigh,gsize%nsckt,gsize%npsoc,ndime,gpart)

!                                  if(gpart>0.and.gpart.ne.lpart) then

!                                     ! Update boundary neighbour counter
!                                     pextn_cnt(2) = pextn_cnt(2) + 1

!                                  end if
!                               end if
!                               if(ndime==2) exit
!                            end do
!                         end do
!                      end do

!                      ! Fill pextn
!                      pextn(pextn_cnt(1)+1) = pextn_cnt(2) + 1             

!                   end do

!                elseif(case==do_list) then

!                   ! Element identifier (inside the part)
!                   ijkvert(auxv(pdime,:)) = ijedge(1:ndime-1)
!                   ijkvert(pdime) = 0
!                   do i=1,ndime
!                      ijkelem(i) = ijkvert(i)*(gsize%nedom(i)-1) + 1
!                   end do

!                   do i=2,gsize%nedom(pdime)-1

!                      ! Element identifier (inside the part)
!                      call kronec(pdime,ndime,kron)
!                      ijkelem(1:ndime) = ijkelem(1:ndime) + kron
!                      call globalid(ijkelem,gsize%nedom,ndime,glnum)

!                      do kedge = -1,1
!                         do medge = -1,1
!                            if(kedge==iaux.or.medge==jaux) then
!                               
!                               iaux2 = kedge
!                               jaux2 = medge
!                               if(kedge.ne.iaux) iaux2 = 0
!                               if(medge.ne.jaux) jaux2 = 0

!                               ! Neighbour part
!                               ijkpart_neigh(pdime) = ijkpart(pdime)
!                               if(ndime==2) then
!                                  ijkpart_neigh(auxv(pdime,:)) = ijkpart(auxv(pdime,:)) + iaux2
!                               else
!                                  ijkpart_neigh(auxv(pdime,:)) = ijkpart(auxv(pdime,:)) + (/iaux2,jaux2/)
!                               end if

!                               ! Global part numeration
!                               do idime=1,ndime
!                                  subgl(idime) = (ijkpart(idime)-1)*gsize%nedom(idime)
!                                  if(ijkpart_neigh(idime)<1) then
!                                     subgl(idime) = gsize%npdir(idime)*gsize%nedom(idime)
!                                  elseif(ijkpart_neigh(idime)>gsize%npdir(idime)) then
!                                     subgl(idime) = -gsize%nedom(idime)
!                                  end if
!                               end do

!                               ! Check part boundary
!                               call check_part_boundary(ijkpart_neigh,gsize%npdir,ndime,isper)
!                               call globalid_2l(ijkpart_neigh,gsize%nsckt,gsize%npsoc,ndime,gpart)

!                               if(gpart>0.and.gpart.ne.lpart) then

!                                  do ielem = -1,1

!                                     ! Update boundary neighbour counter
!                                     pextn_cnt(2) = pextn_cnt(2) + 1

!                                     ! Neighbour element identifier (local)
!                                     ijkneigh(pdime) = ijkelem(pdime) + ielem
!                                     if(ndime==2) then
!                                        ijkneigh(auxv(pdime,:)) = ijkelem(auxv(pdime,:)) + kedge
!                                     else
!                                        ijkneigh(auxv(pdime,:)) = ijkelem(auxv(pdime,:)) +(/kedge,medge/)
!                                     end if
!                                     subgl_ijk = subgl + ijkneigh
!                                     call globalid(subgl_ijk,gsize%nedir,ndime,gneigh)

!                                     ! Fill adjacencies info
!                                     lextn(pextn_cnt(2)) = gneigh
!                                     lextp(pextn_cnt(2)) = gpart

!                                  end do
!                               end if
!                            end if
!                            if(ndime==2) exit
!                         end do
!                      end do
!                   end do
!                end if
!             end if
!             if(ndime==2) exit
!          end do
!       end do
!    end do

!    ! Deallocate auxiliar vector
!    call memfree(auxv,__FILE__,__LINE__)

!  end subroutine edge_loop_adj

!  !================================================================================================
!  subroutine corn_loop(ijkpart,ndime,gsize,tsize,gdata,isper,poin,line,surf,npnumg,npnumt, & 
!       &               nenum,coord,cnt,case,nodes)
!    !-----------------------------------------------------------------------
!    ! 
!    !-----------------------------------------------------------------------
!    implicit none
!    integer(ip)                    , intent(in)    :: ijkpart(3),isper(3),ndime,case
!    type(geom_size_t)              , intent(in)    :: gsize
!    type(topo_size_t)              , intent(in)    :: tsize
!    type(uniform_mesh_descriptor_t), intent(in)    :: gdata
!    type(conditions_t)             , intent(in)    :: poin,line,surf
!    integer(ip)                    , intent(inout) :: npnumg(:),npnumt(:),nenum(:),cnt(3)
!    real(rp)                       , intent(inout) :: coord(:,:)
!    type(conditions_t)             , intent(inout) :: nodes
!    
!    integer(ip) :: ivert,jvert,kvert,ijkpoin(3),ijkvert(3),glnum,flag,neigh(2**ndime)
!    integer(ip) :: i,j,k,lcorn(2,2,2),ijkelem(3)
!    integer(ip) :: poin_cnt,gline,gsurf
!    integer(ip) :: aux_line(2*(ndime-1)*ndime,2),aux_surf(2*ndime,4)
!    integer(ip), allocatable :: auxv(:,:)

!    ! Auxiliar vector of components
!    if(ndime==2) then
!       call memalloc(2,1,auxv,__FILE__,__LINE__)
!       auxv(1,1) = 2
!       auxv(2,1) = 1
!       aux_line(:,1) = (/3,1,2,1/)
!       aux_line(:,2) = (/4,2,4,3/)
!    else if(ndime==3) then
!       call memalloc(3,2,auxv,__FILE__,__LINE__)
!       auxv(1,:) = (/2,3/)
!       auxv(2,:) = (/1,3/)
!       auxv(3,:) = (/1,2/)
!       aux_line(:,1) = (/7,5,3,1,6,5,2,1,4,3,2,1/)
!       aux_line(:,2) = (/8,6,4,2,8,7,4,3,8,7,6,5/)
!       aux_surf(:,1) = (/2,1,3,1,5,1/)
!       aux_surf(:,2) = (/4,3,4,2,6,2/)
!       aux_surf(:,3) = (/6,5,7,5,7,3/)
!       aux_surf(:,4) = (/8,7,8,6,8,4/)
!    end if

!    ! Corner nodes loop
!    poin_cnt = 1
!    do ivert = 0,1
!       do jvert = 0,1
!          do kvert =0,1

!             ! Check if it is a boundary
!             ijkvert = (/ivert,jvert,kvert/)
!             call check_corn_boundary(ijkpart,ijkvert,gsize%npdir,gsize%nsckt,gsize%npsoc, &
!                  &                   ndime,isper,auxv,flag,neigh,lcorn)

!             ! Check case
!             if(flag==case) then

!                ! Point identifier (inside the part)
!                do i=1,ndime
!                   ijkpoin(i) = ijkvert(i)*(gsize%npdom(i)-1) + 1
!                end do
!                call globalid(ijkpoin,gsize%npdom,ndime,glnum)

!                ! Global to local numeration (inside the part)
!                npnumg(glnum) = cnt(3)
!                npnumt(glnum) = cnt(1)
!              
!                ! Coordinates
!                call coord_ijk(ijkpoin,ijkpart,gsize%npdir,gsize%nedom,gsize%nedir, &
!                     &         ndime,gdata,coord(:,cnt(3)))

!                ! Boundary conditions
!                if(case==inter) then
!                   nodes%code(:,cnt(1)) = poin%code(:,poin_cnt)
!                   nodes%valu(:,cnt(1)) = poin%valu(:,poin_cnt)
!                end if

!                ! Update point counter
!                cnt(1) = cnt(1) + 1
!                cnt(3) = cnt(3) + 1

!                if(case==bound) then

!                   ! Boundary conditions for non-global corners
!                   k = 0
!                   do i=1,2**ndime 
!                      if(neigh(i)>0) then
!                         k=k+1
!                      end if
!                   end do
!                   if(k==2) then
!                      do i=1,2*(ndime-1)*ndime
!                         if(neigh(aux_line(i,1))>0.and.neigh(aux_line(i,2))>0) then
!                            gline = i
!                         end if
!                      end do
!                      nodes%code(:,cnt(1)-1) = line%code(:,gline)
!                      nodes%valu(:,cnt(1)-1) = line%valu(:,gline)
!                   else if(k>2.and.k<2**ndime.and.ndime==3) then
!                      do i=1,2*ndime
!                         if(neigh(aux_surf(i,1))>0.and.neigh(aux_surf(i,2))>0.and. &
!                              neigh(aux_surf(i,3))>0.and.neigh(aux_surf(i,4))>0) then
!                            gsurf = i
!                         end if
!                      end do
!                      nodes%code(:,cnt(1)-1) = surf%code(:,gsurf)
!                      nodes%valu(:,cnt(1)-1) = surf%valu(:,gsurf)
!                   end if

!                end if

!                if(ivert.gt.min(1,abs(gsize%nedom(1)-1))) cycle
!                if(jvert.gt.min(1,abs(gsize%nedom(2)-1))) cycle
!                if(kvert.gt.min(1,abs(gsize%nedom(3)-1))) cycle

!                ! Element identifier (inside the part)
!                do i=1,ndime
!                   ijkelem(i) = ijkvert(i)*(gsize%nedom(i)-1) + 1
!                end do
!                call globalid(ijkelem,gsize%nedom,ndime,glnum)

!                ! Global to local numeration (inside the part)
!                nenum(glnum) = cnt(2)

!                ! Update element counter
!                cnt(2) = cnt(2) + 1

!             end if

!             ! Update corner counter
!             poin_cnt = poin_cnt + 1

!             if(ndime==2) exit
!          end do
!       end do
!    end do

!    ! Deallocate auxiliar vector
!    call memfree(auxv,__FILE__,__LINE__)

!  end subroutine corn_loop

!  !==================================================================================================
!  subroutine corn_loop_adj(ijkpart,ndime,gsize,tsize,isper,nb,pextn_cnt,pextn,case,lextn,lextp, &
!       &                   mcase,nelbl)
!    !-----------------------------------------------------------------------
!    ! 
!    !-----------------------------------------------------------------------
!    implicit none
!    integer(ip)           , intent(in)    :: ijkpart(3),isper(3),ndime,case,nb
!    type(geom_size_t)       , intent(in)    :: gsize
!    type(topo_size_t)       , intent(in)    :: tsize
!    integer(ip)           , intent(inout) :: pextn_cnt(2),pextn(nb+1)
!    integer(ip) , optional, intent(inout) :: lextp(:)
!    integer(igp), optional, intent(inout) :: lextn(:)
!    integer(ip) , optional, intent(in)    :: mcase,nelbl(3)
!    
!    integer(ip) :: ijkelem(3),neigh(2**ndime),ijkvert(3),ijkneigh(3),subgl(3),subgl_ijk(3)
!    integer(ip) :: i,j,k,ivert,jvert,kvert,glnum,flag,gpart,lcorn(2,2,2),idime,ijkpart_neigh(3)
!    integer(ip) :: iaux,jaux,kaux,iaux2,jaux2,kaux2,lpart
!    integer(ip), allocatable :: auxv(:,:)
!    integer(igp) :: gneigh

!    ! Auxiliar vector of components
!    if(ndime==2) then
!       call memalloc(2,1,auxv,__FILE__,__LINE__)
!       auxv(1,1) = 2
!       auxv(2,1) = 1
!    else if(ndime==3) then
!       call memalloc(3,2,auxv,__FILE__,__LINE__)
!       auxv(1,:) = (/2,3/)
!       auxv(2,:) = (/1,3/)
!       auxv(3,:) = (/1,2/)
!    end if       

!    ! Global numbering
!    call globalid_2l(ijkpart,gsize%nsckt,gsize%npsoc,ndime,lpart)

!    !Corner nodes loop
!    do ivert = 0,1
!       iaux = 2*ivert-1
!       do jvert = 0,1
!          jaux = 2*jvert-1
!          do kvert =0,1
!             kaux = 2*kvert-1
!    
!             ! Check if it is a boundary
!             ijkvert = (/ivert,jvert,kvert/)
!             call check_corn_boundary(ijkpart,ijkvert,gsize%npdir,gsize%nsckt,gsize%npsoc, &
!                  &                   ndime,isper,auxv,flag,neigh,lcorn)

!             ! Only boundary elements
!             if(flag==bound) then

!                if(case==do_count) then
!                   
!                   ! Update boundary element counter
!                   pextn_cnt(1) = pextn_cnt(1) + 1

!                   do i = -1,1
!                      do j = -1,1
!                         do k = -1,1
!                            if(i==iaux.or.j==jaux.or.k==kaux) then

!                               iaux2 = i
!                               jaux2 = j
!                               kaux2 = k
!                               if(i.ne.iaux) iaux2 = 0
!                               if(j.ne.jaux) jaux2 = 0
!                               if(k.ne.kaux) kaux2 = 0
!                               
!                               ! Neighbour part
!                               ijkpart_neigh = ijkpart + (/iaux2,jaux2,kaux2/)
!                               call check_part_boundary(ijkpart_neigh,gsize%npdir,ndime,isper)
!                               call globalid_2l(ijkpart_neigh,gsize%nsckt,gsize%npsoc,ndime,gpart)
!                               
!                               if(gpart>0.and.gpart.ne.lpart) then
!                                  
!                                  ! Update boundary neighbour counter
!                                  pextn_cnt(2) = pextn_cnt(2) + 1

!                               end if

!                            end if
!                            if(ndime==2) exit
!                         end do
!                      end do
!                   end do

!                   ! Fill pextn
!                   pextn(pextn_cnt(1)+1) = pextn_cnt(2) + 1

!                elseif(case==do_list) then

!                   ! Element identifier (inside the part)
!                   do i=1,ndime
!                      ijkelem(i) = ijkvert(i)*(gsize%nedom(i)-1) + 1
!                   end do

!                   ! Element identifier (inside the part)
!                   call globalid(ijkelem,gsize%nedom,ndime,glnum)

!                   do i = -1,1
!                      do j = -1,1
!                         do k = -1,1
!                            if(i==iaux.or.j==jaux.or.k==kaux) then

!                               iaux2 = i
!                               jaux2 = j
!                               kaux2 = k
!                               if(i.ne.iaux) iaux2 = 0
!                               if(j.ne.jaux) jaux2 = 0
!                               if(k.ne.kaux) kaux2 = 0
!                               
!                               ! Neighbour part
!                               ijkpart_neigh = ijkpart + (/iaux2,jaux2,kaux2/)

!                               ! Global part numeration
!                               do idime=1,ndime
!                                  subgl(idime) = (ijkpart(idime)-1)*gsize%nedom(idime)
!                                  if(ijkpart_neigh(idime)<1) then
!                                     subgl(idime) = gsize%npdir(idime)*gsize%nedom(idime)
!                                  elseif(ijkpart_neigh(idime)>gsize%npdir(idime)) then
!                                     subgl(idime) = -gsize%nedom(idime)
!                                  end if
!                               end do

!                               ! Check boundary partition
!                               call check_part_boundary(ijkpart_neigh,gsize%npdir,ndime,isper)
!                               call globalid_2l(ijkpart_neigh,gsize%nsckt,gsize%npsoc,ndime,gpart)
!                               
!                               if(gpart>0.and.gpart.ne.lpart) then

!                                  ! Update boundary neighbour counter
!                                  pextn_cnt(2) = pextn_cnt(2) + 1

!                                  ! Neighbour element identifier (local)
!                                  ijkneigh = ijkelem + (/i,j,k/)
!                                  subgl_ijk = subgl + ijkneigh
!                                  call globalid(subgl_ijk,gsize%nedir,ndime,gneigh)

!                                  ! Fill adjacencies info
!                                  lextn(pextn_cnt(2)) = gneigh
!                                  lextp(pextn_cnt(2)) = gpart

!                               end if
!                            end if
!                            if(ndime==2) exit
!                         end do
!                      end do
!                   end do
!                end if
!             end if
!             if(ndime==2) exit
!          end do
!       end do
!    end do

!    ! Deallocate auxiliar vector
!    call memfree(auxv,__FILE__,__LINE__)

!  end subroutine corn_loop_adj

!  !================================================================================================
!  subroutine coord_ijk(ijkpoin,ijkpart,npdir,nedom,nedir,ndime,msize,coord)
!    !-----------------------------------------------------------------------
!    ! 
!    !-----------------------------------------------------------------------
!    implicit none
!    integer(ip)                    , intent(in)    :: ijkpoin(3),ijkpart(3),npdir(3),nedom(3),nedir(3)
!    integer(ip)                    , intent(in)    :: ndime
!    type(uniform_mesh_descriptor_t), intent(in)    :: msize
!    real(rp)                       , intent(inout) :: coord(:)

!    integer(ip) :: i,ntdis(3),pdegr,nebl(3)
!    real(rp)    :: leng(3),coor0(3),stret(3),istret,lengbl(3)
!    integer(ip) :: newijkpoin,newnedir,newnedom,info
!    real(rp)    :: newleng

!    ! Unpack uniform_mesh_descriptor
!    pdegr   = msize%pdegr
!    leng(1) = msize%xleng
!    leng(2) = msize%yleng
!    leng(3) = msize%zleng
!    ntdis(1) = msize%ntdix
!    ntdis(2) = msize%ntdiy
!    ntdis(3) = msize%ntdiz
!    stret(1) = msize%xstret
!    stret(2) = msize%ystret
!    stret(3) = msize%zstret
!    nebl(1) = msize%neblx
!    nebl(2) = msize%nebly
!    nebl(3) = msize%neblz
!    lengbl(1) = msize%xlengbl
!    lengbl(2) = msize%ylengbl
!    lengbl(3) = msize%zlengbl

!    ! Set origin of coordinates
!    coor0(1) = msize%x0; coor0(2) = msize%y0; coor0(3) = msize%z0

!    ! Set origin of coordinates
!    coor0(1) = msize%x0; coor0(2) = msize%y0; coor0(3) = msize%z0

!    do i=1,ndime 
!       assert ( ntdis(i).ne.1 ) ! This case is not implemented
!       if(ntdis(i)==0) then
!          coord(i) = coor0(i) + leng(i)*((ijkpart(i)-1)*nedom(i)+(ijkpoin(i)-1)/real(pdegr))/nedir(i)
!       elseif(ntdis(i)==2) then
!          istret=stret(i)
!          coord(i) = coor0(i) - leng(i)*tanh(istret*(1.0_rp-2.0_rp*((ijkpart(i)-1)*nedom(i)+(ijkpoin(i)-1)/real(pdegr))/nedir(i)))/tanh(istret)
!       elseif(ntdis(i)==3.or.ntdis(i)==4) then ! boundary layers (Hunt's case: solid uniform + fluid unif==3 or tanh==4)
!          istret=stret(i)
!          
!!!$          !TO DO
!!!$          if((nebl(i)>0).and.(nedom(i)<nebl(i))) then
!!!$             write(*,*) 'boundary layer cannot be splitted over several sbds, dir ',i,'nedom',nedom,'nebl',nebl
!!!$             call mpi_barrier (1, info)
!!!$             call mpi_finalize( info )
!!!$             stop
!!!$          end if

!          !first boundary layer
!          if((ijkpart(i)==1).and.(ijkpoin(i)<=nebl(i))) then 
!             if(ntdis(i)==3) then
!                coord(i) = coor0(i) - lengbl(i) + lengbl(i)*(ijkpoin(i)-1)/nebl(i)
!             elseif(ntdis(i)==4) then
!                coord(i) = coor0(i) - leng(i) - lengbl(i) + lengbl(i)*(ijkpoin(i)-1)/nebl(i)
!             end if
!          !latest boundary layer
!          elseif((ijkpart(i)==npdir(i)).and.(ijkpoin(i)>(nedom(i)-nebl(i)+1))) then 
!             newijkpoin=ijkpoin(i)-(nedom(i)-nebl(i))
!             coord(i) = coor0(i) + leng(i) + lengbl(i)*(newijkpoin-1)/nebl(i)
!          else !core region
!             if(ntdis(i)==3) then
!                coord(i) = coor0(i) + leng(i)*((ijkpart(i)-1)*nedom(i)+(ijkpoin(i)-nebl(i)-1)/real(pdegr))/(nedir(i)-2*nebl(i))
!             elseif(ntdis(i)==4) then
!                coord(i) = coor0(i) - leng(i)*tanh(istret*(1.0_rp-2.0_rp*((ijkpart(i)-1)*nedom(i)+&
!                     (ijkpoin(i)-nebl(i)-1)/real(pdegr))/(nedir(i)-2*nebl(i))))/tanh(istret)
!             end if
!          end if
!       end if
!    end do
!   
!  end subroutine coord_ijk

!  !================================================================================================
!  subroutine check_part_boundary(ijk,npdir,ndime,isper)
!    !-----------------------------------------------------------------------
!    ! 
!    !-----------------------------------------------------------------------
!    implicit none
!    integer(ip), intent(in)    :: npdir(3),ndime,isper(3)
!    integer(ip), intent(inout) :: ijk(3)

!    integer(ip) :: i

!    ! Periodic mesh with respect to all directions
!    do i=1,ndime
!       if(isper(i)==1) then
!          if(ijk(i)<1) then
!             ijk(i) = npdir(i)
!          else if(ijk(i)>npdir(i)) then
!             ijk(i) = 1!npdir(i)
!          end if
!       else
!          if(ijk(i)<1.or.ijk(i)>npdir(i)) then
!             ijk = 0*ijk
!          end if
!       end if
!    end do

!  end subroutine check_part_boundary

!  !================================================================================================
!  subroutine check_face_boundary(ijkpart,pdime,iface,npdir,nsckt,npsoc,ndime,isper,auxv,flag, &
!       &                         neigh,lface)
!    !-----------------------------------------------------------------------
!    ! 
!    !-----------------------------------------------------------------------
!    implicit none
!    integer(ip), intent(in)  :: ijkpart(3),pdime,iface,npdir(3),ndime,isper(3),auxv(3,2)
!    integer(ip), intent(in)  :: nsckt(3),npsoc(3)
!    integer(ip), intent(out) :: flag,neigh(2),lface(2)

!    integer(ip)              :: i,ijkneig(3),newpo,cnt

!    cnt=1
!    flag = -1
!    lface = 0

!    do i=-1,0

!       ! Auxiliar variables
!       newpo = iface + i

!       ! Check neighbors
!       call neigh_face(newpo,pdime,auxv,ijkpart,ijkneig)

!       ! Check part boundary
!       call check_part_boundary(ijkneig,npdir,ndime,isper)

!       ! Global numbering
!       call globalid_2l(ijkneig,nsckt,npsoc,ndime,neigh(cnt))
!       
!       ! Set flag
!       if(neigh(cnt)>0) then
!          flag = min(flag+1,1)
!          lface(i+2) = 1 
!       end if

!       ! Update counter
!       cnt = cnt + 1

!    end do

!  end subroutine check_face_boundary

!  !================================================================================================
!  subroutine check_edge_boundary(ijkpart,pdime,ijedge,npdir,nsckt,npsoc,ndime,isper,auxv,flag, &
!       &                         neigh,ledge)
!    !-----------------------------------------------------------------------
!    ! 
!    !-----------------------------------------------------------------------
!    implicit none
!    integer(ip), intent(in)  :: ijkpart(3),pdime,ijedge(2),npdir(3),ndime,isper(3),auxv(:,:)
!    integer(ip), intent(in)  :: nsckt(3),npsoc(3)
!    integer(ip), intent(out) :: flag,neigh(2*(ndime-1)),ledge(2,2)

!    integer(ip)              :: i,j,ijkneig(3),newpo(2),cnt

!    cnt=1
!    flag = -1
!    ledge = 0

!    do i=-1,0
!       do j=-1,0

!          ! Auxiliar variables
!          newpo = ijedge + (/i,j/)

!          ! Check neighbors
!          call neigh_edge(newpo,pdime,auxv,ijkpart,ijkneig)

!          ! Check part boundary
!          call check_part_boundary(ijkneig,npdir,ndime,isper)

!          ! Global numbering
!          call globalid_2l(ijkneig,nsckt,npsoc,ndime,neigh(cnt))

!          ! Set flag
!          if(neigh(cnt)>0) then
!             flag = min(flag+1,1)
!             ledge(i+2,j+2) = 1
!          end if

!          ! Update counter
!          cnt = cnt + 1

!          if(ndime==2) exit
!       end do
!    end do

!  end subroutine check_edge_boundary

!  !================================================================================================
!  subroutine check_corn_boundary(ijkpart,ijkvert,npdir,nsckt,npsoc,ndime,isper,auxv,flag,neigh, &
!       &                         lcorn)
!    !-----------------------------------------------------------------------
!    ! 
!    !-----------------------------------------------------------------------
!    implicit none
!    integer(ip), intent(in)  :: ijkpart(3),ijkvert(3),npdir(3),ndime,isper(3),auxv(:,:)
!    integer(ip), intent(in)  :: nsckt(3),npsoc(3)
!    integer(ip), intent(out) :: flag,neigh(2**ndime),lcorn(2,2,2)

!    integer(ip)              :: i,j,k,ijkneig(3),newpo(3),cnt

!    cnt = 1
!    flag = -1
!    lcorn = 0

!    do i=-1,0
!       do j=-1,0
!          do k=-1,0

!             ! Auxiliar variables
!             newpo = ijkvert + (/i,j,k/)

!             ! Check neighbors
!             call neigh_corn(newpo,ijkpart,ijkneig)

!             ! Check part boundary
!             call check_part_boundary(ijkneig,npdir,ndime,isper)

!             ! Global numbering
!             call globalid_2l(ijkneig,nsckt,npsoc,ndime,neigh(cnt))

!             ! Set flag
!             if(neigh(cnt)>0) then
!                flag = min(flag+1,1)
!                lcorn(i+2,j+2,k+2) = 1
!             end if
!             ! Update counter
!             cnt = cnt + 1

!             if(ndime==2) exit
!          end do
!       end do
!    end do

!  end subroutine check_corn_boundary

!  !================================================================================================
!  subroutine neigh_face(newpo,pdime,auxv,ijkpart,ijkneig)
!    !-----------------------------------------------------------------------
!    ! 
!    !-----------------------------------------------------------------------
!    implicit none
!    integer(ip), intent(in)  :: ijkpart(3),newpo,pdime,auxv(3,2)
!    integer(ip), intent(out) :: ijkneig(3)

!    ijkneig(auxv(pdime,:)) = ijkpart(auxv(pdime,:))
!    ijkneig(pdime) = ijkpart(pdime) + newpo

!  end subroutine neigh_face

!  !================================================================================================
!  subroutine neigh_edge(newpo,pdime,auxv,ijkpart,ijkneig)
!    !-----------------------------------------------------------------------
!    ! 
!    !-----------------------------------------------------------------------
!    implicit none
!    integer(ip), intent(in)  :: ijkpart(3),newpo(2),pdime,auxv(:,:)
!    integer(ip), intent(out) :: ijkneig(3)

!    ijkneig(auxv(pdime,:)) = ijkpart(auxv(pdime,:)) + newpo(1:size(auxv(pdime,:)))
!    ijkneig(pdime) = ijkpart(pdime)

!  end subroutine neigh_edge

!  !================================================================================================
!  subroutine neigh_corn(newpo,ijkpart,ijkneig)
!    !-----------------------------------------------------------------------
!    ! 
!    !-----------------------------------------------------------------------
!    implicit none
!    integer(ip), intent(in)  :: ijkpart(3),newpo(3)
!    integer(ip), intent(out) :: ijkneig(3)

!    ijkneig = ijkpart + newpo

!  end subroutine neigh_corn

!  !==================================================================================================
!  subroutine generic_l2g(subgl,npnumg,npnumt,nenum,ndime,pdegr,gsize,tsize,gdata,l2ge,l2gp,trian, &
!       &                 coord)
!    !-----------------------------------------------------------------------
!    ! 
!    !-----------------------------------------------------------------------
!    implicit none
!    integer(ip)                    , intent(in)    :: subgl(ndime),ndime,pdegr
!    integer(ip)                    , intent(in)    :: npnumg(:),npnumt(:),nenum(:)
!    type(geom_size_t)              , intent(in)    :: gsize
!    type(topo_size_t)              , intent(in)    :: tsize
!    type(uniform_mesh_descriptor_t), intent(in)    :: gdata
!    integer(igp)                   , intent(out)   :: l2ge(:),l2gp(:)
!    type(triangulation_t)          , intent(inout) :: trian
!    real(rp)                       , intent(in)    :: coord(:,:)

!    integer(ip)               :: i,j,k,m,n,l,ijk(3),ijklnods(3),num,count,auxva((pdegr+1)**ndime)
!    integer(ip)               :: ne_aux(3),np_aux(3),mcase,auxnum,nedir_l2g(ndime),nvef,aux_glb
!    integer(ip)               :: subgl_aux(3),subgl_ijk(3),nelbl(3),ncorn
!    integer(ip)               :: cnt,ipoin,jpoin,ielem,pdime,iedge,jedge,iface,nddomk,idime
!    integer(ip)               :: auxvc(2**ndime),auxvf(2*ndime),auxvd(2*(ndime-1)*ndime),aux_cnt
!    integer(ip) , allocatable :: auxv(:,:)
!    integer(igp), allocatable :: po_l2g(:),el_l2g(:)

!    ! Allocate
!    call memalloc(gsize%nedomt,el_l2g,__FILE__,__LINE__)
!    call memalloc(tsize%notot ,po_l2g,__FILE__,__LINE__)

!    ! Auxiliar vectro
!    check( pdegr == 1 )  ! Only working for linear geometrical elements
!    if(ndime==2) then
!       call memalloc(2,1,auxv,__FILE__,__LINE__)
!       auxv(1,1) = 2
!       auxv(2,1) = 1
!       if(pdegr==1) then
!          auxva = (/1,3,2,4/)
!       else if(pdegr==2) then
!          auxva = (/1,8,4,5,9,7,2,6,3/)
!       else if(pdegr==3) then
!          auxva = (/1,12,11,4,5,13,16,10,6,14,15,9,2,7,8,3/)
!       end if
!       ne_aux(1)=1
!       ne_aux(2)=gsize%nedom(1)
!       ne_aux(3)=gsize%nedom(2)
!       np_aux(1)=1
!       np_aux(2)=gsize%npdom(1)
!       np_aux(3)=gsize%npdom(2)
!       subgl_aux(1:2) = subgl
!       subgl_aux(3) = 0
!       auxvc = (/1,3,2,4/)
!       auxvd = (/7,8,5,6/)
!       nvef = 8
!       ncorn = 4
!    else if(ndime==3) then
!       call memalloc(3,2,auxv,__FILE__,__LINE__)
!       auxv(1,:) = (/2,3/)
!       auxv(2,:) = (/1,3/)
!       auxv(3,:) = (/1,2/)
!       if(pdegr==1) then
!          auxva = (/1,5,3,7,2,6,4,8/)
!       else if(pdegr==2) then
!          auxva = (/1,13,5,12,25,20,4,16,8,9,22,17,21,27,26,11,24,19,2,14,6,10,23,18,3,15,7/)
!       else if(pdegr==3) then
!          auxva = (/1,17,21,5,16,44,52,32,15,43,51,31,4,20,24,8,9,37,45,25,33,57,61,53,36,60,64,56,14,42, &
!               &    50,30,10,38,46,26,34,58,62,54,35,59,63,55,13,41,49,29,2,18,22,6,11,39,47,27,12,40,48, &
!               &    28,3,19,23,7/)
!       end if
!       ne_aux=gsize%nedom
!       np_aux=gsize%npdom
!       subgl_aux = subgl
!       auxvc = (/1,5,3,7,2,6,4,8/)
!       auxvd = (/17,19,18,20,13,15,14,16,9,11,10,12/)
!       auxvf = (/21,22,23,24,25,26/)
!       nvef = 26
!       ncorn = 8
!    end if
!    mcase = gdata%mater
!    nelbl(1) = gdata%neblx
!    nelbl(2) = gdata%nebly
!    nelbl(3) = gdata%neblz

!    ! Elements
!    do i=1,ne_aux(1)
!       do j=1,ne_aux(2)
!          do k=1,ne_aux(3)

!             ! Identifier
!             if(ndime==2) then
!                ijk=(/j,k,i/)
!             else if(ndime==3) then
!                ijk=(/i,j,k/)
!             end if
!             call globalid(ijk,gsize%nedom,ndime,num)
!             
!             ! Local to global vector
!             subgl_ijk = subgl_aux+ijk
!             call globalid(subgl_ijk,gsize%nedir,ndime,el_l2g(num))

!             ! Fill material
!             if(mcase>0) then
!                call materialid(mcase,subgl_ijk,gsize%nedir,nelbl,trian%elems(nenum(num))%subset_id)
!             end if

!             ! Allocate triangulation elemental objects
!             trian%elems(nenum(num))%num_vefs = nvef
!             call memalloc(trian%elems(nenum(num))%num_vefs,trian%elems(nenum(num))%vefs, &
!                  &        __FILE__,__LINE__)
!             call memalloc(trian%num_dims,ncorn,trian%elems(nenum(num))%coordinates, __FILE__, __LINE__ )
!             call put_topology_element_triangulation(nenum(num),trian)           
!             count = 1
!             ! Elemental corners
!             do m=0,1
!                do n=0,1
!                   do l=0,1

!                      ! Corner identifier in the element
!                      ijklnods= ijk + (/m,n,l/)
!                      call globalid(ijklnods,gsize%npdom,ndime,auxnum)

!                      ! Generate lnods
!                      trian%elems(nenum(num))%vefs(auxvc(count)) = npnumt(auxnum)

!                      ! Generate coords
!                      trian%elems(nenum(num))%coordinates(:,auxva(count)) = coord(:,npnumg(auxnum))

!                      ! Update counter
!                      count = count +1

!                      if(ndime==2) exit
!                   end do
!                end do
!             end do
!                
!             ! Elemental edges
!             count = 1
!             aux_cnt = tsize%nctot
!             do pdime=ndime,1,-1
!                do iedge=0,1
!                   do jedge=0,1

!                      ! Edge identifier in the element
!                      ijklnods(pdime) = ijk(pdime)
!                      ijklnods(auxv(pdime,1)) = ijk(auxv(pdime,1)) + iedge
!                      if(ndime==3) ijklnods(auxv(pdime,2)) = ijk(auxv(pdime,2)) + jedge
!                      call globalid(ijklnods,tsize%nddom(:,pdime),ndime,auxnum)

!                      ! Generate lnods
!                      trian%elems(nenum(num))%vefs(auxvd(count)) = npnumt(aux_cnt + auxnum)

!                      ! Update counter
!                      count = count +1

!                      if(ndime==2) exit
!                   end do
!                end do
!                aux_cnt = aux_cnt + tsize%nddir(pdime)
!             end do

!             ! Elemental faces
!             count = 1
!             if(ndime==3) then
!                aux_cnt = tsize%nctot + tsize%ndtot
!                do pdime=ndime,1,-1
!                   do iface = 0,1

!                      ! Face identifier in the element
!                      ijklnods(pdime) = ijk(pdime) + iface
!                      ijklnods(auxv(pdime,:)) = ijk(auxv(pdime,:))
!                      call globalid(ijklnods,tsize%nfdom(:,pdime),ndime,auxnum)

!                      ! Generate lnods
!                      trian%elems(nenum(num))%vefs(auxvf(count)) = npnumt(aux_cnt + auxnum)

!                      ! Update counter
!                      count = count + 1

!                   end do
!                   aux_cnt = aux_cnt + tsize%nfdir(pdime)
!                end do
!             end if

!          end do
!       end do
!       if(ndime==2) exit
!    end do

!    ! Construct l2g emap (ordered by objects)
!    do ielem=1,gsize%nedomt
!       l2ge(nenum(ielem)) = el_l2g(ielem)
!    end do

!    ! Nodes
!    do i=1,np_aux(1)
!       do j=1,np_aux(2)
!          do k=1,np_aux(3)

!             ! Identifier
!             if(ndime==2) then
!                ijk=(/j,k,i/)
!             else if(ndime==3) then
!                ijk=(/i,j,k/)
!             end if
!             call globalid(ijk,gsize%npdom,ndime,num)

!             ! Local to global vector
!             nedir_l2g = gsize%nedir(1:ndime)*pdegr + 1
!             subgl_ijk = subgl_aux*pdegr+ijk
!             do idime=1,ndime
!                ! Modify global ID in periodic boundary
!                if(subgl_ijk(idime)==(gsize%nedir(idime)+1).and.(gdata%isper(idime)==1)) subgl_ijk(idime) = 1
!             end do
!             call globalid(subgl_ijk,nedir_l2g,ndime,po_l2g(num))
!    
!          end do
!       end do
!       if(ndime==2) exit
!    end do

!    ! Edges
!    aux_cnt = tsize%nctot
!    aux_glb = tsize%ncglb
!    do pdime=ndime,1,-1
!       if(ndime==2) then
!          nddomk=1
!       elseif(ndime==3) then
!          nddomk=tsize%nddom(auxv(pdime,2),pdime)
!       end if
!       do i=1,tsize%nddom(pdime,pdime)
!          do j=1,tsize%nddom(auxv(pdime,1),pdime)
!             do k=1,nddomk

!                ! Identifier
!                ijk(pdime) = i
!                ijk(auxv(pdime,1)) = j
!                if(ndime==3) ijk(auxv(pdime,2)) = k
!                call globalid(ijk,tsize%nddom(:,pdime),ndime,num)

!                ! Local to global vector
!                subgl_ijk = subgl_aux + ijk
!                do idime=1,ndime-1
!                   ! Modify global ID in periodic boundary
!                   if(subgl_ijk(auxv(pdime,idime))==(gsize%nedir(auxv(pdime,idime))+1).and. &
!                        & (gdata%isper(auxv(pdime,idime))==1)) subgl_ijk(auxv(pdime,idime)) = 1
!                end do
!                call globalid(subgl_ijk,tsize%ndglb(:,pdime),ndime,po_l2g(num+aux_cnt))
!                po_l2g(num+aux_cnt) = po_l2g(num+aux_cnt) + aux_glb

!             end do
!          end do
!       end do
!       aux_cnt = aux_cnt + tsize%nddir(pdime)
!       aux_glb = aux_glb + tsize%ndsum(pdime)
!    end do

!    ! Faces
!    if(ndime==3) then
!       aux_cnt = tsize%nctot + tsize%ndtot
!       do pdime=ndime,1,-1
!          do i=1,tsize%nfdom(pdime,pdime)
!             do j=1,tsize%nfdom(auxv(pdime,1),pdime)
!                do k=1,tsize%nfdom(auxv(pdime,2),pdime)

!                   ! Identifier
!                   ijk(pdime) = i
!                   ijk(auxv(pdime,:)) = (/j,k/)
!                   call globalid(ijk,tsize%nfdom(:,pdime),ndime,num)

!                   ! Local to global vector
!                   subgl_ijk = subgl_aux + ijk
!                   ! Modify global ID in periodic boundary
!                   if(subgl_ijk(pdime)==(gsize%nedir(pdime)+1).and. &
!                        & (gdata%isper(pdime)==1)) subgl_ijk(pdime) = 1
!                   call globalid(subgl_ijk,tsize%nfglb(:,pdime),ndime,po_l2g(num+aux_cnt))
!                   po_l2g(num+aux_cnt) = po_l2g(num+aux_cnt) + aux_glb

!                end do
!             end do
!          end do
!          aux_cnt = aux_cnt + tsize%nfdir(pdime)
!          aux_glb = aux_glb + tsize%nfsum(pdime)
!       end do
!    end if

!    ! Construct l2g nmap (ordered by objects)
!    ! Elemental corners
!    do ipoin=1,gsize%npdomt
!       !l2gp(npnumg(ipoin)) = po_l2g(ipoin
!       l2gp(npnumt(ipoin)) = po_l2g(ipoin)
!    end do
!    ! Elemental edges
!    aux_cnt = tsize%nctot
!    do pdime=ndime,1,-1
!       do i=1,tsize%nddir(pdime)
!          l2gp(npnumt(i+aux_cnt)) = po_l2g(i+aux_cnt)
!       end do
!       aux_cnt = aux_cnt + tsize%nddir(pdime)
!    end do
!    ! Elemental faces
!    aux_cnt = tsize%nctot + tsize%ndtot
!    do pdime=ndime,1,-1
!       do i=1,tsize%nfdir(pdime)
!          l2gp(npnumt(i+aux_cnt)) = po_l2g(i+aux_cnt)
!       end do
!       aux_cnt = aux_cnt + tsize%nfdir(pdime)
!    end do

!    ! Deallocate
!    call memfree(el_l2g,__FILE__,__LINE__)
!    call memfree(po_l2g,__FILE__,__LINE__)
!    call memfree(auxv,__FILE__,__LINE__)

!  end subroutine generic_l2g

!  !==================================================================================================
!  subroutine global_to_ijk(lpart,nsckt,npsoc,ndime,ijk)
!    implicit none
!    integer(ip), intent(in)  :: lpart,ndime,nsckt(ndime),npsoc(ndime)
!    integer(ip), intent(out) :: ijk(3)

!    integer(ip) :: isckt,jsckt,ksckt,lsckt
!    integer(ip) :: ipart_aux,jpart_aux,kpart_aux,lpart_aux
!    integer(ip) :: ipart,jpart,kpart
!    real(rp) :: aux1,aux2,aux3,aux4,aux5,aux6,aux7

!    if(ndime==2) then
!       ! First level (Socket level)
!       aux1 = (lpart-1)/(npsoc(1)*npsoc(2))
!       lsckt = floor(aux1) + 1
!       aux2 = (lsckt-1)/nsckt(2)
!       jsckt = lsckt - floor(aux2)*nsckt(2)
!       isckt = floor(aux2) + 1

!       ! Second level (Part inside the socket)
!       lpart_aux = lpart - (lsckt-1)*npsoc(1)*npsoc(2)
!       aux3 = (lpart_aux-1)/npsoc(2)
!       jpart_aux = lpart_aux - floor(aux3)*npsoc(2)
!       ipart_aux = floor(aux3) + 1

!       ! ijk numeration
!       ipart = ipart_aux + (isckt-1)*npsoc(1)
!       jpart = jpart_aux + (jsckt-1)*npsoc(2)
!       kpart = 1
!       ijk = (/ipart,jpart,kpart/)
!    else if(ndime==3) then
!       ! First level (Socket level)
!       aux1 = (lpart-1)/(npsoc(1)*npsoc(2)*npsoc(3))
!       lsckt = floor(aux1) + 1
!       aux2 = (lsckt-1)/(nsckt(2)*nsckt(3))
!       isckt = floor(aux2) + 1
!       aux3 = (lsckt - (isckt-1)*nsckt(2)*nsckt(3) - 1)/nsckt(3)
!       jsckt = floor(aux3) + 1
!       aux4 = lsckt - (isckt-1)*nsckt(2)*nsckt(3) - (jsckt-1)*nsckt(3) - 1
!       ksckt = floor(aux4) + 1

!       ! Second level (Part inside the socket)
!       lpart_aux = lpart - (lsckt-1)*npsoc(1)*npsoc(2)*npsoc(3)
!       aux5 = (lpart_aux-1)/(npsoc(2)*npsoc(3))
!       ipart_aux = floor(aux5) + 1
!       aux6 = (lpart_aux - (ipart_aux-1)*npsoc(2)*npsoc(3) - 1)/npsoc(3)
!       jpart_aux = floor(aux6) + 1
!       aux7 = lpart_aux - (ipart_aux-1)*npsoc(2)*npsoc(3) - (jpart_aux-1)*npsoc(3) - 1
!       kpart_aux = floor(aux7) + 1

!       ! ijk numeration
!       ipart = ipart_aux + (isckt-1)*npsoc(1)
!       jpart = jpart_aux + (jsckt-1)*npsoc(2)
!       kpart = kpart_aux + (ksckt-1)*npsoc(3)
!       ijk = (/ipart,jpart,kpart/)       
!    end if
!  
!  end subroutine global_to_ijk

!  !================================================================================================
!  subroutine aux_surf_num(pdime,aux_surf,gsurf_aux)
!    implicit none
!    integer(ip), intent(in)  :: pdime
!    integer(ip), intent(out) :: aux_surf(4,2),gsurf_aux(4)

!    if(pdime==3) then
!       gsurf_aux = (/3,4,5,6/)
!    else if(pdime==2) then
!       gsurf_aux = (/1,2,5,6/)
!    else if(pdime==1) then
!       gsurf_aux = (/1,2,3,4/)
!    end if
!    aux_surf(1,:) = (/2,4/)
!    aux_surf(2,:) = (/1,3/)
!    aux_surf(3,:) = (/3,4/)
!    aux_surf(4,:) = (/1,2/)    

!  end subroutine aux_surf_num

!  !==================================================================================================
!  subroutine kronec(i,n,kron)
!    !-----------------------------------------------------------------------
!    ! 
!    !-----------------------------------------------------------------------
!    implicit none
!    integer(ip), intent(in)  :: i,n
!    integer(ip), intent(out) :: kron(n)

!    kron = 0
!    kron(i) = 1

!  end subroutine kronec

!  !==================================================================================================
!  subroutine globalid_ip(ijk,nd,ndime,gl)
!    implicit none
!    integer(ip), intent(in)  :: ijk(ndime),nd(ndime),ndime
!    integer(ip), intent(out) :: gl

!    if(ndime==1) then
!       gl = ijk(1)
!    elseif(ndime==2) then
!       gl = (ijk(1)-1)*nd(2)+ijk(2)
!    else if(ndime==3) then
!       gl = (ijk(1)-1)*nd(2)*nd(3)+(ijk(2)-1)*nd(3)+ijk(3)
!    end if

!  end subroutine globalid_ip

!  !==================================================================================================
!  subroutine globalid_igp(ijk,nd,ndime,gl)
!    implicit none
!    integer(ip), intent(in)  :: ijk(ndime),nd(ndime),ndime
!    integer(igp), intent(out) :: gl
!    
!    if(ndime==1) then
!       gl = int(ijk(1),igp)
!    elseif(ndime==2) then
!       gl = (int(ijk(1),igp)-1)*int(nd(2),igp)+int(ijk(2),igp)
!    else if(ndime==3) then
!       gl = (int(ijk(1),igp)-1)*int(nd(2),igp)*int(nd(3),igp)+(int(ijk(2),igp)-1)*int(nd(3),igp)+int(ijk(3),igp)
!    end if
!    
!  end subroutine globalid_igp
!  
!  !==================================================================================================
!  subroutine globalid_2l(ijk,nsckt,npsoc,ndime,gl_2l)
!    implicit none
!    integer(ip), intent(in)  :: ijk(3),nsckt(3),npsoc(3),ndime
!    integer(ip), intent(out) :: gl_2l

!    integer(ip) :: aux1(ndime),aux2(ndime),aux3,aux4,i,gl,local_ijk(ndime)
!    real(rp)    :: work

!    ! Auxiliar variables
!    aux3=1
!    do i=1,ndime
!       work = (ijk(i)-1)/npsoc(i)
!       aux1(i) = floor(work) + 1;
!       aux2(i) = floor(work)*npsoc(i)
!       aux3 = aux3*npsoc(i)
!    end do

!    ! 2-level global numeration
!    call globalid(aux1,nsckt,ndime,aux4)
!    aux4 = (aux4-1)*aux3
!    local_ijk = ijk(1:ndime) - aux2
!    call globalid(local_ijk,npsoc,ndime,gl)
!    gl_2l = aux4 + gl
!    
!  end subroutine globalid_2l

!  !==================================================================================================
!  subroutine materialid(mcase,ijkelem,nedir,nelbl,mater)
!    implicit none
!    integer(ip), intent(in)    :: mcase,ijkelem(3),nedir(3),nelbl(3)
!    integer(ip), intent(inout) :: mater
!    
!    if(mcase==1) then
!       ! (x<leng(x)/2 --> mat == 1)
!       ! (x>leng(x)/2 --> mat == 2)
!       assert(mod(nedir(1),2)==0)
!       if(ijkelem(1).le.nedir(1)/2) then
!          mater = 1
!       else
!          mater = 2
!       end if
!    elseif(mcase==2) then
!       ! (y<leng(y)/2 --> mat == 1)
!       ! (y>leng(y)/2 --> mat == 2)
!       assert(mod(nedir(2),2)==0)
!       if(ijkelem(2).le.nedir(2)/2) then
!          mater = 1
!       else
!          mater = 2
!       end if
!    elseif(mcase==3) then   ! Hunt's case
!       ! (z >= -1.0 and z <= 1.0 --> mat == 1) IMH
!       ! (z < -1.0 or z > 1.0 --> mat == 2) DCY
!       if(ijkelem(3)>nelbl(3) .and. ijkelem(3)<=nedir(3)-nelbl(3)) then  ! We have nelbl elements in the layer (DCY-solid)
!          mater = 1
!       else
!          mater = 2
!       end if
!    elseif(mcase==4) then   ! Backward-facing step
!       ! (x >= 0.0 or y >= 0.0 --> mat == 1) NSI
!       ! (x < 0.0 and y < 0.0 --> mat == 2) Boundary
!       if(ijkelem(1)>nelbl(1) .or. ijkelem(2)>nelbl(2)) then 
!          mater = 1
!       else
!          mater = 2
!       end if
!    end if
!    
!  end subroutine materialid
  
end module uniform_hex_mesh_generator_names
