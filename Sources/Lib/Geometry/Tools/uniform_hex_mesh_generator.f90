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
  use stdio_names
  use reference_fe_names
  use uniform_hex_mesh_generator_parameters_names
  use FPL
  implicit none
# include "debug.i90"
  private 

  type uniform_hex_mesh_t 
     private 
     integer(ip)              :: num_dims
     integer(ip)              :: interpolation_order
     integer(ip), allocatable :: num_cells_x_dim(:)
     integer(ip), allocatable :: num_parts_x_dim_x_level(:)
     integer(ip), allocatable :: is_dir_periodic(:)
     real(rp),    allocatable :: domain_limits(:)
   contains   
     procedure         , non_overridable :: create                                         => uhm_create
     procedure, private, non_overridable :: process_parameters                             => uhm_process_parameters
     procedure         , non_overridable :: generate_part_local_mesh                       => uhm_generate_part_local_mesh
     procedure, private, non_overridable :: check_part_local_mesh_parameters               => uhm_check_part_local_mesh_parameters
     procedure         , non_overridable :: generate_part_aggregation_among_levels         => uhm_generate_part_aggregation_among_levels
     procedure, private, non_overridable :: check_part_aggregation_among_levels_parameters => uhm_check_part_aggregation_among_levels_parameters 
     procedure         , non_overridable :: free                                           => uhm_free
     procedure         , non_overridable :: get_num_dims                                   => uhm_get_num_dims
  end type uniform_hex_mesh_t
  
  interface ijk_to_spatial_numbering
     module procedure ijk_to_spatial_numbering_ip
  end interface ijk_to_spatial_numbering

  public :: uniform_hex_mesh_t

contains

  subroutine uhm_create(this, parameters)
    implicit none
    class(uniform_hex_mesh_t), intent(inout) :: this
    type(ParameterList_t)    , intent(in)    :: parameters
    call this%free()
    call this%process_parameters(parameters)
  end subroutine uhm_create

  subroutine uhm_free(this)
    implicit none
    class(uniform_hex_mesh_t), intent(inout) :: this
    if(allocated(this%num_cells_x_dim)) call memfree(this%num_cells_x_dim,__FILE__,__LINE__)
    if(allocated(this%num_parts_x_dim_x_level)) call memfree(this%num_parts_x_dim_x_level,__FILE__,__LINE__)
    if(allocated(this%is_dir_periodic)) call memfree(this%is_dir_periodic,__FILE__,__LINE__)
    if(allocated(this%domain_limits)) call memfree(this%domain_limits,__FILE__,__LINE__)
  end subroutine uhm_free
  
  pure function uhm_get_num_dims(this)
    implicit none
    class(uniform_hex_mesh_t), intent(in) :: this
    integer(ip) :: uhm_get_num_dims
    uhm_get_num_dims = this%num_dims
  end function uhm_get_num_dims

  subroutine uhm_process_parameters(this,parameters)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine generates geometry data to construct a structured mesh                      !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(uniform_hex_mesh_t), intent(inout) :: this
    type(ParameterList_t)    , intent(in)    :: parameters
    ! Locals
    integer(ip) :: istat, idime
    logical     :: is_present
    integer(ip), allocatable :: array_size(:)
    real(rp), allocatable :: domain_limits(:)
    
    ! Mandatory
    assert(parameters%isPresent(struct_hex_mesh_generator_num_dims_key))
    assert(parameters%isAssignable(struct_hex_mesh_generator_num_dims_key, this%num_dims))
    istat = parameters%get(key = struct_hex_mesh_generator_num_dims_key, value = this%num_dims); check(istat==0)
    
    ! Mandatory (array)
    assert ( parameters%isPresent(key = struct_hex_mesh_generator_num_cells_x_dim_key ) )
    istat = parameters%GetShape(key = struct_hex_mesh_generator_num_cells_x_dim_key, shape = array_size); check(istat==0)
    call memalloc(array_size(1), this%num_cells_x_dim,__FILE__,__LINE__)
    assert(parameters%isAssignable(struct_hex_mesh_generator_num_cells_x_dim_key, this%num_cells_x_dim))
    istat = parameters%get(key = struct_hex_mesh_generator_num_cells_x_dim_key, value = this%num_cells_x_dim); check(istat==0)

    ! Mandatory (array)
    assert (parameters%isPresent(key = struct_hex_mesh_generator_is_dir_periodic_key ) )
    istat = parameters%GetShape(key = struct_hex_mesh_generator_is_dir_periodic_key, shape = array_size); check(istat==0);
    call memalloc(array_size(1), this%is_dir_periodic,__FILE__,__LINE__)
    assert(parameters%isAssignable(struct_hex_mesh_generator_is_dir_periodic_key, this%is_dir_periodic))
    istat = parameters%get(key = struct_hex_mesh_generator_is_dir_periodic_key, value = this%is_dir_periodic); check(istat==0)

    ! Mandatory (array)
    assert( parameters%isPresent(key = struct_hex_mesh_generator_num_parts_x_dim_x_level_key) )
    istat = parameters%GetShape(key = struct_hex_mesh_generator_num_parts_x_dim_x_level_key, shape = array_size); check(istat==0)
    call memalloc(array_size(1), this%num_parts_x_dim_x_level,__FILE__,__LINE__)
    assert(parameters%isAssignable(struct_hex_mesh_generator_num_parts_x_dim_x_level_key, this%num_parts_x_dim_x_level))
    istat = parameters%get(key = struct_hex_mesh_generator_num_parts_x_dim_x_level_key , value = this%num_parts_x_dim_x_level); check(istat==0)
    
    ! Mandatory (array)
    assert( parameters%isPresent(key = struct_hex_mesh_generator_domain_limits_key) )
    istat = parameters%GetShape(key = struct_hex_mesh_generator_domain_limits_key, shape = array_size); check(istat==0)
    call memalloc(array_size(1), this%domain_limits,__FILE__,__LINE__)
    assert(parameters%isAssignable(struct_hex_mesh_generator_domain_limits_key, this%domain_limits))
    istat = parameters%get(key = struct_hex_mesh_generator_domain_limits_key , value = this%domain_limits); check(istat==0)

    ! Here we do not use our memfree because array_size was allocated inside FPL 
    ! (without calling memalloc)
    if(allocated(array_size)) deallocate(array_size)
  end subroutine uhm_process_parameters
    
  subroutine uhm_generate_part_local_mesh(this, &
                                          num_local_cells, &
                                          num_local_vefs, &
                                          num_vertices, &
                                          num_edges, &
                                          num_faces, &
                                          ptr_vefs_x_cell, &
                                          lst_vefs_lids, &
                                          boundary_id, &
                                          coordinates, &
                                          num_ghost_cells, &
                                          cells_gids, &
                                          cells_mypart, &
                                          vefs_gids, &
                                          num_itfc_cells, &
                                          lst_itfc_cells, &
                                          ptr_ext_neighs_x_itfc_cell, &
                                          lst_ext_neighs_gids, &
                                          lst_ext_neighs_part_ids, &
                                          part_id)
    implicit none
    class(uniform_hex_mesh_t) , intent(in)    :: this
    integer(ip)               , intent(out)   :: num_local_cells
    integer(ip)               , intent(out)   :: num_local_vefs
    integer(ip)               , intent(out)   :: num_vertices
    integer(ip)               , intent(out)   :: num_edges
    integer(ip)               , intent(out)   :: num_faces
    integer(ip)  , allocatable, intent(inout) :: ptr_vefs_x_cell(:)            ! Size = num_local_cells + 1
    integer(ip)  , allocatable, intent(inout) :: lst_vefs_lids(:)                ! Size = ptr_vefs_x_cell(num_local_cells+1)-1
    integer(ip)  , allocatable, intent(inout) :: boundary_id(:)                  ! Size = num_local_vefs
    real(rp)     , allocatable, intent(inout) :: coordinates(:,:)

    integer(ip)               , optional, intent(out)   :: num_ghost_cells
    integer(ip)               , optional, intent(out)   :: num_itfc_cells
    integer(igp) , allocatable, optional, intent(inout) :: cells_gids(:)                   ! Size = num_local_cells + num_ghost_cells
    integer(ip)  , allocatable, optional, intent(inout) :: cells_mypart(:)                 ! Size = num_local_cells + num_ghost_cells
    integer(igp) , allocatable, optional, intent(inout) :: vefs_gids(:)                    ! Size = num_local_vefs
    integer(ip)  , allocatable, optional, intent(inout) :: lst_itfc_cells(:)              
    integer(ip)  , allocatable, optional, intent(inout) :: ptr_ext_neighs_x_itfc_cell(:)
    integer(igp) , allocatable, optional, intent(inout) :: lst_ext_neighs_gids(:)         
    integer(ip)  , allocatable, optional, intent(inout) :: lst_ext_neighs_part_ids(:)
    integer(ip)               , optional, intent(in)    :: part_id

    integer(ip), allocatable :: cell_permutation(:)

    integer(ip) :: part_ijk(0:SPACE_DIM-1)
    integer(ip) :: mypart_ijk(0:SPACE_DIM-1)
    integer(ip) :: cell_ijk(0:SPACE_DIM-1)
    integer(ip) :: neighbor_ijk(0:SPACE_DIM-1)
    integer(ip) :: neighbor_part_ijk(0:SPACE_DIM-1)
    integer(ip) :: nface_ijk(0:SPACE_DIM-1)
    integer(ip) :: first_cell_ijk(0:SPACE_DIM-1)

    ! Here total=local+ghost (if any) refers to the things I have,
    ! whereas global to the whole distributed mesh.
    integer(ip) :: num_total_cells_x_dir(0:SPACE_DIM-1)
    integer(ip) :: num_local_cells_x_dir(0:SPACE_DIM-1)
    integer(ip) :: num_left_parts_x_dir(0:SPACE_DIM-1)
    integer(ip) :: num_right_parts_x_dir(0:SPACE_DIM-1)

    integer(igp), allocatable  :: num_global_n_faces(:)
    integer(ip) , allocatable  :: num_total_n_faces(:)
    integer(ip) , allocatable  :: num_global_nfaces_x_dir(:,:)
    integer(ip) , allocatable  :: num_total_nfaces_x_dir(:,:)

    integer(ip)               :: topology, num_nface_types, partial_count
    integer(ip)               :: ighost_cell, ilocal_cell
    integer(ip)               :: num_ghost_cells_ , has_left_ghost, has_right_ghost
    integer(ip)               :: idime, jdime, icell, iface, iface_of_itype, index, itype, itfc_cells

    type(polytope_t)     :: polytope
    !type(node_array_t)        :: node_array
    integer(ip)               :: ones(SPACE_DIM)
    logical                   :: count_it

    integer(ip)               :: num_dims
    integer(ip)               :: num_cells_x_dim(0:SPACE_DIM-1)
    integer(ip)               :: is_dir_periodic(0:SPACE_DIM-1)
    real(rp)                  :: domain_limits(1:SPACE_DIM,2)
    integer(ip)               :: num_parts_x_dim_x_level(0:SPACE_DIM-1)

    if(present(num_ghost_cells)) then
       assert(present(num_itfc_cells))
       assert(present(cells_gids))
       assert(present(cells_mypart))
       assert(present(vefs_gids))
       assert(present(lst_itfc_cells))
       assert(present(ptr_ext_neighs_x_itfc_cell))
       assert(present(lst_ext_neighs_gids))  
       assert(present(lst_ext_neighs_part_ids))
    end if
    
    call this%check_part_local_mesh_parameters()

    ! Copy member to local variables with lowerbound 0
    num_dims                               = this%num_dims
    num_cells_x_dim(0:SPACE_DIM-1)         = this%num_cells_x_dim(1:SPACE_DIM)
    is_dir_periodic(0:SPACE_DIM-1)         = this%is_dir_periodic(1:SPACE_DIM)
    domain_limits(0:SPACE_DIM-1)           = this%domain_limits(1:SPACE_DIM)
    num_parts_x_dim_x_level(0:SPACE_DIM-1) = this%num_parts_x_dim_x_level(1:SPACE_DIM)

    ones = 1
    topology = 2**num_dims-1  ! Hexahedral
    call polytope%create( num_dims, topology ) 
    !call node_array%create ( polytope, ones*this%interpolation_order )
    
    ! PARTS
    ! =====
    ! Get my part coordinates (make it 0-based, assuming part_id is 1-based) and the number of parts I have around (if any)
    if(present(part_id)) then
       call spatial_to_ijk_numbering(num_dims, num_parts_x_dim_x_level, part_id, part_ijk)
    else
       part_ijk = 0
       num_parts_x_dim_x_level = 1
    end if
    num_left_parts_x_dir=1
    num_right_parts_x_dir=1
    do idime = 0, num_dims - 1 
       if(is_dir_periodic(idime)==0 .or. num_parts_x_dim_x_level(idime)==1) then ! Not periodic
          if(part_ijk(idime)==0) num_left_parts_x_dir(idime)=0 
          if(part_ijk(idime)==num_parts_x_dim_x_level(idime)-1) num_right_parts_x_dir(idime)=0 
       end if
    end do

    ! CELLS
    ! =====
    ! Global and local number of cells (per direction and total; local, ghost and global)
    do idime = 0, num_dims - 1 
       num_local_cells_x_dir(idime) = num_cells_x_dim(idime) / num_parts_x_dim_x_level(idime) 
       first_cell_ijk(idime) =  part_ijk(idime) * num_local_cells_x_dir(idime)
    end do
    num_local_cells = 1
    do idime = 0, num_dims - 1
       num_local_cells = num_local_cells * num_local_cells_x_dir(idime)
    end do
    num_total_cells_x_dir = num_local_cells_x_dir + num_left_parts_x_dir + num_right_parts_x_dir
    first_cell_ijk = first_cell_ijk - num_left_parts_x_dir

    if(present(num_ghost_cells)) then

       ! Count ghost cells
       num_ghost_cells_ = 1
       do idime = 0, num_dims - 1
          num_ghost_cells_ = num_ghost_cells_ * num_total_cells_x_dir(idime)
       end do
       num_ghost_cells_ = num_ghost_cells_ - num_local_cells
       num_ghost_cells  = num_ghost_cells_

       ! Number of interface cells (=local-interior)
       num_itfc_cells = 1
       do idime = 0, num_dims - 1
          if(num_local_cells_x_dir(idime) > &
               & num_left_parts_x_dir(idime) + num_right_parts_x_dir(idime)) then
             num_itfc_cells = num_itfc_cells * &
                  & (num_local_cells_x_dir(idime) - num_left_parts_x_dir(idime) - num_right_parts_x_dir(idime))
          else
             num_itfc_cells = 0
             exit
          end if
       end do
       num_itfc_cells  = num_local_cells - num_itfc_cells

       ! Generate permutation vector to store ghost cells at the end
       ! List ghost cells (a permutation array to reorder cells with ghost at the end)
       call memalloc(num_local_cells+num_ghost_cells_, cell_permutation, __FILE__,__LINE__)
       cell_permutation=0
       ighost_cell = 0
       ilocal_cell = 0
       itfc_cells = 0
       do icell = 1, num_local_cells+num_ghost_cells_
          call spatial_to_ijk_numbering(num_dims, num_total_cells_x_dir, icell, cell_ijk)
          do idime = 0, num_dims - 1
             if(is_dir_periodic(idime)==0) then
                if(    (num_left_parts_x_dir(idime)==1 .and.cell_ijk(idime)==0).or. &
                     & (num_right_parts_x_dir(idime)==1.and.cell_ijk(idime)==num_total_cells_x_dir(idime)-1) ) then ! cell is ghost
                   cell_permutation(icell) = num_local_cells + num_ghost_cells_- ighost_cell
                   ighost_cell = ighost_cell + 1
                   exit
                end if
             end if
          end do
          if(cell_permutation(icell)==0) then
             do idime = 0, num_dims - 1
                if(is_dir_periodic(idime)==0) then
                   if((num_left_parts_x_dir(idime)==1 .and.cell_ijk(idime)==1).or. &
                        &  (num_right_parts_x_dir(idime)==1.and.cell_ijk(idime)==num_total_cells_x_dir(idime)-2) ) then ! cell is interface
                      cell_permutation(icell) = num_local_cells - itfc_cells
                      itfc_cells = itfc_cells + 1
                      exit
                   end if
                end if
             end do
          end if
          if(cell_permutation(icell)==0) then
             ilocal_cell = ilocal_cell + 1
             cell_permutation(icell) = ilocal_cell
          end if
       end do
       assert(ilocal_cell == num_local_cells-num_itfc_cells)
       assert(itfc_cells == num_itfc_cells)
       assert(ighost_cell == num_ghost_cells_)

    else

       num_ghost_cells_ = 0
       call memalloc(num_local_cells, cell_permutation, __FILE__,__LINE__)
       do icell = 1, num_local_cells
          cell_permutation(icell)=icell ! = Id
       end do

    end if

    ! The following paragraph is not needed but I keep it because the loop can be useful
    ! num_ghost_cells_ = 0
    ! do iface=1,polytope%get_num_n_faces()
    !    if(polytope%get_n_face_dim(iface)<this%num_dims) then ! do not include the polytope itself
    !       partial_count=1
    !       count_it = .true.
    !       do idime = 0, this%num_dims - 1 
    !          if(polytope%n_face_dir_is_fixed(iface,idime)==1) then
    !             partial_count=partial_count*num_local_cells_x_dir(idime)
    !          else
    !             if( (polytope%n_face_dir_coordinate(iface,idime)==0.and.num_left_parts_x_dir(idime)==0) .or. &
    !                 (polytope%n_face_dir_coordinate(iface,idime)==1.and.num_right_parts_x_dir(idime)==0) ) count_it = .false.
    !          end if
    !       end do
    !       if(count_it) num_ghost_cells_ = num_ghost_cells_ + partial_count
    !    end if
    ! end do

    ! N_FACES
    ! =======
    ! Global and local number of n_faces (per direction and total)
    num_nface_types = 0
    do iface=1,polytope%get_num_n_faces()
       if(polytope%get_n_face_dim(iface)<num_dims.and. &
            & polytope%n_face_coordinate(iface)==0) num_nface_types = num_nface_types + 1
    end do

    num_edges = 0
    num_faces = 0
    call memalloc( num_nface_types+1, num_global_n_faces, __FILE__,__LINE__,lb1=0)
    call memalloc( num_nface_types+1, num_total_n_faces, __FILE__,__LINE__,lb1=0)
    call memalloc( num_dims, num_nface_types, num_global_nfaces_x_dir, __FILE__,__LINE__,lb1=0,lb2=0)
    call memalloc( num_dims, num_nface_types, num_total_nfaces_x_dir, __FILE__,__LINE__,lb1=0,lb2=0)
    itype = -1
    do iface=1,polytope%get_num_n_faces()
       if(polytope%get_n_face_dim(iface)<num_dims.and. &
            & polytope%n_face_coordinate(iface)==0) then
          itype = itype + 1 
          !itype = polytope%n_face_type(iface)
          do idime = 0, num_dims - 1
             num_global_nfaces_x_dir(idime,itype) = &
                  & num_cells_x_dim(idime)  + &
                  & 1 - max(polytope%n_face_dir_is_fixed(iface,idime),is_dir_periodic(idime))
             num_total_nfaces_x_dir(idime,itype) =  &
                  & num_total_cells_x_dir(idime) + &
                  & 1 - max(polytope%n_face_dir_is_fixed(iface,idime),is_dir_periodic(idime)/num_parts_x_dim_x_level(idime))
          end do
          num_global_n_faces(itype+1) = 1
          num_total_n_faces(itype+1) = 1
          do idime = 0, num_dims - 1
             num_global_n_faces(itype+1) = num_global_n_faces(itype+1) * num_global_nfaces_x_dir(idime,itype)
             num_total_n_faces(itype+1) = num_total_n_faces(itype+1) * num_total_nfaces_x_dir(idime,itype)
          end do
          if(polytope%get_n_face_dim(iface)==num_dims-1) then
             num_faces = num_faces + num_total_n_faces(itype+1)
          else if(polytope%get_n_face_dim(iface)>0) then
             num_edges = num_edges + num_total_n_faces(itype+1)
          end if
       end if
    end do
    num_global_n_faces(0) = 1
    num_total_n_faces(0) = 1
    itype = -1
    do iface=1,polytope%get_num_n_faces()
       if(polytope%get_n_face_dim(iface)<num_dims.and. &
            & polytope%n_face_coordinate(iface)==0) then
          !itype = polytope%n_face_type(iface)
          itype = itype + 1 
          num_global_n_faces(itype+1) = num_global_n_faces(itype+1) + num_global_n_faces(itype)
          num_total_n_faces(itype+1) = num_total_n_faces(itype+1) + num_total_n_faces(itype)
       end if
    end do
    num_local_vefs = num_total_n_faces(num_nface_types) - 1
    num_vertices = num_total_n_faces(1) - 1

    ! FILL ARRAYS
    ! Construct local numbering (ptr_vefs_x_cell does not require permutation because all cells have the same number
    ! and, further, the accumulation does not work if it is applied)
    call memalloc(num_local_cells+num_ghost_cells_+1, ptr_vefs_x_cell, __FILE__,__LINE__)
    ptr_vefs_x_cell = polytope%get_num_n_faces() - 1 ! the cell itself does not count
    ptr_vefs_x_cell(1) = 1
    do icell = 1, num_local_cells+num_ghost_cells_
       ptr_vefs_x_cell(icell+1) = ptr_vefs_x_cell(icell+1) + ptr_vefs_x_cell(icell)
    end do
    call memalloc( ptr_vefs_x_cell(num_local_cells+num_ghost_cells_+1)-1, lst_vefs_lids, __FILE__,__LINE__)

    do icell = 1, num_local_cells+num_ghost_cells_
       call spatial_to_ijk_numbering(num_dims, num_total_cells_x_dir, icell, cell_ijk)
       itype = -1
       do iface=1,polytope%get_num_n_faces()
          if(polytope%get_n_face_dim(iface)<num_dims) then ! do not include the polytope itself
             if(polytope%n_face_coordinate(iface)==0) itype = itype + 1 
             do idime = 0, num_dims - 1
                nface_ijk(idime) = mod(cell_ijk(idime) + polytope%n_face_dir_coordinate(iface,idime),num_total_nfaces_x_dir(idime,itype))
             end do
             !itype = polytope%n_face_type(iface)
             lst_vefs_lids(ptr_vefs_x_cell(cell_permutation(icell))+iface-1) = num_total_n_faces(itype) + &
                  &  ijk_to_spatial_numbering(num_dims, num_total_nfaces_x_dir(:,itype), nface_ijk)
          end if
       end do
    end do

    if(present(num_ghost_cells)) then
       ! Cells global numbering and part
       call memalloc(num_local_cells+num_ghost_cells_,cells_gids,__FILE__,__LINE__)
       call memalloc(num_local_cells+num_ghost_cells_,cells_mypart,__FILE__,__LINE__)
       do icell = 1, num_local_cells+num_ghost_cells_
          call spatial_to_ijk_numbering(num_dims, num_total_cells_x_dir, icell, cell_ijk)
          do idime = 0, num_dims - 1
             if( (num_left_parts_x_dir(idime)==1.and.cell_ijk(idime)==0)) then
                mypart_ijk(idime)=part_ijk(idime)-1
             else if( (num_right_parts_x_dir(idime)==1.and.cell_ijk(idime)==num_total_cells_x_dir(idime)-1)) then
                mypart_ijk(idime)=part_ijk(idime)+1
             else
                mypart_ijk(idime)=part_ijk(idime)
             end if
          end do
          cells_mypart(cell_permutation(icell)) = 1 + &
               &   ijk_to_spatial_numbering( num_dims, &
               &                             num_parts_x_dim_x_level, mypart_ijk)
          cell_ijk = first_cell_ijk + cell_ijk
          cells_gids(cell_permutation(icell)) = 1 + ijk_to_spatial_numbering(num_dims, num_cells_x_dim, cell_ijk)
       end do

       ! List ghost cells and compute interface cells pointers
       call memalloc(num_itfc_cells,lst_itfc_cells, __FILE__,__LINE__)
       call memalloc(num_itfc_cells+1,ptr_ext_neighs_x_itfc_cell,__FILE__,__LINE__)
       itfc_cells = 0
       do icell = 1, num_local_cells+num_ghost_cells_
          if(cell_permutation(icell)>(num_local_cells-num_itfc_cells).and.cell_permutation(icell)<=num_local_cells) then ! cell is interface
             call spatial_to_ijk_numbering(num_dims, num_total_cells_x_dir, icell, cell_ijk)
             index = 0
             do iface=1,polytope%get_num_n_faces()
                if(polytope%get_n_face_dim(iface)<num_dims) then
                   count_it = .false.
                   do idime = 0, num_dims - 1
                      neighbor_ijk(idime) = cell_ijk(idime) - 1 + &
                           & 2 * polytope%n_face_dir_coordinate(iface,idime) + &
                           & polytope%n_face_dir_is_fixed(iface,idime)
                      if(neighbor_ijk(idime)<0.or.neighbor_ijk(idime)>num_total_cells_x_dir(idime)-1) then
                         count_it = .false.  ! the neighbor is out of the domain.
                         exit
                      else if( (num_left_parts_x_dir(idime)==1.and.neighbor_ijk(idime)==0) .or. &
                           &   (num_right_parts_x_dir(idime)==1.and.neighbor_ijk(idime)==num_total_cells_x_dir(idime)-1)) then
                         count_it = .true.
                      end if
                   end do
                   if(count_it) index = index + 1
                end if
             end do
             if(index>0) then
                itfc_cells = itfc_cells + 1
                lst_itfc_cells(itfc_cells) = cell_permutation(icell)
                ptr_ext_neighs_x_itfc_cell(itfc_cells+1)=index
             end if
          end if
       end do
       assert(itfc_cells==num_itfc_cells)

       ! Point to head
       ptr_ext_neighs_x_itfc_cell(1)=1
       do itfc_cells=1,num_itfc_cells
          ptr_ext_neighs_x_itfc_cell(itfc_cells+1) = ptr_ext_neighs_x_itfc_cell(itfc_cells+1) + ptr_ext_neighs_x_itfc_cell(itfc_cells)
       end do

       ! List interface cells neighbors and neighbor parts
       call memalloc(ptr_ext_neighs_x_itfc_cell(num_itfc_cells+1)-1, lst_ext_neighs_gids ,__FILE__,__LINE__)
       call memalloc(ptr_ext_neighs_x_itfc_cell(num_itfc_cells+1)-1, lst_ext_neighs_part_ids ,__FILE__,__LINE__)
       itfc_cells = 1
       do icell = 1, num_local_cells+num_ghost_cells_
          if(cell_permutation(icell)>(num_local_cells-num_itfc_cells).and.cell_permutation(icell)<=num_local_cells) then ! cell is interface
             call spatial_to_ijk_numbering(num_dims, num_total_cells_x_dir, icell, cell_ijk)
             index = 0
             do iface=1,polytope%get_num_n_faces()
                if(polytope%get_n_face_dim(iface)<num_dims) then
                   count_it = .false.
                   do idime = 0, num_dims - 1
                      neighbor_ijk(idime) = cell_ijk(idime) - 1 + &
                           & 2 * polytope%n_face_dir_coordinate(iface,idime) + &
                           & polytope%n_face_dir_is_fixed(iface,idime)
                      if(neighbor_ijk(idime)<0.or.neighbor_ijk(idime)>num_total_cells_x_dir(idime)-1) then
                         count_it = .false.  ! the neighbor is out of the domain.
                         exit
                      else if( (num_left_parts_x_dir(idime)==1.and.neighbor_ijk(idime)==0)) then
                         neighbor_part_ijk(idime)=part_ijk(idime)-1
                         count_it = .true.
                      else if( (num_right_parts_x_dir(idime)==1.and.neighbor_ijk(idime)==num_total_cells_x_dir(idime)-1)) then
                         neighbor_part_ijk(idime)=part_ijk(idime)+1
                         count_it = .true.
                      else
                         neighbor_part_ijk(idime)=part_ijk(idime)
                      end if
                   end do
                   if(count_it) then
                      lst_ext_neighs_gids(ptr_ext_neighs_x_itfc_cell(itfc_cells)+index)= cells_gids(cell_permutation(1 + &
                           &   ijk_to_spatial_numbering( num_dims, num_total_cells_x_dir, neighbor_ijk)))

                      ! This should work too (test it!):
                      ! neighbor_ijk = first_cell_ijk + neighbor_ijk
                      ! lst_ext_neighs_gids(ptr_ext_neighs_x_itfc_cell(itfc_cells)+index)= 1 + &
                      !      &   ijk_to_spatial_numbering( num_dims, &
                      !      &                             num_cells_x_dir, neighbor_ijk)

                      lst_ext_neighs_part_ids(ptr_ext_neighs_x_itfc_cell(itfc_cells)+index)= 1 + &
                           &   ijk_to_spatial_numbering( num_dims, &
                           &                             num_parts_x_dim_x_level, neighbor_part_ijk)
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
    ! In case of 2D domains (num_dims=2) when SPACE_DIM is used, it is necessary to initialize the
    ! coordinates array to zero in order to guarantee that the third component is initialized to zero.
    ! The use of SPACE_DIM instead num_dims is based on the fact that this variable is known in
    ! compilation time, allowing the compiler to perform additional optimizations.
    coordinates = 0.0_rp
    call memalloc(num_local_vefs,boundary_id,__FILE__,__LINE__)
    boundary_id=-1
    itype = -1
    do iface=1,polytope%get_num_n_faces()
       if(polytope%get_n_face_dim(iface)<num_dims.and. &
            & polytope%n_face_coordinate(iface)==0) then
          itype = itype + 1 
          !itype = polytope%n_face_type(iface)
          do iface_of_itype = num_total_n_faces(itype), num_total_n_faces(itype+1) - 1
             call spatial_to_ijk_numbering(num_dims, num_total_nfaces_x_dir(:,itype), &
                  &                        iface_of_itype + 1 - num_total_n_faces(itype), nface_ijk)
             nface_ijk = first_cell_ijk + nface_ijk
             if(present(num_ghost_cells)) &
             vefs_gids(iface_of_itype) = num_global_n_faces(itype) + &
                  &                      ijk_to_spatial_numbering( num_dims, &
                  &                                                num_global_nfaces_x_dir(:,itype), &
                  &                                                nface_ijk )
             if(itype==0) then
                do idime = 0, num_dims - 1 
                   coordinates(idime+1,iface_of_itype) = real(nface_ijk(idime),rp) / real(num_cells_x_dim(idime),rp)
                end do
             end if
             index = 0
             do idime = 0, num_dims - 1 
                if(is_dir_periodic(idime)==0) then ! Not periodic
                   if(  polytope%n_face_dir_is_fixed(iface,idime)==0.and.nface_ijk(idime)==0) then 
                      ! idime bit is already 0
                   else if(polytope%n_face_dir_is_fixed(iface,idime)==0.and.nface_ijk(idime)==num_global_nfaces_x_dir(idime,itype)-1) then 
                      index = ibset( index, idime )
                   else
                      index = ibset( index, num_dims + idime ) ! Fix this coordinate
                   end if
                else
                   index = ibset( index, num_dims + idime ) ! Fix this coordinate
                end if
             end do
             boundary_id(iface_of_itype) = polytope%get_ijk_to_index(index)
          end do
       end if
    end do
    
    ! Map coordinates from [0,1]x[0,1]x[0,1] to [xi,xe]x[yi,ye]x[zi,ze]
    do idime = 1, num_dims
      coordinates(idime,:) = (domain_limits(idime,2)-domain_limits(idime,1))*coordinates(idime,:) + domain_limits(idime,1)
    end do

    call memfree( num_global_n_faces, __FILE__,__LINE__)
    call memfree( num_total_n_faces, __FILE__,__LINE__)
    call memfree( num_global_nfaces_x_dir, __FILE__,__LINE__)
    call memfree( num_total_nfaces_x_dir, __FILE__,__LINE__)
    call memfree( num_parts_x_dim_x_level, __FILE__, __LINE__ )
    call memfree( num_cells_x_dim, __FILE__, __LINE__)
    call memfree( is_dir_periodic, __FILE__, __LINE__)
    call memfree( domain_limits, __FILE__, __LINE__)


    call memfree(cell_permutation, __FILE__,__LINE__)

    !call node_array%free()
    call polytope%free()

  end subroutine uhm_generate_part_local_mesh

  
  subroutine uhm_check_part_local_mesh_parameters(this)
    implicit none
    class(uniform_hex_mesh_t) , intent(in) :: this
    character(len=:), allocatable :: error_message
#ifdef DEBUG
    error_message = struct_hex_mesh_generator_num_dims_key //  & 
                    ' (' // ch(this%num_dims) // ') ' //  & 
                    ' must be either 2 or 3'
    massert ( this%num_dims == 2 .or. this%num_dims == 3, error_message ) 
    
    error_message = struct_hex_mesh_generator_num_dims_key //  & 
                    ' (' // ch(this%num_dims) // ') ' //  & 
                    ' must be smaller or equal than the SPACE_DIM' // &
                    ' (' // ch(SPACE_DIM) // ') ' //  & 
                    'FEMPAR library-level constant parameter'
    massert ( this%num_dims <= SPACE_DIM, error_message )

    massert ( size(this%num_cells_x_dim, dim=1) == this%num_dims, error_message )
    massert ( size(this%is_periodic, dim=1) == this%num_dims, error_message )
    massert ( size(this%num_parts_x_dim_x_level, dim=1) >= this%num_dims, error_message )
    massert ( size(this%domain_limits, dim=1) == 2*this%num_dims, error_message )
   
#endif
  end subroutine uhm_check_part_local_mesh_parameters
  
  subroutine uhm_generate_part_aggregation_among_levels(this, &
                                                        task_id, & 
                                                        num_levels, &
                                                        num_parts_x_level, &
                                                        part_aggregation_among_levels)

    implicit none
    class(uniform_hex_mesh_t) , intent(in)    :: this
    integer(ip)               , intent(in)    :: task_id
    integer(ip)               , intent(in)    :: num_levels
    integer(ip)               , intent(in)    :: num_parts_x_level(num_levels)
    integer(ip)  , allocatable, intent(inout) :: part_aggregation_among_levels(:)
    
    integer(ip) :: ilevel,idime,ipart,first,last
    integer(ip) :: part_ijk(0:SPACE_DIM-1)

    assert(task_id<sum(num_parts_x_level))

    call this%check_part_aggregation_among_levels_parameters(num_levels, num_parts_x_level)

    if ( allocated(part_aggregation_among_levels) ) & 
        call memfree(part_aggregation_among_levels, __FILE__,__LINE__)
        
    call memalloc(num_levels, part_aggregation_among_levels, __FILE__,__LINE__)
    part_aggregation_among_levels = -1
    ilevel=1
    ipart=task_id + 1 
    do while(ipart>num_parts_x_level(ilevel))
       ipart  = ipart - num_parts_x_level(ilevel)
       ilevel = ilevel + 1
    end do
    part_aggregation_among_levels(ilevel) = ipart
    do while(ilevel<=num_levels-1)
       first = (ilevel-1)*SPACE_DIM
       last  = first + this%num_dims-1
       call spatial_to_ijk_numbering(this%num_dims, this%num_parts_x_dim_x_level(first:last), ipart, part_ijk)
       do idime = 0, this%num_dims - 1 
          part_ijk(idime) = part_ijk(idime)*this%num_parts_x_dim_x_level(ilevel*SPACE_DIM+idime)/this%num_parts_x_dim_x_level((ilevel-1)*SPACE_DIM+idime)
       end do
       first = ilevel*SPACE_DIM
       last  = first + this%num_dims-1
       ipart = ijk_to_spatial_numbering(this%num_dims,this%num_parts_x_dim_x_level(first:last), part_ijk)+1
       ilevel = ilevel +1
       part_aggregation_among_levels(ilevel) = ipart
    end do
  end subroutine uhm_generate_part_aggregation_among_levels

  subroutine uhm_check_part_aggregation_among_levels_parameters(this, num_levels, num_parts_x_level)
    implicit none
    class(uniform_hex_mesh_t) , intent(in) :: this
    integer(ip)               , intent(in) :: num_levels
    integer(ip)               , intent(in) :: num_parts_x_level(0:num_levels-1)
    integer(ip)                            :: dims, i
    character(len=:), allocatable          :: error_message
#ifdef DEBUG

    error_message = 'size('//struct_hex_mesh_generator_num_parts_x_dim_x_level_key//')'//  & 
                    ' (' // ch(size(this%num_parts_x_dim_x_level)) // ') ' //  & 
                    ' must be ' // ch(SPACE_DIM*num_levels)
    massert( size(this%num_parts_x_dim_x_level) == SPACE_DIM*num_levels, error_message )


    do i = 0, num_levels-1
        error_message = struct_hex_mesh_generator_num_parts_x_dim_x_level_key//  & 
                        ' parts at level (' // ch(num_levels-i) // ') ' //  & 
                        ' must be ' // ch(num_parts_x_level(i))
        massert(product(this%num_parts_x_dim_x_level(SPACE_DIM*i:SPACE_DIM*i+this%num_dims-1)) == num_parts_x_level(i), error_message)
    enddo
    
#endif
  end subroutine uhm_check_part_aggregation_among_levels_parameters
  
  
  pure function ijk_to_spatial_numbering_ip(num_dims, num_x_dim, ijk)
    implicit none
    integer(ip)           , intent(in) :: num_dims
    integer(ip)           , intent(in) :: num_x_dim(0:SPACE_DIM-1) 
    integer(ip)           , intent(in) :: ijk(0:SPACE_DIM-1) 
    integer(ip) :: ijk_to_spatial_numbering_ip
    integer(ip) :: idime, jdime
    integer(ip) :: previous
    ijk_to_spatial_numbering_ip = 0
    do idime = 0, num_dims - 1
       previous = 1
       do jdime = 0, idime - 1 
          previous = previous * num_x_dim(jdime)
       end do
       ijk_to_spatial_numbering_ip = ijk_to_spatial_numbering_ip + previous*ijk(idime)
    end do
  end function ijk_to_spatial_numbering_ip

  pure subroutine spatial_to_ijk_numbering(num_dims, num_x_dim, spatial_numbering, ijk)
    implicit none
    integer(ip)           , intent(in)  :: num_dims
    integer(ip)           , intent(in)  :: num_x_dim(0:SPACE_DIM-1) 
    integer(ip)           , intent(in)  :: spatial_numbering
    integer(ip)           , intent(out) :: ijk(0:SPACE_DIM-1) 
    integer(ip) :: idime,j

    j = spatial_numbering - 1          ! To make it 0-based (assuming spatial_numbering is 1-based)
    do idime = 0, num_dims - 1
       ijk(idime) = mod(j,num_x_dim(idime))
       j = j / num_x_dim(idime)
    end do

  end subroutine spatial_to_ijk_numbering
  
end module uniform_hex_mesh_generator_names
