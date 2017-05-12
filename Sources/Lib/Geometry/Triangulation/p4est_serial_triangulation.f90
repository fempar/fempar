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
module p4est_serial_triangulation_names
  use, intrinsic :: iso_c_binding
  
  use types_names
  use stdio_names
  use memor_names
  use p4est_bindings_names
  
  use FPL

  implicit none
# include "debug.i90"
  private

  integer(ip), parameter :: NUM_CORNERS_2D           = 4
  integer(ip), parameter :: NUM_FACES_2D             = 4
  integer(ip), parameter :: NUM_SUBFACES_FACE_2D     = 2
  integer(ip), parameter :: NUM_FACE_CORNERS_2D      = 2
  integer(ip), parameter :: NUM_VEFS_2D              = NUM_CORNERS_2D+NUM_FACES_2D
  integer(ip), parameter :: P4EST_FACE_CORNERS_2D(0:NUM_FACE_CORNERS_2D-1,0:NUM_FACES_2D-1) = & 
                                                  reshape([0, 2,&
                                                           1, 3,&  
                                                           0, 1,&
                                                           2, 3], [NUM_FACE_CORNERS_2D,NUM_FACES_2D])
                                                  
                                                  
  
  
  
  ! TODO: this data type should extend an abstract triangulation,
  !       and implement its corresponding accessors
  type p4est_serial_triangulation_t
    private
    integer(ip) :: num_cells          = -1
    integer(ip) :: num_dimensions     = -1
    integer(ip) :: num_vefs           = -1
    integer(ip) :: num_proper_vefs    = -1 
    integer(ip) :: num_improper_vefs  = -1 
    type(c_ptr) :: p4est_connectivity = c_null_ptr
    type(c_ptr) :: p4est              = c_null_ptr
    type(c_ptr) :: p4est_mesh         = c_null_ptr
    ! TODO: I am pretty sure that a type(c_ptr) :: p4est_ghost
    !       member variable will be needed (at least in the parallel realization)
    
    ! p4est quadrant connectivity (1:NUM_FACES_2D/3D,1:nQuads) => neighbor quadrant
    integer(P4EST_F90_LOCIDX),pointer     :: QuadToQuad(:,:) => NULL()
    ! p4est face connectivity (1:NUM_FACES_2D/3D,1:nQuads) => neighbor faceId + orientation + non-conform info
    integer(P4EST_F90_QLEVEL),pointer     :: QuadToFace(:,:) => NULL()   
    ! p4est face connectivity for mortars NUM_SUBFACES_FACE_2D/3D,1:nHalfFaces), (~small sides)
    integer(P4EST_F90_LOCIDX),pointer     :: QuadToHalf(:,:) => NULL()
    
    ! TODO: The following 4x member variables should be replaced by our F200X implementation of "std::vector<T>" 
    ! p4est Integer coordinates of first quadrant node (xy/xyz,nQuads)
    integer(P4EST_F90_LOCIDX), allocatable :: QuadCoords(:,:)
    ! p4est Integer Level of quadrant
    integer(P4EST_F90_QLEVEL), allocatable :: QuadLevel(:)
    integer(ip)              , allocatable :: ptr_vefs_per_cell(:)
    integer(ip)              , allocatable :: lst_vefs_lids(:)    
  contains
    procedure, non_overridable          :: create                          => p4est_serial_triangulation_create
    procedure, non_overridable          :: free                            => p4est_serial_triangulation_free
    procedure, non_overridable          :: refine_and_coarsen              => p4est_serial_triangulation_refine_and_coarsen
    procedure, private, non_overridable :: update_p4est_mesh               => p4est_serial_triangulation_update_p4est_mesh
    procedure, private, non_overridable :: update_topology_from_p4est_mesh => p4est_serial_triangulation_update_topology_from_p4est_mesh
    
    procedure, private, non_overridable :: update_ptr_vefs_per_cell        => p4est_serial_triangulation_update_ptr_vefs_per_cell
    procedure, private, non_overridable :: update_lst_vefs_lids            => p4est_serial_triangulation_update_ptr_vefs_per_cell
    procedure, private, non_overridable :: fill_lst_vefs_lids              => p4est_serial_triangulation_fill_lst_vefs_lids
    procedure, private, non_overridable :: free_ptr_vefs_per_cell          => p4est_serial_triangulation_free_ptr_vefs_per_cell
    procedure, private, non_overridable :: free_lst_vefs_lids              => p4est_serial_triangulation_free_lst_vefs_lids
#ifndef ENABLE_P4EST
    procedure, non_overridable :: not_enabled_error => p4est_serial_triangulation_not_enabled_error
#endif
  end type p4est_serial_triangulation_t
  
  public :: p4est_serial_triangulation_t
  
contains

subroutine p4est_serial_triangulation_create (this, parameters)
  implicit none
  class(p4est_serial_triangulation_t), intent(inout) :: this
  type(ParameterList_t)              , intent(in)    :: parameters
  
#ifdef ENABLE_P4EST
  call this%free()
  this%num_cells = 1
  
  ! TODO: Extract num_dimensions out of parameters
  this%num_dimensions = 2
  
  if ( this%num_dimensions == 2 ) then
    call F90_p4est_connectivity_new_unitsquare(this%p4est_connectivity)
    call F90_p4est_new(this%p4est_connectivity, this%p4est)
    call this%update_p4est_mesh()
    call this%update_topology_from_p4est_mesh()
    call this%update_ptr_vefs_per_cell()
    call this%update_lst_vefs_lids()
  else if ( this%num_dimensions == 3 ) then
    check(.false.)
  end if  
#else
  call this%not_enabled_error()
#endif
end subroutine p4est_serial_triangulation_create  

subroutine p4est_serial_triangulation_refine_and_coarsen(this)
  implicit none
  class(p4est_serial_triangulation_t), intent(inout) :: this
  
#ifdef ENABLE_P4EST
  if ( this%num_dimensions == 2 ) then
    call F90_p4est_refine(this%p4est)
    call this%update_p4est_mesh()
    call this%update_topology_from_p4est_mesh()
    call this%update_ptr_vefs_per_cell()
    call this%update_lst_vefs_lids()
  else if ( this%num_dimensions == 3 ) then
    check(.false.)
  end if
  
  ! Update the number of triangulation cells
  this%num_cells = size(this%QuadLevel)
  
#else
  call this%not_enabled_error()
#endif  
  
end subroutine p4est_serial_triangulation_refine_and_coarsen

subroutine p4est_serial_triangulation_update_ptr_vefs_per_cell(this)
  implicit none
  class(p4est_serial_triangulation_t), intent(inout) :: this
  integer(ip) :: icell
  integer(ip) :: num_vefs_per_cell

#ifdef ENABLE_P4EST
  call this%free_ptr_vefs_per_cell()
  
  if ( this%number_dimensions == 2 ) then
    num_vefs_per_cell = NUM_VEFS_2D
  else if ( this%number_dimensions == 3 ) then
    num_vefs_per_cell = 
  end if
  
  
  call memalloc(this%num_cells+1, this%ptr_vefs_per_cell, __FILE__, __LINE__)
  this%ptr_vefs_per_cell(1)=1
  do icell=1, this%num_cells
    if ( this%number_dimensions == 2 ) then
       this%ptr_vefs_per_cell(icell+1) = this%ptr_vefs_per_cell(icell) + 
    end if
  end do
#else
  call this%not_enabled_error()
#endif  
  
end subroutine p4est_serial_triangulation_update_ptr_vefs_per_cell

subroutine p4est_serial_triangulation_update_lst_vefs_per_cell(this)
  implicit none
  class(p4est_serial_triangulation_t), intent(inout) :: this
end subroutine p4est_serial_triangulation_update_lst_vefs_per_cell

subroutine p4est_serial_triangulation_fill_lst_vefs_lids(this)
  implicit none
  class(p4est_serial_triangulation_t), intent(inout) :: this
end subroutine p4est_serial_triangulation_fill_lst_vefs_lids

subroutine p4est_serial_triangulation_free_ptr_vefs_per_cell(this)
  implicit none
  class(p4est_serial_triangulation_t), intent(inout) :: this  
  if (allocated(this%ptr_vefs_per_cell)) &
    call memfree(this%ptr_vefs_per_cell, __FILE__, __LINE__)
end subroutine p4est_serial_triangulation_free_ptr_vefs_per_cell

subroutine p4est_serial_triangulation_free_lst_vefs_lids(this)
  implicit none
  class(p4est_serial_triangulation_t), intent(inout) :: this
  if (allocated(this%lst_vefs_lids)) &
    call memfree(this%lst_vefs_lids, __FILE__, __LINE__)
end subroutine p4est_serial_triangulation_free_lst_vefs_lids

subroutine p4est_serial_triangulation_update_p4est_mesh(this)
  implicit none
  class(p4est_serial_triangulation_t), intent(inout) :: this
  
#ifdef ENABLE_P4EST
  if ( this%num_dimensions == 2 ) then
    call F90_p4est_mesh_new(this%p4est, this%p4est_mesh)
  else if ( this%num_dimensions == 3 ) then
    check(.false.)
  end if
#else
  call this%not_enabled_error()
#endif   
end subroutine p4est_serial_triangulation_update_p4est_mesh

subroutine p4est_serial_triangulation_update_topology_from_p4est_mesh(this)
 implicit none 
 class(p4est_serial_triangulation_t), intent(inout) :: this
 integer(P4EST_F90_LOCIDX) :: local_num_quadrants
 integer(P4EST_F90_GLOIDX) :: global_num_quadrants
 integer(P4EST_F90_GLOIDX) :: global_first_quadrant
 integer(P4EST_F90_LOCIDX) :: num_half_faces
 type(c_ptr) :: QQ, QF, QH
 
#ifdef ENABLE_P4EST
 if ( this%num_dimensions == 2 ) then
  call F90_p4est_get_mesh_info(this%p4est, &
                               this%p4est_mesh, &
                               local_num_quadrants, &
                               global_num_quadrants, &
                               global_first_quadrant, &
                               num_half_faces)
 
  if (allocated(this%QuadCoords)) &
     call memfree(this%QuadCoords, __FILE__, __LINE__)
  
  if (allocated(this%QuadLevel)) &
    call memfree(this%QuadLevel, __FILE__, __LINE__)
  
  call memalloc(2, local_num_quadrants, this%QuadCoords, __FILE__, __LINE__)
  call memalloc(local_num_quadrants, this%QuadLevel, __FILE__, __LINE__ )
  
  call F90_p4est_get_mesh_topology_arrays(this%p4est, &
                                          this%p4est_mesh, &
                                          QQ, &
                                          QF, &
                                          QH, &
                                          this%QuadCoords, &
                                          this%QuadLevel)
  
  call c_f_pointer(qq,this%quadtoquad,[NUM_FACES_2D,local_num_quadrants])
  call c_f_pointer(qf,this%quadtoface,[NUM_FACES_2D,local_num_quadrants])
  if(num_half_faces>0) call c_f_pointer(qh,this%quadtohalf,[NUM_SUBFACES_FACE_2D,num_half_faces])
 else if ( this%num_dimensions == 3 ) then
   check(.false.)
 end if 
  
#else
  call this%not_enabled_error()
#endif
end subroutine p4est_serial_triangulation_update_topology_from_p4est_mesh

subroutine p4est_serial_triangulation_free ( this)
  implicit none
  class(p4est_serial_triangulation_t), intent(inout) :: this

#ifdef ENABLE_P4EST
  if ( this%num_dimensions == 2 ) then
    call F90_p4est_destroy(this%p4est)
    call F90_p4est_connectivity_destroy(this%p4est_connectivity)
    call F90_p4est_mesh_destroy(this%p4est_mesh)
  
    this%p4est_connectivity = c_null_ptr
    this%p4est              = c_null_ptr
    this%p4est_mesh         = c_null_ptr
  else if ( this%num_dimensions == 3 ) then
    check(.false.)
  end if
  
  nullify(this%QuadToQuad)
  nullify(this%QuadToFace)
  nullify(this%QuadToHalf)
  
  if (allocated(this%QuadCoords)) &
     call memfree(this%QuadCoords, __FILE__, __LINE__)
  
  if (allocated(this%QuadLevel)) &
    call memfree(this%QuadLevel, __FILE__, __LINE__)
  
  this%num_dimensions  = -1
  this%num_cells = -1
  this%num_vefs = -1
  this%num_proper_vefs = -1
  this%num_improper_vefs = -1
#else
  call this%not_enabled_error()
#endif     
end subroutine p4est_serial_triangulation_free


#ifndef ENABLE_P4EST
  subroutine p4est_serial_triangulation_not_enabled_error(this)
    class(p4est_serial_triangulation_t), intent(inout) :: this
    write (stderr,*) 'Error: FEMPAR was not compiled with -DENABLE_P4EST.'
    write (stderr,*) "Error: You must activate this CPP macro in order to use P4EST"
    check(.false.)
  end subroutine p4est_serial_triangulation_not_enabled_error
#endif

end module p4est_serial_triangulation_names
