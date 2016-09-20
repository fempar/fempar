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
module hex_boundary_set_ids_descriptor_names
  ! Serial modules
  use types_names
  use memor_names
  implicit none
# include "debug.i90"
  private
   
  ! Convention used within this uniform hex mesh generator module 
  ! for numbering objects within an hexaedron (lexicographic, ZYX)
  ! 2D:
  ! Points code: 1=bottom left; 2=top left;   3=bottom right; 4=top right
  ! Lines code:  1=left line;   2=right line; 3=bottom line;  4=top line 

  ! 3D:
  ! Points code:  1=(0,0,0), 2=(0,0,1), 3=(0,1,0), 4=(0,1,1), 5=(1,0,0),
  !               6=(1,0,1), 7=(1,1,0), 8=(1,1,1)
  ! Lines code:   1=(0,0,z), 2=(0,1,z), 3=(1,0,z), 4=(1,1,z), 5=(0,y,0),
  !               6=(0,y,1), 7=(1,y,0), 8=(1,y,1), 9=(x,0,0), 10=(x,0,1),
  !               11=(x,1,0), 12=(x,1,1)
  ! Surface code: 1=(x,y,0), 2=(x,y,1), 3=(x,0,z), 4=(x,1,z), 5=(0,y,z),
  !               6=(1,y,z)
  ! Formely uniform_conditions_descriptor_t
  type hex_boundary_set_ids_descriptor_t 
    private
    integer(ip)              :: number_dimensions = -1
    integer(ip), allocatable :: boundary_set_ids_vertices(:)
    integer(ip), allocatable :: boundary_set_ids_edges(:)
    integer(ip), allocatable :: boundary_set_ids_faces(:)
  contains
    procedure, non_overridable          :: create                       => hex_boundary_set_ids_descriptor_create
    procedure, non_overridable          :: free                         => hex_boundary_set_ids_descriptor_free
    procedure, non_overridable          :: set_vertex_boundary_set_id   => hex_boundary_set_ids_descriptor_set_vertex_boundary_set_id
    procedure, non_overridable          :: set_edge_boundary_set_id     => hex_boundary_set_ids_descriptor_set_edge_boundary_set_id
    procedure, non_overridable          :: set_face_boundary_set_id     => hex_boundary_set_ids_descriptor_set_face_boundary_set_id
    procedure, non_overridable          :: get_vertex_boundary_set_id   => hex_boundary_set_ids_descriptor_get_vertex_boundary_set_id
    procedure, non_overridable          :: get_edge_boundary_set_id     => hex_boundary_set_ids_descriptor_get_edge_boundary_set_id
    procedure, non_overridable          :: get_face_boundary_set_id     => hex_boundary_set_ids_descriptor_get_face_boundary_set_id
    procedure, non_overridable, private :: get_number_vertices          => hex_boundary_set_ids_descriptor_get_number_vertices
    procedure, non_overridable, private :: get_number_edges             => hex_boundary_set_ids_descriptor_get_number_edges
    procedure, non_overridable, private :: get_number_faces             => hex_boundary_set_ids_descriptor_get_number_faces
  end type hex_boundary_set_ids_descriptor_t
  
  public :: hex_boundary_set_ids_descriptor_t
  
contains

subroutine hex_boundary_set_ids_descriptor_create(this, number_dimensions)
  implicit none
  class(hex_boundary_set_ids_descriptor_t), intent(inout) :: this
  integer(ip)                             , intent(in)    :: number_dimensions
  assert ( number_dimensions == 2 .or. number_dimensions == 3 )
  call this%free()
  this%number_dimensions = number_dimensions
  call memalloc (this%get_number_vertices(), this%boundary_set_ids_vertices, __FILE__, __LINE__ )
  this%boundary_set_ids_vertices = 0
  call memalloc (this%get_number_edges(), this%boundary_set_ids_edges, __FILE__, __LINE__ )
  this%boundary_set_ids_edges = 0
  call memalloc (this%get_number_faces(), this%boundary_set_ids_faces, __FILE__, __LINE__ )
  this%boundary_set_ids_faces = 0
end subroutine hex_boundary_set_ids_descriptor_create

subroutine hex_boundary_set_ids_descriptor_free(this)
  implicit none
  class(hex_boundary_set_ids_descriptor_t), intent(inout) :: this
  this%number_dimensions = -1
  if ( allocated(this%boundary_set_ids_vertices) )  call memfree (this%boundary_set_ids_vertices, __FILE__, __LINE__ )
  if ( allocated(this%boundary_set_ids_edges) )  call memfree (this%boundary_set_ids_edges, __FILE__, __LINE__ )
  if ( allocated(this%boundary_set_ids_faces) )  call memfree (this%boundary_set_ids_faces, __FILE__, __LINE__ )
end subroutine hex_boundary_set_ids_descriptor_free

subroutine hex_boundary_set_ids_descriptor_set_vertex_boundary_set_id ( this, vertex_id, boundary_set_id )
  implicit none
  class(hex_boundary_set_ids_descriptor_t), intent(inout) :: this
  integer(ip)                             , intent(in)    :: vertex_id
  integer(ip)                             , intent(in)    :: boundary_set_id
  assert ( vertex_id >=1 .and. vertex_id <= this%get_number_vertices() )
  this%boundary_set_ids_vertices(vertex_id) = boundary_set_id
end subroutine hex_boundary_set_ids_descriptor_set_vertex_boundary_set_id

subroutine hex_boundary_set_ids_descriptor_set_edge_boundary_set_id ( this, edge_id, boundary_set_id )
  implicit none
  class(hex_boundary_set_ids_descriptor_t), intent(inout) :: this
  integer(ip)                             , intent(in)    :: edge_id
  integer(ip)                             , intent(in)    :: boundary_set_id
  assert ( edge_id >=1 .and. edge_id <= this%get_number_edges() )
  this%boundary_set_ids_edges(edge_id) = boundary_set_id
end subroutine hex_boundary_set_ids_descriptor_set_edge_boundary_set_id

subroutine hex_boundary_set_ids_descriptor_set_face_boundary_set_id ( this, face_id, boundary_set_id )
  implicit none
  class(hex_boundary_set_ids_descriptor_t), intent(inout) :: this
  integer(ip)                             , intent(in)    :: face_id
  integer(ip)                             , intent(in)    :: boundary_set_id
  assert ( face_id >=1 .and. face_id <= this%get_number_faces() )
  this%boundary_set_ids_faces(face_id) = boundary_set_id
end subroutine hex_boundary_set_ids_descriptor_set_face_boundary_set_id

function hex_boundary_set_ids_descriptor_get_vertex_boundary_set_id ( this, vertex_id )
  implicit none
  class(hex_boundary_set_ids_descriptor_t), intent(inout) :: this
  integer(ip)                             , intent(in)    :: vertex_id
  integer(ip) :: hex_boundary_set_ids_descriptor_get_vertex_boundary_set_id
  assert ( vertex_id >=1 .and. vertex_id <= this%get_number_vertices() )
  hex_boundary_set_ids_descriptor_get_vertex_boundary_set_id = this%boundary_set_ids_vertices(vertex_id)
end function hex_boundary_set_ids_descriptor_get_vertex_boundary_set_id

function hex_boundary_set_ids_descriptor_get_edge_boundary_set_id ( this, edge_id )
  implicit none
  class(hex_boundary_set_ids_descriptor_t), intent(inout) :: this
  integer(ip)                             , intent(in)    :: edge_id
  integer(ip) :: hex_boundary_set_ids_descriptor_get_edge_boundary_set_id
  assert ( edge_id >=1 .and. edge_id <= this%get_number_vertices() )
  hex_boundary_set_ids_descriptor_get_edge_boundary_set_id = this%boundary_set_ids_edges(edge_id)
end function hex_boundary_set_ids_descriptor_get_edge_boundary_set_id

function hex_boundary_set_ids_descriptor_get_face_boundary_set_id ( this, face_id )
  implicit none
  class(hex_boundary_set_ids_descriptor_t), intent(inout) :: this
  integer(ip)                             , intent(in)    :: face_id
  integer(ip) :: hex_boundary_set_ids_descriptor_get_face_boundary_set_id
  assert ( face_id >=1 .and. face_id <= this%get_number_vertices() )
  hex_boundary_set_ids_descriptor_get_face_boundary_set_id = this%boundary_set_ids_faces(face_id)
end function hex_boundary_set_ids_descriptor_get_face_boundary_set_id

function hex_boundary_set_ids_descriptor_get_number_vertices(this)
  implicit none
  class(hex_boundary_set_ids_descriptor_t), intent(in) :: this
  integer(ip) :: hex_boundary_set_ids_descriptor_get_number_vertices
  hex_boundary_set_ids_descriptor_get_number_vertices = ISHFT(1_ip, this%number_dimensions) ! 2**number_dimensions
end function hex_boundary_set_ids_descriptor_get_number_vertices

function hex_boundary_set_ids_descriptor_get_number_edges(this)
  implicit none
  class(hex_boundary_set_ids_descriptor_t), intent(in) :: this
  integer(ip) :: hex_boundary_set_ids_descriptor_get_number_edges
  hex_boundary_set_ids_descriptor_get_number_edges = 2*(this%number_dimensions-1)*this%number_dimensions
end function hex_boundary_set_ids_descriptor_get_number_edges

function hex_boundary_set_ids_descriptor_get_number_faces(this)
  implicit none
  class(hex_boundary_set_ids_descriptor_t), intent(in) :: this
  integer(ip) :: hex_boundary_set_ids_descriptor_get_number_faces
  
  ! Following the convention of the uniform_hex_mesh_generator code,
  ! in 2D there are no faces, only edges. This may be re-considered in
  ! the future, such that it becomes consistent with FEMPAR triangulations,
  ! where in 2D there are actually no edges
  hex_boundary_set_ids_descriptor_get_number_faces = 0
  if ( this%number_dimensions == 3 ) then
    hex_boundary_set_ids_descriptor_get_number_faces = 2 * this%number_dimensions
  end if
end function hex_boundary_set_ids_descriptor_get_number_faces

end module hex_boundary_set_ids_descriptor_names

