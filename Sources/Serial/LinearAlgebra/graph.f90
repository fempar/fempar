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
module graph_names
  use types_names
  use memor_names
  implicit none
# include "debug.i90"

  private

  ! Graph (CSR format, may exploit symmetry if desired)
  type graph_t
     integer(ip)                :: &
          nv,                      &    ! Number of vertices  (rows)
          nv2                           ! Number of vertices2 (columns)

     logical :: symmetric_storage       ! .True.   Implicitly assumes that G=(V,E) is such that 
                                        !          (i,j) \belongs E <=> (j,i) \belongs E, forall i,j \belongs V.
                                        !          Only edges (i,j) with j>=i are stored.
                                        ! .False.  All (i,j) \belongs E are stored.  

     integer(ip), allocatable   :: &
          ia(:),                     &    ! Indices to adjacencies        
          ja(:)                           ! Adjacencies
   contains
     procedure :: create => graph_create
     procedure :: print  => graph_print
     procedure :: copy   => graph_copy
     procedure :: free   => graph_free
  end type graph_t

  ! Types
  public :: graph_t


contains
  
  !=============================================================================
  subroutine graph_create (this,symmetric_storage)
    implicit none
    class(graph_t), intent(out) :: this
    logical      , intent(in)  :: symmetric_storage
    this%symmetric_storage = symmetric_storage
  end subroutine graph_create

  !=============================================================================
  subroutine graph_copy (igraph, ograph)
    implicit none
    class(graph_t), intent(in)    :: igraph
    type(graph_t), intent(out) :: ograph
    ograph%symmetric_storage = igraph%symmetric_storage
    ograph%nv   = igraph%nv
    ograph%nv2  = igraph%nv2
    ! Alloc/copy ia array
    call memalloc ( igraph%nv+1, ograph%ia, __FILE__,__LINE__ )
    ograph%ia = igraph%ia
    ! Allocate space/copy ja array 
    call memalloc ( igraph%ia(igraph%nv+1)-1, ograph%ja,__FILE__,__LINE__) 
    ograph%ja = igraph%ja
  end subroutine graph_copy
  
  !=============================================================================
  subroutine graph_free(this)
    implicit none
    class(graph_t), intent(inout)  :: this
    call memfree (this%ia,__FILE__,__LINE__)
    call memfree (this%ja,__FILE__,__LINE__)
  end subroutine graph_free

  !=============================================================================
  subroutine graph_print(graph,lunou)
    implicit none
    class(graph_t),  intent(in) :: graph
    integer(ip)    ,  intent(in) :: lunou

    ! Local variables
    integer(ip) :: i,j

    write (lunou, '(a)')     '*** begin graph data structure ***'
    write (lunou, '(a,i10)') 'Number of vertices (rows):', graph%nv
    write (lunou, '(a,i10)') 'Number of vertices (cols):', graph%nv2

    write (lunou, '(a)')     'Pointers to adjacency lists (ia):'
    write (lunou, '(10i10)')    graph%ia(1:(graph%nv+1))

    write (lunou, '(a)')      'Adjacency lists of each vertex:'
    do i=1,graph%nv
       write(lunou,'(10i10)') i, graph%ja(graph%ia(i):graph%ia(i+1)-1)
    end do
    write (lunou, '(a)')     '*** end graph data structure ***'
  end subroutine graph_print

end module graph_names
