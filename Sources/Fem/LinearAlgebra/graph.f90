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
module fem_graph_class
  use types
  use memor
  implicit none
# include "debug.i90"

  private

  ! Graph
  type fem_graph
     integer(ip)                :: &
        type=0,                    &         ! Type of fem_graph (csr_symm, csr, css, part)
        nv=0,                      &         ! Number of vertices
        nv2=0,                     &         ! Number of vertices2 (bipartite fem_graphs)
        
        ! The following fields are only 
        ! used for fem_graphs of type css
        nzs=0,                     &         ! Half size
        nzt=0                                ! Total size
     integer(ip), allocatable   :: &
        ia(:),                     &         ! Indices of adjacencies
        
        ! The following field is only 
        ! used for fem_graphs of type css
        is(:),                     &         ! Indices to access lower/upper adjacencies
        
        ja(:)                                ! Adjacencies

  end type fem_graph

  ! See comments on mesh_graph.f90 ...
  integer(ip), parameter :: csr_symm = 10    ! Compress storage row for SQUARE symmetric
                                             ! matrices, diagonal entries stored, column 
                                             ! indices in ascending order (as required, 
                                             ! e.g., by PARDISO)

  integer(ip), parameter :: csr  = 20        ! Compress storage row for square (nv==nv2)
                                             ! or rectangular (nv!=nv2) matrices, diagonal
                                             ! entries stored, column indices in ascending 
                                             ! order (as required, e.g., by PARDISO or 
                                             ! Trilinos, respectively)

  integer(ip), parameter :: css   = 30       ! ia and ja as in csr, wo diagonal entries 
                                             ! stored. "is" is equal to ia in csr_symm but 
                                             ! wo diagonal entries stored. column indices 
                                             ! not necessarily in ascending order

  integer(ip), parameter :: part  = 40       ! Compress storage row for graphs, wo 
                                             ! diagonal entries stored, column indices not
                                             ! necessarily in ascending order (e.g., to 
                                             ! perform graph partitioning with METIS)

  ! csr, csr_symm and css fem_graph types are intended to address sparse
  ! matrices, while part type is used for sequential graph 
  ! partitioning. The separation of the mesh_to_graph interface is
  ! done accordingly to this criteria       

  ! Types
  public :: fem_graph

  ! Functions
  public :: fem_graph_copy, fem_graph_free, fem_graph_print, fem_graph_read, &
       &    fem_graph_ja_one_to_zero_indexing, fem_graph_ja_zero_to_one_indexing
  ! Constants
  public :: csr_symm, csr, part, css

contains
  !=============================================================================
  !
  !
  !=============================================================================
  subroutine fem_graph_copy (igraph, ograph)
    implicit none
    type(fem_graph), intent(in)    :: igraph
    type(fem_graph), intent(inout) :: ograph

    ! This routine is only provided for graphs stored in CSR/CSR_SYMM
    ! format. If you need to copy a graph in a different format,
    ! please extend this routine
    assert ( igraph%type == csr_symm .or. igraph%type == csr )

    ograph%type = igraph%type
    ograph%nv   = igraph%nv
    ograph%nv2  = igraph%nv2

    ! Alloc/copy ia array
    call memalloc ( igraph%nv+1, ograph%ia, __FILE__,__LINE__ )
    ograph%ia = igraph%ia
   
    ! Allocate space/copy ja array 
    call memalloc ( igraph%ia(igraph%nv+1)-1, ograph%ja,__FILE__,__LINE__) 
    ograph%ja = igraph%ja

  end subroutine fem_graph_copy

  !=============================================================================
  !
  !
  !=============================================================================
  subroutine fem_graph_free(g)
    !-----------------------------------------------------------------------
    ! This routine 
    !-----------------------------------------------------------------------
    implicit none
    type(fem_graph), intent(inout)  :: g

    call memfree (g%ia,__FILE__,__LINE__)
    call memfree (g%ja,__FILE__,__LINE__)

    if ( g%type == css ) then
       call memfree(g%is,__FILE__,__LINE__)
    end if

  end subroutine fem_graph_free

  !=============================================================================
  subroutine fem_graph_print(lunou, g)
    implicit none
    integer(ip)    ,  intent(in) :: lunou
    type(fem_graph),  intent(in) :: g

    ! Local variables
    integer(ip) :: i,j

    write (lunou, '(a)')     '*** begin fem_graph data structure ***'
    write (lunou, '(a,i10)') 'Graph type:', g%type    
    write (lunou, '(a,i10)') 'Number of vertices (rows):', g%nv
    write (lunou, '(a,i10)') 'Number of vertices (cols):', g%nv2

    if ( g%type == css ) then
      write (lunou, '(a,i10)') 'Number of edges (i,j) s.t. i != j:',  g%nzt
      write (lunou, '(a,i10)') 'Number of edges (i,j) s.t. i  > j):', g%nzs
    end if

    write (lunou, '(a)')     'Pointers to adjacency lists (ia):'
    write (lunou, '(10i10)')    g%ia(1:(g%nv+1))

    if ( g%type == css ) then
      write (lunou, '(a)')     'Pointers to adjacency lists (is):'
      write (lunou, '(10i10)')    g%is(1:(g%nv+1))
    end if


    write (lunou, '(a)')      'Adjacency lists of each vertex:'
    do i=1,g%nv
       write(lunou,'(10i10)') i, g%ja(g%ia(i):g%ia(i+1)-1)
    end do
    write (lunou, '(a)')     '*** end fem_graph data structure ***'

  end subroutine fem_graph_print

  !=============================================================================
  subroutine fem_graph_read(lunin, g)
    implicit none
    integer(ip)    ,  intent(in)    :: lunin
    type(fem_graph),  intent(inout) :: g
    character(256) :: text
    ! Local variables
    integer(ip) :: i,j

    read (lunin, '(a)')       text
    read (lunin, '(a,i10)')   text(1:11), g%type    
    read (lunin, '(a,i10)')   text(1:26), g%nv
    read (lunin, '(a,i10)')   text(1:26), g%nv2

    if ( g%type == css ) then
      read (lunin, '(a,i10)') text(1:34), g%nzt
      read (lunin, '(a,i10)') text(1:35), g%nzs
    end if

    call memalloc(g%nv+1,g%ia,__FILE__,__LINE__)

    read (lunin, '(a)')       text
    do i=1,g%nv/10
       read (lunin, '(10i10)') g%ia( (i-1)*10+1 : i*10 )
    end do
    if(mod(g%nv,10)>0)  read (lunin, '(10i10)') g%ia( (i-1)*10+1 : g%nv+1 )

    if ( g%type == css ) then
      read (lunin, '(a)')     text
      do i=1,g%nv/10
         read (lunin, '(10i10)') g%is( (i-1)*10+1 : i*10 )
      end do
      if(mod(g%nv,10)>0)  read (lunin, '(10i10)') g%is( (i-1)*10+1 : g%nv+1 )
    end if

    call memalloc(g%ia(g%nv+1)-1,g%ja,__FILE__,__LINE__)

    read (lunin, '(a)')       text
    ! Only works when g%ia(i+1)-g%ia(i) < 10 !
    do i=1,g%nv
       read(lunin,'(10i10)') j, g%ja(g%ia(i):g%ia(i+1)-1)
       if(j/=i) then
          write(*,*) 'ERROR reading graph'
          stop
       end if
    end do
    read (lunin, '(a)')       text

  end subroutine fem_graph_read

  !=============================================================================
  subroutine fem_graph_ja_one_to_zero_indexing (g)
    implicit none
    ! Parameters
    type(fem_graph),  intent(inout) :: g
    integer(ip)                 :: j

    do j=g%ia(1), g%ia(g%nv+1)-1
      g%ja(j) = g%ja(j) - 1
    end do

  end subroutine fem_graph_ja_one_to_zero_indexing

  !=============================================================================
  subroutine fem_graph_ja_zero_to_one_indexing (g)
    implicit none
    ! Local Variables
    type(fem_graph),  intent(inout) :: g
    integer(ip)                 :: j

    do j=g%ia(1), g%ia(g%nv+1)-1
      g%ja(j) = g%ja(j) + 1
    end do

  end subroutine fem_graph_ja_zero_to_one_indexing

!!$  ! Inspired on http://en.wikipedia.org/wiki/Breadth-first_search
!!$  subroutine fem_graph_compute_connected_components (g, lconn)
!!$    implicit none
!!$
!!$    ! Parameters
!!$    type(fem_graph), intent(in)   :: g
!!$    type(list)     , intent(out)  :: lconn
!!$ 
!!$    ! Locals
!!$    integer(ip), allocatable :: aux(:)
!!$    integer(ip), allocatable :: marked(:)
!!$    integer(ip), allocatable :: q(:)
!!$    integer(ip)              :: head, tail, current, i, size, j, l, k
!!$
!!$    call memalloc( g%nv, aux     , __FILE__,__LINE__)
!!$    call memalloc( g%nv, q       , __FILE__,__LINE__)
!!$    call memalloc( g%nv, lconn%l , __FILE__,__LINE__)
!!$    call memalloc( g%nv, marked  , __FILE__,__LINE__)
!!$
!!$    lconn%n  = 0
!!$    marked   = 0
!!$    current  = 1 
!!$
!!$    do i=1, g%nv
!!$       if (marked(i) == 0) then
!!$          ! New connected component
!!$          lconn%n = lconn%n +1
!!$          size = 0
!!$1  procedure BFS(G,v):
!!$2      create a queue Q
!!$          head=1
!!$          tail=1
!!$3      enqueue v onto Q
!!$          q(tail)=i
!!$          tail=tail+1
!!$4      mark v
!!$          marked(i)=1
!!$          lconn%l(current)=i
!!$          size  = size + 1
!!$          current = current + 1  
!!$5      while Q is not empty:
!!$          do while (head/=tail)
!!$6         t ← Q.dequeue()
!!$             j=q(head)
!!$             head = head + 1
!!$9         for all edges e in G.adjacentEdges(t) do
!!$             do k=g%ia(j), g%ia(j+1)-1
!!$12           u ← G.adjacentVertex(t,e)
!!$                l=g%ja(k)
!!$13           if u is not marked:
!!$                if (marked(l)==0) then
!!$14              mark u
!!$                   marked(l)=1
!!$                   lconn%l(current)=l
!!$                   size  = size + 1
!!$                   current = current + 1  
!!$15              enqueue u onto Q
!!$                   q(tail)=l
!!$                   tail=tail+1
!!$                end if
!!$             end do
!!$          end do
!!$          aux(lconn%n) = size
!!$       end if
!!$    end do
!!$    
!!$    call memalloc( lconn%n+1, lconn%p, __FILE__,__LINE__)
!!$
!!$
!!$    lconn%p(1) = 1
!!$    do i=1, lconn%n
!!$        lconn%p(i+1) = lconn%p(i) + aux(i)
!!$     end do
!!$
!!$     write(*,*) 'ZZ', g%nv
!!$     write(*,*) 'XX', lconn%n
!!$     write(*,*) 'YY', lconn%p
!!$       
!!$    call memfree( aux,__FILE__,__LINE__)
!!$    call memfree( q  ,__FILE__,__LINE__)
!!$    call memfree( marked,__FILE__,__LINE__)
!!$
!!$  end subroutine fem_graph_compute_connected_components

  !=============================================================================
  !subroutine permut_graph(np,np2,nz,ia,ja,is,lrord,lrord2)
  !  !-----------------------------------------------------------------------
  !  ! This routine 
  !  !-----------------------------------------------------------------------
  !  implicit none
  !  integer(ip), intent(in)    :: np,np2,nz
  !  integer(ip), intent(in)    :: lrord(np),lrord2(np2)
  !  integer(ip), intent(inout) :: ia(np+1),is(np+1),ja(nz)
  !  integer(ip), allocatable   :: iaw(:),jaw(:)
  !  integer(ip)                :: i,j,k,iz,istat

  !  ! Reorder ia
  !  allocate(iaw(np+1),stat=istat)
  !  do i=1,np
  !     iaw(lrord(i)+1)=ia(i+1)-ia(i)
  !  end do
  !  iaw(1)=1
  !  do i=1,np
  !     iaw(i+1)=iaw(i+1)+iaw(i)
  !  end do

  !  ! Save ja
  !  allocate(jaw(nz),stat=istat)
  !  do i=1,np
  !     jaw(iaw(lrord(i)):iaw(lrord(i)+1)-1)=lrord2(ja(ia(i):ia(i+1)-1))
  !  end do

  !  ! Copy iaw and deallocate
  !  ia=iaw
  !  deallocate(iaw)

  !  ! Reorder ja and rebuild is
  !  is=0
  !  do i=1,np
  !     k=0
  !     do iz=ia(i),ia(i+1)-1
  !        if(jaw(iz)<=i) then
  !           ja(ia(i)+k)=jaw(iz)
  !           k=k+1
  !        else
  !           is(i+1)=is(i+1)+1
  !           ja(ia(i+1)-is(i+1))=jaw(iz)
  !        end if
  !     end do
  !  end do
  !  deallocate(jaw)

  !  ! Compress is
  !  is(1)=1
  !  do i=1,np
  !     is(i+1)=is(i+1)+is(i)
  !  end do

  !end subroutine permut_graph

end module fem_graph_class
