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
module mesh_graph_partition
  use types
  use memor
  use fem_mesh_names
  use fem_graph_names
  use mesh_graph
  use fem_partition_names
  implicit none
# include "debug.i90"
  private

  interface mesh_to_graph  !_part
    module procedure mesh_to_graph_matrix
  end interface

  ! Functions
  public :: mesh_to_graph

contains

 ! I am not sure that the computation of all possible graph
 ! types addressing sparse matrices can be encapsulated in one
 ! routine/interface. I will assume that this is possible. 
 ! If not, we will need separate subroutines for each graph type.
  subroutine mesh_to_graph_matrix ( storage, gtype, ndof1, ndof2, part, primal_mesh, primal_graph  )
    implicit none

    ! Parameters
    integer(ip)         , intent(in)  :: storage, gtype, ndof1, ndof2
    type(fem_partition) , intent(in)  :: part
    type(fem_mesh)      , intent(in)  :: primal_mesh
    type(fem_graph)     , intent(out) :: primal_graph

    ! Local variables
    type (fem_mesh)                   :: dual_mesh
    integer(ip), allocatable          :: iwork(:)       ! Integer ip working array
    integer(ip)                       :: pwork(3)       ! Pointers to work space

    ! Valid storage layouts for fem_matrix are blk and scal
    assert  ( storage == blk .or. storage == scal )

    assert  ( ndof1 >=1 .and. ndof2 >= 1 )

    ! Valid fem_graph types for addressing sparse matrices 
    ! are csr_symm or csr or css
    assert( gtype == csr_symm .or. gtype == csr .or. gtype == css )

    ! For vertex-based distributions, the only valid format
    ! is csr, as the local portion of each processor has to
    ! be rectangular (of course if nparts > 1)
    assert( part%ptype == vertex_based .or.  part%ptype == element_based )
    assert( .not. (part%ptype == vertex_based) .or. gtype==csr )

    ! Crs_symm is valid only for square matrices 
    assert  ( (.not. (gtype==csr_symm))  .or. (ndof1 == ndof2) )

    ! Compute dual_mesh
    call mesh_to_dual ( primal_mesh, dual_mesh )

    ! Allocate space for ia on the primal graph
    primal_graph%type = gtype

    if ( storage == blk ) then 
       ! Parallel (info about data-distribution extracted from partition object)
       if ( primal_graph%type == csr .or. primal_graph%type == csr_symm ) then
          if ( part%ptype == vertex_based ) then
             primal_graph%nv  = part%nmap%ni + part%nmap%nb
             primal_graph%nv2 = part%nmap%ni + part%nmap%nb + part%nmap%ne
          else 
             ! Element-based distribution
             primal_graph%nv  = part%nmap%nl
             primal_graph%nv2 = primal_graph%nv
          end if
          call memalloc ( primal_graph%nv+1, primal_graph%ia, __FILE__,__LINE__ )
       else 
          if  ( primal_graph%type == css ) then
             ! Element-based distribution
             primal_graph%nv  = part%nmap%nl
             primal_graph%nv2 = primal_graph%nv
             call memalloc ( primal_graph%nv+1, primal_graph%ia, __FILE__,__LINE__ )
             call memalloc ( primal_graph%nv+1, primal_graph%is, __FILE__,__LINE__ )
          end if
          ! other graph types may go here in an else construct (if required)
       end if

       ! Allocate working space for count_primal_graph and list_primal_graph routines
       ! (TOTAL WS SIZE = primal mesh npoin + maximum number of neighbours of any primal graph node)
       pwork(1) = 1
       pwork(2) = pwork(1) + primal_mesh%npoin
       pwork(3) = pwork(2) + dual_mesh%nnode*primal_mesh%nnode
       call memalloc ( pwork(3), iwork, __FILE__,__LINE__ )


       ! I am not sure that the computation of the number of neighbours
       ! can be encapsulated in one routine comprising all possible graph
       ! types addressing sparse matrices (e.g., do they require the same amount
       ! of workspace ?) I will assume that this is possible. If not, we will need 
       ! separate subroutines for each graph type.

       ! After a lot of thinking, I have decided to provide separate routines
       ! for count and list for each kind of graph. The body of the (old)
       ! common routine (see comment above) became very dirty. Besides it had
       ! a lot of if,else statements and conditional logic within the most deep
       ! loop (so that potentially a lot of overhead). I am aware that this
       ! decision will cause some controverse, but I am not sure which is the
       ! best solution at the moment ...

       if ( primal_graph%type == csr ) then
          call count_primal_graph_csr ( primal_mesh, dual_mesh, primal_graph, &  
               iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
       else
          if ( primal_graph%type == csr_symm ) then
             call count_primal_graph_csr_symm ( primal_mesh, dual_mesh, primal_graph, &  
                  iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
          else
             if ( primal_graph%type == css ) then
                call count_primal_graph_css ( primal_mesh, dual_mesh, primal_graph, &  
                     iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
             end if
          end if
       end if

       ! Allocate space for ja on the primal graph 
       call memalloc (primal_graph%ia(primal_graph%nv+1)-1, primal_graph%ja,         __FILE__,__LINE__)    

       ! I am not sure that the computation of the list of neighbours
       ! can be encapsulated in one routine comprising all possible graph
       ! types addressing sparse matrices (e.g., do they require the same amount
       ! of workspace ?). I will assume that this is possible. If not, we will need
       ! separate subroutines for each graph type.

       ! After a lot of thinking, I have decided to provide separate routines
       ! for count and list for each kind of graph. The body of the (old)
       ! common routine (see comment above) became very dirty. Besides it had
       ! a lot of if,else statements and conditional logic within the most deep
       ! loop (so that potentially a lot of overhead). I am aware that this
       ! decision will cause some controverse, but I am not sure which is the
       ! best solution at the moment ...
       if ( primal_graph%type == csr ) then
          call list_primal_graph_csr  ( primal_mesh, dual_mesh, primal_graph, &
               iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
       else
          if ( primal_graph%type == csr_symm ) then
             call list_primal_graph_csr_symm  ( primal_mesh, dual_mesh, primal_graph, &
                  iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
          else
             if ( primal_graph%type == css ) then
                call list_primal_graph_css  ( primal_mesh, dual_mesh, primal_graph, &
                     iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
             end if
          end if
       end if
    else if ( storage == scal ) then

       ! Parallel (info about data-distribution extracted from partition object)
       if ( primal_graph%type == csr .or. primal_graph%type == csr_symm ) then
          if ( part%ptype == vertex_based ) then
             primal_graph%nv  = (part%nmap%ni + part%nmap%nb)                * ndof1
             primal_graph%nv2 = (part%nmap%ni + part%nmap%nb + part%nmap%ne) * ndof2
          else 
             ! Element-based distribution
             primal_graph%nv  = part%nmap%nl    * ndof1
             primal_graph%nv2 = primal_graph%nv * ndof2
          end if
          call memalloc ( primal_graph%nv+1, primal_graph%ia, __FILE__,__LINE__ )
       else 
          if  ( primal_graph%type == css ) then
             ! Element-based distribution
             primal_graph%nv  = part%nmap%nl    * ndof1
             primal_graph%nv2 = primal_graph%nv * ndof2 
             call memalloc ( primal_graph%nv+1, primal_graph%ia, __FILE__,__LINE__ )
             call memalloc ( primal_graph%nv+1, primal_graph%is, __FILE__,__LINE__ )
          end if
          ! other graph types may go here in an else construct (if required)
       end if

       ! Allocate working space for count_primal_graph and list_primal_graph routines
       ! (TOTAL WS SIZE = primal mesh npoin + maximum number of neighbours of any primal graph node)
       pwork(1) = 1
       pwork(2) = pwork(1) + primal_mesh%npoin
       pwork(3) = pwork(2) + dual_mesh%nnode*primal_mesh%nnode
       call memalloc ( pwork(3), iwork, __FILE__,__LINE__ )


       ! I am not sure that the computation of the number of neighbours
       ! can be encapsulated in one routine comprising all possible graph
       ! types addressing sparse matrices (e.g., do they require the same amount
       ! of workspace ?) I will assume that this is possible. If not, we will need 
       ! separate subroutines for each graph type.

       ! After a lot of thinking, I have decided to provide separate routines
       ! for count and list for each kind of graph. The body of the (old)
       ! common routine (see comment above) became very dirty. Besides it had
       ! a lot of if,else statements and conditional logic within the most deep
       ! loop (so that potentially a lot of overhead). I am aware that this
       ! decision will cause some controverse, but I am not sure which is the
       ! best solution at the moment ...

       if ( primal_graph%type == csr ) then
          call count_primal_graph_csr_scal ( ndof1, ndof2, primal_mesh, dual_mesh, primal_graph, &  
                                             iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
       else
          if ( primal_graph%type == csr_symm ) then
             call count_primal_graph_csr_symm_scal ( ndof1, ndof2, primal_mesh, dual_mesh, primal_graph, &  
                                                     iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
          else
             if ( primal_graph%type == css ) then
                call count_primal_graph_css_scal ( ndof1, ndof2, primal_mesh, dual_mesh, primal_graph, &  
                                                   iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
             end if
          end if
       end if

       ! Allocate space for ja on the primal graph 
       call memalloc (primal_graph%ia(primal_graph%nv+1)-1, primal_graph%ja,      __FILE__,__LINE__)    

       ! I am not sure that the computation of the list of neighbours
       ! can be encapsulated in one routine comprising all possible graph
       ! types addressing sparse matrices (e.g., do they require the same amount
       ! of workspace ?). I will assume that this is possible. If not, we will need
       ! separate subroutines for each graph type.

       ! After a lot of thinking, I have decided to provide separate routines
       ! for count and list for each kind of graph. The body of the (old)
       ! common routine (see comment above) became very dirty. Besides it had
       ! a lot of if,else statements and conditional logic within the most deep
       ! loop (so that potentially a lot of overhead). I am aware that this
       ! decision will cause some controverse, but I am not sure which is the
       ! best solution at the moment ...

       if ( primal_graph%type == csr ) then
          call list_primal_graph_csr_scal  ( ndof1, ndof2, primal_mesh, dual_mesh, primal_graph, &
                                             iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
       else
          if ( primal_graph%type == csr_symm ) then
             call list_primal_graph_csr_symm_scal  ( ndof1, ndof2, primal_mesh, dual_mesh, primal_graph, &
                                                     iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
          else
             if ( primal_graph%type == css ) then
                call list_primal_graph_css_scal  ( ndof1, ndof2, primal_mesh, dual_mesh, primal_graph, &
                                                   iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
             end if
          end if
       end if

    end if

    ! Free dual_mesh
    call fem_mesh_free ( dual_mesh )

  end subroutine mesh_to_graph_matrix

end module mesh_graph_partition
