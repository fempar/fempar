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
module mesh_graph
  use types
  use memor
  use fem_mesh_class
  use fem_graph_class
  use renum_class
  implicit none
# include "debug.i90"
  private

  interface mesh_to_graph
     module procedure mesh_to_graph_part, mesh_to_graph_part_new, mesh_to_graph_share, mesh_to_graph_matrix, mesh_to_graph_matrix_sz
  end interface mesh_to_graph

  ! Functions
  public :: mesh_to_graph,                                                     &
            count_primal_graph_csr, list_primal_graph_csr,                     &
            count_primal_graph_csr_symm, list_primal_graph_csr_symm,           &
            count_primal_graph_css, list_primal_graph_css,                     & 
            count_primal_graph_csr_scal, list_primal_graph_csr_scal,           &
            count_primal_graph_csr_symm_scal, list_primal_graph_csr_symm_scal, &
            count_primal_graph_css_scal, list_primal_graph_css_scal          , &
            mesh_graph_compute_connected_components

! *** IMPORTANT NOTE: These cpp macros should go to a
! common include file or they should be a program 
! subroutine otherwise
#define blk2scal(iv,idof,ndof) (((iv)-1)*(ndof)+(idof))
#define scal2blk(id,idof,ndof) (((id)-1)/(ndof)+1)

contains
  !============================================================================================
  !
  !============================================================================================
  subroutine mesh_to_graph_matrix ( storage, gtype, ndof1, ndof2, primal_mesh, primal_graph,stabi)
    implicit none

    ! Parameters
    integer(ip)          , intent(in)     :: storage, gtype, ndof1, ndof2
    integer(ip), optional, intent(in)     :: stabi
    type(fem_mesh)       , intent(in)     :: primal_mesh
    type(fem_graph)      , intent(out)    :: primal_graph


    ! Local variables
    type (fem_mesh)                        :: dual_mesh
    integer(ip), allocatable               :: iwork(:)       ! Integer ip working array
    integer(ip)                            :: pwork(3)       ! Pointers to work space

    ! Valid storage layouts for fem_matrix are blk and scal
    assert  ( storage == blk .or. storage == scal )

    assert  ( ndof1 >= 1 .and. ndof2 >= 1 )

    ! Valid fem_graph types for addressing sparse matrices are csr_symm or csr or csc
    assert( gtype==csr_symm .or. gtype==csr .or. gtype == css )

    ! Crs_symm is valid only for square matrices 
    assert  ( (.not. (gtype==csr_symm))  .or. (ndof1 == ndof2) )

    ! Compute dual_mesh
    call mesh_to_dual ( primal_mesh, dual_mesh )

    ! Allocate space for ia on the primal graph
    primal_graph%type = gtype

    if ( storage == blk ) then
       if ( primal_graph%type == csr .or. primal_graph%type == csr_symm ) then
          primal_graph%nv  = primal_mesh%npoin
          primal_graph%nv2 = primal_graph%nv
          call memalloc ( primal_graph%nv+1, primal_graph%ia, __FILE__,__LINE__ )
       else
          if ( primal_graph%type == css ) then
             primal_graph%nv  = primal_mesh%npoin
             primal_graph%nv2 = primal_graph%nv
             call memalloc ( primal_graph%nv+1, primal_graph%ia, __FILE__,__LINE__ )
             call memalloc ( primal_graph%nv+1, primal_graph%is, __FILE__,__LINE__ )
          end if
          ! other fem_graph types may go here in an else construct (if required)
       end if

       ! Allocate working space for count_primal_graph and list_primal_graph routines
       ! (TOTAL WS SIZE = primal mesh npoin + maximum number of neighbours of any primal 
       ! graph node)
       pwork(1) = 1
       pwork(2) = pwork(1) + primal_mesh%npoin
!!$       if (present(stabi) .and. stabi == 3) then
!!$          pwork(3) = pwork(2) + (dual_mesh%nnode*primal_mesh%nnode)**3
!!$       else
!!$          pwork(3) = pwork(2) + dual_mesh%nnode*primal_mesh%nnode
!!$       end if
       if (present(stabi)) then
          if (stabi ==3) then
            pwork(3) = pwork(2) + dual_mesh%nnode*primal_mesh%nnode*dual_mesh%nnode*primal_mesh%nnode
          else
            pwork(3) = pwork(2) + dual_mesh%nnode*primal_mesh%nnode
          end if 
       else
          pwork(3) = pwork(2) + dual_mesh%nnode*primal_mesh%nnode
       end if

       call memalloc ( pwork(3), iwork, __FILE__,__LINE__ )


       ! I am not sure that the computation of the number of neighbours
       ! can be encapsulated in one routine comprising all possible fem_graph
       ! types addressing sparse matrices (e.g., do they require the same amount
       ! of workspace ?) I will assume that this is possible. If not, we will need 
       ! separate subroutines for each fem_graph type.

       ! After a lot of thinking, I have decided to provide separate routines
       ! for count and list for each kind of fem_graph. The body of the (old)
       ! common routine (see comment above) became very dirty. Besides it had
       ! a lot of if,else statements and conditional logic within the most deep
       ! loop (so that potentially a lot of overhead). I am aware that this
       ! decision will cause some controverse, but I am not sure which is the
       ! best solution at the moment ...
       if ( primal_graph%type == csr ) then
          if (present(stabi)) then
             if ( stabi == 3 ) then
                call count_primal_graph_extended_stencil_csr ( primal_mesh, dual_mesh, primal_graph, &  
                     iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
             else
                call count_primal_graph_csr ( primal_mesh, dual_mesh, primal_graph, &  
                     iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
             end if
          else
             call count_primal_graph_csr ( primal_mesh, dual_mesh, primal_graph, &  
                  iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
          end if
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
       call memalloc (primal_graph%ia(primal_graph%nv+1)-1, primal_graph%ja,          __FILE__,__LINE__)    

       ! I am not sure that the computation of the list of neighbours
       ! can be encapsulated in one routine comprising all possible fem_graph
       ! types addressing sparse matrices (e.g., do they require the same amount
       ! of workspace ?). I will assume that this is possible. If not, we will need
       ! separate subroutines for each fem_graph type.

       ! After a lot of thinking, I have decided to provide separate routines
       ! for count and list for each kind of fem_graph. The body of the (old)
       ! common routine (see comment above) became very dirty. Besides it had
       ! a lot of if,else statements and conditional logic within the most deep
       ! loop (so that potentially a lot of overhead). I am aware that this
       ! decision will cause some controverse, but I am not sure which is the
       ! best solution at the moment ...
       if ( primal_graph%type == csr ) then
          if (present(stabi)) then
             if (stabi == 3) then
                call list_primal_graph_extended_stencil_csr  ( primal_mesh, dual_mesh, primal_graph, &
                     iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
             else
                call list_primal_graph_csr  ( primal_mesh, dual_mesh, primal_graph, &
                     iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
             end if
          else
             call list_primal_graph_csr  ( primal_mesh, dual_mesh, primal_graph, &
                  iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
          end if
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

       if ( primal_graph%type == csr .or. primal_graph%type == csr_symm ) then
          primal_graph%nv  = primal_mesh%npoin * ndof1
          primal_graph%nv2 = primal_mesh%npoin * ndof2
          call memalloc ( primal_graph%nv+1, primal_graph%ia, __FILE__,__LINE__ )
       else
          if ( primal_graph%type == css ) then
             primal_graph%nv  = primal_mesh%npoin * ndof1
             primal_graph%nv2 = primal_mesh%npoin * ndof2
             call memalloc ( primal_graph%nv+1, primal_graph%ia, __FILE__,__LINE__ )
             call memalloc ( primal_graph%nv+1, primal_graph%is, __FILE__,__LINE__ )
          end if
          ! other fem_graph types may go here in an else construct (if required)
       end if

       ! Allocate working space for count_primal_graph and list_primal_graph routines
       ! (TOTAL WS SIZE = primal mesh npoin + maximum number of neighbours of any primal 
       ! graph node)
       pwork(1) = 1
       pwork(2) = pwork(1) + primal_mesh%npoin
!!$       if (present(stabi) .and. stabi == 3) then
!!$          pwork(3) = pwork(2) + dual_mesh%nnode*primal_mesh%nnode*dual_mesh%nnode*primal_mesh%nnode
!!$       else
!!$          pwork(3) = pwork(2) + dual_mesh%nnode*primal_mesh%nnode
!!$       end if
       if (present(stabi)) then
          if (stabi ==3) then
            pwork(3) = pwork(2) + dual_mesh%nnode*primal_mesh%nnode*dual_mesh%nnode*primal_mesh%nnode
          else
            pwork(3) = pwork(2) + dual_mesh%nnode*primal_mesh%nnode
          end if 
       else
          pwork(3) = pwork(2) + dual_mesh%nnode*primal_mesh%nnode
       end if

       call memalloc ( pwork(3), iwork, __FILE__,__LINE__ )

       ! I am not sure that the computation of the number of neighbours
       ! can be encapsulated in one routine comprising all possible fem_graph
       ! types addressing sparse matrices (e.g., do they require the same amount
       ! of workspace ?) I will assume that this is possible. If not, we will need 
       ! separate subroutines for each fem_graph type.

       ! After a lot of thinking, I have decided to provide separate routines
       ! for count and list for each kind of fem_graph. The body of the (old)
       ! common routine (see comment above) became very dirty. Besides it had
       ! a lot of if,else statements and conditional logic within the most deep
       ! loop (so that potentially a lot of overhead). I am aware that this
       ! decision will cause some controverse, but I am not sure which is the
       ! best solution at the moment ...

       if ( primal_graph%type == csr ) then
          if (present(stabi)) then
             if (stabi == 3) then
                call count_primal_graph_extended_stencil_csr_scal ( ndof1, ndof2, primal_mesh, dual_mesh, &
                     primal_graph,iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
             else
                call count_primal_graph_csr_scal ( ndof1, ndof2, primal_mesh, dual_mesh, primal_graph, &  
                     iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
             end if
           else
              call count_primal_graph_csr_scal ( ndof1, ndof2, primal_mesh, dual_mesh, primal_graph, &  
                                             iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
           end if
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

!!$    write (*,*) primal_graph%ia

       ! Allocate space for ja on the primal graph 
       call memalloc ( primal_graph%ia(primal_graph%nv+1)-1, primal_graph%ja,          __FILE__,__LINE__)    

       ! I am not sure that the computation of the list of neighbours
       ! can be encapsulated in one routine comprising all possible fem_graph
       ! types addressing sparse matrices (e.g., do they require the same amount
       ! of workspace ?). I will assume that this is possible. If not, we will need
       ! separate subroutines for each fem_graph type.

       ! After a lot of thinking, I have decided to provide separate routines
       ! for count and list for each kind of fem_graph. The body of the (old)
       ! common routine (see comment above) became very dirty. Besides it had
       ! a lot of if,else statements and conditional logic within the most deep
       ! loop (so that potentially a lot of overhead). I am aware that this
       ! decision will cause some controverse, but I am not sure which is the
       ! best solution at the moment ...
       if ( primal_graph%type == csr ) then
          if (present(stabi)) then
             if (stabi == 3) then
                call list_primal_graph_extended_stencil_csr_scal  ( ndof1, ndof2, primal_mesh, dual_mesh, & 
                     primal_graph,iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
             else
                call list_primal_graph_csr_scal  ( ndof1, ndof2, primal_mesh, dual_mesh, primal_graph, &
                  iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
             end if
          else
             call list_primal_graph_csr_scal  ( ndof1, ndof2, primal_mesh, dual_mesh, primal_graph, &
                  iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
          end if
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
    call memfree ( iwork,__FILE__,__LINE__)

  end subroutine mesh_to_graph_matrix

!============================================================================================
  subroutine mesh_to_graph_matrix_sz ( storage, gtype, ndof1, ndof2, primal_mesh, primal_graph,szmap)
    implicit none

    ! Parameters
    integer(ip)          , intent(in)     :: storage, gtype, ndof1, ndof2
    integer(ip)          , intent(in)     :: szmap(:)
    type(fem_mesh)       , intent(in)     :: primal_mesh
    type(fem_graph)      , intent(out)    :: primal_graph


    ! Local variables
    type (fem_mesh)                        :: dual_mesh
    integer(ip), allocatable               :: iwork(:)       ! Integer ip working array
    integer(ip)                            :: pwork(3)       ! Pointers to work space
    integer(ip)                            :: szdual_ia(primal_mesh%nelem+1)
    integer(ip), allocatable               :: aux_dual(:),szdual_ja(:)
    
    !SZ is only implemented for the CSR case
    assert ( gtype == csr )

    ! Valid storage layouts for fem_matrix are blk and scal
    assert  ( storage == blk .or. storage == scal )
    assert  ( ndof1 >= 1 .and. ndof2 >= 1 )

    ! Valid fem_graph types for addressing sparse matrices are csr_symm or csr or csc
    assert( gtype==csr_symm .or. gtype==csr .or. gtype == css )

    ! Crs_symm is valid only for square matrices 
    assert  ( (.not. (gtype==csr_symm))  .or. (ndof1 == ndof2) )

    ! Compute dual_mesh
    call mesh_to_dual ( primal_mesh, dual_mesh )

    ! Compute dual SZ map
    call memalloc(primal_mesh%nelem*primal_mesh%nnode,aux_dual,__FILE__,__LINE__)
    call szmap_to_dual(szmap,primal_mesh%npoin,primal_mesh%nelem,szdual_ia,aux_dual)
    pwork(1) = szdual_ia(primal_mesh%nelem+1)
    call memalloc(pwork(1),szdual_ja,__FILE__,__LINE__ )
    szdual_ja = aux_dual(1:pwork(1))
    call memfree(aux_dual,__FILE__,__LINE__)

    ! Allocate space for ia on the primal graph
    primal_graph%type = gtype

    if ( storage == blk ) then
       primal_graph%nv  = primal_mesh%npoin
       primal_graph%nv2 = primal_graph%nv
       call memalloc ( primal_graph%nv+1, primal_graph%ia, __FILE__,__LINE__ )

       ! Allocate working space for count_primal_graph and list_primal_graph routines
       ! (TOTAL WS SIZE = primal mesh npoin + maximum number of neighbours of any primal 
       ! graph node)
       pwork(1) = 1
       pwork(2) = pwork(1) + primal_mesh%npoin
       pwork(3) = pwork(2) + (dual_mesh%nnode*primal_mesh%nnode)**3
       call memalloc ( pwork(3), iwork, __FILE__,__LINE__ )
 
       ! Store number of neighbour nodes in primal graph
       call count_primal_graph_sz_csr ( primal_mesh, dual_mesh, primal_graph,szmap, &  
            szdual_ia,szdual_ja,iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
       
       ! Allocate space for ja on the primal graph 
       call memalloc (primal_graph%ia(primal_graph%nv+1)-1, primal_graph%ja,           __FILE__,__LINE__)    

       call list_primal_graph_sz_csr  ( primal_mesh, dual_mesh, primal_graph, szmap, &
            szdual_ia,szdual_ja,iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
      
    else if ( storage == scal ) then
       
       primal_graph%nv  = primal_mesh%npoin * ndof1
       primal_graph%nv2 = primal_mesh%npoin * ndof2
       call memalloc ( primal_graph%nv+1, primal_graph%ia, __FILE__,__LINE__ )
       
       ! Allocate working space for count_primal_graph and list_primal_graph routines
       ! (TOTAL WS SIZE = primal mesh npoin + maximum number of neighbours of any primal 
       ! graph node)
       pwork(1) = 1
       pwork(2) = pwork(1) + primal_mesh%npoin
       pwork(3) = pwork(2) + dual_mesh%nnode*primal_mesh%nnode*dual_mesh%nnode*primal_mesh%nnode
       
       call memalloc ( pwork(3), iwork, __FILE__,__LINE__ )

       call count_primal_graph_sz_csr_scal ( ndof1, ndof2, primal_mesh, dual_mesh, &
            primal_graph,szmap,szdual_ia,szdual_ja,iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
      
       ! Allocate space for ja on the primal graph 
       call memalloc ( primal_graph%ia(primal_graph%nv+1)-1, primal_graph%ja,           __FILE__,__LINE__)    

       call list_primal_graph_sz_csr_scal  ( ndof1, ndof2, primal_mesh, dual_mesh, &
            primal_graph, szmap,szdual_ia,szdual_ja,iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
    end if

    ! Free dual_mesh
    call fem_mesh_free ( dual_mesh )


  end subroutine mesh_to_graph_matrix_sz

  !============================================================================================
  subroutine  count_primal_graph_csr ( primal_mesh, dual_mesh, primal_graph,  &
       &                                 ws_position, ws_neighbors )
    implicit none

    ! Parameters
    type(fem_mesh) , intent(in)     :: primal_mesh, dual_mesh
    type(fem_graph), intent(inout)  :: primal_graph
    integer(ip), intent(out)        :: ws_position  (primal_mesh%npoin)
    integer(ip), intent(out)        :: ws_neighbors (primal_mesh%nnode*dual_mesh%nnode)

    ! Local variables
    integer(ip)                     :: ineigh  
    integer(ip)                     :: ipoindm ! Dual   mesh point identifier
    integer(ip)                     :: ipoinpm ! Primal mesh point identifier
    integer(ip)                     :: ipoinpg ! Primal graph point identifier 
    integer(ip)                     :: p_ipoindm
    integer(ip)                     :: p_ipoinpm
    integer(ip)                     :: inods1,inods2   

    integer(ip)                     :: first_free_pos
    integer(ip)                     :: num_neighbors

    integer(ip)                     :: npoin

    assert(dual_mesh%nelty/=1)

    ! Initialize work space of filled-up 
    ! positions
    ws_position   = 0

    ! Initialize work space of neighbour 
    ! identifiers 
    ws_neighbors  = 0

    !primal_graph%ia    = -1                                                ! DBG:
    !write (*,*) size(primal_graph%ia)                                      ! DBG: 
    !write (*,*) size(primal_mesh%pnods)                                    ! DBG: 
    !write (*,*) size(primal_mesh%lnods)                                    ! DBG: 
    !write (*,*) size(dual_mesh%pnods)                                      ! DBG: 
    !write (*,*) size(dual_mesh%lnods)                                      ! DBG: 
    !write (*,*) size(ws_position)                                          ! DBG: 
    !write (*,*) size(ws_freq)                                              ! DBG: 
    !write (*,*) size(ws_neighbors)                                         ! DBG: 
    !write (*, '(10i10)')    primal_mesh%pnods  (1:(primal_mesh%nelem+1))   ! DBG: 
    !write (*, '(10i10)')    primal_graph%ia    (1:(primal_graph%nv+1))     ! DBG:  
    !primal_graph%ia    = -1 
    primal_graph%ia(1) = 1

    do ipoinpg=1, primal_graph%nv
       primal_graph%ia(ipoinpg+1) = primal_graph%ia(ipoinpg)
       first_free_pos = 1
       num_neighbors  = 0

       ! Traverse the dual nodes of the dual element number ipoinpg
       ! (i.e., primal elements around primal node number ipoinpg)
       do p_ipoindm = dual_mesh%pnods(ipoinpg), dual_mesh%pnods(ipoinpg+1)-1
          ipoindm = dual_mesh%lnods(p_ipoindm)

          ! Traverse the primal nodes of the primal element number ipoindm
          ! (i.e., dual elements around dual node number ipoindm)
          if(primal_mesh%nelty==1) then
             inods1=(ipoindm-1)*primal_mesh%nnode+1
             inods2=ipoindm*primal_mesh%nnode
          else
             inods1=primal_mesh%pnods(ipoindm)
             inods2=primal_mesh%pnods(ipoindm+1)-1
          end if
          do p_ipoinpm = inods1,inods2
             ipoinpm = primal_mesh%lnods(p_ipoinpm)
             ! write (*,*) ipoinpg, ipoindm, ipoinpm ! DBG: 
             ! If ipoinpm not visited yet
             if ( ws_position(ipoinpm) == 0 ) then
                ws_position(ipoinpm) = first_free_pos
                ws_neighbors(first_free_pos) = ipoinpm
                first_free_pos = first_free_pos + 1
                num_neighbors = num_neighbors + 1 
                ! If ipoinpm already visited
             end if
          end do
       end do

       ! write (*, '(10i10)')    ws_position (1:primal_graph%nv)  ! DBG:
       ! write (*, '(10i10)')    ws_neighbors(1:num_neighbors)    ! DBG: 
       ! write (*, '(10i10)')    ws_freq     (1:num_neighbors)    ! DBG: 

       ! Extract neighbours while restoring working space to initial state
       do ineigh=1, num_neighbors
          primal_graph%ia(ipoinpg+1) = primal_graph%ia(ipoinpg+1) + 1
          ws_position ( ws_neighbors(ineigh) ) = 0
          ws_neighbors(ineigh) = 0
       end do
    end do

    primal_graph%nzt  = primal_graph%ia(primal_graph%nv+1) - primal_graph%ia(1)

  end subroutine count_primal_graph_csr


!============================================================================================
  subroutine  count_primal_graph_extended_stencil_csr ( primal_mesh, dual_mesh, primal_graph,  &
       &                                 ws_position, ws_neighbors )
    implicit none

    ! Parameters
    type(fem_mesh) , intent(in)     :: primal_mesh, dual_mesh
    type(fem_graph), intent(inout)  :: primal_graph
    integer(ip), intent(out)    :: ws_position  (primal_mesh%npoin)
    integer(ip), intent(out)    :: ws_neighbors(primal_mesh%nnode*dual_mesh%nnode*primal_mesh%nnode*dual_mesh%nnode)

    ! Local variables
    integer(ip)                     :: ineigh  
    integer(ip)                     :: ipoindm, jpoindm, kpoindm ! Dual   mesh point identifier
    integer(ip)                     :: ipoinpm, jpoinpm, kpoinpm ! Primal mesh point identifier
    integer(ip)                     :: ipoinpg ! Primal graph point identifier 
    integer(ip)                     :: p_ipoindm, p_jpoindm,p_kpoindm
    integer(ip)                     :: p_ipoinpm, p_jpoinpm,p_kpoinpm
    integer(ip)                     :: inods1,inods2,jnods1,jnods2,knods1,knods2

    integer(ip)                     :: first_free_pos
    integer(ip)                     :: num_neighbors

    integer(ip)                     :: npoin

    assert(dual_mesh%nelty/=1)

    ! Initialize work space of filled-up 
    ! positions
    ws_position   = 0

    ! Initialize work space of neighbour 
    ! identifiers 
    ws_neighbors  = 0

    !primal_graph%ia    = -1                                                ! DBG:
    !write (*,*) size(primal_graph%ia)                                      ! DBG: 
    !write (*,*) size(primal_mesh%pnods)                                    ! DBG: 
    !write (*,*) size(primal_mesh%lnods)                                    ! DBG: 
    !write (*,*) size(dual_mesh%pnods)                                      ! DBG: 
    !write (*,*) size(dual_mesh%lnods)                                      ! DBG: 
    !write (*,*) size(ws_position)                                          ! DBG: 
    !write (*,*) size(ws_freq)                                              ! DBG: 
    !write (*,*) size(ws_neighbors)                                         ! DBG: 
    !write (*, '(10i10)')    primal_mesh%pnods  (1:(primal_mesh%nelem+1))   ! DBG: 
    !write (*, '(10i10)')    primal_graph%ia    (1:(primal_graph%nv+1))     ! DBG:  
    !primal_graph%ia    = -1 
    primal_graph%ia(1) = 1

    do ipoinpg=1, primal_graph%nv
       primal_graph%ia(ipoinpg+1) = primal_graph%ia(ipoinpg)
       first_free_pos = 1
       num_neighbors  = 0

       ! Traverse the dual nodes of the dual element number ipoinpg
       ! (i.e., primal elements around primal node number ipoinpg)
       do p_ipoindm = dual_mesh%pnods(ipoinpg), dual_mesh%pnods(ipoinpg+1)-1
          ipoindm = dual_mesh%lnods(p_ipoindm)

          ! Traverse the primal nodes of the primal element number ipoindm
          ! (i.e., dual elements around dual node number ipoindm)
          if(primal_mesh%nelty==1) then
             inods1=(ipoindm-1)*primal_mesh%nnode+1
             inods2=ipoindm*primal_mesh%nnode
          else
             inods1=primal_mesh%pnods(ipoindm)
             inods2=primal_mesh%pnods(ipoindm+1)-1
          end if
          do p_ipoinpm = inods1,inods2
             ipoinpm = primal_mesh%lnods(p_ipoinpm)

             do p_jpoindm = dual_mesh%pnods(ipoinpm), dual_mesh%pnods(ipoinpm+1)-1
                jpoindm = dual_mesh%lnods(p_jpoindm)

                ! Traverse the primal nodes of the primal element number ipoindm
                ! (i.e., dual elements around dual node number ipoindm)
                if(primal_mesh%nelty==1) then
                   jnods1=(jpoindm-1)*primal_mesh%nnode+1
                   jnods2=jpoindm*primal_mesh%nnode
                else
                   jnods1=primal_mesh%pnods(jpoindm)
                   jnods2=primal_mesh%pnods(jpoindm+1)-1
                end if
                do p_jpoinpm = jnods1,jnods2
                   jpoinpm = primal_mesh%lnods(p_jpoinpm)

                   do p_kpoindm = dual_mesh%pnods(jpoinpm), dual_mesh%pnods(jpoinpm+1)-1
                      kpoindm = dual_mesh%lnods(p_kpoindm)

                      ! Traverse the primal nodes of the primal element number ipoindm
                      ! (i.e., dual elements around dual node number ipoindm)
                      if(primal_mesh%nelty==1) then
                         knods1=(kpoindm-1)*primal_mesh%nnode+1
                         knods2=kpoindm*primal_mesh%nnode
                      else
                         knods1=primal_mesh%pnods(kpoindm)
                         knods2=primal_mesh%pnods(kpoindm+1)-1
                      end if
                      do p_kpoinpm = knods1,knods2
                         kpoinpm = primal_mesh%lnods(p_kpoinpm)
                         ! write (*,*) ipoinpg, ipoindm, ipoinpm ! DBG: 
                         ! If ipoinpm not visited yet
                         if ( ws_position(kpoinpm) == 0 ) then
                            ws_position(kpoinpm) = first_free_pos
                            ws_neighbors(first_free_pos) = kpoinpm
                            first_free_pos = first_free_pos + 1
                            num_neighbors = num_neighbors + 1 
                         end if
                      end do
                   end do
                   
                end do
             end do
          end do
       end do

       ! write (*, '(10i10)')    ws_position (1:primal_graph%nv)  ! DBG:
       ! write (*, '(10i10)')    ws_neighbors(1:num_neighbors)    ! DBG: 
       ! write (*, '(10i10)')    ws_freq     (1:num_neighbors)    ! DBG: 

       ! Extract neighbours while restoring working space to initial state
       do ineigh=1, num_neighbors
          primal_graph%ia(ipoinpg+1) = primal_graph%ia(ipoinpg+1) + 1
          ws_position ( ws_neighbors(ineigh) ) = 0
          ws_neighbors(ineigh) = 0
       end do
    end do

    primal_graph%nzt  = primal_graph%ia(primal_graph%nv+1) - primal_graph%ia(1)

  end subroutine count_primal_graph_extended_stencil_csr

 !============================================================================================
  subroutine  count_primal_graph_sz_csr ( primal_mesh, dual_mesh, primal_graph,szmap,  &
       &                                 szdual_ia,szdual_ja,ws_position, ws_neighbors )
    implicit none

    ! Parameters
    type(fem_mesh) , intent(in)     :: primal_mesh, dual_mesh
    integer(ip)    , intent(in)     :: szmap(:),szdual_ia(:),szdual_ja(:)
    type(fem_graph), intent(inout)  :: primal_graph
    integer(ip), intent(out)    :: ws_position  (primal_mesh%npoin)
    integer(ip), intent(out)    :: ws_neighbors(primal_mesh%nnode*dual_mesh%nnode*primal_mesh%nnode*dual_mesh%nnode)

    ! Local variables
    integer(ip)                     :: ineigh  
    integer(ip)                     :: ipoindm, jpoindm, kpoindm ! Dual   mesh point identifier
    integer(ip)                     :: ipoinpm, jpoinpm, kpoinpm ! Primal mesh point identifier
    integer(ip)                     :: ipoinpg ! Primal graph point identifier 
    integer(ip)                     :: p_ipoindm, p_jpoindm,p_kpoindm
    integer(ip)                     :: p_ipoinpm, p_jpoinpm,p_kpoinpm
    integer(ip)                     :: inods1,inods2,jnods1,jnods2,knods1,knods2

    integer(ip)                     :: first_free_pos
    integer(ip)                     :: num_neighbors

    integer(ip)                     :: npoin

    assert(dual_mesh%nelty/=1)

    ! Initialize work space of filled-up positions
    ws_position   = 0

    ! Initialize work space of neighbour identifiers 
    ws_neighbors  = 0

    primal_graph%ia(1) = 1

    do ipoinpg=1, primal_graph%nv
       primal_graph%ia(ipoinpg+1) = primal_graph%ia(ipoinpg)
       first_free_pos = 1
       num_neighbors  = 0

       ! Projection vs Projection terms
       do p_ipoindm = dual_mesh%pnods(ipoinpg), dual_mesh%pnods(ipoinpg+1)-1
          ipoindm = dual_mesh%lnods(p_ipoindm) 

          !Nodes ipoinpm such that szmap(ipoinpm) = ipoindm
          inods1=szdual_ia(ipoindm)
          inods2=szdual_ia(ipoindm+1)-1
          do p_ipoinpm = inods1,inods2
             ipoinpm = szdual_ja(p_ipoinpm)
             
             !Elems jpoindm such that ipoinpm is in jpoindm
             do p_jpoindm = dual_mesh%pnods(ipoinpm), dual_mesh%pnods(ipoinpm+1)-1
                jpoindm = dual_mesh%lnods(p_jpoindm)

                !Nodes jpoinpm such that jpoinpm is in jpoindm
                if(primal_mesh%nelty==1) then
                   jnods1=(jpoindm-1)*primal_mesh%nnode+1
                   jnods2=jpoindm*primal_mesh%nnode
                else
                   jnods1=primal_mesh%pnods(jpoindm)
                   jnods2=primal_mesh%pnods(jpoindm+1)-1
                end if
                do p_jpoinpm = jnods1,jnods2
                   jpoinpm = primal_mesh%lnods(p_jpoinpm)

                   !Elems kpoindm = szmap(jpoinpm)
                   kpoindm = szmap(jpoinpm)

                   !Nodes kpoinpm such that kpoinpm is in kpoindm
                   if(primal_mesh%nelty==1) then
                      knods1=(kpoindm-1)*primal_mesh%nnode+1
                      knods2=kpoindm*primal_mesh%nnode
                   else
                      knods1=primal_mesh%pnods(kpoindm)
                      knods2=primal_mesh%pnods(kpoindm+1)-1
                   end if
                   do p_kpoinpm = knods1,knods2
                      kpoinpm = primal_mesh%lnods(p_kpoinpm)
                      ! write (*,*) ipoinpg, ipoindm, ipoinpm ! DBG: 
                      ! If ipoinpm not visited yet
                      if ( ws_position(kpoinpm) == 0 ) then
                         ws_position(kpoinpm) = first_free_pos
                         ws_neighbors(first_free_pos) = kpoinpm
                         first_free_pos = first_free_pos + 1
                         num_neighbors = num_neighbors + 1 
                      end if
                   end do
                end do
             end do
          end do
       end do

 
       ! Extract neighbours while restoring working space to initial state
       do ineigh=1, num_neighbors
          primal_graph%ia(ipoinpg+1) = primal_graph%ia(ipoinpg+1) + 1
          ws_position ( ws_neighbors(ineigh) ) = 0
          ws_neighbors(ineigh) = 0
       end do
    end do

    primal_graph%nzt  = primal_graph%ia(primal_graph%nv+1) - primal_graph%ia(1)

  end subroutine count_primal_graph_sz_csr

  !============================================================================================
  subroutine  count_primal_graph_csr_scal ( ndof1, ndof2, primal_mesh, dual_mesh, primal_graph,  &
       &                                    ws_position, ws_neighbors )
    implicit none

    ! Parameters
    integer(ip)    , intent(in)     :: ndof1, ndof2
    type(fem_mesh) , intent(in)     :: primal_mesh, dual_mesh
    type(fem_graph), intent(inout)  :: primal_graph
    integer(ip)    , intent(out)    :: ws_position  (primal_mesh%npoin)
    integer(ip)    , intent(out)    :: ws_neighbors (primal_mesh%nnode*dual_mesh%nnode)

    ! Local variables
    integer(ip)                     :: ineigh  
    integer(ip)                     :: ipoindm ! Dual   mesh point identifier
    integer(ip)                     :: ipoinpm ! Primal mesh point identifier
    integer(ip)                     :: ipoinpg ! Primal graph point identifier 
    integer(ip)                     :: p_ipoindm
    integer(ip)                     :: p_ipoinpm
    integer(ip)                     :: inods1,inods2   

    integer(ip)                     :: first_free_pos
    integer(ip)                     :: num_neighbors

    integer(ip)                     :: npoin, idof1

    assert(dual_mesh%nelty/=1)

    ! Initialize work space of filled-up 
    ! positions
    ws_position   = 0

    ! Initialize work space of neighbour 
    ! identifiers 
    ws_neighbors  = 0

    !primal_graph%ia    = -1                                                ! DBG:
    !write (*,*) size(primal_graph%ia)                                      ! DBG: 
    !write (*,*) size(primal_mesh%pnods)                                    ! DBG: 
    !write (*,*) size(primal_mesh%lnods)                                    ! DBG: 
    !write (*,*) size(dual_mesh%pnods)                                      ! DBG: 
    !write (*,*) size(dual_mesh%lnods)                                      ! DBG: 
    !write (*,*) size(ws_position)                                          ! DBG: 
    !write (*,*) size(ws_freq)                                              ! DBG: 
    !write (*,*) size(ws_neighbors)                                         ! DBG: 
    !write (*, '(10i10)')    primal_mesh%pnods  (1:(primal_mesh%nelem+1))   ! DBG: 
    !write (*, '(10i10)')    primal_graph%ia    (1:(primal_graph%nv+1))     ! DBG:  
    !primal_graph%ia    = -1 
    primal_graph%ia(1) = 1

    do ipoinpg=1, primal_graph%nv, ndof1
       primal_graph%ia(ipoinpg+1) = primal_graph%ia(ipoinpg)
       ! write (*,*) 'ZZ', ipoinpg, scal2blk(ipoinpg,1,ndof1), primal_graph%ia(ipoinpg+1)

       first_free_pos = 1
       num_neighbors  = 0

       ! Traverse the dual nodes of the dual element number ipoinpg
       ! (i.e., primal elements around primal node number ipoinpg)
       do p_ipoindm = dual_mesh%pnods(scal2blk(ipoinpg,1,ndof1)), dual_mesh%pnods(scal2blk(ipoinpg,1,ndof1)+1)-1
          ipoindm = dual_mesh%lnods(p_ipoindm)

          ! Traverse the primal nodes of the primal element number ipoindm
          ! (i.e., dual elements around dual node number ipoindm)
          if(primal_mesh%nelty==1) then
             inods1=(ipoindm-1)*primal_mesh%nnode+1
             inods2=ipoindm*primal_mesh%nnode
          else
             inods1=primal_mesh%pnods(ipoindm)
             inods2=primal_mesh%pnods(ipoindm+1)-1
          end if
          do p_ipoinpm = inods1,inods2
             ipoinpm = primal_mesh%lnods(p_ipoinpm)
             ! write (*,*) ipoinpg, ipoindm, ipoinpm ! DBG: 
             ! If ipoinpm not visited yet
             if ( ws_position(ipoinpm) == 0 ) then
                ws_position(ipoinpm) = first_free_pos
                ws_neighbors(first_free_pos) = ipoinpm
                first_free_pos = first_free_pos + 1
                num_neighbors = num_neighbors + 1 
                ! If ipoinpm already visited
             end if
          end do
       end do

       !write (*, '(a10i10)')    'XX', ws_position (1:primal_graph%nv/ndof1)  ! DBG:
       !write (*, '(a10i10)')    'XX', ws_neighbors(1:num_neighbors)    ! DBG: 

       ! write (*,*) 'XX', ipoinpg, primal_graph%ia(ipoinpg+1), ndof2, num_neighbors

       ! Extract neighbours while restoring working space to initial state
       do ineigh=1, num_neighbors
          primal_graph%ia(ipoinpg+1) = primal_graph%ia(ipoinpg+1) + ndof2
          ws_position ( ws_neighbors(ineigh) ) = 0
          ws_neighbors(ineigh) = 0
       end do

       do idof1 = 2, ndof1
          primal_graph%ia(ipoinpg+idof1) = primal_graph%ia(ipoinpg+idof1-1)                         &  ! Base 
                                         & + (primal_graph%ia(ipoinpg+1) - primal_graph%ia(ipoinpg))   ! Length 
       end do
    end do

    primal_graph%nzt  = primal_graph%ia(primal_graph%nv+1) - primal_graph%ia(1)

  end subroutine count_primal_graph_csr_scal





  !============================================================================================
  subroutine  count_primal_graph_extended_stencil_csr_scal ( ndof1, ndof2, primal_mesh, dual_mesh, primal_graph,  &
       &                                    ws_position, ws_neighbors )
    implicit none

    ! Parameters
    integer(ip)    , intent(in)     :: ndof1, ndof2
    type(fem_mesh) , intent(in)     :: primal_mesh, dual_mesh
    type(fem_graph), intent(inout)  :: primal_graph
    integer(ip)    , intent(out)    :: ws_position  (primal_mesh%npoin)
    integer(ip)    , intent(out)    :: ws_neighbors (primal_mesh%nnode**2*dual_mesh%nnode**2)

    ! Local variables
    integer(ip)                     :: ineigh  
    integer(ip)                     :: ipoindm,jpoindm,kpoindm ! Dual   mesh point identifier
    integer(ip)                     :: ipoinpm,jpoinpm,kpoinpm ! Primal mesh point identifier
    integer(ip)                     :: ipoinpg ! Primal graph point identifier 
    integer(ip)                     :: p_ipoindm, p_jpoindm, p_kpoindm
    integer(ip)                     :: p_ipoinpm, p_jpoinpm, p_kpoinpm
    integer(ip)                     :: inods1,inods2,jnods1,jnods2 ,knods1,knods2   

    integer(ip)                     :: first_free_pos
    integer(ip)                     :: num_neighbors

    integer(ip)                     :: npoin, idof1

    assert(dual_mesh%nelty/=1)

    ! Initialize work space of filled-up 
    ! positions
    ws_position   = 0

    ! Initialize work space of neighbour 
    ! identifiers 
    ws_neighbors  = 0

    !primal_graph%ia    = -1                                                ! DBG:
    !write (*,*) size(primal_graph%ia)                                      ! DBG: 
    !write (*,*) size(primal_mesh%pnods)                                    ! DBG: 
    !write (*,*) size(primal_mesh%lnods)                                    ! DBG: 
    !write (*,*) size(dual_mesh%pnods)                                      ! DBG: 
    !write (*,*) size(dual_mesh%lnods)                                      ! DBG: 
    !write (*,*) size(ws_position)                                          ! DBG: 
    !write (*,*) size(ws_freq)                                              ! DBG: 
    !write (*,*) size(ws_neighbors)                                         ! DBG: 
    !write (*, '(10i10)')    primal_mesh%pnods  (1:(primal_mesh%nelem+1))   ! DBG: 
    !write (*, '(10i10)')    primal_graph%ia    (1:(primal_graph%nv+1))     ! DBG:  
    !primal_graph%ia    = -1 
    primal_graph%ia(1) = 1

    do ipoinpg=1, primal_graph%nv, ndof1
       primal_graph%ia(ipoinpg+1) = primal_graph%ia(ipoinpg)
       ! write (*,*) 'ZZ', ipoinpg, scal2blk(ipoinpg,1,ndof1), primal_graph%ia(ipoinpg+1)

       first_free_pos = 1
       num_neighbors  = 0

       ! Traverse the dual nodes of the dual element number ipoinpg
       ! (i.e., primal elements around primal node number ipoinpg)
       do p_ipoindm = dual_mesh%pnods(scal2blk(ipoinpg,1,ndof1)), dual_mesh%pnods(scal2blk(ipoinpg,1,ndof1)+1)-1
          ipoindm = dual_mesh%lnods(p_ipoindm)

          ! Traverse the primal nodes of the primal element number ipoindm
          ! (i.e., dual elements around dual node number ipoindm)
          if(primal_mesh%nelty==1) then
             inods1=(ipoindm-1)*primal_mesh%nnode+1
             inods2=ipoindm*primal_mesh%nnode
          else
             inods1=primal_mesh%pnods(ipoindm)
             inods2=primal_mesh%pnods(ipoindm+1)-1
          end if
          do p_ipoinpm = inods1,inods2
             ipoinpm = primal_mesh%lnods(p_ipoinpm)

             do p_jpoindm = dual_mesh%pnods(scal2blk(ipoinpm,1,ndof1)), &
                  dual_mesh%pnods(scal2blk(ipoinpm,1,ndof1)+1)-1
                jpoindm = dual_mesh%lnods(p_jpoindm)

                ! Traverse the primal nodes of the primal element number ipoindm
                ! (i.e., dual elements around dual node number ipoindm)
                if(primal_mesh%nelty==1) then
                   jnods1=(jpoindm-1)*primal_mesh%nnode+1
                   jnods2=jpoindm*primal_mesh%nnode
                else
                   jnods1=primal_mesh%pnods(jpoindm)
                   jnods2=primal_mesh%pnods(jpoindm+1)-1
                end if
                do p_jpoinpm = jnods1,jnods2
                   jpoinpm = primal_mesh%lnods(p_jpoinpm)
                   
                   do p_kpoindm = dual_mesh%pnods(jpoinpm), dual_mesh%pnods(jpoinpm+1)-1
                      kpoindm = dual_mesh%lnods(p_kpoindm)

                      ! Traverse the primal nodes of the primal element number ipoindm
                      ! (i.e., dual elements around dual node number ipoindm)
                      if(primal_mesh%nelty==1) then
                         knods1=(kpoindm-1)*primal_mesh%nnode+1
                         knods2=kpoindm*primal_mesh%nnode
                      else
                         knods1=primal_mesh%pnods(kpoindm)
                         knods2=primal_mesh%pnods(kpoindm+1)-1
                      end if
                      do p_kpoinpm = knods1,knods2
                         kpoinpm = primal_mesh%lnods(p_kpoinpm)
                         ! write (*,*) ipoinpg, ipoindm, ipoinpm ! DBG: 
                         ! If ipoinpm not visited yet
                         if ( ws_position(kpoinpm) == 0 ) then
                            ws_position(kpoinpm) = first_free_pos
                            ws_neighbors(first_free_pos) = kpoinpm
                            first_free_pos = first_free_pos + 1
                            num_neighbors = num_neighbors + 1 
                         end if
                      end do
                   end do
                end do
             end do
          end do
       end do

       !write (*, '(a10i10)')    'XX', ws_position (1:primal_graph%nv/ndof1)  ! DBG:
       !write (*, '(a10i10)')    'XX', ws_neighbors(1:num_neighbors)    ! DBG: 

       ! write (*,*) 'XX', ipoinpg, primal_graph%ia(ipoinpg+1), ndof2, num_neighbors

       ! Extract neighbours while restoring working space to initial state
       do ineigh=1, num_neighbors
          primal_graph%ia(ipoinpg+1) = primal_graph%ia(ipoinpg+1) + ndof2
          ws_position ( ws_neighbors(ineigh) ) = 0
          ws_neighbors(ineigh) = 0
       end do

       do idof1 = 2, ndof1
          primal_graph%ia(ipoinpg+idof1) = primal_graph%ia(ipoinpg+idof1-1)                         &  ! Base 
                                         & + (primal_graph%ia(ipoinpg+1) - primal_graph%ia(ipoinpg))   ! Length 
       end do
    end do

    primal_graph%nzt  = primal_graph%ia(primal_graph%nv+1) - primal_graph%ia(1)

  end subroutine count_primal_graph_extended_stencil_csr_scal

 !============================================================================================
  subroutine  count_primal_graph_sz_csr_scal ( ndof1,ndof2,primal_mesh, dual_mesh, primal_graph,szmap,  &
       &                                 szdual_ia,szdual_ja,ws_position, ws_neighbors )
    implicit none

    ! Parameters
    integer(ip)    , intent(in)     :: ndof1, ndof2
    type(fem_mesh) , intent(in)     :: primal_mesh, dual_mesh
    integer(ip)    , intent(in)     :: szmap(:),szdual_ia(:),szdual_ja(:)
    type(fem_graph), intent(inout)  :: primal_graph
    integer(ip), intent(out)    :: ws_position  (primal_mesh%npoin)
    integer(ip), intent(out)    :: ws_neighbors(primal_mesh%nnode*dual_mesh%nnode*primal_mesh%nnode*dual_mesh%nnode)

    ! Local variables
    integer(ip)                     :: ineigh  
    integer(ip)                     :: ipoindm, jpoindm, kpoindm ! Dual   mesh point identifier
    integer(ip)                     :: ipoinpm, jpoinpm, kpoinpm ! Primal mesh point identifier
    integer(ip)                     :: ipoinpg ! Primal graph point identifier 
    integer(ip)                     :: p_ipoindm, p_jpoindm,p_kpoindm
    integer(ip)                     :: p_ipoinpm, p_jpoinpm,p_kpoinpm
    integer(ip)                     :: inods1,inods2,jnods1,jnods2,knods1,knods2

    integer(ip)                     :: first_free_pos
    integer(ip)                     :: num_neighbors

    integer(ip)                     :: npoin, idof1


    assert(dual_mesh%nelty/=1)

    ! Initialize work space of filled-up positions
    ws_position   = 0

    ! Initialize work space of neighbour identifiers 
    ws_neighbors  = 0

    primal_graph%ia(1) = 1

    do ipoinpg=1, primal_graph%nv, ndof1
       primal_graph%ia(ipoinpg+1) = primal_graph%ia(ipoinpg)
       first_free_pos = 1
       num_neighbors  = 0

       ! Projection vs Projection terms
       do p_ipoindm = dual_mesh%pnods(scal2blk(ipoinpg,1,ndof1)), dual_mesh%pnods(scal2blk(ipoinpg,1,ndof1)+1)-1
          ipoindm = dual_mesh%lnods(p_ipoindm) 

          !Nodes ipoinpm such that szmap(ipoinpm) = ipoindm
          inods1=szdual_ia(ipoindm)
          inods2=szdual_ia(ipoindm+1)-1
          do p_ipoinpm = inods1,inods2
             ipoinpm = szdual_ja(p_ipoinpm)
             
             !Elems jpoindm such that ipoinpm is in jpoindm
             do p_jpoindm = dual_mesh%pnods(scal2blk(ipoinpm,1,ndof1)), dual_mesh%pnods(scal2blk(ipoinpm,1,ndof1)+1)-1
                jpoindm = dual_mesh%lnods(p_jpoindm)

                !Nodes jpoinpm such that jpoinpm is in jpoindm
                if(primal_mesh%nelty==1) then
                   jnods1=(jpoindm-1)*primal_mesh%nnode+1
                   jnods2=jpoindm*primal_mesh%nnode
                else
                   jnods1=primal_mesh%pnods(jpoindm)
                   jnods2=primal_mesh%pnods(jpoindm+1)-1
                end if
                do p_jpoinpm = jnods1,jnods2
                   jpoinpm = primal_mesh%lnods(p_jpoinpm)

                   !Elems kpoindm = szmap(jpoinpm)
                   kpoindm = szmap(jpoinpm)

                   !Nodes kpoinpm such that kpoinpm is in kpoindm
                   if(primal_mesh%nelty==1) then
                      knods1=(kpoindm-1)*primal_mesh%nnode+1
                      knods2=kpoindm*primal_mesh%nnode
                   else
                      knods1=primal_mesh%pnods(kpoindm)
                      knods2=primal_mesh%pnods(kpoindm+1)-1
                   end if
                   do p_kpoinpm = knods1,knods2
                      kpoinpm = primal_mesh%lnods(p_kpoinpm)
                      ! write (*,*) ipoinpg, ipoindm, ipoinpm ! DBG: 
                      ! If ipoinpm not visited yet
                      if ( ws_position(kpoinpm) == 0 ) then
                         ws_position(kpoinpm) = first_free_pos
                         ws_neighbors(first_free_pos) = kpoinpm
                         first_free_pos = first_free_pos + 1
                         num_neighbors = num_neighbors + 1 
                      end if
                   end do
                end do
             end do
          end do
       end do


 
       ! Extract neighbours while restoring working space to initial state
       do ineigh=1, num_neighbors
          primal_graph%ia(ipoinpg+1) = primal_graph%ia(ipoinpg+1) + ndof2
          ws_position ( ws_neighbors(ineigh) ) = 0
          ws_neighbors(ineigh) = 0
       end do

       do idof1 = 2, ndof1
          primal_graph%ia(ipoinpg+idof1) = primal_graph%ia(ipoinpg+idof1-1)                         &  ! Base 
                                         & + (primal_graph%ia(ipoinpg+1) - primal_graph%ia(ipoinpg))   ! Length 
       end do
    end do

    primal_graph%nzt  = primal_graph%ia(primal_graph%nv+1) - primal_graph%ia(1)

  end subroutine count_primal_graph_sz_csr_scal

  !============================================================================================
  subroutine  count_primal_graph_csr_symm ( primal_mesh, dual_mesh, primal_graph,  &
       &                                    ws_position, ws_neighbors)
    implicit none

    ! Parameters
    type(fem_mesh) , intent(in)            :: primal_mesh, dual_mesh
    type(fem_graph), intent(inout)         :: primal_graph
    integer(ip), intent(out)               :: ws_position  (primal_mesh%npoin)
    integer(ip), intent(out)               :: ws_neighbors (primal_mesh%nnode*dual_mesh%nnode)

    ! Local variables
    integer(ip)                     :: ineigh  
    integer(ip)                     :: ipoindm ! Dual   mesh point identifier
    integer(ip)                     :: ipoinpm ! Primal mesh point identifier
    integer(ip)                     :: ipoinpg ! Primal graph point identifier 
    integer(ip)                     :: p_ipoindm   
    integer(ip)                     :: p_ipoinpm   
    integer(ip)                     :: inods1,inods2   

    integer(ip)                     :: first_free_pos
    integer(ip)                     :: num_neighbors

    integer(ip)                     :: npoin

    assert(dual_mesh%nelty/=1)

    ! Initialize work space of filled-up 
    ! positions
    ws_position   = 0

    ! Initialize work space of neighbour 
    ! identifiers 
    ws_neighbors  = 0

    !primal_graph%ia    = -1                                                ! DBG:
    !write (*,*) size(primal_graph%ia)                                      ! DBG: 
    !write (*,*) size(primal_mesh%pnods)                                    ! DBG: 
    !write (*,*) size(primal_mesh%lnods)                                    ! DBG: 
    !write (*,*) size(dual_mesh%pnods)                                      ! DBG: 
    !write (*,*) size(dual_mesh%lnods)                                      ! DBG: 
    !write (*,*) size(ws_position)                                          ! DBG: 
    !write (*,*) size(ws_freq)                                              ! DBG: 
    !write (*,*) size(ws_neighbors)                                         ! DBG: 
    !write (*, '(10i10)')    primal_mesh%pnods  (1:(primal_mesh%nelem+1))   ! DBG: 
    !write (*, '(10i10)')    primal_graph%ia    (1:(primal_graph%nv+1))     ! DBG:  
    !primal_graph%ia    = -1 
    primal_graph%ia(1) = 1

    do ipoinpg=1, primal_graph%nv
       primal_graph%ia(ipoinpg+1) = primal_graph%ia(ipoinpg)
       first_free_pos = 1
       num_neighbors  = 0

       ! Traverse the dual nodes of the dual element number ipoinpg
       ! (i.e., primal elements around primal node number ipoinpg)
       do p_ipoindm = dual_mesh%pnods(ipoinpg), dual_mesh%pnods(ipoinpg+1)-1
          ipoindm = dual_mesh%lnods(p_ipoindm)

          ! Traverse the primal nodes of the primal element number ipoindm
          ! (i.e., dual elements around dual node number ipoindm)
          if(primal_mesh%nelty==1) then
             inods1=(ipoindm-1)*primal_mesh%nnode+1
             inods2=ipoindm*primal_mesh%nnode
          else
             inods1=primal_mesh%pnods(ipoindm)
             inods2=primal_mesh%pnods(ipoindm+1)-1
          end if
          do p_ipoinpm = inods1,inods2
             ipoinpm = primal_mesh%lnods(p_ipoinpm)
             ! Only count edges (i,j) s.t., i <= j  
             if ( ipoinpg <= ipoinpm ) then
                ! write (*,*) ipoinpg, ipoindm, ipoinpm ! DBG: 
                ! If ipoinpm not visited yet
                if ( ws_position(ipoinpm) == 0 ) then
                   ws_position(ipoinpm) = first_free_pos
                   ws_neighbors(first_free_pos) = ipoinpm
                   first_free_pos = first_free_pos + 1
                   num_neighbors = num_neighbors + 1 
                   ! If ipoinpm already visited
                end if
             end if
          end do
       end do

       ! write (*, '(10i10)')    ws_position (1:primal_graph%nv)  ! DBG:
       ! write (*, '(10i10)')    ws_neighbors(1:num_neighbors)    ! DBG: 
       ! write (*, '(10i10)')    ws_freq     (1:num_neighbors)    ! DBG: 

       ! Extract neighbours while restoring working space to initial state
       do ineigh=1, num_neighbors
          primal_graph%ia(ipoinpg+1) = primal_graph%ia(ipoinpg+1) + 1
          ws_position ( ws_neighbors(ineigh) ) = 0
          ws_neighbors(ineigh) = 0
       end do
    end do

    primal_graph%nzt  = primal_graph%ia(primal_graph%nv+1) - primal_graph%ia(1)

  end subroutine count_primal_graph_csr_symm

  !============================================================================================
  subroutine  count_primal_graph_csr_symm_scal ( ndof1, ndof2, primal_mesh, dual_mesh, primal_graph,  &
       &                                         ws_position, ws_neighbors)
    implicit none

    ! Parameters
    integer(ip)                            :: ndof1, ndof2
    type(fem_mesh) , intent(in)            :: primal_mesh, dual_mesh
    type(fem_graph), intent(inout)         :: primal_graph
    integer(ip), intent(out)               :: ws_position  (primal_mesh%npoin)
    integer(ip), intent(out)               :: ws_neighbors (primal_mesh%nnode*dual_mesh%nnode)

    ! Local variables
    integer(ip)                     :: ineigh  
    integer(ip)                     :: ipoindm ! Dual   mesh point identifier
    integer(ip)                     :: ipoinpm ! Primal mesh point identifier
    integer(ip)                     :: ipoinpg ! Primal graph point identifier 
    integer(ip)                     :: p_ipoindm   
    integer(ip)                     :: p_ipoinpm   
    integer(ip)                     :: inods1,inods2   

    integer(ip)                     :: first_free_pos
    integer(ip)                     :: num_neighbors

    integer(ip)                     :: npoin, idof1

    assert( dual_mesh%nelty/=1 )

    ! Initialize work space of filled-up 
    ! positions
    ws_position   = 0

    ! Initialize work space of neighbour 
    ! identifiers 
    ws_neighbors  = 0

    !primal_graph%ia    = -1                                                ! DBG:
    !write (*,*) size(primal_graph%ia)                                      ! DBG: 
    !write (*,*) size(primal_mesh%pnods)                                    ! DBG: 
    !write (*,*) size(primal_mesh%lnods)                                    ! DBG: 
    !write (*,*) size(dual_mesh%pnods)                                      ! DBG: 
    !write (*,*) size(dual_mesh%lnods)                                      ! DBG: 
    !write (*,*) size(ws_position)                                          ! DBG: 
    !write (*,*) size(ws_freq)                                              ! DBG: 
    !write (*,*) size(ws_neighbors)                                         ! DBG: 
    !write (*, '(10i10)')    primal_mesh%pnods  (1:(primal_mesh%nelem+1))   ! DBG: 
    !write (*, '(10i10)')    primal_graph%ia    (1:(primal_graph%nv+1))     ! DBG:  
    !primal_graph%ia    = -1 
    primal_graph%ia(1) = 1


    do ipoinpg=1, primal_graph%nv, ndof1
       primal_graph%ia(ipoinpg+1) = primal_graph%ia(ipoinpg)
       first_free_pos = 1
       num_neighbors  = 0

       ! Traverse the dual nodes of the dual element number ipoinpg
       ! (i.e., primal elements around primal node number ipoinpg)
       do p_ipoindm = dual_mesh%pnods(scal2blk(ipoinpg,1,ndof1)), dual_mesh%pnods(scal2blk(ipoinpg,1,ndof1)+1)-1
          ipoindm = dual_mesh%lnods(p_ipoindm)

          ! Traverse the primal nodes of the primal element number ipoindm
          ! (i.e., dual elements around dual node number ipoindm)
          if(primal_mesh%nelty==1) then
             inods1=(ipoindm-1)*primal_mesh%nnode+1
             inods2=ipoindm*primal_mesh%nnode
          else
             inods1=primal_mesh%pnods(ipoindm)
             inods2=primal_mesh%pnods(ipoindm+1)-1
          end if
          do p_ipoinpm = inods1,inods2
             ipoinpm = primal_mesh%lnods(p_ipoinpm)
             ! Only count edges (i,j) s.t., i <= j
             if ( scal2blk(ipoinpg,1,ndof1) <= ipoinpm ) then
                ! write (*,*) ipoinpg, ipoindm, ipoinpm ! DBG: 
                ! If ipoinpm not visited yet
                if ( ws_position(ipoinpm) == 0 ) then
                   ws_position(ipoinpm) = first_free_pos
                   ws_neighbors(first_free_pos) = ipoinpm
                   first_free_pos = first_free_pos + 1
                   num_neighbors = num_neighbors + 1 
                   ! If ipoinpm already visited
                end if
             end if
          end do
       end do

       ! write (*, '(10i10)')    ws_position (1:primal_graph%nv)  ! DBG:
       ! write (*, '(10i10)')    ws_neighbors(1:num_neighbors)    ! DBG: 
       ! write (*, '(10i10)')    ws_freq     (1:num_neighbors)    ! DBG: 

       ! Extract neighbours while restoring working space to initial state
       do ineigh=1, num_neighbors
          primal_graph%ia(ipoinpg+1) = primal_graph%ia(ipoinpg+1) + ndof2
          ws_position ( ws_neighbors(ineigh) ) = 0
          ws_neighbors(ineigh) = 0
       end do

       do idof1 = 2, ndof1
          primal_graph%ia(ipoinpg+idof1) = primal_graph%ia(ipoinpg+idof1-1)                            &  ! Base 
                                         & + (primal_graph%ia(ipoinpg+1) - primal_graph%ia(ipoinpg))   &  ! Length 
                                         & + (1-idof1)
       end do
    end do

    primal_graph%nzt  = primal_graph%ia(primal_graph%nv+1) - primal_graph%ia(1)

  end subroutine count_primal_graph_csr_symm_scal

  !============================================================================================
  subroutine  count_primal_graph_css ( primal_mesh, dual_mesh, primal_graph,  &
       &                               ws_position, ws_neighbors)
    implicit none

    ! Parameters
    type(fem_mesh) , intent(in)            :: primal_mesh, dual_mesh
    type(fem_graph), intent(inout)         :: primal_graph
    integer(ip), intent(out)               :: ws_position  (primal_mesh%npoin)
    integer(ip), intent(out)               :: ws_neighbors (primal_mesh%nnode*dual_mesh%nnode)

    ! Local variables
    integer(ip)                     :: ineigh  
    integer(ip)                     :: ipoindm ! Dual   mesh point identifier
    integer(ip)                     :: ipoinpm ! Primal mesh point identifier
    integer(ip)                     :: ipoinpg ! Primal graph point identifier 
    integer(ip)                     :: p_ipoindm   
    integer(ip)                     :: p_ipoinpm   
    integer(ip)                     :: inods1,inods2   

    integer(ip)                     :: first_free_pos
    integer(ip)                     :: num_neighbors_low
    integer(ip)                     :: num_neighbors_upp

    integer(ip)                     :: npoin

    assert(dual_mesh%nelty/=1)

    ! Initialize work space of filled-up 
    ! positions
    ws_position   = 0

    ! Initialize work space of neighbour 
    ! identifiers 
    ws_neighbors  = 0

    !primal_graph%ia    = -1                                                ! DBG:
    !write (*,*) size(primal_graph%ia)                                      ! DBG: 
    !write (*,*) size(primal_mesh%pnods)                                    ! DBG: 
    !write (*,*) size(primal_mesh%lnods)                                    ! DBG: 
    !write (*,*) size(dual_mesh%pnods)                                      ! DBG: 
    !write (*,*) size(dual_mesh%lnods)                                      ! DBG: 
    !write (*,*) size(ws_position)                                          ! DBG: 
    !write (*,*) size(ws_freq)                                              ! DBG: 
    !write (*,*) size(ws_neighbors)                                         ! DBG: 
    !write (*, '(10i10)')    primal_mesh%pnods  (1:(primal_mesh%nelem+1))   ! DBG: 
    !write (*, '(10i10)')    primal_graph%ia    (1:(primal_graph%nv+1))     ! DBG:  
    !primal_graph%ia    = -1 
    primal_graph%ia(1) = 1
    primal_graph%is(1) = 1


    do ipoinpg=1, primal_graph%nv
       primal_graph%ia(ipoinpg+1) = primal_graph%ia(ipoinpg)
       first_free_pos     = 1
       num_neighbors_low  = 0
       num_neighbors_upp  = 0

       ! Traverse the dual nodes of the dual element number ipoinpg
       ! (i.e., primal elements around primal node number ipoinpg)
       do p_ipoindm = dual_mesh%pnods(ipoinpg), dual_mesh%pnods(ipoinpg+1)-1
          ipoindm = dual_mesh%lnods(p_ipoindm)

          ! Traverse the primal nodes of the primal element number ipoindm
          ! (i.e., dual elements around dual node number ipoindm)
          if(primal_mesh%nelty==1) then
             inods1=(ipoindm-1)*primal_mesh%nnode+1
             inods2=ipoindm*primal_mesh%nnode
          else
             inods1=primal_mesh%pnods(ipoindm)
             inods2=primal_mesh%pnods(ipoindm+1)-1
          end if
          do p_ipoinpm = inods1,inods2
             ipoinpm = primal_mesh%lnods(p_ipoinpm)
             ! Only count edges (i,j) s.t., i > j (lower)  or edges s.t. i < j (upper)
             if ( (ipoinpg > ipoinpm) ) then ! Lower triangle entries
                ! If ipoinpm not visited yet
                if ( ws_position(ipoinpm) == 0 ) then
                   ws_position(ipoinpm) = first_free_pos
                   ws_neighbors(first_free_pos) = ipoinpm
                   first_free_pos = first_free_pos + 1
                   num_neighbors_low = num_neighbors_low + 1 
                   ! If ipoinpm already visited
                end if
             else
                if ( (ipoinpg < ipoinpm) ) then ! Upper triangle entries
                   ! If ipoinpm not visited yet
                   if ( ws_position(ipoinpm) == 0 ) then
                      ws_position(ipoinpm) = first_free_pos
                      ws_neighbors(first_free_pos) = ipoinpm
                      first_free_pos = first_free_pos + 1
                      num_neighbors_upp = num_neighbors_upp + 1 
                      ! If ipoinpm already visited
                   end if
                end if
             end if
          end do
       end do

       ! write (*, '(10i10)')    ws_position (1:primal_graph%nv)  ! DBG:
       ! write (*, '(10i10)')    ws_neighbors(1:num_neighbors)    ! DBG: 
       ! write (*, '(10i10)')    ws_freq     (1:num_neighbors)    ! DBG: 

       ! Extract those neighbours which are listed at least min_freq_neig times
       ! while restoring working space to initial state
       do ineigh=1, num_neighbors_low + num_neighbors_upp
          primal_graph%ia(ipoinpg+1) = primal_graph%ia(ipoinpg+1) + 1
          ws_position ( ws_neighbors(ineigh) ) = 0
          ws_neighbors(ineigh) = 0
       end do
       primal_graph%is(ipoinpg+1) = primal_graph%is(ipoinpg) + num_neighbors_upp
    end do

    primal_graph%nzt  = primal_graph%ia(primal_graph%nv+1) - primal_graph%ia(1)
    primal_graph%nzs = primal_graph%is(primal_graph%nv+1) - primal_graph%is(1)

  end subroutine count_primal_graph_css

  !============================================================================================
  subroutine  count_primal_graph_css_scal ( ndof1, ndof2, primal_mesh, dual_mesh, primal_graph,  &
       &                                    ws_position, ws_neighbors)
    implicit none

   
    ! Parameters
    integer(ip)    , intent(in)            :: ndof1, ndof2
    type(fem_mesh) , intent(in)            :: primal_mesh, dual_mesh
    type(fem_graph), intent(inout)         :: primal_graph
    integer(ip), intent(out)               :: ws_position  (primal_mesh%npoin)
    integer(ip), intent(out)               :: ws_neighbors (primal_mesh%nnode*dual_mesh%nnode)

    ! Local variables
    integer(ip)                     :: ineigh  
    integer(ip)                     :: ipoindm ! Dual   mesh point identifier
    integer(ip)                     :: ipoinpm ! Primal mesh point identifier
    integer(ip)                     :: ipoinpg ! Primal graph point identifier 
    integer(ip)                     :: p_ipoindm   
    integer(ip)                     :: p_ipoinpm   
    integer(ip)                     :: inods1,inods2   

    integer(ip)                     :: first_free_pos
    integer(ip)                     :: num_neighbors_low
    integer(ip)                     :: num_neighbors_upp

    integer(ip)                     :: npoin

    assert(dual_mesh%nelty/=1)

    write (0,*) 'Error: the body of count_primal_graph_css_scal in mesh_graph.f90 still to be written'
    write (0,*) 'Error: the body of count_primal_graph_css in mesh_graph.f90 is a good point to start'
    write (0,*) 'Error: volunteers are welcome !!!'
    
    stop

  end subroutine count_primal_graph_css_scal


  !============================================================================================
  subroutine  list_primal_graph_csr ( primal_mesh, dual_mesh, primal_graph,  &
       &                              ws_position, ws_neighbors )    
    use psb_sort_mod
    implicit none

    ! Parameters
    type(fem_mesh) , intent(in)    :: primal_mesh, dual_mesh
    type(fem_graph), intent(inout) :: primal_graph
    integer(ip), intent(out)       :: ws_position  (primal_mesh%npoin)
    integer(ip), intent(out)       :: ws_neighbors (primal_mesh%nnode*dual_mesh%nnode)

    ! Local variables
    integer(ip)                     :: ineigh  
    integer(ip)                     :: ipoindm ! Dual   mesh point identifier
    integer(ip)                     :: ipoinpm ! Primal mesh point identifier
    integer(ip)                     :: ipoinpg ! Primal graph point identifier 
    integer(ip)                     :: p_ipoindm   
    integer(ip)                     :: p_ipoinpm   
    integer(ip)                     :: inods1,inods2   

    integer(ip)                     :: first_free_ja_pos
    integer(ip)                     :: first_free_pos
    integer(ip)                     :: num_neighbors

    assert(dual_mesh%nelty/=1)

    ! Initialize work space of filled-up 
    ! positions to zeros
    ws_position   = 0

    ! Initialize working space of neighbour 
    ! identifiers 
    ws_neighbors = 0

    first_free_ja_pos = 1

    do ipoinpg=1, primal_graph%nv       
       first_free_pos = 1
       num_neighbors  = 0

       ! Traverse the dual nodes of the dual element number ipoinpg
       ! (i.e., primal elements around primal node number ipoinpg)
       do p_ipoindm = dual_mesh%pnods(ipoinpg), dual_mesh%pnods(ipoinpg+1)-1
          ipoindm = dual_mesh%lnods(p_ipoindm)

          ! Traverse the primal nodes of the primal element number ipoindm
          ! (i.e., dual elements around dual node number ipoindm)
          if(primal_mesh%nelty==1) then
             inods1=(ipoindm-1)*primal_mesh%nnode+1
             inods2=ipoindm*primal_mesh%nnode
          else
             inods1=primal_mesh%pnods(ipoindm)
             inods2=primal_mesh%pnods(ipoindm+1)-1
          end if
          do p_ipoinpm = inods1,inods2
             ipoinpm = primal_mesh%lnods(p_ipoinpm)
             ! If ipoinpm not visited yet
             if ( ws_position(ipoinpm) == 0 ) then
                ws_position(ipoinpm) = first_free_pos
                ws_neighbors(first_free_pos) = ipoinpm
                first_free_pos = first_free_pos + 1
                num_neighbors = num_neighbors + 1 
                ! If ipoinpm already visited
             end if
          end do
       end do

       ! Extract neighbours while restoring working space to initial state
       do ineigh=1, num_neighbors
          primal_graph%ja(first_free_ja_pos) = ws_neighbors(ineigh)
          first_free_ja_pos = first_free_ja_pos + 1
          ws_position ( ws_neighbors(ineigh) ) = 0
          ws_neighbors( ineigh               ) = 0
       end do

       ! Order increasingly column identifiers of current row 
       ! using heap sort algorithm
       call psb_hsort(primal_graph%ja( primal_graph%ia(ipoinpg):primal_graph%ia(ipoinpg+1)-1 ))

    end do

  end subroutine list_primal_graph_csr


!============================================================================================
  subroutine  list_primal_graph_extended_stencil_csr ( primal_mesh, dual_mesh, primal_graph,  &
       &                              ws_position, ws_neighbors )    
    use psb_sort_mod
    implicit none

    ! Parameters
    type(fem_mesh) , intent(in)    :: primal_mesh, dual_mesh
    type(fem_graph), intent(inout) :: primal_graph
    integer(ip), intent(out)       :: ws_position  (primal_mesh%npoin)
    integer(ip), intent(out)       :: ws_neighbors (primal_mesh%nnode**2*dual_mesh%nnode**2)

    ! Local variables
    integer(ip)                     :: ineigh  
    integer(ip)                     :: ipoindm, jpoindm, kpoindm ! Dual   mesh point identifier
    integer(ip)                     :: ipoinpm, jpoinpm, kpoinpm ! Primal mesh point identifier
    integer(ip)                     :: ipoinpg ! Primal graph point identifier 
    integer(ip)                     :: p_ipoindm, p_jpoindm, p_kpoindm   
    integer(ip)                     :: p_ipoinpm, p_jpoinpm, p_kpoinpm   
    integer(ip)                     :: inods1,inods2,jnods1,jnods2,knods1,knods2  

    integer(ip)                     :: first_free_ja_pos
    integer(ip)                     :: first_free_pos
    integer(ip)                     :: num_neighbors

    assert(dual_mesh%nelty/=1)

    ! Initialize work space of filled-up 
    ! positions to zeros
    ws_position   = 0

    ! Initialize working space of neighbour 
    ! identifiers 
    ws_neighbors = 0

    first_free_ja_pos = 1

    do ipoinpg=1, primal_graph%nv       
       first_free_pos = 1
       num_neighbors  = 0

       ! Traverse the dual nodes of the dual element number ipoinpg
       ! (i.e., primal elements around primal node number ipoinpg)
       do p_ipoindm = dual_mesh%pnods(ipoinpg), dual_mesh%pnods(ipoinpg+1)-1
          ipoindm = dual_mesh%lnods(p_ipoindm)

          ! Traverse the primal nodes of the primal element number ipoindm
          ! (i.e., dual elements around dual node number ipoindm)
          if(primal_mesh%nelty==1) then
             inods1=(ipoindm-1)*primal_mesh%nnode+1
             inods2=ipoindm*primal_mesh%nnode
          else
             inods1=primal_mesh%pnods(ipoindm)
             inods2=primal_mesh%pnods(ipoindm+1)-1
          end if
          do p_ipoinpm = inods1,inods2
             ipoinpm = primal_mesh%lnods(p_ipoinpm)

             do p_jpoindm = dual_mesh%pnods(ipoinpm), dual_mesh%pnods(ipoinpm+1)-1
               jpoindm = dual_mesh%lnods(p_jpoindm)

               ! Traverse the primal nodes of the primal element number ipoindm
               ! (i.e., dual elements around dual node number ipoindm)
               if(primal_mesh%nelty==1) then
                  jnods1=(jpoindm-1)*primal_mesh%nnode+1
                  jnods2=jpoindm*primal_mesh%nnode
               else
                  jnods1=primal_mesh%pnods(jpoindm)
                  jnods2=primal_mesh%pnods(jpoindm+1)-1
               end if
               do p_jpoinpm = jnods1,jnods2
                  jpoinpm = primal_mesh%lnods(p_jpoinpm)
                  ! If ipoinpm not visited yet
                  do p_kpoindm = dual_mesh%pnods(jpoinpm), dual_mesh%pnods(jpoinpm+1)-1
                     kpoindm = dual_mesh%lnods(p_kpoindm)

                     ! Traverse the primal nodes of the primal element number ipoindm
                     ! (i.e., dual elements around dual node number ipoindm)
                     if(primal_mesh%nelty==1) then
                        knods1=(kpoindm-1)*primal_mesh%nnode+1
                        knods2=kpoindm*primal_mesh%nnode
                     else
                        knods1=primal_mesh%pnods(kpoindm)
                        knods2=primal_mesh%pnods(kpoindm+1)-1
                     end if
                     do p_kpoinpm = knods1,knods2
                        kpoinpm = primal_mesh%lnods(p_kpoinpm)
                        ! If ipoinpm not visited yet
                        if ( ws_position(kpoinpm) == 0 ) then
                           ws_position(kpoinpm) = first_free_pos
                           ws_neighbors(first_free_pos) = kpoinpm
                           first_free_pos = first_free_pos + 1
                           num_neighbors = num_neighbors + 1 
                           ! If ipoinpm already visited
                        end if
                     end do
                  end do
                  
               end do
            end do
          end do
       end do

       ! Extract neighbours while restoring working space to initial state
       do ineigh=1, num_neighbors
          primal_graph%ja(first_free_ja_pos) = ws_neighbors(ineigh)
          first_free_ja_pos = first_free_ja_pos + 1
          ws_position ( ws_neighbors(ineigh) ) = 0
          ws_neighbors( ineigh               ) = 0
       end do

       ! Order increasingly column identifiers of current row 
       ! using heap sort algorithm
       call psb_hsort(primal_graph%ja( primal_graph%ia(ipoinpg):primal_graph%ia(ipoinpg+1)-1 ))

    end do

  end subroutine list_primal_graph_extended_stencil_csr


!============================================================================================
  subroutine  list_primal_graph_sz_csr ( primal_mesh, dual_mesh, primal_graph,szmap,  &
       &                              szdual_ia,szdual_ja,ws_position, ws_neighbors )    
    use psb_sort_mod
    implicit none

    ! Parameters
    type(fem_mesh) , intent(in)    :: primal_mesh, dual_mesh
    integer(ip)    , intent(in)     :: szmap(:),szdual_ia(:),szdual_ja(:)
    type(fem_graph), intent(inout) :: primal_graph
    integer(ip), intent(out)       :: ws_position  (primal_mesh%npoin)
    integer(ip), intent(out)       :: ws_neighbors (primal_mesh%nnode**2*dual_mesh%nnode**2)

    ! Local variables
    integer(ip)                     :: ineigh  
    integer(ip)                     :: ipoindm, jpoindm, kpoindm ! Dual   mesh point identifier
    integer(ip)                     :: ipoinpm, jpoinpm, kpoinpm ! Primal mesh point identifier
    integer(ip)                     :: ipoinpg ! Primal graph point identifier 
    integer(ip)                     :: p_ipoindm, p_jpoindm, p_kpoindm   
    integer(ip)                     :: p_ipoinpm, p_jpoinpm, p_kpoinpm   
    integer(ip)                     :: inods1,inods2,jnods1,jnods2,knods1,knods2  

    integer(ip)                     :: first_free_ja_pos
    integer(ip)                     :: first_free_pos
    integer(ip)                     :: num_neighbors

    assert(dual_mesh%nelty/=1)

    ! Initialize work space of filled-up 
    ! positions to zeros
    ws_position   = 0

    ! Initialize working space of neighbour 
    ! identifiers 
    ws_neighbors = 0

    first_free_ja_pos = 1

    do ipoinpg=1, primal_graph%nv
       first_free_pos = 1
       num_neighbors  = 0

       ! Projection vs Projection terms
       do p_ipoindm = dual_mesh%pnods(ipoinpg), dual_mesh%pnods(ipoinpg+1)-1
          ipoindm = dual_mesh%lnods(p_ipoindm) 

          !Nodes ipoinpm such that szmap(ipoinpm) = ipoindm
          inods1=szdual_ia(ipoindm)
          inods2=szdual_ia(ipoindm+1)-1
          do p_ipoinpm = inods1,inods2
             ipoinpm = szdual_ja(p_ipoinpm)
             
             !Elems jpoindm such that ipoinpm is in jpoindm
             do p_jpoindm = dual_mesh%pnods(ipoinpm), dual_mesh%pnods(ipoinpm+1)-1
                jpoindm = dual_mesh%lnods(p_jpoindm)

                !Nodes jpoinpm such that jpoinpm is in jpoindm
                if(primal_mesh%nelty==1) then
                   jnods1=(jpoindm-1)*primal_mesh%nnode+1
                   jnods2=jpoindm*primal_mesh%nnode
                else
                   jnods1=primal_mesh%pnods(jpoindm)
                   jnods2=primal_mesh%pnods(jpoindm+1)-1
                end if
                do p_jpoinpm = jnods1,jnods2
                   jpoinpm = primal_mesh%lnods(p_jpoinpm)

                   !Elems kpoindm = szmap(jpoinpm)
                   kpoindm = szmap(jpoinpm)

                   !Nodes kpoinpm such that kpoinpm is in kpoindm
                   if(primal_mesh%nelty==1) then
                      knods1=(kpoindm-1)*primal_mesh%nnode+1
                      knods2=kpoindm*primal_mesh%nnode
                   else
                      knods1=primal_mesh%pnods(kpoindm)
                      knods2=primal_mesh%pnods(kpoindm+1)-1
                   end if
                   do p_kpoinpm = knods1,knods2
                      kpoinpm = primal_mesh%lnods(p_kpoinpm)
                      ! write (*,*) ipoinpg, ipoindm, ipoinpm ! DBG: 
                      ! If ipoinpm not visited yet
                      if ( ws_position(kpoinpm) == 0 ) then
                         ws_position(kpoinpm) = first_free_pos
                         ws_neighbors(first_free_pos) = kpoinpm
                         first_free_pos = first_free_pos + 1
                         num_neighbors = num_neighbors + 1 
                      end if
                   end do
                end do
             end do
          end do
       end do
!!$
!!$       ! Cross terms
!!$       do p_ipoindm = dual_mesh%pnods(ipoinpg), dual_mesh%pnods(ipoinpg+1)-1
!!$          ipoindm = dual_mesh%lnods(p_ipoindm) 
!!$
!!$          ! Nodes ipoinpm such that ipoinpm is in ipoindm
!!$          if(primal_mesh%nelty==1) then
!!$             inods1=(ipoindm-1)*primal_mesh%nnode+1
!!$             inods2=ipoindm*primal_mesh%nnode
!!$          else
!!$             inods1=primal_mesh%pnods(ipoindm)
!!$             inods2=primal_mesh%pnods(ipoindm+1)-1
!!$          end if
!!$          do p_ipoinpm = inods1,inods2
!!$             ipoinpm = primal_mesh%lnods(p_ipoinpm)
!!$
!!$             !Elems jpoindm such that jpoindm=szmap(ipoinpm)
!!$             jpoindm = szmap(ipoinpm)
!!$
!!$             !Nodes jpoinpm such that jpoinpm is in jpoindm
!!$             if(primal_mesh%nelty==1) then
!!$                jnods1=(jpoindm-1)*primal_mesh%nnode+1
!!$                jnods2=jpoindm*primal_mesh%nnode
!!$             else
!!$                jnods1=primal_mesh%pnods(jpoindm)
!!$                jnods2=primal_mesh%pnods(jpoindm+1)-1
!!$             end if
!!$             do p_jpoinpm = jnods1,jnods2
!!$                jpoinpm = primal_mesh%lnods(p_jpoinpm)
!!$
!!$                if ( ws_position(jpoinpm) == 0 ) then
!!$                   ws_position(jpoinpm) = first_free_pos
!!$                   ws_neighbors(first_free_pos) = jpoinpm
!!$                   first_free_pos = first_free_pos + 1
!!$                   num_neighbors = num_neighbors + 1 
!!$                end if
!!$             end do
!!$          end do
!!$       end do
!!$       do p_ipoindm = dual_mesh%pnods(ipoinpg), dual_mesh%pnods(ipoinpg+1)-1
!!$          ipoindm = dual_mesh%lnods(p_ipoindm) 
!!$
!!$          !Nodes ipoinpm such that szmap(ipoinpm) = ipoindm
!!$          inods1=szdual_ia(ipoindm)
!!$          inods2=szdual_ia(ipoindm+1)-1
!!$          do p_ipoinpm = inods1,inods2
!!$             ipoinpm = szdual_ja(p_ipoinpm)
!!$
!!$             !Elems jpoindm such that ipoinpm is in jpoindm
!!$             do p_jpoindm = dual_mesh%pnods(ipoinpm), dual_mesh%pnods(ipoinpm+1)-1
!!$                jpoindm = dual_mesh%lnods(p_jpoindm)
!!$
!!$                !Nodes jpoinpm such that jpoinpm is in jpoindm
!!$                if(primal_mesh%nelty==1) then
!!$                   jnods1=(jpoindm-1)*primal_mesh%nnode+1
!!$                   jnods2=jpoindm*primal_mesh%nnode
!!$                else
!!$                   jnods1=primal_mesh%pnods(jpoindm)
!!$                   jnods2=primal_mesh%pnods(jpoindm+1)-1
!!$                end if
!!$                do p_jpoinpm = jnods1,jnods2
!!$                   jpoinpm = primal_mesh%lnods(p_jpoinpm)
!!$
!!$                   if ( ws_position(jpoinpm) == 0 ) then
!!$                      ws_position(jpoinpm) = first_free_pos
!!$                      ws_neighbors(first_free_pos) = jpoinpm
!!$                      first_free_pos = first_free_pos + 1
!!$                      num_neighbors = num_neighbors + 1 
!!$                   end if
!!$                end do
!!$             end do
!!$          end do
!!$       end do

       ! Extract neighbours while restoring working space to initial state
       do ineigh=1, num_neighbors
          primal_graph%ja(first_free_ja_pos) = ws_neighbors(ineigh)
          first_free_ja_pos = first_free_ja_pos + 1
          ws_position ( ws_neighbors(ineigh) ) = 0
          ws_neighbors( ineigh               ) = 0
       end do

       ! Order increasingly column identifiers of current row 
       ! using heap sort algorithm
       call psb_hsort(primal_graph%ja( primal_graph%ia(ipoinpg):primal_graph%ia(ipoinpg+1)-1 ))

    end do

  end subroutine list_primal_graph_sz_csr

  !============================================================================================
  subroutine  list_primal_graph_csr_scal ( ndof1, ndof2, primal_mesh, dual_mesh, primal_graph,  &
       &                                   ws_position, ws_neighbors )    
    use psb_sort_mod
    implicit none

    ! Parameters
    integer(ip)    , intent(in)    :: ndof1, ndof2
    type(fem_mesh) , intent(in)    :: primal_mesh, dual_mesh
    type(fem_graph), intent(inout) :: primal_graph
    integer(ip)    , intent(out)   :: ws_position  (primal_mesh%npoin)
    integer(ip)    , intent(out)   :: ws_neighbors (primal_mesh%nnode*dual_mesh%nnode)

    ! Local variables
    integer(ip)                    :: ineigh  
    integer(ip)                    :: ipoindm ! Dual   mesh point identifier
    integer(ip)                    :: ipoinpm ! Primal mesh point identifier
    integer(ip)                    :: ipoinpg ! Primal graph point identifier 
    integer(ip)                    :: p_ipoindm   
    integer(ip)                    :: p_ipoinpm   
    integer(ip)                    :: inods1,inods2   

    integer(ip)                    :: first_free_ja_pos
    integer(ip)                    :: first_free_pos
    integer(ip)                    :: num_neighbors, idof1, idof2

    assert(dual_mesh%nelty/=1)

    ! Initialize work space of filled-up 
    ! positions to zeros
    ws_position   = 0

    ! Initialize working space of neighbour 
    ! identifiers 
    ws_neighbors = 0

    first_free_ja_pos = 1

    do ipoinpg=1, primal_graph%nv, ndof1       
       first_free_pos = 1
       num_neighbors  = 0

       ! Traverse the dual nodes of the dual element number ipoinpg
       ! (i.e., primal elements around primal node number ipoinpg)
       do p_ipoindm = dual_mesh%pnods(scal2blk(ipoinpg,1,ndof1)), dual_mesh%pnods(scal2blk(ipoinpg,1,ndof1)+1)-1
          ipoindm = dual_mesh%lnods(p_ipoindm)

          ! Traverse the primal nodes of the primal element number ipoindm
          ! (i.e., dual elements around dual node number ipoindm)
          if(primal_mesh%nelty==1) then
             inods1=(ipoindm-1)*primal_mesh%nnode+1
             inods2=ipoindm*primal_mesh%nnode
          else
             inods1=primal_mesh%pnods(ipoindm)
             inods2=primal_mesh%pnods(ipoindm+1)-1
          end if
          do p_ipoinpm = inods1,inods2
             ipoinpm = primal_mesh%lnods(p_ipoinpm)
             ! If ipoinpm not visited yet
             if ( ws_position(ipoinpm) == 0 ) then
                ws_position(ipoinpm) = first_free_pos
                ws_neighbors(first_free_pos) = ipoinpm
                first_free_pos = first_free_pos + 1
                num_neighbors = num_neighbors + 1 
                ! If ipoinpm already visited
             end if
          end do
       end do

       ! Extract neighbours while restoring working space to initial state
       do ineigh=1, num_neighbors
          do idof2 = 1, ndof2
            primal_graph%ja(first_free_ja_pos) = blk2scal(ws_neighbors(ineigh),idof2,ndof2)
            first_free_ja_pos = first_free_ja_pos + 1
          end do

          ws_position ( ws_neighbors(ineigh) ) = 0
          ws_neighbors( ineigh               ) = 0
       end do


       ! Order increasingly column identifiers of current row 
       ! using heap sort algorithm
       ! write (*,*) 'A', primal_graph%ja( primal_graph%ia(ipoinpg):primal_graph%ia(ipoinpg+1)-1 )
       call psb_hsort(primal_graph%ja( primal_graph%ia(ipoinpg):primal_graph%ia(ipoinpg+1)-1 ))
       ! write (*,*) 'D', primal_graph%ja( primal_graph%ia(ipoinpg):primal_graph%ia(ipoinpg+1)-1 )


       do idof1 = 2, ndof1
           first_free_ja_pos = first_free_ja_pos + primal_graph%ia(ipoinpg+idof1) - primal_graph%ia(ipoinpg+idof1-1)
           primal_graph%ja( primal_graph%ia(ipoinpg+idof1-1):primal_graph%ia(ipoinpg+idof1)-1 ) = &
               & primal_graph%ja( primal_graph%ia(ipoinpg):primal_graph%ia(ipoinpg+1)-1 ) 
       end do

    end do

  end subroutine list_primal_graph_csr_scal

 !============================================================================================
  subroutine  list_primal_graph_extended_stencil_csr_scal ( ndof1, ndof2, primal_mesh, dual_mesh, primal_graph,  &
       &                                   ws_position, ws_neighbors )    
    use psb_sort_mod
    implicit none

    ! Parameters
    integer(ip)    , intent(in)    :: ndof1, ndof2
    type(fem_mesh) , intent(in)    :: primal_mesh, dual_mesh
    type(fem_graph), intent(inout) :: primal_graph
    integer(ip)    , intent(out)   :: ws_position  (primal_mesh%npoin)
    integer(ip)    , intent(out)   :: ws_neighbors (primal_mesh%nnode**2*dual_mesh%nnode**2)

    ! Local variables
    integer(ip)                    :: ineigh  
    integer(ip)                    :: ipoindm,jpoindm,kpoindm ! Dual   mesh point identifier
    integer(ip)                    :: ipoinpm,jpoinpm,kpoinpm ! Primal mesh point identifier
    integer(ip)                    :: ipoinpg ! Primal graph point identifier 
    integer(ip)                    :: p_ipoindm,p_jpoindm ,p_kpoindm  
    integer(ip)                    :: p_ipoinpm,p_jpoinpm,p_kpoinpm   
    integer(ip)                    :: inods1,inods2,jnods1,jnods2,knods1,knods2   

    integer(ip)                    :: first_free_ja_pos
    integer(ip)                    :: first_free_pos
    integer(ip)                    :: num_neighbors, idof1, idof2

    assert(dual_mesh%nelty/=1)

    ! Initialize work space of filled-up 
    ! positions to zeros
    ws_position   = 0

    ! Initialize working space of neighbour 
    ! identifiers 
    ws_neighbors = 0

    first_free_ja_pos = 1

    do ipoinpg=1, primal_graph%nv, ndof1       
       first_free_pos = 1
       num_neighbors  = 0

       ! Traverse the dual nodes of the dual element number ipoinpg
       ! (i.e., primal elements around primal node number ipoinpg)
       do p_ipoindm = dual_mesh%pnods(scal2blk(ipoinpg,1,ndof1)), dual_mesh%pnods(scal2blk(ipoinpg,1,ndof1)+1)-1
          ipoindm = dual_mesh%lnods(p_ipoindm)

          ! Traverse the primal nodes of the primal element number ipoindm
          ! (i.e., dual elements around dual node number ipoindm)
          if(primal_mesh%nelty==1) then
             inods1=(ipoindm-1)*primal_mesh%nnode+1
             inods2=ipoindm*primal_mesh%nnode
          else
             inods1=primal_mesh%pnods(ipoindm)
             inods2=primal_mesh%pnods(ipoindm+1)-1
          end if
          do p_ipoinpm = inods1,inods2
             ipoinpm = primal_mesh%lnods(p_ipoinpm)

             do p_jpoindm = dual_mesh%pnods(scal2blk(ipoinpm,1,ndof1)), &
                  dual_mesh%pnods(scal2blk(ipoinpm,1,ndof1)+1)-1
                jpoindm = dual_mesh%lnods(p_jpoindm)
                
                ! Traverse the primal nodes of the primal element number ipoindm
                ! (i.e., dual elements around dual node number ipoindm)
                if(primal_mesh%nelty==1) then
                   jnods1=(jpoindm-1)*primal_mesh%nnode+1
                   jnods2=jpoindm*primal_mesh%nnode
                else
                   jnods1=primal_mesh%pnods(jpoindm)
                   jnods2=primal_mesh%pnods(jpoindm+1)-1
                end if
                do p_jpoinpm = jnods1,jnods2
                   jpoinpm = primal_mesh%lnods(p_jpoinpm)
                   do p_kpoindm = dual_mesh%pnods(jpoinpm), dual_mesh%pnods(jpoinpm+1)-1
                     kpoindm = dual_mesh%lnods(p_kpoindm)

                     ! Traverse the primal nodes of the primal element number ipoindm
                     ! (i.e., dual elements around dual node number ipoindm)
                     if(primal_mesh%nelty==1) then
                        knods1=(kpoindm-1)*primal_mesh%nnode+1
                        knods2=kpoindm*primal_mesh%nnode
                     else
                        knods1=primal_mesh%pnods(kpoindm)
                        knods2=primal_mesh%pnods(kpoindm+1)-1
                     end if
                     do p_kpoinpm = knods1,knods2
                        kpoinpm = primal_mesh%lnods(p_kpoinpm)
                        ! If ipoinpm not visited yet
                        if ( ws_position(kpoinpm) == 0 ) then
                           ws_position(kpoinpm) = first_free_pos
                           ws_neighbors(first_free_pos) = kpoinpm
                           first_free_pos = first_free_pos + 1
                           num_neighbors = num_neighbors + 1 
                           ! If ipoinpm already visited
                        end if
                     end do
                  end do
                end do
             end do
          end do
       end do

       ! Extract neighbours while restoring working space to initial state
       do ineigh=1, num_neighbors
          do idof2 = 1, ndof2
            primal_graph%ja(first_free_ja_pos) = blk2scal(ws_neighbors(ineigh),idof2,ndof2)
            first_free_ja_pos = first_free_ja_pos + 1
          end do

          ws_position ( ws_neighbors(ineigh) ) = 0
          ws_neighbors( ineigh               ) = 0
       end do


       ! Order increasingly column identifiers of current row 
       ! using heap sort algorithm
       ! write (*,*) 'A', primal_graph%ja( primal_graph%ia(ipoinpg):primal_graph%ia(ipoinpg+1)-1 )
       call psb_hsort(primal_graph%ja( primal_graph%ia(ipoinpg):primal_graph%ia(ipoinpg+1)-1 ))
       ! write (*,*) 'D', primal_graph%ja( primal_graph%ia(ipoinpg):primal_graph%ia(ipoinpg+1)-1 )


       do idof1 = 2, ndof1
           first_free_ja_pos = first_free_ja_pos + primal_graph%ia(ipoinpg+idof1) - primal_graph%ia(ipoinpg+idof1-1)
           primal_graph%ja( primal_graph%ia(ipoinpg+idof1-1):primal_graph%ia(ipoinpg+idof1)-1 ) = &
               & primal_graph%ja( primal_graph%ia(ipoinpg):primal_graph%ia(ipoinpg+1)-1 ) 
       end do

    end do

  end subroutine list_primal_graph_extended_stencil_csr_scal

  !============================================================================================
  subroutine  list_primal_graph_sz_csr_scal ( ndof1, ndof2, primal_mesh, dual_mesh, primal_graph,szmap,  &
       &                                   szdual_ia,szdual_ja,ws_position, ws_neighbors )    
    use psb_sort_mod
    implicit none

    ! Parameters
    integer(ip)    , intent(in)    :: ndof1, ndof2
    type(fem_mesh) , intent(in)    :: primal_mesh, dual_mesh
    integer(ip)    , intent(in)    :: szmap(:),szdual_ia(:),szdual_ja(:)
    type(fem_graph), intent(inout) :: primal_graph
    integer(ip)    , intent(out)   :: ws_position  (primal_mesh%npoin)
    integer(ip)    , intent(out)   :: ws_neighbors (primal_mesh%nnode**2*dual_mesh%nnode**2)

    ! Local variables
    integer(ip)                    :: ineigh  
    integer(ip)                    :: ipoindm,jpoindm,kpoindm ! Dual   mesh point identifier
    integer(ip)                    :: ipoinpm,jpoinpm,kpoinpm ! Primal mesh point identifier
    integer(ip)                    :: ipoinpg ! Primal graph point identifier 
    integer(ip)                    :: p_ipoindm,p_jpoindm ,p_kpoindm  
    integer(ip)                    :: p_ipoinpm,p_jpoinpm,p_kpoinpm   
    integer(ip)                    :: inods1,inods2,jnods1,jnods2,knods1,knods2   

    integer(ip)                    :: first_free_ja_pos
    integer(ip)                    :: first_free_pos
    integer(ip)                    :: num_neighbors, idof1, idof2

    assert(dual_mesh%nelty/=1)

    ! Initialize work space of filled-up 
    ! positions to zeros
    ws_position   = 0

    ! Initialize working space of neighbour 
    ! identifiers 
    ws_neighbors = 0

    first_free_ja_pos = 1

    do ipoinpg=1, primal_graph%nv, ndof1       
       first_free_pos = 1
       num_neighbors  = 0

       ! Projection vs Projection terms
       do p_ipoindm = dual_mesh%pnods(scal2blk(ipoinpg,1,ndof1)), dual_mesh%pnods(scal2blk(ipoinpg,1,ndof1)+1)-1
          ipoindm = dual_mesh%lnods(p_ipoindm) 

          !Nodes ipoinpm such that szmap(ipoinpm) = ipoindm
          inods1=szdual_ia(ipoindm)
          inods2=szdual_ia(ipoindm+1)-1
          do p_ipoinpm = inods1,inods2
             ipoinpm = szdual_ja(p_ipoinpm)
             
             !Elems jpoindm such that ipoinpm is in jpoindm
             do p_jpoindm = dual_mesh%pnods(scal2blk(ipoinpm,1,ndof1)), dual_mesh%pnods(scal2blk(ipoinpm,1,ndof1)+1)-1
                jpoindm = dual_mesh%lnods(p_jpoindm)

                !Nodes jpoinpm such that jpoinpm is in jpoindm
                if(primal_mesh%nelty==1) then
                   jnods1=(jpoindm-1)*primal_mesh%nnode+1
                   jnods2=jpoindm*primal_mesh%nnode
                else
                   jnods1=primal_mesh%pnods(jpoindm)
                   jnods2=primal_mesh%pnods(jpoindm+1)-1
                end if
                do p_jpoinpm = jnods1,jnods2
                   jpoinpm = primal_mesh%lnods(p_jpoinpm)

                   !Elems kpoindm = szmap(jpoinpm)
                   kpoindm = szmap(jpoinpm)

                   !Nodes kpoinpm such that kpoinpm is in kpoindm
                   if(primal_mesh%nelty==1) then
                      knods1=(kpoindm-1)*primal_mesh%nnode+1
                      knods2=kpoindm*primal_mesh%nnode
                   else
                      knods1=primal_mesh%pnods(kpoindm)
                      knods2=primal_mesh%pnods(kpoindm+1)-1
                   end if
                   do p_kpoinpm = knods1,knods2
                      kpoinpm = primal_mesh%lnods(p_kpoinpm)
                      ! write (*,*) ipoinpg, ipoindm, ipoinpm ! DBG: 
                      ! If ipoinpm not visited yet
                      if ( ws_position(kpoinpm) == 0 ) then
                         ws_position(kpoinpm) = first_free_pos
                         ws_neighbors(first_free_pos) = kpoinpm
                         first_free_pos = first_free_pos + 1
                         num_neighbors = num_neighbors + 1 
                      end if
                   end do
                end do
             end do
          end do
       end do
!!$
!!$        ! Cross terms
!!$       do p_ipoindm = dual_mesh%pnods(scal2blk(ipoinpg,1,ndof1)), dual_mesh%pnods(scal2blk(ipoinpg,1,ndof1)+1)-1
!!$          ipoindm = dual_mesh%lnods(p_ipoindm) 
!!$
!!$          ! Nodes ipoinpm such that ipoinpm is in ipoindm
!!$          if(primal_mesh%nelty==1) then
!!$             inods1=(ipoindm-1)*primal_mesh%nnode+1
!!$             inods2=ipoindm*primal_mesh%nnode
!!$          else
!!$             inods1=primal_mesh%pnods(ipoindm)
!!$             inods2=primal_mesh%pnods(ipoindm+1)-1
!!$          end if
!!$          do p_ipoinpm = inods1,inods2
!!$             ipoinpm = primal_mesh%lnods(p_ipoinpm)
!!$
!!$             !Elems jpoindm such that jpoindm=szmap(ipoinpm)
!!$             jpoindm = szmap(ipoinpm)
!!$
!!$             !Nodes jpoinpm such that jpoinpm is in jpoindm
!!$             if(primal_mesh%nelty==1) then
!!$                jnods1=(jpoindm-1)*primal_mesh%nnode+1
!!$                jnods2=jpoindm*primal_mesh%nnode
!!$             else
!!$                jnods1=primal_mesh%pnods(jpoindm)
!!$                jnods2=primal_mesh%pnods(jpoindm+1)-1
!!$             end if
!!$             do p_jpoinpm = jnods1,jnods2
!!$                jpoinpm = primal_mesh%lnods(p_jpoinpm)
!!$
!!$                if ( ws_position(jpoinpm) == 0 ) then
!!$                   write(*,*) "jpoin[3]", jpoinpm
!!$                   ws_position(jpoinpm) = first_free_pos
!!$                   ws_neighbors(first_free_pos) = jpoinpm
!!$                   first_free_pos = first_free_pos + 1
!!$                   num_neighbors = num_neighbors + 1 
!!$                end if
!!$             end do
!!$          end do
!!$       end do
!!$       do p_ipoindm = dual_mesh%pnods(scal2blk(ipoinpg,1,ndof1)), dual_mesh%pnods(scal2blk(ipoinpg,1,ndof1)+1)-1
!!$          ipoindm = dual_mesh%lnods(p_ipoindm) 
!!$
!!$          !Nodes ipoinpm such that szmap(ipoinpm) = ipoindm
!!$          inods1=szdual_ia(ipoindm)
!!$          inods2=szdual_ia(ipoindm+1)-1
!!$          do p_ipoinpm = inods1,inods2
!!$             ipoinpm = szdual_ja(p_ipoinpm)
!!$
!!$             !Elems jpoindm such that ipoinpm is in jpoindm
!!$             do p_jpoindm =dual_mesh%pnods(scal2blk(ipoinpm,1,ndof1)),dual_mesh%pnods(scal2blk(ipoinpm,1,ndof1)+1)-1
!!$                jpoindm = dual_mesh%lnods(p_jpoindm)
!!$
!!$                !Nodes jpoinpm such that jpoinpm is in jpoindm
!!$                if(primal_mesh%nelty==1) then
!!$                   jnods1=(jpoindm-1)*primal_mesh%nnode+1
!!$                   jnods2=jpoindm*primal_mesh%nnode
!!$                else
!!$                   jnods1=primal_mesh%pnods(jpoindm)
!!$                   jnods2=primal_mesh%pnods(jpoindm+1)-1
!!$                end if
!!$                do p_jpoinpm = jnods1,jnods2
!!$                   jpoinpm = primal_mesh%lnods(p_jpoinpm)
!!$
!!$                   if ( ws_position(jpoinpm) == 0 ) then
!!$                      write(*,*) "jpoin[2]", jpoinpm
!!$                      ws_position(jpoinpm) = first_free_pos
!!$                      ws_neighbors(first_free_pos) = jpoinpm
!!$                      first_free_pos = first_free_pos + 1
!!$                      num_neighbors = num_neighbors + 1 
!!$                   end if
!!$                end do
!!$             end do
!!$          end do
!!$       end do

       ! Extract neighbours while restoring working space to initial state
       do ineigh=1, num_neighbors
          do idof2 = 1, ndof2
             primal_graph%ja(first_free_ja_pos) = blk2scal(ws_neighbors(ineigh),idof2,ndof2)
             first_free_ja_pos = first_free_ja_pos + 1
          end do

          ws_position ( ws_neighbors(ineigh) ) = 0
          ws_neighbors( ineigh               ) = 0
       end do


       ! Order increasingly column identifiers of current row 
       ! using heap sort algorithm
       ! write (*,*) 'A', primal_graph%ja( primal_graph%ia(ipoinpg):primal_graph%ia(ipoinpg+1)-1 )
       call psb_hsort(primal_graph%ja( primal_graph%ia(ipoinpg):primal_graph%ia(ipoinpg+1)-1 ))
       ! write (*,*) 'D', primal_graph%ja( primal_graph%ia(ipoinpg):primal_graph%ia(ipoinpg+1)-1 )


       do idof1 = 2, ndof1
          first_free_ja_pos = first_free_ja_pos + primal_graph%ia(ipoinpg+idof1) - primal_graph%ia(ipoinpg+idof1-1)
          primal_graph%ja( primal_graph%ia(ipoinpg+idof1-1):primal_graph%ia(ipoinpg+idof1)-1 ) = &
               & primal_graph%ja( primal_graph%ia(ipoinpg):primal_graph%ia(ipoinpg+1)-1 ) 
       end do

    end do

  end subroutine list_primal_graph_sz_csr_scal

  !============================================================================================
  subroutine  list_primal_graph_csr_symm ( primal_mesh, dual_mesh, primal_graph,  &
       &                                   ws_position, ws_neighbors )    
    use psb_sort_mod
    implicit none

    ! Parameters
    type(fem_mesh) , intent(in)    :: primal_mesh, dual_mesh
    type(fem_graph), intent(inout) :: primal_graph
    integer(ip), intent(out)       :: ws_position  (primal_mesh%npoin)
    integer(ip), intent(out)       :: ws_neighbors (primal_mesh%nnode*dual_mesh%nnode)

    ! Local variables
    integer(ip)                     :: ineigh  
    integer(ip)                     :: ipoindm ! Dual   mesh point identifier
    integer(ip)                     :: ipoinpm ! Primal mesh point identifier
    integer(ip)                     :: ipoinpg ! Primal graph point identifier 
    integer(ip)                     :: p_ipoindm   
    integer(ip)                     :: p_ipoinpm   
    integer(ip)                     :: inods1,inods2   

    integer(ip)                     :: first_free_ja_pos
    integer(ip)                     :: first_free_pos
    integer(ip)                     :: num_neighbors

    assert(dual_mesh%nelty/=1)

    ! Initialize work space of filled-up 
    ! positions to zeros
    ws_position   = 0

    ! Initialize working space of neighbour 
    ! identifiers 
    ws_neighbors = 0

    first_free_ja_pos = 1

    do ipoinpg=1, primal_graph%nv       
       first_free_pos = 1
       num_neighbors  = 0

       ! Traverse the dual nodes of the dual element number ipoinpg
       ! (i.e., primal elements around primal node number ipoinpg)
       do p_ipoindm = dual_mesh%pnods(ipoinpg), dual_mesh%pnods(ipoinpg+1)-1
          ipoindm = dual_mesh%lnods(p_ipoindm)

          ! Traverse the primal nodes of the primal element number ipoindm
          ! (i.e., dual elements around dual node number ipoindm)
          if(primal_mesh%nelty==1) then
             inods1=(ipoindm-1)*primal_mesh%nnode+1
             inods2=ipoindm*primal_mesh%nnode
          else
             inods1=primal_mesh%pnods(ipoindm)
             inods2=primal_mesh%pnods(ipoindm+1)-1
          end if
          do p_ipoinpm = inods1,inods2
             ipoinpm = primal_mesh%lnods(p_ipoinpm)
             ! Only list edges (i,j) s.t., i <= j  
             if ( ipoinpg <= ipoinpm ) then
                ! If ipoinpm not visited yet
                if ( ws_position(ipoinpm) == 0 ) then
                   ws_position(ipoinpm) = first_free_pos
                   ws_neighbors(first_free_pos) = ipoinpm
                   first_free_pos = first_free_pos + 1
                   num_neighbors = num_neighbors + 1 
                   ! If ipoinpm already visited
                end if
             end if
          end do
       end do

       ! Extract neighbours while restoring working space to initial state
       do ineigh=1, num_neighbors
          primal_graph%ja(first_free_ja_pos) = ws_neighbors(ineigh)
          first_free_ja_pos = first_free_ja_pos + 1
          ws_position ( ws_neighbors(ineigh) ) = 0
          ws_neighbors( ineigh               ) = 0
       end do

       ! Order increasingly column identifiers of current row 
       ! using heap sort algorithm
       call psb_hsort(primal_graph%ja( primal_graph%ia(ipoinpg):primal_graph%ia(ipoinpg+1)-1 ))

    end do

  end subroutine list_primal_graph_csr_symm

  !============================================================================================
  subroutine  list_primal_graph_csr_symm_scal ( ndof1, ndof2, primal_mesh, dual_mesh, primal_graph,  &
       &                                        ws_position, ws_neighbors )    
    use psb_sort_mod
    implicit none

    ! Parameters
    integer(ip)    , intent(in)    :: ndof1, ndof2
    type(fem_mesh) , intent(in)    :: primal_mesh, dual_mesh
    type(fem_graph), intent(inout) :: primal_graph
    integer(ip)    , intent(out)   :: ws_position  (primal_mesh%npoin)
    integer(ip)    , intent(out)   :: ws_neighbors (primal_mesh%nnode*dual_mesh%nnode)

    ! Local variables
    integer(ip)                    :: ineigh  
    integer(ip)                    :: ipoindm ! Dual   mesh point identifier
    integer(ip)                    :: ipoinpm ! Primal mesh point identifier
    integer(ip)                    :: ipoinpg ! Primal graph point identifier 
    integer(ip)                    :: p_ipoindm   
    integer(ip)                    :: p_ipoinpm   
    integer(ip)                    :: inods1,inods2   

    integer(ip)                    :: first_free_ja_pos
    integer(ip)                    :: first_free_pos
    integer(ip)                    :: num_neighbors, idof1, idof2

    assert(dual_mesh%nelty/=1)

    ! Initialize work space of filled-up 
    ! positions to zeros
    ws_position   = 0

    ! Initialize working space of neighbour 
    ! identifiers 
    ws_neighbors = 0

    first_free_ja_pos = 1

    do ipoinpg=1, primal_graph%nv, ndof1      
       first_free_pos = 1
       num_neighbors  = 0

       ! Traverse the dual nodes of the dual element number ipoinpg
       ! (i.e., primal elements around primal node number ipoinpg)
       do p_ipoindm = dual_mesh%pnods(scal2blk(ipoinpg,1,ndof1)), dual_mesh%pnods(scal2blk(ipoinpg,1,ndof1)+1)-1
          ipoindm = dual_mesh%lnods(p_ipoindm)

          ! Traverse the primal nodes of the primal element number ipoindm
          ! (i.e., dual elements around dual node number ipoindm)
          if(primal_mesh%nelty==1) then
             inods1=(ipoindm-1)*primal_mesh%nnode+1
             inods2=ipoindm*primal_mesh%nnode
          else
             inods1=primal_mesh%pnods(ipoindm)
             inods2=primal_mesh%pnods(ipoindm+1)-1
          end if
          do p_ipoinpm = inods1,inods2
             ipoinpm = primal_mesh%lnods(p_ipoinpm)
             ! Only list edges (i,j) s.t., i <= j  
             if ( scal2blk(ipoinpg,1,ndof1) <= ipoinpm ) then
                ! If ipoinpm not visited yet
                if ( ws_position(ipoinpm) == 0 ) then
                   ws_position(ipoinpm) = first_free_pos
                   ws_neighbors(first_free_pos) = ipoinpm
                   first_free_pos = first_free_pos + 1
                   num_neighbors = num_neighbors + 1 
                   ! If ipoinpm already visited
                end if
             end if
          end do
       end do

       ! Extract neighbours while restoring working space to initial state
       do ineigh=1, num_neighbors
          do idof2 = 1, ndof2
            primal_graph%ja(first_free_ja_pos) = blk2scal(ws_neighbors(ineigh),idof2,ndof2)
            first_free_ja_pos = first_free_ja_pos + 1
          end do

          ws_position ( ws_neighbors(ineigh) ) = 0
          ws_neighbors( ineigh               ) = 0
       end do

       ! Order increasingly column identifiers of current row 
       ! using heap sort algorithm
       ! write (*,*) 'A', primal_graph%ja( primal_graph%ia(ipoinpg):primal_graph%ia(ipoinpg+1)-1 ) ! DBG:
       call psb_hsort(primal_graph%ja( primal_graph%ia(ipoinpg):primal_graph%ia(ipoinpg+1)-1 ))
       ! write (*,*) 'D', primal_graph%ja( primal_graph%ia(ipoinpg):primal_graph%ia(ipoinpg+1)-1 ) ! DBG:

 
       do idof1 = 2, ndof1
           first_free_ja_pos = first_free_ja_pos + primal_graph%ia(ipoinpg+idof1) - primal_graph%ia(ipoinpg+idof1-1)
           ! write(*,*) 'ZZZ', primal_graph%ia(ipoinpg+idof1) - primal_graph%ia(ipoinpg+idof1-1), size(primal_graph%ja), first_free_ja_pos
           primal_graph%ja( primal_graph%ia(ipoinpg+idof1-1):primal_graph%ia(ipoinpg+idof1)-1 ) = &
               & primal_graph%ja( (idof1-1+primal_graph%ia(ipoinpg)):primal_graph%ia(ipoinpg+1)-1 ) 
       end do

    end do

  end subroutine list_primal_graph_csr_symm_scal

  !============================================================================================
  subroutine  list_primal_graph_css ( primal_mesh, dual_mesh, primal_graph,  &
       &                                ws_position, ws_neighbors )    
    ! use psb_sort_mod
    implicit none

    ! Parameters
    type(fem_mesh) , intent(in)    :: primal_mesh, dual_mesh
    type(fem_graph), intent(inout) :: primal_graph
    integer(ip), intent(out)       :: ws_position  (primal_mesh%npoin)
    integer(ip), intent(out)       :: ws_neighbors (primal_mesh%nnode*dual_mesh%nnode)

    ! Local variables
    integer(ip)                     :: ineigh  
    integer(ip)                     :: ipoindm ! Dual   mesh point identifier
    integer(ip)                     :: ipoinpm ! Primal mesh point identifier
    integer(ip)                     :: ipoinpg ! Primal graph point identifier 
    integer(ip)                     :: p_ipoindm   
    integer(ip)                     :: p_ipoinpm   
    integer(ip)                     :: inods1,inods2   

    integer(ip)                     :: first_free_ja_low, first_free_ja_upp
    integer(ip)                     :: first_free_pos
    integer(ip)                     :: num_neighbors_low, num_neighbors_upp

    assert(dual_mesh%nelty/=1)

    ! Initialize work space of filled-up 
    ! positions to zeros
    ws_position   = 0

    ! Initialize working space of neighbour 
    ! identifiers 
    ws_neighbors = 0

    do ipoinpg=1, primal_graph%nv       
       first_free_pos = 1
       num_neighbors_low  = 0
       num_neighbors_upp  = 0

       ! Traverse the dual nodes of the dual element number ipoinpg
       ! (i.e., primal elements around primal node number ipoinpg)
       do p_ipoindm = dual_mesh%pnods(ipoinpg), dual_mesh%pnods(ipoinpg+1)-1
          ipoindm = dual_mesh%lnods(p_ipoindm)

          ! Traverse the primal nodes of the primal element number ipoindm
          ! (i.e., dual elements around dual node number ipoindm)
          if(primal_mesh%nelty==1) then
             inods1=(ipoindm-1)*primal_mesh%nnode+1
             inods2=ipoindm*primal_mesh%nnode
          else
             inods1=primal_mesh%pnods(ipoindm)
             inods2=primal_mesh%pnods(ipoindm+1)-1
          end if
          do p_ipoinpm = inods1,inods2
             ipoinpm = primal_mesh%lnods(p_ipoinpm)
             if ( (ipoinpg > ipoinpm) ) then ! Lower triangle entries
                ! If ipoinpm not visited yet
                if ( ws_position(ipoinpm) == 0 ) then
                   ws_position(ipoinpm) = first_free_pos
                   ws_neighbors(first_free_pos) = ipoinpm
                   first_free_pos = first_free_pos + 1
                   num_neighbors_low = num_neighbors_low + 1 
                   ! If ipoinpm already visited
                end if
             else
                if ( (ipoinpg < ipoinpm) ) then ! Upper triangle entries
                   ! If ipoinpm not visited yet
                   if ( ws_position(ipoinpm) == 0 ) then
                      ws_position(ipoinpm) = first_free_pos
                      ws_neighbors(first_free_pos) = ipoinpm
                      first_free_pos = first_free_pos + 1
                      num_neighbors_upp = num_neighbors_upp + 1 
                      ! If ipoinpm already visited
                   end if
                end if
             end if
          end do
       end do

       first_free_ja_low = primal_graph%ia (ipoinpg)
       first_free_ja_upp = primal_graph%ia (ipoinpg+1) &
            &  - (primal_graph%is (ipoinpg+1)- primal_graph%is (ipoinpg))

       ! write (*,*) first_free_ja_low, first_free_ja_upp ! DBG:

       ! Extract neighbours while restoring working space to initial state
       do ineigh=1, num_neighbors_low + num_neighbors_upp
          if ( ipoinpg > ws_neighbors(ineigh)  ) then ! Lower triangle entries
             primal_graph%ja(first_free_ja_low) = ws_neighbors(ineigh)
             first_free_ja_low = first_free_ja_low + 1 
             ws_position ( ws_neighbors(ineigh) ) = 0
             ws_neighbors( ineigh               ) = 0
          else
             if ( ipoinpg < ws_neighbors(ineigh) ) then ! Upper triangle entries
                primal_graph%ja(first_free_ja_upp) = ws_neighbors(ineigh)
                first_free_ja_upp = first_free_ja_upp + 1
                ws_position ( ws_neighbors(ineigh) ) = 0
                ws_neighbors( ineigh               ) = 0
             end if
          end if
       end do
       ! write (*,*) 'XXX', ipoinpg, &
       ! & primal_graph%ja (primal_graph%ia (ipoinpg):primal_graph%ia (ipoinpg+1)-1)

    end do

  end subroutine list_primal_graph_css

  !============================================================================================
  subroutine  list_primal_graph_css_scal ( ndof1, ndof2, primal_mesh, dual_mesh, primal_graph,  &
       &                                ws_position, ws_neighbors )    
    ! use psb_sort_mod
    implicit none

    ! Parameters
    integer(ip)    , intent(in)    :: ndof1, ndof2
    type(fem_mesh) , intent(in)    :: primal_mesh, dual_mesh
    type(fem_graph), intent(inout) :: primal_graph
    integer(ip)    , intent(out)   :: ws_position  (primal_mesh%npoin)
    integer(ip)    , intent(out)   :: ws_neighbors (primal_mesh%nnode*dual_mesh%nnode)

    ! Local variables
    integer(ip)                    :: ineigh  
    integer(ip)                    :: ipoindm ! Dual   mesh point identifier
    integer(ip)                    :: ipoinpm ! Primal mesh point identifier
    integer(ip)                    :: ipoinpg ! Primal graph point identifier 
    integer(ip)                    :: p_ipoindm   
    integer(ip)                    :: p_ipoinpm   
    integer(ip)                    :: inods1,inods2 
  
    integer(ip)                    :: first_free_ja_low, first_free_ja_upp
    integer(ip)                    :: first_free_pos
    integer(ip)                    :: num_neighbors_low, num_neighbors_upp

    assert(dual_mesh%nelty/=1)

    write (0,*) 'Error: the body of list_primal_graph_css_scal in mesh_graph.f90 still to be written'
    write (0,*) 'Error: the body of list_primal_graph_css in mesh_graph.f90 is a good point to start'
    write (0,*) 'Error: volunteers are welcome !!!'
    
    stop

  end subroutine list_primal_graph_css_scal


  !================================================================================================
  subroutine mesh_to_graph_part (primal_mesh, dual_mesh, min_freq_neig, primal_graph)
    !---------------------------------------------------------------------------
    ! This routine generates a graph of a given primal mesh (this graph 
    ! is a list of primal points around primal points)
    !
    ! The routine exploits the duality among primal/dual meshes. In order to 
    ! work properly pnods has to be an array of pointers to address lnods
    ! (!!! this does not currently happen with FE meshes with eltyp==1 !!!)
    !
    ! min_freq_neig is an input integer which determines which are
    ! the neighbours of a particular vertex of the graph. In particular, 
    ! neighbours are those vertices which are listed at least min_freq_neig times 
    ! when traversing the list of primal elements around a primal point
    !
    ! IMPORTANT NOTE:  This routine DOES NOT generate self-edges for both primal 
    !                  and dual graphs. However, mesh_to_graph_dual
    !                  generates self-edges. METIS graphs do NOT have self-edges.
    !---------------------------------------------------------------------------
    implicit none

    ! Parameters
    type(fem_mesh) , intent(in)  :: primal_mesh, dual_mesh
    integer(ip), intent(in)      :: min_freq_neig
    type(fem_graph), intent(out) :: primal_graph

    ! Local variables
    integer(ip), allocatable :: iwork(:)            ! Integer ip working array
    integer(ip)              :: pwork(4)            ! Pointers to work space

    ! Allocate space for ia on the primal graph
    primal_graph%type = part
    primal_graph%nv   = primal_mesh%npoin
    call memalloc (primal_graph%nv+1, primal_graph%ia, __FILE__,__LINE__)

    ! Allocate working space for count_primal_graph and list_primal_graph routines
    ! (TOTAL WS SIZE = primal mesh npoin + 2*maximum number of neighbours of any primal
    ! graph node)
    pwork(1) = 1
    pwork(2) = pwork(1) + primal_mesh%npoin
    pwork(3) = pwork(2) + dual_mesh%nnode*primal_mesh%nnode
    pwork(4) = pwork(3) + dual_mesh%nnode*primal_mesh%nnode
    call memalloc (pwork(4), iwork, __FILE__,__LINE__)

    ! write(*,*) 'begin count_primal_graph' ! DBG:
    ! Calculate ia array
    call count_primal_graph_part ( primal_mesh, dual_mesh, min_freq_neig, primal_graph,  &
         &                           iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)), &
         &                           iwork(pwork(3):pwork(4)) ) 
    ! write(*,*) 'end count_primal_graph'   ! DBG:

    ! Allocate space for ja on the primal graph 
    call memalloc (primal_graph%ia(primal_graph%nv+1)-1, primal_graph%ja,           __FILE__,__LINE__)

    ! write(*,*) 'begin list_primal_graph' ! DBG:
    ! Calculate ja array
    call list_primal_graph_part  (primal_mesh, dual_mesh, min_freq_neig, primal_graph, & 
         &                        iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)),  &
         &                        iwork(pwork(3):pwork(4)) )
    ! write(*,*) 'end list_primal_graph'   ! DBG: 

    ! Deallocate working space
    call memfree (iwork,__FILE__,__LINE__)
    return
  end subroutine mesh_to_graph_part

  !============================================================================================
  subroutine  count_primal_graph_part ( primal_mesh, dual_mesh, min_freq_neig, primal_graph,  &
       &                                  ws_position, ws_freq, ws_neighbors)
    implicit none

    ! Parameters
    type(fem_mesh) , intent(in)     :: primal_mesh, dual_mesh
    integer(ip), intent(in)         :: min_freq_neig
    type(fem_graph), intent(inout)  :: primal_graph
    integer(ip), intent(out)        :: ws_position  (primal_mesh%npoin)
    integer(ip), intent(out)        :: ws_freq      (primal_mesh%nnode*dual_mesh%nnode)
    integer(ip), intent(out)        :: ws_neighbors (primal_mesh%nnode*dual_mesh%nnode)

    ! Local variables
    integer(ip)                 :: ineigh  
    integer(ip)                 :: ipoindm ! Dual   mesh point identifier
    integer(ip)                 :: ipoinpm ! Primal mesh point identifier
    integer(ip)                 :: ipoinpg ! Primal graph point identifier 
    integer(ip)                 :: p_ipoindm   
    integer(ip)                 :: p_ipoinpm   
    integer(ip)                 :: inods1 , inods2   
    integer(ip)                 :: inods1d, inods2d

    integer(ip)                 :: first_free_pos
    integer(ip)                 :: num_neighbors

    
    ! I do not agree with the following
    ! assertion. dual_mesh can actually be
    ! a primal finite element mesh when this
    ! routine has to compute dual graphs (i.e., elements 
    ! around elements). This last  situation was not taken 
    ! into account in the code below.
    ! I have modified the code below so that it
    ! is also able to compute dual graphs. 
    
    ! !!!! IMPORTANT NOTE: INCORRECT ASSERTION !!!
    ! assert(dual_mesh%nelty/=1)

    ! CORRECT ASSERTION: Either primal_mesh or dual_mesh has
    !                    to be a dual_mesh
    assert ( primal_mesh%nelty /= 1 .or. dual_mesh%nelty /= 1 ) 

    ! Initialize work space of filled-up 
    ! positions
    ws_position   = 0

    ! Initialize work space of freqs. 
    ! associated to neighbours
    ws_freq       = 0

    ! Initialize work space of neighbour 
    ! identifiers 
    ws_neighbors  = 0

    !primal_graph%ia    = -1                                                ! DBG:
    !write (*,*) size(primal_graph%ia)                                      ! DBG: 
    !write (*,*) size(primal_mesh%pnods)                                    ! DBG: 
    !write (*,*) size(primal_mesh%lnods)                                    ! DBG: 
    !write (*,*) size(dual_mesh%pnods)                                      ! DBG: 
    !write (*,*) size(dual_mesh%lnods)                                      ! DBG: 
    !write (*,*) size(ws_position)                                          ! DBG: 
    !write (*,*) size(ws_freq)                                              ! DBG: 
    !write (*,*) size(ws_neighbors)                                         ! DBG: 
    !write (*, '(10i10)')    primal_mesh%pnods  (1:(primal_mesh%nelem+1))   ! DBG: 
    !write (*, '(10i10)')    primal_graph%ia    (1:(primal_graph%nv+1))     ! DBG:  
    !primal_graph%ia    = -1 
    primal_graph%ia(1) = 1

    do ipoinpg=1, primal_mesh%npoin
       primal_graph%ia(ipoinpg+1) = primal_graph%ia(ipoinpg)
       first_free_pos = 1
       num_neighbors  = 0
       ! Traverse the dual nodes of the dual element number ipoinpg
       ! (i.e., primal elements around primal node number ipoinpg)
       if(dual_mesh%nelty==1) then
          inods1d=(ipoinpg-1)*dual_mesh%nnode+1
          inods2d=ipoinpg*dual_mesh%nnode
       else
          inods1d=dual_mesh%pnods(ipoinpg)
          inods2d=dual_mesh%pnods(ipoinpg+1)-1
       end if
       do p_ipoindm = inods1d, inods2d
       !do p_ipoindm = dual_mesh%pnods(ipoinpg), dual_mesh%pnods(ipoinpg+1)-1
          ipoindm = dual_mesh%lnods(p_ipoindm)

          ! Traverse the primal nodes of the primal element number ipoindm
          ! (i.e., dual elements around dual node number ipoindm)
          if(primal_mesh%nelty==1) then
             inods1=(ipoindm-1)*primal_mesh%nnode+1
             inods2=ipoindm*primal_mesh%nnode
          else
             inods1=primal_mesh%pnods(ipoindm)
             inods2=primal_mesh%pnods(ipoindm+1)-1
          end if
          do p_ipoinpm = inods1,inods2
             ipoinpm = primal_mesh%lnods(p_ipoinpm) 

             ! write (*,*) ipoinpg, ipoindm, ipoinpm ! DBG: 
             ! If ipoinpm not visited yet
             if ( ws_position(ipoinpm) == 0 ) then
                ws_position(ipoinpm) = first_free_pos
                ws_neighbors(first_free_pos) = ipoinpm
                ws_freq(first_free_pos) = 1
                first_free_pos = first_free_pos + 1
                num_neighbors = num_neighbors + 1 
                ! If ipoinpm already visited
             else
                ws_freq (ws_position(ipoinpm)) = ws_freq( ws_position(ipoinpm) ) + 1
             end if

          end do
       end do

       ! write (*, '(10i10)')    ws_position (1:primal_graph%nv)  ! DBG:
       ! write (*, '(10i10)')    ws_neighbors(1:num_neighbors)    ! DBG: 
       ! write (*, '(10i10)')    ws_freq     (1:num_neighbors)    ! DBG: 

       ! Extract those neighbours which are listed at least min_freq_neig times
       ! while restoring working space to initial state
       do ineigh=1, num_neighbors
          if (ws_freq(ineigh).ge.min_freq_neig) then
             if (.not.ws_neighbors(ineigh)==ipoinpg) then ! Exclude self-edges
                primal_graph%ia(ipoinpg+1) = primal_graph%ia(ipoinpg+1) + 1
             end if
          end if
          ws_freq(ineigh) = 0
          ws_position ( ws_neighbors(ineigh) ) = 0
          ws_neighbors(ineigh) = 0
       end do
    end do

  end subroutine count_primal_graph_part

  !============================================================================================
  subroutine  list_primal_graph_part ( primal_mesh, dual_mesh, min_freq_neig, primal_graph,  &
       &                                 ws_position, ws_freq, ws_neighbors )
    implicit none
    ! Parameters
    type(fem_mesh) , intent(in)    :: primal_mesh, dual_mesh
    type(fem_graph), intent(inout) :: primal_graph
    integer(ip), intent(in)        :: min_freq_neig
    integer(ip), intent(out)       :: ws_position  (primal_mesh%npoin)
    integer(ip), intent(out)       :: ws_freq      (primal_mesh%nnode*dual_mesh%nnode)
    integer(ip), intent(out)       :: ws_neighbors (primal_mesh%nnode*dual_mesh%nnode)

    ! Local variables
    integer(ip)                     :: ineigh  
    integer(ip)                     :: ipoindm ! Dual   mesh point identifier
    integer(ip)                     :: ipoinpm ! Primal mesh point identifier
    integer(ip)                     :: ipoinpg ! Primal graph point identifier 
    integer(ip)                     :: p_ipoindm   
    integer(ip)                     :: p_ipoinpm   
    integer(ip)                     :: inods1,inods2
    integer(ip)                     :: inods1d, inods2d  
  

    integer(ip)                     :: first_free_ja_pos
    integer(ip)                     :: first_free_pos
    integer(ip)                     :: num_neighbors

    ! I do not agree with the following
    ! assertion. dual_mesh can actually be
    ! a primal finite element mesh when this
    ! routine has to compute dual graphs (i.e., elements 
    ! around elements). This last  situation was not taken 
    ! into account in the code below.
    ! I have modified the code below so that it
    ! is also able to compute dual graphs. 
    
    ! !!!! IMPORTANT NOTE: INCORRECT ASSERTION !!!
    ! assert(dual_mesh%nelty/=1)

    ! CORRECT ASSERTION: Either primal_mesh or dual_mesh has
    !                    to be a dual_mesh
    assert ( primal_mesh%nelty /= 1 .or. dual_mesh%nelty /= 1 ) 


    ! Initialize work space of filled-up 
    ! positions to zeros
    ws_position   = 0

    ! Initialize work space of freqs. 
    ! associated to neighbours
    ws_freq       = 0

    ! Initialize working space of neighbour 
    ! identifiers 
    ws_neighbors = 0

    first_free_ja_pos = 1

    do ipoinpg=1, primal_mesh%npoin       
       first_free_pos = 1
       num_neighbors  = 0

       ! Traverse the dual nodes of the dual element number ipoinpg
       ! (i.e., primal elements around primal node number ipoinpg)
       if(dual_mesh%nelty==1) then
          inods1d=(ipoinpg-1)*dual_mesh%nnode+1
          inods2d=ipoinpg*dual_mesh%nnode
       else
          inods1d=dual_mesh%pnods(ipoinpg)
          inods2d=dual_mesh%pnods(ipoinpg+1)-1
       end if

       do p_ipoindm = inods1d, inods2d
       ! do p_ipoindm = dual_mesh%pnods(ipoinpg), dual_mesh%pnods(ipoinpg+1)-1
          ipoindm = dual_mesh%lnods(p_ipoindm)

          ! Traverse the primal nodes of the primal element number ipoindm
          ! (i.e., dual elements around dual node number ipoindm)
          if(primal_mesh%nelty==1) then
             inods1=(ipoindm-1)*primal_mesh%nnode+1
             inods2=ipoindm*primal_mesh%nnode
          else
             inods1=primal_mesh%pnods(ipoindm)
             inods2=primal_mesh%pnods(ipoindm+1)-1
          end if
          do p_ipoinpm = inods1,inods2
             ipoinpm = primal_mesh%lnods(p_ipoinpm)
             ! If ipoinpm not visited yet
             if ( ws_position(ipoinpm) == 0 ) then
                ws_position(ipoinpm) = first_free_pos
                ws_neighbors(first_free_pos) = ipoinpm
                ws_freq(first_free_pos) = 1
                first_free_pos = first_free_pos + 1
                num_neighbors = num_neighbors + 1 
                ! If ipoinpm already visited
             else
                ws_freq ( ws_position(ipoinpm) ) = ws_freq( ws_position(ipoinpm) ) + 1
             end if
          end do
       end do

       ! Extract those neighbours which are listed at least min_freq_neig times
       ! while restoring working space to initial state
       do ineigh=1, num_neighbors
          if (ws_freq(ineigh).ge.min_freq_neig) then
             if ( .not.ws_neighbors(ineigh)==ipoinpg ) then ! Exclude self-edges
                primal_graph%ja(first_free_ja_pos) = ws_neighbors(ineigh)
                first_free_ja_pos = first_free_ja_pos + 1
             end if
          end if
          ws_freq(ineigh) = 0
          ws_position ( ws_neighbors(ineigh) ) = 0
          ws_neighbors( ineigh               ) = 0
       end do
    end do

  end subroutine list_primal_graph_part


!!$  !================================================================================================
!!$  subroutine mesh_to_graph_part_new (primal_mesh, dual_mesh, primal_graph, only_faces)
!!$    !---------------------------------------------------------------------------
!!$    ! This routine generates a graph of a given primal mesh (this graph 
!!$    ! is a list of primal points around primal points)
!!$    !
!!$    ! The routine exploits the duality among primal/dual meshes. In order to 
!!$    ! work properly pnods has to be an array of pointers to address lnods
!!$    ! (!!! this does not currently happen with FE meshes with eltyp==1 !!!)
!!$    !
!!$    ! min_freq_neig is an input integer which determines which are
!!$    ! the neighbours of a particular vertex of the graph. In particular, 
!!$    ! neighbours are those vertices which are listed at least min_freq_neig times 
!!$    ! when traversing the list of primal elements around a primal point
!!$    !
!!$    ! IMPORTANT NOTE:  This routine DOES NOT generate self-edges for both primal 
!!$    !                  and dual graphs. However, mesh_to_graph_dual
!!$    !                  generates self-edges. METIS graphs do NOT have self-edges.
!!$    !---------------------------------------------------------------------------
!!$    implicit none
!!$
!!$    ! Parameters
!!$    type(fem_mesh) , intent(in)  :: primal_mesh, dual_mesh
!!$    type(fem_graph), intent(out) :: primal_graph
!!$    integer(ip), intent(in)      :: only_faces 
!!$
!!$    ! Local variables
!!$    integer(ip), allocatable :: iwork(:)            ! Integer ip working array
!!$    integer(ip)              :: pwork(3)            ! Pointers to work space
!!$
!!$    ! Allocate space for ia on the primal graph
!!$    primal_graph%nv = primal_mesh%npoin
!!$    call memalloc (primal_graph%nv+1, primal_graph%ia, __FILE__,__LINE__)
!!$
!!$    ! Allocate working space for count_primal_graph and list_primal_graph routines
!!$    ! (TOTAL WS SIZE = primal mesh npoin + 2*maximum number of neighbours of any primal
!!$    ! graph node)
!!$    pwork(1) = 1
!!$    pwork(2) = pwork(1) + primal_mesh%npoin
!!$    pwork(3) = pwork(2) + dual_mesh%nnode*primal_mesh%nnode
!!$    call memalloc (pwork(3), iwork, __FILE__,__LINE__)
!!$
!!$    ! write(*,*) 'begin count_primal_graph' ! DBG:
!!$    ! Calculate ia array
!!$    call count_primal_graph_part_new ( primal_mesh, dual_mesh, only_faces, primal_graph,  &
!!$         &                           iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)))
!!$    ! write(*,*) 'end count_primal_graph'   ! DBG:
!!$
!!$    ! Allocate space for ja on the primal graph 
!!$    call memalloc (primal_graph%ia(primal_graph%nv+1)-1, primal_graph%ja,            __FILE__,__LINE__)
!!$
!!$    ! write(*,*) 'begin list_primal_graph' ! DBG:
!!$    ! Calculate ja array
!!$    call list_primal_graph_part_new  (primal_mesh, dual_mesh, only_faces, primal_graph, & 
!!$         &                        iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)))
!!$    ! write(*,*) 'end list_primal_graph'   ! DBG: 
!!$
!!$    ! Deallocate working space
!!$    call memfree (iwork,__FILE__,__LINE__)
!!$    return
!!$  end subroutine mesh_to_graph_part_new
!!$
!!$  !============================================================================================
!!$  subroutine  count_primal_graph_part_new ( primal_mesh, dual_mesh, only_faces, primal_graph,  &
!!$       &                                  ws_position, ws_neighbors)
!!$    implicit none
!!$
!!$    ! Parameters
!!$    type(fem_mesh) , intent(in)     :: primal_mesh, dual_mesh
!!$    integer(ip), intent(in)         :: only_faces 
!!$    type(fem_graph), intent(inout)  :: primal_graph
!!$    integer(ip), intent(out)        :: ws_position  (primal_mesh%npoin)
!!$    integer(ip), intent(out)        :: ws_neighbors (primal_mesh%nnode*dual_mesh%nnode)
!!$
!!$    ! Local variables
!!$    integer(ip)                 :: ineigh  
!!$    integer(ip)                 :: ipoindm ! Dual   mesh point identifier
!!$    integer(ip)                 :: ipoinpm ! Primal mesh point identifier
!!$    integer(ip)                 :: ipoinpg ! Primal graph point identifier 
!!$    integer(ip)                 :: p_ipoindm   
!!$    integer(ip)                 :: p_ipoinpm   
!!$    integer(ip)                 :: inods1 , inods2   
!!$    integer(ip)                 :: inods1d, inods2d
!!$
!!$    integer(ip)                 :: first_free_pos
!!$    integer(ip)                 :: num_neighbors
!!$
!!$    
!!$    assert ( only_faces == 0 .or. only_faces == 1 ) 
!!$
!!$    ! Initialize work space of filled-up 
!!$    ! positions
!!$    ws_position   = 0
!!$
!!$    ! Initialize work space of neighbour 
!!$    ! identifiers 
!!$    ws_neighbors  = 0
!!$
!!$    !primal_graph%ia    = -1                                                ! DBG:
!!$    !write (*,*) size(primal_graph%ia)                                      ! DBG: 
!!$    !write (*,*) size(primal_mesh%pnods)                                    ! DBG: 
!!$    !write (*,*) size(primal_mesh%lnods)                                    ! DBG: 
!!$    !write (*,*) size(dual_mesh%pnods)                                      ! DBG: 
!!$    !write (*,*) size(dual_mesh%lnods)                                      ! DBG: 
!!$    !write (*,*) size(ws_position)                                          ! DBG: 
!!$    !write (*,*) size(ws_neighbors)                                         ! DBG: 
!!$    !write (*, '(10i10)')    primal_mesh%pnods  (1:(primal_mesh%nelem+1))   ! DBG: 
!!$    !write (*, '(10i10)')    primal_graph%ia    (1:(primal_graph%nv+1))     ! DBG:  
!!$    !primal_graph%ia    = -1 
!!$    primal_graph%ia(1) = 1
!!$
!!$    do ipoinpg=1, primal_mesh%npoin
!!$       primal_graph%ia(ipoinpg+1) = primal_graph%ia(ipoinpg)
!!$       first_free_pos = 1
!!$       num_neighbors  = 0
!!$       ! Traverse the dual nodes of the dual element number ipoinpg
!!$       ! (i.e., primal elements around primal node number ipoinpg)
!!$       inods1d=dual_mesh%pnods(ipoinpg)
!!$       inods2d=dual_mesh%pnods(ipoinpg+1)-1
!!$       do p_ipoindm = inods1d, inods2d
!!$       !do p_ipoindm = dual_mesh%pnods(ipoinpg), dual_mesh%pnods(ipoinpg+1)-1
!!$          ipoindm = dual_mesh%lnods(p_ipoindm)
!!$
!!$          if (only_faces==0) then 
!!$             ! Traverse the primal nodes of the primal element number ipoindm
!!$             ! (i.e., dual elements around dual node number ipoindm)
!!$             inods1=primal_mesh%pnods(ipoindm)
!!$             inods2=primal_mesh%pnods(ipoindm+1)-1
!!$          else if ( only_faces == 1 ) then
!!$             inods1=primal_mesh%p_face(ipoindm)
!!$             inods2=primal_mesh%pnods(ipoindm+1)-1
!!$          end if
!!$
!!$          do p_ipoinpm = inods1,inods2
!!$             ipoinpm = primal_mesh%lnods(p_ipoinpm) 
!!$
!!$             ! write (*,*) ipoinpg, ipoindm, ipoinpm ! DBG: 
!!$             ! If ipoinpm not visited yet
!!$             if ( ws_position(ipoinpm) == 0 ) then
!!$                ws_position(ipoinpm) = first_free_pos
!!$                ws_neighbors(first_free_pos) = ipoinpm
!!$                first_free_pos = first_free_pos + 1
!!$                num_neighbors = num_neighbors + 1 
!!$                ! If ipoinpm already visited
!!$             end if
!!$
!!$          end do
!!$       end do
!!$
!!$       ! write (*, '(10i10)')    ws_position (1:primal_graph%nv)  ! DBG:
!!$       ! write (*, '(10i10)')    ws_neighbors(1:num_neighbors)    ! DBG: 
!!$
!!$       ! Extract those neighbours which are listed at least min_freq_neig times
!!$       ! while restoring working space to initial state
!!$       do ineigh=1, num_neighbors
!!$          if (.not.ws_neighbors(ineigh)==ipoinpg) then ! Exclude self-edges
!!$              primal_graph%ia(ipoinpg+1) = primal_graph%ia(ipoinpg+1) + 1
!!$          end if
!!$          ws_position ( ws_neighbors(ineigh) ) = 0
!!$          ws_neighbors(ineigh) = 0
!!$       end do
!!$    end do
!!$
!!$  end subroutine count_primal_graph_part_new
!!$
!!$  !============================================================================================
!!$  subroutine  list_primal_graph_part_new ( primal_mesh, dual_mesh, only_faces, primal_graph,  &
!!$       &                                 ws_position, ws_neighbors )
!!$    implicit none
!!$    ! Parameters
!!$    type(fem_mesh) , intent(in)    :: primal_mesh, dual_mesh
!!$    type(fem_graph), intent(inout) :: primal_graph
!!$    integer(ip), intent(in)        :: only_faces 
!!$    integer(ip), intent(out)       :: ws_position  (primal_mesh%npoin)
!!$    integer(ip), intent(out)       :: ws_neighbors (primal_mesh%nnode*dual_mesh%nnode)
!!$
!!$    ! Local variables
!!$    integer(ip)                     :: ineigh  
!!$    integer(ip)                     :: ipoindm ! Dual   mesh point identifier
!!$    integer(ip)                     :: ipoinpm ! Primal mesh point identifier
!!$    integer(ip)                     :: ipoinpg ! Primal graph point identifier 
!!$    integer(ip)                     :: p_ipoindm   
!!$    integer(ip)                     :: p_ipoinpm   
!!$    integer(ip)                     :: inods1,inods2
!!$    integer(ip)                     :: inods1d, inods2d  
!!$  
!!$
!!$    integer(ip)                     :: first_free_ja_pos
!!$    integer(ip)                     :: first_free_pos
!!$    integer(ip)                     :: num_neighbors
!!$
!!$    assert ( only_faces == 0 .or. only_faces == 1 ) 
!!$
!!$    ! Initialize work space of filled-up 
!!$    ! positions to zeros
!!$    ws_position   = 0
!!$
!!$    ! Initialize working space of neighbour 
!!$    ! identifiers 
!!$    ws_neighbors = 0
!!$
!!$    first_free_ja_pos = 1
!!$
!!$    do ipoinpg=1, primal_mesh%npoin       
!!$       first_free_pos = 1
!!$       num_neighbors  = 0
!!$
!!$       ! Traverse the dual nodes of the dual element number ipoinpg
!!$       ! (i.e., primal elements around primal node number ipoinpg)
!!$       inods1d=dual_mesh%pnods(ipoinpg)
!!$       inods2d=dual_mesh%pnods(ipoinpg+1)-1
!!$
!!$       do p_ipoindm = inods1d, inods2d
!!$       ! do p_ipoindm = dual_mesh%pnods(ipoinpg), dual_mesh%pnods(ipoinpg+1)-1
!!$          ipoindm = dual_mesh%lnods(p_ipoindm)
!!$
!!$          ! Traverse the primal nodes of the primal element number ipoindm
!!$          ! (i.e., dual elements around dual node number ipoindm)
!!$          if (only_faces==0) then
!!$             ! Traverse the primal nodes of the primal element number ipoindm
!!$             ! (i.e., dual elements around dual node number ipoindm)
!!$             inods1=primal_mesh%pnods(ipoindm)
!!$             inods2=primal_mesh%pnods(ipoindm+1)-1
!!$          else if ( only_faces == 1 ) then
!!$             inods1=primal_mesh%p_face(ipoindm)
!!$             inods2=primal_mesh%pnods(ipoindm+1)-1
!!$          end if
!!$
!!$          do p_ipoinpm = inods1,inods2
!!$             ipoinpm = primal_mesh%lnods(p_ipoinpm)
!!$             ! If ipoinpm not visited yet
!!$             if ( ws_position(ipoinpm) == 0 ) then
!!$                ws_position(ipoinpm) = first_free_pos
!!$                ws_neighbors(first_free_pos) = ipoinpm
!!$                first_free_pos = first_free_pos + 1
!!$                num_neighbors = num_neighbors + 1 
!!$                ! If ipoinpm already visited
!!$             end if
!!$          end do
!!$       end do
!!$
!!$       ! Extract those neighbours which are listed at least min_freq_neig times
!!$       ! while restoring working space to initial state
!!$       do ineigh=1, num_neighbors
!!$          if ( .not.ws_neighbors(ineigh)==ipoinpg ) then ! Exclude self-edges
!!$             primal_graph%ja(first_free_ja_pos) = ws_neighbors(ineigh)
!!$             first_free_ja_pos = first_free_ja_pos + 1
!!$          end if
!!$          ws_position ( ws_neighbors(ineigh) ) = 0
!!$          ws_neighbors( ineigh               ) = 0
!!$       end do
!!$    end do
!!$
!!$  end subroutine list_primal_graph_part_new

  !============================================================================================
  subroutine mesh_to_graph_part_new (primal_mesh, dual_mesh, dual_graph, only_faces)
    !---------------------------------------------------------------------------
    ! This routine generates a graph of a given primal mesh (this graph 
    ! is a list of primal points around primal points)
    !
    ! The routine exploits the duality among primal/dual meshes. In order to 
    ! work properly pnods has to be an array of pointers to address lnods
    ! (!!! this does not currently happen with FE meshes with eltyp==1 !!!)
    !
    ! min_freq_neig is an input integer which determines which are
    ! the neighbours of a particular vertex of the graph. In particular, 
    ! neighbours are those vertices which are listed at least min_freq_neig times 
    ! when traversing the list of primal elements around a primal point
    !
    ! IMPORTANT NOTE:  This routine DOES NOT generate self-edges for both primal 
    !                  and dual graphs. However, mesh_to_graph_dual
    !                  generates self-edges. METIS graphs do NOT have self-edges.
    !---------------------------------------------------------------------------
    implicit none

    ! Parameters
    type(fem_mesh) , intent(in)  :: primal_mesh, dual_mesh
    type(fem_graph), intent(out) :: dual_graph
    integer(ip), intent(in)      :: only_faces 

    ! Local variables
    integer(ip), allocatable :: iwork(:)            ! Integer ip working array
    integer(ip)              :: pwork(3)            ! Pointers to work space

    assert(only_faces==1)

    ! Allocate space for ia on the dual graph
    dual_graph%nv = dual_mesh%npoin
    call memalloc (dual_graph%nv+1, dual_graph%ia, __FILE__,__LINE__)

    ! Allocate working space for count_dual_graph and list_dual_graph routines
    ! (TOTAL WS SIZE = dual mesh npoin + 2*maximum number of neighbours of any dual
    ! graph node)
    pwork(1) = 1
    pwork(2) = pwork(1) + dual_mesh%npoin
    pwork(3) = pwork(2) + dual_mesh%nnode*primal_mesh%nnode
    call memalloc (pwork(3), iwork, __FILE__,__LINE__)

    ! write(*,*) 'begin count_dual_graph' ! DBG:
    ! Calculate ia array
    call count_dual_graph_part_new ( primal_mesh, dual_mesh, dual_graph,iwork(pwork(1):pwork(2)), &
         &                           iwork(pwork(2):pwork(3)))
    ! write(*,*) 'end count_dual_graph'   ! DBG:

    ! Allocate space for ja on the dual graph 
    call memalloc (dual_graph%ia(dual_graph%nv+1)-1, dual_graph%ja,__FILE__,__LINE__)

    ! write(*,*) 'begin list_dual_graph' ! DBG:
    ! Calculate ja array
    call list_dual_graph_part_new (primal_mesh, dual_mesh, dual_graph, iwork(pwork(1):pwork(2)), & 
         &                         iwork(pwork(2):pwork(3)))
    ! write(*,*) 'end list_dual_graph'   ! DBG: 

    ! Deallocate working space
    call memfree (iwork,__FILE__,__LINE__)
    return
  end subroutine mesh_to_graph_part_new

  !============================================================================================
  subroutine  count_dual_graph_part_new ( primal_mesh, dual_mesh, dual_graph, ws_position, &
       &                                  ws_neighbors)
    implicit none

    ! Parameters
    type(fem_mesh) , intent(in)     :: primal_mesh, dual_mesh
    type(fem_graph), intent(inout)  :: dual_graph
    integer(ip), intent(out)        :: ws_position  (dual_mesh%npoin)
    integer(ip), intent(out)        :: ws_neighbors (primal_mesh%nnode*dual_mesh%nnode)

    ! Local variables
    integer(ip)                 :: ineigh  
    integer(ip)                 :: ipoinpm ! Primal mesh point identifier
    integer(ip)                 :: ipoindm ! Dual mesh point identifier
    integer(ip)                 :: ipoindg ! Dual graph point identifier 
    integer(ip)                 :: p_ipoindm   
    integer(ip)                 :: p_ipoinpm   
    integer(ip)                 :: inods1 , inods2   
    integer(ip)                 :: inods1d, inods2d

    integer(ip)                 :: first_free_pos
    integer(ip)                 :: num_neighbors

    ! Initialize work space of filled-up 
    ! positions
    ws_position   = 0

    ! Initialize work space of neighbour 
    ! identifiers 
    ws_neighbors  = 0

    !primal_graph%ia    = -1                                                ! DBG:
    !write (*,*) size(primal_graph%ia)                                      ! DBG: 
    !write (*,*) size(primal_mesh%pnods)                                    ! DBG: 
    !write (*,*) size(primal_mesh%lnods)                                    ! DBG: 
    !write (*,*) size(dual_mesh%pnods)                                      ! DBG: 
    !write (*,*) size(dual_mesh%lnods)                                      ! DBG: 
    !write (*,*) size(ws_position)                                          ! DBG: 
    !write (*,*) size(ws_neighbors)                                         ! DBG: 
    !write (*, '(10i10)')    primal_mesh%pnods  (1:(primal_mesh%nelem+1))   ! DBG: 
    !write (*, '(10i10)')    primal_graph%ia    (1:(primal_graph%nv+1))     ! DBG:  
    !primal_graph%ia    = -1 
    dual_graph%ia(1) = 1

    do ipoindg=1, dual_mesh%npoin
       dual_graph%ia(ipoindg+1) = dual_graph%ia(ipoindg)
       first_free_pos = 1
       num_neighbors  = 0
       inods1=primal_mesh%p_face(ipoindg)
       inods2=primal_mesh%pnods(ipoindg+1)-1
       do p_ipoinpm = inods1,inods2
          ipoinpm = primal_mesh%lnods(p_ipoinpm) 

          inods1d=dual_mesh%pnods(ipoinpm)
          inods2d=dual_mesh%pnods(ipoinpm+1)-1
          do p_ipoindm = inods1d, inods2d
             ipoindm = dual_mesh%lnods(p_ipoindm)

             ! write (*,*) ipoinpg, ipoindm, ipoinpm ! DBG: 
             ! If ipoinpm not visited yet
             if ( ws_position(ipoindm) == 0 ) then
                ws_position(ipoindm) = first_free_pos
                ws_neighbors(first_free_pos) = ipoindm
                first_free_pos = first_free_pos + 1
                num_neighbors = num_neighbors + 1 
                ! If ipoinpm already visited
             end if

          end do
       end do

       ! write (*, '(10i10)')    ws_position (1:primal_graph%nv)  ! DBG:
       ! write (*, '(10i10)')    ws_neighbors(1:num_neighbors)    ! DBG: 

       ! Extract those neighbours which are listed at least min_freq_neig times
       ! while restoring working space to initial state
       do ineigh=1, num_neighbors
          if (.not.ws_neighbors(ineigh)==ipoindg) then ! Exclude self-edges
              dual_graph%ia(ipoindg+1) = dual_graph%ia(ipoindg+1) + 1
          end if
          ws_position ( ws_neighbors(ineigh) ) = 0
          ws_neighbors(ineigh) = 0
       end do
    end do

  end subroutine count_dual_graph_part_new

  !============================================================================================
  subroutine  list_dual_graph_part_new ( primal_mesh, dual_mesh, dual_graph, ws_position,  &
       &                                 ws_neighbors )
    implicit none
    ! Parameters
    type(fem_mesh) , intent(in)    :: primal_mesh, dual_mesh
    type(fem_graph), intent(inout) :: dual_graph
    integer(ip), intent(out)       :: ws_position  (dual_mesh%npoin)
    integer(ip), intent(out)       :: ws_neighbors (primal_mesh%nnode*dual_mesh%nnode)

    ! Local variables
    integer(ip)                     :: ineigh  
    integer(ip)                     :: ipoinpm ! Primal mesh point identifier
    integer(ip)                     :: ipoindm ! Dual mesh point identifier
    integer(ip)                     :: ipoindg ! Dual graph point identifier 
    integer(ip)                     :: p_ipoindm   
    integer(ip)                     :: p_ipoinpm   
    integer(ip)                     :: inods1,inods2
    integer(ip)                     :: inods1d, inods2d  
  

    integer(ip)                     :: first_free_ja_pos
    integer(ip)                     :: first_free_pos
    integer(ip)                     :: num_neighbors

    ! Initialize work space of filled-up 
    ! positions to zeros
    ws_position   = 0

    ! Initialize working space of neighbour 
    ! identifiers 
    ws_neighbors = 0

    first_free_ja_pos = 1

    do ipoindg=1, dual_mesh%npoin       
       first_free_pos = 1
       num_neighbors  = 0

       inods1=primal_mesh%p_face(ipoindg)
       inods2=primal_mesh%pnods(ipoindg+1)-1
       do p_ipoinpm = inods1,inods2
          ipoinpm = primal_mesh%lnods(p_ipoinpm)

          inods1d=dual_mesh%pnods(ipoinpm)
          inods2d=dual_mesh%pnods(ipoinpm+1)-1
          do p_ipoindm = inods1d, inods2d
             ipoindm = dual_mesh%lnods(p_ipoindm)

             ! If ipoinpm not visited yet
             if ( ws_position(ipoindm) == 0 ) then
                ws_position(ipoindm) = first_free_pos
                ws_neighbors(first_free_pos) = ipoindm
                first_free_pos = first_free_pos + 1
                num_neighbors = num_neighbors + 1 
             end if
          end do
       end do

       ! Extract those neighbours which are listed at least min_freq_neig times
       ! while restoring working space to initial state
       do ineigh=1, num_neighbors
          if ( .not.ws_neighbors(ineigh)==ipoindg ) then ! Exclude self-edges
             dual_graph%ja(first_free_ja_pos) = ws_neighbors(ineigh)
             first_free_ja_pos = first_free_ja_pos + 1
          end if
          ws_position ( ws_neighbors(ineigh) ) = 0
          ws_neighbors( ineigh               ) = 0
       end do
    end do

  end subroutine list_dual_graph_part_new


  ! Inspired on http://en.wikipedia.org/wiki/Breadth-first_search.
  ! Given a mesh (m) and its dual graph (g), it computes the list 
  ! of nodes (lconn) of each connected component in m
  subroutine mesh_graph_compute_connected_components (m, g, lconn)
    implicit none

    ! Parameters
    type(fem_mesh) , intent(in)   :: m   
    type(fem_graph), intent(in)   :: g
    type(list)     , intent(out)  :: lconn
 
    ! Locals
    integer(ip), allocatable :: auxv(:), auxe(:), e(:)
    integer(ip), allocatable :: emarked(:), vmarked(:)
    integer(ip), allocatable :: q(:)
    integer(ip)              :: head, tail, i, esize, vsize, current, & 
                                j, l, k, inods1d, inods2d, p_ipoin, ipoin

    call memalloc ( g%nv   , auxe     , __FILE__,__LINE__)
    call memalloc ( g%nv   , auxv     , __FILE__,__LINE__)
    call memalloc ( g%nv   , q        , __FILE__,__LINE__)
    call memalloc ( g%nv   , emarked  , __FILE__,__LINE__)
    call memalloc ( m%npoin, vmarked  , __FILE__,__LINE__)
    call memalloc ( g%nv  ,  e        , __FILE__,__LINE__)

    lconn%n  = 0
    emarked  = 0
    current  = 1 

    do i=1, g%nv
       if (emarked(i) == 0) then
          ! New connected component
          lconn%n = lconn%n +1
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

             ! Traverse the nodes of the element number j
             if( m%nelty == 1 ) then
                inods1d = (j-1) * m%nnode+1
                inods2d = j     * m%nnode
             else
                inods1d = m%pnods(j)
                inods2d = m%pnods(j+1)-1
             end if

             do p_ipoin = inods1d, inods2d
                ipoin = m%lnods(p_ipoin)
                if (vmarked(ipoin)==0) then
                   vmarked(ipoin)=1
                   vsize = vsize+1
                end if
             end do

!!$9         for all edges e in G.adjacentEdges(t) do
             do k=g%ia(j), g%ia(j+1)-1
!!$12           u ← G.adjacentVertex(t,e)
                l=g%ja(k)
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
             end do
          end do
          auxe(lconn%n) = esize
          auxv(lconn%n) = vsize
       end if
    end do
    
    call memalloc( lconn%n+1, lconn%p, __FILE__,__LINE__)

    lconn%p(1) = 1
    do i=1, lconn%n
       lconn%p(i+1) = lconn%p(i) + auxv(i)
    end do

    call memfree( auxv   ,__FILE__,__LINE__)
    call memfree( q      ,__FILE__,__LINE__)
    call memfree( emarked,__FILE__,__LINE__)

    call memalloc( lconn%p(lconn%n+1)-1, lconn%l, __FILE__,__LINE__)

    current=1
    l=1
    do i=1, lconn%n
       vmarked = 0
       ! Traverse elements of current connected component  
       do current=current,current+auxe(i)-1
          j=e(current)

          ! Traverse the nodes of the element number j
          if( m%nelty == 1 ) then
             inods1d = (j-1) * m%nnode+1
             inods2d = j     * m%nnode
          else
             inods1d = m%pnods(j)
             inods2d = m%pnods(j+1)-1
          end if

          do p_ipoin = inods1d, inods2d
             ipoin = m%lnods(p_ipoin)
             if (vmarked(ipoin)==0) then
                vmarked(ipoin)=1
                lconn%l(l)=ipoin
                l=l+1
             end if
          end do
          
       end do
    end do

    ! write(*,*) 'ZZ', g%nv, m%npoin, l
    ! write(*,*) 'XX', lconn%n
    ! write(*,*) 'YY', lconn%p
    ! write(*,*) 'PP', lconn%l

    call memfree( auxe,__FILE__,__LINE__)
    call memfree( e,__FILE__,__LINE__)
    call memfree( vmarked,__FILE__,__LINE__)

  end subroutine mesh_graph_compute_connected_components

  !================================================================================================
  subroutine mesh_to_graph_share (primal_mesh, dual_mesh, primal_graph, edges)
    !---------------------------------------------------------------------------
    ! This routine generates a graph of a given primal mesh (this graph 
    ! is a list of primal points around primal points)
    !
    ! Unlike previous routines min_freq_neig is not used, this function only
    ! makes sense for min_freq_neig = 1. 
    !
    ! IMPORTANT NOTE:  This routine DOES NOT generate self-edges for both primal 
    !                  and dual graphs. However, mesh_to_graph_dual
    !                  generates self-edges. METIS graphs do NOT have self-edges.
    !---------------------------------------------------------------------------
    implicit none

    ! Parameters
    type(fem_mesh) , intent(in)  :: primal_mesh, dual_mesh
    !integer(ip)    , intent(in)  :: min_freq_neig
    integer(ip)    , parameter   :: min_freq_neig = 1
    type(fem_graph), intent(out) :: primal_graph
    integer(ip)    , allocatable, intent(out) :: edges(:)

    ! Local variables
    integer(ip), allocatable :: iwork(:)            ! Integer ip working array
    integer(ip)              :: pwork(4)            ! Pointers to work space

    ! Allocate space for ia on the primal graph
    primal_graph%nv = primal_mesh%npoin
    call memalloc (primal_graph%nv+1, primal_graph%ia, __FILE__,__LINE__)

    ! Allocate working space for count_primal_graph and list_primal_graph routines
    ! (TOTAL WS SIZE = primal mesh npoin + 2*maximum number of neighbours of any primal
    ! graph node)
    pwork(1) = 1
    pwork(2) = pwork(1) + primal_mesh%npoin
    pwork(3) = pwork(2) + dual_mesh%nnode*primal_mesh%nnode
    pwork(4) = pwork(3) + dual_mesh%nnode*primal_mesh%nnode
    call memalloc (pwork(4), iwork, __FILE__,__LINE__)

    ! write(*,*) 'begin count_primal_graph' ! DBG:
    ! Calculate ia array
    call count_primal_graph_part ( primal_mesh, dual_mesh, min_freq_neig, primal_graph,  &
         &                         iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)), &
         &                         iwork(pwork(3):pwork(4)) ) 
    ! write(*,*) 'end count_primal_graph'   ! DBG:

    ! Allocate space for ja on the primal graph 
    call memalloc (primal_graph%ia(primal_graph%nv+1)-1, primal_graph%ja,            __FILE__,__LINE__)
    call memalloc (primal_graph%ia(primal_graph%nv+1)-1, edges,            __FILE__,__LINE__)

    ! write(*,*) 'begin list_primal_graph' ! DBG:
    ! Calculate ja array
    call list_primal_graph_share  (primal_mesh, dual_mesh, min_freq_neig, primal_graph, edges, & 
         &                        iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)),  &
         &                        iwork(pwork(3):pwork(4)) )
    ! write(*,*) 'end list_primal_graph'   ! DBG: 

    ! Deallocate working space
    call memfree (iwork,__FILE__,__LINE__)
    return
  end subroutine mesh_to_graph_share

  !============================================================================================
  subroutine  list_primal_graph_share ( primal_mesh, dual_mesh, min_freq_neig, primal_graph, edges, &
       &                                 ws_position, ws_freq, ws_neighbors )
    implicit none
    ! Parameters
    type(fem_mesh) , intent(in)    :: primal_mesh, dual_mesh
    type(fem_graph), intent(inout) :: primal_graph
    integer(ip), intent(in)        :: min_freq_neig
    integer(ip), intent(out)       :: edges        (primal_graph%ia(primal_graph%nv+1)-1)
    integer(ip), intent(out)       :: ws_position  (primal_mesh%npoin)
    integer(ip), intent(out)       :: ws_freq      (primal_mesh%nnode*dual_mesh%nnode)
    integer(ip), intent(out)       :: ws_neighbors (primal_mesh%nnode*dual_mesh%nnode)

    ! Local variables
    integer(ip)                     :: ineigh  
    integer(ip)                     :: ipoindm ! Dual   mesh point identifier
    integer(ip)                     :: ipoinpm ! Primal mesh point identifier
    integer(ip)                     :: ipoinpg ! Primal graph point identifier 
    integer(ip)                     :: p_ipoindm   
    integer(ip)                     :: p_ipoinpm   
    integer(ip)                     :: inods1,inods2
    integer(ip)                     :: inods1d, inods2d  
  

    integer(ip)                     :: first_free_ja_pos
    integer(ip)                     :: first_free_pos
    integer(ip)                     :: num_neighbors
    integer(ip)                     :: freq_neig

 
    assert ( primal_mesh%nelty /= 1 .or. dual_mesh%nelty /= 1 ) 

    ! Initialize work space of filled-up positions to zeros
    ws_position   = 0
    ! Initialize work space of freqs. associated to neighbours
    ws_freq       = 0
    ! Initialize working space of neighbour identifiers 
    ws_neighbors = 0

    first_free_ja_pos = 1

    do ipoinpg=1, primal_mesh%npoin   ! do ielem 
       first_free_pos = 1
       num_neighbors  = 0

       if(dual_mesh%nelty==1) then
          inods1d=(ipoinpg-1)*dual_mesh%nnode+1
          inods2d=ipoinpg*dual_mesh%nnode
       else
          inods1d=dual_mesh%pnods(ipoinpg)
          inods2d=dual_mesh%pnods(ipoinpg+1)-1
       end if

       do p_ipoindm = inods1d, inods2d   ! do inode in ielem

          ipoindm = dual_mesh%lnods(p_ipoindm)

          if(primal_mesh%nelty==1) then
             inods1=(ipoindm-1)*primal_mesh%nnode+1
             inods2=ipoindm*primal_mesh%nnode
          else
             inods1=primal_mesh%pnods(ipoindm)
             inods2=primal_mesh%pnods(ipoindm+1)-1
          end if
          do p_ipoinpm = inods1,inods2     ! do elems around inode
             ipoinpm = primal_mesh%lnods(p_ipoinpm)   ! jelem
             ! If ipoinpm not visited yet
             if ( ws_position(ipoinpm) == 0 ) then
                ws_position(ipoinpm) = first_free_pos
                ws_neighbors(first_free_pos) = ipoinpm
                ! This neighbor is visited through p_ipoindm
                ws_freq(first_free_pos) = ibset(ws_freq(first_free_pos),p_ipoindm-inods1d)
                !ws_freq(first_free_pos) = 1
                first_free_pos = first_free_pos + 1
                num_neighbors = num_neighbors + 1 
                ! If ipoinpm already visited
             else
                ! This neighbor is (also) visited through p_ipoindm
                ws_freq(ws_position(ipoinpm)) = ibset(ws_freq(ws_position(ipoinpm)),p_ipoindm-inods1d)
                !ws_freq ( ws_position(ipoinpm) ) = ws_freq( ws_position(ipoinpm) ) + 1
             end if
          end do
       end do

       ! Extract those neighbours which are listed at least min_freq_neig times
       ! while restoring working space to initial state
       do ineigh=1, num_neighbors
          !! In a general case of min_freq_neig/=1 we need to count how many times
          !! this neighbor has been visited, i.e sum ws_freq bits set to 1.
          !freq_neig = 0
          !do p_ipoindm = inods1d, inods2d
          !   freq_neig = freq_neig + btest(ws_freq(ineigh),p_ipoindm)
          !end do
          !if (freq_neig.ge.min_freq_neig) then
          if (ws_freq(ineigh) /= 0 ) then  ! It has been visited
             if ( .not.ws_neighbors(ineigh)==ipoinpg ) then ! Exclude self-edges
                primal_graph%ja(first_free_ja_pos) = ws_neighbors(ineigh)
                edges(first_free_ja_pos) = ws_freq(ineigh)
                first_free_ja_pos = first_free_ja_pos + 1
             end if
          end if
          ws_freq(ineigh) = 0
          ws_position ( ws_neighbors(ineigh) ) = 0
          ws_neighbors( ineigh               ) = 0
       end do
    end do

  end subroutine list_primal_graph_share

 !=============================================================================
  subroutine szmap_to_dual(primal_lnods,npoin,nelem,dual_pnods,dual_lnods)
    !-----------------------------------------------------------------------
    ! This routine generates the dual szmap (list of nodes that point to
    ! a element) of a given primal szmap. 
    !-----------------------------------------------------------------------
    implicit none
    integer(ip)   , intent(in)  :: primal_lnods(npoin)
    integer(ip)   , intent(in)  :: npoin,nelem
    integer(ip)   , intent(out) :: dual_pnods(nelem+1),dual_lnods(:)

    integer(ip)   :: primal_pnods(nelem),ielem,nnode

    primal_pnods = (/(ielem, ielem=1,nelem)/)

    ! Count elements around points and generate pointers to store them (pelpo)
    call count_elements_around_points(1,nelem,npoin,primal_pnods,primal_lnods,1,nnode,dual_pnods, NULL())


    ! List elements around points (lelpo)
    call list_elements_around_points(1,nelem,npoin,primal_pnods,primal_lnods,1,nnode,dual_pnods,dual_lnods, NULL())

    return

  end subroutine szmap_to_dual

end module mesh_graph
