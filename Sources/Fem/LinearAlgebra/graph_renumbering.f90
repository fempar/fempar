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
module graph_renumbering_names
  use types_names
  use memor_names
  use renumbering_names
  use graph_names
  use metis_interface_names
  use mesh_partition_base_names
  use rcm_renumbering_names
#include "debug.i90"
  implicit none
  private
  
  ! Functions
  public :: graph_pt_renumbering, graph_nd_renumbering
  
contains

  !=================================================================================================
  subroutine graph_nd_renumbering(prt_parts, gp, renumbering)
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    implicit none
    type(part_params_t), intent(in)        :: prt_parts
    type(graph_t), target, intent(in)  :: gp
    type(renumbering_t), target,  intent(inout)  :: renumbering
    
    assert(renumbering%n==gp%nv)
    
    if ( gp%nv == 1 ) then
       renumbering%lperm(1) = 1
       renumbering%iperm(1) = 1
    else
#ifdef ENABLE_METIS
       ierr = metis_setdefaultoptions(c_loc(options))
       assert(ierr == METIS_OK) 
       
       options(METIS_OPTION_NUMBERING) = 1
       options(METIS_OPTION_DBGLVL)    = prt_parts%metis_option_debug
       
       ierr = metis_nodend ( c_loc(gp%nv),c_loc(gp%ia),c_loc(gp%ja),C_NULL_PTR,c_loc(options), &
            &                c_loc(renumbering%iperm),c_loc(renumbering%lperm))
       
       assert(ierr == METIS_OK)
#else
       call enable_metis_error_message
#endif
    end if
  end subroutine graph_nd_renumbering

  !=================================================================================================
  subroutine graph_pt_renumbering(prt_parts,gp,ldomn)
    !-----------------------------------------------------------------------
    ! This routine computes a nparts-way-partitioning of the input graph gp
    !-----------------------------------------------------------------------
    implicit none
    type(part_params_t), target, intent(in)    :: prt_parts
    type(graph_t)  , target, intent(inout) :: gp
    integer(ip)      , target, intent(out)     :: ldomn(gp%nv)

    ! Local variables 
    integer(ip), target      :: kedge
    integer(ip)              :: idumm,iv
    integer(ip), allocatable :: lwork(:)
    integer(ip)              :: i, j, m, k, ipart
    integer(ip), allocatable :: iperm(:)

   
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

    if ( prt_parts%strat == part_kway ) then
       options(METIS_OPTION_NUMBERING) = 1
       options(METIS_OPTION_DBGLVL)    = prt_parts%metis_option_debug
       
       ! Enforce contiguous partititions
       options(METIS_OPTION_CONTIG)    = prt_parts%metis_option_contig
       
       ! Explicitly minimize the maximum degree of the subdomain graph
       options(METIS_OPTION_MINCONN)   = prt_parts%metis_option_minconn
       options(METIS_OPTION_UFACTOR)   = prt_parts%metis_option_ufactor
       
       ncon = 1 
       ierr = metis_partgraphkway( c_loc(gp%nv), c_loc(ncon), c_loc(gp%ia)  , c_loc(gp%ja) , & 
                                   C_NULL_PTR  , C_NULL_PTR , C_NULL_PTR    , c_loc(prt_parts%nparts), &
                                   C_NULL_PTR  , C_NULL_PTR , c_loc(options), c_loc(kedge), c_loc(ldomn) )
       
       assert(ierr == METIS_OK) 
       
    else if ( prt_parts%strat == part_recursive ) then
       options(METIS_OPTION_NUMBERING) = 1
       options(METIS_OPTION_DBGLVL)    = prt_parts%metis_option_debug
       options(METIS_OPTION_UFACTOR)   = prt_parts%metis_option_ufactor

       ncon = 1 
       ierr = metis_partgraphrecursive( c_loc(gp%nv), c_loc(ncon), c_loc(gp%ia)  , c_loc(gp%ja) , & 
                                        C_NULL_PTR  , C_NULL_PTR , C_NULL_PTR    , c_loc(prt_parts%nparts), &
                                        C_NULL_PTR  , C_NULL_PTR , c_loc(options), c_loc(kedge), c_loc(ldomn) )
    end if    
#else
    call enable_metis_error_message
#endif

    if ( prt_parts%strat == part_strip ) then
       j = gp%nv
       m = 0
       do ipart=1,prt_parts%nparts
          k = j / (prt_parts%nparts-ipart+1)
          do i = 1, k
             ldomn(m+i) = ipart
          end do
          m = m + k
          j = j - k
       end do
    else if ( prt_parts%strat == part_rcm_strip ) then
       call memalloc ( gp%nv, iperm, __FILE__,__LINE__ )
       call genrcm ( gp%nv, gp%ia(gp%nv+1)-1, gp%ia, gp%ja, iperm )
       j = gp%nv
       m = 0
       do ipart=1,prt_parts%nparts
          k = j / (prt_parts%nparts-ipart+1)
          do i = 1, k
             ldomn(iperm(m+i)) = ipart
          end do
          m = m + k
          j = j - k
       end do
       call memfree ( iperm,__FILE__,__LINE__)
    end if

  end subroutine graph_pt_renumbering

end module graph_renumbering_names
