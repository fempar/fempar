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
module metis_interface
  use types
  use iso_c_binding
  implicit none

#ifdef ENABLE_METIS5

!!$/* The maximum length of the options[] array */
!!$#define METIS_NOPTIONS          40
     integer(c_int), parameter ::  METIS_NOPTIONS    =  40
 

!!$/*------------------------------------------------------------------------
!!$* Enum type definitions 
!!$*-------------------------------------------------------------------------*/
!!$/*! Return codes */
!!$typedef enum {
!!$  METIS_OK              = 1,    /*!< Returned normally */
!!$  METIS_ERROR_INPUT     = -2,   /*!< Returned due to erroneous inputs and/or options */
!!$  METIS_ERROR_MEMORY    = -3,   /*!< Returned due to insufficient memory */
!!$  METIS_ERROR           = -4    /*!< Some other errors */
!!$} rstatus_et; 
!!$
     integer(c_int), parameter ::  METIS_OK    =  1
     integer(c_int), parameter ::  METIS_ERROR_INPUT  = -2
     integer(c_int), parameter ::  METIS_ERROR_MEMORY = -3
     integer(c_int), parameter ::  METIS_ERROR = -4

!!$
!!$/*! Operation type codes */
!!$typedef enum {
!!$  METIS_OP_PMETIS,       
!!$  METIS_OP_KMETIS,
!!$  METIS_OP_OMETIS
!!$} moptype_et;
!!$
     integer(c_int), parameter :: METIS_OP_PMETIS = 0       
     integer(c_int), parameter :: METIS_OP_KMETIS = 1 
     integer(c_int), parameter :: METIS_OP_OMETIS = 2
!!$
!!$/*! Options codes (i.e., options[]) */
!!$typedef enum {
!!$  METIS_OPTION_PTYPE,
!!$  METIS_OPTION_OBJTYPE,
!!$  METIS_OPTION_CTYPE,
!!$  METIS_OPTION_IPTYPE,
!!$  METIS_OPTION_RTYPE,
!!$  METIS_OPTION_DBGLVL,
!!$  METIS_OPTION_NITER,
!!$  METIS_OPTION_NCUTS,
!!$  METIS_OPTION_SEED,
!!$  METIS_OPTION_NO2HOP,
!!$  METIS_OPTION_MINCONN,
!!$  METIS_OPTION_CONTIG,
!!$  METIS_OPTION_COMPRESS,
!!$  METIS_OPTION_CCORDER,
!!$  METIS_OPTION_PFACTOR,
!!$  METIS_OPTION_NSEPS,
!!$  METIS_OPTION_UFACTOR,
!!$  METIS_OPTION_NUMBERING,
!!$
!!$  /* Used for command-line parameter purposes */
!!$  METIS_OPTION_HELP,
!!$  METIS_OPTION_TPWGTS,
!!$  METIS_OPTION_NCOMMON,
!!$  METIS_OPTION_NOOUTPUT,
!!$  METIS_OPTION_BALANCE,
!!$  METIS_OPTION_GTYPE,
!!$  METIS_OPTION_UBVEC
!!$} moptions_et;
!!$
     integer(c_int), parameter :: METIS_OPTION_PTYPE = 0 
     integer(c_int), parameter :: METIS_OPTION_OBJTYPE = 1
     integer(c_int), parameter :: METIS_OPTION_CTYPE = 2
     integer(c_int), parameter :: METIS_OPTION_IPTYPE = 3
     integer(c_int), parameter :: METIS_OPTION_RTYPE = 4
     integer(c_int), parameter :: METIS_OPTION_DBGLVL = 5
     integer(c_int), parameter :: METIS_OPTION_NITER = 6
     integer(c_int), parameter :: METIS_OPTION_NCUTS = 7
     integer(c_int), parameter :: METIS_OPTION_SEED = 8
     integer(c_int), parameter :: METIS_OPTION_NO2HOP = 9
     integer(c_int), parameter :: METIS_OPTION_MINCONN = 10
     integer(c_int), parameter :: METIS_OPTION_CONTIG = 11
     integer(c_int), parameter :: METIS_OPTION_COMPRESS = 12
     integer(c_int), parameter :: METIS_OPTION_CCORDER = 13
     integer(c_int), parameter :: METIS_OPTION_PFACTOR = 14
     integer(c_int), parameter :: METIS_OPTION_NSEPS = 15
     integer(c_int), parameter :: METIS_OPTION_UFACTOR = 16
     integer(c_int), parameter :: METIS_OPTION_NUMBERING = 17
     
     integer(c_int), parameter :: METIS_OPTION_HELP = 18
     integer(c_int), parameter :: METIS_OPTION_TPWGTS = 19
     integer(c_int), parameter :: METIS_OPTION_NCOMMON = 20
     integer(c_int), parameter :: METIS_OPTION_NOOUTPUT = 21
     integer(c_int), parameter :: METIS_OPTION_BALANCE = 22
     integer(c_int), parameter :: METIS_OPTION_GTYPE = 23
     integer(c_int), parameter :: METIS_OPTION_UBVEC = 24

!!$
!!$/*! Partitioning Schemes */
!!$typedef enum {
!!$  METIS_PTYPE_RB, 
!!$  METIS_PTYPE_KWAY                
!!$} mptype_et;
     
     integer(c_int), parameter :: METIS_PTYPE_RB = 0
     integer(c_int), parameter :: METIS_PTYPE_KWAY = 1

!!$
!!$/*! Graph types for meshes */
!!$typedef enum {
!!$  METIS_GTYPE_DUAL,
!!$  METIS_GTYPE_NODAL               
!!$} mgtype_et;
!!$

     integer(c_int), parameter :: METIS_GTYPE_DUAL  = 0
     integer(c_int), parameter :: METIS_GTYPE_NODAL = 1
     

!!$/*! Coarsening Schemes */
!!$typedef enum {
!!$  METIS_CTYPE_RM,
!!$  METIS_CTYPE_SHEM
!!$} mctype_et;
!!$

     integer(c_int), parameter :: METIS_CTYPE_RM = 0
     integer(c_int), parameter :: METIS_CTYPE_SHEM = 1

!!$/*! Initial partitioning schemes */
!!$typedef enum {
!!$  METIS_IPTYPE_GROW,
!!$  METIS_IPTYPE_RANDOM,
!!$  METIS_IPTYPE_EDGE,
!!$  METIS_IPTYPE_NODE,
!!$  METIS_IPTYPE_METISRB
!!$} miptype_et;
!!$
     integer(c_int), parameter :: METIS_IPTYPE_GROW = 0
     integer(c_int), parameter :: METIS_IPTYPE_RANDOM = 1
     integer(c_int), parameter :: METIS_IPTYPE_EDGE = 2
     integer(c_int), parameter :: METIS_IPTYPE_NODE = 3
     integer(c_int), parameter :: METIS_IPTYPE_METISRB = 4

!!$
!!$/*! Refinement schemes */
!!$typedef enum {
!!$  METIS_RTYPE_FM,
!!$  METIS_RTYPE_GREEDY,
!!$  METIS_RTYPE_SEP2SIDED,
!!$  METIS_RTYPE_SEP1SIDED
!!$} mrtype_et;
!!$
     
     integer(c_int), parameter :: METIS_RTYPE_FM = 0
     integer(c_int), parameter :: METIS_RTYPE_GREEDY = 1
     integer(c_int), parameter :: METIS_RTYPE_SEP2SIDED = 2
     integer(c_int), parameter :: METIS_RTYPE_SEP1SIDED = 3 
     
!!$
!!$/*! Debug Levels */
!!$typedef enum {
!!$  METIS_DBG_INFO       = 1,       /*!< Shows various diagnostic messages */
!!$  METIS_DBG_TIME       = 2,       /*!< Perform timing analysis */
!!$  METIS_DBG_COARSEN    = 4,	  /*!< Show the coarsening progress */
!!$  METIS_DBG_REFINE     = 8,	  /*!< Show the refinement progress */
!!$  METIS_DBG_IPART      = 16, 	  /*!< Show info on initial partitioning */
!!$  METIS_DBG_MOVEINFO   = 32, 	  /*!< Show info on vertex moves during refinement */
!!$  METIS_DBG_SEPINFO    = 64, 	  /*!< Show info on vertex moves during sep refinement */
!!$  METIS_DBG_CONNINFO   = 128,     /*!< Show info on minimization of subdomain connectivity */
!!$  METIS_DBG_CONTIGINFO = 256,     /*!< Show info on elimination of connected components */ 
!!$  METIS_DBG_MEMORY     = 2048,    /*!< Show info related to wspace allocation */
!!$} mdbglvl_et;
!!$
  integer(c_int), parameter :: METIS_DBG_INFO       = 1    ! /*!< Shows various diagnostic messages */
  integer(c_int), parameter :: METIS_DBG_TIME       = 2    ! /*!< Perform timing analysis */
  integer(c_int), parameter :: METIS_DBG_COARSEN    = 4	   ! /*!< Show the coarsening progress */
  integer(c_int), parameter :: METIS_DBG_REFINE     = 8	   ! /*!< Show the refinement progress */
  integer(c_int), parameter :: METIS_DBG_IPART      = 16   ! /*!< Show info on initial partitioning */
  integer(c_int), parameter :: METIS_DBG_MOVEINFO   = 32   ! /*!< Show info on vertex moves during refinement */
  integer(c_int), parameter :: METIS_DBG_SEPINFO    = 64   ! /*!< Show info on vertex moves during sep refinement */
  integer(c_int), parameter :: METIS_DBG_CONNINFO   = 128  ! /*!< Show info on minimization of subdomain connectivity */
  integer(c_int), parameter :: METIS_DBG_CONTIGINFO = 256  ! /*!< Show info on elimination of connected components */ 
  integer(c_int), parameter :: METIS_DBG_MEMORY     = 2048 ! /*!< Show info related to wspace allocation */

!!$
!!$/* Types of objectives */
!!$typedef enum {
!!$  METIS_OBJTYPE_CUT,
!!$  METIS_OBJTYPE_VOL,
!!$  METIS_OBJTYPE_NODE
!!$} mobjtype_et;

  integer(c_int), parameter :: METIS_OBJTYPE_CUT =  0 
  integer(c_int), parameter :: METIS_OBJTYPE_VOL =  1
  integer(c_int), parameter :: METIS_OBJTYPE_NODE = 2 

#ifndef METIS_LONG_INTEGERS

  integer(c_int),target :: options(0:METIS_NOPTIONS-1)
  integer(c_int),target :: ncon 
  integer(c_int) :: ierr

  interface
     function fp_metis_nodendextractseparatortree(nvtxs,xadj,adjncy,vwgt,options,perm,iperm,nlevel,graphvert2nodesep) & 
        & bind(c,NAME='FP_METIS_NodeNDExtractSeparatorTree')
       use iso_c_binding
       integer(c_int) :: fp_metis_nodendextractseparatortree
       type(c_ptr), value :: nvtxs, nlevel
       type(c_ptr), value :: xadj, adjncy, vwgt, options, perm, iperm, graphvert2nodesep
     end function fp_metis_nodendextractseparatortree
     function fp_metis_partgraphkway(nvtxs,ncon,xadj,adjncy,vwgt,vsize,adjwgt,nparts,tptwgts,ubvec,options,objval,part) &
        & bind(c,NAME='FP_METIS_PartGraphKway')
       use iso_c_binding
       integer(c_int) :: fp_metis_partgraphkway
       type(c_ptr), value :: nvtxs, ncon, nparts, objval, part
       type(c_ptr), value :: xadj, adjncy, vwgt, vsize, adjwgt, options
       type(c_ptr), value :: tptwgts, ubvec 
       ! WARNING: metis.h, #define REALTYPEWIDTH 64 REQUIRED !!!
     end function fp_metis_partgraphkway
     function fp_metis_partgraphrecursive(nvtxs,ncon,xadj,adjncy,vwgt,vsize,adjwgt,nparts,tptwgts,ubvec,options,objval,part) &
        & bind(c,NAME='FP_METIS_PartGraphRecursive')
       use iso_c_binding
       integer(c_int) :: fp_metis_partgraphrecursive
       type(c_ptr), value :: nvtxs, ncon, nparts, objval, part
       type(c_ptr), value :: xadj, adjncy, vwgt, vsize, adjwgt, options
       type(c_ptr), value :: tptwgts, ubvec 
       ! WARNING: metis.h, #define REALTYPEWIDTH 64 REQUIRED !!!
     end function fp_metis_partgraphrecursive
     function fp_metis_setdefaultoptions(options) bind(c,NAME='FP_METIS_SetDefaultOptions')
       use iso_c_binding
       integer(c_int)      :: fp_metis_setdefaultoptions
       type(c_ptr), value  :: options
     end function fp_metis_setdefaultoptions
  end interface

#else
  ! TODO (done by Alberto in his private 64-bit version of Fempar).
#endif

#else
#ifdef ENABLE_METIS4

#ifndef METIS_LONG_INTEGERS
  ! Metis-4.0 default options except 6 (5 in metis C style).
  integer(c_int)                :: optmt_nd(8)=(/1,3,1,2,0,0,0,1/)
  !integer(ip)                :: optmt_nd(8)=(/1,3,1,2,1+2+4+8+16+32+64+128+256+512,0,0,1/)

  ! Metis-4-0 default options
  integer(c_int)                :: optmt_pr(5)=(/0,0,0,0,0/)

  interface
     subroutine fp_metis_nodendextractseparatortree(n,ia,ja,m,no,lo,li,nl,gv2ns) & 
        & bind(c,NAME='fp_metis_nodendextractseparatortree')
       use iso_c_binding
       integer(c_int) :: n, m, no(8), nl
       integer(c_int) :: ia(*), ja(*), lo(*), li(*), gv2ns(*)
     end subroutine fp_metis_nodendextractseparatortree
     subroutine fp_metis_nodend(n,ia,ja,m,no,lo,li) bind(c,NAME='fp_metis_nodend')
       use iso_c_binding
       integer(c_int) :: n,m,no(8)
       integer(c_int) :: ia(*),ja(*),lo(*),li(*)
     end subroutine fp_metis_nodend
     subroutine fp_metis_nodewnd(n,ia,ja,iw,m,no,lo,li) bind(c,NAME='fp_metis_nodewnd')
       use iso_c_binding
       integer(c_int) :: n,m,no(8)
       integer(c_int) :: ia(*),ja(*),iw(*),lo(*),li(*)
     end subroutine fp_metis_nodewnd
     subroutine fp_metis_partgraphkway(n,ie,je,i1,i2,kw,nf,np,no,nc,led) &
        & bind(c,NAME='fp_metis_partgraphkway')
       use iso_c_binding
       integer(c_int) :: n,i1,i2,kw,nf,np,no(5),nc
       integer(c_int) :: led(*),ie(*),je(*)
     end subroutine fp_metis_partgraphkway
     subroutine fp_metis_partgraphrecursive(n,ie,je,i1,i2,kw,nf,np,no,nc,led) &
        & bind(c,NAME='fp_metis_partgraphrecursive')
       use iso_c_binding
       integer(c_int) :: n,i1,i2,kw,nf,np,no(5),nc
       integer(c_int) :: led(*),ie(*),je(*)
     end subroutine 
  end interface

#else

  ! Metis-4.0 default options except 6 (5 in metis C style).
  integer(c_long)                :: optmt_nd(8)=(/1,3,1,2,0,0,0,1/)
  !integer(ip)                :: optmt_nd(8)=(/1,3,1,2,1+2+4+8+16+32+64+128+256+512,0,0,1/)

  ! Metis-4-0 default options
  integer(c_long)                :: optmt_pr(5)=(/1,3,1,3,0/)

     function fp_metis_nodendextractseparatortree(n,ia,ja,m,no,lo,li,nl,gv2ns) & 
        & bind(c,NAME='fp_metis_nodendextractseparatortree')
       use iso_c_binding
       integer(c_long)  :: n, m, no(8), nl
       integer(c_long) :: ia(*), ja(*), lo(*), li(*), gv2ns(*)
     end subroutine fp_metis_nodendextractseparatortree
     subroutine fp_metis_nodend(n,ia,ja,m,no,lo,li) bind(c,NAME='fp_metis_nodend')
       use iso_c_binding
       integer(c_long)  :: n, m, no(8)
       integer(c_long) :: ia(*), ja(*), lo(*), li(*)
     end subroutine fp_metis_nodend
     subroutine fp_metis_nodewnd(n,ia,ja,iw,m,no,lo,li) bind(c,NAME='fp_metis_nodewnd')
       use iso_c_binding
       integer(c_long)  :: n,m,no(8)
       integer(c_long) :: ia(*),ja(*),iw(*),lo(*),li(*)
     end subroutine fp_metis_nodewnd
     subroutine fp_metis_partgraphkway(n,ie,je,i1,i2,kw,nf,np,no,nc,led) &
        & bind(c,NAME='fp_metis_partgraphkway')
       use iso_c_binding
       integer(c_long)  :: n,kw,nf,np,no(5),nc
       integer(c_long) :: ie(*),je(*),i1,i2,led(*)
     end subroutine fp_metis_partgraphkway
     subroutine fp_metis_partgraphrecursive(n,ie,je,i1,i2,kw,nf,np,no,nc,led) &
        & bind(c,NAME='fp_metis_partgraphrecursive')
       use iso_c_binding
       integer(c_long)  :: n,kw,nf,np,no(5),nc
       integer(c_long) :: ie(*),je(*),i1,i2,led(*)
     end subroutine fp_metis_partgraphrecursive

#endif ! Long integers

#else 

contains 
  ! Public by default
  subroutine enable_metis_error_message
    implicit none
    write (0,*) 'Error: FemPar was not compiled with -DENABLE_METIS5 nor  -DENABLE_METIS4.'
    write (0,*) "Error: You must activate any of these cpp macro in order to use Fempar's interface to Metis"
    call runend
  end subroutine enable_metis_error_message

#endif ! metis4

#endif ! metis5

end module metis_interface
