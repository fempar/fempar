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
module par_context_names
  ! Serial modules
  use types_names

  ! Parallel modules
  use psb_penv_mod_names

#ifdef MPI_MOD
  use mpi
#endif
  implicit none 
#ifdef MPI_H
  include 'mpif.h'
#endif
    
#include "debug.i90"
  private

  ! Parallel context
  type par_context_t
     private 
     ! The following member is an integer which identifies an initialized 
     ! parallel context which is created and manipulated using 
     ! our own wrappers to the MPI library (adapted from PSBLAS 2.4.0 library)

     ! ***IMPORTANT NOTE***: parallel contexts are always of type 
     ! integer: the kind parameter must NOT be specified. This requirement is 
     ! imposed by the underlying message-passing library, i.e., MPI. 
     integer :: icontxt = mpi_comm_null
     integer :: rank    = -1
     integer :: size    = -1
  contains
     procedure, non_overridable, private :: par_context_create_by_comm_world_dup
     procedure, non_overridable, private :: par_context_create_from_base_context_and_task_ids
     procedure, non_overridable, private :: par_context_create_intercomm
     generic                             :: create => par_context_create_by_comm_world_dup, &
                                                      par_context_create_from_base_context_and_task_ids, &
                                                      par_context_create_intercomm
     
     procedure, non_overridable, private :: par_context_split
     procedure, non_overridable, private :: par_context_split_into_partition_with_two_subcontexts
     generic                             :: split => par_context_split, &
                                                     par_context_split_into_partition_with_two_subcontexts
                                                      
     procedure, non_overridable :: free        => par_context_free
     procedure, non_overridable :: nullify     => par_context_nullify
     procedure, non_overridable :: info        => par_context_info
     procedure, non_overridable :: get_icontxt => par_context_get_icontxt
     procedure, non_overridable :: get_rank    => par_context_get_rank
     procedure, non_overridable :: get_size    => par_context_get_size
  end type par_context_t

  ! Types
  public :: par_context_t

contains

  subroutine par_context_create_by_comm_world_dup ( this )
    implicit none 
    class(par_context_t), intent(inout) :: this
    call this%free(finalize=.false.)
    call psb_init ( this%icontxt )
    call psb_info ( this%icontxt, this%rank, this%size )
  end subroutine par_context_create_by_comm_world_dup
  
  subroutine par_context_create_from_base_context_and_task_ids ( this, base_context, ntasks, task_ids )
    implicit none 
    class(par_context_t)          , intent(inout) :: this
    type(par_context_t)           , intent(in)    :: base_context
    integer(ip)                   , intent(in)    :: ntasks
    integer(ip)         , optional, intent(in)    :: task_ids(ntasks)  
    call this%free(finalize=.false.)
    call psb_init ( this%icontxt, ntasks, base_context%icontxt, task_ids)
    call psb_info ( this%icontxt, this%rank, this%size )
  end subroutine par_context_create_from_base_context_and_task_ids
  
  subroutine par_context_create_intercomm ( this, world_context, sub_context_1, sub_context_2)
    implicit none 
    class(par_context_t)          , intent(inout) :: this
    type(par_context_t)           , intent(in)    :: world_context
    type(par_context_t)           , intent(in)    :: sub_context_1
    type(par_context_t)           , intent(in)    :: sub_context_2
    integer                                       :: info
    
    call this%free(finalize=.false.)
    
    ! {{sub_context_1} U {sub_context_2}} MUST be a partition of world_context (precondition). 
    ! Further we assume that the ranks ids' (in world_context) of tasks in sub_context_1 go from 
    ! 0 to size(sub_context_1)-1 and ranks ids' of tasks in sub_context_2 from that number onwards
    assert (world_context%rank >= 0)

    ! Create intracomm
    if(sub_context_1%rank>=0) then
       assert(sub_context_2%rank<0)
       !    mpi_intercomm_create (local_comm, local_leader, peer_comm, remote_leader, tag, newintercomm, info)
       call mpi_intercomm_create(sub_context_1%icontxt, 0, world_context%icontxt, sub_context_1%size, 0, this%icontxt, info)
    else 
       assert(sub_context_2%rank>=0)
       call mpi_intercomm_create(sub_context_2%icontxt, 0, world_context%icontxt, 0 , 0, this%icontxt, info)
    end if
    assert(info==mpi_success)
    
    ! and store current info (who am I and how many we are in these contexts).
    call psb_info ( this%icontxt, this%rank, this%size )
  end subroutine par_context_create_intercomm
  
  subroutine par_context_split ( this, color, new_subcontext )
    implicit none 
    class(par_context_t)          , intent(in)    :: this
    integer                       , intent(in)    :: color
    type(par_context_t)           , intent(inout) :: new_subcontext
    integer                       , parameter     :: key=0
    integer                                       :: info

    call new_subcontext%free(finalize=.false.)
    call mpi_comm_split(this%icontxt, color, key, new_subcontext%icontxt, info)
    assert ( info == mpi_success )
    call psb_info ( new_subcontext%icontxt, new_subcontext%rank, new_subcontext%size )
  end subroutine par_context_split
  
  subroutine par_context_split_into_partition_with_two_subcontexts ( this, in_subcontext1, subcontext1, subcontext2 )
    implicit none 
    class(par_context_t)          , intent(in)    :: this
    logical                       , intent(in)    :: in_subcontext1
    type(par_context_t)           , intent(inout) :: subcontext1
    type(par_context_t)           , intent(inout) :: subcontext2
    integer                                       :: info
    integer                       , parameter     :: key=0
    
    call subcontext1%free(finalize=.false.)
    call subcontext2%free(finalize=.false.)
    
    if ( in_subcontext1 ) then
       call mpi_comm_split(this%icontxt, 1, key, subcontext1%icontxt, info)
       call subcontext2%nullify()
    else
       call mpi_comm_split(this%icontxt, 2, key, subcontext2%icontxt, info)
       call subcontext1%nullify()
    end if
    assert ( info == mpi_success )
    call psb_info ( subcontext1%icontxt, subcontext1%rank, subcontext1%size )
    call psb_info ( subcontext2%icontxt, subcontext2%rank, subcontext2%size )
  end subroutine par_context_split_into_partition_with_two_subcontexts
  
  ! Frees the memory related to the communication object underlying "this"
  ! and finalizes the underlying message-passing library (i.e., MPI) if 
  ! finalize == .true. 
  subroutine par_context_free ( p_context, finalize  )
    implicit none 
    ! Parameters
    class(par_context_t)            , intent(inout) :: p_context
    logical                         , intent(in)    :: finalize
    call psb_exit ( p_context%icontxt, finalize )
  end subroutine par_context_free
  
  subroutine par_context_nullify ( this )
    implicit none 
    class(par_context_t), intent(inout) :: this
    call this%free(finalize=.false.)
    this%icontxt = mpi_comm_null
    call psb_info ( this%icontxt, this%rank, this%size )
  end subroutine par_context_nullify

  subroutine par_context_info ( p_context, rank, size )
    ! Parameters
    class(par_context_t), intent(in)  :: p_context
    integer             , intent(out) :: rank
    integer             , intent(out) :: size
    rank = p_context%rank
    size = p_context%size
  end subroutine par_context_info
  
  function par_context_get_icontxt (this)
    implicit none
    class(par_context_t), intent(in) :: this
    integer :: par_context_get_icontxt
    par_context_get_icontxt = this%icontxt
  end function par_context_get_icontxt
  
  function par_context_get_rank (this)
    implicit none
    class(par_context_t), intent(in) :: this
    integer :: par_context_get_rank
    par_context_get_rank = this%rank
  end function par_context_get_rank
  
  function par_context_get_size (this)
    implicit none
    class(par_context_t), intent(in) :: this
    integer :: par_context_get_size
    par_context_get_size = this%size
  end function par_context_get_size

  subroutine par_context_create_by_color ( my_color, p_context, q_context, b_context )
    implicit none 

    ! Parameters
    integer(ip)        , intent(in)    :: my_color
    type(par_context_t), intent(inout) :: p_context
    type(par_context_t), intent(inout) :: q_context
    type(par_context_t), intent(in)    :: b_context 

    ! ***IMPORTANT NOTE***: for the following two variables, 
    ! the kind parameter must NOT be specified. This requirement is 
    ! imposed by the underlying message-passing library, i.e., MPI 
    integer                         :: info
    integer             , parameter :: key=0
    logical                         :: initialized

    call mpi_comm_split(b_context%icontxt, my_color, key , p_context%icontxt, info)
    q_context%icontxt = mpi_comm_null
    assert(info==mpi_success)

    ! and store current info (who am I and how many we are in these contexts).
    call psb_info ( p_context%icontxt, p_context%rank, p_context%size )
    call psb_info ( q_context%icontxt, q_context%rank, q_context%size )
  end subroutine par_context_create_by_color
  
end module par_context_names
