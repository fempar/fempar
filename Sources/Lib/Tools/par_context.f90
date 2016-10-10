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
  use memor_names

  ! Parallel modules
  use psb_const_mod_names
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

  integer(ip), parameter :: undefined_color=-1
  public :: undefined_color

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
          &                                           par_context_create_from_base_context_and_task_ids, &
          &                                           par_context_create_intercomm

     procedure, non_overridable, private :: par_context_assign
     generic                             :: assignment(=)  => par_context_assign

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
     procedure, non_overridable :: am_i_member => par_context_am_i_member

     procedure, non_overridable :: barrier     => par_context_barrier

     procedure, non_overridable :: sum_scalar_rp            => par_context_sum_scalar_rp
     procedure, non_overridable :: sum_vector_rp            => par_context_sum_vector_rp
     procedure, non_overridable :: max_scalar_rp            => par_context_max_scalar_rp
     procedure, non_overridable :: max_vector_rp            => par_context_max_vector_rp


     procedure :: scatter => par_context_scatter_scalar_ip
     procedure :: gather  => par_context_gather_scalar_ip
     procedure :: bcast   => par_context_bcast_scalar_ip

     procedure :: bcast_subcontext   => par_context_bcast_subcontext


     procedure, private :: par_context_neighbours_exchange_rp
     procedure, private :: par_context_neighbours_exchange_ip
     procedure, private :: par_context_neighbours_exchange_igp
     procedure, private :: par_context_neighbours_exchange_single_ip
     procedure, private :: par_context_neighbours_exchange_wo_pack_unpack_ieep
     generic   :: neighbours_exchange      => par_context_neighbours_exchange_rp, &
          par_context_neighbours_exchange_ip,&
          par_context_neighbours_exchange_igp,&
          par_context_neighbours_exchange_single_ip, &
          par_context_neighbours_exchange_wo_pack_unpack_ieep

     procedure, private :: transfer_to_master_ip          => par_context_transfer_to_master_ip
     procedure, private :: transfer_to_master_ip_1D_array => par_context_transfer_to_master_ip_1D_array
     generic :: transfer_to_master => transfer_to_master_ip, transfer_to_master_ip_1D_array

     procedure, private :: par_context_gather_to_master_ip
     procedure, private :: par_context_gather_to_master_igp
     procedure, private :: par_context_gather_to_master_ip_1D_array
     procedure, private :: par_context_gather_to_masterv_ip_1D_array
     procedure, private :: par_context_gather_to_masterv_igp_1D_array
     procedure, private :: par_context_gather_to_masterv_rp_1D_array
     procedure, private :: par_context_gather_to_masterv_rp_2D_array
     generic   :: gather_to_master => par_context_gather_to_master_ip, &
          par_context_gather_to_master_igp, &
          par_context_gather_to_master_ip_1D_array, &
          par_context_gather_to_masterv_ip_1D_array, &
          par_context_gather_to_masterv_igp_1D_array, &
          par_context_gather_to_masterv_rp_1D_array, &
          par_context_gather_to_masterv_rp_2D_array

     procedure, private :: par_context_scatter_from_masterv_rp_1D_array                                  
     generic   :: scatter_from_master => par_context_scatter_from_masterv_rp_1D_array

  end type par_context_t

  ! Types
  public :: par_context_t

contains

  !=============================================================================
  subroutine par_context_assign(this, that)
    implicit none 
    class(par_context_t), intent(inout) :: this
    class(par_context_t), intent(in)    :: that
    call this%free(finalize=.false.)
    this%icontxt = that%icontxt
    this%rank = that%rank
    this%size = that%size
  end subroutine par_context_assign

  !=============================================================================
  subroutine par_context_create_by_comm_world_dup ( this )
    implicit none 
    class(par_context_t), intent(inout) :: this
    call this%free(finalize=.false.)
    call psb_init ( this%icontxt )
    call psb_info ( this%icontxt, this%rank, this%size )
  end subroutine par_context_create_by_comm_world_dup

  !=============================================================================
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

  !=============================================================================
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

  !=============================================================================
  subroutine par_context_split ( this, color, new_subcontext )
    implicit none 
    class(par_context_t), intent(in)    :: this
    integer             , intent(in)    :: color
    type(par_context_t) , intent(inout) :: new_subcontext
    integer             , parameter     :: key=0
    integer  :: info,my_color

    if(color==undefined_color) then
       my_color = mpi_undefined
    else
       my_color = color
    end if

    call new_subcontext%free(finalize=.false.)
    call mpi_comm_split(this%icontxt, my_color, key, new_subcontext%icontxt, info)
    assert ( info == mpi_success )
    call psb_info ( new_subcontext%icontxt, new_subcontext%rank, new_subcontext%size )
  end subroutine par_context_split

  !=============================================================================
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

  !=============================================================================
  ! Frees the memory related to the communication object underlying "this"
  ! and finalizes the underlying message-passing library (i.e., MPI) if 
  ! finalize == .true. 
  subroutine par_context_free ( p_context, finalize  )
    implicit none 
    ! Parameters
    class(par_context_t)            , intent(inout) :: p_context
    logical                         , intent(in)    :: finalize

    if(p_context%icontxt/=mpi_comm_null) call psb_exit ( p_context%icontxt, finalize )

    p_context%icontxt=mpi_comm_null
    p_context%rank=-1
    p_context%size=-1

  end subroutine par_context_free

  !=============================================================================
  subroutine par_context_nullify ( this )
    implicit none 
    class(par_context_t), intent(inout) :: this
    call this%free(finalize=.false.)
    this%icontxt = mpi_comm_null
    call psb_info ( this%icontxt, this%rank, this%size )
  end subroutine par_context_nullify

  !=============================================================================
  subroutine par_context_info ( p_context, rank, size )
    ! Parameters
    class(par_context_t), intent(in)  :: p_context
    integer             , intent(out) :: rank
    integer             , intent(out) :: size
    rank = p_context%rank
    size = p_context%size
  end subroutine par_context_info

  !=============================================================================
  pure function par_context_get_icontxt (this)
    implicit none
    class(par_context_t), intent(in) :: this
    integer :: par_context_get_icontxt
    par_context_get_icontxt = this%icontxt
  end function par_context_get_icontxt

  !=============================================================================
  pure function par_context_get_rank (this)
    implicit none
    class(par_context_t), intent(in) :: this
    integer :: par_context_get_rank
    par_context_get_rank = this%rank
  end function par_context_get_rank

  !=============================================================================
  pure function par_context_get_size (this)
    implicit none
    class(par_context_t), intent(in) :: this
    integer :: par_context_get_size
    par_context_get_size = this%size
  end function par_context_get_size

  !=============================================================================
  pure function par_context_am_i_member(this)
    implicit none
    class(par_context_t), intent(in) :: this
    logical                          :: par_context_am_i_member
    par_context_am_i_member = (this%rank>=0)
  end function par_context_am_i_member

  !=============================================================================
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

  !=============================================================================
  subroutine par_context_barrier(this)
    implicit none 
    class(par_context_t), intent(in) :: this
    integer :: mpi_comm_p, ierr
    call psb_get_mpicomm (this%icontxt, mpi_comm_p)
    call mpi_barrier ( mpi_comm_p, ierr)
    check ( ierr == mpi_success )
  end subroutine par_context_barrier

  !=============================================================================
  subroutine par_context_sum_scalar_rp (this,alpha)
    implicit none
    class(par_context_t) , intent(in)    :: this
    real(rp)             , intent(inout) :: alpha
    call psb_sum(this%icontxt, alpha)
  end subroutine par_context_sum_scalar_rp

  !=============================================================================
  subroutine par_context_sum_vector_rp(this,alpha)
    implicit none
    class(par_context_t) , intent(in)    :: this
    real(rp)             , intent(inout) :: alpha(:) 
    call psb_sum(this%icontxt, alpha)
  end subroutine par_context_sum_vector_rp

  !=============================================================================
  subroutine par_context_max_scalar_rp (this,alpha)
    implicit none
    class(par_context_t) , intent(in)    :: this
    real(rp)             , intent(inout) :: alpha
    call psb_max(this%icontxt, alpha)
  end subroutine par_context_max_scalar_rp

  !=============================================================================
  subroutine par_context_max_vector_rp(this,alpha)
    implicit none
    class(par_context_t) , intent(in)    :: this
    real(rp)             , intent(inout) :: alpha(:) 
    call psb_max(this%icontxt, alpha)
  end subroutine par_context_max_vector_rp

  !=============================================================================
  subroutine par_context_bcast_subcontext(this,subcontxt,condition)
    implicit none
    class(par_context_t) , intent(in)    :: this
    class(par_context_t) , intent(in)    :: subcontxt
    logical              , intent(inout) :: condition

    integer(ip) :: recv_rank
    integer(ip) :: send_rank
    integer :: mpi_comm, info

    send_rank = psb_root_
    recv_rank = this%size

    if(this%icontxt==send_rank) then
       call psb_snd(this%icontxt, condition, recv_rank)
    else if(this%icontxt==recv_rank) then
       call psb_rcv(this%icontxt, condition, send_rank)
    end if
    if(subcontxt%am_i_member()) then
       call psb_get_mpicomm (subcontxt%icontxt, mpi_comm)
       call mpi_bcast(condition,1,MPI_LOGICAL,psb_root_,mpi_comm,info)
       check( info == mpi_success )
    end if

  end subroutine par_context_bcast_subcontext

  !=============================================================================
  ! When packing   (gathering) ,    buffer <- alpha * x
  ! When unpacking (scattering),    x <- beta*x + buffer
  subroutine par_context_neighbours_exchange_rp ( this, & 
       &                                          num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       &                                          num_snd, list_snd, snd_ptrs, pack_idx,   &
       &                                          alpha, beta, x)
    ! use psb_const_mod_names
    ! use psb_penv_mod_names
    implicit none
    class(par_context_t), intent(in) :: this

    ! Control info to receive
    integer(ip)             , intent(in) :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(in) :: unpack_idx (rcv_ptrs(num_rcv+1)-1)

    ! Control info to send
    integer(ip)             , intent(in) :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in) :: pack_idx (snd_ptrs(num_snd+1)-1)

    ! Floating point data
    real(rp), intent(in)    :: alpha, beta
    real(rp), intent(inout) :: x(:)

    ! Communication related locals 
    integer :: my_pid, i, proc_to_comm, sizmsg
    integer :: the_mpi_comm,  iret
    integer :: p2pstat(mpi_status_size)
    integer :: icontxt

    ! Request handlers for non-blocking receives
    integer, allocatable :: rcvhd(:)

    ! Request handlers for non-blocking receives
    integer, allocatable :: sndhd(:)

    real(rp), allocatable :: sndbuf(:) 
    real(rp), allocatable :: rcvbuf(:)

    call psb_get_mpicomm (this%icontxt, the_mpi_comm)

    call memalloc (num_rcv, rcvhd, __FILE__,__LINE__)
    call memalloc (num_snd, sndhd, __FILE__,__LINE__)

    call memalloc ((snd_ptrs(num_snd+1)-snd_ptrs(1)), sndbuf, __FILE__,__LINE__)
    call memalloc ((rcv_ptrs(num_rcv+1)-rcv_ptrs(1)), rcvbuf, __FILE__,__LINE__)

    ! Pack send buffers
    call pack_rp ( snd_ptrs(num_snd+1)-snd_ptrs(1), pack_idx, alpha, x, sndbuf )

    ! First post all the non blocking receives   
    do i=1, num_rcv
       proc_to_comm = list_rcv(i)

       ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
       call psb_get_rank (proc_to_comm, this%icontxt, proc_to_comm-1)

       ! Message size to be received
       sizmsg = rcv_ptrs(i+1)-rcv_ptrs(i)

       if ( (sizmsg > 0) .and. (list_rcv(i)-1 /= this%rank) ) then
          call mpi_irecv(  rcvbuf(rcv_ptrs(i)), sizmsg,        &
               &  psb_mpi_real, proc_to_comm, &
               &  psb_double_swap_tag, the_mpi_comm, rcvhd(i), iret)
          check ( iret == mpi_success )
       end if
    end do

    ! Secondly post all non-blocking sends
    do i=1, num_snd
       proc_to_comm = list_snd(i)

       ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
       call psb_get_rank (proc_to_comm, this%icontxt, proc_to_comm-1)

       ! Message size to be sent
       sizmsg = snd_ptrs(i+1)-snd_ptrs(i)

       if ( (sizmsg > 0) .and. (list_snd(i)-1 /= this%rank) ) then 
          call mpi_isend(sndbuf(snd_ptrs(i)), sizmsg, &
               & psb_mpi_real, proc_to_comm,    &
               & psb_double_swap_tag, the_mpi_comm, sndhd(i), iret)
          check ( iret == mpi_success )
       end if
    end do

    ! Wait on all non-blocking receives
    do i=1, num_rcv
       proc_to_comm = list_rcv(i)

       ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
       call psb_get_rank (proc_to_comm, this%icontxt, proc_to_comm-1)

       ! Message size to be received
       sizmsg = rcv_ptrs(i+1)-rcv_ptrs(i)

       if ( (sizmsg > 0) .and. (list_rcv(i)-1 /= this%rank) ) then
          call mpi_wait(rcvhd(i), p2pstat, iret)

       else if ( list_rcv(i)-1 == this%rank ) then
          if ( sizmsg /= snd_ptrs(i+1)-snd_ptrs(i) ) then 
             write(0,*) 'Fatal error in single_exchange: mismatch on self sendf', & 
                  & sizmsg, snd_ptrs(i+1)-snd_ptrs(i) 
          end if

          rcvbuf( rcv_ptrs(i):rcv_ptrs(i)+sizmsg-1) = &
               sndbuf( snd_ptrs(i): snd_ptrs(i)+sizmsg-1 )
       end if
    end do

    ! Finally wait on all non-blocking sends
    do i=1, num_snd
       proc_to_comm = list_snd(i)

       ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
       call psb_get_rank (proc_to_comm, this%icontxt, proc_to_comm-1)

       ! Message size to be received
       sizmsg = snd_ptrs(i+1)-snd_ptrs(i)

       if ( (sizmsg > 0) .and. (list_snd(i)-1 /= this%rank) ) then
          call mpi_wait(sndhd(i), p2pstat, iret)
          check ( iret == mpi_success )
       end if
    end do

    ! Unpack recv buffers
    call unpack_rp (rcv_ptrs(num_rcv+1)-rcv_ptrs(1), unpack_idx, beta, rcvbuf, x )

    call memfree (rcvhd,__FILE__,__LINE__) 
    call memfree (sndhd,__FILE__,__LINE__)

    call memfree (sndbuf,__FILE__,__LINE__)
    call memfree (rcvbuf,__FILE__,__LINE__)

  end subroutine par_context_neighbours_exchange_rp

  !=============================================================================
  ! When packing   (gathering) ,    buffer <- alpha * x
  ! When unpacking (scattering),    x <- beta*x + buffer
  subroutine par_context_neighbours_exchange_ip ( this, & 
       &                                          num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       &                                          num_snd, list_snd, snd_ptrs, pack_idx,   &
       &                                          x,chunk_size)
    ! use psb_const_mod_names
    ! use psb_penv_mod_names
    implicit none
    class(par_context_t), intent(in)    :: this
    ! Control info to receive
    integer(ip)             , intent(in)    :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(in)    :: unpack_idx (rcv_ptrs(num_rcv+1)-1)
    ! Control info to send
    integer(ip)             , intent(in)    :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in)    :: pack_idx (snd_ptrs(num_snd+1)-1)
    ! Raw data to be exchanged
    integer(ip)             , intent(inout) :: x(:)
    integer(ip)   , optional, intent(in)    :: chunk_size

    ! Communication related locals 
    integer :: my_pid, i, proc_to_comm, sizmsg
    integer :: the_mpi_comm,  iret
    integer :: p2pstat(mpi_status_size)
    integer :: icontxt

    ! Request handlers for non-blocking receives
    integer, allocatable :: rcvhd(:)

    ! Request handlers for non-blocking receives
    integer, allocatable :: sndhd(:)

    integer(ip), allocatable :: sndbuf(:) 
    integer(ip), allocatable :: rcvbuf(:)

    integer(ip) :: chunk_size_


    if ( present(chunk_size) ) then
       chunk_size_ = chunk_size
    else
       chunk_size_ = 1
    end if

    call psb_get_mpicomm (this%icontxt, the_mpi_comm)

    call memalloc (num_rcv, rcvhd, __FILE__,__LINE__)
    call memalloc (num_snd, sndhd, __FILE__,__LINE__)

    call memalloc ((snd_ptrs(num_snd+1)-snd_ptrs(1))*chunk_size_, sndbuf, __FILE__,__LINE__)
    call memalloc ((rcv_ptrs(num_rcv+1)-rcv_ptrs(1))*chunk_size_, rcvbuf, __FILE__,__LINE__)

    ! Pack send buffers
    call pack_ip ( snd_ptrs(num_snd+1)-snd_ptrs(1), chunk_size_, pack_idx, x, sndbuf )

    ! First post all the non blocking receives   
    do i=1, num_rcv
       proc_to_comm = list_rcv(i)

       ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
       call psb_get_rank (proc_to_comm, this%icontxt, proc_to_comm-1)

       ! Message size to be received
       sizmsg = (rcv_ptrs(i+1)-rcv_ptrs(i))*chunk_size_

       if ( (sizmsg > 0) .and. (list_rcv(i)-1 /= this%rank) ) then
          call mpi_irecv(  rcvbuf((rcv_ptrs(i)-1)*chunk_size_+1), sizmsg,        &
               &  psb_mpi_integer, proc_to_comm, &
               &  psb_double_swap_tag, the_mpi_comm, rcvhd(i), iret)
          check ( iret == mpi_success )
       end if
    end do

    ! Secondly post all non-blocking sends
    do i=1, num_snd
       proc_to_comm = list_snd(i)

       ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
       call psb_get_rank (proc_to_comm, this%icontxt, proc_to_comm-1)

       ! Message size to be sent
       sizmsg = (snd_ptrs(i+1)-snd_ptrs(i))*chunk_size_

       if ( (sizmsg > 0) .and. (list_snd(i)-1 /= this%rank) ) then 
          call mpi_isend(sndbuf((snd_ptrs(i)-1)*chunk_size_+1), sizmsg, &
               & psb_mpi_integer, proc_to_comm,    &
               & psb_double_swap_tag, the_mpi_comm, sndhd(i), iret)
          check ( iret == mpi_success )
       end if
    end do

    ! Wait on all non-blocking receives
    do i=1, num_rcv
       proc_to_comm = list_rcv(i)

       ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
       call psb_get_rank (proc_to_comm, this%icontxt, proc_to_comm-1)

       ! Message size to be received
       sizmsg = (rcv_ptrs(i+1)-rcv_ptrs(i))*chunk_size_

       if ( (sizmsg > 0) .and. (list_rcv(i)-1 /= this%rank) ) then
          call mpi_wait(rcvhd(i), p2pstat, iret)
          check ( iret == mpi_success )
       else if ( list_rcv(i)-1 == this%rank ) then
          if ( sizmsg /= (snd_ptrs(i+1)-snd_ptrs(i))*chunk_size_ ) then 
             write(0,*) 'Fatal error in single_exchange: mismatch on self sendf', & 
                  & sizmsg, snd_ptrs(i+1)-snd_ptrs(i) 
             check(.false.)
          end if
          rcvbuf((rcv_ptrs(i)-1)*chunk_size_+1:(rcv_ptrs(i)-1)*chunk_size_+sizmsg) = &
               sndbuf( (snd_ptrs(i)-1)*chunk_size_+1:(snd_ptrs(i)-1)*chunk_size_+sizmsg )
       end if
    end do

    ! Finally wait on all non-blocking sends
    do i=1, num_snd
       proc_to_comm = list_snd(i)

       ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
       call psb_get_rank (proc_to_comm, this%icontxt, proc_to_comm-1)

       ! Message size to be received
       sizmsg = (snd_ptrs(i+1)-snd_ptrs(i))*chunk_size_

       if ( (sizmsg > 0) .and. (list_snd(i)-1 /= this%rank) ) then
          call mpi_wait(sndhd(i), p2pstat, iret)
          check ( iret == mpi_success )
       end if
    end do

    ! Unpack recv buffers
    call unpack_ip (rcv_ptrs(num_rcv+1)-rcv_ptrs(1), chunk_size_, unpack_idx, rcvbuf, x )

    call memfree (rcvhd,__FILE__,__LINE__) 
    call memfree (sndhd,__FILE__,__LINE__)

    call memfree (sndbuf,__FILE__,__LINE__)
    call memfree (rcvbuf,__FILE__,__LINE__)
  end subroutine par_context_neighbours_exchange_ip

  !=============================================================================
  subroutine par_context_neighbours_exchange_igp ( this, & 
       &                                              num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       &                                              num_snd, list_snd, snd_ptrs, pack_idx,   &
       &                                              x, chunk_size)
    !use psb_const_mod_names
    !use psb_penv_mod_names
    implicit none
    class(par_context_t), intent(in)    :: this
    ! Control info to receive
    integer(ip)             , intent(in)    :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(in)    :: unpack_idx (rcv_ptrs(num_rcv+1)-1)
    ! Control info to send
    integer(ip)             , intent(in)    :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in)    :: pack_idx (snd_ptrs(num_snd+1)-1)
    ! Raw data to be exchanged
    integer(igp)            , intent(inout) :: x(:)
    integer(ip)   , optional, intent(in)    :: chunk_size

    ! Communication related locals 
    integer :: my_pid, i, proc_to_comm, sizmsg
    integer :: the_mpi_comm,  iret
    integer :: p2pstat(mpi_status_size)
    integer :: icontxt

    ! Request handlers for non-blocking receives
    integer, allocatable :: rcvhd(:)

    ! Request handlers for non-blocking receives
    integer, allocatable :: sndhd(:)

    integer(igp), allocatable :: sndbuf(:) 
    integer(igp), allocatable :: rcvbuf(:)

    integer(ip) :: chunk_size_

    if ( present(chunk_size) ) then
       chunk_size_ = chunk_size
    else
       chunk_size_ = 1
    end if

    call psb_get_mpicomm (this%icontxt, the_mpi_comm)

    call memalloc (num_rcv, rcvhd, __FILE__,__LINE__)
    call memalloc (num_snd, sndhd, __FILE__,__LINE__)

    call memalloc ((snd_ptrs(num_snd+1)-snd_ptrs(1))*chunk_size_, sndbuf, __FILE__,__LINE__)
    call memalloc ((rcv_ptrs(num_rcv+1)-rcv_ptrs(1))*chunk_size_, rcvbuf, __FILE__,__LINE__)

    ! Pack send buffers
    call pack_igp ( snd_ptrs(num_snd+1)-snd_ptrs(1), chunk_size_, pack_idx, x, sndbuf )

    ! First post all the non blocking receives   
    do i=1, num_rcv
       proc_to_comm = list_rcv(i)

       ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
       call psb_get_rank (proc_to_comm, this%icontxt, proc_to_comm-1)

       ! Message size to be received
       sizmsg = (rcv_ptrs(i+1)-rcv_ptrs(i))*chunk_size_

       if ( (sizmsg > 0) .and. (list_rcv(i)-1 /= this%rank) ) then
          call mpi_irecv(  rcvbuf((rcv_ptrs(i)-1)*chunk_size_+1), sizmsg,        &
               &  psb_mpi_long_integer, proc_to_comm, &
               &  psb_double_swap_tag, the_mpi_comm, rcvhd(i), iret)
          check ( iret == mpi_success )
       end if
    end do

    ! Secondly post all non-blocking sends
    do i=1, num_snd
       proc_to_comm = list_snd(i)

       ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
       call psb_get_rank (proc_to_comm, this%icontxt, proc_to_comm-1)

       ! Message size to be sent
       sizmsg = (snd_ptrs(i+1)-snd_ptrs(i))*chunk_size_

       if ( (sizmsg > 0) .and. (list_snd(i)-1 /= this%rank) ) then 
          call mpi_isend(sndbuf((snd_ptrs(i)-1)*chunk_size_+1), sizmsg, &
               & psb_mpi_long_integer, proc_to_comm,    &
               & psb_double_swap_tag, the_mpi_comm, sndhd(i), iret)
          check ( iret == mpi_success )
       end if
    end do

    ! Wait on all non-blocking receives
    do i=1, num_rcv
       proc_to_comm = list_rcv(i)

       ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
       call psb_get_rank (proc_to_comm, this%icontxt, proc_to_comm-1)

       ! Message size to be received
       sizmsg = (rcv_ptrs(i+1)-rcv_ptrs(i))*chunk_size_

       if ( (sizmsg > 0) .and. (list_rcv(i)-1 /= this%rank) ) then
          call mpi_wait(rcvhd(i), p2pstat, iret)
          check ( iret == mpi_success )
       else if ( list_rcv(i)-1 == this%rank ) then
          if ( sizmsg /= (snd_ptrs(i+1)-snd_ptrs(i))*chunk_size_ ) then 
             write(0,*) 'Fatal error in single_exchange: mismatch on self sendf', & 
                  & sizmsg, snd_ptrs(i+1)-snd_ptrs(i) 
             check(.false.)     
          end if
          rcvbuf( (rcv_ptrs(i)-1)*chunk_size_+1:(rcv_ptrs(i)-1)*chunk_size_+sizmsg) = &
               sndbuf( (snd_ptrs(i)-1)*chunk_size_+1:(snd_ptrs(i)-1)*chunk_size_+sizmsg )
       end if
    end do

    ! Finally wait on all non-blocking sends
    do i=1, num_snd
       proc_to_comm = list_snd(i)

       ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
       call psb_get_rank (proc_to_comm, this%icontxt, proc_to_comm-1)

       ! Message size to be received
       sizmsg = (snd_ptrs(i+1)-snd_ptrs(i))*chunk_size_

       if ( (sizmsg > 0) .and. (list_snd(i)-1 /= this%rank) ) then
          call mpi_wait(sndhd(i), p2pstat, iret)
          check ( iret == mpi_success )
       end if
    end do

    ! Unpack recv buffers
    call unpack_igp (rcv_ptrs(num_rcv+1)-rcv_ptrs(1), chunk_size_, unpack_idx, rcvbuf, x )

    call memfree (rcvhd,__FILE__,__LINE__) 
    call memfree (sndhd,__FILE__,__LINE__)

    call memfree (sndbuf,__FILE__,__LINE__)
    call memfree (rcvbuf,__FILE__,__LINE__)
  end subroutine par_context_neighbours_exchange_igp

  !=============================================================================
  subroutine par_context_neighbours_exchange_single_ip ( this, & 
       &                                                    num_neighbours, &
       &                                                    list_neighbours, &
       &                                                    input_data,&
       &                                                    output_data)
    !use psb_const_mod_names
    !use psb_penv_mod_names
    implicit none
    class(par_context_t), intent(in) :: this

    integer                 , intent(in)    :: num_neighbours
    integer(ip)             , intent(in)    :: list_neighbours (num_neighbours)
    integer(ip)             , intent(in)    :: input_data
    integer(ip)             , intent(inout) :: output_data(num_neighbours)

    integer(ip), allocatable :: ptrs(:)        ! How much data does the part send/recv to/from each neighbour?
    integer(ip), allocatable :: unpack_idx(:)  ! Where the data received from each neighbour is copied/added 
    ! on the local vectors of the part ?
    integer(ip), allocatable :: pack_idx(:)    ! Where is located the data to be sent to 
    ! each neighbour on the local vectors of the part ?

    integer(ip), allocatable :: buffer(:)  
    integer(ip)              :: i 

    call memalloc ( num_neighbours+1, ptrs, __FILE__, __LINE__ )
    ptrs(1)=1
    do i=2, num_neighbours+1
       ptrs(i)=ptrs(i-1)+1
    end do

    call memalloc ( ptrs(num_neighbours+1)-1, pack_idx, __FILE__, __LINE__ )
    pack_idx = 1 

    call memalloc ( ptrs(num_neighbours+1)-1, unpack_idx, __FILE__, __LINE__ )
    do i=1, ptrs(num_neighbours+1)-1
       unpack_idx(i) = i + 1
    end do

    call memalloc ( num_neighbours+1, buffer, __FILE__, __LINE__ )
    buffer(1) = input_data

    call this%neighbours_exchange ( num_neighbours,    &
         list_neighbours,   &
         ptrs,              &
         unpack_idx,        &  
         num_neighbours,    &
         list_neighbours,   &
         ptrs,              &
         pack_idx,          &
         buffer )

    output_data = buffer(2:)

    call memfree (buffer    , __FILE__, __LINE__ )
    call memfree (pack_idx  , __FILE__, __LINE__ )
    call memfree (unpack_idx, __FILE__, __LINE__ )
    call memfree (ptrs      , __FILE__, __LINE__ )
  end subroutine par_context_neighbours_exchange_single_ip

  !=============================================================================
  subroutine par_context_neighbours_exchange_wo_pack_unpack_ieep ( this, &
       &                                                              number_neighbours, &
       &                                                              neighbour_ids, &
       &                                                              snd_ptrs, &
       &                                                              snd_buf, & 
       &                                                              rcv_ptrs, &
       &                                                              rcv_buf )
    !     use psb_const_mod_names
    !     use psb_penv_mod_names
    ! #ifdef MPI_MOD
    !     use mpi
    ! #endif
    !      implicit none
    ! #ifdef MPI_H
    !      include 'mpif.h'
    ! #endif
    ! Parameters
    class(par_context_t)  , intent(in)    :: this 
    integer(ip)           , intent(in)    :: number_neighbours
    integer(ip)           , intent(in)    :: neighbour_ids(number_neighbours)
    integer(ip)           , intent(in)    :: snd_ptrs(number_neighbours+1)
    integer(ieep)         , intent(in)    :: snd_buf(snd_ptrs(number_neighbours+1)-1)   
    integer(ip)           , intent(in)    :: rcv_ptrs(number_neighbours+1)
    integer(ieep)         , intent(out)   :: rcv_buf(rcv_ptrs(number_neighbours+1)-1)

    ! Communication related locals 
    integer(ip) :: icontxt 
    integer     :: my_pid, proc_to_comm, sizmsg
    integer     :: the_mpi_comm,  iret
    integer     :: p2pstat(mpi_status_size)

    ! Request handlers for non-blocking receives
    integer, allocatable, dimension(:) :: rcvhd

    ! Request handlers for non-blocking receives
    integer, allocatable, dimension(:) :: sndhd

    integer(ip) :: i

    call memalloc (number_neighbours, rcvhd, __FILE__,__LINE__)
    call memalloc (number_neighbours, sndhd, __FILE__,__LINE__)

    call psb_get_mpicomm (this%icontxt, the_mpi_comm)

    ! First post all the non blocking receives   
    do i=1, number_neighbours
       proc_to_comm = neighbour_ids(i)

       ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
       call psb_get_rank (proc_to_comm, this%icontxt, proc_to_comm-1)

       ! Message size to be received
       sizmsg = rcv_ptrs(i+1)-rcv_ptrs(i)

       if ( (sizmsg > 0) .and. (neighbour_ids(i)-1 /= this%rank) ) then
          call mpi_irecv(  rcv_buf(rcv_ptrs(i)), sizmsg, &
               &  psb_mpi_integer1, proc_to_comm, &
               &  psb_double_swap_tag, the_mpi_comm, rcvhd(i), iret)
          check ( iret == mpi_success )
       end if
    end do

    ! Secondly post all non-blocking sends
    do i=1, number_neighbours
       proc_to_comm = neighbour_ids(i)

       ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
       call psb_get_rank (proc_to_comm, this%icontxt, proc_to_comm-1)

       ! Message size to be sent
       sizmsg = snd_ptrs(i+1)-snd_ptrs(i)

       if ( (sizmsg > 0) .and. (neighbour_ids(i)-1 /= this%rank) ) then 
          call mpi_isend(snd_buf(snd_ptrs(i)), sizmsg, &
               & psb_mpi_integer1, proc_to_comm, &
               & psb_double_swap_tag, the_mpi_comm, sndhd(i), iret)
          check ( iret == mpi_success )
       end if
    end do

    ! Wait on all non-blocking receives
    do i=1, number_neighbours
       proc_to_comm = neighbour_ids(i)

       ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
       call psb_get_rank (proc_to_comm, this%icontxt, proc_to_comm-1)

       ! Message size to be received
       sizmsg = rcv_ptrs(i+1)-rcv_ptrs(i)

       if ( (sizmsg > 0) .and. (neighbour_ids(i)-1 /= this%rank) ) then
          call mpi_wait(rcvhd(i), p2pstat, iret)
          check (iret == mpi_success)
       else if ( neighbour_ids(i)-1 == this%rank ) then
          if ( sizmsg /= snd_ptrs(i+1)-snd_ptrs(i) ) then 
             write(0,*) 'Fatal error in single_exchange: mismatch on self sendf', & 
                  & sizmsg, snd_ptrs(i+1)-snd_ptrs(i)
          end if
          rcv_buf( rcv_ptrs(i):rcv_ptrs(i+1)-1 ) = &
               snd_buf( snd_ptrs(i):snd_ptrs(i+1)-1 )
       end if
    end do

    ! Finally wait on all non-blocking sends
    do i=1, number_neighbours
       proc_to_comm = neighbour_ids(i)

       ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
       call psb_get_rank (proc_to_comm, this%icontxt, proc_to_comm-1)

       ! Message size to be received
       sizmsg = snd_ptrs(i+1)-snd_ptrs(i)

       if ( (sizmsg > 0) .and. (neighbour_ids(i)-1 /= this%rank) ) then
          call mpi_wait(sndhd(i), p2pstat, iret)
          check ( iret == mpi_success )
       end if
    end do

    call memfree (rcvhd ,__FILE__,__LINE__) 
    call memfree (sndhd ,__FILE__,__LINE__)
  end subroutine par_context_neighbours_exchange_wo_pack_unpack_ieep

  !=============================================================================
  !=============================================================================
  subroutine par_context_gather_scalar_ip ( this, root, input_data, output_data )
    implicit none
    class(par_context_t), intent(in)   :: this
    integer             , intent(in)   :: root
    integer(ip)         , intent(in)   :: input_data
    integer(ip)         , intent(out)  :: output_data(this%size)
    integer     :: the_mpi_comm, iret
    call psb_get_mpicomm (this%icontxt, the_mpi_comm)
    call mpi_gather( input_data, 1, psb_mpi_integer, output_data, 1, psb_mpi_integer, root, the_mpi_comm, iret)
    check( iret == mpi_success )
  end subroutine par_context_gather_scalar_ip

  !=============================================================================
  subroutine par_context_scatter_scalar_ip ( this, root, input_data, output_data )
    implicit none
    class(par_context_t), intent(in)   :: this
    integer                 , intent(in)   :: root
    integer(ip)             , intent(in)   :: input_data(this%size)
    integer(ip)             , intent(out)  :: output_data
    integer     :: the_mpi_comm, iret
    call psb_get_mpicomm (this%icontxt, the_mpi_comm)
    call mpi_scatter( input_data, 1, psb_mpi_integer, output_data, 1, psb_mpi_integer, root, the_mpi_comm, iret)
    check( iret == mpi_success )
  end subroutine par_context_scatter_scalar_ip

  !=============================================================================
  subroutine par_context_bcast_scalar_ip ( this, root, data )
    implicit none
    class(par_context_t), intent(in)    :: this
    integer                 , intent(in)    :: root
    integer(ip)             , intent(inout) :: data
    call psb_bcast ( this%icontxt, data, root=root)
  end subroutine par_context_bcast_scalar_ip

  !=============================================================================
  !=============================================================================
  subroutine pack_rp ( n, pack_idx, alpha, x, y )
    implicit none

    !Parameters
    integer (ip), intent(in)   :: n
    integer (ip), intent(in)   :: pack_idx(n)
    real    (rp), intent(in)   :: alpha
    real    (rp), intent(in)   :: x(*)
    real    (rp), intent(inout):: y(*)

    !Locals
    integer(ip) :: i

    if (alpha == 0.0_rp) then 
       !do nothing
    else if (alpha == 1.0_rp) then 
       do i=1,n
          y(i) = x(pack_idx(i))
       end do
    else if (alpha == -1.0_rp) then 
       do i=1,n
          y(i) = x(pack_idx(i))
       end do
    else  
       do i=1,n
          y(i) = alpha*x(pack_idx(i))
       end do
    end if

  end subroutine pack_rp

  !=============================================================================
  subroutine unpack_rp ( n, unpack_idx, beta, x, y )
    implicit none

    !Parameters
    integer(ip), intent(in)    :: n
    integer(ip), intent(in)    :: unpack_idx(n)
    real(rp)   , intent(in)    :: beta
    real(rp)   , intent(in)    :: x(*)
    real(rp)   , intent(inout) :: y(*)

    !Locals
    integer(ip) :: i

    if (beta == 0.0_rp) then
       do i=1,n
          y(unpack_idx(i)) = x(i)
       end do
    else if (beta == 1.0_rp) then
       do i=1,n
          y(unpack_idx(i)) = y(unpack_idx(i)) + x(i)
       end do
    else
       do i=1,n
          y(unpack_idx(i)) = beta*y(unpack_idx(i)) + x(i)
       end do
    end if
  end subroutine unpack_rp

  !=============================================================================
  subroutine pack_ip ( n, chunk_size, pack_idx, x, y )
    implicit none

    ! Parameters
    integer (ip), intent(in)     :: n
    integer (ip), intent(in)     :: chunk_size
    integer (ip), intent(in)     :: pack_idx(n)
    integer (ip), intent(in)    :: x(*)
    integer (ip), intent(inout) :: y(*)

    ! Locals
    integer(ip) :: i, j, startx, endx
    integer(ip) :: current
    current=1
    do i=1,n
       startx = (pack_idx(i)-1)*chunk_size + 1
       endx   = startx + chunk_size - 1
       do j=startx, endx
          y(current) = x(j)
          current = current + 1
       end do
    end do
  end subroutine pack_ip

  !=============================================================================
  subroutine unpack_ip ( n, chunk_size, unpack_idx, x, y )
    implicit none

    ! Parameters
    integer (ip), intent(in)     :: n
    integer (ip), intent(in)     :: chunk_size
    integer (ip), intent(in)     :: unpack_idx(n)
    integer (ip), intent(in)    :: x(*)
    integer (ip), intent(inout) :: y(*)

    ! Locals
    integer(ip) :: i, j, starty, endy, current
    current = 1
    do i=1,n
       starty = (unpack_idx(i)-1)*chunk_size + 1
       endy   = starty + chunk_size - 1
       do j=starty, endy
          y(j) = x(current)
          current = current + 1
       end do
    end do
  end subroutine unpack_ip

  !=============================================================================
  subroutine pack_igp ( n, chunk_size, pack_idx, x, y )
    implicit none

    !Parameters
    integer (ip), intent(in)     :: n
    integer (ip), intent(in)     :: chunk_size
    integer (ip), intent(in)     :: pack_idx(n)
    integer (igp), intent(in)    :: x(*)
    integer (igp), intent(inout) :: y(*)

    !Locals
    integer(ip) :: i, j, startx, endx
    integer(ip) :: current
    current=1
    do i=1,n
       startx = (pack_idx(i)-1)*chunk_size + 1
       endx   = startx + chunk_size - 1
       do j=startx, endx
          y(current) = x(j)
          current = current + 1
       end do
    end do
  end subroutine pack_igp

  !=============================================================================
  subroutine unpack_igp ( n, chunk_size, unpack_idx, x, y )
    implicit none

    !Parameters
    integer (ip), intent(in)     :: n
    integer (ip), intent(in)     :: chunk_size
    integer (ip), intent(in)     :: unpack_idx(n)
    integer (igp), intent(in)    :: x(*)
    integer (igp), intent(inout) :: y(*)

    !Locals
    integer(ip) :: i, j, starty, endy, current
    current = 1
    do i=1,n
       starty = (unpack_idx(i)-1)*chunk_size + 1
       endy   = starty + chunk_size - 1
       do j=starty, endy
          y(j) = x(current)
          current = current + 1
       end do
    end do
  end subroutine unpack_igp
  !=============================================================================
  !=============================================================================
  subroutine par_context_transfer_to_master_ip ( this, input_data, output_data )
    implicit none
    class(par_context_t), intent(in)      :: this
    integer(ip)         , intent(in)      :: input_data
    integer(ip)         , intent(inout)   :: output_data
    integer(ip) :: recv_rank
    integer(ip) :: send_rank
    send_rank = 0
    recv_rank = this%size-1
    if ( this%rank == send_rank ) then
       call psb_snd(this%icontxt, input_data, recv_rank)
    else if (this%rank == recv_rank) then
       call psb_rcv(this%icontxt, output_data, send_rank)
    end if
  end subroutine par_context_transfer_to_master_ip

  !=============================================================================
  subroutine par_context_transfer_to_master_ip_1D_array ( this, input_data, output_data )
    implicit none
    class(par_context_t), intent(in)      :: this
    integer(ip)         , intent(in)      :: input_data(:)
    integer(ip)         , intent(inout)   :: output_data(:)
    integer(ip) :: recv_rank
    integer(ip) :: send_rank
    send_rank = 0
    recv_rank = this%size-1 
    if ( this%rank == send_rank ) then
       call psb_snd(this%icontxt, input_data, recv_rank)
    else if (this%rank == recv_rank) then
       call psb_rcv(this%icontxt, output_data, send_rank)
    end if
  end subroutine par_context_transfer_to_master_ip_1D_array

  !=============================================================================
  !=============================================================================
  subroutine par_context_gather_to_master_ip ( this, input_data, output_data )
    implicit none
    class(par_context_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data
    integer(ip)             , intent(out)  :: output_data(this%size)
    integer                :: the_mpi_comm, iret, root
    root    = this%size - 1 
    call psb_get_mpicomm (this%icontxt, the_mpi_comm)
    call mpi_gather( input_data, 1, psb_mpi_integer, output_data, 1, psb_mpi_integer, root, the_mpi_comm, iret)
    check( iret == mpi_success )
  end subroutine par_context_gather_to_master_ip

  !=============================================================================
  subroutine par_context_gather_to_master_igp ( this, input_data, output_data )
    implicit none
    class(par_context_t), intent(in)   :: this
    integer(igp)            , intent(in)   :: input_data
    integer(igp)            , intent(out)  :: output_data(this%size)
    integer                :: the_mpi_comm, iret, root
    root    = this%size - 1 
    call psb_get_mpicomm (this%icontxt, the_mpi_comm)
    call mpi_gather( input_data, 1, psb_mpi_long_integer, output_data, 1, psb_mpi_long_integer, root, the_mpi_comm, iret)
    check( iret == mpi_success )
  end subroutine par_context_gather_to_master_igp

  !=============================================================================
  subroutine par_context_gather_to_master_ip_1D_array ( this, input_data_size, input_data, output_data )
    implicit none
    class(par_context_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data_size
    integer(ip)             , intent(in)   :: input_data(input_data_size)
    integer(ip)             , intent(out)  :: output_data(:)
    integer                :: the_mpi_comm, iret, root
    root    = this%size - 1 
    call psb_get_mpicomm (this%icontxt, the_mpi_comm)
    call mpi_gather( input_data,  input_data_size, psb_mpi_integer, &
         output_data, input_data_size, psb_mpi_integer, &
         root, the_mpi_comm, iret)
    check( iret == mpi_success )
  end subroutine par_context_gather_to_master_ip_1D_array

  !=============================================================================
  subroutine par_context_gather_to_masterv_ip_1D_array ( this, input_data_size, input_data, recv_counts, displs, output_data )
    implicit none
    class(par_context_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data_size
    integer(ip)             , intent(in)   :: input_data(input_data_size)
    integer(ip)             , intent(in)   :: recv_counts(this%size)
    integer(ip)             , intent(in)   :: displs(this%size)
    integer(ip)             , intent(out)  :: output_data(:)
    integer                :: the_mpi_comm, iret, root
    root    = this%size - 1 
    call psb_get_mpicomm (this%icontxt, the_mpi_comm)
    call mpi_gatherv( input_data, input_data_size, psb_mpi_integer, &
         output_data, recv_counts, displs, psb_mpi_integer, &
         root, the_mpi_comm, iret)
    check( iret == mpi_success )
  end subroutine par_context_gather_to_masterv_ip_1D_array

  !=============================================================================
  subroutine par_context_gather_to_masterv_igp_1D_array ( this, input_data_size, input_data, recv_counts, displs, output_data )
    implicit none
    class(par_context_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data_size
    integer(igp)            , intent(in)   :: input_data(input_data_size)
    integer(ip)             , intent(in)   :: recv_counts(this%size)
    integer(ip)             , intent(in)   :: displs(this%size)
    integer(igp)            , intent(out)  :: output_data(:)
    integer                :: the_mpi_comm, iret, root
    root    = this%size - 1 
    call psb_get_mpicomm (this%icontxt, the_mpi_comm)
    call mpi_gatherv( input_data, input_data_size, psb_mpi_long_integer, &
         output_data, recv_counts, displs, psb_mpi_long_integer, &
         root, the_mpi_comm, iret)
    check( iret == mpi_success )
  end subroutine par_context_gather_to_masterv_igp_1D_array

  !=============================================================================
  subroutine par_context_gather_to_masterv_rp_1D_array ( this, input_data_size, input_data, recv_counts, displs, output_data )
    implicit none
    class(par_context_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data_size
    real(rp)                , intent(in)   :: input_data(input_data_size)
    integer(ip)             , intent(in)   :: recv_counts(this%size)
    integer(ip)             , intent(in)   :: displs(this%size)
    real(rp)                , intent(out)  :: output_data(:)
    integer                :: the_mpi_comm, iret, root
    root    = this%size - 1 
    call psb_get_mpicomm (this%icontxt, the_mpi_comm)
    call mpi_gatherv( input_data , input_data_size, psb_mpi_real, &
         output_data, recv_counts, displs, psb_mpi_real, &
         root, the_mpi_comm, iret)
    check( iret == mpi_success )
  end subroutine par_context_gather_to_masterv_rp_1D_array

  !=============================================================================
  subroutine par_context_gather_to_masterv_rp_2D_array ( this, input_data, recv_counts, displs, output_data )
    implicit none
    class(par_context_t), intent(in)   :: this
    real(rp)                , intent(in)   :: input_data(:,:)
    integer(ip)             , intent(in)   :: recv_counts(this%size)
    integer(ip)             , intent(in)   :: displs(this%size)
    real(rp)                , intent(out)  :: output_data(:)
    integer                :: the_mpi_comm, iret, root
    root    = this%size - 1 
    call psb_get_mpicomm (this%icontxt, the_mpi_comm)
    call mpi_gatherv( input_data , size(input_data,1)*size(input_data,2), psb_mpi_real, &
         output_data, recv_counts, displs, psb_mpi_real, &
         root, the_mpi_comm, iret)
    check( iret == mpi_success )
  end subroutine par_context_gather_to_masterv_rp_2D_array

  !=============================================================================
  subroutine par_context_scatter_from_masterv_rp_1D_array ( this, input_data, send_counts, displs, output_data_size, output_data )
    implicit none
    class(par_context_t), intent(in)   :: this
    real(rp)                , intent(in)   :: input_data(:)
    integer(ip)             , intent(in)   :: send_counts(this%size)
    integer(ip)             , intent(in)   :: displs(this%size)
    integer(ip)             , intent(in)   :: output_data_size
    real(rp)                , intent(out)  :: output_data(output_data_size)
    integer                :: the_mpi_comm, iret, root
    root    = this%size - 1 
    call psb_get_mpicomm (this%icontxt, the_mpi_comm)
    call mpi_scatterv( input_data, send_counts, displs, psb_mpi_real, &
         output_data, output_data_size, psb_mpi_real, &
         root, the_mpi_comm, iret)
    check( iret == mpi_success )
  end subroutine par_context_scatter_from_masterv_rp_1D_array

end module par_context_names
