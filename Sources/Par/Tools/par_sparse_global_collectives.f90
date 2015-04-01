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
module par_sparse_global_collectives
   ! Serial modules
   use types
   use memor

   implicit none
   private

  !=============================================
  ! TODO: extend module to support block vectors
  !       (already DONE!!!)
  !=============================================

#include "debug.i90"

   ! Alternative implementations of sparse all to all
   integer(ip), parameter :: sp_all_to_all_all_to_all    = 0
   integer(ip), parameter :: sp_all_to_all_psnd_prcv     = 1     
   integer(ip), parameter :: sp_all_to_all_only_psnd     = 2 
   integer(ip), parameter :: sp_all_to_all_only_prcv     = 3
   integer(ip), parameter :: sp_all_to_all_ircv_and_rsnd = 4
   integer(ip), parameter :: sp_all_to_all_ircv_and_snd  = 5
   integer(ip), parameter :: sp_all_to_all_ircv_and_isnd = 6


   ! Default implementation of sparse_all_to_all
   integer(ip), parameter :: sp_all_to_all_default = & 
                          &  sp_all_to_all_ircv_and_isnd   

   ! Constants
   public :: sp_all_to_all_default,                                   &  
          &  sp_all_to_all_all_to_all, sp_all_to_all_psnd_prcv,       &
          &  sp_all_to_all_only_psnd , sp_all_to_all_only_prcv,       &
          &  sp_all_to_all_ircv_and_rsnd, sp_all_to_all_ircv_and_snd, &
          &  sp_all_to_all_ircv_and_isnd
   ! Types
   ! public ::
   
   ! Subroutines
   public :: single_exchange
contains

   ! When packing   (gathering) ,    buffer <- alpha * x
   ! When unpacking (scattering),    x <- beta*x + buffer
   subroutine single_exchange ( icontxt, & 
                                num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
                                num_snd, list_snd, snd_ptrs, pack_idx,   &
                                alpha, beta, x, work, mode)
     use psb_const_mod
     use psb_penv_mod
#ifdef MPI_MOD
     use mpi
#endif
     implicit none
#ifdef MPI_H
     include 'mpif.h'
#endif
     
     ! Parameters
     integer, intent(in) :: icontxt 
     
     ! **IMPORTANT NOTE**: I will assume that both 
     ! list_rcv and list_snd hold process identifiers 
     ! starting from one. However, the underlying
     ! message-passing library identifies processes 
     ! starting from zero ...

     ! Control info to receive
     integer    , intent(in) :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
     integer(ip), intent(in) :: unpack_idx (rcv_ptrs(num_rcv+1)-1)

     ! Control info to send
     integer    , intent(in) :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
     integer(ip), intent(in) :: pack_idx (snd_ptrs(num_snd+1)-1)

     ! Floating point data
     real(rp), intent(in)    :: alpha, beta

     ! **IMPORTANT NOTE**: x must be of size 
     ! size(x_I + x_G) or of size size(x_G) depending 
     ! on whether pack/unpack arrays contain indices relative
     ! to the beginning of the whole local array or
     ! to the beginning of the interfaces (i.e., x_G)
     ! I have declared x has an assumed-size because 
     ! it accomodates both options
     real(rp), intent(inout) :: x(*)

     ! Optional work array to store snd/rcv
     ! buffers. If work does not have enough room, then
     ! the routine allocates room for it
     real(rp), target, optional, intent(in) :: work(:)
     
     ! Concrete implementation of sparse all to all
     ! (see module header)
     integer(ip), optional, intent(in) :: mode
      
     ! Local variables 
     logical do_pack, do_unpack, valid_mode, albf
     
     ! Communication related locals 
     integer :: my_pid, num_procs, i, proc_to_comm, sizmsg
     integer :: mpi_comm,  iret, info
     integer :: p2pstat(mpi_status_size)


     ! Arrays required by mpi_all_to_all
     integer, allocatable, dimension(:) :: sndidx, rcvidx, &
                                        &  sndsiz, rcvsiz

     ! Request handlers for non-blocking receives
     integer, allocatable, dimension(:) :: rcvhd

     ! Request handlers for non-blocking receives
     integer, allocatable, dimension(:) :: sndhd

     integer(ip)       :: mode_
     real(rp), pointer :: sndbuf(:) => NULL(), rcvbuf(:) => NULL()
   
     if ( present(mode) ) then
       mode_ = mode
     else
       mode_ = sp_all_to_all_default
     end if

     ! Check the set of (currently) supported implementations
     valid_mode = mode_ == sp_all_to_all_all_to_all     .or.  &
               &  mode_ == sp_all_to_all_psnd_prcv      .or.  & 
               &  mode_ == sp_all_to_all_only_psnd      .or.  & 
               &  mode_ == sp_all_to_all_only_prcv      .or.  & 
               &  mode_ == sp_all_to_all_ircv_and_rsnd  .or.  &
               &  mode_ == sp_all_to_all_ircv_and_snd   .or.  & 
               &  mode_ == sp_all_to_all_ircv_and_isnd

     assert ( valid_mode )

     do_pack = mode_ == sp_all_to_all_all_to_all     .or.  & 
            &  mode_ == sp_all_to_all_psnd_prcv      .or.  &
            &  mode_ == sp_all_to_all_only_psnd      .or.  &
            &  mode_ == sp_all_to_all_ircv_and_rsnd  .or.  &
            &  mode_ == sp_all_to_all_ircv_and_snd   .or.  &
            &  mode_ == sp_all_to_all_ircv_and_isnd

     do_unpack =  mode_ == sp_all_to_all_all_to_all     .or.  &
               &  mode_ == sp_all_to_all_psnd_prcv      .or.  &
               &  mode_ == sp_all_to_all_only_prcv      .or.  & 
               &  mode_ == sp_all_to_all_ircv_and_rsnd  .or.  &
               &  mode_ == sp_all_to_all_ircv_and_snd   .or.  & 
               &  mode_ == sp_all_to_all_ircv_and_isnd
     
    
     ! Get my process identifier and number of
     ! of processors available in the parallel
     ! context
     call psb_info(icontxt, my_pid, num_procs)
     
     ! Get MPI communicator associated to icontxt (in
     ! the current implementation of our wrappers
     ! to the MPI library icontxt and mpi_comm are actually 
     ! the same)
     call psb_get_mpicomm (icontxt, mpi_comm)

     if (mode_ == sp_all_to_all_all_to_all) then 
       call memalloc (num_procs, sndidx, __FILE__,__LINE__)
       call memalloc (num_procs, rcvidx, __FILE__,__LINE__)
       call memalloc (num_procs, sndsiz, __FILE__,__LINE__)
       call memalloc (num_procs, rcvsiz, __FILE__,__LINE__)
       
       sndsiz(:) = 0
       sndidx(:) = 0 
       rcvsiz(:) = 0
       rcvidx(:) = 0
       
       ! Prepare rcv info.
       do i=1, num_rcv
           proc_to_comm = list_rcv(i)

           ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
           call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)
           
           ! rcvidx and rcvsiz start from one ! 
           rcvidx( proc_to_comm+1 ) = rcv_ptrs(i)-1
           rcvsiz( proc_to_comm+1 ) = rcv_ptrs (i+1)-rcv_ptrs(i)
       end do 
    
       ! Prepare snd info.
       do i=1, num_snd
           proc_to_comm = list_snd(i)

           ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
           call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)
           
           ! sndsiz and sndsiz start from one ! 
           sndidx( proc_to_comm+1 ) = snd_ptrs(i)-1
           sndsiz( proc_to_comm+1 ) = snd_ptrs (i+1)-snd_ptrs(i)
       end do

       ! write (*,*) 'RCVIDX', rcvidx ! DBG: 
       ! write (*,*) 'RCVSIZ', rcvsiz ! DBG:

       ! write (*,*) 'SNDIDX', sndidx ! DBG:
       ! write (*,*) 'SNDSIZ', sndsiz ! DBG:
    else
       call memalloc (num_rcv, rcvhd, __FILE__,__LINE__)
       if (  mode_ == sp_all_to_all_ircv_and_isnd ) call memalloc (num_snd, sndhd, __FILE__,__LINE__)
    end if

    ! Prepare room for sndbuf and rcvbuf
    if (present(work)) then 
       if ( (snd_ptrs (num_snd+1)-snd_ptrs (1) + & 
             rcv_ptrs (num_rcv+1)-rcv_ptrs (1)) <= size(work,1) ) then
         sndbuf => work(1:snd_ptrs(num_snd+1)-1)
         rcvbuf => work(snd_ptrs(num_snd+1):snd_ptrs(num_snd+1)+rcv_ptrs(num_rcv+1)-2)
         albf = .false.
       else
         call memallocp ((snd_ptrs(num_snd+1)-snd_ptrs(1)), sndbuf, __FILE__,__LINE__)
         call memallocp ((rcv_ptrs(num_rcv+1)-rcv_ptrs(1)), rcvbuf, __FILE__,__LINE__)
         albf=.true.             
       endif 
    else
      call memallocp ((snd_ptrs(num_snd+1)-snd_ptrs(1)), sndbuf, __FILE__,__LINE__)
      call memallocp ((rcv_ptrs(num_rcv+1)-rcv_ptrs(1)), rcvbuf, __FILE__,__LINE__)
      albf=.true. 
    end if

    if (do_pack) then
      ! Pack send buffers
      call pack ( snd_ptrs(num_snd+1)-snd_ptrs(1), pack_idx, alpha, x, sndbuf )
      ! write (*,*) 'P', sndbuf       ! DBG:
      ! write (*,*) 'PIDX', pack_idx  ! DBG:
    end if

     if (mode_ == sp_all_to_all_all_to_all) then
       ! exchange data using mpi_alltoallv
       call mpi_alltoallv( sndbuf, sndsiz, sndidx, psb_mpi_real, &
                           rcvbuf, rcvsiz, rcvidx, psb_mpi_real, mpi_comm, iret)
       
       if ( iret /= mpi_success ) then
         write (0,*) 'Error: mpi_alltoallv returned != mpi_success'
         call psb_abort (icontxt)    
       end if

     else if ( mode_ == sp_all_to_all_psnd_prcv .or. &
             & mode_ == sp_all_to_all_only_psnd .or. &
             & mode_ == sp_all_to_all_only_prcv )  then

       if ( mode_ == sp_all_to_all_psnd_prcv .or. &
          & mode_ == sp_all_to_all_only_psnd ) then
         ! Firstly, post all (locally) blocking sends    
         do i=1, num_snd
           proc_to_comm = list_snd(i)-1           
           ! Message size to be sent
           sizmsg = snd_ptrs(i+1)-snd_ptrs(i)
           if ( (sizmsg > 0) .and. (proc_to_comm /= my_pid) ) then 
              call psb_snd( icontxt, sndbuf(snd_ptrs(i):snd_ptrs(i)+sizmsg-1), & 
                          & proc_to_comm )
           end if
         end do
       end if 

       if ( mode_ == sp_all_to_all_psnd_prcv .or. & 
          & mode_ == sp_all_to_all_only_prcv ) then
         ! write (*,*) 'RCVDDD !!!!'
         ! Secondly, post all blocking receives
         do i=1, num_rcv
           proc_to_comm = list_rcv(i)-1
           ! Message size to be received
           sizmsg = rcv_ptrs(i+1)-rcv_ptrs(i)
           ! write (*,*) 'PSBRCV',  sizmsg, proc_to_comm, my_pid ! DBG:
           if ( (sizmsg > 0) .and. (proc_to_comm /= my_pid) ) then
              call psb_rcv( icontxt, rcvbuf(rcv_ptrs(i):rcv_ptrs(i)+sizmsg-1), & 
                          & proc_to_comm )               
           else if ( proc_to_comm == my_pid ) then
             if ( sizmsg /= rcv_ptrs(i+1)-rcv_ptrs(i) ) then 
                write(0,*) 'Fatal error in single_exchange: mismatch on self sendf', & 
                          & sizmsg, rcv_ptrs(i+1)-rcv_ptrs(i)
             end if
             rcvbuf( rcv_ptrs(i):rcv_ptrs(i)+sizmsg-1 ) = &
                 sndbuf( snd_ptrs(i) : snd_ptrs(i)+sizmsg-1 )
           end if
          end do
        end if  

     else if ( mode_ == sp_all_to_all_ircv_and_rsnd .or. &
             & mode_ == sp_all_to_all_ircv_and_snd  .or. &
             & mode_ == sp_all_to_all_ircv_and_isnd ) then

        ! First post all the non blocking receives   
        do i=1, num_rcv
          proc_to_comm = list_rcv(i)
          
          ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
          call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)
          
          ! Message size to be received
          sizmsg = rcv_ptrs(i+1)-rcv_ptrs(i)
      
          if ( (sizmsg > 0) .and. (list_rcv(i)-1 /= my_pid) ) then
             call mpi_irecv(  rcvbuf(rcv_ptrs(i)), sizmsg,        &
                           &  psb_mpi_real, proc_to_comm, &
                           &  psb_double_swap_tag, mpi_comm, rcvhd(i), iret)

             if ( iret /= mpi_success ) then
                write (0,*) 'Error: mpi_irecv returned != mpi_success'
                call psb_abort (icontxt)    
             end if
          end if
        end do
    
        ! The following barrier is required by MPI_Rsend (see below). 
        ! MPI_Rsend requires that the matching processor has already 
        ! issued the matching receive. MPI_Rsend reduces overhead 
        ! associated to the hand-shaking protocol among sender and 
        ! receiver. 
        if ( mode_ == sp_all_to_all_ircv_and_rsnd  ) then
          call mpi_barrier(mpi_comm, info)
        end if 

        if ( mode_ == sp_all_to_all_ircv_and_rsnd .or. &
             & mode_ == sp_all_to_all_ircv_and_snd ) then

           ! Secondly, post all blocking sends    
           do i=1, num_snd
            proc_to_comm = list_snd(i)
          
            ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
            call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)
          
            ! Message size to be sent
            sizmsg = snd_ptrs(i+1)-snd_ptrs(i)
      
            if ( (sizmsg > 0) .and. (list_snd(i)-1 /= my_pid) ) then 
               if ( mode_ == sp_all_to_all_ircv_and_rsnd ) then 
                  call mpi_rsend(sndbuf(snd_ptrs(i)), sizmsg,  &
                     & psb_mpi_real, proc_to_comm,      &
                     & psb_double_swap_tag, mpi_comm, iret)
               else
                  call mpi_send(sndbuf(snd_ptrs(i)), sizmsg, &
                     & psb_mpi_real, proc_to_comm,    &
                     & psb_double_swap_tag, mpi_comm, iret)
               end if

               if ( iret /= mpi_success ) then
                  write (0,*) 'Error: mpi_rsend/mpi_send returned != mpi_success'
                  call psb_abort (icontxt)    
               end if
            end if
           end do
        else if ( mode_ == sp_all_to_all_ircv_and_isnd ) then 
           ! Secondly post all non-blocking sends
           do i=1, num_snd
            proc_to_comm = list_snd(i)
          
            ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
            call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)
          
            ! Message size to be sent
            sizmsg = snd_ptrs(i+1)-snd_ptrs(i)
      
            if ( (sizmsg > 0) .and. (list_snd(i)-1 /= my_pid) ) then 
                  call mpi_isend(sndbuf(snd_ptrs(i)), sizmsg, &
                     & psb_mpi_real, proc_to_comm,    &
                     & psb_double_swap_tag, mpi_comm, sndhd(i), iret)

               if ( iret /= mpi_success ) then
                  write (0,*) 'Error: mpi_isend returned != mpi_success'
                  call psb_abort (icontxt)    
               end if
            end if
           end do
        end if

        ! Wait on all non-blocking receives
        do i=1, num_rcv
          proc_to_comm = list_rcv(i)
          
          ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
          call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)
          
          ! Message size to be received
          sizmsg = rcv_ptrs(i+1)-rcv_ptrs(i)
      
          if ( (sizmsg > 0) .and. (list_rcv(i)-1 /= my_pid) ) then
             call mpi_wait(rcvhd(i), p2pstat, iret)

             if ( iret /= mpi_success ) then
                write (0,*) 'Error: mpi_wait returned != mpi_success'
                call psb_abort (icontxt)    
             end if
          else if ( list_rcv(i)-1 == my_pid ) then
            if ( sizmsg /= snd_ptrs(i+1)-snd_ptrs(i) ) then 
               write(0,*) 'Fatal error in single_exchange: mismatch on self sendf', & 
                        & sizmsg, snd_ptrs(i+1)-snd_ptrs(i) 
            end if

             rcvbuf( rcv_ptrs(i):rcv_ptrs(i)+sizmsg-1) = &
                 sndbuf( snd_ptrs(i): snd_ptrs(i)+sizmsg-1 )
          end if
        end do

        if ( mode_ == sp_all_to_all_ircv_and_isnd ) then 
          ! Finally wait on all non-blocking sends
          do i=1, num_snd
            proc_to_comm = list_snd(i)
          
            ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
            call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)
          
            ! Message size to be received
            sizmsg = snd_ptrs(i+1)-snd_ptrs(i)
      
            if ( (sizmsg > 0) .and. (list_snd(i)-1 /= my_pid) ) then
               call mpi_wait(sndhd(i), p2pstat, iret)
               if ( iret /= mpi_success ) then
                  write (0,*) 'Error: mpi_wait returned != mpi_success'
                  call psb_abort (icontxt)    
               end if
             end if
          end do
        end if
     end if

    if (do_unpack) then
      ! Unpack recv buffers
      call unpack (rcv_ptrs(num_rcv+1)-rcv_ptrs(1), unpack_idx, beta, rcvbuf, x )
      ! write(*,*)  'RCVPTRS', rcv_ptrs   ! DBG:  
      ! write (*,*) 'UIDX', unpack_idx    ! DBG:
      ! write (*,*) 'U', rcvbuf           ! DBG: 
    end if

    
    ! Free allocated memory ...
    if (mode_ == sp_all_to_all_all_to_all) then 
      call memfree (sndidx,__FILE__,__LINE__)
      call memfree (rcvidx,__FILE__,__LINE__)
      call memfree (sndsiz,__FILE__,__LINE__)
      call memfree (rcvsiz,__FILE__,__LINE__)
    else
      call memfree (rcvhd,__FILE__,__LINE__) 
      if (  mode_ == sp_all_to_all_ircv_and_isnd ) call memfree (sndhd,__FILE__,__LINE__)
    end if

    if (albf) then
      call memfreep (sndbuf,__FILE__,__LINE__)
      call memfreep (rcvbuf,__FILE__,__LINE__)  
    end if
   
   end subroutine single_exchange

   subroutine pack ( neq, pack_idx, alpha, x, y )
     implicit none

     ! Parameters
     integer (ip)  :: neq, pack_idx(:)
     real    (rp)  :: x(*), y(*), alpha

     ! Locals
     integer(ip) :: i, idof, base_l, base_r

     if (alpha == 0.0_rp) then 
        ! do nothing
     else if (alpha == 1.0_rp) then 
        do i=1,neq
           base_l = i
           base_r = pack_idx(i)
           y(base_l) = x(base_r)
           ! write (*,*) 'PACK', base_l, base_r, y(base_l), x(base_r)
           base_l = base_l + 1
           base_r = base_r + 1
        end do
     else if (alpha == -1.0_rp) then 
        do i=1,neq
           base_l = i
           base_r = pack_idx(i)
           y(base_l) = x( base_r )
           base_l = base_l + 1
           base_r = base_r + 1
        end do
     else  
        do i=1,neq
           base_l = i
           base_r = pack_idx(i)
           y(base_l) = alpha*x(base_r)
           base_l = base_l + 1
           base_r = base_r + 1
        end do
     end if

   end subroutine pack

   subroutine unpack ( neq, unpack_idx, beta, x, y )
     implicit none

     ! Parameters
     integer(ip) :: neq, unpack_idx(*)
     real(rp)    :: beta, x(*), y(*)

     ! Locals
     integer(ip) :: i, idof, base_l, base_r

     if (beta == 0.0_rp) then
        do i=1,neq
           base_r = i
           base_l = unpack_idx(i)
           y(base_l) = x(base_r)
           ! write (*,*) 'UNPACK', base_l, base_r, y(base_l), x(base_r)
           base_l = base_l + 1
           base_r = base_r + 1
        end do
     else if (beta == 1.0_rp) then
        do i=1,neq
           base_r = i
           base_l = unpack_idx(i)
           y(base_l) = y(base_l) + x(base_r)
           ! write (*,*) 'UNPACK', base_l, base_r, y(base_l), x(base_r)
           base_l = base_l + 1
           base_r = base_r + 1
        end do
     else
        do i=1,neq
           base_r = i
           base_l = unpack_idx(i)
           y(base_l) = beta*y(base_l) + x(base_r)
           base_l = base_l + 1
           base_r = base_r + 1
        end do
     end if
   end subroutine unpack

end module par_sparse_global_collectives
