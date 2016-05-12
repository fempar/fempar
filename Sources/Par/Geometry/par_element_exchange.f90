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
module par_element_exchange_names
   ! Serial modules
   use types_names
   use memor_names
   
   use par_context_names
   use par_environment_names
   use element_import_names
   use migratory_element_names
   implicit none
#include "debug.i90"
   private

   ! Subroutines
   public :: ghost_elements_exchange
contains

  subroutine ghost_elements_exchange ( par_environment, element_import, data )
    implicit none
    type(par_environment_t)   , intent(in)    :: par_environment
    type(element_import_t)    , intent(in)    :: element_import
    class(migratory_element_t), intent(inout) :: data(:)
    
    if ( par_environment%am_i_l1_task() ) then
      call plain_ghost_element_exchange ( par_environment, &
                                          element_import%get_number_neighbours(), &
                                          element_import%get_neighbours_ids(), &
                                          element_import%get_rcv_ptrs(), &
                                          element_import%get_rcv_leids(), &
                                          element_import%get_snd_ptrs(), &
                                          element_import%get_snd_leids(), &
                                          data)
    end if
  end subroutine ghost_elements_exchange

  subroutine plain_ghost_element_exchange ( par_environment, &
                                            number_neighbours, &
                                            neighbour_ids, &
                                            rcv_ptrs, &
                                            rcv_leids,&
                                            snd_ptrs, &
                                            snd_leids, &
                                            data)
                                             
    use psb_const_mod_names
    use psb_penv_mod_names
#ifdef MPI_MOD
    use mpi
#endif
     implicit none
#ifdef MPI_H
     include 'mpif.h'
#endif
     ! Parameters
     type(par_environment_t)   , intent(in)    :: par_environment 
     integer(ip)               , intent(in)    :: number_neighbours
     integer(ip)               , intent(in)    :: neighbour_ids(number_neighbours)
     integer(ip)               , intent(in)    :: rcv_ptrs(number_neighbours+1)
     integer(ip)               , intent(in)    :: rcv_leids(rcv_ptrs(number_neighbours+1)-1)
     integer(ip)               , intent(in)    :: snd_ptrs(number_neighbours+1)
     integer(ip)               , intent(in)    :: snd_leids(snd_ptrs(number_neighbours+1)-1)   
     class(migratory_element_t), intent(inout) :: data(:) 
     
     ! Communication related locals 
     integer(ip) :: icontxt 
     integer     :: my_pid, num_procs, proc_to_comm, sizmsg
     integer     :: the_mpi_comm,  iret, info
     integer     :: p2pstat(mpi_status_size)

     ! Request handlers for non-blocking receives
     integer, allocatable, dimension(:) :: rcvhd

     ! Request handlers for non-blocking receives
     integer, allocatable, dimension(:) :: sndhd

     integer(ip)  , allocatable :: elemsizes(:)
     integer(ip)  , allocatable :: rcv_ptrs_buf(:)
     integer(ip)  , allocatable :: snd_ptrs_buf(:)
     
     integer(ieep), allocatable :: sndbuf(:)  
     integer(ieep), allocatable :: rcvbuf(:)

     type(par_context_t), pointer :: l1_context
     
     integer(ip) :: current, i, j
     
     assert ( par_environment%am_i_l1_task() )
     
     l1_context => par_environment%get_l1_context()
     icontxt    = l1_context%get_icontxt()
     my_pid     = l1_context%get_rank()
     num_procs  = l1_context%get_size()

     call psb_get_mpicomm (icontxt, the_mpi_comm)

     call memalloc ( size(data), elemsizes, __FILE__, __LINE__ )
     call memalloc ( number_neighbours+1, snd_ptrs_buf, __FILE__, __LINE__ )
     call memalloc ( number_neighbours+1, rcv_ptrs_buf, __FILE__, __LINE__ )
    
     ! Set-up snd_ptrs_bufs using local information
     elemsizes = - 1
     snd_ptrs_buf = 0
     do i=1, number_neighbours
      do j=snd_ptrs(i),snd_ptrs(i+1)-1
       call data(snd_leids(j))%size(elemsizes(snd_leids(j)))  
       snd_ptrs_buf(i+1) = snd_ptrs_buf(i+1) + elemsizes(snd_leids(j))
      end do
     end do
     
     snd_ptrs_buf(1) = 1
     do i=1, number_neighbours
       snd_ptrs_buf(i+1) = snd_ptrs_buf(i) + snd_ptrs_buf(i+1)
     end do
     
     ! Set-up rcv_ptrs using data fetched from my neighbours
     call par_environment%l1_neighbours_exchange ( number_neighbours, &
                                                   neighbour_ids,&
                                                   rcv_ptrs,&
                                                   rcv_leids,&
                                                   number_neighbours,&
                                                   neighbour_ids,&
                                                   snd_ptrs,&
                                                   snd_leids,&
                                                   elemsizes )
     
     rcv_ptrs_buf = 0
     do i=1, number_neighbours
      do j=rcv_ptrs(i),rcv_ptrs(i+1)-1
       rcv_ptrs_buf(i+1) = rcv_ptrs_buf(i+1) + elemsizes(rcv_leids(j))
      end do
     end do
     
     rcv_ptrs_buf(1) = 1
     do i=1, number_neighbours
       rcv_ptrs_buf(i+1) = rcv_ptrs_buf(i) + rcv_ptrs_buf(i+1)
     end do
     
     ! Prepare room for sndbuf
     call memalloc (snd_ptrs_buf(number_neighbours+1)-1, sndbuf, __FILE__,__LINE__)

     ! Prepare room for rcvbuf
     call memalloc (rcv_ptrs_buf(number_neighbours+1)-1, rcvbuf, __FILE__,__LINE__)
     
     call memalloc (number_neighbours, rcvhd, __FILE__,__LINE__)
     call memalloc (number_neighbours, sndhd, __FILE__,__LINE__)

     ! Pack data items into send buffer
     do i=1, number_neighbours
        current = snd_ptrs_buf(i)
        do j=snd_ptrs(i),snd_ptrs(i+1)-1
          call data(snd_leids(j))%pack( elemsizes(snd_leids(j)), sndbuf(current) )
          current = current + elemsizes(snd_leids(j))
        end do
     end do
     
     ! First post all the non blocking receives   
     do i=1, number_neighbours
       proc_to_comm = neighbour_ids(i)
         
       ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
       call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)
          
       ! Message size to be received
       sizmsg = rcv_ptrs_buf(i+1)-rcv_ptrs_buf(i)
      
       if ( (sizmsg > 0) .and. (neighbour_ids(i)-1 /= my_pid) ) then
          call mpi_irecv(  rcvbuf(rcv_ptrs_buf(i)), sizmsg, &
                        &  psb_mpi_integer1, proc_to_comm, &
                        &  psb_double_swap_tag, the_mpi_comm, rcvhd(i), iret)
          check ( iret == mpi_success )
       end if
     end do

     ! Secondly post all non-blocking sends
     do i=1, number_neighbours
        proc_to_comm = neighbour_ids(i)
          
        ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
        call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)
          
        ! Message size to be sent
        sizmsg = snd_ptrs_buf(i+1)-snd_ptrs_buf(i)
    
        if ( (sizmsg > 0) .and. (neighbour_ids(i)-1 /= my_pid) ) then 
             call mpi_isend(sndbuf(snd_ptrs_buf(i)), sizmsg, &
                     & psb_mpi_integer1, proc_to_comm, &
                     & psb_double_swap_tag, the_mpi_comm, sndhd(i), iret)
             check ( iret == mpi_success )
        end if
     end do

     ! Wait on all non-blocking receives
     do i=1, number_neighbours
       proc_to_comm = neighbour_ids(i)
         
       ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
       call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)
          
       ! Message size to be received
       sizmsg = rcv_ptrs_buf(i+1)-rcv_ptrs_buf(i)
      
       if ( (sizmsg > 0) .and. (neighbour_ids(i)-1 /= my_pid) ) then
          call mpi_wait(rcvhd(i), p2pstat, iret)
          check (iret == mpi_success)
       else if ( neighbour_ids(i)-1 == my_pid ) then
          if ( sizmsg /= snd_ptrs_buf(i+1)-snd_ptrs_buf(i) ) then 
             write(0,*) 'Fatal error in single_exchange: mismatch on self sendf', & 
                     & sizmsg, snd_ptrs_buf(i+1)-snd_ptrs_buf(i)
          end if
          rcvbuf( rcv_ptrs_buf(i):rcv_ptrs_buf(i+1)-1 ) = &
                 sndbuf( snd_ptrs_buf(i):snd_ptrs_buf(i+1)-1 )
       end if
     end do

     ! Finally wait on all non-blocking sends
     do i=1, number_neighbours
        proc_to_comm = neighbour_ids(i)
          
        ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
        call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)
          
        ! Message size to be received
        sizmsg = snd_ptrs_buf(i+1)-snd_ptrs_buf(i)
      
        if ( (sizmsg > 0) .and. (neighbour_ids(i)-1 /= my_pid) ) then
          call mpi_wait(sndhd(i), p2pstat, iret)
          check ( iret == mpi_success )
        end if
     end do
     
     ! Unpack data items from recv buffer
     current = rcv_ptrs_buf(1)
     do i=1, number_neighbours
        do j=rcv_ptrs(i),rcv_ptrs(i+1)-1
          call data(rcv_leids(j))%unpack( elemsizes(rcv_leids(j)), rcvbuf(current) )
          current = current + elemsizes(rcv_leids(j))
        end do
     end do

     call memfree (rcvhd ,__FILE__,__LINE__) 
     call memfree (sndhd ,__FILE__,__LINE__)
     call memfree (sndbuf,__FILE__,__LINE__)
     call memfree (rcvbuf,__FILE__,__LINE__)
     call memfree (elemsizes, __FILE__, __LINE__ )
     call memfree (snd_ptrs_buf, __FILE__, __LINE__ )
     call memfree (rcv_ptrs_buf, __FILE__, __LINE__ )
   end subroutine plain_ghost_element_exchange 

end module par_element_exchange_names
