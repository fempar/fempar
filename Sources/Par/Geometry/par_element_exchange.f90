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
   use element_import_names
   use migratory_element_names

   implicit none
#include "debug.i90"
   private

   ! Subroutines
   public :: ghost_elements_exchange
contains

  subroutine ghost_elements_exchange ( icontxt, element_import, data )
    implicit none
    integer(ip)               , intent(in)    :: icontxt
    type(element_import_t)    , intent(in)    :: element_import
    class(migratory_element_t), intent(inout) :: data(:)
   
    call plain_ghost_element_exchange ( icontxt, &
                                        element_import%get_number_neighbours(), &
                                        element_import%get_neighbours_ids(), &
                                        element_import%get_rcv_ptrs(), &
                                        element_import%get_rcv_leids(), &
                                        element_import%get_snd_ptrs(), &
                                        element_import%get_snd_leids(), &
                                        data)
  end subroutine ghost_elements_exchange

  subroutine plain_ghost_element_exchange ( icontxt, &
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
     integer(ip)               , intent(in)    :: icontxt 
     integer(ip)               , intent(in)    :: number_neighbours
     integer(ip)               , intent(in)    :: neighbour_ids(number_neighbours)
     integer(ip)               , intent(in)    :: rcv_ptrs(number_neighbours+1)
     integer(ip)               , intent(in)    :: rcv_leids(rcv_ptrs(number_neighbours+1)-1)
     integer(ip)               , intent(in)    :: snd_ptrs(number_neighbours+1)
     integer(ip)               , intent(in)    :: snd_leids(snd_ptrs(number_neighbours+1)-1)   
     class(migratory_element_t), intent(inout) :: data(:) 
     
     ! Communication related locals 
     integer :: my_pid, num_procs, i, proc_to_comm, sizmsg
     integer :: mpi_comm,  iret, info
     integer :: p2pstat(mpi_status_size)

     ! Request handlers for non-blocking receives
     integer, allocatable, dimension(:) :: rcvhd

     ! Request handlers for non-blocking receives
     integer, allocatable, dimension(:) :: sndhd

     integer(ieep), allocatable :: sndbuf(:)  
     integer(ieep), allocatable :: rcvbuf(:)

     integer(ip) :: elemsize
 
     ! Get my process identifier and number of
     ! of processors available in the parallel
     ! context
     call psb_info(icontxt, my_pid, num_procs)
     
     ! Get MPI communicator associated to icontxt (in
     ! the current implementation of our wrappers
     ! to the MPI library icontxt and mpi_comm are actually 
     ! the same)
     call psb_get_mpicomm (icontxt, mpi_comm)

     ! Get element size
     call data(snd_leids(1))%size(elemsize)

     ! Prepare room for sndbuf
     call memalloc ((snd_ptrs(number_neighbours+1)-snd_ptrs(1))*elemsize, sndbuf, __FILE__,__LINE__)

     ! Prepare room for rcvbuf
     call memalloc ((rcv_ptrs(number_neighbours+1)-rcv_ptrs(1))*elemsize, rcvbuf, __FILE__,__LINE__)

     ! Pack data items into send buffer
     do i=1, snd_ptrs(number_neighbours+1)-1
        call data(snd_leids(i))%pack(elemsize,sndbuf((i-1)*elemsize+1)) 
     end do

     call memalloc (number_neighbours, rcvhd, __FILE__,__LINE__)
     call memalloc (number_neighbours, sndhd, __FILE__,__LINE__)

     ! First post all the non blocking receives   
     do i=1, number_neighbours
       proc_to_comm = neighbour_ids(i)
         
       ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
       call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)
          
       ! Message size to be received
       sizmsg = (rcv_ptrs(i+1)-rcv_ptrs(i))*elemsize
      
       if ( (sizmsg > 0) .and. (neighbour_ids(i)-1 /= my_pid) ) then
          call mpi_irecv(  rcvbuf((rcv_ptrs(i)-1)*elemsize+1), sizmsg, &
                        &  psb_mpi_integer1, proc_to_comm, &
                        &  psb_double_swap_tag, mpi_comm, rcvhd(i), iret)

          if ( iret /= mpi_success ) then
             write (0,*) 'Error: mpi_irecv returned != mpi_success'
             call psb_abort (icontxt)    
          end if
       end if
     end do

     ! Secondly post all non-blocking sends
     do i=1, number_neighbours
        proc_to_comm = neighbour_ids(i)
          
        ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
        call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)
          
        ! Message size to be sent
        sizmsg = (snd_ptrs(i+1)-snd_ptrs(i))*elemsize
    
        if ( (sizmsg > 0) .and. (neighbour_ids(i)-1 /= my_pid) ) then 
             call mpi_isend(sndbuf((snd_ptrs(i)-1)*elemsize+1), sizmsg, &
                     & psb_mpi_integer1, proc_to_comm, &
                     & psb_double_swap_tag, mpi_comm, sndhd(i), iret)

             if ( iret /= mpi_success ) then
                write (0,*) 'Error: mpi_isend returned != mpi_success'
                call psb_abort (icontxt)    
             end if
          end if
       end do

     ! Wait on all non-blocking receives
     do i=1, number_neighbours
       proc_to_comm = neighbour_ids(i)
         
       ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
       call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)
          
       ! Message size to be received
       sizmsg = (rcv_ptrs(i+1)-rcv_ptrs(i))*elemsize
      
       if ( (sizmsg > 0) .and. (neighbour_ids(i)-1 /= my_pid) ) then
          call mpi_wait(rcvhd(i), p2pstat, iret)

          if ( iret /= mpi_success ) then
             write (0,*) 'Error: mpi_wait returned != mpi_success'
             call psb_abort (icontxt)    
          end if
       else if ( neighbour_ids(i)-1 == my_pid ) then
          if ( sizmsg /= (snd_ptrs(i+1)-snd_ptrs(i))*elemsize ) then 
             write(0,*) 'Fatal error in single_exchange: mismatch on self sendf', & 
                     & sizmsg, (snd_ptrs(i+1)-snd_ptrs(i) )*elemsize
          end if

          rcvbuf( (rcv_ptrs(i)-1)*elemsize+1:(rcv_ptrs(i)-1)*elemsize+sizmsg) = &
                 sndbuf( (snd_ptrs(i)-1)*elemsize+1:(snd_ptrs(i)-1)*elemsize+sizmsg )
       end if
    end do

     ! Finally wait on all non-blocking sends
     do i=1, number_neighbours
        proc_to_comm = neighbour_ids(i)
          
        ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
        call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)
          
        ! Message size to be received
        sizmsg = (snd_ptrs(i+1)-snd_ptrs(i))*elemsize
      
        if ( (sizmsg > 0) .and. (neighbour_ids(i)-1 /= my_pid) ) then
          call mpi_wait(sndhd(i), p2pstat, iret)
          if ( iret /= mpi_success ) then
              write (0,*) 'Error: mpi_wait returned != mpi_success'
              call psb_abort (icontxt)    
          end if
        end if
     end do

     ! Unpack data items into send buffer
     do i=1, rcv_ptrs(number_neighbours+1)-1
        call data(rcv_leids(i))%unpack(elemsize,rcvbuf((i-1)*elemsize+1)) 
     end do

     call memfree (rcvhd ,__FILE__,__LINE__) 
     call memfree (sndhd ,__FILE__,__LINE__)
     call memfree (sndbuf,__FILE__,__LINE__)
     call memfree (rcvbuf,__FILE__,__LINE__)
   
   end subroutine plain_ghost_element_exchange 

end module par_element_exchange_names
