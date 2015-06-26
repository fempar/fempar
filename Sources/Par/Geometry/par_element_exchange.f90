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
   use fem_element_import_names
   use migratory_element_names

   implicit none
   private

#include "debug.i90"
   
   ! Subroutines
   public :: ghost_elements_exchange
contains

  !*** Abstract implementation of ghost element exchange ***!
  subroutine ghost_elements_exchange ( icontxt, f_el_import, data )
    implicit none
    ! Locals 
    integer(ip)              :: icontxt
    type(fem_element_import_t) :: f_el_import
    ! Data is an array of polymorphic entries     
    class(migratory_element) :: data(f_el_import%nelem + f_el_import%nghost)
   
    call plain_ghost_element_exchange ( icontxt, f_el_import%npadj, f_el_import%lpadj, &
                                        f_el_import%rcv_ptrs, f_el_import%snd_ptrs, f_el_import%snd_leids, &
                                        f_el_import%nelem, f_el_import%nghost, data)

  end subroutine ghost_elements_exchange

  subroutine plain_ghost_element_exchange ( icontxt, npadj, lpadj, &
                                            rcv_ptrs, snd_ptrs, snd_leids, &
                                            nelem, nghost, data)
                                             
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
     integer, intent(in) :: icontxt 
     
     ! **IMPORTANT NOTE**: I will assume that both 
     ! list_rcv and list_snd hold process identifiers 
     ! starting from one. However, the underlying
     ! message-passing library identifies processes 
     ! starting from zero ...

     ! Control info to receive
     integer    , intent(in)    :: npadj, lpadj(npadj), rcv_ptrs(npadj+1)

     ! Control info to send
     integer    , intent(in)    :: snd_ptrs(npadj+1)
     integer(ip), intent(in)    :: snd_leids(snd_ptrs(npadj+1)-1)

     integer(ip), intent(in)    :: nelem, nghost
     class(migratory_element), intent(inout) :: data(nelem+nghost) 
     
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
     call data(1)%size(elemsize)

     ! Prepare room for sndbuf
     call memalloc ((snd_ptrs(npadj+1)-snd_ptrs(1))*elemsize, sndbuf, __FILE__,__LINE__)

     ! Prepare room for rcvbuf
     call memalloc ((rcv_ptrs(npadj+1)-rcv_ptrs(1))*elemsize, rcvbuf, __FILE__,__LINE__)

     ! Pack data items into send buffer
     do i=1, snd_ptrs(npadj+1)-1
        call data(snd_leids(i))%pack(elemsize,sndbuf((i-1)*elemsize+1)) 
     end do

     call memalloc (npadj, rcvhd, __FILE__,__LINE__)
     call memalloc (npadj, sndhd, __FILE__,__LINE__)

     ! First post all the non blocking receives   
     do i=1, npadj
       proc_to_comm = lpadj(i)
         
       ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
       call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)
          
       ! Message size to be received
       sizmsg = (rcv_ptrs(i+1)-rcv_ptrs(i))*elemsize
      
       if ( (sizmsg > 0) .and. (lpadj(i)-1 /= my_pid) ) then
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
     do i=1, npadj
        proc_to_comm = lpadj(i)
          
        ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
        call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)
          
        ! Message size to be sent
        sizmsg = (snd_ptrs(i+1)-snd_ptrs(i))*elemsize
    
        if ( (sizmsg > 0) .and. (lpadj(i)-1 /= my_pid) ) then 
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
     do i=1, npadj
       proc_to_comm = lpadj(i)
         
       ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
       call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)
          
       ! Message size to be received
       sizmsg = (rcv_ptrs(i+1)-rcv_ptrs(i))*elemsize
      
       if ( (sizmsg > 0) .and. (lpadj(i)-1 /= my_pid) ) then
          call mpi_wait(rcvhd(i), p2pstat, iret)

          if ( iret /= mpi_success ) then
             write (0,*) 'Error: mpi_wait returned != mpi_success'
             call psb_abort (icontxt)    
          end if
       else if ( lpadj(i)-1 == my_pid ) then
          if ( sizmsg /= (snd_ptrs(i+1)-snd_ptrs(i))*elemsize ) then 
             write(0,*) 'Fatal error in single_exchange: mismatch on self sendf', & 
                     & sizmsg, (snd_ptrs(i+1)-snd_ptrs(i) )*elemsize
          end if

          rcvbuf( (rcv_ptrs(i)-1)*elemsize+1:(rcv_ptrs(i)-1)*elemsize+sizmsg) = &
                 sndbuf( (snd_ptrs(i)-1)*elemsize+1:(snd_ptrs(i)-1)*elemsize+sizmsg )
       end if
    end do

     ! Finally wait on all non-blocking sends
     do i=1, npadj
        proc_to_comm = lpadj(i)
          
        ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
        call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)
          
        ! Message size to be received
        sizmsg = (snd_ptrs(i+1)-snd_ptrs(i))*elemsize
      
        if ( (sizmsg > 0) .and. (lpadj(i)-1 /= my_pid) ) then
          call mpi_wait(sndhd(i), p2pstat, iret)
          if ( iret /= mpi_success ) then
              write (0,*) 'Error: mpi_wait returned != mpi_success'
              call psb_abort (icontxt)    
          end if
        end if
     end do

     ! Unpack data items into send buffer
     do i=1, rcv_ptrs(npadj+1)-1
        call data(nelem+i)%unpack(elemsize,rcvbuf((i-1)*elemsize+1)) 
     end do

     call memfree (rcvhd ,__FILE__,__LINE__) 
     call memfree (sndhd ,__FILE__,__LINE__)
     call memfree (sndbuf,__FILE__,__LINE__)
     call memfree (rcvbuf,__FILE__,__LINE__)
   
   end subroutine plain_ghost_element_exchange 

end module par_element_exchange_names
