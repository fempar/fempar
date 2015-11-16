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

  implicit none
# include "debug.i90"
  private

  ! Parallel context
  type par_context_t
     ! Store that it was created
     logical              :: created = .false.

     ! The following member is an integer which identifies an initialized 
     ! parallel context which is created and manipulated using 
     ! our own wrappers to the MPI library (adapted from PSBLAS 2.4.0 library)

     ! ***IMPORTANT NOTE***: parallel contexts are always of type 
     ! integer: the kind parameter must NOT be specified. This requirement is 
     ! imposed by the underlying message-passing library, i.e., MPI. 
     integer              :: icontxt
     integer              :: iam
     integer              :: np

  end type par_context_t

  interface par_context_create
     module procedure par_context_create_by_list, par_context_create_by_color, par_context_create_bridge
  end interface par_context_create

  ! Types
  public :: par_context_t

  ! Functions
  public :: par_context_create, par_context_null, par_context_info, par_context_free

contains

  !=============================================================================

  subroutine par_context_null ( p_context )
#ifdef MPI_MOD
use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(par_context_t), intent(out) :: p_context

    p_context%created = .true.

    p_context%icontxt = mpi_comm_null
    ! Store current info (who am I and how many we are in these contexts).
    call psb_info ( p_context%icontxt, p_context%iam, p_context%np )
  end subroutine par_context_null

  subroutine par_context_create_by_list ( p_context, nproc, b_context, lproc )
#ifdef MPI_MOD
use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif

    type(par_context_t), intent(out) :: p_context

    integer(ip)      , intent(in), optional :: nproc
    type(par_context_t), intent(in), optional :: b_context  ! Base context
    integer(ip)      , intent(in), optional :: lproc(:)

    ! ***IMPORTANT NOTE***: for the following two variables, 
    ! the kind parameter must NOT be specified. This requirement is 
    ! imposed by the underlying message-passing library, i.e., MPI 
    integer                         :: info
    logical                         :: initialized

    if(p_context%created) then
       write(0,*) 'Error in par_context create,'
       write(0,*) 'it has been already created'
       stop 
    end if
    p_context%created = .true.

    ! Using the complete interface initialize parallel environment
    if(present(b_context)) then
       call psb_init ( p_context%icontxt, nproc, b_context%icontxt, lproc)
    else
       call psb_init ( p_context%icontxt, np=nproc, ids=lproc)
    end if
    ! and store current info (who am I and how many we are in this context).
    call psb_info ( p_context%icontxt, p_context%iam, p_context%np )

  end subroutine par_context_create_by_list

  subroutine par_context_create_by_color ( my_color, p_context, q_context, b_context )
#ifdef MPI_MOD
use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif

    ! Parameters
    integer(ip)      , intent(in)  :: my_color
    type(par_context_t), intent(out) :: p_context
    type(par_context_t), intent(out) :: q_context
    type(par_context_t), intent(in), optional :: b_context  ! Base context

    ! ***IMPORTANT NOTE***: for the following two variables, 
    ! the kind parameter must NOT be specified. This requirement is 
    ! imposed by the underlying message-passing library, i.e., MPI 
    integer                         :: info, my_color_
    integer             , parameter :: key=0
    logical                         :: initialized
    type(par_context_t)               :: b_context_  ! Base context

    if(p_context%created) then
       write(0,*) 'Error in par_context create,'
       write(0,*) 'it has been already created'
       check(.false.)
    end if
    p_context%created = .true.
    q_context%created = .true.


    ! If not given initialize parallel environment
    if(present(b_context)) then
       assert (b_context%created .eqv. .true.)
       b_context_%icontxt = b_context%icontxt
    else
       call psb_init ( b_context_%icontxt)
    end if

    ! Create comm by split
    if(my_color<0) then
       my_color_ = mpi_undefined
       !write(*,*) 'COLOR:',my_color, b_context_%icontxt
       call mpi_comm_split(b_context_%icontxt, my_color_, key , q_context%icontxt, info)
       p_context%icontxt = mpi_comm_null
    else
       my_color_ = my_color
       call mpi_comm_split(b_context_%icontxt, my_color_, key , p_context%icontxt, info)
       q_context%icontxt = mpi_comm_null
    end if
    assert(info==MPI_SUCCESS)

    ! and store current info (who am I and how many we are in these contexts).
    call psb_info ( p_context%icontxt, p_context%iam, p_context%np )
    call psb_info ( q_context%icontxt, q_context%iam, q_context%np )
    !write(*,*) 'I AM:',my_color, p_context%iam, q_context%iam

  end subroutine par_context_create_by_color

  subroutine par_context_create_bridge ( w_context, p_context, q_context, b_context )
#ifdef MPI_MOD
use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif

    ! Parameters
    type(par_context_t), intent(in)  :: w_context
    type(par_context_t), intent(in)  :: p_context
    type(par_context_t), intent(in)  :: q_context
    type(par_context_t), intent(out) :: b_context

    ! ***IMPORTANT NOTE***: for the following two variables, 
    ! the kind parameter must NOT be specified. This requirement is 
    ! imposed by the underlying message-passing library, i.e., MPI 
    integer                         :: info, my_color_
    logical                         :: initialized
    type(par_context_t)               :: b_context_  ! Base context

    if(b_context%created) then
       write(0,*) 'Error in par_context create,'
       write(0,*) 'it has been already created'
       stop 
    end if
    ! We assume w_context to be the union of p_context and q_context. Further
    ! we assume that the ranks ids' (in w_context) of tasks in p_context go from 
    ! 0 to np-1 and ranks ids' of tasks in q context from np onwards
    assert ( w_context%created .eqv. .true.)
    assert ( w_context%iam >= 0)
    assert ( p_context%created .eqv. .true.)
    assert ( q_context%created .eqv. .true.)
    b_context%created = .true.


    ! Create intracomm
    if(p_context%iam>=0) then
       assert(q_context%iam<0)
       !    mpi_intercomm_create (local_comm, local_leader, peer_comm, remote_leader, tag, newintercomm, info)
       call mpi_intercomm_create(p_context%icontxt, 0, w_context%icontxt, p_context%np, 0, b_context%icontxt, info)
    else 
       assert(q_context%iam>=0)
       call mpi_intercomm_create(q_context%icontxt, 0, w_context%icontxt, 0 , 0, b_context%icontxt, info)
    end if
    assert(info==MPI_SUCCESS)
    
    ! and store current info (who am I and how many we are in these contexts).
    call psb_info ( b_context%icontxt, b_context%iam, b_context%np )

  end subroutine par_context_create_bridge

  subroutine par_context_info ( p_context, my_pid, num_procs )
    ! Parameters
    type(par_context_t), intent(in)  :: p_context
    integer          , intent(out) :: my_pid, num_procs

    assert( p_context%created )
    my_pid = p_context%iam
    num_procs = p_context%np
  end subroutine par_context_info

  subroutine par_context_free ( p_context, close  )
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    ! Parameters
    type(par_context_t), intent(inout) :: p_context
    logical, intent(in), optional    :: close


    ! ***IMPORTANT NOTE***: for the following  variable, 
    ! the kind parameter must NOT be specified. This requirement is 
    ! imposed by the underlying message-passing library, i.e., MPI 
    integer :: info
    
    assert( p_context%created )
    
    call psb_exit ( p_context%icontxt, close )
    p_context%created = .false.
    
  end subroutine par_context_free

end module par_context_names
