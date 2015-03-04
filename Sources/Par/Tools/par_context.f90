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
  use types

  ! Trilinos-Interfaces modules 
  !use for_trilinos_shadow_interfaces

  ! Parallel modules
  use psb_penv_mod

# include "debug.i90"
  implicit none
  private

  ! ***IMPORTANT NOTE***: I am assuming that the
  ! constructor of par_context_names is responsible
  ! for creating/destroying all the members/objects it 
  ! contains, INCLUDING initialization/finalization of
  ! the parallel environment
  integer(ip), parameter :: inhouse  = 0
  integer(ip), parameter :: trilinos = 1
  integer(ip), parameter :: petsc    = 2

  ! Parallel context
  type par_context
     ! Store that it was created
     logical              :: created = .false.

     ! Epetra Communicator. 
     ! This communicator handles the parallel environment and 
     ! is required for creating epetra_maps
     !type(epetra_mpicomm) :: epcomm

     ! The following member is an integer which identifies an initialized 
     ! parallel context which is created and manipulated using 
     ! our own wrappers to the MPI library (adapted from PSBLAS 2.4.0 library)

     ! ***IMPORTANT NOTE***: parallel contexts are always of type 
     ! integer: the kind parameter must NOT be specified. This requirement is 
     ! imposed by the underlying message-passing library, i.e., MPI. 
     integer              :: icontxt
     integer              :: iam
     integer              :: np

     ! Parallel Software responsible for performing 
     ! parallel computations
     integer(ip)          :: handler
  end type par_context

  interface par_context_create
     module procedure par_context_create_by_list, par_context_create_by_color, par_context_create_bridge
  end interface par_context_create

  ! Types
  public :: par_context

  ! Constants
  public :: inhouse, trilinos, petsc

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
    type(par_context), intent(out) :: p_context

    p_context%created = .true.
    p_context%handler = inhouse

    p_context%icontxt = mpi_comm_null
    ! Store current info (who am I and how many we are in these contexts).
    call psb_info ( p_context%icontxt, p_context%iam, p_context%np )
  end subroutine par_context_null

  subroutine par_context_create_by_list ( handler, p_context, nproc, b_context, lproc )
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif

    ! Parameters
    integer(ip)      , intent(in)  :: handler
    type(par_context), intent(out) :: p_context

    integer(ip)      , intent(in), optional :: nproc
    type(par_context), intent(in), optional :: b_context  ! Base context
    integer(ip)      , intent(in), optional :: lproc(:)

    ! ***IMPORTANT NOTE***: for the following two variables, 
    ! the kind parameter must NOT be specified. This requirement is 
    ! imposed by the underlying message-passing library, i.e., MPI 
    integer                         :: info
    logical                         :: initialized

    assert ( handler == inhouse .or. handler == trilinos )
    if(p_context%created) then
       write(0,*) 'Error in par_context create,'
       write(0,*) 'it has been already created'
       stop 
    end if
    p_context%created = .true.
    p_context%handler = handler

    if ( handler == inhouse ) then
       ! Initialize psb parallel environment and 
       ! we are done
       !call psb_init ( p_context%icontxt )

       ! Using the complete interface initialize parallel environment
       if(present(b_context)) then
          call psb_init ( p_context%icontxt, nproc, b_context%icontxt, lproc)
       else
          call psb_init ( p_context%icontxt, np=nproc, ids=lproc)
       end if
       ! and store current info (who am I and how many we are in this context).
       call psb_info ( p_context%icontxt, p_context%iam, p_context%np )

    else 
       if ( handler == trilinos ) then

          ! Not ready to handle multiple contexts
          if(present(nproc).or.present(b_context).or.present(lproc)) then
             write(0,*) 'Error in par_context create:'
             write(0,*) 'not ready to handle multiple contexts with trilinos'
             stop 
          end if

          ! Initialize mpi parallel environment
          call mpi_initialized(initialized, info)
          if ((.not.initialized).or.(info /= mpi_success)) then 
             !call mpi_init(info) 

             ! *** VERY IMPORTANT NOTE ***: The following two lines are 
             ! required instead of mpi_init above in order to properly
             ! use par_timer when using Trilinos handler. This is
             ! dirty as PSB communication layer is intended to use
             ! only with inhouse handler. A cleaner solution is required here,
             ! with an implementation of par_timer which takes into
             ! account handler and uses corresponding communication routines
             ! accordingly
             call psb_init ( p_context%icontxt )
             info = mpi_success

             if (info /= mpi_success) then
                write(0,*) 'Error in initalizing MPI, bailing out',info 
                stop 
             end if
             ! Construct Epetra_mpicomm object
             !call epetra_mpicomm_construct ( p_context%epcomm )

             ! Store who am I and how many we are
             !call epetra_mpicomm_my_pid   (p_context%epcomm, p_context%iam)
             !call epetra_mpicomm_num_proc (p_context%epcomm, p_context%np)

          end if
       end if
    end if

  end subroutine par_context_create_by_list

  subroutine par_context_create_by_color ( handler, my_color, p_context, q_context, b_context )
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif

    ! Parameters
    integer(ip)      , intent(in)  :: handler
    integer(ip)      , intent(in)  :: my_color
    type(par_context), intent(out) :: p_context
    type(par_context), intent(out) :: q_context
    type(par_context), intent(in), optional :: b_context  ! Base context

    ! ***IMPORTANT NOTE***: for the following two variables, 
    ! the kind parameter must NOT be specified. This requirement is 
    ! imposed by the underlying message-passing library, i.e., MPI 
    integer                         :: info, my_color_
    integer             , parameter :: key=0
    logical                         :: initialized
    type(par_context)               :: b_context_  ! Base context

    assert ( handler == inhouse .or. handler == trilinos )
    if(p_context%created) then
       write(0,*) 'Error in par_context create,'
       write(0,*) 'it has been already created'
       stop 
    end if
    p_context%created = .true.
    p_context%handler = handler
    q_context%created = .true.
    q_context%handler = handler

    if ( handler == inhouse ) then

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
          !write(*,*) 'COLOR:',my_color, b_context_%icontxt
          call mpi_comm_split(b_context_%icontxt, my_color_, key , p_context%icontxt, info)
          q_context%icontxt = mpi_comm_null
       end if
       assert(info==MPI_SUCCESS)

       ! and store current info (who am I and how many we are in these contexts).
       call psb_info ( p_context%icontxt, p_context%iam, p_context%np )
       call psb_info ( q_context%icontxt, q_context%iam, q_context%np )
       !write(*,*) 'I AM:',my_color, p_context%iam, q_context%iam

    else if ( handler == trilinos ) then

       ! Not ready to handle multiple contexts
       write(0,*) 'Error in par_context create:'
       write(0,*) 'not ready to handle multiple contexts with trilinos'
       stop 

    end if

  end subroutine par_context_create_by_color

  subroutine par_context_create_bridge ( handler, w_context, p_context, q_context, b_context )
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif

    ! Parameters
    integer(ip)      , intent(in)  :: handler
    type(par_context), intent(in)  :: w_context
    type(par_context), intent(in)  :: p_context
    type(par_context), intent(in)  :: q_context
    type(par_context), intent(out) :: b_context

    ! ***IMPORTANT NOTE***: for the following two variables, 
    ! the kind parameter must NOT be specified. This requirement is 
    ! imposed by the underlying message-passing library, i.e., MPI 
    integer                         :: info, my_color_
    logical                         :: initialized
    type(par_context)               :: b_context_  ! Base context

    assert ( handler == inhouse .or. handler == trilinos )
    if(b_context%created) then
       write(0,*) 'Error in par_context create,'
       write(0,*) 'it has been already created'
       stop 
    end if
    ! We assume w_context to be the union of p_context and q_context. Further
    ! we assume that the ranks ids' (in w_context) of tasks in p_context go from 
    ! 0 to np-1 and ranks ids' of tasks in q context from np onwards
    assert ( w_context%created .eqv. .true.)
    assert ( w_context%handler == inhouse )
    assert ( w_context%iam >= 0)
    assert ( p_context%created .eqv. .true.)
    assert ( p_context%handler == inhouse )
    assert ( q_context%created .eqv. .true.)
    assert ( q_context%handler == inhouse )
    b_context%created = .true.
    b_context%handler = handler

    if ( handler == inhouse ) then

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
 
    else if ( handler == trilinos ) then

       ! Not ready to handle multiple contexts
       write(0,*) 'Error in par_context create:'
       write(0,*) 'not ready to handle multiple contexts with trilinos'
       stop 

    end if

  end subroutine par_context_create_bridge

  subroutine par_context_info ( p_context, my_pid, num_procs )
    ! Parameters
    type(par_context), intent(in)  :: p_context
    integer          , intent(out) :: my_pid, num_procs

    assert(p_context%created .eqv. .true.)
    my_pid = p_context%iam
    num_procs = p_context%np

    ! assert ( p_context%handler == inhouse .or. p_context%handler == trilinos )

    ! if ( p_context%handler == inhouse ) then
    !    call psb_info ( p_context%icontxt, my_pid, num_procs )
    ! else 
    !    if ( p_context%handler == trilinos ) then
    !       call epetra_mpicomm_my_pid   (p_context%epcomm, my_pid)
    !       call epetra_mpicomm_num_proc (p_context%epcomm, num_procs)
    !    end if
    ! end if

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
    type(par_context), intent(inout) :: p_context
    logical, intent(in), optional    :: close


    ! ***IMPORTANT NOTE***: for the following  variable, 
    ! the kind parameter must NOT be specified. This requirement is 
    ! imposed by the underlying message-passing library, i.e., MPI 
    integer :: info

    assert(p_context%created .eqv. .true.)
    assert(p_context%handler==inhouse.or.p_context%handler==trilinos)

    if ( p_context%handler == inhouse ) then
       ! Finalize psb parallel environment and we are done
       call psb_exit ( p_context%icontxt, close )
    else 
       if ( p_context%handler == trilinos ) then
          ! Destruct epetra_mpicomm object
          !call epetra_mpicomm_destruct ( p_context%epcomm )

          ! *** VERY IMPORTANT NOTE ***: The following two lines are 
          ! required instead of mpi_init above in order to properly
          ! use par_timer when using Trilinos handler. This is
          ! dirty as PSB communication layer is intended to use
          ! only with inhouse handler. A cleaner solution is required here,
          ! with an implementation of par_timer which takes into
          ! account handler and uses corresponding communication routines
          ! accordingly

          ! Finalize mpi parallel environment
          !call mpi_finalize(info)
          call psb_exit ( p_context%icontxt, .true. )

       end if
    end if
    p_context%created = .false.

  end subroutine par_context_free

end module par_context_names
