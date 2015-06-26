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
module par_environment_names
  ! Serial modules
  use types_names
  use memor_names
  use maps_names
  use abstract_environment_names

  ! Parallel modules
  use psb_penv_mod_names
  use par_context_names

# include "debug.i90"
  implicit none
  private

  type, extends(abstract_environment) ::  par_environment_t
     logical                      :: created             ! Has the parallel environment been created?
     type (par_context_t), pointer  :: p_context => NULL() ! Fine process
     type (par_context_t), pointer  :: q_context => NULL() ! Available (unused) processes 
     type (par_context_t), pointer  :: b_context => NULL() ! Intercommunicator betwen p_context and q_context (bcast and recursive call)
     type (par_context_t), pointer  :: w_context => NULL() ! World communicator (all process involved).
     
     ! Number of levels in the multilevel hierarchy of MPI tasks
     integer(ip) :: num_levels = -1 
     
     integer(ip), allocatable ::     &
          id_parts(:),               &    ! Part identifier
          num_parts(:)                    ! Number of parts

   contains
     
     procedure :: par_environment_create_single_level 
     procedure :: par_environment_create_multilevel
     generic :: create =>  par_environment_create_single_level, &
                           par_environment_create_multilevel
     
     procedure :: par_environment_bcast_logical

     
     procedure :: info                => par_environment_info
     procedure :: am_i_fine_task      => par_environment_am_i_fine_task
     procedure :: bcast               => par_environment_bcast_logical
     procedure :: first_level_barrier => par_environment_first_level_barrier
  end type par_environment_t
  

  interface par_environment_create
     module procedure par_environment_create_single_level, par_environment_create_multilevel
  end interface par_environment_create

  ! Types
  public :: par_environment_t, par_environment_create, par_environment_free

contains

  subroutine par_environment_first_level_barrier(env) 
    implicit none
    ! Dummy arguments
    class(par_environment_t),intent(in)  :: env

    ! Local variables
    integer :: mpi_comm_p, ierr

    ! Parallel environment MUST BE already created
    assert ( env%created )
    assert ( associated(env%p_context) )

    ! Get MPI communicator associated to icontxt_b (in
    ! the current implementation of our wrappers
    ! to the MPI library icontxt and mpi_comm are actually 
    ! the same)
    call psb_get_mpicomm (env%p_context%icontxt, mpi_comm_p)
    call mpi_barrier ( mpi_comm_p, ierr)
    check ( ierr == 0 )
  end subroutine par_environment_first_level_barrier

  subroutine par_environment_info(env,me,np) 
    implicit none
    class(par_environment_t),intent(in)  :: env
    integer(ip)           ,intent(out) :: me
    integer(ip)           ,intent(out) :: np

    ! Parallel environment MUST BE already created
    assert ( env%created )

    me = env%p_context%iam 
    np = env%p_context%np

  end subroutine par_environment_info
  
  function par_environment_am_i_fine_task(env) 
    implicit none
    class(par_environment_t) ,intent(in)  :: env
    logical                             :: par_environment_am_i_fine_task 

    ! Parallel environment MUST BE already created
    assert ( env%created )

    par_environment_am_i_fine_task = (env%p_context%iam >= 0)

  end function par_environment_am_i_fine_task

  subroutine par_environment_bcast_logical(env,condition)
#ifdef MPI_MOD
use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    ! Parameters
    class(par_environment_t), intent(in)    :: env
    logical               , intent(inout) :: condition

    ! Locals
    integer :: mpi_comm_b, info

    assert ( env%created )

    if ( env%num_levels > 1 ) then
       ! b_context is an intercomm among p_context & q_context
       ! Therefore the semantics of the mpi_bcast subroutine slightly changes
       ! P_0 in p_context is responsible for bcasting condition to all the processes
       ! in q_context

       ! Get MPI communicator associated to icontxt_b (in
       ! the current implementation of our wrappers
       ! to the MPI library icontxt and mpi_comm are actually 
       ! the same)
       call psb_get_mpicomm (env%b_context%icontxt, mpi_comm_b)
       
       if (env%p_context%iam >=0) then
          if ( env%p_context%iam == psb_root_ ) then
             call mpi_bcast(condition,1,MPI_LOGICAL,MPI_ROOT,mpi_comm_b,info)
             check( info == mpi_success )
          else
             call mpi_bcast(condition,1,MPI_LOGICAL,MPI_PROC_NULL,mpi_comm_b,info)
             check( info == mpi_success )
          end if
       else if (env%q_context%iam >=0) then
          call mpi_bcast(condition,1,MPI_LOGICAL,psb_root_,mpi_comm_b,info)
          check( info == mpi_success )
       end if
    end if

  end subroutine par_environment_bcast_logical

  !=============================================================================
  subroutine par_environment_create_single_level ( p_env, p_context )
    implicit none 
    ! Parameters
    class(par_environment_t)        , intent(out) :: p_env
    type(par_context_t)    , target , intent(in)  :: p_context

    assert(p_context%created)
    p_env%p_context => p_context

    nullify(p_env%w_context)
    nullify(p_env%q_context)
    nullify(p_env%b_context)

    p_env%num_levels = 1 
    call memalloc(p_env%num_levels, p_env%id_parts, __FILE__, __LINE__ )
    p_env%id_parts(1) = p_context%iam + 1 

    call memalloc(p_env%num_levels, p_env%num_parts, __FILE__, __LINE__ )
    p_env%num_parts(1) = p_context%np

    p_env%created = .true. 
  end subroutine par_environment_create_single_level
  
  !=============================================================================
  subroutine par_environment_create_multilevel ( p_env,      & ! Parallel environment
                                                 w_context,  & ! Intracomm including p & q
                                                 p_context,  & ! Fine-tasks intracomm
                                                 q_context,  & ! Rest-of-tasks intracomm
                                                 b_context,  & ! Intercomm p<=>q
                                                 num_levels, &
                                                 id_parts,   &
                                                 num_parts )
    implicit none 
    ! Parameters
    class(par_environment_t)   , intent(out)   :: p_env
    type(par_context_t), target, intent(in)    :: w_context
    type(par_context_t), target, intent(in)    :: p_context
    type(par_context_t), target, intent(in)    :: q_context
    type(par_context_t), target, intent(in)    :: b_context
    integer(ip)              , intent(in)    :: num_levels
    integer(ip)              , intent(in)    :: id_parts(:)
    integer(ip)              , intent(in)    :: num_parts(:)

    assert(p_context%created .eqv. .true.)
    p_env%p_context => p_context

    assert(w_context%created .eqv. .true.)
    p_env%w_context => w_context

    assert(q_context%created .eqv. .true.)
    p_env%q_context => q_context

    assert(b_context%created .eqv. .true.)
    p_env%b_context => b_context

    assert(size(id_parts) == num_levels) 
    assert(size(num_parts) == num_levels)
    p_env%num_levels = num_levels
    call memalloc(p_env%num_levels, p_env%id_parts,__FILE__,__LINE__ )
    call memalloc(p_env%num_levels, p_env%num_parts,__FILE__,__LINE__ )
    p_env%id_parts = id_parts
    p_env%num_parts = num_parts

    p_env%created = .true. 
  end subroutine par_environment_create_multilevel

  !=============================================================================
  subroutine par_environment_free ( p_env )
    implicit none 
    ! Parameters
    type(par_environment_t), intent(inout) :: p_env

    ! Parallel environment MUST BE already created
    assert ( p_env%created .eqv. .true.)

    call memfree(p_env%id_parts , __FILE__, __LINE__ )
    call memfree(p_env%num_parts, __FILE__, __LINE__ )

    nullify ( p_env%w_context ) 
    nullify ( p_env%p_context )
    nullify ( p_env%q_context ) 
    nullify ( p_env%b_context )

    p_env%created = .false. 
  end subroutine par_environment_free

end module par_environment_names
