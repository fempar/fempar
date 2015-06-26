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
module par_mesh_partition_conditions_names

  ! Serial modules
use types_names
  use fem_mesh_names
use fem_mesh_gen_names
  use fem_conditions_names
  !use fem_partition_names
  use fem_materials_names
  ! use partition_import
use fem_mesh_gen_names
  !use fem_mesh_io
  !use stdio
  ! use fem_import_names

  ! Parallel modules
  !use par_partition_names
  use par_context_names
  use par_mesh_names
  

# include "debug.i90"
  implicit none
  private

  ! public :: par_structured_create, par_structured_gen_partition

contains

!!$  !=============================================================================
!!$  subroutine par_structured_create (lpart,ndime,isper,nedir,npdir,nsckt,msize,poin,    &
!!$       &                            line,surf,p_context,p_mesh,p_part,nodes,g_context, &
!!$       &                            c_context,mater,mtype)
!!$    !-----------------------------------------------------------------------
!!$    ! 
!!$    !-----------------------------------------------------------------------
!!$    implicit none
!!$    ! Parameters 
!!$    integer(ip),          intent(in)         :: isper(3),nedir(3),npdir(3),nsckt(3)
!!$    integer(ip),          intent(in)         :: ndime,lpart
!!$    type(mesh_size_t),      intent(in)         :: msize
!!$    type(fem_conditions_t), intent(in)         :: poin,line,surf
!!$    type(par_context_t)   , target, intent(in) :: p_context
!!$    type(par_partition) , target, intent(out):: p_part
!!$    type(par_mesh_t),       intent(out)        :: p_mesh
!!$    type(fem_conditions_t), intent(out)        :: nodes
!!$    type(par_context_t)  , target, intent(in), optional  :: g_context
!!$    type(par_context_t)  , target, intent(in), optional  :: c_context
!!$    type(fem_materials_t), optional, intent(inout) :: mater
!!$    integer(ip),         optional, intent(in)    :: mtype
!!$
!!$    ! Locals
!!$    integer         :: iam, num_procs
!!$
!!$    ! Finalize par_partition generation
!!$    assert ( p_context%handler == inhouse )
!!$    p_part%p_context => p_context
!!$    if(present(g_context)) p_part%g_context => g_context
!!$    if(present(c_context)) p_part%c_context => c_context
!!$
!!$    if(p_context%iam>=0) then
!!$       call fem_mesh_gen_partition_create( lpart,ndime,isper,nedir,npdir,nsckt,msize,poin,line,surf, &
!!$            &                                p_mesh%f_mesh,p_part%f_part,nodes,mater,mtype)
!!$
!!$       call par_context_info ( p_context, iam, num_procs )
!!$
!!$       assert ( p_context%handler==inhouse.and.p_part%f_part%ptype==element_based )
!!$       assert ( p_part%f_part%nparts == num_procs )
!!$       
!!$       call partition_to_import ( p_part%f_part, p_part%f_import )
!!$
!!$       ! call fem_import_print (6,p_part%f_import)
!!$    end if
!!$
!!$    assert ( associated(p_part%p_context) )
!!$    p_mesh%p_part => p_part
!!$
!!$    ! Transfer ndime from fine-grid tasks to coarse-grid task
!!$    call par_mesh_bcast (p_mesh, p_mesh%f_mesh%ndime)
!!$
!!$  end subroutine par_structured_create
!!$
!!$   !=============================================================================
!!$   subroutine par_structured_gen_partition (lpart,ndime,isper,nedir,npdir,nsckt,msize,poin,    &
!!$       &                                    line,surf,p_mesh,nodes,mater,mtype)
!!$    !------------------------------------------------------------------------------------------
!!$    ! Here we assume that (pre-conditions):
!!$    !    x p_mesh%p_part has already been fully created via par_partition_create in advance
!!$    !    x p_mesh has been (partially) been created via par_mesh_create. This is why it is
!!$    !      an intent(inout) instead of an intent(out) argument
!!$    !-----------------------------------------------------------------------------------------
!!$    implicit none
!!$    ! Parameters
!!$    integer(ip),          intent(in)              :: lpart, ndime
!!$    integer(ip),          intent(in)              :: isper(3),nedir(3),npdir(3),nsckt(3)
!!$    type(mesh_size_t),      intent(in)              :: msize
!!$    type(fem_conditions_t), intent(in)              :: poin,line,surf
!!$    type(par_mesh_t),       intent(inout)           :: p_mesh
!!$    type(fem_conditions_t), intent(out)             :: nodes
!!$    type(fem_materials_t) , optional, intent(inout) :: mater
!!$    integer(ip),          optional, intent(in)    :: mtype
!!$
!!$
!!$    ! This subroutine requires w_context (all tasks involved here), p_context (fine tasks, which
!!$    ! are in charge of integration in the finner level) and q_context(unused tasks) to be 
!!$    ! created by a call to par_context_create_by_split.
!!$    assert ( associated(p_mesh%p_part%w_context) )
!!$    assert ( p_mesh%p_part%w_context%created .eqv. .true. )
!!$    assert ( associated(p_mesh%p_part%p_context) )
!!$    assert ( p_mesh%p_part%p_context%created .eqv. .true. )
!!$    assert ( associated(p_mesh%p_part%q_context) )
!!$    assert ( p_mesh%p_part%q_context%created .eqv. .true. )
!!$    ! Check appropriate assignment to context: 1) I'm in w_context and 2) I'm in p_context OR in q_context but not in both
!!$    assert ( p_mesh%p_part%w_context%iam >= 0)
!!$    assert ( (p_mesh%p_part%p_context%iam >=0 .and. p_mesh%p_part%q_context%iam < 0) .or. (p_mesh%p_part%p_context%iam < 0 .and. p_mesh%p_part%q_context%iam >= 0))
!!$
!!$    ! Element-based partitioning/inhouse handler required 
!!$    assert ( p_mesh%p_part%f_part%ptype == element_based )
!!$    assert ( p_mesh%p_part%p_context%handler == inhouse  )
!!$
!!$    if(p_mesh%p_part%p_context%iam>=0) then
!!$       call fem_mesh_gen_partition_create( lpart,ndime,isper,nedir,npdir,nsckt,msize,poin,line,surf, &
!!$            &                              p_mesh%f_mesh,p_mesh%p_part%f_part,nodes,mater,mtype)       
!!$       call partition_to_import ( p_mesh%p_part%f_part, p_mesh%p_part%f_import )
!!$    end if
!!$
!!$    ! *** IMPORTANT NOTE: PENDING !!!!
!!$    ! AFM: How do we solve this with the new design of communicators?
!!$    !      In other words, which tasks and for what ndime is needed? 
!!$    ! Transfer ndime from fine-grid tasks to coarse-grid task
!!$    ! call par_mesh_bcast (p_mesh, p_mesh%f_mesh%ndime)
!!$
!!$    ! AFM: TEMPORARILY I WILL FOLLOW THE FOLLOWING DIRTY SOLUTION.
!!$    ! IT WORKS AS ndime is an intent(in) argument that has been
!!$    ! init by all MPI tasks within the calling driver
!!$    p_mesh%f_mesh%ndime = ndime
!!$    
!!$  end subroutine par_structured_gen_partition

end module par_mesh_partition_conditions_names
