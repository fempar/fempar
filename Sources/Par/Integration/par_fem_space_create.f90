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
# include "debug.i90"
module par_fem_space_names

  ! Fem Modules
  use types
  use memor
  use fem_space_names
  use fem_space_types
  use dof_handler_names
  use hash_table_names

  ! Par Modules
  use par_triangulation_names
  use par_element_exchange

#ifdef memcheck
  use iso_c_binding
#endif
  implicit none
# include "debug.i90"
  private

  public :: par_fem_space_create 

contains

  subroutine par_fem_space_create ( p_trian, dhand, femsp, problem, continuity, order, material, &
       & time_steps_to_store, hierarchical_basis, static_condensation, num_materials  )
    implicit none
    type(par_triangulation), intent(inout) :: p_trian
    type(dof_handler), intent(in)       :: dhand
    type(fem_space), intent(inout)      :: femsp  
    integer(ip),     intent(in)       :: material(:), order(:,:), problem(:)
    logical(lg),     intent(in)       :: continuity(:,:)
    integer(ip), optional, intent(in) :: time_steps_to_store
    logical(lg), optional, intent(in) :: hierarchical_basis
    logical(lg), optional, intent(in) :: static_condensation
    integer(ip), optional, intent(in) :: num_materials 

    integer(ip) :: num_elems, num_ghosts, ielem

    ! Create local fem space
    num_elems = p_trian%f_trian%num_elems
    num_ghosts = p_trian%num_ghosts

    p_trian%f_trian%num_elems = num_elems + num_ghosts
    write(*,*) '***** CREATE FEM SPACE STRUCTURES *****'
    call fem_space_allocate_structures( p_trian%f_trian, dhand, femsp, &
         time_steps_to_store = time_steps_to_store, hierarchical_basis = hierarchical_basis, &
         static_condensation = static_condensation, num_materials = num_materials )
    p_trian%f_trian%num_elems = num_elems
    write(*,*) '***** FILL FEM SPACE STRUCTURES *****'
    call fem_space_fe_list_create( femsp, problem, continuity, order, material )
    

    ! Communicate problem, continuity, order, and material
    write(*,*) '***** EXCHANGE GHOST INFO *****'
    call ghost_elements_exchange ( p_trian%p_context%icontxt, p_trian%f_el_import, femsp%lelem )

    ! Create ghost fem space (only partially, i.e., previous info)
    write(*,*) '***** FILL GHOST ELEMENTS  *****'

    
    !write(*,*) 'num_elems+1', femsp%g_trian%num_elems+1
    !write(*,*) 'num_ghosts', num_ghosts
    !do ielem = femsp%g_trian%num_elems+1, femsp%g_trian%num_elems+num_ghosts
    !   call fem_element_fixed_info_write( p_trian%f_trian%elems(ielem)%topology )
    !end do

    !write(*,*) 'num_elems+1', femsp%g_trian%num_elems+1
    !write(*,*) 'num_ghosts', num_ghosts
    !do ielem = femsp%g_trian%num_elems+1, femsp%g_trian%num_elems+num_ghosts
    !   call fem_element_fixed_info_write( femsp%g_trian%elems(ielem)%topology )
    !end do

    
    call ghost_fe_list_create( femsp, num_ghosts )

  end subroutine par_fem_space_create

  subroutine ghost_fe_list_create( femsp, num_ghosts )
    implicit none
    type(fem_space), intent(inout), target :: femsp
    integer(ip), intent(in)        :: num_ghosts

    integer(ip) :: ielem, nvars, f_type, ivar, f_order, istat, pos_elinf, v_key
    logical(lg) :: created
    
    do ielem = femsp%g_trian%num_elems+1, femsp%g_trian%num_elems+num_ghosts
       write (*,*) '************* GHOST ELEMENT *************',ielem
       nvars = femsp%dof_handler%problems(femsp%lelem(ielem)%problem)%nvars
       femsp%lelem(ielem)%num_vars = nvars
       f_type = femsp%g_trian%elems(ielem)%topology%ftype
       write(*,*) 'f_type ghosts',f_type
       write(*,*) 'nvars',nvars
       assert ( f_type > 0)
       call memalloc(nvars, femsp%lelem(ielem)%f_inf, __FILE__, __LINE__ )
       do ivar=1,nvars
          f_order = femsp%lelem(ielem)%order(ivar)
          write(*,*) 'f_order',f_order
          v_key = femsp%g_trian%num_dims + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)*f_order
          call femsp%ht_elem_info%put(key=v_key,val=femsp%cur_elinf,stat=istat)
          if ( istat == now_stored) then 
             write (*,*) ' FIXED INFO NEW'
             call fem_element_fixed_info_create(femsp%lelem_info(femsp%cur_elinf),f_type,              &
                  &                             f_order,femsp%g_trian%num_dims,created)
             assert(created)
             pos_elinf = femsp%cur_elinf
             femsp%cur_elinf = femsp%cur_elinf + 1
          else if ( istat == was_stored ) then
             write (*,*) ' FIXED INFO ALREADY STORED '
             call femsp%ht_elem_info%get(key=v_key,val=pos_elinf,stat=istat)
             assert ( istat == key_found )
          end if
          femsp%lelem(ielem)%f_inf(ivar)%p => femsp%lelem_info(pos_elinf)
       end do
    end do
  end subroutine ghost_fe_list_create
end module par_fem_space_names

