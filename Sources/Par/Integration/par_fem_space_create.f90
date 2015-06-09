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
  use problem_names
  use dof_handler_names
  use hash_table_names
  use fem_conditions_names

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


  !*********************************************************************************
  ! This subroutine is intended to be called from a parallel driver, having as input
  ! a set of arrays of size number of local elements times number of global variables
  ! for continuity and order, and number of local elements for material and problem,
  ! together with some optional flags. The output of this subroutine is a fem_space
  ! with the required info on ghost elements, together with the dof generation and the
  ! dof graph distribution.
  !*********************************************************************************
  subroutine par_fem_space_create ( p_trian, dhand, femsp, problem, approximations, bcond, continuity, order, material, &
       & which_approx, num_approximations, time_steps_to_store, hierarchical_basis, static_condensation, num_continuity  )
    implicit none
    type(par_triangulation), intent(inout) :: p_trian
    type(dof_handler), intent(in)       :: dhand
    type(fem_space), intent(inout)      :: femsp  
    type(discrete_problem_pointer) , intent(in)    :: approximations(:)
    type(fem_conditions)           , intent(in)    :: bcond
    integer(ip),     intent(in)       :: material(:), order(:,:), problem(:)
    integer(ip),     intent(in)       :: continuity(:,:), which_approx(:)
    integer(ip), intent(in) :: num_approximations
    integer(ip), optional, intent(in) :: time_steps_to_store
    logical(lg), optional, intent(in) :: hierarchical_basis
    logical(lg), optional, intent(in) :: static_condensation
    integer(ip), optional, intent(in) :: num_continuity 

    integer(ip) :: num_elems, ielem

    ! Create local fem space
    ! Note: When this subroutine is called from par_fem_space_create.f90, we have
    ! provided num_ghosts just in order to allocate element lists in 
    ! fem_space_allocate_structures that will also include ghost elements. 
    ! call fem_space_create( p_trian%f_trian, dhand, femsp, &
    !      & problem, approximations, bcond, continuity, order, material, which_approx, num_approximations, &
    !      & time_steps_to_store = time_steps_to_store, &
    !      & hierarchical_basis = hierarchical_basis, static_condensation = static_condensation, &
    !      & num_continuity = num_continuity, num_ghosts = p_trian%num_ghosts )


    call fem_space_allocate_structures(  p_trian%f_trian, dhand, femsp, num_approximations=num_approximations,&
         time_steps_to_store = time_steps_to_store, hierarchical_basis = hierarchical_basis, &
         static_condensation = static_condensation, num_continuity = num_continuity, &
         num_ghosts = p_trian%num_ghosts ) 

    assert(size(approximations)==num_approximations)
    femsp%approximations = approximations

    call fem_space_fe_list_create ( femsp, problem, which_approx, continuity, order, material, bcond )

    ! Communicate problem, continuity, order, and material
    !write(*,*) '***** EXCHANGE GHOST INFO *****'
    call ghost_elements_exchange ( p_trian%p_env%p_context%icontxt, p_trian%f_el_import, femsp%lelem )

    ! Create ghost fem space (only partially, i.e., previous info)
    ! write(*,*) '***** FILL GHOST ELEMENTS  *****'
    call ghost_fe_list_create ( femsp, p_trian%num_ghosts ) 

    call integration_faces_list( femsp )

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

  end subroutine par_fem_space_create


  !*********************************************************************************
  ! This subroutine takes the local finite element space, extended with ghost 
  ! elements that include (after the element exchange) values for order, continuity,
  ! material, problem. With this info, we create the fixed info pointers in ghost
  ! elements, which are needed in order to sort dofs on different processors the 
  ! same way.
  !*********************************************************************************
  subroutine ghost_fe_list_create( femsp, num_ghosts )
    implicit none
    type(fem_space), intent(inout), target :: femsp
    integer(ip), intent(in)        :: num_ghosts

    integer(ip) :: ielem, nvars, f_type, ivar, f_order, istat, pos_elinf, v_key
    logical(lg) :: created
    integer(ip) :: aux_val

    do ielem = femsp%g_trian%num_elems+1, femsp%g_trian%num_elems+num_ghosts
       !write (*,*) '************* GHOST ELEMENT *************',ielem
       nvars = femsp%dof_handler%problems(femsp%lelem(ielem)%problem)%p%nvars
       femsp%lelem(ielem)%num_vars = nvars
       f_type = femsp%g_trian%elems(ielem)%topology%ftype
       !write(*,*) 'f_type ghosts',f_type
       !write(*,*) 'nvars',nvars
       assert ( f_type > 0)
       call memalloc(nvars, femsp%lelem(ielem)%f_inf, __FILE__, __LINE__ )
       do ivar=1,nvars
          f_order = femsp%lelem(ielem)%order(ivar)
          !write(*,*) 'f_order',f_order
          v_key = femsp%g_trian%num_dims + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)*f_order

          call femsp%pos_elem_info%get(key=v_key,val=pos_elinf,stat=istat)
          if ( istat == new_index) then 
             !write (*,*) ' FIXED INFO NEW'
             call fem_element_fixed_info_create(femsp%lelem_info(pos_elinf),f_type,              &
                  &                             f_order,femsp%g_trian%num_dims,created)
             assert(created)
          end if
          femsp%lelem(ielem)%f_inf(ivar)%p => femsp%lelem_info(pos_elinf)

          ! aux_val = femsp%cur_elinf
          ! call femsp%ht_elem_info%put(key=v_key,val=aux_val,stat=istat)
          ! if ( istat == now_stored) then 
          !    write (*,*) ' FIXED INFO NEW'
          !    call fem_element_fixed_info_create(femsp%lelem_info(femsp%cur_elinf),f_type,              &
          !         &                             f_order,femsp%g_trian%num_dims,created)
          !    assert(created)
          !    pos_elinf = femsp%cur_elinf
          !    femsp%cur_elinf = femsp%cur_elinf + 1
          ! else if ( istat == was_stored ) then
          !    write (*,*) ' FIXED INFO ALREADY STORED '
          !    call femsp%ht_elem_info%get(key=v_key,val=pos_elinf,stat=istat)
          !    assert ( istat == key_found )
          ! end if
          ! femsp%lelem(ielem)%f_inf(ivar)%p => femsp%lelem_info(pos_elinf)
       end do
    end do
  end subroutine ghost_fe_list_create

  subroutine interface_faces_list( p_trian, femsp ) 
    implicit none
    ! Parameters
    type(par_triangulation), intent(in)       :: p_trian 
    type(fem_space), intent(inout)            :: femsp

    integer(ip) :: count_int, iobje, ielem, istat

    ! interface faces (among subdomains)

    count_int = 0
    do iobje = 1, p_trian%f_trian%num_objects
       if ( p_trian%f_trian%objects(iobje)%dimension == p_trian%f_trian%num_dims-1 ) then
          if (p_trian%objects(iobje)%interface == 1 ) then
             assert( p_trian%f_trian%objects(iobje)%num_elems_around == 2 )
             ielem = p_trian%f_trian%objects(iobje)%elems_around(1)
             count_int = count_int + 1
             !femsp%interface_faces(count_int) = iobje
          end if
       end if
    end do

    allocate( femsp%interface_faces(count_int), stat=istat)
    check ( istat == 0 )

    count_int = 0
    do iobje = 1, p_trian%f_trian%num_objects
       if ( p_trian%f_trian%objects(iobje)%dimension == p_trian%f_trian%num_dims-1 ) then
          if (p_trian%objects(iobje)%interface == 1 ) then
             assert( p_trian%f_trian%objects(iobje)%num_elems_around == 2 )
             ielem = p_trian%f_trian%objects(iobje)%elems_around(1)
             count_int = count_int + 1
             femsp%interface_faces(count_int)%face_object = iobje
          end if
       end if
    end do

    femsp%num_interface_faces = count_int

  end subroutine interface_faces_list

end module par_fem_space_names

