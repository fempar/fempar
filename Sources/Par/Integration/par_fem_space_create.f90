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
module par_fem_space_names
  ! Fem Modules
  use types_names
  use memor_names
  use fem_space_names
  use fem_space_types_names
  use problem_names
  use dof_handler_names
  use hash_table_names
  use fem_element_names

  ! Par Modules
  use par_environment_names
  use par_triangulation_names
  use par_element_exchange_names
  use par_conditions_names

#ifdef memcheck
use iso_c_binding
#endif
  implicit none
# include "debug.i90"
  private

  type par_fem_space_t
     ! Data structure which stores the local part
     ! of the BC's mapped to the current subdomain
     type(fem_space_t)              :: f_space
     
     integer(ip)                  :: num_interface_faces
     type(fem_face_t), allocatable  :: interface_faces(:)

     ! Pointer to parallel triangulation
     type(par_triangulation_t), pointer :: p_trian => NULL()
  end type par_fem_space_t

  ! Types
  public :: par_fem_space_t

  ! Methods
  public :: par_fem_space_create,  par_fem_space_print, par_fem_space_free

contains


  !*********************************************************************************
  ! This subroutine is intended to be called from a parallel driver, having as input
  ! a set of arrays of size number of local elements times number of global variables
  ! for continuity and order, and number of local elements for material and problem,
  ! together with some optional flags. The output of this subroutine is a fem_space
  ! with the required info on ghost elements.
  !*********************************************************************************
  subroutine par_fem_space_create ( p_trian, dhand, p_femsp, problem, p_cond, continuity, order,     &
                                    material, which_approx, time_steps_to_store, hierarchical_basis, &
                                    static_condensation, num_continuity )
    implicit none
    ! Dummy arguments
    type(par_triangulation_t), target, intent(in)    :: p_trian
    type(dof_handler_t)              , intent(in)    :: dhand
    type(par_fem_space_t)            , intent(inout) :: p_femsp  
    integer(ip)                    , intent(in)    :: problem(:)
    type(par_conditions_t)           , intent(in)    :: p_cond
    integer(ip)                    , intent(in)    :: continuity(:,:)
    integer(ip)                    , intent(in)    :: order(:,:)
    integer(ip)                    , intent(in)    :: material(:)
    integer(ip)                    , intent(in)    :: which_approx(:)
    integer(ip)          , optional, intent(in)    :: time_steps_to_store
    logical          , optional, intent(in)    :: hierarchical_basis
    logical          , optional, intent(in)    :: static_condensation
    integer(ip)          , optional, intent(in)    :: num_continuity   
    
    ! Local variables
    integer(ip) :: istat

    ! Parallel environment MUST BE already created
    assert ( p_trian%p_env%created )
    p_femsp%p_trian => p_trian 

    if( p_femsp%p_trian%p_env%p_context%iam >= 0 ) then

       call fem_space_allocate_structures(  p_trian%f_trian, dhand, p_femsp%f_space,            &
            time_steps_to_store = time_steps_to_store, hierarchical_basis = hierarchical_basis, &
            static_condensation = static_condensation, num_continuity = num_continuity, &
            num_ghosts = p_trian%num_ghosts ) 

!!$    AFM This allocate statement fails at runtime, i.e., istat/=0. Why ????
!!$    allocate ( fspac%approximations(num_approximations), stat=istat)
!!$    write (*,*) 'XXX', num_approximations, istat
!!$    check (istat == 0)    

       call fem_space_fe_list_create ( p_femsp%f_space, problem, which_approx, continuity, order, material, p_cond%f_conditions )

       ! Communicate problem, continuity, order, and material
       !write(*,*) '***** EXCHANGE GHOST INFO *****'
       call ghost_elements_exchange ( p_trian%p_env%p_context%icontxt, p_trian%f_el_import, p_femsp%f_space%lelem )

       ! Create ghost fem space (only partially, i.e., previous info)
       ! write(*,*) '***** FILL GHOST ELEMENTS  *****'
       call ghost_fe_list_create ( p_femsp ) 

       call integration_faces_list( p_femsp%f_space )

       !write(*,*) 'num_elems+1', p_femsp%g_trian%num_elems+1
       !write(*,*) 'num_ghosts', num_ghosts
       !do ielem = p_femsp%g_trian%num_elems+1, p_femsp%g_trian%num_elems+num_ghosts
       !   call fem_element_fixed_info_write( p_trian%f_trian%elems(ielem)%topology )
       !end do

       !write(*,*) 'num_elems+1', p_femsp%g_trian%num_elems+1
       !write(*,*) 'num_ghosts', num_ghosts
       !do ielem = p_femsp%g_trian%num_elems+1, p_femsp%g_trian%num_elems+num_ghosts
       !   call fem_element_fixed_info_write( p_femsp%g_trian%elems(ielem)%topology )
       !end do
    end if

  end subroutine par_fem_space_create
  
  subroutine par_fem_space_free ( p_femsp )
    implicit none
    ! Dummy arguments
    type(par_fem_space_t), intent(inout) :: p_femsp  
    
    ! Local variables
    integer(ip) :: ielem
    
    ! Parallel environment MUST BE already created
    assert ( associated(p_femsp%p_trian) )
    assert ( p_femsp%p_trian%p_env%created )
    
    if( p_femsp%p_trian%p_env%p_context%iam >= 0 ) then
       ! Deallocate type(fem_element_ts) associated to ghost elements
       do ielem = p_femsp%p_trian%f_trian%num_elems+1, p_femsp%p_trian%f_trian%num_elems+p_femsp%p_trian%num_ghosts
          if(allocated(p_femsp%f_space%lelem(ielem)%f_inf)) call memfree(p_femsp%f_space%lelem(ielem)%f_inf,__FILE__,__LINE__)
          if(allocated(p_femsp%f_space%lelem(ielem)%elem2dof)) call memfree(p_femsp%f_space%lelem(ielem)%elem2dof,__FILE__,__LINE__)
          if(allocated(p_femsp%f_space%lelem(ielem)%unkno)) call memfree(p_femsp%f_space%lelem(ielem)%unkno,__FILE__,__LINE__)
          call fem_element_free_unpacked(p_femsp%f_space%lelem(ielem))
       end do
       call fem_space_free(p_femsp%f_space)
    end if
    
    nullify ( p_femsp%p_trian )
  end subroutine par_fem_space_free

  subroutine par_fem_space_print ( p_femsp )
    implicit none
    type(par_fem_space_t), intent(in) :: p_femsp  

    
  end subroutine par_fem_space_print


  !*********************************************************************************
  ! This subroutine takes the local finite element space, extended with ghost 
  ! elements that include (after the element exchange) values for order, continuity,
  ! material, problem. With this info, we create the fixed info pointers in ghost
  ! elements, which are needed in order to sort dofs on different processors the 
  ! same way.
  !*********************************************************************************
  subroutine ghost_fe_list_create( p_femsp )
    implicit none
    type(par_fem_space_t), target, intent(inout) :: p_femsp

    integer(ip) :: ielem, nvars, f_type, ivar, f_order, istat, pos_elinf, v_key
    logical :: created
    integer(ip) :: aux_val, max_num_nodes, lndof, nnode

    do ielem = p_femsp%p_trian%f_trian%num_elems+1, p_femsp%p_trian%f_trian%num_elems+p_femsp%p_trian%num_ghosts
       !write (*,*) '************* GHOST ELEMENT *************',ielem
       nvars = p_femsp%f_space%dof_handler%problems(p_femsp%f_space%lelem(ielem)%problem)%p%nvars
       p_femsp%f_space%lelem(ielem)%num_vars = nvars
       f_type = p_femsp%p_trian%f_trian%elems(ielem)%topology%ftype
       !write(*,*) 'f_type ghosts',f_type
       !write(*,*) 'nvars',nvars
       assert ( f_type > 0)
       call memalloc(nvars, p_femsp%f_space%lelem(ielem)%f_inf, __FILE__, __LINE__ )
       allocate(p_femsp%f_space%lelem(ielem)%nodes_object(nvars), stat=istat )
       do ivar=1,nvars
          f_order = p_femsp%f_space%lelem(ielem)%order(ivar)
          !write(*,*) 'f_order',f_order
          v_key = p_femsp%p_trian%f_trian%num_dims + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)*f_order
          call p_femsp%f_space%pos_elem_info%get(key=v_key,val=pos_elinf,stat=istat)
          if ( istat == new_index) then 
             !write (*,*) ' FIXED INFO NEW'
             call fem_element_fixed_info_create(p_femsp%f_space%lelem_info(pos_elinf),f_type,              &
                  &                             f_order,p_femsp%p_trian%f_trian%num_dims,created)
             assert(created)
          end if
          p_femsp%f_space%lelem(ielem)%f_inf(ivar)%p => p_femsp%f_space%lelem_info(pos_elinf)
          if ( p_femsp%f_space%lelem(ielem)%continuity(ivar) /= 0 ) then
             p_femsp%f_space%lelem(ielem)%nodes_object(ivar)%p => p_femsp%f_space%lelem_info(pos_elinf)%ndxob
          else 
             p_femsp%f_space%lelem(ielem)%nodes_object(ivar)%p => p_femsp%f_space%lelem_info(pos_elinf)%ndxob_int
          end if
       end do

       ! Compute number of DOFs in the elemental matrix associated to ielem
       lndof = 0
       max_num_nodes = 0
       do ivar=1,nvars
          if ( p_femsp%f_space%static_condensation ) then
             nnode = p_femsp%f_space%lelem(ielem)%f_inf(ivar)%p%nnode -                                   &
                  &  p_femsp%f_space%lelem(ielem)%f_inf(ivar)%p%nodes_obj(p_femsp%f_space%g_trian%num_dims+1) ! SB.alert : do not use nodes_obj
          else
             nnode = p_femsp%f_space%lelem(ielem)%f_inf(ivar)%p%nnode 
          end if
          lndof = lndof + nnode
          max_num_nodes = max(max_num_nodes,nnode)
       end do

       call memalloc( max_num_nodes, nvars, p_femsp%f_space%lelem(ielem)%elem2dof, __FILE__,__LINE__ )
       p_femsp%f_space%lelem(ielem)%elem2dof = 0
       call memalloc( max_num_nodes, nvars, p_femsp%f_space%time_steps_to_store, p_femsp%f_space%lelem(ielem)%unkno, __FILE__,__LINE__)
       p_femsp%f_space%lelem(ielem)%unkno = 0.0_rp
    end do

  end subroutine ghost_fe_list_create

  subroutine interface_faces_list( p_trian, p_femsp ) 
    implicit none
    ! Parameters
    type(par_triangulation_t), intent(in)     :: p_trian 
    type(par_fem_space_t)    , intent(inout)  :: p_femsp

    integer(ip) :: count_int, iobje, ielem, istat

    ! interface faces (among subdomains)

    count_int = 0
    do iobje = 1, p_trian%f_trian%num_objects
       if ( p_trian%f_trian%objects(iobje)%dimension == p_trian%f_trian%num_dims-1 ) then
          if (p_trian%objects(iobje)%interface == 1 ) then
             assert( p_trian%f_trian%objects(iobje)%num_elems_around == 2 )
             ielem = p_trian%f_trian%objects(iobje)%elems_around(1)
             count_int = count_int + 1
          end if
       end if
    end do

    allocate( p_femsp%interface_faces(count_int), stat=istat)
    check ( istat == 0 )

    count_int = 0
    do iobje = 1, p_trian%f_trian%num_objects
       if ( p_trian%f_trian%objects(iobje)%dimension == p_trian%f_trian%num_dims-1 ) then
          if (p_trian%objects(iobje)%interface == 1 ) then
             assert( p_trian%f_trian%objects(iobje)%num_elems_around == 2 )
             ielem = p_trian%f_trian%objects(iobje)%elems_around(1)
             count_int = count_int + 1
             p_femsp%interface_faces(count_int)%face_object = iobje
          end if
       end if
    end do

   p_femsp%num_interface_faces = count_int

  end subroutine interface_faces_list

end module par_fem_space_names

