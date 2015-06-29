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
module par_fe_space_names
  ! Fem Modules
  use types_names
  use memor_names
  use fe_space_names
  use fe_space_types_names
  use problem_names
  use dof_descriptor_names
  use hash_table_names
  use finite_element_names

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

  type par_fe_space_t
     ! Data structure which stores the local part
     ! of the BC's mapped to the current subdomain
     type(fe_space_t)              :: fe_space
     
     integer(ip)                  :: num_interface_faces
     type(fe_face_t), allocatable  :: interface_faces(:)

     ! Pointer to parallel triangulation
     type(par_triangulation_t), pointer :: p_trian => NULL()
  end type par_fe_space_t

  ! Types
  public :: par_fe_space_t

  ! Methods
  public :: par_fe_space_create,  par_fe_space_print, par_fe_space_free

contains


  !*********************************************************************************
  ! This subroutine is intended to be called from a parallel driver, having as input
  ! a set of arrays of size number of local elements times number of global variables
  ! for continuity and order, and number of local elements for material and problem,
  ! together with some optional flags. The output of this subroutine is a fe_space
  ! with the required info on ghost elements.
  !*********************************************************************************
  subroutine par_fe_space_create ( p_trian, dof_descriptor, p_fe_space, problem, p_cond, continuity, order,     &
                                    material, which_approx, time_steps_to_store, hierarchical_basis, &
                                    static_condensation, num_continuity )
    implicit none
    ! Dummy arguments
    type(par_triangulation_t), target, intent(in)    :: p_trian
    type(dof_descriptor_t)              , intent(in)    :: dof_descriptor
    type(par_fe_space_t)            , intent(inout) :: p_fe_space  
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
    p_fe_space%p_trian => p_trian 

    if( p_fe_space%p_trian%p_env%p_context%iam >= 0 ) then

       call fe_space_allocate_structures(  p_trian%f_trian, dof_descriptor, p_fe_space%fe_space,            &
            time_steps_to_store = time_steps_to_store, hierarchical_basis = hierarchical_basis, &
            static_condensation = static_condensation, num_continuity = num_continuity, &
            num_ghosts = p_trian%num_ghosts ) 

!!$    AFM This allocate statement fails at runtime, i.e., istat/=0. Why ????
!!$    allocate ( fe_space%approximations(num_approximations), stat=istat)
!!$    write (*,*) 'XXX', num_approximations, istat
!!$    check (istat == 0)    

       call fe_space_fe_list_create ( p_fe_space%fe_space, problem, which_approx, continuity, order, material, p_cond%f_conditions )

       ! Communicate problem, continuity, order, and material
       !write(*,*) '***** EXCHANGE GHOST INFO *****'
       call ghost_elements_exchange ( p_trian%p_env%p_context%icontxt, p_trian%f_el_import, p_fe_space%fe_space%finite_elements )

       ! Create ghost fem space (only partially, i.e., previous info)
       ! write(*,*) '***** FILL GHOST ELEMENTS  *****'
       call ghost_fe_list_create ( p_fe_space ) 

       call integration_faces_list( p_fe_space%fe_space )

       !write(*,*) 'num_elems+1', p_fe_space%g_trian%num_elems+1
       !write(*,*) 'num_ghosts', num_ghosts
       !do ielem = p_fe_space%g_trian%num_elems+1, p_fe_space%g_trian%num_elems+num_ghosts
       !   call finite_element_fixed_info_write( p_trian%f_trian%elems(ielem)%geo_reference_element )
       !end do

       !write(*,*) 'num_elems+1', p_fe_space%g_trian%num_elems+1
       !write(*,*) 'num_ghosts', num_ghosts
       !do ielem = p_fe_space%g_trian%num_elems+1, p_fe_space%g_trian%num_elems+num_ghosts
       !   call finite_element_fixed_info_write( p_fe_space%g_trian%elems(ielem)%geo_reference_element )
       !end do
    end if

  end subroutine par_fe_space_create
  
  subroutine par_fe_space_free ( p_fe_space )
    implicit none
    ! Dummy arguments
    type(par_fe_space_t), intent(inout) :: p_fe_space  
    
    ! Local variables
    integer(ip) :: ielem
    
    ! Parallel environment MUST BE already created
    assert ( associated(p_fe_space%p_trian) )
    assert ( p_fe_space%p_trian%p_env%created )
    
    if( p_fe_space%p_trian%p_env%p_context%iam >= 0 ) then
       ! Deallocate type(finite_element_ts) associated to ghost elements
       do ielem = p_fe_space%p_trian%f_trian%num_elems+1, p_fe_space%p_trian%f_trian%num_elems+p_fe_space%p_trian%num_ghosts
          if(allocated(p_fe_space%fe_space%finite_elements(ielem)%reference_element_vars)) call memfree(p_fe_space%fe_space%finite_elements(ielem)%reference_element_vars,__FILE__,__LINE__)
          if(allocated(p_fe_space%fe_space%finite_elements(ielem)%elem2dof)) call memfree(p_fe_space%fe_space%finite_elements(ielem)%elem2dof,__FILE__,__LINE__)
          if(allocated(p_fe_space%fe_space%finite_elements(ielem)%unkno)) call memfree(p_fe_space%fe_space%finite_elements(ielem)%unkno,__FILE__,__LINE__)
          call finite_element_free_unpacked(p_fe_space%fe_space%finite_elements(ielem))
       end do
       call fe_space_free(p_fe_space%fe_space)
    end if
    
    nullify ( p_fe_space%p_trian )
  end subroutine par_fe_space_free

  subroutine par_fe_space_print ( p_fe_space )
    implicit none
    type(par_fe_space_t), intent(in) :: p_fe_space  

    
  end subroutine par_fe_space_print


  !*********************************************************************************
  ! This subroutine takes the local finite element space, extended with ghost 
  ! elements that include (after the element exchange) values for order, continuity,
  ! material, problem. With this info, we create the fixed info pointers in ghost
  ! elements, which are needed in order to sort dofs on different processors the 
  ! same way.
  !*********************************************************************************
  subroutine ghost_fe_list_create( p_fe_space )
    implicit none
    type(par_fe_space_t), target, intent(inout) :: p_fe_space

    integer(ip) :: ielem, nvars, f_type, ivar, f_order, istat, pos_elinf, v_key
    logical :: created
    integer(ip) :: aux_val, max_num_nodes, lndof, nnode

    do ielem = p_fe_space%p_trian%f_trian%num_elems+1, p_fe_space%p_trian%f_trian%num_elems+p_fe_space%p_trian%num_ghosts
       !write (*,*) '************* GHOST ELEMENT *************',ielem
       nvars = p_fe_space%fe_space%dof_descriptor%problems(p_fe_space%fe_space%finite_elements(ielem)%problem)%p%nvars
       p_fe_space%fe_space%finite_elements(ielem)%num_vars = nvars
       f_type = p_fe_space%p_trian%f_trian%elems(ielem)%geo_reference_element%ftype
       !write(*,*) 'f_type ghosts',f_type
       !write(*,*) 'nvars',nvars
       assert ( f_type > 0)
       call memalloc(nvars, p_fe_space%fe_space%finite_elements(ielem)%reference_element_vars, __FILE__, __LINE__ )
       allocate(p_fe_space%fe_space%finite_elements(ielem)%nodes_per_vef(nvars), stat=istat )
       do ivar=1,nvars
          f_order = p_fe_space%fe_space%finite_elements(ielem)%order(ivar)
          !write(*,*) 'f_order',f_order
          v_key = p_fe_space%p_trian%f_trian%num_dims + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)*f_order
          call p_fe_space%fe_space%pos_elem_info%get(key=v_key,val=pos_elinf,stat=istat)
          if ( istat == new_index) then 
             !write (*,*) ' FIXED INFO NEW'
             call finite_element_fixed_info_create(p_fe_space%fe_space%finite_elements_info(pos_elinf),f_type,              &
                  &                             f_order,p_fe_space%p_trian%f_trian%num_dims)
          end if
          p_fe_space%fe_space%finite_elements(ielem)%reference_element_vars(ivar)%p => p_fe_space%fe_space%finite_elements_info(pos_elinf)
          if ( p_fe_space%fe_space%finite_elements(ielem)%continuity(ivar) /= 0 ) then
             p_fe_space%fe_space%finite_elements(ielem)%nodes_per_vef(ivar)%p => p_fe_space%fe_space%finite_elements_info(pos_elinf)%ndxob
          else 
             p_fe_space%fe_space%finite_elements(ielem)%nodes_per_vef(ivar)%p => p_fe_space%fe_space%finite_elements_info(pos_elinf)%ndxob_int
          end if
       end do

       ! Compute number of DOFs in the elemental matrix associated to ielem
       lndof = 0
       max_num_nodes = 0
       do ivar=1,nvars
          if ( p_fe_space%fe_space%static_condensation ) then
             nnode = p_fe_space%fe_space%finite_elements(ielem)%reference_element_vars(ivar)%p%nnode -                                   &
                  &  p_fe_space%fe_space%finite_elements(ielem)%reference_element_vars(ivar)%p%nodes_vef(p_fe_space%fe_space%g_trian%num_dims+1) ! SB.alert : do not use nodes_vef
          else
             nnode = p_fe_space%fe_space%finite_elements(ielem)%reference_element_vars(ivar)%p%nnode 
          end if
          lndof = lndof + nnode
          max_num_nodes = max(max_num_nodes,nnode)
       end do

       call memalloc( max_num_nodes, nvars, p_fe_space%fe_space%finite_elements(ielem)%elem2dof, __FILE__,__LINE__ )
       p_fe_space%fe_space%finite_elements(ielem)%elem2dof = 0
       call memalloc( max_num_nodes, nvars, p_fe_space%fe_space%time_steps_to_store, p_fe_space%fe_space%finite_elements(ielem)%unkno, __FILE__,__LINE__)
       p_fe_space%fe_space%finite_elements(ielem)%unkno = 0.0_rp
    end do

  end subroutine ghost_fe_list_create

  subroutine interface_faces_list( p_trian, p_fe_space ) 
    implicit none
    ! Parameters
    type(par_triangulation_t), intent(in)     :: p_trian 
    type(par_fe_space_t)    , intent(inout)  :: p_fe_space

    integer(ip) :: count_int, iobje, ielem, istat

    ! interface faces (among subdomains)

    count_int = 0
    do iobje = 1, p_trian%f_trian%num_vefs
       if ( p_trian%f_trian%vefs(iobje)%dimension == p_trian%f_trian%num_dims-1 ) then
          if (p_trian%vefs(iobje)%interface == 1 ) then
             assert( p_trian%f_trian%vefs(iobje)%num_elems_around == 2 )
             ielem = p_trian%f_trian%vefs(iobje)%elems_around(1)
             count_int = count_int + 1
          end if
       end if
    end do

    allocate( p_fe_space%interface_faces(count_int), stat=istat)
    check ( istat == 0 )

    count_int = 0
    do iobje = 1, p_trian%f_trian%num_vefs
       if ( p_trian%f_trian%vefs(iobje)%dimension == p_trian%f_trian%num_dims-1 ) then
          if (p_trian%vefs(iobje)%interface == 1 ) then
             assert( p_trian%f_trian%vefs(iobje)%num_elems_around == 2 )
             ielem = p_trian%f_trian%vefs(iobje)%elems_around(1)
             count_int = count_int + 1
             p_fe_space%interface_faces(count_int)%face_vef = iobje
          end if
       end if
    end do

   p_fe_space%num_interface_faces = count_int

  end subroutine interface_faces_list

end module par_fe_space_names

