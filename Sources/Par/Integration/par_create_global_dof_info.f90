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
module par_create_global_dof_info_names

  ! Fem Modules
  use types_names
  use memor_names
  use sort_names
  use dof_handler_names
  use fem_triangulation_names
  use fem_graph_names
  use create_global_dof_info_names
  use fe_space_names
  use hash_table_names

  use finite_element_names

  ! Par Modules
  use par_triangulation_names
  use par_fe_space_names
  use par_graph_names
  use par_block_graph_names
  use block_dof_distribution_names
  use block_dof_distribution_create_names 

  implicit none
# include "debug.i90"
  private

  public :: par_create_distributed_dof_info

contains

  !*********************************************************************************
  ! This subroutine is intended to generate the dof generation and the dof graph 
  ! distribution. Here, the element2dof array is created as in serial, but next the
  ! dofs in ghost elements on interface faces coupled via interface discontinuous
  ! Galerkin terms are also included. Next, object2dof and the dof_graph are generated
  ! using the same serial functions. Finally, the gluing info for dofs is generated,
  ! based on VEFs objects, the parallel triangulation info, and the ghost element info,
  ! relying on the fact that, with this info, we know how to put together interior and
  ! ghost elements, and so, number the dofs in a conforming way.
  !*********************************************************************************
  subroutine par_create_distributed_dof_info ( dhand, p_trian, p_femsp, &
       blk_dof_dist, p_blk_graph, gtype ) 
    implicit none
    ! Paramters
    type(dof_handler_t)                  , intent(in)     :: dhand
    type(par_triangulation_t)            , intent(in)     :: p_trian
    type(par_fe_space_t)                , intent(inout)  :: p_femsp
    type(block_dof_distribution_t)       , intent(inout)  :: blk_dof_dist 
    type(par_block_graph_t)              , intent(inout)  :: p_blk_graph
    integer(ip)              , optional, intent(in)     :: gtype(dhand%nblocks) 

    integer(ip)               :: iblock, jblock
    type (par_graph_t), pointer :: p_graph

    ! Parallel environment MUST BE already created
    assert ( associated(p_femsp%p_trian) )
    assert ( p_femsp%p_trian%p_env%created )

    if( p_femsp%p_trian%p_env%p_context%iam >= 0 ) then
       call create_element_to_dof_and_ndofs( dhand, p_trian%f_trian, p_femsp%f_space )
       call ghost_dofs_by_integration( dhand, p_trian, p_femsp%f_space )
       call create_object2dof( dhand, p_trian%f_trian, p_femsp%f_space )
    end if

    ! Allocate block_dof_distribution
    call blk_dof_dist%alloc(p_trian%p_env, dhand%nblocks )

    ! Fill block_dof_distribution
    call block_dof_distribution_create ( p_trian, p_femsp, blk_dof_dist )

    ! Allocate par_block_graph
    call p_blk_graph%alloc(blk_dof_dist)

    ! Fill par_block_graph
    if( p_femsp%p_trian%p_env%p_context%iam >= 0 ) then
       do iblock = 1, dhand%nblocks
          do jblock = 1, dhand%nblocks
             p_graph => p_blk_graph%get_block(iblock,jblock)
             if ( iblock == jblock .and. present(gtype) ) then
                call create_dof_graph_block ( iblock, jblock, dhand, p_trian%f_trian, & 
                     p_femsp%f_space, p_graph%f_graph, gtype(iblock) )
             else
                call create_dof_graph_block ( iblock, jblock, dhand, p_trian%f_trian, &
                     p_femsp%f_space, p_graph%f_graph )
             end if
          end do
       end do

       ! call fem_graph_print ( 6, p_graph%f_graph )
    end if
  end subroutine par_create_distributed_dof_info

  !*********************************************************************************
  ! This subroutine takes the triangulation, the dof handler, and the finite element 
  ! space, and puts the dofs that are local but due to ghost elements when using 
  ! discontinuous Galerkin methods. In this case, the face dofs on the ghost element
  ! side are considered local elements (replicated among processors). 
  !*********************************************************************************
  subroutine ghost_dofs_by_integration( dhand, p_trian, femsp )
    implicit none
    type(dof_handler_t)      , intent(in)       :: dhand
    type(par_triangulation_t), intent(in)       :: p_trian
    type(fe_space_t)        , intent(inout)    :: femsp

    integer(ip) :: iobje, i, ielem, l_faci, iprob, nvapb, ivars, l_var, g_var, inode, l_node, count
    integer(ip) :: iblock, lobje

    do iblock = 1, dhand%nblocks

       count = femsp%ndofs(iblock)    

       do iobje = 1,femsp%num_interior_faces
          lobje = femsp%interior_faces(iobje)%face_object
          !do iobje = 1, p_trian%num_itfc_objs
          !lobje = p_trian%lst_itfc_objs(iobje)
          if ( p_trian%objects(lobje)%interface /= -1 ) then
             if ( p_trian%f_trian%objects(lobje)%dimension == p_trian%f_trian%num_dims -1 ) then 
                assert ( p_trian%f_trian%objects(lobje)%num_elems_around == 2 )
                do i=1,2
                   ielem = p_trian%f_trian%objects(lobje)%elems_around(i)
                   if ( ielem > p_trian%f_trian%num_elems ) exit
                end do
                assert ( i < 3 )            
                l_faci = local_position(lobje,p_trian%f_trian%elems(ielem)%objects, &
                     & p_trian%f_trian%elems(ielem)%num_objects )
                iprob = femsp%lelem(ielem)%problem
                nvapb = dhand%prob_block(iblock,iprob)%nd1 
                do ivars = 1, nvapb
                   l_var = dhand%prob_block(iblock,iprob)%a(ivars)
                   g_var = dhand%problems(iprob)%p%l2g_var(l_var)
                   if ( femsp%lelem(ielem)%continuity(g_var) == 0 ) then
                      do inode = femsp%lelem(ielem)%f_inf(l_var)%p%ntxob%p(l_faci), &
                           & femsp%lelem(ielem)%f_inf(l_var)%p%ntxob%p(l_faci+1)-1
                         l_node = femsp%lelem(ielem)%f_inf(l_var)%p%ntxob%l(inode)
                         assert ( femsp%lelem(ielem)%elem2dof(l_node,l_var) == 0  )
                         count = count + 1
                         femsp%lelem(ielem)%elem2dof(l_node,l_var) = count
                      end do
                   end if
                end do
             end if
          end if
       end do

       femsp%ndofs(iblock) = count
    end do

  end subroutine ghost_dofs_by_integration

  !*********************************************************************************
  ! Auxiliary function that says if key is in list
  !*********************************************************************************
  logical function is_object_in_element(key,list,size)
    implicit none
    integer(ip) :: key, size, list(size), i

    is_object_in_element = .false.
    do i = 1,size
       if ( list(i) == key)  then
          is_object_in_element = .true.
          exit
       end if
    end do

  end function is_object_in_element

  !*********************************************************************************
  ! Auxiliary function that says if key is in list
  !*********************************************************************************
  logical function is_object_touched(key,list,size)
    implicit none
    integer(ip) :: key, size, list(size), i

    is_object_touched = .false.
    do i = 1,size
       if ( list(i) == key)  then
          is_object_touched = .true.
          exit
       end if
    end do

  end function is_object_touched


  !*********************************************************************************
  ! Auxiliary function that returns the position of key in list
  !*********************************************************************************
  integer(ip) function local_position(key,list,size)
    implicit none
    integer(ip) :: key, size, list(size)

    do local_position = 1,size
       if ( list(local_position) == key) exit
    end do
    assert ( local_position < size + 1 )

  end function local_position


end module par_create_global_dof_info_names
