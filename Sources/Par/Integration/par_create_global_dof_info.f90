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
  use types
  use memor
  use fem_space_names
  use dof_handler_names
  use fem_triangulation_names
  use fem_graph_names
  use create_global_dof_info_names

  ! Par Modules
  use par_triangulation_names

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
  ! using the same serial functions. Finally, he gluing info for dofs is generated,
  ! based on vefs objects, the parallel triangulation info, and the ghost element info,
  ! relying on the fact that, with this info, we know how to put together interior and
  ! ghost elements, and so, number the dofs in a conforming way.
  !*********************************************************************************
  subroutine par_create_distributed_dof_info ( dhand, p_trian, femsp, dof_graph ) 
    implicit none
    type(par_triangulation), intent(inout) :: p_trian
    type(dof_handler), intent(in)          :: dhand
    type(fem_space), intent(inout)         :: femsp  
    type(fem_graph), allocatable, intent(out) :: dof_graph(:,:)

    call create_element_to_dof_and_ndofs( dhand, p_trian%f_trian, femsp )

    call ghost_discontinuous_Galerkin_dofs( dhand, p_trian, femsp  )

    call create_object2dof( dhand, p_trian%f_trian, femsp )

    call create_dof_graph( dhand, p_trian%f_trian, femsp, dof_graph )

    !call graph_distribution_create

  end subroutine par_create_distributed_dof_info

  !*********************************************************************************
  ! This subroutine takes the triangulation, the dof handler, and the finite element 
  ! space, and puts the dofs that are local but due to ghost elements when using 
  ! discontinuous Galerkin methods. In this case, the face dofs on the ghost element
  ! side are considered local elements (replicated among processors). 
  !*********************************************************************************
  subroutine ghost_discontinuous_Galerkin_dofs( dhand, p_trian, femsp )
    implicit none
    type(par_triangulation), intent(in)       :: p_trian
    type(dof_handler), intent(in)             :: dhand
    type(fem_space), intent(inout)            :: femsp

    integer(ip) :: iobje, i, ielem, l_faci, iprob, nvapb, ivars, l_var, g_var, inode, l_node, count
    integer(ip) :: iblock

    do iblock = 1, dhand%nblocks

       count = femsp%ndofs(iblock)
       do iobje = 1,1! interior faces!!!!trian%num_objects
          if ( p_trian%objects(iobje)%interface == -1 ) then
             do i=1,2
                ielem = p_trian%f_trian%objects(iobje)%elems_around(i)
                if ( ielem > p_trian%f_trian%num_elems ) exit
             end do
             assert ( i < 3 )
             l_faci = local_position(iobje,p_trian%f_trian%elems(ielem)%objects, &
                  & p_trian%f_trian%elems(ielem)%num_objects )
             iprob = femsp%lelem(ielem)%problem
             nvapb = dhand%prob_block(iblock,iprob)%nd1 
             do ivars = 1, nvapb
                l_var = dhand%prob_block(iblock,iprob)%a(ivars)
                g_var = dhand%problems(iprob)%p%l2g_var(l_var)
                do inode = femsp%lelem(ielem)%f_inf(l_var)%p%ntxob%p(l_faci), &
                     & femsp%lelem(ielem)%f_inf(l_var)%p%ntxob%p(l_faci+1)-1
                   l_node = femsp%lelem(ielem)%f_inf(l_var)%p%ntxob%l(inode)
                   if ( femsp%lelem(ielem)%elem2dof(l_node,l_var) == 0 ) then
                      count = count + 1
                      femsp%lelem(ielem)%elem2dof(l_node,l_var) = count
                   end if
                end do
             end do
          end if
       end do

    end do

  end subroutine ghost_discontinuous_Galerkin_dofs


  integer(ip) function local_position(key,list,size)
    implicit none
    integer(ip) :: key, size, list(size)
    
    do local_position = 1,size
       if ( list(local_position) == key) exit
    end do
    assert ( 0 == 1 )
    
  end function local_position

end module par_create_global_dof_info_names
