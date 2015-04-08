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
module assembly_names

  use types
  use fem_space_names
  use dof_handler_names
  use fem_block_matrix_names
  use fem_matrix_names
  use fem_graph_names

  implicit none
# include "debug.i90"
  private

  public :: assembly_element_matrix_global_system_matrix

  interface assembly_element_matrix_global_system_matrix
     module procedure assembly_element_matrix_global_system_block_matrix, &
          & assembly_element_matrix_global_system_monolithic_matrix
  end interface assembly_element_matrix_global_system_matrix
contains

  subroutine assembly_element_matrix_global_system_block_matrix(  elem, dhand, a ) 
    implicit none
    type(dof_handler), intent(in)             :: dhand
    type(fem_element), intent(in)             :: elem
    type(fem_block_matrix), intent(inout)     :: a

    integer(ip) :: ivar, start(dhand%problems(elem%problem)%nvars+1)
    ! JP: end not needed?
    integer(ip) :: end(dhand%problems(elem%problem)%nvars+1), iblock, jblock
    
    ! JP: I also need start when writing problems to acces the element matrix...
    ! Assuming a monolithic problem matrix
    do ivar = 1,dhand%problems(elem%problem)%nvars
      start(ivar+1) = start(ivar+1) + elem%f_inf(ivar)%p%nnode
    end do

    start(1) = 1
    do ivar = 2, dhand%problems(elem%problem)%nvars+1
       start(ivar) = start(ivar) + start(ivar-1)
    end do

    do iblock = 1, dhand%nblocks
       do jblock = 1, dhand%nblocks
          call element_matrix_block_assembly( dhand, elem, start, a%blocks(iblock,jblock)%p_f_matrix, iblock, jblock )
       end do
    end do

  end subroutine assembly_element_matrix_global_system_block_matrix

  subroutine assembly_element_matrix_global_system_monolithic_matrix(  elem, dhand, a ) 
    implicit none
    type(dof_handler), intent(in)             :: dhand
    type(fem_element), intent(in)             :: elem
    type(fem_matrix), intent(inout)           :: a

    integer(ip) :: ivar, start(dhand%problems(elem%problem)%nvars+1)
    ! JP: end not needed?
    integer(ip) :: end(dhand%problems(elem%problem)%nvars+1), iblock, jblock
    

    ! Assuming a monolithic problem matrix
    do ivar = 1,dhand%problems(elem%problem)%nvars
      start(ivar+1) = start(ivar+1) + elem%f_inf(ivar)%p%nnode
    end do

    start(1) = 1
    do ivar = 2, dhand%problems(elem%problem)%nvars+1
       start(ivar) = start(ivar) + start(ivar-1)
    end do

    call element_matrix_block_assembly( dhand, elem, start, a )

  end subroutine assembly_element_matrix_global_system_monolithic_matrix
  
  ! OTHER ALTERNATIVES TO BE IMPLEMENTED :

  ! subroutine assembly_monolithic_element_matrix_one_block_system_matrix
  ! end subroutine assembly_monolithic_element_matrix_one_block_system_matrix

  ! subroutine assembly_one_block_element_matrix_one_block_system_matrix
  ! end subroutine assembly_one_block_element_matrix_one_block_system_matrix

  ! subroutine assembly_one_block_element_matrix_global_system_matrix
  ! end subroutine assembly_one_block_element_matrix_global_system_matrix

  subroutine element_matrix_block_assembly( dhand, elem,start, a, iblock, jblock )
    implicit none
    ! Parameters
    type(dof_handler), intent(in)             :: dhand
    type(fem_element), intent(in)             :: elem
    integer(ip), intent(in)                   :: start(dhand%problems(elem%problem)%nvars+1)
    type(fem_matrix), intent(inout)           :: a
    integer(ip), intent(in), optional         :: iblock, jblock

    integer(ip) :: gtype, iprob, nvapb_i, nvapb_j, ivars, jvars, l_var, m_var, g_var, k_var
    integer(ip) :: inode, jnode, idof, jdof, k, nvapb, iblock_, jblock_
    
    iblock_ = 1
    jblock_ = 1
    
    if ( present(iblock) ) iblock_ = iblock
    if ( present(jblock) ) jblock_ = jblock
    
    gtype = a%gr%type

    iprob = elem%problem

    nvapb_i = dhand%prob_block(iblock,iprob)%nd1
    do ivars = 1, nvapb
       l_var = dhand%prob_block(iblock,iprob)%a(ivars)
       g_var = dhand%problems(iprob)%l2g_var(l_var)
       nvapb_j = dhand%prob_block(jblock,iprob)%nd1
       do jvars = 1, nvapb
          m_var = dhand%prob_block(iblock,iprob)%a(jvars)
          k_var = dhand%problems(iprob)%l2g_var(m_var)
          do inode = 1,elem%f_inf(l_var)%p%nnode
             idof = elem%elem2dof(inode,l_var)
             if ( idof  > 0 ) then
                do jnode = 1,elem%f_inf(m_var)%p%nnode
                   jdof = elem%elem2dof(jnode,m_var)
                   if (  gtype == csr .and. jdof > 0 ) then
                      do k = a%gr%ia(idof),a%gr%ia(idof+1)-1
                         if ( a%gr%ja(k) == jdof ) exit
                      end do
                      assert ( k < a%gr%ia(idof+1) )
                      a%a(k) = elem%p_mat%a(start(l_var)+inode,start(m_var)+jnode)
                   else if ( jdof >= idof ) then! gtype == csr_symm 
                      do k = a%gr%ia(idof),a%gr%ia(idof+1)-1
                         if ( a%gr%ja(k) == jdof ) exit
                      end do
                      assert ( k < a%gr%ia(idof+1) )
                      a%a(k) = elem%p_mat%a(start(l_var)+inode,start(m_var)+jnode)
                   end if
                end do
             end if
          end do
       end do
    end do

  end subroutine element_matrix_block_assembly

end module assembly_names



