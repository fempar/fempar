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

use types_names
  use array_ip1_names
  use fem_element_names
  !use fem_space_names
  use integrable_names
  use dof_handler_names
  use fem_block_matrix_names
  use fem_matrix_names
  use fem_block_vector_names
  use fem_vector_names
  use fem_graph_names
  use plain_vector_names

  implicit none
# include "debug.i90"
  private

  ! Functions
  public :: assembly

contains

  subroutine assembly(elem, dhand, a ) 
    implicit none
    ! Parameters
    type(dof_handler_t), intent(in)    :: dhand
    type(fem_element_t), intent(in)    :: elem
    class(integrable_t), intent(inout) :: a

    select type(a)
    class is(fem_matrix_t)
       call assembly_element_matrix_mono(elem, dhand, a) 
    class is(fem_vector_t)
       call assembly_element_vector_mono(elem, dhand, a)
    class is(plain_vector_t)
       call assembly_element_plain_vector(elem, dhand, a)
       !class is(fem_block_matrix_t)
       !    call assembly_element_matrix_block(elem, dhand, a)
       ! class is(fem_block_vector_t)
       !    call assembly_element_vector_block(elem, dhand, a)
    class default
       ! class not yet implemented
       check(.false.)
    end select
  end subroutine assembly

  subroutine assembly_element_matrix_block(  elem, dhand, a ) 
    implicit none
    type(dof_handler_t), intent(in)             :: dhand
    type(fem_element_t), intent(in)             :: elem
    type(fem_block_matrix_t), intent(inout)     :: a
    
    integer(ip) :: ivar, iblock, jblock

    type(fem_matrix_t), pointer :: f_matrix

    do iblock = 1, dhand%nblocks
       do jblock = 1, dhand%nblocks
          f_matrix => a%get_block(iblock,jblock)
          if ( associated(f_matrix) ) then
             call element_matrix_assembly( dhand, elem, f_matrix, iblock, jblock )
          end if 
       end do
    end do

  end subroutine assembly_element_matrix_block

  subroutine assembly_element_matrix_mono(  elem, dhand, a ) 
    implicit none
    type(dof_handler_t), intent(in)             :: dhand
    type(fem_element_t), intent(in)             :: elem
    type(fem_matrix_t), intent(inout)           :: a

    !integer(ip) :: start(dhand%problems(elem%problem)%p%nvars+1)
    call element_matrix_assembly( dhand, elem, a )

  end subroutine assembly_element_matrix_mono

  subroutine assembly_face_element_matrix_block(  face, elem, dhand, a ) 
    implicit none
    type(dof_handler_t), intent(in)             :: dhand
    type(fem_face_t)   , intent(in)             :: face
    type(fem_element_t), intent(in)             :: elem(2)
    type(fem_block_matrix_t), intent(inout)     :: a

    integer(ip) :: iblock, jblock, i
    type(array_ip1_t) :: start(2)
    type(fem_matrix_t), pointer :: f_matrix

!!$    do i=1,2
!!$       call pointer_variable(  elem(i), dhand, start(i)%a )
!!$    end do
!!$    start(2)%a = start(2)%a + start(1)%a(dhand%problems(elem(1)%problem)%p%nvars+1) - 1
    elem(2)%start%a = elem(2)%start%a + elem(1)%start%a(dhand%problems(elem(1)%problem)%p%nvars+1) - 1

    do iblock = 1, dhand%nblocks
       do jblock = 1, dhand%nblocks
          f_matrix => a%get_block(iblock,jblock)
          if ( associated(f_matrix) ) then
            do i = 1,2
               call element_matrix_assembly( dhand, elem(i), f_matrix, iblock, jblock )
            end do
            call face_element_matrix_assembly( dhand, elem, face, f_matrix, iblock, jblock )
          end if
       end do
    end do

  end subroutine assembly_face_element_matrix_block

  subroutine assembly_face_element_matrix_mono(  face, elem, dhand, a ) 
    implicit none
    type(dof_handler_t), intent(in)             :: dhand
    type(fem_face_t)   , intent(in)             :: face
    type(fem_element_t), intent(in)             :: elem(2)
    type(fem_matrix_t), intent(inout)     :: a

    integer(ip) :: i
    type(array_ip1_t) :: start(2)

!!$    do i=1,2
!!$       call pointer_variable(  elem(i), dhand, start(i)%a )
!!$    end do
!!$    start(2)%a = start(2)%a + start(1)%a(dhand%problems(elem(1)%problem)%p%nvars+1) - 1
    elem(2)%start%a = elem(2)%start%a + elem(1)%start%a(dhand%problems(elem(1)%problem)%p%nvars+1) - 1
    do i = 1,2
       call element_matrix_assembly( dhand, elem(i), a )
    end do
    call face_element_matrix_assembly( dhand, elem, face, a )

  end subroutine assembly_face_element_matrix_mono

  subroutine assembly_element_vector_block(  elem, dhand, a ) 
    implicit none
    type(dof_handler_t), intent(in)             :: dhand
    type(fem_element_t), intent(in)             :: elem
    type(fem_block_vector_t), intent(inout)     :: a

    integer(ip) :: iblock

    do iblock = 1, dhand%nblocks
       call element_vector_assembly( dhand, elem, a%blocks(iblock), &
            & iblock )
    end do

  end subroutine assembly_element_vector_block


  subroutine assembly_element_vector_mono(  elem, dhand, a ) 
    implicit none
    type(dof_handler_t), intent(in)             :: dhand
    type(fem_element_t), intent(in)             :: elem
    type(fem_vector_t), intent(inout)           :: a

    call element_vector_assembly( dhand, elem, a )

  end subroutine assembly_element_vector_mono

  subroutine assembly_element_plain_vector (  elem, dhand, a ) 
    implicit none
    type(dof_handler_t) , intent(in)    :: dhand
    type(fem_element_t) , intent(in)    :: elem
    type(plain_vector_t), intent(inout) :: a
    
    a%b = a%b + elem%p_plain_vector%a

  end subroutine assembly_element_plain_vector

  subroutine assembly_face_vector_block(  face, elem, dhand, a ) 
    implicit none
    type(dof_handler_t), intent(in)             :: dhand
    type(fem_face_t)   , intent(in)             :: face
    type(fem_element_t), intent(in)             :: elem
    type(fem_block_vector_t), intent(inout)     :: a

    integer(ip) :: iblock

    ! Note: This subroutine only has sense on interface / boundary faces, only
    ! related to ONE element.
    do iblock = 1, dhand%nblocks
       call face_vector_assembly( dhand, elem, face, a%blocks(iblock), iblock )
    end do

  end subroutine assembly_face_vector_block

  subroutine assembly_face_vector_mono(  face, elem, dhand, a ) 
    implicit none
    type(dof_handler_t), intent(in)             :: dhand
    type(fem_face_t)   , intent(in)             :: face
    type(fem_element_t), intent(in)             :: elem
    type(fem_vector_t), intent(inout)           :: a


    ! Note: This subroutine only has sense on interface / boundary faces, only
    ! related to ONE element.
    call face_vector_assembly( dhand, elem, face, a )

  end subroutine assembly_face_vector_mono

  subroutine element_matrix_assembly( dhand, elem, a, iblock, jblock )
    implicit none
    ! Parameters
    type(dof_handler_t), intent(in)             :: dhand
    type(fem_element_t), intent(in)             :: elem
    type(fem_matrix_t), intent(inout)           :: a
    integer(ip), intent(in), optional         :: iblock, jblock

    integer(ip) :: gtype, iprob, nvapb_i, nvapb_j, ivars, jvars, l_var, m_var, g_var, k_var
    integer(ip) :: inode, jnode, idof, jdof, k, iblock_, jblock_

    iblock_ = 1
    jblock_ = 1
    if ( present(iblock) ) iblock_ = iblock
    if ( present(jblock) ) jblock_ = jblock
    
    
    gtype = a%gr%type
    iprob = elem%problem


    nvapb_i = dhand%prob_block(iblock_,iprob)%nd1
    nvapb_j = dhand%prob_block(jblock_,iprob)%nd1

    !write (*,*) 'nvapb_i:',nvapb_i
    !write (*,*) 'nvapb_j:',nvapb_j
    !write (*,*) 'start:',elem%start%a

    !write(*,*) 'local matrix',elem%p_mat%nd1, elem%p_mat%nd2
    !write(*,*) 'local matrix',elem%p_mat%a
    
    

    do ivars = 1, nvapb_i
       l_var = dhand%prob_block(iblock_,iprob)%a(ivars)
       g_var = dhand%problems(iprob)%p%l2g_var(l_var)
       do jvars = 1, nvapb_j
          m_var = dhand%prob_block(jblock_,iprob)%a(jvars)
          k_var = dhand%problems(iprob)%p%l2g_var(m_var)
          !write (*,*) 'l_var:',l_var
          !write (*,*) 'm_var:',m_var
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
                      a%a(k) = a%a(k) + elem%p_mat%a(elem%start%a(l_var)+inode-1,elem%start%a(m_var)+jnode-1)
                   else if ( jdof >= idof ) then! gtype == csr_symm 
                      do k = a%gr%ia(idof),a%gr%ia(idof+1)-1
                         if ( a%gr%ja(k) == jdof ) exit
                      end do
                      assert ( k < a%gr%ia(idof+1) )
                      a%a(k) = a%a(k) + elem%p_mat%a(elem%start%a(l_var)+inode-1,elem%start%a(m_var)+jnode-1)
                   end if
                end do
             end if
          end do
       end do
    end do
    
    !write(*,*) 'system_matrix'
    !call fem_matrix_print(6,a)

  end subroutine element_matrix_assembly

  subroutine face_element_matrix_assembly( dhand, elem, face, a, iblock, jblock )
    implicit none
    ! Parameters
    type(dof_handler_t), intent(in)             :: dhand
    type(fem_element_t), intent(in)             :: elem(2)
    type(fem_face_t)   , intent(in)             :: face
    type(fem_matrix_t), intent(inout)           :: a
    integer(ip), intent(in), optional         :: iblock, jblock

    integer(ip) :: gtype, iprob, jprob, nvapb_i, nvapb_j, ivars, jvars, l_var, m_var, g_var, k_var
    integer(ip) :: inode, jnode, idof, jdof, k, iblock_, jblock_, iobje, i, j, ndime

    iblock_ = 1
    jblock_ = 1
    if ( present(iblock) ) iblock_ = iblock
    if ( present(jblock) ) jblock_ = jblock
    gtype = a%gr%type
    do i = 1, 2
       j = 3 - i
       iprob = elem(i)%problem
       jprob = elem(j)%problem
       nvapb_i = dhand%prob_block(iblock_,iprob)%nd1
       nvapb_j = dhand%prob_block(jblock_,jprob)%nd1
       do ivars = 1, nvapb_i
          l_var = dhand%prob_block(iblock_,iprob)%a(ivars)
          g_var = dhand%problems(iprob)%p%l2g_var(l_var)
          do jvars = 1, nvapb_j
             m_var = dhand%prob_block(jblock_,jprob)%a(jvars)
             k_var = dhand%problems(jprob)%p%l2g_var(m_var)
             do inode = 1,elem(i)%f_inf(l_var)%p%nnode
                idof = elem(i)%elem2dof(inode,l_var)
                if ( idof  > 0 ) then
                   ndime = elem(j)%p_geo_info%ndime
                   iobje = face%face_object + elem(j)%p_geo_info%nobje_dim(ndime) - 1
                   do jnode = elem(j)%f_inf(m_var)%p%ntxob%p(iobje),elem(j)%f_inf(m_var)%p%ntxob%p(iobje+1)-1
                      jdof = elem(j)%elem2dof(elem(j)%f_inf(m_var)%p%ntxob%l(jnode),m_var)
                      if (  gtype == csr .and. jdof > 0 ) then
                         do k = a%gr%ia(idof),a%gr%ia(idof+1)-1
                            if ( a%gr%ja(k) == jdof ) exit
                         end do
                         assert ( k < a%gr%ia(idof+1) )
                         a%a(k) = a%a(k) + face%p_mat%a(elem(i)%start%a(l_var)+inode,elem(j)%start%a(m_var)+jnode-1)
                      else if ( jdof >= idof ) then! gtype == csr_symm 
                         do k = a%gr%ia(idof),a%gr%ia(idof+1)-1
                            if ( a%gr%ja(k) == jdof ) exit
                         end do
                         assert ( k < a%gr%ia(idof+1) )
                         a%a(k) = a%a(k) + face%p_mat%a(elem(i)%start%a(l_var)+inode,elem(j)%start%a(m_var)+jnode-1)
                      end if
                   end do
                end if
             end do
          end do
       end do
    end do

  end subroutine face_element_matrix_assembly

  subroutine element_vector_assembly( dhand, elem, a, iblock )
    implicit none
    ! Parameters
    type(dof_handler_t), intent(in)             :: dhand
    type(fem_element_t), intent(in)             :: elem
    type(fem_vector_t), intent(inout)           :: a
    integer(ip), intent(in), optional         :: iblock

    integer(ip) :: iprob, nvapb_i, ivars, l_var, m_var
    integer(ip) :: inode, idof, iblock_, g_var

    iblock_ = 1
    if ( present(iblock) ) iblock_ = iblock
    iprob = elem%problem
    nvapb_i = dhand%prob_block(iblock_,iprob)%nd1
    do ivars = 1, nvapb_i
       l_var = dhand%prob_block(iblock_,iprob)%a(ivars)
       g_var = dhand%problems(iprob)%p%l2g_var(l_var)
       do inode = 1,elem%f_inf(l_var)%p%nnode
          idof = elem%elem2dof(inode,l_var)
          if ( idof  > 0 ) then
             a%b(idof) =  a%b(idof) + elem%p_vec%a(elem%start%a(l_var)+inode-1)
          end if
       end do
    end do

  end subroutine element_vector_assembly


  subroutine face_vector_assembly( dhand, elem, face, a, iblock )
    implicit none
    ! Parameters
    type(dof_handler_t), intent(in)             :: dhand
    type(fem_element_t), intent(in)             :: elem
    type(fem_face_t)   , intent(in)             :: face
    type(fem_vector_t), intent(inout)           :: a
    integer(ip), intent(in), optional         :: iblock

    integer(ip) :: iprob, nvapb_i, ivars, l_var, g_var
    integer(ip) :: inode, idof, iblock_, iobje, ndime

    iblock_ = 1
    ndime = elem%p_geo_info%ndime
    if ( present(iblock) ) iblock_ = iblock
    iprob = elem%problem
    nvapb_i = dhand%prob_block(iblock_,iprob)%nd1
    do ivars = 1, nvapb_i
       l_var = dhand%prob_block(iblock_,iprob)%a(ivars)
       g_var = dhand%problems(iprob)%p%l2g_var(l_var)
       iobje = face%face_object + elem%p_geo_info%nobje_dim(ndime) - 1
       do inode = elem%f_inf(l_var)%p%ntxob%p(iobje),elem%f_inf(l_var)%p%ntxob%p(iobje+1)-1
          idof = elem%elem2dof(elem%f_inf(l_var)%p%ntxob%l(inode),l_var)
          if ( idof  > 0 ) then
             a%b(idof) = a%b(idof) + face%p_vec%a(elem%start%a(l_var)+inode-1)
          end if
       end do
    end do

  end subroutine face_vector_assembly
 
end module assembly_names



