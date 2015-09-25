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
  use allocatable_array_names
  use finite_element_names
  !use fe_space_names
  use integrable_names
  use dof_descriptor_names
  use serial_block_matrix_names
  use serial_scalar_matrix_names
  use serial_block_array_names
  use serial_scalar_array_names
  use serial_scalar_names
  use graph_names

  implicit none
# include "debug.i90"
  private

  ! Functions
  public :: assembly, assembly_face, element_serial_scalar_matrix_assembly, element_serial_scalar_array_assembly

contains

  subroutine assembly(finite_element, dof_descriptor, a ) 
    implicit none
    ! Parameters
    type(dof_descriptor_t), intent(in)    :: dof_descriptor
    type(finite_element_t), intent(in)    :: finite_element
    class(integrable_t)   , intent(inout) :: a

    select type(a)
    class is(serial_scalar_matrix_t)
       call assembly_element_matrix_mono(finite_element, dof_descriptor, a) 
    class is(serial_scalar_array_t)
       call assembly_element_vector_mono(finite_element, dof_descriptor, a)
    class is(serial_block_matrix_t)
       call assembly_element_matrix_block(finite_element, dof_descriptor, a)
    class is(serial_block_array_t)
       call assembly_element_vector_block(finite_element, dof_descriptor, a)
    class default
       ! class not yet implemented
       check(.false.)
    end select
  end subroutine assembly

  subroutine assembly_face(fe_face, finite_elements, dof_descriptor, a ) 
    implicit none
    ! Parameters
    type(dof_descriptor_t), intent(in)    :: dof_descriptor
    type(fe_face_t)       , intent(in)    :: fe_face
    type(finite_element_pointer_t), intent(in)    :: finite_elements(2)
    class(integrable_t)   , intent(inout) :: a

    select type(a)
    class is(serial_scalar_matrix_t)
       call assembly_face_element_matrix_mono(fe_face, finite_elements, dof_descriptor, a) 
    !class is(vector_t)
    !   call assembly_face_element_vector_mono(fe_face, dof_descriptor, a)
    class is(serial_block_matrix_t)
       call assembly_face_element_matrix_block(fe_face, finite_elements, dof_descriptor, a)
    !class is(block_vector_t)
    !   call assembly_face_element_vector_block(fe_face, dof_descriptor, a)
    class default
       ! class not yet implemented
       check(.false.)
    end select
  end subroutine assembly_face

  subroutine assembly_element_matrix_block(  finite_element, dof_descriptor, a ) 
    implicit none
    type(dof_descriptor_t), intent(in)    :: dof_descriptor
    type(finite_element_t), intent(in)    :: finite_element
    type(serial_block_matrix_t)  , intent(inout) :: a
    
    integer(ip) :: ivar, iblock, jblock

    type(serial_scalar_matrix_t), pointer :: f_matrix

    do iblock = 1, dof_descriptor%nblocks
       do jblock = 1, dof_descriptor%nblocks
          f_matrix => a%blocks(iblock,jblock)%serial_scalar_matrix
          if ( associated(f_matrix) ) then
             call element_serial_scalar_matrix_assembly( dof_descriptor, finite_element, f_matrix, iblock, jblock )
          end if 
       end do
    end do

  end subroutine assembly_element_matrix_block

  subroutine assembly_element_matrix_mono(  finite_element, dof_descriptor, a ) 
    implicit none
    type(dof_descriptor_t), intent(in)    :: dof_descriptor
    type(finite_element_t), intent(in)    :: finite_element
    type(serial_scalar_matrix_t)        , intent(inout) :: a

    call element_serial_scalar_matrix_assembly( dof_descriptor, finite_element, a )

  end subroutine assembly_element_matrix_mono

  subroutine assembly_face_element_matrix_block(  fe_face, finite_element, dof_descriptor, a ) 
    implicit none
    type(dof_descriptor_t), intent(in)    :: dof_descriptor
    type(fe_face_t)       , intent(in)    :: fe_face
    type(finite_element_pointer_t), intent(in)    :: finite_element(2)
    type(serial_block_matrix_t)  , intent(inout) :: a

    integer(ip) :: iblock, jblock, i
    type(allocatable_array_ip1_t) :: start_aux
    type(serial_scalar_matrix_t), pointer :: f_matrix

    ! Auxiliar local start array to store start(2)
    call allocatable_array_create(dof_descriptor%problems(finite_element(2)%p%problem)%p%nvars+1,start_aux)
    start_aux%a = finite_element(2)%p%start%a

    finite_element(2)%p%start%a = finite_element(2)%p%start%a + finite_element(1)%p%start%a(dof_descriptor%problems(finite_element(1)%p%problem)%p%nvars+1) - 1

    do iblock = 1, dof_descriptor%nblocks
       do jblock = 1, dof_descriptor%nblocks
          f_matrix => a%blocks(iblock,jblock)%serial_scalar_matrix
          if ( associated(f_matrix) ) then
            do i = 1,2
               call element_serial_scalar_matrix_assembly( dof_descriptor, finite_element(i)%p, f_matrix, iblock, jblock )
            end do
            call face_element_matrix_assembly( dof_descriptor, finite_element, fe_face, f_matrix, iblock, jblock )
          end if
       end do
    end do

    ! Restore start(2)
    finite_element(2)%p%start%a = start_aux%a

  end subroutine assembly_face_element_matrix_block

  subroutine assembly_face_element_matrix_mono(  fe_face, finite_element, dof_descriptor, a ) 
    implicit none
    type(dof_descriptor_t), intent(in)    :: dof_descriptor
    type(fe_face_t)       , intent(in)    :: fe_face
    type(finite_element_pointer_t), intent(in)    :: finite_element(2)
    type(serial_scalar_matrix_t)        , intent(inout) :: a

    integer(ip) :: i
    type(allocatable_array_ip1_t) :: start_aux

    ! Auxiliar local start array to store start(2)
    call allocatable_array_create(dof_descriptor%problems(finite_element(2)%p%problem)%p%nvars+1,start_aux)
    start_aux%a = finite_element(2)%p%start%a

    finite_element(2)%p%start%a = finite_element(2)%p%start%a + finite_element(1)%p%start%a(dof_descriptor%problems(finite_element(1)%p%problem)%p%nvars+1) - 1
    do i = 1,2
       call element_serial_scalar_matrix_assembly( dof_descriptor, finite_element(i)%p, a )
    end do
    call face_element_matrix_assembly( dof_descriptor, finite_element, fe_face, a )

    ! Restore start(2)
    finite_element(2)%p%start%a = start_aux%a

  end subroutine assembly_face_element_matrix_mono

  subroutine assembly_element_vector_block(  finite_element, dof_descriptor, a ) 
    implicit none
    type(dof_descriptor_t), intent(in)    :: dof_descriptor
    type(finite_element_t), intent(in)    :: finite_element
    type(serial_block_array_t)  , intent(inout) :: a

    integer(ip) :: iblock

    do iblock = 1, dof_descriptor%nblocks
       call element_serial_scalar_array_assembly( dof_descriptor, finite_element, a%blocks(iblock), &
            & iblock )
    end do

  end subroutine assembly_element_vector_block

  subroutine assembly_element_vector_mono(  finite_element, dof_descriptor, a ) 
    implicit none
    type(dof_descriptor_t), intent(in)    :: dof_descriptor
    type(finite_element_t), intent(in)    :: finite_element
    type(serial_scalar_array_t)        , intent(inout) :: a
    
    call element_serial_scalar_array_assembly( dof_descriptor, finite_element, a )

  end subroutine assembly_element_vector_mono

  subroutine assembly_face_vector_block(  fe_face, finite_element, dof_descriptor, a ) 
    implicit none
    type(dof_descriptor_t), intent(in)    :: dof_descriptor
    type(fe_face_t)       , intent(in)    :: fe_face
    type(finite_element_t), intent(in)    :: finite_element
    type(serial_block_array_t)  , intent(inout) :: a

    integer(ip) :: iblock

    ! Note: This subroutine only has sense on interface / boundary faces, only
    ! related to ONE element.
    do iblock = 1, dof_descriptor%nblocks
       call face_vector_assembly( dof_descriptor, finite_element, fe_face, a%blocks(iblock), iblock )
    end do

  end subroutine assembly_face_vector_block

  subroutine assembly_face_vector_mono(  fe_face, finite_element, dof_descriptor, a ) 
    implicit none
    type(dof_descriptor_t), intent(in)    :: dof_descriptor
    type(fe_face_t)       , intent(in)    :: fe_face
    type(finite_element_t), intent(in)    :: finite_element
    type(serial_scalar_array_t)        , intent(inout) :: a

    ! Note: This subroutine only has sense on interface / boundary faces, only
    ! related to ONE element.
    call face_vector_assembly( dof_descriptor, finite_element, fe_face, a )

  end subroutine assembly_face_vector_mono

  subroutine element_serial_scalar_matrix_assembly( dof_descriptor, finite_element, a, iblock, jblock )
    implicit none
    ! Parameters
    type(dof_descriptor_t)    , intent(in)    :: dof_descriptor
    type(finite_element_t) , intent(in)    :: finite_element
    type(serial_scalar_matrix_t)         , intent(inout) :: a
    integer(ip), intent(in), optional      :: iblock, jblock

    integer(ip) :: iprob, nvapb_i, nvapb_j, ivars, jvars, l_var, m_var, g_var, k_var
    integer(ip) :: inode, jnode, idof, jdof, k, iblock_, jblock_

    iblock_ = 1
    jblock_ = 1
    if ( present(iblock) ) iblock_ = iblock
    if ( present(jblock) ) jblock_ = jblock
        
    iprob = finite_element%problem
	
    nvapb_i = dof_descriptor%prob_block(iblock_,iprob)%nd1
    nvapb_j = dof_descriptor%prob_block(jblock_,iprob)%nd1

    !write (*,*) 'nvapb_i:',nvapb_i
    !write (*,*) 'nvapb_j:',nvapb_j
    !write (*,*) 'start:',elem%start%a

    !write(*,*) 'local matrix',finite_element%p_mat%nd1, finite_element%p_mat%nd2
    !write(*,*) 'local matrix',finite_element%p_mat%a
    


    do ivars = 1, nvapb_i
       l_var = dof_descriptor%prob_block(iblock_,iprob)%a(ivars)
       g_var = dof_descriptor%problems(iprob)%p%l2g_var(l_var)
       do jvars = 1, nvapb_j
          m_var = dof_descriptor%prob_block(jblock_,iprob)%a(jvars)
          k_var = dof_descriptor%problems(iprob)%p%l2g_var(m_var)
          !write (*,*) 'l_var:',l_var
          !write (*,*) 'm_var:',m_var
          if( dof_descriptor%dof_coupl(g_var,k_var) == 1) then
             do inode = 1,finite_element%reference_element_vars(l_var)%p%nnode
                idof = finite_element%elem2dof(inode,l_var)
                if ( idof  > 0 ) then
                   do jnode = 1,finite_element%reference_element_vars(m_var)%p%nnode
                      jdof = finite_element%elem2dof(jnode,m_var)
                      if (  (.not. a%graph%symmetric_storage) .and. jdof > 0 ) then
                         do k = a%graph%ia(idof),a%graph%ia(idof+1)-1
                            if ( a%graph%ja(k) == jdof ) exit
                         end do
                         assert ( k < a%graph%ia(idof+1) )
                         a%a(k) = a%a(k) + finite_element%p_mat%a(finite_element%start%a(l_var)+inode-1,finite_element%start%a(m_var)+jnode-1)
                      else if ( jdof >= idof ) then 
                         do k = a%graph%ia(idof),a%graph%ia(idof+1)-1
                            if ( a%graph%ja(k) == jdof ) exit
                         end do
                         assert ( k < a%graph%ia(idof+1) )
                         a%a(k) = a%a(k) + finite_element%p_mat%a(finite_element%start%a(l_var)+inode-1,finite_element%start%a(m_var)+jnode-1)
                      end if
                   end do
                end if
             end do
          end if
       end do
    end do
    
    !write(*,*) 'system_matrix'
    !call matrix_print(6,a)

  end subroutine element_serial_scalar_matrix_assembly

  subroutine face_element_matrix_assembly( dof_descriptor, finite_element, fe_face, a, iblock, jblock )
    implicit none
    ! Parameters
    type(dof_descriptor_t), intent(in)    :: dof_descriptor
    type(finite_element_pointer_t), intent(in)    :: finite_element(2)
    type(fe_face_t)       , intent(in)    :: fe_face
    type(serial_scalar_matrix_t)        , intent(inout) :: a
    integer(ip), intent(in), optional :: iblock, jblock

    integer(ip) :: iprob, jprob, nvapb_i, nvapb_j, ivars, jvars, l_var, m_var, g_var, k_var
    integer(ip) :: inode, jnode, idof, jdof, k, iblock_, jblock_, iobje, i, j, ndime

    iblock_ = 1
    jblock_ = 1
    if ( present(iblock) ) iblock_ = iblock
    if ( present(jblock) ) jblock_ = jblock
    do i = 1, 2
       j = 3 - i
       iprob = finite_element(i)%p%problem
       jprob = finite_element(j)%p%problem
       nvapb_i = dof_descriptor%prob_block(iblock_,iprob)%nd1
       nvapb_j = dof_descriptor%prob_block(jblock_,jprob)%nd1
       do ivars = 1, nvapb_i
          l_var = dof_descriptor%prob_block(iblock_,iprob)%a(ivars)
          g_var = dof_descriptor%problems(iprob)%p%l2g_var(l_var)
          do jvars = 1, nvapb_j
             m_var = dof_descriptor%prob_block(jblock_,jprob)%a(jvars)
             k_var = dof_descriptor%problems(jprob)%p%l2g_var(m_var)
             do inode = 1,finite_element(i)%p%reference_element_vars(l_var)%p%nnode
                idof = finite_element(i)%p%elem2dof(inode,l_var)
                if ( idof  > 0 ) then
                   ndime = finite_element(j)%p%p_geo_reference_element%ndime
                   iobje = fe_face%face_vef + finite_element(j)%p%p_geo_reference_element%nvef_dim(ndime) - 1
                   do jnode = finite_element(j)%p%reference_element_vars(m_var)%p%ntxob%p(iobje),finite_element(j)%p%reference_element_vars(m_var)%p%ntxob%p(iobje+1)-1
                      jdof = finite_element(j)%p%elem2dof(finite_element(j)%p%reference_element_vars(m_var)%p%ntxob%l(jnode),m_var)
                      if (  (.not. a%graph%symmetric_storage) .and. jdof > 0 ) then
                         do k = a%graph%ia(idof),a%graph%ia(idof+1)-1
                            if ( a%graph%ja(k) == jdof ) exit
                         end do
                         assert ( k < a%graph%ia(idof+1) )
                         a%a(k) = a%a(k) + fe_face%p_mat%a(finite_element(i)%p%start%a(l_var)+inode,finite_element(j)%p%start%a(m_var)+jnode-1)
                      else if ( jdof >= idof ) then 
                         do k = a%graph%ia(idof),a%graph%ia(idof+1)-1
                            if ( a%graph%ja(k) == jdof ) exit
                         end do
                         assert ( k < a%graph%ia(idof+1) )
                         a%a(k) = a%a(k) + fe_face%p_mat%a(finite_element(i)%p%start%a(l_var)+inode,finite_element(j)%p%start%a(m_var)+jnode-1)
                      end if
                   end do
                end if
             end do
          end do
       end do
    end do

  end subroutine face_element_matrix_assembly

  subroutine element_serial_scalar_array_assembly( dof_descriptor, finite_element, a, iblock )
    implicit none
    ! Parameters
    type(dof_descriptor_t), intent(in)    :: dof_descriptor
    type(finite_element_t), intent(in)    :: finite_element
    type(serial_scalar_array_t)        , intent(inout) :: a
    integer(ip), intent(in), optional :: iblock

    integer(ip) :: iprob, nvapb_i, ivars, l_var, m_var
    integer(ip) :: inode, idof, iblock_, g_var

    iblock_ = 1
    if ( present(iblock) ) iblock_ = iblock
    iprob = finite_element%problem
    nvapb_i = dof_descriptor%prob_block(iblock_,iprob)%nd1
    do ivars = 1, nvapb_i
       l_var = dof_descriptor%prob_block(iblock_,iprob)%a(ivars)
       g_var = dof_descriptor%problems(iprob)%p%l2g_var(l_var)
       do inode = 1,finite_element%reference_element_vars(l_var)%p%nnode
          idof = finite_element%elem2dof(inode,l_var)
          if ( idof  > 0 ) then
             a%b(idof) =  a%b(idof) + finite_element%p_vec%a(finite_element%start%a(l_var)+inode-1)
          end if
       end do
    end do

  end subroutine element_serial_scalar_array_assembly

  subroutine face_vector_assembly( dof_descriptor, finite_element, fe_face, a, iblock )
    implicit none
    ! Parameters
    type(dof_descriptor_t), intent(in)    :: dof_descriptor
    type(finite_element_t), intent(in)    :: finite_element
    type(fe_face_t)       , intent(in)    :: fe_face
    type(serial_scalar_array_t)        , intent(inout) :: a
    integer(ip), intent(in), optional :: iblock

    integer(ip) :: iprob, nvapb_i, ivars, l_var, g_var
    integer(ip) :: inode, idof, iblock_, iobje, ndime

    iblock_ = 1
    ndime = finite_element%p_geo_reference_element%ndime
    if ( present(iblock) ) iblock_ = iblock
    iprob = finite_element%problem
    nvapb_i = dof_descriptor%prob_block(iblock_,iprob)%nd1
    do ivars = 1, nvapb_i
       l_var = dof_descriptor%prob_block(iblock_,iprob)%a(ivars)
       g_var = dof_descriptor%problems(iprob)%p%l2g_var(l_var)
       iobje = fe_face%face_vef + finite_element%p_geo_reference_element%nvef_dim(ndime) - 1
       do inode = finite_element%reference_element_vars(l_var)%p%ntxob%p(iobje),finite_element%reference_element_vars(l_var)%p%ntxob%p(iobje+1)-1
          idof = finite_element%elem2dof(finite_element%reference_element_vars(l_var)%p%ntxob%l(inode),l_var)
          if ( idof  > 0 ) then
             a%b(idof) = a%b(idof) + fe_face%p_vec%a(finite_element%start%a(l_var)+inode-1)
          end if
       end do
    end do

  end subroutine face_vector_assembly
 
end module assembly_names



