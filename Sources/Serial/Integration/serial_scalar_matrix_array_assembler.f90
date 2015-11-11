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
module serial_scalar_matrix_array_assembler_names
  use types_names
  use finite_element_names
  use dof_descriptor_names
  use allocatable_array_names
 
  ! Abstract modules
  use matrix_array_assembler_names
  use matrix_names
  use array_names
  
  ! Concrete implementations
  use serial_scalar_matrix_names
  use serial_scalar_array_names
  
  implicit none
# include "debug.i90"
  private

  type, extends(matrix_array_assembler_t) :: serial_scalar_matrix_array_assembler_t
  contains
	procedure :: assembly => serial_scalar_matrix_array_assembler_assembly
  end type
	 
  ! Data types
  public :: serial_scalar_matrix_array_assembler_t
  public :: element_serial_scalar_matrix_assembly, element_serial_scalar_array_assembly, & 
            face_element_serial_scalar_matrix_assembly, face_element_serial_scalar_array_assembly
  
contains
  subroutine serial_scalar_matrix_array_assembler_assembly(this,dof_descriptor,finite_element) 
    implicit none
    class(serial_scalar_matrix_array_assembler_t), intent(inout) :: this
    type(dof_descriptor_t)                       , intent(in)    :: dof_descriptor
    type(finite_element_t)                       , intent(in)    :: finite_element
	
	class(matrix_t), pointer :: matrix
	class(array_t) , pointer :: array
	
	matrix => this%get_matrix()
	array  => this%get_array()
	
	select type(matrix)
    class is(serial_scalar_matrix_t)
	   call element_serial_scalar_matrix_assembly( dof_descriptor, finite_element, matrix )
	 class default
       check(.false.)
    end select  
	
    select type(array)
    class is(serial_scalar_array_t)
	   call element_serial_scalar_array_assembly( dof_descriptor, finite_element, array )
	 class default
       check(.false.)
    end select 
  end subroutine serial_scalar_matrix_array_assembler_assembly
  
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

    do ivars = 1, nvapb_i
       l_var = dof_descriptor%prob_block(iblock_,iprob)%a(ivars)
       g_var = dof_descriptor%problems(iprob)%p%l2g_var(l_var)
       do jvars = 1, nvapb_j
          m_var = dof_descriptor%prob_block(jblock_,iprob)%a(jvars)
          k_var = dof_descriptor%problems(iprob)%p%l2g_var(m_var)
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
  end subroutine element_serial_scalar_matrix_assembly
  
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
  
  subroutine assembly_face_element_matrix_mono(  fe_face, finite_element, dof_descriptor, a ) 
    implicit none
    type(dof_descriptor_t), intent(in)    :: dof_descriptor
    type(fe_face_t)       , intent(in)    :: fe_face
    type(finite_element_pointer_t), intent(in)    :: finite_element(2)
    type(serial_scalar_matrix_t)        , intent(inout) :: a

    integer(ip) :: i
    type(allocatable_array_ip1_t) :: start_aux

    ! Auxiliar local start array to store start(2)
    call start_aux%create(dof_descriptor%problems(finite_element(2)%p%problem)%p%nvars+1)
    start_aux%a = finite_element(2)%p%start%a

    finite_element(2)%p%start%a = finite_element(2)%p%start%a + finite_element(1)%p%start%a(dof_descriptor%problems(finite_element(1)%p%problem)%p%nvars+1) - 1
    do i = 1,2
       call element_serial_scalar_matrix_assembly( dof_descriptor, finite_element(i)%p, a )
    end do
    call face_element_serial_scalar_matrix_assembly( dof_descriptor, finite_element, fe_face, a )

    ! Restore start(2)
    finite_element(2)%p%start%a = start_aux%a
  end subroutine assembly_face_element_matrix_mono
  
  subroutine face_element_serial_scalar_matrix_assembly( dof_descriptor, finite_element, fe_face, a, iblock, jblock )
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

  end subroutine face_element_serial_scalar_matrix_assembly
  
  subroutine face_element_serial_scalar_array_assembly( dof_descriptor, finite_element, fe_face, a, iblock )
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

  end subroutine face_element_serial_scalar_array_assembly
end module serial_scalar_matrix_array_assembler_names
