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
module fem_matrix_vector
  use types
  use fem_matrix_names
  use fem_vector_names
  implicit none
# include "debug.i90"

  private

  public :: fem_matvec, fem_matvec_trans, fem_matmat, fem_matmat_trans

contains

  subroutine fem_matvec (a,x,y)
    implicit none
    type(fem_matrix) , intent(in)    :: a
    type(fem_vector) , intent(in)    :: x
    type(fem_vector) , intent(inout) :: y
    real(rp) :: aux
    ! write (*,*) a%nd1, a%nd2, x%nd, y%nd 

    assert ( a%nd1 == y%nd )
    assert ( a%nd2 == x%nd )
    assert ( x%storage == y%storage )
    assert ( a%storage == x%storage )

    if (a%storage == blk) then ! gr must describe block sparsity pattern of mat
       if(a%type==css_mat) then
          call matvec_css(a%nd1,a%nd2,a%gr%nv,a%gr%ia,a%gr%is,a%gr%ja,a%d,a%l,a%u,x%b,y%b)
       else if(a%type==csr_mat) then
          if (a%symm == symm_false) then
             call matvec_csr(a%nd1,a%nd2,a%gr%nv,a%gr%nv2,a%gr%ia,a%gr%ja,a%a,x%b,y%b)
          else if (a%symm == symm_true) then
             call matvec_csr_symm(a%nd1,a%nd2,a%gr%nv,a%gr%nv,a%gr%ia,a%gr%ja,a%a,x%b,y%b)          
          end if
       else if(a%type==csc_mat) then
          call matvec_csc(a%nd1,a%nd2,a%gr%nv,a%gr%nv2,a%gr%ia,a%gr%ja,a%a,x%b,y%b)
       end if
    else if (a%storage==scal) then
       if(a%type==css_mat) then
          call matvec_css_scal(a%nd1,a%nd2,a%gr%nv,a%gr%ia,a%gr%is,a%gr%ja,a%d,a%l,a%u,x%b,y%b)
       else if(a%type==csr_mat) then
          if (a%symm == symm_false) then
             call matvec_csr_scal(a%nd1,a%nd2,a%gr%nv,a%gr%nv2,a%gr%ia,a%gr%ja,a%a,x%b,y%b)
          else if (a%symm == symm_true) then
             call matvec_csr_symm_scal(a%nd1,a%nd2,a%gr%nv,a%gr%nv,a%gr%ia,a%gr%ja,a%a,x%b,y%b)          
          end if
       else if(a%type==csc_mat) then
          call matvec_csc_scal(a%nd1,a%nd2,a%gr%nv,a%gr%nv2,a%gr%ia,a%gr%ja,a%a,x%b,y%b)
       end if
    end if

  end subroutine fem_matvec

  subroutine fem_matmat (a, n, ldX, x, ldY, y)
    implicit none

    ! Parameters
    type(fem_matrix) , intent(in)    :: a
    integer(ip)      , intent(in)    :: n
    integer(ip)      , intent(in)    :: ldX
    real(rp)         , intent(in)    :: x(ldX, n)
    integer(ip)      , intent(in)    :: ldY
    real(rp)         , intent(inout) :: y(ldY, n)

    ! Locals 
    integer (ip) :: i
 
    ! Version for blk storage layout still to be developped
    assert ( a%storage == scal ) 
!!$    if (a%storage == blk) then ! gr must describe block sparsity pattern of mat
!!$       if(a%type==css_mat) then
!!$          call matvec_css(a%nd1,a%nd2,a%gr%nv,a%gr%ia,a%gr%is,a%gr%ja,a%d,a%l,a%u,x%b,y%b)
!!$       else if(a%type==csr_mat) then
!!$          if (a%symm == symm_false) then
!!$             call matvec_csr(a%nd1,a%nd2,a%gr%nv,a%gr%nv2,a%gr%ia,a%gr%ja,a%a,x%b,y%b)
!!$          else if (a%symm == symm_true) then
!!$             call matvec_csr_symm(a%nd1,a%nd2,a%gr%nv,a%gr%nv,a%gr%ia,a%gr%ja,a%a,x%b,y%b)          
!!$          end if
!!$       else if(a%type==csc_mat) then
!!$          call matvec_csc(a%nd1,a%nd2,a%gr%nv,a%gr%nv2,a%gr%ia,a%gr%ja,a%a,x%b,y%b)
!!$       end if
!!$    else if (a%storage==scal) then
    do i=1,n
       if(a%type==css_mat) then
          call matvec_css_scal(a%nd1,a%nd2,a%gr%nv,a%gr%ia,a%gr%is,a%gr%ja,a%d,a%l,a%u,x(1:a%gr%nv2,i),y(1:a%gr%nv,i))
       else if(a%type==csr_mat) then
          if (a%symm == symm_false) then
             call matvec_csr_scal(a%nd1,a%nd2,a%gr%nv,a%gr%nv2,a%gr%ia,a%gr%ja,a%a,x(1:a%gr%nv2,i),y(1:a%gr%nv,i))
          else if (a%symm == symm_true) then
             call matvec_csr_symm_scal(a%nd1,a%nd2,a%gr%nv,a%gr%nv,a%gr%ia,a%gr%ja,a%a,x(1:a%gr%nv2,i),y(1:a%gr%nv,i))          
          end if
       else if(a%type==csc_mat) then
          call matvec_csc_scal(a%nd1,a%nd2,a%gr%nv,a%gr%nv2,a%gr%ia,a%gr%ja,a%a,x(1:a%gr%nv2,i),y(1:a%gr%nv,i))
       end if
    end do
!!$    end if



  end subroutine fem_matmat

  subroutine fem_matvec_trans (a,x,y)
    implicit none
    type(fem_matrix) , intent(in)    :: a
    type(fem_vector) , intent(in)    :: x
    type(fem_vector) , intent(inout) :: y

    assert ( a%nd1 == y%nd)
    assert ( a%nd2 == x%nd)
    assert ( x%storage == y%storage )
    assert ( a%storage == x%storage )

    if (a%storage == blk) then ! gr must describe block sparsity pattern of mat
       if(a%type==css_mat) then
          ! call matvec_css_trans(a%nd1,a%nd2,a%gr%nv,a%gr%ia,a%gr%is,a%gr%ja,a%d,a%l,a%u,x%b,y%b)
          write (0,*) 'Error: the body of matvec_css_trans in matvec_dof.i90 still to be written'
          write (0,*) 'Error: volunteers are welcome !!!'
          stop
       else if(a%type==csr_mat) then
          if (a%symm == symm_false) then
             call matvec_csr_trans(a%nd1,a%nd2,a%gr%nv,a%gr%nv2,a%gr%ia,a%gr%ja,a%a,x%b,y%b)
          else if (a%symm == symm_true) then
             ! call matvec_csr_symm_trans(a%nd1,a%nd2,a%gr%nv,a%gr%nv,a%gr%ia,a%gr%ja,a%a,x%b,y%b)
             write (0,*) 'Error: the body of matvec_csr_symm_trans in matvec_dof.i90 still to be written'
             write (0,*) 'Error: volunteers are welcome !!!'
             stop
          end if
       else if(a%type==csc_mat) then
          ! call matvec_csc_trans (a%nd1,a%nd2,a%gr%nv,a%gr%nv,a%gr%ia,a%gr%ja,a%a,x%b,y%b)
          write (0,*) 'Error: the body of matvec_csc_trans in matvec_dof.i90 still to be written'
          write (0,*) 'Error: volunteers are welcome !!!'
          stop
       end if
    else if (a%storage==scal) then
       if(a%type==css_mat) then
          ! call matvec_css_scal_trans(a%nd1,a%nd2,a%gr%nv,a%gr%ia,a%gr%is,a%gr%ja,a%d,a%l,a%u,x%b,y%b)
          write (0,*) 'Error: the body of matvec_css_scal_trans in matvec_dof.i90 still to be written'
          write (0,*) 'Error: volunteers are welcome !!!'
          stop          
       else if(a%type==csr_mat) then
          if (a%symm == symm_false) then
             call matvec_csr_trans_scal(a%nd1,a%nd2,a%gr%nv,a%gr%nv2,a%gr%ia,a%gr%ja,a%a,x%b,y%b)
          else if (a%symm == symm_true) then
             call matvec_csr_symm_trans_scal(a%nd1,a%nd2,a%gr%nv,a%gr%nv,a%gr%ia,a%gr%ja,a%a,x%b,y%b)          
          end if
       else if(a%type==csc_mat) then
          ! call matvec_csc_trans_scal(a%nd1,a%nd2,a%gr%nv,a%gr%nv,a%gr%ia,a%gr%ja,a%a,x%b,y%b)
          write (0,*) 'Error: the body of matvec_csc_trans_scal in matvec_dof.i90 still to be written'
          write (0,*) 'Error: volunteers are welcome !!!'
          stop
       end if
    end if
  end subroutine fem_matvec_trans

  subroutine fem_matmat_trans (a, n, ldX, x, ldY, y)
    implicit none

    ! Parameters
    type(fem_matrix) , intent(in)    :: a
    integer(ip)      , intent(in)    :: n
    integer(ip)      , intent(in)    :: ldX
    real(rp)         , intent(in)    :: x(ldX, n)
    integer(ip)      , intent(in)    :: ldY
    real(rp)         , intent(inout) :: y(ldY, n)

    ! Locals 
    integer (ip) :: i

    ! Version for blk storage layout still to be developped
    assert ( a%storage == scal ) 

!!$    if (a%storage == blk) then ! gr must describe block sparsity pattern of mat
!!$       if(a%type==css_mat) then
!!$          ! call matvec_css_trans(a%nd1,a%nd2,a%gr%nv,a%gr%ia,a%gr%is,a%gr%ja,a%d,a%l,a%u,x%b,y%b)
!!$          write (0,*) 'Error: the body of matvec_css_trans in matvec_dof.i90 still to be written'
!!$          write (0,*) 'Error: volunteers are welcome !!!'
!!$          stop
!!$       else if(a%type==csr_mat) then
!!$          if (a%symm == symm_false) then
!!$             call matvec_csr_trans(a%nd1,a%nd2,a%gr%nv,a%gr%nv2,a%gr%ia,a%gr%ja,a%a,x%b,y%b)
!!$          else if (a%symm == symm_true) then
!!$             ! call matvec_csr_symm_trans(a%nd1,a%nd2,a%gr%nv,a%gr%nv,a%gr%ia,a%gr%ja,a%a,x%b,y%b)
!!$             write (0,*) 'Error: the body of matvec_csr_symm_trans in matvec_dof.i90 still to be written'
!!$             write (0,*) 'Error: volunteers are welcome !!!'
!!$             stop
!!$          end if
!!$       else if(a%type==csc_mat) then
!!$          ! call matvec_csc_trans (a%nd1,a%nd2,a%gr%nv,a%gr%nv,a%gr%ia,a%gr%ja,a%a,x%b,y%b)
!!$          write (0,*) 'Error: the body of matvec_csc_trans in matvec_dof.i90 still to be written'
!!$          write (0,*) 'Error: volunteers are welcome !!!'
!!$          stop
!!$       end if
!!$    else if (a%storage==scal) then

    do i=1,n
       if(a%type==css_mat) then
          ! call matvec_css_scal_trans(a%nd1,a%nd2,a%gr%nv,a%gr%ia,a%gr%is,a%gr%ja,a%d,a%l,a%u,x%b,y%b)
          write (0,*) 'Error: the body of matvec_css_scal_trans in matvec_dof.i90 still to be written'
          write (0,*) 'Error: volunteers are welcome !!!'
          stop          
       else if(a%type==csr_mat) then
          if (a%symm == symm_false) then
             call matvec_csr_trans_scal(a%nd1,a%nd2,a%gr%nv,a%gr%nv2,a%gr%ia,a%gr%ja,a%a,x(1:a%gr%nv2,i),y(1:a%gr%nv,i) )
          else if (a%symm == symm_true) then
             call matvec_csr_symm_trans_scal(a%nd1,a%nd2,a%gr%nv,a%gr%nv,a%gr%ia,a%gr%ja,a%a,x(1:a%gr%nv,i),y(1:a%gr%nv2,i))          
          end if
       else if(a%type==csc_mat) then
          ! call matvec_csc_trans_scal(a%nd1,a%nd2,a%gr%nv,a%gr%nv,a%gr%ia,a%gr%ja,a%a,x%b,y%b)
          write (0,*) 'Error: the body of matvec_csc_trans_scal in matvec_dof.i90 still to be written'
          write (0,*) 'Error: volunteers are welcome !!!'
          stop
       end if
    end do

!!$    end if

  end subroutine fem_matmat_trans

end module fem_matrix_vector
