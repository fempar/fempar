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
!=============================================================================
! Auxiliar modules to handle ndofs by a template like technique
!=============================================================================

module matvec_names
  use types_names
  implicit none

# include "debug.i90"
  private

  public :: matvec, matvec_trans, matvec_symmetric_storage, matvec_symmetric_storage_trans
  
contains

  ! Debugged
  subroutine matvec (nv,nv2,ia,ja,a,x,y)
    implicit none
    integer(ip), intent(in)  :: nv,nv2,ia(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: a(ia(nv+1)-1),x(nv2)
    real(rp)   , intent(out) :: y(nv)
    integer(ip)              :: iv,iz,jv

    y = 0.0_rp
    do iv = 1, nv
       do iz = ia(iv), ia(iv+1)-1
          jv   = ja(iz)
          y(iv) = y(iv) + x(jv)*a(iz)
       end do ! iz
    end do ! iv

  end subroutine matvec

  ! Debugged
  subroutine matvec_symmetric_storage (nv,nv2,ia,ja,a,x,y)
    implicit none
    integer(ip), intent(in)  :: nv,nv2,ia(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: a(ia(nv+1)-1),x(nv2)
    real(rp)   , intent(out) :: y(nv)
    integer(ip)              :: iv,iz,jv

    assert(nv==nv2)

    y = 0.0_rp
    do iv = 1, nv
       y(iv) = y(iv) + x(ja(ia(iv)))*a(ia(iv))
       do iz = ia(iv)+1, ia(iv+1)-1
          jv = ja(iz)
          y(iv) = y(iv) + x(jv)*a(iz)
          y(jv) = y(jv) + x(iv)*a(iz)
       end do ! iz
    end do ! iv

  end subroutine matvec_symmetric_storage

  ! Debugged
  subroutine matvec_trans (nv,nv2,ia,ja,a,x,y)
    implicit none
    integer(ip), intent(in)  :: nv,nv2,ia(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: a(ia(nv+1)-1), x(nv)
    real(rp)   , intent(out) :: y(nv2)
    integer(ip)              :: iv,iz,jv
    
    y = 0.0_rp
    do iv = 1, nv
       do iz = ia(iv), ia(iv+1)-1
          jv = ja(iz)
          y(jv) = y(jv) + x(iv)*a(iz)
       end do ! iz
    end do ! iv
    
  end subroutine matvec_trans

  subroutine matvec_symmetric_storage_trans (nv,nv2,ia,ja,a,x,y)
    implicit none
    integer(ip), intent(in)  :: nv,nv2,ia(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: a(ia(nv+1)-1), x(nv)
    real(rp)   , intent(out) :: y(nv2)
    integer(ip)              :: iv,iz,jv,id,jd
    integer(ip)              :: of, ivc, izc, jvc

    write (0,*) 'Error: the body of matvec_symmetric_storage_trans in matvec.f90 still to be written'
    write (0,*) 'Error: volunteers are welcome !!!'
    stop 
  end subroutine matvec_symmetric_storage_trans



end module matvec_names




       
