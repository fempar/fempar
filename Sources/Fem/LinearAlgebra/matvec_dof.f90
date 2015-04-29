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

module matvec_dof

  use types

  integer(ip), parameter :: ndof1 = 1, ndof2 = 1

  contains


subroutine matvec_css(nv, ia, is, ja, da, la, ua, x, y)
  implicit none
  integer(ip), intent(in)  :: nv,ia(nv+1),is(nv+1),ja(ia(nv+1)-1)
  real(rp)   , intent(in)  :: x(nv), da(nv)
  real(rp)   , intent(in)  :: la(is(nv+1)-1), ua(is(nv+1)-1)
  real(rp)   , intent(out) :: y(nv)

  write (0,*) 'Error: the body of matvec_css in matvec_dof.i90 still to be written'
  write (0,*) 'Error: volunteers are welcome !!!'
  stop

end subroutine matvec_css


! Debugged
subroutine matvec_csr (nv,nv2,ia,ja,a,x,y)
  implicit none
  integer(ip), intent(in)  :: nv,nv2,ia(nv+1),ja(ia(nv+1)-1)
  real(rp)   , intent(in)  :: a(ia(nv+1)-1),x(nv2)
  real(rp)   , intent(out) :: y(nv)
  integer(ip)              :: iv,iz,jv,id,jd
  integer(ip)              :: of, ivc, izc, jvc
  y = 0.0_rp
  do iv = 1, nv, ndof1
     do iz = ia(iv), ia(iv+1)-1, ndof2
        of   = iz - ia(iv)
        jv   = ja(iz)
        ! write (*,*) iv, jv, of !  DBG:
        do ivc = iv, iv + ndof1 - 1
           jvc   = jv
           do izc = ia(ivc)+of, ia(ivc)+of+ndof2-1
              ! write (*,*) 'XXX', ivc, jvc, ia(ivc)+of, ia(ivc)+of+ndof2-1 ! DBG:
              y(ivc) = y(ivc) + x(jvc)*a(izc)
              jvc      = jvc + 1
           end do ! izc
        end do ! ivc
     end do ! iz
  end do ! iv
end subroutine matvec_csr



! Debugged
subroutine matvec_csr_symm (nv,nv2,ia,ja,a,x,y)
  implicit none
  integer(ip), intent(in)  :: nv,nv2,ia(nv+1),ja(ia(nv+1)-1)
  real(rp)   , intent(in)  :: a(ia(nv+1)-1),x(nv2)
  real(rp)   , intent(out) :: y(nv)
  integer(ip)              :: iv,iz,jv,id,jd
  integer(ip)              :: of, ivc, izc, jvc

  y = 0.0_rp
  do iv = 1, nv, ndof1
     
     ! Diagonal block entry
     iz  = ia(iv)
     do ivc = iv, iv + ndof1 - 1
        jvc   = ivc
        izc   = ia(ivc)
        
        ! Diagonal scalar entries
        y(ivc) = y(ivc) + x(jvc)*a(izc)
        
        jvc = ivc+1
        ! Off-diagonal scalar entries
        do izc = ia(ivc)+1, ia(ivc)+ndof2-(ivc-iv)-1
           y(ivc) = y(ivc) + x(jvc)*a(izc)
           y(jvc) = y(jvc) + x(ivc)*a(izc)
           ! write (*,*) 'XXX', ivc, jvc, ia(ivc), ia(ivc)+ndof2-(ivc-iv)-1   ! DBG:
           jvc      = jvc + 1 
        end do ! izc
     end do ! ivc

     ! Strict upper-triangular block entries
     do iz = ia(iv)+ndof2, ia(iv+1)-1, ndof2
        of   = iz - ia(iv)
        jv = ja(iz)
        do ivc = iv, iv + ndof1 - 1
           jvc   = jv
           do izc = ia(ivc)+of, ia(ivc)+of+ndof2-1
              y(ivc) = y(ivc) + x(jvc)*a(izc)
              y(jvc) = y(jvc) + x(ivc)*a(izc)
              jvc      = jvc + 1 
           end do ! izc
           of = of - 1        ! Next global row has an entry less
        end do ! ivc
     end do ! iz
  end do ! iv

end subroutine matvec_csr_symm

! Debugged
subroutine matvec_csr_trans (nv,nv2,ia,ja,a,x,y)
  implicit none
  integer(ip), intent(in)  :: nv,nv2,ia(nv+1),ja(ia(nv+1)-1)
  real(rp)   , intent(in)  :: a(ia(nv+1)-1), x(nv)
  real(rp)   , intent(out) :: y(nv2)
  integer(ip)              :: iv,iz,jv,id,jd
  integer(ip)              :: of, ivc, izc, jvc

  y = 0.0_rp
  do iv = 1, nv, ndof1
     do iz = ia(iv), ia(iv+1)-1, ndof2
        of   = iz - ia(iv)
        ivc  = iv
        jv = ja(iz)
        do ivc = ivc, ivc + ndof1 - 1
           izc   = ia(ivc) + of
           jvc   = jv
           do izc = izc, izc + ndof2 - 1
              y(jvc) = y(jvc) + x(ivc)*a(izc)
              jvc      = jvc + 1 
           end do ! izc
        end do ! ivc
     end do ! iz
  end do ! iv

end subroutine matvec_csr_trans

subroutine matvec_csr_symm_trans (nv,nv2,ia,ja,a,x,y)
  implicit none
  integer(ip), intent(in)  :: nv,nv2,ia(nv+1),ja(ia(nv+1)-1)
  real(rp)   , intent(in)  :: a(ia(nv+1)-1), x(nv)
  real(rp)   , intent(out) :: y(nv2)
  integer(ip)              :: iv,iz,jv,id,jd
  integer(ip)              :: of, ivc, izc, jvc

  write (0,*) 'Error: the body of matvec_csr_symm_trans in matvec_dof.i90 still to be written'
  write (0,*) 'Error: volunteers are welcome !!!'
  stop 
end subroutine matvec_csr_symm_trans


! Not tested
subroutine matvec_csc (nv,nv2,ia,ja,a,x,y)
    implicit none
    integer(ip), intent(in)  :: nv,nv2,ia(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: a(ia(nv+1)-1), x(nv2)
    real(rp)   , intent(out) :: y(nv)
    integer(ip)              :: iv,iz,jv,id,jd
    integer(ip)              :: of, ivc, izc, jvc

    y = 0.0_rp
    do iv = 1, nv2, ndof2
       do iz = ia(iv), ia(iv+1)-1, ndof1
          of = iz - ia(iv)
          jv = ja(iz)
          do ivc = iv, iv + ndof2 - 1
             jvc   = jv
             do izc = ia(ivc)+of, ia(ivc)+of+ndof1-1
                y(jvc) = y(jvc) + x(ivc)*a(izc)
                jvc      = jvc + 1
             end do ! izc
          end do ! ivc
       end do ! iz
    end do ! iv
end subroutine matvec_csc


end module matvec_dof




       
