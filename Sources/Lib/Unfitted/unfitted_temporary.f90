
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

module unfitted_temporary_names
  use types_names
  use memor_names
  use reference_fe_names

  implicit none
# include "debug.i90"

contains

!========================================================================================
! TODO @fverdugo FEMPAR PRIORITY MEDIUM EFFORT HIGH
! A better way to do this?
! This subroutine assumes ref elem with hex topology and linear order
! Where to put this info? At the reference element?
subroutine evaluate_monomials(points,monomials,degree)
  implicit none
  type(quadrature_t),    intent(in)    :: points
  real(rp), allocatable, intent(inout) :: monomials(:,:)
  integer(ip),           intent(in)    :: degree

  integer(ip) :: q_point
  real(rp), pointer :: quad_coords(:,:)
  integer(ip) :: imo, px, py, pz

  assert(allocated(monomials))
  assert(size(monomials,1)==points%get_num_quadrature_points())

  quad_coords => points%get_coordinates()
  select case(points%get_num_dims())
    case(1)
      assert(size(monomials,2)==(degree+1))
      do q_point = 1, points%get_num_quadrature_points()
        imo = 1
        do px = 0, degree
          monomials(q_point,imo) = quad_coords(1,q_point)**px
          imo = imo + 1
        end do
      end do
    case(2)
      assert( size(monomials,2)== ((degree+1)**2) )
      do q_point = 1, points%get_num_quadrature_points()
        imo = 1
        do px = 0, degree
          do py = 0, degree
            monomials(q_point,imo) = (quad_coords(1,q_point)**px)&
                                   * (quad_coords(2,q_point)**py)
            imo = imo + 1
          end do
        end do
      end do
    case(3)
      assert(size(monomials,2)==((degree+1)**3))
      do q_point = 1, points%get_num_quadrature_points()
        imo = 1
        do px = 0, degree
          do py = 0, degree
            do pz = 0, degree
              monomials(q_point,imo) = (quad_coords(1,q_point)**px)&
                                     * (quad_coords(2,q_point)**py)&
                                     * (quad_coords(3,q_point)**pz)
              imo = imo + 1
            end do
          end do
        end do
      end do
    case default
      check(.false.)
  end select

end subroutine evaluate_monomials

!========================================================================================
  subroutine At_times_B_times_A(A,B,C)
  !TODO @fverdugo FEMPAR PRIORITY HIGH EFFORT LOW
  ! Is there a way in FEMPAR to multiply two dense matrices?
  ! This can be optimized
    implicit none
    real(rp), intent(in)    :: A(:,:)
    real(rp), intent(in)    :: B(:,:)
    real(rp), intent(inout) :: C(:,:)
    integer(ip) :: i,j,m,n
    assert( lbound(C,1)==lbound(A,2) .and. ubound(C,1)==ubound(A,2) )
    assert( lbound(A,1)==lbound(B,1) .and. ubound(A,1)==ubound(B,1) )
    assert( lbound(B,2)==lbound(A,1) .and. ubound(B,2)==ubound(A,1) )
    assert( lbound(C,2)==lbound(A,2) .and. ubound(C,2)==ubound(A,2) )
    do i = lbound(C,1), ubound(C,1)
      do j = lbound(C,2), ubound(C,2)
        C(i,j) = 0
        do m = lbound(A,1), ubound(A,1)
          do n = lbound(A,1), ubound(A,1)
            C(i,j) = C(i,j) + A(m,i)*B(m,n)*A(n,j)
          end do
        end do
      end do
    end do
  end subroutine At_times_B_times_A

  subroutine find_facets_neighbors_in_mesh(T,nbel,facets,nbef,neigs)
    implicit none
    integer(ip), intent(in)    :: T(:,:)
    integer(ip), intent(in)    :: nbel
    integer(ip), intent(in)    :: facets(:,:)
    integer(ip), intent(in)    :: nbef
    integer(ip), intent(inout) :: neigs(:,:)

    integer(ip) :: iel
    integer(ip) :: ifa

    assert(size(T,2)>=nbel)
    assert(size(facets,2)==nbef)
    assert(size(neigs,1)==nbef)
    assert(size(neigs,2)>=nbel)

    do iel = 1, nbel
      do ifa = 1, nbef
      neigs(ifa,iel) = find_neigbour(T,iel,nbel,facets,ifa)
      end do
    end do

    contains

      function find_neigbour(T,iel,nbel,facets,ifa)
        implicit none
        integer(ip), intent(in)    :: T(:,:)
        integer(ip), intent(in)    :: iel
        integer(ip), intent(in)    :: nbel
        integer(ip), intent(in)    :: facets(:,:)
        integer(ip), intent(in)    :: ifa
        integer(ip) :: find_neigbour

        integer(ip) :: jel
        integer(ip) :: jen
        integer(ip) :: ifn
        logical :: has_all_facets_nodes
        logical :: node_found

        find_neigbour = 0
        do jel = 1, nbel
          if (iel /= jel) then

            has_all_facets_nodes = .true.
            do ifn = 1, size(facets,1)

              node_found = .false.
              do jen = 1,size(T,1)
                if ( T(jen,jel) == T(facets(ifn,ifa),iel) ) then
                  node_found = .true.
                  exit
                end if
              end do

              has_all_facets_nodes = has_all_facets_nodes .and. node_found
            end do

            if (has_all_facets_nodes) then
              assert(find_neigbour == 0)
              find_neigbour = jel
            end if

          end if
        end do
      end function find_neigbour
  end subroutine find_facets_neighbors_in_mesh

end module unfitted_temporary_names
