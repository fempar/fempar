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
module quadrature_faces_names
  use types_names
  use memor_names
  use fe_space_types_names
  use quadrature_names

# include "debug.i90"
  implicit none
  private

  type face_quadrature_t
     integer(ip)              :: &
        ndime,                   &    ! Number of space dimensions
        etype,                   &    ! Element type ID
        nface,                   &    ! Number of faces of the element
        ngausxface                    ! Number of integration points

     real(rp)   , allocatable :: &
        pos(:,:,:)                      ! Quadrature points position
  end type face_quadrature_t

  ! Types
  public :: face_quadrature_t

  ! Functions
  public :: face_quadrature_create, face_quadrature_free

contains

  !==================================================================================================
  subroutine face_quadrature_create(ndime,nface,quad,fi,face_quad,subface)
    implicit none
    ! Parameters
    integer(ip)          , intent(in)  :: ndime,nface
    type(fem_fixed_info_t) , intent(in)  :: fi
    type(quadrature_t)     , intent(in)  :: quad
    type(face_quadrature_t), intent(out) :: face_quad
    integer(ip)          , intent(in), optional :: subface

    ! Local variables
    integer(ip)                   :: iface,igaus,ngaus,i
    integer(ip)                   :: subface_pos(ndime),idime, refinement_level
    real(rp)                      :: center, length
    real(rp)                      :: nopos(ndime,fi%nnode)
    real(rp)                      :: v1(ndime),M(ndime,ndime-1),posgp(ndime-1,quad%ngaus)

    ! Assign integer parameters
    face_quad%ndime      = ndime
    face_quad%etype      = fi%ftype
    face_quad%nface      = nface
    face_quad%ngausxface = quad%ngaus

    ngaus                = quad%ngaus
    posgp                = quad%pos
    
    if (present(subface)) then
       if (subface > 0) then
          ! TODO: Use refinement level as a variable. For now it is assumed to be at most 1.
          refinement_level = 1
          if (fi%ftype == Q_type_id) then
             call Q_ijkg(subface_pos,subface,ndime,2**refinement_level-1)
             length = 2.0_rp/real(2**refinement_level)
             do idime =1, ndime-1
                center = -1.0_rp + length*(real(subface_pos(idime))+0.5_rp)
                do igaus = 1,ngaus
                   posgp(idime,igaus) = center + length/2.0*posgp(idime,igaus)
                end do
             end do
          else
             write(*,*) __FILE__,__LINE__, 'ERROR:: This element type subface had still not been implemented'
          end if
       end if
    end if
    call memalloc(ndime,quad%ngaus,nface,face_quad%pos,__FILE__,__LINE__)

    if (fi%ftype == P_type_id) then
       call P_refcoord (nopos, ndime, fi%order, fi%nnode)
    elseif (fi%ftype == Q_type_id) then
       call Q_refcoord (nopos, ndime, fi%order, fi%nnode)
    else
       write(*,*) 'quadrature_t_faces:: ERROR! Unknown element type.'
       check (1 == 0)
    end if

    if (ndime == 1) then
       face_quad%pos(1,:,1) = (/-1,1/)
    else if (ndime == 2) then
       do iface = 1,nface
          v1 = nopos(:,fi%crxob%l(fi%crxob%p(fi%nobje_dim(ndime)+iface-1)))
          M(:,1) = 0.5*(nopos(:,fi%crxob%l(fi%crxob%p(fi%nobje_dim(ndime)+iface-1)+1)) -         &
               &        nopos(:,fi%crxob%l(fi%crxob%p(fi%nobje_dim(ndime)+iface-1)  )))
          do igaus = 1,ngaus
             face_quad%pos(:,igaus,iface) = v1 + M(:,1)*(posgp(1,igaus)+1)
          end do
       end do
    else if (ndime == 3) then
       if (fi%ftype == P_type_id) then
          do iface = 1,nface
             v1 = nopos(:,fi%crxob%l(fi%crxob%p(fi%nobje_dim(ndime)+iface-1)))
             M(:,1) = nopos(:,fi%crxob%l(fi%crxob%p(fi%nobje_dim(ndime)+iface-1)+1)) -              &
                  &   nopos(:,fi%crxob%l(fi%crxob%p(fi%nobje_dim(ndime)+iface-1)))
             M(:,2) = nopos(:,fi%crxob%l(fi%crxob%p(fi%nobje_dim(ndime)+iface-1)+2)) -              &
                  &   nopos(:,fi%crxob%l(fi%crxob%p(fi%nobje_dim(ndime)+iface-1)))
             do igaus = 1,ngaus
                face_quad%pos(:,igaus,iface) = v1 + matmul(M,posgp(:,igaus))
             end do
          end do
       elseif (fi%ftype == Q_type_id) then
          do iface = 1,nface
             v1 = 0.0_rp
             do i=1,4
                v1 = v1 + nopos(:,fi%crxob%l(fi%crxob%p(fi%nobje_dim(ndime)+iface-1)+i-1))
             end do
             v1 = v1/4.0_rp

             M(:,1) = 0.5*(nopos(:,fi%crxob%l(fi%crxob%p(fi%nobje_dim(ndime)+iface-1)+1)) -         &
                  &        nopos(:,fi%crxob%l(fi%crxob%p(fi%nobje_dim(ndime)+iface-1))))
             M(:,2) = 0.5*(nopos(:,fi%crxob%l(fi%crxob%p(fi%nobje_dim(ndime)+iface-1)+2)) -         &
                  &        nopos(:,fi%crxob%l(fi%crxob%p(fi%nobje_dim(ndime)+iface-1))))
             do igaus = 1,ngaus
                face_quad%pos(:,igaus,iface) = v1 + matmul(M,posgp(:,igaus))
             end do
          end do
       end if
    else
       write(*,*) __FILE__,__LINE__, 'Face quadrature not implemented for ndime>3'
    end if

  end subroutine face_quadrature_create

  !==================================================================================================
  subroutine face_quadrature_free(rule)
    !-----------------------------------------------------------------------
    ! 
    ! This routine frees the integration rule
    !
    !-----------------------------------------------------------------------
    implicit none
    type(face_quadrature_t), intent(inout) :: rule

    call memfree(rule%pos  ,__FILE__,__LINE__)
  end subroutine face_quadrature_free

end module quadrature_faces_names
