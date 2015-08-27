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
module statistics_names
  use memor_names
  use types_names
  use generate_uniform_triangulation_names
  use fe_space_names
  use interpolation_tools_names
  implicit none
# include "debug.i90"
  private

  type line_statistics_t
     real(rp), allocatable :: a(:,:)        ! Values
     integer(ip)           :: lunio         ! Output unit
     integer(ip)           :: direction     ! Direction in which we want to store statistics
     integer(ip)           :: order         ! Order of interpolation of the field
     integer(ip)           :: npoints       ! Number of points where values are stored
     integer(ip)           :: ncomp         ! Field components
     integer(ip)           :: iblock        ! Field block
     integer(ip)           :: iprob         ! Problem ID
     type(uniform_mesh_descriptor_t), pointer :: gdata => NULL()    ! Pointer to uniform mesh descriptor
     type(fe_space_t)               , pointer :: fe_space => NULL() ! Pointer to FE space
   contains
     procedure :: initialize => initialize_line_statistics
     procedure :: compute    => compute_line_statistics
     procedure :: finalize   => finalize_line_statistics
  end type line_statistics_t

  ! line statistics:
  ! %a(1:ncomp,:)                        -->  sum(field(i))         , for i=1,...,ncomp
  ! %a(ncomp+1:ncomp+(ncomp+1)*ncomp/2)  -->  sum(field(i)*field(j)), for i=1,...,ncomp, j=i,...,ncomp

  ! Types
  public :: line_statistics_t
  
  ! Functions
  public :: globalid_to_ijk

contains

  !==================================================================================================
  subroutine initialize_line_statistics(statistics,gdata,fe_space,direction,order,lunio,iblock,iprob)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine initialize the variables needed to compute the statistics of a field        !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(line_statistics_t)              , intent(inout) :: statistics
    type(uniform_mesh_descriptor_t),target, intent(in)    :: gdata
    type(fe_space_t)               ,target, intent(in)    :: fe_space
    integer(ip)                           , intent(in)    :: direction,order,lunio
    integer(ip), optional                 , intent(in)    :: iblock,iprob
    ! Locals
    integer(ip) :: aux

    ! Fill statistics
    statistics%lunio = lunio
    statistics%direction = direction
    statistics%order = order
    statistics%iblock = 1
    statistics%iprob  = 1
    if(present(iblock)) statistics%iblock = iblock
    if(present(iblock)) statistics%iprob  = iprob
    
    ! Compute npoints
    statistics%npoints = gdata%nedir(direction)*order + 1
    statistics%ncomp   = fe_space%dof_descriptor%prob_block(statistics%iblock,statistics%iprob)%nd1
    
    ! Allocate values array
    aux = statistics%ncomp + statistics%ncomp*(statistics%ncomp+1)/2
    call memalloc(aux,statistics%npoints,statistics%a,__FILE__,__LINE__)
    statistics%a = 0.0_rp

    ! Pointers
    statistics%gdata => gdata
    statistics%fe_space => fe_space

  end subroutine initialize_line_statistics

  !==================================================================================================
  subroutine compute_line_statistics(statistics)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine compute the statistics of a field                                           !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(line_statistics_t)        , intent(inout) :: statistics
    ! Locals
    integer(ip) :: ndime,pdegr,nvapb,count
    integer(ip) :: ijkaux(3),npoin(3),ijk(3),ijkelem(3),ijkpoin(3),ijknode(3)
    integer(ip) :: idime,ielem,iprob,ivar,lvar,inode,jnode,i,j,k,icomp,jcomp
    integer(ip), allocatable :: permu(:)
    real(rp), allocatable :: field(:,:,:,:)

    ! Checks
    check(associated(statistics%gdata))
    check(associated(statistics%fe_space))

    ! Auxiliar variables
    ndime  = statistics%gdata%ndime
    pdegr  = statistics%order
    npoin(1) = statistics%gdata%nedir(1)*pdegr+1
    npoin(2) = statistics%gdata%nedir(2)*pdegr+1
    npoin(3) = statistics%gdata%nedir(3)*pdegr+1
    call memalloc(3**ndime,permu,__FILE__,__LINE__)
    if(ndime==2) then
       permu = (/1,4,7,2,5,8,3,6,9/)
    elseif(ndime==3) then
       permu = (/1,10,19,4,13,22,7,16,25,2,11,20,5,14,23,8,17,26,3,12,21,6,15,24,9,18,27/)
    end if

    ! Auxiliar ijk array
    ijkaux = 1
    do idime=1,ndime
       ijkaux(idime) = 3
    end do

    ! Allocate ijk field array
    call memalloc(npoin(1),npoin(2),npoin(3),statistics%ncomp,field,__FILE__,__LINE__)
    field = 0.0_rp

    ! Loop over elements
    do ielem = 1,statistics%fe_space%g_trian%num_elems
       iprob = statistics%fe_space%finite_elements(ielem)%problem
       nvapb = statistics%fe_space%dof_descriptor%prob_block(statistics%iblock,iprob)%nd1

       if(iprob/=statistics%iprob) cycle

       ! Global ID to ijk
       call globalid_to_ijk(ndime,statistics%gdata%nedir,ielem,ijkelem)

       ! Loop over problem and block variables
       do ivar = 1, nvapb
          lvar = statistics%fe_space%dof_descriptor%prob_block(statistics%iblock,iprob)%a(ivar)

          ! Loop over elemental nodes
          do inode = 1,statistics%fe_space%finite_elements(ielem)%reference_element_vars(lvar)%p%nnode
             jnode = permu(inode)

             ! Local ID to ijk
             call globalid_to_ijk(ndime,ijkaux,jnode,ijknode)

             ! Global point ID
             ijkpoin = (ijkelem - 1)*pdegr + ijknode

             ! Fill field
             field(ijkpoin(1),ijkpoin(2),ijkpoin(3),lvar) = statistics%fe_space%finite_elements(ielem)%unkno(inode,lvar,1)

          end do
       end do
    end do

    ! Compute statistics
    do i=1,npoin(1)
       do j=1,npoin(2)
          do k=1,npoin(3)
             ijk = (/i,j,k/)
             count = 0
             do icomp=1,statistics%ncomp
                statistics%a(icomp,ijk(statistics%direction)) = statistics%a(icomp,ijk(statistics%direction)) + &   ! U_i
                     &                                          field(i,j,k,icomp)
                do jcomp=icomp,statistics%ncomp
                   count = count + 1
                   statistics%a(statistics%ncomp+count,ijk(statistics%direction)) =        &                        ! U_i*U_j
                        & statistics%a(statistics%ncomp+count,ijk(statistics%direction)) + &
                        & field(i,j,k,icomp)*field(i,j,k,jcomp)
                end do
             end do
          end do
       end do
    end do
    do i=1,3
       if(i==statistics%direction) cycle
       statistics%a = statistics%a/npoin(i)
    end do

    ! Deallocate
    call memfree(permu,__FILE__,__LINE__)
    call memfree(field,__FILE__,__LINE__)
    
  end subroutine compute_line_statistics

  !==================================================================================================
  subroutine finalize_line_statistics(statistics)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine prints the statistics of a field and destroy its variables.                 !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(line_statistics_t)        , intent(inout) :: statistics    
    ! Locals
    integer(ip) :: ndime,pdegr,nvapb
    integer(ip) :: ijkaux(3),npoin(3),ijk(3),ijkelem(3),ijkpoin(3),ijknode(3)
    integer(ip) :: idime,ielem,iprob,ivar,lvar,inode,jnode,i,j,k,icomp,jcomp
    integer(ip), allocatable :: permu(:)
    real(rp), allocatable :: ycoord_ijk(:,:,:),elcoord(:,:),ycoord(:)

    ! Checks
    check(associated(statistics%gdata))
    check(associated(statistics%fe_space))

    ! Auxiliar variables
    ndime  = statistics%gdata%ndime
    pdegr  = statistics%order
    npoin(1) = statistics%gdata%nedir(1)*pdegr+1
    npoin(2) = statistics%gdata%nedir(2)*pdegr+1
    npoin(3) = statistics%gdata%nedir(3)*pdegr+1
    call memalloc(3**ndime,permu,__FILE__,__LINE__)
    if(ndime==2) then
       permu = (/1,4,7,2,5,8,3,6,9/)
    elseif(ndime==3) then
       permu = (/1,10,19,4,13,22,7,16,25,2,11,20,5,14,23,8,17,26,3,12,21,6,15,24,9,18,27/)
    end if

    ! Auxiliar ijk array
    ijkaux = 1
    do idime=1,ndime
       ijkaux(idime) = 3
    end do

    ! Allocate ijk coordinates array
    call memalloc(npoin(1),npoin(2),npoin(3),ycoord_ijk,__FILE__,__LINE__)
    call memalloc(npoin(statistics%direction),ycoord,__FILE__,__LINE__)
    ycoord_ijk = 0.0_rp

    ! Loop over elements
    do ielem = 1,statistics%fe_space%g_trian%num_elems
       iprob = statistics%fe_space%finite_elements(ielem)%problem
       nvapb = statistics%fe_space%dof_descriptor%prob_block(statistics%iblock,iprob)%nd1

       if(iprob/=statistics%iprob) cycle

       ! Global ID to ijk
       call globalid_to_ijk(ndime,statistics%gdata%nedir,ielem,ijkelem)

       ! Direction variable
       ivar = statistics%direction
       lvar = statistics%fe_space%dof_descriptor%prob_block(statistics%iblock,iprob)%a(ivar)

       ! Get coordinates
       call memalloc(ndime,statistics%fe_space%finite_elements(ielem)%reference_element_vars(lvar)%p%nnode, &
            &        elcoord,__FILE__,__LINE__)
       call interpolate(ndime, statistics%fe_space%finite_elements(ielem)%p_geo_reference_element%nnode, &
            &           statistics%fe_space%finite_elements(ielem)%reference_element_vars(lvar)%p%nnode,       &
            &           statistics%fe_space%finite_elements(ielem)%inter(lvar)%p,                              &
            &           statistics%fe_space%g_trian%elems(ielem)%coordinates,elcoord)

       ! Loop over elemental nodes
       do inode = 1,statistics%fe_space%finite_elements(ielem)%reference_element_vars(lvar)%p%nnode
          jnode = permu(inode)

          ! Local ID to ijk
          call globalid_to_ijk(ndime,ijkaux,jnode,ijknode)

          ! Global point ID
          ijkpoin = (ijkelem - 1)*pdegr + ijknode

          ! Fill coord
          ycoord_ijk(ijkpoin(1),ijkpoin(2),ijkpoin(3)) = elcoord(statistics%direction,inode)

       end do

       ! Deallocate elcoord
       call memfree(elcoord,__FILE__,__LINE__)
    end do

    ! Select coordinates
    if(statistics%direction==1) then
       ycoord = ycoord_ijk(:,1,1)
    elseif(statistics%direction==2) then
       ycoord = ycoord_ijk(1,:,1)
    elseif(statistics%direction==3) then
       ycoord = ycoord_ijk(1,1,:)
    end if

    ! Write file
    if(ndime==2)then
       do j=1,npoin(statistics%direction)
          write(statistics%lunio,1) j,ycoord(j),statistics%a(:,j)
       end do
    elseif(ndime==3) then
       do j=1,npoin(statistics%direction)
          write(statistics%lunio,2) j,ycoord(j),statistics%a(:,j)
       end do
    end if

    ! Deallocate
    call memfree(permu,__FILE__,__LINE__)
    call memfree(ycoord,__FILE__,__LINE__)
    call memfree(ycoord_ijk,__FILE__,__LINE__)
    call memfree(statistics%a,__FILE__,__LINE__)

    ! Nullify pointers
    statistics%gdata => null()
    statistics%fe_space => null()

1   format(i10,6(1x,e16.8e3))
2   format(i10,10(1x,e16.8e3))
    
  end subroutine finalize_line_statistics

  !==================================================================================================
  subroutine globalid_to_ijk(ndime,nd,gl,ijk)
    implicit none
    integer(ip), intent(in)  :: gl,nd(3),ndime
    integer(ip), intent(out) :: ijk(3)
    ! Locals
    integer(ip) :: aux

    if(ndime==1) then
       ijk(3) = 1
       ijk(2) = 1
       ijk(1) = gl
    elseif(ndime==2) then
       ijk(1) = floor(real(gl-1)/nd(2)) + 1
       ijk(2) = gl - (ijk(1)-1)*nd(2)
       ijk(3) = 1
    else if(ndime==3) then
       ijk(1) = floor(real(gl-1)/(nd(2)*nd(3)))+1
       aux = gl - (ijk(1)-1)*(nd(2)*nd(3))
       ijk(2) = floor(real(aux-1)/nd(3))+1
       ijk(3) = aux - (ijk(2)-1)* nd(3)
    end if

  end subroutine globalid_to_ijk
    

end module statistics_names
