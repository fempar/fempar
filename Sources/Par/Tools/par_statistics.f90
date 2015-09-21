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
module par_statistics_names
  use memor_names
  use types_names
  use statistics_names
  use generate_uniform_triangulation_names
  use par_fe_space_names
  use par_environment_names
  use blocks_dof_distribution_names
  use interpolation_tools_names
  use psi_reduce_mod_names
  implicit none
# include "debug.i90"
  private

  type par_line_statistics_t
     type(line_statistics_t)          :: f_stats
     type(par_fe_space_t)   , pointer :: p_fe_space => NULL()
     type(par_environment_t), pointer :: p_env => NULL()
     integer(ip), allocatable         :: nparts_per_dof(:)
   contains
     procedure :: initialize => initialize_par_line_statistics
     procedure :: compute    => compute_par_line_statistics
     procedure :: finalize   => finalize_par_line_statistics
  end type par_line_statistics_t

  ! Types
  public :: par_line_statistics_t
  
contains

  !==================================================================================================
  subroutine initialize_par_line_statistics(p_stats,gdata,p_fe_space,p_env,blocks_dof_distribution,direction, &
       &                                    order,lunio,iblock,iprob)
    implicit none
    class(par_line_statistics_t)    , intent(inout) :: p_stats
    type(uniform_mesh_descriptor_t), intent(in)    :: gdata
    type(par_fe_space_t)    ,target, intent(in)    :: p_fe_space
    type(par_environment_t) ,target, intent(in)    :: p_env
    type(blocks_dof_distribution_t) , intent(in)   :: blocks_dof_distribution
    integer(ip)                    , intent(in)    :: direction,order,lunio
    integer(ip), optional          , intent(in)    :: iblock,iprob
    ! Locals
    integer(ip) :: iobje,iniobj,endobj

    ! Fill parallel types
    p_stats%p_fe_space => p_fe_space
    p_stats%p_env      => p_env
    
    if(p_env%am_i_fine_task()) then

       ! Fill serial type
       call p_stats%f_stats%initialize(gdata,p_fe_space%serial_fe_space,direction,order,lunio,iblock,iprob)

       ! Construct npats_per_dof
       call memalloc(blocks_dof_distribution%blocks(p_stats%f_stats%iblock)%nl,p_stats%nparts_per_dof,__FILE__,__LINE__)
       p_stats%nparts_per_dof = 0
       do iobje=1,blocks_dof_distribution%blocks(p_stats%f_stats%iblock)%nobjs
          iniobj = blocks_dof_distribution%blocks(p_stats%f_stats%iblock)%lobjs(2,iobje)
          endobj = blocks_dof_distribution%blocks(p_stats%f_stats%iblock)%lobjs(3,iobje)
          p_stats%nparts_per_dof(iniobj:endobj) = blocks_dof_distribution%blocks(p_stats%f_stats%iblock)%lobjs(4,iobje)
       end do

    end if
    
  end subroutine initialize_par_line_statistics

  !==================================================================================================
  subroutine compute_par_line_statistics(p_stats)
    implicit none
    class(par_line_statistics_t), intent(inout) ::p_stats
    ! Locals
    integer(ip) :: ndime,pdegr,nvapb,count
    integer(ip) :: ijkaux(3),npoin(3),ijk(3),ijkelem(3),ijkpoin(3),ijknode(3)
    integer(ip) :: idime,ielem,iprob,ivar,lvar,inode,jnode,i,j,k,icomp,jcomp,jelem,idof
    integer(ip), allocatable :: permu(:)
    real(rp), allocatable :: field(:,:,:,:)

    ! Checks
    check(associated(p_stats%p_fe_space))
    check(associated(p_stats%p_env))
    
    if(p_stats%p_env%am_i_fine_task()) then

       ! Checks
       check(associated(p_stats%f_stats%gdata))
       check(associated(p_stats%f_stats%fe_space))

       ! Auxiliar variables
       ndime  = p_stats%f_stats%gdata%ndime
       pdegr  = p_stats%f_stats%order
       npoin(1) = p_stats%f_stats%gdata%nedir(1)*pdegr+1
       npoin(2) = p_stats%f_stats%gdata%nedir(2)*pdegr+1
       npoin(3) = p_stats%f_stats%gdata%nedir(3)*pdegr+1
       call memalloc(3**ndime,permu,__FILE__,__LINE__)
       if(ndime==2) then
          permu = (/1,4,7,2,5,8,3,6,9/)
       elseif(ndime==3) then
          permu = (/1,10,19,4,13,22,7,16,25,2,11,20,5,14,23,8,17,26,3,12,21,6,15,24,9,18,27/)
       end if

       ! Construct auxiliar array
       ijkaux = 1
       do idime=1,ndime
          ijkaux(idime) = 3
       end do

       ! Allocate ijk field array
       call memalloc(npoin(1),npoin(2),npoin(3),p_stats%f_stats%ncomp,field,__FILE__,__LINE__)
       field = 0.0_rp

       ! Loop over elements
       do ielem = 1, p_stats%p_fe_space%serial_fe_space%g_trian%num_elems
          jelem = p_stats%p_fe_space%p_trian%elems(ielem)%globalID
          iprob = p_stats%p_fe_space%serial_fe_space%finite_elements(ielem)%problem
          nvapb = p_stats%p_fe_space%serial_fe_space%dof_descriptor%prob_block(p_stats%f_stats%iblock,iprob)%nd1

          if(iprob/=p_stats%f_stats%iprob) cycle

          ! Global ID to ijk
          call globalid_to_ijk(ndime,p_stats%f_stats%gdata%nedir,jelem,ijkelem)

          ! Loop over problem and block variables
          do ivar = 1, nvapb
             lvar = p_stats%p_fe_space%serial_fe_space%dof_descriptor%prob_block(p_stats%f_stats%iblock,iprob)%a(ivar)

             ! Loop over elemental nodes
             do inode = 1,p_stats%p_fe_space%serial_fe_space%finite_elements(ielem)%reference_element_vars(lvar)%p%nnode
                idof  =  p_stats%p_fe_space%serial_fe_space%finite_elements(ielem)%elem2dof(inode,lvar)
                jnode = permu(inode)

                ! Local ID to ijk
                call globalid_to_ijk(ndime,ijkaux,jnode,ijknode)

                ! Global point ID
                ijkpoin = (ijkelem - 1)*pdegr + ijknode

                ! Fill field
                if(idof>0) then
                   field(ijkpoin(1),ijkpoin(2),ijkpoin(3),lvar) = &
                        & p_stats%p_fe_space%serial_fe_space%finite_elements(ielem)%unkno(inode,lvar,1)/p_stats%nparts_per_dof(idof)
                else
                   field(ijkpoin(1),ijkpoin(2),ijkpoin(3),lvar) = &
                        & p_stats%p_fe_space%serial_fe_space%finite_elements(ielem)%unkno(inode,lvar,1)
                end if

             end do
          end do
       end do

       ! Compute statistics
       do i=1,npoin(1)
          do j=1,npoin(2)
             do k=1,npoin(3)
                ijk = (/i,j,k/)
                count = 0
                do icomp=1,p_stats%f_stats%ncomp
                   p_stats%f_stats%a(icomp,ijk(p_stats%f_stats%direction)) = &
                        &  p_stats%f_stats%a(icomp,ijk(p_stats%f_stats%direction)) + &                           ! U_i
                        &  field(i,j,k,icomp)
                   do jcomp=icomp,p_stats%f_stats%ncomp
                      count = count + 1
                      p_stats%f_stats%a(p_stats%f_stats%ncomp+count,ijk(p_stats%f_stats%direction)) =        &   ! U_i*U_j
                           & p_stats%f_stats%a(p_stats%f_stats%ncomp+count,ijk(p_stats%f_stats%direction)) + &
                           & field(i,j,k,icomp)*field(i,j,k,jcomp)
                   end do
                end do
             end do
          end do
       end do
       do i=1,3
          if(i==p_stats%f_stats%direction) cycle
          p_stats%f_stats%a = p_stats%f_stats%a/npoin(i)
       end do
       
       ! Deallocate
       call memfree(permu,__FILE__,__LINE__)
       call memfree(field,__FILE__,__LINE__)

    end if
    
  end subroutine compute_par_line_statistics

  !==================================================================================================
  subroutine finalize_par_line_statistics(p_stats)
    implicit none   
    class(par_line_statistics_t), intent(inout) ::p_stats
    ! Locals
    integer(ip) :: ndime,pdegr,nvapb,me,np
    integer(ip) :: ijkaux(3),npoin(3),ijk(3),ijkelem(3),ijkpoin(3),ijknode(3)
    integer(ip) :: idime,ielem,iprob,ivar,lvar,inode,jnode,i,j,k,icomp,jcomp,jelem,idof
    integer(ip), allocatable :: permu(:)
    real(rp), allocatable :: ycoord_ijk(:,:,:),elcoord(:,:),ycoord(:)

    ! Checks
    check(associated(p_stats%p_fe_space))
    check(associated(p_stats%p_env))
    
    if(p_stats%p_env%am_i_fine_task()) then

       ! Checks
       check(associated(p_stats%f_stats%gdata))
       check(associated(p_stats%f_stats%fe_space))

       ! Auxiliar variables
       ndime  = p_stats%f_stats%gdata%ndime
       pdegr  = p_stats%f_stats%order
       npoin(1) = p_stats%f_stats%gdata%nedir(1)*pdegr+1
       npoin(2) = p_stats%f_stats%gdata%nedir(2)*pdegr+1
       npoin(3) = p_stats%f_stats%gdata%nedir(3)*pdegr+1
       call memalloc(3**ndime,permu,__FILE__,__LINE__)
       if(ndime==2) then
          permu = (/1,4,7,2,5,8,3,6,9/)
       elseif(ndime==3) then
          permu = (/1,10,19,4,13,22,7,16,25,2,11,20,5,14,23,8,17,26,3,12,21,6,15,24,9,18,27/)
       end if

       ! Construct auxiliar array
       ijkaux = 1
       do idime=1,ndime
          ijkaux(idime) = 3
       end do

       ! Allocate ijk array
       call memalloc(npoin(1),npoin(2),npoin(3),ycoord_ijk,__FILE__,__LINE__)
       call memalloc(npoin(p_stats%f_stats%direction),ycoord,__FILE__,__LINE__)
       ycoord_ijk = 0.0_rp

       ! Loop over elements
       do ielem = 1, p_stats%p_fe_space%serial_fe_space%g_trian%num_elems
          jelem = p_stats%p_fe_space%p_trian%elems(ielem)%globalID
          iprob = p_stats%p_fe_space%serial_fe_space%finite_elements(ielem)%problem
          nvapb = p_stats%p_fe_space%serial_fe_space%dof_descriptor%prob_block(1,iprob)%nd1

          if(iprob/=p_stats%f_stats%iprob) cycle

          ! Global ID to ijk
          call globalid_to_ijk(ndime,p_stats%f_stats%gdata%nedir,jelem,ijkelem)

          ! Direction variable
          ivar = p_stats%f_stats%direction
          lvar = p_stats%p_fe_space%serial_fe_space%dof_descriptor%prob_block(p_stats%f_stats%iblock,iprob)%a(ivar)

          ! Get coordinates
          call memalloc(ndime,p_stats%p_fe_space%serial_fe_space%finite_elements(ielem)%reference_element_vars(lvar)%p%nnode, &
               &        elcoord,__FILE__,__LINE__)
          call interpolate(ndime, p_stats%p_fe_space%serial_fe_space%finite_elements(ielem)%p_geo_reference_element%nnode, &
               &           p_stats%p_fe_space%serial_fe_space%finite_elements(ielem)%reference_element_vars(lvar)%p%nnode, &
               &           p_stats%p_fe_space%serial_fe_space%finite_elements(ielem)%inter(lvar)%p,                        &
               &           p_stats%p_fe_space%p_trian%f_trian%elems(ielem)%coordinates,elcoord)

          ! Loop over elemental nodes
          do inode = 1,p_stats%p_fe_space%serial_fe_space%finite_elements(ielem)%reference_element_vars(lvar)%p%nnode
             idof  = p_stats%p_fe_space%serial_fe_space%finite_elements(ielem)%elem2dof(inode,lvar)
             jnode = permu(inode)

             ! Local ID to ijk
             call globalid_to_ijk(ndime,ijkaux,jnode,ijknode)

             ! Global point ID
             ijkpoin = (ijkelem - 1)*pdegr + ijknode

             ! Fill veloc
             ycoord_ijk(ijkpoin(1),ijkpoin(2),ijkpoin(3)) = elcoord(p_stats%f_stats%direction,inode)

          end do

          ! Deallocate elcoord
          call memfree(elcoord,__FILE__,__LINE__)
       end do

       ! Select coordinates
       if(p_stats%f_stats%direction==1) then
          ycoord = ycoord_ijk(:,1,1)
       elseif(p_stats%f_stats%direction==2) then
          ycoord = ycoord_ijk(1,:,1)
       elseif(p_stats%f_stats%direction==3) then
          ycoord = ycoord_ijk(1,1,:)
       end if
       
       ! Reduce stats to root
       call psb_sum(p_stats%p_env%p_context%icontxt,p_stats%f_stats%a)
       call psb_amx(p_stats%p_env%p_context%icontxt,ycoord)

       ! Get process info
       call p_stats%p_env%info(me,np)

       ! Write file
       if(me==0) then
          if(ndime==2)then
             do j=1,npoin(p_stats%f_stats%direction)
                write(p_stats%f_stats%lunio,1) j,ycoord(j),p_stats%f_stats%a(:,j)
             end do
          elseif(ndime==3) then
             do j=1,npoin(p_stats%f_stats%direction)
                write(p_stats%f_stats%lunio,2) j,ycoord(j),p_stats%f_stats%a(:,j)
             end do
          end if
       end if

       ! Deallocate
       call memfree(permu,__FILE__,__LINE__)
       call memfree(ycoord,__FILE__,__LINE__)
       call memfree(ycoord_ijk,__FILE__,__LINE__)
       call memfree(p_stats%nparts_per_dof,__FILE__,__LINE__)
       call memfree(p_stats%f_stats%a,__FILE__,__LINE__)

       ! Nullify pointers
       p_stats%f_stats%gdata => null()
       p_stats%f_stats%fe_space => null()

    end if

    ! Nullify pointers
    p_stats%p_fe_space => null()
    p_stats%p_env      => null()

1   format(i10,6(1x,e16.8e3))
2   format(i10,10(1x,e16.8e3))

  end subroutine finalize_par_line_statistics 

end module par_statistics_names
