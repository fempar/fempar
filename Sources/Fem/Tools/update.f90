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
module fem_update_names
  use types_names
  use memor_names
  use fem_space_names
  use fem_vector_names
  use fem_block_vector_names
  use fem_conditions_names
  use analytical_names
  use interpolation_tools_names
  implicit none
# include "debug.i90"
  private

  interface fem_update_solution
     module procedure fem_update_solution_mono, fem_update_solution_block
  end interface fem_update_solution
     

  ! Functions
  public :: fem_update_strong_dirichlet_bcond, fem_update_analytical_bcond, fem_update_solution
 
contains
  
  !==================================================================================================
  subroutine fem_update_strong_dirichlet_bcond( fspac, fcond )
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine updates Dirichlet boundary conditions in unkno from fem_conditions values.  !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    type(fem_space_t)     , intent(inout) :: fspac
    type(fem_conditions_t), intent(in)    :: fcond
    ! Locals
    integer(ip) :: ielem, iobje, ivar, inode, l_node, gvar, lobje, prob

    do ielem = 1, fspac%g_trian%num_elems
       prob = fspac%lelem(ielem)%problem
       do ivar=1, fspac%dof_handler%problems(prob)%p%nvars
          gvar=fspac%dof_handler%problems(prob)%p%l2g_var(ivar)
          do iobje = 1,fspac%lelem(ielem)%p_geo_info%nobje
             lobje = fspac%g_trian%elems(ielem)%objects(iobje)
             do inode = fspac%lelem(ielem)%nodes_object(ivar)%p%p(iobje), &
                  &     fspac%lelem(ielem)%nodes_object(ivar)%p%p(iobje+1)-1 
                l_node = fspac%lelem(ielem)%nodes_object(ivar)%p%l(inode)
                if ( fspac%lelem(ielem)%bc_code(ivar,iobje) /= 0 ) then
                   fspac%lelem(ielem)%unkno(l_node,ivar,1) = fcond%valu(gvar,lobje)
                end if
             end do
          end do
       end do
    end do

  end subroutine fem_update_strong_dirichlet_bcond

  !==================================================================================================
  subroutine fem_update_analytical_bcond(vars_of_unk,case,ctime,fspac,caset,t)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine updates Dirichlet boundary conditions in unkno from an analytical solution. !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    integer(ip)          , intent(in)    :: vars_of_unk(:)
    integer(ip)          , intent(in)    :: case
    real(rp)             , intent(in)    :: ctime
    type(fem_space_t)      , intent(inout) :: fspac
    integer(ip), optional, intent(in)    :: caset,t
    ! Locals
    integer(ip) :: ielem,prob,ndime,iobje,lobje,inode,lnode
    integer(ip) :: nvars,ivar,gvar,gnode,unode,cnt
    real(rp)    :: part(3)
    real(rp), allocatable :: coord(:,:),param(:)

    if(case>0) then

       nvars = size(vars_of_unk,1)
       ndime = fspac%g_trian%num_dims

       ! Allocate parameters
       if(nvars==1) then
          call memalloc(10,param,__FILE__,__LINE__)
       else
          call memalloc(30,param,__FILE__,__LINE__)
       end if


       do ielem = 1, fspac%g_trian%num_elems
          prob  = fspac%lelem(ielem)%problem
          gnode = fspac%lelem(ielem)%p_geo_info%nnode
          cnt   = 0

          do ivar=vars_of_unk(1),vars_of_unk(nvars)

             ! Global variable
             cnt = cnt+1
             gvar=fspac%dof_handler%problems(prob)%p%l2g_var(ivar)

             ! Interpolate coordinates
             unode = fspac%lelem(ielem)%f_inf(ivar)%p%nnode
             call memalloc(ndime,unode,coord,__FILE__,__LINE__)
             call interpolate(ndime,gnode,unode,fspac%lelem(ielem)%inter(ivar)%p, &
                  &           fspac%g_trian%elems(ielem)%coordinates,coord)

             do iobje = 1,fspac%lelem(ielem)%p_geo_info%nobje
                lobje = fspac%g_trian%elems(ielem)%objects(iobje)

                if ( fspac%lelem(ielem)%bc_code(ivar,iobje) /= 0 ) then

                   do inode = fspac%lelem(ielem)%nodes_object(ivar)%p%p(iobje), &
                        &     fspac%lelem(ielem)%nodes_object(ivar)%p%p(iobje+1)-1 
                      lnode = fspac%lelem(ielem)%nodes_object(ivar)%p%l(inode)

                      call analytical_field(case,ndime,coord(:,lnode),ctime,param)

                      if(present(caset)) then
                         if(caset>0) then
                            call analytical_field(caset,ndime,coord(:,lnode),ctime,part)
                            if(present(t)) then
                               fspac%lelem(ielem)%unkno(lnode,ivar,1) =  param(cnt)*part(t)
                            else
                               fspac%lelem(ielem)%unkno(lnode,ivar,1) =  param(cnt)*part(1)
                            end if
                         else
                            fspac%lelem(ielem)%unkno(lnode,ivar,1) =  param(cnt)
                         end if
                      else
                         fspac%lelem(ielem)%unkno(lnode,ivar,1) =  param(cnt)
                      end if

                   end do
                end if
             end do

             ! Deallocate auxiliar coordinates
             call memfree(coord,__FILE__,__LINE__)

          end do
       end do

       ! Deallocate params
       call memfree(param,__FILE__,__LINE__)

    end if
    
  end subroutine fem_update_analytical_bcond

  !==================================================================================================
  subroutine fem_update_solution_mono(fevec,fspac,iblock)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine stores the solution from a fem_vector into unkno.                           !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    type(fem_vector_t)     , intent(in)    :: fevec   
    type(fem_space_t)      , intent(inout) :: fspac
    integer(ip), optional, intent(in)    :: iblock
    ! Locals
    integer(ip) :: ielem,iblock_,iprob,nvapb,ivar,lvar,inode,idof

    iblock_ = 1
    if ( present(iblock) ) iblock_ = iblock
    
    ! Loop over elements
    do ielem = 1, fspac%g_trian%num_elems
       iprob = fspac%lelem(ielem)%problem
       nvapb = fspac%dof_handler%prob_block(iblock_,iprob)%nd1
       
       ! Loop over problem and block variables
       do ivar = 1, nvapb
          lvar = fspac%dof_handler%prob_block(iblock_,iprob)%a(ivar)

          ! Loop over elemental nodes
          do inode = 1,fspac%lelem(ielem)%f_inf(lvar)%p%nnode
             idof = fspac%lelem(ielem)%elem2dof(inode,lvar)
             
             if(idof/=0) then

                ! Update unkno
                fspac%lelem(ielem)%unkno(inode,lvar,1) = fevec%b(idof)

             end if

          end do
       end do
    end do
    
  end subroutine fem_update_solution_mono

  !==================================================================================================
  subroutine fem_update_solution_block(blvec,fspac)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine stores the solution from a fem_vector into unkno.                           !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    type(fem_block_vector_t), intent(in)    :: blvec   
    type(fem_space_t)       , intent(inout) :: fspac
    ! Locals
    integer(ip) :: iblock

    ! Loop over blocks
    do iblock = 1,blvec%nblocks
       
       ! Call monolithic update
       call fem_update_solution_mono(blvec%blocks(iblock),fspac,iblock)

    end do
    
  end subroutine fem_update_solution_block

end module fem_update_names
