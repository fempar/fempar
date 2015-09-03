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
module update_names
  use types_names
  use memor_names
  use fe_space_names
  use vector_names
  use block_vector_names
  use conditions_names
  use analytical_function_names
  use interpolation_tools_names
  use abstract_vector_names
  implicit none
# include "debug.i90"
  private     

  ! Functions
  public :: update_strong_dirichlet_bcond, update_analytical_bcond, update_solution, &
       &    update_solution_mono, update_nonlinear_solution, update_analytical_initial, &
       &    update_initialize, update_initialize_mono
 
contains
  
  !==================================================================================================
  subroutine update_strong_dirichlet_bcond( fe_space, fcond )
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine updates Dirichlet boundary conditions in unkno from conditions values.  !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    type(fe_space_t)     , intent(inout) :: fe_space
    type(conditions_t), intent(in)    :: fcond
    ! Locals
    integer(ip) :: ielem, iobje, ivar, inode, l_node, gvar, lobje, prob

    do ielem = 1, fe_space%g_trian%num_elems
       prob = fe_space%finite_elements(ielem)%problem
       do ivar=1, fe_space%dof_descriptor%problems(prob)%p%nvars
          gvar=fe_space%dof_descriptor%problems(prob)%p%l2g_var(ivar)
          do iobje = 1,fe_space%finite_elements(ielem)%p_geo_reference_element%nvef
             lobje = fe_space%g_trian%elems(ielem)%vefs(iobje)
             do inode = fe_space%finite_elements(ielem)%nodes_per_vef(ivar)%p%p(iobje), &
                  &     fe_space%finite_elements(ielem)%nodes_per_vef(ivar)%p%p(iobje+1)-1 
                l_node = fe_space%finite_elements(ielem)%nodes_per_vef(ivar)%p%l(inode)
                if ( fe_space%finite_elements(ielem)%bc_code(ivar,iobje) /= 0 ) then
                   fe_space%finite_elements(ielem)%unkno(l_node,ivar,1) = fcond%valu(gvar,lobje)
                end if
             end do
          end do
       end do
    end do

  end subroutine update_strong_dirichlet_bcond

  !==================================================================================================
  subroutine update_analytical_bcond(vars_of_unk,ctime,fe_space,tvar)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine updates Dirichlet boundary conditions in unkno from an analytical solution. !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    integer(ip)          , intent(in)    :: vars_of_unk(:)
    real(rp)             , intent(in)    :: ctime
    type(fe_space_t)     , intent(inout) :: fe_space
    integer(ip), optional, intent(in)    :: tvar
    ! Locals
    integer(ip) :: ielem,prob,ndime,iobje,lobje,inode,lnode
    integer(ip) :: nvars,ivar,gvar,gnode,unode,cnt
    real(rp)    :: param(11)
    real(rp), allocatable :: coord(:,:)

    nvars = size(vars_of_unk,1)
    ndime = fe_space%g_trian%num_dims

    do ielem = 1, fe_space%g_trian%num_elems
       prob  = fe_space%finite_elements(ielem)%problem
       gnode = fe_space%finite_elements(ielem)%p_geo_reference_element%nnode
       cnt   = 0

       do ivar=vars_of_unk(1),vars_of_unk(nvars)

          if(fe_space%finite_elements(ielem)%p_analytical_code%a(ivar,1)>0) then

             param = 0.0_rp

             ! Global variable
             cnt = cnt+1
             gvar=fe_space%dof_descriptor%problems(prob)%p%l2g_var(ivar)

             ! Interpolate coordinates
             unode = fe_space%finite_elements(ielem)%reference_element_vars(ivar)%p%nnode
             call memalloc(ndime,unode,coord,__FILE__,__LINE__)
             call interpolate(ndime,gnode,unode,fe_space%finite_elements(ielem)%inter(ivar)%p, &
                  &           fe_space%g_trian%elems(ielem)%coordinates,coord)

             do iobje = 1,fe_space%finite_elements(ielem)%p_geo_reference_element%nvef
                lobje = fe_space%g_trian%elems(ielem)%vefs(iobje)

                if ( fe_space%finite_elements(ielem)%bc_code(ivar,iobje) /= 0 ) then

                   do inode = fe_space%finite_elements(ielem)%nodes_per_vef(ivar)%p%p(iobje), &
                        &     fe_space%finite_elements(ielem)%nodes_per_vef(ivar)%p%p(iobje+1)-1 
                      lnode = fe_space%finite_elements(ielem)%nodes_per_vef(ivar)%p%l(inode)

                      ! Evaluate analytical unknown
                      call evaluate_analytical(fe_space%finite_elements(ielem)%p_analytical_code%a(ivar,1), &
                           &                   fe_space%finite_elements(ielem)%p_analytical_code%a(ivar,2), &
                           &                   ndime,coord(:,lnode),ctime,param,tvar=tvar)

                      fe_space%finite_elements(ielem)%unkno(lnode,ivar,1) = param(1)

                   end do
                end if
             end do

             ! Deallocate auxiliar coordinates
             call memfree(coord,__FILE__,__LINE__)

          end if

       end do
    end do
    
  end subroutine update_analytical_bcond

  !==================================================================================================
  subroutine update_analytical_initial(vars_of_unk,ctime,fe_space,tvar)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine updates Dirichlet boundary conditions in unkno from an analytical solution. !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    integer(ip)          , intent(in)    :: vars_of_unk(:)
    real(rp)             , intent(in)    :: ctime
    type(fe_space_t)     , intent(inout) :: fe_space
    integer(ip), optional, intent(in)    :: tvar
    ! Locals
    integer(ip) :: ielem,prob,ndime,iobje,lobje,inode,lnode
    integer(ip) :: nvars,ivar,gvar,gnode,unode,cnt
    real(rp)    :: param(11)
    real(rp), allocatable :: coord(:,:)

    nvars = size(vars_of_unk,1)
    ndime = fe_space%g_trian%num_dims

    do ielem = 1, fe_space%g_trian%num_elems
       prob  = fe_space%finite_elements(ielem)%problem
       gnode = fe_space%finite_elements(ielem)%p_geo_reference_element%nnode
       cnt   = 0

       do ivar=vars_of_unk(1),vars_of_unk(nvars)

          if(fe_space%finite_elements(ielem)%p_analytical_code%a(ivar,1)>0) then

             param = 0.0_rp

             ! Global variable
             cnt = cnt+1
             gvar=fe_space%dof_descriptor%problems(prob)%p%l2g_var(ivar)

             ! Interpolate coordinates
             unode = fe_space%finite_elements(ielem)%reference_element_vars(ivar)%p%nnode
             call memalloc(ndime,unode,coord,__FILE__,__LINE__)
             call interpolate(ndime,gnode,unode,fe_space%finite_elements(ielem)%inter(ivar)%p, &
                  &           fe_space%g_trian%elems(ielem)%coordinates,coord)

             do iobje = 1,fe_space%finite_elements(ielem)%reference_element_vars(ivar)%p%nvef+1

                do inode = fe_space%finite_elements(ielem)%nodes_per_vef(ivar)%p%p(iobje), &
                     &     fe_space%finite_elements(ielem)%nodes_per_vef(ivar)%p%p(iobje+1)-1 
                   lnode = fe_space%finite_elements(ielem)%nodes_per_vef(ivar)%p%l(inode)

                   ! Evaluate analytical unknown
                   call evaluate_analytical(fe_space%finite_elements(ielem)%p_analytical_code%a(ivar,1), &
                        &                   fe_space%finite_elements(ielem)%p_analytical_code%a(ivar,2), &
                        &                   ndime,coord(:,lnode),ctime,param,tvar=tvar)

                   fe_space%finite_elements(ielem)%unkno(lnode,ivar,1) = param(1)

                end do
             end do

             ! Deallocate auxiliar coordinates
             call memfree(coord,__FILE__,__LINE__)

          end if

       end do
    end do
    
  end subroutine update_analytical_initial

  !==================================================================================================
  subroutine update_solution(vec,fe_space)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine stores the solution from a abstract_vector into unkno.                         !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(abstract_vector_t), intent(in)    :: vec   
    type(fe_space_t)     , intent(inout) :: fe_space

    select type(vec)
    class is(vector_t)
       call update_solution_mono(vec,fe_space)
    class is(block_vector_t)
       call update_solution_block(vec,fe_space)
    class default
       write(*,*) 'update_solution:: vec type not supported'
       check(.false.)
    end select
       
  end subroutine update_solution

  !==================================================================================================
  subroutine update_solution_mono(fevec,fe_space,iblock)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine stores the solution from a vector into unkno.                           !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    type(vector_t)     , intent(in)    :: fevec   
    type(fe_space_t)      , intent(inout) :: fe_space
    integer(ip), optional, intent(in)    :: iblock
    ! Locals
    integer(ip) :: ielem,iblock_,iprob,nvapb,ivar,lvar,inode,idof

    iblock_ = 1
    if ( present(iblock) ) iblock_ = iblock
    
    ! Loop over elements
    do ielem = 1, fe_space%g_trian%num_elems
       iprob = fe_space%finite_elements(ielem)%problem
       nvapb = fe_space%dof_descriptor%prob_block(iblock_,iprob)%nd1
       
       ! Loop over problem and block variables
       do ivar = 1, nvapb
          lvar = fe_space%dof_descriptor%prob_block(iblock_,iprob)%a(ivar)

          ! Loop over elemental nodes
          do inode = 1,fe_space%finite_elements(ielem)%reference_element_vars(lvar)%p%nnode
             idof = fe_space%finite_elements(ielem)%elem2dof(inode,lvar)
             
             if(idof/=0) then

                ! Update unkno
                fe_space%finite_elements(ielem)%unkno(inode,lvar,1) = fevec%b(idof)

             end if

          end do
       end do
    end do
    
  end subroutine update_solution_mono

  !==================================================================================================
  subroutine update_solution_block(blvec,fe_space)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine stores the solution from a vector into unkno.                               !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    type(block_vector_t), intent(in)    :: blvec   
    type(fe_space_t)       , intent(inout) :: fe_space
    ! Locals
    integer(ip) :: iblock

    ! Loop over blocks
    do iblock = 1,blvec%nblocks
       
       ! Call monolithic update
       call update_solution_mono(blvec%blocks(iblock),fe_space,iblock)

    end do
    
  end subroutine update_solution_block

  !==================================================================================================
  subroutine update_nonlinear_solution(fe_space,working_vars,origin,current)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine stores the previous nonlinear solution.                                     !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    type(fe_space_t)     , intent(inout) :: fe_space
    integer(ip)          , intent(in)    :: working_vars(:)
    integer(ip), optional, intent(in)    :: origin,current
    ! Locals
    integer(ip) :: ielem,ivar,origin_,current_,nvars

    nvars = size(working_vars,1)

    ! Set stored steps
    if(present(origin)) then
       origin_ = origin
    else
       origin_ = 1
    end if
    if(present(current)) then
       current_ = current
    else
       current_ = 2
    end if
    
    ! Update unkno
    do ielem = 1, fe_space%g_trian%num_elems
       do ivar = 1,nvars
          fe_space%finite_elements(ielem)%unkno(:,working_vars(ivar),current_) = &
               &    fe_space%finite_elements(ielem)%unkno(:,working_vars(ivar),origin_)
       end do
    end do
    
  end subroutine update_nonlinear_solution

  !==================================================================================================
  subroutine update_initialize(vec,fe_space)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine stores the solution from a abstract_vector into unkno.                         !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(abstract_vector_t), intent(inout) :: vec   
    type(fe_space_t)     , intent(in)    :: fe_space

    select type(vec)
    class is(vector_t)
       call update_initialize_mono(vec,fe_space)
    class is(block_vector_t)
       call update_initialize_block(vec,fe_space)
    class default
       write(*,*) 'update_solution:: vec type not supported'
       check(.false.)
    end select
       
  end subroutine update_initialize

  !==================================================================================================
  subroutine update_initialize_mono(fevec,fe_space,iblock)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine stores the solution from a vector into unkno.                           !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    type(vector_t)       , intent(inout) :: fevec   
    type(fe_space_t)     , intent(in)    :: fe_space
    integer(ip), optional, intent(in)    :: iblock
    ! Locals
    integer(ip) :: ielem,iblock_,iprob,nvapb,ivar,lvar,inode,idof

    iblock_ = 1
    if ( present(iblock) ) iblock_ = iblock
    
    ! Loop over elements
    do ielem = 1, fe_space%g_trian%num_elems
       iprob = fe_space%finite_elements(ielem)%problem
       nvapb = fe_space%dof_descriptor%prob_block(iblock_,iprob)%nd1
       
       ! Loop over problem and block variables
       do ivar = 1, nvapb
          lvar = fe_space%dof_descriptor%prob_block(iblock_,iprob)%a(ivar)

          ! Loop over elemental nodes
          do inode = 1,fe_space%finite_elements(ielem)%reference_element_vars(lvar)%p%nnode
             idof = fe_space%finite_elements(ielem)%elem2dof(inode,lvar)
             
             if(idof/=0) then

                ! Update unkno
                fevec%b(idof) = fe_space%finite_elements(ielem)%unkno(inode,lvar,1)

             end if

          end do
       end do
    end do
    
  end subroutine update_initialize_mono

  !==================================================================================================
  subroutine update_initialize_block(blvec,fe_space)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine stores the solution from a vector into unkno.                               !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    type(block_vector_t), intent(inout) :: blvec   
    type(fe_space_t)    , intent(in)    :: fe_space
    ! Locals
    integer(ip) :: iblock

    ! Loop over blocks
    do iblock = 1,blvec%nblocks
       
       ! Call monolithic update
       call update_initialize_mono(blvec%blocks(iblock),fe_space,iblock)

    end do
    
  end subroutine update_initialize_block

end module update_names
