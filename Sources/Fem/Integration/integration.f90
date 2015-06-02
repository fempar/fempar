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
# include "debug.i90"
module integration_names
  use types
  use assembly_names
  use integrable_names
  use problem_names
  use integration_tools_names
  use femap_interp
  use fem_space_names
  use assembly_names
  use fem_block_matrix_names
  use fem_matrix_names
  use fem_block_vector_names
  use fem_vector_names
  implicit none
  private

  public :: volume_integral

contains

  subroutine volume_integral(femsp,res1,res2)
    implicit none
    ! Parameters
    type(fem_space)  , intent(inout) :: femsp
    class(integrable), intent(inout) :: res1
    class(integrable), optional, intent(inout) :: res2

    ! Locals
    integer(ip) :: ielem,ivar,nvars 
    class(discrete_problem) , pointer :: discrete
    integer(ip) :: start(femsp%dof_handler%nvars_global)

    ! Main element loop
    do ielem=1,femsp%g_trian%num_elems

       nvars = femsp%lelem(ielem)%num_vars
       ! Compute integration tools on ielem for each ivar (they all share the quadrature inside integ)
       do ivar=1,nvars
          call volume_integrator_update(femsp%lelem(ielem)%integ(ivar)%p,femsp%g_trian%elems(ielem)%coordinates)
       end do

       ! Starting position for each dof
       call pointer_variable(femsp%lelem(ielem),femsp%dof_handler, start(1:nvars+1) )

       ! Compute element matrix and rhs
       discrete => femsp%approximations(femsp%lelem(ielem)%approximation)%p
       call discrete%matvec(femsp%lelem(ielem)%integ,femsp%lelem(ielem)%unkno,start(1:nvars+1), &
            &               femsp%lelem(ielem)%p_mat,femsp%lelem(ielem)%p_vec)

       ! Apply boundary conditions
       call impose_strong_dirichlet_data (femsp%lelem(ielem),femsp%dof_handler) 

       ! Assembly first contribution
       select type(res1)
       class is(fem_matrix)
          call assembly_element_matrix(femsp%lelem(ielem),femsp%dof_handler,start(1:nvars+1),res1) 
       class is(fem_vector)
          call assembly_element_vector(femsp%lelem(ielem),femsp%dof_handler,start(1:nvars+1),res1) 
       ! class is(fem_block_matrix)
       !    call assembly_element_matrix_block(femsp%lelem(ielem),femsp%dof_handler,start(1:nvars+1),res1) 
       ! class is(fem_block_vector)
       !    call assembly_element_vector_block(femsp%lelem(ielem),femsp%dof_handler,start(1:nvars+1),res1) 
       class default
          ! class not yet implemented
          check(.false.)
       end select

       ! Assembly second contribution if present
       if(present(res2)) then
          select type(res2)
             class is(fem_matrix)
             call assembly_element_matrix(femsp%lelem(ielem),femsp%dof_handler,start(1:nvars+1),res2) 
             class is(fem_vector)
             call assembly_element_vector(femsp%lelem(ielem),femsp%dof_handler,start(1:nvars+1),res2) 
             ! class is(fem_block_matrix)
             ! call assembly_element_matrix_block(femsp%lelem(ielem),femsp%dof_handler,start(1:nvars+1),res2) 
             ! class is(fem_block_vector)
             ! call assembly_element_vector_block(femsp%lelem(ielem),femsp%dof_handler,start(1:nvars+1),res2) 
             class default
                ! class not yet implemented 
             check(.false.)
          end select
       end if

    end do

  end subroutine volume_integral

end module integration_names
