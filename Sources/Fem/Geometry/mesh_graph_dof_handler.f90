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
module mesh_graph_dof_handler
  ! Serial modules
  use types
  use fem_graph_class
  use fem_mesh_class
  use mesh_graph
  use dof_handler_class
  use fem_blocks_class
  use fem_space_class

# include "debug.i90"

  implicit none
  private

  ! Function
  public :: fem_mesh_graph_dof_handler_create_block

contains

  !=============================================================================
  subroutine fem_mesh_graph_dof_handler_create_block(gtype,mesh,dhand,graph,femsp)
    ! Parameters
    type(fem_mesh)     , intent(in)         :: mesh
    type(dof_handler)  , intent(in)         :: dhand
    integer(ip)        , intent(in)         :: gtype(dhand%blocks%nb) 
    type(fem_graph)    , intent(out)        :: graph(dhand%blocks%nb,dhand%blocks%nb)
    type(fem_space), optional, intent(in)   :: femsp

    ! Locals
    integer(ip)     :: ibloc,jbloc
    type(fem_graph) :: mesh_graph
    
    ! Point to mesh's parallel partition object
    do ibloc=1,dhand%blocks%nb
       do jbloc=1,dhand%blocks%nb
          if ( ibloc == jbloc ) then
             graph(ibloc,jbloc)%type = gtype(ibloc)
          else ! ibloc /= jbloc
             graph(ibloc,jbloc)%type = csr
          end if
       end do
    end do

    ! Construct mesh graph
    call mesh_to_graph( scal, csr, 1, 1, mesh, mesh_graph)

    ! Fill p_graph%f_graph
    do ibloc=1,dhand%blocks%nb
       do jbloc=1,dhand%blocks%nb
          if ( ibloc == jbloc ) then
             call graph_to_dof_graph(gtype(ibloc),ibloc,jbloc,dhand,mesh,mesh_graph, &
                  &                  graph(ibloc,jbloc),femsp) 
          else ! ibloc /= jbloc
             call graph_to_dof_graph(csr,ibloc,jbloc,dhand,mesh,mesh_graph, &
                  &                  graph(ibloc,jbloc),femsp)
          end if
       end do
    end do

    ! Destruct graph
    call fem_graph_free(mesh_graph)

  end subroutine fem_mesh_graph_dof_handler_create_block

end module mesh_graph_dof_handler
