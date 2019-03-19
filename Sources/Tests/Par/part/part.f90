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
module partitioner_input_names
  use fempar_names
# include "debug.i90"
  implicit none
  private

  type, extends(parameter_handler_t) ::  partitioner_input_t 
     private 
   contains
     procedure :: define_parameters    => partitioner_input_define_parameters
  end type partitioner_input_t

  public :: partitioner_input_t

contains

  subroutine partitioner_input_define_parameters(this)
    implicit none
    class(partitioner_input_t), intent(inout) :: this
    integer(ip)                               :: error

    call this%add(dir_path_key, '--dir-path', '.', 'Directory of the source files', switch_ab='-d')
    call this%add(prefix_key, '--prefix', 'square', 'Name of the GiD files', switch_ab='-p')
    call this%add(dir_path_out_key, '--dir-path-out', '.', 'Output Directory', switch_ab='-o')
    call this%add(num_parts_key, '--num_parts', 16, 'Number of parts of the mesh', switch_ab='-n')
    call this%add(num_levels_distribution_key, '--num_levels', 1, 'Number of levels', switch_ab='-l')
    call this%add(num_parts_x_level_key, '--num_parts_x_level', [16,4,1,0,0], 'Number of parts per level (array of fixed size 5)', switch_ab='-npl')

    call this%add(strategy_key, '--strategy', part_kway, 'Partitioning strategy')
    call this%add(debug_key, '--debug',  0, 'Debug mode')
    call this%add(metis_option_debug_key, '--metis_debug',  2, 'Metis debug option')
    call this%add(metis_option_ufactor_key, '--metis_ufactor', 30, 'Metis ufactor option')
    call this%add(metis_option_minconn_key, '--metis_minconn', 0, 'Metis minconn option')
    call this%add(metis_option_contig_key, '--metis_contig', 1, 'Metis contig option')
    call this%add(metis_option_ctype_key, '--metis_ctype', METIS_CTYPE_SHEM, 'Metis ctype option')
    call this%add(metis_option_iptype_key, '--metis_iptype', METIS_IPTYPE_EDGE, 'Metis iptype option')
  end subroutine partitioner_input_define_parameters

end module partitioner_input_names
!==================================================================================================
!==================================================================================================

program partitioner
  use fempar_names
  use partitioner_input_names
  implicit none
  type(partitioner_input_t)              :: input
  type(ParameterList_t)    , pointer     :: parameters
  type(mesh_t)                           :: gmesh
  type(mesh_distribution_t), allocatable :: distr(:)
  type(environment_t)  , allocatable :: env(:)
  type(mesh_t)             , allocatable :: lmesh(:)
  integer(ip) :: ipart, ienv

  call fempar_init()
  call input%create()
  parameters => input%get_values()

  ! Read and partition gmesh into lmesh
  call gmesh%read(parameters)
  call gmesh%write_file_for_postprocess(parameters)
  call gmesh%create_distribution (parameters, distr, env, lmesh)

  ! Write environments
  call environment_write_files             ( parameters, env )

  ! Write partition info
  call mesh_distribution_write_files           ( parameters, distr )
  call mesh_distribution_write_for_postprocess ( parameters, gmesh, distr )

  ! Write local meshes
  call mesh_write_files                 ( parameters, lmesh )
  call mesh_write_files_for_postprocess ( parameters, lmesh )

  ! Deallocate partition objects and envs
  do ipart=1,size(distr)
     call distr(ipart)%free()
     call lmesh(ipart)%free
  end do
  deallocate (distr)
  deallocate (lmesh)

  do ienv=1,size(env)
     call env(ienv)%free()
  end do
  deallocate (env)
  
  call gmesh%free()

  call input%free()
  call fempar_finalize()

end program partitioner
