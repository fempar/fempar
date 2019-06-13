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

program partitioner
  use fempar_names
  implicit none
  type(mesh_t) :: gmesh
  type(mesh_partitioner_t) :: mesh_partitioner
  type(parameterlist_t), pointer :: parameters
  call fempar_init()
  call parameter_handler%process_parameters()
  parameters => parameter_handler%get_values()
  call gmesh%read_fempar_gid_problem_type_format(parameters)
  call gmesh%write_gid_postprocess_format(parameters)
  call mesh_partitioner%create(gmesh, parameters)
  call mesh_partitioner%partition_mesh()
  call mesh_partitioner%write_mesh_parts_fempar_gid_problem_type_format(parameters)
  call mesh_partitioner%write_mesh_partition_gid_postprocess_format(parameters)
  call mesh_partitioner%free()
  call gmesh%free()
  call fempar_finalize()
end program partitioner
