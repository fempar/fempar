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
    type(ParameterList_t), pointer :: list, switches, switches_ab, helpers, required
    integer(ip) :: error

    list        => this%get_values()
    switches    => this%get_switches()
    switches_ab => this%get_switches_ab()
    helpers     => this%get_helpers()
    required    => this%get_required()

    error = list%set(key = dir_path_key            , value = '.') ; check(error==0)
    error = list%set(key = prefix_key              , value = 'square') ; check(error==0)
    error = list%set(key = dir_path_out_key        , value = '.') ; check(error==0)
    error = list%set(key = num_parts_key           , value =  16)              ; check(error==0)
    error = list%set(key = num_levels_distribution_key          , value =  1)               ; check(error==0)
    error = list%set(key = num_parts_per_level_key , value =  [16,4,1,0,0])    ; check(error==0)
    error = list%set(key = strategy_key            , value = part_kway)        ; check(error==0)
    error = list%set(key = debug_key               , value =  0)               ; check(error==0)
    error = list%set(key = metis_option_debug_key  , value =  2)               ; check(error==0)
    error = list%set(key = metis_option_ufactor_key, value = 30)               ; check(error==0)
    error = list%set(key = metis_option_minconn_key, value =  0)               ; check(error==0)
    error = list%set(key = metis_option_contig_key , value =  1)               ; check(error==0)
    error = list%set(key = metis_option_ctype_key  , value = METIS_CTYPE_SHEM) ; check(error==0)
    error = list%set(key = metis_option_iptype_key , value = METIS_IPTYPE_EDGE); check(error==0)

    ! Only some of them are controlled from cli
    error = switches%set(key = dir_path_key    , value = '--dir-path')                  ; check(error==0)
    error = switches%set(key = prefix_key      , value = '--prefix')                    ; check(error==0)
    error = switches%set(key = dir_path_out_key, value = '--dir-path-out')              ; check(error==0)
    error = switches%set(key = num_parts_key   , value = '--num_parts')                 ; check(error==0)
    error = switches%set(key = num_levels_distribution_key  , value = '--num_levels')                ; check(error==0)
    error = switches%set(key = num_parts_per_level_key, value = '--num_parts_per_level'); check(error==0)

    error = switches_ab%set(key = dir_path_key    , value = '-d')             ; check(error==0)
    error = switches_ab%set(key = prefix_key      , value = '-p')             ; check(error==0)
    error = switches_ab%set(key = dir_path_out_key, value = '-o')             ; check(error==0)
    error = switches_ab%set(key = num_parts_key   , value = '-n')             ; check(error==0)
    error = switches_ab%set(key = num_levels_distribution_key  , value = '-l')             ; check(error==0)
    error = switches_ab%set(key = num_parts_per_level_key   , value = '-npl') ; check(error==0)

    error = helpers%set(key = dir_path_key    , value = 'Directory of the source files')                               ; check(error==0)
    error = helpers%set(key = prefix_key      , value = 'Name of the GiD files')                                       ; check(error==0)
    error = helpers%set(key = dir_path_out_key, value = 'Output Directory')                                            ; check(error==0)
    error = helpers%set(key = num_parts_key   , value = 'Number of parts of the mesh')                                 ; check(error==0)
    error = helpers%set(key = num_levels_distribution_key  , value = 'Number of levels')                                            ; check(error==0)
    error = helpers%set(key = num_parts_per_level_key   , value = 'Number of parts per level (array of fixed size 5)') ; check(error==0)

    error = required%set(key = dir_path_key    , value = .false.)           ; check(error==0)
    error = required%set(key = prefix_key      , value = .false.)           ; check(error==0)
    error = required%set(key = dir_path_out_key, value = .false.)           ; check(error==0)
    error = required%set(key = num_parts_key   , value = .false.)           ; check(error==0)
    error = required%set(key = num_levels_distribution_key  , value = .false.)           ; check(error==0)
    error = required%set(key = num_parts_per_level_key  , value = .false.)  ; check(error==0)

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
