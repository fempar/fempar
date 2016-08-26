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
module test_mesh_gen_input_names
  use serial_names
# include "debug.i90"
  implicit none
  private

  type test_mesh_gen_input_t 
     private 
     type(Command_Line_Interface)  :: cli 

     type(ParameterList_t)         :: list
     type(ParameterList_t)         :: switches
     type(ParameterList_t)         :: switches_ab
     type(ParameterList_t)         :: helpers
     type(ParameterList_t)         :: required

   contains
     procedure, non_overridable             :: create         => test_mesh_gen_input_create
     procedure, non_overridable, private    :: set_default    => test_mesh_gen_input_set_default
     procedure, non_overridable, private    :: add_to_cli     => test_mesh_gen_input_add_to_cli
     procedure, non_overridable, private    :: parse          => test_mesh_gen_input_parse
     procedure, non_overridable             :: get_parameters => test_mesh_gen_input_get_parameters
     procedure, non_overridable             :: free           => test_mesh_gen_input_free
  end type test_mesh_gen_input_t

  public :: test_mesh_gen_input_t

contains

  subroutine test_mesh_gen_input_create(this)
    implicit none
    class(test_mesh_gen_input_t), intent(inout) :: this
    call this%free()
     ! Initialize Command Line Interface
    call this%cli%init(progname    = 'test_mesh_genpart',                                                     &
         &        version     = '',                                                                 &
         &        authors     = '',                                                                 &
         &        license     = '',                                                                 &
         &        description =  'FEMPAR driver to generate a uniform mesh.', &
         &        examples    = ['test_mesh_gen -h  ', 'test_mesh_gen -n  ' ])
    call this%set_default()
    call this%add_to_cli()
    call this%parse()
  end subroutine test_mesh_gen_input_create

  !==================================================================================================
  subroutine test_mesh_gen_input_set_default(this)
    implicit none
    class(test_mesh_gen_input_t), intent(inout) :: this
    integer(ip) :: error
    integer(ip) :: tmp(SPACE_DIM) ! To circumbent gfortran error

    call this%list%init()
    error = 0
    error = error + this%list%set(key = dir_path_key                  , value = '.')
    error = error + this%list%set(key = prefix_key                    , value = 'uniform')
    error = error + this%list%set(key = dir_path_out_key              , value = '.')
    error = error + this%list%set(key = number_of_dimensions_key      , value =  2)
    error = error + this%list%set(key = interpolation_order_key       , value =  1)
    tmp = [12,6,0]; error = error + this%list%set(key = number_of_cells_per_dir_key   , value = tmp)
    tmp = [3,3,1]; error = error + this%list%set(key = number_of_parts_per_dir_key   , value = tmp)
    tmp = [0,0,0]; error = error + this%list%set(key = is_dir_periodic_key           , value = tmp)

    ! tmp = [4,2,0]; error = error + this%list%set(key = number_of_cells_per_dir_key   , value = tmp)
    ! tmp = [1,1,1]; error = error + this%list%set(key = number_of_parts_per_dir_key   , value = tmp)
    ! tmp = [1,1,0]; error = error + this%list%set(key = is_dir_periodic_key           , value = tmp)

    !error = error + this%list%set(key = number_of_cells_per_dir_key   , value =  [4,2,0])
    !error = error + this%list%set(key = number_of_parts_per_dir_key   , value =  [1,1,1])
    !error = error + this%list%set(key = is_dir_periodic_key           , value =  [0,0,0])
    check(error==0)

    ! Only some of them are controlled from cli
    call this%switches%init()
    error = error + this%switches%set(key = dir_path_key    , value = '--dir-path')
    error = error + this%switches%set(key = prefix_key      , value = '--prefix')
    error = error + this%switches%set(key = dir_path_out_key, value = '--dir-path-out')
    check(error==0)

    call this%switches_ab%init()
    error = error + this%switches_ab%set(key = dir_path_key    , value = '-d')
    error = error + this%switches_ab%set(key = prefix_key      , value = '-p')
    error = error + this%switches_ab%set(key = dir_path_out_key, value = '-o')
    check(error==0)

    call this%helpers%init()
    error = error + this%helpers%set(key = dir_path_key    , value = 'Directory of the source files')
    error = error + this%helpers%set(key = prefix_key      , value = 'Name of the GiD files')
    error = error + this%helpers%set(key = dir_path_out_key, value = 'Output Directory')
    check(error==0)

    call this%required%init()
    error = error + this%required%set(key = dir_path_key    , value = .false.)
    error = error + this%required%set(key = prefix_key      , value = .false.)
    error = error + this%required%set(key = dir_path_out_key, value = .false.)
    check(error==0)

  end subroutine test_mesh_gen_input_set_default

  !==================================================================================================
  subroutine test_mesh_gen_input_free(this)
    implicit none
    class(test_mesh_gen_input_t), intent(inout) :: this
    call this%list%free()
    call this%switches%free()
    call this%switches_ab%free()
    call this%required%free()
    call this%cli%free()
   end subroutine test_mesh_gen_input_free

  !==================================================================================================
  function test_mesh_gen_input_get_parameters(this)
    implicit none
    class(test_mesh_gen_input_t), target , intent(in) :: this
    type(ParameterList_t), pointer  :: test_mesh_gen_input_get_parameters
    test_mesh_gen_input_get_parameters => this%list
  end function test_mesh_gen_input_get_parameters

  !==================================================================================================
  !
  ! The following methods can be programmed in the library looping over the entries in, e.g. switch.  
  ! To do that we need to manage data types conversions automatically. Here I'm exploiting the knowledge
  ! of the data type of each entry. We could ask fpl...
  !
  !==================================================================================================
  subroutine test_mesh_gen_input_add_to_cli(this)
    implicit none
    class(test_mesh_gen_input_t) , intent(inout) :: this
    integer(ip)        :: error
    character(len=512) :: switch, switch_ab, help, cvalue
    logical            :: required
    integer(ip)        :: ivalue

    ! IO parameters
    error = 0
    error = error + this%list%get       (key = dir_path_key , value = cvalue)
    error = error + this%switches%get   (key = dir_path_key , value = switch)
    error = error + this%switches_ab%get(key = dir_path_key , value = switch_ab)
    error = error + this%helpers%get    (key = dir_path_key , value = help)
    error = error + this%required%get   (key = dir_path_key , value = required)
    call this%cli%add(switch=trim(switch),switch_ab=trim(switch_ab), help=trim(help), &
         &            required=required,act='store',def=trim(cvalue),error=error)
    check(error==0)

    error = 0
    error = error + this%list%get       (key = prefix_key , value = cvalue)
    error = error + this%switches%get   (key = prefix_key , value = switch)
    error = error + this%switches_ab%get(key = prefix_key , value = switch_ab)
    error = error + this%helpers%get    (key = prefix_key , value = help)
    error = error + this%required%get   (key = prefix_key , value = required)
    check(error==0)
    call this%cli%add(switch=trim(switch),switch_ab=trim(switch_ab), help=trim(help), &
         &            required=required,act='store',def=trim(cvalue),error=error)
    check(error==0)

    error = 0
    error = error + this%list%get       (key = dir_path_out_key , value = cvalue)
    error = error + this%switches%get   (key = dir_path_out_key , value = switch)
    error = error + this%switches_ab%get(key = dir_path_out_key , value = switch_ab)
    error = error + this%helpers%get    (key = dir_path_out_key , value = help)
    error = error + this%required%get   (key = dir_path_out_key , value = required)
    check(error==0)
    call this%cli%add(switch=trim(switch),switch_ab=trim(switch_ab), help=trim(help), &
         &            required=required,act='store',def=trim(cvalue),error=error)
    check(error==0)

  end subroutine test_mesh_gen_input_add_to_cli

  subroutine test_mesh_gen_input_parse(this)
    implicit none
    class(test_mesh_gen_input_t), intent(inout) :: this
    integer(ip)    :: istat
    character(512) :: switch, cvalue
    integer(ip)    :: ivalue

    call this%cli%parse(error=istat); check(istat==0)

    istat = this%switches%get(key = dir_path_key , value = switch)
    check(istat==0)
    if (this%cli%is_passed(switch=switch)) then
       call this%cli%get(switch=switch, val=cvalue, error=istat); check(istat==0)
       istat = this%list%set(key = dir_path_key, value=cvalue)
    end if

    istat = this%switches%get(key = prefix_key , value = switch)
    check(istat==0)
    if (this%cli%is_passed(switch=switch)) then
       call this%cli%get(switch=switch, val=cvalue, error=istat); check(istat==0)
       istat = this%list%set(key = prefix_key, value=cvalue)
    end if

    istat = this%switches%get(key = dir_path_out_key , value = switch)
    check(istat==0)
    if (this%cli%is_passed(switch=switch)) then
       call this%cli%get(switch=switch, val=cvalue, error=istat); check(istat==0)
       istat = this%list%set(key = dir_path_out_key, value=cvalue)
    end if

  end subroutine test_mesh_gen_input_parse  

end module test_mesh_gen_input_names 

!==================================================================================================
!==================================================================================================
!==================================================================================================
!==================================================================================================

program partitioner
  use serial_names
  use test_mesh_gen_input_names
  implicit none
  type(test_mesh_gen_input_t)    :: input
  type(ParameterList_t), pointer :: parameters
    integer(ip)                  :: num_local_cells
    integer(ip)                  :: num_ghost_cells
    integer(ip)                  :: num_local_vefs
    integer(ip) , allocatable    :: ptr_vefs_per_cell(:)            ! Size = num_local_cells + 1
    integer(ip) , allocatable    :: lst_vefs_lids(:)                ! Size = ptr_vefs_per_cell(num_local_cells+1)-1
    integer(ip) , allocatable    :: vef_lid_to_vertex_lid(:)        ! Size = num_local_vefs (-1 if vef_lid is not a vertex)
    integer(igp), allocatable    :: cell_gids(:)                    ! Size = num_local_cells 
    integer(igp), allocatable    :: vefs_gids(:)                    ! Size = num_local_vefs
    integer(ip) :: icell, error
    integer(ip) :: n(SPACE_DIM)

  call fempar_init()
  call input%create()
  parameters => input%get_parameters()

  error = parameters%get(key = number_of_cells_per_dir_key, value = n)
  write(*,*) n

  call uniform_hex_mesh_generator_generate_connectivities(parameters,        &
       num_local_cells,       &
       num_ghost_cells,       &
       num_local_vefs,        &
       ptr_vefs_per_cell,     &
       lst_vefs_lids,         &
       vef_lid_to_vertex_lid, &
       cell_gids,             &
       vefs_gids,5)
 
  write(*,*) num_local_cells, num_ghost_cells, num_local_vefs
  do icell=1,num_local_cells+num_ghost_cells
     write(*,*) icell, ptr_vefs_per_cell(icell), ptr_vefs_per_cell(icell+1)-1, lst_vefs_lids(ptr_vefs_per_cell(icell):ptr_vefs_per_cell(icell+1)-1)
  end do

  call memfree(ptr_vefs_per_cell,     __FILE__,__LINE__)
  call memfree(lst_vefs_lids,         __FILE__,__LINE__)
  !call memfree(vef_lid_to_vertex_lid, __FILE__,__LINE__)
  !call memfree(cell_gids,             __FILE__,__LINE__)
  !call memfree(vefs_gids,             __FILE__,__LINE__) 

  call input%free()
  call fempar_finalize()

end program partitioner
