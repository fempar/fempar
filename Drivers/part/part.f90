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

  type partitioner_input_t 
     private 
     type(Command_Line_Interface)  :: cli 

     type(ParameterList_t)         :: list
     type(ParameterList_t)         :: switches
     type(ParameterList_t)         :: switches_ab
     type(ParameterList_t)         :: helpers
     type(ParameterList_t)         :: required

   contains
     procedure, non_overridable             :: create         => partitioner_input_create
     procedure, non_overridable, private    :: set_default    => partitioner_input_set_default
     procedure, non_overridable, private    :: add_to_cli     => partitioner_input_add_to_cli
     procedure, non_overridable, private    :: parse          => partitioner_input_parse
     procedure, non_overridable             :: get_parameters => partitioner_input_get_parameters
     procedure, non_overridable             :: free           => partitioner_input_free
  end type partitioner_input_t

  public :: partitioner_input_t

contains

  subroutine partitioner_input_create(this)
    implicit none
    class(partitioner_input_t), intent(inout) :: this
    call this%free()
     ! Initialize Command Line Interface
    call this%cli%init(progname    = 'part',                                                     &
         &        version     = '',                                                                 &
         &        authors     = '',                                                                 &
         &        license     = '',                                                                 &
         &        description =  'FEMPAR driver to part a GiD mesh.', &
         &        examples    = ['part -h  ', 'part -n  ' ])
    call this%set_default()
    call this%add_to_cli()
    call this%parse()
  end subroutine partitioner_input_create

  !==================================================================================================
  subroutine partitioner_input_set_default(this)
    implicit none
    class(partitioner_input_t), intent(inout) :: this
    integer(ip) :: error

    call this%list%init()
    error = 0
    error = error + this%list%set(key = dir_path_key            , value = '.')
    error = error + this%list%set(key = prefix_key              , value = 'square')
    error = error + this%list%set(key = dir_path_out_key        , value = '.')
    error = error + this%list%set(key = num_parts_key           , value =  16)
    error = error + this%list%set(key = num_levels_key          , value =  1)
    error = error + this%list%set(key = num_parts_per_level_key , value =  [16,4,1,0,0])
    error = error + this%list%set(key = strategy_key            , value = part_kway)
    error = error + this%list%set(key = debug_key               , value =  0)
    error = error + this%list%set(key = metis_option_debug_key  , value =  2)
    error = error + this%list%set(key = metis_option_ufactor_key, value = 30)
    error = error + this%list%set(key = metis_option_minconn_key, value =  0)
    error = error + this%list%set(key = metis_option_contig_key , value =  1)
    error = error + this%list%set(key = metis_option_ctype_key  , value = METIS_CTYPE_SHEM) ! METIS_CTYPE_RM
    error = error + this%list%set(key = metis_option_iptype_key , value = METIS_IPTYPE_EDGE)
    check(error==0)

    ! Only some of them are controlled from cli
    call this%switches%init()
    error = error + this%switches%set(key = dir_path_key    , value = '--dir-path')
    error = error + this%switches%set(key = prefix_key      , value = '--prefix')
    error = error + this%switches%set(key = dir_path_out_key, value = '--dir-path-out')
    error = error + this%switches%set(key = num_parts_key   , value = '--num_parts')
    error = error + this%switches%set(key = num_levels_key  , value = '--num_levels')
    error = error + this%switches%set(key = num_parts_per_level_key, value = '--num_parts_per_level')

    check(error==0)

    call this%switches_ab%init()
    error = error + this%switches_ab%set(key = dir_path_key    , value = '-d')
    error = error + this%switches_ab%set(key = prefix_key      , value = '-p')
    error = error + this%switches_ab%set(key = dir_path_out_key, value = '-o')
    error = error + this%switches_ab%set(key = num_parts_key   , value = '-n')
    error = error + this%switches_ab%set(key = num_levels_key  , value = '-l')
    error = error + this%switches_ab%set(key = num_parts_per_level_key   , value = '-npl')
    check(error==0)

    call this%helpers%init()
    error = error + this%helpers%set(key = dir_path_key    , value = 'Directory of the source files')
    error = error + this%helpers%set(key = prefix_key      , value = 'Name of the GiD files')
    error = error + this%helpers%set(key = dir_path_out_key, value = 'Output Directory')
    error = error + this%helpers%set(key = num_parts_key   , value = 'Number of parts of the mesh')
    error = error + this%helpers%set(key = num_levels_key  , value = 'Number of levels')
    error = error + this%helpers%set(key = num_parts_per_level_key   , value = 'Number of parts per level (array of fixed size 5)')
    check(error==0)

    call this%required%init()
    error = error + this%required%set(key = dir_path_key    , value = .false.)
    error = error + this%required%set(key = prefix_key      , value = .false.)
    error = error + this%required%set(key = dir_path_out_key, value = .false.)
    error = error + this%required%set(key = num_parts_key   , value = .false.)
    error = error + this%required%set(key = num_levels_key  , value = .false.)
    error = error + this%required%set(key = num_parts_per_level_key  , value = .false.)
    check(error==0)

  end subroutine partitioner_input_set_default

  !==================================================================================================
  subroutine partitioner_input_free(this)
    implicit none
    class(partitioner_input_t), intent(inout) :: this
    call this%list%free()
    call this%switches%free()
    call this%switches_ab%free()
    call this%required%free()
    call this%cli%free()
   end subroutine partitioner_input_free

  !==================================================================================================
  function partitioner_input_get_parameters(this)
    implicit none
    class(partitioner_input_t), target , intent(in) :: this
    type(ParameterList_t), pointer  :: partitioner_input_get_parameters
    partitioner_input_get_parameters => this%list
  end function partitioner_input_get_parameters

  !==================================================================================================
  !
  ! The following methods can be programmed in the library looping over the entries in, e.g. switch.  
  ! To do that we need to manage data types conversions automatically. Here I'm exploiting the knowledge
  ! of the data type of each entry. We could ask fpl...
  !
  !==================================================================================================
  subroutine partitioner_input_add_to_cli(this)
    implicit none
    class(partitioner_input_t) , intent(inout) :: this
    integer(ip)        :: error
    character(len=512) :: switch, switch_ab, help ! , cvalue
    logical            :: required
    integer(ip)        :: ivalue
    character(len=:), allocatable :: key, cvalue !, switch, switch_ab, help
    type(ParameterListIterator_t) :: Iterator

    error = 0
    Iterator = this%switches%GetIterator()
    do while (.not. Iterator%HasFinished())
       key = Iterator%GetKey()
       error = error + Iterator%Get(switch)
       error = error + this%switches_ab%get  (key = key , value = switch_ab)
       error = error + this%helpers%get      (key = key , value = help)
       error = error + this%required%get     (key = key , value = required)
       error = error + this%list%GetAsString (key = key , string = cvalue, separator=" ")

       if(this%list%GetDimensions(Key=Iterator%GetKey()) == 0) then 
          call this%cli%add(switch=trim(switch),switch_ab=trim(switch_ab), help=trim(help), &
            &               required=required,act='store',def=trim(cvalue),error=error)
       else if(this%list%GetDimensions(Key=Iterator%GetKey()) == 1) then 
          call this%cli%add(switch=trim(switch),switch_ab=trim(switch_ab), help=trim(help), &
            &               required=required,act='store',def=trim(cvalue),error=error,nargs='5')
       else
          write(*,*) 'Rank >1 arrays not supported by CLI'
          check(.false.)
       end if
          check(error==0)
       call Iterator%Next()
    enddo

    !! IO parameters
    !error = 0
    !error = error + this%list%get       (key = dir_path_key , value = cvalue)
    !error = error + this%switches%get   (key = dir_path_key , value = switch)
    !error = error + this%switches_ab%get(key = dir_path_key , value = switch_ab)
    !error = error + this%helpers%get    (key = dir_path_key , value = help)
    !error = error + this%required%get   (key = dir_path_key , value = required)
    !call this%cli%add(switch=trim(switch),switch_ab=trim(switch_ab), help=trim(help), &
    !     &            required=required,act='store',def=trim(cvalue),error=error)
    !check(error==0)

    !error = 0
    !error = error + this%list%get       (key = prefix_key , value = cvalue)
    !error = error + this%switches%get   (key = prefix_key , value = switch)
    !error = error + this%switches_ab%get(key = prefix_key , value = switch_ab)
    !error = error + this%helpers%get    (key = prefix_key , value = help)
    !error = error + this%required%get   (key = prefix_key , value = required)
    !check(error==0)
    !call this%cli%add(switch=trim(switch),switch_ab=trim(switch_ab), help=trim(help), &
    !     &            required=required,act='store',def=trim(cvalue),error=error)
    !check(error==0)

    !error = 0
    !error = error + this%list%get       (key = dir_path_out_key , value = cvalue)
    !error = error + this%switches%get   (key = dir_path_out_key , value = switch)
    !error = error + this%switches_ab%get(key = dir_path_out_key , value = switch_ab)
    !error = error + this%helpers%get    (key = dir_path_out_key , value = help)
    !error = error + this%required%get   (key = dir_path_out_key , value = required)
    !check(error==0)
    !call this%cli%add(switch=trim(switch),switch_ab=trim(switch_ab), help=trim(help), &
    !     &            required=required,act='store',def=trim(cvalue),error=error)
    !check(error==0)

    !error = 0
    !error = error + this%list%get       (key = num_parts_key , value = ivalue)
    !error = error + this%switches%get   (key = num_parts_key , value = switch)
    !error = error + this%switches_ab%get(key = num_parts_key , value = switch_ab)
    !error = error + this%helpers%get    (key = num_parts_key , value = help)
    !error = error + this%required%get   (key = num_parts_key , value = required)
    !!write(*,*) ivalue
    !write(cvalue,*) ivalue
    !check(error==0)
    !call this%cli%add(switch=trim(switch),switch_ab=trim(switch_ab), help=trim(help), &
    !     &            required=required,act='store',def=trim(cvalue),error=error)
    !check(error==0)

  end subroutine partitioner_input_add_to_cli

  subroutine partitioner_input_parse(this)
    implicit none
    class(partitioner_input_t), intent(inout) :: this
    integer(ip)    :: istat, error
    character(512) :: switch ! , cvalue
    integer(ip)    :: ivalue

    character(len=:), allocatable :: key
    type(ParameterListIterator_t) :: Iterator
    class(*), pointer :: val0
    class(*), pointer :: val1(:)
    
    call this%cli%parse(error=error); check(error==0)

    error = 0
    Iterator = this%switches%GetIterator()
    do while (.not. Iterator%HasFinished())
       key = Iterator%GetKey()
       error = error + Iterator%Get(switch)
       if (this%cli%is_passed(switch=switch)) then
          if(this%list%GetDimensions(key = key)==0) then
             error = error + this%list%GetPointer(key = key, value=val0)
             call this%cli%get(switch=switch, val=val0, error=error)
          else if(this%list%GetDimensions(key = key)==1) then
             error = error + this%list%GetPointer(key = key, value=val1)
             call this%cli%get(switch=switch, val=val1, error=error)
          end if
       end if
       check(error==0)
       call Iterator%Next()
    enddo

    !istat = this%switches%get(key = dir_path_key , value = switch)
    !check(istat==0)
    !if (this%cli%is_passed(switch=switch)) then
    !   call this%cli%get(switch=switch, val=cvalue, error=istat); check(istat==0)
    !   istat = this%list%set(key = dir_path_key, value=cvalue)
    !end if

    !istat = this%switches%get(key = prefix_key , value = switch)
    !check(istat==0)
    !if (this%cli%is_passed(switch=switch)) then
    !   call this%cli%get(switch=switch, val=cvalue, error=istat); check(istat==0)
    !   istat = this%list%set(key = prefix_key, value=cvalue)
    !end if

    !istat = this%switches%get(key = dir_path_out_key , value = switch)
    !check(istat==0)
    !if (this%cli%is_passed(switch=switch)) then
    !   call this%cli%get(switch=switch, val=cvalue, error=istat); check(istat==0)
    !   istat = this%list%set(key = dir_path_out_key, value=cvalue)
    !end if

    !istat = this%switches%get(key = num_parts_key , value = switch)
    !check(istat==0)
    !if (this%cli%is_passed(switch=switch)) then
    !   call this%cli%get(switch=switch, val=ivalue, error=istat); check(istat==0)
    !   istat = this%list%set(key = num_parts_key, value=ivalue)
    !end if

  end subroutine partitioner_input_parse  

end module partitioner_input_names 

!==================================================================================================
!==================================================================================================
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
  type(par_environment_t)  , allocatable :: env(:)
  type(mesh_t)             , allocatable :: lmesh(:)
  integer(ip) :: ipart, ienv

  call fempar_init()
  call input%create()
  parameters => input%get_parameters()

  ! Read and partition gmesh into lmesh
  call gmesh%read(parameters)
  call gmesh%write_file_for_postprocess(parameters)
  call gmesh%create_distribution (parameters, distr, env, lmesh)

  ! Write environments
  call par_environment_write_files             ( parameters, env )

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
