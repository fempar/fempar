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
module parameter_generator_names
  use types_names
  use flap, only : Command_Line_Interface
  use FPL
# include "debug.i90"
  implicit none
  private

  ! This type implements the coupling between FPL and the cli. From a user point
  ! of view it is only necessary to extend it implementing set_default, where the 
  ! parameters required by the user have to be registered with a default value
  ! and for those that could be read from the command line register the switches,
  ! abbreviated_switches, helpers and whether they are mandatory or not.
  type parameter_generator_t 
     private 
     type(Command_Line_Interface)  :: cli 
     type(ParameterList_t)         :: list
     type(ParameterList_t)         :: switches
     type(ParameterList_t)         :: switches_ab
     type(ParameterList_t)         :: helpers
     type(ParameterList_t)         :: required
   contains
     procedure, non_overridable    :: create          => parameter_generator_create
     procedure                     :: set_default     => parameter_generator_set_default
     procedure, non_overridable    :: add_to_cli      => parameter_generator_add_to_cli
     procedure, non_overridable    :: parse           => parameter_generator_parse
     procedure, non_overridable    :: free            => parameter_generator_free
     procedure, non_overridable    :: get_parameters  => parameter_generator_get_parameters 
     procedure, non_overridable    :: get_switches    => parameter_generator_get_switches   
     procedure, non_overridable    :: get_switches_ab => parameter_generator_get_switches_ab
     procedure, non_overridable    :: get_helpers     => parameter_generator_get_helpers    
     procedure, non_overridable    :: get_required    => parameter_generator_get_required   
  end type parameter_generator_t

  public :: parameter_generator_t

contains

  subroutine parameter_generator_create(this)
    implicit none
    class(parameter_generator_t), intent(inout) :: this
    call this%free()
     ! Initialize Command Line Interface
    call this%cli%init(progname    = 'part',                                                     &
         &        version     = '',                                                                 &
         &        authors     = '',                                                                 &
         &        license     = '',                                                                 &
         &        description =  'FEMPAR driver to part a GiD mesh.', &
         &        examples    = ['part -h  ', 'part -n  ' ])

    call this%list%init()
    call this%switches%init()
    call this%switches_ab%init()
    call this%helpers%init()
    call this%required%init()

    call this%set_default()
    call this%add_to_cli()
    call this%parse()
  end subroutine parameter_generator_create

  !==================================================================================================
  subroutine parameter_generator_set_default(this)
    implicit none
    class(parameter_generator_t), intent(inout) :: this
    integer(ip) :: error
    ! This is necessary in derived classes implemted in the user space (here we could access
    ! member variables directly.
    type(ParameterList_t), pointer :: list, switches, switches_ab, helpers, required

    list        => this%get_parameters()
    switches    => this%get_switches()
    switches_ab => this%get_switches_ab()
    helpers     => this%get_helpers()
    required    => this%get_required()

    error = list%set(key = dir_path_key       , value = '.')      ; check(error==0);
    error = list%set(key = prefix_key         , value = 'problem'); check(error==0);

    ! Not all of them need to be controlled from cli
    error = switches%set(key = dir_path_key   , value = '--dir-path'); check(error==0);
    error = switches%set(key = prefix_key     , value = '--prefix')  ; check(error==0);

    error = switches_ab%set(key = dir_path_key, value = '-d'); check(error==0);
    error = switches_ab%set(key = prefix_key  , value = '-p'); check(error==0);

    error = helpers%set(key = dir_path_key    , value = 'Directory of the source files'); check(error==0);
    error = helpers%set(key = prefix_key      , value = 'Name of the project')          ; check(error==0);

    error = required%set(key = dir_path_key   , value = .false.); check(error==0);
    error = required%set(key = prefix_key     , value = .false.); check(error==0);
  end subroutine parameter_generator_set_default

  !==================================================================================================
  function parameter_generator_get_parameters(this)
    implicit none
    class(parameter_generator_t), target , intent(in) :: this
    type(ParameterList_t)  , pointer    :: parameter_generator_get_parameters
    parameter_generator_get_parameters => this%list
  end function parameter_generator_get_parameters

  !==================================================================================================
  function parameter_generator_get_switches(this)
    implicit none
    class(parameter_generator_t), target , intent(in) :: this
    type(ParameterList_t)  , pointer    :: parameter_generator_get_switches
    parameter_generator_get_switches => this%switches
  end function parameter_generator_get_switches
  !==================================================================================================
  function parameter_generator_get_switches_ab(this)
    implicit none
    class(parameter_generator_t), target , intent(in) :: this
    type(ParameterList_t)  , pointer    :: parameter_generator_get_switches_ab
    parameter_generator_get_switches_ab => this%switches_ab
  end function parameter_generator_get_switches_ab
  !==================================================================================================
  function parameter_generator_get_helpers(this)
    implicit none
    class(parameter_generator_t), target , intent(in) :: this
    type(ParameterList_t)  , pointer    :: parameter_generator_get_helpers
    parameter_generator_get_helpers => this%helpers
  end function parameter_generator_get_helpers

  !==================================================================================================
  function parameter_generator_get_required(this)
    implicit none
    class(parameter_generator_t), target , intent(in) :: this
    type(ParameterList_t)  , pointer    :: parameter_generator_get_required
    parameter_generator_get_required => this%required
  end function parameter_generator_get_required

  !==================================================================================================
  subroutine parameter_generator_free(this)
    implicit none
    class(parameter_generator_t), intent(inout) :: this
    call this%list%free()
    call this%switches%free()
    call this%switches_ab%free()
    call this%required%free()
    call this%cli%free()
   end subroutine parameter_generator_free
  !==================================================================================================
  subroutine parameter_generator_add_to_cli(this)
    implicit none
    class(parameter_generator_t) , intent(inout) :: this
    integer(ip)                   :: error
    character(len=:), allocatable :: switch, switch_ab, help ! , cvalue
    logical                       :: required
    integer(ip)                   :: ivalue
    character(len=:), allocatable :: key, cvalue !, switch, switch_ab, help
    type(ParameterListIterator_t) :: Iterator

    error = 0
    Iterator = this%switches%GetIterator()
    do while (.not. Iterator%HasFinished())
       key = Iterator%GetKey()
       error = error + Iterator%GetAsString (switch)
       error = error + this%switches_ab%GetAsString (key = key , String = switch_ab)
       error = error + this%helpers%GetAsString     (key = key , String = help)
       error = error + this%required%Get            (key = key , value = required)
       error = error + this%list%GetAsString        (key = key , string = cvalue, separator=" ")

       if(this%list%GetDimensions(Key=Iterator%GetKey()) == 0) then 
          call this%cli%add(switch=trim(switch),switch_ab=trim(switch_ab), help=trim(help), &
            &               required=required,act='store',def=trim(cvalue),error=error)
       else if(this%list%GetDimensions(Key=Iterator%GetKey()) == 1) then 
          call this%cli%add(switch=trim(switch),switch_ab=trim(switch_ab), help=trim(help), &
            &               required=required,act='store',def=trim(cvalue),error=error,nargs='+')
       else
          write(*,*) 'Rank >1 arrays not supported by CLI'
          check(.false.)
       end if
          check(error==0)
       call Iterator%Next()
    enddo

  end subroutine parameter_generator_add_to_cli

  !==================================================================================================
  subroutine parameter_generator_parse(this)
    implicit none
    class(parameter_generator_t), intent(inout) :: this
    integer(ip)                :: istat, error
    character(len=str_cla_len) :: switch ! , cvalue
    integer(ip)                :: ivalue

    character(len=:), allocatable :: key
    type(ParameterListIterator_t) :: Iterator
    class(*), pointer :: val0
    class(*), pointer :: val1(:)
    integer(ip), allocatable :: val_ip(:)
    real(rp)   , allocatable :: val_rp(:)
    character(512)           :: val_ch
    
    call this%cli%parse(error=error); check(error==0)

    Iterator = this%switches%GetIterator()
    do while (.not. Iterator%HasFinished())
       key = Iterator%GetKey()
       error = Iterator%Get(switch); check(error==0)
       if (this%cli%is_passed(switch=switch)) then
          if(this%list%GetDimensions(key = key)==0) then
             error = this%list%GetPointer(key = key, value=val0); check(error==0)
             select type(val0)
             type is(character(*))
                 call this%cli%get(switch=switch, val=val_ch, error=error)
                 error = this%list%Set(key = key, value=trim(val_ch)); check(error==0)
             class default
                 call this%cli%get(switch=switch, val=val0, error=error)
             end select
          else if(this%list%GetDimensions(key = key)==1) then
             error = this%list%GetPointer(key = key, value=val1); check(error==0)
             select type(val1)
             type is(integer(ip))
                call this%cli%get_varying(switch=switch, val=val_ip, error=error); check(error==0)
                error = this%list%set(key = key, value = val_ip); check(error==0)
             type is(real(rp))
                call this%cli%get_varying(switch=switch, val=val_rp, error=error); check(error==0)
                error = this%list%set(key = key, value = val_rp); check(error==0)
             class default
                check(.false.)
             end select
          end if
       end if
       call Iterator%Next()
    enddo

  end subroutine parameter_generator_parse  

end module parameter_generator_names 
