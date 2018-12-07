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
module parameter_handler_names
  use types_names
  use flap, only : Command_Line_Interface
  use flap_utils_m
  use FPL
# include "debug.i90"
  implicit none
  private

  type, abstract :: parameter_handler_t 
  !------------------------------------------------------------------
  !< This type implements the coupling between FPL and the cli. From a user point
  !< of view it is only necessary to extend it implementing set_default(), where the 
  !< parameters required by the user have to be registered with a default value
  !< and for those that could be read from the command line register the switches,
  !< abbreviated_switches, helpers and whether they are mandatory or not.
  !------------------------------------------------------------------
     private 
     type(Command_Line_Interface)  :: cli 
     type(ParameterList_t)         :: values
     type(ParameterList_t)         :: switches
     type(ParameterList_t)         :: switches_ab
     type(ParameterList_t)         :: helpers
     type(ParameterList_t)         :: required
   contains
     procedure                                  :: create                   => parameter_handler_create
     procedure(define_parameters_interface), deferred :: define_parameters 
     procedure, non_overridable, private        :: assert_lists_consistency => parameter_handler_assert_lists_consistency
     procedure, non_overridable, private        :: add_to_cli               => parameter_handler_add_to_cli
     procedure, non_overridable, private        :: add_to_cli_group         => parameter_handler_add_to_cli_group
     procedure, non_overridable, private        :: parse                    => parameter_handler_parse
     procedure, non_overridable, private        :: parse_group              => parameter_handler_parse_group
     procedure                                  :: free                     => parameter_handler_free
     procedure, non_overridable                 :: get_values               => parameter_handler_get_values 
     procedure, non_overridable                 :: get_switches             => parameter_handler_get_switches   
     procedure, non_overridable                 :: get_switches_ab          => parameter_handler_get_switches_ab
     procedure, non_overridable                 :: get_helpers              => parameter_handler_get_helpers    
     procedure, non_overridable                 :: get_required             => parameter_handler_get_required   
     procedure, non_overridable, private        :: initialize_lists         => parameter_handler_initialize_lists  
  end type parameter_handler_t

  public :: parameter_handler_t
  public :: parameter_handler_free
  public :: count_tokens_cla_string_list 
  public :: process_tokens_cla_string_list
  
  
  abstract interface
    subroutine define_parameters_interface(this)
    !------------------------------------------------------------------
    !< Deferred binding that all subclasses are forced to implement.
    !< Subclasses have to (consistently) set-up the type(ParameterList_t)
    !< member variables above s.t. they can be transferred to the cli
    !< instance by means of the parse binding (see the code of 
    !< parameter_generator_assert_lists_consistency)
    !------------------------------------------------------------------
      import :: parameter_handler_t
      implicit none
      class(parameter_handler_t), intent(inout) :: this
    end subroutine define_parameters_interface
  end interface

contains

  subroutine parameter_handler_create(this, parse_cla, progname, version, help, description, &
                                        license, authors, examples, epilog, disable_hv, &
                                        usage_lun, error_lun, version_lun)
    implicit none
    class(parameter_handler_t), intent(inout) :: this
    logical,            optional, intent(in)    :: parse_cla         !< Parse command line arguments
    character(*),       optional, intent(in)    :: progname          !< Program name.
    character(*),       optional, intent(in)    :: version           !< Program version.
    character(*),       optional, intent(in)    :: help              !< Help message introducing the CLI usage.
    character(*),       optional, intent(in)    :: description       !< Detailed description message introducing the program.
    character(*),       optional, intent(in)    :: license           !< License description.
    character(*),       optional, intent(in)    :: authors           !< Authors list.
    character(*),       optional, intent(in)    :: examples(1:)      !< Examples of correct usage.
    character(*),       optional, intent(in)    :: epilog            !< Epilog message.
    logical,            optional, intent(in)    :: disable_hv        !< Disable automatic insert of 'help' and 'version' CLAs.
    integer(ip),        optional, intent(in)    :: usage_lun         !< Unit number to print usage/help
    integer(ip),        optional, intent(in)    :: version_lun       !< Unit number to print version/license info
    integer(ip),        optional, intent(in)    :: error_lun         !< Unit number to print error info
    
    logical :: parse_cla_
    parse_cla_ = .true.
    if ( present(parse_cla) ) parse_cla_ = parse_cla

    call this%free()

    if ( parse_cla_ ) then
      call this%cli%init(progname, version, help, description, &
                         license, authors, examples, epilog, disable_hv, &
                         usage_lun, error_lun, version_lun)
    end if
    
    call this%initialize_lists()
    call this%define_parameters()
    
#ifdef DEBUG
    call this%assert_lists_consistency()
#endif
    if ( parse_cla_ ) then 
      call this%add_to_cli()
      call this%parse()
    end if
  end subroutine parameter_handler_create

  !==================================================================================================
  function parameter_handler_get_values(this)
    implicit none
    class(parameter_handler_t), target , intent(in) :: this
    type(ParameterList_t)  , pointer    :: parameter_handler_get_values
    parameter_handler_get_values => this%values
  end function parameter_handler_get_values

  !==================================================================================================
  function parameter_handler_get_switches(this)
    implicit none
    class(parameter_handler_t), target , intent(in) :: this
    type(ParameterList_t)  , pointer    :: parameter_handler_get_switches
    parameter_handler_get_switches => this%switches
  end function parameter_handler_get_switches
  !==================================================================================================
  function parameter_handler_get_switches_ab(this)
    implicit none
    class(parameter_handler_t), target , intent(in) :: this
    type(ParameterList_t)  , pointer    :: parameter_handler_get_switches_ab
    parameter_handler_get_switches_ab => this%switches_ab
  end function parameter_handler_get_switches_ab
  !==================================================================================================
  function parameter_handler_get_helpers(this)
    implicit none
    class(parameter_handler_t), target , intent(in) :: this
    type(ParameterList_t)  , pointer    :: parameter_handler_get_helpers
    parameter_handler_get_helpers => this%helpers
  end function parameter_handler_get_helpers

  !==================================================================================================
  function parameter_handler_get_required(this)
    implicit none
    class(parameter_handler_t), target , intent(in) :: this
    type(ParameterList_t)  , pointer    :: parameter_handler_get_required
    parameter_handler_get_required => this%required
  end function parameter_handler_get_required

  !==================================================================================================
  subroutine parameter_handler_free(this)
    implicit none
    class(parameter_handler_t), intent(inout) :: this
    call this%values%free()
    call this%switches%free()
    call this%switches_ab%free()
    call this%helpers%free()
    call this%required%free()
    call this%helpers%free()
    call this%cli%free()
   end subroutine parameter_handler_free
   
   subroutine parameter_handler_assert_lists_consistency(this)
     implicit none
     class(parameter_handler_t), intent(in) :: this
     logical :: has_groups
     character(len=:), allocatable  :: key
     type(ParameterListIterator_t)  :: Iterator
     type(ParameterListIterator_t)  :: Sublist_Iterator
     type(ParameterList_t), pointer :: switches_sublist
     type(ParameterList_t), pointer :: switches_ab_sublist
     type(ParameterList_t), pointer :: helpers_sublist
     type(ParameterList_t), pointer :: required_sublist
     type(ParameterList_t), pointer :: values_sublist
     
!     assert(this%switches%length() == this%switches_ab%length())
     assert(this%switches%length() == this%helpers%length())
!     assert(this%switches%length() == this%required%length())
!     assert(this%switches%length() == this%values%length())
     
     Iterator = this%switches%GetIterator()
     do while (.not. Iterator%HasFinished())
          key = Iterator%GetKey()
          !assert(this%switches_ab%isPresent(key))
          assert(this%helpers%isPresent(key))
          !assert(this%required%isPresent(key))
          assert(this%values%isPresent(key))
          if (Iterator%isSubList()) then
            assert(Iterator%GetSublist(switches_sublist)==0)
            assert(this%switches_ab%GetSubList(key,switches_ab_sublist)==0)
            assert(this%helpers%GetSubList(key,helpers_sublist)==0)
            assert(this%required%GetSubList(key,required_sublist)==0)
            assert(this%values%GetSubList(key,values_sublist)==0)
            assert(switches_sublist%length() == switches_ab_sublist%length())
            assert(switches_sublist%length() == helpers_sublist%length())
            assert(switches_sublist%length() == required_sublist%length())
            assert(switches_sublist%length() == values_sublist%length())
            sublist_iterator = switches_sublist%GetIterator()
            do while (.not. sublist_iterator%HasFinished())
                key = sublist_iterator%GetKey()       
                assert(switches_ab_sublist%isPresent(key))
                assert(helpers_sublist%isPresent(key))
                assert(required_sublist%isPresent(key))
                assert(values_sublist%isPresent(key))
                assert(.not. sublist_iterator%isSubList())
                assert(.not. switches_ab_sublist%isSubList(key))
                assert(.not. helpers_sublist%isSubList(key))
                assert(.not. required_sublist%isSubList(key))
                assert(.not. values_sublist%isSubList(key))
                call sublist_iterator%Next()  
            end do    
          end if
          call Iterator%Next()
     end do
   end subroutine parameter_handler_assert_lists_consistency
   
   
  !==================================================================================================
  subroutine parameter_handler_add_to_cli(this)
    implicit none
    class(parameter_handler_t) , intent(inout) :: this
    integer(ip)                   :: error
    character(len=:), allocatable :: switch, switch_ab, help ! , cvalue
    logical                       :: required
    integer(ip)                   :: ivalue
    character(len=:), allocatable :: key, cvalue !, switch, switch_ab, help
    type(ParameterListIterator_t) :: Iterator

    type(ParameterList_t), pointer :: switches_sublist
    type(ParameterList_t), pointer :: switches_ab_sublist
    type(ParameterList_t), pointer :: helpers_sublist
    type(ParameterList_t), pointer :: required_sublist
    type(ParameterList_t), pointer :: values_sublist
    
    error = 0
    call this%add_to_cli_group(this%switches,this%switches_ab,this%helpers,this%required,this%values)
    Iterator = this%switches%GetIterator()
    do while (.not. Iterator%HasFinished())
       if(.not. Iterator%isSublist()) then
          call Iterator%Next()
          cycle
       endif
       key = Iterator%GetKey()
       error = Iterator%GetSublist(switches_sublist); assert(error==0);
       error = this%switches_ab%GetSubList(key,switches_ab_sublist); assert(error==0);
       error = this%helpers%GetSubList(key,helpers_sublist); assert(error==0);
       error = this%required%GetSubList(key,required_sublist); assert(error==0);
       error = this%values%GetSubList(key,values_sublist); assert(error==0);
       call this%add_to_cli_group(switches_sublist, &
                                  switches_ab_sublist,&
                                  helpers_sublist, &
                                  required_sublist, &
                                  values_sublist, &
                                  key)
       call Iterator%Next()
    enddo

  end subroutine parameter_handler_add_to_cli
  
  !==================================================================================================
  subroutine parameter_handler_add_to_cli_group(this,switches,switches_ab,helpers,required,values,group)
    implicit none
    class(parameter_handler_t) , intent(inout) :: this
    type(parameterlist_t)       , intent(in)   :: switches
    type(parameterlist_t)       , intent(in)   :: switches_ab
    type(parameterlist_t)       , intent(in)   :: helpers
    type(parameterlist_t)        , intent(in)  :: required
    type(parameterlist_t)       , intent(in)   :: values
    character(*), optional       , intent(in)  :: group
    
    integer(ip)                   :: error
    character(len=:), allocatable :: switch, switch_ab, help
    logical                       :: is_required
    integer(ip)                   :: ivalue
    character(len=:), allocatable :: key, cvalue
    type(ParameterListIterator_t) :: Iterator

    error = 0
    Iterator = switches%GetIterator()
    do while (.not. Iterator%HasFinished())
       if(Iterator%isSublist()) then
          call Iterator%Next()
          cycle
       endif
       key = Iterator%GetKey()
       error = error + Iterator%GetAsString   (switch)
       if ( switches_ab%isPresent(key = key) ) then
         error = error + switches_ab%GetAsString(key = key , String = switch_ab)
       end if
       error = error + helpers%GetAsString    (key = key , String = help)
       !error = error + required%Get           (key = key , value = is_required)
       is_required = .false.
       error = error + values%GetAsString     (key = key , string = cvalue, separator=" ")
       if(values%GetDimensions(Key=Iterator%GetKey()) == 0) then 
         if ( switches_ab%isPresent(key = key) ) then
           call this%cli%add(group=group,switch=switch,switch_ab=switch_ab, help=help, &
             &               required=is_required,act='store',def=cvalue,error=error)
         else 
           call this%cli%add(group=group,switch=switch,help=help, &
             &               required=is_required,act='store',def=cvalue,error=error)
         end if  
       else if(values%GetDimensions(Key=Iterator%GetKey()) == 1) then 
         if ( switches_ab%isPresent(key = key) ) then
           call this%cli%add(group=group,switch=switch,switch_ab=switch_ab, help=help, &
             &               required=is_required,act='store',def=cvalue,error=error,nargs='+')
         else 
           call this%cli%add(group=group,switch=switch, help=help, &
             &               required=is_required,act='store',def=cvalue,error=error,nargs='+')
         end if
       else
          write(*,*) 'Rank >1 arrays not supported by CLI'
          assert(.false.)
       end if
       assert(error==0)
       call Iterator%Next()
    end do
  end subroutine parameter_handler_add_to_cli_group

  !==================================================================================================
  subroutine parameter_handler_parse(this)
    implicit none
    class(parameter_handler_t), intent(inout) :: this
    type(parameterlist_t), pointer       :: switches_sublist
    type(parameterlist_t), pointer       :: values_sublist
    integer(ip)                :: istat, error
    character(len=str_cla_len) :: switch ! , cvalue
    integer(ip)                :: ivalue
    character(len=:), allocatable :: key
    type(ParameterListIterator_t) :: Iterator
    
    call this%cli%parse(error=error); assert(error==0)

    call this%parse_group(this%switches, this%values)
    Iterator = this%switches%GetIterator()
    do while (.not. Iterator%HasFinished())
        if(.not. Iterator%isSublist()) then
           call Iterator%Next()
           cycle
        endif
        key = Iterator%GetKey()
        if(this%cli%run_command(key)) then
           error = Iterator%GetSublist(switches_sublist); assert(error==0);
           error = this%values%GetSubList(key,values_sublist); assert(error==0);
           call this%parse_group(switches_sublist, values_sublist, key)
        endif
          call Iterator%Next()
    enddo

  end subroutine parameter_handler_parse  

  
  subroutine parameter_handler_parse_group(this,switches,values,group)
    implicit none
    class(parameter_handler_t) , intent(inout) :: this
    type(parameterlist_t)       , intent(inout)    :: switches
    type(parameterlist_t)       , intent(inout)    :: values
    character(*), optional, intent(in)    :: group
    integer(ip)                :: istat, error
    character(len=str_cla_len) :: switch ! , cvalue
    integer(ip)                :: ivalue

    character(len=:), allocatable :: key
    type(ParameterListIterator_t) :: Iterator
    class(*), pointer :: val0
    class(*), pointer :: val1(:)
    integer(ip), allocatable :: val_ip(:)
    real(rp)   , allocatable :: val_rp(:)
    character(str_cla_len)   :: val_ch

    Iterator = switches%GetIterator()
    do while (.not. Iterator%HasFinished())
       if(Iterator%isSublist()) then
          call Iterator%Next()
          cycle
       endif
       key = Iterator%GetKey()
       error = Iterator%Get(switch); assert(error==0)
       if (this%cli%is_passed(group=group,switch=switch)) then
          if(values%GetDimensions(key = key)==0) then
             error = values%GetPointer(key = key, value=val0); assert(error==0)
             select type(val0)
             type is(character(*))
                 call this%cli%get(group=group,switch=switch, val=val_ch, error=error)
                 error = values%Set(key = key, value=val_ch); assert(error==0)
             class default
                call this%cli%get(group=group,switch=switch, val=val0, error=error)
             end select
          else if(values%GetDimensions(key = key)==1) then
             error = values%GetPointer(key = key, value=val1); assert(error==0)
             select type(val1)
             type is(integer(ip))
                call this%cli%get_varying(group=group,switch=switch, val=val_ip, error=error); assert(error==0)
                error = values%set(key = key, value = val_ip); assert(error==0)
             type is(real(rp))
                call this%cli%get_varying(group=group,switch=switch, val=val_rp, error=error); assert(error==0)
                error = values%set(key = key, value = val_rp); assert(error==0)
             class default
                assert(.false.)
             end select
          end if
       end if
       call Iterator%Next()
    enddo

  end subroutine parameter_handler_parse_group

  !==================================================================================================
  subroutine parameter_handler_initialize_lists(this)
    implicit none
    class(parameter_handler_t), intent(inout) :: this
    call this%values%init()
    call this%switches%init()
    call this%switches_ab%init()
    call this%helpers%init()
    call this%required%init()
  end subroutine parameter_handler_initialize_lists
  
  function count_tokens_cla_string_list (string_list) 
    implicit none
    character(*)       , intent(in) :: string_list
    integer(ip) :: count_tokens_cla_string_list
    character(len=len(string_list)), allocatable :: vals(:) !< String array of values based on buffer value.
    call tokenize(strin=string_list, delimiter=' ', toks=vals, Nt=count_tokens_cla_string_list)
    deallocate(vals)
  end function count_tokens_cla_string_list
  
  subroutine process_tokens_cla_string_list ( string_list, values ) 
    implicit none 
    character(*)       , intent(in)    :: string_list
    character(*)       , intent(inout) :: values(:)
    character(len=len(string_list)), allocatable :: vals(:) !< String array of values based on buffer value.
    integer(ip) :: Nt, i
    call tokenize(strin=string_list, delimiter=' ', toks=vals, Nt=Nt)
    assert ( size(values) == Nt) 
    do i=1, size(vals) 
      values(i)=vals(i)
    end do 
    deallocate(vals)
  end subroutine process_tokens_cla_string_list
  

end module parameter_handler_names 
