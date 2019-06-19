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
module fempar_parameter_handler_names
    use types_names
    use iterative_linear_solver_parameters_names
    use direct_solver_parameters_names
    use nonlinear_solver_parameters_names
    use environment_parameters_names
    use mesh_parameters_names
    use mesh_partitioner_parameters_names
    use metis_names
    use triangulation_parameters_names
    use uniform_hex_mesh_generator_parameters_names
    use p4est_triangulation_parameters_names
    use fe_space_parameters_names
    use reference_fe_names
    use output_handler_parameters_names
    use FPL
    use flap, only : Command_Line_Interface
# include "debug.i90"
    implicit none
    private

    character(*), parameter :: FLAP_HELP_MESSAGE_TABULATOR = "    "
    character(*), parameter :: BRK_LINE = NEW_LINE(FLAP_HELP_MESSAGE_TABULATOR) // FLAP_HELP_MESSAGE_TABULATOR

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
        type(ParameterList_t)         :: choices
    contains
        procedure, non_overridable, private        :: add0D                    => parameter_handler_add0D
        procedure, non_overridable, private        :: add1D                    => parameter_handler_add1D
        procedure, non_overridable, private        :: update0D                 => parameter_handler_update0D
        procedure, non_overridable, private        :: update1D                 => parameter_handler_update1D
        procedure, non_overridable                 :: add_to_cli               => parameter_handler_add_to_cli
        procedure, non_overridable, private        :: add_to_cli_group         => parameter_handler_add_to_cli_group
        procedure, non_overridable                 :: init_cli                 => parameter_handler_init_cli
        procedure, non_overridable                 :: parse                    => parameter_handler_parse
        procedure, non_overridable, private        :: parse_group              => parameter_handler_parse_group
        procedure                                  :: free                     => parameter_handler_free
        procedure, non_overridable                 :: get_values               => parameter_handler_get_values 
        procedure, non_overridable, private        :: get_switches             => parameter_handler_get_switches   
        procedure, non_overridable, private        :: get_switches_ab          => parameter_handler_get_switches_ab
        procedure, non_overridable, private        :: get_helpers              => parameter_handler_get_helpers    
        procedure, non_overridable, private        :: get_required             => parameter_handler_get_required 
        procedure, non_overridable, private        :: get_choices              => parameter_handler_get_choices 
        procedure, non_overridable                 :: initialize_lists         => parameter_handler_initialize_lists  
        generic,   public                          :: add                      => add0D, add1D
        generic,   public                          :: update                   => update0D, update1D 
    end type parameter_handler_t
 
    type, extends(parameter_handler_t) :: fempar_parameter_handler_t 
    private 
        procedure(define_user_parameters), nopass, pointer :: define_user_parameters => NULL()
    contains
        procedure, private, non_overridable :: fph_Get_0D
        procedure, private, non_overridable :: fph_Get_1D
        procedure, public                   :: getasstring              => fph_GetAsString
        
        procedure, private, non_overridable :: define_fempar_parameters => fph_define_fempar_parameters
        procedure, private, non_overridable :: fph_mesh_define_parameters
        procedure, private, non_overridable :: fph_mesh_partitioner_define_parameters
        procedure, private, non_overridable :: fph_static_triang_define_parameters
        procedure, private, non_overridable :: fph_struct_hex_mesh_generator_define_parameters
        procedure, private, non_overridable :: fph_p4est_triang_define_parameters
        procedure, private, non_overridable :: fph_fes_define_parameters
        procedure, private, non_overridable :: fph_coarse_fe_handler_define_parameters
        procedure, private, non_overridable :: fph_ils_define_parameters
        procedure, private, non_overridable :: fph_nls_define_parameters
        procedure, private, non_overridable :: fph_dls_define_parameters
        procedure, private, non_overridable :: fph_output_handler_define_parameters
        
        procedure                           :: process_parameters       => fph_process_parameters
        procedure                           :: free                     => fph_free
        procedure                           :: get_dir_path             => fph_get_dir_path
        procedure                           :: get_prefix               => fph_get_prefix
        generic,   public                   :: get                      => fph_Get_0D, fph_Get_1D
    end type fempar_parameter_handler_t
  
    interface
        subroutine define_user_parameters()
        end subroutine define_user_parameters
    end interface

    type(fempar_parameter_handler_t), save :: parameter_handler

    public :: parameter_handler, define_user_parameters

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PARAMETER_HANDLER_T PROCEDURES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function parameter_handler_get_values(this) result(values)
    !------------------------------------------------------------------
    !< Return pointer to values dictionary
    !------------------------------------------------------------------
        implicit none
        class(parameter_handler_t), target, intent(in) :: this
        type(ParameterList_t),      pointer            :: values
    !------------------------------------------------------------------
        values => this%values
    end function parameter_handler_get_values


    function parameter_handler_get_switches(this) result(switches)
    !------------------------------------------------------------------
    !< Return pointer to switches dictionary
    !------------------------------------------------------------------
        implicit none
        class(parameter_handler_t), target, intent(in) :: this
        type(ParameterList_t),      pointer            :: switches
    !------------------------------------------------------------------
        switches => this%switches
    end function parameter_handler_get_switches


    function parameter_handler_get_switches_ab(this) result(switches_ab)
    !------------------------------------------------------------------
    !< Return pointer to switches_ab dictionary
    !------------------------------------------------------------------
        implicit none
        class(parameter_handler_t), target, intent(in) :: this
        type(ParameterList_t),      pointer            :: switches_ab
    !------------------------------------------------------------------
        switches_ab => this%switches_ab
    end function parameter_handler_get_switches_ab


    function parameter_handler_get_helpers(this) result(helpers)
    !------------------------------------------------------------------
    !< Return pointer to helpers dictionary
    !------------------------------------------------------------------
        implicit none
        class(parameter_handler_t), target, intent(in) :: this
        type(ParameterList_t),      pointer            :: helpers
    !------------------------------------------------------------------
        helpers => this%helpers
    end function parameter_handler_get_helpers


    function parameter_handler_get_required(this) result(required)
    !------------------------------------------------------------------
    !< Return pointer to required dictionary
    !------------------------------------------------------------------
        implicit none
        class(parameter_handler_t), target, intent(in) :: this
        type(ParameterList_t),      pointer            :: required
    !------------------------------------------------------------------
        required => this%required
    end function parameter_handler_get_required


    function parameter_handler_get_choices(this) result(choices)
    !------------------------------------------------------------------
    !< Return pointer to choices dictionary
    !------------------------------------------------------------------
        implicit none
        class(parameter_handler_t), target, intent(in) :: this
        type(ParameterList_t),      pointer            :: choices
    !------------------------------------------------------------------
        choices => this%choices
    end function parameter_handler_get_choices


    subroutine parameter_handler_free(this)
    !------------------------------------------------------------------
    !< Free the parameter_handler_t
    !------------------------------------------------------------------
        implicit none
        class(parameter_handler_t), intent(inout) :: this
    !------------------------------------------------------------------
        call this%values%free()
        call this%switches%free()
        call this%switches_ab%free()
        call this%helpers%free()
        call this%required%free()
        call this%helpers%free()
        call this%choices%free()
        call this%cli%free()
    end subroutine parameter_handler_free


    subroutine parameter_handler_add0D(this, key, switch, value, help, switch_ab, required, choices, group)
    !------------------------------------------------------------------
    !< Add a new argument to command line interface
    !------------------------------------------------------------------
        implicit none
        class(parameter_handler_t), intent(inout) :: this
        character(len=*),           intent(in)    :: key
        character(len=*),           intent(in)    :: switch
        class(*),                   intent(in)    :: value 
        character(len=*),           intent(in)    :: help
        character(len=*), optional, intent(in)    :: switch_ab
        logical,          optional, intent(in)    :: required
        character(len=*), optional, intent(in)    :: choices
        character(len=*), optional, intent(in)    :: group
        type(ParameterList_t),      pointer       :: switches
        type(ParameterList_t),      pointer       :: values
        type(ParameterList_t),      pointer       :: helpers
        type(ParameterList_t),      pointer       :: switches_ab
        type(ParameterList_t),      pointer       :: requires
        type(ParameterList_t),      pointer       :: choice
        integer                                   :: error
    !------------------------------------------------------------------
        if(present(group)) then
            if(this%switches%isSublist(key=group)) then
                error = this%switches%getSublist(key=group, Sublist=switches); assert(error == 0)
            else
                switches => this%switches%NewSublist(key=group)
            endif
            if(this%values%isSublist(key=group)) then
                error = this%values%getSublist(key=group, Sublist=values); assert(error == 0)
            else
                values => this%values%NewSublist(key=group)
            endif
            if(this%helpers%isSublist(key=group)) then
                error = this%helpers%getSublist(key=group, Sublist=helpers); assert(error == 0)
            else
                helpers => this%helpers%NewSublist(key=group)
            endif
            if(present(switch_ab) .and. this%switches_ab%isSublist(key=group)) then
                error = this%switches_ab%getSublist(key=group, Sublist=switches_ab); assert(error == 0)
            else
                switches_ab => this%switches_Ab%NewSublist(key=group)
            endif
            if(present(required) .and. this%required%isSublist(key=group)) then
                error = this%required%getSublist(key=group, Sublist=requires); assert(error == 0)
            else
                requires => this%required%NewSublist(key=group)
            endif
            if(present(choices) .and. this%choices%isSublist(key=group)) then
                error = this%choices%getSublist(key=group, Sublist=choice); assert(error == 0)
            else
                choice => this%choices%NewSublist(key=group)
            endif
        else
            switches    => this%get_switches()
            values      => this%get_values()
            helpers     => this%get_helpers()
            if(present(switch_ab)) switches_ab => this%get_switches_ab()
            if(present(required))  requires    => this%get_required()
            if(present(choices))   choice      => this%get_choices()
        endif
        massert( .not. switches%isPresent(key = key), 'Duplicated Parameter: key "'//key//'" already defined in parameter handler' )
        error = switches%Set(key=key, value=switch); assert(error==0)
        error = helpers%Set(key=key,  value=help);   assert(error==0)
        error = values%Set(key=key,   value=value);  assert(error==0)
        if(present(switch_ab)) then
            error = switches_ab%Set(key=key, value=switch_ab); assert(error==0)
        endif
        if(present(required)) then
            error = requires%Set(key=key, value=required); assert(error==0)
        endif
        if(present(choices)) then
            error = choice%Set(key=key, value=choices); assert(error==0)
        endif
    end subroutine parameter_handler_add0D


    subroutine parameter_handler_add1D(this, key, switch, value, help, switch_ab, required, choices, group)
    !------------------------------------------------------------------
    !< Add a new argument to command line interface
    !------------------------------------------------------------------
        implicit none
        class(parameter_handler_t), intent(inout) :: this
        character(len=*),           intent(in)    :: key
        character(len=*),           intent(in)    :: switch
        class(*),                   intent(in)    :: value(:) 
        character(len=*),           intent(in)    :: help
        character(len=*), optional, intent(in)    :: switch_ab
        logical,          optional, intent(in)    :: required
        character(len=*), optional, intent(in)    :: choices
        character(len=*), optional, intent(in)    :: group
        type(ParameterList_t),      pointer       :: switches
        type(ParameterList_t),      pointer       :: values
        type(ParameterList_t),      pointer       :: helpers
        type(ParameterList_t),      pointer       :: switches_ab
        type(ParameterList_t),      pointer       :: requires
        type(ParameterList_t),      pointer       :: choice
        integer                                   :: error
    !------------------------------------------------------------------
        if(present(group)) then
            if(this%switches%isSublist(key=group)) then
                error = this%switches%getSublist(key=group, Sublist=switches); assert(error == 0)
            else
                switches => this%switches%NewSublist(key=group)
            endif
            if(this%values%isSublist(key=group)) then
                error = this%values%getSublist(key=group, Sublist=values);assert(error == 0)
            else
                values => this%values%NewSublist(key=group)
            endif
            if(this%helpers%isSublist(key=group)) then
                error = this%helpers%getSublist(key=group, Sublist=helpers);assert(error == 0)
            else
                helpers => this%helpers%NewSublist(key=group)
            endif
            if(present(switch_ab) .and. this%switches_ab%isSublist(key=group)) then
                error = this%switches_ab%getSublist(key=group, Sublist=switches_ab);assert(error == 0)
            else
                switches_ab => this%switches_Ab%NewSublist(key=group)
            endif
            if(present(required) .and. this%required%isSublist(key=group)) then
                error = this%required%getSublist(key=group, Sublist=requires);assert(error == 0)
            else
                requires => this%required%NewSublist(key=group)
            endif
            if(present(choices)  .and. this%choices%isSublist(key=group)) then
                error = this%choices%getSublist(key=group, Sublist=choice); assert(error == 0)
            else
                choice => this%choices%NewSublist(key=group)
            endif
        else
            switches    => this%get_switches()
            values      => this%get_values()
            helpers     => this%get_helpers()
            if(present(switch_ab)) switches_ab => this%get_switches_ab()
            if(present(required))  requires    => this%get_required()
            if(present(choices))   choice      => this%get_choices()
        endif
        massert( .not. switches%isPresent(key = key), 'Duplicated Parameter: key "'//key//'" already defined in parameter handler' )
        error = switches%Set(key=key, value=switch); assert(error==0)
        error = helpers%Set(key=key,  value=help);   assert(error==0)
        error = values%Set(key=key,   value=value);  assert(error==0)
        if(present(switch_ab)) then
            error = switches_ab%Set(key=key, value=switch_ab); assert(error==0)
        endif
        if(present(required)) then
            error = requires%Set(key=key, value=required); assert(error==0)
        endif
        if(present(choices)) then
            error = choice%Set(key=key, value=choices); assert(error==0)
        endif
    end subroutine parameter_handler_add1D



    subroutine parameter_handler_update0D(this, key, value, group)
    !------------------------------------------------------------------
    !< Update a scalar value given an argument key
    !------------------------------------------------------------------
        implicit none
        class(parameter_handler_t), intent(inout) :: this
        character(len=*),           intent(in)    :: key
        class(*),                   intent(in)    :: value 
        character(len=*), optional, intent(in)    :: group
        type(ParameterList_t),      pointer       :: values
        integer                                   :: error
    !------------------------------------------------------------------
        if(present(group)) then
            if(this%values%isSublist(key=group)) then
                error = this%values%getSublist(key=group, Sublist=values); assert(error == 0)
            else
                values => this%values%NewSublist(key=group)
            endif
        else
            values => this%get_values()
        endif
        error = values%Set(key=key,   value=value);  assert(error==0)
    end subroutine parameter_handler_update0D


    subroutine parameter_handler_update1D(this, key, value, group)
    !------------------------------------------------------------------
    !< Update a vector value given an argument key
    !------------------------------------------------------------------
        implicit none
        class(parameter_handler_t), intent(inout) :: this
        character(len=*),           intent(in)    :: key
        class(*),                   intent(in)    :: value(:) 
        character(len=*), optional, intent(in)    :: group
        type(ParameterList_t),      pointer       :: values
        integer                                   :: error
    !------------------------------------------------------------------
        if(present(group)) then
            if(this%values%isSublist(key=group)) then
                error = this%values%getSublist(key=group, Sublist=values);assert(error == 0)
            else
                values => this%values%NewSublist(key=group)
            endif
        else
            values => this%get_values()
        endif
        error = values%Set(key=key,   value=value);  assert(error==0)
    end subroutine parameter_handler_update1D

  
   
    subroutine parameter_handler_add_to_cli(this)
    !------------------------------------------------------------------
    !< Add to command line interface procedure
    !------------------------------------------------------------------
        implicit none
        class(parameter_handler_t) , intent(inout) :: this
        integer(ip)                                :: error
        character(len=:), allocatable              :: switch, switch_ab, help ! , cvalue
        logical                                    :: required
        integer(ip)                                :: ivalue
        character(len=:), allocatable              :: key, cvalue !, switch, switch_ab, help
        type(ParameterListIterator_t)              :: Iterator
        type(ParameterList_t), pointer             :: switches_sublist
        type(ParameterList_t), pointer             :: switches_ab_sublist
        type(ParameterList_t), pointer             :: helpers_sublist
        type(ParameterList_t), pointer             :: required_sublist
        type(ParameterList_t), pointer             :: values_sublist
        type(ParameterList_t), pointer             :: choices_sublist
    !------------------------------------------------------------------
        error = 0
        call this%add_to_cli_group(this%switches,this%switches_ab,this%helpers,this%required,this%values,this%choices)
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
            error = this%choices%GetSubList(key,choices_sublist); assert(error==0);
            call this%add_to_cli_group(switches_sublist,    &
                                       switches_ab_sublist, &
                                       helpers_sublist,     &
                                       required_sublist,    &
                                       values_sublist,      &
                                       choices_sublist,     &
                                       key)
            call Iterator%Next()
        enddo
    end subroutine parameter_handler_add_to_cli


    subroutine parameter_handler_add_to_cli_group(this,switches,switches_ab,helpers,required,values,choices,group)
    !------------------------------------------------------------------
    !< Add to command line interface group procedure
    !------------------------------------------------------------------
        implicit none
        class(parameter_handler_t), intent(inout) :: this
        type(parameterlist_t),      intent(in)    :: switches
        type(parameterlist_t),      intent(in)    :: switches_ab
        type(parameterlist_t),      intent(in)    :: helpers
        type(parameterlist_t),      intent(in)    :: required
        type(parameterlist_t),      intent(in)    :: values
        type(parameterlist_t),      intent(in)    :: choices
        character(*), optional,     intent(in)    :: group
        integer(ip)                               :: error
        character(len=:), allocatable             :: switch, switch_ab, help
        logical                                   :: is_required
        integer(ip)                               :: ivalue
        character(len=:), allocatable             :: key, cvalue, choice
        type(ParameterListIterator_t)             :: Iterator
    !------------------------------------------------------------------
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
            if ( choices%isPresent(key = key) ) then
                error = error + choices%GetAsString(key = key , String = choice)
            end if
            error = error + helpers%GetAsString    (key = key , String = help)
            !error = error + required%Get           (key = key , value = is_required)
            is_required = .false.
            error = error + values%GetAsString     (key = key , string = cvalue, separator=" ")
            if(values%GetDimensions(Key=Iterator%GetKey()) == 0) then 
                if ( switches_ab%isPresent(key = key) ) then
                    if ( choices%isPresent(key = key) ) then
                        call this%cli%add(group=group,switch=switch,switch_ab=switch_ab, help=help, &
                                     required=is_required,choices=choice,act='store',def=cvalue,error=error)
                    else
                        call this%cli%add(group=group,switch=switch,switch_ab=switch_ab, help=help, &
                                     required=is_required,act='store',def=cvalue,error=error)
                    end if
                else 
                    if ( choices%isPresent(key = key) ) then
                        call this%cli%add(group=group,switch=switch,help=help, &
                                     required=is_required,choices=choice,act='store',def=cvalue,error=error)
                    else
                        call this%cli%add(group=group,switch=switch,help=help, &
                                     required=is_required,act='store',def=cvalue,error=error)
                    endif
                end if  
            else if(values%GetDimensions(Key=Iterator%GetKey()) == 1) then 
                if ( switches_ab%isPresent(key = key) ) then
                    if ( choices%isPresent(key = key) ) then
                        call this%cli%add(group=group,switch=switch,switch_ab=switch_ab, help=help, &
                                     required=is_required,choices=choice,act='store',def=cvalue,error=error,nargs='+')
                    else
                        call this%cli%add(group=group,switch=switch,switch_ab=switch_ab, help=help, &
                                     required=is_required,act='store',def=cvalue,error=error,nargs='+')
                    endif
                else 
                    if ( choices%isPresent(key = key) ) then
                        call this%cli%add(group=group,switch=switch, help=help, &
                                     required=is_required,choices=choice,act='store',def=cvalue,error=error,nargs='+')
                    else
                        call this%cli%add(group=group,switch=switch, help=help, &
                                     required=is_required,act='store',def=cvalue,error=error,nargs='+')
                    endif
                end if
                else
                write(*,*) 'Rank >1 arrays not supported by CLI'
                assert(.false.)
            end if
            assert(error==0)
            call Iterator%Next()
        end do
    end subroutine parameter_handler_add_to_cli_group
  

    subroutine parameter_handler_parse(this)
    !------------------------------------------------------------------
    !< Parse command line arguments and fill parameter lists
    !------------------------------------------------------------------
        class(parameter_handler_t), intent(inout) :: this
        type(parameterlist_t), pointer            :: switches_sublist
        type(parameterlist_t), pointer            :: values_sublist
        integer(ip)                               :: istat, error
        character(len=str_cla_len)                :: switch ! , cvalue
        integer(ip)                               :: ivalue
        character(len=:), allocatable             :: key
        type(ParameterListIterator_t)             :: Iterator
    !------------------------------------------------------------------
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
  

    subroutine parameter_handler_init_cli(this,                     &
                    progname, version, help, description,           &
                    license, authors, examples, epilog, disable_hv, &
                    usage_lun, error_lun, version_lun)
    !------------------------------------------------------------------
    !< Initialize command line interface
    !------------------------------------------------------------------
        implicit none
        class(parameter_handler_t),   intent(inout) :: this
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

        call this%cli%init(progname, version, help, description, &
                         license, authors, examples, epilog, disable_hv, &
                         usage_lun, error_lun, version_lun)
    end subroutine parameter_handler_init_cli  

  
    subroutine parameter_handler_parse_group(this,switches,values,group)
    !------------------------------------------------------------------
    !< Parser group of arguments from command line and update parameter lists
    !------------------------------------------------------------------
        implicit none
        class(parameter_handler_t), intent(inout) :: this
        type(parameterlist_t),      intent(inout) :: switches
        type(parameterlist_t),      intent(inout) :: values
        character(*), optional,     intent(in)    :: group
        integer(ip)                               :: istat, error
        character(len=str_cla_len)                :: switch ! , cvalue
        integer(ip)                               :: ivalue
        character(len=:), allocatable             :: key
        type(ParameterListIterator_t)             :: Iterator
        class(*), pointer                         :: val0
        class(*), pointer                         :: val1(:)
        integer(ip), allocatable                  :: val_ip(:)
        real(rp)   , allocatable                  :: val_rp(:)
        character(str_cla_len)                    :: val_ch
    !------------------------------------------------------------------
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


    subroutine parameter_handler_initialize_lists(this)
    !------------------------------------------------------------------
    !< Initialize parameter lists
    !------------------------------------------------------------------
        implicit none
        class(parameter_handler_t), intent(inout) :: this
        call this%values%init()
        call this%switches%init()
        call this%switches_ab%init()
        call this%helpers%init()
        call this%required%init()
        call this%choices%init()
    end subroutine parameter_handler_initialize_lists

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FEMPAR_PARAMETER_HANDLER_T PROCEDURES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    subroutine fph_process_parameters(this,define_user_parameters_procedure,parse_cla, &
                                      progname, version, help, description,            &
                                      license, authors, examples, epilog, disable_hv,  &
                                      usage_lun, error_lun, version_lun)
    !------------------------------------------------------------------
    !< Initialize lists, publish parameters and fill values from CLI
    !------------------------------------------------------------------
        implicit none
        class(fempar_parameter_handler_t),           intent(inout) :: this
        procedure(define_user_parameters), optional                :: define_user_parameters_procedure
        logical,                           optional, intent(in)    :: parse_cla
        character(*),                      optional, intent(in)    :: progname          !< Program name.
        character(*),                      optional, intent(in)    :: version           !< Program version.
        character(*),                      optional, intent(in)    :: help              !< Help message introducing the CLI usage.
        character(*),                      optional, intent(in)    :: description       !< Detailed description message introducing the program.
        character(*),                      optional, intent(in)    :: license           !< License description.
        character(*),                      optional, intent(in)    :: authors           !< Authors list.
        character(*),                      optional, intent(in)    :: examples(1:)      !< Examples of correct usage.
        character(*),                      optional, intent(in)    :: epilog            !< Epilog message.
        logical,                           optional, intent(in)    :: disable_hv        !< Disable automatic insert of 'help' and 'version' CLAs.
        integer(ip),                       optional, intent(in)    :: usage_lun         !< Unit number to print usage/help
        integer(ip),                       optional, intent(in)    :: version_lun       !< Unit number to print version/license info
        integer(ip),                       optional, intent(in)    :: error_lun         !< Unit number to print error info
        logical :: parse_cla_
    !------------------------------------------------------------------
        call this%free()
        parse_cla_ = .true.; if ( present(parse_cla) ) parse_cla_ = parse_cla
        if ( parse_cla_ ) call this%init_cli(progname, version, help, description,     &
                                      license, authors, examples, epilog, disable_hv,  &
                                      usage_lun, error_lun, version_lun)
        call this%initialize_lists()
        call this%define_fempar_parameters()
        if (present(define_user_parameters_procedure)) then
            this%define_user_parameters => define_user_parameters_procedure
            call this%define_user_parameters()
            else
            nullify(this%define_user_parameters)
        end if

        if ( parse_cla_ ) then 
            call this%add_to_cli()
            call this%parse()
        end if
    end subroutine fph_process_parameters


    subroutine fph_free(this)
    !------------------------------------------------------------------
    !< Free fempar_parameter_handler_t derived type
    !------------------------------------------------------------------
        implicit none
        class(fempar_parameter_handler_t), intent(inout) :: this
    !------------------------------------------------------------------
        call parameter_handler_free(this) 
        nullify(this%define_user_parameters)
    end subroutine fph_free  
  
    subroutine fph_define_fempar_parameters(this)
    !------------------------------------------------------------------
    !< Publish FEMPAR parameters
    !------------------------------------------------------------------
        implicit none
        class(fempar_parameter_handler_t), intent(inout) :: this
       
        ! Environment
        call this%add(environment_type_key, environment_type_cla_name, structured, 'Type of environment')

        call this%fph_static_triang_define_parameters()
        call this%fph_struct_hex_mesh_generator_define_parameters()
        call this%fph_p4est_triang_define_parameters()
        call this%fph_mesh_define_parameters()
        call this%fph_mesh_partitioner_define_parameters()
        call this%fph_fes_define_parameters()
        call this%fph_ils_define_parameters()
        call this%fph_nls_define_parameters()
        call this%fph_dls_define_parameters()
        call this%fph_coarse_fe_handler_define_parameters()
        call this%fph_output_handler_define_parameters()        
    end subroutine fph_define_fempar_parameters
    
    subroutine fph_static_triang_define_parameters(this)
      implicit none
      class(fempar_parameter_handler_t), intent(inout) :: this
      call this%add(static_triang_generate_from_key, &
           static_triang_generate_cla_name, &
           static_triang_generate_from_struct_hex_mesh_generator, &
           'Way to generate static triangulation') !, &
           ! choices = static_triang_generate_cla_choices)
      
      call this%add(static_triang_geometric_interpolation_order_key, &
           static_triang_geometric_interpolation_order_cla_name, &
           1, 'Interpolation order static triangulation geometrical mapping (only partial/unstable support for order>1)' )
      
      call this%add(static_triang_mesh_data_files_dir_path_key, &
           static_triang_mesh_data_files_dir_path_cla_name, &
           '.', &
           'The relative or full file system path to the folder where the mesh data files are located') !, &
           ! choices = static_triang_generate_cla_choices)
      
      call this%add(static_triang_mesh_data_files_prefix_key, &
           static_triang_mesh_data_files_prefix_cla_name, &
           'mesh', &
           'Token string which is used as a prefix to compose the names of the mesh data files as prefix.*') !, &
           ! choices = static_triang_generate_cla_choices)
    end subroutine fph_static_triang_define_parameters
    
    
    subroutine fph_struct_hex_mesh_generator_define_parameters(this)
      implicit none
      class(fempar_parameter_handler_t), intent(inout) :: this
      
      call this%add(struct_hex_mesh_generator_num_dims_key, &
           struct_hex_mesh_generator_num_dims_cla_name, &
           2, &
           'Number of space dimensions') ! , &
           ! choices = struct_hex_mesh_generator_num_dims_cla_choices )

      call this%add(struct_hex_mesh_generator_domain_limits_key, &
           struct_hex_mesh_generator_domain_limits_cla_name, &
           [0.0,1.0,0.0,1.0,0.0,1.0], &
           'Domain interval per dimension')

      call this%add(struct_hex_mesh_generator_is_dir_periodic_key, &
           struct_hex_mesh_generator_is_dir_periodic_cla_name, &
           [0,0,0], &
           'Is the mesh periodic per dimension')

      call this%add(struct_hex_mesh_generator_num_levels_key, &
           struct_hex_mesh_generator_num_levels_cla_name, &
           1, &
           'Number of levels in the triangulation hierarchy required for the MLBDDC preconditioner')

      call this%add(struct_hex_mesh_generator_num_cells_x_dim_key, &
           struct_hex_mesh_generator_num_cells_x_dim_cla_name, &
           [10,10,10], &
           'Number of cells per dimension per level')  

      call this%add(struct_hex_mesh_generator_num_parts_x_dim_key, &
           struct_hex_mesh_generator_num_parts_x_dim_cla_name, &
           [1,1,1], &
           'Number of parts per dimension per level')

    end subroutine fph_struct_hex_mesh_generator_define_parameters

    subroutine fph_p4est_triang_define_parameters(this)
      implicit none
      class(fempar_parameter_handler_t), intent(inout) :: this
      call this%add(p4est_triang_num_dims_key, &
           p4est_triang_num_dims_cla_name, &
           default_p4est_triang_num_dims, &
           'p4est triangulation number of space dimensions')
      
      call this%add(p4est_triang_num_levels_key, &
           p4est_triang_num_levels_cla_name, &
           default_p4est_triang_num_levels, &
           'Number of levels in the triangulation hierarchy required for the MLBDDC preconditioner')

      call this%add(p4est_triang_domain_limits_key, &
           p4est_triang_domain_limits_cla_name, &
           default_p4est_triang_domain_limits, &
           'p4est triangulation domain interval per dimension')

      call this%add(p4est_triang_log_level_key, &
           p4est_triang_log_level_cla_name, &
           default_p4est_triang_log_level, &
           'p4est library level of logging output')

      call this%add(p4est_triang_2_1_k_balance_key, &
           p4est_triang_2_1_k_balance_cla_name, &
           default_p4est_triang_2_1_k_balance,  &
           'value of k for 2:1 k-balanced forest-of-octrees &
           (use with care, at present, only k={0,1} supported/tested)')

      call this%add(p4est_triang_k_ghost_cells_key, &
           p4est_triang_k_ghost_cells_cla_name, &
           default_p4est_triang_k_ghost_cells, &
           'value of k for the k-ghost cells set of each processor &
           (k=0 works for any FE space; k>0 should work depending on the FE space, although NOT tested, use with care)')
    end subroutine fph_p4est_triang_define_parameters

    subroutine fph_mesh_define_parameters(this)
      implicit none
      class(fempar_parameter_handler_t), intent(inout) :: this

      call this%add(mesh_dir_path_key, &
           mesh_dir_path_cla_name, &
           mesh_default_dir_path_output, &
           'The relative or full file system path to the & 
            folder where the mesh data files')

      call this%add(mesh_prefix_key, &
           mesh_prefix_cla_name, &
           mesh_default_prefix_output, &
           'Token string which is used as a prefix to compose & 
            the names of the mesh data files')
    end subroutine fph_mesh_define_parameters
    
    subroutine fph_mesh_partitioner_define_parameters(this)
      implicit none
      class(fempar_parameter_handler_t), intent(inout) :: this
      character(len=:), allocatable                    :: help_string
      help_string = "Directory of the source files" // BRK_LINE // & 
           "To read from GiD mesh" // BRK_LINE // &
           "To read from GiD mesh" 

      call this%add(mesh_partitioner_dir_path_key, &
           mesh_partitioner_dir_path_cla_name, &
           mesh_partitioner_default_dir_path_output, &
           'The relative or full file system path to the & 
            folder where the data files of the mesh partitioner are generated')

      call this%add(mesh_partitioner_prefix_key, &
           mesh_partitioner_prefix_cla_name, &
           mesh_partitioner_default_prefix_output, &
           'Token string which is used as a prefix to compose & 
            the names of the data files generated by the mesh partitioner as prefix.*')
      
      call this%add(dir_path_key, &
           dir_path_cla_name, &
           '.', &
           help_string)
            
      call this%add(prefix_key, &
           prefix_cla_name, &
           'A', &
           'Name of the GiD files')

      call this%add(mesh_partitioner_num_levels_key, &
           mesh_partitioner_num_levels_cla_name, &
           1, &
           'Number of levels in the triangulation hierarchy required for the MLBDDC preconditioner')

      call this%add(mesh_partitioner_num_parts_x_level_key, &
           mesh_partitioner_num_parts_x_level_cla_name, &
           [1], &
           'Number of parts per level') 

      call this%add(mesh_partitioner_strategy_key, &
           mesh_partitioner_strategy_cla_name, &
           metis_part_kway, &
           'Strategy key for partitioner') 

      call this%add(mesh_partitioner_metis_option_debug_key, &
           mesh_partitioner_metis_option_debug_cla_name, &
           2, &
           'METIS debug key') 

      call this%add(mesh_partitioner_metis_option_ufactor_key, &
           mesh_partitioner_metis_option_ufactor_cla_name, &
           30, &
           'METIS option ufactor') 

      call this%add(mesh_partitioner_metis_option_minconn_key, &
           mesh_partitioner_metis_option_minconn_cla_name, &
           0, &
           'METIS option minconn') 

      call this%add(mesh_partitioner_metis_option_contig_key, &
           mesh_partitioner_metis_option_contig_cla_name, &
           1, &
           'METIS option config') 

      call this%add(mesh_partitioner_metis_option_ctype_key, &
           mesh_partitioner_metis_option_ctype_cla_name, &
           METIS_CTYPE_SHEM, &
           'METIS option ctype')

      call this%add(metis_option_iptype_key, &
           mesh_partitioner_metis_option_iptype_cla_name, &
           METIS_IPTYPE_EDGE, &
           'METIS option iptype')
    end subroutine fph_mesh_partitioner_define_parameters
    
    subroutine fph_fes_define_parameters(this)
      implicit none
      class(fempar_parameter_handler_t), intent(inout) :: this

      ! Number of fields in the FE space
      call this%add(fes_num_fields_key, &
           fes_num_fields_cla_name, 1, &
           'Finite element space number of fields')

      ! Number of reference FEs
      call this%add(fes_num_ref_fes_key, &
           fes_num_ref_fes_cla_name, &
           1, &
           'Finite element space number of fields')

      ! set_ids_to_reference_fes: Given a field ID and cell set ID, 
      ! returns the desired reference FE ID
      call this%add(fes_set_ids_ref_fes_key, &
           fes_set_ids_ref_fes_cla_name, &
           [ 1 ], &
           'Set IDs to reference FEs for every field')

      ! Reference FE IDs    
      call this%add(fes_ref_fe_types_key, &
           fes_ref_fe_types_cla_name, &
           fe_type_lagrangian, &
           'Reference finite element types')

      ! Reference FE order (0,1,2,...) fe_space_orders_key
      call this%add(fes_ref_fe_orders_key, &
           fes_ref_fe_orders_cla_name, &
           [ 1 ], &
           'Reference finite element orders')

      ! FE space conformities 
      ! ( true = no face integration needed, false = face integration required )
      call this%add(fes_ref_fe_conformities_key, &
           fes_ref_fe_conformities_cla_name, &
           [ .true. ], &
           'Finite element space conformities')

      ! FE space continuities 
      ! ( true = continuous FE space, false = otherwise )
      call this%add(fes_ref_fe_continuities_key, &
           fes_ref_fe_continuities_cla_name, &
           [ .true. ], &
           'Finite element space continuities')

      ! FE field types (scalar, vector, tensor)
      call this%add(fes_field_types_key, &
           fes_field_types_cla_name, &
           field_type_scalar, &
           'Finite element space field types')

      ! FE field blocks (scalar, vector, tensor)
      call this%add(fes_field_blocks_key, &
           fes_field_blocks_cla_name, &
           [ 1 ], &
           'Finite element space map field_id to block_id')

      ! FE space construction type homogeneous/heterogeneous 
      ! ( .true. = all cells same reference fe, .false. = otherwise ) 
      call this%add(fes_same_ref_fes_all_cells_key, &
           fes_same_ref_fes_all_cells_cla_name, &
           .true., &
           'Finite element space fixed reference fe logical')

    end subroutine fph_fes_define_parameters
     
    subroutine fph_coarse_fe_handler_define_parameters(this)
      implicit none
      class(fempar_parameter_handler_t), intent(inout) :: this

      call this%add(coarse_space_use_vertices_key, &
           coarse_space_use_vertices_cla_name, &
           .true., &
           'Shape functions on vertices')

      call this%add(coarse_space_use_edges_key, &
           coarse_space_use_edges_cla_name, &
           .true., &
           'Shape functions on edges')

      call this%add(coarse_space_use_faces_key, &
           coarse_space_use_faces_cla_name, &
           .true., &
           'Shape functions on faces')

      ! H-curl COARSE FE HANDLER BDDC
      call this%add(bddc_scaling_function_case_key, &
           '--BDDC_SCALING_FUNCTION_CASE', &
           cardinality, &
           'Scaling type for the BDDC weighting operator')
    end subroutine fph_coarse_fe_handler_define_parameters

    subroutine fph_ils_define_parameters(this)
      implicit none
      class(fempar_parameter_handler_t), intent(inout) :: this

      ! Iterative linear solver keys
      call this%add(ils_type_key, & 
           ils_type_cla_name, & 
           rgmres_name,  & 
           'Iterative linear solver type')

      call this%add(ils_rtol_key, & 
           ils_rtol_cla_name, & 
           default_rtol, & 
           'Relative tolerance for the stopping criteria')

      call this%add(ils_atol_key, & 
           ils_atol_cla_name, & 
           default_atol, & 
           'Absolute tolerance for the stopping criteria')

      call this%add(ils_stopping_criterium_key, & 
           ils_stopping_criterium_cla_name, & 
           default_rgmres_stopping_criteria, & 
           help='Stopping criterium type')

      call this%add(ils_output_frequency_key, & 
           ils_output_frequency_cla_name, & 
           default_output_frequency, & 
           'Frequency for output printing') 

      call this%add(ils_max_num_iterations_key, & 
           ils_max_num_iterations_cla_name, & 
           default_max_num_iterations, & 
           'Maximum number of iterations')

      call this%add(ils_track_convergence_history_key, & 
           ils_track_convergence_history_cla_name, & 
           default_track_convergence_history,'Track convergence history')

      call this%add(ils_max_dim_krylov_basis_key, & 
           ils_max_dim_krylov_basis_cla_name, & 
           default_dkrymax, & 
           'Maximum dimension of Krylov basis') 

      call this%add(ils_orthonorm_strategy_key, & 
           ils_orthonorm_strategy_cla_name, & 
           default_orthonorm_strat, & 
           'Orthonormalization strategy (for GMRES)')

      call this%add(ils_relaxation_key, & 
           ils_relaxation_cla_name, & 
           default_richardson_relaxation, & 
           'Relaxation value (for Richardson)')
    end subroutine fph_ils_define_parameters
     

    subroutine fph_dls_define_parameters(this)
      implicit none
      class(fempar_parameter_handler_t), intent(inout) :: this
      integer(ip) :: i
#ifdef UMFPACK
      integer(ip) :: umfpack_control(UMFPACK_INFO)
#endif 

      call this%add(dls_type_key, & 
           dls_type_cla_name, & 
           pardiso_mkl, & 
           'Direct solver type')

      call this%add(pardiso_mkl_message_level, & 
           pardiso_mkl_message_level_cla_name, & 
           pardiso_mkl_default_message_level, & 
           'PARDISO message level')

      call this%add(pardiso_mkl_iparm, & 
           pardiso_mkl_iparm_cla_name, & 
           [ (pardiso_mkl_default_iparm, i=1,64) ], & 
           'PARDISO parameters')     

      call this%add(pardiso_mkl_matrix_type, &
           pardiso_mkl_matrix_type_cla_name, & 
           pardiso_mkl_default_matrix_type, & 
           'PARDISO matrix type')

#ifdef UMFPACK
      call umfpack_di_defaults(umfpack_control)
      call this%add(umfpack_control_params, &
           umfpack_control_params_cla_name, &
           umfpack_control, &
           'UMFPACK control parameters')
#endif  
    end subroutine fph_dls_define_parameters

    subroutine fph_nls_define_parameters(this)
      implicit none
      class(fempar_parameter_handler_t), intent(inout) :: this
      call this%add(nls_rtol_key, &
           nls_rtol_cla_name, &
           default_nls_rtol, &
           'Relative tolerance for the nonlinear solvers stopping criteria')

      call this%add(nls_atol_key, &
           nls_atol_cla_name, &
           default_nls_atol, &
           'Absolute tolerance for the nonlinear solvers stopping criteria')

      call this%add(nls_stopping_criterium_key, &
           nls_stopping_criterium_cla_name, &
           default_nls_stopping_criterium, &
           'Nonlinear solvers stopping criterium type')

      call this%add(nls_max_num_iterations_key, &
           nls_max_num_iterations_cla_name, &
           default_nls_max_iter, &
           'Nonlinear solvers maximum number of iterations')

      call this%add(nls_print_iteration_output_key, &
           nls_print_iteration_output_cla_name, &
           default_nls_print_iteration_output, &
           'Print output per nonlinear solver iteration')
    end subroutine fph_nls_define_parameters
     
    subroutine fph_output_handler_define_parameters(this)
      implicit none
      class(fempar_parameter_handler_t), intent(inout) :: this
      ! General
      call this%add(output_handler_dir_path_key, &
           output_handler_dir_path_cla_name, &
           '.', &
           'Output Directory')
      
      call this%add(output_handler_prefix_key, &
           output_handler_prefix_cla_name, &
           '.', &
           'Prefix to generate the names of the generated files')
      
      call this%add(output_handler_static_grid_key, &
           output_handler_static_grid_cla_name, &
           output_handler_static_grid_default, &
           'Grid does not change along all time steps')
      
      call this%add(output_handler_vtk_format_key, &
           output_handler_vtk_format_cla_name, &
           output_handler_vtk_format_default, &
           'Prefix to generate the names of the generated files')
      
      call this%add(output_handler_xh5_strategy_key, &
           output_handler_xh5_strategy_cla_name, &
           output_handler_xh5_strategy_default, &
           'Parallel IO strategy (800: Contiguous hyperslabs, 801: Dataset per process). See XH5For docs.', &
           choices = output_handler_xh5_strategy_choices)
      
      call this%add(output_handler_xh5_info_key, &
           output_handler_xh5_info_cla_name, &
           output_handler_xh5_info_default, &
           'HDF5 I/O low level configuration based on MPI info')
      
    end subroutine fph_output_handler_define_parameters
     
    subroutine fph_Get_0D(this, key, value)
    !------------------------------------------------------------------
    !< Get 0D value given its key
    !------------------------------------------------------------------
        implicit none
        class(fempar_parameter_handler_t) , intent(in)    :: this
        character(len=*),                   intent(in)    :: key
        class(*),                           intent(inout) :: value
        integer(ip)                                       :: error
    !------------------------------------------------------------------
        assert(this%values%isAssignable(key, value))
        error = this%values%Get(key=key, value=value)
        assert(error==0)
    end subroutine fph_Get_0D


    subroutine fph_Get_1D(this, key, value)
    !------------------------------------------------------------------
    !< Get 1D value given its key
    !------------------------------------------------------------------
        implicit none
        class(fempar_parameter_handler_t) , intent(in)    :: this
        character(len=*),                   intent(in)    :: key
        class(*),                           intent(inout) :: value(:)
        integer(ip)                                       :: error
    !------------------------------------------------------------------
        assert(this%values%isAssignable(key, value))
        error = this%values%Get(key=key, value=value)
        assert(error==0)
    end subroutine fph_Get_1D


    subroutine fph_GetAsString(this, key, string)
    !------------------------------------------------------------------
    !< Get 0D value given its key
    !------------------------------------------------------------------
        implicit none
        class(fempar_parameter_handler_t) , intent(in)    :: this
        character(len=*),                   intent(in)    :: key
        character(len=:),      allocatable, intent(inout) :: string
        integer(ip)                                       :: error
    !------------------------------------------------------------------
        assert(this%values%isAssignable(key, 'string'))
        error = this%values%GetAsString(key=key, string=string)
        assert(error==0)
    end subroutine fph_GetAsString


    function fph_get_dir_path(this)
    !------------------------------------------------------------------
    !< Get dir_path
    !------------------------------------------------------------------
        implicit none
        class(fempar_parameter_handler_t) , intent(in) :: this
        character(len=:),      allocatable             :: fph_get_dir_path
    !------------------------------------------------------------------
        call this%GetAsString(key=dir_path_key, String = fph_get_dir_path)
    end function fph_get_dir_path 


    function fph_get_prefix(this)
    !------------------------------------------------------------------
    !< Get prefix
    !------------------------------------------------------------------
        implicit none
        class(fempar_parameter_handler_t) , intent(in) :: this
        character(len=:),      allocatable             :: fph_get_prefix
    !------------------------------------------------------------------
        call this%GetAsString(key=prefix_key, String = fph_get_prefix)
    end function fph_get_prefix

end module fempar_parameter_handler_names 
