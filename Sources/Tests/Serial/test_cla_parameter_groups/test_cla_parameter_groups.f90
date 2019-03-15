module custom_parameter_generator_names

    use fempar_names

implicit none
#include "debug.i90"

    type, extends(parameter_handler_t) :: cla_parameter_groups_t
    contains
        procedure, public :: define_parameters  => cla_parameter_groups_define_parameters
        procedure, public :: add_to_group       => cla_parameter_groups_add_to_group
        procedure, public :: print              => cla_parameter_groups_print
    end type

contains

    subroutine cla_parameter_groups_define_parameters(this)
        class(cla_parameter_groups_t), intent(inout) :: this
        type(ParameterList_t), pointer               :: switches
        type(ParameterList_t), pointer               :: switches_ab
        type(ParameterList_t), pointer               :: helpers
        type(ParameterList_t), pointer               :: required
        type(ParameterList_t), pointer               :: values
        type(ParameterList_t), pointer               :: choices
        integer(ip)                                  :: error

        switches    => this%get_switches()
        switches_ab => this%get_switches_ab()
        helpers     => this%get_helpers()
        required    => this%get_required()
        values      => this%get_values()
        choices     => this%get_choices()

        call this%add('commonswitch', '--commonswitch', 1, 'common switch help', switch_ab='-cs', required=.true., choices=[1,2,3])

        call this%add_to_group('group1')
        call this%add_to_group('group2')
        call this%add_to_group('group3')

    end subroutine


    subroutine cla_parameter_groups_add_to_group(this, group)
        class(cla_parameter_groups_t), intent(inout) :: this
        character(len=*),              intent(in)    :: group
        integer(ip)                                  :: error

        call this%add('switch1', '--switch1', 1, 'switch1 help', switch_ab='-s1', required=.false., group=group)
        call this%add('switch2', '--switch2', 2, 'switch2 help', switch_ab='-s2', required=.false., group=group)
        call this%add('switch3', '--switch3', 3, 'switch3 help', switch_ab='-s3', required=.false., group=group)
    end subroutine



    subroutine cla_parameter_groups_print(this)
        class(cla_parameter_groups_t), intent(inout) :: this
        type(ParameterList_t), pointer               :: values
        type(ParameterList_t), pointer               :: switches
        type(ParameterList_t), pointer               :: switches_ab
        type(ParameterList_t), pointer               :: helpers
        type(ParameterList_t), pointer               :: required
        type(ParameterList_t), pointer               :: choices
        integer(ip)                                  :: error

        switches    => this%get_switches()
        switches_ab => this%get_switches_ab()
        helpers     => this%get_helpers()
        required    => this%get_required()
        values      => this%get_values()
        choices     => this%get_choices()


        call switches%print(prefix='| Switches | ')
        call switches_ab%print(prefix='| Abreviated Switches | ')
        call helpers%print(prefix='| Helpers | ')
        call required%print(prefix='| Required | ')
        call values%print(prefix='| Values | ')
        call choices%print(prefix='| Choices | ')
    end subroutine

end module custom_parameter_generator_names


program test_cla_parameter_groups

    use custom_parameter_generator_names
    use fempar_names

implicit none

    type(cla_parameter_groups_t) :: my_parameters

    call FEMPAR_INIT()

    call my_parameters%create(examples=['test_cla_parameter_groups -cs 7 group1 -s1 8 group2 -s2 9'])
    call my_parameters%print()
    call my_parameters%free()


    call FEMPAR_FINALIZE()


end program
