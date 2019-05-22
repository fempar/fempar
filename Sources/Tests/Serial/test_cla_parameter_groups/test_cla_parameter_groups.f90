module custom_parameter_generator_names

    use fempar_names

implicit none
#include "debug.i90"

    type :: cla_parameter_groups_t
    contains
        procedure, public :: process_parameters => cla_parameter_groups_process_parameters
        procedure, public :: print              => cla_parameter_groups_print
    end type

contains

    subroutine add_to_group( group)
        character(len=*),              intent(in)    :: group
        call parameter_handler%add('switch1', '--switch1', 1, 'switch1 help', switch_ab='-s1', required=.false., group=group)
        call parameter_handler%add('switch2', '--switch2', 2, 'switch2 help', switch_ab='-s2', required=.false., group=group)
        call parameter_handler%add('switch3', '--switch3', 3, 'switch3 help', switch_ab='-s3', required=.false., group=group)
    end subroutine


    subroutine define_parameters()
        call parameter_handler%add('commonswitch', '--commonswitch', 1, 'common switch help', switch_ab='-cs', required=.true., choices=[1,2,3])
        call add_to_group('group1')
        call add_to_group('group2')
        call add_to_group('group3')
    end subroutine


    subroutine cla_parameter_groups_process_parameters(this)
        class(cla_parameter_groups_t), intent(inout) :: this
        call parameter_handler%create(examples=['test_cla_parameter_groups -cs 7 group1 -s1 8 group2 -s2 9'])
        call parameter_handler%process_parameters(define_parameters)
    end subroutine


    subroutine cla_parameter_groups_print(this)
        class(cla_parameter_groups_t), intent(inout) :: this
        type(ParameterList_t), pointer               :: values
        type(ParameterList_t), pointer               :: switches
        type(ParameterList_t), pointer               :: switches_ab
        type(ParameterList_t), pointer               :: helpers
        type(ParameterList_t), pointer               :: required
        type(ParameterList_t), pointer               :: choices

        switches    => parameter_handler%get_switches()
        switches_ab => parameter_handler%get_switches_ab()
        helpers     => parameter_handler%get_helpers()
        required    => parameter_handler%get_required()
        values      => parameter_handler%get_values()
        choices     => parameter_handler%get_choices()

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

    call my_parameters%process_parameters()
    call my_parameters%print()

    call FEMPAR_FINALIZE()

end program
