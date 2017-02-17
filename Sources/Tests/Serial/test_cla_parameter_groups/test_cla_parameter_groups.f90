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
        integer(ip)                                  :: error

        switches    => this%get_switches()
        switches_ab => this%get_switches_ab()
        helpers     => this%get_helpers()
        required    => this%get_required()
        values      => this%get_values()

        error = switches%set(   key='commonswitch', value='--commonswitch'); assert(error == 0)
        error = switches_ab%set(key='commonswitch', value='-cs');            assert(error == 0)
        error = helpers%set(    key='commonswitch', value='common switch');  assert(error == 0)
        error = required%set(   key='commonswitch', value=.false.);          assert(error == 0)
        error = values%set(     key='commonswitch', value=1);                assert(error == 0)

        call this%add_to_group('group1')
        call this%add_to_group('group2')
        call this%add_to_group('group3')

    end subroutine


    subroutine cla_parameter_groups_add_to_group(this, group)
        class(cla_parameter_groups_t), intent(inout) :: this
        character(len=*),              intent(in)    :: group

        type(ParameterList_t), pointer               :: switches
        type(ParameterList_t), pointer               :: switches_ab
        type(ParameterList_t), pointer               :: helpers
        type(ParameterList_t), pointer               :: required
        type(ParameterList_t), pointer               :: values

        type(ParameterList_t), pointer               :: switches_group
        type(ParameterList_t), pointer               :: switches_ab_group
        type(ParameterList_t), pointer               :: helpers_group
        type(ParameterList_t), pointer               :: required_group
        type(ParameterList_t), pointer               :: values_group

        integer(ip)                                  :: error

        switches    => this%get_switches()
        switches_ab => this%get_switches_ab()
        helpers     => this%get_helpers()
        required    => this%get_required()
        values      => this%get_values()


        switches_group => switches%NewSubList(group)
        switches_ab_group => switches_ab%NewSubList(group)
        helpers_group => helpers%NewSubList(group)
        required_group => required%NewSubList(group)
        values_group => values%NewSubList(group)


        error = switches_group%set(key='switch1', value='--switch1'); assert(error == 0)
        error = switches_group%set(key='switch2', value='--switch2'); assert(error == 0)
        error = switches_group%set(key='switch3', value='--switch3'); assert(error == 0)

        error = switches_ab_group%set(key='switch1', value='-s1'); assert(error == 0)
        error = switches_ab_group%set(key='switch2', value='-s2'); assert(error == 0)
        error = switches_ab_group%set(key='switch3', value='-s3'); assert(error == 0)

        error = helpers_group%set(key='switch1', value='switch1'); assert(error == 0)
        error = helpers_group%set(key='switch2', value='switch2'); assert(error == 0)
        error = helpers_group%set(key='switch3', value='switch3'); assert(error == 0)

        error = required_group%set(key='switch1', value=.false.); assert(error == 0)
        error = required_group%set(key='switch2', value=.false.); assert(error == 0)
        error = required_group%set(key='switch3', value=.false.); assert(error == 0)

        error = values_group%set(key='switch1', value=1); assert(error == 0)
        error = values_group%set(key='switch2', value=2); assert(error == 0)
        error = values_group%set(key='switch3', value=3); assert(error == 0)

    end subroutine



    subroutine cla_parameter_groups_print(this)
        class(cla_parameter_groups_t), intent(inout) :: this
        type(ParameterList_t), pointer               :: values
        type(ParameterList_t), pointer               :: switches
        type(ParameterList_t), pointer               :: switches_ab
        type(ParameterList_t), pointer               :: helpers
        type(ParameterList_t), pointer               :: required
        integer(ip)                                  :: error

        switches    => this%get_switches()
        switches_ab => this%get_switches_ab()
        helpers     => this%get_helpers()
        required    => this%get_required()
        values      => this%get_values()


        call switches%print(prefix='| Switches | ')
        call switches_ab%print(prefix='| Abreviated Switches | ')
        call helpers%print(prefix='| Helpers | ')
        call required%print(prefix='| Required | ')
        call values%print(prefix='| Values | ')
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


    call FEMPAR_FINALIZE()


end program
