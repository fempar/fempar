module command_line_parameters_names
  use serial_names
  use par_names
  implicit none
#include "debug.i90" 
  private
  
  ! Types
  type par_test_reference_fe_parameters_t
     ! Input problem location
     character(len=256) :: dir_path
     character(len=256) :: prefix
     
     ! Number of parts in which the problem was split
     integer(ip)                       :: nparts
     type(Type_Command_Line_Interface) :: cli
  contains
     procedure :: set_default_params => par_test_reference_fe_parameters_set_default_params
     procedure :: set_cli            => par_test_reference_fe_parameters_set_cli
     procedure :: parse              => par_test_reference_fe_parameters_parse
  end type par_test_reference_fe_parameters_t

  ! Types
  public :: par_test_reference_fe_parameters_t

contains
  !==================================================================================================
  subroutine par_test_reference_fe_parameters_set_default_params(params)
    implicit none
    class(par_test_reference_fe_parameters_t), intent(inout) :: params
    params%dir_path = ''
    params%prefix = ''
  end subroutine par_test_reference_fe_parameters_set_default_params

  !==================================================================================================
  subroutine par_test_reference_fe_parameters_set_cli(params)
    implicit none
    class(par_test_reference_fe_parameters_t), intent(inout) :: params
    integer(ip) :: error
    
    ! Initialize Command Line Interface
    call params%cli%init(progname    = 'par_test_reference_fe', &
         &               version     = '',                               &
         &               authors     = '',                               &
         &               license     = '',                               &
         &               description = "FEMPAR parallel test driver", &
         &               examples    = ['par_test_reference_fe -h'] )

    ! Set Command Line Arguments
    call params%cli%add(switch='--dir-path',switch_ab='-dir-path',help='Absolute or relative PATH to the partitioned&
                       & problem. Must end with /',required=.true., act='store', error=error)
    check(error==0)

    call params%cli%add(switch='--prefix',switch_ab='-prefix',help='Prefix for all input files (mesh, conditions, etc.).& 
                       & E.g., if these files were generated from square.gid GiD project, then --prefix square.',& 
                       & required=.true., act='store', error=error)
    check(error==0)

    call params%cli%add(switch='--nparts',switch_ab='-nparts',help='Number of parts in which the problem was split.',& 
                       & required=.true., act='store', error=error)
    check(error==0)
  end subroutine par_test_reference_fe_parameters_set_cli

  !==================================================================================================
  subroutine par_test_reference_fe_parameters_parse(params)
    implicit none
    class(par_test_reference_fe_parameters_t), intent(inout) :: params
    integer(ip)                                                        :: error
    character(len=256)                                                 :: aux_string
    logical                                                            :: aux_logical

    call params%cli%parse(error=error)
    check(error==0)

    call params%cli%get(switch='-dir-path',val=params%dir_path,error=error); check(error==0)
    call params%cli%get(switch='-prefix',val=params%prefix,error=error); check(error==0)
    call params%cli%get(switch='-nparts',val=params%nparts,error=error); check(error==0)
  end subroutine par_test_reference_fe_parameters_parse
end module command_line_parameters_names



program par_test_reference_fe
  !----------------------------------------------------------
  ! Parallel partitioner test
  !----------------------------------------------------------
  use serial_names
  use par_names
  use command_line_parameters_names

  implicit none
#include "debug.i90" 

  ! Our data
  type(par_context_t)                             :: w_context, p_context, q_context, b_context
  type(par_environment_t)                         :: p_env
  type(par_mesh_t)                                :: p_mesh
  type(par_conditions_t)                          :: p_cond
  type(par_triangulation_t)                       :: p_trian
  type(par_fe_space_t)                            :: p_fe_space
  type(p_reference_fe_t)                          :: reference_fe_array_one(1)
  

  integer(ip)              :: num_levels
  integer(ip), allocatable :: id_parts(:), num_parts(:)
  
  type(par_test_reference_fe_parameters_t) :: test_params


  call meminit

  ! Start parallel execution
  call par_context_create (w_context)

  call test_params%set_default_params()
  call test_params%set_cli()
  call test_params%parse()

  ! This test only works with two levels.
  num_levels = 2
  call memalloc(num_levels, id_parts , __FILE__, __LINE__)
  call memalloc(num_levels, num_parts, __FILE__, __LINE__)

  num_parts = (/test_params%nparts, 1/)
  id_parts = (/w_context%iam+1, 1/)


  ! Create p_context and q_context splitting w_context
  if(w_context%iam < num_parts(1)) then
     call par_context_create ( 1, p_context, q_context, w_context )
  else
     call par_context_create ( 2, q_context, p_context, w_context )
  end if
  assert ( (p_context%iam >=0 .and. q_context%iam < 0) .or. (p_context%iam < 0 .and. q_context%iam >= 0))
  
  ! Create b_context as an intercommunicator among p_context <=> q_context 
  call par_context_create ( w_context, p_context, q_context, b_context )

  ! Create parallel environment
  call par_environment_create( p_env,& 
                               w_context,& 
                               p_context,& 
                               q_context,&
                               b_context,&
                               num_levels,&
                               id_parts, & 
                               num_parts )

  ! Read mesh
  call par_mesh_read ( test_params%dir_path, test_params%prefix, p_env, p_mesh )

  ! Read boundary conditions
  call par_conditions_read(test_params%dir_path, test_params%prefix, p_mesh%f_mesh%npoin, p_env, p_cond)

  ! Generate triangulation
  call par_mesh_to_triangulation (p_mesh, p_trian, p_cond)
  
  ! Simple case
  reference_fe_array_one(1) =  make_reference_fe ( topology = topology_quad, &
                                                   fe_type = fe_type_lagrangian, &
                                                   number_dimensions = 2, &
                                                   order = 1, &
                                                   field_type = field_type_scalar, &
                                                   continuity = .true. )
  
  call p_fe_space%par_fe_space_create( par_triangulation = p_trian, &
                                       par_boundary_conditions = p_cond, &
                                       reference_fe_phy = reference_fe_array_one )

  call p_fe_space%par_fe_space_fill_dof_info()
  
  call p_fe_space%par_fe_space_compute_dof_import()
  !call p_fe_space%par_fe_space_print()
  
  call p_fe_space%par_fe_space_free()
  call reference_fe_array_one(1)%free()
  
  call par_triangulation_free(p_trian)
  call par_conditions_free (p_cond)
  call par_mesh_free (p_mesh)

  call memfree(id_parts , __FILE__, __LINE__)
  call memfree(num_parts, __FILE__, __LINE__)

  call par_environment_free (p_env)
  call par_context_free ( b_context, .false. )
  call par_context_free ( p_context, .false. )
  call par_context_free ( q_context, .false. )
  call par_context_free ( w_context )

  call memstatus

end program par_test_reference_fe
