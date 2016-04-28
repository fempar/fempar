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
module par_command_line_parameters_names
  use types_names
  use Data_Type_Command_Line_Interface
# include "debug.i90"

  implicit none
  private

  type par_test_triangulations_params_t
     character(len=:), allocatable :: default_dir_path
     character(len=:), allocatable :: default_prefix
     character(len=:), allocatable :: default_dir_path_out
   contains
     procedure :: set_default_params => par_test_triangulations_set_par_default_params
  end type par_test_triangulations_params_t

  ! Types
  public :: par_test_triangulations_params_t

  ! Functions
  public :: cli_add_par_params

contains

  subroutine par_test_triangulations_set_par_default_params(params)
    use serial_names
    implicit none
    class(par_test_triangulations_params_t), intent(inout) :: params
    ! IO parameters
    params%default_dir_path     = 'data/'
    params%default_prefix       = 'square'
    params%default_dir_path_out = 'output/'
  end subroutine par_test_triangulations_set_par_default_params

  !==================================================================================================
  subroutine cli_add_par_params(cli,params)
    implicit none
    type(Type_Command_Line_Interface)    , intent(inout) :: cli
    type(par_test_triangulations_params_t)          , intent(in)    :: params
    !class(par_test_triangulations_parallel_params_t), intent(inout) :: par_params
    ! Locals
    integer(ip) :: error
    character   :: aux_string

    ! IO parameters
    call cli%add(switch='--dir_path',switch_ab='-d',                              &
         &       help='Directory of the source files',required=.false., act='store',                &
         &       def=trim(params%default_dir_path),error=error)
    check(error==0)
    call cli%add(switch='--prefix',switch_ab='-pr',help='Name of the GiD files',  &
         &       required=.false.,act='store',def=trim(params%default_prefix),error=error) 
    check(error==0)
    call cli%add(switch='--dir_path_out',switch_ab='-out',help='Output Directory',&
         &       required=.false.,act='store',def=trim(params%default_dir_path_out),error=error)
    check(error==0)
    
  end subroutine cli_add_par_params


end module par_command_line_parameters_names

!****************************************************************************************************
!****************************************************************************************************

program par_test_triangulations
  use serial_names
  use par_names
  !use JP_par_triangulation_names
  !use JP_par_mesh_to_triangulation_names
  use Data_Type_Command_Line_Interface
  use par_command_line_parameters_names
  use mpi
  
  implicit none
#include "debug.i90"
  ! Our data
  type(par_context_t)                    :: p_context
  type(par_environment_t)                :: p_env
  type(par_mesh_t)                       :: p_mesh
  !type(JP_par_triangulation_t)           :: p_trian
  type(par_conditions_t)                 :: p_cond
  type(par_test_triangulations_params_t) :: test_params

  ! Arguments
  character(len=256)       :: dir_path, dir_path_out
  character(len=256)       :: prefix, filename

  integer(ip) :: lunio, istat

  type(Type_Command_Line_Interface):: cli 
 
  !call meminit

  ! Start parallel execution
  !call par_context_create (p_context)
  !call par_environment_create(p_env,p_context)

  ! Read IO parameters
  !call read_flap_cli_par_test_triangulations(cli,test_params)
 
  ! Read mesh
  !call cli%get(switch='-d'  ,val=dir_path    ,error=istat); check(istat==0)
  !call cli%get(switch='-pr' ,val=prefix      ,error=istat); check(istat==0)
  !call cli%get(switch='-out',val=dir_path_out,error=istat); check(istat==0)
  !call par_mesh_read (dir_path, prefix, p_env, p_mesh)

  !call mpi_barrier(p_context%icontxt,istat)
  
  ! Read conditions 
  !call par_conditions_read (dir_path, prefix, p_mesh%f_mesh%npoin, p_env, p_cond)

  ! Generate efs and its boundary conditions
  !call generate_efs(p_mesh%f_mesh, p_cond%f_conditions)

  ! Construct triangulation
  !call JP_par_mesh_to_triangulation ( p_mesh, p_trian)

  ! Print triangulation
  !call p_trian%print(stderr)

  !call p_trian%free()
  !call par_conditions_free ( p_cond )
  !call par_mesh_free (p_mesh)

  !call par_environment_free (p_env)
  !call par_context_free ( p_context, .false. )
  call memstatus

contains

  !==================================================================================================
  subroutine read_flap_cli_par_test_triangulations(cli,test_params)
    use Data_Type_Command_Line_Interface
    use par_command_line_parameters_names
    use serial_names
    implicit none
    type(Type_Command_Line_Interface), intent(out)   :: cli
    type(par_test_triangulations_params_t)      , intent(inout) :: test_params
    
    ! Locals
    integer(ip)                 :: istat

    ! Initialize Command Line Interface
    call cli%init(progname    = 'par_test_triangulations',                                                     &
         &        version     = '',                                                                 &
         &        authors     = '',                                                                 &
         &        license     = '',                                                                 &
         &        description =  'Parallel FEMPAR driver to trinagulations.', &
         &        examples    = ['par_test_triangulations -h  ', 'par_test_triangulations -h  ' ])
    
    ! Set Parallel parameters
    call test_params%set_default_params()
    call cli_add_par_params(cli,test_params) 
    
    call cli%parse(error=istat)
    check(istat == 0)
    
  end subroutine read_flap_cli_par_test_triangulations

end program
