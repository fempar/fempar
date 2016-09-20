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
module command_line_parameters_names
  use types_names
  !use Data_Type_Command_Line_Interface
  use flap, only : command_line_interface
# include "debug.i90"

  implicit none
  private

  type test_geometry_params_t
     character(len=:), allocatable :: default_dir_path
     character(len=:), allocatable :: default_prefix
     character(len=:), allocatable :: default_dir_path_out
   contains
     procedure :: set_default_params => test_geometry_set_default_params
  end type test_geometry_params_t

  ! Types
  public :: test_geometry_params_t

  ! Functions
  public :: cli_add_params

contains

  subroutine test_geometry_set_default_params(params)
    use fempar_names
    implicit none
    class(test_geometry_params_t), intent(inout) :: params
    ! IO parameters
    params%default_dir_path     = '.'
    params%default_prefix       = 'cyl'
    params%default_dir_path_out = '.'
  end subroutine test_geometry_set_default_params

  !==================================================================================================
  subroutine cli_add_params(cli,params)
    implicit none
    type(Command_Line_Interface)    , intent(inout) :: cli
    type(test_geometry_params_t)          , intent(in)    :: params
    !class(test_geometry_parallel_params_t), intent(inout) :: params
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
    
  end subroutine cli_add_params


end module command_line_parameters_names

!****************************************************************************************************
!****************************************************************************************************
program test_geometry
  use types_names
  use iso_c_binding
  use sisl_names
  use geometry_names
  use fempar_names
  !use Data_Type_Command_Line_Interface
  use flap, only : command_line_interface
  use command_line_parameters_names
  implicit none
  type(geometry_t)      :: geometry
  type(line_t), pointer :: line
  type(point_t)         :: point
  real(rp)              :: param
  !real(rp)              :: point(3),param

  character(len=256)       :: dir_path, dir_path_out
  character(len=256)       :: prefix, filename
  integer(ip) :: lunio, istat
  type(Command_Line_Interface):: cli 
  type(test_geometry_params_t)     :: test_params

  call meminit

  ! Read IO parameters
  call read_flap_cli_test_geometry(cli,test_params)
  call cli%get(switch='-d'  ,val=dir_path    ,error=istat); check(istat==0)
  call cli%get(switch='-pr' ,val=prefix      ,error=istat); check(istat==0)

  ! Get line 3 of geometry and test its TBPs
  !call geometry%read(dir_path,prefix)
  line => geometry%get_line(3)

  !point(1) = -10.0_rp
  !point(2) = 0.0_rp
  !point(3) = 0.0_rp
  call point%init(0.0_rp)
  call point%set(1,-10.0_rp)
  param = line%get_parameter(point,1.0e-8_rp)
  write(*,*) 'Line parameter', param

  param = 0.75_rp
  call line%evaluate(param,point)
  write(*,*) 'Line evaluation', point%get_value()

  ! Get line 1 of geometry and test its TBPs
  !call geometry%read(dir_path,prefix)
  line => geometry%get_line(1)

  !point(1) = 0.0_rp
  !point(2) = 10.0_rp
  !point(3) = 20.0_rp
  call point%init( (/0.0_rp,10.0_rp,20.0_rp/) )

  param = line%get_parameter(point,1.0e-8_rp)
  write(*,*) 'Line parameter', param

  param = 0.6_rp
  call line%evaluate(param,point)
  write(*,*) 'Line evaluation', point%get_value()

contains

  !==================================================================================================
  subroutine read_flap_cli_test_geometry(cli,test_params)
    !use Data_Type_Command_Line_Interface
    use command_line_parameters_names
    use fempar_names
    implicit none
    type(Command_Line_Interface) , intent(out)   :: cli
    type(test_geometry_params_t), intent(inout) :: test_params
    
    ! Locals
    integer(ip)                 :: istat

    ! Initialize Command Line Interface
    call cli%init(progname    = 'test_geometry',                                                    &
         &        version     = '',                                                                 &
         &        authors     = '',                                                                 &
         &        license     = '',                                                                 &
         &        description =  'Parallel FEMPAR driver to test geometry.', &
         &        examples    = ['test_geometry -h  ', 'test_geometry -h  ' ])
    
    ! Set Parallel parameters
    call test_params%set_default_params()
    call cli_add_params(cli,test_params) 
    
    call cli%parse(error=istat)
    check(istat == 0)
    
  end subroutine read_flap_cli_test_geometry

end program test_geometry
