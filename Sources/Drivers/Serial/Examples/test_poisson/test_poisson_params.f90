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
module test_poisson_params_names
  use serial_names
# include "debug.i90"

  implicit none
  private

  type test_poisson_params_t 
     private 
     ! IO parameters
     character(len=:), allocatable :: default_dir_path
     character(len=:), allocatable :: default_prefix
     character(len=:), allocatable :: default_dir_path_out
     character(len=:), allocatable :: default_fe_formulation
     
     type(Command_Line_Interface):: cli 

     ! IO parameters
     character(len=256)            :: dir_path
     character(len=256)            :: prefix
     character(len=256)            :: dir_path_out
     character(len=256)            :: fe_formulation

     
   contains
     procedure, non_overridable             :: create       => test_poisson_create
     procedure, non_overridable, private    :: set_default  => test_poisson_set_default
     procedure, non_overridable, private    :: add_to_cli   => test_poisson_add_to_cli
     procedure, non_overridable             :: parse        => test_poisson_parse 
     procedure, non_overridable             :: free         => test_poisson_free
     procedure, non_overridable             :: get_dir_path
     procedure, non_overridable             :: get_prefix
     procedure, non_overridable             :: get_dir_path_out
     procedure, non_overridable             :: get_fe_formulation
  end type test_poisson_params_t

  ! Types
  public :: test_poisson_params_t

contains

  subroutine test_poisson_create(this)
    implicit none
    class(test_poisson_params_t), intent(inout) :: this
    
    call this%free()
    
     ! Initialize Command Line Interface
    call this%cli%init(progname    = 'test_poisson',                                                     &
         &        version     = '',                                                                 &
         &        authors     = '',                                                                 &
         &        license     = '',                                                                 &
         &        description =  'FEMPAR driver to test the new serial fe space.', &
         &        examples    = ['test_poisson -h  ', 'test_poisson -h  ' ])
    
    call this%set_default()
    call this%add_to_cli()
  end subroutine test_poisson_create
  
  subroutine test_poisson_set_default(this)
    implicit none
    class(test_poisson_params_t), intent(inout) :: this
    ! IO parameters
    this%default_dir_path       = 'data/'
    this%default_prefix         = 'square'
    this%default_dir_path_out   = 'output/'
    this%default_fe_formulation = 'cG'
  end subroutine test_poisson_set_default
  
  !==================================================================================================
  subroutine test_poisson_add_to_cli(this)
    implicit none
    class(test_poisson_params_t) , intent(inout) :: this

    ! Locals
    integer(ip) :: error

    ! IO parameters
    call this%cli%add(switch='--dir-path',switch_ab='-d',                              &
         &            help='Directory of the source files',required=.false., act='store',                &
         &            def=trim(this%default_dir_path),error=error)
    check(error==0)
    call this%cli%add(switch='--prefix',switch_ab='-p',help='Name of the GiD files',  &
         &            required=.false.,act='store',def=trim(this%default_prefix),error=error) 
    check(error==0)
    call this%cli%add(switch='--dir-path-out',switch_ab='-o',help='Output Directory',&
         &            required=.false.,act='store',def=trim(this%default_dir_path_out),error=error)
    check(error==0)  
    call this%cli%add(switch='--fe-formulation',switch_ab='-f',help='cG or dG FE formulation for Poisson problem',&
         &            required=.false.,act='store',def=trim(this%default_fe_formulation), choices='cG,dG', error=error)
    check(error==0)  
    
  end subroutine test_poisson_add_to_cli
  
  subroutine test_poisson_parse(this)
    implicit none
    class(test_poisson_params_t), intent(inout) :: this
    integer(ip) :: istat
    
    call this%cli%parse(error=istat); check(istat==0)
    
    ! IO parameters
    call this%cli%get(switch='-d'  ,val=this%dir_path    ,error=istat); check(istat==0)
    call this%cli%get(switch='-p' ,val=this%prefix      ,error=istat); check(istat==0)
    call this%cli%get(switch='-o',val=this%dir_path_out,error=istat); check(istat==0)
    call this%cli%get(switch='-f',val=this%fe_formulation,error=istat); check(istat==0)
  end subroutine test_poisson_parse  

  subroutine test_poisson_free(this)
    implicit none
    class(test_poisson_params_t), intent(inout) :: this
    if(allocated(this%default_dir_path)) deallocate(this%default_dir_path)              
    if(allocated(this%default_prefix)) deallocate(this%default_prefix)                    
    if(allocated(this%default_dir_path_out)) deallocate(this%default_dir_path_out)
    call this%cli%free()
  end subroutine test_poisson_free

  ! GETTERS *****************************************************************************************
  function get_dir_path(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    character(len=256) :: get_dir_path
    get_dir_path = this%dir_path
  end function get_dir_path

  !==================================================================================================
  function get_prefix(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    character(len=256) :: get_prefix
    get_prefix = this%prefix
  end function get_prefix

  !==================================================================================================
  function get_dir_path_out(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    character(len=256) :: get_dir_path_out
    get_dir_path_out = this%dir_path_out
  end function get_dir_path_out
  
  !==================================================================================================
  function get_fe_formulation(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    character(len=256) :: get_fe_formulation
    get_fe_formulation = this%fe_formulation
  end function get_fe_formulation
  
end module test_poisson_params_names
