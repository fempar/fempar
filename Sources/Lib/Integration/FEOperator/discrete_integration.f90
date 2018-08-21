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
module discrete_integration_names
  use field_names
  use reference_fe_names
  use types_names
  use assembler_names
  use fe_space_names
  use memor_names
  use vector_names
  use serial_scalar_array_names

  implicit none
# include "debug.i90"
  private

  
  ! sbadia :  I would keep a pointer to the FE spaces in discrete_integration, and only one integrate
  type, abstract :: discrete_integration_t
   contains
     procedure (integrate_galerkin_interface), deferred :: integrate_galerkin
     procedure                                          :: integrate_petrov_galerkin
     procedure                                          :: integrate_residual
     procedure                                          :: integrate_petrov_galerkin_residual
     procedure                                          :: integrate_tangent
     procedure                                          :: integrate_petrov_galerkin_tangent
     procedure                                          :: set_evaluation_point
     procedure                                          :: set_boundary_data
     procedure                                          :: set_current_time
     procedure                                          :: get_boundary_data
     generic   :: integrate => integrate_galerkin, integrate_petrov_galerkin
  end type discrete_integration_t

  type p_discrete_integration_t
     class(discrete_integration_t), pointer :: p => NULL()  
  end type p_discrete_integration_t

  public :: discrete_integration_t, p_discrete_integration_t

  abstract interface
     subroutine integrate_galerkin_interface ( this, fe_space, assembler )
       import :: discrete_integration_t, serial_fe_space_t, assembler_t
       implicit none
       class(discrete_integration_t)  ,    intent(in)    :: this
       class(serial_fe_space_t)       ,    intent(inout) :: fe_space
       class(assembler_t),    intent(inout) :: assembler
     end subroutine integrate_galerkin_interface
  end interface
    
  contains
     subroutine integrate_petrov_galerkin ( this, fe_space_trial, fe_space_test, assembler )
       implicit none
       class(discrete_integration_t)  ,    intent(in)    :: this
       class(serial_fe_space_t)       ,    intent(inout) :: fe_space_trial
       class(serial_fe_space_t)       ,    intent(inout) :: fe_space_test
       class(assembler_t),    intent(inout) :: assembler
       mcheck(.false.,"You must implement integrate Petrov-Galerkin if you want to use it")       
     end subroutine  integrate_petrov_galerkin 
  
     subroutine integrate_residual ( this, fe_space, assembler )
       implicit none
       class(discrete_integration_t)  ,    intent(in)    :: this
       class(serial_fe_space_t)       ,    intent(inout) :: fe_space
       class(assembler_t),    intent(inout) :: assembler
       mcheck(.false.,"You must implement integrate Galerkin residual if you want to use it")       
     end subroutine  integrate_residual
  
  subroutine integrate_petrov_galerkin_residual ( this, fe_space_trial, fe_space_test, assembler )
       implicit none
       class(discrete_integration_t)  ,    intent(in)    :: this
       class(serial_fe_space_t)       ,    intent(inout) :: fe_space_trial
       class(serial_fe_space_t)       ,    intent(inout) :: fe_space_test
       class(assembler_t),    intent(inout) :: assembler
       mcheck(.false.,"You must implement integrate Petrov-Galerkin residual if you want to use it")       
     end subroutine  integrate_petrov_galerkin_residual
  
     subroutine integrate_tangent ( this, fe_space, assembler )
       implicit none
       class(discrete_integration_t)  ,    intent(in)    :: this
       class(serial_fe_space_t)       ,    intent(inout) :: fe_space
       class(assembler_t),    intent(inout) :: assembler
       mcheck(.false.,"You must implement integrate Galerkin tangent if you want to use it")       
     end subroutine  integrate_tangent
  
  subroutine integrate_petrov_galerkin_tangent ( this, fe_space_trial, fe_space_test, assembler )
       implicit none
       class(discrete_integration_t)  ,    intent(in)    :: this
       class(serial_fe_space_t)       ,    intent(inout) :: fe_space_trial
       class(serial_fe_space_t)       ,    intent(inout) :: fe_space_test
       class(assembler_t),    intent(inout) :: assembler
       mcheck(.false.,"You must implement integrate Petrov-Galerkin tangent if you want to use it")       
     end subroutine  integrate_petrov_galerkin_tangent
  
     subroutine set_evaluation_point ( this, evaluation_point )
       implicit none
       class(discrete_integration_t)  ,    intent(inout)    :: this
       class(vector_t)                ,    intent(in)       :: evaluation_point     
       mcheck(.false.,"You must implement set_evaluation_point if you want to use it")
     end subroutine  set_evaluation_point 
  
  subroutine set_boundary_data ( this, boundary_data )
       implicit none
       class(discrete_integration_t)  ,    intent(inout)    :: this
       type(serial_scalar_array_t)    ,    intent(in)       :: boundary_data     
       mcheck(.false.,"You must implement set_boundary_data if you want to use it")
     end subroutine  set_boundary_data 
          
     subroutine set_current_time ( this, fe_space, current_time)
       implicit none
       class(discrete_integration_t)  ,    intent(inout)    :: this
       class(serial_fe_space_t)       ,    intent(in)       :: fe_space
       real(rp)                       ,    intent(in)       :: current_time     
       mcheck(.false.,"You must implement set_current_time if you want to use it")
     end subroutine  set_current_time
     
     function get_boundary_data ( this )
       implicit none
       class(discrete_integration_t)  ,    intent(inout)    :: this
       type(serial_scalar_array_t), pointer :: get_boundary_data
       mcheck(.false.,"You must implement set_current_time if you want to use it")
     end function get_boundary_data
       
end module discrete_integration_names
