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
module linear_solver_names
  use types_names
  use stdio_names
  
  ! Concrete modules
  use richardson_names
  use cg_names
  
  ! Abstract modules
  use vector_names
  use operator_names
  use base_linear_solver_names
  use environment_names
  
  implicit none
# include "debug.i90"
  private

  ! Set of possible values for the %state variable member
  integer(ip), parameter :: not_created         = 0  ! Fresh void object
  integer(ip), parameter :: environment_set     = 1  ! Pointer to environment set
  integer(ip), parameter :: solver_type_set     = 2  ! Dynamic type of solver set
  
  ! State transition diagram for type(linear_solver_t)
  ! --------------------------------------------------------
  ! Input State      | Action               | Output State 
  ! --------------------------------------------------------
  ! not_created      | create               | environment_set
  ! not_created      | free                 | not_created
  ! environment_set  | set_type_from_pl     | solver_type_set
  ! environment_set  | free                 | not_created
  ! solver_type_set  | set_type_from_pl     | solver_type_set
  ! solver_type_set  | free                 | not_created
  type :: linear_solver_t
    private
    class(environment_t)       , pointer  :: environment
    class(base_linear_solver_t), pointer  :: base_linear_solver
    integer(ip)                           :: state = not_created
  contains
    ! Concrete TBPs
    procedure :: create                          => linear_solver_create
    procedure :: free                            => linear_solver_free
    procedure :: solve                           => linear_solver_solve
    procedure :: set_type_from_pl                => linear_solver_set_type_from_pl
    procedure :: set_parameters_from_pl          => linear_solver_set_parameters_from_pl
    procedure :: set_type_and_parameters_from_pl => linear_solver_set_type_and_parameters_from_pl
    procedure :: set_operators                   => linear_solver_set_operators
    procedure :: set_initial_solution            => linear_solver_set_initial_solution
    procedure :: set_rhs                         => linear_solver_set_rhs
  end type
  
  ! Data types
  public :: linear_solver_t
  
contains    
   subroutine linear_solver_create ( this, environment )
     implicit none
     class(linear_solver_t)      , intent(inout) :: this
     class(environment_t), target, intent(in)    :: environment
     assert ( this%state == not_created )
     this%environment => environment
     this%state = environment_set
   end subroutine linear_solver_create
   
   subroutine linear_solver_free ( this )
     implicit none
     class(linear_solver_t), intent(inout) :: this
     if ( this%state == solver_type_set ) then
       call this%base_linear_solver%free()
       deallocate (this%base_linear_solver)
     end if
     nullify(this%environment)
     this%state = not_created
   end subroutine linear_solver_free
   
   subroutine linear_solver_solve ( this, x )
     implicit none
     class(linear_solver_t), intent(inout) :: this
     class(vector_t)       , intent(inout) :: x 
     assert ( this%state == solver_type_set )
     call this%base_linear_solver%solve(x)
   end subroutine linear_solver_solve

   subroutine linear_solver_set_type_from_pl ( this )
     implicit none
     class(linear_solver_t), intent(inout) :: this
     character(len=:)      , allocatable   :: linear_solver_type
     
     assert ( this%state == environment_set .or. this%state == solver_type_set )
     
     if ( this%state == solver_type_set ) then
       ! PENDING: ONLY FREE IF THE TYPE SELECTED DOES NOT MATCH THE EXISTING ONE
       call this%base_linear_solver%free()
       deallocate ( this%base_linear_solver )
     end if
     
     ! PENDING
     ! 1. Get val associated to key="linear_solver_type" from type(ParameterList)
     ! 2. Select Factory Method associated to val from "global" (and dynamically built) dictionary of Factory Methods
     ! 3. Only create if this%state == environment_set or if base_linear_solver was freed in the block of code above
     !this%base_linear_solver => create_richardson(this%environment)
     this%base_linear_solver => create_cg(this%environment)
     assert ( this%base_linear_solver%get_state() == start )
     this%state = solver_type_set
   end subroutine linear_solver_set_type_from_pl
   
   subroutine linear_solver_set_parameters_from_pl ( this )
     implicit none
     class(linear_solver_t), intent(inout) :: this
     character(len=:)      , allocatable   :: linear_solver_type
     assert ( this%state == solver_type_set )
     call this%base_linear_solver%set_parameters_from_pl()
   end subroutine linear_solver_set_parameters_from_pl
   
   subroutine linear_solver_set_type_and_parameters_from_pl ( this )
     implicit none
     class(linear_solver_t), intent(inout) :: this
     call this%set_type_from_pl()
     call this%set_parameters_from_pl()
   end subroutine linear_solver_set_type_and_parameters_from_pl
   
   subroutine linear_solver_set_operators ( this, A, M )
     implicit none
     class(linear_solver_t), intent(inout) :: this
     class(operator_t)     , intent(in)    :: A, M
     assert ( this%state == solver_type_set )
     call this%base_linear_solver%set_operators(A,M)
   end subroutine linear_solver_set_operators
   
   subroutine linear_solver_set_rhs ( this, b )
     implicit none
     class(linear_solver_t), intent(inout) :: this
     class(vector_t)       , intent(in)    :: b
     assert ( this%state == solver_type_set )
     call this%base_linear_solver%set_rhs(b)
   end subroutine linear_solver_set_rhs
   
   subroutine linear_solver_set_initial_solution( this, initial_solution )
     implicit none
     class(linear_solver_t), intent(inout) :: this
     class(vector_t)       , intent(in)    :: initial_solution
     assert ( this%state == solver_type_set )
     call this%base_linear_solver%set_initial_solution(initial_solution)
   end subroutine linear_solver_set_initial_solution
   
end module linear_solver_names
