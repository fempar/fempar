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
module iterative_linear_solver_names
  use types_names
  use stdio_names
  
  ! Abstract modules
  use vector_names
  use operator_names
  use base_iterative_linear_solver_names
  use iterative_linear_solver_parameters_names
  use iterative_linear_solver_creational_methods_dictionary_names
  use environment_names
  use ParameterList
  
  implicit none
# include "debug.i90"
  private

  ! Set of possible values for the %state variable member
  integer(ip), parameter :: not_created         = 0  ! Fresh void object
  integer(ip), parameter :: environment_set     = 1  ! Pointer to environment set
  integer(ip), parameter :: solver_type_set     = 2  ! Dynamic type of solver set
  
  ! State transition diagram for type(iterative_linear_solver_t)
  ! --------------------------------------------------------
  ! Input State      | Action               | Output State 
  ! --------------------------------------------------------
  ! not_created      | create               | environment_set
  ! not_created      | free                 | not_created
  ! environment_set  | set_type_from_pl     | solver_type_set
  ! environment_set  | free                 | not_created
  ! solver_type_set  | set_type_from_pl     | solver_type_set
  ! solver_type_set  | free                 | not_created
  type :: iterative_linear_solver_t
     private
     class(environment_t)       ,           pointer  :: environment                  => NULL()
     class(base_iterative_linear_solver_t), pointer  :: base_iterative_linear_solver => NULL()
     integer(ip)                                     :: state = not_created
   contains
     ! Concrete TBPs
     procedure :: create                          => iterative_linear_solver_create
     procedure :: free                            => iterative_linear_solver_free
     procedure :: solve                           => iterative_linear_solver_solve
     procedure :: print_convergence_history       => iterative_linear_solver_print_convergence_history
     procedure :: set_type_from_pl                => iterative_linear_solver_set_type_from_pl
     procedure :: set_parameters_from_pl          => iterative_linear_solver_set_parameters_from_pl
     procedure :: set_type_and_parameters_from_pl => iterative_linear_solver_set_type_and_parameters_from_pl
     procedure :: set_operators                   => iterative_linear_solver_set_operators
     procedure :: set_initial_solution            => iterative_linear_solver_set_initial_solution
     procedure :: set_type_from_string            => iterative_linear_solver_set_type_from_string
  end type iterative_linear_solver_t
  
  public :: iterative_linear_solver_t
  
contains    
   subroutine iterative_linear_solver_create ( this, environment )
     implicit none
     class(iterative_linear_solver_t)      , intent(inout) :: this
     class(environment_t), target, intent(in)    :: environment
     assert ( this%state == not_created )
     this%environment => environment
     this%state = environment_set
   end subroutine iterative_linear_solver_create
   
   subroutine iterative_linear_solver_free ( this )
     implicit none
     class(iterative_linear_solver_t), intent(inout) :: this
     if ( this%state == solver_type_set ) then
       call this%base_iterative_linear_solver%free()
       deallocate (this%base_iterative_linear_solver)
     end if
     nullify(this%environment)
     this%state = not_created
   end subroutine iterative_linear_solver_free
   
   subroutine iterative_linear_solver_print_convergence_history ( this, file_path )
     implicit none
     class(iterative_linear_solver_t), intent(in) :: this
     character(len=*)      , intent(in) :: file_path
     assert ( this%state == solver_type_set )
     call this%base_iterative_linear_solver%print_convergence_history(file_path)
   end subroutine iterative_linear_solver_print_convergence_history
   
   subroutine iterative_linear_solver_solve ( this, b, x )
     implicit none
     class(iterative_linear_solver_t), intent(inout) :: this
     class(vector_t)       , intent(in) :: b 
     class(vector_t)       , intent(inout) :: x 
     assert ( this%state == solver_type_set )
     call this%base_iterative_linear_solver%solve(b,x)
   end subroutine iterative_linear_solver_solve

   subroutine iterative_linear_solver_set_type_from_pl ( this, parameter_list )
     implicit none
     class(iterative_linear_solver_t), intent(inout) :: this
     type(ParameterList_t),            intent(in)    :: parameter_list
     character(len=:)      , allocatable             :: iterative_linear_solver_type
     integer(ip)                                     :: DataSizeInBytes
     integer                                         :: FPLError
#ifdef DEBUG
     logical                                         :: is_present
     logical                                         :: is_string
     integer(ip), allocatable                        :: shape(:)
     is_present = parameter_list%isPresent   (Key = ils_type)
     is_string  = parameter_list%isOfDataType(Key = ils_type, mold  = 'string')
     FPLError   = parameter_list%getshape    (Key = ils_type, shape = shape)
     ! check if ITERATIVE_LINEAR_SOLVER_TYPE is present and is a scalar string
     ! in the given parameter list,
     assert(is_present .and. is_string .and. size(shape) == 0) 
#endif
     DataSizeInBytes = parameter_list%DataSizeInBytes(Key=ils_type)
     allocate(character(len=DataSizeInBytes) :: iterative_linear_solver_type, stat=FPLError)
     assert(FPLError == 0)
     FPLError = parameter_list%Get(Key=ils_type, Value=iterative_linear_solver_type)
     assert(FPLError == 0)
     call this%set_type_from_string (iterative_linear_solver_type)
   end subroutine iterative_linear_solver_set_type_from_pl
   
   subroutine iterative_linear_solver_set_parameters_from_pl ( this, parameter_list )
     implicit none
     class(iterative_linear_solver_t), intent(inout) :: this
     type(ParameterList_t),            intent(in)    :: parameter_list
     assert ( this%state == solver_type_set )
     call this%base_iterative_linear_solver%set_parameters_from_pl(parameter_list)
   end subroutine iterative_linear_solver_set_parameters_from_pl
   
   subroutine iterative_linear_solver_set_type_and_parameters_from_pl ( this, parameter_list )
     implicit none
     class(iterative_linear_solver_t), intent(inout) :: this
     type(ParameterList_t),            intent(in)    :: parameter_list
     call this%set_type_from_pl( parameter_list )
     call this%set_parameters_from_pl(parameter_list)
   end subroutine iterative_linear_solver_set_type_and_parameters_from_pl
   
   subroutine iterative_linear_solver_set_operators ( this, A, M )
     implicit none
     class(iterative_linear_solver_t), intent(inout) :: this
     class(operator_t)     , intent(in)    :: A, M
     assert ( this%state == solver_type_set )
     call this%base_iterative_linear_solver%set_operators(A,M)
   end subroutine iterative_linear_solver_set_operators
      
   subroutine iterative_linear_solver_set_initial_solution( this, initial_solution )
     implicit none
     class(iterative_linear_solver_t), intent(inout) :: this
     class(vector_t)       , intent(in)    :: initial_solution
     assert ( this%state == solver_type_set )
     call this%base_iterative_linear_solver%set_initial_solution(initial_solution)
   end subroutine iterative_linear_solver_set_initial_solution

   subroutine iterative_linear_solver_set_type_from_string (this, iterative_linear_solver_type)
     implicit none
     class(iterative_linear_solver_t),                    intent(inout) :: this
     character(len=*)                ,                    intent(in)    :: iterative_linear_solver_type
     procedure(create_iterative_linear_solver_interface), pointer       :: create

     assert ( this%state == environment_set .or. this%state == solver_type_set )
     
     if ( this%state == solver_type_set ) then
       ! PENDING: ONLY FREE IF THE TYPE SELECTED DOES NOT MATCH THE EXISTING ONE
       call this%base_iterative_linear_solver%free()
       deallocate ( this%base_iterative_linear_solver )
     end if

     nullify(create)
     assert(the_iterative_linear_solver_creational_methods_dictionary%isInitialized())
     assert(the_iterative_linear_solver_creational_methods_dictionary%isPresent(Key=iterative_linear_solver_type))
     call the_iterative_linear_solver_creational_methods_dictionary%Get(Key=iterative_linear_solver_type,Proc=create)
     if(associated(create)) call create(this%environment, this%base_iterative_linear_solver)
     
     assert ( this%base_iterative_linear_solver%get_state() == start )
     this%state = solver_type_set
   end subroutine iterative_linear_solver_set_type_from_string
   
end module iterative_linear_solver_names
