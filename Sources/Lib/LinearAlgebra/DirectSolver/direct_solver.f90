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
! licensors of this Program grant you additional permission to convey the 
! resulting work. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module direct_solver_names
    ! Serial modules
    USE types_names
    USE memor_names
    USE operator_names
    USE linear_solver_names
    USE vector_space_names
    USE base_direct_solver_names
    USE direct_solver_parameters_names
    USE sparse_matrix_names
    USE vector_names
    USE serial_scalar_array_names
    USE direct_solver_creational_methods_dictionary_names
    USE FPL
    
implicit none
# include "debug.i90"

private

    type, extends(linear_solver_t) :: direct_solver_t
    private
        class(base_direct_solver_t), pointer :: base_direct_solver  => NULL()
    contains
    private
        procedure, non_overridable, public :: set_type                => direct_solver_set_type
        procedure, non_overridable, public :: set_type_from_pl        => direct_solver_set_type_from_pl
        procedure, non_overridable, public :: set_parameters_from_pl  => direct_solver_set_parameters_from_pl
        procedure, non_overridable, public :: set_matrix              => direct_solver_set_matrix
        procedure, non_overridable, public :: replace_matrix          => direct_solver_replace_matrix
        procedure,                  public :: update_matrix           => direct_solver_update_matrix
        procedure, non_overridable         :: create_vector_spaces    => direct_solver_create_vector_spaces
        procedure, non_overridable, public :: symbolic_setup          => direct_solver_symbolic_setup
        procedure, non_overridable, public :: numerical_setup         => direct_solver_numerical_setup
        procedure, non_overridable, public :: log_info                => direct_solver_log_info
        procedure, non_overridable         :: solve_single_rhs        => direct_solver_solve_single_rhs
        procedure, non_overridable         :: solve_several_rhs       => direct_solver_solve_several_rhs
        procedure,                  public :: apply                   => direct_solver_apply
        procedure,                  public :: apply_add               => direct_solver_apply_add
        procedure, non_overridable, public :: free_in_stages          => direct_solver_free_in_stages
        procedure, non_overridable, public :: free                    => direct_solver_free
        generic,                    public :: solve                   => solve_single_rhs, solve_several_rhs
    end type

public :: direct_solver_t

contains

    subroutine direct_solver_set_type(this, name)
    !-----------------------------------------------------------------
    !< Allocate the concrete direct solver from a given solver name
    !-----------------------------------------------------------------
        class(direct_solver_t),                    intent(inout) :: this
        character(len=*),                          intent(in)    :: name
        procedure(create_direct_solver_interface), pointer       :: create
        logical                                                  :: isPresent
    !-----------------------------------------------------------------
        if(associated(this%base_direct_solver)) then
            call this%base_direct_solver%free_clean()
            deallocate(this%base_direct_solver)
        endif
        nullify(create)
        assert(the_direct_solver_creational_methods_dictionary%isInitialized())
        assert(the_direct_solver_creational_methods_dictionary%isPresent(Key=name))
        call the_direct_solver_creational_methods_dictionary%Get(Key=name,Proc=create)
        if(associated(create)) call create(this%base_direct_solver)
    end subroutine direct_solver_set_type


    subroutine direct_solver_set_type_from_pl(this, parameter_list)
    !-----------------------------------------------------------------
    !< Allocate the concrete direct solver from a solver name stored 
    !< in the parameter list
    !-----------------------------------------------------------------
        class(direct_solver_t), intent(inout) :: this
        type(ParameterList_t),  intent(in)    :: parameter_list
        character(len=:), allocatable         :: name
        integer                               :: DataSizeInBytes
        integer                               :: FPLError
    !-----------------------------------------------------------------
        ! check if DIRECT_SOLVER_TYPE is present and is a scalar string
        ! in the given parameter list,
        assert(parameter_list%isAssignable(dls_type_key, 'string')) 
        FPLError = parameter_list%GetAsString(Key=dls_type_key, String=name)
        assert(FPLError == 0)
        call this%set_type(name)
    end subroutine direct_solver_set_type_from_pl


    subroutine direct_solver_set_parameters_from_pl(this, parameter_list)
    !-----------------------------------------------------------------
    !< Set parameters from values stored in a given parameter list
    !-----------------------------------------------------------------
        class(direct_solver_t), intent(inout) :: this
        type(ParameterList_t),  intent(in)    :: parameter_list
    !-----------------------------------------------------------------
        assert(associated(this%base_direct_solver))
        call this%base_direct_solver%set_parameters_from_pl(parameter_list)
    end subroutine direct_solver_set_parameters_from_pl
    
    subroutine direct_solver_set_matrix(this, matrix)
    !-----------------------------------------------------------------
    !< Associate the concrete direct solver with an matrix
    !-----------------------------------------------------------------
        class(direct_solver_t),         intent(inout) :: this
        class(operator_t),              intent(in)    :: matrix
    !-----------------------------------------------------------------
        assert(associated(this%base_direct_solver))
        select type(matrix)
        class is (sparse_matrix_t)
          call this%base_direct_solver%set_matrix(matrix)
        class is (lvalue_operator_t)
          select type ( op => matrix%get_operator() ) 
          class is (sparse_matrix_t)
            call this%base_direct_solver%set_matrix(op)
          class default
            mcheck(.false., "direct_solver_set_matrix:: matrix MUST be of dynamic type sparse_matrix_t")
          end select 
        class default
          mcheck(.false., "direct_solver_set_matrix:: matrix MUST be of dynamic type sparse_matrix_t")
        end select
    end subroutine direct_solver_set_matrix
    
   subroutine direct_solver_replace_matrix(this, matrix, same_nonzero_pattern)
    !-----------------------------------------------------------------
    !< Update matrix pointer 
    !< If same_nonzero_pattern numerical_setup has to be performed
    !< If not same_nonzero_pattern symbolic_setup has to be performed
    !-----------------------------------------------------------------
        class(direct_solver_t),        intent(inout) :: this
        class(operator_t),             intent(in)    :: matrix
        logical,                       intent(in)    :: same_nonzero_pattern
    !-----------------------------------------------------------------
        assert(associated(this%base_direct_solver))
        select type(matrix)
        class is (sparse_matrix_t)
          call this%base_direct_solver%replace_matrix(matrix, same_nonzero_pattern)  
        class is (lvalue_operator_t)
          select type ( op => matrix%get_operator() ) 
          class is (sparse_matrix_t)
            call this%base_direct_solver%replace_matrix(op, same_nonzero_pattern)  
          class default
            mcheck(.false., "direct_solver_update_matrix:: matrix MUST be of dynamic type sparse_matrix_t")
          end select 
        class default
          mcheck(.false., "direct_solver_update_matrix:: matrix MUST be of dynamic type sparse_matrix_t")
        end select
        if(.not. same_nonzero_pattern) call this%free_vector_spaces()
    end subroutine direct_solver_replace_matrix
    
    subroutine direct_solver_update_matrix(this, same_nonzero_pattern)
    !-----------------------------------------------------------------
    !< Update matrix pointer 
    !< If same_nonzero_pattern numerical_setup has to be performed
    !< If not same_nonzero_pattern symbolic_setup has to be performed
    !-----------------------------------------------------------------
        class(direct_solver_t),        intent(inout) :: this
        logical,                       intent(in)    :: same_nonzero_pattern
    !-----------------------------------------------------------------
      call this%base_direct_solver%update_matrix(same_nonzero_pattern)  
      if(.not. same_nonzero_pattern) call this%free_vector_spaces()
    end subroutine direct_solver_update_matrix
    
    subroutine direct_solver_create_vector_spaces ( this ) 
    !-----------------------------------------------------------------
    !< Clone vector spaces from matrix vector spaces
    !-----------------------------------------------------------------
        class(direct_solver_t),          intent(inout) :: this
        type(sparse_matrix_t),  pointer                :: matrix
        type(vector_space_t),   pointer                :: matrix_domain_vector_space
        type(vector_space_t),   pointer                :: matrix_range_vector_space
        type(vector_space_t),   pointer                :: direct_solver_domain_vector_space
        type(vector_space_t),   pointer                :: direct_solver_range_vector_space
    !-----------------------------------------------------------------
        assert(.not. this%vector_spaces_are_created())
        assert(this%base_direct_solver%matrix_is_set())
        matrix => this%base_direct_solver%get_matrix()
        matrix_domain_vector_space        => matrix%get_domain_vector_space()
        matrix_range_vector_space         => matrix%get_range_vector_space()
        direct_solver_domain_vector_space => this%get_domain_vector_space()
        direct_solver_range_vector_space  => this%get_range_vector_space()
        call matrix_domain_vector_space%clone(direct_solver_domain_vector_space)
        call matrix_range_vector_space%clone(direct_solver_range_vector_space)
    end subroutine direct_solver_create_vector_spaces


    subroutine direct_solver_symbolic_setup(this)
    !-----------------------------------------------------------------
    !< Concrete direct solver performs analysis phase
    !< Check vector spaces status and create or update
    !< vector spaces if needed
    !-----------------------------------------------------------------
        class(direct_solver_t), intent(inout) :: this
    !-----------------------------------------------------------------
        type(sparse_matrix_t),  pointer                :: matrix
        assert(associated(this%base_direct_solver))
        assert(this%base_direct_solver%matrix_is_set())
        if(.not. this%vector_spaces_are_created()) &
            call this%create_vector_spaces()
        matrix => this%base_direct_solver%get_matrix()
        if(matrix%get_num_rows()==0 .or. matrix%get_num_cols()==0) return
        call this%base_direct_solver%symbolic_setup()
    end subroutine direct_solver_symbolic_setup


    subroutine direct_solver_numerical_setup(this)
    !-----------------------------------------------------------------
    !< Concrete direct solver performs factorization phase
    !-----------------------------------------------------------------
        class(direct_solver_t), intent(inout) :: this
    !-----------------------------------------------------------------
        type(sparse_matrix_t),  pointer                :: matrix
        assert(associated(this%base_direct_solver))
        assert(this%base_direct_solver%matrix_is_set())
        if(.not. this%vector_spaces_are_created()) &
            call this%create_vector_spaces()
        matrix => this%base_direct_solver%get_matrix()
        if(matrix%get_num_rows()==0 .or. matrix%get_num_cols()==0) return
        call this%base_direct_solver%numerical_setup()
    end subroutine direct_solver_numerical_setup


    subroutine direct_solver_log_info(this)
    !-----------------------------------------------------------------
    !< Print direct solver log info
    !-----------------------------------------------------------------
        class(direct_solver_t), intent(inout) :: this
    !-----------------------------------------------------------------
        assert(associated(this%base_direct_solver))
        call this%base_direct_solver%log_info()
    end subroutine direct_solver_log_info


    subroutine direct_solver_solve_single_rhs(op, x, y)
      !-----------------------------------------------------------------
      !< Computes y <- A^-1 * x
      !-----------------------------------------------------------------
      class(direct_solver_t)                :: op
      class(vector_t),        intent(in)    :: x
      class(vector_t),        intent(inout) :: y
      !-----------------------------------------------------------------
      type(sparse_matrix_t),  pointer                :: matrix

#ifdef DEBUG
      type(serial_scalar_array_t)     :: r
      real(rp),               pointer :: r_real(:)
      real(rp),               pointer :: x_real(:)
      real(rp)                        :: err
      character(len=10)               :: serr
      real(rp),             parameter :: tol = 1.0e-10
#endif

        call x%GuardTemp()
        assert(associated(op%base_direct_solver))
        assert(op%base_direct_solver%matrix_is_set())
        if(.not. op%vector_spaces_are_created()) &
            call op%create_vector_spaces()
        matrix => op%base_direct_solver%get_matrix()
        if(matrix%get_num_rows()==0 .or. matrix%get_num_cols()==0) return
        select type (x)
            type is (serial_scalar_array_t)
                select type (y)
                    type is (serial_scalar_array_t)
                        call op%base_direct_solver%solve(x,y)
#ifdef DEBUG
                        call op%base_direct_solver%evaluate_precision(x,y)
#endif
                class DEFAULT 
                    check(.false.) 
               end select
         class DEFAULT
            check(.false.)
        end select
        call x%CleanTemp()

    end subroutine direct_solver_solve_single_rhs

    subroutine direct_solver_solve_several_rhs(op, x, y)
      !-----------------------------------------------------------------
      !< Computes y <- A^-1 * x
      !-----------------------------------------------------------------
      class(direct_solver_t), intent(inout) :: op
      real(rp),               intent(inout) :: x(:, :)
      real(rp),               intent(inout) :: y(:, :)
      !-----------------------------------------------------------------
      type(sparse_matrix_t),  pointer                :: matrix

#ifdef DEBUG
      type(serial_scalar_array_t)     :: xarr
      type(serial_scalar_array_t)     :: yarr
      type(serial_scalar_array_t)     :: rarr
      real(rp),               pointer :: x_real(:)
      real(rp),               pointer :: y_real(:)
      real(rp),               pointer :: r_real(:)
      integer(ip)                     :: nrhs
      integer(ip)                     :: i
      real(rp)                        :: err
      character(len=10)               :: serr
      real(rp),             parameter :: tol = 1.0e-10
#endif

      assert(associated(op%base_direct_solver))
      assert(op%base_direct_solver%matrix_is_set())
      if(.not. op%vector_spaces_are_created()) &
           call op%create_vector_spaces()
      matrix => op%base_direct_solver%get_matrix()
      if(matrix%get_num_rows()==0 .or. matrix%get_num_cols()==0) return
      call op%base_direct_solver%solve(x, y)

#ifdef DEBUG
      call op%base_direct_solver%evaluate_precision(x,y) 
#endif

    end subroutine direct_solver_solve_several_rhs


    subroutine direct_solver_apply(this, x, y)
    !-----------------------------------------------------------------
    !< Call to Solve (Computes y <- A^-1 * x)
    !-----------------------------------------------------------------
        class(direct_solver_t), intent(inout)    :: this
        class(vector_t),        intent(in)    :: x
        class(vector_t),        intent(inout) :: y
    !-----------------------------------------------------------------
        call this%solve(x,y)
    end subroutine direct_solver_apply
    
    
    subroutine direct_solver_apply_add(this, x, y)
    !-----------------------------------------------------------------
    !< Call to Solve (Computes y <- A^-1 * x + y)
    !-----------------------------------------------------------------
        class(direct_solver_t), intent(inout)    :: this
        class(vector_t),        intent(in)    :: x
        class(vector_t),        intent(inout) :: y
        class(vector_t), allocatable          :: w
        type(vector_space_t), pointer         :: range_vector_space
        integer(ip)                           :: istat
    !-----------------------------------------------------------------
        range_vector_space => this%get_range_vector_space()
        call range_vector_space%create_vector(w)
        call this%solve(x,w)
        call y%axpby(1.0, w, 1.0)
        call w%free()
        deallocate(w, stat=istat); check(istat==0)
    end subroutine direct_solver_apply_add
    
    subroutine direct_solver_free_in_stages(this, action)
    !-----------------------------------------------------------------
    !< Free direct solver in stages
    !-----------------------------------------------------------------
        class(direct_solver_t), intent(inout) :: this
        integer(ip),            intent(in)    :: action
    !-----------------------------------------------------------------
        if(associated(this%base_direct_solver)) then
            select case (action)
                case (free_numerical_setup)
                    call this%base_direct_solver%free_numerical()
                case (free_symbolic_setup)
                    call this%base_direct_solver%free_symbolic()
                case (free_clean)
                    call this%free_vector_spaces()
                    call this%base_direct_solver%free_clean()
                    deallocate(this%base_direct_solver)
                    nullify(this%base_direct_solver)
                case DEFAULT
                    assert(.false.)
            end select
        endif
    end subroutine direct_solver_free_in_stages


    subroutine direct_solver_free(this)
    !-----------------------------------------------------------------
    !< Free direct solver
    !-----------------------------------------------------------------
        class(direct_solver_t), intent(inout) :: this
    !-----------------------------------------------------------------
        call this%free_in_stages(free_numerical_setup)
        call this%free_in_stages(free_symbolic_setup)
        call this%free_in_stages(free_clean)
    end subroutine direct_solver_free

end module direct_solver_names
