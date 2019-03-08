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

module base_direct_solver_names
    ! Serial modules
    USE types_names
    USE memor_names
    USE sparse_matrix_names
    USE serial_scalar_array_names
    USE FPL

    implicit none
# include "debug.i90"
  
    private

    integer(ip), public, parameter :: BASE_DIRECT_SOLVER_STATE_START    = 0
    integer(ip), public, parameter :: BASE_DIRECT_SOLVER_STATE_SYMBOLIC = 1 ! Symbolic data already computed
    integer(ip), public, parameter :: BASE_DIRECT_SOLVER_STATE_NUMERIC  = 2 ! Numerical data already computed

  !-----------------------------------------------------------------
  ! State transition diagram for type(base_direct_solver_t)
  !-----------------------------------------------------------------
  ! Note: all subclasses must implement this state diagram
  !-----------------------------------------------------------------
  ! Input State         | Action                | Output State 
  !-----------------------------------------------------------------
  ! Start               | symbolic_setup        | Symbolic
  ! Start               | numerical_setup       | Numeric          ! perform symbolic_setup()
  ! Start               | solve                 | Numeric          ! perform numerical_setup()
  ! Start               | free_clean            | Start
  ! Start               | free_symbolic         | Start            ! it does nothing
  ! Start               | free_numeric          | Start            ! it does nothing
  ! Start               | set_matrix            | Start
  ! Start               | update_matrix + *     | Start 

  ! Symbolic            | symbolic_setup                        | Symbolic         ! it does nothing
  ! Symbolic            | numerical_setup                       | Numeric           
  ! Symbolic            | solve                                 | Numeric          ! perform numerical_setup()
  ! Symbolic            | free_clean                            | Start
  ! Symbolic            | free_symbolic                         | Start
  ! Symbolic            | free_numeric                          | Symbolic         ! it does nothing
  ! Symbolic            | update_matrix + same_nonzero_pattern  | Symbolic         ! Re-assigns matrix pointer
  ! Symbolic            | update_matrix + !same_nonzero_pattern | Start            ! Re-assigns matrix pointer + performs free_symbolic()
    
    
  ! Numeric             | symbolic_setup                        | Numeric          ! it does nothing
  ! Numeric             | numeric_setup                         | Numeric          ! performs numerical_setup() iif numerical_setup_pending
  ! Numeric             | solve                                 | Numeric          ! performs numerical_setup() iif numerical_setup_pending
  ! Numeric             | free_numeric                          | Symbolic
  ! Numeric             | free_symbolic                         | Start
  ! Numeric             | free_clean                            | Start
  ! Numeric             | update_matrix + same_nonzero_pattern  | Numeric          ! Re-assigns matrix pointer + activates numerical_setup_pending 
  ! Numeric             | update_matrix + !same_nonzero_pattern | Start            ! Re-assigns matrix pointer + performs free_symbolic()  

    ! A relative tolerance used in evaluate_precision TBPs of base_direct_solver_t
    ! In particular, we consider the solution provided by the direct solver to be
    ! "acceptable" if ||b-Ax*||_inf / ||x*||_inf < evaluate_precision_required_precision.
    ! Obviously, this is only checked in DEBUG mode.
    real(rp), parameter :: evaluate_precision_required_precision = 1.0e-10_rp
    public :: evaluate_precision_required_precision

    type, abstract :: base_direct_solver_t
    private
        character(len=:), allocatable          :: name
        integer(ip)                            :: state  = BASE_DIRECT_SOLVER_STATE_START
        type(sparse_matrix_t), public, pointer :: matrix => NULL()
        ! Direct solvers info
        integer(ip)                            :: mem_peak_symb
        integer(ip)                            :: mem_perm_symb
        integer(ip)                            :: nz_factors
        integer(ip)                            :: mem_peak_num 
        real(rp)                               :: Mflops
    contains
    private
        procedure(base_direct_solver_free_clean_body),         public, deferred :: free_clean_body
        procedure(base_direct_solver_free_symbolic_body),      public, deferred :: free_symbolic_body
        procedure(base_direct_solver_free_numerical_body),     public, deferred :: free_numerical_body
        procedure(base_direct_solver_symbolic_setup_body),     public, deferred :: symbolic_setup_body
        procedure(base_direct_solver_numerical_setup_body),    public, deferred :: numerical_setup_body
        procedure(base_direct_solver_solve_single_rhs_body),   public, deferred :: solve_single_rhs_body
        procedure(base_direct_solver_solve_several_rhs_body),  public, deferred :: solve_several_rhs_body
        procedure(base_direct_solver_set_parameters_from_pl),  public, deferred :: set_parameters_from_pl
        procedure, non_overridable, public :: free_clean                   => base_direct_solver_free_clean
        procedure, non_overridable, public :: free_symbolic                => base_direct_solver_free_symbolic
        procedure, non_overridable, public :: free_numerical               => base_direct_solver_free_numerical
        procedure, non_overridable, public :: symbolic_setup               => base_direct_solver_symbolic_setup
        procedure, non_overridable, public :: numerical_setup              => base_direct_solver_numerical_setup
        procedure, non_overridable         :: solve_single_rhs             => base_direct_solver_solve_single_rhs
        procedure, non_overridable         :: solve_several_rhs            => base_direct_solver_solve_several_rhs
        procedure                          :: evaluate_precision_single_rhs  => base_direct_solver_evaluate_precision_single_rhs 
        procedure                          :: evaluate_precision_several_rhs => base_direct_solver_evaluate_precision_several_rhs 
        procedure, non_overridable, public :: reset                        => base_direct_solver_reset
        procedure, non_overridable, public :: set_name                     => base_direct_solver_set_name
        procedure, non_overridable, public :: set_matrix                   => base_direct_solver_set_matrix
        procedure, non_overridable, public :: get_matrix                   => base_direct_solver_get_matrix
        procedure, non_overridable, public :: replace_matrix               => base_direct_solver_replace_matrix
        procedure, non_overridable, public :: update_matrix                => base_direct_solver_update_matrix
        procedure, non_overridable, public :: matrix_is_set                => base_direct_solver_matrix_is_set
        procedure, non_overridable, public :: set_state_start              => base_direct_solver_set_state_start
        procedure, non_overridable, public :: set_state_symbolic           => base_direct_solver_set_state_symbolic
        procedure, non_overridable, public :: set_state_numeric            => base_direct_solver_set_state_numeric
        procedure, non_overridable, public :: state_is_start               => base_direct_solver_state_is_start
        procedure, non_overridable, public :: state_is_symbolic            => base_direct_solver_state_is_symbolic
        procedure, non_overridable, public :: state_is_numeric             => base_direct_solver_state_is_numeric
        procedure, non_overridable, public :: set_mem_peak_symb            => base_direct_solver_set_mem_peak_symb
        procedure, non_overridable, public :: set_mem_perm_symb            => base_direct_solver_set_mem_perm_symb
        procedure, non_overridable, public :: set_mem_peak_num             => base_direct_solver_set_mem_peak_num
        procedure, non_overridable, public :: set_nz_factors               => base_direct_solver_set_nz_factors
        procedure, non_overridable, public :: set_Mflops                   => base_direct_solver_set_Mflops
        procedure, non_overridable, public :: get_mem_peak_symb            => base_direct_solver_get_mem_peak_symb
        procedure, non_overridable, public :: get_mem_perm_symb            => base_direct_solver_get_mem_perm_symb
        procedure, non_overridable, public :: get_mem_peak_num             => base_direct_solver_get_mem_peak_num
        procedure, non_overridable, public :: get_nz_factors               => base_direct_solver_get_nz_factors
        procedure, non_overridable, public :: get_Mflops                   => base_direct_solver_get_Mflops
        procedure, non_overridable, public :: log_info                     => base_direct_solver_log_info
        generic,                    public :: solve                        => solve_single_rhs, solve_several_rhs
        generic,                    public :: evaluate_precision           => evaluate_precision_single_rhs, evaluate_precision_several_rhs
    end type

    abstract interface
        subroutine create_direct_solver_interface(base_direct_solver)
            import base_direct_solver_t
            class(base_direct_solver_t), pointer, intent(inout) :: base_direct_solver
        end subroutine

        subroutine base_direct_solver_set_parameters_from_pl(this, parameter_list)
            import base_direct_solver_t
            import ParameterList_t
            class(base_direct_solver_t),   intent(inout) :: this
            type(ParameterList_t),         intent(in)    :: parameter_list
        end subroutine base_direct_solver_set_parameters_from_pl

        ! This function returns .true. if the subclass implementor could actually
        ! perform the symbolic factorization provided the status of the input matrix to be 
        ! factorized. Some sparse direct solvers may (depending on the parameter values
        ! selected by the user) require the entries of the matrix to be factorized
        ! to be available at this stage. In such a case, this method would return .false.
        ! if the matrix entries are not avaiable at this stage, deferring the actual 
        ! symbolic factorization to the numerical factorization stage (where the entries
        ! of the matrix MUST be actually ready) 
        function base_direct_solver_symbolic_setup_body(this)
            import base_direct_solver_t
            class(base_direct_solver_t), intent(inout) :: this
            logical :: base_direct_solver_symbolic_setup_body
        end function base_direct_solver_symbolic_setup_body

        subroutine base_direct_solver_numerical_setup_body(this)
            import base_direct_solver_t
            class(base_direct_solver_t), intent(inout) :: this
        end subroutine base_direct_solver_numerical_setup_body

        subroutine base_direct_solver_solve_single_rhs_body(op, x, y)
            import base_direct_solver_t
            import serial_scalar_array_t
            class(base_direct_solver_t), intent(inout) :: op
            type(serial_scalar_array_t), intent(in)    :: x
            type(serial_scalar_array_t), intent(inout) :: y
        end subroutine base_direct_solver_solve_single_rhs_body

        subroutine base_direct_solver_solve_several_rhs_body(op, x, y)
            import base_direct_solver_t
            import ip
            import rp
            class(base_direct_solver_t), intent(inout) :: op
            real(rp),                    intent(inout) :: x(:, :)
            real(rp),                    intent(inout) :: y(:, :)
        end subroutine base_direct_solver_solve_several_rhs_body

        subroutine  base_direct_solver_free_clean_body(this)
            import base_direct_solver_t
            class(base_direct_solver_t), intent(inout) :: this
        end subroutine base_direct_solver_free_clean_body

        subroutine  base_direct_solver_free_symbolic_body(this)
            import base_direct_solver_t
            class(base_direct_solver_t), intent(inout) :: this
        end subroutine base_direct_solver_free_symbolic_body

        subroutine  base_direct_solver_free_numerical_body(this)
            import base_direct_solver_t
            class(base_direct_solver_t), intent(inout) :: this
        end subroutine base_direct_solver_free_numerical_body
    end interface

public :: base_direct_solver_t, create_direct_solver_interface

contains

    subroutine base_direct_solver_symbolic_setup(this)
        class(base_direct_solver_t), intent(inout) :: this
        ! Check pre-conditions
        assert(this%matrix_is_set())
        if(this%state_is_symbolic() .or. this%state_is_numeric()) return
        if (this%symbolic_setup_body()) call this%set_state_symbolic()
    end subroutine base_direct_solver_symbolic_setup

    subroutine base_direct_solver_numerical_setup(this)
        class(base_direct_solver_t), intent(inout) :: this
        ! Check pre-conditions
        if(this%state_is_start()) call this%symbolic_setup()
        if(this%state_is_symbolic()) then
          call this%numerical_setup_body()
          ! post-conditions
          call this%set_state_numeric()
        end if   
    end subroutine base_direct_solver_numerical_setup

    subroutine base_direct_solver_solve_single_rhs(op, x, y)
        class(base_direct_solver_t),  intent(inout) :: op
        type(serial_scalar_array_t), intent(in)    :: x
        type(serial_scalar_array_t), intent(inout) :: y
        ! Check pre-conditions
        if(.not. op%state_is_numeric()) call op%numerical_setup()
        call x%GuardTemp()
        call op%solve_single_rhs_body(x, y)
        call x%CleanTemp()
        ! post-conditions
    end subroutine base_direct_solver_solve_single_rhs
    
    subroutine base_direct_solver_solve_several_rhs(op, x, y)
        class(base_direct_solver_t),  intent(inout) :: op
        real(rp),                     intent(inout) :: x(:, :)
        real(rp),                     intent(inout) :: y(:, :)
        ! Check pre-conditions
        if(.not. op%state_is_numeric()) call op%numerical_setup()
        call op%solve_several_rhs_body(x, y)
        ! post-conditions
    end subroutine base_direct_solver_solve_several_rhs

    subroutine base_direct_solver_evaluate_precision_single_rhs(this, x, y)
      class(base_direct_solver_t),  intent(inout) :: this
      type(serial_scalar_array_t), intent(in)     :: x
      type(serial_scalar_array_t), intent(in)     :: y

      type(sparse_matrix_t),  pointer :: matrix
      type(serial_scalar_array_t)     :: r
      real(rp),               pointer :: r_real(:)
      real(rp),               pointer :: x_real(:)
      real(rp)                        :: err
      character(len=10)               :: serr
      matrix => this%get_matrix()
      call r%clone(x)
      call r%scal(-1.0,x)
      call matrix%apply_add(y,r)
      r_real => r%get_entries()
      x_real => x%get_entries()
      err = 0.0
      if (maxval(abs(x_real))>0) then
         err = maxval(abs(r_real))/maxval(abs(x_real))
      end if
      write(serr,'(e10.3)') err
      wassert( maxval(abs(r_real))<=evaluate_precision_required_precision*maxval(abs(x_real)), 'Direct solver (single rhs): returned solution is not accurate. |b-Ax|_inf/|b|_inf = '//serr )
      call r%free()
    end subroutine base_direct_solver_evaluate_precision_single_rhs
 
    subroutine base_direct_solver_evaluate_precision_several_rhs(this, x, y)
      class(base_direct_solver_t),  intent(inout) :: this
      real(rp),                     intent(in)    :: x(:, :)
      real(rp),                     intent(in)    :: y(:, :)

      type(sparse_matrix_t),  pointer                :: matrix
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
      matrix => this%get_matrix()
      nrhs = size(x,2)
      call xarr%create_and_allocate(size(x,1))
      call yarr%clone(xarr)
      call rarr%clone(xarr)
      x_real =>xarr%get_entries()
      y_real =>yarr%get_entries()
      r_real =>rarr%get_entries()
      do i =1,nrhs
         x_real(:) = x(:,i)
         y_real(:) = y(:,i)
         call rarr%scal(-1.0,xarr)
         call matrix%apply_add(yarr,rarr)
         err = 0.0
         if (maxval(abs(x_real))>0) then
            err = maxval(abs(r_real))/maxval(abs(x_real))
         end if
         write(serr,'(e10.3)') err
         wassert( maxval(abs(r_real))<=evaluate_precision_required_precision*maxval(abs(x_real)),'Direct solver (several rhs): returned solution is not accurate. |b-Ax|_inf/|b|_inf = '//serr )
      end do
      call xarr%free()
      call yarr%free()
      call rarr%free()
    end subroutine base_direct_solver_evaluate_precision_several_rhs

    subroutine  base_direct_solver_free_clean(this)
        class(base_direct_solver_t), intent(inout) :: this
        ! Check pre-conditions
        if(this%state_is_symbolic() .or. this%state_is_numeric()) call this%free_symbolic()
        call this%free_clean_body()
        ! post-conditions
        call this%set_state_start()
    end subroutine base_direct_solver_free_clean

    subroutine  base_direct_solver_free_symbolic(this)
        class(base_direct_solver_t), intent(inout) :: this
        ! Check pre-conditions
        if(this%state_is_start()) return
        if(this%state_is_numeric()) call this%free_numerical()
        call this%free_symbolic_body()
        ! post-conditions
        call this%set_mem_peak_symb(0_ip)
        call this%set_mem_perm_symb(0_ip)
        call this%set_nz_factors(0_ip)
        call this%set_state_start()
    end subroutine base_direct_solver_free_symbolic

    subroutine  base_direct_solver_free_numerical(this)
        class(base_direct_solver_t), intent(inout) :: this
                
        ! Check pre-conditions
        if(.not. this%state_is_numeric()) return
        call this%free_numerical_body
        ! post-conditions
        call this%set_mem_peak_num(0_ip)
        call this%set_Mflops(0._rp)
        call this%set_state_symbolic()
    end subroutine base_direct_solver_free_numerical

    subroutine base_direct_solver_reset(this)
        class(base_direct_solver_t), intent(inout) :: this
        this%state  = BASE_DIRECT_SOLVER_STATE_START
        nullify(this%matrix)
        this%mem_peak_symb = 0
        this%mem_perm_symb = 0
        this%nz_factors    = 0
        this%mem_peak_num  = 0
        this%Mflops        = 0._rp
    end subroutine base_direct_solver_reset

    subroutine base_direct_solver_set_name(this, name)
        class(base_direct_solver_t),   intent(inout) :: this
        character(len=*),              intent(in)    :: name
        assert(this%state_is_start())
        this%name = name
    end subroutine base_direct_solver_set_name

    subroutine base_direct_solver_set_matrix(this, matrix)
        class(base_direct_solver_t),   intent(inout) :: this
        type(sparse_matrix_t), target, intent(in)    :: matrix
        assert(.not. this%state_is_symbolic())
        assert(.not. this%state_is_numeric())
        this%matrix => matrix
    end subroutine base_direct_solver_set_matrix

    function base_direct_solver_get_matrix(this) result(matrix)
        class(base_direct_solver_t),    intent(in) :: this
        type(sparse_matrix_t), pointer             :: matrix
        assert(this%matrix_is_set())
        matrix => this%matrix
    end function base_direct_solver_get_matrix

    subroutine base_direct_solver_replace_matrix(this, matrix, same_nonzero_pattern)
    !-----------------------------------------------------------------
    !< Update matrix pointer 
    !< If same_nonzero_pattern (re-)performs numerical factorization
    !< If not same_nonzero_pattern frees its contents and forces numerical factorization
    !-----------------------------------------------------------------
        class(base_direct_solver_t),   intent(inout) :: this
        type(sparse_matrix_t), target, intent(in)    :: matrix
        logical,                       intent(in)    :: same_nonzero_pattern
    !-----------------------------------------------------------------
        this%matrix => matrix
        if(same_nonzero_pattern .and. (this%state == BASE_DIRECT_SOLVER_STATE_NUMERIC .or. &
                                       this%state == BASE_DIRECT_SOLVER_STATE_SYMBOLIC)) then
            call this%numerical_setup_body()
            call this%set_state_numeric()
        else
            call this%free_symbolic()
            call this%numerical_setup()
        endif
    end subroutine base_direct_solver_replace_matrix
    
    
    subroutine base_direct_solver_update_matrix(this, same_nonzero_pattern)
    !-----------------------------------------------------------------
    !< Signals that the matrix it was associated to has changed
    !< If same_nonzero_pattern (re-)performs numerical factorization
    !< If not same_nonzero_pattern frees its contents and forces numerical factorization
    !-----------------------------------------------------------------
        class(base_direct_solver_t),   intent(inout) :: this
        logical,                       intent(in)    :: same_nonzero_pattern
    !-----------------------------------------------------------------
        if(same_nonzero_pattern .and. (this%state == BASE_DIRECT_SOLVER_STATE_NUMERIC .or. &
                                       this%state == BASE_DIRECT_SOLVER_STATE_SYMBOLIC)) then
            call this%numerical_setup_body()
            call this%set_state_numeric()
        else
            call this%free_symbolic()
            call this%numerical_setup()
        endif
    end subroutine base_direct_solver_update_matrix

    function base_direct_solver_matrix_is_set(this) result(matrix_is_set)
        class(base_direct_solver_t), intent(in) :: this
        logical                                 :: matrix_is_set
        matrix_is_set = associated(this%matrix)
    end function base_direct_solver_matrix_is_set

    subroutine base_direct_solver_set_state_start(this)
        class(base_direct_solver_t), intent(inout) :: this
        this%state = BASE_DIRECT_SOLVER_STATE_START
    end subroutine base_direct_solver_set_state_start

    subroutine base_direct_solver_set_state_symbolic(this)
        class(base_direct_solver_t), intent(inout) :: this
        this%state = BASE_DIRECT_SOLVER_STATE_SYMBOLIC
    end subroutine base_direct_solver_set_state_symbolic

    subroutine base_direct_solver_set_state_numeric(this)
        class(base_direct_solver_t), intent(inout) :: this
        this%state = BASE_DIRECT_SOLVER_STATE_NUMERIC
    end subroutine base_direct_solver_set_state_numeric

    function base_direct_solver_state_is_start(this) result(is_start)
        class(base_direct_solver_t), intent(in) :: this
        logical                                 :: is_start
        is_start = this%state == BASE_DIRECT_SOLVER_STATE_START
    end function base_direct_solver_state_is_start

    function base_direct_solver_state_is_symbolic(this) result(is_symbolic_setup)
        class(base_direct_solver_t), intent(in) :: this
        logical                                 :: is_symbolic_setup
        is_symbolic_setup = this%state == BASE_DIRECT_SOLVER_STATE_SYMBOLIC
    end function base_direct_solver_state_is_symbolic

    function base_direct_solver_state_is_numeric(this) result(is_numerical_setup)
        class(base_direct_solver_t), intent(in) :: this
        logical                                 :: is_numerical_setup
        is_numerical_setup= this%state == BASE_DIRECT_SOLVER_STATE_NUMERIC
    end function base_direct_solver_state_is_numeric

    subroutine base_direct_solver_set_mem_peak_symb(this, mem_peak_symb)
        class(base_direct_solver_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: mem_peak_symb
        this%mem_peak_symb = mem_peak_symb
    end subroutine base_direct_solver_set_mem_peak_symb

    subroutine base_direct_solver_set_mem_perm_symb(this, mem_perm_symb)
        class(base_direct_solver_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: mem_perm_symb
        this%mem_perm_symb = mem_perm_symb
    end subroutine base_direct_solver_set_mem_perm_symb

    subroutine base_direct_solver_set_mem_peak_num(this, mem_peak_num)
        class(base_direct_solver_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: mem_peak_num
        this%mem_peak_num = mem_peak_num
    end subroutine base_direct_solver_set_mem_peak_num

    subroutine base_direct_solver_set_nz_factors(this, nz_factors)
        class(base_direct_solver_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: nz_factors
        this%nz_factors = nz_factors
    end subroutine base_direct_solver_set_nz_factors

    subroutine base_direct_solver_set_Mflops(this, Mflops)
        class(base_direct_solver_t), intent(inout) :: this
        real(rp),                    intent(in)    :: Mflops
        this%Mflops = Mflops
    end subroutine base_direct_solver_set_Mflops

    function base_direct_solver_get_mem_peak_symb(this) result(mem_peak_symb)
        class(base_direct_solver_t), intent(in) :: this
        integer(ip)                             :: mem_peak_symb
        mem_peak_symb = this%mem_peak_symb
    end function base_direct_solver_get_mem_peak_symb

    function base_direct_solver_get_mem_perm_symb(this) result(mem_perm_symb)
        class(base_direct_solver_t), intent(in) :: this
        integer(ip)                             :: mem_perm_symb
        mem_perm_symb = this%mem_perm_symb
    end function base_direct_solver_get_mem_perm_symb

    function base_direct_solver_get_mem_peak_num(this) result(mem_peak_num)
        class(base_direct_solver_t), intent(in) :: this
        integer(ip)                             :: mem_peak_num
        mem_peak_num = this%mem_peak_num
    end function base_direct_solver_get_mem_peak_num

    function base_direct_solver_get_nz_factors(this) result(nz_factors)
        class(base_direct_solver_t), intent(in) :: this
        integer(ip)                             :: nz_factors
        nz_factors = this%nz_factors
    end function base_direct_solver_get_nz_factors

    function base_direct_solver_get_Mflops(this) result(Mflops)
        class(base_direct_solver_t), intent(in) :: this
        real(rp)                                :: Mflops
        Mflops = this%Mflops
    end function base_direct_solver_get_Mflops


    subroutine base_direct_solver_log_info(this)
    !-----------------------------------------------------------------
    !< Print info
    !-----------------------------------------------------------------
        class(base_direct_solver_t), intent(in) :: this
    !-----------------------------------------------------------------
        write (*,'(a)') ''
        write (*,'(a)') '---------------------------------------'
        write (*,'(a)') ' '//this%name//' DIRECT SOLVER LOG INFO'
        write (*,'(a)') '---------------------------------------'
        write (*,'(a,i10)') 'Peak mem.      in KBytes (symb fact) = ', this%get_mem_peak_symb()
        write (*,'(a,i10)') 'Permanent mem. in KBytes (symb fact) = ', this%get_mem_perm_symb()
        write (*,'(a,i10)') 'Size of factors (thousands)          = ', this%get_nz_factors()
        write (*,'(a,i10)') 'Peak mem.      in KBytes (num fact)  = ', this%get_mem_peak_num()
        write (*,'(a,f10.2)') 'MFlops for factorization             = ', this%get_Mflops() 
        write (*,'(a)') ''
    end subroutine base_direct_solver_log_info

end module base_direct_solver_names
