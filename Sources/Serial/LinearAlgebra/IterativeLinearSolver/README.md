# Iterative Linear Solvers: State transition diagram

```fortran
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
```

# Iterative Linear Solvers: Basic usage

```fortran
...
    type(ParameterList_t)                :: parameter_list
    type(ParameterList_t), pointer       :: iterative_linear_solver_params
    type(fe_affine_operator_t)           :: fe_affine_operator
    type(poisson_discrete_integration_t) :: poisson_integration
    type(iterative_linear_solver_t)      :: iterative_linear_solver
    type(serial_environment_t)           :: senv
    class(vector_t), allocatable         :: computed_solution_vector
...
    call FPL_Init()                                                                        !
    call the_iterative_linear_solver_creational_methods_dictionary%init()                  ! Initialization. In a FEMPAR_INIT subroutine in the future
    call parameter_list%Init()                                                             !
...
    iterative_linear_solver_params => parameter_list%NewSubList(Key=cg_name)               ! Create a new parameter sublist to store CG iterative linear solver parameters
    call iterative_linear_solver_params%set(key = ils_type, Value = cg_name)               ! Set the iterative linear solver type parameter into the list
...
     ! Create iterative linear solver, set operators and solve linear system
     call iterative_linear_solver%create(senv)                                                          ! Create iterative linear solver and set environment
     call iterative_linear_solver%set_type_and_parameters_from_pl(parameter_list)                       ! Set concrete type of iterative linear solver and it parameters
     call iterative_linear_solver%set_operators(fe_affine_operator, .identity. fe_affine_operator)      ! Set operators
     call iterative_linear_solver%solve(fe_affine_operator%get_translation(), computed_solution_vector) ! Solve
     call iterative_linear_solver%free()                                                                ! Free
...

```

# Iterative Linear Solvers: Public Parameters

```fortran

  !-------------------------------------------------------------------
  ! List of convergence criteria available for iterative solvers 
  !-------------------------------------------------------------------
  integer(ip), parameter :: res_nrmgiven_rhs_nrmgiven  = 1  ! ||  r(i) ||g <= rtol*||  b    ||g + atol 
  integer(ip), parameter :: res_nrmgiven_res_nrmgiven  = 2  ! ||  r(i) ||g <= rtol*||  r(0) ||g + atol   
  integer(ip), parameter :: delta_rhs                  = 3  ! || dx(i) ||  <= rtol*||  b  || + atol
  integer(ip), parameter :: delta_delta                = 4  ! || dx(i) ||  <= rtol*||dx(1)|| + atol
  integer(ip), parameter :: res_res                    = 5  ! ||  r(i) ||  <= rtol*|| r(0)|| + atol
  integer(ip), parameter :: res_rhs                    = 6  ! ||  r(i) ||  <= rtol*||  b  || + atol
  integer(ip), parameter :: delta_rhs_and_res_res      = 7  ! delta_rhs    AND res_res
  integer(ip), parameter :: delta_rhs_and_res_rhs      = 8  ! delta_rhs    AND res_rhs
  integer(ip), parameter :: delta_delta_and_res_res    = 9  ! delta_delta  AND res_res
  integer(ip), parameter :: delta_delta_and_res_rhs    = 10 ! delta_delta  AND res_rhs 
                                                            ! ||.|| is the 2-norm, dx(i) = x(i) - x(i-1),
                                                            ! r(i) is the residual at the i-th iteration

  !-------------------------------------------------------------------
  ! String parameters with the names of the parameters for iterative linear solvers
  !-------------------------------------------------------------------
  character(len=*), parameter :: ils_type                      = 'iterative_linear_solver_type'
  character(len=*), parameter :: ils_rtol                      = 'iterative_linear_solver_rtol'
  character(len=*), parameter :: ils_atol                      = 'iterative_linear_solver_atol'
  character(len=*), parameter :: ils_stopping_criteria         = 'iterative_linear_solver_stopping_criteria'
  character(len=*), parameter :: ils_output_frequency          = 'iterative_linear_solver_output_frequency'
  character(len=*), parameter :: ils_max_num_iterations        = 'iterative_linear_solver_max_num_iterations'
  character(len=*), parameter :: ils_track_convergence_history = 'iterative_linear_solver_track_convergence_history'

  !-----------------------------------------------------------------
  ! Iterative linear solver names
  !-----------------------------------------------------------------
  character(len=*), parameter :: cg_name         = 'CG'         ! CG iterative linear solver
  character(len=*), parameter :: fgmres_name     = 'FGMRES'     ! FGMREs iterative linear solver
  character(len=*), parameter :: icg_name        = 'ICG'        ! ICG iterative linear solver
  character(len=*), parameter :: lfom_name       = 'LFOM'       ! LFOM iterative linear solver
  character(len=*), parameter :: lgmres_name     = 'LGMRES'     ! LGMRES iterative linear solver
  character(len=*), parameter :: minres_name     = 'MINRES'     ! MINRES iterative linear solver
  character(len=*), parameter :: rgmres_name     = 'RGMRES'     ! RGMRES iterative linear solver
  character(len=*), parameter :: richardson_name = 'RICHARDSON' ! RICHARDSON iterative linear solver


  !-----------------------------------------------------------------
  ! Some common parameters to FGMRES, LFOM, LGMRES and RGMRES iterative linear solvers
  !-----------------------------------------------------------------
  character(len=*), parameter :: ils_dkrymax             = 'iterative_linear_solver_dkrymax'
  character(len=*), parameter :: ils_orthonorm_strat     = 'iterative_linear_solver_orthonorm_strat'
  character(len=*), parameter :: orthonorm_strat_icgsro  = 'ICGSRO' 
  character(len=*), parameter :: orthonorm_strat_mgsro   = 'MGSRO'
  
  integer (ip), parameter :: mgsro  = 1 ! mgs : Modified Gram-Schmidt (appropriate for serial GMRES)
  integer (ip), parameter :: icgsro = 2 ! icgs: Iterative Classical Gram-Schmidt (appropriate for distributed GMRES)

  !-----------------------------------------------------------------
  ! Parameters used in RICHARDSON iterative linear solvers
  !-----------------------------------------------------------------
  character(len=*), parameter :: ils_relaxation = 'iterative_linear_solver_relaxation'
```
