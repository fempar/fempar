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
module par_preconditioner_dd_mlevel_bddc_names
  ! Serial modules
  use types_names
  use memor_names
  use hash_table_names
  use map_names
  use map_apply_names
  use sort_names
  use renumbering_names
  use serial_environment_names

#ifdef ENABLE_BLAS
  use blas77_interfaces_names
#endif

#ifdef ENABLE_LAPACK
  use lapack77_interfaces_names
#endif

  use abstract_solver_names
  use mesh_names
  use matrix_names
  use preconditioner_names
  use graph_names
  use vector_names
  use operator_dd_names
  use postpro_names
  use stdio_names
 
  ! Parallel modules
  use par_dd_base_names
  use psb_penv_mod_names
  use par_environment_names
  use dof_distribution_names
  use block_dof_distribution_create_names
  use par_matrix_names
  use par_vector_names
  use par_context_names
  use par_mesh_names
  use par_graph_names
  use par_timer_names

  ! Abstract modules
  use base_operand_names
  use base_operator_names

# include "debug.i90"
  implicit none
  private

  integer (ip), parameter :: default_unknowns            = interface_unknowns
  integer (ip), parameter :: default_internal_problems   = handled_by_bddc_module
  integer (ip), parameter :: default_kind_coarse_dofs    = corners
  integer (ip), parameter :: default_projection          = galerkin
  integer (ip), parameter :: default_ndime               = 2
  integer (ip), parameter :: default_co_sys_sol_strat    = serial_gather
  integer (ip), parameter :: default_nn_sys_sol_strat    = corners_rest_part_solve_expl_schur 
  integer (ip), parameter :: default_correction_mode     = additive_symmetric
  integer (ip), parameter :: default_pad_collectives     = pad
  integer (ip), parameter :: default_schur_edge_lag_mult = reuse_from_phis
  integer (ip), parameter :: default_subd_elmat_calc     = phit_minus_c_i_t_lambda  

  interface solve
     module procedure matrix_preconditioner_vector_solve, &
                      matrix_preconditioner_r1_solve, & 
                      matrix_preconditioner_r2_solve
  end interface solve


  type par_preconditioner_dd_mlevel_bddc_params_t
     ! Preconditioner params
     integer(ip) :: unknowns           =  default_unknowns            ! Ax=b (all_unknowns) .or. Sg=y (interface_unknowns) ? 
     integer(ip) :: internal_problems  =  default_internal_problems
     integer(ip) :: kind_coarse_dofs   =  default_kind_coarse_dofs    ! Corners or corners_and_edges ?
     integer(ip) :: ndime              =  default_ndime               ! Number of dimensions
     integer(ip) :: projection         =  default_projection          ! Galerkin / Petrov-Galerkin projection 
     integer(ip) :: co_sys_sol_strat   =  default_co_sys_sol_strat    ! Serial gather, serial all_gather, distributed ? 
     integer(ip) :: nn_sys_sol_strat   =  default_nn_sys_sol_strat    ! Strategy to solve the neumann problem (i.e., fine-grid correction) ?
     integer(ip) :: correction_mode    =  default_correction_mode     ! How to combine coarse and fine grid correction ?
     integer(ip) :: pad_collectives    =  default_pad_collectives     ! Pad collectives ? 
     integer(ip) :: schur_edge_lag_mult=  default_schur_edge_lag_mult ! Re-use Schur complement edge-lagrange multipliers and  A_rr^-1 C_r^T
     integer(ip) :: subd_elmat_calc    =  default_subd_elmat_calc     ! Compute A_c as \Phi_t A_i \Phi or as \Phi_t (-C_i^t \Lambda) ? 
     logical                  :: enable_constraint_weights = .false. 
     real(rp), allocatable    :: C_weights(:) 

     ! preconditioner_params and solver_control have their own defaults
     type(preconditioner_params_t)  :: ppars_harm
     type(solver_control_t)      :: spars_harm
     type(solver_control_t)      :: spars_neumann
     type(preconditioner_params_t)  :: ppars_dirichlet 
     type(solver_control_t)      :: spars_dirichlet 
     type(solver_control_t)      :: spars_coarse
     type(preconditioner_params_t)  :: ppars_coarse_serial ! co_sys_sol_strat = serial
     type (par_preconditioner_dd_mlevel_bddc_params_t ), pointer :: ppars_coarse_bddc => NULL()  ! co_sys_sol_strat = recursive

     real(rp), allocatable    :: weight(:)

     !integer(ip) :: num_neumann_solves   = 0
     !integer(ip) :: num_dirichlet_solves = 0
     !integer(ip) :: num_coarse_solves    = 0
     !integer(ip) :: neumann_its   (2000)      ! A much more sophisticated (dynamic) treatment is required here
     !integer(ip) :: dirichlet_its (2000)      ! This is already done in the new slover params type, defined as allocatable and
     !integer(ip) :: coarse_its    (2000)      ! allocated to max_it
     !integer(ip) :: harm_its_agg  = 0 

  end type par_preconditioner_dd_mlevel_bddc_params_t

  type, extends(base_operator_t) :: par_preconditioner_dd_mlevel_bddc_t
     real(rp), pointer :: weight(:) => NULL()

     ! Preconditioner params
     integer (ip)              :: unknowns           ! Ax=b (all_unknowns) .or. Sg=y (interface_unknowns) ? 
     integer (ip)              :: internal_problems
     integer (ip)              :: kind_coarse_dofs   ! Corners, corners_and_edges or other combinations ?
     integer (ip)              :: ndime              ! Number of dimensions
     integer (ip)              :: projection         ! Galerkin / Petrov-Galerkin projection 

     ! Serial gather, serial all_gather, distributed ? 
     integer (ip)              :: co_sys_sol_strat

     ! Strategy to solve the neumann problem (i.e., fine-grid correction) ?
     integer (ip)              :: nn_sys_sol_strat

     ! How to combine coarse and fine grid correction ?
     integer (ip)              :: correction_mode

     ! Pad collectives ? 
     integer (ip)              :: pad_collectives

     ! Re-use Schur complement edge-lagrange multipliers and  A_rr^-1 C_r^T
     ! for the computation of the neumann problem or compute them
     ! from scratch using spars_neumann ?
     integer(ip)               :: schur_edge_lag_mult

     ! How to calculate subdomain element matrices (subd_elmat_calc)?
     ! As \Phi_t A_i \Phi or \Phi_t (-C_i^t \Lambda) ? 
     integer(ip)               :: subd_elmat_calc

     ! BEGIN. Local info preconditioner (local data to each processor)
     integer (ip)              :: nl_coarse          ! Number of local coarse dofs

     ! (nl_coarse = nl_corners + nl_edges)
     integer (ip)              :: nl_corners      ! Number of local corners
     integer (ip)              :: nl_corners_dofs ! Sum of corner dofs
     integer (ip)              :: nl_edges        ! Number of local edges
     integer (ip)              :: nl_edges_dofs   ! Sum of edges dofs

     integer (ip), allocatable :: coarse_dofs (:)      ! List of local coarse dofs.

     integer (ip), allocatable :: perm  (:)       ! Permutation from A_i to matrix P^T A_i P 
     integer (ip), allocatable :: iperm (:)       ! Permutation from P^T A_i P to A_i

     ! Let A_i refer to p_matrix%f_matrix, then the BDDC
     ! algorithm requires a symmetric permutation of A_i, 
     ! P, s.t. P^T A_i P = [ A_cc A_cr ]
     !                     [ A_rc A_rr ],
     ! where corners are numbered first, and then the complement
     ! (edges+the rest of dofs of the subdomain). The arrays perm
     ! and iperm store the correspondence among A_i and P^T A_i P
     type ( graph_t )      :: A_rr_gr
     type ( graph_t )      :: A_rr_trans_gr
     type ( matrix_t )     :: A_rr, A_rr_trans
     real(rp) , allocatable  :: A_rc(:,:), A_cr(:,:), A_cc(:,:)
     real(rp) , allocatable  :: A_cr_trans(:,:), A_rc_trans(:,:), A_cc_trans(:,:)
     logical                 :: enable_constraint_weights
     real(rp) , allocatable  :: C_weights(:) 
          
     ! Explicitly stores the Schur complement corresponding to edge-lagrange multipliers
     ! i.e., S_rr = C_r k_rr^-1 C_r^T
     real    (rp), allocatable :: S_rr (:,:), S_rr_trans(:,:)
     integer (ip), allocatable :: ipiv (:), ipiv_trans(:)     ! pivots from DGETRF's in LAPACK

     ! Explicitly stores the product A_rr^-1 C_r^T
     real (rp), allocatable :: A_rr_inv_C_r_T (:,:), A_rr_trans_inv_C_r_T (:,:)

     ! Explicitly stores the Schur complement corresponding to edge-lagrange multipliers
     ! i.e., S_rr = C_r k_rr^-1 C_r^T
     real    (rp), allocatable :: S_rr_neumann (:,:), S_rr_trans_neumann (:,:)
     integer (ip), allocatable :: ipiv_neumann (:), ipiv_trans_neumann(:)   ! pivots from DGETRF's in LAPACK

     ! Explicitly stores the product A_rr^-1 C_r^T
     real (rp), allocatable :: A_rr_inv_C_r_T_neumann (:,:), A_rr_trans_inv_C_r_T_neumann (:,:)

     real(rp), allocatable :: rPhi (:,:), lPhi(:,:)                ! Harmonic extensions, right and left
                                                                   ! PHI basis regarding the A_c assembly

     real(rp), allocatable :: blk_lag_mul (:,:)  ! Stores the lagrange multipliers corresponding
                                                 ! to the block constrained linear system which  
                                                 ! is required to be solved for the computation of 
                                                 ! the harmonic extensions 
 
   ! Ac global matrix is the assembly of locals ( lPhi · A_i · rPhi ) 

     ! Recall that rPhi is the solution of the following sparse
     ! linear system with several dense RHS
     ! [ A_cc A_cr  I     0] [ rPhi ] = [0]  
     ! [ A_rc A_rr  0 C_r^t] [ rPhi ] = [0] 
     ! [    I    0  0     0] [lambda] = [I] 
     ! [    0   C_r 0     0] [lambda] = [I] 

     ! Then lPhi is the solution of the following sparse
     ! linear system with several dense RHS
     ! [ A_cc_trans A_rc_trans  I     0] [ lPhi ] = [0]  
     ! [ A_cr_trans A_rr_trans  0 C_r^t] [ lPhi ] = [0] 
     ! [      I         0       0     0] [lambda] = [I] 
     ! [      0        C_r      0     0] [lambda] = [I]  

     ! END. Local info preconditioner (local data to each processor)

     ! BEGIN. Global info preconditioner 
     ! (global data only in processor/s in charge of the coarse grid system)
     integer (ip)              :: ng_coarse            ! Number of global coarse dofs
                                                       ! (i.e., size of coarse grid system)

     ! AFM: Although at first glance we though that it would be convenient to replace ptr_coarse_dofs(:) 
     !      by p_mesh_c%f_mesh%pnods, on my experience it is very appropiate to preserve the former to
     !      keep the current codes clean. There are two differences among p_mesh_c%f_mesh%pnods
     !      and ptr_coarse_dofs:
     !       1) p_mesh_c%f_mesh%pnods is of size np_g and ptr_coarse_dofs is of size np_g+1 
     !       2) p_mesh_c%f_mesh%pnods uses F-based indexing and ptr_coarse_dofs C-based indexing
     !      If we use p_mesh_c%f_mesh%pnods as a replacement for ptr_coarse_dofs, then many of
     !      highly stable subroutines that we have in our code that use ptr_coarse_dofs will have
     !      to be adapted to the new definition of ptr_coarse_dofs. Besides ptr_coarse_dofs can
     !      be directly passed to v collectives (mpi_gatherv, mpi_scatterv, etc.), while with p_mesh_c%f_mesh%pnods
     !      we should transform among both kind of indexing back and forth for each call to mpi collectives.
     !      Therefore, my decision will be to keep both ptr_coarse_dofs and p_mesh_c%f_mesh%pnods separately.
     integer (ip), allocatable :: ptr_coarse_dofs (:)  ! Pointers to the list of coarse
                                                       ! dofs of each processor

     integer (ip), allocatable :: lst_coarse_dofs (:)  ! List of coarse dofs of 
                                                       ! each processor

     integer (ip), allocatable :: vars        (:)  ! Physical variables
                                                   ! associated to coarse_dofs

     ! The maximum number of coarse dofs among all subdomains.
     ! This is required ONLY when pad_collectives == pad
     integer (ip)              :: max_coarse_dofs

     ! Coarse grid system coefficient matrix
     type ( graph_t )         :: A_c_gr     ! co_sys_sol_strat = serial
     type ( matrix_t )        :: A_c 
     type ( mesh_t )          :: f_mesh_c

     type ( par_environment_t )   :: p_env_c      ! Parallel environment for the coarse-grid problem
     type ( dof_distribution_t )  :: dof_dist_c   ! co_sys_sol_strat = recursive (or distributed, not implemented yet)
     type ( renumbering_t )             :: erenumbering_c       ! Element renumbering_tbering required to pass from c_mesh to p_mesh_c
                                                ! (this is required for the assembly of the coarse-grid matrix)
     type ( par_mesh_t )          :: p_mesh_c
     type ( par_graph_t )         :: p_graph_c
     type ( par_matrix_t )        :: p_mat_c

     ! END. Global info preconditioner 
     ! (global data only in processor/s in charge of the coarse grid system)

     integer (ip)             :: symm
     integer (ip)             :: sign
     type ( preconditioner_t )     :: M_rr ,M_rr_trans
     type ( preconditioner_t )     :: M_c                                  ! co_sys_sol_strat = serial
   

     type ( par_preconditioner_dd_mlevel_bddc_t ), pointer :: p_M_c        ! co_sys_sol_strat = recursive

     type ( operator_dd_t ) :: A_II_inv ! Only required if unknowns == all_unknowns

     type ( par_matrix_t  )     , pointer  :: p_mat    => NULL()

     type (par_context_t)  :: g_context ! Fine to coarse comm
     type (par_context_t)  :: c_context ! Coarse process
     type (par_context_t)  :: d_context ! Available (coarse unused) process
     type (par_context_t)  :: b_context ! Intercommunicator betwen c_context and d_context (bcast and recursive call)

     ! par_timer objects associated to f_tasks_c_tasks_w_coarse_duties timing
     type(par_timer_t), pointer :: timer_assprec_ov_coarse
     type(par_timer_t), pointer :: timer_assprec_ov_fine
     type(par_timer_t), pointer :: timer_assprec_ov_coarse_header
     type(par_timer_t), pointer :: timer_assprec_ov_fine_header

     type(par_timer_t), pointer :: timer_fillprec_ov_coarse
     type(par_timer_t), pointer :: timer_fillprec_ov_fine
     type(par_timer_t), pointer :: timer_fillprec_ov_coarse_header
     type(par_timer_t), pointer :: timer_fillprec_ov_fine_header

     type(par_timer_t), pointer :: timer_applyprec_ov_coarse
     type(par_timer_t), pointer :: timer_applyprec_ov_fine
     type(par_timer_t), pointer :: timer_applyprec_ov_coarse_header
     type(par_timer_t), pointer :: timer_applyprec_ov_coarse_tail

     type(par_timer_t), pointer :: timer_coll_assprec
     type(par_timer_t), pointer :: timer_coll_fillprec
     type(par_timer_t), pointer :: timer_coll_applyprec

     type(preconditioner_params_t)        :: ppars_harm 
     type(solver_control_t), pointer   :: spars_harm
     type(solver_control_t), pointer   :: spars_neumann
     type(preconditioner_params_t)        :: ppars_dirichlet 
     type(solver_control_t), pointer   :: spars_dirichlet 
     type(solver_control_t), pointer   :: spars_coarse
     type(preconditioner_params_t)        :: ppars_coarse_serial    ! co_sys_sol_strat = serial_gather
     type (par_preconditioner_dd_mlevel_bddc_params_t ) :: ppars_coarse_bddc      ! co_sys_sol_strat = recursive_bddc

     integer(ip) :: num_neumann_solves   = 0
     integer(ip) :: num_dirichlet_solves = 0
     integer(ip) :: num_coarse_solves    = 0

     integer(ip) :: neumann_its   (2000)  ! A much more sophisticated (dynamic) treatment is required here
     integer(ip) :: dirichlet_its (2000)
     integer(ip) :: coarse_its    (2000)

     integer(ip) :: harm_its_agg         = 0 
   contains
     procedure :: apply => par_preconditioner_dd_mlevel_bddc_apply_tbp
     procedure :: apply_fun => par_preconditioner_dd_mlevel_bddc_apply_fun_tbp
     procedure :: fill_values => par_preconditioner_dd_mlevel_bddc_fill_values_tbp
     procedure :: free => par_preconditioner_dd_mlevel_bddc_free_tbp
  end type par_preconditioner_dd_mlevel_bddc_t

  ! public :: default_kind_coarse_dofs, default_co_sys_sol_strat, default_ndime

  ! Types
  public :: par_preconditioner_dd_mlevel_bddc_t, par_preconditioner_dd_mlevel_bddc_params_t

  ! Functions
  public :: par_preconditioner_dd_mlevel_bddc_create, par_preconditioner_dd_mlevel_bddc_ass_struct, &
            par_preconditioner_dd_mlevel_bddc_fill_val, par_preconditioner_dd_mlevel_bddc_free, &
            par_preconditioner_dd_mlevel_bddc_apply_all_unk, &
            par_preconditioner_dd_mlevel_bddc_static_condensation, &
            par_preconditioner_dd_mlevel_bddc_report, &
            par_preconditioner_dd_mlevel_bddc_time_report, &
            par_preconditioner_dd_mlevel_bddc_fill_val_phase_1, &
            par_preconditioner_dd_mlevel_bddc_fill_val_phase_2, &
            assemble_A_c

  logical, parameter :: debug_verbose_level_1 = .false. 
  logical, parameter :: debug_verbose_level_2 = .false.
  logical, parameter :: debug_verbose_level_3 = .false.  ! Prints harmonic extensions 
                                                              ! to gid files, and coarse grid system 
                                                              ! coefficient matrix to matrix market 

contains

  recursive subroutine par_preconditioner_dd_mlevel_bddc_create( p_mat, mlbddc, mlbddc_params )
#ifdef MPI_MOD
use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    !implicit none
    ! Parameters
    type(par_matrix_t)                       , target ,intent(in)  :: p_mat
    type(par_preconditioner_dd_mlevel_bddc_t)                ,intent(out) :: mlbddc
    type(par_preconditioner_dd_mlevel_bddc_params_t), target ,intent(in)  :: mlbddc_params

    ! Locals
    integer            :: i, istat
    logical            :: i_am_coarse_task, i_am_fine_task, i_am_higher_level_task
    type (vector_t)  :: dum 
    type (matrix_t)  :: mat_dum
    type (par_context_t) :: dum_context


    ! This module requires w_context (all tasks involved here), p_context (fine tasks, which
    ! are in charge of integration in the finer level) and q_context(unused tasks) to be 
    ! created by a call to par_context_create_by_split.
    ! b_context should also be created in advance, as it might be required by other preconditioners
    ! that do not have a multilevel structure (say, e.g., the diagonal preconditioner)
    assert ( associated(p_mat%p_env) )
    assert ( p_mat%p_env%created )
    assert ( p_mat%p_env%num_levels > 1 )
    assert ( associated(p_mat%dof_dist) )


    assert ( p_mat%p_env%w_context%created .eqv. .true. )
    assert ( associated(p_mat%p_env%p_context) )
    assert ( p_mat%p_env%p_context%created .eqv. .true. )
    assert ( associated(p_mat%p_env%q_context) )
    assert ( p_mat%p_env%q_context%created .eqv. .true. )
    ! Check appropriate assignment to context: 1) I'm in w_context and 2) I'm in p_context OR in q_context but not in both
    assert ( p_mat%p_env%w_context%iam >= 0)
    assert ( (p_mat%p_env%p_context%iam >=0 .and. p_mat%p_env%q_context%iam < 0) .or. (p_mat%p_env%p_context%iam < 0 .and. p_mat%p_env%q_context%iam >= 0))
    assert ( associated(p_mat%p_env%b_context) )
    assert ( p_mat%p_env%b_context%created .eqv. .true. )

    mlbddc%symm    = p_mat%f_matrix%symm
    mlbddc%sign    = p_mat%f_matrix%sign

    !assert ( mlbddc_params%unknowns == interface_unknowns .or. mlbddc_params%unknowns == all_unknowns )
    assert ( mlbddc_params%unknowns == all_unknowns ) ! Only global unknowns accepted in the multilevel?
    mlbddc%unknowns  = mlbddc_params%unknowns

    assert ( mlbddc_params%internal_problems == handled_by_bddc_module )  ! .or. mlbddc_params%internal_problems == handled_by_user_driver  )
    mlbddc%internal_problems  = mlbddc_params%internal_problems

    assert ( mlbddc_params%ndime == 2 .or. mlbddc_params%ndime == 3 )
    mlbddc%ndime = mlbddc_params%ndime    

    if (mlbddc%ndime==2) then 
       ! Only three kinds of coarse dofs: 1) only corners; 2) both corners and edges; 3) only edges 
       assert(mlbddc_params%kind_coarse_dofs == corners.or.mlbddc_params%kind_coarse_dofs==corners_and_edges.or.mlbddc_params%kind_coarse_dofs==edges)
    end if

    if (mlbddc%ndime==3) then 
       ! Seven kinds of coarse dofs
       assert(mlbddc_params%kind_coarse_dofs==corners.or.mlbddc_params%kind_coarse_dofs==corners_and_edges.or.mlbddc_params%kind_coarse_dofs==corners_edges_and_faces.or.mlbddc_params%kind_coarse_dofs== edges.or.mlbddc_params%kind_coarse_dofs==faces.or.mlbddc_params%kind_coarse_dofs == edges_and_faces.or.mlbddc_params%kind_coarse_dofs == corners_and_faces)
    end if
    mlbddc%kind_coarse_dofs  = mlbddc_params%kind_coarse_dofs

    assert( mlbddc_params%projection == galerkin .or. mlbddc_params%projection == petrov_galerkin )
    mlbddc%projection = mlbddc_params%projection 

    assert( mlbddc_params%enable_constraint_weights .or. .not. mlbddc_params%enable_constraint_weights  )
    mlbddc%enable_constraint_weights = mlbddc_params%enable_constraint_weights    

    ! Up to the moment, two strategies for the solution of the coarse grid system are implemented
    assert ( mlbddc_params%co_sys_sol_strat == serial_gather .or. mlbddc_params%co_sys_sol_strat == recursive_bddc)
    mlbddc%co_sys_sol_strat  = mlbddc_params%co_sys_sol_strat
    
    assert (mlbddc_params%nn_sys_sol_strat==corners_rest_part_solve_expl_schur.or.mlbddc_params%nn_sys_sol_strat==direct_solve_constrained_problem )
    mlbddc%nn_sys_sol_strat = mlbddc_params%nn_sys_sol_strat

    mlbddc%correction_mode = mlbddc_params%correction_mode

    assert ( mlbddc_params%pad_collectives == pad .or. mlbddc_params%pad_collectives == nopad )
    mlbddc%pad_collectives = mlbddc_params%pad_collectives

    assert ( mlbddc_params%schur_edge_lag_mult == reuse_from_phis .or. mlbddc_params%schur_edge_lag_mult == compute_from_scratch )
    mlbddc%schur_edge_lag_mult = mlbddc_params%schur_edge_lag_mult

    assert ( mlbddc_params%subd_elmat_calc == phit_a_i_phi .or. mlbddc_params%subd_elmat_calc == phit_minus_c_i_t_lambda )
    mlbddc%subd_elmat_calc = mlbddc_params%subd_elmat_calc

    mlbddc%p_mat    => p_mat

    mlbddc%num_dirichlet_solves = 0
    mlbddc%num_neumann_solves   = 0
    mlbddc%num_coarse_solves    = 0

   if(allocated(mlbddc_params%weight)) then
       if (mlbddc%p_mat%p_env%p_context%iam>=0) then
          mlbddc%weight => mlbddc_params%weight(1:p_mat%dof_dist%nb)
       else
          nullify(mlbddc%weight)
       end if
    else
       nullify(mlbddc%weight)
    end if

    ! Now create c_context (the coarse context) and d_context (the set of unused tasks) by spliting q_context
    ! allocate(mlbddc%p_mat%dof_dist%c_context, stat = istat)
    ! assert(istat==0)
    ! allocate(mlbddc%p_mat%dof_dist%d_context, stat = istat)
    ! assert(istat==0)
    ! if(mlbddc%p_mat%dof_dist%q_context%iam >= 0) then
    !    if(mlbddc%p_mat%dof_dist%q_context%iam < mlbddc%p_mat%dof_dist%num_parts(2)) then
    !       call par_context_create ( inhouse, 1, mlbddc%p_mat%dof_dist%c_context, mlbddc%p_mat%dof_dist%d_context, mlbddc%p_mat%dof_dist%w_context )
    !    else
    !       call par_context_create ( inhouse, 2, mlbddc%p_mat%dof_dist%d_context, mlbddc%p_mat%dof_dist%c_context, mlbddc%p_mat%dof_dist%w_context )
    !    end if
    ! else
    !    call par_context_create ( inhouse, -1, mlbddc%p_mat%dof_dist%c_context, mlbddc%p_mat%dof_dist%d_context, mlbddc%p_mat%dof_dist%w_context )
    ! end if
    ! write(*,*) 'c_context and d_context created'
    if(mlbddc%p_mat%p_env%q_context%iam >= 0) then
       if(mlbddc%p_mat%p_env%q_context%iam < mlbddc%p_mat%p_env%num_parts(2)) then
          call par_context_create ( 1, mlbddc%c_context, mlbddc%d_context, mlbddc%p_mat%p_env%q_context )
       else
          call par_context_create ( 2, mlbddc%d_context, mlbddc%c_context, mlbddc%p_mat%p_env%q_context )
       end if
    else
       call par_context_null (mlbddc%c_context)
       call par_context_null (mlbddc%d_context)
    end if

    ! Now create g_context, where communication actually occurs
    if(mlbddc%p_mat%p_env%p_context%iam >= 0) then
       ! my_color is id_parts(this_level+1)
       call par_context_create ( mlbddc%p_mat%p_env%id_parts(2)    , mlbddc%g_context, dum_context, mlbddc%p_mat%p_env%w_context )
    else if(mlbddc%c_context%iam >= 0) then
       call par_context_create ( mlbddc%c_context%iam+1      , mlbddc%g_context, dum_context, mlbddc%p_mat%p_env%w_context )
    else
       call par_context_create (                           -1, mlbddc%g_context, dum_context, mlbddc%p_mat%p_env%w_context )
    end if

    ! Which duties I have?
    i_am_fine_task = (mlbddc%p_mat%p_env%p_context%iam >= 0)
    i_am_coarse_task = (mlbddc%c_context%iam >= 0)
    i_am_higher_level_task = (mlbddc%d_context%iam >= 0)

    if ( debug_verbose_level_1 ) then
      write(*,*) 'LEVELS: ', mlbddc%p_mat%p_env%num_levels, 'w_context: ',mlbddc%p_mat%p_env%w_context%iam, mlbddc%p_mat%p_env%w_context%np
      write(*,*) 'LEVELS: ', mlbddc%p_mat%p_env%num_levels, 'p_context: ',mlbddc%p_mat%p_env%p_context%iam, mlbddc%p_mat%p_env%p_context%np
      write(*,*) 'LEVELS: ', mlbddc%p_mat%p_env%num_levels, 'q_context: ',mlbddc%p_mat%p_env%q_context%iam, mlbddc%p_mat%p_env%q_context%np
      write(*,*) 'LEVELS: ', mlbddc%p_mat%p_env%num_levels, 'c_context: ',mlbddc%c_context%iam, mlbddc%c_context%np
      write(*,*) 'LEVELS: ', mlbddc%p_mat%p_env%num_levels, 'd_context: ',mlbddc%d_context%iam, mlbddc%d_context%np
      write(*,*) 'LEVELS: ', mlbddc%p_mat%p_env%num_levels, 'g_context: ',mlbddc%g_context%iam, mlbddc%g_context%np
      write(*,*) 'LEVELS: ', mlbddc%p_mat%p_env%num_levels, 'b_context: ',mlbddc%p_mat%p_env%b_context%iam, mlbddc%p_mat%p_env%b_context%np
    end if


    ! Both the global number of objects and the maximum number
    ! of parts surrounding any node on the interface are required
    ! during preconditioner set-up in all tasks
    ! call par_partition_bcast(mlbddc%p_part,mlbddc%p_mat%p_env%f_part%omap%ng)
    ! call par_partition_bcast(mlbddc%p_part,mlbddc%p_mat%p_env%f_part%max_nparts)
    ! call par_partition_bcast(mlbddc%p_part,mlbddc%p_mat%p_env%num_levels)
    ! Are these arrays allocated in the coarse process?
    !call par_partition_bcast(mlbddc%p_part,mlbddc%p_mat%p_env%num_parts)
    !call par_partition_bcast(mlbddc%p_part,mlbddc%p_mat%p_env%id_parts)

    mlbddc%ppars_harm      = mlbddc_params%ppars_harm

    allocate(mlbddc%spars_harm, stat=istat)
    check( istat == 0 ) 
    mlbddc%spars_harm      = mlbddc_params%spars_harm

    allocate(mlbddc%spars_neumann, stat=istat)
    check( istat == 0 ) 
    mlbddc%spars_neumann   = mlbddc_params%spars_neumann

    mlbddc%ppars_dirichlet = mlbddc_params%ppars_dirichlet

    allocate(mlbddc%spars_dirichlet, stat=istat)
    check( istat == 0 ) 
    mlbddc%spars_dirichlet = mlbddc_params%spars_dirichlet

    allocate(mlbddc%spars_coarse, stat=istat)
    check( istat == 0 ) 
    mlbddc%spars_coarse    = mlbddc_params%spars_coarse

    if(mlbddc%co_sys_sol_strat == serial_gather) then
       mlbddc%ppars_coarse_serial = mlbddc_params%ppars_coarse_serial
    else if(mlbddc%co_sys_sol_strat == recursive_bddc) then
       mlbddc%ppars_coarse_bddc = mlbddc_params%ppars_coarse_bddc
    end if

   if ( temp_fine_coarse_grid_overlap ) then 
       if ( i_am_coarse_task ) then
          allocate(mlbddc%timer_assprec_ov_coarse, stat=istat)
          check(istat == 0)
          allocate(mlbddc%timer_assprec_ov_coarse_header, stat=istat)
          check(istat == 0)
          allocate(mlbddc%timer_fillprec_ov_coarse, stat=istat)
          check(istat == 0)
          allocate(mlbddc%timer_applyprec_ov_coarse, stat=istat)
          check(istat == 0)
          call par_timer_create ( mlbddc%timer_assprec_ov_coarse,  'Symb. Prec. region C',  mlbddc%c_context%icontxt )
          call par_timer_create ( mlbddc%timer_assprec_ov_coarse_header,  'Symb. Prec. region C header',  mlbddc%c_context%icontxt )
          call par_timer_create ( mlbddc%timer_fillprec_ov_coarse, 'Fill. Prec. region C',  mlbddc%c_context%icontxt )
          call par_timer_create ( mlbddc%timer_applyprec_ov_coarse,  'Apply Prec. region C',  mlbddc%c_context%icontxt )
       else if ( i_am_fine_task ) then
          allocate ( mlbddc%timer_assprec_ov_fine, stat=istat)
          check(istat == 0)
          allocate ( mlbddc%timer_assprec_ov_fine_header,stat=istat)
          check(istat == 0)
          allocate ( mlbddc%timer_fillprec_ov_fine, stat=istat)
          check(istat == 0)
          allocate ( mlbddc%timer_fillprec_ov_fine_header, stat=istat)
          check(istat == 0)
          allocate ( mlbddc%timer_fillprec_ov_coarse_header,stat=istat)
          check(istat == 0)
          allocate ( mlbddc%timer_applyprec_ov_fine, stat=istat)
          check(istat == 0)
          allocate ( mlbddc%timer_applyprec_ov_coarse_header,stat=istat)
          check(istat == 0)
          allocate ( mlbddc%timer_applyprec_ov_coarse_tail,stat=istat)
          check(istat == 0)
          call par_timer_create ( mlbddc%timer_assprec_ov_fine, 'Symb. Prec. region F',  mlbddc%p_mat%p_env%p_context%icontxt )
          call par_timer_create ( mlbddc%timer_assprec_ov_fine_header, 'Symb. Prec. region F header',  mlbddc%p_mat%p_env%p_context%icontxt )
          call par_timer_create ( mlbddc%timer_fillprec_ov_fine, 'Fill. Prec. region F',  mlbddc%p_mat%p_env%p_context%icontxt )
          call par_timer_create ( mlbddc%timer_fillprec_ov_fine_header, 'Fill. Prec. region F header',  mlbddc%p_mat%p_env%p_context%icontxt )
          call par_timer_create ( mlbddc%timer_fillprec_ov_coarse_header, 'Fill. Prec. region C header',  mlbddc%p_mat%p_env%p_context%icontxt )
          call par_timer_create ( mlbddc%timer_applyprec_ov_fine, 'Apply Prec. region F',  mlbddc%p_mat%p_env%p_context%icontxt )
          call par_timer_create ( mlbddc%timer_applyprec_ov_coarse_header,  'Apply Prec. region C header',  mlbddc%p_mat%p_env%p_context%icontxt )
          call par_timer_create ( mlbddc%timer_applyprec_ov_coarse_tail, 'Apply Prec. region C tail',  mlbddc%p_mat%p_env%p_context%icontxt )
       end if
       if ( i_am_coarse_task .or. i_am_fine_task ) then
         allocate ( mlbddc%timer_coll_assprec  , stat=istat)
         check(istat == 0)
         allocate ( mlbddc%timer_coll_fillprec , stat=istat)
         check(istat == 0)
         allocate ( mlbddc%timer_coll_applyprec, stat=istat)
         check(istat == 0)
         call par_timer_create ( mlbddc%timer_coll_assprec  , 'Collective Assembly preconditioner',  mlbddc%g_context%icontxt )
         call par_timer_create ( mlbddc%timer_coll_fillprec , 'Collective Fill preconditioner'    ,  mlbddc%g_context%icontxt )
         call par_timer_create ( mlbddc%timer_coll_applyprec, 'Collective Apply preconditioner'   ,  mlbddc%g_context%icontxt )
       end if          
    end if

    if ( i_am_fine_task ) then     
 
       if (mlbddc%enable_constraint_weights) then 
          if(  allocated(mlbddc_params%C_weights )) then 
             call memalloc ( p_mat%dof_dist%nb, mlbddc%C_weights, __FILE__, __LINE__ )
             mlbddc%C_weights = mlbddc_params%C_weights 
          end if
       end if

       ! BEG. FINE-GRID PROBLEM DUTIES
       if ( mlbddc%nn_sys_sol_strat == corners_rest_part_solve_expl_schur ) then

          call matrix_create(p_mat%f_matrix%type, p_mat%f_matrix%symm, mlbddc%A_rr, p_mat%f_matrix%sign)

          if ( mlbddc%projection == petrov_galerkin ) then 
             call matrix_create(p_mat%f_matrix%type, p_mat%f_matrix%symm, mlbddc%A_rr_trans, p_mat%f_matrix%sign)
          end if

          if ( mlbddc%internal_problems == handled_by_bddc_module) then
             call preconditioner_create( mlbddc%A_rr, mlbddc%M_rr, mlbddc%ppars_harm)

             if ( mlbddc%projection == petrov_galerkin ) then 
                call preconditioner_create( mlbddc%A_rr, mlbddc%M_rr_trans, mlbddc%ppars_harm)
             end if
          else
             check(.false.)
          end if

       else if (  mlbddc%nn_sys_sol_strat == direct_solve_constrained_problem ) then

          call matrix_create( p_mat%f_matrix%type, p_mat%f_matrix%symm, mlbddc%A_rr, indefinite)
          if ( mlbddc%projection == petrov_galerkin ) then 
             call matrix_create( p_mat%f_matrix%type, p_mat%f_matrix%symm, mlbddc%A_rr_trans, p_mat%f_matrix%sign)
          end if
          if ( mlbddc%internal_problems == handled_by_bddc_module) then
             call preconditioner_create( mlbddc%A_rr, mlbddc%M_rr, mlbddc%ppars_harm)
             if ( mlbddc%projection==petrov_galerkin ) then 
                call preconditioner_create( mlbddc%A_rr, mlbddc%M_rr_trans, mlbddc%ppars_harm)
             end if
          else
             check(.false.)
          end if
       end if


       if ( mlbddc%unknowns == all_unknowns ) then
          if ( mlbddc%internal_problems == handled_by_bddc_module) then
              call operator_dd_create ( p_mat%f_matrix, & 
                                            p_mat%dof_dist, &
                                            mlbddc%A_II_inv, & 
                                            spars=mlbddc%spars_dirichlet, &
                                            ppars=mlbddc%ppars_dirichlet)
                                           ! symm=symm_true, &
                                           ! sign=positive_definite )
          else
             check (.false.)
!!$             call operator_dd_create ( p_mat%f_matrix, & 
!!$                                           p_mat%dof_dist%f_part, &
!!$                                           mlbddc%A_II_inv, & 
!!$                                           dirichlet_internally=.false., &
!!$                                           spars=mlbddc%spars_dirichlet, &
!!$                                           ppars=mlbddc%ppars_dirichlet)
!!$                                           ! symm=symm_true, &
!!$                                           ! sign=positive_definite )
          end if
       end if
       ! END FINE-GRID PROBLEM DUTIES    
    end if

    if ( i_am_coarse_task .or. i_am_higher_level_task ) then
       if ( mlbddc%internal_problems == handled_by_bddc_module) then
          ! BEG. COARSE-GRID PROBLEM DUTIES

          if(mlbddc%co_sys_sol_strat == serial_gather) then ! There are only coarse tasks
             call matrix_create( p_mat%f_matrix%type, p_mat%f_matrix%symm, mlbddc%A_c)
             call preconditioner_create( mlbddc%A_c, mlbddc%M_c, mlbddc%ppars_coarse_serial)
          
          !else if(mlbddc%co_sys_sol_strat == distributed) then

          else if(mlbddc%co_sys_sol_strat == recursive_bddc) then
             assert(mlbddc%p_mat%p_env%num_levels>2) ! Assuming last level serial

             ! Now create b_context (intercomm)
             assert(istat==0)
             call par_context_create ( mlbddc%p_mat%p_env%q_context, mlbddc%c_context, mlbddc%d_context, mlbddc%b_context )

             ! Create coarse partition (assign c_context -> p_context, q_context -> w_context, d_context -> q_context) 
             call par_environment_create ( mlbddc%p_env_c,         &
                                           mlbddc%p_mat%p_env%q_context, &  ! Intracomm including c & d
                                           mlbddc%c_context, &  ! Fine-tasks intracomm
                                           mlbddc%d_context, &  ! Rest-of-tasks intracomm
                                           mlbddc%b_context,       &  ! Intercomm c<=>d
                                           mlbddc%p_mat%p_env%num_levels-1, &  
                                           mlbddc%p_mat%p_env%id_parts(2:mlbddc%p_mat%p_env%num_levels), &
                                           mlbddc%p_mat%p_env%num_parts(2:mlbddc%p_mat%p_env%num_levels) )

             ! Create coarse mesh
             call par_mesh_create ( mlbddc%p_env_c, mlbddc%p_mesh_c )

             ! Create coarse graph 
             call par_graph_create ( mlbddc%dof_dist_c, mlbddc%p_env_c, mlbddc%p_graph_c )

             ! Create coarse matrix
             call par_matrix_create( p_mat%f_matrix%type, & 
                                     p_mat%f_matrix%symm, & 
                                     mlbddc%dof_dist_c, &
                                     mlbddc%dof_dist_c, &
                                     mlbddc%p_env_c, &
                                     mlbddc%p_mat_c)

             ! Allocate inverse
             allocate(mlbddc%p_M_c, stat=istat)
             check(istat==0)

             ! Recursively call bddc_create
             call par_preconditioner_dd_mlevel_bddc_create(mlbddc%p_mat_c, mlbddc%p_M_c, mlbddc%ppars_coarse_bddc )

          end if
       else
          check(.false.)
       end if

       ! END COARSE-GRID PROBLEM DUTIES    
    end if

  end subroutine par_preconditioner_dd_mlevel_bddc_create

  !=================================================================================================

  recursive subroutine par_preconditioner_dd_mlevel_bddc_free ( mlbddc, mode )
    implicit none
    ! Parameters
    type(par_preconditioner_dd_mlevel_bddc_t), intent(inout) :: mlbddc
    integer(ip)                     , intent(in)    :: mode

    ! Locals
    type (vector_t) :: dum 
    integer           :: iam, istat
    logical           :: i_am_fine_task, i_am_coarse_task, i_am_higher_level_task

    ! The routine requires the partition/context info
    assert ( associated(mlbddc%p_mat%p_env) )
    assert ( mlbddc%p_mat%p_env%created )
    assert ( mlbddc%p_mat%p_env%num_levels > 1 ) 

    assert ( associated(mlbddc%p_mat%p_env%w_context) )
    assert ( mlbddc%p_mat%p_env%w_context%created .eqv. .true. )
    assert ( associated(mlbddc%p_mat%p_env%p_context) )
    assert ( mlbddc%p_mat%p_env%p_context%created .eqv. .true. )
    assert ( associated(mlbddc%p_mat%p_env%q_context) )
    assert ( mlbddc%p_mat%p_env%q_context%created .eqv. .true. )
    ! Check appropriate assignment to context: 1) I'm in w_context and 2) I'm in p_context OR in q_context but not in both
    assert ( mlbddc%p_mat%p_env%w_context%iam >= 0)
    assert ( (mlbddc%p_mat%p_env%p_context%iam >=0 .and. mlbddc%p_mat%p_env%q_context%iam < 0) .or. (mlbddc%p_mat%p_env%p_context%iam < 0 .and. mlbddc%p_mat%p_env%q_context%iam >= 0))
    assert ( mlbddc%c_context%created .eqv. .true. )
    assert ( mlbddc%d_context%created .eqv. .true. )
    assert ( mlbddc%g_context%created .eqv. .true. )
    assert ( associated(mlbddc%p_mat%p_env%b_context) )
    assert ( mlbddc%p_mat%p_env%b_context%created .eqv. .true. )

    ! Which duties I have?
    i_am_fine_task = (mlbddc%p_mat%p_env%p_context%iam >= 0)
    i_am_coarse_task = (mlbddc%c_context%iam >= 0)
    i_am_higher_level_task = (mlbddc%d_context%iam >= 0)
    if(i_am_fine_task) iam = mlbddc%p_mat%p_env%p_context%iam
    if(i_am_coarse_task) iam = mlbddc%c_context%iam
    if(i_am_higher_level_task) iam = mlbddc%d_context%iam

    if ( mode == free_clean ) then
       if ( i_am_fine_task ) then
          ! BEG. FINE-GRID PROBLEM DUTIES
          if ( mlbddc%internal_problems == handled_by_bddc_module) then
             call preconditioner_free ( preconditioner_free_clean  , mlbddc%M_rr )
             if (mlbddc%projection == petrov_galerkin )  then 
                call preconditioner_free ( preconditioner_free_clean, mlbddc%M_rr_trans ) 
             end if
          else
             check(.false.)
          end if
          if ( mlbddc%unknowns == all_unknowns ) then
             call operator_dd_free ( mlbddc%A_II_inv, free_clean )
          end if
          
          if (allocated (mlbddc%C_weights)) call memfree ( mlbddc%C_weights, __FILE__, __LINE__ )
          ! END FINE-GRID PROBLEM DUTIES
       end if

       ! BEG. COARSE-GRID PROBLEM DUTIES
       if ( mlbddc%internal_problems == handled_by_bddc_module) then
          if(mlbddc%co_sys_sol_strat == serial_gather) then

             if ( i_am_coarse_task ) then
                call preconditioner_free ( preconditioner_free_clean  , mlbddc%M_c )
             end if
             
          else if(mlbddc%co_sys_sol_strat == recursive_bddc) then

             assert(mlbddc%p_mat%p_env%num_levels>2) ! Assuming last level direct
             
             if ( i_am_coarse_task .or. i_am_higher_level_task ) then
                ! Recursively call bddc_free
                call par_preconditioner_dd_mlevel_bddc_free(mlbddc%p_M_c, free_clean )
                
                ! These lines should be uncommented when the structures are actually filled
                ! In fact we should add flags to indicate that in each free routine
                
                ! Free coarse matrix
                call par_matrix_free( mlbddc%p_mat_c, free_clean)
                
                ! Free coarse graph
                call par_graph_free ( mlbddc%p_graph_c, free_clean )
                
                ! Free coarse mesh
                call par_mesh_free ( mlbddc%p_mesh_c, free_clean )
                
                ! Free coarse environment
                call par_environment_free ( mlbddc%p_env_c )
                
                call par_context_free (mlbddc%c_context,.false.)
                call par_context_free (mlbddc%d_context,.false.)
                call par_context_free (mlbddc%g_context,.false.)
                call par_context_free (mlbddc%b_context,.false.)
                
                ! Deallocate inverse
                deallocate(mlbddc%p_M_c, stat=istat)
                assert(istat==0)
                
             end if
             
             ! deallocate(mlbddc%p_mat%p_env%c_context, stat = istat)
             ! assert(istat==0)
             ! deallocate(mlbddc%d_context, stat = istat)
             ! assert(istat==0)
             
          end if
       else
          check(.false.)
       end if
       ! END COARSE-GRID PROBLEM DUTIES

       if ( temp_fine_coarse_grid_overlap ) then 
         if ( i_am_coarse_task ) then
           deallocate ( mlbddc%timer_assprec_ov_coarse)
           deallocate ( mlbddc%timer_assprec_ov_coarse_header)
           deallocate ( mlbddc%timer_fillprec_ov_coarse)
           deallocate ( mlbddc%timer_applyprec_ov_coarse)
         else if ( i_am_fine_task ) then
           deallocate ( mlbddc%timer_assprec_ov_fine)
           deallocate ( mlbddc%timer_assprec_ov_fine_header)
           deallocate ( mlbddc%timer_fillprec_ov_fine)
           deallocate ( mlbddc%timer_fillprec_ov_fine_header)
           deallocate ( mlbddc%timer_fillprec_ov_coarse_header)
           deallocate ( mlbddc%timer_applyprec_ov_fine)
           deallocate ( mlbddc%timer_applyprec_ov_coarse_header)
           deallocate ( mlbddc%timer_applyprec_ov_coarse_tail)
         end if
         if ( i_am_coarse_task .or. i_am_fine_task ) then
           deallocate ( mlbddc%timer_coll_assprec)
           deallocate ( mlbddc%timer_coll_fillprec)
           deallocate ( mlbddc%timer_coll_applyprec)
         end if          
       end if
       
       nullify ( mlbddc%p_mat  )
       deallocate(mlbddc%spars_harm)
       deallocate(mlbddc%spars_neumann)
       deallocate(mlbddc%spars_dirichlet)
       deallocate(mlbddc%spars_coarse)
       return  
    end if

    ! TODO: code below needs to be modified after structures are filled...
    if ( mode == free_only_struct  ) then
       if ( i_am_fine_task ) then
          ! BEG. FINE-GRID PROBLEM DUTIES
          if ( mlbddc%internal_problems == handled_by_bddc_module) then
             call preconditioner_free ( preconditioner_free_struct, mlbddc%M_rr )
            if ( mlbddc%projection == petrov_galerkin ) then 
             call preconditioner_free ( preconditioner_free_struct, mlbddc%M_rr_trans )
            end if
          end if

          if ( mlbddc%unknowns == all_unknowns ) then
             call operator_dd_free ( mlbddc%A_II_inv, free_only_struct )
          end if

          call matrix_free ( mlbddc%A_rr, free_only_struct )
          call graph_free ( mlbddc%A_rr_gr )
 
          if ( mlbddc%projection == petrov_galerkin ) then 
          call matrix_free ( mlbddc%A_rr_trans, free_only_struct )
            end if
   
          call memfree ( mlbddc%coarse_dofs,__FILE__,__LINE__)

          if (  mlbddc%nn_sys_sol_strat == corners_rest_part_solve_expl_schur ) then
             call memfree ( mlbddc%perm,__FILE__,__LINE__)
             call memfree ( mlbddc%iperm,__FILE__,__LINE__)
          end if

          ! END FINE-GRID PROBLEM DUTIES
       end if
       
       ! BEG. COARSE-GRID PROBLEM DUTIES
       if(mlbddc%co_sys_sol_strat == serial_gather) then

          if ( i_am_coarse_task ) then
             if ( mlbddc%internal_problems == handled_by_bddc_module) then
                call preconditioner_free ( preconditioner_free_struct, mlbddc%M_c )
             end if
             call matrix_free ( mlbddc%A_c, free_only_struct )
             call graph_free ( mlbddc%A_c_gr )
             call mesh_free ( mlbddc%f_mesh_c )
             call memfree ( mlbddc%vars, __FILE__, __LINE__)
             call memfree ( mlbddc%ptr_coarse_dofs, __FILE__, __LINE__)
          end if

       else if(mlbddc%co_sys_sol_strat == recursive_bddc) then

          assert(mlbddc%p_mat%p_env%num_levels>2) 
          
          if ( i_am_coarse_task .or. i_am_higher_level_task ) then
             ! Recursively call bddc_free
             call par_preconditioner_dd_mlevel_bddc_free(mlbddc%p_M_c, free_only_struct )
             
             if ( i_am_coarse_task ) then
                call renumbering_free ( mlbddc%erenumbering_c )
                call memfree ( mlbddc%vars, __FILE__, __LINE__)
                call memfree ( mlbddc%ptr_coarse_dofs, __FILE__, __LINE__)
             end if
             
             ! Free coarse partition
             call dof_distribution_free ( mlbddc%dof_dist_c )
             
             ! Free coarse mesh
             call par_mesh_free ( mlbddc%p_mesh_c, free_only_struct )
             
             ! Free coarse graph
             call par_graph_free ( mlbddc%p_graph_c, free_only_struct )
             
             ! Free coarse matrix
             call par_matrix_free( mlbddc%p_mat_c, free_only_struct)
          end if
          
       end if
       ! END COARSE-GRID PROBLEM DUTIES

    else if ( mode == free_only_values ) then

       if ( i_am_fine_task ) then
          ! BEG. FINE-GRID PROBLEM DUTIES
          if ( mlbddc%internal_problems == handled_by_bddc_module) then
             call preconditioner_free ( preconditioner_free_values  , mlbddc%M_rr )
             if (mlbddc%projection == petrov_galerkin )  then 
             call preconditioner_free ( preconditioner_free_values  , mlbddc%M_rr_trans )
             end if
          end if
          if ( mlbddc%unknowns == all_unknowns ) then
             call operator_dd_free ( mlbddc%A_II_inv, free_only_values )
          end if

          call matrix_free ( mlbddc%A_rr, free_only_values )

          call memfree ( mlbddc%rPhi,__FILE__,__LINE__)
          call memfree ( mlbddc%lPhi,__FILE__,__LINE__)
          call memfree ( mlbddc%blk_lag_mul,__FILE__,__LINE__) 

             if (mlbddc%projection == petrov_galerkin )  then 
          call matrix_free ( mlbddc%A_rr_trans, free_only_values )
             end if

          if ( mlbddc%nn_sys_sol_strat == corners_rest_part_solve_expl_schur ) then
             if ( mlbddc%kind_coarse_dofs == corners_and_edges .or. & 
                  mlbddc%kind_coarse_dofs == corners_edges_and_faces .or. &
                  mlbddc%kind_coarse_dofs == edges .or. & 
                  mlbddc%kind_coarse_dofs == edges_and_faces .or. & 
                  mlbddc%kind_coarse_dofs == faces .or. & 
                  mlbddc%kind_coarse_dofs == corners_and_faces ) then
                call memfree ( mlbddc%S_rr,__FILE__,__LINE__)
                
                call memfree ( mlbddc%A_rr_inv_C_r_T,__FILE__,__LINE__)

                if (mlbddc%projection == petrov_galerkin )  then 
                   call memfree ( mlbddc%S_rr_trans,__FILE__,__LINE__)

                   call memfree ( mlbddc%A_rr_trans_inv_C_r_T,__FILE__,__LINE__)
                end if

                if (mlbddc%symm == symm_false .or. & 
                     (mlbddc%symm == symm_true .and. (mlbddc%sign == indefinite .or. mlbddc%sign == unknown)) ) then
                   call memfree ( mlbddc%ipiv,__FILE__,__LINE__)
                   if (mlbddc%projection == petrov_galerkin )  then 
                      call memfree ( mlbddc%ipiv_trans,__FILE__,__LINE__)
                   end if

                end if

                if (mlbddc%schur_edge_lag_mult == compute_from_scratch ) then
                   call memfree ( mlbddc%S_rr_neumann, & 
                        & 'par_preconditioner_dd_bddc_free::mlbddc%S_rr_neumann' )
                   
                   call memfree ( mlbddc%A_rr_inv_C_r_T_neumann, & 
                        & 'par_preconditioner_dd_bddc_free::mlbddc%A_rr_inv_C_r_T_neumann' )
                   
                   if (mlbddc%projection == petrov_galerkin )  then 
                      call memfree ( mlbddc%S_rr_trans_neumann, & 
                           & 'par_preconditioner_dd_bddc_free::mlbddc%S_rr_neumann' )

                      call memfree ( mlbddc%A_rr_trans_inv_C_r_T_neumann, & 
                           & 'par_preconditioner_dd_bddc_free::mlbddc%A_rr_inv_C_r_T_neumann' )   
                   end if

                   if (mlbddc%symm == symm_false .or. & 
                        (mlbddc%symm == symm_true .and. & 
                        (mlbddc%sign == indefinite .or. mlbddc%sign == unknown)) ) then
                      call memfree ( mlbddc%ipiv_neumann, & 
                           & 'par_preconditioner_dd_bddc_free::mlbddc%ipiv' )
                      if (mlbddc%projection == petrov_galerkin )  then 
                         call memfree ( mlbddc%ipiv_trans_neumann, & 
                              & 'par_preconditioner_dd_bddc_free::mlbddc%ipiv' )
                      end if
                   end if
                end if

             end if
          end if
          ! END FINE-GRID PROBLEM DUTIES
       end if


       if(mlbddc%co_sys_sol_strat == serial_gather) then
          if ( i_am_coarse_task ) then
             if ( mlbddc%internal_problems == handled_by_bddc_module) then
                call preconditioner_free ( preconditioner_free_values, mlbddc%M_c )
             end if
             call matrix_free ( mlbddc%A_c, free_only_values )
          end if
       else if(mlbddc%co_sys_sol_strat == recursive_bddc) then
          assert(mlbddc%p_mat%p_env%num_levels>2) 
          if ( i_am_coarse_task .or. i_am_higher_level_task ) then
             ! Recursively call bddc_free
             call par_preconditioner_dd_mlevel_bddc_free(mlbddc%p_M_c, free_only_values )

             ! Free coarse matrix
             call par_matrix_free( mlbddc%p_mat_c, free_only_values)
          end if
       end if

    end if

  end subroutine par_preconditioner_dd_mlevel_bddc_free

  !=================================================================================================
  recursive subroutine par_preconditioner_dd_mlevel_bddc_ass_struct ( p_mat, mlbddc ) 
use mpi
    implicit none

    ! Parameters 
    type(par_matrix_t)                , target, intent(in)    :: p_mat
    type(par_preconditioner_dd_mlevel_bddc_t), target, intent(inout) :: mlbddc
 
    type (mesh_t)  :: c_mesh
    logical          :: i_am_coarse_task, i_am_fine_task, i_am_higher_level_task
    !integer(ip)      :: iam

    ! The routine requires the partition/context info
    assert ( associated(p_mat%p_env) )
    assert ( p_mat%p_env%created )
    assert ( p_mat%p_env%num_levels > 1 ) 

    assert ( associated(mlbddc%p_mat%p_env%w_context) )
    assert ( mlbddc%p_mat%p_env%w_context%created .eqv. .true. )
    assert ( associated(mlbddc%p_mat%p_env%p_context) )
    assert ( mlbddc%p_mat%p_env%p_context%created .eqv. .true. )
    assert ( associated(mlbddc%p_mat%p_env%q_context) )
    assert ( mlbddc%p_mat%p_env%q_context%created .eqv. .true. )
    ! Check appropriate assignment to context: 1) I'm in w_context and 2) I'm in p_context OR in q_context but not in both
    assert ( mlbddc%p_mat%p_env%w_context%iam >= 0)
    assert ( (mlbddc%p_mat%p_env%p_context%iam >=0 .and. mlbddc%p_mat%p_env%q_context%iam < 0) .or. (mlbddc%p_mat%p_env%p_context%iam < 0 .and. mlbddc%p_mat%p_env%q_context%iam >= 0))
    assert ( mlbddc%c_context%created .eqv. .true. )
    assert ( mlbddc%d_context%created .eqv. .true. )
    assert ( mlbddc%g_context%created .eqv. .true. )
    assert ( associated(mlbddc%p_mat%p_env%b_context) )
    assert ( mlbddc%p_mat%p_env%b_context%created .eqv. .true. )

    mlbddc%p_mat => p_mat
    
    ! Which duties I have?
    i_am_fine_task = (mlbddc%p_mat%p_env%p_context%iam >= 0)
    i_am_coarse_task = (mlbddc%c_context%iam >= 0)
    i_am_higher_level_task = (mlbddc%d_context%iam >= 0)
    ! if(i_am_fine_task) iam = mlbddc%p_mat%p_env%p_context%iam
    ! if(i_am_coarse_task) iam = mlbddc%c_context%iam
    ! if(i_am_higher_level_task) iam = mlbddc%d_context%iam

    if ( i_am_fine_task ) then
       ! BEG. FINE-GRID PROBLEM DUTIES
       ! Local Computations (no communication at all)
       call identify_coarse_dofs_and_compute_perm_iperm_naive (mlbddc)
       ! END FINE-GRID PROBLEM DUTIES
    else if ( i_am_coarse_task ) then  ! TODO: needed if i_am_higher_label_task too?
       mlbddc%nl_coarse=0           !       but if done, then calls to count_global_coarse
       mlbddc%nl_corners=0          !       can be unified...
       mlbddc%nl_corners_dofs=0
       mlbddc%nl_edges=0
       mlbddc%nl_edges_dofs=0
    end if

    ! TODO: check timers everywhere (already DONE by AFMH)
    if ( temp_fine_coarse_grid_overlap ) then 
       if ( .not. i_am_higher_level_task ) then
          call par_timer_start ( mlbddc%timer_coll_assprec  )
       end if
    end if

    if( .not. i_am_higher_level_task ) then
       if(mlbddc%co_sys_sol_strat == serial_gather) then
          ! BEG. COARSE-GRID PROBLEM DUTIES
          ! Global Computations related to the coarse grid system (communication)
          call gather_coarse_mesh (mlbddc, mlbddc%f_mesh_c)
          
          if ( i_am_coarse_task ) call generate_coarse_graph (mlbddc)
          ! END COARSE-GRID PROBLEM DUTIES
       else if (mlbddc%co_sys_sol_strat == recursive_bddc) then
          ! BEG. COARSE-GRID PROBLEM DUTIES
          ! Global Computations related to the coarse grid system (communication)
          call gather_coarse_mesh (mlbddc, c_mesh)
          
          call generate_coarse_partition_and_graph (c_mesh, mlbddc)
          ! END COARSE-GRID PROBLEM DUTIES
       end if
    end if

    
!!$    if ( temp_fine_coarse_grid_overlap ) then
!!$       if ( .not. i_am_higher_level_task ) then
!!$          call par_timer_stop ( mlbddc%timer_coll_assprec  )
!!$          call psb_barrier ( mlbddc%g_context%icontxt )
!!$          if ( i_am_fine_task ) then
!!$             call par_timer_start ( mlbddc%timer_assprec_ov_fine  ) 
!!$          else if ( i_am_coarse_task ) then
!!$             call par_timer_start ( mlbddc%timer_assprec_ov_coarse ) 
!!$          end if
!!$       end if
!!$    end if
    
    if ( i_am_fine_task ) then

       if ( temp_fine_coarse_grid_overlap ) then
          call par_timer_start ( mlbddc%timer_assprec_ov_fine_header )
       end if

       ! BEG. FINE-GRID PROBLEM DUTIES
       if (mlbddc%nn_sys_sol_strat == corners_rest_part_solve_expl_schur ) then 
          call extract_graph_A_rr (p_mat, mlbddc)
       else if (  mlbddc%nn_sys_sol_strat == direct_solve_constrained_problem ) then
          call augment_graph_with_constraints ( p_mat, mlbddc )
       end if
       
       call matrix_graph ( mlbddc%A_rr_gr, mlbddc%A_rr)

          if (mlbddc%projection == petrov_galerkin) then
             call matrix_graph ( mlbddc%A_rr_gr, mlbddc%A_rr_trans)
          end if

       if ( mlbddc%internal_problems == handled_by_bddc_module) then
             call preconditioner_symbolic( mlbddc%A_rr , mlbddc%M_rr )

          if (mlbddc%projection == petrov_galerkin) then
             call preconditioner_symbolic( mlbddc%A_rr_trans , mlbddc%M_rr_trans )
          end if

       end if
       if ( mlbddc%unknowns == all_unknowns ) then
          call operator_dd_ass_struct (p_mat%f_matrix, mlbddc%A_II_inv ) 
       end if
       ! END FINE-GRID PROBLEM DUTIES

       if ( temp_fine_coarse_grid_overlap ) then
          call par_timer_stop ( mlbddc%timer_assprec_ov_fine_header )
       end if

    else  ! if ( i_am_coarse_task .or. i_am_higher_level_task ) then

       ! BEG. COARSE-GRID PROBLEM DUTIES
       if (  mlbddc%co_sys_sol_strat == serial_gather ) then  ! There are only coarse tasks
          call matrix_graph ( mlbddc%A_c_gr, mlbddc%A_c)
          if ( mlbddc%internal_problems == handled_by_bddc_module) then
             call preconditioner_symbolic( mlbddc%A_c , mlbddc%M_c )
          end if
       else if (  mlbddc%co_sys_sol_strat == recursive_bddc ) then
          ! Assign coarse matrix graph (the partition)
          call par_matrix_graph ( mlbddc%p_graph_c, mlbddc%p_mat_c)
          call par_preconditioner_dd_mlevel_bddc_ass_struct ( mlbddc%p_mat_c, mlbddc%p_M_c ) 
       end if
       ! END COARSE-GRID PROBLEM DUTIES
    end if

  end subroutine par_preconditioner_dd_mlevel_bddc_ass_struct

  subroutine augment_graph_with_constraints (p_mat, mlbddc)
    implicit none

    ! Parameters 
    type(par_matrix_t)         , intent(in)     :: p_mat
    type(par_preconditioner_dd_mlevel_bddc_t), intent(inout)  :: mlbddc

    integer :: info

    mlbddc%A_rr_gr%nv    = mlbddc%p_mat%dof_dist%nl  + &
                                 & mlbddc%nl_corners             + &
                                 & mlbddc%nl_edges              
 
    mlbddc%A_rr_gr%nv2   = mlbddc%A_rr_gr%nv
    mlbddc%A_rr_gr%type  = p_mat%f_matrix%gr%type

    call memalloc ( mlbddc%A_rr_gr%nv+1, mlbddc%A_rr_gr%ia, __FILE__,__LINE__ )

    ! Count neighbours
    call count_graph_augment_graph_with_constraints ( p_mat%f_matrix%gr%type ,  & 
                                                      p_mat%f_matrix%gr%nv   ,  & 
                                                      p_mat%f_matrix%gr%ia   ,  & 
                                                      mlbddc%A_rr_gr%nv  , &
                                                      mlbddc%p_mat%dof_dist%nl, &
                                                      mlbddc%nl_coarse, &
                                                      mlbddc%coarse_dofs, &
                                                      mlbddc%p_mat%dof_dist%max_nparts, &
                                                      mlbddc%p_mat%dof_dist%omap%nl, &
                                                      mlbddc%p_mat%dof_dist%lobjs, &
                                                      mlbddc%A_rr_gr%ia  )

    call memalloc ( mlbddc%A_rr_gr%ia(mlbddc%A_rr_gr%nv+1)-1, mlbddc%A_rr_gr%ja,  __FILE__,__LINE__ )

    ! List neighbours 
    call list_graph_augment_graph_with_constraints ( p_mat%f_matrix%gr%type    , & 
                                                     p_mat%f_matrix%gr%nv      , & 
                                                     p_mat%f_matrix%gr%ia      , & 
                                                     p_mat%f_matrix%gr%ja      , &
                                                     mlbddc%A_rr_gr%nv     , & 
                                                     mlbddc%A_rr_gr%ia     , &
                                                     mlbddc%p_mat%dof_dist%nl, &
                                                     mlbddc%nl_coarse, &
                                                     mlbddc%coarse_dofs, &
                                                     mlbddc%p_mat%dof_dist%max_nparts, &
                                                     mlbddc%p_mat%dof_dist%omap%nl, &
                                                     mlbddc%p_mat%dof_dist%lobjs, &
                                                     mlbddc%A_rr_gr%ja  )

!!$    if ( debug_verbose_level_2 ) then
!!$       call graph_print (6, p_mat%f_matrix%gr )
!!$       call graph_print (6, mlbddc%A_rr_gr)
!!$    end if
  end subroutine augment_graph_with_constraints

  subroutine count_graph_augment_graph_with_constraints ( gtype, anv, aia, arr_nv,  &
                                                          nl, nl_coarse, & 
                                                          coarse_dofs, max_nparts, nobjs, lobjs, & 
                                                          arr_ia )
    implicit none
    ! Parameters
    integer (ip), intent(in)  :: gtype, anv, arr_nv
    integer (ip), intent(in)  :: aia (anv+1)
    integer (ip), intent(in)  :: nl
    integer (ip), intent(in)  :: nl_coarse
    integer (ip), intent(in)  :: coarse_dofs(nl_coarse) 
    integer (ip), intent(in)  :: max_nparts, nobjs
    integer (ip), intent(in)  :: lobjs(max_nparts+4, nobjs)
    integer (ip), intent(out) :: arr_ia (arr_nv+1)
    
    ! Locals
    integer(ip)  :: idarr, ida
    integer (ip) :: i, j, iobj, jd, j1, j2

    
    ! Count # of elements on each 
    ! row of a and map them to the
    ! corresponding rows of a_rr
    arr_ia = 0
    do idarr=1, anv
       arr_ia(idarr+1) = aia(idarr+1)-aia(idarr)
    end do

    if ( gtype == csr ) then
       do i=1, nl_coarse 
          iobj   = coarse_dofs (i)
          ! Corner as a degenerated edge or
          ! edge/face
          j1     = lobjs(2,iobj)
          j2     = lobjs(3,iobj)
         
          do jd = j1, j2
             ! new entry (i+nl,jd)
             arr_ia(jd+1) = arr_ia(jd+1)+1

             ! new entry (jd,i+nl)
             arr_ia(i+nl+1) = arr_ia(i+nl+1)+1
          end do
        end do
    else if ( gtype == csr_symm ) then
       do i=1, nl_coarse 
          iobj   = coarse_dofs (i)

          ! Corner as a degenerated edge or
          ! edge/face
          j1     = lobjs(2,iobj)
          j2     = lobjs(3,iobj)
          
          do jd = j1, j2
             ! new entry (id+nl,jd)
             arr_ia(jd+1) = arr_ia(jd+1)+1
          end do

       end do
    end if

    do idarr=anv+1,arr_nv
       arr_ia(idarr+1) = arr_ia(idarr+1) + 1  
       ! Add explicit zeros on the diagonal
    end do

    ! Transform length to header
    arr_ia(1)=1
    do idarr=1, arr_nv
       arr_ia(idarr+1)=arr_ia(idarr)+arr_ia(idarr+1)
    end do

    ! write(*,*) 'XXX', arr_ia
    
  end subroutine count_graph_augment_graph_with_constraints

  subroutine list_graph_augment_graph_with_constraints ( gtype, anv, aia, aja, arr_nv, arr_ia, &
                                                         nl, nl_coarse, & 
                                                         coarse_dofs, max_nparts, nobjs, lobjs, & 
                                                         arr_ja )
    implicit none
    ! Parameters
    integer (ip) , intent(in)    :: gtype, anv, arr_nv
    integer (ip) , intent(in)    :: aia (anv+1)
    integer (ip) , intent(inout) :: arr_ia(arr_nv+1)
    integer (ip) , intent(in)    :: aja (aia(anv+1)-1)
    integer (ip) , intent(in)    :: nl
    integer (ip) , intent(in)    :: nl_coarse
    integer (ip) , intent(in)    :: coarse_dofs (nl_coarse) 
    integer (ip) , intent(in)    :: max_nparts, nobjs
    integer (ip) , intent(in)    :: lobjs(max_nparts+4, nobjs)
    integer (ip) , intent(out)   :: arr_ja (arr_ia(arr_nv+1)-1)

    ! Locals
    integer(ip)  :: idarr, nz
    integer (ip) :: i, j, iobj, jd, j1, j2

    ! List elements on each row of 
    ! a and map them to the corresponding 
    ! rows of a_rr
    do idarr=1, anv
       nz = aia(idarr+1) - aia(idarr)  
       arr_ja ( arr_ia(idarr):arr_ia(idarr)+nz-1 ) = aja ( aia(idarr):aia(idarr+1)-1 )

       arr_ia ( idarr ) = arr_ia( idarr ) + nz
    end do

    
    if ( gtype == csr ) then
       do i=1, nl_coarse
          iobj = coarse_dofs (i)

          ! Corner as a degenerated edge or
          ! edge/face
          j1 = lobjs(2,iobj)
          j2 = lobjs(3,iobj)
          
          do jd = j1, j2
             ! new entry (jd,i+n)
             arr_ja( arr_ia(jd) ) = i + nl 
             arr_ia( jd )         = arr_ia(jd) + 1

             ! new entry (i+nl,jd)
             arr_ja ( arr_ia(i+nl) ) = jd
             arr_ia ( i+nl )         = arr_ia ( i+nl ) + 1 
          end do

       end do
    else if ( gtype == csr_symm ) then ! Add explicit zeros on the diagonal
       do i=1, nl_coarse 
          iobj   = coarse_dofs (i)

          ! Corner as a degenerated edge or
          ! edge/face
          j1     = lobjs(2,iobj)
          j2     = lobjs(3,iobj)  

          do jd = j1, j2
             ! new entry (jd,i+nl)
             arr_ja( arr_ia(jd) ) = i + nl
             arr_ia(jd)           = arr_ia(jd) + 1
          end do
       end do
    end if

    ! Add explicit zeros on the diagonal
    do idarr=anv+1,arr_nv
       arr_ja( arr_ia(idarr) ) = idarr
       arr_ia(idarr) = arr_ia(idarr) + 1  
    end do

    ! Recover arr_ia
    do idarr = arr_nv + 1, 2, -1
       arr_ia (idarr) = arr_ia (idarr-1)
    end do
    arr_ia (idarr) = 1
 
    ! write (*,*) 'YYY', arr_ia    

  end subroutine list_graph_augment_graph_with_constraints

  subroutine identify_coarse_dofs_and_compute_perm_iperm_naive (mlbddc)
    implicit none

    ! Parameters
    type(par_preconditioner_dd_mlevel_bddc_t), intent(inout) :: mlbddc

    ! Locals
    integer (ip) :: iobj, icorner, iedge

    call count_coarse_dofs ( mlbddc%p_mat%dof_dist%max_nparts  , &
                             mlbddc%p_mat%dof_dist%omap%nl     , &
                             mlbddc%p_mat%dof_dist%lobjs       , &
                             mlbddc%kind_coarse_dofs          , &
                             mlbddc%ndime                     , &
                             mlbddc%nl_corners                , &
                             mlbddc%nl_corners_dofs           , &
                             mlbddc%nl_edges                  , &
                             mlbddc%nl_edges_dofs           ) 

    mlbddc%nl_coarse = mlbddc%nl_corners + & 
         mlbddc%nl_edges

    ! List coarse dofs
    call memalloc (  mlbddc%nl_coarse,     mlbddc%coarse_dofs   ,   __FILE__,__LINE__ )
    
    if ( mlbddc%nn_sys_sol_strat == corners_rest_part_solve_expl_schur ) then 
       
       call memalloc ( mlbddc%p_mat%dof_dist%nl,  mlbddc%perm,__FILE__,__LINE__ )
       call memalloc ( mlbddc%p_mat%dof_dist%nl, mlbddc%iperm,__FILE__,__LINE__ )
       call list_coarse_dofs ( mlbddc%p_mat%dof_dist%max_nparts, &
                               mlbddc%p_mat%dof_dist%omap%nl   , &
                               mlbddc%p_mat%dof_dist%lobjs     , &
                               mlbddc%kind_coarse_dofs        , &
                               mlbddc%ndime                   , &
                               mlbddc%p_mat%dof_dist%nl   , &
                               mlbddc%nl_corners              , &
                               mlbddc%nl_corners_dofs         , &
                               mlbddc%nl_edges                , &
                               mlbddc%nl_edges_dofs           , &
                               mlbddc%coarse_dofs             , &  
                               mlbddc%iperm                   , &
                               mlbddc%perm                    )

    else if ( mlbddc%nn_sys_sol_strat == direct_solve_constrained_problem ) then
       call list_coarse_dofs ( mlbddc%p_mat%dof_dist%max_nparts, &
                               mlbddc%p_mat%dof_dist%omap%nl   , &
                               mlbddc%p_mat%dof_dist%lobjs     , &
                               mlbddc%kind_coarse_dofs        , &
                               mlbddc%ndime                   , &
                               mlbddc%p_mat%dof_dist%nl   , &
                               mlbddc%nl_corners              , &
                               mlbddc%nl_corners_dofs         , &
                               mlbddc%nl_edges                , &
                               mlbddc%nl_edges_dofs           , &
                               mlbddc%coarse_dofs  )
    end if

    if ( debug_verbose_level_1 ) then
       write (*,*)   'Local Number of coarse nodes:', mlbddc%nl_coarse
       write (*,*)   'Local Number of corners:', mlbddc%nl_corners
       write (*,*)   'Local Sum of corners dofs:', mlbddc%nl_corners_dofs
       write (*,*)   'Local Number of edges:', mlbddc%nl_edges
       write (*,*)   'Local Sum of edges dofs:', mlbddc%nl_edges_dofs
       write (*,*)   'Local List of coarse dofs:', mlbddc%coarse_dofs
       ! write (*, *)  'perm                  :', mlbddc%perm
       ! write (*, *)  'iperm                 :', mlbddc%iperm
    end if

  end subroutine identify_coarse_dofs_and_compute_perm_iperm_naive

  subroutine count_coarse_dofs ( max_nparts, nobjs, lobjs, &
                               & kind_coarse_dofs, ndime,  &
                               & ncorners, ncornersdofs,   &
                               & nedges  , nedgesdofs )
    implicit none

    ! Parameters
    integer (ip), intent(in)  :: max_nparts, nobjs
    integer (ip), intent(in)  :: lobjs (max_nparts+4, nobjs)
    integer (ip), intent(in)  :: kind_coarse_dofs, ndime
    integer (ip), intent(out) :: ncorners, ncornersdofs
    integer (ip), intent(out) :: nedges  , nedgesdofs

    ! Locals
    integer (ip) :: iobj, idof, inode
    integer (ip) :: d 

    ncorners     = 0
    ncornersdofs = 0
    nedges       = 0
    nedgesdofs   = 0

    ! The code does not work when the interior object
    ! has only one node !!!!
    assert (lobjs(3,1)-lobjs(2,1)/=0) 

    do iobj=1,nobjs
       ! write (*,*) 'PPP', iobj, lobjs(1,iobj), lobjs(2,iobj), lobjs(3,iobj)
       if ( lobjs(3,iobj)-lobjs(2,iobj) == 0 ) then
          ! write (*,*) 'PPP', iobj, lobjs(1,iobj), lobjs(3,iobj), lobjs(4,iobj), lobjs(5:,iobj)
          ! write (*,*) 'XXX', iobj, dirichlet
          if ( kind_coarse_dofs == corners                  .or. & 
               kind_coarse_dofs == corners_and_edges        .or. & 
               kind_coarse_dofs == corners_edges_and_faces  .or. & 
               kind_coarse_dofs == corners_and_faces ) then
             ncorners = ncorners + 1 
             ncornersdofs = ncornersdofs + 1
          end if
       else if ( kind_coarse_dofs == corners_and_edges       .or. & 
                 kind_coarse_dofs == corners_edges_and_faces .or. & 
                 kind_coarse_dofs == edges                   .or. &
                 kind_coarse_dofs == faces                   .or. &
                 kind_coarse_dofs == edges_and_faces         .or. &
                 kind_coarse_dofs == corners_and_faces ) then
          if ( kind_coarse_dofs == corners_and_edges .or. kind_coarse_dofs == edges ) then
             d = ndime-1
          else if ( kind_coarse_dofs == corners_edges_and_faces .or. kind_coarse_dofs == edges_and_faces ) then
             d = ndime-2 
          end if

          if ( (kind_coarse_dofs /= faces .and. kind_coarse_dofs /= corners_and_faces .and. lobjs(4,iobj)  > d) .or. &
               ((kind_coarse_dofs == faces.or.kind_coarse_dofs == corners_and_faces) .and. lobjs(4,iobj) == 2) ) then
                nedges = nedges + 1
                nedgesdofs = nedgesdofs + (lobjs(3,iobj)-lobjs(2,iobj)+1)
             end if
          end if
    end do

  end subroutine count_coarse_dofs

  subroutine list_coarse_dofs ( max_nparts, nobjs, lobjs, &
                              & kind_coarse_dofs, ndime , &
                              & n,                        &
                              & ncorners, ncornersdofs  , &
                              & nedges  , nedgesdofs    , &
                              & coarse_dofs             , &
                              & iperm, perm )                
    implicit none

    ! Parameters
    integer (ip), intent(in)            :: max_nparts, nobjs
    integer (ip), intent(in)            :: lobjs (max_nparts+4, nobjs)
    integer (ip), intent(in)            :: kind_coarse_dofs, ndime
    integer (ip), intent(in)            :: n
    integer (ip), intent(in)            :: ncorners, ncornersdofs
    integer (ip), intent(in)            :: nedges  , nedgesdofs
    integer (ip), intent(out)           :: coarse_dofs(ncorners+nedges)
    integer (ip), optional, intent(out) :: iperm(n), perm(n)

    ! Locals
    integer (ip) :: iobj, icorner, iedge
    integer (ip) :: idcorners, idedges, idrest
    integer (ip) :: ipoin, id, base, inode
    integer (ip) :: d


    if (present(iperm)) then
       icorner = 1
       iedge   = ncorners + 1

       idcorners = 1
       idedges   = ncornersdofs + 1
       idrest    = ncornersdofs + nedgesdofs + 1


       do iobj=1,nobjs
          if ( lobjs(3,iobj) - lobjs(2,iobj) == 0 ) then
             if ( kind_coarse_dofs == corners                  .or. & 
                  kind_coarse_dofs == corners_and_edges        .or. & 
                  kind_coarse_dofs == corners_edges_and_faces  .or. & 
                  kind_coarse_dofs == corners_and_faces ) then

                coarse_dofs(icorner) = iobj 
                icorner = icorner + 1
                do ipoin=lobjs(2,iobj),lobjs(3,iobj)
                   base = ipoin
                   iperm ( idcorners ) = base
                   base      = base      + 1
                   idcorners = idcorners + 1
                end do
              else
                do ipoin=lobjs(2,iobj),lobjs(3,iobj)
                   base = ipoin
                   iperm ( idrest ) = base
                   base   = base    + 1
                   idrest = idrest  + 1
                end do
             end if
          else if ( kind_coarse_dofs == corners_and_edges       .or. & 
                    kind_coarse_dofs == corners_edges_and_faces .or. & 
                    kind_coarse_dofs == edges                   .or. &
                    kind_coarse_dofs == faces                   .or. &
                    kind_coarse_dofs == edges_and_faces         .or. &
                    kind_coarse_dofs == corners_and_faces ) then

             if ( kind_coarse_dofs == corners_and_edges  .or. kind_coarse_dofs == edges ) then
                d = ndime-1
             else if ( kind_coarse_dofs == corners_edges_and_faces .or. kind_coarse_dofs == edges_and_faces ) then
                d = ndime-2 
             end if

             if ( (kind_coarse_dofs /= faces .and. kind_coarse_dofs /= corners_and_faces .and. lobjs(4,iobj)  > d) .or. &
                 ((kind_coarse_dofs == faces.or.kind_coarse_dofs == corners_and_faces) .and. lobjs(4,iobj) == 2) ) then
                  coarse_dofs(iedge) = iobj 
                  iedge = iedge + 1
                  do ipoin=lobjs(2,iobj),lobjs(3,iobj)
                     base = ipoin
                     iperm ( idedges ) = base
                     base    = base    + 1
                     idedges = idedges + 1
                  end do
             else
                do ipoin=lobjs(2,iobj),lobjs(3,iobj)
                   base = ipoin
                   iperm ( idrest ) = base
                   base   = base    + 1
                   idrest = idrest  + 1
                end do
             end if
          else 
             do ipoin=lobjs(2,iobj),lobjs(3,iobj)
                base = ipoin
                iperm ( idrest ) = base
                base   = base    + 1
                idrest = idrest  + 1
             end do
          end if
       end do

    else
       icorner = 1
       iedge   = ncorners + 1

       do iobj=1,nobjs
          if ( lobjs(3,iobj) - lobjs(2,iobj) == 0 ) then
             if ( kind_coarse_dofs == corners                  .or. & 
                  kind_coarse_dofs == corners_and_edges        .or. & 
                  kind_coarse_dofs == corners_edges_and_faces  .or. & 
                  kind_coarse_dofs == corners_and_faces ) then
                coarse_dofs(icorner) = iobj 
                icorner = icorner + 1
             end if
          else if ( kind_coarse_dofs == corners_and_edges       .or. & 
                    kind_coarse_dofs == corners_edges_and_faces .or. &
                    kind_coarse_dofs == edges                   .or. &
                    kind_coarse_dofs == faces                   .or. &
                    kind_coarse_dofs == edges_and_faces         .or. &
                    kind_coarse_dofs == corners_and_faces ) then
             if ( kind_coarse_dofs == corners_and_edges .or. kind_coarse_dofs == edges  ) then
                d = ndime-1
             else if ( kind_coarse_dofs == corners_edges_and_faces .or. kind_coarse_dofs == edges_and_faces ) then
                d = ndime-2 
             end if
             
             if ( (kind_coarse_dofs /= faces .and. kind_coarse_dofs /= corners_and_faces .and. lobjs(4,iobj)  > d) .or. &
                 ((kind_coarse_dofs == faces.or.kind_coarse_dofs == corners_and_faces) .and. lobjs(4,iobj) == 2) ) then
                coarse_dofs(iedge) = iobj 
                iedge = iedge + 1
             end if
          end if
       end do
    end if

    if (present(perm)) then
       ! Compute perm for its inverse
       do id=1, n
          perm( iperm(id) ) = id
       end do
    end if

  end subroutine list_coarse_dofs

  subroutine extract_graph_A_rr (p_mat, mlbddc)
    implicit none

    ! Parameters 
    type(par_matrix_t)         , intent(in)                    :: p_mat
    type(par_preconditioner_dd_mlevel_bddc_t), intent(inout)                 :: mlbddc

    mlbddc%A_rr_gr%nv    = mlbddc%p_mat%dof_dist%nl - mlbddc%nl_corners_dofs
    mlbddc%A_rr_gr%nv2   = mlbddc%A_rr_gr%nv
    mlbddc%A_rr_gr%type  = p_mat%f_matrix%gr%type

    call memalloc ( mlbddc%A_rr_gr%nv+1, mlbddc%A_rr_gr%ia,  __FILE__,__LINE__ )

    ! Count neighbours
    call count_graph_A_rr ( mlbddc%nl_corners_dofs, &
         p_mat%f_matrix%gr%type     , & 
         p_mat%f_matrix%gr%nv       , & 
         p_mat%f_matrix%gr%ia       , & 
         p_mat%f_matrix%gr%ja       , &
         mlbddc%A_rr_gr%nv     , & 
         mlbddc%perm           , &
         mlbddc%iperm          , &   
         mlbddc%A_rr_gr%ia  )

    call memalloc ( mlbddc%A_rr_gr%ia(mlbddc%A_rr_gr%nv+1)-1, mlbddc%A_rr_gr%ja,  __FILE__,__LINE__ )

    ! List neighbours 
    call list_graph_A_rr  ( mlbddc%nl_corners_dofs, & 
         p_mat%f_matrix%gr%type     , & 
         p_mat%f_matrix%gr%nv       , & 
         p_mat%f_matrix%gr%ia       , & 
         p_mat%f_matrix%gr%ja       , &
         mlbddc%A_rr_gr%nv     , & 
         mlbddc%A_rr_gr%ia     , & 
         mlbddc%perm           , &
         mlbddc%iperm          , &       
         mlbddc%A_rr_gr%ja  )

!!$    if ( debug_verbose_level_2 ) then
!!$       ! call graph_print (6, p_mat%f_matrix%gr )
!!$       call graph_print (6, mlbddc%A_rr_gr)
!!$    end if

  end subroutine extract_graph_A_rr

  subroutine count_graph_A_rr ( nl_corners_dofs, gtype, anv, aia, aja, arrnv, perm, iperm, arria )
    implicit none
    ! Parameters
    integer (ip) , intent(in)  :: nl_corners_dofs
    integer (ip) , intent(in)  :: gtype, anv, arrnv
    integer (ip) , intent(in)  :: aia   (anv+1)
    integer (ip) , intent(in)  :: aja   (aia(anv+1)-1)
    integer (ip) , intent(in)  :: perm  (anv)
    integer (ip) , intent(in)  :: iperm (anv)
    integer (ip), intent(out)  :: arria(arrnv+1)

    ! Locals
    integer(ip) :: idarr, ida, jdarr, jda, jdapos


    arria = 0
    do idarr=1, arrnv
       ida = iperm(idarr+nl_corners_dofs)
       do jdapos = aia(ida), aia(ida+1)-1
          jdarr = perm( aja(jdapos) )
          if ( jdarr > nl_corners_dofs ) then
             jdarr = jdarr -nl_corners_dofs
             if ( gtype == csr ) then
                arria(idarr+1) = arria(idarr+1) + 1 
             else if ( gtype == csr_symm ) then
                if ( jdarr >= idarr ) then
                   arria(idarr+1) = arria(idarr+1) + 1 
                else ! jdarr < idarr
                   arria(jdarr+1) = arria(jdarr+1) + 1                     
                end if
             end if
          end if
       end do
    end do

    arria(1)=1
    do idarr=1, arrnv
       arria(idarr+1)=arria(idarr)+arria(idarr+1)
    end do

  end subroutine count_graph_A_rr

  subroutine list_graph_A_rr ( nl_corners_dofs, gtype, anv, aia, aja, arrnv, arria, perm, iperm, arrja)
    implicit none
    ! Parameters
    integer (ip) , intent(in)    :: nl_corners_dofs
    integer (ip) , intent(in)    :: gtype, anv, arrnv
    integer (ip) , intent(in)    :: aia   (anv+1)
    integer (ip) , intent(inout) :: arria(arrnv+1)
    integer (ip) , intent(in)    :: aja   (aia(anv+1)-1)
    integer (ip) , intent(in)    :: perm  (anv)
    integer (ip) , intent(in)    :: iperm (anv)
    integer (ip) , intent(out)   :: arrja (arria(arrnv+1)-1)

    ! Locals
    integer(ip)                  :: idarr, ida, jdarr, jda, jdapos
    integer(ip)                  :: k1, k2

    do idarr=1, arrnv
       ida = iperm(idarr+nl_corners_dofs)
       k1  = arria (idarr)
       k2  = arria (idarr+1)-1
       do jdapos = aia(ida), aia(ida+1)-1
          jdarr = perm( aja(jdapos) )
          if ( jdarr > nl_corners_dofs ) then
             jdarr = jdarr -nl_corners_dofs
             if ( gtype == csr ) then
                arrja ( arria(idarr) ) = jdarr
                arria (idarr) = arria(idarr) + 1
             else if ( gtype == csr_symm ) then
                if ( jdarr >= idarr ) then
                   arrja ( arria(idarr) ) = jdarr
                   arria(idarr) = arria(idarr) + 1 
                else ! jdarr < idarr
                   arrja ( arria(jdarr) ) = idarr
                   arria(jdarr) = arria(jdarr) + 1            
                end if
             end if
          end if
       end do

       if ( gtype == csr ) then
          ! Order increasingly column identifiers of current row 
          ! using heap sort algorithm
          call sort( k2-k1+1, arrja(k1:k2) )
       end if
    end do

    do idarr=arrnv+1, 2, -1
       arria(idarr) = arria(idarr-1)
    end do
    arria(idarr) = 1

    if ( gtype == csr_symm ) then
       do idarr=1, arrnv
          k1  = arria (idarr)
          k2  = arria (idarr+1)-1
          ! Order increasingly column identifiers of current row 
          ! using heap sort algorithm
          call sort( k2-k1+1, arrja(k1:k2) )
       end do
    end if

  end subroutine list_graph_A_rr

  subroutine gather_coarse_mesh(mlbddc, c_mesh)
    ! Formerly compute_global_ordering_coarse_grid_system (mlbddc)
    implicit none

    ! Parameters 
    type(par_preconditioner_dd_mlevel_bddc_t), intent(inout) :: mlbddc
    type(mesh_t), intent(out) :: c_mesh
    integer(ip), allocatable  :: idum (:)
    logical                   :: i_am_fine_task, i_am_coarse_task, i_am_higher_level_task
    integer                   :: info, iam
    integer(ip)               :: i

    ! Which duties I have?
    i_am_fine_task = (mlbddc%p_mat%p_env%p_context%iam >= 0)
    i_am_coarse_task = (mlbddc%c_context%iam >= 0)
    i_am_higher_level_task = (mlbddc%d_context%iam >= 0)
    if(i_am_fine_task) iam = mlbddc%p_mat%p_env%p_context%iam
    if(i_am_coarse_task) iam = mlbddc%c_context%iam
    if(i_am_higher_level_task) iam = mlbddc%d_context%iam

    if( i_am_fine_task) then
       
       ! Both omap%ng and max_nparts are required next in the coarse-grid
       ! tasks. TRANSFER these data!
       call transfer_coarse_int ( mlbddc%g_context%icontxt, &
                                  mlbddc%g_context%iam, &
                                  mlbddc%g_context%np, &
                                  mlbddc%p_mat%dof_dist%omap%ng )

       call transfer_coarse_int ( mlbddc%g_context%icontxt, &
                                  mlbddc%g_context%iam, &
                                  mlbddc%g_context%np, &
                                  mlbddc%p_mat%dof_dist%max_nparts )

       ! write (*,*) 'I am fine task', mlbddc%g_context%iam, mlbddc%g_context%np
       ! call psb_barrier ( mlbddc%p_mat%p_env%p_context%icontxt)

       call snd_coarse_int_ip ( mlbddc%g_context%icontxt, & 
                                mlbddc%g_context%np, &
                                mlbddc%nl_coarse)

       ! if( mlbddc%pad_collectives == pad ) then
       ! AFM: IMPORTANT NOTE !!!! You need to call all_reduce_coarse_mesh_nnode
       ! EVEN with nopad mode, as c_mesh%nnode is later required for the generation 
       ! of the primal graph in create_graph_from_mesh_matrix (mesh_graph.f90)
       call allreduce_coarse_mesh_nnode(  mlbddc%g_context%icontxt, & 
                                          mlbddc%g_context%np, &
                                          mlbddc%nl_coarse, &
                                          mlbddc%max_coarse_dofs)
       ! end if

       call snd_coarse_mesh_lnods ( mlbddc%g_context%icontxt, & 
                                    mlbddc%g_context%np,      &
                                    mlbddc%pad_collectives,          &  
                                    mlbddc%nl_coarse,                &
                                    mlbddc%coarse_dofs,              &
                                    mlbddc%max_coarse_dofs,          &
                                    mlbddc%p_mat%dof_dist%omap%ng,    &
                                    mlbddc%p_mat%dof_dist%max_nparts, &
                                    mlbddc%p_mat%dof_dist%omap%nl,    &
                                    mlbddc%p_mat%dof_dist%omap%l2g,   &
                                    mlbddc%p_mat%dof_dist%lobjs)      
                                    ! AFM: ng_coarse is not ACTUALLY needed
                                    !      in the fine-grid tasks for the algorithm
                                    !      we are using for assembling coarse-grid graph
                                    ! mlbddc%ng_coarse )

 
    else if( i_am_coarse_task ) then
       ! Both omap%ng and max_nparts are required next in the coarse-grid
       ! tasks. TRANSFER these data!
       call transfer_coarse_int ( mlbddc%g_context%icontxt, &
                                  mlbddc%g_context%iam, &
                                  mlbddc%g_context%np, &
                                  mlbddc%p_mat%dof_dist%omap%ng )

       call transfer_coarse_int ( mlbddc%g_context%icontxt, &
                                  mlbddc%g_context%iam, &
                                  mlbddc%g_context%np, &
                                  mlbddc%p_mat%dof_dist%max_nparts )

       c_mesh%nelem = mlbddc%g_context%np-1
       call memalloc ( c_mesh%nelem+1, c_mesh%pnods,__FILE__,__LINE__ ) 
       call memalloc ( mlbddc%g_context%np+1, mlbddc%ptr_coarse_dofs,__FILE__,__LINE__ ) 

       call rcv_coarse_int_ip ( mlbddc%g_context%icontxt, & 
                                mlbddc%g_context%np, &
                                mlbddc%ptr_coarse_dofs(2) )


       ! AFM: IMPORTANT NOTE !!!! You need to call all_reduce_coarse_mesh_nnode
       ! EVEN with nopad mode, as c_mesh%nnode is later required for the generation 
       ! of the primal graph in create_graph_from_mesh_matrix (mesh_graph.f90)
       ! if( mlbddc%pad_collectives == pad )  then
       call allreduce_coarse_mesh_nnode(  mlbddc%g_context%icontxt, & 
                                          mlbddc%g_context%np, &
                                          mlbddc%nl_coarse, &
                                          mlbddc%max_coarse_dofs)
       ! end if

       c_mesh%nnode = mlbddc%max_coarse_dofs
   
       ! Transform length to header
       mlbddc%ptr_coarse_dofs(1)=0
       c_mesh%pnods(1)=1
       do i = 1, mlbddc%g_context%np-1
          mlbddc%ptr_coarse_dofs(i+1) = mlbddc%ptr_coarse_dofs(i+1) + mlbddc%ptr_coarse_dofs(i)
          c_mesh%pnods(i+1) = mlbddc%ptr_coarse_dofs(i+1) + 1
       end do
       mlbddc%ptr_coarse_dofs(i+1) = mlbddc%ptr_coarse_dofs(i+1) + mlbddc%ptr_coarse_dofs(i)
       
       if ( debug_verbose_level_2 ) then
          write (*,*)  'c_mesh%pnods:', c_mesh%pnods
          call psb_barrier ( mlbddc%c_context%icontxt )
       end if

       call memalloc ( c_mesh%pnods(c_mesh%nelem+1)-1, c_mesh%lnods,__FILE__,__LINE__ ) 

       ! lnods is given in local numbering, the global number of (coarse) nodes being
       ! mlbddc%p_mat%dof_dist%omap%ng which is used to allocate a work array
       ! TODO: use a hash table to avoid it.
       ! If we were interested in the global numbering of nodes, it could be stored
       ! at this point, where it is lost.
       ! AFM: DONE
       call rcv_coarse_mesh_lnods ( mlbddc%g_context%icontxt, & 
                                    mlbddc%g_context%np, &
                                    mlbddc%c_context%np, &
                                    mlbddc%pad_collectives,  &  
                                    mlbddc%max_coarse_dofs,  &
                                    mlbddc%p_mat%dof_dist%omap%ng,  &
                                    mlbddc%ptr_coarse_dofs,   & 
                                    c_mesh%lnods,   &
                                    mlbddc%vars,   &
                                    mlbddc%ng_coarse )

       c_mesh%npoin = mlbddc%ng_coarse

       if ( debug_verbose_level_2 ) then
          write (*,*)  'cmesh%lnods:', c_mesh%lnods
          call psb_barrier ( mlbddc%c_context%icontxt )
       end if

    end if
    
    if ( temp_fine_coarse_grid_overlap ) then
       if ( i_am_fine_task .or. i_am_coarse_task ) then
          if (mlbddc%co_sys_sol_strat == serial_gather) then
             call par_timer_stop ( mlbddc%timer_coll_assprec  )
             call psb_barrier ( mlbddc%g_context%icontxt )
             if ( i_am_fine_task ) then
                call par_timer_start ( mlbddc%timer_assprec_ov_fine )
             else
                call par_timer_start ( mlbddc%timer_assprec_ov_coarse )
             end if
          end if
       end if
    end if

  end subroutine gather_coarse_mesh

  subroutine transfer_coarse_int (icontxt_g, iam_g, np_g, data)
    ! Formerly count_global_coarse_dofs_f_tasks_c_tasks_w_coarse_duties
use psb_const_mod_names
use psb_penv_mod_names
#ifdef MPI_MOD
use mpi
#endif
    implicit none
#ifdef MPI_H
    include 'mpif.h'
#endif

    ! Parameters
    integer    , intent(in)     :: icontxt_g, iam_g, np_g
    integer(ip), intent(inout)  :: data

    ! Locals
    integer     :: mpi_comm_g

    ! Get MPI communicator associated to icontxt (in
    ! the current implementation of our wrappers
    ! to the MPI library icontxt and mpi_comm are actually 
    ! the same)
    call psb_get_mpicomm (icontxt_g, mpi_comm_g)
    
    if (iam_g == 0) then
       call psb_snd(icontxt_g, data, np_g-1)
    else if (iam_g == np_g-1) then
       call psb_rcv(icontxt_g, data, 0)
    end if
    
  end subroutine transfer_coarse_int


  subroutine snd_coarse_int_ip (icontxt_g, np_g, nl_coarse)
    ! Formerly count_global_coarse_dofs_f_tasks_c_tasks_w_coarse_duties
use psb_const_mod_names
use psb_penv_mod_names
#ifdef MPI_MOD
use mpi
#endif
    implicit none
#ifdef MPI_H
    include 'mpif.h'
#endif

    ! Parameters
    integer    , intent(in)  :: icontxt_g, np_g
    integer(ip), intent(in)  :: nl_coarse

    ! Locals
    integer     :: mpi_comm_g, root_g, iret, rcv_dum

    ! Get MPI communicator associated to icontxt (in
    ! the current implementation of our wrappers
    ! to the MPI library icontxt and mpi_comm are actually 
    ! the same)
    call psb_get_mpicomm (icontxt_g, mpi_comm_g)
    root_g = np_g -1 

    call mpi_gather( nl_coarse, 1, psb_mpi_integer, rcv_dum, 1, psb_mpi_integer, root_g, mpi_comm_g, iret)
    check( iret == mpi_success )

  end subroutine snd_coarse_int_ip

  subroutine snd_coarse_int_igp (icontxt_g, np_g, nl_coarse)
    ! Formerly count_global_coarse_dofs_f_tasks_c_tasks_w_coarse_duties
use psb_const_mod_names
use psb_penv_mod_names
#ifdef MPI_MOD
use mpi
#endif
    implicit none
#ifdef MPI_H
    include 'mpif.h'
#endif

    ! Parameters
    integer     , intent(in)  :: icontxt_g, np_g
    integer(igp), intent(in)  :: nl_coarse

    ! Locals
    integer      :: mpi_comm_g, root_g, iret
    integer(igp) :: rcv_dum

    ! Get MPI communicator associated to icontxt (in
    ! the current implementation of our wrappers
    ! to the MPI library icontxt and mpi_comm are actually 
    ! the same)
    call psb_get_mpicomm (icontxt_g, mpi_comm_g)
    root_g = np_g -1 

    call mpi_gather( nl_coarse, 1, psb_mpi_long_integer, rcv_dum, 1, psb_mpi_long_integer, root_g, mpi_comm_g, iret)
    check( iret == mpi_success )

  end subroutine snd_coarse_int_igp

  subroutine rcv_coarse_int_ip ( icontxt_g, np_g, vec_int)
    ! Formerly count_global_coarse_dofs_f_tasks_c_tasks_w_coarse_duties
use psb_const_mod_names
use psb_penv_mod_names
#ifdef MPI_MOD
use mpi
#endif
    implicit none

#ifdef MPI_H
    include 'mpif.h'
#endif

    ! Parameters
    integer    , intent(in)  :: icontxt_g, np_g
    integer(ip), intent(out) :: vec_int(*)

    ! Locals
    integer     :: mpi_comm_g, root_g, int_zero, iret

    integer(ip) :: i

    ! Get MPI communicator associated to icontxt (in
    ! the current implementation of our wrappers
    ! to the MPI library icontxt and mpi_comm are actually 
    ! the same)
    call psb_get_mpicomm (icontxt_g, mpi_comm_g)
    root_g = np_g -1 
    int_zero = 0

    call mpi_gather( int_zero, 1, psb_mpi_integer, vec_int, 1, psb_mpi_integer, root_g, mpi_comm_g, iret)

    check( iret == mpi_success )

  end subroutine rcv_coarse_int_ip

    subroutine rcv_coarse_int_igp ( icontxt_g, np_g, vec_int)
    ! Formerly count_global_coarse_dofs_f_tasks_c_tasks_w_coarse_duties
use psb_const_mod_names
use psb_penv_mod_names
#ifdef MPI_MOD
use mpi
#endif
    implicit none

#ifdef MPI_H
    include 'mpif.h'
#endif

    ! Parameters
    integer     , intent(in)  :: icontxt_g, np_g
    integer(igp), intent(out) :: vec_int(*)

    ! Locals
    integer     :: mpi_comm_g, root_g, iret
    integer(igp):: int_zero

    integer(ip) :: i

    ! Get MPI communicator associated to icontxt (in
    ! the current implementation of our wrappers
    ! to the MPI library icontxt and mpi_comm are actually 
    ! the same)
    call psb_get_mpicomm (icontxt_g, mpi_comm_g)
    root_g = np_g -1 
    int_zero = 0_igp

    
    call mpi_gather( int_zero, 1, psb_mpi_long_integer, vec_int, 1, psb_mpi_long_integer, root_g, mpi_comm_g, iret)
    check( iret == mpi_success )
    
  end subroutine rcv_coarse_int_igp

  subroutine allreduce_coarse_mesh_nnode( icontxt_g, np_g, nl_coarse, max_coarse_dofs)
use psb_const_mod_names
use psb_penv_mod_names
#ifdef MPI_MOD
use mpi
#endif
    implicit none

#ifdef MPI_H
    include 'mpif.h'
#endif

    ! Parameters
    integer    , intent(in)  :: icontxt_g, np_g
    integer(ip), intent(in)  :: nl_coarse
    integer(ip), intent(out) :: max_coarse_dofs

    ! Locals
    integer     :: mpi_comm_g, iret

    ! Get MPI communicator associated to icontxt (in
    ! the current implementation of our wrappers
    ! to the MPI library icontxt and mpi_comm are actually 
    ! the same)
    call psb_get_mpicomm (icontxt_g, mpi_comm_g)
    call mpi_allreduce (nl_coarse, max_coarse_dofs, 1, psb_mpi_integer, mpi_max, mpi_comm_g, iret)
    check( iret == mpi_success )

  end subroutine allreduce_coarse_mesh_nnode

  subroutine snd_coarse_mesh_lnods ( icontxt_g, np_g, pad_collectives, nl_coarse, coarse_dofs, &
                                     max_coarse_dofs, gnobjs, max_nparts, lnobjs, l2g, lobjs )!, &  
                                     ! ng_coarse)
    ! Formerly list_global_coarse_dofs_f_tasks_c_tasks_w_coarse_duties
use psb_const_mod_names
use psb_penv_mod_names

#ifdef MPI_MOD
use mpi
#endif

    implicit none

#ifdef MPI_H
    include 'mpif.h'
#endif

    ! Parameters
    integer    , intent(in)  :: icontxt_g, np_g
    integer(ip), intent(in)  :: pad_collectives, nl_coarse
    integer(ip), intent(in)  :: coarse_dofs (nl_coarse), max_coarse_dofs
    integer(ip), intent(in)  :: gnobjs, max_nparts, lnobjs, l2g(lnobjs), lobjs(max_nparts+4,lnobjs)
    ! integer(ip), intent(out) :: ng_coarse

    ! Locals
    integer      :: mpi_comm_g, root_g, iret, rcv_dum, rcv_dumvec(1)
   
    ! This declaration is just to avoid compiler
    ! errors in the SERIAL-ALL-GATHER itinerary
    ! integer      :: np, mpi_comm, icontxt
    integer(ip)  :: i, j, sz, off 

    integer(ip), allocatable :: coarse_dofs_vars_global (:,:)
    integer(ip), allocatable :: work1 (:)
    integer(ip), allocatable :: work2 (:)
    integer(ip), allocatable :: work3 (:,:) 
    integer(ip), allocatable :: pad_buf_snd (:,:), pad_buf_rcv(:,:)
    integer(ip)              :: ng_corners_on_edges_faces

    ! Get MPI communicator associated to icontxt (in
    ! the current implementation of our wrappers
    ! to the MPI library icontxt and mpi_comm are actually 
    ! the same)
    call psb_get_mpicomm (icontxt_g, mpi_comm_g)
    root_g = np_g -1

    call memalloc ( 2, nl_coarse, coarse_dofs_vars_global, __FILE__,__LINE__ )
    do i=1, nl_coarse
       coarse_dofs_vars_global(1,i) = l2g(coarse_dofs(i))
       coarse_dofs_vars_global(2,i) = lobjs(1,coarse_dofs(i))
       !write(*,*) 'YYY', me_p, l2g(coarse_dofs(i)), lobjs(1,coarse_dofs(i))
    end do

    if  ( pad_collectives == pad ) then
       ! Pack lpadj_self into send buffer
       call memalloc ( 2, max_coarse_dofs, pad_buf_snd, __FILE__,__LINE__ )
       pad_buf_snd (1, 1:nl_coarse) = coarse_dofs_vars_global(1,:)
       pad_buf_snd (2, 1:nl_coarse) = coarse_dofs_vars_global(2,:)

       call mpi_gather( pad_buf_snd, 2*max_coarse_dofs, psb_mpi_integer,  &
                        rcv_dum, rcv_dum, psb_mpi_integer, root_g, mpi_comm_g, iret)
       check ( iret == mpi_success )

       call memfree ( pad_buf_snd,__FILE__,__LINE__)

    else if (pad_collectives == nopad ) then 
       call mpi_gatherv( coarse_dofs_vars_global, nl_coarse*2, psb_mpi_integer, &
                         rcv_dum, rcv_dumvec, rcv_dumvec, psb_mpi_integer, &
                         root_g, mpi_comm_g, iret)
       check ( iret == mpi_success )
    end if

    call memfree ( coarse_dofs_vars_global,__FILE__,__LINE__)
    ! call psb_bcast (icontxt_g, ng_coarse, root=root_g)

  end subroutine snd_coarse_mesh_lnods

  subroutine rcv_coarse_mesh_lnods ( icontxt_g, np_g, np_c, pad_collectives, &
                                     max_coarse_dofs, gnobjs, ptr_coarse_dofs,  lst_coarse_dofs,  &  
                                     vars, ng_coarse)
    ! Formerly list_global_coarse_dofs_f_tasks_c_tasks_w_coarse_duties
    ! TODO: avoid allocating work2(gnobjs) using a hash_table
    ! AFM: DONE
use psb_const_mod_names
use psb_penv_mod_names

#ifdef MPI_MOD
use mpi
#endif

    implicit none

#ifdef MPI_H
    include 'mpif.h'
#endif

    ! Parameters
    integer    , intent(in)               :: icontxt_g, np_g, np_c
    integer(ip), intent(in)               :: pad_collectives
    integer(ip), intent(in)               :: max_coarse_dofs
    integer(ip), intent(in)               :: gnobjs
    integer(ip), intent(in)               :: ptr_coarse_dofs(np_g+1)
    integer(ip), intent(out)              :: lst_coarse_dofs(ptr_coarse_dofs(np_g+1))
    integer(ip), allocatable, intent(out) :: vars(:)
    integer(ip), intent(out) :: ng_coarse

    ! Locals
    integer      :: mpi_comm_g, root_g, iret, snd_dum
    integer(ip)  :: aux_val 
   
    ! This declaration is just to avoid compiler
    ! errors in the SERIAL-ALL-GATHER itinerary
    ! integer      :: np, mpi_comm, icontxt
    integer(ip)  :: i, j, sz, off, istat 

    integer(ip), allocatable :: recvcounts (:)
    integer(ip), allocatable :: work1 (:)
    integer(ip), allocatable :: work2 (:)
    integer(ip), allocatable :: work3 (:,:) 
    integer(ip), allocatable :: pad_buf_snd (:,:), pad_buf_rcv(:,:)
    type(hash_table_ip_ip_t)   :: objs_visited, objs_positions 

    ! Set size of the hash table as 10% of the average number of coarse
    ! DoFs per subdomain    
    call objs_visited%init(max(int(real(gnobjs,rp)/real(np_c,rp)*0.1_rp,ip),5))
    call objs_positions%init(max(int(real(gnobjs,rp)/real(np_c,rp)*0.1_rp,ip),5))

    ! Get MPI communicator associated to icontxt (in
    ! the current implementation of our wrappers
    ! to the MPI library icontxt and mpi_comm are actually 
    ! the same)
    call psb_get_mpicomm (icontxt_g, mpi_comm_g)
    root_g = np_g -1

    if  ( pad_collectives == pad ) then
       ! Pack lpadj_self into send buffer
       call memalloc ( 2, max_coarse_dofs, pad_buf_snd, __FILE__,__LINE__ )

       call memalloc ( 2, max_coarse_dofs*np_g, pad_buf_rcv, __FILE__,__LINE__ )

       call mpi_gather ( pad_buf_snd, 2*max_coarse_dofs, psb_mpi_integer, &
                         pad_buf_rcv, 2*max_coarse_dofs, psb_mpi_integer, &
                         root_g, mpi_comm_g, iret)

       check( iret == mpi_success )

       call memalloc ( ptr_coarse_dofs(np_g+1), work1,  __FILE__,__LINE__ )

       ! Unpack padded data
       off = 1
       do i=1, np_g
          sz = ptr_coarse_dofs(i+1) - ptr_coarse_dofs(i)
          work1 ( ptr_coarse_dofs(i)+1:ptr_coarse_dofs(i+1) ) = pad_buf_rcv ( 1, off:off+sz-1 ) 
          ! lst_vars ( ptr_coarse_dofs(i)+1:ptr_coarse_dofs(i+1) ) = pad_buf_rcv ( 2, off:off+sz-1 ) 
          off = off + max_coarse_dofs
       end do
       
       ! write (*,*) 'XXX', ptr_coarse_dofs(1:np_g+1)
       ! write (*,*) 'YYY', work1 (:,1:ptr_coarse_dofs(np+1)-1)
          
       ! Free padded data
       ! call memalloc ( gnobjs, work2, __FILE__,__LINE__ ) 
       ! work2                     = 0 
       ng_coarse                 = 0
       off = 1
       do i=1, np_g 
         do j=ptr_coarse_dofs(i)+1, ptr_coarse_dofs(i+1) 
           aux_val = ng_coarse+1
           call objs_visited%put(key=work1(j),val=aux_val, stat=istat) 
           if ( istat == now_stored ) then
               aux_val = off+j-ptr_coarse_dofs(i)-1
               call objs_positions%put(key=ng_coarse+1,val=aux_val, stat=istat) 
               assert ( istat == now_stored )
               ng_coarse = ng_coarse + 1
               ! work2(work1(i+1))  = ng_coarse
               lst_coarse_dofs(j) = ng_coarse
           else
               call objs_visited%get(key=work1(j),val=lst_coarse_dofs(j), stat=istat) 
               ! lst_coarse_dofs(i+1) = work2(work1(i+1))
           end if
         end do
         off = off + max_coarse_dofs
      end do
      call memfree ( work1,__FILE__,__LINE__)
      call objs_visited%free 

      call memalloc ( ng_coarse, vars, __FILE__, __LINE__)
      do i=1, ng_coarse
         call objs_positions%get(key=i,val=j, stat=istat)
         assert ( istat == key_found )
         vars(i) = pad_buf_rcv ( 2, j ) 
      end do

      call objs_positions%free 
      ! call memfree ( work2,__FILE__,__LINE__)
      call memfree ( pad_buf_snd,__FILE__,__LINE__)
      call memfree ( pad_buf_rcv,__FILE__,__LINE__)

    else if (pad_collectives == nopad ) then

       call memalloc ( np_g, recvcounts,  __FILE__,__LINE__ )

       do i=1, np_g
          recvcounts(i) = 2*(ptr_coarse_dofs(i+1)-ptr_coarse_dofs(i))
       end do

       call memalloc ( 2, ptr_coarse_dofs(np_g+1), work3,  __FILE__,__LINE__ )

       call memalloc ( np_g+1, work2, __FILE__,__LINE__ )

       work2(1) = 0
       do i=1, np_g
          work2(i+1) = work2(i) + recvcounts(i)
       end do

       call mpi_gatherv( snd_dum, 0, psb_mpi_integer, &
                         work3, recvcounts,  work2, psb_mpi_integer, &
                         root_g, mpi_comm_g, iret)

       check ( iret == mpi_success )

       call memfree  ( recvcounts,__FILE__,__LINE__)
       call memfree ( work2,__FILE__,__LINE__)

       ! call memalloc ( gnobjs, work2, __FILE__,__LINE__ ) 
       ! work2                     = 0 
       ng_coarse                 = 0
       do i=0, ptr_coarse_dofs(np_g+1)-1
          ! lst_vars(i+1) = work3(2,i+1)
          aux_val = ng_coarse + 1
          call objs_visited%put(key=work3(1,i+1),val=aux_val, stat=istat) 
          if ( istat == now_stored ) then
             aux_val = i + 1 
             call objs_positions%put(key=ng_coarse+1,val=aux_val, stat=istat) 

             ng_coarse = ng_coarse + 1
             ! work2(work3(1,i+1))  = ng_coarse
             lst_coarse_dofs(i+1) = ng_coarse
          else
             ! lst_coarse_dofs(i+1) = work2(work3(1,i+1))
             call objs_visited%get(key=work3(1,i+1),val=lst_coarse_dofs(i+1), stat=istat)
          end if
       end do
       ! call memfree ( work2,__FILE__,__LINE__) 
       call objs_visited%free 

       call memalloc ( ng_coarse, vars, __FILE__, __LINE__)
       do i=1, ng_coarse
          call objs_positions%get(key=i,val=j, stat=istat)
          assert ( istat == key_found )
          vars(i) = work3 ( 2, j ) 
       end do
       call memfree ( work3,__FILE__,__LINE__)
       call objs_positions%free 
    end if

    ! call psb_bcast (icontxt_g, ng_coarse, root=root_g)

  end subroutine rcv_coarse_mesh_lnods

  subroutine generate_coarse_graph (mlbddc)
use mpi
    implicit none
    ! Parameters 
    ! On input , mlbddc%f_mesh_c%pnods in C-based indexing
    ! On output, mlbddc%f_mesh_c%pnods in F-based indexing 
    type(par_preconditioner_dd_mlevel_bddc_t), intent(inout) :: mlbddc

    ! Locals
    integer(ip)    :: gtype

    ! AFM: Derived from the fact that f_graph within p_graph_c was not ACTUALLY created
    if (mlbddc%A_c%symm == symm_true) then
       gtype = csr_symm
    else
       gtype = csr
    end if

    call mesh_to_graph_matrix ( gtype, mlbddc%f_mesh_c, mlbddc%A_c_gr )

    if ( debug_verbose_level_2 ) then 
       call graph_print ( 6, mlbddc%A_c_gr )
    end if

  end subroutine generate_coarse_graph

  subroutine generate_coarse_partition_and_graph (c_mesh, mlbddc)
use mpi
    implicit none
    ! Parameters 
    type(mesh_t), intent(inout) :: c_mesh
    type(par_preconditioner_dd_mlevel_bddc_t), intent(inout) :: mlbddc

    ! Locals
    integer(ip)              :: i, sum_nnz, max_sum_nnz
    integer(ip), allocatable :: local_lobjs_coarse (:)
    integer(ip), allocatable :: nl_coarse_neigh(:)   ! Number of local coarse dofs of my neighbours
    logical                  :: i_am_coarse_task, i_am_fine_task, i_am_higher_level_task

    type(mesh_t)           :: dual_f_mesh
    integer(ip), allocatable :: neighbours_coarse_parts (:)
    integer(ip), allocatable :: dual_parts (:)
    integer(ip), allocatable :: l2ge(:)
    ! integer(ip), allocatable :: vars (:)

    integer(ip)               ::  nebou       ! number of boundary elements
    integer(ip), allocatable  ::  pextn(:)    ! Pointers to the lextn
    integer(igp), allocatable ::  lextn(:)    ! List of (GID of) external neighbors
    integer(ip), allocatable  ::  lextp(:)    ! List of parts of external neighbors
    integer(ip), allocatable  ::  lexte(:)    ! Edge information of external neighbors
    integer(ip), allocatable  ::  lextm(:)    ! Edge information of external neighbors
    integer(ip)  :: gtype

    ! Which duties I have?
    i_am_fine_task = (mlbddc%p_mat%p_env%p_context%iam >= 0)
    i_am_coarse_task = (mlbddc%c_context%iam >= 0)
    i_am_higher_level_task = (mlbddc%d_context%iam >= 0)


    if(mlbddc%co_sys_sol_strat == recursive_bddc) then

       if ( i_am_coarse_task ) then
          ! TODO: all rcv_ rounties are wrong because ptr_coarse_dofs starts from 0.
          ! DONE: we now keep mlbddc%ptr_coarse_dofs  & c_mesh%lnods separately
          call rcv_subdomains_surrounding_coarse_dofs ( mlbddc%g_context%icontxt, &
                                                        mlbddc%g_context%iam,     &
                                                        mlbddc%g_context%np,      &
                                                        mlbddc%pad_collectives,          &
                                                        mlbddc%p_mat%dof_dist%max_nparts, &  
                                                        c_mesh%npoin,    &
                                                        c_mesh%nnode,    &
                                                        mlbddc%ptr_coarse_dofs,    &
                                                        c_mesh%lnods,    &
                                                        dual_f_mesh, &
                                                        dual_parts )

          if ( debug_verbose_level_2 ) then
             write (*,*)  'dual_f_mesh%pnods:', dual_f_mesh%pnods
             write (*,*)  'dual_f_mesh%lnods:', dual_f_mesh%lnods
             write (*,*)  'dual_parts:', dual_parts
             call psb_barrier ( mlbddc%c_context%icontxt )
          end if


          !call memalloc (mlbddc%g_context%np-1, l2ge, __FILE__,__LINE__) 
          call memalloc (c_mesh%nelem, l2ge, __FILE__,__LINE__) 
          
          call rcv_coarse_int_ip(mlbddc%g_context%icontxt, &
                                 mlbddc%g_context%np,      &
                                 l2ge)

          if ( debug_verbose_level_2 ) then
             write (*,*)  'l2ge:', l2ge
             call psb_barrier ( mlbddc%c_context%icontxt )
          end if

          if ( temp_fine_coarse_grid_overlap ) then
             call par_timer_stop ( mlbddc%timer_coll_assprec  )
             call psb_barrier ( mlbddc%g_context%icontxt )
             call par_timer_start ( mlbddc%timer_assprec_ov_coarse )
             call par_timer_start ( mlbddc%timer_assprec_ov_coarse_header )
          end if

!!$          call build_partition_adjacency (mlbddc%c_context%iam+1, &
!!$               &                          c_mesh,& 
!!$               &                          l2ge, &
!!$               &                          dual_f_mesh, & 
!!$               &                          dual_parts, &
!!$               &                          nebou, &
!!$               &                          lebou, &
!!$               &                          pextn, &
!!$               &                          lextn, &
!!$               &                          lextp, &
!!$               &                          lexte)
!!$          if ( debug_verbose_level_2 ) then 
!!$             write (*,*)  'nebou:', nebou
!!$             write (*,*)  'pextn:', pextn
!!$             write (*,*)  'lextn:', lextn
!!$             write (*,*)  'lextp:', lextp
!!$             write (*,*)  'lexte:', lexte
!!$             call psb_barrier ( mlbddc%c_context%icontxt )
!!$          end if



          ! AFM: The extraction of vars(:) could be integrated into rcv_coarse_mesh_lnods.
          ! This would be more efficient than extracting it here. However, this would require to 
          ! have vars(:) as a member of par_preconditioner_dd_mlevel_bddc. As far as I know, this is 
          ! temporal data only required here. Right? 
          ! call memalloc (c_mesh%npoin, vars, __FILE__,__LINE__) 
          ! DONE
          ! do i=1, c_mesh%pnods(c_mesh%nelem+1)-1
          !   vars(c_mesh%lnods(i)) = mlbddc%lst_vars(i)
          ! end do
          
          if ( debug_verbose_level_2 ) then 
             write (*,*)  'vars:', mlbddc%vars
             call psb_barrier ( mlbddc%c_context%icontxt )
          end if


          call dof_distribution_coarse_create ( mlbddc%c_context%icontxt, & ! Communication context
                                                mlbddc%c_context%iam    , &
                                                mlbddc%c_context%np     , &
                                                c_mesh                        , & ! Local mesh in an initial ARBITRARY local node/element ordering
                                                l2ge                          , & ! local2global element correspondence
                                                dual_f_mesh                   , & ! Associated dual_mesh with external elements also listed for boundary DoFs
                                                dual_parts                    , & ! Parts associated to each element in the dual mesh
                                                mlbddc%vars                   , & ! Physical unknown corresponding to each DoF in f_mesh
                                                mlbddc%p_mesh_c%f_mesh        , & ! New mesh object
                                                mlbddc%dof_dist_c             , & ! New distributed DoF object conformal to new local mesh
                                                mlbddc%erenumbering_c                  )   ! Resulting Element renumbering


          ! Free dual_f_mesh with (external) adjacency data built-in
          call mesh_free ( dual_f_mesh )

          if ( debug_verbose_level_2 ) then 
             write (*,*)  'mlbddc%p_mesh_c%f_mesh:', mlbddc%p_mesh_c%f_mesh%pnods
             write (*,*)  'mlbddc%p_mesh_c%f_mesh:', mlbddc%p_mesh_c%f_mesh%lnods
             write (*,*)  'erenumbering_c%lperm:', mlbddc%erenumbering_c%lperm 
             write (*,*)  'erenumbering_c%iperm:', mlbddc%erenumbering_c%iperm 
             write (*,*)  'mlbddc%p_mesh_c%f_mesh:', mlbddc%p_mesh_c%f_mesh%lnods
             call dof_distribution_print ( 6, mlbddc%dof_dist_c )
             call psb_barrier ( mlbddc%c_context%icontxt )
          end if

          ! AFM: Derived from the fact that f_graph within p_graph_c was not ACTUALLY created
          if (mlbddc%p_mat_c%f_matrix%symm == symm_true) then
             gtype = csr_symm
          else
             gtype = csr
          end if
          call mesh_to_graph_matrix ( gtype, mlbddc%p_mesh_c%f_mesh, mlbddc%p_graph_c%f_graph)

          if ( debug_verbose_level_2 ) then 
             call graph_print ( 6,  mlbddc%p_graph_c%f_graph )
             call psb_barrier ( mlbddc%c_context%icontxt )
          end if

          if ( temp_fine_coarse_grid_overlap ) then
             call par_timer_stop ( mlbddc%timer_assprec_ov_coarse_header )
          end if

          call memfree (l2ge, __FILE__,__LINE__) 
          ! call memfree (vars, __FILE__,__LINE__) 
          call memfree (dual_parts, __FILE__,__LINE__)
          call mesh_free ( c_mesh )  

       else if ( i_am_fine_task ) then

          call memalloc ( mlbddc%p_mat%dof_dist%npadj, neighbours_coarse_parts, __FILE__,__LINE__) 

          call fetch_parts_coarse_neighbours ( mlbddc%p_mat%p_env%p_context%icontxt, & 
                                               mlbddc%p_mat%dof_dist%npadj, & 
                                               mlbddc%p_mat%dof_dist%lpadj, & 
                                               mlbddc%p_mat%p_env%id_parts(2), &
                                               neighbours_coarse_parts)

          call snd_subdomains_surrounding_coarse_dofs ( mlbddc%g_context%icontxt, &
                                                        mlbddc%g_context%iam,     &
                                                        mlbddc%g_context%np,      &
                                                        mlbddc%pad_collectives,          &
                                                        mlbddc%p_mat%dof_dist%max_nparts, &  
                                                        mlbddc%nl_coarse,                &
                                                        mlbddc%coarse_dofs,              &
                                                        mlbddc%p_mat%dof_dist%omap%nl,    &
                                                        mlbddc%p_mat%dof_dist%lobjs,      &
                                                        mlbddc%p_mat%dof_dist%npadj,      &
                                                        mlbddc%p_mat%dof_dist%lpadj,      &
                                                        mlbddc%p_mat%p_env%id_parts,          &
                                                        neighbours_coarse_parts,            &
                                                        mlbddc%max_coarse_dofs)

          call memfree (neighbours_coarse_parts, __FILE__,__LINE__) 

          ! TODO my_part = iam + 1 everywhere, right? just check it...
          call snd_coarse_int_ip( mlbddc%g_context%icontxt, &
                                   mlbddc%g_context%np,      &
                                   mlbddc%p_mat%p_env%p_context%iam+1 ) 
                                   ! snd_igp ) 

          if ( temp_fine_coarse_grid_overlap ) then
             call par_timer_stop ( mlbddc%timer_coll_assprec  )
             call psb_barrier ( mlbddc%g_context%icontxt )
             call par_timer_start ( mlbddc%timer_assprec_ov_fine )
          end if

       end if

    end if

  end subroutine generate_coarse_partition_and_graph

  subroutine fetch_parts_coarse_neighbours ( icontxt, & 
                                            npadj, & 
                                            lpadj, & 
                                            my_part_coarse, &
                                            parts_coarse_neigh)
#ifdef MPI_MOD
use mpi
#endif 
    implicit none

#ifdef MPI_H
    include 'mpif.h'
#endif

    ! Parameters
    integer(ip), intent(in)  :: icontxt, npadj, my_part_coarse
    integer(ip), intent(in)  :: lpadj(npadj)
    integer(ip), intent(out) :: parts_coarse_neigh(npadj)

    ! Locals
    integer     :: proc_to_comm, sizmsg
    integer     :: mpi_comm,  iret
    integer     :: p2pstat(mpi_status_size)
    integer(ip) :: iposadj

    ! Request handlers for non-blocking receives
    integer, allocatable, dimension(:) :: rcvhd

    ! Request handlers for non-blocking receives
    integer, allocatable, dimension(:) :: sndhd
    
    ! Get MPI communicator associated to icontxt (in
    ! the current implementation of our wrappers
    ! to the MPI library icontxt and mpi_comm are actually 
    ! the same)
    call psb_get_mpicomm (icontxt, mpi_comm)

    call memalloc ( npadj, rcvhd, __FILE__,__LINE__ )
    call memalloc ( npadj, sndhd, __FILE__,__LINE__ )

    ! Determine extra_corners_rcv from extra_corners_snd by means of a nearest neighbours exchange

    ! First post all the non-blocking receives   
    do iposadj=1, npadj

       proc_to_comm = lpadj (iposadj)

       ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
       call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)

       ! Message size to be received
       sizmsg = 1

       call mpi_irecv(  parts_coarse_neigh(iposadj), sizmsg, &
            &           psb_mpi_integer, proc_to_comm, &
            &           psb_int_swap_tag, mpi_comm, rcvhd(iposadj), iret)

       check ( iret == mpi_success )

    end do

    ! Second post all the non-blocking sends   
    do iposadj=1, npadj
       proc_to_comm = lpadj (iposadj)

       ! Get MPI rank id associated to proc_to_comm-1 in proc_to_comm
       ! (In our MPI wrapper, this is not actually required. It is
       !  a formal use of the PSB MPI wrappers.).
       call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)

       ! Message size to be sent
       sizmsg = 1

       call mpi_isend(  my_part_coarse, sizmsg,        &
            &           psb_mpi_integer, proc_to_comm, &
            &           psb_int_swap_tag, mpi_comm, sndhd(iposadj), iret)

       check ( iret == mpi_success )

    end do

    ! Then, wait on all non-blocking receives
    do iposadj=1, npadj
       proc_to_comm = lpadj (iposadj)
       call mpi_wait(rcvhd(iposadj), p2pstat, iret)
       check ( iret == mpi_success )
    end do
    
    ! Finally, wait on all non-blocking sends
    do iposadj=1, npadj
       call mpi_wait(sndhd(iposadj), p2pstat, iret)
       check ( iret == mpi_success )
    end do

    call memfree ( rcvhd,__FILE__,__LINE__)
    call memfree ( sndhd,__FILE__,__LINE__)

  end subroutine fetch_parts_coarse_neighbours

  subroutine snd_subdomains_surrounding_coarse_dofs ( icontxt_g, me_g, np_g, pad_collectives, &
                                                              max_nparts, nl_coarse, coarse_dofs, &  
                                                              lnobjs, lobjs, npadj, lpadj, id_parts, &
                                                              parts_coarse_neigh, max_coarse_dofs)

use psb_const_mod_names
use psb_penv_mod_names

#ifdef MPI_MOD
use mpi
#endif
    implicit none
#ifdef MPI_H
    include 'mpif.h'
#endif

    ! Parameters
    integer    , intent(in)  :: icontxt_g, me_g, np_g
    integer(ip), intent(in)  :: pad_collectives, max_nparts, nl_coarse
    integer(ip), intent(in)  :: coarse_dofs (nl_coarse)
    integer(ip), intent(in)  :: lnobjs
    integer(ip), intent(in)  :: lobjs(max_nparts+4,lnobjs)
    integer(ip), intent(in)  :: npadj
    integer(ip), intent(in)  :: lpadj(npadj)
    integer(ip), intent(in)  :: id_parts(2)
    integer(ip), intent(in)  :: parts_coarse_neigh(npadj)
    integer(ip), intent(in)  :: max_coarse_dofs

    ! Locals    
    integer     :: mpi_comm_g,  iret, root_g, rcv_dum, rcv_dumvec(1)
    integer(ip), allocatable :: local_lobjs_coarse (:,:)
    integer(ip)              :: i, j, k, offl, offg, ipart, ipadj

    ! Get MPI communicator associated to icontxt (in
    ! the current implementation of our wrappers
    ! to the MPI library icontxt and mpi_comm are actually 
    ! the same)
    call psb_get_mpicomm (icontxt_g, mpi_comm_g)
    root_g = np_g - 1

    if ( pad_collectives == nopad ) then 
       call memalloc ( 2*max_nparts+1,       nl_coarse,  local_lobjs_coarse,  __FILE__,__LINE__)  
    else if (pad_collectives == pad) then
       call memalloc ( 2*max_nparts+1, max_coarse_dofs,  local_lobjs_coarse,  __FILE__,__LINE__)  
    end if

    do i=1, nl_coarse
       ! Copy object from the local list of objects to the list of coarse objs
       local_lobjs_coarse(             1,i) = lobjs(4, coarse_dofs(i))              ! Number of subdomains around coarse_dofs(i) 
       local_lobjs_coarse(2:max_nparts+1,i) = lobjs(5:max_nparts+4, coarse_dofs(i)) ! List of subdomains around coarse_dofs(i) 
       do ipart = 1, lobjs(4, coarse_dofs(i)) 
          if( lobjs(4+ipart, coarse_dofs(i)) == id_parts(1) ) then
             local_lobjs_coarse(max_nparts+1+ipart,i) = id_parts(2)
             cycle
          end if
          do ipadj = 1,npadj
             if( lobjs(4+ipart, coarse_dofs(i)) == lpadj(ipadj) ) exit
          end do
          local_lobjs_coarse(max_nparts+1+ipart,i) = parts_coarse_neigh(ipadj)
       end do
    end do

    if ( pad_collectives == nopad ) then
       call mpi_gatherv( local_lobjs_coarse,  nl_coarse*(2*max_nparts+1), psb_mpi_integer, &
                         rcv_dum, rcv_dumvec, rcv_dumvec, psb_mpi_integer, &
                         root_g, mpi_comm_g, iret)
       check ( iret == mpi_success )
    else if (pad_collectives == pad) then
       call mpi_gather ( local_lobjs_coarse,  max_coarse_dofs*(2*max_nparts+1), psb_mpi_integer, &
                         rcv_dum, rcv_dum,  psb_mpi_integer, &
                         root_g, mpi_comm_g,  iret)
       check ( iret == mpi_success )
    end if
   
    call memfree ( local_lobjs_coarse,__FILE__,__LINE__) 
   
  end subroutine snd_subdomains_surrounding_coarse_dofs

  subroutine rcv_subdomains_surrounding_coarse_dofs ( icontxt_g, me_g, np_g, & 
                                                      pad_collectives, &
                                                      max_nparts, ng_coarse, &  
                                                      max_coarse_dofs, &
                                                      ptr_coarse_dofs, &
                                                      lst_coarse_dofs, & 
                                                      dual_f_mesh, &
                                                      dual_f_parts )

use psb_const_mod_names
use psb_penv_mod_names

#ifdef MPI_MOD
use mpi
#endif
    implicit none
#ifdef MPI_H
    include 'mpif.h'
#endif

    ! Parameters
    integer    , intent(in)  :: icontxt_g, me_g, np_g

    integer(ip), intent(in)  :: pad_collectives, max_nparts, ng_coarse
    integer(ip), intent(in)  :: max_coarse_dofs
    integer(ip), intent(in)  :: ptr_coarse_dofs(np_g+1)
    integer(ip), intent(in)  :: lst_coarse_dofs(ptr_coarse_dofs(np_g+1))
    type(mesh_t) , intent(out) :: dual_f_mesh
    integer(ip)    , allocatable, intent(out) :: dual_f_parts(:)

    ! Locals    
    integer     :: mpi_comm_g,  iret, root_g, snd_dum

    integer(ip), allocatable :: local_lobjs_coarse  (:,:), work(:,:), work2(:)
    integer(ip), allocatable :: pad_buf_rcv(:,:)

    integer(ip)              :: i, j, k, offl, offg
    integer(ip), allocatable :: recvcounts (:), displs (:)
    integer(ip), allocatable :: global_lobjs_coarse (:,:)

    ! Get MPI communicator associated to icontxt (in
    ! the current implementation of our wrappers
    ! to the MPI library icontxt and mpi_comm are actually
    ! the same)
    call psb_get_mpicomm (icontxt_g, mpi_comm_g)
    root_g = np_g - 1

    ! global_lbojs_coarse contains a list of subdomains and coarse parts around coarse dofs
    ! i.e. a list of elements and parts (global numbering) around points (local numbering)
    call memalloc ( 2*max_nparts + 1, ng_coarse, global_lobjs_coarse, __FILE__,__LINE__) 

    if ( pad_collectives == nopad ) then

       call memalloc ( np_g, recvcounts,   __FILE__,__LINE__ )
       call memalloc ( np_g+1, displs,   __FILE__,__LINE__ )

       displs(1) = 0 
       do i=1, np_g
          recvcounts(i) = (ptr_coarse_dofs(i+1)-ptr_coarse_dofs(i))*(2*max_nparts+1)
          displs(i+1)   = displs(i) + recvcounts(i)
       end do

       call memalloc ( 2*max_nparts+1, displs(np_g+1), work, __FILE__,__LINE__)

       call mpi_gatherv( snd_dum, 0 ,  psb_mpi_integer, &
                         work ,  recvcounts,  displs, psb_mpi_integer, &
                         root_g, mpi_comm_g, iret)
       check ( iret == mpi_success )

       call memfree ( recvcounts,__FILE__,__LINE__)
       call memfree ( displs    ,__FILE__,__LINE__)

       call memalloc ( ng_coarse, work2, __FILE__,__LINE__ ) 
       work2 = 0 
       k     = 0
       do i=0, ptr_coarse_dofs(np_g+1)-1
          if ( work2(lst_coarse_dofs(i+1)) == 0 ) then
             k = k + 1
             work2(lst_coarse_dofs(i+1)) = 1
             global_lobjs_coarse (:,lst_coarse_dofs(i+1)) = work(:,i+1)
             if (k == ng_coarse) exit 
          end if
       end do
       call memfree ( work,__FILE__,__LINE__) 
       call memfree ( work2,__FILE__,__LINE__) 

    else if (pad_collectives == pad) then

       call memalloc ( 2*max_nparts+1, max_coarse_dofs, local_lobjs_coarse,  __FILE__,__LINE__)  
       call memalloc ( 2*max_nparts+1, np_g*max_coarse_dofs*(2*max_nparts+1), pad_buf_rcv, __FILE__,__LINE__)

       call mpi_gather ( local_lobjs_coarse,  max_coarse_dofs*(2*max_nparts+1), psb_mpi_integer, &
                          pad_buf_rcv,  max_coarse_dofs*(2*max_nparts+1), psb_mpi_integer, &
                          root_g,                      mpi_comm_g,  iret)
       check ( iret == mpi_success )

       ! Unpad work          
       call memalloc ( ng_coarse, work2, __FILE__,__LINE__ ) 
       work2 = 0 
       k     = 0
       do j=1, np_g
          offg = (j-1)*max_coarse_dofs
          offl = 1
          do i=ptr_coarse_dofs(j),ptr_coarse_dofs(j+1)-1
             if ( work2(lst_coarse_dofs(i+1)) == 0 ) then
                k = k + 1
                work2(lst_coarse_dofs(i+1)) = 1
                global_lobjs_coarse (:,lst_coarse_dofs(i+1)) = pad_buf_rcv(:,offg+offl)
                if (k == ng_coarse) exit 
             end if
             offl=offl+1
          end do
       end do
       
       call memfree ( pad_buf_rcv,__FILE__,__LINE__) 
       call memfree ( work2,__FILE__,__LINE__) 

       call memfree ( local_lobjs_coarse,__FILE__,__LINE__) 
   end if

   ! AFM: Here ONLY the members of dual_f_mesh that you know that are
   !      used later by build_partition_adjancency are filled up.
   !      This is DIRTY, as this data structure may not work with other
   !      subroutines that operate on mesh, or even with build_partition_adjacency
   !      if its implementation changes in the future. This is because we are
   !      breaking here (and in many parts of the code I have to admit) the
   !      data encapsulation principle. There should be one common subroutine or interface
   !      for the creation of meshes that is used everywhere a mesh is created and that
   !      is responsible for the proper initialization of ALL of its members, and not
   !      a subset of them (as it is done HERE). I will init at least npoin, nelty here
   !      as a temporal solution (what about ndime, nnode, nboun, etc.?)
   dual_f_mesh%npoin = np_g - 1
   dual_f_mesh%nelem = ng_coarse
   call memalloc (dual_f_mesh%nelem+1, dual_f_mesh%pnods, __FILE__,__LINE__) 
   dual_f_mesh%pnods(1) = 1
   do i=1,dual_f_mesh%nelem
      dual_f_mesh%pnods(i+1) = dual_f_mesh%pnods(i) + global_lobjs_coarse(1,i)
   end do
   call memalloc (dual_f_mesh%pnods(dual_f_mesh%nelem+1)-1,dual_f_mesh%lnods, __FILE__,__LINE__) 
   call memalloc (dual_f_mesh%pnods(dual_f_mesh%nelem+1)-1,dual_f_parts, __FILE__,__LINE__) 
   do i=1,dual_f_mesh%nelem
      k=1
      do j=dual_f_mesh%pnods(i), dual_f_mesh%pnods(i+1)-1 
         k = k + 1
         dual_f_mesh%lnods(j) = global_lobjs_coarse(k,i)
         dual_f_parts(j) = global_lobjs_coarse(max_nparts+k,i)
      end do
   end do

   call memfree ( global_lobjs_coarse, __FILE__,__LINE__) 

 end subroutine rcv_subdomains_surrounding_coarse_dofs

  !=================================================================================================
  recursive subroutine par_preconditioner_dd_mlevel_bddc_fill_val ( mlbddc )
                                           
    implicit none
    ! Parameters 
    type(par_preconditioner_dd_mlevel_bddc_t), target, intent(inout) :: mlbddc

    ! Locals
    type(par_matrix_t), pointer :: p_mat
    logical :: i_am_coarse_task, i_am_fine_task, i_am_higher_level_task

    p_mat => mlbddc%p_mat

    ! Which duties do I have?
    i_am_fine_task = (mlbddc%p_mat%p_env%p_context%iam >= 0)
    i_am_coarse_task = (mlbddc%c_context%iam >= 0)
    i_am_higher_level_task = (mlbddc%d_context%iam >= 0)

    ! "Numerical factorization" of Neumann internal problem at the end of phase1
    call par_preconditioner_dd_mlevel_bddc_fill_val_phase_1 (p_mat, mlbddc)

    ! "Numerical factorization of Coarse/Dirichlet internal problems at the end of phase 2
    call par_preconditioner_dd_mlevel_bddc_fill_val_phase_2 (p_mat, mlbddc) 

    if ( i_am_coarse_task .or. i_am_higher_level_task ) then
       if (  mlbddc%co_sys_sol_strat == recursive_bddc ) then
          !call par_preconditioner_dd_mlevel_bddc_fill_val ( mlbddc%p_mat_c, mlbddc%p_M_c ) 
          call par_preconditioner_dd_mlevel_bddc_fill_val ( mlbddc%p_M_c ) 
       end if
    end if

  end subroutine par_preconditioner_dd_mlevel_bddc_fill_val

  !=================================================================================================
  subroutine par_preconditioner_dd_mlevel_bddc_fill_val_phase_1 ( p_mat, mlbddc )
                                           
    implicit none
    ! Parameters 
    type(par_matrix_t)         , intent(in)                    :: p_mat
    type(par_preconditioner_dd_mlevel_bddc_t), intent(inout), target  :: mlbddc

    ! Locals
    integer(ip) :: me, np, lunou
    logical     :: i_am_fine_task
    type(par_matrix_t) :: p_mat_trans


    ! The routine requires the partition/context info
    assert ( associated(p_mat%p_env) )
    assert ( p_mat%p_env%created )
    assert ( p_mat%p_env%num_levels > 1 ) 

    assert ( associated(mlbddc%p_mat%p_env%w_context) )
    assert ( mlbddc%p_mat%p_env%w_context%created .eqv. .true. )
    assert ( associated(mlbddc%p_mat%p_env%p_context) )
    assert ( mlbddc%p_mat%p_env%p_context%created .eqv. .true. )
    assert ( associated(mlbddc%p_mat%p_env%q_context) )
    assert ( mlbddc%p_mat%p_env%q_context%created .eqv. .true. )
    ! Check appropriate assignment to context: 1) I'm in w_context and 2) I'm in p_context OR in q_context but not in both
    assert ( mlbddc%p_mat%p_env%w_context%iam >= 0)
    assert ( (mlbddc%p_mat%p_env%p_context%iam >=0 .and. mlbddc%p_mat%p_env%q_context%iam < 0) .or. (mlbddc%p_mat%p_env%p_context%iam < 0 .and. mlbddc%p_mat%p_env%q_context%iam >= 0))
    assert ( mlbddc%c_context%created .eqv. .true. )
    assert ( mlbddc%d_context%created .eqv. .true. )
    assert ( mlbddc%g_context%created .eqv. .true. )
    assert ( associated(mlbddc%p_mat%p_env%b_context) )
    assert ( mlbddc%p_mat%p_env%b_context%created .eqv. .true. )

    ! Which duties I have?
    i_am_fine_task = (mlbddc%p_mat%p_env%p_context%iam >= 0)

    if (i_am_fine_task) then
       ! BEG. FINE-GRID PROBLEM DUTIES
       call matrix_fill_val ( mlbddc%A_rr )
      
       if (mlbddc%projection == petrov_galerkin ) then 
       call matrix_fill_val ( mlbddc%A_rr_trans )
       end if

       if ( mlbddc%nn_sys_sol_strat == corners_rest_part_solve_expl_schur ) then

          call memalloc ( mlbddc%A_rr_gr%nv,      &
                          mlbddc%nl_corners_dofs, &
                          mlbddc%A_rc, __FILE__,__LINE__ )

          call memalloc ( mlbddc%nl_corners_dofs, &
                          mlbddc%A_rr_gr%nv,      & 
                          mlbddc%A_cr, __FILE__,__LINE__ ) 

          call memalloc ( mlbddc%nl_corners_dofs, &
                          mlbddc%nl_corners_dofs, &
                          mlbddc%A_cc, __FILE__,__LINE__ )

          if (mlbddc%projection == petrov_galerkin ) then 
          call memalloc ( mlbddc%A_rr_gr%nv,      &
                          mlbddc%nl_corners_dofs, &
                          mlbddc%A_cr_trans, __FILE__,__LINE__ )

          call memalloc ( mlbddc%nl_corners_dofs, &
                          mlbddc%A_rr_gr%nv,      & 
                          mlbddc%A_rc_trans, __FILE__,__LINE__ ) 

          call memalloc ( mlbddc%nl_corners_dofs, &
                          mlbddc%nl_corners_dofs, &
                          mlbddc%A_cc_trans, __FILE__,__LINE__ )
          end if
       

          ! Extract a_rc/a_cr/a_rr/a_cc 
          call extract_values_A_rr_A_cr_A_rc_A_cc  ( mlbddc%nl_corners_dofs, & 
                                                     p_mat%f_matrix%gr%type   , & 
                                                     p_mat%f_matrix%gr%nv     , & 
                                                     p_mat%f_matrix%gr%ia     , & 
                                                     p_mat%f_matrix%gr%ja     , &
                                                     p_mat%f_matrix%a         , &
                                                     mlbddc%A_rr_gr%nv   , & 
                                                     mlbddc%A_rr_gr%ia   , &
                                                     mlbddc%A_rr_gr%ja   , &
                                                     mlbddc%perm         , &
                                                     mlbddc%iperm        , &       
                                                     mlbddc%A_rr%a       , &
                                                     mlbddc%A_cr, & 
                                                     mlbddc%A_rc, &              
                                                     mlbddc%A_cc )

       
          if ( mlbddc%projection == petrov_galerkin ) then 

           call matrix_transpose( p_mat%f_matrix, p_mat_trans%f_matrix )       

      call extract_values_A_rr_A_cr_A_rc_A_cc  ( mlbddc%nl_corners_dofs, & 
                                                 p_mat_trans%f_matrix%gr%type   , & 
                                                 p_mat_trans%f_matrix%gr%nv     , & 
                                                 p_mat_trans%f_matrix%gr%ia     , & 
                                                 p_mat_trans%f_matrix%gr%ja     , &
                                                 p_mat_trans%f_matrix%a         , &
                                                 mlbddc%A_rr_gr%nv   , & 
                                                 mlbddc%A_rr_gr%ia   , &
                                                 mlbddc%A_rr_gr%ja   , &
                                                 mlbddc%perm         , &
                                                 mlbddc%iperm        , &       
                                                 mlbddc%A_rr_trans%a       , &
                                                 mlbddc%A_rc_trans, & 
                                                 mlbddc%A_cr_trans, &              
                                                 mlbddc%A_cc_trans )

           call matrix_free(p_mat_trans%f_matrix)
          
         end if

       else if ( mlbddc%nn_sys_sol_strat == direct_solve_constrained_problem) then 

          call augment_matrix_with_constraints ( p_mat%f_matrix%gr%type, & 
                                                 p_mat%f_matrix%gr%nv   , & 
                                                 p_mat%f_matrix%gr%ia    , & 
                                                 p_mat%f_matrix%gr%ja    , &
                                                 p_mat%f_matrix%a    , &
                                                 mlbddc%A_rr_gr%nv   , & 
                                                 mlbddc%A_rr_gr%ia   , &
                                                 mlbddc%A_rr_gr%ja, &
                                                 mlbddc%p_mat%dof_dist%nl, &
                                                 mlbddc%nl_corners, &
                                                 mlbddc%nl_edges, &
                                                 mlbddc%nl_coarse, &
                                                 mlbddc%coarse_dofs, &
                                                 mlbddc%p_mat%dof_dist%max_nparts, &
                                                 mlbddc%p_mat%dof_dist%omap%nl, &
                                                 mlbddc%p_mat%dof_dist%lobjs, &
                                                 mlbddc%A_rr%a,  &
                                                 mlbddc%p_mat%dof_dist%nb, &
                                                 mlbddc%C_weights )

         if (mlbddc%projection == petrov_galerkin ) then 
              
          call matrix_transpose( p_mat%f_matrix, p_mat_trans%f_matrix )  

          call augment_matrix_with_constraints ( p_mat_trans%f_matrix%gr%type, & 
                                                 p_mat_trans%f_matrix%gr%nv   , & 
                                                 p_mat_trans%f_matrix%gr%ia    , & 
                                                 p_mat_trans%f_matrix%gr%ja    , &
                                                 p_mat_trans%f_matrix%a    , &
                                                 mlbddc%A_rr_gr%nv   , & 
                                                 mlbddc%A_rr_gr%ia   , &
                                                 mlbddc%A_rr_gr%ja, &
                                                 mlbddc%p_mat%dof_dist%nl, &
                                                 mlbddc%nl_corners, &
                                                 mlbddc%nl_edges, &
                                                 mlbddc%nl_coarse, &
                                                 mlbddc%coarse_dofs, &
                                                 mlbddc%p_mat%dof_dist%max_nparts, &
                                                 mlbddc%p_mat%dof_dist%omap%nl, &
                                                 mlbddc%p_mat%dof_dist%lobjs, &
                                                 mlbddc%A_rr_trans%a, &
                                                 mlbddc%p_mat%dof_dist%nb, &
                                                 mlbddc%C_weights  )
        
            call matrix_free(p_mat_trans%f_matrix)

         end if

       end if

       if ( debug_verbose_level_3 ) then
          call par_context_info ( mlbddc%p_mat%p_env%p_context, me, np )
          lunou = io_open ( trim('fine_matrix_mlevel_bddc_' //  trim(ch(mlbddc%p_mat%p_env%num_levels)) // trim('+') // trim(ch(me+1)) // trim('.') // 'mtx' ), 'write')
          call matrix_print_matrix_market ( lunou, mlbddc%A_rr )
          call io_close ( lunou )
       end if

       if ( mlbddc%internal_problems == handled_by_bddc_module) then
          ! call preconditioner_numeric( mlbddc%A_rr , mlbddc%M_rr )
          call preconditioner_numeric( mlbddc%M_rr )
          
          if (mlbddc%projection == petrov_galerkin ) then 
             !call preconditioner_numeric( mlbddc%A_rr_trans, mlbddc%M_rr_trans ) 
             call preconditioner_numeric( mlbddc%M_rr_trans ) 
          end if

       else
          check(.false.)
!!$          call operator_mat_create (mlbddc%A_rr, mlbddc%A_rr_mat_op )
!!$          call abs_operator_convert(mlbddc%A_rr_mat_op, mlbddc%A_rr_op)
       end if
    end if


  end subroutine par_preconditioner_dd_mlevel_bddc_fill_val_phase_1

  !=================================================================================================
  subroutine par_preconditioner_dd_mlevel_bddc_fill_val_phase_2 ( p_mat, mlbddc, realloc_harm_extensions )
                                           
    implicit none
    ! Parameters 
    type(par_matrix_t)         , intent(in)                    :: p_mat
    type(par_preconditioner_dd_mlevel_bddc_t), intent(inout), target  :: mlbddc
    logical, optional, intent(in)                            :: realloc_harm_extensions

    ! Locals 
    real(rp) , allocatable  :: rhs(:,:), lambda_r(:,:), a_ci(:,:)
    real(rp) , allocatable  :: rhs_t(:,:), lambda_r_t(:,:), a_ci_t(:,:)
    real (rp), allocatable  :: work (:,:), work_t(:,:)

    integer          :: iam, me, np
    integer          :: lunou
    logical          :: i_am_coarse_task, i_am_fine_task, i_am_higher_level_task

    ! The routine requires the partition/context info
    assert ( associated(p_mat%p_env) )
    assert ( p_mat%p_env%created )
    assert ( p_mat%p_env%num_levels > 1 ) 

    assert ( mlbddc%p_mat%p_env%w_context%created .eqv. .true. )
    assert ( associated(mlbddc%p_mat%p_env%p_context) )
    assert ( mlbddc%p_mat%p_env%p_context%created .eqv. .true. )
    assert ( associated(mlbddc%p_mat%p_env%q_context) )
    assert ( mlbddc%p_mat%p_env%q_context%created .eqv. .true. )
    ! Check appropriate assignment to context: 1) I'm in w_context and 2) I'm in p_context OR in q_context but not in both
    assert ( mlbddc%p_mat%p_env%w_context%iam >= 0)
    assert ( (mlbddc%p_mat%p_env%p_context%iam >=0 .and. mlbddc%p_mat%p_env%q_context%iam < 0) .or. (mlbddc%p_mat%p_env%p_context%iam < 0 .and. mlbddc%p_mat%p_env%q_context%iam >= 0))
    assert ( mlbddc%c_context%created .eqv. .true. )
    assert ( mlbddc%d_context%created .eqv. .true. )
    assert ( mlbddc%g_context%created .eqv. .true. )
    assert ( associated(mlbddc%p_mat%p_env%b_context) )
    assert ( mlbddc%p_mat%p_env%b_context%created .eqv. .true. )
    
    ! Which duties I have?
    i_am_fine_task = (mlbddc%p_mat%p_env%p_context%iam >= 0)
    i_am_coarse_task = (mlbddc%c_context%iam >= 0)
    i_am_higher_level_task = (mlbddc%d_context%iam >= 0)

    if (i_am_fine_task) then
   
       if ( mlbddc%nn_sys_sol_strat == corners_rest_part_solve_expl_schur ) then

          if ( mlbddc%kind_coarse_dofs == corners_and_edges .or. & 
               mlbddc%kind_coarse_dofs == corners_edges_and_faces .or. &
               mlbddc%kind_coarse_dofs == edges .or. & 
               mlbddc%kind_coarse_dofs == edges_and_faces .or. & 
               mlbddc%kind_coarse_dofs == faces .or. & 
               mlbddc%kind_coarse_dofs == corners_and_faces ) then

             call memalloc ( mlbddc%nl_edges, & 
                             mlbddc%nl_edges+mlbddc%nl_corners, & 
                             rhs, __FILE__,__LINE__ )

             call memalloc ( mlbddc%A_rr_gr%nv, & 
                             mlbddc%nl_corners, & 
                             work, __FILE__,__LINE__ )

             if ( mlbddc%projection == petrov_galerkin ) then 

                call memalloc ( mlbddc%nl_edges, & 
                                mlbddc%nl_edges+mlbddc%nl_corners, & 
                                rhs_t, __FILE__,__LINE__ )

                call memalloc ( mlbddc%A_rr_gr%nv, & 
                                mlbddc%nl_corners, & 
                                work_t, __FILE__,__LINE__ )
             end if


              call compute_edge_lagrange_multipliers_rhs (mlbddc, & 
                                                            mlbddc%A_rr, &
                                                            mlbddc%M_rr, & 
                                                            mlbddc%a_rc, & 
                                                            rhs, work)

             if ( mlbddc%projection == petrov_galerkin ) then 
                call compute_edge_lagrange_multipliers_rhs (mlbddc, & 
                                                            mlbddc%A_rr_trans, &
                                                            mlbddc%M_rr_trans, & 
                                                            mlbddc%a_cr_trans, & 
                                                            rhs_t, work_t)
             end if
          
             call memalloc ( mlbddc%nl_edges, &
                             mlbddc%nl_edges+mlbddc%nl_corners, &
                             lambda_r, __FILE__,__LINE__ )

             if ( mlbddc%projection == petrov_galerkin ) then 
             call memalloc ( mlbddc%nl_edges, &
                             mlbddc%nl_edges+mlbddc%nl_corners, &
                             lambda_r_t, __FILE__,__LINE__ )
             end if
 
             call compute_edge_lagrange_multipliers_schur_complement ( mlbddc, & 
                                                                       mlbddc%A_rr, &
                                                                       mlbddc%M_rr, &
                                                                       mlbddc%S_rr, &
                                                                       mlbddc%S_rr_neumann, &
                                                                       mlbddc%A_rr_inv_C_r_T, &
                                                                       mlbddc%A_rr_inv_C_r_T_neumann, &
                                                                       mlbddc%ipiv, &
                                                                       mlbddc%ipiv_neumann ) 

                if (mlbddc%projection == petrov_galerkin ) then 
             call compute_edge_lagrange_multipliers_schur_complement ( mlbddc, & 
                                                                       mlbddc%A_rr_trans, &
                                                                       mlbddc%M_rr_trans, &
                                                                       mlbddc%S_rr_trans, &
                                                                       mlbddc%S_rr_trans_neumann, &
                                                                       mlbddc%A_rr_trans_inv_C_r_T, &
                                                                       mlbddc%A_rr_trans_inv_C_r_T_neumann, &
                                                                       mlbddc%ipiv_trans, &
                                                                       mlbddc%ipiv_trans_neumann ) 
                end if
             

              call solve_edge_lagrange_multipliers_explicit_schur ( mlbddc, &
                                                                    mlbddc%S_rr, &
                                                                    mlbddc%S_rr_neumann, &
                                                                    mlbddc%ipiv, &
                                                                    mlbddc%ipiv_neumann, &
                                                                    reuse_from_phis, &
                                                                    (mlbddc%nl_edges+mlbddc%nl_corners),&
                                                                    rhs, lambda_r )
                 call memfree ( rhs, __FILE__, __LINE__)

              if (mlbddc%projection == petrov_galerkin ) then 
              
              call solve_edge_lagrange_multipliers_explicit_schur ( mlbddc, &
                                                                      mlbddc%S_rr_trans, &
                                                                      mlbddc%S_rr_trans_neumann, &
                                                                      mlbddc%ipiv_trans, &
                                                                      mlbddc%ipiv_trans_neumann, &
                                                                      reuse_from_phis, &
                                                                     (mlbddc%nl_edges+mlbddc%nl_corners),&
                                                                      rhs_t, lambda_r_t )
                call memfree( rhs_t, __FILE__, __LINE__ )

              end if

          else
             call memalloc ( mlbddc%nl_edges, 0, lambda_r, __FILE__,__LINE__ )   

             call memalloc ( mlbddc%A_rr_gr%nv, 0, work, __FILE__,__LINE__ )    

           if (mlbddc%projection == petrov_galerkin ) then 

             call memalloc ( mlbddc%nl_edges, 0, lambda_r_t, __FILE__,__LINE__ )   

             call memalloc ( mlbddc%A_rr_gr%nv, 0, work_t, __FILE__,__LINE__ ) 
           
           end if
               
          end if

          call memalloc ( mlbddc%p_mat%dof_dist%nl, & 
                          mlbddc%nl_edges+mlbddc%nl_corners, & 
                          mlbddc%rPhi, __FILE__, __LINE__ )

          call memalloc ( mlbddc%nl_edges+mlbddc%nl_corners, &
                          mlbddc%nl_edges+mlbddc%nl_corners, & 
                          mlbddc%blk_lag_mul, __FILE__,__LINE__ )

          call memalloc ( mlbddc%p_mat%dof_dist%nl, & 
                          mlbddc%nl_edges+mlbddc%nl_corners, & 
                          mlbddc%lPhi, __FILE__, __LINE__ )



         call compute_harmonic_extensions_corner_edges_partitioning  ( mlbddc, &
                                                                       mlbddc%A_cr, &
                                                                       mlbddc%A_rc, &
                                                                       mlbddc%A_cc, &
                                                                       mlbddc%A_rr, &
                                                                       mlbddc%M_rr, &
                                                                       mlbddc%A_rr_inv_C_r_T, &
                                                                       mlbddc%rPhi, &
                                                                       lambda_r, work, 'N' )

           if (mlbddc%projection == petrov_galerkin ) then 

          call compute_harmonic_extensions_corner_edges_partitioning ( mlbddc, &
                                                                       mlbddc%A_rc_trans, &
                                                                       mlbddc%A_cr_trans, &
                                                                       mlbddc%A_cc_trans, &
                                                                       mlbddc%A_rr_trans, &
                                                                       mlbddc%M_rr_trans, &
                                                                       mlbddc%A_rr_trans_inv_C_r_T, &
                                                                       mlbddc%lPhi, &
                                                                       lambda_r_t, work_t, 'T' )

           end if


          call memfree ( lambda_r,__FILE__,__LINE__)

          call memfree ( mlbddc%A_rc,__FILE__,__LINE__)
          call memfree ( mlbddc%A_cr,__FILE__,__LINE__)
          call memfree ( mlbddc%A_cc,__FILE__,__LINE__)

          call memfree ( work,__FILE__,__LINE__)
          
          if (mlbddc%projection == petrov_galerkin )  then 
          call memfree ( lambda_r_t,__FILE__,__LINE__)

          call memfree ( mlbddc%A_rc_trans,__FILE__,__LINE__)
          call memfree ( mlbddc%A_cr_trans,__FILE__,__LINE__)
          call memfree ( mlbddc%A_cc_trans,__FILE__,__LINE__)

          call memfree ( work_t,__FILE__,__LINE__)

          end if

       else if ( mlbddc%nn_sys_sol_strat == direct_solve_constrained_problem ) then
          call memalloc ( mlbddc%p_mat%dof_dist%nl, &
                          mlbddc%nl_edges+mlbddc%nl_corners, & 
                          mlbddc%rPhi, __FILE__,__LINE__ )

          call memalloc ( mlbddc%p_mat%dof_dist%nl, &
                          mlbddc%nl_edges+mlbddc%nl_corners, & 
                          mlbddc%lPhi, __FILE__,__LINE__ )

          call memalloc ( mlbddc%nl_edges+mlbddc%nl_corners, & 
                          mlbddc%nl_edges+mlbddc%nl_corners, & 
                          mlbddc%blk_lag_mul, __FILE__,__LINE__ )


          call compute_harmonic_extensions_with_constraints (mlbddc, & 
                                                             mlbddc%A_rr, &
                                                             mlbddc%M_rr, &
                                                             mlbddc%rPhi, & 
                                                             'N' ) 

         if (mlbddc%projection == petrov_galerkin ) then 

         call compute_harmonic_extensions_with_constraints (mlbddc, & 
                                                            mlbddc%A_rr_trans, &
                                                            mlbddc%M_rr_trans, &
                                                            mlbddc%lPhi, & 
                                                            'T' ) 

         end if

       end if
       ! END FINE-GRID PROBLEM DUTIES
       mlbddc%harm_its_agg = mlbddc%spars_harm%it
    end if

!!$    if ( temp_fine_coarse_grid_overlap ) then 
!!$       if ( .not. i_am_higher_level_task ) then
!!$          if (i_am_fine_task) then
!!$             call par_timer_stop ( mlbddc%timer_assprec_ov_fine  ) 
!!$          else if ( i_am_coarse_task ) then
!!$             call par_timer_stop ( mlbddc%timer_assprec_ov_coarse ) 
!!$          end if
!!$       end if
!!$    end if

!!$    if ( temp_fine_coarse_grid_overlap ) then 
!!$       if (.not. i_am_higher_level_task) then
!!$          call par_timer_start ( mlbddc%timer_coll_fillprec  )
!!$       end if
!!$    end if
    
    ! BEG. COARSE-GRID PROBLEM DUTIES

! If the projection is symmetric (Galerkin), lPhi is created but not filled, and now is assigned to rPhi
    if (mlbddc%p_mat%p_env%p_context%iam >= 0) then
       if (mlbddc%projection == galerkin) then 
          mlbddc%lPhi = mlbddc%rPhi 
       end if
    end if
 
    call assemble_A_c ( p_mat, mlbddc, realloc_harm_extensions=realloc_harm_extensions )
   
   ! If recursive BDDC and enable C_weights: 
   ! Project and assemble constraints weights and put it up on the next level.   
    if ( mlbddc%co_sys_sol_strat == recursive_bddc ) then
       if (mlbddc%enable_constraint_weights) then  
         call assemble_C_weights ( mlbddc )
       end if 
    end if
   
!!$    if ( temp_fine_coarse_grid_overlap ) then 
!!$       if (.not. i_am_higher_level_task) then
!!$          call par_timer_stop ( mlbddc%timer_coll_fillprec  )
!!$       end if
!!$    end if
    
    if ( .not. i_am_fine_task ) then ! ( i_am_coarse_task .or. i_am_higher_level_task )
       if (  mlbddc%co_sys_sol_strat == serial_gather ) then  ! There are only coarse tasks
          if ( debug_verbose_level_3 ) then
             lunou =  io_open ( trim('A_mlevel_bddc_c.mtx'), 'write')  
             call matrix_print_matrix_market (lunou,  mlbddc%A_c)
             call io_close (lunou)
          end if
          if ( mlbddc%internal_problems == handled_by_bddc_module) then
             ! call preconditioner_numeric( mlbddc%A_c , mlbddc%M_c )
             call preconditioner_numeric( mlbddc%M_c )
          else
             check(.false.)
!!$             call operator_mat_create (mlbddc%A_c, mlbddc%A_c_mat_op )
!!$             call abs_operator_convert(mlbddc%A_c_mat_op, mlbddc%A_c_op)
          end if
       end if
    end if

    ! END COARSE-GRID PROBLEM DUTIES

    if (i_am_fine_task) then

       if ( temp_fine_coarse_grid_overlap ) then
          call par_timer_start ( mlbddc%timer_fillprec_ov_fine_header )
       end if
       
       if ( mlbddc%unknowns == all_unknowns ) then
          ! BEG. FINE-GRID PROBLEM DUTIES
          call operator_dd_fill_val (p_mat%f_matrix, mlbddc%A_II_inv) 
          ! END. FINE-GRID PROBLEM DUTIES 
          
          if ( debug_verbose_level_3 ) then
             call par_context_info ( mlbddc%p_mat%p_env%p_context, me, np )
             lunou = io_open ( trim('dirichlet_matrix_mlevel_bddc_' // trim(ch(mlbddc%p_mat%p_env%num_levels)) // trim('+') //  trim(ch(me+1)) // trim('.') // 'mtx' ), 'write')
             call matrix_print_matrix_market ( lunou, mlbddc%A_II_inv%A_II )
             call io_close ( lunou )
          end if
       end if
       
       if (temp_fine_coarse_grid_overlap) then
          call par_timer_stop ( mlbddc%timer_fillprec_ov_fine_header  )
       end if

    end if

  end subroutine par_preconditioner_dd_mlevel_bddc_fill_val_phase_2


  subroutine extract_values_A_rr_A_cr_A_rc_A_cc ( nl_corners_dofs, gtype,  & 
       anv, aia, aja, a,        &  
       arrnv, arria, arrja,     &
       perm, iperm, arr, acr, arc, acc  )
    implicit none

    ! Parameters
    integer (ip) , intent(in)    :: nl_corners_dofs
    integer (ip) , intent(in)    :: gtype, anv, arrnv
    integer (ip) , intent(in)    :: aia (anv+1)
    integer (ip) , intent(in)    :: aja (aia(anv+1)-1)
    real    (rp) , intent(in)    :: a   (aia(anv+1)-1)  
    integer (ip) , intent(inout) :: arria (arrnv+1)
    integer (ip) , intent(inout) :: arrja (arria(arrnv+1)-1)
    integer (ip) , intent(in)    :: perm  (anv)
    integer (ip) , intent(in)    :: iperm (anv)
    real    (rp) , intent(out)   :: arr (arria(arrnv+1)-1)
    real    (rp) , intent(out)   :: acr (nl_corners_dofs, arrnv)
    real    (rp) , intent(out)   :: arc (arrnv, nl_corners_dofs) 
    real    (rp) , intent(out)   :: acc (nl_corners_dofs, nl_corners_dofs)

    ! Locals
    integer(ip) :: id, jd, idarr, ida, jdarr, jda, jdapos, i
    integer(ip) :: k1, k2, k, j,l
    integer(ip), allocatable :: work(:)
    real   (rp), allocatable :: awork(:)

    acr = 0.0_rp
    arc = 0.0_rp
    acc = 0.0_rp

    call memalloc ( arrnv, work, __FILE__,__LINE__ )
    call memalloc ( arrnv, awork, __FILE__,__LINE__ )


    ! Traverse first nl_corners_dofs of P^T A_i P 
    do id = 1, nl_corners_dofs
       ida = iperm(id)
       do jdapos = aia(ida), aia(ida+1)-1
          jd = perm( aja(jdapos) )
          ! write (*,*) 'ZZZ1', id, jd, ida, aja(jdapos), aia(ida), aia(ida+1)-1  
          if (jd > nl_corners_dofs .and. gtype == csr_symm) then
             ! write (*,*) 'XXX', idarr, jdarr
             acr (id, jd-nl_corners_dofs) = a(jdapos)
             arc (jd-nl_corners_dofs, id) = a(jdapos)
             !   elseif (jd > nl_corners_dofs) then
             !      acr (id,jd-nl_corners_dofs) = a(jdapos) ! may be obtained here, well defined down
          else ! jd <= nl_corners_dofs .or. gtype != csr_symm
             if ( jd <= nl_corners_dofs ) then
                if (gtype == csr_symm) then
                   acc(id,jd) = a(jdapos)
                   acc(jd,id) = a(jdapos)
                else
                   acc(id,jd) = a(jdapos)
                end if
             end if
          end if
       end do
    end do

    ! write (*,*) acc

    do idarr=1, arrnv
       ida = iperm(idarr+nl_corners_dofs)
       k1  = arria (idarr)
       k2  = arria (idarr+1)-1
       do jdapos = aia(ida), aia(ida+1)-1
          jdarr = perm( aja(jdapos) )
          ! write (*,*) 'ZZZ2', idarr+nl_corners_dofs, jdarr, nl_corners_dofs
          if ( jdarr > nl_corners_dofs ) then
             jdarr = jdarr - nl_corners_dofs
             if ( gtype == csr ) then
                arrja ( arria(idarr) ) = jdarr
                arr   ( arria(idarr) ) = a(jdapos)
                arria (idarr) = arria(idarr) + 1
             else if ( gtype == csr_symm ) then
                if ( jdarr >= idarr ) then
                   arrja ( arria(idarr) ) = jdarr
                   arr   ( arria(idarr) ) = a(jdapos) 
                   arria(idarr) = arria(idarr) + 1
                else ! jdarr < idarr
                   arrja ( arria(jdarr) ) = idarr
                   arr  ( arria(jdarr) ) = a(jdapos)
                   arria(jdarr) = arria(jdarr) + 1
                end if
             end if
          else
             ! write (*,*) 'XXX', idarr, jdarr
             arc (idarr, jdarr) =  a(jdapos) 
             acr (jdarr, idarr) =  a(jdapos)
          end if
       end do

       if ( gtype == csr ) then
          ! Order increasingly column identifiers of current row 
          ! using heap sort algorithm
          k = k2-k1+1
          do i=1,k
             work(i) = i
          end do

          call sort( k2-k1+1, arrja(k1:k2), work(1:k) )

          awork (1:k) = arr (k1:k2)
          do i=0, k-1
             arr ( k1+i ) = awork(work(i+1))
          end do
       end if
    end do

    do idarr=arrnv+1, 2, -1
       arria(idarr) = arria(idarr-1)
    end do
    arria(idarr) = 1

    if ( gtype == csr_symm ) then
       do idarr=1, arrnv
          k1  = arria (idarr)
          k2  = arria (idarr+1)-1
          k   = k2-k1+1
          ! Order increasingly column identifiers of current row 
          ! using heap sort algorithm
          do i=1,k
             work(i) = i
          end do
          ! write (*,*) 'XXA', arrja(k1:k2)
          ! write (*,*) 'YYA', work (1:k)              
          call sort( k2-k1+1, arrja(k1:k2), work(1:k) )
          awork (1:k) = arr (k1:k2)
          ! write (*,*) 'sdasd', awork(1:k) 
          do i=0, k-1
             arr ( k1+i ) = awork( work(i+1) )
          end do
          ! write (*,*) 'sdasd', arr (k1:k2)
          ! write (*,*) 'XXD', arrja(k1:k2)
          ! write (*,*) 'YYD', work (1:k) 
       end do
    end if

    !   Acr extraction 
    if (gtype == csr) then
       do idarr= 1,nl_corners_dofs
          ida = iperm(idarr)
          do jdapos = aia(ida), aia(ida+1)-1
             jdarr = perm( aja(jdapos) )

             if (jdarr > nl_corners_dofs) then 
                acr(idarr,jdarr-nl_corners_dofs) = a(jdapos)
             end if
          end do
       end do
    end if

    call memfree  ( work,__FILE__,__LINE__)
    call memfree  ( awork,__FILE__,__LINE__)
  end subroutine extract_values_A_rr_A_cr_A_rc_A_cc



  subroutine extract_trans_matrix(A, gr_a, A_t, gr_t)

    type(matrix_t), intent(in)     :: A         ! Input matrix
    type(graph_t),  target, intent(in)     :: gr_a      ! Graph of the matrix
    type(matrix_t), intent(out)  :: A_t       ! Output matrix
    type(graph_t), pointer,  intent(out)  :: gr_t      ! Graph of the transpose matrix

    ! Locals 
    type(graph_t) :: aux_graph
    integer :: k,i,j

    if ( gr_a%type == csr ) then
       call matrix_alloc ( csr_mat, symm_false, gr_a, A_t )
    else if (gr_a%type == csr_symm) then
       call matrix_alloc ( csr_mat, symm_true, gr_a, A_t )
    end if

    aux_graph = gr_a

    if (gr_a%type == csr_symm) then 
       A_t%a = A%a
    elseif (gr_a%type == csr ) then
       k = 0    
       do i = 1,gr_t%nv
          do j=1, (gr_t%ia(i+1) - gr_t%ia(i) )
             k = k+1
             A_t%a(k) = A%a(aux_graph%ia( gr_t%ja(k) ) )
             aux_graph%ia(gr_t%ja(k)) = aux_graph%ia(gr_t%ja(k)) + 1
          end do
       end do

    end if

  end subroutine extract_trans_matrix

  subroutine augment_matrix_with_constraints ( gtype, & 
                                               anv, aia, aja, a, &  
                                               arr_nv, arr_ia, arr_ja, &
                                               nl, nl_corners, nl_edges, nl_coarse, & 
                                               coarse_dofs, max_nparts, nobjs, lobjs, &
                                               arr, nb, C_weights )
    implicit none

    ! Parameters
    integer (ip) , intent(in)   :: gtype, anv, arr_nv
    integer (ip) , intent(in)   :: aia (anv+1)
    integer (ip) , intent(in)   :: aja (aia(anv+1)-1)
    real    (rp) , intent(in)   :: a   (aia(anv+1)-1)  
    integer (ip) , intent(inout):: arr_ia (arr_nv+1)
    integer (ip) , intent(in)   :: nl
    integer (ip) , intent(in)   :: nl_corners, nl_edges, nl_coarse
    integer (ip) , intent(in)   :: coarse_dofs (nl_coarse) 
    integer (ip) , intent(in)   :: max_nparts, nobjs
    integer (ip) , intent(in)   :: lobjs(max_nparts+4, nobjs)
    integer (ip) , intent(in)   :: arr_ja (arr_ia(arr_nv+1)-1)
    real    (rp) , intent(out)  :: arr (arr_ia(arr_nv+1)-1)
    integer (ip) ,  intent(in), optional :: nb
    real    (rp) ,  intent(in), allocatable, optional :: C_weights(:)

    ! Locals
    integer(ip) :: idarr, nz, size, iposedge, icoarse, nid
    integer(ip) :: i, j, k, iobj, jd, j1, j2
    integer(ip) :: js, ji1, ji2, je
    real (rp)   :: weight, w_size
    logical     :: aux


    ! Number of interior nodes
    if ( present (C_weights) )  nid = nl - nb 

    ! List elements on each row of 
    ! a and map them to the corresponding 
    ! rows of a_rr
    do idarr=1, anv
       nz = aia(idarr+1) - aia(idarr)  
       arr ( arr_ia(idarr):arr_ia(idarr)+nz-1 ) = a ( aia(idarr):aia(idarr+1)-1 )
       arr_ia ( idarr ) = arr_ia( idarr ) + nz
    end do

    if ( gtype == csr ) then
       do i=1, nl_coarse 
          iobj     = coarse_dofs (i)
          iposedge = i-nl_corners

            ! Corner as a degenerated edge or
            ! face without extra_corners
            j1 = lobjs(2,iobj)
            j2 = lobjs(3,iobj)

            size = 0 
            if ( present (C_weights) ) then 
               if (allocated (C_weights)) then 
                  w_size = 0
               end if
            end if

            do jd = j1, j2
               size = size + 1
               if ( present (C_weights) ) then   
                 if (allocated (C_weights) ) then                
                   w_size = w_size +  C_weights(jd-nid)
                end if
             end if
            end do

            weight = 1.0_rp/size

            do jd = j1, j2
               ! new entry (jd,i+nl)
               if ( present (C_weights) ) then 
                  if ( allocated (C_weights) ) then 
                     weight = C_weights(jd-nid)/w_size
                  end if
               end if
               arr( arr_ia(jd) ) = weight
               arr_ia(jd)        = arr_ia(jd) + 1

               ! new entry (id+nl,jd)
               arr ( arr_ia(i+nl) ) = weight 
               arr_ia ( i+nl )      = arr_ia ( i+nl ) + 1 
            end do

       end do
    else if ( gtype == csr_symm ) then ! Add explicit zeros on the diagonal
       do i=1, nl_coarse
          iobj   = coarse_dofs (i)
          iposedge = i - nl_corners
             
          ! Corner as a degenerated edge or
          ! face without extra_corners
          j1     = lobjs(2,iobj)
          j2     = lobjs(3,iobj)

          size = 0
          if ( present (C_weights) ) w_size = 0

          do jd = j1, j2
             size = size + 1
             if ( present (C_weights) ) then 
               if (allocated (C_weights) ) then 
                w_size = w_size + C_weights(jd-nid)
             end if
             end if
          end do
                
          weight = 1.0_rp/size

          do jd = j1, j2
             ! new entry (jd,id+nl)
             if ( present (C_weights) ) then 
                if ( allocated (C_weights) ) then 
                   weight = C_weights(jd-nid)/w_size
                end if
             end if
             arr( arr_ia(jd) ) = weight
             arr_ia(jd)        = arr_ia(jd) + 1
          end do
                       
       end do
    end if


    ! Add explicit zeros on the diagonal
    do idarr=anv+1,arr_nv
       arr( arr_ia(idarr) ) = 0.0_rp
       arr_ia(idarr) = arr_ia(idarr) + 1 
    end do

    
    ! Recover arr_ia
    do idarr = arr_nv + 1, 2, -1
       arr_ia (idarr) = arr_ia (idarr-1)
    end do
    arr_ia (idarr) = 1
        
  end subroutine augment_matrix_with_constraints


 ! Computes the several RHS's of the edge lagrange multipliers Schur complement
  ! 'N' system : (C_r A_rr^-1 C_r^T) lamba_r = - ( [0 I_r] + C_r * A_rr ^-1 [A_rc 0] )
  ! 'T' system : (C_r A_rr_trans^-1 C_r^T) lamba_r_t = - ( [0 I_r] + C_r * A_rr_trans ^-1 [A_cr_trans 0] )
  subroutine compute_edge_lagrange_multipliers_rhs (mlbddc, A_rr, M_rr,  A_rc, rhs, work)
    implicit none
    ! Parameters 
    type(par_preconditioner_dd_mlevel_bddc_t), intent(inout), target :: mlbddc

    type(matrix_t)       , intent(inout)            :: A_rr
    type(preconditioner_t)      , intent(inout)            :: M_rr

    real(rp)               , intent(inout)         :: A_rc ( mlbddc%A_rr_gr%nv  , &
                                                               mlbddc%nl_corners) 

    real(rp)               , intent(out)           :: rhs ( mlbddc%nl_edges, & 
                                                              (mlbddc%nl_edges+mlbddc%nl_corners))

    real(rp)               , intent(out)           :: work ( mlbddc%A_rr_gr%nv, & 
                                                               mlbddc%nl_corners)
    ! Locals
    integer(ip)             :: base, i, start, icoarse

   ! type(solver_control_t)      :: spars
 

    if ( mlbddc%nl_corners > 0 ) then
       mlbddc%spars_harm%nrhs=mlbddc%nl_corners
       ! mlbddc%spars%method=direct
       work = 0.0_rp
       if ( mlbddc%internal_problems == handled_by_bddc_module) then
             call solve(A_rr, M_rr, &
                        A_rc, mlbddc%A_rr_gr%nv,  &
                        work, mlbddc%A_rr_gr%nv,  &
                        mlbddc%spars_harm)  
       else
          check(.false.)
!!$          call solve( mlbddc%A_rr_op, mlbddc%M_inv_op_rr, &
!!$               A_rc, mlbddc%A_rr_gr%nv,  &
!!$               work, mlbddc%A_rr_gr%nv,  &
!!$               mlbddc%spars_harm)
       end if


    end if

    rhs = 0.0_rp
    base = mlbddc%nl_corners 

    ! rhs = - ( [0 I_r] ) 
    do i=1, mlbddc%nl_edges
        rhs(i, base+i) = -1.0_rp
    end do
    

    if ( mlbddc%nl_corners > 0 ) then
       ! Y = alpha * Y + beta * C_r * X
       ! C_r => nl_edges x n(A_rr)

       call update_with_c_r   ( mlbddc%p_mat%dof_dist%nl, &
                                mlbddc%perm, &
                                mlbddc%iperm, &
                                mlbddc%nl_corners, &
                                mlbddc%nl_edges, &
                                mlbddc%coarse_dofs, &
                                mlbddc%p_mat%dof_dist%max_nparts, &
                                mlbddc%p_mat%dof_dist%omap%nl, &
                                mlbddc%p_mat%dof_dist%lobjs, &
                                mlbddc%nl_edges, & 
                                mlbddc%A_rr_gr%nv, &
                                mlbddc%nl_corners, & 
                                work, & 
                                1.0_rp, &
                                -1.0_rp, & 
                                rhs, &
                                mlbddc%p_mat%dof_dist%nb, &
                                mlbddc%C_weights )
    end if
    
     ! write (*,*) 'XXX', rhs
  end subroutine compute_edge_lagrange_multipliers_rhs


  ! Y = alpha * Y + beta * C_r * X
  ! C_r => nl_edges x n(A_rr)
  subroutine update_with_c_r   ( npoin, perm, iperm,           &
                                 nl_corners, nl_edges, coarse_dofs,  & 
                                 max_nparts, nobjs, lobjs,           &
                                 m, n, k, X, alpha, beta, Y, nb, C_weights)
    implicit none
    integer (ip), intent(in)    :: npoin
    integer (ip), intent(in)    :: perm(npoin), iperm(npoin)
    integer (ip), intent(in)    :: nl_corners, nl_edges
    integer (ip), intent(in)    :: coarse_dofs(nl_corners+nl_edges) 
    integer (ip), intent(in)    :: max_nparts, nobjs
    integer (ip), intent(in)    :: lobjs(max_nparts+4, nobjs)
    integer (ip), intent(in)    :: m, k, n
    real(rp)  , intent(in)    :: X(n,k)
    real(rp)  , intent(in)    :: alpha, beta
    real(rp)  , intent(inout) :: Y(m,k)   
    integer(ip) , intent(in), optional :: nb
    real(rp)  ,   intent(in), allocatable, optional :: C_weights(:)

    ! Locals
    integer (ip) :: r, i, j, iobj, size, jd, l, ld, nid
    integer (ip) :: js, je, p
    real (rp)    :: weight, sum, w_size

!!$    OPTIMIZATION !!!
!!$    ! Y = alpha * Y
!!$    if (alpha == 0.0_rp) then
!!$
!!$    else if (alpha == 1.0_rp) then
!!$
!!$    else ! alpha /= 0.0_rp and alpha /= 1.0_rp
!!$
!!$    end if

    if (present(C_weights)) nid = npoin - nb

    do r = 1, k
       do i=1, nl_edges 
          iobj   = coarse_dofs (nl_corners + i)

          size = 0 
          if (present(C_weights)) w_size = 0

          sum  = 0.0_rp

          js   = lobjs(2,iobj) ! Start
          je   = lobjs(3,iobj) ! End

          
          do jd = js, je
             size = size + 1
             if (present(C_weights)) then
              if (allocated(C_weights)) then  
                w_size = w_size + C_weights(jd-nid)
             end if
             end if
          end do


          do jd = js, je
              weight = 1.0

             if (present(C_weights)) then 
                if (allocated(C_weights)) then 
                   weight = C_weights(jd-nid)/w_size
                end if
             end if

             ld = perm (jd) - nl_corners
             sum = sum + weight*X(ld,r)
          end do

              weight = 1.0_rp/size
          if (present(C_weights)) then 
             if (allocated(C_weights)) then 
                weight = 1.0
             end if
          end if
          Y(i,r) = alpha * Y(i,r) + beta * weight * sum

       end do
    end do

  end subroutine update_with_c_r

  ! Y = alpha * Y + beta * C_r^T * X
  subroutine update_with_c_r_trans  ( npoin, perm, iperm,           &
                                      nl_corners, nl_edges, coarse_dofs,  & 
                                      max_nparts, nobjs, lobjs,           &
                                      m, n, k, X, alpha, beta, Y, nb, C_weights)
    implicit none

    ! Parameters
    integer (ip), intent(in)    :: npoin
    integer (ip), intent(in)    :: perm(npoin), iperm(npoin)
    integer (ip), intent(in)    :: nl_corners, nl_edges
    integer (ip), intent(in)    :: coarse_dofs(nl_corners+nl_edges) 
    integer (ip), intent(in)    :: max_nparts, nobjs
    integer (ip), intent(in)    :: lobjs(max_nparts+4, nobjs)
    integer (ip), intent(in)    :: m, n, k
    real(rp)  , intent(in)    :: X(n,k)
    real(rp)  , intent(in)    :: alpha, beta
    real(rp)  , intent(inout) :: Y(m,k)
    integer(ip) , intent(in), optional :: nb
    real(rp)  ,   intent(in), allocatable, optional :: C_weights(:)

    ! Locals
    integer (ip) :: r, i, j, iobj, size, w_size, id, jd, jd1, jd2, j1, j2, jdof, ld, l, nid
    integer (ip) :: js, je, ji1, ji2, iposedge, p
    real (rp)    :: weight, sum, xirw    

!!$    OPTIMIZATION !!!
!!$    ! Y = alpha * Y
!!$    if (alpha == 0.0_rp) then
!!$
!!$    else if (alpha == 1.0_rp) then
!!$
!!$    else ! alpha /= 0.0_rp and alpha /= 1.0_rp
!!$
!!$    end if

! Number of interior DOFs
  if (present(C_weights)) nid = npoin - nb

    do r = 1, k
      ! write (*,*) 'XXX', r, n/ndof   
      do i=1, nl_edges 
         iobj   = coarse_dofs (nl_corners + i)

         size = 0 
         if (present(C_weights)) w_size = 0

         js   = lobjs(2,iobj) ! Start
         je   = lobjs(3,iobj) ! End
 
         do jd = js, je
            size = size + 1
         if (present(C_weights)) then
          if (allocated(C_weights)) then  
            w_size = w_size + C_weights(jd-nid)
          end if
         end if
         end do 
            
         weight = 1.0_rp/size
         xirw = X(id,r)*weight

         do jd = js, je
            ld = perm (jd) - nl_corners
            if (present(C_weights)) then 
               if (allocated(C_weights)) then 
                  weight = C_weights(jd-nid)/w_size
                  xirw = X(id,r)*weight
               end if
            end if
            Y(ld,r) = alpha*Y(ld,r) + beta * xirw
         end do

       end do
    end do

  end subroutine update_with_c_r_trans

  ! Y = alpha * Y + beta * C_i^T * X
  subroutine update_with_c_i_trans  ( npoin,                        &
                                      nl_corners, nl_edges, coarse_dofs,  & 
                                      max_nparts, nobjs, lobjs,           &
                                      m, n, k, X, alpha, beta, Y, C_weights)
    implicit none

    ! Parameters
    integer (ip), intent(in)    :: npoin
    integer (ip), intent(in)    :: nl_corners, nl_edges
    integer (ip), intent(in)    :: coarse_dofs(nl_corners+nl_edges) 
    integer (ip), intent(in)    :: max_nparts, nobjs
    integer (ip), intent(in)    :: lobjs(max_nparts+4, nobjs)
    integer (ip), intent(in)    :: m, n, k
    real(rp)  , intent(in)    :: X(n,k)
    real(rp)  , intent(in)    :: alpha, beta
    real(rp)  , intent(inout) :: Y(m,k)
    real(rp)  ,   intent(in), allocatable, optional :: C_weights(:)

    ! Locals
    integer (ip) :: r, i, j, iobj, size, id, jd, j1, j2, ld, l, nid
    integer (ip) :: js, je, ji1, ji2, iposedge, p
    real (rp)    :: weight, sum, xirw, w_size    

    ! Compute interior dofs
    nid = npoin - m

    do r = 1, k

      do i=1, nl_corners
          iobj   = coarse_dofs (i)
          ! Corner as a degenerated edge or
          ! edge/face
          j1     = lobjs(2,iobj)
          j2     = lobjs(3,iobj)
          ! size   = j2-j1+1
          ! weight = (1.0_rp/size)

          size = 0    

          do jd = j1, j2
             size = size + 1
          end do

          weight = 1.0_rp/size

          xirw = X(i,r)*weight
          ! write (*,*) jd1, jd2, id, beta*X(id,r)
          do jd = j1, j2
             Y(jd-nid,r) = alpha*Y(jd-nid,r) + beta * xirw 
          end do

       end do


      ! write (*,*) 'XXX', r, n/ndof   
      do i=1, nl_edges 
         iobj   = coarse_dofs (nl_corners + i)
         id     = nl_corners+i
         js   = lobjs(2,iobj) ! Start
         je   = lobjs(3,iobj) ! End

         size = 0 
         if (present(C_weights)) w_size = 0

         do jd = js, je
            size = size + 1
            if ( present (C_weights) ) then 
             if (allocated(C_weights)) then 
                w_size = w_size + C_weights(jd-nid)
             end if  
            end if
         end do

         weight = 1.0_rp/size
         xirw = X(id,r)*weight

         ! write (*,*) 'AAA', iposedge, js, je, ji2, je , size, weight, X(id,r)
         do jd = js, je
            ld = jd
            if (present(C_weights)) then 
             if (allocated(C_weights)) then 
               weight = C_weights(jd-nid)/w_size
               xirw   = X(id,r)*weight
             end if 
            end if
            Y(ld-nid,r) = alpha*Y(ld-nid,r) + beta * xirw
         end do

     !    id  = id + 1

      end do
    end do


  end subroutine update_with_c_i_trans

!!$      SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
!!$*
!!$*  -- LAPACK routine (version 3.3.1) --
!!$*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!!$*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!!$*  -- April 2011                                                      --
!!$*
!!$*     .. Scalar Arguments ..
!!$      CHARACTER          UPLO
!!$      INTEGER            INFO, LDA, N
!!$*     ..
!!$*     .. Array Arguments ..
!!$      DOUBLE PRECISION   A( LDA, * )
!!$*     ..
!!$*
!!$*  Purpose
!!$*  =======
!!$*
!!$*  DPOTRF computes the Cholesky factorization of a real symmetric
!!$*  positive definite matrix A.
!!$*
!!$*  The factorization has the form
!!$*     A = U**T * U,  if UPLO = 'U', or
!!$*     A = L  * L**T,  if UPLO = 'L',
!!$*  where U is an upper triangular matrix and L is lower triangular.
!!$*
!!$*  This is the block version of the algorithm, calling Level 3 BLAS.
!!$*
!!$*  Arguments
!!$*  =========
!!$*
!!$*  UPLO    (input) CHARACTER*1
!!$*          = 'U':  Upper triangle of A is stored;
!!$*          = 'L':  Lower triangle of A is stored.
!!$*
!!$*  N       (input) INTEGER
!!$*          The order of the matrix A.  N >= 0.
!!$*
!!$*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!!$*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!!$*          N-by-N upper triangular part of A contains the upper
!!$*          triangular part of the matrix A, and the strictly lower
!!$*          triangular part of A is not referenced.  If UPLO = 'L', the
!!$*          leading N-by-N lower triangular part of A contains the lower
!!$*          triangular part of the matrix A, and the strictly upper
!!$*          triangular part of A is not referenced.
!!$*
!!$*          On exit, if INFO = 0, the factor U or L from the Cholesky
!!$*          factorization A = U**T*U or A = L*L**T.
!!$*
!!$*  LDA     (input) INTEGER
!!$*          The leading dimension of the array A.  LDA >= max(1,N).
!!$*
!!$*  INFO    (output) INTEGER
!!$*          = 0:  successful exit
!!$*          < 0:  if INFO = -i, the i-th argument had an illegal value
!!$*          > 0:  if INFO = i, the leading minor of order i is not
!!$*                positive definite, and the factorization could not be
!!$*                completed.
!!$*
!!$*  =====================================================================

!!$      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
!!$*
!!$*  -- LAPACK routine (version 3.2) --
!!$*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!!$*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!!$*     November 2006
!!$*
!!$*     .. Scalar Arguments ..
!!$      INTEGER            INFO, LDA, M, N
!!$*     ..
!!$*     .. Array Arguments ..
!!$      INTEGER            IPIV( * )
!!$      DOUBLE PRECISION   A( LDA, * )
!!$*     ..
!!$*
!!$*  Purpose
!!$*  =======
!!$*
!!$*  DGETRF computes an LU factorization of a general M-by-N matrix A
!!$*  using partial pivoting with row interchanges.
!!$*
!!$*  The factorization has the form
!!$*     A = P * L * U
!!$*  where P is a permutation matrix, L is lower triangular with unit
!!$*  diagonal elements (lower trapezoidal if m > n), and U is upper
!!$*  triangular (upper trapezoidal if m < n).
!!$*
!!$*  This is the right-looking Level 3 BLAS version of the algorithm.
!!$*
!!$*  Arguments
!!$*  =========
!!$*
!!$*  M       (input) INTEGER
!!$*          The number of rows of the matrix A.  M >= 0.
!!$*
!!$*  N       (input) INTEGER
!!$*          The number of columns of the matrix A.  N >= 0.
!!$*
!!$*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!!$*          On entry, the M-by-N matrix to be factored.
!!$*          On exit, the factors L and U from the factorization
!!$*          A = P*L*U; the unit diagonal elements of L are not stored.
!!$*
!!$*  LDA     (input) INTEGER
!!$*          The leading dimension of the array A.  LDA >= max(1,M).
!!$*
!!$*  IPIV    (output) INTEGER array, dimension (min(M,N))
!!$*          The pivot indices; for 1 <= i <= min(M,N), row i of the
!!$*          matrix was interchanged with row IPIV(i).
!!$*
!!$*  INFO    (output) INTEGER
!!$*          = 0:  successful exit
!!$*          < 0:  if INFO = -i, the i-th argument had an illegal value
!!$*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!!$*                has been completed, but the factor U is exactly
!!$*                singular, and division by zero will occur if it is used
!!$*                to solve a system of equations.
!!$*
!!$*  =====================================================================


!!$      INTEGER          FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3,
!!$     $                 N4 )
!!$*
!!$*  -- LAPACK auxiliary routine (preliminary version) --
!!$*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!!$*     Courant Institute, Argonne National Lab, and Rice University
!!$*     February 20, 1992
!!$*
!!$*     .. Scalar Arguments ..
!!$      CHARACTER*( * )    NAME, OPTS
!!$      INTEGER            ISPEC, N1, N2, N3, N4
!!$*     ..
!!$*
!!$*  Purpose
!!$*  =======
!!$*
!!$*  ILAENV is called from the LAPACK routines to choose problem-dependent
!!$*  parameters for the local environment.  See ISPEC for a description of
!!$*  the parameters.
!!$*
!!$*  This version provides a set of parameters which should give good,
!!$*  but not optimal, performance on many of the currently available
!!$*  computers.  Users are encouraged to modify this subroutine to set
!!$*  the tuning parameters for their particular machine using the option
!!$*  and problem size information in the arguments.
!!$*
!!$*  This routine will not function correctly if it is converted to all
!!$*  lower case.  Converting it to all upper case is allowed.
!!$*
!!$*  Arguments
!!$*  =========
!!$*
!!$*  ISPEC   (input) INTEGER
!!$*          Specifies the parameter to be returned as the value of
!!$*          ILAENV.
!!$*          = 1: the optimal blocksize; if this value is 1, an unblocked
!!$*               algorithm will give the best performance.
!!$*          = 2: the minimum block size for which the block routine
!!$*               should be used; if the usable block size is less than
!!$*               this value, an unblocked routine should be used.
!!$*          = 3: the crossover point (in a block routine, for N less
!!$*               than this value, an unblocked routine should be used)
!!$*          = 4: the number of shifts, used in the nonsymmetric
!!$*               eigenvalue routines
!!$*          = 5: the minimum column dimension for blocking to be used;
!!$*               rectangular blocks must have dimension at least k by m,
!!$*               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!!$*          = 6: the crossover point for the SVD (when reducing an m by n
!!$*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!!$*               this value, a QR factorization is used first to reduce
!!$*               the matrix to a triangular form.)
!!$*          = 7: the number of processors
!!$*          = 8: the crossover point for the multishift QR and QZ methods
!!$*               for nonsymmetric eigenvalue problems.
!!$*
!!$*  NAME    (input) CHARACTER*(*)
!!$*          The name of the calling subroutine, in either upper case or
!!$*          lower case.
!!$*
!!$*  OPTS    (input) CHARACTER*(*)
!!$*          The character options to the subroutine NAME, concatenated
!!$*          into a single character string.  For example, UPLO = 'U',
!!$*          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!!$*          be specified as OPTS = 'UTN'.
!!$*
!!$*  N1      (input) INTEGER
!!$*  N2      (input) INTEGER
!!$*  N3      (input) INTEGER
!!$*  N4      (input) INTEGER
!!$*          Problem dimensions for the subroutine NAME; these may not all
!!$*          be required.
!!$*
!!$* (ILAENV) (output) INTEGER
!!$*          >= 0: the value of the parameter specified by ISPEC
!!$*          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!!$*
!!$*  Further Details
!!$*  ===============
!!$*
!!$*  The following conventions have been used when calling ILAENV from the
!!$*  LAPACK routines:
!!$*  1)  OPTS is a concatenation of all of the character options to
!!$*      subroutine NAME, in the same order that they appear in the
!!$*      argument list for NAME, even if they are not used in determining
!!$*      the value of the parameter specified by ISPEC.
!!$*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!!$*      that they appear in the argument list for NAME.  N1 is used
!!$*      first, N2 second, and so on, and unused problem dimensions are
!!$*      passed a value of -1.
!!$*  3)  The parameter value returned by ILAENV is checked for validity in
!!$*      the calling subroutine.  For example, ILAENV is used to retrieve
!!$*      the optimal blocksize for STRTRI as follows:
!!$*
!!$*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!!$*      IF( NB.LE.1 ) NB = MAX( 1, N )
!!$*
!!$*  ====================================================================

!!$      SUBROUTINE DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
!!$*
!!$*  -- LAPACK routine (version 3.1) --
!!$*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!!$*     November 2006
!!$*
!!$*     .. Scalar Arguments ..
!!$      CHARACTER          UPLO
!!$      INTEGER            INFO, LDA, LWORK, N
!!$*     ..
!!$*     .. Array Arguments ..
!!$      INTEGER            IPIV( * )
!!$      DOUBLE PRECISION   A( LDA, * ), WORK( * )
!!$*     ..
!!$*
!!$*  Purpose
!!$*  =======
!!$*
!!$*  DSYTRF computes the factorization of a real symmetric matrix A using
!!$*  the Bunch-Kaufman diagonal pivoting method.  The form of the
!!$*  factorization is
!!$*
!!$*     A = U*D*U**T  or  A = L*D*L**T
!!$*
!!$*  where U (or L) is a product of permutation and unit upper (lower)
!!$*  triangular matrices, and D is symmetric and block diagonal with
!!$*  1-by-1 and 2-by-2 diagonal blocks.
!!$*
!!$*  This is the blocked version of the algorithm, calling Level 3 BLAS.
!!$*
!!$*  Arguments
!!$*  =========
!!$*
!!$*  UPLO    (input) CHARACTER*1
!!$*          = 'U':  Upper triangle of A is stored;
!!$*          = 'L':  Lower triangle of A is stored.
!!$*
!!$*  N       (input) INTEGER
!!$*          The order of the matrix A.  N >= 0.
!!$*
!!$*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!!$*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!!$*          N-by-N upper triangular part of A contains the upper
!!$*          triangular part of the matrix A, and the strictly lower
!!$*          triangular part of A is not referenced.  If UPLO = 'L', the
!!$*          leading N-by-N lower triangular part of A contains the lower
!!$*          triangular part of the matrix A, and the strictly upper
!!$*          triangular part of A is not referenced.
!!$*
!!$*          On exit, the block diagonal matrix D and the multipliers used
!!$*          to obtain the factor U or L (see below for further details).
!!$*
!!$*  LDA     (input) INTEGER
!!$*          The leading dimension of the array A.  LDA >= max(1,N).
!!$*
!!$*  IPIV    (output) INTEGER array, dimension (N)
!!$*          Details of the interchanges and the block structure of D.
!!$*          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!!$*          interchanged and D(k,k) is a 1-by-1 diagonal block.
!!$*          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
!!$*          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
!!$*          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
!!$*          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
!!$*          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!!$*
!!$*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!!$*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!!$*
!!$*  LWORK   (input) INTEGER
!!$*          The length of WORK.  LWORK >=1.  For best performance
!!$*          LWORK >= N*NB, where NB is the block size returned by ILAENV.
!!$*
!!$*          If LWORK = -1, then a workspace query is assumed; the routine
!!$*          only calculates the optimal size of the WORK array, returns
!!$*          this value as the first entry of the WORK array, and no error
!!$*          message related to LWORK is issued by XERBLA.
!!$*
!!$*  INFO    (output) INTEGER
!!$*          = 0:  successful exit
!!$*          < 0:  if INFO = -i, the i-th argument had an illegal value
!!$*          > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization
!!$*                has been completed, but the block diagonal matrix D is
!!$*                exactly singular, and division by zero will occur if it
!!$*                is used to solve a system of equations.
!!$*
!!$*  Further Details
!!$*  ===============
!!$*
!!$*  If UPLO = 'U', then A = U*D*U', where
!!$*     U = P(n)*U(n)* ... *P(k)U(k)* ...,
!!$*  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
!!$*  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!!$*  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!!$*  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
!!$*  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!!$*
!!$*             (   I    v    0   )   k-s
!!$*     U(k) =  (   0    I    0   )   s
!!$*             (   0    0    I   )   n-k
!!$*                k-s   s   n-k
!!$*
!!$*  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
!!$*  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
!!$*  and A(k,k), and v overwrites A(1:k-2,k-1:k).
!!$*
!!$*  If UPLO = 'L', then A = L*D*L', where
!!$*     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
!!$*  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
!!$*  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!!$*  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!!$*  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
!!$*  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!!$*
!!$*             (   I    0     0   )  k-1
!!$*     L(k) =  (   0    I     0   )  s
!!$*             (   0    v     I   )  n-k-s+1
!!$*                k-1   s  n-k-s+1
!!$*
!!$*  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
!!$*  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
!!$*  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
!!$*
!!$*  =====================================================================

  ! Explicitly computes the Schur complement corresponding to edge-lagrange multipliers
  ! i.e., S_rr = C_r k_rr^-1 C_r^T, as well as A_rr^-1 C_r^T on A_rr_inv_C_r_T (:,:).
  ! Then, the schur complement it is factorized/overwritten on the same memory
  ! space using LAPACK.
  subroutine compute_edge_lagrange_multipliers_schur_complement ( mlbddc, A_rr, M_rr,  S_rr, S_rr_neumann, A_rr_inv_C_r_T, A_rr_inv_C_r_T_neumann, ipiv, ipiv_neumann,  iparm, msglvl )
use mpi
    implicit none
    ! Parameters 
    type(par_preconditioner_dd_mlevel_bddc_t), intent(inout) :: mlbddc
    type(matrix_t)          , intent(inout)   :: A_rr
    type(preconditioner_t)         , intent(inout)   :: M_rr

    real(rp), allocatable     , intent(inout)   :: A_rr_inv_C_r_T(:,:) , A_rr_inv_C_r_T_neumann(:,:)
    real(rp), allocatable     , intent(inout)   :: S_rr(:,:), S_rr_neumann(:,:)
    integer (ip), allocatable , intent(inout)   :: ipiv(:), ipiv_neumann(:) 
    
    integer                , intent(in), target, optional :: iparm(64)
    integer                , intent(in), optional         :: msglvl

    ! Locals
    integer (ip) :: error, start, icoarse, j
    integer :: info
    real (rp), allocatable :: C_r_T (:,:)
    real (rp), allocatable :: work(:)
    integer (ip)           :: lwork
   ! type(solver_control_t)      :: spars

#ifdef ENABLE_LAPACK
    call memalloc (  mlbddc%nl_edges, mlbddc%nl_edges,  S_rr,  __FILE__,__LINE__ )

    call memalloc ( mlbddc%A_rr_gr%nv,   mlbddc%nl_edges, A_rr_inv_C_r_T,    __FILE__,__LINE__ )

    if (mlbddc%schur_edge_lag_mult == compute_from_scratch ) then
       call memalloc (  mlbddc%nl_edges, mlbddc%nl_edges,  S_rr_neumann,  __FILE__,__LINE__ )
       call memalloc (  mlbddc%A_rr_gr%nv,   mlbddc%nl_edges,  A_rr_inv_C_r_T_neumann,    __FILE__,__LINE__ )
    end if

    call memalloc (  mlbddc%A_rr_gr%nv,    mlbddc%nl_edges, C_r_T, __FILE__,__LINE__ )


    ! Compute C_r_T
    call compute_c_r_trans ( mlbddc%p_mat%dof_dist%nl, &
                             mlbddc%perm, &
                             mlbddc%iperm, &
                             mlbddc%nl_corners, &
                             mlbddc%nl_edges, &
                             mlbddc%coarse_dofs, &
                             mlbddc%p_mat%dof_dist%max_nparts, &
                             mlbddc%p_mat%dof_dist%omap%nl, &
                             mlbddc%p_mat%dof_dist%lobjs, &
                             C_r_T, &  
                             mlbddc%p_mat%dof_dist%nb, &
                             mlbddc%C_weights )


    !!$ call mpi_barrier (mpi_comm_world, info)
    !!$ check(1==0)

    if ( mlbddc%nl_edges > 0 ) then
       ! Compute A_rr^-1 * C_r^T
       mlbddc%spars_harm%nrhs=mlbddc%nl_edges
       ! mlbddc%spars%method=direct

       A_rr_inv_C_r_T = 0.0_rp
       if ( mlbddc%internal_problems == handled_by_bddc_module) then
          call solve( A_rr, M_rr, &
                      C_r_T, mlbddc%A_rr_gr%nv, &
                      A_rr_inv_C_r_T, mlbddc%A_rr_gr%nv, &
                      mlbddc%spars_harm)
       else 
          check(.false.)
!!$          call solve( mlbddc%A_rr_op, mlbddc%M_inv_op_rr, &
!!$               C_r_T, mlbddc%A_rr_gr%nv, &
!!$               A_rr_inv_C_r_T, mlbddc%A_rr_gr%nv, &
!!$               mlbddc%spars_harm)
       end if

       if (mlbddc%schur_edge_lag_mult == compute_from_scratch ) then
          ! Compute A_rr^-1 * C_r^T
          mlbddc%spars_neumann%nrhs = mlbddc%nl_edges
          ! mlbddc%spars%method=direct

          A_rr_inv_C_r_T_neumann = 0.0_rp
          if ( mlbddc%internal_problems == handled_by_bddc_module) then
             call solve( A_rr, M_rr, & 
                         C_r_T, mlbddc%A_rr_gr%nv, &
                         A_rr_inv_C_r_T_neumann, mlbddc%A_rr_gr%nv, &
                         mlbddc%spars_neumann)
          else
             check(.false.)
!!$             call solve( mlbddc%A_rr_op, mlbddc%M_inv_op_rr, & 
!!$                  C_r_T, mlbddc%A_rr_gr%nv, &
!!$                  A_rr_inv_C_r_T_neumann, mlbddc%A_rr_gr%nv, &
!!$                  mlbddc%spars_neumann)
          end if
       end if

    end if

    ! Compute S_rr <-  C_r A_rr^-1 C_r^T
    ! Y = alpha * Y + beta * C_r * X
    ! C_r => nl_edges x n(A_rr)
    S_rr = 0.0

    call update_with_c_r   ( mlbddc%p_mat%dof_dist%nl, &
                             mlbddc%perm, &
                             mlbddc%iperm, &
                             mlbddc%nl_corners, &
                             mlbddc%nl_edges, &
                             mlbddc%coarse_dofs, &
                             mlbddc%p_mat%dof_dist%max_nparts, &
                             mlbddc%p_mat%dof_dist%omap%nl, &
                             mlbddc%p_mat%dof_dist%lobjs, &
                             mlbddc%nl_edges, & 
                             mlbddc%A_rr_gr%nv, &
                             mlbddc%nl_edges, & 
                             A_rr_inv_C_r_T, & 
                             0.0_rp, &
                             1.0_rp, & 
                             S_rr,  & 
                             mlbddc%p_mat%dof_dist%nb, &
                             mlbddc%C_weights )


    if (mlbddc%schur_edge_lag_mult == compute_from_scratch ) then
       ! Compute S_rr <-  C_r A_rr^-1 C_r^T
       ! Y = alpha * Y + beta * C_r * X
       ! C_r => nl_edges x n(A_rr)
       S_rr_neumann = 0.0

       call update_with_c_r   ( mlbddc%p_mat%dof_dist%nl, &
                                mlbddc%perm, &
                                mlbddc%iperm, &
                                mlbddc%nl_corners, &
                                mlbddc%nl_edges, &
                                mlbddc%coarse_dofs, &
                                mlbddc%p_mat%dof_dist%max_nparts, &
                                mlbddc%p_mat%dof_dist%omap%nl, &
                                mlbddc%p_mat%dof_dist%lobjs, &
                                mlbddc%nl_edges, & 
                                mlbddc%A_rr_gr%nv, &
                                mlbddc%nl_edges, & 
                                A_rr_inv_C_r_T_neumann, & 
                                0.0_rp, &
                                1.0_rp, & 
                                S_rr_neumann, &
                                mlbddc%p_mat%dof_dist%nb, &
                                mlbddc%C_weights ) 

    end if


    call memfree ( C_r_T,__FILE__,__LINE__)

    if (mlbddc%symm == symm_true) then
       
       if ( mlbddc%sign == positive_definite ) then
          
          if ( mlbddc%nl_edges > 0 ) then 
             call DPOTRF ( 'U', &
                  mlbddc%nl_edges, &
                  S_rr, &
                  mlbddc%nl_edges, &
                  error )
             
             if (error /= 0) then
                write (0,*) 'Error, DOPTRF: the following ERROR was detected: ', error
                check(1==0)
             end if
          end if
          
       else if ( mlbddc%sign == indefinite .or. mlbddc%sign == unknown ) then
          call memalloc (  mlbddc%nl_edges, ipiv, __FILE__,__LINE__ )
          
          if ( mlbddc%nl_edges > 0 ) then 
             lwork = ILAENV ( 1, 'DSYTRF', 'U', mlbddc%nl_edges, -1, -1, -1 )
             if ( lwork < 0 ) then
                write (0,*) 'Error, ILAENV: illegal value in parameter number: ', lwork
                check(1==0)
             end if
             lwork = mlbddc%nl_edges
             call memalloc (  lwork, work,  __FILE__,__LINE__ )
             
             call DSYTRF ( 'U', &
                  mlbddc%nl_edges, &
                  S_rr, &
                  mlbddc%nl_edges, &
                  ipiv,&
                  work, &
                  lwork, &
                  error )
             
             if (error /= 0) then
                write (0,*) 'Error, DSYTRF: the following ERROR was detected: ', error
                check(1==0)
             end if
             
             call memfree ( work,__FILE__,__LINE__)
          end if
       end if

    else if (mlbddc%symm == symm_false) then

       call memalloc (  mlbddc%nl_edges, ipiv, __FILE__,__LINE__ )

       if ( mlbddc%nl_edges > 0 ) then 

          ! write (*,*) mlbddc%nl_edges  
          ! write (*,*) mlbddc%S_rr

          call DGETRF(  mlbddc%nl_edges, & 
               mlbddc%nl_edges, &
               S_rr, &
               mlbddc%nl_edges, &
               ipiv, &
               error )

          ! write (*,*) mlbddc%ipiv

          if (error /= 0) then
             write (0,*) 'Error, DGETRF: the following ERROR was detected: ', error
             check(1==0)
          end if
       end if
    end if

    if (mlbddc%schur_edge_lag_mult == compute_from_scratch ) then
       if (mlbddc%symm == symm_true) then
          if ( mlbddc%sign == positive_definite ) then
              if ( mlbddc%nl_edges > 0 ) then 
                call DPOTRF ( 'U', &
                    mlbddc%nl_edges, &
                     S_rr_neumann, &
                     mlbddc%nl_edges, &
                     error )

                if (error /= 0) then
                   write (0,*) 'Error, DOPTRF: the following ERROR was detected: ', error
                   check(1==0)
                end if
             end if
          else if ( mlbddc%sign == indefinite .or. mlbddc%sign == unknown ) then
             call memalloc (  mlbddc%nl_edges, & 
                  &  ipiv_neumann, & 
                  &  'compute_edge_lagrange_multipliers_schur_complement::mlbddc%ipiv_neumann' )

             if ( mlbddc%nl_edges > 0 ) then 
                lwork = ILAENV ( 1, 'DSYTRF', 'U', mlbddc%nl_edges, -1, -1, -1 )
                if ( lwork < 0 ) then
                   write (0,*) 'Error, ILAENV: illegal value in parameter number: ', lwork
                   check(1==0)
                end if
                lwork = mlbddc%nl_edges * lwork
                call memalloc (  lwork, & 
                     &  work,  & 
                     &  'compute_edge_lagrange_multipliers_schur_complement::work' )

                call DSYTRF ( 'U', &
                     mlbddc%nl_edges, &
                     S_rr_neumann, &
                     mlbddc%nl_edges, &
                     ipiv_neumann,&
                     work, &
                     lwork, &
                     error )

                if (error /= 0) then
                   write (0,*) 'Error, DSYTRF: the following ERROR was detected: ', error
                   check(1==0)
                end if

                call memfree ( work, & 
                     &  'compute_edge_lagrange_multipliers_schur_complement::work' )
             end if
          else
             write (0,*) 'Error, no dense factorization for symmetric positive_semidefinite matrices available'
             check(1==0)
          end if
       else if (mlbddc%symm == symm_false) then

          call memalloc (  mlbddc%nl_edges, & 
               &  ipiv_neumann, & 
               &  'compute_edge_lagrange_multipliers_schur_complement::mlbddc%ipiv_neumann' )
          
          if ( mlbddc%nl_edges > 0 ) then 

             ! write (*,*) mlbddc%nl_edges * mlbddc%nd1 
             ! write (*,*) mlbddc%S_rr_neumann

             call DGETRF(  mlbddc%nl_edges, & 
                  mlbddc%nl_edges, &
                  S_rr_neumann, &
                  mlbddc%nl_edges, &
                  ipiv_neumann, &
                  error )

             ! write (*,*) mlbddc%ipiv_neumann

             if (error /= 0) then
                write (0,*) 'Error, DGETRF: the following ERROR was detected: ', error
                check(1==0)
             end if
          end if

       end if
    end if
#else
    write (0,*) 'Error: par_preconditioner_dd_mlevel_bddc was not compiled with -DENABLE_LAPACK.'
    write (0,*) 'Error: You must activate this cpp macro in order to use the LAPACK'
    check(1==0)    
#endif
  end subroutine compute_edge_lagrange_multipliers_schur_complement


  subroutine compute_c_r_trans  ( npoin, perm, iperm,           &
       nl_corners, nl_edges, coarse_dofs,  & 
       max_nparts, nobjs, lobjs,           &
       Y, nb, C_weights)
    implicit none

    ! Parameters
    integer (ip), intent(in)    :: npoin
    integer (ip), intent(in)    :: perm(npoin), iperm(npoin)
    integer (ip), intent(in)    :: nl_corners, nl_edges
    integer (ip), intent(in)    :: coarse_dofs(nl_corners+nl_edges) 
    integer (ip), intent(in)    :: max_nparts, nobjs
    integer (ip), intent(in)    :: lobjs(max_nparts+4, nobjs)
    real(rp)  , intent(inout) :: Y( (npoin-nl_corners), nl_edges)
    integer(ip), intent(in), optional :: nb
    real(rp)   , intent(in), allocatable, optional :: C_weights(:)

    ! Locals
    integer (ip) :: i, j, iobj, size, id, jd, jd1, jd2, j1, j2, jdof, ld, l, nid
    integer (ip) :: js, je, ji1, ji2, iposedge, k
    real (rp)    :: weight, w_size 

    Y = 0.0

    if (present(C_weights)) nid = npoin - nb

    ! write (*,*) 'XXX', r, n/ndof   
    do i=1, nl_edges
       iobj   = coarse_dofs (nl_corners + i)

          size = 0
          if (present(C_weights)) w_size = 0

          js   = lobjs(2,iobj) ! Start
          je   = lobjs(3,iobj) ! End

          do jd = js, je
             size = size + 1
             if (present(C_weights)) then 
              if (allocated(C_weights)) then 
             w_size = w_size + C_weights(jd-nid)
              end if 
             end if
          end do

          weight = 1.0_rp/size

          do jd = js, je
             ld = perm (jd) - nl_corners
             if (present(C_weights)) then 
                if (allocated(C_weights)) then 
                   weight = C_weights(jd-nid)/w_size
                end if
             end if
             Y(ld,i) = weight
          end do
    end do

  end subroutine compute_c_r_trans

!!$      SUBROUTINE DSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!!$*
!!$*  -- LAPACK routine (version 3.1) --
!!$*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!!$*     November 2006
!!$*
!!$*     .. Scalar Arguments ..
!!$      CHARACTER          UPLO
!!$      INTEGER            INFO, LDA, LDB, N, NRHS
!!$*     ..
!!$*     .. Array Arguments ..
!!$      INTEGER            IPIV( * )
!!$      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!!$*     ..
!!$*
!!$*  Purpose
!!$*  =======
!!$*
!!$*  DSYTRS solves a system of linear equations A*X = B with a real
!!$*  symmetric matrix A using the factorization A = U*D*U**T or
!!$*  A = L*D*L**T computed by DSYTRF.
!!$*
!!$*  Arguments
!!$*  =========
!!$*
!!$*  UPLO    (input) CHARACTER*1
!!$*          Specifies whether the details of the factorization are stored
!!$*          as an upper or lower triangular matrix.
!!$*          = 'U':  Upper triangular, form is A = U*D*U**T;
!!$*          = 'L':  Lower triangular, form is A = L*D*L**T.
!!$*
!!$*  N       (input) INTEGER
!!$*          The order of the matrix A.  N >= 0.
!!$*
!!$*  NRHS    (input) INTEGER
!!$*          The number of right hand sides, i.e., the number of columns
!!$*          of the matrix B.  NRHS >= 0.
!!$*
!!$*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!!$*          The block diagonal matrix D and the multipliers used to
!!$*          obtain the factor U or L as computed by DSYTRF.
!!$*
!!$*  LDA     (input) INTEGER
!!$*          The leading dimension of the array A.  LDA >= max(1,N).
!!$*
!!$*  IPIV    (input) INTEGER array, dimension (N)
!!$*          Details of the interchanges and the block structure of D
!!$*          as determined by DSYTRF.
!!$*
!!$*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!!$*          On entry, the right hand side matrix B.
!!$*          On exit, the solution matrix X.
!!$*
!!$*  LDB     (input) INTEGER
!!$*          The leading dimension of the array B.  LDB >= max(1,N).
!!$*
!!$*  INFO    (output) INTEGER
!!$*          = 0:  successful exit
!!$*          < 0:  if INFO = -i, the i-th argument had an illegal value
!!$*
!!$*  =====================================================================


  ! Solves explicitly the edge lagrange multipliers Schur complement system:
  ! (C_r A_rr^-1 C_r^T) lamba_r = - ( [0 I_r] + C_r * A_rr ^-1 [A_rc 0] )
  subroutine solve_edge_lagrange_multipliers_explicit_schur (mlbddc, S_rr, S_rr_neumann, ipiv, ipiv_neumann, schur_edge_lag_mult, nrhs, rhs, lambda_r)
    implicit none

    ! Parameters 
    type(par_preconditioner_dd_mlevel_bddc_t), intent(in) :: mlbddc 

    real(rp)                  , intent(in)   :: S_rr(:,:), S_rr_neumann(:,:)
    integer (ip)              , intent(in)   :: ipiv(:), ipiv_neumann(:) 
    integer (ip)              , intent(in)   :: schur_edge_lag_mult
    integer (ip)              , intent(in)   :: nrhs
    real(rp)                  , intent(in), target  :: rhs (mlbddc%nl_edges, nrhs)
    real(rp)                  , intent(out), target :: lambda_r (mlbddc%nl_edges, nrhs)

    ! Locals
    real    (rp), allocatable :: aux (:)
    integer (ip)              :: i, error

#ifdef ENABLE_BLAS
    if ( mlbddc%nl_edges > 0 ) then
       if ( schur_edge_lag_mult == compute_from_scratch ) then
          if (mlbddc%symm == symm_true) then

             lambda_r = rhs 

             if ( mlbddc%sign == positive_definite ) then
                ! write (*,*) lambda_r(1:10,1)

                ! Fwd. substitution with U^T
                call DTRSM ( 'L', &
                     'U', & ! U is in the upper triangle
                     'T', & ! Tranpose U
                     'N', & ! U is not unit diagonal
                     mlbddc%nl_edges,&
                     nrhs,&
                     1.0,&
                     S_rr_neumann,&
                     mlbddc%nl_edges,&
                     lambda_r,&
                     mlbddc%nl_edges)

                ! write (*,*) lambda_r(1:10,1)

                ! Bck. substitution with U
                call DTRSM ( 'L', &
                     'U', & ! U is in the upper triangle
                     'N', & ! Do not transpose U
                     'N', & ! U is not unit diagonal 
                     mlbddc%nl_edges,&
                     nrhs,&
                     1.0,&
                     S_rr_neumann,&
                     mlbddc%nl_edges,&
                     lambda_r,&
                     mlbddc%nl_edges)

                ! write (*,*) lambda_r(1:10,1)
             else if ( mlbddc%sign == indefinite .or. mlbddc%sign == unknown ) then
              ! Bck. substitution with U
                call DSYTRS ( 'U', &
                     mlbddc%nl_edges,&
                     nrhs,&
                     S_rr_neumann,&
                     mlbddc%nl_edges,&
                     ipiv_neumann, &
                     lambda_r,&
                     mlbddc%nl_edges, &
                     error)
                if ( error < 0 ) then
                   write (0,*) 'Error, DSYTRS: illegal value in parameter number: ', error
                   check(1==0)
                end if
             else
                write (0,*) 'Error, no dense factorization for symmetric positive_semidefinite matrices available'
                check(1==0)
             end if

          else if (mlbddc%symm == symm_false) then

             call memalloc (  nrhs,       & 
                  &           aux ,       & 
                  &           'solve_edge_lagrange_multipliers_explicit_schur::aux' )

             ! Apply row interchanges to RHS
             lambda_r = rhs
             do i=1, mlbddc%nl_edges
                ! Swap lambda_r(i) and lambda_r (ipiv_neumann(i))
                aux = lambda_r (ipiv_neumann(i), : )
                lambda_r ( ipiv_neumann(i), : ) = lambda_r(i, :)
                lambda_r (i,:) = aux
             end do

             ! Fwd. substitution with L
             call DTRSM ( 'L', &
                  'L', &
                  'N', &  ! No transpose
                  'U', &  ! L is unit diagonal ...
                  mlbddc%nl_edges,&
                  nrhs,&
                  1.0,&
                  S_rr_neumann,&
                  mlbddc%nl_edges,&
                  lambda_r,&
                  mlbddc%nl_edges)

             ! Bck. substitution with U
             call DTRSM ( 'L', &
                  'U', &
                  'N', & ! Do not transpose
                  'N', & ! U is not unit diagonal
                  mlbddc%nl_edges,&
                  nrhs,&
                  1.0,&
                  S_rr_neumann,&
                  mlbddc%nl_edges,&
                  lambda_r,&
                  mlbddc%nl_edges)

             call memfree ( aux , 'solve_edge_lagrange_multipliers_explicit_schur::aux' )
          end if 
       else

          if (mlbddc%symm == symm_true) then

             lambda_r = rhs 


             if ( mlbddc%sign == positive_definite ) then

                ! write (*,*) lambda_r(1:10,1)

                ! Fwd. substitution with U^T
                call DTRSM ( 'L', &
                     'U', & ! U is in the upper triangle
                     'T', & ! Tranpose U
                     'N', & ! U is not unit diagonal
                     mlbddc%nl_edges,&
                     nrhs,&
                     1.0,&
                     S_rr,&
                     mlbddc%nl_edges,&
                     lambda_r,&
                     mlbddc%nl_edges)

                ! write (*,*) lambda_r(1:10,1)


                ! Bck. substitution with U
                call DTRSM ( 'L', &
                     'U', & ! U is in the upper triangle
                     'N', & ! Do not transpose U
                     'N', & ! U is not unit diagonal 
                     mlbddc%nl_edges  ,&
                     nrhs,&
                     1.0,&
                     S_rr,&
                     mlbddc%nl_edges ,&
                     lambda_r,&
                     mlbddc%nl_edges )

                ! write (*,*) lambda_r(1:10,1)
             else if ( mlbddc%sign == indefinite .or. mlbddc%sign == unknown ) then
                ! Bck. substitution with U
                call DSYTRS ( 'U', &
                     mlbddc%nl_edges,&
                     nrhs,&
                     S_rr,&
                     mlbddc%nl_edges,&
                     ipiv, &
                     lambda_r,&
                     mlbddc%nl_edges, &
                     error)
                if ( error < 0 ) then
                   write (0,*) 'Error, DSYTRS: illegal value in parameter number: ', error
                   check(1==0)
                end if
             else
                write (0,*) 'Error, no dense factorization for symmetric positive_semidefinite matrices available'
                check(1==0)
             end if

          else if (mlbddc%symm == symm_false) then

             call memalloc (  nrhs,       & 
                  &           aux ,       & 
                  &           'solve_edge_lagrange_multipliers_explicit_schur::aux' )

             ! Apply row interchanges to RHS
             lambda_r = rhs
             do i=1, mlbddc%nl_edges 
                ! Swap lambda_r(i) and lambda_r (ipiv(i))
                aux = lambda_r ( ipiv(i), : )
                lambda_r ( ipiv(i), : ) = lambda_r(i, :)
                lambda_r (i,:) = aux
             end do

             ! Fwd. substitution with L
             call DTRSM ( 'L', &
                  'L', &
                  'N', &  ! No transpose
                  'U', &  ! L is unit diagonal ...
                  mlbddc%nl_edges ,&
                  nrhs,&
                  1.0,&
                  S_rr,&
                  mlbddc%nl_edges  ,&
                  lambda_r,&
                  mlbddc%nl_edges  )

             ! Bck. substitution with U
             call DTRSM ( 'L', &
                  'U', &
                  'N', & ! Do not transpose
                  'N', & ! U is not unit diagonal
                  mlbddc%nl_edges  ,&
                  nrhs,&
                  1.0,&
                  S_rr,&
                  mlbddc%nl_edges  ,&
                  lambda_r,&
                  mlbddc%nl_edges  )

             call memfree ( aux , 'solve_edge_lagrange_multipliers_explicit_schur::aux' )
          end if
       end if
    end if
#else
      write (0,*) 'Error: par_preconditioner_dd_mlevel_bddc was not compiled with -DENABLE_BLAS.'
      write (0,*) 'Error: You must activate this cpp macro in order to use the BLAS'
      check(1==0)    
#endif 
    end subroutine solve_edge_lagrange_multipliers_explicit_schur
 


 ! Computes the harmonic extensions, i.e., p_prec_dd_mlevel_bddc%rPhi(:,:)
  ! Given the following partition of harm into corners & edges dofs:       
  ! harm = [rPhi_c]
  !        [rPhi_r],
  ! rPhi_c is given by [I_c 0]
  ! and rPhi_r by the solution of the following linear system
  ! A_rr rPhi_r = - ([A_rc 0] + C_r^T lamba_r) 
  subroutine compute_harmonic_extensions_corner_edges_partitioning (mlbddc, A_cr, A_rc, A_cc, A_rr, M_rr, A_rr_inv_C_r_T, rPhi, lambda_r, work, system, iparm, msglvl)
    implicit none
    type(par_preconditioner_dd_mlevel_bddc_t), intent(inout), target :: mlbddc
    real(rp)                 , intent(in)   :: A_cr( mlbddc%nl_corners, &
                                                     mlbddc%A_rr_gr%nv  )

    real(rp)                 , intent(in)   :: A_rc( mlbddc%A_rr_gr%nv, &
                                                     mlbddc%nl_corners  )

    real(rp)                 , intent(in)   :: A_cc( mlbddc%nl_corners  , &
                                                     mlbddc%nl_corners  )

    type(matrix_t)         , intent(inout)   :: A_rr
    type(preconditioner_t)        , intent(inout)   :: M_rr

    real (rp)                , intent(in)      :: A_rr_inv_C_r_T (:,:) 
    real(rp)                 , intent(inout)   :: rPhi (mlbddc%p_mat%dof_dist%nl,&
                                                        mlbddc%nl_edges+mlbddc%nl_corners)
    real(rp)                 , intent(in)   :: lambda_r ( mlbddc%nl_edges  , & 
                                                          (mlbddc%nl_edges+mlbddc%nl_corners) )

    real(rp)                 , intent(in)  :: work ( mlbddc%A_rr_gr%nv, & 
                                                     mlbddc%nl_corners   )

    character(len=1)         , intent(in)  :: system ! 'N' for rPhi computation, 'T' for lPhi computation

    integer                  , intent(in), target, optional :: iparm(64)
    integer                  , intent(in), optional         :: msglvl

    ! Locals
    real(rp), allocatable :: work1(:,:), work2(:,:), work3(:,:), aux(:)
    integer(ip)           :: j, k, base, iadj, jdof

    type(post_file_t)      :: lupos
    integer              :: me, np
    ! type(solver_control_t)      :: spars

    ! harm(:, 1:(mlbddc%nl_edges+mlbddc%nl_corners) ) = 0.0
    ! Set harm_c
    ! do i=1, mlbddc%nl_corners 
    !    rPhi(i, i) = 1.0_rp
    ! end do

   if ( mlbddc%nl_coarse > 0 ) then 

      call memalloc ( mlbddc%A_rr_gr%nv, (mlbddc%nl_corners+mlbddc%nl_edges) ,                       work2, __FILE__,__LINE__ )

      if ( mlbddc%nn_sys_sol_strat == corners_rest_part_solve_expl_schur ) then

         if ( mlbddc%kind_coarse_dofs == corners_and_edges .or. & 
              mlbddc%kind_coarse_dofs == corners_edges_and_faces .or. &
              mlbddc%kind_coarse_dofs == edges .or. & 
              mlbddc%kind_coarse_dofs == edges_and_faces .or. & 
              mlbddc%kind_coarse_dofs == faces .or. & 
              mlbddc%kind_coarse_dofs == corners_and_faces ) then

            work2 = 0.0_rp
            work2 ( :, 1:mlbddc%nl_corners ) = -work(:,:)

#ifdef ENABLE_BLAS
            if ( mlbddc%nl_edges > 0 ) then
               call DGEMM( 'N', &
                    'N', &
                    mlbddc%A_rr_gr%nv, &
                    (mlbddc%nl_corners+mlbddc%nl_edges) , &
                    mlbddc%nl_edges , &
                    -1.0, &
                    A_rr_inv_C_r_T, &
                    mlbddc%A_rr_gr%nv, &
                    lambda_r, &
                    mlbddc%nl_edges , &
                    1.0, &
                    work2, &
                    mlbddc%A_rr_gr%nv)
             ! write (*,*) a_ci ! DBG:
            end if
#else
            write (0,*) 'Error: par_preconditioner_dd_mlevel_bddc was not compiled with -DENABLE_BLAS.'
            write (0,*) 'Error: You must activate this cpp macro in order to use the BLAS'
            check(.false.)    
#endif
         else if (mlbddc%kind_coarse_dofs == corners) then
            call memalloc ( mlbddc%A_rr_gr%nv, (mlbddc%nl_corners+mlbddc%nl_edges), work1, __FILE__,__LINE__ )

            work1 = 0.0_rp
            work1 ( :, 1:mlbddc%nl_corners ) = -A_rc(:,:)

            mlbddc%spars_harm%nrhs=(mlbddc%nl_edges+mlbddc%nl_corners) 
            ! mlbddc%spars%method=direct
            work2 = 0.0_rp

            if ( mlbddc%internal_problems == handled_by_bddc_module) then
               call solve( A_rr, M_rr, &
                           work1, mlbddc%A_rr_gr%nv, &
                           work2, mlbddc%A_rr_gr%nv, &
                           mlbddc%spars_harm)
            else
               check(.false.)
!!$               call solve( mlbddc%A_rr_op, mlbddc%M_inv_op_rr, &
!!$                           work1, mlbddc%A_rr_gr%nv, &
!!$                           work2, mlbddc%A_rr_gr%nv, &
!!$                           mlbddc%spars_harm)
            end if
            call memfree ( work1,__FILE__,__LINE__) 
         end if
      end if
  
      if ( system == 'N') then 
         if ( mlbddc%subd_elmat_calc == phit_minus_c_i_t_lambda ) then
            if ( mlbddc%nl_corners > 0 ) then
               ! Compute block lagrange multipliers
               mlbddc%blk_lag_mul(1:mlbddc%nl_corners ,& 
                    1:mlbddc%nl_corners ) = -A_cc
               
               
               mlbddc%blk_lag_mul(1:mlbddc%nl_corners ,& 
                    mlbddc%nl_corners +1:) = 0.0_rp
               
#ifdef ENABLE_BLAS
               call DGEMM( 'N', &
                    'N', &
                    mlbddc%nl_corners  , &
                    mlbddc%nl_coarse   , &
                    mlbddc%A_rr_gr%nv, &
                    -1.0, &
                    A_cr, &
                    mlbddc%nl_corners, &
                    work2, &
                    mlbddc%A_rr_gr%nv, &
                    1.0, &
                    mlbddc%blk_lag_mul, &
                    mlbddc%nl_coarse  )
#else
               write (0,*) 'Error: par_preconditioner_dd_mlevel_bddc was not compiled with -DENABLE_BLAS.'
               write (0,*) 'Error: You must activate this cpp macro in order to use the BLAS'
               check(.false.)    
#endif
               
            end if
            mlbddc%blk_lag_mul(mlbddc%nl_corners +1:,:) = lambda_r
         end if
      end if

      base = mlbddc%nl_corners   
      do j=1, mlbddc%p_mat%dof_dist%nl  
         k =  mlbddc%iperm(j)
         ! write (*,*) j, k, mlbddc%p_mat%dof_dist%nl  , base DBG:
         if (j <= base) then
            ! harm_c [I_c 0] 
            rPhi(k,:) = 0.0_rp
            rPhi(k,j) = 1.0_rp
         else 
            ! rPhi_r
            rPhi(k,:) =  work2 (j-base,:)           
         end if
      end do

      if ( debug_verbose_level_3 ) then

            call memalloc ( mlbddc%p_mat%dof_dist%nl, aux, __FILE__,__LINE__ )

            call par_context_info ( mlbddc%p_mat%p_env%p_context, me, np )

            do iadj=1, mlbddc%nl_coarse 

               do j=1,mlbddc%p_mat%dof_dist%nl 
                  aux(j) = rPhi(j,iadj)
               end do

               if (system == 'N' ) then 
                  call postpro_open_file(1, trim('rPhi_mlevel_bddc_' // trim(ch(iadj)) // trim('.') // trim(ch(me+1))), lupos)
               elseif (system == 'T') then 
                  call postpro_open_file(1, trim('lPhi_mlevel_bddc_' // trim(ch(iadj)) // trim('.') // trim(ch(me+1))), lupos)
               end if
               call postpro(lupos,aux,'UNKNO',1,1.0) ! CDR
               call postpro_close_file(lupos)
            end do

            call memfree ( aux,__FILE__,__LINE__) 


      end if

      call memfree ( work2,__FILE__,__LINE__) 

   end if
  end subroutine compute_harmonic_extensions_corner_edges_partitioning


  ! Computes the harmonic extensions, i.e., mlbddc%rPhi(:,:)
  subroutine compute_harmonic_extensions_with_constraints (mlbddc, A_rr, M_rr, rPhi, system, iparm, msglvl)
    implicit none
    type(par_preconditioner_dd_mlevel_bddc_t), intent(inout), target        :: mlbddc
    integer                  , intent(in), target, optional :: iparm(64)
    integer                  , intent(in), optional         :: msglvl

    type(matrix_t)      , intent(inout)   :: A_rr
    type(preconditioner_t)     , intent(inout)   :: M_rr
    real(rp)              , intent(inout)   :: rPhi( mlbddc%p_mat%dof_dist%nl, &
                                                     mlbddc%nl_edges+mlbddc%nl_corners)    
    character(len=1)      , intent(in)      :: system ! Information about the regular system or the system assoc                                                        iated to the transpose matrix

    ! Locals
    real(rp), allocatable :: work1(:,:), work2(:,:), work3(:,:)
    integer(ip)           :: i, j, k, base, iadj, jdof, icoarse

    type(post_file_t)      :: lupos
    integer              :: me, np
    ! type(solver_control_t)      :: spars
    ! MATRIX MARKET          ! Write files with specified matrices
  integer(ip)              :: lunou    ! Creation of the file
  character(len=256)       :: filename ! Name of the file
    ! rPhi(:, 1:(mlbddc%nl_edges+mlbddc%nl_corners) ) = 0.0
    ! Set rPhi_c
    ! do i=1, mlbddc%nl_corners 
    !    rPhi(i, i) = 1.0_rp
    ! end do

    if ( mlbddc%nl_corners+mlbddc%nl_edges > 0 ) then
       call memalloc ( mlbddc%A_rr_gr%nv, (mlbddc%nl_corners+mlbddc%nl_edges) ,             work1, __FILE__,__LINE__ )

       call memalloc ( mlbddc%A_rr_gr%nv, (mlbddc%nl_corners+mlbddc%nl_edges) ,             work2, __FILE__,__LINE__ )



       work1 = 0.0_rp
       j = 1
       do i = mlbddc%p_mat%dof_dist%nl   + 1, mlbddc%A_rr_gr%nv
         work1 (i,j) = 1.0_rp
         j = j + 1 
       end do
      
       mlbddc%spars_harm%nrhs=(mlbddc%nl_corners+mlbddc%nl_edges) 
       ! mlbddc%spars%method=direct
       work2 = 0.0_rp

       if ( mlbddc%internal_problems == handled_by_bddc_module) then
          call solve( A_rr, M_rr, &
                      work1, mlbddc%A_rr_gr%nv, &
                      work2,  mlbddc%A_rr_gr%nv, &
                      mlbddc%spars_harm)
       else 
          check(.false.)
!!$          call solve( mlbddc%A_rr_op, mlbddc%M_inv_op_rr, &
!!$                      work1, mlbddc%A_rr_gr%nv, &
!!$                      work2,  mlbddc%A_rr_gr%nv, &
!!$                      mlbddc%spars_harm)
       end if

             rPhi = work2 (1:mlbddc%p_mat%dof_dist%nl ,:) 

             if ( system == 'N' ) then 
                if ( mlbddc%subd_elmat_calc == phit_minus_c_i_t_lambda ) then
                   mlbddc%blk_lag_mul = work2 (mlbddc%p_mat%dof_dist%nl +1:,:)
                end if
             end if

       if ( debug_verbose_level_3 ) then
          call memalloc ( 1, mlbddc%p_mat%dof_dist%nl, work3, __FILE__,__LINE__ )
          call par_context_info ( mlbddc%p_mat%p_env%p_context, me, np )


          do iadj=1, mlbddc%nl_coarse 
             do j=1,mlbddc%p_mat%dof_dist%nl 
                work3(1,j) = rPhi(j,iadj)
             end do

if (system == 'N') then 
    call postpro_open_file(1, trim('rPhi_mlevel_bddc_' // trim(ch(iadj)) // trim('.') // trim(ch(me+1))), lupos)
elseif (system == 'T') then
    call postpro_open_file(1, trim('lPhi_mlevel_bddc_' // trim(ch(iadj)) // trim('.') // trim(ch(me+1))), lupos)
end if
             call postpro(lupos,work3,'UNKNO',1,1.0) ! CDR
             call postpro_close_file(lupos)
          end do

          call memfree ( work3,__FILE__,__LINE__) 
       end if

       call memfree ( work1,__FILE__,__LINE__) 
       call memfree ( work2,__FILE__,__LINE__) 
    end if
  end subroutine compute_harmonic_extensions_with_constraints



  subroutine assemble_A_c (p_mat, mlbddc, realloc_harm_extensions)
    implicit none
    ! Parameters  
    type(par_matrix_t)         , intent(in)     :: p_mat
    type(par_preconditioner_dd_mlevel_bddc_t), intent(inout)  :: mlbddc
    logical, intent(in), optional :: realloc_harm_extensions

    ! Locals
    real(rp), allocatable :: a_ci (:,:), a_ci_gathered (:)
    real(rp), allocatable :: work (:,:)
    logical               :: realloc_harm_extensions_
    integer               :: info, iam
    integer (ip)          :: i, sum, j
    logical               :: i_am_fine_task, i_am_coarse_task, i_am_higher_level_task


    if ( present(realloc_harm_extensions) ) then
       realloc_harm_extensions_ = realloc_harm_extensions
    else
       realloc_harm_extensions_ = .true.
    end if

    ! Which duties I have?
    i_am_fine_task = (mlbddc%p_mat%p_env%p_context%iam >= 0)
    i_am_coarse_task = (mlbddc%c_context%iam >= 0)
    i_am_higher_level_task = (mlbddc%d_context%iam >= 0)
    if(i_am_fine_task) iam = mlbddc%p_mat%p_env%p_context%iam
    if(i_am_coarse_task) iam = mlbddc%c_context%iam
    if(i_am_higher_level_task) iam = mlbddc%d_context%iam

    if( i_am_coarse_task ) then
       sum = 0
       do i=1, mlbddc%g_context%np
          sum = sum + (mlbddc%ptr_coarse_dofs(i+1)-mlbddc%ptr_coarse_dofs(i)) * & 
               & (mlbddc%ptr_coarse_dofs(i+1)-mlbddc%ptr_coarse_dofs(i))  
       end do
       call memalloc ( sum, a_ci_gathered, __FILE__,__LINE__ )
    end if
 
    if ( i_am_fine_task ) then
       ! BEGIN FINE-GRID RELATED DUTIES
       call memalloc ( mlbddc%nl_coarse, mlbddc%nl_coarse, a_ci,__FILE__,__LINE__ )
      if ( mlbddc%nl_coarse > 0 ) then

          if ( mlbddc%subd_elmat_calc == phit_minus_c_i_t_lambda ) then

             call memalloc ( mlbddc%p_mat%dof_dist%nb, (mlbddc%nl_edges+mlbddc%nl_corners), work, __FILE__,__LINE__ )

             work = 0.0_rp

             ! Compute S \PHI_G == -C_i^t * blk_lag_mul (see art002)

             call update_with_c_i_trans ( mlbddc%p_mat%dof_dist%nl, &
                  mlbddc%nl_corners, &
                  mlbddc%nl_edges, &
                  mlbddc%coarse_dofs, &
                  mlbddc%p_mat%dof_dist%max_nparts, &
                  mlbddc%p_mat%dof_dist%omap%nl, &
                  mlbddc%p_mat%dof_dist%lobjs, &
                  mlbddc%p_mat%dof_dist%nb  , & 
                  (mlbddc%nl_edges+mlbddc%nl_corners) , & 
                  (mlbddc%nl_edges+mlbddc%nl_corners) , &
                  mlbddc%blk_lag_mul, & 
                  0.0_rp, &
                  -1.0_rp, & 
                  work, & 
                  mlbddc%C_weights)


             a_ci = 0.0_rp 

#ifdef ENABLE_BLAS
             call DGEMM( 'T', &
                  'N', &
                  mlbddc%nl_coarse  , &
                  mlbddc%nl_coarse  , &
                  mlbddc%p_mat%dof_dist%nb  , &
                  1.0, &
                  mlbddc%lPhi(mlbddc%p_mat%dof_dist%ni  +1,1), &
                  mlbddc%p_mat%dof_dist%nl   , &
                  work, &
                  mlbddc%p_mat%dof_dist%nb  , &
                  0.0, &
                  a_ci, &
                  mlbddc%nl_coarse  )

#else
             write (0,*) 'Error: par_preconditioner_dd_mlevel_bddc was not compiled with -DENABLE_BLAS.'
             write (0,*) 'Error: You must activate this cpp macro in order to use the BLAS'
             check(.false.)    
#endif

             call memfree ( work, 'assemble_a_c::work' )
          else

             ! HERE a_ci = \Phi^t A_i \Phi 
             call memalloc ( mlbddc%p_mat%dof_dist%nl, &
                  (mlbddc%nl_edges+mlbddc%nl_corners),  &
                  work, 'assemble_a_c::work' )

             work = 0.0_rp
#ifdef ENABLE_MKL
             ! work2 <- 0.0 * work2 + 1.0 * A * work1
             if (p_mat%f_matrix%symm == symm_true) then
                call mkl_dcsrmm ('N', &
                     p_mat%f_matrix%gr%nv, &
                     (mlbddc%nl_edges+mlbddc%nl_corners), &
                     p_mat%f_matrix%gr%nv, &
                     1.0, &
                     'SUNF', &
                     p_mat%f_matrix%a, &
                     p_mat%f_matrix%gr%ja, &
                     p_mat%f_matrix%gr%ia(1), &
                     p_mat%f_matrix%gr%ia(2), &
                     mlbddc%rPhi, &
                     p_mat%f_matrix%gr%nv, &
                     0.0, &
                     work, &
                     p_mat%f_matrix%gr%nv)

             else if (p_mat%f_matrix%symm == symm_false) then
                call mkl_dcsrmm ( 'N', &
                     p_mat%f_matrix%gr%nv, &
                     (mlbddc%nl_edges+mlbddc%nl_corners), &
                     p_mat%f_matrix%gr%nv, &
                     1.0, &
                     'GXXF', &
                     p_mat%f_matrix%a, &
                     p_mat%f_matrix%gr%ja, &
                     p_mat%f_matrix%gr%ia(1), &
                     p_mat%f_matrix%gr%ia(2), &
                     mlbddc%rPhi, &
                     p_mat%f_matrix%gr%nv, &
                     0.0, &
                     work, &
                     p_mat%f_matrix%gr%nv)
             end if
             ! write (*,*) work ! DBG:
#else
             call matmat ( p_mat%f_matrix, & 
                  (mlbddc%nl_edges+mlbddc%nl_corners), &
                  mlbddc%p_mat%dof_dist%nl, &
                  mlbddc%rPhi, &
                  mlbddc%p_mat%dof_dist%nl, &
                  work )
#endif

             a_ci = 0.0_rp 

#ifdef ENABLE_BLAS
             call DGEMM( 'T', &
                  'N', &
                  mlbddc%nl_coarse, &
                  mlbddc%nl_coarse, &
                  mlbddc%p_mat%dof_dist%nl, &
                  1.0, &
                  mlbddc%lPhi, &
                  mlbddc%p_mat%dof_dist%nl , &
                  work, &
                  mlbddc%p_mat%dof_dist%nl, &
                  0.0, &
                  a_ci, &
                  mlbddc%nl_coarse)

#else
             write (0,*) 'Error: par_preconditioner_dd_mlevel_bddc was not compiled with -DENABLE_BLAS.'
             write (0,*) 'Error: You must activate this cpp macro in order to use the BLAS'
             check(.false.)    
#endif
             call memfree ( work,__FILE__,__LINE__)

          end if
       end if

       if (realloc_harm_extensions_) then

          ! if ( mlbddc%unknowns == interface_unknowns ) then
          ! Re-alloc right side harmonic extensions (rPhi)
          call memalloc ( mlbddc%p_mat%dof_dist%nb, mlbddc%nl_edges+mlbddc%nl_corners,  work, __FILE__,__LINE__ )

          work = mlbddc%rPhi( &
               mlbddc%p_mat%dof_dist%ni+1:,:)

          call memfree ( mlbddc%rPhi,__FILE__,__LINE__)

          call memalloc ( mlbddc%p_mat%dof_dist%nb, mlbddc%nl_edges+mlbddc%nl_corners,  mlbddc%rPhi, __FILE__,__LINE__ )

          mlbddc%rPhi = work
          call memfree ( work,__FILE__,__LINE__)

          ! Re-alloc left side harmonic extensions (lPhi)
          call memalloc ( mlbddc%p_mat%dof_dist%nb, mlbddc%nl_edges+mlbddc%nl_corners,  work, __FILE__,__LINE__ )

          work = mlbddc%lPhi( &
               mlbddc%p_mat%dof_dist%ni+1:,:)

          call memfree ( mlbddc%lPhi,__FILE__,__LINE__)

          call memalloc ( mlbddc%p_mat%dof_dist%nb, mlbddc%nl_edges+mlbddc%nl_corners,  mlbddc%lPhi, __FILE__,__LINE__ )

          mlbddc%lPhi = work
          call memfree ( work,__FILE__,__LINE__)

          ! end if
       end if
    end if


    if ( temp_fine_coarse_grid_overlap ) then 
       if (.not. i_am_higher_level_task) then
          if (i_am_fine_task) then
             call par_timer_stop ( mlbddc%timer_assprec_ov_fine  ) 
             if ( mlbddc%co_sys_sol_strat == serial_gather ) then
                call par_timer_stop ( mlbddc%timer_fillprec_ov_coarse_header )
             end if
          else if ( i_am_coarse_task ) then
             call par_timer_stop ( mlbddc%timer_assprec_ov_coarse ) 
          end if
          
          call par_timer_start ( mlbddc%timer_coll_fillprec  )
       end if
    end if

    if ( i_am_coarse_task ) then
       call rcv_coarse_stiffness_matrices ( mlbddc%g_context%icontxt, &
                                            mlbddc%g_context%np,      &
                                            mlbddc%pad_collectives,          &
                                            mlbddc%max_coarse_dofs,          &
                                            mlbddc%ptr_coarse_dofs,          &
                                            sum, &
                                            a_ci_gathered ) 

       if ( temp_fine_coarse_grid_overlap ) then
          call par_timer_stop ( mlbddc%timer_coll_fillprec  )
          call psb_barrier ( mlbddc%g_context%icontxt )
          call par_timer_start ( mlbddc%timer_fillprec_ov_coarse )
          if ( mlbddc%co_sys_sol_strat == recursive_bddc ) then
             call par_timer_start ( mlbddc%p_M_c%timer_fillprec_ov_coarse_header )
          end if
       end if

       if (  mlbddc%co_sys_sol_strat == serial_gather ) then  
          call matrix_fill_val (mlbddc%A_c)

          call sum_coarse_stiffness_matrices ( mlbddc%A_c, &
                                               mlbddc%g_context%np, & 
                                               mlbddc%ptr_coarse_dofs, &
                                               mlbddc%f_mesh_c%lnods, &
                                               sum, &
                                               a_ci_gathered)
          
       else if ( mlbddc%co_sys_sol_strat == recursive_bddc ) then
          assert(mlbddc%p_mat%p_env%num_levels>2) ! Assuming last level direct
          
          call matrix_fill_val (mlbddc%p_mat_c%f_matrix)

          call sum_coarse_stiffness_matrices ( mlbddc%p_mat_c%f_matrix, &
                                               mlbddc%g_context%np, & 
                                               mlbddc%ptr_coarse_dofs, &
                                               mlbddc%p_mesh_c%f_mesh%lnods, &
                                               sum, &
                                               a_ci_gathered, &
                                               mlbddc%erenumbering_c%iperm)
       end if
       call memfree ( a_ci_gathered,__FILE__,__LINE__)


    else if ( i_am_fine_task ) then
       call snd_coarse_stiffness_matrices ( mlbddc%g_context%icontxt, &
                                            mlbddc%g_context%np,      &
                                            mlbddc%pad_collectives,          &
                                            mlbddc%nl_coarse,                &
                                            mlbddc%max_coarse_dofs,          &
                                            a_ci )
       call memfree ( a_ci,__FILE__,__LINE__)     

       if ( temp_fine_coarse_grid_overlap ) then
          call par_timer_stop ( mlbddc%timer_coll_fillprec  )
          call psb_barrier ( mlbddc%g_context%icontxt )
          call par_timer_start ( mlbddc%timer_fillprec_ov_fine  )
       end if

    end if

  end subroutine assemble_A_c

  subroutine sum_coarse_stiffness_matrices  ( A_c,              &
                                              np,               &
                                              ptr_coarse_dofs,  &
                                              lst_coarse_dofs,  & 
                                              sz_a_ci_gathered, &
                                              a_ci_gathered,    &
                                              perm )
    implicit none

    ! Parameters
    type(matrix_t), intent(inout) :: A_c
    integer(ip)     , intent(in)    :: np, sz_a_ci_gathered
    integer(ip)     , intent(in)    :: ptr_coarse_dofs (np+1)
    integer(ip)     , intent(in)    :: lst_coarse_dofs (ptr_coarse_dofs(np+1))
    real   (rp)     , intent(in)    :: a_ci_gathered (sz_a_ci_gathered)
    integer(ip)     , optional, intent(in)  :: perm(np-1) ! Element re-numbering

    ! Locals
    real(rp), allocatable  :: elm(:,:)
    integer (ip) :: i, nnode, idof, jdof, inode, jnode
    integer (ip) :: base, base_old, base_new
    integer      :: lunou
    integer (ip), allocatable :: base_displs(:)


    ! write (*,*) 'ZZZZZZ', a_ci_gathered

    if ( .not. present(perm) ) then
       base = 0 
       do i=1, np-1
          nnode = ptr_coarse_dofs(i+1) - ptr_coarse_dofs(i)

          call memalloc (  nnode, nnode, elm, __FILE__, __LINE__ )

          ! base + (j-1)*lda + i       
          do jnode=1, nnode
             do inode=1, nnode
                elm(inode, jnode) = a_ci_gathered ( base + (inode-1) + (jnode-1)*nnode + 1)
             end do
          end do

          call matrix_assembly ( nnode, &
                                     lst_coarse_dofs(ptr_coarse_dofs(i)+1:ptr_coarse_dofs(i+1)), &
                                     elm, &
                                     A_c ) 

          call memfree ( elm, __FILE__, __LINE__ )

          base = base + nnode*nnode
       end do

    else
       call memalloc ( np, base_displs, __FILE__, __LINE__  )

       base_displs(1) = 0
       do i=1, np-1
          base_displs(i+1) = base_displs(i)   + (ptr_coarse_dofs(i+1) - ptr_coarse_dofs(i)) * & 
                                                (ptr_coarse_dofs(i+1) - ptr_coarse_dofs(i))
       end do

       base_new = 0

       do i=1, np-1
          nnode = ptr_coarse_dofs(perm(i)+1) - ptr_coarse_dofs(perm(i))

          base_old  = base_displs(perm(i))

          call memalloc (  nnode, nnode, elm, __FILE__, __LINE__ )


          ! base + (j-1)*lda + i       
          do jnode=1, nnode
             do inode=1, nnode                                 
                elm(inode, jnode) = a_ci_gathered ( base_old + (inode-1) + (jnode-1)*nnode + 1 )
             end do
          end do
                    
          call matrix_assembly ( nnode, &
                                     lst_coarse_dofs(base_new+1:base_new+nnode), &
                                     elm, &
                                     A_c ) 
           
          call memfree ( elm, __FILE__, __LINE__ )
                  
          base_new = base_new + nnode
          
       end do

       call memfree ( base_displs, __FILE__, __LINE__  )
       
    end if

    if (debug_verbose_level_2) then 
       call matrix_print ( 6, A_c ) ! DBG:
    end if


  end subroutine sum_coarse_stiffness_matrices


  subroutine snd_coarse_stiffness_matrices ( icontxt_g, np_g, & 
                                             pad_collectives, & 
                                             nl_coarse,  &
                                             max_coarse_dofs, &
                                             a_ci )
use psb_const_mod_names
use psb_penv_mod_names

#ifdef MPI_MOD
use mpi
#endif

    implicit none

#ifdef MPI_H
    include 'mpif.h'
#endif

    ! Parameters
    integer    , intent(in)  :: icontxt_g, np_g
    integer(ip), intent(in)  :: pad_collectives, nl_coarse
    integer(ip), intent(in)  :: max_coarse_dofs
    real   (rp), intent(in)  :: a_ci (nl_coarse,nl_coarse)

    ! Locals
    integer                  :: mpi_comm_g, iret, root_g
    integer(ip)              :: j, off
    real   (rp), allocatable :: pad_buf_snd (:)
    real(rp)                 :: rcv_dumvec(1), rcv_dum
    integer                  :: int_rcv_dumvec(1), int_rcv_dum

    ! Get MPI communicator associated to icontxt (in
    ! the current implementation of our wrappers
    ! to the MPI library icontxt and mpi_comm are actually 
    ! the same)
    call psb_get_mpicomm (icontxt_g, mpi_comm_g)
    root_g = np_g -1

    if ( pad_collectives == pad ) then 
       call memalloc ( max_coarse_dofs*max_coarse_dofs, pad_buf_snd, __FILE__,__LINE__ )

       ! Why this is needed ? Commented it out 
       ! pad_buf_snd = -1.0 

       ! Pack a_ci into send buffer
       off = 0
       do j = 1, nl_coarse 
          pad_buf_snd ( off+1:off+nl_coarse ) = a_ci(:,j)
          off = off + max_coarse_dofs
       end do

       call mpi_gather ( pad_buf_snd, max_coarse_dofs*max_coarse_dofs, psb_mpi_real, &
                         rcv_dum, int_rcv_dum, psb_mpi_real, &
                         root_g, mpi_comm_g, iret)
       check ( iret == mpi_success )

       call memfree ( pad_buf_snd,__FILE__,__LINE__)

    else if ( pad_collectives == nopad ) then
       call mpi_gatherv( a_ci,  nl_coarse * nl_coarse, psb_mpi_real, &
                         rcv_dum, int_rcv_dumvec, int_rcv_dumvec, psb_mpi_real, &
                         root_g, mpi_comm_g, iret)
       check ( iret == mpi_success )
    end if

  end subroutine snd_coarse_stiffness_matrices


  subroutine rcv_coarse_stiffness_matrices ( icontxt_g, np_g, & 
                                             pad_collectives, & 
                                             max_coarse_dofs, &
                                             ptr_coarse_dofs, &
                                             sz_a_ci_gathered, &
                                             a_ci_gathered  )
use psb_const_mod_names
use psb_penv_mod_names

#ifdef MPI_MOD
use mpi
#endif

    implicit none

#ifdef MPI_H
    include 'mpif.h'
#endif

    ! Parameters
    integer    , intent(in)  :: icontxt_g,np_g
    integer(ip), intent(in)  :: pad_collectives
    integer(ip), intent(in)  :: ptr_coarse_dofs (np_g+1)
    integer(ip), intent(in)  :: max_coarse_dofs
    integer(ip), intent(in)  :: sz_a_ci_gathered
    real   (rp), intent(out) :: a_ci_gathered (*)

    ! Locals
    integer                  :: mpi_comm_g, iret, root_g
    integer(ip)              :: i,j, off, offrcv, ncoarse
    integer(ip), allocatable :: recvcounts (:), displs (:)
    real   (rp), allocatable :: pad_buf_snd (:), pad_buf_rcv(:)
    real (rp)                :: snd_dum

    ! Get MPI communicator associated to icontxt (in
    ! the current implementation of our wrappers
    ! to the MPI library icontxt and mpi_comm are actually 
    ! the same)
    call psb_get_mpicomm (icontxt_g, mpi_comm_g)
    root_g = np_g -1

    if ( pad_collectives == pad ) then 

       call memalloc ( max_coarse_dofs*max_coarse_dofs, pad_buf_snd, __FILE__,__LINE__ )
       call memalloc ( max_coarse_dofs*max_coarse_dofs*np_g, pad_buf_rcv, __FILE__,__LINE__ )

       call mpi_gather ( pad_buf_snd, max_coarse_dofs*max_coarse_dofs, psb_mpi_real, &
                         pad_buf_rcv, max_coarse_dofs*max_coarse_dofs, psb_mpi_real, &
                         root_g, mpi_comm_g, iret)
       check ( iret == mpi_success )
  
       ! Unpack padded data
       offrcv = 0
       off    = 0 
       do i=1, np_g
          ncoarse = ptr_coarse_dofs(i+1)-ptr_coarse_dofs(i)  
          do j = 1, ncoarse
             ! write (*,*) 'ZZZ', pad_buf_rcv ( offrcv+1:offrcv+ncoarse )
             a_ci_gathered ( off+1: off+ncoarse ) = pad_buf_rcv ( offrcv+1:offrcv+ncoarse ) 
             offrcv = offrcv + max_coarse_dofs
             off    = off + ncoarse
          end do
          offrcv = i*max_coarse_dofs*max_coarse_dofs
       end do

       call memfree ( pad_buf_rcv,__FILE__,__LINE__)
       call memfree ( pad_buf_snd, __FILE__,__LINE__ )

    else if ( pad_collectives == nopad ) then

       call memalloc ( np_g, recvcounts,  __FILE__,__LINE__ )
       call memalloc ( np_g+1, displs,  __FILE__,__LINE__ )

       displs(1) = 0 
       do i=1, np_g
          recvcounts(i) = (ptr_coarse_dofs(i+1)-ptr_coarse_dofs(i)) * (ptr_coarse_dofs(i+1)-ptr_coarse_dofs(i))
          displs(i+1)   = displs(i) + recvcounts(i)
       end do

       call mpi_gatherv( snd_dum,  0, psb_mpi_real, &
                         a_ci_gathered, recvcounts,  displs, psb_mpi_real, &
                         root_g, mpi_comm_g, iret)
       check ( iret == mpi_success )

       call memfree ( recvcounts,__FILE__,__LINE__)
       call memfree ( displs    ,__FILE__,__LINE__)

    end if

  end subroutine rcv_coarse_stiffness_matrices

  !=================================================================================================
  recursive subroutine par_preconditioner_dd_mlevel_bddc_apply_all_unk (mlbddc, x, y)
    implicit none
    ! Parameters
    type(par_preconditioner_dd_mlevel_bddc_t), intent(in) :: mlbddc 
    type(par_vector_t)                , intent(in)        :: x
    type(par_vector_t)                , intent(inout)     :: y

    ! Locals
    type(par_vector_t)      :: r, dx, v1, v2, v3, aux
    type(par_vector_t)      :: r_I, r_G, x_I, x_G, r_v1, y_I, y_G, dx_I, dx_G
    type(par_vector_t)      :: p_r_c, p_z_c
    type(vector_t)      :: r_c, z_c, dum_vec
    logical               :: i_am_coarse_task, i_am_fine_task, i_am_higher_level_task

    ! The routine requires the partition/context info
    assert ( associated(mlbddc%p_mat) )
    assert ( associated(mlbddc%p_mat%p_env) )
    assert ( mlbddc%p_mat%p_env%created )
    assert ( mlbddc%p_mat%p_env%num_levels > 1 ) 

    assert ( associated(mlbddc%p_mat%p_env%w_context) )
    assert ( mlbddc%p_mat%p_env%w_context%created .eqv. .true. )
    assert ( associated(mlbddc%p_mat%p_env%p_context) )
    assert ( mlbddc%p_mat%p_env%p_context%created .eqv. .true. )
    assert ( associated(mlbddc%p_mat%p_env%q_context) )
    assert ( mlbddc%p_mat%p_env%q_context%created .eqv. .true. )
    ! Check appropriate assignment to context: 1) I'm in w_context and 2) I'm in p_context OR in q_context but not in both
    assert ( mlbddc%p_mat%p_env%w_context%iam >= 0)
    assert ( (mlbddc%p_mat%p_env%p_context%iam >=0 .and. mlbddc%p_mat%p_env%q_context%iam < 0) .or. (mlbddc%p_mat%p_env%p_context%iam < 0 .and. mlbddc%p_mat%p_env%q_context%iam >= 0))
    assert ( mlbddc%c_context%created .eqv. .true. )
    assert ( mlbddc%d_context%created .eqv. .true. )
    assert ( mlbddc%g_context%created .eqv. .true. )
    assert ( associated(mlbddc%p_mat%p_env%b_context) )
    assert ( mlbddc%p_mat%p_env%b_context%created .eqv. .true. )

    ! Which duties I have?
    i_am_fine_task = (mlbddc%p_mat%p_env%p_context%iam >= 0)
    i_am_coarse_task = (mlbddc%c_context%iam >= 0)
    i_am_higher_level_task = (mlbddc%d_context%iam >= 0)

!!$    if ( temp_fine_coarse_grid_overlap ) then 
!!$       if (i_am_fine_task) then
!!$          call par_timer_stop ( mlbddc%timer_fillprec_ov_fine  ) 
!!$       else if ( i_am_coarse_task  ) then
!!$          call par_timer_stop ( mlbddc%timer_fillprec_ov_coarse ) 
!!$       end if
!!$    end if

    assert ( mlbddc%correction_mode == additive .or. mlbddc%correction_mode == additive_symmetric )

    if ( i_am_fine_task ) then

       if ( mlbddc%correction_mode == additive ) then

          ! Create temporary vector
          call par_vector_alloc ( mlbddc%p_mat%dof_dist, mlbddc%p_mat%p_env, r )

          ! par_vector_create view calls next
          ! requires a state to be assigned to r !!!
          r%state = x%state
          call par_vector_create_view ( r,                                      1, mlbddc%p_mat%dof_dist%ni, r_I )
          call par_vector_create_view ( r, mlbddc%p_mat%dof_dist%ni+1, mlbddc%p_mat%dof_dist%nl, r_G )
          call par_vector_create_view ( x, mlbddc%p_mat%dof_dist%ni+1, mlbddc%p_mat%dof_dist%nl, x_G )

          ! Init interior vertices to zero
          call vector_zero ( r_I%f_vector )

          ! Phase 1: compute coarse-grid correction v1
          call par_vector_clone ( x_G, v1 )

          ! Copy residual into temporal vector
          call par_vector_copy  ( x_G, r_G  )

          if ( r%state == part_summed ) then 
             ! Summ  residual
             call par_vector_comm ( r_G )
          end if

          ! 1.2. ws_vec_G <- W_i R_i r
          ! Weight residual
          if ( associated (mlbddc%weight) ) then
             call par_vector_weight ( r_G, mlbddc%weight ) 
          else
             call par_vector_weight ( r_G ) 
          end if

          call par_preconditioner_dd_mlevel_bddc_compute_c_g_corr_ass_r_c ( mlbddc, r_G, dum_vec )

          call par_vector_clone ( x_G, v2 )

          call par_preconditioner_dd_mlevel_bddc_compute_fine_grid_correction ( mlbddc, r, v2 )

!!$          if ( temp_fine_coarse_grid_overlap ) then 
!!$             call par_timer_stop ( mlbddc%timer_applyprec_ov_fine  ) 
!!$          end if

          call par_preconditioner_dd_mlevel_bddc_compute_c_g_corr_scatter ( mlbddc, dum_vec, v1 )

          call par_vector_create_view ( y, mlbddc%p_mat%dof_dist%ni+1, mlbddc%p_mat%dof_dist%nl, y_G )
          call par_vector_create_view ( x,                                 1, mlbddc%p_mat%dof_dist%ni, x_I )

          call par_vector_copy ( v1, y_G )
          call par_vector_pxpy ( v2, y_G )
          call par_vector_free ( v1 )
          call par_vector_free ( v2 )

          if ( associated (mlbddc%weight) ) then
             call par_vector_weight ( y_G, mlbddc%weight ) 
          else
             call par_vector_weight ( y_G ) 
          end if
          call par_vector_comm   ( y_G )

          ! r_I <- A_IG * y_G
          call operator_dd_apply_A_IG ( mlbddc%A_II_inv, y_G%f_vector, r_I%f_vector )

          ! r <- x_I - r_I <- -r_I + x_I
          call par_vector_pxmy  ( x_I, r_I )

          call par_vector_create_view (  y, 1, mlbddc%p_mat%dof_dist%ni, y_I  )
          call operator_dd_solve_A_II ( mlbddc%A_II_inv, r_I%f_vector, y_I%f_vector )

          !mlbddc%num_dirichlet_solves = mlbddc%num_dirichlet_solves + 1
          !mlbddc%dirichlet_its (mlbddc%num_dirichlet_solves) = mlbddc%spars_dirichlet%it

          call par_vector_free   (r)

       else if (mlbddc%correction_mode == additive_symmetric) then
          ! Create temporary vectors
          call par_vector_clone ( x, r  )
          call par_vector_clone ( x, aux  )
          call par_vector_clone ( y, dx )

          call par_vector_copy ( x, aux )
          if ( aux%state == full_summed ) then
             call par_vector_weight ( aux )
          end if

          ! par_vector_create view calls next
          ! requires a state to be assigned to r !!!
          r%state = x%state
          call par_vector_create_view ( r,                                      1, mlbddc%p_mat%dof_dist%ni, r_I )
          call par_vector_create_view ( r, mlbddc%p_mat%dof_dist%ni+1, mlbddc%p_mat%dof_dist%nl, r_G )

          call par_vector_create_view ( aux,                                      1, mlbddc%p_mat%dof_dist%ni, x_I )
          call par_vector_create_view ( aux, mlbddc%p_mat%dof_dist%ni+1, mlbddc%p_mat%dof_dist%nl, x_G )

          call par_vector_create_view ( y,                                      1, mlbddc%p_mat%dof_dist%ni, y_I )
          call par_vector_create_view ( y, mlbddc%p_mat%dof_dist%ni+1, mlbddc%p_mat%dof_dist%nl, y_G )

          call par_vector_create_view ( dx,                                      1, mlbddc%p_mat%dof_dist%ni, dx_I )
          call par_vector_create_view ( dx, mlbddc%p_mat%dof_dist%ni+1, mlbddc%p_mat%dof_dist%nl, dx_G )


          ! 1) Compute dx_I = A_II^-1 r_I,   dx_G = 0
          call operator_dd_solve_A_II ( mlbddc%A_II_inv, x_I%f_vector, dx_I%f_vector )
          call vector_zero            ( dx_G%f_vector )
          !mlbddc%num_dirichlet_solves = mlbddc%num_dirichlet_solves + 1
          !mlbddc%dirichlet_its (mlbddc%num_dirichlet_solves) = mlbddc%spars_dirichlet%it

          ! 2) Compute x = x + dx     (idem x_I = x_I + dx_I,    x_G = x_G)
          call par_vector_copy ( dx, y )

          ! 3) Compute r = r - A dx
          ! r <- Adx
          call par_matvec ( mlbddc%p_mat, dx, r ) 
          ! r <- -r + x
          call par_vector_pxmy  ( aux, r )


          ! 4) Compute dx = A_{BDDC}^{-1} r

          ! Phase 1: compute coarse-grid correction v1
          call par_vector_clone ( x_G, v1 )

          if ( r%state == part_summed ) then 
             ! Summ  residual
             call par_vector_comm ( r_G )
          end if

          ! 1.2. ws_vec_G <- W_i R_i r
          ! Weight residual
          if ( associated (mlbddc%weight) ) then
             call par_vector_weight ( r_G, mlbddc%weight ) 
          else
             call par_vector_weight ( r_G ) 
          end if

          call par_preconditioner_dd_mlevel_bddc_compute_c_g_corr_ass_r_c ( mlbddc, r_G, dum_vec )

          call par_vector_clone ( x_G, v2 )

          call par_preconditioner_dd_mlevel_bddc_compute_fine_grid_correction ( mlbddc, r, v2 )

!!$          if ( temp_fine_coarse_grid_overlap ) then 
!!$             call par_timer_stop ( mlbddc%timer_applyprec_ov_fine  ) 
!!$          end if

          call par_preconditioner_dd_mlevel_bddc_compute_c_g_corr_scatter ( mlbddc, z_c, v1 )

          call par_vector_copy ( v1, dx_G )
          call par_vector_pxpy ( v2, dx_G )
          call par_vector_free ( v1 )
          call par_vector_free ( v2 )

          if ( associated (mlbddc%weight) ) then
             call par_vector_weight ( dx_G, mlbddc%weight ) 
          else
             call par_vector_weight ( dx_G ) 
          end if
          call par_vector_comm   ( dx_G )

          call par_vector_zero ( dx_I )

          ! 5) Do 2) and 3)
          ! 2) Compute x = x + dx
          call par_vector_pxpy ( dx, y )

          ! 3) Compute r = r - A dx 
          call par_vector_copy ( r, aux )
          ! r <- Adx
          call par_matvec ( mlbddc%p_mat, dx, r ) 
          ! r <- -r + x
          call par_vector_pxmy  ( aux, r )

          ! 6) Compute 1), 2)
          call operator_dd_solve_A_II ( mlbddc%A_II_inv, r_I%f_vector, dx_I%f_vector )
          call vector_zero            ( dx_G%f_vector )
          !mlbddc%num_dirichlet_solves = mlbddc%num_dirichlet_solves + 1
          !mlbddc%dirichlet_its (mlbddc%num_dirichlet_solves) = mlbddc%spars_dirichlet%it

          ! 2) Compute x = x + dx
          call par_vector_pxpy ( dx, y )

          call par_vector_free   (r)
          call par_vector_free   (aux)
          call par_vector_free   (dx)

          if ( temp_fine_coarse_grid_overlap ) then 
            if (i_am_fine_task) then
               call par_timer_stop ( mlbddc%timer_applyprec_ov_coarse_tail )
            end if
          end if
       end if
    else ! if ( i_am_coarse_task .or. i_am_higher_level_task ) then 
       if(mlbddc%co_sys_sol_strat == serial_gather) then
          if ( i_am_coarse_task ) then
             ! Assemble coarse-grid residual
             call vector_alloc ( mlbddc%A_c%gr%nv, r_c )    
             call par_preconditioner_dd_mlevel_bddc_compute_c_g_corr_ass_r_c ( mlbddc, r, r_c )

             ! Solve coarse-grid problem serially
             call vector_alloc ( mlbddc%A_c%gr%nv, z_c )    
             mlbddc%spars_coarse%nrhs=1
             call vector_zero(z_c)
             if ( mlbddc%internal_problems == handled_by_bddc_module) then
                call solve( mlbddc%A_c, mlbddc%M_c, r_c, z_c, mlbddc%spars_coarse)
             else
                check(.false.)
!!$             call solve( mlbddc%A_c_op, mlbddc%M_inv_op_c, r_c, z_c, mlbddc%spars_coarse)
             end if
             !mlbddc%num_coarse_solves = mlbddc%num_coarse_solves + 1
             !mlbddc%coarse_its (mlbddc%num_coarse_solves) = mlbddc%spars_coarse%it

             ! Scatter solution of coarse-grid problem 
             call par_preconditioner_dd_mlevel_bddc_compute_c_g_corr_scatter ( mlbddc, z_c, v1 )

             call vector_free (z_c)
             call vector_free (r_c)
          end if

       else if(mlbddc%co_sys_sol_strat == recursive_bddc) then
          assert(mlbddc%p_mat%p_env%num_levels>2)
          if ( i_am_coarse_task .or. i_am_higher_level_task ) then
             ! Assemble coarse-grid residual in parallel
             call par_vector_alloc ( mlbddc%dof_dist_c, mlbddc%p_env_c, p_r_c)
             p_r_c%state = part_summed
             call par_preconditioner_dd_mlevel_bddc_compute_c_g_corr_ass_r_c ( mlbddc, r, p_r_c%f_vector )
             ! call vector_print_matrix_market (6, p_r_c%f_vector)
             ! Solve coarse-grid problem via recursive bddc
             ! AFM: In the future, the following call should be replaced by a call that allows
             ! the solution of the coarse-grid problem via a krylov subspace solver  
             ! A simple solve specialization may help here ??? 
             call par_vector_alloc ( mlbddc%dof_dist_c, mlbddc%p_env_c, p_z_c)
             p_z_c%state = full_summed

             ! call par_preconditioner_dd_mlevel_bddc_apply_all_unk ( mlbddc%p_mat_c, mlbddc%p_M_c, p_r_c, p_z_c )
             mlbddc%spars_coarse%nrhs=1
             call abstract_solve( mlbddc%p_mat_c, mlbddc%p_M_c, p_r_c, p_z_c, mlbddc%spars_coarse, mlbddc%p_env_c)

             ! Scatter solution of coarse-grid problem 
             call par_preconditioner_dd_mlevel_bddc_compute_c_g_corr_scatter ( mlbddc, p_z_c%f_vector, v1 )

             call par_vector_free (p_z_c)
             call par_vector_free (p_r_c)
          end if
       end if
!!$       if ( temp_fine_coarse_grid_overlap ) then 
!!$          if ( i_am_coarse_task ) then
!!$             call par_timer_stop ( mlbddc%timer_applyprec_ov_coarse  ) 
!!$          end if
!!$       end if
    end if
    
  end subroutine par_preconditioner_dd_mlevel_bddc_apply_all_unk


  !=============================================================================
  subroutine par_preconditioner_dd_mlevel_bddc_static_condensation (p_mat, mlbddc, b, x)
    implicit none
    ! Parameters
    type(par_matrix_t)                , intent(in)    :: p_mat
    type(par_preconditioner_dd_mlevel_bddc_t), intent(inout) :: mlbddc
    type(par_vector_t)                , intent(in)    :: b
    type(par_vector_t)                , intent(inout) :: x

    ! Locals
    type(par_vector_t)      :: r
    type(par_vector_t)      :: r_I, b_I, x_G, x_I

    ! Pointer to part/context object is required
    assert ( associated(p_mat%p_env) )
    assert ( p_mat%p_env%created )
    assert ( p_mat%p_env%num_levels > 1 ) 

    if ( mlbddc%p_mat%p_env%p_context%iam >= 0 ) then
       call par_vector_clone ( b, r )

       call par_vector_create_view ( r, 1, mlbddc%p_mat%dof_dist%ni, r_I )
       call par_vector_create_view ( b, 1, mlbddc%p_mat%dof_dist%ni, b_I )
       call par_vector_create_view ( x, 1, mlbddc%p_mat%dof_dist%ni, x_I )
       call par_vector_create_view ( x, mlbddc%p_mat%dof_dist%ni+1, mlbddc%p_mat%dof_dist%nl, x_G )
       
       ! r_I <- A_IG * y_G
       call operator_dd_apply_A_IG ( mlbddc%A_II_inv, x_G%f_vector, r_I%f_vector )
       
       ! r <- x_I - r_I <- -r_I + x_I
       call par_vector_pxmy ( b_I, r_I )
       
       call operator_dd_solve_A_II ( mlbddc%A_II_inv, r_I%f_vector, x_I%f_vector )

       mlbddc%num_dirichlet_solves = mlbddc%num_dirichlet_solves + 1
       mlbddc%dirichlet_its (mlbddc%num_dirichlet_solves) = mlbddc%spars_dirichlet%it
       
       call par_vector_free (r)
    end if

  end subroutine par_preconditioner_dd_mlevel_bddc_static_condensation

  subroutine par_preconditioner_dd_mlevel_bddc_compute_c_g_corr_ass_r_c ( mlbddc, r_G, r_c )
   implicit none
   ! Parameters
   type(par_preconditioner_dd_mlevel_bddc_t) ,intent(in) :: mlbddc
   type(par_vector_t)                 ,intent(in) :: r_G
   type(vector_t)                 ,intent(inout) :: r_c 

   ! Locals
   integer :: me, np
   integer(ip) :: i, j
   real(rp), allocatable :: r_ci(:)
   logical :: i_am_fine_task, i_am_coarse_task, i_am_higher_level_task

   real(rp), allocatable :: r_ci_gathered(:)
   integer(ip) :: sum

   ! Which duties I have?
   i_am_fine_task = (mlbddc%p_mat%p_env%p_context%iam >= 0)
   i_am_coarse_task = (mlbddc%c_context%iam >= 0)
   i_am_higher_level_task = (mlbddc%d_context%iam >= 0)

   if ( i_am_fine_task ) then
     call memalloc ( mlbddc%nl_coarse, r_ci, __FILE__,__LINE__ )
     
     ! 1.3. rci <- rPhi^T * W_i R_i r | G
     call apply_harm_trans ( mlbddc, r_G, r_ci )
   end if
   
   if ( temp_fine_coarse_grid_overlap ) then 
      if ( .not. i_am_higher_level_task ) then
         if (i_am_fine_task) then
            call par_timer_stop ( mlbddc%timer_fillprec_ov_fine  ) 
            if ( mlbddc%co_sys_sol_strat == serial_gather ) then
               call par_timer_stop ( mlbddc%timer_applyprec_ov_coarse_header )
            end if
         else if ( i_am_coarse_task ) then
            call par_timer_stop ( mlbddc%timer_fillprec_ov_coarse ) 
         end if
         call par_timer_start ( mlbddc%timer_coll_applyprec  )
      end if
   end if
   
   ! 1.4. r_c <- comm r_ci
   if ( i_am_coarse_task ) then
      sum = 0
      do i=1, mlbddc%g_context%np
         sum = sum + (mlbddc%ptr_coarse_dofs(i+1)- &
              mlbddc%ptr_coarse_dofs(i))  
      end do
      call memalloc (  sum, r_ci_gathered, __FILE__,__LINE__ )

      call rcv_coarse_stiffness_vectors ( mlbddc%g_context%icontxt, &
           mlbddc%g_context%np,      &
           mlbddc%pad_collectives,          &
           mlbddc%max_coarse_dofs,          &
           mlbddc%ptr_coarse_dofs,          &
           sum,                                &      
           r_ci_gathered )

      if ( temp_fine_coarse_grid_overlap ) then
         call par_timer_stop ( mlbddc%timer_coll_applyprec  )
         call psb_barrier ( mlbddc%g_context%icontxt )
         call par_timer_start ( mlbddc%timer_applyprec_ov_coarse )
         if ( mlbddc%co_sys_sol_strat == recursive_bddc ) then
            call par_timer_start ( mlbddc%p_M_c%timer_applyprec_ov_coarse_header )
         end if
      end if

      if (  mlbddc%co_sys_sol_strat == serial_gather ) then  
         call sum_coarse_stiffness_vectors ( r_c, &
                                             mlbddc%g_context%np, & 
                                             mlbddc%ng_coarse, & 
                                             mlbddc%ptr_coarse_dofs, &
                                             mlbddc%f_mesh_c%lnods, &
                                             sum, &
                                             r_ci_gathered)
      else if ( mlbddc%co_sys_sol_strat == recursive_bddc ) then
         assert(mlbddc%p_mat%p_env%num_levels>2) ! Assuming last level direct
         call sum_coarse_stiffness_vectors ( r_c, &
                                             mlbddc%g_context%np, & 
                                             mlbddc%ng_coarse, & 
                                             mlbddc%ptr_coarse_dofs, &
                                             mlbddc%p_mesh_c%f_mesh%lnods, &
                                             sum, &
                                             r_ci_gathered, & 
                                             mlbddc%erenumbering_c%iperm)

      end if
      call memfree (  r_ci_gathered,__FILE__,__LINE__)

   else if (i_am_fine_task) then
      call snd_coarse_stiffness_vectors ( mlbddc%g_context%icontxt, &
           mlbddc%g_context%np,      &
           mlbddc%pad_collectives,          &
           mlbddc%nl_coarse,                &
           mlbddc%max_coarse_dofs,          &
           r_ci )

      if ( temp_fine_coarse_grid_overlap ) then
         call par_timer_stop ( mlbddc%timer_coll_applyprec  )
         call psb_barrier ( mlbddc%g_context%icontxt )
         call par_timer_start ( mlbddc%timer_applyprec_ov_fine  )
      end if

   end if


   if ( temp_fine_coarse_grid_overlap ) then 
      if ( .not. i_am_higher_level_task ) then
         call par_timer_stop ( mlbddc%timer_coll_applyprec  )
      end if
   end if

   if ( i_am_fine_task ) then
     call memfree ( r_ci,__FILE__,__LINE__)
   end if
 end subroutine par_preconditioner_dd_mlevel_bddc_compute_c_g_corr_ass_r_c


 subroutine par_preconditioner_dd_mlevel_bddc_compute_c_g_corr_scatter ( mlbddc, z_c, v_G )
   implicit none
   ! Parameters
   type(par_preconditioner_dd_mlevel_bddc_t), intent(in) :: mlbddc
   type(vector_t)                , intent(inout) :: z_c
   type(par_vector_t)                , intent(inout) :: v_G 

   ! Locals
   real(rp), allocatable :: z_ci(:)
   logical :: i_am_fine_task, i_am_coarse_task, i_am_higher_level_task

   ! Which duties I have?
   i_am_fine_task = (mlbddc%p_mat%p_env%p_context%iam >= 0)
   i_am_coarse_task = (mlbddc%c_context%iam >= 0)

   if ( i_am_fine_task ) then
      call memalloc ( mlbddc%nl_coarse, z_ci, __FILE__,__LINE__ )
   end if

   if ( temp_fine_coarse_grid_overlap ) then 
      if (i_am_fine_task) then
         call par_timer_stop ( mlbddc%timer_applyprec_ov_fine  ) 
      else if ( i_am_coarse_task ) then
         call par_timer_stop ( mlbddc%timer_applyprec_ov_coarse ) 
      end if
   end if

   ! 1.6. z_ci <- scatter z_c
   if ( i_am_fine_task ) then
      call rcv_z_c( mlbddc%g_context%icontxt, & 
                    mlbddc%g_context%np,      &
                    mlbddc%pad_collectives,          &
                    mlbddc%nl_coarse,                &
                    mlbddc%max_coarse_dofs,          &
                    z_ci )
   else if ( i_am_coarse_task ) then
      if (  mlbddc%co_sys_sol_strat == serial_gather ) then  
         call snd_z_c (mlbddc%g_context%icontxt, & 
                       mlbddc%g_context%np,      &
                       mlbddc%pad_collectives,          &
                       mlbddc%max_coarse_dofs,          &
                       mlbddc%ptr_coarse_dofs,          &
                       mlbddc%f_mesh_c%lnods,           &
                       z_c )
      else if ( mlbddc%co_sys_sol_strat == recursive_bddc ) then
         assert(mlbddc%p_mat%p_env%num_levels>2) ! Assuming last level direct
         call snd_z_c (mlbddc%g_context%icontxt, & 
                       mlbddc%g_context%np,      &
                       mlbddc%pad_collectives,          &
                       mlbddc%max_coarse_dofs,          &
                       mlbddc%ptr_coarse_dofs,          &
                       mlbddc%p_mesh_c%f_mesh%lnods,    &
                       z_c,                                &
                       mlbddc%erenumbering_c%iperm)
      end if
   end if

   if ( temp_fine_coarse_grid_overlap ) then 
      if (i_am_fine_task) then
         call par_timer_start ( mlbddc%timer_applyprec_ov_coarse_tail )
      end if
   end if

   if ( i_am_fine_task ) then
      ! 1.7. v1_G <- harm * z_ci | G
      call apply_harm ( mlbddc, z_ci, v_G )
      call memfree ( z_ci,__FILE__,__LINE__) 
   end if

 end subroutine par_preconditioner_dd_mlevel_bddc_compute_c_g_corr_scatter


  subroutine par_preconditioner_dd_mlevel_bddc_compute_fine_grid_correction ( mlbddc, r, v_G )
    implicit none
    ! Parameters
    type(par_preconditioner_dd_mlevel_bddc_t), intent(in) :: mlbddc
    type(par_vector_t)                , intent(in)    :: r       ! Local residual
    type(par_vector_t)                , intent(inout) :: v_G 

    ! Locals
    real(rp), allocatable :: r_r (:), z_r (:)
    real(rp), allocatable :: rhs (:), lambda_r (:)
    real(rp), allocatable :: work(:)
    real(rp), allocatable :: aug_r(:), aug_v(:)
    ! type(solver_control_t)      :: spars

    ! Phase 2: compute substructure correction v_G
    if ( mlbddc%nn_sys_sol_strat == corners_rest_part_solve_expl_schur ) then

       call memalloc ( mlbddc%A_rr_gr%nv, r_r, __FILE__,__LINE__ )

       call apply_perm_to_residual ( mlbddc%p_mat%dof_dist%nl , &
                                     mlbddc%nl_corners , & 
                                     mlbddc%perm, &                                  
                                     r%f_vector%b, &
                                     r_r )

       if ( mlbddc%kind_coarse_dofs == corners_and_edges .or. & 
            mlbddc%kind_coarse_dofs == corners_edges_and_faces .or. &
            mlbddc%kind_coarse_dofs == edges .or. & 
            mlbddc%kind_coarse_dofs == edges_and_faces .or. & 
            mlbddc%kind_coarse_dofs == faces .or. & 
            mlbddc%kind_coarse_dofs == corners_and_faces ) then
          call memalloc ( mlbddc%nl_edges , rhs, __FILE__,__LINE__ )
          call memalloc ( mlbddc%A_rr_gr%nv, work, __FILE__,__LINE__ )

          call compute_neumann_edge_lagrange_multipliers_rhs (mlbddc, r_r, rhs, work)

          ! write(*,*) rhs(1:10)

          call memalloc ( mlbddc%nl_edges  , lambda_r, __FILE__,__LINE__ )

              call solve_edge_lagrange_multipliers_explicit_schur   ( mlbddc, &
                                                                      mlbddc%S_rr, &
                                                                      mlbddc%S_rr_neumann, &
                                                                      mlbddc%ipiv, &
                                                                      mlbddc%ipiv_neumann, &
                                                                      mlbddc%schur_edge_lag_mult, &
                                                                      1, &
                                                                      rhs, lambda_r )

          call memfree ( rhs,__FILE__,__LINE__)
       else
          call memalloc ( 0,                         lambda_r, __FILE__,__LINE__ )

          call memalloc ( 0,                          work, __FILE__,__LINE__ )
       end if

       call memalloc ( mlbddc%A_rr_gr%nv, z_r, __FILE__,__LINE__ )

       call solve_neumann_problem (mlbddc, r_r, lambda_r, work, z_r)

       call memfree ( work,__FILE__,__LINE__)
       

       call apply_iperm_to_z_r ( mlbddc%unknowns, &  
                                 mlbddc%p_mat%dof_dist%nl , &
                                 mlbddc%p_mat%dof_dist%ni , &
                                 mlbddc%nl_corners , &
                                 mlbddc%iperm, &
                                 z_r, &
                                 v_G%f_vector%b )

       call memfree ( z_r,__FILE__,__LINE__)
       call memfree ( lambda_r,__FILE__,__LINE__) 
       call memfree ( r_r,__FILE__,__LINE__)

    else if ( mlbddc%nn_sys_sol_strat == direct_solve_constrained_problem ) then

       call memalloc ( (mlbddc%p_mat%dof_dist%nl+mlbddc%nl_coarse), aug_r, __FILE__, __LINE__ )
       call memalloc ( (mlbddc%p_mat%dof_dist%nl+mlbddc%nl_coarse), aug_v, __FILE__, __LINE__ )

       aug_r(1:mlbddc%p_mat%dof_dist%nl )  = r%f_vector%b
       aug_r(mlbddc%p_mat%dof_dist%nl+1:) = 0.0_rp

       mlbddc%spars_neumann%nrhs=1
       ! mlbddc%spars%method=direct
       aug_v = 0.0_rp
       if ( mlbddc%internal_problems == handled_by_bddc_module) then
          call solve( mlbddc%A_rr, mlbddc%M_rr, aug_r, aug_v, mlbddc%spars_neumann)
       else
          check(.false.)
!!$            call solve( mlbddc%A_rr_op, mlbddc%M_inv_op_rr, aug_r, aug_v, mlbddc%spars_neumann)
       end if


       v_G%f_vector%b = aug_v ( mlbddc%p_mat%dof_dist%ni+1:mlbddc%p_mat%dof_dist%nl   ) 

       call memfree ( aug_r,__FILE__,__LINE__)
       call memfree ( aug_v,__FILE__,__LINE__)
    end if

    !mlbddc%num_neumann_solves = mlbddc%num_neumann_solves + 1
    !mlbddc%neumann_its (mlbddc%num_neumann_solves) = mlbddc%spars%it

  end subroutine par_preconditioner_dd_mlevel_bddc_compute_fine_grid_correction


  !=============================================================================
  ! r_ci <- rPhi^T * r
  subroutine apply_harm_trans  (mlbddc, r, r_ci)
    implicit none

    ! Parameters
    type(par_preconditioner_dd_mlevel_bddc_t), intent(in)    :: mlbddc
    type(par_vector_t)                , intent(in)    :: r
    real(rp)                        , intent(inout) :: r_ci(mlbddc%nl_coarse  )

#ifdef ENABLE_BLAS
      r_ci = 0.0_rp
      ! r_ci <- 1.0 * rPhi^T * r + 0.0 * r_ci 
      ! if ( mlbddc%unknowns == all_unknowns ) then
      !    call DGEMV(  'T', & 
      !                mlbddc%p_mat%dof_dist%nl  , &
      !                mlbddc%nl_coarse  , &
      !                1.0, &
      !                mlbddc%rPhi, &
      !                mlbddc%p_mat%dof_dist%nl  , &
      !                r%f_vector%b, &
      !                1, &
      !                0.0, & 
      !                r_ci, & 
      !                1)
      !else if ( mlbddc%unknowns == interface_unknowns ) then
         call DGEMV(  'T', & 
                      mlbddc%p_mat%dof_dist%nb  , &
                      mlbddc%nl_coarse  , &
                      1.0, &
                      mlbddc%lPhi, &
                      mlbddc%p_mat%dof_dist%nb  , &
                      r%f_vector%b, &
                      1, &
                      0.0, & 
                      r_ci, & 
                      1)
      ! end if
#else
     write (0,*) 'Error: par_preconditioner_dd_mlevel_bddc was not compiled with -DENABLE_BLAS.'
     write (0,*) 'Error: You must activate this cpp macro in order to use the BLAS'
     check(.false.)    
#endif
   end subroutine apply_harm_trans

  !=============================================================================
  ! r <- rPhi * r_ci
  subroutine apply_harm  (mlbddc, r_ci, r)
    implicit none

    ! Parameters
    type(par_preconditioner_dd_mlevel_bddc_t), intent(in)    :: mlbddc
    real(rp)                        , intent(in)    :: r_ci(mlbddc%nl_coarse  )
    type(par_vector_t)                , intent(inout) :: r

#ifdef ENABLE_BLAS
      ! r <- 1.0 * rPhi * r_ci + 0.0 * r 
      r%f_vector%b = 0.0
      ! if ( mlbddc%unknowns == all_unknowns ) then
      !   call DGEMV(  'N', & 
      !                mlbddc%p_mat%dof_dist%nl  , &
      !                mlbddc%nl_coarse  , &
      !                1.0, &
      !                mlbddc%rPhi, &
      !                mlbddc%p_mat%dof_dist%nl  , &
      !                r_ci, &
      !                1,    &
      !                0.0,  & 
      !                r%f_vector%b, & 
      !                1)         
      !else if ( mlbddc%unknowns == interface_unknowns ) then
         call DGEMV(  'N', & 
                      mlbddc%p_mat%dof_dist%nb  , &
                      mlbddc%nl_coarse  , &
                      1.0, &
                      mlbddc%rPhi, &
                      mlbddc%p_mat%dof_dist%nb  , &
                      r_ci, &
                      1,    &
                      0.0,  & 
                      r%f_vector%b, & 
                      1)
      ! end if
#else
     write (0,*) 'Error: par_preconditioner_dd_mlevel_bddc was not compiled with -DENABLE_BLAS.'
     write (0,*) 'Error: You must activate this cpp macro in order to use the BLAS'
     check(.false.)    
#endif
   end subroutine apply_harm

  subroutine snd_coarse_stiffness_vectors (  icontxt_g, np_g, & 
                                             pad_collectives, & 
                                             nl_coarse, &
                                             max_coarse_dofs, &  
                                             r_ci )
                                                     
use psb_const_mod_names
use psb_penv_mod_names

#ifdef MPI_MOD
use mpi
#endif

    implicit none

#ifdef MPI_H
    include 'mpif.h'
#endif

    ! Parameters
    integer    , intent(in)  :: icontxt_g, np_g
    integer(ip), intent(in)  :: pad_collectives, nl_coarse
    integer(ip), intent(in)  :: max_coarse_dofs
    real   (rp), intent(in)  :: r_ci (nl_coarse)

    ! Locals
    integer                  :: mpi_comm_g, iret, root_g
    real   (rp), allocatable :: pad_buf_snd (:)
    real   (rp)              :: rcv_dum, rcv_dumvec(1)
    integer(ip)              :: int_rcv_dum, int_rcv_dumvec(1)

    ! Get MPI communicator associated to icontxt (in
    ! the current implementation of our wrappers
    ! to the MPI library icontxt and mpi_comm are actually 
    ! the same)
    call psb_get_mpicomm (icontxt_g, mpi_comm_g)
    root_g = np_g -1

    if ( pad_collectives == pad ) then
       call memalloc ( max_coarse_dofs, pad_buf_snd, __FILE__,__LINE__ )
 
       ! Pack r_ci into send buffer
       pad_buf_snd ( 1:nl_coarse ) = r_ci(:)

       ! write (*,*) 'XXXX', r_ci
       ! write (*,*) 'YYYY', pad_buf_snd
      
       call mpi_gather ( pad_buf_snd, max_coarse_dofs, psb_mpi_real, &
                         rcv_dum, int_rcv_dum, psb_mpi_real, &
                         root_g, mpi_comm_g, iret)
       check ( iret == mpi_success )

       call memfree ( pad_buf_snd,__FILE__,__LINE__)
    else if ( pad_collectives == nopad ) then
       call mpi_gatherv( r_ci,  nl_coarse, psb_mpi_real, &
                         rcv_dum, int_rcv_dumvec, int_rcv_dumvec, psb_mpi_real, &
                         root_g, mpi_comm_g, iret)
       check ( iret == mpi_success )
    end if 

  end subroutine snd_coarse_stiffness_vectors

  subroutine rcv_coarse_stiffness_vectors ( icontxt_g, np_g, & 
                                            pad_collectives, & 
                                            max_coarse_dofs, &  
                                            ptr_coarse_dofs, &
                                            sz_r_ci_gathered, &
                                            r_ci_gathered  )
use psb_const_mod_names
use psb_penv_mod_names

#ifdef MPI_MOD
use mpi
#endif

    implicit none

#ifdef MPI_H
    include 'mpif.h'
#endif

    ! Parameters
    integer    , intent(in)  :: icontxt_g, np_g
    integer(ip), intent(in)  :: pad_collectives
    integer(ip), intent(in)  :: ptr_coarse_dofs (np_g+1)
    integer(ip), intent(in)  :: max_coarse_dofs
    integer(ip), intent(in)  :: sz_r_ci_gathered 
    real   (rp), intent(out) :: r_ci_gathered (sz_r_ci_gathered)

    ! Locals
    integer                  :: mpi_comm_g, iret, root_g
    integer(ip)              :: i, k, off, offrcv, ncoarse
    integer(ip), allocatable :: recvcounts (:), displs (:)
    real   (rp), allocatable :: pad_buf_snd(:), pad_buf_rcv(:)
    real   (rp)              :: snd_dum

    ! Get MPI communicator associated to icontxt (in
    ! the current implementation of our wrappers
    ! to the MPI library icontxt and mpi_comm are actually 
    ! the same)
    call psb_get_mpicomm (icontxt_g, mpi_comm_g)
    root_g = np_g -1

    if ( pad_collectives == pad ) then
       call memalloc ( max_coarse_dofs*np_g, pad_buf_rcv, __FILE__,__LINE__ )
       call memalloc ( max_coarse_dofs, pad_buf_snd, __FILE__,__LINE__ )

       call mpi_gather ( pad_buf_snd, max_coarse_dofs, psb_mpi_real, &
                         pad_buf_rcv, max_coarse_dofs, psb_mpi_real, &
                         root_g, mpi_comm_g, iret)
       check ( iret == mpi_success )
          
       ! Unpack padded data
       offrcv = 0
       off    = 0 
       do i=1, np_g
         ncoarse = ptr_coarse_dofs(i+1)-ptr_coarse_dofs(i)
         r_ci_gathered ( off+1: off+ncoarse ) = pad_buf_rcv ( offrcv+1:offrcv+ncoarse ) 
         offrcv = offrcv + max_coarse_dofs
         off    = off    + ncoarse
       end do

       call memfree ( pad_buf_snd,__FILE__,__LINE__)
       call memfree ( pad_buf_rcv,__FILE__,__LINE__)
    else if ( pad_collectives == nopad ) then
       call memalloc ( np_g, recvcounts,   __FILE__,__LINE__ )
       ! call memalloc ( np_g+1, displs,  __FILE__,__LINE__ )

       ! displs(1) = 0 
       do i=1, np_g
         recvcounts(i) = (ptr_coarse_dofs(i+1)-ptr_coarse_dofs(i))
         ! displs(i+1)   = displs(i) + recvcounts(i)
       end do

       call mpi_gatherv( snd_dum, 0, psb_mpi_real, &
                         r_ci_gathered, recvcounts, ptr_coarse_dofs, psb_mpi_real, &
                         root_g, mpi_comm_g, iret)
       check ( iret == mpi_success )

       call memfree ( recvcounts,__FILE__,__LINE__)
       ! call memfree ( displs    ,__FILE__,__LINE__)
    end if

  end subroutine rcv_coarse_stiffness_vectors

  subroutine sum_coarse_stiffness_vectors  ( r_c,              &
                                             np, ng_coarse,    &
                                             ptr_coarse_dofs,  &
                                             lst_coarse_dofs,  & 
                                             sz_r_ci_gathered, &
                                             r_ci_gathered,    &
                                             perm )
    implicit none

    ! Parameters
    type(vector_t), intent(inout)         :: r_c
    integer(ip)     , intent(in)            :: np, ng_coarse, sz_r_ci_gathered
    integer(ip)     , intent(in)            :: ptr_coarse_dofs (np+1)
    integer(ip)     , intent(in)            :: lst_coarse_dofs (ptr_coarse_dofs(np+1))
    real   (rp)     , intent(in)            :: r_ci_gathered (sz_r_ci_gathered)
    integer(ip)     , optional, intent(in)  :: perm(np-1) ! Element re-numbering

    ! Locals
    integer (ip) :: i, nnode, inode, base, & 
                 &  ignode, base_ptr, base_dof
 
    ! write (*,*) ptr_coarse_dofs
    ! write (*,*) lst_coarse_dofs
    ! write (*,*) sz_r_ci_gathered
    ! write (*,*) ng_coarse   

     ! AFM: r_c should have been allocated in advance by the caller subroutine.
     !      As vector_alloc already initializes r_c to zero, the following
     !      sentence might not be needed. I leave it for safety reasons.
     call vector_zero(r_c)
     if ( .not. present(perm) ) then
        base_ptr = 0
        base_dof = 1
        do i=1, np-1
          nnode = ptr_coarse_dofs(i+1) - ptr_coarse_dofs(i)
          do inode=1, nnode
            ignode  = lst_coarse_dofs ( base_ptr + inode )
            ! base + (j-1)*lda + i
            r_c%b(ignode) = r_c%b(ignode) + r_ci_gathered(base_dof) 
            base_dof=base_dof+1 
          end do
          base_ptr = base_ptr + nnode
        end do
     else
        base_ptr = 1
        do i=1, np-1
           do inode=ptr_coarse_dofs(perm(i))+1, ptr_coarse_dofs(perm(i)+1)
              ignode  = lst_coarse_dofs ( base_ptr )
              r_c%b(ignode) = r_c%b(ignode) + r_ci_gathered(inode)
              base_ptr = base_ptr + 1
           end do
        end do
     end if
     
    ! write (*,*) r_c
  end subroutine sum_coarse_stiffness_vectors

  subroutine snd_z_c ( icontxt_g, np_g, &  
                       pad_collectives,  &
                       max_coarse_dofs,  &  
                       ptr_coarse_dofs,  &
                       lst_coarse_dofs,  & 
                       z_c, &
                       perm)
use psb_const_mod_names
use psb_penv_mod_names

#ifdef MPI_MOD
use mpi
#endif
    implicit none

#ifdef MPI_H
    include 'mpif.h'
#endif

    ! Parameters
    integer    , intent(in)      :: icontxt_g, np_g
    integer(ip), intent(in)      :: pad_collectives
    integer(ip), intent(in)      :: max_coarse_dofs
    integer(ip), intent(in)      :: ptr_coarse_dofs (np_g+1)
    integer(ip), intent(in)      :: lst_coarse_dofs (ptr_coarse_dofs (np_g+1))
    type(vector_t), intent(in) :: z_c
    integer(ip), optional, intent(in) :: perm(np_g-1)

    ! Locals
    integer (ip)              :: i, k, m, j, nnode, base
    integer (ip), allocatable :: sendcounts (:), displs (:)
    real (rp)   , allocatable :: z_c_scattered (:)
    real (rp)   , allocatable :: pad_buf_rcv(:)
    real (rp)                 :: rcv_dum

    ! Locals    
    integer     :: mpi_comm_g,  iret, root_g

    ! Get MPI communicator associated to icontxt (in
    ! the current implementation of our wrappers
    ! to the MPI library icontxt and mpi_comm are actually 
    ! the same)
    call psb_get_mpicomm (icontxt_g, mpi_comm_g)
    root_g = np_g -1

    if ( pad_collectives == pad ) then

       call memalloc ( max_coarse_dofs*np_g, z_c_scattered, __FILE__,__LINE__ ) 
       call memalloc ( max_coarse_dofs, pad_buf_rcv, __FILE__,__LINE__ )

       if ( .not. present(perm) ) then
          m = 1
          do j=1, np_g-1
            do i=ptr_coarse_dofs(j)+1, ptr_coarse_dofs(j+1)
               k                 = lst_coarse_dofs(i)
               z_c_scattered (m) = z_c%b(k)
               m                 = m + 1
            end do
            m = m + (max_coarse_dofs - (ptr_coarse_dofs(j+1)-ptr_coarse_dofs(j)))
          end do
       else
          m = 1
          do j=1, np_g-1
             nnode = ptr_coarse_dofs(perm(j)+1)-ptr_coarse_dofs(perm(j))
             base = (perm(j)-1)*max_coarse_dofs + 1
             do i=1, nnode
                k                    = lst_coarse_dofs(m)
                z_c_scattered (base) = z_c%b(k)
                m                    = m + 1
                base                 = base + 1
             end do
          end do
       end if

       call mpi_scatter  ( z_c_scattered, max_coarse_dofs,  psb_mpi_real, &
                           pad_buf_rcv  , max_coarse_dofs,  psb_mpi_real, &
                           root_g, mpi_comm_g, iret)
       check ( iret == mpi_success )

       call memfree ( pad_buf_rcv  , __FILE__, __LINE__)
       call memfree ( z_c_scattered, __FILE__, __LINE__)

    else if ( pad_collectives == nopad ) then 

       call memalloc ( np_g, sendcounts,  __FILE__,__LINE__ )
       ! call memalloc ( np_g+1, displs,   __FILE__,__LINE__ )

       ! displs(1) = 0 
       do i=1, np_g
          sendcounts(i) = (ptr_coarse_dofs(i+1)-ptr_coarse_dofs(i))
          ! displs(i+1)   = displs(i) + sendcounts(i)
       end do

       call memalloc ( ptr_coarse_dofs(np_g+1), z_c_scattered, __FILE__,__LINE__ ) 
 
       if ( .not. present(perm) ) then 
         m = 1
         do j=1, np_g-1
            do i=ptr_coarse_dofs(j)+1, ptr_coarse_dofs(j+1)
               k    = lst_coarse_dofs(i)
               z_c_scattered (m) = z_c%b(k)
               m                 = m   + 1
            end do
         end do
       else
         m = 1
         do j=1, np_g-1
            do i=ptr_coarse_dofs(perm(j))+1, ptr_coarse_dofs(perm(j)+1)
               k                 = lst_coarse_dofs(m)
               z_c_scattered (i) = z_c%b(k)
               m                 = m   + 1
            end do
         end do
      end if

       call mpi_scatterv ( z_c_scattered , sendcounts, ptr_coarse_dofs, psb_mpi_real, &
                           rcv_dum, 0,     psb_mpi_real, &
                           root_g, mpi_comm_g, iret)
       check ( iret == mpi_success )

       call memfree ( sendcounts,__FILE__,__LINE__)
       call memfree ( z_c_scattered,__FILE__,__LINE__)

    end if
          
  end subroutine snd_z_c

  subroutine rcv_z_c ( icontxt_g, np_g,  &  
                       pad_collectives,  &
                       nl_coarse,        &
                       max_coarse_dofs,  &  
                       z_ci  )
use psb_const_mod_names
use psb_penv_mod_names

#ifdef MPI_MOD
use mpi
#endif
    implicit none

#ifdef MPI_H
    include 'mpif.h'
#endif

    ! Parameters
    integer    , intent(in)  :: icontxt_g, np_g
    integer(ip), intent(in)  :: pad_collectives
    integer(ip), intent(in)  :: nl_coarse
    integer(ip), intent(in)  :: max_coarse_dofs
    real   (rp), intent(out) :: z_ci (nl_coarse)

    ! Locals
    integer (ip)              :: i, k, m, j
    real (rp)   , allocatable :: pad_buf_rcv(:)
    real (rp)                 :: snd_dum, snd_dumvec(1)
    integer (ip)              :: int_snd_dum, int_snd_dumvec(1)

    ! Locals    
    integer     :: mpi_comm_g, iret, root_g
    
    ! Get MPI communicator associated to icontxt (in
    ! the current implementation of our wrappers
    ! to the MPI library icontxt and mpi_comm are actually 
    ! the same)
    call psb_get_mpicomm (icontxt_g, mpi_comm_g)
    root_g = np_g -1

    if ( pad_collectives == pad ) then

       call memalloc ( max_coarse_dofs, pad_buf_rcv, __FILE__,__LINE__ )

       call mpi_scatter  ( snd_dum, int_snd_dum,  psb_mpi_real, &
                           pad_buf_rcv  , max_coarse_dofs,  psb_mpi_real, &
                           root_g, mpi_comm_g, iret)
       check ( iret == mpi_success )

       z_ci = pad_buf_rcv ( 1 : nl_coarse )
       call memfree ( pad_buf_rcv,__FILE__,__LINE__)

    else if ( pad_collectives == nopad ) then 

       call mpi_scatterv ( snd_dumvec, int_snd_dumvec, int_snd_dumvec, psb_mpi_real, &
                           z_ci, nl_coarse,     psb_mpi_real, &
                           root_g, mpi_comm_g, iret)
       check ( iret == mpi_success )
    end if
  end subroutine rcv_z_c


  ! r_r <- eliminate_corners ( P^T r )
  ! nl  <- number of local dofs
  ! nc  <- number of corner dofs
  subroutine apply_perm_to_residual ( nl, nc, perm, r, r_r ) 
    implicit none
    ! Parameters 
    integer (ip), intent(in)    :: nl, nc
    integer (ip), intent(in)    :: perm (nl)
    real (rp)   , intent(in)    :: r    (nl)
    real (rp)   , intent(out)   :: r_r  (nl-nc)

    ! Locals  
    integer (ip) :: i, j

    do i=1, nl
       j = perm(i)
       if ( j > nc ) then
           r_r ( j-nc ) = r(i) 
      end if
    end do

  end subroutine apply_perm_to_residual
  
  ! z_G <-  restrict_to_interface ( P  [z_c = 0])
  !                                    [z_r    ]
  ! nl  <- number of local dofs
  ! ni  <- number of interior dofs
  ! nc  <- number of corner dofs
  subroutine apply_iperm_to_z_r ( unk, nl, ni, nc, iperm, z_r, z_G ) 
    implicit none
    ! Parameters 
    integer (ip), intent(in)  :: unk, nl, ni, nc
    integer (ip), intent(in)  :: iperm (nl)
    real (rp)   , intent(in)  :: z_r (nl-nc)
    real (rp)   , intent(out) :: z_G (*)

    ! Locals  
    integer (ip) :: i, j

    ! if ( unk == interface_unknowns ) then
       do i=1, nl
          j = iperm(i)
          if ( j > ni ) then
             if ( i > nc ) then 
                z_G ( j-ni ) = z_r (i-nc)
             else ! i <= nc
                z_G ( j-ni ) = 0.0_rp
             end if
          end if
       end do
    ! else if ( unk == all_unknowns ) then
    !   do i=1, nl
    !      j = iperm(i)
    !      if ( i > nc ) then 
    !         z_G ( j ) = z_r (i-nc)
    !      else ! i <= nc
    !         z_G ( j ) = 0.0_rp
    !      end if
    !   end do
    ! end if

  end subroutine apply_iperm_to_z_r

  ! Computes the RHS of the substructure correction edge lagrange multipliers Schur complement
  ! system : (C_r A_rr^-1 C_r^T) lamba_r = C_r A_rr ^-1 r_r
  subroutine compute_neumann_edge_lagrange_multipliers_rhs (mlbddc, r_r, rhs, work, iparm, msglvl)
    implicit none
    ! Parameters 
    type(par_preconditioner_dd_mlevel_bddc_t), intent(in), target :: mlbddc

    real(rp)                 , intent(inout)         :: r_r  (mlbddc%A_rr_gr%nv)
    real(rp)                 , intent(out)           :: rhs  (mlbddc%nl_edges )
    real(rp)                 , intent(out)           :: work (mlbddc%A_rr_gr%nv) 

    integer           , intent(in), target, optional :: iparm(64)
    integer           , intent(in), optional         :: msglvl

    ! Locals
    integer(ip)             :: base, i

   ! type(solver_control_t)      :: spars

    ! work <- A_rr^-1 * r_r
   mlbddc%spars_neumann%nrhs=1
   ! mlbddc%spars%method=direct
   work = 0.0_rp
   if ( mlbddc%internal_problems == handled_by_bddc_module) then
      call solve( mlbddc%A_rr, mlbddc%M_rr, r_r, work, mlbddc%spars_neumann)
   else
      check(.false.)
!!$      call solve( mlbddc%A_rr_op, mlbddc%M_inv_op_rr, r_r, work, mlbddc%spars_neumann)
   end if


     rhs = 0.0_rp

    ! Y = alpha * Y + beta * C_r * X
    ! C_r => nl_edges x n(A_rr)

    call update_with_c_r   ( mlbddc%p_mat%dof_dist%nl, &
                             mlbddc%perm, &
                             mlbddc%iperm, &
                             mlbddc%nl_corners, &
                             mlbddc%nl_edges, &
                             mlbddc%coarse_dofs, &
                             mlbddc%p_mat%dof_dist%max_nparts, &
                             mlbddc%p_mat%dof_dist%omap%nl, &
                             mlbddc%p_mat%dof_dist%lobjs, &
                             mlbddc%nl_edges   , & 
                             mlbddc%A_rr_gr%nv, &
                             1, & 
                             work, & 
                             0.0_rp, &
                             1.0_rp, & 
                             rhs, &
                             mlbddc%p_mat%dof_dist%nb, &
                             mlbddc%C_weights )

  end subroutine compute_neumann_edge_lagrange_multipliers_rhs

  ! Solves the Neumann problem given by 
  ! A_rr z_r = r_r - C_r^T lamba_r 
  subroutine solve_neumann_problem (mlbddc, r_r, lambda_r, work, z_r)
    implicit none
    type(par_preconditioner_dd_mlevel_bddc_t), intent(in), target        :: mlbddc
    real(rp)                 , intent(inout)                :: r_r (mlbddc%A_rr_gr%nv)
    real(rp)                 , intent(in)                   :: lambda_r (mlbddc%nl_edges  )
    real(rp)                 , intent(in)                   :: work (mlbddc%A_rr_gr%nv) 
    real(rp)                 , intent(out)                  :: z_r  (mlbddc%A_rr_gr%nv)

    ! Locals
    real(rp), allocatable :: work2(:)
    integer(ip)           :: j, k, base
    ! type(solver_control_t)      :: spars

    if ( mlbddc%nn_sys_sol_strat == corners_rest_part_solve_expl_schur ) then
       if ( mlbddc%kind_coarse_dofs == corners_and_edges .or. & 
            mlbddc%kind_coarse_dofs == corners_edges_and_faces .or. &
            mlbddc%kind_coarse_dofs == edges .or. & 
            mlbddc%kind_coarse_dofs == edges_and_faces .or. & 
            mlbddc%kind_coarse_dofs == faces .or. & 
            mlbddc%kind_coarse_dofs == corners_and_faces ) then
         ! z_r <- A_rr^-1 r_r
         z_r = work
         
          ! write (*,*) lambda_r(1:10)

         ! z_r <- z_r - A_rr^-1 * C_r^T * lambda_r
#ifdef ENABLE_BLAS

         if (mlbddc%schur_edge_lag_mult == compute_from_scratch ) then
            call DGEMV( 'N', & 
                        mlbddc%A_rr_gr%nv, &
                        mlbddc%nl_edges, &
                        -1.0, &
                        mlbddc%A_rr_inv_C_r_T_neumann, &
                        mlbddc%A_rr_gr%nv, &
                        lambda_r, &
                        1, &
                        1.0, & 
                        z_r, & 
                        1)
         else
            call DGEMV( 'N', & 
                        mlbddc%A_rr_gr%nv, &
                        mlbddc%nl_edges , &
                        -1.0, &
                        mlbddc%A_rr_inv_C_r_T, &
                        mlbddc%A_rr_gr%nv, &
                        lambda_r, &
                        1, &
                        1.0, & 
                        z_r, & 
                        1)
            ! write (*,*) z_r(1:10)
            ! write (*,*) a_ci ! DBG:
         end if
#else
          write (0,*) 'Error: par_preconditioner_dd_mlevel_bddc was not compiled with -DENABLE_BLAS.'
          write (0,*) 'Error: You must activate this cpp macro in order to use the BLAS'
          check(.false.)    
#endif
       else if ( mlbddc%kind_coarse_dofs == corners ) then 
          mlbddc%spars_neumann%nrhs=1
          ! mlbddc%spars%method=direct
          z_r = 0.0_rp
          if ( mlbddc%internal_problems == handled_by_bddc_module) then
             call solve( mlbddc%A_rr, mlbddc%M_rr, r_r, z_r,mlbddc%spars_neumann)
          else
             check(.false.)
!!$             call solve( mlbddc%A_rr_op, mlbddc%M_inv_op_rr, r_r, z_r,mlbddc%spars_neumann)
          end if
       end if
    end if
  end subroutine solve_neumann_problem 

    subroutine par_preconditioner_dd_mlevel_bddc_report ( mlbddc, prefix )
use psb_const_mod_names
use psb_penv_mod_names
      
#ifdef MPI_MOD
use mpi
#endif  
      
      implicit none

#ifdef MPI_H
      include 'mpif.h'
#endif

      ! Parameters
      type(par_preconditioner_dd_mlevel_bddc_t), intent(inout) :: mlbddc
      character *(*)           , intent(in)    :: prefix 

      ! Locals
      integer(ip)           :: lunio, i, cur, j
      real(rp), allocatable :: gathered_data (:)
      real(rp), allocatable :: snd_buffer    (:)
      real(rp)              :: avg_max_dirichlet_its, avg_max_neumann_its, avg_coarse_solves, max_harm_its_agg, aux
      integer               :: iret 
      
      ! Get Rid of coarse-grid tasks 
      if (mlbddc%p_mat%p_env%p_context%iam<0) return

      avg_max_dirichlet_its = 0.0_rp
      avg_max_neumann_its   = 0.0_rp 

      if (mlbddc%p_mat%p_env%p_context%iam == 0) then
         lunio = io_open('dirichlet_report_' // trim(prefix) // '.' // trim(ch(mlbddc%p_mat%p_env%p_context%np)))
         call memalloc ( (4+mlbddc%num_dirichlet_solves)*mlbddc%p_mat%p_env%p_context%np, gathered_data, __FILE__,__LINE__ )
      else
         call memalloc ( 0, gathered_data, __FILE__,__LINE__ )
      end if
      
      call memalloc ( 4+mlbddc%num_dirichlet_solves, snd_buffer, __FILE__,__LINE__ )

      snd_buffer(1) = mlbddc%A_II_inv%M_II%cg
      snd_buffer(2) = mlbddc%A_II_inv%M_II%ca
      snd_buffer(3) = mlbddc%A_II_inv%M_II%cs
      snd_buffer(4) = mlbddc%A_II_inv%M_II%lev

      do j=1, mlbddc%num_dirichlet_solves
         snd_buffer(4+j) = mlbddc%dirichlet_its(j)
      end do

      
      call mpi_gather( snd_buffer   ,  4+mlbddc%num_dirichlet_solves, psb_mpi_real,  &
                       gathered_data,  4+mlbddc%num_dirichlet_solves, psb_mpi_real, 0, mlbddc%p_mat%p_env%p_context%icontxt, iret)

      if (mlbddc%p_mat%p_env%p_context%iam == 0) then
         cur=1
         aux = 0.0_rp
         do i=1, mlbddc%p_mat%p_env%p_context%np
            write (lunio, '(i10, 4f10.3)', advance='no') i, gathered_data(cur:cur+4-1)
            cur=cur+4
            do j=1, mlbddc%num_dirichlet_solves-1
               write (lunio, '(f10.3)',advance='no') gathered_data(cur)
               aux = aux + gathered_data(cur)
               cur=cur+1
            end do
            write (lunio, '(f10.3)') gathered_data(cur)
            aux = aux + gathered_data(cur)
            aux = aux/mlbddc%num_dirichlet_solves
            
            if ( aux > avg_max_dirichlet_its ) then
               avg_max_dirichlet_its = aux 
            end if

            cur=cur+1
         end do
      
      end if

      call memfree ( snd_buffer,__FILE__,__LINE__)
      call memfree ( gathered_data,__FILE__,__LINE__)  

      if (mlbddc%p_mat%p_env%p_context%iam == 0) then
         call io_close (lunio)
      end if

      if (mlbddc%p_mat%p_env%p_context%iam == 0) then
         lunio = io_open('neumann_report_' // trim(prefix) // '.' // trim(ch(mlbddc%p_mat%p_env%p_context%np)))
         call memalloc ( (6+mlbddc%num_neumann_solves)*mlbddc%p_mat%p_env%p_context%np, gathered_data, __FILE__,__LINE__ )
      else
         call memalloc ( 0, gathered_data, __FILE__,__LINE__ )
      end if

      call memalloc ( 6+mlbddc%num_neumann_solves, snd_buffer, __FILE__,__LINE__ )

      snd_buffer(1) = mlbddc%M_rr%cg
      snd_buffer(2) = mlbddc%M_rr%ca
      snd_buffer(3) = mlbddc%M_rr%cs
      snd_buffer(4) = mlbddc%M_rr%lev
      snd_buffer(5) = mlbddc%nl_coarse
      snd_buffer(6) = mlbddc%harm_its_agg

      do j=1, mlbddc%num_neumann_solves
         snd_buffer(6+j) = mlbddc%neumann_its(j)
      end do


      call mpi_gather( snd_buffer   ,  6+mlbddc%num_neumann_solves, psb_mpi_real,  &
                       gathered_data,  6+mlbddc%num_neumann_solves, psb_mpi_real, 0, mlbddc%p_mat%p_env%p_context%icontxt, iret)

      if (mlbddc%p_mat%p_env%p_context%iam == 0) then
         cur=1
         aux = 0.0_rp
         max_harm_its_agg = 0.0_rp
         do i=1, mlbddc%p_mat%p_env%p_context%np
            write (lunio, '(i10, 6f10.3)', advance='no') i, gathered_data(cur:cur+6-1)
            if ( gathered_data(cur+6-1) > max_harm_its_agg ) then
               max_harm_its_agg = gathered_data(cur+6-1)
            end if
            cur=cur+6
            do j=1, mlbddc%num_neumann_solves-1
               write (lunio, '(f10.3)',advance='no') gathered_data(cur)
               aux = aux + gathered_data(cur)
               cur=cur+1
            end do
            write (lunio, '(f10.3)') gathered_data(cur)
            
            aux = aux + gathered_data(cur)
            aux = aux/mlbddc%num_neumann_solves
            
            if ( aux > avg_max_neumann_its ) then
               avg_max_neumann_its = aux 
            end if
            
            cur=cur+1
         end do

         avg_coarse_solves = 0.0_rp
         do i=1, mlbddc%num_coarse_solves
            avg_coarse_solves = avg_coarse_solves + mlbddc%coarse_its(i)
         end do
         avg_coarse_solves = avg_coarse_solves/mlbddc%num_coarse_solves
         
      end if

      call memfree ( snd_buffer,__FILE__,__LINE__)
      call memfree ( gathered_data,__FILE__,__LINE__)

      if (mlbddc%p_mat%p_env%p_context%iam == 0) then
         call io_close (lunio)
         lunio = io_open('max_report_' // trim(prefix) // '.' // trim(ch(mlbddc%p_mat%p_env%p_context%np)))
         write (lunio, '(f10.3)') max_harm_its_agg, avg_max_dirichlet_its, avg_max_neumann_its, avg_coarse_solves
         call io_close (lunio)
      end if

    end subroutine par_preconditioner_dd_mlevel_bddc_report

    recursive subroutine par_preconditioner_dd_mlevel_bddc_time_report ( mlbddc )
      implicit none
      ! Parameters
      type(par_preconditioner_dd_mlevel_bddc_t), target, intent(inout) :: mlbddc

      ! Locals
      logical :: i_am_coarse_task, i_am_fine_task, i_am_higher_level_task

      assert ( associated(mlbddc%p_mat%p_env) )
      assert ( mlbddc%p_mat%p_env%created )
      assert ( mlbddc%p_mat%p_env%num_levels > 1 ) 

      assert ( associated(mlbddc%p_mat%p_env%w_context) )
      assert ( mlbddc%p_mat%p_env%w_context%created .eqv. .true. )
      assert ( associated(mlbddc%p_mat%p_env%p_context) )
      assert ( mlbddc%p_mat%p_env%p_context%created .eqv. .true. )
      assert ( associated(mlbddc%p_mat%p_env%q_context) )
      assert ( mlbddc%p_mat%p_env%q_context%created .eqv. .true. )
      ! Check appropriate assignment to context: 1) I'm in w_context and 2) I'm in p_context OR in q_context but not in both
      assert ( mlbddc%p_mat%p_env%w_context%iam >= 0)
      assert ( (mlbddc%p_mat%p_env%p_context%iam >=0 .and. mlbddc%p_mat%p_env%q_context%iam < 0) .or. (mlbddc%p_mat%p_env%p_context%iam < 0 .and. mlbddc%p_mat%p_env%q_context%iam >= 0))
      assert ( mlbddc%c_context%created .eqv. .true. )
      assert ( mlbddc%d_context%created .eqv. .true. )
      assert ( mlbddc%g_context%created .eqv. .true. )
      assert ( associated(mlbddc%p_mat%p_env%b_context) )
      assert ( mlbddc%p_mat%p_env%b_context%created .eqv. .true. )

      ! Which duties I have?
      i_am_fine_task = (mlbddc%p_mat%p_env%p_context%iam >= 0)
      i_am_coarse_task = (mlbddc%c_context%iam >= 0)
      i_am_higher_level_task = (mlbddc%d_context%iam >= 0)

      if ( temp_fine_coarse_grid_overlap ) then
        if (i_am_fine_task) then
          call par_timer_report (mlbddc%timer_assprec_ov_fine,   .false. )

          if (mlbddc%co_sys_sol_strat == serial_gather) then
             call par_timer_report (mlbddc%timer_assprec_ov_fine_header,   .false. )
          end if

          call par_timer_report (mlbddc%timer_fillprec_ov_fine,   .false. )

          if (mlbddc%co_sys_sol_strat == serial_gather) then
             call par_timer_report (mlbddc%timer_fillprec_ov_fine_header,   .false. )
             call par_timer_report (mlbddc%timer_fillprec_ov_coarse_header,   .false. )
          end if

          call par_timer_report (mlbddc%timer_applyprec_ov_fine,   .false. )

          if (mlbddc%co_sys_sol_strat == serial_gather) then
             call par_timer_report (mlbddc%timer_applyprec_ov_coarse_header,   .false. )
             call par_timer_report (mlbddc%timer_applyprec_ov_coarse_tail,   .false. )
          end if
          

        else if (i_am_coarse_task) then
          call par_timer_report (mlbddc%timer_assprec_ov_coarse, .false. )

          if (mlbddc%co_sys_sol_strat == recursive_bddc) then
             call par_timer_report (mlbddc%timer_assprec_ov_coarse_header, .false. )
          end if

          call par_timer_report (mlbddc%timer_fillprec_ov_coarse, .false. )
          call par_timer_report (mlbddc%timer_applyprec_ov_coarse, .false. )

        end if

        if ( i_am_fine_task .or. i_am_coarse_task ) then
          call par_timer_report ( mlbddc%timer_coll_assprec, .false.  )
          call par_timer_report ( mlbddc%timer_coll_fillprec, .false.  )
          call par_timer_report ( mlbddc%timer_coll_applyprec, .false.  )
        end if

        if ( i_am_coarse_task .or. i_am_higher_level_task ) then
          if (mlbddc%co_sys_sol_strat == recursive_bddc) then
            call par_preconditioner_dd_mlevel_bddc_time_report(mlbddc%p_M_c) 
          end if
        end if
      end if

    end subroutine par_preconditioner_dd_mlevel_bddc_time_report

! UPDATE CONSTRAINTS WEIGHTS VECTOR  ==========================================================

 subroutine assemble_C_weights ( mlbddc )
    implicit none

    ! Parameters  
    type(par_preconditioner_dd_mlevel_bddc_t), intent(inout)  :: mlbddc

    ! Locals
    real(rp), allocatable :: C_weights_i(:), C_weights_i_gathered(:)
    type(vector_t)      :: C_weights_next_level
    type(par_vector_t)      :: constraint_weights
    integer(ip)           :: i, j, sum
    logical               :: i_am_fine_task, i_am_coarse_task, i_am_higher_level_task
    integer(ip)           :: nl_coarse_dofs, ni, nb
    integer(ip)           :: nid, iobj, jd, j1, j2      

    ! Which duties I have  
    i_am_fine_task = (mlbddc%p_mat%p_env%p_context%iam >= 0)
    i_am_coarse_task = (mlbddc%c_context%iam >= 0)
    i_am_higher_level_task = (mlbddc%d_context%iam >= 0)

  if( i_am_coarse_task ) then
       sum = 0
       do i=1, mlbddc%g_context%np
          sum = sum + ( mlbddc%ptr_coarse_dofs(i+1)-mlbddc%ptr_coarse_dofs(i) )  
       end do
       call memalloc ( sum, C_weights_i_gathered, __FILE__,__LINE__ )
    end if
 
    if ( i_am_fine_task ) then
       ! BEGIN FINE-GRID RELATED DUTIES
       ! if (allocated(mlbddc%C_weights) ) then 
       !    !First of all compute local CORNER DOF over next level edge value 
       !    nid = mlbddc%p_mat%dof_dist%nl -  mlbddc%p_mat%dof_dist%nb
       !    do i=1, mlbddc%nl_coarse 
       !       iobj  = mlbddc%coarse_dofs (i)
       !       ! Each COARSE DOF is given its size
       !       j1 = mlbddc%p_mat%dof_dist%lobjs(2,iobj)
       !       j2 = mlbddc%p_mat%dof_dist%lobjs(3,iobj)

       !       do jd = j1, j2         
       !          if (j2-j1 .eq. 0) then 
       !             ! Dividing into 2 we are substracting the interior edge contribution 
       !             mlbddc%C_weights(jd-nid) = mlbddc%C_weights(jd-nid)/2
       !          end if
       !       end do
       !    end do
       ! end if
    
       call memalloc ( mlbddc%nl_coarse, C_weights_i, __FILE__,__LINE__ )
       C_weights_i = 0.0_rp

       if ( mlbddc%nl_coarse > 0 ) then           
              call vector_alloc(mlbddc%p_mat%dof_dist%nb, constraint_weights%f_vector)
              constraint_weights%f_vector%b = mlbddc%C_weights
              call apply_harm_trans( mlbddc, constraint_weights, C_weights_i ) 
              call vector_free ( constraint_weights%f_vector ) 
       end if
    end if

    if ( i_am_coarse_task ) then

              ! Assemble COARSE constraint weights
              C_weights_i_gathered = 0.0_rp
              call rcv_coarse_stiffness_vectors  ( mlbddc%g_context%icontxt, &
                                            mlbddc%g_context%np,      &
                                            mlbddc%pad_collectives,          &
                                            mlbddc%max_coarse_dofs,          &
                                            mlbddc%ptr_coarse_dofs,          &
                                            sum,                                &
                                            C_weights_i_gathered ) 

              if (temp_fine_coarse_grid_overlap ) then 
                 call psb_barrier ( mlbddc%g_context%icontxt )
              end if
      
         call vector_alloc ( mlbddc%ng_coarse, C_weights_next_level )
         call sum_coarse_stiffness_vectors ( C_weights_next_level, &
                                             mlbddc%g_context%np, & 
                                             mlbddc%ng_coarse, & 
                                             mlbddc%ptr_coarse_dofs, &
                                             mlbddc%p_mesh_c%f_mesh%lnods, &
                                             sum, &
                                             C_weights_i_gathered, & 
                                             mlbddc%erenumbering_c%iperm)

        ! Take only the interface weights from the full constraints weights vector 
        nl_coarse_dofs = mlbddc%p_M_c%nl_corners_dofs + mlbddc%p_M_c%nl_edges_dofs

        ni = mlbddc%p_mat_c%dof_dist%ni
        nb  = mlbddc%p_mat_c%dof_dist%nb

        call memalloc (  nb ,  mlbddc%p_M_c%C_weights , __FILE__, __LINE__) 
        mlbddc%p_M_c%C_weights = C_weights_next_level%b(ni+1:ni+nb)
        call vector_free ( C_weights_next_level ) 

        call memfree ( C_weights_i_gathered,__FILE__,__LINE__)

    elseif( i_am_fine_task ) then 
 
       call snd_coarse_stiffness_vectors  ( mlbddc%g_context%icontxt, &
                                            mlbddc%g_context%np,      &
                                            mlbddc%pad_collectives,          &
                                            mlbddc%nl_coarse,                &
                                            mlbddc%max_coarse_dofs,          &
                                            C_weights_i )
   
       if ( temp_fine_coarse_grid_overlap ) then
          call psb_barrier ( mlbddc%g_context%icontxt )
       end if

       call memfree (C_weights_i, __FILE__ , __LINE__ )
       
    end if

  end subroutine assemble_C_weights


  subroutine dof_distribution_coarse_create ( icontxt    , & ! Communication context
                                              me         , &
                                              np         , &
                                              f_mesh     , & ! Local mesh in an initial ARBITRARY local node/element ordering
                                              l2ge       , & ! local2global element correspondence
                                              df_mesh    , & ! Associated dual_mesh with external elements also listed for boundary DoFs
                                              dual_parts , & ! Parts associated to each element in the dual mesh
                                              vars       , & ! Physical unknown corresponding to each DoF in f_mesh
                                              new_f_mesh , & ! New mesh object
                                              dof_dist   , & ! New dof_distribution object conformal to new local mesh
                                              erenumbering       )   ! Resulting Element renumbering
                                                                
                                              
    implicit none
    integer(ip)           , intent(in)    :: icontxt, me, np 
    type(mesh_t)        , intent(in)    :: f_mesh
    integer(ip)           , intent(in)    :: l2ge(f_mesh%nelem)
    type(mesh_t)        , intent(in)    :: df_mesh
    integer(ip)           , intent(in)    :: dual_parts(df_mesh%pnods(df_mesh%nelem+1)-1)
    integer(ip)           , intent(in)    :: vars(f_mesh%npoin)
    type(mesh_t)        , intent(out)   :: new_f_mesh
    type(dof_distribution_t), intent(out)   :: dof_dist
    type(renumbering_t)           , intent(out)   :: erenumbering    


    ! Locals
    integer(ip)                :: nboun, eboun, max_nelem
    integer(igp)               :: ngelem
    integer(ip) , allocatable  :: nlboun(:), elboun(:)
    type(renumbering_t)                :: nrenumbering
    type(par_timer_t)            :: t1, t2, t3, t4, t5
    logical, parameter         :: temporize_phases = .false. 

    dof_dist%nparts = np
    dof_dist%ipart  = me+1
    
    if ( temporize_phases ) then
      call par_timer_create ( t1, 'count_nboun_list_nlboun', ictxt=icontxt )
      call par_timer_start  ( t1 )
    end if
    
    call count_nboun_eboun_list_nlboun_elboun_mlevel_bddc ( me+1      , &
                                                            f_mesh    , &
                                                            l2ge      , &
                                                            df_mesh   , &
                                                            dual_parts, &
                                                            nboun     , &
                                                            eboun     , &
                                                            nlboun    , &
                                                            elboun ) 
    if ( temporize_phases ) then
       call par_timer_stop ( t1 )
    end if
    

    nrenumbering%n = f_mesh%npoin
    erenumbering%n = f_mesh%nelem

    call memalloc (nrenumbering%n, nrenumbering%iperm, __FILE__,__LINE__)
    call memalloc (erenumbering%n, erenumbering%iperm, __FILE__,__LINE__)

    if ( temporize_phases ) then
      call par_timer_create ( t2, 'create_lpadj_lobjs', ictxt=icontxt )
      call par_timer_start  ( t2 )
    end if

    call create_lpadj_lobjs_with_vars_mlevel_bddc ( f_mesh, &
                                                    df_mesh, &
                                                    dual_parts, &
                                                    icontxt, &
                                                    me+1, &
                                                    np, &
                                                    nboun, &
                                                    eboun, &
                                                    elboun, &
                                                    vars, &
                                                    nlboun, &
                                                    dof_dist%npadj      , &
                                                    dof_dist%lpadj      , & 
                                                    dof_dist%max_nparts , &
                                                    max_nelem             , &
                                                    dof_dist%nobjs      , & 
                                                    dof_dist%lobjs      , &   
                                                    nrenumbering%iperm            , &
                                                    erenumbering%iperm )

    
    if ( temporize_phases ) then
      call par_timer_stop ( t2 )
    end if

    call memfree (nlboun,__FILE__,__LINE__)
    call memfree (elboun,__FILE__,__LINE__)

    call memalloc (nrenumbering%n, nrenumbering%lperm, __FILE__,__LINE__)
    call memalloc (erenumbering%n, erenumbering%lperm, __FILE__,__LINE__)

    call renumbering_inverse (nrenumbering%n, nrenumbering%iperm, nrenumbering%lperm)
    call renumbering_inverse (erenumbering%n, erenumbering%iperm, erenumbering%lperm)

    if ( temporize_phases ) then
      call par_timer_create ( t3, 'create_int_objs', ictxt=icontxt )
      call par_timer_start  ( t3 )
    end if

    call create_int_objs ( me+1, &
                           dof_dist%npadj, &
                           dof_dist%lpadj, &
                           dof_dist%max_nparts , &
                           dof_dist%nobjs      , &
                           dof_dist%lobjs      , &
                           dof_dist%int_objs%n , &
                           dof_dist%int_objs%p , &
                           dof_dist%int_objs%l ) 
    if ( temporize_phases ) then
      call par_timer_stop ( t3 )
    end if

    if ( temporize_phases ) then
      call par_timer_create ( t4, 'create_omap', ictxt=icontxt )
      call par_timer_start  ( t4 )
    end if
    call create_omap ( icontxt    , & ! Communication context
                       me         , &
                       np         , &
                       dof_dist%npadj, &
                       dof_dist%lpadj, & 
                       dof_dist%int_objs%p, &  
                       dof_dist%int_objs%l, &
                       dof_dist%max_nparts , & 
                       dof_dist%nobjs      , & 
                       dof_dist%lobjs      , &
                       dof_dist%omap%nl,     &
                       dof_dist%omap%ng,     &
                       dof_dist%omap%ni,     &
                       dof_dist%omap%nb,     &
                       dof_dist%omap%ne,     &
                       dof_dist%omap%l2g )
    if ( temporize_phases ) then
      call par_timer_stop ( t4 )
    end if

    if ( temporize_phases ) then
      ! Apply nrenumbering/erenumbering to input mesh into output mesh
      call par_timer_create ( t5, 'mesh_l2l', ictxt=icontxt )
      call par_timer_start  ( t5 )
    end if

    call mesh_l2l ( nrenumbering, erenumbering, f_mesh, new_f_mesh )
 
    if ( temporize_phases ) then
      call par_timer_stop ( t5 )
    end if

    if ( temporize_phases ) then
      call par_timer_report (t1)
      call par_timer_report (t2, .false.)
      call par_timer_report (t3, .false.)
      call par_timer_report (t4, .false.)
      call par_timer_report (t5, .false.)
    end if

    call renumbering_free(nrenumbering)

    ! Compute dof_import_t instance such that DoF nearest neighbour exchanges
    ! can be performed among subdomains on any intermmediate coarse-grid mesh
    call dof_distribution_compute_import(dof_dist)

  end subroutine dof_distribution_coarse_create


  subroutine count_nboun_eboun_list_nlboun_elboun_mlevel_bddc ( ipart     , &
                                                                f_mesh    , &
                                                                l2ge      , &
                                                                df_mesh   , &
                                                                dual_parts, &
                                                                nboun     , &
                                                                eboun     , &
                                                                nlboun    , &
                                                                elboun )
    implicit none

    ! Parameters
    integer(ip)   , intent(in)                :: ipart
    type(mesh_t), intent(in)                :: f_mesh
    integer(ip)   , intent(in)                :: l2ge(f_mesh%nelem)
    type(mesh_t), intent(in)                :: df_mesh
    integer(ip)   , intent(in)                :: dual_parts(df_mesh%pnods(df_mesh%nelem+1)-1)
    integer(ip)   , intent(out)               :: nboun
    integer(ip)   , intent(out)               :: eboun
    integer(ip)   , intent(out), allocatable  :: nlboun(:)
    integer(ip)   , intent(out), allocatable  :: elboun(:)


    ! Locals
    integer(ip), allocatable :: auxn(:), auxe(:)
    integer(ip)              :: i, j, istat, aux_val
    logical                  :: is_boun_poin
    type(hash_table_ip_ip_t)   :: boun_elems 
    type(hash_table_ip_ip_t)   :: g2le

    call boun_elems%init(max(int(real(f_mesh%npoin,rp)*0.05_rp,ip),5))
    call g2le%init(max(int(real(f_mesh%nelem,rp)*0.25_rp,ip),5))

    call memalloc (f_mesh%npoin, auxn, __FILE__, __LINE__)
    call memalloc (f_mesh%nelem, auxe, __FILE__, __LINE__)

    do i=1, f_mesh%nelem
       aux_val = i
       ! write(*,*) 'XXX', i, l2ge(i)
       call g2le%put(key=l2ge(i), val=aux_val, stat=istat)
    end do

    nboun=0
    eboun=0
    do i=1, f_mesh%npoin ! For each DoF of the primal mesh
       is_boun_poin = .false.
       do j=df_mesh%pnods(i), df_mesh%pnods(i+1)-1 
          if (.not. is_boun_poin .and. dual_parts(j) /= ipart) then
             nboun      = nboun + 1
             auxn(nboun) = i 
             is_boun_poin = .true.
             exit
          end if
       end do

       if ( is_boun_poin ) then
          do j=df_mesh%pnods(i), df_mesh%pnods(i+1)-1 
             ! write(*,*) 'XXX', i, dual_parts(j), ipart, df_mesh%lnods(j)
             if (dual_parts(j) == ipart) then
                aux_val = 1
                call boun_elems%put(key=df_mesh%lnods(j), val=aux_val, stat=istat)
                if (istat == now_stored) then
                   eboun       = eboun + 1
                   ! auxe(eboun) = df_mesh%lnods(j); df_mesh%lnods(j) holds a global EID
                   ! write(*,*) 'YYY', df_mesh%lnods(j)
                   ! call g2le%print
                   call g2le%get(key=df_mesh%lnods(j),val=auxe(eboun), stat=istat) 
                   ! write(*,*) 'YYY', auxe(eboun)
                end if
             end if
          end do
       end if
    end do

    call memalloc (nboun, nlboun, __FILE__,__LINE__)
    nlboun = auxn(1:nboun)
    call memfree  (auxn, __FILE__, __LINE__)

    call memalloc (eboun, elboun, __FILE__,__LINE__)
    elboun = auxe(1:eboun)
    call memfree  (auxe, __FILE__, __LINE__)

    call boun_elems%free
    call g2le%free

  end subroutine count_nboun_eboun_list_nlboun_elboun_mlevel_bddc


  subroutine create_lpadj_lobjs_with_vars_mlevel_bddc ( f_mesh     , &
                                                        df_mesh    , &
                                                        dual_parts , &
                                                        icontxt    , &
                                                        ipart      , &
                                                        nparts     , &
                                                        nboun      , &
                                                        eboun      , &
                                                        elboun     , &
                                                        vars       , &
                                                        nlboun     , &
                                                        npadj      , &
                                                        lpadj      , & 
                                                        max_nparts , &
                                                        max_nelem  , & 
                                                        nobjs      , & 
                                                        lobjs      , &
                                                        nl2ln2o    , &
                                                        el2ln2o    )
    implicit none

    type(mesh_t), intent(in) :: f_mesh
    type(mesh_t), intent(in) :: df_mesh
    integer(ip)   , intent(in) :: dual_parts(df_mesh%pnods(df_mesh%nelem+1)-1)
    integer(ip)   , intent(in) :: icontxt
    integer(ip)   , intent(in) :: ipart
    integer(ip)   , intent(in) :: nparts 
    integer(ip)   , intent(in) :: nboun
    integer(ip)   , intent(in) :: eboun
    integer(ip)   , intent(in) :: elboun(eboun)

    integer(ip) , intent(in)    :: vars(f_mesh%npoin)
    integer(ip) , intent(inout) :: nlboun(nboun)

    integer(ip), intent(out)    :: npadj 
    integer(ip), intent(out), allocatable :: lpadj(:)
    integer(ip), intent(out)    :: max_nparts
    integer(ip), intent(out)    :: max_nelem
    integer(ip), intent(out)    :: nobjs
    integer(ip), intent(out), allocatable :: lobjs(:,:)
    integer(ip), intent(out)    :: nl2ln2o(f_mesh%npoin)
    integer(ip), intent(out)    :: el2ln2o(f_mesh%nelem)


    ! Locals
    integer(ip)               :: est_max_nparts, i, j, k, jpart, count, pos, count_nelem, pos_jelem
    integer(ip)               :: start, end, ni, ielpo, jelem
    integer(igp), allocatable :: ws_parts_list_sep (:,:)
    type(hash_table_ip_ip_t)    :: ws_parts_visited, ws_parts_visited_all
    type(hash_table_igp_ip_t)   :: ws_elems_visited
    ! integer(ip), allocatable :: ws_parts_visited (:)
    ! integer(ip), allocatable :: ws_parts_visited_all (:)
    integer(ip), allocatable  :: ws_parts_visited_list_all (:)
    integer(igp), allocatable :: ws_sort_l1 (:), ws_sort_l2 (:)
    integer(ip), allocatable  :: ws_lobjs_temp (:,:)
    integer(ip), allocatable  :: is_boun (:), eboun_pos(:)
    integer(ip)               :: cur_internal_elem, cur_boundary_elem
    integer(ip)               :: istat
    integer(ip), parameter    :: tbl_length = 100

    type(par_timer_t)    :: t1, t2, t3, t4, t5
    logical, parameter :: temporize_phases = .false. 
    integer(ip) :: aux_val 


    
    ! A. Martin: This is dirty and deprecated. We should declare an interface 
    ! to icomp somewhere !!! I guess that the same has to
    ! be applied for intsort and sortix sort subroutines (TODO)
    !integer icomp
    !external icomp   

    call memalloc ( f_mesh%nelem, eboun_pos, __FILE__, __LINE__ )

    if ( temporize_phases ) then
      call par_timer_create ( t1, 'lpadj_lobjs_vars:initial set-up', ictxt=icontxt )
      call par_timer_start  ( t1 )
    end if

    eboun_pos = -1
    do i=1, eboun
       eboun_pos(elboun(i))=i
    end do

    cur_internal_elem = 1 
    cur_boundary_elem = f_mesh%nelem - eboun + 1

    do i=1, f_mesh%nelem
       if (eboun_pos(i)==-1) then
          el2ln2o(cur_internal_elem) = i 
          cur_internal_elem = cur_internal_elem + 1
       else
          el2ln2o (cur_boundary_elem) = i 
          cur_boundary_elem = cur_boundary_elem + 1
       end if
    end do

    call memalloc ( f_mesh%npoin, is_boun, __FILE__, __LINE__ )

    ! Compute an estimation for the maximum number 
    ! of parts around any boundary node
    is_boun        = 0
    est_max_nparts = 0
    do i=1, nboun
       is_boun( nlboun(i) ) = 1
       count = df_mesh%pnods(nlboun(i)+1)-df_mesh%pnods(nlboun(i))
       if ( count > est_max_nparts) then
          est_max_nparts=count
       end if
    end do

    ! write(*,*) 'is_boun', is_boun
    ! write(*,*) 'est_max', est_max_nparts 

    ! Re-number internal nodes first
    k=1
    do i=1, f_mesh%npoin
       if ( is_boun(i) == 0 ) then
          nl2ln2o(k)=i
          k=k+1
       end if
    end do

    if ( temporize_phases ) then
      call par_timer_stop  ( t1 )
    end if

    call memfree ( is_boun,__FILE__,__LINE__) 


    call memalloc ( 2*est_max_nparts+2, nboun, ws_parts_list_sep,                     __FILE__,__LINE__ )

    ! call memalloc ( nparts, ws_parts_visited,                __FILE__,__LINE__ )

    ! call memalloc ( nparts, ws_parts_visited_all,                __FILE__,__LINE__ )

    ! The size of the following array does not scale 
    ! properly with nparts. A hash table (in its current form)
    ! does not solve the problem as there is no subroutine to
    ! extract all the elements in one-shot, or one by one iteratively.
    ! Anyway, if these subroutines would exist, I do not consider
    ! that an associative table is the best option to implement
    ! in an efficient way this data structure. A dynamic linked list
    ! may be sufficient. Volunteers?
    call memalloc ( nparts, ws_parts_visited_list_all,                     __FILE__,__LINE__ )


    if ( temporize_phases ) then
      call par_timer_create ( t2, 'lpadj_lobjs_vars:subds_surrounding', ictxt=icontxt )
      call par_timer_start  ( t2 )
    end if

    aux_val              = 1 
    max_nparts           =0
    max_nelem            =0
    npadj                =0
    ! ws_parts_visited_all =0
    call ws_parts_visited_all%init(tbl_length)
    ws_parts_list_sep    =0
    ! ws_parts_visited     =0
    do i=1,nboun

       call ws_parts_visited%init(tbl_length)
       call ws_elems_visited%init(tbl_length)

       ! write (*,*) 'XXX', ws_parts_visited
       ! ws_parts_visited(ipart) = 1
       call ws_parts_visited%put(key=ipart,val=aux_val,stat=istat)
       ws_parts_list_sep (1,i) = vars(nlboun(i))
       ws_parts_list_sep (3,i) = ipart
       count                   = 1 

       pos_jelem   = est_max_nparts+3
       count_nelem = 0

       ! Traverse elements around nlboun(i)
       do ielpo=df_mesh%pnods(nlboun(i)), df_mesh%pnods(nlboun(i)+1)-1
          ws_parts_list_sep (pos_jelem, i) = df_mesh%lnods(ielpo)
          pos_jelem   = pos_jelem   + 1 
          count_nelem = count_nelem + 1

          ! YES: extract part of neighbour
          jpart = dual_parts(ielpo)

          call ws_parts_visited%put(key=jpart,val=aux_val,stat=istat)
          if ( istat == now_stored ) then
             count = count + 1
             ws_parts_list_sep (count+2,i) = jpart
             ! ws_parts_visited(jpart) = 1 
          end if

          call ws_parts_visited_all%put(key=jpart,val=aux_val,stat=istat)
          if ( istat == now_stored ) then
             npadj = npadj + 1
             ws_parts_visited_list_all(npadj) = jpart
             ! ws_parts_visited_all(jpart) = 1
          end if

       end do

       ws_parts_list_sep(2,i) = count

       ! Re-set ws_parts_visited
       ! do j=2, count+1
       !   ws_parts_visited (ws_parts_list_sep(j,i)) = 0 
       ! end do
       call ws_parts_visited%free

       ! Sort list of parts in increasing order by part identifiers
       ! This is required by the call to icomp subroutine below 
       ! call psb_hsort ( ws_parts_list_sep( 3:(count+2), i) )
       call sort ( (count+2)-3+1, ws_parts_list_sep( 3:(count+2), i) )
       
       ! call psb_hsort ( ws_parts_list_sep( est_max_nparts+3:pos_jelem-1, i) )
       call sort ( (pos_jelem-1)-(est_max_nparts+3)+1, ws_parts_list_sep( est_max_nparts+3:pos_jelem-1, i) )


       ! if ( count == 4 ) then
       ! write (*,*) 'PARTS', ws_parts_list_sep( 2:(count+1), i)
       ! write (*,*) 'ELEMS', ws_parts_list_sep( est_max_nparts+2:pos_jelem-1, i)
       !end if


       ! write (*,*) 'ZZZ', i, ws_parts_list_sep( 3:(count+2), i)
       ! write (*,*) 'ZZZ', i, ws_parts_list_sep( est_max_nparts+3:pos_jelem-1, i)
       ! write (*,*) i, count_nelem

       max_nparts = max(max_nparts, count)
       max_nelem  = max(max_nelem , count_nelem)
    end do

    if ( temporize_phases ) then
      call par_timer_stop  ( t2 )
    end if

    call memfree( eboun_pos,__FILE__,__LINE__)

    ! call memfree ( ws_parts_visited,__FILE__,__LINE__)

    ! call memfree ( ws_parts_visited_all,__FILE__,__LINE__)

    call ws_parts_visited_all%free

    call memalloc ( npadj, lpadj,                     __FILE__,__LINE__ )

    lpadj = ws_parts_visited_list_all(1:npadj)

    ! write (*,*) 'XXX', npadj, lpadj

    call memfree ( ws_parts_visited_list_all,__FILE__,__LINE__)

    call memalloc ( 2*est_max_nparts+2, ws_sort_l1,                     __FILE__,__LINE__ )

    call memalloc ( 2*est_max_nparts+2, ws_sort_l2,                     __FILE__,__LINE__ )


    ! write(*,*) 'XXX',  ws_parts_list_sep(1:max_nparts+1, 10486-(f_mesh%npoin-nboun))
    ! write(*,*) 'XXX',  ws_parts_list_sep(1:max_nparts+1, 24243-(f_mesh%npoin-nboun))


    if ( temporize_phases ) then
      call par_timer_create ( t3, 'lpadj_lobjs_vars:intsort', ictxt=icontxt )
      call par_timer_start  ( t3 )
    end if

    ! Re-number boundary nodes in increasing order by the 
    ! number of parts they belong and, among vertices sharing the same parts,
    ! in increasing order by the list of parts shared by each vertex
    !call intsort( 2*est_max_nparts+2,       & ! Rows of ws_parts_list_sep
    !   &          2*est_max_nparts+2,   & ! LD of ws_parts_list_sep
    !   &          nboun,              & ! Cols of ws_parts_list_sep
    !   &          ws_parts_list_sep,  &
    !   &          nlboun, &
    !   &          ws_sort_l1, ws_sort_l2)
    call sort_array_cols_by_row_section( 2*est_max_nparts+2,       & ! Rows of ws_parts_list_sep
       &          2*est_max_nparts+2,   & ! LD of ws_parts_list_sep
       &          nboun,              & ! Cols of ws_parts_list_sep
       &          ws_parts_list_sep,  &
       &          nlboun, &
       &          ws_sort_l1, ws_sort_l2)

!!$    do i=1,nboun
!!$       write(*,*) 'WSXXX', ws_parts_list_sep(:,i)
!!$    end do


    if ( temporize_phases ) then
      call par_timer_stop  ( t3 )
    end if

    ! write(*,*) nlboun
    ! write(*,*) 'XXX', ws_parts_list_sep(1:max_nparts+1, 1)


    call memfree ( ws_sort_l1,__FILE__,__LINE__)

    call memfree ( ws_sort_l2,__FILE__,__LINE__)

    ! Sincronize tasks
    call psb_max ( icontxt, max_nparts ) 

    ! Identify interface objects 
    ! (faces, edges and vertices) 
    call memalloc ( max_nparts+4, nboun, ws_lobjs_temp,                     __FILE__,__LINE__ )

    if ( temporize_phases ) then
      call par_timer_create ( t4, 'lpadj_lobjs_vars:pack_objs', ictxt=icontxt )
      call par_timer_start  ( t4 )
    end if

    nobjs=1

    ! Prepare first object
    ws_lobjs_temp(1,nobjs) = ws_parts_list_sep(1,1)  ! Extract physical unknown 
    ws_lobjs_temp(2,nobjs) = f_mesh%npoin-nboun+1    ! begin obj
    ws_lobjs_temp(4,nobjs) = ws_parts_list_sep(2,1)
    ws_lobjs_temp(5:4+ws_parts_list_sep(2,1),nobjs) = &
       &          ws_parts_list_sep(3:2+ws_parts_list_sep(2,1),1)

    ni = f_mesh%npoin-nboun
    k  = ni + 1

    ! Loop over vertices of the current separator (i.e., inode)
    do i=1,nboun-1

       if(icomp(max_nparts+2,ws_parts_list_sep(:,i),ws_parts_list_sep(:,i+1)) /= 0) then

          ! Complete current object
          ws_lobjs_temp(3,nobjs)=k ! end obj
          ws_lobjs_temp(4,nobjs)=ws_parts_list_sep(2,i)
          ws_lobjs_temp(5:4+ws_parts_list_sep(2,i),nobjs) = &
               & ws_parts_list_sep(3:2+ws_parts_list_sep(2,i),i)

          ! Prepare next object
          nobjs=nobjs+1
          ws_lobjs_temp(1,nobjs)=ws_parts_list_sep(1,i) ! Extract physical unknown
          ws_lobjs_temp(2,nobjs)=k+1                    ! begin obj
       end if
       k=k+1
    end do
  
    ! Complete last object
    ws_lobjs_temp(3,nobjs)=k        ! end obj
    ws_lobjs_temp(4,nobjs)=ws_parts_list_sep(2,i)
    ws_lobjs_temp(5:4+ws_parts_list_sep(2,i),nobjs) = &
         & ws_parts_list_sep(3:2+ws_parts_list_sep(2,i),i)    

    nl2ln2o(ni+1:) = nlboun

    nobjs = nobjs +1 

    if ( temporize_phases ) then
      call par_timer_stop  ( t4 )
    end if

    ! Reallocate lobjs and add internal object first
    call memalloc (max_nparts+4, nobjs, lobjs,__FILE__,__LINE__)
  

    if ( temporize_phases ) then
      call par_timer_create ( t5, 'lpadj_lobjs_vars:generate_lobjs', ictxt=icontxt )
      call par_timer_start  ( t5 )
    end if
    
    ! Internal object
    lobjs (1,1) = -1
    lobjs (2,1) = 1
    lobjs (3,1) = ni
    lobjs (4,1) = 1
    lobjs (5,1) = ipart
    lobjs (6:max_nparts+4,1) = 0

    ! Copy ws_lobjs_temp to lobjs ...
    do i=1,nobjs-1
       lobjs(:,i+1)=ws_lobjs_temp(1:(max_nparts+4), i)
    end do

    if ( temporize_phases ) then
      call par_timer_stop  ( t5 )
    end if

    if ( temporize_phases ) then
      call par_timer_report (t1, .false.)
      call par_timer_report (t2, .false.)
      call par_timer_report (t3, .false.)
      call par_timer_report (t4, .false.)
      call par_timer_report (t5, .false.)
    end if

    call memfree ( ws_parts_list_sep,__FILE__,__LINE__)

    call memfree ( ws_lobjs_temp,__FILE__,__LINE__)
  end subroutine create_lpadj_lobjs_with_vars_mlevel_bddc


  !=============================================================================
  subroutine matrix_assembly(nn,ln,ea,mat)
    implicit none
    ! Parameters
    integer(ip)      ,intent(in)    :: nn
    integer(ip)      ,intent(in)    :: ln(:)
    real(rp)         ,intent(in)    :: ea(:,:)
    type(matrix_t) ,intent(inout) :: mat
    
    ! Locals
    integer(ip) :: nl, ne
    
    nl=size(ln)
    assert(nl==nn)
    
    ne=size(ea,dim=1)
    assert(ne==nn)
    
    ne=size(ea,dim=2)
    assert(ne==nn)
    
    assert(mat%type==csr_mat)
    
    call ass_csr_mat_scal(mat%symm,nn,ln,ea,mat%gr%nv, &
         &                mat%gr%ia,mat%gr%ja,mat%a)

  end subroutine matrix_assembly

  subroutine ass_csr_mat_scal (ks,nn,ln,ea,nv,ia,ja,a)
    implicit none
    ! Parameters
    integer(ip) ,intent(in)    :: ks, nn, nv
    integer(ip) ,intent(in)    :: ln(nn), ia(nv+1), ja(ia(nv+1)-1)
    real(rp)    ,intent(in)    :: ea(nn,nn)
    real(rp)    ,intent(inout) :: a(ia(nv+1)-1)

    ! Locals
    integer(ip)                :: in,jn,iv,jv,iz
    integer(ip)                :: distance
    
    if(ks==1) then ! Unsymmetric 
       do in=1,nn
          iv = ln(in)
          do jn = 1,nn
             jv = ln(jn)

             iz   = ia(iv)
             do while( ja(iz) /= jv )
                iz = iz + 1
             end do

             a(iz) = a(iz) + ea(in, jn)
          end do ! end jn
       end do ! end in
    else ! Symmetric (upper triangle)
       do in=1,nn
          iv=ln(in)
          do jn = 1,nn
             jv = ln(jn)

             ! Note that: 
             ! If (distance>0)      entry belongs to upper triangle
             ! else if (distance<0) entry belongs to lower triangle
             ! else                 entry belongs to main diagonal
             distance = jv - iv

             if ( distance >= 0 ) then
                ! All entries are in the upper triangule
                iz   = ia(iv)
                do while( ja(iz) /= jv )
                   iz = iz + 1
                end do
                a(iz) = a(iz) + ea(in,jn)
             end if
          end do
       end do
    end if

  end subroutine ass_csr_mat_scal

  !============================================================================================
  subroutine mesh_to_graph_matrix ( gtype, primal_mesh, primal_graph )
    implicit none

    ! Parameters
    integer(ip)          , intent(in)     :: gtype
    type(mesh_t)       , intent(in)     :: primal_mesh
    type(graph_t)      , intent(out)    :: primal_graph


    ! Local variables
    type (mesh_t)                        :: dual_mesh
    integer(ip), allocatable               :: iwork(:)       ! Integer ip working array
    integer(ip)                            :: pwork(3)       ! Pointers to work space

    ! Valid graph_t types for addressing sparse matrices are csr_symm or csr
    assert( gtype==csr_symm .or. gtype==csr )


    ! Compute dual_mesh
    call mesh_to_dual ( primal_mesh, dual_mesh )

    ! Allocate space for ia on the primal graph
    primal_graph%type = gtype

    primal_graph%nv  = primal_mesh%npoin
    primal_graph%nv2 = primal_mesh%npoin
    call memalloc ( primal_graph%nv+1, primal_graph%ia, __FILE__,__LINE__ )
       
    ! Allocate working space for count_primal_graph and list_primal_graph routines
    ! (TOTAL WS SIZE = primal mesh npoin + maximum number of neighbours of any primal 
    ! graph node)
    pwork(1) = 1
    pwork(2) = pwork(1) + primal_mesh%npoin
    pwork(3) = pwork(2) + dual_mesh%nnode*primal_mesh%nnode

    call memalloc ( pwork(3), iwork, __FILE__,__LINE__ )
    
    call count_primal_graph_csr_scal ( primal_graph%type, primal_mesh, dual_mesh, primal_graph, &  
                                       iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )

    ! Allocate space for ja on the primal graph 
    call memalloc ( primal_graph%ia(primal_graph%nv+1)-1, primal_graph%ja, __FILE__,__LINE__)
     
    call list_primal_graph_csr_scal  ( primal_graph%type, primal_mesh, dual_mesh, primal_graph, &
                                       iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)) )
    
    ! Free dual_mesh
    call mesh_free ( dual_mesh )
    call memfree ( iwork,__FILE__,__LINE__)
    
  end subroutine mesh_to_graph_matrix

  !============================================================================================
  subroutine  count_primal_graph_csr_scal ( gtype, primal_mesh, dual_mesh, primal_graph,  &
       &                                    ws_position, ws_neighbors )
    implicit none

    ! Parameters
    integer(ip)    , intent(in)     :: gtype
    type(mesh_t) , intent(in)     :: primal_mesh, dual_mesh
    type(graph_t), intent(inout)  :: primal_graph
    integer(ip)    , intent(out)    :: ws_position  (primal_mesh%npoin)
    integer(ip)    , intent(out)    :: ws_neighbors (primal_mesh%nnode*dual_mesh%nnode)

    ! Local variables
    integer(ip)                     :: ineigh  
    integer(ip)                     :: ipoindm ! Dual   mesh point identifier
    integer(ip)                     :: ipoinpm ! Primal mesh point identifier
    integer(ip)                     :: ipoinpg ! Primal graph point identifier 
    integer(ip)                     :: p_ipoindm
    integer(ip)                     :: p_ipoinpm
    integer(ip)                     :: inods1,inods2   

    integer(ip)                     :: first_free_pos
    integer(ip)                     :: num_neighbors

    integer(ip)                     :: npoin, idof1

    ! Initialize work space of filled-up positions
    ws_position   = 0

    ! Initialize work space of neighbour identifiers 
    ws_neighbors  = 0

    primal_graph%ia(1) = 1

    do ipoinpg=1, primal_graph%nv
       primal_graph%ia(ipoinpg+1) = primal_graph%ia(ipoinpg)

       first_free_pos = 1
       num_neighbors  = 0

       ! Traverse the dual nodes of the dual element number ipoinpg
       ! (i.e., primal elements around primal node number ipoinpg)
       do p_ipoindm = dual_mesh%pnods(ipoinpg), dual_mesh%pnods(ipoinpg+1)-1
          ipoindm = dual_mesh%lnods(p_ipoindm)

          ! Traverse the primal nodes of the primal element number ipoindm
          ! (i.e., dual elements around dual node number ipoindm)
          inods1=primal_mesh%pnods(ipoindm)
          inods2=primal_mesh%pnods(ipoindm+1)-1
          do p_ipoinpm = inods1,inods2
             ipoinpm = primal_mesh%lnods(p_ipoinpm)
             if ( gtype == csr .or. ( gtype == csr_symm .and. ipoinpg <= ipoinpm ) ) then
                ! write (*,*) ipoinpg, ipoindm, ipoinpm ! DBG: 
                ! If ipoinpm not visited yet
                if ( ws_position(ipoinpm) == 0 ) then
                   ws_position(ipoinpm) = first_free_pos
                   ws_neighbors(first_free_pos) = ipoinpm
                   first_free_pos = first_free_pos + 1
                   num_neighbors = num_neighbors + 1 
                   ! If ipoinpm already visited
                end if
             end if
          end do
       end do
          
       ! Extract neighbours while restoring working space to initial state
       do ineigh=1, num_neighbors
          primal_graph%ia(ipoinpg+1) = primal_graph%ia(ipoinpg+1) + 1
             ws_position ( ws_neighbors(ineigh) ) = 0
             ws_neighbors(ineigh) = 0
          end do
       end do

     end subroutine count_primal_graph_csr_scal

     !============================================================================================
     subroutine  list_primal_graph_csr_scal ( gtype, primal_mesh, dual_mesh, primal_graph,  &
          &                                   ws_position, ws_neighbors )    
       implicit none

       ! Parameters
       integer(ip)    , intent(in)    :: gtype
       type(mesh_t) , intent(in)    :: primal_mesh, dual_mesh
       type(graph_t), intent(inout) :: primal_graph
       integer(ip)    , intent(out)   :: ws_position  (primal_mesh%npoin)
       integer(ip)    , intent(out)   :: ws_neighbors (primal_mesh%nnode*dual_mesh%nnode)

       ! Local variables
       integer(ip)                    :: ineigh  
       integer(ip)                    :: ipoindm ! Dual   mesh point identifier
       integer(ip)                    :: ipoinpm ! Primal mesh point identifier
       integer(ip)                    :: ipoinpg ! Primal graph point identifier 
       integer(ip)                    :: p_ipoindm   
       integer(ip)                    :: p_ipoinpm   
       integer(ip)                    :: inods1,inods2   

       integer(ip)                    :: first_free_ja_pos
       integer(ip)                    :: first_free_pos
       integer(ip)                    :: num_neighbors, idof1, idof2

       ! Initialize work space of filled-up positions to zeros
       ws_position   = 0

       ! Initialize working space of neighbour identifiers 
       ws_neighbors = 0

       first_free_ja_pos = 1

       do ipoinpg=1, primal_graph%nv     
          first_free_pos = 1
          num_neighbors  = 0

          ! Traverse the dual nodes of the dual element number ipoinpg
          ! (i.e., primal elements around primal node number ipoinpg)
          do p_ipoindm = dual_mesh%pnods(ipoinpg), dual_mesh%pnods(ipoinpg+1)-1
             ipoindm = dual_mesh%lnods(p_ipoindm)

             ! Traverse the primal nodes of the primal element number ipoindm
             ! (i.e., dual elements around dual node number ipoindm)
             inods1=primal_mesh%pnods(ipoindm)
             inods2=primal_mesh%pnods(ipoindm+1)-1
             do p_ipoinpm = inods1,inods2
                ipoinpm = primal_mesh%lnods(p_ipoinpm)
                ! Only list edges (i,j) s.t., i <= j  
                if ( gtype == csr .or. ( gtype == csr_symm .and. ipoinpg <= ipoinpm ) ) then
                   ! If ipoinpm not visited yet
                   if ( ws_position(ipoinpm) == 0 ) then
                      ws_position(ipoinpm) = first_free_pos
                      ws_neighbors(first_free_pos) = ipoinpm
                      first_free_pos = first_free_pos + 1
                      num_neighbors = num_neighbors + 1 
                      ! If ipoinpm already visited
                   end if
                end if
             end do
          end do

          ! Extract neighbours while restoring working space to initial state
          do ineigh=1, num_neighbors
             primal_graph%ja(first_free_ja_pos) = ws_neighbors(ineigh)
             first_free_ja_pos = first_free_ja_pos + 1
             
             ws_position ( ws_neighbors(ineigh) ) = 0
             ws_neighbors( ineigh               ) = 0
          end do

          ! Order increasingly column identifiers of current row 
          ! using heap sort algorithm
          ! write (*,*) 'A', primal_graph%ja( primal_graph%ia(ipoinpg):primal_graph%ia(ipoinpg+1)-1 ) ! DBG:
          call sort(primal_graph%ia(ipoinpg+1)-primal_graph%ia(ipoinpg), & 
                    primal_graph%ja(primal_graph%ia(ipoinpg):primal_graph%ia(ipoinpg+1)-1))
          ! write (*,*) 'D', primal_graph%ja( primal_graph%ia(ipoinpg):primal_graph%ia(ipoinpg+1)-1 ) ! DBG:

       end do

     end subroutine list_primal_graph_csr_scal

     !=============================================================================
     subroutine par_preconditioner_dd_mlevel_bddc_apply_tbp (op, x, y)
       implicit none
       ! Parameters
       class(par_preconditioner_dd_mlevel_bddc_t)    , intent(in)    :: op
       class(base_operand_t)   , intent(in)    :: x
       class(base_operand_t)   , intent(inout) :: y

!!$       assert (associated(op%mat))

       call x%GuardTemp()

       select type(x)
          class is (par_vector_t)
          select type(y)
             class is(par_vector_t)
             call par_preconditioner_dd_mlevel_bddc_apply_all_unk ( op, x, y )
             class default
             write(0,'(a)') 'matrix_t%apply: unsupported y class'
             check(1==0)
          end select
          class default
          write(0,'(a)') 'par_preconditioner_dd_mlevel_bddc_t%apply: unsupported x class'
          check(1==0)
       end select

       call x%CleanTemp()
     end subroutine par_preconditioner_dd_mlevel_bddc_apply_tbp


     !=============================================================================
     function par_preconditioner_dd_mlevel_bddc_apply_fun_tbp (op, x) result(y)
       implicit none
       ! Parameters
       class(par_preconditioner_dd_mlevel_bddc_t), intent(in)   :: op
       class(base_operand_t), intent(in)  :: x
       class(base_operand_t), allocatable :: y
       type(par_vector_t), allocatable :: local_y

       call x%GuardTemp()

       select type(x)
          class is (par_vector_t)
          allocate(local_y)
          call par_vector_alloc ( x%dof_dist, x%p_env, local_y)
          call par_preconditioner_dd_mlevel_bddc_apply_all_unk ( op, x, local_y )
          call move_alloc(local_y, y)
          call y%SetTemp()
          class default
          write(0,'(a)') 'par_preconditioner_dd_mlevel_bddc_t%apply_fun: unsupported x class'
          check(1==0)
       end select

       call x%CleanTemp()
     end function par_preconditioner_dd_mlevel_bddc_apply_fun_tbp

     !=============================================================================
     subroutine par_preconditioner_dd_mlevel_bddc_fill_values_tbp (op)
       implicit none
       ! Parameters
       class(par_preconditioner_dd_mlevel_bddc_t)    , intent(inout) :: op

       assert (associated(op%p_mat))

       if(op%do_fill_values) call par_preconditioner_dd_mlevel_bddc_fill_val ( op )

     end subroutine par_preconditioner_dd_mlevel_bddc_fill_values_tbp

        subroutine par_preconditioner_dd_mlevel_bddc_free_tbp(this)
          implicit none
          class(par_preconditioner_dd_mlevel_bddc_t), intent(inout) :: this
        end subroutine par_preconditioner_dd_mlevel_bddc_free_tbp

  subroutine matrix_preconditioner_vector_solve (A,M,b,x,pars)
    implicit none
    type(matrix_t)    ,intent(in)    :: A     ! Matrix
    type(preconditioner_t)   ,intent(in)    :: M     ! Preconditioner
    type(vector_t)    ,intent(in)    :: b     ! RHS
    type(vector_t)    ,intent(inout) :: x     ! Approximate solution
    type(solver_control_t),intent(inout) :: pars  ! Solver parameters

    ! Locals
    type(serial_environment_t) :: senv

    call abstract_solve (A, M, b, x, pars, senv)

  end subroutine matrix_preconditioner_vector_solve

  subroutine matrix_preconditioner_r1_solve (A,M,b,x,pars)
    implicit none
    type(matrix_t)            ,intent(in)    :: A          ! Matrix
    type(preconditioner_t)           ,intent(in)    :: M          ! Preconditioner
    real(rp)           , target ,intent(in)    :: b(A%gr%nv) ! RHS
    real(rp)           , target ,intent(inout) :: x(A%gr%nv) ! Approximate solution
    type(solver_control_t)        ,intent(inout) :: pars       ! Solver parameters

    ! Locals
    type(vector_t)         :: vector_b
    type(vector_t)         :: vector_x
    type(serial_environment_t) :: senv

    ! fill vector_b members
    vector_b%neq     =  A%gr%nv
    vector_b%mode    =  reference
    vector_b%b => b   

    ! fill vector_x members
    vector_x%neq     = A%gr%nv 
    vector_x%mode    = reference 
    vector_x%b       => x 

    call abstract_solve (A, M, vector_b, vector_x, pars, senv)

  end subroutine matrix_preconditioner_r1_solve

  subroutine matrix_preconditioner_r2_solve (A,M,b,ldb,x,ldx,pars)
    implicit none
    type(matrix_t)            ,intent(in)    :: A                 ! Matrix
    type(preconditioner_t)           ,intent(in)    :: M                 ! Preconditioner
    type(solver_control_t)        ,intent(inout) :: pars              ! Solver parameters
    integer(ip)                 ,intent(in)    :: ldb, ldx
    real(rp)           , target ,intent(in)    :: b(ldb, pars%nrhs) ! RHS
    real(rp)           , target ,intent(inout) :: x(ldx, pars%nrhs) ! Approximate solution

    ! Locals
    type(vector_t)         :: vector_b
    type(vector_t)         :: vector_x
    type(serial_environment_t) :: senv
    integer(ip)              :: k
    integer(ip)              :: tot_its

!!$    if (pars%method == direct) then
!!$       call preconditioner_apply (A, M, pars%nrhs, vector_b, ldb, vector_x, ldx)
!!$    else
    tot_its = 0 
    do k=1, pars%nrhs
       ! fill b members
       vector_b%neq     =  A%gr%nv
       vector_b%mode    =  reference
       vector_b%b       => b(:,k)
       
       ! fill vector_x members
       vector_x%neq     =  A%gr%nv
       vector_x%mode    =  reference 
       vector_x%b       => x(:,k)

       call abstract_solve (A, M, vector_b, vector_x, pars, senv)
       
       tot_its = tot_its + pars%it
    end do
    pars%it = tot_its
!!$    end if
  end subroutine matrix_preconditioner_r2_solve

end module par_preconditioner_dd_mlevel_bddc_names



