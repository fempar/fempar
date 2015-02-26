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
module par_dd_base
  ! Module that includes the parameters used in the several par_precond_dd's
  use types

  implicit none

  ! Common to BDDC/BNN
  integer(ip), parameter :: f_tasks_w_coarse_duties = 0 ! (1) All MPI Tasks are in p_context communicator .or.
                                                        ! (2) MPI Tasks separated into fine_tasks (p_context) 
                                                        !     and coarse_tasks (c_contxt)
                                                        ! In both cases tasks in p_context communicator responsible
                                                        ! for the solution of the coarse-grid problem

  integer(ip), parameter :: f_tasks_c_tasks_w_coarse_duties = 1 ! MPI Tasks MUST BE separated into  
                                                                ! fine_tasks (p_context) 
                                                                ! and coarse_tasks (c_contxt). The first ones 
                                                                ! are responsible for fine-grid duties, the second ones
                                                                ! for coarse-grid duties 


  integer(ip), parameter :: all_unknowns          = 0  ! Precondition the global linear system Ax=b
  integer(ip), parameter :: interface_unknowns    = 1  ! Precondition interface problem Sy=g


  ! BDDC
  ! Strategy for the computation of the fine-grid correction 
  integer(ip), parameter  :: corners_rest_part_solve_impl_schur = 0
  integer(ip), parameter  :: corners_rest_part_solve_expl_schur = 1
  integer(ip), parameter  :: direct_solve_constrained_problem   = 2

  ! How to combine coarse and fine grid correction ?
  integer(ip), parameter :: only_f_correction  = 0
  integer(ip), parameter :: additive           = 1
  integer(ip), parameter :: multiplicative_cfc = 2 ! Coarse+Fine+Coarse correction on each iteration
  integer(ip), parameter :: multiplicative_fc  = 3 ! Fine+Coarse correction on each iteration. 
                                                   ! Coarse correction at the beginning of Krylov subspace solver
                                                   ! (implementation pending for Galerkin BDDC)
  integer(ip), parameter :: multiplicative_fc_enhanced  = 4 ! Fine+Coarse correction on each iteration. 
                                                            ! Coarse correction at the beginning Krylov subspace solver
                                                            ! Saves one Dirichlet problem
  integer(ip), parameter :: additive_symmetric  = 5 ! Symmetrized version of the global BDDC method required for inexact solvers.

  ! Type of coarse dofs to be considered for
  ! the coarse grid system (corners are sufficient
  ! for 2D problems, while for 3D, both corners and
  ! edges are required)
  integer (ip), parameter :: corners           = 0
  integer (ip), parameter :: corners_and_edges = 1
  integer (ip), parameter :: corners_edges_and_faces = 2 


  integer (ip), parameter :: naive     =  0 ! Corners are degenerated edges. Does not allow 
                                            ! to select corners as nodes belonging to edges/faces
                                            ! No comms. are required. May fail if there are not
                                            ! enough corners to properly get rid of the kernel of the 
                                            ! Neumman operator

  integer (ip), parameter :: laplacian =  1 ! Corners are degenerated edges or nodes belonging to
                                            ! edges/faces in such a way that there is at least 1
                                            ! corner per edge/face in 2D/3D, resp.
                                            ! Comms. among nearest neighbours are required.

  integer (ip), parameter :: elasticity = 2 ! Corners are degenerated edges and/or nodes belonging to
                                            ! edges/faces in such a way that there are at least 2/3
                                            ! corners per edge/face in 2D/3D, resp.
                                            ! Comms. among nearest neighbours are required.

  ! Strategy for the solution of the coarse grid system
  integer (ip), parameter :: serial_gather     = 0 ! Only one processor gathers and solves
                                                   ! serially the coarse grid system by means
                                                   ! of a direct sparse solver (PARDISO)

  integer (ip), parameter :: serial_all_gather = 1 ! All processors in the communicator gather
                                                   ! and solve serially the coarse grid system
                                                   ! by means of a direct sparse solver (PARDISO)

  integer (ip), parameter :: distributed       = 2 ! All processors in the communicatior solve
                                                   ! the coarse grid system by means of a distributed
                                                   ! direct sparse solver (e.g., MUMPS; not implemented 
                                                   ! yet at all)

  integer (ip), parameter :: recursive_bddc   = 3  ! A recursive call to bddc is performed (multilevel)

  integer    , parameter :: root_pid           = 0 ! In case of serial_gather, the processor id in charge
                                                   ! of the global information (internal/private constant)

  integer (ip), parameter :: petrov_galerkin    = 0 ! Petrov-galerkin projection is considered(nonsymmetric)
  integer (ip), parameter :: galerkin           = 1 ! Galerkin projection for both test and solution function

  ! BNN

  ! Which rigid body motions to include in the coarse space ?
  integer(ip), parameter :: translation              = 0
  integer(ip), parameter :: translation_and_rotation = 1

  ! Strategy for the solution of the neumann problem
  integer(ip), parameter  :: impl_sol_aug_with_z_in_kernel = 0
  integer(ip), parameter  :: impl_sol_aug_with_z_mapped    = 1
  integer(ip), parameter  :: expl_sol_with_z_mapped        = 2

  ! Contrain all unks (i.e., interior+interface) or only interface
  ! unks in the computation of the fine grid correction 
  integer(ip), parameter :: constrain_all_unks            = 0
  integer(ip), parameter :: constrain_only_interface_unks = 1

  ! How to treat non-floating subdomains fine-grid correction ?
  integer (ip), parameter :: treat_as_floating    = 0 
  integer (ip), parameter :: treat_wo_constrain   = 1

  ! Pad collectives ?
  integer(ip), parameter :: nopad  = 0
  integer(ip), parameter :: pad    = 1

  ! Re-use Schur complement edge-lagrange multipliers and  A_rr^-1 C_r^T
  ! for the computation of the neumann problem or compute them
  ! from scratch using spars_neumann
  integer(ip), parameter  :: reuse_from_phis      = 0
  integer(ip), parameter  :: compute_from_scratch = 1
  integer(ip), parameter  :: compute_as_Y_T_A_Y   = 2 ! Where Y=A_rr^{-1} C_r^T
                                                       ! Javier's manuscript notes
                                                       ! to symmetrize Var. 3 (art008)
                                                       ! wo computing from scratch 

  ! How to calculate subdomain element matrices (subd_elmat_calc)?
  integer(ip), parameter :: phit_a_i_phi            = 0   ! \Phi_t A_i \Phi
  integer(ip), parameter :: phit_minus_c_i_t_lambda = 1   ! \Phi_t (-C_i^t \Lambda)

  ! Who is responsible for the set-up preconditioners
  ! for the internal problems?
  integer(ip), parameter :: handled_by_bddc_module = 0 
  integer(ip), parameter :: handled_by_user_driver = 1 

  ! Temporize coarse-grid computations in isolation for f_tasks_w_coarse_duties
  ! model of execution in case of the BDDC. Temporize coarse-grid computations
  ! in isolation in case of the BNN method 
  logical(lg) , parameter :: temp_coarse_grid_computations = .false.

  ! Temporize fine/coarse-grid overlap for f_tasks_c_tasks_w_coarse_duties model
  ! of execution in case of the BDDC. This does not apply to the BNN method as
  ! we have not still implemented f_tasks_c_tasks_w_coarse_duties model of execution
  ! for the BNN method
  logical(lg) , parameter :: temp_fine_coarse_grid_overlap  = .false.


  integer (ip), parameter :: nokernel  =  0 ! The local Neumann operator has no kernel
 
  ! (Already defined above for the BDDC method)
!!$  integer (ip), parameter :: laplacian =  1 ! The local Neumann operator has the kernel
!!$                                            ! corresponding to the Laplacian PDE 
!!$
!!$  integer (ip), parameter :: elasticity = 2 ! The local Neumann operator has the kernel
!!$                                            ! corresponding to the Eslasticity PDE


  ! Domain decomposition Preconditioners supported
  integer(ip), parameter :: dd_prec_none        = 0 ! Identity preconditioner
  integer(ip), parameter :: dd_prec_diag        = 1 ! D^-1, where D=diag(A)
  integer(ip), parameter :: dd_prec_nn          = 2
  integer(ip), parameter :: dd_prec_mlevel_bddc = 5


  integer(ip), parameter :: dd_default_preconditioner = dd_prec_none

end module par_dd_base
