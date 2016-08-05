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

module mlbddc_names
 ! Tools
 use types_names
 use list_types_names
 use allocatable_array_names
 
 ! Integration related modules
 use base_static_triangulation_names
 use new_serial_fe_space_names
 use new_fe_affine_operator_names
 
 ! Linear Algebra related modules
 use operator_names
 use par_scalar_array_names
 use serial_scalar_array_names
 use vector_space_names
 use vector_names
 
 use matrix_names
 use base_sparse_matrix_names
 use sparse_matrix_parameters_names
 use sparse_matrix_names
 use par_sparse_matrix_names
 use direct_solver_names
 use direct_solver_parameters_names

#ifdef ENABLE_BLAS
 use blas77_interfaces_names
#endif

#ifdef ENABLE_LAPACK
 use lapack77_interfaces_names
#endif
  
 ! Parallel communication-related data structures
 use par_environment_names
 
 implicit none
# include "debug.i90"
 private

 type, extends(fe_object_accessor_t) :: dof_object_accessor_t
   private
   type(mlbddc_t), pointer :: mlbddc
   integer(ip)             :: field_id
   integer(ip)             :: dof_object_lid
 contains
    procedure  :: dof_object_accessor_create
    generic    :: create                      => dof_object_accessor_create
    procedure  :: free                        => dof_object_accessor_free
    procedure  :: next                        => dof_object_accessor_next
    procedure  :: set_lid                     => dof_object_accessor_set_lid
    procedure  :: get_number_dofs_on_object   => dof_object_accessor_get_number_dofs_on_object
    procedure  :: get_dofs_on_object_iterator => dof_object_accessor_get_dofs_on_object_iterator
 end type dof_object_accessor_t
 
 type dof_object_iterator_t
    private
    type(dof_object_accessor_t) :: current_dof_object_accessor
  contains
     procedure, non_overridable          :: create       => dof_object_iterator_create
     procedure, non_overridable          :: free         => dof_object_iterator_free
     procedure, non_overridable          :: init         => dof_object_iterator_init
     procedure, non_overridable          :: next         => dof_object_iterator_next
     procedure, non_overridable          :: has_finished => dof_object_iterator_has_finished
     procedure, non_overridable          :: current      => dof_object_iterator_current
  end type dof_object_iterator_t
 
 type, extends(operator_t) :: mlbddc_t
   private
 
   ! Pointer to the fe_affine_operator_t this mlbddc_t instance has been created from
   type(new_fe_affine_operator_t)    , pointer :: fe_affine_operator => NULL()
   
   !**** BEGIN member variables which control coarse DoFs on top of coarse VEFs ****!
   !********************************************************************************!
   
   ! IMPORTANT NOTE: The member variables encompassed within this commented region 
   !                 were originally member variables of type(par_fe_space_t). I 
   !                 think that having these out of type(par_fe_space_t) leads to a
   !                 much more extensible sw design. Recall that these member variables 
   !                 are created/filled as a result of a process that is to be customized 
   !                 by the user. Besides, we might have different type(mlbddc_t) instances 
   !                 living together within the same program s.t. each of these instances 
   !                 need a different configuration for these member variables. Having these 
   !                 within type(par_fe_space_t) clearly prevents such use cases. Other 
   !                 issues that have to be addressed and for which I currently do not have a 
   !                 clear answer are: (1) who is in charge of filling these member variables? 
   !                 (2) Are these member variables in their current state sufficient/general 
   !                 enough to represent in the computer the result of any coarse DoF 
   !                 identification process? (3) Is this the most appropriate place to place 
   !                 them?
   
   ! Aggregation of fine DoFs into DoF objects
   integer(ip)                   , allocatable :: num_dofs_objects_per_field(:)
   type(list_t)                  , allocatable :: dofs_objects_per_field(:) 
   
   ! GIDs of the DoFs objects. For their generation, we though to be a good idea 
   ! to take as basis the GIDs of the VEFs objects they are built from (instead of
   ! generating them from scratch, which in turn would imply further communication).
   type(allocatable_array_igp1_t), allocatable :: dofs_objects_gids_per_field(:)
     
   ! GIDs of the VEFs objects the DoFs objects are put on top of
   type(allocatable_array_ip1_t) , allocatable :: vefs_lids_dofs_objects_per_field(:)
   
   type(coo_sparse_matrix_t)                   :: constraint_matrix
   !********************************************************************************!
   !**** END member variables which control coarse DoFs on top of coarse VEFs ****!
 
   ! Constrained Neumann problem-related member variables
   ! B = [ A C^T ]
   !     [ C   0 ]
   type(direct_solver_t)                       :: constrained_neumann_solver
   type(sparse_matrix_t)                       :: constrained_neumann_matrix
   
   ! Dirichlet problem-related member variables
   ! A => [ A_II A_IG ]
   !      [ A_GI A_GG ]
   type(direct_solver_t)                       :: dirichlet_solver
   type(sparse_matrix_t)                       :: A_II
   type(sparse_matrix_t)                       :: A_IG
   type(sparse_matrix_t)                       :: A_GI
   type(sparse_matrix_t)                       :: A_GG
   
   ! Coarse-grid problem related member variables
   real(rp), allocatable                       :: phi(:,:)
   
   ! Pointer to data structure which is in charge of coarse DoF handling.
   ! It will be a nullified pointer on L1 tasks, and associated via target 
   ! allocation in the case of L2-Ln tasks.
   type(coarse_fe_space_t)       , pointer     :: coarse_fe_space => NULL()
   
   ! Coarse-grid matrix. It is temporarily stored in a type(par_sparse_matrix_t)
   ! data structure, although, in my opinion, in the seek of extensibility, 
   ! some sort of operator is required here that plays the same role as
   ! type(fe_affine_operator_t) on L1 tasks. It will be a nullified pointer on 
   ! L1 tasks, and associated via target allocation in the case of L2-Ln tasks.
   type(direct_solver_t)                       :: coarse_solver
   type(par_sparse_matrix_t)     , pointer     :: coarse_grid_matrix => NULL()
   
   ! Next level in the preconditioning hierarchy. It will be a nullified pointer on 
   ! L1 tasks, and associated via target allocation in the case of L2-Ln tasks.
   type(mlbddc_coarse_t)         , pointer     :: mlbddc_coarse   => NULL()
 contains
    procedure, non_overridable          :: create                                          => mlbddc_create
    procedure, non_overridable, private :: create_vector_spaces                            => mlbddc_create_vector_spaces
    
    ! Symbolic setup-related TBPs
    procedure, non_overridable          :: symbolic_setup                                  => mlbddc_symbolic_setup
    procedure, non_overridable, private :: setup_dofs_objects_and_constraint_matrix        => mlbddc_setup_dofs_objects_and_constraint_matrix
    procedure, non_overridable, private :: setup_coarse_fe_space                           => mlbddc_setup_coarse_fe_space
    procedure, non_overridable, private :: transfer_number_fields                          => mlbddc_transfer_number_fields    
    procedure, non_overridable, private :: transfer_field_type                             => mlbddc_transfer_field_type
    procedure, non_overridable, private :: gather_ptr_dofs_per_fe_and_field                => mlbddc_gather_ptr_dofs_per_fe_and_field
    procedure, non_overridable, private :: gather_coarse_dofs_gids_rcv_counts_and_displs   => mlbddc_gather_coarse_dofs_gids_rcv_counts_and_displs
    procedure, non_overridable, private :: gather_coarse_dofs_gids                         => mlbddc_gather_coarse_dofs_gids
    procedure, non_overridable, private :: gather_vefs_gids_dofs_objects                   => mlbddc_gather_vefs_gids_dofs_objects
    procedure, non_overridable, private :: symbolic_setup_dirichlet_problem                => mlbddc_symbolic_setup_dirichlet_problem
    procedure, non_overridable, private :: symbolic_setup_dirichlet_solver                 => mlbddc_symbolic_setup_dirichlet_solver
    procedure, non_overridable, private :: symbolic_setup_constrained_neumann_problem      => mlbddc_symbolic_setup_constrained_neumann_problem
    procedure, non_overridable, private :: symbolic_setup_constrained_neumann_solver       => mlbddc_symbolic_setup_constrained_neumann_solver
    procedure, non_overridable, private :: symbolic_setup_coarse_grid_matrix               => mlbddc_symbolic_setup_coarse_grid_matrix
    procedure, non_overridable, private :: coarse_grid_matrix_symbolic_assembly            => mlbddc_coarse_grid_matrix_symbolic_assembly 
    procedure, non_overridable, private :: symbolic_setup_mlbddc_coarse                    => mlbddc_symbolic_setup_mlbddc_coarse
    procedure, non_overridable, private :: symbolic_setup_coarse_solver                    => mlbddc_symbolic_setup_coarse_solver
    
    ! Numerical setup-related TBPs
    procedure, non_overridable          :: numerical_setup                                 => mlbddc_numerical_setup
    procedure, non_overridable, private :: numerical_setup_dirichlet_problem               => mlbddc_numerical_setup_dirichlet_problem
    procedure, non_overridable, private :: numerical_setup_dirichlet_solver                => mlbddc_numerical_setup_dirichlet_solver
    procedure, non_overridable, private :: numerical_setup_constrained_neumann_problem     => mlbddc_numerical_constrained_neumann_problem
    procedure, non_overridable, private :: numerical_setup_constrained_neumann_solver      => mlbddc_numerical_setup_constrained_neumann_solver
    procedure, non_overridable, private :: allocate_coarse_grid_basis                      => mlbddc_allocate_coarse_grid_basis
    procedure, non_overridable, private :: setup_coarse_grid_basis                         => mlbddc_setup_coarse_grid_basis
    procedure, non_overridable, private :: compute_subdomain_elmat                         => mlbddc_compute_subdomain_elmat
    procedure, non_overridable, private :: compute_and_gather_subdomain_elmat              => mlbddc_compute_and_gather_subdomain_elmat
    procedure, non_overridable, private :: compute_subdomain_elmat_counts_and_displs       => mlbddc_compute_subdomain_elmat_counts_and_displs
    procedure, non_overridable, private :: numerical_setup_coarse_grid_matrix              => mlbddc_numerical_setup_coarse_grid_matrix
    procedure, non_overridable, private :: coarse_grid_matrix_numerical_assembly           => mlbddc_coarse_grid_matrix_numerical_assembly 
    procedure, non_overridable, private :: numerical_setup_mlbddc_coarse                   => mlbddc_numerical_setup_mlbddc_coarse
    procedure, non_overridable, private :: numerical_setup_coarse_solver                   => mlbddc_numerical_setup_coarse_solver
    
    ! Apply related TBPs
    procedure                           :: apply                                           => mlbddc_apply
    procedure, non_overridable, private :: apply_par_scalar_array                          => mlbddc_apply_par_scalar_array
    procedure, non_overridable, private :: solve_coarse_problem                            => mlbddc_solve_coarse_problem
    procedure, non_overridable, private :: compute_coarse_correction                       => mlbddc_compute_coarse_correction
    procedure, non_overridable, private :: setup_coarse_grid_residual                      => mlbddc_setup_coarse_grid_residual
    procedure, non_overridable, private :: compute_coarse_dofs_values                      => mlbddc_compute_coarse_dofs_values
    procedure, non_overridable, private :: compute_and_gather_coarse_dofs_values           => mlbddc_compute_and_gather_coarse_dofs_values
    procedure, non_overridable, private :: compute_coarse_dofs_values_counts_and_displs    => mlbddc_compute_coarse_dofs_values_counts_and_displs
    procedure, non_overridable, private :: coarse_grid_residual_assembly                   => mlbddc_coarse_grid_residual_assembly
    procedure, non_overridable, private :: scatter_and_interpolate_coarse_grid_correction  => mlbddc_scatter_and_interpolate_coarse_grid_correction
    procedure, non_overridable, private :: scatter_coarse_grid_correction                  => mlbddc_scatter_coarse_grid_correction
    procedure, non_overridable, private :: fill_coarse_dofs_values_scattered               => mlbddc_fill_coarse_dofs_values_scattered
    procedure, non_overridable, private :: interpolate_coarse_grid_correction              => mlbddc_interpolate_coarse_grid_correction
    procedure, non_overridable, private :: solve_dirichlet_problem                         => mlbddc_solve_dirichlet_problem
    procedure, non_overridable, private :: apply_A_GI                                      => mlbddc_apply_A_GI
    procedure, non_overridable, private :: apply_A_IG                                      => mlbddc_apply_A_IG
    procedure, non_overridable, private :: solve_constrained_neumann_problem               => mlbddc_solve_constrained_neumann_problem
    procedure, non_overridable, private :: apply_weight_operator                           => mlbddc_apply_weight_operator
    procedure, non_overridable, private :: create_interior_interface_views                 => mlbddc_create_interior_interface_views
    

    ! Free-related TBPs
    procedure, non_overridable          :: free                                             => mlbddc_free
    procedure, non_overridable          :: free_clean                                       => mlbddc_free_clean
    procedure, non_overridable          :: free_symbolic_setup                              => mlbddc_free_symbolic_setup
    procedure, non_overridable          :: free_dofs_objects_and_constraint_matrix          => mlbddc_free_dofs_objects_and_constraint_matrix
    procedure, non_overridable          :: free_symbolic_setup_dirichlet_problem            => mlbddc_free_symbolic_setup_dirichlet_problem
    procedure, non_overridable          :: free_symbolic_setup_dirichlet_solver             => mlbddc_free_symbolic_setup_dirichlet_solver
    procedure, non_overridable          :: free_symbolic_setup_constrained_neumann_problem  => mlbddc_free_symbolic_setup_constrained_neumann_problem
    procedure, non_overridable          :: free_symbolic_setup_constrained_neumann_solver   => mlbddc_free_symbolic_setup_constrained_neumann_solver
    procedure, non_overridable          :: free_symbolic_setup_coarse_solver                => mlbddc_free_symbolic_setup_coarse_solver
    
    procedure, non_overridable          :: free_numerical_setup                             => mlbddc_free_numerical_setup
    procedure, non_overridable          :: free_numerical_setup_dirichlet_problem           => mlbddc_free_numerical_setup_dirichlet_problem
    procedure, non_overridable          :: free_numerical_setup_dirichlet_solver            => mlbddc_free_numerical_setup_dirichlet_solver
    procedure, non_overridable          :: free_numerical_setup_constrained_neumann_problem => mlbddc_free_numerical_setup_constrained_neumann_problem
    procedure, non_overridable          :: free_numerical_setup_constrained_neumann_solver  => mlbddc_free_numerical_setup_constrained_neumann_solver
    procedure, non_overridable          :: free_coarse_grid_basis                           => mlbddc_free_coarse_grid_basis
    procedure, non_overridable          :: free_numerical_setup_coarse_solver               => mlbddc_free_numerical_setup_coarse_solver
    
    ! Miscellaneous 
    procedure, non_overridable, private :: get_total_number_coarse_dofs                     => mlbddc_get_total_number_coarse_dofs
    procedure, non_overridable, private :: get_block_number_coarse_dofs                     => mlbddc_get_block_number_coarse_dofs
    procedure, non_overridable          :: create_field_dofs_object_iterator                => mlbddc_create_field_dofs_object_iterator
    procedure, non_overridable, private :: get_par_sparse_matrix                            => mlbddc_get_par_sparse_matrix
    procedure, non_overridable, private :: get_fe_space                                     => mlbddc_get_fe_space
    procedure, non_overridable, private :: get_par_environment                              => mlbddc_get_par_environment
    procedure, non_overridable, private :: am_i_l1_task                                     => mlbddc_am_i_l1_task
    procedure                           :: is_linear                                        => mlbddc_is_linear
 end type mlbddc_t
 
 type, extends(coarse_fe_object_accessor_t) :: coarse_dof_object_accessor_t
   private
   type(mlbddc_coarse_t), pointer :: mlbddc_coarse
   integer(ip)                    :: field_id
   integer(ip)                    :: coarse_dof_object_lid
 contains
    procedure  :: coarse_dof_object_accessor_create
    generic    :: create                      => coarse_dof_object_accessor_create
    procedure  :: free                        => coarse_dof_object_accessor_free
    procedure  :: next                        => coarse_dof_object_accessor_next
    procedure  :: set_lid                     => coarse_dof_object_accessor_set_lid
    procedure  :: get_number_dofs_on_object   => coarse_dof_object_accessor_get_number_dofs_on_object
    procedure  :: get_dofs_on_object_iterator => coarse_dof_object_accessor_get_dofs_on_object_iterator
 end type coarse_dof_object_accessor_t
 
 type coarse_dof_object_iterator_t
    private
    type(coarse_dof_object_accessor_t) :: current_coarse_dof_object_accessor
  contains
     procedure, non_overridable          :: create       => coarse_dof_object_iterator_create
     procedure, non_overridable          :: free         => coarse_dof_object_iterator_free
     procedure, non_overridable          :: init         => coarse_dof_object_iterator_init
     procedure, non_overridable          :: next         => coarse_dof_object_iterator_next
     procedure, non_overridable          :: has_finished => coarse_dof_object_iterator_has_finished
     procedure, non_overridable          :: current      => coarse_dof_object_iterator_current
  end type coarse_dof_object_iterator_t
 
 
 
 type :: mlbddc_coarse_t 
   private
   ! Some sort of operator is required here that plays the role of type(mlbddc_t)%fe_affine_operator.
   ! This operator should be built on the previous level and passed here. Let us use the 
   ! coarse_fe_space built on the previous level in the mean time.
   type(coarse_fe_space_t)       , pointer     :: fe_space          => NULL()
   type(par_sparse_matrix_t)     , pointer     :: par_sparse_matrix => NULL()
   
   !**** BEGIN member variables which control coarse DoFs on top of coarse VEFs ****!
   !********************************************************************************!
   ! IMPORTANT NOTE: The member variables encompassed within this commented region 
   !                 were originally member variables of type(coarse_fe_space_t).
   ! Aggregation of fine DoFs into DoF objects
   integer(ip)                   , allocatable :: num_dofs_objects_per_field(:)
   type(list_t)                  , allocatable :: dofs_objects_per_field(:) 
   ! GIDs of the DoFs objects. For their generation, we though to be a good idea 
   ! to take as basis the GIDs of the VEFs objects they are built from (instead of
   ! generating them from scratch, which in turn would imply further communication).
   type(allocatable_array_igp1_t), allocatable :: dofs_objects_gids_per_field(:)
   ! GIDs of the VEFs objects the DoFs objects are put on top of
   type(allocatable_array_ip1_t) , allocatable :: vefs_lids_dofs_objects_per_field(:)
   type(coo_sparse_matrix_t)                   :: constraint_matrix
   !********************************************************************************!
   !**** END member variables which control coarse DoFs on top of coarse VEFs ****!
   
   
   ! Constrained Neumann problem-related member variables
   ! B = [ A C^T ]
   !     [ C   0 ]
   type(direct_solver_t)                       :: constrained_neumann_solver
   type(sparse_matrix_t)                       :: constrained_neumann_matrix
   
   ! Dirichlet problem-related member variables
   ! A => [ A_II A_IG ]
   !      [ A_GI A_GG ]
   type(direct_solver_t)                       :: dirichlet_solver
   type(sparse_matrix_t)                       :: A_II
   type(sparse_matrix_t)                       :: A_IG
   type(sparse_matrix_t)                       :: A_GI
   type(sparse_matrix_t)                       :: A_GG
   
   ! Coarse-grid problem related member variables
   real(rp), allocatable                       :: phi(:,:)
   
   ! Pointer to data structure which is in charge of coarse DoF handling.
   ! It will be a nullified pointer on L1 tasks, and associated via target 
   ! allocation in the case of L2-Ln tasks.
   type(coarse_fe_space_t)       , pointer     :: coarse_fe_space => NULL()
   
   ! Coarse-grid matrix. It is temporarily stored in a type(par_sparse_matrix_t)
   ! data structure, although, in my opinion, in the seek of extensibility, 
   ! some sort of operator is required here that plays the same role as
   ! type(fe_affine_operator_t) on L1 tasks. It will be a nullified pointer on 
   ! L1 tasks, and associated via target allocation in the case of L2-Ln tasks.
   type(direct_solver_t)                       :: coarse_solver
   type(par_sparse_matrix_t)     , pointer     :: coarse_grid_matrix => NULL()
   
   ! Next level in the preconditioning hierarchy. It will be a nullified pointer on 
   ! L1 tasks, and associated via target allocation in the case of L2-Ln tasks.
   type(mlbddc_coarse_t)         , pointer     :: mlbddc_coarse   => NULL()
 contains
 
   procedure, non_overridable          :: create                                            => mlbddc_coarse_create
   !
   ! Symbolic setup-related TBPs
   procedure, non_overridable          :: symbolic_setup                                    => mlbddc_coarse_symbolic_setup
   procedure, non_overridable, private :: setup_dofs_objects_and_constraint_matrix          => mlbddc_coarse_setup_dofs_objects_and_constraint_matrix
   procedure, non_overridable, private :: setup_coarse_fe_space                             => mlbddc_coarse_setup_coarse_fe_space
   procedure, non_overridable, private :: transfer_number_fields                            => mlbddc_coarse_transfer_number_fields
   procedure, non_overridable, private :: transfer_field_type                               => mlbddc_coarse_transfer_field_type
   procedure, non_overridable, private :: gather_ptr_dofs_per_fe_and_field                  => mlbddc_coarse_gather_ptr_dofs_per_fe_and_field
   procedure, non_overridable, private :: gather_coarse_dofs_gids_rcv_counts_and_displs     => mlbddc_coarse_gather_coarse_dofs_gids_rcv_counts_and_displs
   procedure, non_overridable, private :: gather_coarse_dofs_gids                           => mlbddc_coarse_gather_coarse_dofs_gids
   procedure, non_overridable, private :: gather_vefs_gids_dofs_objects                     => mlbddc_coarse_gather_vefs_gids_dofs_objects
   procedure, non_overridable, private :: symbolic_setup_dirichlet_problem                  => mlbddc_coarse_symbolic_setup_dirichlet_problem
   procedure, non_overridable, private :: symbolic_setup_dirichlet_solver                   => mlbddc_coarse_symbolic_setup_dirichlet_solver
   procedure, non_overridable, private :: symbolic_setup_constrained_neumann_problem        => mlbddc_coarse_symbolic_setup_constrained_neumann_problem
   procedure, non_overridable, private :: symbolic_setup_constrained_neumann_solver         => mlbddc_coarse_symbolic_setup_constrained_neumann_solver
   procedure, non_overridable, private :: symbolic_setup_coarse_grid_matrix                 => mlbddc_coarse_symbolic_setup_coarse_grid_matrix
   procedure, non_overridable, private :: coarse_grid_matrix_symbolic_assembly              => mlbddc_coarse_coarse_grid_matrix_symbolic_assembly 
   procedure, non_overridable, private :: symbolic_setup_mlbddc_coarse                      => mlbddc_coarse_symbolic_setup_mlbddc_coarse
   procedure, non_overridable, private :: symbolic_setup_coarse_solver                      => mlbddc_coarse_symbolic_setup_coarse_solver

   
   ! Numerical setup-related TBPS
   procedure, non_overridable          :: numerical_setup                                   => mlbddc_coarse_numerical_setup
   procedure, non_overridable, private :: numerical_setup_dirichlet_problem                 => mlbddc_coarse_numerical_setup_dirichlet_problem
   procedure, non_overridable, private :: numerical_setup_dirichlet_solver                  => mlbddc_coarse_numerical_setup_dirichlet_solver
   procedure, non_overridable, private :: numerical_setup_constrained_neumann_problem       => mlbddc_coarse_numerical_constrained_neumann_problem
   procedure, non_overridable, private :: numerical_setup_constrained_neumann_solver        => mlbddc_coarse_numerical_setup_constrained_neumann_solver
   procedure, non_overridable, private :: allocate_coarse_grid_basis                        => mlbddc_coarse_allocate_coarse_grid_basis
   procedure, non_overridable, private :: setup_coarse_grid_basis                           => mlbddc_coarse_setup_coarse_grid_basis
   procedure, non_overridable, private :: compute_subdomain_elmat                           => mlbddc_coarse_compute_subdomain_elmat
   procedure, non_overridable, private :: compute_and_gather_subdomain_elmat                => mlbddc_coarse_compute_and_gather_subdomain_elmat
   procedure, non_overridable, private :: compute_subdomain_elmat_counts_and_displs         => mlbddc_coarse_compute_subdomain_elmat_counts_and_displs
   procedure, non_overridable, private :: numerical_setup_coarse_grid_matrix                => mlbddc_coarse_numerical_setup_coarse_grid_matrix
   procedure, non_overridable, private :: coarse_grid_matrix_numerical_assembly             => mlbddc_coarse_coarse_grid_matrix_numerical_assembly
   procedure, non_overridable, private :: numerical_setup_mlbddc_coarse                     => mlbddc_coarse_numerical_setup_mlbddc_coarse
   procedure, non_overridable, private :: numerical_setup_coarse_solver                     => mlbddc_coarse_numerical_setup_coarse_solver
   
   ! Apply related TBPs
    procedure                           :: apply                                            => mlbddc_coarse_apply
    procedure, non_overridable, private :: apply_par_scalar_array                           => mlbddc_coarse_apply_par_scalar_array
    procedure, non_overridable, private :: solve_coarse_problem                             => mlbddc_coarse_solve_coarse_problem
    procedure, non_overridable, private :: compute_coarse_correction                        => mlbddc_coarse_compute_coarse_correction
    procedure, non_overridable, private :: setup_coarse_grid_residual                       => mlbddc_coarse_setup_coarse_grid_residual
    procedure, non_overridable, private :: compute_coarse_dofs_values                       => mlbddc_coarse_compute_coarse_dofs_values
    procedure, non_overridable, private :: compute_and_gather_coarse_dofs_values            => mlbddc_coarse_compute_and_gather_coarse_dofs_values
    procedure, non_overridable, private :: compute_coarse_dofs_values_counts_and_displs     => mlbddc_coarse_compute_coarse_dofs_values_counts_and_displs
    procedure, non_overridable, private :: coarse_grid_residual_assembly                    => mlbddc_coarse_coarse_grid_residual_assembly
    procedure, non_overridable, private :: scatter_and_interpolate_coarse_grid_correction   => mlbddc_coarse_scatter_and_interpolate_coarse_grid_correction
    procedure, non_overridable, private :: scatter_coarse_grid_correction                   => mlbddc_coarse_scatter_coarse_grid_correction
    procedure, non_overridable, private :: fill_coarse_dofs_values_scattered                => mlbddc_coarse_fill_coarse_dofs_values_scattered
    procedure, non_overridable, private :: interpolate_coarse_grid_correction               => mlbddc_coarse_interpolate_coarse_grid_correction
    procedure, non_overridable, private :: solve_dirichlet_problem                          => mlbddc_coarse_solve_dirichlet_problem
    procedure, non_overridable, private :: apply_A_GI                                       => mlbddc_coarse_apply_A_GI
    procedure, non_overridable, private :: apply_A_IG                                       => mlbddc_coarse_apply_A_IG
    procedure, non_overridable, private :: solve_constrained_neumann_problem                => mlbddc_coarse_solve_constrained_neumann_problem
    procedure, non_overridable, private :: apply_weight_operator                            => mlbddc_coarse_apply_weight_operator
    procedure, non_overridable, private :: create_interior_interface_views                  => mlbddc_coarse_create_interior_interface_views
   
   
   procedure, non_overridable          :: free                                              => mlbddc_coarse_free
   procedure, non_overridable          :: free_clean                                        => mlbddc_coarse_free_clean
   procedure, non_overridable          :: free_symbolic_setup                               => mlbddc_coarse_free_symbolic_setup
   procedure, non_overridable, private :: free_dofs_objects_and_constraint_matrix           => mlbddc_coarse_free_dofs_objects_and_constraint_matrix
   procedure, non_overridable, private :: free_symbolic_setup_dirichlet_problem             => mlbddc_coarse_free_symbolic_setup_dirichlet_problem
   procedure, non_overridable, private :: free_symbolic_setup_dirichlet_solver              => mlbddc_coarse_free_symbolic_setup_dirichlet_solver
   procedure, non_overridable, private :: free_symbolic_setup_constrained_neumann_problem   => mlbddc_coarse_free_symbolic_setup_constrained_neumann_problem
   procedure, non_overridable, private :: free_symbolic_setup_constrained_neumann_solver    => mlbddc_coarse_free_symbolic_setup_constrained_neumann_solver
   procedure, non_overridable, private :: free_symbolic_setup_coarse_solver                 => mlbddc_coarse_free_symbolic_setup_coarse_solver

    procedure, non_overridable          :: free_numerical_setup                             => mlbddc_coarse_free_numerical_setup
    procedure, non_overridable, private :: free_numerical_setup_dirichlet_problem           => mlbddc_coarse_free_numerical_setup_dirichlet_problem
    procedure, non_overridable, private :: free_numerical_setup_dirichlet_solver            => mlbddc_coarse_free_numerical_setup_dirichlet_solver
    procedure, non_overridable, private :: free_numerical_setup_constrained_neumann_problem => mlbddc_coarse_free_numerical_setup_constrained_neumann_problem
    procedure, non_overridable, private :: free_numerical_setup_constrained_neumann_solver  => mlbddc_coarse_free_numerical_setup_constrained_neumann_solver
    procedure, non_overridable, private :: free_coarse_grid_basis                           => mlbddc_coarse_free_coarse_grid_basis
    procedure, non_overridable, private :: free_numerical_setup_coarse_solver               => mlbddc_coarse_free_numerical_setup_coarse_solver
 
   
   procedure, non_overridable           :: get_total_number_coarse_dofs                     => mlbddc_coarse_get_total_number_coarse_dofs
   procedure, non_overridable           :: get_block_number_coarse_dofs                     => mlbddc_coarse_get_block_number_coarse_dofs
   procedure, non_overridable           :: create_field_dofs_object_iterator                => mlbddc_coarse_create_field_dofs_object_iterator
   
   procedure, non_overridable, private  :: get_par_sparse_matrix                            => mlbddc_coarse_get_par_sparse_matrix
   procedure, non_overridable, private  :: get_fe_space                                     => mlbddc_coarse_get_fe_space
   procedure, non_overridable, private  :: get_par_environment                              => mlbddc_coarse_get_par_environment
   procedure, non_overridable, private  :: am_i_l1_task                                     => mlbddc_coarse_am_i_l1_task
   end type mlbddc_coarse_t
 
 public :: mlbddc_t

contains

#include "sbm_dof_object_accessor.i90"
#include "sbm_dof_object_iterator.i90"
#include "sbm_coarse_dof_object_accessor.i90"
#include "sbm_coarse_dof_object_iterator.i90"


  subroutine mlbddc_create ( this, fe_affine_operator )
    implicit none
    class(mlbddc_t)                       , intent(inout) :: this
    type(new_fe_affine_operator_t), target, intent(in)    :: fe_affine_operator
    type(new_par_fe_space_t), pointer :: fe_space
    
    call this%free()
    
    this%fe_affine_operator => fe_affine_operator
    fe_space                => this%get_fe_space()
    call this%create_vector_spaces()
    
    
    assert ( fe_space%get_number_blocks() == 1 )
    assert ( fe_space%get_number_fields() == 1 )
  end subroutine mlbddc_create
  
  subroutine mlbddc_create_vector_spaces (mlbddc)
    implicit none
    class(mlbddc_t), intent(inout)  :: mlbddc
    type(vector_space_t), pointer :: fe_affine_operator_domain_vector_space
    type(vector_space_t), pointer :: fe_affine_operator_range_vector_space
    type(vector_space_t), pointer :: mlbddc_domain_vector_space
    type(vector_space_t), pointer :: mlbddc_range_vector_space
    fe_affine_operator_domain_vector_space => mlbddc%fe_affine_operator%get_domain_vector_space()
    fe_affine_operator_range_vector_space => mlbddc%fe_affine_operator%get_range_vector_space()
    mlbddc_domain_vector_space => mlbddc%get_domain_vector_space()
    mlbddc_range_vector_space => mlbddc%get_range_vector_space()
    call fe_affine_operator_domain_vector_space%clone(mlbddc_domain_vector_space)
    call fe_affine_operator_range_vector_space%clone(mlbddc_range_vector_space)
  end subroutine mlbddc_create_vector_spaces
  
  
  subroutine mlbddc_symbolic_setup ( this )
    implicit none
    class(mlbddc_t), intent(inout) :: this
    type(par_environment_t), pointer :: par_environment  
    par_environment => this%get_par_environment()
    
    if( par_environment%am_i_l1_task() ) then
       if ( par_environment%get_l1_size() > 1 ) then
         call this%setup_dofs_objects_and_constraint_matrix()
       end if
    end if
       
    call this%setup_coarse_fe_space()
    call this%symbolic_setup_coarse_grid_matrix()
    call this%symbolic_setup_mlbddc_coarse()
    
    if ( par_environment%am_i_l1_task() ) then
      if ( par_environment%get_l1_size() > 1 ) then      
        call this%symbolic_setup_dirichlet_problem()
        call this%symbolic_setup_dirichlet_solver()
        call this%symbolic_setup_constrained_neumann_problem()
        call this%symbolic_setup_constrained_neumann_solver()
      else
        call this%symbolic_setup_coarse_solver()
      end if
    end if
  end subroutine mlbddc_symbolic_setup
  
  subroutine mlbddc_setup_dofs_objects_and_constraint_matrix (this)
    implicit none
    class(mlbddc_t)                   , intent(inout) :: this
    type(par_environment_t), pointer :: par_environment
    type(new_par_fe_space_t)   , pointer :: fe_space
    integer(ip)                      :: i 
    
    assert ( this%am_i_l1_task() )
    fe_space => this%get_fe_space()
    call fe_space%setup_dofs_objects_and_constraint_matrix(this%num_dofs_objects_per_field, &
                                                           this%dofs_objects_per_field, &
                                                           this%dofs_objects_gids_per_field, &
                                                           this%vefs_lids_dofs_objects_per_field, &
                                                           this%constraint_matrix)
    
    !write (*,'(a)') '****print mlbddc_setup_dofs_objects_and_constraint_matrix****'
    !write (*,'(a,i10)'  ) 'num_dofs_objects_per_field:', this%num_dofs_objects_per_field 
    !do i=1, fe_space%get_number_fields()
    !  write (*,'(a,i5)') '****FIELD ', i, '****'
    !  write (*,'(a)') 'dofs_objects_per_field'
    !  call this%dofs_objects_per_field(i)%print(6)
    !  call this%constraint_matrix%print(6)
    !end do
    !write (*,'(a)') '****print mlbddc_setup_dofs_objects_and_constraint_matrix****'
  end subroutine mlbddc_setup_dofs_objects_and_constraint_matrix
  
  subroutine mlbddc_setup_coarse_fe_space(this)
  implicit none
  class(mlbddc_t), intent(inout) :: this
  integer(ip)                          :: istat
  integer(ip)                          :: number_fields
  integer(ip), allocatable             :: field_type(:)
  integer(ip), allocatable             :: ptr_dofs_per_fe_and_field(:)
  integer(ip), allocatable             :: coarse_dofs_gids_recv_counts(:)
  integer(ip), allocatable             :: coarse_dofs_gids_displs(:)
  integer(igp), allocatable            :: lst_dofs_gids(:)
  integer(igp), allocatable            :: lst_vefs_gids_dofs_objects(:)
  
  type(par_environment_t)  , pointer     :: par_environment
  type(new_par_fe_space_t)     , pointer     :: fe_space
  type(new_par_triangulation_t), pointer     :: triangulation
  
  par_environment => this%get_par_environment()
  fe_space        => this%get_fe_space()

   ! All MPI tasks (even if they are not involved in the L2 from L1 gather) should also allocate the
   ! allocatable arrays due to the fact that non-allocated allocatable arrays cannot
   ! be passed as actual arguments of dummy arguments that do not have the allocatable attribute 
   ! Otherwise, the code crashes with a segmentation fault.
   call memalloc (0, field_type, __FILE__, __LINE__)
   call memalloc (0, ptr_dofs_per_fe_and_field, __FILE__, __LINE__)
   call memalloc (0, lst_dofs_gids, __FILE__, __LINE__)
   call memalloc (0, lst_vefs_gids_dofs_objects, __FILE__, __LINE__)
   call memalloc (0, coarse_dofs_gids_recv_counts, __FILE__, __LINE__)
   call memalloc (0, coarse_dofs_gids_displs, __FILE__, __LINE__)

  ! L2 tasks gather from L1 tasks all raw data required to set-up the coarse fe space on L2 tasks
  if ( par_environment%am_i_l1_to_l2_task() ) then
     call this%transfer_number_fields(number_fields) 
     call this%transfer_field_type(number_fields, field_type)
     call this%gather_ptr_dofs_per_fe_and_field(number_fields, ptr_dofs_per_fe_and_field)
     call this%gather_coarse_dofs_gids_rcv_counts_and_displs (coarse_dofs_gids_recv_counts, coarse_dofs_gids_displs)
     call this%gather_coarse_dofs_gids(coarse_dofs_gids_recv_counts, coarse_dofs_gids_displs, lst_dofs_gids)
     call this%gather_vefs_gids_dofs_objects(coarse_dofs_gids_recv_counts, coarse_dofs_gids_displs, lst_vefs_gids_dofs_objects)
  end if

  if ( par_environment%am_i_lgt1_task() ) then
     ! lgt1 MPI tasks (recursively) build coarse triangulation
     allocate  ( this%coarse_fe_space, stat = istat )
     check( istat == 0 )
     triangulation => fe_space%get_par_triangulation()
     call this%coarse_fe_space%create (triangulation%get_coarse_triangulation(), &
                                       number_fields, &
                                       field_type, &
                                       ptr_dofs_per_fe_and_field, &
                                       lst_dofs_gids, &
                                       lst_vefs_gids_dofs_objects)
  else
     ! L1 tasks do not hold any piece of the coarse triangulation
     nullify(this%coarse_fe_space)
  end if

  ! All tasks free raw data (see actual reason on the top part of this subroutine)
  call memfree (field_type, __FILE__, __LINE__)
  call memfree (ptr_dofs_per_fe_and_field, __FILE__, __LINE__)
  call memfree (lst_dofs_gids, __FILE__, __LINE__)
  call memfree (lst_vefs_gids_dofs_objects, __FILE__, __LINE__)
  call memfree (coarse_dofs_gids_recv_counts, __FILE__, __LINE__)
  call memfree (coarse_dofs_gids_displs, __FILE__, __LINE__)
end subroutine mlbddc_setup_coarse_fe_space

subroutine mlbddc_transfer_number_fields ( this, number_fields )
  implicit none
  class(mlbddc_t)   , intent(in)             :: this
  integer(ip)       , intent(out)            :: number_fields
  integer(ip)                                :: dummy_integer_ip
  type(par_environment_t), pointer           :: par_environment
  type(new_par_fe_space_t), pointer              :: fe_space

  par_environment => this%get_par_environment()
  assert ( par_environment%am_i_l1_to_l2_task() )
  fe_space        => this%get_fe_space()
  if ( par_environment%am_i_l1_to_l2_root() ) then
     call par_environment%l1_to_l2_transfer(input_data=dummy_integer_ip, &
                                            output_data=number_fields)
  else
     number_fields = fe_space%get_number_fields()
     call par_environment%l1_to_l2_transfer(input_data=number_fields, &
                                            output_data=dummy_integer_ip) 
  end if
end subroutine mlbddc_transfer_number_fields


subroutine mlbddc_transfer_field_type ( this, number_fields, field_type )
  implicit none
  class(mlbddc_t)         , intent(in)    :: this
  integer(ip)             , intent(in)    :: number_fields
  integer(ip), allocatable, intent(inout) :: field_type(:)
  integer(ip)                             :: dummy_integer_array_ip(0)
  type(par_environment_t)  , pointer      :: par_environment
  type(new_par_fe_space_t)     , pointer      :: fe_space
  
  par_environment => this%get_par_environment()
  fe_space        => this%get_fe_space()

  assert ( par_environment%am_i_l1_to_l2_task() )
  if ( par_environment%am_i_l1_to_l2_root() ) then
     if ( allocated (field_type) ) call memfree ( field_type, __FILE__, __LINE__ )
     call memalloc ( number_fields, field_type, __FILE__, __LINE__ )
     call par_environment%l1_to_l2_transfer(input_data=dummy_integer_array_ip, &
                                            output_data=field_type)
  else
     call par_environment%l1_to_l2_transfer(input_data=fe_space%get_field_type(), &
                                            output_data=dummy_integer_array_ip) 
  end if
end subroutine mlbddc_transfer_field_type

subroutine mlbddc_gather_ptr_dofs_per_fe_and_field( this, number_fields, ptr_dofs_per_fe_and_field )
  implicit none
  class(mlbddc_t)   , intent(in)    :: this
  integer(ip)             , intent(in)    :: number_fields
  integer(ip), allocatable, intent(inout) :: ptr_dofs_per_fe_and_field(:)
  integer(ip)                             :: i, num_local_cells
  integer(ip)                             :: dummy_integer_array(0)
  type(par_environment_t)     , pointer   :: par_environment
  type(new_par_fe_space_t)        , pointer   :: fe_space
  type(new_par_triangulation_t)   , pointer   :: triangulation
  type(coarse_triangulation_t), pointer   :: coarse_triangulation
  
  par_environment => this%get_par_environment()  
  fe_space    => this%get_fe_space()
  assert ( par_environment%am_i_l1_to_l2_task() )
  if ( par_environment%am_i_l1_to_l2_root() ) then
     triangulation    => fe_space%get_par_triangulation()
     coarse_triangulation => triangulation%get_coarse_triangulation()
     num_local_cells = coarse_triangulation%get_num_local_cells()
     if (allocated(ptr_dofs_per_fe_and_field)) call memfree ( ptr_dofs_per_fe_and_field, __FILE__, __LINE__ )
     call memalloc (num_local_cells*number_fields+1, ptr_dofs_per_fe_and_field, __FILE__, __LINE__ )
     call par_environment%l2_from_l1_gather( input_data_size = number_fields, &
                                                          input_data      = dummy_integer_array, &
                                                          output_data     = ptr_dofs_per_fe_and_field(2:))
     ptr_dofs_per_fe_and_field(1) = 1
     do i=1, num_local_cells*number_fields
       ptr_dofs_per_fe_and_field(i+1) = ptr_dofs_per_fe_and_field(i) + ptr_dofs_per_fe_and_field(i+1)
     end do
  else
     call par_environment%l2_from_l1_gather( input_data_size = fe_space%get_number_fields(), &
                                                          input_data      = this%num_dofs_objects_per_field, &
                                                          output_data     = dummy_integer_array )
  end if
end subroutine mlbddc_gather_ptr_dofs_per_fe_and_field

  subroutine mlbddc_gather_coarse_dofs_gids_rcv_counts_and_displs( this, recv_counts, displs )
    implicit none
    class(mlbddc_t)     , intent(in)    :: this
    integer(ip) , allocatable , intent(inout) :: recv_counts(:) 
    integer(ip) , allocatable , intent(inout) :: displs(:)
    integer(ip)                               :: i
    integer(ip)                               :: l1_to_l2_size
    integer(ip)                               :: dummy_integer_array(0)
    type(par_environment_t)  , pointer        :: par_environment
    
    par_environment => this%get_par_environment()
    assert ( par_environment%am_i_l1_to_l2_task() )
    if ( par_environment%am_i_l1_to_l2_root() ) then
      l1_to_l2_size = par_environment%get_l1_to_l2_size()
      if ( allocated (recv_counts) ) call memfree ( recv_counts, __FILE__, __LINE__ )
      if ( allocated (displs) ) call memfree ( displs, __FILE__, __LINE__ )
      call memalloc ( l1_to_l2_size, recv_counts, __FILE__, __LINE__ )
      call memalloc ( l1_to_l2_size, displs, __FILE__, __LINE__ )
      call par_environment%l2_from_l1_gather( input_data = 0, &
                                                           output_data = recv_counts ) 
      displs(1) = 0
      do i=2, l1_to_l2_size
        displs(i) = displs(i-1) + recv_counts(i-1)
      end do
    else
      call par_environment%l2_from_l1_gather( input_data  = sum(this%num_dofs_objects_per_field), &
                                                           output_data = dummy_integer_array ) 
    end if
  end subroutine mlbddc_gather_coarse_dofs_gids_rcv_counts_and_displs
  
  subroutine mlbddc_gather_coarse_dofs_gids ( this, recv_counts, displs, lst_gids )
    implicit none
    class(mlbddc_t)     , intent(in)    :: this
    integer(ip)               , intent(in)    :: recv_counts(:)
    integer(ip)               , intent(in)    :: displs(:)
    integer(igp), allocatable , intent(inout) :: lst_gids(:)
    integer(ip)                               :: l1_to_l2_size
    integer(igp)                              :: dummy_integer_array_igp(0)
    integer(ip)                               :: dummy_integer_array_ip(0)
    integer(ip)                               :: i, spos, epos
    integer(igp), allocatable                 :: buffer(:)
    type(par_environment_t)  , pointer        :: par_environment
    type(new_par_fe_space_t)     , pointer        :: fe_space
    
    par_environment => this%get_par_environment()
    assert ( par_environment%am_i_l1_to_l2_task() )
    if ( par_environment%am_i_l1_to_l2_root() ) then
      l1_to_l2_size = par_environment%get_l1_to_l2_size()
      if (allocated(lst_gids)) call memfree ( lst_gids, __FILE__, __LINE__ )
      call memalloc ( displs(l1_to_l2_size), lst_gids, __FILE__, __LINE__ )
      call par_environment%l2_from_l1_gather( input_data_size = 0, &
                                                           input_data      = dummy_integer_array_igp, &
                                                           recv_counts     = recv_counts, &
                                                           displs          = displs, &
                                                           output_data     = lst_gids )
    else
      fe_space    => this%get_fe_space()
      ! Pack dofs_objects_gids_per_field(:) into plain buffer for further data exchange
      call memalloc ( sum(this%num_dofs_objects_per_field), buffer, __FILE__, __LINE__ ) 
      spos = 1
      do i=1, fe_space%get_number_fields()
        epos = spos + this%num_dofs_objects_per_field(i)-1
        buffer(spos:epos) = this%dofs_objects_gids_per_field(i)%a
        spos = epos +1 
      end do
      
      call par_environment%l2_from_l1_gather( input_data_size = size(buffer), &
                                                           input_data      = buffer, &
                                                           recv_counts     = dummy_integer_array_ip, &
                                                           displs          = dummy_integer_array_ip, &
                                                           output_data     = dummy_integer_array_igp )
      call memfree ( buffer, __FILE__, __LINE__ )
    end if    
  end subroutine mlbddc_gather_coarse_dofs_gids

  
  subroutine mlbddc_gather_vefs_gids_dofs_objects ( this, recv_counts, displs, vef_gids )
    implicit none
    class(mlbddc_t)     , intent(in)    :: this
    integer(ip)               , intent(in)    :: recv_counts(:)
    integer(ip)               , intent(in)    :: displs(:)
    integer(igp), allocatable , intent(inout) :: vef_gids(:)
    integer(ip)                               :: l1_to_l2_size
    integer(igp)                              :: dummy_integer_array_igp(0)
    integer(ip)                               :: dummy_integer_array_ip(0)
    integer(ip)                               :: i, j, spos, epos
    integer(igp), allocatable                 :: buffer(:)
    type(par_environment_t)  , pointer        :: par_environment
    type(new_par_fe_space_t)     , pointer    :: fe_space
    type(new_par_triangulation_t), pointer    :: triangulation 
    
    type(dof_object_iterator_t)        :: dofs_object_iterator
    type(dof_object_accessor_t)        :: dof_object
    
    par_environment => this%get_par_environment()
    assert ( par_environment%am_i_l1_to_l2_task() )
    if ( par_environment%am_i_l1_to_l2_root() ) then
      l1_to_l2_size = par_environment%get_l1_to_l2_size()
      if (allocated(vef_gids)) call memfree ( vef_gids, __FILE__, __LINE__ )
      call memalloc ( displs(l1_to_l2_size), vef_gids, __FILE__, __LINE__ )
      call par_environment%l2_from_l1_gather( input_data_size = 0, &
                                              input_data      = dummy_integer_array_igp, &
                                              recv_counts     = recv_counts, &
                                              displs          = displs, &
                                              output_data     = vef_gids )
    else
      fe_space      => this%get_fe_space()
      triangulation => fe_space%get_par_triangulation()
      ! Pack vefs_lids_dofs_objects_per_field(:) into plain buffer for further data exchange
      call memalloc ( sum(this%num_dofs_objects_per_field), buffer, __FILE__, __LINE__ ) 
      spos = 1
      do i=1, fe_space%get_number_fields()
        epos = spos + this%num_dofs_objects_per_field(i)-1
        dofs_object_iterator = this%create_field_dofs_object_iterator(i)
        do j=spos, epos
          call dofs_object_iterator%current(dof_object)
          buffer(j) = dof_object%get_gid()
          call dofs_object_iterator%next()
        end do  
        spos = epos +1 
      end do
      call par_environment%l2_from_l1_gather( input_data_size = size(buffer), &
                                              input_data      = buffer, &
                                              recv_counts     = dummy_integer_array_ip, &
                                              displs          = dummy_integer_array_ip, &
                                              output_data     = dummy_integer_array_igp )
      call memfree ( buffer, __FILE__, __LINE__ )
    end if    
  end subroutine mlbddc_gather_vefs_gids_dofs_objects
  
  subroutine mlbddc_symbolic_setup_dirichlet_problem ( this) 
    implicit none
    class(mlbddc_t)           , intent(inout) :: this
    type(new_par_fe_space_t)      , pointer       :: fe_space
    type(par_sparse_matrix_t) , pointer       :: par_sparse_matrix
    type(sparse_matrix_t)     , pointer       :: A
    
    assert ( this%am_i_l1_task() )
    fe_space => this%get_fe_space()
    par_sparse_matrix => this%get_par_sparse_matrix()
    A => par_sparse_matrix%get_sparse_matrix()
    if ( A%get_symmetric_storage() ) then
       call A%split_2x2_symbolic( num_row = fe_space%get_block_number_interior_dofs(1), &
                                  num_col = fe_space%get_block_number_interior_dofs(1), &
                                  A_II    = this%A_II, &
                                  A_IG    = this%A_IG, &
                                  A_GG    = this%A_GG)
    else
       call A%split_2x2_symbolic( num_row = fe_space%get_block_number_interior_dofs(1), &
                                  num_col = fe_space%get_block_number_interior_dofs(1), &
                                  A_II    = this%A_II, &
                                  A_IG    = this%A_IG, &
                                  A_GI    = this%A_GI, &
                                  A_GG    = this%A_GG)
    end if
  end subroutine mlbddc_symbolic_setup_dirichlet_problem
  
  subroutine mlbddc_symbolic_setup_dirichlet_solver ( this) 
    implicit none
    class(mlbddc_t), intent(inout) :: this
    assert ( this%am_i_l1_task() )
    call this%dirichlet_solver%set_type(name   = pardiso_mkl)
    call this%dirichlet_solver%set_matrix(matrix = this%A_II)
    call this%dirichlet_solver%symbolic_setup() 
  end subroutine mlbddc_symbolic_setup_dirichlet_solver
  
  subroutine mlbddc_symbolic_setup_constrained_neumann_problem(this)
    implicit none
    class(mlbddc_t)           , intent(inout) :: this
    type(par_environment_t)   , pointer       :: par_environment
    type(par_sparse_matrix_t) , pointer       :: par_sparse_matrix
    type(sparse_matrix_t)     , pointer       :: A
    
    assert ( this%am_i_l1_task() )
    par_sparse_matrix => this%get_par_sparse_matrix()
    A => par_sparse_matrix%get_sparse_matrix()
    
    call A%expand_matrix_symbolic(C_T               = this%constraint_matrix, &
                                  to                = this%constrained_neumann_matrix, &
                                  symmetric_storage = A%get_symmetric_storage(), & 
                                  symmetric         = A%is_symmetric(), &
                                  sign              = SPARSE_MATRIX_SIGN_INDEFINITE)
    
  end subroutine mlbddc_symbolic_setup_constrained_neumann_problem
  
  subroutine mlbddc_symbolic_setup_constrained_neumann_solver(this)
    implicit none
    class(mlbddc_t)           , intent(inout) :: this
    assert ( this%am_i_l1_task() )
    call this%constrained_neumann_solver%set_type(name   = pardiso_mkl)
    call this%constrained_neumann_solver%set_matrix(matrix = this%constrained_neumann_matrix)
    call this%constrained_neumann_solver%symbolic_setup() 
  end subroutine mlbddc_symbolic_setup_constrained_neumann_solver
  
  subroutine mlbddc_symbolic_setup_coarse_grid_matrix ( this )
    implicit none
    class(mlbddc_t)           , intent(inout) :: this
    type(par_environment_t)   , pointer       :: L1_environment
    type(par_environment_t)   , pointer       :: L2_environment
    integer(ip)                               :: istat
    
    L1_environment => this%get_par_environment()
    if ( L1_environment%am_i_lgt1_task() ) then
       ! lgt1 MPI tasks symbolically set-up this%coarse_grid_matrix
       allocate  ( this%coarse_grid_matrix, stat = istat )
       check( istat == 0 )

       L2_environment => L1_environment%get_next_level()
       call this%coarse_grid_matrix%create( p_env             = L2_environment, &
                                            dof_import        = this%coarse_fe_space%get_block_dof_import(1), &
                                            symmetric_storage = .true., &
                                            is_symmetric      = .true., &
                                            sign              = SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE )
       
       if ( L2_environment%am_i_l1_task() ) then
          call this%coarse_grid_matrix_symbolic_assembly()
       end if
    else       
       ! L1 tasks do not hold any piece of the coarse triangulation
       nullify(this%coarse_grid_matrix)
    end if
  end subroutine mlbddc_symbolic_setup_coarse_grid_matrix
  
  subroutine mlbddc_coarse_grid_matrix_symbolic_assembly ( this )
    implicit none
    class(mlbddc_t)           , intent(inout) :: this
    type(par_environment_t)   , pointer       :: L1_environment
    type(par_environment_t)   , pointer       :: L2_environment
    type(coarse_fe_iterator_t) :: iterator
    type(coarse_fe_accessor_t) :: coarse_fe
    type(i1p_t), allocatable :: elem2dof(:)
    logical, pointer :: field_coupling(:,:)
    integer(ip) :: ifield, jfield, i, j, istat
    
    L1_environment => this%get_par_environment()
    L2_environment => L1_environment%get_next_level()
    assert ( associated (L2_environment) )
    assert ( L2_environment%am_i_l1_task() )
    
    allocate ( elem2dof(this%coarse_fe_space%get_number_fields()), stat=istat)
    check(istat==0)     
    
    field_coupling => this%coarse_fe_space%get_field_coupling()
    
    iterator  = this%coarse_fe_space%create_coarse_fe_iterator()
    do while ( .not. iterator%has_finished() )
       coarse_fe = iterator%current()
       call coarse_fe%get_elem2dof(elem2dof)
       do ifield=1, this%coarse_fe_space%get_number_fields()
          do jfield=1, this%coarse_fe_space%get_number_fields()
             if ((field_coupling(ifield,jfield))) then
               do j=1, size(elem2dof(jfield)%p)
                 do i=1, size(elem2dof(ifield)%p)
                    call this%coarse_grid_matrix%insert(ia=elem2dof(ifield)%p(i), &
                                                        ja=elem2dof(jfield)%p(j) )
                 end do
               end do  
             end if
          end do
       end do
       call iterator%next()
    end do
    
    deallocate ( elem2dof, stat=istat )
    check(istat==0)     
    
    call this%coarse_grid_matrix%convert(csr_format)
    
    !call this%coarse_grid_matrix%print(6)
  end subroutine mlbddc_coarse_grid_matrix_symbolic_assembly
  
  
  subroutine mlbddc_symbolic_setup_mlbddc_coarse(this)
    implicit none
    class(mlbddc_t), intent(inout)     :: this
    type(par_environment_t)  , pointer :: par_environment
    integer(ip)                        :: istat
    
    par_environment => this%get_par_environment()
    if ( par_environment%am_i_lgt1_task() ) then
     ! lgt1 MPI tasks symbolically setup mlbddc coarse
     allocate  ( this%mlbddc_coarse, stat = istat )
     check( istat == 0 )
     call this%mlbddc_coarse%create(this%coarse_fe_space, &
                                    this%coarse_grid_matrix)
     call this%mlbddc_coarse%symbolic_setup()
    else
      ! L1 tasks do not hold any piece of the coarse triangulation
      nullify(this%mlbddc_coarse)
     end if
  end subroutine mlbddc_symbolic_setup_mlbddc_coarse
  
  subroutine mlbddc_symbolic_setup_coarse_solver(this)
    implicit none
    class(mlbddc_t), intent(inout)     :: this
    type(par_environment_t)  , pointer :: par_environment
    type(par_sparse_matrix_t), pointer :: par_sparse_matrix
    par_environment => this%get_par_environment()
    assert ( par_environment%get_l1_size() == 1 )
    assert ( par_environment%am_i_l1_task() )
    par_sparse_matrix => this%get_par_sparse_matrix()
    call this%coarse_solver%set_type(name   = pardiso_mkl)
    call this%coarse_solver%set_matrix(matrix = par_sparse_matrix%get_sparse_matrix())
    call this%coarse_solver%symbolic_setup() 
  end subroutine mlbddc_symbolic_setup_coarse_solver
   
  subroutine mlbddc_numerical_setup ( this )
    implicit none
    class(mlbddc_t)        , intent(inout)   :: this
    type(par_environment_t), pointer         :: par_environment
    
    par_environment => this%get_par_environment()
    if ( par_environment%am_i_l1_task() ) then
      if ( par_environment%get_l1_size() > 1 ) then  
        call this%numerical_setup_constrained_neumann_problem()
        call this%numerical_setup_constrained_neumann_solver()
        call this%setup_coarse_grid_basis()
      end if
    end if  
    
    call this%numerical_setup_coarse_grid_matrix()
    call this%numerical_setup_mlbddc_coarse()
    
    if ( par_environment%am_i_l1_task() ) then
      if ( par_environment%get_l1_size() > 1 ) then  
        call this%numerical_setup_dirichlet_problem()
        call this%numerical_setup_dirichlet_solver()
      else if ( par_environment%get_l1_size() == 1 ) then
        call this%numerical_setup_coarse_solver()
      end if
    end if  
  end subroutine mlbddc_numerical_setup
  
  
  subroutine mlbddc_numerical_setup_dirichlet_problem (this) 
    implicit none
    class(mlbddc_t)           , intent(inout) :: this
    type(new_par_fe_space_t)      , pointer       :: fe_space
    type(par_sparse_matrix_t) , pointer       :: par_sparse_matrix
    type(sparse_matrix_t)     , pointer       :: A
    
    assert ( this%am_i_l1_task() )
    fe_space => this%get_fe_space()
    par_sparse_matrix => this%get_par_sparse_matrix()
    A => par_sparse_matrix%get_sparse_matrix()
    if ( A%get_symmetric_storage() ) then
       call A%split_2x2_numeric(  num_row = fe_space%get_block_number_interior_dofs(1), &
                                  num_col = fe_space%get_block_number_interior_dofs(1), &
                                  A_II    = this%A_II, &
                                  A_IG    = this%A_IG, &
                                  A_GG    = this%A_GG)
    else
       call A%split_2x2_numeric(  num_row = fe_space%get_block_number_interior_dofs(1), &
                                  num_col = fe_space%get_block_number_interior_dofs(1), &
                                  A_II    = this%A_II, &
                                  A_IG    = this%A_IG, &
                                  A_GI    = this%A_GI, &
                                  A_GG    = this%A_GG)
    end if
    
    !call A%print(6)
    !call this%A_II%print(6)
    !call this%A_IG%print(6)
    !call this%A_GG%print(6)
  end subroutine mlbddc_numerical_setup_dirichlet_problem
  
  subroutine mlbddc_numerical_setup_dirichlet_solver(this)
    implicit none
    class(mlbddc_t)           , intent(inout) :: this
    assert ( this%am_i_l1_task() )
    call this%dirichlet_solver%numerical_setup() 
    !call this%dirichlet_solver%log_info()
  end subroutine mlbddc_numerical_setup_dirichlet_solver
  
  subroutine mlbddc_numerical_constrained_neumann_problem(this)
   implicit none
   class(mlbddc_t)           , intent(inout) :: this
   type(new_par_fe_space_t)      , pointer       :: fe_space
   type(par_sparse_matrix_t) , pointer       :: par_sparse_matrix
   type(sparse_matrix_t)     , pointer       :: A
    
    assert ( this%am_i_l1_task() )
    fe_space => this%get_fe_space()
    par_sparse_matrix => this%get_par_sparse_matrix()
    A => par_sparse_matrix%get_sparse_matrix()
    call A%expand_matrix_numeric(C_T = this%constraint_matrix, &
                                 to  = this%constrained_neumann_matrix)
    
    !call this%constrained_neumann_matrix%print(6)
  end subroutine  mlbddc_numerical_constrained_neumann_problem
  
  subroutine mlbddc_numerical_setup_constrained_neumann_solver(this)
    implicit none
    class(mlbddc_t)           , intent(inout) :: this
    assert ( this%am_i_l1_task() )
    call this%constrained_neumann_solver%numerical_setup() 
    !call this%constrained_neumann_solver%log_info()
  end subroutine mlbddc_numerical_setup_constrained_neumann_solver
  
  subroutine mlbddc_allocate_coarse_grid_basis ( this )
    implicit none
    class(mlbddc_t)           , intent(inout) :: this
    type(new_par_fe_space_t)      , pointer       :: fe_space
    assert ( this%am_i_l1_task() )
    fe_space => this%get_fe_space() 
    call this%free_coarse_grid_basis()
    call memalloc (fe_space%get_block_number_dofs(1), &
                   this%get_block_number_coarse_dofs(1), &
                   this%Phi, &
                   __FILE__, &
                   __LINE__) 
  end subroutine mlbddc_allocate_coarse_grid_basis
  
  subroutine mlbddc_setup_coarse_grid_basis ( this )
     implicit none
     class(mlbddc_t)           , intent(inout) :: this
     type(new_par_fe_space_t)      , pointer       :: fe_space
     real(rp), allocatable                     :: work1(:,:)
     real(rp), allocatable                     :: work2(:,:)
     integer(ip)                               :: i, j
     
     assert ( this%am_i_l1_task() )
     fe_space => this%get_fe_space() 
     
     call memalloc ( this%constrained_neumann_matrix%get_num_rows(), &
                     this%get_block_number_coarse_dofs(1), &
                     work1, __FILE__,__LINE__ )
     
     call memalloc ( this%constrained_neumann_matrix%get_num_rows(), &
                     this%get_block_number_coarse_dofs(1), &
                     work2, __FILE__,__LINE__ )
     
     work1 = 0.0_rp
     work2 = 0.0_rp
     j=1
     do i = fe_space%get_block_number_dofs(1)+1, this%constrained_neumann_matrix%get_num_rows()
       work1 (i,j) = 1.0_rp
       j = j + 1 
     end do
     
     call this%constrained_neumann_solver%solve(work1, &
                                                work2)
     
     call this%allocate_coarse_grid_basis()
     this%Phi = work2 (1:fe_space%get_block_number_dofs(1),:) 
     
     !do, i=1,this%constrained_neumann_matrix%get_num_rows()
     !  write(*,"(100g15.5)") ( work1(i,j), j=1,this%get_block_number_coarse_dofs(1) )
     !enddo
     !do, i=1,this%constrained_neumann_matrix%get_num_rows()
     !  write(*,"(100g15.5)") ( work2(i,j), j=1,this%get_block_number_coarse_dofs(1) )
     !enddo
     call memfree ( work1, __FILE__,__LINE__ )
     call memfree ( work2, __FILE__,__LINE__ )
  end subroutine mlbddc_setup_coarse_grid_basis 
  
  ! Computes subdomain_elmat = \Phi^t A_i \Phi 
  subroutine mlbddc_compute_subdomain_elmat ( this, subdomain_elmat )
    implicit none
    class(mlbddc_t)      , intent(in)    :: this
    real(rp), allocatable, intent(inout) :: subdomain_elmat(:,:)
    
    real(rp), allocatable :: work(:,:)
    type(par_sparse_matrix_t), pointer :: par_sparse_matrix
    type(sparse_matrix_t), pointer :: A
    
    assert ( this%am_i_l1_task() )
    par_sparse_matrix => this%get_par_sparse_matrix()
    A => par_sparse_matrix%get_sparse_matrix()
    
    if ( allocated(subdomain_elmat) ) then
      call memfree ( subdomain_elmat, __FILE__, __LINE__ )
    end if
    
    call memalloc ( this%get_block_number_coarse_dofs(1), &
                    this%get_block_number_coarse_dofs(1), &
                    subdomain_elmat, &
                    __FILE__, &
                    __LINE__ );
    
    call memalloc ( A%get_num_rows(), & 
                    this%get_block_number_coarse_dofs(1), &
                    work, &
                    __FILE__, & 
                    __LINE__ )

    work = 0.0_rp
    call A%apply_to_dense_matrix ( n     = this%get_block_number_coarse_dofs(1), &
                                   alpha = 1.0_rp, &
                                   ldb   = A%get_num_rows(), &
                                   b     = this%Phi, &
                                   beta  = 0.0_rp, &
                                   ldc   = A%get_num_rows(), &
                                   c     = work )        
    subdomain_elmat = 0.0_rp
#ifdef ENABLE_BLAS
    call DGEMM( 'T', &
                'N', &
                this%get_block_number_coarse_dofs(1), &
                this%get_block_number_coarse_dofs(1), &
                A%get_num_rows(), &
                1.0, &
                this%Phi, &
                A%get_num_rows() , &
                work, &
                A%get_num_rows(), &
                0.0, &
                subdomain_elmat, &
                this%get_block_number_coarse_dofs(1))
#else
    write (0,*) 'Error: mlbddc.f90 was not compiled with -DENABLE_BLAS.'
    write (0,*) 'Error: You must activate this cpp macro in order to use the BLAS'
    check(.false.)    
#endif
    call memfree ( work, __FILE__, __LINE__)
  end subroutine mlbddc_compute_subdomain_elmat
 
  subroutine mlbddc_compute_subdomain_elmat_counts_and_displs ( this, counts, displs )
    implicit none
    class(mlbddc_t)           , intent(inout) :: this
    integer(ip), allocatable  , intent(inout) :: counts(:)
    integer(ip), allocatable  , intent(inout) :: displs(:)
    type(par_environment_t), pointer :: par_environment
    integer(ip) :: i, l1_to_l2_size
    type(coarse_fe_iterator_t) :: iterator
    type(coarse_fe_accessor_t) :: coarse_fe
    
    par_environment => this%get_par_environment()
    assert (par_environment%am_i_l1_to_l2_root())
    l1_to_l2_size = par_environment%get_l1_to_l2_size()
    if ( allocated(counts) ) call memfree ( counts, __FILE__, __LINE__ )
    if ( allocated(displs) ) call memfree ( displs, __FILE__, __LINE__ )
    call memalloc ( l1_to_l2_size, counts, __FILE__, __LINE__ )
    call memalloc ( l1_to_l2_size, displs, __FILE__, __LINE__ )
    
    i = 1
    counts(l1_to_l2_size) = 0
    iterator = this%coarse_fe_space%create_coarse_fe_iterator()
    do while ( .not. iterator%has_finished() )
      coarse_fe = iterator%current()
      if ( coarse_fe%is_local() ) then
        counts(i) = coarse_fe%get_number_dofs()**2    
        i = i +1
      end if
      call iterator%next()
    end do
    
    displs(1) = 0
    do i=2, l1_to_l2_size
      displs(i) = displs(i-1) + counts(i-1)
    end do
  end subroutine mlbddc_compute_subdomain_elmat_counts_and_displs
    
  
  subroutine mlbddc_compute_and_gather_subdomain_elmat ( this, subdomain_elmat_gathered )
     implicit none
     class(mlbddc_t)      , intent(inout) :: this
     real(rp), allocatable, intent(inout) :: subdomain_elmat_gathered(:)
     type(par_environment_t), pointer :: par_environment
     integer(ip), allocatable :: counts(:)
     integer(ip), allocatable :: displs(:)
     real(rp), allocatable :: subdomain_elmat(:,:)
     real(rp) :: dummy_real_array_rp_2D(0,0)
     real(rp) :: dummy_real_array_rp_1D(0)
     integer(ip) :: dummy_integer_array_ip(0)
     integer(ip) :: l1_to_l2_size
     par_environment => this%get_par_environment()
     assert (par_environment%am_i_l1_to_l2_task())
     
     if ( par_environment%am_i_l1_to_l2_root() ) then
       call this%compute_subdomain_elmat_counts_and_displs(counts, displs)
       l1_to_l2_size = par_environment%get_l1_to_l2_size()
       
       if ( allocated(subdomain_elmat_gathered) ) & 
         call memfree (subdomain_elmat_gathered, __FILE__, __LINE__)
         
       call memalloc ( displs(l1_to_l2_size), & 
                       subdomain_elmat_gathered, & 
                       __FILE__, __LINE__ ) 
       
       call par_environment%l2_from_l1_gather( input_data      = dummy_real_array_rp_2D, &
                                               recv_counts     = counts, &
                                               displs          = displs, &
                                               output_data     = subdomain_elmat_gathered )
       
       call memfree ( counts, __FILE__, __LINE__ )
       call memfree ( displs, __FILE__, __LINE__ )
     else
       call this%compute_subdomain_elmat(subdomain_elmat)
       call par_environment%l2_from_l1_gather( input_data      = subdomain_elmat, &
                                               recv_counts     = dummy_integer_array_ip, &
                                               displs          = dummy_integer_array_ip, &
                                               output_data     = dummy_real_array_rp_1D )
       call memfree ( subdomain_elmat, __FILE__, __LINE__ )
     end if     
  end subroutine mlbddc_compute_and_gather_subdomain_elmat
  
  subroutine mlbddc_numerical_setup_coarse_grid_matrix ( this )
    implicit none
    class(mlbddc_t)           , intent(inout) :: this
    type(par_environment_t)   , pointer       :: L1_environment
    type(par_environment_t)   , pointer       :: L2_environment
    integer(ip)                               :: istat
    real(rp), allocatable                     :: subdomain_elmat_gathered(:)
    L1_environment => this%get_par_environment()
    if ( L1_environment%am_i_l1_to_l2_task() ) then
      call this%compute_and_gather_subdomain_elmat(subdomain_elmat_gathered)
    end if
    if ( L1_environment%am_i_lgt1_task() ) then
      L2_environment => L1_environment%get_next_level() 
      if ( L2_environment%am_i_l1_task() ) then
          call this%coarse_grid_matrix_numerical_assembly(subdomain_elmat_gathered)
       end if
    end if
    ! subdomain_elmat_gathered is only allocated on L2 MPI tasks
    if (allocated(subdomain_elmat_gathered) ) call memfree(subdomain_elmat_gathered, __FILE__, __LINE__ )
  end subroutine mlbddc_numerical_setup_coarse_grid_matrix
  
  subroutine mlbddc_coarse_grid_matrix_numerical_assembly ( this, subdomain_elmat_gathered )
    implicit none
    class(mlbddc_t)           , intent(inout) :: this
    real(rp)                  , intent(in)    :: subdomain_elmat_gathered(*)
    
    type(par_environment_t)   , pointer       :: L1_environment
    type(par_environment_t)   , pointer       :: L2_environment
    type(coarse_fe_iterator_t) :: iterator
    type(coarse_fe_accessor_t) :: coarse_fe
    type(i1p_t), allocatable :: elem2dof(:)
    logical, pointer :: field_coupling(:,:)
    integer(ip) :: ifield, jfield, i, j, istat, current
    
    L1_environment => this%get_par_environment()
    L2_environment => L1_environment%get_next_level()
    assert ( associated (L2_environment) )
    assert ( L2_environment%am_i_l1_task() )
    allocate ( elem2dof(this%coarse_fe_space%get_number_fields()), stat=istat)
    check(istat==0)     
    field_coupling => this%coarse_fe_space%get_field_coupling()
    current = 1
    iterator  = this%coarse_fe_space%create_coarse_fe_iterator()
    do while ( .not. iterator%has_finished() )
       coarse_fe = iterator%current()
       call coarse_fe%get_elem2dof(elem2dof)
       if ( coarse_fe%is_local() ) then
          do ifield=1, this%coarse_fe_space%get_number_fields()
             do jfield=1, this%coarse_fe_space%get_number_fields()
                if ((field_coupling(ifield,jfield))) then
                   do j=1, size(elem2dof(jfield)%p)
                      do i=1, size(elem2dof(ifield)%p)
                         call this%coarse_grid_matrix%insert(ia =  elem2dof(ifield)%p(i), &
                                                             ja  = elem2dof(jfield)%p(j), &
                                                             val = subdomain_elmat_gathered(current) )
                         current = current + 1
                      end do
                   end do
                end if
             end do
          end do
       end if
       call iterator%next()
    end do
    deallocate ( elem2dof, stat=istat )
    check(istat==0)     
    call this%coarse_grid_matrix%convert(csr_format)
    !call this%coarse_grid_matrix%print(6)
  end subroutine mlbddc_coarse_grid_matrix_numerical_assembly
  
  subroutine mlbddc_numerical_setup_mlbddc_coarse ( this )
   implicit none
   class(mlbddc_t)           , intent(inout) :: this
   type(par_environment_t)  , pointer :: par_environment
    
    par_environment => this%get_par_environment()
    if ( par_environment%am_i_lgt1_task() ) then
     call this%mlbddc_coarse%numerical_setup()
    end if
  end subroutine mlbddc_numerical_setup_mlbddc_coarse 
  
  subroutine mlbddc_numerical_setup_coarse_solver(this)
    implicit none
    class(mlbddc_t)           , intent(inout) :: this
    type(par_environment_t)  , pointer :: par_environment
    par_environment => this%get_par_environment()
    assert ( par_environment%get_l1_size() == 1 )
    assert ( par_environment%am_i_l1_task() )
    call this%coarse_solver%numerical_setup() 
    !call this%coarse_solver%log_info()
  end subroutine mlbddc_numerical_setup_coarse_solver
  
  !=============================================================================
  subroutine mlbddc_apply (op, x, y)
    implicit none
    ! Parameters
    class(mlbddc_t)   , intent(in)    :: op
    class(vector_t)   , intent(in)    :: x
    class(vector_t)   , intent(inout) :: y
    call op%abort_if_not_in_domain(x)
    call op%abort_if_not_in_range(y)
    call x%GuardTemp()
    select type(x)
       class is (par_scalar_array_t)
       select type(y)
          class is(par_scalar_array_t)
          call op%apply_par_scalar_array( x, y )
       end select
    end select
    call x%CleanTemp()
  end subroutine mlbddc_apply
  
  subroutine mlbddc_apply_par_scalar_array ( this, x, y )
    implicit none
    class(mlbddc_t)         , intent(in)    :: this
    type(par_scalar_array_t), intent(in)    :: x
    type(par_scalar_array_t), intent(inout) :: y
       
    type(par_scalar_array_t) :: y_I, y_G
    type(par_scalar_array_t) :: residual, residual_I, residual_G
    type(par_scalar_array_t) :: delta, delta_I, delta_G
    type(par_scalar_array_t) :: constrained_neumann_correction, &
                                constrained_neumann_correction_I, &
                                constrained_neumann_correction_G
    type(par_scalar_array_t) :: coarse_correction, &
                                coarse_correction_I, &
                                coarse_correction_G
    type(par_environment_t), pointer :: par_environment

    par_environment => this%get_par_environment()
    if ( par_environment%get_l1_size() == 1 ) then
      call this%solve_coarse_problem(x, y)
    else ! l1_size > 1 .or. l1_size < 1    
      call this%create_interior_interface_views(y, &
                                                y_I, & 
                                                y_G)
 
      ! Clone, copy, and create views for input residual
      call residual%clone(x)
      call residual%copy(x)
      call this%create_interior_interface_views(residual, &
                                                residual_I, & 
                                                residual_G)
    
      ! Clone, init, and create views for delta
      call delta%clone(x)
      call this%create_interior_interface_views(delta, &
                                                delta_I, & 
                                                delta_G)
    
      ! Clone and create views for constrained_neumann_correction
      call constrained_neumann_correction%clone(x)
      call this%create_interior_interface_views(constrained_neumann_correction, &
                                                constrained_neumann_correction_I, & 
                                                constrained_neumann_correction_G)
    
      ! Clone and create views for coarse_correction
      call coarse_correction%clone(x)
      call this%create_interior_interface_views(coarse_correction, &
                                                coarse_correction_I, & 
                                                coarse_correction_G)
    
    
      ! 1) Compute delta_I = A_II^-1 residual_I,   delta_G = 0
      call this%solve_dirichlet_problem( residual_I, delta_I )
    
      ! 2) y_I = delta_I
      call y_I%copy(delta_I)
    
      ! 3) ! I-A P_D = [ 0                    0 ]
           !           [ -A_GI A_II^{-1}      I ]
      !delta_G = A_GI * A_II^{-1}*residual_I
      call this%apply_A_GI(delta_I, delta_G)
      call delta%comm()
      call residual_I%init(0.0_rp)
      ! residual_G = residual_G - delta_G =
      !              residual_G - sum(A_GI*A_II^{-1}*residual_I,i=1..P)
      call residual_G%axpby(-1.0_rp, delta_G, 1.0_rp)

      ! 4) Compute delta_G = A_{BDDC}^{-1} r
      call this%apply_weight_operator(residual,residual) 
      call this%compute_coarse_correction(residual, coarse_correction)
      call this%solve_constrained_neumann_problem ( residual, constrained_neumann_correction )
      call delta_G%copy (coarse_correction_G)
      call delta_G%axpby(1.0_rp, constrained_neumann_correction_G,1.0_rp)
      call this%apply_weight_operator(delta,delta) 
      call delta%comm()
      ! y_G = delta_G
      call y_G%copy(delta_G)

      ! 5) I-P_D A = [ 0      -A_II^{-1} A_IG ]
      !              [ 0                   I  ]
      call this%apply_A_IG(delta_G, residual_I)
      call this%solve_dirichlet_problem(residual_I, delta_I)
      ! y_I = y_I - delta_I
      call y_I%axpby(-1.0_rp, delta_I, 1.0_rp)

      call constrained_neumann_correction_I%free()
      call constrained_neumann_correction_G%free()
      call constrained_neumann_correction%free()
      call coarse_correction_I%free()
      call coarse_correction_G%free()
      call coarse_correction%free()
      call delta_I%free()
      call delta_G%free()
      call delta%free()
      call residual_I%free()
      call residual_G%free()
      call residual%free()
      call y_I%free()
      call y_G%free()
    end if
  end subroutine mlbddc_apply_par_scalar_array
  
  subroutine mlbddc_solve_coarse_problem(this,x,y)
    implicit none
    class(mlbddc_t)            , intent(in)    :: this
    type(par_scalar_array_t)   , intent(in)    :: x
    type(par_scalar_array_t)   , intent(inout) :: y
    type(par_environment_t)    , pointer       :: par_environment
    type(serial_scalar_array_t), pointer       :: serial_y 

    par_environment => this%get_par_environment()
    assert ( par_environment%get_l1_size() == 1 )
    assert ( par_environment%am_i_l1_task() )

    serial_y => y%get_serial_scalar_array()
    call this%coarse_solver%solve(x%get_serial_scalar_array(), &
                                  serial_y)
  end subroutine mlbddc_solve_coarse_problem
  
  subroutine mlbddc_compute_coarse_correction (this, residual, coarse_correction)
    implicit none
    class(mlbddc_t)         , intent(in)    :: this
    type(par_scalar_array_t), intent(in)    :: residual
    type(par_scalar_array_t), intent(inout) :: coarse_correction
    
    type(par_scalar_array_t) :: coarse_residual
    type(par_scalar_array_t) :: coarse_coarse_correction
    type(par_environment_t), pointer :: L1_environment
    type(par_environment_t), pointer :: L2_environment
    
    L1_environment => this%get_par_environment()
    if ( L1_environment%am_i_lgt1_task() ) then
       L2_environment => L1_environment%get_next_level()
       call coarse_residual%create_and_allocate( p_env      = L2_environment, &
                                                 dof_import = this%coarse_fe_space%get_block_dof_import(1) )
       call coarse_coarse_correction%clone(coarse_residual)
    end if   
    
    ! Transfer residual from L1 to L2 tasks
    call this%setup_coarse_grid_residual ( residual, coarse_residual )
    
    ! Solve coarse problem on > L1 tasks
    if ( L1_environment%am_i_lgt1_task() ) then
       call this%mlbddc_coarse%apply(coarse_residual, coarse_coarse_correction)
    end if
    call this%scatter_and_interpolate_coarse_grid_correction ( coarse_coarse_correction, coarse_correction )
    if ( L1_environment%am_i_lgt1_task() ) then
       call coarse_coarse_correction%free()
       call coarse_residual%free()
    end if
    
  end subroutine mlbddc_compute_coarse_correction
  
  subroutine mlbddc_setup_coarse_grid_residual ( this, vector, coarse_grid_vector )
    implicit none
    class(mlbddc_t)         , intent(in)    :: this
    type(par_scalar_array_t), intent(in)    :: vector
    type(par_scalar_array_t), intent(inout) :: coarse_grid_vector
    type(par_environment_t), pointer        :: L1_environment
    type(par_environment_t), pointer        :: L2_environment
    real(rp), allocatable                   :: coarse_dofs_values_gathered(:)
    
    L1_environment => this%get_par_environment()
    if ( L1_environment%am_i_l1_to_l2_task() ) then
      call this%compute_and_gather_coarse_dofs_values(vector, coarse_dofs_values_gathered)
    end if
    
    if ( L1_environment%am_i_lgt1_task() ) then
       call coarse_grid_vector%init(0.0_rp) 
       L2_environment => L1_environment%get_next_level()
       if ( L2_environment%am_i_l1_task() ) then
          call this%coarse_grid_residual_assembly(coarse_dofs_values_gathered, coarse_grid_vector)
       end if
    else
       call coarse_grid_vector%free()
    end if
    if ( allocated(coarse_dofs_values_gathered) ) call memfree(coarse_dofs_values_gathered, __FILE__, __LINE__ )
  end subroutine mlbddc_setup_coarse_grid_residual
  
  
  ! Computes subdomain_elvec = \Phi^t v_i
  subroutine mlbddc_compute_coarse_dofs_values ( this, vector, coarse_dofs_values )
    implicit none
    class(mlbddc_t)         , intent(in) :: this
    type(par_scalar_array_t), intent(in) :: vector
    real(rp), allocatable, intent(inout) :: coarse_dofs_values(:)
    type(serial_scalar_array_t), pointer :: serial_scalar_array
    
    assert ( this%am_i_l1_task() )
    serial_scalar_array => vector%get_serial_scalar_array()
    
    if ( allocated(coarse_dofs_values) ) then
      call memfree ( coarse_dofs_values, __FILE__, __LINE__ )
    end if
    
    call memalloc ( this%get_block_number_coarse_dofs(1), &
                    coarse_dofs_values, &
                    __FILE__, &
                    __LINE__ );
    
    coarse_dofs_values = 0.0_rp
#ifdef ENABLE_BLAS
    call DGEMV(  'T', & 
                 serial_scalar_array%get_size(), &
                 this%get_block_number_coarse_dofs(1), &
                 1.0_rp, &
                 this%Phi, &
                 serial_scalar_array%get_size(), &
                 serial_scalar_array%get_entries(), &
                 1, &
                 0.0_rp, & 
                 coarse_dofs_values, & 
                 1)
#else
    write (0,*) 'Error: mlbddc.f90 was not compiled with -DENABLE_BLAS.'
    write (0,*) 'Error: You must activate this cpp macro in order to use the BLAS'
    check(.false.)    
#endif
  end subroutine mlbddc_compute_coarse_dofs_values
 
  subroutine mlbddc_compute_coarse_dofs_values_counts_and_displs ( this, counts, displs )
    implicit none
    class(mlbddc_t)           , intent(in)    :: this
    integer(ip), allocatable  , intent(inout) :: counts(:)
    integer(ip), allocatable  , intent(inout) :: displs(:)
    type(par_environment_t), pointer :: par_environment
    integer(ip) :: i, l1_to_l2_size
    type(coarse_fe_iterator_t) :: iterator
    type(coarse_fe_accessor_t) :: coarse_fe
    
    par_environment => this%get_par_environment()
    assert (par_environment%am_i_l1_to_l2_root())
    l1_to_l2_size = par_environment%get_l1_to_l2_size()
    if ( allocated(counts) ) call memfree ( counts, __FILE__, __LINE__ )
    if ( allocated(displs) ) call memfree ( displs, __FILE__, __LINE__ )
    call memalloc ( l1_to_l2_size, counts, __FILE__, __LINE__ )
    call memalloc ( l1_to_l2_size, displs, __FILE__, __LINE__ )
    
    i=1
    counts(l1_to_l2_size) = 0
    iterator = this%coarse_fe_space%create_coarse_fe_iterator()
    do while ( .not. iterator%has_finished() )
      coarse_fe = iterator%current()
      if ( coarse_fe%is_local() ) then
        counts(i) = coarse_fe%get_number_dofs()
        i=i+1
      end if
      call iterator%next()
    end do
    
    displs(1) = 0
    do i=2, l1_to_l2_size
      displs(i) = displs(i-1) + counts(i-1)
    end do
  end subroutine mlbddc_compute_coarse_dofs_values_counts_and_displs
  
  subroutine mlbddc_compute_and_gather_coarse_dofs_values ( this, vector, coarse_dofs_values_gathered )
     implicit none
     class(mlbddc_t)         , intent(in)    :: this
     type(par_scalar_array_t), intent(in)    :: vector
     real(rp)   , allocatable, intent(inout) :: coarse_dofs_values_gathered(:)
     type(par_environment_t), pointer :: par_environment
     integer(ip), allocatable :: counts(:)
     integer(ip), allocatable :: displs(:)
     real(rp), allocatable :: coarse_dofs_values(:)
     real(rp) :: dummy_real_array_rp(0)
     integer(ip) :: dummy_integer_array_ip(0)
     integer(ip) :: l1_to_l2_size
     par_environment => this%get_par_environment()
     assert (par_environment%am_i_l1_to_l2_task())
     
     if ( par_environment%am_i_l1_to_l2_root() ) then
       call this%compute_coarse_dofs_values_counts_and_displs(counts, displs)
       l1_to_l2_size = par_environment%get_l1_to_l2_size()
       
       if ( allocated(coarse_dofs_values_gathered) ) & 
         call memfree (coarse_dofs_values_gathered, __FILE__, __LINE__)
         
       call memalloc ( displs(l1_to_l2_size), & 
                       coarse_dofs_values_gathered, & 
                       __FILE__, __LINE__ ) 
       
       call par_environment%l2_from_l1_gather( input_data_size = 0, &
                                               input_data      = dummy_real_array_rp, &
                                               recv_counts     = counts, &
                                               displs          = displs, &
                                               output_data     = coarse_dofs_values_gathered )
       
       call memfree ( counts, __FILE__, __LINE__ )
       call memfree ( displs, __FILE__, __LINE__ )
     else
       call this%compute_coarse_dofs_values(vector, coarse_dofs_values)
       call par_environment%l2_from_l1_gather( input_data_size = size(coarse_dofs_values), &
                                               input_data      = coarse_dofs_values, &
                                               recv_counts     = dummy_integer_array_ip, &
                                               displs          = dummy_integer_array_ip, &
                                               output_data     = dummy_real_array_rp )
       call memfree ( coarse_dofs_values, __FILE__, __LINE__ )
     end if     
  end subroutine mlbddc_compute_and_gather_coarse_dofs_values
    
  ! This subroutine assumes that coarse_grid_vector has been already created+allocated, and initialized to zero.
  subroutine mlbddc_coarse_grid_residual_assembly ( this, coarse_dofs_values_gathered, coarse_grid_vector ) 
    implicit none
    class(mlbddc_t)         , intent(in)    :: this
    real(rp)                , intent(in)    :: coarse_dofs_values_gathered(*)
    type(par_scalar_array_t), intent(inout) :: coarse_grid_vector
     
    type(par_environment_t)   , pointer       :: L1_environment
    type(par_environment_t)   , pointer       :: L2_environment
    type(coarse_fe_iterator_t)                :: iterator
    type(coarse_fe_accessor_t)                :: coarse_fe
    
    type(i1p_t), allocatable :: elem2dof(:)
    integer(ip) :: ifield, i, istat, current
    
    L1_environment => this%get_par_environment()
    L2_environment => L1_environment%get_next_level()
    assert ( associated (L2_environment) )
    assert ( L2_environment%am_i_l1_task() )

    allocate ( elem2dof(this%coarse_fe_space%get_number_fields()), stat=istat)
    check(istat==0)     
    current = 1
    iterator  = this%coarse_fe_space%create_coarse_fe_iterator()
    do while ( .not. iterator%has_finished() )
       coarse_fe = iterator%current()
       call coarse_fe%get_elem2dof(elem2dof)
       if ( coarse_fe%is_local() ) then
          do ifield=1, this%coarse_fe_space%get_number_fields()
            do i=1, size(elem2dof(ifield)%p)
              call coarse_grid_vector%add(i   =  elem2dof(ifield)%p(i), &
                                          val = coarse_dofs_values_gathered(current) )
              current = current + 1
            end do
          end do
       end if
       call iterator%next()
    end do
    deallocate ( elem2dof, stat=istat )
    check(istat==0)
    call coarse_grid_vector%comm()
  end subroutine mlbddc_coarse_grid_residual_assembly
  
  subroutine mlbddc_scatter_and_interpolate_coarse_grid_correction (this, coarse_grid_vector, vector)
    implicit none
    class(mlbddc_t)           , intent(in)    :: this
    type(par_scalar_array_t)  , intent(in)    :: coarse_grid_vector
    type(par_scalar_array_t)  , intent(inout) :: vector
    type(par_environment_t)   , pointer :: L1_environment
    real(rp), allocatable :: coarse_dofs_values(:)
    
    L1_environment => this%get_par_environment()
    if ( L1_environment%am_i_l1_to_l2_task() ) then
      call this%scatter_coarse_grid_correction(coarse_grid_vector, coarse_dofs_values)
    end if
    
    if ( L1_environment%am_i_l1_task() ) then
      call this%interpolate_coarse_grid_correction(coarse_dofs_values, vector)
    end if
    
    if ( allocated(coarse_dofs_values) ) call memfree(coarse_dofs_values, __FILE__, __LINE__ )
  end subroutine mlbddc_scatter_and_interpolate_coarse_grid_correction
  
  subroutine mlbddc_scatter_coarse_grid_correction ( this, coarse_grid_vector, coarse_dofs_values )
    implicit none
    class(mlbddc_t)         , intent(in)    :: this
    type(par_scalar_array_t), intent(in)    :: coarse_grid_vector
    real(rp), allocatable   , intent(inout) :: coarse_dofs_values(:)
    
    integer(ip) :: l1_to_l2_size
    real(rp), allocatable :: coarse_dofs_values_scattered(:)
    integer(ip), allocatable :: counts(:), displs(:)
    type(par_environment_t), pointer :: par_environment
    real(rp) :: dummy_real_array_rp(0)
    integer(ip) :: dummy_integer_array_ip(0)
    
    par_environment => this%get_par_environment()
    assert (par_environment%am_i_l1_to_l2_task())
    if ( par_environment%am_i_l1_to_l2_root() ) then
       
       call this%compute_coarse_dofs_values_counts_and_displs(counts, displs)
       l1_to_l2_size = par_environment%get_l1_to_l2_size()
       
       call memalloc ( displs(l1_to_l2_size), & 
                       coarse_dofs_values_scattered, & 
                       __FILE__, __LINE__ ) 
       
       call this%fill_coarse_dofs_values_scattered ( coarse_grid_vector, &
                                                     coarse_dofs_values_scattered)
       
       call par_environment%l2_to_l1_scatter( input_data       = coarse_dofs_values_scattered, &
                                              send_counts      = counts, &
                                              displs           = displs, &
                                              output_data_size = 0, &
                                              output_data      = dummy_real_array_rp )
       
       call memfree ( coarse_dofs_values_scattered, & 
                       __FILE__, __LINE__ ) 
       
       call memfree ( counts, __FILE__, __LINE__ )
       call memfree ( displs, __FILE__, __LINE__ )
       
    else 
       if ( allocated (coarse_dofs_values) ) call memfree( coarse_dofs_values, __FILE__, __LINE__ )
       call memalloc ( this%get_block_number_coarse_dofs(1), &
                       coarse_dofs_values, &
                       __FILE__, __LINE__ );
       
       call par_environment%l2_to_l1_scatter( input_data       = dummy_real_array_rp, &
                                              send_counts      = dummy_integer_array_ip, &
                                              displs           = dummy_integer_array_ip, &
                                              output_data_size = size(coarse_dofs_values), &
                                              output_data      = coarse_dofs_values )
    end if
  end subroutine mlbddc_scatter_coarse_grid_correction
  
  subroutine mlbddc_fill_coarse_dofs_values_scattered ( this, coarse_grid_vector, coarse_dofs_values_scattered ) 
    implicit none
    class(mlbddc_t)          , intent(in)     :: this
    type(par_scalar_array_t) , intent(in)     :: coarse_grid_vector
    real(rp)                 , intent(inout)  :: coarse_dofs_values_scattered(*)
    type(par_environment_t), pointer    :: par_environment
    type(i1p_t)           , allocatable :: elem2dof(:)
    type(coarse_fe_iterator_t)          :: iterator
    type(coarse_fe_accessor_t)          :: coarse_fe
    integer(ip)                         :: istat, current, ifield
    
    par_environment => this%get_par_environment()
    assert (par_environment%am_i_l1_to_l2_root())
    
    allocate ( elem2dof(this%coarse_fe_space%get_number_fields()), stat=istat)
    check(istat==0)
    
    current = 1
    iterator  = this%coarse_fe_space%create_coarse_fe_iterator()
    do while ( .not. iterator%has_finished() )
       coarse_fe = iterator%current()
       call coarse_fe%get_elem2dof(elem2dof)
       if ( coarse_fe%is_local() ) then
          do ifield=1, this%coarse_fe_space%get_number_fields()
              call coarse_grid_vector%extract_subvector ( iblock       = 1, &
                                                          size_indices = size(elem2dof(ifield)%p), &
                                                          indices     = elem2dof(ifield)%p, &
                                                          values      = coarse_dofs_values_scattered(current) )
              current = current + size(elem2dof(ifield)%p)
          end do
       end if
       call iterator%next()
    end do
    
    deallocate ( elem2dof, stat=istat )
    check(istat==0)
  end subroutine mlbddc_fill_coarse_dofs_values_scattered 
    
  subroutine mlbddc_interpolate_coarse_grid_correction (this, coarse_dofs_values, vector)
    implicit none
    class(mlbddc_t)            , intent(in)    :: this
    real(rp)                   , intent(in)    :: coarse_dofs_values(*)
    type(par_scalar_array_t)   , intent(inout) :: vector
    type(par_environment_t)    , pointer       :: L1_environment
    type(serial_scalar_array_t), pointer       :: serial_scalar_array 
    real(rp)                   , pointer       :: serial_scalar_array_entries(:)
    
    L1_environment => this%get_par_environment()
    assert ( L1_environment%am_i_l1_task() )
    
    call vector%init(0.0_rp)
    serial_scalar_array         => vector%get_serial_scalar_array()
    serial_scalar_array_entries => serial_scalar_array%get_entries()
    
#ifdef ENABLE_BLAS
    call DGEMV(  'N', & 
                  serial_scalar_array%get_size(), &
                  this%get_block_number_coarse_dofs(1), &
                  1.0_rp, &
                  this%Phi, &
                  serial_scalar_array%get_size(), &
                  coarse_dofs_values, &
                  1,    &
                  0.0_rp,  & 
                  serial_scalar_array_entries, & 
                  1)
#else
     write (0,*) 'Error: mlbddc.f90 was not compiled with -DENABLE_BLAS.'
     write (0,*) 'Error: You must activate this cpp macro in order to use the BLAS'
     check(.false.)    
#endif
     
  end subroutine mlbddc_interpolate_coarse_grid_correction
  
  
  subroutine mlbddc_solve_dirichlet_problem(this, x_I, y_I)
    implicit none
    class(mlbddc_t)           , intent(in)    :: this
    type(par_scalar_array_t)  , intent(in)    :: x_I
    type(par_scalar_array_t)  , intent(inout) :: y_I
    type(serial_scalar_array_t), pointer      :: y_I_serial
    if ( this%am_i_l1_task() ) then
       y_I_serial => y_I%get_serial_scalar_array()
       call this%dirichlet_solver%solve(x_I%get_serial_scalar_array(), &
                                        y_I_serial)
    end if   
  end subroutine mlbddc_solve_dirichlet_problem  
  
  subroutine mlbddc_apply_A_GI(this, x_I, y_G)
    implicit none
    class(mlbddc_t)           , intent(in)    :: this
    type(par_scalar_array_t)  , intent(in)    :: x_I
    type(par_scalar_array_t)  , intent(inout) :: y_G
    type(par_sparse_matrix_t), pointer :: par_sparse_matrix
    type(sparse_matrix_t), pointer :: A
    type(serial_scalar_array_t), pointer    :: y_G_serial
    
    if ( this%am_i_l1_task() ) then
      par_sparse_matrix => this%get_par_sparse_matrix()
      A => par_sparse_matrix%get_sparse_matrix()
      y_G_serial => y_G%get_serial_scalar_array()
      if ( A%get_symmetric_storage() ) then
         call this%A_IG%apply_transpose(x_I%get_serial_scalar_array(), &
                                        y_G_serial)
      else
         call this%A_GI%apply(x_I%get_serial_scalar_array(), &
                              y_G_serial)
      end if
    end if
    
  end subroutine mlbddc_apply_A_GI
  
  subroutine mlbddc_apply_A_IG(this, x_G, y_I)
    implicit none
    class(mlbddc_t)           , intent(in)    :: this
    type(par_scalar_array_t)  , intent(in)    :: x_G
    type(par_scalar_array_t)  , intent(inout) :: y_I
    type(par_sparse_matrix_t), pointer :: par_sparse_matrix
    type(sparse_matrix_t), pointer :: A
    type(serial_scalar_array_t), pointer    :: y_I_serial
    if ( this%am_i_l1_task() ) then
      par_sparse_matrix => this%get_par_sparse_matrix()
      A => par_sparse_matrix%get_sparse_matrix()
      y_I_serial => y_I%get_serial_scalar_array()
      if ( A%get_symmetric_storage() ) then
         call this%A_IG%apply(x_G%get_serial_scalar_array(), &
                              y_I_serial)
      end if
    end if 
  end subroutine mlbddc_apply_A_IG
  
  subroutine mlbddc_solve_constrained_neumann_problem(this, x, y)
    implicit none
    class(mlbddc_t)            , intent(in)    :: this
    type(par_scalar_array_t)   , intent(in)    :: x
    type(par_scalar_array_t)   , intent(inout) :: y
    type(new_par_fe_space_t)       , pointer :: par_fe_space
    type(serial_scalar_array_t), pointer :: x_serial
    type(serial_scalar_array_t), pointer :: y_serial
    real(rp)                   , pointer :: y_serial_entries(:)
    type(serial_scalar_array_t)          :: augmented_x
    type(serial_scalar_array_t)          :: augmented_y
    real(rp), allocatable                :: augmented_x_entries(:)
    real(rp), allocatable                :: augmented_y_entries(:)
    integer(ip)                          :: block_number_dofs
    integer(ip)                          :: block_number_coarse_dofs

    if ( this%am_i_l1_task() ) then
      par_fe_space => this%get_fe_space()
      block_number_dofs        = par_fe_space%get_block_number_dofs(1)
      block_number_coarse_dofs = this%get_block_number_coarse_dofs(1)
      
      x_serial => x%get_serial_scalar_array()
    
      ! Set-up augmented_x from x_serial
      call augmented_x%create(block_number_dofs + block_number_coarse_dofs)
      call memalloc ( block_number_dofs + block_number_coarse_dofs, & 
                      augmented_x_entries, __FILE__,__LINE__)  
      augmented_x_entries(1:block_number_dofs)  = x_serial%get_entries() 
      augmented_x_entries(block_number_dofs+1:) = 0.0_rp 
      call augmented_x%set_view_entries(augmented_x_entries)

      ! Set-up augmented_y
      call augmented_y%create(block_number_dofs + block_number_coarse_dofs)
      call memalloc ( block_number_dofs + block_number_coarse_dofs, & 
                      augmented_y_entries, __FILE__,__LINE__)  
      call augmented_y%set_view_entries(augmented_y_entries)
      call this%constrained_neumann_solver%solve(augmented_x, &
                                                 augmented_y)

      ! Set-up y from augmented_y
      y_serial         => y%get_serial_scalar_array()
      y_serial_entries => y_serial%get_entries()
      y_serial_entries = augmented_y_entries(1:block_number_dofs)
      
      call memfree ( augmented_x_entries, __FILE__,__LINE__)  
      call memfree ( augmented_y_entries, __FILE__,__LINE__)  
    end if  
  end subroutine mlbddc_solve_constrained_neumann_problem
       
  subroutine mlbddc_apply_weight_operator(this, x, y)
    implicit none
    class(mlbddc_t)    , intent(in)    :: this
    type(par_scalar_array_t)  , intent(in)    :: x
    type(par_scalar_array_t)  , intent(inout) :: y
    type(serial_scalar_array_t), pointer :: x_local
    type(serial_scalar_array_t), pointer :: y_local
    real(rp), pointer :: x_local_entries(:)
    real(rp), pointer :: y_local_entries(:)
    type(dof_object_iterator_t) :: dofs_object_iterator
    type(dof_object_accessor_t) :: dof_object
    type(list_iterator_t) :: dofs_on_object
 
    if ( this%am_i_l1_task() ) then
      x_local         => x%get_serial_scalar_array()
      x_local_entries => x_local%get_entries()
      y_local         => y%get_serial_scalar_array()
      y_local_entries => y_local%get_entries()
      dofs_object_iterator = this%create_field_dofs_object_iterator(1)
      do while ( .not. dofs_object_iterator%has_finished() ) 
         call dofs_object_iterator%current(dof_object)
         dofs_on_object = dof_object%get_dofs_on_object_iterator()
         do while ( .not. dofs_on_object%is_upper_bound() )
           y_local_entries(dofs_on_object%get_current()) = &
             x_local_entries(dofs_on_object%get_current())/dof_object%get_number_parts_around()
           call dofs_on_object%next()
         end do
         call dofs_object_iterator%next()
      end do
    end if
  end subroutine mlbddc_apply_weight_operator

  subroutine mlbddc_create_interior_interface_views ( this, x, x_I, X_G )
    implicit none
    class(mlbddc_t)         , intent(in)       :: this
    type(par_scalar_array_t), intent(in)       :: x
    type(par_scalar_array_t), intent(inout)    :: x_I
    type(par_scalar_array_t), intent(inout)    :: x_G
    type(new_par_fe_space_t), pointer :: par_fe_space
    par_fe_space => this%get_fe_space()
    call x%create_view(1, &
                       par_fe_space%get_block_number_interior_dofs(1), &
                       x_I)
    call x%create_view(par_fe_space%get_block_number_interior_dofs(1)+1, &
                       par_fe_space%get_block_number_dofs(1), &
                       x_G)
  end subroutine mlbddc_create_interior_interface_views 
  
  
  subroutine mlbddc_free(this)
    implicit none
    class(mlbddc_t)           , intent(inout) :: this
    call this%free_numerical_setup()
    call this%free_symbolic_setup()
    call this%free_clean()
  end subroutine mlbddc_free
  
  subroutine mlbddc_free_clean(this)
    implicit none
    class(mlbddc_t)           , intent(inout) :: this
    nullify(this%fe_affine_operator)
    call this%free_vector_spaces()
  end subroutine mlbddc_free_clean
  
  subroutine mlbddc_free_symbolic_setup(this)
    implicit none
    class(mlbddc_t)           , intent(inout) :: this
    type(par_environment_t)   , pointer       :: par_environment
    integer(ip)                               :: istat
    
    ! This will be replaced by an apropriate check on the state diagram
    if ( associated(this%fe_affine_operator) ) then 
      par_environment => this%get_par_environment()
      if ( par_environment%am_i_l1_task() ) then
        if ( par_environment%get_l1_size() > 1 ) then  
          call this%free_dofs_objects_and_constraint_matrix()
          call this%free_symbolic_setup_dirichlet_solver()
          call this%free_symbolic_setup_dirichlet_problem()
          call this%free_symbolic_setup_constrained_neumann_solver()
          call this%free_symbolic_setup_constrained_neumann_problem()
          nullify (this%mlbddc_coarse)
          nullify (this%coarse_grid_matrix)
          nullify (this%coarse_fe_space)
        else
          call this%free_symbolic_setup_coarse_solver()
        end if    
      else
          ! mlbddc_coarse should be freed before coarse_fe_space
          ! (as the former was created from the latter)
          call this%mlbddc_coarse%free_symbolic_setup()
          call this%mlbddc_coarse%free_clean()
          deallocate (this%mlbddc_coarse, stat=istat)
          check (istat==0)
      
          call this%coarse_grid_matrix%free()
          deallocate  ( this%coarse_grid_matrix, stat = istat )
          check( istat == 0 )
        
          call this%coarse_fe_space%free()
          deallocate (this%coarse_fe_space, stat=istat)
          check (istat==0)
      end if
    end if  
  end subroutine mlbddc_free_symbolic_setup 
  
  subroutine mlbddc_free_dofs_objects_and_constraint_matrix(this)
    implicit none
    class(mlbddc_t)           , intent(inout) :: this
    integer(ip)                               :: i, istat 
    
    assert ( this%am_i_l1_task() )
    
    if (allocated(this%num_dofs_objects_per_field)) then 
      call memfree ( this%num_dofs_objects_per_field, __FILE__, __LINE__ )
    end if
  
    if (allocated(this%dofs_objects_per_field)) then
      do i=1, size(this%dofs_objects_per_field)
        call this%dofs_objects_per_field(i)%free()
      end do  
      deallocate ( this%dofs_objects_per_field, stat=istat )
      check (istat == 0)
    end if
  
    if (allocated(this%dofs_objects_gids_per_field)) then
      do i=1, size(this%dofs_objects_gids_per_field)
        call this%dofs_objects_gids_per_field(i)%free()
      end do  
      deallocate ( this%dofs_objects_gids_per_field, stat=istat )
      check (istat == 0)  
    end if
  
    if (allocated(this%vefs_lids_dofs_objects_per_field)) then
      do i=1, size(this%vefs_lids_dofs_objects_per_field)
        call this%vefs_lids_dofs_objects_per_field(i)%free()
      end do  
      deallocate ( this%vefs_lids_dofs_objects_per_field, stat=istat )
      check (istat == 0)  
    end if
    
    call this%constraint_matrix%free()
  end subroutine mlbddc_free_dofs_objects_and_constraint_matrix
  
  subroutine mlbddc_free_symbolic_setup_dirichlet_problem(this)
    implicit none
    class(mlbddc_t)           , intent(inout) :: this
    assert ( this%am_i_l1_task() )
    call this%A_II%free()
    call this%A_IG%free()
    call this%A_GI%free()
    call this%A_GG%free()
  end subroutine mlbddc_free_symbolic_setup_dirichlet_problem
  
  subroutine mlbddc_free_symbolic_setup_dirichlet_solver(this)
    implicit none
    class(mlbddc_t)           , intent(inout) :: this
    assert ( this%am_i_l1_task() )
    call this%dirichlet_solver%free()
  end subroutine mlbddc_free_symbolic_setup_dirichlet_solver
  
  subroutine mlbddc_free_symbolic_setup_constrained_neumann_problem(this)
    implicit none
    class(mlbddc_t), intent(inout) :: this
    assert ( this%am_i_l1_task() )
    call this%constrained_neumann_matrix%free()
  end subroutine mlbddc_free_symbolic_setup_constrained_neumann_problem
  
  subroutine mlbddc_free_symbolic_setup_constrained_neumann_solver(this)
    implicit none
    class(mlbddc_t)           , intent(inout) :: this
    assert ( this%am_i_l1_task() )
    call this%constrained_neumann_solver%free()
  end subroutine mlbddc_free_symbolic_setup_constrained_neumann_solver
  
  subroutine mlbddc_free_symbolic_setup_coarse_solver(this)
    implicit none
    class(mlbddc_t)           , intent(inout) :: this
    assert ( this%am_i_l1_task() )
    call this%coarse_solver%free()
  end subroutine mlbddc_free_symbolic_setup_coarse_solver
  
  subroutine mlbddc_free_numerical_setup(this)
    implicit none
    class(mlbddc_t)           , intent(inout) :: this
    type(par_environment_t)   , pointer       :: par_environment
    ! This will be replaced by an apropriate check on the state diagram
    if ( associated(this%fe_affine_operator) ) then 
      par_environment => this%get_par_environment()
      if ( par_environment%am_i_l1_task() ) then
        if ( par_environment%get_l1_size() > 1 ) then
          call this%free_numerical_setup_dirichlet_solver()
          call this%free_numerical_setup_dirichlet_problem()
          call this%free_numerical_setup_constrained_neumann_solver()
          call this%free_numerical_setup_constrained_neumann_problem()
          call this%free_coarse_grid_basis()
        else if ( par_environment%get_l1_size() == 1 ) then
          call this%free_numerical_setup_coarse_solver()
        end if
      else
        call this%mlbddc_coarse%free_numerical_setup()
        call this%coarse_grid_matrix%free_in_stages(free_numerical_setup)
      end if
    end if
  end subroutine mlbddc_free_numerical_setup
  
  subroutine mlbddc_free_numerical_setup_dirichlet_problem(this)
    implicit none
    class(mlbddc_t)           , intent(inout) :: this
    assert ( this%am_i_l1_task() )
    call this%A_II%free_in_stages(free_numerical_setup)
    call this%A_IG%free_in_stages(free_numerical_setup)
    call this%A_GI%free_in_stages(free_numerical_setup)
    call this%A_GG%free_in_stages(free_numerical_setup)
  end subroutine mlbddc_free_numerical_setup_dirichlet_problem
  
  subroutine mlbddc_free_numerical_setup_dirichlet_solver(this)
    implicit none
    class(mlbddc_t)           , intent(inout) :: this
    assert ( this%am_i_l1_task() )
    call this%dirichlet_solver%free_in_stages(free_numerical_setup)
  end subroutine mlbddc_free_numerical_setup_dirichlet_solver
  
  subroutine mlbddc_free_numerical_setup_constrained_neumann_problem(this)
    implicit none
    class(mlbddc_t)           , intent(inout) :: this
    assert ( this%am_i_l1_task() )
    call this%constrained_neumann_matrix%free_in_stages(free_numerical_setup)
  end subroutine mlbddc_free_numerical_setup_constrained_neumann_problem
  
  subroutine mlbddc_free_numerical_setup_constrained_neumann_solver(this)
    implicit none
    class(mlbddc_t)           , intent(inout) :: this
    assert ( this%am_i_l1_task() )
    call this%constrained_neumann_solver%free_in_stages(free_numerical_setup)
  end subroutine mlbddc_free_numerical_setup_constrained_neumann_solver
  
  subroutine mlbddc_free_coarse_grid_basis ( this ) 
    implicit none
    class(mlbddc_t)           , intent(inout) :: this
    assert ( this%am_i_l1_task() )
    if ( allocated ( this%Phi ) ) then
      call memfree ( this%Phi, __FILE__, __LINE__ )
    end if
  end subroutine mlbddc_free_coarse_grid_basis
  
  subroutine mlbddc_free_numerical_setup_coarse_solver ( this )
    implicit none
    class(mlbddc_t)           , intent(inout) :: this
    type(par_environment_t)  , pointer :: par_environment
    par_environment => this%get_par_environment()
    assert ( par_environment%get_l1_size() == 1 )
    assert ( par_environment%am_i_l1_task() )
    call this%coarse_solver%free_in_stages(free_numerical_setup)
  end subroutine mlbddc_free_numerical_setup_coarse_solver
  
  function mlbddc_get_total_number_coarse_dofs ( this )
    implicit none
    class(mlbddc_t)           , intent(in) :: this
    integer(ip)                               :: mlbddc_get_total_number_coarse_dofs
    type(new_par_fe_space_t)     , pointer        :: fe_space
    integer(ip)                               :: field_id
    
    assert ( this%am_i_l1_task() )
    
    mlbddc_get_total_number_coarse_dofs = 0
    fe_space => this%get_fe_space()
    do field_id = 1, fe_space%get_number_fields()
       mlbddc_get_total_number_coarse_dofs = mlbddc_get_total_number_coarse_dofs + & 
                                             this%num_dofs_objects_per_field(field_id)
    end do
  end function mlbddc_get_total_number_coarse_dofs
  
  function mlbddc_get_block_number_coarse_dofs ( this, block_id )
    implicit none
    class(mlbddc_t)           , intent(in)    :: this
    integer(ip)               , intent(in)    :: block_id
    integer(ip)                               :: mlbddc_get_block_number_coarse_dofs 
    type(new_par_fe_space_t)  , pointer       :: fe_space
    integer(ip)                               :: field_id
    integer(ip)               , pointer       :: field_blocks(:)
    assert ( this%am_i_l1_task() )
    assert ( block_id == 1 )
    
    fe_space     => this%get_fe_space()
    field_blocks => fe_space%get_field_blocks()
    
    mlbddc_get_block_number_coarse_dofs = 0
    do field_id = 1, fe_space%get_number_fields()
       if ( field_blocks(field_id) == block_id ) then
         mlbddc_get_block_number_coarse_dofs = mlbddc_get_block_number_coarse_dofs + & 
                                               this%num_dofs_objects_per_field(field_id)
       end if                                      
    end do
  end function mlbddc_get_block_number_coarse_dofs
  
  function mlbddc_create_field_dofs_object_iterator(this, field_id)
    implicit none
    class(mlbddc_t), intent(in) :: this
    integer(ip)    , intent(in) :: field_id
    type(dof_object_iterator_t) :: mlbddc_create_field_dofs_object_iterator
    call mlbddc_create_field_dofs_object_iterator%create(1, field_id, this)
  end function mlbddc_create_field_dofs_object_iterator
  
  
  
  ! Helper function that extracts a run-time polymorphic class(matrix_t)
  ! from the fe_affine_operator, and dynamically casts it into  
  ! type(par_sparse_matrix_t). If the dynamic cast cannot be performed 
  ! [because class(matrix_t) is NOT of type(par_sparse_matrix_t)], then it 
  ! aborts the execution of the program.
  function mlbddc_get_par_sparse_matrix(this)
    implicit none
    class(mlbddc_t)          , intent(in) :: this
    type(par_sparse_matrix_t), pointer    :: mlbddc_get_par_sparse_matrix
    class(matrix_t)          , pointer    :: matrix
    matrix => this%fe_affine_operator%get_matrix()
    select type(matrix)
    type is (par_sparse_matrix_t)
      mlbddc_get_par_sparse_matrix => matrix
    class default
      check(.false.)
    end select
  end function mlbddc_get_par_sparse_matrix
  
  ! Helper function that extracts a run-time polymorphic class(serial_fe_space_t)
  ! from the fe_affine_operator, and dynamically casts it into  
  ! type(par_fe_space_t). If the dynamic cast cannot be performed 
  ! [because class(serial_fe_space_t) is NOT of type(par_fe_space_t)], then it 
  ! aborts the execution of the program.
  function mlbddc_get_fe_space(this)
    implicit none
    class(mlbddc_t)          , intent(in) :: this
    type(new_par_fe_space_t)     , pointer    :: mlbddc_get_fe_space
    class(new_serial_fe_space_t) , pointer    :: fe_space
    fe_space => this%fe_affine_operator%get_fe_space()
    select type(fe_space)
    type is (new_par_fe_space_t)
      mlbddc_get_fe_space => fe_space
    class default
      check(.false.)
    end select
  end function mlbddc_get_fe_space
  
  function mlbddc_get_par_environment(this)
    implicit none
    class(mlbddc_t)          , intent(in) :: this
    type(par_environment_t)  , pointer    :: mlbddc_get_par_environment
    type(new_par_fe_space_t)     , pointer    :: coarse_fe_space
    coarse_fe_space => this%get_fe_space()
    mlbddc_get_par_environment => coarse_fe_space%get_par_environment()
  end function mlbddc_get_par_environment
  
  function mlbddc_am_i_l1_task(this)
    implicit none
    class(mlbddc_t)          , intent(in) :: this
    logical                               :: mlbddc_am_i_l1_task
    type(par_environment_t)   , pointer   :: par_environment
    par_environment => this%get_par_environment()
    mlbddc_am_i_l1_task = par_environment%am_i_l1_task()
  end function mlbddc_am_i_l1_task
  
  function mlbddc_is_linear( op )
    implicit none
    class(mlbddc_t)          , intent(in) :: op
    logical :: mlbddc_is_linear
    mlbddc_is_linear = .true.
  end function mlbddc_is_linear
  
  subroutine mlbddc_coarse_create ( this, fe_space, par_sparse_matrix )
    implicit none
    class(mlbddc_coarse_t)           , intent(inout) :: this
    type(coarse_fe_space_t)  , target, intent(in)    :: fe_space
    type(par_sparse_matrix_t), target, intent(in)    :: par_sparse_matrix
    call this%free()
    this%fe_space => fe_space
    this%par_sparse_matrix => par_sparse_matrix
  end subroutine mlbddc_coarse_create
  
  subroutine mlbddc_coarse_symbolic_setup ( this )
    implicit none
    class(mlbddc_coarse_t), intent(inout) :: this
    type(par_environment_t), pointer :: par_environment
    integer(ip)                               :: istat
    par_environment => this%get_par_environment()
    
    if( par_environment%am_i_l1_task() ) then
       if ( par_environment%get_l1_size() > 1 ) then
         call this%setup_dofs_objects_and_constraint_matrix()
       end if
    end if
       
    call this%setup_coarse_fe_space()
    call this%symbolic_setup_coarse_grid_matrix()
    call this%symbolic_setup_mlbddc_coarse()
    
    if ( par_environment%am_i_l1_task() ) then
      if ( par_environment%get_l1_size() > 1 ) then      
        call this%symbolic_setup_dirichlet_problem()
        call this%symbolic_setup_dirichlet_solver()
        call this%symbolic_setup_constrained_neumann_problem()
        call this%symbolic_setup_constrained_neumann_solver()
      else 
        call this%symbolic_setup_coarse_solver()
      end if
    end if
  end subroutine mlbddc_coarse_symbolic_setup
  
  subroutine mlbddc_coarse_setup_dofs_objects_and_constraint_matrix (this)
    implicit none
    class(mlbddc_coarse_t), intent(inout) :: this
    type(par_environment_t), pointer :: par_environment
    type(coarse_fe_space_t)   , pointer :: fe_space
    integer(ip)                      :: i 
    
    assert ( this%am_i_l1_task() )
    fe_space => this%get_fe_space()
    call fe_space%setup_dofs_objects_and_constraint_matrix(this%num_dofs_objects_per_field, &
                                                           this%dofs_objects_per_field, &
                                                           this%dofs_objects_gids_per_field, &
                                                           this%vefs_lids_dofs_objects_per_field, &
                                                           this%constraint_matrix)
    
    !write (*,'(a)') '****print mlbddc_coarse_setup_dofs_objects_and_constraint_matrix****'
    !write (*,'(a,i10)'  ) 'num_dofs_objects_per_field:', this%num_dofs_objects_per_field 
    !do i=1, fe_space%get_number_fields()
    !  write (*,'(a,i5)') '****FIELD ', i, '****'
    !  write (*,'(a)') 'dofs_objects_per_field'
    !  call this%dofs_objects_per_field(i)%print(6)
    !  call this%constraint_matrix%print(6)
    !end do
    !write (*,'(a)') '****print mlbddc_coarse_setup_dofs_objects_and_constraint_matrix****'
  end subroutine mlbddc_coarse_setup_dofs_objects_and_constraint_matrix
  
  subroutine mlbddc_coarse_setup_coarse_fe_space(this)
  implicit none
  class(mlbddc_coarse_t), intent(inout) :: this
  integer(ip)                           :: istat
  integer(ip)                           :: number_fields
  integer(ip), allocatable              :: field_type(:)
  integer(ip), allocatable              :: ptr_dofs_per_fe_and_field(:)
  integer(ip), allocatable              :: coarse_dofs_gids_recv_counts(:)
  integer(ip), allocatable              :: coarse_dofs_gids_displs(:)
  integer(igp), allocatable             :: lst_dofs_gids(:)
  integer(igp), allocatable             :: lst_vefs_gids_dofs_objects(:)
  type(par_environment_t), pointer      :: par_environment
  type(coarse_fe_space_t), pointer      :: fe_space
  type(coarse_triangulation_t), pointer :: coarse_triangulation
  
  
  par_environment => this%get_par_environment()
  fe_space        => this%get_fe_space()
  
   ! All MPI tasks (even if they are not involved in the L2 from L1 gather) should also allocate the
   ! allocatable arrays due to the fact that non-allocated allocatable arrays cannot
   ! be passed as actual arguments of dummy arguments that do not have the allocatable attribute 
   ! Otherwise, the code crashes with a segmentation fault.
   call memalloc (0, field_type, __FILE__, __LINE__)
   call memalloc (0, ptr_dofs_per_fe_and_field, __FILE__, __LINE__)
   call memalloc (0, lst_dofs_gids, __FILE__, __LINE__)
   call memalloc (0, lst_vefs_gids_dofs_objects, __FILE__, __LINE__)
   call memalloc (0, coarse_dofs_gids_recv_counts, __FILE__, __LINE__)
   call memalloc (0, coarse_dofs_gids_displs, __FILE__, __LINE__)

   
  ! L2 tasks gather from L1 tasks all raw data required to set-up the coarse triangulation on L2 tasks
  if ( par_environment%am_i_l1_to_l2_task() ) then
     call this%transfer_number_fields(number_fields) 
     call this%transfer_field_type(number_fields, field_type)
     call this%gather_ptr_dofs_per_fe_and_field(number_fields, ptr_dofs_per_fe_and_field)
     call this%gather_coarse_dofs_gids_rcv_counts_and_displs (coarse_dofs_gids_recv_counts, coarse_dofs_gids_displs)
     call this%gather_coarse_dofs_gids(coarse_dofs_gids_recv_counts, coarse_dofs_gids_displs, lst_dofs_gids)
     call this%gather_vefs_gids_dofs_objects(coarse_dofs_gids_recv_counts, coarse_dofs_gids_displs, lst_vefs_gids_dofs_objects)
  end if
  

  if ( par_environment%am_i_lgt1_task() ) then
     ! lgt1 MPI tasks (recursively) build coarse triangulation
     allocate  ( this%coarse_fe_space, stat = istat )
     check( istat == 0 )
     coarse_triangulation => fe_space%get_triangulation()
     call this%coarse_fe_space%create (coarse_triangulation%get_coarse_triangulation(), &
                                       number_fields, &
                                       field_type, &
                                       ptr_dofs_per_fe_and_field, &
                                       lst_dofs_gids, &
                                       lst_vefs_gids_dofs_objects)
  else
     ! L1 tasks do not hold any piece of the coarse triangulation
     nullify(this%coarse_fe_space)
  end if

  ! All tasks free raw data (see actual reason on the top part of this subroutine)
  call memfree (field_type, __FILE__, __LINE__)
  call memfree (ptr_dofs_per_fe_and_field, __FILE__, __LINE__)
  call memfree (lst_dofs_gids, __FILE__, __LINE__)
  call memfree (lst_vefs_gids_dofs_objects, __FILE__, __LINE__)
  call memfree (coarse_dofs_gids_recv_counts, __FILE__, __LINE__)
  call memfree (coarse_dofs_gids_displs, __FILE__, __LINE__)
end subroutine mlbddc_coarse_setup_coarse_fe_space

subroutine mlbddc_coarse_transfer_number_fields ( this, number_fields )
  implicit none
  class(mlbddc_coarse_t)   , intent(in)      :: this
  integer(ip)              , intent(out)     :: number_fields
  integer(ip)                                :: dummy_integer_ip
  type(par_environment_t), pointer           :: par_environment
  type(coarse_fe_space_t), pointer           :: fe_space

  par_environment => this%get_par_environment()
  assert ( par_environment%am_i_l1_to_l2_task() )
  fe_space        => this%get_fe_space()
  if ( par_environment%am_i_l1_to_l2_root() ) then
     call par_environment%l1_to_l2_transfer(input_data=dummy_integer_ip, &
                                            output_data=number_fields)
  else
     number_fields = fe_space%get_number_fields()
     call par_environment%l1_to_l2_transfer(input_data=number_fields, &
                                            output_data=dummy_integer_ip) 
  end if
end subroutine mlbddc_coarse_transfer_number_fields

subroutine mlbddc_coarse_transfer_field_type ( this, number_fields, field_type )
  implicit none
  class(mlbddc_coarse_t)   , intent(in)    :: this
  integer(ip)             , intent(in)       :: number_fields
  integer(ip), allocatable, intent(inout)    :: field_type(:)
  integer(ip)                                :: dummy_integer_array_ip(0)
  type(par_environment_t), pointer           :: par_environment
  type(coarse_fe_space_t), pointer           :: fe_space

  par_environment => this%get_par_environment()
  assert ( par_environment%am_i_l1_to_l2_task() )
  fe_space        => this%get_fe_space()
  if ( par_environment%am_i_l1_to_l2_root() ) then
     if ( allocated (field_type) ) call memfree ( field_type, __FILE__, __LINE__ )
     call memalloc ( number_fields, field_type, __FILE__, __LINE__ )
     call par_environment%l1_to_l2_transfer(input_data=dummy_integer_array_ip, &
                                            output_data=field_type)
  else
     call par_environment%l1_to_l2_transfer(input_data=fe_space%get_field_type(), &
                                            output_data=dummy_integer_array_ip) 
  end if
end subroutine mlbddc_coarse_transfer_field_type

subroutine mlbddc_coarse_gather_ptr_dofs_per_fe_and_field( this, number_fields, ptr_dofs_per_fe_and_field )
  implicit none
  class(mlbddc_coarse_t)   , intent(in)      :: this
  integer(ip)             , intent(in)       :: number_fields
  integer(ip), allocatable, intent(inout)    :: ptr_dofs_per_fe_and_field(:)
  integer(ip)                                :: i, num_local_cells
  integer(ip)                                :: dummy_integer_array(0)
  type(par_environment_t), pointer           :: par_environment
  type(coarse_fe_space_t), pointer           :: fe_space
  type(coarse_triangulation_t), pointer      :: triangulation
  type(coarse_triangulation_t), pointer      :: coarse_triangulation
  
  par_environment => this%get_par_environment()
  assert ( par_environment%am_i_l1_to_l2_task() )
  fe_space => this%get_fe_space()
  if ( par_environment%am_i_l1_to_l2_root() ) then
     triangulation => fe_space%get_triangulation()
     coarse_triangulation => triangulation%get_coarse_triangulation()
     num_local_cells = coarse_triangulation%get_num_local_cells()
     if (allocated(ptr_dofs_per_fe_and_field)) call memfree ( ptr_dofs_per_fe_and_field, __FILE__, __LINE__ )
     call memalloc (num_local_cells*number_fields+1, ptr_dofs_per_fe_and_field, __FILE__, __LINE__ )
     call par_environment%l2_from_l1_gather( input_data_size = number_fields, &
                                             input_data      = dummy_integer_array, &
                                             output_data     = ptr_dofs_per_fe_and_field(2:))
     ptr_dofs_per_fe_and_field(1) = 1
     do i=1, num_local_cells*number_fields
       ptr_dofs_per_fe_and_field(i+1) = ptr_dofs_per_fe_and_field(i) + ptr_dofs_per_fe_and_field(i+1)
     end do
  else
     call par_environment%l2_from_l1_gather( input_data_size = fe_space%get_number_fields(), &
                                             input_data      = this%num_dofs_objects_per_field, &
                                             output_data     = dummy_integer_array )
  end if
end subroutine mlbddc_coarse_gather_ptr_dofs_per_fe_and_field

  subroutine mlbddc_coarse_gather_coarse_dofs_gids_rcv_counts_and_displs( this, recv_counts, displs )
    implicit none
    class(mlbddc_coarse_t)     , intent(in)    :: this
    integer(ip) , allocatable , intent(inout) :: recv_counts(:) 
    integer(ip) , allocatable , intent(inout) :: displs(:)
    integer(ip)                               :: i
    integer(ip)                               :: l1_to_l2_size
    integer(ip)                               :: dummy_integer_array(0)
    type(par_environment_t), pointer          :: par_environment
    type(coarse_fe_space_t), pointer          :: fe_space

    par_environment => this%get_par_environment()
    assert ( par_environment%am_i_l1_to_l2_task() )
    fe_space => this%get_fe_space()
    if ( par_environment%am_i_l1_to_l2_root() ) then
      l1_to_l2_size = par_environment%get_l1_to_l2_size()
      if ( allocated (recv_counts) ) call memfree ( recv_counts, __FILE__, __LINE__ )
      if ( allocated (displs) ) call memfree ( displs, __FILE__, __LINE__ )
      call memalloc ( l1_to_l2_size, recv_counts, __FILE__, __LINE__ )
      call memalloc ( l1_to_l2_size, displs, __FILE__, __LINE__ )
      call par_environment%l2_from_l1_gather( input_data = 0, &
                                                           output_data = recv_counts ) 
      displs(1) = 0
      do i=2, l1_to_l2_size
        displs(i) = displs(i-1) + recv_counts(i-1)
      end do
    else
      call par_environment%l2_from_l1_gather( input_data  = sum(this%num_dofs_objects_per_field), &
                                                           output_data = dummy_integer_array ) 
    end if
  end subroutine mlbddc_coarse_gather_coarse_dofs_gids_rcv_counts_and_displs
  
  subroutine mlbddc_coarse_gather_coarse_dofs_gids ( this, recv_counts, displs, lst_gids )
    implicit none
    class(mlbddc_coarse_t)     , intent(in)    :: this
    integer(ip)               , intent(in)    :: recv_counts(:)
    integer(ip)               , intent(in)    :: displs(:)
    integer(igp), allocatable , intent(inout) :: lst_gids(:)
    integer(ip)                               :: l1_to_l2_size
    integer(igp)                              :: dummy_integer_array_igp(0)
    integer(ip)                               :: dummy_integer_array_ip(0)
    type(par_environment_t), pointer          :: par_environment
    type(coarse_fe_space_t), pointer          :: fe_space
    integer(ip)                               :: i, spos, epos
    integer(igp), allocatable                 :: buffer(:)
    
    par_environment => this%get_par_environment()
    assert ( par_environment%am_i_l1_to_l2_task() )
    fe_space => this%get_fe_space()
    if ( par_environment%am_i_l1_to_l2_root() ) then
      l1_to_l2_size = par_environment%get_l1_to_l2_size()
      if (allocated(lst_gids)) call memfree ( lst_gids, __FILE__, __LINE__ )
      call memalloc ( displs(l1_to_l2_size), lst_gids, __FILE__, __LINE__ )
      call par_environment%l2_from_l1_gather( input_data_size = 0, &
                                              input_data      = dummy_integer_array_igp, &
                                              recv_counts     = recv_counts, &
                                              displs          = displs, &
                                              output_data     = lst_gids )
    else
      ! Pack dofs_objects_gids_per_field(:) into plain buffer for further data exchange
      call memalloc ( sum(this%num_dofs_objects_per_field), buffer, __FILE__, __LINE__ ) 
      spos = 1
      do i=1, fe_space%get_number_fields()
        epos = spos + this%num_dofs_objects_per_field(i)-1
        buffer(spos:epos) = this%dofs_objects_gids_per_field(i)%a
        spos = epos +1 
      end do
    
      call par_environment%l2_from_l1_gather( input_data_size = size(buffer), &
                                              input_data      = buffer, &
                                              recv_counts     = dummy_integer_array_ip, &
                                              displs          = dummy_integer_array_ip, &
                                              output_data     = dummy_integer_array_igp )
      call memfree ( buffer, __FILE__, __LINE__ )
    end if    
  end subroutine mlbddc_coarse_gather_coarse_dofs_gids

  subroutine mlbddc_coarse_gather_vefs_gids_dofs_objects ( this, recv_counts, displs, vef_gids )
    implicit none
    class(mlbddc_coarse_t)     , intent(in)   :: this
    integer(ip)               , intent(in)    :: recv_counts(:)
    integer(ip)               , intent(in)    :: displs(:)
    integer(igp), allocatable , intent(inout) :: vef_gids(:)
    integer(ip)                               :: l1_to_l2_size
    integer(igp)                              :: dummy_integer_array_igp(0)
    integer(ip)                               :: dummy_integer_array_ip(0)
    type(par_environment_t), pointer          :: par_environment
    type(coarse_fe_space_t), pointer          :: fe_space
    integer(ip)                               :: i, j, spos, epos
    integer(igp), allocatable                 :: buffer(:)
    type(coarse_dof_object_iterator_t)        :: coarse_dofs_object_iterator
    type(coarse_dof_object_accessor_t)        :: coarse_dof_object
    
    par_environment => this%get_par_environment()
    assert ( par_environment%am_i_l1_to_l2_task() )
    fe_space => this%get_fe_space()
    if ( par_environment%am_i_l1_to_l2_root() ) then
      l1_to_l2_size = par_environment%get_l1_to_l2_size()
      if (allocated(vef_gids)) call memfree ( vef_gids, __FILE__, __LINE__ )
      call memalloc ( displs(l1_to_l2_size), vef_gids, __FILE__, __LINE__ )
      call par_environment%l2_from_l1_gather( input_data_size = 0, &
                                                           input_data      = dummy_integer_array_igp, &
                                                           recv_counts     = recv_counts, &
                                                           displs          = displs, &
                                                           output_data     = vef_gids )
    else
      ! Pack vefs_lids_dofs_objects_per_field(:) into plain buffer for further data exchange
      call memalloc ( sum(this%num_dofs_objects_per_field), buffer, __FILE__, __LINE__ ) 
      spos = 1
      do i=1, fe_space%get_number_fields()
        epos = spos + this%num_dofs_objects_per_field(i)-1
        coarse_dofs_object_iterator = this%create_field_dofs_object_iterator(i)
        do j=spos, epos
          coarse_dof_object = coarse_dofs_object_iterator%current()
          buffer(j) = coarse_dof_object%get_gid()
          call coarse_dofs_object_iterator%next()
        end do  
        spos = epos +1 
      end do
    
      call par_environment%l2_from_l1_gather( input_data_size = size(buffer), &
                                              input_data      = buffer, &
                                              recv_counts     = dummy_integer_array_ip, &
                                              displs          = dummy_integer_array_ip, &
                                              output_data     = dummy_integer_array_igp )
      
      call memfree ( buffer, __FILE__, __LINE__ )
    end if    
  end subroutine mlbddc_coarse_gather_vefs_gids_dofs_objects
  
  subroutine mlbddc_coarse_symbolic_setup_dirichlet_problem ( this) 
    implicit none
    class(mlbddc_coarse_t)    , intent(inout) :: this
    type(par_environment_t)   , pointer       :: par_environment
    type(coarse_fe_space_t)   , pointer       :: fe_space
    type(par_sparse_matrix_t) , pointer       :: par_sparse_matrix
    type(sparse_matrix_t)     , pointer       :: A
    
    assert ( this%am_i_l1_task() )
    fe_space => this%get_fe_space()
    par_sparse_matrix => this%get_par_sparse_matrix()
    A => par_sparse_matrix%get_sparse_matrix()
    if ( A%get_symmetric_storage() ) then
       call A%split_2x2_symbolic( num_row = fe_space%get_block_number_interior_dofs(1), &
                                  num_col = fe_space%get_block_number_interior_dofs(1), &
                                  A_II    = this%A_II, &
                                  A_IG    = this%A_IG, &
                                  A_GG    = this%A_GG)
    else
       call A%split_2x2_symbolic( num_row = fe_space%get_block_number_interior_dofs(1), &
                                  num_col = fe_space%get_block_number_interior_dofs(1), &
                                  A_II    = this%A_II, &
                                  A_IG    = this%A_IG, &
                                  A_GI    = this%A_GI, &
                                  A_GG    = this%A_GG)
    end if
  end subroutine mlbddc_coarse_symbolic_setup_dirichlet_problem
  
    subroutine mlbddc_coarse_symbolic_setup_dirichlet_solver ( this) 
    implicit none
    class(mlbddc_coarse_t), intent(inout) :: this
    assert ( this%am_i_l1_task() )
    call this%dirichlet_solver%set_type(name   = pardiso_mkl)
    call this%dirichlet_solver%set_matrix(matrix = this%A_II)
    call this%dirichlet_solver%symbolic_setup() 
  end subroutine mlbddc_coarse_symbolic_setup_dirichlet_solver
  
  subroutine mlbddc_coarse_symbolic_setup_constrained_neumann_problem(this)
    implicit none
    class(mlbddc_coarse_t)           , intent(inout) :: this
    type(par_environment_t)   , pointer       :: par_environment
    type(par_sparse_matrix_t) , pointer       :: par_sparse_matrix
    type(sparse_matrix_t)     , pointer       :: A
    
    assert ( this%am_i_l1_task() )
    par_sparse_matrix => this%get_par_sparse_matrix()
    A => par_sparse_matrix%get_sparse_matrix()
    
    call A%expand_matrix_symbolic(C_T               = this%constraint_matrix, &
                                  to                = this%constrained_neumann_matrix, &
                                  symmetric_storage = A%get_symmetric_storage(), & 
                                  symmetric         = A%is_symmetric(), &
                                  sign              = SPARSE_MATRIX_SIGN_INDEFINITE)
    
  end subroutine mlbddc_coarse_symbolic_setup_constrained_neumann_problem
  
  subroutine mlbddc_coarse_symbolic_setup_constrained_neumann_solver(this)
    implicit none
    class(mlbddc_coarse_t)           , intent(inout) :: this
    assert ( this%am_i_l1_task() )
    call this%constrained_neumann_solver%set_type(name   = pardiso_mkl)
    call this%constrained_neumann_solver%set_matrix(matrix = this%constrained_neumann_matrix)
    call this%constrained_neumann_solver%symbolic_setup() 
  end subroutine mlbddc_coarse_symbolic_setup_constrained_neumann_solver
  
  subroutine mlbddc_coarse_symbolic_setup_coarse_grid_matrix ( this )
    implicit none
    class(mlbddc_coarse_t)           , intent(inout) :: this
    type(par_environment_t)   , pointer       :: L1_environment
    type(par_environment_t)   , pointer       :: L2_environment
    integer(ip)                               :: istat
    
    L1_environment => this%get_par_environment()
    if ( L1_environment%am_i_lgt1_task() ) then
       ! lgt1 MPI tasks symbolically set-up this%coarse_grid_matrix
       allocate  ( this%coarse_grid_matrix, stat = istat )
       check( istat == 0 )

       L2_environment => L1_environment%get_next_level()
       call this%coarse_grid_matrix%create( p_env             = L2_environment, &
                                            dof_import        = this%coarse_fe_space%get_block_dof_import(1), &
                                            symmetric_storage = .true., &
                                            is_symmetric      = .true., &
                                            sign              = SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE )
       
       if ( L2_environment%am_i_l1_task() ) then
          call this%coarse_grid_matrix_symbolic_assembly()
       end if
    else       
       ! L1 tasks do not hold any piece of the coarse triangulation
       nullify(this%coarse_grid_matrix)
    end if
  end subroutine mlbddc_coarse_symbolic_setup_coarse_grid_matrix
  
  subroutine mlbddc_coarse_coarse_grid_matrix_symbolic_assembly ( this )
    implicit none
    class(mlbddc_coarse_t)           , intent(inout) :: this
    type(par_environment_t)   , pointer       :: L1_environment
    type(par_environment_t)   , pointer       :: L2_environment
    type(coarse_fe_iterator_t) :: iterator
    type(coarse_fe_accessor_t) :: coarse_fe
    type(i1p_t), allocatable :: elem2dof(:)
    logical, pointer :: field_coupling(:,:)
    integer(ip) :: ifield, jfield, i, j, istat
    
    L1_environment => this%get_par_environment()
    L2_environment => L1_environment%get_next_level()
    assert ( associated (L2_environment) )
    assert ( L2_environment%am_i_l1_task() )
    
    allocate ( elem2dof(this%coarse_fe_space%get_number_fields()), stat=istat)
    check(istat==0)     
    
    field_coupling => this%coarse_fe_space%get_field_coupling()
    
    iterator  = this%coarse_fe_space%create_coarse_fe_iterator()
    do while ( .not. iterator%has_finished() )
       coarse_fe = iterator%current()
       call coarse_fe%get_elem2dof(elem2dof)
       do ifield=1, this%coarse_fe_space%get_number_fields()
          do jfield=1, this%coarse_fe_space%get_number_fields()
             if ((field_coupling(ifield,jfield))) then
               do j=1, size(elem2dof(jfield)%p)
                 do i=1, size(elem2dof(ifield)%p)
                    call this%coarse_grid_matrix%insert(ia=elem2dof(ifield)%p(i), &
                                                        ja=elem2dof(jfield)%p(j) )
                 end do
               end do  
             end if
          end do
       end do
       call iterator%next()
    end do
    
    deallocate ( elem2dof, stat=istat )
    check(istat==0)     
    
    call this%coarse_grid_matrix%convert(csr_format)
    
    !call this%coarse_grid_matrix%print(6)
  end subroutine mlbddc_coarse_coarse_grid_matrix_symbolic_assembly
  
  subroutine mlbddc_coarse_symbolic_setup_mlbddc_coarse(this)
    implicit none
    class(mlbddc_coarse_t), intent(inout)     :: this
    type(par_environment_t)  , pointer :: par_environment
    integer(ip)                        :: istat
    
    par_environment => this%get_par_environment()
    if ( par_environment%am_i_lgt1_task() ) then
     ! lgt1 MPI tasks symbolically setup mlbddc coarse
     allocate  ( this%mlbddc_coarse, stat = istat )
     check( istat == 0 )
     call this%mlbddc_coarse%create(this%coarse_fe_space, &
                                    this%coarse_grid_matrix)
     call this%mlbddc_coarse%symbolic_setup()
    else
      ! L1 tasks do not hold any piece of mlbddc_coarse
      nullify(this%mlbddc_coarse)
     end if
  end subroutine mlbddc_coarse_symbolic_setup_mlbddc_coarse

  subroutine mlbddc_coarse_symbolic_setup_coarse_solver(this)
    implicit none
    class(mlbddc_coarse_t), intent(inout) :: this
    type(par_environment_t)  , pointer    :: par_environment
    type(par_sparse_matrix_t), pointer :: par_sparse_matrix
    par_environment => this%get_par_environment()
    assert ( par_environment%get_l1_size() == 1 )
    assert ( par_environment%am_i_l1_task() )
    par_sparse_matrix => this%get_par_sparse_matrix()
    call this%coarse_solver%set_type(name   = pardiso_mkl)
    call this%coarse_solver%set_matrix(matrix = par_sparse_matrix%get_sparse_matrix())
    call this%coarse_solver%symbolic_setup() 
  end subroutine mlbddc_coarse_symbolic_setup_coarse_solver

  
  subroutine mlbddc_coarse_numerical_setup ( this )
    implicit none
    class(mlbddc_coarse_t)        , intent(inout)   :: this
    type(par_environment_t), pointer         :: par_environment
    type(par_scalar_array_t) :: dum
    
    par_environment => this%get_par_environment()
    if ( par_environment%am_i_l1_task() ) then
      if ( par_environment%get_l1_size() > 1 ) then  
        call this%numerical_setup_constrained_neumann_problem()
        call this%numerical_setup_constrained_neumann_solver()
        call this%setup_coarse_grid_basis()
      end if
    end if  
    
    call this%numerical_setup_coarse_grid_matrix()
    call this%numerical_setup_mlbddc_coarse()
    
    if ( par_environment%am_i_l1_task() ) then
      if ( par_environment%get_l1_size() > 1 ) then  
        call this%numerical_setup_dirichlet_problem()
        call this%numerical_setup_dirichlet_solver()
      else if ( par_environment%get_l1_size() == 1 ) then
        call this%numerical_setup_coarse_solver()
      end if
    end if 
    
  end subroutine mlbddc_coarse_numerical_setup
  
  
  subroutine mlbddc_coarse_numerical_setup_dirichlet_problem (this) 
    implicit none
    class(mlbddc_coarse_t)           , intent(inout) :: this
    type(coarse_fe_space_t)      , pointer       :: fe_space
    type(par_sparse_matrix_t) , pointer       :: par_sparse_matrix
    type(sparse_matrix_t)     , pointer       :: A
    
    assert ( this%am_i_l1_task() )
    fe_space => this%get_fe_space()
    par_sparse_matrix => this%get_par_sparse_matrix()
    A => par_sparse_matrix%get_sparse_matrix()
    if ( A%get_symmetric_storage() ) then
       call A%split_2x2_numeric(  num_row = fe_space%get_block_number_interior_dofs(1), &
                                  num_col = fe_space%get_block_number_interior_dofs(1), &
                                  A_II    = this%A_II, &
                                  A_IG    = this%A_IG, &
                                  A_GG    = this%A_GG)
    else
       call A%split_2x2_numeric(  num_row = fe_space%get_block_number_interior_dofs(1), &
                                  num_col = fe_space%get_block_number_interior_dofs(1), &
                                  A_II    = this%A_II, &
                                  A_IG    = this%A_IG, &
                                  A_GI    = this%A_GI, &
                                  A_GG    = this%A_GG)
    end if
    
    !call A%print(6)
    !call this%A_II%print(6)
    !call this%A_IG%print(6)
    !call this%A_GG%print(6)
  end subroutine mlbddc_coarse_numerical_setup_dirichlet_problem
  
  subroutine mlbddc_coarse_numerical_setup_dirichlet_solver(this)
    implicit none
    class(mlbddc_coarse_t)           , intent(inout) :: this
    assert ( this%am_i_l1_task() )
    call this%dirichlet_solver%numerical_setup() 
    !call this%dirichlet_solver%log_info()
  end subroutine mlbddc_coarse_numerical_setup_dirichlet_solver
  
  subroutine mlbddc_coarse_numerical_constrained_neumann_problem(this)
   implicit none
   class(mlbddc_coarse_t)           , intent(inout) :: this
   type(coarse_fe_space_t)      , pointer       :: fe_space
   type(par_sparse_matrix_t) , pointer       :: par_sparse_matrix
   type(sparse_matrix_t)     , pointer       :: A
    
    assert ( this%am_i_l1_task() )
    fe_space => this%get_fe_space()
    par_sparse_matrix => this%get_par_sparse_matrix()
    A => par_sparse_matrix%get_sparse_matrix()
    call A%expand_matrix_numeric(C_T = this%constraint_matrix, &
                                 to  = this%constrained_neumann_matrix)
    
    !call this%constrained_neumann_matrix%print(6)
  end subroutine  mlbddc_coarse_numerical_constrained_neumann_problem
  
  subroutine mlbddc_coarse_numerical_setup_constrained_neumann_solver(this)
    implicit none
    class(mlbddc_coarse_t)           , intent(inout) :: this
    assert ( this%am_i_l1_task() )
    call this%constrained_neumann_solver%numerical_setup() 
    !call this%constrained_neumann_solver%log_info()
  end subroutine mlbddc_coarse_numerical_setup_constrained_neumann_solver
  
  subroutine mlbddc_coarse_allocate_coarse_grid_basis ( this )
    implicit none
    class(mlbddc_coarse_t)           , intent(inout) :: this
    type(coarse_fe_space_t), pointer       :: fe_space
    assert ( this%am_i_l1_task() )
    fe_space => this%get_fe_space() 
    call this%free_coarse_grid_basis()
    call memalloc (fe_space%get_block_number_dofs(1), &
                   this%get_block_number_coarse_dofs(1), &
                   this%Phi, &
                   __FILE__, &
                   __LINE__) 
  end subroutine mlbddc_coarse_allocate_coarse_grid_basis
  
  subroutine mlbddc_coarse_setup_coarse_grid_basis ( this )
     implicit none
     class(mlbddc_coarse_t)           , intent(inout) :: this
     type(coarse_fe_space_t), pointer :: fe_space
     real(rp), allocatable :: work1(:,:)
     real(rp), allocatable :: work2(:,:)
     integer(ip) :: i, j
     
     assert ( this%am_i_l1_task() )
     fe_space => this%get_fe_space() 
     
     call memalloc ( this%constrained_neumann_matrix%get_num_rows(), &
                     this%get_block_number_coarse_dofs(1), &
                     work1, __FILE__,__LINE__ )
     
     call memalloc ( this%constrained_neumann_matrix%get_num_rows(), &
                     this%get_block_number_coarse_dofs(1), &
                     work2, __FILE__,__LINE__ )
     
     work1 = 0.0_rp
     work2 = 0.0_rp
     j=1
     do i = fe_space%get_block_number_dofs(1)+1, this%constrained_neumann_matrix%get_num_rows()
       work1 (i,j) = 1.0_rp
       j = j + 1 
     end do
     
     call this%constrained_neumann_solver%solve(work1, &
                                                work2)
     
     call this%allocate_coarse_grid_basis()
     this%Phi = work2 (1:fe_space%get_block_number_dofs(1),:) 
     
     !do, i=1,this%constrained_neumann_matrix%get_num_rows()
     !  write(*,"(100g15.5)") ( work1(i,j), j=1,this%get_block_number_coarse_dofs(1) )
     !enddo
     !do, i=1,this%constrained_neumann_matrix%get_num_rows()
     !  write(*,"(100g15.5)") ( work2(i,j), j=1,this%get_block_number_coarse_dofs(1) )
     !enddo
     call memfree ( work1, __FILE__,__LINE__ )
     call memfree ( work2, __FILE__,__LINE__ )
  end subroutine mlbddc_coarse_setup_coarse_grid_basis 
  
  ! Computes subdomain_elmat = \Phi^t A_i \Phi 
  subroutine mlbddc_coarse_compute_subdomain_elmat ( this, subdomain_elmat )
    implicit none
    class(mlbddc_coarse_t)      , intent(in)    :: this
    real(rp), allocatable, intent(inout) :: subdomain_elmat(:,:)
    
    real(rp), allocatable :: work(:,:)
    type(par_sparse_matrix_t), pointer :: par_sparse_matrix
    type(sparse_matrix_t), pointer :: A
    
    assert ( this%am_i_l1_task() )
    par_sparse_matrix => this%get_par_sparse_matrix()
    A => par_sparse_matrix%get_sparse_matrix()
    
    if ( allocated(subdomain_elmat) ) then
      call memfree ( subdomain_elmat, __FILE__, __LINE__ )
    end if
    
    call memalloc ( this%get_block_number_coarse_dofs(1), &
                    this%get_block_number_coarse_dofs(1), &
                    subdomain_elmat, &
                    __FILE__, &
                    __LINE__ );
    
    call memalloc ( A%get_num_rows(), & 
                    this%get_block_number_coarse_dofs(1), &
                    work, &
                    __FILE__, & 
                    __LINE__ )

    work = 0.0_rp
    call A%apply_to_dense_matrix ( n     = this%get_block_number_coarse_dofs(1), &
                                   alpha = 1.0_rp, &
                                   ldb   = A%get_num_rows(), &
                                   b     = this%Phi, &
                                   beta  = 0.0_rp, &
                                   ldc   = A%get_num_rows(), &
                                   c     = work )        
    subdomain_elmat = 0.0_rp
#ifdef ENABLE_BLAS
    call DGEMM( 'T', &
                'N', &
                this%get_block_number_coarse_dofs(1), &
                this%get_block_number_coarse_dofs(1), &
                A%get_num_rows(), &
                1.0, &
                this%Phi, &
                A%get_num_rows() , &
                work, &
                A%get_num_rows(), &
                0.0, &
                subdomain_elmat, &
                this%get_block_number_coarse_dofs(1))
#else
    write (0,*) 'Error: mlbddc.f90 was not compiled with -DENABLE_BLAS.'
    write (0,*) 'Error: You must activate this cpp macro in order to use the BLAS'
    check(.false.)    
#endif
    call memfree ( work, __FILE__, __LINE__)
  end subroutine mlbddc_coarse_compute_subdomain_elmat
 
  subroutine mlbddc_coarse_compute_subdomain_elmat_counts_and_displs ( this, counts, displs )
    implicit none
    class(mlbddc_coarse_t)           , intent(inout) :: this
    integer(ip), allocatable  , intent(inout) :: counts(:)
    integer(ip), allocatable  , intent(inout) :: displs(:)
    type(par_environment_t), pointer :: par_environment
    integer(ip) :: i, l1_to_l2_size
    type(coarse_fe_iterator_t) :: iterator
    type(coarse_fe_accessor_t) :: coarse_fe
    
    par_environment => this%get_par_environment()
    assert (par_environment%am_i_l1_to_l2_root())
    l1_to_l2_size = par_environment%get_l1_to_l2_size()
    if ( allocated(counts) ) call memfree ( counts, __FILE__, __LINE__ )
    if ( allocated(displs) ) call memfree ( displs, __FILE__, __LINE__ )
    call memalloc ( l1_to_l2_size, counts, __FILE__, __LINE__ )
    call memalloc ( l1_to_l2_size, displs, __FILE__, __LINE__ )
    
    i = 1
    counts(l1_to_l2_size) = 0
    iterator = this%coarse_fe_space%create_coarse_fe_iterator()
    do while ( .not. iterator%has_finished() )
      coarse_fe = iterator%current()
      if ( coarse_fe%is_local() ) then
        counts(i) = coarse_fe%get_number_dofs()**2    
        i = i +1
      end if
      call iterator%next()
    end do
    
    displs(1) = 0
    do i=2, l1_to_l2_size
      displs(i) = displs(i-1) + counts(i-1)
    end do
  end subroutine mlbddc_coarse_compute_subdomain_elmat_counts_and_displs
    
  
  subroutine mlbddc_coarse_compute_and_gather_subdomain_elmat ( this, subdomain_elmat_gathered )
     implicit none
     class(mlbddc_coarse_t)      , intent(inout) :: this
     real(rp), allocatable, intent(inout) :: subdomain_elmat_gathered(:)
     type(par_environment_t), pointer :: par_environment
     integer(ip), allocatable :: counts(:)
     integer(ip), allocatable :: displs(:)
     real(rp), allocatable :: subdomain_elmat(:,:)
     real(rp) :: dummy_real_array_rp_2D(0,0)
     real(rp) :: dummy_real_array_rp_1D(0)
     integer(ip) :: dummy_integer_array_ip(0)
     integer(ip) :: l1_to_l2_size
     par_environment => this%get_par_environment()
     assert (par_environment%am_i_l1_to_l2_task())
     
     if ( par_environment%am_i_l1_to_l2_root() ) then
       call this%compute_subdomain_elmat_counts_and_displs(counts, displs)
       l1_to_l2_size = par_environment%get_l1_to_l2_size()
       
       if ( allocated(subdomain_elmat_gathered) ) & 
         call memfree (subdomain_elmat_gathered, __FILE__, __LINE__)
         
       call memalloc ( displs(l1_to_l2_size), & 
                       subdomain_elmat_gathered, & 
                       __FILE__, __LINE__ ) 
       
       call par_environment%l2_from_l1_gather( input_data      = dummy_real_array_rp_2D, &
                                               recv_counts     = counts, &
                                               displs          = displs, &
                                               output_data     = subdomain_elmat_gathered )
       
       call memfree ( counts, __FILE__, __LINE__ )
       call memfree ( displs, __FILE__, __LINE__ )
     else
       call this%compute_subdomain_elmat(subdomain_elmat)
       call par_environment%l2_from_l1_gather( input_data      = subdomain_elmat, &
                                               recv_counts     = dummy_integer_array_ip, &
                                               displs          = dummy_integer_array_ip, &
                                               output_data     = dummy_real_array_rp_1D )
       call memfree ( subdomain_elmat, __FILE__, __LINE__ )
     end if     
  end subroutine mlbddc_coarse_compute_and_gather_subdomain_elmat
  
  subroutine mlbddc_coarse_numerical_setup_coarse_grid_matrix ( this )
    implicit none
    class(mlbddc_coarse_t)           , intent(inout) :: this
    type(par_environment_t)   , pointer       :: L1_environment
    type(par_environment_t)   , pointer       :: L2_environment
    integer(ip)                               :: istat
    real(rp), allocatable                     :: subdomain_elmat_gathered(:)
    L1_environment => this%get_par_environment()
    if ( L1_environment%am_i_l1_to_l2_task() ) then
      call this%compute_and_gather_subdomain_elmat(subdomain_elmat_gathered)
    end if
    if ( L1_environment%am_i_lgt1_task() ) then
      L2_environment => L1_environment%get_next_level() 
      if ( L2_environment%am_i_l1_task() ) then
          call this%coarse_grid_matrix_numerical_assembly(subdomain_elmat_gathered)
       end if
    end if
    ! subdomain_elmat_gathered is only allocated on L2 MPI tasks
    if (allocated(subdomain_elmat_gathered) ) call memfree(subdomain_elmat_gathered, __FILE__, __LINE__ )
  end subroutine mlbddc_coarse_numerical_setup_coarse_grid_matrix
  
  subroutine mlbddc_coarse_coarse_grid_matrix_numerical_assembly ( this, subdomain_elmat_gathered )
    implicit none
    class(mlbddc_coarse_t)           , intent(inout) :: this
    real(rp)                  , intent(in)    :: subdomain_elmat_gathered(*)
    
    type(par_environment_t)   , pointer       :: L1_environment
    type(par_environment_t)   , pointer       :: L2_environment
    type(coarse_fe_iterator_t) :: iterator
    type(coarse_fe_accessor_t) :: coarse_fe
    type(i1p_t), allocatable :: elem2dof(:)
    logical, pointer :: field_coupling(:,:)
    integer(ip) :: ifield, jfield, i, j, istat, current
    
    L1_environment => this%get_par_environment()
    L2_environment => L1_environment%get_next_level()
    assert ( associated (L2_environment) )
    assert ( L2_environment%am_i_l1_task() )
    allocate ( elem2dof(this%coarse_fe_space%get_number_fields()), stat=istat)
    check(istat==0)     
    field_coupling => this%coarse_fe_space%get_field_coupling()
    current = 1
    iterator  = this%coarse_fe_space%create_coarse_fe_iterator()
    do while ( .not. iterator%has_finished() )
       coarse_fe = iterator%current()
       call coarse_fe%get_elem2dof(elem2dof)
       if ( coarse_fe%is_local() ) then
          do ifield=1, this%coarse_fe_space%get_number_fields()
             do jfield=1, this%coarse_fe_space%get_number_fields()
                if ((field_coupling(ifield,jfield))) then
                   do j=1, size(elem2dof(jfield)%p)
                      do i=1, size(elem2dof(ifield)%p)
                         call this%coarse_grid_matrix%insert(ia =  elem2dof(ifield)%p(i), &
                                                             ja  = elem2dof(jfield)%p(j), &
                                                             val = subdomain_elmat_gathered(current) )
                         current = current + 1
                      end do
                   end do
                end if
             end do
          end do
       end if
       call iterator%next()
    end do
    deallocate ( elem2dof, stat=istat )
    check(istat==0)     
    call this%coarse_grid_matrix%convert(csr_format)
    !call this%coarse_grid_matrix%print(6)
  end subroutine mlbddc_coarse_coarse_grid_matrix_numerical_assembly
    
  subroutine mlbddc_coarse_numerical_setup_mlbddc_coarse ( this )
   implicit none
   class(mlbddc_coarse_t)           , intent(inout) :: this
   type(par_environment_t)  , pointer :: par_environment
    
    par_environment => this%get_par_environment()
    if ( par_environment%am_i_lgt1_task() ) then
     call this%mlbddc_coarse%numerical_setup()
    end if
  end subroutine mlbddc_coarse_numerical_setup_mlbddc_coarse 

  subroutine mlbddc_coarse_numerical_setup_coarse_solver(this)
    implicit none
    class(mlbddc_coarse_t)   , intent(inout) :: this
    type(par_environment_t)  , pointer :: par_environment
    par_environment => this%get_par_environment()
    assert ( par_environment%get_l1_size() == 1 )
    assert ( par_environment%am_i_l1_task() )
    call this%coarse_solver%numerical_setup() 
    !call this%coarse_solver%log_info()
  end subroutine mlbddc_coarse_numerical_setup_coarse_solver
  
  !=============================================================================
  subroutine mlbddc_coarse_apply (op, x, y)
    implicit none
    ! Parameters
    class(mlbddc_coarse_t)   , intent(in)    :: op
    class(vector_t)   , intent(in)    :: x
    class(vector_t)   , intent(inout) :: y
    !call op%abort_if_not_in_domain(x)
    !call op%abort_if_not_in_range(y)
    call x%GuardTemp()
    select type(x)
       class is (par_scalar_array_t)
       select type(y)
          class is(par_scalar_array_t)
          call op%apply_par_scalar_array( x, y )
       end select
    end select
    call x%CleanTemp()
  end subroutine mlbddc_coarse_apply
  
  subroutine mlbddc_coarse_apply_par_scalar_array ( this, x, y )
    implicit none
    class(mlbddc_coarse_t)         , intent(in)    :: this
    type(par_scalar_array_t), intent(in)    :: x
    type(par_scalar_array_t), intent(inout) :: y
       
    type(par_scalar_array_t) :: y_I, y_G
    type(par_scalar_array_t) :: residual, residual_I, residual_G
    type(par_scalar_array_t) :: delta, delta_I, delta_G
    type(par_scalar_array_t) :: constrained_neumann_correction, &
                                constrained_neumann_correction_I, &
                                constrained_neumann_correction_G
    type(par_scalar_array_t) :: coarse_correction, &
                                coarse_correction_I, &
                                coarse_correction_G
    type(par_environment_t), pointer :: par_environment

    par_environment => this%get_par_environment()
    if ( par_environment%get_l1_size() == 1 ) then
      call this%solve_coarse_problem(x, y)
    else ! l1_size > 1 .or. l1_size < 1    
      call this%create_interior_interface_views(y, &
                                                y_I, & 
                                                y_G)
 
      ! Clone, copy, and create views for input residual
      call residual%clone(x)
      call residual%copy(x)
      call this%create_interior_interface_views(residual, &
                                                residual_I, & 
                                                residual_G)
    
      ! Clone, init, and create views for delta
      call delta%clone(x)
      call this%create_interior_interface_views(delta, &
                                                delta_I, & 
                                                delta_G)
    
      ! Clone and create views for constrained_neumann_correction
      call constrained_neumann_correction%clone(x)
      call this%create_interior_interface_views(constrained_neumann_correction, &
                                                constrained_neumann_correction_I, & 
                                                constrained_neumann_correction_G)
    
      ! Clone and create views for coarse_correction
      call coarse_correction%clone(x)
      call this%create_interior_interface_views(coarse_correction, &
                                                coarse_correction_I, & 
                                                coarse_correction_G)
    
    
      ! 1) Compute delta_I = A_II^-1 residual_I,   delta_G = 0
      call this%solve_dirichlet_problem( residual_I, delta_I )
    
      ! 2) y_I = delta_I
      call y_I%copy(delta_I)
    
      ! 3) ! I-A P_D = [ 0                    0 ]
           !           [ -A_GI A_II^{-1}      I ]
      !delta_G = A_GI * A_II^{-1}*residual_I
      call this%apply_A_GI(delta_I, delta_G)
      call delta%comm()
      call residual_I%init(0.0_rp)
      ! residual_G = residual_G - delta_G =
      !              residual_G - sum(A_GI*A_II^{-1}*residual_I,i=1..P)
      call residual_G%axpby(-1.0_rp, delta_G, 1.0_rp)

      ! 4) Compute delta_G = A_{BDDC}^{-1} r
      call this%apply_weight_operator(residual,residual) 
      call this%compute_coarse_correction(residual, coarse_correction)
      call this%solve_constrained_neumann_problem ( residual, constrained_neumann_correction )
      call delta_G%copy (coarse_correction_G)
      call delta_G%axpby(1.0_rp, constrained_neumann_correction_G,1.0_rp)
      call this%apply_weight_operator(delta,delta) 
      call delta%comm()
      ! y_G = delta_G
      call y_G%copy(delta_G)

      ! 5) I-P_D A = [ 0      -A_II^{-1} A_IG ]
      !              [ 0                   I  ]
      call this%apply_A_IG(delta_G, residual_I)
      call this%solve_dirichlet_problem(residual_I, delta_I)
      ! y_I = y_I - delta_I
      call y_I%axpby(-1.0_rp, delta_I, 1.0_rp)

      call constrained_neumann_correction_I%free()
      call constrained_neumann_correction_G%free()
      call constrained_neumann_correction%free()
      call coarse_correction_I%free()
      call coarse_correction_G%free()
      call coarse_correction%free()
      call delta_I%free()
      call delta_G%free()
      call delta%free()
      call residual_I%free()
      call residual_G%free()
      call residual%free()
      call y_I%free()
      call y_G%free()
    end if
  end subroutine mlbddc_coarse_apply_par_scalar_array
  
  subroutine mlbddc_coarse_solve_coarse_problem(this,x,y)
    implicit none
    class(mlbddc_coarse_t)            , intent(in)    :: this
    type(par_scalar_array_t)   , intent(in)    :: x
    type(par_scalar_array_t)   , intent(inout) :: y
    type(par_environment_t)    , pointer       :: par_environment
    type(serial_scalar_array_t), pointer       :: serial_y 

    par_environment => this%get_par_environment()
    assert ( par_environment%get_l1_size() == 1 )
    assert ( par_environment%am_i_l1_task() )

    serial_y => y%get_serial_scalar_array()
    call this%coarse_solver%solve(x%get_serial_scalar_array(), &
                                  serial_y)
  end subroutine mlbddc_coarse_solve_coarse_problem
  
  subroutine mlbddc_coarse_compute_coarse_correction (this, residual, coarse_correction)
    implicit none
    class(mlbddc_coarse_t)         , intent(in)    :: this
    type(par_scalar_array_t), intent(in)    :: residual
    type(par_scalar_array_t), intent(inout) :: coarse_correction
    
    type(par_scalar_array_t) :: coarse_residual
    type(par_scalar_array_t) :: coarse_coarse_correction
    type(par_environment_t), pointer :: L1_environment
    type(par_environment_t), pointer :: L2_environment
    
    L1_environment => this%get_par_environment()
    if ( L1_environment%am_i_lgt1_task() ) then
       L2_environment => L1_environment%get_next_level()
       call coarse_residual%create_and_allocate( p_env      = L2_environment, &
                                                 dof_import = this%coarse_fe_space%get_block_dof_import(1) )
       call coarse_coarse_correction%clone(coarse_residual)
    end if   
    
    ! Transfer residual from L1 to L2 tasks
    call this%setup_coarse_grid_residual ( residual, coarse_residual )
    
    ! Solve coarse problem on > L1 tasks
    if ( L1_environment%am_i_lgt1_task() ) then
       ! call this%mlbddc_coarse_coarse%apply(coarse_residual, coarse_coarse_correction)
       call coarse_coarse_correction%copy(coarse_residual)
    end if
    call this%scatter_and_interpolate_coarse_grid_correction ( coarse_coarse_correction, coarse_correction )
    if ( L1_environment%am_i_lgt1_task() ) then
       call coarse_coarse_correction%free()
       call coarse_residual%free()
    end if
    
  end subroutine mlbddc_coarse_compute_coarse_correction
  
  subroutine mlbddc_coarse_setup_coarse_grid_residual ( this, vector, coarse_grid_vector )
    implicit none
    class(mlbddc_coarse_t)         , intent(in)    :: this
    type(par_scalar_array_t), intent(in)    :: vector
    type(par_scalar_array_t), intent(inout) :: coarse_grid_vector
    type(par_environment_t), pointer        :: L1_environment
    type(par_environment_t), pointer        :: L2_environment
    real(rp), allocatable                   :: coarse_dofs_values_gathered(:)
    
    L1_environment => this%get_par_environment()
    if ( L1_environment%am_i_l1_to_l2_task() ) then
      call this%compute_and_gather_coarse_dofs_values(vector, coarse_dofs_values_gathered)
    end if
    
    if ( L1_environment%am_i_lgt1_task() ) then
       call coarse_grid_vector%init(0.0_rp) 
       L2_environment => L1_environment%get_next_level()
       if ( L2_environment%am_i_l1_task() ) then
          call this%coarse_grid_residual_assembly(coarse_dofs_values_gathered, coarse_grid_vector)
       end if
    else
       call coarse_grid_vector%free()
    end if
    if ( allocated(coarse_dofs_values_gathered) ) call memfree(coarse_dofs_values_gathered, __FILE__, __LINE__ )
  end subroutine mlbddc_coarse_setup_coarse_grid_residual
  
  
  ! Computes subdomain_elvec = \Phi^t v_i
  subroutine mlbddc_coarse_compute_coarse_dofs_values ( this, vector, coarse_dofs_values )
    implicit none
    class(mlbddc_coarse_t)         , intent(in) :: this
    type(par_scalar_array_t), intent(in) :: vector
    real(rp), allocatable, intent(inout) :: coarse_dofs_values(:)
    type(serial_scalar_array_t), pointer :: serial_scalar_array
    
    assert ( this%am_i_l1_task() )
    serial_scalar_array => vector%get_serial_scalar_array()
    
    if ( allocated(coarse_dofs_values) ) then
      call memfree ( coarse_dofs_values, __FILE__, __LINE__ )
    end if
    
    call memalloc ( this%get_block_number_coarse_dofs(1), &
                    coarse_dofs_values, &
                    __FILE__, &
                    __LINE__ );
    
    coarse_dofs_values = 0.0_rp
#ifdef ENABLE_BLAS
    call DGEMV(  'T', & 
                 serial_scalar_array%get_size(), &
                 this%get_block_number_coarse_dofs(1), &
                 1.0_rp, &
                 this%Phi, &
                 serial_scalar_array%get_size(), &
                 serial_scalar_array%get_entries(), &
                 1, &
                 0.0_rp, & 
                 coarse_dofs_values, & 
                 1)
#else
    write (0,*) 'Error: mlbddc_coarse.f90 was not compiled with -DENABLE_BLAS.'
    write (0,*) 'Error: You must activate this cpp macro in order to use the BLAS'
    check(.false.)    
#endif
  end subroutine mlbddc_coarse_compute_coarse_dofs_values
 
  subroutine mlbddc_coarse_compute_coarse_dofs_values_counts_and_displs ( this, counts, displs )
    implicit none
    class(mlbddc_coarse_t)           , intent(in)    :: this
    integer(ip), allocatable  , intent(inout) :: counts(:)
    integer(ip), allocatable  , intent(inout) :: displs(:)
    type(par_environment_t), pointer :: par_environment
    integer(ip) :: i, l1_to_l2_size
    type(coarse_fe_iterator_t) :: iterator
    type(coarse_fe_accessor_t) :: coarse_fe
    
    par_environment => this%get_par_environment()
    assert (par_environment%am_i_l1_to_l2_root())
    l1_to_l2_size = par_environment%get_l1_to_l2_size()
    if ( allocated(counts) ) call memfree ( counts, __FILE__, __LINE__ )
    if ( allocated(displs) ) call memfree ( displs, __FILE__, __LINE__ )
    call memalloc ( l1_to_l2_size, counts, __FILE__, __LINE__ )
    call memalloc ( l1_to_l2_size, displs, __FILE__, __LINE__ )
    
    i=1
    counts(l1_to_l2_size) = 0
    iterator = this%coarse_fe_space%create_coarse_fe_iterator()
    do while ( .not. iterator%has_finished() )
      coarse_fe = iterator%current()
      if ( coarse_fe%is_local() ) then
        counts(i) = coarse_fe%get_number_dofs()
        i=i+1
      end if
      call iterator%next()
    end do
    
    displs(1) = 0
    do i=2, l1_to_l2_size
      displs(i) = displs(i-1) + counts(i-1)
    end do
  end subroutine mlbddc_coarse_compute_coarse_dofs_values_counts_and_displs
  
  subroutine mlbddc_coarse_compute_and_gather_coarse_dofs_values ( this, vector, coarse_dofs_values_gathered )
     implicit none
     class(mlbddc_coarse_t)         , intent(in)    :: this
     type(par_scalar_array_t), intent(in)    :: vector
     real(rp)   , allocatable, intent(inout) :: coarse_dofs_values_gathered(:)
     type(par_environment_t), pointer :: par_environment
     integer(ip), allocatable :: counts(:)
     integer(ip), allocatable :: displs(:)
     real(rp), allocatable :: coarse_dofs_values(:)
     real(rp) :: dummy_real_array_rp(0)
     integer(ip) :: dummy_integer_array_ip(0)
     integer(ip) :: l1_to_l2_size
     par_environment => this%get_par_environment()
     assert (par_environment%am_i_l1_to_l2_task())
     
     if ( par_environment%am_i_l1_to_l2_root() ) then
       call this%compute_coarse_dofs_values_counts_and_displs(counts, displs)
       l1_to_l2_size = par_environment%get_l1_to_l2_size()
       
       if ( allocated(coarse_dofs_values_gathered) ) & 
         call memfree (coarse_dofs_values_gathered, __FILE__, __LINE__)
         
       call memalloc ( displs(l1_to_l2_size), & 
                       coarse_dofs_values_gathered, & 
                       __FILE__, __LINE__ ) 
       
       call par_environment%l2_from_l1_gather( input_data_size = 0, &
                                               input_data      = dummy_real_array_rp, &
                                               recv_counts     = counts, &
                                               displs          = displs, &
                                               output_data     = coarse_dofs_values_gathered )
       
       call memfree ( counts, __FILE__, __LINE__ )
       call memfree ( displs, __FILE__, __LINE__ )
     else
       call this%compute_coarse_dofs_values(vector, coarse_dofs_values)
       call par_environment%l2_from_l1_gather( input_data_size = size(coarse_dofs_values), &
                                               input_data      = coarse_dofs_values, &
                                               recv_counts     = dummy_integer_array_ip, &
                                               displs          = dummy_integer_array_ip, &
                                               output_data     = dummy_real_array_rp )
       call memfree ( coarse_dofs_values, __FILE__, __LINE__ )
     end if     
  end subroutine mlbddc_coarse_compute_and_gather_coarse_dofs_values
    
  ! This subroutine assumes that coarse_grid_vector has been already created+allocated, and initialized to zero.
  subroutine mlbddc_coarse_coarse_grid_residual_assembly ( this, coarse_dofs_values_gathered, coarse_grid_vector ) 
    implicit none
    class(mlbddc_coarse_t)         , intent(in)    :: this
    real(rp)                , intent(in)    :: coarse_dofs_values_gathered(*)
    type(par_scalar_array_t), intent(inout) :: coarse_grid_vector
     
    type(par_environment_t)   , pointer       :: L1_environment
    type(par_environment_t)   , pointer       :: L2_environment
    type(coarse_fe_iterator_t)                :: iterator
    type(coarse_fe_accessor_t)                :: coarse_fe
    
    type(i1p_t), allocatable :: elem2dof(:)
    integer(ip) :: ifield, i, istat, current
    
    L1_environment => this%get_par_environment()
    L2_environment => L1_environment%get_next_level()
    assert ( associated (L2_environment) )
    assert ( L2_environment%am_i_l1_task() )

    allocate ( elem2dof(this%coarse_fe_space%get_number_fields()), stat=istat)
    check(istat==0)     
    current = 1
    iterator  = this%coarse_fe_space%create_coarse_fe_iterator()
    do while ( .not. iterator%has_finished() )
       coarse_fe = iterator%current()
       call coarse_fe%get_elem2dof(elem2dof)
       if ( coarse_fe%is_local() ) then
          do ifield=1, this%coarse_fe_space%get_number_fields()
            do i=1, size(elem2dof(ifield)%p)
              call coarse_grid_vector%add(i   =  elem2dof(ifield)%p(i), &
                                          val = coarse_dofs_values_gathered(current) )
              current = current + 1
            end do
          end do
       end if
       call iterator%next()
    end do
    deallocate ( elem2dof, stat=istat )
    check(istat==0)
    call coarse_grid_vector%comm()
  end subroutine mlbddc_coarse_coarse_grid_residual_assembly
  
  subroutine mlbddc_coarse_scatter_and_interpolate_coarse_grid_correction (this, coarse_grid_vector, vector)
    implicit none
    class(mlbddc_coarse_t)           , intent(in)    :: this
    type(par_scalar_array_t)  , intent(in)    :: coarse_grid_vector
    type(par_scalar_array_t)  , intent(inout) :: vector
    type(par_environment_t)   , pointer :: L1_environment
    real(rp), allocatable :: coarse_dofs_values(:)
    
    L1_environment => this%get_par_environment()
    if ( L1_environment%am_i_l1_to_l2_task() ) then
      call this%scatter_coarse_grid_correction(coarse_grid_vector, coarse_dofs_values)
    end if
    
    if ( L1_environment%am_i_l1_task() ) then
      call this%interpolate_coarse_grid_correction(coarse_dofs_values, vector)
    end if
    
    if ( allocated(coarse_dofs_values) ) call memfree(coarse_dofs_values, __FILE__, __LINE__ )
  end subroutine mlbddc_coarse_scatter_and_interpolate_coarse_grid_correction
  
  subroutine mlbddc_coarse_scatter_coarse_grid_correction ( this, coarse_grid_vector, coarse_dofs_values )
    implicit none
    class(mlbddc_coarse_t)         , intent(in)    :: this
    type(par_scalar_array_t), intent(in)    :: coarse_grid_vector
    real(rp), allocatable   , intent(inout) :: coarse_dofs_values(:)
    
    integer(ip) :: l1_to_l2_size
    real(rp), allocatable :: coarse_dofs_values_scattered(:)
    integer(ip), allocatable :: counts(:), displs(:)
    type(par_environment_t), pointer :: par_environment
    real(rp) :: dummy_real_array_rp(0)
    integer(ip) :: dummy_integer_array_ip(0)
    
    par_environment => this%get_par_environment()
    assert (par_environment%am_i_l1_to_l2_task())
    if ( par_environment%am_i_l1_to_l2_root() ) then
       
       call this%compute_coarse_dofs_values_counts_and_displs(counts, displs)
       l1_to_l2_size = par_environment%get_l1_to_l2_size()
       
       call memalloc ( displs(l1_to_l2_size), & 
                       coarse_dofs_values_scattered, & 
                       __FILE__, __LINE__ ) 
       
       call this%fill_coarse_dofs_values_scattered ( coarse_grid_vector, &
                                                     coarse_dofs_values_scattered)
       
       call par_environment%l2_to_l1_scatter( input_data       = coarse_dofs_values_scattered, &
                                              send_counts      = counts, &
                                              displs           = displs, &
                                              output_data_size = 0, &
                                              output_data      = dummy_real_array_rp )
       
       call memfree ( coarse_dofs_values_scattered, & 
                       __FILE__, __LINE__ ) 
       
       call memfree ( counts, __FILE__, __LINE__ )
       call memfree ( displs, __FILE__, __LINE__ )
       
    else 
       if ( allocated (coarse_dofs_values) ) call memfree( coarse_dofs_values, __FILE__, __LINE__ )
       call memalloc ( this%get_block_number_coarse_dofs(1), &
                       coarse_dofs_values, &
                       __FILE__, __LINE__ );
       
       call par_environment%l2_to_l1_scatter( input_data       = dummy_real_array_rp, &
                                              send_counts      = dummy_integer_array_ip, &
                                              displs           = dummy_integer_array_ip, &
                                              output_data_size = size(coarse_dofs_values), &
                                              output_data      = coarse_dofs_values )
    end if
  end subroutine mlbddc_coarse_scatter_coarse_grid_correction
  
  subroutine mlbddc_coarse_fill_coarse_dofs_values_scattered ( this, coarse_grid_vector, coarse_dofs_values_scattered ) 
    implicit none
    class(mlbddc_coarse_t)          , intent(in)     :: this
    type(par_scalar_array_t) , intent(in)     :: coarse_grid_vector
    real(rp)                 , intent(inout)  :: coarse_dofs_values_scattered(*)
    type(par_environment_t), pointer    :: par_environment
    type(i1p_t)           , allocatable :: elem2dof(:)
    type(coarse_fe_iterator_t)          :: iterator
    type(coarse_fe_accessor_t)          :: coarse_fe
    integer(ip)                         :: istat, current, ifield
    
    par_environment => this%get_par_environment()
    assert (par_environment%am_i_l1_to_l2_root())
    
    allocate ( elem2dof(this%coarse_fe_space%get_number_fields()), stat=istat)
    check(istat==0)
    
    current = 1
    iterator  = this%coarse_fe_space%create_coarse_fe_iterator()
    do while ( .not. iterator%has_finished() )
       coarse_fe = iterator%current()
       call coarse_fe%get_elem2dof(elem2dof)
       if ( coarse_fe%is_local() ) then
          do ifield=1, this%coarse_fe_space%get_number_fields()
              call coarse_grid_vector%extract_subvector ( iblock       = 1, &
                                                          size_indices = size(elem2dof(ifield)%p), &
                                                          indices     = elem2dof(ifield)%p, &
                                                          values      = coarse_dofs_values_scattered(current) )
              current = current + size(elem2dof(ifield)%p)
          end do
       end if
       call iterator%next()
    end do
    
    deallocate ( elem2dof, stat=istat )
    check(istat==0)
  end subroutine mlbddc_coarse_fill_coarse_dofs_values_scattered 
    
  subroutine mlbddc_coarse_interpolate_coarse_grid_correction (this, coarse_dofs_values, vector)
    implicit none
    class(mlbddc_coarse_t)            , intent(in)    :: this
    real(rp)                   , intent(in)    :: coarse_dofs_values(*)
    type(par_scalar_array_t)   , intent(inout) :: vector
    type(par_environment_t)    , pointer       :: L1_environment
    type(serial_scalar_array_t), pointer       :: serial_scalar_array 
    real(rp)                   , pointer       :: serial_scalar_array_entries(:)
    
    L1_environment => this%get_par_environment()
    assert ( L1_environment%am_i_l1_task() )
    
    call vector%init(0.0_rp)
    serial_scalar_array         => vector%get_serial_scalar_array()
    serial_scalar_array_entries => serial_scalar_array%get_entries()
    
#ifdef ENABLE_BLAS
    call DGEMV(  'N', & 
                  serial_scalar_array%get_size(), &
                  this%get_block_number_coarse_dofs(1), &
                  1.0_rp, &
                  this%Phi, &
                  serial_scalar_array%get_size(), &
                  coarse_dofs_values, &
                  1,    &
                  0.0_rp,  & 
                  serial_scalar_array_entries, & 
                  1)
#else
     write (0,*) 'Error: mlbddc_coarse.f90 was not compiled with -DENABLE_BLAS.'
     write (0,*) 'Error: You must activate this cpp macro in order to use the BLAS'
     check(.false.)    
#endif
     
  end subroutine mlbddc_coarse_interpolate_coarse_grid_correction
  
  
  subroutine mlbddc_coarse_solve_dirichlet_problem(this, x_I, y_I)
    implicit none
    class(mlbddc_coarse_t)           , intent(in)    :: this
    type(par_scalar_array_t)  , intent(in)    :: x_I
    type(par_scalar_array_t)  , intent(inout) :: y_I
    type(serial_scalar_array_t), pointer      :: y_I_serial
    if ( this%am_i_l1_task() ) then
       y_I_serial => y_I%get_serial_scalar_array()
       call this%dirichlet_solver%solve(x_I%get_serial_scalar_array(), &
                                        y_I_serial)
    end if   
  end subroutine mlbddc_coarse_solve_dirichlet_problem  
  
  subroutine mlbddc_coarse_apply_A_GI(this, x_I, y_G)
    implicit none
    class(mlbddc_coarse_t)           , intent(in)    :: this
    type(par_scalar_array_t)  , intent(in)    :: x_I
    type(par_scalar_array_t)  , intent(inout) :: y_G
    type(par_sparse_matrix_t), pointer :: par_sparse_matrix
    type(sparse_matrix_t), pointer :: A
    type(serial_scalar_array_t), pointer    :: y_G_serial
    
    if ( this%am_i_l1_task() ) then
      par_sparse_matrix => this%get_par_sparse_matrix()
      A => par_sparse_matrix%get_sparse_matrix()
      y_G_serial => y_G%get_serial_scalar_array()
      if ( A%get_symmetric_storage() ) then
         call this%A_IG%apply_transpose(x_I%get_serial_scalar_array(), &
                                        y_G_serial)
      else
         call this%A_GI%apply(x_I%get_serial_scalar_array(), &
                              y_G_serial)
      end if
    end if
    
  end subroutine mlbddc_coarse_apply_A_GI
  
  subroutine mlbddc_coarse_apply_A_IG(this, x_G, y_I)
    implicit none
    class(mlbddc_coarse_t)           , intent(in)    :: this
    type(par_scalar_array_t)  , intent(in)    :: x_G
    type(par_scalar_array_t)  , intent(inout) :: y_I
    type(par_sparse_matrix_t), pointer :: par_sparse_matrix
    type(sparse_matrix_t), pointer :: A
    type(serial_scalar_array_t), pointer    :: y_I_serial
    if ( this%am_i_l1_task() ) then
      par_sparse_matrix => this%get_par_sparse_matrix()
      A => par_sparse_matrix%get_sparse_matrix()
      y_I_serial => y_I%get_serial_scalar_array()
      if ( A%get_symmetric_storage() ) then
         call this%A_IG%apply(x_G%get_serial_scalar_array(), &
                              y_I_serial)
      end if
    end if 
  end subroutine mlbddc_coarse_apply_A_IG
  
  subroutine mlbddc_coarse_solve_constrained_neumann_problem(this, x, y)
    implicit none
    class(mlbddc_coarse_t)            , intent(in)    :: this
    type(par_scalar_array_t)   , intent(in)    :: x
    type(par_scalar_array_t)   , intent(inout) :: y
    type(coarse_fe_space_t)    , pointer :: coarse_fe_space
    type(serial_scalar_array_t), pointer :: x_serial
    type(serial_scalar_array_t), pointer :: y_serial
    real(rp)                   , pointer :: y_serial_entries(:)
    type(serial_scalar_array_t)          :: augmented_x
    type(serial_scalar_array_t)          :: augmented_y
    real(rp), allocatable                :: augmented_x_entries(:)
    real(rp), allocatable                :: augmented_y_entries(:)
    integer(ip)                          :: block_number_dofs
    integer(ip)                          :: block_number_coarse_dofs

    if ( this%am_i_l1_task() ) then
      coarse_fe_space => this%get_fe_space()
      block_number_dofs        = coarse_fe_space%get_block_number_dofs(1)
      block_number_coarse_dofs = this%get_block_number_coarse_dofs(1)
      
      x_serial => x%get_serial_scalar_array()
    
      ! Set-up augmented_x from x_serial
      call augmented_x%create(block_number_dofs + block_number_coarse_dofs)
      call memalloc ( block_number_dofs + block_number_coarse_dofs, & 
                      augmented_x_entries, __FILE__,__LINE__)  
      augmented_x_entries(1:block_number_dofs)  = x_serial%get_entries() 
      augmented_x_entries(block_number_dofs+1:) = 0.0_rp 
      call augmented_x%set_view_entries(augmented_x_entries)

      ! Set-up augmented_y
      call augmented_y%create(block_number_dofs + block_number_coarse_dofs)
      call memalloc ( block_number_dofs + block_number_coarse_dofs, & 
                      augmented_y_entries, __FILE__,__LINE__)  
      call augmented_y%set_view_entries(augmented_y_entries)
      call this%constrained_neumann_solver%solve(augmented_x, &
                                                 augmented_y)

      ! Set-up y from augmented_y
      y_serial         => y%get_serial_scalar_array()
      y_serial_entries => y_serial%get_entries()
      y_serial_entries = augmented_y_entries(1:block_number_dofs)
      
      call memfree ( augmented_x_entries, __FILE__,__LINE__)  
      call memfree ( augmented_y_entries, __FILE__,__LINE__)  
    end if  
  end subroutine mlbddc_coarse_solve_constrained_neumann_problem
       
  subroutine mlbddc_coarse_apply_weight_operator(this, x, y)
    implicit none
    class(mlbddc_coarse_t)    , intent(in)    :: this
    type(par_scalar_array_t)  , intent(in)    :: x
    type(par_scalar_array_t)  , intent(inout) :: y
    type(serial_scalar_array_t), pointer :: x_local
    type(serial_scalar_array_t), pointer :: y_local
    real(rp), pointer :: x_local_entries(:)
    real(rp), pointer :: y_local_entries(:)
    type(coarse_dof_object_iterator_t) :: coarse_dofs_object_iterator
    type(coarse_dof_object_accessor_t) :: coarse_dof_object
    type(list_iterator_t) :: dofs_on_object
 
    if ( this%am_i_l1_task() ) then
      x_local         => x%get_serial_scalar_array()
      x_local_entries => x_local%get_entries()
      y_local         => y%get_serial_scalar_array()
      y_local_entries => y_local%get_entries()
      coarse_dofs_object_iterator = this%create_field_dofs_object_iterator(1)
      do while ( .not. coarse_dofs_object_iterator%has_finished() ) 
         coarse_dof_object = coarse_dofs_object_iterator%current()
         dofs_on_object = coarse_dof_object%get_dofs_on_object_iterator()
         do while ( .not. dofs_on_object%is_upper_bound() )
           y_local_entries(dofs_on_object%get_current()) = &
             x_local_entries(dofs_on_object%get_current())/coarse_dof_object%get_number_parts_around()
           call dofs_on_object%next()
         end do
         call coarse_dofs_object_iterator%next()
      end do
    end if
  end subroutine mlbddc_coarse_apply_weight_operator

  subroutine mlbddc_coarse_create_interior_interface_views ( this, x, x_I, X_G )
    implicit none
    class(mlbddc_coarse_t)         , intent(in)       :: this
    type(par_scalar_array_t), intent(in)       :: x
    type(par_scalar_array_t), intent(inout)    :: x_I
    type(par_scalar_array_t), intent(inout)    :: x_G
    type(coarse_fe_space_t), pointer :: coarse_fe_space
    coarse_fe_space => this%get_fe_space()
    call x%create_view(1, &
                       coarse_fe_space%get_block_number_interior_dofs(1), &
                       x_I)
    call x%create_view(coarse_fe_space%get_block_number_interior_dofs(1)+1, &
                       coarse_fe_space%get_block_number_dofs(1), &
                       x_G)
  end subroutine mlbddc_coarse_create_interior_interface_views 
  
  
  subroutine mlbddc_coarse_free(this)
    implicit none
    class(mlbddc_coarse_t)           , intent(inout) :: this
    call this%free_numerical_setup()
    call this%free_symbolic_setup()
    call this%free_clean()
  end subroutine mlbddc_coarse_free
  
  subroutine mlbddc_coarse_free_clean(this)
    implicit none
    class(mlbddc_coarse_t)           , intent(inout) :: this
    nullify(this%coarse_fe_space)
    nullify(this%par_sparse_matrix)
  end subroutine mlbddc_coarse_free_clean
  
  subroutine mlbddc_coarse_free_symbolic_setup(this)
    implicit none
    class(mlbddc_coarse_t)           , intent(inout) :: this
    type(par_environment_t)   , pointer       :: par_environment
    integer(ip)                               :: istat
    
    ! This will be replaced by an apropriate check on the state diagram
    if ( associated(this%fe_space) ) then 
      par_environment => this%get_par_environment()
      if ( par_environment%am_i_l1_task() ) then
        if ( par_environment%get_l1_size() > 1 ) then  
          call this%free_dofs_objects_and_constraint_matrix()
          call this%free_symbolic_setup_dirichlet_solver()
          call this%free_symbolic_setup_dirichlet_problem()
          call this%free_symbolic_setup_constrained_neumann_solver()
          call this%free_symbolic_setup_constrained_neumann_problem()
          nullify (this%mlbddc_coarse)
          nullify (this%coarse_grid_matrix)
          nullify (this%coarse_fe_space)
        else
          call this%free_symbolic_setup_coarse_solver()
        end if    
      else
          ! mlbddc_coarse should be freed before coarse_fe_space
          ! (as the former was created from the latter)
          call this%mlbddc_coarse%free_symbolic_setup()
          call this%mlbddc_coarse%free_clean()
          deallocate (this%mlbddc_coarse, stat=istat)
          check (istat==0)
      
          call this%coarse_grid_matrix%free()
          deallocate  ( this%coarse_grid_matrix, stat = istat )
          check( istat == 0 )
        
          call this%coarse_fe_space%free()
          deallocate (this%coarse_fe_space, stat=istat)
          check (istat==0)
      end if
    end if   
  end subroutine mlbddc_coarse_free_symbolic_setup 
  
  subroutine mlbddc_coarse_free_dofs_objects_and_constraint_matrix(this)
    implicit none
    class(mlbddc_coarse_t)           , intent(inout) :: this
    integer(ip)                               :: i, istat 
    
    assert ( this%am_i_l1_task() )
    
    if (allocated(this%num_dofs_objects_per_field)) then 
      call memfree ( this%num_dofs_objects_per_field, __FILE__, __LINE__ )
    end if
  
    if (allocated(this%dofs_objects_per_field)) then
      do i=1, size(this%dofs_objects_per_field)
        call this%dofs_objects_per_field(i)%free()
      end do  
      deallocate ( this%dofs_objects_per_field, stat=istat )
      check (istat == 0)
    end if
  
    if (allocated(this%dofs_objects_gids_per_field)) then
      do i=1, size(this%dofs_objects_gids_per_field)
        call this%dofs_objects_gids_per_field(i)%free()
      end do  
      deallocate ( this%dofs_objects_gids_per_field, stat=istat )
      check (istat == 0)  
    end if
  
    if (allocated(this%vefs_lids_dofs_objects_per_field)) then
      do i=1, size(this%vefs_lids_dofs_objects_per_field)
        call this%vefs_lids_dofs_objects_per_field(i)%free()
      end do  
      deallocate ( this%vefs_lids_dofs_objects_per_field, stat=istat )
      check (istat == 0)  
    end if
    
    call this%constraint_matrix%free()
  end subroutine mlbddc_coarse_free_dofs_objects_and_constraint_matrix
  

  subroutine mlbddc_coarse_free_symbolic_setup_dirichlet_problem(this)
    implicit none
    class(mlbddc_coarse_t)           , intent(inout) :: this
    assert ( this%am_i_l1_task() )
    call this%A_II%free()
    call this%A_IG%free()
    call this%A_GI%free()
    call this%A_GG%free()
  end subroutine mlbddc_coarse_free_symbolic_setup_dirichlet_problem
  
  subroutine mlbddc_coarse_free_symbolic_setup_dirichlet_solver(this)
    implicit none
    class(mlbddc_coarse_t)           , intent(inout) :: this
    assert ( this%am_i_l1_task() )
    call this%dirichlet_solver%free()
  end subroutine mlbddc_coarse_free_symbolic_setup_dirichlet_solver
  
  subroutine mlbddc_coarse_free_symbolic_setup_coarse_solver(this)
    implicit none
    class(mlbddc_coarse_t)           , intent(inout) :: this
    assert ( this%am_i_l1_task() )
    call this%coarse_solver%free()
  end subroutine mlbddc_coarse_free_symbolic_setup_coarse_solver
  
  subroutine mlbddc_coarse_free_symbolic_setup_constrained_neumann_problem(this)
    implicit none
    class(mlbddc_coarse_t), intent(inout) :: this
    assert ( this%am_i_l1_task() )
    call this%constrained_neumann_matrix%free()
  end subroutine mlbddc_coarse_free_symbolic_setup_constrained_neumann_problem
  
  subroutine mlbddc_coarse_free_symbolic_setup_constrained_neumann_solver(this)
    implicit none
    class(mlbddc_coarse_t)           , intent(inout) :: this
    assert ( this%am_i_l1_task() )
    call this%constrained_neumann_solver%free()
  end subroutine mlbddc_coarse_free_symbolic_setup_constrained_neumann_solver
  
  subroutine mlbddc_coarse_free_numerical_setup(this)
    implicit none
    class(mlbddc_coarse_t)           , intent(inout) :: this
    type(par_environment_t)   , pointer       :: par_environment
    ! This will be replaced by an apropriate check on the state diagram
    
    if ( associated(this%fe_space) ) then 
      par_environment => this%get_par_environment()
      if ( par_environment%am_i_l1_task() ) then
        if ( par_environment%get_l1_size() > 1 ) then
          call this%free_numerical_setup_dirichlet_solver()
          call this%free_numerical_setup_dirichlet_problem()
          call this%free_numerical_setup_constrained_neumann_solver()
          call this%free_numerical_setup_constrained_neumann_problem()
          call this%free_coarse_grid_basis()
        else if ( par_environment%get_l1_size() == 1 ) then
          call this%free_numerical_setup_coarse_solver()
        end if
      else
        call this%mlbddc_coarse%free_numerical_setup()
        call this%coarse_grid_matrix%free_in_stages(free_numerical_setup)
      end if
    end if
  end subroutine mlbddc_coarse_free_numerical_setup
  
  subroutine mlbddc_coarse_free_numerical_setup_dirichlet_problem(this)
    implicit none
    class(mlbddc_coarse_t)           , intent(inout) :: this
    assert ( this%am_i_l1_task() )
    call this%A_II%free_in_stages(free_numerical_setup)
    call this%A_IG%free_in_stages(free_numerical_setup)
    call this%A_GI%free_in_stages(free_numerical_setup)
    call this%A_GG%free_in_stages(free_numerical_setup)
  end subroutine mlbddc_coarse_free_numerical_setup_dirichlet_problem
  
  subroutine mlbddc_coarse_free_numerical_setup_dirichlet_solver(this)
    implicit none
    class(mlbddc_coarse_t)           , intent(inout) :: this
    assert ( this%am_i_l1_task() )
    call this%dirichlet_solver%free_in_stages(free_numerical_setup)
  end subroutine mlbddc_coarse_free_numerical_setup_dirichlet_solver
  
  subroutine mlbddc_coarse_free_numerical_setup_constrained_neumann_problem(this)
    implicit none
    class(mlbddc_coarse_t)           , intent(inout) :: this
    assert ( this%am_i_l1_task() )
    call this%constrained_neumann_matrix%free_in_stages(free_numerical_setup)
  end subroutine mlbddc_coarse_free_numerical_setup_constrained_neumann_problem
  
  subroutine mlbddc_coarse_free_numerical_setup_constrained_neumann_solver(this)
    implicit none
    class(mlbddc_coarse_t)           , intent(inout) :: this
    assert ( this%am_i_l1_task() )
    call this%constrained_neumann_solver%free_in_stages(free_numerical_setup)
  end subroutine mlbddc_coarse_free_numerical_setup_constrained_neumann_solver
  
  subroutine mlbddc_coarse_free_coarse_grid_basis ( this ) 
    implicit none
    class(mlbddc_coarse_t)           , intent(inout) :: this
    assert ( this%am_i_l1_task() )
    if ( allocated ( this%Phi ) ) then
      call memfree ( this%Phi, __FILE__, __LINE__ )
    end if
  end subroutine mlbddc_coarse_free_coarse_grid_basis
  
  subroutine mlbddc_coarse_free_numerical_setup_coarse_solver ( this )
    implicit none
    class(mlbddc_coarse_t), intent(inout) :: this
    type(par_environment_t), pointer :: par_environment
    par_environment => this%get_par_environment()
    assert ( par_environment%get_l1_size() == 1 )
    assert ( par_environment%am_i_l1_task() )
    call this%coarse_solver%free_in_stages(free_numerical_setup)
  end subroutine mlbddc_coarse_free_numerical_setup_coarse_solver
  
  
  function mlbddc_coarse_get_total_number_coarse_dofs ( this )
    implicit none
    class(mlbddc_coarse_t)           , intent(in) :: this
    integer(ip)                                   :: mlbddc_coarse_get_total_number_coarse_dofs
    type(coarse_fe_space_t)     , pointer         :: fe_space
    integer(ip)                                   :: field_id
    
    assert ( this%am_i_l1_task() )
    
    mlbddc_coarse_get_total_number_coarse_dofs = 0
    fe_space => this%get_fe_space()
    do field_id = 1, fe_space%get_number_fields()
       mlbddc_coarse_get_total_number_coarse_dofs = mlbddc_coarse_get_total_number_coarse_dofs + & 
                                             this%num_dofs_objects_per_field(field_id)
    end do
  end function mlbddc_coarse_get_total_number_coarse_dofs
  
  function mlbddc_coarse_get_block_number_coarse_dofs ( this, block_id )
    implicit none
    class(mlbddc_coarse_t)           , intent(in)    :: this
    integer(ip)               , intent(in)    :: block_id
    integer(ip)                               :: mlbddc_coarse_get_block_number_coarse_dofs 
    type(coarse_fe_space_t)      , pointer    :: fe_space
    integer(ip)                               :: field_id
    integer(ip)               , pointer       :: field_blocks(:)
    assert ( this%am_i_l1_task() )
    assert ( block_id == 1 )
    
    fe_space     => this%get_fe_space()
    field_blocks => fe_space%get_field_blocks()
    
    mlbddc_coarse_get_block_number_coarse_dofs = 0
    do field_id = 1, fe_space%get_number_fields()
       if ( field_blocks(field_id) == block_id ) then
         mlbddc_coarse_get_block_number_coarse_dofs = mlbddc_coarse_get_block_number_coarse_dofs + & 
                                               this%num_dofs_objects_per_field(field_id)
       end if                                      
    end do
  end function mlbddc_coarse_get_block_number_coarse_dofs
  
  function mlbddc_coarse_create_field_dofs_object_iterator(this, field_id)
    implicit none
    class(mlbddc_coarse_t)   , intent(in) :: this
    integer(ip)              , intent(in) :: field_id
    type(coarse_dof_object_iterator_t) :: mlbddc_coarse_create_field_dofs_object_iterator
    assert ( field_id >=1 .and. field_id <= this%fe_space%get_number_fields() )
    call mlbddc_coarse_create_field_dofs_object_iterator%create(1, field_id, this)
  end function mlbddc_coarse_create_field_dofs_object_iterator
  
  ! Helper function that extracts a run-time polymorphic class(matrix_t)
  ! from XXX, and dynamically casts it into  
  ! type(par_sparse_matrix_t). If the dynamic cast cannot be performed 
  ! [because class(matrix_t) is NOT of type(par_sparse_matrix_t)], then it 
  ! aborts the execution of the program.
  function mlbddc_coarse_get_par_sparse_matrix(this)
    implicit none
    class(mlbddc_coarse_t)   , intent(in) :: this
    type(par_sparse_matrix_t), pointer    :: mlbddc_coarse_get_par_sparse_matrix
    mlbddc_coarse_get_par_sparse_matrix => this%par_sparse_matrix
  end function mlbddc_coarse_get_par_sparse_matrix
  
  ! Helper function that extracts type(coarse_fe_space_t) from XXX
  function mlbddc_coarse_get_fe_space(this)
    implicit none
    class(mlbddc_coarse_t)   , intent(in) :: this
    type(coarse_fe_space_t)  , pointer    :: mlbddc_coarse_get_fe_space
    mlbddc_coarse_get_fe_space => this%fe_space
  end function mlbddc_coarse_get_fe_space
  
  function mlbddc_coarse_get_par_environment(this)
    implicit none
    class(mlbddc_coarse_t)   , intent(in) :: this
    type(par_environment_t)  , pointer    :: mlbddc_coarse_get_par_environment
    type(coarse_fe_space_t)  , pointer    :: fe_space
    fe_space => this%get_fe_space()
    mlbddc_coarse_get_par_environment => fe_space%get_par_environment()
  end function mlbddc_coarse_get_par_environment
  
  function mlbddc_coarse_am_i_l1_task(this)
    implicit none
    class(mlbddc_coarse_t)   , intent(in) :: this
    logical                               :: mlbddc_coarse_am_i_l1_task
    type(par_environment_t)   , pointer   :: par_environment
    par_environment => this%get_par_environment()
    mlbddc_coarse_am_i_l1_task = par_environment%am_i_l1_task()
  end function mlbddc_coarse_am_i_l1_task
  
end module mlbddc_names
