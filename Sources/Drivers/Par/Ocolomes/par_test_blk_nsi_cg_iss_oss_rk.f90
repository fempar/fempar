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
module my_nonlinear_operator_names
  use serial_names
  use par_names
  implicit none
# include "debug.i90"
  private

  type, extends(par_picard_nonlinear_operator_t) :: momentum_operator_t 
     type(par_block_matrix_t)                                  :: p_block_matrix_l
     type(par_block_matrix_t)                                  :: p_block_matrix_nl
     type(par_block_matrix_t)                                  :: p_block_matrix_m
     type(par_block_matrix_t)                                  :: p_block_u_matrix
     type(par_block_matrix_t)                                  :: p_mass_x_matrix
     type(par_block_vector_t)                                  :: p_block_vector_l
     type(par_block_vector_t)                                  :: p_block_vector_nl
     type(par_block_vector_t)                                  :: p_block_vector_m
     type(par_block_vector_t)                                  :: p_block_vector
     type(par_block_vector_t)                                  :: p_block_unknown
     type(block_operator_t)                                    :: Am, Al, Anl
     type(abs_operator_t)                                      :: A_uu, scal_Am, scal_Al, scal_Anl
     type(block_operator_t)                                    :: A_momentum
     type(block_operand_t)                                     :: b_momentum, bm, bl, bnl, brhs
     type(block_operand_t)                                     :: x_momentum
     type(block_preconditioner_l_t)                            :: M_momentum
     type(par_preconditioner_dd_mlevel_bddc_t)                 :: p_mlevel_bddc_u
     type(par_preconditioner_dd_diagonal_t)                    :: p_diagonal_x
     type(par_preconditioner_dd_mlevel_bddc_params_t)          :: p_mlevel_bddc_pars_u
     type(par_preconditioner_dd_mlevel_bddc_params_t), pointer :: point_to_p_mlevel_bddc_pars_u 
     type(inverse_operator_t)                                  :: inv_A_uu
     type(solver_control_t)                                    :: sctrl_A_uu
     real(rp)                                                  :: dt, aii
   contains
     procedure :: build => build_momentum_operator
     procedure :: free  => free_momentum_operator
     procedure :: compute_vector => compute_vector_momentum_operator
  end type momentum_operator_t
  
  type, extends(par_picard_nonlinear_operator_t) :: pressure_operator_t 
     type(par_block_matrix_t)                                  :: p_mass_u_matrix
     type(par_block_matrix_t)                                  :: p_lapla_p_matrix
     type(par_block_matrix_t)                                  :: p_block_matrix
     type(par_block_vector_t)                                  :: p_block_vector
     type(par_block_vector_t)                                  :: p_block_unknown
     type(block_operator_t)                                    :: A_pressure
     type(block_operand_t)                                     :: b_pressure
     type(block_operand_t)                                     :: x_pressure
     type(block_preconditioner_l_t)                            :: M_pressure
     type(par_preconditioner_dd_mlevel_bddc_t)                 :: p_mlevel_bddc_p
     type(par_preconditioner_dd_diagonal_t)                    :: p_diagonal_w
     type(par_preconditioner_dd_mlevel_bddc_params_t)          :: p_mlevel_bddc_pars_p
     type(par_preconditioner_dd_mlevel_bddc_params_t), pointer :: point_to_p_mlevel_bddc_pars_p  
   contains
     procedure :: build => build_pressure_operator
     procedure :: free  => free_pressure_operator
     procedure :: compute_vector => compute_vector_pressure_operator
  end type pressure_operator_t
  
  type, extends(par_picard_nonlinear_operator_t) :: momentum_update_operator_t 
     type(par_block_matrix_t)               :: p_mass_u_matrix
     type(par_block_matrix_t)               :: p_block_matrix
     type(par_block_vector_t)               :: p_block_vector
     type(par_block_vector_t)               :: p_block_unknown
     type(block_operator_t)                 :: A_update_u
     type(block_operand_t)                  :: b_update_u
     type(block_operand_t)                  :: x_update_u
     type(block_preconditioner_l_t)         :: M_update_u
     type(par_preconditioner_dd_diagonal_t) :: p_diagonal_u
   contains
     procedure :: build => build_momentum_update_operator
     procedure :: free  => free_momentum_update_operator
     procedure :: compute_vector => compute_vector_momentum_update_operator
  end type momentum_update_operator_t
  
  type, extends(par_picard_nonlinear_operator_t) :: projection_update_operator_t 
     type(par_block_matrix_t)               :: p_mass_x_matrix
     type(par_block_matrix_t)               :: p_block_matrix
     type(par_block_vector_t)               :: p_block_vector
     type(par_block_vector_t)               :: p_block_unknown
     type(block_operator_t)                 :: A_update_x
     type(block_operand_t)                  :: b_update_x
     type(block_operand_t)                  :: x_update_x
     type(block_preconditioner_l_t)         :: M_update_x
     type(par_preconditioner_dd_diagonal_t) :: p_diagonal_x
   contains
     procedure :: build => build_projection_update_operator
     procedure :: free  => free_projection_update_operator
     procedure :: compute_vector => compute_vector_projection_update_operator
  end type projection_update_operator_t

 ! Types
  public :: momentum_operator_t, pressure_operator_t, momentum_update_operator_t, &
       &    projection_update_operator_t
  
contains

  !==================================================================================================
  subroutine build_momentum_operator(nlop,p_blk_graph,p_env,ndime,max_iter,nltol)
    implicit none
    class(momentum_operator_t), target, intent(inout) :: nlop
    type(par_block_graph_t)           , intent(in)    :: p_blk_graph
    class(par_environment_t)          , intent(in)    :: p_env
    integer(ip)                       , intent(in)    :: ndime,max_iter
    real(rp)                          , intent(in)    :: nltol
    ! Locals
    integer(ip) :: i,num_levels,istat
    type(par_graph_t), pointer :: p_p_graph

    ! Set number of levels
    p_p_graph => p_blk_graph%get_block(1,1)
    num_levels = p_p_graph%p_env%num_levels

    ! Modify default nonlinear parameters
    nlop%max_iter = max_iter
    nlop%nltol = nltol

    ! Modify default inverse_operator solver parameters
    nlop%sctrl_A_uu%method = direct
    
    ! Initialize parameters
    nlop%dt  = 0.0_rp
    nlop%aii = 0.0_rp

    ! Allocate matrices and vectors
    ! Mass matrix
    call nlop%p_block_matrix_m%alloc(p_blk_graph)
    call nlop%p_block_matrix_m%set_block_to_zero(1,2)
    call nlop%p_block_matrix_m%set_block_to_zero(1,3)
    call nlop%p_block_matrix_m%set_block_to_zero(2,1)
    call nlop%p_block_matrix_m%set_block_to_zero(2,2)
    call nlop%p_block_matrix_m%set_block_to_zero(2,3)
    call nlop%p_block_matrix_m%set_block_to_zero(3,1)
    call nlop%p_block_matrix_m%set_block_to_zero(3,2)
    call nlop%p_block_matrix_m%set_block_to_zero(3,3)
    nlop%p_block_matrix_m%fill_values_stage = update_transient
    ! Linear matrix
    call nlop%p_block_matrix_l%alloc(p_blk_graph)
    call nlop%p_block_matrix_l%set_block_to_zero(1,2)
    call nlop%p_block_matrix_l%set_block_to_zero(1,3)
    call nlop%p_block_matrix_l%set_block_to_zero(2,1)
    call nlop%p_block_matrix_l%set_block_to_zero(2,2)
    call nlop%p_block_matrix_l%set_block_to_zero(2,3)
    call nlop%p_block_matrix_l%set_block_to_zero(3,1)
    call nlop%p_block_matrix_l%set_block_to_zero(3,2)
    call nlop%p_block_matrix_l%set_block_to_zero(3,3)
    nlop%p_block_matrix_l%fill_values_stage = update_constant
    ! Nonlinear matrix
    call nlop%p_block_matrix_nl%alloc(p_blk_graph)
    call nlop%p_block_matrix_nl%set_block_to_zero(1,2)
    call nlop%p_block_matrix_nl%set_block_to_zero(2,1)
    call nlop%p_block_matrix_nl%set_block_to_zero(2,2)
    call nlop%p_block_matrix_nl%set_block_to_zero(2,3)
    call nlop%p_block_matrix_nl%set_block_to_zero(3,2)
    nlop%p_block_matrix_nl%fill_values_stage = update_nonlinear
    ! Auxiliar matrix (Block-U preconditioner)
    call nlop%p_block_u_matrix%alloc(p_blk_graph)
    call nlop%p_block_u_matrix%set_block_to_zero(1,2)
    call nlop%p_block_u_matrix%set_block_to_zero(1,3)
    call nlop%p_block_u_matrix%set_block_to_zero(2,1)
    call nlop%p_block_u_matrix%set_block_to_zero(2,2)
    call nlop%p_block_u_matrix%set_block_to_zero(2,3)
    call nlop%p_block_u_matrix%set_block_to_zero(3,1)
    call nlop%p_block_u_matrix%set_block_to_zero(3,2)
    call nlop%p_block_u_matrix%set_block_to_zero(3,3)
    nlop%p_block_u_matrix%fill_values_stage = update_nonlinear
    ! Auxiliar matrix (Block-X preconditioner)
    call nlop%p_mass_x_matrix%alloc(p_blk_graph)
    call nlop%p_mass_x_matrix%set_block_to_zero(1,1)
    call nlop%p_mass_x_matrix%set_block_to_zero(1,2)
    call nlop%p_mass_x_matrix%set_block_to_zero(1,3)
    call nlop%p_mass_x_matrix%set_block_to_zero(2,1)
    call nlop%p_mass_x_matrix%set_block_to_zero(2,2)
    call nlop%p_mass_x_matrix%set_block_to_zero(2,3)
    call nlop%p_mass_x_matrix%set_block_to_zero(3,1)
    call nlop%p_mass_x_matrix%set_block_to_zero(3,2)
    nlop%p_mass_x_matrix%fill_values_stage = update_nonlinear
    ! Mass BCs vector
    call nlop%p_block_vector_m%alloc(p_blk_graph)
    nlop%p_block_vector_m%blocks(1)%state = part_summed
    nlop%p_block_vector_m%blocks(2)%state = part_summed
    nlop%p_block_vector_m%blocks(3)%state = part_summed
    ! Linear BCs vector
    call nlop%p_block_vector_l%alloc(p_blk_graph)
    nlop%p_block_vector_l%blocks(1)%state = part_summed
    nlop%p_block_vector_l%blocks(2)%state = part_summed
    nlop%p_block_vector_l%blocks(3)%state = part_summed
    ! Nonlinear BCs vector
    call nlop%p_block_vector_nl%alloc(p_blk_graph)
    nlop%p_block_vector_nl%blocks(1)%state = part_summed
    nlop%p_block_vector_nl%blocks(2)%state = part_summed
    nlop%p_block_vector_nl%blocks(3)%state = part_summed
    ! RHS vector
    call nlop%p_block_vector%alloc(p_blk_graph)
    nlop%p_block_vector%blocks(1)%state = part_summed
    nlop%p_block_vector%blocks(2)%state = part_summed
    nlop%p_block_vector%blocks(3)%state = part_summed
    call nlop%p_block_vector%init(0.0_rp)
    ! Unknown vector
    call nlop%p_block_unknown%alloc(p_blk_graph)
    nlop%p_block_unknown%blocks(1)%state = full_summed
    nlop%p_block_unknown%blocks(2)%state = full_summed
    nlop%p_block_unknown%blocks(3)%state = full_summed

    ! Define (recursive) parameters for K^-1
    nlop%point_to_p_mlevel_bddc_pars_u => nlop%p_mlevel_bddc_pars_u
    do i=1, num_levels-1
       nlop%point_to_p_mlevel_bddc_pars_u%ndime            = ndime
       nlop%point_to_p_mlevel_bddc_pars_u%unknowns         = all_unknowns
       nlop%point_to_p_mlevel_bddc_pars_u%pad_collectives  = pad
       nlop%point_to_p_mlevel_bddc_pars_u%projection       = galerkin                     
       nlop%point_to_p_mlevel_bddc_pars_u%subd_elmat_calc  = phit_minus_c_i_t_lambda            !default  
       nlop%point_to_p_mlevel_bddc_pars_u%correction_mode  = additive_symmetric                 !default 
       nlop%point_to_p_mlevel_bddc_pars_u%nn_sys_sol_strat = corners_rest_part_solve_expl_schur ! default 
       if(ndime==3) then
          nlop%point_to_p_mlevel_bddc_pars_u%kind_coarse_dofs = corners_edges_and_faces
       else
          nlop%point_to_p_mlevel_bddc_pars_u%kind_coarse_dofs = corners_and_edges
       end if
       if ( i < num_levels-1 ) then
          nlop%point_to_p_mlevel_bddc_pars_u%co_sys_sol_strat     = recursive_bddc
          nlop%point_to_p_mlevel_bddc_pars_u%ppars_harm%type      = pardiso_mkl_prec !umfpack_prec 
          nlop%point_to_p_mlevel_bddc_pars_u%ppars_dirichlet%type = pardiso_mkl_prec !umfpack_prec   
          if ( i == 1 ) then
             nlop%point_to_p_mlevel_bddc_pars_u%spars_coarse%method = direct
             nlop%point_to_p_mlevel_bddc_pars_u%spars_coarse%itmax  = 200
             nlop%point_to_p_mlevel_bddc_pars_u%spars_coarse%rtol   = 1.0e-20
             nlop%point_to_p_mlevel_bddc_pars_u%spars_coarse%trace  = 1
             nlop%point_to_p_mlevel_bddc_pars_u%correction_mode     = additive
          end if
          allocate(nlop%point_to_p_mlevel_bddc_pars_u%ppars_coarse_bddc, stat = istat)
          check(istat==0)
          nlop%point_to_p_mlevel_bddc_pars_u => nlop%point_to_p_mlevel_bddc_pars_u%ppars_coarse_bddc
       else
          nlop%point_to_p_mlevel_bddc_pars_u%co_sys_sol_strat         = serial_gather
          nlop%point_to_p_mlevel_bddc_pars_u%ppars_harm%type          = pardiso_mkl_prec !umfpack_prec  
          nlop%point_to_p_mlevel_bddc_pars_u%ppars_dirichlet%type     = pardiso_mkl_prec !umfpack_prec  
          nlop%point_to_p_mlevel_bddc_pars_u%ppars_coarse_serial%type = pardiso_mkl_prec !umfpack_prec  
          nullify ( nlop%point_to_p_mlevel_bddc_pars_u%ppars_coarse_bddc )
       end if
    end do
    nlop%point_to_p_mlevel_bddc_pars_u => nlop%p_mlevel_bddc_pars_u
    do i=1, num_levels-1
       nlop%point_to_p_mlevel_bddc_pars_u => nlop%point_to_p_mlevel_bddc_pars_u%ppars_coarse_bddc
    end do
    
    ! Construct U-preconditioner (K^-1)
    call par_preconditioner_dd_mlevel_bddc_create(nlop%p_block_u_matrix%get_block(1,1),nlop%p_mlevel_bddc_u, &
         &                                        nlop%p_mlevel_bddc_pars_u)
    call par_preconditioner_dd_mlevel_bddc_ass_struct(nlop%p_block_u_matrix%get_block(1,1),nlop%p_mlevel_bddc_u)
    nlop%p_mlevel_bddc_u%fill_values_stage = update_transient
    nlop%p_mlevel_bddc_u%free_values_stage = update_transient

    ! Construct X-preconditioner (Mx^-1)
    call par_preconditioner_dd_diagonal_create(nlop%p_mass_x_matrix%get_block(3,3),nlop%p_diagonal_x)
    call par_preconditioner_dd_diagonal_ass_struct(nlop%p_mass_x_matrix%get_block(3,3),nlop%p_diagonal_x)
    nlop%p_diagonal_x%fill_values_stage = update_transient!nonlinear
    nlop%p_diagonal_x%free_values_stage = update_transient!nonlinear

    ! Create Mass operator
    call nlop%Am%create(2,2)
    call nlop%Am%set_block(1,1,nlop%p_block_matrix_m%get_block(1,1))
    call nlop%Am%set_block_to_zero(1,2)
    call nlop%Am%set_block_to_zero(2,1)
    call nlop%Am%set_block_to_zero(2,2)

    ! Create Linear operator
    call nlop%Al%create(2,2)
    call nlop%Al%set_block(1,1,nlop%p_block_matrix_l%get_block(1,1))
    call nlop%Al%set_block_to_zero(1,2)
    call nlop%Al%set_block_to_zero(2,1)
    call nlop%Al%set_block_to_zero(2,2)

    ! Create Nonlinear operator
    call nlop%Anl%create(2,2)
    call nlop%Anl%set_block(1,1,nlop%p_block_matrix_nl%get_block(1,1))
    call nlop%Anl%set_block(1,2,nlop%p_block_matrix_nl%get_block(1,3))
    call nlop%Anl%set_block(2,1,nlop%p_block_matrix_nl%get_block(3,1))
    call nlop%Anl%set_block(2,2,nlop%p_block_matrix_nl%get_block(3,3))

    ! Construct scalar products
    nlop%scal_Am = nlop%dt * nlop%Am%get_block(1,1)
    nlop%scal_Al = nlop%aii * nlop%Al%get_block(1,1)
    nlop%scal_Anl = nlop%aii * nlop%Anl%get_block(1,1)

    ! Construct abstract operator for UU block
    nlop%A_uu = nlop%scal_Am + nlop%scal_Al + nlop%scal_Anl

    ! Create inverse operator for Block UU
    call nlop%inv_A_uu%create(nlop%A_uu,nlop%p_mlevel_bddc_u,nlop%sctrl_A_uu,p_env)
    
    ! Construct global operator
    call nlop%A_momentum%create(2,2)
    call nlop%A_momentum%set_block(1,1,nlop%A_uu)
    call nlop%A_momentum%set_block(1,2,nlop%aii*nlop%p_block_matrix_nl%get_block(1,3))
!!$    call nlop%A_momentum%set_block(2,1,nlop%p_block_matrix_nl%get_block(3,1))
!!$    call nlop%A_momentum%set_block(2,2,nlop%p_block_matrix_nl%get_block(3,3))
    call nlop%A_momentum%set_block(2,1,nlop%aii*nlop%p_block_matrix_nl%get_block(3,1))
    call nlop%A_momentum%set_block(2,2,nlop%aii*nlop%p_block_matrix_nl%get_block(3,3))  
    
    ! Create Mass BCs operand
    call nlop%bm%create(2)
    call nlop%bm%set_block(1,nlop%p_block_vector_m%blocks(1))
    call nlop%bm%set_block(2,nlop%p_block_vector_m%blocks(3))

    ! Create Linear BCs operand
    call nlop%bl%create(2)
    call nlop%bl%set_block(1,nlop%p_block_vector_l%blocks(1))
    call nlop%bl%set_block(2,nlop%p_block_vector_l%blocks(3))

    ! Create Nonlinear BCs operand
    call nlop%bnl%create(2)
    call nlop%bnl%set_block(1,nlop%p_block_vector_nl%blocks(1))
    call nlop%bnl%set_block(2,nlop%p_block_vector_nl%blocks(3))

    ! Create RHS operand
    call nlop%brhs%create(2)
    call nlop%brhs%set_block(1,nlop%p_block_vector%blocks(1))
    call nlop%brhs%set_block(2,nlop%p_block_vector%blocks(3))

    ! Construct solution operand
    call nlop%x_momentum%create(2)
    call nlop%x_momentum%set_block(1,nlop%p_block_unknown%blocks(1))
    call nlop%x_momentum%set_block(2,nlop%p_block_unknown%blocks(3))

    ! Construct global block preconditioner
    call nlop%M_momentum%create(2)
    call nlop%M_momentum%set_block(1,1,nlop%inv_A_uu)
    call nlop%M_momentum%set_block(2,1,nlop%aii*nlop%p_block_matrix_nl%get_block(3,1))
!!$    call nlop%M_momentum%set_block(2,1,nlop%p_block_matrix_nl%get_block(3,1))
    call nlop%M_momentum%set_block(2,2,nlop%p_diagonal_x)
    
    ! Fill picard nonlinear operator pointers
    call nlop%create(nlop%A_momentum, nlop%M_momentum, nlop%b_momentum, nlop%x_momentum, &
         &           nlop%p_block_unknown,                                               &
         &           A_int_c = nlop%p_block_matrix_l, b_int_c = nlop%p_block_vector_l,   &
         &           A_int_t = nlop%p_block_matrix_m, b_int_t = nlop%p_block_vector_m,   &
         &           A_int_n = nlop%p_block_matrix_nl,b_int_n = nlop%p_block_vector_nl)
    
  end subroutine build_momentum_operator

  !==================================================================================================
  subroutine compute_vector_momentum_operator(nlop)
    implicit none
    class(momentum_operator_t), intent(inout) :: nlop
    
    ! Construct global operand
    nlop%b_momentum = nlop%dt * nlop%bm + nlop%aii * nlop%bl + nlop%aii * nlop%bnl + nlop%brhs
    
  end subroutine compute_vector_momentum_operator

  !==================================================================================================
  subroutine build_pressure_operator(nlop,p_blk_graph,p_env,ndime)
    implicit none
    class(pressure_operator_t), target, intent(inout) :: nlop
    type(par_block_graph_t)           , intent(in)    :: p_blk_graph
    class(par_environment_t)          , intent(in)    :: p_env
    integer(ip)                       , intent(in)    :: ndime
    ! Locals
    integer(ip)                 :: i,num_levels,istat
    type(par_graph_t) , pointer :: p_p_graph
    type(par_matrix_t), pointer :: aux

    ! Set number of levels
    p_p_graph => p_blk_graph%get_block(1,1)
    num_levels = p_p_graph%p_env%num_levels

    ! Allocate matrices and vectors
    call nlop%p_block_matrix%alloc(p_blk_graph)
    call nlop%p_block_matrix%set_block_to_zero(1,3)
    call nlop%p_block_matrix%set_block_to_zero(2,2)
    call nlop%p_block_matrix%set_block_to_zero(2,3)
    call nlop%p_block_matrix%set_block_to_zero(3,1)
    call nlop%p_block_matrix%set_block_to_zero(3,2)
    call nlop%p_block_matrix%set_block_to_zero(3,3)
    nlop%p_block_matrix%fill_values_stage = update_transient
    ! RHS vector
    call nlop%p_block_vector%alloc(p_blk_graph)
    nlop%p_block_vector%blocks(1)%state = part_summed
    nlop%p_block_vector%blocks(2)%state = part_summed
    nlop%p_block_vector%blocks(3)%state = part_summed
    call nlop%p_block_vector%init(0.0_rp)
    ! Unknown vector
    call nlop%p_block_unknown%alloc(p_blk_graph)
    nlop%p_block_unknown%blocks(1)%state = full_summed
    nlop%p_block_unknown%blocks(2)%state = full_summed
    nlop%p_block_unknown%blocks(3)%state = full_summed

    ! Auxiliar matrices
    call nlop%p_mass_u_matrix%alloc(p_blk_graph)
    call nlop%p_mass_u_matrix%set_block_to_zero(1,2)
    call nlop%p_mass_u_matrix%set_block_to_zero(1,3)
    call nlop%p_mass_u_matrix%set_block_to_zero(2,1)
    call nlop%p_mass_u_matrix%set_block_to_zero(2,2)
    call nlop%p_mass_u_matrix%set_block_to_zero(2,3)
    call nlop%p_mass_u_matrix%set_block_to_zero(3,1)
    call nlop%p_mass_u_matrix%set_block_to_zero(3,2)
    call nlop%p_mass_u_matrix%set_block_to_zero(3,3)
    call nlop%p_lapla_p_matrix%alloc(p_blk_graph)
    call nlop%p_lapla_p_matrix%set_block_to_zero(1,1)
    call nlop%p_lapla_p_matrix%set_block_to_zero(1,2)
    call nlop%p_lapla_p_matrix%set_block_to_zero(1,3)
    call nlop%p_lapla_p_matrix%set_block_to_zero(2,1)
    call nlop%p_lapla_p_matrix%set_block_to_zero(2,3)
    call nlop%p_lapla_p_matrix%set_block_to_zero(3,1)
    call nlop%p_lapla_p_matrix%set_block_to_zero(3,2)
    call nlop%p_lapla_p_matrix%set_block_to_zero(3,3)

    ! Define (recursive) parameters for Mp^-1
    nlop%point_to_p_mlevel_bddc_pars_p => nlop%p_mlevel_bddc_pars_p
    do i=1, num_levels-1
       nlop%point_to_p_mlevel_bddc_pars_p%ndime            = ndime
       nlop%point_to_p_mlevel_bddc_pars_p%unknowns         = all_unknowns
       nlop%point_to_p_mlevel_bddc_pars_p%pad_collectives  = pad
       nlop%point_to_p_mlevel_bddc_pars_p%projection       = galerkin                     
       nlop%point_to_p_mlevel_bddc_pars_p%subd_elmat_calc  = phit_minus_c_i_t_lambda            !default  
       nlop%point_to_p_mlevel_bddc_pars_p%correction_mode  = additive_symmetric                 !default 
       nlop%point_to_p_mlevel_bddc_pars_p%nn_sys_sol_strat = corners_rest_part_solve_expl_schur ! default 
       if(ndime==3) then
          nlop%point_to_p_mlevel_bddc_pars_p%kind_coarse_dofs = corners_edges_and_faces
       else
          nlop%point_to_p_mlevel_bddc_pars_p%kind_coarse_dofs = corners_and_edges
       end if
       if ( i < num_levels-1 ) then
          nlop%point_to_p_mlevel_bddc_pars_p%co_sys_sol_strat     = recursive_bddc
          nlop%point_to_p_mlevel_bddc_pars_p%ppars_harm%type      = pardiso_mkl_prec !umfpack_prec 
          nlop%point_to_p_mlevel_bddc_pars_p%ppars_dirichlet%type = pardiso_mkl_prec !umfpack_prec   
          if ( i == 1 ) then
             nlop%point_to_p_mlevel_bddc_pars_p%spars_coarse%method = direct
             nlop%point_to_p_mlevel_bddc_pars_p%spars_coarse%itmax  = 200
             nlop%point_to_p_mlevel_bddc_pars_p%spars_coarse%rtol   = 1.0e-20
             nlop%point_to_p_mlevel_bddc_pars_p%spars_coarse%trace  = 1
             nlop%point_to_p_mlevel_bddc_pars_p%correction_mode     = additive
          end if
          allocate(nlop%point_to_p_mlevel_bddc_pars_p%ppars_coarse_bddc, stat = istat)
          check(istat==0)
          nlop%point_to_p_mlevel_bddc_pars_p => nlop%point_to_p_mlevel_bddc_pars_p%ppars_coarse_bddc
       else
          nlop%point_to_p_mlevel_bddc_pars_p%co_sys_sol_strat         = serial_gather
          nlop%point_to_p_mlevel_bddc_pars_p%ppars_harm%type          = pardiso_mkl_prec !umfpack_prec  
          nlop%point_to_p_mlevel_bddc_pars_p%ppars_dirichlet%type     = pardiso_mkl_prec !umfpack_prec  
          nlop%point_to_p_mlevel_bddc_pars_p%ppars_coarse_serial%type = pardiso_mkl_prec !umfpack_prec  
          nullify ( nlop%point_to_p_mlevel_bddc_pars_p%ppars_coarse_bddc )
       end if
    end do
    nlop%point_to_p_mlevel_bddc_pars_p => nlop%p_mlevel_bddc_pars_p
    do i=1, num_levels-1
       nlop%point_to_p_mlevel_bddc_pars_p => nlop%point_to_p_mlevel_bddc_pars_p%ppars_coarse_bddc
    end do

    ! Construct W-preconditioner (Mu^-1)
    call par_preconditioner_dd_diagonal_create(nlop%p_mass_u_matrix%get_block(1,1),nlop%p_diagonal_w)
    call par_preconditioner_dd_diagonal_ass_struct(nlop%p_mass_u_matrix%get_block(1,1),nlop%p_diagonal_w)
    nlop%p_diagonal_w%fill_values_stage = update_constant
    nlop%p_diagonal_w%free_values_stage = update_constant

    ! Construct P-preconditioner (Lp^-1)
    call par_preconditioner_dd_mlevel_bddc_create(nlop%p_lapla_p_matrix%get_block(2,2),nlop%p_mlevel_bddc_p, &
         &                                        nlop%p_mlevel_bddc_pars_p)
    call par_preconditioner_dd_mlevel_bddc_ass_struct(nlop%p_lapla_p_matrix%get_block(2,2),nlop%p_mlevel_bddc_p)
    nlop%p_mlevel_bddc_p%fill_values_stage = update_constant
    nlop%p_mlevel_bddc_p%free_values_stage = update_constant
    
    ! Create W-P Block operator
    call nlop%A_pressure%create(2,2)
    call nlop%A_pressure%set_block(1,1,nlop%p_block_matrix%get_block(1,1))
    call nlop%A_pressure%set_block(1,2,nlop%p_block_matrix%get_block(1,2))
    call nlop%A_pressure%set_block(2,1,nlop%p_block_matrix%get_block(2,1))
    call nlop%A_pressure%set_block_to_zero(2,2)

    ! Create W-P Block operand
    call nlop%b_pressure%create(2)
    call nlop%b_pressure%set_block(1,nlop%p_block_vector%blocks(1))
    call nlop%b_pressure%set_block(2,nlop%p_block_vector%blocks(2))
    call nlop%x_pressure%create(2)
    call nlop%x_pressure%set_block(1,nlop%p_block_unknown%blocks(1))
    call nlop%x_pressure%set_block(2,nlop%p_block_unknown%blocks(2))

    ! Create W-P Block preconditioner
    call nlop%M_pressure%create(2)
    call nlop%M_pressure%set_block(1,1,nlop%p_diagonal_w)
    call nlop%M_pressure%set_block(2,1,nlop%p_block_matrix%get_block(2,1))
    call nlop%M_pressure%set_block(2,2,nlop%p_mlevel_bddc_p)

    ! Fill picard nonlinear operator pointers
    call nlop%create(nlop%A_pressure, nlop%M_pressure, nlop%b_pressure, nlop%x_pressure, &
         &           nlop%p_block_unknown, A_int_t=nlop%p_block_matrix, b_int_t=nlop%p_block_vector)
    
  end subroutine build_pressure_operator

  !==================================================================================================
  subroutine compute_vector_pressure_operator(nlop)
    implicit none
    class(pressure_operator_t), intent(inout) :: nlop
    
    ! Dummy
    
  end subroutine compute_vector_pressure_operator

  !==================================================================================================
  subroutine build_momentum_update_operator(nlop,p_blk_graph)
    implicit none
    class(momentum_update_operator_t), intent(inout) :: nlop
    type(par_block_graph_t)          , intent(in)    :: p_blk_graph

    ! Allocate integration matrices and vectors
    call nlop%p_block_matrix%alloc(p_blk_graph)
    call nlop%p_block_matrix%set_block_to_zero(1,2)
    call nlop%p_block_matrix%set_block_to_zero(1,3)
    call nlop%p_block_matrix%set_block_to_zero(2,1)
    call nlop%p_block_matrix%set_block_to_zero(2,2)
    call nlop%p_block_matrix%set_block_to_zero(2,3)
    call nlop%p_block_matrix%set_block_to_zero(3,1)
    call nlop%p_block_matrix%set_block_to_zero(3,2)
    call nlop%p_block_matrix%set_block_to_zero(3,3)
    nlop%p_block_matrix%fill_values_stage = update_transient
    ! RHS vector
    call nlop%p_block_vector%alloc(p_blk_graph)
    nlop%p_block_vector%blocks(1)%state = part_summed
    nlop%p_block_vector%blocks(2)%state = part_summed
    nlop%p_block_vector%blocks(3)%state = part_summed
    call nlop%p_block_vector%init(0.0_rp)
    ! Unknown vector
    call nlop%p_block_unknown%alloc(p_blk_graph)
    nlop%p_block_unknown%blocks(1)%state = full_summed
    nlop%p_block_unknown%blocks(2)%state = full_summed
    nlop%p_block_unknown%blocks(3)%state = full_summed

    ! Auxiliar matrices
    call nlop%p_mass_u_matrix%alloc(p_blk_graph)
    call nlop%p_mass_u_matrix%set_block_to_zero(1,2)
    call nlop%p_mass_u_matrix%set_block_to_zero(1,3)
    call nlop%p_mass_u_matrix%set_block_to_zero(2,1)
    call nlop%p_mass_u_matrix%set_block_to_zero(2,2)
    call nlop%p_mass_u_matrix%set_block_to_zero(2,3)
    call nlop%p_mass_u_matrix%set_block_to_zero(3,1)
    call nlop%p_mass_u_matrix%set_block_to_zero(3,2)
    call nlop%p_mass_u_matrix%set_block_to_zero(3,3)
    
    ! Construct U-preconditioner (Mu^-1)
    call par_preconditioner_dd_diagonal_create(nlop%p_mass_u_matrix%get_block(1,1),nlop%p_diagonal_u)
    call par_preconditioner_dd_diagonal_ass_struct(nlop%p_mass_u_matrix%get_block(1,1),nlop%p_diagonal_u)
    nlop%p_diagonal_u%fill_values_stage = update_transient
    nlop%p_diagonal_u%free_values_stage = update_transient
    
    ! Create U Block operator
    call nlop%A_update_u%create(1,1)
    call nlop%A_update_u%set_block(1,1,nlop%p_block_matrix%get_block(1,1))

    ! Create U Block operand
    call nlop%b_update_u%create(1)
    call nlop%b_update_u%set_block(1,nlop%p_block_vector%blocks(1))
    call nlop%x_update_u%create(1)
    call nlop%x_update_u%set_block(1,nlop%p_block_unknown%blocks(1))

    ! Create U Block preconditioner
    call nlop%M_update_u%create(1)
    call nlop%M_update_u%set_block(1,1,nlop%p_diagonal_u)

    ! Fill picard nonlinear operator pointers
    call nlop%create(nlop%A_update_u, nlop%M_update_u, nlop%b_update_u, nlop%x_update_u, &
         &           nlop%p_block_unknown, A_int_t=nlop%p_block_matrix, b_int_t=nlop%p_block_vector)
    
  end subroutine build_momentum_update_operator

  !==================================================================================================
  subroutine compute_vector_momentum_update_operator(nlop)
    implicit none
    class(momentum_update_operator_t), intent(inout) :: nlop
    
    ! Dummy
    
  end subroutine compute_vector_momentum_update_operator

 !==================================================================================================
  subroutine build_projection_update_operator(nlop,p_blk_graph)
    implicit none
    class(projection_update_operator_t), intent(inout) :: nlop
    type(par_block_graph_t)            , intent(in)    :: p_blk_graph

    ! Allocate matrices and vectors
    call nlop%p_block_matrix%alloc(p_blk_graph)
    call nlop%p_block_matrix%set_block_to_zero(1,1)
    call nlop%p_block_matrix%set_block_to_zero(1,2)
    call nlop%p_block_matrix%set_block_to_zero(1,3)
    call nlop%p_block_matrix%set_block_to_zero(2,1)
    call nlop%p_block_matrix%set_block_to_zero(2,2)
    call nlop%p_block_matrix%set_block_to_zero(2,3)
    call nlop%p_block_matrix%set_block_to_zero(3,1)
    call nlop%p_block_matrix%set_block_to_zero(3,2)
    nlop%p_block_matrix%fill_values_stage = update_transient
    ! RHS vector
    call nlop%p_block_vector%alloc(p_blk_graph)
    nlop%p_block_vector%blocks(1)%state = part_summed
    nlop%p_block_vector%blocks(2)%state = part_summed
    nlop%p_block_vector%blocks(3)%state = part_summed
    call nlop%p_block_vector%init(0.0_rp)
    ! Unknown vector
    call nlop%p_block_unknown%alloc(p_blk_graph)
    nlop%p_block_unknown%blocks(1)%state = full_summed
    nlop%p_block_unknown%blocks(2)%state = full_summed
    nlop%p_block_unknown%blocks(3)%state = full_summed

    ! Auxiliar matrices
    call nlop%p_mass_x_matrix%alloc(p_blk_graph)
    call nlop%p_mass_x_matrix%set_block_to_zero(1,1)
    call nlop%p_mass_x_matrix%set_block_to_zero(1,2)
    call nlop%p_mass_x_matrix%set_block_to_zero(1,3)
    call nlop%p_mass_x_matrix%set_block_to_zero(2,1)
    call nlop%p_mass_x_matrix%set_block_to_zero(2,2)
    call nlop%p_mass_x_matrix%set_block_to_zero(2,3)
    call nlop%p_mass_x_matrix%set_block_to_zero(3,1)
    call nlop%p_mass_x_matrix%set_block_to_zero(3,2)

    ! Construct X-preconditioner (Mx^-1)
    call par_preconditioner_dd_diagonal_create(nlop%p_mass_x_matrix%get_block(3,3),nlop%p_diagonal_x)
    call par_preconditioner_dd_diagonal_ass_struct(nlop%p_mass_x_matrix%get_block(3,3),nlop%p_diagonal_x)
    nlop%p_diagonal_x%fill_values_stage = update_transient
    nlop%p_diagonal_x%free_values_stage = update_transient
    
    ! Create U Block operator
    call nlop%A_update_x%create(1,1)
    call nlop%A_update_x%set_block(1,1,nlop%p_block_matrix%get_block(3,3))

    ! Create U Block operand
    call nlop%b_update_x%create(1)
    call nlop%b_update_x%set_block(1,nlop%p_block_vector%blocks(3))
    call nlop%x_update_x%create(1)
    call nlop%x_update_x%set_block(1,nlop%p_block_unknown%blocks(3))

    ! Create U Block preconditioner
    call nlop%M_update_x%create(1)
    call nlop%M_update_x%set_block(1,1,nlop%p_diagonal_x)
    nlop%M_update_x%fill_values_stage = update_nonlinear!transient

    ! Fill picard nonlinear operator pointers
    call nlop%create(nlop%A_update_x, nlop%M_update_x, nlop%b_update_x, nlop%x_update_x, &
         &           nlop%p_block_unknown, A_int_t=nlop%p_block_matrix, b_int_t=nlop%p_block_vector)
    
  end subroutine build_projection_update_operator

  !==================================================================================================
  subroutine compute_vector_projection_update_operator(nlop)
    implicit none
    class(projection_update_operator_t), intent(inout) :: nlop
    
    ! Dummy
    
  end subroutine compute_vector_projection_update_operator

  !==================================================================================================
  subroutine free_momentum_operator(nlop)
    implicit none
    class(momentum_operator_t), intent(inout) :: nlop

    ! Unassign picard nonlinear operator
    call nlop%unassign()

    ! Destroy global Block preconditioner
    call nlop%M_momentum%destroy()

    ! Destroy global Block operand
    call nlop%b_momentum%destroy()
    call nlop%x_momentum%destroy()

    ! Destroy intermediate operands
    call nlop%bm%destroy()
    call nlop%bl%destroy()
    call nlop%bnl%destroy()
    call nlop%brhs%destroy()

    ! Destroy global Block operator
    call nlop%A_momentum%free()

    ! Destroy intermediate operators
    call nlop%Am%destroy()
    call nlop%Al%destroy()
    call nlop%Anl%destroy()

    ! Destroy X-preconditioner
    call par_preconditioner_dd_diagonal_free(nlop%p_diagonal_x,free_struct)
    call par_preconditioner_dd_diagonal_free(nlop%p_diagonal_x,free_clean)

    ! Destroy U-preconditioner
    call par_preconditioner_dd_mlevel_bddc_free(nlop%p_mlevel_bddc_u,free_struct)
    call par_preconditioner_dd_mlevel_bddc_free(nlop%p_mlevel_bddc_u,free_clean)

    ! Deallocate vectors
    call nlop%p_block_unknown%free()
    call nlop%p_block_vector%free()
    call nlop%p_block_vector_m%free()
    call nlop%p_block_vector_l%free()
    call nlop%p_block_vector_nl%free()
    
    ! Deallocate matrices
    call nlop%p_mass_x_matrix%free()
    call nlop%p_block_u_matrix%free()
    call nlop%p_block_matrix_m%free()
    call nlop%p_block_matrix_l%free()
    call nlop%p_block_matrix_nl%free()

  end subroutine free_momentum_operator

  !==================================================================================================
  subroutine free_pressure_operator(nlop)
    implicit none
    class(pressure_operator_t), intent(inout) :: nlop

    ! Unassign picard nonlinear operator
    call nlop%unassign()

    ! Destroy global Block preconditioner
    call nlop%M_pressure%destroy()

    ! Destroy global Block operand
    call nlop%b_pressure%destroy()
    call nlop%x_pressure%destroy()

    ! Destroy global Block operator
    call nlop%A_pressure%destroy()

    ! Destroy P-preconditioner
    call par_preconditioner_dd_mlevel_bddc_free(nlop%p_mlevel_bddc_p,free_struct)
    call par_preconditioner_dd_mlevel_bddc_free(nlop%p_mlevel_bddc_p,free_clean)

    ! Destroy W-preconditioner
    call par_preconditioner_dd_diagonal_free(nlop%p_diagonal_w,free_struct)
    call par_preconditioner_dd_diagonal_free(nlop%p_diagonal_w,free_clean)

    ! Deallocate vectors
    call nlop%p_block_unknown%free()
    call nlop%p_block_vector%free()
    
    ! Deallocate matrices
    call nlop%p_mass_u_matrix%free()
    call nlop%p_lapla_p_matrix%free()
    call nlop%p_block_matrix%free()

  end subroutine free_pressure_operator

  !==================================================================================================
  subroutine free_momentum_update_operator(nlop)
    implicit none
    class(momentum_update_operator_t), intent(inout) :: nlop

    ! Unassign picard nonlinear operator
    call nlop%unassign()

    ! Destroy global Block preconditioner
    call nlop%M_update_u%destroy()

    ! Destroy global Block operand
    call nlop%b_update_u%destroy()
    call nlop%x_update_u%destroy()

    ! Destroy global Block operator
    call nlop%A_update_u%destroy()

    ! Destroy U-preconditioner
    call par_preconditioner_dd_diagonal_free(nlop%p_diagonal_u,free_struct)
    call par_preconditioner_dd_diagonal_free(nlop%p_diagonal_u,free_clean)

    ! Deallocate vectors
    call nlop%p_block_unknown%free()
    call nlop%p_block_vector%free()
    
    ! Deallocate matrices
    call nlop%p_mass_u_matrix%free()
    call nlop%p_block_matrix%free()

  end subroutine free_momentum_update_operator

  !==================================================================================================
  subroutine free_projection_update_operator(nlop)
    implicit none
    class(projection_update_operator_t), intent(inout) :: nlop

    ! Unassign picard nonlinear operator
    call nlop%unassign()

    ! Destroy global Block preconditioner
    call nlop%M_update_x%destroy()

    ! Destroy global Block operand
    call nlop%b_update_x%destroy()
    call nlop%x_update_x%destroy()

    ! Destroy global Block operator
    call nlop%A_update_x%destroy()

    ! Destroy X-preconditioner
    call par_preconditioner_dd_diagonal_free(nlop%p_diagonal_x,free_struct)
    call par_preconditioner_dd_diagonal_free(nlop%p_diagonal_x,free_clean)

    ! Deallocate vectors
    call nlop%p_block_unknown%free()
    call nlop%p_block_vector%free()
    
    ! Deallocate matrices
    call nlop%p_mass_x_matrix%free()
    call nlop%p_block_matrix%free()

  end subroutine free_projection_update_operator

end module my_nonlinear_operator_names

module command_line_parameters_names
  use types_names
  use Data_Type_Command_Line_Interface
  implicit none
  private
  
  ! Types
  type par_test_blk_nsi_cg_iss_oss_rk_params_t
     ! Names
     character(len=:), allocatable :: default_authors   
     character(len=:), allocatable :: default_pwd       
     character(len=:), allocatable :: default_prefix    
     ! Geometry
     character(len=:), allocatable :: default_nelems    
     character(len=:), allocatable :: default_nparts    
     character(len=:), allocatable :: default_nsockets  
     character(len=:), allocatable :: default_ndiscret  
     character(len=:), allocatable :: default_length    
     character(len=:), allocatable :: default_periodic  
     character(len=:), allocatable :: default_origin    
     character(len=:), allocatable :: default_stretch
     character(len=:), allocatable :: default_nelbound  
     character(len=:), allocatable :: default_sizebound 
     character(len=:), allocatable :: default_matercase 
     ! Time integration
     character(len=:), allocatable :: default_initime   
     character(len=:), allocatable :: default_finaltime 
     character(len=:), allocatable :: default_timestep  
     character(len=:), allocatable :: default_rkstage   
     character(len=:), allocatable :: default_rkorder   
     character(len=:), allocatable :: default_rkflag   
     character(len=:), allocatable :: default_rkimex
     ! Solution
     character(len=:), allocatable :: default_analytical_v
     character(len=:), allocatable :: default_analytical_p
     character(len=:), allocatable :: default_analytical_t
     character(len=:), allocatable :: default_initial_cond
     ! Problem
     character(len=:), allocatable :: default_kfl_conv 
     character(len=:), allocatable :: default_kfl_skew 
     character(len=:), allocatable :: default_diffu    
     ! Discretization
     character(len=:), allocatable :: default_kfl_proj 
     character(len=:), allocatable :: default_ktauc    
     character(len=:), allocatable :: default_vel_order
     character(len=:), allocatable :: default_pre_order
     ! Solver
     character(len=:), allocatable :: default_rtol
     ! Nonlinearity
     character(len=:), allocatable :: default_nltol    
     character(len=:), allocatable :: default_max_iter 
     ! Outputs
     character(len=:), allocatable :: default_write_unkno
     character(len=:), allocatable :: default_trace_unkno
     character(len=:), allocatable :: default_write_dissipation
     character(len=:), allocatable :: default_trace_dissipation
  end type par_test_blk_nsi_cg_iss_oss_rk_params_t

  ! Types
  public :: par_test_blk_nsi_cg_iss_oss_rk_params_t
  
  ! Functions
  public :: cli_add_params,set_default_params,set_default_params_analytical,set_default_params_taylor_green

contains

  !==================================================================================================
  subroutine set_default_params(params)
    implicit none
    type(par_test_blk_nsi_cg_iss_oss_rk_params_t), intent(inout) :: params

    ! Names
    params%default_authors   = '.false.'
    params%default_pwd       = './'
    params%default_prefix    = 'par_test_blk_nsi_cg_iss_oss_rk'
    ! Geometry
    params%default_nelems    = '4 4'
    params%default_nparts    = '2 2'
    params%default_nsockets  = '1 1'
    params%default_ndiscret  = '0 0'
    params%default_length    = '1.0 1.0'
    params%default_periodic  = '0 0'
    params%default_origin    = '0.0 0.0'
    params%default_stretch   = '2.75 2.75'
    params%default_nelbound  = '0 0'
    params%default_sizebound = '0.0 0.0'
    params%default_matercase = '0'
    ! Time integration
    params%default_initime   = '0.0'
    params%default_finaltime = '1.0'
    params%default_timestep  = '1.0'
    params%default_rkstage   = '2'
    params%default_rkorder   = '1'
    params%default_rkflag    = '0'
    params%default_rkimex    = '1'
    ! Solution
    params%default_analytical_v = '0 0'
    params%default_analytical_p = '0'
    params%default_analytical_t = '0'
    params%default_initial_cond = '.false.'
    ! Problem
    params%default_kfl_conv = '1'
    params%default_kfl_skew = '0'
    params%default_diffu    = '1.0'
    ! Discretization
    params%default_kfl_proj  = '1'
    params%default_ktauc     = '4.0'
    params%default_vel_order = '2'
    params%default_pre_order = '1'
    ! Solver
    params%default_rtol    = '1.0e-14'
    ! Nonlinearity
    params%default_nltol    = '1.0e-10'
    params%default_max_iter = '20'    
    ! Outputs
    params%default_write_unkno = '.true.'
    params%default_trace_unkno = '0.0'
    params%default_write_dissipation = '.true.'
    params%default_trace_dissipation = '0.0'

  end subroutine set_default_params

  !==================================================================================================
  subroutine set_default_params_analytical(params)
    implicit none
    type(par_test_blk_nsi_cg_iss_oss_rk_params_t), intent(inout) :: params

    ! Names
    params%default_prefix = 'analytical_test'
    ! Solution
    params%default_analytical_v = '1 2'
    params%default_analytical_p = '3'
    params%default_analytical_t = '1'
    ! Discretization
    params%default_ktauc = '0.0' 
    ! Outputs
     params%default_write_dissipation = '.false.'

  end subroutine set_default_params_analytical

  !==================================================================================================
  subroutine set_default_params_taylor_green(params)
    implicit none
    type(par_test_blk_nsi_cg_iss_oss_rk_params_t), intent(inout) :: params

    ! Names
    params%default_prefix    = 'tgv'
    ! Geometry
    params%default_nelems    = '9 9 9'
    params%default_nparts    = '3 3 3'
    params%default_nsockets  = '1 1 1'
    params%default_ndiscret  = '0 0 0'
    params%default_length    = '6.283185307179586476925286766559005768394 6.283185307179586476925286766559005768394 6.283185307179586476925286766559005768394'
    params%default_periodic  = '1 1 1'
    params%default_origin    = '0.0 0.0 0.0'
    params%default_stretch   = '2.75 2.75 2.75'
    params%default_nelbound  = '0 0 0'
    params%default_sizebound = '0.0 0.0 0.0'
    ! Solution
    params% default_analytical_v = '7 8 0'
    params%default_initial_cond  = '.true.'
    ! Problem
    params%default_kfl_conv = '1'
    params%default_kfl_skew = '0'
    params%default_diffu    = '1.0/1600.0'
    ! Solver
    params%default_rtol     = '1.0e-07'
    ! Nonlinearity
    params%default_nltol    = '1.0e-05'
    ! Outputs
    params%default_write_unkno = '.true.'
    params%default_trace_unkno = '0.01'
    params%default_write_dissipation = '.true.'
    params%default_trace_dissipation = '0.01'

  end subroutine set_default_params_taylor_green

  !==================================================================================================
  subroutine cli_add_params(cli,params,group)
    implicit none
    type(Type_Command_Line_Interface)            , intent(inout) :: cli
    type(par_test_blk_nsi_cg_iss_oss_rk_params_t), intent(in)    :: params
    character(*)                                 , intent(in)    :: group
    ! Locals
    integer(ip) :: error
    
    ! Set Command Line Arguments
    ! Names
    call cli%add(group=trim(group),switch='--authors',switch_ab='-a',help='Print authors names',required=.false., &
         &       act='store_true',def=trim(params%default_authors),error=error)
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--path_working_directory',switch_ab='-pwd',help='Working directory path', &
         &       required=.false.,act='store',def=trim(params%default_pwd),error=error)
     if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--prefix',switch_ab='-prf',help='Test prefix', &
         &       required=.false.,act='store',def=trim(params%default_prefix),error=error)
     if(error/=0) then; check(.false.); end if

    ! Geometry
    call cli%add(group=trim(group),switch='--num_elements',switch_ab='-ne',help='Number of elements per direction', &
         &       required=.false.,act='store',nargs='*',def=trim(params%default_nelems),error=error)
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--num_parts',switch_ab='-np',help='Number of parts per direction', &
         &       required=.false.,act='store',nargs='*',def=trim(params%default_nparts),error=error)
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--num_sockets',switch_ab='-ns',help='Number of sockets per direction', &
         &       required=.false.,act='store',nargs='*',def=trim(params%default_nsockets),error=error)
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--discretization_type',help='Discretization type per direction (0=uniform, 1=cubic, 2=tanh, 3=imh+unif, 4:imh+tanh)', &
         &       required=.false.,act='store',nargs='*',def=trim(params%default_ndiscret),error=error)
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--domain_length',switch_ab='-dl',help='Domain length per direction', &
         &       required=.false.,act='store',nargs='*',def=trim(params%default_length),error=error)
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--periodicity',switch_ab='-peri',help='Periodicity per direction (0: No, 1: Yes)', &
         &       required=.false.,act='store',nargs='*',def=trim(params%default_periodic),error=error)
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--origin',switch_ab='-O',help='Origin coordinates', &
         &       required=.false.,act='store',nargs='*',def=trim(params%default_origin),error=error)
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--stretching',switch_ab='-st',help='Stretching parameter', &
         &       required=.false.,act='store',nargs='*',def=trim(params%default_stretch),error=error)
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--num_elements_bound_layer',switch_ab='-nb',help='Number of elements in the boundary layer per direction', &
         &       required=.false.,act='store',nargs='*',def=trim(params%default_nelbound),error=error)
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--size_bound_layer',switch_ab='-sb',help='Size of the boundary layer per direction', &
         &       required=.false.,act='store',nargs='*',def=trim(params%default_sizebound),error=error)
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--material_case',switch_ab='-mc',help='Material case', &
         &       required=.false.,act='store',def=trim(params%default_matercase),error=error)
    if(error/=0) then; check(.false.); end if

    ! Time integration
    call cli%add(group=trim(group),switch='--initial_time',switch_ab='-it',help='Initial time', &
         &       required=.false.,act='store',def=trim(params%default_initime),error=error)
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--final_time',switch_ab='-ft',help='Initial time', &
         &       required=.false.,act='store',def=trim(params%default_finaltime),error=error)
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--time_step',switch_ab='-dt',help='Time step size', &
         &       required=.false.,act='store',def=trim(params%default_timestep),error=error)
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--SRK_num_stages',switch_ab='-rks',help='Number of Runge-Kutta stages', &
         &       required=.false.,act='store',def=trim(params%default_rkstage),error=error)
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--SRK_temporal_order',switch_ab='-rko',help='Order of Runge-Kutta scheme', &
         &       required=.false.,act='store',def=trim(params%default_rkorder),error=error)
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--SRK_flag',switch_ab='-rkf',help='Flag for Runge-Kutta scheme', &
         &       required=.false.,act='store',def=trim(params%default_rkflag),error=error)
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--SRK_imex',switch_ab='-rkimex',help='0: explicit convection, 1: implicit convection', &
         &       required=.false.,act='store',choices='0,1',def=trim(params%default_rkimex),error=error)
    if(error/=0) then; check(.false.); end if

    ! Solution
    call cli%add(group=trim(group),switch='--analytical_v',switch_ab='-av',help='Analytical solution for velocity components', &
         &       required=.false.,act='store',nargs='*',def=trim(params%default_analytical_v),error=error)
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--analytical_p',switch_ab='-ap',help='Analytical solution for pressure', &
         &       required=.false.,act='store',def=trim(params%default_analytical_p),error=error)
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--analytical_t',switch_ab='-at',help='Analytical temporal component', &
         &       required=.false.,act='store',def=trim(params%default_analytical_t),error=error)
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--initial_condition',switch_ab='-insol',help='Give initial solution', &
         &       required=.false.,act='store',def=trim(params%default_initial_cond),error=error)
    if(error/=0) then; check(.false.); end if

    ! Problem
    call cli%add(group=trim(group),switch='--kfl_conv',switch_ab='-cnv',help='Convection flag', &
         &       required=.false.,act='store',choices='0,1',def=trim(params%default_kfl_conv),error=error)
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--kfl_skew',switch_ab='-skw',help='Skewness flag', &
         &       required=.false.,act='store',choices='0,1',def=trim(params%default_kfl_skew),error=error)
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--diffusion',switch_ab='-nu',help='Diffusion value', &
         &       required=.false.,act='store',def=trim(params%default_diffu),error=error)
    if(error/=0) then; check(.false.); end if

    ! Discretization
    call cli%add(group=trim(group),switch='--kfl_proj',switch_ab='-prj',help='Projection matrix flag (0: L2, 1: tau)', &
         &       required=.false.,act='store',choices='0,1',def=trim(params%default_kfl_proj),error=error)
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--k_tau_c',switch_ab='-ktc',help='Constant of tau_c', &
         &       required=.false.,act='store',def=trim(params%default_ktauc),error=error)
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--velocity_order',switch_ab='-vo',help='Velocity interpolation order', &
         &       required=.false.,act='store',def=trim(params%default_vel_order),error=error)
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--pressure_order',switch_ab='-po',help='Pressure interpolation order', &
         &       required=.false.,act='store',def=trim(params%default_pre_order),error=error)
    if(error/=0) then; check(.false.); end if

    ! Solver
    call cli%add(group=trim(group),switch='--solver_tolerance',switch_ab='-rtol',help='Solver tolerance', &
         &       required=.false.,act='store',def=trim(params%default_rtol),error=error)
    if(error/=0) then; check(.false.); end if

    ! Nonlinearity
    call cli%add(group=trim(group),switch='--nonlinear_tolerance',switch_ab='-nltol',help='Nonlinear tolerance', &
         &       required=.false.,act='store',def=trim(params%default_nltol),error=error)    
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--max_nonlinear_iter',switch_ab='-nlit',help='Maximum nonlinear iterations', &
         &       required=.false.,act='store',def=trim(params%default_max_iter),error=error)
    if(error/=0) then; check(.false.); end if

    ! Outputs
    call cli%add(group=trim(group),switch='--write_unkno_VTK',switch_ab='-wunk',help='Write solution in VTK format', &
         &       required=.false.,act='store',def=trim(params%default_write_unkno),error=error)
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--trace_unkno_VTK',switch_ab='-tunk',help='Intervals between two consecutive solution writes', &
         &       required=.false.,act='store',def=trim(params%default_trace_unkno),error=error)
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--write_dissipation',switch_ab='-wdis',help='Write dissipation outputs', &
         &       required=.false.,act='store',def=trim(params%default_write_dissipation),error=error)
    if(error/=0) then; check(.false.); end if
    call cli%add(group=trim(group),switch='--trace_dissipation',switch_ab='-tdis',help='Intervals between two consecutive dissipation writes', &
         &       required=.false.,act='store',def=trim(params%default_trace_dissipation),error=error)
    if(error/=0) then; check(.false.); end if

  end subroutine cli_add_params
  
end module command_line_parameters_names

program par_test_blk_nsi_cg_iss_oss_rk
  use serial_names
  use par_names
  use nsi_names
  use nsi_cg_iss_oss_names
  use norm_names
  use lib_vtk_io_interface_names
  use my_nonlinear_operator_names
  use Data_Type_Command_Line_Interface
  use command_line_parameters_names
  use ir_precision
  use dof_distribution_names
  implicit none
# include "debug.i90"

  ! Geometry & Integration Types
  type(uniform_mesh_descriptor_t)       :: gdata
  type(uniform_conditions_descriptor_t) :: bdata
  type(reference_element_t)             :: geo_reference_element
  type(par_triangulation_t)             :: p_trian
  type(par_conditions_t)                :: p_cond
  type(dof_descriptor_t)                :: dof_descriptor
  type(par_fe_space_t)                  :: p_fe_space  
  type(par_block_graph_t)               :: p_blk_graph
  type(block_dof_distribution_t)        :: blk_dof_dist
  type(par_scalar_t)                    :: enorm_u, enorm_p, dummy

  ! Problem types
  type(nsi_problem_t)                                 :: myprob
  type(nsi_cg_iss_oss_discrete_t)                     :: mydisc
  type(nsi_cg_iss_oss_rk_momentum_t)         , target :: cg_iss_oss_rk_momentum
  type(nsi_cg_iss_oss_rk_momentum_rhs_t)     , target :: cg_iss_oss_rk_momentum_rhs
  type(nsi_cg_iss_oss_rk_pressure_t)         , target :: cg_iss_oss_rk_pressure
  type(nsi_cg_iss_oss_rk_momentum_update_t)  , target :: cg_iss_oss_rk_momentum_update
  type(nsi_cg_iss_oss_rk_projection_update_t), target :: cg_iss_oss_rk_projection_update
  type(nsi_cg_iss_oss_lapla_p_t)             , target :: lapla_p_integration
  type(nsi_cg_iss_oss_massu_t)               , target :: mass_u_integration
  type(nsi_cg_iss_oss_massx_t)               , target :: mass_x_integration
  type(nsi_cg_iss_oss_dissipation_t)         , target :: dissipation_integration
  type(error_norm_t)                         , target :: error_compute
  type(rungekutta_integrator_t)              , target :: rkinteg
  type(time_integration_t)                   , target :: tinteg
  type(discrete_integration_pointer_t)                :: approx(1)
  type(momentum_operator_t)                           :: momentum_operator
  type(pressure_operator_t)                           :: pressure_operator
  type(momentum_update_operator_t)                    :: momentum_update_operator
  type(projection_update_operator_t)                  :: projection_update_operator

  ! Solver types
  type(solver_control_t)        :: sctrl
  type(par_context_t)           :: w_context
  type(par_context_t)           :: p_context
  type(par_context_t)           :: q_context
  type(par_context_t)           :: b_context
  type(par_environment_t)       :: p_env
  type(par_timer_t)             :: p_timer

  ! Postproces types
  type(vtk_t) :: fevtk

  ! Integers
  integer(ip) :: num_levels,me,np
  integer(ip) :: gtype(3) = (/ csr, csr, csr /)
  integer(ip) :: ibloc,jbloc,istat,i
  integer(ip) :: num_approximations = 1
  integer(ip) :: setterms(6,2),settable(3)
  integer(ip) :: order_v,order_p
  integer(ip) :: nstage,rk_order,rk_flag,rk_imex
  integer(ip) :: analytical_vx,analytical_vy,analytical_vz,analytical_p,analytical_t
  integer(ip) :: lunio_dis
  integer(ip) :: max_iter

  ! Reals
  real(rp) :: initime,finaltime,dt
  real(rp) :: dissipations(10)
  real(rp) :: trace_unkno,trace_dissipation
  real(rp) :: nltol

  ! Logicals
  logical  :: initial_condition,write_unkno_VTK,write_dissipation

  ! Parameters
  integer(ip), parameter :: velocity=1, pressure=2
  integer(ip), parameter :: nonlinear=0, linear=1
  integer(ip), parameter :: explicit=0, implicit=1

  ! Allocatable
  integer(ip), allocatable :: id_parts(:)
  integer(ip), allocatable :: num_parts(:)
  integer(ip), allocatable :: continuity(:,:)
  integer(ip), allocatable :: order(:,:)
  integer(ip), allocatable :: material(:)
  integer(ip), allocatable :: problem(:)
  integer(ip), allocatable :: vars_block(:)
  integer(ip), allocatable :: dof_coupling(:,:)
  integer(ip), allocatable :: analytical_v(:)

  ! Arguments
  type(Type_Command_Line_Interface):: cli 
  character(len=:), allocatable :: group
  character(len=256)            :: dir_path_out,prefix

  call meminit

  ! Read parameters from command-line
  call read_flap_cli_par_test_blk_nsi_cg_iss_oss_rk(cli)
  call cli%parse(error=istat)
  if(cli%run_command('analytical')) then
     group = 'analytical'
  elseif(cli%run_command('Taylor_Green')) then
     group = 'Taylor_Green'
  end if

  ! Get names
  call cli%get(group=trim(group),switch='-pwd',val=dir_path_out,error=istat); if(istat/=0) then; check(.false.); end if
  call cli%get(group=trim(group),switch='-prf',val=prefix,error=istat); if(istat/=0) then; check(.false.); end if

  ! Generate geometry data 
  call uniform_mesh_descriptor_create(gdata,cli,group)

  ! Generate boundary data
  call uniform_conditions_descriptor_create(2*gdata%ndime+1,2*gdata%ndime+1,gdata%ndime,bdata)
  bdata%poin%code(gdata%ndime+1,1:2**gdata%ndime-1) = 0
  bdata%line%code(gdata%ndime+1,:) = 0
  bdata%surf%code(gdata%ndime+1,:) = 0
  bdata%poin%code(gdata%ndime+1,2**gdata%ndime) = 1
  bdata%poin%valu(gdata%ndime+1,2**gdata%ndime) = 0.0_rp
  bdata%poin%valu(1:gdata%ndime,:) = 1.0_rp
  bdata%line%valu(1:gdata%ndime,:) = 1.0_rp
  bdata%poin%code(gdata%ndime+2:2*gdata%ndime+1,:) = 0
  bdata%line%code(gdata%ndime+2:2*gdata%ndime+1,:) = 0
  bdata%surf%code(gdata%ndime+2:2*gdata%ndime+1,:) = 0

  ! Generate element geometrical fixed info
  call reference_element_create(geo_reference_element,Q_type_id,1,gdata%ndime)

  ! Set levels
  num_levels = 2
  call memalloc(num_levels, id_parts , __FILE__, __LINE__)
  call memalloc(num_levels, num_parts, __FILE__, __LINE__)
  num_parts = (/gdata%nparts, 1/)
  id_parts = (/w_context%iam+1, 1/)

  ! Start parallel execution
  call par_context_create (w_context)

  ! Create p_context and q_context splitting w_context
  if(w_context%iam < num_parts(1)) then
     call par_context_create(1,p_context,q_context,w_context)
  else
     call par_context_create(2,q_context,p_context,w_context)
  end if
  check((p_context%iam>=0 .and. q_context%iam<0) .or. (p_context%iam<0 .and. q_context%iam>= 0))

  ! Create b_context as an intercommunicator among p_context <=> q_context 
  call par_context_create(w_context,p_context,q_context,b_context)

  ! Create parallel environment
  call par_environment_create(p_env,w_context,p_context,q_context,b_context,num_levels,id_parts,num_parts)
  call p_env%info(me,np)

  ! Generate par triangulation
  call par_generate_uniform_triangulation(p_env,gdata,bdata,geo_reference_element,p_trian,p_cond,material)
  if(gdata%ndime==2.and.gdata%isper(1)==1.and.gdata%isper(2)==1) then
     if(me==0) p_cond%f_conditions%code(gdata%ndime+1,1) = 1
  elseif(gdata%ndime==3.and.gdata%isper(1)==1.and.gdata%isper(2)==1.and.gdata%isper(3)==1) then
     if(me==0) p_cond%f_conditions%code(gdata%ndime+1,1) = 1
  end if

  ! Define Runge-Kutta method
  call cli%get(group=trim(group),switch='-rks',val=nstage,error=istat); if(istat/=0) then; check(.false.); end if
  call cli%get(group=trim(group),switch='-rko',val=rk_order,error=istat); if(istat/=0) then; check(.false.); end if
  call cli%get(group=trim(group),switch='-rkf',val=rk_flag,error=istat); if(istat/=0) then; check(.false.); end if
  call cli%get(group=trim(group),switch='-rkimex',val=rk_imex,error=istat); if(istat/=0) then; check(.false.); end if
  settable      = (/nstage,rk_order,rk_flag/)
  if(rk_imex==0) then ! Explicit convection
     setterms(1,:) = (/update_constant,implicit/)   ! Diffusion
     setterms(2,:) = (/update_transient,explicit/)  ! Convection
     setterms(3,:) = (/update_transient,explicit/)  ! Pressure Gradient
     setterms(4,:) = (/update_transient,explicit/)  ! OSS_vu
     setterms(5,:) = (/update_transient,explicit/)  ! OSS_vx
     setterms(6,:) = (/update_transient,explicit/)  ! Force
  elseif(rk_imex==1) then ! Implicit convection
     setterms(1,:) = (/update_constant,implicit/)   ! Diffusion
     setterms(2,:) = (/update_nonlinear,implicit/)  ! Convection
     setterms(3,:) = (/update_transient,explicit/)  ! Pressure Gradient
     setterms(4,:) = (/update_nonlinear,implicit/)  ! OSS_vu
     setterms(5,:) = (/update_nonlinear,implicit/)  ! OSS_vx
     setterms(6,:) = (/update_transient,implicit/)  ! Force
  end if
  call rkinteg%create(setterms,settable)

  ! Create problem
  call myprob%create(gdata%ndime)
  call mydisc%create(myprob)
  call mydisc%vars_block(myprob,vars_block)
  call mydisc%dof_coupling(myprob,dof_coupling)
  call cg_iss_oss_rk_momentum%create(myprob,mydisc)
  call cg_iss_oss_rk_momentum_rhs%create(myprob,mydisc)
  call cg_iss_oss_rk_pressure%create(myprob,mydisc)
  call cg_iss_oss_rk_momentum_update%create(myprob,mydisc)
  call cg_iss_oss_rk_projection_update%create(myprob,mydisc)
  call lapla_p_integration%create(myprob,mydisc)
  call mass_u_integration%create(myprob,mydisc)
  call mass_x_integration%create(myprob,mydisc)
  call dissipation_integration%create(myprob,mydisc)
  cg_iss_oss_rk_momentum%rkinteg => rkinteg
  cg_iss_oss_rk_momentum_rhs%rkinteg => rkinteg
  cg_iss_oss_rk_pressure%tinteg  => tinteg
  cg_iss_oss_rk_momentum_update%rkinteg => rkinteg
  dissipation_integration%tinteg  => tinteg
  ! Temporal Command Line variables
  call cli%get(group=trim(group),switch='-dt',val=dt,error=istat); if(istat/=0) then; check(.false.); end if  
  call cli%get(group=trim(group),switch='-it',val=initime,error=istat); if(istat/=0) then; check(.false.); end if  
  call cli%get(group=trim(group),switch='-ft',val=finaltime,error=istat); if(istat/=0) then; check(.false.); end if  
  rkinteg%dtinv   = 1.0_rp/dt
  tinteg%dtinv    = 1.0_rp/dt
  rkinteg%ftime   = finaltime
  rkinteg%itime   = initime
  ! Problem Command Line variables
  call cli%get(group=trim(group),switch='-cnv',val=myprob%kfl_conv,error=istat); if(istat/=0) then; check(.false.); end if  
  call cli%get(group=trim(group),switch='-skw',val=myprob%kfl_skew,error=istat); if(istat/=0) then; check(.false.); end if  
  call cli%get(group=trim(group),switch='-nu',val=myprob%diffu,error=istat); if(istat/=0) then; check(.false.); end if
  ! Discretization Command Line variables
  call cli%get(group=trim(group),switch='-prj',val=mydisc%kfl_proj,error=istat); if(istat/=0) then; check(.false.); end if  
  call cli%get(group=trim(group),switch='-ktc',val=mydisc%ktauc,error=istat); if(istat/=0) then; check(.false.); end if  
  call cli%get(group=trim(group),switch='-vo',val=order_v,error=istat); if(istat/=0) then; check(.false.); end if  
  call cli%get(group=trim(group),switch='-po',val=order_p,error=istat); if(istat/=0) then; check(.false.); end if  

  ! Create dof_descriptor
  call dof_descriptor%create(3,1,mydisc%nvars,vars_block,dof_coupling)
  call dof_descriptor%set_problem(1,mydisc)

  ! Allocate auxiliar elemental arrays
  call memalloc(p_trian%f_trian%num_elems,dof_descriptor%nvars_global,continuity, __FILE__,__LINE__)
  call memalloc(p_trian%f_trian%num_elems,dof_descriptor%nvars_global,order,__FILE__,__LINE__)
  call memalloc(p_trian%f_trian%num_elems,problem,__FILE__,__LINE__)
  continuity             = 1
  order                  = order_v
  order(:,gdata%ndime+1) = order_p
  problem                = 1

  ! Create par_fe_space
  call par_fe_space_create(p_trian,dof_descriptor,p_fe_space,problem,p_cond,continuity,order,material, &
       &                   time_steps_to_store=3+rkinteg%rk_table_implicit%stage,                      &
       &                   hierarchical_basis=.false.,static_condensation=.false.,num_continuity=1)

  ! Initialize VTK output
  call cli%get(group=trim(group),switch='-wunk',val=write_unkno_VTK,error=istat); if(istat/=0) then; check(.false.); end if
  call cli%get(group=trim(group),switch='-tunk',val=trace_unkno,error=istat); if(istat/=0) then; check(.false.); end if
  if(write_unkno_VTK) then
     call fevtk%initialize(p_trian%f_trian,p_fe_space%fe_space,myprob,p_env,dir_path_out,prefix, &
          &                nparts=gdata%nparts,linear_order=.true.)
  end if

  ! Initialize dissipation output
  call cli%get(group=trim(group),switch='-wdis',val=write_dissipation,error=istat); if(istat/=0) then; check(.false.); end if
  call cli%get(group=trim(group),switch='-tdis',val=trace_dissipation,error=istat); if(istat/=0) then; check(.false.); end if
  call dummy%create(p_env)
  if(me==0.and.write_dissipation) then
     lunio_dis = io_open(trim(dir_path_out)//trim(prefix)//'.dis',position='append')
  end if

  ! Create dof info
  call par_create_distributed_dof_info(dof_descriptor,p_trian,p_fe_space,blk_dof_dist,p_blk_graph,gtype)  

  ! Assign analytical solution
  call cli%get_varying(group=trim(group),switch='-av',val=analytical_v,error=istat); if(istat/=0) then; check(.false.); end if   
  call cli%get(group=trim(group),switch='-ap',val=analytical_p,error=istat); if(istat/=0) then; check(.false.); end if  
  call cli%get(group=trim(group),switch='-at',val=analytical_t,error=istat); if(istat/=0) then; check(.false.); end if   
  analytical_vx=analytical_v(1)
  analytical_vy=analytical_v(2)
  if(gdata%ndime==2) then
     call par_fe_space_set_analytical_code(p_fe_space,(/analytical_vx,analytical_vy,analytical_p,0,0/), &
          &                            (/analytical_t,analytical_t,0,0,0/))
  else
     analytical_vz=analytical_v(3)
     call par_fe_space_set_analytical_code(p_fe_space,(/analytical_vx,analytical_vy,analytical_vz,analytical_p,0,0,0/), &
          &                            (/analytical_t,analytical_t,analytical_t,0,0,0,0/))
  end if

  ! Create picard nonlinear operators
  call cli%get(group=trim(group),switch='-nltol',val=nltol,error=istat); if(istat/=0) then; check(.false.); end if  
  call cli%get(group=trim(group),switch='-nlit',val=max_iter,error=istat); if(istat/=0) then; check(.false.); end if    
  call momentum_operator%build(p_blk_graph,p_env,gdata%ndime,max_iter,nltol)
  call pressure_operator%build(p_blk_graph,p_env,gdata%ndime)
  call momentum_update_operator%build(p_blk_graph)
  call projection_update_operator%build(p_blk_graph)

  ! Solver control parameters 
  sctrl%method  = lgmres
  sctrl%trace   = 100
  sctrl%itmax   = 800
  sctrl%dkrymax = 800
  sctrl%stopc   = res_nrmgiven_res_nrmgiven
  sctrl%orto    = icgs
  call cli%get(group=trim(group),switch='-rtol',val=sctrl%rtol,error=istat); if(istat/=0) then; check(.false.); end if  
  sctrl%track_conv_his = .false.

  ! Assign initial condition
  call cli%get(group=trim(group),switch='-insol',val=initial_condition,error=istat); if(istat/=0) then; check(.false.); end if  
  if(initial_condition) then
     call par_update_analytical_initial((/(i,i=1,gdata%ndime+1)/),rkinteg%itime,p_fe_space)
  end if

  ! Do time steps
  call par_timer_create(p_timer,'TEMPORAL_LOOP', w_context%icontxt)
  call par_timer_init(p_timer)
  call par_timer_start(p_timer)   
  call do_time_steps_rk_nsi(rkinteg,sctrl,1.0e-7_rp,100,p_env,gdata,p_fe_space,momentum_operator, &
       &                    cg_iss_oss_rk_momentum,pressure_operator,cg_iss_oss_rk_pressure,      &
       &                    momentum_update_operator,cg_iss_oss_rk_momentum_update,               &
       &                    projection_update_operator,cg_iss_oss_rk_projection_update,           &
       &                    cg_iss_oss_rk_momentum_rhs,tinteg,fevtk,write_unkno_VTK,trace_unkno,  &
       &                    write_dissipation,dissipation_integration,trace_dissipation,dummy,    &
       &                    lunio_dis)
  call par_timer_stop(p_timer)   
  call par_timer_report(p_timer) 

  ! Print PVD file
  if(write_unkno_VTK) istat = fevtk%write_PVD() 

  ! Compute error norm
  call error_compute%create(myprob,mydisc)
  approx(1)%p => error_compute
  error_compute%ctime = rkinteg%ctime
  error_compute%unknown_id = velocity
  call enorm_u%create(p_env)
  call enorm_u%init()
  call par_volume_integral(approx,p_fe_space,enorm_u)
  call enorm_u%reduce()
  error_compute%unknown_id = pressure
  call enorm_p%create(p_env)
  call enorm_p%init()
  call par_volume_integral(approx,p_fe_space,enorm_p)
  call enorm_p%reduce()
  if(me==0) then
     write(*,*) 'Velocity error norm: ', sqrt(enorm_u%get())
     write(*,*) 'Pressure error norm: ', sqrt(enorm_p%get()) 
  end if

  ! Deallocate
  deallocate(analytical_v)
  call memfree(continuity,__FILE__,__LINE__)
  call memfree(order,__FILE__,__LINE__)
  call memfree(material,__FILE__,__LINE__)
  call memfree(problem,__FILE__,__LINE__)
  call memfree(vars_block,__FILE__,__LINE__)
  call memfree(dof_coupling,__FILE__,__LINE__)
  call memfree(id_parts , __FILE__, __LINE__)
  call memfree(num_parts, __FILE__, __LINE__)
  call momentum_operator%free()
  call pressure_operator%free()
  call momentum_update_operator%free()
  call projection_update_operator%free()
  if(write_unkno_VTK) call fevtk%free
  if(write_dissipation.and.me==0) call io_close(lunio_dis)
  call p_blk_graph%free
  call blk_dof_dist%free
  call par_fe_space_free(p_fe_space) 
  call myprob%free
  call mydisc%free
  call error_compute%free
  call cg_iss_oss_rk_momentum%free
  call cg_iss_oss_rk_momentum_rhs%free
  call cg_iss_oss_rk_pressure%free
  call cg_iss_oss_rk_momentum_update%free
  call cg_iss_oss_rk_projection_update%free
  call lapla_p_integration%free
  call mass_u_integration%free
  call mass_x_integration%free
  call dissipation_integration%free
  call rkinteg%free
  call dof_descriptor_free(dof_descriptor)
  call par_triangulation_free(p_trian)
  call par_conditions_free (p_cond)
  call reference_element_free(geo_reference_element)
  call uniform_conditions_descriptor_free(bdata)
  call par_environment_free (p_env)
  call par_context_free ( b_context, .false. )
  call par_context_free ( p_context, .false. )
  call par_context_free ( q_context, .false. )
  call par_context_free ( w_context )

  call memstatus

contains

  !==================================================================================================
  subroutine read_flap_cli_par_test_blk_nsi_cg_iss_oss_rk(cli)
    implicit none
    type(Type_Command_Line_Interface), intent(out) :: cli
    ! Locals
    type(par_test_blk_nsi_cg_iss_oss_rk_params_t) :: analytical_params,taylor_green_params
    logical     :: authors_print
    integer(ip) :: error
    
    authors_print = .false.
    
    ! Initialize Command Line Interface
    call cli%init(progname    = 'par_test_blk_nsi_cg_iss_oss_rk', &
         &        version     = '',                               &
         &        authors     = '',                               &
         &        license     = '',                               &
         &        description = 'Parallel FEMPAR driver to solve transient Navier-Stokes problems using Continuous-Galerkin Inf-Sup stable elements with Orthogonal Subscales stabilization for the spatial discretization and Segregated Runge-Kutta schemes for the temporal discretization.', &
         &        examples    = ['par_test_blk_nsi_cg_iss_oss_rk -h         ', &
         &                       'par_test_blk_nsi_cg_iss_oss_rk analytical ', &
         &                       'par_test_blk_nsi_cg_iss_oss_rk Taylor_Green'])
    
    ! Set Command Line Arguments Groups, i.e. commands
    call cli%add_group(group='analytical',description='solve a problem with an analytical solution')
    call cli%add_group(group='Taylor_Green',description='solve the Taylor-Green Vortex flow problem')
    
    ! Set Command Line Arguments for each group
    call set_default_params(analytical_params)
    call set_default_params_analytical(analytical_params)
    call cli_add_params(cli,analytical_params,'analytical')
    call set_default_params(taylor_green_params)
    call set_default_params_taylor_green(taylor_green_params)
    call cli_add_params(cli,taylor_green_params,'Taylor_Green')
    
  end subroutine read_flap_cli_par_test_blk_nsi_cg_iss_oss_rk
  
  !==================================================================================================
  subroutine do_time_steps_rk_nsi(rkinteg,sctrl,sttol,maxst,p_env,gdata,p_fe_space,momentum_operator, &
       &                          momentum_integration,pressure_operator,pressure_integration,        &
       &                          momentum_update_operator,momentum_update_integration,               &
       &                          projection_update_operator,projection_update_integration,           &
       &                          momentum_rhs_integration,tinteg,fevtk,write_unkno_VTK,trace_unkno,  &
       &                          write_dissipation,dissipation_integration,trace_dissipation,dummy,  &
       &                          lunio_dis)
    implicit none
    type(rungekutta_integrator_t)        , intent(inout) :: rkinteg
    type(solver_control_t)               , intent(inout) :: sctrl
    real(rp)                             , intent(in)    :: sttol
    integer(ip)                          , intent(in)    :: maxst
    class(abstract_environment_t)        , intent(in)    :: p_env
    type(uniform_mesh_descriptor_t)      , intent(in)    :: gdata
    type(par_fe_space_t)                 , intent(inout) :: p_fe_space
    class(momentum_operator_t)           , intent(inout) :: momentum_operator
    class(pressure_operator_t)           , intent(inout) :: pressure_operator
    class(momentum_update_operator_t)    , intent(inout) :: momentum_update_operator
    class(projection_update_operator_t)  , intent(inout) :: projection_update_operator
    class(discrete_integration_t), target, intent(inout) :: momentum_integration
    class(discrete_integration_t), target, intent(inout) :: pressure_integration
    class(discrete_integration_t), target, intent(inout) :: momentum_update_integration
    class(discrete_integration_t), target, intent(inout) :: projection_update_integration
    class(discrete_integration_t), target, intent(inout) :: momentum_rhs_integration
    class(time_integration_t)            , intent(inout) :: tinteg
    type(vtk_t)                          , intent(inout) :: fevtk
    logical                              , intent(in)    :: write_unkno_VTK,write_dissipation
    real(rp)                             , intent(in)    :: trace_unkno,trace_dissipation
    class(discrete_integration_t), target, intent(inout) :: dissipation_integration
    type(par_scalar_t)                   , intent(inout) :: dummy
    integer(ip)                          , intent(in)    :: lunio_dis
    ! Locals
    type(discrete_integration_pointer_t) :: approx(1)
    integer(ip) :: istage,nstage,istep,ielem,me,np
    real(rp)    :: rtime,ctime,prevtime,prevtime_vtk,prevtime_dis
    logical     :: status

    ! Initialize time steps
    rtime = 1.0_rp
    istep = 1
    prevtime_vtk = 0.0_rp
    prevtime_dis = 0.0_rp

    ! Get process info
    call p_env%info(me,np)

    ! Update steady boundary conditions
    call par_update_strong_dirichlet_bcond(p_fe_space,p_cond)
          
    ! Update analytical/time dependent boundary conditions
    call par_update_analytical_bcond((/(i,i=1,gdata%ndime+1)/),rkinteg%ctime,p_fe_space)

    ! Initialize vector
    call par_update_initialize(momentum_operator%x_sol,p_fe_space)
    call par_update_initialize(pressure_operator%x_sol,p_fe_space)
    call par_update_initialize(momentum_update_operator%x_sol,p_fe_space)
    call par_update_initialize(projection_update_operator%x_sol,p_fe_space)

    !**************************** NSI specific tasks *****************************************!
    ! Preconditioner matrices
    approx(1)%p => lapla_p_integration
    call par_volume_integral(approx,p_fe_space,pressure_operator%p_lapla_p_matrix)
    approx(1)%p => mass_u_integration
    call par_volume_integral(approx,p_fe_space,pressure_operator%p_mass_u_matrix)
    call par_volume_integral(approx,p_fe_space,momentum_update_operator%p_mass_u_matrix)
    approx(1)%p => mass_x_integration
    call par_volume_integral(approx,p_fe_space,momentum_operator%p_mass_x_matrix)
    if(p_env%am_i_fine_task()) then
       projection_update_operator%p_mass_x_matrix%blocks(3,3)%p_p_matrix%f_matrix%a = &
            & momentum_operator%p_mass_x_matrix%blocks(3,3)%p_p_matrix%f_matrix%a
    end if
    !*****************************************************************************************!

    ! Compute constant operators
    ! Momentum
    approx(1)%p => momentum_integration
    momentum_integration%integration_stage = update_constant
    call momentum_operator%fill_constant(approx,p_fe_space) 
    ! Momentum update
    approx(1)%p => momentum_update_integration
    momentum_update_integration%integration_stage = update_constant
    call momentum_update_operator%fill_constant(approx,p_fe_space)
    !**************************** NSI specific tasks *****************************************!
    ! Pressure
    approx(1)%p => pressure_integration
    pressure_integration%integration_stage = update_constant
    call pressure_operator%fill_constant(approx,p_fe_space)
    ! Projection update
    approx(1)%p => projection_update_integration
    projection_update_integration%integration_stage = update_constant
    call projection_update_operator%fill_constant(approx,p_fe_space)
    !*****************************************************************************************!

    ! Time steps loop
    step: do while (rkinteg%ctime<rkinteg%ftime.and.rtime>sttol.and.istep<=maxst)

       ! Current time
       prevtime = rkinteg%ctime
       rkinteg%ctime = rkinteg%ctime + 1.0_rp/rkinteg%dtinv
       if(me==0) then
          write(*,'(a)') '============================================================'
          write(*,'(a15,i8,a,a20,e15.8)') 'Time step:     ',istep,', ','Current time: ',rkinteg%ctime
       end if
       
       ! Loop over stages
       stage: do istage=1,rkinteg%rk_table(1)%p%stage

          ! Set current time
          ctime = prevtime + 1.0_rp/rkinteg%dtinv*rkinteg%rk_table(1)%p%c(istage)
          rkinteg%istage = istage
          if(me==0) then
             write(*,'(a)') '------------------------------------------------------------'
             write(*,'(a20,i3,a,a21,e15.8)') 'Runge-Kutta stage:  ',istage,',','Current time: ',ctime
          end if

          ! Set momentum operator scalars
          momentum_operator%dt  = rkinteg%dtinv
          momentum_operator%aii = rkinteg%rk_table_implicit%A(istage,istage)

          ! Update operators with scalars
          momentum_operator%scal_Am = momentum_operator%dt * momentum_operator%Am%get_block(1,1)
          momentum_operator%scal_Al = momentum_operator%aii * momentum_operator%Al%get_block(1,1)
          momentum_operator%scal_Anl = momentum_operator%aii * momentum_operator%Anl%get_block(1,1)
          momentum_operator%A_uu = momentum_operator%scal_Am + momentum_operator%scal_Al + momentum_operator%scal_Anl
          call momentum_operator%A_momentum%set_block(1,1,momentum_operator%A_uu)
          call momentum_operator%A_momentum%set_block(1,2, &
               & momentum_operator%aii*momentum_operator%p_block_matrix_nl%get_block(1,3))
          call momentum_operator%A_momentum%set_block(2,1, &
               & momentum_operator%aii*momentum_operator%p_block_matrix_nl%get_block(3,1))
          call momentum_operator%A_momentum%set_block(2,2, &
               & momentum_operator%aii*momentum_operator%p_block_matrix_nl%get_block(3,3))  
          call momentum_operator%M_momentum%set_block(2,1, &
               & momentum_operator%aii*momentum_operator%p_block_matrix_nl%get_block(3,1))
         
          ! Skip 1st stage
          if(istage==1.and.momentum_operator%aii==0) then
             if(me==0) then
                write(*,*) 'First stage skipped for momentum equation.'
             end if
             ! (U_n --> U_1)
             approx(1)%p => momentum_integration
             call par_update_nonlinear_solution(p_fe_space,approx(1)%p%working_vars,3,3+1)
             ! (U_0 --> U_k)
             call par_update_nonlinear_solution(p_fe_space,approx(1)%p%working_vars,1,2)
          else

             ! Update analytical/time dependent boundary conditions
             call par_update_analytical_bcond((/(i,i=1,gdata%ndime+1)/),ctime,p_fe_space)

             ! Compute momentum transient operators
             !********************************** This is DIRTY *************************************!
             approx(1)%p => momentum_integration
             momentum_integration%integration_stage = update_constant  ! Liniar operator
             call momentum_operator%fill_constant(approx,p_fe_space) 
             !**************************************************************************************! 
             approx(1)%p => mass_u_integration
             momentum_integration%integration_stage = update_transient ! Mass matrix
             call momentum_operator%fill_transient(approx,p_fe_space) 

             ! Compute momentum operand (RHS)
             approx(1)%p => momentum_rhs_integration
             call momentum_operator%p_block_vector%init(0.0_rp)
             call par_volume_integral(approx,p_fe_space,momentum_operator%p_block_vector)

             !********************************** This is DIRTY *************************************!
             if(p_env%am_i_fine_task()) then
                momentum_operator%p_block_u_matrix%blocks(1,1)%p_p_matrix%f_matrix%a = &
                     & momentum_operator%dt*momentum_operator%p_block_matrix_m%blocks(1,1)%p_p_matrix%f_matrix%a + &
                     & momentum_operator%aii*momentum_operator%p_block_matrix_l%blocks(1,1)%p_p_matrix%f_matrix%a
             end if
             call momentum_operator%M%fill_values(update_transient) ! Recompute M with updated matrix values 
             !**************************************************************************************! 

             ! Momentum equation solution
             approx(1)%p => momentum_integration
             call momentum_operator%apply(sctrl,p_env,approx,p_fe_space)

             ! Free Transient
             call momentum_operator%M%free_values(update_transient)

             ! Store unkno to istage position (U --> U_i)
             call par_update_nonlinear_solution(p_fe_space,approx(1)%p%working_vars,1,3+istage)   

          end if

          !**************************** NSI specific tasks *****************************************!
          ! Update boundary conditions (velocity derivative)
          call par_update_analytical_bcond((/(i,i=1,gdata%ndime)/),ctime,p_fe_space,2)

          ! Compute pressure transient operators
          approx(1)%p => pressure_integration
          pressure_integration%integration_stage = update_transient
          tinteg%ctime = ctime
          call pressure_operator%fill_transient(approx,p_fe_space) 

          ! Pressure equation solution
          approx(1)%p => pressure_integration
          call pressure_operator%apply(sctrl,p_env,approx,p_fe_space)

          ! Store pressure (P --> P_i)
          call par_update_nonlinear_solution(p_fe_space,pressure_integration%working_vars,1,3+istage) 

          ! Restore velocity (U_i --> U)
          call par_update_nonlinear_solution(p_fe_space,momentum_integration%working_vars,3+istage,1) 
          !*****************************************************************************************! 

       end do stage

       ! Set current time
       ctime = rkinteg%ctime
       if(me==0) then
          write(*,'(a)') '------------------------------------------------------------'
          write(*,'(a24,a21,e15.8)') 'Runge-Kutta update      ','Current time: ',ctime
       end if

       ! Update analytical/time dependent boundary conditions
       call par_update_analytical_bcond((/(i,i=1,gdata%ndime+1)/),ctime,p_fe_space)           

       ! Compute momentum update transient operators
       approx(1)%p => momentum_update_integration
       momentum_update_integration%integration_stage = update_transient
       call momentum_update_operator%fill_transient(approx,p_fe_space) 
       call momentum_update_operator%M%fill_values(update_transient) ! Recompute M with updated matrix values 
       
       ! Momentum update equation solution
       approx(1)%p => momentum_update_integration
       call momentum_update_operator%apply(sctrl,p_env,approx,p_fe_space)

       ! Free preconditioner 
       call momentum_update_operator%M%free_values(update_transient)

       ! Store unkno to previous step position (Un+1 --> U_n)
       call par_update_nonlinear_solution(p_fe_space,approx(1)%p%working_vars,1,3)    
       
       !**************************** NSI specific tasks *****************************************!   
       ! Compute projection transient operators
       approx(1)%p => projection_update_integration
       projection_update_integration%integration_stage = update_transient
       call projection_update_operator%fill_transient(approx,p_fe_space) 
       call projection_update_operator%M%fill_values(update_transient) ! Recompute M with updated matrix values 

       ! Projection update equation solution
       approx(1)%p => projection_update_integration
       call projection_update_operator%apply(sctrl,p_env,approx,p_fe_space)

       ! Free preconditioner 
       call projection_update_operator%M%free_values(update_transient)

       ! Store projection to previous step position (Un+1 --> U_n)
       call par_update_nonlinear_solution(p_fe_space,approx(1)%p%working_vars,1,3)   

       ! Update boundary conditions (velocity derivative)
       call par_update_analytical_bcond((/(i,i=1,gdata%ndime)/),ctime,p_fe_space,2)

       ! Compute pressure transient operators
       approx(1)%p => pressure_integration
       pressure_integration%integration_stage = update_transient
       tinteg%ctime = ctime
       call pressure_operator%fill_transient(approx,p_fe_space) 
       
       ! Pressure equation solution
       approx(1)%p => pressure_integration
       call pressure_operator%apply(sctrl,p_env,approx,p_fe_space)
       
       ! Restore velocity (U_n+1 --> U)
       call par_update_nonlinear_solution(p_fe_space,momentum_integration%working_vars,3,1) 

       ! Update analytical/time dependent boundary conditions
       call par_update_analytical_bcond((/(i,i=1,gdata%ndime+1)/),ctime,p_fe_space)
       !*****************************************************************************************!  

       ! Check steady state

       ! Update counter
       istep = istep + 1

       ! Postprocess
       ! ===========
       ! Compute dissipations
       status = ((ctime-prevtime_dis).ge.trace_dissipation)
       if(write_dissipation.and.status) then
          prevtime_dis = ctime
          approx(1)%p => dissipation_integration
          call par_volume_integral(approx,p_fe_space,dummy)
          select type(dissipation_integration)
          class is(nsi_cg_iss_oss_dissipation_t) 
             dissipations = dissipation_integration%values
          end select
          if(p_env%am_i_fine_task()) then
             select type(p_env)
             class is(par_environment_t) 
                call psb_sum(p_env%p_context%icontxt,dissipations)
             end select
          end if
          if(me==0) then
             ! Volume average
             dissipations(:) = dissipations(:)/(gdata%xleng*gdata%yleng*gdata%zleng)
             ! Write on a file
             write(lunio_dis,*) 'TIME__:',ctime
             write(lunio_dis,*) 'DISNUM:',dissipations(1)                                                
             write(lunio_dis,*) 'DISVIS:',dissipations(2)                                                
             write(lunio_dis,*) 'DISCON:',dissipations(3)                                                
             write(lunio_dis,*) 'DISTIM:',dissipations(4)                                                
             write(lunio_dis,*) 'DISFEX:',dissipations(5)                                         
             write(lunio_dis,*) 'DISTOT:',sum(dissipations(2:5))
             write(lunio_dis,*) 'ENERGY:',dissipations(6)                                                
             write(lunio_dis,*) 'ENSTRO:',dissipations(7)
             write(lunio_dis,*) 'SKEWNE:',dissipations(8)/sqrt(dissipations(9)**3)                       
             write(lunio_dis,*) 'DISTUH:',dissipations(2)+dissipations(3)+dissipations(4)+dissipations(5)
          end if
       end if

       ! Print solution to VTK file
       status = ((ctime-prevtime_vtk).ge.trace_unkno)
       if(write_unkno_VTK.and.status) then 
          prevtime_vtk = ctime
          select type(p_env)
          class is(par_environment_t) 
             istat = fevtk%write_VTK(n_part=p_env%p_context%iam,t_step=ctime)
          end select
          istat = fevtk%write_PVTK(t_step=ctime)
       end if

    end do step

    ! Free constant
    call pressure_operator%M%free_values(update_constant)

  end subroutine do_time_steps_rk_nsi
  
end program par_test_blk_nsi_cg_iss_oss_rk
