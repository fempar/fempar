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
module my_par_linear_algebra_names
  use serial_names
  use par_names
  implicit none
# include "debug.i90"
  private
  
  type my_par_linear_algebra_t
     type(par_block_matrix_t)                                  :: p_blk_matrix
     type(par_block_vector_t)                                  :: p_blk_vector
     type(par_block_vector_t)                                  :: p_blk_unknown
     type(par_block_matrix_t)                                  :: mass_p_matrix
     type(block_operator_t)                                    :: block_operator
     type(block_operand_t)                                     :: block_operand_vec
     type(block_operand_t)                                     :: block_operand_unk
     type(block_preconditioner_l_t)                            :: block_preconditioner
     type(par_preconditioner_dd_mlevel_bddc_t)                 :: p_mlevel_bddc_u
     type(par_preconditioner_dd_mlevel_bddc_t)                 :: p_mlevel_bddc_p
     type(par_preconditioner_dd_mlevel_bddc_t)                 :: p_mlevel_bddc_x
     type(par_preconditioner_dd_mlevel_bddc_params_t)          :: p_mlevel_bddc_pars_u
     type(par_preconditioner_dd_mlevel_bddc_params_t)          :: p_mlevel_bddc_pars_p
     type(par_preconditioner_dd_mlevel_bddc_params_t)          :: p_mlevel_bddc_pars_x
     type(par_preconditioner_dd_mlevel_bddc_params_t), pointer :: point_to_p_mlevel_bddc_pars_u 
     type(par_preconditioner_dd_mlevel_bddc_params_t), pointer :: point_to_p_mlevel_bddc_pars_p 
     type(par_preconditioner_dd_mlevel_bddc_params_t), pointer :: point_to_p_mlevel_bddc_pars_x
     type(block_operator_t)                                    :: block_up_operator
     type(block_operator_t)                                    :: block_up_x_operator
     type(block_operator_t)                                    :: block_x_up_operator
     type(block_operator_t)                                    :: block_x_operator
     type(block_operand_t)                                     :: block_up_operand_vec
     type(block_operand_t)                                     :: block_up_operand_unk
     type(block_operand_t)                                     :: block_x_operand_vec
     type(block_operand_t)                                     :: block_x_operand_unk
     type(block_preconditioner_l_t)                            :: block_up_preconditioner
     type(block_preconditioner_l_t)                            :: block_x_preconditioner
   contains
     procedure :: build => build_par_linear_algebra
     procedure :: free  => free_par_linear_algebra
     
  end type my_par_linear_algebra_t

  ! Types
  public :: my_par_linear_algebra_t

contains
  
  !==================================================================================================
  subroutine build_par_linear_algebra(la,p_blk_graph,ndime)
    implicit none
    class(my_par_linear_algebra_t), target, intent(inout) :: la
    type(par_block_graph_t)               , intent(in)    :: p_blk_graph
    integer(ip)                           , intent(in)    :: ndime
    ! Locals
    integer(ip)                :: i,num_levels,istat
    type(par_graph_t), pointer :: p_p_graph
    type(par_matrix_t), pointer :: aux

    ! Set number of levels
    p_p_graph => p_blk_graph%get_block(1,1)
    num_levels = p_p_graph%p_env%num_levels
    
    ! Allocate matrices and vectors
    call la%p_blk_matrix%alloc(p_blk_graph)
    call la%p_blk_matrix%set_block_to_zero(2,2)
    call la%p_blk_matrix%set_block_to_zero(2,3)
    call la%p_blk_matrix%set_block_to_zero(3,2)
    call la%p_blk_vector%alloc(p_blk_graph)
    call la%p_blk_unknown%alloc(p_blk_graph)
    la%p_blk_vector%blocks(1)%state = part_summed
    la%p_blk_vector%blocks(2)%state = part_summed
    la%p_blk_vector%blocks(3)%state = part_summed
    la%p_blk_unknown%blocks(1)%state = full_summed
    la%p_blk_unknown%blocks(2)%state = full_summed
    la%p_blk_unknown%blocks(3)%state = full_summed    
    call par_block_vector_fill_complete(la%p_blk_vector)
    call par_block_vector_fill_complete(la%p_blk_unknown)
    call la%p_blk_vector%init(0.0_rp)

    ! Auxiliar matrices
    call la%mass_p_matrix%alloc(p_blk_graph)
    call la%mass_p_matrix%set_block_to_zero(1,1)
    call la%mass_p_matrix%set_block_to_zero(1,2)
    call la%mass_p_matrix%set_block_to_zero(1,3)
    call la%mass_p_matrix%set_block_to_zero(2,1)
    call la%mass_p_matrix%set_block_to_zero(2,3)
    call la%mass_p_matrix%set_block_to_zero(3,1)
    call la%mass_p_matrix%set_block_to_zero(3,2)
    call la%mass_p_matrix%set_block_to_zero(3,3)

    ! Define (recursive) parameters for K^-1
    la%point_to_p_mlevel_bddc_pars_u => la%p_mlevel_bddc_pars_u
    do i=1, num_levels-1
       la%point_to_p_mlevel_bddc_pars_u%ndime            = ndime
       la%point_to_p_mlevel_bddc_pars_u%unknowns         = all_unknowns
       la%point_to_p_mlevel_bddc_pars_u%pad_collectives  = pad
       la%point_to_p_mlevel_bddc_pars_u%projection       = petrov_galerkin                     
       la%point_to_p_mlevel_bddc_pars_u%subd_elmat_calc  = phit_minus_c_i_t_lambda            !default  
       la%point_to_p_mlevel_bddc_pars_u%correction_mode  = additive_symmetric                 !default 
       la%point_to_p_mlevel_bddc_pars_u%nn_sys_sol_strat = corners_rest_part_solve_expl_schur ! default 
       if(ndime==3) then
          la%point_to_p_mlevel_bddc_pars_u%kind_coarse_dofs = corners_edges_and_faces
       else
          la%point_to_p_mlevel_bddc_pars_u%kind_coarse_dofs = corners_and_edges
       end if
       if ( i < num_levels-1 ) then
          la%point_to_p_mlevel_bddc_pars_u%co_sys_sol_strat     = recursive_bddc
          la%point_to_p_mlevel_bddc_pars_u%ppars_harm%type      = pardiso_mkl_prec !umfpack_prec 
          la%point_to_p_mlevel_bddc_pars_u%ppars_dirichlet%type = pardiso_mkl_prec !umfpack_prec   
          if ( i == 1 ) then
             la%point_to_p_mlevel_bddc_pars_u%spars_coarse%method = direct
             la%point_to_p_mlevel_bddc_pars_u%spars_coarse%itmax  = 200
             la%point_to_p_mlevel_bddc_pars_u%spars_coarse%rtol   = 1.0e-20
             la%point_to_p_mlevel_bddc_pars_u%spars_coarse%trace  = 1
             la%point_to_p_mlevel_bddc_pars_u%correction_mode     = additive
          end if
          allocate(la%point_to_p_mlevel_bddc_pars_u%ppars_coarse_bddc, stat = istat)
          check(istat==0)
          la%point_to_p_mlevel_bddc_pars_u => la%point_to_p_mlevel_bddc_pars_u%ppars_coarse_bddc
       else
          la%point_to_p_mlevel_bddc_pars_u%co_sys_sol_strat         = serial_gather
          la%point_to_p_mlevel_bddc_pars_u%ppars_harm%type          = pardiso_mkl_prec !umfpack_prec  
          la%point_to_p_mlevel_bddc_pars_u%ppars_dirichlet%type     = pardiso_mkl_prec !umfpack_prec  
          la%point_to_p_mlevel_bddc_pars_u%ppars_coarse_serial%type = pardiso_mkl_prec !umfpack_prec  
          nullify ( la%point_to_p_mlevel_bddc_pars_u%ppars_coarse_bddc )
       end if
    end do
    la%point_to_p_mlevel_bddc_pars_u => la%p_mlevel_bddc_pars_u
    do i=1, num_levels-1
       la%point_to_p_mlevel_bddc_pars_u => la%point_to_p_mlevel_bddc_pars_u%ppars_coarse_bddc
    end do

    ! Define (recursive) parameters for Mp^-1
    la%point_to_p_mlevel_bddc_pars_p => la%p_mlevel_bddc_pars_p
    do i=1, num_levels-1
       la%point_to_p_mlevel_bddc_pars_p%ndime            = ndime
       la%point_to_p_mlevel_bddc_pars_p%unknowns         = all_unknowns
       la%point_to_p_mlevel_bddc_pars_p%pad_collectives  = pad
       la%point_to_p_mlevel_bddc_pars_p%projection       = galerkin                     
       la%point_to_p_mlevel_bddc_pars_p%subd_elmat_calc  = phit_minus_c_i_t_lambda            !default  
       la%point_to_p_mlevel_bddc_pars_p%correction_mode  = additive_symmetric                 !default 
       la%point_to_p_mlevel_bddc_pars_p%nn_sys_sol_strat = corners_rest_part_solve_expl_schur ! default 
       if(ndime==3) then
          la%point_to_p_mlevel_bddc_pars_p%kind_coarse_dofs = corners_edges_and_faces
       else
          la%point_to_p_mlevel_bddc_pars_p%kind_coarse_dofs = corners_and_edges
       end if
       if ( i < num_levels-1 ) then
          la%point_to_p_mlevel_bddc_pars_p%co_sys_sol_strat     = recursive_bddc
          la%point_to_p_mlevel_bddc_pars_p%ppars_harm%type      = pardiso_mkl_prec !umfpack_prec 
          la%point_to_p_mlevel_bddc_pars_p%ppars_dirichlet%type = pardiso_mkl_prec !umfpack_prec   
          if ( i == 1 ) then
             la%point_to_p_mlevel_bddc_pars_p%spars_coarse%method = direct
             la%point_to_p_mlevel_bddc_pars_p%spars_coarse%itmax  = 200
             la%point_to_p_mlevel_bddc_pars_p%spars_coarse%rtol   = 1.0e-20
             la%point_to_p_mlevel_bddc_pars_p%spars_coarse%trace  = 1
             la%point_to_p_mlevel_bddc_pars_p%correction_mode     = additive
          end if
          allocate(la%point_to_p_mlevel_bddc_pars_p%ppars_coarse_bddc, stat = istat)
          check(istat==0)
          la%point_to_p_mlevel_bddc_pars_p => la%point_to_p_mlevel_bddc_pars_p%ppars_coarse_bddc
       else
          la%point_to_p_mlevel_bddc_pars_p%co_sys_sol_strat         = serial_gather
          la%point_to_p_mlevel_bddc_pars_p%ppars_harm%type          = pardiso_mkl_prec !umfpack_prec  
          la%point_to_p_mlevel_bddc_pars_p%ppars_dirichlet%type     = pardiso_mkl_prec !umfpack_prec  
          la%point_to_p_mlevel_bddc_pars_p%ppars_coarse_serial%type = pardiso_mkl_prec !umfpack_prec  
          nullify ( la%point_to_p_mlevel_bddc_pars_p%ppars_coarse_bddc )
       end if
    end do
    la%point_to_p_mlevel_bddc_pars_p => la%p_mlevel_bddc_pars_p
    do i=1, num_levels-1
       la%point_to_p_mlevel_bddc_pars_p => la%point_to_p_mlevel_bddc_pars_p%ppars_coarse_bddc
    end do

    ! Define (recursive) parameters for Mx^-1
    la%point_to_p_mlevel_bddc_pars_x => la%p_mlevel_bddc_pars_x
    do i=1, num_levels-1
       la%point_to_p_mlevel_bddc_pars_x%ndime            = ndime
       la%point_to_p_mlevel_bddc_pars_x%unknowns         = all_unknowns
       la%point_to_p_mlevel_bddc_pars_x%pad_collectives  = pad
       la%point_to_p_mlevel_bddc_pars_x%projection       = galerkin                     
       la%point_to_p_mlevel_bddc_pars_x%subd_elmat_calc  = phit_minus_c_i_t_lambda            !default  
       la%point_to_p_mlevel_bddc_pars_x%correction_mode  = additive_symmetric                 !default 
       la%point_to_p_mlevel_bddc_pars_x%nn_sys_sol_strat = corners_rest_part_solve_expl_schur ! default 
       if(ndime==3) then
          la%point_to_p_mlevel_bddc_pars_x%kind_coarse_dofs = corners_edges_and_faces
       else
          la%point_to_p_mlevel_bddc_pars_x%kind_coarse_dofs = corners_and_edges
       end if
       if ( i < num_levels-1 ) then
          la%point_to_p_mlevel_bddc_pars_x%co_sys_sol_strat     = recursive_bddc
          la%point_to_p_mlevel_bddc_pars_x%ppars_harm%type      = pardiso_mkl_prec !umfpack_prec 
          la%point_to_p_mlevel_bddc_pars_x%ppars_dirichlet%type = pardiso_mkl_prec !umfpack_prec   
          if ( i == 1 ) then
             la%point_to_p_mlevel_bddc_pars_x%spars_coarse%method = direct
             la%point_to_p_mlevel_bddc_pars_x%spars_coarse%itmax  = 200
             la%point_to_p_mlevel_bddc_pars_x%spars_coarse%rtol   = 1.0e-20
             la%point_to_p_mlevel_bddc_pars_x%spars_coarse%trace  = 1
             la%point_to_p_mlevel_bddc_pars_x%correction_mode     = additive
          end if
          allocate(la%point_to_p_mlevel_bddc_pars_x%ppars_coarse_bddc, stat = istat)
          check(istat==0)
          la%point_to_p_mlevel_bddc_pars_x => la%point_to_p_mlevel_bddc_pars_x%ppars_coarse_bddc
       else
          la%point_to_p_mlevel_bddc_pars_x%co_sys_sol_strat         = serial_gather
          la%point_to_p_mlevel_bddc_pars_x%ppars_harm%type          = pardiso_mkl_prec !umfpack_prec  
          la%point_to_p_mlevel_bddc_pars_x%ppars_dirichlet%type     = pardiso_mkl_prec !umfpack_prec  
          la%point_to_p_mlevel_bddc_pars_x%ppars_coarse_serial%type = pardiso_mkl_prec !umfpack_prec  
          nullify ( la%point_to_p_mlevel_bddc_pars_x%ppars_coarse_bddc )
       end if
    end do
    la%point_to_p_mlevel_bddc_pars_x => la%p_mlevel_bddc_pars_x
    do i=1, num_levels-1
       la%point_to_p_mlevel_bddc_pars_x => la%point_to_p_mlevel_bddc_pars_x%ppars_coarse_bddc
    end do

    ! Create U-Preconditioner (K^-1)
    call par_preconditioner_dd_mlevel_bddc_create(la%p_blk_matrix%get_block(1,1),la%p_mlevel_bddc_u, &
         &                                        la%p_mlevel_bddc_pars_u)
    call par_preconditioner_dd_mlevel_bddc_ass_struct(la%p_blk_matrix%get_block(1,1),la%p_mlevel_bddc_u)

    ! Create P-Preconditioner (Mp^-1)
    call par_preconditioner_dd_mlevel_bddc_create(la%mass_p_matrix%get_block(2,2),la%p_mlevel_bddc_p, &
         &                                        la%p_mlevel_bddc_pars_p)
    call par_preconditioner_dd_mlevel_bddc_ass_struct(la%mass_p_matrix%get_block(2,2),la%p_mlevel_bddc_p)

    ! Create U-Preconditioner (Mx^-1)
    call par_preconditioner_dd_mlevel_bddc_create(la%p_blk_matrix%get_block(3,3),la%p_mlevel_bddc_x, &
         &                                        la%p_mlevel_bddc_pars_x)
    call par_preconditioner_dd_mlevel_bddc_ass_struct(la%p_blk_matrix%get_block(3,3),la%p_mlevel_bddc_x)
    
    ! Create U-P Block operator
    call la%block_up_operator%create(2,2)
    call la%block_up_operator%set_block(1,1,la%p_blk_matrix%get_block(1,1))
    call la%block_up_operator%set_block(1,2,la%p_blk_matrix%get_block(1,2))
    call la%block_up_operator%set_block(2,1,la%p_blk_matrix%get_block(2,1))
    call la%block_up_operator%set_block_to_zero(2,2)

    ! Create U-P Block operand
    call la%block_up_operand_vec%create(2)
    call la%block_up_operand_vec%set_block(1,la%p_blk_vector%blocks(1))
    call la%block_up_operand_vec%set_block(2,la%p_blk_vector%blocks(2))
    call la%block_up_operand_unk%create(2)
    call la%block_up_operand_unk%set_block(1,la%p_blk_unknown%blocks(1))
    call la%block_up_operand_unk%set_block(2,la%p_blk_unknown%blocks(2))

    ! Create U-P Block preconditioner
    call la%block_up_preconditioner%create(2)
    call la%block_up_preconditioner%set_block(1,1,la%p_mlevel_bddc_u)
    call la%block_up_preconditioner%set_block(2,1,la%p_blk_matrix%get_block(2,1))
    call la%block_up_preconditioner%set_block(2,2,la%p_mlevel_bddc_p)

    ! Create UP-X Block operator
    call la%block_up_x_operator%create(2,1)
    call la%block_up_x_operator%set_block(1,1,la%p_blk_matrix%get_block(1,3))
    call la%block_up_x_operator%set_block_to_zero(2,1) 

    ! Create X-UP Block operator
    call la%block_x_up_operator%create(1,2)
    call la%block_x_up_operator%set_block(1,1,la%p_blk_matrix%get_block(3,1))
    call la%block_x_up_operator%set_block_to_zero(1,2)

    ! Create X Block operator
    call la%block_x_operator%create(1,1)
    call la%block_x_operator%set_block(1,1,la%p_blk_matrix%get_block(3,3))

    ! Create X Block operand
    call la%block_x_operand_vec%create(1)
    call la%block_x_operand_vec%set_block(1,la%p_blk_vector%blocks(3))
    call la%block_x_operand_unk%create(1)
    call la%block_x_operand_unk%set_block(1,la%p_blk_unknown%blocks(3))

    ! Create X Block preconditioner
    call la%block_x_preconditioner%create(1)
    call la%block_x_preconditioner%set_block(1,1,la%p_mlevel_bddc_x)

    ! Create Global Block operator
    call la%block_operator%create(2,2)
    call la%block_operator%set_block(1,1,la%block_up_operator)
    call la%block_operator%set_block(1,2,la%block_up_x_operator)
    call la%block_operator%set_block(2,1,la%block_x_up_operator)
    call la%block_operator%set_block(2,2,la%block_x_operator)

    ! Create Global Block operand
    call la%block_operand_vec%create(2)
    call la%block_operand_vec%set_block(1,la%block_up_operand_vec)
    call la%block_operand_vec%set_block(2,la%block_x_operand_vec)
    call la%block_operand_unk%create(2)
    call la%block_operand_unk%set_block(1,la%block_up_operand_unk)
    call la%block_operand_unk%set_block(2,la%block_x_operand_unk)

    ! Create Global Block preconditioner
    call la%block_preconditioner%create(2)
    call la%block_preconditioner%set_block(1,1,la%block_up_preconditioner)
    call la%block_preconditioner%set_block(2,1,la%block_x_up_operator)
    call la%block_preconditioner%set_block(2,2,la%block_x_preconditioner)

  end subroutine build_par_linear_algebra

  !==================================================================================================
  subroutine free_par_linear_algebra(la)
    implicit none
    class(my_par_linear_algebra_t), intent(inout) :: la

    ! Destroy Global Block preconditioner
    call la%block_preconditioner%destroy()

    ! Destroy Global Block operand
    call la%block_operand_vec%destroy()
    call la%block_operand_unk%destroy()
    
    ! Destroy Global Block operator
    call la%block_operator%destroy()

    ! Destroy X Block preconditioner
    call la%block_x_preconditioner%destroy()

    ! Destroy X Block operand
    call la%block_x_operand_vec%destroy()
    call la%block_x_operand_unk%destroy()

    ! Destroy X Block operator
    call la%block_x_operator%destroy()

    ! Destroy X-UP Block operator
    call la%block_x_up_operator%destroy()

    ! Destroy UP-X Block operator
    call la%block_up_x_operator%destroy()

    ! Destroy U-P Block preconditioner
    call la%block_up_preconditioner%destroy()

    ! Destroy U-P Block operand
    call la%block_up_operand_vec%destroy()
    call la%block_up_operand_unk%destroy()

    ! Destroy U-P Block operator
    call la%block_up_operator%destroy()

    ! Destroy X-preconditioner
    call par_preconditioner_dd_mlevel_bddc_free(la%p_mlevel_bddc_x,free_struct)
    call par_preconditioner_dd_mlevel_bddc_free(la%p_mlevel_bddc_x,free_clean)

    ! Destroy P-preconditioner
    call par_preconditioner_dd_mlevel_bddc_free(la%p_mlevel_bddc_p,free_struct)
    call par_preconditioner_dd_mlevel_bddc_free(la%p_mlevel_bddc_p,free_clean)

    ! Destroy U-preconditioner
    call par_preconditioner_dd_mlevel_bddc_free(la%p_mlevel_bddc_u,free_struct)
    call par_preconditioner_dd_mlevel_bddc_free(la%p_mlevel_bddc_u,free_clean)

    ! Deallocate auxiliar matrices
    call la%mass_p_matrix%free()
          
    ! Deallocate problem matrix
    call la%p_blk_matrix%free()

    ! Deallocate problem vectors
    call la%p_blk_vector%free()
    call la%p_blk_unknown%free()
    
  end subroutine free_par_linear_algebra
  
end module my_par_linear_algebra_names

program par_test_blk_nsi_cg_iss_oss
  use serial_names
  use par_names
  use my_par_linear_algebra_names
  use nsi_names
  use nsi_cg_iss_oss_names
  use norm_names
  use lib_vtk_io_interface_names
  implicit none
# include "debug.i90"

  ! Types
  type(uniform_mesh_descriptor_t)         :: gdata
  type(uniform_conditions_descriptor_t)   :: bdata
  type(reference_element_t)               :: geo_reference_element
  type(par_context_t)                     :: w_context
  type(par_context_t)                     :: p_context
  type(par_context_t)                     :: q_context
  type(par_context_t)                     :: b_context
  type(par_environment_t)                 :: p_env
  type(par_triangulation_t)               :: p_trian
  type(par_conditions_t)                  :: p_cond
  type(dof_descriptor_t)                  :: dof_descriptor
  type(par_fe_space_t)                    :: p_fe_space  
  type(nsi_problem_t)                     :: myprob
  type(nsi_cg_iss_oss_discrete_t)         :: mydisc
  type(time_integration_t)       , target :: tinteg
  type(nsi_cg_iss_oss_matvec_t)  , target :: cg_iss_oss_matvec
  type(discrete_integration_pointer_t)    :: approx(1)
  type(vtk_t)                             :: fevtk
  type(par_block_graph_t)                 :: p_blk_graph
  type(block_dof_distribution_t)          :: blk_dof_dist
  type(par_scalar_t)                      :: enorm_u, enorm_p
  type(error_norm_t)             , target :: error_compute
  type(nsi_cg_iss_oss_massp_t)   , target :: mass_p_integration
  type(my_par_linear_algebra_t)           :: p_linear_algebra
  type(solver_control_t)                  :: sctrl

  ! Integers
  integer(ip) :: num_levels,me,np
  integer(ip) :: gtype(3) = (/ csr, csr, csr /)
  integer(ip) :: ibloc,jbloc,istat,i,j
  integer(ip) :: num_approximations = 1

  ! Parameters
  integer(ip), parameter :: velocity=1, pressure=2

  ! Allocatables
  integer(ip), allocatable :: id_parts(:)
  integer(ip), allocatable :: num_parts(:)
  integer(ip), allocatable :: continuity(:,:)
  integer(ip), allocatable :: order(:,:)
  integer(ip), allocatable :: material(:)
  integer(ip), allocatable :: problem(:)
  integer(ip), allocatable :: vars_block(:)
  integer(ip), allocatable :: dof_coupling(:,:)

  ! Arguments
  character(len=256) :: dir_path_out,prefix
  integer(ip)        :: nex,ney,nez,npx,npy,npz

  call meminit

  ! Read parameters from command-line
  call read_pars_cl_par_test_blk_nsi_cg_iss_oss(prefix,dir_path_out,nex,ney,nez,npx,npy,npz)

  ! Generate geometry data
  call uniform_mesh_descriptor_create(gdata,nex,ney,nez,npx,npy,npz)

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

  ! Create problem
  call myprob%create(gdata%ndime)
  call mydisc%create(myprob)
  call mydisc%vars_block(myprob,vars_block)
  call mydisc%dof_coupling(myprob,dof_coupling)
  call cg_iss_oss_matvec%create(myprob,mydisc)
  cg_iss_oss_matvec%tinteg => tinteg
  tinteg%dtinv      = 0.0_rp
  mydisc%kfl_proj   = 1
  mydisc%kfl_lump   = 1
  myprob%kfl_conv   = 0
  myprob%diffu      = 1.0_rp

  ! Create dof_descriptor
  call dof_descriptor%create(3,1,mydisc%nvars,vars_block,dof_coupling)
  call dof_descriptor%set_problem(1,mydisc)

  ! Allocate auxiliar elemental arrays
  call memalloc(p_trian%f_trian%num_elems,dof_descriptor%nvars_global,continuity, __FILE__,__LINE__)
  call memalloc(p_trian%f_trian%num_elems,dof_descriptor%nvars_global,order,__FILE__,__LINE__)
  call memalloc(p_trian%f_trian%num_elems,problem,__FILE__,__LINE__)
  continuity             = 1
  order                  = 2
  order(:,gdata%ndime+1) = 1
  problem                = 1

  ! Create par_fe_space
  call par_fe_space_create(p_trian,dof_descriptor,p_fe_space,problem,p_cond,continuity,order,material, &
       &                   time_steps_to_store=3,hierarchical_basis=.false.,                           &
       &                   static_condensation=.false.,num_continuity=1)

  ! Initialize VTK output
  call fevtk%initialize(p_trian%f_trian,p_fe_space%fe_space,myprob,p_env,dir_path_out,prefix, &
       &                nparts=gdata%nparts,linear_order=.true.)

  ! Create dof info
  call par_create_distributed_dof_info(dof_descriptor,p_trian,p_fe_space,blk_dof_dist,p_blk_graph,gtype)  

  ! Assign analytical solution
  if(gdata%ndime==2) then
     call par_fe_space_set_analytical_code(p_fe_space,(/4,5,3,0,0/),(/0,0,0,0,0/))
  else
     write(*,*) 'analytical function not ready for 3D'
  end if

  ! Create linear algebra structures
  call p_linear_algebra%build(p_blk_graph,gdata%ndime)

  ! Compute auxiliar matrix
  call mass_p_integration%create(myprob,mydisc)
  approx(1)%p => mass_p_integration
  call par_volume_integral(approx,p_fe_space,p_linear_algebra%mass_p_matrix)
  
  ! Apply boundary conditions to unkno
  call par_update_strong_dirichlet_bcond(p_fe_space,p_cond)
  call par_update_analytical_bcond((/(i,i=1,gdata%ndime+1)/),0.0_rp,p_fe_space)

  ! Solver control parameters
  sctrl%method  = lgmres
  sctrl%trace   = 100
  sctrl%itmax   = 800
  sctrl%dkrymax = 800
  sctrl%stopc   = res_nrmgiven_res_nrmgiven
  sctrl%orto    = icgs
  sctrl%rtol    = 1.0e-14_rp
  sctrl%track_conv_his = .false.

  ! Point discrete integration
  approx(1)%p => cg_iss_oss_matvec

  ! Do nonlinear iterations
  call nonlinear_iteration(sctrl,sctrl%rtol*1.0e2_rp,10,p_env,approx,p_fe_space,p_linear_algebra)

  ! Print solution to VTK file
  istat = fevtk%write_VTK(n_part=p_env%p_context%iam)
  if(p_env%am_i_fine_task()) istat = fevtk%write_PVTK()

  ! Compute error norm
  call error_compute%create(myprob,mydisc)
  approx(1)%p => error_compute
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

  ! Free linear algebra structures
  call p_linear_algebra%free()

  ! Deallocate
  call memfree(continuity,__FILE__,__LINE__)
  call memfree(order,__FILE__,__LINE__)
  call memfree(material,__FILE__,__LINE__)
  call memfree(problem,__FILE__,__LINE__)
  call memfree(vars_block,__FILE__,__LINE__)
  call memfree(dof_coupling,__FILE__,__LINE__)
  call memfree(id_parts , __FILE__, __LINE__)
  call memfree(num_parts, __FILE__, __LINE__)
  call fevtk%free
  call p_blk_graph%free
  call blk_dof_dist%free
  call par_fe_space_free(p_fe_space) 
  call myprob%free
  call mydisc%free
  call error_compute%free
  call cg_iss_oss_matvec%free
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
  subroutine read_pars_cl_par_test_blk_nsi_cg_iss_oss(prefix,dir_path_out,nex,ney,nez,npx,npy,npz)
    implicit none
    character*(*), intent(out) :: prefix, dir_path_out
    integer(ip)  , intent(out) :: nex,ney,nez,npx,npy,npz
    character(len=256)         :: program_name
    character(len=256)         :: argument 
    integer                    :: numargs,iargc

    numargs = iargc()
    call getarg(0, program_name)
    if (.not. (numargs==8) ) then
       write (6,*) 'Usage: ', trim(program_name), ' prefix dir_path_out nex ney nez npx npy npz'
       stop
    end if

    call getarg(1, argument)
    prefix = trim(argument)

    call getarg(2,argument)
    dir_path_out = trim(argument)

    call getarg(3, argument)
    read (argument,*) nex

    call getarg(4, argument)
    read (argument,*) ney

    call getarg(5, argument)
    read (argument,*) nez

    call getarg(6, argument)
    read (argument,*) npx

    call getarg(7, argument)
    read (argument,*) npy

    call getarg(8, argument)
    read (argument,*) npz

  end subroutine read_pars_cl_par_test_blk_nsi_cg_iss_oss

  !==================================================================================================
  subroutine compute_preconditioner(la)
    implicit none
    type(my_par_linear_algebra_t), intent(inout) :: la

    ! precond numeric (K^-1)
    !call par_preconditioner_dd_mlevel_bddc_fill_val(la%p_blk_matrix%get_block(1,1),la%p_mlevel_bddc_u)
    call par_preconditioner_dd_mlevel_bddc_fill_val(la%p_mlevel_bddc_u)

    ! precond numeric (Mp^-1)
    !call par_preconditioner_dd_mlevel_bddc_fill_val(la%mass_p_matrix%get_block(2,2),la%p_mlevel_bddc_p)
    call par_preconditioner_dd_mlevel_bddc_fill_val(la%p_mlevel_bddc_p)

    ! precond numeric (Mx^-1)
    !call par_preconditioner_dd_mlevel_bddc_fill_val(la%p_blk_matrix%get_block(3,3),la%p_mlevel_bddc_x)
    call par_preconditioner_dd_mlevel_bddc_fill_val(la%p_mlevel_bddc_x)
          
  end subroutine compute_preconditioner
  
  !==================================================================================================
  subroutine free_preconditioner(la)
    implicit none
    type(my_par_linear_algebra_t), intent(inout) :: la

    ! precond free (K^-1)
    call par_preconditioner_dd_mlevel_bddc_free(la%p_mlevel_bddc_u,free_values)

    ! precond free (Mp^-1)
    call par_preconditioner_dd_mlevel_bddc_free(la%p_mlevel_bddc_p,free_values)

    ! precond free (Mx^-1)
    call par_preconditioner_dd_mlevel_bddc_free(la%p_mlevel_bddc_x,free_values)
          
  end subroutine free_preconditioner

  !==================================================================================================
  subroutine nonlinear_iteration( sctrl, nltol, maxit, env, approx, p_fe_space, la )
    implicit none
    type(solver_control_t)               , intent(inout) :: sctrl
    real(rp)                             , intent(in)    :: nltol
    integer(ip)                          , intent(in)    :: maxit    
    class(abstract_environment_t)        , intent(in)    :: env
    type(discrete_integration_pointer_t) , intent(inout) :: approx(:)
    type(par_fe_space_t)                 , intent(inout) :: p_fe_space
    type(my_par_linear_algebra_t), target, intent(inout) :: la
    ! Locals
    integer(ip) :: iiter,me,np
    real(rp)    :: resnorm,ininorm
    class(base_operator_t), pointer :: A, M
    class(base_operand_t) , pointer :: x, b
    type(block_operand_t)           :: y
    logical                         :: exit_loop

!!$    ! Assign operators
!!$    A => la%block_operator
!!$    M => la%block_preconditioner
!!$    b => la%block_operand_vec
!!$    x => la%block_operand_unk
!!$
!!$    call p_env%info(me,np)
!!$    
!!$    iiter = 0
!!$    do while( iiter < maxit )
!!$
!!$       ! Update counter
!!$       iiter = iiter+1
!!$
!!$       ! Initialize Matrix and vector
!!$       ! ***************** Abstract procedure to initialize a base_operator ************************!
!!$       call par_block_matrix_zero(la%p_blk_matrix)
!!$       !********************************************************************************************!
!!$       call b%init(0.0_rp)
!!$
!!$       ! Integrate system
!!$       !call volume_integral(approx,fe_space,A,b)
!!$       ! ************ Abstract integration of a base_operator and/or base_operand*******************!
!!$       call par_volume_integral(approx,p_fe_space,la%p_blk_matrix,la%p_blk_vector)
!!$       !********************************************************************************************!
!!$
!!$       ! Check convergence
!!$       if(iiter==1) ininorm = b%nrm2()  
!!$       y = b - A*x 
!!$       resnorm = y%nrm2()
!!$       exit_loop = ( resnorm < nltol*ininorm)
!!$       call p_env%bcast(exit_loop)
!!$       if( exit_loop ) then
!!$          if(me==0) then
!!$             write(*,*) 'Nonlinear iterations: ', iiter
!!$             write(*,*) 'Nonlinear error norm: ', resnorm
!!$          end if
!!$          exit
!!$       end if
!!$
!!$       ! Compute Numeric preconditioner
!!$       ! ***************** Abstract procedure to compute precond numeric ***************************!
!!$       !call compute_preconditioner(la)
!!$       call M%fill_values()
!!$       !********************************************************************************************!
!!$
!!$       ! Solve system
!!$       call abstract_solve(A,M,b,x,sctrl,env)
!!$       call solver_control_log_conv_his(sctrl)
!!$       call solver_control_free_conv_his(sctrl)
!!$
!!$       ! Free Numeric preconditioner
!!$       ! ******************** Abstract procedure to free precond numeric ***************************!
!!$       call free_preconditioner(la)
!!$       !********************************************************************************************!
!!$       
!!$       ! Store solution to unkno
!!$       ! ***************** Abstract procedure to update from a base operant ************************!
!!$       call par_update_solution(la%p_blk_unknown,p_fe_space)
!!$       !********************************************************************************************!
!!$       
!!$       ! Store nonlinear iteration ( k+1 --> k )
!!$       call par_update_nonlinear_solution(p_fe_space)
!!$       
!!$    end do
!!$
!!$    ! Deallocate
!!$    call y%free()

  end subroutine nonlinear_iteration


end program par_test_blk_nsi_cg_iss_oss
