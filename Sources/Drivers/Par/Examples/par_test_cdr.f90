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
program par_test_cdr
  !----------------------------------------------------------
  ! Parallel partitioner test
  !----------------------------------------------------------
  use fem
  use par
  use cdr_names
  use cdr_stabilized_continuous_Galerkin_names 
  use mpi
  
  implicit none
#include "debug.i90" 
  ! Our data
  type(par_context)       :: w_context, p_context, q_context, b_context
  type(par_environment)   :: p_env
  type(par_mesh)          :: p_mesh
  type(par_triangulation) :: p_trian
  type(par_fem_space)     :: p_fspac

  type(par_matrix), target                :: p_mat
  type(par_vector), target                :: p_vec, p_unk
  class(base_operand) , pointer           :: x, y
  class(base_operator), pointer           :: A

  ! Preconditioner-related data structures
  type(par_precond_dd_diagonal)           :: p_prec_dd_diag
  type(par_precond_dd_mlevel_bddc), target :: p_mlevel_bddc
  ! type(par_precond_dd_mlevel_bddc), pointer  :: point_to_p_mlevel_bddc
  type(par_precond_dd_mlevel_bddc_params), target  :: p_mlevel_bddc_pars
  type(par_precond_dd_mlevel_bddc_params), pointer :: point_to_p_mlevel_bddc_pars
  integer(ip), allocatable :: kind_coarse_dofs(:)


  type(solver_control)     :: sctrl

  type(block_dof_distribution)    :: blk_dof_dist
  type(dof_handler)               :: dhand
  type(par_block_graph)           :: p_blk_graph
  integer(ip)                     :: gtype(1) = (/ csr_symm /)
  type(par_conditions)            :: p_cond

  type(cdr_problem)               :: my_problem
  type(cdr_approximation), target :: my_approximation
  integer(ip)                     :: num_approximations
  type(discrete_problem_pointer)  :: approximations(1)

  integer(ip)              :: num_levels, nparts, ndime
  integer(ip), allocatable :: id_parts(:), num_parts(:)

  ! Arguments
  integer(ip)                   :: lunio
  character(len=256)            :: dir_path, dir_path_out
  character(len=256)            :: prefix
  character(len=:), allocatable :: name
  integer(ip)              :: i, j, ierror, iblock

  integer(ip), allocatable :: order(:,:), material(:), problem(:), which_approx(:)
  integer(ip), allocatable :: continuity(:,:)

  call meminit

  ! Read parameters from command-line
  call read_pars_cl ( dir_path, prefix, dir_path_out, nparts, ndime )

  !! This is a hack. It only works for two levels.
  num_levels = 2
  call memalloc(num_levels, id_parts , __FILE__, __LINE__)
  call memalloc(num_levels, num_parts, __FILE__, __LINE__)

  num_parts = (/nparts, 1/)
  id_parts = (/w_context%iam+1, 1/)


  ! Start parallel execution
  call par_context_create (w_context)

  ! Create p_context and q_context splitting w_context
  if(w_context%iam < num_parts(1)) then
     call par_context_create ( 1, p_context, q_context, w_context )
  else
     call par_context_create ( 2, q_context, p_context, w_context )
  end if
  assert ( (p_context%iam >=0 .and. q_context%iam < 0) .or. (p_context%iam < 0 .and. q_context%iam >= 0))
  
  ! Create b_context as an intercommunicator among p_context <=> q_context 
  call par_context_create ( w_context, p_context, q_context, b_context )

  ! Create parallel environment
  call par_environment_create( p_env,& 
                               w_context,& 
                               p_context,& 
                               q_context,&
                               b_context,&
                               num_levels,&
                               id_parts, & 
                               num_parts )

  ! Read mesh
  call par_mesh_read ( dir_path, prefix, p_env, p_mesh )

  ! Read boundary conditions
  call par_conditions_read(dir_path, prefix, p_mesh%f_mesh%npoin, p_env, p_cond)
  if ( p_env%am_i_fine_task() ) p_cond%f_conditions%code = 0 !(dG)

  call par_mesh_to_triangulation (p_mesh, p_trian, p_cond)

  !write (*,*) '********** CREATE DOF HANDLER**************'
  call dhand%create( 1, 1, 1 )


  call my_problem%create( p_trian%f_trian%num_dims )
  call my_approximation%create(my_problem)
  num_approximations=1
  approximations(1)%p => my_approximation
  
  call dhand%set_problem( 1, my_approximation )
  ! ... for as many problems as we have

  call memalloc( p_trian%f_trian%num_elems, dhand%nvars_global, continuity, __FILE__, __LINE__)
  continuity = 0 !(dG)
  call memalloc( p_trian%f_trian%num_elems, dhand%nvars_global, order, __FILE__, __LINE__)
  order = 1
  call memalloc( p_trian%f_trian%num_elems, material, __FILE__, __LINE__)
  material = 1
  call memalloc( p_trian%f_trian%num_elems, problem, __FILE__, __LINE__)
  problem = 1
  call memalloc( p_trian%f_trian%num_elems, which_approx, __FILE__, __LINE__)
  which_approx = 1


  ! Continuity
  ! write(*,*) 'Continuity', continuity
  call par_fem_space_create ( p_trian, dhand, p_fspac, problem, &
                              num_approximations, approximations, &
                              p_cond, continuity, order, material, &
                              which_approx, time_steps_to_store = 1, &
                              hierarchical_basis = logical(.false.,lg), &
                              & static_condensation = logical(.false.,lg), num_continuity = 1 )

  if ( p_env%am_i_fine_task() ) then
     call update_strong_dirichlet_boundary_conditions( p_fspac%f_space )
  end if

  call par_create_distributed_dof_info ( dhand, p_trian, p_fspac, blk_dof_dist, p_blk_graph, gtype )  

  call par_matrix_alloc ( csr_mat, symm_true, p_blk_graph%get_block(1,1), p_mat, positive_definite )

  call par_vector_alloc ( blk_dof_dist%get_block(1), p_env, p_vec )
  p_vec%state = part_summed

  call par_vector_alloc ( blk_dof_dist%get_block(1), p_env, p_unk )
  p_unk%state = full_summed

  if ( p_env%am_i_fine_task() ) then
     call volume_integral( approximations, p_fspac%f_space, p_mat%f_matrix, p_vec%f_vector)
     !call fem_matrix_print ( 6, p_mat%f_matrix )
     !call fem_vector_print ( 6, p_vec%f_vector )
  end if

  call p_unk%init(1.0_rp)

  A => p_mat
  x => p_vec
  y => p_unk
  y = x - A*y
  write(*,*) 'XXX error norm XXX', y%nrm2()
  ! AFM: I had to re-assign the state of punk as the expression
  ! y = x - A*y changed its state to part_summed!!! 
  p_unk%state = full_summed

  ! Define (recursive) parameters
  point_to_p_mlevel_bddc_pars => p_mlevel_bddc_pars
  do i=1, num_levels-1
     point_to_p_mlevel_bddc_pars%ndime            = ndime
     point_to_p_mlevel_bddc_pars%unknowns         = all_unknowns
     point_to_p_mlevel_bddc_pars%pad_collectives  = pad
     point_to_p_mlevel_bddc_pars%projection       = galerkin                           !default
     point_to_p_mlevel_bddc_pars%subd_elmat_calc  = phit_minus_c_i_t_lambda            !default  
     point_to_p_mlevel_bddc_pars%correction_mode  = additive_symmetric                 !default 
     point_to_p_mlevel_bddc_pars%nn_sys_sol_strat = corners_rest_part_solve_expl_schur ! default 

     if ( i < num_levels-1 ) then
        point_to_p_mlevel_bddc_pars%co_sys_sol_strat = recursive_bddc
        point_to_p_mlevel_bddc_pars%ppars_harm%type      = pardiso_mkl_prec !umfpack_prec 
        point_to_p_mlevel_bddc_pars%ppars_dirichlet%type = pardiso_mkl_prec !umfpack_prec 
        
        if ( i == 1 ) then
           point_to_p_mlevel_bddc_pars%spars_coarse%method = direct
           point_to_p_mlevel_bddc_pars%spars_coarse%itmax  = 200
           point_to_p_mlevel_bddc_pars%spars_coarse%rtol   = 1.0e-08
           point_to_p_mlevel_bddc_pars%spars_coarse%trace  = 1
           point_to_p_mlevel_bddc_pars%correction_mode  = additive
        end if
        allocate(point_to_p_mlevel_bddc_pars%ppars_coarse_bddc, stat = ierror)
        check(ierror==0)
        point_to_p_mlevel_bddc_pars => point_to_p_mlevel_bddc_pars%ppars_coarse_bddc
     else
        point_to_p_mlevel_bddc_pars%co_sys_sol_strat = serial_gather
        point_to_p_mlevel_bddc_pars%ppars_harm%type          =pardiso_mkl_prec !umfpack_prec  
        point_to_p_mlevel_bddc_pars%ppars_dirichlet%type     =pardiso_mkl_prec !umfpack_prec  
        point_to_p_mlevel_bddc_pars%ppars_coarse_serial%type =pardiso_mkl_prec !umfpack_prec  
        nullify ( point_to_p_mlevel_bddc_pars%ppars_coarse_bddc )
     end if
  end do

  call memalloc ( ndime, kind_coarse_dofs, __FILE__, __LINE__ )
  kind_coarse_dofs(1) = corners
  kind_coarse_dofs(2) = corners_and_edges
  if ( ndime == 3 ) then
     kind_coarse_dofs(3) = corners_edges_and_faces
  end if

  sctrl%method=cg
  sctrl%trace=1
  sctrl%itmax=800
  sctrl%dkrymax=800
  sctrl%stopc=res_res
  sctrl%orto=icgs
  sctrl%rtol=1.0e-06

  do j=1,ndime

     point_to_p_mlevel_bddc_pars => p_mlevel_bddc_pars
     do i=1, num_levels-1
        point_to_p_mlevel_bddc_pars%kind_coarse_dofs = kind_coarse_dofs(j)
        point_to_p_mlevel_bddc_pars => point_to_p_mlevel_bddc_pars%ppars_coarse_bddc
     end do


     call p_unk%init(0.0_rp)

     ! Create multilevel bddc inverse 
     call par_precond_dd_mlevel_bddc_create( p_mat, p_mlevel_bddc, p_mlevel_bddc_pars )

     ! Ass struct
     call par_precond_dd_mlevel_bddc_ass_struct ( p_mat, p_mlevel_bddc )

     ! Fill val
     call par_precond_dd_mlevel_bddc_fill_val ( p_mat, p_mlevel_bddc )

     call abstract_solve(p_mat,p_mlevel_bddc,p_vec,p_unk,sctrl,p_env)

     ! Free bddc inverse
     call par_precond_dd_mlevel_bddc_free( p_mlevel_bddc, free_only_values)
     call par_precond_dd_mlevel_bddc_free( p_mlevel_bddc, free_only_struct)
     call par_precond_dd_mlevel_bddc_free( p_mlevel_bddc, free_clean)


  end do

  call memfree ( kind_coarse_dofs, __FILE__, __LINE__ )

!!$  call par_precond_dd_diagonal_create ( p_mat, p_prec_dd_diag )
!!$  call par_precond_dd_diagonal_ass_struct ( p_mat, p_prec_dd_diag )
!!$  call par_precond_dd_diagonal_fill_val ( p_mat, p_prec_dd_diag )
!!$  
!!$
!!$  call abstract_solve(p_mat,p_prec_dd_diag,p_vec,p_unk,sctrl,p_env)
!!$
!!$  call par_precond_dd_diagonal_free ( p_prec_dd_diag, free_only_values )
!!$  call par_precond_dd_diagonal_free ( p_prec_dd_diag, free_only_struct )
!!$  call par_precond_dd_diagonal_free ( p_prec_dd_diag, free_clean )



  call par_matrix_free (p_mat)
  call par_vector_free (p_vec)
  call par_vector_free (p_unk)

  call memfree( continuity, __FILE__, __LINE__)
  call memfree( order, __FILE__, __LINE__)
  call memfree( material, __FILE__, __LINE__)
  call memfree( problem, __FILE__, __LINE__)
  call memfree( which_approx, __FILE__, __LINE__)

  call p_blk_graph%free
  call blk_dof_dist%free
  call par_fem_space_free(p_fspac) 
  call my_problem%free
  call my_approximation%free
  call dof_handler_free (dhand)
  call par_triangulation_free(p_trian)
  call par_conditions_free (p_cond)
  call par_mesh_free (p_mesh)

  call memfree(id_parts , __FILE__, __LINE__)
  call memfree(num_parts, __FILE__, __LINE__)

  call par_environment_free (p_env)
  call par_context_free ( b_context, .false. )
  call par_context_free ( p_context, .false. )
  call par_context_free ( q_context, .false. )
  call par_context_free ( w_context )

  call memstatus

contains
  subroutine read_pars_cl (dir_path, prefix, dir_path_out, nparts, ndime)
    implicit none
    character*(*), intent(out)   :: dir_path, prefix, dir_path_out
    integer(ip)  , intent(out)   :: nparts, ndime

    character(len=256)           :: program_name
    character(len=256)           :: argument 
    integer                      :: numargs,iargc

    numargs = iargc()
    call getarg(0, program_name)
    if (.not. (numargs==5) ) then
       write (6,*) 'Usage: ', trim(program_name), ' dir_path prefix dir_path_out nparts ndime'
       stop
    end if

    call getarg(1, argument)
    dir_path = trim(argument)

    call getarg(2, argument)
    prefix = trim(argument)

    call getarg(3,argument)
    dir_path_out = trim(argument)

    call getarg(4,argument)
    read(argument,*) nparts

    call getarg(5,argument)
    read(argument,*) ndime

  end subroutine read_pars_cl

  subroutine update_strong_dirichlet_boundary_conditions( fspac )
    implicit none
    type(fem_space), intent(inout)    :: fspac
    
    integer(ip) :: ielem, iobje, ivar, inode, l_node

    do ielem = 1, fspac%g_trian%num_elems
       do iobje = 1,fspac%lelem(ielem)%p_geo_info%nobje
          do ivar=1, fspac%dof_handler%problems(problem(ielem))%p%nvars
             
             do inode = fspac%lelem(ielem)%nodes_object(ivar)%p%p(iobje), &
                  &     fspac%lelem(ielem)%nodes_object(ivar)%p%p(iobje+1)-1 
                l_node = fspac%lelem(ielem)%nodes_object(ivar)%p%l(inode)
                if ( fspac%lelem(ielem)%bc_code(ivar,iobje) /= 0 ) then
                   fspac%lelem(ielem)%unkno(l_node,ivar,1) = 1.0_rp
                end if
             end do
          end do
       end do
    end do
  end subroutine update_strong_dirichlet_boundary_conditions

end program par_test_cdr
