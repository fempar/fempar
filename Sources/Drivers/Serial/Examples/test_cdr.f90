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
program test_cdr
  use serial_names
  use cdr_names
  use cdr_stabilized_continuous_Galerkin_names
  implicit none
#include "debug.i90"
  ! Our data
  type(mesh_t)           :: f_mesh
  type(triangulation_t)  :: f_trian
  type(matrix_t)         :: f_mat
  type(conditions_t)     :: f_cond
  type(dof_descriptor_t) :: dof_descriptor
  type(fe_space_t)       :: fe_space
  type(graph_t), pointer :: f_graph
  type(block_graph_t)    :: f_blk_graph
  logical                :: symmetric_storage(1) = (/ .true. /)

  type(cdr_problem_t)                   :: my_problem
  type(cdr_discrete_t)                  :: my_discrete
  type(cdr_approximation_t), target     :: my_approximation
  integer(ip)                         :: num_approximations
  type(discrete_integration_pointer_t)  :: approximations(1)

  type(matrix_t), target             :: my_matrix
  type(vector_t), target             :: my_vector, feunk
  class(abstract_vector_t) , pointer :: x, y
  class(abstract_operator_t), pointer :: A

  type(preconditioner_t)        :: feprec
  type(preconditioner_params_t) :: ppars
  type(solver_control_t)     :: sctrl
  type(serial_environment_t) :: senv

  ! Arguments
  character(len=256)       :: dir_path, dir_path_out
  character(len=256)       :: prefix, filename
  integer(ip)              :: i, j, vars_prob(1) = 1, ierror, iblock

  integer(ip), allocatable :: order(:,:), material(:), problem(:)

  integer(ip), allocatable :: continuity(:,:)

  integer(ip) :: lunio, istat

  call meminit

  ! Read parameters from command-line
  call  read_pars_cl_test_cdr ( dir_path, prefix, dir_path_out )

  ! Read mesh
  call mesh_read (dir_path, prefix, f_mesh, permute_c2z=.true.)

  ! Read conditions 
  call conditions_read (dir_path, prefix, f_mesh%npoin, f_cond)

  !call mesh_write ( 6, f_mesh, 'square_uniform' )

  call mesh_to_triangulation ( f_mesh, f_trian, gcond = f_cond )

  ! write(*,*) 'conditions%ncode', f_cond%ncode
  ! write(*,*) 'conditions%nvalu', f_cond%nvalu
  ! write(*,*) 'conditions%ncond', f_cond%ncond

  ! write(*,*) 'conditions%code', f_cond%code
  ! write(*,*) 'conditions%valu', f_cond%valu
  f_cond%code = 0 ! (dG)
  !call triangulation_print( 6 , f_trian )

  vars_prob = 1
  call dof_descriptor%create( 1, 1, 1 ) !, vars_block, dof_coupl )
  !                      ( dof_descriptor, nblocks, nprobs, nvars_global, vars_block, dof_coupl )

  !write (6,*) '*** physical problem  ***'
  !write (6,*) 'Number of variables of problem: ',  my_problem%nvars
  !write (6,*) 'Local to global (of variables) for problem: ' ,  my_problem%l2g_var

  call my_problem%create( f_trian%num_dims )
  call my_discrete%create( my_problem)
  call my_approximation%create(my_problem,my_discrete)
  num_approximations=1
  approximations(1)%p => my_approximation

  call dof_descriptor%set_problem( 1, my_discrete )
  !                     ( ndime, dof_descriptor, l2g_vars, iprob ) 
  ! ... for as many problems as we have

  !call dof_descriptor_print ( dof_descriptor, 6 )


  call memalloc( f_trian%num_elems, dof_descriptor%nvars_global, continuity, __FILE__, __LINE__)
   continuity = 0 ! (dG)
  ! continuity = 1
  call memalloc( f_trian%num_elems, dof_descriptor%nvars_global, order, __FILE__, __LINE__)
  order = 5
  call memalloc( f_trian%num_elems, material, __FILE__, __LINE__)
  material = 1
  call memalloc( f_trian%num_elems, problem, __FILE__, __LINE__)
  problem = 1

  ! Continuity
  !write(*,*) 'Continuity', continuity

  call fe_space_create ( f_trian, dof_descriptor, fe_space, problem, f_cond, continuity, order, material, &
       & time_steps_to_store = 1, hierarchical_basis = .false., & 
       & static_condensation = .false., num_continuity = 1 )

  f_cond%valu = 1.0_rp
  call update_strong_dirichlet_bcond( fe_space, f_cond )

  call create_dof_info( dof_descriptor, f_trian, fe_space, f_blk_graph, symmetric_storage )

  f_graph => f_blk_graph%get_block(1,1)
  call matrix_alloc( .true., f_graph, my_matrix, positive_definite )

  call vector_alloc( f_graph%nv, my_vector )
  
  call volume_integral( approximations, fe_space, my_matrix, my_vector)

  !sctrl%method=direct
  !ppars%type = pardiso_mkl_prec
  !call preconditioner_create  (my_matrix, feprec, ppars)
  !call preconditioner_symbolic(my_matrix, feprec)
  !call preconditioner_numeric (my_matrix, feprec)
  !call preconditioner_log_info(feprec)

  write (*,*) '********** STARTING RES COMP **********,dof_graph(1,1)%nv',f_graph%nv
  call vector_alloc( f_graph%nv, feunk )
  call feunk%init(1.0_rp)



  ! feunk = my_vector - my_matrix*feunk 
  !y = x - A*y 

  !call vector_print( 6, feunk)
  write(*,*) 'XXX error vs exact norm XXX', feunk%nrm2()
  
  !call abstract_solve(my_matrix,feprec,my_vector,feunk,sctrl,senv)

  !call solver_control_free_conv_his(sctrl)

  A => my_matrix
  x => my_vector
  y => feunk
  call vector_print( 6, feunk)
  call matrix_print( 6, my_matrix)

  ! feunk = my_vector - my_matrix*feunk 
  y = x - A*y 

  write(*,*) 'XXX error solver norm XXX', feunk%nrm2()

  !call vector_print( 6, feunk)

  !call preconditioner_free ( preconditioner_free_values, feprec)
  !call preconditioner_free ( preconditioner_free_struct, feprec)
  !call preconditioner_free ( preconditioner_free_clean, feprec)

  !write (*,*) '********** FINISHED ASSEMBLY **********'


  ! call matrix_print( 6, my_matrix)
  ! call preconditioner_dd_mlevel_bddc_create ( f_mat, mlbddc, mlbddc_params )

  call memfree( continuity, __FILE__, __LINE__)
  call memfree( order, __FILE__, __LINE__)
  call memfree( material, __FILE__, __LINE__)
  call memfree( problem, __FILE__, __LINE__)

  call f_blk_graph%free()
  call vector_free( feunk )
  call vector_free( my_vector )
  call matrix_free( my_matrix) 
  call fe_space_free(fe_space) 
  call my_problem%free
  call my_discrete%free
  call my_approximation%free
  call dof_descriptor_free ( dof_descriptor )
  call triangulation_free ( f_trian )
  call conditions_free ( f_cond )
  call mesh_free (f_mesh)

  call memstatus

contains

  subroutine read_pars_cl_test_cdr (dir_path, prefix, dir_path_out)
    implicit none
    character*(*), intent(out)   :: dir_path, prefix, dir_path_out
    character(len=256)           :: program_name
    character(len=256)           :: argument 
    integer                      :: numargs,iargc

    numargs = iargc()
    call getarg(0, program_name)
    if (.not. (numargs==3) ) then
       write (6,*) 'Usage: ', trim(program_name), ' dir_path prefix dir_path_out'
       stop
    end if

    call getarg(1, argument)
    dir_path = trim(argument)

    call getarg(2, argument)
    prefix = trim(argument)

    call getarg(3,argument)
    dir_path_out = trim(argument)

  end subroutine read_pars_cl_test_cdr

!!$  subroutine update_strong_dirichlet_boundary_conditions( fe_space )
!!$    implicit none
!!$
!!$    type(fe_space_t), intent(inout)    :: fe_space
!!$
!!$    integer(ip) :: ielem, iobje, ivar, inode, l_node
!!$
!!$    do ielem = 1, fe_space%g_trian%num_elems
!!$       do ivar=1, fe_space%dof_descriptor%problems(problem(ielem))%p%nvars
!!$          !write (*,*) 'ielem',ielem
!!$          !write (*,*) 'ivar',ivar
!!$          !write (*,*) 'KKKKKKKKKKKKKKKKKKKKK'
!!$          !write (*,*) 'fe_space%lelem(ielem)%nodes_per_vef(ivar)%p%p',fe_space%lelem(ielem)%nodes_per_vef(ivar)%p%p
!!$          !write (*,*) 'fe_space%lelem(ielem)%nodes_per_vef(ivar)%p%l',fe_space%lelem(ielem)%nodes_per_vef(ivar)%p%l
!!$          do iobje = 1,fe_space%lelem(ielem)%p_geo_reference_element%nvef
!!$
!!$             do inode = fe_space%lelem(ielem)%nodes_per_vef(ivar)%p%p(iobje), &
!!$                  &     fe_space%lelem(ielem)%nodes_per_vef(ivar)%p%p(iobje+1)-1 
!!$                l_node = fe_space%lelem(ielem)%nodes_per_vef(ivar)%p%l(inode)
!!$                if ( fe_space%lelem(ielem)%bc_code(ivar,iobje) /= 0 ) then
!!$                   fe_space%lelem(ielem)%unkno(l_node,ivar,1) = 1.0_rp
!!$                end if
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$
!!$  end subroutine update_strong_dirichlet_boundary_conditions

end program test_cdr
