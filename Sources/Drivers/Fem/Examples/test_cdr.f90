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
  use fem 
  use cdr_names
  use cdr_stabilized_continuous_Galerkin_names
  implicit none
#include "debug.i90"
  ! Our data
  type(fem_mesh)           :: f_mesh
  type(fem_triangulation)  :: f_trian
  type(fem_matrix)         :: f_mat
  type(fem_conditions)     :: f_cond
  type(dof_handler)        :: dhand
  type(fem_space)          :: fspac
  type(fem_graph), pointer :: f_graph
  type(fem_block_graph)    :: f_blk_graph
  integer(ip)              :: gtype(1) = (/ csr_symm /)

  type(cdr_problem)               :: my_problem
  type(cdr_approximation), target :: my_approximation
  integer(ip)                     :: num_approximations
  type(discrete_problem_pointer)  :: approximations(1)

  type(fem_matrix), target             :: my_matrix
  type(fem_vector), target             :: my_vector, feunk
  class(base_operand) , pointer :: x, y
  class(base_operator), pointer :: A

  type(fem_precond)        :: feprec
  type(fem_precond_params) :: ppars
  type(solver_control)     :: sctrl
  type(serial_environment) :: senv

  ! Arguments
  character(len=256)       :: dir_path, dir_path_out
  character(len=256)       :: prefix, filename
  integer(ip)              :: i, j, vars_prob(1) = 1, ierror, iblock

  integer(ip), allocatable :: order(:,:), material(:), problem(:), which_approx(:)

  integer(ip), allocatable :: continuity(:,:)

  integer(ip) :: lunio, istat

  call meminit

  ! Read parameters from command-line
  call  read_pars_cl_test_cdr ( dir_path, prefix, dir_path_out )

  ! Read mesh
  call fem_mesh_read (dir_path, prefix, f_mesh, permute_c2z=.true.)

  ! Read conditions 
  call fem_conditions_read (dir_path, prefix, f_mesh%npoin, f_cond)

  !call fem_mesh_write ( 6, f_mesh, 'square_uniform' )

  call mesh_to_triangulation ( f_mesh, f_trian, gcond = f_cond )

  ! write(*,*) 'conditions%ncode', f_cond%ncode
  ! write(*,*) 'conditions%nvalu', f_cond%nvalu
  ! write(*,*) 'conditions%ncond', f_cond%ncond

  ! write(*,*) 'conditions%code', f_cond%code
  ! write(*,*) 'conditions%valu', f_cond%valu
  f_cond%code = 0 ! (dG)
  !call triangulation_print( 6 , f_trian )

  vars_prob = 1
  call dhand%create( 1, 1, 1 ) !, vars_block, dof_coupl )
  !                      ( dhand, nblocks, nprobs, nvars_global, vars_block, dof_coupl )

  !write (6,*) '*** physical problem  ***'
  !write (6,*) 'Number of variables of problem: ',  my_problem%nvars
  !write (6,*) 'Local to global (of variables) for problem: ' ,  my_problem%l2g_var

  call my_problem%create( f_trian%num_dims )
  call my_approximation%create(my_problem)
  num_approximations=1
  approximations(1)%p => my_approximation

  call dhand%set_problem( 1, my_approximation )
  !                     ( ndime, dhand, l2g_vars, iprob ) 
  ! ... for as many problems as we have

  !call dof_handler_print ( dhand, 6 )


  call memalloc( f_trian%num_elems, dhand%nvars_global, continuity, __FILE__, __LINE__)
   continuity = 0 ! (dG)
  ! continuity = 1
  call memalloc( f_trian%num_elems, dhand%nvars_global, order, __FILE__, __LINE__)
  order = 5
  call memalloc( f_trian%num_elems, material, __FILE__, __LINE__)
  material = 1
  call memalloc( f_trian%num_elems, problem, __FILE__, __LINE__)
  problem = 1
  call memalloc( f_trian%num_elems, which_approx, __FILE__, __LINE__)
  which_approx = 1 

  ! Continuity
  !write(*,*) 'Continuity', continuity

  call fem_space_create ( f_trian, dhand, fspac, problem, num_approximations, approximations, f_cond, continuity, order, material, &
       & which_approx=which_approx,  time_steps_to_store = 1, hierarchical_basis = logical(.false.,lg), & 
       & static_condensation = logical(.false.,lg), num_continuity = 1 )

  call update_strong_dirichlet_boundary_conditions( fspac )

  call create_dof_info( dhand, f_trian, fspac, f_blk_graph, gtype )

  f_graph => f_blk_graph%get_block(1,1)
  call fem_matrix_alloc( csr_mat, symm_true, f_graph, my_matrix, positive_definite )

  call fem_vector_alloc( f_graph%nv, my_vector )
  
  call volume_integral( approximations, fspac, my_matrix, my_vector)

  !sctrl%method=direct
  !ppars%type = pardiso_mkl_prec
  !call fem_precond_create  (my_matrix, feprec, ppars)
  !call fem_precond_symbolic(my_matrix, feprec)
  !call fem_precond_numeric (my_matrix, feprec)
  !call fem_precond_log_info(feprec)

  write (*,*) '********** STARTING RES COMP **********,dof_graph(1,1)%nv',f_graph%nv
  call fem_vector_alloc( f_graph%nv, feunk )
  call feunk%init(1.0_rp)



  ! feunk = my_vector - my_matrix*feunk 
  !y = x - A*y 

  !call fem_vector_print( 6, feunk)
  write(*,*) 'XXX error vs exact norm XXX', feunk%nrm2()
  
  !call abstract_solve(my_matrix,feprec,my_vector,feunk,sctrl,senv)

  !call solver_control_free_conv_his(sctrl)

  A => my_matrix
  x => my_vector
  y => feunk
  call fem_vector_print( 6, feunk)
  call fem_matrix_print( 6, my_matrix)

  ! feunk = my_vector - my_matrix*feunk 
  y = x - A*y 

  write(*,*) 'XXX error solver norm XXX', feunk%nrm2()

  !call fem_vector_print( 6, feunk)

  !call fem_precond_free ( precond_free_values, feprec)
  !call fem_precond_free ( precond_free_struct, feprec)
  !call fem_precond_free ( precond_free_clean, feprec)

  !write (*,*) '********** FINISHED ASSEMBLY **********'


  ! call fem_matrix_print( 6, my_matrix)
  ! call fem_precond_dd_mlevel_bddc_create ( f_mat, mlbddc, mlbddc_params )

  call memfree( continuity, __FILE__, __LINE__)
  call memfree( order, __FILE__, __LINE__)
  call memfree( material, __FILE__, __LINE__)
  call memfree( problem, __FILE__, __LINE__)
  call memfree( which_approx, __FILE__, __LINE__)

  call f_blk_graph%free()
  call fem_vector_free( feunk )
  call fem_vector_free( my_vector )
  call fem_matrix_free( my_matrix) 
  call fem_space_free(fspac) 
  call my_problem%free
  call my_approximation%free
  call dof_handler_free ( dhand )
  call fem_triangulation_free ( f_trian )
  call fem_conditions_free ( f_cond )
  call fem_mesh_free (f_mesh)

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

  subroutine update_strong_dirichlet_boundary_conditions( fspac )
    implicit none

    type(fem_space), intent(inout)    :: fspac

    integer(ip) :: ielem, iobje, ivar, inode, l_node

    do ielem = 1, fspac%g_trian%num_elems
       do ivar=1, fspac%dof_handler%problems(problem(ielem))%p%nvars
          !write (*,*) 'ielem',ielem
          !write (*,*) 'ivar',ivar
          !write (*,*) 'KKKKKKKKKKKKKKKKKKKKK'
          !write (*,*) 'fspac%lelem(ielem)%nodes_object(ivar)%p%p',fspac%lelem(ielem)%nodes_object(ivar)%p%p
          !write (*,*) 'fspac%lelem(ielem)%nodes_object(ivar)%p%l',fspac%lelem(ielem)%nodes_object(ivar)%p%l
          do iobje = 1,fspac%lelem(ielem)%p_geo_info%nobje

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

end program test_cdr
