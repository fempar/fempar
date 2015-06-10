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
  type(par_context)                       :: context
  type(par_environment)                   :: p_env
  type(par_mesh)                          :: p_mesh
  type(par_triangulation)                 :: p_trian
  type(par_matrix)                        :: p_mat
  type(par_vector)                        :: p_vec, p_unk

  type(dof_distribution) , allocatable :: dof_dist(:)

  type(dof_handler)  :: dhand
  type(fem_space)    :: fspac

  type(par_graph), allocatable    :: dof_graph(:,:)
  type(fem_conditions)  :: f_cond

  type(cdr_problem)               :: my_problem
  type(cdr_approximation), target :: my_approximation

  ! Arguments
  integer(ip)                   :: lunio
  character(len=256)            :: dir_path, dir_path_out
  character(len=256)            :: prefix
  character(len=:), allocatable :: name
  integer(ip)              :: i, j, ierror, iblock

  integer(ip), allocatable :: order(:,:), material(:), problem(:), which_approx(:)

  integer(ip), allocatable :: continuity(:,:)

  type(discrete_problem_pointer) :: approximations(1)

  call meminit

  ! Start parallel execution
  call par_context_create (context)

  ! Create parallel environment
  call par_environment_create( p_env, context )

  ! Read parameters from command-line
  call read_pars_cl ( dir_path, prefix, dir_path_out )

  ! Read mesh
  call par_mesh_read ( dir_path, prefix, p_env, p_mesh )

  ! Read conditions (next 5 lines should go in a subroutine...should they?)
  call fem_conditions_compose_name(prefix,name) 
  call par_filename(context,name)
  lunio = io_open(trim(dir_path) // '/' // trim(name),status='old')
  call fem_conditions_read(lunio,p_mesh%f_mesh%npoin,f_cond)

  call par_mesh_to_triangulation (p_mesh, p_trian, f_cond)

  !write (*,*) '********** CREATE DOF HANDLER**************'
  call dhand%create( 1, 1, 1 )
  !call dof_handler_print ( dhand, 6 )


  call my_problem%create( p_trian%f_trian%num_dims )
  call my_approximation%create(my_problem)
  approximations(1)%p => my_approximation

  call dhand%set_problem( 1, my_approximation )
  ! ... for as many problems as we have

  call memalloc( p_trian%f_trian%num_elems, dhand%nvars_global, continuity, __FILE__, __LINE__)
  continuity = 1
  call memalloc( p_trian%f_trian%num_elems, dhand%nvars_global, order, __FILE__, __LINE__)
  order = 1
  call memalloc( p_trian%f_trian%num_elems, material, __FILE__, __LINE__)
  material = 1
  call memalloc( p_trian%f_trian%num_elems, problem, __FILE__, __LINE__)
  problem = 1
  
  call memalloc( p_trian%f_trian%num_elems, which_approx, __FILE__, __LINE__)
  which_approx = 1


  ! if ( context%iam > 0 ) then
  !    !pause
  !    do while ( 1 > 0)
  !       i = i + 1
  !    end do
  ! else
  !    write (*,*) 'Processor 0 not stopped'
  !    !i = 1
  !    !do while ( 1 > 0)
  !    !   i = i + 1
  !    !end do
  ! end if  


  ! Continuity
  !write(*,*) 'Continuity', continuity

  call par_fem_space_create ( p_trian, dhand, fspac, problem, approximations, &
                              f_cond, continuity, order, material, &
                              which_approx, num_approximations=1, time_steps_to_store = 1, &
                              hierarchical_basis = logical(.false.,lg), &
                              & static_condensation = logical(.false.,lg), num_continuity = 1 )


  call par_create_distributed_dof_info ( dhand, p_trian, fspac, dof_dist, dof_graph, (/ csr_symm /) )  

  call par_matrix_alloc ( csr_mat, symm_true, dof_graph(1,1), p_mat, positive_definite )
  call par_vector_alloc ( dof_dist(1), p_env, p_vec )
  call par_vector_alloc ( dof_dist(1), p_env, p_unk )

  call volume_integral( fspac, p_mat%f_matrix, p_vec%f_vector)

  !call dof_distribution_print ( 6, dof_dist(1) )

  call par_matrix_free (p_mat)
  call par_vector_free (p_vec)
  call par_vector_free (p_unk)

  call memfree( continuity, __FILE__, __LINE__)
  call memfree( order, __FILE__, __LINE__)
  call memfree( material, __FILE__, __LINE__)
  call memfree( problem, __FILE__, __LINE__)
  call memfree( which_approx, __FILE__, __LINE__)

  do i = 1, dhand%nblocks
     do j = 1, dhand%nblocks
        call par_graph_free( dof_graph(i,j) )
     end do
  end do
  call memfree ( dof_graph, __FILE__, __LINE__ )

  do iblock=1, dhand%nblocks
     call dof_distribution_free(dof_dist(iblock))
  end do
  deallocate(dof_dist)

  call fem_space_free(fspac) 
  call my_problem%free
  call dof_handler_free (dhand)
  call par_triangulation_free(p_trian)
  call fem_conditions_free (f_cond)
  call par_mesh_free (p_mesh)
  call par_environment_free (p_env)
  call par_context_free ( context )

  ! call memstatus

contains
  subroutine read_pars_cl (dir_path, prefix, dir_path_out)
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

  end subroutine read_pars_cl

end program par_test_cdr
