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
  type(fem_mesh)                          :: f_mesh
  type(fem_triangulation)                 :: f_trian
  type(fem_matrix)                        :: f_mat
  type(fem_conditions)                    :: f_cond
  type(dof_handler)  :: dhand
  type(fem_space)    :: fspac
  type(fem_graph), allocatable    :: dof_graph(:,:)

  type(cdr_problem)               :: my_problem
  type(cdr_approximation), target :: my_approximation
  type(discrete_problem_pointer)  :: approximations(1)

  type(fem_matrix)             :: my_matrix

  ! Arguments
  character(len=256)       :: dir_path, dir_path_out
  character(len=256)       :: prefix, filename
  integer(ip)              :: i, j, vars_prob(1) = 1, ierror, iblock

  integer(ip), allocatable :: order(:,:), material(:), problem(:), which_approx(:)

  integer(ip), allocatable :: continuity(:,:)

  integer(ip) :: lunio

  call meminit

  ! Read parameters from command-line
  call  read_pars_cl_test_cdr ( dir_path, prefix, dir_path_out )

  ! Read mesh
  filename = trim(dir_path)//'data/'//trim(prefix)//'.msh'
  lunio = io_open(filename,status='old')
  call fem_mesh_read(lunio,f_mesh,permute_c2z=.true.)
  call io_close(lunio)

  ! Read conditions 
  filename = trim(dir_path)//'data/'//trim(prefix)//'.cnd'
  lunio = io_open(filename,status='old')
  call fem_conditions_read(lunio,f_mesh%npoin,f_cond)
  call io_close(lunio)

  !call fem_mesh_write ( 6, f_mesh, 'square_uniform' )

  call mesh_to_triangulation ( f_mesh, f_trian, gcond = f_cond )

  write(*,*) 'conditions%ncode', f_cond%ncode
  write(*,*) 'conditions%nvalu', f_cond%nvalu
  write(*,*) 'conditions%ncond', f_cond%ncond

  write(*,*) 'conditions%code', f_cond%code
  write(*,*) 'conditions%valu', f_cond%valu

  call triangulation_print( 6 , f_trian )

  vars_prob = 1
  call dhand%create( 1, 1, 1 ) !, vars_block, dof_coupl )
  !                      ( dhand, nblocks, nprobs, nvars_global, vars_block, dof_coupl )
  !call dof_handler_print ( dhand, 6 )

  call my_problem%create( f_trian%num_dims )
  call dhand%set_problem( 1, my_problem )
  !                     ( ndime, dhand, l2g_vars, iprob ) 

  ! ... for as many problems as we have

  call my_approximation%create(my_problem)
  approximations(1)%p => my_approximation


  call memalloc( f_trian%num_elems, dhand%nvars_global, continuity, __FILE__, __LINE__)
  continuity = 1
  call memalloc( f_trian%num_elems, dhand%nvars_global, order, __FILE__, __LINE__)
  order = 4
  call memalloc( f_trian%num_elems, material, __FILE__, __LINE__)
  material = 1
  call memalloc( f_trian%num_elems, problem, __FILE__, __LINE__)
  problem = 1
  call memalloc( f_trian%num_elems, which_approx, __FILE__, __LINE__)
  which_approx = 1














  ! Continuity
  !write(*,*) 'Continuity', continuity

  call fem_space_create ( f_trian, dhand, fspac, problem, approximations, f_cond, continuity, order, material, &
       & which_approx, num_approximations=1, time_steps_to_store = 1, hierarchical_basis = logical(.false.,lg), & 
       & static_condensation = logical(.false.,lg), num_continuity = 1 )

  call create_dof_info( dhand, f_trian, fspac, dof_graph )

  call fem_space_print( 6, fspac )


  call fem_matrix_alloc(csr_mat,symm_false, dof_graph(1,1), my_matrix)

  !write (*,*) '********** STARTING ASSEMBLY **********'
  call volume_integral(fspac,my_matrix)
  !write (*,*) '********** FINISHED ASSEMBLY **********'






  ! call fem_precond_dd_mlevel_bddc_create ( f_mat, mlbddc, mlbddc_params )

  call memfree( continuity, __FILE__, __LINE__)
  call memfree( order, __FILE__, __LINE__)
  call memfree( material, __FILE__, __LINE__)
  call memfree( problem, __FILE__, __LINE__)
  call memfree( which_approx, __FILE__, __LINE__)

  do i = 1, dhand%nblocks
     do j = 1, dhand%nblocks
        call fem_graph_free( dof_graph(i,j) )
     end do
  end do
     

  call fem_space_free(fspac) 

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

end program test_cdr
