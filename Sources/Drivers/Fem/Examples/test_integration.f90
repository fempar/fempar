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
program test_integration
  !-----------------------------------------------------------------------
  ! Test assembly of a matrix
  !-----------------------------------------------------------------------
  use fem
  !use problems
  use nsi_names
  use nsi_cg_names
  type(fem_mesh)                          :: f_mesh
  type(fem_triangulation)                 :: f_trian
  type(fem_matrix)                        :: f_mat
  type(fem_conditions)                    :: f_cond

  integer(ip)        :: i, j, ierror, iblock
  type(dof_handler)  :: dhand
  type(fem_space)    :: my_fem_space

  type(fem_graph), allocatable :: my_dof_graph(:,:)
  type(fem_matrix)             :: my_matrix
  integer(ip)    , allocatable :: order(:,:), material(:), problem(:), which_approx(:), continuity(:,:)

  type(discrete_problem_pointer)  :: approximations(1)
  type(nsi_problem)               :: my_problem
  type(nsi_approximation), target :: my_aproximation

  character(len=:), allocatable  :: name 
  character(len=256)       :: dir_path, dir_path_out
  character(len=256)       :: prefix
  integer(ip)              :: lunio

  call  read_pars_cl ( dir_path, prefix, dir_path_out )

  !write (*,*) '********** READ MESH **********'

  call fem_mesh_compose_name ( prefix, name )
  lunio = io_open( trim(dir_path) // '/' // trim(name), 'read' )
  call fem_mesh_read ( lunio, f_mesh )
  call io_close(lunio)

  !write (*,*) '********** CREATE TRIANGULATION **********'

  call fem_triangulation_create ( f_mesh%nelem , f_trian )
  call mesh_to_triangulation    ( f_mesh       , f_trian )
  call fem_triangulation_to_dual( f_trian )

  !write (*,*) '********** CREATE DOF HANDLER **********'
  call dhand%create( 1, 1, 1 )
  !call dof_handler_print ( dhand, 6 )

  call memalloc( f_trian%num_elems, dhand%nvars_global, continuity, __FILE__, __LINE__)
  continuity = 1
  call memalloc( f_trian%num_elems, dhand%nvars_global, order, __FILE__, __LINE__)
  order = 3
  call memalloc( f_trian%num_elems, material, __FILE__, __LINE__)
  material = 1
  call memalloc( f_trian%num_elems, problem, __FILE__, __LINE__)
  problem = 1
  call memalloc( f_trian%num_elems, which_approx, __FILE__, __LINE__)
  which_approx = 1

  !write (*,*) '********** CREATE PROBLEM AND DISCRETIZATION **********'
  call my_problem%create(f_trian%num_dims)
  call my_aproximation%create(my_problem)
  approximations(1)%p => my_aproximation

  !write (*,*) '********** CREATE FEM SPACE AND DOF GRAPH **********'
  call fem_space_create ( f_trian, dhand, my_fem_space, problem, approximations, f_cond, continuity, order, material, &
       & which_approx, num_approximations=1, time_steps_to_store = 1, hierarchical_basis = logical(.false.,lg), &
       & static_condensation = logical(.false.,lg), num_continuity = 1 )
  call create_dof_info( dhand, f_trian, my_fem_space, my_dof_graph )
  call fem_matrix_alloc(csr_mat,symm_false, my_dof_graph(1,1), my_matrix)

  !write (*,*) '********** STARTING ASSEMBLY **********'
  call volume_integral(my_fem_space,my_matrix)
  !write (*,*) '********** FINISHED ASSEMBLY **********'


  !write (*,*) '********** FREE STRUCTURES **********'
  call fem_triangulation_free(f_trian)
  call fem_mesh_free ( f_mesh )

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


end program test_integration
