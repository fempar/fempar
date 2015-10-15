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

module mypart_names
  use serial_names
  implicit none
  
  type, extends(migratory_element_t) :: mypart_t
     integer(ip) :: mypart
   contains
     procedure :: size   => mypart_size
     procedure :: pack   => mypart_pack
     procedure :: unpack => mypart_unpack
  end type mypart_t

contains

  subroutine mypart_size (my, n)
    implicit none
    class(mypart_t), intent(in)  :: my
    integer(ip)    , intent(out) :: n
    ! Locals
    integer(ieep) :: mold(1)
    n  = size(transfer(1_ip ,mold))
  end subroutine mypart_size
  
  subroutine mypart_pack (my, n, buffer)
    implicit none
    class(mypart_t), intent(in)   :: my
    integer(ip)    , intent(in)   :: n
    integer(ieep)  , intent(out)  :: buffer(n)
    
    ! Locals
    integer(ieep) :: mold(1)
    integer(ip)   :: size_of_ip
    integer(ip)   :: start, end
    
    size_of_ip   = size(transfer(1_ip ,mold))
        
    start = 1
    end   = start + size_of_ip -1
    buffer(start:end) = transfer(my%mypart,mold)
  end subroutine mypart_pack
    
  subroutine mypart_unpack(my, n, buffer)
    implicit none
    class(mypart_t) , intent(inout)  :: my
    integer(ip)     , intent(in)     :: n
    integer(ieep)   , intent(in)     :: buffer(n)
        
    ! Locals
    integer(ieep) :: mold(1)
    integer(ip)   :: size_of_ip
    integer(ip)   :: start, end
        
    size_of_ip   = size(transfer(1_ip ,mold))
    
    start = 1
    end   = start + size_of_ip -1
    my%mypart  = transfer(buffer(start:end), my%mypart)
  end subroutine mypart_unpack 

end module mypart_names


program par_test_cdr_unstructured
  !----------------------------------------------------------
  ! Parallel partitioner test
  !----------------------------------------------------------
  use serial_names
  use par_names
  use prob_names
  use mypart_names
  
  implicit none
#include "debug.i90" 
  ! Our data
  type(par_context_t)       :: w_context, p_context, q_context, p_p_context, p_q_context, b_context
  type(par_environment_t)   :: p_env
  type(par_mesh_t)          :: p_mesh
  type(par_triangulation_t) :: p_trian
  type(par_fe_space_t)      :: p_fe_space

  class(matrix_t)          , pointer  :: matrix
  type(par_scalar_matrix_t), pointer  :: p_mat
  class(array_t)           , pointer  :: array
  type(par_scalar_array_t) , target   :: p_unk
  type(par_scalar_array_t) , pointer  :: p_vec
  type(fe_affine_operator_t)          :: fe_affine_operator
  class(vector_t) , pointer           :: x, y
  class(operator_t), pointer  :: A

  ! Preconditioner-related data structures
  type(par_preconditioner_dd_diagonal_t)           :: p_prec_dd_diag
  type(par_preconditioner_dd_mlevel_bddc_t), target :: p_mlevel_bddc
  ! type(par_preconditioner_dd_mlevel_bddc_t), pointer  :: point_to_p_mlevel_bddc
  type(par_preconditioner_dd_mlevel_bddc_params_t), target  :: p_mlevel_bddc_pars
  type(par_preconditioner_dd_mlevel_bddc_params_t), pointer :: point_to_p_mlevel_bddc_pars
  integer(ip), allocatable :: kind_coarse_dofs(:)

  type(solver_control_t)                :: sctrl
  type(dof_descriptor_t)                :: dof_descriptor
  logical                               :: symmetric_storage(1) = (/ .true. /)
  type(par_conditions_t)                :: p_cond

  type(cdr_problem_t)                   :: my_problem
  type(cdr_discrete_t)                  :: my_discrete
  type(cdr_approximation_t), target     :: my_approximation
  integer(ip)                           :: num_approximations
  type(p_discrete_integration_t)  :: approximations(1)

  integer(ip)              :: num_levels, nparts, ndime
  integer(ip), allocatable :: id_parts(:), num_parts(:)

  ! Arguments
  integer(ip)                   :: lunio
  character(len=256)            :: dir_path, dir_path_out
  character(len=256)            :: prefix
  character(len=:), allocatable :: name
  integer(ip)                   :: i, j, k, ierror, iblock, num_uniform_refinement_steps
  type(par_timer_t)             :: par_uniform_refinement_timer, par_mesh_to_triangulation_timer, par_fe_space_create_timer
  type(par_timer_t)             :: total_time


  integer(ip), allocatable :: order(:,:), material(:), problem(:)
  integer(ip), allocatable :: continuity(:,:)
  logical    , allocatable :: enable_face_integration(:,:)


  interface
     subroutine malloc_stats() bind(c,name='malloc_stats')
       use iso_c_binding
       implicit none
     end subroutine malloc_stats
  end interface

  call meminit

  ! Read parameters from command-line
  call read_pars_cl ( dir_path, prefix, dir_path_out, num_levels, num_parts, ndime, num_uniform_refinement_steps )

  call memalloc(num_levels, id_parts, __FILE__, __LINE__)
  call mesh_conditions_read_redistribute_on_first_level_tasks ( num_levels, num_parts, id_parts, &
                                                                w_context, p_context, q_context, b_context, p_env, &
                                                                p_mesh, p_cond )

  call par_timer_create ( par_mesh_to_triangulation_timer, 'PAR_MESH_TO_TRIANGULATION', w_context%icontxt )
  call par_timer_create ( par_fe_space_create_timer, 'PAR_FE_SPACE_CREATE', w_context%icontxt )
  call par_timer_create ( par_uniform_refinement_timer, 'PAR_UNIFORM_REFINEMENT', w_context%icontxt )

  do i=1, num_uniform_refinement_steps
     call par_timer_init (par_mesh_to_triangulation_timer)
     call par_timer_start (par_mesh_to_triangulation_timer)   
     call par_mesh_to_triangulation (p_mesh, p_trian, p_cond)
     call par_timer_stop (par_mesh_to_triangulation_timer)   
     call par_timer_report (par_mesh_to_triangulation_timer)   
     call par_mesh_free(p_mesh)

     call par_timer_init (par_uniform_refinement_timer)
     call par_timer_start (par_uniform_refinement_timer) 
     call par_uniform_refinement ( p_trian, p_mesh, p_cond )
     call par_timer_stop (par_uniform_refinement_timer)  
     call par_timer_report (par_uniform_refinement_timer)
     call par_triangulation_free(p_trian)

     if ( p_env%p_context%iam == 0 ) then
        call malloc_stats()
     end if
  end do
  call par_timer_init (par_mesh_to_triangulation_timer)
  call par_timer_start (par_mesh_to_triangulation_timer)   
  call par_mesh_to_triangulation (p_mesh, p_trian, p_cond)
  call par_timer_stop (par_mesh_to_triangulation_timer)   
  call par_timer_report (par_mesh_to_triangulation_timer) 
  
  if ( p_env%p_context%iam == 0 ) then
     call malloc_stats()
  end if

  !write (*,*) '********** CREATE DOF HANDLER**************'
  call dof_descriptor%create( 1, 1, 1 )


  call my_problem%create( p_trian%f_trian%num_dims )
  call my_discrete%create( my_problem )
  call my_approximation%create(my_problem,my_discrete)
  num_approximations=1
  approximations(1)%discrete_integration => my_approximation
  
  call dof_descriptor%set_problem( 1, my_discrete )
  ! ... for as many problems as we have

  call memalloc( p_trian%f_trian%num_elems, dof_descriptor%nvars_global, continuity, __FILE__, __LINE__)
  continuity = 1
  call memalloc( p_trian%f_trian%num_elems, dof_descriptor%nvars_global, enable_face_integration, __FILE__, __LINE__)
  enable_face_integration = .false.
  call memalloc( p_trian%f_trian%num_elems, dof_descriptor%nvars_global, order, __FILE__, __LINE__)
  order = 1
  call memalloc( p_trian%f_trian%num_elems, material, __FILE__, __LINE__)
  material = 1
  call memalloc( p_trian%f_trian%num_elems, problem, __FILE__, __LINE__)
  problem = 1

  ! Continuity
  ! write(*,*) 'Continuity', continuity
  call p_fe_space%create ( p_trian, dof_descriptor, problem, &
                              p_cond, continuity, enable_face_integration, order, material, &
                              time_steps_to_store = 1, &
                              hierarchical_basis = .false., &
                              & static_condensation = .false., num_continuity = 1 )
  
  call par_create_distributed_dof_info ( p_fe_space )

  ! if ( p_env%am_i_fine_task() ) p_cond%f_conditions%valu=1.0_rp
  call par_update_strong_dirichlet_bcond( p_fe_space, p_cond )

  call fe_affine_operator%create ( (/.true./), &
								   (/.true./), & 
								   (/positive_definite/), &
								   p_fe_space, &
								   approximations)
  
  call fe_affine_operator%symbolic_setup()
  call fe_affine_operator%numerical_setup()

  matrix => fe_affine_operator%get_matrix()
  select type(matrix)
  class is(par_scalar_matrix_t)
    p_mat => matrix
  class default
    check(.false.)
  end select 

  array => fe_affine_operator%get_array()
  select type(array)
  class is(par_scalar_array_t)
    p_vec => array
  class default
    check(.false.)
  end select 
  
  call p_unk%clone(p_vec)
  call p_unk%init(1.0_rp)

  A => p_mat
  x => p_vec
  y => p_unk
  y = x - A*y

  ! Define (recursive) parameters
  point_to_p_mlevel_bddc_pars => p_mlevel_bddc_pars
  do i=1, num_levels-1
     point_to_p_mlevel_bddc_pars%ndime            = ndime
     point_to_p_mlevel_bddc_pars%unknowns         = all_unknowns
     point_to_p_mlevel_bddc_pars%pad_collectives  = pad
     point_to_p_mlevel_bddc_pars%projection       = galerkin                           !default
     point_to_p_mlevel_bddc_pars%subd_elmat_calc  = phit_minus_c_i_t_lambda            !default  
     point_to_p_mlevel_bddc_pars%correction_mode  = additive                 !default 
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

  call par_timer_create ( total_time, 'Total time', w_context%icontxt )

  do j=2,ndime
     do k=1,4
        point_to_p_mlevel_bddc_pars => p_mlevel_bddc_pars
        do i=1, num_levels-1
           point_to_p_mlevel_bddc_pars%kind_coarse_dofs = kind_coarse_dofs(j)
           point_to_p_mlevel_bddc_pars => point_to_p_mlevel_bddc_pars%ppars_coarse_bddc
        end do

        call p_unk%init(0.0_rp)

        ! Create multilevel bddc inverse 
        call par_preconditioner_dd_mlevel_bddc_create( p_mat, p_mlevel_bddc, p_mlevel_bddc_pars )

        call par_timer_start (total_time) 

        ! Ass struct
        call par_preconditioner_dd_mlevel_bddc_ass_struct ( p_mat, p_mlevel_bddc )

        ! Fill val
        call par_preconditioner_dd_mlevel_bddc_fill_val ( p_mlevel_bddc )
        call par_preconditioner_dd_mlevel_bddc_static_condensation (p_mat, p_mlevel_bddc, p_vec, p_unk)

        call abstract_solve(p_mat,p_mlevel_bddc,p_vec,p_unk,sctrl,p_env)

        call par_timer_stop (total_time)   
        call par_timer_report (total_time) 

        if ( p_env%p_context%iam == 0 .or. p_mlevel_bddc%c_context%iam == 0 .or. p_mlevel_bddc%d_context%iam == 0  ) then
           call malloc_stats()
        end if

        ! Free bddc inverse
        call par_preconditioner_dd_mlevel_bddc_free_in_stages( p_mlevel_bddc, free_values)
        call par_preconditioner_dd_mlevel_bddc_free_in_stages( p_mlevel_bddc, free_struct)
        call par_preconditioner_dd_mlevel_bddc_free_in_stages( p_mlevel_bddc, free_clean)

     end do
  end do

  call memfree ( kind_coarse_dofs, __FILE__, __LINE__ )

  call p_unk%free()
  call fe_affine_operator%free()

  call memfree( continuity, __FILE__, __LINE__)
  call memfree( enable_face_integration, __FILE__, __LINE__)
  call memfree( order, __FILE__, __LINE__)
  call memfree( material, __FILE__, __LINE__)
  call memfree( problem, __FILE__, __LINE__)

  call p_fe_space%free() 
  call my_problem%free
  call my_discrete%free
  call my_approximation%free
  call dof_descriptor_free (dof_descriptor)
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

contains
  subroutine read_pars_cl (dir_path, prefix, dir_path_out, num_levels, num_parts, ndime, num_uniform_refinement_steps)
    implicit none
    character(*)              , intent(out) :: dir_path, prefix, dir_path_out
    integer(ip)               , intent(out) :: num_levels
    integer(ip)  , allocatable, intent(out) :: num_parts(:) 
    integer(ip)               , intent(out) :: ndime
    integer(ip)               , intent(out) :: num_uniform_refinement_steps



    character(len=256)           :: program_name
    character(len=256)           :: argument 
    integer                      :: numargs,iargc

    numargs = iargc()
    call getarg(0, program_name)
    if (.not. (numargs==8) ) then
       write (6,*) 'Usage: ', trim(program_name), ' dir_path prefix dir_path_out num_levels p1 p2 ndime num_uniform_refinement_steps'
       stop
    end if

    call getarg(1, argument)
    dir_path = trim(argument)

    call getarg(2, argument)
    prefix = trim(argument)

    call getarg(3,argument)
    dir_path_out = trim(argument)

    call getarg(4,argument)
    read(argument,*) num_levels
    assert(num_levels==3)
    
    call memalloc(num_levels, num_parts, __FILE__, __LINE__)

    call getarg(5,argument)
    read(argument,*) num_parts(1)

    call getarg(6,argument)
    read(argument,*) num_parts(2)
    num_parts(3) = 1

    call getarg(7,argument)
    read(argument,*) ndime

    call getarg(8,argument)
    read(argument,*) num_uniform_refinement_steps

  end subroutine read_pars_cl


  ! This subroutine reads par_mesh and par_conditions from the distributed file system on a subset of 1st level MPI tasks
  ! with size equivalent to the number of 2nd level MPI tasks. Then, it redistributes the data structures from this subset
  ! to the whole set of 1st level MPI tasks such that each of these get a single element. Besides, it sets the mapping governing
  ! aggregation among 1st level and 2nd level subdomains (i.e., id_parts), and the whole parallel environment (including parallel
  ! contexts)
  subroutine mesh_conditions_read_redistribute_on_first_level_tasks ( num_levels, num_parts, id_parts, &
                                                                      w_context, p_context, q_context, b_context, p_env, &
                                                                      p_mesh, p_cond )
    implicit none

    ! Dummy arguments
    integer(ip)                    , intent(in)  :: num_levels
    integer(ip)                    , intent(in)  :: num_parts(num_levels)
    integer(ip)                    , intent(out) :: id_parts(num_levels)
    type(par_context_t)            , intent(out) :: w_context, p_context, q_context, b_context
    type(par_environment_t), target, intent(out) :: p_env
    type(par_mesh_t)               , intent(out) :: p_mesh
    type(par_conditions_t)         , intent(out) :: p_cond

    ! Locals
    type(par_context_t) :: p_p_context, p_q_context

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

    id_parts(1) = p_context%iam+1
    id_parts(2) = q_context%iam+1
    id_parts(3) = 1

    if ( p_context%iam >= 0  ) then
       if ( p_context%iam < num_parts(2) ) then
          call par_context_create ( 1, p_p_context, p_q_context, p_context )
       else
          call par_context_create ( 2, p_q_context, p_p_context, p_context )
       end if
       
       if ( p_p_context%iam >= 0 ) then
          ! Create parallel environment
          call par_environment_create( p_env, p_p_context )
       
          ! Read mesh
          call par_mesh_read ( dir_path, prefix, p_env, p_mesh )
          
          ! Read boundary conditions
          call par_conditions_read(dir_path, prefix, p_mesh%f_mesh%npoin, p_env, p_cond)
       end if
       
       call scatter_mesh_and_conditions(num_levels, id_parts, p_context, p_p_context, p_q_context, p_mesh, p_cond)

       if ( p_p_context%iam >= 0 ) then
          ! Free parallel environmnet
          call par_environment_free ( p_env )
       end if
       
       call par_context_free ( p_p_context, .false. )
       call par_context_free ( p_q_context, .false. )
    end if

    ! Create parallel environment
    call par_environment_create( p_env,& 
                                 w_context,& 
                                 p_context,& 
                                 q_context,&
                                 b_context,&
                                 num_levels,&
                                 id_parts, & 
                                 num_parts )

    ! AFM: **IMPORTANT NOTE: This is very UGLY, and violates the OO encapsulation principle.
    ! I required to do such a thing as type(par_mesh_t) and type(par_conditions_t) do
    ! not provide sufficient expressiveness to re-set the corresponding 
    ! type(par_environment_t) instance.
    p_mesh%p_env => p_env
    p_cond%p_env => p_env

    
  end subroutine mesh_conditions_read_redistribute_on_first_level_tasks



  subroutine scatter_mesh_and_conditions(num_levels, id_parts, p_context, p_p_context, p_q_context, p_mesh, p_conditions)
    use mpi
    implicit none
    ! Dummy arguments
    integer(ip)           , intent(in)    :: num_levels
    integer(ip)           , intent(inout) :: id_parts(num_levels)
    type(par_context_t)   , intent(in)    :: p_context, p_p_context, p_q_context
    type(par_mesh_t)      , intent(inout) :: p_mesh
    type(par_conditions_t), intent(inout) :: p_conditions

    integer(ip)           , allocatable   :: num_elems_per_subdomain(:)
    integer(ip)           , allocatable   :: proc_map(:) 
    integer(ip)                           :: proc_map_local
    integer(ip)                           :: i, k, offset, from_proc
    integer                               :: ierr
    type(par_context_t)                   :: g_context, dum_context
    type(mesh_t)                          :: dual_mesh
    type(graph_t)                         :: dual_graph

    type(mesh_t), allocatable               :: lmesh(:)
    type(mesh_distribution_t), allocatable  :: ldist(:)
    type(conditions_t), allocatable         :: lconditions(:)
    integer(ip), allocatable                :: sizes_per_subdomain(:,:)
    integer(ip), allocatable                :: local_sizes(:)
    integer(ip), allocatable                :: sendcounts(:), displs(:)
    integer(ip)                             :: dummy_displs(1)
    integer(ip)                             :: dummy_sendcounts(1)
    integer(ieep)                           :: dummy_buffer_to_be_scattered(1)
    integer(ieep), allocatable              :: buffer_to_be_scattered(:)
    integer(ieep), allocatable              :: local_buffer(:)
    integer(ip)                             :: istat


    ! All gather number of elements per subdomain on p_p_context
    if ( p_p_context%iam >= 0 ) then
       call memalloc(p_p_context%np+1, num_elems_per_subdomain, __FILE__, __LINE__)
       num_elems_per_subdomain(1) = 1
       call mpi_allgather(p_mesh%f_mesh_dist%emap%nl, 1, MPI_INTEGER, &
                          num_elems_per_subdomain(2), 1, MPI_INTEGER, & 
                          p_p_context%icontxt, ierr)
       check(ierr==0)
       do i=1, p_p_context%np
          num_elems_per_subdomain(i+1)=num_elems_per_subdomain(i)+num_elems_per_subdomain(i+1) 
       end do
    end if

    if ( p_p_context%iam == 0 ) then
       call memalloc( num_elems_per_subdomain(p_p_context%np+1)-1, proc_map, __FILE__, __LINE__)
       do i=1,p_p_context%np
          proc_map(i)=i
          offset = p_p_context%np-i
          do j=num_elems_per_subdomain(i)+1, num_elems_per_subdomain(i+1)-1
             proc_map(j+offset)=i
          end do
       end do
       call mpi_scatter(proc_map      , 1, MPI_INTEGER, &
                        proc_map_local, 1, MPI_INTEGER, 0, p_context%icontxt, ierr)
       check(ierr==0)
       call memfree( proc_map, __FILE__, __LINE__)
    else
       call mpi_scatter(            -1, 1, MPI_INTEGER, &
                        proc_map_local, 1, MPI_INTEGER, 0, p_context%icontxt, ierr)
       check(ierr==0)
    end if


    ! Set 1st_to_2nd_level_mapping on 1st level MPI tasks
    id_parts(2) = proc_map_local !+ p_context%np

    if ( p_p_context%iam >= 0 ) then
       ! Generate dual mesh (i.e., list of elements around points)
       call mesh_to_dual(p_mesh%f_mesh, dual_mesh)

       ! Generate dual graph. Two elements will be neighbours it they share at least one vertex
       call create_graph_from_mesh(dual_mesh, p_mesh%f_mesh, 1, dual_graph)

       call split_mesh_into_one_element_per_subdomain( p_p_context, &
                                                       num_elems_per_subdomain, &
                                                       p_mesh%f_mesh, &
                                                       p_mesh%f_mesh_dist, & 
                                                       p_conditions%f_conditions, & 
                                                       dual_graph, & 
                                                       lmesh, & 
                                                       ldist, & 
                                                       lconditions )
 
       ! Deallocate (old) local portions of p_mesh, and p_conditions
       call mesh_free(p_mesh%f_mesh)
       call mesh_distribution_free(p_mesh%f_mesh_dist)
       call conditions_free(p_conditions%f_conditions)

       ! Deallocate temporary data
       call dual_graph%free()
       call mesh_free(dual_mesh)
       call memfree(num_elems_per_subdomain, __FILE__, __LINE__)
    end if

    call par_context_create ( proc_map_local, g_context, dum_context, p_context )

    if ( g_context%iam == 0 ) then
       assert ( g_context%np == size(lmesh) )
       call memalloc ( g_context%np  , sendcounts, __FILE__, __LINE__)
       call memalloc ( g_context%np+1,     displs, __FILE__, __LINE__)
       ! Store the size of lmesh, ldist, lconditions in a 2d array of size 3xnelem
       call memalloc ( 3, g_context%np, sizes_per_subdomain, __FILE__, __LINE__)
       displs(1)=0
       do i=1, size(lmesh)
          sizes_per_subdomain(1,i) = mesh_size(lmesh(i))
          sizes_per_subdomain(2,i) = mesh_distribution_size(ldist(i))
          sizes_per_subdomain(3,i) = conditions_size(lconditions(i))
          sendcounts(i) = sizes_per_subdomain(1,i) + sizes_per_subdomain(2,i) + sizes_per_subdomain(3,i)
          displs(i+1) = displs(i) + sendcounts(i)
       end do
       call memalloc ( 3, local_sizes, __FILE__, __LINE__)
       ! Scatter this 2d array to processors in g_context
       call mpi_scatter(sizes_per_subdomain, 3, MPI_INTEGER, &
                        local_sizes        , 3, MPI_INTEGER, 0, g_context%icontxt, ierr)
       check(ierr==0)

       ! Determine the size of the buffer to be scattered
       call memalloc ( displs(g_context%np+1)-1, buffer_to_be_scattered, __FILE__, __LINE__)
       
       ! Pack lmesh, ldist, lconditions into large buffer to be scattered
       k=1
       do i=1, size(lmesh)
          call mesh_pack ( lmesh(i), sizes_per_subdomain(1,i), buffer_to_be_scattered(k) )
          k = k + sizes_per_subdomain(1,i)
          call mesh_distribution_pack ( ldist(i), sizes_per_subdomain(2,i), buffer_to_be_scattered(k) )
          k = k + sizes_per_subdomain(2,i)
          call conditions_pack ( lconditions(i), sizes_per_subdomain(3,i), buffer_to_be_scattered(k) )
          k = k + sizes_per_subdomain(3,i)
       end do

       ! Free lmesh, ldist, lconditions
       do i=1, size(lmesh)
          call mesh_free( lmesh(i) )
          call mesh_distribution_free ( ldist(i) )
          call conditions_free ( lconditions(i) )
       end do
       
       deallocate(lmesh, stat=istat)
       check(istat==0)

       deallocate(ldist, stat=istat)
       check(istat==0)
       
       deallocate(lconditions, stat=istat)
       check(istat==0)

       ! Scatter large buffer
       call memalloc ( local_sizes(1)+local_sizes(2)+local_sizes(3), local_buffer, __FILE__, __LINE__)
       call mpi_scatterv(buffer_to_be_scattered, sendcounts, displs,&
                         MPI_INTEGER1, local_buffer, size(local_buffer), MPI_INTEGER1,&
                         0, g_context%icontxt, ierr)
       check(ierr==0)

       call memfree ( sendcounts, __FILE__, __LINE__)
       call memfree ( displs    , __FILE__, __LINE__)
       call memfree ( sizes_per_subdomain, __FILE__, __LINE__)
       call memfree ( buffer_to_be_scattered, __FILE__, __LINE__)

    else
       call memalloc ( 3, local_sizes, __FILE__, __LINE__)
       call mpi_scatter(         -1, 3, MPI_INTEGER, &
                        local_sizes, 3, MPI_INTEGER, 0, g_context%icontxt, ierr)
       check(ierr==0)
       
       ! Scatter large buffer
       call memalloc ( local_sizes(1)+local_sizes(2)+local_sizes(3), local_buffer, __FILE__, __LINE__)
       call mpi_scatterv(dummy_buffer_to_be_scattered, dummy_sendcounts, dummy_displs,&
                         MPI_INTEGER1, local_buffer, size(local_buffer), MPI_INTEGER1,&
                         0, g_context%icontxt, ierr)
       check(ierr==0)
    end if

    k=1
    call mesh_unpack ( p_mesh%f_mesh, local_sizes(1), local_buffer(k) )
    k = k + local_sizes(1)
    call mesh_distribution_unpack ( p_mesh%f_mesh_dist, local_sizes(2), local_buffer(k) )
    k = k + local_sizes(2)
    call conditions_unpack ( p_cond%f_conditions, local_sizes(3), local_buffer(k) )

    ! call mesh_distribution_print ( 6, p_mesh%f_mesh_dist )

    call memfree ( local_sizes, __FILE__, __LINE__)
    call memfree ( local_buffer, __FILE__, __LINE__)
    call par_context_free ( g_context, .false. )

  end subroutine scatter_mesh_and_conditions


  subroutine split_mesh_into_one_element_per_subdomain ( p_p_context, &
                                                         num_elems_per_subdomain, &
                                                         mesh, &
                                                         mesh_dist, &
                                                         conditions, &
                                                         dual_graph, &
                                                         lmesh, &
                                                         ldist, &
                                                         lconditions ) 
    implicit none

    ! Dummy arguments
    type(par_context_t)                   , intent(in)  :: p_p_context
    integer(ip)                           , intent(in)  :: num_elems_per_subdomain(p_p_context%np+1)
    type(mesh_t)                          , intent(in)  :: mesh
    type(mesh_distribution_t)             , intent(in)  :: mesh_dist
    type(conditions_t)                    , intent(in)  :: conditions
    type(graph_t)                         , intent(in)  :: dual_graph
    type(mesh_t)             , allocatable, intent(out) :: lmesh(:)
    type(mesh_distribution_t), allocatable, intent(out) :: ldist(:)
    type(conditions_t)       , allocatable, intent(out) :: lconditions(:)

    ! Local variables
    integer(ip)                  :: istat, ielem, i, j, k, l
    type(element_import_t)       :: elem_import
    type(mypart_t), allocatable  :: data(:) 
    type(hash_table_ip_ip_t)     :: is_boundary_elem
    type(hash_table_igp_ip_t)    :: ghost_g2l
    type(mesh_distribution_t)    :: dum_mesh_distribution 


    allocate(lmesh(mesh%nelem), stat=istat)
    check(istat==0)

    allocate(ldist(mesh%nelem), stat=istat)
    check(istat==0)

    allocate(lconditions(mesh%nelem), stat=istat)
    check(istat==0)

    call element_import_create ( mesh_dist, elem_import )

    ! call element_import_print ( 6, elem_import ) 

    allocate ( data(elem_import%nelem + elem_import%nghost), stat=istat )
    check(istat==0)

    data(1)%mypart = p_p_context%iam + 1 
    do ielem=2, mesh%nelem
       data(ielem)%mypart = num_elems_per_subdomain(p_p_context%iam+1)+(ielem-1)+(p_p_context%np-(p_p_context%iam+1))
    end do

    call ghost_elements_exchange( p_p_context%icontxt, elem_import, data )

    call is_boundary_elem%init( max(5,int(sqrt(real(mesh%nelem,rp))*0.1_rp,ip)) )
    do i = 1, mesh_dist%nebou
       call is_boundary_elem%put(key=mesh_dist%lebou(i), val=i, stat=istat)
       assert(istat==now_stored)
    end do

    call ghost_g2l%init( max(5,int(sqrt(real(mesh%nelem,rp))*0.1_rp,ip)) )
    assert ( elem_import%nghost == size(elem_import%rcv_geids) )
    do i = 1, elem_import%nghost
       call ghost_g2l%put(key=elem_import%rcv_geids(i), val=i+mesh%nelem, stat=istat)
       assert(istat==now_stored)
    end do

    do ielem = 1, mesh%nelem
       ! Generate mesh
       lmesh(ielem)%nelem = 1
       lmesh(ielem)%nnode = mesh%pnods(ielem+1) - mesh%pnods(ielem)
       lmesh(ielem)%ndime = mesh%ndime
       lmesh(ielem)%npoin = lmesh(ielem)%nnode 

       call memalloc ( 2, lmesh(ielem)%pnods, __FILE__, __LINE__ )
       lmesh(ielem)%pnods(1) = 1
       lmesh(ielem)%pnods(2) = 1 + lmesh(ielem)%nnode

       call memalloc ( lmesh(ielem)%nnode, lmesh(ielem)%lnods, __FILE__, __LINE__ )
       do j=1, lmesh(ielem)%nnode
          lmesh(ielem)%lnods(j)=j
       end do

       call memalloc (mesh%ndime, lmesh(ielem)%nnode, lmesh(ielem)%coord, __FILE__, __LINE__ )

       i=1
       do j=mesh%pnods(ielem),mesh%pnods(ielem+1)-1 
          lmesh(ielem)%coord(:,i) = mesh%coord(:,mesh%lnods(j))
          i=i+1
       end do

       ! Generate conditions
       call conditions_create ( conditions%ncode, conditions%nvalu, lmesh(ielem)%nnode, lconditions(ielem))
       i=1
       do j=mesh%pnods(ielem),mesh%pnods(ielem+1)-1 
          lconditions(ielem)%code(:,i) = conditions%code(:,mesh%lnods(j))
          lconditions(ielem)%valu(:,i) = conditions%valu(:,mesh%lnods(j))
          i=i+1
       end do

       ! Generate mesh_distribution
       ldist(ielem)%ipart  = data(ielem)%mypart
       ldist(ielem)%nparts = num_elems_per_subdomain(p_p_context%np+1)-1
       ldist(ielem)%nebou  = 1
       ldist(ielem)%nnbou  = lmesh(ielem)%nnode
       call memalloc ( ldist(ielem)%nebou, ldist(ielem)%lebou, __FILE__, __LINE__ )
       ldist(ielem)%lebou = 1 

       call memalloc ( ldist(ielem)%nnbou, ldist(ielem)%lnbou, __FILE__, __LINE__ )
       ldist(ielem)%lnbou = lmesh(ielem)%lnods

       call map_alloc(1,mesh_dist%emap%ng,ldist(ielem)%emap)
       ldist(ielem)%emap%l2g(1) = mesh_dist%emap%l2g(ielem)

       call map_alloc(ldist(ielem)%nnbou, mesh_dist%nmap%ng, ldist(ielem)%nmap)
       i=1
       do j=mesh%pnods(ielem), mesh%pnods(ielem+1)-1 
          ldist(ielem)%nmap%l2g(i) = mesh_dist%nmap%l2g(mesh%lnods(j))
          i=i+1
       end do

       call memalloc (2, ldist(ielem)%pextn, __FILE__, __LINE__)
       ldist(ielem)%pextn(1)=1

       call is_boundary_elem%get(key=ielem, val=i, stat=istat)
       if ( istat==key_found ) then ! Boundary element
          ldist(ielem)%pextn(2) = ldist(ielem)%pextn(1) + & 
               (dual_graph%ia(ielem+1)-dual_graph%ia(ielem)) + & 
               (mesh_dist%pextn(i+1)-mesh_dist%pextn(i))

          call memalloc ( ldist(ielem)%pextn(2)-1, ldist(ielem)%lextn, __FILE__, __LINE__)
          call memalloc ( ldist(ielem)%pextn(2)-1, ldist(ielem)%lextp, __FILE__, __LINE__)

          k=1
          do j=dual_graph%ia(ielem), dual_graph%ia(ielem+1)-1
             ldist(ielem)%lextn(k) = mesh_dist%emap%l2g(dual_graph%ja(j))
             ldist(ielem)%lextp(k) = data(dual_graph%ja(j))%mypart
             k=k+1
          end do

          do j=mesh_dist%pextn(i), mesh_dist%pextn(i+1)-1
             ldist(ielem)%lextn(k) = mesh_dist%lextn(j)
             call ghost_g2l%get(key=mesh_dist%lextn(j), val=l, stat=istat)
             assert (istat==key_found)
             assert ( mesh_dist%lextn(j) == elem_import%rcv_geids(l-mesh%nelem) )
             ldist(ielem)%lextp(k) = data(l)%mypart
             k=k+1
          end do
       else ! Interior element
          ldist(ielem)%pextn(2) = ldist(ielem)%pextn(1) + (dual_graph%ia(ielem+1)-dual_graph%ia(ielem))
          call memalloc ( ldist(ielem)%pextn(2)-1, ldist(ielem)%lextn, __FILE__, __LINE__)
          call memalloc ( ldist(ielem)%pextn(2)-1, ldist(ielem)%lextp, __FILE__, __LINE__)
          k=1
          do j=dual_graph%ia(ielem), dual_graph%ia(ielem+1)-1
             ldist(ielem)%lextn(k) = mesh_dist%emap%l2g(dual_graph%ja(j))
             ldist(ielem)%lextp(k) = data(dual_graph%ja(j))%mypart
             k=k+1
          end do
       end if
       ! call mesh_distribution_print (6, ldist(ielem))
    end do
    call is_boundary_elem%free()
    call ghost_g2l%free()
    call element_import_free ( elem_import )
    deallocate ( data, stat=istat )
    check(istat==0)

  end subroutine split_mesh_into_one_element_per_subdomain

    function mesh_size (mesh)
      implicit none
      type(mesh_t), intent(in) :: mesh
      integer(ip) :: mesh_size

      ! Locals
      integer(ieep) :: mold(1)
      integer(ip)   :: size_of_ip, size_of_rp
      
      size_of_ip   = size(transfer(1_ip ,mold))
      size_of_rp   = size(transfer(1.0_rp,mold))
      
      mesh_size = size_of_ip*4 + &  
                  size_of_ip * size(mesh%pnods) + &
                  size_of_ip * size(mesh%lnods) + &
                  size_of_rp * size(mesh%coord,1) * size(mesh%coord,2) 
    end function mesh_size

    subroutine mesh_pack (mesh, n, buffer)
      implicit none
      type(mesh_t)  , intent(in)  :: mesh
      integer(ip)   , intent(in)  :: n
      integer(ieep) , intent(out) :: buffer(n)

      ! Locals
      integer(ieep) :: mold(1)
      integer(ip)   :: size_of_ip, size_of_igp, size_of_rp
      integer(ip)   :: start, end
      
      size_of_ip   = size(transfer(1_ip ,mold))
      size_of_igp  = size(transfer(1_igp,mold))
      size_of_rp   = size(transfer(1.0_rp,mold))

      start = 1
      end   = start + size_of_ip - 1
      buffer(start:end) = transfer(mesh%nelem,mold)

      start = end + 1
      end   = start + size_of_ip - 1
      buffer(start:end) = transfer(mesh%nnode,mold)

      start = end + 1
      end   = start + size_of_ip - 1
      buffer(start:end) = transfer(mesh%ndime,mold)

      start = end + 1
      end   = start + size_of_ip - 1
      buffer(start:end) = transfer(mesh%npoin,mold)

      start = end + 1
      end   = start + size(mesh%pnods)*size_of_ip - 1
      buffer(start:end) = transfer(mesh%pnods,mold)

      start = end + 1
      end   = start + size(mesh%lnods)*size_of_ip - 1
      buffer(start:end) = transfer(mesh%lnods,mold)

      start = end + 1
      end   = start + size(mesh%coord,1)*size(mesh%coord,2)*size_of_rp - 1
      buffer(start:end) = transfer(mesh%coord,mold)
    end subroutine mesh_pack

    subroutine mesh_unpack (mesh, n, buffer)
      implicit none
      type(mesh_t)              , intent(inout)  :: mesh
      integer(ip)               , intent(in)  :: n
      integer(ieep)             , intent(in) :: buffer(n)

      ! Locals
      integer(ieep)         :: mold(1)
      integer(ip)           :: size_of_ip, size_of_igp, size_of_rp
      integer(ip)           :: start, end
      real(rp), allocatable :: aux(:)
      integer(ip)           :: i,j,k
      
      size_of_ip   = size(transfer(1_ip ,mold))
      size_of_igp  = size(transfer(1_igp,mold))
      size_of_rp   = size(transfer(1.0_rp,mold))

      start = 1
      end   = start + size_of_ip - 1
      mesh%nelem = transfer(buffer(start:end),mesh%nelem)

      start = end + 1
      end   = start + size_of_ip - 1
      mesh%nnode = transfer(buffer(start:end),mesh%nnode)

      start = end + 1
      end   = start + size_of_ip - 1
      mesh%ndime = transfer(buffer(start:end),mesh%ndime)

      start = end + 1
      end   = start + size_of_ip - 1
      mesh%npoin = transfer(buffer(start:end),mesh%npoin)

      call memalloc(mesh%nelem+1, mesh%pnods,__FILE__,__LINE__)
      start = end + 1
      end   = start + size(mesh%pnods)*size_of_ip - 1
      mesh%pnods = transfer(buffer(start:end),mesh%pnods)

      call memalloc(mesh%pnods(mesh%nelem+1)-1, mesh%lnods,__FILE__,__LINE__)
      start = end + 1
      end   = start + size(mesh%lnods)*size_of_ip - 1
      mesh%lnods = transfer(buffer(start:end),mesh%lnods)

      call memalloc(mesh%ndime*mesh%npoin, aux,__FILE__,__LINE__)
      start = end + 1
      end   = start + size(aux)*size_of_rp - 1
      aux   = transfer(buffer(start:end),aux)

      call memalloc(mesh%ndime,mesh%npoin, mesh%coord,__FILE__,__LINE__)
      k=1
      do j=1,mesh%npoin
       do i=1,mesh%ndime
         mesh%coord(i,j) = aux(k)
         k=k+1
       end do
      end do
      call memfree(aux,__FILE__,__LINE__)
         
    end subroutine mesh_unpack

    function mesh_distribution_size (mesh_distribution)
      implicit none
      type(mesh_distribution_t), intent(in) :: mesh_distribution
      integer(ip) :: mesh_distribution_size

      ! Locals
      integer(ieep) :: mold(1)
      integer(ip)   :: size_of_ip, size_of_igp
      
      size_of_ip   = size(transfer(1_ip ,mold))
      size_of_igp  = size(transfer(1_igp,mold))
      
      mesh_distribution_size = size_of_ip*4 + &  
                               size_of_ip  * size(mesh_distribution%pextn) + &
                               size_of_ip  * size(mesh_distribution%lextp) + &
                               size_of_igp * size(mesh_distribution%lextn) + &
                               size_of_ip  * size(mesh_distribution%lebou) + &
                               size_of_ip  * size(mesh_distribution%lnbou) + &
                               size_of_igp + & ! emap%ng
                               size_of_ip  + & ! emap%nl
                               size_of_igp * size(mesh_distribution%emap%l2g) + &
                               size_of_igp + & ! nmap%ng
                               size_of_ip  + & ! nmap%nl
                               size_of_igp * size(mesh_distribution%nmap%l2g)

    end function mesh_distribution_size

    subroutine mesh_distribution_pack (mesh_distribution, n, buffer)
      implicit none
      type(mesh_distribution_t) , intent(in)  :: mesh_distribution
      integer(ip)               , intent(in)  :: n
      integer(ieep)             , intent(out) :: buffer(n)

      ! Locals
      integer(ieep) :: mold(1)
      integer(ip)   :: size_of_ip, size_of_igp
      integer(ip)   :: start, end
      
      size_of_ip   = size(transfer(1_ip ,mold))
      size_of_igp  = size(transfer(1_igp,mold))

      start = 1
      end   = start + size_of_ip - 1
      buffer(start:end) = transfer(mesh_distribution%ipart,mold)

      start = end + 1
      end   = start + size_of_ip - 1
      buffer(start:end) = transfer(mesh_distribution%nparts,mold)

      start = end + 1
      end   = start + size_of_ip - 1
      buffer(start:end) = transfer(mesh_distribution%nebou,mold)

      start = end + 1
      end   = start + size_of_ip - 1
      buffer(start:end) = transfer(mesh_distribution%nnbou,mold)

      start = end + 1
      end   = start + size(mesh_distribution%pextn)*size_of_ip - 1
      buffer(start:end) = transfer(mesh_distribution%pextn,mold)

      start = end + 1
      end   = start + size(mesh_distribution%lextp)*size_of_ip - 1
      buffer(start:end) = transfer(mesh_distribution%lextp,mold)

      start = end + 1
      end   = start + size(mesh_distribution%lextn)*size_of_igp - 1
      buffer(start:end) = transfer(mesh_distribution%lextn,mold)

      start = end + 1
      end   = start + size(mesh_distribution%lebou)*size_of_ip - 1
      buffer(start:end) = transfer(mesh_distribution%lebou,mold)

      start = end + 1
      end   = start + size(mesh_distribution%lnbou)*size_of_ip - 1
      buffer(start:end) = transfer(mesh_distribution%lnbou,mold)

      start = end + 1
      end   = start + size_of_igp - 1
      buffer(start:end) = transfer(mesh_distribution%emap%ng,mold)

      start = end + 1
      end   = start + size_of_ip - 1
      buffer(start:end) = transfer(mesh_distribution%emap%nl,mold)
      
      start = end + 1
      end   = start + size(mesh_distribution%emap%l2g)*size_of_igp - 1
      buffer(start:end) = transfer(mesh_distribution%emap%l2g,mold)

      start = end + 1
      end   = start + size_of_igp - 1
      buffer(start:end) = transfer(mesh_distribution%nmap%ng,mold)
 
      start = end + 1
      end   = start + size_of_ip - 1
      buffer(start:end) = transfer(mesh_distribution%nmap%nl,mold)
      
      start = end + 1
      end   = start + size(mesh_distribution%nmap%l2g)*size_of_igp - 1
      buffer(start:end) = transfer(mesh_distribution%nmap%l2g,mold)
    end subroutine mesh_distribution_pack

    subroutine mesh_distribution_unpack (mesh_distribution, n, buffer)
      implicit none
      type(mesh_distribution_t) , intent(inout)  :: mesh_distribution
      integer(ip)               , intent(in)     :: n
      integer(ieep)             , intent(in)     :: buffer(n)

      ! Locals
      integer(ieep) :: mold(1)
      integer(ip)   :: size_of_ip, size_of_igp
      integer(ip)   :: start, end
      
      size_of_ip   = size(transfer(1_ip ,mold))
      size_of_igp  = size(transfer(1_igp,mold))

      start = 1
      end   = start + size_of_ip - 1
      mesh_distribution%ipart = transfer(buffer(start:end),mesh_distribution%ipart)

      start = end + 1
      end   = start + size_of_ip - 1
      mesh_distribution%nparts = transfer(buffer(start:end),mesh_distribution%nparts)

      start = end + 1
      end   = start + size_of_ip - 1
      mesh_distribution%nebou = transfer(buffer(start:end),mesh_distribution%nebou)

      start = end + 1
      end   = start + size_of_ip - 1
      mesh_distribution%nnbou = transfer(buffer(start:end),mesh_distribution%nnbou)

      call memalloc( mesh_distribution%nebou+1, mesh_distribution%pextn, __FILE__, __LINE__)
      start = end + 1
      end   = start + size(mesh_distribution%pextn)*size_of_ip - 1
      mesh_distribution%pextn = transfer(buffer(start:end),mesh_distribution%pextn)

      call memalloc( mesh_distribution%pextn(mesh_distribution%nebou+1)-1, & 
                     mesh_distribution%lextp, __FILE__, __LINE__)
      start = end + 1
      end   = start + size(mesh_distribution%lextp)*size_of_ip - 1
      mesh_distribution%lextp = transfer(buffer(start:end),mesh_distribution%lextp)

      call memalloc( mesh_distribution%pextn(mesh_distribution%nebou+1)-1, & 
                     mesh_distribution%lextn, __FILE__, __LINE__)
      start = end + 1
      end   = start + size(mesh_distribution%lextn)*size_of_igp - 1
      mesh_distribution%lextn = transfer(buffer(start:end),mesh_distribution%lextn)

      call memalloc( mesh_distribution%nebou, mesh_distribution%lebou, __FILE__, __LINE__)
      start = end + 1
      end   = start + size(mesh_distribution%lebou)*size_of_ip - 1
      mesh_distribution%lebou = transfer(buffer(start:end),mesh_distribution%lebou)

      call memalloc( mesh_distribution%nnbou, mesh_distribution%lnbou, __FILE__, __LINE__)
      start = end + 1
      end   = start + size(mesh_distribution%lnbou)*size_of_ip - 1
      mesh_distribution%lnbou = transfer(buffer(start:end),mesh_distribution%lnbou)

      start = end + 1
      end   = start + size_of_igp - 1
      mesh_distribution%emap%ng = transfer(buffer(start:end),mesh_distribution%emap%ng)

      start = end + 1
      end   = start + size_of_ip - 1
      mesh_distribution%emap%nl = transfer(buffer(start:end),mesh_distribution%emap%nl)
 
      call memalloc( mesh_distribution%emap%nl, mesh_distribution%emap%l2g, __FILE__, __LINE__)
      start = end + 1
      end   = start + size_of_igp*size(mesh_distribution%emap%l2g) - 1
      mesh_distribution%emap%l2g = transfer(buffer(start:end),mesh_distribution%emap%l2g)

      start = end + 1
      end   = start + size_of_igp - 1
      mesh_distribution%nmap%ng = transfer(buffer(start:end),mesh_distribution%nmap%ng)

      start = end + 1
      end   = start + size_of_ip - 1
      mesh_distribution%nmap%nl = transfer(buffer(start:end),mesh_distribution%nmap%nl)
     
      call memalloc( mesh_distribution%nmap%nl, mesh_distribution%nmap%l2g, __FILE__, __LINE__)
      start = end + 1
      end   = start + size(mesh_distribution%nmap%l2g)*size_of_igp - 1
      mesh_distribution%nmap%l2g = transfer(buffer(start:end),mesh_distribution%nmap%l2g)

    end subroutine mesh_distribution_unpack

    function conditions_size (conditions)
      implicit none
      type(conditions_t), intent(in) :: conditions
      integer(ip) :: conditions_size

      ! Locals
      integer(ieep) :: mold(1)
      integer(ip)   :: size_of_ip, size_of_rp
      
      size_of_ip = size(transfer(1_ip ,mold))
      size_of_rp = size(transfer(1.0_rp,mold))

      conditions_size = size_of_ip*3 + &  
                        size_of_ip  * size(conditions%code,1) * size(conditions%code,2) + &
                        size_of_rp  * size(conditions%valu,1) * size(conditions%valu,2) 
    end function conditions_size


    subroutine conditions_pack (conditions, n, buffer)
      implicit none
      type(conditions_t) , intent(in)  :: conditions
      integer(ip)        , intent(in)  :: n
      integer(ieep)      , intent(out) :: buffer(n)

      ! Locals
      integer(ieep) :: mold(1)
      integer(ip)   :: size_of_ip, size_of_rp
      integer(ip)   :: start, end
      
      size_of_ip   = size(transfer(1_ip ,mold))
      size_of_rp = size(transfer(1.0_rp,mold))

      start = 1
      end   = start + size_of_ip - 1
      buffer(start:end) = transfer(conditions%ncode,mold)

      start = end + 1
      end   = start + size_of_ip - 1
      buffer(start:end) = transfer(conditions%nvalu,mold)

      start = end + 1
      end   = start + size_of_ip - 1
      buffer(start:end) = transfer(conditions%ncond,mold)

      start = end + 1
      end   = start + size(conditions%code,1)*size(conditions%code,2)*size_of_ip - 1
      buffer(start:end) = transfer(conditions%code,mold)

      start = end + 1
      end   = start + size(conditions%valu,1)*size(conditions%valu,2)*size_of_rp - 1
      buffer(start:end) = transfer(conditions%valu,mold)

    end subroutine conditions_pack

    subroutine conditions_unpack (conditions, n, buffer)
      implicit none
      type(conditions_t) , intent(inout)  :: conditions
      integer(ip)        , intent(in)     :: n
      integer(ieep)      , intent(in)     :: buffer(n)

      ! Locals
      integer(ieep) :: mold(1)
      integer(ip)   :: size_of_ip, size_of_rp
      integer(ip)   :: start, end, i, j, k
      real(rp), allocatable :: aux_rp(:)
      integer(ip), allocatable :: aux_ip(:)
      
      size_of_ip   = size(transfer(1_ip ,mold))
      size_of_rp = size(transfer(1.0_rp,mold))

      start = 1
      end   = start + size_of_ip - 1
      conditions%ncode = transfer(buffer(start:end),conditions%ncode)

      start = end + 1
      end   = start + size_of_ip - 1
      conditions%nvalu = transfer(buffer(start:end),conditions%nvalu)

      start = end + 1
      end   = start + size_of_ip - 1
      conditions%ncond = transfer(buffer(start:end),conditions%ncond)

      call memalloc(conditions%ncode*conditions%ncond, aux_ip,__FILE__,__LINE__)
      start    = end + 1
      end      = start + size(aux_ip)*size_of_ip - 1
      aux_ip   = transfer(buffer(start:end),aux_ip)

      call memalloc(conditions%ncode, conditions%ncond, conditions%code,__FILE__,__LINE__)
      k=1
      do j=1,conditions%ncond
       do i=1,conditions%ncode
         conditions%code(i,j) = aux_ip(k)
         k=k+1
       end do
      end do
      call memfree(aux_ip,__FILE__,__LINE__)

      call memalloc(conditions%nvalu*conditions%ncond, aux_rp,__FILE__,__LINE__)
      start  = end + 1
      end    = start + size(aux_rp)*size_of_rp - 1
      aux_rp = transfer(buffer(start:end),aux_rp)
      call memalloc(conditions%nvalu, conditions%ncond, conditions%valu,__FILE__,__LINE__)
      k=1
      do j=1,conditions%ncond
       do i=1,conditions%nvalu
         conditions%valu(i,j) = aux_rp(k)
         k=k+1
       end do
      end do
      call memfree(aux_rp,__FILE__,__LINE__)

    end subroutine conditions_unpack

end program par_test_cdr_unstructured
