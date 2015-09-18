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
module postprocess_field_nsi_names
  use types_names
  use problem_names
  use postprocess_field_names
# include "debug.i90"
  implicit none

  type, extends(postprocess_field_t) :: postprocess_field_velocity_t
   contains
     procedure :: compute_field => nsi_compute_field_vel
  end type postprocess_field_velocity_t

  type, extends(postprocess_field_t) :: postprocess_field_pressure_t
   contains
     procedure :: compute_field => nsi_compute_field_pre
  end type postprocess_field_pressure_t

  public :: postprocess_field_velocity_t, postprocess_field_pressure_t
contains
  subroutine nsi_compute_field_vel(postprocess_field)
    !******** READ ME: *********************************************************!
    ! This subroutine nsi_compute_field has been created only for testing       !
    ! purposes and it is left here as an easy example of how to use the         !
    ! postprocess_fields_names module.                                          !
    !---------------------------------------------------------------------------!
    implicit none
    class(postprocess_field_velocity_t), intent(inout) :: postprocess_field
    ! Locals
    integer(ip)        :: nelem, ndime, i

    ! Extract some parameter form fe_space
    nelem = postprocess_field%fe_space%g_trian%num_elems
    ndime = postprocess_field%fe_space%g_trian%num_dims

    do i=1,nelem
       postprocess_field%fe_postprocess_field(i)%nodal_properties =                              &
            &          postprocess_field%fe_space%finite_elements(i)%unkno(:,1:ndime,1)
    end do
    

  end subroutine nsi_compute_field_vel

  subroutine nsi_compute_field_pre(postprocess_field)
    !******** READ ME: *********************************************************!
    ! This subroutine nsi_compute_field has been created only for testing       !
    ! purposes and it is left here as an easy example of how to use the         !
    ! postprocess_fields_names module.                                          !
    !---------------------------------------------------------------------------!
    implicit none
    class(postprocess_field_pressure_t), intent(inout) :: postprocess_field
    ! Locals
    integer(ip)        :: nelem, ndime, i

    ! Extract some parameter form fe_space
    nelem = postprocess_field%fe_space%g_trian%num_elems
    ndime = postprocess_field%fe_space%g_trian%num_dims

    do i=1,nelem
       postprocess_field%fe_postprocess_field(i)%nodal_properties =                              &
            &          postprocess_field%fe_space%finite_elements(i)%unkno(:,1+ndime:1+ndime,1)
    end do
    
  end subroutine nsi_compute_field_pre
end module postprocess_field_nsi_names

program par_test_nsi_iss
  use serial_names
  use par_names
  use nsi_names
  use nsi_cg_iss_names
  use lib_vtk_io_interface_names
  use postprocess_field_nsi_names
  implicit none
# include "debug.i90"

  ! Types
  type(uniform_mesh_descriptor_t)                    :: gdata
  type(uniform_conditions_descriptor_t)              :: bdata
  type(reference_element_t)                          :: geo_reference_element
  type(par_context_t)                                :: w_context
  type(par_context_t)                                :: p_context
  type(par_context_t)                                :: q_context
  type(par_context_t)                                :: b_context
  type(par_environment_t)                            :: p_env
  type(par_triangulation_t)                          :: p_trian
  type(par_conditions_t)                             :: p_cond
  type(dof_descriptor_t)                             :: dof_descriptor
  type(par_fe_space_t)                               :: p_fe_space  
  type(nsi_problem_t)                                :: myprob
  type(nsi_cg_iss_discrete_t)               , target :: mydisc
  type(nsi_cg_iss_matvec_t)                 , target :: cg_iss_matvec
  type(discrete_integration_pointer_t)               :: approx(1)
  type(vtk_t)                                        :: fevtk
  type(par_block_graph_t)                            :: p_blk_graph
  type(block_dof_distribution_t)                     :: blk_dof_dist
  type(par_preconditioner_dd_mlevel_bddc_t)       , target  :: p_mlevel_bddc
  type(par_preconditioner_dd_mlevel_bddc_params_t), target  :: p_mlevel_bddc_pars
  type(par_preconditioner_dd_identity_t)                    :: p_prec_dd_diag
  type(par_preconditioner_dd_mlevel_bddc_params_t), pointer :: point_to_p_mlevel_bddc_pars 
  type(par_matrix_t), target                         :: p_mat
  type(par_vector_t), target                         :: p_vec
  type(par_vector_t), target                         :: p_unk
  type(solver_control_t)                             :: sctrl
  class(base_operand_t) , pointer           :: x, y
  class(base_operator_t), pointer           :: A
  type(postprocess_field_velocity_t)    :: postprocess_vel
  type(postprocess_field_pressure_t)    :: postprocess_pre

  ! Integers
  integer(ip) :: num_levels
  integer(ip) :: gtype(1) = (/ csr /)
  integer(ip) :: ibloc,jbloc,istat,i,j
  integer(ip) :: num_approximations = 1

  ! Allocatables
  integer(ip), allocatable :: id_parts(:)
  integer(ip), allocatable :: num_parts(:)
  integer(ip), allocatable :: continuity(:,:)
  integer(ip), allocatable :: face_coupling(:,:)
  integer(ip), allocatable :: order(:,:)
  integer(ip), allocatable :: material(:)
  integer(ip), allocatable :: problem(:)
  integer(ip), allocatable :: which_approx(:)

  ! Arguments
  character(len=256) :: dir_path_out,prefix
  integer(ip)        :: nex,ney,nez,npx,npy,npz

  call meminit

  ! Read parameters from command-line
  call read_pars_cl_par_test_nsi(prefix,dir_path_out,nex,ney,nez,npx,npy,npz)

  ! Generate geometry data
  call uniform_mesh_descriptor_create(gdata,nex,ney,nez,npx,npy,npz)

  ! Generate boundary data
  call uniform_conditions_descriptor_create(gdata%ndime+1,gdata%ndime+1,gdata%ndime,bdata)
  bdata%poin%code(gdata%ndime+1,1:2**gdata%ndime-1) = 0
  bdata%line%code(gdata%ndime+1,:) = 0
  bdata%surf%code(gdata%ndime+1,:) = 0
  bdata%poin%code(gdata%ndime+1,2**gdata%ndime) = 1
  bdata%poin%valu(gdata%ndime+1,2**gdata%ndime) = 0.0_rp
  bdata%poin%valu(1:gdata%ndime,:) = 1.0_rp
  bdata%line%valu(1:gdata%ndime,:) = 1.0_rp

  ! Generate element geometrical fixed info
  call reference_element_create(geo_reference_element,Q_type_id,1,gdata%ndime)

  ! Set levels
  num_levels = 2
  call memalloc(num_levels, id_parts , __FILE__, __LINE__)
  call memalloc(num_levels, num_parts, __FILE__, __LINE__)
  num_parts = (/gdata%nparts, 1/)
  id_parts = (/w_context%iam+1, 1/)

  ! Start parallel execution
  call par_context_create(w_context)

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

  ! Generate par triangulation
  call par_generate_uniform_triangulation(p_env,gdata,bdata,geo_reference_element,p_trian,p_cond,material)

  ! Create dof_descriptor
  call dof_descriptor%create(1,1,gdata%ndime+1)

  ! Create problem
  call myprob%create(gdata%ndime)
  call mydisc%create(myprob)
  call cg_iss_matvec%create(myprob,mydisc)
  call dof_descriptor%set_problem(1,mydisc)
  approx(1)%p       => cg_iss_matvec
  mydisc%dtinv      = 0.0_rp
  myprob%kfl_conv   = 1
  myprob%diffu      = 1.0_rp
  myprob%case_veloc = 1
  myprob%case_press = 1

  ! Allocate auxiliar elemental arrays
  call memalloc(p_trian%f_trian%num_elems,dof_descriptor%nvars_global,continuity, __FILE__,__LINE__)
  call memalloc(p_trian%f_trian%num_elems,dof_descriptor%nvars_global,face_coupling, __FILE__,__LINE__)
  call memalloc(p_trian%f_trian%num_elems,dof_descriptor%nvars_global,order,__FILE__,__LINE__)
  call memalloc(p_trian%f_trian%num_elems,problem,__FILE__,__LINE__)
  call memalloc(p_trian%f_trian%num_elems,which_approx,__FILE__,__LINE__)
  continuity             = 1
  face_coupling             = 1
  order(:,1:gdata%ndime) = 2
  order(:,gdata%ndime+1) = 1
  problem                = 1
  which_approx           = 1 

  ! Create par_fe_space
  call par_fe_space_create(p_trian,dof_descriptor,p_fe_space,problem,p_cond,continuity,face_coupling, &
       &                    order,material,which_approx,time_steps_to_store=3,                        &
       &                    hierarchical_basis=.false.,                         &
       &                    static_condensation=.false.,num_continuity=1)

  ! Initialize VTK output
  call fevtk%initialize(p_trian%f_trian,p_fe_space%fe_space,myprob,p_env,dir_path_out,prefix, &
       &                nparts=gdata%nparts)!,linear_order=.false.)

  ! Create dof info
  call par_create_distributed_dof_info(dof_descriptor,p_trian,p_fe_space,blk_dof_dist,p_blk_graph,gtype)  

  !if(p_env%am_i_fine_task()) call par_graph_print(6,p_blk_graph%get_block(1,1))

  ! Allocate matrices and vectors
  call par_matrix_alloc(csr_mat,symm_false,p_blk_graph%get_block(1,1),p_mat)
  call par_vector_alloc(blk_dof_dist%get_block(1),p_env,p_vec)
  call par_vector_alloc(blk_dof_dist%get_block(1),p_env,p_unk)
  p_vec%state = part_summed
  p_unk%state = full_summed
  call p_vec%init(0.0_rp)

  !call p_unk%init(1.0_rp)
  !call par_vector_weight(p_unk)
  !call p_unk%comm()
  !call par_vector_print(6,p_unk)
  
  ! Apply boundary conditions to unkno
  if ( p_env%am_i_fine_task() ) p_cond%f_conditions%valu = 1.0_rp
  call par_update_strong_dirichlet_bcond(p_fe_space,p_cond)
  !call par_update_analytical_bcond((/1:gdata%ndime/),myprob%case_veloc,0.0_rp,p_fe_space
  call par_update_analytical_bcond((/(i, i=1,gdata%ndime)/),myprob%case_veloc,0.0_rp,p_fe_space)
  call par_update_analytical_bcond((/gdata%ndime+1/),myprob%case_press,0.0_rp,p_fe_space)

  ! Integrate
  if(p_env%am_i_fine_task()) then
     call volume_integral(approx,p_fe_space%fe_space,p_mat%f_matrix,p_vec%f_vector)
  end if

  ! Define (recursive) parameters
  point_to_p_mlevel_bddc_pars => p_mlevel_bddc_pars
  do i=1, num_levels-1
     point_to_p_mlevel_bddc_pars%ndime            = gdata%ndime
     point_to_p_mlevel_bddc_pars%unknowns         = all_unknowns
     point_to_p_mlevel_bddc_pars%pad_collectives  = pad
     point_to_p_mlevel_bddc_pars%projection       = petrov_galerkin                     
     point_to_p_mlevel_bddc_pars%subd_elmat_calc  = phit_minus_c_i_t_lambda            !default  
     point_to_p_mlevel_bddc_pars%correction_mode  = additive_symmetric                 !default 
     point_to_p_mlevel_bddc_pars%nn_sys_sol_strat = corners_rest_part_solve_expl_schur ! default 
     if(gdata%ndime==3) then
        point_to_p_mlevel_bddc_pars%kind_coarse_dofs = corners_edges_and_faces
     else
        point_to_p_mlevel_bddc_pars%kind_coarse_dofs = corners_and_edges
     end if
     if ( i < num_levels-1 ) then
        point_to_p_mlevel_bddc_pars%co_sys_sol_strat     = recursive_bddc
        point_to_p_mlevel_bddc_pars%ppars_harm%type      = pardiso_mkl_prec !umfpack_prec 
        point_to_p_mlevel_bddc_pars%ppars_dirichlet%type = pardiso_mkl_prec !umfpack_prec   
        if ( i == 1 ) then
           point_to_p_mlevel_bddc_pars%spars_coarse%method = direct
           point_to_p_mlevel_bddc_pars%spars_coarse%itmax  = 200
           point_to_p_mlevel_bddc_pars%spars_coarse%rtol   = 1.0e-20
           point_to_p_mlevel_bddc_pars%spars_coarse%trace  = 1
           point_to_p_mlevel_bddc_pars%correction_mode     = additive
        end if
        allocate(point_to_p_mlevel_bddc_pars%ppars_coarse_bddc, stat = istat)
        check(istat==0)
        point_to_p_mlevel_bddc_pars => point_to_p_mlevel_bddc_pars%ppars_coarse_bddc
     else
        point_to_p_mlevel_bddc_pars%co_sys_sol_strat         = serial_gather
        point_to_p_mlevel_bddc_pars%ppars_harm%type          = pardiso_mkl_prec !umfpack_prec  
        point_to_p_mlevel_bddc_pars%ppars_dirichlet%type     = pardiso_mkl_prec !umfpack_prec  
        point_to_p_mlevel_bddc_pars%ppars_coarse_serial%type = pardiso_mkl_prec !umfpack_prec  
        nullify ( point_to_p_mlevel_bddc_pars%ppars_coarse_bddc )
     end if
  end do
  
  point_to_p_mlevel_bddc_pars => p_mlevel_bddc_pars
  do i=1, num_levels-1
     point_to_p_mlevel_bddc_pars => point_to_p_mlevel_bddc_pars%ppars_coarse_bddc
  end do

  ! Define solver control parameters
  sctrl%method=rgmres
  sctrl%trace=100
  sctrl%itmax=800
  sctrl%dkrymax=800
  sctrl%stopc=res_nrmgiven_res_nrmgiven
  sctrl%orto=icgs
  sctrl%rtol=1.0e-14_rp

!!$  call p_unk%init(1.0_rp)
!!$  A => p_mat
!!$  x => p_vec
!!$  y => p_unk
!!$  y = x - A*y
!!$  write(*,*) 'XXX error norm XXX', y%nrm2()
!!$  p_unk%state = full_summed

  ! Create Preconditioner 
  call par_preconditioner_dd_mlevel_bddc_create(p_mat,p_mlevel_bddc,p_mlevel_bddc_pars)
  call par_preconditioner_dd_mlevel_bddc_ass_struct(p_mat,p_mlevel_bddc)
  call par_preconditioner_dd_mlevel_bddc_fill_val(p_mat,p_mlevel_bddc)

!!$  call par_preconditioner_dd_identity_create ( p_mat, p_prec_dd_diag )
!!$  call par_preconditioner_dd_identity_ass_struct ( p_mat, p_prec_dd_diag )
!!$  call par_preconditioner_dd_identity_fill_val ( p_mat, p_prec_dd_diag )
!!$
!!$  !call par_vector_print(6,p_vec)
!!$  call abstract_solve(p_mat,p_prec_dd_diag,p_vec,p_unk,sctrl,p_env)
!!$  call par_vector_print(6,p_unk)
!!$
!!$  call par_preconditioner_dd_identity_free ( p_prec_dd_diag, free_values )
!!$  call par_preconditioner_dd_identity_free ( p_prec_dd_diag, free_struct )
!!$  call par_preconditioner_dd_identity_free ( p_prec_dd_diag, free_clean )

  ! Solve
  call abstract_solve(p_mat,p_mlevel_bddc,p_vec,p_unk,sctrl,p_env)
  call solver_control_log_conv_his(sctrl)
  call solver_control_free_conv_his(sctrl)  
  !call par_vector_print(6,p_unk)

  ! Store solution to unkno
  call par_update_solution(p_unk,p_fe_space)

  ! Compute postprocess field
  call postprocess_vel%create('velocity',gdata%ndime,p_fe_space%fe_space,p_env,                     &
       &                      use_interpolation_order_from_variable,variable_identifier=1)
  call postprocess_pre%create('pressure',1,p_fe_space%fe_space,p_env,                               &
       &                      use_interpolation_order_from_variable,variable_identifier=2)
  call postprocess_vel%compute_and_finalize_field()
  call postprocess_pre%compute_and_finalize_field()
     
  ! Print solution to VTK file
  istat = fevtk%write_VTK_start(n_part=p_env%p_context%iam,o_fmt='ascii')
  istat = fevtk%write_VTK_unknowns()
  istat = fevtk%write_VTK_field(postprocess_vel)
  istat = fevtk%write_VTK_field(postprocess_pre)
  istat = fevtk%write_VTK_end() 
  istat = fevtk%write_PVTK()
  call postprocess_vel%free
  call postprocess_pre%free
  
  ! Free preconditioner
  call par_preconditioner_dd_mlevel_bddc_free(p_mlevel_bddc,free_values)
  call par_preconditioner_dd_mlevel_bddc_free(p_mlevel_bddc,free_struct)
  call par_preconditioner_dd_mlevel_bddc_free(p_mlevel_bddc,free_clean)

  ! Deallocate
  call memfree(id_parts , __FILE__, __LINE__)
  call memfree(num_parts, __FILE__, __LINE__)
  call memfree(continuity,__FILE__,__LINE__)
  call memfree(face_coupling,__FILE__,__LINE__)
  call memfree(order,__FILE__,__LINE__)
  call memfree(material,__FILE__,__LINE__)
  call memfree(problem,__FILE__,__LINE__)
  call memfree(which_approx,__FILE__,__LINE__)
  call fevtk%free
  call par_matrix_free (p_mat)
  call par_vector_free (p_vec)
  call par_vector_free (p_unk)
  call p_blk_graph%free
  call blk_dof_dist%free
  call par_fe_space_free(p_fe_space) 
  call myprob%free
  call mydisc%free
  call cg_iss_matvec%free
  call dof_descriptor_free (dof_descriptor)
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
  subroutine read_pars_cl_par_test_nsi(prefix,dir_path_out,nex,ney,nez,npx,npy,npz)
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

  end subroutine read_pars_cl_par_test_nsi

end program par_test_nsi_iss
