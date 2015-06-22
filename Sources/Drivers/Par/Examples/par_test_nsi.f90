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
program par_test_nsi_iss
  use fem
  use par 
  use nsi_names
  use nsi_cg_iss_names
  use lib_vtk_io_interface
  implicit none
# include "debug.i90"

  ! Types
  type(geom_data)                                  :: gdata
  type(bound_data)                                 :: bdata
  type(fem_fixed_info)                             :: ginfo
  type(par_context)                                :: w_context
  type(par_context)                                :: p_context
  type(par_context)                                :: q_context
  type(par_context)                                :: b_context
  type(par_environment)                            :: p_env
  type(par_triangulation)                          :: p_trian
  type(par_conditions)                             :: p_cond
  type(dof_handler)                                :: dhand
  type(par_fem_space)                              :: p_fspac  
  type(nsi_problem)                                :: myprob
  type(nsi_cg_iss_discrete)               , target :: mydisc
  type(nsi_cg_iss_matvec)                 , target :: matvec
  type(discrete_integration_pointer)               :: approx(1)
  type(fem_vtk)                                    :: fevtk
  type(par_block_graph)                            :: p_blk_graph
  type(block_dof_distribution)                     :: blk_dof_dist
  type(par_precond_dd_mlevel_bddc)       , target  :: p_mlevel_bddc
  type(par_precond_dd_mlevel_bddc_params), target  :: p_mlevel_bddc_pars
  type(par_precond_dd_mlevel_bddc_params), pointer :: point_to_p_mlevel_bddc_pars 
  type(par_matrix), target                         :: p_mat
  type(par_vector), target                         :: p_vec
  type(par_vector), target                         :: p_unk
  type(solver_control)                             :: sctrl

  ! Logicals
  logical(lg) :: ginfo_state

  ! Integers
  integer(ip) :: num_levels
  integer(ip) :: gtype(1) = (/ csr /)
  integer(ip) :: ibloc,jbloc,istat,i,j
  integer(ip) :: num_approximations = 1

  ! Allocatables
  integer(ip), allocatable :: id_parts(:)
  integer(ip), allocatable :: num_parts(:)
  integer(ip), allocatable :: continuity(:,:)
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
  call geom_data_create(gdata,nex,ney,nez,npx,npy,npz)

  ! Generate boundary data
  call bound_data_create(gdata%ndime+1,gdata%ndime+1,gdata%ndime,bdata)
  bdata%poin%code(gdata%ndime+1,1:2**gdata%ndime-1) = 0
  bdata%line%code(gdata%ndime+1,:) = 0
  bdata%surf%code(gdata%ndime+1,:) = 0
  bdata%poin%code(gdata%ndime+1,2**gdata%ndime) = 1
  bdata%poin%valu(gdata%ndime+1,2**gdata%ndime) = 0.0_rp
  bdata%poin%valu(1:gdata%ndime,:) = 1.0_rp
  bdata%line%valu(1:gdata%ndime,:) = 1.0_rp

  ! Generate element geometrical fixed info
  call fem_element_fixed_info_create(ginfo,Q_type_id,1,gdata%ndime,ginfo_state)

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
  call par_gen_triangulation(p_env,gdata,bdata,ginfo,p_trian,p_cond,material)

  ! Create dof_handler
  call dhand%create(1,1,gdata%ndime+1)

  ! Create problem
  call myprob%create(gdata%ndime)
  call mydisc%create(myprob)
  call matvec%create(myprob,mydisc)
  call dhand%set_problem(1,mydisc)
  approx(1)%p     => matvec
  mydisc%dtinv    = 0.0_rp
  myprob%kfl_conv = 1
  myprob%diffu    = 1.0_rp

  ! Allocate auxiliar elemental arrays
  call memalloc(p_trian%f_trian%num_elems,dhand%nvars_global,continuity, __FILE__,__LINE__)
  call memalloc(p_trian%f_trian%num_elems,dhand%nvars_global,order,__FILE__,__LINE__)
  call memalloc(p_trian%f_trian%num_elems,problem,__FILE__,__LINE__)
  call memalloc(p_trian%f_trian%num_elems,which_approx,__FILE__,__LINE__)
  continuity             = 1
  order(:,1:gdata%ndime) = 2
  order(:,gdata%ndime+1) = 1
  problem                = 1
  which_approx           = 1 

  ! Create par_fem_space
  call par_fem_space_create(p_trian,dhand,p_fspac,problem,p_cond,continuity,order,material, &
       &                    which_approx,time_steps_to_store=3,                             &
       &                    hierarchical_basis=logical(.false.,lg),                         &
       &                    static_condensation=logical(.false.,lg),num_continuity=1)

  ! Initialize VTK output
  if(p_env%am_i_fine_task()) then
     call fevtk%initialize(p_trian%f_trian,p_fspac%f_space,myprob,dir_path_out,prefix, &
          &                nparts=gdata%nparts)
  end if

  ! Create dof info
  call par_create_distributed_dof_info(dhand,p_trian,p_fspac,blk_dof_dist,p_blk_graph,gtype)  

  ! Allocate matrices and vectors
  call par_matrix_alloc(csr_mat,symm_false,p_blk_graph%get_block(1,1),p_mat)
  call par_vector_alloc(blk_dof_dist%get_block(1),p_env,p_vec)
  call par_vector_alloc(blk_dof_dist%get_block(1),p_env,p_unk)
  p_vec%state = part_summed
  p_unk%state = full_summed
  call p_vec%init(0.0_rp)
  
  ! Apply boundary conditions to unkno
  if(p_env%am_i_fine_task()) then
     call update_strong_dirichlet_boundary_conditions(p_fspac%f_space)
  end if

  ! Integrate
  if(p_env%am_i_fine_task()) then
     call volume_integral(approx,p_fspac%f_space,p_mat%f_matrix,p_vec%f_vector)
  end if

  ! Define (recursive) parameters
  point_to_p_mlevel_bddc_pars => p_mlevel_bddc_pars
  do i=1, num_levels-1
     point_to_p_mlevel_bddc_pars%ndime            = gdata%ndime
     point_to_p_mlevel_bddc_pars%unknowns         = all_unknowns
     point_to_p_mlevel_bddc_pars%pad_collectives  = pad
     point_to_p_mlevel_bddc_pars%projection       = galerkin                           !default
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
  sctrl%stopc=res_nrmgiven_rhs_nrmgiven
  sctrl%orto=icgs
  sctrl%rtol=1.0e-20

  ! Create Preconditioner 
  call par_precond_dd_mlevel_bddc_create(p_mat,p_mlevel_bddc,p_mlevel_bddc_pars)
  call par_precond_dd_mlevel_bddc_ass_struct(p_mat,p_mlevel_bddc)
  call par_precond_dd_mlevel_bddc_fill_val(p_mat,p_mlevel_bddc)

  ! Solve
  call abstract_solve(p_mat,p_mlevel_bddc,p_vec,p_unk,sctrl,p_env)
  call solver_control_log_conv_his(sctrl)
  call solver_control_free_conv_his(sctrl)  
  if(p_env%am_i_fine_task()) call fem_vector_print(6,p_unk%f_vector)

  ! Print solution to VTK file
  if(p_env%am_i_fine_task()) istat = fevtk%write_VTK(n_part=p_env%p_context%iam,o_fmt='ascii')
  if(p_env%p_context%iam==0) istat = fevtk%write_PVTK()

  ! Free preconditioner
  call par_precond_dd_mlevel_bddc_free(p_mlevel_bddc,free_only_values)
  call par_precond_dd_mlevel_bddc_free(p_mlevel_bddc,free_only_struct)
  call par_precond_dd_mlevel_bddc_free(p_mlevel_bddc,free_clean)

  ! Deallocate
  call memfree(id_parts , __FILE__, __LINE__)
  call memfree(num_parts, __FILE__, __LINE__)
  call memfree(continuity,__FILE__,__LINE__)
  call memfree(order,__FILE__,__LINE__)
  call memfree(material,__FILE__,__LINE__)
  call memfree(problem,__FILE__,__LINE__)
  call memfree(which_approx,__FILE__,__LINE__)
  if(p_env%am_i_fine_task()) call fevtk%free
  call par_matrix_free (p_mat)
  call par_vector_free (p_vec)
  call par_vector_free (p_unk)
  call p_blk_graph%free
  call blk_dof_dist%free
  call par_fem_space_free(p_fspac) 
  call myprob%free
  call mydisc%free
  call matvec%free
  call dof_handler_free (dhand)
  call par_triangulation_free(p_trian)
  call par_conditions_free (p_cond)
  call fem_element_fixed_info_free(ginfo)
  call bound_data_free(bdata)
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

subroutine update_strong_dirichlet_boundary_conditions( fspac )
    implicit none

    type(fem_space), intent(inout)    :: fspac

    integer(ip) :: ielem, iobje, ivar, inode, l_node

    do ielem = 1, fspac%g_trian%num_elems
       do ivar=1, fspac%dof_handler%problems(problem(ielem))%p%nvars-1
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

  end subroutine update_strong_dirichlet_boundary_conditions

end program par_test_nsi_iss
