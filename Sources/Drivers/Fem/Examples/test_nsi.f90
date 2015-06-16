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
program test_nsi_iss
  use fem
  use nsi_names
  use nsi_cg_iss_names
  implicit none
# include "debug.i90"
  
  ! Types
  type(geom_data)                        :: gdata
  type(bound_data)                       :: bdata
  type(fem_fixed_info)                   :: ginfo
  type(fem_triangulation)                :: f_trian
  type(fem_conditions)                   :: f_cond
  type(dof_handler)                      :: dhand
  type(fem_space)                        :: fspac  
  type(nsi_problem)                      :: myprob
  type(nsi_cg_iss_approximation), target :: myapprox
  type(discrete_problem_pointer)         :: approximations(1)
  type(fem_matrix)              , target :: femat
  type(fem_vector)              , target :: fevec,feunk
  type(fem_precond)                      :: feprec
  type(fem_precond_params)               :: ppars
  type(solver_control)                   :: sctrl
  type(serial_environment)               :: senv
  class(base_operand)          , pointer :: x, b
  class(base_operator)         , pointer :: A, M

  ! Logicals
  logical(lg) :: ginfo_state

  ! Integers
  integer(ip) :: gtype(1) = (/ csr /)
  integer(ip) :: ibloc,jbloc,istat
  integer(ip) :: num_approximations = 1

  ! Allocatable
  integer(ip), allocatable :: continuity(:,:)
  integer(ip), allocatable :: order(:,:)
  integer(ip), allocatable :: material(:)
  integer(ip), allocatable :: problem(:)
  integer(ip), allocatable :: which_approx(:)
  type(fem_graph), pointer :: f_graph
  type(fem_block_graph)    :: f_blk_graph

  ! Arguments
  character(len=256) :: dir_path_out,prefix
  integer(ip)        :: nex,ney,nez
  
  call meminit

  ! Read parameters from command-line
  call read_pars_cl_test_nsi(prefix,dir_path_out,nex,ney,nez)

  ! Generate geometry data
  call geom_data_create(gdata,nex,ney,nez)

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

  ! Generate triangulation
  call gen_triangulation(1,gdata,bdata,ginfo,f_trian,f_cond,material)

  ! Create dof_handler
  call dhand%create(1,1,gdata%ndime+1)

  ! Create problem
  call myprob%create(gdata%ndime)
  call myapprox%create(myprob)
  approximations(1)%p => myapprox
  call dhand%set_problem(1,myapprox)
  myapprox%dtinv  = 0.0_rp
  myprob%kfl_conv = 1
  myprob%diffu    = 1.0_rp

  ! Allocate auxiliar elemental arrays
  call memalloc(f_trian%num_elems,dhand%nvars_global,continuity, __FILE__,__LINE__)
  call memalloc(f_trian%num_elems,dhand%nvars_global,order,__FILE__,__LINE__)
  call memalloc(f_trian%num_elems,problem,__FILE__,__LINE__)
  call memalloc(f_trian%num_elems,which_approx,__FILE__,__LINE__)
  continuity             = 1
  order(:,1:gdata%ndime) = 2
  order(:,gdata%ndime+1) = 1
  problem                = 1
  which_approx           = 1 
  
  ! Create fem_space
  call fem_space_create(f_trian,dhand,fspac,problem,num_approximations,approximations,f_cond, &
       &                continuity,order,material,which_approx,time_steps_to_store=3,         &
       &                hierarchical_basis=logical(.false.,lg),                               &
       &                static_condensation=logical(.false.,lg),num_continuity=1)

  ! Create dof info
  call create_dof_info(dhand,f_trian,fspac,f_blk_graph,gtype)

  f_graph => f_blk_graph%get_block(1,1)

  ! Allocate matrices and vectors
  call fem_matrix_alloc(csr_mat,symm_false,f_graph,femat)
  call fem_vector_alloc(f_graph%nv,fevec)
  call fem_vector_alloc(f_graph%nv,feunk)
  call fevec%init(0.0_rp)

  ! Update boundary conditions
  call update_strong_dirichlet_boundary_conditions(fspac)

  ! Integrate
  call volume_integral(fspac,femat,fevec)

  ! Construct preconditioner
  sctrl%method = direct
  ppars%type   = pardiso_mkl_prec
  call fem_precond_create(femat,feprec,ppars)
  call fem_precond_symbolic(femat,feprec)
  call fem_precond_numeric(femat,feprec)
  call fem_precond_log_info(feprec)

  ! Define operators
  A => femat
  b => fevec
  x => feunk

  ! Solve
  call abstract_solve(femat,feprec,fevec,feunk,sctrl,senv)
  !call abstract_solve(A,M,b,x,sctrl,senv)
  call solver_control_log_conv_his(sctrl)
  call solver_control_free_conv_his(sctrl)

  call fem_vector_print(6,feunk)

  ! Free preconditioner
  call fem_precond_free(precond_free_values,feprec)
  call fem_precond_free(precond_free_struct,feprec)
  call fem_precond_free(precond_free_clean,feprec)
  
  ! Deallocate
  call memfree(continuity,__FILE__,__LINE__)
  call memfree(order,__FILE__,__LINE__)
  call memfree(material,__FILE__,__LINE__)
  call memfree(problem,__FILE__,__LINE__)
  call memfree(which_approx,__FILE__,__LINE__)
  call f_blk_graph%free()
  call fem_vector_free(feunk)
  call fem_vector_free(fevec)
  call fem_matrix_free(femat) 
  call fem_space_free(fspac) 
  call myprob%free
  call myapprox%free
  call dof_handler_free(dhand)
  call fem_triangulation_free(f_trian)
  call fem_conditions_free(f_cond)
  call fem_element_fixed_info_free(ginfo)
  call bound_data_free(bdata)

  call memstatus

contains
  
  !==================================================================================================
  subroutine read_pars_cl_test_nsi(prefix,dir_path_out,nex,ney,nez)
    implicit none
    character*(*), intent(out) :: prefix, dir_path_out
    integer(ip)  , intent(out) :: nex,ney,nez
    character(len=256)         :: program_name
    character(len=256)         :: argument 
    integer                    :: numargs,iargc

    numargs = iargc()
    call getarg(0, program_name)
    if (.not. (numargs==5) ) then
       write (6,*) 'Usage: ', trim(program_name), ' prefix dir_path_out nex ney nez'
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

  end subroutine read_pars_cl_test_nsi

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
  
end program test_nsi_iss
