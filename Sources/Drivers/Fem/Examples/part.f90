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
program partitioner
  !-----------------------------------------------------------------------
  ! A finite element preprocessor
  !-----------------------------------------------------------------------
use fem_names
  !use mesh_graph_partition
  implicit none
  integer(ip)                       :: power, type, use_graph, nparts, nstr
  type(mesh_size_t)                   :: msize
  type(fem_conditions_t)              :: poin,line,surf
  type(fem_conditions_t)              :: gnodes, gbouns
  type(fem_mesh_t)                    :: gmesh
  type(fem_materials_t)               :: gmat
  type(fem_mesh_distribution_t) , allocatable :: distr(:)
  type(fem_mesh_t)      , allocatable :: lmesh(:)
  type(fem_conditions_t), allocatable :: lnodes(:),lbouns(:)
  type(fem_materials_t),  allocatable :: lmater(:)
  type(post_file_t)                   :: lupos
  integer(ip)                       :: ipart,idof,ndofn,ndime,prob,gid_sq
  integer(ip), parameter            :: cdr=1, sto=2, mhd=3
  integer(ip), parameter            :: gid=1, square=2
  integer(ip)                       :: isper(3),nedir(3)

  character(len=256)                :: dir_path, dir_path_out
  character(len=256)                :: name, prefix
  integer(ip), allocatable          :: ldome(:)
  integer(ip)                       :: i, j

  type(part_params_t) :: prt_pars

  call meminit

  call read_mesh_part_pars_cl(nparts,dir_path,prefix,dir_path_out)
  
  ! Read mesh
  call fem_mesh_read(dir_path, prefix, gmesh, permute_c2z=.true.)

  ! Read conditions
  call fem_conditions_read(dir_path, prefix, gmesh%npoin, gnodes)

  ! Read materials
  call fem_materials_read(dir_path, prefix, gmat)
  gmat%list = 1 
  
  ! Write original mesh for postprocess
  call fem_mesh_write(dir_path_out, prefix, gmesh)

!!$  ! Write original conditions
!!$  call fem_conditions_compose_name ( comp_prefix, name ) 
!!$  lunio = io_open(name)
!!$  if(gmesh%nboun>0) then
!!$     call fem_conditions_write(lunio,gnodes,gbouns)
!!$  else
!!$     call fem_conditions_write(lunio,gnodes)
!!$  end if
!!$  call io_close(lunio)

  ! Create partition
  prt_pars%nparts                = nparts
  prt_pars%strat                 = part_kway
  prt_pars%metis_option_ufactor  = 1
  prt_pars%debug                 = 0
  prt_pars%metis_option_minconn  = 1
  prt_pars%metis_option_debug    = 2

  call fem_mesh_distribution_create (prt_pars, gmesh, distr, lmesh)

  ! Output domain partition to GiD file
  call memalloc (gmesh%nelem, ldome, __FILE__,__LINE__)
  do i=1, nparts
     do j=1, distr(i)%emap%nl
        ldome(distr(i)%emap%l2g(j)) = i
     end do
  end do

  name = trim(dir_path) // trim(prefix) // '.res'
  call postpro_open_file(1,name,lupos)
  call postpro_gp_init(lupos,1,gmesh%nnode,gmesh%ndime)
  call postpro_gp(lupos,gmesh%ndime,gmesh%nnode,ldome,'EDOMS',1,1.0)
  call postpro_close_file(lupos)
  call memfree (ldome,__FILE__,__LINE__)

  ! Write partition info
  call fem_mesh_distribution_write_files ( dir_path_out, prefix, nparts, distr )

  ! Write local meshes
  call fem_mesh_write_files ( dir_path_out, prefix, nparts, lmesh )

  ! Create local conditions
  allocate(lnodes(nparts))
  do ipart=1,nparts
     call fem_conditions_create(gnodes%ncode,gnodes%nvalu,lmesh(ipart)%npoin,lnodes(ipart))
     call map_apply_g2l(distr(ipart)%nmap,gnodes%ncode,gnodes%code,lnodes(ipart)%code) ! ,nren)
     call map_apply_g2l(distr(ipart)%nmap,gnodes%nvalu,gnodes%valu,lnodes(ipart)%valu) ! ,nren) 
  end do

!!$  ! Conditions on elements (Neumann bc's)
!!$  if(gmesh%nboun.gt.0) then
!!$     allocate(lbouns(nparts))
!!$     do ipart=1,nparts
!!$        call fem_conditions_create(gbouns%ncode,gbouns%nvalu,lmesh(ipart)%nboun,lbouns(ipart))
!!$        call map_apply_g2l(distr(ipart)%bmap,gbouns%ncode,gbouns%code,lbouns(ipart)%code)
!!$        call map_apply_g2l(distr(ipart)%bmap,gbouns%nvalu,gbouns%valu,lbouns(ipart)%valu)
!!$     end do
!!$     ! Write local conditions
!!$     call fem_conditions_write_files(dir_path_out, prefix, nparts, lnodes, lbouns )
!!$  else
  ! Write local conditions
  call fem_conditions_write_files(dir_path_out, prefix, nparts, lnodes )
!!$  end if

  ! ! Create material  
  ! allocate(lmater(nparts))
  ! do ipart=1,nparts
  !    call fem_materials_create(lmesh(ipart)%nelem,lmater(ipart))
  !    call map_apply_g2l(parts(ipart)%emap,gmat%list,lmater(ipart)%list,eren)
  ! end do

  ! ! Write material
  ! call fem_materials_write_files(dir_path_out, prefix, nparts, lmater )

  ! Deallocate partition objects
  do ipart=1,nparts
     call fem_mesh_distribution_free (distr(ipart))
     call fem_mesh_free (lmesh(ipart))
     call fem_conditions_free (lnodes(ipart))
!!$     call fem_materials_free (lmater(ipart))
  end do
  deallocate (distr)
  deallocate (lmesh)
  deallocate (lnodes)
!!$  deallocate (lmater)

  ! Deallocate mesh_renum objects
  !call renum_free (nren)
  !call renum_free (eren)

  ! Deallocate mesh generated by geom.f90
!!$  call fem_materials_free(gmat)
  call fem_conditions_free(gnodes)
  call fem_mesh_free(gmesh)

  ! call mem_report

contains

  ! *****************************************************************************!
  ! Read mesh partition params from command-line options.                        ! 
  ! Command-line options processing for f90 is discussed, e.g.,                  !
  ! on the following URL: http://people.sc.fsu.edu/~jburkardt/f_src/args/args.f90!
  ! Still to confirm whether this support is standard in f90 or depends          !
  ! on the compiler (i.e., INTEL, GNU, etc.)                                     !
  ! *****************************************************************************!
  subroutine read_mesh_part_pars_cl(nparts,dir_path,prefix,dir_path_out)
    implicit none 
    character*(*), intent(out)  :: dir_path, prefix, dir_path_out
    character(len=256)          :: program_name
    character(len=256)          :: argument 
    integer                     :: numargs, iargc, nparts

    numargs = iargc()
    call getarg(0, program_name)
    if (.not. (numargs.eq.4)) then
       write(*,*) 'Usage: ', trim(program_name), ' nparts dir_path_data prefix dir_path_out'
       stop
    end if

    ! Error-check is required here !
    call getarg(1, argument)
    read (argument,*) nparts

    call getarg(2, argument)
    dir_path = trim(argument)

    call getarg(3, argument)
    prefix = trim(argument)

    call getarg(4, argument)
    dir_path_out = trim(argument)

  end subroutine read_mesh_part_pars_cl

end program partitioner
