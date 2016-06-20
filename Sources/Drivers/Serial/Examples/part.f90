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
  use serial_names
  implicit none
  integer(ip)                     :: nparts
  type(mesh_t)                    :: gmesh
  type(materials_t)               :: gmat
  type(mesh_distribution_t) , allocatable :: distr(:)
  type(mesh_t)      , allocatable :: lmesh(:)
  type(conditions_t), allocatable :: lnodes(:),lbouns(:)
  type(materials_t),  allocatable :: lmater(:)
  type(post_file_t)                 :: lupos
  integer(ip)                       :: ipart,idof,ndofn,ndime,prob,gid_sq
  integer(ip), parameter            :: cdr=1, sto=2, mhd=3
  integer(ip), parameter            :: gid=1, square=2
  integer(ip)                       :: isper(3),nedir(3)

  character(len=256)                :: dir_path, dir_path_out
  character(len=256)                :: name, prefix
  integer(ip), allocatable          :: ldome(:)
  integer(ip)                       :: i, j

  type(mesh_distribution_params_t) :: prt_pars

  call meminit

  call read_mesh_part_pars_cl(nparts,dir_path,prefix,dir_path_out)

  !call gmesh%read(dir_path, prefix, permute_c2z=.true.)
  call gmesh%read(dir_path, prefix)

  ! To debug
  ! call gmesh%to_dual()
  ! do i=1,gmesh%npoin
  !    write(*,*) i,gmesh%lelpo(gmesh%pelpo(i):gmesh%pelpo(i+1)-1)
  ! end do
  ! write(*,*) '----------------------------------------------------'
  ! call gmesh%generate_vefs()
  ! write(*,*) gmesh%nvefs
  ! do i=1,gmesh%nelem
  !    write(*,*) gmesh%lvefs(gmesh%pvefs(i):gmesh%pvefs(i+1)-1)
  ! end do
  ! write(*,*) '----------------------------------------------------'
  ! j=0
  ! do i=1,gmesh%nvefs
  !    if(gmesh%lvef_geo(i)/=0) then
  !       j=j+1
  !       write(*,*) i,gmesh%lvef_geo(i), gmesh%lvef_set(i)
  !    end if
  ! end do
  ! write(*,*) 'found boundaries:', j

  ! Write original mesh for postprocess
  call gmesh%write_file_for_postprocess(dir_path_out, prefix)

  ! Create partition
  prt_pars%nparts                = nparts
  prt_pars%strat                 = part_kway
  prt_pars%metis_option_ufactor  = 30
  prt_pars%debug                 = 0
  prt_pars%metis_option_minconn  = 0 
  prt_pars%metis_option_contig   = 1 
  prt_pars%metis_option_debug    = 2
  !prt_pars%metis_option_ctype    = METIS_CTYPE_SHEM ! METIS_CTYPE_RM ! Random matching
  !prt_pars%metis_option_iptype   = METIS_IPTYPE_EDGE

  call gmesh%create_distribution (prt_pars, distr, lmesh)

  ! Output domain partition to GiD file
  call memalloc (gmesh%nelem, ldome, __FILE__,__LINE__)
  do i=1, nparts
     do j=1, distr(i)%num_local_cells
        ldome(distr(i)%l2g_cells(j)) = i
     end do
  end do

  name = trim(dir_path)// '/' // trim(prefix) // '.post.res'
  call postpro_open_file(1,name,lupos)
  call postpro_gp_init(lupos,1,gmesh%nnode,gmesh%ndime)
  call postpro_gp(lupos,gmesh%ndime,gmesh%nnode,ldome,'EDOMS',1,1.0)
  call postpro_close_file(lupos)
  call memfree (ldome,__FILE__,__LINE__)

  ! Write partition info
  call mesh_distribution_write_files    ( dir_path_out, prefix, nparts, distr )

  ! Write local meshes
  call mesh_write_files                 ( dir_path_out, prefix, nparts, lmesh )
  call mesh_write_files_for_postprocess ( dir_path_out, prefix, nparts, lmesh )

  ! Deallocate partition objects
  do ipart=1,nparts
     call distr(ipart)%free()
     call lmesh(ipart)%free
  end do
  deallocate (distr)
  deallocate (lmesh)

  call gmesh%free()

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
