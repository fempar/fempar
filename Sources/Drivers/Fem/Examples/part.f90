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
  use fem
  !use mesh_graph_partition
  implicit none
  integer(ip)                       :: power, type, use_graph, nparts, nstr
  type(mesh_size)                   :: msize
  type(fem_conditions)              :: poin,line,surf
  type(fem_conditions)              :: gnodes, gbouns
  type(fem_mesh)                    :: gmesh
  type(fem_materials)               :: gmat
  type(fem_mesh_distribution) , allocatable :: distr(:)
  type(fem_mesh)      , allocatable :: lmesh(:)
  type(fem_conditions), allocatable :: lnodes(:),lbouns(:)
  type(fem_materials),  allocatable :: lmater(:)
  integer(ip)                       :: lunio
  type(post_file)                   :: lupos
  integer(ip)                       :: ipart,idof,ndofn,ndime,prob,gid_sq
  integer(ip), parameter            :: cdr=1, sto=2, mhd=3
  integer(ip), parameter            :: gid=1, square=2
  integer(ip)                       :: isper(3),nedir(3)

  character(len=256)                :: dir_path,dir_path_out
  character(len=256)                :: prefix,comp_prefix
  character(len=256)                :: name,preord
  character(len=:), allocatable     :: name_mesh
  character(len=256)                :: filemesh,filecond
  integer(ip), allocatable          :: ldome(:)
  integer(ip)                       :: i, j

  type(part_params) :: prt_pars

  call meminit

  call read_mesh_part_pars_cl(power,nparts,dir_path,dir_path_out,prefix,prob,gid_sq)

  if(gid_sq==gid) then  ! Read GiD project
     ! Problem name ( without .gid/ )
     nstr = len_trim(dir_path)
     prefix = trim(dir_path(1:nstr-5))
     dir_path = trim(dir_path)//'data/'
     comp_prefix = trim(dir_path)//trim(prefix)

     ! Read mesh
     call fem_mesh_compose_name ( comp_prefix, name_mesh ) 
     lunio = io_open(name_mesh)
     call fem_mesh_read(lunio,gmesh)
     call io_close(lunio)

     ! Read conditions
     call fem_conditions_compose_name ( comp_prefix, name ) 
     lunio = io_open(name)
     call fem_conditions_read(lunio,gmesh%npoin,gnodes,gmesh%nboun,gbouns)
     call io_close(lunio)

     ! Read materials
     call fem_materials_compose_name ( comp_prefix, name ) 
     lunio = io_open(name)
     call fem_materials_read(lunio,gmat)
     gmat%list = 1 
     call io_close(lunio)

     ! Update comp_prefix for output folder
     comp_prefix = trim(dir_path_out)//trim(prefix)

  elseif(gid_sq==square) then
     comp_prefix = trim(dir_path_out)//trim(prefix)

     ! Mesh generation data
     isper=0
     nedir(1)=2**power
     nedir(2)=2**power
     call fem_mesh_alloc(2,4,nedir,isper,gmesh) ! (ndime,nnode,nedir,isper,msh,wcoor)

     msize%ntdix = 0        ! Type of discretization in x (0=uniform, 1=refined near walls)
     msize%ntdiy = 0        ! Type of discretization in y (0=uniform, 1=refined near walls)

     msize%xleng = 1.0_rp   ! Size of the domain in x
     msize%yleng = 1.0_rp   ! Size of the domain in y
     msize%zleng = 1.0_rp   ! Size of the domain in z
     msize%zx1   = 0.1_rp   ! size of the elements at x=0   (left)  
     msize%zx2   = 0.1_rp   ! size of the elements at x=a/2 (center)
     msize%zy1   = 0.1_rp   ! size of the elements at y=0   (bottom)
     msize%zy2   = 0.1_rp   ! size of the elements at y=b/2 (center)

     call fem_mesh_box(msize,gmesh)

     ndime = 2             ! Problem dimension

     ! Degrees of freedom
     if(prob==1) then       ! CDR
        ndofn=1
     elseif(prob==2) then   ! Stokes
        ndofn=ndime+1
     elseif(prob==3) then   ! MHD
        ndofn=2*(ndime+1)
     end if
     call fem_conditions_create(ndofn,ndofn,4,poin)
     call fem_conditions_create(ndofn,ndofn,4,line)
     call fem_conditions_create(ndofn,ndofn,1,surf)

     ! Points code: 1=bottom left; 2=bottom right; 3=top right; 4=top left
     ! Lines code:  1=bottom line; 2=right line;   3=top line;  4=left line 

     if(prob==1) then            ! CDR
        ! Point conditions
        poin%code(ndofn,1)=1; poin%valu(:,1)=0.0_rp
        poin%code(ndofn,2)=1; poin%valu(:,2)=0.0_rp
        poin%code(ndofn,3)=1; poin%valu(:,3)=0.0_rp
        poin%code(ndofn,4)=1; poin%valu(:,4)=0.0_rp

        ! Line conditions
        line%code(ndofn,1)=1; line%valu(:,1)=0.0_rp
        line%code(ndofn,2)=1; line%valu(:,2)=0.0_rp
        line%code(ndofn,3)=1; line%valu(:,3)=0.0_rp
        line%code(ndofn,4)=1; line%valu(:,4)=0.0_rp

     elseif(prob==2) then          ! Stokes
        ! Point conditions
        do idof=1, ndime
           poin%code(idof,1)=1
           poin%code(idof,2)=1
           poin%code(idof,3)=1
           poin%code(idof,4)=1 
        end do
        poin%code(ndofn,1)=1  ! hay que fijar la presion en un punto en flujo cerrado
        poin%code(ndofn,2)=0
        poin%code(ndofn,3)=0
        poin%code(ndofn,4)=0
        ! Point values
        poin%valu(:,1)=0.0_rp
        poin%valu(:,2)=0.0_rp
        poin%valu(:,3)=0.0_rp
        poin%valu(:,4)=0.0_rp

        ! Line conditions
        do idof=1, ndime
           line%code(idof,1)=1
           line%code(idof,2)=1
           line%code(idof,3)=1
           line%code(idof,4)=1
        end do
        line%code(ndofn,1)=0
        line%code(ndofn,2)=0
        line%code(ndofn,3)=0
        line%code(ndofn,4)=0
        ! Line values
        line%valu(:,1)=0.0_rp
        line%valu(:,2)=0.0_rp
        line%valu(:,3)=0.0_rp
        line%valu(:,4)=0.0_rp

        ! Stokes driven cavity flow
        line%valu(1,3)=1.0_rp

     elseif(prob==3) then          ! MHD
        ! Point conditions
        do idof=1, ndime        ! Velocity
           poin%code(idof,1)=1
           poin%code(idof,2)=1
           poin%code(idof,3)=1
           poin%code(idof,4)=1 
        end do
        poin%code(ndime+1,1)=1  ! Pressure (hay que fijar la presion en un punto en flujo cerrado)
        poin%code(ndime+1,2)=0
        poin%code(ndime+1,3)=0
        poin%code(ndime+1,4)=0
        do idof=1, ndime        ! Induction
           poin%code(ndime+1+idof,1)=1
           poin%code(ndime+1+idof,2)=1
           poin%code(ndime+1+idof,3)=1
           poin%code(ndime+1+idof,4)=1 
        end do
        poin%code(ndofn,1)=1    ! Magnetic pressure
        poin%code(ndofn,2)=1
        poin%code(ndofn,3)=1
        poin%code(ndofn,4)=1
        ! Point values
        poin%valu(:,1)=0.0_rp
        poin%valu(:,2)=0.0_rp
        poin%valu(:,3)=0.0_rp
        poin%valu(:,4)=0.0_rp

        ! Line conditions
        do idof=1, ndime        ! Velocity
           line%code(idof,1)=1
           line%code(idof,2)=1
           line%code(idof,3)=1
           line%code(idof,4)=1
        end do
        line%code(ndime+1,1)=0  ! Pressure
        line%code(ndime+1,2)=0
        line%code(ndime+1,3)=0
        line%code(ndime+1,4)=0
        if(ndime==2) then
           line%code(ndime+2,1)=1; line%code(ndime+3,1)=0  ! Induction (bx; by) we impose the tangential comp.
           line%code(ndime+2,2)=0; line%code(ndime+3,2)=1
           line%code(ndime+2,3)=1; line%code(ndime+3,3)=0
           line%code(ndime+2,4)=0; line%code(ndime+3,4)=1
           line%code(ndofn,1)=1 ! Magnetic pressure
           line%code(ndofn,2)=1
           line%code(ndofn,3)=1
           line%code(ndofn,4)=1
        end if
        ! Line values
        line%valu(:,1)=0.0_rp
        line%valu(:,2)=0.0_rp
        line%valu(:,3)=0.0_rp
        line%valu(:,4)=0.0_rp

     end if

     surf%code(:,1)=0
     surf%valu(:,1)=0.0_rp

     call fem_conditions_box(gmesh,poin,line,surf,gnodes) 
  end if
  
  ! Write original mesh
  call fem_mesh_compose_name ( comp_prefix, name_mesh ) 
  lunio = io_open(name_mesh)
  call fem_mesh_write(lunio,gmesh)
  call io_close(lunio)

!!$  ! Write original conditions
!!$  call fem_conditions_compose_name ( comp_prefix, name ) 
!!$  lunio = io_open(name)
!!$  if(gmesh%nboun>0) then
!!$     call fem_conditions_write(lunio,gnodes,gbouns)
!!$  else
!!$     call fem_conditions_write(lunio,gnodes)
!!$  end if
!!$  call io_close(lunio)
!!$
!!$  ! Write original mesh for postprocess
!!$  call postpro_compose_mesh_name( comp_prefix, name ) 
!!$  lunio = io_open(name)
!!$  call fem_mesh_write(lunio,gmesh)
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
  name = trim(comp_prefix) // '.res'
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
     if(gmesh%nboun.gt.0) call fem_conditions_free (lbouns(ipart))
  end do
  deallocate (distr)
  deallocate (lmesh)
  deallocate (lnodes)
!!$  deallocate (lmater)
  if(gmesh%nboun.gt.0) deallocate (lbouns)

  ! Deallocate mesh_renum objects
  !call renum_free (nren)
  !call renum_free (eren)

  ! Deallocate mesh generated by geom.f90
!!$  call fem_materials_free(gmat)
  call fem_conditions_free(gnodes)
  if(gmesh%nboun>0) call fem_conditions_free(gbouns)
  call fem_mesh_free(gmesh)

  ! call mem_report

contains

  ! subroutine mem_report 
  !  integer(imp)  :: mmax, mcur
  !  call memmax (mmax)
  !  call memcur (mcur)
  !  write (*,*) '** [Fempar Mem Report] ** Max. Dyn. Mem.:', mmax, &
  !     &        ' bytes. Cur. Dyn. Mem.:', mcur, ' bytes'  
  ! end subroutine  

  ! *****************************************************************************!
  ! Read mesh partition params from command-line options.                        ! 
  ! Command-line options processing for f90 is discussed, e.g.,                  !
  ! on the following URL: http://people.sc.fsu.edu/~jburkardt/f_src/args/args.f90!
  ! Still to confirm whether this support is standard in f90 or depends          !
  ! on the compiler (i.e., INTEL, GNU, etc.)                                     !
  ! *****************************************************************************!
  subroutine read_mesh_part_pars_cl(power,nparts,dir_path,dir_path_out,prefix,problem,gid_sq)
    implicit none 
    character*(*), intent(out)  :: dir_path, prefix, dir_path_out
    integer(ip), intent(out)    :: power,nparts,problem,gid_sq
    character(len=256)          :: program_name
    character(len=256)          :: argument 
    integer                     :: numargs, iargc

    numargs = iargc()
    call getarg(0, program_name)
    if (.not. (numargs.eq.7.or.numargs.eq.4)) then
       write(*,*) 'Usage: ', trim(program_name), ' [gid|square] nparts dir_path_data dir_path_out [ power prefix [cdr|sto|mhd] ]'
       stop
    end if

    ! Error-check is required here !
    call getarg(1,argument)
    if (trim(argument)=='gid') then
       gid_sq=gid
    else
       if (trim(argument)=='square') then
          gid_sq=square
       else
          write(*,*) 'Usage: ', trim(program_name), ' [gid|square] nparts dir_path_data dir_path_out [ power prefix [cdr|sto|mhd] ]'
          stop
       end if
    end if

    ! Error-check is required here !
    call getarg(2, argument)
    read (argument,*) nparts

    call getarg(3, argument)
    dir_path = trim(argument)

    call getarg(4, argument)
    dir_path_out = trim(argument)

    ! Default partitioner options (square domain)
    if(gid_sq==square) then
       call getarg(5, argument)
       read (argument,*) power     

       call getarg(6, argument)
       prefix = trim(argument)

       call getarg(7, argument)
       if(trim(argument)=='cdr') then
          problem=cdr
       elseif(trim(argument)=='sto') then
          problem=sto
       else
          if (trim(argument)=='mhd') then
             problem=mhd
          else
             write(*,*) 'Usage: ', trim(program_name), ' [gid|square] nparts dir_path_data dir_path_out [ power prefix [cdr|sto|mhd] ]'
             stop
          end if
       end if
    end if

  end subroutine read_mesh_part_pars_cl

end program partitioner
