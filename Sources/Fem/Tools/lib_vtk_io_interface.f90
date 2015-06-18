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

module lib_vtk_io_interface

  use fem, only: fem_mesh, fem_triangulation, fem_mesh_print, triangulation_print
  use memor
  use stdio
  use fem_space_types
  use fem_space_names
  use problem_names
  use types
  use Lib_VTK_IO
  use ISO_C_BINDING

  implicit none

  interface
    function mkdir_recursive(path) bind(c,name="mkdir_recursive")
      use iso_c_binding
      logical :: mkdir_recursive
      character(kind=c_char,len=1), intent(IN) :: path(*)
    end function mkdir_recursive
  end interface

  ! Type for storing field descriptors
  type vtk_field
     character(len=:), allocatable :: var_location   ! 'Node' or 'Cell' field
     character(len=:), allocatable :: var_name       ! Name of the field
     character(len=:), allocatable :: field_type     ! Field data type 'Float32', 'Float64', 'Int32', etc.
     integer(ip)                   :: num_comp = 0   ! Number of components
     logical                       :: filled = .False.
  end type vtk_field

  ! Type for storing mesh data
  type vtk_mesh
     character(len=:), allocatable :: dir_path       ! Directory where the results are going to be stored
     character(len=:), allocatable :: prefix         ! Name prefix of the VTK files
     real(rp)   , allocatable      :: X(:),Y(:),Z(:) ! Coordinates of the mesh
     integer(ip), allocatable      :: connec(:)      ! Connectivity matrix
     integer(ip), allocatable      :: offset(:)      ! VTK element offset
     integer(1) , allocatable      :: ctype(:)       ! VTK element type
     integer(ip)                   :: nnods          ! Number of nodes
     integer(ip)                   :: ndim           ! Dimensions of the mesh
     type(vtk_field), allocatable  :: fields(:)      ! Array storing fields info
     logical                       :: filled = .False.
!     integer(ip), allocatable      :: elem2subelem_i(:),elem2subelem_j(:)
!     integer(ip)                   :: num_sub_elems
!     type(array_ip2)               :: nodes_subelem(max_order)
  end type vtk_mesh

  ! Type for storing several mesh data with its field descriptors
  ! It also contains information about the number of parts (PVTK) and time steps (PVD)
  ! It stores the directory path and the prefix where to write in disk
  type fem_vtk
     type(vtk_mesh), allocatable   :: mesh(:)         ! VTK mesh data and field descriptors
     type(fem_space), pointer      :: p_f_space => NULL()  ! Poins to fem_space
     integer(ip)                   :: num_meshes = 0  ! Number of VTK meshes stored
     integer(ip)                   :: num_steps = 1   ! Number of time steps
     integer(ip)                   :: num_parts = 1   ! Number of parts
     contains
        procedure          :: initialize              ! Initialize vtk_mesh derived type. Convert mesh to VTK
        procedure          :: free                    ! Deallocate
        procedure          :: write_VTK               ! Write a VTU file
        procedure          :: write_PVTK              ! Write a PVTU file
        procedure          :: write_PVD               ! Write a PVD file
        procedure, private :: fill_mesh_from_triangulation
        procedure, private :: fill_fields_from_physical_problem
        procedure, private :: get_VTK_time_output_path
        procedure, private :: get_PVD_time_output_path
        procedure, private :: get_vtk_filename
        procedure, private :: get_pvtk_filename
        procedure, private :: set_num_steps
        procedure, private :: set_num_parts
        procedure, private :: set_dir_path
        procedure, private :: set_prefix
  end type fem_vtk

  character(len=5) :: mesh_prefix = 'mesh_'
  character(len=5) :: part_prefix = 'part_'
  character(len=5) :: time_prefix = 'time_'
  character(len=4) :: vtk_ext = '.vtu'
  character(len=4) :: pvd_ext = '.pvd'
  character(len=5) :: pvtk_ext = '.pvtu'

  integer(ip) :: ftype_conn(max_FE_types,max_nobje) ! Face connectivities for P and Q P1 elements
  integer(1)  :: celltypes(max_ndime,max_FE_types)  ! VTK cell type: (dimensions,P/Q_type_id) 
!  integer(ip) :: type_id(12)                        ! Reverse VTK cell type

  public :: fem_vtk


contains

  ! Subroutine to initialize fem_vtk derived type
  subroutine initialize(f_vtk, f_trian, f_space, phys_prob, dir_path, prefix, nparts, nsteps, nmesh)
  ! ----------------------------------------------------------------------------------
    class(fem_vtk),          intent(INOUT) :: f_vtk
    type(fem_triangulation), intent(IN)    :: f_trian
    type(fem_space), target, intent(IN)    :: f_space
    class(physical_problem), intent(IN)    :: phys_prob
    character(len=*),        intent(IN)    :: dir_path
    character(len=*),        intent(IN)    :: prefix  
    integer(ip), optional,   intent(IN)    :: nparts
    integer(ip), optional,   intent(IN)    :: nsteps
    integer(ip), optional,   intent(OUT)   :: nmesh
    integer(ip)                            :: nm
  ! ----------------------------------------------------------------------------------

    ! Face connectivities for P and Q P1 elements
    ftype_conn(:,:) = 0_ip
    ftype_conn(P_type_id,1:4) = (/1,2,3,4/)
    ftype_conn(Q_type_id,1:8) = (/1,2,4,3,5,6,8,7/)

    ! VTK cell type: (dimensions,P/Q_type_id) 
    celltypes(1,:) = 3          ! VTK_LINE
    celltypes(2,P_type_id) = 5  ! VTK_TRIANGLE
    celltypes(2,Q_type_id) = 9  ! VTK_QUAD
    celltypes(3,P_type_id) = 10 ! VTK_TETRA
    celltypes(3,Q_type_id) = 12 ! VTK_HEXAHEDRON

!    ! Reverse VTK type
!    type_id(:) = 0
!    type_id(5) =  P_type_id; type_id(10) =  P_type_id
!    type_id(9) =  Q_type_id; type_id(12) =  Q_type_id


    call f_vtk%fill_mesh_from_triangulation(f_trian, nm)
    call f_vtk%fill_fields_from_physical_problem(phys_prob, nm)
    f_vtk%p_f_space => f_space
    call f_vtk%set_dir_path(dir_path,nm)
    call f_vtk%set_prefix(prefix,nm)
    if(present(nparts)) call f_vtk%set_num_parts(nparts)
    if(present(nsteps)) call f_vtk%set_num_steps(nsteps)
    if(present(nmesh)) nmesh = nm

  ! ----------------------------------------------------------------------------------
  end subroutine initialize


  ! Subroutine to store several meshes in a fem_vtk derived type
  subroutine fill_mesh_from_triangulation(f_vtk, f_trian, nmesh)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(fem_vtk),          intent(INOUT) :: f_vtk
    type(fem_triangulation), intent(IN)    :: f_trian
    integer(ip), optional,   intent(OUT)   :: nmesh
    type(vtk_mesh), allocatable            :: f_vtk_tmp(:)
    integer(ip) :: i, j, tnnod, counter
  ! ----------------------------------------------------------------------------------

    ! Meshes allocation
    if(f_vtk%num_meshes == 0) then 
        f_vtk%num_meshes = 1
        if(allocated(f_vtk%mesh)) deallocate(f_vtk%mesh)
        allocate(f_vtk%mesh(f_vtk%num_meshes))
    else
        f_vtk%num_meshes = f_vtk%num_meshes + 1 
        call move_alloc(from=f_vtk%mesh, to=f_vtk_tmp)
        allocate(f_vtk%mesh(f_vtk%num_meshes))
        f_vtk%mesh(1:size(f_vtk_tmp,dim=1)) = f_vtk_tmp(:)
        deallocate(f_vtk_tmp)
    endif
    if(present(nmesh)) nmesh = f_vtk%num_meshes

    call memalloc ( f_trian%num_elems, f_vtk%mesh(f_vtk%num_meshes)%ctype, __FILE__,__LINE__)
    call memalloc ( f_trian%num_elems, f_vtk%mesh(f_vtk%num_meshes)%offset, __FILE__,__LINE__)
    f_vtk%mesh(f_vtk%num_meshes)%ctype = 0; f_vtk%mesh(f_vtk%num_meshes)%offset = 0

    ! Fill VTK cell type and offset arrays and and count nodes
    tnnod = 0
    do i=1, f_trian%num_elems
        f_vtk%mesh(f_vtk%num_meshes)%ctype(i) = celltypes(f_trian%num_dims,f_trian%elems(i)%topology%ftype)
        tnnod = tnnod + f_trian%elems(i)%topology%nnode
        f_vtk%mesh(f_vtk%num_meshes)%offset(i:) = tnnod
    enddo
    f_vtk%mesh(f_vtk%num_meshes)%nnods = tnnod
    f_vtk%mesh(f_vtk%num_meshes)%ndim = f_trian%num_dims

    call memalloc ( tnnod, f_vtk%mesh(f_vtk%num_meshes)%connec, __FILE__,__LINE__)
    call memalloc ( tnnod, f_vtk%mesh(f_vtk%num_meshes)%X, __FILE__,__LINE__)
    call memalloc ( tnnod, f_vtk%mesh(f_vtk%num_meshes)%Y, __FILE__,__LINE__)
    call memalloc ( tnnod, f_vtk%mesh(f_vtk%num_meshes)%Z, __FILE__,__LINE__)
    f_vtk%mesh(f_vtk%num_meshes)%X = 0._rp; f_vtk%mesh(f_vtk%num_meshes)%Y = 0._rp; f_vtk%mesh(f_vtk%num_meshes)%Z = 0._rp
    counter = 1
    tnnod = 0

    ! Fill VTK coordinate arrays
    do i=1, f_trian%num_elems
        do j=1, f_trian%elems(i)%topology%nnode
            f_vtk%mesh(f_vtk%num_meshes)%connec(tnnod+j) = ftype_conn( f_trian%elems(i)%topology%ftype,mod(j-1,f_trian%elems(i)%topology%nnode)+1) - 1 + tnnod 
            if (f_trian%num_dims >=1) f_vtk%mesh(f_vtk%num_meshes)%X(counter) = f_trian%elems(i)%coordinates(1,j)
            if (f_trian%num_dims >=2) f_vtk%mesh(f_vtk%num_meshes)%Y(counter) = f_trian%elems(i)%coordinates(2,j)
            if (f_trian%num_dims >=3) f_vtk%mesh(f_vtk%num_meshes)%Z(counter) = f_trian%elems(i)%coordinates(3,j)
            counter = counter + 1
        enddo
        tnnod = tnnod + f_trian%elems(i)%topology%nnode
    enddo
    f_vtk%mesh(f_vtk%num_meshes)%filled = .True.
  ! ----------------------------------------------------------------------------------
  end subroutine fill_mesh_from_triangulation


  ! Fill fields information into fem_vtk derivet type
  subroutine fill_fields_from_physical_problem(f_vtk, phys_prob, nmesh)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(fem_vtk),          intent(INOUT) :: f_vtk
    class(physical_problem), intent(IN)    :: phys_prob
    integer(ip), optional,   intent(IN)    :: nmesh
    type(vtk_field), allocatable           :: vtk_f_tmp(:)
    integer(ip)                            :: nm
    integer(ip)                            :: i
  ! ----------------------------------------------------------------------------------

    nm = f_vtk%num_meshes
    if(present(nmesh)) nm = nmesh

    if(.not. allocated (f_vtk%mesh(nm)%fields)) then
        allocate(f_vtk%mesh(nm)%fields(phys_prob%nunks))
    else
        deallocate(f_vtk%mesh(nm)%fields)
        allocate(f_vtk%mesh(nm)%fields(phys_prob%nunks))
    endif

    do i=1, phys_prob%nunks
        f_vtk%mesh(nm)%fields(i)%var_location = 'Node'
        f_vtk%mesh(nm)%fields(i)%field_type = 'Float64'
        f_vtk%mesh(nm)%fields(i)%var_name = 'Unknown_'//trim(adjustl(ch(i)))
        if(allocated(phys_prob%vars_of_unk)) then
            if(size(phys_prob%vars_of_unk, dim=1) >= i) f_vtk%mesh(nm)%fields(i)%num_comp = phys_prob%vars_of_unk(i)
        endif
        if(allocated(phys_prob%unkno_names)) then
            if(size(phys_prob%unkno_names, dim=1) >= i) f_vtk%mesh(nm)%fields(i)%var_name = trim(adjustl(phys_prob%unkno_names(i)))
        endif
        f_vtk%mesh(nm)%fields(i)%filled = .True.
    enddo
  ! ----------------------------------------------------------------------------------
  end subroutine fill_fields_from_physical_problem


  ! Write a single VTK file to disk
  function write_VTK(f_vtk, f_name, n_part, t_step, n_mesh, o_fmt) result(E_IO)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(fem_vtk),             intent(INOUT) :: f_vtk
    character(len=*), optional, intent(IN)    :: f_name
    integer(ip),      optional, intent(IN)    :: n_part
    integer(ip),      optional, intent(IN)    :: t_step
    integer(ip),      optional, intent(IN)    :: n_mesh
    character(len=*), optional, intent(IN)    :: o_fmt
    character(len=:), allocatable             :: fn,dp
    integer(ip)                               :: nm, np, ts
    character(len=:), allocatable             :: of
    integer(ip)                               :: fid, nnods, nels, E_IO
    logical                                   :: ok
    integer(ip)                               :: i, j, f, idx, tnnod, curr_nvar, tncomp, nnode, elnnod
    real(rp), allocatable                     :: field(:,:)
  ! ----------------------------------------------------------------------------------

    nm = f_vtk%num_meshes
    if(present(n_mesh)) nm = n_mesh

    np = f_vtk%num_parts
    if(present(n_part)) np = n_part

    ts = f_vtk%num_steps
    if(present(t_step)) ts = t_step 

    dp = f_vtk%get_VTK_time_output_path(f_path=f_vtk%mesh(nm)%dir_path, t_step=ts, n_mesh=nm)
    fn = f_vtk%get_VTK_filename(f_prefix=f_vtk%mesh(nm)%prefix, n_part=np, n_mesh=nm)
    fn = dp//fn
    if(present(f_name)) fn = f_name

    ok = mkdir_recursive(dp//C_NULL_CHAR)
!    print*, '-->', stat(dp//C_NULL_CHAR, int(o'775',c_int16_t))

    of = 'raw'
    if(present(o_fmt)) of = trim(adjustl(o_fmt))

    nnods = f_vtk%mesh(nm)%nnods
    nels = size(f_vtk%mesh(nm)%ctype, dim=1)

    E_IO = VTK_INI_XML(output_format = trim(adjustl(of)), filename = trim(adjustl(fn)), mesh_topology = 'UnstructuredGrid', cf=fid)
    E_IO = VTK_GEO_XML(NN = nnods, NC = nels,  X=f_vtk%mesh(nm)%X, Y=f_vtk%mesh(nm)%Y, Z=f_vtk%mesh(nm)%Z, cf=fid)
    E_IO = VTK_CON_XML(NC= nels, connect   = f_vtk%mesh(nm)%connec, offset    = f_vtk%mesh(nm)%offset, cell_type = f_vtk%mesh(nm)%ctype, cf=fid)
 

    if(allocated(f_vtk%mesh(nm)%fields)) then

        E_IO = VTK_DAT_XML(var_location='node',var_block_action='open', cf=fid)

        tncomp = 0
        do f=1, size(f_vtk%mesh(nm)%fields, dim=1) 
            if(f_vtk%mesh(nm)%fields(f)%filled) then
                tnnod = 0
                curr_nvar = tncomp + f_vtk%mesh(nm)%fields(f)%num_comp
!                allocate(field(f_vtk%mesh(nm)%fields(f)%num_comp,nnods))
                call memalloc( f_vtk%mesh(nm)%fields(f)%num_comp, nnods, field, __FILE__,__LINE__)
    
                do i=1, nels
                    elnnod = f_vtk%p_f_space%lelem(i)%f_inf(1)%p%nobje_dim(2)-1 !Num nodes (dim=2 -> vertex)
                    do j=1, elnnod
                        nnode = f_vtk%p_f_space%lelem(i)%f_inf(curr_nvar)%p%ntxob%p(j)
                        idx = f_vtk%p_f_space%lelem(i)%f_inf(curr_nvar)%p%ntxob%l(nnode)
                        field(1:f_vtk%mesh(nm)%fields(f)%num_comp,j+tnnod) = &
                             f_vtk%p_f_space%lelem(i)%unkno(idx, tncomp+1:tncomp+f_vtk%mesh(nm)%fields(f)%num_comp, ts)
                    enddo
                    tnnod = tnnod + elnnod 
                enddo
        
                tncomp = curr_nvar
                E_IO = VTK_VAR_XML(NC_NN=nnods,N_COL=f_vtk%mesh(nm)%fields(f)%num_comp, varname=f_vtk%mesh(nm)%fields(f)%var_name,var=field, cf=fid)
                if(allocated(field)) call memfree(field)
            endif
        enddo

        E_IO = VTK_DAT_XML(var_location='node',var_block_action='close', cf=fid)
    endif

    E_IO = VTK_GEO_XML(cf=fid)
    E_IO = VTK_END_XML()


  ! ----------------------------------------------------------------------------------
  end function write_VTK


  ! Write the PVTK file referencing serveral parts of the partitioned mesh
  function write_PVTK(f_vtk, f_name, n_mesh, t_step) result(E_IO)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(fem_vtk),             intent(INOUT) :: f_vtk
    character(len=*), optional, intent(IN)    :: f_name
    integer(ip),      optional, intent(IN)    :: n_mesh
    integer(ip),      optional, intent(IN)    :: t_step
    integer(ip)                               :: nm
    character(len=:),allocatable              :: var_name
    character(len=:),allocatable              :: fn ,dp
    integer(ip)                               :: i, fid, nnods, nels, ts, E_IO
  ! ----------------------------------------------------------------------------------

    nm = f_vtk%num_meshes
    if(present(n_mesh)) nm = n_mesh

    ts = f_vtk%num_steps
    if(present(t_step)) ts = t_step 


    dp = f_vtk%get_VTK_time_output_path(f_path=f_vtk%mesh(nm)%dir_path, t_step=ts, n_mesh=nm)
    fn = f_vtk%get_PVTK_filename(f_prefix=f_vtk%mesh(nm)%prefix, n_mesh=nm)
    fn = dp//fn
    if(present(f_name)) fn = f_name

    nnods = f_vtk%mesh(nm)%nnods
    nels = size(f_vtk%mesh(nm)%ctype, dim=1)

    ! pvtu
    E_IO = PVTK_INI_XML(filename = trim(adjustl(fn)), mesh_topology = 'PUnstructuredGrid', tp='Float64')
    do i=1, f_vtk%num_parts
        E_IO = PVTK_GEO_XML(source=trim(adjustl(f_vtk%get_VTK_filename(f_prefix=f_vtk%mesh(nm)%prefix, n_part=f_vtk%num_parts, n_mesh=nm))))
    enddo
    do i=1, size(f_vtk%mesh(nm)%fields, dim=1)
        if(f_vtk%mesh(nm)%fields(i)%filled) then
            E_IO = PVTK_DAT_XML(var_location = trim(adjustl(f_vtk%mesh(nm)%fields(i)%var_location)), var_block_action = 'OPEN')
            if(allocated(f_vtk%mesh(nm)%fields(i)%var_name)) then
                var_name = f_vtk%mesh(nm)%fields(i)%var_name
            else
                var_name = 'Unknown_'//trim(adjustl(ch(i)))
            endif
            E_IO = PVTK_VAR_XML(varname = trim(adjustl(var_name)), tp=trim(adjustl(f_vtk%mesh(nm)%fields(i)%field_type)))
        endif
        E_IO = PVTK_DAT_XML(var_location = trim(adjustl(f_vtk%mesh(nm)%fields(i)%var_location)), var_block_action = 'CLOSE')
    enddo
    E_IO = PVTK_END_XML()

  ! ----------------------------------------------------------------------------------
  end function write_PVTK


  ! Write the PVD file referencing several PVTK files in a timeline
  function write_PVD(f_vtk, f_name, n_mesh) result(E_IO)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(fem_vtk),             intent(INOUT) :: f_vtk
    character(len=*), optional, intent(IN)    :: f_name
    integer(ip),      optional, intent(IN)    :: n_mesh
    integer(ip)                               :: nm, rf
    character(len=:),allocatable              :: var_name
    character(len=:),allocatable              :: pvdfn, pvtkfn ,dp
    integer(ip)                               :: i, fid, nnods, nels, ts, E_IO
  ! ----------------------------------------------------------------------------------

    nm = f_vtk%num_meshes
    if(present(n_mesh)) nm = n_mesh

    pvdfn = trim(adjustl(f_vtk%mesh(nm)%dir_path))//'/'//trim(adjustl(f_vtk%mesh(nm)%prefix))//pvd_ext
    if(present(f_name)) pvdfn = f_name


    E_IO = PVD_INI_XML(filename=trim(adjustl(pvdfn)),cf=rf)
    do ts=1, f_vtk%num_steps
      dp = f_vtk%get_PVD_time_output_path(f_path=f_vtk%mesh(nm)%dir_path, t_step=ts)
      pvtkfn = f_vtk%get_PVTK_filename(f_prefix=f_vtk%mesh(nm)%prefix, n_mesh=nm)
      pvtkfn = dp//pvtkfn
      E_IO = PVD_DAT_XML(filename=trim(adjustl(pvtkfn)),timestep=ts, cf=rf)
    enddo
    E_IO = PVD_END_XML(cf=rf)

  ! ----------------------------------------------------------------------------------
  end function write_PVD



  ! Build time output dir path for the vtk files in each timestep
  function get_VTK_time_output_path(f_vtk, f_path, t_step, n_mesh) result(dp)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(fem_vtk),             intent(INOUT) :: f_vtk
    character(len=*), optional, intent(IN)    :: f_path
    integer(ip),      optional, intent(IN)    :: t_step
    integer(ip),      optional, intent(IN)    :: n_mesh
    character(len=:), allocatable             :: dp
    character(len=:), allocatable             :: fp
    integer(ip)                               :: nm, ts
  ! ----------------------------------------------------------------------------------

    nm = f_vtk%num_meshes
    if(present(n_mesh)) nm = n_mesh

    ts = f_vtk%num_steps
    if(present(t_step)) ts = t_step 

    fp = f_vtk%mesh(nm)%dir_path
    if(present(f_path)) fp = f_path

    dp = trim(adjustl(fp))//'/'//time_prefix//trim(adjustl(ch(ts)))//'/'

  end function get_VTK_time_output_path


  ! Build time output dir path for the vtk files in each timestep
  function get_PVD_time_output_path(f_vtk, f_path, t_step) result(dp)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(fem_vtk),             intent(INOUT) :: f_vtk
    character(len=*), optional, intent(IN)    :: f_path
    integer(ip),      optional, intent(IN)    :: t_step
    character(len=:), allocatable             :: dp
    character(len=:), allocatable             :: fp
    integer(ip)                               :: ts
  ! ----------------------------------------------------------------------------------

    ts = f_vtk%num_steps
    if(present(t_step)) ts = t_step 

    dp = time_prefix//trim(adjustl(ch(ts)))//'/'

  end function get_PVD_time_output_path


  ! Build VTK filename
  function get_VTK_filename(f_vtk, f_prefix, n_part, n_mesh) result(fn)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(fem_vtk),             intent(INOUT) :: f_vtk
    character(len=*), optional, intent(IN)    :: f_prefix
    integer(ip),      optional, intent(IN)    :: n_part
    integer(ip),      optional, intent(IN)    :: n_mesh
    character(len=:), allocatable             :: fn
    character(len=:), allocatable             :: fp
    integer(ip)                               :: nm, np
  ! ----------------------------------------------------------------------------------

    nm = f_vtk%num_meshes
    if(present(n_mesh)) nm = n_mesh

    np = f_vtk%num_parts
    if(present(n_part)) np = n_part

    fp = f_vtk%mesh(nm)%prefix
    if(present(f_prefix)) fp = f_prefix

    fn = trim(adjustl(fp))//'_'//mesh_prefix//trim(adjustl(ch(nm)))//'_'//part_prefix//trim(adjustl(ch(np)))//vtk_ext

  end function get_VTK_filename


  ! Build VTK filename
  function get_PVTK_filename(f_vtk, f_prefix, n_mesh) result(fn)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(fem_vtk),             intent(INOUT) :: f_vtk
    character(len=*), optional, intent(IN)    :: f_prefix
    integer(ip),      optional, intent(IN)    :: n_mesh
    character(len=:), allocatable             :: fn
    character(len=:), allocatable             :: fp
    integer(ip)                               :: nm
  ! ----------------------------------------------------------------------------------

    nm = f_vtk%num_meshes
    if(present(n_mesh)) nm = n_mesh

    fp = f_vtk%mesh(nm)%prefix
    if(present(f_prefix)) fp = f_prefix

    fn = trim(adjustl(fp))//'_'//mesh_prefix//trim(adjustl(ch(nm)))//pvtk_ext

  end function get_PVTK_filename




  ! Set the number of time steps of the simulation to be writen in the PVD
  subroutine set_num_steps(f_vtk, t_steps)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(fem_vtk), intent(INOUT) :: f_vtk
    integer(ip),    intent(IN)    :: t_steps
  ! ----------------------------------------------------------------------------------

    f_vtk%num_steps = t_steps
  ! ----------------------------------------------------------------------------------
  end subroutine set_num_steps


  ! Set the number of parts of the partitioned mesh to be writen in the PVTK
  subroutine set_num_parts(f_vtk, n_parts)
  ! ----------------------------------------------------------------------------------
    class(fem_vtk), intent(INOUT) :: f_vtk
    integer(ip),    intent(IN)    :: n_parts
  ! ----------------------------------------------------------------------------------

    f_vtk%num_parts = n_parts
  ! ----------------------------------------------------------------------------------
  end subroutine set_num_parts


  ! Set the name of the output directory
  subroutine set_dir_path(f_vtk, dir_path, nmesh)
  ! ----------------------------------------------------------------------------------
    class(fem_vtk),        intent(INOUT) :: f_vtk
    character(len=*),      intent(IN)    :: dir_path
    integer(ip), optional, intent(IN)    :: nmesh
    integer(ip)                          :: nm = 1
  ! ----------------------------------------------------------------------------------

    if(present(nmesh)) nm = nmesh
    f_vtk%mesh(nm)%dir_path = dir_path
  ! ----------------------------------------------------------------------------------
  end subroutine set_dir_path


  ! Set the name of the output directory
  subroutine set_prefix(f_vtk, prefix, nmesh)
  ! ----------------------------------------------------------------------------------
    class(fem_vtk),        intent(INOUT) :: f_vtk
    character(len=*),      intent(IN)    :: prefix
    integer(ip), optional, intent(IN)    :: nmesh
    integer(ip)                          :: nm = 1
  ! ----------------------------------------------------------------------------------

    if(present(nmesh)) nm = nmesh
    f_vtk%mesh(nm)%prefix = prefix
  ! ----------------------------------------------------------------------------------
  end subroutine set_prefix


  ! Free the fem_vtk derived type
  subroutine free (f_vtk)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(fem_vtk), intent(inout) :: f_vtk
    integer(ip)                   :: i, j
  ! ----------------------------------------------------------------------------------


    if(allocated(f_vtk%mesh)) then
        do i=1, size(f_vtk%mesh,1)
            call memfree(f_vtk%mesh(i)%X)
            call memfree(f_vtk%mesh(i)%Y)
            call memfree(f_vtk%mesh(i)%Z)
            call memfree(f_vtk%mesh(i)%connec)
            call memfree(f_vtk%mesh(i)%offset)
            call memfree(f_vtk%mesh(i)%ctype)
            if(allocated(f_vtk%mesh(i)%fields)) then
                do j=1, size(f_vtk%mesh(i)%fields)
                    if(allocated(f_vtk%mesh(i)%fields(j)%var_location)) deallocate(f_vtk%mesh(i)%fields(j)%var_location)
                    if(allocated(f_vtk%mesh(i)%fields(j)%var_name)) deallocate(f_vtk%mesh(i)%fields(j)%var_name)
                    if(allocated(f_vtk%mesh(i)%fields(j)%field_type)) deallocate(f_vtk%mesh(i)%fields(j)%field_type)
                    f_vtk%mesh(i)%fields(j)%filled = .False.
                enddo
                deallocate(f_vtk%mesh(i)%fields)
            endif
            f_vtk%mesh(i)%filled = .False.
        enddo
        deallocate(f_vtk%mesh)
    endif
    f_vtk%num_meshes = 0
  ! ----------------------------------------------------------------------------------
  end subroutine free


end module lib_vtk_io_interface
