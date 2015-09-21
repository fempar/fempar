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

module lib_vtk_io_interface_names

  use mesh_names
  use triangulation_names
  use types_names
  use memor_names
  use stdio_names
  use allocatable_array_names
  use interpolation_names
  use fe_space_types_names
  use serial_fe_space_names
  use problem_names
  use interpolation_tools_names
  use abstract_environment_names
  use Lib_VTK_IO
  use ISO_C_BINDING
  use postprocess_field_names
  use hash_table_names

  implicit none
# include "debug.i90"

  private
  
  interface
    function mkdir_recursive(path) bind(c,name="mkdir_recursive")
      use iso_c_binding
      integer(kind=c_int) :: mkdir_recursive
      character(kind=c_char,len=1), intent(IN) :: path(*)
    end function mkdir_recursive
  end interface

  ! Write process status parameters
  integer(ip), parameter :: unknown = 0
  integer(ip), parameter :: started = 1
  integer(ip), parameter :: pointdata_opened = 2
  integer(ip), parameter :: pointdata_closed = 3
  integer(ip), parameter :: ended  = 4
  integer(ip), parameter :: started_read = 5

  ! Type for storing field descriptors
  type vtk_field_t
     character(len=:), allocatable :: var_location   ! 'Node' or 'Cell' field
     character(len=:), allocatable :: var_name       ! Name of the field
     character(len=:), allocatable :: field_type     ! Field data type 'Float32', 'Float64', 'Int32', etc.
     integer(ip)                   :: num_comp = 0   ! Number of components
     logical                       :: filled = .False.
  end type vtk_field_t

  ! Type for storing mesh data
  type vtk_mesh_t
     character(len=:), allocatable :: dir_path       ! Directory where the results are going to be stored
     character(len=:), allocatable :: prefix         ! Name prefix of the VTK files
     real(rp)   , allocatable      :: X(:),Y(:),Z(:) ! Coordinates of the mesh
     integer(ip), allocatable      :: connec(:)      ! Connectivity matrix
     integer(ip), allocatable      :: offset(:)      ! VTK element offset
     integer(1) , allocatable      :: ctype(:)       ! VTK element type
     integer(ip)                   :: nnods          ! Number of nodes
     integer(ip)                   :: ndim           ! Dimensions of the mesh
     type(vtk_field_t), allocatable:: unknowns(:)    ! Array storing field_ts info
     type(vtk_field_t), allocatable:: postprocess_fields(:) ! Array storing postprocess field_ts info
     logical                       :: linear_order = .False.
     logical                       :: filled = .False.
!     integer(ip), allocatable     :: elem2subelem_i(:),elem2subelem_j(:)
     integer(ip)                   :: num_sub_elems
     type(allocatable_array_ip2_t)             :: nodes_subelem(max_order)
     integer(ip)                   :: status = unknown ! Status of the write process
  end type vtk_mesh_t

  ! Type for storing several mesh data with its field descriptors
  ! It also contains information about the number of parts (PVTK) and time steps (PVD)
  ! It stores the directory path and the prefix where to write in disk
  type vtk_t
     type(vtk_mesh_t), allocatable :: mesh(:)         ! VTK mesh data and field_t descriptors
     type(serial_fe_space_t), pointer     :: fe_space => NULL()  ! Poins to fe_space_t
     class(physical_problem_t), pointer     :: p_phys_prob => NULL()  ! Poins to physical_problem_t
     class(abstract_environment_t), pointer      :: env => NULL()  ! Poins to fe_space_t
     real(rp), allocatable         :: steps(:)        ! Array of parameters (time, eigenvalues,etc.)
     integer(ip)                   :: steps_counter=0 ! time steps counter
     integer(ip)                   :: num_meshes = 0  ! Number of VTK meshes stored
     integer(ip)                   :: num_steps = 0   ! Number of time steps
     integer(ip)                   :: num_parts = 0   ! Number of parts
     integer(ip)                   :: root_proc = 0   ! Root processor
     contains
        procedure          :: initialize              ! Initialize mesh_t derived type
        procedure          :: write_VTK_start         ! Start VTU file writing
        procedure          :: write_VTK_unknowns      ! Write the unknowns into the VTU file
        procedure          :: write_VTK_field         ! Write the field into the VTU file
        procedure          :: write_VTK_end           ! Ends VTU file writing
        procedure          :: write_VTK               ! Write a VTU file
        procedure          :: write_PVTK              ! Write a PVTU file
        procedure          :: write_PVD               ! Write a PVD file
        procedure          :: read_VTK_start          ! Start VTU file writing
        procedure          :: read_VTK_unknowns       ! Read the unknowns into the VTU file
        procedure          :: read_VTK_end            ! Ends VTU file writing
        procedure          :: read_VTK                ! Read a VTU file
        procedure          :: free                    ! Deallocate
        procedure, private :: initialize_linear_order
        procedure, private :: initialize_superlinear_order
        procedure, private :: fill_mesh_from_triangulation
        procedure, private :: fill_mesh_superlinear_order
        procedure, private :: fill_unknowns_from_physical_problem
        procedure, private :: add_new_postprocess_field
        procedure, private :: create_dir_hierarchy_on_root_process
        procedure, private :: get_VTK_time_output_path
        procedure, private :: get_PVD_time_output_path
        procedure, private :: get_vtk_filename
        procedure, private :: get_pvtk_filename
        procedure, private :: set_root_proc
        procedure, private :: set_num_steps
        procedure, private :: set_num_parts
        procedure, private :: set_dir_path
        procedure, private :: set_prefix
        procedure, private :: append_step
  end type vtk_t

!  character(len=5) :: mesh_prefix = 'mesh_'
!  character(len=5) :: part_prefix = 'part_'
  character(len=5) :: time_prefix = 'time_'
  character(len=4) :: vtk_ext = '.vtu'
  character(len=4) :: pvd_ext = '.pvd'
  character(len=5) :: pvtk_ext = '.pvtu'

!  integer(ip) :: ftype_conn(max_FE_types,max_nvef) ! Connectivities for P and Q P1 elements
  integer(1)  :: celltypes(max_ndime,max_FE_types)  ! VTK cell type: (dimensions,P/Q_type_id) 
!  integer(ip) :: type_id(12)                        ! Reverse VTK cell type

  public :: vtk_t

contains

  ! Subroutine to initialize vtk_t derived type
  subroutine initialize(f_vtk, f_trian, fe_space, phys_prob, env, dir_path, prefix, root_proc, nparts, nsteps, nmesh, linear_order)
    class(vtk_t),          intent(INOUT)   :: f_vtk
    type(triangulation_t), intent(IN)      :: f_trian
    type(serial_fe_space_t), target, intent(INOUT) :: fe_space
    class(physical_problem_t), target, intent(IN)  :: phys_prob
    class(abstract_environment_t), target, intent(IN)    :: env
    character(len=*),        intent(IN)    :: dir_path
    character(len=*),        intent(IN)    :: prefix  
    integer(ip), optional,   intent(IN)    :: root_proc
    integer(ip), optional,   intent(IN)    :: nparts
    integer(ip), optional,   intent(IN)    :: nsteps
    logical,     optional,   intent(IN)    :: linear_order
    integer(ip), optional,   intent(OUT)   :: nmesh
    integer(ip)                            :: nm = 1
    logical                                :: lo = .False., ft = .False.
    integer(ip)                            :: me, np, st, rp
  ! ----------------------------------------------------------------------------------

!! Node permutation
!    ! Face connectivities for P and Q P1 elements
!    ftype_conn(:,:) = 0_ip
!    ftype_conn(P_type_id,1:4) = (/1,2,3,4/)
!    ftype_conn(Q_type_id,1:8) = (/1,2,4,3,5,6,8,7/)
!
!!    ! Reverse VTK type
!!    type_id(:) = 0
!!    type_id(5) =  P_type_id; type_id(10) =  P_type_id
!!    type_id(9) =  Q_type_id; type_id(12) =  Q_type_id

    ! VTK cell type: (dimensions,P/Q_type_id) 
    celltypes(1,:) = 3          ! VTK_LINE
    celltypes(2,P_type_id) = 5  ! VTK_TRIANGLE
    celltypes(2,Q_type_id) = 8  ! VTK_PIXEL
    celltypes(3,P_type_id) = 10 ! VTK_TETRA
    celltypes(3,Q_type_id) = 11 ! VTK_VOXEL

    if(present(linear_order)) lo = linear_order

    f_vtk%fe_space => fe_space
    f_vtk%p_phys_prob => phys_prob
    f_vtk%env => env

    me = 0; np = 1; rp = 0
    if(associated(f_vtk%env)) then 
        call f_vtk%env%info(me,np) 
        ft =  f_vtk%env%am_i_fine_task() 
    endif
    if(present(root_proc)) rp = root_proc
    call f_vtk%set_root_proc(rp)

    if(ft) then
        if(lo) then 
          call f_vtk%initialize_linear_order(f_trian, fe_space, phys_prob, dir_path, prefix, nm)
        else
          call f_vtk%initialize_superlinear_order(fe_space, phys_prob, dir_path, prefix, nm)
        endif
        if(present(nmesh)) nmesh = nm

        call f_vtk%set_dir_path(dir_path,nm)
        call f_vtk%set_prefix(prefix,nm)
        if(present(nparts)) np = nparts
        call f_vtk%set_num_parts(np)
        st = 1
        if(present(nsteps)) st = nsteps
        call f_vtk%set_num_steps(st)

    endif

  ! ----------------------------------------------------------------------------------
  end subroutine initialize


  ! Subroutine to initialize vtk_t derived type
  subroutine initialize_linear_order(f_vtk, f_trian, fe_space, phys_prob, dir_path, prefix, nmesh)
  ! ----------------------------------------------------------------------------------
    class(vtk_t),          intent(INOUT)   :: f_vtk
    type(triangulation_t), intent(IN)      :: f_trian
    type(serial_fe_space_t), target, intent(IN)   :: fe_space
    class(physical_problem_t), intent(IN)  :: phys_prob
    character(len=*),        intent(IN)    :: dir_path
    character(len=*),        intent(IN)    :: prefix  
    integer(ip), optional,   intent(OUT)   :: nmesh
    integer(ip)                            :: nm
  ! ----------------------------------------------------------------------------------

    call f_vtk%fill_mesh_from_triangulation(f_trian, nm)

    if(present(nmesh)) nmesh = nm

  ! ----------------------------------------------------------------------------------
  end subroutine initialize_linear_order


  ! Subroutine to initialize vtk_t derived type with high order mesh
  subroutine initialize_superlinear_order(f_vtk, fe_space, phys_prob, dir_path, prefix, nmesh)
  ! ----------------------------------------------------------------------------------
    class(vtk_t),          intent(INOUT)    :: f_vtk
    type(serial_fe_space_t), target, intent(INout) :: fe_space
    class(physical_problem_t), intent(IN)   :: phys_prob
    character(len=*),        intent(IN)     :: dir_path
    character(len=*),        intent(IN)     :: prefix  
    integer(ip), optional,   intent(OUT)    :: nmesh
    integer(ip)                             :: nm
  ! ----------------------------------------------------------------------------------

    call f_vtk%fill_mesh_superlinear_order(fe_space, nm)
!    call f_vtk%fill_unknowns_from_physical_problem(phys_prob, nm)

    if(present(nmesh)) nmesh = nm

  ! ----------------------------------------------------------------------------------
  end subroutine initialize_superlinear_order


  ! Subroutine to store several meshes in a vtk_t derived type
  subroutine fill_mesh_from_triangulation(f_vtk, f_trian, nmesh)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(vtk_t),          intent(INOUT) :: f_vtk
    type(triangulation_t), intent(IN)    :: f_trian
    integer(ip), optional,   intent(OUT)   :: nmesh
    type(vtk_mesh_t), allocatable            :: f_vtk_tmp(:)
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
    f_vtk%mesh(f_vtk%num_meshes)%linear_order = .True.

    call memalloc ( f_trian%num_elems, f_vtk%mesh(f_vtk%num_meshes)%ctype, __FILE__,__LINE__)
    call memalloc ( f_trian%num_elems, f_vtk%mesh(f_vtk%num_meshes)%offset, __FILE__,__LINE__)
    f_vtk%mesh(f_vtk%num_meshes)%ctype = 0; f_vtk%mesh(f_vtk%num_meshes)%offset = 0

    ! Fill VTK cell type and offset arrays and and count nodes
    tnnod = 0
    do i=1, f_trian%num_elems
        f_vtk%mesh(f_vtk%num_meshes)%ctype(i) = celltypes(f_trian%num_dims,f_trian%elems(i)%geo_reference_element%ftype)
        tnnod = tnnod + f_trian%elems(i)%geo_reference_element%nnode
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
        do j=1, f_trian%elems(i)%geo_reference_element%nnode
! node permutation (only VTK elements 9 and 12)
!            f_vtk%mesh(f_vtk%num_meshes)%connec(tnnod+j) = ftype_conn( f_trian%elems(i)%geo_reference_element%ftype,mod(j-1,f_trian%elems(i)%geo_reference_element%nnode)+1) - 1 + tnnod
            f_vtk%mesh(f_vtk%num_meshes)%connec(tnnod+j) = j + tnnod - 1
            if (f_trian%num_dims >=1) f_vtk%mesh(f_vtk%num_meshes)%X(counter) = f_trian%elems(i)%coordinates(1,j)
            if (f_trian%num_dims >=2) f_vtk%mesh(f_vtk%num_meshes)%Y(counter) = f_trian%elems(i)%coordinates(2,j)
            if (f_trian%num_dims >=3) f_vtk%mesh(f_vtk%num_meshes)%Z(counter) = f_trian%elems(i)%coordinates(3,j)
            counter = counter + 1
        enddo
        tnnod = tnnod + f_trian%elems(i)%geo_reference_element%nnode
    enddo
    f_vtk%mesh(f_vtk%num_meshes)%filled = .True.
  ! ----------------------------------------------------------------------------------
  end subroutine fill_mesh_from_triangulation

  subroutine fill_mesh_superlinear_order(f_vtk,fe_space, nmesh)
    implicit none
    ! Parmeter
    class(vtk_t)    , intent(inout) :: f_vtk
    type(serial_fe_space_t)  , intent(inout)    :: fe_space
    integer(ip),optional      , intent(out)   :: nmesh
  
    ! Local variables
    type(vtk_mesh_t), allocatable            :: f_vtk_tmp(:)
    integer(ip)           :: ielem, subelem, order, ndime, nnode, geo_nnode, f_type, pos_voint, istat
    integer(ip)           :: count_poinsX, count_subelem, lnode, gnode, num_subelems, v_key, ltype(2)
    integer(ip)           :: first_coord(fe_space%g_trian%num_dims), g_coord(fe_space%g_trian%num_dims)
    integer(ip)           :: l_coord(fe_space%g_trian%num_dims)
    type(interpolation_t)   :: interp(max_order)
    type(allocatable_array_rp2_t)       :: coords(max_order)
    
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
    f_vtk%mesh(f_vtk%num_meshes)%linear_order = .False.

    ndime = fe_space%g_trian%num_dims
    f_type = Q_type_id

    ! Construct of the interpolation and the nodes mapping for each order
    do order = 1, max_order
       ! Number of interpolation and Geometrical nodes (assumed to be linear)
       nnode = Q_nnods(ndime,order)
       geo_nnode = Q_nnods(ndime,1)

       ! Construct the matrix of the coordinates of the node in the reference element
       call allocatable_array_create(ndime,nnode,coords(order))

       ! Construct the mapping of the nodes of the subelems
       num_subelems = Q_nnods(ndime,order-1)
       call allocatable_array_create(geo_nnode,num_subelems,f_vtk%mesh(f_vtk%num_meshes)%nodes_subelem(order))
       do subelem = 1, num_subelems
          call Q_ijkg(first_coord,subelem,ndime,order-1)
          do lnode = 1, geo_nnode
             call Q_ijkg(l_coord,lnode,ndime,1)
             g_coord = first_coord + l_coord
             gnode = Q_gijk(g_coord,ndime,order)
             f_vtk%mesh(f_vtk%num_meshes)%nodes_subelem(order)%a(lnode,subelem) = gnode
          end do
       end do
    end do


    ! Count the number of subelems and points for the postprocess
    count_poinsX = 0
    count_subelem = 0
    do ielem = 1, fe_space%g_trian%num_elems
       order = maxval(fe_space%finite_elements(ielem)%order)
       num_subelems = Q_nnods(ndime,order-1)
       geo_nnode = Q_nnods(ndime,1)
       count_subelem = count_subelem + num_subelems
       count_poinsX = count_poinsX + num_subelems*geo_nnode
    end do

    ! Store the info
    f_vtk%mesh(f_vtk%num_meshes)%nnods = count_poinsX
    f_vtk%mesh(f_vtk%num_meshes)%num_sub_elems = count_subelem


    ! Allocate the vectors
    call memalloc ( f_vtk%mesh(f_vtk%num_meshes)%nnods, f_vtk%mesh(f_vtk%num_meshes)%X, __FILE__,__LINE__)
    call memalloc ( f_vtk%mesh(f_vtk%num_meshes)%nnods, f_vtk%mesh(f_vtk%num_meshes)%Y, __FILE__,__LINE__)
    call memalloc ( f_vtk%mesh(f_vtk%num_meshes)%nnods, f_vtk%mesh(f_vtk%num_meshes)%Z, __FILE__,__LINE__)    
    call memalloc ( f_vtk%mesh(f_vtk%num_meshes)%nnods, f_vtk%mesh(f_vtk%num_meshes)%connec, __FILE__,__LINE__)    
    call memalloc ( f_vtk%mesh(f_vtk%num_meshes)%num_sub_elems, f_vtk%mesh(f_vtk%num_meshes)%ctype, __FILE__,__LINE__)    
    call memalloc ( f_vtk%mesh(f_vtk%num_meshes)%num_sub_elems, f_vtk%mesh(f_vtk%num_meshes)%offset, __FILE__,__LINE__)    


    count_poinsX = 0
    count_subelem = 0
    f_vtk%mesh(f_vtk%num_meshes)%X = 0._rp; f_vtk%mesh(f_vtk%num_meshes)%Y = 0._rp; f_vtk%mesh(f_vtk%num_meshes)%Z = 0._rp
    do ielem = 1, fe_space%g_trian%num_elems
       order = maxval(fe_space%finite_elements(ielem)%order)
       num_subelems = Q_nnods(ndime,order-1)
       geo_nnode    = Q_nnods(ndime,1)
       nnode        = Q_nnods(ndime,order)

       ! Take the coordinates from the geometry mesh
       do lnode = 1, geo_nnode
!          gnode = fe_space%g_mesh%lnods(fe_space%g_mesh%pnods(ielem)+lnode-1)
!          coords(1)%a(:,lnode) = fe_space%g_mesh%coord(:,gnode)
            coords(1)%a(:,lnode) = fe_space%g_trian%elems(ielem)%coordinates(:,lnode)
       end do

       ! Interpolate to the coordinate of all the nodes
       if (order>1) then
          ltype(2) = ndime + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)
          ltype(1) = ndime + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)*order
          v_key    = (max_ndime+1)*(max_FE_types+1)*(max_order) * ltype(1) + ltype(2)
          call fe_space%pos_volume_integrator%get(key=v_key, val=pos_voint, stat = istat)
          assert(istat.ne.new_index)
          call interpolate(ndime,geo_nnode,nnode,fe_space%linter(pos_voint),coords(1)%a,coords(order)%a)
       end if
       ! Store the coordinates of each subelem
       do subelem = 1, num_subelems
         
          count_subelem = count_subelem + 1
          do lnode = 1, geo_nnode
             gnode = f_vtk%mesh(f_vtk%num_meshes)%nodes_subelem(order)%a(lnode,subelem)
             count_poinsX = count_poinsX + 1
             if (ndime >= 1) f_vtk%mesh(f_vtk%num_meshes)%X(count_poinsX) = coords(order)%a(1,gnode)
             if (ndime >= 2) f_vtk%mesh(f_vtk%num_meshes)%Y(count_poinsX) = coords(order)%a(2,gnode)
             if (ndime >= 3) f_vtk%mesh(f_vtk%num_meshes)%Z(count_poinsX) = coords(order)%a(3,gnode)
             f_vtk%mesh(f_vtk%num_meshes)%connec(count_poinsX) = count_poinsX-1
          end do

          ! Store the type of element
          assert(fe_space%finite_elements(ielem)%p_geo_reference_element%ftype == Q_type_id)
          f_vtk%mesh(f_vtk%num_meshes)%ctype(count_subelem) = celltypes(ndime,fe_space%finite_elements(ielem)%p_geo_reference_element%ftype)          

           ! Fill offset
           f_vtk%mesh(f_vtk%num_meshes)%offset(count_subelem) = count_poinsX
        end do

     end do
     
     ! Free memory
     do order = 1, max_order
        call allocatable_array_free(coords(order))
     end do

  end subroutine fill_mesh_superlinear_order

  ! Fill field_ts information into vtk_t derivet type
  subroutine fill_unknowns_from_physical_problem(f_vtk, nmesh)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(vtk_t),          intent(INOUT) :: f_vtk
    integer(ip), optional,   intent(IN)    :: nmesh
    class(physical_problem_t), pointer     :: phys_prob
    type(vtk_field_t), allocatable         :: vtk_f_tmp(:)
    integer(ip)                            :: nm
    integer(ip)                            :: i
  ! ----------------------------------------------------------------------------------

    nm = f_vtk%num_meshes
    if(present(nmesh)) nm = nmesh

    ! check that f_vtk%phys_prob has been initialized f_vtk%initilize()
    check(associated(f_vtk%p_phys_prob)) 

    phys_prob => f_vtk%p_phys_prob

    if(.not. allocated (f_vtk%mesh(nm)%unknowns)) then
        allocate(f_vtk%mesh(nm)%unknowns(phys_prob%nunks))
    else
        deallocate(f_vtk%mesh(nm)%unknowns)
        allocate(f_vtk%mesh(nm)%unknowns(phys_prob%nunks))
    endif

    do i=1, phys_prob%nunks
        f_vtk%mesh(nm)%unknowns(i)%var_location = 'node'
        f_vtk%mesh(nm)%unknowns(i)%field_type = 'Float64'
        f_vtk%mesh(nm)%unknowns(i)%var_name = 'Unknown_'//trim(adjustl(ch(i)))
        if(allocated(phys_prob%vars_of_unk)) then
            if(size(phys_prob%vars_of_unk, dim=1) >= i) f_vtk%mesh(nm)%unknowns(i)%num_comp = phys_prob%vars_of_unk(i)
        endif
        if(allocated(phys_prob%unkno_names)) then
            if(size(phys_prob%unkno_names, dim=1) >= i) f_vtk%mesh(nm)%unknowns(i)%var_name = trim(adjustl(phys_prob%unkno_names(i)))
        endif
        f_vtk%mesh(nm)%unknowns(i)%filled = .True.
    enddo
  ! ----------------------------------------------------------------------------------
  end subroutine fill_unknowns_from_physical_problem

  ! Add and fill a new postprocess_field in f_vtk%mesh(nmesh)%postprocess_field(:)
  subroutine add_new_postprocess_field(f_vtk, postprocess_field, nmesh)
  ! ----------------------------------------------------------------------------------
    class(vtk_t),              intent(INOUT) :: f_vtk
    class(postprocess_field_t), intent(in)    :: postprocess_field
    integer(ip), optional,     intent(IN)    :: nmesh
    type(vtk_field_t), allocatable           :: tmp_postprocess_fields(:) !< Temporal postprocess fields
    integer(ip)                              :: n_postprocess_fields   !< Number of temporal fiels stored
    integer(ip)                              :: nm             !< Real Number of Mesh
  ! ----------------------------------------------------------------------------------

    nm = f_vtk%num_meshes
    if(present(nmesh)) nm = nmesh
    
    if (allocated(f_vtk%mesh(nm)%postprocess_fields)) then
       n_postprocess_fields = size(f_vtk%mesh(nm)%postprocess_fields)
       allocate(tmp_postprocess_fields(n_postprocess_fields))
       tmp_postprocess_fields = f_vtk%mesh(nm)%postprocess_fields
       deallocate(f_vtk%mesh(nm)%postprocess_fields)
       allocate(f_vtk%mesh(nm)%postprocess_fields(n_postprocess_fields+1))
       f_vtk%mesh(nm)%postprocess_fields(1:n_postprocess_fields) = tmp_postprocess_fields
       deallocate(tmp_postprocess_fields)
    else
       allocate(f_vtk%mesh(nm)%postprocess_fields(1))
       n_postprocess_fields = 0
    end if
    f_vtk%mesh(nm)%postprocess_fields(n_postprocess_fields+1)%var_location = 'node'
    f_vtk%mesh(nm)%postprocess_fields(n_postprocess_fields+1)%var_name     = postprocess_field%name
    f_vtk%mesh(nm)%postprocess_fields(n_postprocess_fields+1)%field_type   = 'Float64'
    f_vtk%mesh(nm)%postprocess_fields(n_postprocess_fields+1)%num_comp     = postprocess_field%nvars
    f_vtk%mesh(nm)%postprocess_fields(n_postprocess_fields+1)%filled       = .true.
  ! ----------------------------------------------------------------------------------
  end subroutine add_new_postprocess_field

  ! Start the write of a single VTK file to disk (if I am fine MPI task)
  function write_VTK_start(f_vtk, f_name, n_part, t_step, n_mesh, o_fmt, f_id) result(E_IO)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(vtk_t),             intent(INOUT)   :: f_vtk   !< VTK_t derived type
    character(len=*), optional, intent(IN)    :: f_name  !< VTK File NAME
    integer(ip),      optional, intent(IN)    :: n_part  !< Number of the PART
    real(rp),         optional, intent(IN)    :: t_step  !< Time STEP value
    integer(ip),      optional, intent(IN)    :: n_mesh  !< Number of the MESH
    character(len=*), optional, intent(IN)    :: o_fmt   !< Ouput ForMaT
    integer(ip),      optional, intent(OUT)   :: f_id    !< File ID
    character(len=:), allocatable             :: fn      !< Real File Name
    character(len=:), allocatable             :: dp      !< Real Directory Path
    character(len=:), allocatable             :: of      !< Real Output Format
    real(rp)                                  :: ts      !< Real Time Step
    logical                                   :: ft      !< Fine Task
    integer(ip)                               :: nm      !< Real Number of the Mesh
    integer(ip)                               :: np      !< Real Number of the Part
    integer(ip)                               :: me      !< Task identifier
    integer(ip)                               :: fid     !< Real File ID
    integer(ip)                               :: nnods   !< Number of NODeS
    integer(ip)                               :: nels    !< Number of ELementS
    integer(ip)                               :: E_IO    !< Error IO
  ! ----------------------------------------------------------------------------------

    check(associated(f_vtk%env))
 
    ft =  f_vtk%env%am_i_fine_task() 

    E_IO = 0
    fid = -1
    
    if(ft) then
        me = 0; np = 1
        call f_vtk%env%info(me,np) 
        np = me
        if(present(n_part)) np = n_part

        nm = f_vtk%num_meshes
        if(present(n_mesh)) nm = n_mesh
    
        ts = 0._rp
        if(present(t_step)) ts = t_step 
        call f_vtk%append_step(ts)


        if(f_vtk%mesh(nm)%status == unknown .or. f_vtk%mesh(nm)%status == ended) then
    
            dp = f_vtk%get_VTK_time_output_path(f_path=f_vtk%mesh(nm)%dir_path, t_step=ts, n_mesh=nm)
            fn = f_vtk%get_VTK_filename(f_prefix=f_vtk%mesh(nm)%prefix, n_part=np, n_mesh=nm)
            fn = dp//fn
            if(present(f_name)) fn = f_name

            if( f_vtk%create_dir_hierarchy_on_root_process(dp,issue_final_barrier=.True.) == 0) then    
                of = 'raw'
                if(present(o_fmt)) of = trim(adjustl(o_fmt))
        
                nnods = f_vtk%mesh(nm)%nnods
                nels = size(f_vtk%mesh(nm)%ctype, dim=1)
        
                E_IO = VTK_INI_XML(output_format = trim(adjustl(of)), filename = trim(adjustl(fn)), mesh_topology = 'UnstructuredGrid', cf=fid)
                E_IO = VTK_GEO_XML(NN = nnods, NC = nels,  X=f_vtk%mesh(nm)%X, Y=f_vtk%mesh(nm)%Y, Z=f_vtk%mesh(nm)%Z, cf=fid)
                E_IO = VTK_CON_XML(NC= nels, connect   = f_vtk%mesh(nm)%connec, offset    = f_vtk%mesh(nm)%offset, cell_type = f_vtk%mesh(nm)%ctype, cf=fid)
                f_vtk%mesh(nm)%status = started
            endif
        endif
        if(present(f_id)) f_id = fid

    endif

  ! ----------------------------------------------------------------------------------
  end function write_VTK_start

  ! Write a single VTK file to disk (if I am fine MPI task)
  function write_VTK_unknowns(f_vtk, n_mesh, f_id) result(E_IO)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(vtk_t),             intent(INOUT)   :: f_vtk          !< VTK_t derived type
    integer(ip),      optional, intent(IN)    :: n_mesh         !< Number of MESH
    integer(ip),      optional, intent(IN)    :: f_id           !< File ID
    real(rp), allocatable                     :: field(:,:)     !< FIELD(ncomp,nnod)
    integer(ip)                               :: nm             !< Real Number of Mesh
    integer(ip)                               :: E_IO           !< IO Error
    logical                                   :: ft             !< Fine Task

    check(associated(f_vtk%env))
    ft =  f_vtk%env%am_i_fine_task() 

    E_IO = 0
    
    if(ft) then        

       nm = f_vtk%num_meshes
       if(present(n_mesh)) nm = n_mesh

       if(f_vtk%mesh(f_vtk%num_meshes)%linear_order .eqv. .true.) then
          E_IO = write_VTK_unknowns_linear(f_vtk, nm, f_id)
       else
          E_IO = write_VTK_unknowns_superlinear(f_vtk, nm, f_id)
       end if
          

    end if

  end function write_VTK_unknowns

  ! Write a single VTK file to disk (if I am fine MPI task) for linear interpolation
  function write_VTK_unknowns_linear(f_vtk, nm, f_id) result(E_IO)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(vtk_t),             intent(INOUT)   :: f_vtk          !< VTK_t derived type
    integer(ip),                intent(IN)    :: nm             !< Real Number of Mesh
    integer(ip),      optional, intent(IN)    :: f_id           !< File ID
    real(rp), allocatable                     :: field(:,:)     !< FIELD(ncomp,nnod)
    integer(ip)                               :: tidx           !< Actual Time InDeX
    integer(ip)                               :: nnods          !< Number of NODes
    integer(ip)                               :: nels           !< Number of ELementS
    integer(ip)                               :: tnnod          !< Total Number of NODes
    integer(ip)                               :: curr_ncomp     !< CURRent Number of COMPonent
    integer(ip)                               :: tncomp         !< Total Number of COMPonents
    integer(ip)                               :: nnode          !< NODE Number
    integer(ip)                               :: elnnod         !< Number of NODes per ELement
    integer(ip)                               :: i, j, f, idx   !< Indices
    integer(ip)                               :: E_IO           !< IO Error
    logical                                   :: ft             !< Fine Task
  ! ----------------------------------------------------------------------------------
    
    tidx = 1

    if(.not. allocated(f_vtk%mesh(nm)%unknowns)) call f_vtk%fill_unknowns_from_physical_problem(nm)

    if(allocated(f_vtk%mesh(nm)%unknowns) .and. (f_vtk%mesh(nm)%status >= started) .and. &
         (f_vtk%mesh(nm)%status < pointdata_closed)) then

       if(f_vtk%mesh(nm)%status < pointdata_opened) then
          if(present(f_id)) then
             E_IO = VTK_DAT_XML(var_location='node',var_block_action='open', cf=f_id)
          else
             E_IO = VTK_DAT_XML(var_location='node',var_block_action='open')
          endif
          f_vtk%mesh(nm)%status = pointdata_opened
       endif

       nnods = f_vtk%mesh(nm)%nnods
       nels = size(f_vtk%mesh(nm)%ctype, dim=1)

       tncomp = 0
       do f=1, size(f_vtk%mesh(nm)%unknowns, dim=1) 
          if(f_vtk%mesh(nm)%unknowns(f)%filled) then
             tnnod = 0
             curr_ncomp = tncomp + f_vtk%mesh(nm)%unknowns(f)%num_comp
             !                   allocate(field(f_vtk%mesh(nm)%unknowns(f)%num_comp,nnods))
             call memalloc( f_vtk%mesh(nm)%unknowns(f)%num_comp, nnods, field, __FILE__,__LINE__)

             do i=1, nels
                elnnod = f_vtk%fe_space%finite_elements(i)%reference_element_vars(1)%p%nvef_dim(2)-1 !Num nodes (dim=2 -> vertex)
                do j=1, elnnod
                   nnode = f_vtk%fe_space%finite_elements(i)%reference_element_vars(curr_ncomp)%p%ntxob%p(j)
                   idx = f_vtk%fe_space%finite_elements(i)%reference_element_vars(curr_ncomp)%p%ntxob%l(nnode)
                   field(1:f_vtk%mesh(nm)%unknowns(f)%num_comp,j+tnnod) = &
                        f_vtk%fe_space%finite_elements(i)%unkno(idx, tncomp+1:tncomp+f_vtk%mesh(nm)%unknowns(f)%num_comp, tidx)
                enddo
                tnnod = tnnod + elnnod 
             enddo

             tncomp = curr_ncomp
             if(present(f_id)) then 
                E_IO = VTK_VAR_XML(NC_NN=nnods,N_COL=f_vtk%mesh(nm)%unknowns(f)%num_comp, varname=f_vtk%mesh(nm)%unknowns(f)%var_name,var=field, cf=f_id)
             else
                E_IO = VTK_VAR_XML(NC_NN=nnods,N_COL=f_vtk%mesh(nm)%unknowns(f)%num_comp, varname=f_vtk%mesh(nm)%unknowns(f)%var_name,var=field)
             endif
             if(allocated(field)) call memfree(field, __FILE__,__LINE__)
          endif
       enddo

    endif

  ! ----------------------------------------------------------------------------------
  end function write_VTK_unknowns_linear

  ! Write a single VTK file to disk (if I am fine MPI task) for high order interpolation
  function write_VTK_unknowns_superlinear(f_vtk, nm, f_id) result(E_IO)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(vtk_t),             intent(INOUT)   :: f_vtk          !< VTK_t derived type
    integer(ip),                intent(IN)    :: nm             !< Real Number of Mesh
    integer(ip),      optional, intent(IN)    :: f_id           !< File ID
    real(rp), allocatable                     :: field(:,:)     !< FIELD(ncomp,nnod)
    integer(ip)                               :: tidx           !< Actual Time InDeX
    integer(ip)                               :: nnods          !< Number of NODes
    integer(ip)                               :: nels           !< Number of ELementS
    integer(ip)                               :: nsubels        !< Number of sub_ELementS
    integer(ip)                               :: tnnod          !< Total Number of NODes
    integer(ip)                               :: curr_ncomp     !< CURRent Number of COMPonent
    integer(ip)                               :: tncomp         !< Total Number of COMPonents
    integer(ip)                               :: nnode, gnode   !< NODE Number
    integer(ip)                               :: ndime          !< Physical dimensions
    integer(ip)                               :: elnnod         !< Number of NODes per ELement
    integer(ip)                               :: i, j, f, idx   !< Indices
    integer(ip)                               :: subelem, icomp !< Indices
    integer(ip)                               :: order          !< Order of interpolation
    integer(ip)                               :: E_IO           !< IO Error
    logical                                   :: ft             !< Fine Task
    integer(ip)                               :: f_type, pos_voint, istat, v_key, ltype(2)
    real(rp), allocatable                     :: origin_field(:,:), target_field(:,:)
  ! ----------------------------------------------------------------------------------
    
    tidx = 1

    if(.not. allocated(f_vtk%mesh(nm)%unknowns)) call f_vtk%fill_unknowns_from_physical_problem(nm)

    if(allocated(f_vtk%mesh(nm)%unknowns) .and. (f_vtk%mesh(nm)%status >= started) .and. &
         (f_vtk%mesh(nm)%status < pointdata_closed)) then

       if(f_vtk%mesh(nm)%status < pointdata_opened) then
          if(present(f_id)) then
             E_IO = VTK_DAT_XML(var_location='node',var_block_action='open', cf=f_id)
          else
             E_IO = VTK_DAT_XML(var_location='node',var_block_action='open')
          endif
          f_vtk%mesh(nm)%status = pointdata_opened
       endif

       nnods = f_vtk%mesh(nm)%nnods
       nels = f_vtk%fe_space%g_trian%num_elems
       ndime = f_vtk%fe_space%g_trian%num_dims
       f_type = Q_type_id

       tncomp = 0
       do f=1, size(f_vtk%mesh(nm)%unknowns, dim=1) 
          if(f_vtk%mesh(nm)%unknowns(f)%filled) then
             tnnod = 0
             curr_ncomp = tncomp + f_vtk%mesh(nm)%unknowns(f)%num_comp
             !                   allocate(field(f_vtk%mesh(nm)%unknowns(f)%num_comp,nnods))
             call memalloc( f_vtk%mesh(nm)%unknowns(f)%num_comp, nnods, field, __FILE__,__LINE__)

             do i=1, nels

                ! Check (only working for Q-type elements)
                check(f_vtk%fe_space%finite_elements(i)%p_geo_reference_element%ftype == Q_type_id)

                ! Set order of interpolation
                order = maxval(f_vtk%fe_space%finite_elements(i)%order)
                nsubels = Q_nnods(ndime,order-1)
                gnode   = Q_nnods(ndime,1)
                nnode   = Q_nnods(ndime,order)

                ! Interpolate unknowns with different order of interpolation
                if(f_vtk%fe_space%finite_elements(i)%order(curr_ncomp).ne.order) then

                   ! Allocate auxiliar target and origin fields
                   call memalloc(f_vtk%mesh(nm)%unknowns(f)%num_comp,nnode,target_field,__FILE__,__LINE__)
                   call memalloc(f_vtk%mesh(nm)%unknowns(f)%num_comp,gnode,origin_field,__FILE__,__LINE__)
                   do icomp=1,f_vtk%mesh(nm)%unknowns(f)%num_comp
                      origin_field(icomp,:) = f_vtk%fe_space%finite_elements(i)%unkno(1:gnode,tncomp+icomp,tidx)
                   end do
                   target_field = 0.0_rp

                   ! Interpolate fields
                   ltype(2) = ndime + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)
                   ltype(1) = ndime + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)*order
                   v_key    = (max_ndime+1)*(max_FE_types+1)*(max_order) * ltype(1) + ltype(2)
                   call f_vtk%fe_space%pos_volume_integrator%get(key=v_key, val=pos_voint, stat = istat)
                   assert(istat.ne.new_index)
                   call interpolate(f_vtk%mesh(nm)%unknowns(f)%num_comp,gnode,nnode,f_vtk%fe_space%linter(pos_voint), &
                        &           origin_field,target_field)

                   do icomp=1,f_vtk%mesh(nm)%unknowns(f)%num_comp
                      f_vtk%fe_space%finite_elements(i)%unkno(:,tncomp+icomp,tidx) = target_field(icomp,:)
                   end do

                   ! Deallocate target and origin fields
                   call memfree(origin_field,__FILE__,__LINE__)
                   call memfree(target_field,__FILE__,__LINE__)

                end if

                ! Loop over subelements
                do subelem = 1,nsubels

                   ! Loop over geometrical nodes in subelement
                   do j=1, gnode
                      idx = f_vtk%mesh(nm)%nodes_subelem(order)%a(j,subelem)
                      field(1:f_vtk%mesh(nm)%unknowns(f)%num_comp,j+tnnod) = &
                           f_vtk%fe_space%finite_elements(i)%unkno(idx, tncomp+1:tncomp+f_vtk%mesh(nm)%unknowns(f)%num_comp, tidx)
                   enddo
                   tnnod = tnnod + gnode
                end do
             enddo

             tncomp = curr_ncomp
             if(present(f_id)) then 
                E_IO = VTK_VAR_XML(NC_NN=nnods,N_COL=f_vtk%mesh(nm)%unknowns(f)%num_comp, varname=f_vtk%mesh(nm)%unknowns(f)%var_name,var=field, cf=f_id)
             else
                E_IO = VTK_VAR_XML(NC_NN=nnods,N_COL=f_vtk%mesh(nm)%unknowns(f)%num_comp, varname=f_vtk%mesh(nm)%unknowns(f)%var_name,var=field)
             endif
             if(allocated(field)) call memfree(field, __FILE__,__LINE__)
          endif
       enddo

    endif

  ! ----------------------------------------------------------------------------------
  end function write_VTK_unknowns_superlinear

  ! Write a single VTK file to disk (if I am fine MPI task)
  function write_VTK_field(f_vtk, postproc_field, n_mesh, f_id) result(E_IO)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(vtk_t),             intent(INOUT)   :: f_vtk          !< VTK_t derived type
    class(postprocess_field_t), intent(inout) :: postproc_field !< Postprocess field structure to be written
    integer(ip),      optional, intent(IN)    :: n_mesh         !< Number of MESH
    integer(ip),      optional, intent(IN)    :: f_id           !< File ID
    real(rp), allocatable                     :: field(:,:)     !< FIELD(ncomp,nnod)
    integer(ip)                               :: nm             !< Real Number of Mesh
    integer(ip)                               :: E_IO           !< IO Error
    logical                                   :: ft             !< Fine Task

    check(associated(f_vtk%env))
    ft =  f_vtk%env%am_i_fine_task() 

    E_IO = 0
    
    if(ft) then        

       nm = f_vtk%num_meshes
       if(present(n_mesh)) nm = n_mesh

       if(f_vtk%mesh(f_vtk%num_meshes)%linear_order .eqv. .true.) then
          E_IO = write_VTK_field_linear(f_vtk, postproc_field, nm, f_id)
       else
          E_IO = write_VTK_field_superlinear(f_vtk, postproc_field, nm, f_id)
       end if
          

    end if

  end function write_VTK_field

  ! Write a single VTK file to disk (if I am fine MPI task)
  function write_VTK_field_linear(f_vtk,postproc_field, nm, f_id) result(E_IO)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(vtk_t),             intent(INOUT)   :: f_vtk          !< VTK_t derived type
    class(postprocess_field_t), intent(inout) :: postproc_field !< Postprocess field structure to be written
    integer(ip),                intent(IN)    :: nm             !< Number of MESH
    integer(ip),      optional, intent(IN)    :: f_id           !< File ID
    real(rp), allocatable                     :: field(:,:)     !< FIELD(ncomp,nnod)
    integer(ip)                               :: tidx           !< Actual Time InDeX
    integer(ip)                               :: nnods          !< Number of NODes
    integer(ip)                               :: nels           !< Number of ELementS
    integer(ip)                               :: tnnod          !< Total Number of NODes
    integer(ip)                               :: curr_ncomp     !< CURRent Number of COMPonent
    integer(ip)                               :: tncomp         !< Total Number of COMPonents
    integer(ip)                               :: inode          !< NODE Number
    integer(ip)                               :: elnnod         !< Number of NODes per ELement
    integer(ip)                               :: i, j, f, idx   !< Indices
    integer(ip)                               :: E_IO           !< IO Error
    logical                                   :: ft             !< Fine Task
  ! ----------------------------------------------------------------------------------

    check(associated(f_vtk%env))
    ft =  f_vtk%env%am_i_fine_task() 

    E_IO = 0
    
    if(ft) then
                
        tidx = 1
            
        if (f_vtk%mesh(nm)%status >= started .and. f_vtk%mesh(nm)%status < pointdata_closed) then
        
           call f_vtk%add_new_postprocess_field(postproc_field,nm)

           if (f_vtk%mesh(nm)%status < pointdata_opened) then
              if(present(f_id)) then
                 E_IO = VTK_DAT_XML(var_location='node',var_block_action='open', cf=f_id)
              else
                 E_IO = VTK_DAT_XML(var_location='node',var_block_action='open')
              endif
              f_vtk%mesh(nm)%status = pointdata_opened
           end if
           
           tnnod = 0
           nnods = f_vtk%mesh(nm)%nnods
           nels = size(f_vtk%mesh(nm)%ctype, dim=1)
           
           check(postproc_field%is_finalized())
           
           call memalloc( postproc_field%nvars, nnods, field, __FILE__,__LINE__)
           
           do i=1, nels
              ! Search a component of element i with the same interpolation as the postprocess field
              curr_ncomp = 0
              do f=1, f_vtk%p_phys_prob%nunks
                 curr_ncomp = curr_ncomp + f_vtk%mesh(nm)%unknowns(f)%num_comp
                 if (f_vtk%fe_space%finite_elements(i)%reference_element_vars(f)%p%nnode.eq.postproc_field%fe_postprocess_field(i)%nnode) exit
              end do
              ! Set the elemental nodes in the VTK interpolation
              elnnod = f_vtk%fe_space%finite_elements(i)%reference_element_vars(1)%p%nvef_dim(2)-1 !Num nodes (dim=2 -> vertex)
              do j=1, elnnod
                 inode = f_vtk%fe_space%finite_elements(i)%reference_element_vars(curr_ncomp)%p%ntxob%p(j)
                 idx = f_vtk%fe_space%finite_elements(i)%reference_element_vars(curr_ncomp)%p%ntxob%l(inode)

                 field(1:postproc_field%nvars,j+tnnod) = &
                      postproc_field%fe_postprocess_field(i)%nodal_properties(idx, 1:postproc_field%nvars)
              enddo
              tnnod = tnnod + elnnod 
           enddo
           
           if(present(f_id)) then 
              E_IO = VTK_VAR_XML(NC_NN=nnods,N_COL=postproc_field%nvars, varname=postproc_field%name,var=field, cf=f_id)
           else
              E_IO = VTK_VAR_XML(NC_NN=nnods,N_COL=postproc_field%nvars, varname=postproc_field%name,var=field)
           endif
           if(allocated(field)) call memfree(field)

        endif
        
     endif

  ! ----------------------------------------------------------------------------------
   end function write_VTK_field_linear


  ! Write a single VTK file to disk (if I am fine MPI task) for high order interpolation
  function write_VTK_field_superlinear(f_vtk,postproc_field, nm, f_id) result(E_IO)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(vtk_t),             intent(INOUT)   :: f_vtk          !< VTK_t derived type
    class(postprocess_field_t), intent(inout) :: postproc_field !< Postprocess field structure to be written
    integer(ip),                intent(IN)    :: nm             !< Real Number of Mesh
    integer(ip),      optional, intent(IN)    :: f_id           !< File ID
    real(rp), allocatable                     :: field(:,:)     !< FIELD(ncomp,nnod)
    integer(ip)                               :: tidx           !< Actual Time InDeX
    integer(ip)                               :: nnods          !< Number of NODes
    integer(ip)                               :: nels           !< Number of ELementS
    integer(ip)                               :: nsubels        !< Number of sub_ELementS
    integer(ip)                               :: tnnod          !< Total Number of NODes
    integer(ip)                               :: curr_ncomp     !< CURRent Number of COMPonent
    integer(ip)                               :: tncomp         !< Total Number of COMPonents
    integer(ip)                               :: nnode, gnode   !< NODE Number
    integer(ip)                               :: ndime          !< Physical dimensions
    integer(ip)                               :: elnnod         !< Number of NODes per ELement
    integer(ip)                               :: i, j, f, idx   !< Indices
    integer(ip)                               :: subelem, icomp !< Indices
    integer(ip)                               :: order          !< Order of interpolation
    integer(ip)                               :: E_IO           !< IO Error
    logical                                   :: ft             !< Fine Task
    integer(ip)                               :: f_type, pos_voint, istat, v_key, ltype(2)
    real(rp), allocatable                     :: origin_field(:,:), target_field(:,:)
  ! ----------------------------------------------------------------------------------
    
    tidx = 1

    if((f_vtk%mesh(nm)%status >= started) .and. (f_vtk%mesh(nm)%status < pointdata_closed)) then

       call f_vtk%add_new_postprocess_field(postproc_field,nm)

       if(f_vtk%mesh(nm)%status < pointdata_opened) then
          if(present(f_id)) then
             E_IO = VTK_DAT_XML(var_location='node',var_block_action='open', cf=f_id)
          else
             E_IO = VTK_DAT_XML(var_location='node',var_block_action='open')
          endif
          f_vtk%mesh(nm)%status = pointdata_opened
       endif

       nnods = f_vtk%mesh(nm)%nnods
       nels = f_vtk%fe_space%g_trian%num_elems
       ndime = f_vtk%fe_space%g_trian%num_dims
       f_type = Q_type_id

       tncomp = 0

       check(postproc_field%is_finalized())

       tnnod = 0
       
       f = size(f_vtk%mesh(nm)%postprocess_fields,dim=1)

       call memalloc( f_vtk%mesh(nm)%postprocess_fields(f)%num_comp, nnods, field, __FILE__,__LINE__)
       
       do i=1, nels
          ! Search a component of element i with the same interpolation as the postprocess field
          curr_ncomp = 0
          do j=1, f_vtk%p_phys_prob%nunks
             curr_ncomp = curr_ncomp + f_vtk%mesh(nm)%unknowns(j)%num_comp
             if (f_vtk%fe_space%finite_elements(i)%reference_element_vars(j)%p%nnode.eq.postproc_field%fe_postprocess_field(i)%nnode) exit
          end do
          
          ! Check (only working for Q-type elements)
          check(f_vtk%fe_space%finite_elements(i)%p_geo_reference_element%ftype == Q_type_id)
          
          ! Set order of interpolation
          order = maxval(f_vtk%fe_space%finite_elements(i)%order)
          nsubels = Q_nnods(ndime,order-1)
          gnode   = Q_nnods(ndime,1)
          nnode   = Q_nnods(ndime,order)
          
          ! Interpolate unknowns with different order of interpolation
          if(f_vtk%fe_space%finite_elements(i)%order(curr_ncomp).ne.order) then
             
             ! Allocate auxiliar target and origin fields
             call memalloc(f_vtk%mesh(nm)%postprocess_fields(f)%num_comp,nnode,target_field,__FILE__,__LINE__)
             call memalloc(f_vtk%mesh(nm)%postprocess_fields(f)%num_comp,gnode,origin_field,__FILE__,__LINE__)
             do icomp=1,f_vtk%mesh(nm)%postprocess_fields(f)%num_comp
                origin_field(icomp,:) = postproc_field%fe_postprocess_field(i)%nodal_properties(1:gnode,tncomp+icomp)
             end do
             target_field = 0.0_rp
             
             ! Interpolate fields
             ltype(2) = ndime + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)
             ltype(1) = ndime + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)*order
             v_key    = (max_ndime+1)*(max_FE_types+1)*(max_order) * ltype(1) + ltype(2)
             call f_vtk%fe_space%pos_volume_integrator%get(key=v_key, val=pos_voint, stat = istat)
             assert(istat.ne.new_index)
             call interpolate(f_vtk%mesh(nm)%postprocess_fields(f)%num_comp,gnode,nnode,f_vtk%fe_space%linter(pos_voint), &
                  &           origin_field,target_field)
             
             do icomp=1,f_vtk%mesh(nm)%postprocess_fields(f)%num_comp
                postproc_field%fe_postprocess_field(i)%nodal_properties(:,tncomp+icomp) = target_field(icomp,:)
             end do
             
             ! Deallocate target and origin fields
             call memfree(origin_field,__FILE__,__LINE__)
             call memfree(target_field,__FILE__,__LINE__)

          end if
          
          ! Loop over subelements
          do subelem = 1,nsubels
             
             ! Loop over geometrical nodes in subelement
             do j=1, gnode
                idx = f_vtk%mesh(nm)%nodes_subelem(order)%a(j,subelem)
                field(1:f_vtk%mesh(nm)%postprocess_fields(f)%num_comp,j+tnnod) = &
                     postproc_field%fe_postprocess_field(i)%nodal_properties(idx, tncomp+1:tncomp+f_vtk%mesh(nm)%postprocess_fields(f)%num_comp)
             enddo
             tnnod = tnnod + gnode
          end do
       enddo
       
       tncomp = curr_ncomp
       if(present(f_id)) then 
          E_IO = VTK_VAR_XML(NC_NN=nnods,N_COL=f_vtk%mesh(nm)%postprocess_fields(f)%num_comp, varname=f_vtk%mesh(nm)%postprocess_fields(f)%var_name,var=field, cf=f_id)
       else
          E_IO = VTK_VAR_XML(NC_NN=nnods,N_COL=f_vtk%mesh(nm)%postprocess_fields(f)%num_comp, varname=f_vtk%mesh(nm)%postprocess_fields(f)%var_name,var=field)
       endif
       if(allocated(field)) call memfree(field, __FILE__,__LINE__)
    endif


  ! ----------------------------------------------------------------------------------
  end function write_VTK_field_superlinear

  ! Write a single VTK file to disk (if I am fine MPI task)
  function write_VTK_end(f_vtk, n_mesh, f_id) result(E_IO)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(vtk_t),               intent(INOUT) :: f_vtk   !< VTK_t derived type
    integer(ip),      optional, intent(IN)    :: n_mesh  !< Number of MESH
    integer(ip),      optional, intent(INOUT) :: f_id    !< File ID
    integer(ip)                               :: nm      !< Real Number of Mesh
    integer(ip)                               :: E_IO    !< IO Error
    logical                                   :: ft      !< Fine Task
  ! ----------------------------------------------------------------------------------

    check(associated(f_vtk%env))
    ft =  f_vtk%env%am_i_fine_task() 

    E_IO = 0

    if (ft) then
       
       nm = f_vtk%num_meshes
       if(present(n_mesh)) nm = n_mesh
       
       if ((f_vtk%mesh(nm)%status >= started) .and. (f_vtk%mesh(nm)%status /= ended)) then
          
          if(f_vtk%mesh(nm)%status == pointdata_opened) then
             if(present(f_id)) then
                E_IO = VTK_DAT_XML(var_location='node',var_block_action='close', cf=f_id)
             else
                E_IO = VTK_DAT_XML(var_location='node',var_block_action='close')
             endif
             f_vtk%mesh(nm)%status = pointdata_closed
          endif
          
          if(present(f_id)) then       
             E_IO = VTK_GEO_XML(cf=f_id)
             E_IO = VTK_END_XML(cf=f_id)
          else
             E_IO = VTK_GEO_XML()
             E_IO = VTK_END_XML()
          endif
          
          f_vtk%mesh(nm)%status = ended
       endif
    endif

  ! ----------------------------------------------------------------------------------
  end function write_VTK_end





  ! Write a single VTK file to disk (if I am fine MPI task)
  function write_VTK(f_vtk, f_name, n_part, t_step, n_mesh, o_fmt) result(E_IO)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(vtk_t),               intent(INOUT) :: f_vtk   !< VTK_t derived type
    character(len=*), optional, intent(IN)    :: f_name  !< VTK File NAME
    integer(ip),      optional, intent(IN)    :: n_part  !< Number of the PART
    real(rp),         optional, intent(IN)    :: t_step  !< Time STEP value
    integer(ip),      optional, intent(IN)    :: n_mesh  !< Number of the MESH
    character(len=*), optional, intent(IN)    :: o_fmt   !< Ouput ForMaT
    character(len=:), allocatable             :: fn      !< Real File Name
    character(len=:), allocatable             :: dp      !< Real Directory Path
    character(len=:), allocatable             :: of      !< Real Output Format
    real(rp)                                  :: ts      !< Real Time Step
    logical                                   :: ft      !< Fine Task
    integer(ip)                               :: nm      !< Real Number of the Mesh
    integer(ip)                               :: np      !< Real Number of the Part
    integer(ip)                               :: me      !< Task identifier
    integer(ip)                               :: fid     !< Real File ID
    integer(ip)                               :: E_IO    !< Error IO
  ! ----------------------------------------------------------------------------------

    me = 0; np = 1
    call f_vtk%env%info(me,np) 
    np = me
    if(present(n_part)) np = n_part

    nm = f_vtk%num_meshes
    if(present(n_mesh)) nm = n_mesh
    
    ts = 0._rp
    if(present(t_step)) ts = t_step 


    ft =  f_vtk%env%am_i_fine_task()
    E_IO = 0
    
    if(ft) then
       dp = f_vtk%get_VTK_time_output_path(f_path=f_vtk%mesh(nm)%dir_path, t_step=ts, n_mesh=nm)
       fn = f_vtk%get_VTK_filename(f_prefix=f_vtk%mesh(nm)%prefix, n_part=np, n_mesh=nm)
       fn = dp//fn
       if(present(f_name)) fn = f_name

       of = 'raw'
       if(present(o_fmt)) of = trim(adjustl(o_fmt))

       E_IO = f_vtk%write_VTK_start(fn, np, ts, nm, of, fid) 
       E_IO = f_vtk%write_VTK_unknowns(nm, fid)
       E_IO = f_vtk%write_VTK_end(nm, fid)
    end if

  ! ----------------------------------------------------------------------------------
  end function write_VTK


  ! Write the PVTK file referencing serveral parts of the partitioned mesh (if I am root_proc)
  function write_PVTK(f_vtk, f_name, n_mesh, t_step) result(E_IO)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(vtk_t),             intent(INOUT) :: f_vtk
    character(len=*), optional, intent(IN)    :: f_name
    integer(ip),      optional, intent(IN)    :: n_mesh
    real(rp),         optional, intent(IN)    :: t_step
    integer(ip)                               :: nm, rf
    character(len=:),allocatable              :: var_name
    character(len=:),allocatable              :: fn ,dp
    real(rp)                                  :: ts
    integer(ip)                               :: i, fid, nnods, nels, E_IO
    integer(ip)                               :: me, np
    logical                                   :: isDir
  ! ----------------------------------------------------------------------------------

    me = 0; np = 1
    check(associated(f_vtk%env))
    call f_vtk%env%info(me,np) 

    E_IO = 0

    if( f_vtk%env%am_i_fine_task() .and. me == f_vtk%root_proc) then

        nm = f_vtk%num_meshes
        if(present(n_mesh)) nm = n_mesh
        ts = 0_rp
        if(allocated(f_vtk%steps)) then
            if(f_vtk%steps_counter >0 .and. f_vtk%steps_counter <= size(f_vtk%steps,1)) &
                ts = f_vtk%steps(f_vtk%steps_counter)
        endif
        if(present(t_step)) ts = t_step 

        dp = f_vtk%get_VTK_time_output_path(f_path=f_vtk%mesh(nm)%dir_path, t_step=ts, n_mesh=nm)
        fn = f_vtk%get_PVTK_filename(f_prefix=f_vtk%mesh(nm)%prefix, n_mesh=nm, t_step=ts)
        fn = dp//fn
        if(present(f_name)) fn = f_name

        nnods = f_vtk%mesh(nm)%nnods
        nels = size(f_vtk%mesh(nm)%ctype, dim=1)

!        inquire( file=trim(dp)//'/.', exist=isDir ) 
        if(f_vtk%create_dir_hierarchy_on_root_process(trim(dp), issue_final_barrier=.False.) == 0) then
            ! pvtu
            E_IO = PVTK_INI_XML(filename = trim(adjustl(fn)), mesh_topology = 'PUnstructuredGrid', tp='Float64', cf=rf)
            do i=0, f_vtk%num_parts-1
                E_IO = PVTK_GEO_XML(source=trim(adjustl(f_vtk%get_VTK_filename(f_prefix=f_vtk%mesh(nm)%prefix, n_part=i, n_mesh=nm))), cf=rf)
            enddo
            if(allocated(f_vtk%mesh(nm)%unknowns).or.allocated(f_vtk%mesh(nm)%postprocess_fields)) then
               E_IO = PVTK_DAT_XML(var_location = 'Node', var_block_action = 'OPEN', cf=rf)
               ! Write unknowns point data
               if (allocated(f_vtk%mesh(nm)%unknowns)) then
                  do i=1, size(f_vtk%mesh(nm)%unknowns, dim=1)
                     if(f_vtk%mesh(nm)%unknowns(i)%filled .and. trim(adjustl(f_vtk%mesh(nm)%unknowns(i)%var_location)) == 'node') then
                        if(allocated(f_vtk%mesh(nm)%unknowns(i)%var_name)) then
                           var_name = f_vtk%mesh(nm)%unknowns(i)%var_name
                        else
                           var_name = 'Unknown_'//trim(adjustl(ch(i)))
                        endif
                        E_IO = PVTK_VAR_XML(varname = trim(adjustl(var_name)), tp=trim(adjustl(f_vtk%mesh(nm)%unknowns(i)%field_type)), Nc=f_vtk%mesh(nm)%unknowns(i)%num_comp , cf=rf)
                     endif
                  enddo
               end if
               ! Write postprocess fields point data
               if (allocated(f_vtk%mesh(nm)%postprocess_fields)) then
                  do i=1, size(f_vtk%mesh(nm)%postprocess_fields, dim=1)
                     if(f_vtk%mesh(nm)%postprocess_fields(i)%filled .and. trim(adjustl(f_vtk%mesh(nm)%postprocess_fields(i)%var_location)) == 'node') then
                        if(allocated(f_vtk%mesh(nm)%postprocess_fields(i)%var_name)) then
                           var_name = f_vtk%mesh(nm)%postprocess_fields(i)%var_name
                        else
                           var_name = 'Postprocess_field_'//trim(adjustl(ch(i)))
                        endif
                        E_IO = PVTK_VAR_XML(varname = trim(adjustl(var_name)), tp=trim(adjustl(f_vtk%mesh(nm)%postprocess_fields(i)%field_type)), Nc=f_vtk%mesh(nm)%postprocess_fields(i)%num_comp , cf=rf)
                     endif
                  enddo
               end if
               E_IO = PVTK_DAT_XML(var_location = 'Node', var_block_action = 'CLOSE', cf=rf)
               E_IO = PVTK_DAT_XML(var_location = 'Cell', var_block_action = 'OPEN', cf=rf)
               ! Write unknowns cell data
               if (allocated(f_vtk%mesh(nm)%unknowns)) then
                  do i=1, size(f_vtk%mesh(nm)%unknowns, dim=1)
                     if(f_vtk%mesh(nm)%unknowns(i)%filled .and. trim(adjustl(f_vtk%mesh(nm)%unknowns(i)%var_location)) == 'cell') then
                        if(allocated(f_vtk%mesh(nm)%unknowns(i)%var_name)) then
                           var_name = f_vtk%mesh(nm)%unknowns(i)%var_name
                        else
                           var_name = 'Unknown_'//trim(adjustl(ch(i)))
                        endif
                        E_IO = PVTK_VAR_XML(varname = trim(adjustl(var_name)), tp=trim(adjustl(f_vtk%mesh(nm)%unknowns(i)%field_type)), cf=rf)
                     endif
                  enddo
               end if
               ! Write postprocess fields cell data
               if (allocated(f_vtk%mesh(nm)%postprocess_fields)) then
                  do i=1, size(f_vtk%mesh(nm)%postprocess_fields, dim=1)
                     if(f_vtk%mesh(nm)%postprocess_fields(i)%filled .and. trim(adjustl(f_vtk%mesh(nm)%postprocess_fields(i)%var_location)) == 'cell') then
                        if(allocated(f_vtk%mesh(nm)%postprocess_fields(i)%var_name)) then
                           var_name = f_vtk%mesh(nm)%postprocess_fields(i)%var_name
                        else
                           var_name = 'Postprocess_fields_'//trim(adjustl(ch(i)))
                        endif
                        E_IO = PVTK_VAR_XML(varname = trim(adjustl(var_name)), tp=trim(adjustl(f_vtk%mesh(nm)%postprocess_fields(i)%field_type)), cf=rf)
                     endif
                  enddo
               end if
               E_IO = PVTK_DAT_XML(var_location = 'Cell', var_block_action = 'CLOSE', cf=rf)
            endif
            E_IO = PVTK_END_XML(cf=rf)
        endif

    endif

  ! ----------------------------------------------------------------------------------
  end function write_PVTK


  ! Write the PVD file referencing several PVTK files in a timeline (if I am root proc)
  function write_PVD(f_vtk, f_name, n_mesh) result(E_IO)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(vtk_t),             intent(INOUT) :: f_vtk
    character(len=*), optional, intent(IN)    :: f_name
    integer(ip),      optional, intent(IN)    :: n_mesh
    integer(ip)                               :: nm, rf
    character(len=:),allocatable              :: var_name
    character(len=:),allocatable              :: pvdfn, pvtkfn ,dp
    integer(ip)                               :: i, fid, nnods, nels, ts, E_IO
    integer(ip)                               :: me, np
    logical                                   :: isDir
  ! ----------------------------------------------------------------------------------

    me = 0; np = 1
    check(associated(f_vtk%env))
    call f_vtk%env%info(me,np) 

    E_IO = 0

    if(f_vtk%env%am_i_fine_task() .and. me == f_vtk%root_proc) then

        nm = f_vtk%num_meshes
        if(present(n_mesh)) nm = n_mesh
    
        pvdfn = trim(adjustl(f_vtk%mesh(nm)%dir_path))//'/'//trim(adjustl(f_vtk%mesh(nm)%prefix))//'_'//trim(adjustl(ch(nm)))//pvd_ext
        if(present(f_name)) pvdfn = f_name

!        inquire( file=trim(adjustl(f_vtk%mesh(nm)%dir_path))//'/.', exist=isDir )
        if(f_vtk%create_dir_hierarchy_on_root_process(trim(trim(adjustl(f_vtk%mesh(nm)%dir_path))), issue_final_barrier=.False.) == 0) then
            if(allocated(f_vtk%steps)) then
                if(size(f_vtk%steps,1) >= min(f_vtk%num_steps,f_vtk%steps_counter)) then
                    E_IO = PVD_INI_XML(filename=trim(adjustl(pvdfn)),cf=rf)
                    do ts=1, min(f_vtk%num_steps,f_vtk%steps_counter)
                        dp = f_vtk%get_PVD_time_output_path(f_path=f_vtk%mesh(nm)%dir_path, t_step=f_vtk%steps(ts))
                        pvtkfn = f_vtk%get_PVTK_filename(f_prefix=f_vtk%mesh(nm)%prefix, n_mesh=nm, t_step=f_vtk%steps(ts))
                        pvtkfn = dp//pvtkfn
                        E_IO = PVD_DAT_XML(filename=trim(adjustl(pvtkfn)),timestep=ts, cf=rf)
                    enddo
                    E_IO = PVD_END_XML(cf=rf)
                endif
            endif
        endif
    
    endif

  ! ----------------------------------------------------------------------------------
  end function write_PVD

  ! Read a single VTK file from disk (if I am fine MPI task)
  function read_VTK(f_vtk, f_name, n_part, t_step, n_mesh, o_fmt) result(E_IO)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(vtk_t),               intent(INOUT) :: f_vtk   !< VTK_t derived type
    character(len=*), optional, intent(IN)    :: f_name  !< VTK File NAME
    integer(ip),      optional, intent(IN)    :: n_part  !< Number of the PART
    real(rp),         optional, intent(IN)    :: t_step  !< Time STEP value
    integer(ip),      optional, intent(IN)    :: n_mesh  !< Number of the MESH
    character(len=*), optional, intent(IN)    :: o_fmt   !< Ouput ForMaT
    character(len=:), allocatable             :: fn      !< Real File Name
    character(len=:), allocatable             :: dp      !< Real Directory Path
    character(len=:), allocatable             :: of      !< Real Output Format
    real(rp)                                  :: ts      !< Real Time Step
    logical                                   :: ft      !< Fine Task
    integer(ip)                               :: nm      !< Real Number of the Mesh
    integer(ip)                               :: np      !< Real Number of the Part
    integer(ip)                               :: me      !< Task identifier
    integer(ip)                               :: fid     !< Real File ID
    integer(ip)                               :: E_IO    !< Error IO
  ! ----------------------------------------------------------------------------------

    me = 0; np = 1
    call f_vtk%env%info(me,np) 
    np = me
    if(present(n_part)) np = n_part

    nm = f_vtk%num_meshes
    if(present(n_mesh)) nm = n_mesh
    
    ts = 0._rp
    if(present(t_step)) ts = t_step 

    ft =  f_vtk%env%am_i_fine_task()

    E_IO = 0

    if(ft) then
        dp = f_vtk%get_VTK_time_output_path(f_path=f_vtk%mesh(nm)%dir_path, t_step=ts, n_mesh=nm)
        fn = f_vtk%get_VTK_filename(f_prefix=f_vtk%mesh(nm)%prefix, n_part=np, n_mesh=nm)
        fn = dp//fn
        if(present(f_name)) fn = f_name

        of = 'raw'
        if(present(o_fmt)) of = trim(adjustl(o_fmt))

        E_IO = f_vtk%read_VTK_start(fn, np, ts, nm, of, fid) 
        E_IO = f_vtk%read_VTK_unknowns(nm, fid)
        E_IO = f_vtk%read_VTK_end(nm, fid)
    endif

  ! ----------------------------------------------------------------------------------
  end function read_VTK

  ! Start the read of a single VTK file from disk (if I am fine MPI task)
  function read_VTK_start(f_vtk, f_name, n_part, t_step, n_mesh, o_fmt, f_id) result(E_IO)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(vtk_t),             intent(INOUT)   :: f_vtk   !< VTK_t derived type
    character(len=*), optional, intent(IN)    :: f_name  !< VTK File NAME
    integer(ip),      optional, intent(IN)    :: n_part  !< Number of the PART
    real(rp),         optional, intent(IN)    :: t_step  !< Time STEP value
    integer(ip),      optional, intent(IN)    :: n_mesh  !< Number of the MESH
    character(len=*), optional, intent(IN)    :: o_fmt   !< Ouput ForMaT
    integer(ip),      optional, intent(OUT)   :: f_id    !< File ID
    character(len=:), allocatable             :: fn      !< Real File Name
    character(len=:), allocatable             :: dp      !< Real Directory Path
    character(len=:), allocatable             :: of      !< Real Output Format
    real(rp)                                  :: ts      !< Real Time Step
    logical                                   :: ft      !< Fine Task
    integer(ip)                               :: nm      !< Real Number of the Mesh
    integer(ip)                               :: np      !< Real Number of the Part
    integer(ip)                               :: me      !< Task identifier
    integer(ip)                               :: fid     !< Real File ID
    integer(ip)                               :: nnods   !< Number of NODeS
    integer(ip)                               :: nels    !< Number of ELementS
    integer(ip)                               :: E_IO    !< Error IO
  ! ----------------------------------------------------------------------------------

    check(associated(f_vtk%env))
 
    ft =  f_vtk%env%am_i_fine_task() 

    E_IO = 0
    fid = -1
    
    if(ft) then
        me = 0; np = 1
        call f_vtk%env%info(me,np) 
        np = me
        if(present(n_part)) np = n_part

        nm = f_vtk%num_meshes
        if(present(n_mesh)) nm = n_mesh
    
        ts = 0._rp
        if(present(t_step)) ts = t_step 
        call f_vtk%append_step(ts)


        if(f_vtk%mesh(nm)%status == unknown .or. f_vtk%mesh(nm)%status == ended) then
    
            dp = f_vtk%get_VTK_time_output_path(f_path=f_vtk%mesh(nm)%dir_path, t_step=ts, n_mesh=nm)
            fn = f_vtk%get_VTK_filename(f_prefix=f_vtk%mesh(nm)%prefix, n_part=np, n_mesh=nm)
            fn = dp//fn
            if(present(f_name)) fn = f_name

            if( f_vtk%create_dir_hierarchy_on_root_process(dp,issue_final_barrier=.True.) == 0) then    
                of = 'raw'
                if(present(o_fmt)) of = trim(adjustl(o_fmt))
        
                nnods = f_vtk%mesh(nm)%nnods
                nels = size(f_vtk%mesh(nm)%ctype, dim=1)
        
                E_IO = VTK_INI_XML_READ(input_format = trim(adjustl(of)), filename = trim(adjustl(fn)), mesh_topology = 'UnstructuredGrid', cf=fid)
                f_vtk%mesh(nm)%status = started_read
            endif
        endif
        if(present(f_id)) f_id = fid

    endif

  ! ----------------------------------------------------------------------------------
  end function read_VTK_start

  ! Read a single VTK file from disk (if I am fine MPI task)
  function read_VTK_unknowns(f_vtk, n_mesh, f_id) result(E_IO)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(vtk_t),             intent(INOUT)   :: f_vtk          !< VTK_t derived type
    integer(ip),      optional, intent(IN)    :: n_mesh         !< Number of MESH
    integer(ip),      optional, intent(IN)    :: f_id           !< File ID
    real(rp), allocatable                     :: field(:,:)     !< FIELD(ncomp,nnod)
    integer(ip)                               :: nm             !< Real Number of Mesh
    integer(ip)                               :: E_IO           !< IO Error
    logical                                   :: ft             !< Fine Task

    check(associated(f_vtk%env))
    ft =  f_vtk%env%am_i_fine_task() 

    E_IO = 0
    
    if(ft) then        

       nm = f_vtk%num_meshes
       if(present(n_mesh)) nm = n_mesh

       if(f_vtk%mesh(nm)%status==started_read) then
          
          if(f_vtk%mesh(f_vtk%num_meshes)%linear_order .eqv. .true.) then
             E_IO = read_VTK_unknowns_linear(f_vtk, nm, f_id)
          else
             E_IO = read_VTK_unknowns_superlinear(f_vtk, nm, f_id)
          end if
       
       end if          

    end if

  end function read_VTK_unknowns

  ! Read a single VTK file to disk (if I am fine MPI task) for linear interpolation
  function read_VTK_unknowns_linear(f_vtk, nm, f_id) result(E_IO)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(vtk_t),             intent(INOUT)   :: f_vtk          !< VTK_t derived type
    integer(ip),                intent(IN)    :: nm             !< Real Number of Mesh
    integer(ip),      optional, intent(IN)    :: f_id           !< File ID
    real(rp), allocatable                     :: field(:,:)     !< FIELD(ncomp,nnod)
    integer(ip)                               :: tidx           !< Actual Time InDeX
    integer(ip)                               :: nnods          !< Number of NODes
    integer(ip)                               :: nels           !< Number of ELementS
    integer(ip)                               :: tnnod          !< Total Number of NODes
    integer(ip)                               :: curr_ncomp     !< CURRent Number of COMPonent
    integer(ip)                               :: tncomp         !< Total Number of COMPonents
    integer(ip)                               :: nnode          !< NODE Number
    integer(ip)                               :: elnnod         !< Number of NODes per ELement
    integer(ip)                               :: i, j, f, idx   !< Indices
    integer(ip)                               :: E_IO           !< IO Error
    logical                                   :: ft             !< Fine Task
  ! ----------------------------------------------------------------------------------
    
    tidx = 1

    if(.not. allocated(f_vtk%mesh(nm)%unknowns)) call f_vtk%fill_unknowns_from_physical_problem(nm)

    if(allocated(f_vtk%mesh(nm)%unknowns)) then

       nnods = f_vtk%mesh(nm)%nnods
       nels = size(f_vtk%mesh(nm)%ctype, dim=1)

       tncomp = 0
       do f=1, size(f_vtk%mesh(nm)%unknowns, dim=1) 
          if(f_vtk%mesh(nm)%unknowns(f)%filled) then
             tnnod = 0
             curr_ncomp = tncomp + f_vtk%mesh(nm)%unknowns(f)%num_comp
             !                   allocate(field(f_vtk%mesh(nm)%unknowns(f)%num_comp,nnods))
             call memalloc( f_vtk%mesh(nm)%unknowns(f)%num_comp, nnods, field, __FILE__,__LINE__)

             if(present(f_id)) then 
                E_IO = VTK_VAR_XML_READ(var_location='node',varname=f_vtk%mesh(nm)%unknowns(f)%var_name,NC_NN=nnods, &
                     &                  NCOMP=f_vtk%mesh(nm)%unknowns(f)%num_comp,var=field,cf=f_id)
             else
                E_IO = VTK_VAR_XML_READ(var_location='node',varname=f_vtk%mesh(nm)%unknowns(f)%var_name,NC_NN=nnods, &
                     &                  NCOMP=f_vtk%mesh(nm)%unknowns(f)%num_comp,var=field)
             end if

             do i=1, nels
                elnnod = f_vtk%fe_space%finite_elements(i)%reference_element_vars(1)%p%nvef_dim(2)-1 !Num nodes (dim=2 -> vertex)
                do j=1, elnnod
                   nnode = f_vtk%fe_space%finite_elements(i)%reference_element_vars(curr_ncomp)%p%ntxob%p(j)
                   idx = f_vtk%fe_space%finite_elements(i)%reference_element_vars(curr_ncomp)%p%ntxob%l(nnode)
                   f_vtk%fe_space%finite_elements(i)%unkno(idx, tncomp+1:tncomp+f_vtk%mesh(nm)%unknowns(f)%num_comp, tidx) = &
                        &  field(1:f_vtk%mesh(nm)%unknowns(f)%num_comp,j+tnnod) 
                        
                enddo
                tnnod = tnnod + elnnod 
             enddo

             tncomp = curr_ncomp
             if(allocated(field)) call memfree(field, __FILE__,__LINE__)
          endif
       enddo

    endif

  ! ----------------------------------------------------------------------------------
  end function read_VTK_unknowns_linear

  ! Write a single VTK file to disk (if I am fine MPI task) for high order interpolation
  function read_VTK_unknowns_superlinear(f_vtk, nm, f_id) result(E_IO)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(vtk_t),             intent(INOUT)   :: f_vtk          !< VTK_t derived type
    integer(ip),                intent(IN)    :: nm             !< Real Number of Mesh
    integer(ip),      optional, intent(IN)    :: f_id           !< File ID
    real(rp), allocatable                     :: field(:,:)     !< FIELD(ncomp,nnod)
    integer(ip)                               :: tidx           !< Actual Time InDeX
    integer(ip)                               :: nnods          !< Number of NODes
    integer(ip)                               :: nels           !< Number of ELementS
    integer(ip)                               :: nsubels        !< Number of sub_ELementS
    integer(ip)                               :: tnnod          !< Total Number of NODes
    integer(ip)                               :: curr_ncomp     !< CURRent Number of COMPonent
    integer(ip)                               :: tncomp         !< Total Number of COMPonents
    integer(ip)                               :: nnode, gnode   !< NODE Number
    integer(ip)                               :: ndime          !< Physical dimensions
    integer(ip)                               :: elnnod         !< Number of NODes per ELement
    integer(ip)                               :: i, j, f, idx   !< Indices
    integer(ip)                               :: subelem, icomp !< Indices
    integer(ip)                               :: order          !< Order of interpolation
    integer(ip)                               :: E_IO           !< IO Error
    logical                                   :: ft             !< Fine Task
    integer(ip)                               :: f_type, pos_elinf, istat, v_key, ltype(2), inode, jnode
    real(rp), allocatable                     :: target_field(:,:)
  ! ----------------------------------------------------------------------------------
    
    tidx = 1

    if(.not. allocated(f_vtk%mesh(nm)%unknowns)) call f_vtk%fill_unknowns_from_physical_problem(nm)

    if(allocated(f_vtk%mesh(nm)%unknowns) ) then

       nnods = f_vtk%mesh(nm)%nnods
       nels = f_vtk%fe_space%g_trian%num_elems
       ndime = f_vtk%fe_space%g_trian%num_dims
       f_type = Q_type_id

       tncomp = 0
       do f=1, size(f_vtk%mesh(nm)%unknowns, dim=1) 
          if(f_vtk%mesh(nm)%unknowns(f)%filled) then
             tnnod = 0
             curr_ncomp = tncomp + f_vtk%mesh(nm)%unknowns(f)%num_comp
             !                   allocate(field(f_vtk%mesh(nm)%unknowns(f)%num_comp,nnods))
             call memalloc( f_vtk%mesh(nm)%unknowns(f)%num_comp, nnods, field, __FILE__,__LINE__)

             if(present(f_id)) then 
                E_IO = VTK_VAR_XML_READ(var_location='node',varname=f_vtk%mesh(nm)%unknowns(f)%var_name,NC_NN=nnods, &
                     &                  NCOMP=f_vtk%mesh(nm)%unknowns(f)%num_comp,var=field,cf=f_id)
             else
                E_IO = VTK_VAR_XML_READ(var_location='node',varname=f_vtk%mesh(nm)%unknowns(f)%var_name,NC_NN=nnods, &
                     &                  NCOMP=f_vtk%mesh(nm)%unknowns(f)%num_comp,var=field)
             end if

             do i=1, nels

                ! Check (only working for Q-type elements)
                check(f_vtk%fe_space%finite_elements(i)%p_geo_reference_element%ftype == Q_type_id)

                ! Set order of interpolation
                order = maxval(f_vtk%fe_space%finite_elements(i)%order)
                nsubels = Q_nnods(ndime,order-1)
                gnode   = Q_nnods(ndime,1)
                nnode   = Q_nnods(ndime,order)

                ! Loop over subelements
                do subelem = 1,nsubels

                   ! Loop over geometrical nodes in subelement
                   do j=1, gnode
                      idx = f_vtk%mesh(nm)%nodes_subelem(order)%a(j,subelem)
                       f_vtk%fe_space%finite_elements(i)%unkno(idx, tncomp+1:tncomp+f_vtk%mesh(nm)%unknowns(f)%num_comp, tidx) = &
                           field(1:f_vtk%mesh(nm)%unknowns(f)%num_comp,j+tnnod)
                   enddo
                   tnnod = tnnod + gnode
                end do
                
                ! Interpolate unknowns with different order of interpolation
                if(f_vtk%fe_space%finite_elements(i)%order(curr_ncomp).ne.order) then

                   v_key = ndime + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)*order
                   call f_vtk%fe_space%pos_elem_info%get(key=v_key,val=pos_elinf,stat=istat)
                   call memalloc(f_vtk%mesh(nm)%unknowns(f)%num_comp,nnode,target_field,__FILE__,__LINE__)

                   do icomp=1,f_vtk%mesh(nm)%unknowns(f)%num_comp
                      target_field(icomp,:) = f_vtk%fe_space%finite_elements(i)%unkno(:,tncomp+icomp,tidx)
                   end do

                   do icomp=1,f_vtk%mesh(nm)%unknowns(f)%num_comp
                      do inode=1,f_vtk%fe_space%finite_elements(i)%reference_element_vars(curr_ncomp)%p%nnode
                         jnode = f_vtk%fe_space%finite_elements_info(pos_elinf)%ntxob%p(inode)
                         idx   = f_vtk%fe_space%finite_elements_info(pos_elinf)%ntxob%l(jnode)
                         f_vtk%fe_space%finite_elements(i)%unkno(inode,tncomp+icomp,tidx) = target_field(icomp,idx)
                      end do
                   end do

                   call memfree(target_field,__FILE__,__LINE__)

                end if


             enddo

             tncomp = curr_ncomp
             if(allocated(field)) call memfree(field, __FILE__,__LINE__)
          endif
       enddo

    endif

  ! ----------------------------------------------------------------------------------
  end function read_VTK_unknowns_superlinear

  ! Read a single VTK file from disk (if I am fine MPI task)
  function read_VTK_end(f_vtk, n_mesh, f_id) result(E_IO)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(vtk_t),               intent(INOUT) :: f_vtk   !< VTK_t derived type
    integer(ip),      optional, intent(IN)    :: n_mesh  !< Number of MESH
    integer(ip),      optional, intent(INOUT) :: f_id    !< File ID
    integer(ip)                               :: nm      !< Real Number of Mesh
    integer(ip)                               :: E_IO    !< IO Error
    logical                                   :: ft      !< Fine Task
  ! ----------------------------------------------------------------------------------

    check(associated(f_vtk%env))
    ft =  f_vtk%env%am_i_fine_task() 

    E_IO = 0

    nm = f_vtk%num_meshes
    if(present(n_mesh)) nm = n_mesh
    
    if(ft .and. (f_vtk%mesh(nm)%status >= started_read) .and. (f_vtk%mesh(nm)%status /= ended)) then

        if(present(f_id)) then       
            E_IO = VTK_END_XML_READ(cf=f_id)
        else
            E_IO = VTK_END_XML_READ()
        endif

        f_vtk%mesh(nm)%status = ended
    endif

  ! ----------------------------------------------------------------------------------
  end function read_VTK_end


  function create_dir_hierarchy_on_root_process(f_vtk, path, issue_final_barrier) result(res)
  ! ----------------------------------------------------------------------------------
    class(vtk_t),   intent(INOUT)   :: f_vtk
    character(len=*),  intent(IN)    :: path
    logical, optional, intent(IN)   :: issue_final_barrier
    logical                         :: ft, ifb
    integer(kind=c_int)             :: res
    integer(ip)                     :: me, np
  ! ----------------------------------------------------------------------------------
    me = 0; np = 1; ft = .False.; ifb = .False.
    if(present(issue_final_barrier)) ifb = issue_final_barrier
    check(associated(f_vtk%env))
    call f_vtk%env%info(me,np) 
    check(f_vtk%root_proc <= np-1)

    res=0
    if(me == f_vtk%root_proc) then
       res = mkdir_recursive(path//C_NULL_CHAR)
       check ( res == 0 ) 
    end if

    if(ifb) call f_vtk%env%first_level_barrier()

  ! ----------------------------------------------------------------------------------
  end function create_dir_hierarchy_on_root_process


  ! Build time output dir path for the vtk files in each timestep
  function get_VTK_time_output_path(f_vtk, f_path, t_step, n_mesh) result(dp)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(vtk_t),             intent(INOUT) :: f_vtk
    character(len=*), optional, intent(IN)    :: f_path
    real(rp),         optional, intent(IN)    :: t_step
    integer(ip),      optional, intent(IN)    :: n_mesh
    character(len=:), allocatable             :: dp
    character(len=:), allocatable             :: fp
    integer(ip)                               :: nm
    real(rp)                                  :: ts
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
    class(vtk_t),             intent(INOUT)   :: f_vtk
    character(len=*), optional, intent(IN)    :: f_path
    real(RP),         optional, intent(IN)    :: t_step
    character(len=:), allocatable             :: dp
    character(len=:), allocatable             :: fp
    real(rp)                                  :: ts
  ! ----------------------------------------------------------------------------------

    ts = f_vtk%num_steps
    if(present(t_step)) ts = t_step 

    dp = time_prefix//trim(adjustl(ch(ts)))//'/'

  end function get_PVD_time_output_path


  ! Build VTK filename
  function get_VTK_filename(f_vtk, f_prefix, n_part, n_mesh) result(fn)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(vtk_t),             intent(INOUT) :: f_vtk
    character(len=*), optional, intent(IN)    :: f_prefix
    integer(ip),      optional, intent(IN)    :: n_part
    integer(ip),      optional, intent(IN)    :: n_mesh
    character(len=:), allocatable             :: fn
    character(len=:), allocatable             :: fp
    integer(ip)                               :: nm, me, np
  ! ----------------------------------------------------------------------------------

    nm = f_vtk%num_meshes
    if(present(n_mesh)) nm = n_mesh

    me = 0; np = 1
    call f_vtk%env%info(me, np)
    np = me
    if(present(n_part)) np = n_part

    fp = f_vtk%mesh(nm)%prefix
    if(present(f_prefix)) fp = f_prefix

    fn = trim(adjustl(fp))//'_'//trim(adjustl(ch(nm)))//'_'//trim(adjustl(ch(np)))//vtk_ext

  end function get_VTK_filename


  ! Build VTK filename
  function get_PVTK_filename(f_vtk, f_prefix, n_mesh, t_step) result(fn)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(vtk_t),             intent(INOUT) :: f_vtk
    character(len=*), optional, intent(IN)    :: f_prefix
    integer(ip),      optional, intent(IN)    :: n_mesh
    real(rp),         optional, intent(IN)    :: t_step
    character(len=:), allocatable             :: fn
    character(len=:), allocatable             :: fp
    integer(ip)                               :: nm
    real(rp)                                  :: ts
  ! ----------------------------------------------------------------------------------

    nm = f_vtk%num_meshes
    if(present(n_mesh)) nm = n_mesh

    ts = 0._rp

    if(allocated(f_vtk%steps)) then
        if(f_vtk%steps_counter >0 .and. f_vtk%steps_counter <= size(f_vtk%steps,1)) &
            ts = f_vtk%steps(f_vtk%steps_counter)
    endif
    if(present(t_step)) ts = t_step

    fp = f_vtk%mesh(nm)%prefix
    if(present(f_prefix)) fp = f_prefix

    fn = trim(adjustl(fp))//'_'//trim(adjustl(ch(nm)))//'_'//trim(adjustl(ch(ts)))//pvtk_ext

  end function get_PVTK_filename


  ! Set the number of time steps of the simulation to be writen in the PVD
  subroutine set_root_proc(f_vtk, root)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(vtk_t), intent(INOUT) :: f_vtk
    integer(ip),    intent(IN)    :: root
  ! ----------------------------------------------------------------------------------

    f_vtk%root_proc = root
  ! ----------------------------------------------------------------------------------
  end subroutine set_root_proc


  ! Set the number of time steps of the simulation to be writen in the PVD
  subroutine set_num_steps(f_vtk, t_steps)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(vtk_t), intent(INOUT) :: f_vtk
    integer(ip),    intent(IN)  :: t_steps
    real(rp), allocatable       :: aux_steps(:)
  ! ----------------------------------------------------------------------------------

    f_vtk%num_steps = t_steps
    if(.not.allocated(f_vtk%steps)) then
        call memalloc ( f_vtk%num_steps, f_vtk%steps, __FILE__,__LINE__)
    elseif(size(f_vtk%steps)<f_vtk%num_steps) then
        call memalloc ( size(f_vtk%steps,1), aux_steps, __FILE__,__LINE__)
        aux_steps(1:size(f_vtk%steps,1)) = f_vtk%steps(1:size(f_vtk%steps,1))
        if (allocated(f_vtk%steps)) call memfree(f_vtk%steps, __FILE__,__LINE__)
        call memalloc ( f_vtk%num_steps, f_vtk%steps, __FILE__,__LINE__)
        f_vtk%steps(1:size(aux_steps,1)) = aux_steps(1:size(aux_steps,1))
        if (allocated(aux_steps)) call memfree(aux_steps, __FILE__,__LINE__)
    endif        

  ! ----------------------------------------------------------------------------------
  end subroutine set_num_steps

  ! Set the number of time steps of the simulation to be writen in the PVD
  subroutine append_step(f_vtk, c_step)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(vtk_t), intent(INOUT) :: f_vtk
    real(rp), intent(IN)        :: c_step !Current time step
    real(rp), allocatable       :: aux_steps(:)
  ! ----------------------------------------------------------------------------------

    f_vtk%steps_counter = f_vtk%steps_counter + 1
    if(.not.allocated(f_vtk%steps)) then
        call memalloc ( f_vtk%num_steps, f_vtk%steps, __FILE__,__LINE__)
    elseif(size(f_vtk%steps)<max(f_vtk%num_steps,f_vtk%steps_counter)) then
        call memalloc ( size(f_vtk%steps,1), aux_steps, __FILE__,__LINE__)
        aux_steps(1:size(f_vtk%steps,1)) = f_vtk%steps(1:size(f_vtk%steps,1))
        if (allocated(f_vtk%steps)) call memfree(f_vtk%steps, __FILE__,__LINE__)
        call memalloc ( max(f_vtk%num_steps,f_vtk%steps_counter), f_vtk%steps, __FILE__,__LINE__)
        f_vtk%steps(1:size(aux_steps,1)) = aux_steps(1:size(aux_steps,1))
        if (allocated(aux_steps)) call memfree(aux_steps, __FILE__,__LINE__)
    endif        

    if (f_vtk%num_steps < f_vtk%steps_counter) f_vtk%num_steps = f_vtk%steps_counter

    f_vtk%steps(f_vtk%steps_counter) = c_step

  ! ----------------------------------------------------------------------------------
  end subroutine append_step

  ! Set the number of parts of the partitioned mesh to be writen in the PVTK
  subroutine set_num_parts(f_vtk, n_parts)
  ! ----------------------------------------------------------------------------------
    class(vtk_t), intent(INOUT) :: f_vtk
    integer(ip),    intent(IN)    :: n_parts
  ! ----------------------------------------------------------------------------------

    f_vtk%num_parts = n_parts
  ! ----------------------------------------------------------------------------------
  end subroutine set_num_parts


  ! Set the name of the output directory
  subroutine set_dir_path(f_vtk, dir_path, nmesh)
  ! ----------------------------------------------------------------------------------
    class(vtk_t),        intent(INOUT) :: f_vtk
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
    class(vtk_t),        intent(INOUT) :: f_vtk
    character(len=*),      intent(IN)    :: prefix
    integer(ip), optional, intent(IN)    :: nmesh
    integer(ip)                          :: nm = 1
  ! ----------------------------------------------------------------------------------

    if(present(nmesh)) nm = nmesh
    f_vtk%mesh(nm)%prefix = prefix
  ! ----------------------------------------------------------------------------------
  end subroutine set_prefix


  ! Free the vtk_t derived type
  subroutine free (f_vtk)
  ! ----------------------------------------------------------------------------------
    implicit none
    class(vtk_t), intent(inout) :: f_vtk
    integer(ip)                   :: i, j
    logical                       :: ft
  ! ----------------------------------------------------------------------------------

    check(associated(f_vtk%env))
    ft = f_vtk%env%am_i_fine_task() 

    if(ft) then

        if(allocated(f_vtk%mesh)) then
            do i=1, size(f_vtk%mesh,1)
                call memfree(f_vtk%mesh(i)%X, __FILE__,__LINE__)
                call memfree(f_vtk%mesh(i)%Y, __FILE__,__LINE__)
                call memfree(f_vtk%mesh(i)%Z, __FILE__,__LINE__)
                call memfree(f_vtk%mesh(i)%connec, __FILE__,__LINE__)
                call memfree(f_vtk%mesh(i)%offset, __FILE__,__LINE__)
                call memfree(f_vtk%mesh(i)%ctype, __FILE__,__LINE__)
                if (.not. f_vtk%mesh(i)%linear_order) then
                    do j=1, max_order
                        call allocatable_array_free(f_vtk%mesh(i)%nodes_subelem(j))
                    enddo
                endif
                if(allocated(f_vtk%mesh(i)%unknowns)) then
                    do j=1, size(f_vtk%mesh(i)%unknowns)
                        if(allocated(f_vtk%mesh(i)%unknowns(j)%var_location)) deallocate(f_vtk%mesh(i)%unknowns(j)%var_location)
                        if(allocated(f_vtk%mesh(i)%unknowns(j)%var_name)) deallocate(f_vtk%mesh(i)%unknowns(j)%var_name)
                        if(allocated(f_vtk%mesh(i)%unknowns(j)%field_type)) deallocate(f_vtk%mesh(i)%unknowns(j)%field_type)
                        f_vtk%mesh(i)%unknowns(j)%filled = .False.
                    enddo
                    deallocate(f_vtk%mesh(i)%unknowns)
                endif
                f_vtk%mesh(i)%filled = .False.
            enddo
            deallocate(f_vtk%mesh)
        endif

        if (allocated(f_vtk%steps)) call memfree(f_vtk%steps, __FILE__,__LINE__)

    endif
    f_vtk%num_meshes = 0
    f_vtk%num_steps = 0
    f_vtk%num_parts = 0
    f_vtk%root_proc = 0
    f_vtk%fe_space => NULL()
    f_vtk%env => NULL()
  ! ----------------------------------------------------------------------------------
  end subroutine free


end module lib_vtk_io_interface_names
