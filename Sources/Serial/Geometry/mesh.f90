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
module mesh_names
  use types_names
  use memor_names
  use list_types_names
  use hash_table_names
  !use materials_names
  use stdio_names
  use reference_fe_names
  use reference_fe_factory_names
  !use conditions_names
  use mesh_distribution_names
  use metis_interface_names
  use rcm_renumbering_names
  implicit none
# include "debug.i90"
  private

  integer(ip), parameter :: c_order = 0
  integer(ip), parameter :: z_order = 1
  integer(ip), parameter :: max_num_elem_types = 3

  integer(ip), target :: permu_2DP1(3) = (/ 1, 2, 3/)
  integer(ip), target :: permu_2DQ1(4) = (/ 1, 2, 4, 3/)
  integer(ip), target :: permu_3DP1(4) = (/ 1, 2, 3, 4/)
  integer(ip), target :: permu_3DPR(6) = (/ 1, 2, 3, 4, 5, 6/)
  integer(ip), target :: permu_3DQ1(8) = (/ 1, 2, 4, 3, 5, 6, 8, 7/)
  integer(ip), target :: permu_id  (8) = (/ 1, 2, 3, 4, 5, 6, 7, 8/)

  type mesh_t
     ! Sizes
     integer(ip)                :: &
          order=c_order,           &         ! GiD element order (c)
          nelty=1,                 &         ! Number of element types
          ndime,                   &         ! Number of space dimensions
          npoin,                   &         ! Number of nodes (vertices)
          nvefs,                   &         ! Number of vefs
          nelem,                   &         ! Number of elements
          nnode,                   &         ! Maximum number of nodes per element
          nboun,                   &         ! Number of boundary elements
          nnodb                              ! Maximum number of nodes per boundary element

     ! Elements
     integer(ip), allocatable ::  &
          pnods(:),               &         ! pointers to the lnods
          lnods(:),               &         ! list of vefs of each element
          legeo(:),               &         ! List of geometry (volume) each element lies in
          leset(:),               &         ! List of sets associated to each element
          pvefs(:),               &         ! pointers to the lvefs
          lvefs(:),               &         ! list of vefs of each element
          lvef_geo(:),            &         ! List of geometric entities (volume, surface, point) each vef lies in
          lvef_set(:)                       ! List of sets associated to each vef

     ! Booundary (understood as a subset of vefs on which the definition of a set and a geometry is relevant,
     ! e.g. an internal interface separating materials over which a force has to be computed)
     type(list_t)             ::  &
          bound                             ! boundary elements (vefs)      
     integer(ip), allocatable ::  &
          lbgeo(:),               &         ! List of geometric entities (volume, surface, point) each boundary lies in
          lbset(:)                          ! List of sets associated to each boundary

     ! Dual mesh (elements around vertices)
     integer(ip)              ::  &
          nelpo = 0                         ! Nonzero when created
     integer(ip), allocatable ::  &
          pelpo(:),               &
          lelpo(:)

     real(rp), allocatable ::     &
          coord(:,:)                         ! Vertex coordinates

     type(p_reference_fe_t)   ::  ref_fe_list(max_num_elem_types)
    
    contains
     ! JP-TODO: program get and set for variables.
     procedure, non_overridable :: to_dual              => mesh_to_dual_new
     procedure, non_overridable :: create_distribution  => create_mesh_distribution
     procedure, non_overridable :: generate_vefs        => mesh_generate_vefs
     procedure, non_overridable :: get_sizes            => mesh_get_sizes
     procedure, non_overridable :: move_cells           => mesh_move_cells
     procedure, non_overridable :: move_coordinates     => mesh_move_coordinates
     procedure, non_overridable :: get_coordinates      => mesh_get_coordinates
     procedure, non_overridable :: get_boundary         => mesh_get_boundary
     procedure, non_overridable :: free                 => mesh_free
     procedure, non_overridable :: read_from_unit       => mesh_read_from_unit
     procedure, non_overridable :: read_from_file       => mesh_read_from_file
     generic :: read => read_from_file, read_from_unit
     procedure, non_overridable, nopass :: compose_name => mesh_compose_name
     procedure, non_overridable :: write_file_for_postprocess => mesh_write_file_for_postprocess
  end type mesh_t

  ! Types
  public :: mesh_t

  ! Constants
  public :: c_order, z_order

  ! Functions
  public :: mesh_write_file, mesh_write_files, mesh_write_files_for_postprocess

contains

  !=============================================================================
  subroutine mesh_get_sizes(this,ndime,npoin,nnode,nelem)
    class(mesh_t), intent(inout) :: this
    integer(ip), intent(inout) :: ndime,npoin,nelem,nnode
    ndime=this%ndime  ! Number of space dimensions
    npoin=this%npoin  ! Number of nodes (vertices)
    nelem=this%nelem  ! Number of elements
    nnode=this%nnode  ! Maximum number of nodes per element
  end subroutine mesh_get_sizes
  !=============================================================================
  subroutine mesh_move_cells(this,pvefs,lvefs)
    class(mesh_t)           , intent(inout) :: this
    integer(ip), allocatable, intent(inout) :: pvefs(:), lvefs(:)
    call memmovealloc(this%pnods,pvefs,__FILE__,__LINE__)
    call memmovealloc(this%lnods,lvefs,__FILE__,__LINE__)
  end subroutine mesh_move_cells
  !=============================================================================
  subroutine mesh_move_coordinates(this,coord)
    class(mesh_t)        , intent(inout) :: this
    real(rp), allocatable, intent(inout) :: coord(:,:)
    call memmovealloc(this%coord,coord,__FILE__,__LINE__)
  end subroutine mesh_move_coordinates
  !=============================================================================
  function mesh_get_coordinates(this)
    class(mesh_t), target, intent(inout) :: this
    real(rp)     , pointer       :: mesh_get_coordinates(:,:)
    mesh_get_coordinates => this%coord
  end function mesh_get_coordinates
  !=============================================================================
  subroutine mesh_get_boundary(this,boundary,lbgeo,lbset)
    class(mesh_t), target   , intent(inout) :: this
    type(list_t), pointer   , intent(inout) :: boundary
    integer(ip) , pointer   , intent(inout) :: lbgeo(:), lbset(:)
    boundary => this%bound
    lbgeo => this%lbgeo
    lbset => this%lbset
    !call memmovealloc(this%lbgeo,lbgeo,__FILE__,__LINE__)
    !call memmovealloc(this%lbset,lbset,__FILE__,__LINE__)
  end subroutine mesh_get_boundary
  !=============================================================================
  ! subroutine mesh_move_ref_fe(this,ref_fe_index,ref_fe_list)
  !   class(mesh_t), intent(inout) :: this
  !   type(p_reference_fe_t), allocatable, intent(inout) :: ref_fe_list(:)
  !   integer(ip)           , allocatable, intent(inout) :: ref_fe_index(:)
  !   call memmovealloc(this%ref_fe_list, ref_fe_list ,__FILE__,__LINE__)
  !   call memmovealloc(this%ref_fe_index,ref_fe_index,__FILE__,__LINE__)
  ! end subroutine mesh_move_ref_fe

!=============================================================================
  subroutine mesh_to_dual_new(this)
    class(mesh_t), intent(inout)     :: this
    
    ! Local variables
    integer(ip)              :: inode, ipoin, ielem, size_lnods

    if(this%nelpo>0) return

    call memalloc (this%npoin+1, this%pelpo, __FILE__,__LINE__)
    size_lnods = this%pnods(this%nelem+1)-1 
    
    ! Compute the number of elements around each point
    this%pelpo=0
    do inode=1, size_lnods
       ipoin=this%lnods(inode)
       this%pelpo(ipoin+1)=this%pelpo(ipoin+1)+1
    end do
    
    ! Find the maximum number of elements around a point
    this%nelpo=0
    do ipoin=1,this%npoin
       this%nelpo=max(this%nelpo,this%pelpo(ipoin+1))
    end do
    
    ! Compute pointers to the starting position of the list
    ! of elements around each point
    this%pelpo(1)=1
    do ipoin=1,this%npoin
       this%pelpo(ipoin+1)=this%pelpo(ipoin+1)+this%pelpo(ipoin)
    end do

    ! Allocate lelpo and fill it
    call memalloc (this%pelpo(this%npoin+1), this%lelpo, __FILE__,__LINE__)

    ! Compute the list of elements around each point.
    ! pelpo is used instead of auxiliary work space.
    do ielem=1,this%nelem 
       do inode=this%pnods(ielem),this%pnods(ielem+1)-1 
          ipoin=this%lnods(inode)
          this%lelpo(this%pelpo(ipoin))=ielem
          this%pelpo(ipoin)=this%pelpo(ipoin)+1
       end do
    end do
    
    ! Recover pelpo
    do ipoin=this%npoin+1, 2, -1
       this%pelpo(ipoin)=this%pelpo(ipoin-1)
    end do
    this%pelpo(1) = 1

  end subroutine mesh_to_dual_new

  !=============================================================================
  subroutine mesh_free (msh)
    !-----------------------------------------------------------------------
    ! This routine generates deallocates a mesh
    !-----------------------------------------------------------------------
    implicit none
    class(mesh_t), intent(inout)  :: msh

    if (allocated(msh%pnods)) call memfree (msh%pnods,__FILE__,__LINE__)
    if (allocated(msh%lnods)) call memfree (msh%lnods,__FILE__,__LINE__)
    if (allocated(msh%legeo)) call memfree (msh%legeo,__FILE__,__LINE__)
    if (allocated(msh%leset)) call memfree (msh%leset,__FILE__,__LINE__)
    if (allocated(msh%pvefs)) call memfree (msh%pvefs,__FILE__,__LINE__)
    if (allocated(msh%lvefs)) call memfree (msh%lvefs,__FILE__,__LINE__)
    if (allocated(msh%lvef_geo)) call memfree (msh%lvef_geo,__FILE__,__LINE__)
    if (allocated(msh%lvef_set)) call memfree (msh%lvef_set,__FILE__,__LINE__)
    if (allocated(msh%coord))    call memfree (msh%coord,__FILE__,__LINE__)

    msh%ndime=0
    msh%npoin=0
    msh%nvefs=0
    msh%nelem=0
    msh%nnode=0
    msh%order=c_order
    msh%nelty=1

    call msh%bound%free()
    if (allocated(msh%lbgeo)) call memfree (msh%lbgeo,__FILE__,__LINE__)
    if (allocated(msh%lbset)) call memfree (msh%lbset,__FILE__,__LINE__)
    msh%nnodb=0

    if (allocated(msh%pelpo)) call memfree (msh%pelpo,__FILE__,__LINE__)
    if (allocated(msh%lelpo)) call memfree (msh%lelpo,__FILE__,__LINE__)
    msh%nelpo = 0

  end subroutine mesh_free

  !===============================================================================================
  subroutine mesh_copy(msh_old,msh_new)
    implicit none
    type(mesh_t), intent(in)    :: msh_old
    type(mesh_t), intent(inout) :: msh_new

    msh_new%ndime = msh_old%ndime
    msh_new%npoin = msh_old%npoin
    msh_new%nelem = msh_old%nelem
    msh_new%nnode = msh_old%nnode
    
    call memalloc(msh_new%nelem+1,msh_new%pnods,__FILE__,__LINE__)
    msh_new%pnods = msh_old%pnods

    call memalloc(msh_new%pnods(msh_new%nelem+1),msh_new%lnods,__FILE__,__LINE__)
    msh_new%lnods = msh_old%lnods

    if (allocated(msh_old%coord)) then
       call memalloc(number_space_dimensions,msh_new%npoin,msh_new%coord,__FILE__,__LINE__)
       msh_new%coord = msh_old%coord
    end if

  end subroutine mesh_copy
  
  !=============================================================================
  subroutine mesh_read_from_unit(msh,lunio)
    !------------------------------------------------------------------------
    !
    ! This routine reads a mesh writen by GiD according to fempar problem type.
    !
    !------------------------------------------------------------------------
    implicit none
    integer(ip)      , intent(in)  :: lunio
    class(mesh_t)     , intent(out) :: msh
    !logical, optional, intent(in)  :: permute_c2z
    integer(ip)     :: idime,ipoin,inode,ielem,nboun,iboun,vbound,nnode,nnodb ! ,istat,jpoin,jelem,iboun
    character(14)   :: dum1
    character(7)    :: dum2
    character(7)    :: dum3
    character(10)   :: dum4
    character(10)   :: dum5
    character(6)   :: dum6
    character(1000) :: tel
    integer(ip), allocatable :: lnods_aux(:)
    integer(ip), allocatable :: bound_list_aux(:)
    type(list_iterator_t)    :: bound_iterator
    integer(ip), pointer     :: permu(:)
    logical                  :: permute_c2z_

    ! if(present(permute_c2z)) then
    !    permute_c2z_ = permute_c2z
    ! else
    !    permute_c2z_ = .false.
    ! end if
    ! write(*,*) 'Permuting c2z:',permute_c2z_

    ! Read first line: "MESH dimension  3  order  0  types  1  elements        352  nodes        100  boundaries        144"
    ! Read first line: "MESH dimension  2  order  0  types  1  elements        100  nodes        121  boundaries         40"
    ! Read first line: "MESH dimension  2  order  0  types  1  elements          1  vertices          4  vefs          8
    read(lunio,'(a14,1x,i2, a7,1x,i2, a7,1x,i2, a10,1x,i10, a10,1x,i10, a6,1x,i10)') &
         & dum1,msh%ndime,dum2,msh%order,dum3,msh%nelty,dum4,msh%nelem, dum5,msh%npoin,dum6,nboun

    write(*,*) 'Read mesh with parameters:',msh%ndime,msh%order,msh%nelty,msh%nelem,msh%npoin,nboun

    ! Read nodes
    call memalloc(number_space_dimensions,msh%npoin,msh%coord,__FILE__,__LINE__)
    msh%coord = 0.0_rp
    do while(tel(1:5).ne.'coord')
       read(lunio,'(a)') tel
    end do
    read(lunio,'(a)') tel
    do while(tel(1:5).ne.'end c')
       read(tel,*) ipoin,(msh%coord(idime,ipoin),idime=1,msh%ndime)
       read(lunio,'(a)') tel
    end do

    ! Read elements' size (pnods)
    call memalloc(msh%nelem+1,msh%pnods,__FILE__,__LINE__)
    do while(tel(1:5).ne.'eleme')
       read(lunio,'(a)') tel
    end do
    read(lunio,'(a)') tel
    do while(tel(1:5).ne.'end e')
       read(tel,*) ielem,msh%pnods(ielem+1)
       read(lunio,'(a)') tel
    end do
    ! Transform length to header and get mesh%nnode
    msh%pnods(1) = 1
    msh%nnode    = 0
    do ielem = 2, msh%nelem+1
       msh%nnode = max(msh%nnode,msh%pnods(ielem))
       msh%pnods(ielem) = msh%pnods(ielem)+msh%pnods(ielem-1)
    end do

    ! Read elements
    call memalloc(msh%pnods(msh%nelem+1),msh%lnods,__FILE__,__LINE__)
    call memalloc(msh%nelem,msh%legeo,__FILE__,__LINE__)
    call memalloc(msh%nelem,msh%leset,__FILE__,__LINE__)
    call io_rewind(lunio)
    do while(tel(1:5).ne.'eleme')
       read(lunio,'(a)') tel
    end do
    read(lunio,'(a)') tel
    do while(tel(1:5).ne.'end e')
       read(tel,*) ielem,nnode,(msh%lnods(msh%pnods(ielem)-1+inode),inode=1,nnode),msh%leset(ielem),msh%legeo(ielem)
       read(lunio,'(a)') tel
    end do

    ! Read boundary elements' size (pnodb)
    call msh%bound%create(nboun)
    do while(tel(1:5).ne.'vefs')
       read(lunio,'(a)') tel
    end do
    read(lunio,'(a)') tel
    do while(tel(1:5).ne.'end v')
       read(tel,*) iboun, vbound
       call msh%bound%sum_to_pointer_index(iboun, vbound)
       read(lunio,'(a)') tel
    end do
    ! Transform length to header and get mesh%nnodb
    call msh%bound%calculate_header()
    msh%nnodb    = 0
    do iboun = 2, msh%bound%get_num_pointers()+1
       msh%nnodb = max(msh%nnodb,msh%bound%get_sublist_size(iboun))
    end do

    ! Read boundary elements
    call msh%bound%allocate_list_from_pointer()
    call memalloc(msh%bound%get_num_pointers(),msh%lbgeo,__FILE__,__LINE__)
    call memalloc(msh%bound%get_num_pointers(),msh%lbset,__FILE__,__LINE__)
    call io_rewind(lunio)
    do while(tel(1:5).ne.'vefs')
       read(lunio,'(a)') tel
    end do
    read(lunio,'(a)') tel
    do while(tel(1:5).ne.'end v')
       read(tel,*) iboun,nnodb
       allocate(bound_list_aux(msh%bound%get_sublist_size(iboun)))
       read(tel,*) iboun,nnodb, (bound_list_aux(inode),inode=1,nnodb),msh%lbset(iboun),msh%lbgeo(iboun)
       bound_iterator = msh%bound%create_iterator(iboun)
       do inode=1, nnodb
          call bound_iterator%set_current(bound_list_aux(inode))
          call bound_iterator%next()
       enddo
       deallocate(bound_list_aux)
       read(lunio,'(a)') tel
    end do

    ! Reordering (c to z) the nodes of the mesh, if needed
    !if(permute_c2z_) then
    if(msh%order==c_order) then
       msh%order= z_order
       call memalloc(msh%nnode, lnods_aux, __FILE__, __LINE__)
       do ielem = 1,msh%nelem
          nnode = msh%pnods(ielem+1) - msh%pnods(ielem)
          lnods_aux(1:nnode) = msh%lnods(msh%pnods(ielem):msh%pnods(ielem+1)-1)
          if(msh%ndime == 2) then
             if(nnode == 3)  then    ! Linear triangles (2DP1)
                permu => permu_2DP1
             elseif(nnode == 4) then ! Linear quadrilaterals(2DQ1)
                permu => permu_2DQ1
             end if
          elseif(msh%ndime == 3) then
             if(nnode == 4) then     ! Linear tetrahedra (3DP1)
                permu => permu_3DP1
             elseif(nnode == 8) then ! Linear hexahedra (3DQ1)
                permu => permu_3DQ1
             end if
          end if
          do inode = 1, nnode
             msh%lnods(msh%pnods(ielem)+inode-1) = lnods_aux(permu(inode))
          end do
       end do
       call memfree(lnods_aux,__FILE__,__LINE__)
    end if

  end subroutine mesh_read_from_unit

  !=============================================================================
  subroutine mesh_write_file (msh,lunio,title)
    !------------------------------------------------------------------------
    !
    ! This routine writes a mesh in the format defined by GiD fempar problem type.
    !
    !------------------------------------------------------------------------
    implicit none
    integer(ip)  , intent(in)           :: lunio
    class(mesh_t) , intent(in)          :: msh
    character(*) , intent(in), optional :: title

    integer(ip)                         :: ielem, idime, ipoin, inode, iboun, jboun
    integer(ip), allocatable            :: bound_list_aux(:)
    type(list_iterator_t)               :: bound_iterator


    ! Read first line: "MESH dimension  2  order  0  types  2  elements         86  nodes         63  boundaries         28"
    write(lunio,'(a14,1x,i2,a7,1x,i2,a7,1x,i2,a10,1x,i10,a7,1x,i10,a12,1x,i10)') &
         & 'MESH dimension',msh%ndime,'  order',msh%order,'  types',msh%nelty,'  elements', &
         & msh%nelem,'  nodes',msh%npoin,'  boundaries',msh%bound%get_num_pointers()

    ! Coordinates
    write(lunio,'(a)')'coordinates'
    assert(allocated(msh%coord))
    do ipoin=1,msh%npoin
       write(lunio,'(i10,3(1x,e16.8e3))') ipoin,(msh%coord(idime,ipoin),idime=1,msh%ndime)
    end do
    write(lunio,'(a)')'end coordinates'

    ! Elements
    write(lunio,'(a)')'elements'
    do ielem=1,msh%nelem
       write(lunio,'(i10,65(1x,i10))') ielem, msh%pnods(ielem+1)-msh%pnods(ielem),&
            &  msh%lnods(msh%pnods(ielem):msh%pnods(ielem+1)-1),msh%leset(ielem),msh%legeo(ielem)
    end do
    write(lunio,'(a)') 'end elements'

    ! Boundary elements
    write(lunio,'(a)')'boundaries'
    do iboun=1,msh%bound%get_num_pointers()
       allocate(bound_list_aux(msh%bound%get_sublist_size(iboun)))
       bound_iterator = msh%bound%create_iterator(iboun)
       do jboun=1, msh%bound%get_sublist_size(iboun)
          bound_list_aux(jboun) = bound_iterator%get_current()
          call bound_iterator%next()
       enddo
       write(lunio,'(i10,65(1x,i10))') iboun, bound_iterator%get_size(), &
            &  bound_list_aux ,msh%lbset(iboun),msh%lbgeo(iboun)
       deallocate(bound_list_aux)
    end do
    write(lunio,'(a)') 'end boundaries'

  end subroutine mesh_write_file

  !=============================================================================
  subroutine mesh_write_post_file (msh,lunio,title)
    !------------------------------------------------------------------------
    !
    ! This routine writes a mesh in GiD format (only works for linear elements).
    !
    !------------------------------------------------------------------------
    implicit none
    integer(ip)      , intent(in)           :: lunio
    class(mesh_t)     , intent(in)           :: msh
    character(*)     , intent(in), optional :: title

    integer(ip)                    :: ielem, idime, ipoin, inode, nnode
    character(13)                  :: elemt
    character(len=:), allocatable  :: title_

    integer(ip)     , pointer      :: permu(:)

    permu => permu_id
    if(msh%ndime==2) then
       if(msh%nnode==3) then
          elemt='Triangle' 
          if(msh%order==z_order) permu => permu_2DP1
       else
          elemt='Quadrilateral'
          if(msh%order==z_order) permu => permu_2DQ1
       end if
    else
       if(msh%nnode==4) then 
          elemt='Tetrahedra'
          if(msh%order==z_order) permu => permu_3DP1
       else if(msh%nnode==6) then 
          elemt='Prism'
          if(msh%order==z_order) permu => permu_3DPR
       else if(msh%nnode==8) then 
          elemt='Hexahedra'
          if(msh%order==z_order) permu => permu_3DQ1
       end if
    end if

    ! Header
    title_ = 'TITLE'
    if(present(title)) title_=title
    write(lunio,1) adjustl(trim(title_)),msh%ndime,adjustl(trim(elemt)),msh%nnode

    ! Coordinates
    write(lunio,2)'coordinates'
    if (allocated(msh%coord)) then
       do ipoin=1,msh%npoin
          write(lunio,3) ipoin,(msh%coord(idime,ipoin),idime=1,msh%ndime)
       end do
    end if
    write(lunio,2)'end coordinates'

    ! Connectivity
    if(msh%nelty==1) then
       write(lunio,2)'elements'
       do ielem=1,msh%nelem
          nnode = msh%pnods(ielem+1)-msh%pnods(ielem)
          write(lunio,4) ielem, (msh%lnods(msh%pnods(ielem)-1+permu(inode)),inode=1,nnode),1
       end do
       write(lunio,2) 'end elements'
    else
       ! Write hexahedra or prismas (3D) or quads(2)
       write(lunio,2)'elements'
       do ielem=1,msh%nelem
          nnode = msh%pnods(ielem+1)-msh%pnods(ielem)
          if(nnode == msh%nnode) &
               write(lunio,4) ielem, (msh%lnods(msh%pnods(ielem)-1+permu(inode)),inode=1,nnode),1
       end do
       write(lunio,2) 'end elements'
       ! Now write tetrahedra (3D) or triangles (2D)
       if(msh%ndime==2) then
          nnode = 3
          elemt = 'Triangle' 
          if(msh%order==z_order) permu => permu_2DP1
       else if(msh%ndime==3) then
          nnode = 4
          elemt='Tetrahedra'
          if(msh%order==z_order) permu => permu_3DP1
       end if
       write(lunio,1) adjustl(trim(title_)),msh%ndime,adjustl(trim(elemt)),nnode
       write(lunio,2)'coordinates'
       write(lunio,2)'end coordinates'
       write(lunio,2)'elements'
       do ielem=1,msh%nelem
          if(msh%pnods(ielem+1)-msh%pnods(ielem) == nnode) &
               write(lunio,4) ielem, (msh%lnods(msh%pnods(ielem)-1+permu(inode)),inode=1,nnode),1
       end do
       write(lunio,2) 'end elements'
       ! Eventually write prismas (3D)
       if(msh%ndime==3.and.msh%nnode==8) then
          nnode = 4
          elemt='Prism'
          if(msh%order==z_order) permu => permu_3DPR
          write(lunio,1) adjustl(trim(title_)),msh%ndime,adjustl(trim(elemt)),nnode
          write(lunio,2)'coordinates'
          write(lunio,2)'end coordinates'
          write(lunio,2)'elements'
          do ielem=1,msh%nelem
             if(msh%pnods(ielem+1)-msh%pnods(ielem) == nnode) &
                  write(lunio,4) ielem, (msh%lnods(msh%pnods(ielem)-1+permu(inode)),inode=1,nnode),1
          end do
          write(lunio,2) 'end elements'
       end if
    end if

1   format('MESH ',a,' dimension ',i1,' Elemtype ',a,' Nnode ',i2)
2   format(a)
3   format(i10, 3(1x,e16.8e3))
4   format(i10,65(1x,i10))
5   format('BOUNDARY ',a,' Nnodb ',i2)
6   format(i6,10(1x,i6))

  end subroutine mesh_write_post_file

  !=============================================================================
  subroutine mesh_compose_name ( prefix, name ) 
    implicit none
    character(len=*)             , intent(in)    :: prefix 
    character(len=:), allocatable, intent(inout) :: name
    name = trim(prefix) // '.mesh'
  end subroutine mesh_compose_name
  !=============================================================================
  subroutine mesh_compose_post_name ( prefix, name ) 
    implicit none
    character(len=*)             , intent(in)    :: prefix 
    character(len=:), allocatable, intent(inout) :: name
    name = trim(prefix) // '.post.msh'
  end subroutine mesh_compose_post_name

  !=============================================================================
  subroutine mesh_write_files ( dir_path, prefix, nparts, lmesh )
     implicit none
     ! Parameters 
     character(*)   , intent(in)  :: dir_path 
     character(*)   , intent(in)  :: prefix
     integer(ip)     , intent(in)  :: nparts
     type(mesh_t)  , intent(in)  :: lmesh (nparts)

     character(len=:), allocatable :: name, rename ! Deferred-length allocatable character arrays

     ! Locals 
     integer (ip) :: i,lunio

     call mesh_compose_name ( prefix, name )

     do i=nparts, 1, -1  
        rename=name
        call numbered_filename_compose(i,nparts,rename)
        lunio = io_open( trim(dir_path) // '/' // trim(rename), 'write' )
        call mesh_write_file(lmesh(i),lunio)
        call io_close(lunio)
     end do
     
     ! name, and rename should be automatically deallocated by the compiler when they
     ! go out of scope. Should we deallocate them explicitly for safety reasons?
   end subroutine mesh_write_files

  !=============================================================================
  subroutine mesh_write_files_for_postprocess ( dir_path, prefix, nparts, lmesh )
     implicit none
     ! Parameters 
     character(*)    , intent(in)  :: dir_path 
     character(*)    , intent(in)  :: prefix
     integer(ip)     , intent(in)  :: nparts
     type(mesh_t)    , intent(in)  :: lmesh (nparts)
     character(len=:), allocatable :: name, rename ! Deferred-length allocatable character arrays

     ! Locals 
     integer (ip) :: i,lunio

     do i=nparts, 1, -1  
        name=prefix
        call numbered_filename_compose(i,nparts,name)
        call mesh_compose_post_name (name, rename)
        lunio = io_open( trim(dir_path) // '/' // trim(rename), 'write' )
        call mesh_write_post_file(lmesh(i),lunio)
        call io_close(lunio)
     end do
     
   end subroutine mesh_write_files_for_postprocess

  !=============================================================================
   subroutine mesh_read_files ( dir_path, prefix, nparts, lmesh )
     implicit none
     ! Parameters 
     character(*), intent(in)       :: dir_path 
     character(*), intent(in)       :: prefix
     integer(ip)   , intent(in)     :: nparts
     type(mesh_t), intent(out)    :: lmesh (nparts)
     character(len=:), allocatable  :: name, rename ! Deferred-length allocatable character arrays
     
     ! Locals 
     integer (ip)                     :: i,lunio
     
     call mesh_compose_name ( prefix, name )
     do i=nparts, 1, -1  
        rename=name
        call numbered_filename_compose(i,nparts,rename)
        lunio = io_open( trim(dir_path) // '/' // trim(rename), 'read' )
        call lmesh(i)%read(lunio)
        call io_close(lunio)
     end do
     
   end subroutine mesh_read_files

  !=============================================================================
   subroutine mesh_read_from_file (f_mesh,  dir_path, prefix ) !, permute_c2z
     implicit none 
     ! Parameters
     character (*)                , intent(in)  :: dir_path
     character (*)                , intent(in)  :: prefix
     class(mesh_t)                , intent(out) :: f_mesh
     !logical, optional, intent(in)  :: permute_c2z
     ! Locals
     integer(ip)                    :: lunio
     character(len=:), allocatable  :: name

     ! Read mesh
     call mesh_compose_name ( prefix, name )
     lunio = io_open( trim(dir_path)//'/'//trim(name), 'read', status='old' )
     call f_mesh%read(lunio)  !, permute_c2z
     call io_close(lunio)
   end subroutine mesh_read_from_file

   !=============================================================================
   subroutine mesh_write_file_for_postprocess ( f_mesh, dir_path, prefix)
     implicit none 
     ! Parameters
     character (*)                , intent(in)  :: dir_path
     character (*)                , intent(in)  :: prefix
     class(mesh_t)                , intent(in)  :: f_mesh

     ! Locals
     integer(ip)                    :: lunio
     character(len=:), allocatable  :: name

     call mesh_compose_post_name ( prefix, name )
     lunio = io_open( trim(dir_path)//'/'//trim(name), 'write' )
     call mesh_write_post_file(f_mesh,lunio)
     call io_close(lunio)

   end subroutine mesh_write_file_for_postprocess

  !==================================================================================================
  ! subroutine mesh_generate_vefs_gid_and_dim(mesh,num_global_cells,num_global_verices,cells_gid,vertx_gid,vef_gid)
  !   implicit none
  !   ! Parameters
  !   class(mesh_t), intent(inout) :: mesh
  !   integer(igp) , intent(inout) :: num_global_cells,num_global_verices
  !   integer(igp) , intent(inout) :: cells_gid(:)
  !   integer(igp) , intent(inout) :: vertx_gid(:)
  !   integer(igp) , allocatable, intent(inout) :: vef_gid(:)
  !   integer(ip) :: ielem, ivef, num_vertices

  !   call memalloc(mesh%pvefs(nelem+1),lst_vefs_gid,__FILE__,__LINE__)

  !   do ielem=1,mesh%nelem
  !      num_vertices = mesh%pnods(ielem+1)-mesh%pnods(ielem)
  !      !num_vefs     = mesh%pvefs(ielem+1)-mesh%pvefs(ielem)
  !      call pos_ref_fe%get(key=num_vertices,val=ielem_type,stat=istat)
  !      assert(istat==old_index)
  !      num_edges = mesh%ref_fe_list(ielem_type)%p%get_number_vefs_of_dimension(1)
  !      num_faces = mesh%ref_fe_list(ielem_type)%p%get_number_vefs_of_dimension(2)
  !      index=mesh%pvefs(ielem)-1
  !      do ivef=1,num_vertices
  !         index=index+1
  !         vef_lid = mesh%lvefs(index)
  !         vef_gid(index) = vertx_gid(vef_lid)
  !         vef_dim(index) = 0
  !      end do
  !      do ivef=num_vertices+1,num_vertices+num_edges
  !         vef_lid = mesh%lvefs(index)
  !         vef_gid(index) = ishft(int(ipart,igp),int(32,igp)) + int(vef_lid, igp) + ishft(int(1,igp),int(60,igp))
  !         vef_dim(index) = 1
  !      end do
  !      do ivef=num_vertices+num_edges+1,num_vertices+num_edges+num_faces
  !         vef_lid = mesh%lvefs(index)
  !         vef_gid(index) = ishft(int(ipart,igp),int(32,igp)) + int(vef_lid, igp) + ishft(int(1,igp),int(60,igp))
  !         vef_dim(index) = 2
  !      end do
  !   end do

  ! end subroutine mesh_generate_vefs_gid_and_dim
  !==================================================================================================
  subroutine mesh_generate_vefs(mesh)
    implicit none
    ! Parameters
    class(mesh_t), intent(inout) :: mesh
    ! integer(igp) , optional, intent(inout) :: nelem_global, npoin_global 
    ! integer(igp) , optional, intent(inout) :: cells_gid(:)
    ! integer(igp) , optional, intent(inout) :: vertx_gid(:)
    ! integer(igp) , optional, intent(inout) :: vef_gid(:)

    type(list_t), pointer    :: vertices_ivef
    type(list_t), pointer    :: vertices_jvef
    logical     :: equal
    integer(ip) :: istat, count, iboun, ivert, jvert, nnodb
    integer(ip) :: ielem, ivef, ielem_type, ielem_num_nodes
    integer(ip) :: ielem_num_vefs, ielem_first_vef_id, ielem_num_vef_verts
    integer(ip) :: vertex_of_ivef(4)
    integer(ip) :: jelpo, jelem, jvef, jelem_type, jelem_num_nodes
    integer(ip) :: jelem_num_vefs, jelem_first_vef_id, jelem_num_vef_verts
    integer(ip) :: vertex_of_jvef(4)
    type(list_iterator_t) :: vertices_ivef_iterator
    type(list_iterator_t) :: vertices_jvef_iterator
    type(list_iterator_t) :: bound_iterator

    ! We only require linear continuous elements so we could only have
    ! tetrahedra, hexahedra and prisms combined.
    type(position_hash_table_t) :: pos_ref_fe

    ! if(present(nelem_global)) then
    !    assert(present(npoin_global))
    !    assert(present(cells_gid))
    !    assert(present(vertx_gid))
    !    assert(present(vef_gid))
    !    assert(size(cells_gid)==mesh%nelem)
    !    assert(size(vertx_gid)==mesh%npoin)
    !    assert(size(vef_gid)==mesh%pvefs(mesh%nelem+1))
    ! end if

    ! Build reference_fe list
    call pos_ref_fe%init(max_num_elem_types)
    do ielem=1,mesh%nelem
       ielem_num_nodes=mesh%pnods(ielem+1)-mesh%pnods(ielem)
       call pos_ref_fe%get(key=ielem_num_nodes,val=ielem_type,stat=istat)
       if ( istat == new_index ) then
          if(mesh%ndime==2) then
             if(ielem_num_nodes==3) then ! Triangle
                assert(.false.)
             else ! Quadrilateral
                mesh%ref_fe_list(ielem_type) = &
                     &    make_reference_fe ( topology = topology_quad, fe_type = fe_type_lagrangian, &
                     &                        number_dimensions = mesh%ndime, order = 1,          &
                     &                        field_type = field_type_vector, continuity = .true. )
             end if
          else
             if(ielem_num_nodes==4) then ! Tetrahedra
                assert(.false.)
             else if(ielem_num_nodes==6) then ! Prism
                assert(.false.)
             else if(ielem_num_nodes==8) then ! Hexahedra
                assert(.false.)
             end if
          end if
       end if
    end do

    ! Compute mesh%pvefs and allocate mesh%lvefs
    call memalloc(mesh%nelem+1, mesh%pvefs, __FILE__, __LINE__ )
    mesh%pvefs(1)=1
    do ielem=1,mesh%nelem
       ielem_num_nodes=mesh%pnods(ielem+1)-mesh%pnods(ielem)
       call pos_ref_fe%get(key=ielem_num_nodes,val=ielem_type,stat=istat)
       assert(istat==old_index)
       ielem_num_vefs = ielem_num_nodes &
            &   + mesh%ref_fe_list(ielem_type)%p%get_number_vefs_of_dimension(1) 
       if(mesh%ndime==3) ielem_num_vefs = ielem_num_vefs &
            &   + mesh%ref_fe_list(ielem_type)%p%get_number_vefs_of_dimension(2)
       mesh%pvefs(ielem+1)=mesh%pvefs(ielem)+ielem_num_vefs
    end do
    call memalloc(mesh%pvefs(mesh%nelem+1), mesh%lvefs , __FILE__, __LINE__ )
    mesh%lvefs=0
    !write(*,*) mesh%pvefs

    ! Fill vefs
    mesh%nvefs = mesh%npoin
    do ielem=1,mesh%nelem
       ielem_num_nodes = mesh%pnods(ielem+1)-mesh%pnods(ielem)
       ! Fill vertices
       mesh%lvefs(mesh%pvefs(ielem):mesh%pvefs(ielem)+ielem_num_nodes-1)=mesh%lnods(mesh%pnods(ielem):mesh%pnods(ielem+1)-1)
       call pos_ref_fe%get(key=ielem_num_nodes,val=ielem_type,stat=istat)
       assert(istat==old_index)
       ! Fill edges
       ielem_num_vefs     = mesh%ref_fe_list(ielem_type)%p%get_number_vefs_of_dimension(1)
       ielem_first_vef_id = mesh%ref_fe_list(ielem_type)%p%get_first_vef_id_of_dimension(1)
       vertices_ivef => mesh%ref_fe_list(ielem_type)%p%get_vertices_vef()
       do ivef=1,ielem_num_vefs
          if(mesh%lvefs(mesh%pvefs(ielem)-1+ielem_first_vef_id-1+ivef)==0) then ! Not filled yet
             mesh%nvefs=mesh%nvefs+1                                               ! Count it
             mesh%lvefs(mesh%pvefs(ielem)-1+ielem_first_vef_id-1+ivef)=mesh%nvefs  ! Fill it
             vertices_ivef_iterator = vertices_ivef%create_iterator(ielem_first_vef_id+ivef-1)
             vertex_of_ivef(1) = mesh%lnods(mesh%pnods(ielem)-1+vertices_ivef_iterator%reach_from_current(0))
             vertex_of_ivef(2) = mesh%lnods(mesh%pnods(ielem)-1+vertices_ivef_iterator%reach_from_current(1))
             do jelpo=mesh%pelpo(vertex_of_ivef(1)),mesh%pelpo(vertex_of_ivef(1)+1)-1
                jelem=mesh%lelpo(jelpo)
                if(jelem>ielem) then
                   jelem_num_nodes=mesh%pnods(jelem+1)-mesh%pnods(jelem)
                   call pos_ref_fe%get(key=jelem_num_nodes,val=jelem_type,stat=istat)
                   assert(istat==old_index)
                   jelem_num_vefs     = mesh%ref_fe_list(jelem_type)%p%get_number_vefs_of_dimension(1)
                   jelem_first_vef_id = mesh%ref_fe_list(jelem_type)%p%get_first_vef_id_of_dimension(1)
                   vertices_jvef => mesh%ref_fe_list(jelem_type)%p%get_vertices_vef()
                   vertices_jvef_iterator = vertices_jvef%create_iterator(jelem_first_vef_id+ivef-1)
                   do jvef=1,jelem_num_vefs
                      vertex_of_jvef(1) = mesh%lnods(mesh%pnods(jelem)-1+vertices_jvef_iterator%reach_from_current(0))
                      vertex_of_jvef(2) = mesh%lnods(mesh%pnods(jelem)-1+vertices_jvef_iterator%reach_from_current(1))
                      ! Compare, here we are using that edges have two vertices, hard coded
                      equal = (vertex_of_ivef(1)==vertex_of_jvef(1).and.vertex_of_ivef(2)==vertex_of_jvef(2)).or. &
                           &  (vertex_of_ivef(1)==vertex_of_jvef(2).and.vertex_of_ivef(2)==vertex_of_jvef(1))
                      if(equal) then ! Fill it
                         mesh%lvefs(mesh%pvefs(jelem)-1+jelem_first_vef_id-1+jvef)=mesh%nvefs
                         exit
                      end if
                   end do
                end if
             end do
          end if
       end do
       ! Fill faces (similar code except for the number of vertices of each face that is variable)
       if(mesh%ndime==3) then
          ielem_num_vefs      = mesh%ref_fe_list(ielem_type)%p%get_number_vefs_of_dimension(2)
          ielem_first_vef_id  = mesh%ref_fe_list(ielem_type)%p%get_first_vef_id_of_dimension(2)
          do ivef=1,ielem_num_vefs
             if(mesh%lvefs(mesh%pvefs(ielem)-1+ielem_first_vef_id-1+ivef)==0) then ! Not filled yet
                mesh%nvefs=mesh%nvefs+1                                       ! Count it
                mesh%lvefs(mesh%pvefs(ielem)-1+ielem_first_vef_id-1+ivef)=mesh%nvefs ! Fill it
                ielem_num_vef_verts = mesh%ref_fe_list(ielem_type)%p%get_number_vertices_vef(ielem_first_vef_id+ivef-1)
                vertices_ivef => mesh%ref_fe_list(ielem_type)%p%get_vertices_vef()
                vertex_of_ivef = 0
                vertices_ivef_iterator = vertices_ivef%create_iterator(ielem_first_vef_id+ivef-1)
                do ivert=1,ielem_num_vef_verts
                   vertex_of_ivef(ivert)=mesh%lnods( mesh%pnods(ielem)-1+vertices_ivef_iterator%get_current())
                   call vertices_ivef_iterator%next()
                end do
                do jelpo=mesh%pelpo(vertex_of_ivef(1)),mesh%pelpo(vertex_of_ivef(1)+1)-1
                   jelem=mesh%lelpo(jelpo)
                   if(jelem>ielem) then
                      jelem_num_nodes=mesh%pnods(jelem+1)-mesh%pnods(jelem)
                      call pos_ref_fe%get(key=jelem_num_nodes,val=jelem_type,stat=istat)
                      assert(istat==old_index)
                      jelem_num_vefs      = mesh%ref_fe_list(jelem_type)%p%get_number_vefs_of_dimension(2)
                      jelem_first_vef_id  = mesh%ref_fe_list(jelem_type)%p%get_first_vef_id_of_dimension(2)
                      vertices_jvef => mesh%ref_fe_list(jelem_type)%p%get_vertices_vef()
                      do jvef=1,jelem_num_vefs
                         jelem_num_vef_verts = mesh%ref_fe_list(jelem_type)%p%get_number_vertices_vef(jelem_first_vef_id+jvef-1)
                         vertices_jvef_iterator = vertices_jvef%create_iterator(jelem_first_vef_id+jvef-1)
                         if(jelem_num_vef_verts==ielem_num_vef_verts) then
                            vertex_of_jvef = 0
                            do jvert=1,jelem_num_vef_verts
                               vertex_of_jvef(jvert)=mesh%lnods( mesh%pnods(jelem)-1+vertices_jvef_iterator%get_current())
                               call vertices_jvef_iterator%next()
                            end do
                            count=0
                            do ivert=1,ielem_num_vef_verts
                               do jvert=1,jelem_num_vef_verts
                                  if(vertex_of_ivef(ivert)==vertex_of_jvef(jvert)) then
                                     count=count+1
                                     exit
                                  end if
                               end do
                            end do
                            equal=(count==ielem_num_vef_verts)
                            if(equal) then ! Fill it
                               mesh%lvefs(mesh%pvefs(jelem)-1+jelem_first_vef_id-1+jvef)=mesh%nvefs
                               exit
                            end if
                         end if
                      end do
                   end if
                end do
             end if
          end do
       end if
    end do

    ! Identify boundary faces and assign set and geometry to vefs
    call memalloc(mesh%nvefs, mesh%lvef_geo, __FILE__, __LINE__ )
    call memalloc(mesh%nvefs, mesh%lvef_set, __FILE__, __LINE__ )
    do iboun=1,mesh%bound%get_num_pointers()
       bound_iterator = mesh%bound%create_iterator(iboun)
       nnodb=bound_iterator%get_size()
       if(nnodb==1) then      ! Vertex
          ivert=bound_iterator%reach_from_current(0)
          mesh%lvef_geo(ivert)=mesh%lbgeo(iboun)
          mesh%lvef_set(ivert)=mesh%lbset(iboun)
       else if(nnodb==2) then ! Edge
          vertex_of_ivef(1) = bound_iterator%reach_from_current(0)
          vertex_of_ivef(2) = bound_iterator%reach_from_current(1)
          elems1: do jelpo=mesh%pelpo(vertex_of_ivef(1)),mesh%pelpo(vertex_of_ivef(1)+1)-1
             jelem=mesh%lelpo(jelpo)
             jelem_num_nodes=mesh%pnods(jelem+1)-mesh%pnods(jelem)
             call pos_ref_fe%get(key=jelem_num_nodes,val=jelem_type,stat=istat)
             assert(istat==old_index)
             jelem_num_vefs     = mesh%ref_fe_list(jelem_type)%p%get_number_vefs_of_dimension(1)
             jelem_first_vef_id = mesh%ref_fe_list(jelem_type)%p%get_first_vef_id_of_dimension(1)
             vertices_jvef => mesh%ref_fe_list(jelem_type)%p%get_vertices_vef()
             do jvef=1,jelem_num_vefs
                vertices_jvef_iterator = vertices_jvef%create_iterator(jelem_first_vef_id+jvef-1)
                vertex_of_jvef(1) = mesh%lnods(mesh%pnods(jelem)-1+vertices_jvef_iterator%reach_from_current(0))
                vertex_of_jvef(2) = mesh%lnods(mesh%pnods(jelem)-1+vertices_jvef_iterator%reach_from_current(1))
                ! Compare, here we are using that edges have two vertices, hard coded
                equal = (vertex_of_ivef(1)==vertex_of_jvef(1).and.vertex_of_ivef(2)==vertex_of_jvef(2)).or. &
                     &  (vertex_of_ivef(1)==vertex_of_jvef(2).and.vertex_of_ivef(2)==vertex_of_jvef(1))
                if(equal) then ! Fill it
                   mesh%lvef_geo( mesh%lvefs(mesh%pvefs(jelem)-1+jelem_first_vef_id-1+jvef)  ) = mesh%lbgeo(iboun)
                   mesh%lvef_set( mesh%lvefs(mesh%pvefs(jelem)-1+jelem_first_vef_id-1+jvef)  ) = mesh%lbset(iboun)
                   exit elems1
                end if
             end do
          end do elems1
       else                   ! Face
          vertex_of_ivef = 0
          do ivert=1,nnodb
             vertex_of_ivef(ivert)= bound_iterator%get_current()
             call bound_iterator%next()
          end do
          elems2: do jelpo=mesh%pelpo(vertex_of_ivef(1)),mesh%pelpo(vertex_of_ivef(1)+1)-1
             jelem=mesh%lelpo(jelpo)
             jelem_num_nodes=mesh%pnods(jelem+1)-mesh%pnods(jelem)
             call pos_ref_fe%get(key=jelem_num_nodes,val=jelem_type,stat=istat)
             assert(istat==old_index)
             jelem_num_vefs     = mesh%ref_fe_list(jelem_type)%p%get_number_vefs_of_dimension(2)
             jelem_first_vef_id = mesh%ref_fe_list(jelem_type)%p%get_first_vef_id_of_dimension(2)
             vertices_jvef => mesh%ref_fe_list(jelem_type)%p%get_vertices_vef()
             do jvef=1,jelem_num_vefs
                jelem_num_vef_verts = mesh%ref_fe_list(jelem_type)%p%get_number_vertices_vef(jelem_first_vef_id+jvef-1)
                vertices_jvef_iterator = vertices_jvef%create_iterator(jelem_first_vef_id+jvef-1)
                if(jelem_num_vef_verts==nnodb) then
                   vertex_of_jvef = 0
                   do jvert=1,jelem_num_vef_verts
                      vertex_of_jvef(jvert)=mesh%lnods( mesh%pnods(jelem)-1+vertices_jvef_iterator%get_current())
                      call vertices_jvef_iterator%next()
                   end do
                   count=0
                   do ivert=1,ielem_num_vef_verts
                      do jvert=1,jelem_num_vef_verts
                         if(vertex_of_ivef(ivert)==vertex_of_jvef(jvert)) then
                            count=count+1
                            exit
                         end if
                      end do
                   end do
                   equal=(count==ielem_num_vef_verts)
                   if(equal) then ! Fill it
                      mesh%lvef_geo( mesh%lvefs(mesh%pvefs(jelem)-1+jelem_first_vef_id-1+jvef)  ) = mesh%lbgeo(iboun)
                      mesh%lvef_set( mesh%lvefs(mesh%pvefs(jelem)-1+jelem_first_vef_id-1+jvef)  ) = mesh%lbset(iboun)
                      exit elems2
                   end if
                end if
             end do
          end do elems2
       end if
    end do

  end subroutine mesh_generate_vefs

  
  
  subroutine create_mesh_distribution( femesh, prt_pars, distr, lmesh)
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    implicit none

    ! Parameters
    class(mesh_t)                   , intent(inout)      :: femesh
    type(mesh_distribution_params_t), intent(in)         :: prt_pars
    type(mesh_distribution_t) , allocatable, intent(out) :: distr(:) ! Mesh distribution instances
    type(mesh_t)              , allocatable, intent(out) :: lmesh(:) ! Local mesh instances

    ! Local variables
    type(list_t)                 :: fe_graph    ! Dual graph (to be partitioned)
    integer(ip)   , allocatable  :: ldome(:)    ! Part of each element
    integer(ip)                  :: ipart
    integer                      :: istat

    integer(ip) :: ielem,jelem,iedge,inode,ipoin,jpoin
    real(rp) :: cnorm,vnorm
    integer(ip)   , allocatable  :: weight(:)
    real(rp)      , allocatable  :: coord_i(:),coord_j(:),veloc(:)

    ! Generate dual mesh (i.e., list of elements around points)
    call femesh%to_dual()

    ! Create dual (i.e. list of elements around elements)
    call create_dual_graph(femesh,fe_graph)

    ! ! Create a weigth according to advection
    ! call memalloc (femesh%ndime, coord_i, __FILE__,__LINE__)
    ! call memalloc (femesh%ndime, coord_j, __FILE__,__LINE__)
    ! call memalloc (femesh%ndime, veloc, __FILE__,__LINE__)
    ! call memalloc (fe_graph%p(fe_graph%n+1)-1, weight, __FILE__,__LINE__)
    ! do ielem=1,femesh%nelem
    !    coord_i=0.0_rp
    !    do inode=femesh%pnods(ielem),femesh%pnods(ielem+1)-1
    !       ipoin=femesh%lnods(inode)
    !       coord_i(:)=coord_i(:)+femesh%coord(:,ipoin)
    !    end do
    !    coord_i=coord_i/(femesh%pnods(ielem+1)-femesh%pnods(ielem))
    !    veloc(1) = 1.0_rp !  coord_i(2)
    !    veloc(2) = 0.0_rp !  -coord_i(1)
    !    vnorm=sqrt(veloc(1)*veloc(1)+veloc(2)*veloc(2))
    !    do iedge=fe_graph%p(ielem),fe_graph%p(ielem+1)-1
    !       jelem=fe_graph%l(iedge)
    !       coord_j=0.0_rp
    !       do inode=femesh%pnods(jelem),femesh%pnods(jelem+1)-1
    !          jpoin=femesh%lnods(inode)
    !          coord_j(:)=coord_j(:)+femesh%coord(:,jpoin)
    !       end do
    !       coord_j=coord_j/(femesh%pnods(jelem+1)-femesh%pnods(jelem))
    !       coord_j(1)=coord_j(1)-coord_i(1)
    !       coord_j(2)=coord_j(2)-coord_i(2)
    !       cnorm=sqrt(coord_j(1)*coord_j(1)+coord_j(2)*coord_j(2))
    !       if(vnorm*cnorm>1.0e-8) then
    !          weight(iedge) = int(100.0_rp * ((veloc(1)*coord_j(1)+veloc(2)*coord_j(2)) &
    !               &                       / (vnorm*cnorm))**10 + 0)
    !       else
    !          weight(iedge) = 0
    !       end if
    !    end do
    ! end do

    ! Partition dual graph to assign a domain to each element (in ldome)
    call memalloc (femesh%nelem, ldome, __FILE__,__LINE__)   
    ! ! write(*,*) weight
    ! call graph_pt_renumbering(prt_pars,fe_graph,ldome,weight)
    ! call memfree ( coord_i, __FILE__,__LINE__)
    ! call memfree ( coord_j, __FILE__,__LINE__)
    ! call memfree ( veloc, __FILE__,__LINE__)
    ! call memfree ( weight, __FILE__,__LINE__)
    call graph_pt_renumbering(prt_pars,fe_graph,ldome)

    ! Now free fe_graph, not needed anymore?
    allocate(distr(prt_pars%nparts), stat=istat)
    check(istat==0)
    allocate(lmesh(prt_pars%nparts), stat=istat)
    check(istat==0) 

    do ipart=1, prt_pars%nparts
       distr(ipart)%ipart  = ipart
       distr(ipart)%nparts = prt_pars%nparts
    end do
    call build_maps(prt_pars%nparts, ldome, femesh, distr)

    ! Build local meshes and their duals and generate partition adjacency
    do ipart=1,prt_pars%nparts

       ! Generate Local mesh
       call mesh_g2l(distr(ipart)%num_local_vertices,  &
                     distr(ipart)%l2g_vertices,        &
                     distr(ipart)%num_local_cells,     &
                     distr(ipart)%l2g_cells,           &
                     femesh,                           &
                     lmesh(ipart))

       call build_adjacency_new (femesh, ldome,             &
            &                    ipart,                     &
            &                    lmesh(ipart),              &
            &                    distr(ipart)%l2g_vertices, &
            &                    distr(ipart)%l2g_cells,    &
            &                    distr(ipart)%nebou,        &
            &                    distr(ipart)%nnbou,        &
            &                    distr(ipart)%lebou,        &
            &                    distr(ipart)%lnbou,        &
            &                    distr(ipart)%pextn,        &
            &                    distr(ipart)%lextn,        &
            &                    distr(ipart)%lextp )
    end do
    call fe_graph%free()
    call memfree(ldome,__FILE__,__LINE__)
  end subroutine create_mesh_distribution

  !================================================================================================
  subroutine create_dual_graph(mesh,graph)
    ! Parameters
    type(mesh_t) , intent(in)  :: mesh
    type(list_t),  intent(out) :: graph
    ! Locals
    integer(ip), allocatable :: lelem(:)
    integer(ip), allocatable :: keadj(:)

    call memalloc(           mesh%nelem,lelem,__FILE__,__LINE__)
    call memalloc(mesh%nelpo*mesh%nnode,keadj,__FILE__,__LINE__)
    lelem=0
    call graph%create(mesh%nelem)

    call count_elemental_graph(mesh%ndime,mesh%npoin,mesh%nelem, &
    !call count_elemental_graph(1,mesh%npoin,mesh%nelem, &
         &                     mesh%pnods,mesh%lnods,mesh%nnode,mesh%nelpo, &
         &                     mesh%pelpo,mesh%lelpo,lelem,keadj,graph)

    call graph%allocate_list_from_pointer()

    call list_elemental_graph(mesh%ndime,mesh%npoin,mesh%nelem, &
    !call list_elemental_graph(1,mesh%npoin,mesh%nelem, &
         &                    mesh%pnods,mesh%lnods,mesh%nnode,mesh%nelpo, &
         &                    mesh%pelpo,mesh%lelpo,lelem,keadj,graph)

    call memfree(lelem,__FILE__,__LINE__)
    call memfree(keadj,__FILE__,__LINE__)

  end subroutine create_dual_graph

  !-----------------------------------------------------------------------
  subroutine count_elemental_graph(ncomm,npoin,nelem,pnods,lnods,nnode, &
       &                           nelpo,pelpo,lelpo,lelem,keadj,graph)
    implicit none
    integer(ip),  intent(in)    :: ncomm,npoin,nelem
    integer(ip),  intent(in)    :: pnods(nelem+1),lnods(pnods(nelem+1))
    integer(ip),  intent(in)    :: nnode             ! Number of nodes of each element (max.)
    integer(ip),  intent(in)    :: nelpo             ! Number of elements around points (max.)
    integer(ip),  intent(in)    :: pelpo(npoin+1)    ! Number of elements around points
    integer(ip),  intent(in)    :: lelpo(pelpo(npoin+1))    ! List of elements around points
    type(list_t), intent(inout) :: graph                    ! Number of edges on each element (list_t)
    integer(ip),  intent(out)   :: lelem(nelem)             ! Auxiliar array
    integer(ip),  intent(out)   :: keadj(nelpo*nnode)       ! Auxiliar array
    integer(ip)                 :: ielem,jelem,inode,knode  ! Indices
    integer(ip)                 :: ipoin,ielpo,index        ! Indices
    integer(ip)                 :: neadj,ielel,jelel,nelel  ! Indices

    lelem=0
    neadj=1
    knode=nnode
    call graph%create(nelem)
    do ielem=1,nelem
       call graph%sum_to_pointer_index(ielem-1, neadj)
       ! Loop over nodes and their surrounding elements and count
       ! how many times they are repeated as neighbors of ielem
       nelel=0
       keadj=0
       index=pnods(ielem)-1
       knode=pnods(ielem+1)-pnods(ielem)
       do inode=1,knode
          ipoin=lnods(index+inode)
          do ielpo=pelpo(ipoin),pelpo(ipoin+1)-1
             jelem=lelpo(ielpo)
             if(lelem(jelem)==0) then
                nelel=nelel+1
                keadj(nelel)=jelem
             end if
             lelem(jelem)=lelem(jelem)+1
          end do
       end do

       ! Now we loop over the elements around ielem and define neighbors
       ! as those sharing ncomm nodes. The smaller this number is the
       ! higher the connectivity of elemental graph is. Note that prisms,
       ! for example, could share 3 or 4 nodes depending on the face, so
       ! ndime (the number of space dimensions) is a good choice.
       jelel=0
       do ielel=1,nelel
          jelem=keadj(ielel)
          if(lelem(jelem)>=ncomm) jelel=jelel+1
       end do
       neadj=neadj+jelel

       ! Reset lelem
       do ielel=1,nelel
          jelem=keadj(ielel)
          lelem(jelem)=0
       end do
    end do

    call graph%sum_to_pointer_index(nelem, neadj)
    call graph%calculate_header()

  end subroutine count_elemental_graph
  !-----------------------------------------------------------------------
  subroutine list_elemental_graph(ncomm,npoin,nelem,pnods,lnods,nnode, &
       &                          nelpo,pelpo,lelpo,lelem,keadj,graph)
    implicit none
    integer(ip),  intent(in)    :: ncomm,npoin,nelem
    integer(ip),  intent(in)    :: pnods(nelem+1),lnods(pnods(nelem+1))
    integer(ip),  intent(in)    :: nnode             ! Number of nodes of each element (max.)
    integer(ip),  intent(in)    :: nelpo             ! Number of elements around points (max.)
    integer(ip),  intent(in)    :: pelpo(npoin+1)    ! Number of elements around points
    integer(ip),  intent(in)    :: lelpo(pelpo(npoin+1))    ! List of elements around points
    type(list_t), intent(inout) :: graph                    ! List of edges on each element (list_t)
    integer(ip),  intent(out)   :: lelem(nelem)             ! Auxiliar array
    integer(ip),  intent(out)   :: keadj(nelpo*nnode)       ! Auxiliar array
    integer(ip)                 :: ielem,jelem,inode,knode  ! Indices
    integer(ip)                 :: ipoin,ielpo,index        ! Indices
    integer(ip)                 :: neadj,ielel,jelel,nelel  ! Indices
    type(list_iterator_t)       :: graph_iterator

    lelem=0
    knode=nnode
    do ielem=1,nelem
       ! Loop over nodes and their surrounding elements and count
       ! how many times they are repeated as neighbors of ielem
       nelel=0
       keadj=0
       index=pnods(ielem)-1
       knode=pnods(ielem+1)-pnods(ielem)
       do inode=1,knode
          ipoin=lnods(index+inode)
          do ielpo=pelpo(ipoin),pelpo(ipoin+1)-1
             jelem=lelpo(ielpo)
             if(lelem(jelem)==0) then
                nelel=nelel+1
                keadj(nelel)=jelem
             end if
             lelem(jelem)=lelem(jelem)+1
          end do
       end do

       ! Now we loop over the elements around ielem and define neighbors
       call graph%allocate_list_from_pointer()
       graph_iterator = graph%create_iterator(ielem)
       do ielel=1,nelel
          jelem=keadj(ielel)
          if(lelem(jelem)>=ncomm) then
             call graph_iterator%set_current(jelem)
             call graph_iterator%next()
          end if
       end do

       ! Reset lelem
       do ielel=1,nelel
          jelem=keadj(ielel)
          lelem(jelem)=0
       end do
    end do

  end subroutine list_elemental_graph

  !================================================================================================
  subroutine build_adjacency_new ( gmesh, ldome, my_part, lmesh, l2gn, l2ge, &
       &                           nebou, nnbou, lebou, lnbou, pextn, lextn, lextp)
    implicit none
    integer(ip)   , intent(in)  :: my_part
    type(mesh_t)  , intent(in)  :: gmesh,lmesh
    !type(mesh_t)  , intent(in)  :: dual_lmesh
    integer(ip)   , intent(in)  :: ldome(gmesh%nelem)
    integer(igp)  , intent(in)  :: l2gn(lmesh%npoin)
    integer(igp)  , intent(in)  :: l2ge(lmesh%nelem)
    !integer(ip)   , intent(in)  :: dual_parts( dual_lmesh%pnods(dual_lmesh%nelem+1)-1)
    integer(ip)   , intent(out) :: nebou
    integer(ip)   , intent(out) :: nnbou
    integer(ip)   , allocatable, intent(out) ::  lebou(:)    ! List of boundary elements
    integer(ip)   , allocatable, intent(out) ::  lnbou(:)    ! List of boundary nodes
    integer(ip)   , allocatable, intent(out) ::  pextn(:)    ! Pointers to the lextn
    integer(igp)  , allocatable, intent(out) ::  lextn(:)    ! List of (GID of) external neighbors
    integer(ip)   , allocatable, intent(out) ::  lextp(:)    ! List of parts of external neighbors

    integer(ip) :: lelem, ielem, jelem, pelem, pnode, inode1, inode2, ipoin, lpoin, jpart, iebou, istat, touch
    integer(ip) :: nextn, nexte, nepos
    integer(ip), allocatable :: local_visited(:)
    type(hash_table_ip_ip_t)   :: external_visited

    if(my_part==0) then
       write(*,*)  'Parts:'
       do ielem=1,gmesh%nelem
          write(*,*)  ielem, ldome(ielem)
       end do
       write(*,*)  'Global mesh:',gmesh%npoin,gmesh%nelem
       do ielem=1,gmesh%nelem
          write(*,*)  ielem, gmesh%lnods(gmesh%pnods(ielem):gmesh%pnods(ielem+1)-1)
       end do
       write(*,*)  'Global dual mesh:',gmesh%nelpo
       do ipoin=1,gmesh%npoin
          write(*,*)  ipoin, gmesh%lelpo(gmesh%pelpo(ipoin):gmesh%pelpo(ipoin+1)-1)
       end do
       write(*,*)  'Local mesh:',lmesh%npoin,lmesh%nelem
       do lelem=1,lmesh%nelem
          write(*,*)  lelem, l2ge(lelem),lmesh%lnods(lmesh%pnods(lelem):lmesh%pnods(lelem+1)-1)
       end do
       write(*,*)  'Local2Global (nodes)'
       do lpoin=1,lmesh%npoin
          write(*,*)  lpoin, l2gn(lpoin)
       end do
    end if

    ! Count boundary nodes
    nnbou = 0 
    do lpoin=1, lmesh%npoin
       ipoin = l2gn(lpoin)
       do pelem = gmesh%pelpo(ipoin), gmesh%pelpo(ipoin+1) - 1
          ielem = gmesh%lelpo(pelem)
          jpart = ldome(ielem)
          if ( jpart /= my_part ) then 
             nnbou = nnbou +1
             exit
          end if
       end do
    end do

    ! List boundary nodes
    call memalloc ( nnbou, lnbou, __FILE__, __LINE__ ) 
    nnbou = 0
    do lpoin=1, lmesh%npoin
       ipoin = l2gn(lpoin)
       do pelem = gmesh%pelpo(ipoin), gmesh%pelpo(ipoin+1) - 1
          ielem = gmesh%lelpo(pelem)
          jpart = ldome(ielem)
          if ( jpart /= my_part ) then 
             lnbou(nnbou+1) = ipoin
             nnbou = nnbou +1
             exit
          end if
       end do
    end do

    ! As the dual mesh is given with global IDs we need a hash table to do the touch.
    call memalloc(lmesh%nelem, local_visited,__FILE__,__LINE__)
    local_visited = 0
    call external_visited%init(20)

    ! 1) Count boundary elements and external edges
    touch = 1
    nebou = 0 ! number of boundary elements
    nextn = 0 ! number of external edges
    do lelem = 1, lmesh%nelem
       nexte = 0   ! number of external neighbours of this element
       ielem = l2ge(lelem)
       inode1 = gmesh%pnods(ielem)
       inode2 = gmesh%pnods(ielem+1)-1
       do pnode = inode1, inode2
          ipoin = gmesh%lnods(pnode)
          do pelem = gmesh%pelpo(ipoin), gmesh%pelpo(ipoin+1) - 1
             jelem = gmesh%lelpo(pelem)
             if(jelem/=ielem) then
                jpart = ldome(jelem)
                if(jpart/=my_part) then                                   ! This is an external element
                   if(local_visited(lelem) == 0 ) nebou = nebou +1        ! Count it
                   !call external_visited%put(key=jelem,val=1, stat=istat) ! Touch jelem as external neighbor of lelem.
                   call external_visited%put(key=jelem,val=touch, stat=istat) ! Touch jelem as external neighbor of lelem.
                   if(istat==now_stored) nexte = nexte + 1                ! Count external neighbours of lelem
                   local_visited(lelem) = nexte                           ! Touch lelem also storing the number
                end if                                                    ! of external neighbours it has
             end if
          end do
       end do
       nextn = nextn + nexte
       ! Clean hash table
       if(local_visited(lelem) /= 0 ) then 
          do pnode = inode1, inode2
             ipoin = gmesh%lnods(pnode)
             do pelem = gmesh%pelpo(ipoin), gmesh%pelpo(ipoin+1) - 1
                jelem = gmesh%lelpo(pelem)
                if(jelem/=ielem) then
                   jpart = ldome(jelem)
                   if(jpart/=my_part) then
                      call external_visited%del(key=jelem, stat=istat)
                   end if
                end if
             end do
          end do
       end if
       call external_visited%print
    end do

    if(my_part==0) then
       write(*,*)  'Visited (boundary) elements:'
       do lelem=1,lmesh%nelem
          write(*,*)  local_visited(lelem)
       end do
    end if

    ! 2) Allocate arrays and store list and pointers to externals
    call memalloc(nebou  , lebou,__FILE__,__LINE__)
    call memalloc(nebou+1, pextn,__FILE__,__LINE__)
    call memalloc(nextn  , lextn,__FILE__,__LINE__)
    call memalloc(nextn  , lextp,__FILE__,__LINE__)

    iebou = 0
    pextn(1) = 1
    do lelem = 1, lmesh%nelem
       if(local_visited(lelem) /= 0 ) then
          iebou = iebou +1
          lebou(iebou) = lelem
          pextn(iebou+1) = local_visited(lelem) + pextn(iebou)
       end if
    end do

    if(my_part==0) then
       write(*,*)  'Boundary elements:'
       do iebou=1,nebou
          write(*,*)  lebou(iebou)
       end do
    end if

    ! 3) Store boundary elements and external edges
    !do lelem = 1, lmesh%nelem
    do iebou = 1, nebou
       lelem = lebou(iebou)
       ielem = l2ge(lelem)
       nexte = 0   ! number of external neighbours of this element
       inode1 = gmesh%pnods(ielem)
       inode2 = gmesh%pnods(ielem+1)-1
       do pnode = inode1, inode2
          ipoin = gmesh%lnods(pnode)
          do pelem = gmesh%pelpo(ipoin), gmesh%pelpo(ipoin+1) - 1
             jelem = gmesh%lelpo(pelem)
             if(jelem/=ielem) then
                jpart = ldome(jelem)
                if(jpart/=my_part) then                                   ! This is an external element
                   call external_visited%put(key=jelem,val=touch, stat=istat) ! Touch jelem as external neighbor of lelem.
                   if(istat==now_stored) then
                      lextn(pextn(iebou)+nexte) = jelem
                      lextp(pextn(iebou)+nexte) = jpart
                      nexte = nexte + 1
                   end if
                end if
             end if
          end do
       end do
       ! Clean hash table
       do pnode = inode1, inode2
          ipoin = gmesh%lnods(pnode)
          do pelem = gmesh%pelpo(ipoin), gmesh%pelpo(ipoin+1) - 1
             jelem = gmesh%lelpo(pelem)
             if(jelem/=ielem) then
                jpart = ldome(jelem)
                if(jpart/=my_part) then                                   ! This is an external element
                   call external_visited%del(key=jelem, stat=istat)
                end if
             end if
          end do
       end do
    end do

    call external_visited%free
    call memfree(local_visited,__FILE__,__LINE__)

  end subroutine build_adjacency_new

  !================================================================================================
  subroutine build_maps( nparts, ldome, femesh, distr )
    ! This routine builds (node and element) partition maps without using the objects
    ! and (unlike parts_sizes, parts_maps, etc.) does not generate a new global numbering.
    implicit none
    integer(ip)                , intent(in)    :: nparts
    type(mesh_t)             , intent(in)    :: femesh
    integer(ip)                , intent(in)    :: ldome(femesh%nelem)
    type(mesh_distribution_t), intent(inout) :: distr(nparts)

    integer(ip)   , allocatable  :: nedom(:) ! Number of points per part (here is not header!)
    integer(ip)   , allocatable  :: npdom(:) ! Number of elements per part (here is not header!)
    integer(ip)   , allocatable  :: work1(:)
    integer(ip)   , allocatable  :: work2(:)
    integer(ip) :: ielem, ipart, inode, iboun

    ! Number of elements of each part and global to local element map (is one to one)
    call memalloc (nparts, nedom,__FILE__,__LINE__)
    nedom=0
    do ielem=1,femesh%nelem
       ipart = ldome(ielem)
       nedom(ipart)=nedom(ipart)+1
    end do
    ! Allocate local to global maps
    do ipart=1,nparts
       distr(ipart)%num_local_cells  = nedom(ipart)
       distr(ipart)%num_global_cells = int(femesh%nelem,igp)
       call memalloc(distr(ipart)%num_local_cells, distr(ipart)%l2g_cells, __FILE__, __LINE__)
    end do
    nedom = 0
    do ielem=1,femesh%nelem
       ipart = ldome(ielem)
       nedom(ipart)=nedom(ipart)+1
       distr(ipart)%l2g_cells(nedom(ipart)) = ielem
    end do

    call memfree ( nedom,__FILE__,__LINE__)

    ! Number of nodes of each part and global to local node map (is NOT one to one)
    call memalloc ( nparts, npdom,__FILE__,__LINE__)
    call memalloc ( femesh%npoin, work1,__FILE__,__LINE__)
    call memalloc ( femesh%npoin, work2,__FILE__,__LINE__)
    npdom=0
    do ipart = 1, nparts
       work1 = 0
       work2 = 0
       do ielem=1,femesh%nelem
          if(ldome(ielem)==ipart) then
             do inode = femesh%pnods(ielem), femesh%pnods(ielem+1) - 1 
                if(work1(femesh%lnods(inode)) == 0 ) then
                   npdom(ipart) = npdom(ipart)+1
                   work1(femesh%lnods(inode)) = 1
                   work2(npdom(ipart)) = femesh%lnods(inode)
                end if
             end do
          end if
       end do
       distr(ipart)%num_local_vertices  = npdom(ipart)
       distr(ipart)%num_global_vertices = int(femesh%npoin,igp)
       call memalloc(distr(ipart)%num_local_vertices, distr(ipart)%l2g_vertices, __FILE__, __LINE__)
       distr(ipart)%l2g_vertices = work2(1:npdom(ipart))
    end do
    call memfree ( work1,__FILE__,__LINE__)
    call memfree ( work2,__FILE__,__LINE__)
    call memfree ( npdom,__FILE__,__LINE__)
  end subroutine build_maps

  ! Inspired on http://en.wikipedia.org/wiki/Breadth-first_search.
  ! Given a mesh (m) and its dual graph (g), it computes the list 
  ! of nodes (lconn) of each connected component in m. Can be very
  ! useful as a tool to determine whether the mesh partitioning process
  ! leads to disconnected subdomains or not.
  subroutine mesh_graph_compute_connected_components (m, g, lconn)
    implicit none

    ! Parameters
    type(mesh_t) , intent(in)   :: m   
    type(list_t),  intent(in)   :: g
    type(list_t),  intent(out)  :: lconn

    ! Locals
    integer(ip), allocatable :: auxv(:), auxe(:), e(:)
    integer(ip), allocatable :: emarked(:), vmarked(:)
    integer(ip), allocatable :: q(:)
    integer(ip)              :: head, tail, i, esize, vsize, current, & 
         j, l, k, inods1d, inods2d, p_ipoin, ipoin, graph_num_rows, lconnn
    type(list_iterator_t)    :: graph_column_iterator
    type(list_iterator_t)    :: lconn_iterator

    graph_num_rows = g%get_num_pointers()
    call memalloc ( graph_num_rows   , auxe     , __FILE__,__LINE__)
    call memalloc ( graph_num_rows   , auxv     , __FILE__,__LINE__)
    call memalloc ( graph_num_rows   , q        , __FILE__,__LINE__)
    call memalloc ( graph_num_rows   , emarked  , __FILE__,__LINE__)
    call memalloc ( m%npoin          , vmarked  , __FILE__,__LINE__)
    call memalloc ( graph_num_rows   ,  e       , __FILE__,__LINE__)

    lconnn  = 0
    emarked  = 0
    current  = 1 

    do i=1, graph_num_rows
       if (emarked(i) == 0) then
          ! New connected component
          lconnn = lconnn +1
          esize   = 0
          vsize   = 0
          vmarked = 0 
!!$1  procedure BFS(G,v):
!!$2      create a queue Q
          head=1
          tail=1
!!$3      enqueue v onto Q
          q(tail)=i
          tail=tail+1
!!$4      mark v
          emarked(i)=1
          e(current)=i
          esize  = esize + 1
          current = current + 1  

!!$5      while Q is not empty:
          do while (head/=tail)
!!$6         t ← Q.dequeue()
             j=q(head)
             head = head + 1

             ! Traverse the nodes of the element number j
             inods1d = m%pnods(j)
             inods2d = m%pnods(j+1)-1

             do p_ipoin = inods1d, inods2d
                ipoin = m%lnods(p_ipoin)
                if (vmarked(ipoin)==0) then
                   vmarked(ipoin)=1
                   vsize = vsize+1
                end if
             end do

!!$9         for all edges e in G.adjacentEdges(t) do
             graph_column_iterator = g%create_iterator(j)
             do while(.not. graph_column_iterator%is_upper_bound())
!!$12           u ← G.adjacentVertex(t,e)
                l=graph_column_iterator%get_current()
!!$13           if u is not emarked:
                if (emarked(l)==0) then
!!$14              mark u
                   emarked(l)=1
                   e(current)=l
                   esize  = esize + 1
                   current = current + 1  

!!$15              enqueue u onto Q
                   q(tail)=l
                   tail=tail+1
                end if
                call graph_column_iterator%next()
             end do
          end do
          auxe(lconnn) = esize
          auxv(lconnn) = vsize
       end if
    end do

    call lconn%create(lconnn)

    do i=1, lconnn
       call lconn%sum_to_pointer_index(i, auxv(i))
    end do

    call memfree( auxv   ,__FILE__,__LINE__)
    call memfree( q      ,__FILE__,__LINE__)
    call memfree( emarked,__FILE__,__LINE__)

    call lconn%calculate_header()
    call lconn%allocate_list_from_pointer()

    current=1
    lconn_iterator = lconn%create_iterator()
    do i=1, lconn_iterator%get_size()
       vmarked = 0
       ! Traverse elements of current connected component  
       do current=current,current+auxe(i)-1
          j=e(current)

          ! Traverse the nodes of the element number j
          inods1d = m%pnods(j)
          inods2d = m%pnods(j+1)-1

          do p_ipoin = inods1d, inods2d
             ipoin = m%lnods(p_ipoin)
             if (vmarked(ipoin)==0) then
                vmarked(ipoin)=1
                call lconn_iterator%set_current(ipoin)
                call lconn_iterator%next()
             end if
          end do

       end do
    end do

    ! write(*,*) 'ZZ', g%nv, m%npoin, l
    ! write(*,*) 'XX', lconn%n
    ! write(*,*) 'YY', lconn%p
    ! write(*,*) 'PP', lconn%l

    call memfree( auxe,__FILE__,__LINE__)
    call memfree( e,__FILE__,__LINE__)
    call memfree( vmarked,__FILE__,__LINE__)

  end subroutine mesh_graph_compute_connected_components

  !=================================================================================================
  subroutine graph_nd_renumbering(prt_parts, gp, iperm, lperm)
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    implicit none
    type(mesh_distribution_params_t), intent(in)         :: prt_parts
    type(list_t)                    , target, intent(inout) :: gp
    integer(ip)                     , target, intent(out):: iperm(gp%get_size())
    integer(ip)                     , target, intent(out):: lperm(gp%get_size())
    
    if ( gp%get_num_pointers() == 1 ) then
       lperm(1) = 1
       iperm(1) = 1
    else
#ifdef ENABLE_METIS
       ierr = metis_setdefaultoptions(c_loc(options))
       assert(ierr == METIS_OK) 
       
       options(METIS_OPTION_NUMBERING) = 1
       options(METIS_OPTION_DBGLVL)    = prt_parts%metis_option_debug
       
       ierr = metis_nodend ( gp%get_num_pointers_c_loc() ,gp%get_pointers_c_loc() , gp%get_list_c_loc(), &
            &                C_NULL_PTR, c_loc(options), c_loc(iperm),c_loc(lperm))
       
       assert(ierr == METIS_OK)
#else
       call enable_metis_error_message
#endif
    end if
  end subroutine graph_nd_renumbering

  !=================================================================================================
  subroutine graph_pt_renumbering(prt_parts,gp,ldomn,weight)
    !-----------------------------------------------------------------------
    ! This routine computes a nparts-way-partitioning of the input graph gp
    !-----------------------------------------------------------------------
    implicit none
    type(mesh_distribution_params_t), target, intent(in)    :: prt_parts
    type(list_t)                    , target, intent(inout) :: gp
    integer(ip)                     , target, intent(out)   :: ldomn(gp%get_num_pointers())
    integer(ip)                     , target, optional, intent(in)  :: weight(gp%get_size())

    ! Local variables 
    integer(ip), target      :: kedge
    integer(ip)              :: idumm,iv
    integer(ip), allocatable :: lwork(:)
    integer(ip)              :: i, j, m, k, ipart
    integer(ip), allocatable :: iperm(:)
   
#ifdef ENABLE_METIS
    ierr = metis_setdefaultoptions(c_loc(options))
    assert(ierr == METIS_OK) 

!!$      From METIS 5.0 manual:
!!$
!!$      The following options are valid for METIS PartGraphRecursive:
!!$      
!!$      METIS_OPTION_CTYPE, METIS_OPTION_IPTYPE, METIS_OPTION_RTYPE,
!!$      METIS_OPTION_NO2HOP, METIS_OPTION_NCUTS, METIS_OPTION_NITER,
!!$      METIS_OPTION_SEED, METIS_OPTION_UFACTOR, METIS_OPTION_NUMBERING,
!!$      METIS_OPTION_DBGLVL
!!$     
!!$      The following options are valid for METIS PartGraphKway:
!!$ 
!!$      METIS_OPTION_OBJTYPE, METIS_OPTION_CTYPE, METIS_OPTION_IPTYPE,
!!$      METIS_OPTION_RTYPE, METIS_OPTION_NO2HOP, METIS_OPTION_NCUTS,
!!$      METIS_OPTION_NITER, METIS_OPTION_UFACTOR, METIS_OPTION_MINCONN,
!!$      METIS_OPTION_CONTIG, METIS_OPTION_SEED, METIS_OPTION_NUMBERING,
!!$      METIS_OPTION_DBGLVL

    if ( prt_parts%strat == part_kway ) then
       options(METIS_OPTION_NUMBERING) = 1
       options(METIS_OPTION_DBGLVL)    = prt_parts%metis_option_debug
       
       ! Enforce contiguous partititions
       options(METIS_OPTION_CONTIG)    = prt_parts%metis_option_contig
       
       ! Explicitly minimize the maximum degree of the subdomain graph
       options(METIS_OPTION_MINCONN)   = prt_parts%metis_option_minconn
       options(METIS_OPTION_UFACTOR)   = prt_parts%metis_option_ufactor

       ! Select random (default) or sorted heavy edge matching
       options(METIS_OPTION_CTYPE)     = prt_parts%metis_option_ctype
       options(METIS_OPTION_IPTYPE)     = prt_parts%metis_option_iptype
       
       ncon = 1 
       
       write(*,*) 'k_way',present(weight)
       if(present(weight)) then
          write(*,*) 'calling metis',options(METIS_OPTION_CTYPE)
          options(METIS_OPTION_NITER) = 100

          ierr = metis_partgraphkway( gp%get_num_pointers_c_loc(), c_loc(ncon), gp%get_pointers_c_loc(), gp%get_list_c_loc() , & 
!                                    vw             vsize       adjw
                                      C_NULL_PTR  , C_NULL_PTR , c_loc(weight) , c_loc(prt_parts%nparts), &
                                      C_NULL_PTR  , C_NULL_PTR , c_loc(options), c_loc(kedge), c_loc(ldomn) )
       else
          ierr = metis_partgraphkway( gp%get_num_pointers_c_loc(), c_loc(ncon), gp%get_pointers_c_loc(), gp%get_list_c_loc() , & 
                                      C_NULL_PTR  , C_NULL_PTR , C_NULL_PTR    , c_loc(prt_parts%nparts), &
                                      C_NULL_PTR  , C_NULL_PTR , c_loc(options), c_loc(kedge), c_loc(ldomn) )
       end if

       assert(ierr == METIS_OK) 
       
    else if ( prt_parts%strat == part_recursive ) then
       write(*,*) 'part_recursive',present(weight)

       options(METIS_OPTION_NUMBERING) = 1
       options(METIS_OPTION_DBGLVL)    = prt_parts%metis_option_debug
       options(METIS_OPTION_UFACTOR)   = prt_parts%metis_option_ufactor

       ncon = 1 
       ierr = metis_partgraphrecursive( gp%get_num_pointers_c_loc(), c_loc(ncon), gp%get_pointers_c_loc(), gp%get_list_c_loc() , & 
                                        C_NULL_PTR  , C_NULL_PTR , C_NULL_PTR    , c_loc(prt_parts%nparts), &
                                        C_NULL_PTR  , C_NULL_PTR , c_loc(options), c_loc(kedge), c_loc(ldomn) )
    end if    
#else
    call enable_metis_error_message
#endif

    if ( prt_parts%strat == part_strip ) then
       j = gp%get_num_pointers()
       m = 0
       do ipart=1,prt_parts%nparts
          k = j / (prt_parts%nparts-ipart+1)
          do i = 1, k
             ldomn(m+i) = ipart
          end do
          m = m + k
          j = j - k
       end do
    else if ( prt_parts%strat == part_rcm_strip ) then
       call memalloc ( gp%get_num_pointers(), iperm, __FILE__,__LINE__ )
       call genrcm ( gp, iperm )
       j = gp%get_num_pointers()
       m = 0
       do ipart=1,prt_parts%nparts
          k = j / (prt_parts%nparts-ipart+1)
          do i = 1, k
             ldomn(iperm(m+i)) = ipart
          end do
          m = m + k
          j = j - k
       end do
       call memfree ( iperm,__FILE__,__LINE__)
    end if

  end subroutine graph_pt_renumbering


  !================================================================================================
  subroutine mesh_g2l(num_local_vertices, l2g_vertices, num_local_cells, l2g_cells, gmesh, lmesh)
    implicit none
    integer(ip),     intent(in)    :: num_local_vertices
    integer(igp),    intent(in)    :: l2g_vertices(num_local_vertices)
    integer(ip),     intent(in)    :: num_local_cells
    integer(igp),    intent(in)    :: l2g_cells(num_local_cells)
    type(mesh_t)   , intent(in)    :: gmesh
    type(mesh_t)   , intent(inout) :: lmesh
    type(hash_table_igp_ip_t)      :: ws_inmap
    type(hash_table_igp_ip_t)      :: el_inmap
    integer(ip)    , allocatable   :: node_list(:)
    integer(ip)                    :: aux, ipoin,inode,inodb,knode,knodb,lnodb_size,istat
    integer(ip)                    :: ielem_lmesh,ielem_gmesh,iboun_lmesh,iboun_gmesh
    integer(ip)                    :: p_ielem_gmesh,p_ipoin_lmesh,p_ipoin_gmesh
    type(list_iterator_t)          :: bound_iterator
    logical :: count_it


    lmesh%order=gmesh%order
    lmesh%nelty=gmesh%nelty
    lmesh%ndime=gmesh%ndime
    lmesh%npoin=num_local_vertices
    lmesh%nelem=num_local_cells

    call ws_inmap%init(max(int(num_local_vertices*0.25,ip),10))
    do ipoin=1,num_local_vertices
       ! aux is used to avoid compiler warning related to val being an intent(inout) argument
       aux = ipoin
       call ws_inmap%put(key=l2g_vertices(ipoin),val=aux,stat=istat) 
    end do

    call el_inmap%init(max(int(num_local_cells*0.25,ip),10))
    do ipoin=1,num_local_cells
       ! aux is used to avoid compiler warning related to val being an intent(inout) argument
       aux = ipoin
       call el_inmap%put(key=l2g_cells(ipoin),val=aux,stat=istat) 
    end do

    ! Elements
    call memalloc(lmesh%nelem+1, lmesh%pnods, __FILE__,__LINE__)
    call memalloc(lmesh%nelem  , lmesh%legeo, __FILE__,__LINE__)
    call memalloc(lmesh%nelem  , lmesh%leset, __FILE__,__LINE__)
    lmesh%nnode=0
    lmesh%pnods=0
    lmesh%pnods(1)=1
    do ielem_lmesh=1,lmesh%nelem
       ielem_gmesh = l2g_cells(ielem_lmesh)
       knode = gmesh%pnods(ielem_gmesh+1)-gmesh%pnods(ielem_gmesh)
       lmesh%pnods(ielem_lmesh+1)=lmesh%pnods(ielem_lmesh)+knode
       lmesh%nnode=max(lmesh%nnode,knode)
       lmesh%legeo(ielem_lmesh)=gmesh%legeo(ielem_gmesh)
       lmesh%leset(ielem_lmesh)=gmesh%leset(ielem_gmesh)
    end do
    call memalloc (lmesh%pnods(lmesh%nelem+1), lmesh%lnods, __FILE__,__LINE__)
    do ielem_lmesh=1,lmesh%nelem
       ielem_gmesh = l2g_cells(ielem_lmesh)
       p_ipoin_gmesh = gmesh%pnods(ielem_gmesh)-1
       p_ipoin_lmesh = lmesh%pnods(ielem_lmesh)-1
       knode = gmesh%pnods(ielem_gmesh+1)-gmesh%pnods(ielem_gmesh)
       do inode=1,knode
          call ws_inmap%get(key=int(gmesh%lnods(p_ipoin_gmesh+inode),igp),val=lmesh%lnods(p_ipoin_lmesh+inode),stat=istat) 
       end do
    end do

    ! Boundary elements
    iboun_lmesh=0
    lmesh%nnodb=0
    lnodb_size=0
    do iboun_gmesh=1,gmesh%bound%get_num_pointers()
       bound_iterator = gmesh%bound%create_iterator(iboun_gmesh)
       knodb = bound_iterator%get_size()
       count_it=.true.
       do while(.not. bound_iterator%is_upper_bound())
          call ws_inmap%get(key=int(bound_iterator%get_current(),igp),val=knode,stat=istat)
          call bound_iterator%next()
          if(istat==key_not_found) then
             count_it=.false.
             exit
          end if
       end do
       if(count_it) then
          lnodb_size=lnodb_size+knodb
          lmesh%nnodb=max(lmesh%nnodb,knodb)
          iboun_lmesh=iboun_lmesh+1
       end if
    end do

    if(iboun_lmesh>0) then

       call memalloc (  lmesh%nnodb,   node_list, __FILE__,__LINE__)
       call memalloc(   iboun_lmesh, lmesh%lbgeo, __FILE__,__LINE__)
       call memalloc(   iboun_lmesh, lmesh%lbset, __FILE__,__LINE__)

       call lmesh%bound%create(iboun_lmesh)

       iboun_lmesh=2
       do iboun_gmesh=1,gmesh%bound%get_num_pointers()
          bound_iterator = gmesh%bound%create_iterator(iboun_gmesh)
          knodb = bound_iterator%get_size()
          count_it=.true.
          do inode=1,knodb
             call ws_inmap%get(key=int(bound_iterator%get_current(),igp),val=node_list(inode),stat=istat)
             call bound_iterator%next()
             if(istat==key_not_found) then
                count_it=.false.
                exit
             end if
          end do
          if(count_it) then
             call lmesh%bound%sum_to_pointer_index(iboun_lmesh, knodb)
             lmesh%lbgeo(iboun_lmesh)=gmesh%lbgeo(iboun_gmesh)
             lmesh%lbset(iboun_lmesh)=gmesh%lbset(iboun_gmesh)
             iboun_lmesh=iboun_lmesh+1
          end if
       end do

       call lmesh%bound%calculate_header()
       call lmesh%bound%allocate_list_from_pointer()
       do iboun_lmesh=1,lmesh%bound%get_num_pointers()
          bound_iterator = lmesh%bound%create_iterator(iboun_lmesh)
          knodb = bound_iterator%get_size()
          do inode=1,knodb
             call ws_inmap%get(key=int(bound_iterator%get_current(),igp),val=node_list(inode),stat=istat)
			 call bound_iterator%set_current(node_list(inode))
			 call bound_iterator%next()
          enddo
       enddo

       call memfree (node_list, __FILE__,__LINE__)
    end if
    
    call ws_inmap%free
    call el_inmap%free

    call memalloc(number_space_dimensions, lmesh%npoin, lmesh%coord, __FILE__,__LINE__)
    !call map_apply_g2l(nmap, gmesh%ndime, gmesh%coord, lmesh%coord)
    do ipoin=1,num_local_vertices
       lmesh%coord(:,ipoin)=gmesh%coord(:,l2g_vertices(ipoin))
    end do

  end subroutine mesh_g2l
  
  
end module mesh_names
