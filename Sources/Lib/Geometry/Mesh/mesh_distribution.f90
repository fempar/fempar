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
module mesh_distribution_names
  use types_names
  use memor_names
  use stdio_names
  use metis_names
  use mesh_partitioner_parameters_names
  use FPL
  implicit none
# include "debug.i90"
  private

  !> Derived data type which describes the local subdomain corresponding to each
  !> MPI task and its interface with its neighbouring subdomains
  type mesh_distribution_t
     ! private
     integer(ip) ::                &
        ipart,                     & ! Part identifier
        nparts                       ! Number of parts

     integer(ip), allocatable ::   &
        pextn(:),                  & ! Pointers to the lext*
        lextp(:)                     ! List of parts of external neighbors
     
     integer(igp), allocatable ::  &
        lextn(:)                     ! List of (GIDs of) external neighbors

     integer(ip) ::                &
        nebou,                     & ! Number of boundary elements
        nnbou                        ! Number of boundary nodes

     integer(ip), allocatable  ::  & 
        lebou(:),                  &  ! List of boundary elements
        lnbou(:)                      ! List of boundary nodes
        
     integer(ip)               :: num_local_vertices
     integer(igp)              :: num_global_vertices
     integer(igp), allocatable :: l2g_vertices(:)        ! Local2global array of vertices
     
     integer(ip)               :: num_local_cells        ! Number of local cells
     integer(igp)              :: num_global_cells       ! Number of global cells
     integer(igp), allocatable :: l2g_cells(:)           ! Local2global array of cells

   contains
     procedure, non_overridable         :: create => mesh_distribution_create
     procedure, non_overridable         :: free   => mesh_distribution_free

     procedure, non_overridable         :: mesh_distribution_read_dir_path_prefix 
     procedure, non_overridable         :: mesh_distribution_read_file_unit 
     procedure, non_overridable, nopass :: mesh_distribution_compose_name
     generic                            :: read => mesh_distribution_read_dir_path_prefix, &
                                                   mesh_distribution_read_file_unit

     procedure, non_overridable         :: mesh_distribution_write_dir_path_prefix 
     procedure, non_overridable         :: mesh_distribution_write_file_unit 
     generic                            :: write => mesh_distribution_write_dir_path_prefix, &
                                                    mesh_distribution_write_file_unit
    ! Getters
     procedure, non_overridable         :: get_ipart               => mesh_distribution_get_ipart
     procedure, non_overridable         :: get_nparts              => mesh_distribution_get_nparts
     procedure, non_overridable         :: get_pextn               => mesh_distribution_get_pextn
     procedure, non_overridable         :: get_lextp               => mesh_distribution_get_lextp
     procedure, non_overridable         :: get_lextn               => mesh_distribution_get_lextn
     procedure, non_overridable         :: get_nebou               => mesh_distribution_get_nebou
     procedure, non_overridable         :: get_nnbou               => mesh_distribution_get_nnbou
     procedure, non_overridable         :: get_lebou               => mesh_distribution_get_lebou
     procedure, non_overridable         :: get_lnbou               => mesh_distribution_get_lnbou
     procedure, non_overridable         :: get_num_local_vertices  => mesh_distribution_get_num_local_vertices
     procedure, non_overridable         :: get_num_global_vertices => mesh_distribution_get_num_global_vertices
     procedure, non_overridable         :: get_l2g_vertices        => mesh_distribution_get_l2g_vertices
     procedure, non_overridable         :: get_num_local_cells     => mesh_distribution_get_num_local_cells
     procedure, non_overridable         :: get_num_global_cells    => mesh_distribution_get_num_global_cells
     procedure, non_overridable         :: get_l2g_cells           => mesh_distribution_get_l2g_cells

     procedure, non_overridable :: print  => mesh_distribution_print
     procedure, non_overridable :: create_empty => mesh_distribution_create_empty
     procedure, non_overridable :: get_sizes    => mesh_distribution_get_sizes
     procedure, non_overridable :: move_gids    => mesh_distribution_move_gids
     procedure, non_overridable :: move_external_elements_info => mesh_distribution_move_external_elements_info
  end type mesh_distribution_t


  type mesh_distribution_params_t
     integer(ip) :: nparts         = 2    ! nparts
     integer(ip) :: num_levels     = 1    ! nlevels
     integer(ip), allocatable :: num_parts_x_level (:)

     integer(ip) :: debug       = 1    ! Print info partition

     integer(ip) :: strat = part_kway  ! Partitioning algorithm (part_kway,
                                       ! part_recursive,part_strip,part_rcm_strip)

     ! Only applicable to metis 5.0 for both part_kway and part_recursive
     ! Use METIS defaults (i.e., == -1) 30 for part_kway, and 1 for part_recursive
     integer(ip) :: metis_option_ufactor = -1 ! Imbalance tol of x/1000 + 1

     ! Only applicable to metis 5.0 and part_kway
     integer(ip) :: metis_option_minconn = 1 ! (Try to) Minimize maximum degree 
                                             ! of subdomain graph
     integer(ip) :: metis_option_contig  = 1 ! (Try to) Produce partitions 
                                             ! that are contiguous
     
     integer(ip) :: metis_option_ctype  = METIS_CTYPE_RM    ! Random matching
     integer(ip) :: metis_option_iptype = METIS_IPTYPE_GROW ! Grow bisection greedy

     ! Applicable to both metis 4.0 and metis 5.0
     integer(ip) :: metis_option_debug  =  0 
     contains
       procedure, non_overridable :: get_parameters_from_fpl =>  mesh_distribution_get_parameters_from_fpl
       procedure, non_overridable :: free => mesh_distribution_parameters_free
  end type mesh_distribution_params_t

  ! Types
  public :: mesh_distribution_t, mesh_distribution_params_t

  ! Functions
  public :: mesh_distribution_write_files
  public :: mesh_distribution_read_files

contains

  subroutine mesh_distribution_create (this, &
                                       ipart, &
                                       nparts, &
                                       num_local_cells, &
                                       num_global_cells, &
                                       nebou, &
                                       num_local_vertices, &
                                       num_global_vertices, &
                                       nnbou, &
                                       l2g_cells, &
                                       l2g_vertices, &
                                       lebou, &
                                       lnbou, &
                                       pextn, &
                                       lextn, &
                                       lextp)
    implicit none
    class(mesh_distribution_t), intent(inout) :: this
    integer(ip)               , intent(in)    :: ipart
    integer(ip)               , intent(in)    :: nparts
    integer(ip)               , intent(in)    :: num_local_cells
    integer(igp)              , intent(in)    :: num_global_cells
    integer(ip)               , intent(in)    :: nebou
    integer(ip)               , intent(in)    :: num_local_vertices
    integer(igp)              , intent(in)    :: num_global_vertices
    integer(ip)               , intent(in)    :: nnbou
    integer(igp)              , intent(in)    :: l2g_cells(num_local_cells)
    integer(igp)              , intent(in)    :: l2g_vertices(num_local_vertices)
    integer(ip)               , intent(in)    :: lebou(nebou)
    integer(ip)               , intent(in)    :: lnbou(nnbou)
    integer(ip)               , intent(in)    :: pextn(nebou+1)
    integer(igp)              , intent(in)    :: lextn(pextn(nebou+1)-1)
    integer(ip)               , intent(in)    :: lextp(pextn(nebou+1)-1)
    call this%free()
    
    this%ipart               = ipart 
    this%nparts              = nparts 
    this%num_local_cells     = num_local_cells
    this%num_global_cells    = num_global_cells
    this%nebou               = nebou
    this%num_local_vertices  = num_local_vertices
    this%num_global_vertices = num_global_vertices
    this%nnbou               = nnbou
    call memalloc(size(l2g_cells)   , this%l2g_cells,__FILE__,__LINE__)
    this%l2g_cells = l2g_cells
    call memalloc(size(l2g_vertices), this%l2g_vertices,__FILE__,__LINE__)
    this%l2g_vertices = l2g_vertices
    call memalloc(size(lebou)       , this%lebou,__FILE__,__LINE__)
    this%lebou = lebou
    call memalloc(size(lnbou)       , this%lnbou,__FILE__,__LINE__)
    this%lnbou = lnbou
    call memalloc(size(pextn)       , this%pextn,__FILE__,__LINE__)
    this%pextn = pextn
    call memalloc(size(lextn)       , this%lextn,__FILE__,__LINE__)
    this%lextn = lextn
    call memalloc(size(lextp)       , this%lextp,__FILE__,__LINE__)
    this%lextp = lextp
  end subroutine mesh_distribution_create 

  !=============================================================================
  subroutine mesh_distribution_free (this)
    implicit none

    ! Parameters
    class(mesh_distribution_t), intent(inout)  :: this
    if(allocated(this%lebou)) call memfree ( this%lebou,__FILE__,__LINE__)
    if(allocated(this%lnbou)) call memfree ( this%lnbou,__FILE__,__LINE__)
    if(allocated(this%pextn)) call memfree ( this%pextn ,__FILE__,__LINE__)
    if(allocated(this%lextn)) call memfree ( this%lextn ,__FILE__,__LINE__)
    if(allocated(this%lextp)) call memfree ( this%lextp ,__FILE__,__LINE__)
    if(allocated(this%l2g_vertices)) call memfree ( this%l2g_vertices, __FILE__,__LINE__)
    if(allocated(this%l2g_cells)) call memfree ( this%l2g_cells, __FILE__,__LINE__)
    this%ipart               = 0 
    this%nparts              = 0 
    this%num_local_cells     = 0
    this%num_global_cells    = 0
    this%nebou               = 0
    this%num_local_vertices  = 0
    this%num_global_vertices = 0
    this%nnbou               = 0   
  end subroutine mesh_distribution_free

  subroutine mesh_distribution_read_dir_path_prefix (this,  dir_path, prefix)
     implicit none 
     ! Parameters
     class(mesh_distribution_t), intent(inout) :: this
     character(*)              , intent(in)    :: dir_path
     character(*)              , intent(in)    :: prefix
     ! Locals
     integer(ip)                    :: lunio
     character(len=:), allocatable  :: name

     ! Read mesh
     call this%mesh_distribution_compose_name ( prefix, name )
     lunio = io_open( trim(dir_path)//'/'//trim(name), 'read', status='old' ); check(lunio>0)
     call this%read(lunio)
     call io_close(lunio)
  end subroutine mesh_distribution_read_dir_path_prefix

  subroutine mesh_distribution_read_file_unit (this, lunio)
    implicit none 
    ! Parameters
    class(mesh_distribution_t), intent(inout) :: this
    integer(ip)               , intent(in)    :: lunio
    !-----------------------------------------------------------------------
    ! This subroutine reads a mesh_distribution object
    !-----------------------------------------------------------------------

    call this%free()

    read ( lunio, '(10i10)' ) this%ipart, this%nparts

    read ( lunio, '(10i10)' ) this%nebou       
    call memalloc ( this%nebou, this%lebou,__FILE__,__LINE__  )
    read ( lunio, '(10i10)' ) this%lebou
        
    read ( lunio, '(10i10)' ) this%nnbou      
    call memalloc ( this%nnbou, this%lnbou,__FILE__,__LINE__  )
    read ( lunio, '(10i10)' ) this%lnbou
    
    
    call memalloc ( this%nebou+1, this%pextn ,__FILE__,__LINE__  )
    read ( lunio, '(10i10)' ) this%pextn
    
    call memalloc ( this%pextn(this%nebou+1)-1, this%lextn ,__FILE__,__LINE__  )
    call memalloc ( this%pextn(this%nebou+1)-1, this%lextp ,__FILE__,__LINE__  )
    read ( lunio, '(10i10)' ) this%lextn
    read ( lunio, '(10i10)' ) this%lextp

    read ( lunio, '(10i10)' ) this%num_local_vertices, &
                              this%num_global_vertices
    if(this%num_local_vertices>0) then
       call memalloc(this%num_local_vertices, this%l2g_vertices, __FILE__, __LINE__)
       read ( lunio,'(10i10)') this%l2g_vertices
    end if

    read ( lunio, '(10i10)' ) this%num_local_cells, &
                              this%num_global_cells
    if(this%num_local_cells>0) then
       call memalloc(this%num_local_cells, this%l2g_cells, __FILE__, __LINE__)
       read ( lunio,'(10i10)') this%l2g_cells
    end if
  end subroutine mesh_distribution_read_file_unit


  !=============================================================================
    pure function mesh_distribution_get_ipart(this) result(ipart)
        class(mesh_distribution_t), intent(in) :: this
        integer(ip)                            :: ipart
        ipart = this%ipart
    end function mesh_distribution_get_ipart

  !=============================================================================
    pure function mesh_distribution_get_nparts(this) result(nparts)
        class(mesh_distribution_t), intent(in) :: this
        integer(ip)                            :: nparts
        nparts = this%nparts
    end function mesh_distribution_get_nparts


  !=============================================================================
    function mesh_distribution_get_pextn(this) result(pextn)
        class(mesh_distribution_t), target, intent(in) :: this
        integer(ip),                pointer            :: pextn(:)
        pextn => this%pextn
    end function mesh_distribution_get_pextn


  !=============================================================================
    function mesh_distribution_get_lextp(this) result(lextp)
        class(mesh_distribution_t), target, intent(in) :: this
        integer(ip),                pointer            :: lextp(:)
        lextp => this%lextp
    end function mesh_distribution_get_lextp


  !=============================================================================
    function mesh_distribution_get_lextn(this) result(lextn)
        class(mesh_distribution_t), target, intent(in) :: this
        integer(igp),               pointer            :: lextn(:)
        lextn => this%lextn
    end function mesh_distribution_get_lextn


  !=============================================================================
    pure function mesh_distribution_get_nebou(this) result(nebou)
        class(mesh_distribution_t), intent(in) :: this
        integer(ip)                            :: nebou
        nebou = this%nebou
    end function mesh_distribution_get_nebou


  !=============================================================================
    pure function mesh_distribution_get_nnbou(this) result(nnbou)
        class(mesh_distribution_t), intent(in) :: this
        integer(ip)                            :: nnbou
        nnbou = this%nnbou
    end function mesh_distribution_get_nnbou


  !=============================================================================
    function mesh_distribution_get_lebou(this) result(lebou)
        class(mesh_distribution_t), target, intent(in) :: this
        integer(ip),                pointer            :: lebou(:)
        lebou => this%lebou
    end function mesh_distribution_get_lebou


  !=============================================================================
    function mesh_distribution_get_lnbou(this) result(lnbou)
        class(mesh_distribution_t), target, intent(in) :: this
        integer(ip),                pointer            :: lnbou(:)
        lnbou => this%lnbou
    end function mesh_distribution_get_lnbou


  !=============================================================================
    pure function mesh_distribution_get_num_local_vertices(this) result(num_local_vertices)
        class(mesh_distribution_t), intent(in) :: this
        integer(ip)                            :: num_local_vertices
        num_local_vertices = this%num_local_vertices
    end function mesh_distribution_get_num_local_vertices


  !=============================================================================
    pure function mesh_distribution_get_num_global_vertices(this) result(num_global_vertices)
        class(mesh_distribution_t), intent(in) :: this
        integer(igp)                           :: num_global_vertices
        num_global_vertices = this%num_global_vertices
    end function mesh_distribution_get_num_global_vertices


  !=============================================================================
    function mesh_distribution_get_l2g_vertices(this) result(l2g_vertices)
        class(mesh_distribution_t), target, intent(in) :: this
        integer(igp),               pointer            :: l2g_vertices(:)
        l2g_vertices = this%l2g_vertices
    end function mesh_distribution_get_l2g_vertices


  !=============================================================================
    pure function mesh_distribution_get_num_local_cells(this) result(num_local_cells)
        class(mesh_distribution_t), intent(in) :: this
        integer(ip)                            :: num_local_cells
        num_local_cells = this%num_local_cells
    end function mesh_distribution_get_num_local_cells


  !=============================================================================
    pure function mesh_distribution_get_num_global_cells(this) result(num_global_cells)
        class(mesh_distribution_t), intent(in) :: this
        integer(igp)                           :: num_global_cells
        num_global_cells = this%num_global_cells
    end function mesh_distribution_get_num_global_cells


  !=============================================================================
    function mesh_distribution_get_l2g_cells(this) result(l2g_cells)
        class(mesh_distribution_t), target, intent(in) :: this
        integer(igp),               pointer            :: l2g_cells(:)
        l2g_cells => this%l2g_cells
    end function mesh_distribution_get_l2g_cells


  !=============================================================================
  subroutine mesh_distribution_compose_name ( prefix, name ) 
    implicit none
    character (len=*)             , intent(in)    :: prefix 
    character (len=:), allocatable, intent(inout) :: name
    name = trim(prefix) // '.prt'
  end subroutine mesh_distribution_compose_name 


  !=============================================================================
  subroutine mesh_distribution_get_sizes(this,ipart,nparts)
    class(mesh_distribution_t), intent(inout) :: this
    integer(ip), intent(inout) :: ipart,nparts
    ipart=this%ipart
    nparts=this%nparts
  end subroutine mesh_distribution_get_sizes
  !=============================================================================
  subroutine mesh_distribution_move_gids(this,cells_gid,vefs_gid)
    class(mesh_distribution_t), intent(inout) :: this
    integer(igp), intent(inout), allocatable :: vefs_gid(:)
    integer(igp), intent(inout), allocatable :: cells_gid(:)
    call memmovealloc(this%l2g_vertices,vefs_gid,__FILE__,__LINE__)
    call memmovealloc(this%l2g_cells,cells_gid,__FILE__,__LINE__)
  end subroutine mesh_distribution_move_gids
  !=============================================================================
  subroutine mesh_distribution_move_external_elements_info(this,nebou,lebou,pextn,lextn,lextp)
    class(mesh_distribution_t), intent(inout) :: this
    integer(ip), intent(inout)   :: nebou
    integer(ip), intent(inout), allocatable :: lebou(:)
    integer(ip), intent(inout), allocatable :: pextn(:)
    integer(igp), intent(inout), allocatable :: lextn(:)
    integer(ip), intent(inout), allocatable :: lextp(:)
    nebou=this%nebou
    call memmovealloc(this%lebou,lebou,__FILE__,__LINE__)
    call memmovealloc(this%pextn,pextn,__FILE__,__LINE__)
    call memmovealloc(this%lextn,lextn,__FILE__,__LINE__)
    call memmovealloc(this%lextp,lextp,__FILE__,__LINE__)
  end subroutine mesh_distribution_move_external_elements_info
  !=============================================================================
  subroutine mesh_distribution_create_empty(this)
    class(mesh_distribution_t), intent(inout) :: this
    call memalloc ( 1, this%pextn ,__FILE__,__LINE__  )
    this%pextn(1) = 1
  end subroutine mesh_distribution_create_empty
  
  !=============================================================================
  subroutine mesh_distribution_print (msh_dist, lu_out)
    !-----------------------------------------------------------------------
    ! This subroutine prints a mesh_distribution object
    !-----------------------------------------------------------------------
    implicit none

    ! Parameters
    integer(ip)              , intent(in)  :: lu_out
    class(mesh_distribution_t), intent(in)  :: msh_dist

    ! Local variables
    integer (ip) :: i, j

    if(lu_out>0) then
       write(lu_out,'(a)') '*** begin mesh_distribution data structure ***'

       write(lu_out,'(a,i10)') 'Number of parts:', &
           &  msh_dist%nparts

       write(lu_out,'(a,i10)') 'Number of elements on the boundary:', &
          &  msh_dist%nebou

       write(lu_out,'(a,i10)') 'Number of neighbours:', &
          &  msh_dist%pextn(msh_dist%nebou+1)-msh_dist%pextn(1)

       write(lu_out,'(a)') 'GEIDs of boundary elements:'
       do i=1,msh_dist%nebou
          write(lu_out,'(10i10)') msh_dist%l2g_cells(msh_dist%lebou(i))
       end do

       write(lu_out,'(a)') 'GEIDs of neighbors:'
       do i=1,msh_dist%nebou
          write(lu_out,'(10i10)') (msh_dist%lextn(j),j=msh_dist%pextn(i),msh_dist%pextn(i+1)-1)
       end do

       write(lu_out,'(a)') 'Parts of neighbours:'
       do i=1,msh_dist%nebou
          write(lu_out,'(10i10)') (msh_dist%lextp(j),j=msh_dist%pextn(i),msh_dist%pextn(i+1)-1)
       end do

       write(lu_out,'(a)') '*** end mesh_distribution data structure ***'
    end if
  end subroutine mesh_distribution_print

  !=============================================================================
  subroutine mesh_distribution_write_dir_path_prefix (this, dir_path, prefix)
     implicit none 
     ! Parameters
     class(mesh_distribution_t), intent(inout) :: this
     character(*)              , intent(in)    :: dir_path
     character(*)              , intent(in)    :: prefix
     ! Locals
     integer(ip) :: lunio
     character(len=:), allocatable  :: name

     ! Read mesh
     call this%mesh_distribution_compose_name ( prefix, name )
     lunio = io_open( trim(dir_path)//'/'//trim(name), 'write' ); check(lunio>0)
     call this%write(lunio)
     call io_close(lunio)
  end subroutine mesh_distribution_write_dir_path_prefix

  !=============================================================================
  subroutine mesh_distribution_write_file_unit (this, lunio)
    implicit none
    ! Parameters
    class(mesh_distribution_t), intent(in) :: this
    integer                   , intent(in) :: lunio
    write ( lunio, '(10i10)' ) this%ipart, this%nparts
    write ( lunio, '(10i10)' ) this%nebou
    write ( lunio, '(10i10)' ) this%lebou
    write ( lunio, '(10i10)' ) this%nnbou
    write ( lunio, '(10i10)' ) this%lnbou
    write ( lunio, '(10i10)' ) this%pextn
    write ( lunio, '(10i10)' ) this%lextn
    write ( lunio, '(10i10)' ) this%lextp
    write ( lunio, '(10i10)' ) this%num_local_vertices, &
                               this%num_global_vertices
    if(this%num_local_vertices>0) write ( lunio,'(10i10)') this%l2g_vertices
    
    write ( lunio, '(10i10)' ) this%num_local_cells, &
                               this%num_global_cells
    if(this%num_local_cells>0) write ( lunio,'(10i10)') this%l2g_cells
  end subroutine mesh_distribution_write_file_unit

  !=============================================================================
  subroutine mesh_distribution_write_files ( parameter_list, parts )
    implicit none
    ! Parameters
    type(ParameterList_t)    , intent(in) :: parameter_list
    type(mesh_distribution_t), intent(in)  :: parts(:)

    ! Locals
    integer(ip)                   :: nparts
    integer(ip)                   :: istat
    logical                       :: is_present
    character(len=:), allocatable :: dir_path
    character(len=:), allocatable :: prefix
    character(len=:), allocatable :: name, rename
    integer(ip)                   :: lunio
    integer(ip)                   :: i

    nparts = size(parts)

    ! Mandatory parameters
    assert(parameter_list%isAssignable(dir_path_out_key, 'string'))
    istat = parameter_list%GetAsString(key = dir_path_out_key, String = dir_path)
    assert(istat == 0)
    
    assert(parameter_list%isAssignable(prefix_key, 'string'))
    istat = istat + parameter_list%GetAsString(key = prefix_key, String = prefix)
    assert(istat==0)

    call mesh_distribution_compose_name ( prefix, name )
    
    do i=1,nparts
       rename=name
       call numbered_filename_compose(i,nparts,rename)
       lunio = io_open (trim(dir_path) // '/' // trim(rename)); check(lunio>0)
       call parts(i)%write (lunio)
       call io_close (lunio)
    end do
  end subroutine  mesh_distribution_write_files

  !=============================================================================
  subroutine mesh_distribution_read_files ( dir_path, prefix, nparts, parts )
    implicit none
    ! Parameters 
    character(*), intent(in)    :: dir_path 
    character(*), intent(in)    :: prefix
    integer(ip)     , intent(in)    :: nparts
    type(mesh_distribution_t), intent(inout)  :: parts(nparts)

    ! Locals 
    integer (ip)                        :: i
    character(len=:), allocatable       :: name,rename ! Deferred-length allocatable character arrays
    integer(ip)                         :: lunio

    call mesh_distribution_compose_name ( prefix, name )

    do i=1,nparts
       rename=name
       call numbered_filename_compose(i,nparts,rename)
       lunio = io_open (trim(dir_path) // '/' // trim(rename)); check(lunio>0)
       call parts(i)%read (lunio)
       call io_close (lunio)
    end do
  end subroutine  mesh_distribution_read_files

  subroutine mesh_distribution_get_parameters_from_fpl(this,parameter_list)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine generates geometry data to construct a structured mesh                      !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(mesh_distribution_params_t), intent(inout) :: this
    type(ParameterList_t)            , intent(in)    :: parameter_list
    ! Locals
    integer(ip)              :: istat
    integer(ip), allocatable :: param_size(:), param(:)

    ! Mandatory parameters: either nparts or num_levels
    assert(parameter_list%isPresent(key = num_parts_key).or.parameter_list%isPresent(key = num_levels_distribution_key))
    if( parameter_list%isPresent(num_parts_key)) then
       assert(parameter_list%isAssignable(num_parts_key, this%nparts))
       istat = parameter_list%get(key = num_parts_key , value = this%nparts)
       assert(istat==0)
    end if
    if( parameter_list%isPresent(num_levels_distribution_key) ) then
       assert(parameter_list%isAssignable(num_levels_distribution_key, this%num_levels))
       istat = parameter_list%get(key = num_levels_distribution_key  , value = this%num_levels)
       assert(istat==0)
       
       assert(parameter_list%isPresent(key = num_parts_x_level_key ))
       assert( parameter_list%GetDimensions(key = num_parts_x_level_key) == 1)

       ! Get the array using the local variable
       istat =  parameter_list%GetShape(key = num_parts_x_level_key, shape = param_size ); check(istat==0)
       call memalloc(param_size(1), param,__FILE__,__LINE__)
       assert(parameter_list%isAssignable(num_parts_x_level_key, param))
       istat = parameter_list%get(key = num_parts_x_level_key, value = param)
       assert(istat==0)

       call memalloc(this%num_levels, this%num_parts_x_level,__FILE__,__LINE__)
       this%num_parts_x_level = param(1:this%num_levels)
       call memfree(param,__FILE__,__LINE__)

       this%nparts = this%num_parts_x_level(1)
    else
       this%num_levels=1
       call memalloc(this%num_levels, this%num_parts_x_level,__FILE__,__LINE__)
       this%num_parts_x_level(1)=this%nparts
    end if

    ! Optional paramters
    if( parameter_list%isPresent(debug_key) ) then
       assert(parameter_list%isAssignable(debug_key, this%debug))
       istat = parameter_list%get(key = debug_key  , value = this%debug)
       assert(istat==0)
    end if

    if( parameter_list%isPresent(strategy_key) ) then
       assert(parameter_list%isAssignable(strategy_key, this%strat))
       istat = parameter_list%get(key = strategy_key  , value = this%strat)
       assert(istat==0)
       assert(this%strat==part_kway.or.this%strat==part_recursive.or.this%strat==part_strip.or.this%strat==part_rcm_strip)
    end if

    if( parameter_list%isPresent(metis_option_debug_key) ) then
       assert(parameter_list%isAssignable(metis_option_debug_key, this%metis_option_debug))
       istat = parameter_list%get(key = metis_option_debug_key  , value = this%metis_option_debug)
       check(istat==0)
    end if

    if( parameter_list%isPresent(metis_option_ufactor_key) ) then
       assert(parameter_list%isAssignable(metis_option_ufactor_key, this%metis_option_ufactor))
       istat = parameter_list%get(key = metis_option_ufactor_key, value = this%metis_option_ufactor)
       assert(istat==0)
    end if

    if( parameter_list%isPresent(metis_option_minconn_key) ) then
       assert(parameter_list%isAssignable(metis_option_minconn_key, this%metis_option_minconn))
       istat = parameter_list%get(key = metis_option_minconn_key, value = this%metis_option_minconn)
       check(istat==0)
    end if

    if( parameter_list%isPresent(metis_option_contig_key) ) then
       assert(parameter_list%isAssignable(metis_option_contig_key, this%metis_option_contig))
       istat = parameter_list%get(key = metis_option_contig_key , value = this%metis_option_contig)
       assert(istat==0)
    end if

    if( parameter_list%isPresent(metis_option_ctype_key) ) then
       assert(parameter_list%isAssignable(metis_option_ctype_key, this%metis_option_ctype))
       istat = parameter_list%get(key = metis_option_ctype_key  , value = this%metis_option_ctype)
       assert(istat==0)
    end if

  end subroutine mesh_distribution_get_parameters_from_fpl

  !=============================================================================
  subroutine mesh_distribution_parameters_free(this)
    implicit none
    class(mesh_distribution_params_t), intent(inout) :: this
    call memfree(this%num_parts_x_level,__FILE__,__LINE__)
  end subroutine mesh_distribution_parameters_free


end module mesh_distribution_names
