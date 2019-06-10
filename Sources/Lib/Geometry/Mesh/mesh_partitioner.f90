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
module mesh_partitioner_names
  use types_names
  use list_types_names
  use hash_table_names
  use rcm_renumbering_names
  use memor_names
  use stdio_names
  use metis_names
  use mesh_names
  use mesh_distribution_names
  use mesh_partitioner_parameters_names
  use FPL
  implicit none
# include "debug.i90"
  private
  
  type mesh_partitioner_t
     private
     
     type(mesh_t)             , allocatable :: lmesh(:)
     type(mesh_distribution_t), allocatable :: distr(:)
     
     integer(ip)              :: nparts
     integer(ip)              :: num_levels
     integer(ip), allocatable :: num_parts_x_level (:)
     
     integer(ip)              :: debug                   ! Print (on screen?) info partition
     integer(ip)              :: strat                   ! Partitioning algorithm (part_kway,part_recursive,part_strip,part_rcm_strip)
     integer(ip)              :: metis_option_ufactor    ! Imbalance tol of metis_option_ufactor/1000 + 1
     integer(ip)              :: metis_option_minconn 
     integer(ip)              :: metis_option_contig  
     integer(ip)              :: metis_option_ctype 
     integer(ip)              :: metis_option_iptype
     integer(ip)              :: metis_option_debug
  contains 
     procedure, non_overridable          :: partition_mesh               => mesh_partitioner_partition_mesh
     procedure, non_overridable, private :: set_default_parameter_values => mesh_partitioner_set_default_parameter_values
     procedure, non_overridable, private :: set_parameters_from_pl       => mesh_partitioner_set_parameters_from_pl
     procedure, non_overridable, private :: graph_pt_renumbering
     procedure, non_overridable, private :: mesh_partitioner_write_mesh_parts_dir_path_prefix
     procedure, non_overridable, private :: mesh_partitioner_write_mesh_parts_pl
     generic                             :: write_mesh_parts             => mesh_partitioner_write_mesh_parts_dir_path_prefix, &
                                                                            mesh_partitioner_write_mesh_parts_pl
     procedure, non_overridable          :: free                         => mesh_partitioner_free
  end type mesh_partitioner_t
  
contains

  !=============================================================================
  subroutine mesh_partitioner_partition_mesh ( this, mesh, parameter_list ) 
     implicit none
     class(mesh_partitioner_t), intent(inout) :: this
     type(mesh_t)             , intent(inout) :: mesh
     type(ParameterList_t)    , intent(in)    :: parameter_list
    
     ! Local variables
    type(list_t)                     :: fe_graph         ! Dual graph (to be partitioned)
    type(list_t)                     :: parts_graph      ! Parts graph (to be partitioned)
    integer(ip), allocatable, target :: ldome(:)         ! Part of each element
    type(i1p_t), allocatable         :: ldomp(:)         ! Part of each part (recursively)
    integer(ip), allocatable         :: parts_mapping(:) ! Part of each element
    integer(ip) :: istat, ilevel, jlevel, ipart, itask, num_tasks
     
     
     call this%free()
     call this%set_parameters_from_pl(parameter_list)
     

    ! Generate dual mesh (i.e., list of elements around points)
    call mesh%to_dual()

    ! Create dual (i.e. list of elements around elements)
    call mesh%build_dual_graph(fe_graph)
   
    ! Partition dual graph to assign a domain to each element (in ldome)
    call memalloc (mesh%get_num_cells(), ldome, __FILE__,__LINE__)   
    
    call this%graph_pt_renumbering(fe_graph,ldome)

    allocate(ldomp(this%num_levels), stat=istat); check(istat==0);
    ldomp(1)%p => ldome
    do ilevel=1,this%num_levels-1
       call memallocp(this%num_parts_x_level(ilevel),ldomp(ilevel+1)%p, __FILE__,__LINE__)
       if(this%num_parts_x_level(ilevel+1)>1) then  ! Typically in the last level there is onle one part
          call build_parts_graph (this%num_parts_x_level(ilevel), ldomp(ilevel)%p, fe_graph, parts_graph)
          call fe_graph%free()
          fe_graph = parts_graph
          this%nparts = this%num_parts_x_level(ilevel+1)
          call this%graph_pt_renumbering(parts_graph,ldomp(ilevel+1)%p)
       else
          ldomp(ilevel+1)%p = 1
       end if
       call parts_graph%free()
    end do
    this%nparts = this%num_parts_x_level(1)
    call fe_graph%free()

    num_tasks = 0
    do ilevel=1,this%num_levels
       num_tasks = num_tasks + this%num_parts_x_level(ilevel)
    end do
    !allocate(env(num_tasks), stat=istat); check(istat==0) 
    itask = 0
    call memalloc(this%num_levels,parts_mapping,__FILE__,__LINE__)
    do ilevel=1,this%num_levels
       do ipart = 1, this%num_parts_x_level(ilevel)
          itask = itask+1
          do jlevel = 1 , ilevel - 1 
             parts_mapping(jlevel) = 0
          end do
          parts_mapping(ilevel) = ipart
          do jlevel = ilevel+1 , this%num_levels
             parts_mapping(jlevel) = ldomp(jlevel)%p( parts_mapping(jlevel-1) )
          end do
          !call env(itask)%assign_parts_to_tasks(this%num_levels,this%num_parts_x_level,parts_mapping)
       end do
    end do
    call memfree(parts_mapping,__FILE__,__LINE__)

    this%nparts = this%num_parts_x_level(1)
    do ilevel=1,this%num_levels-1
       call memfreep(ldomp(ilevel+1)%p, __FILE__,__LINE__)
    end do
    deallocate(ldomp, stat=istat); check(istat==0);

    allocate(this%distr(this%nparts), stat=istat); check(istat==0)
    allocate(this%lmesh(this%nparts), stat=istat); check(istat==0) 

    call build_maps(this%nparts, ldome, mesh, this%distr)

    ! Build local meshes and their duals and generate partition adjacency
    do ipart=1,this%nparts
       ! Generate Local mesh
       call mesh_g2l(this%distr(ipart)%num_local_vertices,  &
                     this%distr(ipart)%l2g_vertices,        &
                     this%distr(ipart)%num_local_cells,     &
                     this%distr(ipart)%l2g_cells,           &
                     mesh,                           &
                     this%lmesh(ipart))
       call build_adjacency (mesh, ldome,             &
            &                    ipart,                     &
            &                    this%lmesh(ipart),              &
            &                    this%distr(ipart)%l2g_vertices, &
            &                    this%distr(ipart)%l2g_cells,    &
            &                    this%distr(ipart)%nebou,        &
            &                    this%distr(ipart)%nnbou,        &
            &                    this%distr(ipart)%lebou,        &
            &                    this%distr(ipart)%lnbou,        &
            &                    this%distr(ipart)%pextn,        &
            &                    this%distr(ipart)%lextn,        &
            &                    this%distr(ipart)%lextp )
    end do
    call memfree(ldome,__FILE__,__LINE__)     
  end subroutine mesh_partitioner_partition_mesh
  
  !=============================================================================
  subroutine mesh_partitioner_set_parameters_from_pl(this,parameter_list)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    type(ParameterList_t)    , intent(in)    :: parameter_list
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
  end subroutine mesh_partitioner_set_parameters_from_pl
  
  !=============================================================================
  subroutine mesh_partitioner_set_default_parameter_values(this)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    this%debug = mesh_partitioner_default_debug
    this%strat = mesh_partitioner_default_strat
    this%metis_option_ufactor = mesh_partitioner_default_metis_option_ufactor
    this%metis_option_minconn = mesh_partitioner_default_metis_option_minconn
    this%metis_option_contig = mesh_partitioner_default_metis_option_contig
    this%metis_option_ctype = mesh_partitioner_default_metis_option_ctype
    this%metis_option_iptype = mesh_partitioner_default_metis_option_iptype
    this%metis_option_debug = mesh_partitioner_default_metis_option_debug
  end subroutine mesh_partitioner_set_default_parameter_values
  
  subroutine mesh_partitioner_write_mesh_parts_pl(this, parameter_list)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    type(ParameterList_t)    , intent(in)    :: parameter_list
    
  end subroutine mesh_partitioner_write_mesh_parts_pl
  
  subroutine mesh_partitioner_write_mesh_parts_dir_path_prefix(this, dir_path, prefix)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    character(len=*)         , intent(in)    :: dir_path
    character(len=*)         , intent(in)    :: prefix
    
  end subroutine mesh_partitioner_write_mesh_parts_dir_path_prefix
  
  !=============================================================================
  subroutine mesh_partitioner_free(this)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    if ( allocated(this%num_parts_x_level) ) call memfree(this%num_parts_x_level,__FILE__,__LINE__)
    call this%set_default_parameter_values()
  end subroutine mesh_partitioner_free
  
  
    !================================================================================================
  subroutine build_adjacency ( gmesh, ldome, my_part, lmesh, l2gn, l2ge, &
       &                           nebou, nnbou, lebou, lnbou, pextn, lextn, lextp)
    implicit none
    integer(ip)   , intent(in)  :: my_part
    type(mesh_t)  , intent(in)  :: gmesh,lmesh
    integer(ip)   , intent(in)  :: ldome(gmesh%nelem)
    integer(igp)  , intent(in)  :: l2gn(lmesh%npoin)
    integer(igp)  , intent(in)  :: l2ge(lmesh%nelem)
    integer(ip)   , intent(out) :: nebou
    integer(ip)   , intent(out) :: nnbou
    integer(ip)   , allocatable, intent(inout) ::  lebou(:)    ! List of boundary elements
    integer(ip)   , allocatable, intent(inout) ::  lnbou(:)    ! List of boundary vertices
    integer(ip)   , allocatable, intent(inout) ::  pextn(:)    ! Pointers to the lextn
    integer(igp)  , allocatable, intent(inout) ::  lextn(:)    ! List of (GID of) external neighbors
    integer(ip)   , allocatable, intent(inout) ::  lextp(:)    ! List of parts of external neighbors

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
       write(*,*)  'Local2Global (vertices)'
       do lpoin=1,lmesh%npoin
          write(*,*)  lpoin, l2gn(lpoin)
       end do
    end if

    ! Count boundary vertices
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

    ! List boundary vertices
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
    call external_visited%free()
    call memfree(local_visited,__FILE__,__LINE__)
  end subroutine build_adjacency

  !================================================================================================
  subroutine build_maps( nparts, ldome, femesh, distr )
    ! This routine builds (node and element) partition maps without using the objects
    ! and (unlike parts_sizes, parts_maps, etc.) does not generate a new global numbering.
    implicit none
    integer(ip)                , intent(in)    :: nparts
    type(mesh_t)               , intent(in)    :: femesh
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

    ! Number of vertices of each part and global to local node map (is NOT one to one)
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
  ! of vertices (lconn) of each connected component in m. Can be very
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

             ! Traverse the vertices of the element number j
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

          ! Traverse the vertices of the element number j
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

    call memfree( auxe,__FILE__,__LINE__)
    call memfree( e,__FILE__,__LINE__)
    call memfree( vmarked,__FILE__,__LINE__)
  end subroutine mesh_graph_compute_connected_components

  !=================================================================================================
  subroutine graph_pt_renumbering(this,gp,ldomn,weight)
    !-----------------------------------------------------------------------
    ! This routine computes a nparts-way-partitioning of the input graph gp
    !-----------------------------------------------------------------------
    implicit none
    class(mesh_partitioner_t)       , target, intent(in)    :: this
    type(list_t)                    , target, intent(inout) :: gp
    integer(ip)                     , target, intent(inout) :: ldomn(gp%get_num_pointers())
    integer(ip),            optional, target, intent(in)    :: weight(gp%get_size())

    ! Local variables 
    integer(ip), target      :: kedge
    integer(ip)              :: idumm,iv
    integer(ip), allocatable :: lwork(:)
    integer(ip)              :: i, j, m, k, ipart
    integer(ip), allocatable :: iperm(:)
#ifdef ENABLE_METIS
    integer(c_int),target :: options(0:METIS_NOPTIONS-1)
    integer(c_int),target :: ncon 
    integer(c_int)        :: ierr
#endif    
   
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

    if ( this%strat == part_kway ) then
       options(METIS_OPTION_NUMBERING) = 1
       options(METIS_OPTION_DBGLVL)    = this%metis_option_debug
       
       ! Enforce contiguous partititions
       options(METIS_OPTION_CONTIG)    = this%metis_option_contig
       
       ! Explicitly minimize the maximum degree of the subdomain graph
       options(METIS_OPTION_MINCONN)   = this%metis_option_minconn
       options(METIS_OPTION_UFACTOR)   = this%metis_option_ufactor

       ! Select random (default) or sorted heavy edge matching
       options(METIS_OPTION_CTYPE)     = this%metis_option_ctype
       options(METIS_OPTION_IPTYPE)    = this%metis_option_iptype
       
       ncon = 1
       if(present(weight)) then
          options(METIS_OPTION_NITER) = 100
          ierr = metis_partgraphkway( gp%get_num_pointers_c_loc(), c_loc(ncon), gp%get_pointers_c_loc(), gp%get_list_c_loc() , & 
                                      C_NULL_PTR  , C_NULL_PTR , c_loc(weight) , c_loc(this%nparts), &
                                      C_NULL_PTR  , C_NULL_PTR , c_loc(options), c_loc(kedge), c_loc(ldomn) )
       else
          ierr = metis_partgraphkway( gp%get_num_pointers_c_loc(), c_loc(ncon), gp%get_pointers_c_loc(), gp%get_list_c_loc() , & 
                                      C_NULL_PTR  , C_NULL_PTR , C_NULL_PTR    , c_loc(this%nparts), &
                                      C_NULL_PTR  , C_NULL_PTR , c_loc(options), c_loc(kedge), c_loc(ldomn) )
       end if

       assert(ierr == METIS_OK) 
       
    else if ( this%strat == part_recursive ) then
       options(METIS_OPTION_NUMBERING) = 1
       options(METIS_OPTION_DBGLVL)    = this%metis_option_debug
       options(METIS_OPTION_UFACTOR)   = this%metis_option_ufactor

       ncon = 1 
       ierr = metis_partgraphrecursive( gp%get_num_pointers_c_loc(), c_loc(ncon), gp%get_pointers_c_loc(), gp%get_list_c_loc() , & 
                                        C_NULL_PTR  , C_NULL_PTR , C_NULL_PTR    , c_loc(this%nparts), &
                                        C_NULL_PTR  , C_NULL_PTR , c_loc(options), c_loc(kedge), c_loc(ldomn) )
    end if    
#else
    call enable_metis_error_message
#endif

    if ( this%strat == part_strip ) then
       j = gp%get_num_pointers()
       m = 0
       do ipart=1,this%nparts
          k = j / (this%nparts-ipart+1)
          do i = 1, k
             ldomn(m+i) = ipart
          end do
          m = m + k
          j = j - k
       end do
    else if ( this%strat == part_rcm_strip ) then
       call memalloc ( gp%get_num_pointers(), iperm, __FILE__,__LINE__ )
       call genrcm ( gp, iperm )
       j = gp%get_num_pointers()
       m = 0
       do ipart=1,this%nparts
          k = j / (this%nparts-ipart+1)
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
    integer(ip)                    :: aux, ipoin,inode,knode,kvef_size,lvef_size,istat
    integer(ip)                    :: ielem_lmesh,ielem_gmesh,ivef_lmesh,ivef_gmesh
    integer(ip)                    :: p_ielem_gmesh,p_ipoin_lmesh,p_ipoin_gmesh
    type(list_iterator_t)          :: given_vefs_iterator
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
    ivef_lmesh=0
    lmesh%nnodb=0
    lvef_size=0
    do ivef_gmesh=1,gmesh%given_vefs%get_num_pointers()
       given_vefs_iterator = gmesh%given_vefs%create_iterator(ivef_gmesh)
       kvef_size = given_vefs_iterator%get_size()
       count_it=.true.
       do while(.not. given_vefs_iterator%is_upper_bound())
          call ws_inmap%get(key=int(given_vefs_iterator%get_current(),igp),val=knode,stat=istat)
          call given_vefs_iterator%next()
          if(istat==key_not_found) then
             count_it=.false.
             exit
          end if
       end do
       if(count_it) then
          lvef_size=lvef_size+kvef_size
          lmesh%nnodb=max(lmesh%nnodb,kvef_size)
          ivef_lmesh=ivef_lmesh+1
       end if
    end do

    if(ivef_lmesh>0) then

       call memalloc (  lmesh%nnodb,   node_list, __FILE__,__LINE__)
       call memalloc(   ivef_lmesh, lmesh%lst_vefs_geo, __FILE__,__LINE__)
       call memalloc(   ivef_lmesh, lmesh%lst_vefs_set, __FILE__,__LINE__)

       call lmesh%given_vefs%create(ivef_lmesh)

       ivef_lmesh=1
       do ivef_gmesh=1,gmesh%given_vefs%get_num_pointers()
          given_vefs_iterator = gmesh%given_vefs%create_iterator(ivef_gmesh)
          kvef_size = given_vefs_iterator%get_size()
          count_it=.true.
          do inode=1,kvef_size
             call ws_inmap%get(key=int(given_vefs_iterator%get_current(),igp),val=node_list(inode),stat=istat)
             call given_vefs_iterator%next()
             if(istat==key_not_found) then
                count_it=.false.
                exit
             end if
          end do
          if(count_it) then
             call lmesh%given_vefs%sum_to_pointer_index(ivef_lmesh, kvef_size)
             lmesh%lst_vefs_geo(ivef_lmesh)=gmesh%lst_vefs_geo(ivef_gmesh)
             lmesh%lst_vefs_set(ivef_lmesh)=gmesh%lst_vefs_set(ivef_gmesh)
             ivef_lmesh=ivef_lmesh+1
          end if
       end do

       call lmesh%given_vefs%calculate_header()
       call lmesh%given_vefs%allocate_list_from_pointer()

       ivef_lmesh=1
       do ivef_gmesh=1,gmesh%given_vefs%get_num_pointers()
          given_vefs_iterator = gmesh%given_vefs%create_iterator(ivef_gmesh)
          kvef_size = given_vefs_iterator%get_size()
          count_it=.true.
          do inode=1,kvef_size
             call ws_inmap%get(key=int(given_vefs_iterator%get_current(),igp),val=node_list(inode),stat=istat)
             call given_vefs_iterator%next()
             if(istat==key_not_found) then
                count_it=.false.
                exit
             end if
          end do
          if(count_it) then
             given_vefs_iterator = lmesh%given_vefs%create_iterator(ivef_lmesh)
             do inode=1,kvef_size
                call given_vefs_iterator%set_current(node_list(inode))
                call given_vefs_iterator%next()
             enddo
             ivef_lmesh=ivef_lmesh+1
          end if
       end do
       call memfree (node_list, __FILE__,__LINE__)
    end if
    
    call ws_inmap%free
    call el_inmap%free

    call memalloc(SPACE_DIM, lmesh%npoin, lmesh%coord, __FILE__,__LINE__)
    do ipoin=1,num_local_vertices
       lmesh%coord(:,ipoin)=gmesh%coord(:,l2g_vertices(ipoin))
    end do

  end subroutine mesh_g2l

  subroutine build_parts_graph (nparts, ldome, fe_graph, parts_graph)
    implicit none
    integer(ip)             , intent(in)  :: nparts
    type(list_t)            , intent(in)  :: fe_graph
    integer(ip)             , intent(in)  :: ldome(:)
    type(list_t)            , intent(out) :: parts_graph

    integer(ip)              :: istat,ielem,jelem,ipart,jpart
    integer(ip)              :: num_parts_around, touched
    type(list_iterator_t)                    :: fe_graph_iterator
    type(list_iterator_t)                    :: parts_graph_iterator
    type(position_hash_table_t), allocatable :: visited_parts_touched(:)
    type(hash_table_ip_ip_t)   , allocatable :: visited_parts_numbers(:)

    call parts_graph%create(nparts)

    ! The maximum number of parts around a part can be estimated from the
    ! maximum number of elements connected to an element, that is, by the 
    ! maximum degree of elements graph. Note, however, that it can be bigger.
    num_parts_around=0
    do ielem=1,fe_graph%get_num_pointers()
       fe_graph_iterator = fe_graph%create_iterator(ielem)
       num_parts_around = max(num_parts_around,fe_graph_iterator%get_size())
    end do
    allocate(visited_parts_touched(nparts),stat=istat); assert(istat==0);
    allocate(visited_parts_numbers(nparts),stat=istat); assert(istat==0);
    do ipart=1,nparts
       call visited_parts_touched(ipart)%init(num_parts_around)
       call visited_parts_numbers(ipart)%init(num_parts_around)
    end do

    ! Now compute graph pointers and fill tables
    do ielem=1,fe_graph%get_num_pointers()
       ipart = ldome(ielem)
       fe_graph_iterator = fe_graph%create_iterator(ielem)
       do while(.not.fe_graph_iterator%is_upper_bound())
          jelem = fe_graph_iterator%get_current()
          jpart = ldome(jelem)
          call visited_parts_touched(ipart)%get(key=jpart,val=num_parts_around,stat=istat) ! Touch it (jpart is around ipart)
          if(istat==new_index) then
             call visited_parts_numbers(ipart)%put(key=num_parts_around,val=jpart,stat=istat) ! Store it
             assert(istat==now_stored)
          end if
          call fe_graph_iterator%next()
       end do
    end do
    do ipart=1,nparts
       call parts_graph%sum_to_pointer_index(ipart,visited_parts_touched(ipart)%last())
    end do

    ! Fill graph from tables
    call parts_graph%calculate_header()
    call parts_graph%allocate_list_from_pointer()
    do ipart=1,nparts
       num_parts_around = 0
       parts_graph_iterator = parts_graph%create_iterator(ipart)
       do while(.not.parts_graph_iterator%is_upper_bound())
          num_parts_around = num_parts_around + 1
          call visited_parts_numbers(ipart)%get(key=num_parts_around,val=jpart,stat=istat) 
          assert(istat==key_found)
          call parts_graph_iterator%set_current(jpart)
          call parts_graph_iterator%next()
       end do
       assert(num_parts_around==visited_parts_touched(ipart)%last())
       call visited_parts_touched(ipart)%free() ! This could be done before, as far as the assert is eliminated
       call visited_parts_numbers(ipart)%free()
    end do

    deallocate(visited_parts_touched,stat=istat)
    deallocate(visited_parts_numbers,stat=istat)
  end subroutine build_parts_graph
  
end module mesh_partitioner_names
