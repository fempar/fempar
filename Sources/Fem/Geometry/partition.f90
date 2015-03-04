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
module fem_partition_names
  use types
  use maps_names
  use memor
  use stdio
  implicit none
  private

  !=================================================================
  ! TODO: improve I/O subroutines fem_partition_write and 
  !       fem_partition_read (in particular, use appropiate
  !       formats, instead of *)
  !=================================================================

  ! Some parameters for confortable reading 
  integer(ip), parameter :: element_based=0
  integer(ip), parameter :: vertex_based=1

  ! Alternative implementations of the algorithm which maps
  ! shared points on the interfaces to owner parts
  integer(ip), parameter :: owner_bal_two_parts_min_sev_parts  = 0 
  ! Balance shared points among two parts.
  ! Map arbitrarily shared points among several
  ! parts to the part which min. identifer

  ! Partition data available
  integer(ip), parameter :: basic=10
  integer(ip), parameter :: adjacencies=20
  integer(ip), parameter :: interfaces=30
  integer(ip), parameter :: extended_adjacencies=40

  ! Partition
  type fem_partition

     ! Control data
     integer(ip) ::                &
        ptype = element_based,     &    ! Vertex based or Element based partition
        pinfo = basic                   ! Information available

     ! Basic data
     ! The minimum amount of information needed to describe a distributed
     ! mesh. The rest of data in this type can be generated from this
     ! basic data but it is quite difficult/expensive to do so in parallel.
     integer(ip) ::                &
        ipart  = 1,                &    ! Part identifier
        nparts = 1                      ! Number of parts
     type(map)   ::                &
        nmap,                      &    ! Nodes local to global map
        bmap                            ! Boundaries local to global map
     
     type(map_igp) :: emap

     ! Adjacency data
     ! We also store the list of (global IDs of) external neighbors of near-boundary
     ! entities and the part they belong to. That is, in an element based distribution
     ! this list contains external neighbors of those elements having nodes on the
     ! boundary whereas in a vertex based distribution it contains neighbors of
     ! those vertices belonging to elements on the boundary. 
     integer(ip), allocatable ::   &
        pextn(:),                  &    ! Pointers to the lextr
        lextp(:),                  &    ! List of parts of external neighbors
        lexte(:),                  &    ! Edge information of external neighbors (SB: To be removed)
        lextm(:)                        ! Material of external neighbors         (SB: To be removed)
     
     integer(igp), allocatable ::  &
        lextn(:)                        ! List of (GID of) external neighbors

     ! Together with Adjacency data, minimum set of data required to generate
     ! Interface data (see below) in parallel
     integer(ip) ::  nebou,                     &    ! Number of boundary elements
                     nnbou,                     &    ! Number of boundary nodes
                     nelem,                     &    ! Number of local elements 
                     npoin                           ! Number of local nodes  
     integer(ip), allocatable ::   & 
        lebou(:),                  &  ! List of boundary elements 
        lnbou(:)                      ! List of boundary nodes

     integer(igp), allocatable ::  & 
        l2ge (:),                  &  ! Local2Global for elements 
        l2gn (:)                      ! Local2Global for vertices   

     ! Interface data
     ! Apart from the (local view of the) graph of parts in npadj and lpadj, these data 
     ! structures store the information needed to easily build communication
     ! maps and to build coarse grids. It consists of a (local) list of nodes defining
     ! interface objects (corners, edges and faces), the number and a list of subdomains 
     ! they are shared by and a global numering of them.
     integer(ip) ::                &
        npadj                           ! Number of adjacent parts
     integer(ip), allocatable ::   &
        lpadj(:)                        ! List of adjacent parts

     integer(ip) ::                &
        max_nparts,                &    ! Maximum number of parts around a local object
        nobjs                           ! Number of local objects
     integer(ip), allocatable  ::  &
        lobjs(:,:)                      ! List of local objects
     type(list)  ::                &    ! List of objects on each edge to 
        int_objs                        ! an adjacent part / Interface_objects
     type(map)   ::                &
        omap                            ! Objects local to global map

     ! To add to the documentation, the contents of lobjs are as follows:
     ! lobjs(1,*) : node of the separator tree/physical variable
     ! lobjs(2,*) : first node of the object
     ! lobjs(3,*) : last node of the object
     ! lobjs(4,*) : number of parts around this object
     ! lobjs(5:,*) : list of parts around this object
  end type fem_partition

  interface fem_partition_read
     module procedure fem_partition_read_new, &
          &           fem_partition_read_old
  end interface fem_partition_read

  interface fem_partition_write
     module procedure fem_partition_write_new, &
          &           fem_partition_write_old
  end interface fem_partition_write

  ! Constants
  public :: element_based, vertex_based, & 
            basic, adjacencies, interfaces, extended_adjacencies

  ! Types
  public :: fem_partition

  ! Functions
  public :: fem_partition_free, fem_partition_print,               & 
         &  fem_partition_read, fem_partition_write,               &
         &  fem_partition_compose_name, fem_partition_write_files, &
         &  fem_partition_write_files_new,                         &
         &  fem_partition_read_files,                              &
         &  fem_partition_read_files_new
contains

  !=============================================================================
  subroutine fem_partition_free (f_part)
    !-----------------------------------------------------------------------
    ! This routine deallocates a partition object
    !-----------------------------------------------------------------------
    implicit none

    ! Parameters
    type(fem_partition), intent(inout)  :: f_part

    if (f_part%pinfo==extended_adjacencies) then
       call memfree ( f_part%lebou,__FILE__,__LINE__)
       call memfree ( f_part%lnbou,__FILE__,__LINE__)
       call memfree ( f_part%l2ge,__FILE__,__LINE__)
       call memfree ( f_part%l2gn,__FILE__,__LINE__)
       call memfree ( f_part%pextn ,__FILE__,__LINE__)
       call memfree ( f_part%lextn ,__FILE__,__LINE__)
       call memfree ( f_part%lextp ,__FILE__,__LINE__)
       call memfree ( f_part%lexte ,__FILE__,__LINE__)
       call memfree ( f_part%lextm ,__FILE__,__LINE__)
   else
       call map_free(f_part%nmap)
       call map_free(f_part%emap)
       call map_free(f_part%bmap)
       if (f_part%pinfo==basic) return

       ! Adjacency data
       call memfree ( f_part%pextn ,__FILE__,__LINE__)
       call memfree ( f_part%lextn ,__FILE__,__LINE__)
       call memfree ( f_part%lextp ,__FILE__,__LINE__)
       call memfree ( f_part%lexte ,__FILE__,__LINE__)
       call memfree ( f_part%lextm ,__FILE__,__LINE__)
       if (f_part%pinfo==adjacencies) return

       ! Interface data
       call memfree ( f_part%lpadj,__FILE__,__LINE__)
       call memfree ( f_part%lobjs,__FILE__,__LINE__)
       call memfree ( f_part%int_objs%p ,__FILE__,__LINE__)
       call memfree ( f_part%int_objs%l ,__FILE__,__LINE__)
       call map_free(f_part%omap)
   end if

!    call memfree ( f_part%lpadj,__FILE__,__LINE__)
!    call memfree ( f_part%lobjs,__FILE__,__LINE__)
!    call memfree ( f_part%int_objs%p,__FILE__,__LINE__)
!    call memfree ( f_part%int_objs%l,__FILE__,__LINE__)

!    call memfree ( f_part%pextn ,__FILE__,__LINE__)
!    call memfree ( f_part%lextn ,__FILE__,__LINE__)
!    call memfree ( f_part%lextp ,__FILE__,__LINE__)
!    call memfree ( f_part%lexte ,__FILE__,__LINE__)
!    call memfree ( f_part%lextm ,__FILE__,__LINE__)

!    call map_free (f_part%nmap)
!    call map_free (f_part%emap)
!    call map_free (f_part%omap)
!    call map_free (f_part%bmap)
 
    return

  end subroutine fem_partition_free

  !=============================================================================
  subroutine fem_partition_print (lu_out, part)
    !-----------------------------------------------------------------------
    ! This routine prints a partition object
    !-----------------------------------------------------------------------
    use types
    implicit none


    ! Parameters
    type(fem_partition), intent(in)  :: part
    integer(ip)        , intent(in)  :: lu_out

    ! Local variables
    integer (ip) :: i, j

    if(lu_out>0) then

       write(lu_out,'(a)') '*** begin fem_partition data structure ***'

       write(lu_out,'(a,i10)') 'Number of parts:', &
          &  part%nparts

       write(lu_out,'(a,i10)') 'Number of adjacent parts:', &
          &  part%npadj
       
       write(lu_out,'(10i10)') part%lpadj(1:part%npadj)

       write(lu_out,'(a,i10)') 'Number of elements on the boundary:', &
          &  part%emap%nb

!!$       write(lu_out,'(a)') 'Neighbors and their parts:'
!!$       do i=1,part%emap%nb
!!$          write(lu_out,'(10i10)') part%emap%l2g(part%emap%ni+i), &
!!$               &                  (part%lextn(j),part%lextp(j),j=part%pextn(i),part%pextn(i+1)-1)
!!$       end do

       write(lu_out,'(a,i10)') 'Number of local objects:', &
          &  part%nobjs

       write(lu_out,'(a)') 'List of local objects:'
       do i=1,part%omap%nl
          write(lu_out,'(10i10)') part%omap%l2g(i),part%lobjs(:,i)
       end do

       write(lu_out,'(a)') 'List of interface objects:'
       do i=1,part%npadj
          write(lu_out,'(10i10)') i, &
             & (part%int_objs%l(j),j=part%int_objs%p(i),part%int_objs%p(i+1)-1)
       end do

       write(lu_out,'(a)') '*** end fem_partition data structure ***'

    end if
 
    return

  end subroutine fem_partition_print

  !=============================================================================
  subroutine fem_partition_write_old (file_path, f_part)
    ! Parameters
    character *(*)     , intent(in)  :: file_path
    type(fem_partition), intent(in)  :: f_part
    !-----------------------------------------------------------------------
    ! This routine writes a partition object to file file_path
    !-----------------------------------------------------------------------
    ! Locals 
    integer :: lunio
    lunio = io_open (file_path, 'write')

    write ( lunio, * ) f_part%ptype
    write ( lunio, * ) f_part%ipart
    write ( lunio, * ) f_part%nparts
    write ( lunio, * ) f_part%max_nparts 
    write ( lunio, * ) f_part%nobjs      
    write ( lunio, * ) f_part%npadj       

    write ( lunio, * ) f_part%int_objs%n 
    write ( lunio, * ) f_part%nmap%nl
    write ( lunio, * ) f_part%nmap%ng
    write ( lunio, * ) f_part%nmap%ni
    write ( lunio, * ) f_part%nmap%nb
    write ( lunio, * ) f_part%nmap%ne

    write ( lunio, * ) f_part%emap%nl
    write ( lunio, * ) f_part%emap%ng
    write ( lunio, * ) f_part%emap%ni
    write ( lunio, * ) f_part%emap%nb
    write ( lunio, * ) f_part%emap%ne

    write ( lunio, * ) f_part%omap%nl
    write ( lunio, * ) f_part%omap%ng
    write ( lunio, * ) f_part%omap%ni
    write ( lunio, * ) f_part%omap%nb
    write ( lunio, * ) f_part%omap%ne

    write ( lunio, * ) f_part%nmap%l2g
    write ( lunio, * ) f_part%emap%l2g
    write ( lunio, * ) f_part%omap%l2g

    write ( lunio, * ) f_part%lpadj
    write ( lunio, * ) f_part%lobjs
    write ( lunio, * ) f_part%int_objs%p

    write ( lunio, * ) f_part%int_objs%l

    call io_close (lunio)
  end subroutine fem_partition_write_old

  subroutine fem_partition_write_new (lunio, f_part)
    ! Parameters
    integer            , intent(in) :: lunio
    type(fem_partition), intent(in) :: f_part
    !-----------------------------------------------------------------------
    ! This routine writes a partition to lunio
    !-----------------------------------------------------------------------

    ! Basic data always present
    write ( lunio, '(10i10)' ) f_part%ptype, f_part%pinfo, & 
                               f_part%ipart, f_part%nparts
  
 
    if ( f_part%pinfo==extended_adjacencies) then
       write ( lunio, '(10i10)' ) f_part%nebou
       write ( lunio, '(10i10)' ) f_part%nelem
       write ( lunio, '(10i10)' ) f_part%lebou
       write ( lunio, '(10i10)' ) f_part%l2ge
       write ( lunio, '(10i10)' ) f_part%nnbou
       write ( lunio, '(10i10)' ) f_part%npoin
       write ( lunio, '(10i10)' ) f_part%lnbou
       write ( lunio, '(10i10)' ) f_part%l2gn
       write ( lunio, '(10i10)' ) f_part%pextn
       write ( lunio, '(10i10)' ) f_part%lextn
       write ( lunio, '(10i10)' ) f_part%lextp
       write ( lunio, '(10i10)' ) f_part%lexte
       write ( lunio, '(10i10)' ) f_part%lextm
    else 
       call map_write(lunio,f_part%nmap)
       call map_write(lunio,f_part%emap)
       call map_write(lunio,f_part%bmap)
       if (f_part%pinfo==basic) return

       ! Adjacency data
       write ( lunio, '(10i10)' ) f_part%pextn
       write ( lunio, '(10i10)' ) f_part%lextn
       write ( lunio, '(10i10)' ) f_part%lextp
       write ( lunio, '(10i10)' ) f_part%lexte
       write ( lunio, '(10i10)' ) f_part%lextm
       if (f_part%pinfo==adjacencies) return

       ! Interface data
       write ( lunio, '(10i10)' ) f_part%npadj       
       write ( lunio, '(10i10)' ) f_part%lpadj
       write ( lunio, '(10i10)' ) f_part%max_nparts, f_part%nobjs      
       write ( lunio, '(10i10)' ) f_part%lobjs
       write ( lunio, '(10i10)' ) f_part%int_objs%n 
       write ( lunio, '(10i10)' ) f_part%int_objs%p
       write ( lunio, '(10i10)' ) f_part%int_objs%l
       call map_write(lunio,f_part%omap)
    end if

  end subroutine fem_partition_write_new


  !=============================================================================
  subroutine fem_partition_read_old (file_path, f_part)
    ! Parameters
    character *(*)     , intent(in)     :: file_path
    type(fem_partition), intent(out)    :: f_part
    !-----------------------------------------------------------------------
    ! This routine writes a partition object to file file_path
    !-----------------------------------------------------------------------
    ! Locals 
    integer :: lunio
    lunio = io_open (file_path, 'read', status='old')

    read ( lunio, * ) f_part%ptype
    read ( lunio, * ) f_part%ipart
    read ( lunio, * ) f_part%nparts
    read ( lunio, * ) f_part%max_nparts 
    read ( lunio, * ) f_part%nobjs      
    read ( lunio, * ) f_part%npadj

    read ( lunio, * ) f_part%int_objs%n 

    read ( lunio, * ) f_part%nmap%nl
    read ( lunio, * ) f_part%nmap%ng
    read ( lunio, * ) f_part%nmap%ni
    read ( lunio, * ) f_part%nmap%nb
    read ( lunio, * ) f_part%nmap%ne

    read ( lunio, * ) f_part%emap%nl
    read ( lunio, * ) f_part%emap%ng
    read ( lunio, * ) f_part%emap%ni
    read ( lunio, * ) f_part%emap%nb
    read ( lunio, * ) f_part%emap%ne

    read ( lunio, * ) f_part%omap%nl
    read ( lunio, * ) f_part%omap%ng
    read ( lunio, * ) f_part%omap%ni
    read ( lunio, * ) f_part%omap%nb
    read ( lunio, * ) f_part%omap%ne

    call map_alloc(f_part%nmap%nl, f_part%nmap%ng, f_part%nmap)
    call map_alloc(f_part%emap%nl, f_part%emap%ng, f_part%emap)
    call map_alloc(f_part%omap%nl, f_part%omap%ng, f_part%omap)

    read ( lunio, * ) f_part%nmap%l2g
    read ( lunio, * ) f_part%emap%l2g
    read ( lunio, * ) f_part%omap%l2g

    call memalloc ( f_part%npadj       , f_part%lpadj , __FILE__,__LINE__  )
    call memalloc ( f_part%max_nparts+4, f_part%nobjs ,  f_part%lobjs, __FILE__,__LINE__  )
    call memalloc ( f_part%npadj+1     , f_part%int_objs%p , __FILE__,__LINE__  )

    read ( lunio, * ) f_part%lpadj
    read ( lunio, * ) f_part%lobjs
    read ( lunio, * ) f_part%int_objs%p  

    call memalloc ( f_part%int_objs%p(f_part%npadj+1)-1,f_part%int_objs%l , __FILE__,__LINE__ )

    read ( lunio, * ) f_part%int_objs%l
   
    call io_close (lunio)
  end subroutine fem_partition_read_old

  subroutine fem_partition_read_new (lunio, f_part)
    ! Parameters
    integer(ip)        , intent(in)    :: lunio
    type(fem_partition), intent(inout) :: f_part
    !-----------------------------------------------------------------------
    ! This routine reads a partition
    !-----------------------------------------------------------------------

    ! Basic data always present
    read ( lunio, '(10i10)' ) f_part%ptype, f_part%pinfo, & 
                              f_part%ipart, f_part%nparts

    if (f_part%pinfo==extended_adjacencies) then
       read ( lunio, '(10i10)' ) f_part%nebou       
       read ( lunio, '(10i10)' ) f_part%nelem       
       call memalloc ( f_part%nebou, f_part%lebou,__FILE__,__LINE__  )
       read ( lunio, '(10i10)' ) f_part%lebou

       call memalloc ( f_part%nelem, f_part%l2ge,__FILE__,__LINE__  )
       read ( lunio, '(10i10)' ) f_part%l2ge

       read ( lunio, '(10i10)' ) f_part%nnbou      
       read ( lunio, '(10i10)' ) f_part%npoin   
       call memalloc ( f_part%nnbou, f_part%lnbou,__FILE__,__LINE__  )
       read ( lunio, '(10i10)' ) f_part%lnbou

       call memalloc ( f_part%npoin, f_part%l2gn,__FILE__,__LINE__  )
       read ( lunio, '(10i10)' ) f_part%l2gn

       call memalloc ( f_part%nebou+1, f_part%pextn ,__FILE__,__LINE__  )
       read ( lunio, '(10i10)' ) f_part%pextn

       call memalloc ( f_part%pextn(f_part%nebou+1)-1, f_part%lextn ,__FILE__,__LINE__  )
       call memalloc ( f_part%pextn(f_part%nebou+1)-1, f_part%lextp ,__FILE__,__LINE__  )
       call memalloc ( f_part%pextn(f_part%nebou+1)-1, f_part%lexte , __FILE__,__LINE__  )
       call memalloc ( f_part%pextn(f_part%nebou+1)-1, f_part%lextm , __FILE__,__LINE__  )
       read ( lunio, '(10i10)' ) f_part%lextn
       read ( lunio, '(10i10)' ) f_part%lextp
       read ( lunio, '(10i10)' ) f_part%lexte
       read ( lunio, '(10i10)' ) f_part%lextm
    else 
       call map_read(lunio,f_part%nmap)
       call map_read(lunio,f_part%emap)
       call map_read(lunio,f_part%bmap)
       if (f_part%ptype==basic) return

       ! Adjacency data
       call memalloc ( f_part%emap%nb+1, f_part%pextn , __FILE__,__LINE__  )
       read ( lunio, '(10i10)' ) f_part%pextn

       call memalloc ( f_part%pextn(f_part%emap%nb+1)-1, f_part%lextn ,__FILE__,__LINE__  )
       call memalloc ( f_part%pextn(f_part%emap%nb+1)-1, f_part%lextp ,__FILE__,__LINE__  )
       call memalloc ( f_part%pextn(f_part%emap%nb+1)-1, f_part%lexte ,__FILE__,__LINE__  )
       call memalloc ( f_part%pextn(f_part%emap%nb+1)-1, f_part%lextm , __FILE__,__LINE__  )

       read ( lunio, '(10i10)' ) f_part%lextn
       read ( lunio, '(10i10)' ) f_part%lextp
       read ( lunio, '(10i10)' ) f_part%lexte
       read ( lunio, '(10i10)' ) f_part%lextm
       if (f_part%ptype==adjacencies) return

       ! Interface data
       read ( lunio, '(10i10)' ) f_part%npadj
       call memalloc ( f_part%npadj, f_part%lpadj, __FILE__,__LINE__  )
       read ( lunio, '(10i10)' ) f_part%lpadj

       read ( lunio, '(10i10)' ) f_part%max_nparts, f_part%nobjs      
       call memalloc ( f_part%max_nparts+4, f_part%nobjs ,  f_part%lobjs,__FILE__,__LINE__  )
       read ( lunio, '(10i10)' ) f_part%lobjs
       read ( lunio, '(10i10)' ) f_part%int_objs%n 
       call memalloc ( f_part%int_objs%n+1, f_part%int_objs%p , __FILE__,__LINE__  )
       read ( lunio, '(10i10)' ) f_part%int_objs%p
       call memalloc ( f_part%int_objs%p(f_part%npadj+1)-1,                        f_part%int_objs%l , __FILE__,__LINE__ )
       read ( lunio, '(10i10)' ) f_part%int_objs%l
       call map_read(lunio,f_part%omap)
    end if 
  end subroutine fem_partition_read_new

  !=============================================================================
  subroutine fem_partition_compose_name ( prefix, name ) 
    implicit none
    character *(*), intent(in)        :: prefix 
    character *(*), intent(out)       :: name
    name = trim(prefix) // '.prt'
  end subroutine 

  !=============================================================================

  subroutine fem_partition_write_files ( dir_path, prefix, nparts, parts )
    implicit none
    ! Parameters 
    character *(*), intent(in)       :: dir_path 
    character *(*), intent(in)       :: prefix
    integer(ip)   , intent(in)       :: nparts
    type(fem_partition), intent(in)  :: parts(nparts)
    ! Locals 
    integer (ip)                     :: i
    character(256)                   :: name, rename 

    call fem_partition_compose_name ( prefix, name )

    do i=1,nparts
       rename=name
       call numbered_filename_compose(i,nparts,rename)
       call fem_partition_write ( trim(dir_path) // '/' // trim(rename), parts(i) )
    end do

  end subroutine  fem_partition_write_files

subroutine fem_partition_write_files_new ( dir_path, prefix, nparts, parts )
    implicit none
    ! Parameters
    character *(*), intent(in)       :: dir_path
    character *(*), intent(in)       :: prefix
    integer(ip)   , intent(in)       :: nparts
    type(fem_partition), intent(in)  :: parts(nparts)
    ! Locals
    integer (ip)                     :: i
    character(256)                   :: name, rename

    integer :: lunio

    call fem_partition_compose_name ( prefix, name )

    do i=1,nparts
       rename=name
       call numbered_filename_compose(i,nparts,rename)
       lunio = io_open (trim(dir_path) // '/' // trim(rename))
       call fem_partition_write ( lunio, parts(i) )
       call io_close (lunio)
    end do

  end subroutine  fem_partition_write_files_new

  !=============================================================================
  subroutine fem_nested_partition_read_files ( dir_path, prefix, nparts, num_nested, parts )
    implicit none
    ! Parameters 
    character *(*), intent(in)          :: dir_path 
    character *(*), intent(in)          :: prefix
    integer(ip)   , intent(in)          :: nparts,num_nested
    type(fem_partition), intent(inout)  :: parts(nparts,num_nested)

    ! Locals 
    integer (ip)                        :: i,k
    character(256)                     :: name,rename 

    call fem_partition_compose_name ( prefix, name )

    do k=1,num_nested
       do i=1,nparts
          rename=name
          call numbered_filename_compose(k,num_nested,rename)
          call numbered_filename_compose(i,nparts,rename)
          call fem_partition_read ( trim(dir_path) // '/' // trim(rename), parts(i,k) )
       end do
    end do

  end subroutine  fem_nested_partition_read_files

  subroutine fem_partition_read_files ( dir_path, prefix, nparts, parts )
    implicit none
    ! Parameters 
    character *(*), intent(in)          :: dir_path 
    character *(*), intent(in)          :: prefix
    integer(ip)   , intent(in)          :: nparts
    type(fem_partition), intent(inout)  :: parts(nparts)
    integer         :: lunio

    ! Locals 
    integer (ip)                        :: i
    character(256)                      :: name,rename 

    call fem_partition_compose_name ( prefix, name )

    do i=1,nparts
       rename=name
       call numbered_filename_compose(i,nparts,rename)
       lunio = io_open (trim(dir_path) // '/' // trim(rename))
       call fem_partition_read ( lunio, parts(i) )
       call io_close (lunio)
       !call fem_partition_read ( trim(dir_path) // '/' // trim(rename), parts(i) )
    end do

  end subroutine  fem_partition_read_files
!!$  
!!$  subroutine fem_partition_read_files_new ( dir_path, prefix, nparts, parts )
!!$    implicit none
!!$    ! Parameters 
!!$    character *(*), intent(in)          :: dir_path 
!!$    character *(*), intent(in)          :: prefix
!!$    integer(ip)   , intent(in)          :: nparts
!!$    type(fem_partition), intent(inout)  :: parts(nparts)
!!$
!!$    ! Locals 
!!$    integer (ip)                        :: i,lunio
!!$    character(256)                      :: name,rename 
!!$
!!$    call fem_partition_compose_name ( prefix, name )
!!$
!!$    do i=1,nparts
!!$       rename=name
!!$       call numbered_filename_compose(i,nparts,rename)
!!$       lunio = io_open (trim(dir_path) // '/' // trim(rename))
!!$       call fem_partition_read ( lunio, parts(i) )
!!$    end do
!!$
!!$  end subroutine  fem_partition_read_files_new

  subroutine fem_partition_read_files_new ( dir_path, prefix, nparts, parts )
    implicit none
    ! Parameters 
    character *(*), intent(in)          :: dir_path 
    character *(*), intent(in)          :: prefix
    integer(ip)   , intent(in)          :: nparts
    type(fem_partition), intent(inout)  :: parts(nparts)

    ! Locals 
    integer (ip)                        :: i
    character(256)                      :: name,rename 
    integer(ip)                         :: lunio

    call fem_partition_compose_name ( prefix, name )

    do i=1,nparts
       rename=name
       call numbered_filename_compose(i,nparts,rename)
       lunio = io_open (trim(dir_path) // '/' // trim(rename))
       call fem_partition_read ( lunio, parts(i) )
       call io_close (lunio)
    end do

  end subroutine  fem_partition_read_files_new

end module fem_partition_names
