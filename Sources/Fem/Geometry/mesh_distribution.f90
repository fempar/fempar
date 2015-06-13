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
module fem_mesh_distribution_names
  use types
  use memor
  use stdio
  use maps_names
  implicit none
  private

  ! Data required to describe on each MPI task the distribution of the mesh
  type fem_mesh_distribution
     integer(ip) ::                &
        ipart  = 1,                &    ! Part identifier
        nparts = 1                      ! Number of parts

     integer(ip), allocatable ::   &
        pextn(:),                  &    ! Pointers to the lext*
        lextp(:)                        ! List of parts of external neighbors
     
     integer(igp), allocatable ::   &
        lextn(:)                        ! List of (GIDs of) external neighbors

     integer(ip) ::  nebou,        &    ! Number of boundary elements
                     nnbou              ! Number of boundary nodes 

     integer(ip), allocatable ::   & 
        lebou(:),                  &  ! List of boundary elements 
        lnbou(:)                      ! List of boundary nodes

     type(map_igp) ::  & 
        emap,                  &  ! Local2Global for elements 
        nmap                      ! Local2Global for vertices
  end type fem_mesh_distribution

  ! Types
  public :: fem_mesh_distribution

  ! Functions
  public :: fem_mesh_distribution_free, fem_mesh_distribution_print,               & 
         &  fem_mesh_distribution_read, fem_mesh_distribution_write,               &
         &  fem_mesh_distribution_compose_name, fem_mesh_distribution_write_files, &
         &  fem_mesh_distribution_read_files
contains

  !=============================================================================
  subroutine fem_mesh_distribution_free (f_msh_dist)
    !-----------------------------------------------------------------------
    ! This subroutine deallocates a mesh_distribution object
    !-----------------------------------------------------------------------
    implicit none

    ! Parameters
    type(fem_mesh_distribution), intent(inout)  :: f_msh_dist

    call memfree ( f_msh_dist%lebou,__FILE__,__LINE__)
    call memfree ( f_msh_dist%lnbou,__FILE__,__LINE__)
    call memfree ( f_msh_dist%pextn ,__FILE__,__LINE__)
    call memfree ( f_msh_dist%lextn ,__FILE__,__LINE__)
    call memfree ( f_msh_dist%lextp ,__FILE__,__LINE__)

    call map_free (f_msh_dist%nmap)
    call map_free (f_msh_dist%emap)

  end subroutine fem_mesh_distribution_free

  !=============================================================================
  subroutine fem_mesh_distribution_print (lu_out, msh_dist)
    !-----------------------------------------------------------------------
    ! This subroutine prints a mesh_distribution object
    !-----------------------------------------------------------------------
    use types
    implicit none

    ! Parameters
    integer(ip)                , intent(in)  :: lu_out
    type(fem_mesh_distribution), intent(in)  :: msh_dist

    ! Local variables
    integer (ip) :: i, j

    if(lu_out>0) then

       write(lu_out,'(a)') '*** begin fem_mesh_distribution data structure ***'

       write(lu_out,'(a,i10)') 'Number of parts:', &
          &  msh_dist%nparts

       write(lu_out,'(a,i10)') 'Number of elements on the boundary:', &
          &  msh_dist%nebou

       write(lu_out,'(a)') 'GEIDs of boundary elements, neighbors and their parts:'
       do i=1,msh_dist%nebou
          write(lu_out,'(10i10)') msh_dist%emap%l2g(msh_dist%lebou(i)), &
               &                  (msh_dist%lextn(j),msh_dist%lextp(j),j=msh_dist%pextn(i),msh_dist%pextn(i+1)-1)
       end do

       write(lu_out,'(a)') '*** end fem_mesh_distribution data structure ***'

    end if
 
  end subroutine fem_mesh_distribution_print

  subroutine fem_mesh_distribution_write (lunio, f_msh_dist)
    ! Parameters
    integer            , intent(in) :: lunio
    type(fem_mesh_distribution), intent(in) :: f_msh_dist
    !-----------------------------------------------------------------------
    ! This subroutine writes a mesh_distribution to lunio
    !-----------------------------------------------------------------------

    write ( lunio, '(10i10)' ) f_msh_dist%ipart, f_msh_dist%nparts
    write ( lunio, '(10i10)' ) f_msh_dist%nebou
    write ( lunio, '(10i10)' ) f_msh_dist%lebou
    write ( lunio, '(10i10)' ) f_msh_dist%nnbou
    write ( lunio, '(10i10)' ) f_msh_dist%lnbou
    write ( lunio, '(10i10)' ) f_msh_dist%pextn
    write ( lunio, '(10i10)' ) f_msh_dist%lextn
    write ( lunio, '(10i10)' ) f_msh_dist%lextp

    call map_write (lunio, f_msh_dist%nmap)
    call map_write (lunio, f_msh_dist%emap)

  end subroutine fem_mesh_distribution_write

  subroutine fem_mesh_distribution_read (lunio, f_msh_dist)
    ! Parameters
    integer(ip)        , intent(in)            :: lunio
    type(fem_mesh_distribution), intent(inout) :: f_msh_dist
    !-----------------------------------------------------------------------
    ! This subroutine reads a mesh_distribution object
    !-----------------------------------------------------------------------

    read ( lunio, '(10i10)' ) f_msh_dist%ipart, f_msh_dist%nparts

    read ( lunio, '(10i10)' ) f_msh_dist%nebou       
    call memalloc ( f_msh_dist%nebou, f_msh_dist%lebou,__FILE__,__LINE__  )
    read ( lunio, '(10i10)' ) f_msh_dist%lebou
        
    read ( lunio, '(10i10)' ) f_msh_dist%nnbou      
    call memalloc ( f_msh_dist%nnbou, f_msh_dist%lnbou,__FILE__,__LINE__  )
    read ( lunio, '(10i10)' ) f_msh_dist%lnbou
    
    
    call memalloc ( f_msh_dist%nebou+1, f_msh_dist%pextn ,__FILE__,__LINE__  )
    read ( lunio, '(10i10)' ) f_msh_dist%pextn
    
    call memalloc ( f_msh_dist%pextn(f_msh_dist%nebou+1)-1, f_msh_dist%lextn ,__FILE__,__LINE__  )
    call memalloc ( f_msh_dist%pextn(f_msh_dist%nebou+1)-1, f_msh_dist%lextp ,__FILE__,__LINE__  )
    read ( lunio, '(10i10)' ) f_msh_dist%lextn
    read ( lunio, '(10i10)' ) f_msh_dist%lextp

    call map_read (lunio, f_msh_dist%nmap)
    call map_read (lunio, f_msh_dist%emap)

  end subroutine fem_mesh_distribution_read

  !=============================================================================
  subroutine fem_mesh_distribution_compose_name ( prefix, name ) 
    implicit none
    character (len=*), intent(in)                 :: prefix 
    character (len=:), allocatable, intent(inout) :: name
    name = trim(prefix) // '.prt'
  end subroutine 

  subroutine fem_mesh_distribution_write_files ( dir_path, prefix, nparts, parts )
    implicit none
    ! Parameters
    character (len=*), intent(in)   :: dir_path
    character (len=*), intent(in)   :: prefix
    integer(ip)  , intent(in)       :: nparts
    type(fem_mesh_distribution), intent(in)  :: parts(nparts)
    ! Locals
    integer (ip)                     :: i
    character(len=:), allocatable    :: name, rename ! Deferred-length allocatable character arrays
    
    integer :: lunio
    
    call fem_mesh_distribution_compose_name ( prefix, name )
    
    do i=1,nparts
       rename=name
       call numbered_filename_compose(i,nparts,rename)
       lunio = io_open (trim(dir_path) // '/' // trim(rename))
       call fem_mesh_distribution_write ( lunio, parts(i) )
       call io_close (lunio)
    end do

    ! name, and rename should be automatically deallocated by the compiler when they
    ! go out of scope. Should we deallocate them explicitly for safety reasons?
  end subroutine  fem_mesh_distribution_write_files

  subroutine fem_mesh_distribution_read_files ( dir_path, prefix, nparts, parts )
    implicit none
    ! Parameters 
    character (*), intent(in)    :: dir_path 
    character (*), intent(in)    :: prefix
    integer(ip)  , intent(in)    :: nparts
    type(fem_mesh_distribution), intent(inout)  :: parts(nparts)

    ! Locals 
    integer (ip)                        :: i
    character(len=:), allocatable       :: name,rename ! Deferred-length allocatable character arrays
    integer(ip)                         :: lunio

    call fem_mesh_distribution_compose_name ( prefix, name )

    do i=1,nparts
       rename=name
       call numbered_filename_compose(i,nparts,rename)
       lunio = io_open (trim(dir_path) // '/' // trim(rename))
       call fem_mesh_distribution_read ( lunio, parts(i) )
       call io_close (lunio)
    end do

    ! name, and rename should be automatically deallocated by the compiler when they
    ! go out of scope. Should we deallocate them explicitly for safety reasons?
  end subroutine  fem_mesh_distribution_read_files

end module fem_mesh_distribution_names
