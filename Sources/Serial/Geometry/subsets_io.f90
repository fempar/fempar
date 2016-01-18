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
module subsets_io_names
  use types_names
  use stdio_names
  use memor_names
  use subsets_names
  implicit none
  private

  ! Functions
  public :: subsets_read_file, subsets_read,     &
            subsets_compose_name, subsets_write, &
            subsets_write_file, subsets_write_files

contains

  !=============================================================================
  ! TODO:
  !  * use stdio functions instead of line to read from files (thus allowing
  !    the posibility of using comments, checking errors, etc.)
  !  * implement other formats (when needed)
  !
  !=============================================================================
  subroutine subsets_read_file(lunio,subset)
    !------------------------------------------------------------------------
    !
    ! This routine reads subsets
    !
    !------------------------------------------------------------------------
    implicit none
    integer(ip)    , intent(in)            :: lunio 
    type(subsets_t), intent(out)           :: subset
    integer(ip)                            :: i,nelem,ielem
    character(1024)                        :: tel

    ! Look for starting
    read(lunio,'(a)') tel
    do while(tel(1:5).ne.'eleme')
       read(lunio,'(a)') tel
    end do

    ! Read first lines to get nelem
    nelem = 0
    do while(tel(1:9).ne.'end eleme')
        read(lunio,'(a)') tel
	       nelem = nelem + 1
    end do
    nelem = nelem - 1

    call subsets_create(nelem,subset)

    ! Look for starting nodes subsets
  !  read(lunio,'(a)') tel
  !  do while(tel(1:5).ne.'eleme')
  !    read(lunio,'(a)') tel
  !  end do
    call io_rewind(lunio)
    read(lunio,'(a)') tel         ! Read "elements"
    read(lunio,'(a)') tel         ! Read 1st element listes
    do while(tel(1:9).ne.'end eleme')
       read(tel,*) ielem,subset%list(ielem)
       read(lunio,'(a)') tel
    end do

    return

  end subroutine subsets_read_file


  !=============================================================================
  subroutine subsets_write_file(lunio,subset)
    !------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------
    implicit none
    ! Parameters
    integer(ip)    , intent(in) :: lunio
    type(subsets_t), intent(in) :: subset
    
    ! Locals
    integer(ip)                 :: ielem
    
    write(lunio,1) 'elements'
    do ielem=1,subset%nelem
       write(lunio,2) ielem,subset%list(ielem)
    end do
    write(lunio,1) 'end elements'
    
1   format(a)
2   format(i10,i10)
    
  end subroutine subsets_write_file

  !=============================================================================
  subroutine subsets_write ( dir_path, prefix, f_subset )
    implicit none 
    ! Parameters
    character (*)  , intent(in)    :: dir_path
    character (*)  , intent(in)    :: prefix
    type(subsets_t), intent(inout) :: f_subset
    
    ! Locals
    character(len=:), allocatable  :: name
    integer(ip) :: lunio

    ! Read subsets
    call subsets_compose_name ( prefix, name )
    lunio = io_open( trim(dir_path)//'/'//trim(name), 'write' )
    call subsets_write_file(lunio, f_subset)
    call io_close(lunio)
    
  end subroutine subsets_write

  !=============================================================================
  subroutine subsets_compose_name ( prefix, name ) 
    implicit none
    character (len=*)             , intent(in)    :: prefix 
    character (len=:), allocatable, intent(inout) :: name
    name = trim(prefix) // '.sbs'
  end subroutine subsets_compose_name

  !=============================================================================
  subroutine subsets_write_files ( dir_path, prefix, nparts, lsubset )
    implicit none
    ! Parameters 
    character(*)    , intent(in)    :: dir_path 
    character(*)    , intent(in)    :: prefix
    integer(ip)     , intent(in)    :: nparts
    type(subsets_t) , intent(in)    :: lsubset (nparts)
    character(len=:), allocatable   :: name, rename ! Deferred-length allocatable character arrays

    ! Locals 
    integer (ip)                    :: i, lunio

    call subsets_compose_name ( prefix, name )

    do i=1,nparts
       rename = name
       call numbered_filename_compose(i,nparts,rename)
       lunio = io_open(trim(dir_path)//'/'//trim(rename),'write' )
       call subsets_write_file(lunio,lsubset(i))
       call io_close(lunio)
    end do

  end subroutine subsets_write_files

  !=============================================================================
  subroutine subsets_read ( dir_path, prefix, f_subset )
    implicit none 
    ! Parameters
    character (*)    , intent(in)  :: dir_path
    character (*)    , intent(in)  :: prefix
    type(subsets_t)  , intent(out) :: f_subset
    
    ! Locals
    character(len=:), allocatable  :: name
    integer(ip)                    :: lunio

    ! Read subsets
    call subsets_compose_name ( prefix, name )
    lunio = io_open( trim(dir_path)//'/'//trim(name), 'read', status='old' )
    call subsets_read_file(lunio, f_subset)
    call io_close(lunio)
    
  end subroutine subsets_read

end module subsets_io_names
