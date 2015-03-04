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
module fem_materials_io
  use types
  use stdio
  use memor
  use fem_materials_names
  implicit none
  private

  ! Functions
  public :: fem_materials_read, fem_materials_compose_name, fem_materials_write
  public :: fem_materials_write_files

contains

  !=============================================================================
  ! TODO:
  !  * use stdio functions instead of line to read from files (thus allowing
  !    the posibility of using comments, checking errors, etc.)
  !  * implement other formats (when needed)
  !
  !=============================================================================
  subroutine fem_materials_read(lunio,mat)
    !------------------------------------------------------------------------
    !
    ! This routine reads materials
    !
    !------------------------------------------------------------------------
    implicit none
    integer(ip)         , intent(in)            :: lunio 
    type(fem_materials) , intent(out)           :: mat
    integer(ip)                                 :: i,nelem,ielem
    character(1024)                             :: tel

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

    call fem_materials_create(nelem,mat)

    ! Look for starting nodes materials
  !  read(lunio,'(a)') tel
  !  do while(tel(1:5).ne.'eleme')
  !    read(lunio,'(a)') tel
  !  end do
    call io_rewind(lunio)
    read(lunio,'(a)') tel         ! Read "elements"
    read(lunio,'(a)') tel         ! Read 1st element listes
    do while(tel(1:9).ne.'end eleme')
       read(tel,*) ielem,mat%list(ielem)
       read(lunio,'(a)') tel
    end do

    return

  end subroutine fem_materials_read


!=============================================================================
  subroutine fem_materials_write(lunio,mater)
    !------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------
    implicit none
    integer(ip)         , intent(in) :: lunio
    type(fem_materials), intent(in) :: mater
    character(80) :: fmt
    integer(ip)   :: ielem

    !fmt='(2x,i10)','//trim(ch(nodes%ncode))//'(1x,i2),'//trim(ch(nodes%nvalu))//'(1x,e14.7))'
    !write(*,*) fmt

    write(lunio,1) 'elements'
    do ielem=1,mater%nelem
       write(lunio,2) ielem,mater%list(ielem)
    end do
    write(lunio,1) 'end elements'

1   format(a)
2   format(i10,i10)

  end subroutine fem_materials_write

 !=============================================================================
  subroutine fem_materials_compose_name ( prefix, name ) 
    implicit none
    character *(*), intent(in)        :: prefix 
    character *(*), intent(out)       :: name
    name = trim(prefix) // '.mat'
  end subroutine fem_materials_compose_name



!=============================================================================
  subroutine fem_materials_write_files ( dir_path, prefix, nparts, lmater )
    implicit none
    ! Parameters 
    character *(*), intent(in)       :: dir_path 
    character *(*), intent(in)       :: prefix
    integer(ip)   , intent(in)       :: nparts
    type(fem_materials), intent(in) :: lmater (nparts)
    character(256)                   :: name, rename

    ! Locals 
    integer (ip)                     :: i, lunio

    call fem_materials_compose_name ( prefix, name )

    do i=1,nparts
       rename = name
       call numbered_filename_compose(i,nparts,rename)
       lunio = io_open(trim(dir_path)//'/'//trim(rename),'write' )
       call fem_materials_write(lunio,lmater(i))
       call io_close(lunio)
    end do

  end subroutine fem_materials_write_files

end module fem_materials_io
