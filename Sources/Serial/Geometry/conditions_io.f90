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
module conditions_io_names
use types_names
use stdio_names
use memor_names
  use conditions_names
  implicit none
  private

  ! Functions
  public :: conditions_read_file, conditions_read, conditions_write, &
       &    conditions_write_file, conditions_compose_name, & 
       &    conditions_write_files

contains

  !=============================================================================
  ! TODO:
  !  * use stdio functions instead of line to read from files (thus allowing
  !    the posibility of using comments, checking errors, etc.)
  !  * implement other formats (when needed)
  !
  !=============================================================================
  subroutine conditions_read_file(lunio,npoin,nodes)
    !------------------------------------------------------------------------
    !
    ! This routine reads conditions
    !
    !------------------------------------------------------------------------
    implicit none
    integer(ip)         , intent(in)            :: lunio,npoin
    type(conditions_t), intent(out)           :: nodes
    integer(ip)                                 :: i,ncode,nvalu,icode,ivalu,ipoin,iboun
    character(1024)                             :: tel

    ! Look for starting
    read(lunio,'(a)') tel
    do while(tel(1:5).ne.'CONDI')
       read(lunio,'(a)') tel
    end do

    ! Read first lines to get ncode and nvalu
    read(lunio,'(a)') tel
    do i=1,75
       if(tel(i:i+4)=='NCODE') then
          ! read(tel(i+9:i+10),*) ncode
          read(tel(i+5:80),*) ncode
       end if
    end do

    read(lunio,'(a)') tel
    do i=1,75
       if(tel(i:i+4)=='NVALU') then
          read(tel(i+5:80),*) nvalu
       end if
    end do

    call conditions_create(ncode,nvalu,npoin,nodes)

    ! Look for starting nodes conditions
    read(lunio,'(a)') tel
    do while(tel(1:5).ne.'NODES')
       read(lunio,'(a)') tel
    end do

    read(lunio,'(a)') tel
    do while(tel(1:9).ne.'END NODES')
       read(tel,*) ipoin,(nodes%code(icode,ipoin),icode=1,ncode), &
            &            (nodes%valu(ivalu,ipoin),ivalu=1,nvalu)
       read(lunio,'(a)') tel
    end do

  end subroutine conditions_read_file

  !=============================================================================
  subroutine conditions_write_file(lunio, nodes, bouns)
    !------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------
    implicit none
    integer(ip)                   , intent(in) :: lunio
    type(conditions_t)          , intent(in) :: nodes
    type(conditions_t), optional, intent(in) :: bouns
    character(80) :: fmt
    integer(ip)   :: ipoin,iboun,icode,ivalu

    fmt='(i10,'//trim(ch(nodes%ncode))//'(1x,i2),'//trim(ch(nodes%nvalu))//'(1x,e14.7))'

    write(lunio,1) 'CONDITIONS'
    write(lunio,2) 'NCODE', nodes%ncode
    write(lunio,2) 'NVALU', nodes%nvalu
    write(lunio,1) 'NODES'
    do ipoin=1,nodes%ncond
       write(lunio,fmt) ipoin,(nodes%code(icode,ipoin),icode=1,nodes%ncode), &
            &               (nodes%valu(ivalu,ipoin),ivalu=1,nodes%nvalu)
    end do
    write(lunio,1) 'END NODES'

1   format(a)
2   format(a5,i10)

  end subroutine conditions_write_file

 !=============================================================================
  subroutine conditions_compose_name ( prefix, name ) 
    implicit none
    character (len=*)             , intent(in)  :: prefix 
    character (len=:), allocatable, intent(out) :: name
    name = trim(prefix) // '.cnd'
  end subroutine conditions_compose_name

 !=============================================================================
  subroutine conditions_write_files ( dir_path, prefix, nparts, lnodes, lbouns )
    implicit none
    ! Parameters 
    character *(*), intent(in)       :: dir_path 
    character *(*), intent(in)       :: prefix
    integer(ip)   , intent(in)       :: nparts
    type(conditions_t), intent(in) :: lnodes (nparts)
    type(conditions_t), optional, intent(in) :: lbouns (nparts)
    character(len=:), allocatable    :: name, rename

    ! Locals 
    integer (ip)                     :: i, lunio

    call conditions_compose_name ( prefix, name )

    do i=1,nparts
       rename = name
       call numbered_filename_compose(i,nparts,rename)
       lunio = io_open(trim(dir_path)//'/'//trim(rename),'write' )
       if(present(lbouns)) then
          call conditions_write_file(lunio,lnodes(i),lbouns(i))
       else
          call conditions_write_file(lunio,lnodes(i))
       end if
       call io_close(lunio)
    end do

  end subroutine conditions_write_files

  !=============================================================================
  subroutine conditions_read ( dir_path, prefix, npoin, f_conditions )
    implicit none 
    ! Parameters
    character (*)       , intent(in)  :: dir_path
    character (*)       , intent(in)  :: prefix
    integer(ip)         , intent(in)  :: npoin
    type(conditions_t), intent(out) :: f_conditions
    
    ! Locals
    character(len=:), allocatable  :: name
    integer(ip) :: lunio

    ! Read conditions
    call conditions_compose_name ( prefix, name )
    lunio = io_open( trim(dir_path)//'/'//trim(name), 'read', status='old' )
    call conditions_read_file(lunio, npoin, f_conditions)
    call io_close(lunio)
    
  end subroutine conditions_read

  !=============================================================================
  subroutine conditions_write ( dir_path, prefix, f_conditions )
    implicit none 
    ! Parameters
    character (*)       , intent(in) :: dir_path
    character (*)       , intent(in) :: prefix
    type(conditions_t), intent(in) :: f_conditions
    
    ! Locals
    character(len=:), allocatable  :: name
    integer(ip) :: lunio

    ! Read conditions
    call conditions_compose_name ( prefix, name )
    lunio = io_open( trim(dir_path)//'/'//trim(name), 'read', status='old' )
    call conditions_write_file(lunio, f_conditions)
    call io_close(lunio)
    
  end subroutine conditions_write

end module conditions_io_names
