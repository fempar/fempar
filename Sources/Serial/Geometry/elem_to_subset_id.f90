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
module elem_to_subset_id_names
  use types_names
  use memor_names
  use stdio_names
  implicit none
# include "debug.i90"
  
  private

  type elem_to_subset_id_t
     private
     integer(ip)                :: nelem=0                          
     integer(ip), allocatable   :: elem_to_subset_id(:)
   contains
     procedure, non_overridable :: create        => elem_to_subset_id_create
     procedure, non_overridable :: free          => elem_to_subset_id_free
     procedure, non_overridable :: get_subset_id => elem_to_subset_id_get_subset_id 
     procedure, nopass          :: compose_name  => elem_to_subset_id_compose_name
     procedure, non_overridable :: read          => elem_to_subset_id_read
     procedure, non_overridable :: read_file     => elem_to_subset_id_read_file
     procedure, non_overridable :: write         => elem_to_subset_id_read
     procedure, non_overridable :: write_file    => elem_to_subset_id_write_file
  end type elem_to_subset_id_t

  public :: elem_to_subset_id_t

contains

  !===============================================================================================
  subroutine elem_to_subset_id_create( this, nelem )
     implicit none
     class(elem_to_subset_id_t) , intent(inout) :: this
     integer(ip)                , intent(in)    :: nelem
     call this%free()
     this%nelem=nelem
     call memalloc (this%nelem,this%elem_to_subset_id, __FILE__,__LINE__)
     this%elem_to_subset_id=0 
  end subroutine elem_to_subset_id_create

  !===============================================================================================
  subroutine elem_to_subset_id_free( this )
     implicit none
     class(elem_to_subset_id_t), intent(inout) :: this
     this%nelem= 0
     if (allocated(this%elem_to_subset_id)) call memfree (this%elem_to_subset_id,__FILE__,__LINE__)
  end subroutine elem_to_subset_id_free
  
  !===============================================================================================
  function elem_to_subset_id_get_subset_id ( this, ielem ) result ( elem_to_subset_id )
     implicit none
     class(elem_to_subset_id_t), intent(in) :: this
     integer(ip)               , intent(in) :: ielem  
     integer(ip)                            :: elem_to_subset_id
     elem_to_subset_id = this%elem_to_subset_id(ielem)
  end function elem_to_subset_id_get_subset_id
  
  !===============================================================================================
  subroutine elem_to_subset_id_compose_name ( prefix, name ) 
     implicit none
     character (len=*)             , intent(in)    :: prefix 
     character (len=:), allocatable, intent(inout) :: name
     name = trim(prefix) // '.sbs'
  end subroutine elem_to_subset_id_compose_name
  
  !===============================================================================================
  subroutine elem_to_subset_id_read( this, dir_path, prefix )
     implicit none
     ! Parameters
     class(elem_to_subset_id_t) , intent(inout) :: this
     character (*)              , intent(in)    :: dir_path
     character (*)              , intent(in)    :: prefix
     ! Locals
     character(len=:), allocatable              :: name
     integer(ip)                                :: lunio
        
     call this%compose_name( prefix, name )
     lunio = io_open( trim(dir_path)//'/'//trim(name), 'read', status='old' )
     call this%read_file(lunio)
     call io_close(lunio)
  end subroutine elem_to_subset_id_read
  
  !===============================================================================================
  subroutine elem_to_subset_id_read_file( this, lunio )
     implicit none
     ! Parameters
     class(elem_to_subset_id_t) , intent(inout) :: this
     integer(ip)                                :: lunio
     ! Locals
     integer(ip)                                :: i,ielem,nelem
     character(1024)                            :: tel
     
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
     
     call this%create(nelem)
     
     call io_rewind(lunio)
     read(lunio,'(a)') tel   ! Read "elements"
     read(lunio,'(a)') tel   ! Read 1st element lists
     do while(tel(1:9).ne.'end eleme')
        read(tel,*) ielem,this%elem_to_subset_id(ielem)
        read(lunio,'(a)') tel
     end do
     
  end subroutine elem_to_subset_id_read_file
  
  !===============================================================================================
  subroutine elem_to_subset_id_write( this, dir_path, prefix )
     implicit none
     ! Parameters
     class(elem_to_subset_id_t) , intent(inout) :: this
     character (*)              , intent(in)    :: dir_path
     character (*)              , intent(in)    :: prefix
     ! Locals
     character(len=:), allocatable              :: name
     integer(ip)                                :: lunio
     
     call this%compose_name( prefix, name )
     lunio = io_open( trim(dir_path)//'/'//trim(name), 'write' )
     call this%write_file(lunio)
     call io_close(lunio)
  end subroutine elem_to_subset_id_write

  !=============================================================================
  subroutine elem_to_subset_id_write_file(this,lunio)
     implicit none
     ! Parameters
     class(elem_to_subset_id_t) , intent(inout) :: this
     integer(ip)                , intent(in)    :: lunio
     ! Locals
     integer(ip)                                :: ielem
     
     write(lunio,1) 'elements'
     do ielem=1,this%nelem
        write(lunio,2) ielem,this%elem_to_subset_id(ielem)
     end do
     write(lunio,1) 'end elements'
    
1    format(a)
2    format(i10,i10)

  end subroutine elem_to_subset_id_write_file
  
end module elem_to_subset_id_names
