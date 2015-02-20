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
module template_element_class
  use types
  use migratory_element_class

  implicit none
  private

  type, extends(migratory_element) ::template_element
     integer(igp) :: globalID
     real(rp)     :: area 
   contains
     procedure :: size   => template_element_size
     procedure :: pack   => template_element_pack
     procedure :: unpack => template_element_unpack
     procedure :: setglobalID => template_element_setglobalid
     procedure :: setarea     => template_element_setarea
  end type template_element

  public :: template_element

contains
  subroutine template_element_size (my, n)
    implicit none
    class(template_element), intent(in)  :: my
    integer(ip)            , intent(out) :: n
    
    ! Locals
    integer(ip) :: mold(1)
    integer(ip) :: size_of_igp, size_of_rp

    size_of_igp = size(transfer(1_igp,mold))
    size_of_rp  = size(transfer(1_rp ,mold))

    n = size_of_igp + size_of_rp

  end subroutine template_element_size

  subroutine template_element_pack (my, n, buffer)
    implicit none
    class(template_element), intent(in)  :: my
    integer(ip)            , intent(in)  :: n
    integer(ip)            , intent(out) :: buffer(n)

    ! Locals
    integer(ip) :: mold(1)
    integer(ip) :: size_of_igp, size_of_rp
    integer(ip) :: start, end

    size_of_igp = size(transfer(1_igp, mold))
    size_of_rp  = size(transfer(1_rp , mold))

    start = 1
    end   = start + size_of_igp -1

    buffer(start:end) = transfer(my%globalID,mold)

    start = end + 1
    end   = start + size_of_rp - 1 

    buffer(start:end) = transfer(my%area,mold)
    
  end subroutine template_element_pack

  subroutine template_element_unpack(my, n, buffer)
    implicit none
    class(template_element), intent(inout) :: my
    integer(ip)            , intent(in)    :: n
    integer(ip)            , intent(in)    :: buffer(n)

    ! Locals
    integer(ip) :: mold(1)
    integer(ip) :: size_of_igp, size_of_rp
    integer(ip) :: start, end

    size_of_igp = size(transfer(1_igp,mold))
    size_of_rp  = size(transfer(1_rp ,mold))

    start = 1
    end   = start + size_of_igp -1

    my%globalID = transfer(buffer(start:end), my%globalID)
    
    start = end + 1
    end   = start + size_of_rp - 1
    
    my%area = transfer(buffer(start:end), my%area)

  end subroutine template_element_unpack

  subroutine template_element_setglobalid (my, globalid)
    implicit none
    class(template_element), intent(inout)  :: my
    integer(igp)           , intent(in)     :: globalid
    
    my%globalID = globalid
    
  end subroutine template_element_setglobalid
  
  subroutine template_element_setarea (my, area)
    implicit none
    class(template_element), intent(inout) :: my
    real(rp)               , intent(in)    :: area
    
    my%area = area
    
  end subroutine template_element_setarea

end module template_element_class
