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
!***********************************************************************
! All allocatable arrays
!***********************************************************************
# define var_attr allocatable, target
# define point(a,b) call move_alloc(a,b)
# define generic_status_test             allocated
# define generic_memalloc_interface      memalloc
# define generic_memrealloc_interface    memrealloc
# define generic_memfree_interface       memfree
# define generic_memmovealloc_interface  memmovealloc
# include "debug.i90"
!***********************************************************************
module interpolation_tools_names
  use types
  use memor
#ifdef memcheck
  use iso_c_binding
#endif
  use fem_space_types
  use interpolation_names
  use array_names
  implicit none
  private

  type interpolator_pointer
     type(array_rp2), pointer :: p => NULL()
  end type interpolator_pointer

  ! Types
  public :: interpolator_pointer

# define var_type type(interpolator_pointer)
# define var_size 8
# define bound_kind ip
# include "mem_header.i90"

  public :: memalloc,  memrealloc,  memfree, memmovealloc

  ! Functions
  public :: interpolator_create, interpolator_free

contains

  !==================================================================================================
# include "mem_body.i90"

  !==================================================================================================
    subroutine interpolator_create(gtype,utype,ndime,g_ord,u_ord,gnode,unode,int_array,khie)
    implicit none
    integer(ip)          , intent(in)  :: gtype, utype
    integer(ip)          , intent(in)  :: ndime, g_ord, u_ord, gnode, unode
    type(array_rp2)      , intent(out) :: int_array
    logical(lg), optional, intent(in)  :: khie
    ! Locals
    integer(ip)           :: nlocs, i
    type(interpolation)   :: inter
    real(rp), allocatable :: ref_coord(:,:),auxpo(:)

    call memalloc(ndime,unode,ref_coord,__FILE__,__LINE__)
    if(utype==P_type_id) then    
       call P_refcoord(ref_coord,ndime,u_ord,unode)
    elseif (utype == Q_type_id) then
       call Q_refcoord(ref_coord,ndime,u_ord,unode)
    end if
    
    ! Construct interpolation
    if (utype == P_type_id) then
       call interpolation_create(1,1,1,ndime,gnode,unode,inter)
       call interpolation_local(ref_coord,inter)
    elseif (utype == Q_type_id) then
       nlocs=nint(real(unode,rp)**(1.0_rp/real(ndime,rp)))
       call memalloc(nlocs,auxpo,__FILE__,__LINE__)
       do i=1,nlocs
          auxpo(i) = ref_coord(1,i)
       end do
       call interpolation_create(1,1,1,ndime,gnode,unode,inter)
       call interpolation_local(auxpo,inter,ndime,g_ord,nlocs)
       call memfree(auxpo,__FILE__,__LINE__)
    end if

    ! Store shape
    call array_create(gnode,unode,int_array)
    int_array%a = inter%shape

    ! Deallocate
    call interpolation_free(inter)
    call memfree(ref_coord,__FILE__,__LINE__)

  end subroutine interpolator_create

  !==================================================================================================
  subroutine interpolator_free(int_array)
    implicit none
    type(array_rp2), intent(inout) :: int_array

    call array_free(int_array)

  end subroutine interpolator_free

end module interpolation_tools_names
