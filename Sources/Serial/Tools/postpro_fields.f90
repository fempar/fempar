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

module postpro_fields_names

  use types_names
  use memor_names
  use fe_space_names
  
  implicit none

  ! Elemental postprocess field type
  type fe_postprocess_field_t
     integer(ip)           :: nnode
!     integer(ip)           :: nvars
     real(rp), allocatable :: nodal_properties(:,:)
  end type fe_postprocess_field_t

  ! Postprocess field type
  type postprocess_field_t
     logical                       :: filled = .false.
     integer(ip)                   :: nelem                 ! Number of elements
     integer(ip)                   :: nvars                 ! Number of variables 
     integer(ip)                   :: iunk                  ! Unknown index of fe_space
     integer(ip),      allocatable :: interpolation_list(:) ! ?
     character(len=:), allocatable :: name                  ! Name of the postprocess field
     type(fe_space_t), pointer     :: fe_space => NULL()    ! Points to fe_space_t
     type(fe_postprocess_field_t), allocatable :: fe_postprocess_field(:)
   contains
     procedure    :: create => postprocess_field_create
     procedure    :: free   => postprocess_field_free
  end type postprocess_field_t
  
  ! Types
  public :: postprocess_field_t
  
contains
  
  !==================================================================================================
  subroutine postprocess_field_create(postprocess_field,name,nvars,nelem,iunk,fe_space)
    implicit none
    class(postprocess_field_t), intent(inout)  :: postprocess_field
    integer(ip),                intent(in)     :: nvars, nelem, iunk
    character(len=:), allocatable, intent(in)  :: name
    type(fe_space_t),  target,  intent(in)     :: fe_space
    ! Locals
    integer(ip)      :: nnode, ielem
    
    postprocess_field%nvars = nvars
    postprocess_field%nelem = nelem
    postprocess_field%name  = name
    postprocess_field%iunk  = iunk
    postprocess_field%fe_space => fe_space
    call memalloc(nvars, postprocess_field%interpolation_list, __FILE__, __LINE__, valin=0)
    allocate(postprocess_field%fe_postprocess_field(1:nelem))
    do ielem=1,nelem
       nnode = fe_space%finite_elements(ielem)%reference_element_vars(iunk)%p%nnode
       call memalloc(nnode, nvars, postprocess_field%fe_postprocess_field(ielem)%nodal_properties,&
            &        __FILE__, __LINE__, valin=0.0_rp)
       postprocess_field%fe_postprocess_field%nnode = nnode
    end do
    
  end subroutine postprocess_field_create

  !==================================================================================================
  subroutine postprocess_field_free(postprocess_field)
    implicit none
    class(postprocess_field_t), intent(inout) :: postprocess_field
    ! Locals
    integer(ip)     :: ielem
    
    do ielem=1,postprocess_field%nelem
       call memfree(postprocess_field%fe_postprocess_field(ielem)%nodal_properties, __FILE__,     &
            &       __LINE__)
    end do

    postprocess_field%fe_space => null()
    deallocate(postprocess_field%fe_postprocess_field)
    call memfree(postprocess_field%interpolation_list)

  end subroutine postprocess_field_free

end module postpro_fields_names
