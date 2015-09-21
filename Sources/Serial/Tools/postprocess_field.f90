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

module postprocess_field_names

  use types_names
  use memor_names
  use serial_fe_space_names
  use abstract_environment_names
  implicit none
# include "debug.i90"
  private

  integer(ip), parameter :: use_max_interpolation_order = 0
  integer(ip), parameter :: use_interpolation_order_from_variable = 1

  ! Elemental postprocess field type
  type fe_postprocess_field_t
     integer(ip)           :: nnode
     real(rp), allocatable :: nodal_properties(:,:)
   contains
     procedure    :: create => fe_postprocess_field_create
     procedure    :: free   => fe_postprocess_field_free
  end type fe_postprocess_field_t

  ! Postprocess field type
  type, abstract :: postprocess_field_t
     logical                       :: filled = .false.
     integer(ip)                   :: nvars                 ! Number of variables 
     character(len=:), allocatable :: name                  ! Name of the postprocess field
     type(serial_fe_space_t), pointer     :: fe_space => NULL()  ! Points to fe_space_t
     class(abstract_environment_t), pointer :: env => NULL() ! Points to environment
     type(fe_postprocess_field_t), allocatable :: fe_postprocess_field(:)
   contains
     generic      :: create       => postprocess_field_create
     procedure, private :: postprocess_field_create
     procedure    :: free         => postprocess_field_free
     procedure    :: is_finalized => postprocess_field_is_finalized
     procedure    :: compute_and_finalize_field => postprocess_field_compute_and_finalize_field
     procedure(postprocess_field_compute_interface), deferred :: compute_field
  end type postprocess_field_t

  abstract interface
     subroutine postprocess_field_compute_interface(postprocess_field)
       import :: postprocess_field_t
       implicit none
       class(postprocess_field_t), intent(inout)  :: postprocess_field
     end subroutine postprocess_field_compute_interface
  end interface

  ! Types
  public :: postprocess_field_t, use_max_interpolation_order, use_interpolation_order_from_variable
  
contains
  
  !==================================================================================================
  subroutine postprocess_field_create(postprocess_field,name,nvars,fe_space,env,                    &
       &                              interpolation_order_mode,variable_identifier)
    implicit none
    class(postprocess_field_t), intent(inout)  :: postprocess_field
    integer(ip),                intent(in)     :: nvars    ! Number of components of the postprocess_field
    character(len=*),           intent(in)     :: name     ! Name of the postprocess field
    type(serial_fe_space_t),  target,  intent(in)     :: fe_space 
    integer(ip),                intent(in)     :: interpolation_order_mode
    integer(ip),   optional,    intent(in)     :: variable_identifier
    class(abstract_environment_t), target, intent(in) :: env
    ! Locals
    integer(ip)      :: ndime, ielem

    assert((interpolation_order_mode.eq.use_max_interpolation_order).or.(interpolation_order_mode.eq.use_interpolation_order_from_variable))
    if (interpolation_order_mode.eq.use_interpolation_order_from_variable) then 
       check(present(variable_identifier))
    end if

    postprocess_field%nvars = nvars
    postprocess_field%name  = name
    postprocess_field%fe_space => fe_space
    postprocess_field%env => env
    
    if (postprocess_field%env%am_I_fine_task()) then
       allocate(postprocess_field%fe_postprocess_field(1:fe_space%g_trian%num_elems))
       
       do ielem=1,fe_space%g_trian%num_elems
          call postprocess_field%fe_postprocess_field(ielem)%create(postprocess_field%fe_space,ielem,&
               &                                   nvars, interpolation_order_mode, variable_identifier)
       end do
    end if
    
  end subroutine postprocess_field_create

  !==================================================================================================
  subroutine postprocess_field_free(postprocess_field)
    implicit none
    class(postprocess_field_t), intent(inout) :: postprocess_field
    ! Locals
    integer(ip)     :: ielem, nelem
    
    if (postprocess_field%env%am_I_fine_task()) then
       nelem = postprocess_field%fe_space%g_trian%num_elems
       do ielem=1,nelem
          call postprocess_field%fe_postprocess_field(ielem)%free
       end do
       deallocate(postprocess_field%fe_postprocess_field)
    end if

    postprocess_field%fe_space => null()

  end subroutine postprocess_field_free

  !==================================================================================================
  subroutine fe_postprocess_field_create(fe_postprocess_field, fe_space, ielem, nvars,              &
       &                                 interpolation_order_mode, variable_identifier)
    implicit none
    class(fe_postprocess_field_t), intent(inout) :: fe_postprocess_field
    type(serial_fe_space_t),              intent(inout) :: fe_space
    integer(ip),                   intent(in)    :: ielem, nvars, interpolation_order_mode
    integer(ip),   optional,       intent(in)    :: variable_identifier
    ! Locals
    integer(ip)   :: i, nnode, istat
    
    ! Get the maximum interpolation order among all varibles in ielem
    nnode = -1
    do i=1,fe_space%finite_elements(ielem)%num_vars
       if (nnode.lt.fe_space%finite_elements(ielem)%reference_element_vars(i)%p%nnode) then
          nnode = fe_space%finite_elements(ielem)%reference_element_vars(i)%p%nnode
       end if
    end do
    
    select case(interpolation_order_mode)
    case(use_max_interpolation_order)
       fe_postprocess_field%nnode = nnode
    case(use_interpolation_order_from_variable)
       fe_postprocess_field%nnode = fe_space%finite_elements_info(variable_identifier)%nnode
    end select

    call memalloc(nnode, nvars, fe_postprocess_field%nodal_properties, __FILE__, __LINE__,          &
         &        valin=0.0_rp)
    
  end subroutine fe_postprocess_field_create

  !==================================================================================================
  subroutine fe_postprocess_field_free(fe_postprocess_field)
    implicit none
    class(fe_postprocess_field_t), intent(inout) :: fe_postprocess_field

    call memfree(fe_postprocess_field%nodal_properties, __FILE__, __LINE__)
  end subroutine fe_postprocess_field_free

  !==================================================================================================
  subroutine postprocess_field_compute_and_finalize_field(postprocess_field)
    implicit none
    class(postprocess_field_t), intent(inout) :: postprocess_field
    
    if (postprocess_field%env%am_I_fine_task()) then
       call postprocess_field%compute_field()
       postprocess_field%filled = .true.
    end if
  end subroutine postprocess_field_compute_and_finalize_field

  !==================================================================================================
  function postprocess_field_is_finalized(postprocess_field) result(status)
    implicit none
    class(postprocess_field_t), intent(inout) :: postprocess_field
    logical                                   :: status

    status = postprocess_field%filled
  end function postprocess_field_is_finalized

end module postprocess_field_names
