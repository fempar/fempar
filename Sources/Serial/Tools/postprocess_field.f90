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
  use fe_space_names
  use abstract_environment_names
  use fe_space_types_names, only : max_FE_types, max_ndime
  use problem_names,        only : physical_problem_t
  use hash_table_names,     only : new_index
  implicit none
# include "debug.i90"
  private

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
     integer(ip)                   :: interpolation_order   ! Interpolation order of the field
     character(len=:), allocatable :: name                  ! Name of the postprocess field
     type(fe_space_t), pointer     :: p_fe_space => NULL()  ! Points to fe_space_t
     class(abstract_environment_t), pointer :: p_env => NULL() ! Points to environment
     type(fe_postprocess_field_t), allocatable :: fe_postprocess_field(:)
   contains
     procedure    :: create       => postprocess_field_create
     procedure    :: free         => postprocess_field_free
     procedure    :: is_finalized => postprocess_field_is_finalized
     procedure    :: compute_and_finalize_field => postprocess_field_compute_and_finalize_field
     procedure(postprocess_field_compute_interface), deferred :: compute_field
  end type postprocess_field_t

  abstract interface
     subroutine postprocess_field_compute_interface(postprocess_field,prob,fill_process)
       import :: postprocess_field_t, physical_problem_t, ip
       implicit none
       class(postprocess_field_t), intent(inout)  :: postprocess_field
       class(physical_problem_t),   intent(in)     :: prob
       integer(ip), optional,      intent(in)     :: fill_process
     end subroutine postprocess_field_compute_interface
  end interface

  ! Types
  public :: postprocess_field_t
  
contains
  
  !==================================================================================================
  subroutine postprocess_field_create(postprocess_field,name,nvars,interpolation_order,fe_space,env)
    implicit none
    class(postprocess_field_t), intent(inout)  :: postprocess_field
    integer(ip),                intent(in)     :: nvars, interpolation_order
    character(len=*),           intent(in)     :: name
    type(fe_space_t),  target,  intent(in)     :: fe_space
    class(abstract_environment_t), target, intent(in) :: env
    ! Locals
    integer(ip)      :: ndime, ielem

    postprocess_field%nvars = nvars
    postprocess_field%name  = name
    postprocess_field%interpolation_order = interpolation_order
    postprocess_field%p_fe_space => fe_space
    postprocess_field%p_env => env
    
    if (postprocess_field%p_env%am_I_fine_task()) then
       allocate(postprocess_field%fe_postprocess_field(1:fe_space%g_trian%num_elems))
       
       do ielem=1,fe_space%g_trian%num_elems
          call postprocess_field%fe_postprocess_field(ielem)%create(postprocess_field%p_fe_space,ielem,&
               &                                                    nvars, interpolation_order)
       end do
    end if
    
  end subroutine postprocess_field_create

  !==================================================================================================
  subroutine postprocess_field_free(postprocess_field)
    implicit none
    class(postprocess_field_t), intent(inout) :: postprocess_field
    ! Locals
    integer(ip)     :: ielem, nelem
    
    if (postprocess_field%p_env%am_I_fine_task()) then
       nelem = postprocess_field%p_fe_space%g_trian%num_elems
       do ielem=1,nelem
          call postprocess_field%fe_postprocess_field(ielem)%free
       end do
       deallocate(postprocess_field%fe_postprocess_field)
    end if

    postprocess_field%p_fe_space => null()

  end subroutine postprocess_field_free

  !==================================================================================================
  subroutine fe_postprocess_field_create(fe_postprocess_field, fe_space, ielem, nvars,              &
       &                                 interpolation_order)
    implicit none
    class(fe_postprocess_field_t), intent(inout) :: fe_postprocess_field
    type(fe_space_t),              intent(inout) :: fe_space
    integer(ip),                   intent(in)    :: ielem, nvars, interpolation_order
    ! Locals
    integer(ip)   :: f_type, v_key, ndime, pos_elinf, istat
    
    ndime = fe_space%g_trian%num_dims
    f_type = fe_space%g_trian%elems(ielem)%geo_reference_element%ftype
    v_key = ndime + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)*interpolation_order
    call fe_space%pos_elem_info%get(key=v_key,val=pos_elinf,stat=istat)
    ! Check that the interpolation choosen is already in use
    check( istat /= new_index) 
    fe_postprocess_field%nnode = fe_space%finite_elements_info(pos_elinf)%nnode

    call memalloc(fe_postprocess_field%nnode, nvars, fe_postprocess_field%nodal_properties,         &
         &        __FILE__, __LINE__, valin=0.0_rp)
    
  end subroutine fe_postprocess_field_create

  !==================================================================================================
  subroutine fe_postprocess_field_free(fe_postprocess_field)
    implicit none
    class(fe_postprocess_field_t), intent(inout) :: fe_postprocess_field

    call memfree(fe_postprocess_field%nodal_properties, __FILE__, __LINE__)
  end subroutine fe_postprocess_field_free

  !==================================================================================================
  subroutine postprocess_field_compute_and_finalize_field(postprocess_field,prob,fill_process)
    implicit none
    class(postprocess_field_t), intent(inout) :: postprocess_field
    class(physical_problem_t),  intent(in)    :: prob
    integer(ip), optional,      intent(in)    :: fill_process
    
    if (postprocess_field%p_env%am_I_fine_task()) then
       call postprocess_field%compute_field(prob,fill_process)
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
