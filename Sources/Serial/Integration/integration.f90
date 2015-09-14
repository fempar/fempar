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
# include "debug.i90"
module integration_names
  use types_names
  use assembly_names
  use integrable_names
  use problem_names
  use integration_tools_names
  use femap_interp_names
  use finite_element_names
  use fe_space_names
  use assembly_names
  use block_matrix_names
  use matrix_names
  use block_vector_names
  use serial_scalar_array_names
  use finite_element_names
  use dof_descriptor_names
  implicit none
  private

  ! Abstract assembly interface
  abstract interface
     subroutine assembly_interface(finite_element, dof_descriptor, a)
       import :: finite_element_t, dof_descriptor_t, integrable_t
       implicit none
       type(finite_element_t), intent(in)    :: finite_element
       type(dof_descriptor_t), intent(in)    :: dof_descriptor
       class(integrable_t)   , intent(inout) :: a
     end subroutine assembly_interface
     subroutine assembly_face_interface(fe_face, finite_element, dof_descriptor, a)
       import :: fe_face_t, finite_element_pointer_t, dof_descriptor_t, integrable_t
       implicit none
       type(fe_face_t)               , intent(in)    :: fe_face
       type(finite_element_pointer_t), intent(in)    :: finite_element(2)
       type(dof_descriptor_t), intent(in)    :: dof_descriptor
       class(integrable_t)   , intent(inout) :: a
     end subroutine assembly_face_interface
  end interface

  public :: volume_integral

contains

  subroutine volume_integral(approx,fe_space,res1,res2,alternative_assembly,alternative_assembly_face)
    implicit none
    ! Parameters
    type(fe_space_t)                    , intent(inout) :: fe_space
    class(integrable_t)                 , intent(inout) :: res1
    class(integrable_t), optional       , intent(inout) :: res2
    type(discrete_integration_pointer_t), intent(inout) :: approx(:)
    procedure(assembly_interface)      , optional :: alternative_assembly
    procedure(assembly_face_interface) , optional :: alternative_assembly_face

    ! Locals
    integer(ip) :: i,ielem,iface,ivar,nvars
    type(finite_element_pointer_t) :: finite_elements(2)

    do i=1,size(approx)

       if(approx(i)%p%domain_dimension==3) then

          do ielem=1,fe_space%g_trian%num_elems
             if(associated(approx(i)%p%domain)) then
                if(approx(i)%p%domain(ielem)==0) cycle
             end if

             nvars = fe_space%finite_elements(ielem)%num_vars
             ! Compute integration tools on ielem for each ivar (they all share the quadrature inside integ)
             do ivar=1,nvars
                call volume_integrator_update(fe_space%finite_elements(ielem)%integ(ivar)%p,fe_space%g_trian%elems(ielem)%coordinates)
             end do

             call approx(i)%p%compute(fe_space%finite_elements(ielem))

             ! Assembly first contribution
             if(present(alternative_assembly)) then
                call alternative_assembly(fe_space%finite_elements(ielem),fe_space%dof_descriptor,res1) 
             else
                call assembly(fe_space%finite_elements(ielem),fe_space%dof_descriptor,res1) 
             end if

             if(present(res2)) then
                if(present(alternative_assembly)) then
                   call alternative_assembly(fe_space%finite_elements(ielem),fe_space%dof_descriptor,res2) 
                else
                   call assembly(fe_space%finite_elements(ielem),fe_space%dof_descriptor,res2)
                end if
             end if

          end do

       else if(approx(i)%p%domain_dimension==2) then

          check(associated(approx(i)%p%domain))
          do iface=1,size(approx(i)%p%domain)

             ! TO DO: develop face_integrator_update in the same line of volume_integrator_update
             ! adapting the old code commented in integration_tools.f90 line 352 (old subroutine
             ! integ_faces).
             !
             ! nvars = fe_space%fe_faces(iface)%num_vars
             ! ! Compute integration tools on iface for each ivar (they all share the quadrature inside integ)
             ! do ivar=1,nvars
             !    call volume_integrator_update(fe_space%fe_faces(iface)%integ(ivar)%p,fe_space%g_trian%elems(ielem)%coordinates)
             ! end do

             call approx(i)%p%compute_face(fe_space%fe_faces(iface))

             finite_elements(1)%p => fe_space%finite_elements(fe_space%fe_faces(iface)%neighbor_element(1))
             finite_elements(2)%p => fe_space%finite_elements(fe_space%fe_faces(iface)%neighbor_element(2))

             ! Assembly first contribution
             if(present(alternative_assembly)) then
                call alternative_assembly_face(fe_space%fe_faces(iface), finite_elements, fe_space%dof_descriptor,res1) 
             else
                call assembly_face(fe_space%fe_faces(iface), finite_elements, fe_space%dof_descriptor,res1) 
             end if

             if(present(res2)) then
                if(present(alternative_assembly)) then
                   call alternative_assembly_face(fe_space%fe_faces(iface), finite_elements, fe_space%dof_descriptor,res2) 
                else
                   call assembly_face(fe_space%fe_faces(iface), finite_elements, fe_space%dof_descriptor,res2)
                end if
             end if

          end do

       end if
    end do
  end subroutine volume_integral

end module integration_names
