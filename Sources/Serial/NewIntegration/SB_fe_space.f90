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
module SB_fe_space_names
  use memor_names
  use allocatable_array_names
  use matrix_names
  use array_names
  use serial_scalar_matrix_names
  use serial_scalar_array_names
  use serial_block_matrix_names
  use serial_block_array_names
  use SB_matrix_array_assembler_names
  use SB_serial_scalar_matrix_array_assembler_names
  use SB_serial_block_matrix_array_assembler_names
  use types_names
  use reference_fe_names
  use triangulation_names
  use reference_fe_factory_names
  use integration_tools_names
  use migratory_element_names
  use conditions_names
  use hash_table_names
  use graph_names
  use sort_names
  implicit none
# include "debug.i90"
  private

  ! This module includes the physical FE and FE space structure. It includes the
  ! following types:
  !
  ! * finite_element_t: It includes pointers to the unknown/geometrical reference_fe_t
  !   and the corresponding volume integrator, the local-to-global mapping for DOFs in
  !   elem2dof and the boundary conditions restricted to the element
  ! * fe_space_t: The abstract fe_space_t, which only includes some deferred methods 
  !   that must be provided by any fe_space_t concretization
  ! * serial_fe_space_t: The serial version of the fe_space_t, which includes an array
  !   of finite_element_t, the array of all possible reference_fe_t and volume_integrator_t
  !   in the FE space, the total number of DOFs and an array that provides the DOFs
  !   in a VEF.

  type, abstract :: SB_finite_element_t
     private   
     type(elem_topology_t), pointer :: cell 
     class(reference_fe_t), pointer :: geometry_reference_fe  
     type(SB_quadrature_t), pointer :: quadrature
     type(fe_map_t), pointer :: fe_map
     type(SB_p_volume_integrator_t), pointer :: volume_integrator(:)
   contains
     procedure :: get_volume_integrator
     procedure :: get_quadrature	
     procedure :: get_fe_map

     procedure(fe_get_elem2dof_interface), deferred :: get_elem2dof 
     procedure(fe_get_bc_code_interface), deferred :: get_bc_code 
     procedure(fe_get_bc_value_interface), deferred :: get_bc_value 
     procedure(fe_get_number_nodes_interface), deferred :: get_number_nodes 
     procedure(fe_get_number_nodes_field_interface), deferred :: get_number_nodes_field 
  end type SB_finite_element_t

  abstract interface
     subroutine fe_get_number_nodes_field_interface ( this, number_nodes, n )
       import :: SB_finite_element_t, ip
       implicit none
       class(SB_finite_element_t), intent(in) :: this
       integer(ip), intent(in)                :: n
       integer(ip), intent(inout)             :: number_nodes(:)
     end subroutine fe_get_number_nodes_field_interface
  end interface

  abstract interface
     function fe_get_number_nodes_interface ( this )
       import :: SB_finite_element_t, ip
       implicit none
       class(SB_finite_element_t), intent(in) :: this
       integer(ip)                 :: fe_get_number_nodes_interface
     end function fe_get_number_nodes_interface
  end interface

  abstract interface
     function fe_get_volume_integrator_interface ( this )
       import :: SB_finite_element_t, SB_p_volume_integrator_t
       implicit none
       class(SB_finite_element_t), target, intent(in)  :: this
       type(SB_p_volume_integrator_t), pointer :: fe_get_volume_integrator_interface(:)
     end function fe_get_volume_integrator_interface
  end interface
  abstract interface
     function fe_get_elem2dof_interface( this )
       import :: SB_finite_element_t, i1p_t
       implicit none
       class(SB_finite_element_t), target, intent(in) :: this
       type(i1p_t), pointer :: fe_get_elem2dof_interface(:)
     end function fe_get_elem2dof_interface
  end interface
  abstract interface
     function fe_get_bc_code_interface( this )
       import :: SB_finite_element_t, i1p_t
       implicit none
       class(SB_finite_element_t), target, intent(in) :: this
       type(i1p_t), pointer :: fe_get_bc_code_interface(:)
     end function fe_get_bc_code_interface
  end interface
  abstract interface
     function fe_get_bc_value_interface( this )
       import :: SB_finite_element_t, r1p_t
       implicit none
       class(SB_finite_element_t), target, intent(in) :: this
       type(r1p_t), pointer :: fe_get_bc_value_interface(:)
     end function fe_get_bc_value_interface
  end interface
  abstract interface
     subroutine fe_get_quadrature_interface ( this, quadrature, n )
       import :: SB_finite_element_t, ip, SB_quadrature_t
       implicit none
       class(SB_finite_element_t), target, intent(in)  :: this
       integer(ip), intent(in) :: n
       type(SB_quadrature_t), pointer, intent(inout) :: quadrature
     end subroutine fe_get_quadrature_interface
  end interface
  abstract interface
     subroutine fe_get_fe_map_interface ( this, fe_map, n )
       import :: SB_finite_element_t, ip, fe_map_t
       implicit none
       class(SB_finite_element_t), target, intent(in)  :: this
       integer(ip), intent(in) :: n
       type(fe_map_t), pointer, intent(inout) :: fe_map
     end subroutine fe_get_fe_map_interface
  end interface

  type :: SB_p_finite_element_t
     class(SB_finite_element_t), pointer :: p => NULL ()
  end type SB_p_finite_element_t

!   type, extends(SB_finite_element_t) :: SB_simple_finite_element_t 
!   private
!   type(p_reference_fe_t), pointer :: reference_fe(:)
!   !class(reference_fe_t), pointer :: reference_fe
!   !type(SB_volume_integrator_t), pointer :: volume_integrator
!   integer(ip)     , allocatable   :: elem2dof(:)
!   integer(ip), allocatable :: bc_code(:)
!   real(rp), allocatable :: bc_value(:)
! contains
!   procedure :: create => simple_fe_create
!   procedure :: free   => simple_fe_free
!   procedure :: print  => simple_fe_print

!   !procedure :: get_volume_integrator => simple_fe_get_volume_integrator

!   procedure :: get_elem2dof => simple_fe_get_elem2dof
!   procedure :: get_bc_code => simple_fe_get_bc_code
!   procedure :: get_bc_value => simple_fe_get_bc_value
!   procedure :: get_number_nodes => simple_fe_get_number_nodes
!   procedure :: get_number_nodes_field => simple_fe_get_number_nodes_field

!   ! Legacy part not touched
!   !procedure :: size   => simple_finite_element_size
!   !procedure :: pack   => simple_finite_element_pack
!   !procedure :: unpack => simple_finite_element_unpack    
! end type SB_simple_finite_element_t

!type :: SB_p_simple_finite_element_t
!  class(SB_simple_finite_element_t), pointer :: p => NULL ()
!end type SB_p_simple_finite_element_t

public :: SB_finite_element_t!, SB_simple_finite_element_t, SB_p_simple_finite_element_t

type, extends(SB_finite_element_t) :: SB_composite_finite_element_t
private
!type(SB_p_simple_finite_element_t), allocatable :: field_fe(:) 
integer(ip) :: number_nodes
type(i1p_t), pointer :: elem2dof(:)
type(i1p_t), pointer :: bc_code(:)
type(r1p_t), pointer :: bc_value(:)
type(p_reference_fe_t), pointer :: reference_fe(:)
contains
  !procedure :: create => composite_fe_create
  !procedure :: free   => composite_fe_free
  !procedure :: print  => composite_fe_print

  !procedure :: get_volume_integrator => composite_fe_get_volume_integrator

procedure :: get_elem2dof => composite_fe_get_elem2dof
procedure :: get_bc_code => composite_fe_get_bc_code
procedure :: get_bc_value => composite_fe_get_bc_value
procedure :: get_number_nodes => composite_fe_get_number_nodes
procedure :: get_number_nodes_field => composite_fe_get_number_nodes_field

! Legacy part not touched
!procedure :: size   => composite_finite_element_size
!procedure :: pack   => composite_finite_element_pack
!procedure :: unpack => composite_finite_element_unpack  
end type SB_composite_finite_element_t

type, abstract :: SB_fe_space_t
type(triangulation_t), pointer :: triangulation
type(p_reference_fe_t), allocatable :: reference_fe_geo_list(:)
type(p_reference_fe_t), allocatable :: reference_fe_phy_list(:)
type(SB_p_quadrature_t), allocatable :: quadrature(:)
type(p_fe_map_t), allocatable :: fe_map(:)! Integrator
type(SB_p_volume_integrator_t), allocatable :: volume_integrator(:)	
contains
procedure (get_fe_interface), deferred :: get_fe
procedure (initialize_volume_integrator_interface), deferred :: initialize_volume_integrator
procedure :: initialize_quadrature
procedure :: initialize_fe_map
procedure (create_assembler_interface)        , deferred :: create_assembler
procedure (symbolic_setup_assembler_interface), deferred :: symbolic_setup_assembler
procedure (get_blocks_interface), deferred :: get_blocks
end type SB_fe_space_t

abstract interface
   ! Gives a pointer of the i (global index) FE of the FE space
function get_fe_interface( this, i )
 import :: SB_fe_space_t, SB_finite_element_t, ip
 implicit none
 class(SB_fe_space_t), target, intent(in) :: this
 integer(ip) :: i
 class(SB_finite_element_t), pointer :: get_fe_interface
end function get_fe_interface

function get_blocks_interface( this )
 import :: SB_fe_space_t, ip
 implicit none
 class(SB_fe_space_t), target, intent(in) :: this
 integer(ip), pointer :: get_blocks_interface(:)
end function get_blocks_interface

! Computes the fe_mapping_t and updates the volume integrator
subroutine initialize_volume_integrator_interface( this )
 import :: SB_fe_space_t, ip
 implicit none
 ! Parameters
 class(SB_fe_space_t), intent(inout) :: this  
end subroutine initialize_volume_integrator_interface

! Selects the dynamic type of class(matrix_array_assembler_t), and that of
! the class(matrix_t) and class(array_t) data members that it contains.
! It also creates these latter two data members. This subroutine follows
! the FACTORY METHOD design pattern.
function create_assembler_interface(this,& 
    diagonal_blocks_symmetric_storage,&
    diagonal_blocks_symmetric,&
    diagonal_blocks_sign)
 import :: SB_fe_space_t, SB_matrix_array_assembler_t, ip
 implicit none
 class(SB_fe_space_t)           , intent(in) :: this
 logical                        , intent(in) :: diagonal_blocks_symmetric_storage(:)
 logical                        , intent(in) :: diagonal_blocks_symmetric(:)
 integer(ip)                    , intent(in) :: diagonal_blocks_sign(:)
 class(SB_matrix_array_assembler_t), pointer :: create_assembler_interface
end function create_assembler_interface

! Set-ups symbolically the class(matrix_t) data member that it contains.
! Essentially, this amounts for computing the coupling among DoFs in the 
! sparsity pattern of class(matrix_t).
subroutine symbolic_setup_assembler_interface(this,matrix_array_assembler)
 import :: SB_fe_space_t, SB_matrix_array_assembler_t
 implicit none
 class(SB_fe_space_t)              , intent(in)    :: this
 class(SB_matrix_array_assembler_t), intent(inout) :: matrix_array_assembler
end subroutine symbolic_setup_assembler_interface
end interface

type :: SB_p_fe_space_t
class(SB_fe_space_t), pointer :: p => NULL ()
end type SB_p_fe_space_t

public :: SB_fe_space_t, SB_p_fe_space_t

! type, extends(SB_fe_space_t) :: SB_serial_fe_space_t
! private

! type(SB_composite_finite_element_t), allocatable :: fe_array(:)
! !type(vector_space_t)  :: fe_function_space


! integer(ip) :: number_dofs     

! ! Acceleration arrays
! type(list_2d_t)       :: vef2dof       ! An auxiliary array to accelerate some parts of the code

! contains

!   ! TBPs
! procedure :: create => serial_fe_space_create
! procedure :: free   => serial_fe_space_free
! procedure :: print  => serial_fe_space_print

! procedure :: get_fe => serial_fe_space_get_fe
! procedure :: initialize_volume_integrator  => serial_fe_space_initialize_volume_integrator
! procedure :: get_max_number_nodes => serial_fe_space_get_max_number_nodes
! procedure :: get_blocks => serial_fe_space_get_blocks

! procedure :: fill_dof_info => serial_fe_space_fill_dof_info
! procedure :: create_assembler => serial_fe_space_create_assembler

! procedure :: symbolic_setup_assembler => serial_fe_space_symbolic_setup_assembler

! ! procedure :: fe_iterator ! Work to be done Alberto / Javier

! end type SB_serial_fe_space_t

! type :: SB_p_serial_fe_space_t
! type(SB_serial_fe_space_t), pointer :: p => NULL ()
! end type SB_p_serial_fe_space_t

! public :: SB_serial_fe_space_t, SB_p_serial_fe_space_t

type, extends(SB_fe_space_t) :: SB_composite_fe_space_t  

type(SB_composite_finite_element_t), allocatable :: fe_array(:)
integer(ip) :: number_fe_spaces   
!class(SB_p_serial_fe_space_t), allocatable :: fe_space_array(:)  
integer(ip), allocatable :: field_blocks(:)
integer(ip) :: number_blocks
logical, allocatable :: fields_coupling(:,:)
integer(ip), allocatable :: number_dofs(:)
! Acceleration arrays
type(list_2d_t), allocatable       :: vef2dof(:)       ! An auxiliary array to accelerate some parts of the code

contains
  ! TBPs
procedure :: get_number_blocks => composite_fe_space_get_number_blocks
procedure :: get_number_fe_spaces => composite_fe_space_get_number_fe_spaces

procedure :: create => composite_fe_space_create
!procedure :: free   => composite_fe_space_free
!procedure :: print  => composite_fe_space_print

procedure :: get_fe => composite_fe_space_get_fe
procedure :: initialize_volume_integrator  => composite_fe_space_initialize_volume_integrator
procedure :: get_max_number_nodes => composite_fe_space_get_max_number_nodes
procedure :: get_blocks => composite_fe_space_get_blocks
procedure :: get_fields_coupling => composite_fe_space_get_fields_coupling
procedure :: fill_dof_info => composite_fe_space_fill_dof_info
procedure :: create_assembler => composite_fe_space_create_assembler

procedure :: symbolic_setup_assembler => composite_fe_space_symbolic_setup_assembler

end type SB_composite_fe_space_t

public :: SB_composite_fe_space_t
contains

  ! Includes with all the TBP and supporting subroutines for the types above.
  ! In a future, we would like to use the submodule features of FORTRAN 2008.

#include "sbm_composite_finite_element.i90"

!#include "sbm_simple_finite_element.i90"

!#include "sbm_serial_fe_space.i90"

#include "sbm_composite_fe_space.i90"

end module SB_fe_space_names
