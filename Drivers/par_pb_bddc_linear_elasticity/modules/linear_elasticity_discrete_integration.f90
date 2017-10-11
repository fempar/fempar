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
module linear_elasticity_discrete_integration_names
  use fempar_names
  use linear_elasticity_constitutive_models_names
  use linear_elasticity_analytical_functions_names
  implicit none
# include "debug.i90"
  private
  real(rp), parameter, public :: nu_heterogeneous = 0.2_rp

  type, extends(discrete_integration_t), abstract :: linear_elasticity_discrete_integration_t
  private
  integer(ip) :: number_dimensions
  integer(ip) :: number_fields
  integer(ip) :: number_components
  real(rp)    :: mass_coefficient
  real(rp)    :: residual_coefficient
  real(rp)    :: currrent_time
  character(len=:)  , allocatable :: terms_to_integrate
  character(len=256), allocatable :: fe_type(:)
  character(len=256), allocatable :: field_type(:)
  character(len=256), allocatable :: field_name(:)
  integer(ip), allocatable :: field_blocks(:)
  logical    , allocatable :: field_coupling(:,:)
  class(vector_function_t), pointer :: source_term  => null()
  ! The solution needs to be stored here to avoid an aliasing problem
  ! in nonlinear_operator_apply
  class(fe_function_t), pointer :: fe_function => null()
  class(fe_function_t), pointer :: fe_function_old => null()
  type(linear_elasticity_analytical_functions_t), pointer :: analytical_functions => null()
  real(rp) :: current_time
contains
  procedure :: get_number_fields
  procedure :: get_number_components
  procedure :: get_field_blocks
  procedure :: get_field_coupling
  procedure :: get_field_type
  procedure :: get_fe_type
  procedure :: get_field_name
  procedure :: get_fe_function
  procedure :: set_number_fields
  procedure :: set_number_components
  procedure :: set_field_blocks
  procedure :: set_field_coupling
  procedure :: set_analytical_functions
  procedure :: set_fe_function
  procedure :: set_old_fe_function
  procedure :: set_terms_to_integrate
  procedure :: set_mass_coefficient
  procedure :: set_residual_coefficient
  procedure :: set_current_time
  procedure(create_interface)          , deferred :: create
  procedure(is_symmetric_interface)    , deferred :: is_symmetric
  procedure(is_coercive_interface)     , deferred :: is_coercive
  procedure :: free => linear_elasticity_discrete_integration_free
end type linear_elasticity_discrete_integration_t

abstract interface
  subroutine  create_interface(this,number_dimensions,analytical_functions)
    import :: linear_elasticity_discrete_integration_t, linear_elasticity_analytical_functions_t, ip
    implicit none
    class(linear_elasticity_discrete_integration_t)       , intent(inout) :: this
    integer(ip)                                   , intent(in)    :: number_dimensions
    type(linear_elasticity_analytical_functions_t), target, intent(in)    :: analytical_functions
  end subroutine create_interface
  function is_symmetric_interface(this)
    import :: linear_elasticity_discrete_integration_t
    implicit none
    class(linear_elasticity_discrete_integration_t)       , intent(inout) :: this
    logical :: is_symmetric_interface
  end function is_symmetric_interface
  function is_coercive_interface(this)
    import :: linear_elasticity_discrete_integration_t
    implicit none
    class(linear_elasticity_discrete_integration_t)       , intent(inout) :: this
    logical :: is_coercive_interface
  end function is_coercive_interface
  subroutine init_solution_interface(this,fe_space,solution)
    import :: linear_elasticity_discrete_integration_t, serial_fe_space_t, fe_function_t
    class(linear_elasticity_discrete_integration_t)       , intent(inout) :: this
    class(serial_fe_space_t)                      , intent(in)    :: fe_space
    type(fe_function_t), target                   , intent(inout) :: solution
  end subroutine init_solution_interface
  !subroutine check_solution_interface(this,fe_space)
  !  import :: linear_elasticity_discrete_integration_t, serial_fe_space_t
  !  class(linear_elasticity_discrete_integration_t)       , intent(inout) :: this
  !  class(serial_fe_space_t)                      , intent(in)    :: fe_space
  !end subroutine check_solution_interface
end interface

character(*), parameter :: tangent_terms                 = 'tangent'
character(*), parameter :: translation_terms             = 'translation'
character(*), parameter :: tangent_and_translation_terms = 'both'
public :: tangent_terms, translation_terms, tangent_and_translation_terms

character(*), parameter :: implicit_terms  = 'implicit'
character(*), parameter :: explicit_terms  = 'explicit'
character(*), parameter :: transient_terms = 'transient' ! + implicit
public :: implicit_terms, explicit_terms, transient_terms

type, extends(linear_elasticity_discrete_integration_t) :: irreducible_discrete_integration_t
private
contains
procedure :: create             => irreducible_discrete_integration_create
procedure :: integrate_galerkin => irreducible_discrete_integration_integrate
procedure :: is_symmetric       => irreducible_discrete_integration_is_symmetric
procedure :: is_coercive        => irreducible_discrete_integration_is_coercive
!procedure :: init_solution => irreducible_discrete_integration_init_solution
end type irreducible_discrete_integration_t

type, extends(irreducible_discrete_integration_t) :: irreducible_heterogeneous_discrete_integration_t
private
real(rp), public    :: elastic_modulus
contains
  !procedure :: create             => irreducible_discrete_integration_create
procedure :: integrate_galerkin => irreducible_heterogeneous_discrete_integration_integrate
!procedure :: is_symmetric       => irreducible_discrete_integration_is_symmetric
!procedure :: is_coercive        => irreducible_discrete_integration_is_coercive
!procedure :: init_solution => irreducible_discrete_integration_init_solution 
end type irreducible_heterogeneous_discrete_integration_t  

character(*), parameter :: discrete_integration_type_irreducible = 'irreducible'
character(*), parameter :: discrete_integration_type_mixed_u_p   = 'mixed_u_p'

public :: linear_elasticity_discrete_integration_t, irreducible_discrete_integration_t, irreducible_heterogeneous_discrete_integration_t
public :: discrete_integration_type_irreducible

!public :: linear_elasticity_discrete_integration_t, irreducible_discrete_integration_t, mixed_u_p_discrete_integration_t
!public :: discrete_integration_type_irreducible, discrete_integration_type_mixed_u_p

contains

subroutine linear_elasticity_discrete_integration_free ( this )
implicit none
class(linear_elasticity_discrete_integration_t)   ,intent(inout)  :: this
integer(ip) :: istat
deallocate(this%fe_type,stat=istat)  ; check(istat==0)
deallocate(this%field_type,stat=istat); check(istat==0)
deallocate(this%field_name,stat=istat); check(istat==0)
call memfree(this%field_blocks,__FILE__,__LINE__)
call memfree(this%field_coupling,__FILE__,__LINE__)     
end subroutine linear_elasticity_discrete_integration_free

subroutine set_terms_to_integrate(this,terms_to_integrate)
implicit none
class(linear_elasticity_discrete_integration_t), intent(inout)  :: this
character(*)                           , intent(in)     :: terms_to_integrate
this%terms_to_integrate = terms_to_integrate
end subroutine set_terms_to_integrate

subroutine set_analytical_functions ( this, analytical_functions )
implicit none
class(linear_elasticity_discrete_integration_t)   ,intent(inout)  :: this
type(linear_elasticity_analytical_functions_t), target, intent(in)    :: analytical_functions
this%analytical_functions => analytical_functions
end subroutine set_analytical_functions

subroutine set_fe_function (this, fe_function)
implicit none
class(linear_elasticity_discrete_integration_t), intent(inout) :: this
type(fe_function_t) , target           , intent(in)    :: fe_function
this%fe_function => fe_function
end subroutine set_fe_function

subroutine set_old_fe_function (this, old_fe_function)
implicit none
class(linear_elasticity_discrete_integration_t), intent(inout) :: this
class(fe_function_t),            target, intent(in)    :: old_fe_function
this%fe_function_old => old_fe_function
end subroutine set_old_fe_function

subroutine set_mass_coefficient(this,mass_coefficient)
implicit none
class(linear_elasticity_discrete_integration_t), intent(inout) :: this
real(rp)                               , intent(in)    :: mass_coefficient
this%mass_coefficient = mass_coefficient
end subroutine set_mass_coefficient

subroutine set_residual_coefficient(this, residual_coefficient)
implicit none
class(linear_elasticity_discrete_integration_t), intent(inout) :: this
real(rp)                               , intent(in)    :: residual_coefficient
this%residual_coefficient = residual_coefficient
end subroutine set_residual_coefficient

function get_fe_function (this)
implicit none
class(linear_elasticity_discrete_integration_t), target, intent(in) :: this
class(fe_function_t), pointer                               :: get_fe_function
get_fe_function => this%fe_function
end function get_fe_function

!=============================================================================

subroutine set_current_time(this,current_time)
implicit none
class(linear_elasticity_discrete_integration_t), intent(inout) :: this
real(rp)                               , intent(in)    :: current_time
this%current_time = current_time
end subroutine set_current_time

!==============================================================================

function get_number_fields(this)
implicit none
class(linear_elasticity_discrete_integration_t), intent(in)   :: this
integer(ip) :: get_number_fields
get_number_fields = this%number_fields
end function get_number_fields

subroutine set_number_fields(this,number_fields)
implicit none
class(linear_elasticity_discrete_integration_t), intent(inout) :: this
integer(ip)                        , intent(in)    :: number_fields
this%number_fields = number_fields
end subroutine  set_number_fields

function get_number_components(this)
implicit none
class(linear_elasticity_discrete_integration_t), intent(in)    :: this
integer(ip) :: get_number_components
get_number_components = this%number_components
end function get_number_components

subroutine set_number_components(this, number_components)
implicit none
class(linear_elasticity_discrete_integration_t), intent(inout) :: this
integer(ip)                        , intent(in)    :: number_components
this%number_components = number_components
end subroutine set_number_components

function get_field_blocks(this)
implicit none
class(linear_elasticity_discrete_integration_t), target, intent(in) :: this
integer(ip), pointer :: get_field_blocks(:)
get_field_blocks => this%field_blocks
end function get_field_blocks

subroutine set_field_blocks(this,field_blocks)
implicit none
class(linear_elasticity_discrete_integration_t), intent(inout) :: this
integer(ip)                            , intent(in)    :: field_blocks(:)
if(allocated(this%field_blocks)) call memfree(this%field_blocks, __FILE__,__LINE__)
call memalloc(this%number_fields,this%field_blocks, __FILE__,__LINE__)
this%field_blocks = field_blocks
end subroutine set_field_blocks

function get_field_coupling(this)
implicit none
class(linear_elasticity_discrete_integration_t), target, intent(in) :: this
logical , pointer :: get_field_coupling(:,:)
get_field_coupling => this%field_coupling
end function get_field_coupling

subroutine  set_field_coupling(this,field_coupling)
implicit none
class(linear_elasticity_discrete_integration_t), intent(inout) :: this
logical                                , intent(in)    :: field_coupling(:,:)
if(allocated(this%field_coupling)) call memfree(this%field_coupling, __FILE__,__LINE__)
call memalloc(this%number_fields,this%number_fields,this%field_coupling, __FILE__,__LINE__)
this%field_coupling = field_coupling
end subroutine set_field_coupling

function get_field_type(this,field_id)
implicit none
class(linear_elasticity_discrete_integration_t), intent(inout) :: this
integer(ip)                            , intent(in)    :: field_id
character(len=:), allocatable :: get_field_type
get_field_type = trim(this%field_type(field_id))
end function get_field_type

function get_fe_type(this,field_id)
implicit none
class(linear_elasticity_discrete_integration_t), intent(inout) :: this
integer(ip)                            , intent(in)    :: field_id
character(len=:), allocatable :: get_fe_type
get_fe_type = trim(this%fe_type(field_id))
end function get_fe_type

function get_field_name(this,field_id)
implicit none
class(linear_elasticity_discrete_integration_t), intent(inout) :: this
integer(ip)                            , intent(in)    :: field_id
character(len=:), allocatable :: get_field_name
get_field_name = trim(this%field_name(field_id))
end function get_field_name

#include "sbm_irreducible_discrete_integration.i90"
#include "sbm_irreducible_heterogeneous_discrete_integration.i90"
!#include "sbm_mixed_u_p_discrete_integration.i90"

end module linear_elasticity_discrete_integration_names
