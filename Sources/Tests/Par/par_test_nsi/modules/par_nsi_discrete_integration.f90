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
module nsi_discrete_integration_names
  use fempar_names
  use par_nsi_constitutive_models_names
  use par_nsi_analytical_functions_names
  implicit none
# include "debug.i90"
  private

  type, extends(discrete_integration_t) :: nsi_discrete_integration_t
     private
     integer(ip) :: number_dimensions
     integer(ip) :: number_fields
     integer(ip) :: number_components
     real(rp)    :: mass_coefficient
     real(rp)    :: residual_coefficient
     real(rp)    :: currrent_time
     character(len=:), allocatable :: terms_to_integrate
     character(len=:), allocatable :: terms_in_the_residual
     character(len=256), allocatable :: fe_type(:)
     character(len=256), allocatable :: field_type(:)
     character(len=256), allocatable :: field_name(:)
     integer(ip), allocatable :: field_blocks(:)
     logical    , allocatable :: field_coupling(:,:)
     class(vector_function_t), pointer :: source_term  => null()
     ! The solution needs to be stored here to avoid an aliasing problem
     ! in nonlinear_operator_apply
     class(fe_function_t), pointer :: fe_function => null()
     !class(fe_function_t), pointer :: fe_function_old => null()
     type(par_nsi_analytical_functions_t), pointer :: analytical_functions => null()
     real(rp) :: current_time
     real(rp) :: viscosity
   contains
     procedure :: create               => nsi_discrete_integration_create
     procedure :: set_evaluation_point => nsi_discrete_integration_set_evaluation_point
     procedure :: integrate_galerkin   => nsi_discrete_integration_integrate  
     procedure :: integrate_residual   => nsi_discrete_integration_integrate_residual  
     procedure :: integrate_tangent    => nsi_discrete_integration_integrate_tangent      
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
     procedure :: set_terms_in_the_residual
     procedure :: set_mass_coefficient
     procedure :: set_residual_coefficient
     procedure :: set_current_time
     procedure :: free => nsi_discrete_integration_free
  end type nsi_discrete_integration_t
  
  character(*), parameter :: tangent_terms                 = 'tangent'
  character(*), parameter :: translation_terms             = 'translation'
  character(*), parameter :: tangent_and_translation_terms = 'both'
  public :: tangent_terms, translation_terms, tangent_and_translation_terms

  character(*), parameter :: implicit_terms  = 'implicit'
  character(*), parameter :: explicit_terms  = 'explicit'
  character(*), parameter :: transient_terms = 'transient' ! + implicit
  
  character(*), parameter :: discrete_integration_type_galerkin = 'galerkin'
  public :: implicit_terms, explicit_terms, transient_terms

  public :: nsi_discrete_integration_t
  public :: discrete_integration_type_galerkin
  
contains

  subroutine nsi_discrete_integration_free ( this )
     implicit none
     class(nsi_discrete_integration_t)   ,intent(inout)  :: this
     integer(ip) :: istat
     deallocate(this%fe_type,stat=istat)  ; check(istat==0)
     deallocate(this%field_type,stat=istat); check(istat==0)
     deallocate(this%field_name,stat=istat); check(istat==0)
     call memfree(this%field_blocks,__FILE__,__LINE__)
     call memfree(this%field_coupling,__FILE__,__LINE__)     
  end subroutine nsi_discrete_integration_free
  
  subroutine nsi_discrete_integration_set_evaluation_point ( this, evaluation_point )
     implicit none
     class(nsi_discrete_integration_t)   ,intent(inout)  :: this
     class(vector_t)                     ,intent(in)     :: evaluation_point
     call this%fe_function%set_free_dof_values(evaluation_point)
  end subroutine nsi_discrete_integration_set_evaluation_point

  subroutine set_terms_to_integrate(this,terms_to_integrate)
     implicit none
     class(nsi_discrete_integration_t), intent(inout)  :: this
     character(*)                           , intent(in)     :: terms_to_integrate
     this%terms_to_integrate = terms_to_integrate
  end subroutine set_terms_to_integrate  

  subroutine set_terms_in_the_residual(this,terms_in_the_residual)
     implicit none
     class(nsi_discrete_integration_t), intent(inout)  :: this
     character(*)                           , intent(in)     :: terms_in_the_residual
     this%terms_in_the_residual = terms_in_the_residual
  end subroutine set_terms_in_the_residual

  subroutine set_analytical_functions ( this, analytical_functions )
     implicit none
     class(nsi_discrete_integration_t)   ,intent(inout)  :: this
     type(par_nsi_analytical_functions_t), target, intent(in)    :: analytical_functions
     this%analytical_functions => analytical_functions
  end subroutine set_analytical_functions
  
  subroutine set_fe_function (this, fe_function)
    implicit none
    class(nsi_discrete_integration_t), intent(inout) :: this
    type(fe_function_t) , target           , intent(in)    :: fe_function
    this%fe_function => fe_function
  end subroutine set_fe_function

  subroutine set_old_fe_function (this, old_fe_function)
    implicit none
    class(nsi_discrete_integration_t), intent(inout) :: this
    class(fe_function_t),            target, intent(in)    :: old_fe_function
    !this%fe_function_old => old_fe_function
  end subroutine set_old_fe_function

  subroutine set_mass_coefficient(this,mass_coefficient)
    implicit none
    class(nsi_discrete_integration_t), intent(inout) :: this
    real(rp)                               , intent(in)    :: mass_coefficient
    this%mass_coefficient = mass_coefficient
  end subroutine set_mass_coefficient

  subroutine set_residual_coefficient(this, residual_coefficient)
    implicit none
    class(nsi_discrete_integration_t), intent(inout) :: this
    real(rp)                               , intent(in)    :: residual_coefficient
    this%residual_coefficient = residual_coefficient
  end subroutine set_residual_coefficient
  
  function get_fe_function (this)
    implicit none
    class(nsi_discrete_integration_t), target, intent(in) :: this
    class(fe_function_t), pointer                               :: get_fe_function
    get_fe_function => this%fe_function
  end function get_fe_function

  !=============================================================================
  
  subroutine set_current_time(this,fe_space,current_time)
    implicit none
    class(nsi_discrete_integration_t), intent(inout) :: this
    class(serial_fe_space_t)         , intent(in)    :: fe_space
    real(rp)                         , intent(in)    :: current_time
    this%current_time = current_time
  end subroutine set_current_time

!==============================================================================
  
  function get_number_fields(this)
    implicit none
    class(nsi_discrete_integration_t), intent(in)   :: this
    integer(ip) :: get_number_fields
    get_number_fields = this%number_fields
  end function get_number_fields

  subroutine set_number_fields(this,number_fields)
    implicit none
    class(nsi_discrete_integration_t), intent(inout) :: this
    integer(ip)                        , intent(in)    :: number_fields
    this%number_fields = number_fields
  end subroutine  set_number_fields

  function get_number_components(this)
    implicit none
    class(nsi_discrete_integration_t), intent(in)    :: this
    integer(ip) :: get_number_components
    get_number_components = this%number_components
  end function get_number_components

  subroutine set_number_components(this, number_components)
    implicit none
    class(nsi_discrete_integration_t), intent(inout) :: this
    integer(ip)                        , intent(in)    :: number_components
    this%number_components = number_components
  end subroutine set_number_components

  function get_field_blocks(this)
    implicit none
    class(nsi_discrete_integration_t), target, intent(in) :: this
    integer(ip), pointer :: get_field_blocks(:)
    get_field_blocks => this%field_blocks
  end function get_field_blocks

  subroutine set_field_blocks(this,field_blocks)
    implicit none
    class(nsi_discrete_integration_t), intent(inout) :: this
    integer(ip)                            , intent(in)    :: field_blocks(:)
    if(allocated(this%field_blocks)) call memfree(this%field_blocks, __FILE__,__LINE__)
    call memalloc(this%number_fields,this%field_blocks, __FILE__,__LINE__)
    this%field_blocks = field_blocks
  end subroutine set_field_blocks

  function get_field_coupling(this)
    implicit none
    class(nsi_discrete_integration_t), target, intent(in) :: this
    logical , pointer :: get_field_coupling(:,:)
    get_field_coupling => this%field_coupling
  end function get_field_coupling

  subroutine  set_field_coupling(this,field_coupling)
    implicit none
    class(nsi_discrete_integration_t), intent(inout) :: this
    logical                                , intent(in)    :: field_coupling(:,:)
    if(allocated(this%field_coupling)) call memfree(this%field_coupling, __FILE__,__LINE__)
    call memalloc(this%number_fields,this%number_fields,this%field_coupling, __FILE__,__LINE__)
    this%field_coupling = field_coupling
  end subroutine set_field_coupling
  
  function get_field_type(this,field_id)
    implicit none
    class(nsi_discrete_integration_t), intent(inout) :: this
    integer(ip)                            , intent(in)    :: field_id
    character(len=:), allocatable :: get_field_type
    get_field_type = trim(this%field_type(field_id))
  end function get_field_type
  
  function get_fe_type(this,field_id)
    implicit none
    class(nsi_discrete_integration_t), intent(inout) :: this
    integer(ip)                            , intent(in)    :: field_id
    character(len=:), allocatable :: get_fe_type
    get_fe_type = trim(this%fe_type(field_id))
  end function get_fe_type
  
  function get_field_name(this,field_id)
    implicit none
    class(nsi_discrete_integration_t), intent(inout) :: this
    integer(ip)                            , intent(in)    :: field_id
    character(len=:), allocatable :: get_field_name
    get_field_name = trim(this%field_name(field_id))
  end function get_field_name

#include "sbm_nsi_discrete_integration.i90"
  
end module nsi_discrete_integration_names
