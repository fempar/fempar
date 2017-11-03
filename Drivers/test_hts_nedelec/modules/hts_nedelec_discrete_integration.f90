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
module hts_nedelec_discrete_integration_names
  use fempar_names
  use hts_nedelec_params_names
  use hts_theta_method_names 
  
  implicit none
# include "debug.i90"
  private
  
  integer(ip), parameter :: hts = 1
  integer(ip), parameter :: air = 2
  
  type, extends(discrete_integration_t) :: hts_nedelec_discrete_integration_t
     private
     ! Physical parameters 
     real(rp)                          :: air_resistivity
     real(rp)                          :: hts_resistivity
     real(rp)                          :: air_permeability
     real(rp)                          :: hts_permeability 
     real(rp)                          :: critical_electric_field 
     real(rp)                          :: critical_current
     real(rp)                          :: nonlinear_exponent 
     character(len=:), allocatable     :: integration_type  
     
     ! Solution to integrate nonlinear, time dependent 
     type(theta_method_t)    , pointer :: theta_method 
     type(fe_function_t)     , pointer :: H_current         
     type(fe_function_t)     , pointer :: H_previous        
     class(vector_function_t), pointer :: source_term        => NULL()
   contains
     procedure :: create                => hts_nedelec_discrete_integration_create
     procedure :: integrate_galerkin    => hts_nedelec_discrete_integration_integrate 
     procedure :: set_integration_type 
     procedure :: compute_resistivity
     procedure :: compute_tangent_resistivity 
  end type hts_nedelec_discrete_integration_t
  
  public :: hts_nedelec_discrete_integration_t
  
contains
   
  subroutine hts_nedelec_discrete_integration_create(this, theta_method, H_current, H_previous, hts_nedelec_params, source_term) 
  implicit none 
  class(hts_nedelec_discrete_integration_t), intent(inout) :: this
  type(theta_method_t)    ,target          , intent(in)    :: theta_method
  type(fe_function_t)     ,target          , intent(in)    :: H_current 
  type(fe_function_t)     ,target          , intent(in)    :: H_previous 
  type(hts_nedelec_params_t)               , intent(in)    :: hts_nedelec_params 
  class(vector_function_t), target         , intent(in)    :: source_term
    
  this%air_resistivity         = hts_nedelec_params%get_air_resistivity() 
  this%hts_resistivity         = hts_nedelec_params%get_hts_resistivity() 
  this%air_permeability        = hts_nedelec_params%get_air_permeability() 
  this%hts_permeability        = hts_nedelec_params%get_hts_permeability()
  this%critical_electric_field = hts_nedelec_params%get_critical_electric_field()
  this%critical_current        = hts_nedelec_params%get_critical_current()
  this%nonlinear_exponent      = hts_nedelec_params%get_nonlinear_exponent()
  
  ! Set solution and source term 
  this%theta_method  => theta_method 
  this%H_previous    => H_previous 
  this%H_current     => H_current
  this%source_term   => source_term 
  end subroutine hts_nedelec_discrete_integration_create
    
  subroutine hts_nedelec_discrete_integration_integrate ( this, fe_space, assembler )
    implicit none
    class(hts_nedelec_discrete_integration_t)   , intent(in)    :: this
    class(serial_fe_space_t)                    , intent(inout) :: fe_space
    class(assembler_t)             , intent(inout) :: assembler

    ! FE space traversal-related data types
    class(fe_cell_iterator_t), allocatable :: fe
    
    ! FE integration-related data types
    type(quadrature_t)       , pointer :: quad
    type(point_t)            , pointer :: quad_coords(:)    
    type(vector_field_t)               :: H_value_previous
    type(vector_field_t) , allocatable :: H_current_curl_values(:)  
    type(fe_cell_function_vector_t)    :: fe_cell_function_previous, fe_cell_function_current 
    type(vector_field_t), allocatable  :: shape_values_H(:,:)
    type(vector_field_t), allocatable  :: curl_values_H(:,:)    
    
    ! Nonlinear parameters computed in the integrate 
    real(rp)   :: resistivity, tangent_resistivity
    real(rp)   :: permeability 
         
    ! FE matrix and vector i.e., A_K + f_K
    real(rp), allocatable              :: elmat(:,:), elvec(:)

    integer(ip)  :: istat
    integer(ip)  :: qpoint, num_quad_points
    integer(ip)  :: idof, jdof, num_dofs
    real(rp)     :: factor, time_factor 
    type(vector_field_t), allocatable :: source_term_values(:,:)
    real(rp), allocatable             :: current_time(:)
    
    assert ( associated(this%source_term) )
    assert ( associated(this%H_current) )
    assert ( associated(this%H_previous) )
    
    call fe_space%set_up_cell_integration()
    call fe_cell_function_previous%create(fe_space, 1) 
    call fe_cell_function_current%create(fe_space,  1) 

    call fe_space%create_fe_cell_iterator(fe)
    num_dofs = fe%get_num_dofs()
    call memalloc ( num_dofs, num_dofs, elmat, __FILE__, __LINE__ )
    call memalloc ( num_dofs, elvec, __FILE__, __LINE__ )
    quad             => fe%get_quadrature()
    num_quad_points  = quad%get_num_quadrature_points()
    
    time_factor = this%theta_method%get_theta() * this%theta_method%get_time_step() 
    allocate (source_term_values(num_quad_points,1), stat=istat); check(istat==0)
    allocate (H_current_curl_values(num_quad_points), stat=istat); check(istat==0) 
    allocate (current_time(1), stat=istat); check(istat==0) 
    current_time = this%theta_method%get_current_time() 
    
    do while ( .not. fe%has_finished())
       
       ! Update FE-integration related data structures
       call fe%update_integration()
       call fe_cell_function_previous%update(fe, this%H_previous)
       call fe_cell_function_current%update(fe, this%H_current)

       ! Get quadrature coordinates to evaluate boundary value
       quad_coords => fe%get_quadrature_points_coordinates()
       
       ! Evaluate pressure source term at quadrature points
       call this%source_term%get_values_set( quad_coords, current_time, source_term_values)
       
       ! Evaluate current curl values  
       call fe_cell_function_current%compute_quadrature_points_curl_values(H_current_curl_values)
       
       ! Compute element matrix and vector
       elmat = 0.0_rp
       elvec = 0.0_rp
             
       call fe%get_values(shape_values_H,1)
       call fe%get_curls(curl_values_H,1)
        
       do qpoint = 1, num_quad_points
          factor = fe%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
          
          ! Compute nonlinear resistivity 
          resistivity  = this%compute_resistivity( H_current_curl_values(qpoint), fe%get_set_id() )
          permeability = this%air_permeability
          
          ! Previous solution to integrate RHS contribution 
          call fe_cell_function_previous%get_value( qpoint, H_value_previous )
          
          !  mu_0 ( H,v ) + rho( curl(H), curl(v) )  
          do idof=1, fe%get_num_dofs_field(1)
            do jdof=1, fe%get_num_dofs_field(1)  
              elmat(idof,jdof) = elmat(idof,jdof) + &
              (permeability/time_factor*shape_values_H(jdof,qpoint)*shape_values_H(idof,qpoint) + resistivity*curl_values_H(jdof,qpoint)*curl_values_H(idof,qpoint))*factor
                  
                  if ( this%integration_type=='add_tangent_terms' .and. fe%get_set_id()== hts ) then 
                   tangent_resistivity = this%compute_tangent_resistivity(H_current_curl_values(qpoint), curl_values_H(jdof,qpoint)) 
                   elmat(idof,jdof) = elmat(idof,jdof) + tangent_resistivity*curl_values_H(idof,qpoint)*H_current_curl_values(qpoint)*factor                    
                  end if       
                  
            end do            
            elvec(idof) = elvec(idof) + (source_term_values(qpoint,1)+permeability/time_factor*H_value_previous)*shape_values_H(idof,qpoint)*factor
          end do
       
      end do ! Qpoint loop 
          
       call fe%assembly( this%H_current, elmat, elvec, assembler )
       call fe%next()
    end do
    call fe_space%free_fe_cell_iterator(fe)

    deallocate (source_term_values, stat=istat); check(istat==0)
    deallocate (H_current_curl_values, stat=istat); check(istat==0)
    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )
  end subroutine hts_nedelec_discrete_integration_integrate
   
  ! -----------------------------------------------------------------------------------------------
  function compute_resistivity(this, curl_H, material) result(resistivity)
  implicit none 
  class(hts_nedelec_discrete_integration_t)   , intent(in)       :: this
  type(vector_field_t)                        , intent(inout)    :: curl_H 
  integer                                     , intent(in)       :: material
  real(rp)                                                       :: resistivity 
  ! Locals
  real(rp) :: Ec, Jc, n 
  
  Ec = this%critical_electric_field 
  Jc = this%critical_current 
  n  = this%nonlinear_exponent 

  if ( material == hts ) then ! HTS DOMAIN: Nonlinear resistivity = Ec/Jc*|| curl(H) / Jc ||**n 
     if (this%nonlinear_exponent .ge. 1 ) then 
        resistivity = Ec/Jc*(curl_H%nrm2()/Jc)**n + 1e-16_rp
     else
        resistivity = this%hts_resistivity
     end if
  else if ( material == air ) then ! Air domain 
     resistivity = this%air_resistivity 
  else 
     assert(.false.) 
  end if

  end function compute_resistivity
  
  ! -----------------------------------------------------------------------------------------------
  function compute_tangent_resistivity(this, curl_H, curl_shape) result(tangent_resistivity) 
  implicit none 
  class(hts_nedelec_discrete_integration_t)   , intent(in)       :: this
  type(vector_field_t)                        , intent(inout)    :: curl_H 
  type(vector_field_t)                        , intent(in)       :: curl_shape 
  real(rp)                                                       :: tangent_resistivity 
   ! Locals
  real(rp) :: Ec, Jc, n 
  
  Ec = this%critical_electric_field 
  Jc = this%critical_current 
  n  = this%nonlinear_exponent 
  
  ! Tangent Nonlinear resistivity = Ec/Jc * n *|| curl(H) / Jc ||**(n-2) < curl(H)/Jc | curl(Phi)/Jc >
  if (this%nonlinear_exponent .ge. 1 ) then 
      tangent_resistivity = n*Ec/Jc*( curl_H%nrm2()/Jc )**(n-2.0_rp) * ((1.0_rp/Jc)*curl_H) * ((1.0_rp/Jc)*curl_shape) 
  else
      tangent_resistivity = 0.0_rp 
  end if 
  
  end function compute_tangent_resistivity
  
  ! -----------------------------------------------------------------------------------------------
  subroutine set_integration_type(this, integration_type) 
  implicit none 
  class(hts_nedelec_discrete_integration_t), intent(inout) :: this
  character(len=*)                         , intent(in)    :: integration_type                             
  
  assert( (integration_type == 'regular') .or. (integration_type=='add_tangent_terms') ) 
  this%integration_type = integration_type 
  
  end subroutine 
  
  
end module hts_nedelec_discrete_integration_names
