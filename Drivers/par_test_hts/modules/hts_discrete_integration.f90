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
module hts_discrete_integration_names
  use fempar_names
  use hts_theta_method_names 
  use par_test_hts_params_names

  implicit none
# include "debug.i90"
  private

  integer(ip), parameter :: air = 0
  integer(ip), parameter :: hts = 1


  type, extends(discrete_integration_t) :: hts_discrete_integration_t
  class(vector_function_t)        , pointer :: source_term        => NULL() 
  type(fe_function_t)             , pointer :: H_current          => NULL()
  type(fe_function_t)             , pointer :: H_previous         => NULL()
  type(theta_method_t)            , pointer :: theta_method       => NULL()    
  character(len=:)            , allocatable :: terms_to_integrate
  character(len=:)            , allocatable :: terms_in_the_residual
  ! Parameter values 
  real(rp)                          :: air_resistivity
  real(rp)                          :: hts_resistivity
  real(rp)                          :: air_permeability
  real(rp)                          :: hts_permeability 
  real(rp)                          :: critical_electric_field 
  real(rp)                          :: critical_current
  real(rp)                          :: nonlinear_exponent
		logical                           :: is_analytical_solution 
  character(len=:), allocatable     :: hts_device_type 
contains
		procedure :: set_evaluation_point => hts_discrete_integration_set_evaluation_point
  procedure :: integrate_galerkin   => hts_discrete_integration_integrate  
  procedure :: integrate_residual   => hts_discrete_integration_integrate_residual  
  procedure :: integrate_tangent    => hts_discrete_integration_integrate_tangent   
		procedure :: set_terms_to_integrate
		procedure :: set_analytical_functions
  procedure :: set_fe_functions 
  procedure :: set_theta_method 
  procedure :: set_parameter_values 
  procedure :: compute_resistivity 
  procedure :: compute_tangent_resistivity
end type hts_discrete_integration_t

public :: hts_discrete_integration_t, air, hts 

contains

  subroutine hts_discrete_integration_set_evaluation_point ( this, evaluation_point )
     implicit none
     class(hts_discrete_integration_t)   ,intent(inout)  :: this
     class(vector_t)                     ,intent(in)     :: evaluation_point
     call this%H_current%set_free_dof_values(evaluation_point)
  end subroutine hts_discrete_integration_set_evaluation_point
		
		  subroutine set_terms_to_integrate(this,terms_to_integrate)
     implicit none
     class(hts_discrete_integration_t)   ,intent(inout)  :: this
     character(*)                        ,intent(in)     :: terms_to_integrate
     this%terms_to_integrate = terms_to_integrate
  end subroutine set_terms_to_integrate  
		
subroutine set_analytical_functions ( this, source_term )
 implicit none
 class(hts_discrete_integration_t)    ,intent(inout)  :: this
 class(vector_function_t), target,        intent(in)     :: source_term
 this%source_term => source_term 
end subroutine set_analytical_functions

subroutine set_fe_functions (this, H_previous, H_current)
 implicit none
 class(hts_discrete_integration_t), intent(inout)       :: this
 type(fe_function_t)  , target, intent(in)  :: H_previous
 type(fe_function_t)  , target, intent(in)  :: H_current 

 this%H_previous => H_previous 
 this%H_current  => H_current 
end subroutine set_fe_functions

subroutine set_theta_method (this, theta_method)
 implicit none
 class(hts_discrete_integration_t), intent(inout)    :: this
 type(theta_method_t) , target       , intent(in)    :: theta_method
 this%theta_method => theta_method 
end subroutine set_theta_method

subroutine set_parameter_values( this, test_params ) 
 implicit none 
 class(hts_discrete_integration_t), intent(inout)       :: this
 type(par_test_hts_params_t)         , intent(in)       :: test_params

 this%air_resistivity         = test_params%get_air_resistivity() 
 this%hts_resistivity         = test_params%get_hts_resistivity() 
 this%air_permeability        = test_params%get_air_permeability() 
 this%hts_permeability        = test_params%get_hts_permeability()
 this%critical_electric_field = test_params%get_critical_electric_field()
 this%critical_current        = test_params%get_critical_current()
 this%nonlinear_exponent      = test_params%get_nonlinear_exponent()
 this%is_analytical_solution  = test_params%get_is_analytical_solution()
 this%hts_device_type         = test_params%get_hts_device_type() 
end subroutine set_parameter_values

subroutine hts_discrete_integration_integrate ( this, fe_space, assembler )
 implicit none
 class(hts_discrete_integration_t)       , intent(in)    :: this
 class(serial_fe_space_t)                , intent(inout) :: fe_space
 class(assembler_t)                      , intent(inout) :: assembler

  call this%integrate_tangent(fe_space, assembler)
  call this%integrate_residual(fe_space, assembler)
end subroutine hts_discrete_integration_integrate

subroutine hts_discrete_integration_integrate_tangent ( this, fe_space, assembler )
 implicit none
 class(hts_discrete_integration_t)       , intent(in)    :: this
 class(serial_fe_space_t)                , intent(inout) :: fe_space
 class(assembler_t)                      , intent(inout) :: assembler

 ! FE space traversal-related data types
 class(fe_cell_iterator_t), allocatable :: fe

 ! FE integration-related data types
 type(quadrature_t)       , pointer :: quad
 type(point_t)            , pointer :: quad_coords(:)
 type(vector_field_t), allocatable  :: shape_curls(:,:)
 type(vector_field_t), allocatable  :: shape_values(:,:)
 type(vector_field_t), pointer      :: H_previous_values(:)
 type(vector_field_t), pointer      :: H_current_values(:)
 type(vector_field_t) , allocatable :: H_current_curl_values(:)
 type(fe_cell_function_vector_t)    :: fe_cell_function_previous
 type(fe_cell_function_vector_t)    :: fe_cell_function_current

 ! FE Jacobian matrix, A_K 
 real(rp), allocatable              :: elmat(:,:)

 integer(ip)  :: istat
 integer(ip)  :: qpoint, num_quad_points
 integer(ip)  :: idof, jdof, num_dofs, max_num_dofs 
 real(rp)     :: factor, time_factor 
 type(vector_field_t), allocatable :: source_term_values(:,:)
 real(rp)     :: current_time(1)

 ! Parameter values 
 real(rp) :: permeability 
 type(tensor_field_t) :: resistivity
 type(tensor_field_t) :: tangent_resistivity

 integer(ip)  :: number_fields

 assert (associated(this%source_term)) 
 assert (associated(this%H_previous))
 assert (associated(this%H_current))
	
 call fe_space%set_up_cell_integration()
 call fe_space%create_fe_cell_iterator(fe)
 call fe_cell_function_previous%create(fe_space, 1)
 call fe_cell_function_current%create(fe_space, 1) 

 max_num_dofs = fe_space%get_max_num_dofs_on_a_cell()
 call memalloc ( max_num_dofs, max_num_dofs, elmat, __FILE__, __LINE__ )
 current_time(1) = this%theta_method%get_current_time()

 quad            => fe%get_quadrature()
 num_quad_points = quad%get_num_quadrature_points()
 allocate (source_term_values(num_quad_points,1), stat=istat); check(istat==0)
 allocate (H_current_curl_values(num_quad_points), stat=istat); check(istat==0)
 time_factor = this%theta_method%get_theta() * this%theta_method%get_time_step() 

 do while ( .not. fe%has_finished())
    if ( fe%is_local() ) then

       ! Update FE-integration related data structures
       call fe%update_integration()
       call fe_cell_function_previous%update(fe, this%H_previous)
       call fe_cell_function_current%update(fe, this%H_current) 

       quad            => fe%get_quadrature()
       num_quad_points =  quad%get_num_quadrature_points()
       num_dofs        =  fe%get_num_dofs()

       ! Get quadrature coordinates to evaluate source_term
       quad_coords => fe%get_quadrature_points_coordinates()

       ! Evaluate pressure source term at quadrature points
       call this%source_term%get_values_set( quad_coords, current_time, source_term_values)

       ! Evaluate values
       H_previous_values => fe_cell_function_previous%get_quadrature_points_values()
       H_current_values  => fe_cell_function_current%get_quadrature_points_values() 
       call fe_cell_function_current%compute_quadrature_points_curl_values(H_current_curl_values)

       ! Compute element matrix and vector
       elmat = 0.0_rp
       call fe%get_curls(shape_curls)
       call fe%get_values(shape_values)
       do qpoint = 1, num_quad_points
          factor = fe%get_det_jacobian(qpoint) * quad%get_weight(qpoint)

          ! Exact resistivity 
          !resistivity = (current_time(1)**2.0_rp *( quad_coords(qpoint)%get(1)**2.0_rp + quad_coords(qpoint)%get(2)**2.0_rp ))**this%nonlinear_exponent
          resistivity  = this%compute_resistivity( H_current_values(qpoint), H_current_curl_values(qpoint), fe%get_set_id() ) 
          permeability = this%air_permeability

          do idof = 1, num_dofs
             do jdof = 1, num_dofs
                ! A_K(i,j) =  (curl(phi_i),curl(phi_j)) + (phi_i,phi_j)
                elmat(idof,jdof) = elmat(idof,jdof) + factor * ( permeability/time_factor * shape_values(jdof,qpoint) * shape_values(idof,qpoint) + & 
                     resistivity * shape_curls(jdof,qpoint) * shape_curls(idof,qpoint))

                if ( fe%get_set_id() > 0 ) then 
                   tangent_resistivity = this%compute_tangent_resistivity( H_current_values(qpoint), H_current_curl_values(qpoint), & 
                                                                           shape_values(jdof,qpoint), shape_curls(jdof,qpoint)) 

                   elmat(idof,jdof) = elmat(idof,jdof) + tangent_resistivity * shape_curls(idof,qpoint) * H_current_curl_values(qpoint)*factor                    
                end if

             end do
          end do
       end do

       ! Apply boundary conditions
       call fe%assembly( elmat, assembler )
						
    end if
    call fe%next()
 end do
 call fe_space%free_fe_cell_iterator(fe)

 deallocate (shape_values, stat=istat);       check(istat==0)
 deallocate (shape_curls, stat=istat);        check(istat==0)
 deallocate (source_term_values, stat=istat); check(istat==0)
 call memfree ( elmat, __FILE__, __LINE__ )
 call fe_cell_function_previous%free() 
 call fe_cell_function_current%free() 
 deallocate ( H_current_curl_values, stat=istat); check(istat==0)
end subroutine hts_discrete_integration_integrate_tangent 

subroutine hts_discrete_integration_integrate_residual ( this, fe_space, assembler )
 implicit none
 class(hts_discrete_integration_t)    , intent(in)    :: this
 class(serial_fe_space_t)             , intent(inout) :: fe_space
 class(assembler_t)                   , intent(inout) :: assembler 

 ! FE space traversal-related data types
 class(fe_cell_iterator_t), allocatable :: fe

 ! FE integration-related data types
 type(quadrature_t)       , pointer :: quad
 type(point_t)            , pointer :: quad_coords(:)
 type(vector_field_t), allocatable  :: shape_curls(:,:)
 type(vector_field_t), allocatable  :: shape_values(:,:)
 type(vector_field_t) , pointer     :: H_current_values(:)
 type(vector_field_t) , pointer     :: H_previous_values(:)
 type(vector_field_t) , allocatable :: H_current_curl_values(:)
 type(fe_cell_function_vector_t)    :: fe_cell_function_previous
 type(fe_cell_function_vector_t)    :: fe_cell_function_current

 ! FE residual i.e., r_K
 integer(ip)          , pointer     :: fe_dofs(:)
 real(rp), allocatable              :: elvec(:) 

 integer(ip)  :: istat
 integer(ip)  :: qpoint, num_quad_points
 integer(ip)  :: idof, jdof, num_dofs, max_num_dofs 
 real(rp)     :: factor, time_factor 
 type(vector_field_t), allocatable :: source_term_values(:,:)
 real(rp)     :: current_time(1)

 ! Parameter values 
 real(rp) :: permeability 
 type(tensor_field_t) :: resistivity 
 
 integer(ip)  :: number_fields

 assert (associated(this%source_term)) 
 assert (associated(this%H_previous))
 assert (associated(this%H_current))
	
 call fe_space%set_up_cell_integration()
 call fe_space%create_fe_cell_iterator(fe)
 call fe_cell_function_previous%create(fe_space, 1)
 call fe_cell_function_current%create(fe_space, 1) 

 max_num_dofs = fe_space%get_max_num_dofs_on_a_cell()
 call memalloc ( max_num_dofs, elvec, __FILE__, __LINE__ )
 current_time(1) = this%theta_method%get_current_time()

 quad            => fe%get_quadrature()
 num_quad_points = quad%get_num_quadrature_points()
 allocate (source_term_values(num_quad_points,1), stat=istat); check(istat==0)
 allocate (H_current_curl_values(num_quad_points), stat=istat); check(istat==0)
 time_factor = this%theta_method%get_theta() * this%theta_method%get_time_step() 

 do while ( .not. fe%has_finished())
    if ( fe%is_local() ) then

       ! Update FE-integration related data structures
       call fe%update_integration()
       call fe_cell_function_previous%update(fe, this%H_previous)
       call fe_cell_function_current%update(fe, this%H_current) 

       quad            => fe%get_quadrature()
       num_quad_points =  quad%get_num_quadrature_points()
       num_dofs        =  fe%get_num_dofs()

       ! Get quadrature coordinates to evaluate source_term
       quad_coords => fe%get_quadrature_points_coordinates()

       ! Evaluate source term at quadrature points
       call this%source_term%get_values_set( quad_coords, current_time, source_term_values)

       ! Evaluate solution values  
       H_previous_values => fe_cell_function_previous%get_quadrature_points_values()
       H_current_values  => fe_cell_function_current%get_quadrature_points_values()
       call fe_cell_function_current%compute_quadrature_points_curl_values(H_current_curl_values)

       ! Compute element vector
       elvec = 0.0_rp
       call fe%get_curls(shape_curls)
       call fe%get_values(shape_values)
       do qpoint = 1, num_quad_points
          factor = fe%get_det_jacobian(qpoint) * quad%get_weight(qpoint)

          resistivity  = this%compute_resistivity( H_current_values(qpoint), H_current_curl_values(qpoint), fe%get_set_id() )           
          permeability = this%air_permeability
										
          do idof = 1, num_dofs									
             ! R_K(i) = A(u)*u - F  
             elvec(idof) = elvec(idof) + factor * permeability/time_factor * shape_values(idof,qpoint) * H_current_values(qpoint)  &
                                       + factor * resistivity * shape_curls(idof, qpoint) * H_current_curl_values(qpoint)          &  
                                       - factor * (source_term_values(qpoint,1) + permeability/time_factor*H_previous_values(qpoint)) * shape_values(idof,qpoint)
          end do
       end do
							
       call fe%assembly(elvec, assembler)
							
    end if
    call fe%next()
 end do
 call fe_space%free_fe_cell_iterator(fe)

 deallocate (shape_values, stat=istat);       check(istat==0)
 deallocate (shape_curls, stat=istat);        check(istat==0)
 deallocate (source_term_values, stat=istat); check(istat==0)
 call memfree ( elvec, __FILE__, __LINE__ )
 call fe_cell_function_previous%free() 
 call fe_cell_function_current%free() 
 deallocate ( H_current_curl_values, stat=istat); check(istat==0)
end subroutine hts_discrete_integration_integrate_residual 

! -----------------------------------------------------------------------------------------------
function compute_resistivity(this, H_value, curl_H, material) result(resistivity_tensor)
  implicit none 
  class(hts_discrete_integration_t)        , intent(in)    :: this
  type(vector_field_t)                     , intent(in)    :: H_value 
  type(vector_field_t)                     , intent(inout) :: curl_H 
  integer                                  , intent(in)    :: material
  type(tensor_field_t) :: resistivity_tensor

  ! Locals
  real(rp) :: resistivity, Ec, Jc, n 

  Ec = this%critical_electric_field 
  Jc = this%critical_current 
  n  = this%nonlinear_exponent 

  if ( this%is_analytical_solution ) then ! Analytical resistivity = rho * <H|H>**n 
     massert(this%hts_device_type==bulk, 'In analytical solutions device type must be an isotropic bulk') 
     resistivity = this%hts_resistivity * ( H_value*H_value )**n  
  else 
     if ( material == air ) then ! Air domain 
     resistivity = this%air_resistivity
     else if ( material > air ) then ! HTS DOMAIN: Nonlinear resistivity = Ec/Jc*|| curl(H) / Jc ||**n 
	 if ( n>0 ) then 
     resistivity = Ec/Jc*(curl_H%nrm2()/Jc)**n + 1e-14_rp
	 else 
	 resistivity = this%hts_resistivity
	 end if 
     else 
        massert(.false., 'Invalid material' ) 
     end if
  end if

  call resistivity_tensor%init(0.0_rp) 
  call resistivity_tensor%set(1,1,resistivity) 
  call resistivity_tensor%set(2,2,resistivity)
  select case ( this%hts_device_type ) 
  case ( bulk ) 
    call resistivity_tensor%set(3,3,resistivity) 
  case ( stack ) 
    call resistivity_tensor%set(3,3,this%air_resistivity )
  case DEFAULT 
    massert(.false., 'Invalid HTS device type' ) 
  end select 
          
end function compute_resistivity

! -----------------------------------------------------------------------------------------------
function compute_tangent_resistivity(this, H_value, curl_H, shape_value, curl_shape) result(tangent_resistivity_tensor) 
  implicit none 
  class(hts_discrete_integration_t), intent(in)        :: this
  type(vector_field_t)                , intent(in)     :: H_value
  type(vector_field_t)                , intent(inout)  :: curl_H 
  type(vector_field_t)                , intent(in)     :: shape_value 
  type(vector_field_t)                , intent(in)     :: curl_shape 
  type(tensor_field_t)   :: tangent_resistivity_tensor 
  ! Locals
  real(rp) :: tangent_resistivity, Ec, Jc, n 

  Ec = this%critical_electric_field 
  Jc = this%critical_current 
  n  = this%nonlinear_exponent 

  if ( this%is_analytical_solution ) then ! Analytical resistivity = rho * n * <H|H>**(n-1) * <Phi|H>
     massert(this%hts_device_type==bulk, 'In analytical solutions device type must be an isotropic bulk')
     if ( n > 0.0_rp) then 
        tangent_resistivity = this%hts_resistivity * 2.0_rp*n*( H_value*H_value )**(n-1.0_rp)*(shape_value * H_value)
     else 
        tangent_resistivity = 0.0_rp
     end if
  else 
     ! Tangent Nonlinear resistivity = Ec/Jc * n *|| curl(H) / Jc ||**(n-2) < curl(H)/Jc | curl(Phi)/Jc >
     if (this%nonlinear_exponent >= 1 ) then 
        tangent_resistivity = n*Ec/Jc*( curl_H%nrm2()/Jc )**(n-2.0_rp) * ((1.0_rp/Jc)*curl_H) * ((1.0_rp/Jc)*curl_shape) 
     else
        tangent_resistivity = 0.0_rp 
     end if
  end if
  
  call tangent_resistivity_tensor%init(0.0_rp) 
  call tangent_resistivity_tensor%set(1,1,tangent_resistivity) 
  call tangent_resistivity_tensor%set(2,2,tangent_resistivity)
  select case ( this%hts_device_type ) 
  case ( bulk ) 
    call tangent_resistivity_tensor%set(3,3,tangent_resistivity) 
  case ( stack ) 
    call tangent_resistivity_tensor%set(3,3,0.0_rp )
  case DEFAULT 
    massert(.false., 'Invalid HTS device type' ) 
  end select 
  
end function compute_tangent_resistivity

end module hts_discrete_integration_names
