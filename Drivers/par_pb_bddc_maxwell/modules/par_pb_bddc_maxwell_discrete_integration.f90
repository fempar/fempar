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
module par_pb_bddc_maxwell_discrete_integration_names
  use fempar_names
  use par_pb_bddc_maxwell_params_names
  use par_pb_bddc_maxwell_analytical_functions_names
  
  implicit none
# include "debug.i90"
  private
  
  integer(ip), parameter :: preconditioner = 0 
  integer(ip), parameter :: problem = 1  
  public :: preconditioner, problem 
  
  integer(ip), parameter :: BLACK =  0 
  integer(ip), parameter :: WHITE =  1
  public :: BLACK, WHITE 
  
  type, extends(discrete_integration_t) :: maxwell_cG_discrete_integration_t
   type(par_pb_bddc_maxwell_analytical_functions_t), pointer :: analytical_functions => NULL()
   type(fe_function_t)                 , pointer :: fe_function          => NULL()
   integer(ip)                                   :: integration_type
   character(len=:), allocatable                 :: materials_distribution_case 
   contains
     procedure :: set_analytical_functions
     procedure :: set_materials_distribution_case 
     procedure :: set_fe_function 
     procedure :: set_integration_type 
     procedure :: integrate_galerkin 
  end type maxwell_cG_discrete_integration_t
  
  public :: maxwell_cG_discrete_integration_t
  
contains
   
  subroutine set_analytical_functions ( this, analytical_functions )
     implicit none
     class(maxwell_cG_discrete_integration_t)    ,intent(inout)  :: this
     type(par_pb_bddc_maxwell_analytical_functions_t), target, intent(in)    :: analytical_functions
     this%analytical_functions => analytical_functions
  end subroutine set_analytical_functions
  
  subroutine set_fe_function (this, fe_function)
    implicit none
    class(maxwell_cG_discrete_integration_t), intent(inout)       :: this
    type(fe_function_t)                     , target, intent(in)  :: fe_function
    this%fe_function => fe_function
  end subroutine set_fe_function
      
  subroutine set_integration_type (this, integration_type )
    implicit none
    class(maxwell_cG_discrete_integration_t), intent(inout)  :: this
    integer(ip)                             , intent(in)     :: integration_type
    
    this%integration_type = integration_type
  end subroutine set_integration_type
  
  subroutine set_materials_distribution_case (this,  materials_distribution_case )
    implicit none
    class(maxwell_cG_discrete_integration_t), intent(inout)  :: this
    character(len=*)                        , intent(in)     :: materials_distribution_case
    
    this%materials_distribution_case = materials_distribution_case
  end subroutine set_materials_distribution_case

  subroutine integrate_galerkin ( this, fe_space, assembler )
    implicit none
    class(maxwell_cG_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)                , intent(inout) :: fe_space
    class(assembler_t)                      , intent(inout) :: assembler

    ! FE space traversal-related data types
    class(fe_cell_iterator_t) , allocatable :: fe

    ! FE integration-related data types
    type(quadrature_t)       , pointer :: quad
    type(point_t)            , pointer :: quad_coords(:)
    type(vector_field_t), allocatable  :: shape_curls(:,:)
    type(vector_field_t), allocatable  :: shape_values(:,:)

    ! FE matrix and vector i.e., A_K + f_K
    real(rp), allocatable              :: elmat(:,:), elvec(:)

    integer(ip)  :: istat
    integer(ip)  :: qpoint, num_quad_points
    integer(ip)  :: idof, jdof, num_dofs, max_num_dofs 
    real(rp)     :: factor
    type(vector_field_t), allocatable  :: source_term_values(:)

    integer(ip)                             :: number_fields
    class(vector_function_t)  , pointer     :: source_term
    type(curl_curl_coeff_holder_t), pointer :: curl_curl_coeff_holder(:)
    type(mass_coeff_holder_t), pointer      :: mass_coeff_holder(:)
    real(rp), allocatable                   :: mass_coeff(:), curl_curl_coeff(:) 
    integer(ip)                             :: material_id 
    
    assert (associated(this%analytical_functions)) 
    assert (associated(this%fe_function))

    ! Extract analytical functions 
    source_term         => this%analytical_functions%get_source_term()
    curl_curl_coeff_holder  => this%analytical_functions%get_curl_curl_coeff()
    mass_coeff_holder => this%analytical_functions%get_mass_coeff()
    
    call fe_space%set_up_cell_integration()
    call fe_space%create_fe_cell_iterator(fe)
    
    max_num_dofs = fe_space%get_max_num_dofs_on_a_cell()
    call memalloc ( max_num_dofs, max_num_dofs, elmat, __FILE__, __LINE__ )
    call memalloc ( max_num_dofs, elvec, __FILE__, __LINE__ )
 
    quad            => fe%get_quadrature()
    num_quad_points = quad%get_num_quadrature_points()
    allocate (source_term_values(num_quad_points), stat=istat); check(istat==0)
    call memalloc( num_quad_points, curl_curl_coeff, __FILE__, __LINE__ )
    call memalloc( num_quad_points, mass_coeff, __FILE__, __LINE__ )

    do while ( .not. fe%has_finished())
    
       if ( fe%is_local() ) then  
    
          ! Update FE-integration related data structures
          call fe%update_integration()

          ! Very important: this has to be inside the loop, as different FEs can be present! Study the multifield case, what happens with source term? 
          quad            => fe%get_quadrature()
          num_quad_points =  quad%get_num_quadrature_points()
          num_dofs        =  fe%get_num_dofs()
            
          ! Get quadrature coordinates to evaluate source_term
          quad_coords => fe%get_quadrature_points_coordinates()
    
          ! Evaluate analytical functions at quadrature points 
          select case ( this%materials_distribution_case )
          case ( homogeneous, heterogeneous )
          material_id = 1
          case ( checkerboard, channels, radial  ) 
          material_id = 1 + fe%get_set_id()
          end select 
          
          call source_term%get_values_set(quad_coords, source_term_values)
          call curl_curl_coeff_holder(material_id)%p%get_values_set(quad_coords, curl_curl_coeff)
          call mass_coeff_holder(material_id)%p%get_values_set(quad_coords, mass_coeff)  
                    
          ! Compute element matrix and vector
          elmat = 0.0_rp
          elvec = 0.0_rp
          call fe%get_curls(shape_curls)
          call fe%get_values(shape_values)
          do qpoint = 1, num_quad_points
             factor = fe%get_det_jacobian(qpoint) * quad%get_weight(qpoint)                        
             do idof = 1, num_dofs
                do jdof = 1, num_dofs
                   ! A_K(i,j) =  (curl(phi_i),curl(phi_j)) + (phi_i,phi_j)
                   elmat(idof,jdof) = elmat(idof,jdof) + factor * ( curl_curl_coeff(qpoint)* shape_curls(jdof,qpoint)* shape_curls(idof,qpoint) & 
                                                                  + mass_coeff(qpoint)*shape_values(jdof,qpoint)*shape_values(idof,qpoint) )
                end do
             end do
             
             do idof = 1, num_dofs
                ! F_K(i) = (f(i),phi_i)
                elvec(idof) = elvec(idof) + factor * source_term_values(qpoint)*shape_values(idof,qpoint)
             end do
          end do
           
          ! Apply boundary conditions
          call fe%assembly( this%fe_function, elmat, elvec, assembler )
       elseif ( this%integration_type == preconditioner ) then 
       
          ! Update FE-integration related data structures
          call fe%update_integration()
       
          ! Get parameter values 
          select case ( this%materials_distribution_case )
          case ( homogeneous, heterogeneous )
          material_id = 1
          case ( checkerboard, channels, radial ) 
          material_id = 1 + fe%get_set_id()
          end select 
                  
          ! Very important: this has to be inside the loop, as different FEs can be present! Study the multifield case, what happens with source term? 
          quad            => fe%get_quadrature()
          num_quad_points =  quad%get_num_quadrature_points()
          num_dofs        =  fe%get_num_dofs()
            
          ! Get quadrature coordinates to evaluate source_term
          quad_coords => fe%get_quadrature_points_coordinates()
          call mass_coeff_holder(material_id)%p%get_values_set(quad_coords, mass_coeff)   
                            
          ! Compute element matrix and vector
          elmat = 0.0_rp
          call fe%get_curls(shape_curls)
          call fe%get_values(shape_values)
          do qpoint = 1, num_quad_points
             factor = fe%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
          
             do idof = 1, num_dofs
                do jdof = 1, num_dofs
                   ! A_M(i,j) =  (phi_i,phi_j)
                   elmat(idof,jdof) = elmat(idof,jdof) + factor * mass_coeff(qpoint)*shape_values(jdof,qpoint)*shape_values(idof,qpoint)
                end do
             end do
          end do
          
          ! Assemble elmat to active DoFs on ghost cells 
            call fe%assembly( elmat, assembler )
       end if 
       call fe%next()
    end do
    
    call fe_space%free_fe_cell_iterator(fe)

    deallocate (shape_values, stat=istat);       check(istat==0)
    deallocate (shape_curls, stat=istat);        check(istat==0)
    deallocate (source_term_values, stat=istat); check(istat==0)
    call memfree ( curl_curl_coeff, __FILE__, __LINE__ ) 
    call memfree ( mass_coeff, __FILE__, __LINE__ ) 
    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )
  end subroutine integrate_galerkin
  
end module par_pb_bddc_maxwell_discrete_integration_names
