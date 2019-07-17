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

!****************************************************************************************************

module tutorial_03_discrete_integration_names
  use fempar_names
  
  implicit none
# include "debug.i90"
  private
  type, extends(discrete_integration_t) :: poisson_discrete_integration_t
     private
     class(scalar_function_t), pointer :: source_term 
     type(fe_function_t),      pointer :: fe_function => NULL()
   contains
     procedure :: set_source_term
     procedure :: set_fe_function
     procedure :: integrate_galerkin
  end type poisson_discrete_integration_t
  
  public :: poisson_discrete_integration_t
  
contains
  subroutine set_source_term ( this, source_term )
     implicit none
     class(poisson_discrete_integration_t)        , intent(inout) :: this
     class(scalar_function_t),              target, intent(in)    :: source_term
     this%source_term => source_term
  end subroutine set_source_term

  subroutine set_fe_function (this, fe_function)
     implicit none
     class(poisson_discrete_integration_t), intent(inout) :: this
     type(fe_function_t)             , target, intent(in)    :: fe_function
     this%fe_function => fe_function
  end subroutine set_fe_function

  subroutine integrate_galerkin ( this, fe_space, assembler )
    implicit none
    class(poisson_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)         , intent(inout) :: fe_space
    class(assembler_t)      , intent(inout) :: assembler

    ! FE space traversal-related data types
    class(fe_cell_iterator_t), allocatable :: fe

    ! FE integration-related data types
    type(quadrature_t)       , pointer :: quad
    type(point_t)            , pointer :: quad_coords(:)
    type(vector_field_t), allocatable  :: shape_gradients(:,:)
    real(rp)            , allocatable  :: shape_values(:,:)

    ! FE matrix and vector i.e., A_K + f_K
    real(rp), allocatable              :: elmat(:,:), elvec(:)

    integer(ip)  :: istat
    integer(ip)  :: qpoint, num_quad_points
    integer(ip)  :: idof, jdof, num_dofs, max_num_dofs
    real(rp)     :: factor
    real(rp)     :: source_term_value
    
    call fe_space%create_fe_cell_iterator(fe)
    max_num_dofs = fe_space%get_max_num_dofs_on_a_cell()
    call memalloc ( max_num_dofs, max_num_dofs, elmat, __FILE__, __LINE__ )
    call memalloc ( max_num_dofs, elvec, __FILE__, __LINE__ )

    do while ( .not. fe%has_finished() )

       if ( fe%is_local() ) then
          ! Update FE-integration related data structures
          call fe%update_integration()

          ! Very important: this has to be inside the loop, as different FEs can be present!
          quad            => fe%get_quadrature()
          num_quad_points = quad%get_num_quadrature_points()
          num_dofs = fe%get_num_dofs()
          
          ! Get quadrature coordinates to evaluate source_term
          quad_coords => fe%get_quadrature_points_coordinates()
                    
          ! Compute element matrix and vector
          elmat = 0.0_rp
          elvec = 0.0_rp
          call fe%get_gradients(shape_gradients)
          call fe%get_values(shape_values)
          do qpoint = 1, num_quad_points
             factor = fe%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
             do idof = 1, num_dofs
                do jdof = 1, num_dofs
                   ! A_K(i,j) = (grad(phi_i),grad(phi_j))
                   elmat(idof,jdof) = elmat(idof,jdof) + factor * shape_gradients(jdof,qpoint) * shape_gradients(idof,qpoint)
                end do
             end do
             
             ! Source term
             call this%source_term%get_value_space(quad_coords(qpoint),source_term_value)
             do idof = 1, num_dofs
                elvec(idof) = elvec(idof) + factor * source_term_value * shape_values(idof,qpoint) 
             end do
          end do
          
          call fe%assembly( this%fe_function, elmat, elvec, assembler )
       end if
       call fe%next()
    end do
    call fe_space%free_fe_cell_iterator(fe)
    call memfree(shape_values, __FILE__, __LINE__)
    deallocate (shape_gradients, stat=istat); check(istat==0);
    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )
  end subroutine integrate_galerkin
  
end module tutorial_03_discrete_integration_names

module tutorial_03_functions_names
  use fempar_names
  implicit none
# include "debug.i90"
  private
  
  type, extends(scalar_function_t) :: base_scalar_function_t
    private
    real(rp)              :: alpha
    real(rp)              :: circle_radius
    real(rp), allocatable :: circle_center(:)
  contains
    procedure, non_overridable :: create  => base_scalar_function_create
    procedure, non_overridable :: free    => base_scalar_function_free
  end type base_scalar_function_t
  
  type, extends(base_scalar_function_t) :: sharp_circular_wave_source_term_t
    private
  contains
    procedure, non_overridable :: get_value_space => sharp_circular_wave_source_term_get_value_space
  end type sharp_circular_wave_source_term_t
  
  type, extends(base_scalar_function_t) :: sharp_circular_wave_solution_t
    private
  contains
    procedure, non_overridable :: get_value_space    => sharp_circular_wave_solution_get_value_space
    procedure, non_overridable :: get_gradient_space => sharp_circular_wave_solution_get_gradient_space
  end type sharp_circular_wave_solution_t
  
  public :: sharp_circular_wave_source_term_t
  public :: sharp_circular_wave_solution_t
  
contains
  !===============================================================================================
  subroutine base_scalar_function_create ( this, num_dims, alpha, circle_radius, circle_center )
    implicit none
    class(base_scalar_function_t), intent(inout) :: this
    integer(ip)                  , intent(in)    :: num_dims
    real(rp)                     , intent(in)    :: alpha
    real(rp)                     , intent(in)    :: circle_radius
    real(rp)                     , intent(in)    :: circle_center(:)
    call this%free()
    massert ( num_dims == size(circle_center), "base_scalar_function_create_t :: num_dims MUST match size(circle_center)" )
    call this%set_num_dims(num_dims)
    this%alpha = alpha
    this%circle_radius = circle_radius
    call memalloc ( this%get_num_dims(), this%circle_center, __FILE__, __LINE__ ) 
    this%circle_center(:) = circle_center(:)
  end subroutine  base_scalar_function_create 
 
  !===============================================================================================
  subroutine base_scalar_function_free ( this)
    implicit none
    class(base_scalar_function_t), intent(inout) :: this
    call this%set_num_dims(-1_ip)
    this%alpha         = 0.0_rp
    this%circle_radius = 0.0_rp
    if (allocated(this%circle_center) ) call memfree(this%circle_center, __FILE__, __LINE__ )
  end subroutine  base_scalar_function_free 
  
  !===============================================================================================
  subroutine sharp_circular_wave_source_term_get_value_space ( this, point, result )
    implicit none
    class(sharp_circular_wave_source_term_t), intent(in)    :: this
    type(point_t)                           , intent(in)    :: point
    real(rp)                                , intent(inout) :: result
    real(rp)    :: term1, sqrt_term1
    real(rp)    :: term2, term2_2
    real(rp)    :: alpha, alpha_2, alpha_3
    real(rp)    :: x_minus_xc_2, y_minus_yc_2, z_minus_zc_2
    real(rp)    :: xc, yc, zc
    real(rp)    :: x, y, z, r
    real(rp)    :: dx2, dy2, dz2
    integer(ip) :: num_dims
    num_dims = this%get_num_dims()
    alpha = this%alpha
    alpha_2 = alpha*alpha
    alpha_3 = alpha_2*alpha
    xc = this%circle_center(1)
    yc = this%circle_center(2)
    r  = this%circle_radius
    x  = point%get(1); call assert_if_not_within_range(x)
    y  = point%get(2); call assert_if_not_within_range(y)
    if ( num_dims == 2 ) then
      x_minus_xc_2 = (x-xc)*(x-xc)
      y_minus_yc_2 = (y-yc)*(y-yc)
      term1 = x_minus_xc_2+y_minus_yc_2
      sqrt_term1 = sqrt(term1)
      term2   = sqrt_term1-r
      term2_2 = term2*term2
      dx2 = alpha                              /(sqrt_term1*(alpha_2*term2_2+1.0_rp)) - & 
            (alpha*x_minus_xc_2)               /((term1**1.5_rp)*(alpha_2*term2_2+1.0_rp)) - &
            (2.0_rp*alpha_3*x_minus_xc_2*term2)/(term1*(alpha_2*term2_2+1.0_rp)*(alpha_2*term2_2+1.0_rp))
      dy2 = alpha                              /(sqrt_term1*(alpha_2*term2_2+1.0_rp)) - & 
            (alpha*y_minus_yc_2)               /((term1**1.5_rp)*(alpha_2*term2_2+1.0_rp)) - &
            (2.0_rp*alpha_3*y_minus_yc_2*term2)/(term1*(alpha_2*term2_2+1.0_rp)*(alpha_2*term2_2+1.0_rp))
      result = -dx2-dy2
    else
      z  = point%get(3); call assert_if_not_within_range(y)  
      zc = this%circle_center(3)
      x_minus_xc_2 = (x-xc)*(x-xc)
      y_minus_yc_2 = (y-yc)*(y-yc)
      z_minus_zc_2 = (z-zc)*(z-zc)
      term1 = x_minus_xc_2+y_minus_yc_2+z_minus_zc_2
      sqrt_term1 = sqrt(term1)
      term2   = sqrt_term1-r
      term2_2 = term2*term2
      dx2 = alpha                              /(sqrt_term1*(alpha_2*term2_2+1.0_rp)) - & 
            (alpha*x_minus_xc_2)               /((term1**1.5_rp)*(alpha_2*term2_2+1.0_rp)) - &
            (2.0_rp*alpha_3*x_minus_xc_2*term2)/(term1*(alpha_2*term2_2+1.0_rp)*(alpha_2*term2_2+1.0_rp))
      dy2 = alpha                              /(sqrt_term1*(alpha_2*term2_2+1.0_rp)) - & 
            (alpha*y_minus_yc_2)               /((term1**1.5_rp)*(alpha_2*term2_2+1.0_rp)) - &
            (2.0_rp*alpha_3*y_minus_yc_2*term2)/(term1*(alpha_2*term2_2+1.0_rp)*(alpha_2*term2_2+1.0_rp))
      dz2 = alpha                              /(sqrt_term1*(alpha_2*term2_2+1.0_rp)) - & 
            (alpha*z_minus_zc_2)               /((term1**1.5_rp)*(alpha_2*term2_2+1.0_rp)) - &
            (2.0_rp*alpha_3*z_minus_zc_2*term2)/(term1*(alpha_2*term2_2+1.0_rp)*(alpha_2*term2_2+1.0_rp))      
      result = -dx2-dy2-dz2
    end if
  end subroutine sharp_circular_wave_source_term_get_value_space
  
  !===============================================================================================
  subroutine sharp_circular_wave_solution_get_value_space ( this, point, result )
    implicit none
    class(sharp_circular_wave_solution_t), intent(in)    :: this
    type(point_t)                        , intent(in)    :: point
    real(rp)                             , intent(inout) :: result
    real(rp)    :: alpha
    real(rp)    :: xc, yc, zc
    real(rp)    :: x, y, z, r
    integer(ip) :: num_dims
    num_dims = this%get_num_dims()
    alpha = this%alpha
    xc = this%circle_center(1)
    yc = this%circle_center(2)
    r  = this%circle_radius
    x  = point%get(1); call assert_if_not_within_range(x)
    y  = point%get(2); call assert_if_not_within_range(y)
    if ( num_dims == 2 ) then
      result = atan(alpha*(sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc))-r))
    else if (num_dims == 3) then
      zc = this%circle_center(3)
      z  = point%get(3); ; call assert_if_not_within_range(z)
      result = atan(alpha*(sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc))-r))
    end if
  end subroutine sharp_circular_wave_solution_get_value_space
  
  subroutine sharp_circular_wave_solution_get_gradient_space ( this, point, result )
    implicit none
    class(sharp_circular_wave_solution_t), intent(in)    :: this
    type(point_t)                        , intent(in)    :: point
    type(vector_field_t)                 , intent(inout) :: result
    real(rp)    :: term1, sqrt_term1
    real(rp)    :: term2, term2_2
    real(rp)    :: alpha, alpha_2
    real(rp)    :: x_minus_xc_2, y_minus_yc_2, z_minus_zc_2
    real(rp)    :: xc, yc, zc
    real(rp)    :: x, y, z, r
    real(rp)    :: dx, dy, dz
    integer(ip) :: num_dims
    num_dims = this%get_num_dims()
    alpha = this%alpha
    alpha_2 = alpha*alpha
    xc = this%circle_center(1)
    yc = this%circle_center(2)
    r  = this%circle_radius
    x  = point%get(1); call assert_if_not_within_range(x)
    y  = point%get(2); call assert_if_not_within_range(y)
    x_minus_xc_2 = (x-xc)*(x-xc)
    y_minus_yc_2 = (y-yc)*(y-yc)
    if ( num_dims == 2 ) then 
      term1 = x_minus_xc_2+y_minus_yc_2
      sqrt_term1 = sqrt(term1)
      term2 = sqrt_term1-r
      term2_2 = term2*term2
      dx = (alpha*(x-xc))/(sqrt_term1*(alpha_2*term2_2+1.0_rp))
      dy = (alpha*(y-yc))/(sqrt_term1*(alpha_2*term2_2+1.0_rp))
      call result%set(1, dx)
      call result%set(2, dy)
    else
      zc = this%circle_center(3)
      z  = point%get(3); call assert_if_not_within_range(z)
      z_minus_zc_2 = (z-zc)*(z-zc)
      term1 = x_minus_xc_2+y_minus_yc_2+z_minus_zc_2
      sqrt_term1 = sqrt(term1)
      term2 = sqrt_term1-r
      term2_2 = term2*term2
      dx = (alpha*(x-xc))/(sqrt_term1*(alpha_2*term2_2+1.0_rp))
      dy = (alpha*(y-yc))/(sqrt_term1*(alpha_2*term2_2+1.0_rp))
      dz = (alpha*(z-zc))/(sqrt_term1*(alpha_2*term2_2+1.0_rp))
      call result%set(1, dx)
      call result%set(2, dy)
      call result%set(3, dz)
    end if
  end subroutine sharp_circular_wave_solution_get_gradient_space 
  
  subroutine assert_if_not_within_range(x)
    real(rp), intent(in) :: x
    massert (x >= 0.0_rp .and. x <= 1.0_rp, "sharp_circular_wave_solution_t :: this function only works with unit cube/square")
  end subroutine assert_if_not_within_range
  
end module tutorial_03_functions_names

module tutorial_03_error_estimator_names
  use fempar_names
  
  implicit none
# include "debug.i90"
  private
  
  
  type, extends(error_estimator_t) :: poisson_error_estimator_t
     private
     class(scalar_function_t)            , pointer     :: exact_solution    => NULL()
     type(fe_function_t)                 , pointer     :: discrete_solution => NULL()
     type(fe_cell_function_scalar_t)                   :: fe_cell_function
     type(vector_field_t)                , allocatable :: exact_solution_gradients(:) 
   contains
     procedure :: create                    => poisson_error_estimator_create
     procedure :: free                      => poisson_error_estimator_free
     procedure :: set_exact_solution        => poisson_error_estimator_set_exact_solution
     procedure :: set_discrete_solution     => poisson_error_estimator_set_discrete_solution
     procedure :: compute_local_estimates   => poisson_error_estimator_compute_local_estimates
     procedure :: compute_local_true_errors => poisson_error_estimator_compute_local_true_errors
     procedure :: get_error_norm_exponent   => poisson_error_estimator_get_error_norm_exponent
  end type poisson_error_estimator_t
  
  public :: poisson_error_estimator_t

contains

  subroutine poisson_error_estimator_create ( this, fe_space, parameter_list )
    implicit none
    class(poisson_error_estimator_t), intent(inout) :: this
    class(serial_fe_space_t)   , target, intent(in) :: fe_space
    type(parameterlist_t)              , intent(in) :: parameter_list
    class(environment_t), pointer :: environment
    integer(ip) :: istat
    call this%free()
    call ee_create(this,fe_space,parameter_list)
    assert(fe_space%get_num_fields()==1)
    call this%fe_cell_function%create(fe_space,1)
    allocate(this%exact_solution_gradients(fe_space%get_max_num_quadrature_points()),stat=istat)
    check(istat==0)
  end subroutine poisson_error_estimator_create

  subroutine poisson_error_estimator_free ( this )
    implicit none
    class(poisson_error_estimator_t), intent(inout) :: this
    integer(ip)                   :: istat
    call this%fe_cell_function%free()
    if ( allocated(this%exact_solution_gradients) ) then
       deallocate(this%exact_solution_gradients, stat=istat)
       check(istat==0)
    end if
    call ee_free(this)
  end subroutine poisson_error_estimator_free

  subroutine poisson_error_estimator_set_exact_solution ( this, exact_solution )
    implicit none
    class(poisson_error_estimator_t)        , intent(inout) :: this
    class(scalar_function_t)       , target , intent(in)    :: exact_solution
    this%exact_solution => exact_solution
  end subroutine poisson_error_estimator_set_exact_solution

  subroutine poisson_error_estimator_set_discrete_solution ( this, discrete_solution )
    implicit none
    class(poisson_error_estimator_t)         , intent(inout) :: this
    type(fe_function_t)             , target , intent(in)    :: discrete_solution
    this%discrete_solution => discrete_solution
  end subroutine poisson_error_estimator_set_discrete_solution

  subroutine poisson_error_estimator_compute_local_true_errors(this)
    implicit none
    class(poisson_error_estimator_t), intent(inout) :: this
    class(serial_fe_space_t)  , pointer :: fe_space
    class(triangulation_t)    , pointer :: triangulation
    type(vector_field_t)      , pointer :: discrete_solution_gradients(:)
    type(std_vector_real_rp_t), pointer :: sq_local_true_errors
    real(rp)                  , pointer :: sq_local_true_errors_entries(:)
    class(fe_cell_iterator_t) , allocatable :: fe
    type(quadrature_t)        , pointer :: quad
    type(point_t)             , pointer :: quad_coords(:)
    real(rp) :: factor
    integer(ip) :: qpoint, num_quad_points
    real(rp) :: sq_local_estimate
    integer(ip) :: istat 

    assert (associated(this%exact_solution))
    assert (associated(this%discrete_solution))
    fe_space => this%get_fe_space()
    assert ( associated(fe_space) )
    
    triangulation => fe_space%get_triangulation()

    sq_local_true_errors => this%get_sq_local_true_errors()
    call sq_local_true_errors%resize(0)
    call sq_local_true_errors%resize(triangulation%get_num_cells(), 0.0_rp)
    sq_local_true_errors_entries => sq_local_true_errors%get_pointer()

    call fe_space%create_fe_cell_iterator(fe)
    do while(.not. fe%has_finished())
       call fe%update_integration()
       call this%fe_cell_function%update(fe,this%discrete_solution)
       quad => fe%get_quadrature()
       num_quad_points = quad%get_num_quadrature_points()
       quad_coords => fe%get_quadrature_points_coordinates()
       call this%exact_solution%get_gradients_set( quad_coords, &
                                                   this%exact_solution_gradients )
       discrete_solution_gradients => this%fe_cell_function%get_quadrature_points_gradients()
       sq_local_estimate = 0.0_rp
       do qpoint = 1, num_quad_points
          factor = fe%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
          this%exact_solution_gradients(qpoint) = this%exact_solution_gradients(qpoint) - &
                                                  discrete_solution_gradients(qpoint)
          sq_local_estimate = sq_local_estimate + factor * this%exact_solution_gradients(qpoint) * &
                                     this%exact_solution_gradients(qpoint)
       end do
       sq_local_true_errors_entries(fe%get_gid()) = sq_local_estimate
       call fe%next()
    end do
    call fe_space%free_fe_cell_iterator(fe)
  end subroutine poisson_error_estimator_compute_local_true_errors

  subroutine poisson_error_estimator_compute_local_estimates(this)
    class(poisson_error_estimator_t), intent(inout) :: this
  end subroutine poisson_error_estimator_compute_local_estimates

  function poisson_error_estimator_get_error_norm_exponent(this)
    class(poisson_error_estimator_t), intent(in) :: this
    real(rp) :: poisson_error_estimator_get_error_norm_exponent
    poisson_error_estimator_get_error_norm_exponent = 0.5_rp
  end function poisson_error_estimator_get_error_norm_exponent

end module tutorial_03_error_estimator_names

