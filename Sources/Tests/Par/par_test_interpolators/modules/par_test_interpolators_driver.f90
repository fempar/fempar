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
module test_interpolators_driver_names
  use fempar_names
  use par_test_interpolators_params_names
  use interpolators_analytical_functions_names
  use interpolators_conditions_names
# include "debug.i90"

  implicit none
  private

  integer(ip), parameter :: MAGNETIC_FIELD_ID = 1
  integer(ip), parameter :: PRESSURE_FIELD_ID = 2

  type par_test_interpolators_driver_t 
     private 

     ! Place-holder for parameter-value set provided through command-line interface
     type(par_test_interpolators_params_t)     :: test_params
     type(ParameterList_t) , pointer           :: parameter_list

     ! Cells and lower dimension objects container
     type(par_triangulation_t)                 :: triangulation
     integer(ip), allocatable                  :: cell_set_ids(:)

     ! Analytical functions of the problem
     type(interpolators_analytical_functions_t) :: problem_functions

     ! Discrete weak problem integration-related data type instances 
     type(par_fe_space_t)                        :: fe_space 
     type(p_reference_fe_t) , allocatable        :: reference_fes(:) 
     type(standard_l1_coarse_fe_handler_t)       :: coarse_fe_handler 
     type(p_l1_coarse_fe_handler_t), allocatable :: coarse_fe_handlers(:)
     integer(ip)            , allocatable        :: set_ids_to_reference_fes(:,:)
     type(interpolators_conditions_t)            :: interpolators_conditions

     ! Operators related data 
     type(block_layout_t)                        :: block_layout 

     ! Problem solution FE function
     type(fe_function_t)                         :: solution
     type(fe_function_t)                         :: time_solution 
     real(rp)                                    :: time 

     type(environment_t)                    :: par_environment

   contains
     procedure                  :: run_simulation
     procedure                  :: setup_environment
     procedure                  :: parse_command_line_parameters
     procedure        , private :: setup_triangulation
     procedure        , private :: setup_reference_fes
     procedure        , private :: setup_coarse_fe_handlers 
     procedure        , private :: setup_fe_space
     procedure        , private :: interpolate_analytical_functions 
     procedure        , private :: check_solution
     procedure        , private :: check_time_solution
     procedure                  :: free_command_line_parameters
     procedure                  :: free_environment 
     procedure        , private :: free
  end type par_test_interpolators_driver_t

  ! Types
  public :: par_test_interpolators_driver_t

contains

  subroutine parse_command_line_parameters(this)
    implicit none
    class(par_test_interpolators_driver_t), intent(inout) :: this 
    call this%test_params%create()
    this%parameter_list => this%test_params%get_values()
  end subroutine parse_command_line_parameters

  !========================================================================================
  subroutine setup_environment(this)
    implicit none
    class(par_test_interpolators_driver_t), intent(inout) :: this
    integer(ip) :: istat
    if ( this%test_params%get_triangulation_type() == triangulation_generate_structured ) then
       istat = this%parameter_list%set(key = environment_type_key, value = structured) ; check(istat==0)
    else
       istat = this%parameter_list%set(key = environment_type_key, value = unstructured) ; check(istat==0)
    end if
    istat = this%parameter_list%set(key = execution_context_key, value = mpi_context) ; check(istat==0)
    call this%par_environment%create(this%parameter_list)
  end subroutine setup_environment

  subroutine setup_triangulation(this)
    implicit none
    class(par_test_interpolators_driver_t), intent(inout) :: this
    class(cell_iterator_t)            , allocatable :: cell
    integer(ip)                                     :: set_id 

    call this%triangulation%create(this%parameter_list, this%par_environment)

  end subroutine setup_triangulation

  subroutine setup_reference_fes(this)
    implicit none
    class(par_test_interpolators_driver_t), intent(inout) :: this
    integer(ip) :: istat
    class(vef_iterator_t), allocatable  :: vef
    class(cell_iterator_t), allocatable       :: cell
    class(reference_fe_t), pointer :: reference_fe_geo

       allocate(this%reference_fes(2), stat=istat); check(istat==0)

       if ( this%par_environment%am_i_l1_task() ) then
          call this%triangulation%create_cell_iterator(cell)
          reference_fe_geo => cell%get_reference_fe()

          this%reference_fes(MAGNETIC_FIELD_ID) =  make_reference_fe ( topology   = reference_fe_geo%get_topology(),           &
                                                                       fe_type    = fe_type_nedelec,                           &
                                                                       num_dims   = this%triangulation%get_num_dims(),         &
                                                                       order      = this%test_params%get_reference_fe_order(), &
                                                                       field_type = field_type_vector,                         &
                                                                       conformity = .true. ) 

          this%reference_fes(PRESSURE_FIELD_ID) =  make_reference_fe ( topology   = reference_fe_geo%get_topology(),           &
                                                                       fe_type    = fe_type_lagrangian,                        &
                                                                       num_dims   = this%triangulation%get_num_dims(),         &
                                                                       order      = this%test_params%get_reference_fe_order(), &
                                                                       field_type = field_type_scalar,                         &
                                                                       conformity = .true. ) 
          call this%triangulation%free_cell_iterator(cell)
       end if

    ! if ( trim(this%test_params%get_triangulation_type()) == 'structured' ) then
    call this%triangulation%create_vef_iterator(vef)
    do while ( .not. vef%has_finished() )
       if(vef%is_at_boundary()) then
          call vef%set_set_id(1)
       else
          call vef%set_set_id(0)
       end if
       call vef%next()
    end do
    call this%triangulation%free_vef_iterator(vef)
    ! end if    

  end subroutine setup_reference_fes

  subroutine setup_coarse_fe_handlers(this)
    implicit none
    class(par_test_interpolators_driver_t), target, intent(inout) :: this
    integer(ip) :: istat, i 

    allocate(this%coarse_fe_handlers(2), stat=istat); check(istat==0)
    do i=1,2
       this%coarse_fe_handlers(i)%p => this%coarse_fe_handler
    end do

  end subroutine setup_coarse_fe_handlers

  subroutine setup_fe_space(this)
    implicit none
    class(par_test_interpolators_driver_t), intent(inout) :: this 

    call this%interpolators_conditions%set_num_dims(this%triangulation%get_num_dims() + 1)

       call this%fe_space%create( triangulation       = this%triangulation,      &
                                  reference_fes       = this%reference_fes,      &
                                  coarse_fe_handlers  = this%coarse_fe_handlers, &
                                  conditions          = this%interpolators_conditions )
    call this%fe_space%set_up_cell_integration()
    call this%fe_space%set_up_facet_integration()

    ! Operators related data needed for fe_function 
    !call this%block_layout%create( this%fe_space%get_num_fields() )
    !call this%fe_space%generate_global_dof_numbering()        
  end subroutine setup_fe_space

  subroutine interpolate_analytical_functions (this)
    implicit none
    class(par_test_interpolators_driver_t), intent(inout) :: this 
    class(vector_t) , pointer :: dof_values
    real(rp) :: time_ 

    call this%problem_functions%set_num_dims(this%triangulation%get_num_dims())

    ! Set boundary conditions 
    call this%interpolators_conditions%set_boundary_function_Hx(this%problem_functions%get_boundary_function_Hx())
    call this%interpolators_conditions%set_boundary_function_Hy(this%problem_functions%get_boundary_function_Hy())
    if ( this%triangulation%get_num_dims() == 3) then 
       call this%interpolators_conditions%set_boundary_function_Hz(this%problem_functions%get_boundary_function_Hz())
    end if
    call this%interpolators_conditions%set_boundary_function_pressure(this%problem_functions%get_boundary_function_pressure()) 

    ! Project functions 
    call this%solution%create(this%fe_space)
    call this%fe_space%interpolate( MAGNETIC_FIELD_ID, this%problem_functions%get_magnetic_field_solution(), this%solution ) 
    call this%fe_space%interpolate( PRESSURE_FIELD_ID, this%problem_functions%get_pressure_solution(), this%solution )
    call this%fe_space%interpolate_dirichlet_values( this%solution )

    ! Project transient functions 
    call random_number( this%time ) 
    call this%time_solution%create(this%fe_space)
    call this%fe_space%interpolate(  field_id    = MAGNETIC_FIELD_ID,                                    & 
                                     function     = this%problem_functions%get_magnetic_field_solution(), &
                                     fe_function  = this%time_solution ,                                  & 
                                     time         = this%time ) 
    call this%fe_space%interpolate( field_id     = PRESSURE_FIELD_ID,                                    &
                                    function     = this%problem_functions%get_pressure_solution(),       & 
                                    fe_function  = this%time_solution ,                                  &  
                                    time         = this%time)
    call this%fe_space%interpolate_dirichlet_values( fe_function = this%time_solution, & 
                                                     time        = this%time )
  end subroutine interpolate_analytical_functions

  subroutine check_solution(this)
    implicit none
    class(par_test_interpolators_driver_t), intent(inout) :: this
    class(vector_function_t), pointer :: H_exact_function
    class(scalar_function_t), pointer :: p_exact_function
    type(error_norms_vector_t) :: H_error_norm
    type(error_norms_scalar_t) :: p_error_norm 
    real(rp) :: mean, l1, l2, lp, linfty, h1, hcurl, h1_s, w1p_s, w1p, w1infty_s, w1infty
    real(rp) :: error_tolerance

    H_exact_function => this%problem_functions%get_magnetic_field_solution()
    p_exact_function => this%problem_functions%get_pressure_solution() 

    error_tolerance = 1.0e-06 

    call H_error_norm%create(this%fe_space, MAGNETIC_FIELD_ID )
    mean = H_error_norm%compute(H_exact_function, this%solution, mean_norm)   
    l1 = H_error_norm%compute(H_exact_function, this%solution, l1_norm)   
    l2 = H_error_norm%compute(H_exact_function, this%solution, l2_norm)   
    lp = H_error_norm%compute(H_exact_function, this%solution, lp_norm)   
    linfty = H_error_norm%compute(H_exact_function, this%solution, linfty_norm)   
    h1_s = H_error_norm%compute(H_exact_function, this%solution, h1_seminorm) 
    h1 = H_error_norm%compute(H_exact_function, this%solution, h1_norm) 
    hcurl = H_error_norm%compute(H_exact_function, this%solution, hcurl_seminorm) 
    w1p_s = H_error_norm%compute(H_exact_function, this%solution, w1p_seminorm)   
    w1p = H_error_norm%compute(H_exact_function, this%solution, w1p_norm)   
    w1infty_s = H_error_norm%compute(H_exact_function, this%solution, w1infty_seminorm) 
    w1infty = H_error_norm%compute(H_exact_function, this%solution, w1infty_norm)

    if ( this%par_environment%am_i_l1_root() ) then 
    write(*,*) 'PROJECTED MAGNETIC FIELD FUNCTION ERROR NORMS *************'
    write(*,'(a20,e32.25)') 'mean_norm:', mean; check ( abs(mean) < error_tolerance )
    write(*,'(a20,e32.25)') 'l1_norm:', l1; check ( l1 < error_tolerance )
    write(*,'(a20,e32.25)') 'l2_norm:', l2; check ( l2 < error_tolerance )
    write(*,'(a20,e32.25)') 'lp_norm:', lp; check ( lp < error_tolerance )
    write(*,'(a20,e32.25)') 'linfnty_norm:', linfty; check ( linfty < error_tolerance )
    write(*,'(a20,e32.25)') 'h1_seminorm:', h1_s; check ( h1_s < error_tolerance )
    write(*,'(a20,e32.25)') 'h1_norm:', h1; check ( h1 < error_tolerance )
    write(*,'(a20,e32.25)') 'hcurl_norm:', hcurl; check ( hcurl < error_tolerance )
    write(*,'(a20,e32.25)') 'w1p_seminorm:', w1p_s; check ( w1p_s < error_tolerance )
    write(*,'(a20,e32.25)') 'w1p_norm:', w1p; check ( w1p < error_tolerance )
    write(*,'(a20,e32.25)') 'w1infty_seminorm:', w1infty_s; check ( w1infty_s < error_tolerance )
    write(*,'(a20,e32.25)') 'w1infty_norm:', w1infty; check ( w1infty < error_tolerance )
    end if 
    call H_error_norm%free()


    call p_error_norm%create(this%fe_space, PRESSURE_FIELD_ID )
    mean = p_error_norm%compute(p_exact_function, this%solution, mean_norm)   
    l1 = p_error_norm%compute(p_exact_function, this%solution, l1_norm)   
    l2 = p_error_norm%compute(p_exact_function, this%solution, l2_norm)   
    lp = p_error_norm%compute(p_exact_function, this%solution, lp_norm)   
    linfty = p_error_norm%compute(p_exact_function, this%solution, linfty_norm)   
    h1_s = p_error_norm%compute(p_exact_function, this%solution, h1_seminorm) 
    h1 = p_error_norm%compute(p_exact_function, this%solution, h1_norm) 
    w1p_s = p_error_norm%compute(p_exact_function, this%solution, w1p_seminorm)   
    w1p = p_error_norm%compute(p_exact_function, this%solution, w1p_norm)   
    w1infty_s = p_error_norm%compute(p_exact_function, this%solution, w1infty_seminorm) 
    w1infty = p_error_norm%compute(p_exact_function, this%solution, w1infty_norm)

    if (this%par_environment%am_i_l1_root()) then 
    write(*,*) 'PROJECTED SCALAR PRESSURE FUNCTION ERROR NORMS ************************'
    write(*,'(a20,e32.25)') 'mean_norm:', mean; check ( abs(mean) < error_tolerance )
    write(*,'(a20,e32.25)') 'l1_norm:', l1; check ( l1 < error_tolerance )
    write(*,'(a20,e32.25)') 'l2_norm:', l2; check ( l2 < error_tolerance )
    write(*,'(a20,e32.25)') 'lp_norm:', lp; check ( lp < error_tolerance )
    write(*,'(a20,e32.25)') 'linfnty_norm:', linfty; check ( linfty < error_tolerance )
    write(*,'(a20,e32.25)') 'h1_seminorm:', h1_s; check ( h1_s < error_tolerance )
    write(*,'(a20,e32.25)') 'h1_norm:', h1; check ( h1 < error_tolerance )
    write(*,'(a20,e32.25)') 'w1p_seminorm:', w1p_s; check ( w1p_s < error_tolerance )
    write(*,'(a20,e32.25)') 'w1p_norm:', w1p; check ( w1p < error_tolerance )
    write(*,'(a20,e32.25)') 'w1infty_seminorm:', w1infty_s; check ( w1infty_s < error_tolerance )
    write(*,'(a20,e32.25)') 'w1infty_norm:', w1infty; check ( w1infty < error_tolerance )
    end if 
    call p_error_norm%free() 

  end subroutine check_solution

  subroutine check_time_solution(this)
    implicit none
    class(par_test_interpolators_driver_t), intent(inout) :: this
    class(vector_function_t), pointer :: H_exact_function
    class(scalar_function_t), pointer :: p_exact_function
    type(error_norms_vector_t) :: H_error_norm
    type(error_norms_scalar_t) :: p_error_norm 
    real(rp) :: mean, l1, l2, lp, linfty, h1, hcurl, h1_s, w1p_s, w1p, w1infty_s, w1infty
    real(rp) :: error_tolerance, time_

    H_exact_function => this%problem_functions%get_magnetic_field_solution()
    p_exact_function => this%problem_functions%get_pressure_solution() 

    error_tolerance = 1.0e-06 

    call H_error_norm%create(this%fe_space, MAGNETIC_FIELD_ID )

    
    mean = H_error_norm%compute(H_exact_function, this%time_solution, mean_norm, time=this%time)   
    l1 = H_error_norm%compute(H_exact_function, this%time_solution, l1_norm, time=this%time)   
    l2 = H_error_norm%compute(H_exact_function, this%time_solution, l2_norm, time=this%time)   
    lp = H_error_norm%compute(H_exact_function, this%time_solution, lp_norm, time=this%time)   
    linfty = H_error_norm%compute(H_exact_function, this%time_solution, linfty_norm, time=this%time)   
    h1_s = H_error_norm%compute(H_exact_function, this%time_solution, h1_seminorm, time=this%time) 
    h1 = H_error_norm%compute(H_exact_function, this%time_solution, h1_norm, time=this%time) 
    hcurl = H_error_norm%compute(H_exact_function, this%time_solution, hcurl_seminorm, time=this%time) 
    w1p_s = H_error_norm%compute(H_exact_function, this%time_solution, w1p_seminorm, time=this%time)   
    w1p = H_error_norm%compute(H_exact_function, this%time_solution, w1p_norm, time=this%time)   
    w1infty_s = H_error_norm%compute(H_exact_function, this%time_solution, w1infty_seminorm, time=this%time) 
    w1infty = H_error_norm%compute(H_exact_function, this%time_solution, w1infty_norm, time=this%time)

    if (this%par_environment%am_i_l1_root()) then
    write(*,*) 'PROJECTED TRANSIENT MAGNETIC FIELD FUNCTION ERROR NORMS *************'
    write(*,'(a20,e32.25)') 'mean_norm:', mean; check ( abs(mean) < error_tolerance )
    write(*,'(a20,e32.25)') 'l1_norm:', l1; check ( l1 < error_tolerance )
    write(*,'(a20,e32.25)') 'l2_norm:', l2; check ( l2 < error_tolerance )
    write(*,'(a20,e32.25)') 'lp_norm:', lp; check ( lp < error_tolerance )
    write(*,'(a20,e32.25)') 'linfnty_norm:', linfty; check ( linfty < error_tolerance )
    write(*,'(a20,e32.25)') 'h1_seminorm:', h1_s; check ( h1_s < error_tolerance )
    write(*,'(a20,e32.25)') 'h1_norm:', h1; check ( h1 < error_tolerance )
    write(*,'(a20,e32.25)') 'hcurl_norm:', hcurl; check ( hcurl < error_tolerance )
    write(*,'(a20,e32.25)') 'w1p_seminorm:', w1p_s; check ( w1p_s < error_tolerance )
    write(*,'(a20,e32.25)') 'w1p_norm:', w1p; check ( w1p < error_tolerance )
    write(*,'(a20,e32.25)') 'w1infty_seminorm:', w1infty_s; check ( w1infty_s < error_tolerance )
    write(*,'(a20,e32.25)') 'w1infty_norm:', w1infty; check ( w1infty < error_tolerance )
    end if 
    call H_error_norm%free()


    call p_error_norm%create(this%fe_space, PRESSURE_FIELD_ID )
    mean = p_error_norm%compute(p_exact_function, this%time_solution, mean_norm, time=this%time)   
    l1 = p_error_norm%compute(p_exact_function, this%time_solution, l1_norm, time=this%time)   
    l2 = p_error_norm%compute(p_exact_function, this%time_solution, l2_norm, time=this%time)   
    lp = p_error_norm%compute(p_exact_function, this%time_solution, lp_norm, time=this%time)   
    linfty = p_error_norm%compute(p_exact_function, this%time_solution, linfty_norm, time=this%time)   
    h1_s = p_error_norm%compute(p_exact_function, this%time_solution, h1_seminorm, time=this%time) 
    h1 = p_error_norm%compute(p_exact_function, this%time_solution, h1_norm, time=this%time) 
    w1p_s = p_error_norm%compute(p_exact_function, this%time_solution, w1p_seminorm, time=this%time)   
    w1p = p_error_norm%compute(p_exact_function, this%time_solution, w1p_norm, time=this%time)   
    w1infty_s = p_error_norm%compute(p_exact_function, this%time_solution, w1infty_seminorm, time=this%time) 
    w1infty = p_error_norm%compute(p_exact_function, this%time_solution, w1infty_norm, time=this%time)

    if (this%par_environment%am_i_l1_root()) then
    write(*,*) 'PROJECTED TRANSIENT SCALAR PRESSURE FUNCTION ERROR NORMS ************************'
    write(*,'(a20,e32.25)') 'mean_norm:', mean; check ( abs(mean) < error_tolerance )
    write(*,'(a20,e32.25)') 'l1_norm:', l1; check ( l1 < error_tolerance )
    write(*,'(a20,e32.25)') 'l2_norm:', l2; check ( l2 < error_tolerance )
    write(*,'(a20,e32.25)') 'lp_norm:', lp; check ( lp < error_tolerance )
    write(*,'(a20,e32.25)') 'linfnty_norm:', linfty; check ( linfty < error_tolerance )
    write(*,'(a20,e32.25)') 'h1_seminorm:', h1_s; check ( h1_s < error_tolerance )
    write(*,'(a20,e32.25)') 'h1_norm:', h1; check ( h1 < error_tolerance )
    write(*,'(a20,e32.25)') 'w1p_seminorm:', w1p_s; check ( w1p_s < error_tolerance )
    write(*,'(a20,e32.25)') 'w1p_norm:', w1p; check ( w1p < error_tolerance )
    write(*,'(a20,e32.25)') 'w1infty_seminorm:', w1infty_s; check ( w1infty_s < error_tolerance )
    write(*,'(a20,e32.25)') 'w1infty_norm:', w1infty; check ( w1infty < error_tolerance )
    end if 
    call p_error_norm%free() 

  end subroutine check_time_solution

  subroutine run_simulation(this) 
    implicit none
    class(par_test_interpolators_driver_t), intent(inout) :: this
    call this%free()
    call this%setup_triangulation()
    call this%setup_reference_fes()
    call this%setup_coarse_fe_handlers()
    call this%setup_fe_space()
    call this%interpolate_analytical_functions()
    call this%check_solution()
    call this%check_time_solution() 
    call this%free()
  end subroutine run_simulation

  subroutine free_command_line_parameters(this)
    implicit none
    class(par_test_interpolators_driver_t), intent(inout) :: this
    call this%test_params%free()
  end subroutine free_command_line_parameters

  subroutine free_environment(this)
    implicit none
    class(par_test_interpolators_driver_t), intent(inout) :: this
    call this%par_environment%free()
  end subroutine free_environment

  subroutine free(this)
    implicit none
    class(par_test_interpolators_driver_t), intent(inout) :: this
    integer(ip) :: i, istat

    call this%block_layout%free() 
    call this%solution%free() 
    call this%time_solution%free() 
    call this%fe_space%free()

    if ( this%par_environment%am_i_l1_task() ) then 
       if (allocated(this%coarse_fe_handlers))  then 
          deallocate(this%coarse_fe_handlers, stat=istat);  check(istat==0) 
       end if
    end if

    if ( allocated(this%reference_fes) ) then
       do i=1, size(this%reference_fes)
          if ( this%par_environment%am_i_l1_task() ) call this%reference_fes(i)%p%free()
       end do
       deallocate(this%reference_fes, stat=istat); check(istat==0)
    end if

    if ( this%par_environment%am_i_l1_task() ) call this%triangulation%free()
    if ( allocated( this%set_ids_to_reference_fes ) ) then 
       call memfree(this%set_ids_to_reference_fes, __FILE__, __LINE__ )
    end if

  end subroutine free

end module test_interpolators_driver_names
