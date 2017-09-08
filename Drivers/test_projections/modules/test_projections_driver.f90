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
module test_projections_driver_names
  use fempar_names
  use projections_params_names
  use projections_analytical_functions_names
		use projections_discrete_integration_names 
  use projections_conditions_names
# include "debug.i90"

  implicit none
  private

		integer(ip), parameter :: MAGNETIC_FIELD_ID = 1
  integer(ip), parameter :: PRESSURE_FIELD_ID = 2
		
  type test_projections_driver_t 
     private 

     ! Place-holder for parameter-value set provided through command-line interface
     type(projections_params_t)           :: test_params
     type(ParameterList_t)                :: parameter_list

     ! Cells and lower dimension objects container
     type(serial_triangulation_t)                :: triangulation

     ! Analytical functions of the problem
     type(projections_analytical_functions_t) :: problem_functions

     ! Discrete weak problem integration-related data type instances 
     type(serial_fe_space_t)                     :: fe_space 
     type(p_reference_fe_t) , allocatable        :: reference_fes(:) 
     type(projections_conditions_t)              :: projections_conditions
					type(fe_affine_operator_t)                  :: fe_affine_operator
					type(projections_discrete_integration_t)    :: projections_integration

     ! Problem solution FE function
     type(fe_function_t)                         :: solution

   contains
     procedure                  :: run_simulation
     procedure        , private :: parse_command_line_parameters
     procedure        , private :: setup_triangulation
     procedure        , private :: setup_reference_fes
     procedure        , private :: setup_fe_space
     procedure        , private :: project_analytical_function 
     procedure        , private :: check_solution
     procedure        , private :: write_solution
     procedure        , private :: free
  end type test_projections_driver_t

  ! Types
  public :: test_projections_driver_t

contains

  subroutine parse_command_line_parameters(this)
    implicit none
    class(test_projections_driver_t ), intent(inout) :: this
    call this%test_params%create()
    call this%test_params%parse(this%parameter_list)
  end subroutine parse_command_line_parameters

  subroutine setup_triangulation(this)
    implicit none
    class(test_projections_driver_t), intent(inout) :: this
    call this%triangulation%create(this%parameter_list)
  end subroutine setup_triangulation

  subroutine setup_reference_fes(this)
    implicit none
    class(test_projections_driver_t), intent(inout) :: this
    integer(ip) :: istat
    class(vef_iterator_t), allocatable  :: vef
    class(cell_iterator_t), allocatable       :: cell
    class(reference_fe_t), pointer :: reference_fe_geo


    allocate(this%reference_fes(2), stat=istat)
    check(istat==0)

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

  subroutine setup_fe_space(this)
    implicit none
    class(test_projections_driver_t), intent(inout) :: this 
				
    call this%projections_conditions%set_num_dims(this%triangulation%get_num_dims() + 1)
    call this%fe_space%create( triangulation       = this%triangulation, &
                               reference_fes       = this%reference_fes, &
                               conditions          = this%projections_conditions )
    call this%fe_space%set_up_cell_integration()
    call this%fe_space%set_up_facet_integration()
				
				call this%fe_affine_operator%create ( sparse_matrix_storage_format      = csr_format, &
                                          diagonal_blocks_symmetric_storage = [ .false.  ], &
                                          diagonal_blocks_symmetric         = [ .false. ], &
                                          diagonal_blocks_sign              = [ SPARSE_MATRIX_SIGN_UNKNOWN ], &
                                          fe_space                          = this%fe_space,           &
                                          discrete_integration              = this%projections_integration )
								
  end subroutine setup_fe_space

  subroutine project_analytical_function (this)
    implicit none
    class(test_projections_driver_t), intent(inout) :: this 
				class(vector_t) , pointer :: dof_values
									
    call this%problem_functions%set_num_dims(this%triangulation%get_num_dims())

				! Set boundary conditions 
	   call this%projections_conditions%set_boundary_function_Hx(this%problem_functions%get_boundary_function_Hx())
	   call this%projections_conditions%set_boundary_function_Hy(this%problem_functions%get_boundary_function_Hy())
	   if ( this%triangulation%get_num_dims() == 3) then 
	   call this%projections_conditions%set_boundary_function_Hz(this%problem_functions%get_boundary_function_Hz())
	   end if 
				call this%projections_conditions%set_boundary_function_pressure(this%problem_functions%get_boundary_function_pressure()) 

				! Project functions 
				call this%solution%create(this%fe_space)
				call this%fe_space%project_function( this%problem_functions%get_magnetic_field_solution(), this%solution , MAGNETIC_FIELD_ID) 
				call this%fe_space%project_function( this%problem_functions%get_pressure_solution(), this%solution , PRESSURE_FIELD_ID)
				call this%fe_space%project_Dirichlet_boundary_function( this%solution )
				
				!				WRITE(*,*) ' PROJECTED VALUES **************************************' 
				!		dof_values => this%solution%get_fixed_dof_values() 
				!			
				!select type (dof_values)
    !class is (serial_scalar_array_t)  
    !   call dof_values%print_matrix_market(6)
    !class DEFAULT
    !   assert(.false.) 
    !end select
				
  end subroutine project_analytical_function 

  subroutine check_solution(this)
    implicit none
    class(test_projections_driver_t), intent(inout) :: this
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
    write(*,*) 'PROJECTED MAGNETIC FIELD FUNCTION ERROR NORMS *************'
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
				call H_error_norm%free()
								
							
				call p_error_norm%create(this%fe_space, PRESSURE_FIELD_ID )
				write(*,*) 'PROJECTED SCALAR PRESSURE FUNCTION ERROR NORMS ************************'
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
   
				call p_error_norm%free() 
				
  end subroutine check_solution 
  
  subroutine write_solution(this)
    implicit none
    class(test_projections_driver_t), intent(in) :: this
    type(output_handler_t)                           :: oh
    if(this%test_params%get_write_solution()) then
        call oh%create()
        call oh%attach_fe_space(this%fe_space)
        call oh%add_fe_function(this%solution, MAGNETIC_FIELD_ID, 'H-solution')
								call oh%add_fe_function(this%solution, PRESSURE_FIELD_ID, 'p-solution')
        call oh%open(this%test_params%get_dir_path_out(), this%test_params%get_prefix())
        call oh%write()
        call oh%close()
        call oh%free()
    endif
  end subroutine write_solution

  subroutine run_simulation(this) 
    implicit none
    class(test_projections_driver_t), intent(inout) :: this
    call this%free()
    call this%parse_command_line_parameters()
    call this%setup_triangulation()
    call this%setup_reference_fes()
    call this%setup_fe_space()
    call this%project_analytical_function()
    call this%write_solution()
    call this%check_solution()
    call this%free()
  end subroutine run_simulation

  subroutine free(this)
    implicit none
    class(test_projections_driver_t), intent(inout) :: this
    integer(ip) :: i, istat
    call this%solution%free()   
    call this%fe_space%free()
    if ( allocated(this%reference_fes) ) then
       do i=1, size(this%reference_fes)
          call this%reference_fes(i)%p%free()
       end do
       deallocate(this%reference_fes, stat=istat); check(istat==0)
    end if
    call this%triangulation%free()
    call this%test_params%free()
				call this%fe_affine_operator%free() 
  end subroutine free

end module test_projections_driver_names
