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
module test_fe_cell_predicate_driver_names
  use fempar_names
 
# include "debug.i90"

  implicit none
  private

  type test_fe_cell_predicate_driver_t 
     private 
     type(serial_triangulation_t)                 :: triangulation
     type(p_reference_fe_t), allocatable          :: reference_fes(:) 
     type(serial_fe_space_t)                      :: fe_space 
     type(environment_t)                          :: serial_environment
     type(fe_cell_set_id_predicate_t)             :: fe_cell_predicate
   contains

     procedure                     :: run_simulation
     procedure                     :: setup_environment
     procedure                     :: free_environment
     procedure, private            :: setup_triangulation
     procedure, private            :: setup_reference_fe
     procedure, private            :: setup_fe_cell_predicate
     procedure, private            :: setup_fe_space
     procedure, private            :: check_fe_cell_predicate
     procedure, private            :: free
  end type test_fe_cell_predicate_driver_t

  ! Types
  public :: test_fe_cell_predicate_driver_t

contains

  subroutine run_simulation(this)
    implicit none
    class(test_fe_cell_predicate_driver_t), intent(inout) :: this
    call this%setup_triangulation()
    call this%setup_reference_fe()
    call this%setup_fe_cell_predicate()
    call this%setup_fe_space()
    call this%check_fe_cell_predicate
    call this%free()
  end subroutine run_simulation

  subroutine setup_environment(this, world_context)
    implicit none
    class(test_fe_cell_predicate_driver_t), intent(inout) :: this
    class(execution_context_t)            , intent(in)    :: world_context
    type(ParameterList_t)  :: parameter_list
    integer(ip) :: istat
    call parameter_list%init()
    istat = 0
    istat = istat + parameter_list%set(key = environment_type_key, value = structured)
    check(istat==0)
    call this%serial_environment%create(world_context, parameter_list)
    call parameter_list%free()
  end subroutine setup_environment
  
  subroutine free_environment(this)
    implicit none
    class(test_fe_cell_predicate_driver_t), intent(inout) :: this
    call this%serial_environment%free()
  end subroutine free_environment
  
  subroutine setup_triangulation(this)
    implicit none
    class(test_fe_cell_predicate_driver_t), intent(inout) :: this
    type(ParameterList_t)                                          :: parameter_list
    integer(ip)                                                    :: istat
    class(cell_iterator_t), allocatable                            :: cell
    call parameter_list%init()
    istat = 0
    istat = istat + parameter_list%set(key = num_dims_key, value = 2)
    istat = istat + parameter_list%set(key = num_cells_x_dir_key, value = [200,200,0])
    istat = istat + parameter_list%set(key = is_dir_periodic_key, value = [0,0,0])
    istat = istat + parameter_list%set(key = triangulation_generate_key, value = triangulation_generate_structured)
    check(istat==0)
    call this%triangulation%create(this%serial_environment, parameter_list)
    call parameter_list%free()
    call this%triangulation%create_cell_iterator(cell)
    do while ( .not. cell%has_finished() )
       if (mod(cell%get_gid(),2)==0) then
         call cell%set_set_id(2) 
       else
         call cell%set_set_id(1) 
       end if
       call cell%next()
    end do 
    call this%triangulation%free_cell_iterator(cell)
  end subroutine setup_triangulation

  subroutine setup_reference_fe(this)
    implicit none
    class(test_fe_cell_predicate_driver_t), intent(inout)   :: this
    class(cell_iterator_t),           allocatable    :: cell
    class(reference_fe_t),                pointer    :: reference_fe_geo
    integer(ip)                                      :: reference_fe_order
    integer(ip)                                      :: istat
    call this%triangulation%create_cell_iterator(cell)
    reference_fe_geo => cell%get_reference_fe()
    reference_fe_order = 1
    allocate(this%reference_fes(1), stat=istat); check(istat==0)
    this%reference_fes(1) =  make_reference_fe ( topology = reference_fe_geo%get_topology(), &
         fe_type = fe_type_lagrangian, &
         num_dims = this%triangulation%get_num_dims(), &
         order = reference_fe_order, &
         field_type = field_type_scalar, &
         conformity = .true. , &
         continuity = .true.)
    call this%triangulation%free_cell_iterator(cell)
  end subroutine setup_reference_fe

  subroutine setup_fe_cell_predicate(this)
    implicit none
    class(test_fe_cell_predicate_driver_t), intent(inout) :: this
       call this%fe_cell_predicate%set_cell_set_id(2)
  end subroutine setup_fe_cell_predicate  
  
  subroutine setup_fe_space(this)
    implicit none
    class(test_fe_cell_predicate_driver_t), intent(inout)   :: this
    integer(ip)                                      :: set_ids_to_reference_fes
    call this%fe_space%create( triangulation = this%triangulation, &
                               reference_fes = this%reference_fes)
    call this%fe_space%set_up_cell_integration()
  end subroutine setup_fe_space

  subroutine check_fe_cell_predicate(this)
    implicit none
    class(test_fe_cell_predicate_driver_t), intent(inout)   :: this
    class(fe_cell_iterator_t), allocatable                           :: fe
    call this%fe_space%create_fe_cell_iterator(fe,this%fe_cell_predicate)
    do while ( .not. fe%has_finished() )
       check(fe%get_set_id() == 2)
       call fe%next()
    end do 
    call this%fe_space%free_fe_cell_iterator(fe)
  end subroutine check_fe_cell_predicate  
  
  subroutine free(this)
    implicit none
    class(test_fe_cell_predicate_driver_t), intent(inout) :: this
    call this%fe_space%free()
    call this%reference_fes(1)%p%free()
    call this%triangulation%free()
  end subroutine free
  

end module test_fe_cell_predicate_driver_names
