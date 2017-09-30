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
module par_pb_bddc_linear_elasticity_driver_names
  use fempar_names
  use par_pb_bddc_linear_elasticity_params_names

  use linear_elasticity_discrete_integration_names
  use linear_elasticity_conditions_names

  use linear_elasticity_analytical_functions_names
# include "debug.i90"

  implicit none
  private

  type par_pb_bddc_linear_elasticity_fe_driver_t 
     private 

     ! Place-holder for parameter-value set provided through command-line interface
     type(par_pb_bddc_linear_elasticity_params_t)      :: test_params
     type(ParameterList_t), pointer       :: parameter_list

     ! Cells and lower dimension objects container
     type(par_triangulation_t)             :: triangulation
     integer(ip), allocatable              :: cell_set_ids(:)

     ! Discrete weak problem integration-related data type instances 
     type(par_fe_space_t)                        :: fe_space 
     type(p_reference_fe_t), allocatable         :: reference_fes(:) 
     type(standard_l1_coarse_fe_handler_t)       :: standard_coarse_fe_handler
!!! Need to change this
     type(vector_laplacian_pb_bddc_l1_coarse_fe_handler_t):: vector_laplacian_coarse_fe_handler
     type(p_l1_coarse_fe_handler_t), allocatable :: coarse_fe_handlers(:)
!!! Need to change this
     type(irreducible_discrete_integration_t) :: linear_elasticity_integration
     type(linear_elasticity_conditions_t)           :: linear_elasticity_conditions
     type(linear_elasticity_analytical_functions_t) :: linear_elasticity_analytical_functions


     ! Place-holder for the coefficient matrix and RHS of the linear system
     type(fe_affine_operator_t)            :: fe_affine_operator

     ! MLBDDC preconditioner
     type(mlbddc_t)                            :: mlbddc

     ! Iterative linear solvers data type
     type(iterative_linear_solver_t)           :: iterative_linear_solver

     ! Poisson problem solution FE function
     type(fe_function_t)                   :: solution
     type(constant_vector_function_t)      :: zero_vector

     ! Environment required for fe_affine_operator + vtk_handler
     type(environment_t)                    :: par_environment

     ! Timers
     type(timer_t) :: timer_triangulation
     type(timer_t) :: timer_fe_space
     type(timer_t) :: timer_assemply
     type(timer_t) :: timer_solver_setup
     type(timer_t) :: timer_solver_run

   contains
     procedure                  :: parse_command_line_parameters
     procedure                  :: setup_timers
     procedure                  :: report_timers
     procedure                  :: free_timers
     procedure                  :: setup_environment
     procedure        , private :: setup_triangulation
     procedure        , private :: setup_cell_set_ids
     procedure        , private :: setup_discrete_integration     
     procedure        , private :: setup_reference_fes
     procedure        , private :: setup_coarse_fe_handlers
     procedure        , private :: setup_fe_space
     procedure        , private :: setup_system
     procedure        , private :: setup_solver
     procedure        , private :: assemble_system
     procedure        , private :: solve_system
     procedure        , private :: check_solution
     procedure        , private :: write_solution
     procedure                  :: run_simulation
     procedure        , private :: free
     procedure                  :: free_command_line_parameters
     procedure                  :: free_environment
  end type par_pb_bddc_linear_elasticity_fe_driver_t

  ! Types
  public :: par_pb_bddc_linear_elasticity_fe_driver_t

contains

  subroutine parse_command_line_parameters(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_fe_driver_t), intent(inout) :: this
    call this%test_params%create()
    this%parameter_list => this%test_params%get_values()
  end subroutine parse_command_line_parameters

  !========================================================================================
  subroutine setup_timers(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_fe_driver_t), intent(inout) :: this
    class(execution_context_t), pointer :: w_context
    w_context => this%par_environment%get_w_context()
    call this%timer_triangulation%create(w_context,"SETUP TRIANGULATION")
    call this%timer_fe_space%create(     w_context,"SETUP FE SPACE")
    call this%timer_assemply%create(     w_context,"FE INTEGRATION AND ASSEMBLY")
    call this%timer_solver_setup%create( w_context,"SETUP SOLVER AND PRECONDITIONER")
    call this%timer_solver_run%create(   w_context,"SOLVER RUN")
  end subroutine setup_timers

  !========================================================================================
  subroutine report_timers(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_fe_driver_t), intent(inout) :: this
    call this%timer_triangulation%report(.true.)
    call this%timer_fe_space%report(.false.)
    call this%timer_assemply%report(.false.)
    call this%timer_solver_setup%report(.false.)
    call this%timer_solver_run%report(.false.)
    if (this%par_environment%get_l1_rank() == 0) then
       write(*,*)
    end if
  end subroutine report_timers

  !========================================================================================
  subroutine free_timers(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_fe_driver_t), intent(inout) :: this
    call this%timer_triangulation%free()
    call this%timer_fe_space%free()
    call this%timer_assemply%free()
    call this%timer_solver_setup%free()
    call this%timer_solver_run%free()
  end subroutine free_timers

  !========================================================================================
  subroutine setup_environment(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_fe_driver_t), intent(inout) :: this
    integer(ip) :: istat
    if ( this%test_params%get_triangulation_type() == triangulation_generate_structured ) then
       istat = this%parameter_list%set(key = environment_type_key, value = structured) ; check(istat==0)
    else
       istat = this%parameter_list%set(key = environment_type_key, value = unstructured) ; check(istat==0)
    end if
    istat = this%parameter_list%set(key = execution_context_key, value = mpi_context) ; check(istat==0)
    call this%par_environment%create (this%parameter_list)
  end subroutine setup_environment

  subroutine setup_triangulation(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_fe_driver_t), intent(inout) :: this
    class(vef_iterator_t), allocatable  :: vef
!!! Do we really need this  
    logical :: fixed_pressure


    call this%triangulation%create(this%parameter_list, this%par_environment)

    if ( this%test_params%get_triangulation_type() == triangulation_generate_structured ) then
       fixed_pressure = .false.
       call this%triangulation%create_vef_iterator(vef)
       do while ( .not. vef%has_finished() )
          if(vef%is_at_boundary()) then
             call vef%set_set_id(1)
!!! Do we really need this
             if(.not.fixed_pressure) then
                fixed_pressure = .true.
                call vef%set_set_id(2)
             end if
          else
             call vef%set_set_id(0)
          end if
          call vef%next()
       end do
       call this%triangulation%free_vef_iterator(vef)
    end if


    if ( this%test_params%get_coarse_fe_handler_type() == pb_bddc ) then
       call this%setup_cell_set_ids() 
    end if
    call this%triangulation%setup_coarse_triangulation()
    !write(*,*) 'CG: NUMBER OBJECTS', this%triangulation%get_num_objects()
    if ( this%test_params%get_coarse_fe_handler_type() == standard_bddc ) then
       call this%setup_cell_set_ids() 
    end if

  end subroutine setup_triangulation

  subroutine setup_cell_set_ids(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_fe_driver_t), intent(inout) :: this
    class(cell_iterator_t), allocatable       :: cell
    integer(ip)                               :: istat, idummy
    integer(ip), allocatable                  :: cells_set(:)
    type(point_t), allocatable :: cell_coords(:)
    type(point_t) :: grav_center
    integer(ip)   :: inode, l1_rank  
    real(rp), allocatable:: px1(:), px2(:), py1(:), py2(:),  pz1(:), pz2(:)

!!! Need to change this
    !!this%linear_elasticity_integration%diffusion_inclusion = this%test_params%get_jump()    
    !!this%vector_laplacian_coarse_fe_handler%diffusion_inclusion = this%test_params%get_jump()

    if ( this%par_environment%am_i_l1_task() ) then
       l1_rank = this%par_environment%get_l1_rank()
       call memalloc( this%triangulation%get_num_local_cells(), cells_set, __FILE__, __LINE__ ) 
       call this%triangulation%create_cell_iterator(cell)
       allocate (cell_coords(cell%get_num_nodes()),stat=istat)
       do while( .not. cell%has_finished() )
          if ( cell%is_local() ) then
             call cell%get_nodes_coordinates(cell_coords)
             call grav_center%init(0.0_rp)
             do inode = 1, cell%get_num_nodes()
                grav_center = grav_center + cell_coords(inode)
             end do
             grav_center = (1.0_rp/cell%get_num_nodes())*grav_center
             cells_set( cell%get_gid() ) = cell_set_id( grav_center, &
                  this%triangulation%get_num_dims(), &
                  this%test_params%get_jump(), this%test_params%get_inclusion(), &
                  this%test_params%get_nchannel_x_direction(), &
                  this%test_params%get_nparts_with_channels(), &
                  this%test_params%get_nparts())
          end if
          call cell%next()
       end do
       call this%triangulation%fill_cells_set( cells_set )
       call memfree( cells_set, __FILE__, __LINE__ )
       call this%triangulation%free_cell_iterator(cell)
    end if

    if ( allocated(px1) ) call memfree( px1, __FILE__, __LINE__ )
    if ( allocated(px2) ) call memfree( px2, __FILE__, __LINE__ )
    if ( allocated(py1) ) call memfree( py1, __FILE__, __LINE__ )
    if ( allocated(py2) ) call memfree( py2, __FILE__, __LINE__ )
    if ( allocated(pz1) ) call memfree( pz1, __FILE__, __LINE__ )
    if ( allocated(pz2) ) call memfree( pz2, __FILE__, __LINE__ )

  contains

    function cell_set_id( coord, num_dims, jump, inclusion, nchannel_x_direction, nparts_with_channels,nparts)
      implicit none
      type(point_t), intent(in)  :: coord
      integer(ip)  , intent(in)  :: num_dims
      integer(ip)  , intent(in)  :: jump
      integer(ip)  , intent(in)  :: inclusion
      integer(ip)  , intent(in)  :: nchannel_x_direction(3)
      integer(ip)  , intent(in)  :: nparts_with_channels(3)
      integer(ip)  , intent(in)  :: nparts(3)
      type(point_t) :: origin, opposite
      integer(ip) :: cell_set_id
      integer(ip) :: i,j,k, nchannel, nchannel_in_each_direction
      real(rp)    :: y_pos_0, y_pos_1, z_pos_0, z_pos_1
      real(rp)    :: box_width, half_channel_width, center, eps
      real(rp) :: p1(6), p2(6), p1_b(4), p2_b(4)
      real(rp) :: p1_c(256), p2_c(256)


      cell_set_id = 1
      assert(nparts_with_channels(1)<=nparts(1))
      assert(nparts_with_channels(2)<=nparts(2))
      assert(nparts_with_channels(3)<=nparts(3))

      ! Consider one channel : [0,1], [0.25,0.5], [0.25,0.5]
      if ( inclusion == 1 ) then
         call origin%set(1,0.0_rp)  ; call origin%set(2, 0.25_rp) ; call origin%set(3,0.25_rp);
         call opposite%set(1,0.5_rp); call opposite%set(2,0.50_rp); call opposite%set(3,0.50_rp);
         if ( is_point_in_rectangle( origin, opposite, coord, num_dims ) ) cell_set_id = jump
      else if ( inclusion == 2 ) then
         nchannel = 8
         do i = 1, nchannel
            j = 2*i
            y_pos_0 = real(j-1)/real(2*nchannel+1)
            y_pos_1 = j/real(2*nchannel+1)
            z_pos_0 = (j-1)/real(2*nchannel+1)
            z_pos_1 = j/real(2*nchannel+1)
            call origin%set(1,y_pos_0)  ; call origin%set(2, y_pos_0) ; call origin%set(3,z_pos_0);
            call opposite%set(1,1.0_rp); call opposite%set(2,y_pos_1); call opposite%set(3,z_pos_1);
            if ( is_point_in_rectangle( origin, opposite, coord, num_dims ) ) cell_set_id = jump
            ! 
            z_pos_0 = j/real(2*nchannel+1)
            z_pos_1 = (j+1)/real(2*nchannel+1)
            call origin%set(1,y_pos_0)  ; call origin%set(2, y_pos_0) ; call origin%set(3,z_pos_0);
            call opposite%set(1,y_pos_1); call opposite%set(2,1.0_rp); call opposite%set(3,z_pos_1);
            if ( is_point_in_rectangle( origin, opposite, coord, num_dims ) ) cell_set_id = jump
         end do
      else if ( inclusion == 3 ) then
         ! Hieu's test in PB-BDDC article (two channels)
         p1 = [2.0_rp/36.0_rp, 7.0_rp/36.0_rp, 14.0_rp/36.0_rp, 19.0_rp/36.0_rp, 26.0_rp/36.0_rp, 30.0_rp/36.0_rp] ! lower y value
         p2 = [4.0_rp/36.0_rp, 9.0_rp/36.0_rp, 16.0_rp/36.0_rp, 21.0_rp/36.0_rp, 28.0_rp/36.0_rp, 32.0_rp/36.0_rp] ! upper y value
         do i = 1, 6
            call origin%set(1,0.0_rp)  ; call origin%set(2, p1(i)) ; call origin%set(3,p1(7-i));
            call opposite%set(1,1.0_rp); call opposite%set(2,p2(i)); call opposite%set(3,p2(7-i));
            if ( is_point_in_rectangle( origin, opposite, coord, num_dims ) ) cell_set_id = jump + i - 1
         end do
         do i = 7, 12  
            call origin%set(2,0.0_rp)  ; call origin%set(1, p1(i-6)) ; call origin%set(3,p1(i-6));
            call opposite%set(2,1.0_rp); call opposite%set(1,p2(i-6)); call opposite%set(3,p2(i-6));
            if ( is_point_in_rectangle( origin, opposite, coord, num_dims ) ) cell_set_id = jump + i - 1
         end do
      else if ( inclusion == 4 ) then
         ! Hieu's test in PB-BDDC article (two channels)
         p1_b = [4.0_rp/32.0_rp, 12.0_rp/32.0_rp, 20.0_rp/32.0_rp, 28.0_rp/32.0_rp] ! lower y value
         p2_b = [6.0_rp/32.0_rp, 14.0_rp/32.0_rp, 22.0_rp/32.0_rp, 30.0_rp/32.0_rp] ! upper y value
         nchannel = 1
         ! x edges
         do j = 1, 4
            do k = 1,4
               call origin%set(1,0.0_rp)  ; call origin%set(2, p1_b(j)) ; call origin%set(3,p1_b(k));
               call opposite%set(1,1.0_rp); call opposite%set(2,p2_b(j)); call opposite%set(3,p2_b(k));
               nchannel = nchannel + 1
               if ( is_point_in_rectangle( origin, opposite, coord, num_dims ) ) cell_set_id = nchannel
            end do
         end do
         ! y edges
         do j = 1, 4
            do k = 1,4
               call origin%set(2,0.0_rp)  ; call origin%set(1, p1_b(j)) ; call origin%set(3,p1_b(k));
               call opposite%set(2,1.0_rp); call opposite%set(1,p2_b(j)); call opposite%set(3,p2_b(k));
               nchannel = nchannel + 1
               if ( is_point_in_rectangle( origin, opposite, coord, num_dims ) ) cell_set_id = nchannel
            end do
         end do
         ! z edges
         do j = 1, 4
            do k = 1,4
               call origin%set(3,0.0_rp)  ; call origin%set(2, p1_b(j)) ; call origin%set(1,p1_b(k));
               call opposite%set(3,1.0_rp); call opposite%set(2,p2_b(j)); call opposite%set(1,p2_b(k));
               nchannel = nchannel + 1
               if ( is_point_in_rectangle( origin, opposite, coord, num_dims ) ) cell_set_id = nchannel
            end do
         end do
      else if ( inclusion == 5 ) then
         ! Number of channels can be choosen from the command line using option -nc

         ! defining positions of the channels
         box_width = 1.0_rp/nchannel_x_direction(1)
         half_channel_width = box_width/5
         center = box_width/2
         eps = 1e-14_rp
         do j=1, nchannel_x_direction(1)
            p1_c(j)=center - half_channel_width 
            p2_c(j)=center + half_channel_width  
            center = center + box_width
         enddo

         nchannel = 1
         ! x edges
         do j = 1, nchannel_x_direction(1)
            do k = 1,nchannel_x_direction(1)
               call origin%set(1,0.0_rp)  ; call origin%set(2, p1_c(j)) ; call origin%set(3,p1_c(k));
               call opposite%set(1,1.0_rp); call opposite%set(2,p2_c(j)); call opposite%set(3,p2_c(k));
               nchannel = nchannel + 1
               if ( is_point_in_rectangle( origin, opposite, coord, num_dims ) ) cell_set_id = nchannel
            end do
         end do
         ! y edges
         do j = 1, nchannel_x_direction(1)
            do k = 1,nchannel_x_direction(1)
               call origin%set(2,0.0_rp)  ; call origin%set(1, p1_c(j)) ; call origin%set(3,p1_c(k));
               call opposite%set(2,1.0_rp); call opposite%set(1,p2_c(j)); call opposite%set(3,p2_c(k));
               nchannel = nchannel + 1
               if ( is_point_in_rectangle( origin, opposite, coord, num_dims ) ) cell_set_id = nchannel
            end do
         end do
         ! z edges
         do j = 1, nchannel_x_direction(1)
            do k = 1,nchannel_x_direction(1)
               call origin%set(3,0.0_rp)  ; call origin%set(2, p1_c(j)) ; call origin%set(1,p1_c(k));
               call opposite%set(3,1.0_rp); call opposite%set(2,p2_c(j)); call opposite%set(1,p2_c(k));
               nchannel = nchannel + 1
               if ( is_point_in_rectangle( origin, opposite, coord, num_dims ) ) cell_set_id = nchannel
            end do
         end do
      else if ( inclusion == 6 ) then
         ! Number of channels can be choosen from the command line using option -nc
         ! Channels are positioned so that some will touch but not cross the interface
         ! if the nchannel_x_direction is a multiple of the number partitions per direction 

         ! defining positions of the channels
         box_width = 1.0_rp/nchannel_x_direction(1)
         half_channel_width = box_width/5
         center = half_channel_width
         eps = 1e-14_rp
         do j=1, nchannel_x_direction(1)
            p1_c(j)=center - half_channel_width 
            p2_c(j)=center + half_channel_width  
            center = center + box_width
         enddo

         nchannel = 1
         ! x edges
         do j = 1, nchannel_x_direction(1)
            do k = 1,nchannel_x_direction(1)
               call origin%set(1,0.0_rp)  ; call origin%set(2, p1_c(j)) ; call origin%set(3,p1_c(k));
               call opposite%set(1,1.0_rp); call opposite%set(2,p2_c(j)); call opposite%set(3,p2_c(k));
               nchannel = nchannel + 1
               if ( is_point_in_rectangle( origin, opposite, coord, num_dims ) ) cell_set_id = nchannel
            end do
         end do
         ! y edges
         do j = 1, nchannel_x_direction(1)
            do k = 1,nchannel_x_direction(1)
               call origin%set(2,0.0_rp)  ; call origin%set(1, p1_c(j)) ; call origin%set(3,p1_c(k));
               call opposite%set(2,1.0_rp); call opposite%set(1,p2_c(j)); call opposite%set(3,p2_c(k));
               nchannel = nchannel + 1
               if ( is_point_in_rectangle( origin, opposite, coord, num_dims ) ) cell_set_id = nchannel
            end do
         end do
         ! z edges
         do j = 1, nchannel_x_direction(1)
            do k = 1,nchannel_x_direction(1)
               call origin%set(3,0.0_rp)  ; call origin%set(2, p1_c(j)) ; call origin%set(1,p1_c(k));
               call opposite%set(3,1.0_rp); call opposite%set(2,p2_c(j)); call opposite%set(1,p2_c(k));
               nchannel = nchannel + 1
               if ( is_point_in_rectangle( origin, opposite, coord, num_dims ) ) cell_set_id = nchannel
            end do
         end do
      else if ( inclusion == 7 ) then
         ! Number of channels can be choosen from the command line using option -nc
         ! Channels are positioned so that some will touch but not cross the interface
         ! if the nchannel_x_direction is a multiple of the number partitions per direction 

         ! defining positions of the channels
         if (.not.(allocated(px1))) then
            call memalloc(nchannel_x_direction(1), px1, __FILE__, __LINE__ )
            call memalloc(nchannel_x_direction(1), px2, __FILE__, __LINE__ )
            box_width = 1.0*nparts_with_channels(1)/(nchannel_x_direction(1)*nparts(1))
            half_channel_width = box_width/5
            center = half_channel_width
            do j=1, nchannel_x_direction(1)
               px1(j)=center - half_channel_width 
               px2(j)=center + half_channel_width  
               center = center + box_width
            enddo
            call memalloc(nchannel_x_direction(2), py1, __FILE__, __LINE__ )
            call memalloc(nchannel_x_direction(2), py2, __FILE__, __LINE__ )
            box_width = 1.0*nparts_with_channels(2)/(nchannel_x_direction(2)*nparts(2))
            half_channel_width = box_width/5
            center = half_channel_width
            do j=1, nchannel_x_direction(2)
               py1(j)=center - half_channel_width 
               py2(j)=center + half_channel_width  
               center = center + box_width
            enddo
            call memalloc(nchannel_x_direction(3), pz1, __FILE__, __LINE__ )
            call memalloc(nchannel_x_direction(3), pz2, __FILE__, __LINE__ )
            box_width = 1.0*nparts_with_channels(3)/(nchannel_x_direction(3)*nparts(3))
            half_channel_width = box_width/5
            center = half_channel_width
            do j=1, nchannel_x_direction(3)
               pz1(j)=center - half_channel_width 
               pz2(j)=center + half_channel_width  
               center = center + box_width
            enddo
         end if

         nchannel = 1
         ! x edges
         do j = 1, nchannel_x_direction(2)
            do k = 1,nchannel_x_direction(3)
               call origin%set(1,0.0_rp)  ; call origin%set(2, py1(j)) ; call origin%set(3,pz1(k));
               call opposite%set(1,1.0_rp); call opposite%set(2,py2(j)); call opposite%set(3,pz2(k));
               nchannel = nchannel + 1
               if ( is_point_in_rectangle( origin, opposite, coord, num_dims ) ) cell_set_id = nchannel
            end do
         end do
         ! y edges
         do j = 1, nchannel_x_direction(1)
            do k = 1,nchannel_x_direction(3)
               call origin%set(2,0.0_rp)  ; call origin%set(1, px1(j)) ; call origin%set(3,pz1(k));
               call opposite%set(2,1.0_rp); call opposite%set(1,px2(j)); call opposite%set(3,pz2(k));
               nchannel = nchannel + 1
               if ( is_point_in_rectangle( origin, opposite, coord, num_dims ) ) cell_set_id = nchannel
            end do
         end do
         ! z edges
         do j = 1, nchannel_x_direction(1)
            do k = 1,nchannel_x_direction(2)
               call origin%set(3,0.0_rp)  ; call origin%set(2, px1(j)) ; call origin%set(1,py1(k));
               call opposite%set(3,1.0_rp); call opposite%set(2,px2(j)); call opposite%set(1,py2(k));
               nchannel = nchannel + 1
               if ( is_point_in_rectangle( origin, opposite, coord, num_dims ) ) cell_set_id = nchannel
            end do
         end do
         if ( cell_set_id /= 1 ) cell_set_id = 2
      else if ( inclusion == 8 ) then
         i = mod(mod(l1_rank,nparts(1)),2)
         j = mod(mod(l1_rank,nparts(1)*nparts(2))/nparts(1),2)
         k = mod(l1_rank/(nparts(1)*nparts(2)),2)
         if (((i==0) .and. (j==0) .and. (k==0)).or.&
              (i==0 .and. j==1 .and. k==0).or.&
              (i==1 .and. j==0 .and. k==0).or.&
              ((i==0) .and. (j==0) .and. (k==1))) then
            cell_set_id=2
         else 
            cell_set_id=1
         end if

      end if

    end function cell_set_id

    function is_point_in_rectangle( origin, opposite, coord, num_dims )
      implicit none
      type(point_t), intent(in)  :: origin
      type(point_t), intent(in)  :: opposite
      type(point_t), intent(in)  :: coord
      integer(ip)  , intent(in)  :: num_dims
      logical :: is_point_in_rectangle
      integer(ip) :: i
      is_point_in_rectangle = .true.
      do i = 1, num_dims
         if ( coord%get(i) < origin%get(i) .or. coord%get(i) > opposite%get(i) ) then
            is_point_in_rectangle = .false.
            exit
         end if
      end do
    end function is_point_in_rectangle

  end subroutine setup_cell_set_ids

  !========================================================================================

  subroutine setup_discrete_integration(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_fe_driver_t), intent(inout) :: this
    integer(ip) :: istat

    !allocate(irreducible_discrete_integration_t :: this%linear_elasticity_integration, stat=istat); check(istat==0)
    call this%linear_elasticity_integration%create(this%triangulation%get_num_dims(),this%linear_elasticity_analytical_functions)

  end subroutine setup_discrete_integration

  !========================================================================================

  subroutine setup_reference_fes(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_fe_driver_t), intent(inout) :: this
    integer(ip) :: istat
    class(cell_iterator_t), allocatable       :: cell
    class(reference_fe_t), pointer :: reference_fe_geo

    allocate(this%reference_fes(this%linear_elasticity_integration%get_number_fields()), stat=istat)
    check(istat==0)

    if ( this%par_environment%am_i_l1_task() ) then
       call this%triangulation%create_cell_iterator(cell)
       reference_fe_geo => cell%get_reference_fe()
       this%reference_fes(1) =  make_reference_fe ( topology = reference_fe_geo%get_topology(), &
            fe_type = fe_type_lagrangian, &
            num_dims = this%triangulation%get_num_dims(), &
            order = this%test_params%get_reference_fe_order(), &
            field_type = field_type_vector, &
            conformity = .true. )
       call this%triangulation%free_cell_iterator(cell)
    end if
  end subroutine setup_reference_fes

  subroutine setup_coarse_fe_handlers(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_fe_driver_t), target, intent(inout) :: this
    integer(ip) :: istat
    allocate(this%coarse_fe_handlers(1), stat=istat)
    check(istat==0)

    if ( this%test_params%get_coarse_fe_handler_type() == pb_bddc ) then
       this%coarse_fe_handlers(1)%p => this%vector_laplacian_coarse_fe_handler
    else if (this%test_params%get_coarse_fe_handler_type() == standard_bddc) then
       this%coarse_fe_handlers(1)%p => this%standard_coarse_fe_handler
    end if
  end subroutine setup_coarse_fe_handlers

  subroutine setup_fe_space(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_fe_driver_t), intent(inout) :: this

    integer(ip) :: set_ids_to_reference_fes(1,2)

    call this%linear_elasticity_analytical_functions%set_num_dimensions(this%triangulation%get_num_dims())
    call this%linear_elasticity_conditions%set_number_components(this%linear_elasticity_integration%get_number_components())
    call this%linear_elasticity_conditions%set_number_dimensions(this%triangulation%get_num_dims())    
    call this%linear_elasticity_conditions%set_boundary_function(this%linear_elasticity_analytical_functions%get_solution_function_u())

    call this%fe_space%create( triangulation       = this%triangulation,      &
         reference_fes       = this%reference_fes,      &
         coarse_fe_handlers  = this%coarse_fe_handlers, &
         conditions          = this%linear_elasticity_conditions )


    call this%fe_space%set_up_cell_integration()
    call this%fe_space%set_up_facet_integration()
    !call this%fe_space%print()
  end subroutine setup_fe_space

  subroutine setup_system (this)
    implicit none
    class(par_pb_bddc_linear_elasticity_fe_driver_t), intent(inout) :: this
    type(vector_field_t) :: zero_vector_field

!!! Do we really need this
    !call this%linear_elasticity_integration%set_source_term(this%linear_elasticity_analytical_functions%get_source_term())

    ! if (test_single_scalar_valued_reference_fe) then
    call this%fe_affine_operator%create ( sparse_matrix_storage_format      = csr_format, &
         diagonal_blocks_symmetric_storage = [ .true. ], &
         diagonal_blocks_symmetric         = [ .true. ], &
         diagonal_blocks_sign              = [ SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE ], &
         fe_space                          = this%fe_space, &
         discrete_integration              = this%linear_elasticity_integration, &
         field_blocks                      = this%linear_elasticity_integration%get_field_blocks(), &
         field_coupling                    = this%linear_elasticity_integration%get_field_coupling()) 

    call this%solution%create(this%fe_space) 
    call zero_vector_field%init(0.0_rp)
    call this%zero_vector%create(zero_vector_field)
    call this%fe_space%interpolate(1, this%zero_vector, this%solution)
    call this%fe_space%interpolate_dirichlet_values(this%solution)
    call this%linear_elasticity_integration%set_fe_function(this%solution)

  end subroutine setup_system

  subroutine setup_solver (this)
    implicit none

    class(par_pb_bddc_linear_elasticity_fe_driver_t), intent(inout) :: this
    type(parameterlist_t) :: parameter_list
    type(parameterlist_t), pointer :: plist, dirichlet, neumann, coarse
    integer(ip) :: FPLError
    integer(ip) :: ilev
    integer(ip) :: iparm(64)
    logical, parameter :: si_solver = .false.
    logical, parameter :: use_bddc_preconditioner = .true.

    if (.not.(use_bddc_preconditioner)) then
       call this%iterative_linear_solver%create(this%fe_space%get_environment())
       call this%iterative_linear_solver%set_type_from_string(cg_name)
       call parameter_list%init()
       FPLError = parameter_list%set(key = ils_rtol, value = 1.0e-9_rp)
       assert(FPLError == 0)
       FPLError = parameter_list%set(key = ils_max_num_iterations, value = 5000)
       assert(FPLError == 0)
       call this%iterative_linear_solver%set_parameters_from_pl(parameter_list)
       call this%iterative_linear_solver%set_operators(this%fe_affine_operator, .identity. this%fe_affine_operator) 
       call parameter_list%free()
       return
    end if

    if ( this%par_environment%get_l1_rank() == 0 ) then
       if (si_solver) then
          write(*,*) "si_solver:: 1"
       else 
          write(*,*) "si_solver:: 0"
       end if
    end if
    call this%fe_space%setup_coarse_fe_space(this%parameter_list)


    ! Prepare the internal parameter list of pardiso
    ! See https://software.intel.com/en-us/node/470298 for details
    iparm      = 0 ! Init all entries to zero
    ! Default values
    iparm(1)   = 1 ! no solver default
    iparm(2)   = 2 ! fill-in reordering from METIS
    iparm(8)   = 2 ! numbers of iterative refinement steps
    iparm(10)  = 8 ! perturb the pivot elements with 1E-8
    iparm(21)  = 1 ! 1x1 + 2x2 pivots
    ! Customization
    iparm(11)  = 1 ! use scaling (default 0)
    iparm(13)  = 1 ! use maximum weighted matching algorithm (default 0)

    ! Fill the fempar list to be eventually given to bddc
    plist => this%parameter_list 
    if ( this%par_environment%get_l1_size() == 1 ) then
       FPLError = plist%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
       if ( si_solver ) then
          !FPLError = plist%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_sin); assert(FPLError == 0)
          !FPLError = plist%set(key=pardiso_mkl_message_level, value=0); assert(FPLError == 0)
          !FPLError = plist%set(key=pardiso_mkl_iparm, value=iparm); assert(FPLError == 0)
       end if
    end if
    do ilev=1, this%par_environment%get_num_levels()-1
       ! Set current level Dirichlet solver parameters
       dirichlet => plist%NewSubList(key=mlbddc_dirichlet_solver_params)
       FPLError = dirichlet%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
       if ( si_solver ) then
          !FPLError = dirichlet%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_sin); assert(FPLError == 0)
          !FPLError = dirichlet%set(key=pardiso_mkl_message_level, value=0); assert(FPLError == 0)
          !FPLError = dirichlet%set(key=pardiso_mkl_iparm, value=iparm); assert(FPLError == 0)
       end if

       ! Set current level Neumann solver parameters
       neumann => plist%NewSubList(key=mlbddc_neumann_solver_params)
       FPLError = neumann%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
       if ( si_solver ) then
          FPLError = neumann%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_sin); assert(FPLError == 0)
          FPLError = neumann%set(key=pardiso_mkl_message_level, value=0); assert(FPLError == 0)
          FPLError = neumann%set(key=pardiso_mkl_iparm, value=iparm); assert(FPLError == 0)
       end if

       coarse => plist%NewSubList(key=mlbddc_coarse_solver_params) 
       plist  => coarse 
    end do
    ! Set coarsest-grid solver parameters
    FPLError = coarse%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
    if ( si_solver ) then
       !FPLError = coarse%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_sin); assert(FPLError == 0)
       !FPLError = coarse%set(key=pardiso_mkl_message_level, value=0); assert(FPLError == 0)
       !FPLError = coarse%set(key=pardiso_mkl_iparm, value=iparm); assert(FPLError == 0)
    end if


    ! Set-up MLBDDC preconditioner
    call this%mlbddc%create(this%fe_affine_operator, this%parameter_list)
    call this%mlbddc%symbolic_setup()
    call this%mlbddc%numerical_setup()

    call parameter_list%init()
    FPLError = parameter_list%set(key = ils_rtol, value = 1.0e-9_rp)
    FPLError = parameter_list%set(key = ils_max_num_iterations, value = 5000)
    assert(FPLError == 0)

    call this%iterative_linear_solver%create(this%fe_space%get_environment())
    call this%iterative_linear_solver%set_type_from_string(cg_name)
    call this%iterative_linear_solver%set_parameters_from_pl(parameter_list)
    call this%iterative_linear_solver%set_operators(this%fe_affine_operator, this%mlbddc) 
    call parameter_list%free()
  end subroutine setup_solver


  subroutine assemble_system (this)
    implicit none
    class(par_pb_bddc_linear_elasticity_fe_driver_t), intent(inout) :: this
    class(matrix_t)                  , pointer       :: matrix
    class(vector_t)                  , pointer       :: rhs
    call this%fe_affine_operator%numerical_setup()
    rhs                => this%fe_affine_operator%get_translation()
    matrix             => this%fe_affine_operator%get_matrix()

    !select type(matrix)
    !class is (sparse_matrix_t)  
    !   call matrix%print_matrix_market(6) 
    !class DEFAULT
    !   assert(.false.) 
    !end select

    !select type(rhs)
    !class is (serial_scalar_array_t)  
    !   call rhs%print(6) 
    !class DEFAULT
    !   assert(.false.) 
    !end select
  end subroutine assemble_system


  subroutine solve_system(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_fe_driver_t), intent(inout) :: this
    class(matrix_t)                         , pointer       :: matrix
    class(vector_t)                         , pointer       :: rhs
    class(vector_t)                         , pointer       :: dof_values

    matrix     => this%fe_affine_operator%get_matrix()
    rhs        => this%fe_affine_operator%get_translation()
    !write(*,*)'2-norm',rhs%nrm2()
    dof_values => this%solution%get_free_dof_values()
    call this%iterative_linear_solver%solve(-this%fe_affine_operator%get_translation(), &
         dof_values)
    !write(*,*)'solution',dof_values%nrm2()

    !select type (dof_values)
    !class is (serial_scalar_array_t)  
    !   call dof_values%print(6)
    !class DEFAULT
    !   assert(.false.) 
    !end select

    !select type (matrix)
    !class is (sparse_matrix_t)  
    !   call this%direct_solver%update_matrix(matrix, same_nonzero_pattern=.true.)
    !   call this%direct_solver%solve(rhs , dof_values )
    !class DEFAULT
    !   assert(.false.) 
    !end select
  end subroutine solve_system

  subroutine check_solution(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_fe_driver_t), intent(inout) :: this
    type(error_norms_scalar_t) :: scalar_error_norm 
    type(error_norms_vector_t) :: vector_error_norm 
    real(rp)    :: mean, l1, l2, lp, linfty, h1, h1_s, w1p_s, w1p, w1infty_s, w1infty
    integer(ip) :: field_id

    do field_id = 1, this%linear_elasticity_integration%get_number_fields()

       call vector_error_norm%create(this%fe_space,field_id)    
       mean      = vector_error_norm%compute(this%linear_elasticity_analytical_functions%get_solution_function_u(), this%solution, mean_norm)   
       l1        = vector_error_norm%compute(this%linear_elasticity_analytical_functions%get_solution_function_u(), this%solution, l1_norm)   
       l2        = vector_error_norm%compute(this%linear_elasticity_analytical_functions%get_solution_function_u(), this%solution, l2_norm)   
       lp        = vector_error_norm%compute(this%linear_elasticity_analytical_functions%get_solution_function_u(), this%solution, lp_norm)   
       linfty    = vector_error_norm%compute(this%linear_elasticity_analytical_functions%get_solution_function_u(), this%solution, linfty_norm)   
       h1_s      = vector_error_norm%compute(this%linear_elasticity_analytical_functions%get_solution_function_u(), this%solution, h1_seminorm) 
       h1        = vector_error_norm%compute(this%linear_elasticity_analytical_functions%get_solution_function_u(), this%solution, h1_norm) 
       w1p_s     = vector_error_norm%compute(this%linear_elasticity_analytical_functions%get_solution_function_u(), this%solution, w1p_seminorm)   
       w1p       = vector_error_norm%compute(this%linear_elasticity_analytical_functions%get_solution_function_u(), this%solution, w1p_norm)   
       w1infty_s = vector_error_norm%compute(this%linear_elasticity_analytical_functions%get_solution_function_u(), this%solution, w1infty_seminorm) 
       w1infty   = vector_error_norm%compute(this%linear_elasticity_analytical_functions%get_solution_function_u(), this%solution, w1infty_norm)  
       if ( this%par_environment%am_i_l1_root() ) then
          write(*,'(a12,a8)')      this%linear_elasticity_integration%get_field_name(field_id), ' field: '
          write(*,'(a20,e32.25)') 'mean_norm:', mean             ; check ( abs(mean) < 1.0e-04 ) 
          write(*,'(a20,e32.25)') 'l1_norm:', l1                 ; check ( l1 < 1.0e-04 )        
          write(*,'(a20,e32.25)') 'l2_norm:', l2                 ; check ( l2 < 1.0e-04 )        
          write(*,'(a20,e32.25)') 'lp_norm:', lp                 ; check ( lp < 1.0e-04 )        
          write(*,'(a20,e32.25)') 'linfnty_norm:', linfty        ; check ( linfty < 1.0e-04 )    
          write(*,'(a20,e32.25)') 'h1_seminorm:', h1_s           ; check ( h1_s < 1.0e-04 )      
          write(*,'(a20,e32.25)') 'h1_norm:', h1                 ; check ( h1 < 1.0e-04 )        
          write(*,'(a20,e32.25)') 'w1p_seminorm:', w1p_s         ; check ( w1p_s < 1.0e-04 )     
          write(*,'(a20,e32.25)') 'w1p_norm:', w1p               ; check ( w1p < 1.0e-04 )       
          write(*,'(a20,e32.25)') 'w1infty_seminorm:', w1infty_s ; check ( w1infty_s < 1.0e-04 ) 
          write(*,'(a20,e32.25)') 'w1infty_norm:', w1infty       ; check ( w1infty < 1.0e-04 )   
       end if
       call vector_error_norm%free()
    end do
  end subroutine check_solution

  !========================================================================================

  subroutine write_solution(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_fe_driver_t), intent(inout) :: this
    type(output_handler_t)                          :: oh
    real(rp),allocatable :: cell_vector(:)
    real(rp),allocatable :: mypart_vector(:)
    real(rp), allocatable :: set_id_cell_vector(:)

    if(this%test_params%get_write_solution()) then
       if (this%par_environment%am_i_l1_task()) then

          call memalloc(this%triangulation%get_num_local_cells(),mypart_vector,__FILE__,__LINE__)
          mypart_vector(:) = this%par_environment%get_l1_rank()

          call oh%create()
          call oh%attach_fe_space(this%fe_space)
          call oh%add_fe_function(this%solution, 1, 'solution')
          call oh%add_cell_vector(mypart_vector,'l1_rank')
          !call oh%add_cell_vector(set_id_cell_vector, 'set_id')
          call oh%open(this%test_params%get_dir_path(), this%test_params%get_prefix())
          call oh%write()
          call oh%close()
          call oh%free()

          if (allocated(cell_vector)) call memfree(cell_vector,__FILE__,__LINE__)
          call memfree(mypart_vector,__FILE__,__LINE__)

       end if
    endif
  end subroutine write_solution

  subroutine run_simulation(this) 
    implicit none
    class(par_pb_bddc_linear_elasticity_fe_driver_t), intent(inout) :: this

    call this%timer_triangulation%start()
    call this%setup_triangulation()
    call this%timer_triangulation%stop()

    call this%timer_fe_space%start()
    call this%setup_discrete_integration()
    call this%setup_reference_fes()
    call this%setup_coarse_fe_handlers()
    call this%setup_fe_space()
    call this%timer_fe_space%stop()

    call this%timer_assemply%start()
    call this%setup_system()
    call this%assemble_system()
    call this%timer_assemply%stop()

    call this%timer_solver_setup%start()
    call this%setup_solver()
    call this%timer_solver_setup%stop()

    call this%timer_solver_run%start()
    call this%solve_system()
    call this%timer_solver_run%stop()

    call this%check_solution()
    call this%write_solution()
    call this%free()
  end subroutine run_simulation

  subroutine free(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_fe_driver_t), intent(inout) :: this
    integer(ip) :: i, istat

    call this%solution%free() 
    call this%mlbddc%free()
    call this%iterative_linear_solver%free()
    call this%fe_affine_operator%free()
    call this%fe_space%free()
    if ( allocated(this%reference_fes) ) then
       do i=1, size(this%reference_fes)
          call this%reference_fes(i)%free()
       end do
       deallocate(this%reference_fes, stat=istat)
       check(istat==0)
    end if
    call this%triangulation%free()
    if (allocated(this%cell_set_ids)) call memfree(this%cell_set_ids,__FILE__,__LINE__)
    call this%linear_elasticity_integration%free()
  end subroutine free

  !========================================================================================
  subroutine free_environment(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_fe_driver_t), intent(inout) :: this
    call this%par_environment%free()
  end subroutine free_environment

  !========================================================================================
  subroutine free_command_line_parameters(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_fe_driver_t), intent(inout) :: this
    call this%test_params%free()
  end subroutine free_command_line_parameters


end module par_pb_bddc_linear_elasticity_driver_names
