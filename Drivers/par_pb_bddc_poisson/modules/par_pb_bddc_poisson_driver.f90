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
module par_pb_bddc_poisson_driver_names
  use fempar_names
  use par_pb_bddc_poisson_params_names
  use pb_bddc_poisson_cG_discrete_integration_names
  use pb_bddc_poisson_conditions_names
  use pb_bddc_poisson_analytical_functions_names
# include "debug.i90"

  implicit none
  private

  type par_pb_bddc_poisson_fe_driver_t 
     private 

     ! Place-holder for parameter-value set provided through command-line interface
     type(par_pb_bddc_poisson_params_t)      :: test_params
     type(ParameterList_t), pointer          :: parameter_list

     ! Cells and lower dimension objects container
     type(par_triangulation_t)             :: triangulation

     ! Discrete weak problem integration-related data type instances 
     type(par_fe_space_t)                      :: fe_space 
     type(p_reference_fe_t), allocatable       :: reference_fes(:) 
     type(H1_l1_coarse_fe_handler_t)             :: H1_coarse_fe_handler
     type(standard_l1_coarse_fe_handler_t)       :: standard_coarse_fe_handler
     type(p_l1_coarse_fe_handler_t), allocatable :: coarse_fe_handlers(:)
     type(poisson_CG_discrete_integration_t)   :: poisson_integration
     type(poisson_conditions_t)                :: poisson_conditions
     type(poisson_analytical_functions_t)      :: poisson_analytical_functions


     ! Place-holder for the coefficient matrix and RHS of the linear system
     type(fe_affine_operator_t)            :: fe_affine_operator

     ! MLBDDC preconditioner
     type(mlbddc_t)                            :: mlbddc

     ! Iterative linear solvers data type
     type(iterative_linear_solver_t)           :: iterative_linear_solver

     ! Poisson problem solution FE function
     type(fe_function_t)                   :: solution

     ! Environment required for fe_affine_operator + vtk_handler
     type(environment_t)                   :: environment
   contains
     procedure                  :: run_simulation
     procedure                  :: parse_command_line_parameters
     procedure                  :: setup_environment
     !procedure                  :: get_icontxt
     procedure        , private :: setup_triangulation
     procedure        , private :: setup_cell_set_ids
     procedure        , private :: setup_reference_fes
     procedure        , private :: setup_coarse_fe_handlers
     procedure        , private :: setup_fe_space
     procedure        , private :: setup_system
     procedure        , private :: setup_solver
     procedure        , private :: assemble_system
     procedure        , private :: solve_system
     procedure        , private :: check_solution
     procedure        , private :: write_solution
     procedure        , private :: write_matrices
     procedure        , private :: print_info
     procedure        , private :: free
     procedure                  :: free_environment
     procedure                  :: free_command_line_parameters
  end type par_pb_bddc_poisson_fe_driver_t

  ! Types
  public :: par_pb_bddc_poisson_fe_driver_t

contains

  subroutine parse_command_line_parameters(this)
    implicit none
    class(par_pb_bddc_poisson_fe_driver_t), intent(inout) :: this
    call this%test_params%create()
    this%parameter_list => this%test_params%get_values()
  end subroutine parse_command_line_parameters

  subroutine setup_triangulation(this)
    implicit none
    class(par_pb_bddc_poisson_fe_driver_t), intent(inout) :: this
    type(vef_iterator_t)  :: vef
    
    call this%triangulation%create(this%parameter_list, this%environment)

    if ( this%test_params%get_triangulation_type() == triangulation_generate_structured ) then
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
    end if

    if ( this%test_params%get_coarse_fe_handler_type() == pb_bddc ) then
      call this%setup_cell_set_ids() 
    end if  
    call this%triangulation%setup_coarse_triangulation()
    !write(*,*) 'CG: NUMBER OBJECTS', this%triangulation%get_number_objects()
    if ( this%test_params%get_coarse_fe_handler_type() == standard_bddc ) then
      call this%setup_cell_set_ids() 
    end if

  end subroutine setup_triangulation

  subroutine setup_cell_set_ids(this)
    implicit none
    class(par_pb_bddc_poisson_fe_driver_t), intent(inout) :: this
    class(cell_iterator_t), allocatable       :: cell
    integer(ip)                               :: istat, idummy
    integer(ip), allocatable                  :: cells_set(:)
    type(point_t), allocatable :: cell_coords(:)
    type(point_t) :: grav_center
    integer(ip)   :: inode, l1_rank  
    real(rp), allocatable:: px1(:), px2(:), py1(:), py2(:),  pz1(:), pz2(:)

    this%poisson_integration%diffusion_inclusion = this%test_params%get_jump()    
    this%H1_coarse_fe_handler%diffusion_inclusion = this%test_params%get_jump()

    if ( this%environment%am_i_l1_task() ) then
       l1_rank = this%environment%get_l1_rank()
       call memalloc( this%triangulation%get_num_local_cells(), cells_set, __FILE__, __LINE__ ) 
       call this%triangulation%create_cell_iterator(cell)
       allocate (cell_coords(cell%get_num_nodes()),stat=istat)
       do while( .not. cell%has_finished() )
          if ( cell%is_local() ) then
             call cell%get_coordinates(cell_coords)
             call grav_center%init(0.0_rp)
             do inode = 1, cell%get_num_nodes()
                grav_center = grav_center + cell_coords(inode)
             end do
             grav_center = (1.0_rp/cell%get_num_nodes())*grav_center
             cells_set( cell%get_lid() ) = cell_set_id( grav_center, &
                  this%triangulation%get_num_dimensions(), &
                  this%test_params%get_jump(), this%test_params%get_inclusion(), &
                  this%test_params%get_nchannel_per_direction(), &
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
    
    function cell_set_id( coord, num_dimensions, jump, inclusion, nchannel_per_direction, nparts_with_channels,nparts)
      implicit none
      type(point_t), intent(in)  :: coord
      integer(ip)  , intent(in)  :: num_dimensions
      integer(ip)  , intent(in)  :: jump
      integer(ip)  , intent(in)  :: inclusion
      integer(ip)  , intent(in)  :: nchannel_per_direction(3)
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
         if ( is_point_in_rectangle( origin, opposite, coord, num_dimensions ) ) cell_set_id = jump
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
            if ( is_point_in_rectangle( origin, opposite, coord, num_dimensions ) ) cell_set_id = jump
            ! 
            z_pos_0 = j/real(2*nchannel+1)
            z_pos_1 = (j+1)/real(2*nchannel+1)
            call origin%set(1,y_pos_0)  ; call origin%set(2, y_pos_0) ; call origin%set(3,z_pos_0);
            call opposite%set(1,y_pos_1); call opposite%set(2,1.0_rp); call opposite%set(3,z_pos_1);
            if ( is_point_in_rectangle( origin, opposite, coord, num_dimensions ) ) cell_set_id = jump
         end do
      else if ( inclusion == 3 ) then
         ! Hieu's test in PB-BDDC article (two channels)
         p1 = [2.0_rp/36.0_rp, 7.0_rp/36.0_rp, 14.0_rp/36.0_rp, 19.0_rp/36.0_rp, 26.0_rp/36.0_rp, 30.0_rp/36.0_rp] ! lower y value
         p2 = [4.0_rp/36.0_rp, 9.0_rp/36.0_rp, 16.0_rp/36.0_rp, 21.0_rp/36.0_rp, 28.0_rp/36.0_rp, 32.0_rp/36.0_rp] ! upper y value
         do i = 1, 6
            call origin%set(1,0.0_rp)  ; call origin%set(2, p1(i)) ; call origin%set(3,p1(7-i));
            call opposite%set(1,1.0_rp); call opposite%set(2,p2(i)); call opposite%set(3,p2(7-i));
            if ( is_point_in_rectangle( origin, opposite, coord, num_dimensions ) ) cell_set_id = jump + i - 1
         end do
         do i = 7, 12  
            call origin%set(2,0.0_rp)  ; call origin%set(1, p1(i-6)) ; call origin%set(3,p1(i-6));
            call opposite%set(2,1.0_rp); call opposite%set(1,p2(i-6)); call opposite%set(3,p2(i-6));
            if ( is_point_in_rectangle( origin, opposite, coord, num_dimensions ) ) cell_set_id = jump + i - 1
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
               if ( is_point_in_rectangle( origin, opposite, coord, num_dimensions ) ) cell_set_id = nchannel
            end do
         end do
         ! y edges
         do j = 1, 4
            do k = 1,4
               call origin%set(2,0.0_rp)  ; call origin%set(1, p1_b(j)) ; call origin%set(3,p1_b(k));
               call opposite%set(2,1.0_rp); call opposite%set(1,p2_b(j)); call opposite%set(3,p2_b(k));
               nchannel = nchannel + 1
               if ( is_point_in_rectangle( origin, opposite, coord, num_dimensions ) ) cell_set_id = nchannel
            end do
         end do
         ! z edges
         do j = 1, 4
            do k = 1,4
               call origin%set(3,0.0_rp)  ; call origin%set(2, p1_b(j)) ; call origin%set(1,p1_b(k));
               call opposite%set(3,1.0_rp); call opposite%set(2,p2_b(j)); call opposite%set(1,p2_b(k));
               nchannel = nchannel + 1
               if ( is_point_in_rectangle( origin, opposite, coord, num_dimensions ) ) cell_set_id = nchannel
            end do
         end do
      else if ( inclusion == 5 ) then
         ! Number of channels can be choosen from the command line using option -nc

         ! defining positions of the channels
         box_width = 1.0_rp/nchannel_per_direction(1)
         half_channel_width = box_width/5
         center = box_width/2
         eps = 1e-14_rp
         do j=1, nchannel_per_direction(1)
            p1_c(j)=center - half_channel_width 
            p2_c(j)=center + half_channel_width  
            center = center + box_width
         enddo

         nchannel = 1
         ! x edges
         do j = 1, nchannel_per_direction(1)
            do k = 1,nchannel_per_direction(1)
               call origin%set(1,0.0_rp)  ; call origin%set(2, p1_c(j)) ; call origin%set(3,p1_c(k));
               call opposite%set(1,1.0_rp); call opposite%set(2,p2_c(j)); call opposite%set(3,p2_c(k));
               nchannel = nchannel + 1
               if ( is_point_in_rectangle( origin, opposite, coord, num_dimensions ) ) cell_set_id = nchannel
            end do
         end do
         ! y edges
         do j = 1, nchannel_per_direction(1)
            do k = 1,nchannel_per_direction(1)
               call origin%set(2,0.0_rp)  ; call origin%set(1, p1_c(j)) ; call origin%set(3,p1_c(k));
               call opposite%set(2,1.0_rp); call opposite%set(1,p2_c(j)); call opposite%set(3,p2_c(k));
               nchannel = nchannel + 1
               if ( is_point_in_rectangle( origin, opposite, coord, num_dimensions ) ) cell_set_id = nchannel
            end do
         end do
         ! z edges
         do j = 1, nchannel_per_direction(1)
            do k = 1,nchannel_per_direction(1)
               call origin%set(3,0.0_rp)  ; call origin%set(2, p1_c(j)) ; call origin%set(1,p1_c(k));
               call opposite%set(3,1.0_rp); call opposite%set(2,p2_c(j)); call opposite%set(1,p2_c(k));
               nchannel = nchannel + 1
               if ( is_point_in_rectangle( origin, opposite, coord, num_dimensions ) ) cell_set_id = nchannel
            end do
         end do
      else if ( inclusion == 6 ) then
         ! Number of channels can be choosen from the command line using option -nc
         ! Channels are positioned so that some will touch but not cross the interface
         ! if the nchannel_per_direction is a multiple of the number partitions per direction 

         ! defining positions of the channels
         box_width = 1.0_rp/nchannel_per_direction(1)
         half_channel_width = box_width/5
         center = half_channel_width
         eps = 1e-14_rp
         do j=1, nchannel_per_direction(1)
            p1_c(j)=center - half_channel_width 
            p2_c(j)=center + half_channel_width  
            center = center + box_width
         enddo

         nchannel = 1
         ! x edges
         do j = 1, nchannel_per_direction(1)
            do k = 1,nchannel_per_direction(1)
               call origin%set(1,0.0_rp)  ; call origin%set(2, p1_c(j)) ; call origin%set(3,p1_c(k));
               call opposite%set(1,1.0_rp); call opposite%set(2,p2_c(j)); call opposite%set(3,p2_c(k));
               nchannel = nchannel + 1
               if ( is_point_in_rectangle( origin, opposite, coord, num_dimensions ) ) cell_set_id = nchannel
            end do
         end do
         ! y edges
         do j = 1, nchannel_per_direction(1)
            do k = 1,nchannel_per_direction(1)
               call origin%set(2,0.0_rp)  ; call origin%set(1, p1_c(j)) ; call origin%set(3,p1_c(k));
               call opposite%set(2,1.0_rp); call opposite%set(1,p2_c(j)); call opposite%set(3,p2_c(k));
               nchannel = nchannel + 1
               if ( is_point_in_rectangle( origin, opposite, coord, num_dimensions ) ) cell_set_id = nchannel
            end do
         end do
         ! z edges
         do j = 1, nchannel_per_direction(1)
            do k = 1,nchannel_per_direction(1)
               call origin%set(3,0.0_rp)  ; call origin%set(2, p1_c(j)) ; call origin%set(1,p1_c(k));
               call opposite%set(3,1.0_rp); call opposite%set(2,p2_c(j)); call opposite%set(1,p2_c(k));
               nchannel = nchannel + 1
               if ( is_point_in_rectangle( origin, opposite, coord, num_dimensions ) ) cell_set_id = nchannel
            end do
         end do
      else if ( inclusion == 7 ) then
         ! Number of channels can be choosen from the command line using option -nc
         ! Channels are positioned so that some will touch but not cross the interface
         ! if the nchannel_per_direction is a multiple of the number partitions per direction 

         ! defining positions of the channels
         if (.not.(allocated(px1))) then
            call memalloc(nchannel_per_direction(1), px1, __FILE__, __LINE__ )
            call memalloc(nchannel_per_direction(1), px2, __FILE__, __LINE__ )
            box_width = 1.0*nparts_with_channels(1)/(nchannel_per_direction(1)*nparts(1))
            half_channel_width = box_width/5
            center = half_channel_width
            do j=1, nchannel_per_direction(1)
               px1(j)=center - half_channel_width 
               px2(j)=center + half_channel_width  
               center = center + box_width
            enddo
            call memalloc(nchannel_per_direction(2), py1, __FILE__, __LINE__ )
            call memalloc(nchannel_per_direction(2), py2, __FILE__, __LINE__ )
            box_width = 1.0*nparts_with_channels(2)/(nchannel_per_direction(2)*nparts(2))
            half_channel_width = box_width/5
            center = half_channel_width
            do j=1, nchannel_per_direction(2)
               py1(j)=center - half_channel_width 
               py2(j)=center + half_channel_width  
               center = center + box_width
            enddo
            call memalloc(nchannel_per_direction(3), pz1, __FILE__, __LINE__ )
            call memalloc(nchannel_per_direction(3), pz2, __FILE__, __LINE__ )
            box_width = 1.0*nparts_with_channels(3)/(nchannel_per_direction(3)*nparts(3))
            half_channel_width = box_width/5
            center = half_channel_width
            do j=1, nchannel_per_direction(3)
               pz1(j)=center - half_channel_width 
               pz2(j)=center + half_channel_width  
               center = center + box_width
            enddo
         end if

         nchannel = 1
         ! x edges
         do j = 1, nchannel_per_direction(2)
            do k = 1,nchannel_per_direction(3)
               call origin%set(1,0.0_rp)  ; call origin%set(2, py1(j)) ; call origin%set(3,pz1(k));
               call opposite%set(1,1.0_rp); call opposite%set(2,py2(j)); call opposite%set(3,pz2(k));
               nchannel = nchannel + 1
               if ( is_point_in_rectangle( origin, opposite, coord, num_dimensions ) ) cell_set_id = nchannel
            end do
         end do
         ! y edges
         do j = 1, nchannel_per_direction(1)
            do k = 1,nchannel_per_direction(3)
               call origin%set(2,0.0_rp)  ; call origin%set(1, px1(j)) ; call origin%set(3,pz1(k));
               call opposite%set(2,1.0_rp); call opposite%set(1,px2(j)); call opposite%set(3,pz2(k));
               nchannel = nchannel + 1
               if ( is_point_in_rectangle( origin, opposite, coord, num_dimensions ) ) cell_set_id = nchannel
            end do
         end do
         ! z edges
         do j = 1, nchannel_per_direction(1)
            do k = 1,nchannel_per_direction(2)
               call origin%set(3,0.0_rp)  ; call origin%set(2, px1(j)) ; call origin%set(1,py1(k));
               call opposite%set(3,1.0_rp); call opposite%set(2,px2(j)); call opposite%set(1,py2(k));
               nchannel = nchannel + 1
               if ( is_point_in_rectangle( origin, opposite, coord, num_dimensions ) ) cell_set_id = nchannel
            end do
         end do
      else if ( inclusion == 8 ) then
         i = mod(mod(l1_rank,nparts(1)),2)
         j = mod(mod(l1_rank,nparts(1)*nparts(2))/nparts(2),2)
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

    function is_point_in_rectangle( origin, opposite, coord, num_dimensions )
      implicit none
      type(point_t), intent(in)  :: origin
      type(point_t), intent(in)  :: opposite
      type(point_t), intent(in)  :: coord
      integer(ip)  , intent(in)  :: num_dimensions
      logical :: is_point_in_rectangle
      integer(ip) :: i
      is_point_in_rectangle = .true.
      do i = 1, num_dimensions
         if ( coord%get(i) < origin%get(i) .or. coord%get(i) > opposite%get(i) ) then
            is_point_in_rectangle = .false.
            exit
         end if
      end do
    end function is_point_in_rectangle
    
  end subroutine setup_cell_set_ids


  subroutine setup_reference_fes(this)
    implicit none
    class(par_pb_bddc_poisson_fe_driver_t), intent(inout) :: this
    integer(ip) :: istat
    class(cell_iterator_t), allocatable       :: cell
    class(lagrangian_reference_fe_t), pointer :: reference_fe_geo

    allocate(this%reference_fes(1), stat=istat)
    check(istat==0)

    if ( this%environment%am_i_l1_task() ) then
       call this%triangulation%create_cell_iterator(cell)
       reference_fe_geo => cell%get_reference_fe_geo()
       this%reference_fes(1) =  make_reference_fe ( topology = reference_fe_geo%get_topology(), &
            fe_type = fe_type_lagrangian, &
            number_dimensions = this%triangulation%get_num_dimensions(), &
            order = this%test_params%get_reference_fe_order(), &
            field_type = field_type_scalar, &
            conformity = .true. )
       call this%triangulation%free_cell_iterator(cell)
    end if
  end subroutine setup_reference_fes

  subroutine setup_coarse_fe_handlers(this)
    implicit none
    class(par_pb_bddc_poisson_fe_driver_t), target, intent(inout) :: this
    integer(ip) :: istat
    
    allocate(this%coarse_fe_handlers(1), stat=istat)
    check(istat==0)

    if ( this%test_params%get_coarse_fe_handler_type() == pb_bddc ) then
       this%coarse_fe_handlers(1)%p => this%H1_coarse_fe_handler
    else if (this%test_params%get_coarse_fe_handler_type() == standard_bddc) then
       this%coarse_fe_handlers(1)%p => this%standard_coarse_fe_handler
    end if
  end subroutine setup_coarse_fe_handlers
  
  subroutine setup_fe_space(this)
    implicit none
    class(par_pb_bddc_poisson_fe_driver_t), intent(inout) :: this

    call this%fe_space%create( triangulation       = this%triangulation, &
                               reference_fes       = this%reference_fes, &
                               coarse_fe_handlers  = this%coarse_fe_handlers, &
                               conditions          = this%poisson_conditions )

    call this%fe_space%initialize_fe_integration()
    call this%fe_space%initialize_fe_face_integration()
    call this%poisson_analytical_functions%set_num_dimensions(this%triangulation%get_num_dimensions())
    call this%poisson_conditions%set_boundary_function(this%poisson_analytical_functions%get_boundary_function())
    call this%fe_space%interpolate_dirichlet_values(this%poisson_conditions)    
    !call this%fe_space%print()
  end subroutine setup_fe_space

  subroutine setup_system (this)
    implicit none
    class(par_pb_bddc_poisson_fe_driver_t), intent(inout) :: this

    call this%poisson_integration%set_analytical_functions(this%poisson_analytical_functions)

    ! if (test_single_scalar_valued_reference_fe) then
    call this%fe_affine_operator%create ( sparse_matrix_storage_format      = csr_format, &
         diagonal_blocks_symmetric_storage = [ .true. ], &
         diagonal_blocks_symmetric         = [ .true. ], &
         diagonal_blocks_sign              = [ SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE ], &
         fe_space                          = this%fe_space, &
         discrete_integration              = this%poisson_integration )
  end subroutine setup_system

  subroutine setup_solver (this)
    implicit none
    class(par_pb_bddc_poisson_fe_driver_t), intent(inout) :: this
    type(parameterlist_t) :: parameter_list
    type(parameterlist_t), pointer :: plist, dirichlet, neumann, coarse
    integer(ip) :: FPLError
    integer(ip) :: ilev

    call this%fe_space%setup_coarse_fe_space(this%parameter_list)
    
    plist => this%parameter_list 
    if ( this%environment%get_l1_size() == 1 ) then
       FPLError = plist%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
       !FPLError = plist%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_spd); assert(FPLError == 0)
       !FPLError = plist%set(key=pardiso_mkl_message_level, value=0); assert(FPLError == 0)
       !FPLError = plist%set(key=pardiso_mkl_iparm, value=iparm); assert(FPLError == 0)
    end if
    do ilev=1, this%environment%get_num_levels()-1
       ! Set current level Dirichlet solver parameters
       dirichlet => plist%NewSubList(key=mlbddc_dirichlet_solver_params)
       FPLError = dirichlet%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
       !FPLError = dirichlet%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_spd); assert(FPLError == 0)
       !FPLError = dirichlet%set(key=pardiso_mkl_message_level, value=0); assert(FPLError == 0)
       !FPLError = dirichlet%set(key=pardiso_mkl_iparm, value=iparm); assert(FPLError == 0)
       
       ! Set current level Neumann solver parameters
       neumann => plist%NewSubList(key=mlbddc_neumann_solver_params)
       FPLError = neumann%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
       !FPLError = neumann%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_sin); assert(FPLError == 0)
       !FPLError = neumann%set(key=pardiso_mkl_message_level, value=0); assert(FPLError == 0)
       !FPLError = neumann%set(key=pardiso_mkl_iparm, value=iparm); assert(FPLError == 0)
     
       coarse => plist%NewSubList(key=mlbddc_coarse_solver_params) 
       plist  => coarse 
    end do
    ! Set coarsest-grid solver parameters
    FPLError = coarse%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
    !FPLError = coarse%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_spd); assert(FPLError == 0)
    !FPLError = coarse%set(key=pardiso_mkl_message_level, value=0); assert(FPLError == 0)
    !FPLError = coarse%set(key=pardiso_mkl_iparm, value=iparm); assert(FPLError == 0)
    
    
    ! Set-up MLBDDC preconditioner
    call this%mlbddc%create(this%fe_affine_operator, this%parameter_list)
    call this%mlbddc%symbolic_setup()
    call this%mlbddc%numerical_setup()

    call parameter_list%init()
    FPLError = parameter_list%set(key = ils_rtol, value = 1.0e-6_rp)
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
    class(par_pb_bddc_poisson_fe_driver_t), intent(inout) :: this
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
    class(par_pb_bddc_poisson_fe_driver_t), intent(inout) :: this
    class(matrix_t)                         , pointer       :: matrix
    class(vector_t)                         , pointer       :: rhs
    class(vector_t)                         , pointer       :: dof_values

    matrix     => this%fe_affine_operator%get_matrix()
    rhs        => this%fe_affine_operator%get_translation()
    dof_values => this%solution%get_dof_values()
    call this%iterative_linear_solver%solve(this%fe_affine_operator%get_translation(), &
         dof_values)

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
    class(par_pb_bddc_poisson_fe_driver_t), intent(inout) :: this
    type(error_norms_scalar_t) :: error_norm 
    real(rp) :: mean, l1, l2, lp, linfty, h1, h1_s, w1p_s, w1p, w1infty_s, w1infty

    call error_norm%create(this%fe_space,1)    
    mean = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, mean_norm)   
    l1 = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, l1_norm)   
    l2 = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, l2_norm)   
    lp = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, lp_norm)   
    linfty = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, linfty_norm)   
    h1_s = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, h1_seminorm) 
    h1 = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, h1_norm) 
    w1p_s = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, w1p_seminorm)   
    w1p = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, w1p_norm)   
    w1infty_s = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, w1infty_seminorm) 
    w1infty = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, w1infty_norm)  
    if ( this%environment%get_l1_rank() == 0 ) then
       write(*,'(a20,e32.25)') 'mean_norm:', mean; !check ( abs(mean) < 1.0e-04 )
       write(*,'(a20,e32.25)') 'l1_norm:', l1; !check ( l1 < 1.0e-04 )
       write(*,'(a20,e32.25)') 'l2_norm:', l2; !check ( l2 < 1.0e-04 )
       write(*,'(a20,e32.25)') 'lp_norm:', lp; !check ( lp < 1.0e-04 )
       write(*,'(a20,e32.25)') 'linfnty_norm:', linfty; !check ( linfty < 1.0e-04 )
       write(*,'(a20,e32.25)') 'h1_seminorm:', h1_s; !check ( h1_s < 1.0e-04 )
       write(*,'(a20,e32.25)') 'h1_norm:', h1; !check ( h1 < 1.0e-04 )
       write(*,'(a20,e32.25)') 'w1p_seminorm:', w1p_s; !check ( w1p_s < 1.0e-04 )
       write(*,'(a20,e32.25)') 'w1p_norm:', w1p; !check ( w1p < 1.0e-04 )
       write(*,'(a20,e32.25)') 'w1infty_seminorm:', w1infty_s; !check ( w1infty_s < 1.0e-04 )
       write(*,'(a20,e32.25)') 'w1infty_norm:', w1infty; !check ( w1infty < 1.0e-04 )
    end if
    call error_norm%free()
  end subroutine check_solution

  subroutine write_solution(this)
    implicit none
    class(par_pb_bddc_poisson_fe_driver_t), intent(in) :: this
    type(output_handler_t)                             :: oh
    type(parameterlist_t)                              :: parameter_list
    class(base_static_triangulation_t), pointer        :: triangulation
    class(fe_iterator_t), allocatable :: fe
    real(rp), allocatable                              :: set_id_cell_vector(:)
    integer(ip)                                        :: i, istat
    if(this%test_params%get_write_solution()) then
       if ( this%environment%am_i_l1_task() ) then
          call build_set_id_cell_vector()
          call oh%create()
          call oh%attach_fe_space(this%fe_space)
          call oh%add_fe_function(this%solution, 1, 'solution')
          call oh%add_cell_vector(set_id_cell_vector, 'set_id')
          call parameter_list%init()
          !istat = parameter_list%set(key=vtk_format, value='ascii');
          call oh%open(this%test_params%get_dir_path_out(), this%test_params%get_prefix())!, parameter_list=parameter_list)
          call oh%write()
          call oh%close()
          call oh%free()
          call free_set_id_cell_vector()
          !call parameter_list%free()
       end if
    end if
  contains
    subroutine build_set_id_cell_vector()
      triangulation => this%fe_space%get_triangulation()
      call memalloc(triangulation%get_num_local_cells(), set_id_cell_vector, __FILE__, __LINE__)
      i = 1
      call this%fe_space%create_fe_iterator(fe)  
      do while ( .not. fe%has_finished())
         if (fe%is_local()) then
            set_id_cell_vector(i) = fe%get_set_id()
            i = i + 1
         end if
         call fe%next()
      enddo
      call this%fe_space%free_fe_iterator(fe)  
    end subroutine build_set_id_cell_vector

    subroutine free_set_id_cell_vector()
      call memfree(set_id_cell_vector, __FILE__, __LINE__)
    end subroutine free_set_id_cell_vector
  end subroutine write_solution

  subroutine write_matrices(this)
    implicit none
    class(par_pb_bddc_poisson_fe_driver_t), intent(in) :: this
    character(:), allocatable :: matrix_filename
    character(:), allocatable :: mapping_filename
    class(matrix_t), pointer :: matrix
    integer(ip) :: luout
    integer(igp) :: num_global_dofs
    integer(igp), allocatable :: dofs_gids(:)
    type(serial_scalar_array_t) :: mapping
    integer(ip) :: i
    
    if ( this%test_params%get_write_matrices() ) then
      if ( this%environment%am_i_l1_task() ) then
        matrix_filename = this%test_params%get_dir_path_out() // "/" // this%test_params%get_prefix() 
        call numbered_filename_compose(this%environment%get_l1_rank(),this%environment%get_l1_size(),matrix_filename)
        luout = io_open ( matrix_filename, 'write')
        matrix => this%fe_affine_operator%get_matrix()
        select type(matrix)
        class is (par_sparse_matrix_t)  
          call matrix%print_matrix_market(luout) 
        class DEFAULT
          assert(.false.) 
        end select
        call io_close(luout)
        
        call this%fe_space%compute_num_global_dofs_and_their_gids(num_global_dofs, dofs_gids)
        call mapping%create_and_allocate(size(dofs_gids))
        do i=1, size(dofs_gids)
          call mapping%insert(i,real(dofs_gids(i),rp))
        end do
        
        mapping_filename = this%test_params%get_dir_path_out() // "/" // this%test_params%get_prefix() // "_" //  "mapping"
        call numbered_filename_compose(this%environment%get_l1_rank(),this%environment%get_l1_size(),mapping_filename)
        luout = io_open ( mapping_filename, 'write')
        call mapping%print_matrix_market(luout)
        call io_close(luout)
        call mapping%free()
        call memfree(dofs_gids, __FILE__, __LINE__)
        
        if ( this%environment%get_l1_rank() == 0 ) then
          mapping_filename = this%test_params%get_dir_path_out() // "/" // this%test_params%get_prefix() // "_" //  "num_global_dofs"
          luout = io_open ( mapping_filename, 'write')
          write(luout,*) num_global_dofs
          call io_close(luout)
        end if  
        
      end if
   end if
  end subroutine write_matrices
  

  !========================================================================================
  subroutine print_info (this)
    implicit none
    class(par_pb_bddc_poisson_fe_driver_t), intent(in) :: this

    integer(ip) :: num_sub_domains
    real(rp) :: num_total_cells
    real(rp) :: num_dofs
    integer(ip) :: num_coarse_dofs

    class(environment_t), pointer :: environment
    class(coarse_fe_space_t), pointer :: coarse_fe_space

    environment => this%fe_space%get_environment()

    if (environment%am_i_l1_task()) then
      num_total_cells  = real(this%triangulation%get_num_local_cells(),kind=rp)
      num_dofs         = real(this%fe_space%get_field_number_dofs(1),kind=rp)
      call environment%l1_sum(num_total_cells )
      call environment%l1_sum(num_dofs        )
    end if

    if (environment%get_l1_rank() == 0) then
      num_sub_domains = environment%get_l1_size()
      write(*,'(a,i22)') 'num_sub_domains:          ', num_sub_domains
      write(*,'(a,i22)') 'num_total_cells:          ', nint(num_total_cells , kind=ip )
      write(*,'(a,i22)') 'num_dofs (sub-assembled): ', nint(num_dofs        , kind=ip )
    end if

    if (environment%am_i_lgt1_task()) then
      coarse_fe_space => this%fe_space%get_coarse_fe_space()
      num_coarse_dofs = coarse_fe_space%get_field_number_dofs(1)
      write(*,'(a,i22)') 'num_coarse_dofs:  ', num_coarse_dofs
    end if

  end subroutine print_info
  
  
  subroutine run_simulation(this) 
    implicit none
    class(par_pb_bddc_poisson_fe_driver_t), intent(inout) :: this
    type(timer_t) :: t_solve_system
    
    !call this%free()
    !call this%parse_command_line_parameters()
    call this%setup_triangulation()
    call this%setup_reference_fes()
    call this%setup_coarse_fe_handlers()
    call this%setup_fe_space()
    call this%setup_system()
    call this%assemble_system()
    
    call t_solve_system%create(this%environment%get_w_context() , "SOLVE SYSTEM", TIMER_MODE_MIN)
    call t_solve_system%start()
    call this%setup_solver()
    call this%solution%create(this%fe_space) 
    call this%solve_system()
    call t_solve_system%stop()
    call t_solve_system%report()
    
    !call this%check_solution()
    call this%write_solution()
    call this%write_matrices()
    call this%print_info()
    call this%free()
  end subroutine run_simulation

  subroutine free(this)
    implicit none
    class(par_pb_bddc_poisson_fe_driver_t), intent(inout) :: this
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
    if ( allocated(this%coarse_fe_handlers) ) then
       deallocate(this%coarse_fe_handlers, stat=istat)
       check(istat==0)
    end if
    
    call this%triangulation%free()
    !call this%test_params%free()
  end subroutine free


  !function get_icontxt(this)
  !  implicit none
  !  class(par_pb_bddc_poisson_fe_driver_t), intent(in) :: this
  !  integer(ip) :: get_icontxt
  !  class(execution_context_t), pointer :: w_context
  !  w_context => this%environment%get_w_context()
  !  select type(w_context)
  !  type is (mpi_context_t)
  !     get_icontxt = w_context%get_icontxt()
  !  end select
  !end function get_icontxt
  
  subroutine free_command_line_parameters(this)
    implicit none
    class(par_pb_bddc_poisson_fe_driver_t), intent(inout) :: this
    call this%test_params%free()
  end subroutine free_command_line_parameters

  subroutine setup_environment(this)
    implicit none
    class(par_pb_bddc_poisson_fe_driver_t), intent(inout) :: this
    integer(ip) :: istat
    if ( this%test_params%get_triangulation_type() == triangulation_generate_structured ) then
       istat = this%parameter_list%set(key = environment_type_key, value = structured) ; check(istat==0)
    else
       istat = this%parameter_list%set(key = environment_type_key, value = unstructured) ; check(istat==0)
    end if
    istat = this%parameter_list%set(key = execution_context_key, value = mpi_context) ; check(istat==0)
    call this%environment%create (this%parameter_list)
  end subroutine setup_environment

  subroutine free_environment(this)
    implicit none
    class(par_pb_bddc_poisson_fe_driver_t), intent(inout) :: this
    call this%environment%free()
  end subroutine free_environment

end module par_pb_bddc_poisson_driver_names
