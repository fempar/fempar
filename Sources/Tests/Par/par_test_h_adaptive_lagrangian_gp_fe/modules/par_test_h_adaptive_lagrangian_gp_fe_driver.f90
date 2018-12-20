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
module par_test_h_adaptive_lagrangian_gp_fe_driver_names
  use fempar_names
  use par_test_h_adaptive_lagrangian_gp_fe_params_names
  use lagrangian_gp_conditions_names
  use lagrangian_gp_analytical_functions_names
# include "debug.i90"

  implicit none
  private

  integer(ip), parameter :: PAR_TEST_LAGRANGIAN_GP_FULL = 1 ! Has to be == 1
  integer(ip), parameter :: PAR_TEST_LAGRANGIAN_GP_VOID = 2

  integer(ip), parameter :: TEMPERATURE_FIELD_ID  = 1

  type par_test_h_adaptive_lagrangian_gp_fe_driver_t 
     private 
     
     ! Place-holder for parameter-value set provided through command-line interface
     type(par_test_h_adaptive_lagrangian_gp_params_t) :: test_params
     type(ParameterList_t), pointer                   :: parameter_list
     
     ! Cells and lower dimension objects container
     type(p4est_par_triangulation_t)       :: triangulation
     integer(ip), allocatable              :: cell_set_ids(:)
     
     ! Discrete weak problem integration-related data type instances 
     type(par_fe_space_t)                              :: fe_space 
     type(p_reference_fe_t), allocatable               :: reference_fes(:) 
     type(h_adaptive_algebraic_l1_coarse_fe_handler_t) :: coarse_fe_handler
     type(p_l1_coarse_fe_handler_t), allocatable       :: coarse_fe_handlers(:)
     type(lagrangian_gp_conditions_t)                  :: lagrangian_gp_conditions
     type(lagrangian_gp_analytical_functions_t)        :: lagrangian_gp_analytical_functions

     ! Poisson problem solution FE function
     type(fe_function_t)                   :: solution
     
     ! Environment required for fe_affine_operator + vtk_handler
     type(environment_t)                    :: par_environment

     ! Timers
     type(timer_t) :: timer_triangulation
     type(timer_t) :: timer_fe_space
     type(timer_t) :: timer_assemply
     type(timer_t) :: timer_solver

   contains
     procedure                  :: parse_command_line_parameters
     procedure                  :: setup_timers
     procedure                  :: report_timers
     procedure                  :: free_timers
     procedure                  :: setup_environment
     procedure        , private :: setup_triangulation
     procedure        , private :: setup_reference_fes
     procedure        , private :: setup_coarse_fe_handlers
     procedure        , private :: setup_fe_space
     procedure        , private :: interpolate_analytical_functions
     procedure        , private :: check_solution
     procedure        , private :: write_solution
     procedure                  :: run_simulation
     procedure        , private :: print_info 
     procedure        , private :: free
     procedure                  :: free_command_line_parameters
     procedure                  :: free_environment
     procedure, nopass, private :: popcorn_fun => par_test_h_adaptive_lagrangian_gp_driver_popcorn_fun
     procedure                  :: set_cells_for_uniform_refinement
     procedure                  :: set_cells_for_refinement
     procedure                  :: set_cells_for_coarsening
     procedure                  :: refine_and_coarsen
     procedure                  :: set_cells_weights
  end type par_test_h_adaptive_lagrangian_gp_fe_driver_t 

  ! Types
  public :: par_test_h_adaptive_lagrangian_gp_fe_driver_t

contains

  subroutine parse_command_line_parameters(this)
    implicit none
    class(par_test_h_adaptive_lagrangian_gp_fe_driver_t), intent(inout) :: this
    call this%test_params%create()
    this%parameter_list => this%test_params%get_values()
  end subroutine parse_command_line_parameters

!========================================================================================
subroutine setup_timers(this)
    implicit none
    class(par_test_h_adaptive_lagrangian_gp_fe_driver_t), intent(inout) :: this
    class(execution_context_t), pointer :: w_context
    w_context => this%par_environment%get_w_context()
    call this%timer_triangulation%create(w_context,"SETUP TRIANGULATION")
    call this%timer_fe_space%create(     w_context,"SETUP FE SPACE")
    call this%timer_assemply%create(     w_context,"FE INTEGRATION AND ASSEMBLY")
    call this%timer_solver%create(   w_context,"PRECONDITIONER SETUP + SOLVER RUN")
end subroutine setup_timers

!========================================================================================
subroutine report_timers(this)
    implicit none
    class(par_test_h_adaptive_lagrangian_gp_fe_driver_t), intent(inout) :: this
    call this%timer_triangulation%report(.true.)
    call this%timer_fe_space%report(.false.)
    call this%timer_assemply%report(.false.)
    call this%timer_solver%report(.false.)
    if (this%par_environment%get_l1_rank() == 0) then
      write(*,*)
    end if
end subroutine report_timers

!========================================================================================
subroutine free_timers(this)
    implicit none
    class(par_test_h_adaptive_lagrangian_gp_fe_driver_t), intent(inout) :: this
    call this%timer_triangulation%free()
    call this%timer_fe_space%free()
    call this%timer_assemply%free()
    call this%timer_solver%free()
end subroutine free_timers

!========================================================================================
  subroutine setup_environment(this, world_context)
    implicit none
    class(par_test_h_adaptive_lagrangian_gp_fe_driver_t), intent(inout) :: this
    class(execution_context_t)                    , intent(in)    :: world_context
    integer(ip) :: istat
    istat = this%parameter_list%set(key = environment_type_key, value = p4est) ; check(istat==0)
    call this%par_environment%create (world_context, this%parameter_list)
  end subroutine setup_environment
   
  subroutine setup_triangulation(this)
    implicit none
    class(par_test_h_adaptive_lagrangian_gp_fe_driver_t), intent(inout) :: this

    class(cell_iterator_t), allocatable :: cell, coarser_cell
    type(point_t), allocatable :: cell_coords(:)
    integer(ip) :: istat
    integer(ip) :: set_id
    real(rp) :: x, y
    integer(ip) :: num_void_neigs
    integer(ip) :: i

    integer(ip)           :: ivef
    class(vef_iterator_t), allocatable  :: vef, vef_of_vef
    type(list_t), pointer :: vefs_of_vef
    type(list_t), pointer :: vertices_of_line
    type(list_iterator_t) :: vefs_of_vef_iterator
    type(list_iterator_t) :: vertices_of_line_iterator
    class(reference_fe_t), pointer :: reference_fe_geo
    integer(ip) :: ivef_pos_in_cell, vef_of_vef_pos_in_cell
    integer(ip) :: vertex_pos_in_cell, icell_arround
    integer(ip) :: inode, num
    class(environment_t), pointer :: environment
    real(rp)                      :: domain(6)
    character(len=:), allocatable :: subparts_coupling_criteria
    integer(ip) :: num_refs
    integer(ip) :: vef_set_id

    ! Create a p4est mesh with a custom domain 
    domain = this%test_params%get_domain_limits() 
    subparts_coupling_criteria = this%test_params%get_subparts_coupling_criteria()
    istat = this%parameter_list%set(key = struct_hex_triang_domain_limits_key , value = domain); check(istat==0)
    istat = this%parameter_list%set(key = subparts_coupling_criteria_key, value = subparts_coupling_criteria); check(istat==0)
    call this%triangulation%create(this%par_environment,this%parameter_list)
    
    ! Generate initial uniform mesh and set the cell ids to use void fes
    if (this%test_params%get_use_void_fes()) then
    
      do i = 1, 3
        call this%set_cells_for_uniform_refinement()
        call this%triangulation%refine_and_coarsen()
        call this%set_cells_weights()
        call this%triangulation%redistribute()
        call this%triangulation%clear_refinement_and_coarsening_flags()
      end do
    
      if ( this%par_environment%am_i_l1_task() ) then
        call memalloc(this%triangulation%get_num_local_cells(),this%cell_set_ids)
        call this%triangulation%create_cell_iterator(cell)
        allocate(cell_coords(1:cell%get_num_nodes()),stat=istat); check(istat == 0)
        do while( .not. cell%has_finished() )
          if (cell%is_local()) then
            set_id = PAR_TEST_LAGRANGIAN_GP_VOID
            call cell%get_nodes_coordinates(cell_coords)
            select case (trim(this%test_params%get_use_void_fes_case()))
            case ('half')
              y = cell_coords(1)%get(2)
              if (y>=0.5) set_id = PAR_TEST_LAGRANGIAN_GP_FULL
            case ('quarter')
              x = cell_coords(1)%get(1)
              y = cell_coords(1)%get(2)
              if (x>=0.5 .and. y>=0.5) set_id = PAR_TEST_LAGRANGIAN_GP_FULL
            case ('popcorn')
              do inode = 1,cell%get_num_nodes()
                if ( this%popcorn_fun(cell_coords(inode),&
                  this%triangulation%get_num_dims()) < 0.0 ) then
                  set_id = PAR_TEST_LAGRANGIAN_GP_FULL
                  exit
                end if
              end do
            case default
              check(.false.)
            end select
            this%cell_set_ids(cell%get_gid()) = set_id
          end if
          call cell%next()
        end do
        deallocate(cell_coords, stat = istat); check(istat == 0)
        call this%triangulation%fill_cells_set(this%cell_set_ids)
      end if
        call this%triangulation%free_cell_iterator(cell)
    end if

    ! Setup Dirichlet BCs on the boundary of the initial conforming uniform mesh
    environment => this%triangulation%get_environment()
    if (environment%am_i_l1_task()) then
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

    ! For the popcorn example, set all the vefs on the interface between 
    ! full/void if there are void fes on the Dirichlet boundary. 
    ! NOTE: This algorithm is only valid for conforming meshes.
    
    if (environment%am_i_l1_task()) then
    
      if ( this%test_params%get_use_void_fes() .and. &
           trim(this%test_params%get_use_void_fes_case()) == 'popcorn' ) then
        call this%triangulation%create_vef_iterator(vef)
        call this%triangulation%create_vef_iterator(vef_of_vef)
        call this%triangulation%create_cell_iterator(cell)
        call this%triangulation%create_cell_iterator(coarser_cell)
        do while ( .not. vef%has_finished() )

           ! If it is an INTERIOR face
           if( vef%get_dim() == this%triangulation%get_num_dims()-1 & 
               .and. .not. vef%is_at_boundary() ) then
   
             ! Compute number of void neighbors
             num_void_neigs = 0
             do icell_arround = 1,vef%get_num_cells_around()
               call vef%get_cell_around(icell_arround,cell)
               if (cell%get_set_id() == PAR_TEST_LAGRANGIAN_GP_VOID) num_void_neigs = num_void_neigs + 1
             end do
             
             do icell_arround = 1,vef%get_num_improper_cells_around()
               call vef%get_improper_cell_around(icell_arround,cell)
               if (cell%get_set_id() == PAR_TEST_LAGRANGIAN_GP_VOID) num_void_neigs = num_void_neigs + 1
             end do

             if(num_void_neigs==1) then ! If vef (face) is between a full and a void cell
                 
                 ! Set this vef and vefs of this vef at the Dirichlet boundary
                 if ( vef%is_proper() .or. all_coarser_cells_around_vef_are_void() ) call vef%set_set_id(1)
                 
                 ! Do a loop on all edges in 3D (vertex in 2D) of the face
                 call vef%get_cell_around(1,cell) ! There is always one cell around
                 reference_fe_geo => cell%get_reference_fe()
                 ivef_pos_in_cell = cell%get_vef_lid_from_gid(vef%get_gid())
                 vefs_of_vef => reference_fe_geo%get_facets_n_face()
                 vefs_of_vef_iterator = vefs_of_vef%create_iterator(ivef_pos_in_cell)
                 do while( .not. vefs_of_vef_iterator%is_upper_bound() )

                    ! Set edge (resp. vertex) as Dirichlet
                    vef_of_vef_pos_in_cell = vefs_of_vef_iterator%get_current()
                    call cell%get_vef(vef_of_vef_pos_in_cell, vef_of_vef)
                    if ( vef%is_proper() .or. all_coarser_cells_around_vef_are_void() ) call vef_of_vef%set_set_id(1)

                    ! If 3D, traverse vertices of current line
                    if ( this%triangulation%get_num_dims() == 3 ) then
                      vertices_of_line          => reference_fe_geo%get_vertices_n_face()
                      vertices_of_line_iterator = vertices_of_line%create_iterator(vef_of_vef_pos_in_cell)
                      do while( .not. vertices_of_line_iterator%is_upper_bound() )

                        ! Set vertex as Dirichlet
                        vertex_pos_in_cell = vertices_of_line_iterator%get_current()
                        call cell%get_vef(vertex_pos_in_cell, vef_of_vef)
                        if ( vef%is_proper() .or. all_coarser_cells_around_vef_are_void() ) call vef_of_vef%set_set_id(1)

                        call vertices_of_line_iterator%next()
                      end do ! Loop in vertices in 3D only
                    end if

                    call vefs_of_vef_iterator%next()
                 end do ! Loop in edges (resp. vertices)

             end if ! If face on void/full boundary
           end if ! If vef is an interior face

           call vef%next()
        end do ! Loop in vefs
        call this%triangulation%free_cell_iterator(cell)
        call this%triangulation%free_cell_iterator(coarser_cell)
        call this%triangulation%free_vef_iterator(vef)
        call this%triangulation%free_vef_iterator(vef_of_vef)
      end if
    end if
    
    num_refs = this%test_params%get_num_refinements() 
    if ( this%test_params%get_use_void_fes() ) then
      num_refs = num_refs-3
    end if  
    
    do i = 1, num_refs
      call this%set_cells_for_refinement()
      call this%triangulation%refine_and_coarsen()
      call this%set_cells_weights()
      call this%triangulation%redistribute()
      call this%triangulation%clear_refinement_and_coarsening_flags()
    end do
#ifdef ENABLE_MKL    
    call this%triangulation%setup_coarse_triangulation()
#endif
   
   contains
     
     function all_coarser_cells_around_vef_are_void()
       implicit none
       logical :: all_coarser_cells_around_vef_are_void
       all_coarser_cells_around_vef_are_void = .true.
       do icell_arround = 1,vef%get_num_improper_cells_around()
         call vef%get_improper_cell_around(icell_arround,coarser_cell)
         if (coarser_cell%get_set_id() == PAR_TEST_LAGRANGIAN_GP_FULL) all_coarser_cells_around_vef_are_void = .false.
         exit
       end do
     end function all_coarser_cells_around_vef_are_void
     
  end subroutine setup_triangulation
  
  subroutine setup_reference_fes(this)
    implicit none
    class(par_test_h_adaptive_lagrangian_gp_fe_driver_t), intent(inout) :: this
    integer(ip) :: istat
    class(cell_iterator_t), allocatable       :: cell
    class(reference_fe_t), pointer :: reference_fe_geo
    
    if (this%test_params%get_use_void_fes()) then
      allocate(this%reference_fes(2), stat=istat)
    else
      allocate(this%reference_fes(1), stat=istat)
    end if
    check(istat==0)
    
    if ( this%par_environment%am_i_l1_task() ) then
      call this%triangulation%create_cell_iterator(cell)
      reference_fe_geo => cell%get_reference_fe()
      this%reference_fes(PAR_TEST_LAGRANGIAN_GP_FULL) =  make_reference_fe ( topology = reference_fe_geo%get_topology(), &
                                                   fe_type = fe_type_lagrangian_gp, &
                                                   num_dims = this%triangulation%get_num_dims(), &
                                                   order = this%test_params%get_reference_fe_order(), &
                                                   field_type = field_type_scalar, &
                                                   conformity = .true., &
                                                   continuity = .false. )
      if (this%test_params%get_use_void_fes()) then
        this%reference_fes(PAR_TEST_LAGRANGIAN_GP_VOID) =  make_reference_fe ( topology = reference_fe_geo%get_topology(), &
                                                   fe_type = fe_type_void, &
                                                   num_dims = this%triangulation%get_num_dims(), &
                                                   order = -1, &
                                                   field_type = field_type_scalar, &
                                                   conformity = .true. )
      end if
      call this%triangulation%free_cell_iterator(cell)
    end if
  end subroutine setup_reference_fes
  
  subroutine setup_coarse_fe_handlers(this)
    implicit none
    class(par_test_h_adaptive_lagrangian_gp_fe_driver_t), target, intent(inout) :: this
    integer(ip) :: istat
    allocate(this%coarse_fe_handlers(1), stat=istat)
    check(istat==0)
    this%coarse_fe_handlers(1)%p => this%coarse_fe_handler
  end subroutine setup_coarse_fe_handlers

  subroutine setup_fe_space(this)
    implicit none
    class(par_test_h_adaptive_lagrangian_gp_fe_driver_t), intent(inout) :: this

    integer(ip) :: set_ids_to_reference_fes(1,2)

    call this%lagrangian_gp_analytical_functions%set_num_dims(this%triangulation%get_num_dims())
    call this%lagrangian_gp_analytical_functions%set_solution_polynomial_order(this%test_params%get_reference_fe_order())
    call this%lagrangian_gp_conditions%set_boundary_function(this%lagrangian_gp_analytical_functions%get_boundary_function())
 
    if (this%test_params%get_use_void_fes()) then
      set_ids_to_reference_fes(1,PAR_TEST_LAGRANGIAN_GP_FULL) = PAR_TEST_LAGRANGIAN_GP_FULL
      set_ids_to_reference_fes(1,PAR_TEST_LAGRANGIAN_GP_VOID) = PAR_TEST_LAGRANGIAN_GP_VOID
      call this%fe_space%create( triangulation            = this%triangulation,       &
                                 reference_fes            = this%reference_fes,       &
                                 set_ids_to_reference_fes = set_ids_to_reference_fes, &
                                 coarse_fe_handlers       = this%coarse_fe_handlers,  &
                                 conditions               = this%lagrangian_gp_conditions )
    else
      call this%fe_space%create( triangulation       = this%triangulation,      &
                                 reference_fes       = this%reference_fes,      &
                                 coarse_fe_handlers  = this%coarse_fe_handlers, &
                                 conditions          = this%lagrangian_gp_conditions )
    end if
    
    call this%fe_space%set_up_cell_integration()
    call this%fe_space%set_up_facet_integration() 
#ifdef ENABLE_MKL    
    call this%fe_space%setup_coarse_fe_space(this%parameter_list)
#endif    
    !call this%fe_space%print()
  end subroutine setup_fe_space

  subroutine check_solution(this)
    implicit none
    class(par_test_h_adaptive_lagrangian_gp_fe_driver_t), intent(inout) :: this
    type(error_norms_scalar_t) :: error_norm 
    real(rp) :: mean, l1, l2, lp, linfty, h1, h1_s, w1p_s, w1p, w1infty_s, w1infty
    real(rp) :: requested_precision
    
    call error_norm%create(this%fe_space,1)    
    mean = error_norm%compute(this%lagrangian_gp_analytical_functions%get_solution_function(), this%solution, mean_norm)   
    l1 = error_norm%compute(this%lagrangian_gp_analytical_functions%get_solution_function(), this%solution, l1_norm)   
    l2 = error_norm%compute(this%lagrangian_gp_analytical_functions%get_solution_function(), this%solution, l2_norm)   
    lp = error_norm%compute(this%lagrangian_gp_analytical_functions%get_solution_function(), this%solution, lp_norm)   
    linfty = error_norm%compute(this%lagrangian_gp_analytical_functions%get_solution_function(), this%solution, linfty_norm)   
    h1_s = error_norm%compute(this%lagrangian_gp_analytical_functions%get_solution_function(), this%solution, h1_seminorm) 
    h1 = error_norm%compute(this%lagrangian_gp_analytical_functions%get_solution_function(), this%solution, h1_norm) 
    w1p_s = error_norm%compute(this%lagrangian_gp_analytical_functions%get_solution_function(), this%solution, w1p_seminorm)   
    w1p = error_norm%compute(this%lagrangian_gp_analytical_functions%get_solution_function(), this%solution, w1p_norm)   
    w1infty_s = error_norm%compute(this%lagrangian_gp_analytical_functions%get_solution_function(), this%solution, w1infty_seminorm) 
    w1infty = error_norm%compute(this%lagrangian_gp_analytical_functions%get_solution_function(), this%solution, w1infty_norm)  
    if ( this%par_environment%am_i_l1_root() ) then
      requested_precision=1.0e-8_rp
      write(*,'(a20,e32.25)') 'mean_norm:', mean; check ( abs(mean) < requested_precision )
      write(*,'(a20,e32.25)') 'l1_norm:', l1; check ( l1 < requested_precision )
      write(*,'(a20,e32.25)') 'l2_norm:', l2; check ( l2 < requested_precision )
      write(*,'(a20,e32.25)') 'lp_norm:', lp; check ( lp < requested_precision )
      write(*,'(a20,e32.25)') 'linfnty_norm:', linfty; check ( linfty < requested_precision )
      write(*,'(a20,e32.25)') 'h1_seminorm:', h1_s; check ( h1_s < requested_precision )
      write(*,'(a20,e32.25)') 'h1_norm:', h1; check ( h1 < requested_precision )
      write(*,'(a20,e32.25)') 'w1p_seminorm:', w1p_s; check ( w1p_s < requested_precision )
      write(*,'(a20,e32.25)') 'w1p_norm:', w1p; check ( w1p < requested_precision )
      write(*,'(a20,e32.25)') 'w1infty_seminorm:', w1infty_s; check ( w1infty_s < requested_precision )
      write(*,'(a20,e32.25)') 'w1infty_norm:', w1infty; check ( w1infty < requested_precision )
    end if  
    call error_norm%free()
  end subroutine check_solution
  
  subroutine write_solution(this)
    implicit none
    class(par_test_h_adaptive_lagrangian_gp_fe_driver_t), intent(in) :: this
    type(output_handler_t)                          :: oh
    real(rp),allocatable :: cell_vector(:)
    real(rp),allocatable :: mypart_vector(:)
    class(cell_iterator_t), allocatable :: cell 

    if(this%test_params%get_write_solution()) then
      if (this%par_environment%am_i_l1_task()) then

        if (this%test_params%get_use_void_fes()) then
          call memalloc(this%triangulation%get_num_local_cells(),cell_vector,__FILE__,__LINE__)
          call this%triangulation%create_cell_iterator(cell)
          do while ( .not. cell%has_finished() )
            if ( cell%is_local() ) then       
              cell_vector(cell%get_gid()) = cell%get_set_id()
             end if 
             call cell%next()
          end do 
          call this%triangulation%free_cell_iterator(cell)
        end if

        call memalloc(this%triangulation%get_num_local_cells(),mypart_vector,__FILE__,__LINE__)
        mypart_vector(:) = this%par_environment%get_l1_rank()

        call oh%create()
        call oh%attach_fe_space(this%fe_space)
        call oh%add_fe_function(this%solution, 1, 'solution')
        if (this%test_params%get_use_void_fes()) then
          call oh%add_cell_vector(cell_vector,'cell_set_ids')
        end if
        call oh%add_cell_vector(mypart_vector,'l1_rank')
        call oh%open(this%test_params%get_dir_path(), this%test_params%get_prefix())
        call oh%write()
        call oh%close()
        call oh%free()

        if (allocated(cell_vector)) call memfree(cell_vector,__FILE__,__LINE__)
        call memfree(mypart_vector,__FILE__,__LINE__)

      end if
    endif
  end subroutine write_solution
  
  subroutine interpolate_analytical_functions(this)
    implicit none
    class(par_test_h_adaptive_lagrangian_gp_fe_driver_t), intent(inout) :: this    

    ! Create solution function
    call this%solution%create(this%fe_space)

    ! Project analytical function to solution function using the gp-lagrangian fe
    call this%fe_space%interpolate( TEMPERATURE_FIELD_ID, this%lagrangian_gp_analytical_functions%get_solution_function(), this%solution )
    call this%fe_space%interpolate_dirichlet_values( this%solution )  

  end subroutine interpolate_analytical_functions
  
  subroutine run_simulation(this) 
    implicit none
    class(par_test_h_adaptive_lagrangian_gp_fe_driver_t), intent(inout) :: this

    call this%timer_triangulation%start()
    call this%setup_triangulation()
    call this%timer_triangulation%stop()

    
    call this%timer_fe_space%start()
    call this%setup_reference_fes()
    call this%setup_coarse_fe_handlers()
    call this%setup_fe_space()
    call this%timer_fe_space%stop()
    
    call this%interpolate_analytical_functions()
    
    call this%check_solution()
    call this%print_info() 
    
    call this%refine_and_coarsen()
    
    call this%set_cells_weights()
    call this%triangulation%redistribute()
    call this%fe_space%redistribute(this%solution)
    call this%fe_space%set_up_cell_integration()
    
    call this%check_solution()
    call this%write_solution()
    call this%free()
  end subroutine run_simulation
  
  subroutine refine_and_coarsen(this)
    implicit none
    class(par_test_h_adaptive_lagrangian_gp_fe_driver_t), intent(inout) :: this
    integer(ip) :: i
    
    do i=1,3
       
       if ( i == 2 ) then
         call this%set_cells_for_coarsening()
       else
         call this%set_cells_for_refinement()
       end if

       call this%triangulation%refine_and_coarsen()
       call this%fe_space%refine_and_coarsen( this%solution ) 
       call this%fe_space%set_up_cell_integration()

       call this%check_solution()

    end do  
    
  end subroutine refine_and_coarsen  
  
  subroutine free(this)
    implicit none
    class(par_test_h_adaptive_lagrangian_gp_fe_driver_t), intent(inout) :: this
    integer(ip) :: i, istat
    
    call this%solution%free()
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
    if (allocated(this%cell_set_ids)) call memfree(this%cell_set_ids,__FILE__,__LINE__)
  end subroutine free  

  !========================================================================================
  subroutine free_environment(this)
    implicit none
    class(par_test_h_adaptive_lagrangian_gp_fe_driver_t), intent(inout) :: this
    call this%par_environment%free()
  end subroutine free_environment

  !========================================================================================
  subroutine free_command_line_parameters(this)
    implicit none
    class(par_test_h_adaptive_lagrangian_gp_fe_driver_t), intent(inout) :: this
    call this%test_params%free()
  end subroutine free_command_line_parameters

  function par_test_h_adaptive_lagrangian_gp_driver_popcorn_fun(point,num_dim) result (val)
    implicit none
    type(point_t), intent(in) :: point
    integer(ip),   intent(in) :: num_dim
    real(rp) :: val
    type(point_t) :: p
    real(rp) :: x, y, z
    real(rp) :: xk, yk, zk
    real(rp) :: r0, sg, A
    integer(ip) :: k
    p = point
    if (num_dim < 3) call p%set(3,0.62)
    x = ( 2.0*p%get(1) - 1.0 )
    y = ( 2.0*p%get(2) - 1.0 )
    z = ( 2.0*p%get(3) - 1.0 )
    r0 = 0.6
    sg = 0.2
    A  = 2.0
    val = sqrt(x**2 + y**2 + z**2) - r0
    do k = 0,11
        if (0 <= k .and. k <= 4) then
            xk = (r0/sqrt(5.0))*2.0*cos(2.0*k*pi/5.0)
            yk = (r0/sqrt(5.0))*2.0*sin(2.0*k*pi/5.0)
            zk = (r0/sqrt(5.0))
        else if (5 <= k .and. k <= 9) then
            xk = (r0/sqrt(5.0))*2.0*cos((2.0*(k-5)-1.0)*pi/5.0)
            yk = (r0/sqrt(5.0))*2.0*sin((2.0*(k-5)-1.0)*pi/5.0)
            zk =-(r0/sqrt(5.0))
        else if (k == 10) then
            xk = 0
            yk = 0
            zk = r0
        else
            xk = 0
            yk = 0
            zk = -r0
        end if
        val = val - A*exp( -( (x - xk)**2  + (y - yk)**2 + (z - zk)**2 )/(sg**2) )
    end do
  end function par_test_h_adaptive_lagrangian_gp_driver_popcorn_fun

  subroutine set_cells_for_refinement(this)
    implicit none
    class(par_test_h_adaptive_lagrangian_gp_fe_driver_t), intent(inout) :: this
    class(cell_iterator_t), allocatable :: cell
    class(environment_t), pointer :: environment
    type(point_t), allocatable    :: cell_coordinates(:)
    real(rp) :: cx, cy, cz
    integer(ip) :: inode, istat
    character(len=:), allocatable :: refinement_pattern_case
    ! Centered refined pattern 
    real(rp)               :: inner_region_size(0:SPACE_DIM-1)  
    real(rp)               :: domain(6)
    real(rp)               :: domain_length(0:SPACE_DIM-1) 
    logical, allocatable   :: is_node_coord_within_inner_region(:)  
    integer(ip)            :: idime 

    environment => this%triangulation%get_environment()

    if ( environment%am_i_l1_task() ) then
       call this%triangulation%create_cell_iterator(cell)
       allocate(cell_coordinates( cell%get_num_nodes() ) , stat=istat); check(istat==0)

       if ( this%test_params%get_refinement_pattern_case() == inner_region ) then 
          call memalloc(this%triangulation%get_num_dims(), is_node_coord_within_inner_region, __FILE__, __LINE__ ) 
          inner_region_size = this%test_params%get_inner_region_size() 
          domain = this%test_params%get_domain_limits()
          domain_length(0) = domain(2)-domain(1) 
          domain_length(1) = domain(4)-domain(3) 
          domain_length(2) = domain(6)-domain(5) 
       end if

       do while ( .not. cell%has_finished() )
          if ( cell%is_local() ) then       
             if ( cell%get_level() == 0 ) then
                call cell%set_for_refinement() 
                call cell%next(); cycle 
             end if
             call cell%get_nodes_coordinates(cell_coordinates)
             cx = 0.0_rp
             cy = 0.0_rp 
             cz = 0.0_rp 
             select case ( this%test_params%get_refinement_pattern_case() ) 
             case ( uniform )
                call cell%set_for_refinement()
             case ( even_cells ) 
                if ( (mod(cell%get_gid(),2)==0) )then
                   call cell%set_for_refinement()
                end if
             case ( inner_region ) 
                node_loop: do inode=1, cell%get_num_nodes()      
                   is_node_coord_within_inner_region=.false. 
                   do idime=0, this%triangulation%get_num_dims()-1
                      if ( (domain_length(idime)-inner_region_size(idime))/2.0_rp <= cell_coordinates(inode)%get(idime+1) .and. &
                           (domain_length(idime)+inner_region_size(idime))/2.0_rp >= cell_coordinates(inode)%get(idime+1) ) then 
                         is_node_coord_within_inner_region(idime+1)=.true. 
                      else 
                         cycle node_loop   
                      end if
                   end do
                   if ( all(is_node_coord_within_inner_region) ) then 
                   call cell%set_for_refinement(); exit 
                   end if 
                end do node_loop
             case DEFAULT 
                massert(.false., 'Refinement pattern case selected is not among the options provided: even_cells, inner_region') 
             end select
          end if
          call cell%next()
       end do
       call this%triangulation%free_cell_iterator(cell)
       deallocate(cell_coordinates, stat=istat); check(istat==0)
    end if
    
    if (allocated(is_node_coord_within_inner_region)) call memfree(is_node_coord_within_inner_region)
  end subroutine set_cells_for_refinement
  
  subroutine set_cells_for_uniform_refinement(this)
    implicit none
    class(par_test_h_adaptive_lagrangian_gp_fe_driver_t), intent(inout) :: this
    class(cell_iterator_t), allocatable :: cell
    class(environment_t), pointer :: environment
    environment => this%triangulation%get_environment()
    if ( environment%am_i_l1_task() ) then
       call this%triangulation%create_cell_iterator(cell)
       do while ( .not. cell%has_finished() )
          if ( cell%is_local() ) then       
              call cell%set_for_refinement()
          end if
          call cell%next()
       end do
       call this%triangulation%free_cell_iterator(cell)
    end if
  end subroutine set_cells_for_uniform_refinement
  
  
  subroutine set_cells_weights(this)
    implicit none
    class(par_test_h_adaptive_lagrangian_gp_fe_driver_t), intent(inout) :: this
    class(cell_iterator_t), allocatable :: cell
    class(environment_t), pointer :: environment
    environment => this%triangulation%get_environment()
    if ( environment%am_i_l1_task() ) then
      call this%triangulation%create_cell_iterator(cell)
      do while ( .not. cell%has_finished() )
        if ( cell%is_local() ) then
          call cell%set_weight(int(mod(cell%get_ggid(),2)+1,ip))
        end if   
        call cell%next()
      end do
      call this%triangulation%free_cell_iterator(cell)
    end if
  end subroutine set_cells_weights
  
  subroutine set_cells_for_coarsening(this)
    implicit none
    class(par_test_h_adaptive_lagrangian_gp_fe_driver_t), intent(inout) :: this
    class(cell_iterator_t)      , allocatable :: cell
    class(environment_t), pointer :: environment
    environment => this%triangulation%get_environment()
    if ( environment%am_i_l1_task() ) then
      call this%triangulation%create_cell_iterator(cell)
      do while ( .not. cell%has_finished() )
        if ( cell%is_local() ) then
          call cell%set_for_coarsening()
        end if  
        call cell%next()
      end do
      call this%triangulation%free_cell_iterator(cell)
    end if
  end subroutine set_cells_for_coarsening
  
      !========================================================================================
  subroutine print_info (this)
    implicit none
    class(par_test_h_adaptive_lagrangian_gp_fe_driver_t), intent(inout) :: this

    integer(ip) :: num_sub_domains
    real(rp) :: num_dofs
    real(rp) :: num_interface_dofs 
    integer(ip) :: num_coarse_dofs

    class(environment_t), pointer :: environment
    class(coarse_fe_space_t), pointer :: coarse_fe_space
    class(execution_context_t), pointer :: w_context

    environment => this%fe_space%get_environment()

    if (environment%am_i_l1_task()) then
      num_dofs             = real(this%fe_space%get_field_num_dofs(1),kind=rp)
      num_interface_dofs   = real(this%fe_space%get_total_num_interface_dofs(),kind=rp)
      call environment%l1_sum(num_dofs        )
      call environment%l1_sum(num_interface_dofs )
    end if

    w_context => environment%get_w_context()
    call w_context%barrier()
    
    if (environment%get_l1_rank() == 0) then
      num_sub_domains = environment%get_l1_size()
      write(*,'(a35,i22)') 'num_sub_domains:', num_sub_domains
      write(*,'(a35,i22)') 'num_total_cells:', this%triangulation%get_num_global_cells()
      write(*,'(a35,i22)') 'num_dofs (sub-assembled):', nint(num_dofs            , kind=ip )
      write(*,'(a35,i22)') 'num_interface_dofs (sub-assembled):', nint(num_interface_dofs  , kind=ip )
    end if
 
#ifdef ENABLE_MKL    
    if (environment%am_i_lgt1_task()) then
        coarse_fe_space => this%fe_space%get_coarse_fe_space()
        num_coarse_dofs = coarse_fe_space%get_field_num_dofs(1)
        write(*,'(a35,i22)') 'num_coarse_dofs:', num_coarse_dofs
    end if
#endif    

  end subroutine print_info
  
  
end module par_test_h_adaptive_lagrangian_gp_fe_driver_names
