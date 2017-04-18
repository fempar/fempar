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

module unfitted_l1_coarse_fe_handler_names
  use fe_space_names
  use par_sparse_matrix_names
  use FPL
  use types_names
  use matrix_names
  use base_sparse_matrix_names
  use par_scalar_array_names
  use environment_names
  use dof_import_names
  use serial_scalar_array_names

  implicit none
# include "debug.i90"
  private

!========================================================================================
  type, extends(standard_l1_coarse_fe_handler_t) :: unfitted_l1_coarse_fe_handler_t
    private

    ! A `par_sparse_matrix_t` should be sufficient to handle several fields,
    ! but not to handle several blocks
    class(par_sparse_matrix_t), pointer :: matrix => null()
    !real(rp), allocatable :: flag_ldofs(:)
    !logical :: is_set_up = .false.

  contains

    ! Public TBPs
    procedure :: create                   => unfitted_l1_create  
    procedure :: free                     => unfitted_l1_free
    !procedure :: get_num_coarse_dofs      => unfitted_l1_get_num_coarse_dofs
    !procedure :: setup_constraint_matrix  => unfitted_l1_setup_constraint_matrix
    procedure :: setup_weighting_operator => unfitted_l1_setup_weighting_operator

    !! Provate TBPs
    !procedure,private, non_overridable          :: allocate_and_fill_flag_ldofs  => &
    !  unfitted_l1_allocate_and_fill_flag_ldofs  
    !procedure, private, non_overridable, nopass :: fill_is_cut_cell => unfitted_l1_fill_is_cut_cell

  end type unfitted_l1_coarse_fe_handler_t

  public :: unfitted_l1_coarse_fe_handler_t

contains

!========================================================================================
subroutine unfitted_l1_create(this, matrix)
  implicit none
  class(unfitted_l1_coarse_fe_handler_t), intent(inout) :: this
  class(matrix_t), target,     intent(in)    :: matrix
  call this%free()
  select type (matrix)
  class is (par_sparse_matrix_t)
    this%matrix => matrix
  class default
    check(.false.)
  end select
end subroutine unfitted_l1_create

!========================================================================================
subroutine unfitted_l1_free(this)
  implicit none
  class(unfitted_l1_coarse_fe_handler_t), intent(inout) :: this
  !call memfree(this%flag_ldofs,__FILE__,__LINE__)
  this%matrix => null()
  !this%is_set_up = .false.
end subroutine unfitted_l1_free


!!========================================================================================
!subroutine unfitted_l1_get_num_coarse_dofs(this, par_fe_space, parameter_list, num_coarse_dofs)
!  implicit none
!  class(unfitted_l1_coarse_fe_handler_t), intent(in)    :: this
!  type(par_fe_space_t)                  , intent(in)    :: par_fe_space
!  type(parameterlist_t)                 , intent(in)    :: parameter_list
!  integer(ip)                           , intent(inout) :: num_coarse_dofs(:)
!  call this%allocate_and_fill_flag_ldofs(par_fe_space,parameter_list)
!end subroutine unfitted_l1_get_num_coarse_dofs

!!========================================================================================
!subroutine unfitted_l1_setup_constraint_matrix(this, par_fe_space, parameter_list, constraint_matrix) 
!  implicit none
!  class(unfitted_l1_coarse_fe_handler_t), intent(in)    :: this
!  type(par_fe_space_t)                  , intent(in)    :: par_fe_space
!  type(parameterlist_t)                 , intent(in)    :: parameter_list
!  type(coo_sparse_matrix_t)             , intent(inout) :: constraint_matrix
!end subroutine unfitted_l1_setup_constraint_matrix

!========================================================================================
subroutine unfitted_l1_setup_weighting_operator(this, par_fe_space, parameter_list, weighting_operator) 
  implicit none
  class(unfitted_l1_coarse_fe_handler_t), intent(in)    :: this
  type(par_fe_space_t)                  , intent(in)    :: par_fe_space
  type(parameterlist_t)                 , intent(in)    :: parameter_list
  real(rp), allocatable                 , intent(inout) :: weighting_operator(:)

  type(par_scalar_array_t) :: par_array
  type(environment_t), pointer :: p_env
  type(dof_import_t),  pointer :: dof_import
  type(serial_scalar_array_t), pointer :: serial_array
  real(rp), pointer     :: assembled_diag(:)
  real(rp), allocatable :: sub_assembled_diag(:)
  integer(ip) :: istat
  integer(ip) :: block_id

  ! We assume a single block (for the moment)
  block_id = 1

  ! Get the sub-assembled diagonal
  call this%matrix%extract_diagonal(sub_assembled_diag)

  ! Communicate to compute the fully assembled diagonal
  p_env => par_fe_space%get_environment()
  dof_import => par_fe_space%get_block_dof_import(block_id)
  call par_array%create_and_allocate(p_env, dof_import)
  serial_array   => par_array%get_serial_scalar_array()
  assembled_diag => serial_array%get_entries()
  assembled_diag(:) = sub_assembled_diag(:)
  call par_array%comm()

  ! Compute the weighting
  if (allocated(weighting_operator) ) then
    call memfree ( weighting_operator, __FILE__, __LINE__ )
  end if
  call memalloc(size(sub_assembled_diag),weighting_operator, __FILE__, __LINE__ )
  weighting_operator(:) = sub_assembled_diag(:)/assembled_diag(:)

  ! Clean up
  deallocate(sub_assembled_diag,stat=istat); check(istat == 0)
  call par_array%free()

end subroutine unfitted_l1_setup_weighting_operator

!!========================================================================================
!subroutine unfitted_l1_allocate_and_fill_flag_ldofs(this,par_fe_space,parameter_list)
!  implicit none
!  class(unfitted_l1_coarse_fe_handler_t), intent(inout) :: this
!  type(par_fe_space_t)                  , intent(in)    :: par_fe_space
!  type(parameterlist_t)                 , intent(in)    :: parameter_list
!
!  type(environment_t), pointer           :: par_environment
!  integer(ip)                            :: field_id, block_id
!  integer(ip), pointer                   :: field_to_block(:)
!  integer(ip)                            :: num_lfdofs
!  class(base_static_triangulation_t), pointer :: triangulation
!  type(cell_import_t), pointer           :: cell_import
!  integer(ip)                            :: max_num_parts_arround
!  logical, allocatable                   :: visited_lparts(:)
!  integer(ip), allocatable               :: l2g_parts(:)
!  integer(ip)                            :: num_total_cells
!  logical                                :: is_cut_cell(:)
!  logical                                :: use_vertices, use_edges, use_faces
!  integer(ip)                            :: icell_around
!  integer(ip)                            :: ivef
!  integer(ip)                            :: idof, dof_lid
!  logical                                :: dofs_on_vef
!  type(fe_object_iterator_t)             :: object_iterator
!  type(fe_object_accessor_t)             :: object
!  type(fe_vefs_on_object_iterator_t)     :: vefs_on_object_iterator
!  type(fe_vef_accessor_t)                :: vef
!  type(fe_accessor_t)                    :: fe
!  type(list_iterator_t)                  :: own_dofs_on_vef_iterator
!  integer(ip), pointer                   :: elem2dof(:)
!  integer(ip)                            :: mypart_id, part_id, local_part_id
!  logical                                :: is_the_smallest
!  integer(ip)                            :: new_corner_id
!  logical(ip)                            :: vef_touches_cut_cell
!
!  par_environment   => par_fe_space%get_par_environment()
!  assert ( par_environment%am_i_l1_task() )
!  assert ( associated ( par_environment ) )
!
!  ! We assume a single field for the moment
!  field_id = 1
!  assert(par_fe_space%get_number_fields() == 1)
!
!  ! Allocate the local dofs flag array
!  field_to_block => par_fe_space%get_field_block()
!  block_id = field_to_block(field_id)
!  num_lfdofs = par_fe_space%get_block_number_dofs(block_id)
!  call memalloc(num_lfdofs,this%flag_ldofs,__FILE__,__LINE__)
!
!  ! We mark all the local dofs with a negative value
!  this%flag_ldofs(:) = -1 
!
!  ! Allocate work vector for tracking visited parts
!  triangulation     => par_fe_space%get_triangulation()
!  cell_import       => triangulation%get_cell_import()
!  max_num_parts_arround = cell_import%get_number_neighbours()
!  call memalloc(max_num_parts_arround,visited_lparts,__FILE__,__LINE__)
!  call memalloc(max_num_parts_arround,l2g_parts,__FILE__,__LINE__)
!
!  ! Compute a vector of flags identifying cut/non cut cells for both local and ghost cells
!  ! This requires communication
!  num_total_cells = triangulation%get_num_local_cells() + triangulation%get_num_ghost_cells()
!  call memalloc(num_total_cells,is_cut_cell,__FILE__,__LINE__)
!  call this%fill_is_cut_cell(triangulation,is_cut_cell)
!
!  ! Recover options
!  call this%get_coarse_space_use_vertices_edges_faces(parameter_list,& 
!                                                      use_vertices, &
!                                                      use_edges, &
!                                                      use_faces)
!  ! Loop in objects
!  object_iterator = par_fe_space%create_fe_object_iterator()
!  do while ( .not. object_iterator%has_finished() )
!     call object_iterator%current(object)
!
!    ! Check if c, ce, or cef, and skip object accordingly
!    select case ( object%get_dimension() )
!    case (0)
!      if (.not. use_vertices) then
!        call object_iterator%next(); cycle
!      end if  
!    case (1)
!      if (.not. use_edges) then
!        call object_iterator%next(); cycle
!      end if  
!    case (2)
!      if (.not. use_faces) then
!        call object_iterator%next(); cycle
!      end if  
!    end select
!
!    ! First stage ----------------
!    visited_lparts(:) = .false.
!    l2g_parts(:) = 0
!    
!    ! Loop in vefs of the object
!    vefs_on_object_iterator = object%create_fe_vefs_on_object_iterator()
!    do while ( .not. vefs_on_object_iterator%has_finished() )
!      call vefs_on_object_iterator%current(vef)
!
!      ! Loop in ghost cells around the vef
!      do icell_around=1, vef%get_num_cells_around()
!        call vef%get_cell_around(icell_around,fe)
!        if ( fe%is_ghost() ) then 
!
!          part_id = fe%get_my_part()
!          local_part_id = cell_import%get_local_neighbour_id(part_id)
!           
!          ! Loop in own dofs in the vef as seen from the ghost element
!          call fe%get_field_elem2dof(field_id, elem2dof)
!          ivef = fe%find_lpos_vef_lid(vef%get_lid())
!          own_dofs_on_vef_iterator = fe%create_own_dofs_on_vef_iterator(ivef, field_id)
!          do while ( .not. own_dofs_on_vef_iterator%is_upper_bound() )
!            idof    = own_dofs_on_vef_iterator%get_current()
!            dof_lid = elem2dof(idof)
!
!            ! Mark the dofs of the vef with a 0 and mark the current part as visited
!            if ( dof_lid > 0 ) then
!              flag_ldofs(dof_lid) = 0
!              visited_lparts(local_part_id) = .true.
!              l2g_parts(local_part_id) = part_id
!            end if
!
!             call own_dofs_on_vef_iterator%next()
!          end do
!           
!        end if
!      end do
!
!      call vefs_on_object_iterator%next()
!    end do
!
!    ! At this point all the dofs of the current object that are also in another subdomain are marked with 0
!
!    ! Second stage --------------
!
!    ! (only if the object is an edge)
!    if ( object%get_dimension() .ne. 1 )
!      call object_iterator%next(); cycle
!    end if
!
!    ! (and only if the current part id is smaller than all the visited parts)
!   is_the_smallest = .true.
!   mypart_id = par_environment%get_l1_rank() + 1
!   do local_part_id = 1,cell_import%get_number_neighbours()
!     if (visited_lparts(local_part_id)) then
!       if ( l2g_parts(local_part_id) < mypart_id) is_the_smallest = .false.
!     end if
!   end do
!   if (.not. is_the_smallest) then
!      call object_iterator%next(); cycle
!   end if
!
!    ! Initialize to 1 the counter of new corners
!    new_corner_id = 1
!
!    ! Loop in vefs of the object
!    vefs_on_object_iterator = object%create_fe_vefs_on_object_iterator()
!    do while ( .not. vefs_on_object_iterator%has_finished() )
!      call vefs_on_object_iterator%current(vef)
!
!      ! Check if the vef touches a cut cell (local or ghost)
!      vef_touches_cut_cell = .false.
!      do icell_around=1, vef%get_num_cells_around()
!        call vef%get_cell_around(icell_around,fe)
!        if ( is_cut_cell(fe%get_lid()) ) then
!          vef_touches_cut_cell = .true.
!          exit
!        end if
!      end do
!
!      ! Only if the current vef touches a cut cell:
!      if (.not. vef_touches_cut_cell) then
!        call vefs_on_object_iterator%next(); cycle
!      end if
!    
!      ! Loop in ghost cells around the vef
!      do icell_around=1, vef%get_num_cells_around()
!        call vef%get_cell_around(icell_around,fe)
!        if ( fe%is_local() ) then 
!           
!          ! Loop in own dofs in the vef as seen from the ghost cell
!          call fe%get_field_elem2dof(field_id, elem2dof)
!          ivef = fe%find_lpos_vef_lid(vef%get_lid())
!          own_dofs_on_vef_iterator = fe%create_own_dofs_on_vef_iterator(ivef, field_id)
!          do while ( .not. own_dofs_on_vef_iterator%is_upper_bound() )
!            idof    = own_dofs_on_vef_iterator%get_current()
!            dof_lid = elem2dof(idof)
!
!            if ( dof_lid > 0 ) then
!
!              ! If flag corresponding to the current dof is 0
!              ! (i.e. a dof on the interface that is not yet marked as a new corner),
!              if (flag_ldofs(dof_lid)==0) then
!
!                ! Mark this dof with the counter of new corners
!                flag_ldofs(dof_lid) = new_corner_id
!
!                ! Increment the counter of new corners
!                new_corner_id = new_corner_id + 1
!
!              end if
!            end if
!
!             call own_dofs_on_vef_iterator%next()
!          end do
!           
!        end if
!      end do
!
!      call vefs_on_object_iterator%next()
!    end do
!
!    call object_iterator%next()
!  end do
!
!  ! Communicate the value of the flag
!
!  ! At this point we have identified all the coarse dofs inside each object
!  ! A dof is classified in function of the flag as follows:
!  ! flag == -1 the dof is interior (is not in any other subdomain)
!  ! flag == 0 dof is in another subdomain but it is not a new corner
!  ! flag > 0  dof is in another subdomain and it is a new corner
!
!
!  ! Loop in objects
!
!    ! Recover the dofs of the object
!
!
!
!
!
!
!
!  ! Clean up
!  call memfree(visited_lparts,__FILE__,__LINE__)
!  call memfree(is_cut_cell   ,__FILE__,__LINE__)
!  call memfree(l2g_parts     ,__FILE__,__LINE__)
!
!  ! Mark as set up
!  this%is_set_up = .true.
!
!end subroutine unfitted_l1_setup

!!========================================================================================
!subroutine unfitted_l1_fill_is_cut_cell(triangulation,is_cut_cell)
!  implicit none
!  class(base_static_triangulation_t), pointer, intent(in) :: triangulation
!  logical, allocatable, intent(inout) :: is_cut_cell
!  !TODO
!  mcheck(.false.,"Not yet implemented")
!end subroutine unfitted_l1_fill_is_cut_cell

end module unfitted_l1_coarse_fe_handler_names
!***************************************************************************************************
