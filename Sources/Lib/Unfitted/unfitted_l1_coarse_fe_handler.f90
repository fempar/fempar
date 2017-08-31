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
  use fe_affine_operator_names
  use list_types_names
  use triangulation_names
  use unfitted_triangulations_names
  use unfitted_fe_spaces_names
  use cell_import_names
  use stiffness_weighting_l1_coarse_fe_handler_names

  implicit none
# include "debug.i90"
  private

!========================================================================================
  type, extends(stiffness_weighting_l1_coarse_fe_handler_t) :: unfitted_l1_coarse_fe_handler_t
    private

    class(par_fe_space_t),      pointer     :: par_fe_space   => null()
    class(parameterlist_t),     pointer     :: parameter_list => null()

    type(list_t)                            :: object_gid_to_dof_gids
    integer(ip),                allocatable :: dof_gid_to_cdof_id_in_object(:)
    integer(ip),                allocatable :: object_gid_to_min_neighbour(:)

  contains

    ! Creation / Deletion methods
    generic :: create                         => unfitted_l1_create
    procedure, private :: unfitted_l1_create
    procedure, private :: stiffness_l1_create => unfitted_l1_stiffness_l1_create
    procedure :: free                         => unfitted_l1_free

    ! Overwritten TPBs
    procedure :: get_num_coarse_dofs          => unfitted_l1_get_num_coarse_dofs
    procedure :: setup_constraint_matrix      => unfitted_l1_setup_constraint_matrix

    !! Private TBPs
    procedure, private, non_overridable :: setup_object_gid_to_dof_gids       => unfitted_l1_setup_object_gid_to_dof_gids
    procedure, private, non_overridable :: setup_dof_gid_to_cdof_id_in_object => unfitted_l1_setup_dof_gid_to_cdof_id_in_object
    procedure, private, non_overridable :: identify_problematic_dofs          => unfitted_l1_identify_problematic_dofs

  end type unfitted_l1_coarse_fe_handler_t

  public :: unfitted_l1_coarse_fe_handler_t

contains

!========================================================================================
subroutine unfitted_l1_create(this, fe_affine_operator, parameter_list)

  implicit none
  class(unfitted_l1_coarse_fe_handler_t), intent(inout) :: this
  class(fe_affine_operator_t), target,    intent(in)    :: fe_affine_operator
  class(parameterlist_t),      target,    intent(in)    :: parameter_list

  class(matrix_t),            pointer :: matrix
  class(serial_fe_space_t),   pointer :: fe_space
  type(environment_t),        pointer :: environment

  call this%free()

  call this%stiffness_weighting_l1_coarse_fe_handler_t%create(fe_affine_operator)

  fe_space => fe_affine_operator%get_fe_space()
  select type (fe_space)
    class is (par_fe_space_t)
      this%par_fe_space => fe_space
    class default
      check(.false.)
  end select

  this%parameter_list => parameter_list

  environment => this%par_fe_space%get_environment()
  assert (associated(environment))

  if (environment%am_i_l1_task()) then
    call this%setup_object_gid_to_dof_gids()
    call this%setup_dof_gid_to_cdof_id_in_object()
  end if

end subroutine unfitted_l1_create

!========================================================================================
subroutine unfitted_l1_stiffness_l1_create(this, fe_affine_operator)
  implicit none
  class(unfitted_l1_coarse_fe_handler_t), intent(inout) :: this
  class(fe_affine_operator_t), target,    intent(in)    :: fe_affine_operator
  mcheck(.false.,'This method does not make sense for this class')
end subroutine unfitted_l1_stiffness_l1_create

!========================================================================================
subroutine unfitted_l1_free(this)
  implicit none
  class(unfitted_l1_coarse_fe_handler_t), intent(inout) :: this
  call this%stiffness_weighting_l1_coarse_fe_handler_t%free()
  this%par_fe_space   => null()
  this%parameter_list => null()
  if ( allocated(this%dof_gid_to_cdof_id_in_object) ) call memfree(this%dof_gid_to_cdof_id_in_object,__FILE__,__LINE__)
  if ( allocated(this%object_gid_to_min_neighbour) ) call memfree(this%object_gid_to_min_neighbour,__FILE__,__LINE__)
  call this%object_gid_to_dof_gids%free()
end subroutine unfitted_l1_free

!========================================================================================
subroutine unfitted_l1_get_num_coarse_dofs(this,field_id,par_fe_space,parameter_list,num_coarse_dofs)
  implicit none
  class(unfitted_l1_coarse_fe_handler_t), intent(in)    :: this
  integer(ip)                           , intent(in)    :: field_id
  type(par_fe_space_t)                  , intent(in)    :: par_fe_space
  type(parameterlist_t)                 , intent(in)    :: parameter_list
  integer(ip)                           , intent(inout) :: num_coarse_dofs(:)

  type(environment_t), pointer  :: environment
  integer(ip)                   :: max_cdof_gid_in_object
  logical, allocatable          :: visited_cdof_gids_in_object(:)
  type(list_iterator_t)         :: dofs_in_object_iterator
  integer(ip)                   :: dof_gid
  type(fe_object_iterator_t)    :: object

  assert(field_id == 1)

  environment => this%par_fe_space%get_environment()
  assert ( associated ( environment ) )
  assert ( environment%am_i_l1_task() )
  assert ( size(num_coarse_dofs) == this%par_fe_space%get_num_fe_objects() )


  max_cdof_gid_in_object = maxval(this%dof_gid_to_cdof_id_in_object) + 1
  call memalloc(max_cdof_gid_in_object,visited_cdof_gids_in_object,__FILE__,__LINE__)

  ! Loop in objects
  call this%par_fe_space%create_fe_object_iterator(object)
  do while ( .not. object%has_finished() )

    ! Count how many coarse dofs are on this object:
    ! i.e., loop on the local dofs of the object and count how many different numbers are found
    visited_cdof_gids_in_object(:) = .false.
    dofs_in_object_iterator = this%object_gid_to_dof_gids%create_iterator(object%get_gid())
    do while (.not. dofs_in_object_iterator%is_upper_bound())
      dof_gid = dofs_in_object_iterator%get_current()
      visited_cdof_gids_in_object(this%dof_gid_to_cdof_id_in_object(dof_gid)+1) = .true.
      call dofs_in_object_iterator%next()
    end do
    num_coarse_dofs(object%get_gid()) =  count(visited_cdof_gids_in_object)

    call object%next()
  end do

  call memfree(visited_cdof_gids_in_object,__FILE__,__LINE__)
  call this%par_fe_space%free_fe_object_iterator(object)

end subroutine unfitted_l1_get_num_coarse_dofs

!========================================================================================
subroutine unfitted_l1_setup_constraint_matrix(this,field_id,par_fe_space,parameter_list,constraint_matrix)
  implicit none
  class(unfitted_l1_coarse_fe_handler_t), intent(in)    :: this
  integer(ip)                           , intent(in)    :: field_id
  type(par_fe_space_t)                  , intent(in)    :: par_fe_space
  type(parameterlist_t)                 , intent(in)    :: parameter_list
  type(coo_sparse_matrix_t)             , intent(inout) :: constraint_matrix

  type(environment_t), pointer :: environment
  integer(ip)                  :: block_id
  integer(ip),         pointer :: field_to_block(:)
  integer(ip)                  :: num_cols, num_rows
  integer(ip)                  :: max_cdof_gid_in_object
  type(list_iterator_t)        :: dofs_in_object_iterator
  integer(ip)                  :: dof_gid, cdof_gid_in_object
  integer(ip)                  :: num_fdofs_in_cdof
  integer(ip)                  :: cdof_gid
  type(fe_object_iterator_t)   :: object

  environment => this%par_fe_space%get_environment()
  assert (associated(environment))
  assert (environment%am_i_l1_task())

  ! We assume a single field for the moment
  assert(field_id == 1)
  assert(this%par_fe_space%get_num_fields() == 1)

  ! Free any dynamic memory that constraint_matrix may have inside
  call constraint_matrix%free()

  ! Create constraint matrix (transposed)
  field_to_block => this%par_fe_space%get_field_blocks()
  block_id = field_to_block(field_id)
  num_cols = this%par_fe_space%get_block_num_dofs(block_id)
  num_rows = this%par_fe_space%get_block_num_coarse_dofs(block_id)
  call constraint_matrix%create ( num_cols, num_rows )

  cdof_gid = 0
  max_cdof_gid_in_object = maxval(this%dof_gid_to_cdof_id_in_object) + 1

  ! Loop in objects
  call this%par_fe_space%create_fe_object_iterator(object)
  do while ( .not. object%has_finished() )

    ! Loop in all possible values of cdof_gid_in_object
    ! TODO Maybe these are too many loops
    do cdof_gid_in_object = 0,max_cdof_gid_in_object

      ! Count how many fine dofs are on this c dof
      ! i.e., do a loop in the fine dofs on this object and count how many time the current
      ! coarse dof is found
      num_fdofs_in_cdof = 0
      dofs_in_object_iterator = this%object_gid_to_dof_gids%create_iterator(object%get_gid())
      do while (.not. dofs_in_object_iterator%is_upper_bound())
        dof_gid = dofs_in_object_iterator%get_current()
        if ( cdof_gid_in_object == this%dof_gid_to_cdof_id_in_object(dof_gid) ) num_fdofs_in_cdof = num_fdofs_in_cdof + 1
        call dofs_in_object_iterator%next()
      end do

      ! If there are 0 then cycle
      if (num_fdofs_in_cdof == 0) cycle

      ! Increment in 1 the current row
      cdof_gid = cdof_gid + 1

      ! Loop in all fine dofs of this object
        ! If a fine dof is in the current cdof_gid_in_object add 1/num_fdofs_in_cdof in the constaint matrix
      dofs_in_object_iterator = this%object_gid_to_dof_gids%create_iterator(object%get_gid())
      do while (.not. dofs_in_object_iterator%is_upper_bound())
        dof_gid = dofs_in_object_iterator%get_current()
        if ( cdof_gid_in_object == this%dof_gid_to_cdof_id_in_object(dof_gid) ) then
          call constraint_matrix%insert(dof_gid,cdof_gid, 1.0_rp/real(num_fdofs_in_cdof,rp))
        end if
        call dofs_in_object_iterator%next()
      end do

    end do

    call object%next()
  end do
  call constraint_matrix%sort_and_compress()

  call this%par_fe_space%free_fe_object_iterator(object)

  !call this%standard_l1_coarse_fe_handler_t%setup_constraint_matrix(par_fe_space,parameter_list,constraint_matrix)

end subroutine unfitted_l1_setup_constraint_matrix

!========================================================================================
subroutine unfitted_l1_setup_object_gid_to_dof_gids(this)

  implicit none
  class(unfitted_l1_coarse_fe_handler_t), intent(inout) :: this

  integer(ip)                        :: num_objects
  integer(ip), allocatable           :: object_gid_to_num_dofs_in_object(:)
  integer(ip)                        :: icell_around
  integer(ip)                        :: ivef_within_cell
  integer(ip)                        :: ivef_within_object
  integer(ip)                        :: idof, dof_gid
  logical                            :: dofs_on_vef
  type(environment_t), pointer       :: environment
  type(fe_object_iterator_t)         :: object
  type(fe_vef_iterator_t)            :: vef
  class(fe_cell_iterator_t), allocatable  :: fe
  type(list_iterator_t)              :: own_dofs_on_vef_iterator
  integer(ip), pointer               :: fe_dofs(:)
  logical                            :: use_vertices, use_edges, use_faces
  logical, allocatable               :: visited_dofs(:)
  integer(ip)                        :: field_id, block_id
  integer(ip), pointer               :: field_to_block(:)
  integer(ip)                        :: num_dofs, num_dofs_in_object
  integer(ip), allocatable           :: dof_gids(:)
  integer(ip)                        :: icount
  type(list_iterator_t)              :: dofs_in_object_iterator

  ! We assume a single field for the moment
  field_id = 1
  assert(this%par_fe_space%get_num_fields() == 1)

  ! Auxiliary
  field_to_block => this%par_fe_space%get_field_blocks()
  block_id = field_to_block(field_id)
  num_dofs  = this%par_fe_space%get_block_num_dofs(block_id)
  num_objects = this%par_fe_space%get_num_fe_objects()
  call memalloc(num_dofs,visited_dofs,__FILE__,__LINE__)
  call memalloc(num_dofs,dof_gids,__FILE__,__LINE__)
  call memalloc(num_objects,object_gid_to_num_dofs_in_object,__FILE__,__LINE__)
  icount = 0

  ! Allocate array of min neighbors and initialize with the biggest integer
  call memalloc(num_objects,this%object_gid_to_min_neighbour,__FILE__,__LINE__)
  this%object_gid_to_min_neighbour(:) = huge(num_objects)

  call this%get_coarse_space_use_vertices_edges_faces(this%parameter_list,& 
                                                      use_vertices, &
                                                      use_edges, &
                                                      use_faces)

  ! Fill the auxiliary data
  call this%par_fe_space%create_fe_object_iterator(object)
  call this%par_fe_space%create_fe_cell_iterator(fe)
  call this%par_fe_space%create_fe_vef_iterator(vef)
  do while ( .not. object%has_finished() )

    ! Check if c, ce, or cef, and skip object accordingly
    select case ( object%get_dim() )
    case (0)
      if (.not. use_vertices) then
        call object%next(); cycle
      end if
    case (1)
      if (.not. use_edges) then
        call object%next(); cycle
      end if
    case (2)
      if (.not. use_faces) then
        call object%next(); cycle
      end if
    end select

    ! Loop in vefs of the object
    visited_dofs(:) = .false.
    num_dofs_in_object = 0
    do ivef_within_object=1, object%get_num_vefs()
       call object%get_vef(ivef_within_object,vef)

      ! Loop in ghost cells around the vef
      do icell_around=1, vef%get_num_cells_around()
        call vef%get_cell_around(icell_around,fe)
        if ( fe%is_ghost() ) then

          call fe%get_field_fe_dofs(field_id, fe_dofs)
          ivef_within_cell = fe%get_vef_lid_from_gid(vef%get_gid())

          ! Loop in own dofs in the vef as seen from the ghost element
          own_dofs_on_vef_iterator = fe%create_own_dofs_on_vef_iterator(ivef_within_cell, field_id)
          do while ( .not. own_dofs_on_vef_iterator%is_upper_bound() )
            idof    = own_dofs_on_vef_iterator%get_current()
            dof_gid = fe_dofs(idof)
            if ( dof_gid > 0 ) then
              if(.not. visited_dofs(dof_gid)) then

                ! Store the minimum (active) neighbor part
                if ( this%object_gid_to_min_neighbour(object%get_gid()) > fe%get_my_part() ) then
                  this%object_gid_to_min_neighbour(object%get_gid()) = fe%get_my_part()
                end if

                ! Store the dofs on this object
                visited_dofs(dof_gid) = .true.
                icount = icount + 1
                num_dofs_in_object = num_dofs_in_object + 1
                dof_gids(icount) = dof_gid

              end if
            end if
            call own_dofs_on_vef_iterator%next()
          end do

        end if
      end do

    end do
    object_gid_to_num_dofs_in_object(object%get_gid()) = num_dofs_in_object
    call object%next()
  end do

  ! Initialize the list
  call this%object_gid_to_dof_gids%free()
  call this%object_gid_to_dof_gids%create(this%par_fe_space%get_num_fe_objects())
  call object%first()
  do while ( .not. object%has_finished() )
    call this%object_gid_to_dof_gids%sum_to_pointer_index(&
      object%get_gid(),object_gid_to_num_dofs_in_object(object%get_gid()))
    call object%next()
  end do
  call this%object_gid_to_dof_gids%calculate_header()
  call this%object_gid_to_dof_gids%allocate_list_from_pointer()

  ! Fill the list
  icount = 0
  call object%first()
  do while ( .not. object%has_finished() )

    dofs_in_object_iterator = this%object_gid_to_dof_gids%create_iterator(object%get_gid())
    do while (.not. dofs_in_object_iterator%is_upper_bound())
      icount = icount + 1
      call dofs_in_object_iterator%set_current(dof_gids(icount))
      call dofs_in_object_iterator%next()
    end do

    call object%next()
  end do

  ! Clean up
  call memfree(object_gid_to_num_dofs_in_object,__FILE__,__LINE__)
  call memfree(dof_gids,__FILE__,__LINE__)
  call memfree(visited_dofs,__FILE__,__LINE__)
  call this%par_fe_space%free_fe_vef_iterator(vef)
  call this%par_fe_space%free_fe_cell_iterator(fe)
  call this%par_fe_space%free_fe_object_iterator(object)

end subroutine unfitted_l1_setup_object_gid_to_dof_gids

!========================================================================================
subroutine unfitted_l1_setup_dof_gid_to_cdof_id_in_object(this)

  implicit none
  class(unfitted_l1_coarse_fe_handler_t), intent(inout) :: this

  logical, allocatable                 :: is_problematic_dof(:)
  type(fe_object_iterator_t)           :: object
  type(list_iterator_t)                :: dofs_in_object_iterator
  integer(ip)                          :: dof_gid
  integer(ip)                          :: field_id, block_id
  integer(ip), pointer                 :: field_to_block(:)
  integer(ip)                          :: num_dofs
  integer(ip)                          :: new_corner_counter
  integer(ip)                          :: my_part_id
  type(par_scalar_array_t)             :: par_array
  type(environment_t), pointer         :: p_env
  type(dof_import_t),  pointer         :: dof_import
  type(serial_scalar_array_t), pointer :: serial_array
  real(rp), pointer                    :: array_entries(:)

  ! We assume a single field for the moment
  field_id = 1
  assert(this%par_fe_space%get_num_fields() == 1)
  
  ! Get my part id
  p_env => this%par_fe_space%get_environment()
  my_part_id = p_env%get_l1_rank() + 1

  ! Allocate and initialize the member variable
  field_to_block => this%par_fe_space%get_field_blocks()
  block_id = field_to_block(field_id)
  num_dofs = this%par_fe_space%get_block_num_dofs(block_id)
  call memalloc(num_dofs,this%dof_gid_to_cdof_id_in_object,__FILE__,__LINE__)
  this%dof_gid_to_cdof_id_in_object(:) = 0

  ! Identify which dofs are problematic
  call memalloc(num_dofs,is_problematic_dof,__FILE__,__LINE__)
  call this%identify_problematic_dofs(is_problematic_dof)

  ! Loop in objects
  call this%par_fe_space%create_fe_object_iterator(object)
  do while ( .not. object%has_finished() )

    ! Skip objects that are not edges
    if (object%get_dim() .ne. 1) then
      call object%next(); cycle
    end if

    ! Skip if we are not the part with minimum id
    assert( this%object_gid_to_min_neighbour(object%get_gid()) .ne. my_part_id )
    if ( this%object_gid_to_min_neighbour(object%get_gid()) < my_part_id ) then
      call object%next(); cycle
    end if

    ! Loop in fine dofs on this object
      ! If the dof is problematic, mark it as a new corner
    new_corner_counter = 0
    dofs_in_object_iterator = this%object_gid_to_dof_gids%create_iterator(object%get_gid())
    do while (.not. dofs_in_object_iterator%is_upper_bound())
      dof_gid = dofs_in_object_iterator%get_current()
      if ( is_problematic_dof(dof_gid) ) then
        new_corner_counter = new_corner_counter + 1
        this%dof_gid_to_cdof_id_in_object(dof_gid) = new_corner_counter
      end if
      call dofs_in_object_iterator%next()
    end do

    call object%next()
  end do
  call this%par_fe_space%free_fe_object_iterator(object)

  ! Communicate to make it consistent between sub-domains
  dof_import => this%par_fe_space%get_block_dof_import(block_id)
  call par_array%create_and_allocate(p_env, dof_import)
  serial_array   => par_array%get_serial_scalar_array()
  array_entries => serial_array%get_entries()
  array_entries(:) = real(this%dof_gid_to_cdof_id_in_object(:),kind=rp)
  call par_array%comm()
  this%dof_gid_to_cdof_id_in_object(:) = nint(array_entries,kind=ip)

  call memfree(is_problematic_dof,__FILE__,__LINE__)
  call par_array%free()

end subroutine unfitted_l1_setup_dof_gid_to_cdof_id_in_object

!========================================================================================
subroutine unfitted_l1_identify_problematic_dofs(this,is_problematic_dof)

  implicit none
  class(unfitted_l1_coarse_fe_handler_t), intent(in)    :: this
  logical, allocatable,                   intent(inout) :: is_problematic_dof(:)

  integer(ip), allocatable :: is_cut_cell(:)
  integer(ip) :: field_id
  class(par_fe_space_t), pointer :: par_fe_space
  class(par_unfitted_fe_space_t), pointer :: par_unf_fe_space
  class(triangulation_t), pointer :: triangulation
  class(fe_cell_iterator_t), allocatable  :: fe
  integer(ip) :: num_total_cells
  type(environment_t), pointer :: environment
  type(cell_import_t), pointer :: cell_import
  integer(ip), pointer :: fe_dofs(:)

  ! We assume a single field for the moment
  field_id = 1
  assert(this%par_fe_space%get_num_fields() == 1)
  

  ! Recover the unfitted par fe space
  par_fe_space => this%par_fe_space
  select type(par_fe_space)
    class is (par_unfitted_fe_space_t)
      par_unf_fe_space => par_fe_space
    class default
      check(.false.)
  end select

  ! Mark the local cut elements
  triangulation => par_fe_space%get_triangulation()
  num_total_cells = triangulation%get_num_local_cells() + triangulation%get_num_ghost_cells()
  call memalloc(num_total_cells,is_cut_cell,__FILE__,__LINE__)
  is_cut_cell(:) = 0
  call par_unf_fe_space%create_fe_cell_iterator(fe)
  do while ( .not. fe%has_finished() )
    if ( fe%is_ghost() ) then
      call fe%next(); cycle
    end if
    if (fe%is_cut()) is_cut_cell(fe%get_gid()) = 1
    call fe%next()
  end do

  ! Communicate so that the ghost also have this info
  ! TODO this could be avoided if ghost cells have also the coordinates
  ! For the structured triangulation seems to be true, but for the unstructured?
  environment => triangulation%get_environment()
  cell_import => triangulation%get_cell_import()
  if(environment%get_l1_size()>1) &
  call environment%l1_neighbours_exchange ( cell_import%get_num_neighbours(), &
                                                cell_import%get_neighbours_ids(),    &
                                                cell_import%get_rcv_ptrs(),          &
                                                cell_import%get_rcv_leids(),         &
                                                cell_import%get_num_neighbours(), &
                                                cell_import%get_neighbours_ids(),    &
                                                cell_import%get_snd_ptrs(),          &
                                                cell_import%get_snd_leids(),         &
                                                is_cut_cell )

  ! Mark all the dofs belonging to cut elements
  is_problematic_dof(:) = .false.
  call fe%first()
  do while ( .not. fe%has_finished() )
    call fe%get_field_fe_dofs(field_id, fe_dofs)
    if (is_cut_cell(fe%get_gid())==1) is_problematic_dof(pack(fe_dofs,fe_dofs>0)) = .true.
    call fe%next()
  end do


  ! Clean up
  call memfree(is_cut_cell,__FILE__,__LINE__)
  call par_unf_fe_space%free_fe_cell_iterator(fe)

end subroutine unfitted_l1_identify_problematic_dofs

end module unfitted_l1_coarse_fe_handler_names
!***************************************************************************************************
