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
  use base_static_triangulation_names
  use unfitted_triangulations_names
  use unfitted_fe_spaces_names
  use cell_import_names

  implicit none
# include "debug.i90"
  private

!========================================================================================
  type, extends(standard_l1_coarse_fe_handler_t) :: unfitted_l1_coarse_fe_handler_t
    private

    class(par_sparse_matrix_t), pointer     :: matrix         => null()
    class(par_fe_space_t),      pointer     :: par_fe_space   => null()
    class(parameterlist_t),     pointer     :: parameter_list => null()

    real(rp),                   allocatable :: stiffness_weighting(:)
    type(list_t)                            :: object_lid_to_dof_lids
    integer(ip),                allocatable :: dof_lid_to_cdof_id_in_object(:)
    integer(ip),                allocatable :: object_lid_to_min_neigbour(:) 

  contains

    ! Public TBPs
    procedure :: create                   => unfitted_l1_create
    procedure :: free                     => unfitted_l1_free
    procedure :: get_num_coarse_dofs      => unfitted_l1_get_num_coarse_dofs
    procedure :: setup_constraint_matrix  => unfitted_l1_setup_constraint_matrix
    procedure :: setup_weighting_operator => unfitted_l1_setup_weighting_operator

    !! Private TBPs
    procedure, private, non_overridable :: setup_stiffness_weighting          => unfitted_l1_setup_stiffness_weighting
    procedure, private, non_overridable :: setup_object_lid_to_dof_lids       => unfitted_l1_setup_object_lid_to_dof_lids
    procedure, private, non_overridable :: setup_dof_lid_to_cdof_id_in_object => unfitted_l1_setup_dof_lid_to_cdof_id_in_object
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
  type(environment_t),        pointer :: par_environment

  call this%free()

  matrix => fe_affine_operator%get_matrix()
  select type (matrix)
    class is (par_sparse_matrix_t)
      this%matrix => matrix
    class default
      check(.false.)
  end select

  fe_space => fe_affine_operator%get_fe_space()
  select type (fe_space)
    class is (par_fe_space_t)
      this%par_fe_space => fe_space
    class default
      check(.false.)
  end select

  this%parameter_list => parameter_list

  par_environment => this%par_fe_space%get_par_environment()
  assert (associated(par_environment))

  if (par_environment%am_i_l1_task()) then
    call this%setup_stiffness_weighting()
    call this%setup_object_lid_to_dof_lids()
    call this%setup_dof_lid_to_cdof_id_in_object()
  end if

end subroutine unfitted_l1_create

!========================================================================================
subroutine unfitted_l1_free(this)
  implicit none
  class(unfitted_l1_coarse_fe_handler_t), intent(inout) :: this
  this%matrix         => null()
  this%par_fe_space   => null()
  this%parameter_list => null()
  if ( allocated(this%stiffness_weighting) ) call memfree(this%stiffness_weighting,__FILE__,__LINE__)
  if ( allocated(this%dof_lid_to_cdof_id_in_object) ) call memfree(this%dof_lid_to_cdof_id_in_object,__FILE__,__LINE__)
  if ( allocated(this%object_lid_to_min_neigbour) ) call memfree(this%object_lid_to_min_neigbour,__FILE__,__LINE__)
  call this%object_lid_to_dof_lids%free()
end subroutine unfitted_l1_free

!========================================================================================
subroutine unfitted_l1_get_num_coarse_dofs(this,par_fe_space,parameter_list,num_coarse_dofs)
  implicit none
  class(unfitted_l1_coarse_fe_handler_t), intent(in)    :: this
  type(par_fe_space_t)                  , intent(in)    :: par_fe_space
  type(parameterlist_t)                 , intent(in)    :: parameter_list
  integer(ip)                           , intent(inout) :: num_coarse_dofs(:)

  type(environment_t), pointer  :: par_environment
  integer(ip)                   :: max_cdof_lid_in_object
  logical, allocatable          :: visited_cdof_lids_in_object(:)
  type(list_iterator_t)         :: dofs_in_object_iterator
  integer(ip)                   :: dof_lid
  type(fe_object_iterator_t)    :: object_iterator
  type(fe_object_accessor_t)    :: object

  par_environment => this%par_fe_space%get_par_environment()
  assert ( associated ( par_environment ) )
  assert ( par_environment%am_i_l1_task() )
  assert ( size(num_coarse_dofs) == this%par_fe_space%get_number_fe_objects() )


  max_cdof_lid_in_object = maxval(this%dof_lid_to_cdof_id_in_object) + 1
  call memalloc(max_cdof_lid_in_object,visited_cdof_lids_in_object,__FILE__,__LINE__)

  ! Loop in objects
  object_iterator = par_fe_space%create_fe_object_iterator()
  do while ( .not. object_iterator%has_finished() )
    call object_iterator%current(object)

    ! Count how many coarse dofs are on this object:
    ! i.e., loop on the local dofs of the object and count how many different numbers are found
    visited_cdof_lids_in_object(:) = .false.
    dofs_in_object_iterator = this%object_lid_to_dof_lids%create_iterator(object%get_lid())
    do while (.not. dofs_in_object_iterator%is_upper_bound())
      dof_lid = dofs_in_object_iterator%get_current()
      visited_cdof_lids_in_object(this%dof_lid_to_cdof_id_in_object(dof_lid)+1) = .true.
      call dofs_in_object_iterator%next()
    end do
    num_coarse_dofs(object%get_lid()) =  count(visited_cdof_lids_in_object)

    call object_iterator%next()
  end do

  call memfree(visited_cdof_lids_in_object,__FILE__,__LINE__)

  !call this%standard_l1_coarse_fe_handler_t%get_num_coarse_dofs(par_fe_space,parameter_list,num_coarse_dofs)

end subroutine unfitted_l1_get_num_coarse_dofs

!========================================================================================
subroutine unfitted_l1_setup_constraint_matrix(this,par_fe_space,parameter_list,constraint_matrix)
  implicit none
  class(unfitted_l1_coarse_fe_handler_t), intent(in)    :: this
  type(par_fe_space_t)                  , intent(in)    :: par_fe_space
  type(parameterlist_t)                 , intent(in)    :: parameter_list
  type(coo_sparse_matrix_t)             , intent(inout) :: constraint_matrix

  type(environment_t), pointer :: par_environment
  integer(ip)                  :: field_id, block_id
  integer(ip),         pointer :: field_to_block(:)
  integer(ip)                  :: num_cols, num_rows
  integer(ip)                  :: max_cdof_lid_in_object
  type(list_iterator_t)        :: dofs_in_object_iterator
  integer(ip)                  :: dof_lid, cdof_lid_in_object
  integer(ip)                  :: num_fdofs_in_cdof
  integer(ip)                  :: cdof_lid
  type(fe_object_iterator_t)   :: object_iterator
  type(fe_object_accessor_t)   :: object

  par_environment => this%par_fe_space%get_par_environment()
  assert (associated(par_environment))
  assert (par_environment%am_i_l1_task())

  ! We assume a single field for the moment
  field_id = 1
  assert(this%par_fe_space%get_number_fields() == 1)

  ! Free any dynamic memory that constraint_matrix may have inside
  call constraint_matrix%free()

  ! Create constraint matrix (transposed)
  field_to_block => this%par_fe_space%get_field_blocks()
  block_id = field_to_block(field_id)
  num_cols = this%par_fe_space%get_block_number_dofs(block_id)
  num_rows = this%par_fe_space%get_block_number_coarse_dofs(block_id)
  call constraint_matrix%create ( num_cols, num_rows )

  cdof_lid = 0
  max_cdof_lid_in_object = maxval(this%dof_lid_to_cdof_id_in_object) + 1

  ! Loop in objects
  object_iterator = this%par_fe_space%create_fe_object_iterator()
  do while ( .not. object_iterator%has_finished() )
    call object_iterator%current(object)

    ! Loop in all possible values of cdof_lid_in_object
    ! TODO Maybe these are too many loops
    do cdof_lid_in_object = 0,max_cdof_lid_in_object

      ! Count how many fine dofs are on this c dof
      ! i.e., do a loop in the fine dofs on this object and count how many time the current
      ! coarse dof is found
      num_fdofs_in_cdof = 0
      dofs_in_object_iterator = this%object_lid_to_dof_lids%create_iterator(object%get_lid())
      do while (.not. dofs_in_object_iterator%is_upper_bound())
        dof_lid = dofs_in_object_iterator%get_current()
        if ( cdof_lid_in_object == this%dof_lid_to_cdof_id_in_object(dof_lid) ) num_fdofs_in_cdof = num_fdofs_in_cdof + 1
        call dofs_in_object_iterator%next()
      end do

      ! If there are 0 then cycle
      if (num_fdofs_in_cdof == 0) cycle

      ! Increment in 1 the current row
      cdof_lid = cdof_lid + 1

      ! Loop in all fine dofs of this object
        ! If a fine dof is in the current cdof_lid_in_object add 1/num_fdofs_in_cdof in the constaint matrix
      dofs_in_object_iterator = this%object_lid_to_dof_lids%create_iterator(object%get_lid())
      do while (.not. dofs_in_object_iterator%is_upper_bound())
        dof_lid = dofs_in_object_iterator%get_current()
        if ( cdof_lid_in_object == this%dof_lid_to_cdof_id_in_object(dof_lid) ) then
          call constraint_matrix%insert(dof_lid,cdof_lid, 1.0_rp/real(num_fdofs_in_cdof,rp))
        end if
        call dofs_in_object_iterator%next()
      end do

    end do

    call object_iterator%next()
  end do
  call constraint_matrix%sort_and_compress()

  !call this%standard_l1_coarse_fe_handler_t%setup_constraint_matrix(par_fe_space,parameter_list,constraint_matrix)

end subroutine unfitted_l1_setup_constraint_matrix

!========================================================================================
subroutine unfitted_l1_setup_weighting_operator(this,par_fe_space,parameter_list,weighting_operator)
  implicit none
  class(unfitted_l1_coarse_fe_handler_t), intent(in)    :: this
  type(par_fe_space_t)                  , intent(in)    :: par_fe_space
  type(parameterlist_t)                 , intent(in)    :: parameter_list
  real(rp), allocatable                 , intent(inout) :: weighting_operator(:)

  integer(ip)                            :: field_id, block_id
  integer(ip), pointer                   :: field_to_block(:)
  integer(ip)                            :: num_dofs

  ! Clean up
  if (allocated(weighting_operator) ) then
    call memfree ( weighting_operator, __FILE__, __LINE__ )
  end if

  ! We assume a single field for the moment
  field_id = 1
  assert(this%par_fe_space%get_number_fields() == 1)

  ! Allocate the weighting
  field_to_block => this%par_fe_space%get_field_blocks()
  block_id = field_to_block(field_id)
  num_dofs = this%par_fe_space%get_block_number_dofs(block_id)
  call memalloc(num_dofs,weighting_operator,__FILE__,__LINE__)

  ! Set the weighting with the stored value
  weighting_operator(:) = this%stiffness_weighting(:)

end subroutine unfitted_l1_setup_weighting_operator

!========================================================================================
subroutine unfitted_l1_setup_stiffness_weighting(this)

  implicit none
  class(unfitted_l1_coarse_fe_handler_t), intent(inout) :: this

  type(par_scalar_array_t)             :: par_array
  type(environment_t), pointer         :: p_env
  type(dof_import_t),  pointer         :: dof_import
  type(serial_scalar_array_t), pointer :: serial_array
  real(rp), pointer                    :: assembled_diag(:)
  real(rp), allocatable                :: sub_assembled_diag(:)
  integer(ip)                          :: istat
  integer(ip)                          :: block_id

  ! We assume a single block (for the moment)
  block_id = 1

  ! Get the sub-assembled diagonal
  call this%matrix%extract_diagonal(sub_assembled_diag)

  ! Communicate to compute the fully assembled diagonal
  p_env => this%par_fe_space%get_environment()
  dof_import => this%par_fe_space%get_block_dof_import(block_id)
  call par_array%create_and_allocate(p_env, dof_import)
  serial_array   => par_array%get_serial_scalar_array()
  assembled_diag => serial_array%get_entries()
  assembled_diag(:) = sub_assembled_diag(:)
  call par_array%comm()

  ! Compute the weighting
  call memalloc(size(sub_assembled_diag),this%stiffness_weighting, __FILE__, __LINE__ )
  this%stiffness_weighting(:) = sub_assembled_diag(:)/assembled_diag(:)

  ! Clean up
  deallocate(sub_assembled_diag,stat=istat); check(istat == 0)
  call par_array%free()

end subroutine unfitted_l1_setup_stiffness_weighting

!========================================================================================
subroutine unfitted_l1_setup_object_lid_to_dof_lids(this)

  implicit none
  class(unfitted_l1_coarse_fe_handler_t), intent(inout) :: this

  integer(ip)                        :: num_objects
  integer(ip), allocatable           :: object_lid_to_num_dofs_in_object(:)
  integer(ip)                        :: icell_around
  integer(ip)                        :: ivef
  integer(ip)                        :: idof, dof_lid
  logical                            :: dofs_on_vef
  type(environment_t), pointer       :: par_environment
  type(fe_object_iterator_t)         :: object_iterator
  type(fe_object_accessor_t)         :: object
  type(fe_vefs_on_object_iterator_t) :: vefs_on_object_iterator
  type(fe_vef_accessor_t)            :: vef
  type(fe_accessor_t)                :: fe
  type(list_iterator_t)              :: own_dofs_on_vef_iterator
  integer(ip), pointer               :: elem2dof(:)
  logical                            :: use_vertices, use_edges, use_faces
  logical, allocatable               :: visited_dofs(:)
  integer(ip)                        :: field_id, block_id
  integer(ip), pointer               :: field_to_block(:)
  integer(ip)                        :: num_dofs, num_dofs_in_object
  integer(ip), allocatable           :: dof_lids(:)
  integer(ip)                        :: icount
  type(list_iterator_t)              :: dofs_in_object_iterator

  ! We assume a single field for the moment
  field_id = 1
  assert(this%par_fe_space%get_number_fields() == 1)

  ! Auxiliary
  field_to_block => this%par_fe_space%get_field_blocks()
  block_id = field_to_block(field_id)
  num_dofs  = this%par_fe_space%get_block_number_dofs(block_id)
  num_objects = this%par_fe_space%get_number_fe_objects()
  call memalloc(num_dofs,visited_dofs,__FILE__,__LINE__)
  call memalloc(num_dofs,dof_lids,__FILE__,__LINE__)
  call memalloc(num_objects,object_lid_to_num_dofs_in_object,__FILE__,__LINE__)
  icount = 0

  ! Allocate array of min neighbors and initialize with the biggest integer
  call memalloc(num_objects,this%object_lid_to_min_neigbour,__FILE__,__LINE__)
  this%object_lid_to_min_neigbour(:) = huge(num_objects)

  call this%get_coarse_space_use_vertices_edges_faces(this%parameter_list,& 
                                                      use_vertices, &
                                                      use_edges, &
                                                      use_faces)

  ! Fill the auxiliary data
  object_iterator = this%par_fe_space%create_fe_object_iterator()
  do while ( .not. object_iterator%has_finished() )
    call object_iterator%current(object)

    ! Check if c, ce, or cef, and skip object accordingly
    select case ( object%get_dimension() )
    case (0)
      if (.not. use_vertices) then
        call object_iterator%next(); cycle
      end if
    case (1)
      if (.not. use_edges) then
        call object_iterator%next(); cycle
      end if
    case (2)
      if (.not. use_faces) then
        call object_iterator%next(); cycle
      end if
    end select

    ! Loop in vefs of the object
    vefs_on_object_iterator = object%create_fe_vefs_on_object_iterator()
    visited_dofs(:) = .false.
    num_dofs_in_object = 0
    do while ( .not. vefs_on_object_iterator%has_finished() )
      call vefs_on_object_iterator%current(vef)

      ! Loop in ghost cells around the vef
      do icell_around=1, vef%get_num_cells_around()
        call vef%get_cell_around(icell_around,fe)
        if ( fe%is_ghost() ) then

          call fe%get_field_elem2dof(field_id, elem2dof)
          ivef = fe%find_lpos_vef_lid(vef%get_lid())

          ! Loop in own dofs in the vef as seen from the ghost element
          own_dofs_on_vef_iterator = fe%create_own_dofs_on_vef_iterator(ivef, field_id)
          do while ( .not. own_dofs_on_vef_iterator%is_upper_bound() )
            idof    = own_dofs_on_vef_iterator%get_current()
            dof_lid = elem2dof(idof)
            if ( dof_lid > 0 ) then
              if(.not. visited_dofs(dof_lid)) then

                ! Store the minimum (active) neighbor part
                if ( this%object_lid_to_min_neigbour(object%get_lid()) > fe%get_my_part() ) then
                  this%object_lid_to_min_neigbour(object%get_lid()) = fe%get_my_part()
                end if

                ! Store the dofs on this object
                visited_dofs(dof_lid) = .true.
                icount = icount + 1
                num_dofs_in_object = num_dofs_in_object + 1
                dof_lids(icount) = dof_lid

              end if
            end if
            call own_dofs_on_vef_iterator%next()
          end do

        end if
      end do

      call vefs_on_object_iterator%next()
    end do
    object_lid_to_num_dofs_in_object(object%get_lid()) = num_dofs_in_object
    call object_iterator%next()
  end do

  ! Initialize the list
  call this%object_lid_to_dof_lids%free()
  call this%object_lid_to_dof_lids%create(this%par_fe_space%get_number_fe_objects())
  object_iterator = this%par_fe_space%create_fe_object_iterator()
  do while ( .not. object_iterator%has_finished() )
    call object_iterator%current(object)
    call this%object_lid_to_dof_lids%sum_to_pointer_index(&
      object%get_lid(),object_lid_to_num_dofs_in_object(object%get_lid()))
    call object_iterator%next()
  end do
  call this%object_lid_to_dof_lids%calculate_header()
  call this%object_lid_to_dof_lids%allocate_list_from_pointer()

  ! Fill the list
  icount = 0
  object_iterator = this%par_fe_space%create_fe_object_iterator()
  do while ( .not. object_iterator%has_finished() )
    call object_iterator%current(object)

    dofs_in_object_iterator = this%object_lid_to_dof_lids%create_iterator(object%get_lid())
    do while (.not. dofs_in_object_iterator%is_upper_bound())
      icount = icount + 1
      call dofs_in_object_iterator%set_current(dof_lids(icount))
      call dofs_in_object_iterator%next()
    end do

    call object_iterator%next()
  end do

  ! Clean up
  call memfree(object_lid_to_num_dofs_in_object,__FILE__,__LINE__)
  call memfree(dof_lids,__FILE__,__LINE__)
  call memfree(visited_dofs,__FILE__,__LINE__)

end subroutine unfitted_l1_setup_object_lid_to_dof_lids

!========================================================================================
subroutine unfitted_l1_setup_dof_lid_to_cdof_id_in_object(this)

  implicit none
  class(unfitted_l1_coarse_fe_handler_t), intent(inout) :: this

  logical, allocatable                 :: is_problematic_dof(:)
  type(fe_object_iterator_t)           :: object_iterator
  type(fe_object_accessor_t)           :: object
  type(list_iterator_t)                :: dofs_in_object_iterator
  integer(ip)                          :: dof_lid
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
  assert(this%par_fe_space%get_number_fields() == 1)
  
  ! Get my part id
  p_env => this%par_fe_space%get_par_environment()
  my_part_id = p_env%get_l1_rank() + 1

  ! Allocate and initialize the member variable
  field_to_block => this%par_fe_space%get_field_blocks()
  block_id = field_to_block(field_id)
  num_dofs = this%par_fe_space%get_block_number_dofs(block_id)
  call memalloc(num_dofs,this%dof_lid_to_cdof_id_in_object,__FILE__,__LINE__)
  this%dof_lid_to_cdof_id_in_object(:) = 0

  ! Identify which dofs are problematic
  call memalloc(num_dofs,is_problematic_dof,__FILE__,__LINE__)
  call this%identify_problematic_dofs(is_problematic_dof)

  ! Loop in objects
  object_iterator = this%par_fe_space%create_fe_object_iterator()
  do while ( .not. object_iterator%has_finished() )
    call object_iterator%current(object)

    ! Skip objects that are not edges
    if (object%get_dimension() .ne. 1) then
      call object_iterator%next(); cycle
    end if

    ! Skip if we are not the part with minimum id
    assert( this%object_lid_to_min_neigbour(object%get_lid()) .ne. my_part_id )
    if ( this%object_lid_to_min_neigbour(object%get_lid()) < my_part_id ) then
      call object_iterator%next(); cycle
    end if

    ! Loop in fine dofs on this object
      ! If the dof is problematic, mark it as a new corner
    new_corner_counter = 0
    dofs_in_object_iterator = this%object_lid_to_dof_lids%create_iterator(object%get_lid())
    do while (.not. dofs_in_object_iterator%is_upper_bound())
      dof_lid = dofs_in_object_iterator%get_current()
      if ( is_problematic_dof(dof_lid) ) then
        new_corner_counter = new_corner_counter + 1
        this%dof_lid_to_cdof_id_in_object(dof_lid) = new_corner_counter
      end if
      call dofs_in_object_iterator%next()
    end do

    call object_iterator%next()
  end do

  ! Communicate to make it consistent between sub-domains
  dof_import => this%par_fe_space%get_block_dof_import(block_id)
  call par_array%create_and_allocate(p_env, dof_import)
  serial_array   => par_array%get_serial_scalar_array()
  array_entries => serial_array%get_entries()
  array_entries(:) = real(this%dof_lid_to_cdof_id_in_object(:),kind=rp)
  call par_array%comm()
  this%dof_lid_to_cdof_id_in_object(:) = nint(array_entries,kind=ip)

  call memfree(is_problematic_dof,__FILE__,__LINE__)
  call par_array%free()

end subroutine unfitted_l1_setup_dof_lid_to_cdof_id_in_object

!========================================================================================
subroutine unfitted_l1_identify_problematic_dofs(this,is_problematic_dof)

  implicit none
  class(unfitted_l1_coarse_fe_handler_t), intent(in)    :: this
  logical, allocatable,                   intent(inout) :: is_problematic_dof(:)

  integer(ip), allocatable :: is_cut_cell(:)
  integer(ip) :: field_id
  class(par_fe_space_t), pointer :: par_fe_space
  class(par_unfitted_fe_space_t), pointer :: par_unf_fe_space
  class(base_static_triangulation_t), pointer :: triangulation
  type(unfitted_fe_iterator_t) :: fe_iterator
  type(unfitted_fe_accessor_t) :: fe
  type(unfitted_cell_accessor_t), pointer :: cell
  integer(ip) :: num_total_cells
  type(environment_t), pointer :: par_environment
  type(cell_import_t), pointer :: cell_import
  integer(ip), pointer :: elem2dof(:)

  ! We assume a single field for the moment
  field_id = 1
  assert(this%par_fe_space%get_number_fields() == 1)

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
  fe_iterator = par_unf_fe_space%create_unfitted_fe_iterator()
  do while ( .not. fe_iterator%has_finished() )
    call fe_iterator%current(fe)
    if ( fe%is_ghost() ) then
      call fe_iterator%next(); cycle
    end if
    cell => fe%get_unfitted_cell_accessor()
    if (cell%is_cut()) is_cut_cell(cell%get_lid()) = 1
    call fe_iterator%next()
  end do

  ! Communicate so that the ghost also have this info
  ! TODO this could be avoided if ghost cells have also the coordinates
  ! For the structured triangulation seems to be true, but for the unstructured?
  par_environment => triangulation%get_par_environment()
  cell_import => triangulation%get_cell_import()
  if(par_environment%get_l1_size()>1) &
  call par_environment%l1_neighbours_exchange ( cell_import%get_number_neighbours(), &
                                                cell_import%get_neighbours_ids(),    &
                                                cell_import%get_rcv_ptrs(),          &
                                                cell_import%get_rcv_leids(),         &
                                                cell_import%get_number_neighbours(), &
                                                cell_import%get_neighbours_ids(),    &
                                                cell_import%get_snd_ptrs(),          &
                                                cell_import%get_snd_leids(),         &
                                                is_cut_cell )

  ! Mark all the dofs belonging to cut elements
  is_problematic_dof(:) = .false.
  fe_iterator = par_unf_fe_space%create_unfitted_fe_iterator()
  do while ( .not. fe_iterator%has_finished() )
    call fe_iterator%current(fe)
    call fe%get_field_elem2dof(field_id, elem2dof)
    if (is_cut_cell(fe%get_lid())==1) is_problematic_dof(pack(elem2dof,elem2dof>0)) = .true.
    call fe_iterator%next()
  end do


  ! Clean up
  call memfree(is_cut_cell,__FILE__,__LINE__)

end subroutine unfitted_l1_identify_problematic_dofs

end module unfitted_l1_coarse_fe_handler_names
!***************************************************************************************************
