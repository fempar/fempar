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
module environment_names
  use types_names
  use memor_names
  use stdio_names
  use FPL
  use par_io_names
  use uniform_hex_mesh_generator_names
  use timer_names
  ! Parallel modules
  use execution_context_names
  use mpi_context_names
  use mpi_omp_context_names
  use serial_context_names
  implicit none

# include "debug.i90"
  private

  integer(ip) , parameter :: not_created = 0            ! when contexts are not created
  integer(ip) , parameter :: created     = 1            ! when contexts have been created from lower levels
  integer(ip) , parameter :: created_from_scratch = 2   ! when contexts have been created from scratch

  integer(ip)     , parameter :: mpi_context    = 0
  integer(ip)     , parameter :: serial_context = 1
  integer(ip)     , parameter :: mpi_omp_context = 2
  character(len=*), parameter :: execution_context_key  = 'execution_context'
  public :: mpi_context
  public :: serial_context
  public :: mpi_omp_context
  public :: execution_context_key

  integer(ip)     , parameter :: structured   = 0
  integer(ip)     , parameter :: unstructured = 1
  character(len=*), parameter :: environment_type_key = 'environment_type'
  public :: structured
  public :: unstructured
  public :: environment_type_key

  ! This type manages the assigment of tasks to different levels as well as
  ! the dutties assigned to each of them. The array parts_mapping gives the
  ! part assigned to this task and the parts of coarser levels it belongs to. 
  ! Different levels in the hierarchy are managed recursively.
  type ::  environment_t
     private 
     integer(ip)                       :: state = not_created
     integer(ip)                       :: num_levels = 0         ! =0 when parts are not assigned to tasks
     integer(ip)         , allocatable :: num_parts_per_level(:)
     integer(ip)         , allocatable :: parts_mapping(:)
     class(execution_context_t), allocatable :: world_context          ! All tasks context (=l1+lgt1 tasks)
     class(execution_context_t), allocatable :: l1_context             ! 1st lev tasks context
     class(execution_context_t), allocatable :: lgt1_context           ! > 1st lev tasks context
     class(execution_context_t), allocatable :: l1_to_l2_context       ! Mixed l1 and l2 tasks context 
     ! to exchange data between levels
     type(environment_t), pointer :: next_level 
   contains
     procedure :: read                           => environment_read_file
     procedure :: write                          => environment_write_file

     procedure          :: create                => environment_create
     procedure          :: assign_parts_to_tasks => environment_assign_parts_to_tasks
     procedure, private :: fill_contexts         => environment_fill_contexts

     procedure :: free                           => environment_free
     procedure :: print                          => environment_print
     procedure :: created                        => environment_created
     ! Getters
     procedure :: get_num_tasks                  => environment_get_num_tasks
     procedure :: get_next_level                 => environment_get_next_level
     procedure :: get_l1_context                 => environment_get_l1_context
     procedure :: get_l1_rank                    => environment_get_l1_rank
     procedure :: get_l1_size                    => environment_get_l1_size
     procedure :: get_l1_to_l2_rank              => environment_get_l1_to_l2_rank
     procedure :: get_l1_to_l2_size              => environment_get_l1_to_l2_size
     procedure :: get_lgt1_rank                  => environment_get_lgt1_rank
     ! Who am I?
     procedure :: am_i_lgt1_task                 => environment_am_i_lgt1_task
     procedure :: am_i_l1_to_l2_task             => environment_am_i_l1_to_l2_task
     procedure :: am_i_l1_to_l2_root             => environment_am_i_l1_to_l2_root
     procedure :: am_i_l1_root                   => environment_am_i_l1_root

     procedure :: create_l1_timer                => environment_create_l1_timer
     
     procedure :: get_l2_part_id_l1_task_is_mapped_to => environment_get_l2_part_id_l1_task_is_mapped_to


     procedure, private :: environment_l1_neighbours_exchange_rp
     procedure, private :: environment_l1_neighbours_exchange_ip
     procedure, private :: environment_l1_neighbours_exchange_igp
     procedure, private :: environment_l1_neighbours_exchange_single_ip
     procedure, private :: environment_l1_neighbours_exchange_wo_pack_unpack_ieep
     generic   :: l1_neighbours_exchange      => environment_l1_neighbours_exchange_rp, &
          environment_l1_neighbours_exchange_ip,&
          environment_l1_neighbours_exchange_igp,&
          environment_l1_neighbours_exchange_single_ip, &
          environment_l1_neighbours_exchange_wo_pack_unpack_ieep

     procedure, private :: environment_l1_scatter_scalar_ip
     generic   :: l1_scatter => environment_l1_scatter_scalar_ip                                                 

     procedure, private :: environment_l1_gather_scalar_ip
     generic   :: l1_gather => environment_l1_gather_scalar_ip 

     procedure, private :: environment_l1_bcast_scalar_ip
     generic   :: l1_bcast => environment_l1_bcast_scalar_ip 

     procedure, private :: environment_l2_from_l1_gather_ip
     procedure, private :: environment_l2_from_l1_gather_igp
     procedure, private :: environment_l2_from_l1_gather_ip_1D_array
     procedure, private :: environment_l2_from_l1_gatherv_ip_1D_array
     procedure, private :: environment_l2_from_l1_gatherv_igp_1D_array
     procedure, private :: environment_l2_from_l1_gatherv_rp_1D_array
     procedure, private :: environment_l2_from_l1_gatherv_rp_2D_array
     generic   :: l2_from_l1_gather => environment_l2_from_l1_gather_ip, &
          environment_l2_from_l1_gather_igp, &
          environment_l2_from_l1_gather_ip_1D_array, &
          environment_l2_from_l1_gatherv_ip_1D_array, &
          environment_l2_from_l1_gatherv_igp_1D_array, &
          environment_l2_from_l1_gatherv_rp_1D_array, &
          environment_l2_from_l1_gatherv_rp_2D_array

     procedure, private :: environment_l2_to_l1_scatterv_rp_1D_array                                  
     generic   :: l2_to_l1_scatter => environment_l2_to_l1_scatterv_rp_1D_array


     procedure, private :: environment_l1_to_l2_transfer_ip
     procedure, private :: environment_l1_to_l2_transfer_ip_1D_array
     generic  :: l1_to_l2_transfer => environment_l1_to_l2_transfer_ip, &
          environment_l1_to_l2_transfer_ip_1D_array             


     ! Deferred TBPs inherited from class(environment_t)
     !procedure :: info                        => environment_info
     procedure :: am_i_l1_task                => environment_am_i_l1_task
     procedure :: l1_lgt1_bcast               => environment_l1_lgt1_bcast
     procedure :: l1_barrier                  => environment_l1_barrier
     procedure :: l1_sum_scalar_rp            => environment_l1_sum_scalar_rp
     procedure :: l1_sum_vector_rp            => environment_l1_sum_vector_rp
     procedure :: l1_max_scalar_rp            => environment_l1_max_scalar_rp
     procedure :: l1_max_vector_rp            => environment_l1_max_vector_rp
     generic  :: l1_sum                       => l1_sum_scalar_rp, l1_sum_vector_rp
     generic  :: l1_max                       => l1_max_scalar_rp, l1_max_vector_rp

  end type environment_t

  ! Types
  public :: environment_t, environment_compose_name, environment_write_files

contains

  !=============================================================================
  subroutine environment_compose_name ( prefix, name ) 
    implicit none
    character (len=*), intent(in)    :: prefix 
    character (len=:), allocatable, intent(inout) :: name
    name = trim(prefix) // '.env'
  end subroutine environment_compose_name

  !=============================================================================
  subroutine environment_write_files ( parameter_list, envs )
    implicit none
    ! Parameters
    type(ParameterList_t)  , intent(in) :: parameter_list
    type(environment_t), intent(in)  :: envs(:)

    ! Locals
    integer(ip)          :: nenvs
    integer(ip)          :: istat
    logical              :: is_present
    character(len=:), allocatable :: dir_path
    character(len=:), allocatable :: prefix
    character(len=:), allocatable :: name, rename
    integer(ip)          :: lunio
    integer(ip)          :: i
    logical              :: same_data_type
    integer(ip), allocatable :: shape(:)

    nenvs = size(envs)

     ! Mandatory parameters
    is_present         = parameter_list%isPresent(Key=dir_path_key)
    if(is_present) then
#ifdef DEBUG
        same_data_type = parameter_list%isOfDataType(Key = dir_path_key, mold = dir_path)
        istat          = parameter_list%getshape(Key=dir_path_key, shape=shape)
        assert(istat == 0)
        if(same_data_type .and. size(shape) == 0) then
#endif
            istat = parameter_list%GetAsString(key = dir_path_key, string = dir_path)
            check(istat==0)
#ifdef DEBUG
        else
            write(*,'(a)') ' Warning! '//trim(dir_path_key)//' ignored. Wrong data type or shape. '
        endif
#endif
    endif


    is_present         = parameter_list%isPresent(Key=prefix_key)
    if(is_present) then
#ifdef DEBUG
        same_data_type = parameter_list%isOfDataType(Key = prefix_key, mold = prefix)
        istat          = parameter_list%getshape(Key=prefix_key, shape=shape)
        assert(istat == 0)
        if(same_data_type .and. size(shape) == 0) then
#endif
            istat = parameter_list%GetAsString(key = prefix_key, string = prefix)
            check(istat==0)
#ifdef DEBUG
        else
            write(*,'(a)') ' Warning! '//trim(dir_path_key)//' ignored. Wrong data type or shape. '
        endif
#endif
    endif

    call environment_compose_name ( prefix, name )

    do i=1,nenvs
       rename=name
       call numbered_filename_compose(i,nenvs,rename)
       lunio = io_open (trim(dir_path) // '/' // trim(rename)); check(lunio>0)
       call envs(i)%write (lunio)
       call io_close (lunio)
    end do

    ! name, and rename should be automatically deallocated by the compiler when they
    ! go out of scope. Should we deallocate them explicitly for safety reasons?
  end subroutine  environment_write_files

  !=============================================================================
  subroutine environment_read_file ( this,lunio)
    implicit none 
    class(environment_t), intent(inout) :: this
    integer(ip)             , intent(in)    :: lunio
    call this%free()
    read ( lunio, '(10i10)' ) this%num_levels
    call memalloc ( this%num_levels, this%num_parts_per_level,__FILE__,__LINE__  )
    call memalloc ( this%num_levels, this%parts_mapping,__FILE__,__LINE__  )
    read ( lunio, '(10i10)' ) this%num_parts_per_level
    read ( lunio, '(10i10)' ) this%parts_mapping
  end subroutine environment_read_file

  !=============================================================================
  subroutine environment_write_file ( this,lunio)
    implicit none 
    class(environment_t), intent(in) :: this
    integer(ip)             , intent(in) :: lunio
    assert( this%num_levels>0)
    write ( lunio, '(10i10)' ) this%num_levels
    write ( lunio, '(10i10)' ) this%num_parts_per_level
    write ( lunio, '(10i10)' ) this%parts_mapping    
  end subroutine environment_write_file

  !=============================================================================
  subroutine environment_create ( this, parameters)
    !$ use omp_lib
    implicit none 
    class(environment_t), intent(inout) :: this
    type(ParameterList_t)   , intent(in)    :: parameters

    ! Some refactoring is needed here separating number_of_parts_per_level (and dir)
    ! from the rest of the mesh information.
    type(uniform_hex_mesh_t) :: uniform_hex_mesh

    integer(ip)          :: istat
    logical              :: is_present
    integer(ip)          :: environment_type
    integer(ip)          :: execution_context
    character(len=256)   :: dir_path
    character(len=256)   :: prefix
    character(len=:), allocatable   :: name
    integer(ip)                     :: lunio
    integer(ip)               :: num_levels
    integer(ip) , allocatable :: num_parts_per_level(:)
    integer(ip) , allocatable :: parts_mapping(:)

    call this%free()

    ! Optional parameters
    if(parameters%isPresent(key = execution_context_key)) then
       istat = parameters%get(key = execution_context_key, value = execution_context); check(istat==0)
    else
       execution_context = serial_context
    end if
    if( parameters%isPresent(key = environment_type_key) ) then
       istat = parameters%get(key = environment_type_key, value = environment_type); check(istat==0)
    else
       environment_type = unstructured
    end if

    if(execution_context==serial_context) then
       allocate(serial_context_t :: this%world_context,stat=istat); check(istat==0)
    else if(execution_context==mpi_context) then
       allocate(mpi_context_t :: this%world_context,stat=istat); check(istat==0)
    else if(execution_context==mpi_omp_context) then
       allocate(mpi_omp_context_t :: this%world_context,stat=istat); check(istat==0)
    end if
    call this%world_context%create()

    if(environment_type==unstructured) then

       if(this%world_context%get_num_tasks()>1) then
          ! Mandatory parameters
          is_present =  parameters%isPresent(key = dir_path_key)      ; assert(is_present)
          is_present =  parameters%isPresent(key = prefix_key)        ; assert(is_present)
          istat = parameters%get(key = dir_path_key, value = dir_path); check(istat==0)
          istat = parameters%get(key = prefix_key  , value = prefix)  ; check(istat==0)

          call environment_compose_name(prefix, name )  
          call par_filename( this%world_context%get_current_task()+1, this%world_context%get_num_tasks() , name )
          lunio = io_open( trim(dir_path) // '/' // trim(name), 'read' ); check(lunio>0)
          !$ write(*,*) omp_get_thread_num(),lunio 

          ! Read parts assignment to tasks and verify that the multilevel environment 
          ! can be executed in current context (long-lasting Alberto's concern). 
          call this%read(lunio)
          check(this%get_num_tasks() <= this%world_context%get_num_tasks())

          ! Recursively create multilevel environment
          call this%fill_contexts()

       else
          allocate(this%l1_context,mold=this%world_context,stat=istat);check(istat==0)
          this%l1_context = this%world_context
          allocate(this%lgt1_context,mold=this%world_context,stat=istat);check(istat==0)
          call this%lgt1_context%nullify()
          allocate(this%l1_to_l2_context,mold=this%world_context,stat=istat);check(istat==0)
          call this%l1_to_l2_context%nullify()
       end if

    else if(environment_type==structured) then

       call uniform_hex_mesh%get_data_from_parameter_list(parameters)

       ! Generate parts, assign them to tasks and verify that the multilevel environment 
       ! can be executed in current context (long-lasting Alberto's concern). 
       call uniform_hex_mesh%generate_levels_and_parts(this%world_context%get_current_task(), &
            &                                          num_levels, &
            &                                          num_parts_per_level, &
            &                                          parts_mapping)
       call this%assign_parts_to_tasks(num_levels, num_parts_per_level, parts_mapping)
       call memfree(num_parts_per_level,__FILE__,__LINE__)
       call memfree(parts_mapping,__FILE__,__LINE__)
       check(this%get_num_tasks() <= this%world_context%get_num_tasks())

       ! Recursively create multilevel environment
       call this%fill_contexts()

       call uniform_hex_mesh%free()

    end if
    this%state = created_from_scratch

  end subroutine environment_create

  !=============================================================================
  recursive subroutine environment_fill_contexts (this)
    implicit none 
    ! Parameters
    class(environment_t), intent(inout) :: this
    integer                             :: my_color
    integer(ip)                         :: istat

    assert ( this%num_levels >= 1 )
    assert ( allocated(this%world_context))
    assert ( this%world_context%get_current_task() >= 0 )

    ! Create this%l1_context and this%lgt1_context by splitting world_context
    call this%world_context%split_by_condition ( this%world_context%get_current_task() < this%num_parts_per_level(1), this%l1_context, this%lgt1_context )

    ! Create l1_to_l2_context, where inter-level data transfers actually occur
    if ( this%num_levels > 1 ) then
       if(this%l1_context%get_current_task() >= 0) then
          my_color = this%parts_mapping(2)
       else if( this%lgt1_context%get_current_task() < this%num_parts_per_level(2)  ) then
          my_color = this%lgt1_context%get_current_task()+1
       else
          my_color = undefined_color
       end if
       call this%world_context%split_by_color ( my_color, this%l1_to_l2_context )
    else
       allocate(this%l1_to_l2_context,mold=this%world_context,stat=istat);check(istat==0)
       call this%l1_to_l2_context%nullify()
    end if

    if ( this%num_levels > 1 .and. this%lgt1_context%get_current_task() >= 0 ) then
       allocate(this%next_level, stat=istat);check(istat == 0)
       allocate(this%next_level%world_context,mold=this%world_context,stat=istat);check(istat==0)
       this%next_level%world_context = this%lgt1_context
       call this%next_level%assign_parts_to_tasks(this%num_levels-1, this%num_parts_per_level(2:),this%parts_mapping(2:))
       call this%next_level%fill_contexts()
    else
       nullify(this%next_level)
    end if
    this%state = created

  end subroutine environment_fill_contexts

  !=============================================================================
  subroutine environment_assign_parts_to_tasks ( this, num_levels, num_parts_per_level, parts_mapping)
    implicit none 
    ! Parameters
    class(environment_t), intent(inout) :: this
    integer(ip)         , intent(in)    :: num_levels
    integer(ip)         , intent(in)    :: num_parts_per_level(num_levels)
    integer(ip)         , intent(in)    :: parts_mapping(num_levels)
    integer(ip)                         :: istat

    assert ( num_levels >= 1 )
    call this%free()
    this%num_levels = num_levels
    call memalloc(this%num_levels, this%parts_mapping,__FILE__,__LINE__ )
    call memalloc(this%num_levels, this%num_parts_per_level,__FILE__,__LINE__ )
    this%parts_mapping = parts_mapping
    this%num_parts_per_level = num_parts_per_level

  end subroutine environment_assign_parts_to_tasks

  !=============================================================================
  recursive subroutine environment_free ( this )
    implicit none 
    ! Parameters
    class(environment_t), intent(inout) :: this
    integer(ip)                             :: istat

    if(this%state >= created) then
       if ( this%num_levels > 1 .and. this%lgt1_context%get_current_task() >= 0 ) then
          call this%next_level%free()
          deallocate ( this%next_level, stat = istat ); assert ( istat == 0 )
       end if
       call this%l1_context%free(finalize=.false.)
       call this%lgt1_context%free(finalize=.false.)
       call this%l1_to_l2_context%free(finalize=.false.)
       call this%world_context%free(finalize=(this%state == created_from_scratch))
       this%state = not_created
    end if
    if(this%num_levels > 0) then
       this%num_levels = 0
       call memfree(this%parts_mapping , __FILE__, __LINE__ )
       call memfree(this%num_parts_per_level, __FILE__, __LINE__ )
    end if

  end subroutine environment_free

  !=============================================================================
  recursive subroutine environment_print ( this )
    implicit none 
    ! Parameters
    class(environment_t), intent(in) :: this
    integer(ip)                          :: istat

    if (this%state>=created) then
       write(*,*) 'LEVELS: ', this%num_levels, 'l1_context      : ',this%l1_context%get_current_task(), this%l1_context%get_num_tasks()
       write(*,*) 'LEVELS: ', this%num_levels, 'lgt1_context    : ',this%lgt1_context%get_current_task(), this%lgt1_context%get_num_tasks()
       !write(*,*) 'LEVELS: ', this%num_levels, 'l1_lgt1_context : ',this%l1_lgt1_context%get_rank(), this%l1_lgt1_context%get_size()
       write(*,*) 'LEVELS: ', this%num_levels, 'l1_to_l2_context: ',this%l1_to_l2_context%get_current_task(), this%l1_to_l2_context%get_num_tasks()
       if ( this%num_levels > 1 .and. this%lgt1_context%get_current_task() >= 0 ) then
          call this%next_level%print()
       end if
    end if
  end subroutine environment_print

  !=============================================================================
  function environment_created ( this )
    implicit none 
    class(environment_t), intent(in) :: this
    logical                              :: environment_created
    environment_created =  (this%state>=created)
  end function environment_created

  !=============================================================================

  function environment_get_num_tasks (this)
    implicit none 
    class(environment_t), intent(in) :: this
    integer(ip) :: ilevel, environment_get_num_tasks
    assert(this%num_levels>0)
    environment_get_num_tasks = 0
    do ilevel=1,this%num_levels
       environment_get_num_tasks = environment_get_num_tasks + this%num_parts_per_level(ilevel)
    end do
  end function environment_get_num_tasks

  !=============================================================================
  function environment_get_next_level( this )
    implicit none 
    ! Parameters
    class(environment_t), target, intent(in) :: this
    type(environment_t)         , pointer    :: environment_get_next_level

    environment_get_next_level => this%next_level
  end function environment_get_next_level
  
  !=============================================================================
  function environment_get_l1_context ( this ) result(l1_context)
    implicit none 
    ! Parameters
    class(environment_t),       target, intent(in) :: this
    class(execution_context_t), pointer            :: l1_context
    l1_context => this%l1_context
  end function environment_get_l1_context

  !=============================================================================
  function environment_get_l1_rank ( this )
    implicit none 
    ! Parameters
    class(environment_t), intent(in) :: this
    integer                              :: environment_get_l1_rank
    environment_get_l1_rank = this%l1_context%get_current_task()
  end function environment_get_l1_rank

  !=============================================================================
  function environment_get_l1_size ( this )
    implicit none 
    ! Parameters
    class(environment_t), intent(in) :: this
    integer                              :: environment_get_l1_size
    environment_get_l1_size = this%l1_context%get_num_tasks()
  end function environment_get_l1_size

  !=============================================================================
  function environment_get_l1_to_l2_rank ( this )
    implicit none 
    ! Parameters
    class(environment_t), intent(in) :: this
    integer                              :: environment_get_l1_to_l2_rank
    environment_get_l1_to_l2_rank = this%l1_to_l2_context%get_current_task()
  end function environment_get_l1_to_l2_rank

  !=============================================================================
  pure function environment_get_l1_to_l2_size ( this )
    implicit none 
    ! Parameters
    class(environment_t), intent(in) :: this
    integer                              :: environment_get_l1_to_l2_size
    environment_get_l1_to_l2_size = this%l1_to_l2_context%get_num_tasks()
  end function environment_get_l1_to_l2_size

  !=============================================================================
  function environment_get_lgt1_rank ( this )
    implicit none 
    ! Parameters
    class(environment_t), intent(in) :: this
    integer                              :: environment_get_lgt1_rank
    environment_get_lgt1_rank = this%lgt1_context%get_current_task()
  end function environment_get_lgt1_rank

  !=============================================================================
  function environment_am_i_lgt1_task(this) 
    implicit none
    class(environment_t) ,intent(in)  :: this
    logical                               :: environment_am_i_lgt1_task
    environment_am_i_lgt1_task = (this%lgt1_context%get_current_task() >= 0)
  end function environment_am_i_lgt1_task

  !=============================================================================
  function environment_am_i_l1_to_l2_task(this)
    implicit none
    class(environment_t), intent(in) :: this
    logical                              :: environment_am_i_l1_to_l2_task 
    environment_am_i_l1_to_l2_task = (this%l1_to_l2_context%get_current_task() >= 0)
  end function environment_am_i_l1_to_l2_task

  !=============================================================================
  function environment_am_i_l1_to_l2_root(this)
    implicit none
    class(environment_t), intent(in) :: this
    logical                              :: environment_am_i_l1_to_l2_root
    environment_am_i_l1_to_l2_root = (this%l1_to_l2_context%get_current_task() == this%l1_to_l2_context%get_num_tasks()-1)
  end function environment_am_i_l1_to_l2_root

  !=============================================================================
  function environment_am_i_l1_root(this)
    implicit none
    class(environment_t), intent(in) :: this
    logical                              :: environment_am_i_l1_root
    !environment_am_i_l1_root = (this%l1_context%get_current_task() == 0)
    environment_am_i_l1_root = this%l1_context%am_i_root()
  end function environment_am_i_l1_root

  !=============================================================================
  function environment_get_l2_part_id_l1_task_is_mapped_to (this)
    implicit none
    class(environment_t), intent(in) :: this
    integer(ip) :: environment_get_l2_part_id_l1_task_is_mapped_to 
    environment_get_l2_part_id_l1_task_is_mapped_to = this%parts_mapping(2) 
  end function environment_get_l2_part_id_l1_task_is_mapped_to

  !=============================================================================
  subroutine environment_l1_barrier(this) 
    implicit none
    class(environment_t),intent(in)  :: this
    call this%l1_context%barrier()
  end subroutine environment_l1_barrier

  !=============================================================================
  subroutine environment_info(this,me,np) 
    implicit none
    class(environment_t),intent(in)  :: this
    integer(ip)           ,intent(out) :: me
    integer(ip)           ,intent(out) :: np
    me = this%l1_context%get_current_task()
    np = this%l1_context%get_num_tasks()
  end subroutine environment_info

  !=============================================================================
  function environment_am_i_l1_task(this) 
    implicit none
    class(environment_t) ,intent(in)  :: this
    logical                             :: environment_am_i_l1_task 
    environment_am_i_l1_task = (this%l1_context%get_current_task() >= 0)
  end function environment_am_i_l1_task

  !=============================================================================
  subroutine environment_l1_lgt1_bcast(this,condition)
    implicit none 
    ! Parameters
    class(environment_t), intent(in)    :: this
    logical                 , intent(inout) :: condition
    assert ( this%created() )
    call this%world_context%bcast_subcontext(this%l1_context,this%lgt1_context,condition)
  end subroutine environment_l1_lgt1_bcast

  !=============================================================================
  subroutine environment_create_l1_timer(this, message, mode, timer)
    class(environment_t)      , intent(in)    :: this
    character(len=*)          , intent(in)    :: message
    character(len=*), optional, intent(in)    :: mode
    type(timer_t)             , intent(out)   :: timer
    call timer%create(this%l1_context, message, mode)
  end subroutine environment_create_l1_timer  
  
  !=============================================================================
  subroutine environment_l1_sum_scalar_rp (this,alpha)
    implicit none
    class(environment_t) , intent(in)    :: this
    real(rp)                 , intent(inout) :: alpha
    assert ( this%am_i_l1_task() )
    call this%l1_context%sum_scalar_rp(alpha)
  end subroutine environment_l1_sum_scalar_rp

  !=============================================================================
  subroutine environment_l1_sum_vector_rp(this,alpha)
    implicit none
    class(environment_t) , intent(in)    :: this
    real(rp)                 , intent(inout) :: alpha(:) 
    assert (this%am_i_l1_task())
    call this%l1_context%sum_vector_rp(alpha)
  end subroutine environment_l1_sum_vector_rp

  !=============================================================================
  subroutine environment_l1_max_scalar_rp (this,alpha)
    implicit none
    class(environment_t) , intent(in)    :: this
    real(rp)                 , intent(inout) :: alpha
    assert ( this%am_i_l1_task() )
    call this%l1_context%max_scalar_rp(alpha)
  end subroutine environment_l1_max_scalar_rp

  !=============================================================================
  subroutine environment_l1_max_vector_rp(this,alpha)
    implicit none
    class(environment_t) , intent(in)    :: this
    real(rp)                 , intent(inout) :: alpha(:) 
    assert (this%am_i_l1_task())
    call this%l1_context%max_vector_rp(alpha)
  end subroutine environment_l1_max_vector_rp

  !=============================================================================
  !=============================================================================
  subroutine environment_l1_neighbours_exchange_rp ( this, & 
       num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       num_snd, list_snd, snd_ptrs, pack_idx,   &
       alpha, beta, x)
    implicit none
    class(environment_t), intent(in) :: this

    ! Control info to receive
    integer(ip)             , intent(in) :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(in) :: unpack_idx (rcv_ptrs(num_rcv+1)-1)

    ! Control info to send
    integer(ip)             , intent(in) :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in) :: pack_idx (snd_ptrs(num_snd+1)-1)

    ! Floating point data
    real(rp), intent(in)    :: alpha, beta
    real(rp), intent(inout) :: x(:)

    assert (this%am_i_l1_task())
    call this%l1_context%neighbours_exchange ( num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
         &                                      num_snd, list_snd, snd_ptrs, pack_idx,   &
         &                                      alpha, beta, x)

  end subroutine environment_l1_neighbours_exchange_rp
  !=============================================================================
  subroutine environment_l1_neighbours_exchange_ip ( this, & 
       &                                                 num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       &                                                 num_snd, list_snd, snd_ptrs, pack_idx,   &
       &                                                 x, chunk_size)
    implicit none
    class(environment_t), intent(in)    :: this
    ! Control info to receive
    integer(ip)             , intent(in)    :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(in)    :: unpack_idx (rcv_ptrs(num_rcv+1)-1)
    ! Control info to send
    integer(ip)             , intent(in)    :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in)    :: pack_idx (snd_ptrs(num_snd+1)-1)
    ! Raw data to be exchanged
    integer(ip)             , intent(inout) :: x(:)
    integer(ip)   , optional, intent(in)    :: chunk_size

    assert( this%am_i_l1_task() )
    call this%l1_context%neighbours_exchange ( num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
         &                                     num_snd, list_snd, snd_ptrs, pack_idx,   &
         &                                     x,chunk_size)

  end subroutine environment_l1_neighbours_exchange_ip

  !=============================================================================
  subroutine environment_l1_neighbours_exchange_igp ( this, & 
       num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       num_snd, list_snd, snd_ptrs, pack_idx,   &
       x, chunk_size)
    implicit none
    class(environment_t), intent(in)    :: this
    ! Control info to receive
    integer(ip)             , intent(in)    :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(in)    :: unpack_idx (rcv_ptrs(num_rcv+1)-1)
    ! Control info to send
    integer(ip)             , intent(in)    :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in)    :: pack_idx (snd_ptrs(num_snd+1)-1)
    ! Raw data to be exchanged
    integer(igp)            , intent(inout) :: x(:)
    integer(ip)   , optional, intent(in)    :: chunk_size

    assert( this%am_i_l1_task() )

    call this%l1_context%neighbours_exchange ( num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
         &                                     num_snd, list_snd, snd_ptrs, pack_idx,   &
         &                                     x, chunk_size)

  end subroutine environment_l1_neighbours_exchange_igp

  !=============================================================================
  subroutine environment_l1_neighbours_exchange_single_ip ( this, & 
       num_neighbours, &
       list_neighbours, &
       input_data,&
       output_data)
    implicit none
    class(environment_t), intent(in)    :: this

    integer                 , intent(in)    :: num_neighbours
    integer(ip)             , intent(in)    :: list_neighbours (num_neighbours)
    integer(ip)             , intent(in)    :: input_data
    integer(ip)             , intent(inout) :: output_data(num_neighbours)

    assert  (this%am_i_l1_task())
    call this%l1_context%neighbours_exchange ( num_neighbours, &
         &                                     list_neighbours, &
         &                                     input_data,&
         &                                     output_data)

  end subroutine environment_l1_neighbours_exchange_single_ip

  !=============================================================================
  subroutine environment_l1_neighbours_exchange_wo_pack_unpack_ieep ( this, &
       &                                                                  number_neighbours, &
       &                                                                  neighbour_ids, &
       &                                                                  snd_ptrs, &
       &                                                                  snd_buf, & 
       &                                                                  rcv_ptrs, &
       &                                                                  rcv_buf )
    ! Parameters
    class(environment_t)  , intent(in)    :: this 
    integer(ip)               , intent(in)    :: number_neighbours
    integer(ip)               , intent(in)    :: neighbour_ids(number_neighbours)
    integer(ip)               , intent(in)    :: snd_ptrs(number_neighbours+1)
    integer(ieep)             , intent(in)    :: snd_buf(snd_ptrs(number_neighbours+1)-1)   
    integer(ip)               , intent(in)    :: rcv_ptrs(number_neighbours+1)
    integer(ieep)             , intent(out)   :: rcv_buf(rcv_ptrs(number_neighbours+1)-1)

    assert ( this%am_i_l1_task() )
    call this%l1_context%neighbours_exchange ( number_neighbours, &
         &                                     neighbour_ids, &
         &                                     snd_ptrs, &
         &                                     snd_buf, & 
         &                                     rcv_ptrs, &
         &                                     rcv_buf )

  end subroutine environment_l1_neighbours_exchange_wo_pack_unpack_ieep

  !=============================================================================
  subroutine environment_l1_gather_scalar_ip ( this, input_data, output_data )
    implicit none
    class(environment_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data
    integer(ip)             , intent(out)  :: output_data(:) ! ( this%l1_context%get_num_tasks())
    assert ( this%am_i_l1_task() )
    call this%l1_context%gather ( input_data, output_data )
  end subroutine environment_l1_gather_scalar_ip

  subroutine environment_l1_scatter_scalar_ip ( this, input_data, output_data )
    implicit none
    class(environment_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data(:) ! ( this%l1_context%get_num_tasks())
    integer(ip)             , intent(out)  :: output_data
    assert( this%am_i_l1_task() )
    call this%l1_context%scatter ( input_data, output_data )
  end subroutine environment_l1_scatter_scalar_ip

  subroutine environment_l1_bcast_scalar_ip ( this, data )
    implicit none
    class(environment_t), intent(in)    :: this
    integer(ip)             , intent(inout) :: data
    integer(ip) :: icontxt
    assert ( this%am_i_l1_task() )
    call this%l1_context%bcast ( data )
  end subroutine environment_l1_bcast_scalar_ip

  !=============================================================================
  !=============================================================================
  subroutine environment_l1_to_l2_transfer_ip ( this, input_data, output_data )
    implicit none
    class(environment_t), intent(in)      :: this
    integer(ip)             , intent(in)      :: input_data
    integer(ip)             , intent(inout)   :: output_data
    assert ( this%am_i_l1_to_l2_task() )
    call this%l1_to_l2_context%root_send_master_rcv(input_data, output_data )
  end subroutine environment_l1_to_l2_transfer_ip

  !=============================================================================
  subroutine environment_l1_to_l2_transfer_ip_1D_array ( this, input_data, output_data )
    implicit none
    class(environment_t), intent(in)      :: this
    integer(ip)             , intent(in)      :: input_data(:)
    integer(ip)             , intent(inout)   :: output_data(:)
    assert ( this%am_i_l1_to_l2_task() )
    call this%l1_to_l2_context%root_send_master_rcv(input_data, output_data )
  end subroutine environment_l1_to_l2_transfer_ip_1D_array

  !=============================================================================
  subroutine environment_l2_from_l1_gather_ip ( this, input_data, output_data )
    implicit none
    class(environment_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data
    integer(ip)             , intent(out)  :: output_data(:) ! ( this%l1_to_l2_context%get_num_tasks())
    assert ( this%am_i_l1_to_l2_task() )
    call this%l1_to_l2_context%gather_to_master(input_data,output_data )
  end subroutine environment_l2_from_l1_gather_ip

  !=============================================================================
  subroutine environment_l2_from_l1_gather_igp ( this, input_data, output_data )
    implicit none
    class(environment_t), intent(in)   :: this
    integer(igp)            , intent(in)   :: input_data
    integer(igp)            , intent(out)  :: output_data(:) ! ( this%l1_to_l2_context%get_num_tasks())
    assert ( this%am_i_l1_to_l2_task() )
    call this%l1_to_l2_context%gather_to_master(input_data,output_data )
  end subroutine environment_l2_from_l1_gather_igp

  !=============================================================================
  subroutine environment_l2_from_l1_gather_ip_1D_array ( this, input_data_size, input_data, output_data )
    implicit none
    class(environment_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data_size
    integer(ip)             , intent(in)   :: input_data(input_data_size)
    integer(ip)             , intent(out)  :: output_data(:)
    assert ( this%am_i_l1_to_l2_task() )
    call this%l1_to_l2_context%gather_to_master(input_data_size, input_data,output_data )
  end subroutine environment_l2_from_l1_gather_ip_1D_array

  !=============================================================================
  subroutine environment_l2_from_l1_gatherv_ip_1D_array ( this, input_data_size, input_data, recv_counts, displs, output_data )
    implicit none
    class(environment_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data_size
    integer(ip)             , intent(in)   :: input_data(input_data_size)
    integer(ip)             , intent(in)   :: recv_counts(:) ! ( this%l1_to_l2_context%get_num_tasks())
    integer(ip)             , intent(in)   :: displs(:) ! ( this%l1_to_l2_context%get_num_tasks())
    integer(ip)             , intent(out)  :: output_data(:)
    assert ( this%am_i_l1_to_l2_task() )
    call this%l1_to_l2_context%gather_to_master(input_data_size, input_data, recv_counts, displs, output_data )
  end subroutine environment_l2_from_l1_gatherv_ip_1D_array

  !=============================================================================
  subroutine environment_l2_from_l1_gatherv_igp_1D_array ( this, input_data_size, input_data, recv_counts, displs, output_data )
    implicit none
    class(environment_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data_size
    integer(igp)            , intent(in)   :: input_data(input_data_size)
    integer(ip)             , intent(in)   :: recv_counts(:) ! ( this%l1_to_l2_context%get_num_tasks())
    integer(ip)             , intent(in)   :: displs(:) ! ( this%l1_to_l2_context%get_num_tasks())
    integer(igp)            , intent(out)  :: output_data(:)
    assert ( this%am_i_l1_to_l2_task() )
    call this%l1_to_l2_context%gather_to_master(input_data_size, input_data, recv_counts, displs, output_data )
  end subroutine environment_l2_from_l1_gatherv_igp_1D_array

  !=============================================================================
  subroutine environment_l2_from_l1_gatherv_rp_1D_array ( this, input_data_size, input_data, recv_counts, displs, output_data )
    implicit none
    class(environment_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data_size
    real(rp)                , intent(in)   :: input_data(input_data_size)
    integer(ip)             , intent(in)   :: recv_counts(:) ! ( this%l1_to_l2_context%get_num_tasks())
    integer(ip)             , intent(in)   :: displs(:) ! ( this%l1_to_l2_context%get_num_tasks())
    real(rp)                , intent(out)  :: output_data(:)
    assert ( this%am_i_l1_to_l2_task() )
    call this%l1_to_l2_context%gather_to_master(input_data_size, input_data, recv_counts, displs, output_data )
  end subroutine environment_l2_from_l1_gatherv_rp_1D_array

  !=============================================================================
  subroutine environment_l2_from_l1_gatherv_rp_2D_array ( this, input_data, recv_counts, displs, output_data )
    implicit none
    class(environment_t), intent(in)   :: this
    real(rp)                , intent(in)   :: input_data(:,:)
    integer(ip)             , intent(in)   :: recv_counts(:) ! ( this%l1_to_l2_context%get_num_tasks())
    integer(ip)             , intent(in)   :: displs(:) ! ( this%l1_to_l2_context%get_num_tasks())
    real(rp)                , intent(out)  :: output_data(:)
    assert ( this%am_i_l1_to_l2_task() )
    call this%l1_to_l2_context%gather_to_master(input_data, recv_counts, displs, output_data )
  end subroutine environment_l2_from_l1_gatherv_rp_2D_array

  !=============================================================================
  subroutine environment_l2_to_l1_scatterv_rp_1D_array ( this, input_data, send_counts, displs, output_data_size, output_data )
    implicit none
    class(environment_t), intent(in)   :: this
    real(rp)                , intent(in)   :: input_data(:)
    integer(ip)             , intent(in)   :: send_counts(:) ! ( this%l1_to_l2_context%get_num_tasks())
    integer(ip)             , intent(in)   :: displs(:) ! ( this%l1_to_l2_context%get_num_tasks())
    integer(ip)             , intent(in)   :: output_data_size
    real(rp)                , intent(out)  :: output_data(output_data_size)
    assert( this%am_i_l1_to_l2_task() )
    call this%l1_to_l2_context%scatter_from_master (input_data, send_counts, displs, output_data_size, output_data )
  end subroutine environment_l2_to_l1_scatterv_rp_1D_array

end module environment_names
