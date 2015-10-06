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
module par_scalar_array_names
  ! Serial modules
  use types_names
  use memor_names
  use stdio_names
  use serial_scalar_array_names
  use array_names
  use vector_names

#ifdef ENABLE_BLAS       
  use blas77_interfaces_names
#endif

#ifdef memcheck       
  use iso_c_binding
#endif

  ! Parallel modules
  use par_context_names
  use par_environment_names
  use dof_distribution_names
  use psb_penv_mod_names

  ! Abstract types
  use vector_names

# include "debug.i90"

  implicit none
  private

  integer(ip), parameter :: undefined     = -1 ! Undefined. State to be assigned by the user or externally.
  integer(ip), parameter :: part_summed   = 0  ! partially summed element-based vector
  integer(ip), parameter :: full_summed   = 1  ! fully     summed element-based vector

  ! Distributed Vector
  type, extends(array_t) :: par_scalar_array_t
     ! Data structure which stores the local part 
     ! of the vector mapped to the current processor.
     type( serial_scalar_array_t ) :: serial_scalar_array
     ! Partially or fully summed
     integer(ip)  :: state 
     ! Parallel DoF distribution control info.
     type ( dof_distribution_t ), pointer  :: dof_dist => NULL()
     type ( par_environment_t ) , pointer  :: p_env => NULL()
   contains
     procedure :: create_and_allocate => par_scalar_array_create_and_allocate
     procedure :: create              => par_scalar_array_create
     procedure :: allocate            => par_scalar_array_allocate
     procedure :: create_view         => par_scalar_array_create_view
     procedure :: weight              => par_scalar_array_weight
     procedure :: print               => par_scalar_array_print
     procedure :: print_market_market => par_scalar_array_print_matrix_market

     ! Provide type bound procedures (tbp) implementors
     procedure :: dot  => par_scalar_array_dot
     procedure :: copy => par_scalar_array_copy
     procedure :: init => par_scalar_array_init
     procedure :: scal => par_scalar_array_scal
     procedure :: axpby => par_scalar_array_axpby
     procedure :: nrm2 => par_scalar_array_nrm2
     procedure :: clone => par_scalar_array_clone
     procedure :: comm  => par_scalar_array_comm
     procedure :: same_vector_space => par_scalar_array_same_vector_space
     procedure :: free_in_stages  => par_scalar_array_free_in_stages
     procedure :: default_initialization  => par_scalar_array_default_init
  end type par_scalar_array_t


  ! Types
  public :: par_scalar_array_t

  ! Constants 
  public :: undefined, part_summed, full_summed


!***********************************************************************
! Allocatable arrays of type(par_scalar_array_t)
!***********************************************************************
# define var_attr allocatable, target
# define point(a,b) call move_alloc(a,b)
# define generic_status_test             allocated
# define generic_memalloc_interface      memalloc
# define generic_memrealloc_interface    memrealloc
# define generic_memfree_interface       memfree
# define generic_memmovealloc_interface  memmovealloc

# define var_type type(par_scalar_array_t)
# define var_size 52
# define bound_kind ip
# include "mem_header.i90"

  public :: memalloc,  memrealloc,  memfree, memmovealloc

contains

# include "mem_body.i90"

  !=============================================================================
  subroutine par_scalar_array_default_init (this)
    implicit none
    class(par_scalar_array_t), intent(inout) :: this

    call this%serial_scalar_array%default_initialization()
    nullify(this%dof_dist)
    nullify(this%p_env)
    call this%NullifyTemporary()
  end subroutine par_scalar_array_default_init

  !=============================================================================
  subroutine par_scalar_array_create_and_allocate (this, dof_dist, p_env)
    implicit none
    ! Parameters
    class(par_scalar_array_t), intent(out) :: this
    type(dof_distribution_t) , intent(in)  :: dof_dist
    type(par_environment_t)  , intent(in)  :: p_env

    call this%create(dof_dist, p_env)
    call this%allocate()
  end subroutine par_scalar_array_create_and_allocate
  
  !=============================================================================
  subroutine par_scalar_array_create (this, dof_dist, p_env)
    implicit none
    ! Parameters
    class(par_scalar_array_t)       , intent(out) :: this
    type(dof_distribution_t), target, intent(in)  :: dof_dist
    type(par_environment_t) , target, intent(in)  :: p_env

    ! p_env%p_context is required within this subroutine
    assert ( associated(p_env%p_context) )
    assert ( p_env%p_context%created .eqv. .true.)
    this%dof_dist => dof_dist
    this%p_env    => p_env 
    this%state = undefined
    if(this%p_env%p_context%iam<0) return
    call this%serial_scalar_array%create (dof_dist%nl)
  end subroutine par_scalar_array_create
  
  !=============================================================================
  subroutine par_scalar_array_allocate (this)
    implicit none
    ! Parameters
    class(par_scalar_array_t), intent(inout) :: this
    ! p_env%p_context is required within this subroutine
    assert ( associated(this%p_env%p_context) )
    assert ( this%p_env%p_context%created .eqv. .true.) 
    assert ( associated(this%dof_dist) )
    if(this%p_env%p_context%iam<0) return
    call this%serial_scalar_array%allocate ()
  end subroutine par_scalar_array_allocate

  !=============================================================================
  subroutine par_scalar_array_create_view (this, start, end, t_p_vec)
    implicit none
    class(par_scalar_array_t), intent(in), target :: this
    integer(ip)     , intent(in)         :: start
    integer(ip)     , intent(in)         :: end
    type(par_scalar_array_t), intent(out)        :: t_p_vec

    ! The routine requires the partition/context info
    assert ( associated( this%dof_dist ) )
    assert ( associated( this%p_env%p_context ) )
    assert ( this%p_env%p_context%created .eqv. .true.)

    ! Associate dof distribution and parallel environment 
    t_p_vec%dof_dist => this%dof_dist
    t_p_vec%p_env    => this%p_env

    assert ( this%state /= undefined )
    t_p_vec%state = this%state

    if(this%p_env%p_context%iam<0) return

    ! Call vector_create_view
    call this%serial_scalar_array%create_view (start, end, t_p_vec%serial_scalar_array) 

  end subroutine par_scalar_array_create_view

  !=============================================================================
  ! VERY IMPORTANT: comm_interface is well-defined if and only if p_vec is a 
  !                 vector on the interface
  !=============================================================================
  subroutine comm_interface (p_vec)
    use par_sparse_global_collectives_names
    implicit none

    ! Parameters
    type(par_scalar_array_t), intent(inout)         :: p_vec

    ! Pointer to part/context object is required
    assert ( associated(p_vec%dof_dist) )
    assert ( associated(p_vec%p_env%p_context) )
    assert ( p_vec%p_env%p_context%created .eqv. .true.)
    if(p_vec%p_env%p_context%iam<0) return


    ! call import_print (6, p_vec%dof_dist%dof_import)
    ! write(*,*) 'SE 1', size(p_vec%serial_scalar_array%b)
    ! call vector_print ( 6, p_vec%serial_scalar_array )

    ! First stage: owners receive/reduce, non-owners send
    call single_exchange ( p_vec%p_env%p_context%icontxt, &
         p_vec%dof_dist%dof_import%num_rcv,    &
         p_vec%dof_dist%dof_import%list_rcv,   &
         p_vec%dof_dist%dof_import%rcv_ptrs,   &
         p_vec%dof_dist%dof_import%unpack_idx, &
         p_vec%dof_dist%dof_import%num_snd,    &
         p_vec%dof_dist%dof_import%list_snd,   &
         p_vec%dof_dist%dof_import%snd_ptrs,   &
         p_vec%dof_dist%dof_import%pack_idx,   &
         1.0_rp,                         &
         1.0_rp,                         &
         p_vec%serial_scalar_array%b ) 

    ! write(*,*) 'SE 2'

    ! Second stage: owners send, non-owners receive/insert
    call single_exchange ( p_vec%p_env%p_context%icontxt, &
         p_vec%dof_dist%dof_import%num_snd,    &
         p_vec%dof_dist%dof_import%list_snd,   &
         p_vec%dof_dist%dof_import%snd_ptrs,   &
         p_vec%dof_dist%dof_import%pack_idx,   &
         p_vec%dof_dist%dof_import%num_rcv,    &
         p_vec%dof_dist%dof_import%list_rcv,   &
         p_vec%dof_dist%dof_import%rcv_ptrs,   &
         p_vec%dof_dist%dof_import%unpack_idx, &
         1.0_rp,                         &
         0.0_rp,                         &
         p_vec%serial_scalar_array%b )

  end subroutine comm_interface

  subroutine par_scalar_array_weight ( this, weight )
    implicit none
    ! Parameters
    class(par_scalar_array_t), intent(inout) :: this
    real(rp)        , intent(in), target, optional :: weight(*)

    ! Local variables
    integer(ip)      :: ni 
    type(par_scalar_array_t) :: p_vec_G

    ! Pointer to part/context object is required
    assert ( associated(this%dof_dist) )
    assert ( associated(this%p_env%p_context) )
    assert ( this%p_env%p_context%created .eqv. .true.)
    if(this%p_env%p_context%iam<0) return


    ni = this%serial_scalar_array%size - this%dof_dist%nb
    call this%create_view ( ni+1, this%serial_scalar_array%size, p_vec_G )
    call weight_interface ( p_vec_G, weight )

    this%state = part_summed
  end subroutine par_scalar_array_weight

  !=============================================================================
  ! VERY IMPORTANT: weight_interface is well-defined if and only if p_vec is a 
  !                 vector on the interface
  !=============================================================================
  subroutine weight_interface ( p_vec, weight )
    implicit none
    ! Parameters
    type(par_scalar_array_t), intent(inout)                :: p_vec
    real(rp)        , intent(in), target, optional :: weight(*)

    ! Local variables
    integer(ip) :: iobj, i, i1, i2 
    real(rp)    :: weigt

    ! Pointer to part/context object is required
    assert ( associated(p_vec%dof_dist) )
    assert ( associated(p_vec%p_env%p_context) )
    assert ( p_vec%p_env%p_context%created .eqv. .true.)
    if(p_vec%p_env%p_context%iam<0) return


    if ( present(weight) ) then
       i1 = 1
       i2 = p_vec%dof_dist%nb
       do i=i1,i2
          p_vec%serial_scalar_array%b(i) =  p_vec%serial_scalar_array%b(i) * weight(i)
       end do
    else 

       do iobj=2, p_vec%dof_dist%nobjs
          weigt=1.0_rp/real(p_vec%dof_dist%lobjs(4,iobj))
          i1 = p_vec%dof_dist%lobjs(2,iobj) - p_vec%dof_dist%ni
          i2 = p_vec%dof_dist%lobjs(3,iobj) - p_vec%dof_dist%ni
#ifdef ENABLE_BLAS
          call dscal ( (i2-i1+1), weigt, p_vec%serial_scalar_array%b(i1:i2), 1 )
#else
          p_vec%serial_scalar_array%b(i1:i2) = weigt * p_vec%serial_scalar_array%b(i1:i2)
#endif 
       end do

    end if
  end subroutine weight_interface

  !=============================================================================
  ! VERY IMPORTANT: dot_interface is well-defined if and only if x/y are 
  !                 vectors on the interface
  !============================================================================= 
  subroutine dot_interface(x,y,t)
    implicit none
    ! Parameters  
    type(par_scalar_array_t), intent(in)  :: x
    type(par_scalar_array_t), intent(in)  :: y
    real(rp)    , intent(out)     :: t

    ! Locals 
    integer(ip)                :: ierrc
    type(par_scalar_array_t)              :: ws_vec

    ! Pointer to part/context object is required
    assert ( associated(x%dof_dist   ) )
    assert ( associated(x%p_env%p_context) )
    assert ( associated(y%dof_dist   ) )
    assert ( associated(y%p_env%p_context) )

    assert ( x%state /= undefined .and. y%state /= undefined  )

    if ( (x%state == part_summed .and. y%state == full_summed) .or. (x%state == full_summed .and. y%state == part_summed) ) then
       ! Perform local dot products
       t=x%serial_scalar_array%dot(y%serial_scalar_array)
    else if ( (x%state == full_summed .and. y%state == full_summed) ) then
       ! Perform local weighted dot products
       call weighted_dot (x, y, t)
    else if ( (x%state == part_summed .and. y%state == part_summed) ) then
       ! Allocate space for ws_vec
       call ws_vec%clone(x)
       ! ws_vec <- x  
       call ws_vec%copy(x)
       ! Transform ws_vec from partially summed to fully summed
       call comm_interface ( ws_vec )
       t=x%serial_scalar_array%dot(ws_vec%serial_scalar_array)
       call ws_vec%free()
    end if
  end subroutine dot_interface

  subroutine par_scalar_array_print ( this, luout )
    implicit none
    class(par_scalar_array_t), intent(in) :: this
    integer(ip)              , intent(in) :: luout

    ! Pointer to part/context object is required
    assert ( associated(this%dof_dist   ) )
    assert ( associated(this%p_env%p_context) )
    assert ( this%p_env%p_context%created .eqv. .true.)
    if(this%p_env%p_context%iam<0) return
    write(luout,'(a)') '*** begin par_scalar_array_t data structure ***'
    call  dof_distribution_print (luout, this%dof_dist)
    call  this%serial_scalar_array%print(luout)
    write(luout,'(a)') '*** end par_scalar_array_t data structure ***'
  end subroutine par_scalar_array_print

  subroutine par_scalar_array_print_matrix_market ( this, dir_path, prefix )
    implicit none
    ! Parameters
    class(par_scalar_array_t), intent(in)  :: this
    character (*)            , intent(in)  :: dir_path
    character (*)            , intent(in)  :: prefix

    ! Locals
    integer         :: iam, num_procs, lunou
    integer(ip)     :: j, ndigs_iam, ndigs_num_procs, id_map
    character(256)  :: name 
    character(256)  :: zeros
    character(256)  :: part_id

    assert ( associated(this%dof_dist   ) )
    assert ( associated(this%p_env%p_context) )
    assert ( this%p_env%p_context%created .eqv. .true.)
    if(this%p_env%p_context%iam<0) return

    name = trim(prefix) // '.par_vector' // '.mtx'

    ! Get context info
    call par_context_info ( this%p_env%p_context, iam, num_procs )

    ! Form the file_path of the partition object to be read
    iam = iam + 1 ! Partition identifers start from 1 !!

    ndigs_num_procs = count_digits_par_vector (num_procs)
    zeros = ' '   
    ndigs_iam = count_digits_par_vector ( iam )

    ! write(*,*) ndgs_num_procs, ndigs_iam DBG

    do j=1,  ndigs_num_procs - ndigs_iam
       zeros (j:j) = '0'
    end do
    part_id = ch(iam)

    ! Read partition data from path_file file
    lunou =  io_open (trim(dir_path) // '/' // trim(name) // '.' // trim(zeros) // trim(part_id), 'write')       
    call this%serial_scalar_array%print_matrix_market (lunou)
    call io_close (lunou)
  end subroutine par_scalar_array_print_matrix_market

  function count_digits_par_vector ( i )
    implicit none
    ! Parameters
    integer(ip), intent(in) :: i 
    integer(ip)             :: count_digits_par_vector
    ! Locals   
    integer(ip)             :: x 
    x = i 
    if (x < 0) x = -x;
    count_digits_par_vector = 1;
    x = x/10;
    do while( x > 0)
       count_digits_par_vector = count_digits_par_vector + 1
       x = x/10;
    end do
  end function count_digits_par_vector


  ! Auxiliary (module private) routine
  ! for computing dot products and euclidean
  ! norms for element-based data distributions 
  subroutine weighted_dot (x,y,t)
    implicit none
    ! Parameters  
    type(par_scalar_array_t), intent(in)  :: x
    type(par_scalar_array_t), intent(in)  :: y
    real(rp)    , intent(out)     :: t

    ! Local variables
    integer(ip)                   :: iobj, i, id, i1, i2 
    real(rp)                      :: weigt

    ! Pointer to part/context object is required
    assert ( associated(x%dof_dist   ) )
    assert ( associated(x%p_env%p_context) )
    assert ( associated(y%dof_dist   ) )
    assert ( associated(y%p_env%p_context) )

    t = 0.0_rp 
    do iobj=2, x%dof_dist%nobjs
       weigt=1.0_rp/real(x%dof_dist%lobjs(4,iobj))
       i1=x%dof_dist%lobjs(2,iobj) - x%dof_dist%ni
       i2=x%dof_dist%lobjs(3,iobj) - x%dof_dist%ni
       call flat_weighted_dot ( i1, i2, & 
            & x%serial_scalar_array%b, y%serial_scalar_array%b, weigt, t) 
    end do
  end subroutine weighted_dot

  ! Auxiliary (module private) routine
  ! for computing dot products and euclidean
  ! norms for element-based data distributions 
  subroutine flat_weighted_dot ( i1, i2, x, y, weigt, t)
    implicit none
    ! Parameters  
    integer(ip), intent(in)    :: i1, i2
    real(rp),    intent(in)    :: x(*), y(*), weigt
    real(rp),    intent(inout) :: t

    ! Locals   
    integer(ip)                :: i, id

    do i=i1,i2
       t = t + weigt * x(i) * y(i)
    end do
  end subroutine flat_weighted_dot

  ! alpha <- op1^T * op2
  function par_scalar_array_dot(op1,op2) result(alpha)
    implicit none
    class(par_scalar_array_t), intent(in)    :: op1
    class(vector_t), intent(in)  :: op2
    real(rp) :: alpha

    ! Locals 
    integer(ip)                   :: ni 
    type(par_scalar_array_t)              :: x_I, x_G, y_I, y_G
    real(rp)                      :: s

    ! Alpha should be defined in all tasks (not only in coarse-grid ones)
    alpha = 0.0_rp

    call op1%GuardTemp()
    call op2%GuardTemp()
    select type(op2)
       class is (par_scalar_array_t)
          ! Pointer to part/context object is required
       assert ( associated(op1%dof_dist   ) )
       assert ( associated(op1%p_env%p_context) )
       assert ( associated(op2%dof_dist   ) )
       assert ( associated(op2%p_env%p_context) )
       assert ( op1%p_env%p_context%created .eqv. .true.)
       if(op1%p_env%p_context%iam<0) return

       assert ( op1%state /= undefined .and. op2%state /= undefined  )

       ni = op1%serial_scalar_array%size - op1%dof_dist%nb
       if ( ni > 0 ) then
          call op1%create_view (1, ni, x_I)
          call op2%create_view (1, ni, y_I)
          alpha = x_I%serial_scalar_array%dot ( y_I%serial_scalar_array )
       else
          alpha = 0.0_rp
       end if

       call op1%create_view(ni+1, op1%serial_scalar_array%size, x_G )
       call op2%create_view(ni+1, op2%serial_scalar_array%size, y_G )
       call dot_interface          ( x_G, y_G, s )

       alpha = alpha + s

       ! Reduce-sum local dot products on all processes
       call psb_sum ( op2%p_env%p_context%icontxt, alpha )
       class default
       write(0,'(a)') 'par_scalar_array_t%dot: unsupported op2 class'
       check(1==0)
    end select
    call op1%CleanTemp()
    call op2%CleanTemp()
  end function par_scalar_array_dot

  ! op1 <- op2 
  subroutine par_scalar_array_copy(op1,op2)
    implicit none
    class(par_scalar_array_t), intent(inout) :: op1
    class(vector_t), intent(in)  :: op2

    call op2%GuardTemp()
    select type(op2)
       class is (par_scalar_array_t)
          ! Pointer to part/context object is required
       assert ( associated(op2%dof_dist   ) )
       assert ( associated(op2%p_env%p_context) )
       assert ( associated(op1%dof_dist   ) )
       assert ( associated(op1%p_env%p_context) )
       assert ( op2%p_env%p_context%created .eqv. .true.)
       if(op2%p_env%p_context%iam<0) return

       assert ( op2%state /= undefined  )
       ! Perform local copy
       call op1%serial_scalar_array%copy ( op2%serial_scalar_array )
       op1%state = op2%state
       class default
       write(0,'(a)') 'par_scalar_array_t%copy: unsupported op2 class'
       check(1==0)
    end select
    call op2%CleanTemp()
  end subroutine par_scalar_array_copy

  ! op1 <- alpha * op2
  subroutine par_scalar_array_scal(op1,alpha,op2)
    implicit none
    class(par_scalar_array_t), intent(inout) :: op1
    real(rp), intent(in) :: alpha
    class(vector_t), intent(in) :: op2

    call op2%GuardTemp()
    select type(op2)
       class is (par_scalar_array_t)
          ! Pointer to part/context object is required
       assert ( associated(op2%dof_dist   ) )
       assert ( associated(op2%p_env%p_context) )
       assert ( associated(op1%dof_dist   ) )
       assert ( associated(op1%p_env%p_context) )
       assert ( op1%p_env%p_context%created .eqv. .true.)
       if(op1%p_env%p_context%iam<0) return

       ! Check matching partition/handler
       assert ( op2%state /= undefined )

       ! Scal local copy
       call op1%serial_scalar_array%scal ( alpha, op2%serial_scalar_array )
       op1%state = op2%state
       class default
       write(0,'(a)') 'par_scalar_array_t%scal: unsupported op2 class'
       check(1==0)
    end select
    call op2%CleanTemp()
  end subroutine par_scalar_array_scal
  ! op <- alpha
  subroutine par_scalar_array_init(op,alpha)
    implicit none
    class(par_scalar_array_t), intent(inout) :: op
    real(rp), intent(in) :: alpha

    ! Pointer to part/context object is required
    assert ( associated(op%dof_dist   ) )
    assert ( associated(op%p_env%p_context) )
    assert ( op%p_env%p_context%created .eqv. .true.)
    if(op%p_env%p_context%iam<0) return

    ! Init local copy
    call op%serial_scalar_array%init(alpha)
    ! op%state = full_summed
  end subroutine par_scalar_array_init

  ! op1 <- alpha*op2 + beta*op1
  subroutine par_scalar_array_axpby(op1, alpha, op2, beta)
    implicit none
    class(par_scalar_array_t), intent(inout) :: op1
    real(rp), intent(in) :: alpha
    class(vector_t), intent(in) :: op2
    real(rp), intent(in) :: beta

    call op2%GuardTemp()
    select type(op2)
       class is (par_scalar_array_t)
          ! Pointer to part/context object is required
       assert ( associated(op2%dof_dist   ) )
       assert ( associated(op2%p_env%p_context) )
       assert ( associated(op1%dof_dist   ) )
       assert ( associated(op1%p_env%p_context) )
       assert ( op1%p_env%p_context%created )
       if(op1%p_env%p_context%iam<0) return

       ! Check matching partition/handler
       assert ( op1%state /= undefined .and. op2%state /= undefined  )
       assert ( op1%state == op2%state )

       call op1%serial_scalar_array%axpby( alpha, op2%serial_scalar_array, beta )
       class default
       write(0,'(a)') 'par_scalar_array_t%axpby: unsupported op2 class'
       check(1==0)
    end select
    call op2%CleanTemp()
  end subroutine par_scalar_array_axpby

  ! alpha <- nrm2(op)
  function par_scalar_array_nrm2(op) result(alpha)
    implicit none
    class(par_scalar_array_t), intent(in)  :: op
    real(rp) :: alpha

    ! Alpha should be defined in all tasks (not only in coarse-grid ones)
    alpha = 0.0_rp

    call op%GuardTemp()
    select type(op)
       class is (par_scalar_array_t)
          ! p_env%p_context is required within this subroutine 
       assert ( associated(op%dof_dist) )
       assert ( associated(op%p_env%p_context) )
       assert ( op%p_env%p_context%created )
       if(op%p_env%p_context%iam<0) return
       assert ( op%state /= undefined )
       alpha = op%dot(op)
       alpha = sqrt(alpha)
       class default
       write(0,'(a)') 'par_scalar_array_t%nrm2: unsupported op2 class'
       check(1==0)
    end select
    call op%CleanTemp()
  end function par_scalar_array_nrm2

  ! op1 <- clone(op2) 
  subroutine par_scalar_array_clone(op1,op2)
    implicit none
    class(par_scalar_array_t)          , intent(inout) :: op1
    class(vector_t), target, intent(in)    :: op2

    call op2%GuardTemp()
    select type(op2)
       class is (par_scalar_array_t)
          ! p_env%p_context is required within this subroutine 
       assert ( associated(op2%dof_dist) )
       assert ( associated(op2%p_env%p_context) )
       assert ( op2%p_env%p_context%created )
       op1%dof_dist => op2%dof_dist
       op1%p_env    => op2%p_env
       op1%state    =  op2%state
       if(op2%p_env%p_context%iam<0) return
       call op1%serial_scalar_array%clone ( op2%serial_scalar_array )
       class default
       write(0,'(a)') 'par_scalar_array_t%clone: unsupported op2 class'
       check(1==0)
    end select
    call op2%CleanTemp()
  end subroutine par_scalar_array_clone

  ! op <- comm(op)
  subroutine par_scalar_array_comm(op)
    implicit none
    class(par_scalar_array_t), intent(inout) :: op

    ! Local variables
    integer(ip)      :: ni 
    type(par_scalar_array_t) :: op_G

    ! Pointer to part/context object is required
    assert ( associated(op%dof_dist) )
    assert ( associated(op%p_env%p_context) )
    assert ( op%p_env%p_context%created .eqv. .true.)
    if(op%p_env%p_context%iam<0) return
    ni = op%serial_scalar_array%size - op%dof_dist%nb
    call op%create_view(ni+1, op%serial_scalar_array%size, op_G)
    call comm_interface ( op_G )
    op%state = full_summed
  end subroutine par_scalar_array_comm

  subroutine par_scalar_array_free_in_stages(this,action)
    implicit none
    class(par_scalar_array_t), intent(inout) :: this
    integer(ip)              , intent(in)    :: action

    ! The routine requires the partition/context info
    assert ( associated( this%dof_dist ) )
    assert ( associated( this%p_env%p_context ) )
    assert ( this%p_env%p_context%created .eqv. .true.)
    assert ( action == free_clean .or. action == free_struct .or. action == free_values )	 

    if(this%p_env%p_context%iam<0) then 
       if ( action == free_clean ) then
          this%state = undefined
          nullify ( this%dof_dist )
          nullify ( this%p_env )
       end if
       return
    end if

    ! Free local part
    call this%serial_scalar_array%free_in_stages(action)

    if ( action == free_clean ) then
       this%state = undefined
       nullify ( this%dof_dist )
       nullify ( this%p_env )
    end if
  end subroutine par_scalar_array_free_in_stages
  
  function par_scalar_array_same_vector_space(this,vector)
   implicit none
   class(par_scalar_array_t), intent(in) :: this
   class(vector_t)             , intent(in) :: vector
   logical :: par_scalar_array_same_vector_space
   par_scalar_array_same_vector_space = .false.
   select type(vector)
   class is (par_scalar_array_t)
     par_scalar_array_same_vector_space = (this%dof_dist%nl == vector%dof_dist%nl)
   end select
 end function par_scalar_array_same_vector_space
  
  
end module par_scalar_array_names
