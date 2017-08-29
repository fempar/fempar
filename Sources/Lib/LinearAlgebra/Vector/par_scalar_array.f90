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

#ifdef ENABLE_BLAS       
  use blas77_interfaces_names
#endif

#ifdef memcheck       
  use iso_c_binding
#endif

  ! Parallel modules
  use environment_names
  use dof_import_names

  ! Abstract types
  use vector_names
  use array_names

  implicit none
# include "debug.i90"
  private

  ! Distributed-memory scalar array
  type, extends(array_t) :: par_scalar_array_t
     private
     type ( serial_scalar_array_t )        :: serial_scalar_array
     type ( environment_t ) , pointer  :: p_env      => NULL()     
     type ( dof_import_t )      , pointer  :: dof_import => NULL()
   contains
     procedure :: create_and_allocate     => par_scalar_array_create_and_allocate
     procedure :: create                  => par_scalar_array_create
     procedure :: allocate                => par_scalar_array_allocate
     procedure :: create_view             => par_scalar_array_create_view
     procedure :: set_view_entries        => par_scalar_array_set_view_entries
     procedure :: print                   => par_scalar_array_print
     procedure :: print_matrix_market     => par_scalar_array_print_matrix_market
     procedure :: get_par_environment     => par_scalar_array_get_par_environment
     procedure :: get_serial_scalar_array => par_scalar_array_get_serial_scalar_array
    
     procedure, private :: par_scalar_array_insert_single_entry
     procedure, private :: par_scalar_array_insert_multiple_entries
     procedure, private :: par_scalar_array_add_single_entry
     procedure, private :: par_scalar_array_add_multiple_entries                      
    
     generic :: insert                  => par_scalar_array_insert_single_entry, &
                                             par_scalar_array_insert_multiple_entries
     generic :: add                     => par_scalar_array_add_single_entry, &
                                             par_scalar_array_add_multiple_entries
     
     ! Provide type bound procedures (tbp) implementors
     procedure :: dot                    => par_scalar_array_dot
     procedure :: local_dot              => par_scalar_array_local_dot
     procedure :: copy                   => par_scalar_array_copy
     procedure :: init                   => par_scalar_array_init
     procedure :: scal                   => par_scalar_array_scal
     procedure :: axpby                  => par_scalar_array_axpby
     procedure :: nrm2                   => par_scalar_array_nrm2
     procedure :: clone                  => par_scalar_array_clone
     procedure :: comm                   => par_scalar_array_comm
     procedure :: same_vector_space      => par_scalar_array_same_vector_space
     procedure :: free_in_stages         => par_scalar_array_free_in_stages
     procedure :: default_initialization => par_scalar_array_default_init
     procedure :: get_num_blocks      => par_scalar_array_get_num_blocks
     procedure :: extract_subvector      => par_scalar_array_extract_subvector
     procedure :: insert_subvector       => par_scalar_array_insert_subvector
  end type par_scalar_array_t


  ! Types
  public :: par_scalar_array_t


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
    nullify(this%dof_import)
    nullify(this%p_env)
    call this%NullifyTemporary()
  end subroutine par_scalar_array_default_init

  !=============================================================================
  subroutine par_scalar_array_create_and_allocate (this, p_env, dof_import)
    implicit none
    ! Parameters
    class(par_scalar_array_t), intent(inout) :: this
    type(environment_t)  , intent(in)    :: p_env
    type(dof_import_t)       , intent(in)    :: dof_import
    call this%free()
    call this%create(p_env, dof_import)
    call this%allocate()
  end subroutine par_scalar_array_create_and_allocate
  
  !=============================================================================
  subroutine par_scalar_array_create (this, p_env, dof_import)
    implicit none
    ! Parameters
    class(par_scalar_array_t)       , intent(inout) :: this
    type(environment_t) , target, intent(in)    :: p_env
    type(dof_import_t)      , target, intent(in)    :: dof_import
    call this%free()
    assert ( p_env%created() ) 
    this%p_env      => p_env 
    this%dof_import => dof_import
    if(.not. this%p_env%am_i_l1_task()) return
    call this%serial_scalar_array%create (dof_import%get_num_dofs())
  end subroutine par_scalar_array_create
  
  !=============================================================================
  subroutine par_scalar_array_allocate (this)
    implicit none
    ! Parameters
    class(par_scalar_array_t), intent(inout) :: this
    if(.not. this%p_env%am_i_l1_task()) return
    call this%serial_scalar_array%allocate ()
  end subroutine par_scalar_array_allocate

  !=============================================================================
  subroutine par_scalar_array_create_view (this, start, end, t_p_vec)
    implicit none
    class(par_scalar_array_t),  target , intent(in)     :: this
    integer(ip)                        , intent(in)     :: start
    integer(ip)                        , intent(in)     :: end
    type(par_scalar_array_t)           , intent(inout)  :: t_p_vec
    
    assert ( associated(this%p_env) )
    assert ( this%p_env%created() )
    assert ( associated(this%dof_import) )
    
    call t_p_vec%free()
    
    t_p_vec%dof_import => this%dof_import
    t_p_vec%p_env    => this%p_env
    if(.not. this%p_env%am_i_l1_task()) return
    ! Call vector_create_view
    call this%serial_scalar_array%create_view (start, end, t_p_vec%serial_scalar_array) 
  end subroutine par_scalar_array_create_view
  
  !=============================================================================
  subroutine par_scalar_array_set_view_entries ( this, entries ) 
    implicit none
    class(par_scalar_array_t), intent(inout) ::  this
    real(rp)                 ,  intent(in)   ::  entries(:)
    if(.not. this%p_env%am_i_l1_task()) return
    call this%serial_scalar_array%set_view_entries (entries) 
  end subroutine par_scalar_array_set_view_entries

  !=============================================================================
  function par_scalar_array_get_par_environment( this )
   implicit none
   class(par_scalar_array_t), target, intent(in) :: this 
   type(environment_t)  , pointer            :: par_scalar_array_get_par_environment
   par_scalar_array_get_par_environment => this%p_env
  end function par_scalar_array_get_par_environment
  
  !=============================================================================
  function par_scalar_array_get_serial_scalar_array( this )
   implicit none
   class(par_scalar_array_t), target, intent(in) :: this 
   type(serial_scalar_array_t)  , pointer        :: par_scalar_array_get_serial_scalar_array
   par_scalar_array_get_serial_scalar_array => this%serial_scalar_array
  end function par_scalar_array_get_serial_scalar_array
  
  subroutine par_scalar_array_insert_single_entry (this, i, val)
    implicit none
    class(par_scalar_array_t), intent(inout) :: this
    integer(ip)                 , intent(in)    :: i
    real(rp)                    , intent(in)    :: val
    if(.not. this%p_env%am_i_l1_task()) return
    call this%serial_scalar_array%insert(i,val)
  end subroutine par_scalar_array_insert_single_entry
  
  subroutine par_scalar_array_insert_multiple_entries (this, num_entries, ia, ioffset, val)
    implicit none
    class(par_scalar_array_t), intent(inout) :: this
    integer(ip)                 , intent(in) :: num_entries
    integer(ip)                 , intent(in) :: ia(num_entries)
    integer(ip)                 , intent(in) :: ioffset
    real(rp)                    , intent(in) :: val(:)
    if(.not. this%p_env%am_i_l1_task()) return
    call this%serial_scalar_array%insert(num_entries,ia,ioffset,val)
  end subroutine par_scalar_array_insert_multiple_entries
  
  subroutine par_scalar_array_add_single_entry (this, i, val)
    implicit none
    class(par_scalar_array_t), intent(inout) :: this
    integer(ip)                 , intent(in)    :: i
    real(rp)                    , intent(in)    :: val
    if(.not. this%p_env%am_i_l1_task()) return
    call this%serial_scalar_array%add(i,val)
  end subroutine par_scalar_array_add_single_entry
  
  subroutine par_scalar_array_add_multiple_entries (this, num_entries, ia, ioffset, val)
    implicit none
    class(par_scalar_array_t), intent(inout) :: this
    integer(ip)                 , intent(in) :: num_entries
    integer(ip)                 , intent(in) :: ia(num_entries)
    integer(ip)                 , intent(in) :: ioffset
    real(rp)                    , intent(in) :: val(:)
    if(.not. this%p_env%am_i_l1_task()) return
    call this%serial_scalar_array%add(num_entries,ia,ioffset,val)
  end subroutine par_scalar_array_add_multiple_entries
  
  !=============================================================================
  subroutine par_scalar_array_print ( this, luout )
    implicit none
    class(par_scalar_array_t), intent(in) :: this
    integer(ip)              , intent(in) :: luout
    if(.not. this%p_env%am_i_l1_task()) return
    write(luout,'(a)') '*** begin par_scalar_array_t data structure ***'
    call  this%dof_import%print(luout)
    call  this%serial_scalar_array%print(luout)
    write(luout,'(a)') '*** end par_scalar_array_t data structure ***'
  end subroutine par_scalar_array_print
  
  !=============================================================================
  subroutine par_scalar_array_print_matrix_market ( this, luout )
    implicit none
    ! Parameters
    class(par_scalar_array_t), intent(in)  :: this
    integer(ip)              , intent(in)  :: luout
    if(.not. this%p_env%am_i_l1_task()) return
    call this%serial_scalar_array%print_matrix_market (luout)
  end subroutine par_scalar_array_print_matrix_market

  !=============================================================================
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



  ! alpha <- op1^T * op2
  function par_scalar_array_dot(op1,op2) result(alpha)
    implicit none
    class(par_scalar_array_t), intent(in)  :: op1
    class(vector_t)          , intent(in)  :: op2
    real(rp)                               :: alpha
    
    ! Alpha should be defined in all tasks (not only in coarse-grid ones)
    alpha = 0.0_rp

    call op1%GuardTemp()
    call op2%GuardTemp()
    select type(op2)
       class is (par_scalar_array_t)
       if(.not. op1%p_env%am_i_l1_task()) return
       
       alpha = flat_weighted_dot ( op1%serial_scalar_array%get_entries(), &
                                   op2%serial_scalar_array%get_entries(), &
                                   op1%dof_import%get_weight() )
    
       ! Reduce-sum local dot products on all L1 processes
       call op2%p_env%l1_sum(alpha)
       
       class default
       write(0,'(a)') 'par_scalar_array_t%dot: unsupported op2 class'
       check(1==0)
    end select
    call op1%CleanTemp()
    call op2%CleanTemp()
  end function par_scalar_array_dot
  
  ! Auxiliary (module private) routine
  ! for computing dot products and euclidean
  ! norms for element-based data distributions 
  function flat_weighted_dot ( x, y, weigt) result(t)
    implicit none
    ! Parameters
    real(rp), intent(in)    :: x(:), y(:), weigt(:)
    real(rp) :: t

    ! Locals   
    integer(ip)                :: i
    
    t=0.0_rp
    do i=1, size(x)
       t = t + weigt(i) * x(i) * y(i)
    end do
    
  end function flat_weighted_dot
  
    ! alpha <- op1^T * op2 without final allreduce
  function par_scalar_array_local_dot(op1,op2) result(alpha)
    implicit none
    class(par_scalar_array_t), intent(in)    :: op1
    class(vector_t), intent(in)  :: op2
    real(rp) :: alpha

    ! Alpha should be defined in all tasks (not only in coarse-grid ones)
    alpha = 0.0_rp

    call op1%GuardTemp()
    call op2%GuardTemp()
    select type(op2)
       class is (par_scalar_array_t)
       if(.not. op1%p_env%am_i_l1_task()) return
       alpha = flat_weighted_dot ( op1%serial_scalar_array%get_entries(), &
                                   op2%serial_scalar_array%get_entries(), &
                                   op1%dof_import%get_weight() )
       class default
       write(0,'(a)') 'par_scalar_array_t%dot: unsupported op2 class'
       check(1==0)
    end select
    call op1%CleanTemp()
    call op2%CleanTemp()
  end function par_scalar_array_local_dot
  
  ! op1 <- op2 
  subroutine par_scalar_array_copy(op1,op2)
    implicit none
    class(par_scalar_array_t), intent(inout) :: op1
    class(vector_t), intent(in)  :: op2

    call op2%GuardTemp()
    select type(op2)
       class is (par_scalar_array_t)
       if(.not. op2%p_env%am_i_l1_task()) return

       ! Perform local copy
       call op1%serial_scalar_array%copy ( op2%serial_scalar_array )
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
       if(.not. op1%p_env%am_i_l1_task()) return
       
       ! Scal local copy
       call op1%serial_scalar_array%scal ( alpha, op2%serial_scalar_array )
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

    if(.not. op%p_env%am_i_l1_task()) return

    ! Init local copy
    call op%serial_scalar_array%init(alpha)
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
       if(.not. op1%p_env%am_i_l1_task()) return
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
       if(.not. op%p_env%am_i_l1_task()) return
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
    class(par_scalar_array_t), target, intent(inout) :: op1
    class(vector_t)          , target, intent(in)    :: op2
    class(vector_t), pointer :: p
    p => op1
    if(associated(p,op2)) return ! It's aliasing
    
    call op2%GuardTemp()
    select type(op2)
       class is (par_scalar_array_t)
       call op1%free()
       call op1%create(op2%p_env, op2%dof_import)
       if(.not. op2%p_env%am_i_l1_task()) return
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
    real(rp)           , pointer             :: data(:)
    
    if(.not. op%p_env%am_i_l1_task()) return
    
    data => op%serial_scalar_array%get_entries()

    call op%p_env%l1_neighbours_exchange ( op%dof_import%get_num_rcv(),    &
                                           op%dof_import%get_list_rcv(),   &
                                           op%dof_import%get_rcv_ptrs(),   &
                                           op%dof_import%get_unpack_idx(), &
                                           op%dof_import%get_num_snd(),    &
                                           op%dof_import%get_list_snd(),   &
                                           op%dof_import%get_snd_ptrs(),   &
                                           op%dof_import%get_pack_idx(),   &
                                           1.0_rp,                         &
                                           1.0_rp,                         &
                                           data ) 

    ! Second stage: owners send, non-owners receive/insert
    call op%p_env%l1_neighbours_exchange ( op%dof_import%get_num_snd(),    &
                                           op%dof_import%get_list_snd(),   &
                                           op%dof_import%get_snd_ptrs(),   &
                                           op%dof_import%get_pack_idx(),   &
                                           op%dof_import%get_num_rcv(),    &
                                           op%dof_import%get_list_rcv(),   &
                                           op%dof_import%get_rcv_ptrs(),   &
                                           op%dof_import%get_unpack_idx(), &
                                           1.0_rp,                         &
                                           0.0_rp,                         &
                                           data )
  end subroutine par_scalar_array_comm

  !=============================================================================
  subroutine par_scalar_array_free_in_stages(this,action)
    implicit none
    class(par_scalar_array_t), intent(inout) :: this
    integer(ip)              , intent(in)    :: action

    assert ( action == free_clean .or. action == free_symbolic_setup .or. action == free_numerical_setup )	 

    if ( associated ( this%p_env ) ) then
      if(.not. this%p_env%am_i_l1_task()) then
         if ( action == free_clean ) then
            nullify ( this%dof_import )
            nullify ( this%p_env )
         end if
         return
      end if
    end if

    ! Free local part
    call this%serial_scalar_array%free_in_stages(action)
    if ( associated ( this%p_env ) ) then
      if ( action == free_clean ) then
         nullify ( this%dof_import )
         nullify ( this%p_env )
      end if
    end if
  end subroutine par_scalar_array_free_in_stages
  
  !=============================================================================
  function par_scalar_array_same_vector_space(this,vector)
   implicit none
   class(par_scalar_array_t), intent(in) :: this
   class(vector_t)             , intent(in) :: vector
   logical :: par_scalar_array_same_vector_space
   par_scalar_array_same_vector_space = .false.
   select type(vector)
   class is (par_scalar_array_t)
     par_scalar_array_same_vector_space = (associated(this%p_env,vector%p_env)) 
     if(this%p_env%am_i_l1_task()) then
        par_scalar_array_same_vector_space = par_scalar_array_same_vector_space .and. (associated(this%dof_import,vector%dof_import))
     end if   
   end select
  end function par_scalar_array_same_vector_space
  
  !=============================================================================
  function par_scalar_array_get_num_blocks(this) result(res)
   implicit none 
   class(par_scalar_array_t), intent(in)   :: this
   integer(ip) :: res
   res = 1
  end function par_scalar_array_get_num_blocks
 
  !=============================================================================
  subroutine par_scalar_array_extract_subvector( this, &
                                               & iblock, &
                                               & size_indices, &
                                               & indices, &
                                               & values )
   implicit none
   class(par_scalar_array_t), intent(in)    :: this 
   integer(ip)              , intent(in)    :: iblock
   integer(ip)              , intent(in)    :: size_indices
   integer(ip)              , intent(in)    :: indices(size_indices)
   real(rp)                 , intent(inout) :: values(*)
   if(this%p_env%am_i_l1_task()) then
      call this%serial_scalar_array%extract_subvector(iblock, &
                                                    & size_indices, &
                                                    & indices, &
                                                    & values )
   end if
  end subroutine par_scalar_array_extract_subvector

  !=============================================================================
  subroutine par_scalar_array_insert_subvector( this, &
                                              & iblock, &
                                              & size_indices, &
                                              & indices, &
                                              & values )
   implicit none
   class(par_scalar_array_t), intent(inout) :: this 
   integer(ip)              , intent(in)    :: iblock
   integer(ip)              , intent(in)    :: size_indices
   integer(ip)              , intent(in)    :: indices(size_indices)
   real(rp)                 , intent(in)    :: values(*)
   if(this%p_env%am_i_l1_task()) then
      call this%serial_scalar_array%insert_subvector(iblock, &
                                                    & size_indices, &
                                                    & indices, &
                                                    & values )
   end if
   end subroutine par_scalar_array_insert_subvector
  

end module par_scalar_array_names
