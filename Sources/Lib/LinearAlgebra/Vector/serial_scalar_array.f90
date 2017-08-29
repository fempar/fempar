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
module serial_scalar_array_names
  use types_names
  use memor_names
#ifdef ENABLE_BLAS
  use blas77_interfaces_names
#endif
  use array_names
  use vector_names

  implicit none
# include "debug.i90"
  
  integer(ip), parameter :: not_created   = 0 ! Initial state
  integer(ip), parameter :: created       = 1 ! size/state/is_a_view already set
  integer(ip), parameter :: entries_ready = 2 ! entries ready

  private
  
  ! State transition diagram for type(serial_scalar_array_t)
  ! -------------------------------------------------
  ! Input State | Action               | Output State 
  ! -------------------------------------------------
  ! not_created | free_numerical_setup | not_created
  ! not_created | free_clean           | not_created
  ! not_created | free                 | not_created
  ! not_created | create               | created
  ! not_created | create_and_allocate  | entries_ready
  ! not_created | clone                | same status as cloned source
  ! not_created | create_view          | entries_ready
  
  ! created     | free_clean           | not_created
  ! created     | free_numerical_setup | created 
  ! created     | allocate             | entries_ready
  ! created     | free                 | not_created
  ! created     | set_view_entries     | entries_ready
  ! created     | clone                | same status as cloned source
  
  ! entries_ready | clone                | same status as cloned source 
  ! entries_ready | free_numerical_setup | created
  ! entries_ready | free                 | not_created
  type, extends(array_t) :: serial_scalar_array_t
     private
     integer(ip)                :: size  = 0
     integer(ip)                :: state = not_created 
     logical                    :: is_a_view  = .false.
     real(rp), pointer          :: b(:) => NULL()
   contains
     procedure :: create_and_allocate    => serial_scalar_array_create_and_allocate
     procedure :: create                 => serial_scalar_array_create
     procedure :: allocate               => serial_scalar_array_allocate
     procedure :: create_view            => serial_scalar_array_create_view
     procedure :: set_view_entries       => serial_scalar_array_set_view_entries
     procedure :: print                  => serial_scalar_array_print
     procedure :: print_matrix_market    => serial_scalar_array_print_matrix_market
     procedure :: get_size               => serial_scalar_array_get_size
     procedure :: get_entries            => serial_scalar_array_get_entries
     
     procedure, private :: serial_scalar_array_insert_single_entry
     procedure, private :: serial_scalar_array_insert_multiple_entries
     procedure, private :: serial_scalar_array_add_single_entry
     procedure, private :: serial_scalar_array_add_multiple_entries
     
     generic  :: insert                  => serial_scalar_array_insert_single_entry, &
                                            serial_scalar_array_insert_multiple_entries
     generic  :: add                     => serial_scalar_array_add_single_entry, &
                                            serial_scalar_array_add_multiple_entries

     procedure :: dot                    => serial_scalar_array_dot
     procedure :: local_dot              => serial_scalar_array_dot
     procedure :: copy                   => serial_scalar_array_copy
     procedure :: init                   => serial_scalar_array_init
     procedure :: scal                   => serial_scalar_array_scal
     procedure :: axpby                  => serial_scalar_array_axpby
     procedure :: nrm2                   => serial_scalar_array_nrm2
     procedure :: clone                  => serial_scalar_array_clone
     procedure :: same_vector_space      => serial_scalar_array_same_vector_space
     procedure :: free_in_stages         => serial_scalar_array_free_in_stages
     procedure :: default_initialization => serial_scalar_array_default_init
     procedure :: get_num_blocks      => serial_scalar_array_get_num_blocks
     procedure :: extract_subvector      => serial_scalar_array_extract_subvector
     procedure :: insert_subvector       => serial_scalar_array_insert_subvector
  end type serial_scalar_array_t

  ! Types
  public :: serial_scalar_array_t

contains

  !=============================================================================
  subroutine serial_scalar_array_default_init (this)
    implicit none
    class(serial_scalar_array_t), intent(inout) :: this
    this%size  = 0
    this%state = not_created 
    this%is_a_view  = .false.
    nullify(this%b)
    call this%NullifyTemporary()
  end subroutine serial_scalar_array_default_init

  !=============================================================================
  subroutine serial_scalar_array_create_and_allocate(this,size)
    implicit none
    class(serial_scalar_array_t), intent(inout) :: this
    integer(ip)                 , intent(in)    :: size
    call this%create(size)
    call this%allocate()
  end subroutine serial_scalar_array_create_and_allocate

  !=============================================================================
  subroutine serial_scalar_array_create(this,size)
    implicit none
    class(serial_scalar_array_t), intent(inout) :: this
    integer(ip)                 , intent(in)    :: size
    assert ( this%state == not_created )
    this%size = size  
    this%state = created
  end subroutine serial_scalar_array_create

  !=============================================================================
  subroutine serial_scalar_array_allocate(this)
    implicit none
    class(serial_scalar_array_t), intent(inout) :: this
    assert ( this%state == created )
    call memallocp(this%size,this%b,__FILE__,__LINE__)
    this%state = entries_ready
    this%is_a_view = .false.
  end subroutine serial_scalar_array_allocate

  !=============================================================================
  subroutine serial_scalar_array_create_view (this, start, end, tvec)
    implicit none
    class(serial_scalar_array_t), intent(in), target  :: this
    integer(ip)     , intent(in)                      :: start
    integer(ip)     , intent(in)                      :: end
    type(serial_scalar_array_t), intent(inout)        :: tvec
    
    assert ( this%state == entries_ready )
    assert ( tvec%state == not_created )
    tvec%size =  end-start+1 ! Number of equations
    tvec%b => this%b(start:end)
    tvec%is_a_view = .true.
    tvec%state = entries_ready
  end subroutine serial_scalar_array_create_view
  
  !=============================================================================
  subroutine serial_scalar_array_set_view_entries ( this, entries ) 
    implicit none
    class(serial_scalar_array_t), intent(inout) ::  this
    real(rp)            , target, intent(in)    ::  entries(:)
    assert ( this%state == created .or. (this%state == entries_ready .and. this%is_a_view) )
    this%b         => entries
    this%is_a_view = .true.
    this%state     = entries_ready
  end subroutine serial_scalar_array_set_view_entries
  
  !=============================================================================
  subroutine serial_scalar_array_print (this, luout)
    implicit none
    class(serial_scalar_array_t), intent(in) :: this
    integer(ip)                , intent(in) :: luout
    assert ( this%state == entries_ready )
    write (luout, '(a)')     '*** begin serial_scalar_array data structure ***'
    write(luout,'(a,i10)') 'size', this%size
    write (luout,'(e25.16)') this%b
    write (luout, '(a)')     '*** end serial_scalar_array data structure ***'
  end subroutine serial_scalar_array_print

  !=============================================================================
  subroutine serial_scalar_array_print_matrix_market ( this, luout )
    implicit none
    class(serial_scalar_array_t), intent(in) :: this
    integer(ip)                , intent(in) :: luout
    integer (ip) :: i
    assert ( this%state == entries_ready )
    write (luout,'(a)') '%%MatrixMarket matrix array real general'
    write (luout,*) this%size , 1
    do i=1,this%size 
       write (luout,'(e32.25)') this%b( i )
    end do
  end subroutine serial_scalar_array_print_matrix_market
  
  subroutine serial_scalar_array_insert_single_entry (this, i, val)
    implicit none
    class(serial_scalar_array_t), intent(inout) :: this
    integer(ip)                 , intent(in)    :: i
    real(rp)                    , intent(in)    :: val
    assert ( i<=0 .or. (i>=1 .and. i<= this%size) )
    if ( i > 0 ) this%b(i) = val 
  end subroutine serial_scalar_array_insert_single_entry
  
  subroutine serial_scalar_array_insert_multiple_entries (this, num_entries, ia, ioffset, val)
    implicit none
    class(serial_scalar_array_t), intent(inout) :: this
    integer(ip)                 , intent(in)    :: num_entries
    integer(ip)                 , intent(in)    :: ia(num_entries)
    integer(ip)                 , intent(in)    :: ioffset
    real(rp)                    , intent(in)    :: val(:)
    integer(ip) :: i, j
    
    do i=1, num_entries
      j = ia(i) 
      assert ( j <= 0 .or. (j >=1 .and. j <= this%size) )
      if (j > 0) this%b(j) = val(i+ioffset) 
    end do
  end subroutine serial_scalar_array_insert_multiple_entries
  
  subroutine serial_scalar_array_add_single_entry (this, i, val)
    implicit none
    class(serial_scalar_array_t), intent(inout) :: this
    integer(ip)                 , intent(in)    :: i
    real(rp)                    , intent(in)    :: val
    assert ( i<=0 .or. (i>=1 .and. i<= this%size) )
    if ( i > 0 ) this%b(i) = this%b(i) + val
  end subroutine serial_scalar_array_add_single_entry
  
    subroutine serial_scalar_array_add_multiple_entries (this, num_entries, ia, ioffset, val)
    implicit none
    class(serial_scalar_array_t), intent(inout) :: this
    integer(ip)                 , intent(in)    :: num_entries
    integer(ip)                 , intent(in)    :: ia(num_entries)
    integer(ip)                 , intent(in)    :: ioffset
    real(rp)                    , intent(in)    :: val(:)
    integer(ip) :: i, j
    
    do i=1, num_entries
      j = ia(i) 
      assert ( j <= 0 .or. (j >=1 .and. j <= this%size) )
      if (j > 0) this%b(j) = this%b(j) + val(i+ioffset) 
    end do
    
  end subroutine serial_scalar_array_add_multiple_entries
  
  function serial_scalar_array_get_size ( this )
    implicit none
    class(serial_scalar_array_t), intent(in) :: this
    integer(ip) :: serial_scalar_array_get_size
    assert ( this%state == entries_ready )
    serial_scalar_array_get_size = this%size
  end function serial_scalar_array_get_size
  
  !=============================================================================
  function serial_scalar_array_get_entries ( this )
    implicit none
    class(serial_scalar_array_t), target, intent(in) :: this
    real(rp)                    , pointer            :: serial_scalar_array_get_entries(:)
    assert ( this%state == entries_ready )
    serial_scalar_array_get_entries => this%b
  end function serial_scalar_array_get_entries

  ! alpha <- op1^T * op2
  function serial_scalar_array_dot(op1,op2) result(alpha)
    implicit none
    class(serial_scalar_array_t), intent(in)    :: op1
    class(vector_t), intent(in)  :: op2
    real(rp) :: alpha
    assert ( op1%state == entries_ready )
    
    call op1%GuardTemp()
    call op2%GuardTemp()
    select type(op2)
       class is (serial_scalar_array_t)
       assert ( op2%state == entries_ready )
       assert ( op1%size == op2%size )
#ifdef ENABLE_BLAS
       alpha = ddot( op1%size, op1%b, 1, op2%b, 1 )
#else
       check(1==0)
#endif
       class default
       write(0,'(a)') 'serial_scalar_array_t%dot: unsupported op2 class'
       check(1==0)
    end select
    call op1%CleanTemp()
    call op2%CleanTemp()
  end function serial_scalar_array_dot

  ! op1 <- op2 
  subroutine serial_scalar_array_copy(op1,op2)
    implicit none
    class(serial_scalar_array_t), intent(inout) :: op1
    class(vector_t), intent(in)  :: op2
    assert ( op1%state == entries_ready )
    call op2%GuardTemp()
    select type(op2)
       class is (serial_scalar_array_t)
       assert ( op2%state == entries_ready )
       assert ( op2%size == op1%size )
#ifdef ENABLE_BLAS
       call dcopy ( op2%size, op2%b, 1, op1%b, 1 ) 
#else
       op1%b=op2%b
#endif
       class default
       write(0,'(a)') 'serial_scalar_array_t%copy: unsupported op2 class'
       check(1==0)
    end select
    call op2%CleanTemp()
  end subroutine serial_scalar_array_copy

  ! op1 <- alpha * op2
  subroutine serial_scalar_array_scal(op1,alpha,op2)
    implicit none
    class(serial_scalar_array_t), intent(inout) :: op1
    real(rp), intent(in) :: alpha
    class(vector_t), intent(in) :: op2
    assert ( op1%state == entries_ready )
    
    call op2%GuardTemp()
    select type(op2)
       class is (serial_scalar_array_t)
       assert ( op2%state == entries_ready )
       assert ( op2%size == op1%size )
#ifdef ENABLE_BLAS
       call dcopy ( op2%size, op2%b, 1, op1%b, 1)
       call dscal ( op1%size, alpha, op1%b, 1)
#else
       op1%b=alpha*op2%b
#endif
       class default
       write(0,'(a)') 'serial_scalar_array_t%scal: unsupported op2 class'
       check(1==0)
    end select
    call op2%CleanTemp()
  end subroutine serial_scalar_array_scal
  ! op <- alpha
  subroutine serial_scalar_array_init(op,alpha)
    implicit none
    class(serial_scalar_array_t), intent(inout) :: op
    real(rp), intent(in) :: alpha
    assert ( op%state == entries_ready )
    op%b=alpha
  end subroutine serial_scalar_array_init

  ! op1 <- alpha*op2 + beta*op1
  subroutine serial_scalar_array_axpby(op1, alpha, op2, beta)
    implicit none
    class(serial_scalar_array_t), intent(inout) :: op1
    real(rp), intent(in) :: alpha
    class(vector_t), intent(in) :: op2
    real(rp), intent(in) :: beta

    assert ( op1%state == entries_ready )
    call op2%GuardTemp()
    select type(op2)
       class is (serial_scalar_array_t)
       assert ( op2%size == op1%size )
       assert ( op2%state == entries_ready )
       if ( beta == 0.0_rp ) then
          call op1%scal(alpha, op2)
       else if ( beta == 1.0_rp ) then
          ! AXPY
#ifdef ENABLE_BLAS
          call daxpy ( op2%size, alpha, op2%b, 1, op1%b, 1 )
#else
          op1%b=op1%b+alpha*op2%b
#endif
       else
          ! SCAL + AXPY
          call op1%scal(beta, op1)
#ifdef ENABLE_BLAS
          call daxpy ( op2%size, alpha, op2%b, 1, op1%b, 1 )    
#else
          op1%b=op1%b+alpha*op2%b
#endif  
       end if
       class default
       write(0,'(a)') 'serial_scalar_array_t%axpby: unsupported op2 class'
       check(1==0)
    end select
    call op2%CleanTemp()
  end subroutine serial_scalar_array_axpby

  ! alpha <- nrm2(op)
  function serial_scalar_array_nrm2(op) result(alpha)
    implicit none
    class(serial_scalar_array_t), intent(in)  :: op
    real(rp) :: alpha
    call op%GuardTemp()
    assert ( op%state == entries_ready )
#ifdef ENABLE_BLAS
    alpha = dnrm2( op%size, op%b, 1 )
#else
    alpha = op%dot(op)
    alpha = sqrt(alpha)
#endif

    call op%CleanTemp()
  end function serial_scalar_array_nrm2

  ! op1 <- clone(op2) 
  subroutine serial_scalar_array_clone(op1,op2)
    implicit none
    class(serial_scalar_array_t), target, intent(inout) :: op1
    class(vector_t), target, intent(in)                :: op2
    class(vector_t), pointer :: p
    p => op1
    if(associated(p,op2)) return ! It's aliasing
    
    call op2%GuardTemp()
    select type(op2)
       class is (serial_scalar_array_t)
       assert ( op2%state == created .or. op2%state == entries_ready )
       call op1%free()  
       call op1%create(op2%size)
       if ( op2%state == entries_ready) then
          call op1%allocate()
       end if
       class default
       write(0,'(a)') 'serial_scalar_array_t%clone: unsupported op2 class'
       check(.false.)
    end select
    call op2%CleanTemp()
  end subroutine serial_scalar_array_clone

  !=============================================================================
  subroutine serial_scalar_array_free_in_stages(this,action)
    implicit none
    class(serial_scalar_array_t), intent(inout) :: this
    integer(ip)                 , intent(in)    :: action
    assert ( action == free_clean .or. action == free_symbolic_setup .or. action == free_numerical_setup )
    
    if ( action == free_clean ) then
       ! Undo create
       assert ( .not. this%state == entries_ready )
       if ( this%state == created ) then
         this%size = 0
         this%state = not_created
         this%is_a_view = .false.
       end if
    else if ( action == free_numerical_setup ) then
       ! Undo allocate
       if ( this%state == entries_ready ) then
          if ( this%is_a_view ) then
            this%is_a_view = .false.
            nullify(this%b)
          else 
            call memfreep(this%b,__FILE__,__LINE__)  
            nullify(this%b)
          end if
          this%state = created
       end if
    end if
  end subroutine serial_scalar_array_free_in_stages
  
  !=============================================================================
  function serial_scalar_array_same_vector_space(this,vector)
    implicit none
    class(serial_scalar_array_t), intent(in) :: this
    class(vector_t)             , intent(in) :: vector
    logical :: serial_scalar_array_same_vector_space
    serial_scalar_array_same_vector_space = .false.
    assert (this%state == created .or. this%state == entries_ready)
    select type(vector)
    class is (serial_scalar_array_t)
      serial_scalar_array_same_vector_space = (this%size == vector%size)
    end select
  end function serial_scalar_array_same_vector_space
	
  !=============================================================================
  function serial_scalar_array_get_num_blocks(this) result(res)
    implicit none 
    class(serial_scalar_array_t), intent(in)   :: this
    integer(ip) :: res
    res = 1
  end function serial_scalar_array_get_num_blocks
 
  !=============================================================================
  subroutine serial_scalar_array_extract_subvector( this, &
                                                  & iblock, &
                                                  & size_indices, &
                                                  & indices, &
                                                  & values )
    implicit none
    class(serial_scalar_array_t), intent(in)    :: this 
    integer(ip)                 , intent(in)    :: iblock
    integer(ip)                 , intent(in)    :: size_indices
    integer(ip)                 , intent(in)    :: indices(size_indices)
    real(rp)                    , intent(inout) :: values(*)
    integer(ip)                                 :: i

    assert(this%state == entries_ready)
    assert(iblock == 1)
        
    do i=1, size_indices
      assert ( indices(i) <= this%size )
      if ( indices(i) > 0 ) then
        values(i) = this%b(indices(i))
      end if
    end do 
  end subroutine serial_scalar_array_extract_subvector
 
  !=============================================================================
  subroutine serial_scalar_array_insert_subvector( this, &
                                                 & iblock, &
                                                 & size_indices, &
                                                 & indices, &
                                                 & values )
    implicit none
    class(serial_scalar_array_t), intent(inout) :: this 
    integer(ip)                 , intent(in)    :: iblock
    integer(ip)                 , intent(in)    :: size_indices
    integer(ip)                 , intent(in)    :: indices(size_indices)
    real(rp)                    , intent(in)    :: values(*)
    integer(ip)                                 :: i

    assert(this%state == entries_ready)
    assert(iblock == 1)
        
    do i=1, size_indices
      assert ( indices(i) <= this%size )
      if ( indices(i) > 0 ) then
        this%b(indices(i)) = values(i)
      end if
    end do 
  end subroutine serial_scalar_array_insert_subvector

end module serial_scalar_array_names
