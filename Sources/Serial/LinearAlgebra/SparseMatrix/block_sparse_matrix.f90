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
module block_sparse_matrix_names
  use types_names
  use memor_names
  use sparse_matrix_names
  use serial_block_array_names
  use serial_scalar_array_names
  
  ! Abstract types
  use vector_names
  use operator_names
  use matrix_names
  use vector_space_names

  implicit none
# include "debug.i90"

  private
  
  type p_sparse_matrix_t
    type(sparse_matrix_t), pointer :: sparse_matrix
  end type p_sparse_matrix_t

  type, extends(matrix_t):: block_sparse_matrix_t
     private 
     integer(ip) :: nblocks = -1
     type(p_sparse_matrix_t), allocatable :: blocks(:,:)
   contains
     procedure, private :: block_sparse_matrix_create_only_blocks_container
     procedure, private :: block_sparse_matrix_create_blocks_container_and_all_blocks
     !generic :: create  => block_sparse_matrix_create_only_blocks_container, &
     !                      block_sparse_matrix_create_blocks_container_and_all_blocks
     generic :: create => block_sparse_matrix_create_blocks_container_and_all_blocks
     
     ! This set of methods should be re-thought/developed in the future whenever we have a clearer
     ! idea of which services these objects should provide and how they should behave
     ! (i.e., once we know how to precisely define their state transition diagram)
     ! procedure, private :: block_sparse_matrix_create_diagonal_block
     ! procedure, private :: block_sparse_matrix_create_offdiagonal_block
     ! generic :: create_block => block_sparse_matrix_create_diagonal_block, &
     !                            block_sparse_matrix_create_offdiagonal_block
     procedure :: set_block_to_zero => block_sparse_matrix_set_block_to_zero                     
     
					procedure :: compress_storage              => block_sparse_matrix_compress_storage
     procedure :: allocate                      => block_sparse_matrix_allocate
     procedure :: free_in_stages                => block_sparse_matrix_free_in_stages
     procedure :: get_block                     => block_sparse_matrix_get_block
     procedure :: get_nblocks                   => block_sparse_matrix_get_nblocks
     procedure :: apply                         => block_sparse_matrix_apply
     procedure, private :: create_vector_spaces => block_sparse_matrix_create_vector_spaces
     procedure, non_overridable :: get_iterator => block_sparse_matrix_get_iterator
  end type block_sparse_matrix_t

  type, extends(matrix_iterator_t) :: block_sparse_matrix_iterator_t
      private
      integer(ip) :: i_index
      integer(ip) :: j_index
      type(sparse_matrix_iterator_t) :: sparse_iterator
      type(block_sparse_matrix_t), pointer :: block_sparse_matrix
    contains
       procedure, non_overridable :: free         => block_sparse_matrix_iterator_free
       procedure, non_overridable :: next         => block_sparse_matrix_iterator_next
       procedure, non_overridable :: has_finished => block_sparse_matrix_iterator_has_finished
       procedure, non_overridable :: get_row      => block_sparse_matrix_iterator_get_row
       procedure, non_overridable :: get_column   => block_sparse_matrix_iterator_get_column
       procedure, non_overridable :: get_value    => block_sparse_matrix_iterator_get_value
       procedure, non_overridable :: set_value    => block_sparse_matrix_iterator_set_value
       procedure, non_overridable :: get_iblock   => block_sparse_matrix_iterator_get_iblock
       procedure, non_overridable :: get_jblock   => block_sparse_matrix_iterator_get_jblock
   end type block_sparse_matrix_iterator_t

  ! Types
  public :: block_sparse_matrix_t
  public :: block_sparse_matrix_iterator_t

contains

  !=============================================================================
  subroutine block_sparse_matrix_create_blocks_container_and_all_blocks(this, & 
                                                                        nblocks, & 
                                                                        diagonal_blocks_num_rows, &
                                                                        diagonal_blocks_num_cols, &
                                                                        diagonal_blocks_symmetric_storage,&
                                                                        diagonal_blocks_symmetric,& 
                                                                        diagonal_blocks_sign)
    implicit none
    class(block_sparse_matrix_t), intent(inout) :: this
    integer(ip)                 , intent(in)    :: nblocks
    integer(ip)                 , intent(in)    :: diagonal_blocks_num_rows(nblocks)
    integer(ip)                 , intent(in)    :: diagonal_blocks_num_cols(nblocks)
    logical                     , intent(in)    :: diagonal_blocks_symmetric_storage(nblocks)
    logical                     , intent(in)    :: diagonal_blocks_symmetric(nblocks)
    integer(ip)                 , intent(in)    :: diagonal_blocks_sign(nblocks)

    integer(ip) :: ib,jb,istat

    call this%block_sparse_matrix_create_only_blocks_container(nblocks)
    do ib=1, this%nblocks 
       do jb=1, this%nblocks
          allocate ( this%blocks(ib,jb)%sparse_matrix )
          if ( ib == jb ) then
             assert ( diagonal_blocks_num_rows(ib) == diagonal_blocks_num_cols(ib) )
             call this%blocks(ib,jb)%sparse_matrix%create ( diagonal_blocks_num_rows(ib), &
                                                            diagonal_blocks_symmetric_storage(ib), &
                                                            diagonal_blocks_symmetric(ib), &
                                                            diagonal_blocks_sign(ib) )
          else
             call this%blocks(ib,jb)%sparse_matrix%create( diagonal_blocks_num_rows(ib), &
                                                           diagonal_blocks_num_cols(jb) )
          end if
       end do
    end do
    call this%create_vector_spaces()
  end subroutine block_sparse_matrix_create_blocks_container_and_all_blocks

  !=============================================================================
  subroutine block_sparse_matrix_create_only_blocks_container(this,nblocks)
    implicit none
    ! Parameters
    class(block_sparse_matrix_t), intent(inout) :: this
    integer(ip)                 , intent(in)    :: nblocks
    integer(ip) :: istat
    this%nblocks = nblocks
    allocate ( this%blocks(this%nblocks,this%nblocks), stat=istat )
    check (istat==0)
  end subroutine block_sparse_matrix_create_only_blocks_container
  
  !=============================================================================
  subroutine block_sparse_matrix_compress_storage(this, sparse_matrix_storage_format )
    implicit none
    class(block_sparse_matrix_t), intent(inout) :: this
				character(*)                , intent(in)    :: sparse_matrix_storage_format

    integer(ip) :: ib,jb
    do ib=1, this%nblocks
       do jb=1, this%nblocks
          if ( associated( this%blocks(ib,jb)%sparse_matrix ) ) then
            call this%blocks(ib,jb)%sparse_matrix%convert(sparse_matrix_storage_format)
          end if
       end do
    end do
  end subroutine block_sparse_matrix_compress_storage
		
  !=============================================================================
  subroutine block_sparse_matrix_allocate(this)
    implicit none
    class(block_sparse_matrix_t), intent(inout) :: this
    integer(ip) :: ib,jb
    do ib=1, this%nblocks
       do jb=1, this%nblocks
          if ( associated( this%blocks(ib,jb)%sparse_matrix ) ) then
            call this%blocks(ib,jb)%sparse_matrix%allocate()
          end if
       end do
    end do
  end subroutine block_sparse_matrix_allocate
  
  subroutine block_sparse_matrix_create_vector_spaces(this)
    implicit none
    class(block_sparse_matrix_t), intent(inout) :: this
    integer(ip) :: ib
    
    integer(ip), allocatable :: size_of_blocks_domain(:)
    integer(ip), allocatable :: size_of_blocks_range(:)
    type(serial_block_array_t) :: range_vector
    type(serial_block_array_t) :: domain_vector
    class(vector_space_t), pointer :: range_vector_space
    class(vector_space_t), pointer :: domain_vector_space

    call memalloc(this%nblocks,size_of_blocks_domain,__FILE__,__LINE__)
    call memalloc(this%nblocks,size_of_blocks_range,__FILE__,__LINE__)
    
    do ib=1, this%nblocks
        size_of_blocks_domain(ib) = this%blocks(ib,ib)%sparse_matrix%get_num_cols()
        size_of_blocks_range(ib) = this%blocks(ib,ib)%sparse_matrix%get_num_rows()
    end do
    
    ! Create and set domain and range vector spaces
    call range_vector%create(this%nblocks, size_of_blocks_range)
    call domain_vector%create(this%nblocks, size_of_blocks_domain)
    
    range_vector_space => this%get_range_vector_space()
    call range_vector_space%create(range_vector)
    
    domain_vector_space => this%get_domain_vector_space()
    call domain_vector_space%create(domain_vector)
    
    call range_vector%free()
    call domain_vector%free()
    
    call memfree(size_of_blocks_domain,__FILE__,__LINE__)
    call memfree(size_of_blocks_range,__FILE__,__LINE__)
  end subroutine block_sparse_matrix_create_vector_spaces
  

  function block_sparse_matrix_get_block (bmat,ib,jb)
    implicit none
    ! Parameters
    class(block_sparse_matrix_t), target, intent(in) :: bmat
    integer(ip)                    , intent(in) :: ib,jb
    type(sparse_matrix_t)               , pointer    :: block_sparse_matrix_get_block
    block_sparse_matrix_get_block =>  bmat%blocks(ib,jb)%sparse_matrix
  end function block_sparse_matrix_get_block

  function block_sparse_matrix_get_nblocks (this)
    implicit none
    ! Parameters
    class(block_sparse_matrix_t), target, intent(in) :: this
    integer(ip)                                :: block_sparse_matrix_get_nblocks
    block_sparse_matrix_get_nblocks = this%nblocks
  end function block_sparse_matrix_get_nblocks

  ! op%apply(x,y) <=> y <- op*x
  ! Implicitly assumes that y is already allocated
  subroutine block_sparse_matrix_apply(op,x,y)
    implicit none
    class(block_sparse_matrix_t), intent(in)    :: op
    class(vector_t), intent(in)    :: x
    class(vector_t), intent(inout) :: y
    ! Locals
    integer(ip) :: ib,jb
    type(serial_scalar_array_t)          :: aux
    type(serial_scalar_array_t), pointer :: y_block

    call op%abort_if_not_in_domain(x)
    call op%abort_if_not_in_range(y)
    
    call x%GuardTemp()
    call y%init(0.0_rp)
    select type(x)
       class is (serial_block_array_t)
       select type(y)
          class is(serial_block_array_t)
          do ib=1,op%nblocks
             call aux%clone(y%get_block(ib))
             y_block => y%get_block(ib)
             do jb=1,op%nblocks
                if ( associated(op%blocks(ib,jb)%sparse_matrix) ) then
                   ! aux <- A(ib,jb) * x(jb)
                   call op%blocks(ib,jb)%sparse_matrix%apply(x%get_block(jb),aux)
                   ! y(ib) <- y(ib) + aux
                   call y_block%axpby(1.0_rp,aux,1.0_rp)
                end if
             end do
             call aux%free()
          end do
       end select
    end select
    call x%CleanTemp()
  end subroutine block_sparse_matrix_apply

  subroutine block_sparse_matrix_free_in_stages(this,action)
    implicit none
    class(block_sparse_matrix_t), intent(inout) :: this
    integer(ip)                 , intent(in)    :: action
    integer(ip) :: ib,jb, istat

    do ib=1, this%nblocks 
       do jb=1, this%nblocks
          if ( associated(this%blocks(ib,jb)%sparse_matrix) ) then
             if ( action == free_clean ) then
                call this%blocks(ib,jb)%sparse_matrix%free_in_stages(free_clean)
                deallocate (this%blocks(ib,jb)%sparse_matrix, stat=istat)
                check(istat==0)
             else if ( action == free_symbolic_setup ) then 
                call this%blocks(ib,jb)%sparse_matrix%free_in_stages(free_symbolic_setup)
             else if ( action == free_numerical_setup ) then
                call this%blocks(ib,jb)%sparse_matrix%free_in_stages(free_numerical_setup)
             end if
          end if
       end do
    end do
    if (action == free_clean) then
       this%nblocks = 0
       deallocate (this%blocks, stat=istat)
       check(istat==0)
       call this%free_vector_spaces()
    end if
  end subroutine block_sparse_matrix_free_in_stages
  
  !subroutine block_sparse_matrix_create_diagonal_block (bmat,ib,diagonal_block_num_rows_and_cols,diagonal_block_symmetric_storage,diagonal_block_symmetric,diagonal_block_sign)
  !  implicit none
  !  ! Parameters
  !  class(block_sparse_matrix_t)   , intent(inout) :: bmat
  !  integer(ip)                    , intent(in)    :: ib
  !  integer(ip)                    , intent(in)    :: diagonal_block_num_rows_and_cols
  !  logical                        , intent(in)    :: diagonal_block_symmetric_storage
  !  logical                        , intent(in)    :: diagonal_block_symmetric
  !  integer(ip)                    , intent(in)    :: diagonal_block_sign
  !  integer(ip) :: istat

  !  assert ( .not. associated ( bmat%blocks(ib,ib)%sparse_matrix ) )
  !  allocate ( bmat%blocks(ib,ib)%sparse_matrix, stat=istat )
  !  check ( istat==0 )
  !  call bmat%blocks(ib,ib)%sparse_matrix%create ( diagonal_block_num_rows_and_cols, diagonal_block_symmetric_storage, diagonal_block_symmetric, diagonal_block_sign )
  !end subroutine block_sparse_matrix_create_diagonal_block

  !subroutine block_sparse_matrix_create_offdiagonal_block (bmat,ib,jb,num_rows,num_cols)
  !  implicit none
  !  ! Parameters
  !  class(block_sparse_matrix_t)   , intent(inout) :: bmat
  !  integer(ip)                    , intent(in)    :: ib,jb
  !  integer(ip)                    , intent(in)    :: num_rows, num_cols
  !  integer(ip) :: istat
  !  assert ( ib /= jb ) 
  !  assert ( .not. associated(bmat%blocks(ib,jb)%sparse_matrix ))
  !  allocate ( bmat%blocks(ib,ib)%sparse_matrix, stat=istat )
  !  check ( istat==0 )
  !  call bmat%blocks(ib,jb)%sparse_matrix%create (num_rows,num_cols)
  !end subroutine block_sparse_matrix_create_offdiagonal_block
  
  subroutine block_sparse_matrix_set_block_to_zero (bmat,ib,jb)
   implicit none
   ! Parameters
   class(block_sparse_matrix_t), intent(inout) :: bmat
   integer(ip)                 , intent(in)    :: ib,jb
   integer(ip) :: istat

   if ( associated(bmat%blocks(ib,jb)%sparse_matrix) ) then
      ! Undo create
      call bmat%blocks(ib,jb)%sparse_matrix%free_in_stages(free_clean)
      deallocate (bmat%blocks(ib,jb)%sparse_matrix,stat=istat)
      check(istat==0)
      nullify ( bmat%blocks(ib,jb)%sparse_matrix )
   end if
  end subroutine block_sparse_matrix_set_block_to_zero
  
  subroutine block_sparse_matrix_get_iterator(this, iterator)
    !-----------------------------------------------------------------
    !< Get a pointer to an iterator over the matrix entries
    !-----------------------------------------------------------------
    class(block_sparse_matrix_t), target , intent(in)  :: this
    class(block_sparse_matrix_iterator_t), intent(out) :: iterator
    !-----------------------------------------------------------------
    iterator%i_index = 1
    iterator%j_index = 1
    iterator%block_sparse_matrix => this
    call this%blocks(iterator%i_index,iterator%j_index)%sparse_matrix%get_iterator(iterator%sparse_iterator)
  end subroutine block_sparse_matrix_get_iterator

  !-----------------------------------------------------------------
  !< BLOCK_SPARSE_MATRIX_ITERATOR SUBROUTINES
  !-----------------------------------------------------------------
  ! NOT TESTED!!
  subroutine block_sparse_matrix_iterator_free(this)
    !-----------------------------------------------------------------
    !< Set the pointer to the following entry of the matrix
    !-----------------------------------------------------------------
    class(block_sparse_matrix_iterator_t), intent(inout) :: this
    call this%sparse_iterator%free()
    this%i_index = -1
    this%j_index = -1
    this%block_sparse_matrix => NULL()
  end subroutine block_sparse_matrix_iterator_free

  subroutine block_sparse_matrix_iterator_next(this)
    !-----------------------------------------------------------------
    !< Set the pointer to the following entry of the matrix
    !-----------------------------------------------------------------
    class(block_sparse_matrix_iterator_t), intent(inout) :: this
    call this%sparse_iterator%next()
    if (this%sparse_iterator%has_finished()) then
       this%j_index =  this%j_index + 1
       if ( this%j_index > this%block_sparse_matrix%nblocks) then
          this%j_index = 1
          this%i_index = this%i_index + 1
          if ( this%i_index <= this%block_sparse_matrix%nblocks) then
              call this%block_sparse_matrix%blocks(this%i_index,this%j_index)%sparse_matrix%get_iterator(this%sparse_iterator)
           end if
        end if
     end if
  end subroutine block_sparse_matrix_iterator_next

  function block_sparse_matrix_iterator_has_finished(this)
    !-----------------------------------------------------------------
    !< Check if the pointer of the matrix has reached the end
    !-----------------------------------------------------------------
    class(block_sparse_matrix_iterator_t), intent(in) :: this
    logical :: block_sparse_matrix_iterator_has_finished
    
    block_sparse_matrix_iterator_has_finished = (this%i_index > this%block_sparse_matrix%nblocks)
    
  end function block_sparse_matrix_iterator_has_finished

  function block_sparse_matrix_iterator_get_row(this)
    !-----------------------------------------------------------------
    !< Get the row index of the entry of the matrix
    !-----------------------------------------------------------------
    class(block_sparse_matrix_iterator_t), intent(in) :: this
    integer(ip) :: block_sparse_matrix_iterator_get_row
    block_sparse_matrix_iterator_get_row = this%sparse_iterator%get_row()
  end function block_sparse_matrix_iterator_get_row

  function block_sparse_matrix_iterator_get_column(this)
    !-----------------------------------------------------------------
    !<  Get the column index of the entry of the matrix
    !-----------------------------------------------------------------
    class(block_sparse_matrix_iterator_t), intent(in) :: this
    integer(ip) :: block_sparse_matrix_iterator_get_column
    block_sparse_matrix_iterator_get_column = this%sparse_iterator%get_column()
  end function block_sparse_matrix_iterator_get_column

  function block_sparse_matrix_iterator_get_value(this)
    !-----------------------------------------------------------------
    !< Get the value of the entry of the matrix
    !-----------------------------------------------------------------
    class(block_sparse_matrix_iterator_t), intent(in) :: this
    real(rp) :: block_sparse_matrix_iterator_get_value
    block_sparse_matrix_iterator_get_value = this%sparse_iterator%get_value()
  end function block_sparse_matrix_iterator_get_value

  subroutine block_sparse_matrix_iterator_set_value(this,new_value)
    !-----------------------------------------------------------------
    !< Set the value of the entry of the matrix
    !-----------------------------------------------------------------
    class(block_sparse_matrix_iterator_t), intent(inout) :: this
    real(rp)                            , intent(in)    :: new_value
    call this%sparse_iterator%set_value(new_value)
  end subroutine block_sparse_matrix_iterator_set_value

  function block_sparse_matrix_iterator_get_iblock(this)
    !-----------------------------------------------------------------
    !< Get the i-block index
    !-----------------------------------------------------------------
    class(block_sparse_matrix_iterator_t), intent(in) :: this
    integer(ip) :: block_sparse_matrix_iterator_get_iblock
    block_sparse_matrix_iterator_get_iblock = this%i_index
  end function block_sparse_matrix_iterator_get_iblock

  function block_sparse_matrix_iterator_get_jblock(this)
    !-----------------------------------------------------------------
    !< Get the j-block index
    !-----------------------------------------------------------------
    class(block_sparse_matrix_iterator_t), intent(in) :: this
    integer(ip) :: block_sparse_matrix_iterator_get_jblock
    block_sparse_matrix_iterator_get_jblock = this%j_index
  end function block_sparse_matrix_iterator_get_jblock


end module block_sparse_matrix_names
