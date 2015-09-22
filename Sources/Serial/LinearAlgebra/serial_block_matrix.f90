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
module serial_block_matrix_names
  use types_names
  use memor_names
  use serial_scalar_matrix_names
  use serial_block_array_names
  use serial_scalar_array_names
  
  ! Abstract types
  use vector_names
  use operator_names
  use matrix_names

  implicit none
# include "debug.i90"

  private
  
  type p_serial_scalar_matrix_t
    type(serial_scalar_matrix_t), pointer :: serial_scalar_matrix
  end type p_serial_scalar_matrix_t

  type, extends(matrix_t):: serial_block_matrix_t
    ! private ! IBM XLF 14.1 bug
    integer(ip) :: nblocks = -1
    type(p_serial_scalar_matrix_t), allocatable :: blocks(:,:)
  contains
    procedure, private :: serial_block_matrix_create_only_blocks_container
    procedure, private :: serial_block_matrix_create_blocks_container_and_all_blocks
    generic :: create  => serial_block_matrix_create_only_blocks_container, &
	                      serial_block_matrix_create_blocks_container_and_all_blocks
	
    procedure, private :: serial_block_matrix_create_diagonal_block
	procedure, private :: serial_block_matrix_create_offdiagonal_block
    generic :: create_block => serial_block_matrix_create_diagonal_block, &
	                           serial_block_matrix_create_offdiagonal_block
							 
	procedure :: allocate           => serial_block_matrix_allocate
	procedure :: free_in_stages     => serial_block_matrix_free_in_stages
    procedure :: set_block_to_zero  => serial_block_matrix_set_block_to_zero
    procedure :: get_block          => serial_block_matrix_get_block
    procedure :: get_nblocks        => serial_block_matrix_get_nblocks
    procedure :: apply              => serial_block_matrix_apply
    procedure :: apply_fun          => serial_block_matrix_apply_fun
  end type serial_block_matrix_t

  ! Types
  public :: serial_block_matrix_t

contains

  !=============================================================================
  subroutine serial_block_matrix_create_blocks_container_and_all_blocks(bmat,nblocks,diagonal_blocks_symmetric_storage,diagonal_blocks_symmetric,diagonal_blocks_sign)
    implicit none
    class(serial_block_matrix_t), intent(out) :: bmat
	integer(ip)                 , intent(in)  :: nblocks
	logical                     , intent(in)  :: diagonal_blocks_symmetric_storage(nblocks)
    logical                     , intent(in)  :: diagonal_blocks_symmetric(nblocks)
    integer(ip)                 , intent(in)  :: diagonal_blocks_sign(nblocks)

    integer(ip) :: ib,jb,istat

    bmat%nblocks = nblocks
    allocate ( bmat%blocks(bmat%nblocks,bmat%nblocks), stat=istat )
	check (istat == 0)
    do ib=1, bmat%nblocks 
       do jb=1, bmat%nblocks
         allocate ( bmat%blocks(ib,jb)%serial_scalar_matrix )
         if ( ib == jb ) then
           call bmat%blocks(ib,jb)%serial_scalar_matrix%create (diagonal_blocks_symmetric_storage(ib), diagonal_blocks_symmetric(ib), diagonal_blocks_sign(ib) )
         else
           call bmat%blocks(ib,jb)%serial_scalar_matrix%create()
         end if
       end do
    end do
  end subroutine serial_block_matrix_create_blocks_container_and_all_blocks
  
  !=============================================================================
  subroutine serial_block_matrix_create_only_blocks_container(bmat,nblocks)
    implicit none
    ! Parameters
    class(serial_block_matrix_t), intent(out) :: bmat
	integer(ip)                 , intent(in)  :: nblocks
	integer(ip) :: istat
    bmat%nblocks = nblocks
    allocate ( bmat%blocks(bmat%nblocks,bmat%nblocks), stat=istat )
	check (istat==0)
  end subroutine serial_block_matrix_create_only_blocks_container

  subroutine serial_block_matrix_create_diagonal_block (bmat,ib,diagonal_block_symmetric_storage,diagonal_block_symmetric,diagonal_block_sign)
    implicit none
    ! Parameters
    class(serial_block_matrix_t)   , intent(inout) :: bmat
    integer(ip)                    , intent(in)    :: ib
    logical                        , intent(in)    :: diagonal_block_symmetric_storage
	logical                        , intent(in)    :: diagonal_block_symmetric
    integer(ip)                    , intent(in)    :: diagonal_block_sign
    integer(ip) :: istat

    assert ( associated ( bmat%blocks(ib,ib)%serial_scalar_matrix ) )
    allocate ( bmat%blocks(ib,ib)%serial_scalar_matrix, stat=istat )
	check ( istat==0 )
    call bmat%blocks(ib,ib)%serial_scalar_matrix%create ( diagonal_block_symmetric_storage, diagonal_block_symmetric, diagonal_block_sign )
  end subroutine serial_block_matrix_create_diagonal_block
  
  subroutine serial_block_matrix_create_offdiagonal_block (bmat,ib,jb)
    implicit none
    ! Parameters
    class(serial_block_matrix_t)   , intent(inout) :: bmat
    integer(ip)                    , intent(in)    :: ib,jb
	
    assert ( ib /= jb ) 
    assert ( associated(bmat%blocks(ib,jb)%serial_scalar_matrix ))
    call bmat%blocks(ib,jb)%serial_scalar_matrix%create ()
  end subroutine serial_block_matrix_create_offdiagonal_block
  
  !=============================================================================
  subroutine serial_block_matrix_allocate(this)
    implicit none
    class(serial_block_matrix_t), intent(inout) :: this
    integer(ip) :: ib,jb
	
    do ib=1, this%nblocks
       do jb=1, this%nblocks
         if ( associated( this%blocks(ib,jb)%serial_scalar_matrix ) ) then
           call this%blocks(ib,jb)%serial_scalar_matrix%allocate()
		 end if
       end do
    end do
  end subroutine serial_block_matrix_allocate

  subroutine serial_block_matrix_set_block_to_zero (bmat,ib,jb)
    implicit none
    ! Parameters
    class(serial_block_matrix_t), intent(inout) :: bmat
    integer(ip)                 , intent(in)    :: ib,jb
	integer(ip) :: istat

    if ( associated(bmat%blocks(ib,jb)%serial_scalar_matrix) ) then
	   ! Undo create
       call bmat%blocks(ib,jb)%serial_scalar_matrix%free_in_stages(free_clean)
       deallocate (bmat%blocks(ib,jb)%serial_scalar_matrix,stat=istat)
	   check(istat==0)
       nullify ( bmat%blocks(ib,jb)%serial_scalar_matrix )
    end if
  end subroutine serial_block_matrix_set_block_to_zero

  function serial_block_matrix_get_block (bmat,ib,jb)
    implicit none
    ! Parameters
    class(serial_block_matrix_t), target, intent(in) :: bmat
    integer(ip)                    , intent(in) :: ib,jb
    type(serial_scalar_matrix_t)               , pointer    :: serial_block_matrix_get_block

    serial_block_matrix_get_block =>  bmat%blocks(ib,jb)%serial_scalar_matrix
  end function serial_block_matrix_get_block

  function serial_block_matrix_get_nblocks (this)
    implicit none
    ! Parameters
    class(serial_block_matrix_t), target, intent(in) :: this
    integer(ip)                                :: serial_block_matrix_get_nblocks
    serial_block_matrix_get_nblocks = this%nblocks
  end function serial_block_matrix_get_nblocks


  ! op%apply(x,y) <=> y <- op*x
  ! Implicitly assumes that y is already allocated
  subroutine serial_block_matrix_apply(op,x,y)
    implicit none
    class(serial_block_matrix_t), intent(in)    :: op
    class(vector_t), intent(in)    :: x
    class(vector_t), intent(inouT) :: y
    ! Locals
    integer(ip) :: ib,jb
    type(serial_scalar_array_t) :: aux

    call x%GuardTemp()

    call y%init(0.0_rp)
    select type(x)
    class is (serial_block_array_t)
       select type(y)
       class is(serial_block_array_t)
          do ib=1,op%nblocks
             call aux%clone(y%blocks(ib))
             do jb=1,op%nblocks
                if ( associated(op%blocks(ib,jb)%serial_scalar_matrix) ) then
                   ! aux <- A(ib,jb) * x(jb)
                   call op%blocks(ib,jb)%serial_scalar_matrix%apply(x%blocks(jb),aux)
                   ! y(ib) <- y(ib) + aux
                   call y%blocks(ib)%axpby(1.0_rp,aux,1.0_rp)
                end if
             end do
             call aux%free()
          end do
       class default
          write(0,'(a)') 'serial_block_matrix_t%apply: unsupported y class'
          check(1==0)
       end select
    class default
       write(0,'(a)') 'serial_block_matrix_t%apply: unsupported x class'
       check(1==0)
    end select

    call x%CleanTemp()

  end subroutine serial_block_matrix_apply

  ! op%apply(x)
  ! Allocates room for (temporary) y
  function serial_block_matrix_apply_fun(op,x) result(y)
    implicit none
    class(serial_block_matrix_t), intent(in)  :: op
    class(vector_t), intent(in)  :: x
    class(vector_t), allocatable :: y 
    ! Locals
    integer(ip) :: ib,jb
    type(serial_block_array_t), allocatable :: local_y
    type(serial_scalar_array_t) :: aux

    select type(x)
    class is (serial_block_array_t)
       allocate(local_y)
       call local_y%create(op%nblocks)
       do ib=1,op%nblocks
          call aux%clone(local_y%blocks(ib))
          do jb=1,op%nblocks
             ! aux <- A(ib,jb) * x(jb)
             call op%blocks(ib,jb)%serial_scalar_matrix%apply(x%blocks(jb),aux)
             ! y(ib) <- y(ib) + aux
             call local_y%blocks(ib)%axpby(1.0_rp,aux,1.0_rp)
          end do
          call aux%free()
       end do
       call move_alloc(local_y, y)
       call y%SetTemp()
    class default
       write(0,'(a)') 'serial_block_matrix_t%apply_fun: unsupported x class'
       check(1==0)
    end select
  end function serial_block_matrix_apply_fun
  
    subroutine serial_block_matrix_free_in_stages(this,action)
    implicit none
    class(serial_block_matrix_t), intent(inout) :: this
    integer(ip)                 , intent(in)    :: action
    integer(ip) :: ib,jb, istat

    do ib=1, this%nblocks 
       do jb=1, this%nblocks
          if ( associated(this%blocks(ib,jb)%serial_scalar_matrix) ) then
		    if ( action == free_clean ) then
			  call this%blocks(ib,jb)%serial_scalar_matrix%free_in_stages(free_clean)
			  deallocate (this%blocks(ib,jb)%serial_scalar_matrix, stat=istat)
			  check(istat==0)
            else if ( action == free_struct ) then 
			  call this%blocks(ib,jb)%serial_scalar_matrix%free_in_stages(free_struct)
            else if ( action == free_values ) then
			  call this%blocks(ib,jb)%serial_scalar_matrix%free_in_stages(free_values)
            end if
          end if
       end do
    end do
	if (action == free_clean) then
	   this%nblocks = 0
	   deallocate ( this%blocks, stat=istat)
	   check(istat==0)
	end if
  end subroutine serial_block_matrix_free_in_stages
  
    ! Debugged
  subroutine matvec (nv,nv2,ia,ja,a,x,y)
    implicit none
    integer(ip), intent(in)  :: nv,nv2,ia(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: a(ia(nv+1)-1),x(nv2)
    real(rp)   , intent(out) :: y(nv)
    integer(ip)              :: iv,iz,jv

    y = 0.0_rp
    do iv = 1, nv
       do iz = ia(iv), ia(iv+1)-1
          jv   = ja(iz)
          y(iv) = y(iv) + x(jv)*a(iz)
       end do ! iz
    end do ! iv

  end subroutine matvec

  ! Debugged
  subroutine matvec_symmetric_storage (nv,nv2,ia,ja,a,x,y)
    implicit none
    integer(ip), intent(in)  :: nv,nv2,ia(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: a(ia(nv+1)-1),x(nv2)
    real(rp)   , intent(out) :: y(nv)
    integer(ip)              :: iv,iz,jv

    assert(nv==nv2)

    y = 0.0_rp
    do iv = 1, nv
       y(iv) = y(iv) + x(ja(ia(iv)))*a(ia(iv))
       do iz = ia(iv)+1, ia(iv+1)-1
          jv = ja(iz)
          y(iv) = y(iv) + x(jv)*a(iz)
          y(jv) = y(jv) + x(iv)*a(iz)
       end do ! iz
    end do ! iv

  end subroutine matvec_symmetric_storage

  ! Debugged
  subroutine matvec_trans (nv,nv2,ia,ja,a,x,y)
    implicit none
    integer(ip), intent(in)  :: nv,nv2,ia(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: a(ia(nv+1)-1), x(nv)
    real(rp)   , intent(out) :: y(nv2)
    integer(ip)              :: iv,iz,jv
    
    y = 0.0_rp
    do iv = 1, nv
       do iz = ia(iv), ia(iv+1)-1
          jv = ja(iz)
          y(jv) = y(jv) + x(iv)*a(iz)
       end do ! iz
    end do ! iv
    
  end subroutine matvec_trans

  subroutine matvec_symmetric_storage_trans (nv,nv2,ia,ja,a,x,y)
    implicit none
    integer(ip), intent(in)  :: nv,nv2,ia(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: a(ia(nv+1)-1), x(nv)
    real(rp)   , intent(out) :: y(nv2)
    integer(ip)              :: iv,iz,jv,id,jd
    integer(ip)              :: of, ivc, izc, jvc

    write (0,*) 'Error: the body of matvec_symmetric_storage_trans in matvec.f90 still to be written'
    write (0,*) 'Error: volunteers are welcome !!!'
    stop 
  end subroutine matvec_symmetric_storage_trans
  

end module serial_block_matrix_names
