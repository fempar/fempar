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
module par_block_matrix_names
  ! Serial modules
  use types_names
  use memor_names
  use operator_names
  use matrix_names
  use vector_names

  ! Parallel modules
  use par_environment_names
  use par_scalar_matrix_names
  use par_scalar_array_names
  use par_block_array_names
  use dof_distribution_names
  use blocks_dof_distribution_names

  implicit none
# include "debug.i90"

  private

  type p_par_scalar_matrix_t
    type(par_scalar_matrix_t), pointer :: par_scalar_matrix
  end type p_par_scalar_matrix_t

  type, extends(matrix_t) :: par_block_matrix_t
    ! private ! IBM XLF 14.1 bug
    integer(ip) :: nblocks
    type(p_par_scalar_matrix_t), allocatable :: blocks(:,:)
  contains
  
    procedure, private :: par_block_matrix_create_only_blocks_container
    procedure, private :: par_block_matrix_create_blocks_container_and_all_blocks
    generic :: create  => par_block_matrix_create_only_blocks_container, &
	                      par_block_matrix_create_blocks_container_and_all_blocks
						  
	procedure, private :: par_block_matrix_create_diagonal_block
	procedure, private :: par_block_matrix_create_offdiagonal_block
    generic :: create_block => par_block_matrix_create_diagonal_block, &
	                           par_block_matrix_create_offdiagonal_block
							 
	procedure  :: allocate => par_block_matrix_allocate
  
    procedure :: set_block_to_zero       => par_block_matrix_set_block_to_zero
    procedure :: free                    => par_block_matrix_free_in_one_shot
	procedure :: free_in_stages          => par_block_matrix_free_in_stages
    procedure :: get_block               => par_block_matrix_get_block
    procedure :: get_nblocks             => par_block_matrix_get_nblocks
    procedure :: apply                   => par_block_matrix_apply
    procedure :: apply_fun               => par_block_matrix_apply_fun
  end type par_block_matrix_t

  ! Types
  public :: par_block_matrix_t

contains

  !=============================================================================
  subroutine par_block_matrix_create_blocks_container_and_all_blocks(this, & 
																	 blocks_dof_distribution, & 
																	 diagonal_blocks_symmetric_storage, & 
																	 diagonal_blocks_symmetric,& 
																	 diagonal_blocks_sign)
    implicit none
    class(par_block_matrix_t)      , intent(out) :: this
	type(blocks_dof_distribution_t), intent(in)  :: blocks_dof_distribution
	logical                        , intent(in)  :: diagonal_blocks_symmetric_storage(blocks_dof_distribution%nblocks)
    logical                        , intent(in)  :: diagonal_blocks_symmetric(blocks_dof_distribution%nblocks)
    integer(ip)                    , intent(in)  :: diagonal_blocks_sign(blocks_dof_distribution%nblocks)

    integer(ip) :: ib,jb,istat

    call this%create(blocks_dof_distribution%nblocks)
    do ib=1, this%nblocks 
       do jb=1, this%nblocks
         allocate ( this%blocks(ib,jb)%par_scalar_matrix )
         if ( ib == jb ) then
           call this%blocks(ib,jb)%par_scalar_matrix%create( diagonal_blocks_symmetric_storage(ib), & 
															 diagonal_blocks_symmetric(ib), & 
															 diagonal_blocks_sign(ib), &
															 blocks_dof_distribution%get_block(ib),&
															 blocks_dof_distribution%p_env )
         else
           call this%blocks(ib,jb)%par_scalar_matrix%create(blocks_dof_distribution%get_block(ib), &
															blocks_dof_distribution%get_block(jb), &
														    blocks_dof_distribution%p_env )
         end if
       end do
    end do
  end subroutine par_block_matrix_create_blocks_container_and_all_blocks
  
  !=============================================================================
  subroutine par_block_matrix_create_only_blocks_container(this,nblocks)
    implicit none
    ! Parameters
    class(par_block_matrix_t), intent(out) :: this
	integer(ip)              , intent(in)  :: nblocks
	integer(ip) :: istat
    this%nblocks = nblocks
    allocate ( this%blocks(this%nblocks,this%nblocks), stat=istat )
	check (istat==0)
  end subroutine par_block_matrix_create_only_blocks_container
  
  !=============================================================================
  subroutine par_block_matrix_create_diagonal_block(this, & 
													ib, &
							                        dof_distribution, &
													p_env, &
													symmetric_storage, & 
													is_symmetric,& 
													sign)
    implicit none
    class(par_block_matrix_t)      , intent(out) :: this
	integer(ip)                    , intent(in)  :: ib
	type(dof_distribution_t)       , intent(in)  :: dof_distribution
	type(par_environment_t)        , intent(in)  :: p_env
	logical                        , intent(in)  :: symmetric_storage
    logical                        , intent(in)  :: is_symmetric
    integer(ip)                    , intent(in)  :: sign

    integer(ip) :: istat

	assert ( associated ( this%blocks(ib,ib)%par_scalar_matrix ) )
    allocate ( this%blocks(ib,ib)%par_scalar_matrix, stat=istat )
	check ( istat==0 )
	
    call this%blocks(ib,ib)%par_scalar_matrix%create( symmetric_storage, & 
												      is_symmetric, & 
													  sign, &
													  dof_distribution,&
													  p_env )
	end subroutine par_block_matrix_create_diagonal_block
	
  !=============================================================================
  subroutine par_block_matrix_create_offdiagonal_block(this, & 
													   ib, &
							                           dof_distribution, &
													   dof_distribution_cols, &
													   p_env )
    implicit none
    class(par_block_matrix_t)      , intent(out) :: this
	integer(ip)                    , intent(in)  :: ib
	type(dof_distribution_t)       , intent(in)  :: dof_distribution
	type(dof_distribution_t)       , intent(in)  :: dof_distribution_cols
	type(par_environment_t)        , intent(in)  :: p_env

    integer(ip) :: istat

	assert ( associated ( this%blocks(ib,ib)%par_scalar_matrix ) )
    allocate ( this%blocks(ib,ib)%par_scalar_matrix, stat=istat )
	check ( istat==0 )
	
    call this%blocks(ib,ib)%par_scalar_matrix%create( dof_distribution,&
													  dof_distribution_cols,&
													  p_env )

  end subroutine par_block_matrix_create_offdiagonal_block
  
  !=============================================================================
  subroutine par_block_matrix_allocate(this)
    implicit none
    class(par_block_matrix_t), intent(inout) :: this
    integer(ip) :: ib,jb
	
    do ib=1, this%nblocks
       do jb=1, this%nblocks
         if ( associated( this%blocks(ib,jb)%par_scalar_matrix ) ) then
           call this%blocks(ib,jb)%par_scalar_matrix%allocate()
		 end if
       end do
    end do
  end subroutine par_block_matrix_allocate
  
  !=============================================================================
  subroutine par_block_matrix_set_block_to_zero (this,ib,jb)
    implicit none
    ! Parameters
    class(par_block_matrix_t), intent(inout) :: this
    integer(ip)           , intent(in)  :: ib,jb

    if ( associated(this%blocks(ib,jb)%par_scalar_matrix) ) then
       call this%blocks(ib,jb)%par_scalar_matrix%free_in_stages(free_clean)
       deallocate (this%blocks(ib,jb)%par_scalar_matrix)
       nullify ( this%blocks(ib,jb)%par_scalar_matrix )
    end if

  end subroutine par_block_matrix_set_block_to_zero

  function par_block_matrix_get_block (bmat,ib,jb)
    implicit none
    ! Parameters
    class(par_block_matrix_t), target, intent(in) :: bmat
    integer(ip)                    , intent(in) :: ib,jb
    type(par_scalar_matrix_t)               , pointer    :: par_block_matrix_get_block

    par_block_matrix_get_block =>  bmat%blocks(ib,jb)%par_scalar_matrix
  end function par_block_matrix_get_block

  function par_block_matrix_get_nblocks (bmat)
    implicit none
    ! Parameters
    class(par_block_matrix_t), target, intent(in) :: bmat
    integer(ip)                                :: par_block_matrix_get_nblocks
    par_block_matrix_get_nblocks = bmat%nblocks
  end function par_block_matrix_get_nblocks


  ! op%apply(x,y) <=> y <- op*x
  ! Implicitly assumes that y is already allocated
  subroutine par_block_matrix_apply(op,x,y)
    implicit none
    class(par_block_matrix_t), intent(in)    :: op
    class(vector_t)    , intent(in)    :: x
    class(vector_t)    , intent(inouT) :: y
    ! Locals
    integer(ip)        :: ib,jb
    type(par_scalar_array_t) :: aux

    call x%GuardTemp()

    call y%init(0.0_rp)
    select type(x)
    class is (par_block_array_t)
       select type(y)
       class is(par_block_array_t)
          do ib=1,op%nblocks
             call aux%clone(y%blocks(ib))
             do jb=1,op%nblocks
                if ( associated(op%blocks(ib,jb)%par_scalar_matrix) ) then
                   ! aux <- A(ib,jb) * x(jb)
                   call op%blocks(ib,jb)%par_scalar_matrix%apply(x%blocks(jb),aux)
                   ! y(ib) <- y(ib) + aux
                   call y%blocks(ib)%axpby(1.0_rp,aux,1.0_rp)
                end if
             end do
             call aux%free()
          end do
       class default
          write(0,'(a)') 'par_block_matrix_t%apply: unsupported y class'
          check(1==0)
       end select
    class default
       write(0,'(a)') 'par_block_matrix_t%apply: unsupported x class'
       check(1==0)
    end select

    call x%CleanTemp()

  end subroutine par_block_matrix_apply

  ! op%apply(x)
  ! Allocates room for (temporary) y
  function par_block_matrix_apply_fun(op,x) result(y)
    implicit none
    class(par_block_matrix_t), intent(in)  :: op
    class(vector_t)    , intent(in)  :: x
    class(vector_t)    , allocatable :: y 
    ! Locals
    integer(ip) :: ib,jb
    type(par_block_array_t), allocatable :: local_y
    type(par_scalar_array_t) :: aux

    select type(x)
    class is (par_block_array_t)
       allocate(local_y)
       call local_y%create(op%nblocks)
       do ib=1,op%nblocks
          call aux%clone(local_y%blocks(ib))
          do jb=1,op%nblocks
             ! aux <- A(ib,jb) * x(jb)
             call op%blocks(ib,jb)%par_scalar_matrix%apply(x%blocks(jb),aux)
             ! y(ib) <- y(ib) + aux
             call local_y%blocks(ib)%axpby(1.0_rp,aux,1.0_rp)
          end do
          call aux%free()
       end do
       call move_alloc(local_y, y)
       call y%SetTemp()
    class default
       write(0,'(a)') 'par_block_matrix_t%apply_fun: unsupported x class'
       check(1==0)
    end select
  end function par_block_matrix_apply_fun
  
    subroutine par_block_matrix_free_in_one_shot(this)
    implicit none
    class(par_block_matrix_t), intent(inout) :: this
    integer(ip) :: ib,jb
    call this%free_in_stages(free_values)
	call this%free_in_stages(free_struct)
	call this%free_in_stages(free_clean)
  end subroutine par_block_matrix_free_in_one_shot
  
    subroutine par_block_matrix_free_in_stages(this,action)
    implicit none
    class(par_block_matrix_t), intent(inout) :: this
    integer(ip)                 , intent(in)    :: action
    integer(ip) :: ib,jb, istat

    do ib=1, this%nblocks 
       do jb=1, this%nblocks
          if ( associated(this%blocks(ib,jb)%par_scalar_matrix) ) then
		    if ( action == free_clean ) then
			  call this%blocks(ib,jb)%par_scalar_matrix%free_in_stages(free_clean)
			  deallocate (this%blocks(ib,jb)%par_scalar_matrix, stat=istat)
			  check(istat==0)
            else if ( action == free_struct ) then 
			  call this%blocks(ib,jb)%par_scalar_matrix%free_in_stages(free_struct)
            else if ( action == free_values ) then
			  call this%blocks(ib,jb)%par_scalar_matrix%free_in_stages(free_values)
            end if
          end if
       end do
    end do
	if (action == free_clean) then
	   this%nblocks = 0
	   deallocate ( this%blocks, stat=istat)
	   check(istat==0)
	end if
  end subroutine par_block_matrix_free_in_stages

end module par_block_matrix_names
