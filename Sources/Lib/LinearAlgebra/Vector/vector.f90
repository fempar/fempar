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
module vector_names
  use types_names
  use memory_guard_names
  implicit none

  private
# include "debug.i90"
  
  type, abstract, extends(memory_guard_t) :: vector_t
   contains
     procedure (allocate_interface), deferred :: allocate
     procedure (dot_interface) , deferred     :: dot
     procedure (local_dot_interface) , deferred :: local_dot
     procedure (copy_interface), deferred     :: copy
     procedure (init_interface), deferred     :: init
     procedure (scal_interface), deferred     :: scal
     procedure (axpby_interface), deferred    :: axpby ! y <- a*x + b*y
     procedure (nrm2_interface), deferred     :: nrm2 
     procedure (clone_interface), deferred    :: clone
     procedure (comm_interface), deferred     :: comm
     procedure (same_vector_space_interface), deferred :: same_vector_space
     procedure (get_num_blocks_interface), deferred :: get_num_blocks
     procedure (vector_extract_subvector_interface), deferred :: extract_subvector 
     procedure (vector_insert_subvector_interface) , deferred :: insert_subvector 
     procedure (entrywise_product_interface), deferred :: entrywise_product
     procedure (entrywise_invert_interface), deferred :: entrywise_invert
     
     procedure :: mold 
     procedure :: sum_vector
     procedure :: sub_vector
     procedure :: minus_vector
     procedure, pass(left)  :: scal_left_vector
     procedure, pass(right) :: scal_right_vector
     procedure :: assign_vector

     generic  :: operator(+) => sum_vector
     generic  :: operator(-) => sub_vector, minus_vector
     generic  :: operator(*) => scal_left_vector, scal_right_vector
     generic  :: assignment(=) => assign_vector
  end type vector_t

  abstract interface
     subroutine allocate_interface(this)
       import :: vector_t
       implicit none
       class(vector_t)         ,intent(inout) :: this 
     end subroutine allocate_interface
     ! alpha <- op1^T * op2
     function dot_interface(op1,op2) result(alpha)
       import :: vector_t, rp
       implicit none
       class(vector_t), intent(in)  :: op1,op2
       real(rp) :: alpha
     end function dot_interface
     ! alpha <- op1^T * op2 without final allreduce
     function local_dot_interface(op1,op2) result(alpha)
       import :: vector_t, rp
       implicit none
       class(vector_t), intent(in)  :: op1,op2
       real(rp) :: alpha
     end function local_dot_interface
     ! op1 <- op2 
     subroutine copy_interface(op1,op2)
       import :: vector_t
       implicit none
       class(vector_t), intent(inout) :: op1
       class(vector_t), intent(in)    :: op2
     end subroutine copy_interface
     ! op1 <- alpha * op2
     subroutine scal_interface(op1,alpha,op2)
       import :: vector_t, rp
       implicit none
       class(vector_t), intent(inout) :: op1
       real(rp), intent(in) :: alpha
       class(vector_t), intent(in) :: op2
     end subroutine scal_interface
     ! op <- alpha
     subroutine init_interface(op,alpha)
       import :: vector_t, rp
       implicit none
       class(vector_t), intent(inout) :: op
       real(rp), intent(in) :: alpha
     end subroutine init_interface
     ! op1 <- alpha*op2 + beta*op1
     subroutine axpby_interface(op1, alpha, op2, beta)
       import :: vector_t, rp
       implicit none
       class(vector_t), intent(inout) :: op1
       real(rp), intent(in) :: alpha
       class(vector_t), intent(in) :: op2
       real(rp), intent(in) :: beta
     end subroutine axpby_interface
     ! alpha <- nrm2(op)
     function nrm2_interface(op) result(alpha)
       import :: vector_t, rp
       implicit none
       class(vector_t), intent(in)  :: op
       real(rp) :: alpha
     end function nrm2_interface
     ! op1 <- clone(op2) 
     subroutine clone_interface(op1,op2)
       import :: vector_t
       implicit none
       class(vector_t),target         ,intent(inout) :: op1
       class(vector_t),target         ,intent(in)    :: op2
     end subroutine clone_interface
     ! op <- comm(op), i.e. fully assembled op <- subassembled op 
     subroutine comm_interface(op)
       import :: vector_t
       implicit none
       class(vector_t), intent(inout) :: op
     end subroutine comm_interface
     ! Determines whether this belongs to the same
     ! vector space as vector   
     function same_vector_space_interface(this,vector)
       import :: vector_t
       implicit none
       class(vector_t), intent(in) :: this
       class(vector_t), intent(in) :: vector
       logical :: same_vector_space_interface
     end function
     ! Provide the number of blocks of the vector
     function get_num_blocks_interface(this) result (res)
       import :: vector_t, ip
       implicit none 
       class(vector_t), intent(in)   :: this
       integer(ip) :: res
     end function get_num_blocks_interface
     ! Extract subvector from a subset of indices of a vector
     subroutine vector_extract_subvector_interface( this, &
                                                  & iblock, &
                                                  & size_indices, &
                                                  & indices, &
                                                  & values )
       import :: vector_t, ip, rp
       implicit none
       class(vector_t), intent(in)     :: this 
       integer(ip)    , intent(in)     :: iblock
       integer(ip)    , intent(in)     :: size_indices
       integer(ip)    , intent(in)     :: indices(size_indices)
       real(rp)       , intent(inout)  :: values(*)
     end subroutine vector_extract_subvector_interface
     ! Insert subvector from a subset of indices of a vector
     subroutine vector_insert_subvector_interface( this, &
                                                  & iblock, &
                                                  & size_indices, &
                                                  & indices, &
                                                  & values )
       import :: vector_t, ip, rp
       implicit none
       class(vector_t), intent(inout)  :: this 
       integer(ip)    , intent(in)     :: iblock
       integer(ip)    , intent(in)     :: size_indices
       integer(ip)    , intent(in)     :: indices(size_indices)
       real(rp)       , intent(in)     :: values(*)
     end subroutine vector_insert_subvector_interface
     ! op1 <- op2 ⚪ op3: op1 is the entrywise or Hadamard product of op2 and op3,
     ! i.e. (op1_1,...,op1_n) = (op2_1*op3_1,...,op2_n*op3_n), i = 1:n, n=size(op1)
     subroutine entrywise_product_interface(op1,op2,op3)
       import :: vector_t
       implicit none
       class(vector_t), intent(inout) :: op1
       class(vector_t), intent(in)    :: op2
       class(vector_t), intent(in)    :: op3
     end subroutine entrywise_product_interface
     ! op1_i <- 1 / op1_i, i = 1:n, n=size(op1), i.e. invert each entry of op1
     subroutine entrywise_invert_interface(op1)
       import :: vector_t
       implicit none
       class(vector_t), intent(inout) :: op1
     end subroutine entrywise_invert_interface
  end interface

  public :: vector_t

contains  

  subroutine mold(this,vector) 
    implicit none
    class(vector_t)             , intent(in)    :: this
    class(vector_t), allocatable, intent(inout) :: vector
    integer(ip) :: istat
    
    if ( allocated(vector) ) then
      if ( .not. same_type_as(this,vector) ) then
         call vector%free()
         deallocate(vector, stat=istat); check(istat==0);
         allocate(vector, mold=this, stat=istat); check(istat==0);
      end if
    else
      allocate(vector, mold=this, stat=istat); check(istat==0);
    end if
  end subroutine mold

  ! res <- op1 + op2
  function sum_vector(op1,op2) result (res)
    implicit none
    class(vector_t), intent(in)  :: op1, op2
    class(vector_t), allocatable :: res
    
    call op1%GuardTemp()
    call op2%GuardTemp()

    allocate(res, mold=op1); call res%default_initialization()
    call res%clone(op1)
    call res%copy(op1)
    ! res <- 1.0*op2 + 1.0*res
    call res%axpby(1.0, op2, 1.0)

    call op1%CleanTemp()
    call op2%CleanTemp()
    call res%SetTemp()
  end function sum_vector

  ! res <- op1 - op2
  function sub_vector(op1,op2) result (res)
    implicit none
    class(vector_t), intent(in)  :: op1, op2
    class(vector_t), allocatable :: res

    call op1%GuardTemp()
    call op2%GuardTemp()

    allocate(res, mold=op1); call res%default_initialization()
    call res%clone(op1)
    call res%copy(op1)
    ! res <- -1.0*op2 + 1.0*res
    call res%axpby(-1.0, op2, 1.0)

    call op1%CleanTemp()
    call op2%CleanTemp()
    call res%SetTemp()
  end function sub_vector
  
  ! res <- -op
  function minus_vector(op) result (res)
    implicit none
    class(vector_t), intent(in)   :: op
    class(vector_t), allocatable  :: res

    call op%GuardTemp()
    
    allocate(res, mold=op); call res%default_initialization()
    call res%clone(op)
    call res%scal(-1.0, op)

    call op%CleanTemp()
    call res%SetTemp()
  end function minus_vector

    ! res <- op*alpha
  function scal_left_vector(left,alpha) result (res)
    implicit none
    class(vector_t), intent(in)   :: left
    real (rp), intent(in)             :: alpha
    class(vector_t), allocatable  :: res

    call left%GuardTemp()
    
    allocate(res, mold=left); call res%default_initialization()
    call res%clone(left)
    call res%scal(alpha, left)

    call left%CleanTemp()
    call res%SetTemp()
  end function scal_left_vector

  ! res <- alpha*op
  function scal_right_vector(alpha,right) result (res)
    implicit none
    real (rp), intent(in)             :: alpha
    class(vector_t), intent(in)   :: right
    class(vector_t), allocatable  :: res

    call right%GuardTemp()
    
    allocate(res, mold=right); call res%default_initialization()
    call res%clone(right)
    call res%scal(alpha, right)
    
    call right%CleanTemp()
    call res%SetTemp()
  end function scal_right_vector

  ! op1 <- op2
  subroutine assign_vector(op1,op2) 
    implicit none
    class(vector_t), target, intent(inout):: op1
    class(vector_t), target, intent(in):: op2
    class(vector_t), pointer :: p
    p => op1
    if(associated(p,op2)) return ! It's aliasing

    call op2%GuardTemp()
    
    call op1%clone(op2)
    call op1%copy(op2)
    
    call op2%CleanTemp()
  end subroutine assign_vector
    
end module vector_names
