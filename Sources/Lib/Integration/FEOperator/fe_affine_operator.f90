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
module fe_affine_operator_names
  use types_names
  use memor_names
  
  use vector_space_names
  use reference_fe_names
  use fe_space_names
  use operator_names
  use vector_names
  
  use matrix_array_assembler_names
  use sparse_matrix_array_assembler_names
  use block_sparse_matrix_array_assembler_names  
  use par_sparse_matrix_array_assembler_names

  use sparse_matrix_names, only: sparse_matrix_t
  use block_sparse_matrix_names
  use par_sparse_matrix_names
  
  use serial_scalar_array_names
  use serial_block_array_names
  use par_scalar_array_names
  
  use array_names
  use matrix_names
  use discrete_integration_names
  use environment_names
  use direct_solver_names
  use block_layout_names

  implicit none
# include "debug.i90"

  private

  integer(ip), parameter :: start               = 0 
  integer(ip), parameter :: created             = 1 
  integer(ip), parameter :: symbolically_setup  = 2  
  integer(ip), parameter :: numerically_setup   = 3 

  ! State transition diagram for type(fe_affine_operator_t)
  ! -------------------------------------------------
  ! Input State | Action               | Output State 
  ! -------------------------------------------------
  ! start       | create               | created
  ! start       | free_numerical_setup          | start
  ! start       | free_symbolic_setup  | start
  ! start       | free_clean           | start 
  ! start       | free                 | start 

  ! created     | symbolic_setup       | symbolically_setup
  ! created     | numerical_setup      | numerically_setup
  ! created     | get_tangent+         | numerically_setup
  !               get_translation+               "
  !               get_domain_vector_space+       "
  !               get_range_vector_space+        "
  !               abort_if_not_in_range+         "
  !               abort_if_not_in_domain         "
  ! created     | free_numerical_setup          | created
  ! created     | free_symbolic_setup          | created
  ! created     | free_clean           | start
  ! created     | free                 | start

  ! symbolically_setup    | symbolic_setup       | symbolically_setup
  ! symbolically_setup    | numerical_setup      | numerically_setup
  ! symbolically_setup    | free_numerical_setup          | symbolically_setup
  ! symbolically_setup    | free_symbolic_setup          | created
  ! symbolically_setup    | free_clean           | start
  ! symbolically_setup    | free                 | start
  ! symbolically_setup    | get_tangent+         | numerically_setup
  !                         get_translation+               "
  !                         get_domain_vector_space+       "
  !                         get_range_vector_space+        "
  !                         abort_if_not_in_range+         "
  !                         abort_if_not_in_domain         "

  ! numerically_setup    | symbolic_setup       | numerically_setup
  ! numerically_setup    | numerical_setup      | numerically_setup
  ! numerically_setup    | free_numerical_setup          | symbolically_setup
  ! numerically_setup    | free_symbolic_setup  | created
  ! numerically_setup    | free                 | start
  ! numerically_setup    | free_clean           |  
  ! numerically_setup    | get_tangent+         | numerically_setup
  !                        get_translation+               "
  !                        get_domain_vector_space+       "
  !                        get_range_vector_space+        "
  !                        abort_if_not_in_range+         "
  !                        abort_if_not_in_domain         "


 type, extends(operator_t):: fe_affine_operator_t
  private
  integer(ip)                                     :: state  = start
  character(:)                      , allocatable :: sparse_matrix_storage_format
  type(block_layout_t)                            :: block_layout
  class(serial_fe_space_t)          , pointer     :: fe_space               => NULL() ! test_fe_space
  class(serial_fe_space_t)          , pointer     :: trial_fe_space         => NULL() ! To be used in the future
  class(discrete_integration_t)     , pointer     :: discrete_integration   => NULL()
  class(matrix_array_assembler_t)   , pointer     :: matrix_array_assembler => NULL()
contains
  procedure          :: create                      => fe_affine_operator_create
  procedure          :: symbolic_setup              => fe_affine_operator_symbolic_setup
  procedure          :: numerical_setup             => fe_affine_operator_numerical_setup
  procedure          :: apply                       => fe_affine_operator_apply
  procedure          :: apply_add                   => fe_affine_operator_apply_add
  procedure          :: is_linear                   => fe_affine_operator_is_linear
  procedure          :: get_tangent                 => fe_affine_operator_get_tangent
  procedure          :: get_translation             => fe_affine_operator_get_translation
  procedure          :: get_matrix                  => fe_affine_operator_get_matrix
  procedure          :: get_array                   => fe_affine_operator_get_array
  procedure          :: get_fe_space                => fe_affine_operator_get_fe_space
  procedure          :: free_in_stages              => fe_affine_operator_free_in_stages
  procedure          :: free                        => fe_affine_operator_free
  procedure          :: get_domain_vector_space     => fe_affine_operator_get_domain_vector_space
  procedure          :: get_range_vector_space      => fe_affine_operator_get_range_vector_space
  procedure          :: abort_if_not_in_range       => fe_affine_operator_abort_if_not_in_range
  procedure          :: abort_if_not_in_domain      => fe_affine_operator_abort_if_not_in_domain
  procedure          :: create_direct_solver        => fe_affine_operator_create_direct_solver
  procedure          :: update_direct_solver_matrix => fe_affine_operator_update_direct_solver_matrix
  procedure, private :: create_serial_assembler     => fe_affine_operator_create_serial_assembler
  procedure, private :: create_par_assembler        => fe_affine_operator_create_par_assembler
  procedure, private :: fe_affine_operator_free_numerical_setup
  procedure, private :: fe_affine_operator_free_symbolic_setup
  procedure, private :: fe_affine_operator_free_clean
  procedure, private :: fe_affine_operator_setup
  procedure, private :: fe_affine_operator_fill_values
  procedure, private :: create_vector_spaces        => fe_affine_operator_create_vector_spaces
end type fe_affine_operator_t

! Types
public :: fe_affine_operator_t

contains
subroutine fe_affine_operator_create (this, &
                                      sparse_matrix_storage_format, &
                                      diagonal_blocks_symmetric_storage,&
                                      diagonal_blocks_symmetric,&
                                      diagonal_blocks_sign,&
                                      fe_space,&
                                      discrete_integration, &
                                      field_blocks, &
                                      field_coupling, &
                                      trial_fe_space)
 implicit none
 class(fe_affine_operator_t)              , intent(inout) :: this
 character(*)                             , intent(in)    :: sparse_matrix_storage_format
 logical                                  , intent(in)    :: diagonal_blocks_symmetric_storage(:)
 logical                                  , intent(in)    :: diagonal_blocks_symmetric(:)
 integer(ip)                              , intent(in)    :: diagonal_blocks_sign(:)
 class(serial_fe_space_t)        , target, intent(inout) :: fe_space
 class(discrete_integration_t)   , target, intent(in)    :: discrete_integration
 integer(ip)                     , optional, intent(in)    :: field_blocks(:)
 logical                         , optional, intent(in)    :: field_coupling(:,:)
 class(serial_fe_space_t), target, optional, intent(inout) :: trial_fe_space
 
 call this%free()

 assert(this%state == start)

#ifdef DEBUG
 if ( present(field_blocks) ) then
   assert ( size(field_blocks) == fe_space%get_number_fields() )
 end if
 if ( present(field_coupling) ) then
   assert ( size(field_coupling,1) == fe_space%get_number_fields() )
   assert ( size(field_coupling,2) == fe_space%get_number_fields() )
 end if
#endif
 
 this%sparse_matrix_storage_format = sparse_matrix_storage_format
 this%fe_space                     => fe_space
 this%discrete_integration         => discrete_integration
 
 call this%block_layout%create(fe_space%get_number_fields(), field_blocks, field_coupling )
 call this%fe_space%fill_dof_info(this%block_layout)
 if ( present(trial_fe_space) ) then
    this%trial_fe_space => trial_fe_space
    call this%trial_fe_space%fill_dof_info(this%block_layout)
 end if
 
  select type(fe_space => this%fe_space)
  class is(serial_fe_space_t) 
    this%matrix_array_assembler  => this%create_serial_assembler(diagonal_blocks_symmetric_storage, &
                                                                 diagonal_blocks_symmetric, &
                                                                 diagonal_blocks_sign)
  class is(par_fe_space_t) 
    this%matrix_array_assembler  => this%create_par_assembler(diagonal_blocks_symmetric_storage, &
                                                              diagonal_blocks_symmetric, &
                                                              diagonal_blocks_sign)
  class default
    check(.false.)
  end select
    
 call this%create_vector_spaces()
 this%state = created
end subroutine fe_affine_operator_create

  subroutine fe_affine_operator_create_vector_spaces(this)
    implicit none
    class(fe_affine_operator_t), intent(inout) :: this
    type(vector_space_t), pointer                 :: fe_affine_operator_domain_vector_space
    type(vector_space_t), pointer                 :: fe_affine_operator_range_vector_space
    type(vector_space_t), pointer                 :: matrix_domain_vector_space
    type(vector_space_t), pointer                 :: matrix_range_vector_space
    class(matrix_t)     , pointer :: matrix
    matrix => this%matrix_array_assembler%get_matrix()
    matrix_domain_vector_space => matrix%get_domain_vector_space()
    matrix_range_vector_space => matrix%get_range_vector_space()
    fe_affine_operator_domain_vector_space => operator_get_domain_vector_space(this)
    fe_affine_operator_range_vector_space => operator_get_range_vector_space(this)
    call matrix_domain_vector_space%clone(fe_affine_operator_domain_vector_space)
    call matrix_range_vector_space%clone(fe_affine_operator_range_vector_space)
  end subroutine fe_affine_operator_create_vector_spaces


subroutine fe_affine_operator_symbolic_setup (this)
 implicit none
 class(fe_affine_operator_t), intent(inout) :: this
 assert ( .not. this%state == start )
 if ( this%state == created ) then 
    this%state = symbolically_setup
 end if
end subroutine fe_affine_operator_symbolic_setup

subroutine fe_affine_operator_numerical_setup (this)
 implicit none
 class(fe_affine_operator_t), intent(inout) :: this

 assert ( .not. this%state == start )

 if ( this%state == created ) then
    call this%symbolic_setup()
 end if
 if ( this%state == symbolically_setup ) then
    this%state = numerically_setup
    call this%matrix_array_assembler%allocate_array()
    call this%matrix_array_assembler%init_array(0.0_rp)
 elseif ( this%state == numerically_setup ) then
    call this%matrix_array_assembler%init_array(0.0_rp)
    call this%matrix_array_assembler%init_matrix(0.0_rp)
 end if

 call this%fe_affine_operator_fill_values()
 call this%matrix_array_assembler%compress_storage(this%sparse_matrix_storage_format)
end subroutine fe_affine_operator_numerical_setup

function fe_affine_operator_create_serial_assembler (this, &
                                                     diagonal_blocks_symmetric_storage,&
                                                     diagonal_blocks_symmetric, & 
                                                     diagonal_blocks_sign)
  implicit none
  class(fe_affine_operator_t)     , intent(in) :: this
  logical                         , intent(in) :: diagonal_blocks_symmetric_storage(:)
  logical                         , intent(in) :: diagonal_blocks_symmetric(:)
  integer(ip)                     , intent(in) :: diagonal_blocks_sign(:)
  class(matrix_array_assembler_t) , pointer    :: fe_affine_operator_create_serial_assembler

  ! Locals
  class(matrix_t), pointer :: matrix
  class(array_t) , pointer :: array
  integer(ip)          :: ife_space, jfe_space
  integer(ip)          :: iblock, jblock

  if (this%block_layout%get_num_blocks() == 1) then
     allocate ( sparse_matrix_array_assembler_t :: fe_affine_operator_create_serial_assembler )
     allocate ( sparse_matrix_t :: matrix )
     allocate ( serial_scalar_array_t  :: array )
     select type(matrix)
        class is(sparse_matrix_t)
        call matrix%create(this%block_layout%get_block_num_dofs(1), &
                           diagonal_blocks_symmetric_storage(1),&
                           diagonal_blocks_symmetric(1),&
                           diagonal_blocks_sign(1))
        class default
        check(.false.)
     end select
     select type(array)
        class is(serial_scalar_array_t)
        call array%create(this%block_layout%get_block_num_dofs(1))
        class default
        check(.false.)
     end select
  else
     allocate ( block_sparse_matrix_array_assembler_t :: fe_affine_operator_create_serial_assembler )
     allocate ( block_sparse_matrix_t :: matrix )
     allocate ( serial_block_array_t  :: array )
     select type(matrix)
        class is (block_sparse_matrix_t)
        call matrix%create(this%block_layout%get_num_blocks(), &
             this%block_layout%get_num_dofs_x_block(),&
             this%block_layout%get_num_dofs_x_block(),&
             diagonal_blocks_symmetric_storage,&
             diagonal_blocks_symmetric,&
             diagonal_blocks_sign)

        do jblock=1,this%block_layout%get_num_blocks()
           do iblock=1,this%block_layout%get_num_blocks()
              if (.not. this%block_layout%blocks_coupled(iblock,jblock) ) then
                 call matrix%set_block_to_zero(iblock,jblock)
              end if
           end do
        end do
        class default
        check(.false.)
     end select
     select type(array)
        class is(serial_block_array_t)
        call array%create(this%block_layout%get_num_blocks(),this%block_layout%get_num_dofs_x_block())
        class default
        check(.false.)
     end select
  end if
  call fe_affine_operator_create_serial_assembler%set_matrix(matrix)
  call fe_affine_operator_create_serial_assembler%set_array(array)
end function fe_affine_operator_create_serial_assembler

function fe_affine_operator_create_par_assembler(this, &
                                                 diagonal_blocks_symmetric_storage,&
                                                 diagonal_blocks_symmetric, & 
                                                 diagonal_blocks_sign)
  implicit none
  class(fe_affine_operator_t)       , intent(in) :: this
  logical                           , intent(in) :: diagonal_blocks_symmetric_storage(:)
  logical                           , intent(in) :: diagonal_blocks_symmetric(:)
  integer(ip)                       , intent(in) :: diagonal_blocks_sign(:)
  class(matrix_array_assembler_t)   , pointer    :: fe_affine_operator_create_par_assembler

  ! Locals
  class(matrix_t), pointer :: matrix
  class(array_t) , pointer :: array
  type(environment_t), pointer :: par_environment
  
  
  select type(fe_space => this%fe_space)
  class is(par_fe_space_t)
   par_environment => fe_space%get_par_environment()
   if (this%block_layout%get_num_blocks() == 1) then
     allocate ( par_sparse_matrix_array_assembler_t :: fe_affine_operator_create_par_assembler )
     allocate ( par_sparse_matrix_t :: matrix )
     allocate ( par_scalar_array_t  :: array )
     select type(matrix)
        class is(par_sparse_matrix_t)
        call matrix%create(par_environment, &
                           fe_space%get_block_dof_import(1), &
                           diagonal_blocks_symmetric_storage(1),&
                           diagonal_blocks_symmetric(1),&
                           diagonal_blocks_sign(1))
        class default
        check(.false.)
     end select
     select type(array)
        class is(par_scalar_array_t)
        call array%create(par_environment, &
                          fe_space%get_block_dof_import(1))
        class default
        check(.false.)
     end select
   else
     check(.false.)
   end if
   call fe_affine_operator_create_par_assembler%set_matrix(matrix)
   call fe_affine_operator_create_par_assembler%set_array(array)
  class default
   check(.false.)
  end select
end function fe_affine_operator_create_par_assembler

subroutine fe_affine_operator_free_numerical_setup(this)
 implicit none
 class(fe_affine_operator_t), intent(inout) :: this
 call this%matrix_array_assembler%free_in_stages(free_numerical_setup)
end subroutine fe_affine_operator_free_numerical_setup

subroutine fe_affine_operator_free_symbolic_setup(this)
 implicit none
 class(fe_affine_operator_t), intent(inout) :: this
 call this%matrix_array_assembler%free_in_stages(free_symbolic_setup)
end subroutine fe_affine_operator_free_symbolic_setup

subroutine fe_affine_operator_free_clean(this)
 implicit none
 class(fe_affine_operator_t), intent(inout) :: this
 integer(ip) :: istat
 deallocate(this%sparse_matrix_storage_format)
 nullify(this%fe_space)
 nullify(this%trial_fe_space)
 nullify(this%discrete_integration)
 call this%matrix_array_assembler%free_in_stages(free_clean)
 deallocate(this%matrix_array_assembler, stat=istat )
 check(istat==0)
 nullify(this%matrix_array_assembler)
 call this%free_vector_spaces()
 call this%block_layout%free()
end subroutine fe_affine_operator_free_clean

subroutine fe_affine_operator_free_in_stages(this,action)
 implicit none
 class(fe_affine_operator_t), intent(inout) :: this
 integer(ip)                , intent(in)    :: action
 integer(ip)                                :: istat

 if ( this%state == start ) then
    return
 else if ( this%state == created ) then
    if ( action == free_clean ) then
       call this%fe_affine_operator_free_clean()
       this%state = start
    end if
 else if ( this%state == symbolically_setup ) then
    if ( action == free_symbolic_setup ) then
       call this%fe_affine_operator_free_symbolic_setup()
       this%state = created
    else if ( action == free_clean ) then
       call this%fe_affine_operator_free_symbolic_setup()
       call this%fe_affine_operator_free_clean()
       this%state=start
    end if
 else if ( this%state == numerically_setup ) then
    if ( action == free_numerical_setup ) then
       call this%fe_affine_operator_free_numerical_setup()
       this%state = symbolically_setup
    else if ( action == free_symbolic_setup ) then
       call this%fe_affine_operator_free_numerical_setup()
       call this%fe_affine_operator_free_symbolic_setup()
       this%state = created
    else if ( action == free_clean ) then
       call this%fe_affine_operator_free_numerical_setup()
       call this%fe_affine_operator_free_symbolic_setup()
       call this%fe_affine_operator_free_clean()
       this%state=start
    end if
 end if

end subroutine fe_affine_operator_free_in_stages

subroutine fe_affine_operator_free(this)
 implicit none
 class(fe_affine_operator_t), intent(inout) :: this
 call this%free_in_stages(free_numerical_setup)
 call this%free_in_stages(free_symbolic_setup)
 call this%free_in_stages(free_clean)
end subroutine fe_affine_operator_free

function fe_affine_operator_get_matrix(this)
 implicit none
 class(fe_affine_operator_t), target, intent(in) :: this
 class(matrix_t), pointer :: fe_affine_operator_get_matrix
 call this%fe_affine_operator_setup()
 fe_affine_operator_get_matrix => this%matrix_array_assembler%get_matrix()
end function fe_affine_operator_get_matrix

function fe_affine_operator_get_array(this)
 implicit none
 class(fe_affine_operator_t), target, intent(in) :: this
 class(array_t), pointer :: fe_affine_operator_get_array
 call this%fe_affine_operator_setup()
 fe_affine_operator_get_array => this%matrix_array_assembler%get_array()
end function fe_affine_operator_get_array

function fe_affine_operator_get_fe_space(this)
 implicit none
 class(fe_affine_operator_t), target, intent(in) :: this
 class(serial_fe_space_t), pointer :: fe_affine_operator_get_fe_space
 assert ( .not. this%state == start )
 fe_affine_operator_get_fe_space => this%fe_space
end function fe_affine_operator_get_fe_space

! op%apply(x,y) <=> y <- op*x
! Implicitly assumes that y is already allocated
subroutine fe_affine_operator_apply(this,x,y) 
 implicit none
 class(fe_affine_operator_t), intent(in)    :: this
 class(vector_t) , intent(in)    :: x
 class(vector_t) , intent(inout) :: y 
 class(matrix_t) , pointer       :: matrix
 class(array_t)  , pointer       :: array
 call this%fe_affine_operator_setup()
 call this%abort_if_not_in_domain(x)
 call this%abort_if_not_in_range(y)
 call x%GuardTemp()
 matrix => this%matrix_array_assembler%get_matrix()
 call matrix%apply(x,y)
 array => this%matrix_array_assembler%get_array()
 call y%axpby( -1.0_rp, array, 1.0_rp )
 call x%CleanTemp()
end subroutine fe_affine_operator_apply

! op%apply(x,y) <=> y <- op*x+y
! Implicitly assumes that y is already allocated
subroutine fe_affine_operator_apply_add(this,x,y) 
 implicit none
 class(fe_affine_operator_t), intent(in)    :: this
 class(vector_t) , intent(in)    :: x
 class(vector_t) , intent(inout) :: y 
 class(matrix_t) , pointer       :: matrix
 class(array_t)  , pointer       :: array
 call this%fe_affine_operator_setup()
 call this%abort_if_not_in_domain(x)
 call this%abort_if_not_in_range(y)
 call x%GuardTemp()
 matrix => this%matrix_array_assembler%get_matrix()
 call matrix%apply_add(x,y)
 array => this%matrix_array_assembler%get_array()
 call y%axpby( -1.0_rp, array, 1.0_rp )
 call x%CleanTemp()
end subroutine fe_affine_operator_apply_add

function fe_affine_operator_is_linear(this)
 implicit none
 class(fe_affine_operator_t), intent(in) :: this
 logical :: fe_affine_operator_is_linear
 fe_affine_operator_is_linear = .false.
end function fe_affine_operator_is_linear

function fe_affine_operator_get_tangent(this,x) result(tangent)
 implicit none
 class(fe_affine_operator_t), intent(in) :: this
 class(vector_t), optional  , intent(in) :: x
 type(lvalue_operator_t)          :: tangent
 call this%fe_affine_operator_setup()
 tangent = this%get_matrix()
 call tangent%SetTemp()
end function fe_affine_operator_get_tangent

function fe_affine_operator_get_translation(this) result(translation)
 implicit none
 class(fe_affine_operator_t), intent(in) :: this
 class(vector_t), pointer                :: translation
 call this%fe_affine_operator_setup()
 translation => this%get_array()
end function fe_affine_operator_get_translation

subroutine fe_affine_operator_abort_if_not_in_domain ( this, vector )
 implicit none
 class(fe_affine_operator_t), intent(in)  :: this
 class(vector_t)            , intent(in)  :: vector
 assert ( .not. this%state == start )
 call operator_abort_if_not_in_domain(this,vector)
end subroutine fe_affine_operator_abort_if_not_in_domain

subroutine fe_affine_operator_abort_if_not_in_range ( this, vector )
 implicit none
 class(fe_affine_operator_t), intent(in) :: this
 class(vector_t)            , intent(in) :: vector
 assert ( .not. this%state == start )
 call operator_abort_if_not_in_range(this,vector)
end subroutine fe_affine_operator_abort_if_not_in_range

function fe_affine_operator_get_domain_vector_space ( this )
 implicit none
 class(fe_affine_operator_t), target, intent(in) :: this
 type(vector_space_t)               , pointer    :: fe_affine_operator_get_domain_vector_space
 assert ( .not. this%state == start )
 fe_affine_operator_get_domain_vector_space => operator_get_domain_vector_space(this)
end function fe_affine_operator_get_domain_vector_space

function fe_affine_operator_get_range_vector_space ( this )
 implicit none
 class(fe_affine_operator_t), target, intent(in) :: this
 type(vector_space_t)                  , pointer :: fe_affine_operator_get_range_vector_space
 assert ( .not. this%state == start )
 fe_affine_operator_get_range_vector_space => operator_get_range_vector_space(this)
end function fe_affine_operator_get_range_vector_space

function fe_affine_operator_apply_fun(op,x) result(y)
 implicit none
 class(fe_affine_operator_t)   , intent(in)  :: op
 class(vector_t)     , intent(in)  :: x
 class(vector_t)     , allocatable :: y
 type(vector_space_t), pointer     :: range_vector_space
 range_vector_space => op%get_range_vector_space()
 call op%fe_affine_operator_setup()
 call range_vector_space%create_vector(y)
 call op%apply(x,y)
end function fe_affine_operator_apply_fun

subroutine fe_affine_operator_setup(this)
 implicit none
 class(fe_affine_operator_t)  :: this
  assert (this%state/=start)
  if(this%state==created .or. this%state==symbolically_setup) call this%numerical_setup()
end subroutine fe_affine_operator_setup

subroutine fe_affine_operator_fill_values(this)
  implicit none
  class(fe_affine_operator_t), intent(inout) :: this
  class(environment_t), pointer :: environment
  environment => this%fe_space%get_environment()
  if ( environment%am_i_l1_task() ) then
    if ( associated(this%trial_fe_space) ) then
     call this%discrete_integration%integrate( this%fe_space, this%trial_fe_space, this%matrix_array_assembler )
    else
     call this%discrete_integration%integrate( this%fe_space, this%matrix_array_assembler )
    end if  
  end if  
end subroutine fe_affine_operator_fill_values


subroutine fe_affine_operator_create_direct_solver(this, name, direct_solver)
  implicit none
  class(fe_affine_operator_t), intent(in)    :: this
  character(len=*),            intent(in)    :: name
  type(direct_solver_t),       intent(inout) :: direct_solver
  call direct_solver%set_type(name)
  select type(matrix => this%matrix_array_assembler%get_matrix())
    type is (sparse_matrix_t)
      call direct_solver%set_matrix(matrix)
    class DEFAULT
      check(.false.)
  end select
end subroutine fe_affine_operator_create_direct_solver


subroutine fe_affine_operator_update_direct_solver_matrix(this, same_nonzero_pattern, direct_solver)
  implicit none
  class(fe_affine_operator_t), intent(in)    :: this
  logical,                     intent(in)    :: same_nonzero_pattern
  type(direct_solver_t),       intent(inout) :: direct_solver
  select type(matrix => this%matrix_array_assembler%get_matrix())
    type is (sparse_matrix_t)
      call direct_solver%update_matrix(matrix, same_nonzero_pattern)
    class DEFAULT
      check(.false.)
  end select
end subroutine fe_affine_operator_update_direct_solver_matrix


end module fe_affine_operator_names
