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
  use problem_names
  
  ! Abstract types
  use vector_space_names
  use dof_descriptor_names
  use fe_space_names
  use operator_names
  use vector_names
  use matrix_array_assembler_names
  use array_names
  use matrix_names

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
  ! start       | free_values          | start
  ! start       | free_struct          | start
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
  ! created     | free_values          | created
  ! created     | free_struct          | created
  ! created     | free_clean           | start
  ! created     | free                 | start

  ! symbolically_setup    | symbolic_setup       | symbolically_setup
  ! symbolically_setup    | numerical_setup      | numerically_setup
  ! symbolically_setup    | free_values          | symbolically_setup
  ! symbolically_setup    | free_struct          | created
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
  ! numerically_setup    | free_values          | symbolically_setup
  ! numerically_setup    | free_struct          | created
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
     integer(ip)                                  :: state  = start
     class(fe_space_t)              , pointer     :: fe_space               => NULL()
     type(p_discrete_integration_t) , allocatable :: approximations(:)      
     class(matrix_array_assembler_t), pointer     :: matrix_array_assembler => NULL()
   contains
     procedure          :: create                      => fe_affine_operator_create
     procedure          :: symbolic_setup              => fe_affine_operator_symbolic_setup
     procedure          :: numerical_setup             => fe_affine_operator_numerical_setup
     procedure          :: apply                       => fe_affine_operator_apply
     procedure          :: is_linear                   => fe_affine_operator_is_linear
     procedure          :: get_tangent                 => fe_affine_operator_get_tangent
     procedure          :: get_translation             => fe_affine_operator_get_translation
     procedure          :: get_matrix                  => fe_affine_operator_get_matrix
     procedure          :: get_array                   => fe_affine_operator_get_array
     procedure          :: free_in_stages              => fe_affine_operator_free_in_stages
     procedure          :: free                        => fe_affine_operator_free
     procedure          :: get_domain_vector_space     => fe_affine_operator_get_domain_vector_space
     procedure          :: get_range_vector_space      => fe_affine_operator_get_range_vector_space
     procedure          :: abort_if_not_in_range       => fe_affine_operator_abort_if_not_in_range
     procedure          :: abort_if_not_in_domain      => fe_affine_operator_abort_if_not_in_domain
     procedure, private :: fe_affine_operator_free_values
     procedure, private :: fe_affine_operator_free_struct
     procedure, private :: fe_affine_operator_free_clean
     procedure, private :: fe_affine_operator_setup
  end type fe_affine_operator_t

  ! Types
  public :: fe_affine_operator_t

contains
  subroutine fe_affine_operator_create (this, &
                                        diagonal_blocks_symmetric_storage,&
                                        diagonal_blocks_symmetric,&
                                        diagonal_blocks_sign,&
                                        fe_space,&
                                        approximations)
    implicit none
    class(fe_affine_operator_t)            , intent(out) :: this
    class(fe_space_t)              , target, intent(in)  :: fe_space
    logical                                , intent(in)  :: diagonal_blocks_symmetric_storage(fe_space%dof_descriptor%nblocks)
    logical                                , intent(in)  :: diagonal_blocks_symmetric(fe_space%dof_descriptor%nblocks)
    integer(ip)                            , intent(in)  :: diagonal_blocks_sign(fe_space%dof_descriptor%nblocks)
    type(p_discrete_integration_t)         , intent(in)  :: approximations(:)
    
    integer(ip) :: iapprox

    assert(this%state == start)

    this%fe_space               => fe_space
    allocate ( this%approximations(size(approximations)) )
    do iapprox=1, size(approximations)
      this%approximations(iapprox)%discrete_integration => approximations(iapprox)%discrete_integration
    end do
    this%matrix_array_assembler => fe_space%create_matrix_array_assembler(diagonal_blocks_symmetric_storage, &
                                                                          diagonal_blocks_symmetric, &
                                                                          diagonal_blocks_sign)
    this%state = created
  end subroutine fe_affine_operator_create
  
  subroutine fe_affine_operator_symbolic_setup (this)
    implicit none
    class(fe_affine_operator_t), intent(inout) :: this

    assert ( .not. this%state == start )
    
    if ( this%state == created ) then 
      call this%fe_space%symbolic_setup_matrix_array_assembler(this%matrix_array_assembler)
      this%state = symbolically_setup
    end if

  end subroutine fe_affine_operator_symbolic_setup
  
  subroutine fe_affine_operator_numerical_setup (this)
    implicit none
    class(fe_affine_operator_t), intent(inout) :: this
    class(matrix_t), pointer :: matrix

    type(vector_space_t), pointer :: fe_affine_operator_domain_vector_space
    type(vector_space_t), pointer :: fe_affine_operator_range_vector_space
    type(vector_space_t), pointer :: matrix_domain_vector_space
    type(vector_space_t), pointer :: matrix_range_vector_space

    assert ( .not. this%state == start )

    if ( this%state == created ) then
       call this%symbolic_setup()
    end if

    if ( this%state == symbolically_setup ) then
       call this%matrix_array_assembler%allocate()
       call this%fe_space%volume_integral(this%approximations,this%matrix_array_assembler)
       matrix => this%matrix_array_assembler%get_matrix()
       matrix_domain_vector_space => matrix%get_domain_vector_space()
       matrix_range_vector_space => matrix%get_range_vector_space()
       fe_affine_operator_domain_vector_space => operator_get_domain_vector_space(this)
       fe_affine_operator_range_vector_space => operator_get_range_vector_space(this)
       call matrix_domain_vector_space%clone(fe_affine_operator_domain_vector_space)
       call matrix_range_vector_space%clone(fe_affine_operator_range_vector_space)
       this%state = numerically_setup
    end if

  end subroutine fe_affine_operator_numerical_setup

  subroutine fe_affine_operator_free_values(this)
    implicit none
    class(fe_affine_operator_t), intent(inout) :: this
    call this%matrix_array_assembler%free_in_stages(free_values)
    call this%free_vector_spaces()
  end subroutine fe_affine_operator_free_values

  subroutine fe_affine_operator_free_struct(this)
    implicit none
    class(fe_affine_operator_t), intent(inout) :: this
    call this%matrix_array_assembler%free_in_stages(free_struct)
  end subroutine fe_affine_operator_free_struct

  subroutine fe_affine_operator_free_clean(this)
    implicit none
    class(fe_affine_operator_t), intent(inout) :: this
    integer(ip) :: istat
    call this%matrix_array_assembler%free_in_stages(free_clean)
    deallocate ( this%matrix_array_assembler, stat=istat )
    check(istat==0)
    nullify(this%matrix_array_assembler)
    nullify(this%fe_space)
    deallocate(this%approximations, stat=istat)
    check(istat==0)
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
       if ( action == free_struct ) then
          call this%fe_affine_operator_free_struct()
          this%state = created
       else if ( action == free_clean ) then
          call this%fe_affine_operator_free_struct()
          call this%fe_affine_operator_free_clean()
          this%state=start
       end if
    else if ( this%state == numerically_setup ) then
       if ( action == free_values ) then
          call this%fe_affine_operator_free_values()
          this%state = symbolically_setup
       else if ( action == free_struct ) then
          call this%fe_affine_operator_free_values()
          call this%fe_affine_operator_free_struct()
          this%state = created
       else if ( action == free_clean ) then
          call this%fe_affine_operator_free_values()
          call this%fe_affine_operator_free_struct()
          call this%fe_affine_operator_free_clean()
          this%state=start
       end if
    end if
    
  end subroutine fe_affine_operator_free_in_stages
  
  subroutine fe_affine_operator_free(this)
    implicit none
    class(fe_affine_operator_t), intent(inout) :: this
    call this%free_in_stages(free_values)
    call this%free_in_stages(free_struct)
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
  
  ! op%apply(x,y) <=> y <- op*x
  ! Implicitly assumes that y is already allocated
  subroutine fe_affine_operator_apply(op,x,y) 
    implicit none
    class(fe_affine_operator_t), intent(in)    :: op
    class(vector_t) , intent(in)    :: x
    class(vector_t) , intent(inout) :: y 
    class(matrix_t) , pointer       :: matrix
    class(array_t)  , pointer       :: array
    call op%fe_affine_operator_setup()
    call op%abort_if_not_in_domain(x)
    call op%abort_if_not_in_range(y)
    call x%GuardTemp()
    matrix => op%matrix_array_assembler%get_matrix()
    call matrix%apply(x,y)
    array => op%matrix_array_assembler%get_array()
    call y%axpby( -1.0_rp, array, 1.0_rp )
    call x%CleanTemp()
  end subroutine fe_affine_operator_apply
  
  function fe_affine_operator_is_linear(op)
    implicit none
    class(fe_affine_operator_t), intent(in) :: op
    logical :: fe_affine_operator_is_linear
    fe_affine_operator_is_linear = .false.
  end function fe_affine_operator_is_linear
  
  function fe_affine_operator_get_tangent(op,x) result(tangent)
    implicit none
    class(fe_affine_operator_t), intent(in) :: op
    class(vector_t), optional  , intent(in) :: x
    type(dynamic_state_operator_t)          :: tangent
    call op%fe_affine_operator_setup()
    tangent = op%get_matrix()
    call tangent%SetTemp()
  end function fe_affine_operator_get_tangent
  
  function fe_affine_operator_get_translation(op) result(translation)
    implicit none
    class(fe_affine_operator_t), intent(in) :: op
    class(vector_t), pointer                :: translation
    call op%fe_affine_operator_setup()
    translation => op%get_array()
  end function fe_affine_operator_get_translation

  subroutine fe_affine_operator_abort_if_not_in_domain ( this, vector )
    implicit none
    class(fe_affine_operator_t), intent(in)  :: this
    class(vector_t)            , intent(in)  :: vector
    call this%fe_affine_operator_setup()
    call operator_abort_if_not_in_domain(this,vector)
  end subroutine fe_affine_operator_abort_if_not_in_domain

  subroutine fe_affine_operator_abort_if_not_in_range ( this, vector )
    implicit none
    class(fe_affine_operator_t), intent(in) :: this
    class(vector_t)            , intent(in) :: vector
    call this%fe_affine_operator_setup()
    call operator_abort_if_not_in_range(this,vector)
  end subroutine fe_affine_operator_abort_if_not_in_range

  function fe_affine_operator_get_domain_vector_space ( this )
    implicit none
    class(fe_affine_operator_t), target, intent(in) :: this
    type(vector_space_t)               , pointer    :: fe_affine_operator_get_domain_vector_space
    call this%fe_affine_operator_setup()
    fe_affine_operator_get_domain_vector_space => operator_get_domain_vector_space(this)
  end function fe_affine_operator_get_domain_vector_space
  
  function fe_affine_operator_get_range_vector_space ( this )
    implicit none
    class(fe_affine_operator_t), target, intent(in) :: this
    type(vector_space_t)                  , pointer :: fe_affine_operator_get_range_vector_space
    call this%fe_affine_operator_setup()
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
    call this%symbolic_setup()
    call this%numerical_setup()
  end subroutine fe_affine_operator_setup

end module fe_affine_operator_names
