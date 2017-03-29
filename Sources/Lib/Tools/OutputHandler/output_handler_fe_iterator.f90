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
!---------------------------------------------------------------------
!* Author: Víctor Sande Veiga
! Date: 2016-11-29
! Version: 0.0.1
! Category: IO
!
!--------------------------------------------------------------------- 
!### Iterator over [[serial_fe_space_t(type)]]
!
! Contains the following public entities:
! [[output_handler_fe_iterator_names(module)]]
!
!---------------------------------------------------------------------
module output_handler_fe_iterator_names
!---------------------------------------------------------------------
!* Author: Víctor Sande Veiga
! Date: 2016-11-29
! Version: 0.0.1
! Category: IO
!
!--------------------------------------------------------------------- 
!### Iterator over [[serial_fe_space_t(type)]]
!
! Parent class for all custom iterator over [[serial_fe_space_t(type)]]
!
! Contains the following public entities:
! [[output_handler_fe_iterator_t(type)]]
!
!---------------------------------------------------------------------

use fe_space_names, only: fe_accessor_t, serial_fe_space_t

implicit none
#include "debug.i90"
private

    type output_handler_fe_iterator_t
    !-----------------------------------------------------------------
    !* Author: Víctor Sande Veiga
    ! Date: 2016-11-29
    ! Version: 0.0.1
    ! Category: IO
    !
    !-----------------------------------------------------------------
    !### Iterator over the finite elements in [[serial_fe_space_t(type)]]
    !
    ! Extending this class, allows the definition of custom loops
    ! over the finite elements in [[serial_fe_space_t(type)]]
    !-----------------------------------------------------------------
    private
        !type(fe_iterator_t) :: fe_iterator
        class(fe_accessor_t), allocatable :: fe_accessor
    contains
    private
        procedure, public :: create          => output_handler_fe_iterator_create
        !procedure, public :: get_fe_iterator => output_handler_fe_iterator_get_fe_iterator
        procedure, public :: free            => output_handler_fe_iterator_free
        procedure, public :: init            => output_handler_fe_iterator_init
        procedure, public :: next            => output_handler_fe_iterator_next
        procedure, public :: has_finished    => output_handler_fe_iterator_has_finished
        procedure         ::                    output_handler_fe_iterator_current
        generic,   public :: current         => output_handler_fe_iterator_current
  end type output_handler_fe_iterator_t

public :: output_handler_fe_iterator_t

contains


    subroutine output_handler_fe_iterator_create(this, fe_space)
    !-----------------------------------------------------------------
    !< Create the iterator from a given [[serial_fe_space_t(type)]]
    !-----------------------------------------------------------------
        class(output_handler_fe_iterator_t), intent(inout) :: this
        class(serial_fe_space_t)           , intent(in)    :: fe_space
    !-----------------------------------------------------------------
        call this%free()
        !this%fe_iterator = fe_space%create_fe_iterator()
        call fe_space%create_fe_accessor(this%fe_accessor)
    end subroutine output_handler_fe_iterator_create


    !function output_handler_fe_iterator_get_fe_iterator(this) result(fe_iterator)
    !!-----------------------------------------------------------------
    !!< Return a pointer to [[output_handler_fe_iterator_t(type)]]
    !!-----------------------------------------------------------------
    !    class(output_handler_fe_iterator_t), target, intent(in) :: this
    !    type(fe_iterator_t), pointer                            :: fe_iterator
    !!-----------------------------------------------------------------
    !    fe_iterator => this%fe_iterator 
    !end function output_handler_fe_iterator_get_fe_iterator


    subroutine output_handler_fe_iterator_free(this)
    !-----------------------------------------------------------------
    !< Free procedure
    !-----------------------------------------------------------------
        class(output_handler_fe_iterator_t), intent(inout) :: this
    !-----------------------------------------------------------------
        call this%fe_accessor%free()
    end subroutine output_handler_fe_iterator_free


    subroutine output_handler_fe_iterator_init(this)
    !-----------------------------------------------------------------
    !< Initialize the contained [[fe_iterator_t(type)]]
    !-----------------------------------------------------------------
        class(output_handler_fe_iterator_t), intent(inout) :: this
    !-----------------------------------------------------------------
        call this%fe_accessor%first()
    end subroutine output_handler_fe_iterator_init


    subroutine output_handler_fe_iterator_next(this)
    !-----------------------------------------------------------------
    !< Jump to the next position of [[fe_iterator_t(type)]]
    !-----------------------------------------------------------------
        class(output_handler_fe_iterator_t), intent(inout) :: this
    !-----------------------------------------------------------------
        call this%fe_accessor%next()
    end subroutine output_handler_fe_iterator_next


    function output_handler_fe_iterator_has_finished(this) result(has_finished)
    !-----------------------------------------------------------------
    !< Ask if [[fe_iterator_t(type)]] has reached the last position
    !-----------------------------------------------------------------
        class(output_handler_fe_iterator_t), intent(inout) :: this
        logical                                            :: has_finished
    !-----------------------------------------------------------------
        has_finished = this%fe_accessor%past_the_end()
    end function output_handler_fe_iterator_has_finished


    subroutine output_handler_fe_iterator_current(this, fe_accessor) 
    !-----------------------------------------------------------------
    !< Return the [[fe_accessor_t(type)]] of the current position
    !-----------------------------------------------------------------
        class(output_handler_fe_iterator_t), intent(in)    :: this
        type(fe_accessor_t),                 intent(inout) :: fe_accessor
    !-----------------------------------------------------------------
        !call this%fe_iterator%current(fe_accessor)
        fe_accessor = this%fe_accessor
    end subroutine output_handler_fe_iterator_current

end module output_handler_fe_iterator_names
