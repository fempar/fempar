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

module direct_solver_creational_methods_dictionary_names

USE direct_solver_parameters_names,  only: pardiso_mkl, umfpack
USE ProcedureDictionary
USE base_direct_solver_names,        only: create_direct_solver_interface
USE pardiso_mkl_direct_solver_names, only: create_pardiso_mkl_direct_solver
USE umfpack_direct_solver_names,     only: create_umfpack_direct_solver

implicit none
private

    type :: direct_solver_creational_methods_dictionary_t
    private
        type(ProcedureDictionary_t) :: creational_methods
    contains
    private
        procedure, non_overridable, public :: Init          => direct_solver_creational_methods_dictionary_Init
        procedure, non_overridable, public :: isInitialized => direct_solver_creational_methods_dictionary_isInitialized
        procedure, non_overridable, public :: Set           => direct_solver_creational_methods_dictionary_Set
        procedure, non_overridable, public :: Get           => direct_solver_creational_methods_dictionary_Get
        procedure, non_overridable, public :: isPresent     => direct_solver_creational_methods_dictionary_isPresent
        procedure, non_overridable, public :: Free          => direct_solver_creational_methods_dictionary_Free
    end type

    type(direct_solver_creational_methods_dictionary_t) :: The_direct_solver_creational_methods_dictionary

public :: The_direct_solver_creational_methods_dictionary

contains
    subroutine direct_solver_creational_methods_dictionary_Init(this,Size)
    !-----------------------------------------------------------------
        class(direct_solver_creational_methods_dictionary_t), intent(INOUT) :: this
        integer, optional,                                intent(IN)    :: Size
    !-----------------------------------------------------------------
        call this%creational_methods%Free()
        call this%creational_methods%Init(Size = Size)
        call this%Set(Key=pardiso_mkl, Proc=create_pardiso_mkl_direct_solver)
        call this%Set(Key=umfpack,     Proc=create_umfpack_direct_solver)
    end subroutine

    function direct_solver_creational_methods_dictionary_isInitialized(this) result(isInitialized)
    !-----------------------------------------------------------------
        class(direct_solver_creational_methods_dictionary_t), intent(INOUT) :: this
        logical                                                         :: isInitialized
    !-----------------------------------------------------------------
        isInitialized = this%creational_methods%isInitialized()
    end function

    subroutine direct_solver_creational_methods_dictionary_Set(this,Key,Proc)
    !-----------------------------------------------------------------
        class(direct_solver_creational_methods_dictionary_t), intent(INOUT) :: this
        character(len=*),                                 intent(IN)    :: Key
        procedure(create_direct_solver_interface)                       :: Proc
        procedure(create_direct_solver_interface), pointer              :: ProcPointer
    !-----------------------------------------------------------------
        ProcPointer => Proc
        call this%creational_methods%Set(Key=Key,Proc=ProcPointer)
    end subroutine

    subroutine direct_solver_creational_methods_dictionary_Get(this,Key,Proc)
    !-----------------------------------------------------------------
        class(direct_solver_creational_methods_dictionary_t),  intent(IN) :: this
        character(len=*),                                  intent(IN) :: Key
        procedure(create_direct_solver_interface), pointer            :: Proc
    !-----------------------------------------------------------------
        call this%creational_methods%Get(Key=Key,Proc=Proc)
    end subroutine

    function direct_solver_creational_methods_dictionary_isPresent(this, Key) result(isPresent)
    !-----------------------------------------------------------------
        class(direct_solver_creational_methods_dictionary_t), intent(INOUT) :: this
        character(len=*),                                  intent(IN) :: Key
        logical                                                         :: isPresent
    !-----------------------------------------------------------------
        isPresent = this%creational_methods%isPresent(Key=Key)
    end function

    subroutine direct_solver_creational_methods_dictionary_Free(this)
    !-----------------------------------------------------------------
        class(direct_solver_creational_methods_dictionary_t), intent(INOUT) :: this
    !-----------------------------------------------------------------
        call this%creational_methods%Free()
    end subroutine
end module direct_solver_creational_methods_dictionary_names

