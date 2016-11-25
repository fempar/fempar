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

module iterative_linear_solver_creational_methods_dictionary_names

USE ProcedureDictionary
USE iterative_linear_solver_parameters_names,  only: cg_name,     &
                                                     fgmres_name, &
                                                     icg_name,    &
                                                     lfom_name,   &
                                                     lgmres_name, &
                                                     minres_name, &
                                                     rgmres_name, &
                                                     richardson_name
USE base_iterative_linear_solver_names,        only: create_iterative_linear_solver_interface
USE cg_names,                                  only: create_cg
USE fgmres_names,                              only: create_fgmres
USE icg_names,                                 only: create_icg
USE lfom_names,                                only: create_lfom
USE lgmres_names,                              only: create_lgmres
USE minres_names,                              only: create_minres
USE rgmres_names,                              only: create_rgmres
USE richardson_names,                          only: create_richardson

implicit none
private

    type :: iterative_linear_solver_creational_methods_dictionary_t
    private
        type(ProcedureDictionary_t) :: creational_methods
    contains
    private
        procedure, non_overridable, public :: Init          => creational_methods_dictionary_Init
        procedure, non_overridable, public :: isInitialized => creational_methods_dictionary_isInitialized
        procedure, non_overridable, public :: Set           => creational_methods_dictionary_Set
        procedure, non_overridable, public :: Get           => creational_methods_dictionary_Get
        procedure, non_overridable, public :: isPresent     => creational_methods_dictionary_isPresent
        procedure, non_overridable, public :: Free          => creational_methods_dictionary_Free
    end type

    type(iterative_linear_solver_creational_methods_dictionary_t), save :: The_iterative_linear_solver_creational_methods_dictionary

public :: The_iterative_linear_solver_creational_methods_dictionary

contains
    subroutine creational_methods_dictionary_Init(this,Size)
    !-----------------------------------------------------------------
        class(iterative_linear_solver_creational_methods_dictionary_t), intent(INOUT) :: this
        integer, optional,                                              intent(IN)    :: Size
    !-----------------------------------------------------------------
        call this%creational_methods%Free()
        call this%creational_methods%Init(Size = Size)
        call this%Set(Key=cg_name, Proc=create_cg)
        call this%Set(Key = fgmres_name,     Proc = create_fgmres)
        call this%Set(Key = icg_name,        Proc = create_icg)
        call this%Set(Key = lfom_name,       Proc = create_lfom)
        call this%Set(Key = lgmres_name,     Proc = create_lgmres)
        call this%Set(Key = minres_name,     Proc = create_minres)
        call this%Set(Key = rgmres_name,     Proc = create_rgmres)
        call this%Set(Key = richardson_name, Proc = create_richardson)
    end subroutine

    function creational_methods_dictionary_isInitialized(this) result(isInitialized)
    !-----------------------------------------------------------------
        class(iterative_linear_solver_creational_methods_dictionary_t), intent(INOUT) :: this
        logical                                                         :: isInitialized
    !-----------------------------------------------------------------
        isInitialized = this%creational_methods%isInitialized()
    end function

    subroutine creational_methods_dictionary_Set(this,Key,Proc)
    !-----------------------------------------------------------------
        class(iterative_linear_solver_creational_methods_dictionary_t), intent(INOUT) :: this
        character(len=*),                                 intent(IN)    :: Key
        procedure(create_iterative_linear_solver_interface)                       :: Proc
        procedure(create_iterative_linear_solver_interface), pointer              :: ProcPointer
    !-----------------------------------------------------------------
        ProcPointer => Proc
        call this%creational_methods%Set(Key=Key,Proc=ProcPointer)
    end subroutine

    subroutine creational_methods_dictionary_Get(this,Key,Proc)
    !-----------------------------------------------------------------
        class(iterative_linear_solver_creational_methods_dictionary_t),  intent(IN) :: this
        character(len=*),                                  intent(IN) :: Key
        procedure(create_iterative_linear_solver_interface), pointer            :: Proc
    !-----------------------------------------------------------------
        call this%creational_methods%Get(Key=Key,Proc=Proc)
    end subroutine

    function creational_methods_dictionary_isPresent(this, Key) result(isPresent)
    !-----------------------------------------------------------------
        class(iterative_linear_solver_creational_methods_dictionary_t), intent(INOUT) :: this
        character(len=*),                                  intent(IN) :: Key
        logical                                                         :: isPresent
    !-----------------------------------------------------------------
        isPresent = this%creational_methods%isPresent(Key=Key)
    end function

    subroutine creational_methods_dictionary_Free(this)
    !-----------------------------------------------------------------
        class(iterative_linear_solver_creational_methods_dictionary_t), intent(INOUT) :: this
    !-----------------------------------------------------------------
        call this%creational_methods%Free()
    end subroutine
end module iterative_linear_solver_creational_methods_dictionary_names


