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

module parameters_consistency_names

use ParameterList
use types_names

implicit none
#include "debug.i90"
private

    interface parameter_consistency
        module procedure parameter_consistency_0D
        module procedure parameter_consistency_1D
        module procedure parameter_consistency_2D
        module procedure parameter_consistency_3D
        module procedure parameter_consistency_4D
        module procedure parameter_consistency_5D
        module procedure parameter_consistency_6D
        module procedure parameter_consistency_7D
    end interface

public :: parameter_consistency

contains


   function parameter_consistency_0D (parameter_list, key, value ) result(is_consistent)
        implicit none
        type(ParameterList_t), intent(in) :: parameter_list
        character(len=*),      intent(in) :: key
        class(*),              intent(in) :: value
        logical                           :: is_consistent
#ifdef DEBUG
        integer                           :: error
        integer(ip), allocatable          :: parameter_shape(:)
#endif
        is_consistent = .true.
#ifdef DEBUG
        is_consistent   = parameter_list%isPresent   (Key = key)
        if(is_consistent) then
            is_consistent  = (parameter_list%isOfDataType(Key = key, mold  = value) .and. &
                             (parameter_list%GetShape(Key = key, shape = parameter_shape) == 0))
            is_consistent  = (is_consistent .and. size(parameter_shape) == 0)
        endif 
#endif
   end function parameter_consistency_0D


   function parameter_consistency_1D (parameter_list, key, value ) result(is_consistent)
        implicit none
        type(ParameterList_t), intent(in) :: parameter_list
        character(len=*),      intent(in) :: key
        class(*),              intent(in) :: value(:)
        logical                           :: is_consistent
#ifdef DEBUG
        integer                           :: error
        integer(ip), allocatable          :: parameter_shape(:)
#endif
        is_consistent = .true.
#ifdef DEBUG
        is_consistent   = parameter_list%isPresent   (Key = key)
        if(is_consistent) then
            is_consistent  = (parameter_list%isOfDataType(Key = key, mold  = value) .and. &
                             (parameter_list%GetShape(Key = key, shape = parameter_shape) == 0))
            is_consistent  = (is_consistent .and. all(parameter_shape == shape(Value)))
        endif 
#endif
   end function parameter_consistency_1D


   function parameter_consistency_2D (parameter_list, key, value ) result(is_consistent)
        implicit none
        type(ParameterList_t), intent(in) :: parameter_list
        character(len=*),      intent(in) :: key
        class(*),              intent(in) :: value(:,:)
        logical                           :: is_consistent
#ifdef DEBUG
        integer                           :: error
        integer(ip), allocatable          :: parameter_shape(:)
#endif
        is_consistent = .true.
#ifdef DEBUG
        is_consistent   = parameter_list%isPresent   (Key = key)
        if(is_consistent) then
            is_consistent  = (parameter_list%isOfDataType(Key = key, mold  = value) .and. &
                             (parameter_list%GetShape(Key = key, shape = parameter_shape) == 0))
            is_consistent  = (is_consistent .and. all(parameter_shape == shape(Value)))
        endif 
#endif
   end function parameter_consistency_2D


   function parameter_consistency_3D (parameter_list, key, value ) result(is_consistent)
        implicit none
        type(ParameterList_t), intent(in) :: parameter_list
        character(len=*),      intent(in) :: key
        class(*),              intent(in) :: value(:,:,:)
        logical                           :: is_consistent
#ifdef DEBUG
        integer                           :: error
        integer(ip), allocatable          :: parameter_shape(:)
#endif
        is_consistent = .true.
#ifdef DEBUG
        is_consistent   = parameter_list%isPresent   (Key = key)
        if(is_consistent) then
            is_consistent  = (parameter_list%isOfDataType(Key = key, mold  = value) .and. &
                             (parameter_list%GetShape(Key = key, shape = parameter_shape) == 0))
            is_consistent  = (is_consistent .and. all(parameter_shape == shape(Value)))
        endif 
#endif
   end function parameter_consistency_3D


   function parameter_consistency_4D (parameter_list, key, value ) result(is_consistent)
        implicit none
        type(ParameterList_t), intent(in) :: parameter_list
        character(len=*),      intent(in) :: key
        class(*),              intent(in) :: value(:,:,:,:)
        logical                           :: is_consistent
#ifdef DEBUG
        integer                           :: error
        integer(ip), allocatable          :: parameter_shape(:)
#endif
        is_consistent = .true.
#ifdef DEBUG
        is_consistent   = parameter_list%isPresent   (Key = key)
        if(is_consistent) then
            is_consistent  = (parameter_list%isOfDataType(Key = key, mold  = value) .and. &
                             (parameter_list%GetShape(Key = key, shape = parameter_shape) == 0))
            is_consistent  = (is_consistent .and. all(parameter_shape == shape(Value)))
        endif 
#endif
   end function parameter_consistency_4D


   function parameter_consistency_5D (parameter_list, key, value ) result(is_consistent)
        implicit none
        type(ParameterList_t), intent(in) :: parameter_list
        character(len=*),      intent(in) :: key
        class(*),              intent(in) :: value(:,:,:,:,:)
        logical                           :: is_consistent
#ifdef DEBUG
        integer                           :: error
        integer(ip), allocatable          :: parameter_shape(:)
#endif
        is_consistent = .true.
#ifdef DEBUG
        is_consistent   = parameter_list%isPresent   (Key = key)
        if(is_consistent) then
            is_consistent  = (parameter_list%isOfDataType(Key = key, mold  = value) .and. &
                             (parameter_list%GetShape(Key = key, shape = parameter_shape) == 0))
            is_consistent  = (is_consistent .and. all(parameter_shape == shape(Value)))
        endif 
#endif
   end function parameter_consistency_5D


   function parameter_consistency_6D (parameter_list, key, value ) result(is_consistent)
        implicit none
        type(ParameterList_t), intent(in) :: parameter_list
        character(len=*),      intent(in) :: key
        class(*),              intent(in) :: value(:,:,:,:,:,:)
        logical                           :: is_consistent
#ifdef DEBUG
        integer                           :: error
        integer(ip), allocatable          :: parameter_shape(:)
#endif
        is_consistent = .true.
#ifdef DEBUG
        is_consistent   = parameter_list%isPresent   (Key = key)
        if(is_consistent) then
            is_consistent  = (parameter_list%isOfDataType(Key = key, mold  = value) .and. &
                             (parameter_list%GetShape(Key = key, shape = parameter_shape) == 0))
            is_consistent  = (is_consistent .and. all(parameter_shape == shape(Value)))
        endif 
#endif
   end function parameter_consistency_6D


   function parameter_consistency_7D (parameter_list, key, value ) result(is_consistent)
        implicit none
        type(ParameterList_t), intent(in) :: parameter_list
        character(len=*),      intent(in) :: key
        class(*),              intent(in) :: value(:,:,:,:,:,:,:)
        logical                           :: is_consistent
#ifdef DEBUG
        integer                           :: error
        integer(ip), allocatable          :: parameter_shape(:)
#endif
        is_consistent = .true.
#ifdef DEBUG
        is_consistent   = parameter_list%isPresent   (Key = key)
        if(is_consistent) then
            is_consistent  = (parameter_list%isOfDataType(Key = key, mold  = value) .and. &
                             (parameter_list%GetShape(Key = key, shape = parameter_shape) == 0))
            is_consistent  = (is_consistent .and. all(parameter_shape == shape(Value)))
        endif 
#endif
   end function parameter_consistency_7D

end module parameters_consistency_names
