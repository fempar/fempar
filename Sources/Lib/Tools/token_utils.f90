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

module token_utils_names
    use types_names
    use flap_utils_m
    implicit none
# include "debug.i90"
    private

    public :: count_tokens_cla_string_list 
    public :: process_tokens_cla_string_list

contains

    function count_tokens_cla_string_list (string_list) 
    !------------------------------------------------------------------
    !< Count number of space separated tokens in a given string
    !------------------------------------------------------------------
        implicit none
        character(*),                     intent(in) :: string_list
        integer(ip)                                  :: count_tokens_cla_string_list
        character(len=len(string_list)), allocatable :: vals(:) 
    !------------------------------------------------------------------
        call tokenize(strin=string_list, delimiter=' ', toks=vals, Nt=count_tokens_cla_string_list)
        deallocate(vals)
    end function count_tokens_cla_string_list


    subroutine process_tokens_cla_string_list ( string_list, values )
    !------------------------------------------------------------------
    !< Return space sparated tokens in a given string
    !------------------------------------------------------------------ 
        implicit none 
        character(*),                  intent(in)    :: string_list
        character(*),                  intent(inout) :: values(:)
        character(len=len(string_list)), allocatable :: vals(:)
        integer(ip) :: Nt, i
    !------------------------------------------------------------------
        call tokenize(strin=string_list, delimiter=' ', toks=vals, Nt=Nt)
        assert ( size(values) == Nt) 
        do i=1, size(vals) 
            values(i)=vals(i)
        end do 
        deallocate(vals)
    end subroutine process_tokens_cla_string_list

end module token_utils_names 
