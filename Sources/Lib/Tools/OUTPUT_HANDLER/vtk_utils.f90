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

module vtk_utils

USE IR_Precision, only: str
USE types_names

implicit none
#include "debug.i90"
private

    ! File extensions and time prefix
    character(len=5), parameter :: time_prefix = 'time_'
    character(len=4), parameter :: vtk_ext     = '.vtu'
    character(len=4), parameter :: pvd_ext     = '.pvd'
    character(len=5), parameter :: pvtu_ext    = '.pvtu'


contains

    function get_vtk_output_directory(dir_path, time_step) result(path)
    !-----------------------------------------------------------------
    !< Build time output dir path for the vtk files in each timestep
    !-----------------------------------------------------------------
        character(len=*),  intent(in)    :: dir_path
        real(rp),          intent(in)    :: time_step
        character(len=:), allocatable    :: path
        character(len=:), allocatable    :: fp
    !-----------------------------------------------------------------
        path = trim(adjustl(dir_path))//'/'//time_prefix//trim(adjustl(str(no_sign=.true., n=time_step)))//'/'
    end function get_vtk_output_directory


    function get_pvd_output_directory(dir_path, time_step) result(path)
    !-----------------------------------------------------------------
    !< Build output dir path for the PVD files
    !-----------------------------------------------------------------
        character(len=*), intent(in)    :: dir_path
        real(RP),         intent(in)    :: time_step
        character(len=:), allocatable   :: path
    !-----------------------------------------------------------------
        path = time_prefix//trim(adjustl(str(no_sign=.true., n=time_step)))//'/'
    end function get_pvd_output_directory


    function get_vtk_filename(prefix, part) result(filename)
    !-----------------------------------------------------------------
    !< Build VTK filename
    !-----------------------------------------------------------------
        character(len=*),  intent(in) :: prefix
        integer(ip),       intent(in) :: part
        character(len=:), allocatable :: filename
    !-----------------------------------------------------------------
        filename = trim(adjustl(prefix))//'_'//trim(adjustl(str(no_sign=.true., n=part)))//vtk_ext
    end function get_vtk_filename


    function get_pvtu_filename(prefix, time_step) result(filename)
    !-----------------------------------------------------------------
    !< Build pvtu filename
    !-----------------------------------------------------------------
        character(len=*),   intent(in) :: prefix
        real(rp),           intent(in) :: time_step
        character(len=:), allocatable  :: filename
    !-----------------------------------------------------------------
        filename = trim(adjustl(prefix))//'_'//trim(adjustl(str(no_sign=.true., n=time_step)))//pvtu_ext
    end function get_pvtu_filename

end module vtk_utils

