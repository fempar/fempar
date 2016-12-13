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
program test_geometry

use fempar_names
use test_geometry_params
use gid_geometry_reader_names

implicit none
    type(geometry_t)      :: geometry
    type(line_t), pointer :: line
    type(point_t)         :: point
    real(rp)              :: param
    !real(rp)              :: point(3),param

    character(len=str_cla_len)       :: dir_path, dir_path_out
    character(len=str_cla_len)       :: prefix, filename
    integer(ip)                      :: lunio, istat
    type(parameterlist_t), pointer   :: parameterlist
    type(test_geometry_params_t)     :: test_params

    call FEMPAR_INIT()

    call test_params%create()
    parameterlist => test_params%get_parameters()
!    call geometry%read(parameterlist)
!    call geometry%free()

    !-----------------------------------------------------------------
    !< Create a hexa
    !<     7-------8
    !<    /       /|
    !<   /       / |
    !<  /       /  |
    !< 5---3---6___4       Z   Y
    !< |  /    |  /        |  /
    !  | /     | /         | /
    !< |/      |/          |/ 
    !< 1-------2           +----- X
    !-----------------------------------------------------------------
    call geometry%create(num_points=17, num_lines=18, num_surfaces=7, num_volumes=1)

    call geometry%add_point(point_id=1, coord=[0._rp,0._rp,0._rp])

    call geometry%add_line(line_id=1, coord1=[0._rp,0._rp,0._rp], &
                                      coord2=[1._rp,0._rp,0._rp])

    call geometry%add_line(line_id=2, coord1=[0._rp,0._rp,0._rp], &
                                      coord2=[1._rp,0._rp,0._rp])

    call geometry%add_quad(coord1=[0._rp,0._rp,0._rp], &
                           coord2=[1._rp,0._rp,0._rp], &
                           coord3=[0._rp,0._rp,1._rp], &
                           coord4=[1._rp,0._rp,1._rp])

    call geometry%add_hexa(coord1=[0._rp,0._rp,0._rp], &
                           coord2=[1._rp,0._rp,0._rp], &
                           coord3=[0._rp,0._rp,1._rp], &
                           coord4=[1._rp,0._rp,1._rp], &
                           coord5=[0._rp,1._rp,0._rp], &
                           coord6=[1._rp,1._rp,0._rp], &
                           coord7=[0._rp,1._rp,1._rp], &
                           coord8=[1._rp,1._rp,1._rp])


    call geometry%init()
    call geometry%free()

    
    call the_gid_geometry_reader%fill_geometry(parameterlist, geometry)
    call geometry%free()

    call FEMPAR_FINALIZE()

end program test_geometry
