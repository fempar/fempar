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
program par_test_element_exchange
  !----------------------------------------------------------
  ! Parallel partitioner test
  !----------------------------------------------------------
  use fem
  use par
  implicit none
#include "debug.i90"
  ! Our data
  type(par_context)        :: context
  type(par_partition)      :: p_part
  type(par_mesh)           :: p_mesh
  type(par_triangulation)  :: p_trian

  ! Arguments
  integer(ip)              :: handler
  character(len=256)       :: dir_path, dir_path_out
  character(len=256)       :: prefix
  integer(ip)              :: i, j

  call meminit

  ! Start parallel execution
  handler = inhouse
  call par_context_create (handler, context)

  ! Read parameters from command-line
  call  read_pars_cl_par_test_element_exchange ( dir_path, prefix, dir_path_out )

  !if ( context%iam > 0 ) then
  !   stop
  !end if
  
  write(*,*) ' KK PROC ', context%iam
  
  ! Read partition info. Associate contexts
  call par_partition_create ( dir_path, prefix, context, p_part )

  write(*,*) ' KK PROC 2 ', context%iam
  ! Read mesh
  call par_mesh_create ( dir_path, prefix, p_part, p_mesh )

  call par_mesh_to_triangulation (p_mesh, p_trian)

  ! do i=1,p_trian%num_elems + p_trian%num_ghosts
  !   write(*,'(10i10)') p_trian%elems(i)%objects_GIDs
  ! end do
  ! do i=1,p_trian%num_elems + p_trian%num_ghosts
  !    write(*,'(10i10)') p_trian%f_trian%elems(i)%objects
  ! end do

  call par_triangulation_free(p_trian)
  call par_mesh_free (p_mesh)
  call par_partition_free (p_part)
  call par_context_free ( context )

  call memstatus

contains
  subroutine read_pars_cl_par_test_element_exchange (dir_path, prefix, dir_path_out)
    implicit none
    character*(*), intent(out)   :: dir_path, prefix, dir_path_out
    character(len=256)           :: program_name
    character(len=256)           :: argument 
    integer                      :: numargs,iargc

    numargs = iargc()
    call getarg(0, program_name)
    if (.not. (numargs==3) ) then
       write (6,*) 'Usage: ', trim(program_name), ' dir_path prefix dir_path_out'
       stop
    end if

    call getarg(1, argument)
    dir_path = trim(argument)

    call getarg(2, argument)
    prefix = trim(argument)

    call getarg(3,argument)
    dir_path_out = trim(argument)

  end subroutine read_pars_cl_par_test_element_exchange

end program par_test_element_exchange
