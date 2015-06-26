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
module par_timer_names
  ! Serial modules
  use types_names

  ! Parallel modules
  use psb_penv_mod_names

# include "debug.i90"
  implicit none
  private

  integer(ip), parameter :: max_message = 256
  ! integer(ip), parameter :: fmt_message =  20

  type par_timer_t
      integer                :: ictxt    ! Parallel context
      integer                :: root     ! PID responsible of gathering/reporting timings
      character(max_message) :: message  ! Concept being measured (e.g., assembly)
      real(8)                :: start    ! last call to start
      real(8)                :: stop     ! last call to stop
      real(8)                :: accum    ! sum of all stop-start 
  end type par_timer_t

  ! Public types
  public :: par_timer_t 

  ! Public routines
  public :: par_timer_create, par_timer_init, & 
            par_timer_start,  par_timer_stop, &
            par_timer_report
contains
  
    subroutine par_timer_create ( p_timer, message, ictxt, root )
#ifdef MPI_MOD
      use mpi
#endif
      implicit none 
#ifdef MPI_H
      include 'mpif.h'
#endif
      ! Parameters
      type(par_timer_t), intent(out)          :: p_timer
      character*(*)  , intent(in)           :: message
      integer        , intent(in), optional :: ictxt
      integer        , intent(in), optional :: root
      
      ! Locals
      integer                               :: ictxt_, root_

      ictxt_ = mpi_comm_world
      root_  = 0 
      if ( present(ictxt) ) ictxt_ = ictxt
      if ( present(root)  ) root_  = root 
      
      p_timer%ictxt   = ictxt_
      p_timer%root    = root_ 
      p_timer%message = trim(message)
      p_timer%start   = 0.0 
      p_timer%stop    = 0.0
      p_timer%accum   = 0.0
      p_timer%accum   = 1.79769E+308 ! Largest double precision number
    end subroutine par_timer_create

    subroutine par_timer_init ( p_timer )
      implicit none 
      ! Parameters
      type(par_timer_t), intent(inout)          :: p_timer
      p_timer%start  = 0.0 
      p_timer%stop   = 0.0
      p_timer%accum  = 0.0
      p_timer%accum  = 1.79769E+308 ! Largest double precision number
    end subroutine par_timer_init

    subroutine par_timer_start ( p_timer )
      implicit none 
      ! Parameters
      type(par_timer_t), intent(inout)          :: p_timer
      call psb_barrier (p_timer%ictxt)
      p_timer%start  = psb_wtime()
    end subroutine par_timer_start

    subroutine par_timer_stop ( p_timer )
      implicit none 
      ! Parameters
      type(par_timer_t), intent(inout)          :: p_timer
      p_timer%stop = psb_wtime()

      if ( p_timer%stop - p_timer%start >= 0.0) then
         ! p_timer%accum = p_timer%accum + (p_timer%stop - p_timer%start)
         if ( p_timer%accum >  (p_timer%stop - p_timer%start) ) p_timer%accum = (p_timer%stop - p_timer%start)
         ! p_timer%accum = (p_timer%stop - p_timer%start)
      end if 
      p_timer%start  = 0.0 
      p_timer%stop   = 0.0

    end subroutine par_timer_stop
   
    subroutine par_timer_report ( p_timer, show_header, luout )
      implicit none 
      ! Parameters
      type(par_timer_t), intent(inout)     :: p_timer
      logical, intent(in), optional      :: show_header 
      integer(ip), intent(in), optional  :: luout
      
      ! Locals
      character(len=*), parameter    :: fmt_header = '(a25,1x,3(2x,a15),3(2x,a15))'
      character(len=*), parameter    :: fmt_data   = '(a25,1x,3(2x,es15.9),3(2x,es15.9))'
      real(8)                        :: accum_max, accum_min, accum_sum
      integer                        :: my_id, num_procs
      logical                        :: show_header_

      call psb_info ( p_timer%ictxt, my_id, num_procs )

      accum_max = p_timer%accum
      accum_min = p_timer%accum
      accum_sum = p_timer%accum
      show_header_ = .true.
      if (present(show_header)) show_header_ = show_header 

      if ( show_header_ ) then 
         if (present(luout) ) then
            if ( my_id == p_timer%root ) write(luout,fmt_header) '', 'Min (secs.)', 'Max (secs.)', 'Avg (secs.)'
         else
            if ( my_id == p_timer%root ) write(*,fmt_header) '', 'Min (secs.)', 'Max (secs.)', 'Avg (secs.)'
         end if
      end if         
  
      call psb_max ( p_timer%ictxt, accum_max, p_timer%root )
      call psb_min ( p_timer%ictxt, accum_min, p_timer%root )
      call psb_sum ( p_timer%ictxt, accum_sum, p_timer%root )
      
      if (present(luout) ) then
         if ( my_id == p_timer%root ) write(luout,fmt_data) adjustl(p_timer%message), accum_min, accum_max, accum_sum/num_procs
      else
         if ( my_id == p_timer%root ) write(*,fmt_data) adjustl(p_timer%message), accum_min, accum_max, accum_sum/num_procs
      end if

    end subroutine par_timer_report

end module par_timer_names
