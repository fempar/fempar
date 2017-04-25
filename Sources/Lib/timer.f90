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
module timer_names
  use types_names
  use execution_context_names

  implicit none
# include "debug.i90"
  private
  
  ! What to do in case several pairs of start()-stop() calls are performed?
  character(len=*), parameter :: TIMER_MODE_SUM     = "timer_mode_sum"
  character(len=*), parameter :: TIMER_MODE_MIN     = "timer_mode_min"
  character(len=*), parameter :: TIMER_MODE_LAST    = "timer_mode_last"
  character(len=*), parameter :: DEFAULT_TIMER_MODE = TIMER_MODE_LAST

  type timer_t
     private
     class(execution_context_t), allocatable :: context
     character(:), allocatable :: message  ! Concept being measured (e.g., assembly)
     character(:), allocatable :: mode     ! timer operation mode (see options above)
     real(rp)                  :: t_start  ! last call to start
     real(rp)                  :: t_stop   ! last call to stop
     real(rp)                  :: t_accum  ! sum of all stop-start 
   contains
     procedure, non_overridable :: create => timer_create
     procedure, non_overridable :: free   => timer_free
     procedure, non_overridable :: init   => timer_init
     procedure, non_overridable :: start  => timer_start
     procedure, non_overridable :: stop   => timer_stop
     procedure, non_overridable :: report => timer_report
     procedure, non_overridable :: get_time => timer_get_time
  end type timer_t

  ! Public types
  public :: timer_t 

contains  
    subroutine timer_create ( this, context, message, mode )
      ! Parameters
      class(timer_t)        , intent(inout) :: this
      character(len=*)          , intent(in)    :: message
      class(execution_context_t), intent(in)    :: context
      character(len=*), optional, intent(in)    :: mode
      integer(ip) :: istat

      call this%free()

      if ( present(mode) ) then
        assert ( mode == TIMER_MODE_SUM .or. mode == TIMER_MODE_MIN .or. mode == TIMER_MODE_LAST )
        this%mode = mode
      else
        this%mode = DEFAULT_TIMER_MODE
      end if
      
      allocate(this%context,mold=context,stat=istat);check(istat==0)
      this%context = context
      this%message = message
      call this%init()
    end subroutine timer_create
    
    subroutine timer_free ( this )
      ! Parameters
      class(timer_t), intent(inout) :: this  
      integer(ip) :: istat
      
      if (allocated(this%message) ) then
         deallocate(this%message, stat=istat); check(istat==0);
      end if
      
      if (allocated(this%mode) ) then
         deallocate(this%mode, stat=istat); check(istat==0);
      end if
      
      if (allocated(this%context)) then
        call this%context%free(finalize=.false.)
        deallocate(this%context,stat=istat); check(istat==0);
      end if    
      
      this%t_start   = 0.0_rp
      this%t_stop    = 0.0_rp
      this%t_accum   = 1.79769E+308_rp ! Largest double precision number
    end subroutine timer_free

    subroutine timer_init ( this )
      implicit none 
      ! Parameters
      class(timer_t), intent(inout) :: this
      this%t_start  = 0.0_rp
      this%t_stop   = 0.0_rp
      if ( this%mode == TIMER_MODE_MIN ) then
        this%t_accum   = 1.79769e+308_rp ! Largest double precision number
      else
        this%t_accum   = 0.0_rp
      end if
    end subroutine timer_init

    subroutine timer_start ( this )
      implicit none 
      ! Parameters
      class(timer_t), intent(inout)          :: this
      call this%context%barrier()
      this%t_start  = this%context%time()
    end subroutine timer_start

    subroutine timer_stop ( this )
      implicit none 
      ! Parameters
      class(timer_t), intent(inout)          :: this
      real(rp) :: cur_time

      this%t_stop = this%context%time()
      if ( this%t_stop - this%t_start >= 0.0_rp) then
        cur_time = this%t_stop - this%t_start
      else
        cur_time = 0.0_rp
      end if  
      
      if ( this%mode == TIMER_MODE_MIN ) then
         if ( this%t_accum > cur_time ) this%t_accum = cur_time
      else if ( this%mode == TIMER_MODE_SUM ) then
         this%t_accum = this%t_accum + cur_time
      else if ( this%mode == TIMER_MODE_LAST ) then
         this%t_accum = cur_time
      end if
      
      this%t_start  = 0.0_rp
      this%t_stop   = 0.0_rp
    end subroutine timer_stop
   
    subroutine timer_report ( this, show_header, luout )
      implicit none 
      ! Parameters
      class(timer_t), intent(inout)  :: this
      logical, intent(in), optional      :: show_header 
      integer(ip), intent(in), optional  :: luout
      
      ! Locals
      character(len=*), parameter    :: fmt_header = '(a25,1x,3(2x,a15),3(2x,a15))'
      character(len=*), parameter    :: fmt_data   = '(a25,1x,3(2x,es15.9),3(2x,es15.9))'
      real(rp)                       :: accum_max, accum_min, accum_sum
      logical                        :: show_header_

      accum_max = this%t_accum
      accum_min = this%t_accum
      accum_sum = this%t_accum
      show_header_ = .true.
      if (present(show_header)) show_header_ = show_header 

      if ( show_header_ ) then 
         if (present(luout) ) then
            if ( this%context%am_i_root() ) write(luout,fmt_header) '', 'Min (secs.)', 'Max (secs.)', 'Avg (secs.)'
         else
            if ( this%context%am_i_root() ) write(*,fmt_header) '', 'Min (secs.)', 'Max (secs.)', 'Avg (secs.)'
         end if
      end if         

      call this%context%max_scalar_rp(accum_max)
      call this%context%min_scalar_rp(accum_min)
      call this%context%sum_scalar_rp(accum_sum)
      
      if (present(luout) ) then
         if ( this%context%am_i_root() ) write(luout,fmt_data) adjustl(this%message), accum_min, accum_max, accum_sum/this%context%get_num_tasks()
      else
         if ( this%context%am_i_root() ) write(*,fmt_data) adjustl(this%message), accum_min, accum_max, accum_sum/this%context%get_num_tasks()
      end if

    end subroutine timer_report

    function timer_get_time(this)
      implicit none 
      class(timer_t), intent(in)  :: this
      real(rp) :: timer_get_time
      timer_get_time = this%t_accum
    end function timer_get_time

end module timer_names
