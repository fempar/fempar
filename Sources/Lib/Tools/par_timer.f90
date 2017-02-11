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

  implicit none
# include "debug.i90"
  private
  
  ! What to do in case several pairs of start()-stop() calls are performed?
  character(len=*), parameter :: PAR_TIMER_MODE_SUM     = "par_timer_mode_sum"
  character(len=*), parameter :: PAR_TIMER_MODE_MIN     = "par_timer_mode_min"
  character(len=*), parameter :: PAR_TIMER_MODE_LAST    = "par_timer_mode_last"
  character(len=*), parameter :: DEFAULT_PAR_TIMER_MODE = PAR_TIMER_MODE_LAST

  type par_timer_t
     private
     integer                   :: ictxt    ! Parallel context
     integer                   :: root     ! PID responsible of gathering/reporting timings
     character(:), allocatable :: message  ! Concept being measured (e.g., assembly)
     character(:), allocatable :: mode     ! par_timer operation mode (see options above)
     real(8)                   :: t_start  ! last call to start
     real(8)                   :: t_stop   ! last call to stop
     real(8)                   :: t_accum  ! sum of all stop-start 
   contains
     procedure, non_overridable :: create => par_timer_create
     procedure, non_overridable :: free   => par_timer_free
     procedure, non_overridable :: init   => par_timer_init
     procedure, non_overridable :: start  => par_timer_start
     procedure, non_overridable :: stop   => par_timer_stop
     procedure, non_overridable :: report => par_timer_report
  end type par_timer_t

  ! Public types
  public :: par_timer_t 

contains  
    subroutine par_timer_create ( this, message, ictxt, root, mode )
#ifdef MPI_MOD
      use mpi
#endif
      implicit none 
#ifdef MPI_H
      include 'mpif.h'
#endif
      ! Parameters
      class(par_timer_t), intent(inout)        :: this
      character(len=*) , intent(in)            :: message
      integer          , intent(in), optional  :: ictxt
      integer          , intent(in), optional  :: root
      character(len=*) , intent(in), optional  :: mode
      
      ! Locals
      integer                               :: ictxt_, root_

      call this%free()
      
      ictxt_ = mpi_comm_world
      root_  = 0 
      if ( present(ictxt) ) ictxt_ = ictxt
      if ( present(root)  ) root_  = root 
      if ( present(mode) ) then
        assert ( mode == PAR_TIMER_MODE_SUM .or. mode == PAR_TIMER_MODE_MIN .or. mode == PAR_TIMER_MODE_LAST )
        this%mode = mode
      else
        this%mode = DEFAULT_PAR_TIMER_MODE
      end if
      
      this%ictxt   = ictxt_
      this%root    = root_ 
      this%message = message
      this%t_start   = 0.0 
      this%t_stop    = 0.0
      if ( this%mode == PAR_TIMER_MODE_MIN ) then
        this%t_accum   = 1.79769E+308 ! Largest double precision number
      else
        this%t_accum   = 0.0
      end if
    end subroutine par_timer_create
    
    subroutine par_timer_free ( this )
#ifdef MPI_MOD
      use mpi
#endif
      implicit none 
#ifdef MPI_H
      include 'mpif.h'
#endif
      ! Parameters
      class(par_timer_t), intent(inout) :: this  
      integer(ip) :: istat
      
      if (allocated(this%message) ) then
         deallocate(this%message, stat=istat); check(istat==0);
      end if
      
      if (allocated(this%mode) ) then
         deallocate(this%mode, stat=istat); check(istat==0);
      end if
      
      this%ictxt     = mpi_comm_null
      this%root      = -1
      this%t_start   = 0.0 
      this%t_stop    = 0.0
      this%t_accum   = 1.79769E+308 ! Largest double precision number
    end subroutine par_timer_free

    subroutine par_timer_init ( this )
      implicit none 
      ! Parameters
      class(par_timer_t), intent(inout) :: this
      this%t_start  = 0.0 
      this%t_stop   = 0.0
      if ( this%mode == PAR_TIMER_MODE_MIN ) then
        this%t_accum   = 1.79769E+308 ! Largest double precision number
      else
        this%t_accum   = 0.0
      end if
    end subroutine par_timer_init

    subroutine par_timer_start ( this )
      implicit none 
      ! Parameters
      class(par_timer_t), intent(inout)          :: this
      call psb_barrier (this%ictxt)
      this%t_start  = psb_wtime()
    end subroutine par_timer_start

    subroutine par_timer_stop ( this )
      implicit none 
      ! Parameters
      class(par_timer_t), intent(inout)          :: this
      real(8) :: cur_time
      
      this%t_stop = psb_wtime()
      
      if ( this%t_stop - this%t_start >= 0.0) then
        cur_time = this%t_stop - this%t_start
      else
        cur_time = 0.0
      end if  
      
      if ( this%mode == PAR_TIMER_MODE_MIN ) then
         if ( this%t_accum > cur_time ) this%t_accum = cur_time
      else if ( this%mode == PAR_TIMER_MODE_SUM ) then
         this%t_accum = this%t_accum + cur_time
      else if ( this%mode == PAR_TIMER_MODE_LAST ) then
         this%t_accum = cur_time
      end if
      
      this%t_start  = 0.0 
      this%t_stop   = 0.0
    end subroutine par_timer_stop
   
    subroutine par_timer_report ( this, show_header, luout )
      implicit none 
      ! Parameters
      class(par_timer_t), intent(inout)     :: this
      logical, intent(in), optional      :: show_header 
      integer(ip), intent(in), optional  :: luout
      
      ! Locals
      character(len=*), parameter    :: fmt_header = '(a25,1x,3(2x,a15),3(2x,a15))'
      character(len=*), parameter    :: fmt_data   = '(a25,1x,3(2x,es15.9),3(2x,es15.9))'
      real(8)                        :: accum_max, accum_min, accum_sum
      integer                        :: my_id, num_procs
      logical                        :: show_header_

      call psb_info ( this%ictxt, my_id, num_procs )

      accum_max = this%t_accum
      accum_min = this%t_accum
      accum_sum = this%t_accum
      show_header_ = .true.
      if (present(show_header)) show_header_ = show_header 

      if ( show_header_ ) then 
         if (present(luout) ) then
            if ( my_id == this%root ) write(luout,fmt_header) '', 'Min (secs.)', 'Max (secs.)', 'Avg (secs.)'
         else
            if ( my_id == this%root ) write(*,fmt_header) '', 'Min (secs.)', 'Max (secs.)', 'Avg (secs.)'
         end if
      end if         
  
      call psb_max ( this%ictxt, accum_max, this%root )
      call psb_min ( this%ictxt, accum_min, this%root )
      call psb_sum ( this%ictxt, accum_sum, this%root )
      
      if (present(luout) ) then
         if ( my_id == this%root ) write(luout,fmt_data) adjustl(this%message), accum_min, accum_max, accum_sum/num_procs
      else
         if ( my_id == this%root ) write(*,fmt_data) adjustl(this%message), accum_min, accum_max, accum_sum/num_procs
      end if

    end subroutine par_timer_report

end module par_timer_names
