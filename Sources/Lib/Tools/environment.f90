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
module environment_names
  use types_names
  implicit none
  private

  type, abstract :: environment_t
   contains
     procedure (info_interface)                       , deferred  :: info
     procedure (am_i_l1_task_interface)               , deferred  :: am_i_l1_task
     procedure (am_i_l1_root_interface)               , deferred  :: am_i_l1_root
     procedure (get_l1_rank_interface)                , deferred  :: get_l1_rank
     procedure (get_l1_size_interface)                , deferred  :: get_l1_size
     procedure (bcast_interface)                      , deferred  :: l1_lgt1_bcast
     procedure (l1_barrier_interface)                 , deferred  :: l1_barrier
     procedure (l1_sum_scalar_rp_interface), private  , deferred  :: l1_sum_scalar_rp
     procedure (l1_sum_vector_rp_interface), private  , deferred  :: l1_sum_vector_rp
     procedure (l1_max_scalar_rp_interface), private  , deferred  :: l1_max_scalar_rp
     procedure (l1_max_vector_rp_interface), private  , deferred  :: l1_max_vector_rp
     generic  :: l1_sum                                           => l1_sum_scalar_rp, l1_sum_vector_rp
     generic  :: l1_max                                           => l1_max_scalar_rp, l1_max_vector_rp
  end type environment_t

  ! Abstract interfaces
  abstract interface
     subroutine info_interface(this,me,np) 
       import :: environment_t, ip
       implicit none
       class(environment_t),intent(in)  :: this
       integer(ip)         ,intent(out) :: me
       integer(ip)         ,intent(out) :: np
     end subroutine info_interface

     function am_i_l1_task_interface(this) 
       import :: environment_t, ip
       implicit none
       class(environment_t) ,intent(in)  :: this
       logical                           :: am_i_l1_task_interface 
     end function am_i_l1_task_interface

     function am_i_l1_root_interface(this)
       import :: environment_t
       implicit none
       class(environment_t), intent(in) :: this
       logical                              :: am_i_l1_root_interface
     end function am_i_l1_root_interface

     !=============================================================================
     function get_l1_rank_interface ( this )
       import :: environment_t
       implicit none 
       ! Parameters
       class(environment_t), intent(in) :: this
       integer                          :: get_l1_rank_interface
     end function get_l1_rank_interface

     !=============================================================================
     function get_l1_size_interface ( this )
       import :: environment_t
       implicit none 
       ! Parameters
       class(environment_t), intent(in) :: this
       integer                          :: get_l1_size_interface
     end function get_l1_size_interface
     
     !=============================================================================
     subroutine bcast_interface (this, condition)
       import :: environment_t
       implicit none
       class(environment_t) ,intent(in)    :: this
       logical              ,intent(inout) :: condition
     end subroutine bcast_interface

     subroutine l1_barrier_interface (this)
       import :: environment_t
       implicit none
       class(environment_t) ,intent(in)    :: this
     end subroutine l1_barrier_interface

     subroutine l1_sum_scalar_rp_interface (this,alpha)
       import :: environment_t, rp
       implicit none
       class(environment_t) , intent(in)    :: this
       real(rp)             , intent(inout) :: alpha
     end subroutine l1_sum_scalar_rp_interface

     subroutine l1_sum_vector_rp_interface(this,alpha)
       import :: environment_t, rp
       implicit none
       class(environment_t) , intent(in)    :: this
       real(rp)             , intent(inout) :: alpha(:) 
     end subroutine l1_sum_vector_rp_interface

     subroutine l1_max_scalar_rp_interface (this,alpha)
       import :: environment_t, rp
       implicit none
       class(environment_t) , intent(in)    :: this
       real(rp)             , intent(inout) :: alpha
     end subroutine l1_max_scalar_rp_interface

     subroutine l1_max_vector_rp_interface(this,alpha)
       import :: environment_t, rp
       implicit none
       class(environment_t) , intent(in)    :: this
       real(rp)             , intent(inout) :: alpha(:) 
     end subroutine l1_max_vector_rp_interface
  end interface

  public :: environment_t

end module environment_names
