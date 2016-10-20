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
module serial_context_names
  use types_names
  use memor_names
  use par_context_names
  implicit none 
#include "debug.i90"
  private

  
  type, extends(par_context_t) :: serial_context_t
     private 
   contains
     ! These functions should be non_overridable but there is a bug in gfotran
     procedure :: create             => serial_context_create
     procedure :: assign             => serial_context_assign
     procedure :: split_by_condition => serial_context_split_by_condition
     procedure :: split_by_color     => serial_context_split_by_color
     procedure :: free               => serial_context_free
     procedure :: nullify            => serial_context_nullify
     procedure :: am_i_member        => serial_context_am_i_member
     procedure :: barrier            => serial_context_barrier
     procedure :: sum_scalar_rp      => serial_context_sum_scalar_rp
     procedure :: sum_vector_rp      => serial_context_sum_vector_rp
     procedure :: max_scalar_rp      => serial_context_max_scalar_rp
     procedure :: max_vector_rp      => serial_context_max_vector_rp
     procedure :: scatter            => serial_context_scatter_scalar_ip
     procedure :: gather             => serial_context_gather_scalar_ip
     procedure :: bcast              => serial_context_bcast_scalar_ip
     procedure :: bcast_subcontext   => serial_context_bcast_subcontext
     procedure, private :: neighbours_exchange_rp                  =>  serial_context_neighbours_exchange_rp                 
     procedure, private :: neighbours_exchange_ip                  =>  serial_context_neighbours_exchange_ip                 
     procedure, private :: neighbours_exchange_igp                 =>  serial_context_neighbours_exchange_igp                
     procedure, private :: neighbours_exchange_single_ip           =>  serial_context_neighbours_exchange_single_ip          
     procedure, private :: neighbours_exchange_wo_pack_unpack_ieep =>  serial_context_neighbours_exchange_wo_pack_unpack_ieep
     procedure, private :: root_send_master_rcv_ip          => serial_context_root_send_master_rcv_ip
     procedure, private :: root_send_master_rcv_ip_1D_array => serial_context_root_send_master_rcv_ip_1D_array
     procedure, private :: gather_to_master_ip              => serial_context_gather_to_master_ip           
     procedure, private :: gather_to_master_igp             => serial_context_gather_to_master_igp          
     procedure, private :: gather_to_master_ip_1D_array     => serial_context_gather_to_master_ip_1D_array  
     procedure, private :: gather_to_masterv_ip_1D_array    => serial_context_gather_to_masterv_ip_1D_array 
     procedure, private :: gather_to_masterv_igp_1D_array   => serial_context_gather_to_masterv_igp_1D_array
     procedure, private :: gather_to_masterv_rp_1D_array    => serial_context_gather_to_masterv_rp_1D_array 
     procedure, private :: gather_to_masterv_rp_2D_array    => serial_context_gather_to_masterv_rp_2D_array  
     procedure, private :: scatter_from_masterv_rp_1D_array => serial_context_scatter_from_masterv_rp_1D_array

     ! procedure, non_overridable :: create             => serial_context_create
     ! procedure, non_overridable :: assign             => serial_context_assign
     ! procedure, non_overridable :: split_by_condition => serial_context_split_by_condition
     ! procedure, non_overridable :: split_by_color     => serial_context_split_by_color
     ! procedure, non_overridable :: free               => serial_context_free
     ! procedure, non_overridable :: nullify            => serial_context_nullify
     ! procedure, non_overridable :: am_i_member        => serial_context_am_i_member
     ! procedure, non_overridable :: barrier            => serial_context_barrier
     ! procedure, non_overridable :: sum_scalar_rp      => serial_context_sum_scalar_rp
     ! procedure, non_overridable :: sum_vector_rp      => serial_context_sum_vector_rp
     ! procedure, non_overridable :: max_scalar_rp      => serial_context_max_scalar_rp
     ! procedure, non_overridable :: max_vector_rp      => serial_context_max_vector_rp
     ! procedure, non_overridable :: scatter            => serial_context_scatter_scalar_ip
     ! procedure, non_overridable :: gather             => serial_context_gather_scalar_ip
     ! procedure, non_overridable :: bcast              => serial_context_bcast_scalar_ip
     ! procedure, non_overridable :: bcast_subcontext   => serial_context_bcast_subcontext
     ! procedure, non_overridable, private :: neighbours_exchange_rp                  =>  serial_context_neighbours_exchange_rp                 
     ! procedure, non_overridable, private :: neighbours_exchange_ip                  =>  serial_context_neighbours_exchange_ip                 
     ! procedure, non_overridable, private :: neighbours_exchange_igp                 =>  serial_context_neighbours_exchange_igp                
     ! procedure, non_overridable, private :: neighbours_exchange_single_ip           =>  serial_context_neighbours_exchange_single_ip          
     ! procedure, non_overridable, private :: neighbours_exchange_wo_pack_unpack_ieep =>  serial_context_neighbours_exchange_wo_pack_unpack_ieep
     ! procedure, non_overridable, private :: root_send_master_rcv_ip          => serial_context_root_send_master_rcv_ip
     ! procedure, non_overridable, private :: root_send_master_rcv_ip_1D_array => serial_context_root_send_master_rcv_ip_1D_array
     ! procedure, non_overridable, private :: gather_to_master_ip              => serial_context_gather_to_master_ip           
     ! procedure, non_overridable, private :: gather_to_master_igp             => serial_context_gather_to_master_igp          
     ! procedure, non_overridable, private :: gather_to_master_ip_1D_array     => serial_context_gather_to_master_ip_1D_array  
     ! procedure, non_overridable, private :: gather_to_masterv_ip_1D_array    => serial_context_gather_to_masterv_ip_1D_array 
     ! procedure, non_overridable, private :: gather_to_masterv_igp_1D_array   => serial_context_gather_to_masterv_igp_1D_array
     ! procedure, non_overridable, private :: gather_to_masterv_rp_1D_array    => serial_context_gather_to_masterv_rp_1D_array 
     ! procedure, non_overridable, private :: gather_to_masterv_rp_2D_array    => serial_context_gather_to_masterv_rp_2D_array  
     ! procedure, non_overridable, private :: scatter_from_masterv_rp_1D_array => serial_context_scatter_from_masterv_rp_1D_array
  end type serial_context_t

  ! Types
  public :: serial_context_t

contains

  !=============================================================================
  subroutine serial_context_assign(this, that)
    implicit none 
    class(serial_context_t), intent(inout) :: this
    class(par_context_t)   , intent(in)    :: that
    select type(that)
    type is(serial_context_t)
       assert(that%get_current_task()==0)
       assert(that%get_num_tasks()==1)
       call this%set_current_task(0)
       call this%set_num_tasks(1)
    class default
       check(.false.)
    end select
  end subroutine serial_context_assign

  !=============================================================================
  subroutine serial_context_create ( this )
    implicit none 
    class(serial_context_t), intent(inout) :: this
    call this%set_current_task(0)
    call this%set_num_tasks(1)
  end subroutine serial_context_create

  !=============================================================================
  subroutine serial_context_split_by_color ( this, color, new_subcontext )
    implicit none 
    class(serial_context_t), intent(in)    :: this
    integer                , intent(in)    :: color
    class(par_context_t), allocatable , intent(inout) :: new_subcontext
  end subroutine serial_context_split_by_color

  !=============================================================================
  subroutine serial_context_split_by_condition ( this, in_subcontext1, subcontext1, subcontext2 )
    implicit none 
    class(serial_context_t), intent(in)    :: this
    logical                , intent(in)    :: in_subcontext1
    class(par_context_t), allocatable , intent(inout) :: subcontext1
    class(par_context_t), allocatable , intent(inout) :: subcontext2
  end subroutine serial_context_split_by_condition

  !=============================================================================
  subroutine serial_context_free ( p_context, finalize  )
    implicit none 
    class(serial_context_t), intent(inout) :: p_context
    logical                , intent(in)    :: finalize
    call p_context%set_current_task(-1)
    call p_context%set_num_tasks(-1)
  end subroutine serial_context_free

  !=============================================================================
  subroutine serial_context_nullify ( this )
    implicit none 
    class(serial_context_t), intent(inout) :: this
  end subroutine serial_context_nullify

  !=============================================================================
  ! pure function serial_context_get_current_task (this)
  !   implicit none
  !   class(serial_context_t), intent(in) :: this
  !   integer :: serial_context_get_current_task
  !   serial_context_get_current_task = 0
  ! end function serial_context_get_current_task

  !=============================================================================
  ! pure function serial_context_get_num_tasks (this)
  !   implicit none
  !   class(serial_context_t), intent(in) :: this
  !   integer :: serial_context_get_num_tasks
  !   serial_context_get_num_tasks = 1
  ! end function serial_context_get_num_tasks

  !=============================================================================
  pure function serial_context_am_i_member(this)
    implicit none
    class(serial_context_t), intent(in) :: this
    logical                          :: serial_context_am_i_member
    serial_context_am_i_member = .true.
  end function serial_context_am_i_member

  !=============================================================================
  subroutine serial_context_barrier(this)
    implicit none 
    class(serial_context_t), intent(in) :: this
  end subroutine serial_context_barrier

  !=============================================================================
  subroutine serial_context_sum_scalar_rp (this,alpha)
    implicit none
    class(serial_context_t) , intent(in)    :: this
    real(rp)             , intent(inout) :: alpha
  end subroutine serial_context_sum_scalar_rp

  !=============================================================================
  subroutine serial_context_sum_vector_rp(this,alpha)
    implicit none
    class(serial_context_t) , intent(in)    :: this
    real(rp)             , intent(inout) :: alpha(:) 
  end subroutine serial_context_sum_vector_rp

  !=============================================================================
  subroutine serial_context_max_scalar_rp (this,alpha)
    implicit none
    class(serial_context_t) , intent(in)    :: this
    real(rp)             , intent(inout) :: alpha
  end subroutine serial_context_max_scalar_rp

  !=============================================================================
  subroutine serial_context_max_vector_rp(this,alpha)
    implicit none
    class(serial_context_t) , intent(in)    :: this
    real(rp)             , intent(inout) :: alpha(:) 
  end subroutine serial_context_max_vector_rp

  !=============================================================================
  subroutine serial_context_bcast_subcontext(this,subcontxt,condition)
    implicit none
    class(serial_context_t) , intent(in)    :: this
    class(par_context_t) , intent(in)    :: subcontxt
    logical              , intent(inout) :: condition
  end subroutine serial_context_bcast_subcontext

  !=============================================================================
  ! When packing   (gathering) ,    buffer <- alpha * x
  ! When unpacking (scattering),    x <- beta*x + buffer
  subroutine serial_context_neighbours_exchange_rp ( this, & 
       &                                          num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       &                                          num_snd, list_snd, snd_ptrs, pack_idx,   &
       &                                          alpha, beta, x)
    implicit none
    class(serial_context_t), intent(in) :: this
    integer(ip)             , intent(in) :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(in) :: unpack_idx (rcv_ptrs(num_rcv+1)-1)
    integer(ip)             , intent(in) :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in) :: pack_idx (snd_ptrs(num_snd+1)-1)
    real(rp), intent(in)    :: alpha, beta
    real(rp), intent(inout) :: x(:)

  end subroutine serial_context_neighbours_exchange_rp

  !=============================================================================
  subroutine serial_context_neighbours_exchange_ip ( this, & 
       &                                          num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       &                                          num_snd, list_snd, snd_ptrs, pack_idx,   &
       &                                          x,chunk_size)
    implicit none
    class(serial_context_t), intent(in)    :: this
    integer(ip)             , intent(in)    :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(in)    :: unpack_idx (rcv_ptrs(num_rcv+1)-1)
    integer(ip)             , intent(in)    :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in)    :: pack_idx (snd_ptrs(num_snd+1)-1)
    integer(ip)             , intent(inout) :: x(:)
    integer(ip)   , optional, intent(in)    :: chunk_size
  end subroutine serial_context_neighbours_exchange_ip

  !=============================================================================
  subroutine serial_context_neighbours_exchange_igp ( this, & 
       &                                              num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       &                                              num_snd, list_snd, snd_ptrs, pack_idx,   &
       &                                              x, chunk_size)
    implicit none
    class(serial_context_t), intent(in)    :: this
    integer(ip)             , intent(in)    :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(in)    :: unpack_idx (rcv_ptrs(num_rcv+1)-1)
    integer(ip)             , intent(in)    :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in)    :: pack_idx (snd_ptrs(num_snd+1)-1)
    integer(igp)            , intent(inout) :: x(:)
    integer(ip)   , optional, intent(in)    :: chunk_size

  end subroutine serial_context_neighbours_exchange_igp

  !=============================================================================
  subroutine serial_context_neighbours_exchange_single_ip ( this, & 
       &                                                    num_neighbours, &
       &                                                    list_neighbours, &
       &                                                    input_data,&
       &                                                    output_data)
    implicit none
    class(serial_context_t), intent(in) :: this
    integer                 , intent(in)    :: num_neighbours
    integer(ip)             , intent(in)    :: list_neighbours (num_neighbours)
    integer(ip)             , intent(in)    :: input_data
    integer(ip)             , intent(inout) :: output_data(num_neighbours)
  end subroutine serial_context_neighbours_exchange_single_ip

  !=============================================================================
  subroutine serial_context_neighbours_exchange_wo_pack_unpack_ieep ( this, &
       &                                                              number_neighbours, &
       &                                                              neighbour_ids, &
       &                                                              snd_ptrs, &
       &                                                              snd_buf, & 
       &                                                              rcv_ptrs, &
       &                                                              rcv_buf )
    implicit none
    class(serial_context_t)  , intent(in)    :: this 
    integer(ip)           , intent(in)    :: number_neighbours
    integer(ip)           , intent(in)    :: neighbour_ids(number_neighbours)
    integer(ip)           , intent(in)    :: snd_ptrs(number_neighbours+1)
    integer(ieep)         , intent(in)    :: snd_buf(snd_ptrs(number_neighbours+1)-1)   
    integer(ip)           , intent(in)    :: rcv_ptrs(number_neighbours+1)
    integer(ieep)         , intent(out)   :: rcv_buf(rcv_ptrs(number_neighbours+1)-1)
    rcv_buf = snd_buf ! needed to satisfy the intent
  end subroutine serial_context_neighbours_exchange_wo_pack_unpack_ieep

  !=============================================================================
  !=============================================================================
  subroutine serial_context_gather_scalar_ip ( this, input_data, output_data )
    implicit none
    class(serial_context_t), intent(in)   :: this
    integer(ip)         , intent(in)   :: input_data
    integer(ip)         , intent(out)  :: output_data(:)
    check(.false.)       ! This routine should be never called
  end subroutine serial_context_gather_scalar_ip

  !=============================================================================
  subroutine serial_context_scatter_scalar_ip ( this, input_data, output_data )
    implicit none
    class(serial_context_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data(:)
    integer(ip)             , intent(out)  :: output_data
    check(.false.)       ! This routine should be never called
  end subroutine serial_context_scatter_scalar_ip

  !=============================================================================
  subroutine serial_context_bcast_scalar_ip ( this, data )
    implicit none
    class(serial_context_t), intent(in)    :: this
    integer(ip)             , intent(inout) :: data
    check(.false.)       ! This routine should be never called
  end subroutine serial_context_bcast_scalar_ip

  !=============================================================================
  subroutine serial_context_root_send_master_rcv_ip ( this, input_data, output_data )
    implicit none
    class(serial_context_t), intent(in)      :: this
    integer(ip)         , intent(in)      :: input_data
    integer(ip)         , intent(inout)   :: output_data
    check(.false.)       ! This routine should be never called
  end subroutine serial_context_root_send_master_rcv_ip

  !=============================================================================
  subroutine serial_context_root_send_master_rcv_ip_1D_array ( this, input_data, output_data )
    implicit none
    class(serial_context_t), intent(in)      :: this
    integer(ip)         , intent(in)      :: input_data(:)
    integer(ip)         , intent(inout)   :: output_data(:)
    check(.false.)       ! This routine should be never called
  end subroutine serial_context_root_send_master_rcv_ip_1D_array

  !=============================================================================
  !=============================================================================
  subroutine serial_context_gather_to_master_ip ( this, input_data, output_data )
    implicit none
    class(serial_context_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data
    integer(ip)             , intent(out)  :: output_data(:)
    check(.false.)       ! This routine should be never called
    output_data = 0      ! needed to satisfy the intent
  end subroutine serial_context_gather_to_master_ip

  !=============================================================================
  subroutine serial_context_gather_to_master_igp ( this, input_data, output_data )
    implicit none
    class(serial_context_t), intent(in)   :: this
    integer(igp)            , intent(in)   :: input_data
    integer(igp)            , intent(out)  :: output_data(:)
    check(.false.)       ! This routine should be never called
    output_data = 0      ! needed to satisfy the intent
  end subroutine serial_context_gather_to_master_igp

  !=============================================================================
  subroutine serial_context_gather_to_master_ip_1D_array ( this, input_data_size, input_data, output_data )
    implicit none
    class(serial_context_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data_size
    integer(ip)             , intent(in)   :: input_data(input_data_size)
    integer(ip)             , intent(out)  :: output_data(:)
    check(.false.)       ! This routine should be never called
    output_data = 0      ! needed to satisfy the intent
  end subroutine serial_context_gather_to_master_ip_1D_array

  !=============================================================================
  subroutine serial_context_gather_to_masterv_ip_1D_array ( this, input_data_size, input_data, recv_counts, displs, output_data )
    implicit none
    class(serial_context_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data_size
    integer(ip)             , intent(in)   :: input_data(input_data_size)
    integer(ip)             , intent(in)   :: recv_counts(:)
    integer(ip)             , intent(in)   :: displs(:)
    integer(ip)             , intent(out)  :: output_data(:)
    check(.false.)       ! This routine should be never called
    output_data = 0      ! needed to satisfy the intent
  end subroutine serial_context_gather_to_masterv_ip_1D_array

  !=============================================================================
  subroutine serial_context_gather_to_masterv_igp_1D_array ( this, input_data_size, input_data, recv_counts, displs, output_data )
    implicit none
    class(serial_context_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data_size
    integer(igp)            , intent(in)   :: input_data(input_data_size)
    integer(ip)             , intent(in)   :: recv_counts(:)
    integer(ip)             , intent(in)   :: displs(:)
    integer(igp)            , intent(out)  :: output_data(:)
    check(.false.)       ! This routine should be never called
    output_data = 0.0_rp ! needed to satisfy the intent
  end subroutine serial_context_gather_to_masterv_igp_1D_array

  !=============================================================================
  subroutine serial_context_gather_to_masterv_rp_1D_array ( this, input_data_size, input_data, recv_counts, displs, output_data )
    implicit none
    class(serial_context_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data_size
    real(rp)                , intent(in)   :: input_data(input_data_size)
    integer(ip)             , intent(in)   :: recv_counts(:)
    integer(ip)             , intent(in)   :: displs(:)
    real(rp)                , intent(out)  :: output_data(:)
    check(.false.)       ! This routine should be never called
    output_data = 0.0_rp ! needed to satisfy the intent
  end subroutine serial_context_gather_to_masterv_rp_1D_array

  !=============================================================================
  subroutine serial_context_gather_to_masterv_rp_2D_array ( this, input_data, recv_counts, displs, output_data )
    implicit none
    class(serial_context_t), intent(in)   :: this
    real(rp)                , intent(in)   :: input_data(:,:)
    integer(ip)             , intent(in)   :: recv_counts(:)
    integer(ip)             , intent(in)   :: displs(:)
    real(rp)                , intent(out)  :: output_data(:)
    check(.false.)       ! This routine should be never called
    output_data = 0.0_rp ! needed to satisfy the intent
  end subroutine serial_context_gather_to_masterv_rp_2D_array

  !=============================================================================
  subroutine serial_context_scatter_from_masterv_rp_1D_array ( this, input_data, send_counts, displs, output_data_size, output_data )
    implicit none
    class(serial_context_t), intent(in)   :: this
    real(rp)               , intent(in)   :: input_data(:)
    integer(ip)            , intent(in)   :: send_counts(:)
    integer(ip)            , intent(in)   :: displs(:)
    integer(ip)            , intent(in)   :: output_data_size
    real(rp)               , intent(out)  :: output_data(output_data_size)
    check(.false.)       ! This routine should be never called
    output_data = 0.0_rp ! needed to satisfy the intent
  end subroutine serial_context_scatter_from_masterv_rp_1D_array

end module serial_context_names
