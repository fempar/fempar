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
  use execution_context_names
  implicit none 
#include "debug.i90"
  private

  
  type, extends(execution_context_t) :: serial_context_t
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
     procedure :: am_i_root          => serial_context_am_i_root
     procedure :: barrier            => serial_context_barrier
     procedure :: time               => serial_context_time
     procedure :: sum_scalar_rp      => serial_context_sum_scalar_rp
     procedure :: sum_vector_rp      => serial_context_sum_vector_rp
     procedure :: max_scalar_rp      => serial_context_max_scalar_rp
     procedure :: max_vector_rp      => serial_context_max_vector_rp
     procedure :: min_scalar_rp      => serial_context_min_scalar_rp
     procedure :: max_scalar_ip      => serial_context_max_scalar_ip
     procedure :: sum_scalar_igp     => serial_context_sum_scalar_igp
     procedure :: sum_vector_igp     => serial_context_sum_vector_igp
     procedure :: scatter_ip         => serial_context_scatter_scalar_ip
     procedure :: gather_ip          => serial_context_gather_scalar_ip
     procedure :: bcast_ip           => serial_context_bcast_scalar_ip
     procedure :: bcast_ip_1D_array  => serial_context_bcast_scalar_ip_1D_array
     procedure :: scatter_igp        => serial_context_scatter_scalar_igp
     procedure :: gather_igp         => serial_context_gather_scalar_igp
     procedure :: bcast_igp          => serial_context_bcast_scalar_igp
     procedure :: bcast_subcontext   => serial_context_bcast_subcontext
     procedure :: neighbours_exchange_rp                  =>  serial_context_neighbours_exchange_rp        
     procedure :: neighbours_exchange_wo_alpha_beta_rp    =>  serial_context_neighbours_exchange_wo_alpha_beta_rp
     procedure :: neighbours_exchange_wo_alpha_beta_rp_v  =>  serial_context_neighbours_exchange_wo_alpha_beta_rp_v
     procedure :: neighbours_exchange_ip                  =>  serial_context_neighbours_exchange_ip                 
     procedure :: neighbours_exchange_igp                 =>  serial_context_neighbours_exchange_igp                
     procedure :: neighbours_exchange_single_ip           =>  serial_context_neighbours_exchange_single_ip          
     procedure :: neighbours_exchange_wo_pack_unpack_ieep =>  serial_context_neighbours_exchange_wo_pack_unpack_ieep
     procedure :: neighbours_exchange_wo_unpack_ip        =>  serial_context_neighbours_exchange_wo_unpack_ip
     procedure :: neighbours_exchange_variable_igp        =>  serial_context_neighbours_exchange_variable_igp       
     procedure :: neighbours_exchange_variable_ip         =>  serial_context_neighbours_exchange_variable_ip       
     procedure :: send_ip           => serial_context_send_ip
     procedure :: send_igp          => serial_context_send_igp     
     procedure :: send_rp           => serial_context_send_rp
     procedure :: send_ip_1D_array  => serial_context_send_ip_1D_array
     procedure :: send_igp_1D_array => serial_context_send_igp_1D_array     
     procedure :: send_rp_1D_array  => serial_context_send_rp_1D_array
     procedure :: rcv_ip            => serial_context_rcv_ip
     procedure :: rcv_igp           => serial_context_rcv_igp       
     procedure :: rcv_rp            => serial_context_rcv_rp
     procedure :: rcv_ip_1D_array   => serial_context_rcv_ip_1D_array
     procedure :: rcv_igp_1D_array  => serial_context_rcv_igp_1D_array     
     procedure :: rcv_rp_1D_array   => serial_context_rcv_rp_1D_array
     procedure :: root_send_master_rcv_ip          => serial_context_root_send_master_rcv_ip
     procedure :: root_send_master_rcv_ip_1D_array => serial_context_root_send_master_rcv_ip_1D_array
     procedure :: root_send_master_rcv_rp          => serial_context_root_send_master_rcv_rp
     procedure :: root_send_master_rcv_rp_1D_array => serial_context_root_send_master_rcv_rp_1D_array
     procedure :: root_send_master_rcv_logical     => serial_context_root_send_master_rcv_logical
     procedure :: gather_to_master_ip              => serial_context_gather_to_master_ip           
     procedure :: gather_to_master_igp             => serial_context_gather_to_master_igp          
     procedure :: gather_to_master_ip_1D_array     => serial_context_gather_to_master_ip_1D_array  
     procedure :: gather_to_masterv_ip_1D_array    => serial_context_gather_to_masterv_ip_1D_array 
     procedure :: gather_to_masterv_igp_1D_array   => serial_context_gather_to_masterv_igp_1D_array
     procedure :: gather_to_masterv_rp_1D_array    => serial_context_gather_to_masterv_rp_1D_array 
     procedure :: gather_to_masterv_rp_2D_array    => serial_context_gather_to_masterv_rp_2D_array  
     procedure :: scatter_from_master_ip           => serial_context_scatter_from_master_ip
     procedure :: scatter_from_masterv_ip_1D_array => serial_context_scatter_from_masterv_ip_1D_array
     procedure :: scatter_from_masterv_rp_1D_array => serial_context_scatter_from_masterv_rp_1D_array
  end type serial_context_t

  ! Types
  public :: serial_context_t

contains

  !=============================================================================
  subroutine serial_context_assign(this, that)
    implicit none 
    class(serial_context_t), intent(inout) :: this
    class(execution_context_t)   , intent(in)    :: that
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
    class(execution_context_t), allocatable , intent(inout) :: new_subcontext
  end subroutine serial_context_split_by_color

  !=============================================================================
  subroutine serial_context_split_by_condition ( this, in_subcontext1, subcontext1, subcontext2 )
    implicit none 
    class(serial_context_t), intent(in)    :: this
    logical                , intent(in)    :: in_subcontext1
    class(execution_context_t), allocatable , intent(inout) :: subcontext1
    class(execution_context_t), allocatable , intent(inout) :: subcontext2
  end subroutine serial_context_split_by_condition

  !=============================================================================
  subroutine serial_context_free ( this, finalize  )
    implicit none 
    class(serial_context_t), intent(inout) :: this
    logical                , intent(in)    :: finalize
    call this%set_current_task(-1)
    call this%set_num_tasks(-1)
  end subroutine serial_context_free

  !=============================================================================
  subroutine serial_context_nullify ( this )
    implicit none 
    class(serial_context_t), intent(inout) :: this
  end subroutine serial_context_nullify

  !=============================================================================
  pure function serial_context_am_i_member(this)
    implicit none
    class(serial_context_t), intent(in) :: this
    logical                          :: serial_context_am_i_member
    serial_context_am_i_member = .true.
  end function serial_context_am_i_member

  !=============================================================================
  pure function serial_context_am_i_root(this)
    implicit none
    class(serial_context_t), intent(in) :: this
    logical                          :: serial_context_am_i_root
    serial_context_am_i_root = .true.
  end function serial_context_am_i_root

  !=============================================================================
  subroutine serial_context_barrier(this)
    implicit none 
    class(serial_context_t), intent(in) :: this
  end subroutine serial_context_barrier

  !=============================================================================
  function serial_context_time(this)
    implicit none 
    class(serial_context_t), intent(in) :: this
    real(rp) :: serial_context_time
    integer :: clock_reading, clock_rate, clock_max
    call system_clock ( clock_reading, clock_rate, clock_max ) 
    serial_context_time =  real ( clock_reading, kind = rp ) / real ( clock_rate, kind = rp)     
  end function serial_context_time

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
  subroutine serial_context_min_scalar_rp (this,alpha)
    implicit none
    class(serial_context_t) , intent(in)    :: this
    real(rp)             , intent(inout) :: alpha
  end subroutine serial_context_min_scalar_rp
  
  !=============================================================================
  subroutine serial_context_max_scalar_ip (this,n)
    implicit none
    class(serial_context_t) , intent(in)    :: this
    integer(ip)             , intent(inout) :: n
  end subroutine serial_context_max_scalar_ip
  
  !=============================================================================
  subroutine serial_context_sum_scalar_igp (this,n)
    implicit none
    class(serial_context_t) , intent(in)    :: this
    integer(igp)            , intent(inout) :: n
  end subroutine serial_context_sum_scalar_igp
  
  !=============================================================================
  subroutine serial_context_sum_vector_igp (this,n)
    implicit none
    class(serial_context_t) , intent(in)    :: this
    integer(igp)            , intent(inout) :: n(:)
  end subroutine serial_context_sum_vector_igp

  !=============================================================================
  subroutine serial_context_bcast_subcontext(this,subcontxt1,subcontxt2,condition)
    implicit none
    class(serial_context_t) , intent(in)    :: this
    class(execution_context_t) , intent(in)    :: subcontxt1
    class(execution_context_t) , intent(in)    :: subcontxt2
    logical              , intent(inout) :: condition
  end subroutine serial_context_bcast_subcontext

  !=============================================================================
  ! When packing   (gathering) ,    buffer <- alpha * x
  ! When unpacking (scattering),    x <- beta*x + buffer
  subroutine serial_context_neighbours_exchange_rp ( this, & 
       &                                          num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       &                                          num_snd, list_snd, snd_ptrs, pack_idx,   &
       &                                          alpha, beta, x, y)
    implicit none
    class(serial_context_t), intent(inout) :: this
    integer(ip)             , intent(in) :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(in) :: unpack_idx (rcv_ptrs(num_rcv+1)-1)
    integer(ip)             , intent(in) :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in) :: pack_idx (snd_ptrs(num_snd+1)-1)
    real(rp), intent(in)    :: alpha, beta
    real(rp), intent(in)    :: x(:)
    real(rp), intent(inout) :: y(:)
  end subroutine serial_context_neighbours_exchange_rp
  
  !=============================================================================
  subroutine serial_context_neighbours_exchange_wo_alpha_beta_rp ( this, & 
       &                                                        num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       &                                                        num_snd, list_snd, snd_ptrs, pack_idx,   &
       &                                                        x, y, chunk_size)
    implicit none
    class(serial_context_t) , intent(in)    :: this
    ! Control info to receive
    integer(ip)             , intent(in)    :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(in)    :: unpack_idx (rcv_ptrs(num_rcv+1)-1)
    ! Control info to send
    integer(ip)             , intent(in)    :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in)    :: pack_idx (snd_ptrs(num_snd+1)-1)
    ! Raw data to be exchanged
    real(rp)                , intent(in)    :: x(:)
    real(rp)                , intent(inout) :: y(:)
    integer(ip)   , optional, intent(in)    :: chunk_size  
    check(.false.)
  end subroutine serial_context_neighbours_exchange_wo_alpha_beta_rp
  
  !=============================================================================
  subroutine serial_context_neighbours_exchange_wo_alpha_beta_rp_v ( this, & 
       &                                                        num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       &                                                        num_snd, list_snd, snd_ptrs, pack_idx,   &
       &                                                        x, y, ptr_chunk_size_snd, ptr_chunk_size_rcv )
    implicit none
    class(serial_context_t) , intent(in)    :: this
    ! Control info to receive
    integer(ip)             , intent(in)    :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(in)    :: unpack_idx (rcv_ptrs(num_rcv+1)-1)
    ! Control info to send
    integer(ip)             , intent(in)    :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in)    :: pack_idx (snd_ptrs(num_snd+1)-1)
    ! Raw data to be exchanged
    real(rp)                , intent(in)    :: x(:)
    real(rp)                , intent(inout) :: y(:)
    integer(ip)             , intent(in)    :: ptr_chunk_size_snd(:)
    integer(ip)             , intent(in)    :: ptr_chunk_size_rcv(:)
    mcheck(.false.,'serial_context_neighbours_exchange_wo_alpha_beta_rp_v is not implemented')
  end subroutine serial_context_neighbours_exchange_wo_alpha_beta_rp_v
  
  !=============================================================================
  subroutine serial_context_neighbours_exchange_ip ( this, & 
       &                                          num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       &                                          num_snd, list_snd, snd_ptrs, pack_idx,   &
       &                                          x,y,chunk_size)
    implicit none
    class(serial_context_t), intent(in)    :: this
    integer(ip)             , intent(in)    :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(in)    :: unpack_idx (rcv_ptrs(num_rcv+1)-1)
    integer(ip)             , intent(in)    :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in)    :: pack_idx (snd_ptrs(num_snd+1)-1)
    integer(ip)             , intent(in)    :: x(:)
    integer(ip)             , intent(inout) :: y(:)
    integer(ip)   , optional, intent(in)    :: chunk_size
  end subroutine serial_context_neighbours_exchange_ip

  !=============================================================================
  subroutine serial_context_neighbours_exchange_igp ( this, & 
       &                                              num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       &                                              num_snd, list_snd, snd_ptrs, pack_idx,   &
       &                                              x, y, chunk_size,mask)
    implicit none
    class(serial_context_t), intent(in)    :: this
    integer(ip)             , intent(in)    :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(in)    :: unpack_idx (rcv_ptrs(num_rcv+1)-1)
    integer(ip)             , intent(in)    :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in)    :: pack_idx (snd_ptrs(num_snd+1)-1)
    integer(igp)            , intent(in)    :: x(:)
    integer(igp)            , intent(inout) :: y(:)
    integer(ip)   , optional, intent(in)    :: chunk_size
    integer(igp)  , optional, intent(in)    :: mask
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
       &                                                              num_neighbours, &
       &                                                              neighbour_ids, &
       &                                                              snd_ptrs, &
       &                                                              snd_buf, & 
       &                                                              rcv_ptrs, &
       &                                                              rcv_buf )
    implicit none
    class(serial_context_t)  , intent(in)    :: this 
    integer(ip)           , intent(in)    :: num_neighbours
    integer(ip)           , intent(in)    :: neighbour_ids(num_neighbours)
    integer(ip)           , intent(in)    :: snd_ptrs(num_neighbours+1)
    integer(ieep)         , intent(in)    :: snd_buf(snd_ptrs(num_neighbours+1)-1)   
    integer(ip)           , intent(in)    :: rcv_ptrs(num_neighbours+1)
    integer(ieep)         , intent(out)   :: rcv_buf(rcv_ptrs(num_neighbours+1)-1)
    rcv_buf = snd_buf ! needed to satisfy the intent
  end subroutine serial_context_neighbours_exchange_wo_pack_unpack_ieep

  !=============================================================================
  subroutine serial_context_neighbours_exchange_wo_unpack_ip ( this, &
                                                               num_rcv, list_rcv, rcv_ptrs, rcv_buf, &
                                                               num_snd, list_snd, snd_ptrs, pack_idx,   &
                                                               x, chunk_size)
    implicit none
    class(serial_context_t) , intent(in)    :: this
    ! Control info to receive
    integer(ip)             , intent(in)    :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(out)   :: rcv_buf(:)
    ! Control info to send
    integer(ip)             , intent(in)    :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in)    :: pack_idx (snd_ptrs(num_snd+1)-1)
    ! Raw data to be exchanged
    integer(ip)             , intent(in)    :: x(:)
    integer(ip)   , optional, intent(in)    :: chunk_size
    check(.false.)
  end subroutine serial_context_neighbours_exchange_wo_unpack_ip 

  !=============================================================================
  subroutine serial_context_neighbours_exchange_variable_igp ( this, & 
       &                                              num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       &                                              num_snd, list_snd, snd_ptrs, pack_idx,   &
       &                                              x, y, ptr_chunk_size,mask)
    implicit none
    class(serial_context_t), intent(in)    :: this
    integer(ip)             , intent(in)    :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(in)    :: unpack_idx (rcv_ptrs(num_rcv+1)-1)
    integer(ip)             , intent(in)    :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in)    :: pack_idx (snd_ptrs(num_snd+1)-1)
    integer(igp)            , intent(in)    :: x(:)
    integer(igp)            , intent(inout) :: y(:)
    integer(ip)             , intent(in)    :: ptr_chunk_size(:)
    integer(igp)  , optional, intent(in)    :: mask
    mcheck(.false.,'serial_context_neighbours_exchange_variable_igp is not implemented')
  end subroutine serial_context_neighbours_exchange_variable_igp
  
  !=============================================================================
  subroutine serial_context_neighbours_exchange_variable_ip ( this, & 
       &                                              num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       &                                              num_snd, list_snd, snd_ptrs, pack_idx,   &
       &                                              x, y, ptr_chunk_size,mask)
    implicit none
    class(serial_context_t), intent(in)    :: this
    integer(ip)             , intent(in)    :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(in)    :: unpack_idx (rcv_ptrs(num_rcv+1)-1)
    integer(ip)             , intent(in)    :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in)    :: pack_idx (snd_ptrs(num_snd+1)-1)
    integer(ip)             , intent(in)    :: x(:)
    integer(ip)             , intent(inout) :: y(:)
    integer(ip)             , intent(in)    :: ptr_chunk_size(:)
    integer(ip)   , optional, intent(in)    :: mask
    mcheck(.false.,'serial_context_neighbours_exchange_variable_ip is not implemented')
  end subroutine serial_context_neighbours_exchange_variable_ip
  
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
  subroutine serial_context_bcast_scalar_ip_1D_array ( this, data )
    implicit none
    class(serial_context_t), intent(in)    :: this
    integer(ip)             , intent(inout) :: data(:)
    check(.false.)       ! This routine should be never called
  end subroutine serial_context_bcast_scalar_ip_1D_array

  !=============================================================================
  subroutine serial_context_gather_scalar_igp ( this, input_data, output_data )
    implicit none
    class(serial_context_t), intent(in)   :: this
    integer(igp)         , intent(in)   :: input_data
    integer(igp)         , intent(out)  :: output_data(:)
    check(.false.)       ! This routine should be never called
  end subroutine serial_context_gather_scalar_igp
  
  !=============================================================================
  subroutine serial_context_scatter_scalar_igp ( this, input_data, output_data )
    implicit none
    class(serial_context_t), intent(in)   :: this
    integer(igp)             , intent(in)   :: input_data(:)
    integer(igp)             , intent(out)  :: output_data
    check(.false.)       ! This routine should be never called
  end subroutine serial_context_scatter_scalar_igp

  !=============================================================================
  subroutine serial_context_bcast_scalar_igp ( this, data )
    implicit none
    class(serial_context_t), intent(in)    :: this
    integer(igp)             , intent(inout) :: data
    check(.false.)       ! This routine should be never called
  end subroutine serial_context_bcast_scalar_igp
  
  !=============================================================================
  subroutine serial_context_send_ip ( this, rcv_task, data )
    implicit none
    class(serial_context_t), intent(in)    :: this
    integer(ip)         , intent(in)    :: rcv_task
    integer(ip)         , intent(in)    :: data
    mcheck(.false.,'serial_context_send_ip not implemented')
  end subroutine serial_context_send_ip

  subroutine serial_context_rcv_ip ( this, send_task, data )
    implicit none
    class(serial_context_t) , intent(in)    :: this
    integer(ip)          , intent(inout) :: send_task
    integer(ip)          , intent(inout) :: data
    mcheck(.false.,'serial_context_rcv_ip_ not implemented')
  end subroutine serial_context_rcv_ip

  !=============================================================================
  subroutine serial_context_send_igp ( this, rcv_task, data )
    implicit none
    class(serial_context_t), intent(in)    :: this
    integer(ip)            , intent(in)    :: rcv_task
    integer(igp)           , intent(in)    :: data
    mcheck(.false.,'serial_context_send_ip not implemented')
  end subroutine serial_context_send_igp

  subroutine serial_context_rcv_igp ( this, send_task, data )
    implicit none
    class(serial_context_t) , intent(in)    :: this
    integer(ip)             , intent(inout) :: send_task
    integer(igp)            , intent(inout) :: data
    mcheck(.false.,'serial_context_rcv_igp not implemented')
  end subroutine serial_context_rcv_igp
  
  !=============================================================================
  subroutine serial_context_send_rp ( this, rcv_task, data )
    implicit none
    class(serial_context_t), intent(in)    :: this
    integer(ip)         , intent(in)    :: rcv_task
    real(rp)         , intent(in)    :: data
    mcheck(.false.,'serial_context_send_rp not implemented')
  end subroutine serial_context_send_rp

  subroutine serial_context_rcv_rp ( this, send_task, data )
    implicit none
    class(serial_context_t) , intent(in)    :: this
    integer(ip)          , intent(inout) :: send_task
    real(rp)          , intent(inout) :: data
    mcheck(.false.,'serial_context_rcv_rp_ not implemented')
  end subroutine serial_context_rcv_rp

  !=============================================================================
  subroutine serial_context_send_ip_1D_array ( this, rcv_task, data )
    implicit none
    class(serial_context_t), intent(in)    :: this
    integer(ip)         , intent(in)    :: rcv_task
    integer(ip)         , intent(in)    :: data(:)
    mcheck(.false.,'serial_context_send_ip_1D_array not implemented')
  end subroutine serial_context_send_ip_1D_array

  subroutine serial_context_rcv_ip_1D_array ( this, send_task, data)
    implicit none
    class(serial_context_t) , intent(in)    :: this
    integer(ip)          , intent(inout) :: send_task
    integer(ip)          , intent(inout) :: data(:)
    mcheck(.false.,'serial_context_rcv_ip_1D_array not implemented')
  end subroutine serial_context_rcv_ip_1D_array

  !=============================================================================
  subroutine serial_context_send_igp_1D_array ( this, rcv_task, data )
    implicit none
    class(serial_context_t), intent(in)    :: this
    integer(ip)            , intent(in)    :: rcv_task
    integer(igp)           , intent(in)    :: data(:)
    mcheck(.false.,'serial_context_send_igp_1D_array not implemented')
  end subroutine serial_context_send_igp_1D_array

  subroutine serial_context_rcv_igp_1D_array ( this, send_task, data)
    implicit none
    class(serial_context_t) , intent(in)    :: this
    integer(ip)             , intent(inout) :: send_task
    integer(igp)            , intent(inout) :: data(:)
    mcheck(.false.,'serial_context_rcv_igp_1D_array not implemented')
  end subroutine serial_context_rcv_igp_1D_array
  
  !=============================================================================
  subroutine serial_context_send_rp_1D_array ( this, rcv_task, data )
    implicit none
    class(serial_context_t), intent(in)    :: this
    integer(ip)         , intent(in)    :: rcv_task
    real(rp)         , intent(in)    :: data(:)
    mcheck(.false.,'serial_context_send_rp_1D_array not implemented')
  end subroutine serial_context_send_rp_1D_array

  subroutine serial_context_rcv_rp_1D_array ( this, send_task, data)
    implicit none
    class(serial_context_t) , intent(in)    :: this
    integer(ip)          , intent(inout) :: send_task
    real(rp)          , intent(inout) :: data(:)
    mcheck(.false.,'serial_context_rcv_rp_1D_array not implemented')
  end subroutine serial_context_rcv_rp_1D_array

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
  subroutine serial_context_root_send_master_rcv_rp ( this, input_data, output_data )
    implicit none
    class(serial_context_t), intent(in)      :: this
    real(rp)            , intent(in)      :: input_data
    real(rp)            , intent(inout)   :: output_data
    check(.false.)       ! This routine should be never called
  end subroutine serial_context_root_send_master_rcv_rp

  !=============================================================================
  subroutine serial_context_root_send_master_rcv_rp_1D_array ( this, input_data, output_data )
    implicit none
    class(serial_context_t), intent(in)      :: this
    real(rp)            , intent(in)      :: input_data(:)
    real(rp)            , intent(inout)   :: output_data(:)
    check(.false.)       ! This routine should be never called
  end subroutine serial_context_root_send_master_rcv_rp_1D_array
  
  !=============================================================================
  subroutine serial_context_root_send_master_rcv_logical ( this, input_data, output_data )
    implicit none
    class(serial_context_t), intent(in)      :: this
    logical                , intent(in)      :: input_data
    logical                , intent(inout)   :: output_data
    check(.false.)       ! This routine should be never called
  end subroutine serial_context_root_send_master_rcv_logical 
  
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
  subroutine serial_context_scatter_from_master_ip ( this, input_data, output_data )
    implicit none
    class(serial_context_t), intent(in)   :: this
    integer(ip)            , intent(in)   :: input_data(:)
    integer(ip)            , intent(out)  :: output_data
    check(.false.)       ! This routine should be never called
    output_data = 0_ip  ! needed to satisfy the intent
  end subroutine serial_context_scatter_from_master_ip

  !=============================================================================
  subroutine serial_context_scatter_from_masterv_ip_1D_array ( this, input_data, send_counts, displs, output_data_size, output_data )
    implicit none
    class(serial_context_t), intent(in)   :: this
    integer(ip)            , intent(in)   :: input_data(:)
    integer(ip)            , intent(in)   :: send_counts(:)
    integer(ip)            , intent(in)   :: displs(:)
    integer(ip)            , intent(in)   :: output_data_size
    integer(ip)            , intent(out)  :: output_data(output_data_size)
    check(.false.)       ! This routine should be never called
    output_data = 0.0_rp ! needed to satisfy the intent
  end subroutine serial_context_scatter_from_masterv_ip_1D_array
  
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
