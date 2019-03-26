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
!
! Hybrid execution context: 
! =========================
! This implementation is provided to make hybrid execution portable (to any
! system having mpi and openmp). More sophisticated implementations are 
! posssible, see http://charm.cs.illinois.edu/manuals/html/ampi/manual.html.
! 
! There are two types of mpi_omp_contexts: those that contain all the threads 
! of the current rank, called homogeneous and those that contain some threads
! of the current rank, called heterogeneous. 
!
! We assume that homogeneous contexts are used for neighbor exchanges and 
! collectives not involving the master task (see comments in exectution_context.f90). 
! In this case mpi_omp_context_root_rank is in charge of performing collectives.
! and further inter-rank work may be required.
!
! In the case of collectives involving the master task, we allow for heterogeneous
! context (for the master side). Typically, one rank will contain several master
! tasks mapped to different threads. We therefore need to allow them to call mpi.
!
module mpi_omp_context_names
  ! Serial modules
  use types_names
  use memor_names

  ! Parallel modules
  use execution_context_names
  !$ use omp_lib
#ifdef MPI_MOD
  use mpi
#endif
  implicit none 
#ifdef MPI_H
  include 'mpif.h'
#endif
#include "debug.i90"
  private

  ! Global shared variables to implement communications
  logical     , allocatable :: lg2_buffer(:,:)
  integer(ip) , allocatable :: ip1_buffer(:)
  integer(ip) , allocatable :: ip1_v_buffer(:)
  integer(igp), allocatable :: igp1_buffer(:)
  integer(igp), allocatable :: igp1_v_buffer(:)
  integer(ip) , allocatable :: ip2_buffer(:,:)
  real(rp)    , allocatable :: rp1_buffer(:)
  real(rp)    , allocatable :: rp1_v_buffer(:)
  real(rp)    , allocatable :: rp2_buffer(:,:)
  
  ! Constants that define buffer types. They must be changed if
  ! our definitions in types.f90 changes.
  integer, parameter :: mpi_omp_context_ieep = mpi_integer1
  integer, parameter :: mpi_omp_context_ip   = mpi_integer
  integer, parameter :: mpi_omp_context_igp  = mpi_integer8
  integer, parameter :: mpi_omp_context_rp   = mpi_double_precision
  integer, parameter :: mpi_omp_context_lg   = mpi_logical
  integer, parameter :: mpi_omp_context_root_rank   = 0
  integer, parameter :: mpi_omp_context_default_root_thread = 0
  integer, parameter :: mpi_omp_context_tag  = 1453524 ! which number should go here?

  integer, parameter :: mpi_omp_context_homogeneous   = 0
  integer, parameter :: mpi_omp_context_heterogeneous = 1
  
  
  ! Parallel context
  type, extends(execution_context_t) :: mpi_omp_context_t
     private 
     ! ***IMPORTANT NOTE***: parallel contexts are always of type 
     ! integer: the kind parameter must NOT be specified. This requirement is 
     ! imposed by the underlying message-passing library, i.e., MPI. 
     ! The same comment applies to other integers in mpi interfaces (except
     ! buffers).

     logical :: created_from_mpi = .false.
     integer :: type       = mpi_omp_context_homogeneous
     integer :: icontxt    = mpi_comm_null
     integer :: num_ranks      = -1
     integer :: current_rank   = -1
     !integer :: master_rank    = -1
     integer :: num_threads    = -1
     integer :: current_thread = -1
     integer :: root_thread    = mpi_omp_context_default_root_thread
     integer :: max_num_threads = -1
     integer :: min_num_threads = -1
     !integer, allocatable :: threads_x_rank(:)
     !integer :: master_thread  = -1
     real(rp) :: time_exchange
     real(rp) :: time_bcasts
     real(rp) :: time_reduction
     real(rp) :: time_gather
     real(rp) :: time_scatter
   contains
     ! These functions should be non_overridable but there is a bug in gfotran
     procedure :: create             => mpi_omp_context_create
     procedure :: assign             => mpi_omp_context_assign
     procedure :: split_by_condition => mpi_omp_context_split_by_condition
     procedure :: split_by_color     => mpi_omp_context_split_by_color
     procedure :: free               => mpi_omp_context_free
     procedure :: nullify            => mpi_omp_context_nullify
     procedure :: report_times       => mpi_omp_context_report_times
     procedure :: am_i_member        => mpi_omp_context_am_i_member
     procedure :: am_i_root          => mpi_omp_context_am_i_root
     procedure :: barrier            => mpi_omp_context_barrier
     procedure :: time               => mpi_omp_context_time
     procedure :: sum_scalar_rp      => mpi_omp_context_sum_scalar_rp
     procedure :: sum_vector_rp      => mpi_omp_context_sum_vector_rp
     procedure :: max_scalar_rp      => mpi_omp_context_max_scalar_rp
     procedure :: max_vector_rp      => mpi_omp_context_max_vector_rp
     procedure :: min_scalar_rp      => mpi_omp_context_min_scalar_rp
     procedure :: max_scalar_ip      => mpi_omp_context_max_scalar_ip
     procedure :: sum_scalar_igp     => mpi_omp_context_sum_scalar_igp
     procedure :: sum_vector_igp     => mpi_omp_context_sum_vector_igp
     procedure :: scatter_ip         => mpi_omp_context_scatter_scalar_ip
     procedure :: gather_ip          => mpi_omp_context_gather_scalar_ip
     procedure :: bcast_ip           => mpi_omp_context_bcast_scalar_ip
     procedure :: bcast_ip_1D_array  => mpi_omp_context_bcast_scalar_ip_1D_array
     procedure :: scatter_igp        => mpi_omp_context_scatter_scalar_igp
     procedure :: gather_igp         => mpi_omp_context_gather_scalar_igp
     procedure :: bcast_igp          => mpi_omp_context_bcast_scalar_igp
     procedure :: bcast_subcontext   => mpi_omp_context_bcast_subcontext
     procedure :: neighbours_exchange_rp                   => mpi_omp_context_neighbours_exchange_rp  
     procedure :: neighbours_exchange_wo_alpha_beta_rp     => mpi_omp_context_neighbours_exchange_wo_alpha_beta_rp
     procedure :: neighbours_exchange_wo_alpha_beta_rp_v   => mpi_omp_context_neighbours_exchange_wo_alpha_beta_rp_v
     procedure :: neighbours_exchange_ip                   => mpi_omp_context_neighbours_exchange_ip                 
     procedure :: neighbours_exchange_igp                  => mpi_omp_context_neighbours_exchange_igp                
     procedure :: neighbours_exchange_single_ip            => mpi_omp_context_neighbours_exchange_single_ip          
     procedure :: neighbours_exchange_wo_pack_unpack_ieep  => mpi_omp_context_neighbours_exchange_wo_pack_unpack_ieep
     procedure :: neighbours_exchange_wo_unpack_ip         => mpi_omp_context_neighbours_exchange_wo_unpack_ip
     procedure :: neighbours_exchange_variable_igp         => mpi_omp_context_neighbours_exchange_variable_igp       
     procedure :: neighbours_exchange_variable_ip          => mpi_omp_context_neighbours_exchange_variable_ip 
     procedure :: send_ip           => mpi_omp_context_send_ip
     procedure :: send_igp          => mpi_omp_context_send_igp     
     procedure :: send_rp           => mpi_omp_context_send_rp
     procedure :: send_ip_1D_array  => mpi_omp_context_send_ip_1D_array
     procedure :: send_igp_1D_array => mpi_omp_context_send_igp_1D_array     
     procedure :: send_rp_1D_array  => mpi_omp_context_send_rp_1D_array
     procedure :: rcv_ip            => mpi_omp_context_rcv_ip
     procedure :: rcv_igp           => mpi_omp_context_rcv_igp     
     procedure :: rcv_rp            => mpi_omp_context_rcv_rp     
     procedure :: rcv_ip_1D_array   => mpi_omp_context_rcv_ip_1D_array
     procedure :: rcv_igp_1D_array  => mpi_omp_context_rcv_igp_1D_array     
     procedure :: rcv_rp_1D_array   => mpi_omp_context_rcv_rp_1D_array
     procedure :: root_send_master_rcv_ip          => mpi_omp_context_root_send_master_rcv_ip
     procedure :: root_send_master_rcv_ip_1D_array => mpi_omp_context_root_send_master_rcv_ip_1D_array
     procedure :: root_send_master_rcv_rp          => mpi_omp_context_root_send_master_rcv_rp
     procedure :: root_send_master_rcv_rp_1D_array => mpi_omp_context_root_send_master_rcv_rp_1D_array
     procedure :: root_send_master_rcv_logical     => mpi_omp_context_root_send_master_rcv_logical
     procedure :: gather_to_master_ip              => mpi_omp_context_gather_to_master_ip            
     procedure :: gather_to_master_igp             => mpi_omp_context_gather_to_master_igp           
     procedure :: gather_to_master_ip_1D_array     => mpi_omp_context_gather_to_master_ip_1D_array   
     procedure :: gather_to_masterv_ip_1D_array    => mpi_omp_context_gather_to_masterv_ip_1D_array  
     procedure :: gather_to_masterv_igp_1D_array   => mpi_omp_context_gather_to_masterv_igp_1D_array 
     procedure :: gather_to_masterv_rp_1D_array    => mpi_omp_context_gather_to_masterv_rp_1D_array  
     procedure :: gather_to_masterv_rp_2D_array    => mpi_omp_context_gather_to_masterv_rp_2D_array  
     procedure :: scatter_from_master_ip           => mpi_omp_context_scatter_from_master_ip
     procedure :: scatter_from_masterv_ip_1D_array => mpi_omp_context_scatter_from_masterv_ip_1D_array
     procedure :: scatter_from_masterv_rp_1D_array => mpi_omp_context_scatter_from_masterv_rp_1D_array
  end type mpi_omp_context_t

  ! Types
  public :: mpi_omp_context_t
  
  interface
     subroutine report_bindings(Fcomm) bind(c,name='report_bindings')
       use iso_c_binding
       implicit none
       integer, value, intent(in) :: Fcomm
     end subroutine report_bindings
  end interface

contains

  !=============================================================================
  ! This function creates an homogenous context duplicating mpi_comm_world
  subroutine mpi_omp_context_create ( this )
    implicit none 
    class(mpi_omp_context_t), intent(inout) :: this
    integer :: current_task, num_tasks, istat, i
    logical :: initialized
    integer :: provided_thread_support, comm, num_ranks_with_min_threads
    !integer, allocatable :: threads_x_rank(:)
    
    call this%free(finalize=.false.)

    this%num_threads    = 1
    this%current_thread = 0
    !$ this%num_threads    = omp_get_num_threads()
    !$ this%current_thread = omp_get_thread_num()
    this%root_thread = mpi_omp_context_default_root_thread
    !write(*,*) 'Before mpi_init',this%current_thread,this%num_threads

    if(this%current_thread==this%root_thread) then
       ! Init mpi and dup mpi_comm_world
       !write(*,*) 'Calling mpi_init',this%current_thread,this%num_threads
       call mpi_initialized(initialized,istat); check((.not.initialized).and.(istat == mpi_success))
       call mpi_init_thread(mpi_thread_multiple, provided_thread_support, istat); check(istat == mpi_success)
       check(provided_thread_support==mpi_thread_multiple)
       this%created_from_mpi = .true.
       assert(comm/=mpi_comm_null)
       call mpi_comm_dup(mpi_comm_world,this%icontxt,istat); check(istat == mpi_success)
       
       ! We assume that the number of threads is the same in all ranks except the last one
       call mpi_allreduce(this%num_threads,this%max_num_threads,1,mpi_omp_context_ip,mpi_max,this%icontxt,istat); check(istat == mpi_success)
       call mpi_allreduce(this%num_threads,this%min_num_threads,1,mpi_omp_context_ip,mpi_min,this%icontxt,istat); check(istat == mpi_success)
       !assert(this%min_num_threads==1)
       assert(this%num_threads==this%max_num_threads.or.this%num_threads==this%min_num_threads)
       if(this%num_threads==this%min_num_threads) then
          call mpi_allreduce(1,num_ranks_with_min_threads,1,mpi_omp_context_ip,mpi_sum,this%icontxt,istat); check(istat == mpi_success)
       else
          call mpi_allreduce(0,num_ranks_with_min_threads,1,mpi_omp_context_ip,mpi_sum,this%icontxt,istat); check(istat == mpi_success)
       end if
       assert(num_ranks_with_min_threads==1.or.this%max_num_threads==this%min_num_threads)

       ! Allocate buffers to manage intrarank (interthread) communications
       call memalloc(2,this%max_num_threads,lg2_buffer,__FILE__,__LINE__,lb2=0)
       call memalloc(  this%max_num_threads+1,ip1_buffer,__FILE__,__LINE__,lb1=0)
       call memalloc(  this%max_num_threads,ip1_v_buffer,__FILE__,__LINE__,lb1=0)
       call memalloc(  this%max_num_threads,igp1_buffer,__FILE__,__LINE__,lb1=0)
       call memalloc(  this%max_num_threads,igp1_v_buffer,__FILE__,__LINE__,lb1=0)
       call memalloc(4,this%max_num_threads,ip2_buffer,__FILE__,__LINE__,lb2=0)
       call memalloc(  this%max_num_threads,rp1_buffer,__FILE__,__LINE__,lb1=0)
       call memalloc(  this%max_num_threads,rp1_v_buffer,__FILE__,__LINE__,lb1=0)
       call memalloc(1,this%max_num_threads,rp2_buffer,__FILE__,__LINE__,lb2=0)

       ip1_buffer = this%icontxt
       ip2_buffer(1,:)=this%max_num_threads
       ip2_buffer(2,:)=this%min_num_threads ! This is not needed, I leave it just in case
    else
       this%created_from_mpi = .false.
    end if
    this%type = mpi_omp_context_homogeneous
    !$OMP BARRIER

    this%icontxt = ip1_buffer(this%current_thread)
    call mpi_comm_size(this%icontxt,this%num_ranks,istat)    ; check(istat == mpi_success)
    call mpi_comm_rank(this%icontxt,this%current_rank,istat) ; check(istat == mpi_success)
    this%max_num_threads = ip2_buffer(1,this%current_thread)
    this%min_num_threads = ip2_buffer(2,this%current_thread)

    call this%set_current_task( this%max_num_threads*this%current_rank + this%current_thread )
    call this%set_num_tasks   ( this%max_num_threads*(this%num_ranks-1)+ this%min_num_threads )
    
    !write(*,*) 'After task def',this%get_current_task(),this%get_num_tasks()
!!!#ifdef DEBUG
    call report_bindings(this%icontxt)
!!!#endif
    this%time_exchange  = 0.0_rp
    this%time_bcasts    = 0.0_rp
    this%time_reduction = 0.0_rp
    this%time_gather    = 0.0_rp
    this%time_scatter   = 0.0_rp

  end subroutine mpi_omp_context_create

  !=============================================================================
  ! This function creates an homogenous context by copy
  subroutine mpi_omp_context_assign(this, that)
    implicit none 
    class(mpi_omp_context_t)  , intent(inout) :: this
    class(execution_context_t), intent(in)    :: that
    integer  :: current_task, num_tasks, istat
    call this%free(finalize=.false.)
    select type(that)
    type is(mpi_omp_context_t)
       ! Uncomment the following line for maximum checking
       ! call mpi_initialized(initialized,info); check((initialized).and.(info == mpi_success))
       assert(that%icontxt/=mpi_comm_null)
       this%created_from_mpi = .false.
       this%type = that%type
       this%icontxt = that%icontxt
       this%root_thread = that%root_thread
       this%num_threads    = that%num_threads
       this%current_thread = that%current_thread
       this%max_num_threads = that%max_num_threads
       this%min_num_threads = that%min_num_threads
       assert(this%icontxt/=mpi_comm_null)
       call mpi_comm_size(this%icontxt,this%num_ranks,istat)    ; check(istat == mpi_success)
       call mpi_comm_rank(this%icontxt,this%current_rank,istat) ; check(istat == mpi_success)
       assert(this%num_ranks    == that%num_ranks)
       assert(this%current_rank == that%current_rank)
       call this%set_current_task(that%get_current_task())
       call this%set_num_tasks(that%get_num_tasks())
       this%time_exchange  = that%time_exchange  
       this%time_bcasts    = that%time_bcasts    
       this%time_reduction = that%time_reduction 
       this%time_gather    = that%time_gather    
       this%time_scatter   = that%time_scatter   
    class default
       check(.false.)
    end select
  end subroutine mpi_omp_context_assign

  !=============================================================================
  subroutine mpi_omp_context_split_by_color ( this, color, new_subcontext )
    implicit none 
    class(mpi_omp_context_t)  , intent(in)    :: this
    integer                   , intent(in)    :: color
    class(execution_context_t), allocatable , intent(inout) :: new_subcontext
    integer, parameter :: key=0
    integer :: i,j, istat, my_color, current_rank, num_ranks
    logical :: is_homogeneous, all_are_homogeneous


    if(allocated(new_subcontext)) then
       call new_subcontext%free(finalize=.false.)
    else
       allocate(new_subcontext,mold=this,stat=istat);check(istat==0)
    end if

    if(color==undefined_color) then
       my_color = mpi_undefined
    else
       my_color = color
    end if

    ! When num_threads==max_num_threads, we assume homogeneous (only one color) or 
    ! heterogeneous (all colors different) ranks.
    ! 
    ! When there are heterogeneous ranks we need max_num_threads invocations to 
    ! mpi_comm_split. In the heterogeneous case, these calls will be performed by 
    ! each thread and this rank will belong to max_num_threads comms. In the 
    ! homogeneous case the rank will belong to the communicator that has the maximum 
    ! number of ranks.
    !
    ! When all ranks are homogeneous only one call to mpi_comm_split is enough.
    !
    ! 1) Decide if the ranks are homogeneous or heterogeneous.
    ip1_buffer(this%current_thread)=my_color
    !$OMP BARRIER
    if(this%current_thread==this%root_thread) then
       if(ip1_buffer(this%num_threads-1)/=ip1_buffer(this%root_thread)) then
          is_homogeneous = .false.
          ! Check the rest of the threads do not repeat
          do i=0,this%num_threads-1
             do j=0,this%num_threads-1
                if(i/=j) then
                   check(ip1_buffer(i)/=ip1_buffer(this%root_thread))
                end if
             end do
          end do
       else
          is_homogeneous = .true.
          ! Check the rest of the threads are all equal
          do i=0,this%num_threads-1
             check(ip1_buffer(i)==ip1_buffer(this%root_thread))
          end do
       end if
       call mpi_allreduce(is_homogeneous,all_are_homogeneous,1,mpi_omp_context_lg,mpi_land,this%icontxt,istat); check(istat == mpi_success)
       lg2_buffer(1,:) = is_homogeneous
       lg2_buffer(2,:) = all_are_homogeneous 
    end if
    !$OMP BARRIER
    is_homogeneous      = lg2_buffer(1,this%current_thread)
    all_are_homogeneous = lg2_buffer(2,this%current_thread)

    select type(new_subcontext)
    type is(mpi_omp_context_t)

       new_subcontext%num_threads    = this%num_threads     ! =omp_get_num_threads()
       new_subcontext%current_thread = this%current_thread  ! =omp_get_thread_num()

       !2) Create communicators by split
       if(all_are_homogeneous) then ! only one call to mpi_comm_split and shares the comm
          new_subcontext%created_from_mpi = .false.
          if(this%current_thread==this%root_thread) then
             call mpi_comm_split(this%icontxt, my_color, key, new_subcontext%icontxt, istat);  assert ( istat == mpi_success )
             assert((my_color == mpi_undefined).or.(new_subcontext%icontxt/=mpi_comm_null))
             ip1_buffer(:) = new_subcontext%icontxt
             new_subcontext%created_from_mpi = .true.
          end if
          !$OMP BARRIER
          new_subcontext%icontxt = ip1_buffer(new_subcontext%current_thread)
          if(new_subcontext%icontxt/=mpi_comm_null) then
             call mpi_comm_rank(new_subcontext%icontxt,new_subcontext%current_rank,istat) ; check(istat == mpi_success)
             call mpi_comm_size(new_subcontext%icontxt,new_subcontext%num_ranks,istat)    ; check(istat == mpi_success)
          end if
       else
          !$OMP CRITICAL
          call mpi_comm_split(this%icontxt, my_color, key, new_subcontext%icontxt, istat);  assert ( istat == mpi_success )
          !$OMP END CRITICAL
          assert((my_color == mpi_undefined).or.(new_subcontext%icontxt/=mpi_comm_null))
          new_subcontext%created_from_mpi = .true.

          !3) Select appropriate communicators in homogeneous ranks
          !$OMP BARRIER
          if(is_homogeneous) then
             new_subcontext%type=mpi_omp_context_homogeneous
             ! Select the communicator having the largest set.
             !call mpi_comm_rank(new_subcontext%icontxt,current_rank,istat) ; check(istat == mpi_success)
             call mpi_comm_size(new_subcontext%icontxt,num_ranks,istat)    ; check(istat == mpi_success)
             ip1_buffer(this%current_thread)=num_ranks
             !$OMP BARRIER
             if(this%current_thread==this%root_thread) then
                num_ranks = 0
                do i=0,this%num_threads-1
                   if(num_ranks < ip1_buffer(i)) then
                      num_ranks = ip1_buffer(i)
                      j = i
                   end if
                end do
                ip1_buffer = j
             end if
             !$OMP BARRIER
             if(this%current_thread/=ip1_buffer(this%current_thread)) then
                call mpi_comm_free(new_subcontext%icontxt,istat); check(istat == mpi_success)
                new_subcontext%created_from_mpi = .false.
             else
                call mpi_comm_rank(new_subcontext%icontxt,current_rank,istat) ; check(istat == mpi_success)
                call mpi_comm_size(new_subcontext%icontxt,num_ranks,istat)    ; check(istat == mpi_success)
                ip2_buffer(1,:) = current_rank
                ip2_buffer(2,:) = num_ranks
                !call mpi_bcast(this%master_thread,1,mpi_omp_context_ip,num_ranks-1,this%icontxt,istat); check( istat == mpi_success )
             end if
             !$OMP BARRIER
             new_subcontext%current_rank = ip2_buffer(1,new_subcontext%current_thread)
             new_subcontext%num_ranks    = ip2_buffer(2,new_subcontext%current_thread)
             new_subcontext%root_thread  = ip1_buffer(this%current_thread)
          else
             new_subcontext%type=mpi_omp_context_heterogeneous
             call mpi_comm_rank(new_subcontext%icontxt,new_subcontext%current_rank,istat) ; check(istat == mpi_success)
             call mpi_comm_size(new_subcontext%icontxt,new_subcontext%num_ranks,istat)    ; check(istat == mpi_success)
             new_subcontext%root_thread  = mpi_omp_context_default_root_thread
             ! bcast master rank
             !call mpi_bcast(this%current_thread,1,mpi_omp_context_ip,this%current_rank,this%icontxt,istat); check( istat == mpi_success )
          end if

       end if

       if(new_subcontext%icontxt/=mpi_comm_null) then
          if(this%current_thread==this%root_thread) then
             call mpi_allreduce(new_subcontext%num_threads,new_subcontext%max_num_threads,1,mpi_omp_context_ip,mpi_max,new_subcontext%icontxt,istat)
             call mpi_allreduce(new_subcontext%num_threads,new_subcontext%min_num_threads,1,mpi_omp_context_ip,mpi_min,new_subcontext%icontxt,istat)
             ip2_buffer(1,:) = new_subcontext%max_num_threads
             ip2_buffer(2,:) = new_subcontext%min_num_threads
          end if
          !$OMP BARRIER
          new_subcontext%max_num_threads = ip2_buffer(1,new_subcontext%current_thread)
          new_subcontext%min_num_threads = ip2_buffer(2,new_subcontext%current_thread)

          !write(*,*) 'In mpi_omp', new_subcontext%current_rank, new_subcontext%num_ranks, new_subcontext%current_thread, new_subcontext%num_threads, new_subcontext%max_num_threads, new_subcontext%min_num_threads

          call new_subcontext%set_current_task( new_subcontext%max_num_threads*new_subcontext%current_rank + new_subcontext%current_thread )
          call new_subcontext%set_num_tasks   ( new_subcontext%max_num_threads*(new_subcontext%num_ranks-1)+ new_subcontext%min_num_threads )

          new_subcontext%time_exchange  = 0.0_rp
          new_subcontext%time_bcasts    = 0.0_rp
          new_subcontext%time_reduction = 0.0_rp
          new_subcontext%time_gather    = 0.0_rp
          new_subcontext%time_scatter   = 0.0_rp

       end if
       class default
       check(.false.)
    end select

  end subroutine mpi_omp_context_split_by_color

  !=============================================================================
  subroutine mpi_omp_context_split_by_condition ( this, in_subcontext1, subcontext1, subcontext2 )
    implicit none 
    class(mpi_omp_context_t)  , intent(in)    :: this
    logical                   , intent(in)    :: in_subcontext1
    class(execution_context_t), allocatable, intent(inout) :: subcontext1
    class(execution_context_t), allocatable, intent(inout) :: subcontext2
    integer                :: istat, current_rank, num_ranks, current_thread, num_threads, i
    integer, parameter     :: key=0

    if(allocated(subcontext1)) then
       call subcontext1%free(finalize=.false.)
    else
       allocate(subcontext1,mold=this,stat=istat);check(istat==0)
    end if
    if(allocated(subcontext2)) then
       call subcontext2%free(finalize=.false.)
    else
       allocate(subcontext2,mold=this,stat=istat);check(istat==0)
    end if

    ! Check that all threads in this rank have the same condition (so they are in the same context)
    if(in_subcontext1) then
       ip1_buffer(this%current_thread)=1
    else
       ip1_buffer(this%current_thread)=0
    end if
    !$OMP BARRIER
    do i=0,this%num_threads-1
       check(ip1_buffer(i)==ip1_buffer(this%root_thread))
    end do

    current_thread = this%root_thread
    num_threads    = 1
    !$ current_thread = omp_get_thread_num()
    !$ num_threads    = omp_get_num_threads()
    
    
    select type(subcontext1)
    type is(mpi_omp_context_t)
       select type(subcontext2)
       type is(mpi_omp_context_t)

          subcontext1%created_from_mpi = .false.
          subcontext2%created_from_mpi = .false.

          if( in_subcontext1 ) then
             subcontext1%current_thread = current_thread
             subcontext1%num_threads    = num_threads
             if(this%current_thread==this%root_thread) then
                call mpi_comm_split(this%icontxt, 1, key, subcontext1%icontxt, istat); assert ( istat == mpi_success )
                assert(subcontext1%icontxt/=mpi_comm_null)
                subcontext1%created_from_mpi = .true. 
                ip1_buffer(:) = subcontext1%icontxt
                call mpi_allreduce(subcontext1%num_threads,subcontext1%max_num_threads,1,mpi_omp_context_ip,mpi_max,subcontext1%icontxt,istat)
                call mpi_allreduce(subcontext1%num_threads,subcontext1%min_num_threads,1,mpi_omp_context_ip,mpi_min,subcontext1%icontxt,istat)
                ip2_buffer(1,:) = subcontext1%max_num_threads
                ip2_buffer(2,:) = subcontext1%min_num_threads
             end if
             !$OMP BARRIER
             subcontext1%icontxt = ip1_buffer(current_thread)
             subcontext1%max_num_threads = ip2_buffer(1,current_thread)
             subcontext1%min_num_threads = ip2_buffer(2,current_thread)
             call mpi_comm_rank(subcontext1%icontxt,subcontext1%current_rank,istat) ; check(istat == mpi_success)
             call mpi_comm_size(subcontext1%icontxt,subcontext1%num_ranks,istat)    ; check(istat == mpi_success)
             call subcontext1%set_current_task( subcontext1%max_num_threads*subcontext1%current_rank + subcontext1%current_thread )
             call subcontext1%set_num_tasks   ( subcontext1%max_num_threads*(subcontext1%num_ranks-1)+ subcontext1%min_num_threads )
             subcontext1%time_exchange  = 0.0_rp
             subcontext1%time_bcasts    = 0.0_rp
             subcontext1%time_reduction = 0.0_rp
             subcontext1%time_gather    = 0.0_rp
             subcontext1%time_scatter   = 0.0_rp

             call subcontext2%nullify()             
          else
             subcontext2%current_thread = current_thread
             subcontext2%num_threads    = num_threads
             if(this%current_thread==this%root_thread) then
                call mpi_comm_split(this%icontxt, 2, key, subcontext2%icontxt, istat); assert ( istat == mpi_success )
                subcontext2%created_from_mpi = .true. 
                assert(subcontext2%icontxt/=mpi_comm_null)
                ip1_buffer(:) = subcontext2%icontxt
                call mpi_allreduce(subcontext2%num_threads,subcontext2%max_num_threads,1,mpi_omp_context_ip,mpi_max,subcontext2%icontxt,istat)
                call mpi_allreduce(subcontext2%num_threads,subcontext2%min_num_threads,1,mpi_omp_context_ip,mpi_min,subcontext2%icontxt,istat)
                ip2_buffer(1,:) = subcontext2%max_num_threads
                ip2_buffer(2,:) = subcontext2%min_num_threads
             end if
             !$OMP BARRIER
             subcontext2%icontxt = ip1_buffer(current_thread)
             subcontext2%max_num_threads = ip2_buffer(1,current_thread)
             subcontext2%min_num_threads = ip2_buffer(2,current_thread)
             call mpi_comm_rank(subcontext2%icontxt,subcontext2%current_rank,istat) ; check(istat == mpi_success)
             call mpi_comm_size(subcontext2%icontxt,subcontext2%num_ranks,istat)    ; check(istat == mpi_success)
             call subcontext2%set_current_task( subcontext2%max_num_threads*subcontext2%current_rank + subcontext2%current_thread )
             call subcontext2%set_num_tasks   ( subcontext2%max_num_threads*(subcontext2%num_ranks-1)+ subcontext2%min_num_threads )
             subcontext2%time_exchange  = 0.0_rp
             subcontext2%time_bcasts    = 0.0_rp
             subcontext2%time_reduction = 0.0_rp
             subcontext2%time_gather    = 0.0_rp
             subcontext2%time_scatter   = 0.0_rp

             call subcontext1%nullify()             
          end if

       class default
          check(.false.)
       end select
    class default
       check(.false.)
    end select
    
  end subroutine mpi_omp_context_split_by_condition

  !=============================================================================
  subroutine mpi_omp_context_free ( this, finalize  )
    implicit none 
    class(mpi_omp_context_t), intent(inout) :: this
    logical                 , intent(in)    :: finalize
    integer :: istat

    if(this%icontxt/=mpi_comm_null) call report_bindings(this%icontxt) 

    !$OMP BARRIER
    if(this%created_from_mpi) then
       if(this%icontxt/=mpi_comm_null.and.this%icontxt/=mpi_comm_world) then
          call mpi_comm_free(this%icontxt,istat); check(istat == mpi_success)
       end if
       if(finalize) then
          call memfree(lg2_buffer,__FILE__,__LINE__   )
          call memfree(ip1_buffer,__FILE__,__LINE__   )
          call memfree(ip1_v_buffer,__FILE__,__LINE__ )
          call memfree(igp1_buffer,__FILE__,__LINE__  )
          call memfree(igp1_v_buffer,__FILE__,__LINE__)
          call memfree(ip2_buffer,__FILE__,__LINE__   )
          call memfree(rp1_buffer,__FILE__,__LINE__   )
          call memfree(rp1_v_buffer,__FILE__,__LINE__ )
          call memfree(rp2_buffer,__FILE__,__LINE__   )
          call mpi_finalize(istat); check(istat == mpi_success)
       end if
    end if

    this%created_from_mpi=.false.
    this%icontxt=mpi_comm_null
    this%num_ranks      = -1
    this%current_rank   = -1
    this%num_threads    = -1
    this%current_thread = -1
    call this%set_current_task(-1)
    call this%set_num_tasks(-1)

  end subroutine mpi_omp_context_free

  subroutine mpi_omp_context_report_times ( this, show_header, luout )
    implicit none 
    class(mpi_omp_context_t), intent(inout) :: this
    logical, intent(in), optional      :: show_header 
    integer(ip), intent(in), optional  :: luout
      
      ! Locals
      character(len=*), parameter    :: fmt_header = '(a25,1x,3(2x,a15),3(2x,a15))'
      character(len=*), parameter    :: fmt_data   = '(a25,1x,3(2x,es15.9),3(2x,es15.9))'
      real(rp)                       :: accum_max, accum_min, accum_sum
      logical                        :: show_header_

      accum_max = this%time_exchange
      accum_min = this%time_exchange
      accum_sum = this%time_exchange
      call this%max_scalar_rp(accum_max)
      call this%min_scalar_rp(accum_min)
      call this%sum_scalar_rp(accum_sum)

      show_header_ = .true.
      if (present(show_header)) show_header_ = show_header 

      if ( show_header_ ) then 
         if (present(luout) ) then
            if ( this%am_i_root() ) write(luout,fmt_header) '', 'Min (secs.)', 'Max (secs.)', 'Avg (secs.)'
         else
            if ( this%am_i_root() ) write(*,fmt_header) '', 'Min (secs.)', 'Max (secs.)', 'Avg (secs.)'
         end if
      end if         
      
      if (present(luout) ) then
         if ( this%am_i_root() ) write(luout,fmt_data) 'NEIGHBOURS EXCHANGE', accum_min, accum_max, accum_sum/this%get_num_tasks()
      else
         if ( this%am_i_root() ) write(*,fmt_data) 'NEIGHBOURS EXCHANGE', accum_min, accum_max, accum_sum/this%get_num_tasks()
      end if

      this%time_exchange = 0.0_rp
      
    end subroutine mpi_omp_context_report_times


  !=============================================================================
  subroutine mpi_omp_context_nullify ( this )
    implicit none 
    class(mpi_omp_context_t), intent(inout) :: this
    call this%free(finalize=.false.)
    this%icontxt = mpi_comm_null
    this%num_ranks      = -1
    this%current_rank   = -1
    this%num_threads    = -1
    this%current_thread = -1
    call this%set_current_task(-1)
    call this%set_num_tasks(-1)
    this%time_exchange  = 0.0_rp
    this%time_bcasts    = 0.0_rp
    this%time_reduction = 0.0_rp
    this%time_gather    = 0.0_rp
    this%time_scatter   = 0.0_rp
  end subroutine mpi_omp_context_nullify

  !=============================================================================
  pure function mpi_omp_context_am_i_member(this)
    implicit none
    class(mpi_omp_context_t), intent(in) :: this
    logical                              :: mpi_omp_context_am_i_member
    mpi_omp_context_am_i_member = (this%get_current_task()>=0)
  end function mpi_omp_context_am_i_member

  !=============================================================================
  pure function mpi_omp_context_am_i_root(this)
    implicit none
    class(mpi_omp_context_t), intent(in) :: this
    logical                          :: mpi_omp_context_am_i_root
    mpi_omp_context_am_i_root = (this%get_current_task()==0)
  end function mpi_omp_context_am_i_root

  !=============================================================================
  subroutine mpi_omp_context_barrier(this)
    implicit none 
    class(mpi_omp_context_t), intent(in) :: this
    integer :: istat
    if(this%current_thread==this%root_thread) then
       call mpi_barrier ( this%icontxt, istat); check ( istat == mpi_success )
    end if
    !$OMP BARRIER
  end subroutine mpi_omp_context_barrier

  !=============================================================================
  function mpi_omp_context_time(this)
    implicit none 
    class(mpi_omp_context_t), intent(in) :: this
    real(rp) :: mpi_omp_context_time
    mpi_omp_context_time =  mpi_wtime ()
  end function mpi_omp_context_time

  !=============================================================================
  subroutine mpi_omp_context_sum_scalar_rp (this,alpha)
    implicit none
    class(mpi_omp_context_t) , intent(in)    :: this
    real(rp)                 , intent(inout) :: alpha
    integer  :: istat,i
    rp1_buffer(this%current_thread) = alpha
    !$OMP BARRIER
    if(this%current_thread==this%root_thread) then
       do i=0,this%num_threads-1
          if(i/=this%root_thread) then
             rp1_buffer(this%root_thread) = rp1_buffer(this%root_thread) + rp1_buffer(i)
          end if
       end  do
       call mpi_allreduce(rp1_buffer(this%root_thread),alpha,1,mpi_omp_context_rp,mpi_sum,this%icontxt,istat); check ( istat == mpi_success )
       rp1_buffer=alpha
    end if
    !$OMP BARRIER
    alpha = rp1_buffer(this%current_thread)
  end subroutine mpi_omp_context_sum_scalar_rp

  !=============================================================================
  subroutine mpi_omp_context_sum_vector_rp(this,alpha)
    implicit none
    class(mpi_omp_context_t) , intent(in)    :: this
    real(rp)                 , intent(inout) :: alpha(:) 
    integer :: i, istat

    if(size(alpha)>size(rp2_buffer,1)) then
       !$OMP BARRIER
       if(this%current_thread==this%root_thread) call memrealloc(size(alpha),this%num_threads,rp2_buffer,__FILE__,__LINE__,lb2=0)
       !$OMP BARRIER
    end if
    rp2_buffer(:,this%current_thread)=alpha
    !$OMP BARRIER
    if(this%current_thread==this%root_thread) then
       do i=0,this%num_threads-1
          if(i/=this%root_thread) then
             rp2_buffer(:,this%root_thread) = rp2_buffer(:,this%root_thread) + rp2_buffer(:,i)
          end if
       end  do
       call mpi_allreduce(rp2_buffer(:,this%root_thread),alpha,size(alpha),mpi_omp_context_rp,mpi_sum,this%icontxt,istat)
       check ( istat == mpi_success )
       do i=0,this%num_threads-1
          rp2_buffer(:,i) = alpha
       end  do
    end if
    !$OMP BARRIER
    alpha = rp2_buffer(:,this%current_thread)

  end subroutine mpi_omp_context_sum_vector_rp

  !=============================================================================
  subroutine mpi_omp_context_max_scalar_rp (this,alpha)
    implicit none
    class(mpi_omp_context_t) , intent(in)    :: this
    real(rp)             , intent(inout) :: alpha
    integer  :: i, istat
    rp1_buffer(this%current_thread) = alpha
    !$OMP BARRIER
    if(this%current_thread==this%root_thread) then
       do i=0,this%num_threads-1
          rp1_buffer(this%root_thread) = max(rp1_buffer(this%root_thread),rp1_buffer(i))
       end  do
       call mpi_allreduce(rp1_buffer(this%root_thread),alpha,1,mpi_omp_context_rp,mpi_max,this%icontxt,istat)
       check ( istat == mpi_success )
       rp1_buffer=alpha
    end if
    !$OMP BARRIER
    alpha = rp1_buffer(this%current_thread)
  end subroutine mpi_omp_context_max_scalar_rp

  !=============================================================================
  subroutine mpi_omp_context_max_vector_rp(this,alpha)
    implicit none
    class(mpi_omp_context_t) , intent(in)    :: this
    real(rp)             , intent(inout) :: alpha(:) 
    integer  :: i, istat

    if(size(alpha)>size(rp2_buffer,1)) then
       !$OMP BARRIER
       if(this%current_thread==this%root_thread) &
            & call memrealloc(size(alpha),this%num_threads,rp2_buffer,__FILE__,__LINE__,lb2=0)
       !$OMP BARRIER
    end if
    rp2_buffer(:,this%current_thread)=alpha
    !$OMP BARRIER
    if(this%current_thread==this%root_thread) then
       do i=0,this%num_threads-1
          rp2_buffer(:,this%root_thread) = max(rp2_buffer(:,this%root_thread),rp2_buffer(:,i))
       end  do
       call mpi_allreduce(rp2_buffer(:,this%root_thread),alpha,size(alpha),mpi_omp_context_rp,mpi_max,this%icontxt,istat)
       check ( istat == mpi_success )
       do i=0,this%num_threads-1
          rp2_buffer(:,i) = alpha
       end  do
    end if
    !$OMP BARRIER
    alpha = rp2_buffer(:,this%current_thread)

  end subroutine mpi_omp_context_max_vector_rp

  !=============================================================================
  subroutine mpi_omp_context_min_scalar_rp (this,alpha)
    implicit none
    class(mpi_omp_context_t) , intent(in)    :: this
    real(rp)             , intent(inout) :: alpha
    integer  :: i, istat
    rp1_buffer(this%current_thread) = alpha
    !$OMP BARRIER
    if(this%current_thread==this%root_thread) then
       do i=0,this%num_threads-1
          rp1_buffer(this%root_thread) = min(rp1_buffer(this%root_thread),rp1_buffer(i))
       end  do
       call mpi_allreduce(rp1_buffer,alpha,1,mpi_omp_context_rp,mpi_min,this%icontxt,istat); check ( istat == mpi_success )
       rp1_buffer=alpha
    end if
    !$OMP BARRIER
    alpha = rp1_buffer(this%current_thread)
  end subroutine mpi_omp_context_min_scalar_rp
  
  !=============================================================================
  subroutine mpi_omp_context_max_scalar_ip (this,n)
    implicit none
    class(mpi_omp_context_t) , intent(in)    :: this
    integer(ip)              , intent(inout) :: n
    integer  :: i, istat
    integer(ip) :: dat
    ip1_buffer(this%current_thread) = n
    !$OMP BARRIER
    if(this%current_thread==this%root_thread) then
      do i=0,this%num_threads-1
          ip1_buffer(this%root_thread) = max(ip1_buffer(this%root_thread),ip1_buffer(i))
       end  do
       call mpi_allreduce(ip1_buffer,n,1,mpi_omp_context_ip,mpi_max,this%icontxt,istat); check ( istat == mpi_success )
       ip1_buffer=n
    end if
    !$OMP BARRIER
    n = ip1_buffer(this%current_thread)
  end subroutine mpi_omp_context_max_scalar_ip
  
  !=============================================================================
  subroutine mpi_omp_context_sum_scalar_igp (this,n)
    implicit none
    class(mpi_omp_context_t) , intent(in)    :: this
    integer(igp)             , intent(inout) :: n
    integer  :: i, istat
    integer(ip) :: dat
    igp1_buffer(this%current_thread) = n
    !$OMP BARRIER
    if(this%current_thread==this%root_thread) then
      do i=0,this%num_threads-1
          ip1_buffer(this%root_thread) = max(ip1_buffer(this%root_thread),igp1_buffer(i))
       end  do
       call mpi_allreduce(igp1_buffer,n,1,mpi_omp_context_igp,mpi_sum,this%icontxt,istat); check ( istat == mpi_success )
       ip1_buffer=n
    end if
    !$OMP BARRIER
    n = ip1_buffer(this%current_thread)
  end subroutine mpi_omp_context_sum_scalar_igp
  
  !=============================================================================
  subroutine mpi_omp_context_sum_vector_igp(this,n)
    implicit none
    class(mpi_omp_context_t) , intent(in)    :: this
    integer(igp)             , intent(inout) :: n(:)
    mcheck(.false., "mpi_omp_context_sum_vector_igp :: Implementation pending!")
  end subroutine mpi_omp_context_sum_vector_igp
 
  !=============================================================================
  subroutine mpi_omp_context_bcast_subcontext(this,subcontxt1,subcontxt2,condition)
    implicit none
    class(mpi_omp_context_t)       , intent(in)    :: this
    class(execution_context_t) , intent(in)    :: subcontxt1
    class(execution_context_t) , intent(in)    :: subcontxt2
    logical                    , intent(inout) :: condition
    integer :: recv_rank, send_rank, istat

    ! this%root_thread

    select type(subcontxt2)
    type is(mpi_omp_context_t)
       select type(subcontxt1)
       type is(mpi_omp_context_t)

          if(this%current_thread==this%root_thread) then
             send_rank = mpi_omp_context_root_rank
             if(subcontxt1%am_i_member()) then
                recv_rank = subcontxt1%get_num_tasks()/subcontxt1%max_num_threads
             else if(subcontxt2%am_i_member()) then
                recv_rank = (this%get_num_tasks() - subcontxt2%get_num_tasks())/this%max_num_threads
             end if

             if(this%current_rank==send_rank.and.recv_rank<this%num_ranks) then
                call mpi_send(condition, 1, mpi_omp_context_lg, recv_rank,  &
                     & mpi_omp_context_tag, this%icontxt, istat); check( istat == mpi_success )
             else if(this%current_rank==recv_rank) then
                call mpi_recv(condition, 1, mpi_omp_context_lg, send_rank,  &
                     & mpi_omp_context_tag, this%icontxt, mpi_status_ignore, istat); check( istat == mpi_success )
             end if

             if(subcontxt2%am_i_member()) then
                call mpi_bcast(condition,1,mpi_omp_context_lg,mpi_omp_context_root_rank,subcontxt2%icontxt,istat); check( istat == mpi_success )
             end if

             if(condition) then
                ip1_buffer = 1
             else
                ip1_buffer = 0
             end if

          end if

          !$OMP BARRIER
          if(ip1_buffer(this%current_thread)==1) then
             condition = .true.
          else
             condition = .false.
          end if

          class default
          check(.false.)
       end select
       class default
       check(.false.)
    end select

  end subroutine mpi_omp_context_bcast_subcontext

  !=============================================================================
  ! When packing   (gathering) ,    buffer <- alpha * x
  ! When unpacking (scattering),    x <- beta*x + buffer
  subroutine mpi_omp_context_neighbours_exchange_rp ( this, & 
       &                                          num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       &                                          num_snd, list_snd, snd_ptrs, pack_idx,   &
       &                                          alpha, beta, x, y)
    implicit none
    class(mpi_omp_context_t), intent(inout) :: this

    ! Control info to receive
    integer(ip)             , intent(in) :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(in) :: unpack_idx (rcv_ptrs(num_rcv+1)-1)

    ! Control info to send
    integer(ip)             , intent(in) :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in) :: pack_idx (snd_ptrs(num_snd+1)-1)

    ! Floating point data
    real(rp), intent(in)    :: alpha, beta
    real(rp), intent(in)    :: x(:)
    real(rp), intent(inout) :: y(:)

    ! Communication related locals 
    integer :: i, proc_to_comm, task_to_comm, thread_to_comm, msg_tag, sizmsg, istat
    integer :: p2pstat(mpi_status_size)

    ! Request handlers for non-blocking receives
    integer, allocatable :: rcvhd(:)

    ! Request handlers for non-blocking receives
    integer, allocatable :: sndhd(:)

    real(rp), allocatable :: sndbuf(:) 
    real(rp), allocatable :: rcvbuf(:)

    real(rp) :: t_start, t_stop
    t_start = mpi_wtime ()
    
    call memalloc (num_rcv, rcvhd, __FILE__,__LINE__)
    call memalloc (num_snd, sndhd, __FILE__,__LINE__)

    call memalloc ((snd_ptrs(num_snd+1)-snd_ptrs(1)), sndbuf, __FILE__,__LINE__)
    call memalloc ((rcv_ptrs(num_rcv+1)-rcv_ptrs(1)), rcvbuf, __FILE__,__LINE__)

    ! Pack send buffers
    call pack_rp ( snd_ptrs(num_snd+1)-snd_ptrs(1), pack_idx, alpha, x, sndbuf )

    ! First post all the non blocking receives   
    do i=1, num_rcv
       task_to_comm = list_rcv(i) - 1 
       proc_to_comm = task_to_comm/this%max_num_threads
       thread_to_comm = task_to_comm - proc_to_comm * this%max_num_threads
       msg_tag = this%current_thread + this%max_num_threads * thread_to_comm ! send_thread + num_threads*recv_thread

       ! Message size to be received
       sizmsg = rcv_ptrs(i+1)-rcv_ptrs(i)

       if ( (sizmsg > 0) .and. (list_rcv(i)-1 /= this%get_current_task()) ) then
          call mpi_irecv(  rcvbuf(rcv_ptrs(i)), sizmsg,        &
               &  mpi_omp_context_rp, proc_to_comm, &
               &  msg_tag, this%icontxt, rcvhd(i), istat)
          check ( istat == mpi_success )
       end if
    end do

    ! Secondly post all non-blocking sends
    do i=1, num_snd
       task_to_comm = list_snd(i) - 1
       proc_to_comm = task_to_comm/this%max_num_threads
       thread_to_comm = task_to_comm - proc_to_comm * this%max_num_threads
       msg_tag = this%max_num_threads * this%current_thread + thread_to_comm ! send_thread + num_threads*recv_thread

       ! Message size to be sent
       sizmsg = snd_ptrs(i+1)-snd_ptrs(i)

       if ( (sizmsg > 0) .and. (list_snd(i)-1 /= this%get_current_task()) ) then 
          call mpi_isend(sndbuf(snd_ptrs(i)), sizmsg, &
               & mpi_omp_context_rp, proc_to_comm,    &
               & msg_tag, this%icontxt, sndhd(i), istat)
          check ( istat == mpi_success )
       end if
    end do

    ! Wait on all non-blocking receives
    do i=1, num_rcv
       proc_to_comm = list_rcv(i) - 1

       ! Message size to be received
       sizmsg = rcv_ptrs(i+1)-rcv_ptrs(i)

       if ( (sizmsg > 0) .and. (list_rcv(i)-1 /= this%get_current_task()) ) then
          call mpi_wait(rcvhd(i), p2pstat, istat)

       else if ( list_rcv(i)-1 == this%get_current_task() ) then
          if ( sizmsg /= snd_ptrs(i+1)-snd_ptrs(i) ) then 
             write(0,*) 'Fatal error in single_exchange: mismatch on self sendf', & 
                  & sizmsg, snd_ptrs(i+1)-snd_ptrs(i) 
          end if

          rcvbuf( rcv_ptrs(i):rcv_ptrs(i)+sizmsg-1) = &
               sndbuf( snd_ptrs(i): snd_ptrs(i)+sizmsg-1 )
       end if
    end do

    ! Finally wait on all non-blocking sends
    do i=1, num_snd
       proc_to_comm = list_snd(i) - 1 

       ! Message size to be received
       sizmsg = snd_ptrs(i+1)-snd_ptrs(i)

       if ( (sizmsg > 0) .and. (list_snd(i)-1 /= this%get_current_task()) ) then
          call mpi_wait(sndhd(i), p2pstat, istat)
          check ( istat == mpi_success )
       end if
    end do

    ! Unpack recv buffers
    call unpack_rp (rcv_ptrs(num_rcv+1)-rcv_ptrs(1), unpack_idx, beta, rcvbuf, y )

    call memfree (rcvhd,__FILE__,__LINE__) 
    call memfree (sndhd,__FILE__,__LINE__)

    call memfree (sndbuf,__FILE__,__LINE__)
    call memfree (rcvbuf,__FILE__,__LINE__)

    t_stop = mpi_wtime ()
    this%time_exchange = this%time_exchange + t_stop - t_start
    
  end subroutine mpi_omp_context_neighbours_exchange_rp

  
  !=============================================================================
  subroutine mpi_omp_context_neighbours_exchange_wo_alpha_beta_rp ( this, & 
       &                                                            num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       &                                                            num_snd, list_snd, snd_ptrs, pack_idx,   &
       &                                                            x, y, chunk_size)
    implicit none
    class(mpi_omp_context_t), intent(in)    :: this
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
  end subroutine mpi_omp_context_neighbours_exchange_wo_alpha_beta_rp
  
  !=============================================================================
  subroutine mpi_omp_context_neighbours_exchange_wo_alpha_beta_rp_v ( this, & 
       &                                                              num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       &                                                              num_snd, list_snd, snd_ptrs, pack_idx,   &
       &                                                              x, y, ptr_chunk_size_snd, ptr_chunk_size_rcv )
    implicit none
    class(mpi_omp_context_t), intent(in)    :: this
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
    mcheck(.false.,'mpi_omp_context_neighbours_exchange_wo_alpha_beta_rp_v is not implemented')
  end subroutine mpi_omp_context_neighbours_exchange_wo_alpha_beta_rp_v
  
  !=============================================================================
  ! When packing   (gathering) ,    buffer <- alpha * x
  ! When unpacking (scattering),    x <- beta*x + buffer
  subroutine mpi_omp_context_neighbours_exchange_ip ( this, & 
       &                                          num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       &                                          num_snd, list_snd, snd_ptrs, pack_idx,   &
       &                                          x,y,chunk_size)
    implicit none
    class(mpi_omp_context_t), intent(in)    :: this
    ! Control info to receive
    integer(ip)             , intent(in)    :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(in)    :: unpack_idx (rcv_ptrs(num_rcv+1)-1)
    ! Control info to send
    integer(ip)             , intent(in)    :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in)    :: pack_idx (snd_ptrs(num_snd+1)-1)
    ! Raw data to be exchanged
    integer(ip)             , intent(in)    :: x(:)
    integer(ip)             , intent(inout) :: y(:)
    integer(ip)   , optional, intent(in)    :: chunk_size

    ! Communication related locals 
    integer :: i, proc_to_comm, task_to_comm, thread_to_comm, msg_tag, sizmsg, istat
    integer :: p2pstat(mpi_status_size)

    ! Request handlers for non-blocking receives
    integer, allocatable :: rcvhd(:)

    ! Request handlers for non-blocking receives
    integer, allocatable :: sndhd(:)

    integer(ip), allocatable :: sndbuf(:) 
    integer(ip), allocatable :: rcvbuf(:)

    integer(ip) :: chunk_size_


    if ( present(chunk_size) ) then
       chunk_size_ = chunk_size
    else
       chunk_size_ = 1
    end if

    call memalloc (num_rcv, rcvhd, __FILE__,__LINE__)
    call memalloc (num_snd, sndhd, __FILE__,__LINE__)

    call memalloc ((snd_ptrs(num_snd+1)-snd_ptrs(1))*chunk_size_, sndbuf, __FILE__,__LINE__)
    call memalloc ((rcv_ptrs(num_rcv+1)-rcv_ptrs(1))*chunk_size_, rcvbuf, __FILE__,__LINE__)

    ! Pack send buffers
    call pack_ip ( snd_ptrs(num_snd+1)-snd_ptrs(1), chunk_size_, pack_idx, x, sndbuf )

    ! First post all the non blocking receives   
    do i=1, num_rcv
       task_to_comm = list_rcv(i) - 1 
       proc_to_comm = task_to_comm/this%max_num_threads
       thread_to_comm = task_to_comm - proc_to_comm * this%max_num_threads
       msg_tag = this%current_thread + this%max_num_threads * thread_to_comm ! send_thread + num_threads*recv_thread

       ! Message size to be received
       sizmsg = (rcv_ptrs(i+1)-rcv_ptrs(i))*chunk_size_

       if ( (sizmsg > 0) .and. (list_rcv(i)-1 /= this%get_current_task()) ) then
          call mpi_irecv(  rcvbuf((rcv_ptrs(i)-1)*chunk_size_+1), sizmsg,        &
               &  mpi_omp_context_ip, proc_to_comm, &
               &  msg_tag, this%icontxt, rcvhd(i), istat)
          check ( istat == mpi_success )
       end if
    end do

    ! Secondly post all non-blocking sends
    do i=1, num_snd
       task_to_comm = list_snd(i) - 1
       proc_to_comm = task_to_comm/this%max_num_threads
       thread_to_comm = task_to_comm - proc_to_comm * this%max_num_threads
       msg_tag = this%max_num_threads * this%current_thread + thread_to_comm ! send_thread + num_threads*recv_thread

       ! Message size to be sent
       sizmsg = (snd_ptrs(i+1)-snd_ptrs(i))*chunk_size_

       if ( (sizmsg > 0) .and. (list_snd(i)-1 /= this%get_current_task()) ) then 
          call mpi_isend(sndbuf((snd_ptrs(i)-1)*chunk_size_+1), sizmsg, &
               & mpi_omp_context_ip, proc_to_comm,    &
               & msg_tag, this%icontxt, sndhd(i), istat)
          check ( istat == mpi_success )
       end if
    end do

    ! Wait on all non-blocking receives
    do i=1, num_rcv
       proc_to_comm = list_rcv(i) - 1

       ! Message size to be received
       sizmsg = (rcv_ptrs(i+1)-rcv_ptrs(i))*chunk_size_

       if ( (sizmsg > 0) .and. (list_rcv(i)-1 /= this%get_current_task()) ) then
          call mpi_wait(rcvhd(i), p2pstat, istat)
          check ( istat == mpi_success )
       else if ( list_rcv(i)-1 == this%get_current_task() ) then
          if ( sizmsg /= (snd_ptrs(i+1)-snd_ptrs(i))*chunk_size_ ) then 
             write(0,*) 'Fatal error in single_exchange: mismatch on self sendf', & 
                  & sizmsg, snd_ptrs(i+1)-snd_ptrs(i) 
             check(.false.)
          end if
          rcvbuf((rcv_ptrs(i)-1)*chunk_size_+1:(rcv_ptrs(i)-1)*chunk_size_+sizmsg) = &
               sndbuf( (snd_ptrs(i)-1)*chunk_size_+1:(snd_ptrs(i)-1)*chunk_size_+sizmsg )
       end if
    end do

    ! Finally wait on all non-blocking sends
    do i=1, num_snd
       proc_to_comm = list_snd(i) - 1

       ! Message size to be received
       sizmsg = (snd_ptrs(i+1)-snd_ptrs(i))*chunk_size_

       if ( (sizmsg > 0) .and. (list_snd(i)-1 /= this%get_current_task()) ) then
          call mpi_wait(sndhd(i), p2pstat, istat)
          check ( istat == mpi_success )
       end if
    end do

    ! Unpack recv buffers
    call unpack_ip (rcv_ptrs(num_rcv+1)-rcv_ptrs(1), chunk_size_, unpack_idx, rcvbuf, y )

    call memfree (rcvhd,__FILE__,__LINE__) 
    call memfree (sndhd,__FILE__,__LINE__)

    call memfree (sndbuf,__FILE__,__LINE__)
    call memfree (rcvbuf,__FILE__,__LINE__)
  end subroutine mpi_omp_context_neighbours_exchange_ip

  !=============================================================================
  subroutine mpi_omp_context_neighbours_exchange_igp ( this, & 
       &                                              num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       &                                              num_snd, list_snd, snd_ptrs, pack_idx,   &
       &                                              x, y, chunk_size, mask)
    implicit none
    class(mpi_omp_context_t), intent(in)    :: this
    ! Control info to receive
    integer(ip)             , intent(in)    :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(in)    :: unpack_idx (rcv_ptrs(num_rcv+1)-1)
    ! Control info to send
    integer(ip)             , intent(in)    :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in)    :: pack_idx (snd_ptrs(num_snd+1)-1)
    ! Raw data to be exchanged
    integer(igp)            , intent(in)    :: x(:)
    integer(igp)            , intent(inout) :: y(:)
    integer(ip)   , optional, intent(in)    :: chunk_size
    integer(igp)  , optional, intent(in)    :: mask

    ! Communication related locals 
    integer :: i, proc_to_comm, task_to_comm, thread_to_comm, msg_tag, sizmsg, istat
    integer :: p2pstat(mpi_status_size)

    ! Request handlers for non-blocking receives
    integer, allocatable :: rcvhd(:)

    ! Request handlers for non-blocking receives
    integer, allocatable :: sndhd(:)

    integer(igp), allocatable :: sndbuf(:) 
    integer(igp), allocatable :: rcvbuf(:)

    integer(ip) :: chunk_size_

    if ( present(chunk_size) ) then
       chunk_size_ = chunk_size
    else
       chunk_size_ = 1
    end if

    call memalloc (num_rcv, rcvhd, __FILE__,__LINE__)
    call memalloc (num_snd, sndhd, __FILE__,__LINE__)

    call memalloc ((snd_ptrs(num_snd+1)-snd_ptrs(1))*chunk_size_, sndbuf, __FILE__,__LINE__)
    call memalloc ((rcv_ptrs(num_rcv+1)-rcv_ptrs(1))*chunk_size_, rcvbuf, __FILE__,__LINE__)

    ! Pack send buffers
    call pack_igp ( snd_ptrs(num_snd+1)-snd_ptrs(1), chunk_size_, pack_idx, x, sndbuf )

    ! First post all the non blocking receives   
    do i=1, num_rcv
       task_to_comm = list_rcv(i) - 1
       proc_to_comm = task_to_comm/this%max_num_threads
       thread_to_comm = task_to_comm - proc_to_comm * this%max_num_threads
       msg_tag = this%current_thread + this%max_num_threads * thread_to_comm ! send_thread + num_threads*recv_thread

       ! Message size to be received
       sizmsg = (rcv_ptrs(i+1)-rcv_ptrs(i))*chunk_size_

       if ( (sizmsg > 0) .and. (list_rcv(i)-1 /= this%get_current_task()) ) then
          call mpi_irecv(  rcvbuf((rcv_ptrs(i)-1)*chunk_size_+1), sizmsg,        &
               &  mpi_omp_context_igp, proc_to_comm, &
               &  msg_tag, this%icontxt, rcvhd(i), istat)
          check ( istat == mpi_success )
       end if
    end do

    ! Secondly post all non-blocking sends
    do i=1, num_snd
       task_to_comm = list_snd(i) - 1
       proc_to_comm = task_to_comm/this%max_num_threads
       thread_to_comm = task_to_comm - proc_to_comm * this%max_num_threads
       msg_tag = this%max_num_threads * this%current_thread + thread_to_comm ! send_thread + num_threads*recv_thread

       ! Message size to be sent
       sizmsg = (snd_ptrs(i+1)-snd_ptrs(i))*chunk_size_

       if ( (sizmsg > 0) .and. (list_snd(i)-1 /= this%get_current_task()) ) then 
          call mpi_isend(sndbuf((snd_ptrs(i)-1)*chunk_size_+1), sizmsg, &
               & mpi_omp_context_igp, proc_to_comm,    &
               & msg_tag, this%icontxt, sndhd(i), istat)
          check ( istat == mpi_success )
       end if
    end do

    ! Wait on all non-blocking receives
    do i=1, num_rcv
       proc_to_comm = list_rcv(i) - 1

       ! Message size to be received
       sizmsg = (rcv_ptrs(i+1)-rcv_ptrs(i))*chunk_size_

       if ( (sizmsg > 0) .and. (list_rcv(i)-1 /= this%get_current_task()) ) then
          call mpi_wait(rcvhd(i), p2pstat, istat)
          check ( istat == mpi_success )
       else if ( list_rcv(i)-1 == this%get_current_task() ) then
          if ( sizmsg /= (snd_ptrs(i+1)-snd_ptrs(i))*chunk_size_ ) then 
             write(0,*) 'Fatal error in single_exchange: mismatch on self sendf', & 
                  & sizmsg, snd_ptrs(i+1)-snd_ptrs(i) 
             check(.false.)     
          end if
          rcvbuf( (rcv_ptrs(i)-1)*chunk_size_+1:(rcv_ptrs(i)-1)*chunk_size_+sizmsg) = &
               sndbuf( (snd_ptrs(i)-1)*chunk_size_+1:(snd_ptrs(i)-1)*chunk_size_+sizmsg )
       end if
    end do

    ! Finally wait on all non-blocking sends
    do i=1, num_snd
       proc_to_comm = list_snd(i) - 1

       ! Message size to be received
       sizmsg = (snd_ptrs(i+1)-snd_ptrs(i))*chunk_size_

       if ( (sizmsg > 0) .and. (list_snd(i)-1 /= this%get_current_task()) ) then
          call mpi_wait(sndhd(i), p2pstat, istat)
          check ( istat == mpi_success )
       end if
    end do

    ! Unpack recv buffers
    call unpack_igp (rcv_ptrs(num_rcv+1)-rcv_ptrs(1), chunk_size_, unpack_idx, rcvbuf, y, mask )

    call memfree (rcvhd,__FILE__,__LINE__) 
    call memfree (sndhd,__FILE__,__LINE__)

    call memfree (sndbuf,__FILE__,__LINE__)
    call memfree (rcvbuf,__FILE__,__LINE__)
  end subroutine mpi_omp_context_neighbours_exchange_igp

  !=============================================================================
  subroutine mpi_omp_context_neighbours_exchange_single_ip ( this, & 
       &                                                    num_neighbours, &
       &                                                    list_neighbours, &
       &                                                    input_data,&
       &                                                    output_data)
    implicit none
    class(mpi_omp_context_t), intent(in) :: this

    integer                 , intent(in)    :: num_neighbours
    integer(ip)             , intent(in)    :: list_neighbours (num_neighbours)
    integer(ip)             , intent(in)    :: input_data
    integer(ip)             , intent(inout) :: output_data(num_neighbours)

    integer(ip), allocatable :: ptrs(:)        ! How much data does the part send/recv to/from each neighbour?
    integer(ip), allocatable :: unpack_idx(:)  ! Where the data received from each neighbour is copied/added 
    ! on the local vectors of the part ?
    integer(ip), allocatable :: pack_idx(:)    ! Where is located the data to be sent to 
    ! each neighbour on the local vectors of the part ?

    integer(ip), allocatable :: buffer(:)  
    integer(ip)              :: i 

    call memalloc ( num_neighbours+1, ptrs, __FILE__, __LINE__ )
    ptrs(1)=1
    do i=2, num_neighbours+1
       ptrs(i)=ptrs(i-1)+1
    end do

    call memalloc ( ptrs(num_neighbours+1)-1, pack_idx, __FILE__, __LINE__ )
    pack_idx = 1 

    call memalloc ( ptrs(num_neighbours+1)-1, unpack_idx, __FILE__, __LINE__ )
    do i=1, ptrs(num_neighbours+1)-1
       unpack_idx(i) = i + 1
    end do

    call memalloc ( num_neighbours+1, buffer, __FILE__, __LINE__ )
    buffer(1) = input_data

    call this%neighbours_exchange ( num_neighbours,    &
         list_neighbours,   &
         ptrs,              &
         unpack_idx,        &  
         num_neighbours,    &
         list_neighbours,   &
         ptrs,              &
         pack_idx,          &
         buffer,            &
         buffer )

    output_data = buffer(2:)

    call memfree (buffer    , __FILE__, __LINE__ )
    call memfree (pack_idx  , __FILE__, __LINE__ )
    call memfree (unpack_idx, __FILE__, __LINE__ )
    call memfree (ptrs      , __FILE__, __LINE__ )
  end subroutine mpi_omp_context_neighbours_exchange_single_ip

  !=============================================================================
  subroutine mpi_omp_context_neighbours_exchange_wo_pack_unpack_ieep ( this, &
       &                                                              num_neighbours, &
       &                                                              neighbour_ids, &
       &                                                              snd_ptrs, &
       &                                                              snd_buf, & 
       &                                                              rcv_ptrs, &
       &                                                              rcv_buf )
    implicit none
    class(mpi_omp_context_t)  , intent(in)    :: this 
    integer(ip)           , intent(in)    :: num_neighbours
    integer(ip)           , intent(in)    :: neighbour_ids(num_neighbours)
    integer(ip)           , intent(in)    :: snd_ptrs(num_neighbours+1)
    integer(ieep)         , intent(in)    :: snd_buf(snd_ptrs(num_neighbours+1)-1)   
    integer(ip)           , intent(in)    :: rcv_ptrs(num_neighbours+1)
    integer(ieep)         , intent(out)   :: rcv_buf(rcv_ptrs(num_neighbours+1)-1)

    ! Communication related locals 
    integer :: i, proc_to_comm, task_to_comm, thread_to_comm, msg_tag, sizmsg, istat
    integer :: p2pstat(mpi_status_size)

    ! Request handlers for non-blocking receives
    integer, allocatable, dimension(:) :: rcvhd

    ! Request handlers for non-blocking receives
    integer, allocatable, dimension(:) :: sndhd

    call memalloc (num_neighbours, rcvhd, __FILE__,__LINE__)
    call memalloc (num_neighbours, sndhd, __FILE__,__LINE__)

    ! First post all the non blocking receives   
    do i=1, num_neighbours
       task_to_comm = neighbour_ids(i) - 1
       proc_to_comm = task_to_comm/this%max_num_threads
       thread_to_comm = task_to_comm - proc_to_comm * this%max_num_threads
       msg_tag = this%current_thread + this%max_num_threads * thread_to_comm ! send_thread + num_threads*recv_thread

       ! Message size to be received
       sizmsg = rcv_ptrs(i+1)-rcv_ptrs(i)

       if ( (sizmsg > 0) .and. (neighbour_ids(i)-1 /= this%get_current_task()) ) then
          call mpi_irecv(  rcv_buf(rcv_ptrs(i)), sizmsg, &
               &  mpi_omp_context_ieep, proc_to_comm, &
               &  msg_tag, this%icontxt, rcvhd(i), istat)
          check ( istat == mpi_success )
       end if
    end do

    ! Secondly post all non-blocking sends
    do i=1, num_neighbours
       task_to_comm = neighbour_ids(i) - 1
       proc_to_comm = task_to_comm/this%max_num_threads
       thread_to_comm = task_to_comm - proc_to_comm * this%max_num_threads
       msg_tag = this%max_num_threads * this%current_thread + thread_to_comm ! send_thread + num_threads*recv_thread

       ! Message size to be sent
       sizmsg = snd_ptrs(i+1)-snd_ptrs(i)

       if ( (sizmsg > 0) .and. (neighbour_ids(i)-1 /= this%get_current_task()) ) then 
          call mpi_isend(snd_buf(snd_ptrs(i)), sizmsg, &
               & mpi_omp_context_ieep, proc_to_comm, &
               & msg_tag, this%icontxt, sndhd(i), istat)
          check ( istat == mpi_success )
       end if
    end do

    ! Wait on all non-blocking receives
    do i=1, num_neighbours
       proc_to_comm = neighbour_ids(i) - 1

       ! Message size to be received
       sizmsg = rcv_ptrs(i+1)-rcv_ptrs(i)

       if ( (sizmsg > 0) .and. (neighbour_ids(i)-1 /= this%get_current_task()) ) then
          call mpi_wait(rcvhd(i), p2pstat, istat)
          check (istat == mpi_success)
       else if ( neighbour_ids(i)-1 == this%get_current_task() ) then
          if ( sizmsg /= snd_ptrs(i+1)-snd_ptrs(i) ) then 
             write(0,*) 'Fatal error in single_exchange: mismatch on self sendf', & 
                  & sizmsg, snd_ptrs(i+1)-snd_ptrs(i)
          end if
          rcv_buf( rcv_ptrs(i):rcv_ptrs(i+1)-1 ) = &
               snd_buf( snd_ptrs(i):snd_ptrs(i+1)-1 )
       end if
    end do

    ! Finally wait on all non-blocking sends
    do i=1, num_neighbours
       proc_to_comm = neighbour_ids(i) - 1

       ! Message size to be received
       sizmsg = snd_ptrs(i+1)-snd_ptrs(i)

       if ( (sizmsg > 0) .and. (neighbour_ids(i)-1 /= this%get_current_task()) ) then
          call mpi_wait(sndhd(i), p2pstat, istat)
          check ( istat == mpi_success )
       end if
    end do

    call memfree (rcvhd ,__FILE__,__LINE__) 
    call memfree (sndhd ,__FILE__,__LINE__)
  end subroutine mpi_omp_context_neighbours_exchange_wo_pack_unpack_ieep

  !=============================================================================
  subroutine mpi_omp_context_neighbours_exchange_wo_unpack_ip ( this, &
                                                               num_rcv, list_rcv, rcv_ptrs, rcv_buf, &
                                                               num_snd, list_snd, snd_ptrs, pack_idx,   &
                                                               x, chunk_size)
    implicit none
    class(mpi_omp_context_t) , intent(in)    :: this
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
  end subroutine mpi_omp_context_neighbours_exchange_wo_unpack_ip

  !=============================================================================
  subroutine mpi_omp_context_neighbours_exchange_variable_igp ( this, & 
       &                                              num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       &                                              num_snd, list_snd, snd_ptrs, pack_idx,   &
       &                                              x, y, ptr_chunk_size, mask)
    implicit none
    class(mpi_omp_context_t), intent(in)    :: this
    ! Control info to receive
    integer(ip)             , intent(in)    :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(in)    :: unpack_idx (rcv_ptrs(num_rcv+1)-1)
    ! Control info to send
    integer(ip)             , intent(in)    :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in)    :: pack_idx (snd_ptrs(num_snd+1)-1)
    ! Raw data to be exchanged
    integer(igp)            , intent(in)    :: x(:)
    integer(igp)            , intent(inout) :: y(:)
    integer(ip)             , intent(in)    :: ptr_chunk_size(:)
    integer(igp)  , optional, intent(in)    :: mask
    mcheck(.false.,'mpi_omp_context_neighbours_exchange_variable_igp is not implemented')
  end subroutine mpi_omp_context_neighbours_exchange_variable_igp
  
  !=============================================================================
  subroutine mpi_omp_context_neighbours_exchange_variable_ip ( this, & 
       &                                              num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       &                                              num_snd, list_snd, snd_ptrs, pack_idx,   &
       &                                              x, y, ptr_chunk_size, mask)
    implicit none
    class(mpi_omp_context_t), intent(in)    :: this
    ! Control info to receive
    integer(ip)             , intent(in)    :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(in)    :: unpack_idx (rcv_ptrs(num_rcv+1)-1)
    ! Control info to send
    integer(ip)             , intent(in)    :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in)    :: pack_idx (snd_ptrs(num_snd+1)-1)
    ! Raw data to be exchanged
    integer(ip)             , intent(in)    :: x(:)
    integer(ip)             , intent(inout) :: y(:)
    integer(ip)             , intent(in)    :: ptr_chunk_size(:)
    integer(ip)   , optional, intent(in)    :: mask
    mcheck(.false.,'mpi_omp_context_neighbours_exchange_variable_ip is not implemented')
  end subroutine mpi_omp_context_neighbours_exchange_variable_ip
  
  
  !=============================================================================
  subroutine mpi_omp_context_gather_scalar_ip ( this, input_data, output_data )
    implicit none
    class(mpi_omp_context_t), intent(in)   :: this
    integer(ip)         , intent(in)   :: input_data
    integer(ip)         , intent(out)  :: output_data(:) ! (this%get_num_tasks())
    integer  ::  istat
    
    ip1_buffer(this%current_thread)=input_data
    !$OMP BARRIER
    if(this%current_thread==this%root_thread) then
       call mpi_gather( ip1_buffer, this%max_num_threads, mpi_omp_context_ip, &
         &              output_data, this%max_num_threads, mpi_omp_context_ip, mpi_omp_context_root_rank, this%icontxt, istat)
       check( istat == mpi_success )
    end if
  end subroutine mpi_omp_context_gather_scalar_ip

  !=============================================================================
  subroutine mpi_omp_context_scatter_scalar_ip ( this, input_data, output_data )
    implicit none
    class(mpi_omp_context_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data(:) ! (this%get_num_tasks())
    integer(ip)             , intent(out)  :: output_data
    integer  ::  istat

    if(this%current_thread==this%root_thread) then
       call mpi_scatter( input_data, this%max_num_threads, mpi_omp_context_ip, &
                         ip1_buffer, this%max_num_threads, mpi_omp_context_ip, mpi_omp_context_root_rank, this%icontxt, istat)
       check( istat == mpi_success )
    end if
    !$OMP BARRIER
    output_data = ip1_buffer(this%current_thread)
  end subroutine mpi_omp_context_scatter_scalar_ip

  !=============================================================================
  subroutine mpi_omp_context_bcast_scalar_ip ( this, data )
    implicit none
    class(mpi_omp_context_t), intent(in)    :: this
    integer(ip)             , intent(inout) :: data
    integer  ::  istat
    if(this%current_thread==this%root_thread) then
       call mpi_bcast(data,1,mpi_omp_context_ip,mpi_omp_context_root_rank,this%icontxt,istat); check( istat == mpi_success )
       ip1_buffer = data
    end if
    !$OMP BARRIER
    data = ip1_buffer(this%current_thread)
  end subroutine mpi_omp_context_bcast_scalar_ip
  
  !=============================================================================
  subroutine mpi_omp_context_bcast_scalar_ip_1D_array ( this, data )
    implicit none
    class(mpi_omp_context_t), intent(in)    :: this
    integer(ip)             , intent(inout) :: data(:)
    integer  ::  istat
    if(this%current_thread==this%root_thread) then
       call mpi_bcast(data,size(data),mpi_omp_context_ip,mpi_omp_context_root_rank,this%icontxt,istat); check( istat == mpi_success )
       ip1_buffer = data
    end if
    !$OMP BARRIER
    data = ip1_buffer(this%current_thread)
  end subroutine mpi_omp_context_bcast_scalar_ip_1D_array

  !=============================================================================
  subroutine mpi_omp_context_gather_scalar_igp ( this, input_data, output_data )
    implicit none
    class(mpi_omp_context_t), intent(in)   :: this
    integer(igp)         , intent(in)   :: input_data
    integer(igp)         , intent(out)  :: output_data(:) ! (this%get_num_tasks())
    integer  ::  istat
    
    igp1_buffer(this%current_thread)=input_data
    !$OMP BARRIER
    if(this%current_thread==this%root_thread) then
       call mpi_gather( igp1_buffer, this%max_num_threads, mpi_omp_context_igp, &
         &              output_data, this%max_num_threads, mpi_omp_context_igp, mpi_omp_context_root_rank, this%icontxt, istat)
       check( istat == mpi_success )
    end if
  end subroutine mpi_omp_context_gather_scalar_igp

  !=============================================================================
  subroutine mpi_omp_context_scatter_scalar_igp ( this, input_data, output_data )
    implicit none
    class(mpi_omp_context_t), intent(in)   :: this
    integer(igp)             , intent(in)   :: input_data(:) ! (this%get_num_tasks())
    integer(igp)             , intent(out)  :: output_data
    integer  ::  istat

    if(this%current_thread==this%root_thread) then
       call mpi_scatter( input_data, this%max_num_threads, mpi_omp_context_igp, &
                         igp1_buffer, this%max_num_threads, mpi_omp_context_igp, mpi_omp_context_root_rank, this%icontxt, istat)
       check( istat == mpi_success )
    end if
    !$OMP BARRIER
    output_data = igp1_buffer(this%current_thread)
  end subroutine mpi_omp_context_scatter_scalar_igp

  !=============================================================================
  subroutine mpi_omp_context_bcast_scalar_igp ( this, data )
    implicit none
    class(mpi_omp_context_t), intent(in)    :: this
    integer(igp)             , intent(inout) :: data
    integer  ::  istat
    if(this%current_thread==this%root_thread) then
       call mpi_bcast(data,1,mpi_omp_context_igp,mpi_omp_context_root_rank,this%icontxt,istat); check( istat == mpi_success )
       igp1_buffer = data
    end if
    !$OMP BARRIER
    data = igp1_buffer(this%current_thread)
  end subroutine mpi_omp_context_bcast_scalar_igp

  !=============================================================================
  !=============================================================================
  subroutine pack_rp ( n, pack_idx, alpha, x, y )
    implicit none

    !Parameters
    integer (ip), intent(in)   :: n
    integer (ip), intent(in)   :: pack_idx(n)
    real    (rp), intent(in)   :: alpha
    real    (rp), intent(in)   :: x(*)
    real    (rp), intent(inout):: y(*)

    !Locals
    integer(ip) :: i

    if (alpha == 0.0_rp) then 
       !do nothing
    else if (alpha == 1.0_rp) then 
       do i=1,n
          y(i) = x(pack_idx(i))
       end do
    else if (alpha == -1.0_rp) then 
       do i=1,n
          y(i) = x(pack_idx(i))
       end do
    else  
       do i=1,n
          y(i) = alpha*x(pack_idx(i))
       end do
    end if

  end subroutine pack_rp

  !=============================================================================
  subroutine unpack_rp ( n, unpack_idx, beta, x, y )
    implicit none

    !Parameters
    integer(ip), intent(in)    :: n
    integer(ip), intent(in)    :: unpack_idx(n)
    real(rp)   , intent(in)    :: beta
    real(rp)   , intent(in)    :: x(*)
    real(rp)   , intent(inout) :: y(*)

    !Locals
    integer(ip) :: i

    if (beta == 0.0_rp) then
       do i=1,n
          y(unpack_idx(i)) = x(i)
       end do
    else if (beta == 1.0_rp) then
       do i=1,n
          y(unpack_idx(i)) = y(unpack_idx(i)) + x(i)
       end do
    else
       do i=1,n
          y(unpack_idx(i)) = beta*y(unpack_idx(i)) + x(i)
       end do
    end if
  end subroutine unpack_rp

  !=============================================================================
  subroutine pack_ip ( n, chunk_size, pack_idx, x, y )
    implicit none

    ! Parameters
    integer (ip), intent(in)     :: n
    integer (ip), intent(in)     :: chunk_size
    integer (ip), intent(in)     :: pack_idx(n)
    integer (ip), intent(in)    :: x(*)
    integer (ip), intent(inout) :: y(*)

    ! Locals
    integer(ip) :: i, j, startx, endx
    integer(ip) :: current
    current=1
    do i=1,n
       startx = (pack_idx(i)-1)*chunk_size + 1
       endx   = startx + chunk_size - 1
       do j=startx, endx
          y(current) = x(j)
          current = current + 1
       end do
    end do
  end subroutine pack_ip

  !=============================================================================
  subroutine unpack_ip ( n, chunk_size, unpack_idx, x, y )
    implicit none

    ! Parameters
    integer (ip), intent(in)     :: n
    integer (ip), intent(in)     :: chunk_size
    integer (ip), intent(in)     :: unpack_idx(n)
    integer (ip), intent(in)    :: x(*)
    integer (ip), intent(inout) :: y(*)

    ! Locals
    integer(ip) :: i, j, starty, endy, current
    current = 1
    do i=1,n
       starty = (unpack_idx(i)-1)*chunk_size + 1
       endy   = starty + chunk_size - 1
       do j=starty, endy
          y(j) = x(current)
          current = current + 1
       end do
    end do
  end subroutine unpack_ip

  !=============================================================================
  subroutine pack_igp ( n, chunk_size, pack_idx, x, y )
    implicit none

    !Parameters
    integer (ip), intent(in)     :: n
    integer (ip), intent(in)     :: chunk_size
    integer (ip), intent(in)     :: pack_idx(n)
    integer (igp), intent(in)    :: x(*)
    integer (igp), intent(inout) :: y(*)

    !Locals
    integer(ip) :: i, j, startx, endx
    integer(ip) :: current
    current=1
    do i=1,n
       startx = (pack_idx(i)-1)*chunk_size + 1
       endx   = startx + chunk_size - 1
       do j=startx, endx
          y(current) = x(j)
          current = current + 1
       end do
    end do
  end subroutine pack_igp

  !=============================================================================
  subroutine unpack_igp ( n, chunk_size, unpack_idx, x, y, mask )
    implicit none

    !Parameters
    integer (ip), intent(in)     :: n
    integer (ip), intent(in)     :: chunk_size
    integer (ip), intent(in)     :: unpack_idx(n)
    integer (igp), intent(in)    :: x(*)
    integer (igp), intent(inout) :: y(*)
    integer (igp), optional, intent(in) :: mask
    
    !Locals
    integer(ip) :: i, j, starty, endy, current
    current = 1
    do i=1,n
       starty = (unpack_idx(i)-1)*chunk_size + 1
       endy   = starty + chunk_size - 1
       do j=starty, endy
          if (present(mask)) then
            if ( x(current) /= mask ) then
              y(j) = x(current)
            end if
          else
            y(j) = x(current)
          end if
          current = current + 1
       end do
    end do
  end subroutine unpack_igp

  !=============================================================================
  subroutine mpi_omp_context_send_ip ( this, rcv_task, data )
    implicit none
    class(mpi_omp_context_t), intent(in)    :: this
    integer(ip)         , intent(in)    :: rcv_task
    integer(ip)         , intent(in)    :: data
    mcheck(.false.,'mpi_omp_context_send_ip not implemented')
  end subroutine mpi_omp_context_send_ip

  subroutine mpi_omp_context_rcv_ip ( this, send_task, data )
    implicit none
    class(mpi_omp_context_t) , intent(in)    :: this
    integer(ip)          , intent(inout) :: send_task
    integer(ip)          , intent(inout) :: data
    mcheck(.false.,'mpi_omp_context_rcv_ip_ not implemented')
  end subroutine mpi_omp_context_rcv_ip

  !=============================================================================
  subroutine mpi_omp_context_send_igp ( this, rcv_task, data )
    implicit none
    class(mpi_omp_context_t), intent(in)    :: this
    integer(ip)         , intent(in)    :: rcv_task
    integer(igp)         , intent(in)    :: data
    mcheck(.false.,'mpi_omp_context_send_igp not implemented')
  end subroutine mpi_omp_context_send_igp

  subroutine mpi_omp_context_rcv_igp ( this, send_task, data )
    implicit none
    class(mpi_omp_context_t) , intent(in)    :: this
    integer(ip)          , intent(inout) :: send_task
    integer(igp)          , intent(inout) :: data
    mcheck(.false.,'mpi_omp_context_rcv_igp_ not implemented')
  end subroutine mpi_omp_context_rcv_igp
  
  !=============================================================================
  subroutine mpi_omp_context_send_rp ( this, rcv_task, data )
    implicit none
    class(mpi_omp_context_t), intent(in)    :: this
    integer(ip)         , intent(in)    :: rcv_task
    real(rp)         , intent(in)    :: data
    mcheck(.false.,'mpi_omp_context_send_rp not implemented')
  end subroutine mpi_omp_context_send_rp

  subroutine mpi_omp_context_rcv_rp ( this, send_task, data )
    implicit none
    class(mpi_omp_context_t) , intent(in)    :: this
    integer(ip)          , intent(inout) :: send_task
    real(rp)          , intent(inout) :: data
    mcheck(.false.,'mpi_omp_context_rcv_rp_ not implemented')
  end subroutine mpi_omp_context_rcv_rp

  !=============================================================================

  subroutine mpi_omp_context_send_ip_1D_array ( this, rcv_task, data )
    implicit none
    class(mpi_omp_context_t), intent(in)    :: this
    integer(ip)         , intent(in)    :: rcv_task
    integer(ip)         , intent(in)    :: data(:)
    mcheck(.false.,'mpi_omp_context_send_ip_1D_array not implemented')
  end subroutine mpi_omp_context_send_ip_1D_array

  subroutine mpi_omp_context_rcv_ip_1D_array ( this, send_task, data)
    implicit none
    class(mpi_omp_context_t) , intent(in)    :: this           
    integer(ip)         , intent(inout) :: send_task
    integer(ip)         , intent(inout) :: data(:)
    mcheck(.false.,'mpi_omp_context_rcv_ip_1D_array not implemented')
  end subroutine mpi_omp_context_rcv_ip_1D_array    

  !=============================================================================

  subroutine mpi_omp_context_send_igp_1D_array ( this, rcv_task, data )
    implicit none
    class(mpi_omp_context_t), intent(in)    :: this
    integer(ip)          , intent(in)    :: rcv_task
    integer(igp)         , intent(in)    :: data(:)
    mcheck(.false.,'mpi_omp_context_send_igp_1D_array not implemented')
  end subroutine mpi_omp_context_send_igp_1D_array

  subroutine mpi_omp_context_rcv_igp_1D_array ( this, send_task, data)
    implicit none
    class(mpi_omp_context_t) , intent(in)    :: this
    integer(ip)           , intent(inout) :: send_task
    integer(igp)          , intent(inout) :: data(:)
    mcheck(.false.,'mpi_omp_context_rcv_igp_1D_array not implemented')
  end subroutine mpi_omp_context_rcv_igp_1D_array  
  
  !=============================================================================
  subroutine mpi_omp_context_send_rp_1D_array ( this, rcv_task, data )
    implicit none
    class(mpi_omp_context_t), intent(in)    :: this
    integer(ip)         , intent(in)    :: rcv_task
    real(rp)            , intent(in)    :: data(:)
    mcheck(.false.,'mpi_omp_context_send_rp_1D_array not implemented')
  end subroutine mpi_omp_context_send_rp_1D_array

  subroutine mpi_omp_context_rcv_rp_1D_array ( this, send_task, data)
    implicit none
    class(mpi_omp_context_t) , intent(in)    :: this
    integer(ip)          , intent(inout) :: send_task
    real(rp)          , intent(inout) :: data(:)
    mcheck(.false.,'mpi_omp_context_rcv_rp_1D_array not implemented')
  end subroutine mpi_omp_context_rcv_rp_1D_array

  !=============================================================================
  subroutine mpi_omp_context_root_send_master_rcv_ip ( this, input_data, output_data )
    implicit none
    class(mpi_omp_context_t), intent(in)      :: this
    integer(ip)             , intent(in)      :: input_data
    integer(ip)             , intent(inout)   :: output_data
    integer :: send_rank, recv_rank, istat
    integer :: send_thread
    send_rank = mpi_omp_context_root_rank
    recv_rank = (this%get_num_tasks()-1)/this%max_num_threads
    send_thread = this%root_thread
    if(this%current_rank==send_rank.and.this%current_thread==send_thread) then
       call mpi_send(input_data, 1, mpi_omp_context_ip, recv_rank,  &
               & mpi_omp_context_tag, this%icontxt, istat); check( istat == mpi_success )
    else if(this%current_rank==recv_rank) then
       call mpi_recv(output_data, 1, mpi_omp_context_ip, send_rank,  &
               & mpi_omp_context_tag, this%icontxt, mpi_status_ignore, istat); check( istat == mpi_success )
    end if
    end subroutine mpi_omp_context_root_send_master_rcv_ip

  !=============================================================================
  subroutine mpi_omp_context_root_send_master_rcv_ip_1D_array ( this, input_data, output_data )
    implicit none
    class(mpi_omp_context_t), intent(in)      :: this
    integer(ip)             , intent(in)      :: input_data(:)
    integer(ip)             , intent(inout)   :: output_data(:)
    integer :: send_rank, recv_rank, istat
    integer :: send_thread
    send_rank = mpi_omp_context_root_rank
    recv_rank = (this%get_num_tasks()-1)/this%max_num_threads
    send_thread = this%root_thread
    if(this%current_rank==send_rank.and.this%current_thread==send_thread) then
       call mpi_send(input_data, size(input_data), mpi_omp_context_ip, recv_rank,  &
               & mpi_omp_context_tag, this%icontxt, istat); check( istat == mpi_success )
    else if(this%current_rank==recv_rank) then
       call mpi_recv(output_data, size(output_data), mpi_omp_context_ip, send_rank,  &
               & mpi_omp_context_tag, this%icontxt, mpi_status_ignore, istat); check( istat == mpi_success )
    end if
  end subroutine mpi_omp_context_root_send_master_rcv_ip_1D_array

  !=============================================================================
  subroutine mpi_omp_context_root_send_master_rcv_rp ( this, input_data, output_data )
    implicit none
    class(mpi_omp_context_t), intent(in)      :: this
    real(rp)                , intent(in)      :: input_data
    real(rp)                , intent(inout)   :: output_data
    integer :: send_rank, recv_rank, istat
    integer :: send_thread
    send_rank = mpi_omp_context_root_rank
    recv_rank = (this%get_num_tasks()-1)/this%max_num_threads
    send_thread = this%root_thread
    if(this%current_rank==send_rank.and.this%current_thread==send_thread) then
       call mpi_send(input_data, 1, mpi_omp_context_rp, recv_rank,  &
               & mpi_omp_context_tag, this%icontxt, istat); check( istat == mpi_success )
    else if(this%current_rank==recv_rank) then
       call mpi_recv(output_data, 1, mpi_omp_context_rp, send_rank,  &
               & mpi_omp_context_tag, this%icontxt, mpi_status_ignore, istat); check( istat == mpi_success )
    end if
    end subroutine mpi_omp_context_root_send_master_rcv_rp

  !=============================================================================
  subroutine mpi_omp_context_root_send_master_rcv_rp_1D_array ( this, input_data, output_data )
    implicit none
    class(mpi_omp_context_t), intent(in)      :: this
    real(rp)                , intent(in)      :: input_data(:)
    real(rp)                , intent(inout)   :: output_data(:)
    integer :: send_rank, recv_rank, istat
    integer :: send_thread
    send_rank = mpi_omp_context_root_rank
    recv_rank = (this%get_num_tasks()-1)/this%max_num_threads
    send_thread = this%root_thread
    if(this%current_rank==send_rank.and.this%current_thread==send_thread) then
       call mpi_send(input_data, size(input_data), mpi_omp_context_rp, recv_rank,  &
               & mpi_omp_context_tag, this%icontxt, istat); check( istat == mpi_success )
    else if(this%current_rank==recv_rank) then
       call mpi_recv(output_data, size(output_data), mpi_omp_context_rp, send_rank,  &
               & mpi_omp_context_tag, this%icontxt, mpi_status_ignore, istat); check( istat == mpi_success )
    end if
  end subroutine mpi_omp_context_root_send_master_rcv_rp_1D_array
  
  !=============================================================================
  subroutine mpi_omp_context_root_send_master_rcv_logical ( this, input_data, output_data )
    implicit none
    class(mpi_omp_context_t), intent(in)      :: this
    logical                 , intent(in)      :: input_data
    logical                 , intent(inout)   :: output_data
    integer :: send_rank, recv_rank, istat
    integer :: send_thread
    send_rank = mpi_omp_context_root_rank
    recv_rank = (this%get_num_tasks()-1)/this%max_num_threads
    send_thread = this%root_thread
    if(this%current_rank==send_rank.and.this%current_thread==send_thread) then
       call mpi_send(input_data, 1, mpi_omp_context_lg, recv_rank,  &
               & mpi_omp_context_tag, this%icontxt, istat); check( istat == mpi_success )
    else if(this%current_rank==recv_rank) then
       call mpi_recv(output_data, 1, mpi_omp_context_lg, send_rank,  &
               & mpi_omp_context_tag, this%icontxt, mpi_status_ignore, istat); check( istat == mpi_success )
    end if
  end subroutine mpi_omp_context_root_send_master_rcv_logical
  
  !=============================================================================
  !=============================================================================
  subroutine mpi_omp_context_gather_to_master_ip ( this, input_data, output_data )
    implicit none
    class(mpi_omp_context_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data
    integer(ip)             , intent(out)  :: output_data(:) ! (this%get_num_tasks())
    integer  :: istat, master_rank
    ip1_buffer(this%current_thread)=input_data
    !$OMP BARRIER
    master_rank   = (this%get_num_tasks() - 1)/this%max_num_threads
    assert( (this%type==mpi_omp_context_homogeneous.and.this%current_rank/=master_rank).or.(this%current_rank==master_rank))
    if( (this%current_rank/=master_rank.and.this%current_thread==this%root_thread).or.this%current_rank==master_rank) then
       call mpi_gather( ip1_buffer , this%max_num_threads, mpi_omp_context_ip, &
            &           output_data, this%max_num_threads, mpi_omp_context_ip, &
            &           master_rank, this%icontxt, istat)
       !call mpi_gather( input_data, 1, mpi_omp_context_ip, output_data, 1, mpi_omp_context_ip, master, this%icontxt, istat)
       check( istat == mpi_success )
    end if
  end subroutine mpi_omp_context_gather_to_master_ip

  !=============================================================================
  subroutine mpi_omp_context_gather_to_master_igp ( this, input_data, output_data )
    implicit none
    class(mpi_omp_context_t), intent(in)   :: this
    integer(igp)            , intent(in)   :: input_data
    integer(igp)            , intent(out)  :: output_data(:) ! (this%get_num_tasks())
    integer  :: istat, master_rank
    igp1_buffer(this%current_thread)=input_data
    !$OMP BARRIER
    master_rank   = (this%get_num_tasks() - 1)/this%max_num_threads
    assert( (this%type==mpi_omp_context_homogeneous.and.this%current_rank/=master_rank).or.(this%current_rank==master_rank))
    if( (this%current_rank/=master_rank.and.this%current_thread==this%root_thread).or.this%current_rank==master_rank) then
       call mpi_gather( igp1_buffer, this%max_num_threads, mpi_omp_context_igp, &
            &           output_data, this%max_num_threads, mpi_omp_context_igp, &
            &           master_rank, this%icontxt, istat)
       !call mpi_gather( input_data, 1, mpi_omp_context_igp, output_data, 1, mpi_omp_context_igp, master, this%icontxt, istat)
       check( istat == mpi_success )
    end if

  end subroutine mpi_omp_context_gather_to_master_igp

  !=============================================================================
  subroutine mpi_omp_context_gather_to_master_ip_1D_array ( this, input_data_size, input_data, output_data )
    implicit none
    class(mpi_omp_context_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data_size
    integer(ip)             , intent(in)   :: input_data(input_data_size)
    integer(ip)             , intent(out)  :: output_data(:)
    integer  :: istat, master_rank, send_size
    master_rank   = (this%get_num_tasks()-1)/this%max_num_threads
    assert( (this%type==mpi_omp_context_homogeneous.and.this%current_rank/=master_rank).or.(this%current_rank==master_rank))
    send_size = input_data_size*this%max_num_threads
    if(this%current_rank/=master_rank) then
       if(this%current_thread==this%root_thread) then
          if(size(ip1_v_buffer)<send_size) &
               & call memrealloc(send_size,ip1_v_buffer,__FILE__,__LINE__,lb1=0)
       end if
       !$OMP BARRIER
       ip1_v_buffer(input_data_size*this%current_thread:input_data_size*(this%current_thread+1)-1)=input_data(:)
       !$OMP BARRIER
       if(this%current_thread==this%root_thread) then
          call mpi_gather( ip1_v_buffer, send_size, mpi_omp_context_ip, &
               &           output_data , send_size, mpi_omp_context_ip, &
               &           master_rank, this%icontxt, istat)
          check( istat == mpi_success )
       end if
    else if(this%current_rank==master_rank) then
       if(size(ip1_v_buffer)<send_size) &
            & call memrealloc(send_size,ip1_v_buffer,__FILE__,__LINE__,lb1=0)
       ip1_v_buffer=0
       call mpi_gather( ip1_v_buffer, send_size, mpi_omp_context_ip, &
            &           output_data , send_size, mpi_omp_context_ip, &
            &           master_rank, this%icontxt, istat)
       check( istat == mpi_success )
    end if

  end subroutine mpi_omp_context_gather_to_master_ip_1D_array

  !=============================================================================
  subroutine mpi_omp_context_gather_to_masterv_ip_1D_array ( this, input_data_size, input_data, recv_counts, displs, output_data )
    implicit none
    class(mpi_omp_context_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data_size
    integer(ip)             , intent(in)   :: input_data(input_data_size)
    integer(ip)             , intent(in)   :: recv_counts(:) ! (this%get_num_tasks())
    integer(ip)             , intent(in)   :: displs(:)      ! (this%get_num_tasks())
    integer(ip)             , intent(out)  :: output_data(:)
    integer :: istat, master_rank, send_size,i,j
    integer(ip), allocatable :: recv_counts_(:)
    integer(ip), allocatable :: displs_(:)

    master_rank   = (this%get_num_tasks()-1)/this%max_num_threads
    assert( (this%type==mpi_omp_context_homogeneous.and.this%current_rank/=master_rank).or.(this%current_rank==master_rank))

    if(this%current_rank/=master_rank) then
       ! Compute size of buffer
       ip1_buffer(this%current_thread+1) = input_data_size
       !$OMP BARRIER
       if(this%current_thread==this%root_thread) then
          ip1_buffer(0)=0
          do i=0,this%num_threads-1
             ip1_buffer(i+1)=ip1_buffer(i+1)+ip1_buffer(i)
          end do
          send_size = ip1_buffer(this%num_threads)
          if(size(ip1_v_buffer)<send_size) &
               & call memrealloc(send_size,ip1_v_buffer,__FILE__,__LINE__,lb1=0)
       end if
       !$OMP BARRIER
       ip1_v_buffer(ip1_buffer(this%current_thread):ip1_buffer(this%current_thread+1)-1)=input_data(:)
       !$OMP BARRIER
       if(this%current_thread==this%root_thread) then
          call mpi_gatherv(ip1_v_buffer, send_size, mpi_omp_context_ip, &
               &            output_data, recv_counts, displs, mpi_omp_context_ip, &
               &           master_rank, this%icontxt, istat)
          check( istat == mpi_success )
       end if
    else if(this%current_rank==master_rank) then
       call memalloc(this%num_ranks,recv_counts_,__FILE__,__LINE__)
       call memalloc(this%num_ranks,displs_,__FILE__,__LINE__)
       do i=1,this%num_ranks-1
          recv_counts_(i) = 0
          do j=1,this%max_num_threads
             recv_counts_(i) = recv_counts_(i) + recv_counts((i-1)*this%max_num_threads+j)
          end do
       end do
       recv_counts_(this%num_ranks) = recv_counts( (this%num_ranks-1)*this%max_num_threads+1)
       ! JP: I would like to eliminate displ from the interface and assume that data is always
       ! contiguous in the output vector. Here I assume this is the case and I don't use displ 
       displs_(1) = 0
       do i=2, this%num_ranks
         displs_(i) = displs_(i-1) + recv_counts_(i-1)
       end do
       call mpi_gatherv( input_data, input_data_size, mpi_omp_context_ip, &
            &           output_data, recv_counts_, displs_, mpi_omp_context_ip, &
            &           master_rank, this%icontxt, istat)
       check( istat == mpi_success )
       call memfree(recv_counts_,__FILE__,__LINE__)
       call memfree(displs_,__FILE__,__LINE__)
    end if

  end subroutine mpi_omp_context_gather_to_masterv_ip_1D_array

  subroutine mpi_omp_context_gather_to_masterv_igp_1D_array ( this, input_data_size, input_data, recv_counts, displs, output_data )
    implicit none
    class(mpi_omp_context_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data_size
    integer(igp)            , intent(in)   :: input_data(input_data_size)
    integer(ip)             , intent(in)   :: recv_counts(:) ! (this%get_num_tasks())
    integer(ip)             , intent(in)   :: displs(:)      ! (this%get_num_tasks())
    integer(igp)            , intent(out)  :: output_data(:)
    integer :: istat, master_rank, send_size,i,j
    integer(ip), allocatable :: recv_counts_(:)
    integer(ip), allocatable :: displs_(:)

    master_rank   = (this%get_num_tasks()-1)/this%max_num_threads
    assert( (this%type==mpi_omp_context_homogeneous.and.this%current_rank/=master_rank).or.(this%current_rank==master_rank))

    if(this%current_rank/=master_rank) then
       ! Compute size of buffer
       ip1_buffer(this%current_thread+1) = input_data_size
       !$OMP BARRIER
       if(this%current_thread==this%root_thread) then
          ip1_buffer(0)=0
          do i=0,this%num_threads-1
             ip1_buffer(i+1)=ip1_buffer(i+1)+ip1_buffer(i)
          end do
          send_size = ip1_buffer(this%num_threads)
          if(size(igp1_v_buffer)<send_size) &
               & call memrealloc(send_size,igp1_v_buffer,__FILE__,__LINE__,lb1=0)
       end if
       !$OMP BARRIER
       igp1_v_buffer(ip1_buffer(this%current_thread):ip1_buffer(this%current_thread+1)-1)=input_data(:)
       !$OMP BARRIER
       if(this%current_thread==this%root_thread) then
          call mpi_gatherv(igp1_v_buffer, send_size, mpi_omp_context_igp, &
               &             output_data, recv_counts, displs, mpi_omp_context_igp, &
               &           master_rank, this%icontxt, istat)
          check( istat == mpi_success )
       end if
    else if(this%current_rank==master_rank) then
       call memalloc(this%num_ranks,recv_counts_,__FILE__,__LINE__)
       call memalloc(this%num_ranks,displs_,__FILE__,__LINE__)
       do i=1,this%num_ranks-1
          recv_counts_(i) = 0
          do j=1,this%max_num_threads
             recv_counts_(i) = recv_counts_(i) + recv_counts((i-1)*this%max_num_threads+j)
          end do
       end do
       recv_counts_(this%num_ranks) = recv_counts( (this%num_ranks-1)*this%max_num_threads+1)
       ! JP: I would like to eliminate displ from the interface and assume that data is always
       ! contiguous in the output vector. Here I assume this is the case and I don't use displ 
       displs_(1) = 0
       do i=2, this%num_ranks
         displs_(i) = displs_(i-1) + recv_counts_(i-1)
       end do
       call mpi_gatherv( input_data, input_data_size, mpi_omp_context_igp, &
            &           output_data, recv_counts_, displs_, mpi_omp_context_igp, &
            &           master_rank, this%icontxt, istat)
       check( istat == mpi_success )
       call memfree(recv_counts_,__FILE__,__LINE__)
       call memfree(displs_,__FILE__,__LINE__)
    end if

  end subroutine mpi_omp_context_gather_to_masterv_igp_1D_array

  subroutine mpi_omp_context_gather_to_masterv_rp_1D_array ( this, input_data_size, input_data, recv_counts, displs, output_data )
    implicit none
    class(mpi_omp_context_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data_size
    real(rp)                , intent(in)   :: input_data(input_data_size)
    integer(ip)             , intent(in)   :: recv_counts(:) ! (this%get_num_tasks())
    integer(ip)             , intent(in)   :: displs(:)      ! (this%get_num_tasks())
    real(rp)                , intent(out)  :: output_data(:)
    integer :: istat, master_rank, send_size,i,j
    integer(ip), allocatable :: recv_counts_(:)
    integer(ip), allocatable :: displs_(:)

    master_rank   = (this%get_num_tasks()-1)/this%max_num_threads
    assert( (this%type==mpi_omp_context_homogeneous.and.this%current_rank/=master_rank).or.(this%current_rank==master_rank))

    if(this%current_rank/=master_rank) then
       ! Compute size of buffer
       ip1_buffer(this%current_thread+1) = input_data_size
       !$OMP BARRIER
       if(this%current_thread==this%root_thread) then
          ip1_buffer(0)=0
          do i=0,this%num_threads-1
             ip1_buffer(i+1)=ip1_buffer(i+1)+ip1_buffer(i)
          end do
          send_size = ip1_buffer(this%num_threads)
          if(size(rp1_v_buffer)<send_size) &
               & call memrealloc(send_size,rp1_v_buffer,__FILE__,__LINE__,lb1=0)
       end if
       !$OMP BARRIER
       rp1_v_buffer(ip1_buffer(this%current_thread):ip1_buffer(this%current_thread+1)-1) = input_data(:)
       !$OMP BARRIER
       if(this%current_thread==this%root_thread) then
          call mpi_gatherv(rp1_v_buffer, send_size, mpi_omp_context_rp, &
               &           output_data, recv_counts, displs, mpi_omp_context_rp, &
               &           master_rank, this%icontxt, istat)
          check( istat == mpi_success )
       end if
    else if(this%current_rank==master_rank) then
       call memalloc(this%num_ranks,recv_counts_,__FILE__,__LINE__)
       call memalloc(this%num_ranks,displs_,__FILE__,__LINE__)
       do i=1,this%num_ranks-1
          recv_counts_(i) = 0
          do j=1,this%max_num_threads
             recv_counts_(i) = recv_counts_(i) + recv_counts((i-1)*this%max_num_threads+j)
          end do
       end do
       recv_counts_(this%num_ranks) = recv_counts( (this%num_ranks-1)*this%max_num_threads+1)
       ! JP: I would like to eliminate displ from the interface and assume that data is always
       ! contiguous in the output vector. Here I assume this is the case and I don't use displ 
       displs_(1) = 0
       do i=2, this%num_ranks
         displs_(i) = displs_(i-1) + recv_counts_(i-1)
       end do
       call mpi_gatherv( input_data, input_data_size, mpi_omp_context_rp, &
            &           output_data, recv_counts_, displs_, mpi_omp_context_rp, &
            &           master_rank, this%icontxt, istat)
       check( istat == mpi_success )
       call memfree(recv_counts_,__FILE__,__LINE__)
       call memfree(displs_,__FILE__,__LINE__)
    end if

  end subroutine mpi_omp_context_gather_to_masterv_rp_1D_array

  subroutine mpi_omp_context_gather_to_masterv_rp_2D_array ( this, input_data, recv_counts, displs, output_data )
    implicit none
    class(mpi_omp_context_t), intent(in)   :: this
    real(rp)                , intent(in)   :: input_data(:,:)
    integer(ip)             , intent(in)   :: recv_counts(:) ! (this%get_num_tasks())
    integer(ip)             , intent(in)   :: displs(:)      ! (this%get_num_tasks())
    real(rp)                , intent(out)  :: output_data(:)
    integer :: istat, master_rank, send_size,i,j
    integer(ip), allocatable :: recv_counts_(:)
    integer(ip), allocatable :: displs_(:)

    master_rank   = (this%get_num_tasks()-1)/this%max_num_threads
    assert( (this%type==mpi_omp_context_homogeneous.and.this%current_rank/=master_rank).or.(this%current_rank==master_rank))

    if(this%current_rank/=master_rank) then
       ! Compute size of buffer
       ip1_buffer(this%current_thread+1) = size(input_data,1) * size(input_data,2)
       !$OMP BARRIER
       if(this%current_thread==this%root_thread) then
          ip1_buffer(0)=0
          do i=0,this%num_threads-1
             ip1_buffer(i+1)=ip1_buffer(i+1)+ip1_buffer(i)
          end do
          send_size = ip1_buffer(this%num_threads)
          if(size(rp1_v_buffer)<send_size) &
               & call memrealloc(send_size,rp1_v_buffer,__FILE__,__LINE__,lb1=0)
       end if
       !$OMP BARRIER
       do j=1,size(input_data,2)
          do i=1,size(input_data,1)
             rp1_v_buffer(ip1_buffer(this%current_thread)+(j-1)*size(input_data,1)+i-1) = input_data(i,j)
          end do
       end do
       !$OMP BARRIER
       if(this%current_thread==this%root_thread) then
          call mpi_gatherv(rp1_v_buffer, send_size, mpi_omp_context_rp, &
               &           output_data, recv_counts, displs, mpi_omp_context_rp, &
               &           master_rank, this%icontxt, istat)
          check( istat == mpi_success )
       end if
    else if(this%current_rank==master_rank) then
       call memalloc(this%num_ranks,recv_counts_,__FILE__,__LINE__)
       call memalloc(this%num_ranks,displs_,__FILE__,__LINE__)
       do i=1,this%num_ranks-1
          recv_counts_(i) = 0
          do j=1,this%max_num_threads
             recv_counts_(i) = recv_counts_(i) + recv_counts((i-1)*this%max_num_threads+j)
          end do
       end do
       recv_counts_(this%num_ranks) = recv_counts( (this%num_ranks-1)*this%max_num_threads+1)
       ! JP: I would like to eliminate displ from the interface and assume that data is always
       ! contiguous in the output vector. Here I assume this is the case and I don't use displ 
       displs_(1) = 0
       do i=2, this%num_ranks
         displs_(i) = displs_(i-1) + recv_counts_(i-1)
       end do
       call mpi_gatherv( input_data, size(input_data,1)*size(input_data,2), mpi_omp_context_rp, &
            &           output_data, recv_counts_, displs_, mpi_omp_context_rp, &
            &           master_rank, this%icontxt, istat)
       check( istat == mpi_success )
       call memfree(recv_counts_,__FILE__,__LINE__)
       call memfree(displs_,__FILE__,__LINE__)
    end if

  end subroutine mpi_omp_context_gather_to_masterv_rp_2D_array

  !=============================================================================
  ! subroutine mpi_omp_context_gather_to_masterv_igp_1D_array ( this, input_data_size, input_data, recv_counts, displs, output_data )
  !   implicit none
  !   class(mpi_omp_context_t), intent(in)   :: this
  !   integer(ip)         , intent(in)   :: input_data_size
  !   integer(igp)        , intent(in)   :: input_data(input_data_size)
  !   integer(ip)         , intent(in)   :: recv_counts(:) ! (this%get_num_tasks())
  !   integer(ip)         , intent(in)   :: displs(:)      ! (this%get_num_tasks())
  !   integer(igp)        , intent(out)  :: output_data(:)
  !   integer :: istat, master
  !   master    = this%get_num_tasks() - 1 
  !   call mpi_gatherv( input_data, input_data_size, mpi_omp_context_igp, &
  !        & output_data, recv_counts, displs, mpi_omp_context_igp, master, this%icontxt, istat)
  !   check( istat == mpi_success )
  ! end subroutine mpi_omp_context_gather_to_masterv_igp_1D_array

  !=============================================================================
  ! subroutine mpi_omp_context_gather_to_masterv_rp_1D_array ( this, input_data_size, input_data, recv_counts, displs, output_data )
  !   implicit none
  !   class(mpi_omp_context_t), intent(in)   :: this
  !   integer(ip)         , intent(in)   :: input_data_size
  !   real(rp)            , intent(in)   :: input_data(input_data_size)
  !   integer(ip)         , intent(in)   :: recv_counts(:) ! (this%get_num_tasks())
  !   integer(ip)         , intent(in)   :: displs(:)      ! (this%get_num_tasks())
  !   real(rp)            , intent(out)  :: output_data(:)
  !   integer :: istat, master
  !   master = this%get_num_tasks() - 1 
  !   call mpi_gatherv( input_data , input_data_size, mpi_omp_context_rp, &
  !        & output_data, recv_counts, displs, mpi_omp_context_rp, master, this%icontxt, istat)
  !   check( istat == mpi_success )
  ! end subroutine mpi_omp_context_gather_to_masterv_rp_1D_array

  !=============================================================================
  ! subroutine mpi_omp_context_gather_to_masterv_rp_2D_array ( this, input_data, recv_counts, displs, output_data )
  !   implicit none
  !   class(mpi_omp_context_t), intent(in)   :: this
  !   real(rp)            , intent(in)   :: input_data(:,:)
  !   integer(ip)         , intent(in)   :: recv_counts(:) ! (this%get_num_tasks())
  !   integer(ip)         , intent(in)   :: displs(:)      ! (this%get_num_tasks())
  !   real(rp)            , intent(out)  :: output_data(:)
  !   integer :: istat, master
  !   master = this%get_num_tasks() - 1 
  !   call mpi_gatherv( input_data , size(input_data,1)*size(input_data,2), mpi_omp_context_rp, &
  !        & output_data, recv_counts, displs, mpi_omp_context_rp, master, this%icontxt, istat)
  !   check( istat == mpi_success )
  ! end subroutine mpi_omp_context_gather_to_masterv_rp_2D_array
  
 !=============================================================================
  subroutine mpi_omp_context_scatter_from_master_ip ( this, input_data, output_data )
    implicit none
    class(mpi_omp_context_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data(:) ! (this%get_num_tasks())
    integer(ip)             , intent(out)  :: output_data
    integer  :: istat, master_rank
    mcheck( .false., 'mpi_omp_context_t::scatter_from_master_ip: This TPB has not been tested' )
    !$OMP BARRIER
    master_rank   = (this%get_num_tasks() - 1)/this%max_num_threads
    assert( (this%type==mpi_omp_context_homogeneous.and.this%current_rank/=master_rank).or.(this%current_rank==master_rank))
    if( (this%current_rank/=master_rank.and.this%current_thread==this%root_thread).or.this%current_rank==master_rank) then
       call mpi_scatter( input_data , this%max_num_threads, mpi_omp_context_ip, &
            &            output_data, this%max_num_threads, mpi_omp_context_ip, &
            &            master_rank, this%icontxt, istat)
       check( istat == mpi_success )
    end if
  end subroutine mpi_omp_context_scatter_from_master_ip
  
  !=============================================================================
  subroutine mpi_omp_context_scatter_from_masterv_ip_1D_array ( this, input_data, send_counts, displs, output_data_size, output_data )
    implicit none
    class(mpi_omp_context_t), intent(in)   :: this
    integer(ip)         , intent(in)   :: input_data(:)
    integer(ip)         , intent(in)   :: send_counts(:) ! (this%get_num_tasks())
    integer(ip)         , intent(in)   :: displs(:)      ! (this%get_num_tasks())
    integer(ip)         , intent(in)   :: output_data_size
    integer(ip)         , intent(out)  :: output_data(output_data_size)
    integer :: istat, master_rank, recv_size,i,j
    integer(ip), allocatable :: send_counts_(:)
    integer(ip), allocatable :: displs_(:)
    
    mcheck( .false., 'mpi_omp_context_t::scatter_from_masterv_ip_1D_array: This TPB has not been tested' )
    master_rank   = (this%get_num_tasks()-1)/this%max_num_threads
    assert( (this%type==mpi_omp_context_homogeneous.and.this%current_rank/=master_rank).or.(this%current_rank==master_rank))

    call memalloc(this%num_ranks,send_counts_,__FILE__,__LINE__)
    call memalloc(this%num_ranks,displs_,__FILE__,__LINE__)
    do i=1,this%num_ranks-1
       send_counts_(i) = 0
       do j=1,this%max_num_threads
          send_counts_(i) = send_counts_(i) + send_counts((i-1)*this%max_num_threads+j)
       end do
    end do
    send_counts_(this%num_ranks) = send_counts((this%num_ranks-1)*this%max_num_threads+1)
    ! JP: I would like to eliminate displ from the interface and assume that data is always
    ! contiguous in the output vector. Here I assume this is the case and I don't use displ 
    displs_(1) = 0
    do i=2, this%num_ranks
      displs_(i) = displs_(i-1) + send_counts_(i-1)
    end do
    if(this%current_rank==master_rank) then
       call mpi_scatterv( input_data, send_counts_, displs_, mpi_omp_context_ip, &
         &                output_data, output_data_size, mpi_omp_context_ip, &
         &                master_rank, this%icontxt, istat)
       check( istat == mpi_success )
    else if(this%current_rank/=master_rank) then
       ! Compute size of receiving buffer
       ip1_buffer(this%current_thread+1) = output_data_size
       !$OMP BARRIER
       if(this%current_thread==this%root_thread) then
          ip1_buffer(0)=0
          do i=0,this%num_threads-1
             ip1_buffer(i+1)=ip1_buffer(i+1)+ip1_buffer(i)
          end do
          recv_size = ip1_buffer(this%num_threads)
          if(size(ip1_v_buffer)<recv_size) &
               & call memrealloc(recv_size,ip1_v_buffer,__FILE__,__LINE__,lb1=0)
          call mpi_scatterv( input_data, send_counts_, displs_, mpi_omp_context_ip, &
               &             ip1_v_buffer, recv_size, mpi_omp_context_ip, master_rank, this%icontxt, istat)
          check( istat == mpi_success )
       end if
       !$OMP BARRIER
       output_data = ip1_v_buffer(ip1_buffer(this%current_thread):ip1_buffer(this%current_thread+1)-1)
    end if
    call memfree(send_counts_,__FILE__,__LINE__)
    call memfree(displs_,__FILE__,__LINE__)

  end subroutine mpi_omp_context_scatter_from_masterv_ip_1D_array
  
  
  !=============================================================================
  subroutine mpi_omp_context_scatter_from_masterv_rp_1D_array ( this, input_data, send_counts, displs, output_data_size, output_data )
    implicit none
    class(mpi_omp_context_t), intent(in)   :: this
    real(rp)            , intent(in)   :: input_data(:)
    integer(ip)         , intent(in)   :: send_counts(:) ! (this%get_num_tasks())
    integer(ip)         , intent(in)   :: displs(:)      ! (this%get_num_tasks())
    integer(ip)         , intent(in)   :: output_data_size
    real(rp)            , intent(out)  :: output_data(output_data_size)
    integer :: istat, master_rank, recv_size,i,j
    integer(ip), allocatable :: send_counts_(:)
    integer(ip), allocatable :: displs_(:)

    master_rank   = (this%get_num_tasks()-1)/this%max_num_threads
    assert( (this%type==mpi_omp_context_homogeneous.and.this%current_rank/=master_rank).or.(this%current_rank==master_rank))

    call memalloc(this%num_ranks,send_counts_,__FILE__,__LINE__)
    call memalloc(this%num_ranks,displs_,__FILE__,__LINE__)
    do i=1,this%num_ranks-1
       send_counts_(i) = 0
       do j=1,this%max_num_threads
          send_counts_(i) = send_counts_(i) + send_counts((i-1)*this%max_num_threads+j)
       end do
    end do
    send_counts_(this%num_ranks) = send_counts((this%num_ranks-1)*this%max_num_threads+1)
    ! JP: I would like to eliminate displ from the interface and assume that data is always
    ! contiguous in the output vector. Here I assume this is the case and I don't use displ 
    displs_(1) = 0
    do i=2, this%num_ranks
      displs_(i) = displs_(i-1) + send_counts_(i-1)
    end do
    if(this%current_rank==master_rank) then
       call mpi_scatterv( input_data, send_counts_, displs_, mpi_omp_context_rp, &
         &                output_data, output_data_size, mpi_omp_context_rp, &
         &                master_rank, this%icontxt, istat)
       check( istat == mpi_success )
    else if(this%current_rank/=master_rank) then
       ! Compute size of receiving buffer
       ip1_buffer(this%current_thread+1) = output_data_size
       !$OMP BARRIER
       if(this%current_thread==this%root_thread) then
          ip1_buffer(0)=0
          do i=0,this%num_threads-1
             ip1_buffer(i+1)=ip1_buffer(i+1)+ip1_buffer(i)
          end do
          recv_size = ip1_buffer(this%num_threads)
          if(size(rp1_v_buffer)<recv_size) &
               & call memrealloc(recv_size,rp1_v_buffer,__FILE__,__LINE__,lb1=0)
          call mpi_scatterv( input_data, send_counts_, displs_, mpi_omp_context_rp, &
               &             rp1_v_buffer, recv_size, mpi_omp_context_rp, master_rank, this%icontxt, istat)
          check( istat == mpi_success )
       end if
       !$OMP BARRIER
       output_data = rp1_v_buffer(ip1_buffer(this%current_thread):ip1_buffer(this%current_thread+1)-1)
    end if
    call memfree(send_counts_,__FILE__,__LINE__)
    call memfree(displs_,__FILE__,__LINE__)

  end subroutine mpi_omp_context_scatter_from_masterv_rp_1D_array

end module mpi_omp_context_names
