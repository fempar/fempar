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
module mesh_partitioner_names
  use types_names
  use memor_names
  use stdio_names
  use metis_names
  use mesh_names
  use mesh_distribution_names
  use mesh_partitioner_parameters_names
  use FPL
  implicit none
# include "debug.i90"
  private
  
  type mesh_partitioner_t
     private
     integer(ip)              :: nparts
     integer(ip)              :: num_levels
     integer(ip), allocatable :: num_parts_x_level (:)
     
     integer(ip)              :: debug                   ! Print (on screen?) info partition
     integer(ip)              :: strat                   ! Partitioning algorithm (part_kway,part_recursive,part_strip,part_rcm_strip)
     integer(ip)              :: metis_option_ufactor    ! Imbalance tol of metis_option_ufactor/1000 + 1
     integer(ip)              :: metis_option_minconn 
     integer(ip)              :: metis_option_contig  
     integer(ip)              :: metis_option_ctype 
     integer(ip)              :: metis_option_iptype
     integer(ip)              :: metis_option_debug
  contains 
     procedure, non_overridable          :: partition_mesh               => mesh_partitioner_partition_mesh
     procedure, non_overridable          :: set_default_parameter_values => mesh_partitioner_set_default_parameter_values
     procedure, non_overridable, private :: set_parameters_from_pl       => mesh_partitioner_set_parameters_from_pl
     procedure, non_overridable, private :: mesh_partitioner_write_mesh_parts_dir_path_prefix
     procedure, non_overridable, private :: mesh_partitioner_write_mesh_parts_pl
     generic                             :: write_mesh_parts             => mesh_partitioner_write_mesh_parts_dir_path_prefix, &
                                                                            mesh_partitioner_write_mesh_parts_pl
     procedure, non_overridable          :: free                         => mesh_partitioner_free
  end type mesh_partitioner_t
  
contains

  !=============================================================================
  subroutine mesh_partitioner_partition_mesh ( this, mesh, parameter_list ) 
     implicit none
     class(mesh_partitioner_t), intent(inout) :: this
     type(mesh_t)             , intent(in)    :: mesh
     type(ParameterList_t)   , intent(in)    :: parameter_list
     call this%free()
     call this%set_parameters_from_pl(parameter_list)
  end subroutine mesh_partitioner_partition_mesh
  
  !=============================================================================
  subroutine mesh_partitioner_set_parameters_from_pl(this,parameter_list)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    type(ParameterList_t)    , intent(in)    :: parameter_list
    ! Locals
    integer(ip)              :: istat
    integer(ip), allocatable :: param_size(:), param(:)

    ! Mandatory parameters: either nparts or num_levels
    assert(parameter_list%isPresent(key = num_parts_key).or.parameter_list%isPresent(key = num_levels_distribution_key))
    if( parameter_list%isPresent(num_parts_key)) then
       assert(parameter_list%isAssignable(num_parts_key, this%nparts))
       istat = parameter_list%get(key = num_parts_key , value = this%nparts)
       assert(istat==0)
    end if
    if( parameter_list%isPresent(num_levels_distribution_key) ) then
       assert(parameter_list%isAssignable(num_levels_distribution_key, this%num_levels))
       istat = parameter_list%get(key = num_levels_distribution_key  , value = this%num_levels)
       assert(istat==0)
       
       assert(parameter_list%isPresent(key = num_parts_x_level_key ))
       assert( parameter_list%GetDimensions(key = num_parts_x_level_key) == 1)

       ! Get the array using the local variable
       istat =  parameter_list%GetShape(key = num_parts_x_level_key, shape = param_size ); check(istat==0)
       call memalloc(param_size(1), param,__FILE__,__LINE__)
       assert(parameter_list%isAssignable(num_parts_x_level_key, param))
       istat = parameter_list%get(key = num_parts_x_level_key, value = param)
       assert(istat==0)

       call memalloc(this%num_levels, this%num_parts_x_level,__FILE__,__LINE__)
       this%num_parts_x_level = param(1:this%num_levels)
       call memfree(param,__FILE__,__LINE__)

       this%nparts = this%num_parts_x_level(1)
    else
       this%num_levels=1
       call memalloc(this%num_levels, this%num_parts_x_level,__FILE__,__LINE__)
       this%num_parts_x_level(1)=this%nparts
    end if

    ! Optional paramters
    if( parameter_list%isPresent(debug_key) ) then
       assert(parameter_list%isAssignable(debug_key, this%debug))
       istat = parameter_list%get(key = debug_key  , value = this%debug)
       assert(istat==0)
    end if

    if( parameter_list%isPresent(strategy_key) ) then
       assert(parameter_list%isAssignable(strategy_key, this%strat))
       istat = parameter_list%get(key = strategy_key  , value = this%strat)
       assert(istat==0)
       assert(this%strat==part_kway.or.this%strat==part_recursive.or.this%strat==part_strip.or.this%strat==part_rcm_strip)
    end if

    if( parameter_list%isPresent(metis_option_debug_key) ) then
       assert(parameter_list%isAssignable(metis_option_debug_key, this%metis_option_debug))
       istat = parameter_list%get(key = metis_option_debug_key  , value = this%metis_option_debug)
       check(istat==0)
    end if

    if( parameter_list%isPresent(metis_option_ufactor_key) ) then
       assert(parameter_list%isAssignable(metis_option_ufactor_key, this%metis_option_ufactor))
       istat = parameter_list%get(key = metis_option_ufactor_key, value = this%metis_option_ufactor)
       assert(istat==0)
    end if

    if( parameter_list%isPresent(metis_option_minconn_key) ) then
       assert(parameter_list%isAssignable(metis_option_minconn_key, this%metis_option_minconn))
       istat = parameter_list%get(key = metis_option_minconn_key, value = this%metis_option_minconn)
       check(istat==0)
    end if

    if( parameter_list%isPresent(metis_option_contig_key) ) then
       assert(parameter_list%isAssignable(metis_option_contig_key, this%metis_option_contig))
       istat = parameter_list%get(key = metis_option_contig_key , value = this%metis_option_contig)
       assert(istat==0)
    end if

    if( parameter_list%isPresent(metis_option_ctype_key) ) then
       assert(parameter_list%isAssignable(metis_option_ctype_key, this%metis_option_ctype))
       istat = parameter_list%get(key = metis_option_ctype_key  , value = this%metis_option_ctype)
       assert(istat==0)
    end if
  end subroutine mesh_partitioner_set_parameters_from_pl
  
  !=============================================================================
  subroutine mesh_partitioner_set_default_parameter_values(this)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    this%debug = mesh_partitioner_default_debug
    this%strat = mesh_partitioner_default_strat
    this%metis_option_ufactor = mesh_partitioner_default_metis_option_ufactor
    this%metis_option_minconn = mesh_partitioner_default_metis_option_minconn
    this%metis_option_contig = mesh_partitioner_default_metis_option_contig
    this%metis_option_ctype = mesh_partitioner_default_metis_option_ctype
    this%metis_option_iptype = mesh_partitioner_default_metis_option_iptype
    this%metis_option_debug = mesh_partitioner_default_metis_option_debug
  end subroutine mesh_partitioner_set_default_parameter_values
  
  subroutine mesh_partitioner_write_mesh_parts_pl(this, parameter_list)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    type(ParameterList_t)    , intent(in)    :: parameter_list
    
  end subroutine mesh_partitioner_write_mesh_parts_pl
  
  subroutine mesh_partitioner_write_mesh_parts_dir_path_prefix(this, dir_path, prefix)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    character(len=*)         , intent(in)    :: dir_path
    character(len=*)         , intent(in)    :: prefix
    
  end subroutine mesh_partitioner_write_mesh_parts_dir_path_prefix
  
  !=============================================================================
  subroutine mesh_partitioner_free(this)
    implicit none
    class(mesh_partitioner_t), intent(inout) :: this
    if ( allocated(this%num_parts_x_level) ) call memfree(this%num_parts_x_level,__FILE__,__LINE__)
    call this%set_default_parameter_values()
  end subroutine mesh_partitioner_free
  
end module mesh_partitioner_names
