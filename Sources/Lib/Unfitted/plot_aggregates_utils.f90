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
module plot_aggregates_utils_names
  use types_names
  use list_types_names
  use triangulation_names
  
  implicit none
# include "debug.i90"
  private

  public :: colorize_aggregate_ids
  
contains

  !============================================================================
  subroutine colorize_aggregate_ids(triangulation, aggregate_ids)

    implicit none
    class(triangulation_t), intent(in)    :: triangulation
    integer(ip),            intent(inout) :: aggregate_ids(:)
    
    type(list_t) :: graph
    integer(ip), allocatable :: aggr2color(:)
    integer(ip) :: i

    call make_aggregate_ids_contiguous(aggregate_ids)
    call create_aggregates_dual_graph(triangulation, aggregate_ids,graph)
    call memalloc(graph%get_num_pointers(),aggr2color,__FILE__,__LINE__)
    call colorize_graph(graph,aggr2color)
    do i=1,size(aggregate_ids)
      aggregate_ids(i) = aggr2color(aggregate_ids(i))
    end do
    call memfree(aggr2color,__FILE__,__LINE__)
    call graph%free()

  end subroutine colorize_aggregate_ids

  !============================================================================
  subroutine make_aggregate_ids_contiguous(aggregate_ids)
    implicit none
    integer(ip), intent(inout) :: aggregate_ids(:)

    logical, allocatable :: used_vals(:)
    integer(ip), allocatable :: old2new(:)
    integer(ip) :: current
    integer(ip) :: i

    assert(minval(aggregate_ids)>=0)
    call memalloc(maxval(aggregate_ids)+1,used_vals,__FILE__,__LINE__)
    used_vals(:) = .false.

    do i = 1, size(aggregate_ids)
      used_vals( aggregate_ids(i) + 1 ) = .true.
    end do

    call memalloc(size(used_vals),old2new,__FILE__,__LINE__)
    old2new(:) = -1

    current = 1
    do i = 1, size(used_vals)
      if ( used_vals(i) ) then
        old2new(i) = current
        current = current + 1
      end if
    end do

    do i = 1, size(aggregate_ids)
      assert(old2new( aggregate_ids(i) + 1 )>0)
      aggregate_ids(i) = old2new( aggregate_ids(i) + 1 )
    end do

    call memfree(old2new,__FILE__,__LINE__)
    call memfree(used_vals,__FILE__,__LINE__)

  end subroutine make_aggregate_ids_contiguous

  !============================================================================
  subroutine create_aggregates_dual_graph(triangulation,aggregate_ids,graph)

    implicit none
    class(triangulation_t), intent(in) :: triangulation
    integer(ip), intent(in) :: aggregate_ids(:)
    type(list_t), intent(inout) :: graph

    class(cell_iterator_t), allocatable :: cell
    class(cell_iterator_t), allocatable :: cell_around
    class(vef_iterator_t),allocatable :: vef
    integer(ip) :: ivef
    integer(ip) :: icell_around
    logical, allocatable :: visited_neigs(:)
    type(list_iterator_t) :: iter
    integer(ip) :: aggr_id

    call memalloc(triangulation%get_num_local_cells(),visited_neigs,__FILE__,__LINE__)

    call graph%create(maxval(aggregate_ids))

    call triangulation%create_vef_iterator(vef)
    call triangulation%create_cell_iterator(cell_around)
    call triangulation%create_cell_iterator(cell)

    ! Count number of neighbors
    do aggr_id = 1,maxval(aggregate_ids)
      call cell%first()
      visited_neigs(:) = .false.
      do while (.not. cell%has_finished()) 
        if ( aggregate_ids(cell%get_gid()) == aggr_id ) then
          do ivef = 1,cell%get_num_vefs()
            call cell%get_vef(ivef,vef)
            do icell_around = 1,vef%get_num_cells_around()
              call vef%get_cell_around(icell_around,cell_around)
              if (  (aggregate_ids(cell_around%get_gid()) /=  aggr_id ) .and.&
                    ( .not. visited_neigs(cell_around%get_gid())      )      ) then
                call graph%sum_to_pointer_index(aggr_id,1)
                visited_neigs(cell_around%get_gid()) = .true.
              end if
            end do
          end do
        end if
        call cell%next()
      end do
    end do

    call graph%calculate_header()
    call graph%allocate_list_from_pointer()

    ! Count store the neighbors
    do aggr_id = 1,maxval(aggregate_ids)
      call cell%first()
      visited_neigs(:) = .false.
      iter = graph%create_iterator(aggr_id)
      do while (.not. cell%has_finished()) 
        if ( aggregate_ids(cell%get_gid()) == aggr_id ) then
          do ivef = 1,cell%get_num_vefs()
            call cell%get_vef(ivef,vef)
            do icell_around = 1,vef%get_num_cells_around()
              call vef%get_cell_around(icell_around,cell_around)
              if (  (aggregate_ids(cell_around%get_gid()) /=  aggr_id ) .and.&
                    ( .not. visited_neigs(cell_around%get_gid())      )      ) then
                call iter%set_current(aggregate_ids(cell_around%get_gid()))
                call iter%next()
                visited_neigs(cell_around%get_gid()) = .true.
              end if
            end do
          end do
        end if
        call cell%next()
      end do
    end do

    call triangulation%free_cell_iterator(cell)
    call triangulation%free_cell_iterator(cell_around)
    call triangulation%free_vef_iterator(vef)

    call memfree(visited_neigs,__FILE__,__LINE__)

  end subroutine create_aggregates_dual_graph

  !============================================================================
  subroutine colorize_graph (graph, node2color)

    implicit none
    type(list_t),     intent(in)    :: graph
    integer(ip),      intent(inout) :: node2color(:)

    integer(ip) :: node
    type(list_iterator_t) :: iter
    integer(ip), parameter :: max_num_colors = 20
    logical :: used_colors(max_num_colors)
    integer(ip) :: color

    massert(size(node2color)>=graph%get_num_pointers(),'Not enough space in the given vector')

    ! Initialize the output to a dummy value
    node2color(:) = 0

    ! Loop in all the nodes of the graph
    do node = 1, graph%get_num_pointers()

      ! Loop over all the neighbors of this node and take all the colors used by the neighbors
      used_colors(:) = .false.
      iter = graph%create_iterator(node)
      do while (.not. iter%is_upper_bound())
        color = node2color( iter%get_current() )
        if ( color > 0 ) used_colors( color ) = .true.
        call iter%next()
      end do

      ! Find the minimum not used color
      do color = 1,max_num_colors
        if ( .not. used_colors(color) ) then
          node2color(node) = color
          exit
        end if
      end do

      massert( node2color(node)>0, 'Color has not been assigned. Is max_num_colors large enough?' )

    end do

  end subroutine colorize_graph

end module plot_aggregates_utils_names
