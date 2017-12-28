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
!---------------------------------------------------------------------
!*Author: Víctor Sande
! Date: 2016-11-28
! Version: 0.0.1
! Category: IO
!
!---------------------------------------------------------------------
!### Local view per cell for [[serial_fe_space_t(type)]] and [[fe_function_t(type)]] types.
!
! Contains the following public entities: 
! [[output_handler_fe_cell_function_names(module)]]
!---------------------------------------------------------------------
module output_handler_fe_cell_function_names
!---------------------------------------------------------------------
!* Author: Alberto F. Martín
! Date: 2016-11-29
! Version: 0.0.1
! Category: IO
! 
!---------------------------------------------------------------------
!### Local view per cell for [[serial_fe_space_t(type)]] and [[fe_function_t(type)]] types.
! 
! Contains the following public entities: 
! [[output_handler_fe_cell_function_t(type)]] 
!---------------------------------------------------------------------
    use types_names
    use list_types_names
    use hash_table_names
    use allocatable_array_names
    use std_vector_real_rp_names

    use reference_fe_names
    use triangulation_names
    use fe_space_names
    use environment_names
    use field_names
  
    ! Linear algebra
    use vector_names
    use serial_scalar_array_names

    use output_handler_fe_field_names
    use output_handler_patch_names
    use output_handler_parameters_names
  
implicit none
# include "debug.i90"
private  
  
    type :: output_handler_fe_cell_function_t
    !-----------------------------------------------------------------
    !*Author: Alberto F. Martín
    ! Date: 2016-11-29
    ! Version: 0.0.1
    ! Category: IO
    ! 
    !-----------------------------------------------------------------
    !### Local view per cell delimiter for [[serial_fe_space_t(type)]] and [[fe_function_t(type)]] types.
    !-----------------------------------------------------------------
    private
        logical                                        :: mixed_cell_topologies = .false.
        integer(ip)                                    :: num_dims     = 0
        integer(ip)                                    :: num_nodes          = 0
        integer(ip)                                    :: num_cells          = 0
        class(fe_cell_iterator_t),  pointer            :: fe => NULL()
        type(quadrature_t),        allocatable         :: quadratures(:)
        type(cell_map_t),            allocatable       :: cell_maps(:)
        type(cell_integrator_t), allocatable           :: cell_integrators(:)
        type(hash_table_ip_ip_t)                       :: quadratures_and_maps_position ! Key = max_order_within_fe
        type(hash_table_ip_ip_t)                       :: cell_integrators_position   ! Key = [max_order_within_fe,
                                                                                        !       reference_fe_id]
        contains
        private
            procedure, non_overridable, public :: create                    => output_handler_fe_cell_function_create
            procedure, non_overridable, public :: get_num_nodes             => output_handler_fe_cell_function_get_num_nodes
            procedure, non_overridable, public :: get_num_cells             => output_handler_fe_cell_function_get_num_cells
            procedure, non_overridable, public :: get_num_dims              => output_handler_fe_cell_function_get_num_dims
            procedure, non_overridable, public :: has_mixed_cell_topologies => output_handler_fe_cell_function_has_mixed_cell_topologies
            procedure, non_overridable, public :: free                      => output_handler_fe_cell_function_free
            procedure, non_overridable :: generate_cell_integ_pos_key       => output_handler_fe_cell_function_generate_cell_integ_pos_key
            procedure, non_overridable :: get_num_reference_fes             => output_handler_fe_cell_function_get_num_reference_fes
            procedure, non_overridable, public :: get_quadrature            => output_handler_fe_cell_function_get_quadrature
            procedure, non_overridable, public :: get_cell_map              => output_handler_fe_cell_function_get_cell_map
            procedure, non_overridable, public :: get_cell_integrator       => output_handler_fe_cell_function_get_cell_integrator 
            procedure, non_overridable, public :: get_fe                    => output_handler_get_fe
    end type output_handler_fe_cell_function_t
  
    public :: output_handler_fe_cell_function_t

    
contains

!  ! Includes with all the TBP and supporting subroutines for the types above.
!  ! In a future, we would like to use the submodule features of FORTRAN 2008.

    subroutine output_handler_fe_cell_function_create ( this, fe, num_refinements )
    !-----------------------------------------------------------------
    !< Create output_handler_fe_cell_function. 
    !< This procedure must be called every time the mesh changes.
    !< It loops over the finite elements to precalculate some values
    !< like *num_dims*, *num_cells*, *num_nodes*, 
    !< *quadratures*, *mixed_cell_topologies*, etc.
    !-----------------------------------------------------------------
        class(output_handler_fe_cell_function_t), intent(inout) :: this
        class(fe_cell_iterator_t),     target   , intent(in)    :: fe
        integer(ip), optional,                    intent(in)    :: num_refinements
        class(triangulation_t), pointer                         :: triangulation
        class(reference_fe_t),            pointer               :: reference_fe
        class(reference_fe_t), pointer                          :: reference_fe_geo
        class(reference_fe_t), pointer                          :: previous_reference_fe_geo
        class(environment_t),             pointer               :: environment
        integer(ip)                                             :: current_quadrature_and_map
        integer(ip)                                             :: current_cell_integrator
        integer(ip)                                             :: max_order, max_order_reference_fe_id, max_order_field_id
        integer(ip)                                             :: cell_integ_pos_key
        integer(ip)                                             :: istat, field_id, quadrature_and_map_pos
        integer(ip)                                             :: reference_fe_id
        integer(ip)                                             :: num_field
        character(len=:), allocatable                           :: field_type
        character(len=:), allocatable                           :: diff_operator
        class(serial_fe_space_t),  pointer                      :: fe_space
        
        !-----------------------------------------------------------------
        fe_space    => fe%get_fe_space()
        environment => fe_space%get_environment()
        if (environment%am_i_l1_task()) then
            call this%free()
            this%fe => fe
            triangulation => fe_space%get_triangulation()
            this%num_dims = triangulation%get_num_dims()
            this%num_cells      = 0
            this%num_nodes      = 0

            allocate ( this%quadratures(fe_space%get_num_reference_fes()), stat=istat); check (istat==0)
            allocate ( this%cell_maps(fe_space%get_num_reference_fes()), stat=istat); check (istat==0)
            allocate ( this%cell_integrators(fe_space%get_num_reference_fes()), stat=istat); check (istat==0)

            ! Create quadratures, cell_maps, and cell_integrators
            call this%quadratures_and_maps_position%init()
            call this%cell_integrators_position%init()
            current_quadrature_and_map = 1
            current_cell_integrator  = 1
            
            nullify(previous_reference_fe_geo)
            call this%fe%first()
            do while ( .not. fe%has_finished() ) 
                ! Call to "max()" in order to take into account 
                ! reference_fe_void_t (defined with order == -1)
                max_order = max(fe%get_max_order_all_fields(),1)
                reference_fe_geo => fe%get_reference_fe_geo()
                max_order_reference_fe_id = fe%get_max_order_reference_fe_id()
                call this%quadratures_and_maps_position%put(key = max_order_reference_fe_id, &
                                                            val = current_quadrature_and_map, &
                                                            stat = istat)
                if (istat == now_stored) then
                    ! Create quadrature and cell_map associated to current max_order_within_fe
                    call reference_fe_geo%create_data_out_quadrature(num_refinements = max_order-1, &
                                                                     quadrature      = this%quadratures(current_quadrature_and_map))
                    call this%cell_maps(current_quadrature_and_map)%create(this%quadratures(current_quadrature_and_map),&
                                                                         reference_fe_geo)
                    current_quadrature_and_map = current_quadrature_and_map + 1
                end if
                do field_id=1, fe_space%get_num_fields()
                    cell_integ_pos_key = this%generate_cell_integ_pos_key(fe_space%get_num_reference_fes(), &
                                                                        max_order_reference_fe_id, &
                                                                        fe%get_reference_fe_id(field_id))
                    call this%cell_integrators_position%put(key=cell_integ_pos_key, &
                                                              val=current_cell_integrator, &
                                                              stat=istat)
                    if (istat == now_stored) then
                        call this%quadratures_and_maps_position%get(key = max_order_reference_fe_id, &
                                                                    val = quadrature_and_map_pos, &
                                                                    stat = istat)
                        assert ( istat == key_found )
                        call this%cell_integrators(current_cell_integrator)%create(this%quadratures(quadrature_and_map_pos),&
                                                                                       fe%get_reference_fe(field_id))
                        current_cell_integrator = current_cell_integrator + 1
                    end if
                end do
                if ( this%fe%is_local() ) then
                    ! Call to "max()" in order to take into account 
                    ! reference_fe_void_t (defined with order == -1)
                    max_order = max(fe%get_max_order_all_fields(),1)
                    ! Local cell and node counter
                    this%num_cells = this%num_cells + reference_fe_geo%get_num_subcells(max_order-1)
                    this%num_nodes = this%num_nodes + &
                            (reference_fe_geo%get_num_subcells(max_order-1)*reference_fe_geo%get_num_vertices())

                    ! Check if there are several topology types or a single one
                    if(associated(previous_reference_fe_geo)) then
                       if  (.not. same_type_as(previous_reference_fe_geo, reference_fe_geo)) then
                            this%mixed_cell_topologies = .true.
                       end if
                    end if
                    previous_reference_fe_geo => reference_fe_geo
                 endif
                call this%fe%next()
            end do
            call this%fe%first()           
        end if
    end subroutine output_handler_fe_cell_function_create


    function output_handler_fe_cell_function_get_num_nodes(this) result(num_nodes)
    !-----------------------------------------------------------------
    !< Return the number of nodes
    !-----------------------------------------------------------------
        class(output_handler_fe_cell_function_t),  intent(in) :: this
        integer(ip)                                           :: num_nodes
    !-----------------------------------------------------------------
        num_nodes = this%num_nodes
    end function output_handler_fe_cell_function_get_num_nodes


    function output_handler_fe_cell_function_get_num_cells(this) result(num_cells)
    !-----------------------------------------------------------------
    !< Return the number of cells
    !-----------------------------------------------------------------
        class(output_handler_fe_cell_function_t),  intent(in) :: this
        integer(ip)                                           :: num_cells
    !-----------------------------------------------------------------
        num_cells = this%num_cells
    end function output_handler_fe_cell_function_get_num_cells


    function output_handler_fe_cell_function_get_num_dims(this) result(num_dims)
    !-----------------------------------------------------------------
    !< Return the number of dimensions
    !-----------------------------------------------------------------
        class(output_handler_fe_cell_function_t),  intent(in) :: this
        integer(ip)                                           :: num_dims
    !-----------------------------------------------------------------
        num_dims = this%num_dims
    end function output_handler_fe_cell_function_get_num_dims


    function output_handler_fe_cell_function_has_mixed_cell_topologies(this) result(mixed_cell_topologies)
    !-----------------------------------------------------------------
    !< Return if the topology is composed by mixed cell types
    !-----------------------------------------------------------------
        class(output_handler_fe_cell_function_t),  intent(in) :: this
        logical                                               :: mixed_cell_topologies
    !-----------------------------------------------------------------
        mixed_cell_topologies = this%mixed_cell_topologies
    end function output_handler_fe_cell_function_has_mixed_cell_topologies

    subroutine output_handler_fe_cell_function_free ( this )
    !-----------------------------------------------------------------
    !< Free procedure
    !-----------------------------------------------------------------
        class(output_handler_fe_cell_function_t), intent(inout) :: this
        integer(ip)                                             :: istat, i
    !-----------------------------------------------------------------
        call this%quadratures_and_maps_position%free()

        if(allocated(this%quadratures)) then
            do i=1, size(this%quadratures)
                call this%quadratures(i)%free()
            end do
            deallocate(this%quadratures, stat=istat)
            check(istat==0)
        end if

        if(allocated(this%cell_maps)) then
            do i=1, size(this%cell_maps)
                call this%cell_maps(i)%free()
            end do
            deallocate(this%cell_maps, stat=istat)
            check(istat==0)
        end if

        if(allocated(this%cell_integrators)) then
            do i=1, size(this%cell_integrators)
                call this%cell_integrators(i)%free()
            end do
            deallocate(this%cell_integrators, stat=istat)
            check(istat==0)
        end if
        nullify(this%fe)
        this%num_cells          = 0
        this%num_nodes          = 0
        this%num_dims     = 0
        this%mixed_cell_topologies = .false.
    end subroutine output_handler_fe_cell_function_free


    function output_handler_fe_cell_function_generate_cell_integ_pos_key (this, num_reference_fes, max_order_reference_fe_id, reference_fe_id ) result(cell_integ_pos_key)
    !-----------------------------------------------------------------
    !< Generate cell_integ_pos_key
    !-----------------------------------------------------------------
        class(output_handler_fe_cell_function_t), intent(in) :: this
        integer(ip),                              intent(in) :: num_reference_fes
        integer(ip),                              intent(in) :: max_order_reference_fe_id
        integer(ip),                              intent(in) :: reference_fe_id
        integer(ip)                                          :: cell_integ_pos_key
    !-----------------------------------------------------------------
        cell_integ_pos_key = reference_fe_id + (max_order_reference_fe_id)*num_reference_fes
      end function output_handler_fe_cell_function_generate_cell_integ_pos_key


    function output_handler_fe_cell_function_get_quadrature ( this) result(quadrature)
    !-----------------------------------------------------------------
    !< Return the [[quadrature_t(type)]]
    !-----------------------------------------------------------------
        class(output_handler_fe_cell_function_t), target, intent(in) :: this
        type(quadrature_t),                       pointer            :: quadrature
        integer(ip)                                                  :: quadratures_position
        integer(ip)                                                  :: istat
    !-----------------------------------------------------------------
        call this%quadratures_and_maps_position%get(key=this%fe%get_max_order_reference_fe_id(), &
             val=quadratures_position, &
             stat=istat)
        assert ( .not. istat == key_not_found )
        quadrature => this%quadratures(quadratures_position)
    end function output_handler_fe_cell_function_get_quadrature


    function output_handler_fe_cell_function_get_cell_map ( this) result(cell_map)
    !-----------------------------------------------------------------
    !< Return the [[cell_map_t(type)]]
    !-----------------------------------------------------------------
        class(output_handler_fe_cell_function_t), target, intent(in) :: this
        type(cell_map_t),                           pointer          :: cell_map
        integer(ip)                                                  :: cell_maps_position
        integer(ip)                                                  :: istat
    !-----------------------------------------------------------------
        call this%quadratures_and_maps_position%get(key=this%fe%get_max_order_reference_fe_id(), &
             val=cell_maps_position, &
             stat=istat)
        assert ( .not. istat == key_not_found )
        cell_map => this%cell_maps(cell_maps_position)
    end function output_handler_fe_cell_function_get_cell_map


    function output_handler_fe_cell_function_get_cell_integrator ( this, field_id ) result(cell_integrator)
    !-----------------------------------------------------------------
    !< Return the [[cell_integrator_t(type)]] corresponding with the **field_id**
    !-----------------------------------------------------------------
        class(output_handler_fe_cell_function_t), target, intent(in) :: this
        integer(ip),                                      intent(in) :: field_id
        type(cell_integrator_t),                pointer            :: cell_integrator
        integer(ip)                                                  :: cell_integ_pos_key
        integer(ip)                                                  :: cell_integ_pos
        integer(ip)                                                  :: istat
    !-----------------------------------------------------------------
        cell_integ_pos_key = &
             this%generate_cell_integ_pos_key(this%get_num_reference_fes(), &
             this%fe%get_max_order_reference_fe_id(), &
             this%fe%get_reference_fe_id(field_id))

        call this%cell_integrators_position%get(key=cell_integ_pos_key, &
             val=cell_integ_pos, &
             stat=istat)
        assert ( .not. istat == key_not_found )
        cell_integrator => this%cell_integrators(cell_integ_pos)
    end function output_handler_fe_cell_function_get_cell_integrator
    
    
   function output_handler_get_fe ( this ) result(fe_cell_iterator)
    !-----------------------------------------------------------------
    !< Return the [[cell_integrator_t(type)]] corresponding with the **field_id**
    !-----------------------------------------------------------------
        class(output_handler_fe_cell_function_t), target, intent(in) :: this
        class(fe_cell_iterator_t),         pointer                   :: fe_cell_iterator
        fe_cell_iterator => this%fe
    end function output_handler_get_fe
    
    


    function output_handler_fe_cell_function_get_num_reference_fes ( this ) result(num_reference_fes)
    !-----------------------------------------------------------------
    !< Return the number of [[reference_fe_t(type)]]
    !-----------------------------------------------------------------
        class(output_handler_fe_cell_function_t), intent(in)   :: this
        integer(ip)                                            :: num_reference_fes
        class(serial_fe_space_t), pointer                      :: serial_fe_space
    !-----------------------------------------------------------------
        serial_fe_space => this%fe%get_fe_space()
        num_reference_fes = serial_fe_space%get_num_reference_fes()
    end function output_handler_fe_cell_function_get_num_reference_fes

end module output_handler_fe_cell_function_names
