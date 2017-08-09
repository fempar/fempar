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
! [[output_handler_cell_fe_function_names(module)]]
!---------------------------------------------------------------------
module output_handler_cell_fe_function_names
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
! [[output_handler_cell_fe_function_t(type)]] 
!---------------------------------------------------------------------
    use types_names
    use list_types_names
    use hash_table_names
    use allocatable_array_names

    use reference_fe_names
    use base_static_triangulation_names
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

    type :: fill_patch_field_procedure_t
    !-----------------------------------------------------------------
    !*Author: Víctor Sande
    ! Date: 2016-11-29
    ! Version: 0.0.1
    ! Category: IO
    ! 
    !-----------------------------------------------------------------
    !### Derived type containing a pointer to a procedure interface
    !-----------------------------------------------------------------
        procedure(fill_patch_field_interface), nopass, pointer :: p => NULL()
    end type
  
    type :: output_handler_cell_fe_function_t
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
        class(fe_iterator_t), pointer                  :: current_fe            => NULL()
        type(quadrature_t),        allocatable         :: quadratures(:)
        type(fe_map_t),            allocatable         :: fe_maps(:)
        type(cell_integrator_t), allocatable         :: cell_integrators(:)
        type(hash_table_ip_ip_t)                       :: quadratures_and_maps_position ! Key = max_order_within_fe
        type(hash_table_ip_ip_t)                       :: cell_integrators_position   ! Key = [max_order_within_fe,
                                                                                        !       reference_fe_id]
        type(fill_patch_field_procedure_t), allocatable:: fill_patch_field(:)

        contains
        private
            procedure, non_overridable, public :: create                    => output_handler_cell_fe_function_create
            procedure, non_overridable, public :: get_num_nodes          => output_handler_cell_fe_function_get_num_nodes
            procedure, non_overridable, public :: get_num_cells          => output_handler_cell_fe_function_get_num_cells
            procedure, non_overridable, public :: get_num_dims     => output_handler_cell_fe_function_get_num_dims
            procedure, non_overridable, public :: has_mixed_cell_topologies => output_handler_cell_fe_function_has_mixed_cell_topologies
            procedure, non_overridable, public :: fill_patch                => output_handler_cell_fe_function_fill_patch
            procedure, non_overridable, public :: free                      => output_handler_cell_fe_function_free

            ! Strategy procedures to fill patch field data
            procedure, non_overridable :: apply_fill_patch_field_strategy => &
                                                            output_handler_cell_fe_function_apply_fill_patch_field_strategy
            procedure, non_overridable :: fill_patch_scalar_field_val     => &
                                                            output_handler_cell_fe_function_fill_patch_scalar_field_val
            procedure, non_overridable :: fill_patch_scalar_field_grad    => &
                                                            output_handler_cell_fe_function_fill_patch_scalar_field_grad
            procedure, non_overridable :: fill_patch_vector_field_val     => &
                                                            output_handler_cell_fe_function_fill_patch_vector_field_val
            procedure, non_overridable :: fill_patch_vector_field_grad    => &
                                                            output_handler_cell_fe_function_fill_patch_vector_field_grad
            procedure, non_overridable :: fill_patch_vector_field_div     => &
                                                            output_handler_cell_fe_function_fill_patch_vector_field_div
            procedure, non_overridable :: fill_patch_vector_field_curl    => &
                                                            output_handler_cell_fe_function_fill_patch_vector_field_curl
            procedure, non_overridable :: fill_patch_tensor_field_val     => &
                                                            output_handler_cell_fe_function_fill_patch_tensor_field_val
            procedure, non_overridable :: fill_patch_cell_vector => output_handler_cell_fe_function_fill_patch_cell_vector

            procedure, non_overridable :: generate_cell_integ_pos_key => output_handler_cell_fe_function_generate_cell_integ_pos_key
            procedure, non_overridable :: get_num_reference_fes   => output_handler_cell_fe_function_get_num_reference_fes
            procedure, non_overridable :: get_quadrature             => output_handler_cell_fe_function_get_quadrature
            procedure, non_overridable :: get_fe_map                 => output_handler_cell_fe_function_get_fe_map
            procedure, non_overridable :: get_cell_integrator      => output_handler_cell_fe_function_get_cell_integrator      
    end type output_handler_cell_fe_function_t


    interface 
        subroutine fill_patch_field_interface(this, fe_function, field_id, patch_field)
            import ip
            import fe_function_t
            import output_handler_patch_field_t
            import output_handler_cell_fe_function_t
            class(output_handler_cell_fe_function_t), intent(inout) :: this
            type(fe_function_t),                      intent(in)    :: fe_function
            integer(ip),                              intent(in)    :: field_id
            type(output_handler_patch_field_t),       intent(inout) :: patch_field
        end subroutine fill_patch_field_interface
    end interface
  
public :: output_handler_cell_fe_function_t
  
contains

!  ! Includes with all the TBP and supporting subroutines for the types above.
!  ! In a future, we would like to use the submodule features of FORTRAN 2008.

    subroutine output_handler_cell_fe_function_create ( this, fe, num_fields, fe_fields, num_refinements )
    !-----------------------------------------------------------------
    !< Create output_handler_cell_fe_function. 
    !< This procedure must be called every time the mesh changes.
    !< It loops over the finite elements to precalculate some values
    !< like *num_dims*, *num_cells*, *num_nodes*, 
    !< *quadratures*, *mixed_cell_topologies*, etc.
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t), intent(inout) :: this
        class(fe_iterator_t)                    , intent(inout) :: fe
        integer(ip),                              intent(in)    :: num_fields
        type(output_handler_fe_field_t),          intent(in)    :: fe_fields(1:num_fields)
        integer(ip), optional,                    intent(in)    :: num_refinements
        class(base_static_triangulation_t), pointer             :: triangulation
        class(reference_fe_t),            pointer               :: reference_fe
        class(lagrangian_reference_fe_t), pointer               :: reference_fe_geo
        class(lagrangian_reference_fe_t), pointer               :: previous_reference_fe_geo
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
            triangulation => fe_space%get_triangulation()
            this%num_dims = triangulation%get_num_dims()
            this%num_cells      = 0
            this%num_nodes      = 0

            allocate ( this%quadratures(fe_space%get_num_reference_fes()), stat=istat); check (istat==0)
            allocate ( this%fe_maps(fe_space%get_num_reference_fes()), stat=istat); check (istat==0)
            allocate ( this%cell_integrators(fe_space%get_num_reference_fes()), stat=istat); check (istat==0)

            ! Create quadratures, fe_maps, and cell_integrators
            call this%quadratures_and_maps_position%init()
            call this%cell_integrators_position%init()
            current_quadrature_and_map = 1
            current_cell_integrator  = 1
            
            nullify(previous_reference_fe_geo)
            call fe%first()
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
                    ! Create quadrature and fe_map associated to current max_order_within_fe
                    call reference_fe_geo%create_data_out_quadrature(num_refinements = max_order-1, &
                                                                     quadrature      = this%quadratures(current_quadrature_and_map))
                    call this%fe_maps(current_quadrature_and_map)%create(this%quadratures(current_quadrature_and_map),&
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
                if ( fe%is_local() ) then
                    ! Call to "max()" in order to take into account 
                    ! reference_fe_void_t (defined with order == -1)
                    max_order = max(fe%get_max_order_all_fields(),1)
                    ! Local cell and node counter
                    this%num_cells = this%num_cells + reference_fe_geo%get_num_subcells(max_order-1)
                    this%num_nodes = this%num_nodes + &
                            (reference_fe_geo%get_num_subcells(max_order-1)*reference_fe_geo%get_num_vertices())

                    ! Check if there are several topology types or a single one
                    !if(associated(previous_reference_fe_geo) .and.                       &
                    !    .not. same_type_as(previous_reference_fe_geo, reference_fe_geo)) &
                    !        this%mixed_cell_topologies = .true.
                    if(associated(previous_reference_fe_geo)) then
                        if(.not. same_type_as(previous_reference_fe_geo, reference_fe_geo)) &
                            this%mixed_cell_topologies = .true.
                    end if
                    previous_reference_fe_geo => reference_fe_geo
                endif
                call fe%next()
            end do
            ! Configure fill_patch_field strategy for each field
            if(allocated(this%fill_patch_field)) deallocate(this%fill_patch_field)
            allocate(this%fill_patch_field(num_fields))
            do num_field = 1, num_fields
                field_type    = fe_fields(num_field)%get_field_type()
                diff_operator = fe_fields(num_field)%get_diff_operator()
                call this%apply_fill_patch_field_strategy(field_type, diff_operator, this%fill_patch_field(num_field)%p)
            end do
        end if
    end subroutine output_handler_cell_fe_function_create


    function output_handler_cell_fe_function_get_num_nodes(this) result(num_nodes)
    !-----------------------------------------------------------------
    !< Return the number of nodes
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t),  intent(in) :: this
        integer(ip)                                           :: num_nodes
    !-----------------------------------------------------------------
        num_nodes = this%num_nodes
    end function output_handler_cell_fe_function_get_num_nodes


    function output_handler_cell_fe_function_get_num_cells(this) result(num_cells)
    !-----------------------------------------------------------------
    !< Return the number of cells
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t),  intent(in) :: this
        integer(ip)                                           :: num_cells
    !-----------------------------------------------------------------
        num_cells = this%num_cells
    end function output_handler_cell_fe_function_get_num_cells


    function output_handler_cell_fe_function_get_num_dims(this) result(num_dims)
    !-----------------------------------------------------------------
    !< Return the number of dimensions
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t),  intent(in) :: this
        integer(ip)                                           :: num_dims
    !-----------------------------------------------------------------
        num_dims = this%num_dims
    end function output_handler_cell_fe_function_get_num_dims


    function output_handler_cell_fe_function_has_mixed_cell_topologies(this) result(mixed_cell_topologies)
    !-----------------------------------------------------------------
    !< Return if the topology is composed by mixed cell types
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t),  intent(in) :: this
        logical                                               :: mixed_cell_topologies
    !-----------------------------------------------------------------
        mixed_cell_topologies = this%mixed_cell_topologies
    end function output_handler_cell_fe_function_has_mixed_cell_topologies


    subroutine output_handler_cell_fe_function_fill_patch(this, fe_iterator, num_fields, fe_fields, num_cell_vectors, cell_vectors, patch)
    !-----------------------------------------------------------------
    !< Fill a [[output_handler_patch_t(type)]] from a given [[fe_iterator_t(type)]].
    !< The **pach** contains a local view of the coordinates, connectivities 
    !< and field data per cell.
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t),  intent(inout) :: this
        class(fe_iterator_t),             target,  intent(in)    :: fe_iterator
        integer(ip),                               intent(in)    :: num_fields
        type(output_handler_fe_field_t),           intent(in)    :: fe_fields(1:num_fields)
        integer(ip),                               intent(in)    :: num_cell_vectors
        type(output_handler_cell_vector_t),        intent(in)    :: cell_vectors(1:num_cell_vectors)
        type(output_handler_patch_t),              intent(inout) :: patch
        integer(ip)                                              :: reference_fe_id
        integer(ip)                                              :: idx
        integer(ip)                                              :: field_id
        integer(ip)                                              :: max_order_within_fe
        class(serial_fe_space_t),          pointer               :: fe_space
        type(fe_function_t),               pointer               :: fe_function
        class(lagrangian_reference_fe_t),  pointer               :: reference_fe_geo
        class(environment_t),              pointer               :: environment
        type(point_t),                     pointer               :: coordinates(:)
        type(fe_map_t),                    pointer               :: fe_map
        type(quadrature_t),                pointer               :: quadrature
        type(output_handler_patch_field_t),pointer               :: patch_field
        type(allocatable_array_rp1_t),     pointer               :: patch_cell_vector
        type(allocatable_array_ip2_t),     pointer               :: patch_subcells_connectivity
        character(len=:), allocatable                            :: field_type
        character(len=:), allocatable                            :: diff_operator
    !-----------------------------------------------------------------
        this%current_fe => fe_iterator
        fe_space => fe_iterator%get_fe_space()
        environment => fe_space%get_environment()
        if (environment%am_i_l1_task()) then
            max_order_within_fe = max(fe_iterator%get_max_order_all_fields(),1)
            reference_fe_geo    => fe_iterator%get_reference_fe_geo()
            fe_map              => this%get_fe_map()
            coordinates         => fe_map%get_coordinates()
            call this%current_fe%get_coordinates(coordinates)

            quadrature => this%get_quadrature()
            call fe_map%update(quadrature)

            ! Set subcell information into patch
            call patch%set_cell_type(reference_fe_geo%get_topology())
            call patch%set_num_dims(reference_fe_geo%get_num_dims())
            call patch%set_num_vertices_per_subcell(quadrature%get_num_quadrature_points())
            call patch%set_num_subcells(reference_fe_geo%get_num_subcells(num_refinements=max_order_within_fe-1))
            call patch%set_num_vertices_per_subcell(reference_fe_geo%get_num_vertices())

            ! Set patch coordinates from fe_map
            call patch%set_coordinates(fe_map%get_quadrature_points_coordinates())

            ! Set patch connectivities from reference_fe_geo given num_refinements
            patch_subcells_connectivity => patch%get_subcells_connectivity()
            call patch_subcells_connectivity%create(reference_fe_geo%get_num_vertices(), &
                                                    reference_fe_geo%get_num_subcells(num_refinements=max_order_within_fe-1))
            call reference_fe_geo%get_subcells_connectivity(num_refinements=max_order_within_fe-1, &
                                                            connectivity=patch_subcells_connectivity%a)

            ! Fill patch fe field data
            do idx = 1, num_fields
                fe_function  => fe_fields(idx)%get_fe_function()
                field_id     =  fe_fields(idx)%get_field_id()
                patch_field  => patch%get_field(idx)

                assert(associated(this%fill_patch_field(idx)%p))
                call this%fill_patch_field(idx)%p(this, fe_function, field_id, patch_field)
            end do

            ! Fill patch cell vectors data
            do idx = 1, num_cell_vectors
                patch_cell_vector  => patch%get_cell_vector(idx)
                call this%fill_patch_cell_vector(cell_vectors(idx), patch%get_num_subcells(), patch_cell_vector)
            end do
        end if
    end subroutine output_handler_cell_fe_function_fill_patch


    subroutine output_handler_cell_fe_function_apply_fill_patch_field_strategy(this, field_type, diff_operator, proc)
    !-----------------------------------------------------------------
    !< Choose strategy to fill a patch field.
    !< Patch field calculation is distributed in several procedures
    !< in order to apply some diff operators.
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t),       intent(inout) :: this
        character(len=:), allocatable,                  intent(in)    :: field_type
        character(len=:), allocatable,                  intent(in)    :: diff_operator
        procedure(fill_patch_field_interface), pointer, intent(inout) :: proc
    !-----------------------------------------------------------------
        select case(field_type)
            ! Select procedures to fill patch field from a scalar field  (Value or Grad)
            case ( field_type_scalar )
                select case (diff_operator)
                    case (no_diff_operator)
                        proc => output_handler_cell_fe_function_fill_patch_scalar_field_val
                    case (grad_diff_operator)
                        proc => output_handler_cell_fe_function_fill_patch_scalar_field_grad
                    case DEFAULT
                        check(.false.)
                end select
            ! Select procedures to fill patch field from a vector field (Value or Grad or Div or Curl)
            case ( field_type_vector )
                select case (diff_operator)
                    case (no_diff_operator)
                        proc => output_handler_cell_fe_function_fill_patch_vector_field_val
                    case (grad_diff_operator)
                        proc => output_handler_cell_fe_function_fill_patch_vector_field_grad
                    case (div_diff_operator)
                        proc => output_handler_cell_fe_function_fill_patch_vector_field_div
                    case (curl_diff_operator)
                        proc => output_handler_cell_fe_function_fill_patch_vector_field_curl
                    case DEFAULT
                        check(.false.)
                end select
            ! Select procedures to fill patch field from a tensor field (only Value)
            case ( field_type_tensor )
                select case (diff_operator)
                    case (no_diff_operator)
                        proc => output_handler_cell_fe_function_fill_patch_tensor_field_val
                    case DEFAULT
                        check(.false.)
                end select
        end select
    end subroutine output_handler_cell_fe_function_apply_fill_patch_field_strategy


    subroutine output_handler_cell_fe_function_fill_patch_scalar_field_val(this, fe_function, field_id, patch_field)
    !-----------------------------------------------------------------
    !< Fill the [[output_handler_patch_field_t(type)]] with field values given a scalar fe_field
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t), intent(inout) :: this
        type(fe_function_t),                      intent(in)    :: fe_function
        integer(ip),                              intent(in)    :: field_id
        type(output_handler_patch_field_t),       intent(inout) :: patch_field
        class(reference_fe_t),                    pointer       :: reference_fe
        type(fe_map_t),                           pointer       :: fe_map
        type(cell_integrator_t),                pointer       :: cell_integrator
        type(allocatable_array_rp1_t),            pointer       :: patch_field_nodal_values
        real(rp),              allocatable                      :: scalar_function_values(:)
        type(allocatable_array_rp1_t),            pointer       :: patch_field_scalar_function_values
    !-----------------------------------------------------------------
        ! Get reference_Fe
        reference_fe => this%current_fe%get_reference_fe(field_id)
        assert(reference_fe%get_field_type() == field_type_scalar)

        ! Get and Update volume integrator
        fe_map            => this%get_fe_map()
        cell_integrator => this%get_cell_integrator(field_id) 
        call cell_integrator%update(fe_map)

        ! Gather DoFs of current cell + field_id on nodal_values 
        patch_field_nodal_values => patch_field%get_nodal_values()
        call patch_field_nodal_values%create(reference_fe%get_num_shape_functions())
        call fe_function%gather_nodal_values(this%current_fe, field_id, patch_field_nodal_values%a)

        ! Calculate scalar field values
        call patch_field%set_field_type(field_type_scalar)
        patch_field_scalar_function_values => patch_field%get_scalar_function_values()
        call patch_field_scalar_function_values%move_alloc_out(scalar_function_values) 
        call cell_integrator%evaluate_fe_function(patch_field_nodal_values%a, scalar_function_values)
        call patch_field_scalar_function_values%move_alloc_in(scalar_function_values) 

    end subroutine output_handler_cell_fe_function_fill_patch_scalar_field_val


    subroutine output_handler_cell_fe_function_fill_patch_scalar_field_grad(this, fe_function, field_id, patch_field)
    !-----------------------------------------------------------------
    !< Fill the [[output_handler_patch_field_t(type)]] with field gradients given a scalar fe_field
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t), intent(inout) :: this
        type(fe_function_t),                      intent(in)    :: fe_function
        integer(ip),                              intent(in)    :: field_id
        type(output_handler_patch_field_t),       intent(inout) :: patch_field
        class(reference_fe_t),                    pointer       :: reference_fe
        type(fe_map_t),                           pointer       :: fe_map
        type(cell_integrator_t),                pointer       :: cell_integrator
        type(allocatable_array_rp1_t),            pointer       :: patch_field_nodal_values
        type(vector_field_t),  allocatable                      :: vector_function_values(:)
        type(allocatable_array_vector_field_t),   pointer       :: patch_field_vector_function_values
    !-----------------------------------------------------------------
        ! Get reference_Fe
        reference_fe => this%current_fe%get_reference_fe(field_id)
        assert(reference_fe%get_field_type() == field_type_scalar)

        ! Get and Update volume integrator
        fe_map            => this%get_fe_map()
        cell_integrator => this%get_cell_integrator(field_id) 
        call cell_integrator%update(fe_map)

        ! Gather DoFs of current cell + field_id on nodal_values 
        patch_field_nodal_values => patch_field%get_nodal_values()
        call patch_field_nodal_values%create(reference_fe%get_num_shape_functions())
        call fe_function%gather_nodal_values(this%current_fe, field_id, patch_field_nodal_values%a)

        ! Calculate scalar field gradients
        call patch_field%set_field_type(field_type_vector)
        patch_field_vector_function_values => patch_field%get_vector_function_values()
        call patch_field_vector_function_values%move_alloc_out(vector_function_values) 
        call cell_integrator%evaluate_gradient_fe_function(patch_field_nodal_values%a, vector_function_values)
        call patch_field_vector_function_values%move_alloc_in(vector_function_values) 

    end subroutine output_handler_cell_fe_function_fill_patch_scalar_field_grad


    subroutine output_handler_cell_fe_function_fill_patch_vector_field_val(this, fe_function, field_id, patch_field)
    !-----------------------------------------------------------------
    !< Fill the [[output_handler_patch_field_t(type)]] with field values given a vector fe_field
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t), intent(inout) :: this
        type(fe_function_t),                      intent(in)    :: fe_function
        integer(ip),                              intent(in)    :: field_id
        type(output_handler_patch_field_t),       intent(inout) :: patch_field
        class(reference_fe_t),                    pointer       :: reference_fe
        type(fe_map_t),                           pointer       :: fe_map
        type(cell_integrator_t),                pointer       :: cell_integrator
        type(allocatable_array_rp1_t),            pointer       :: patch_field_nodal_values
        type(vector_field_t),  allocatable                      :: vector_function_values(:)
        type(allocatable_array_vector_field_t),   pointer       :: patch_field_vector_function_values
    !-----------------------------------------------------------------
        ! Get reference_Fe
        reference_fe => this%current_fe%get_reference_fe(field_id)
        assert(reference_fe%get_field_type() == field_type_vector)

        ! Get and Update volume integrator
        fe_map            => this%get_fe_map()
        cell_integrator => this%get_cell_integrator(field_id) 
        call cell_integrator%update(fe_map)

        ! Gather DoFs of current cell + field_id on nodal_values 
        patch_field_nodal_values => patch_field%get_nodal_values()
        call patch_field_nodal_values%create(reference_fe%get_num_shape_functions())
        call fe_function%gather_nodal_values(this%current_fe, field_id, patch_field_nodal_values%a)

        ! Calculate vector field values
        call patch_field%set_field_type(field_type_vector)
        patch_field_vector_function_values => patch_field%get_vector_function_values()
        call patch_field_vector_function_values%move_alloc_out(vector_function_values) 
        call cell_integrator%evaluate_fe_function(patch_field_nodal_values%a, vector_function_values)
        call patch_field_vector_function_values%move_alloc_in(vector_function_values) 

    end subroutine output_handler_cell_fe_function_fill_patch_vector_field_val


    subroutine output_handler_cell_fe_function_fill_patch_vector_field_grad(this, fe_function, field_id, patch_field)
    !-----------------------------------------------------------------
    !< Fill the [[output_handler_patch_field_t(type)]] with field gradients given a vector fe_field
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t), intent(inout) :: this
        type(fe_function_t),                      intent(in)    :: fe_function
        integer(ip),                              intent(in)    :: field_id
        type(output_handler_patch_field_t),       intent(inout) :: patch_field
        class(reference_fe_t),                    pointer       :: reference_fe
        type(fe_map_t),                           pointer       :: fe_map
        type(cell_integrator_t),                pointer       :: cell_integrator
        type(allocatable_array_rp1_t),            pointer       :: patch_field_nodal_values
        type(tensor_field_t),  allocatable                      :: tensor_function_values(:)
        type(allocatable_array_tensor_field_t),   pointer       :: patch_field_tensor_function_values
    !-----------------------------------------------------------------
        ! Get reference_Fe
        reference_fe => this%current_fe%get_reference_fe(field_id)
        assert(reference_fe%get_field_type() == field_type_vector)

        ! Get and Update volume integrator
        fe_map            => this%get_fe_map()
        cell_integrator => this%get_cell_integrator(field_id) 
        call cell_integrator%update(fe_map)

        ! Gather DoFs of current cell + field_id on nodal_values 
        patch_field_nodal_values => patch_field%get_nodal_values()
        call patch_field_nodal_values%create(reference_fe%get_num_shape_functions())
        call fe_function%gather_nodal_values(this%current_fe, field_id, patch_field_nodal_values%a)

        ! Calculate vector field gradients
        call patch_field%set_field_type(field_type_tensor)
        patch_field_tensor_function_values => patch_field%get_tensor_function_values()
        call patch_field_tensor_function_values%move_alloc_out(tensor_function_values) 
        call cell_integrator%evaluate_gradient_fe_function(patch_field_nodal_values%a, tensor_function_values)
        call patch_field_tensor_function_values%move_alloc_in(tensor_function_values) 

    end subroutine output_handler_cell_fe_function_fill_patch_vector_field_grad


    subroutine output_handler_cell_fe_function_fill_patch_vector_field_div(this, fe_function, field_id, patch_field)
    !-----------------------------------------------------------------
    !< Fill the [[output_handler_patch_field_t(type)]] with field divergence given a vector fe_field
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t), intent(inout) :: this
        type(fe_function_t),                      intent(in)    :: fe_function
        integer(ip),                              intent(in)    :: field_id
        type(output_handler_patch_field_t),       intent(inout) :: patch_field
        class(reference_fe_t),                    pointer       :: reference_fe
        type(quadrature_t),                       pointer       :: quadrature
        type(fe_map_t),                           pointer       :: fe_map
        type(cell_integrator_t),                pointer       :: cell_integrator
        type(allocatable_array_rp1_t),            pointer       :: patch_field_nodal_values
        real(rp),              allocatable                      :: scalar_function_values(:)
        type(tensor_field_t),  allocatable                      :: tensor_function_values(:)
        type(allocatable_array_rp1_t),            pointer       :: patch_field_scalar_function_values
        type(allocatable_array_tensor_field_t),   pointer       :: patch_field_tensor_function_values
        integer(ip)                                             :: qpoint
        integer(ip)                                             :: dim
    !-----------------------------------------------------------------
        ! Get reference_Fe
        reference_fe    => this%current_fe%get_reference_fe(field_id)
        quadrature      => this%get_quadrature()
        assert(reference_fe%get_field_type() == field_type_vector)

        ! Get and Update volume integrator
        fe_map            => this%get_fe_map()
        cell_integrator => this%get_cell_integrator(field_id) 
        call cell_integrator%update(fe_map)

        ! Gather DoFs of current cell + field_id on nodal_values 
        patch_field_nodal_values => patch_field%get_nodal_values()
        call patch_field_nodal_values%create(reference_fe%get_num_shape_functions())
        call fe_function%gather_nodal_values(this%current_fe, field_id, patch_field_nodal_values%a)

        call patch_field%set_field_type(field_type_scalar)

        ! get scalar and tensor function values
        patch_field_scalar_function_values => patch_field%get_scalar_function_values()
        patch_field_tensor_function_values => patch_field%get_tensor_function_values()
        call patch_field_scalar_function_values%move_alloc_out(scalar_function_values) 
        call patch_field_tensor_function_values%move_alloc_out(tensor_function_values) 

        ! Calculate gradients
        call cell_integrator%evaluate_gradient_fe_function(patch_field_nodal_values%a, tensor_function_values)

        ! Allocate scalar function values
        if ( allocated(scalar_function_values) ) then
           call memrealloc(quadrature%get_num_quadrature_points(), scalar_function_values, __FILE__, __LINE__)
        else
           call memalloc(quadrature%get_num_quadrature_points(), scalar_function_values, __FILE__, __LINE__)
        end if
        
        ! Calculate divergence
        scalar_function_values = 0._rp
        do qpoint = 1, quadrature%get_num_quadrature_points()
            do dim = 1, this%num_dims
                scalar_function_values(qpoint) = scalar_function_values(qpoint) + tensor_function_values(qpoint)%get(dim, dim)
            enddo
        enddo 

        ! return scalar and tensor function values
        call patch_field_scalar_function_values%move_alloc_in(scalar_function_values)
        call patch_field_tensor_function_values%move_alloc_in(tensor_function_values)

    end subroutine output_handler_cell_fe_function_fill_patch_vector_field_div


    subroutine output_handler_cell_fe_function_fill_patch_vector_field_curl(this, fe_function, field_id, patch_field)
    !-----------------------------------------------------------------
    !< Fill the [[output_handler_patch_field_t(type)]] with field gradients given a vector fe_field
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t), intent(inout) :: this
        type(fe_function_t),                      intent(in)    :: fe_function
        integer(ip),                              intent(in)    :: field_id
        type(output_handler_patch_field_t),       intent(inout) :: patch_field
        class(reference_fe_t),                    pointer       :: reference_fe
        type(quadrature_t),                       pointer       :: quadrature
        type(fe_map_t),                           pointer       :: fe_map
        type(cell_integrator_t),                pointer       :: cell_integrator
        type(allocatable_array_rp1_t),            pointer       :: patch_field_nodal_values
        real(rp),              allocatable                      :: scalar_function_values(:)
        type(vector_field_t),  allocatable                      :: vector_function_values(:)
        type(tensor_field_t),  allocatable                      :: tensor_function_values(:)
        type(allocatable_array_rp1_t),            pointer       :: patch_field_scalar_function_values
        type(allocatable_array_vector_field_t),   pointer       :: patch_field_vector_function_values
        type(allocatable_array_tensor_field_t),   pointer       :: patch_field_tensor_function_values
        integer(ip)                                             :: qpoint
        integer(ip)                                             :: istat
    !-----------------------------------------------------------------
        ! Get reference_Fe
        reference_fe    => this%current_fe%get_reference_fe(field_id)
        quadrature      => this%get_quadrature()
        assert(reference_fe%get_field_type() == field_type_vector)

        ! Get and Update volume integrator
        fe_map            => this%get_fe_map()
        cell_integrator => this%get_cell_integrator(field_id) 
        call cell_integrator%update(fe_map)

        ! Gather DoFs of current cell + field_id on nodal_values 
        patch_field_nodal_values => patch_field%get_nodal_values()
        call patch_field_nodal_values%create(reference_fe%get_num_shape_functions())
        call fe_function%gather_nodal_values(this%current_fe, field_id, patch_field_nodal_values%a)

        ! get tensor function values
        patch_field_tensor_function_values => patch_field%get_tensor_function_values()
        call patch_field_tensor_function_values%move_alloc_out(tensor_function_values) 

        ! Calculate gradients
        call cell_integrator%evaluate_gradient_fe_function(patch_field_nodal_values%a, tensor_function_values)

        if(this%num_dims == 2) then
            call patch_field%set_field_type(field_type_scalar)
            patch_field_scalar_function_values => patch_field%get_scalar_function_values()
            call patch_field_scalar_function_values%move_alloc_out(scalar_function_values) 
            ! Allocate scalar function values
            if ( allocated(scalar_function_values) ) then
               call memrealloc(quadrature%get_num_quadrature_points(), scalar_function_values, __FILE__, __LINE__)
            else
               call memalloc(quadrature%get_num_quadrature_points(), scalar_function_values, __FILE__, __LINE__)
            end if
            ! Calculate curl
            do qpoint = 1, quadrature%get_num_quadrature_points()
                scalar_function_values(qpoint) = tensor_function_values(qpoint)%get(1,2)-tensor_function_values(qpoint)%get(2,1)
            enddo
            call patch_field_scalar_function_values%move_alloc_in(scalar_function_values)     

        elseif(this%num_dims == 3) then
            call patch_field%set_field_type(field_type_vector)
            patch_field_vector_function_values => patch_field%get_vector_function_values()
            call patch_field_vector_function_values%move_alloc_out(vector_function_values) 
            ! Allocate vector function values
            if ( allocated(vector_function_values) ) then
               if ( size(vector_function_values) < quadrature%get_num_quadrature_points() ) then
                  deallocate(vector_function_values, stat=istat); check(istat==0)
                  allocate(vector_function_values(quadrature%get_num_quadrature_points()), stat=istat); check(istat==0)
               endif
            else
               allocate(vector_function_values(quadrature%get_num_quadrature_points()), stat=istat); check(istat==0)
            end if
            ! Calculate curl
            do qpoint = 1, quadrature%get_num_quadrature_points()
                call vector_function_values(qpoint)%set(1, &
                        tensor_function_values(qpoint)%get(2,3)-tensor_function_values(qpoint)%get(3,2))
                call vector_function_values(qpoint)%set(2, &
                        tensor_function_values(qpoint)%get(3,1)-tensor_function_values(qpoint)%get(1,3))
                call vector_function_values(qpoint)%set(3, &
                        tensor_function_values(qpoint)%get(1,2)-tensor_function_values(qpoint)%get(2,1))
            enddo
            call patch_field_vector_function_values%move_alloc_in(vector_function_values) 
        endif

        ! return tensor function values
        call patch_field_tensor_function_values%move_alloc_in(tensor_function_values)

    end subroutine output_handler_cell_fe_function_fill_patch_vector_field_curl


    subroutine output_handler_cell_fe_function_fill_patch_tensor_field_val(this, fe_function, field_id, patch_field)
    !-----------------------------------------------------------------
    !< Fill the [[output_handler_patch_field_t(type)]] with field values given a tensor field
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t), intent(inout) :: this
        type(fe_function_t),                      intent(in)    :: fe_function
        integer(ip),                              intent(in)    :: field_id
        type(output_handler_patch_field_t),       intent(inout) :: patch_field
        class(reference_fe_t),                    pointer       :: reference_fe
        type(fe_map_t),                           pointer       :: fe_map
        type(cell_integrator_t),                pointer       :: cell_integrator
        type(allocatable_array_rp1_t),            pointer       :: patch_field_nodal_values
        type(tensor_field_t),  allocatable                      :: tensor_function_values(:)
        type(allocatable_array_tensor_field_t),   pointer       :: patch_field_tensor_function_values
    !-----------------------------------------------------------------
        ! Get reference_Fe
        reference_fe => this%current_fe%get_reference_fe(field_id)
        assert(reference_fe%get_field_type() == field_type_tensor)

        ! Get and Update volume integrator
        fe_map            => this%get_fe_map()
        cell_integrator => this%get_cell_integrator(field_id) 
        call cell_integrator%update(fe_map)

        ! Gather DoFs of current cell + field_id on nodal_values 
        patch_field_nodal_values => patch_field%get_nodal_values()
        call patch_field_nodal_values%create(reference_fe%get_num_shape_functions())
        call fe_function%gather_nodal_values(this%current_fe, field_id, patch_field_nodal_values%a)

        ! Calculate tensor field values
        call patch_field%set_field_type(field_type_tensor)
        patch_field_tensor_function_values => patch_field%get_tensor_function_values()
        call patch_field_tensor_function_values%move_alloc_out(tensor_function_values) 
        call cell_integrator%evaluate_fe_function(patch_field_nodal_values%a, tensor_function_values )
        call patch_field_tensor_function_values%move_alloc_in(tensor_function_values)

    end subroutine output_handler_cell_fe_function_fill_patch_tensor_field_val


    subroutine output_handler_cell_fe_function_fill_patch_cell_vector(this, cell_vector, num_subcells, patch_cell_vector)
    !-----------------------------------------------------------------
    !< Fill the [[output_handler_patch_field_t(type)]] with field values given a **cell_vector**
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t), intent(inout) :: this
        type(output_handler_cell_vector_t),       intent(in)    :: cell_vector
        integer(ip),                              intent(in)    :: num_subcells
        type(allocatable_array_rp1_t),            intent(inout) :: patch_cell_vector
        real(rp), pointer                                       :: values(:)
    !-----------------------------------------------------------------
        ! Gather DoFs of current cell + field_id on nodal_values 
        values => cell_vector%get_cell_vector()
        call patch_cell_vector%create(num_subcells)
        patch_cell_vector%a = values(this%current_fe%get_lid())
    end subroutine output_handler_cell_fe_function_fill_patch_cell_vector


    subroutine output_handler_cell_fe_function_free ( this )
    !-----------------------------------------------------------------
    !< Free procedure
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t), intent(inout) :: this
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

        if(allocated(this%fe_maps)) then
            do i=1, size(this%fe_maps)
                call this%fe_maps(i)%free()
            end do
            deallocate(this%fe_maps, stat=istat)
            check(istat==0)
        end if

        if(allocated(this%cell_integrators)) then
            do i=1, size(this%cell_integrators)
                call this%cell_integrators(i)%free()
            end do
            deallocate(this%cell_integrators, stat=istat)
            check(istat==0)
        end if

        this%num_cells          = 0
        this%num_nodes          = 0
        this%num_dims     = 0
        this%mixed_cell_topologies = .false.
    end subroutine output_handler_cell_fe_function_free


    function output_handler_cell_fe_function_generate_cell_integ_pos_key (this, num_reference_fes, max_order_reference_fe_id, reference_fe_id ) result(cell_integ_pos_key)
    !-----------------------------------------------------------------
    !< Generate cell_integ_pos_key
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t), intent(in) :: this
        integer(ip),                              intent(in) :: num_reference_fes
        integer(ip),                              intent(in) :: max_order_reference_fe_id
        integer(ip),                              intent(in) :: reference_fe_id
        integer(ip)                                          :: cell_integ_pos_key
    !-----------------------------------------------------------------
        cell_integ_pos_key = reference_fe_id + (max_order_reference_fe_id)*num_reference_fes
      end function output_handler_cell_fe_function_generate_cell_integ_pos_key


    function output_handler_cell_fe_function_get_quadrature ( this ) result(quadrature)
    !-----------------------------------------------------------------
    !< Return the [[quadrature_t(type)]]
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t), target, intent(in) :: this
        type(quadrature_t),                       pointer            :: quadrature
        integer(ip)                                                  :: quadratures_position
        integer(ip)                                                  :: istat
    !-----------------------------------------------------------------
        assert ( associated(this%current_fe) )
        call this%quadratures_and_maps_position%get(key=this%current_fe%get_max_order_reference_fe_id(), &
             val=quadratures_position, &
             stat=istat)
        assert ( .not. istat == key_not_found )
        quadrature => this%quadratures(quadratures_position)
    end function output_handler_cell_fe_function_get_quadrature


    function output_handler_cell_fe_function_get_fe_map ( this ) result(fe_map)
    !-----------------------------------------------------------------
    !< Return the [[fe_map_t(type)]]
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t), target, intent(in) :: this
        type(fe_map_t),                           pointer            :: fe_map
        integer(ip)                                                  :: fe_maps_position
        integer(ip)                                                  :: istat
    !-----------------------------------------------------------------
        assert ( associated(this%current_fe) )
        call this%quadratures_and_maps_position%get(key=this%current_fe%get_max_order_reference_fe_id(), &
             val=fe_maps_position, &
             stat=istat)
        assert ( .not. istat == key_not_found )
        fe_map => this%fe_maps(fe_maps_position)
    end function output_handler_cell_fe_function_get_fe_map


    function output_handler_cell_fe_function_get_cell_integrator ( this, field_id ) result(cell_integrator)
    !-----------------------------------------------------------------
    !< Return the [[cell_integrator_t(type)]] corresponding with the **field_id**
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t), target, intent(in) :: this
        integer(ip),                                      intent(in) :: field_id
        type(cell_integrator_t),                pointer            :: cell_integrator
        integer(ip)                                                  :: cell_integ_pos_key
        integer(ip)                                                  :: cell_integ_pos
        integer(ip)                                                  :: istat
    !-----------------------------------------------------------------
        assert ( associated(this%current_fe) )

        cell_integ_pos_key = &
             this%generate_cell_integ_pos_key(this%get_num_reference_fes(), &
             this%current_fe%get_max_order_reference_fe_id(), &
             this%current_fe%get_reference_fe_id(field_id))

        call this%cell_integrators_position%get(key=cell_integ_pos_key, &
             val=cell_integ_pos, &
             stat=istat)
        assert ( .not. istat == key_not_found )
        cell_integrator => this%cell_integrators(cell_integ_pos)
    end function output_handler_cell_fe_function_get_cell_integrator


    function output_handler_cell_fe_function_get_num_reference_fes ( this ) result(num_reference_fes)
    !-----------------------------------------------------------------
    !< Return the number of [[reference_fe_t(type)]]
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t), intent(in)   :: this
        integer(ip)                                            :: num_reference_fes
        class(serial_fe_space_t), pointer                      :: serial_fe_space
    !-----------------------------------------------------------------
        assert ( associated(this%current_fe) )
        serial_fe_space => this%current_fe%get_fe_space()
        num_reference_fes = serial_fe_space%get_num_reference_fes()
    end function output_handler_cell_fe_function_get_num_reference_fes

end module output_handler_cell_fe_function_names
