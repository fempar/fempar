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
module output_handler_cell_fe_function_names
    use types_names
    use list_types_names
    use hash_table_names
    use allocatable_array_names

    use reference_fe_names
    use fe_space_names
    use fe_function_names
    use environment_names
    use field_names
  
    ! Linear algebra
    use vector_names
    use serial_scalar_array_names

    use output_handler_fe_field_names
    use output_handler_patch_names
  
implicit none
# include "debug.i90"
private
  
  
    type :: output_handler_cell_fe_function_t
    private
        type(fe_accessor_t), pointer           :: current_fe => NULL()
        type(quadrature_t),        allocatable :: quadratures(:)
        type(fe_map_t),            allocatable :: fe_maps(:)
        type(volume_integrator_t), allocatable :: volume_integrators(:)
        type(hash_table_ip_ip_t)               :: quadratures_and_maps_position ! Key = max_order_within_fe
        type(hash_table_ip_ip_t)               :: volume_integrators_position   ! Key = [max_order_within_fe,
                                                                                !       reference_fe_id]
        contains
            procedure, non_overridable :: create                              => output_handler_cell_fe_function_create
            procedure, non_overridable :: fill_patch                          => output_handler_cell_fe_function_fill_patch
            procedure, non_overridable :: free                                => output_handler_cell_fe_function_free
            
            procedure, non_overridable, private :: fill_patch_field           => output_handler_cell_fe_function_fill_patch_field
            procedure, non_overridable, private :: generate_vol_integ_pos_key => output_handler_cell_fe_function_generate_vol_integ_pos_key
            procedure, non_overridable, private :: get_number_reference_fes   => output_handler_cell_fe_function_get_number_reference_fes
            procedure, non_overridable, private :: get_quadrature             => output_handler_cell_fe_function_get_quadrature
            procedure, non_overridable, private :: get_fe_map                 => output_handler_cell_fe_function_get_fe_map
            procedure, non_overridable, private :: get_volume_integrator      => output_handler_cell_fe_function_get_volume_integrator      
    end type output_handler_cell_fe_function_t
  
public :: output_handler_cell_fe_function_t
  
contains

!  ! Includes with all the TBP and supporting subroutines for the types above.
!  ! In a future, we would like to use the submodule features of FORTRAN 2008.

    subroutine output_handler_cell_fe_function_create ( this, fe_space, num_refinements )
    !-----------------------------------------------------------------
    !< Create output_handler_cell_fe_function
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t), intent(inout) :: this
        class(serial_fe_space_t),           intent(in)    :: fe_space
        integer(ip), optional                             :: num_refinements
        type(fe_iterator_t)                               :: fe_iterator
        type(fe_accessor_t)                               :: fe
        class(reference_fe_t),            pointer         :: reference_fe
        class(lagrangian_reference_fe_t), pointer         :: reference_fe_geo
        class(environment_t),             pointer         :: environment
        integer(ip)                                       :: current_quadrature_and_map
        integer(ip)                                       :: current_volume_integrator
        integer(ip)                                       :: max_order_within_fe, max_order_field_id
        integer(ip)                                       :: vol_integ_pos_key
        integer(ip)                                       :: istat, field_id, quadrature_and_map_pos
        integer(ip)                                       :: reference_fe_id
    !-----------------------------------------------------------------
        environment => fe_space%get_environment()
        if (environment%am_i_l1_task()) then

            allocate ( this%quadratures(fe_space%get_number_reference_fes()), stat=istat); check (istat==0)
            allocate ( this%fe_maps(fe_space%get_number_reference_fes()), stat=istat); check (istat==0)
            allocate ( this%volume_integrators(fe_space%get_number_reference_fes()), stat=istat); check (istat==0)

            ! Create quadratures, fe_maps, and volume_integrators
            call this%quadratures_and_maps_position%init()
            call this%volume_integrators_position%init()
            current_quadrature_and_map = 1
            current_volume_integrator  = 1
            fe_iterator = fe_space%create_fe_iterator()
            do while ( .not. fe_iterator%has_finished() ) 
                call fe_iterator%current(fe)
                reference_fe_geo => fe%get_reference_fe_geo()

                max_order_within_fe = fe%get_max_order()

                call this%quadratures_and_maps_position%put(key = max_order_within_fe, &
                                                            val = current_quadrature_and_map, &
                                                            stat = istat)
                if (istat == now_stored) then
                    ! Create quadrature and fe_map associated to current max_order_within_fe
                    call reference_fe_geo%create_data_out_quadrature(num_refinements = max_order_within_fe-1, &
                                                                     quadrature      = this%quadratures(current_quadrature_and_map))
                    call this%fe_maps(current_quadrature_and_map)%create(this%quadratures(current_quadrature_and_map),&
                                                                         reference_fe_geo)
                    current_quadrature_and_map = current_quadrature_and_map + 1
                end if
                do field_id=1, fe_space%get_number_fields()
                    vol_integ_pos_key = this%generate_vol_integ_pos_key(fe_space%get_number_reference_fes(), &
                                                                        max_order_within_fe, &
                                                                        fe%get_reference_fe_id(field_id))
                    call this%volume_integrators_position%put(key=vol_integ_pos_key, &
                                                              val=current_volume_integrator, &
                                                              stat=istat)
                    if (istat == now_stored) then
                        call this%quadratures_and_maps_position%get(key = max_order_within_fe, &
                                                                    val = quadrature_and_map_pos, &
                                                                    stat = istat)
                        assert ( istat == key_found )
                        call this%volume_integrators(current_volume_integrator)%create(this%quadratures(quadrature_and_map_pos),&
                                                                                       fe%get_reference_fe(field_id))
                        current_volume_integrator = current_volume_integrator + 1
                    end if
                end do
                call fe_iterator%next()
            end do
        end if
    end subroutine output_handler_cell_fe_function_create


    subroutine output_handler_cell_fe_function_fill_patch(this, fe_accessor, number_fields, fe_fields, patch)
    !-----------------------------------------------------------------
    !< Fill a patch given a fe_accessor
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t),  intent(inout) :: this
        type(fe_accessor_t),               target, intent(in)    :: fe_accessor
        integer(ip),                               intent(in)    :: number_fields
        type(output_handler_fe_field_t),           intent(in)    :: fe_fields(1:number_fields)
        type(output_handler_patch_t),              intent(inout) :: patch
        integer(ip)                                              :: reference_fe_id
        integer(ip)                                              :: number_field
        integer(ip)                                              :: max_order_within_fe
        class(serial_fe_space_t),          pointer               :: fe_space
        class(lagrangian_reference_fe_t),  pointer               :: reference_fe_geo
        class(environment_t),              pointer               :: environment
        type(point_t),                     pointer               :: coordinates(:)
        type(fe_map_t),                    pointer               :: fe_map
        type(quadrature_t),                pointer               :: quadrature
        type(output_handler_patch_field_t),pointer               :: patch_field
        type(allocatable_array_ip2_t),     pointer               :: patch_subcells_connectivity
    !-----------------------------------------------------------------
        this%current_fe => fe_accessor
        fe_space => fe_accessor%get_fe_space()
        environment => fe_space%get_environment()
        if (environment%am_i_l1_task()) then
            max_order_within_fe =  fe_accessor%get_max_order()
            reference_fe_geo    => fe_accessor%get_reference_fe_geo()
            fe_map              => this%get_fe_map()
            coordinates         => fe_map%get_coordinates()
            call this%current_fe%get_coordinates(coordinates)

            quadrature => this%get_quadrature()
            call fe_map%update(quadrature)

            call patch%set_number_dimensions(reference_fe_geo%get_number_dimensions())
            call patch%set_number_vertices_per_subcell(quadrature%get_number_quadrature_points())
            call patch%set_number_subcells(reference_fe_geo%get_number_subcells(num_refinements=max_order_within_fe-1))
            call patch%set_number_vertices_per_subcell(reference_fe_geo%get_number_vertices())

            fe_map      => this%get_fe_map()
            call patch%set_coordinates(fe_map%get_quadrature_points_coordinates())

            patch_subcells_connectivity => patch%get_subcells_connectivity()
            call patch_subcells_connectivity%create(reference_fe_geo%get_number_vertices(), &
                                                    reference_fe_geo%get_number_subcells(num_refinements=max_order_within_fe-1))
            call reference_fe_geo%get_subcells_connectivity(num_refinements=max_order_within_fe-1, &
                                                            connectivity=patch_subcells_connectivity%a)

            do number_field = 1, number_fields
                patch_field  => patch%get_field(number_field)
                call this%fill_patch_field(fe_fields(number_field), patch_field)
            end do
        end if
    end subroutine output_handler_cell_fe_function_fill_patch


    subroutine output_handler_cell_fe_function_fill_patch_field(this, fe_field, patch_field)
    !-----------------------------------------------------------------
    !< Fill the path field given a fe_field
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t), intent(inout) :: this
        type(output_handler_fe_field_t),          intent(in)    :: fe_field
        type(output_handler_patch_field_t),       intent(inout) :: patch_field
        integer(ip)                                             :: field_id
        integer(ip)                                             :: reference_fe_id
        class(reference_fe_t),                    pointer       :: reference_fe
        type(fe_function_t),                      pointer       :: fe_function
        type(fe_map_t),                           pointer       :: fe_map

        type(volume_integrator_t),                pointer       :: volume_integrator
        type(allocatable_array_rp1_t),            pointer       :: patch_field_nodal_values
        ! Values + gradients for scalar fields
        real(rp),             allocatable                       :: scalar_function_values(:)
        type(vector_field_t), allocatable                       :: scalar_function_gradients(:)
        type(allocatable_array_rp1_t),            pointer       :: patch_field_scalar_function_values
        type(allocatable_array_vector_field_t),   pointer       :: patch_field_scalar_function_gradients
            
        ! Values + gradients for vector fields
        type(vector_field_t), allocatable                       :: vector_function_values(:)
        type(tensor_field_t), allocatable                       :: vector_function_gradients(:)
        type(allocatable_array_vector_field_t),   pointer       :: patch_field_vector_function_values
        type(allocatable_array_tensor_field_t),   pointer       :: patch_field_vector_function_gradients
            
        ! Values for tensor fields (gradients not supported yet)
        type(tensor_field_t) , allocatable                :: tensor_function_values(:)
        type(allocatable_array_tensor_field_t),   pointer :: patch_field_tensor_function_values
    !-----------------------------------------------------------------
        field_id          =  fe_field%get_field_id()
        fe_function       => fe_field%get_fe_function()
        fe_map            => this%get_fe_map()
        volume_integrator => this%get_volume_integrator(field_id) 
        call volume_integrator%update(fe_map)

        reference_fe => this%current_fe%get_reference_fe(field_id)
        reference_fe_id = this%current_fe%get_reference_fe_id(field_id)

        ! Gather DoFs of current cell + field_id on nodal_values 
        patch_field_nodal_values => patch_field%get_nodal_values()
        call patch_field_nodal_values%create(reference_fe%get_number_shape_functions())
        call fe_function%gather_nodal_values(this%current_fe, field_id, patch_field_nodal_values%a)

        call patch_field%set_field_type(this%current_fe%get_field_type(field_id))
        select case(this%current_fe%get_field_type(field_id))
        case ( field_type_scalar )
            patch_field_scalar_function_values    => patch_field%get_scalar_function_values()
            patch_field_scalar_function_gradients => patch_field%get_scalar_function_gradients()
            call patch_field_scalar_function_values%move_alloc_out(scalar_function_values) 
            call patch_field_scalar_function_gradients%move_alloc_out(scalar_function_gradients) 
            ! Evaluate values and gradients at quadrature points
            call volume_integrator%evaluate_fe_function(patch_field_nodal_values%a, scalar_function_values)
            call volume_integrator%evaluate_gradient_fe_function(patch_field_nodal_values%a, scalar_function_gradients)
            call patch_field_scalar_function_values%move_alloc_in(scalar_function_values) 
            call patch_field_scalar_function_gradients%move_alloc_in(scalar_function_gradients) 
        case ( field_type_vector )
            patch_field_vector_function_values    => patch_field%get_vector_function_values()
            patch_field_vector_function_gradients => patch_field%get_vector_function_gradients()
            call patch_field_vector_function_values%move_alloc_out(vector_function_values) 
            call patch_field_vector_function_gradients%move_alloc_out(vector_function_gradients) 
            ! Evaluate values and gradients at quadrature points
            call volume_integrator%evaluate_fe_function(patch_field_nodal_values%a, vector_function_values)
            call volume_integrator%evaluate_gradient_fe_function(patch_field_nodal_values%a, vector_function_gradients)
            call patch_field_vector_function_values%move_alloc_in(vector_function_values) 
            call patch_field_vector_function_gradients%move_alloc_in(vector_function_gradients) 
        case ( field_type_tensor )
            patch_field_vector_function_values    => patch_field%get_vector_function_values()
            call patch_field_tensor_function_values%move_alloc_out(tensor_function_values) 
            ! Evaluate values and gradients at quadrature points
            call volume_integrator%evaluate_fe_function(patch_field_nodal_values%a, tensor_function_values )
            call patch_field_tensor_function_values%move_alloc_in(tensor_function_values) 
        case default
            assert(.false.)
        end select
    end subroutine output_handler_cell_fe_function_fill_patch_field


    subroutine output_handler_cell_fe_function_free ( this )
    !-----------------------------------------------------------------
    !< Free procedure
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t), intent(inout) :: this
        integer(ip)                                             :: istat, i
    !-----------------------------------------------------------------
        call this%quadratures_and_maps_position%free()
        call this%volume_integrators_position%free()

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

        if(allocated(this%volume_integrators)) then
            do i=1, size(this%volume_integrators)
                call this%volume_integrators(i)%free()
            end do
            deallocate(this%volume_integrators, stat=istat)
            check(istat==0)
        end if

    end subroutine output_handler_cell_fe_function_free


    function output_handler_cell_fe_function_generate_vol_integ_pos_key (this, num_reference_fes, max_order_within_fe, reference_fe_id ) result(vol_integ_pos_key)
    !-----------------------------------------------------------------
    !< Generate vol_integ_pos_key
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t), intent(in) :: this
        integer(ip),                              intent(in) :: num_reference_fes
        integer(ip),                              intent(in) :: max_order_within_fe
        integer(ip),                              intent(in) :: reference_fe_id
        integer(ip)                                          :: vol_integ_pos_key
    !-----------------------------------------------------------------
        vol_integ_pos_key = reference_fe_id + (max_order_within_fe)*num_reference_fes
      end function output_handler_cell_fe_function_generate_vol_integ_pos_key


    function output_handler_cell_fe_function_get_quadrature ( this ) result(quadrature)
    !-----------------------------------------------------------------
    !< Return the quadrature
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t), target, intent(in) :: this
        type(quadrature_t),                       pointer            :: quadrature
        integer(ip)                                                  :: quadratures_position
        integer(ip)                                                  :: istat
    !-----------------------------------------------------------------
        assert ( associated(this%current_fe) )
        call this%quadratures_and_maps_position%get(key=this%current_fe%get_max_order(), &
             val=quadratures_position, &
             stat=istat)
        assert ( .not. istat == key_not_found )
        quadrature => this%quadratures(quadratures_position)
    end function output_handler_cell_fe_function_get_quadrature


    function output_handler_cell_fe_function_get_fe_map ( this ) result(fe_map)
    !-----------------------------------------------------------------
    !< Return the fe_map
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t), target, intent(in) :: this
        type(fe_map_t),                           pointer            :: fe_map
        integer(ip)                                                  :: fe_maps_position
        integer(ip)                                                  :: istat
    !-----------------------------------------------------------------
        assert ( associated(this%current_fe) )
        call this%quadratures_and_maps_position%get(key=this%current_fe%get_max_order(), &
             val=fe_maps_position, &
             stat=istat)
        assert ( .not. istat == key_not_found )
        fe_map => this%fe_maps(fe_maps_position)
    end function output_handler_cell_fe_function_get_fe_map


    function output_handler_cell_fe_function_get_volume_integrator ( this, field_id ) result(volume_integrator)
    !-----------------------------------------------------------------
    !< Return the volume integrator corresponding with the field_id
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t), target, intent(in) :: this
        integer(ip),                                      intent(in) :: field_id
        type(volume_integrator_t),                pointer            :: volume_integrator
        integer(ip)                                                  :: vol_integ_pos_key
        integer(ip)                                                  :: vol_integ_pos
        integer(ip)                                                  :: istat
    !-----------------------------------------------------------------
        assert ( associated(this%current_fe) )

        vol_integ_pos_key = &
             this%generate_vol_integ_pos_key(this%get_number_reference_fes(), &
             this%current_fe%get_max_order(), &
             this%current_fe%get_reference_fe_id(field_id))

        call this%volume_integrators_position%get(key=vol_integ_pos_key, &
             val=vol_integ_pos, &
             stat=istat)
        assert ( .not. istat == key_not_found )
        volume_integrator => this%volume_integrators(vol_integ_pos)
    end function output_handler_cell_fe_function_get_volume_integrator


    function output_handler_cell_fe_function_get_number_reference_fes ( this ) result(number_reference_fes)
    !-----------------------------------------------------------------
    !< Return the number of reference fes
    !-----------------------------------------------------------------
        class(output_handler_cell_fe_function_t), intent(in)   :: this
        integer(ip)                                            :: number_reference_fes
        class(serial_fe_space_t), pointer                      :: serial_fe_space
    !-----------------------------------------------------------------
        assert ( associated(this%current_fe) )
        serial_fe_space => this%current_fe%get_fe_space()
        number_reference_fes = serial_fe_space%get_number_reference_fes()
    end function output_handler_cell_fe_function_get_number_reference_fes

end module output_handler_cell_fe_function_names
