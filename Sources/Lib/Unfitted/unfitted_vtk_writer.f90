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

module unfitted_vtk_writer_names
 use types_names
 use memor_names
 use field_names
 
 use reference_fe_names
 use triangulation_names
 use fe_space_names
 use environment_names


  use unfitted_triangulations_names
  use unfitted_fe_spaces_names
  use piecewise_cell_map_names
  use IR_Precision ! VTK_IO
  use Lib_VTK_IO ! VTK_IO

  implicit none
# include "debug.i90"
  private

  type :: unfitted_vtk_writer_t
    private

    ! Used by VTKIO
    integer(ip) :: Nn
    integer(ip) :: Ne
    real(rp),    allocatable, dimension(:) :: x
    real(rp),    allocatable, dimension(:) :: y
    real(rp),    allocatable, dimension(:) :: z
    integer(ip), allocatable, dimension(:) :: cell_type
    integer(ip), allocatable, dimension(:) :: offset
    integer(ip), allocatable, dimension(:) :: connect
    real(rp),    allocatable, dimension(:) :: cell_data
    real(rp),    allocatable, dimension(:) :: point_data
    real(rp),    allocatable, dimension(:) :: v_x
    real(rp),    allocatable, dimension(:) :: v_y
    real(rp),    allocatable, dimension(:) :: v_z
    real(rp),    allocatable, dimension(:) :: pid
    class(environment_t), pointer          :: environment => null()

  contains

    procedure, non_overridable :: attach_triangulation  => uvtkw_attach_triangulation
    procedure, non_overridable :: attach_boundary_faces => uvtkw_attach_boundary_faces
    procedure, non_overridable :: attach_fitted_faces => uvtkw_attach_fitted_faces
    procedure, non_overridable :: attach_boundary_quadrature_points => uvtkw_attach_boundary_quadrature_points
    procedure, non_overridable :: attach_facets_quadrature_points => uvtkw_attach_facets_quadrature_points
    procedure, non_overridable :: attach_fe_function    => uvtkw_attach_fe_function
    procedure, non_overridable :: write_to_vtk_file     => uvtkw_write_to_vtk_file
    procedure, non_overridable :: write_to_pvtk_file    => uvtkw_write_to_pvtk_file
    procedure, non_overridable :: write                 => uvtkw_write
    procedure, non_overridable :: free                  => uvtkw_free

  end type unfitted_vtk_writer_t

  public :: unfitted_vtk_writer_t

contains

!========================================================================================  
  subroutine  uvtkw_attach_triangulation(this,triangulation)

    implicit none
    class(unfitted_vtk_writer_t),   intent(inout) :: this
    class(triangulation_t), intent(in)    :: triangulation

    integer(ip) :: num_cells, num_cell_nodes, num_subcells, num_subcell_nodes, num_dime
    integer(ip) :: istat, icell, inode, ino, isubcell
    class(cell_iterator_t), allocatable  :: cell
    type(point_t), allocatable, dimension(:) :: cell_coords, subcell_coords
    integer(ip) :: the_cell_type, the_subcell_type
    integer(ip), allocatable :: nodes_vtk2fempar(:), nodesids(:)
    integer(ip) :: my_part_id
    
    call this%free()

    this%environment => triangulation%get_environment()
    if ( .not. this%environment%am_i_l1_task() ) return
    my_part_id = this%environment%get_l1_rank() + 1

    call triangulation%create_cell_iterator(cell)

    select type (triangulation)
    class is (serial_unfitted_triangulation_t)
      num_subcells = triangulation%get_total_num_subcells()
      num_subcell_nodes = triangulation%get_max_num_nodes_in_subcell()
    class is (par_unfitted_triangulation_t)
      num_subcells = triangulation%get_total_num_subcells()
      num_subcell_nodes = triangulation%get_max_num_nodes_in_subcell()
    class is (unfitted_p4est_serial_triangulation_t)
      num_subcells = triangulation%get_total_num_subcells()
      num_subcell_nodes = triangulation%get_max_num_nodes_in_subcell()
    class default
      check(.false.)
    end select

    num_dime = triangulation%get_num_dims()
    num_cells = triangulation%get_num_local_cells()
    num_cell_nodes = cell%get_num_nodes()
    this%Ne = num_cells + num_subcells
    this%Nn = num_cell_nodes*num_cells  + num_subcell_nodes*num_subcells

    call memalloc ( this%Nn, this%x, __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%y, __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%z, __FILE__, __LINE__ )
    call memalloc ( this%Ne, this%cell_type, __FILE__, __LINE__ )
    call memalloc ( this%Ne, this%offset   , __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%connect  , __FILE__, __LINE__ )
    call memalloc ( this%Ne, this%cell_data , __FILE__, __LINE__ )
    call memalloc ( this%Ne, this%pid       , __FILE__, __LINE__ )
    allocate ( cell_coords(1:num_cell_nodes), stat = istat ); check(istat == 0)
    allocate ( subcell_coords(1:num_subcell_nodes), stat = istat ); check(istat == 0)
    call memalloc ( num_cell_nodes, nodes_vtk2fempar, __FILE__, __LINE__ )
    call memalloc ( num_cell_nodes, nodesids        , __FILE__, __LINE__ )

    select case (num_dime)
      case(3)
        the_cell_type = 12_I1P
        the_subcell_type = 10_I1P
        nodes_vtk2fempar(:) = [1, 2 , 4, 3, 5, 6, 8, 7]
      case(2)
        the_cell_type = 9_I1P
        the_subcell_type = 5_I1P
        nodes_vtk2fempar(:) = [1, 2 , 4, 3]
      case default
      check(.false.)
    end select

    ! Fill date to be passed to vtkio
    call cell%first()
    icell = 1
    inode = 1
    do while ( .not. cell%has_finished() )

      ! Skip ghost elems
      if (cell%is_ghost()) then
        call cell%next(); cycle
      end if

      call cell%update_sub_triangulation()

      call cell%get_nodes_coordinates( cell_coords )

      do ino = 1, num_cell_nodes !TODO is it possible to avoid loops like this one?
        this%x(inode) = cell_coords(ino)%get(1)
        this%y(inode) = cell_coords(ino)%get(2)
        this%z(inode) = cell_coords(ino)%get(3)
        nodesids(ino) = inode
        inode = inode + 1
      end do

      this%connect(nodesids(:)) = nodesids( nodes_vtk2fempar(:) ) - 1

      this%offset(icell)        = inode - 1
      this%cell_type(icell)     = the_cell_type

      if ( cell%is_interior() ) then
        this%cell_data( icell )  = -1
      else if ( cell%is_cut() ) then
        this%cell_data( icell )  = 0
      else
        this%cell_data( icell )  = 1
      end if

      this%pid(icell) = real(my_part_id,kind=rp)

      icell = icell + 1

      do isubcell = 1, cell%get_num_subcells()
        call cell%get_phys_coords_of_subcell(isubcell,subcell_coords)

        do ino = 1, num_subcell_nodes
          this%x(inode) = subcell_coords(ino)%get(1)
          this%y(inode) = subcell_coords(ino)%get(2)
          this%z(inode) = subcell_coords(ino)%get(3)
          this%connect(inode) = inode - 1
          inode = inode + 1
        end do

        this%offset(icell)        = inode - 1
        this%cell_type(icell)     = the_subcell_type

        if ( cell%is_interior_subcell(isubcell) ) then
          this%cell_data( icell )  = -2
        else
          this%cell_data( icell )  = 2
        end if

        this%pid(icell) = real(my_part_id,kind=rp)

        icell = icell + 1

      end do

      call cell%next()
    end do


    if (num_dime == 2_ip) this%z(:) = 0


    deallocate ( cell_coords, stat = istat ); check(istat == 0)
    deallocate ( subcell_coords, stat = istat ); check(istat == 0)
    call memfree ( nodes_vtk2fempar, __FILE__, __LINE__ )
    call memfree ( nodesids, __FILE__, __LINE__ )
    call triangulation%free_cell_iterator(cell)

  end subroutine  uvtkw_attach_triangulation

!========================================================================================
  subroutine uvtkw_attach_boundary_faces( this, triangulation )
  
    implicit none
    class(unfitted_vtk_writer_t),   intent(inout) :: this
    class(triangulation_t), intent(in)    :: triangulation
  
    integer(ip) :: num_subfacets, num_subfacet_nodes, num_dime
    integer(ip) :: istat, iface, inode, ino, isubfacet
    class(cell_iterator_t), allocatable  :: cell
    type(point_t), allocatable, dimension(:) :: subfacet_coords
    integer(ip) :: the_subfacet_type

    call this%free()
  
    this%environment => triangulation%get_environment()
    if ( .not. this%environment%am_i_l1_task() ) return

    call triangulation%create_cell_iterator(cell)
    num_dime = triangulation%get_num_dims()
    
    select type (triangulation)
    class is (serial_unfitted_triangulation_t)
      num_subfacets = triangulation%get_total_num_subfacets()
      num_subfacet_nodes = triangulation%get_max_num_nodes_in_subfacet()
    class is (unfitted_p4est_serial_triangulation_t)
      num_subfacets = triangulation%get_total_num_subfacets()
      num_subfacet_nodes = triangulation%get_max_num_nodes_in_subfacet()
    class is (par_unfitted_triangulation_t)
      num_subfacets = triangulation%get_total_num_subfacets()
      num_subfacet_nodes = triangulation%get_max_num_nodes_in_subfacet()
    class default
      check(.false.)
    end select
    
    this%Ne = num_subfacets
    this%Nn = num_subfacet_nodes*num_subfacets
  
    call memalloc ( this%Nn, this%x, __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%y, __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%z, __FILE__, __LINE__ )
    call memalloc ( this%Ne, this%cell_type, __FILE__, __LINE__ )
    call memalloc ( this%Ne, this%offset   , __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%connect  , __FILE__, __LINE__ )
    allocate ( subfacet_coords(1:num_subfacet_nodes), stat = istat ); check(istat == 0)
  
    select case (num_dime)
      case(3)
        the_subfacet_type = 5_I1P
      case(2)
        the_subfacet_type = 3_I1P
      case default
      check(.false.)
    end select
  
    ! Fill date to be passed to vtkio
    call cell%first()
    iface = 1
    inode = 1
    do while ( .not. cell%has_finished() )
  
      ! Skip ghost elems
      if (cell%is_ghost()) then
        call cell%next(); cycle
      end if
  
      call cell%update_sub_triangulation()
  
      do isubfacet = 1, cell%get_num_subfacets()
        call cell%get_phys_coords_of_subfacet(isubfacet,subfacet_coords)
  
        do ino = 1, num_subfacet_nodes
          this%x(inode) = subfacet_coords(ino)%get(1)
          this%y(inode) = subfacet_coords(ino)%get(2)
          this%z(inode) = subfacet_coords(ino)%get(3)
          this%connect(inode) = inode - 1
          inode = inode + 1
        end do
  
        this%offset(iface)        = inode - 1
        this%cell_type(iface)     = the_subfacet_type
  
        iface = iface + 1
  
      end do
  
      call cell%next()
    end do
  
  
    if (num_dime == 2_ip) this%z(:) = 0
  
    deallocate ( subfacet_coords, stat = istat ); check(istat == 0)
    call triangulation%free_cell_iterator(cell)
  
  end subroutine uvtkw_attach_boundary_faces


!========================================================================================
  subroutine uvtkw_attach_fitted_faces( this, triangulation )
  
    implicit none
    class(unfitted_vtk_writer_t),   intent(inout) :: this
    class(triangulation_t), intent(in)    :: triangulation

    integer(ip) :: num_facets, num_facet_nodes, num_subfacets, num_subfacet_nodes, num_dime
    integer(ip) :: istat, ifacet, inode, ino, isubfacet
    class(vef_iterator_t), allocatable  :: vef
    type(point_t), allocatable, dimension(:) :: facet_coords, subfacet_coords
    integer(ip) :: the_facet_type, the_subfacet_type
    integer(ip) :: my_part_id
    
    call this%free()

    this%environment => triangulation%get_environment()
    if ( .not. this%environment%am_i_l1_task() ) return
    my_part_id = this%environment%get_l1_rank() + 1


    select type (triangulation)
    class is (serial_unfitted_triangulation_t)
      num_subfacets = triangulation%get_total_num_fitted_sub_facets()
      num_subfacet_nodes = triangulation%get_max_num_nodes_in_subfacet()
    class is (par_unfitted_triangulation_t)
      mcheck(.false.,'Not yet implemented')
      !!num_subfacets = triangulation%get_total_num_fitted_sub_facets()
      !!num_subfacet_nodes = triangulation%get_max_num_nodes_in_subfacet()
    class is (unfitted_p4est_serial_triangulation_t)
      mcheck(.false.,'Not yet implemented')
      !!num_subfacets = triangulation%get_total_num_fitted_sub_facets()
      !!num_subfacet_nodes = triangulation%get_max_num_nodes_in_subfacet()
    class default
      check(.false.)
    end select

    num_dime = triangulation%get_num_dims()
    num_facets = 0
    call triangulation%create_vef_iterator(vef)
    do while (.not. vef%has_finished())
      if (vef%is_facet()) then
        num_facets = num_facets + 1
        num_facet_nodes = vef%get_num_nodes()
      end if
      call vef%next()
    end do
    this%Ne = num_facets + num_subfacets
    this%Nn = num_facet_nodes*num_facets  + num_subfacet_nodes*num_subfacets

    call memalloc ( this%Nn, this%x, __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%y, __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%z, __FILE__, __LINE__ )
    call memalloc ( this%Ne, this%cell_type, __FILE__, __LINE__ )
    call memalloc ( this%Ne, this%offset   , __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%connect  , __FILE__, __LINE__ )
    call memalloc ( this%Ne, this%cell_data , __FILE__, __LINE__ )
    call memalloc ( this%Ne, this%pid       , __FILE__, __LINE__ )
    allocate ( facet_coords(1:num_facet_nodes), stat = istat ); check(istat == 0)
    allocate ( subfacet_coords(1:num_subfacet_nodes), stat = istat ); check(istat == 0)

    select case (num_dime)
      case(3)
        the_facet_type    = 9_I1P
        the_subfacet_type = 5_I1P
      case(2)
        the_facet_type    = 3_I1P
        the_subfacet_type = 3_I1P
      case default
      check(.false.)
    end select

    ! Fill date to be passed to vtkio
    call vef%first()
    ifacet = 1
    inode = 1
    do while ( .not. vef%has_finished() )

      if (.not. vef%is_facet()) then
        call vef%next(); cycle
      end if

      call vef%update_sub_triangulation()

      call vef%get_nodes_coordinates( facet_coords )

      do ino = 1, num_facet_nodes !TODO is it possible to avoid loops like this one?
        this%x(inode) = facet_coords(ino)%get(1)
        this%y(inode) = facet_coords(ino)%get(2)
        this%z(inode) = facet_coords(ino)%get(3)
        this%connect(inode) =  inode -1
        inode = inode + 1
      end do

      this%offset(ifacet)        = inode - 1
      this%cell_type(ifacet)     = the_facet_type

      if ( vef%is_interior() ) then
        this%cell_data( ifacet )  = -1
      else if ( vef%is_cut() ) then
        this%cell_data( ifacet )  = 0
      else
        this%cell_data( ifacet )  = 1
      end if

      this%pid(ifacet) = real(my_part_id,kind=rp)

      ifacet = ifacet + 1

      do isubfacet = 1, vef%get_num_subvefs()
        call vef%get_phys_coords_of_subvef(isubfacet,subfacet_coords)

        do ino = 1, num_subfacet_nodes
          this%x(inode) = subfacet_coords(ino)%get(1)
          this%y(inode) = subfacet_coords(ino)%get(2)
          this%z(inode) = subfacet_coords(ino)%get(3)
          this%connect(inode) = inode - 1
          inode = inode + 1
        end do

        this%offset(ifacet)        = inode - 1
        this%cell_type(ifacet)     = the_subfacet_type

        if ( vef%is_interior_subvef(isubfacet) ) then
          this%cell_data( ifacet )  = -2
        else
          this%cell_data( ifacet )  = 2
        end if

        this%pid(ifacet) = real(my_part_id,kind=rp)

        ifacet = ifacet + 1

      end do

      call vef%next()
    end do

    if (num_dime == 2_ip) this%z(:) = 0

    deallocate ( facet_coords, stat = istat ); check(istat == 0)
    deallocate ( subfacet_coords, stat = istat ); check(istat == 0)
    call triangulation%free_vef_iterator(vef)
  
  end subroutine uvtkw_attach_fitted_faces

!========================================================================================  
  subroutine  uvtkw_attach_fe_function(this,fe_function,fe_space)

    implicit none
    class(unfitted_vtk_writer_t),   intent(inout) :: this
    class(fe_function_t),                   intent(in)    :: fe_function
    class(serial_fe_space_t),      intent(in)    :: fe_space

    class(triangulation_t), pointer :: triangulation
    class(fe_cell_iterator_t), allocatable  :: fe
    class(reference_fe_t), pointer :: ref_fe
    type(quadrature_t) :: subcel_nodal_quad
    type(interpolation_t) :: fe_interpol
    integer(ip) :: ipoint, subcell, inode
    real(rp), allocatable :: nodal_vals(:), subelem_nodal_vals(:)
    real(rp), pointer :: subcell_coords(:,:)
    type(point_t), allocatable :: subcell_points(:)
    integer(ip) :: num_elem_nodes, num_subelem_nodes, num_dime, istat, idime

    call this%free()

    triangulation => fe_space%get_triangulation()

    this%environment => triangulation%get_environment()
    if ( .not. this%environment%am_i_l1_task() ) return

    select type(triangulation)
      class is(serial_unfitted_triangulation_t)
        call this%attach_triangulation(triangulation)
        num_subelem_nodes = triangulation%get_max_num_nodes_in_subcell()
        num_elem_nodes = triangulation%get_max_num_shape_functions()
      class is (par_unfitted_triangulation_t)
        call this%attach_triangulation(triangulation)
        num_subelem_nodes = triangulation%get_max_num_nodes_in_subcell()
        num_elem_nodes = triangulation%get_max_num_shape_functions()
      class is (unfitted_p4est_serial_triangulation_t)
        call this%attach_triangulation(triangulation)
        num_subelem_nodes = triangulation%get_max_num_nodes_in_subcell()
        num_elem_nodes = triangulation%get_max_num_shape_functions()
      class default
        check(.false.)
    end select


    num_dime = triangulation%get_num_dims()

    call memalloc ( num_elem_nodes, nodal_vals, __FILE__, __LINE__ )
    call memalloc ( num_subelem_nodes, subelem_nodal_vals, __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%point_data, __FILE__, __LINE__ )
    allocate(subcell_points(1:num_subelem_nodes),stat = istat); check(istat == 0)

    call subcel_nodal_quad%create(num_dime,num_subelem_nodes)
    subcell_coords => subcel_nodal_quad%get_coordinates()

    call fe_space%create_fe_cell_iterator(fe)

    ipoint = 1
    do while ( .not. fe%has_finished() )

       ! Skip ghost elems
       if (fe%is_ghost()) then
         call fe%next(); cycle
       end if

       ref_fe => fe%get_reference_fe(1) ! TODO we assume a single field

       ! Recover nodal values
       if (fe%is_exterior()) then
         nodal_vals(:) = 0.0
       else
         call fe_function%gather_nodal_values(fe,1,nodal_vals)!TODO we assume a single field
       end if

       do inode = 1,num_elem_nodes
         this%point_data(ipoint) = nodal_vals(inode)
         ipoint = ipoint + 1
       end do

       call fe%update_sub_triangulation()

       do subcell = 1, fe%get_num_subcells()

         ! Get the subcell values
         if (fe%is_interior_subcell(subcell)) then
           call fe%get_ref_coords_of_subcell(subcell,subcell_points)
           ! TODO this is always a nightmare
           do inode = 1, num_subelem_nodes
             do idime = 1, num_dime
               subcell_coords(idime,inode) = subcell_points(inode)%get(idime)
             end do
           end do
           ! TODO this can be done with a volume integrator
           call ref_fe%create_interpolation(subcel_nodal_quad,fe_interpol)
           call ref_fe%evaluate_fe_function_scalar(fe_interpol,nodal_vals,subelem_nodal_vals)
         else
           subelem_nodal_vals(:) = 0.0
         end if

         do inode = 1, num_subelem_nodes
            this%point_data(ipoint) = subelem_nodal_vals(inode)
            ipoint = ipoint + 1
         end do

       end do

       call fe%next()
    end do

    call memfree ( nodal_vals,         __FILE__, __LINE__ )
    call memfree ( subelem_nodal_vals, __FILE__, __LINE__ )
    deallocate(subcell_points,stat = istat); check(istat == 0)
    call subcel_nodal_quad%free()
    call fe_interpol%free()
    call fe_space%free_fe_cell_iterator(fe)


  end subroutine  uvtkw_attach_fe_function

!========================================================================================
  subroutine uvtkw_attach_boundary_quadrature_points( this, fe_space )
  
    implicit none
    class(unfitted_vtk_writer_t),   intent(inout) :: this
    class(serial_unfitted_fe_space_t),      intent(in)    :: fe_space
  
    class(fe_cell_iterator_t),allocatable :: fe
    type(quadrature_t), pointer :: quadrature
    type(point_t), pointer :: quadrature_points_coordinates(:)
    type(piecewise_cell_map_t),     pointer :: cell_map
    integer(ip) :: num_dime, num_subfacets, num_gp_subfacet
    integer(ip) :: qpoint, num_quad_points
    integer(ip) :: ipoint
    type(vector_field_t)  :: normal_vec
    integer(ip), parameter :: the_point_type = 1

    class(triangulation_t), pointer :: triangulation

    call this%free()

    triangulation => fe_space%get_triangulation()

    this%environment => triangulation%get_environment()
    if ( .not. this%environment%am_i_l1_task() ) return

    num_dime       = triangulation%get_num_dims()

    select type(triangulation)
      class is (serial_unfitted_triangulation_t)
        num_subfacets   = triangulation%get_total_num_subfacets()
      class default
      check(.false.)
    end select
  
    call fe_space%create_fe_cell_iterator(fe)

    num_gp_subfacet = 0
    do while ( .not. fe%has_finished() )
       call fe%update_boundary_integration()
       quadrature => fe%get_boundary_quadrature()
       num_gp_subfacet = quadrature%get_num_quadrature_points()
       if (num_gp_subfacet > 0) exit
       call fe%next()
    end do
  
    this%Ne = num_subfacets*num_gp_subfacet
    this%Nn = num_subfacets*num_gp_subfacet
  
    call memalloc ( this%Nn, this%x, __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%y, __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%z, __FILE__, __LINE__ )
    call memalloc ( this%Ne, this%cell_type, __FILE__, __LINE__ )
    call memalloc ( this%Ne, this%offset   , __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%connect  , __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%v_x, __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%v_y, __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%v_z, __FILE__, __LINE__ )
  
    ipoint = 1
    call fe%first()
    do while ( .not. fe%has_finished() )
  
       ! Skip ghost elems
       if (fe%is_ghost()) then
         call fe%next(); cycle
       end if
  
       ! Update FE-integration related data structures
       call fe%update_boundary_integration()
  
       ! As the quadrature changes elem by elem, this has to be inside the loop
       quadrature => fe%get_boundary_quadrature()
       num_quad_points = quadrature%get_num_quadrature_points()
       cell_map => fe%get_boundary_piecewise_cell_map()
  
       ! Physical coordinates of the quadrature points
       quadrature_points_coordinates => cell_map%get_quadrature_points_coordinates()
  
       do qpoint = 1, num_quad_points
         this%x(ipoint) = quadrature_points_coordinates(qpoint)%get(1)
         this%y(ipoint) = quadrature_points_coordinates(qpoint)%get(2)
         this%z(ipoint) = quadrature_points_coordinates(qpoint)%get(3)
         call cell_map%get_normal(qpoint,normal_vec)
         this%v_x(ipoint) = normal_vec%get(1)
         this%v_y(ipoint) = normal_vec%get(2)
         this%v_z(ipoint) = normal_vec%get(3)
         this%cell_type(ipoint) = the_point_type
         this%offset(ipoint) = ipoint
         this%connect(ipoint) = ipoint - 1
         ipoint = ipoint + 1
       end do
  
       call fe%next()
    end do
  
    if (num_dime == 2_ip) this%z(:) = 0
    if (num_dime == 2_ip) this%v_z(:) = 0
    call fe_space%free_fe_cell_iterator(fe)
  
  end subroutine uvtkw_attach_boundary_quadrature_points

!========================================================================================
  subroutine uvtkw_attach_facets_quadrature_points( this, fe_space )
  
    implicit none
    class(unfitted_vtk_writer_t),   intent(inout) :: this
    class(serial_unfitted_fe_space_t),      intent(inout) :: fe_space
  
    class(fe_facet_iterator_t),allocatable :: fe_facet
    type(quadrature_t), pointer :: quadrature
    type(point_t), pointer :: quadrature_points_coordinates(:)
    integer(ip) :: num_dime
    integer(ip) :: qpoint, num_quad_points
    integer(ip) :: ipoint
    type(vector_field_t)  :: normal_vec(2)
    integer(ip), parameter :: the_point_type = 1
    real(rp) :: mystatus

    class(triangulation_t), pointer :: triangulation

    call this%free()

    triangulation => fe_space%get_triangulation()

    this%environment => triangulation%get_environment()
    if ( .not. this%environment%am_i_l1_task() ) return

    num_dime       = triangulation%get_num_dims()

    call fe_space%set_up_facet_integration()
  
    call fe_space%create_fe_facet_iterator(fe_facet)

    this%Ne = 0
    do while ( .not. fe_facet%has_finished() )
      if (fe_facet%is_ghost()) then
        call fe_facet%next(); cycle
      end if
      call fe_facet%update_integration()
      quadrature => fe_facet%get_quadrature()
      this%Ne = this%Ne + quadrature%get_num_quadrature_points()
      call fe_facet%next()
    end do
  
    this%Nn = this%Ne
  
    call memalloc ( this%Nn, this%x, __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%y, __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%z, __FILE__, __LINE__ )
    call memalloc ( this%Ne, this%cell_type, __FILE__, __LINE__ )
    call memalloc ( this%Ne, this%offset   , __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%connect  , __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%v_x, __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%v_y, __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%v_z, __FILE__, __LINE__ )
    call memalloc ( this%Ne, this%cell_data , __FILE__, __LINE__ )
  
    ipoint = 1
    this%cell_data(:) = 0
    call fe_facet%first()
    do while ( .not. fe_facet%has_finished() )
  
       ! Skip ghost elems
       if (fe_facet%is_ghost()) then
         call fe_facet%next(); cycle
       end if
  
       ! Update FE-integration related data structures
       call fe_facet%update_integration()
  
       ! As the quadrature changes elem by elem, this has to be inside the loop
       quadrature => fe_facet%get_quadrature()
       num_quad_points = quadrature%get_num_quadrature_points()
  
       ! Physical coordinates of the quadrature points
       quadrature_points_coordinates => fe_facet%get_quadrature_points_coordinates()

       if (fe_facet%is_cut()) then
         mystatus = 0.0_rp
       else if (fe_facet%is_interior()) then
         mystatus = -1.0_rp
       else if (fe_facet%is_exterior()) then
         mystatus = 1.0_rp
       else
         check(.false.)
       end if
  
       do qpoint = 1, num_quad_points
         this%x(ipoint) = quadrature_points_coordinates(qpoint)%get(1)
         this%y(ipoint) = quadrature_points_coordinates(qpoint)%get(2)
         this%z(ipoint) = quadrature_points_coordinates(qpoint)%get(3)
         call fe_facet%get_normals(qpoint,normal_vec)
         this%v_x(ipoint) = normal_vec(1)%get(1)
         this%v_y(ipoint) = normal_vec(1)%get(2)
         this%v_z(ipoint) = normal_vec(1)%get(3)
         this%cell_data(ipoint) = mystatus
         this%cell_type(ipoint) = the_point_type
         this%offset(ipoint) = ipoint
         this%connect(ipoint) = ipoint - 1
         ipoint = ipoint + 1
       end do
  
       call fe_facet%next()
    end do
  
    if (num_dime == 2_ip) this%z(:) = 0
    if (num_dime == 2_ip) this%v_z(:) = 0
    call fe_space%free_fe_facet_iterator(fe_facet)
  
  end subroutine uvtkw_attach_facets_quadrature_points

!========================================================================================  
  subroutine uvtkw_write_to_vtk_file(this,filename)

    implicit none
    class(unfitted_vtk_writer_t), intent(in) :: this
    character(*), intent(in) :: filename

    integer(ip) :: E_IO

    assert(allocated(this%x        ))
    assert(allocated(this%y        ))
    assert(allocated(this%z        ))
    assert(allocated(this%connect  ))
    assert(allocated(this%offset   ))
    assert(allocated(this%cell_type))

    assert(associated(this%environment))
    assert(this%environment%am_i_l1_task())

    E_IO = VTK_INI_XML(output_format = 'ascii', filename = filename, mesh_topology = 'UnstructuredGrid')
    E_IO = VTK_GEO_XML(NN = this%Nn, NC = this%Ne, X = this%x, Y = this%y, Z = this%z)
    E_IO = VTK_CON_XML(NC = this%Ne, connect = this%connect, offset = this%offset, cell_type = int(this%cell_type,I1P) )

    E_IO = VTK_DAT_XML(var_location = 'cell', var_block_action = 'opeN')
    if (allocated(this%cell_data)) then
      E_IO = VTK_VAR_XML(NC_NN = this%Ne, varname = 'cell_scalars', var = this%cell_data)
    end if
    if (allocated(this%pid)) then
      E_IO = VTK_VAR_XML(NC_NN = this%Ne, varname = 'part_id', var = this%pid)
    end if
    E_IO = VTK_DAT_XML(var_location = 'cell', var_block_action = 'CLOSE')


    if (allocated(this%point_data)) then
      E_IO = VTK_DAT_XML(var_location = 'node', var_block_action = 'opeN')
      E_IO = VTK_VAR_XML(NC_NN = this%Nn, varname = 'point_scalars', var = this%point_data)
      E_IO = VTK_DAT_XML(var_location = 'node', var_block_action = 'CLOSE')
    end if

    if (allocated(this%v_x)) then
      E_IO = VTK_DAT_XML(var_location = 'node', var_block_action = 'opeN')
      E_IO = VTK_VAR_XML(NC_NN = this%Nn, varname = 'point_vectors', varX = this%v_x, varY = this%v_y, varZ = this%v_z )
      E_IO = VTK_DAT_XML(var_location = 'node', var_block_action = 'CLOSE')
    end if

    E_IO = VTK_GEO_XML()
    E_IO = VTK_END_XML()

  end subroutine uvtkw_write_to_vtk_file

!========================================================================================  
  subroutine uvtkw_write_to_pvtk_file(this,file_prefix,num_parts)

    implicit none
    class(unfitted_vtk_writer_t), intent(in) :: this
    character(*),                 intent(in) :: file_prefix
    integer(ip),                  intent(in) :: num_parts

    integer(ip) :: file_id
    integer(ip) :: E_IO
    integer(ip) :: i

    assert(associated(this%environment))
    assert(this%environment%am_i_l1_task())

    E_IO = 0

    E_IO = PVTK_INI_XML(filename      = trim(adjustl(file_prefix//'.pvtu')), &
                        mesh_topology = 'PUnstructuredGrid', &
                        tp            = 'Float64',           &
                        cf            = file_id)
    assert(E_IO == 0)
    do i=0, num_parts-1
        E_IO = PVTK_GEO_XML(source=trim(adjustl( trim(adjustl(file_prefix))//'_'//trim(adjustl(str(no_sign=.true.,n=i)))//'.vtu' )), cf=file_id)
        assert(E_IO == 0)
    enddo

    E_IO = PVTK_DAT_XML(var_location = 'Node', var_block_action = 'OPEN', cf=file_id)
    assert(E_IO == 0)

    if (allocated(this%point_data)) then
        E_IO = PVTK_VAR_XML(varname = 'point_scalars', &
                            tp      = 'Float64',       &
                            Nc      = 1,               &
                            cf      = file_id ); assert(E_IO == 0)
    end if

    if (allocated(this%v_x)) then
        E_IO = PVTK_VAR_XML(varname = 'point_vectors', &
                            tp      = 'Float64',       &
                            Nc      = 3,               &
                            cf      = file_id ); assert(E_IO == 0)
    end if

    E_IO = PVTK_DAT_XML(var_location = 'Node', var_block_action = 'CLOSE', cf=file_id)
    assert(E_IO == 0)

    E_IO = PVTK_DAT_XML(var_location = 'Cell', var_block_action = 'OPEN', cf=file_id)
    assert(E_IO == 0)

    if (allocated(this%cell_data)) then
        E_IO = PVTK_VAR_XML(varname = 'cell_scalars',  &
                            tp      = 'Float64',       &
                            Nc      = 1,               &
                            cf      = file_id ); assert(E_IO == 0)
    end if

    if (allocated(this%pid)) then
        E_IO = PVTK_VAR_XML(varname = 'part_id',  &
                            tp      = 'Float64',       &
                            Nc      = 1,               &
                            cf      = file_id ); assert(E_IO == 0)
    end if

    E_IO = PVTK_DAT_XML(var_location = 'Cell', var_block_action = 'CLOSE', cf=file_id)
    assert(E_IO == 0)

    E_IO = PVTK_END_XML(cf = file_id)
    assert(E_IO == 0)

  end subroutine uvtkw_write_to_pvtk_file

!========================================================================================  
  subroutine uvtkw_write(this,file_prefix)

    implicit none
    class(unfitted_vtk_writer_t), intent(in) :: this
    character(*), intent(in) :: file_prefix

    integer(ip)                                :: me, np

    assert(associated(this%environment))
    me = this%environment%get_l1_rank()
    np = this%environment%get_l1_size()

    if (this%environment%am_i_l1_task()) &
      call this%write_to_vtk_file(trim(adjustl( trim(adjustl(file_prefix))//'_'//trim(adjustl(str(no_sign=.true.,n=me)))//'.vtu' )))
    if (me == 0) &
      call this%write_to_pvtk_file(file_prefix,np)

  end subroutine uvtkw_write

!========================================================================================  
  subroutine uvtkw_free(this)

    implicit none
    class(unfitted_vtk_writer_t), intent(inout) :: this

    if (allocated(this%x         )) call memfree ( this%x          , __FILE__, __LINE__ )
    if (allocated(this%y         )) call memfree ( this%y          , __FILE__, __LINE__ )
    if (allocated(this%z         )) call memfree ( this%z          , __FILE__, __LINE__ )
    if (allocated(this%cell_type )) call memfree ( this%cell_type  , __FILE__, __LINE__ )
    if (allocated(this%offset    )) call memfree ( this%offset     , __FILE__, __LINE__ )
    if (allocated(this%connect   )) call memfree ( this%connect    , __FILE__, __LINE__ )
    if (allocated(this%cell_data )) call memfree ( this%cell_data  , __FILE__, __LINE__ )
    if (allocated(this%point_data)) call memfree ( this%point_data , __FILE__, __LINE__ )
    if (allocated(this%v_x       )) call memfree ( this%v_x        , __FILE__, __LINE__ )
    if (allocated(this%v_y       )) call memfree ( this%v_y        , __FILE__, __LINE__ )
    if (allocated(this%v_z       )) call memfree ( this%v_z        , __FILE__, __LINE__ )
    if (allocated(this%pid       )) call memfree ( this%pid        , __FILE__, __LINE__ )
    this%environment => null()

  end subroutine uvtkw_free

end module unfitted_vtk_writer_names
!***************************************************************************************************


