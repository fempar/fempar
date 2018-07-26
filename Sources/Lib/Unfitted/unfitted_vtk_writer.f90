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
 use stdio_names
 
 use reference_fe_names
 use triangulation_names
 use fe_space_names
 use environment_names
 use list_types_names
 use conditions_names
 use function_names


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
    procedure, non_overridable :: attach_vefs         => uvtkw_attach_vefs
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
    integer(ip), allocatable :: nodes_vtk2fempar(:), nodesids(:)
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
      num_subfacets = triangulation%get_total_num_fitted_sub_facets()
      num_subfacet_nodes = triangulation%get_max_num_nodes_in_subfacet()
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
    call memalloc ( num_facet_nodes, nodes_vtk2fempar, __FILE__, __LINE__ )
    call memalloc ( num_facet_nodes, nodesids        , __FILE__, __LINE__ )
    allocate ( facet_coords(1:num_facet_nodes), stat = istat ); check(istat == 0)
    allocate ( subfacet_coords(1:num_subfacet_nodes), stat = istat ); check(istat == 0)

    select case (num_dime)
      case(3)
        the_facet_type    = 9_I1P
        the_subfacet_type = 5_I1P
        nodes_vtk2fempar(:) = [1, 2 , 4, 3]
      case(2)
        the_facet_type    = 3_I1P
        the_subfacet_type = 3_I1P
        nodes_vtk2fempar(:) = [1, 2]
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
        nodesids(ino) = inode
        inode = inode + 1
      end do

      this%connect(nodesids(:)) = nodesids( nodes_vtk2fempar(:) ) - 1

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

    call memfree ( nodes_vtk2fempar, __FILE__, __LINE__ )
    call memfree ( nodesids, __FILE__, __LINE__ )
    deallocate ( facet_coords, stat = istat ); check(istat == 0)
    deallocate ( subfacet_coords, stat = istat ); check(istat == 0)
    call triangulation%free_vef_iterator(vef)
  
  end subroutine uvtkw_attach_fitted_faces

!========================================================================================  
  subroutine uvtkw_attach_vefs(this,triangulation,conditions)

    implicit none
    class(unfitted_vtk_writer_t),   intent(inout) :: this
    class(triangulation_t), intent(in)    :: triangulation
    class(conditions_t), optional, intent(in) :: conditions

    class(vef_iterator_t),allocatable :: vef
    class(cell_iterator_t), allocatable :: cell

    type(point_t), allocatable :: nodal_coords(:)
    integer(ip), allocatable :: nodal_ids(:)
    integer(ip) :: istat
    integer(ip) :: vef_lid
    class(reference_fe_t), pointer :: reference_fe
    type(list_iterator_t) :: nodal_iter
    integer(ip) :: num_dime
    integer(ip) :: inode, icell, ino
    integer(ip), parameter :: vef_type(0:2) = [1,3,9]
    integer(ip), parameter :: nodes_vtk2fempar_0d(1) = [1]
    integer(ip), parameter :: nodes_vtk2fempar_1d(2) = [1, 2]
    integer(ip), parameter :: nodes_vtk2fempar_2d(4) = [1, 2 , 4, 3]
    type(p_scalar_function_t) :: func
    class(scalar_function_t), pointer :: fun

    call this%free()

    this%environment => triangulation%get_environment()

    num_dime = triangulation%get_num_dims()

    call triangulation%create_vef_iterator(vef)
    call triangulation%create_cell_iterator(cell)
    allocate(nodal_coords(triangulation%get_max_num_shape_functions()),stat=istat); check(istat == 0)
    allocate(nodal_ids(triangulation%get_max_num_shape_functions()),stat=istat); check(istat == 0)

    ! Count number of nodes and number of cells
    this%Nn = 0
    this%Ne = 0
    call vef%first()
    do while ( .not. vef%has_finished() )
      call vef%get_cell_around(1,cell)
      vef_lid = cell%get_vef_lid_from_gid(vef%get_gid())
      reference_fe => cell%get_reference_fe()
      nodal_iter = reference_fe%create_dofs_on_n_face_iterator(vef_lid)
      this%Nn = this%Nn + nodal_iter%get_size()
      this%Ne = this%Ne + 1
      call vef%next()
    end do

    ! Allocate
    call memalloc ( this%Nn, this%x, __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%y, __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%z, __FILE__, __LINE__ )
    call memalloc ( this%Ne, this%cell_type, __FILE__, __LINE__ )
    call memalloc ( this%Ne, this%offset   , __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%connect  , __FILE__, __LINE__ )
    call memalloc ( this%Ne, this%cell_data , __FILE__, __LINE__ )
    if (present(conditions)) then
      call memalloc ( this%Nn, this%point_data  , __FILE__, __LINE__ )
      this%point_data(:) = 0.0
    end if

    call vef%first()
    inode = 1
    icell = 1
    do while ( .not. vef%has_finished() )

      call vef%get_cell_around(1,cell)
      call cell%get_nodes_coordinates(nodal_coords)
      vef_lid = cell%get_vef_lid_from_gid(vef%get_gid())
      reference_fe => cell%get_reference_fe()
      nodal_iter = reference_fe%create_dofs_on_n_face_iterator(vef_lid)
      ino = 1
      nodal_ids(:) = 0
      do while (.not. nodal_iter%is_upper_bound())
        this%x(inode) = nodal_coords(nodal_iter%get_current())%get(1)
        this%y(inode) = nodal_coords(nodal_iter%get_current())%get(2)
        this%z(inode) = nodal_coords(nodal_iter%get_current())%get(3)
        nodal_ids(ino) = inode
        if (present(conditions)) then
          call conditions%get_function(boundary_id=vef%get_set_id(),component_id=1_ip,function=func%p)
          fun => func%p
          if (associated(fun)) then
            call  fun%get_value(nodal_coords(nodal_iter%get_current()), this%point_data(inode))
          end if
        end if
        inode = inode + 1
        ino = ino + 1
        call nodal_iter%next()
      end do

      if (vef%get_dim() == 0) then
        this%connect(nodal_ids(1:1)) = nodal_ids( nodes_vtk2fempar_0d(:) ) - 1
      else if (vef%get_dim() == 1) then
        this%connect(nodal_ids(1:2)) = nodal_ids( nodes_vtk2fempar_1d(:) ) - 1
      else if (vef%get_dim() == 2) then
        this%connect(nodal_ids(1:4)) = nodal_ids( nodes_vtk2fempar_2d(:) ) - 1
      else
        check(.false.)
      end if

      this%offset(icell)        = inode - 1
      this%cell_type(icell)     = vef_type(vef%get_dim())
      this%cell_data( icell )   = vef%get_set_id()

      icell = icell + 1
      call vef%next()
    end do

    deallocate(nodal_coords,stat=istat); check(istat == 0)
    deallocate(nodal_ids,stat=istat); check(istat == 0)
    call triangulation%free_vef_iterator(vef)
    call triangulation%free_cell_iterator(cell)

  end subroutine uvtkw_attach_vefs

!========================================================================================  
  subroutine  uvtkw_attach_fe_function(this,fe_function,fe_space)

    implicit none
    class(unfitted_vtk_writer_t),   intent(inout) :: this
    class(fe_function_t),                   intent(in)    :: fe_function
    class(serial_fe_space_t),      intent(in)    :: fe_space

    class(triangulation_t), pointer :: triangulation
    class(fe_cell_iterator_t), allocatable  :: fe
    class(reference_fe_t), pointer :: reference_fe_cell
    type(quadrature_t) :: subcell_nodal_quad
    type(interpolation_t) :: fe_cell_interpol, fe_subcell_interpol
    integer(ip) :: ipoint, subcell
    real(rp), allocatable :: cell_nodal_vals(:)
    integer(ip) :: num_subcell_nodes, num_dime, istat, idime, num_subcells, icell
    integer(ip) :: num_cells, num_cell_nodes, num_subelems_x_cell, num_subelems_x_subcell
    integer(ip) :: my_part_id
    integer(ip) :: field_id
    type(tet_lagrangian_reference_fe_t)    :: geo_reference_fe_subcell
    type(tet_lagrangian_reference_fe_t)    :: reference_fe_subcell
    type(quadrature_t), pointer            :: nodal_quadrature_subcell
    type(cell_map_t)                       :: cell_map_subcell
    type(cell_map_t)                       :: cell_map_subcell_ref
    class(reference_fe_t), pointer         :: geo_reference_fe_cell
    type(quadrature_t), pointer            :: nodal_quadrature_cell
    type(cell_map_t)                       :: cell_map_cell
    integer(ip) :: isubelem, jno, ino
    integer(ip), allocatable :: cell_subelems_connectivities(:,:)
    integer(ip), allocatable :: subcell_subelems_connectivities(:,:)
    type(point_t), pointer :: mapped_cell_subelem_coords(:)
    type(point_t), pointer :: mapped_subcell_subelem_coords(:)
    type(point_t), pointer :: mapped_subcell_subelem_coords_ref(:)
    integer(ip) :: num_cell_eval_points, max_num_shape_funs
    integer(ip) :: num_subcell_eval_points
    type(point_t), pointer :: cell_coords(:)
    type(point_t), pointer :: subcell_coords(:)
    real(rp), allocatable :: sol_at_cell_eval_points(:)
    real(rp), allocatable :: sol_at_subcell_eval_points(:)
    type(vector_field_t), allocatable :: vector_sol_at_cell_eval_points(:)
    type(vector_field_t), allocatable :: vector_sol_at_subcell_eval_points(:)
    integer(ip), allocatable :: nodes_vtk2fempar(:)
    integer(ip), allocatable :: nodesids(:)
    integer(ip) :: the_cell_type, the_subcell_type
    real(rp), pointer :: subcell_quad_coords(:,:)

    call this%free()

    triangulation => fe_space%get_triangulation()

    this%environment => triangulation%get_environment()
    if ( .not. this%environment%am_i_l1_task() ) return
    my_part_id = this%environment%get_l1_rank() + 1

    num_dime = triangulation%get_num_dims()

    call fe_space%create_fe_cell_iterator(fe)
    call fe%first_local_non_void(field_id=1)
    assert(.not. fe%has_finished())

    call geo_reference_fe_subcell%create( topology = topology_tet,&
                                      num_dims = num_dime,&
                                      order = 1,&
                                      field_type = field_type_scalar,&
                                      conformity = .true., &
                                      continuity = .false. )

    call reference_fe_subcell%create( topology = topology_tet,&
                                      num_dims = num_dime,&
                                      order = fe%get_max_order_all_fields(),&
                                      field_type = field_type_scalar,&
                                      conformity = .true., &
                                      continuity = .false. )

    ! Count things cell level

    num_cells = triangulation%get_num_local_cells()
    num_cell_nodes = fe%get_num_nodes()
    max_num_shape_funs = fe_space%get_max_num_shape_functions()
    reference_fe_cell => fe%get_max_order_reference_fe()
    num_subelems_x_cell = reference_fe_cell%get_num_subcells(fe%get_max_order_all_fields()-1)
    nodal_quadrature_cell => reference_fe_cell%get_nodal_quadrature()
    num_cell_eval_points = nodal_quadrature_cell%get_num_quadrature_points()

    ! Count things subcell level

    nodal_quadrature_subcell => reference_fe_subcell%get_nodal_quadrature()
    num_subcell_eval_points = nodal_quadrature_subcell%get_num_quadrature_points()
    num_subelems_x_subcell = reference_fe_subcell%get_num_subcells(fe%get_max_order_all_fields()-1)
    select type(triangulation)
      class is(serial_unfitted_triangulation_t)
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

    ! Once we have counted, allocate things
    this%Ne = num_cells*num_subelems_x_cell + num_subcells*num_subelems_x_subcell
    this%Nn = num_cell_nodes*num_cells*num_subelems_x_cell  + num_subcell_nodes*num_subcells*num_subelems_x_subcell

    call memalloc ( this%Nn, this%x, __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%y, __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%z, __FILE__, __LINE__ )
    call memalloc ( this%Ne, this%cell_type, __FILE__, __LINE__ )
    call memalloc ( this%Ne, this%offset   , __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%connect  , __FILE__, __LINE__ )
    call memalloc ( this%Ne, this%cell_data , __FILE__, __LINE__ )
    call memalloc ( this%Ne, this%pid       , __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%point_data, __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%v_x, __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%v_y, __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%v_z, __FILE__, __LINE__ )
    this%v_x(:) = 0.0
    this%v_y(:) = 0.0
    this%v_z(:) = 0.0
    this%point_data(:) = 0.0

    call memalloc ( num_cell_nodes, nodes_vtk2fempar, __FILE__, __LINE__ )
    call memalloc ( num_cell_nodes, nodesids        , __FILE__, __LINE__ )
    call memalloc ( max_num_shape_funs, cell_nodal_vals, __FILE__, __LINE__ )
    call memalloc ( num_cell_eval_points, sol_at_cell_eval_points, __FILE__, __LINE__ )
    call memalloc ( num_subcell_eval_points, sol_at_subcell_eval_points, __FILE__, __LINE__ )
    call memalloc ( num_cell_nodes, num_subelems_x_cell, cell_subelems_connectivities, __FILE__,__LINE__ )
    call memalloc ( num_subcell_nodes, num_subelems_x_subcell, subcell_subelems_connectivities,__FILE__, __LINE__ )

    allocate( vector_sol_at_cell_eval_points(num_cell_eval_points), stat=istat); check(istat == 0)
    allocate( vector_sol_at_subcell_eval_points(num_subcell_eval_points), stat=istat); check(istat == 0)

    ! Setup auxiliary things

    call cell_map_subcell%create( nodal_quadrature_subcell, geo_reference_fe_subcell )
    call cell_map_subcell_ref%create( nodal_quadrature_subcell, geo_reference_fe_subcell )
    call reference_fe_subcell%get_subcells_connectivity(fe%get_max_order_all_fields()-1, subcell_subelems_connectivities)

    geo_reference_fe_cell => fe%get_reference_fe_geo()
    reference_fe_cell => fe%get_max_order_reference_fe()

    call cell_map_cell%create( nodal_quadrature_cell, geo_reference_fe_cell )
    call reference_fe_cell%get_subcells_connectivity(fe%get_max_order_all_fields()-1, cell_subelems_connectivities)

    call subcell_nodal_quad%create(num_dime,num_cell_eval_points)

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

    ! Loops to fill data

    ipoint = 1
    icell = 1
    call fe%first()
    do while ( .not. fe%has_finished() )

       ! Skip ghost elems
       if (fe%is_ghost()) then
         call fe%next(); cycle
       end if

       call fe%update_sub_triangulation()

       sol_at_cell_eval_points(:) = 0.0
       do ino = 1, num_cell_eval_points
         call vector_sol_at_cell_eval_points(ino)%init(0.0)
       end do
       if (.not. fe%is_exterior()) then
         do field_id = 1, fe%get_num_fields()
           call fe_function%gather_nodal_values(fe,field_id,cell_nodal_vals)
           reference_fe_cell => fe%get_reference_fe(field_id)
           call reference_fe_cell%create_interpolation(nodal_quadrature_cell,fe_cell_interpol)
           if (fe%get_field_type(field_id)==field_type_scalar) then
             call reference_fe_cell%evaluate_fe_function_scalar(fe_cell_interpol,cell_nodal_vals,sol_at_cell_eval_points)
           else if (fe%get_field_type(field_id)==field_type_vector) then
             call reference_fe_cell%evaluate_fe_function_vector(fe_cell_interpol,cell_nodal_vals,vector_sol_at_cell_eval_points)
           else
             check(.false.)
           end if
         end do
       end if

       cell_coords => cell_map_cell%get_coordinates()
       call fe%get_nodes_coordinates( cell_coords )
       call cell_map_cell%update(nodal_quadrature_cell,no_ressemblance)
       mapped_cell_subelem_coords => cell_map_cell%get_quadrature_points_coordinates()

       do isubelem = 1, num_subelems_x_cell
         do ino = 1, num_cell_nodes
           jno = cell_subelems_connectivities(ino,isubelem)
           this%x(ipoint) = mapped_cell_subelem_coords(jno)%get(1)
           this%y(ipoint) = mapped_cell_subelem_coords(jno)%get(2)
           this%z(ipoint) = mapped_cell_subelem_coords(jno)%get(3)
           nodesids(ino) = ipoint
           do field_id = 1, fe%get_num_fields()
             if (fe%get_field_type(field_id)==field_type_scalar) then
               this%point_data(ipoint) = sol_at_cell_eval_points(jno)
             else if (fe%get_field_type(field_id)==field_type_vector) then
               this%v_x(ipoint) = vector_sol_at_cell_eval_points(jno)%get(1)
               this%v_y(ipoint) = vector_sol_at_cell_eval_points(jno)%get(2)
               this%v_z(ipoint) = vector_sol_at_cell_eval_points(jno)%get(3)
             else
               check(.false.)
             end if
           end do
           ipoint = ipoint + 1
         end do
         this%connect(nodesids(:)) = nodesids( nodes_vtk2fempar(:) ) - 1

         this%offset(icell)        = ipoint - 1
         this%cell_type(icell)     = the_cell_type
         if ( fe%is_interior() ) then
           this%cell_data( icell )  = -1
         else if ( fe%is_cut() ) then
           this%cell_data( icell )  = 0
         else
           this%cell_data( icell )  = 1
         end if
         this%pid(icell) = real(my_part_id,kind=rp)
         icell = icell + 1
       end do

       do subcell = 1, fe%get_num_subcells()

         subcell_coords => cell_map_subcell%get_coordinates()
         call fe%get_phys_coords_of_subcell(subcell,subcell_coords)
         call cell_map_subcell%update(nodal_quadrature_subcell,no_ressemblance)
         mapped_subcell_subelem_coords => cell_map_subcell%get_quadrature_points_coordinates()

         sol_at_subcell_eval_points(:) = 0.0
         do ino = 1, num_subcell_eval_points
           call vector_sol_at_subcell_eval_points(ino)%init(0.0)
         end do

         if (fe%is_interior_subcell(subcell)) then
           subcell_coords => cell_map_subcell_ref%get_coordinates()
           call fe%get_ref_coords_of_subcell(subcell,subcell_coords)
           call cell_map_subcell_ref%update(nodal_quadrature_subcell,no_ressemblance)
           mapped_subcell_subelem_coords_ref => cell_map_subcell_ref%get_quadrature_points_coordinates()
           subcell_quad_coords => subcell_nodal_quad%get_coordinates()
           do ino = 1, num_subcell_eval_points
             do idime = 1, num_dime
               subcell_quad_coords(idime,ino) = mapped_subcell_subelem_coords_ref(ino)%get(idime)
             end do
           end do

           do field_id = 1, fe%get_num_fields()
             reference_fe_cell => fe%get_reference_fe(field_id)
             call reference_fe_cell%create_interpolation(subcell_nodal_quad,fe_subcell_interpol)
             call fe_function%gather_nodal_values(fe,field_id,cell_nodal_vals)
             if (fe%get_field_type(field_id)==field_type_scalar) then
               call reference_fe_cell%evaluate_fe_function_scalar(fe_subcell_interpol,cell_nodal_vals,sol_at_subcell_eval_points)
             else if (fe%get_field_type(field_id)==field_type_vector) then
               call reference_fe_cell%evaluate_fe_function_vector(fe_subcell_interpol,cell_nodal_vals,vector_sol_at_subcell_eval_points)
             else
               check(.false.)
             end if
           end do

         end if

         do isubelem = 1, num_subelems_x_subcell
           do ino = 1, num_subcell_nodes
             jno = subcell_subelems_connectivities(ino,isubelem)
             this%x(ipoint) = mapped_subcell_subelem_coords(jno)%get(1)
             this%y(ipoint) = mapped_subcell_subelem_coords(jno)%get(2)
             this%z(ipoint) = mapped_subcell_subelem_coords(jno)%get(3)

             do field_id = 1, fe%get_num_fields()
               if (fe%get_field_type(field_id)==field_type_scalar) then
                 this%point_data(ipoint) = sol_at_subcell_eval_points(jno)
               else if (fe%get_field_type(field_id)==field_type_vector) then
                 this%v_x(ipoint) = vector_sol_at_subcell_eval_points(jno)%get(1)
                 this%v_y(ipoint) = vector_sol_at_subcell_eval_points(jno)%get(2)
                 this%v_z(ipoint) = vector_sol_at_subcell_eval_points(jno)%get(3)
               else
                 check(.false.)
               end if
             end do

             this%connect(ipoint) = ipoint - 1
             ipoint = ipoint + 1
           end do
           this%offset(icell)        = ipoint - 1
           this%cell_type(icell)     = the_subcell_type
           if ( fe%is_interior_subcell(subcell) ) then
             this%cell_data( icell )  = -2
           else
             this%cell_data( icell )  = 2
           end if
           this%pid(icell) = real(my_part_id,kind=rp)
           icell = icell + 1
         end do

       end do

       call fe%next()
    end do

    call memfree ( cell_subelems_connectivities, __FILE__,__LINE__ )
    call memfree ( subcell_subelems_connectivities,__FILE__, __LINE__ )
    call memfree ( sol_at_cell_eval_points,         __FILE__, __LINE__ )
    call memfree ( sol_at_subcell_eval_points, __FILE__, __LINE__ )
    call memfree ( cell_nodal_vals, __FILE__, __LINE__ )
    deallocate( vector_sol_at_cell_eval_points, stat=istat); check(istat == 0)
    deallocate( vector_sol_at_subcell_eval_points, stat=istat); check(istat == 0)
    call subcell_nodal_quad%free()
    call fe_cell_interpol%free()
    call fe_subcell_interpol%free()
    call geo_reference_fe_subcell%free()
    call reference_fe_subcell%free()
    call cell_map_subcell%free()
    call cell_map_subcell_ref%free()
    call cell_map_cell%free()
    call fe_space%free_fe_cell_iterator(fe)
    call memfree ( nodesids, __FILE__, __LINE__ )
    call memfree ( nodes_vtk2fempar, __FILE__, __LINE__ )

  end subroutine  uvtkw_attach_fe_function

!========================================================================================
  subroutine uvtkw_attach_boundary_quadrature_points( this, fe_space )
  
    implicit none
    class(unfitted_vtk_writer_t),   intent(inout) :: this
    class(serial_fe_space_t)    ,   intent(in)    :: fe_space
  
    class(fe_cell_iterator_t),allocatable :: fe
    type(quadrature_t), pointer :: quadrature
    type(point_t), pointer :: quadrature_points_coordinates(:)
    type(piecewise_cell_map_t),     pointer :: cell_map
    integer(ip) :: num_dime
    integer(ip) :: qpoint, num_quad_points, total_num_quad_points
    integer(ip) :: ipoint
    type(vector_field_t)  :: normal_vec
    integer(ip), parameter :: the_point_type = 1

    class(triangulation_t), pointer :: triangulation

    call this%free()

    triangulation => fe_space%get_triangulation()

    this%environment => triangulation%get_environment()
    if ( .not. this%environment%am_i_l1_task() ) return

    num_dime       = triangulation%get_num_dims()
  
    call fe_space%create_fe_cell_iterator(fe)

    total_num_quad_points = 0
    do while ( .not. fe%has_finished() )
       call fe%update_boundary_integration()
       quadrature => fe%get_boundary_quadrature()
       num_quad_points = quadrature%get_num_quadrature_points()
       total_num_quad_points = total_num_quad_points + num_quad_points
       call fe%next()
    end do
  
    this%Ne = total_num_quad_points
    this%Nn = total_num_quad_points
  
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
    class(serial_fe_space_t),      intent(in) :: fe_space
  
    class(fe_facet_iterator_t),allocatable :: fe_facet
    type(quadrature_t), pointer :: quadrature
    type(point_t), pointer :: quadrature_points_coordinates(:)
    integer(ip) :: num_dime
    integer(ip) :: qpoint, num_quad_points
    integer(ip) :: ipoint
    type(vector_field_t)  :: normal_vec(2)
    integer(ip), parameter :: the_point_type = 1
    real(rp) :: mystatus
    real(rp) :: dS
    type(vector_field_t), allocatable  :: shape_gradients(:,:)
    real(rp)            , allocatable  :: shape_values(:,:)
    integer(ip) :: num_cells_arround, ineigh, istat

    class(triangulation_t), pointer :: triangulation

    call this%free()

    triangulation => fe_space%get_triangulation()

    this%environment => triangulation%get_environment()
    if ( .not. this%environment%am_i_l1_task() ) return

    num_dime       = triangulation%get_num_dims()

    !call fe_space%serial_fe_space_t%set_up_facet_integration()
  
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

       do ineigh = 1, fe_facet%get_num_cells_around()
         call fe_facet%get_gradients(ineigh,shape_gradients,1)
         call fe_facet%get_values   (ineigh,shape_values,1)
       end do
  
       do qpoint = 1, num_quad_points
         dS =  fe_facet%get_det_jacobian(qpoint)*quadrature%get_weight(qpoint)
         this%x(ipoint) = quadrature_points_coordinates(qpoint)%get(1)
         this%y(ipoint) = quadrature_points_coordinates(qpoint)%get(2)
         this%z(ipoint) = quadrature_points_coordinates(qpoint)%get(3)
         call fe_facet%get_normal(qpoint,normal_vec)
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
    call memfree(shape_values, __FILE__, __LINE__)
    deallocate (shape_gradients, stat=istat); check(istat==0);
  
  end subroutine uvtkw_attach_facets_quadrature_points

!========================================================================================  
  subroutine uvtkw_write_to_vtk_file(this,filename)

    implicit none
    class(unfitted_vtk_writer_t), intent(in) :: this
    character(*), intent(in) :: filename

    integer(ip) :: E_IO
    integer(ip) :: iounit

    assert(allocated(this%x        ))
    assert(allocated(this%y        ))
    assert(allocated(this%z        ))
    assert(allocated(this%connect  ))
    assert(allocated(this%offset   ))
    assert(allocated(this%cell_type))

    assert(associated(this%environment))
    assert(this%environment%am_i_l1_task())


    if (this%Nn == 0 .or. this%Ne == 0) then
      iounit = io_open(file=filename,action='write')
      check(iounit>0)
      write(iounit,*) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">                                 '
      write(iounit,*) '  <UnstructuredGrid>                                                                                      '
      write(iounit,*) '    <Piece NumberOfPoints="+0" NumberOfCells="+0">                                                        '
      write(iounit,*) '      <Points>                                                                                            '
      write(iounit,*) '        <DataArray type="Float64" NumberOfComponents="3" Name="Points" format="ascii"></DataArray>        '
      write(iounit,*) '      </Points>                                                                                           '
      write(iounit,*) '      <Cells>                                                                                             '
      write(iounit,*) '        <DataArray type="Int32" Name="connectivity" format="ascii"></DataArray>                           '
      write(iounit,*) '        <DataArray type="Int32" Name="offsets" format="ascii"></DataArray>                                '
      write(iounit,*) '        <DataArray type="Int8" Name="types" format="ascii"></DataArray>                                   '
      write(iounit,*) '      </Cells>                                                                                            '
      write(iounit,*) '      <CellData>                                                                                          '
      write(iounit,*) '        <DataArray type="Float64" Name="cell_scalars" NumberOfComponents="1" format="ascii"></DataArray>  '
      write(iounit,*) '        <DataArray type="Float64" Name="part_id" NumberOfComponents="1" format="ascii"> </DataArray>      '
      write(iounit,*) '      </CellData>                                                                                         '
      write(iounit,*) '      <PointData>                                                                                         '
      write(iounit,*) '        <DataArray type="Float64" Name="point_scalars" NumberOfComponents="1" format="ascii"></DataArray> '
      write(iounit,*) '        <DataArray type="Float64" Name="point_vectors" NumberOfComponents="3" format="ascii"></DataArray> '
      write(iounit,*) '      </PointData>                                                                                        '
      write(iounit,*) '    </Piece>                                                                                              '
      write(iounit,*) '  </UnstructuredGrid>                                                                                     '
      write(iounit,*) '</VTKFile>                                                                                                '
      call io_close(iounit)
      return
    end if

    E_IO = VTK_INI_XML(output_format = 'binary', filename = filename, mesh_topology = 'UnstructuredGrid')
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


    E_IO = VTK_DAT_XML(var_location = 'node', var_block_action = 'opeN')
    if (allocated(this%point_data)) then
      E_IO = VTK_VAR_XML(NC_NN = this%Nn, varname = 'point_scalars', var = this%point_data)
    end if

    if (allocated(this%v_x)) then
      E_IO = VTK_VAR_XML(NC_NN = this%Nn, varname = 'point_vectors', varX = this%v_x, varY = this%v_y, varZ = this%v_z )
    end if
    E_IO = VTK_DAT_XML(var_location = 'node', var_block_action = 'CLOSE')

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


