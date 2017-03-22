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

module poisson_unfitted_vtk_writer_names

  use fempar_names
  use serial_unfitted_triangulation_names
  use serial_unfitted_fe_space_names
  use IR_Precision ! VTK_IO
  use Lib_VTK_IO ! VTK_IO

  implicit none
# include "debug.i90"
  private

  type :: poisson_unfitted_vtk_writer_t
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

  contains

    procedure, non_overridable :: attach_triangulation => puvtk_attach_triangulation
    procedure, non_overridable :: attach_fe_function   => puvtk_attach_fe_function
    procedure, non_overridable :: write_to_vtk_file    => puvtk_write_to_vtk_file
    procedure, non_overridable :: free                 => puvtk_free

  end type poisson_unfitted_vtk_writer_t

  public :: poisson_unfitted_vtk_writer_t

contains

  !------------------------------------------------------------------
  subroutine  puvtk_attach_triangulation(this,triangulation)

    implicit none
    class(poisson_unfitted_vtk_writer_t),   intent(inout) :: this
    class(serial_unfitted_triangulation_t), intent(in)    :: triangulation

    integer(ip) :: num_cells, num_cell_nodes, num_subcells, num_subcell_nodes, num_dime
    integer(ip) :: istat, icell, inode, ino, isubcell
    type(unfitted_cell_iterator_t)  :: cell_iter
    type(unfitted_cell_accessor_t)  :: cell
    type(point_t), allocatable, dimension(:) :: cell_coords, subcell_coords
    integer(ip) :: the_cell_type, the_subcell_type
    integer(ip), allocatable :: nodes_vtk2fempar(:), nodesids(:)
    
    call this%free()

    cell_iter = triangulation%create_unfitted_cell_iterator()
    call cell_iter%current(cell)

    num_dime = triangulation%get_num_dimensions()
    num_cells = triangulation%get_num_cells()
    num_cell_nodes = cell%get_num_nodes()
    num_subcells = triangulation%get_total_num_of_subcells()
    num_subcell_nodes = triangulation%get_max_num_nodes_in_subcell()
    this%Ne = num_cells + num_subcells
    this%Nn = num_cell_nodes*num_cells  + num_subcell_nodes*num_subcells

    call memalloc ( this%Nn, this%x, __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%y, __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%z, __FILE__, __LINE__ )
    call memalloc ( this%Ne, this%cell_type, __FILE__, __LINE__ )
    call memalloc ( this%Ne, this%offset   , __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%connect  , __FILE__, __LINE__ )
    call memalloc ( this%Ne, this%cell_data       , __FILE__, __LINE__ )
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
    call cell_iter%init()
    icell = 1
    inode = 1
    do while ( .not. cell_iter%has_finished() )

      call cell_iter%current(cell)

      call cell%update_sub_triangulation()

      call cell%get_coordinates( cell_coords )

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

      icell = icell + 1

      do isubcell = 1, cell%get_number_of_subcells()
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

        icell = icell + 1

      end do

      call cell_iter%next()
    end do


    if (num_dime == 2_ip) this%z(:) = 0


    deallocate ( cell_coords, stat = istat ); check(istat == 0)
    deallocate ( subcell_coords, stat = istat ); check(istat == 0)
    call memfree ( nodes_vtk2fempar, __FILE__, __LINE__ )
    call memfree ( nodesids, __FILE__, __LINE__ )

  end subroutine  puvtk_attach_triangulation

  !------------------------------------------------------------------
  subroutine  puvtk_attach_fe_function(this,fe_function,fe_space)

    implicit none
    class(poisson_unfitted_vtk_writer_t),   intent(inout) :: this
    class(fe_function_t),                   intent(in)    :: fe_function
    class(serial_unfitted_fe_space_t),      intent(in)    :: fe_space

    class(base_static_triangulation_t), pointer :: triangulation
    type(unfitted_fe_iterator_t) :: fe_iterator
    type(unfitted_fe_accessor_t), target :: fe
    class(fe_accessor_t), pointer :: fe_ptr
    class(unfitted_cell_accessor_t), pointer :: cell
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
    select type(triangulation)
      class is(serial_unfitted_triangulation_t)
        call this%attach_triangulation(triangulation)
        num_subelem_nodes = triangulation%get_max_num_nodes_in_subcell()
      class default
        check(.false.)
    end select

    num_elem_nodes = triangulation%get_max_number_shape_functions()
    num_dime = triangulation%get_num_dimensions()

    call memalloc ( num_elem_nodes, nodal_vals, __FILE__, __LINE__ )
    call memalloc ( num_subelem_nodes, subelem_nodal_vals, __FILE__, __LINE__ )
    call memalloc ( this%Nn, this%point_data, __FILE__, __LINE__ )
    allocate(subcell_points(1:num_subelem_nodes),stat = istat); check(istat == 0)

    call subcel_nodal_quad%create(num_dime,num_subelem_nodes)
    subcell_coords => subcel_nodal_quad%get_coordinates()

    ipoint = 1
    fe_iterator = fe_space%create_unfitted_fe_iterator()
    do while ( .not. fe_iterator%has_finished() )

       call fe_iterator%current(fe)
       fe_ptr => fe
       cell => fe%get_unfitted_cell_accessor()
       ref_fe => fe%get_reference_fe(1) ! TODO we assume a single field

       ! Recover nodal values
       if (cell%is_exterior()) then
         nodal_vals(:) = 0.0
       else
         call fe_function%gather_nodal_values(fe_ptr,1,nodal_vals)!TODO we assume a single field
       end if

       do inode = 1,num_elem_nodes
         this%point_data(ipoint) = nodal_vals(inode)
         ipoint = ipoint + 1
       end do

       call cell%update_sub_triangulation()

       do subcell = 1, cell%get_number_of_subcells()

         ! Get the subcell values
         if (cell%is_interior_subcell(subcell)) then
           call cell%get_ref_coords_of_subcell(subcell,subcell_points)
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

       call fe_iterator%next()
    end do

    call memfree ( nodal_vals,         __FILE__, __LINE__ )
    call memfree ( subelem_nodal_vals, __FILE__, __LINE__ )
    deallocate(subcell_points,stat = istat); check(istat == 0)
    call subcel_nodal_quad%free()
    call fe_interpol%free()


  end subroutine  puvtk_attach_fe_function

  !------------------------------------------------------------------
  subroutine puvtk_write_to_vtk_file(this,filename)

    implicit none
    class(poisson_unfitted_vtk_writer_t), intent(in) :: this
    character(*), intent(in) :: filename

    integer(ip) :: E_IO

    assert(allocated(this%x        ))
    assert(allocated(this%y        ))
    assert(allocated(this%z        ))
    assert(allocated(this%connect  ))
    assert(allocated(this%offset   ))
    assert(allocated(this%cell_type))

    E_IO = VTK_INI_XML(output_format = 'ascii', filename = filename, mesh_topology = 'UnstructuredGrid')
    E_IO = VTK_GEO_XML(NN = this%Nn, NC = this%Ne, X = this%x, Y = this%y, Z = this%z)
    E_IO = VTK_CON_XML(NC = this%Ne, connect = this%connect, offset = this%offset, cell_type = int(this%cell_type,I1P) )

    if (allocated(this%cell_data)) then
      E_IO = VTK_DAT_XML(var_location = 'cell', var_block_action = 'opeN')
      E_IO = VTK_VAR_XML(NC_NN = this%Ne, varname = 'cell_scalars', var = this%cell_data)
      E_IO = VTK_DAT_XML(var_location = 'cell', var_block_action = 'CLOSE')
    end if

    if (allocated(this%point_data)) then
      E_IO = VTK_DAT_XML(var_location = 'node', var_block_action = 'opeN')
      E_IO = VTK_VAR_XML(NC_NN = this%Nn, varname = 'point_scalars', var = this%point_data)
      E_IO = VTK_DAT_XML(var_location = 'node', var_block_action = 'CLOSE')
    end if

    E_IO = VTK_GEO_XML()
    E_IO = VTK_END_XML()

  end subroutine puvtk_write_to_vtk_file

  !------------------------------------------------------------------
  subroutine puvtk_free(this)

    implicit none
    class(poisson_unfitted_vtk_writer_t), intent(inout) :: this

    if (allocated(this%x         )) call memfree ( this%x          , __FILE__, __LINE__ )
    if (allocated(this%y         )) call memfree ( this%y          , __FILE__, __LINE__ )
    if (allocated(this%z         )) call memfree ( this%z          , __FILE__, __LINE__ )
    if (allocated(this%cell_type )) call memfree ( this%cell_type  , __FILE__, __LINE__ )
    if (allocated(this%offset    )) call memfree ( this%offset     , __FILE__, __LINE__ )
    if (allocated(this%connect   )) call memfree ( this%connect    , __FILE__, __LINE__ )
    if (allocated(this%cell_data )) call memfree ( this%cell_data  , __FILE__, __LINE__ )
    if (allocated(this%point_data)) call memfree ( this%point_data , __FILE__, __LINE__ )

  end subroutine puvtk_free

end module poisson_unfitted_vtk_writer_names
!***************************************************************************************************


