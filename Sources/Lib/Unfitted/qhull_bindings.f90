module qhull_bindings_names

  use types_names
  use, intrinsic :: iso_c_binding, only: c_int, c_double

  implicit none

  integer(ip), parameter :: qh_coordT = c_double
  integer(ip), parameter :: qh_int    = c_int

#ifdef ENABLE_QHULL
  interface

    function delaunay_init_and_compute_f90(num_dims, num_points, points, num_cells) bind(C,name="delaunay_init_and_compute")
      import :: qh_int, qh_coordT
      integer(qh_int), value, intent(in)    :: num_dims
      integer(qh_int), value, intent(in)    :: num_points
      real(qh_coordT)       , intent(in)    :: points(num_dims,*)
      integer(qh_int)       , intent(inout) :: num_cells
      integer(qh_int) :: delaunay_init_and_compute_f90
    end function delaunay_init_and_compute_f90

    function delaunay_fill_cells_f90(num_dims, num_cells, cells) bind(C,name="delaunay_fill_cells")
      import :: qh_int
      integer(qh_int), value, intent(in)    :: num_dims
      integer(qh_int), value, intent(in)    :: num_cells
      integer(qh_int)       , intent(inout) :: cells(num_dims+1,*)
      integer(qh_int) :: delaunay_fill_cells_f90
    end function delaunay_fill_cells_f90

    function delaunay_free_f90() bind(C,name="delaunay_free")
      import :: qh_int
      integer(qh_int) :: delaunay_free_f90
    end function delaunay_free_f90

  end interface
#endif

end module qhull_bindings_names

