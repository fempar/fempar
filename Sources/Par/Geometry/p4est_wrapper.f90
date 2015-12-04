module p4est_wrapper
   
   interface
      type(C_PTR) function p4estw_create_unit_cube(fortran_comm) bind(c, name='p4estw_create_unit_cube')
         use iso_c_binding
         implicit none
         integer(C_INT), value, intent(in)  :: fortran_comm !todo: what kind of integer should be specified here?
      end function p4estw_create_unit_cube
      
      !p4est_wrapper_t* p4estw_create_from_geometry(int num_elements, int num_points, int* list_nodes, double coordinates[][NUM_DIMS], sc_MPI_Comm mpicomm);
      type(C_PTR) function p4estw_create_from_geometry(num_elements, num_points, list_nodes, coordinates, fortran_comm) bind(c, name='p4estw_create_from_geometry')
         use iso_c_binding
         implicit none
         integer(C_INT), value, intent(in)  :: num_elements, num_points
         type(C_PTR), value, intent(in) :: list_nodes, coordinates         
         integer(C_INT), value, intent(in)  :: fortran_comm !todo: what kind of integer should be specified here?
      end function p4estw_create_from_geometry
      
      integer(C_INT) function p4estw_refine_all(p4est) bind(c, name='p4estw_refine_all')
         use iso_c_binding
         implicit none
         type(C_PTR), value, intent(in) :: p4est
      end function p4estw_refine_all

      integer(C_INT) function p4estw_refine_selected(p4est, what_refine_size, what_refine) bind(c, name='p4estw_refine_selected')
         use iso_c_binding
         implicit none
         type(C_PTR), value, intent(in) :: p4est
         integer(C_INT), value, intent(in) :: what_refine_size
         type(C_PTR), value, intent(in) :: what_refine
      end function p4estw_refine_selected

      integer(C_INT) function p4estw_partition(p4est) bind(c, name='p4estw_partition')
         use iso_c_binding
         implicit none
         type(C_PTR), value, intent(in) :: p4est
      end function p4estw_partition

      integer(C_INT) function p4estw_destroy(p4est) bind(c, name='p4estw_destroy')
         use iso_c_binding
         implicit none
         type(C_PTR), value, intent(in) :: p4est
      end function p4estw_destroy

      
      integer(C_INT) function p4estw_get_sizes(p4est, num_elements, num_trees, first_gloabl_idx) bind(c, name='p4estw_get_sizes')
         use iso_c_binding
         implicit none
         type(C_PTR), value, intent(in) :: p4est
         integer(C_INT), intent(out) :: num_elements, num_trees
         integer(C_LONG), intent(out) :: first_gloabl_idx
      end function p4estw_get_sizes

      integer(C_INT) function p4estw_get_quadrants(p4est, tree_indices, tree_offsets, quadrant_data) bind(c, name='p4estw_get_quadrants')
         use iso_c_binding
         implicit none
         type(C_PTR), value, intent(in) :: p4est
         type(C_PTR), value :: tree_indices
         type(C_PTR), value :: tree_offsets
         type(C_PTR), value :: quadrant_data
      end function p4estw_get_quadrants
           
!       integer(C_INT) function p4estw_enumerate_quadrants(p4est) bind(c, name='p4estw_enumerate_quadrants')
!          use iso_c_binding
!          implicit none
!          type(C_PTR), value, intent(in) :: p4est
!       end function p4estw_enumerate_quadrants
           
      integer(C_INT) function p4estw_get_ghost_sizes(p4est, num_ghosts, num_ghost_trees, num_ghost_processors) bind(c, name='p4estw_get_ghost_sizes')
         use iso_c_binding
         implicit none
         type(C_PTR), value, intent(in) :: p4est
         integer(C_INT), intent(out) :: num_ghosts, num_ghost_trees, num_ghost_processors
      end function p4estw_get_ghost_sizes
      
      integer(C_INT) function p4estw_get_ghosts(p4est, tree_indices, tree_offsets, processor_indices, processor_offsets, quadrant_data, global_indices) bind(c, name='p4estw_get_ghosts')
         use iso_c_binding
         implicit none
         type(C_PTR), value, intent(in) :: p4est
         type(C_PTR), value :: tree_indices
         type(C_PTR), value :: tree_offsets
         type(C_PTR), value :: processor_indices
         type(C_PTR), value :: processor_offsets
         type(C_PTR), value :: quadrant_data
         type(C_PTR), value :: global_indices
      end function p4estw_get_ghosts
           

      
   end interface

   public :: p4estw_create_from_geometry, p4estw_refine_all, p4estw_refine_selected
   public :: p4estw_get_sizes, p4estw_get_quadrants, p4estw_get_ghost_sizes, p4estw_get_ghosts
   contains
   
!=============================================================================
#ifndef ENABLE_P4EST
  subroutine enable_p4est_error_message
      implicit none
      write (0,*) 'Error: Fempar was not compiled with -DENABLE_P4EST.'
      write (0,*) 'Error: You must activate this cpp macro in order to use p4est'
      stop
  end subroutine
#endif

   
end module p4est_wrapper
      
