! at some point, I will also have to provide the following information:
! for each subdomain, list of neighbouring subdomains and the list of 
! neighbouring elements GIDs. To be used in mesh_distribution.

program create_adapted_mesh
  use serial_names
  use par_names
  use p4est_wrapper
  use migratory_element_names
  !use JP_element_topology_names
  
  use, intrinsic :: ISO_C_BINDING  
  
  implicit none
#include "debug.i90"

  type(par_context_t) :: p_context
  type(c_ptr)         :: p4est
  integer(ip)         :: p4est_error
  integer(ip)         :: num_elements, num_trees, num_ghosts, num_ghost_trees, num_ghost_processors
  integer(ip)         :: elem_idx
  integer(igp)        :: first_global_idx
  
  integer(ip), allocatable, target :: tree_indices(:), tree_offsets(:) 
  integer(igp), allocatable, target :: quadrant_data(:) 
  
  type(cell_t)              :: element
  type(hash_migratory_element_set_t)       :: element_set
  type(hash_migratory_element_iterator_t)  :: element_iterator


  call par_context_create(p_context)
  p4est = p4estw_create_unit_cube(p_context%icontxt)

  p4est_error = p4estw_refine_all(p4est)
  p4est_error = p4estw_refine_all(p4est)
  p4est_error = p4estw_partition(p4est)
  
  p4est_error = p4estw_get_sizes(p4est, num_elements, num_trees, first_global_idx)
  call memalloc (num_trees, tree_indices, __FILE__, __LINE__)
  call memalloc (num_trees, tree_offsets, __FILE__, __LINE__)
  call memalloc (num_elements, quadrant_data, __FILE__, __LINE__)
  p4est_error = p4estw_get_quadrants(p4est, C_LOC(tree_indices), C_LOC(tree_offsets), C_LOC(quadrant_data))
 
  do elem_idx = 1, num_elements
    
  end do


  p4est_error = p4estw_get_ghost_sizes(p4est, num_ghosts, num_ghost_trees, num_ghost_processors)
  
  !call element_set%create(num_elements, element)
  
  
  write(*,*) "rank ", p_context%iam, "num_elem ", num_elements, "first_glob_id ", first_global_idx, "num ghosts", num_ghosts, "num gh proc ", num_ghost_processors
  write(*,*) "rank ", p_context%iam, "quads ", quadrant_data

  call par_context_free ( p_context )
  
end program create_adapted_mesh

