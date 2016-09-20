#include <stdio.h>
#include <assert.h>

#ifdef ENABLE_P4EST

#ifndef P4_TO_P8
#include <p4est_vtk.h>
#include <p4est_bits.h>
#include <p4est_iterate.h>
#include <p4est_lnodes.h>
#include <p4est_ghost.h>
#else
#include <p8est_vtk.h>
#include <p8est_bits.h>
#include <p8est_iterate.h>
#include <p8est_lnodes.h>
#include <p8est_ghost.h>
#endif

// todo: probably provide 2 different variants of p4estw_get_data for 2 and 3 dimmensions
#define NUM_DIMS 2

#define NUM_ELEM_NODES  4
#define NUM_ELEM_FACES  4

#define NUM_NEIGHBOURING_VERTICES 2
#define DOES_NOT_EXIST -1

#define IS_NOT_CORNER  -1
#define UNNUMBERED_CORNER  -2


#define REF_FLAG_NOTHING 0
#define REF_FLAG_REFINE 1
#define REF_FLAG_COARSEN -1
#define REF_FLAG_INVALID -11111

// element data
// todo: rewise
typedef struct element_data
{
  int global_idx;
  int local_idx;

  // should be changed to enum. It is possible to use C enums in Fortran...
  int refinement_flag;
}
element_data_t;


typedef struct p4est_wrapper
{
  p4est_t            *p4est;
  p4est_ghost_t      *ghost;
  element_data_t     *ghost_data;
  p4est_lnodes_t     *lnodes;
  int                mpi_rank, mpi_num_proc;
  int                info_found;
  int                global_elem_idx_offset;

  int num_nodes, num_full_nodes;
  int num_elements;

}
p4est_wrapper_t;

p4est_wrapper_t* p4estw_create_unit_cube(int fortran_mpicomm);
p4est_wrapper_t* p4estw_create_from_geometry(int num_elements, int num_points, int* list_nodes, double coordinates[][NUM_DIMS], int fortran_mpicomm);

int p4estw_get_sizes(p4est_wrapper_t *p4estw, int *num_elements, int *num_trees, int64_t *first_global_idx);
int p4estw_get_quadrants(p4est_wrapper_t *p4estw, int *tree_indices, int *tree_offsets, int64_t *quadrant_morton_ids);
int p4estw_get_ghost_sizes(p4est_wrapper_t *p4estw, int *num_ghosts, int *num_trees, int *num_processors);
int p4estw_get_ghosts(p4est_wrapper_t *p4estw, int *tree_indices, int *tree_offsets, int *proc_indices, int *proc_offsets, int *quadrant_data, int64_t *global_indices);

int p4estw_refine_all(p4est_wrapper_t* p4estw);

#endif
