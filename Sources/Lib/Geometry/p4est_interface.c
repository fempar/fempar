#include "p4est_interface.h"

#ifdef ENABLE_P4EST

void create_ghost_data(p4est_wrapper_t* p4estw);
void destroy_additional_data(p4est_wrapper_t* p4estw);


//*******************************************************************************************
// Retrieving information
//*******************************************************************************************

int p4estw_get_sizes(p4est_wrapper_t *p4estw, int *num_elements, int *num_trees, int64_t *first_global_idx)
{
  *num_elements = p4estw->p4est->local_num_quadrants;
  *num_trees = p4estw->p4est->last_local_tree - p4estw->p4est->first_local_tree + 1;
  *first_global_idx = p4estw->p4est->global_first_quadrant[p4estw->p4est->mpirank];
}

// how many integers per quad. For 2D it is 2 coordinates and level -> 3
// todo: maybe it should be 4 as in 3D
// todo: or should we pass global idx as well?
#define DATA_PER_QUAD  3

// this has been used previously
// we obtained x, y and level
int p4estw_get_quadrants_old(p4est_wrapper_t *p4estw, int *tree_indices, int *tree_offsets, int *quadrant_data)
{
  int tt, q, Q,k, local_tree_num;
  p4est_tree_t       *tree;     /* Pointer to one octree */
  p4est_quadrant_t *quad;
  sc_array_t         *tquadrants;       /* Quadrant array for one tree */
  element_data_t     *elem_data;


  //printf("tree range %d ---- %d, num elements %d\n", p4estw->p4est->first_local_tree, p4estw->p4est->last_local_tree, p4estw->p4est->local_num_quadrants);
  for (tt = p4estw->p4est->first_local_tree, k = 0, local_tree_num = 0;
       tt <= p4estw->p4est->last_local_tree; ++tt, ++local_tree_num)
  {
//    printf(" treee %d\n", tt);
    tree_indices[local_tree_num] = tt + 1;  // enumerate from 1 for Fortran
    tree_offsets[local_tree_num] = k + 1;   // enumerate from 1 for Fortran
    tree = p4est_tree_array_index (p4estw->p4est->trees, tt);   /* Current tree */
    tquadrants = &tree->quadrants;
    Q = (p4est_locidx_t) tquadrants->elem_count;
    for (q = 0; q < Q; ++q, ++k) {
      quad = p4est_quadrant_array_index (tquadrants, q);
//      printf("quadrant %d, x %d, y %d, l %d\n", k, quad->x, quad->y, quad->level);
      quadrant_data[DATA_PER_QUAD * k + 0] = quad->x;
      quadrant_data[DATA_PER_QUAD * k + 1] = quad->y;
      quadrant_data[DATA_PER_QUAD * k + 2] = quad->level;
      elem_data = (element_data_t* ) quad->p.user_data;
      elem_data->local_idx = k;   // enumerate from 0, this is used in p4est. In Fempar, local element numbers from 1.
    }
  }
  tree_offsets[local_tree_num] = k + 1;  // enumerate from 1 for Fortran
  
  return 0;
}

//get morton indices
int p4estw_get_quadrants(p4est_wrapper_t *p4estw, int *tree_indices, int *tree_offsets, int64_t *quadrant_morton_ids)
{
  int tt, q, Q,k, local_tree_num;
  p4est_tree_t       *tree;     /* Pointer to one octree */
  p4est_quadrant_t *quad;
  sc_array_t         *tquadrants;       /* Quadrant array for one tree */
  element_data_t     *elem_data;


  //printf("tree range %d ---- %d, num elements %d\n", p4estw->p4est->first_local_tree, p4estw->p4est->last_local_tree, p4estw->p4est->local_num_quadrants);
  for (tt = p4estw->p4est->first_local_tree, k = 0, local_tree_num = 0;
       tt <= p4estw->p4est->last_local_tree; ++tt, ++local_tree_num) {
//    printf(" treee %d\n", tt);
    tree_indices[local_tree_num] = tt + 1;  // enumerate from 1 for Fortran
    tree_offsets[local_tree_num] = k + 1;   // enumerate from 1 for Fortran
    tree = p4est_tree_array_index (p4estw->p4est->trees, tt);   /* Current tree */
    tquadrants = &tree->quadrants;
    Q = (p4est_locidx_t) tquadrants->elem_count;
    for (q = 0; q < Q; ++q, ++k) {
      quad = p4est_quadrant_array_index (tquadrants, q);
//      printf("quadrant %d, x %d, y %d, l %d\n", k, quad->x, quad->y, quad->level);
//       quadrant_data[DATA_PER_QUAD * k + 0] = quad->x;
//       quadrant_data[DATA_PER_QUAD * k + 1] = quad->y;
//       quadrant_data[DATA_PER_QUAD * k + 2] = quad->level;

      printf("level %d, x %d, y %d\n", quad->level, quad->x, quad->y);
      quadrant_morton_ids[k] = p4est_quadrant_linear_id(quad, P4EST_QMAXLEVEL);
      
      // probably we do not need to keep any special element data
//       elem_data = (element_data_t* ) quad->p.user_data;
//       elem_data->local_idx = k;   // enumerate from 0, this is used in p4est. In Fempar, local element numbers from 1.
    }
  }
  tree_offsets[local_tree_num] = k + 1;  // enumerate from 1 for Fortran
  
 return 0;
}


// // this is probably not needed, enumeration is done in p4estw_get_quadrants as well
// int p4estw_enumerate_quadrants(p4est_wrapper_t *p4estw)
// {
//   int tt, q, Q,k, local_tree_num;
//   p4est_tree_t       *tree;     /* Pointer to one octree */
//   p4est_quadrant_t *quad;
//   sc_array_t         *tquadrants;       /* Quadrant array for one tree */
//   element_data_t     *elem_data;
// 
// 
//   //printf("tree range %d ---- %d, num elements %d\n", p4estw->p4est->first_local_tree, p4estw->p4est->last_local_tree, p4estw->p4est->local_num_quadrants);
//   for (tt = p4estw->p4est->first_local_tree, k = 0, local_tree_num = 0;
//        tt <= p4estw->p4est->last_local_tree; ++tt, ++local_tree_num)
//   {
//     tree = p4est_tree_array_index (p4estw->p4est->trees, tt);   /* Current tree */
//     tquadrants = &tree->quadrants;
//     Q = (p4est_locidx_t) tquadrants->elem_count;
//     for (q = 0; q < Q; ++q, ++k) {
//       quad = p4est_quadrant_array_index (tquadrants, q);
//       elem_data = (element_data_t* ) quad->p.user_data;
//       elem_data->local_idx = k;   // enumerate from 0, this is used in p4est. In Fempar, local element numbers from 1.
// //       printf("proc: %d, local idx %d, global_start %ld\n", p4estw->p4est->mpirank, elem_data->local_idx, p4estw->p4est->global_first_quadrant[p4estw->p4est->mpirank]);
//     }
//   }
//   
//   return 0;
// }


int p4estw_get_ghost_sizes(p4est_wrapper_t *p4estw, int *num_ghosts, int *num_trees, int *num_processors)
{
  int idx;

  create_ghost_data(p4estw);
  p4est_ghost_t    *ghost = p4estw->ghost;

  *num_ghosts = ghost->ghosts.elem_count;

  *num_trees = 0;
  for (idx = 0; idx < ghost->num_trees; idx++)
  {
    if(ghost->tree_offsets[idx+1] > ghost->tree_offsets[idx])
      (*num_trees)++;
  }

  *num_processors = 0;
  for(idx = 0; idx < ghost->mpisize; idx++)
  {
    if(ghost->proc_offsets[idx+1] > ghost->proc_offsets[idx])
      (*num_processors)++;
  }

//  printf(" returning num ghosts %d, num trees %d, num proc %d\n", *num_ghosts, *num_trees, *num_processors);


  //  sprintf(str, "trees: ");
//  for (i = 0; i < p4estw->ghost->num_trees+1; i++)
//    sprintf(str, "%s%d, ", str, p4estw->ghost->tree_offsets[i]);

//  sprintf(str, "%s;  procs: ", str);
//  for(i = 0; i < p4estw->ghost->mpisize + 1; i++)
//    sprintf(str, "%s%d, ", str, p4estw->ghost->proc_offsets[i]);

//  printf("   +++proc %d:    %s\n\n\n",  p4estw->mpi_rank, str);

  return 0;

  
}

int p4estw_get_ghosts(p4est_wrapper_t *p4estw, int *tree_indices, int *tree_offsets, int *proc_indices, int *proc_offsets, int *quadrant_data, int64_t *global_indices)
{
  int idx, p4est_idx, process_id, l_proc_idx;
  p4est_ghost_t    *ghost = p4estw->ghost;
  p4est_quadrant_t *quad;

  tree_offsets[0] = 1;
  p4est_idx = 0;
  for(idx = 0; idx < ghost->num_trees; idx++)
  {
    //printf(" mpi %d: get ghosts idx %d, tree offset %d, next %d\n", p4estw->mpi_rank, idx, ghost->tree_offsets[idx], ghost->tree_offsets[idx+1]);

    if(ghost->tree_offsets[idx+1] > ghost->tree_offsets[idx])
    {
      tree_indices[p4est_idx] = idx + 1;
      tree_offsets[p4est_idx + 1] = ghost->tree_offsets[idx+1] + 1;
      //printf(" mpi %d: get ghosts p4estidx %d, ghost->tree_offsets[idx+1] %d, tree_offsets[p4est_idx] %d\n", p4estw->mpi_rank,  p4est_idx, ghost->tree_offsets[idx+1], tree_offsets[p4est_idx]);
      p4est_idx++;
    }
  }

  proc_offsets[0] = 1;
  p4est_idx = 0;
  for(idx = 0; idx < ghost->mpisize; idx++)
  {
    if(ghost->proc_offsets[idx+1] > ghost->proc_offsets[idx])
    {
      proc_indices[p4est_idx] = idx;
      proc_offsets[p4est_idx + 1] = ghost->proc_offsets[idx+1] + 1;
      p4est_idx++;
    }
  }

  l_proc_idx = 0;
  for(idx = 0; idx < ghost->ghosts.elem_count; idx++)
  {
    if(proc_offsets[l_proc_idx + 1] == idx + 1)
      l_proc_idx += 1;
    process_id = proc_indices[l_proc_idx];
    quad = (p4est_quadrant_t*)sc_array_index(&p4estw->ghost->ghosts, idx);
    quadrant_data[DATA_PER_QUAD*idx + 0] = quad->x;
    quadrant_data[DATA_PER_QUAD*idx + 1] = quad->y;
    quadrant_data[DATA_PER_QUAD*idx + 2] = quad->level;
    global_indices[idx] = quad->p.piggy3.local_num + p4estw->p4est->global_first_quadrant[process_id] + 1; // +1 for Fortran
  }

  return 0;
}

//*******************************************************************************************
// Creating of the p4est wrapper
//*******************************************************************************************

p4est_connectivity_t *create_connectivity(int num_elements, int num_points, int* list_nodes, double coordinates[][NUM_DIMS]);
p4est_wrapper_t* p4estw_create_from_connectivity(p4est_connectivity_t *connectivity, int fortran_mpicomm);


p4est_wrapper_t* p4estw_create_from_geometry(int num_elements, int num_points, int* list_nodes, double coordinates[][NUM_DIMS], int fortran_mpicomm)
{
  //printf(" creating from geometry, %d elements and %d nodes\n", num_elements, num_points);
  p4est_connectivity_t *connectivity = create_connectivity(num_elements, num_points, list_nodes, coordinates);
  return p4estw_create_from_connectivity(connectivity, fortran_mpicomm);
}

p4est_wrapper_t* p4estw_create_unit_cube(int fortran_mpicomm)
{
  //printf(" creating from geometry, %d elements and %d nodes\n", num_elements, num_points);
  p4est_connectivity_t *connectivity;
  
  if(NUM_DIMS == 2)
    connectivity = p4est_connectivity_new_unitsquare ();
  else if(NUM_DIMS == 3) 
    assert(0);
    //connectivity = p8est_connectivity_new_unitcube ();
  else 
    assert(0);
  
  return p4estw_create_from_connectivity(connectivity, fortran_mpicomm);
}


p4est_wrapper_t* p4estw_create_from_connectivity(p4est_connectivity_t *connectivity, int fortran_mpicomm)
{
  p4est_wrapper_t *p4estw = malloc(sizeof(p4est_wrapper_t));
  sc_MPI_Comm mpicomm = MPI_Comm_f2c(fortran_mpicomm);

  p4estw->p4est = p4est_new (mpicomm, connectivity, sizeof(element_data_t), NULL/*init_element_data*/, NULL);
  p4estw->ghost = NULL;
  p4estw->ghost_data = NULL;
  p4estw->lnodes = NULL;
  p4estw->info_found = 0;
  sc_MPI_Comm_rank(p4estw->p4est->mpicomm, &p4estw->mpi_rank);
  sc_MPI_Comm_size(p4estw->p4est->mpicomm, &p4estw->mpi_num_proc);

  return p4estw;
}




//*******************************************************************************************
// Creating of the connectivity
// This currently works only for 2D, has to be extended/rewised
//*******************************************************************************************

typedef struct vertex_data
{
  int num_elements;
  int corner;
  int *element_idx;
}
vertex_data_t;

const int neighbour_vertices[NUM_ELEM_NODES][NUM_NEIGHBOURING_VERTICES] = {{1,2},{0,3},{0,3},{1,2}};


int local_vertex_id(int elem, int vertex, p4est_topidx_t *tree_to_vertex)
{
  int lpoint;

  for(lpoint = 0; lpoint < NUM_ELEM_NODES; lpoint++)
  {
    if(tree_to_vertex[elem * NUM_ELEM_NODES + lpoint] == vertex)
      return lpoint;
  }

  return DOES_NOT_EXIST;
}

const int local_face_id_map[NUM_ELEM_NODES][NUM_ELEM_NODES] = {{-1,2,0,-1},{6,-1,-1,1},{4,-1,-1,3},{-1,5,7,-1}};


int local_face_id(int local_node_1, int local_node_2)
{
  int face_idx = local_face_id_map[local_node_1][local_node_2];
  return face_idx;
}



p4est_connectivity_t *create_connectivity(int num_elements, int num_points, int* list_nodes, double coordinates[][NUM_DIMS])
{
  p4est_connectivity_t *connectivity;
  //const int fempar_to_p4est_elem_vertex_map[NUM_ELEM_NODES] = {0,1,3,2};
  // now fempar uses the same numbering as p4est, so this is no longer needed
  const int fempar_to_p4est_elem_vertex_map[NUM_ELEM_NODES] = {0,1,2,3};
  vertex_data_t vertices[num_points];

  p4est_topidx_t tree_to_vertex[num_elements * NUM_ELEM_NODES];
  p4est_topidx_t tree_to_tree[num_elements * NUM_ELEM_FACES];
  int8_t         tree_to_face[num_elements * NUM_ELEM_FACES];
  p4est_topidx_t *tree_to_corner = NULL;

  p4est_topidx_t num_corners;
  p4est_topidx_t *ctt_offset = NULL;
  p4est_topidx_t *corner_to_tree = NULL;
  int8_t  *corner_to_corner = NULL;

  int lelem, elem, lelem2, elem2, elem2_lpoint1, elem2_lpoint2, lpoint, lpoint2, point, point2, lface, lface2, idx, fempar_idx, neighbour_vertices_idx, ori, corner;

  double *coordinates_p4est;
  // initialize by zero
  for(point = 0; point < num_points; point++)
  {
    vertices[point].num_elements = 0;
    vertices[point].element_idx = NULL;
    vertices[point].corner = UNNUMBERED_CORNER;
  }


  // convert to another local node order and enumerate from 0.
  for(elem = 0; elem < num_elements; elem++)
  {
    for(lpoint = 0; lpoint < NUM_ELEM_NODES; lpoint++)
    {
      idx = NUM_ELEM_NODES * elem + lpoint;
      fempar_idx = NUM_ELEM_NODES * elem + fempar_to_p4est_elem_vertex_map[lpoint];
      tree_to_vertex[idx] = list_nodes[fempar_idx] - 1;
//      printf("idx %d, fempar idx %d, ln %d, ttv %d, ttv[0] %d\n", idx, fempar_idx, list_nodes[fempar_idx], tree_to_vertex[idx], tree_to_vertex[0]);
    }

    for(lface = 0; lface < NUM_ELEM_FACES; lface++)
    {
      idx = NUM_ELEM_FACES * elem + lface;
      tree_to_tree[idx] = elem;
    }
  }



  // first find out, how many adjacent elements each vertex has
  for(idx = 0; idx < num_elements * NUM_ELEM_NODES; idx++)
  {
    //printf("idx %d, ttv %d\n", idx, tree_to_vertex[idx]);
    vertices[tree_to_vertex[idx]].num_elements += 1;
  }
//  printf("\n\n");

  // alocate space for adjacent elements, faces and vertices on oposite ends of those faces
  for(point = 0; point < num_points; point++)
  {
    //printf("point %d, num elems %d\n", point, vertices[point].num_elements);
    assert(vertices[point].num_elements > 0);
    vertices[point].element_idx = malloc(vertices[point].num_elements * sizeof(int));

    // maximal number of faces is elems + 1. if num faces == num_elems, it is interior vertex (and thus corner in p4est sense)
    //vertices[point].faces_idx = malloc((vertices[point].num_elements + 1) * sizeof(int));
    //vertices[point].face_vertices_idx = malloc((vertices[point].num_elements + 1) * sizeof(int));

    // to be used as counter in the next step
    vertices[point].num_elements = 0;
  }

  // fill array of adjacent elements
  for(elem = 0; elem < num_elements; elem++)
  {
    for(lpoint = 0; lpoint < NUM_ELEM_NODES; lpoint++)
    {
      idx = NUM_ELEM_NODES * elem + lpoint;
      vertices[tree_to_vertex[idx]].element_idx[vertices[tree_to_vertex[idx]].num_elements] = elem;
      vertices[tree_to_vertex[idx]].num_elements++;
    }
  }


  // find vertices connected by face with each vertex
  for(point = 0; point < num_points; point++)
  {
    //printf("* point %d\n", point);

    for(lelem = 0; lelem < vertices[point].num_elements; lelem++)
    {
      // element adjacent to given vertex
      elem = vertices[point].element_idx[lelem];
      //printf("** elem %d\n", elem);

      // find local idx of the given vertex within this element
      lpoint = local_vertex_id(elem, point, tree_to_vertex);
      assert(lpoint != DOES_NOT_EXIST);

      // those are its neigbouring nodes within this element. Check if they are allready stored
      for(neighbour_vertices_idx = 0; neighbour_vertices_idx < 2; neighbour_vertices_idx++)
      {
        point2 = tree_to_vertex[elem * NUM_ELEM_NODES + neighbour_vertices[lpoint][neighbour_vertices_idx]];
        //printf("*** point2 %d\n", point2);
        lpoint2 = local_vertex_id(elem, point2, tree_to_vertex);
        assert(lpoint2 != DOES_NOT_EXIST);

        // find the second element
        for(lelem2 = 0; lelem2 < vertices[point].num_elements; lelem2++)
        {
          elem2 = vertices[point].element_idx[lelem2];
          if(elem == elem2)
            continue;

          // check if the second element contains first vertex
          elem2_lpoint1 = local_vertex_id(elem2, point, tree_to_vertex);
          //printf("**** elem2 %d, elem2lpoin1 %d\n", elem2, elem2_lpoint1);
          // if does not contain, try next element
          if(elem2_lpoint1 == DOES_NOT_EXIST)
            continue;

          // check if the second element contains second vertex
          elem2_lpoint2 = local_vertex_id(elem2, point2, tree_to_vertex);
          //printf("               elem2lpoin2 %d\n", elem2_lpoint2);
          // if does not contain, try next element
          if(elem2_lpoint2 == DOES_NOT_EXIST)
            continue;

          // I have found it
          break;
        }

        lface = local_face_id(lpoint, lpoint2);
        if(lface == -1)
        {
          printf("wrong lface from nodes (%d, %d), local  nodes (%d, %d), elements (>%d<)\n", point, point2, lpoint, lpoint2, elem);
          assert(0);
        }

        // check whether I have found neighbouring element
        if(lelem2 == vertices[point].num_elements)
        {
          //printf("***** boundary edge, elem %d, local face id %d\n", elem, lface);
          // I did not, we are on the boundary
          elem2 = -1;
          assert(tree_to_tree[elem * NUM_ELEM_FACES + lface%4] == elem);
          tree_to_face[elem * NUM_ELEM_FACES + lface%4] = lface%4;

          // a vertex connected to boundary edge is not a corner in the sense of p4est
          vertices[point].corner = IS_NOT_CORNER;
        }
        else
        {
          //printf("***** inner edge, elem %d\n", elem);

          // found neighbouing element, face between two elements
          lface2 = local_face_id(elem2_lpoint1, elem2_lpoint2);
          if(lface2 == -1)
          {
            printf("wrong lface from nodes (%d, %d), local  nodes (%d, %d), elements (%d, >%d<)\n", point, point2, elem2_lpoint1, elem2_lpoint2, elem, elem2);
            assert(0);
          }

          if(lface/4 == lface2/4)
            ori = 0;
          else
            ori = 1;
          // set the real orientation flag as is (hopefully) seen by p4est
          lface = 4*ori + lface%4;
          lface2 = 4*ori + lface2%4;

          tree_to_tree[elem * NUM_ELEM_FACES + lface%4] = elem2;
          tree_to_face[elem * NUM_ELEM_FACES + lface%4] = lface2;
        }
      }
    }
  }


  num_corners = 0;
  for(point = 0; point < num_points; point++)
    if(vertices[point].corner == UNNUMBERED_CORNER)
      num_corners++;

  ctt_offset = malloc((num_corners + 1) * sizeof(p4est_topidx_t));
  ctt_offset[0] = 0;

  if(num_corners > 0)
  {
    corner_to_tree = malloc(num_elements * num_corners * sizeof(p4est_topidx_t));
    corner_to_corner = malloc(num_elements * num_corners * sizeof(int8_t));
    tree_to_corner = malloc(num_elements * NUM_ELEM_NODES * sizeof(p4est_topidx_t));

    corner = 0;
    for(point = 0; point < num_points; point++)
    {
      if(vertices[point].corner == UNNUMBERED_CORNER)
      {
        vertices[point].corner = corner;
        assert(vertices[point].num_elements > 2);
        ctt_offset[corner + 1] = ctt_offset[corner] + vertices[point].num_elements;
        for(lelem = 0; lelem < vertices[point].num_elements; lelem++)
        {
          elem = vertices[point].element_idx[lelem];
          idx = ctt_offset[corner] + lelem;
          corner_to_tree[idx] = elem;
          corner_to_corner[idx] = local_vertex_id(elem, point, tree_to_vertex);
        }
        corner++;
      }
    }

    for(elem = 0; elem < num_elements; elem++)
    {
      for(lpoint = 0; lpoint < NUM_ELEM_NODES; lpoint++)
      {
        idx = NUM_ELEM_NODES*elem + lpoint;
        point = tree_to_vertex[idx];
        assert(vertices[point].corner != UNNUMBERED_CORNER);
        if(vertices[point].corner >= 0)
          tree_to_corner[idx] = vertices[point].corner;
        else
          tree_to_corner[idx] = -1;
      }
    }
  }

//  idx = corner_to_tree[2]; corner_to_tree[2] = corner_to_tree[3]; corner_to_tree[3] = idx;
//  idx = corner_to_corner[2]; corner_to_corner[2] = corner_to_corner[3]; corner_to_corner[3] = idx;

//  printf("Num points %d, num elements %d, num corners %d\n", num_points, num_elements, num_corners);

//  printf("tree_to_vertex: \n");
//  for(elem = 0; elem < num_elements; elem++)
//  {
//    for(lpoint = 0; lpoint < NUM_ELEM_NODES; lpoint++)
//    {
//      idx = NUM_ELEM_FACES * elem + lpoint;
//      printf("%d, ", tree_to_vertex[idx]);
//    }
//    printf("\n");
//  }
//  printf("\n");

//  printf("tree_to_tree: \n");
//  for(elem = 0; elem < num_elements; elem++)
//  {
//    for(lpoint = 0; lpoint < NUM_ELEM_NODES; lpoint++)
//    {
//      idx = NUM_ELEM_FACES * elem + lpoint;
//      printf("%d, ", tree_to_tree[idx]);
//    }
//    printf("\n");
//  }
//  printf("\n");

//  printf("tree_to_face: \n");
//  for(elem = 0; elem < num_elements; elem++)
//  {
//    for(lpoint = 0; lpoint < NUM_ELEM_NODES; lpoint++)
//    {
//      idx = NUM_ELEM_FACES * elem + lpoint;
//      printf("%d, ", tree_to_face[idx]);
//    }
//    printf("\n");
//  }
//  printf("\n");

//  if(num_corners > 0)
//  {
//    printf("tree_to_corner: \n");
//    for(elem = 0; elem < num_elements; elem++)
//    {
//      for(lpoint = 0; lpoint < NUM_ELEM_NODES; lpoint++)
//      {
//        idx = NUM_ELEM_FACES * elem + lpoint;
//        printf("%d, ", tree_to_corner[idx]);
//      }
//      printf("\n");
//    }
//    printf("\n");
//  }

//  if(num_corners > 0)
//  {
//    printf("ctt_offset: \n");
//    for(corner = 0; corner < num_corners + 1; corner++)
//    {
//      printf("%d, ",ctt_offset[corner]);
//    }
//    printf("\n\n");

//    printf("corner_to_tree: \n");
//    for(corner = 0; corner < num_corners; corner++)
//    {
//      for(idx = ctt_offset[corner]; idx < ctt_offset[corner+1]; idx++)
//      {
//        printf("%d, ", corner_to_tree[idx]);
//      }
//      printf("\n");
//    }
//    printf("\n");

//    printf("corner_to_corner: \n");
//    for(corner = 0; corner < num_corners; corner++)
//    {
//      for(idx = ctt_offset[corner]; idx < ctt_offset[corner+1]; idx++)
//      {
//        printf("%d, ", corner_to_corner[idx]);
//      }
//      printf("\n");
//    }
//    printf("\n");
//  }

  for(point = 0; point < num_points; point++)
  {
    if(vertices[point].element_idx != NULL)
      free(vertices[point].element_idx);
  }

  coordinates_p4est = malloc(num_points * 3 * sizeof(double));
  memset(coordinates_p4est, 0, num_points * 3 * sizeof(double));

  connectivity =  p4est_connectivity_new_copy (num_points, num_elements, num_corners,
                                      coordinates_p4est, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      tree_to_corner, ctt_offset,
                                      corner_to_tree, corner_to_corner);

  free(coordinates_p4est);


  if(ctt_offset != NULL)
    free(ctt_offset);
  if(corner_to_tree != NULL)
    free(corner_to_tree);
  if(corner_to_corner != NULL)
    free(corner_to_corner);
  if(tree_to_corner != NULL)
    free(tree_to_corner);

  return connectivity;
}



//*******************************************************************************************
// Additional data
//*******************************************************************************************

void create_ghost_data(p4est_wrapper_t* p4estw)
{
  if((p4estw->ghost == NULL) || (p4estw->ghost_data == NULL))
  {
    assert(p4estw->ghost == NULL);
    assert(p4estw->ghost_data == NULL);

    /* Create the ghost layer to learn about parallel neighbors. */
    p4estw->ghost = p4est_ghost_new (p4estw->p4est, P4EST_CONNECT_FULL);
    /* create space for storing the ghost data */
    p4estw->ghost_data = P4EST_ALLOC (element_data_t, p4estw->ghost->ghosts.elem_count);
    /* synchronize the ghost data */
    p4est_ghost_exchange_data (p4estw->p4est, p4estw->ghost, p4estw->ghost_data);
  }
}

void destroy_additional_data(p4est_wrapper_t* p4estw)
{
  if(p4estw->ghost != NULL)
  {
    p4est_ghost_destroy (p4estw->ghost);
    p4estw->ghost = NULL;
  }

  if(p4estw->ghost_data != NULL)
  {
    P4EST_FREE (p4estw->ghost_data);
    p4estw->ghost_data = NULL;
  }

  if(p4estw->lnodes != NULL)
  {
    p4est_lnodes_destroy (p4estw->lnodes);
    p4estw->lnodes = NULL;
  }
}


//*******************************************************************************************
// Refinements
//*******************************************************************************************

static int
refine_all_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  return 1;
}

static void
init_element_data (p4est_t * p4est, p4est_topidx_t which_tree,
                              p4est_quadrant_t * q)
{
  /* the data associated with a quadrant is accessible by p.user_data */
  element_data_t       *data = (element_data_t *) q->p.user_data;

  data->global_idx = -12345;
  data->local_idx = -12345;
  data->refinement_flag = REF_FLAG_INVALID;
}

int p4estw_refine_all(p4est_wrapper_t* p4estw)
{
  destroy_additional_data(p4estw);

  int recursive = 0;

  printf("*********** refine all \n");
  p4est_refine (p4estw->p4est, recursive, refine_all_fn, init_element_data);
  p4est_balance (p4estw->p4est, P4EST_CONNECT_FACE, init_element_data);
  return 0;
}

int p4estw_partition(p4est_wrapper_t* p4estw)
{
  destroy_additional_data(p4estw);

  int partforcoarsen = 0;
  p4est_partition (p4estw->p4est, partforcoarsen, NULL);
  return 0;
}

#endif
