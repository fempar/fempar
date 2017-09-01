/* Module file taken "AS-IS" and adapted to FEMPAR needs from p4est_wrapper.c 
   https://github.com/cburstedde/hopest 4feed803f0c61564203a7bc3f2ca1a6adb63d3cd */

/*
  This file is part of hopest.
  hopest is a Fortran/C library and application for high-order mesh
  preprocessing and interfacing to the p4est apaptive mesh library.

  Copyright (C) 2014 by the developers.

  hopest is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  hopest is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with hopest; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/ 

#ifdef ENABLE_P4EST

#include <sc.h>
#include <p4est.h>
#include <p4est_extended.h>
#include <p4est_mesh.h>
#include <p4est_bits.h>
#include <p4est_base.h>
#include "p4est_wrapper.h"
#include <p8est.h>
#include <p8est_extended.h>
#include <p8est_mesh.h>
#include <p8est_bits.h>
#ifndef SC_ENABLE_MPI
  static int sc_mpi_initialized = 0;
#else
#endif

// Init p4est environment 
// TODO: pass a Fortran communicator and transform it into a C communicator
//       http://stackoverflow.com/questions/42530620/how-to-pass-mpi-communicator-handle-from-fortran-to-c-using-iso-c-binding
void F90_p4est_init()
{
#ifndef SC_ENABLE_MPI
  /* Initialize MPI; see sc_mpi.h.
   * If configure --enable-mpi is given these are true MPI calls.
   * Else these are dummy functions that simulate a single-processor run. */
  int mpiret;
  
  // Assume that we have a dummy main program called p4est_init_environment
  int argc      = 1;
  char ** argv  = (char *[]) {"p4est_init_environment"};
  
  if (!sc_mpi_initialized)
  {
    mpiret = sc_MPI_Init (&argc, &argv);
    SC_CHECK_MPI (mpiret==1);
    
    /* These 2x functions are optional.  If called they store the MPI rank as a
     * static variable so subsequent global p4est log messages are only issued
     * from processor zero.  Here we turn off most of the logging; see sc.h. */
    sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
    p4est_init (NULL, SC_LP_PRODUCTION);
    
    sc_mpi_initialized = 1;
  }
#else
#endif  
}

// Finalize p4est_environment
void F90_p4est_finalize()
{
#ifndef SC_ENABLE_MPI
    int mpiret;
    if ( sc_mpi_initialized )
    {    
        sc_finalize ();   
        mpiret = sc_MPI_Finalize ();
        SC_CHECK_MPI (mpiret);    
    } 
#else  
#endif    
}

void F90_p4est_connectivity_new_unitsquare(p4est_connectivity_t **p4est_connectivity)
{
  F90_p4est_connectivity_destroy(p4est_connectivity);
  *p4est_connectivity = p4est_connectivity_new_unitsquare();
  P4EST_ASSERT (p4est_connectivity_is_valid (*p4est_connectivity));
}

void F90_p8est_connectivity_new_unitcube(p8est_connectivity_t **p8est_connectivity)
{
  F90_p8est_connectivity_destroy(p8est_connectivity);
  *p8est_connectivity = p8est_connectivity_new_unitcube();
  P4EST_ASSERT (p8est_connectivity_is_valid (*p8est_connectivity));
}

void F90_p4est_new ( p4est_connectivity_t *conn,
                     p4est_t              **p4est_out)
{
  p4est_t* p4est;
  P4EST_ASSERT (p4est_connectivity_is_valid (conn));
  
  // TODO: from where are we going to extract the MPI communicator in case
  //       SC_ENABLE_MPI is defined???
  /* Create a forest that is not refined; it consists of the root octant. */                                        
  p4est = p4est_new (sc_MPI_COMM_WORLD, conn, 0, NULL, NULL);
  *p4est_out = p4est;
}

void F90_p8est_new ( p8est_connectivity_t *conn,
                     p8est_t              **p8est_out)
{
  p8est_t* p8est;
  P4EST_ASSERT (p8est_connectivity_is_valid (conn));
  
  // TODO: from where are we going to extract the MPI communicator in case
  //       SC_ENABLE_MPI is defined???
  /* Create a forest that is not refined; it consists of the root octant. */                                        
  p8est = p8est_new (sc_MPI_COMM_WORLD, conn, 0, NULL, NULL);
  *p8est_out = p8est;
}


void init_fn_callback_2d(p4est_t * p4est,p4est_topidx_t which_tree,p4est_quadrant_t * quadrant)
{
    p4est_tree_t       *tree;
    p4est_quadrant_t   *q;
    sc_array_t         *quadrants;
    int                *user_pointer;
    int                *quadrant_data;
    int                 output;
    
    P4EST_ASSERT(which_tree == 0);
    
    // Extract a reference to the first (and uniquely allowed) tree
    tree = p4est_tree_array_index (p4est->trees,0);
    quadrants = &(tree->quadrants);
    q = p4est_quadrant_array_index(quadrants, current_quadrant_index);
    P4EST_ASSERT(p4est_quadrant_compare(q,quadrant) == 0);
    user_pointer  = (int *) p4est->user_pointer;
    quadrant_data = (int *) quadrant->p.user_data;
    *quadrant_data = user_pointer[current_quadrant_index];
    current_quadrant_index = (current_quadrant_index+1) % (quadrants->elem_count);    
}

void init_fn_callback_3d(p8est_t * p8est,p4est_topidx_t which_tree,p8est_quadrant_t * quadrant)
{
    p8est_tree_t       *tree;
    p8est_quadrant_t   *q;
    sc_array_t         *quadrants;
    int                *user_pointer;
    int                *quadrant_data;
    int                 output;
    
    P4EST_ASSERT(which_tree == 0);
    
    // Extract a reference to the first (and uniquely allowed) tree
    tree = p8est_tree_array_index (p8est->trees,0);
    quadrants = &(tree->quadrants);
    q = p8est_quadrant_array_index(quadrants, current_quadrant_index);
    P4EST_ASSERT(p8est_quadrant_compare(q,quadrant) == 0);
    user_pointer  = (int *) p8est->user_pointer;
    quadrant_data = (int *) quadrant->p.user_data;
    *quadrant_data = user_pointer[current_quadrant_index];
    current_quadrant_index = (current_quadrant_index+1) % (quadrants->elem_count);    
}


void F90_p4est_set_user_pointer(int * user_data, p4est_t * p4est) 
{
    p4est_reset_data(p4est,sizeof(int),init_fn_callback_2d,(void *)user_data);
}

void F90_p8est_set_user_pointer(int * user_data, p8est_t * p8est) 
{
    p8est_reset_data(p8est,sizeof(int),init_fn_callback_3d,(void *)user_data);
}

void F90_p4est_mesh_new(p4est_t  *p4est,
                        p4est_mesh_t  **mesh_out )
{
  p4est_mesh_t       *mesh;
  p4est_ghost_t      *ghost;
  
  F90_p4est_mesh_destroy(mesh_out);

  //create ghost layer and mesh
  ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
  mesh = p4est_mesh_new (p4est,ghost,P4EST_CONNECT_FULL);

  p4est_ghost_destroy(ghost);
  
  //return mesh as pointer address;
  *mesh_out=(p4est_mesh_t *)mesh;
}

void F90_p8est_mesh_new(p8est_t  *p8est,
                        p8est_mesh_t  **mesh_out )
{
  p8est_mesh_t       *mesh;
  p8est_ghost_t      *ghost;
  
  F90_p8est_mesh_destroy(mesh_out);

  //create ghost layer and mesh
  ghost = p8est_ghost_new (p8est, P8EST_CONNECT_FULL);
  mesh = p8est_mesh_new (p8est,ghost,P8EST_CONNECT_FULL);

  p8est_ghost_destroy(ghost);
  
  //return mesh as pointer address;
  *mesh_out=(p8est_mesh_t *)mesh;
}

void F90_p4est_connectivity_destroy(p4est_connectivity_t **p4est_connectivity)
{
    if (*p4est_connectivity) p4est_connectivity_destroy(*p4est_connectivity);
}

void F90_p8est_connectivity_destroy(p8est_connectivity_t **p8est_connectivity)
{
    if (*p8est_connectivity) p8est_connectivity_destroy(*p8est_connectivity);
}

void F90_p4est_destroy(p4est_t **p4est)
{
    if (*p4est) p4est_destroy(*p4est);
}

void F90_p8est_destroy(p8est_t **p8est)
{
    if (*p8est) p8est_destroy(*p8est);
}

void F90_p4est_mesh_destroy(p4est_mesh_t **p4est_mesh)
{
    if (*p4est_mesh) 
    {    
      p4est_mesh_destroy(*p4est_mesh);
    }   
}

void F90_p8est_mesh_destroy(p8est_mesh_t **p8est_mesh)
{
    if (*p8est_mesh) 
    {    
      p8est_mesh_destroy(*p8est_mesh);
    }   
}

void F90_p4est_get_mesh_info (p4est_t        *p4est,
                              p4est_mesh_t   *mesh,
                              p4est_locidx_t *local_num_quadrants,
                              p4est_gloidx_t *global_num_quadrants,
                              p4est_gloidx_t *global_first_quadrant,
                              p4est_locidx_t *num_half_faces)
{
    SC_CHECK_ABORTF (mesh->local_num_quadrants == p4est->local_num_quadrants,
                     "mesh->local_num_quadrants [%d] and p4est->local_num_quadrants mismatch [%d]!",
                     mesh->local_num_quadrants,  p4est->local_num_quadrants);
    
    SC_CHECK_ABORTF (p4est->trees->elem_count == 1,
                     "p4est with more [%ld] than one tree!", p4est->trees->elem_count);
    
    *local_num_quadrants   = p4est->local_num_quadrants;
    *global_num_quadrants  = p4est->global_num_quadrants;
    *global_first_quadrant = p4est->global_first_quadrant[p4est->mpirank];
    *num_half_faces        = mesh->quad_to_half->elem_count;
}

void F90_p8est_get_mesh_info (p8est_t        *p8est,
                              p8est_mesh_t   *mesh,
                              p4est_locidx_t *local_num_quadrants,
                              p4est_gloidx_t *global_num_quadrants,
                              p4est_gloidx_t *global_first_quadrant,
                              p4est_locidx_t *num_half_faces)
{
    SC_CHECK_ABORTF (mesh->local_num_quadrants == p8est->local_num_quadrants,
                     "mesh->local_num_quadrants [%d] and p8est->local_num_quadrants mismatch [%d]!",
                     mesh->local_num_quadrants,  p8est->local_num_quadrants);
    
    SC_CHECK_ABORTF (p8est->trees->elem_count == 1,
                     "p8est with more [%ld] than one tree!", p8est->trees->elem_count);
    
    *local_num_quadrants   = p8est->local_num_quadrants;
    *global_num_quadrants  = p8est->global_num_quadrants;
    *global_first_quadrant = p8est->global_first_quadrant[p8est->mpirank];
    *num_half_faces        = mesh->quad_to_half->elem_count;
}


void F90_p4est_get_mesh_topology_arrays( p4est_t        *p4est,
                                         p4est_mesh_t   *mesh,
                                         p4est_locidx_t **quad_to_quad,
                                         int8_t         **quad_to_face, 
                                         p4est_locidx_t **quad_to_half, 
                                         p4est_locidx_t **quad_to_corner,
                                         p4est_qcoord_t *quadcoords,
                                         int8_t         *quadlevel ) 
{
  int iquad,iquadloc;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *q;
  sc_array_t         *quadrants;
  
  // Extract a reference to the first (and uniquely allowed) tree
  tree = p4est_tree_array_index (p4est->trees,0);
  for (iquad = 0; iquad < mesh->local_num_quadrants; iquad++) {  
    quadrants = &(tree->quadrants);
    iquadloc = iquad - tree->quadrants_offset;
    q = p4est_quadrant_array_index(quadrants, iquadloc);
    quadlevel [iquad    ] = q->level;
    quadcoords[iquad*2  ] = q->x;
    quadcoords[iquad*2+1] = q->y;
  }

  *quad_to_quad=mesh->quad_to_quad;
  *quad_to_face=mesh->quad_to_face;
  *quad_to_half = NULL;
  if(mesh->quad_to_half->elem_count>0) *quad_to_half = (p4est_locidx_t *) mesh->quad_to_half->array;
  *quad_to_corner=mesh->quad_to_corner;
}

void F90_p8est_get_mesh_topology_arrays( p8est_t        *p8est,
                                         p8est_mesh_t   *mesh,
                                         p4est_locidx_t **quad_to_quad,
                                         int8_t         **quad_to_face, 
                                         p4est_locidx_t **quad_to_half, 
                                         p4est_locidx_t *quad_to_quad_by_edge,
                                         int8_t         *quad_to_edge,
                                         p4est_locidx_t *quad_to_half_by_edge,
                                         p4est_locidx_t **quad_to_corner,
                                         p4est_qcoord_t *quadcoords,
                                         int8_t         *quadlevel ) 
{
  int iquad,iquadloc;
  p8est_tree_t       *tree;
  p8est_quadrant_t   *q;
  sc_array_t         *quadrants;
  edge_info_t edge_info;
  
  // Extract a reference to the first (and uniquely allowed) tree
  tree = p8est_tree_array_index (p8est->trees,0);
  for (iquad = 0; iquad < mesh->local_num_quadrants; iquad++) {  
    quadrants = &(tree->quadrants);
    iquadloc = iquad - tree->quadrants_offset;
    q = p8est_quadrant_array_index(quadrants, iquadloc);
    quadlevel [iquad    ] = q->level;
    quadcoords[iquad*3  ] = q->x;
    quadcoords[iquad*3+1] = q->y;
    quadcoords[iquad*3+2] = q->z;
  }

  // Extract the neighbor info for edges. Initialize it to -1 (like in quad_to_corner)
  for(int i=0;i<12*(mesh->local_num_quadrants);i++)
  {
    quad_to_quad_by_edge[i] = -1;
    quad_to_edge[i] = -1;
  }
  for(int i=0;i<2*(mesh->local_num_quadrants);i++)
  {
    quad_to_half_by_edge[i] = -1;
  }


  edge_info.quad_to_quad_by_edge = quad_to_quad_by_edge;
  edge_info.quad_to_edge         = quad_to_edge;
  edge_info.quad_to_half_by_edge = quad_to_half_by_edge;
  p8est_iterate(p8est, NULL, &edge_info, NULL, NULL,edge_callback, NULL);

  *quad_to_quad=mesh->quad_to_quad;
  *quad_to_face=mesh->quad_to_face;
  *quad_to_half = NULL;
  if(mesh->quad_to_half->elem_count>0) *quad_to_half = (p4est_locidx_t *) mesh->quad_to_half->array;
  *quad_to_corner=mesh->quad_to_corner;
}

void edge_callback(p8est_iter_edge_info_t * info, void * user_data)
{
  p8est_iter_edge_side_t * cells_around;
  edge_info_t *edge_info;
  p4est_locidx_t *quad_to_quad_by_edge;
  int8_t         *quad_to_edge;
  p4est_locidx_t *quad_to_half_by_edge;

  int k;
  p4est_locidx_t ineig[4], jneig[4];
  int8_t ineig_iedge[4], jneig_jedge[4];
  
  for(int i=0;i<2;i++) ineig[i]   = -1;
  for(int i=0;i<2;i++) jneig[i]   = -1;

  P4EST_ASSERT( (info->sides.elem_count) <= 4 );

  edge_info = (edge_info_t *) user_data;
  quad_to_quad_by_edge = edge_info->quad_to_quad_by_edge;
  quad_to_edge = edge_info->quad_to_edge;
  quad_to_half_by_edge = edge_info->quad_to_half_by_edge;
  
  cells_around = (p8est_iter_edge_side_t *) info->sides.array;
  
  // First treat boundary edges
  if ( info->sides.elem_count == 1 || info->sides.elem_count == 2 )
  {
      for(int i=0;i<(info->sides.elem_count);i++)
      {
          if (cells_around[i].is_hanging)
          {
              ineig[0]       = cells_around[i].is.hanging.quadid[0];
              ineig[1]       = cells_around[i].is.hanging.quadid[1];
              ineig_iedge[0] = cells_around[i].edge;
              quad_to_quad_by_edge[ 12*ineig[0] + ineig_iedge[0] ] = ineig[0];
              quad_to_edge        [ 12*ineig[0] + ineig_iedge[0] ] = ineig_iedge[0];
              quad_to_quad_by_edge[ 12*ineig[1] + ineig_iedge[0] ] = ineig[1];
              quad_to_edge        [ 12*ineig[1] + ineig_iedge[0] ] = ineig_iedge[0];
          }
          else
          {
              ineig[0]       = cells_around[i].is.full.quadid;
              ineig_iedge[0] = cells_around[i].edge;
              quad_to_quad_by_edge[ 12*ineig[0] + ineig_iedge[0] ] = ineig[0];
              quad_to_edge        [ 12*ineig[0] + ineig_iedge[0] ] = ineig_iedge[0];  
          } 
      }
      return;
  }
  
  k=0;
  for(int i=0;i<(info->sides.elem_count);i++)
  {
    if (cells_around[i].is_hanging)
    {
      ineig[2*k]       = cells_around[i].is.hanging.quadid[0];
      ineig[2*k+1]     = cells_around[i].is.hanging.quadid[1];
      ineig_iedge[2*k] = cells_around[i].edge;
    }
    else
    {
      ineig[2*k]       = cells_around[i].is.full.quadid;
      ineig_iedge[2*k] = cells_around[i].edge;
    } 

    for(int j=i;j<(info->sides.elem_count);j++)
    {
        if ( (cells_around[i].faces[0] != cells_around[j].faces[0]) &&
             (cells_around[i].faces[0] != cells_around[j].faces[1]) &&
             (cells_around[i].faces[1] != cells_around[j].faces[0]) &&
             (cells_around[i].faces[1] != cells_around[j].faces[1]) )
        {
          P4EST_ASSERT(k<2);
          if (cells_around[j].is_hanging)
          {
            jneig[2*k]       = cells_around[j].is.hanging.quadid[0];
            jneig[2*k+1]     = cells_around[j].is.hanging.quadid[1];
            jneig_jedge[2*k] = cells_around[j].edge;
          }
          else
          {
            jneig[2*k]       = cells_around[j].is.full.quadid;
            jneig_jedge[2*k] = cells_around[j].edge;
          } 
          
          if (cells_around[i].is_hanging && cells_around[j].is_hanging) 
          {
              // i side
              quad_to_quad_by_edge[ 12*ineig[2*k]   + ineig_iedge[2*k] ] = jneig[2*k];
              quad_to_edge        [ 12*ineig[2*k]   + ineig_iedge[2*k] ] = jneig_jedge[2*k];
              quad_to_quad_by_edge[ 12*ineig[2*k+1] + ineig_iedge[2*k] ] = jneig[2*k+1];
              quad_to_edge        [ 12*ineig[2*k+1] + ineig_iedge[2*k] ] = jneig_jedge[2*k];
              
              //j side
              quad_to_quad_by_edge[ 12*jneig[2*k]   + jneig_jedge[2*k] ] = ineig[2*k];
              quad_to_edge        [ 12*jneig[2*k]   + jneig_jedge[2*k] ] = ineig_iedge[2*k];
              quad_to_quad_by_edge[ 12*jneig[2*k+1] + jneig_jedge[2*k] ] = ineig[2*k+1];
              quad_to_edge        [ 12*jneig[2*k+1] + jneig_jedge[2*k] ] = ineig_iedge[2*k];
          }
          else if (! cells_around[i].is_hanging && cells_around[j].is_hanging) 
          {
              // i side
              quad_to_quad_by_edge[ 12*ineig[2*k] + ineig_iedge[2*k] ] = ineig[2*k];
              quad_to_edge        [ 12*ineig[2*k] + ineig_iedge[2*k] ] = jneig_jedge[2*k]-24;
              quad_to_half_by_edge[ 2*ineig[2*k]                     ] = jneig[2*k];
              quad_to_half_by_edge[ 2*ineig[2*k] +                 1 ] = jneig[2*k+1];
              
              //j side
              quad_to_quad_by_edge[ 12*jneig[2*k]   + jneig_jedge[2*k] ] = ineig[2*k];
              quad_to_edge        [ 12*jneig[2*k]   + jneig_jedge[2*k] ] = 24+ineig_iedge[2*k];
              quad_to_quad_by_edge[ 12*jneig[2*k+1] + jneig_jedge[2*k] ] = ineig[2*k];
              quad_to_edge        [ 12*jneig[2*k+1] + jneig_jedge[2*k] ] = 48+ineig_iedge[2*k];
          }
          else if (cells_around[i].is_hanging && !cells_around[j].is_hanging) 
          {
              // i side
              quad_to_quad_by_edge[ 12*ineig[2*k]   + ineig_iedge[2*k] ] = jneig[2*k];
              quad_to_edge        [ 12*ineig[2*k]   + ineig_iedge[2*k] ] = 24+jneig_jedge[2*k];
              quad_to_quad_by_edge[ 12*ineig[2*k+1] + ineig_iedge[2*k] ] = jneig[2*k];
              quad_to_edge        [ 12*ineig[2*k+1] + ineig_iedge[2*k] ] = 48+jneig_jedge[2*k];
             
              //j side
              quad_to_quad_by_edge[ 12*jneig[2*k] + jneig_jedge[2*k] ] = jneig[2*k];
              quad_to_edge        [ 12*jneig[2*k] + jneig_jedge[2*k] ] = ineig_iedge[2*k]-24;
              quad_to_half_by_edge[ 2*jneig[2*k]                     ] = ineig[2*k];
              quad_to_half_by_edge[ 2*jneig[2*k] +                 1 ] = ineig[2*k+1];
          }   
          else // !cells_around[i].is_hanging && !cells_around[j].is_hanging
          {
              quad_to_quad_by_edge[ 12*ineig[2*k] + ineig_iedge[2*k] ] = jneig[2*k];
              quad_to_edge[ 12*ineig[2*k] + ineig_iedge[2*k] ] = jneig_jedge[2*k];
              quad_to_quad_by_edge[ 12*jneig[2*k] + jneig_jedge[2*k] ] = ineig[2*k];
              quad_to_edge[ 12*jneig[2*k] + jneig_jedge[2*k] ] = ineig_iedge[2*k];
          }
          k++;
          if (k==2) return;
        }
    }
  }
}

int refine_callback_2d(p4est_t * p4est,
                    p4est_topidx_t which_tree,
                    p4est_quadrant_t * quadrant)
{
    P4EST_ASSERT(which_tree == 0);
    return (*((int *)quadrant->p.user_data) ==  FEMPAR_refinement_flag);
}

int refine_callback_3d(p8est_t * p8est,
                    p4est_topidx_t which_tree,
                    p8est_quadrant_t * quadrant)
{
    P4EST_ASSERT(which_tree == 0);
    return (*((int *)quadrant->p.user_data) ==  FEMPAR_refinement_flag);
}


void  refine_replace_callback_2d (p4est_t * p4est,
                               p4est_topidx_t which_tree,
                               int num_outgoing,
                               p4est_quadrant_t * outgoing[],
                               int num_incoming,
                               p4est_quadrant_t * incoming[])
 {
    int quadrant_index;
    int *quadrant_data;
    P4EST_ASSERT(which_tree   == 0);
    P4EST_ASSERT(num_outgoing == 1);
    P4EST_ASSERT(num_incoming == P4EST_CHILDREN);
    for (quadrant_index=0; quadrant_index < P4EST_CHILDREN; quadrant_index++)
    {
      quadrant_data = (int *) incoming[quadrant_index]->p.user_data;
      *quadrant_data = FEMPAR_do_nothing_flag;
    }
 }

void  refine_replace_callback_3d (p8est_t * p8est,
                               p4est_topidx_t which_tree,
                               int num_outgoing,
                               p8est_quadrant_t * outgoing[],
                               int num_incoming,
                               p8est_quadrant_t * incoming[])
 {
    int quadrant_index;
    int *quadrant_data;
    P4EST_ASSERT(which_tree   == 0);
    P4EST_ASSERT(num_outgoing == 1);
    P4EST_ASSERT(num_incoming == P8EST_CHILDREN);
    for (quadrant_index=0; quadrant_index < P8EST_CHILDREN; quadrant_index++)
    {
      quadrant_data = (int *) incoming[quadrant_index]->p.user_data;
      *quadrant_data = FEMPAR_do_nothing_flag;
    }
 }

void F90_p4est_refine( p4est_t * p4est )
{
    p4est_refine_ext(p4est, 0, -1, refine_callback_2d, NULL, refine_replace_callback_2d);
}

void F90_p8est_refine( p8est_t * p8est )
{
    p8est_refine_ext(p8est, 0, -1, refine_callback_3d, NULL, refine_replace_callback_3d);
}

int coarsen_callback_2d (p4est_t * p4est,
                      p4est_topidx_t which_tree,
                      p4est_quadrant_t * quadrants[])
{
    int quadrant_index;
    int coarsen;
    P4EST_ASSERT(which_tree == 0);
    
    coarsen = 1;
    for (quadrant_index=0; quadrant_index < P4EST_CHILDREN; quadrant_index++)
    {
      coarsen = (*((int *)(quadrants[quadrant_index]->p.user_data)) ==  FEMPAR_coarsening_flag);
      if (!coarsen) return coarsen;
    }
    return coarsen;
}

int coarsen_callback_3d (p8est_t * p8est,
                      p4est_topidx_t which_tree,
                      p8est_quadrant_t * quadrants[])
{
    int quadrant_index;
    int coarsen;
    P4EST_ASSERT(which_tree == 0);
    
    coarsen = 1;
    for (quadrant_index=0; quadrant_index < P8EST_CHILDREN; quadrant_index++)
    {
      coarsen = (*((int *)(quadrants[quadrant_index]->p.user_data)) ==  FEMPAR_coarsening_flag);
      if (!coarsen) return coarsen;
    }
    return coarsen;
}

void F90_p4est_coarsen( p4est_t * p4est )
{
    p4est_coarsen(p4est, 0, coarsen_callback_2d, NULL);
}

void F90_p8est_coarsen( p8est_t * p8est )
{
    p8est_coarsen(p8est, 0, coarsen_callback_3d, NULL);
}

void F90_p4est_copy( p4est_t * p4est_input, p4est_t ** p4est_output )
{
   F90_p4est_destroy(p4est_output);
   *p4est_output = p4est_copy(p4est_input,0);
}

void F90_p8est_copy( p8est_t * p8est_input, p8est_t ** p8est_output )
{
   F90_p8est_destroy(p8est_output);
   *p8est_output = p8est_copy(p8est_input,0);
}

void F90_p4est_balance( p4est_t * p4est )
{
  p4est_balance(p4est, P4EST_CONNECT_FULL, NULL);
}

void F90_p8est_balance( p8est_t * p8est )
{
  p8est_balance(p8est, P8EST_CONNECT_FULL, NULL);
}

void F90_p4est_update_refinement_and_coarsening_flags(p4est_t * p4est_old, p4est_t * p4est_new)
{
    p4est_tree_t       *tree_old;
    p4est_quadrant_t   *q_old;
    sc_array_t         *quadrants_old;
    int                old_quadrant_index;
    
    p4est_tree_t       *tree_new;
    p4est_quadrant_t   *q_new;
    sc_array_t         *quadrants_new;
    int                new_quadrant_index;
    
    int * user_pointer;
   
    P4EST_ASSERT(p4est_old->user_pointer == p4est_new->user_pointer);
    
    user_pointer = (int *) p4est_old->user_pointer;
    
    // Extract references to the first (and uniquely allowed) trees
    tree_old = p4est_tree_array_index (p4est_old->trees,0);
    tree_new = p4est_tree_array_index (p4est_new->trees,0);
    
    quadrants_old = &(tree_old->quadrants);
    quadrants_new = &(tree_new->quadrants);
    
    new_quadrant_index = 0;
    for (old_quadrant_index=0; old_quadrant_index < quadrants_old->elem_count;)
    {
       q_old = p4est_quadrant_array_index(quadrants_old, old_quadrant_index);
       q_new = p4est_quadrant_array_index(quadrants_new, new_quadrant_index);
       if ( p4est_quadrant_compare(q_old,q_new) == 0 ) //q_old was not refined nor coarsened
       {
           user_pointer[old_quadrant_index] = FEMPAR_do_nothing_flag;
           old_quadrant_index++;
           new_quadrant_index++;
       }
       else if ( p4est_quadrant_is_parent(q_old,q_new)  )  //q_old was refined
       { 
           user_pointer[old_quadrant_index] = FEMPAR_refinement_flag;
           old_quadrant_index++;
           new_quadrant_index = new_quadrant_index + P4EST_CHILDREN;
       }
       else if ( p4est_quadrant_is_parent(q_new,q_old) ) //q_old and its siblings were coarsened 
       {
           for (int i=0; i < P4EST_CHILDREN; i++)
           {
               user_pointer[old_quadrant_index] = FEMPAR_coarsening_flag;
               old_quadrant_index++;
           }
           new_quadrant_index++;
       }
       else
       {
         P4EST_ASSERT(0);
       }
    }
}

void F90_p8est_update_refinement_and_coarsening_flags(p8est_t * p8est_old, p8est_t * p8est_new)
{
    p8est_tree_t       *tree_old;
    p8est_quadrant_t   *q_old;
    sc_array_t         *quadrants_old;
    int                old_quadrant_index;
    
    p8est_tree_t       *tree_new;
    p8est_quadrant_t   *q_new;
    sc_array_t         *quadrants_new;
    int                new_quadrant_index;
    
    int * user_pointer;
   
    P4EST_ASSERT(p8est_old->user_pointer == p8est_new->user_pointer);
    
    user_pointer = (int *) p8est_old->user_pointer;
    
    // Extract references to the first (and uniquely allowed) trees
    tree_old = p8est_tree_array_index (p8est_old->trees,0);
    tree_new = p8est_tree_array_index (p8est_new->trees,0);
    
    quadrants_old = &(tree_old->quadrants);
    quadrants_new = &(tree_new->quadrants);
    
    new_quadrant_index = 0;
    for (old_quadrant_index=0; old_quadrant_index < quadrants_old->elem_count;)
    {
       q_old = p8est_quadrant_array_index(quadrants_old, old_quadrant_index);
       q_new = p8est_quadrant_array_index(quadrants_new, new_quadrant_index);
       if ( p8est_quadrant_compare(q_old,q_new) == 0 ) //q_old was not refined nor coarsened
       {
           user_pointer[old_quadrant_index] = FEMPAR_do_nothing_flag;
           old_quadrant_index++;
           new_quadrant_index++;
       }
       else if ( p8est_quadrant_is_parent(q_old,q_new)  )  //q_old was refined
       { 
           user_pointer[old_quadrant_index] = FEMPAR_refinement_flag;
           old_quadrant_index++;
           new_quadrant_index = new_quadrant_index + P8EST_CHILDREN;
       }
       else if ( p8est_quadrant_is_parent(q_new,q_old) ) //q_old and its siblings were coarsened 
       {
           for (int i=0; i < P8EST_CHILDREN; i++)
           {
               user_pointer[old_quadrant_index] = FEMPAR_coarsening_flag;
               old_quadrant_index++;
           }
           new_quadrant_index++;
       }
       else
       {
         P4EST_ASSERT(0);
       }
    }

}

void F90_p4est_get_quadrant_vertex_coordinates(p4est_connectivity_t * connectivity,
                                               p4est_topidx_t treeid,
                                               p4est_qcoord_t x,
                                               p4est_qcoord_t y, 
                                               int8_t level,
                                               int   corner,
                                               double vxyz[3])
{
    P4EST_ASSERT(treeid == 0);
    p4est_quadrant_t myself;
    p4est_quadrant_t neighbour;
    
    // Create myself
    myself.x     = x;
    myself.y     = y; 
    myself.level = level;
    
    if ( corner == 0 ) {
          neighbour = myself;
      }      
    else if ( corner == 1 ) {
       p4est_quadrant_face_neighbor(&myself, 
                                    corner, 
                                    &neighbour);
    }
    else if ( corner == 2 ) { 
       p4est_quadrant_face_neighbor(&myself, 
                                    corner+1, 
                                    &neighbour); 
   }
    else if ( corner == 3 ) {   
       p4est_quadrant_corner_neighbor(&myself, 
                                    corner, 
                                    &neighbour); 
   }
   
    // Extract numerical coordinates of lower_left corner of my corner neighbour
    p4est_qcoord_to_vertex (connectivity, 
                            treeid, 
                            neighbour.x, 
                            neighbour.y, 
                            vxyz);
    vxyz[2] = 0.0;
}

void F90_p8est_get_quadrant_vertex_coordinates(p8est_connectivity_t * connectivity,
                                               p4est_topidx_t treeid,
                                               p4est_qcoord_t x,
                                               p4est_qcoord_t y, 
                                               p4est_qcoord_t z, 
                                               int8_t level,
                                               int   corner,
                                               double vxyz[3])
{
    P4EST_ASSERT(treeid == 0);
    p8est_quadrant_t myself;
    p8est_quadrant_t neighbour;
    
    // Create myself
    myself.x     = x;
    myself.y     = y; 
    myself.z     = z; 
    myself.level = level;
    
    if ( corner == 0 ) {
          neighbour = myself;
      }      
    else if ( corner == 1 ) {
       p8est_quadrant_face_neighbor(&myself, 
                                    1, 
                                    &neighbour);
    }
    else if ( corner == 2 ) { 
       p8est_quadrant_face_neighbor(&myself, 
                                    3, 
                                    &neighbour); 
   }
    else if ( corner == 3 ) {   
       p8est_quadrant_edge_neighbor(&myself, 
                                    11, 
                                    &neighbour); 
    }
    else if ( corner == 4 ) {   
       p8est_quadrant_face_neighbor(&myself, 
                                    5, 
                                    &neighbour); 
   }
    else if ( corner == 5 ) {   
       p8est_quadrant_edge_neighbor(&myself, 
                                    7, 
                                    &neighbour); 
   }
    else if ( corner == 6 ) {   
       p8est_quadrant_edge_neighbor(&myself, 
                                    3, 
                                    &neighbour); 
   }
    else if ( corner == 7 ) {   
       p8est_quadrant_corner_neighbor(&myself, 
                                    7, 
                                    &neighbour); 
   }
   
    // Extract numerical coordinates of lower_left corner of my corner neighbour
    p8est_qcoord_to_vertex (connectivity, 
                            treeid, 
                            neighbour.x, 
                            neighbour.y, 
                            neighbour.z, 
                            vxyz);
}

int F90_p4est_is_ancestor ( p4est_qcoord_t q1_x,
                            p4est_qcoord_t q1_y,
                            int8_t q1_level,
                            p4est_qcoord_t q2_x,
                            p4est_qcoord_t q2_y,
                            int8_t q2_level )
{
    p4est_quadrant_t q1;
    p4est_quadrant_t q2;
    
    q1.x = q1_x;
    q1.y = q1_y;
    q1.level = q1_level;
    
    q2.x = q2_x;
    q2.y = q2_y;
    q2.level = q2_level;
    
    return p4est_quadrant_is_ancestor ( &q1, &q2 );
    
}

void F90_p4est_quadrant_set_morton ( int level,
                                     int64_t id,
                                     p4est_qcoord_t *q_x,
                                     p4est_qcoord_t *q_y )
{
    p4est_quadrant_t q;
    p4est_quadrant_set_morton( &q, level, id );
    *q_x = q.x;
    *q_y = q.y;
}

//void p4_savemesh ( char    filename[],
//                   p4est_t *p4est)
//{
//  p4est_t              *p4est2;
//  p4est_connectivity_t *conn2 = NULL;
//  int ip,ic;
//  
//  p4est_save(filename,p4est,0);
//  p4est2=p4est_load_ext(filename,mpicomm,0,0,1,0,NULL,&conn2);
//  // TODO: optional check
//  ic = p4est_connectivity_is_equal(p4est->connectivity,conn2);
//  ip = p4est_is_equal(p4est,p4est2,0);
//  printf("Conn, p4est %i %i \n",ic,ip);
//  p4est_destroy(p4est2);
//  p4est_connectivity_destroy(conn2);
//}

#endif
