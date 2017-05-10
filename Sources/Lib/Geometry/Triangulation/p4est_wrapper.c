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
//#include <sc_io.h>
//#include <p8est_algorithms.h>
//#include <p8est_bits.h>
//#include <p8est_vtk.h>
//#include <p8est_mesh.h>
//#include <p8est_extended.h>
//#include "p4est_HO_geometry.h"
//#include <p4est/p4est_wrapper.h>
//// 3D mode; must be included AFTER every other file
//#include <p4est_to_p8est.h>

//p8est_connectivity_t *global_conn;

//static sc_MPI_Comm           mpicomm;

//// Init p4est 
//void p4_initvars()
//{

//  /* Initialize MPI; see sc_mpi.h.
//   * If configure --enable-mpi is given these are true MPI calls.
//   * Else these are dummy functions that simulate a single-processor run. */
//  mpicomm = sc_MPI_COMM_WORLD;

//  /* These functions are optional.  If called they store the MPI rank as a
//   * static variable so subsequent global p4est log messages are only issued
//   * from processor zero.  Here we turn off most of the logging; see sc.h. */
//  sc_init (mpicomm, 1, 1, NULL, SC_LP_ESSENTIAL);
//  p4est_init (NULL, SC_LP_PRODUCTION);
//}


//// Build p4est data structures by loading existing connectivity file
//void p4_loadmesh(char    filename[],
//                 p4est_t **p4est_out )
//{
//  p4est_t              *p4est;
//  p4est_connectivity_t *conn = NULL;

//  p4est=p4est_load_ext(filename,mpicomm,0,0,1,0,NULL,&conn);
//  P4EST_ASSERT (p4est_is_valid (p4est));
//  P4EST_ASSERT (p4est_connectivity_is_valid (conn));
//  *p4est_out=p4est;
//}



//// Build p4est data structures with existing connectivity from HDF5 mesh
//void p4_connectivity_treevertex (p4est_topidx_t num_vertices,
//                                 p4est_topidx_t num_trees,
//                                 double         *vertices,
//                                 p4est_topidx_t *tree_to_vertex,
//                                 p4est_topidx_t num_periodics,
//                                 p4est_topidx_t *join_faces,
//                                 p4est_connectivity_t **conn_out )
//{
//  p4est_topidx_t        tree;
//  int                   face,i;
//  p4est_connectivity_t *conn = NULL;


//  P4EST_GLOBAL_PRODUCTIONF ("creating connectivity from tree to vertex only...%d %d \n",num_vertices,num_trees);

//  conn = p4est_connectivity_new (num_vertices, num_trees,
//                                 0, 0,
//                                 0, 0);
//  for (i=0; i < 3*num_vertices; ++i){
//    conn->vertices[i] = vertices[i];
//  }
//  for (i=0; i < 8*num_trees; ++i){
//    conn->tree_to_vertex[i] = tree_to_vertex[i];
//  }

//  /*
//   * Fill tree_to_tree and tree_to_face to make sure we have a valid
//   * connectivity.
//   */

//  for (tree = 0; tree < conn->num_trees; ++tree) {
//    for (face = 0; face < P4EST_FACES; ++face) {
//      conn->tree_to_tree[P4EST_FACES * tree + face] = tree;
//      conn->tree_to_face[P4EST_FACES * tree + face] = face;
//    }
//  }

//  P4EST_ASSERT (p4est_connectivity_is_valid (conn));

//  p4est_connectivity_complete (conn);

//  P4EST_ASSERT (p4est_connectivity_is_valid (conn));
//  // Join Faces
//  if ( num_periodics>0) {
//    for (i=0; i < num_periodics; ++i){
//      p8est_connectivity_join_faces(conn,join_faces[5*i],join_faces[5*i+1], 
//                                         join_faces[5*i+2],join_faces[5*i+3],
//                                         join_faces[5*i+4]);
//    } 
//  }

//  P4EST_ASSERT (p4est_connectivity_is_valid (conn));

//  P4EST_GLOBAL_PRODUCTIONF
//    ("New connectivity with %lld trees and %lld vertices\n",
//     (long long) conn->num_trees, (long long) conn->num_vertices);

//  *conn_out=conn;
//}

//void p4_build_p4est ( p4est_connectivity_t *conn,
//                      p4est_t              **p4est_out,
//                      p4est_geometry_t     **geom_out)
//{
//  p4est_t* p4est;
//  p4est_geometry_t   *geom;

//  fflush(stdout);
//  P4EST_ASSERT (p4est_connectivity_is_valid (conn));
//  /* Create a forest that is not refined; it consists of the root octant. */
//  p4est = p4est_new (mpicomm, conn, 0, NULL, NULL);

//  geom = P4EST_ALLOC_ZERO (p4est_geometry_t, 1);
//  geom->name = "hopest_readfromhdf5";
//  geom->X = p4_geometry_X;

//  *p4est_out=p4est;
//  *geom_out=geom;
//}


//void p4_build_bcs(p4est_t        *p4est,
//                  p4est_topidx_t num_trees,
//                  int32_t        *bcelemmap)
//{
//  int itree,iside;

//  p8est_connectivity_t *conn=p4est->connectivity;
//  P4EST_ASSERT (p4est_connectivity_is_valid (conn));

//  P4EST_ASSERT (p4est->trees->elem_count == num_trees);
//  p4est_connectivity_set_attr(conn,6*sizeof(int32_t));
//  P4EST_ASSERT (p4est_connectivity_is_valid (conn));
//  
//  for(itree=0; itree<num_trees; itree++) {
//    for(iside=0; iside<6; iside++) {
//      ((int32_t*) conn->tree_to_attr)[itree*6+iside]=bcelemmap[itree*6+iside];
//    }
//  }
//}


//void p4_build_mesh(p4est_t  *p4est,
//                   p4est_mesh_t  **mesh_out )
//{
//  p4est_mesh_t       *mesh;
//  p4est_ghost_t      *ghost;

//  /* create ghost layer and mesh */
//  ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
//  mesh = p4est_mesh_new_ext (p4est, ghost, 1,1,P4EST_CONNECT_FULL);

//  //return mesh as pointer adress;
//  *mesh_out=(p4est_mesh_t *)mesh;
//}

//void p4_get_bcs(p4est_t        *p4est,
//                int32_t        **bcelemmap)
//{
//  
//  p8est_connectivity_t *conn=p4est->connectivity;
//  *bcelemmap=(int32_t*) conn->tree_to_attr;
//}


//void p4_get_mesh_info ( p4est_t        *p4est,
//                        p4est_mesh_t   *mesh,
//                        p4est_locidx_t *local_num_quadrants,
//                        p4est_gloidx_t *global_num_quadrants,
//                        p4est_gloidx_t *global_first_quadrant,
//                        p4est_locidx_t *num_half_faces,
//                        p4est_topidx_t *num_trees )
//{
//  *local_num_quadrants   = p4est->local_num_quadrants;
//  *global_num_quadrants  = p4est->global_num_quadrants;
//  *global_first_quadrant = p4est->global_first_quadrant[p4est->mpirank];
//  /* big face with 4 small neighbours */
//  *num_half_faces = (p4est_locidx_t) mesh->quad_to_half->elem_count;
//  *num_trees = (p4est_topidx_t) p4est->trees->elem_count;
//  SC_CHECK_ABORTF (mesh->local_num_quadrants == p4est->global_num_quadrants,
//                   "Global quads %lld and local quads %d mismatch ! ",
//                   (long long) p4est->global_num_quadrants, mesh->local_num_quadrants );
//  /* TODO:
//   * In parallel, local and global num_quadrants WILL be different.
//   * This check will then die.
//   */
//  SC_CHECK_ABORTF ((p4est_gloidx_t) mesh->local_num_quadrants
//                   == p4est->global_num_quadrants,
//                   "Global quads %lld and local quads %d mismatch ! ",
//                   (long long) p4est->global_num_quadrants,
//                   mesh->local_num_quadrants );
//}

//void p4_get_quadrants( p4est_t       *p4est,
//                       p4est_mesh_t   *mesh,
//                       p4est_locidx_t local_num_quadrants,
//                       p4est_locidx_t num_half_faces,
//                       p4est_qcoord_t  *intsize,
//                       p4est_topidx_t **quad_to_tree,
//                       p4est_locidx_t **quad_to_quad,
//                       int8_t         **quad_to_face, 
//                       p4est_locidx_t **quad_to_half, 
//                       p4est_qcoord_t *quadcoords,
//                       int8_t         *quadlevel ) 
//{
//  int iquad,iquadloc;
//  p8est_tree_t       *tree;
//  p8est_quadrant_t   *q;
//  sc_array_t         *quadrants;

//  *intsize = P4EST_ROOT_LEN;

//  P4EST_ASSERT (local_num_quadrants == p4est->local_num_quadrants);
//  for (iquad = 0; iquad < mesh->local_num_quadrants; iquad++) {
//    tree = p8est_tree_array_index (p4est->trees,mesh->quad_to_tree[iquad]);
//    quadrants = &(tree->quadrants);
//    iquadloc = iquad - tree->quadrants_offset;
//    q = p8est_quadrant_array_index(quadrants, iquadloc);
//    quadlevel [iquad    ] = q->level;
//    quadcoords[iquad*3  ] = q->x;
//    quadcoords[iquad*3+1] = q->y;
//    quadcoords[iquad*3+2] = q->z;
//  }

//  *quad_to_tree=mesh->quad_to_tree;
//  *quad_to_quad=mesh->quad_to_quad;
//  *quad_to_face=mesh->quad_to_face;

//  *quad_to_half = NULL;
//  if(num_half_faces>0) *quad_to_half = (p4est_locidx_t *) mesh->quad_to_half->array;
//}

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


