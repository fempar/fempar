/* Module file taken "AS-IS" and adapted to FEMPAR needs from p4est_wrapper.h
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

//#include <p8est_mesh.h>
//#include <p8est_geometry.h>

//#ifdef __cplusplus
//extern              "C"         /* prevent C++ name mangling */
//{
//#if 0
//}
//#endif
//#endif

//void p4_initvars();

//void p4_loadmesh(char    filename[],
//                 p8est_t **p4est_out );

//void p4_connectivity_treevertex (p4est_topidx_t num_vertices,
//                                 p4est_topidx_t num_trees,
//                                 double         *vertices,
//                                 p4est_topidx_t *tree_to_vertex,
//                                 p4est_topidx_t num_periodics,
//                                 p4est_topidx_t *join_faces,
//                                 p8est_connectivity_t        **conn_out );
//void p4_build_p4est ( p8est_connectivity_t *conn,
//                      p8est_t              **p4est_out,
//                      p8est_geometry_t     **geom_out);

//void p4_build_bcs(p8est_t        *p4est,
//                  p4est_topidx_t num_trees,
//                  int32_t        *bcelemmap);

//void p4_get_bcs(p8est_t        *p4est,
//                int32_t        **bcelemmap);

//void p4_build_mesh(p8est_t  *p4est,
//                   p8est_mesh_t  **mesh_out );

//void p4_get_mesh_info ( p8est_t        *p4est,
//                        p8est_mesh_t   *mesh,
//                        p4est_locidx_t *local_num_quadrants,
//                        p4est_gloidx_t *global_num_quadrants,
//                        p4est_gloidx_t *global_first_quadrant,
//                        p4est_locidx_t *num_half_faces,
//                        p4est_topidx_t *num_trees );

//void p4_get_quadrants( p8est_t       *p4est,
//                       p8est_mesh_t   *mesh,
//                       p4est_locidx_t local_num_quadrants,
//                       p4est_locidx_t   num_half_faces,
//                       p4est_qcoord_t  *intsize,
//                       p4est_topidx_t **quad_to_tree,
//                       p4est_locidx_t **quad_to_quad,
//                       int8_t         **quad_to_face, 
//                       p4est_locidx_t **quad_to_half, 
//                       p4est_qcoord_t *quadcoords,
//                       int8_t         *quadlevel );

//void p4_savemesh ( char    filename[],
//                   p8est_t *p4est);
//                   
//void p4est_save_all ( char    filename[],
//                      p8est_t *p4est);

//#ifdef __cplusplus
//#if 0
//{
//#endif
//}
//#endif
