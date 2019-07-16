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
module p4est_triangulation_parameters_names
  use types_names
  implicit none
 
  ! Numbers designating the level of logging output.
  !
  ! Priorities TRACE to VERBOSE are appropriate when all parallel processes
  ! contribute log messages.  INFO and above must not clutter the output of
  ! large parallel runs.  STATISTICS can be used for important measurements.
  ! PRODUCTION is meant for rudimentary information on the program flow.
  ! ESSENTIAL can be used for one-time messages, say at program startup.
  ! 
  integer(ip), parameter :: SC_LP_DEFAULT        = -1 !**< this selects the SC default threshold */
  integer(ip), parameter :: SC_LP_ALWAYS         = 0  !**< this will log everything */  
  integer(ip), parameter :: SC_LP_TRACE          = 1  !**< this will prefix file and line number */
  integer(ip), parameter :: SC_LP_DEBUG          = 2  !**< any information on the internal state */
  integer(ip), parameter :: SC_LP_VERBOSE        = 3  !**< information on conditions, decisions */
  integer(ip), parameter :: SC_LP_INFO           = 4  !**< the main things a function is doing */
  integer(ip), parameter :: SC_LP_STATISTICS     = 5  !**< important for consistency/performance */
  integer(ip), parameter :: SC_LP_PRODUCTION     = 6  !**< a few lines for a major api function */
  integer(ip), parameter :: SC_LP_ESSENTIAL      = 7  !**< this logs a few lines max per program */
  integer(ip), parameter :: SC_LP_ERROR          = 8  !**< this logs errors only */
  integer(ip), parameter :: SC_LP_SILENT         = 9  !**< this never logs anything */
  integer(ip), parameter :: FEMPAR_SC_LP_DEFAULT = SC_LP_SILENT
  
  ! Parameter handling
  character(len=*), parameter :: p4est_triang_num_dims_key                     = 'P4EST_TRIANG_NUM_DIMS'
  character(len=*), parameter :: p4est_triang_domain_limits_key                = 'P4EST_TRIANG_DOMAIN_LIMITS'
  character(len=*), parameter :: p4est_triang_geometry_interpolation_order_key = 'P4EST_TRIANG_GEOMETRY_INTERPOLATION_ORDER'
  character(len=*), parameter :: p4est_triang_log_level_key                    = 'P4EST_TRIANG_LOG_LEVEL'
  character(len=*), parameter :: p4est_triang_2_1_k_balance_key                = 'P4EST_TRIANG_2_1_K_BALANCE'
  character(len=*), parameter :: p4est_triang_k_ghost_cells_key                = 'P4EST_TRIANG_K_GHOST_CELLS'

  character(len=*), parameter :: p4est_triang_num_dims_cla_name                     = '--'//p4est_triang_num_dims_key
  character(len=*), parameter :: p4est_triang_domain_limits_cla_name                = '--'//p4est_triang_domain_limits_key
  character(len=*), parameter :: p4est_triang_geometry_interpolation_order_cla_name = '--'//p4est_triang_geometry_interpolation_order_key
  character(len=*), parameter :: p4est_triang_log_level_cla_name                    = '--'//p4est_triang_log_level_key
  character(len=*), parameter :: p4est_triang_2_1_k_balance_cla_name                = '--'//p4est_triang_2_1_k_balance_key
  character(len=*), parameter :: p4est_triang_k_ghost_cells_cla_name                = '--'//p4est_triang_k_ghost_cells_key
  
  character(len=*), parameter :: p4est_triang_log_level_choices         = '-1,0,1,2,3,4,5,6,7,8,9'
  character(len=*), parameter :: p4est_triang_num_dims_cla_choices      = '2,3'
  character(len=*), parameter :: p4est_triang_2_1_k_balance_cla_choices = '0,1,2'
  character(len=*), parameter :: p4est_triang_k_ghost_cells_cla_choices = '0,1,2'

  character(len=*), parameter, public ::p4est_triang_num_dims_cla_help      = 'p4est triangulation number of space dimensions'
  character(len=*), parameter, public ::p4est_triang_domain_limits_cla_help = 'p4est triangulation domain interval per dimension'
  character(len=*), parameter, public ::p4est_triang_geometry_interpolation_order_cla_help = 'p4est triangulation polynomial order of the Lagrangian FE space used to discretize the geometry of the domain'
  character(len=*), parameter, public ::p4est_triang_2_1_k_balance_cla_help = 'Value of k for 2:1 k-balanced forest-of-octrees (use with care, at present,' // BRK_LINE // &
                                                                              'only k={0,1} supported/tested)'
  character(len=*), parameter, public ::p4est_triang_k_ghost_cells_cla_help = 'Value of k for the k-ghost cells set of each processor'                   //BRK_LINE // &
                                                                              '(k=0 works for any FE space; k>0 should work depending on the FE space,' // BRK_LINE // &
                                                                              'although NOT tested, use with care)'

  character(len=*), parameter, public :: p4est_triang_log_level_help    = "P4EST library level of logging output"      // BRK_LINE // & 
                   BULLET_FLAP_HELP_MESSAGE // "-1: SC_LP_DEFAULT    (this selects the SC default threshold)" // BRK_LINE // & 
                   BULLET_FLAP_HELP_MESSAGE // " 0: SC_LP_ALWAYS     (this will log everything)"              // BRK_LINE // &
                   BULLET_FLAP_HELP_MESSAGE // " 1: SC_LP_TRACE      (this will prefix file and line number)" // BRK_LINE // &
                   BULLET_FLAP_HELP_MESSAGE // " 2: SC_LP_DEBUG      (any information on the internal state)" // BRK_LINE // &
                   BULLET_FLAP_HELP_MESSAGE // " 3: SC_LP_VERBOSE    (information on conditions, decisions)"  // BRK_LINE // &
                   BULLET_FLAP_HELP_MESSAGE // " 4: SC_LP_INFO       (the main things a function is doing)"   // BRK_LINE // &
                   BULLET_FLAP_HELP_MESSAGE // " 5: SC_LP_STATISTICS (important for consistency/performance)" // BRK_LINE // &
                   BULLET_FLAP_HELP_MESSAGE // " 6: SC_LP_PRODUCTION (a few lines for a major api function)"  // BRK_LINE // &
                   BULLET_FLAP_HELP_MESSAGE // " 7: SC_LP_ESSENTIAL  (this logs a few lines max per program)" // BRK_LINE // &
                   BULLET_FLAP_HELP_MESSAGE // " 8: SC_LP_ERROR      (this logs errors only)"                 // BRK_LINE // &
                   BULLET_FLAP_HELP_MESSAGE // " 9: SC_LP_SILENT     (this never logs anything)"
  
  ! Unit square/cube by default; all points on and in the domain have space coordinates >= 0.0_rp
  real(rp)   , parameter :: default_p4est_triang_domain_limits (*)            = [0.0_rp,1.0_rp,0.0_rp,1.0_rp,0.0_rp,1.0_rp]
  integer(ip), parameter :: default_p4est_triang_geometry_interpolation_order = 1
  integer(ip), parameter :: default_p4est_triang_log_level                    = FEMPAR_SC_LP_DEFAULT
  integer(ip), parameter :: default_p4est_triang_num_dims                     = 2 
  integer(ip), parameter :: default_p4est_triang_2_1_k_balance                = 0
  integer(ip), parameter :: default_p4est_triang_k_ghost_cells                = 0
    
  ! For 2D
  integer(ip), parameter :: NUM_SUBCELLS_IN_TOUCH_FACE_2D = 2
  integer(ip), parameter :: NUM_CORNERS_2D                = 4
  integer(ip), parameter :: NUM_FACES_2D                  = 4
  integer(ip), parameter :: NUM_SUBFACES_FACE_2D          = 2
  integer(ip), parameter :: NUM_FACE_CORNERS_2D           = 2
  integer(ip), parameter :: NUM_FACES_AT_CORNER_2D        = 2
  integer(ip), parameter :: NUM_VEFS_2D                   = NUM_CORNERS_2D+NUM_FACES_2D
  integer(ip), parameter :: MAX_NUM_CELLS_AROUND_2D       = 4 

  integer(ip), target :: P4EST_FACE_CORNERS_2D(NUM_FACE_CORNERS_2D,NUM_FACES_2D) = & 
                                                  reshape([1, 3,&
                                                           2, 4,&  
                                                           1, 2,&
                                                           3, 4], [NUM_FACE_CORNERS_2D,NUM_FACES_2D])
                                                  
  integer(ip), target :: P4EST_FACES_AT_CORNER_2D(NUM_FACES_AT_CORNER_2D, NUM_CORNERS_2D) = & 
                                                    reshape([1, 3,&
                                                             2, 3,&  
                                                             1, 4,&
                                                             2, 4], [NUM_FACES_AT_CORNER_2D, NUM_CORNERS_2D])
                                                    
  integer(ip), target :: FEMPAR_SUBCELLS_IN_TOUCH_FACE_2D(NUM_SUBCELLS_IN_TOUCH_FACE_2D,NUM_FACES_2D) = &
                                                    reshape([1, 2,&
                                                             3, 4,&  
                                                             1, 3,&
                                                             2, 4], [NUM_SUBCELLS_IN_TOUCH_FACE_2D, NUM_FACES_2D])
                                                    
  integer(ip), target :: P4EST_SUBCELLS_IN_TOUCH_FACE_2D(NUM_SUBCELLS_IN_TOUCH_FACE_2D,NUM_FACES_2D) = &
                                                    reshape([1, 3,&
                                                             2, 4,&  
                                                             1, 2,&
                                                             3, 4], [NUM_SUBCELLS_IN_TOUCH_FACE_2D, NUM_FACES_2D])
                                                    

  integer(ip), target :: P4EST_CORNER_IN_FACE_2D(NUM_FACES_2D,NUM_CORNERS_2D) = & 
                                                  reshape([ 1,-1, 1,-1,&
                                                           -1, 1, 2,-1,&
                                                            2,-1,-1, 1,&
                                                           -1, 2,-1, 2],[NUM_FACES_2D,NUM_CORNERS_2D])
  
  integer(ip), target :: P4EST_OPPOSITE_CORNER_2D(NUM_CORNERS_2D) = [ 4, 3, 2, 1 ]
  integer(ip), target :: P4EST_2_FEMPAR_CORNER_2D(NUM_CORNERS_2D) = [ 1, 2, 3, 4 ]
  integer(ip), target :: P4EST_2_FEMPAR_FACE_2D  (NUM_FACES_2D)   = [ 3, 4, 1, 2 ]
  integer(ip), target :: FEMPAR_2_P4EST_FACE_2D  (NUM_FACES_2D)   = [ 3, 4, 1, 2 ]
  
  integer(ip), target :: P4EST_FACES_SUBFACE_IMPROPER_VERTEX_LID_2D(NUM_FACES_2D,NUM_SUBFACES_FACE_2D) = & 
                                                  reshape([ 4, 3, 4, 2, &
                                                            2, 1, 3, 1 ] ,[NUM_FACES_2D,NUM_SUBFACES_FACE_2D])
                                                  
  integer(ip), target :: P4EST_FACES_IN_TOUCH_2D(NUM_FACES_2D) = [ 2, 1, 4, 3 ]
  
  integer(ip), target :: P4EST_FACES_SUBFACE_FACE_NEIGHBOUR_2D(NUM_FACES_2D/2, NUM_SUBFACES_FACE_2D,1) = &
                                                  reshape([ 4, 2, &
                                                            3, 1 ] ,[NUM_FACES_2D/2,NUM_SUBFACES_FACE_2D,1])
                                                  
  integer(ip), target :: P4EST_FACES_SUBFACE_SUBFACE_NEIGHBOUR_2D(NUM_SUBFACES_FACE_2D,1) = &
                                                  reshape( [ 2, 1 ], [NUM_SUBFACES_FACE_2D,1] )
  
  ! For 3D
  integer(ip), parameter :: NUM_SUBCELLS_IN_TOUCH_FACE_3D = 4
  integer(ip), parameter :: NUM_SUBCELLS_IN_TOUCH_EDGE_3D = 2
  integer(ip), parameter :: NUM_CORNERS_3D                = 8
  integer(ip), parameter :: NUM_FACES_3D                  = 6
  integer(ip), parameter :: NUM_EDGES_3D                  = 12
  integer(ip), parameter :: NUM_SUBFACES_FACE_3D          = 4
  integer(ip), parameter :: NUM_SUBEDGES_EDGE_3D          = 2
  integer(ip), parameter :: NUM_FACE_CORNERS_3D           = 4
  integer(ip), parameter :: NUM_FACE_EDGES_3D             = 4
  integer(ip), parameter :: NUM_EDGE_CORNERS_3D           = 2
  integer(ip), parameter :: NUM_FACES_AT_CORNER_3D        = 3
  integer(ip), parameter :: NUM_FACES_AT_EDGE_3D          = 2
  integer(ip), parameter :: NUM_EDGES_AT_CORNER_3D        = 3
  integer(ip), parameter :: NUM_VEFS_3D                   = NUM_CORNERS_3D+NUM_FACES_3D+NUM_EDGES_3D
  integer(ip), parameter :: MAX_NUM_CELLS_AROUND_3D       = 8 


  integer(ip), target :: P4EST_OPPOSITE_FACE_3D(NUM_FACES_3D) = [ 2, 1, 4, 3, 6, 5 ]

  
  integer(ip), target :: P4EST_FACE_CORNERS_3D(NUM_FACE_CORNERS_3D,NUM_FACES_3D) = & 
                                                  reshape([1, 3, 5, 7,&
                                                           2, 4, 6, 8,&  
                                                           1, 2, 5, 6,&
                                                           3, 4, 7, 8,&
                                                           1, 2, 3, 4,&
                                                           5, 6, 7, 8], [NUM_FACE_CORNERS_3D,NUM_FACES_3D])

  integer(ip), target :: P4EST_EDGE_CORNERS_3D(NUM_EDGE_CORNERS_3D,NUM_EDGES_3D) = & 
                                                  reshape([ 1, 2,&
                                                            3, 4,&
                                                            5, 6,&
                                                            7, 8,&
                                                            1, 3,&
                                                            2, 4,&
                                                            5, 7,&
                                                            6, 8,&
                                                            1, 5,&
                                                            2, 6,&
                                                            3, 7,&
                                                            4, 8], [NUM_EDGE_CORNERS_3D,NUM_EDGES_3D])

  integer(ip), target :: P4EST_FACE_EDGES_3D(NUM_FACE_EDGES_3D,NUM_FACES_3D) = & 
                                                  reshape([ 9, 11,  5,  7,&
                                                           10, 12,  6,  8,&
                                                            9, 10,  1,  3,&
                                                           11, 12,  2,  4,&
                                                            5,  6,  1,  2,&
                                                            7,  8,  3,  4 ], [NUM_FACE_EDGES_3D,NUM_FACES_3D])

  integer(ip), target :: P4EST_FACES_AT_CORNER_3D(NUM_FACES_AT_CORNER_3D, NUM_CORNERS_3D) = & 
                                                    reshape([ 1, 3, 5,&
                                                              2, 3, 5,&
                                                              1, 4, 5,&
                                                              2, 4, 5,&
                                                              1, 3, 6,&
                                                              2, 3, 6,&
                                                              1, 4, 6,&
                                                              2, 4, 6], [NUM_FACES_AT_CORNER_3D, NUM_CORNERS_3D])

  integer(ip), target :: P4EST_FACES_AT_EDGE_3D(NUM_FACES_AT_EDGE_3D, NUM_EDGES_3D) = & 
                                                    reshape([ 3, 5,&
                                                              4, 5,&
                                                              3, 6,&
                                                              4, 6,&
                                                              1, 5,&
                                                              2, 5,&
                                                              1, 6,&
                                                              2, 6,&
                                                              1, 3,&
                                                              2, 3,&
                                                              1, 4,&
                                                              2, 4 ], [NUM_FACES_AT_EDGE_3D, NUM_EDGES_3D])

  integer(ip), target :: P4EST_EDGES_AT_CORNER_3D(NUM_EDGES_AT_CORNER_3D, NUM_CORNERS_3D) = & 
                                                    reshape([ 1,  5,  9,&
                                                              1,  6, 10,&
                                                              2,  5, 11,&
                                                              2,  6, 12,&
                                                              3,  7,  9,&
                                                              3,  8, 10,&
                                                              4,  7, 11,&
                                                              4,  8, 12], [NUM_EDGES_AT_CORNER_3D, NUM_CORNERS_3D])

  integer(ip), target :: P4EST_CORNER_IN_FACE_3D(NUM_FACES_3D,NUM_CORNERS_3D) = &
                                                  reshape([ 1,-1, 1,-1, 1,-1,&
                                                           -1, 1, 2,-1, 2,-1,&
                                                            2,-1,-1, 1, 3,-1,&
                                                           -1, 2,-1, 2, 4,-1,&
                                                            3,-1, 3,-1,-1, 1,&
                                                           -1, 3, 4,-1,-1, 2,&
                                                            4,-1,-1, 3,-1, 3,&
                                                           -1, 4,-1, 4,-1, 4 ],[NUM_FACES_3D,NUM_CORNERS_3D])

  integer(ip), target :: P4EST_EDGE_IN_FACE_3D(NUM_FACES_3D,NUM_EDGES_3D) = &
                                                  reshape([ -1, -1,  3, -1,  3, -1,&
                                                            -1, -1, -1,  3,  4, -1,&
                                                            -1, -1,  4, -1, -1,  3,&
                                                            -1, -1, -1,  4, -1,  4,&
                                                             3, -1, -1, -1,  1, -1,&
                                                            -1,  3, -1, -1,  2, -1,&
                                                             4, -1, -1, -1, -1,  1,&
                                                            -1,  4, -1, -1, -1,  2,&
                                                             1, -1,  1, -1, -1, -1,&
                                                            -1,  1,  2, -1, -1, -1,&
                                                             2, -1, -1,  1, -1, -1,&
                                                            -1,  2, -1,  2, -1, -1 ],[NUM_FACES_3D,NUM_EDGES_3D])

  integer(ip), target :: P4EST_CORNER_IN_EDGE_3D(NUM_EDGES_3D,NUM_CORNERS_3D) = &
                                                  reshape([ 1,-1,-1,-1, 1,-1,-1,-1, 1,-1,-1,-1,&
                                                            2,-1,-1,-1,-1, 1,-1,-1,-1, 1,-1,-1,&
                                                           -1, 1,-1,-1, 2,-1,-1,-1,-1,-1, 1,-1,&
                                                           -1, 2,-1,-1,-1, 2,-1,-1,-1,-1,-1, 1,&
                                                           -1,-1, 1,-1,-1,-1, 1,-1, 2,-1,-1,-1,&
                                                           -1,-1, 2,-1,-1,-1,-1, 1,-1, 2,-1,-1,&
                                                           -1,-1,-1, 1,-1,-1, 2,-1,-1,-1, 2,-1,&
                                                           -1,-1,-1, 2,-1,-1,-1, 2,-1,-1,-1, 2 ],[NUM_EDGES_3D,NUM_CORNERS_3D])

  integer(ip), target :: FEMPAR_SUBCELLS_IN_TOUCH_FACE_3D(NUM_SUBCELLS_IN_TOUCH_FACE_3D,NUM_FACES_3D) = &
                                                    reshape([ 1,  2,  3,  4,&
                                                              5,  6,  7,  8,&
                                                              1,  2,  5,  6,&
                                                              3,  4,  7,  8,&
                                                              1,  3,  5,  7,&
                                                              2,  4,  6,  8 ], [NUM_SUBCELLS_IN_TOUCH_FACE_3D, NUM_FACES_3D])
                                                    
  integer(ip), target :: P4EST_SUBCELLS_IN_TOUCH_FACE_3D(NUM_SUBCELLS_IN_TOUCH_FACE_3D,NUM_FACES_3D) = &
                                                    reshape([ 1,  3,  5,  7,&
                                                              2,  4,  6,  8,&
                                                              1,  2,  5,  6,&
                                                              3,  4,  7,  8,&
                                                              1,  2,  3,  4,&
                                                              5,  6,  7,  8 ], [NUM_SUBCELLS_IN_TOUCH_FACE_3D, NUM_FACES_3D])                                                    
                                                    

  integer(ip), target :: FEMPAR_EDGE_OF_SUBCELLS_IN_TOUCH_FACE_3D(NUM_SUBCELLS_IN_TOUCH_FACE_3D,NUM_FACES_3D) = &
                                                    reshape([ 6,  2,  1,  5,&
                                                              8,  4,  3,  7,&
                                                             10,  3,  1,  9,&
                                                             12,  4,  2, 11,&
                                                             11,  7,  5,  9,&
                                                             12,  8,  6, 10 ], [NUM_SUBCELLS_IN_TOUCH_FACE_3D, NUM_FACES_3D])

  integer(ip), target :: FEMPAR_SUBCELLS_IN_TOUCH_EDGE_3D(NUM_SUBCELLS_IN_TOUCH_EDGE_3D,NUM_EDGES_3D) = &
                                                    reshape([ 1,  2,&
                                                              3,  4,&
                                                              5,  6,&
                                                              7,  8,&
                                                              1,  3,&
                                                              2,  4,&
                                                              5,  7,&
                                                              6,  8,&
                                                              1,  5,&
                                                              2,  6,&
                                                              3,  7,&
                                                              4,  8 ], [NUM_SUBCELLS_IN_TOUCH_EDGE_3D, NUM_EDGES_3D])
                                                    
  integer(ip), target :: P4EST_SUBCELLS_IN_TOUCH_EDGE_3D(NUM_SUBCELLS_IN_TOUCH_EDGE_3D,NUM_EDGES_3D) =  &
                                                    reshape([ 1,  2,&
                                                              3,  4,&
                                                              5,  6,&
                                                              7,  8,&
                                                              1,  3,&
                                                              2,  4,&
                                                              5,  7,&
                                                              6,  8,&
                                                              1,  5,&
                                                              2,  6,&
                                                              3,  7,&
                                                              4,  8 ], [NUM_SUBCELLS_IN_TOUCH_EDGE_3D, NUM_EDGES_3D])                                          
                                                    

  integer(ip), target :: P4EST_OPPOSITE_CORNER_3D(NUM_CORNERS_3D) = [ 8, 7, 6, 5, 4, 3, 2, 1 ]
  integer(ip), target :: P4EST_2_FEMPAR_CORNER_3D(NUM_CORNERS_3D) = [ 1, 2, 3, 4, 5, 6, 7, 8 ]
  integer(ip), target :: P4EST_2_FEMPAR_FACE_3D  (NUM_FACES_3D)   = [ 5, 6, 3, 4, 1, 2 ]
  integer(ip), target :: FEMPAR_2_P4EST_FACE_3D  (NUM_FACES_3D)   = [ 5, 6, 3, 4, 1, 2 ]
  integer(ip), target :: P4EST_2_FEMPAR_EDGE_3D  (NUM_EDGES_3D)   = [ 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12 ]
  
  integer(ip), target :: P4EST_FACES_SUBFACE_FACE_NEIGHBOUR_3D(NUM_FACES_3D/2, NUM_SUBFACES_FACE_3D, 2) = &
                                                  reshape( [ 4, 2, 2, &
                                                             3, 1, 1, &
                                                             4, 2, 2, & 
                                                             3, 1, 1, & 
                                                             6, 6, 4, &
                                                             6, 6, 4, &
                                                             5, 5, 3, &
                                                             5, 5, 3 ], [NUM_FACES_3D/2,NUM_SUBFACES_FACE_3D,2])
                                                  
  integer(ip), target :: P4EST_FACES_SUBFACE_EDGE_NEIGHBOUR_3D(NUM_FACES_3D/2, NUM_SUBFACES_FACE_3D) = &
                                                  reshape( [ 4, 8, 12, &
                                                             3, 7, 11, &
                                                             2, 6, 10, & 
                                                             1, 5,  9 ], [NUM_FACES_3D/2,NUM_SUBFACES_FACE_3D] )                                                
                                                  
  integer(ip), target :: P4EST_FACES_SUBFACE_SUBFACE_NEIGHBOUR_3D(NUM_SUBFACES_FACE_3D,2) = &
                                                  reshape( [ 2, 1, 4, 3, &
                                                             3, 4, 1, 2 ], [NUM_SUBFACES_FACE_3D,2])
                                                  
  integer(ip), target :: P4EST_EDGES_SUBEDGE_FACE_NEIGHBOUR_3D(NUM_EDGES_3D/4, NUM_SUBEDGES_EDGE_3D) = &
                                                  reshape( [ 2, 4, 6, &
                                                             1, 3, 5 ], [NUM_EDGES_3D/4,NUM_SUBEDGES_EDGE_3D] )
  

end module p4est_triangulation_parameters_names
