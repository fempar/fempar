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
module umfpack_interface
  ! use types
  use iso_c_binding
  implicit none

#ifdef ENABLE_UMFPACK

! I did the Fortran2003 bindings for this particular version of UMFPACK.
! In the future, it might be required that these definitions are adapted
! to be conformal with (possible) changes in umfpack.h header file.  

! #define UMFPACK_DATE "Oct 10, 2014"
! #define UMFPACK_VER_CODE(main,sub) ((main) * 1000 + (sub))
! #define UMFPACK_MAIN_VERSION 5
! #define UMFPACK_SUB_VERSION 7
! #define UMFPACK_SUBSUB_VERSION 1
! #define UMFPACK_VER
! #define UMFPACK_VER_CODE(UMFPACK_MAIN_VERSION,UMFPACK_SUB_VERSION)

! /*
!-------------------------------------------------------------------------- */
! /* size of Info and Control arrays */
!-------------------------------------------------------------------------- */
! /* These might be larger in future versions, since there are only 3 unused
!  * entries in Info, and no unused entries in Control. */
!  
! #define UMFPACK_INFO 90
! #define UMFPACK_CONTROL 20

 integer(c_int), parameter :: UMFPACK_INFO = 90
 integer(c_int), parameter :: UMFPACK_CONTROL = 20


! /*
!-------------------------------------------------------------------------- */
! /* contents of Info */
!-------------------------------------------------------------------------- */
 
!/* Note that umfpack_report.m must coincide with these definitions.  S is
! * the submatrix of A after removing row/col singletons and empty rows/cols.
!*/
 
!/* returned by all routines that use Info: */
!#define UMFPACK_STATUS 0        /* UMFPACK_OK, or other result */
!#define UMFPACK_NROW 1          /* n_row input value */
!#define UMFPACK_NCOL 16         /* n_col input value */
!#define UMFPACK_NZ 2            /* # of entries in A */
! 
  integer(c_int), parameter :: UMFPACK_STATUS =  0 
  integer(c_int), parameter :: UMFPACK_NROW   =  1 
  integer(c_int), parameter :: UMFPACK_NCOL   =  16 
  integer(c_int), parameter :: UMFPACK_NZ     =  2 

! /* computed in UMFPACK_*symbolic and UMFPACK_numeric: */
! #define UMFPACK_SIZE_OF_UNIT 3          /* sizeof (Unit) */
  integer(c_int), parameter :: UMFPACK_SIZE_OF_UNIT = 3 

 
! /* computed in UMFPACK_*symbolic: */
!#define UMFPACK_SIZE_OF_INT 4           /* sizeof (int) */
!#define UMFPACK_SIZE_OF_LONG 5          /* sizeof (SuiteSparse_long) */
!#define UMFPACK_SIZE_OF_POINTER 6       /* sizeof (void *) */
!#define UMFPACK_SIZE_OF_ENTRY 7         /* sizeof (Entry), real or complex */
!#define UMFPACK_NDENSE_ROW 8            /* number of dense rows */
!#define UMFPACK_NEMPTY_ROW 9            /* number of empty rows */
!#define UMFPACK_NDENSE_COL 10           /* number of dense rows */
!#define UMFPACK_NEMPTY_COL 11           /* number of empty rows */
!#define UMFPACK_SYMBOLIC_DEFRAG 12      /* # of memory compactions */
!#define UMFPACK_SYMBOLIC_PEAK_MEMORY 13 /* memory used by symbolic analysis*/
!#define UMFPACK_SYMBOLIC_SIZE 14        /* size of Symbolic object, in Units*/
!#define UMFPACK_SYMBOLIC_TIME 15        /* time (sec.) for symbolic analysis*/
!#define UMFPACK_SYMBOLIC_WALLTIME 17    /* wall clock time for sym. analysis*/
!#define UMFPACK_STRATEGY_USED 18        /* strategy used: sym, unsym */
!#define UMFPACK_ORDERING_USED 19        /* ordering used: colamd, amd, given*/
!#define UMFPACK_QFIXED 31               /* whether Q is fixed or refined */
!#define UMFPACK_DIAG_PREFERRED 32       /* whether diagonal pivoting attempted*/
!#define UMFPACK_PATTERN_SYMMETRY 33     /* symmetry of pattern of S */
!#define UMFPACK_NZ_A_PLUS_AT 34         /* nnz (S+S'), excl. diagonal */
!#define UMFPACK_NZDIAG 35               /* nnz (diag (S)) */

  integer(c_int), parameter :: UMFPACK_SIZE_OF_INT     = 4 
  integer(c_int), parameter :: UMFPACK_SIZE_OF_LONG    = 5 
  integer(c_int), parameter :: UMFPACK_SIZE_OF_POINTER = 6 
  integer(c_int), parameter :: UMFPACK_SIZE_OF_ENTRY   = 7 

  integer(c_int), parameter :: UMFPACK_NDENSE_ROW = 8           
  integer(c_int), parameter :: UMFPACK_NEMPTY_ROW = 9 
  integer(c_int), parameter :: UMFPACK_NDENSE_COL = 10
  integer(c_int), parameter :: UMFPACK_NEMPTY_COL = 11 
  integer(c_int), parameter :: UMFPACK_SYMBOLIC_DEFRAG = 12
  integer(c_int), parameter :: UMFPACK_SYMBOLIC_PEAK_MEMORY = 13 

  integer(c_int), parameter :: UMFPACK_SYMBOLIC_SIZE = 14
  integer(c_int), parameter :: UMFPACK_SYMBOLIC_TIME = 15
  integer(c_int), parameter :: UMFPACK_SYMBOLIC_WALLTIME = 17
  integer(c_int), parameter :: UMFPACK_STRATEGY_USED = 18
  integer(c_int), parameter :: UMFPACK_ORDERING_USED = 19
  integer(c_int), parameter :: UMFPACK_QFIXED = 31       
  integer(c_int), parameter :: UMFPACK_DIAG_PREFERRED = 32
  integer(c_int), parameter :: UMFPACK_PATTERN_SYMMETRY = 33
  integer(c_int), parameter :: UMFPACK_NZ_A_PLUS_AT = 34
  integer(c_int), parameter :: UMFPACK_NZDIAG = 35

! /* AMD statistics, computed in UMFPACK_*symbolic: */
! #define UMFPACK_SYMMETRIC_LUNZ 36	/* nz in L+U, if AMD ordering used */
! #define UMFPACK_SYMMETRIC_FLOPS 37	/* flops for LU, if AMD ordering used */
! #define UMFPACK_SYMMETRIC_NDENSE 38	/* # of "dense" rows/cols in S+S' */
! #define UMFPACK_SYMMETRIC_DMAX 39	/* max nz in cols of L, for AMD */

  integer(c_int), parameter :: UMFPACK_SYMMETRIC_LUNZ = 36
  integer(c_int), parameter :: UMFPACK_SYMMETRIC_FLOPS = 37 
  integer(c_int), parameter :: UMFPACK_SYMMETRIC_NDENSE = 38
  integer(c_int), parameter :: UMFPACK_SYMMETRIC_DMAX = 39

! /* 51:55 unused */

! /* statistcs for singleton pruning */
! #define UMFPACK_COL_SINGLETONS 56	/* # of column singletons */
! #define UMFPACK_ROW_SINGLETONS 57	/* # of row singletons */
! #define UMFPACK_N2 58			/* size of S */
! #define UMFPACK_S_SYMMETRIC 59	/* 1 if S square and symmetricly perm.*/

  integer(c_int), parameter :: UMFPACK_COL_SINGLETONS = 56
  integer(c_int), parameter :: UMFPACK_ROW_SINGLETONS = 57
  integer(c_int), parameter :: UMFPACK_N2 = 58
  integer(c_int), parameter :: UMFPACK_S_SYMMETRIC = 59

! /* estimates computed in UMFPACK_*symbolic: */
! #define UMFPACK_NUMERIC_SIZE_ESTIMATE 20    /* final size of Numeric->Memory */
! #define UMFPACK_PEAK_MEMORY_ESTIMATE 21	    /* for symbolic & numeric */
! #define UMFPACK_FLOPS_ESTIMATE 22	    /* flop count */
! #define UMFPACK_LNZ_ESTIMATE 23		    /* nz in L, incl. diagonal */
! #define UMFPACK_UNZ_ESTIMATE 24		    /* nz in U, incl. diagonal */
! #define UMFPACK_VARIABLE_INIT_ESTIMATE 25   /* initial size of Numeric->Memory*/
! #define UMFPACK_VARIABLE_PEAK_ESTIMATE 26   /* peak size of Numeric->Memory */
! #define UMFPACK_VARIABLE_FINAL_ESTIMATE 27  /* final size of Numeric->Memory */
! #define UMFPACK_MAX_FRONT_SIZE_ESTIMATE 28  /* max frontal matrix size */
! #define UMFPACK_MAX_FRONT_NROWS_ESTIMATE 29 /* max # rows in any front */
! #define UMFPACK_MAX_FRONT_NCOLS_ESTIMATE 30 /* max # columns in any front */

  integer(c_int), parameter :: UMFPACK_NUMERIC_SIZE_ESTIMATE = 20
  integer(c_int), parameter :: UMFPACK_PEAK_MEMORY_ESTIMATE = 21
  integer(c_int), parameter :: UMFPACK_FLOPS_ESTIMATE = 22
  integer(c_int), parameter :: UMFPACK_LNZ_ESTIMATE = 23
  integer(c_int), parameter :: UMFPACK_UNZ_ESTIMATE = 24
  integer(c_int), parameter :: UMFPACK_VARIABLE_INIT_ESTIMATE = 25
  integer(c_int), parameter :: UMFPACK_VARIABLE_PEAK_ESTIMATE = 26
  integer(c_int), parameter :: UMFPACK_VARIABLE_FINAL_ESTIMATE = 27
  integer(c_int), parameter :: UMFPACK_MAX_FRONT_SIZE_ESTIMATE = 28
  integer(c_int), parameter :: UMFPACK_MAX_FRONT_NROWS_ESTIMATE = 29
  integer(c_int), parameter :: UMFPACK_MAX_FRONT_NCOLS_ESTIMATE = 30

! /* exact values, (estimates shown above) computed in UMFPACK_numeric: */
! #define UMFPACK_NUMERIC_SIZE 40	    /* final size of Numeric->Memory */
! #define UMFPACK_PEAK_MEMORY 41	    /* for symbolic & numeric */
! #define UMFPACK_FLOPS 42		    /* flop count */
! #define UMFPACK_LNZ 43		    /* nz in L, incl. diagonal */
! #define UMFPACK_UNZ 44		    /* nz in U, incl. diagonal */
! #define UMFPACK_VARIABLE_INIT 45	    /* initial size of Numeric->Memory*/
! #define UMFPACK_VARIABLE_PEAK 46	    /* peak size of Numeric->Memory */
! #define UMFPACK_VARIABLE_FINAL 47	    /* final size of Numeric->Memory */
! #define UMFPACK_MAX_FRONT_SIZE 48	    /* max frontal matrix size */
! #define UMFPACK_MAX_FRONT_NROWS 49	    /* max # rows in any front */
! #define UMFPACK_MAX_FRONT_NCOLS 50	    /* max # columns in any front */

 integer(c_int), parameter :: UMFPACK_NUMERIC_SIZE = 40
 integer(c_int), parameter :: UMFPACK_PEAK_MEMORY = 41
 integer(c_int), parameter :: UMFPACK_FLOPS = 42
 integer(c_int), parameter :: UMFPACK_LNZ = 43
 integer(c_int), parameter :: UMFPACK_UNZ = 44
 integer(c_int), parameter :: UMFPACK_VARIABLE_INIT = 45
 integer(c_int), parameter :: UMFPACK_VARIABLE_PEAK = 46
 integer(c_int), parameter :: UMFPACK_VARIABLE_FINAL = 47
 integer(c_int), parameter :: UMFPACK_MAX_FRONT_SIZE = 48
 integer(c_int), parameter :: UMFPACK_MAX_FRONT_NROWS = 49
 integer(c_int), parameter :: UMFPACK_MAX_FRONT_NCOLS = 50

! /* computed in UMFPACK_numeric: */
! #define UMFPACK_NUMERIC_DEFRAG 60	    /* # of garbage collections */
! #define UMFPACK_NUMERIC_REALLOC 61	    /* # of memory reallocations */
! #define UMFPACK_NUMERIC_COSTLY_REALLOC 62 /* # of costlly memory realloc's */
! #define UMFPACK_COMPRESSED_PATTERN 63	    /* # of integers in LU pattern */
! #define UMFPACK_LU_ENTRIES 64		    /* # of reals in LU factors */
! #define UMFPACK_NUMERIC_TIME 65	    /* numeric factorization time */
! #define UMFPACK_UDIAG_NZ 66		    /* nz on diagonal of U */
! #define UMFPACK_RCOND 67		    /* est. reciprocal condition # */
! #define UMFPACK_WAS_SCALED 68		    /* none, max row, or sum row */
! #define UMFPACK_RSMIN 69		    /* min (max row) or min (sum row) */
! #define UMFPACK_RSMAX 70		    /* max (max row) or max (sum row) */
! #define UMFPACK_UMIN 71		    /* min abs diagonal entry of U */
! #define UMFPACK_UMAX 72		    /* max abs diagonal entry of U */
! #define UMFPACK_ALLOC_INIT_USED 73	    /* alloc_init parameter used */
! #define UMFPACK_FORCED_UPDATES 74	    /* # of forced updates */
! #define UMFPACK_NUMERIC_WALLTIME 75	    /* numeric wall clock time */
! #define UMFPACK_NOFF_DIAG 76		    /* number of off-diagonal pivots */

! #define UMFPACK_ALL_LNZ 77		    /* nz in L, if no dropped entries */
! #define UMFPACK_ALL_UNZ 78		    /* nz in U, if no dropped entries */
! #define UMFPACK_NZDROPPED 79		    /* # of dropped small entries */

  integer(c_int), parameter :: UMFPACK_NUMERIC_DEFRAG = 60
  integer(c_int), parameter :: UMFPACK_NUMERIC_REALLOC = 61
  integer(c_int), parameter :: UMFPACK_NUMERIC_COSTLY_REALLOC = 62
  integer(c_int), parameter :: UMFPACK_COMPRESSED_PATTERN = 63
  integer(c_int), parameter :: UMFPACK_LU_ENTRIES = 64
  integer(c_int), parameter :: UMFPACK_NUMERIC_TIME = 65
  integer(c_int), parameter :: UMFPACK_UDIAG_NZ = 66
  integer(c_int), parameter :: UMFPACK_RCOND = 67
  integer(c_int), parameter :: UMFPACK_WAS_SCALED = 68
  integer(c_int), parameter :: UMFPACK_RSMIN = 69
  integer(c_int), parameter :: UMFPACK_RSMAX = 70
  integer(c_int), parameter :: UMFPACK_UMIN = 71
  integer(c_int), parameter :: UMFPACK_UMAX = 72
  integer(c_int), parameter :: UMFPACK_ALLOC_INIT_USED = 73
  integer(c_int), parameter :: UMFPACK_FORCED_UPDATES = 74
  integer(c_int), parameter :: UMFPACK_NUMERIC_WALLTIME = 75
  integer(c_int), parameter :: UMFPACK_NOFF_DIAG = 76

  integer(c_int), parameter :: UMFPACK_ALL_LNZ = 77
  integer(c_int), parameter :: UMFPACK_ALL_UNZ = 78
  integer(c_int), parameter :: UMFPACK_NZDROPPED = 79

! /* computed in UMFPACK_solve: */
! #define UMFPACK_IR_TAKEN 80	    /* # of iterative refinement steps taken */
! #define UMFPACK_IR_ATTEMPTED 81   /* # of iter. refinement steps attempted */
! #define UMFPACK_OMEGA1 82	    /* omega1, sparse backward error estimate */
! #define UMFPACK_OMEGA2 83	    /* omega2, sparse backward error estimate */
! #define UMFPACK_SOLVE_FLOPS 84    /* flop count for solve */
! #define UMFPACK_SOLVE_TIME 85	    /* solve time (seconds) */
! #define UMFPACK_SOLVE_WALLTIME 86 /* solve time (wall clock, seconds) */

  integer(c_int), parameter :: UMFPACK_IR_TAKEN = 80
  integer(c_int), parameter :: UMFPACK_IR_ATTEMPTED = 81
  integer(c_int), parameter :: UMFPACK_OMEGA1 = 82
  integer(c_int), parameter :: UMFPACK_OMEGA2 = 83
  integer(c_int), parameter :: UMFPACK_SOLVE_FLOPS = 84
  integer(c_int), parameter :: UMFPACK_SOLVE_TIME = 85
  integer(c_int), parameter :: UMFPACK_SOLVE_WALLTIME = 86

! /* Info [87, 88, 89] unused */

! /* Unused parts of Info may be used in future versions of UMFPACK. */

! /* -------------------------------------------------------------------------- */
! /* contents of Control */
! /* -------------------------------------------------------------------------- */

! /* used in all UMFPACK_report_* routines: */
! #define UMFPACK_PRL 0	/* print level */
integer(c_int), parameter :: UMFPACK_PRL = 0

! /* used in UMFPACK_*symbolic only: */
! #define UMFPACK_DENSE_ROW 1      /* dense row parameter */
! #define UMFPACK_DENSE_COL 2      /* dense col parameter */
! #define UMFPACK_BLOCK_SIZE 4     /* BLAS-3 block size */
! #define UMFPACK_STRATEGY 5       /* auto, symmetric, or unsym. */
! #define UMFPACK_ORDERING 10      /* ordering method to use */
! #define UMFPACK_FIXQ 13          /* -1: no fixQ, 0: default, 1: fixQ */
! #define UMFPACK_AMD_DENSE 14     /* for AMD ordering */
! #define UMFPACK_AGGRESSIVE 19	   /* whether or not to use aggressive */
! #define UMFPACK_SINGLETONS 11    /* singleton filter on if true */

 integer(c_int), parameter :: UMFPACK_DENSE_ROW = 1
 integer(c_int), parameter :: UMFPACK_DENSE_COL = 2
 integer(c_int), parameter :: UMFPACK_BLOCK_SIZE = 4
 integer(c_int), parameter :: UMFPACK_STRATEGY = 5
 integer(c_int), parameter :: UMFPACK_ORDERING = 10
 integer(c_int), parameter :: UMFPACK_FIXQ = 13
 integer(c_int), parameter :: UMFPACK_AMD_DENSE = 14
 integer(c_int), parameter :: UMFPACK_AGGRESSIVE = 19
 integer(c_int), parameter :: UMFPACK_SINGLETONS = 11

! /* used in UMFPACK_numeric only: */
! #define UMFPACK_PIVOT_TOLERANCE 3	/* threshold partial pivoting setting */
! #define UMFPACK_ALLOC_INIT 6		/* initial allocation ratio */
! #define UMFPACK_SYM_PIVOT_TOLERANCE 15	/* threshold, only for diag. entries */
! #define UMFPACK_SCALE 16		/* what row scaling to do */
! #define UMFPACK_FRONT_ALLOC_INIT 17	/* frontal matrix allocation ratio */
! #define UMFPACK_DROPTOL 18		/* drop tolerance for entries in L,U */

  integer(c_int), parameter :: UMFPACK_PIVOT_TOLERANCE = 3
  integer(c_int), parameter :: UMFPACK_ALLOC_INIT = 6
  integer(c_int), parameter :: UMFPACK_SYM_PIVOT_TOLERANCE = 15
  integer(c_int), parameter :: UMFPACK_SCALE = 16
  integer(c_int), parameter :: UMFPACK_FRONT_ALLOC_INIT = 17
  integer(c_int), parameter :: UMFPACK_DROPTOL = 18

! /* used in UMFPACK_*solve only: */
! #define UMFPACK_IRSTEP 7		/* max # of iterative refinements */
  integer(c_int), parameter :: UMFPACK_IRSTEP = 7

! /* compile-time settings - Control [8..11] cannot be changed at run time: */
! #define UMFPACK_COMPILED_WITH_BLAS 8	    /* uses the BLAS */
  integer(c_int), parameter :: UMFPACK_COMPILED_WITH_BLAS = 8

! /* 9,12: unused */
!/* -------------------------------------------------------------------------- */

! /* Control [UMFPACK_STRATEGY] is one of the following: */
! #define UMFPACK_STRATEGY_AUTO 0		/* use sym. or unsym. strategy */
! #define UMFPACK_STRATEGY_UNSYMMETRIC 1	/* COLAMD(A), coletree postorder,
!                                            not prefer diag*/
! #define UMFPACK_STRATEGY_OBSOLETE 2     /* 2-by-2 is no longer available */
! #define UMFPACK_STRATEGY_SYMMETRIC 3	/* AMD(A+A'), no coletree postorder,
!                                            prefer diagonal */
  integer(c_int), parameter :: UMFPACK_STRATEGY_AUTO = 0
  integer(c_int), parameter :: UMFPACK_STRATEGY_UNSYMMETRIC = 1
  integer(c_int), parameter :: UMFPACK_STRATEGY_OBSOLETE = 2
  integer(c_int), parameter :: UMFPACK_STRATEGY_SYMMETRIC = 3

! /* Control [UMFPACK_SCALE] is one of the following: */
!#define UMFPACK_SCALE_NONE 0	/* no scaling */
!#define UMFPACK_SCALE_SUM 1	/* default: divide each row by sum (abs (row))*/
!#define UMFPACK_SCALE_MAX 2	/* divide each row by max (abs (row)) */

  integer(c_int), parameter :: UMFPACK_SCALE_NONE = 0
  integer(c_int), parameter :: UMFPACK_SCALE_SUM = 1
  integer(c_int), parameter :: UMFPACK_SCALE_MAX = 2

! /* Control [UMFPACK_ORDERING] and Info [UMFPACK_ORDERING_USED] are one of: */
! #define UMFPACK_ORDERING_CHOLMOD 0      /* use CHOLMOD (AMD/COLAMD then METIS)*/
! #define UMFPACK_ORDERING_AMD 1          /* use AMD/COLAMD */
! #define UMFPACK_ORDERING_GIVEN 2        /* user-provided Qinit */
! #define UMFPACK_ORDERING_METIS 3        /* use METIS */
! #define UMFPACK_ORDERING_BEST 4         /* try many orderings, pick best */
! #define UMFPACK_ORDERING_NONE 5         /* natural ordering */
! #define UMFPACK_ORDERING_USER 6         /* user-provided function */

  integer(c_int), parameter :: UMFPACK_ORDERING_CHOLMOD = 0 
  integer(c_int), parameter :: UMFPACK_ORDERING_AMD = 1
  integer(c_int), parameter :: UMFPACK_ORDERING_GIVEN = 2
  integer(c_int), parameter :: UMFPACK_ORDERING_METIS = 3
  integer(c_int), parameter :: UMFPACK_ORDERING_BEST = 4
  integer(c_int), parameter :: UMFPACK_ORDERING_NONE = 5
  integer(c_int), parameter :: UMFPACK_ORDERING_USER = 6
! /* AMD/COLAMD means: use AMD for symmetric strategy, COLAMD for unsymmetric */

! /* -------------------------------------------------------------------------- */
! /* default values of Control: */
! /* -------------------------------------------------------------------------- */

! #define UMFPACK_DEFAULT_PRL 1
! #define UMFPACK_DEFAULT_DENSE_ROW 0.2
! #define UMFPACK_DEFAULT_DENSE_COL 0.2
! #define UMFPACK_DEFAULT_PIVOT_TOLERANCE 0.1
! #define UMFPACK_DEFAULT_SYM_PIVOT_TOLERANCE 0.001
! #define UMFPACK_DEFAULT_BLOCK_SIZE 32
! #define UMFPACK_DEFAULT_ALLOC_INIT 0.7
! #define UMFPACK_DEFAULT_FRONT_ALLOC_INIT 0.5
! #define UMFPACK_DEFAULT_IRSTEP 2
! #define UMFPACK_DEFAULT_SCALE UMFPACK_SCALE_SUM
! #define UMFPACK_DEFAULT_STRATEGY UMFPACK_STRATEGY_AUTO
! #define UMFPACK_DEFAULT_AMD_DENSE AMD_DEFAULT_DENSE
! #define UMFPACK_DEFAULT_FIXQ 0
! #define UMFPACK_DEFAULT_AGGRESSIVE 1
! #define UMFPACK_DEFAULT_DROPTOL 0
! #define UMFPACK_DEFAULT_ORDERING UMFPACK_ORDERING_AMD
! #define UMFPACK_DEFAULT_SINGLETONS TRUE

!  integer(c_int), parameter :: UMFPACK_DEFAULT_PRL =  1
!  integer(c_int), parameter :: UMFPACK_DEFAULT_DENSE_ROW =  0.2
!  integer(c_int), parameter :: UMFPACK_DEFAULT_DENSE_COL =  0.2
!  integer(c_int), parameter :: UMFPACK_DEFAULT_PIVOT_TOLERANCE =  0.1
!  integer(c_int), parameter :: UMFPACK_DEFAULT_SYM_PIVOT_TOLERANCE =  0.001
!  integer(c_int), parameter :: UMFPACK_DEFAULT_BLOCK_SIZE =  32
!  integer(c_int), parameter :: UMFPACK_DEFAULT_ALLOC_INIT =  0.7
!  integer(c_int), parameter :: UMFPACK_DEFAULT_FRONT_ALLOC_INIT =  0.5
!  integer(c_int), parameter :: UMFPACK_DEFAULT_IRSTEP =  2
!  integer(c_int), parameter :: UMFPACK_DEFAULT_SCALE =  UMFPACK_SCALE_SUM
!  integer(c_int), parameter :: UMFPACK_DEFAULT_STRATEGY = UMFPACK_STRATEGY_AUTO
!  integer(c_int), parameter :: UMFPACK_DEFAULT_AMD_DENSE = AMD_DEFAULT_DENSE
!  integer(c_int), parameter :: UMFPACK_DEFAULT_FIXQ =  0
!  integer(c_int), parameter :: UMFPACK_DEFAULT_AGGRESSIVE =  1
!  integer(c_int), parameter :: UMFPACK_DEFAULT_DROPTOL =  0
!  integer(c_int), parameter :: UMFPACK_DEFAULT_ORDERING = UMFPACK_ORDERING_AMD
!  integer(c_int), parameter :: UMFPACK_DEFAULT_SINGLETONS = TRUE

! /* default values of Control may change in future versions of UMFPACK. */
! /* -------------------------------------------------------------------------- */
! /* status codes */
! /* -------------------------------------------------------------------------- */

!#define UMFPACK_OK (0)
integer(c_int), parameter :: UMFPACK_OK = 0

! /* status > 0 means a warning, but the method was successful anyway. */
! /* A Symbolic or Numeric object was still created. */
! #define UMFPACK_WARNING_singular_matrix (1)
integer(c_int), parameter :: UMFPACK_WARNING_singular_matrix = 1

! /* The following warnings were added in umfpack_*_get_determinant */
! #define UMFPACK_WARNING_determinant_underflow (2)
! #define UMFPACK_WARNING_determinant_overflow (3)
integer(c_int), parameter :: UMFPACK_WARNING_determinant_underflow = 2
integer(c_int), parameter :: UMFPACK_WARNING_determinant_overflow = 3

! /* status < 0 means an error, and the method was not successful. */
! /* No Symbolic of Numeric object was created. */
! #define UMFPACK_ERROR_out_of_memory (-1)
! #define UMFPACK_ERROR_invalid_Numeric_object (-3)
! #define UMFPACK_ERROR_invalid_Symbolic_object (-4)
! #define UMFPACK_ERROR_argument_missing (-5)
! #define UMFPACK_ERROR_n_nonpositive (-6)
! #define UMFPACK_ERROR_invalid_matrix (-8)
! #define UMFPACK_ERROR_different_pattern (-11)
! #define UMFPACK_ERROR_invalid_system (-13)
! #define UMFPACK_ERROR_invalid_permutation (-15)
! #define UMFPACK_ERROR_internal_error (-911) /* yes, call me if you get this! */
! #define UMFPACK_ERROR_file_IO (-17)
! #define UMFPACK_ERROR_ordering_failed (-18)

  integer(c_int), parameter :: UMFPACK_ERROR_out_of_memory = -1
  integer(c_int), parameter :: UMFPACK_ERROR_invalid_Numeric_object = -3
  integer(c_int), parameter :: UMFPACK_ERROR_invalid_Symbolic_object = -4
  integer(c_int), parameter :: UMFPACK_ERROR_argument_missing = -5
  integer(c_int), parameter :: UMFPACK_ERROR_n_nonpositive = -6
  integer(c_int), parameter :: UMFPACK_ERROR_invalid_matrix = -8
  integer(c_int), parameter :: UMFPACK_ERROR_different_pattern = -11
  integer(c_int), parameter :: UMFPACK_ERROR_invalid_system = -13
  integer(c_int), parameter :: UMFPACK_ERROR_invalid_permutation = -15
  integer(c_int), parameter :: UMFPACK_ERROR_internal_error = -911
  integer(c_int), parameter :: UMFPACK_ERROR_file_IO = -17
  integer(c_int), parameter :: UMFPACK_ERROR_ordering_failed = -18

!/* -------------------------------------------------------------------------- */
!/* solve codes */
!/* -------------------------------------------------------------------------- */
!
!/* Solve the system ( )x=b, where ( ) is defined below.  "t" refers to the */
!/* linear algebraic transpose (complex conjugate if A is complex), or the (') */
!/* operator in MATLAB.  "at" refers to the array transpose, or the (.') */
!/* operator in MATLAB. */

!#define UMFPACK_A	(0)	/* Ax=b    */
!#define UMFPACK_At	(1)	/* A'x=b   */
!#define UMFPACK_Aat	(2)	/* A.'x=b  */

!#define UMFPACK_Pt_L	(3)	/* P'Lx=b  */
!#define UMFPACK_L	(4)	/* Lx=b    */
!#define UMFPACK_Lt_P	(5)	/* L'Px=b  */
!#define UMFPACK_Lat_P	(6)	/* L.'Px=b */
!#define UMFPACK_Lt	(7)	/* L'x=b   */
!#define UMFPACK_Lat	(8)	/* L.'x=b  */

!#define UMFPACK_U_Qt	(9)	/* UQ'x=b  */
!#define UMFPACK_U	(10)	/* Ux=b    */
!#define UMFPACK_Q_Ut	(11)	/* QU'x=b  */
!#define UMFPACK_Q_Uat	(12)	/* QU.'x=b */
!#define UMFPACK_Ut	(13)	/* U'x=b   */
!#define UMFPACK_Uat	(14)	/* U.'x=b  */

  integer(c_int), parameter :: UMFPACK_A =  0
  integer(c_int), parameter :: UMFPACK_At =  1
  integer(c_int), parameter :: UMFPACK_Aat =  2

  integer(c_int), parameter :: UMFPACK_Pt_L =  3
  integer(c_int), parameter :: UMFPACK_L =  4
  integer(c_int), parameter :: UMFPACK_Lt_P =  5
  integer(c_int), parameter :: UMFPACK_Lat_P =  6
  integer(c_int), parameter :: UMFPACK_Lt =  7
  integer(c_int), parameter :: UMFPACK_Lat =  8

  integer(c_int), parameter :: UMFPACK_U_Qt =  9
  integer(c_int), parameter :: UMFPACK_U =  10
  integer(c_int), parameter :: UMFPACK_Q_Ut =  11
  integer(c_int), parameter :: UMFPACK_Q_Uat =  12
  integer(c_int), parameter :: UMFPACK_Ut =  13
  integer(c_int), parameter :: UMFPACK_Uat =  14

  interface
      ! int umfpack_di_symbolic
      ! (
      !   int n_row,
      !   int n_col,
      !   const int Ap [ ],
      !   const int Ai [ ],
      !   const double Ax [ ],
      !   void **Symbolic,
      !   const double Control [UMFPACK_CONTROL],
      !   double Info [UMFPACK_INFO]
      ! ) ;
      function umfpack_di_symbolic(n_row, n_col, Ap, Ai, Ax, Symbolic, Control, Info) bind(c,NAME='umfpack_di_symbolic')
         use iso_c_binding
         integer(c_int), value       :: n_row
         integer(c_int), value       :: n_col
         integer(c_int), intent(in)  :: Ap(*)
         integer(c_int), intent(in)  :: Ai(*)
         ! I had to declare Ax as type(c_ptr) as it is required in some contexts
         ! to pass C_NULL_PTR as actual argument to this dummy argument
         type(c_ptr)   , value, intent(in)  :: Ax 
         type(c_ptr)   , intent(out) :: Symbolic 
         real(c_double), intent(in)  :: Control(*) 
         real(c_double), intent(out) :: Info(*) 
         integer(c_int)              :: umfpack_di_symbolic 
      end function umfpack_di_symbolic

      ! int umfpack_di_numeric
      ! (
      !  const int Ap [ ],
      !  const int Ai [ ],
      !  const double Ax [ ],
      !  void *Symbolic,
      !  void **Numeric,
      !  const double Control [UMFPACK_CONTROL],
      !  double Info [UMFPACK_INFO]
      ! ) ;
       function umfpack_di_numeric(Ap, Ai, Ax, Symbolic, Numeric, Control, Info) & 
               & bind(c,NAME='umfpack_di_numeric')
         use iso_c_binding
         integer(c_int), intent(in)     :: Ap(*)
         integer(c_int), intent(in)     :: Ai(*)
         real(c_double), intent(in)     :: Ax(*)
         type(c_ptr), value, intent(in) :: Symbolic 
         type(c_ptr), intent(out)       :: Numeric 
         real(c_double), intent(in)     :: Control(*) 
         real(c_double), intent(out)    :: Info(*) 
         integer(c_int)                 :: umfpack_di_numeric 
       end function umfpack_di_numeric

      ! int umfpack_di_solve
      ! (
      !  int sys,
      !  const int Ap [ ],
      !  const int Ai [ ],
      !  const double Ax [ ],
      !  double X [ ],
      !  const double B [ ],
      !  void *Numeric,
      !  const double Control [UMFPACK_CONTROL],
      !  double Info [UMFPACK_INFO]
      ! ) ;
      function umfpack_di_solve(sys, Ap, Ai, Ax, X, B, Numeric, Control, Info) & 
               & bind(c,NAME='umfpack_di_solve')
         use iso_c_binding
         integer(c_int), value          :: sys 
         integer(c_int), intent(in)     :: Ap(*)
         integer(c_int), intent(in)     :: Ai(*)
         real(c_double), intent(in)     :: Ax(*)
         real(c_double), intent(out)    :: X(*) 
         real(c_double), intent(in)     :: B(*) 
         type(c_ptr), value, intent(in) :: Numeric 
         real(c_double), intent(in)     :: Control(*) 
         real(c_double), intent(out)    :: Info(*) 
         integer(c_int)                 :: umfpack_di_solve 
      end function umfpack_di_solve
      
      ! void umfpack_di_free_symbolic
      ! (
      ! void **Symbolic
      ! ) ;
      subroutine umfpack_di_free_symbolic(Symbolic) &
               & bind(c,NAME='umfpack_di_free_symbolic')
         use iso_c_binding
         type(c_ptr) :: Symbolic 
      end subroutine umfpack_di_free_symbolic

      ! void umfpack_di_free_numeric
      ! (
      ! void **Numeric
      ! ) ;
      subroutine umfpack_di_free_numeric(Numeric) &
                & bind(c,NAME='umfpack_di_free_numeric')
          use iso_c_binding
         type(c_ptr) :: Numeric 
      end subroutine umfpack_di_free_numeric

      ! void umfpack_di_defaults
      ! (
      !  double Control [UMFPACK_CONTROL]
      ! ) ;
      subroutine umfpack_di_defaults(Control) &
                & bind(c,NAME='umfpack_di_defaults')
         use iso_c_binding
         real(c_double), intent(out)  :: Control(*) 
      end subroutine umfpack_di_defaults
  end interface

#endif

end module umfpack_interface
