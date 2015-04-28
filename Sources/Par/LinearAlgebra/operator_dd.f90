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

! Operator required for non-overlapping domain decomposition methods 
module fem_operator_dd_names

  ! Abstract modules
  use solver_base
  use abstract_solver

  ! Serial modules
  use types
  use memor
  use fem_graph_names
  use fem_matrix_names
  use matvec_dof
  use fem_vector_names
  use fem_precond_names
  use serial_environment_names
  
  ! Parallel modules
  use dof_distribution_names
  use fem_graph_distribution
  use fem_matrix_distribution

# include "debug.i90"
  
  implicit none
  private

  ! Interface to Intel MKL mkl_dcsrmm 
  ! (Level 3 Sparse BLAS dense-sparse matrix multiplication)
#ifdef ENABLE_MKL
  interface
     ! http://software.intel.com/sites/products/documentation/hpc/
     ! compilerpro/en-us/cpp/win/mkl/refman/bla/functn_mkl_dcsrmm.html#functn_mkl_dcsrmm
     SUBROUTINE mkl_dcsrmm (transa, m, n, k, alpha, matdescra, val, indx, pntrb, pntre, b, ldb, beta, c, ldc)
       implicit none
       CHARACTER*1 transa
       CHARACTER matdescra(*)
       INTEGER m, n, k, ldb, ldc
       INTEGER indx(*), pntrb(m), pntre(m)
       REAL alpha, beta
       REAL val(*), b(ldb,*), c(ldc,*)
     end SUBROUTINE mkl_dcsrmm

     ! http://software.intel.com/sites/products/documentation/hpc/compilerpro/en-us/
     ! cpp/win/mkl/refman/bla/functn_mkl_dcsrmv.html#functn_mkl_dcsrmv
     SUBROUTINE mkl_dcsrmv(transa, m, k, alpha, matdescra, val, indx, pntrb, pntre, x, beta, y)
       implicit none
       CHARACTER*1 transa
       CHARACTER matdescra(*)
       INTEGER m, k
       INTEGER indx(*), pntrb(m), pntre(m)
       REAL alpha, beta
       REAL val(*), x(*), y(*)
     end SUBROUTINE mkl_dcsrmv

  end interface
#endif 

  type fem_operator_dd
     type ( fem_graph )  :: A_II_gr 
     type ( fem_graph )  :: A_IG_gr
     type ( fem_graph )  :: A_GI_gr
     type ( fem_graph )  :: A_GG_gr

     type ( fem_matrix )  :: A_II 
     type ( fem_matrix )  :: A_IG 
     type ( fem_matrix )  :: A_GI
     type ( fem_matrix )  :: A_GG

     integer (ip)             :: symm
     type ( fem_precond )     :: M_II

     type(solver_control)    , pointer :: spars
     logical                           :: spars_allocated
     type(fem_precond_params), pointer :: ppars
     logical                           :: ppars_allocated

     type(dof_distribution) , pointer  :: dof_dist => NULL()
  end type fem_operator_dd

  ! Types
  public :: fem_operator_dd


  public :: fem_operator_dd_create, fem_operator_dd_free,             &  
            fem_operator_dd_ass_struct, fem_operator_dd_fill_val,     & 
            fem_operator_dd_apply, fem_operator_dd_solve_A_II,        &
            fem_operator_dd_apply_A_IG, fem_operator_dd_apply_A_GI,   &
            fem_operator_dd_apply_A_GG, fem_operator_dd_evaluate_rhs, &
            fem_operator_dd_solve_interior, fem_operator_dd_info,     &
            fem_operator_dd_compute_ut_op_u, fem_operator_dd_apply_A_IG_several_rhs, &
            fem_operator_dd_apply_A_GI_plus_A_GG

contains

  subroutine fem_operator_dd_create ( f_matrix, dof_dist, f_operator, spars, ppars, symm, sign )
    implicit none
    ! Parameters
    type(fem_matrix)        , intent(in), target           :: f_matrix
    type(dof_distribution)  , intent(in), target           :: dof_dist
    type(fem_operator_dd)   , intent(out)                  :: f_operator
    type(solver_control)    , intent(in), target, optional :: spars
    type(fem_precond_params), intent(in), target, optional :: ppars
    ! With this two optional parameters one may select the
    ! structure and sign of the Schur complement operator
    ! independently of f_matrix
    integer (ip)            , intent(in), optional :: symm
    integer (ip)            , intent(in), optional :: sign

    ! Locals
    type    (fem_matrix) :: mat_dum
    integer (ip)         :: symm_, sign_

    if ( present(symm) ) then
       symm_ = symm
    else
       symm_ = f_matrix%symm
    end if

    if ( present(sign) ) then
       sign_ = sign
    else
       sign_ = f_matrix%sign
    end if

    f_operator%dof_dist => dof_dist
    f_operator%symm = symm_

    if ( present(ppars) ) then
       f_operator%ppars => ppars
       f_operator%ppars_allocated = .false.
    else
       allocate(f_operator%ppars)
       f_operator%ppars_allocated = .true.
    end if

    if ( present(spars)  ) then
       f_operator%spars => spars
       f_operator%spars_allocated = .false.
    else
       allocate(f_operator%spars)
       f_operator%spars_allocated = .true.
    end if

    ! AFM: Here we are initializing a precond but using a matrix that
    !      is not the matrix it will be used for. But as I KNOW that the
    !      only info needed at this stage is symmetry and sign 
    !      I will apply the following (DIRTY) patch.
    mat_dum%type    = f_matrix%type
    mat_dum%symm    = symm_
    mat_dum%sign    = sign_
    call fem_precond_create( mat_dum , f_operator%M_II, f_operator%ppars)

  end subroutine fem_operator_dd_create 

  subroutine fem_operator_dd_free ( f_operator, mode )
    implicit none

    ! Parameters
    type(fem_operator_dd), intent(inout)             :: f_operator
    integer(ip)       , intent(in)                   :: mode

    ! Locals
    type (fem_vector) :: dum 

    if ( mode == free_clean ) then
       call fem_precond_free ( precond_free_clean  , f_operator%M_II )
    else if ( mode == free_only_struct  ) then
       call fem_precond_free ( precond_free_struct , f_operator%M_II )
    else if ( mode == free_only_values ) then
       call fem_precond_free ( precond_free_values , f_operator%M_II )
    end if

    ! Free memory associated to the blocks of the operator
    if ( mode == free_only_values ) then
       call fem_matrix_free ( f_operator%A_II, free_only_values )
       call fem_matrix_free ( f_operator%A_IG, free_only_values )
       call fem_matrix_free ( f_operator%A_GG, free_only_values )
    end if

    if ( mode == free_only_struct ) then
       call fem_matrix_free ( f_operator%A_II, free_only_struct )
       call fem_matrix_free ( f_operator%A_IG, free_only_struct )
       call fem_matrix_free ( f_operator%A_GG, free_only_struct )

       call fem_graph_free ( f_operator%A_II_gr )
       call fem_graph_free ( f_operator%A_IG_gr )
       call fem_graph_free ( f_operator%A_GG_gr )
    end if

    if ( f_operator%symm == symm_false ) then
       if ( mode == free_only_values ) then
          call fem_matrix_free( f_operator%A_GI, free_only_values  )
       end if

       if ( mode == free_only_struct ) then
          call fem_matrix_free( f_operator%A_GI, free_only_struct  )
          call fem_graph_free ( f_operator%A_GI_gr ) 
       end if
    end if

    
    if ( mode == free_clean ) then
       if ( f_operator%ppars_allocated  ) then
          deallocate(f_operator%ppars)
       else
          nullify(f_operator%ppars)
       end if

       if ( f_operator%spars_allocated  ) then
          deallocate(f_operator%spars)
       else
          nullify(f_operator%spars)
       end if
    end if

  end subroutine fem_operator_dd_free


  !=============================================================================
  subroutine fem_operator_dd_ass_struct ( f_matrix, f_operator )
    implicit none

    ! Parameters 
    type(fem_matrix)  , intent(in)                    :: f_matrix
    type(fem_operator_dd), intent(inout)              :: f_operator

      ! Split graph of local process into 2x2 block partitioning
    if ( f_operator%symm == symm_false ) then
        call fem_graph_split_2x2_partitioning ( csr, &
                                                f_matrix%gr, & 
                                                f_operator%dof_dist, & 
                                                G_II=f_operator%A_II_gr, G_IG=f_operator%A_IG_gr, &
                                                G_GI=f_operator%A_GI_gr, G_GG=f_operator%A_GG_gr  )
        call fem_matrix_graph(f_operator%A_II_gr, f_operator%A_II)
        call fem_matrix_graph(f_operator%A_IG_gr, f_operator%A_IG)
        call fem_matrix_graph(f_operator%A_GI_gr, f_operator%A_GI)
        call fem_matrix_graph(f_operator%A_GG_gr, f_operator%A_GG)
        ! call fem_graph_print ( 6, f_operator%A_GI%gr )
     else if ( f_operator%symm == symm_true ) then
        call fem_graph_split_2x2_partitioning_symm ( csr_symm, & 
                                                     f_matrix%gr, & 
                                                     f_operator%dof_dist, & 
                                                     G_II=f_operator%A_II_gr, & 
                                                     G_IG=f_operator%A_IG_gr, &
                                                     G_GG=f_operator%A_GG_gr )
        call fem_matrix_graph(f_operator%A_II_gr, f_operator%A_II)
        call fem_matrix_graph(f_operator%A_IG_gr, f_operator%A_IG)
        call fem_matrix_graph(f_operator%A_GG_gr, f_operator%A_GG)
     end if

     call fem_precond_symbolic(f_operator%A_II, f_operator%M_II)

     ! call fem_graph_print ( 6, f_operator%A_II%gr ) ! DBG:
     ! call fem_graph_print ( 6, f_operator%A_IG%gr ) ! DBG:
     ! call fem_graph_print ( 6, f_operator%A_GI%gr ) ! DBG:
     ! call fem_graph_print ( 6, f_operator%A_GG%gr ) ! DBG:

  end subroutine fem_operator_dd_ass_struct

  !=============================================================================
  subroutine fem_operator_dd_fill_val ( f_matrix, f_operator ) !, me )
    use stdio
    implicit none
    
    ! Parameters
    type(fem_matrix)   , intent(in)              :: f_matrix
    type(fem_operator_dd), intent(inout), target :: f_operator
    ! integer(ip)          , intent(in)            :: me

    ! Locals
    ! integer(ip) :: lunou

    ! Split matrix of local process into 2x2 block partitioning
    if ( f_operator%symm == symm_false ) then
        call fem_matrix_split_2x2_partitioning ( f_operator%symm,                   &
                                                 f_matrix,                          &
                                                 f_operator%dof_dist,               & 
                                                 f_operator%A_II, f_operator%A_IG,  &
                                                 f_operator%A_GI, f_operator%A_GG )

     else if ( f_operator%symm == symm_true ) then
        call fem_matrix_split_2x2_partitioning ( f_operator%symm,                            &
                                                 f_matrix,                                   & 
                                                 f_operator%dof_dist,                        & 
                                                 A_II=f_operator%A_II, A_IG=f_operator%A_IG, &  
                                                 A_GG=f_operator%A_GG )
     end if

     ! lunou = io_open ( trim('matrix' // trim(ch(me)) // trim('.') // 'mtx' ), 'write')
     ! call fem_matrix_print_matrix_market ( lunou, f_operator%A_II )
     ! call io_close (lunou)

     call fem_precond_numeric(f_operator%A_II, f_operator%M_II)

     ! call fem_matrix_print ( 6, f_matrix )  ! DBG:
     ! call fem_matrix_print ( 6, f_operator%A_II )    ! DBG:
     ! call fem_matrix_print ( 6, f_operator%A_IG )    ! DBG:
     ! call fem_matrix_print ( 6, f_operator%A_GI )    ! DBG: 
     ! call fem_matrix_print ( 6, f_operator%A_GG )    ! DBG:
  end subroutine fem_operator_dd_fill_val

  !=============================================================================
  ! Computes y_G -< S_G x_G, where S_G the local Schur
  ! complement is given by S_G <- A_GG - A_G_I A_II^-1 A_IG
  ! VERY IMPORTANT: fem_operator_dd_apply is well-defined if and only if x_G and 
  !                 y_G are vectors on the interface
  !=============================================================================  
  subroutine fem_operator_dd_apply ( f_operator, x_G, y_G )
    implicit none

    ! Parameters 
    type(fem_operator_dd), intent(inout),target             :: f_operator
    type(fem_vector)     , intent(in)                       :: x_G
    type(fem_vector)     , intent(inout)                    :: y_G

    ! Locals 
    type(fem_vector)       :: ws_vec       ! Local workspace vector

    type(fem_vector)       :: ws_vec_I     ! View of the interior nodes
    type(fem_vector)       :: ws_vec_I_2   ! Extra ws required by pardiso
    type(fem_vector)       :: ws_vec_G     ! View of the interfaces
                                           
    integer(ip)            :: ni, nl

 
    ni = f_operator%dof_dist%ni   ! # of interior nodes
    nl = f_operator%dof_dist%nl   ! # of local nodes

    ! Allocate space for ws_vec
    call fem_vector_alloc ( nl,   ws_vec)
 
    ! Create a view of ws_vec to point to the interfaces
    call fem_vector_create_view ( ws_vec,    1, ni, ws_vec_I )
    call fem_vector_clone ( ws_vec_I, ws_vec_I_2 )

    ! Create a view of ws_vec to point to the interior nodes
    call fem_vector_create_view ( ws_vec, ni+1, nl, ws_vec_G )

    ! Perform the following computations:
    ! (a) y_G <- A_GG * x_G
    call fem_operator_dd_apply_A_GG ( f_operator, x_G, y_G )
    ! call fem_matrix_print (6, f_operator%A_GG) ! DBG:
    ! call fem_vector_print (6,x_G)     ! DBG:
    ! call fem_vector_print (6,y_G)     ! DBG:

    ! (b) ws_vec_I_2  <- A_IG * x_G
    call fem_operator_dd_apply_A_IG ( f_operator, x_G, ws_vec_I_2 )
    ! call fem_vector_print (6, x_G)         ! DBG:
    ! call fem_vector_print (6, ws_vec_I_2)  ! DBG:

    ! (c) ws_vec_I  <- A_II^-1 *  ws_vec_I_2
    call fem_operator_dd_solve_A_II ( f_operator, ws_vec_I_2, ws_vec_I )
    ! call fem_vector_print (6, ws_vec_I_2)     ! DBG:
    ! call fem_vector_print (6, ws_vec_I)       ! DBG: 

    ! (d) ws_vec_G <- A_GI * ws_vec_I
    call fem_operator_dd_apply_A_GI ( f_operator, ws_vec_I, ws_vec_G )
    ! call fem_matrix_print (6, f_operator%A_IG) ! DBG:
    ! call fem_vector_print (6, ws_vec_I)     ! DBG: 
    ! call fem_vector_print (6, ws_vec_G)     ! DBG:

    ! (e) y_G <- y_G - ws_vec_G 
    call fem_vector_mxpy( ws_vec_G, y_G )

    ! Free memory of local worksize vector
    call fem_vector_free ( ws_vec )
    call fem_vector_free ( ws_vec_I_2 )

  end subroutine fem_operator_dd_apply 

  
  ! Computes y_I <- A_II^-1 * x_I
  subroutine fem_operator_dd_solve_A_II ( f_operator, x_I, y_I )
    implicit none

    ! Parameters 
    type(fem_operator_dd), intent(inout),target :: f_operator
    type(fem_vector)     , intent(in)   ,target :: x_I
    type(fem_vector)     , intent(inout),target :: y_I

    ! Locals
    type(serial_environment) :: senv

    call abstract_solve(f_operator%A_II, f_operator%M_II, x_I, y_I, f_operator%spars, senv) 

  end subroutine fem_operator_dd_solve_A_II

  ! Computes y_I <- A_IG * x_G
  subroutine fem_operator_dd_apply_A_IG ( f_operator, x_G, y_I )
    implicit none

    ! Parameters 
    type(fem_operator_dd), intent(inout)                    :: f_operator
    type(fem_vector   )  , intent(in)                       :: x_G
    type(fem_vector   )  , intent(inout)                    :: y_I
    integer(ip)           :: ni

    ni = f_operator%dof_dist%ni ! # of interior nodes

#ifdef ENABLE_MKL
    call mkl_dcsrmv ( 'N', & 
                      f_operator%A_IG%gr%nv, &
                      f_operator%A_IG%gr%nv2, &
                      1.0, &
                      'GXXF', & 
                      f_operator%A_IG%a, &
                      f_operator%A_IG%gr%ja, & 
                      f_operator%A_IG%gr%ia(1), & 
                      f_operator%A_IG%gr%ia(2), & 
                      x_G%b, &
                      0.0,   &
                      y_I%b)
#else
    call matvec_csr_scal ( f_operator%A_IG%nd1,   &   
                         & f_operator%A_IG%nd2,   &
                         & f_operator%A_IG%gr%nv, &
                         & f_operator%A_IG%gr%nv2,&
                         & f_operator%A_IG%gr%ia, &
                         & f_operator%A_IG%gr%ja, &
                         & f_operator%A_IG%a,     &
                         & x_G%b,        & 
                         & y_I%b  )
#endif

  end subroutine fem_operator_dd_apply_A_IG

  ! Computes y_G < A_GI * x_I
  subroutine fem_operator_dd_apply_A_GI ( f_operator, x_I, y_G )
    implicit none

    ! Parameters 
    type(fem_operator_dd), intent(inout)  :: f_operator
    type(fem_vector)     , intent(in)     :: x_I
    type(fem_vector)     , intent(inout)  :: y_G

    if ( f_operator%A_GG%symm == symm_false ) then
#ifdef ENABLE_MKL
       call mkl_dcsrmv ( 'N', & 
                         f_operator%A_GI%gr%nv, &
                         f_operator%A_GI%gr%nv2, &
                         1.0, &
                         'GXXF', & 
                         f_operator%A_GI%a, &
                         f_operator%A_GI%gr%ja, & 
                         f_operator%A_GI%gr%ia(1), & 
                         f_operator%A_GI%gr%ia(2), & 
                         x_I%b, &
                         0.0,    &
                         y_G%b) 
#else
       call matvec_csr_scal ( f_operator%A_GI%nd1,   &    
                          &   f_operator%A_GI%nd2,   &
                          &   f_operator%A_GI%gr%nv, &
                          &   f_operator%A_GI%gr%nv2,&
                          &   f_operator%A_GI%gr%ia, &
                          &   f_operator%A_GI%gr%ja, &
                          &   f_operator%A_GI%a,     &
                          &   x_I%b,        & ! x_I
                          &   y_G%b         ) ! y_G
#endif
    else if ( f_operator%A_GG%symm == symm_true ) then
#ifdef ENABLE_MKL
       call mkl_dcsrmv ( 'T', & 
                         f_operator%A_IG%gr%nv, &
                         f_operator%A_IG%gr%nv2, &
                         1.0, &
                         'GXXF', & 
                         f_operator%A_IG%a, &
                         f_operator%A_IG%gr%ja, & 
                         f_operator%A_IG%gr%ia(1), & 
                         f_operator%A_IG%gr%ia(2), & 
                         x_I%b, &
                         0.0,   &
                         y_G%b)
#else
       call matvec_csr_trans_scal ( f_operator%A_IG%nd1,   &   
                    &               f_operator%A_IG%nd2,   &
                    &               f_operator%A_IG%gr%nv, &
                    &               f_operator%A_IG%gr%nv2,&
                    &               f_operator%A_IG%gr%ia, &
                    &               f_operator%A_IG%gr%ja, &
                    &               f_operator%A_IG%a,     &
                    &               x_I%b,        & ! x_I 
                    &               y_G%b         ) ! y_G 
#endif
    end if
  end subroutine fem_operator_dd_apply_A_GI


  ! Computes y_G < (A_GI*x_I + A_GG*x_G) 
  subroutine fem_operator_dd_apply_A_GI_plus_A_GG ( f_operator, x, y_G )
    implicit none

    ! Parameters 
    type(fem_operator_dd), intent(inout)  :: f_operator
    type(fem_vector)     , intent(in)     :: x
    type(fem_vector)     , intent(inout)  :: y_G
    type(fem_vector)                      :: x_I, x_G

   call fem_vector_create_view ( x, 1, f_operator%A_II%gr%nv, x_I )
   call fem_vector_create_view ( x, f_operator%A_II%gr%nv+1, f_operator%A_II%gr%nv+f_operator%A_GG%gr%nv, x_G )

   if ( f_operator%A_GG%symm == symm_false ) then
#ifdef ENABLE_MKL
      call mkl_dcsrmv ( 'N', & 
                        f_operator%A_GI%gr%nv, &
                        f_operator%A_GI%gr%nv2, &
                        1.0, &
                        'GXXF', & 
                        f_operator%A_GI%a, &
                        f_operator%A_GI%gr%ja, & 
                        f_operator%A_GI%gr%ia(1), & 
                        f_operator%A_GI%gr%ia(2), & 
                        x_I%b, &
                        0.0, &
                        y_G%b)
#else
      call matvec_csr_scal ( f_operator%A_GI%nd1,   &    
                         &   f_operator%A_GI%nd2,   &
                         &   f_operator%A_GI%gr%nv, &
                         &   f_operator%A_GI%gr%nv2,&
                         &   f_operator%A_GI%gr%ia, &
                         &   f_operator%A_GI%gr%ja, &
                         &   f_operator%A_GI%a,     &
                         &   x_I%b,        & ! x_I
                         &   y_G%b         ) ! y_G
#endif
   else if ( f_operator%A_GG%symm == symm_true ) then
#ifdef ENABLE_MKL
      call mkl_dcsrmv ( 'T', & 
                        f_operator%A_IG%gr%nv, &
                        f_operator%A_IG%gr%nv2, &
                        1.0, &
                        'GXXF', & 
                        f_operator%A_IG%a, &
                        f_operator%A_IG%gr%ja, & 
                        f_operator%A_IG%gr%ia(1), & 
                        f_operator%A_IG%gr%ia(2), & 
                        x_I%b, &
                        0.0,   &
                        y_G%b)
#else
      call matvec_csr_trans_scal ( f_operator%A_IG%nd1,   &   
                   &               f_operator%A_IG%nd2,   &
                   &               f_operator%A_IG%gr%nv, &
                   &               f_operator%A_IG%gr%nv2,&
                   &               f_operator%A_IG%gr%ia, &
                   &               f_operator%A_IG%gr%ja, &
                   &               f_operator%A_IG%a,     &
                   &               x_I%b,        & ! x_I 
                   &               y_G%b         ) ! y_G 
#endif
   end if

   if ( f_operator%A_GG%symm == symm_false ) then
#ifdef ENABLE_MKL
      call mkl_dcsrmv ( 'N', & 
                        f_operator%A_GG%gr%nv, &
                        f_operator%A_GG%gr%nv2, &
                        1.0, &
                        'GXXF', & 
                        f_operator%A_GG%a, &
                        f_operator%A_GG%gr%ja, & 
                        f_operator%A_GG%gr%ia(1), & 
                        f_operator%A_GG%gr%ia(2), & 
                        x_G%b, &
                        1.0,   &
                        y_G%b)
#else
        check(1==0)
!!$           call matvec_csr_scal ( f_operator%A_GG%nd1,    &    
!!$                          &       f_operator%A_GG%nd2,    &
!!$                          &       f_operator%A_GG%gr%nv,  &
!!$                          &       f_operator%A_GG%gr%nv2, &
!!$                          &       f_operator%A_GG%gr%ia,  &
!!$                          &       f_operator%A_GG%gr%ja,  &
!!$                          &       f_operator%A_GG%a,      &
!!$                          &       x_G%b,         & 
!!$                          &       y_G%b          )
#endif
           ! call fem_vector_print ( 6, y_G ) ! DBG:
     else if ( f_operator%A_GG%symm == symm_true ) then
#ifdef ENABLE_MKL
        call mkl_dcsrmv ( 'N', & 
                          f_operator%A_GG%gr%nv, &
                          f_operator%A_GG%gr%nv2, &
                          1.0, &
                          'SUNF', & 
                          f_operator%A_GG%a, &
                          f_operator%A_GG%gr%ja, & 
                          f_operator%A_GG%gr%ia(1), & 
                          f_operator%A_GG%gr%ia(2), & 
                          x_G%b, &
                          1.0,   &
                          y_G%b)
#else
        check(1==0)
!!$           call matvec_csr_symm_scal ( f_operator%A_GG%nd1,    &   
!!$                        &              f_operator%A_GG%nd2,    &
!!$                        &              f_operator%A_GG%gr%nv,  &
!!$                        &              f_operator%A_GG%gr%nv2, &
!!$                        &              f_operator%A_GG%gr%ia,  &
!!$                        &              f_operator%A_GG%gr%ja,  &
!!$                        &              f_operator%A_GG%a,      &
!!$                        &              x_G%b,         & 
!!$                        &              y_G%b          )
#endif
        ! call fem_vector_print ( 6, y_G ) ! DBG:
     end if

  end subroutine fem_operator_dd_apply_A_GI_plus_A_GG




  ! Computes y_G < A_GG * x_G
  subroutine fem_operator_dd_apply_A_GG ( f_operator, x_G, y_G )
    implicit none
    ! Parameters 
    type(fem_operator_dd), intent(inout)                 :: f_operator
    type(fem_vector)     , intent(in)                    :: x_G
    type(fem_vector)     , intent(inout)                 :: y_G

    if ( f_operator%A_GG%symm == symm_false ) then
#ifdef ENABLE_MKL
!!$            call mkl_dcsrgemv( 'N',                   & 
!!$                               f_operator%A_GG%gr%nv, &
!!$                               f_operator%A_GG%a,     &
!!$                               f_operator%A_GG%gr%ia, & 
!!$                               f_operator%A_GG%gr%ja, &
!!$                               x_G%b,                 &
!!$                               y_G%b)
       call mkl_dcsrmv ( 'N', & 
                         f_operator%A_GG%gr%nv, &
                         f_operator%A_GG%gr%nv2, &
                         1.0, &
                         'GXXF', & 
                         f_operator%A_GG%a, &
                         f_operator%A_GG%gr%ja, & 
                         f_operator%A_GG%gr%ia(1), & 
                         f_operator%A_GG%gr%ia(2), & 
                         x_G%b, &
                         0.0,   &
                         y_G%b)
#else
       call matvec_csr_scal ( f_operator%A_GG%nd1,    &    
                      &       f_operator%A_GG%nd2,    &
                      &       f_operator%A_GG%gr%nv,  &
                      &       f_operator%A_GG%gr%nv2, &
                      &       f_operator%A_GG%gr%ia,  &
                      &       f_operator%A_GG%gr%ja,  &
                      &       f_operator%A_GG%a,      &
                      &       x_G%b,         & 
                      &       y_G%b          )
#endif
       ! call fem_vector_print ( 6, y_G ) ! DBG:
    else if ( f_operator%A_GG%symm == symm_true ) then
#ifdef ENABLE_MKL
!!$            call mkl_dcsrsymv( 'U',                   & 
!!$                               f_operator%A_GG%gr%nv, &
!!$                               f_operator%A_GG%a,     &
!!$                               f_operator%A_GG%gr%ia, & 
!!$                               f_operator%A_GG%gr%ja, &
!!$                               x_G%b,                 &
!!$                               y_G%b)
       call mkl_dcsrmv ( 'N', & 
                         f_operator%A_GG%gr%nv, &
                         f_operator%A_GG%gr%nv2, &
                         1.0, &
                         'SUNF', & 
                         f_operator%A_GG%a, &
                         f_operator%A_GG%gr%ja, & 
                         f_operator%A_GG%gr%ia(1), & 
                         f_operator%A_GG%gr%ia(2), & 
                         x_G%b, &
                         0.0,   &
                         y_G%b)
#else
       call matvec_csr_symm_scal ( f_operator%A_GG%nd1,    &   
                    &              f_operator%A_GG%nd2,    &
                    &              f_operator%A_GG%gr%nv,  &
                    &              f_operator%A_GG%gr%nv2, &
                    &              f_operator%A_GG%gr%ia,  &
                    &              f_operator%A_GG%gr%ja,  &
                    &              f_operator%A_GG%a,      &
                    &              x_G%b,         & 
                    &              y_G%b          )
#endif
       ! call fem_vector_print ( 6, y_G ) ! DBG:
    end if
  end subroutine  fem_operator_dd_apply_A_GG

  ! Computes g_G < b_G - A_GI * A_II^-1 b_I
  subroutine fem_operator_dd_evaluate_rhs ( f_operator, b_I, b_G, g_G )
    implicit none
    ! Parameters 
    type(fem_operator_dd), intent(inout) :: f_operator
    type(fem_vector)     , intent(in)    :: b_I
    type(fem_vector)     , intent(in)    :: b_G
    type(fem_vector)     , intent(inout) :: g_G
    
    ! Locals 
    type(fem_vector)                     :: ws_vec_I     ! Work space, interior nodes
    type(fem_vector)                     :: ws_vec_G     ! Work space, interface

    call fem_vector_clone ( b_I, ws_vec_I )
    call fem_vector_clone ( b_G, ws_vec_G )
    
    ! ws_vec_I <- A_II^-1 * b_I   
    call fem_operator_dd_solve_A_II ( f_operator, b_I, ws_vec_I )
    ! call fem_vector_print ( 6, b_I      ) ! DBG:
    ! call fem_vector_print ( 6, ws_vec_I ) ! DBG:
    
    ! ws_vec_G <- A_GI * ws_vec_I  
    call fem_operator_dd_apply_A_GI ( f_operator, ws_vec_I, ws_vec_G )
    
    ! g_G <- b_G 
    call fem_vector_copy ( b_G, g_G )
    
    ! g_G <- g_G - ws_vec_G
    ! call fem_vector_print ( 6, g_G )            ! DBG:
    ! call fem_vector_print ( 6, ws_vec_G )       ! DBG:
    call fem_vector_mxpy ( ws_vec_G, g_G )
    ! call fem_vector_print ( 6, g_G ) ! DBG:
    
    call fem_vector_free ( ws_vec_I )
    call fem_vector_free ( ws_vec_G )
  end subroutine  fem_operator_dd_evaluate_rhs

  ! Computes x_I <- A_II^-1 * (b_I - A_IG * x_G)
  subroutine fem_operator_dd_solve_interior ( f_operator, b_I, x_G, x_I )
    implicit none

    ! Parameters 
    type(fem_operator_dd), intent(inout),target  :: f_operator
    type(fem_vector)  , intent(in)      :: b_I
    type(fem_vector)  , intent(in)      :: x_G
    type(fem_vector)  , intent(inout)   :: x_I

    ! Locals 
    type(fem_vector)  :: ws_vec_I

    call fem_vector_clone ( b_I, ws_vec_I   )

    ! Compute x_I <- A_II^-1 * (b_I - A_IG * x_G)
    ! x_I <- A_IG * x_G  
    call fem_operator_dd_apply_A_IG ( f_operator, x_G, x_I )
    ! call fem_vector_print ( 6, x_G )       ! DBG:
    ! call fem_vector_print ( 6, x_I )       ! DBG:
    ! call fem_matrix_print ( 6, f_operator%A_IG )    ! DBG:
    ! call fem_matrix_print ( 6, f_operator%A_II )
    ! call fem_vector_print ( 6, x_G )       ! DBG:
    ! call fem_vector_print ( 6, ws_vec_I )  ! DBG:
    ! call fem_vector_print ( 6, x_G )        ! DBG:
    ! call fem_vector_print ( 6, ws_vec_I )   ! DBG:
    
    ! Bypass fem_vector* operations and use fem_vector* ones
    ! to avoid dealing with the state member of fem_vector
    ! ws_vec_I <- b_I
    call fem_vector_copy ( b_I, ws_vec_I )
    ! call fem_vector_print ( 6, b_I ) ! DBG:
    ! call fem_vector_print ( 6, ws_vec_I ) ! DBG: 
    
    ! ws_vec_I <- ws_vec_I - x_I  
    ! call fem_vector_print ( 6, ws_vec_I ) ! DBG:
    ! call fem_vector_print ( 6, x_I      ) ! DBG:  
    call fem_vector_mxpy  ( x_I, ws_vec_I )
    ! call fem_vector_print ( 6, ws_vec_I ) ! DBG: 
    ! call fem_vector_print ( 6, x_I ) ! DBG: 
    
    ! (c) x_I  <- A_II^-1 *  ws_vec_I
    call fem_operator_dd_solve_A_II ( f_operator, ws_vec_I, x_I )
    ! call fem_vector_print ( 6, ws_vec_I ) ! DBG:
    ! call fem_vector_print ( 6, x_I )      ! DBG:
    
    call fem_vector_free ( ws_vec_I   )
    
    ! call fem_vector_print ( 6, x_I ) ! DBG: 
    ! call fem_vector_print ( 6, b_I ) ! DBG:

  end subroutine fem_operator_dd_solve_interior

  subroutine fem_operator_dd_info ( f_operator, me, np )
    implicit none

    ! Parameters 
    type(fem_operator_dd), intent(in)    :: f_operator
    integer              , intent(out)   :: me
    integer              , intent(out)   :: np

    me = 0
    np = 1
  end subroutine fem_operator_dd_info

  ! Computes the product v <- u^T op u restricted to the second block row of op 
  ! Documentation:   http://software.intel.com/sites/products/documentation/hpc/
  !                  compilerpro/en-us/cpp/win/mkl/refman/bla/functn_mkl_dcsrmm.html#functn_mkl_dcsrmm
  subroutine fem_operator_dd_compute_ut_op_u (f_operator, n, u, v, opu)  
#ifdef ENABLE_BLAS       
    use blas77_interfaces
#endif
    implicit none 
    ! Parameters  
    type(fem_operator_dd) , intent(in)            :: f_operator
    integer (ip)          , intent(in)            :: n
    real (rp)             , intent(in)            :: u (f_operator%A_II%gr%nv + f_operator%A_GG%gr%nv, n) 
    real (rp)             , intent(out)           :: v (n,n)
    real (rp)             , optional, intent(out) :: opu (f_operator%A_GG%gr%nv,  n) 
 
    ! Locals
    real(rp), allocatable :: work(:,:)

    call memalloc ( f_operator%A_GG%gr%nv, n, work, __FILE__,__LINE__ )

#ifdef ENABLE_MKL
    ! work <- 0.0*work + 1.0 * A_GG * u_G
    work = 0.0_rp
    if (f_operator%A_GG%symm == symm_true) then
        call mkl_dcsrmm ('N', & 
                         f_operator%A_GG%gr%nv, &
                         n, &
                         f_operator%A_GG%gr%nv, &
                         1.0, &
                         'SUNF', & 
                         f_operator%A_GG%a, &
                         f_operator%A_GG%gr%ja, & 
                         f_operator%A_GG%gr%ia(1), & 
                         f_operator%A_GG%gr%ia(2), & 
                         u(f_operator%A_II%gr%nv+1,1), &
                         f_operator%A_II%gr%nv + f_operator%A_GG%gr%nv, &  
                         0.0, & 
                         work, & 
                         f_operator%A_GG%gr%nv)

    else if (f_operator%A_GG%symm == symm_false) then
        call mkl_dcsrmm ( 'N', & 
                         f_operator%A_GG%gr%nv, &
                         n, &
                         f_operator%A_GG%gr%nv, &
                         1.0, &
                         'GXXF', &  
                         f_operator%A_GG%a, &
                         f_operator%A_GG%gr%ja, & 
                         f_operator%A_GG%gr%ia(1), & 
                         f_operator%A_GG%gr%ia(2), & 
                         u(f_operator%A_II%gr%nv+1,1), &
                         f_operator%A_II%gr%nv + f_operator%A_GG%gr%nv, &  
                         0.0, & 
                         work, & 
                         f_operator%A_GG%gr%nv)
     end if
     
#else
     call fem_matmat ( f_operator%A_GG, & 
                       n, & 
                       f_operator%A_II%gr%nv + f_operator%A_GG%gr%nv, &  
                       u(f_operator%A_II%gr%nv+1,1), & 
                       f_operator%A_GG%gr%nv, &
                       work )
#endif

    if (present(opu)) opu = work
    
#ifdef ENABLE_BLAS
      v = 0.0_rp
      call DGEMM( 'T', &
                  'N', &
                  n, &
                  n, &
                  f_operator%A_GG%gr%nv, &
                  1.0, &
                  u(f_operator%A_II%gr%nv+1,1), &
                  f_operator%A_II%gr%nv + f_operator%A_GG%gr%nv, &
                  work, &
                  f_operator%A_GG%gr%nv, &
                  0.0, &
                  v, &
                  n)
#else
     write (0,*) 'Error: fem_operator_dd_compute_ut_op_u was not compiled with -DENABLE_BLAS.'
     write (0,*) 'Error: You must activate this cpp macro in order to use the BLAS'
     check(1==0)    
#endif


#ifdef ENABLE_MKL
    work = 0.0_rp
    if (f_operator%A_II%symm == symm_true) then
        call mkl_dcsrmm ( 'T', & 
                         f_operator%A_IG%gr%nv, &
                         n, &
                         f_operator%A_IG%gr%nv2, &
                         1.0, &
                         'GXXF', & 
                         f_operator%A_IG%a, &
                         f_operator%A_IG%gr%ja, & 
                         f_operator%A_IG%gr%ia(1), & 
                         f_operator%A_IG%gr%ia(2), & 
                         u, &
                         f_operator%A_II%gr%nv + f_operator%A_GG%gr%nv, &  
                         0.0, & 
                         work, & 
                         f_operator%A_GG%gr%nv)

    else if (f_operator%A_II%symm == symm_false) then
        call mkl_dcsrmm ( 'N', & 
                         f_operator%A_GI%gr%nv, &
                         n, &
                         f_operator%A_GI%gr%nv2, &
                         1.0, &
                         'GXXF', &  
                         f_operator%A_GI%a, &
                         f_operator%A_GI%gr%ja, & 
                         f_operator%A_GI%gr%ia(1), & 
                         f_operator%A_GI%gr%ia(2), & 
                         u, &
                         f_operator%A_II%gr%nv + f_operator%A_GG%gr%nv, &  
                         0.0, & 
                         work, & 
                         f_operator%A_GG%gr%nv)
    end if

#else

    if (f_operator%A_II%symm == symm_true) then
       call fem_matmat_trans ( f_operator%A_IG, & 
                               n, & 
                               f_operator%A_II%gr%nv + f_operator%A_GG%gr%nv, &  
                               u, & 
                               f_operator%A_GG%gr%nv, &
                               work )
!!$        call mkl_dcsrmm ( 'T', & 
!!$                         f_operator%A_IG%gr%nv, &
!!$                         n, &
!!$                         f_operator%A_IG%gr%nv2, &
!!$                         1.0, &
!!$                         'GXXF', & 
!!$                         f_operator%A_IG%a, &
!!$                         f_operator%A_IG%gr%ja, & 
!!$                         f_operator%A_IG%gr%ia(1), & 
!!$                         f_operator%A_IG%gr%ia(2), & 
!!$                         u, &
!!$                         f_operator%A_II%gr%nv + f_operator%A_GG%gr%nv, &  
!!$                         0.0, & 
!!$                         work, & 
!!$                         f_operator%A_GG%gr%nv)
    else if (f_operator%A_II%symm == symm_false) then
       call fem_matmat ( f_operator%A_GI, & 
                         n, & 
                         f_operator%A_II%gr%nv + f_operator%A_GG%gr%nv, &  
                         u, & 
                         f_operator%A_GG%gr%nv, &
                         work )
    end if
     ! call fem_matmat ( f_operator%A_GI, n, u(1:f_operator%A_II%gr%nv,:), work )
     ! write (0,*) 'Error: fem_operator_dd_compute_ut_op_u was not compiled with -DENABLE_MKL.'
     ! write (0,*) 'Error: You must activate this cpp macro in order to use the Sparse BLAS in MKL'
     ! check(1==0) 
#endif

    if (present(opu)) opu = opu + work

#ifdef ENABLE_BLAS
      call DGEMM( 'T', &
                  'N', &
                  n, &
                  n, &
                  f_operator%A_GG%gr%nv, &
                  1.0, &
                  u(f_operator%A_II%gr%nv+1,1), &
                  f_operator%A_II%gr%nv + f_operator%A_GG%gr%nv, &
                  work, &
                  f_operator%A_GG%gr%nv, &
                  1.0, &
                  v, &
                  n)
#else
      write (0,*) 'Error: fem_operator_dd_compute_ut_op_u was not compiled with -DENABLE_BLAS.'
      write (0,*) 'Error: You must activate this cpp macro in order to use the BLAS'
      check(1==0)    
#endif

    call memfree ( work,__FILE__,__LINE__)
    ! write (*,*) 'XXX', v
   end subroutine fem_operator_dd_compute_ut_op_u

   ! Computes Y <- beta*Y + alfa*A_IG * X
   subroutine fem_operator_dd_apply_A_IG_several_rhs ( f_operator, nrhs, alpha, X, ldX, beta, Y, ldY )
     implicit none

     ! Parameters 
     type(fem_operator_dd), intent(in)    :: f_operator
     integer(ip), intent(in)              :: nrhs, ldX, ldY
     real(rp)   , intent(in)              :: alpha, beta
     real(rp)   , intent(in)              :: X (ldX, nrhs)
     real(rp)   , intent(inout)           :: Y (ldY , nrhs)

#ifdef ENABLE_MKL
     call mkl_dcsrmm ( 'N', & 
          f_operator%A_IG%gr%nv, &
          nrhs, &
          f_operator%A_IG%gr%nv2, &
          alpha, &
          'GXXF', &  
          f_operator%A_IG%a, &
          f_operator%A_IG%gr%ja, & 
          f_operator%A_IG%gr%ia(1), & 
          f_operator%A_IG%gr%ia(2), & 
          X, &
          ldX, &  
          beta, & 
          Y, & 
          ldY)
#else
     ! write (0,*) 'Error: fem_operator_dd_apply_A_IG_several_rhs was not compiled with -DENABLE_MKL.'
     ! write (0,*) 'Error: You must activate this cpp macro in order to use the Sparse BLAS in MKL'
     ! check(1==0)

     ! The following lines are not efficient at all.
     ! The scaling of Y, i.e., Y = alpha *Y, should be
     ! part of fem_matmat subroutine because now two 
     ! traversals of Y are required. I expect this piece
     ! of code to be cost-negligible so I did not expend
     ! too much time in here ... 
     assert  ( beta == 0.0 )
     call fem_matmat ( f_operator%A_IG, & 
          nrhs, & 
          ldX, &  
          X, & 
          ldY, &
          Y ) 
     Y(1:f_operator%A_IG%gr%nv,:) = alpha*Y(1:f_operator%A_IG%gr%nv,:)
#endif

  end subroutine fem_operator_dd_apply_A_IG_several_rhs

end module fem_operator_dd_names
