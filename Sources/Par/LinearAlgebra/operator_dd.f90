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
module operator_dd_names

  ! Abstract modules
  use abstract_solver_names

  ! Serial modules
  use types_names
  use memor_names
  use graph_names
  use matrix_names
  use matvec_names
  use vector_names
  use preconditioner_names
  use serial_environment_names
  
  ! Parallel modules
  use dof_distribution_names

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

  type operator_dd_t
     type ( graph_t )  :: A_II_gr 
     type ( graph_t )  :: A_IG_gr
     type ( graph_t )  :: A_GI_gr
     type ( graph_t )  :: A_GG_gr

     type ( matrix_t )  :: A_II 
     type ( matrix_t )  :: A_IG 
     type ( matrix_t )  :: A_GI
     type ( matrix_t )  :: A_GG

     integer (ip)             :: symm
     type ( preconditioner_t )     :: M_II

     type(solver_control_t)    , pointer :: spars
     logical                           :: spars_allocated
     type(preconditioner_params_t), pointer :: ppars
     logical                           :: ppars_allocated

     type(dof_distribution_t) , pointer  :: dof_dist => NULL()
  end type operator_dd_t

  ! Types
  public :: operator_dd_t


  public :: operator_dd_create, operator_dd_free,             &  
            operator_dd_ass_struct, operator_dd_fill_val,     & 
            operator_dd_apply, operator_dd_solve_A_II,        &
            operator_dd_apply_A_IG, operator_dd_apply_A_GI,   &
            operator_dd_apply_A_GG, operator_dd_evaluate_rhs, &
            operator_dd_solve_interior, operator_dd_info,     &
            operator_dd_compute_ut_op_u, operator_dd_apply_A_IG_several_rhs, &
            operator_dd_apply_A_GI_plus_A_GG

contains

  subroutine operator_dd_create ( f_matrix, dof_dist, f_operator, spars, ppars, symm, sign )
    implicit none
    ! Parameters
    type(matrix_t)        , intent(in), target           :: f_matrix
    type(dof_distribution_t)  , intent(in), target           :: dof_dist
    type(operator_dd_t)   , intent(out)                  :: f_operator
    type(solver_control_t)    , intent(in), target, optional :: spars
    type(preconditioner_params_t), intent(in), target, optional :: ppars
    ! With this two optional parameters one may select the
    ! structure and sign of the Schur complement operator
    ! independently of f_matrix
    integer (ip)            , intent(in), optional :: symm
    integer (ip)            , intent(in), optional :: sign

    ! Locals
    type    (matrix_t) :: mat_dum
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
    call preconditioner_create( mat_dum , f_operator%M_II, f_operator%ppars)

  end subroutine operator_dd_create 

  subroutine operator_dd_free ( f_operator, mode )
    implicit none

    ! Parameters
    type(operator_dd_t), intent(inout)             :: f_operator
    integer(ip)       , intent(in)                   :: mode

    ! Locals
    type (vector_t) :: dum 

    if ( mode == free_clean ) then
       call preconditioner_free ( preconditioner_free_clean  , f_operator%M_II )
    else if ( mode == free_struct  ) then
       call preconditioner_free ( preconditioner_free_struct , f_operator%M_II )
    else if ( mode == free_values ) then
       call preconditioner_free ( preconditioner_free_values , f_operator%M_II )
    end if

    ! Free memory associated to the blocks of the operator
    if ( mode == free_values ) then
       call matrix_free ( f_operator%A_II, free_values )
       call matrix_free ( f_operator%A_IG, free_values )
       call matrix_free ( f_operator%A_GG, free_values )
    end if

    if ( mode == free_struct ) then
       call matrix_free ( f_operator%A_II, free_struct )
       call matrix_free ( f_operator%A_IG, free_struct )
       call matrix_free ( f_operator%A_GG, free_struct )

       call graph_free ( f_operator%A_II_gr )
       call graph_free ( f_operator%A_IG_gr )
       call graph_free ( f_operator%A_GG_gr )
    end if

    if ( f_operator%symm == symm_false ) then
       if ( mode == free_values ) then
          call matrix_free( f_operator%A_GI, free_values  )
       end if

       if ( mode == free_struct ) then
          call matrix_free( f_operator%A_GI, free_struct  )
          call graph_free ( f_operator%A_GI_gr ) 
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

  end subroutine operator_dd_free


  !=============================================================================
  subroutine operator_dd_ass_struct ( f_matrix, f_operator )
    implicit none

    ! Parameters 
    type(matrix_t)  , intent(in)                    :: f_matrix
    type(operator_dd_t), intent(inout)              :: f_operator

      ! Split graph of local process into 2x2 block partitioning
    if ( f_operator%symm == symm_false ) then
        call split_graph_I_G ( csr, &
                               f_matrix%gr, & 
                               f_operator%dof_dist, & 
                               G_II=f_operator%A_II_gr, G_IG=f_operator%A_IG_gr, &
                               G_GI=f_operator%A_GI_gr, G_GG=f_operator%A_GG_gr  )
        call matrix_graph(f_operator%A_II_gr, f_operator%A_II)
        call matrix_graph(f_operator%A_IG_gr, f_operator%A_IG)
        call matrix_graph(f_operator%A_GI_gr, f_operator%A_GI)
        call matrix_graph(f_operator%A_GG_gr, f_operator%A_GG)
        ! call graph_print ( 6, f_operator%A_GI%gr )
     else if ( f_operator%symm == symm_true ) then
        call split_graph_I_G_symm ( csr_symm, & 
                                    f_matrix%gr, & 
                                    f_operator%dof_dist, & 
                                    G_II=f_operator%A_II_gr, & 
                                    G_IG=f_operator%A_IG_gr, &
                                    G_GG=f_operator%A_GG_gr )
        call matrix_graph(f_operator%A_II_gr, f_operator%A_II)
        call matrix_graph(f_operator%A_IG_gr, f_operator%A_IG)
        call matrix_graph(f_operator%A_GG_gr, f_operator%A_GG)
     end if

     call preconditioner_symbolic(f_operator%A_II, f_operator%M_II)

     ! call graph_print ( 6, f_operator%A_II%gr ) ! DBG:
     ! call graph_print ( 6, f_operator%A_IG%gr ) ! DBG:
     ! call graph_print ( 6, f_operator%A_GI%gr ) ! DBG:
     ! call graph_print ( 6, f_operator%A_GG%gr ) ! DBG:

  end subroutine operator_dd_ass_struct

  !=============================================================================
  subroutine operator_dd_fill_val ( f_matrix, f_operator ) !, me )
    use stdio_names
    implicit none
    
    ! Parameters
    type(matrix_t)   , intent(in)      :: f_matrix
    type(operator_dd_t), intent(inout) :: f_operator
    ! integer(ip)          , intent(in)            :: me

    ! Locals
    ! integer(ip) :: lunou

    ! Split matrix of local process into 2x2 block partitioning
    if ( f_operator%symm == symm_false ) then
        call split_matrix_I_G ( f_operator%symm,                   &
                                f_matrix,                          &
                                f_operator%dof_dist,               & 
                                f_operator%A_II, f_operator%A_IG,  &
                                f_operator%A_GI, f_operator%A_GG )

     else if ( f_operator%symm == symm_true ) then
        call split_matrix_I_G( f_operator%symm,                            &
                               f_matrix,                                   & 
                               f_operator%dof_dist,                        & 
                               A_II=f_operator%A_II, A_IG=f_operator%A_IG, &  
                               A_GG=f_operator%A_GG )
     end if

     ! lunou = io_open ( trim('matrix' // trim(ch(me)) // trim('.') // 'mtx' ), 'write')
     ! call matrix_print_matrix_market ( lunou, f_operator%A_II )
     ! call io_close (lunou)

     ! call preconditioner_numeric(f_operator%A_II, f_operator%M_II)
     call preconditioner_numeric(f_operator%M_II)

     ! call matrix_print ( 6, f_matrix )  ! DBG:
     ! call matrix_print ( 6, f_operator%A_II )    ! DBG:
     ! call matrix_print ( 6, f_operator%A_IG )    ! DBG:
     ! call matrix_print ( 6, f_operator%A_GI )    ! DBG: 
     ! call matrix_print ( 6, f_operator%A_GG )    ! DBG:
  end subroutine operator_dd_fill_val

  !=============================================================================
  ! Computes y_G -< S_G x_G, where S_G the local Schur
  ! complement is given by S_G <- A_GG - A_G_I A_II^-1 A_IG
  ! VERY IMPORTANT: operator_dd_apply is well-defined if and only if x_G and 
  !                 y_G are vectors on the interface
  !=============================================================================  
  subroutine operator_dd_apply ( f_operator, x_G, y_G )
    implicit none

    ! Parameters 
    type(operator_dd_t), intent(in)     :: f_operator
    type(vector_t)     , intent(in)     :: x_G
    type(vector_t)     , intent(inout)  :: y_G

    ! Locals 
    type(vector_t)       :: ws_vec       ! Local workspace vector_t

    type(vector_t)       :: ws_vec_I     ! View of the interior nodes
    type(vector_t)       :: ws_vec_I_2   ! Extra ws required by pardiso
    type(vector_t)       :: ws_vec_G     ! View of the interfaces
                                           
    integer(ip)            :: ni, nl

 
    ni = f_operator%dof_dist%ni   ! # of interior nodes
    nl = f_operator%dof_dist%nl   ! # of local nodes

    ! Allocate space for ws_vec
    call vector_alloc ( nl,   ws_vec)
 
    ! Create a view of ws_vec to point to the interfaces
    call vector_create_view ( ws_vec,    1, ni, ws_vec_I )
    call vector_clone ( ws_vec_I, ws_vec_I_2 )

    ! Create a view of ws_vec to point to the interior nodes
    call vector_create_view ( ws_vec, ni+1, nl, ws_vec_G )

    ! Perform the following computations:
    ! (a) y_G <- A_GG * x_G
    call operator_dd_apply_A_GG ( f_operator, x_G, y_G )
    ! call matrix_print (6, f_operator%A_GG) ! DBG:
    ! call vector_print (6,x_G)     ! DBG:
    ! call vector_print (6,y_G)     ! DBG:

    ! (b) ws_vec_I_2  <- A_IG * x_G
    call operator_dd_apply_A_IG ( f_operator, x_G, ws_vec_I_2 )
    ! call vector_print (6, x_G)         ! DBG:
    ! call vector_print (6, ws_vec_I_2)  ! DBG:

    ! (c) ws_vec_I  <- A_II^-1 *  ws_vec_I_2
    call operator_dd_solve_A_II ( f_operator, ws_vec_I_2, ws_vec_I )
    ! call vector_print (6, ws_vec_I_2)     ! DBG:
    ! call vector_print (6, ws_vec_I)       ! DBG: 

    ! (d) ws_vec_G <- A_GI * ws_vec_I
    call operator_dd_apply_A_GI ( f_operator, ws_vec_I, ws_vec_G )
    ! call matrix_print (6, f_operator%A_IG) ! DBG:
    ! call vector_print (6, ws_vec_I)     ! DBG: 
    ! call vector_print (6, ws_vec_G)     ! DBG:

    ! (e) y_G <- y_G - ws_vec_G 
    call vector_mxpy( ws_vec_G, y_G )

    ! Free memory of local worksize vector
    call vector_free ( ws_vec )
    call vector_free ( ws_vec_I_2 )

  end subroutine operator_dd_apply 

  
  ! Computes y_I <- A_II^-1 * x_I
  subroutine operator_dd_solve_A_II ( f_operator, x_I, y_I )
    implicit none

    ! Parameters 
    type(operator_dd_t), intent(in)    :: f_operator
    type(vector_t)     , intent(in)    :: x_I
    type(vector_t)     , intent(inout) :: y_I

    ! Locals
    type(serial_environment_t) :: senv

    call abstract_solve(f_operator%A_II, f_operator%M_II, x_I, y_I, f_operator%spars, senv) 

  end subroutine operator_dd_solve_A_II

  ! Computes y_I <- A_IG * x_G
  subroutine operator_dd_apply_A_IG ( f_operator, x_G, y_I )
    implicit none

    ! Parameters 
    type(operator_dd_t), intent(in)       :: f_operator
    type(vector_t   )  , intent(in)       :: x_G
    type(vector_t   )  , intent(inout)    :: y_I
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
    call matvec_csr ( f_operator%A_IG%gr%nv, &
                    & f_operator%A_IG%gr%nv2,&
                    & f_operator%A_IG%gr%ia, &
                    & f_operator%A_IG%gr%ja, &
                    & f_operator%A_IG%a,     &
                    & x_G%b,        & 
                    & y_I%b  )
#endif

  end subroutine operator_dd_apply_A_IG

  ! Computes y_G < A_GI * x_I
  subroutine operator_dd_apply_A_GI ( f_operator, x_I, y_G )
    implicit none

    ! Parameters 
    type(operator_dd_t), intent(in)    :: f_operator
    type(vector_t)     , intent(in)    :: x_I
    type(vector_t)     , intent(inout) :: y_G

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
       call matvec_csr ( f_operator%A_GI%gr%nv, &
                         f_operator%A_GI%gr%nv2,&
                         f_operator%A_GI%gr%ia, &
                         f_operator%A_GI%gr%ja, &
                         f_operator%A_GI%a,     &
                         x_I%b,        & ! x_I
                         y_G%b         ) ! y_G
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
       call matvec_csr_trans ( f_operator%A_IG%gr%nv, &
                    &          f_operator%A_IG%gr%nv2,&
                    &          f_operator%A_IG%gr%ia, &
                    &          f_operator%A_IG%gr%ja, &
                    &          f_operator%A_IG%a,     &
                    &          x_I%b,        & ! x_I 
                    &          y_G%b         ) ! y_G 
#endif
    end if
  end subroutine operator_dd_apply_A_GI


  ! Computes y_G < (A_GI*x_I + A_GG*x_G) 
  subroutine operator_dd_apply_A_GI_plus_A_GG ( f_operator, x, y_G )
    implicit none

    ! Parameters 
    type(operator_dd_t), intent(inout)  :: f_operator
    type(vector_t)     , intent(in)     :: x
    type(vector_t)     , intent(inout)  :: y_G
    type(vector_t)                      :: x_I, x_G

   call vector_create_view ( x, 1, f_operator%A_II%gr%nv, x_I )
   call vector_create_view ( x, f_operator%A_II%gr%nv+1, f_operator%A_II%gr%nv+f_operator%A_GG%gr%nv, x_G )

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
      call matvec_csr ( f_operator%A_GI%gr%nv, &
                        f_operator%A_GI%gr%nv2,&
                        f_operator%A_GI%gr%ia, &
                        f_operator%A_GI%gr%ja, &
                        f_operator%A_GI%a,     &
                        x_I%b,        & ! x_I
                        y_G%b         ) ! y_G
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
      call matvec_csr ( f_operator%A_IG%gr%nv, &
                   &    f_operator%A_IG%gr%nv2,&
                   &    f_operator%A_IG%gr%ia, &
                   &    f_operator%A_IG%gr%ja, &
                   &    f_operator%A_IG%a,     &
                   &    x_I%b,        & ! x_I 
                   &    y_G%b         ) ! y_G 
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
           ! call vector_print ( 6, y_G ) ! DBG:
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
        ! call vector_print ( 6, y_G ) ! DBG:
     end if

  end subroutine operator_dd_apply_A_GI_plus_A_GG




  ! Computes y_G < A_GG * x_G
  subroutine operator_dd_apply_A_GG ( f_operator, x_G, y_G )
    implicit none
    ! Parameters 
    type(operator_dd_t), intent(in)    :: f_operator
    type(vector_t)     , intent(in)    :: x_G
    type(vector_t)     , intent(inout) :: y_G

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
       call matvec_csr ( f_operator%A_GG%gr%nv,  &
                         f_operator%A_GG%gr%nv2, &
                         f_operator%A_GG%gr%ia,  &
                         f_operator%A_GG%gr%ja,  &
                         f_operator%A_GG%a,      &
                         x_G%b,         & 
                         y_G%b          )
#endif
       ! call vector_print ( 6, y_G ) ! DBG:
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
       call matvec_csr_symm ( f_operator%A_GG%gr%nv,  &
                              f_operator%A_GG%gr%nv2, &
                              f_operator%A_GG%gr%ia,  &
                              f_operator%A_GG%gr%ja,  &
                              f_operator%A_GG%a,      &
                              x_G%b,         & 
                              y_G%b          )
#endif
       ! call vector_print ( 6, y_G ) ! DBG:
    end if
  end subroutine  operator_dd_apply_A_GG

  ! Computes g_G < b_G - A_GI * A_II^-1 b_I
  subroutine operator_dd_evaluate_rhs ( f_operator, b_I, b_G, g_G )
    implicit none
    ! Parameters 
    type(operator_dd_t), intent(inout) :: f_operator
    type(vector_t)     , intent(in)    :: b_I
    type(vector_t)     , intent(in)    :: b_G
    type(vector_t)     , intent(inout) :: g_G
    
    ! Locals 
    type(vector_t)                     :: ws_vec_I     ! Work space, interior nodes
    type(vector_t)                     :: ws_vec_G     ! Work space, interface

    call vector_clone ( b_I, ws_vec_I )
    call vector_clone ( b_G, ws_vec_G )
    
    ! ws_vec_I <- A_II^-1 * b_I   
    call operator_dd_solve_A_II ( f_operator, b_I, ws_vec_I )
    ! call vector_print ( 6, b_I      ) ! DBG:
    ! call vector_print ( 6, ws_vec_I ) ! DBG:
    
    ! ws_vec_G <- A_GI * ws_vec_I  
    call operator_dd_apply_A_GI ( f_operator, ws_vec_I, ws_vec_G )
    
    ! g_G <- b_G 
    call vector_copy ( b_G, g_G )
    
    ! g_G <- g_G - ws_vec_G
    ! call vector_print ( 6, g_G )            ! DBG:
    ! call vector_print ( 6, ws_vec_G )       ! DBG:
    call vector_mxpy ( ws_vec_G, g_G )
    ! call vector_print ( 6, g_G ) ! DBG:
    
    call vector_free ( ws_vec_I )
    call vector_free ( ws_vec_G )
  end subroutine  operator_dd_evaluate_rhs

  ! Computes x_I <- A_II^-1 * (b_I - A_IG * x_G)
  subroutine operator_dd_solve_interior ( f_operator, b_I, x_G, x_I )
    implicit none

    ! Parameters 
    type(operator_dd_t), intent(inout),target  :: f_operator
    type(vector_t)  , intent(in)      :: b_I
    type(vector_t)  , intent(in)      :: x_G
    type(vector_t)  , intent(inout)   :: x_I

    ! Locals 
    type(vector_t)  :: ws_vec_I

    call vector_clone ( b_I, ws_vec_I   )

    ! Compute x_I <- A_II^-1 * (b_I - A_IG * x_G)
    ! x_I <- A_IG * x_G  
    call operator_dd_apply_A_IG ( f_operator, x_G, x_I )
    ! call vector_print ( 6, x_G )       ! DBG:
    ! call vector_print ( 6, x_I )       ! DBG:
    ! call matrix_print ( 6, f_operator%A_IG )    ! DBG:
    ! call matrix_print ( 6, f_operator%A_II )
    ! call vector_print ( 6, x_G )       ! DBG:
    ! call vector_print ( 6, ws_vec_I )  ! DBG:
    ! call vector_print ( 6, x_G )        ! DBG:
    ! call vector_print ( 6, ws_vec_I )   ! DBG:
    
    ! Bypass vector* operations and use vector* ones
    ! to avoid dealing with the state member of vector
    ! ws_vec_I <- b_I
    call vector_copy ( b_I, ws_vec_I )
    ! call vector_print ( 6, b_I ) ! DBG:
    ! call vector_print ( 6, ws_vec_I ) ! DBG: 
    
    ! ws_vec_I <- ws_vec_I - x_I  
    ! call vector_print ( 6, ws_vec_I ) ! DBG:
    ! call vector_print ( 6, x_I      ) ! DBG:  
    call vector_mxpy  ( x_I, ws_vec_I )
    ! call vector_print ( 6, ws_vec_I ) ! DBG: 
    ! call vector_print ( 6, x_I ) ! DBG: 
    
    ! (c) x_I  <- A_II^-1 *  ws_vec_I
    call operator_dd_solve_A_II ( f_operator, ws_vec_I, x_I )
    ! call vector_print ( 6, ws_vec_I ) ! DBG:
    ! call vector_print ( 6, x_I )      ! DBG:
    
    call vector_free ( ws_vec_I   )
    
    ! call vector_print ( 6, x_I ) ! DBG: 
    ! call vector_print ( 6, b_I ) ! DBG:

  end subroutine operator_dd_solve_interior

  subroutine operator_dd_info ( f_operator, me, np )
    implicit none

    ! Parameters 
    type(operator_dd_t), intent(in)    :: f_operator
    integer              , intent(out)   :: me
    integer              , intent(out)   :: np

    me = 0
    np = 1
  end subroutine operator_dd_info

  ! Computes the product v <- u^T op u restricted to the second block row of op 
  ! Documentation:   http://software.intel.com/sites/products/documentation/hpc/
  !                  compilerpro/en-us/cpp/win/mkl/refman/bla/functn_mkl_dcsrmm.html#functn_mkl_dcsrmm
  subroutine operator_dd_compute_ut_op_u (f_operator, n, u, v, opu)  
#ifdef ENABLE_BLAS       
use blas77_interfaces_names
#endif
    implicit none 
    ! Parameters  
    type(operator_dd_t) , intent(in)            :: f_operator
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
     call matmat ( f_operator%A_GG, & 
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
     write (0,*) 'Error: operator_dd_compute_ut_op_u was not compiled with -DENABLE_BLAS.'
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
       call matmat_trans ( f_operator%A_IG, & 
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
       call matmat ( f_operator%A_GI, & 
                         n, & 
                         f_operator%A_II%gr%nv + f_operator%A_GG%gr%nv, &  
                         u, & 
                         f_operator%A_GG%gr%nv, &
                         work )
    end if
     ! call matmat ( f_operator%A_GI, n, u(1:f_operator%A_II%gr%nv,:), work )
     ! write (0,*) 'Error: operator_dd_compute_ut_op_u was not compiled with -DENABLE_MKL.'
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
      write (0,*) 'Error: operator_dd_compute_ut_op_u was not compiled with -DENABLE_BLAS.'
      write (0,*) 'Error: You must activate this cpp macro in order to use the BLAS'
      check(1==0)    
#endif

    call memfree ( work,__FILE__,__LINE__)
    ! write (*,*) 'XXX', v
   end subroutine operator_dd_compute_ut_op_u

   ! Computes Y <- beta*Y + alfa*A_IG * X
   subroutine operator_dd_apply_A_IG_several_rhs ( f_operator, nrhs, alpha, X, ldX, beta, Y, ldY )
     implicit none

     ! Parameters 
     type(operator_dd_t), intent(in)    :: f_operator
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
     ! write (0,*) 'Error: operator_dd_apply_A_IG_several_rhs was not compiled with -DENABLE_MKL.'
     ! write (0,*) 'Error: You must activate this cpp macro in order to use the Sparse BLAS in MKL'
     ! check(1==0)

     ! The following lines are not efficient at all.
     ! The scaling of Y, i.e., Y = alpha *Y, should be
     ! part of matmat subroutine because now two 
     ! traversals of Y are required. I expect this piece
     ! of code to be cost-negligible so I did not expend
     ! too much time in here ... 
     assert  ( beta == 0.0 )
     call matmat ( f_operator%A_IG, & 
          nrhs, & 
          ldX, &  
          X, & 
          ldY, &
          Y ) 
     Y(1:f_operator%A_IG%gr%nv,:) = alpha*Y(1:f_operator%A_IG%gr%nv,:)
#endif

  end subroutine operator_dd_apply_A_IG_several_rhs

  subroutine split_matrix_I_G (  output_symm, &
                                 A, dof_dist, A_II, A_IG, A_GI, A_GG  )
    !-----------------------------------------------------------------------
    ! Given a 2x2 interior/interface block partitioning described by the
    ! "dof_dist" input parameter: 
    !      
    !  A = [A_II A_IG]
    !      [A_GI A_GG]
    !
    ! this routine computes A_II, A_IG, A_GI and A_GG given the global 
    ! matrix A (see parameter "grph"). Note that A_II, A_IG, A_GI and 
    ! A_GG are all optional. Depending on whether A is of type csr_mat 
    ! + symm==false or csr_mat + symm==true the following output is 
    ! produced:
    !
    !      - csr_mat + symm=false: A_II, A_IG, A_GI and A_GG are stored 
    !                              in csr_mat + symm=false format 
    !
    !      - csr_mat + symm=true:  A_II, A_IG stored in csr_mat + symm=true 
    !                              and A_IG in csr_mat + symm=false. Asking for 
    !                              A_GI is an error in this case. 
    !
    ! * IMPORTANT NOTE: this routine assumes that gr pointer of A_II, A_IG, 
    !                   A_GI and A_GG is already associated. Otherwise, it 
    !                   does not work.
    !
    ! * TO-DO: Unpack member allocatable arrays of matrix
    !          into explicit size arrays and pass them to auxiliary routines 
    !          for performance reasons
    !-----------------------------------------------------------------------
    implicit none

    ! Parameters
    integer(ip)           , intent(in)                :: output_symm
    type(matrix_t)      , intent(in)                :: A
    type(dof_distribution_t), intent(in)                :: dof_dist 

    type(matrix_t)     , intent(inout), optional   :: A_II
    type(matrix_t)     , intent(inout), optional   :: A_IG
    type(matrix_t)     , intent(inout), optional   :: A_GI
    type(matrix_t)     , intent(inout), optional   :: A_GG

    integer :: ipoing, offset

    ! Locals
    logical     :: csr_mat_symm, csr_mat_unsymm

    logical     :: present_a_ii, present_a_ig, & 
         & present_a_gi, present_a_gg

    integer(ip) :: ni_rows, nb_rows, ni_cols, nb_cols

    csr_mat_symm   = (A%type == csr_mat .and. output_symm == symm_true )
    csr_mat_unsymm = (A%type == csr_mat .and. output_symm == symm_false)

    assert ( csr_mat_symm .or. csr_mat_unsymm )
    assert ( .not. present(A_GI) .or. csr_mat_unsymm )

       ni_rows = dof_dist%ni  
       nb_rows = dof_dist%nb  
       ni_cols = dof_dist%ni  
       nb_cols = dof_dist%nb  

    ! If any inout matrix is present, we are done !
    if ( csr_mat_unsymm ) then 
       if ( .not. (present(A_II) .or. present(A_IG) & 
            & .or. present(A_GI) .or. present(A_GG) ) ) return
    else 
       if ( csr_mat_symm ) then
          if ( .not. (present(A_II) .or. present(A_IG) .or. present(A_GG) ) ) return
       end if
    end if

    present_a_ii = present( A_II )
    present_a_ig = present( A_IG )
    present_a_gi = present( A_GI )
    present_a_gg = present( A_GG )

    if ( present_a_ii ) then
       assert ( associated(A_II%gr) ) 

       A_II%type     = A%type
       A_II%symm     = output_symm


       call memalloc ( A_II%gr%ia(A_II%gr%nv+1)-A_II%gr%ia(1), A_II%a,          __FILE__,__LINE__)       

    end if

    if ( present_a_ig ) then
       assert ( associated(A_IG%gr) ) 
       A_IG%type    = A%type
       A_IG%symm    = symm_false

       call memalloc ( A_IG%gr%ia(A_IG%gr%nv+1)-A_IG%gr%ia(1), A_IG%a,__FILE__,__LINE__ )


    end if

    if ( present_a_gi ) then
       assert ( associated(A_GI%gr) ) 
       A_GI%type    = A%type
       A_GI%symm    = symm_false

          call memalloc ( A_GI%gr%ia(A_GI%gr%nv+1)-A_GI%gr%ia(1), A_GI%a,                       __FILE__,__LINE__ )
       ! write (*,*)  A_GI%gr%ia(A_GI%gr%nv+1)-A_GI%gr%ia(1) ! DBG:
    end if

    if ( present_a_gg ) then

       assert ( associated(A_GG%gr) ) 
       A_GG%type    = A%type
       A_GG%symm    = output_symm

          call memalloc ( A_GG%gr%ia(A_GG%gr%nv+1)-A_GG%gr%ia(1), A_GG%a,                       __FILE__,__LINE__ )

    end if

    ! List values on each row of G_II/G_IG
    if ( present_a_ii .or. present_a_ig ) then

          do  ipoing=1, ni_rows 
             if ( output_symm == A%symm ) then
                if (present_a_ii) then 
                   A_II%a(  A_II%gr%ia(ipoing):A_II%gr%ia(ipoing+1)-1 ) = &  
                        A%a(  A%gr%ia(ipoing):A%gr%ia(ipoing)+(A_II%gr%ia(ipoing+1)-A_II%gr%ia(ipoing))-1 )
                end if
             else if ( output_symm == symm_true .and. A%symm == symm_false ) then
                if (present_a_ii) then
                   offset = (A%gr%ia(ipoing+1)- A%gr%ia(ipoing))-(A_IG%gr%ia(ipoing+1)-A_IG%gr%ia(ipoing))-(A_II%gr%ia(ipoing+1)-A_II%gr%ia(ipoing))
                   A_II%a( A_II%gr%ia(ipoing):A_II%gr%ia(ipoing+1)-1 ) = &  
                        A%a( A%gr%ia(ipoing)+offset:A%gr%ia(ipoing)+offset+(A_II%gr%ia(ipoing+1)-A_II%gr%ia(ipoing))-1 )
                end if
             else if ( output_symm == symm_false .and. A%symm == symm_true ) then
                ! Not implemented yet. Trigger an assertion.
                assert ( 1 == 0 )
             end if

             if (present_a_ig) then
                A_IG%a( A_IG%gr%ia(ipoing):A_IG%gr%ia(ipoing+1)-1 ) = &
                     A%a( A%gr%ia(ipoing+1)-(A_IG%gr%ia(ipoing+1)-A_IG%gr%ia(ipoing)):A%gr%ia(ipoing+1)-1)
             end if
          end do
    end if


    ! List values on each row of G_GI/G_GG
    if ( present_a_gi .or. present_a_gg ) then
       do ipoing=ni_rows+1, ni_rows + nb_rows
          if ( present_a_gi ) then
             ! write (*,*) A_GI%gr%ia( ipoing-ni_rows ),  A_GI%gr%ia( ipoing+1-ni_rows )-1                               ! DBG:
             ! write (*,*) A%gr%ia(ipoing), A%gr%ia(ipoing)+(A_GI%gr%ia(ipoing+1-ni_rows)-A_GI%gr%ia(ipoing-ni_rows))-1  ! DBG:
             A_GI%a(  A_GI%gr%ia( ipoing-ni_rows ) : A_GI%gr%ia( ipoing+1-ni_rows )-1 ) = &  
                  A%a ( A%gr%ia(ipoing):A%gr%ia(ipoing)+(A_GI%gr%ia(ipoing+1-ni_rows)-A_GI%gr%ia(ipoing-ni_rows))-1 )
          end if
          if ( present_a_gg ) then
             A_GG%a( A_GG%gr%ia( ipoing-ni_rows ) : A_GG%gr%ia(ipoing+1-ni_rows)-1 )   = &  
                  A%a ( A%gr%ia(ipoing+1)-(A_GG%gr%ia(ipoing+1-ni_rows)-A_GG%gr%ia(ipoing-ni_rows)):A%gr%ia(ipoing+1)-1 )        
          end if
       end do
    end if

  end subroutine split_matrix_I_G

  
  subroutine split_graph_I_G ( output_type, & 
                                              & grph, dof_dist, G_II, G_IG, G_GI, G_GG  )
    !-----------------------------------------------------------------------
    ! Given a 2x2 interior/interface block partitioning described by the
    ! "dof_dist" input parameter: 
    !      
    !  A = [A_II A_IG]
    !      [A_GI A_GG]
    !
    ! this routine computes the graphs associated with A_II, A_IG,
    ! A_GI and A_GG given the graph of the global matrix A (see parameter "grph"). 
    ! Note that G_II, G_IG, G_GI and A_GG are all optional. Depending on whether
    ! "output_type" is of type csr or csr_symm the following output is produced:
    !
    !      - csr:       G_II, G_IG, G_GI and G_GG are stored in csr format 
    !
    !      - csr_symm:  G_II, G_GG stored in csr_symm and G_IG in csr. Asking
    !                   for G_GI is an error in this case. 
    !
    ! * IMPORTANT NOTE: this routine assumes that "grph" has column indices
    !                   listed in increasing order. Otherwise, it does not work.
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    integer(ip)            , intent(in) :: output_type
    type(graph_t)        , intent(in) :: grph
    type(dof_distribution_t) , intent(in) :: dof_dist 

    type(graph_t)     , intent(out), optional   :: G_II
    type(graph_t)     , intent(out), optional   :: G_IG
    type(graph_t)     , intent(out), optional   :: G_GI
    type(graph_t)     , intent(out), optional   :: G_GG


    assert ( grph%type   == csr_symm .or. grph%type   == csr )
    assert ( output_type == csr_symm .or. output_type == csr )
    assert ( .not. present(G_GI) .or. output_type == csr )
    
    ! If any output graph is present, we are done !
    if ( grph%type == csr ) then 
      if ( .not. (present(G_II) .or. present(G_IG) & 
         & .or. present(G_GI) .or. present(G_GG) ) ) return
    else 
        if ( grph%type == csr_symm ) then
          if ( .not. (present(G_II) .or. present(G_IG) .or. present(G_GG) ) ) return
        end if
    end if

    call split_graph_I_G_count_list ( output_type, grph, dof_dist, G_II=G_II, G_IG=G_IG, G_GI=G_GI, G_GG=G_GG  )

  end subroutine split_graph_I_G

  ! Auxiliary routine. TO-DO: Unpack member allocatable arrays of graph
  ! into explicit size arrays and pass them to auxiliary routines for 
  ! performance reasons
  subroutine split_graph_I_G_count_list ( output_type, grph, dof_dist, G_II, G_IG, G_GI, G_GG  )
    implicit none
    ! Parameters
    integer(ip)         , intent(in)              :: output_type
    type(graph_t)     , intent(in)              :: grph
    type(dof_distribution_t) , intent(in)              :: dof_dist 

    type(graph_t)     , intent(out), optional   :: G_II
    type(graph_t)     , intent(out), optional   :: G_IG
    type(graph_t)     , intent(out), optional   :: G_GI
    type(graph_t)     , intent(out), optional   :: G_GG

    integer(ip) :: nz_ii, nz_ig, nz_gi, nz_gg
    integer(ip) :: ni_rows, nb_rows, ni_cols, nb_cols

    integer(ip) :: ipoing, ipoing_neig, pos_neig

    logical     :: present_g_ii, present_g_ig, & 
         & present_g_gi, present_g_gg


    ni_rows = dof_dist%ni 
    nb_rows = dof_dist%nb 
    ni_cols = dof_dist%ni  
    nb_cols = dof_dist%nb  

    present_g_ii = present( G_II )
    present_g_ig = present( G_IG )
    present_g_gi = present( G_GI )
    present_g_gg = present( G_GG )

    if ( present_g_ii ) then
       G_II%nv   = ni_rows
       G_II%nv2  = ni_cols
       G_II%type = output_type
       call memalloc ( G_II%nv+1, G_II%ia,__FILE__,__LINE__)
       G_II%ia(1) = 1
    end if

    if ( present_g_ig ) then
       G_IG%nv   = ni_rows
       G_IG%nv2  = nb_cols
       G_IG%type = csr
       call memalloc ( G_IG%nv+1, G_IG%ia, __FILE__,__LINE__ )
       G_IG%ia(1) = 1
    end if

    if ( present_g_gi ) then
       G_GI%nv   = nb_rows
       G_GI%nv2  = ni_cols
       G_GI%type = csr 
       call memalloc ( G_GI%nv+1, G_GI%ia,__FILE__,__LINE__ )
       G_GI%ia(1) = 1
    end if

    if ( present_g_gg ) then
       G_GG%nv   = nb_rows
       G_GG%nv2  = nb_cols
       G_GG%type = output_type 
       call memalloc ( G_GG%nv+1, G_GG%ia,__FILE__,__LINE__ )
       G_GG%ia(1) = 1
    end if

    nz_ii = 1
    nz_ig = 1
    nz_gi = 1
    nz_gg = 1

    ! List number of nonzeros on each row of G_II/G_IG
    if ( present_g_ii .or. present_g_ig ) then
       do ipoing=1, ni_rows

          if ( grph%type == output_type ) then

             do pos_neig=grph%ia(ipoing), grph%ia(ipoing+1)-1
                ipoing_neig = grph%ja(pos_neig)

                if ( ipoing_neig <= ni_cols ) then
                   nz_ii = nz_ii + 1
                else 
                   nz_ig = nz_ig + 1 
                end if
             end do

          else if ( grph%type == csr      .and. output_type == csr_symm ) then

             do pos_neig=grph%ia(ipoing), grph%ia(ipoing+1)-1
                ipoing_neig = grph%ja(pos_neig)

                if ( ipoing_neig >= ipoing .and. ipoing_neig <= ni_cols ) then
                   nz_ii = nz_ii + 1
                else if ( ipoing_neig > ni_cols ) then
                   nz_ig = nz_ig + 1 
                end if
             end do

          else if ( grph%type == csr_symm .and. output_type == csr ) then
             ! Not implemented yet. Trigger an assertion.
             assert ( 1 == 0 )
          end if

          if ( present_g_ii ) then
             G_II%ia(ipoing+1) = nz_ii
          end if

          if ( present_g_ig ) then
             G_IG%ia(ipoing+1) = nz_ig
          end if

       end do
    end if

    !write (*,*), 'XXX', grph%nv, size(grph%ia), &
    !              dof_dist%ni + dof_dist%nb,  dof_dist%nl ! DBG

    if ( present_g_gi .or. present_g_gg ) then  
       ! List number of nonzeros on each row of G_GI/G_GG)
       do ipoing=ni_rows+1, ni_rows + nb_rows
          if ( grph%type == output_type ) then
             do pos_neig=grph%ia(ipoing), grph%ia(ipoing+1)-1
                ipoing_neig = grph%ja(pos_neig)
                if ( ipoing_neig <= ni_cols ) then
                   nz_gi = nz_gi + 1
                else
                   nz_gg = nz_gg + 1
                end if
             end do
          else if ( grph%type == csr      .and. output_type == csr_symm ) then
             do pos_neig=grph%ia(ipoing), grph%ia(ipoing+1)-1
                ipoing_neig = grph%ja(pos_neig)
                if ( ipoing_neig <= ni_cols ) then
                   nz_gi = nz_gi + 1
                else if ( ipoing_neig >= ipoing ) then
                   nz_gg = nz_gg + 1
                end if
             end do

          else if ( grph%type == csr_symm .and. output_type == csr ) then
             ! Not implemented yet. Trigger an assertion.
             assert ( 1 == 0 )
          end if

          if ( present_g_gi ) then
             G_GI%ia(ipoing+1-ni_rows) = nz_gi
          end if

          if ( present_g_gg ) then
             G_GG%ia(ipoing+1-ni_rows) = nz_gg
          end if

       end do

    end if

    ! write (*,*) nz_ii, nz_ig, nz_gi ! DBG:
    if ( present_g_ii ) then
       ! write (*,'(10i10)') G_II%ia(1:G_II%nv+1)
       call memalloc (nz_ii-1, G_II%ja, __FILE__,__LINE__)
    end if

    if ( present_g_ig ) then
       ! write (*,'(10i10)') G_IG%ia(1:G_IG%nv+1)
       call memalloc (nz_ig-1, G_IG%ja, __FILE__,__LINE__)
    end if

    if ( present_g_gi ) then
       ! write (*,'(10i10)') G_GI%ia(1:G_GI%nv+1)
       call memalloc (nz_gi-1, G_GI%ja, __FILE__,__LINE__)
    end if

    if ( present_g_gg ) then
       ! write (*,'(10i10)') G_GG%ia(1:G_GG%nv+1)
       call memalloc (nz_gg-1, G_GG%ja, __FILE__,__LINE__)
    end if

    nz_ii = 1
    nz_ig = 1
    nz_gi = 1
    nz_gg = 1

    ! List nonzeros on each row of G_II/G_IG
    if ( present_g_ii .or. present_g_ig ) then
       do ipoing=1, ni_rows
          if ( grph%type == output_type ) then
             do pos_neig=grph%ia(ipoing), grph%ia(ipoing+1)-1
                ipoing_neig = grph%ja(pos_neig)
                if ( ipoing_neig <= ni_cols ) then
                   if ( present_g_ii ) then
                      G_II%ja(nz_ii) = ipoing_neig
                      nz_ii = nz_ii + 1 
                   end if
                else 
                   if ( present_g_ig ) then
                      G_IG%ja(nz_ig) = ipoing_neig - ni_cols 
                      nz_ig = nz_ig + 1 
                   end if
                end if
             end do
          else if ( grph%type == csr      .and. output_type == csr_symm ) then
             do pos_neig=grph%ia(ipoing), grph%ia(ipoing+1)-1
                ipoing_neig = grph%ja(pos_neig)
                if ( ipoing_neig >= ipoing .and. ipoing_neig <= ni_cols ) then
                   if ( present_g_ii ) then
                      G_II%ja(nz_ii) = ipoing_neig
                      nz_ii = nz_ii + 1 
                   end if
                else if ( ipoing_neig > ni_cols ) then
                   if ( present_g_ig ) then
                      G_IG%ja(nz_ig) = ipoing_neig - ni_cols 
                      nz_ig = nz_ig + 1 
                   end if
                end if
             end do
          else if ( grph%type == csr_symm .and. output_type == csr ) then
             ! Not implemented yet. Trigger an assertion.
             assert ( 1 == 0 )
          end if
       end do
    end if


    ! List number of nonzeros on each row of G_GI/G_GG
    if ( present_g_gi .or. present_g_gg ) then 
       ! write (*,*) 'XXX', ni_rows, nb_rows, grph%type == output_type
       do ipoing=ni_rows+1, ni_rows + nb_rows
          if ( grph%type == output_type ) then
             ! write (*,*) 'YYY', size(grph%ia), size(grph%ja), ni_cols
             do pos_neig=grph%ia(ipoing), grph%ia(ipoing+1)-1
                ! write (*,*) 'ZZZ', pos_neig
                ipoing_neig = grph%ja(pos_neig)
                if ( ipoing_neig <= ni_cols .and. present_g_gi ) then
                   ! if ( present_g_gi ) then 
                   G_GI%ja( nz_gi ) = ipoing_neig 
                   nz_gi = nz_gi + 1
                   ! end if
                else
                   ! assert ( nz_gg <= size(G_GG%ja) )
                   if ( present_g_gg ) then
                      G_GG%ja( nz_gg ) = ipoing_neig - ni_cols  
                      nz_gg = nz_gg + 1
                   end if
                end if
             end do
          else if ( grph%type == csr      .and. output_type == csr_symm ) then
             do pos_neig=grph%ia(ipoing), grph%ia(ipoing+1)-1
                ipoing_neig = grph%ja(pos_neig)
                if ( ipoing_neig <= ni_cols ) then
                   if ( present_g_gi ) then
                      G_GI%ja( nz_gi ) = ipoing_neig 
                      nz_gi = nz_gi + 1
                   end if
                else if ( ipoing_neig >= ipoing ) then
                   if ( present_g_gg ) then
                      G_GG%ja( nz_gg ) = ipoing_neig - ni_cols  
                      nz_gg = nz_gg + 1
                   end if
                end if
             end do
          else if ( grph%type == csr_symm .and. output_type == csr ) then
             ! Not implemented yet. Trigger an assertion.
             assert ( 1 == 0 )
          end if
       end do
    end if

  end subroutine split_graph_I_G_count_list
  
  subroutine split_graph_I_G_symm ( output_type, & 
                                              & grph, dof_dist, G_II, G_IG, G_GG  )
    !-----------------------------------------------------------------------
    ! Given a 2x2 interior/interface block partitioning described by the
    ! "dof_dist" input parameter: 
    !      
    !  A = [A_II A_IG]
    !      [A_GI A_GG]
    !
    ! this routine computes the graphs associated with A_II, A_IG,
    ! A_GI and A_GG given the graph of the global matrix A (see parameter "grph"). 
    ! Note that G_II, G_IG, G_GI and A_GG are all optional. Depending on whether
    ! "output_type" is of type csr or csr_symm the following output is produced:
    !
    !      - csr:       G_II, G_IG, G_GI and G_GG are stored in csr format 
    !
    !      - csr_symm:  G_II, G_GG stored in csr_symm and G_IG in csr. Asking
    !                   for G_GI is an error in this case. 
    !
    ! * IMPORTANT NOTE: this routine assumes that "grph" has column indices
    !                   listed in increasing order. Otherwise, it does not work.
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    integer(ip)         , intent(in)              :: output_type
    type(graph_t)     , intent(in)              :: grph
    type(dof_distribution_t) , intent(in)              :: dof_dist 

    type(graph_t)     , intent(out) :: G_II
    type(graph_t)     , intent(out) :: G_IG
    type(graph_t)     , intent(out) :: G_GG

    ! assert ( grph%type   == csr_symm )
    assert ( output_type == csr_symm .or. output_type == csr )
    
    call split_graph_I_G_count_list_symm ( output_type, grph, dof_dist, G_II=G_II, G_IG=G_IG, G_GG=G_GG  )

  end subroutine split_graph_I_G_symm

  ! Auxiliary routine. TO-DO: Unpack member allocatable arrays of graph
  ! into explicit size arrays and pass them to auxiliary routines for 
  ! performance reasons
  subroutine split_graph_I_G_count_list_symm ( output_type, grph, dof_dist, G_II, G_IG, G_GG  )
    implicit none
    ! Parameters
    integer(ip)            , intent(in) :: output_type
    type(graph_t)        , intent(in) :: grph
    type(dof_distribution_t) , intent(in) :: dof_dist 

    type(graph_t)     , intent(out)   :: G_II
    type(graph_t)     , intent(out)   :: G_IG
    type(graph_t)     , intent(out)   :: G_GG

    integer(ip) :: nz_ii, nz_ig, nz_gi, nz_gg
    integer(ip) :: ni_rows, nb_rows, ni_cols, nb_cols

    integer(ip) :: ipoing, ipoing_neig, pos_neig

    logical :: tmp

    ni_rows = dof_dist%ni 
    nb_rows = dof_dist%nb 
    ni_cols = dof_dist%ni 
    nb_cols = dof_dist%nb 

    G_II%nv   = ni_rows
    G_II%nv2  = ni_cols
    G_II%type = output_type
    call memalloc ( G_II%nv+1, G_II%ia,__FILE__,__LINE__)
    G_II%ia(1) = 1

    G_IG%nv   = ni_rows
    G_IG%nv2  = nb_cols
    G_IG%type = csr
    call memalloc ( G_IG%nv+1, G_IG%ia, __FILE__,__LINE__ )
    G_IG%ia(1) = 1

    G_GG%nv   = nb_rows
    G_GG%nv2  = nb_cols
    G_GG%type = output_type 
    call memalloc ( G_GG%nv+1, G_GG%ia,__FILE__,__LINE__ )
    G_GG%ia(1) = 1

    nz_ii = 1
    nz_ig = 1
    nz_gi = 1
    nz_gg = 1

    ! List number of nonzeros on each row of G_II/G_IG
    do ipoing=1, ni_rows
       if ( grph%type == output_type ) then
          do pos_neig=grph%ia(ipoing), grph%ia(ipoing+1)-1
             ipoing_neig = grph%ja(pos_neig)

             if ( ipoing_neig <= ni_cols ) then
                nz_ii = nz_ii + 1
             else 
                nz_ig = nz_ig + 1 
             end if
          end do

       else if ( grph%type == csr .and. output_type == csr_symm ) then

          do pos_neig=grph%ia(ipoing), grph%ia(ipoing+1)-1
             ipoing_neig = grph%ja(pos_neig)

             if ( ipoing_neig >= ipoing .and. ipoing_neig <= ni_cols ) then
                nz_ii = nz_ii + 1
             else if ( ipoing_neig > ni_cols ) then
                nz_ig = nz_ig + 1 
             end if
          end do

       else if ( grph%type == csr_symm .and. output_type == csr ) then
          ! Not implemented yet. Trigger an assertion.
          assert ( 1 == 0 )
       end if

       G_II%ia(ipoing+1) = nz_ii
       G_IG%ia(ipoing+1) = nz_ig

    end do

    !write (*,*), 'XXX', grph%nv, size(grph%ia), &
    !              dof_dist%ni + dof_dist%nb,  dof_dist%nl ! DBG

    ! List number of nonzeros on each row of G_GI/G_GG)
    do ipoing=ni_rows+1, ni_rows + nb_rows
       if ( grph%type == output_type ) then
          do pos_neig=grph%ia(ipoing), grph%ia(ipoing+1)-1
             ipoing_neig = grph%ja(pos_neig)
             if ( ipoing_neig <= ni_cols ) then
                nz_gi = nz_gi + 1
             else
                nz_gg = nz_gg + 1
             end if
          end do
       else if ( grph%type == csr      .and. output_type == csr_symm ) then
          do pos_neig=grph%ia(ipoing), grph%ia(ipoing+1)-1
             ipoing_neig = grph%ja(pos_neig)
             if ( ipoing_neig <= ni_cols ) then
                nz_gi = nz_gi + 1
             else if ( ipoing_neig >= ipoing ) then
                nz_gg = nz_gg + 1
             end if
          end do

       else if ( grph%type == csr_symm .and. output_type == csr ) then
          ! Not implemented yet. Trigger an assertion.
          assert ( 1 == 0 )
       end if

       G_GG%ia(ipoing+1-ni_rows) = nz_gg

    end do


    ! write (*,*) nz_ii, nz_ig, nz_gi ! DBG:
    ! write (*,'(10i10)') G_II%ia(1:G_II%nv+1)
    call memalloc (nz_ii-1, G_II%ja, __FILE__,__LINE__)

    ! write (*,'(10i10)') G_IG%ia(1:G_IG%nv+1)
    call memalloc (nz_ig-1, G_IG%ja, __FILE__,__LINE__)

    ! write (*,'(10i10)') G_GG%ia(1:G_GG%nv+1)
    call memalloc (nz_gg-1, G_GG%ja, __FILE__,__LINE__)

    nz_ii = 1
    nz_ig = 1
    nz_gi = 1
    nz_gg = 1

    ! List nonzeros on each row of G_II/G_IG
    do ipoing=1, ni_rows
       if ( grph%type == output_type ) then
          do pos_neig=grph%ia(ipoing), grph%ia(ipoing+1)-1
             ipoing_neig = grph%ja(pos_neig)
             if ( ipoing_neig <= ni_cols ) then
                G_II%ja(nz_ii) = ipoing_neig
                nz_ii = nz_ii + 1 
             else 
                G_IG%ja(nz_ig) = ipoing_neig - ni_cols 
                nz_ig = nz_ig + 1 
             end if
          end do
       else if ( grph%type == csr      .and. output_type == csr_symm ) then
          do pos_neig=grph%ia(ipoing), grph%ia(ipoing+1)-1
             ipoing_neig = grph%ja(pos_neig)
             if ( ipoing_neig >= ipoing .and. ipoing_neig <= ni_cols ) then
                G_II%ja(nz_ii) = ipoing_neig
                nz_ii = nz_ii + 1 
             else if ( ipoing_neig > ni_cols ) then
                G_IG%ja(nz_ig) = ipoing_neig - ni_cols 
                nz_ig = nz_ig + 1 
             end if
          end do
       else if ( grph%type == csr_symm .and. output_type == csr ) then
          ! Not implemented yet. Trigger an assertion.
          assert ( 1 == 0 )
       end if

    end do

    ! call graph_print ( 6, grph )


    ! List number of nonzeros on each row of G_GI/G_GG
    ! write (*,*) 'XXX', ni_rows, nb_rows, grph%type == output_type
    do ipoing=ni_rows+1, ni_rows + nb_rows
       if ( grph%type == output_type ) then
          ! write (*,*) 'YYY', size(grph%ia), size(grph%ja), ni_cols
          do pos_neig=grph%ia(ipoing), grph%ia(ipoing+1)-1
             ! write (*,*) 'ZZZ', pos_neig
             ipoing_neig = grph%ja(pos_neig)
             if ( ipoing_neig > ni_cols ) then
                G_GG%ja( nz_gg ) = ipoing_neig - ni_cols  
                nz_gg = nz_gg + 1
             end if
          end do
       else if ( grph%type == csr      .and. output_type == csr_symm ) then
          do pos_neig=grph%ia(ipoing), grph%ia(ipoing+1)-1
             ipoing_neig = grph%ja(pos_neig)
             if ( ipoing_neig >= ipoing ) then
                G_GG%ja( nz_gg ) = ipoing_neig - ni_cols  
                nz_gg = nz_gg + 1
             end if
          end do
       else if ( grph%type == csr_symm .and. output_type == csr ) then
          ! Not implemented yet. Trigger an assertion.
          assert ( 1 == 0 )
       end if
    end do

  end subroutine split_graph_I_G_count_list_symm


end module operator_dd_names
