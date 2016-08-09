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
module blas77_precision_names
  implicit none
  private

  ! Real precision kind parameters
  integer, parameter  :: DP = 8 ! kind(1.d0) ! Double precision
  integer, parameter  :: SP = 4 ! kind(1.e0) ! Single precision

  public :: SP, DP
end module blas77_precision_names

module blas77_interfaces_names
  implicit none
  ! private

  interface
     SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
       use blas77_precision_names
       implicit none
       ! .. Scalar Arguments ..
       INTEGER , INTENT(IN)    :: INCX, INCY, N
       ! ..
       ! .. Array Arguments ..
       REAL(DP), INTENT(INOUT) :: DX(*), DY(*)
     END SUBROUTINE DSWAP

     SUBROUTINE SSWAP(N,DX,INCX,DY,INCY)
       use blas77_precision_names
       implicit none   
       ! .. Scalar Arguments ..
       INTEGER , INTENT(IN)    :: INCX, INCY, N
       ! ..
       ! .. Array Arguments ..
       REAL(SP), INTENT(INOUT) :: DX(*), DY(*)
     END SUBROUTINE SSWAP

     SUBROUTINE DSCAL(N, DA, DX, INCX)
       use blas77_precision_names
       implicit none
       ! .. Scalar Arguments ..
       REAL(DP), INTENT(IN)    :: DA
       INTEGER , INTENT(IN)    :: INCX, N
       ! ..
       ! .. Array Arguments ..
       REAL(DP), INTENT(INOUT) :: DX(*)
     END SUBROUTINE DSCAL

     SUBROUTINE SSCAL(N, DA, DX, INCX)
       use blas77_precision_names
       implicit none   
       ! .. Scalar Arguments ..
       REAL(SP), INTENT(IN)    :: DA
       INTEGER , INTENT(IN)    :: INCX, N
       ! ..
       ! .. Array Arguments ..
       REAL(SP), INTENT(INOUT) :: DX(*)
     END SUBROUTINE SSCAL

     SUBROUTINE DCOPY(N, DX, INCX, DY, INCY)
       use blas77_precision_names
       implicit none
       ! .. Scalar Arguments ..
       INTEGER , INTENT(IN) :: INCX, INCY, N
       ! ..
       ! .. Array Arguments ..
       REAL(DP), INTENT(IN)    :: DX(*)
       REAL(DP), INTENT(INOUT) :: DY(*)
     END SUBROUTINE DCOPY

     SUBROUTINE SCOPY(N, DX, INCX, DY, INCY)
       use blas77_precision_names
       implicit none   
       ! .. Scalar Arguments ..
       INTEGER , INTENT(IN) :: INCX, INCY, N
       ! ..
       ! .. Array Arguments ..
       REAL(SP), INTENT(IN)    :: DX(*)
       REAL(SP), INTENT(INOUT) :: DY(*)
     END SUBROUTINE SCOPY

     SUBROUTINE DAXPY(N, DA, DX, INCX, DY, INCY)
       use blas77_precision_names
       implicit none   
       ! .. Scalar Arguments ..
       REAL(DP), INTENT(IN) :: DA
       INTEGER , INTENT(IN) :: INCX, INCY, N
       ! ..
       ! .. Array Arguments ..
       REAL(DP), INTENT(IN)    :: DX(*)
       REAL(DP), INTENT(INOUT) :: DY(*)
     END SUBROUTINE DAXPY

     SUBROUTINE SAXPY(N, DA, DX, INCX, DY, INCY)
use blas77_precision_names
       implicit none
       ! .. Scalar Arguments ..
       REAL(SP), INTENT(IN) :: DA
       INTEGER , INTENT(IN) :: INCX, INCY, N
       ! ..
       ! .. Array Arguments ..
       REAL(SP), INTENT(IN)    :: DX(*)
       REAL(SP), INTENT(INOUT) :: DY(*)
     END SUBROUTINE SAXPY

     FUNCTION DDOT(N, DX, INCX, DY, INCY)
       use blas77_precision_names
       implicit none   
       ! .. Scalar Arguments ..
       INTEGER , INTENT(IN) :: INCX, INCY, N
       REAL(DP)             :: DDOT
       ! ..
       ! .. Array Arguments ..
       REAL(DP), INTENT(IN)    :: DX(*)
       REAL(DP), INTENT(IN)    :: DY(*)    
     END FUNCTION DDOT

     FUNCTION SDOT(N, DX, INCX, DY, INCY)
       use blas77_precision_names
       implicit none
       ! .. Scalar Arguments ..
       INTEGER , INTENT(IN) :: INCX, INCY, N
       REAL(SP)             :: SDOT
       ! ..
       ! .. Array Arguments ..
       REAL(SP), INTENT(IN)    :: DX(*)
       REAL(SP), INTENT(IN)    :: DY(*)    
     END FUNCTION SDOT

     FUNCTION DNRM2(N,X,INCX)
       use blas77_precision_names
       implicit none   
       ! .. Scalar Arguments ..
       INTEGER , INTENT(IN) :: INCX, N
       REAL(DP)             :: DNRM2
       ! ..
       ! .. Array Arguments ..
       REAL(DP), INTENT(IN)    :: X(*)
     END FUNCTION DNRM2

     FUNCTION SNRM2(N,X,INCX)
       use blas77_precision_names
       implicit none   

       ! .. Scalar Arguments ..
       INTEGER , INTENT(IN) :: INCX, N
       REAL(DP)             :: SNRM2
       ! ..
       ! .. Array Arguments ..
       REAL(SP), INTENT(IN)    :: X(*)
     END FUNCTION SNRM2

     SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
       use blas77_precision_names
       implicit none  

       ! .. Scalar Arguments ..
       REAL(DP) , INTENT(IN)  :: ALPHA,BETA
       INTEGER  , INTENT(IN)  :: INCX, INCY, LDA, M, N
       CHARACTER, INTENT(IN)  :: TRANS

       ! ..
       ! .. Array Arguments ..
       REAL(DP), INTENT(IN)   :: A(LDA,*), X(*)
       REAL(DP), INTENT(INOUT):: Y(*)
     END SUBROUTINE DGEMV

     SUBROUTINE SGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
       use blas77_precision_names
       implicit none  

       ! .. Scalar Arguments ..
       REAL(SP) , INTENT(IN)  :: ALPHA,BETA
       INTEGER  , INTENT(IN)  :: INCX, INCY, LDA, M, N
       CHARACTER, INTENT(IN)  :: TRANS

       ! ..
       ! .. Array Arguments ..
       REAL(SP), INTENT(IN)   :: A(LDA,*), X(*)
       REAL(SP), INTENT(INOUT):: Y(*)
     END SUBROUTINE SGEMV

     SUBROUTINE DTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
       use blas77_precision_names
       implicit none

       ! .. Scalar Arguments ..
       INTEGER  , INTENT(IN) :: INCX, LDA, N
       CHARACTER, INTENT(IN) :: DIAG, TRANS, UPLO

       ! ..
       ! .. Array Arguments ..
       REAL(DP) , INTENT(IN)   :: A(LDA,*)
       REAL(DP) , INTENT(INOUT):: X(*)
     END SUBROUTINE DTRSV

     SUBROUTINE STRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
       use blas77_precision_names
       implicit none

       ! .. Scalar Arguments ..
       INTEGER  , INTENT(IN) :: INCX, LDA, N
       CHARACTER, INTENT(IN) :: DIAG, TRANS, UPLO

       ! ..
       ! .. Array Arguments ..
       REAL(SP) , INTENT(IN)   :: A(LDA,*)
       REAL(SP) , INTENT(INOUT):: X(*)
     END SUBROUTINE STRSV

     SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
       use blas77_precision_names
       implicit none
       ! .. Scalar Arguments ..
       REAL(DP) , INTENT(IN)   :: ALPHA, BETA
       INTEGER  , INTENT(IN)   :: K, LDA, LDB, LDC, M, N
       CHARACTER, INTENT(IN)   :: TRANSA,TRANSB
       ! ..
       ! .. Array Arguments ..
       REAL(DP), INTENT(IN)    :: A(LDA,*), B(LDB,*)
       REAL(DP), INTENT(INOUT) :: C(LDC,*)
     END SUBROUTINE DGEMM

     SUBROUTINE SGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
       use blas77_precision_names
       implicit none
       ! .. Scalar Arguments ..
       REAL(SP) , INTENT(IN)   :: ALPHA, BETA
       INTEGER  , INTENT(IN)   :: K, LDA, LDB, LDC, M, N
       CHARACTER, INTENT(IN)   :: TRANSA,TRANSB
       ! ..
       ! .. Array Arguments ..
       REAL(SP), INTENT(IN)    :: A(LDA,*), B(LDB,*)
       REAL(SP), INTENT(INOUT) :: C(LDC,*)
     END SUBROUTINE SGEMM

     SUBROUTINE DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
       use blas77_precision_names
       implicit none
       ! *     .. Scalar Arguments ..
       REAL(DP) , INTENT(IN)    :: ALPHA
       INTEGER  , INTENT(IN)    :: LDA, LDB, M, N
       CHARACTER, INTENT(IN)    :: DIAG, SIDE, TRANSA, UPLO
       ! *     ..
       ! *     .. Array Arguments ..
       REAL(DP),  INTENT(IN)    :: A(LDA,*)
       REAL(DP),  INTENT(INOUT) :: B(LDB,*)
     END SUBROUTINE DTRSM

     SUBROUTINE STRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
       use blas77_precision_names
       implicit none
       ! *     .. Scalar Arguments ..
       REAL(SP) , INTENT(IN)    :: ALPHA
       INTEGER  , INTENT(IN)    :: LDA, LDB, M, N
       CHARACTER, INTENT(IN)    :: DIAG, SIDE, TRANSA, UPLO
       ! *     ..
       ! *     .. Array Arguments ..
       REAL(SP),  INTENT(IN)    :: A(LDA,*)
       REAL(SP),  INTENT(INOUT) :: B(LDB,*)
     END SUBROUTINE STRSM
  end interface

end module blas77_interfaces_names

module lapack77_interfaces_names
  implicit none
  ! private

  interface
     SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
       use blas77_precision_names
       implicit none

       ! *     .. Scalar Arguments ..
       INTEGER , INTENT(IN)      :: LDA, M, N
       INTEGER , INTENT(OUT)     :: INFO

       ! *     ..
       ! *     .. Array Arguments ..
       INTEGER , INTENT(OUT)    :: IPIV ( * )
       REAL(DP), INTENT(INOUT)  :: A ( LDA, * )
     END SUBROUTINE DGETRF

     SUBROUTINE DGETRS(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)
       use blas77_precision_names
       implicit none

       ! *     .. Scalar Arguments ..
       CHARACTER, INTENT(IN)     :: TRANS
       INTEGER  , INTENT(IN)     :: LDA, LDB, N, NRHS
       INTEGER  , INTENT(OUT)    :: INFO
    
       ! *     ..
       ! *     .. Array Arguments ..
       INTEGER , INTENT(IN)    :: IPIV ( * )
       REAL(DP), INTENT(IN)    :: A ( LDA, * )
       REAL(DP), INTENT(INOUT) :: B ( LDB, * )
     END SUBROUTINE DGETRS 

     SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
       use blas77_precision_names
       implicit none  
       ! .. Scalar Arguments ..
       CHARACTER, INTENT(IN)     :: UPLO
       INTEGER  , INTENT(IN)     :: LDA, N
       INTEGER  , INTENT(OUT)    :: INFO
       ! *     ..
       !       .. Array Arguments ..
       REAL(DP) , INTENT(INOUT)  :: A( LDA, * ) 
     END SUBROUTINE DPOTRF

     INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
       use blas77_precision_names
       ! .. Scalar Arguments ..
       CHARACTER*(*), INTENT(IN) :: NAME, OPTS
       INTEGER      , INTENT(IN) :: ISPEC, N1, N2, N3, N4  
     END FUNCTION ILAENV

     SUBROUTINE DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
       use blas77_precision_names
       ! .. Scalar Arguments ..
       CHARACTER, INTENT(IN)  :: UPLO
       INTEGER  , INTENT(IN)  :: LDA, LWORK, N
       INTEGER  , INTENT(OUT) :: INFO

       ! .. Array Arguments ..
       INTEGER , INTENT(OUT)   :: IPIV(*)
       REAL(DP), INTENT(INOUT) :: A(LDA,*), WORK(LWORK)
     END SUBROUTINE DSYTRF

     SUBROUTINE DSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
       use blas77_precision_names
       ! .. Scalar Arguments ..
       CHARACTER, INTENT(IN)  :: UPLO
       INTEGER  , INTENT(IN)  :: LDA, LDB, N, NRHS
       INTEGER  , INTENT(OUT) :: INFO

       ! .. Array Arguments ..
       INTEGER , INTENT(IN)    :: IPIV( * )
       REAL(DP), INTENT(IN)    :: A( LDA, * )
       REAL(DP), INTENT(INOUT) :: B( LDB, * )
     END SUBROUTINE DSYTRS

  end interface

end module lapack77_interfaces_names
