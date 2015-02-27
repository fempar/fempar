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
module fem_matrix_partition
  use types
  use memor
  use fem_matrix_names
  use fem_partition_names
  implicit none
# include "debug.i90"
  private

  public :: fem_matrix_split_2x2_partitioning

contains
  subroutine fem_matrix_split_2x2_partitioning (  storage, ndof1, ndof2, output_symm, &
                                                  A, prt, A_II, A_IG, A_GI, A_GG  )
    !-----------------------------------------------------------------------
    ! Given a 2x2 interior/interface block partitioning described by the
    ! "prt" input parameter: 
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
    ! * TO-DO: Unpack member allocatable arrays of fem_matrix
    !          into explicit size arrays and pass them to auxiliary routines 
    !          for performance reasons
    !-----------------------------------------------------------------------
    implicit none

    ! Parameters
    integer(ip)          , intent(in)                :: storage, ndof1, ndof2, output_symm
    type(fem_matrix)     , intent(in)                :: A
    type(fem_partition)  , intent(in)                :: prt 

    type(fem_matrix)     , intent(inout), optional   :: A_II
    type(fem_matrix)     , intent(inout), optional   :: A_IG
    type(fem_matrix)     , intent(inout), optional   :: A_GI
    type(fem_matrix)     , intent(inout), optional   :: A_GG

    integer :: ipoing, offset

    ! Locals
    logical     :: csr_mat_symm, csr_mat_unsymm

    logical     :: present_a_ii, present_a_ig, & 
         & present_a_gi, present_a_gg

    integer(ip) :: ni_rows, nb_rows, ni_cols, nb_cols


    ! Valid storage layouts for fem_matrix are blk and scal
    assert  ( storage == blk .or. storage == scal )
    assert  ( ndof1 >= 1 .and. ndof2 >= 1 )
    assert  ( ndof1 == ndof2 )

    csr_mat_symm   = (A%type == csr_mat .and. output_symm == symm_true )
    csr_mat_unsymm = (A%type == csr_mat .and. output_symm == symm_false)

    assert ( csr_mat_symm .or. csr_mat_unsymm )
    assert ( prt%ptype == element_based )
    assert ( .not. present(A_GI) .or. csr_mat_unsymm )

    if ( storage == blk ) then
       ni_rows = prt%nmap%ni
       nb_rows = prt%nmap%nb
       ni_cols = ni_rows
       nb_cols = nb_rows
    else if ( storage == scal ) then
       ni_rows = prt%nmap%ni * ndof1
       nb_rows = prt%nmap%nb * ndof1
       ni_cols = prt%nmap%ni * ndof2
       nb_cols = prt%nmap%nb * ndof2
    end if

    ! If any inout fem_matrix is present, we are done !
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

       A_II%nd1      = A%nd1
       A_II%nd2      = A%nd2
       A_II%type     = A%type
       A_II%symm     = output_symm
       A_II%storage  = A%storage


       if (A%storage == blk) then 
          call memalloc ( A%nd1, A%nd2, A_II%gr%ia(A_II%gr%nv+1)-A_II%gr%ia(1), A_II%a,       __FILE__,__LINE__)
       else if ( A%storage == scal ) then
          call memalloc ( 1, 1, A_II%gr%ia(A_II%gr%nv+1)-A_II%gr%ia(1), A_II%a,          __FILE__,__LINE__)       
       end if

    end if

    if ( present_a_ig ) then
       assert ( associated(A_IG%gr) ) 
       A_IG%nd1     = A%nd1
       A_IG%nd2     = A%nd2
       A_IG%type    = A%type
       A_IG%symm    = symm_false
       A_IG%storage = A%storage

       if (A%storage == blk) then 
          call memalloc ( A%nd1, A%nd2, A_IG%gr%ia(A_IG%gr%nv+1)-A_IG%gr%ia(1), A_IG%a,          __FILE__,__LINE__ )
       else if (A%storage == scal) then
          call memalloc ( 1, 1, A_IG%gr%ia(A_IG%gr%nv+1)-A_IG%gr%ia(1), A_IG%a,__FILE__,__LINE__ )
       end if


    end if

    if ( present_a_gi ) then
       assert ( associated(A_GI%gr) ) 
       A_GI%nd1     = A%nd1
       A_GI%nd2     = A%nd2
       A_GI%type    = A%type
       A_GI%symm    = symm_false
       A_GI%storage = A%storage

       if (A%storage == blk) then 
          call memalloc ( A%nd1, A%nd2, A_GI%gr%ia(A_GI%gr%nv+1)-A_GI%gr%ia(1), A_GI%a,                         __FILE__,__LINE__ )
       else if (A%storage == scal) then
          call memalloc ( 1, 1, A_GI%gr%ia(A_GI%gr%nv+1)-A_GI%gr%ia(1), A_GI%a,                       __FILE__,__LINE__ )
       end if
       ! write (*,*)  A_GI%gr%ia(A_GI%gr%nv+1)-A_GI%gr%ia(1) ! DBG:
    end if

    if ( present_a_gg ) then

       assert ( associated(A_GG%gr) ) 
       A_GG%nd1     = A%nd1
       A_GG%nd2     = A%nd2
       A_GG%type    = A%type
       A_GG%symm    = output_symm
       A_GG%storage = A%storage

       if (A%storage == blk) then 
          call memalloc ( A%nd1, A%nd2, A_GG%gr%ia(A_GG%gr%nv+1)-A_GG%gr%ia(1), A_GG%a,                         __FILE__,__LINE__ )
       else if (A%storage == scal) then
          call memalloc ( 1, 1, A_GG%gr%ia(A_GG%gr%nv+1)-A_GG%gr%ia(1), A_GG%a,                       __FILE__,__LINE__ )
       end if

    end if

    ! List values on each row of G_II/G_IG
    if ( present_a_ii .or. present_a_ig ) then
       if ( storage == blk ) then
          
          ! For blk storage, output_symm /= A%symm is not implemented yet.
          ! I guess that we actually do not need to implement it, as the usage
          ! of blk storage is almost deprecated
          assert ( output_symm == A%symm )

          do  ipoing=1, ni_rows 

             if (present_a_ii) then 
                A_II%a( 1:A%nd1, 1:A%nd2, A_II%gr%ia(ipoing):A_II%gr%ia(ipoing+1)-1 ) = &  
                     A%a( 1:A%nd1, 1:A%nd2, A%gr%ia(ipoing):A%gr%ia(ipoing)+(A_II%gr%ia(ipoing+1)-A_II%gr%ia(ipoing))-1 )
             end if

             if (present_a_ig) then
                A_IG%a( 1:A%nd1, 1:A%nd2, A_IG%gr%ia(ipoing):A_IG%gr%ia(ipoing+1)-1 ) = &
                     A%a( 1:A%nd1, 1:A%nd2, A%gr%ia(ipoing+1)-(A_IG%gr%ia(ipoing+1)-A_IG%gr%ia(ipoing)):A%gr%ia(ipoing+1)-1)
             end if

          end do
       else if (storage == scal ) then

          do  ipoing=1, ni_rows 
             if ( output_symm == A%symm ) then
                if (present_a_ii) then 
                   A_II%a( 1, 1, A_II%gr%ia(ipoing):A_II%gr%ia(ipoing+1)-1 ) = &  
                        A%a( 1, 1, A%gr%ia(ipoing):A%gr%ia(ipoing)+(A_II%gr%ia(ipoing+1)-A_II%gr%ia(ipoing))-1 )
                end if
             else if ( output_symm == symm_true .and. A%symm == symm_false ) then
                if (present_a_ii) then
                   offset = (A%gr%ia(ipoing+1)- A%gr%ia(ipoing))-(A_IG%gr%ia(ipoing+1)-A_IG%gr%ia(ipoing))-(A_II%gr%ia(ipoing+1)-A_II%gr%ia(ipoing))
                   A_II%a( 1, 1, A_II%gr%ia(ipoing):A_II%gr%ia(ipoing+1)-1 ) = &  
                        A%a( 1, 1, A%gr%ia(ipoing)+offset:A%gr%ia(ipoing)+offset+(A_II%gr%ia(ipoing+1)-A_II%gr%ia(ipoing))-1 )
                end if
             else if ( output_symm == symm_false .and. A%symm == symm_true ) then
                ! Not implemented yet. Trigger an assertion.
                assert ( 1 == 0 )
             end if
             
             if (present_a_ig) then
                A_IG%a( 1, 1, A_IG%gr%ia(ipoing):A_IG%gr%ia(ipoing+1)-1 ) = &
                     A%a( 1, 1, A%gr%ia(ipoing+1)-(A_IG%gr%ia(ipoing+1)-A_IG%gr%ia(ipoing)):A%gr%ia(ipoing+1)-1)
             end if
          end do
          
       end if
    end if


    ! List values on each row of G_GI/G_GG
    if ( present_a_gi .or. present_a_gg ) then
       if ( storage == blk ) then

          ! For blk storage, output_symm /= A%symm is not implemented yet.
          ! I guess that we actually do not need to implement it, as the usage
          ! of blk storage is almost deprecated
          assert ( output_symm == A%symm )
          
          do ipoing=ni_rows+1, ni_rows + nb_rows
             if ( present_a_gi ) then
                ! write (*,*) A_GI%gr%ia( ipoing-prt%nmap%ni ),  A_GI%gr%ia( ipoing+1-prt%nmap%ni )-1                               ! DBG:
                ! write (*,*) A%gr%ia(ipoing), A%gr%ia(ipoing)+(A_GI%gr%ia(ipoing+1-prt%nmap%ni)-A_GI%gr%ia(ipoing-prt%nmap%ni))-1  ! DBG:
                A_GI%a( 1:A%nd1, 1:A%nd2, A_GI%gr%ia( ipoing-prt%nmap%ni ) : A_GI%gr%ia( ipoing+1-prt%nmap%ni )-1 ) = &  
                     A%a ( 1:A%nd1, 1:A%nd2, A%gr%ia(ipoing):A%gr%ia(ipoing)+(A_GI%gr%ia(ipoing+1-prt%nmap%ni)-A_GI%gr%ia(ipoing-prt%nmap%ni))-1 )
             end if
             if ( present_a_gg ) then
                A_GG%a( 1:A%nd1, 1:A%nd2, A_GG%gr%ia( ipoing-prt%nmap%ni ) : A_GG%gr%ia(ipoing+1-prt%nmap%ni)-1 )   = &  
                     A%a ( 1:A%nd1, 1:A%nd2, A%gr%ia(ipoing+1)-(A_GG%gr%ia(ipoing+1-prt%nmap%ni)-A_GG%gr%ia(ipoing-prt%nmap%ni)):A%gr%ia(ipoing+1)-1 )        
             end if
          end do
       else if (storage == scal ) then


          do ipoing=ni_rows+1, ni_rows + nb_rows
             if ( present_a_gi ) then
                ! write (*,*) A_GI%gr%ia( ipoing-ni_rows ),  A_GI%gr%ia( ipoing+1-ni_rows )-1                               ! DBG:
                ! write (*,*) A%gr%ia(ipoing), A%gr%ia(ipoing)+(A_GI%gr%ia(ipoing+1-ni_rows)-A_GI%gr%ia(ipoing-ni_rows))-1  ! DBG:
                A_GI%a( 1, 1, A_GI%gr%ia( ipoing-ni_rows ) : A_GI%gr%ia( ipoing+1-ni_rows )-1 ) = &  
                     A%a ( 1, 1, A%gr%ia(ipoing):A%gr%ia(ipoing)+(A_GI%gr%ia(ipoing+1-ni_rows)-A_GI%gr%ia(ipoing-ni_rows))-1 )
             end if
             if ( present_a_gg ) then
                A_GG%a( 1, 1, A_GG%gr%ia( ipoing-ni_rows ) : A_GG%gr%ia(ipoing+1-ni_rows)-1 )   = &  
                     A%a ( 1, 1, A%gr%ia(ipoing+1)-(A_GG%gr%ia(ipoing+1-ni_rows)-A_GG%gr%ia(ipoing-ni_rows)):A%gr%ia(ipoing+1)-1 )        
             end if
          end do



       end if
     end if

     end subroutine fem_matrix_split_2x2_partitioning

end module fem_matrix_partition

!!$<     if ( present_a_gi .or. present_a_gg ) then  
!!$<       do ipoing=prt%nmap%ni+1, prt%nmap%ni + prt%nmap%nb
!!$<          if ( present_a_gi ) then
!!$<             ! write (*,*) A_GI%gr%ia( ipoing-prt%nmap%ni ),  A_GI%gr%ia( ipoing+1-prt%nmap%ni )-1                               ! DBG:
!!$<             ! write (*,*) A%gr%ia(ipoing), A%gr%ia(ipoing)+(A_GI%gr%ia(ipoing+1-prt%nmap%ni)-A_GI%gr%ia(ipoing-prt%nmap%ni))-1  ! DBG:
!!$<             A_GI%a( 1:A%nd1, 1:A%nd2, A_GI%gr%ia( ipoing-prt%nmap%ni ) : A_GI%gr%ia( ipoing+1-prt%nmap%ni )-1 ) = &  
!!$<               A%a ( 1:A%nd1, 1:A%nd2, A%gr%ia(ipoing):A%gr%ia(ipoing)+(A_GI%gr%ia(ipoing+1-prt%nmap%ni)-A_GI%gr%ia(ipoing-prt%nmap%ni))-1 )
!!$<          end if
!!$<          if ( present_a_gg ) then
!!$<             A_GG%a( 1:A%nd1, 1:A%nd2, A_GG%gr%ia( ipoing-prt%nmap%ni ) : A_GG%gr%ia(ipoing+1-prt%nmap%ni)-1 )   = &  
!!$<               A%a ( 1:A%nd1, 1:A%nd2, A%gr%ia(ipoing+1)-(A_GG%gr%ia(ipoing+1-prt%nmap%ni)-A_GG%gr%ia(ipoing-prt%nmap%ni)):A%gr%ia(ipoing+1)-1 )        
!!$<          end if
!!$<       end do
!!$<     end if
