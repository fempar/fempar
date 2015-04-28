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
module fem_graph_distribution
  use types
  use memor
  use fem_graph_names
  use dof_distribution_names
  implicit none
# include "debug.i90"
  private

  public :: fem_graph_split_2x2_partitioning, fem_graph_split_2x2_partitioning_symm

contains
  subroutine fem_graph_split_2x2_partitioning ( output_type, & 
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
    type(fem_graph)        , intent(in) :: grph
    type(dof_distribution) , intent(in) :: dof_dist 

    type(fem_graph)     , intent(out), optional   :: G_II
    type(fem_graph)     , intent(out), optional   :: G_IG
    type(fem_graph)     , intent(out), optional   :: G_GI
    type(fem_graph)     , intent(out), optional   :: G_GG


    assert ( grph%type   == csr_symm .or. grph%type   == csr )
    assert ( output_type == csr_symm .or. output_type == csr )
    assert ( .not. present(G_GI) .or. output_type == csr )
    
    ! If any output fem_graph is present, we are done !
    if ( grph%type == csr ) then 
      if ( .not. (present(G_II) .or. present(G_IG) & 
         & .or. present(G_GI) .or. present(G_GG) ) ) return
    else 
        if ( grph%type == csr_symm ) then
          if ( .not. (present(G_II) .or. present(G_IG) .or. present(G_GG) ) ) return
        end if
    end if

    call split_2x2_partitioning_count_list ( output_type, grph, dof_dist, G_II=G_II, G_IG=G_IG, G_GI=G_GI, G_GG=G_GG  )

  end subroutine fem_graph_split_2x2_partitioning

  ! Auxiliary routine. TO-DO: Unpack member allocatable arrays of graph
  ! into explicit size arrays and pass them to auxiliary routines for 
  ! performance reasons
  subroutine split_2x2_partitioning_count_list ( output_type, grph, dof_dist, G_II, G_IG, G_GI, G_GG  )
    implicit none
    ! Parameters
    integer(ip)         , intent(in)              :: output_type
    type(fem_graph)     , intent(in)              :: grph
    type(dof_distribution) , intent(in)              :: dof_dist 

    type(fem_graph)     , intent(out), optional   :: G_II
    type(fem_graph)     , intent(out), optional   :: G_IG
    type(fem_graph)     , intent(out), optional   :: G_GI
    type(fem_graph)     , intent(out), optional   :: G_GG

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

  end subroutine split_2x2_partitioning_count_list
  
  subroutine fem_graph_split_2x2_partitioning_symm ( output_type, & 
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
    type(fem_graph)     , intent(in)              :: grph
    type(dof_distribution) , intent(in)              :: dof_dist 

    type(fem_graph)     , intent(out) :: G_II
    type(fem_graph)     , intent(out) :: G_IG
    type(fem_graph)     , intent(out) :: G_GG

    ! assert ( grph%type   == csr_symm )
    assert ( output_type == csr_symm .or. output_type == csr )
    
    call split_2x2_partitioning_count_list_symm ( output_type, grph, dof_dist, G_II=G_II, G_IG=G_IG, G_GG=G_GG  )

  end subroutine fem_graph_split_2x2_partitioning_symm

  ! Auxiliary routine. TO-DO: Unpack member allocatable arrays of graph
  ! into explicit size arrays and pass them to auxiliary routines for 
  ! performance reasons
  subroutine split_2x2_partitioning_count_list_symm ( output_type, grph, dof_dist, G_II, G_IG, G_GG  )
    implicit none
    ! Parameters
    integer(ip)            , intent(in) :: output_type
    type(fem_graph)        , intent(in) :: grph
    type(dof_distribution) , intent(in) :: dof_dist 

    type(fem_graph)     , intent(out)   :: G_II
    type(fem_graph)     , intent(out)   :: G_IG
    type(fem_graph)     , intent(out)   :: G_GG

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

    ! call fem_graph_print ( 6, grph )


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

  end subroutine split_2x2_partitioning_count_list_symm

end module fem_graph_distribution
