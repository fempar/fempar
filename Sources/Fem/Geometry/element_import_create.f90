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
module fem_element_import_create_names
use types_names
use memor_names
  use fem_mesh_distribution_names
  use fem_element_import_names
  use sort_names
  use hash_table_names

  implicit none
# include "debug.i90"
  private

  ! Functions
  public :: fem_element_import_create

contains

  subroutine fem_element_import_create (f_msh_dist, f_element_import)
    implicit none

    ! Parameters
    type(fem_mesh_distribution), intent(in)  :: f_msh_dist
    type(fem_element_import)   , intent(out) :: f_element_import
    
    ! Locals
    f_element_import%nparts = f_msh_dist%nparts
    f_element_import%nelem  = f_msh_dist%emap%nl
   
    ! Count neighbours
    call count_neighbours ( f_msh_dist%nebou, &
                            f_msh_dist%pextn, &
                            f_msh_dist%lextp, &
                            f_element_import%npadj) 

    call memalloc ( f_element_import%npadj, f_element_import%lpadj, __FILE__, __LINE__)

    ! List neighbours
    call list_neighbours ( f_msh_dist%nebou, &
                           f_msh_dist%pextn, &
                           f_msh_dist%lextp, &
                           f_element_import%npadj, &
                           f_element_import%lpadj)

    call memalloc ( f_element_import%npadj+1, f_element_import%rcv_ptrs, __FILE__, __LINE__)
    call memalloc ( f_element_import%npadj+1, f_element_import%snd_ptrs, __FILE__, __LINE__)
    
    call count_elements_snd_rcv (f_msh_dist%nebou, &
                                 f_msh_dist%lebou, &
                                 f_msh_dist%pextn, &
                                 f_msh_dist%lextp, &
                                 f_msh_dist%lextn, &
                                 f_element_import%npadj, &
                                 f_element_import%lpadj, &
                                 f_element_import%snd_ptrs, &
                                 f_element_import%rcv_ptrs)

    f_element_import%nghost = f_element_import%rcv_ptrs(f_element_import%npadj+1)-1

    call memalloc ( f_element_import%rcv_ptrs(f_element_import%npadj+1)-1, &
                    f_element_import%rcv_geids, __FILE__, __LINE__)

    call memalloc ( f_element_import%snd_ptrs(f_element_import%npadj+1)-1, &
                    f_element_import%snd_leids, __FILE__, __LINE__)

    call list_elements_snd_rcv(f_msh_dist%nebou, &
                               f_msh_dist%lebou, &
                               f_msh_dist%emap%nl, &
                               f_msh_dist%emap%l2g,  &
                               f_msh_dist%pextn, &
                               f_msh_dist%lextp, &
                               f_msh_dist%lextn, &
                               f_element_import%npadj, &
                               f_element_import%lpadj, &
                               f_element_import%snd_ptrs, &
                               f_element_import%rcv_ptrs, & 
                               f_element_import%snd_leids, &
                               f_element_import%rcv_geids)
                           
  end subroutine fem_element_import_create

  subroutine count_neighbours ( nebou, pextn, lextp, npadj)
    implicit none
    integer(ip), intent(in)  :: nebou
    integer(ip), intent(in)  :: pextn(nebou+1)
    integer(ip), intent(in)  :: lextp(pextn(nebou+1)-1)
    integer(ip), intent(out) :: npadj

    ! Locals
    type(hash_table_ip_ip)   :: parts_visited
    integer(ip)              :: i, istat, touch
    
    call parts_visited%init(20)

    ! Traverse all external elements. Determine whether the part the element
    ! is mapped to has been already taken into account. If not, increment npadj.
    npadj = 0
    touch = 1
    do i=1, pextn(nebou+1)-1
       !call parts_visited%put(key=lextp(i),val=1,stat=istat)
       call parts_visited%put(key=lextp(i),val=touch,stat=istat)
       if(istat==now_stored) npadj = npadj + 1 
    end do
    call parts_visited%free

  end subroutine count_neighbours

  subroutine list_neighbours ( nebou, pextn, lextp, npadj, lpadj)
    implicit none
    integer(ip), intent(in)  :: nebou
    integer(ip), intent(in)  :: pextn(nebou+1)
    integer(ip), intent(in)  :: lextp(pextn(nebou+1)-1)
    integer(ip), intent(in)  :: npadj 
    integer(ip), intent(out) :: lpadj(npadj)

    ! Locals
    type(hash_table_ip_ip)   :: parts_visited
    integer(ip)              :: i, j, istat, touch
    
    call parts_visited%init(20)

    ! Traverse all external elements. Determine whether the part the element
    ! is mapped to has been already taken into account. If not, list it in lpadj.
    j = 1
    touch = 1
    do i=1, pextn(nebou+1)-1
       !call parts_visited%put(key=lextp(i),val=1,stat=istat)
       call parts_visited%put(key=lextp(i),val=touch,stat=istat)
       if(istat==now_stored) then
          lpadj(j) = lextp(i)
          j = j + 1
       end if
    end do

    call parts_visited%free

  end subroutine list_neighbours

  subroutine count_elements_snd_rcv (nebou, lebou, pextn, lextp, lextn, npadj, lpadj, snd_ptrs, rcv_ptrs)
    implicit none
    integer(ip) , intent(in)  :: nebou
    integer(ip) , intent(in)  :: lebou(nebou)
    integer(ip) , intent(in)  :: pextn(nebou+1)
    integer(ip) , intent(in)  :: lextp(pextn(nebou+1)-1)
    integer(igp), intent(in)  :: lextn(pextn(nebou+1)-1)
    integer(ip) , intent(in)  :: npadj 
    integer(ip) , intent(in)  :: lpadj(npadj)
    integer(ip) , intent(out) :: snd_ptrs(npadj+1)
    integer(ip) , intent(out) :: rcv_ptrs(npadj+1)

    ! Locals
    integer(ip) :: i, j, istat, iedge , touch
    type(hash_table_ip_ip) , allocatable :: snd_lids_per_proc(:) 
    type(hash_table_igp_ip), allocatable :: rcv_gids_per_proc(:)

    snd_ptrs = 0
    rcv_ptrs = 0    

    allocate(snd_lids_per_proc(npadj))
    allocate(rcv_gids_per_proc(npadj))

    do i=1, npadj
       ! Assume that nebou elements will be sent to each processor, 
       ! and take its 10% for the size of the hash table
       call snd_lids_per_proc(i)%init(max( int( real(nebou,rp)*0.1_rp, ip), 5))

       ! Assume that as many elements will be received from each processor
       ! as the number of remote edges communicating elements in this part
       ! and any other remote part, and take its 5%
       call rcv_gids_per_proc(i)%init( max ( int( real(pextn(nebou+1),rp)*0.05_rp,ip), 5) )
    end do

    touch = 1
    do i=1, nebou
       do j=pextn(i), pextn(i+1)-1
          
          iedge = 1
          do while ((lpadj(iedge) /= lextp(j)).and. &
               &      (iedge < npadj) )
             iedge = iedge + 1
          end do

          !call snd_lids_per_proc(iedge)%put(key=lebou(i), val=1, stat=istat)
          call snd_lids_per_proc(iedge)%put(key=lebou(i), val=touch, stat=istat)
          if ( istat == now_stored ) then
             snd_ptrs(iedge+1) =  snd_ptrs(iedge+1) + 1
          end if 
          
          !call rcv_gids_per_proc(iedge)%put(key=lextn(j), val=1, stat=istat) 
          call rcv_gids_per_proc(iedge)%put(key=lextn(j), val=touch, stat=istat) 
          if ( istat == now_stored ) then
             rcv_ptrs(iedge+1) =  rcv_ptrs(iedge+1) + 1
          end if

       end do
    end do

    snd_ptrs(1) = 1
    rcv_ptrs(1) = 1
    do i=1, npadj
       call snd_lids_per_proc(i)%free
       call rcv_gids_per_proc(i)%free

       snd_ptrs(i+1) = snd_ptrs(i) + snd_ptrs(i+1)
       rcv_ptrs(i+1) = rcv_ptrs(i) + rcv_ptrs(i+1)
    end do

  end subroutine count_elements_snd_rcv

  subroutine list_elements_snd_rcv (nebou, lebou, nelem, l2ge, pextn, lextp, &
                                    lextn, npadj, lpadj, snd_ptrs, rcv_ptrs, snd_leids, rcv_geids)
    implicit none
    integer(ip) , intent(in)    :: nebou
    integer(ip) , intent(in)    :: lebou(nebou)
    integer(ip) , intent(in)    :: nelem
    integer(igp), intent(in)    :: l2ge(nelem) 
    integer(ip) , intent(in)    :: pextn(nebou+1)
    integer(ip) , intent(in)    :: lextp(pextn(nebou+1)-1)
    integer(igp), intent(in)    :: lextn(pextn(nebou+1)-1)
    integer(ip) , intent(in)    :: npadj 
    integer(ip) , intent(in)    :: lpadj(npadj)
    integer(ip) , intent(inout) :: snd_ptrs(npadj+1)
    integer(ip) , intent(inout) :: rcv_ptrs(npadj+1)
    integer(ip) , intent(out)   :: snd_leids(snd_ptrs(npadj+1)-1)
    integer(igp), intent(out)   :: rcv_geids(rcv_ptrs(npadj+1)-1)

    ! Locals
    integer(ip) :: i, j, istat, iedge , touch
    type(hash_table_ip_ip) , allocatable :: snd_lids_per_proc(:) 
    type(hash_table_igp_ip), allocatable :: rcv_gids_per_proc(:)
    integer(igp) , allocatable           :: snd_geids(:)    

    allocate(snd_lids_per_proc(npadj))
    allocate(rcv_gids_per_proc(npadj))

    call memalloc ( snd_ptrs(npadj+1)-1, snd_geids, __FILE__, __LINE__)

    do i=1, npadj
       ! Assume that nebou elements will be sent to each processor, 
       ! and take its 10% for the size of the hash table
       call snd_lids_per_proc(i)%init(max(int(real(nebou,rp)*0.1_rp,ip),5))

       ! Assume that as many elements will be received from each processor
       ! as the number of remote edges communicating elements in this part
       ! and any other remote part, and take its 5%
       call rcv_gids_per_proc(i)%init(max(int(real(pextn(nebou+1),rp)*0.05_rp,ip),5))
    end do

    touch = 1
    do i=1, nebou
       do j=pextn(i), pextn(i+1)-1
          
          iedge = 1
          do while ((lpadj(iedge) /= lextp(j)).and. &
               &      (iedge < npadj) )
             iedge = iedge + 1
          end do

          !call snd_lids_per_proc(iedge)%put(key=lebou(i), val=1, stat=istat)
          call snd_lids_per_proc(iedge)%put(key=lebou(i), val=touch, stat=istat)

          if ( istat == now_stored ) then
             ! Add to snd_elids
             snd_leids(snd_ptrs(iedge)) = lebou(i)
             
             ! Add to snd_geids
             snd_geids(snd_ptrs(iedge)) = l2ge(lebou(i))

             snd_ptrs(iedge) = snd_ptrs(iedge) + 1
          end if 
          
          !call rcv_gids_per_proc(iedge)%put(key=lextn(j), val=1, stat=istat) 
          call rcv_gids_per_proc(iedge)%put(key=lextn(j), val=touch, stat=istat) 
          
          if ( istat == now_stored ) then
             ! Add to rcv_geids
             rcv_geids(rcv_ptrs(iedge)) = lextn(j)
             rcv_ptrs(iedge) = rcv_ptrs(iedge) + 1
          end if
          
       end do
    end do

    do i=npadj+1, 2, -1
       call snd_lids_per_proc(i-1)%free
       call rcv_gids_per_proc(i-1)%free
       
       snd_ptrs(i) = snd_ptrs(i-1)
       rcv_ptrs(i) = rcv_ptrs(i-1)
    end do
    snd_ptrs(1) = 1
    rcv_ptrs(1) = 1

    do i=1, npadj
       call sort ( snd_ptrs(i+1)-snd_ptrs(i), & 
                   snd_geids(snd_ptrs(i):snd_ptrs(i+1)-1), &
                   snd_leids(snd_ptrs(i):snd_ptrs(i+1)-1)  )

       call sort ( rcv_ptrs(i+1)-rcv_ptrs(i), & 
                   rcv_geids(rcv_ptrs(i):rcv_ptrs(i+1)-1) )
    end do

    call memfree ( snd_geids, __FILE__, __LINE__)

  end subroutine list_elements_snd_rcv

end module fem_element_import_create_names
