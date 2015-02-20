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
module fem_space_faces
  use types
  use memor
  use hash_table_class
  use fem_mesh_class
  use elmat_class
  use elvec_class
  use new_integration_class
  use fem_mesh_faces
  use fem_space_class
  use fem_space_types
  use fem_partition_class
  use array_class
  use adaptivity

# include "debug.i90"
  implicit none

  private
  integer(ip), parameter :: max_subfaces = 4
  ! Functions
  public :: fem_space_faces_list_create

  interface fill_fem_faces_list
     module procedure fill_fem_faces_list_nostruct, fill_fem_faces_list_struct
  end interface fill_fem_faces_list

contains

 !===================================================================================================
  subroutine fem_space_faces_list_create( fspac, ndime,prob_nunk,prob_list_nunk, fpart )
    implicit none
    ! Parameters
    type(fem_space)      , intent(inout),target  :: fspac
    integer(ip)          , intent(in)            :: ndime
    integer(ip), optional, intent(in)            :: prob_nunk
    integer(ip), optional, intent(in)            :: prob_list_nunk(:)
    type(fem_partition), optional, intent(in)    :: fpart

    ! Local variables
    integer(ip) :: nunkn, fadof, aux, fadof_aux(2),nodfac, iobje
    integer(ip) :: iface, jelem, kelem,  istat, nnode, num_interpolations, iint
    integer(ip) :: pos_famat, pos_aux_famat(2), pos_favec, pos_faint,iunk
    integer(ip) :: v_key, gtype(2), utype(2), g_ord(2), u_ord(2),subface(2)!, nface(2)
    type(fem_fixed_info_pointer) :: gfinf(2),ufinf(2)

    if (fspac%o_mesh%tface == 0) then
       write(*,*) __FILE__,__LINE__,'ERROR tface is 0'
       return
    else
       ! Create the list of faces
       allocate(fspac%lface( fspac%o_mesh%tface + fspac%o_mesh%bou_tface), stat=istat)
       if(istat/=0) then 
          write(6,'(a)') '[Fempar Fatal Error] ***Memory (de)allocation failed.'
          write(6,'(a,i20)') 'Error code: ', istat
          write(6,'(a)') 'Caller routine: fem_space::fem_space_create:: lface'
       end if

       ! Create the list of face integrators
       allocate(fspac%lfaci(max_elmat), stat=istat)
       if(istat/=0) then 
          write(6,'(a)') '[Fempar Fatal Error] ***Memory (de)allocation failed.'
          write(6,'(a,i20)') 'Error code: ', istat
          write(6,'(a)') 'Caller routine: fem_space::fem_space_create::lfaci'
       end if
    end if

    ! Fill the list of faces with the neighbouring elements
    if (associated(fspac%o_mesh%constraint_list)) then
       if (present(fpart)) then
          write(*,*) __FILE__,__LINE__,'The listing of the faces is not implemented'
          stop
       else
          call fill_fem_faces_list_nostruct_adapt( fspac%lface, fspac%lelem, fspac%o_mesh )
       end if
    else
       if (present(fpart)) then
          call fill_fem_faces_list_struct( fspac%lface, fspac%o_mesh, fpart )
       else
          call fill_fem_faces_list_nostruct( fspac%lface, fspac%o_mesh )
       end if
    end if

    ! Loop over interior faces
    do iface = 1, fspac%o_mesh%tface
       num_interpolations = fspac%lelem(fspac%lface(iface)%nei_elem(1))%nint
       call memalloc(num_interpolations,fspac%lface(iface)%integ,__FILE__,__LINE__)
       call memalloc(num_interpolations, fspac%lface(iface)%o2n,__FILE__,__LINE__)
       fadof = 0; fadof_aux = 0
       do jelem =1,2
          kelem = fspac%lface(iface)%nei_elem(jelem)  
          if ( present(prob_nunk) ) then
             nunkn = prob_nunk
          else if ( present(prob_list_nunk) ) then
             nunkn = prob_list_nunk(kelem)
          else
             write(6,*) __FILE__,__LINE__,'ERROR! nunk not provided'
             assert(0==1)
          end if
          do iunk=1,nunkn
             if (fspac%kfl_scond == scond_on) then
                nnode = fspac%lelem(kelem)%f_inf(fspac%lelem(kelem)%iv(iunk))%p%nnode -   &
                     &  fspac%lelem(kelem)%f_inf(fspac%lelem(kelem)%iv(iunk))%p%nodes_obj(ndime+1)
             else
                nnode = fspac%lelem(kelem)%f_inf(fspac%lelem(kelem)%iv(iunk))%p%nnode
             end if
             fadof = fadof + nnode
             fadof_aux(jelem)   =  fadof_aux(jelem)   + nnode
             iobje  = fspac%lelem(fspac%lface(iface)%nei_elem(jelem))%f_inf(fspac%lelem(kelem)%iv(iunk))%p%nobje_dim(ndime) +    &
                  &   fspac%lface(iface)%pos_elem(jelem)-1
             nodfac = fspac%lelem(fspac%lface(iface)%nei_elem(jelem))%f_inf(fspac%lelem(kelem)%iv(iunk))%p%ntxob_i(iobje+1) -    &
                  &   fspac%lelem(fspac%lface(iface)%nei_elem(jelem))%f_inf(fspac%lelem(kelem)%iv(iunk))%p%ntxob_i(iobje)
             fadof_aux(3-jelem) =  fadof_aux(3-jelem) + nodfac
          end do
       end do

       ! Create elemental matrix and vectors
       call fspac%ht_pos_elmat%put(key=fadof,val=fspac%cur_elmat,stat=istat)!
       call fspac%ht_pos_elvec%put(key=fadof,val=fspac%cur_elvec,stat=istat)!
       if ( istat == now_stored ) then
          call elmat_alloc (1,1,fadof,fspac%lelmat(fspac%cur_elmat),scal)!
          pos_famat = fspac%cur_elmat!
          fspac%cur_elmat = fspac%cur_elmat + 1!
          call elvec_alloc ( 1, fadof,fspac%lelvec(fspac%cur_elvec),scal)!
          pos_favec = fspac%cur_elvec!
          fspac%cur_elvec = fspac%cur_elvec + 1!
       else if ( istat == was_stored ) then
          call fspac%ht_pos_elmat%get(key=fadof,val=pos_famat,stat=istat)!
          assert ( istat == key_found )
          call fspac%ht_pos_elvec%get(key=fadof,val=pos_favec,stat=istat)!
          assert ( istat == key_found )
       end if
       fspac%lface(iface)%p_mat        => fspac%lelmat(pos_famat)!
       fspac%lface(iface)%p_vec        => fspac%lelvec(pos_favec)!

       ! Create auxiliar assembling matrix
       do jelem =1, 2
          call fspac%ht_pos_elmat%put(key=fadof_aux(jelem),val=fspac%cur_elmat,stat=istat)
          if ( istat == now_stored ) then
             call elmat_alloc (1,1,fadof_aux(jelem),fspac%lelmat(fspac%cur_elmat),scal)
             pos_aux_famat(jelem) = fspac%cur_elmat
             fspac%cur_elmat = fspac%cur_elmat + 1
          else if ( istat == was_stored ) then
             call fspac%ht_pos_elmat%get(key=fadof_aux(jelem),val=pos_aux_famat(jelem),stat=istat)
             assert ( istat == key_found )
          end if
          fspac%lface(iface)%aux_mat(jelem)%p => fspac%lelmat(pos_aux_famat(jelem))
       end do

       do iint = 1, num_interpolations

          gtype(1) = fspac%lelem(fspac%lface(iface)%nei_elem(1))%p_geo_info%ftype
          gtype(2) = fspac%lelem(fspac%lface(iface)%nei_elem(2))%p_geo_info%ftype
          utype(1) = fspac%lelem(fspac%lface(iface)%nei_elem(1))%f_inf(iint)%p%ftype
          utype(2) = fspac%lelem(fspac%lface(iface)%nei_elem(2))%f_inf(iint)%p%ftype
          g_ord(1) = fspac%lelem(fspac%lface(iface)%nei_elem(1))%p_geo_info%order
          g_ord(2) = fspac%lelem(fspac%lface(iface)%nei_elem(2))%p_geo_info%order
          u_ord(1) = fspac%lelem(fspac%lface(iface)%nei_elem(1))%f_inf(iint)%p%order
          u_ord(2) = fspac%lelem(fspac%lface(iface)%nei_elem(2))%f_inf(iint)%p%order
          subface(1) = fspac%lface(iface)%subface(1)
          subface(2) = fspac%lface(iface)%subface(2)

          ! Compute_key
!!$          aux      = (max_ndime+1)*(max_FE_types+1)*(max_order+1)*(max_subfaces+2)
!!$          v_key    =            ndime+(max_ndime+1)*gtype(1)+(max_ndime+1)*(max_FE_types+1)*g_ord(1)
!!$          v_key    = aux*v_key +ndime+(max_ndime+1)*gtype(2)+(max_ndime+1)*(max_FE_types+1)*g_ord(2)
!!$          v_key    = aux*v_key +ndime+(max_ndime+1)*utype(1)+(max_ndime+1)*(max_FE_types+1)*u_ord(1)&
!!$               &               + (max_ndime+1)*(max_FE_types+1)*(max_order+1)*subface(1)
!!$          v_key    = aux*v_key +ndime+(max_ndime+1)*utype(2)+(max_ndime+1)*(max_FE_types+1)*u_ord(2)&
!!$               &               + (max_ndime+1)*(max_FE_types+1)*(max_order+1)*subface(2) 
          aux   = (max_FE_types+1)*(max_order+1)*(max_subfaces+1)
          v_key =  utype(1) + (max_FE_types+1)*u_ord(1) + (max_FE_types+1)*(max_order+1)*subface(1)
          v_key = aux*v_key + &
               &   utype(2) + (max_FE_types+1)*u_ord(2) + (max_FE_types+1)*(max_order+1)*subface(2)
          aux   = (max_FE_types+1)*(max_order+1)
          v_key = aux*v_key + gtype(1)+(max_FE_types+1)*g_ord(1)
          v_key = aux*v_key + gtype(2)+(max_FE_types+1)*g_ord(2)
          ! Put in hash table
          call fspac%ht_pos_face_integ%put(key=v_key, val=fspac%cur_lfaci, stat = istat)
          if ( istat == now_stored ) then
             do jelem=1, 2
                gfinf(jelem)%p => fspac%lelem(fspac%lface(iface)%nei_elem(jelem))%p_geo_info
                ufinf(jelem)%p => fspac%lelem(fspac%lface(iface)%nei_elem(jelem))%f_inf(iint)%p
             end do
             call integ_create(gfinf, ufinf,ndime, fspac%lfaci(fspac%cur_lfaci),subface)
             pos_faint       = fspac%cur_lfaci
             fspac%cur_lfaci = fspac%cur_lfaci + 1
          else if ( istat == was_stored ) then
             call fspac%ht_pos_face_integ%get(key=v_key,val=pos_faint,stat=istat)
             assert ( istat == key_found )
          end if
          fspac%lface(iface)%integ(iint)%p => fspac%lfaci(pos_faint)

          gtype = fspac%lface(iface)%nei_elem
          call array_create(fspac%lface(iface)%integ(iint)%p%quad%ngaus,fspac%lface(iface)%o2n(iint))
          !call memalloc(fspac%lface(iface)%p_integ%quad%ngaus,fspac%lface(iface)%o2n,__FILE__,__LINE__)
          call permute_nodes_object(fspac%lelem(gtype(2))%f_inf(iint)%p,                            &
               &      fspac%lelem(gtype(1))%f_inf(1)%p,fspac%lface(iface)%o2n(iint)%a,              &
               &      fspac%lelem(gtype(2))%f_inf(iint)%p%nobje_dim(ndime) +                        &
               &      fspac%lface(iface)%pos_elem(2)-1,                                             &
               &      fspac%lelem(gtype(1))%f_inf(iint)%p%nobje_dim(ndime) +                        &
               &      fspac%lface(iface)%pos_elem(1)-1,                                             &
               & fspac%o_mesh%lnods(fspac%o_mesh%pnods(gtype(2)):fspac%o_mesh%pnods(gtype(2)+1)-1), &
               & fspac%o_mesh%lnods(fspac%o_mesh%pnods(gtype(1)):fspac%o_mesh%pnods(gtype(1)+1)-1), &
               & ndime-1,fspac%lface(iface)%integ(iint)%p%ufint_ref(1)%nlocf-1,                     &
               & fspac%lface(iface)%subface(2),  fspac%lface(iface)%subface(1))
       end do
    end do

    gfinf(2)%p => NULL()
    ufinf(2)%p => NULL()
    ! BOUNDARY FACES
    do iface = fspac%o_mesh%tface + 1, fspac%o_mesh%tface + fspac%o_mesh%bou_tface
       jelem = 1
       kelem = fspac%lface(iface)%nei_elem(jelem)
       num_interpolations = fspac%lelem(fspac%lface(iface)%nei_elem(1))%nint
       call memalloc(num_interpolations,fspac%lface(iface)%integ,__FILE__,__LINE__)
       fadof = 0; fadof_aux = 0
       if ( present(prob_nunk) ) then
          nunkn = prob_nunk
       else if ( present(prob_list_nunk) ) then
          nunkn = prob_list_nunk(kelem)
       else
          write(6,*) __FILE__,__LINE__,'ERROR! nunk not provided'
          assert(0==1)
       end if
       do iunk=1,nunkn
          if (fspac%kfl_scond == scond_on) then
             nnode = fspac%lelem(kelem)%f_inf(fspac%lelem(kelem)%iv(iunk))%p%nnode -   &
                  &  fspac%lelem(kelem)%f_inf(fspac%lelem(kelem)%iv(iunk))%p%nodes_obj(ndime+1)
          else
             nnode = fspac%lelem(kelem)%f_inf(fspac%lelem(kelem)%iv(iunk))%p%nnode
          end if
          fadof = fadof + nnode
       end do

       ! Assign elemental matrix and vector pointers
       call fspac%ht_pos_elmat%put(key=fadof,val=fspac%cur_elmat,stat=istat)!
       call fspac%ht_pos_elvec%put(key=fadof,val=fspac%cur_elvec,stat=istat)!
       if ( istat == now_stored ) then
          write(*,*) __FILE__,__LINE__,'A boundary face elmat got istat = now_stored!!'
          call elmat_alloc (1,1,fadof,fspac%lelmat(fspac%cur_elmat),scal)!
          pos_famat = fspac%cur_elmat!
          fspac%cur_elmat = fspac%cur_elmat + 1!
          call elvec_alloc ( 1, fadof,fspac%lelvec(fspac%cur_elvec),scal)!
          pos_favec = fspac%cur_elvec!
          fspac%cur_elvec = fspac%cur_elvec + 1!
       else if ( istat == was_stored ) then
          call fspac%ht_pos_elmat%get(key=fadof,val=pos_famat,stat=istat)!
          assert ( istat == key_found )
          call fspac%ht_pos_elvec%get(key=fadof,val=pos_favec,stat=istat)!
          assert ( istat == key_found )
       end if
       fspac%lface(iface)%p_mat        => fspac%lelmat(pos_famat)!
       fspac%lface(iface)%p_vec        => fspac%lelvec(pos_favec)!

       do iint = 1, num_interpolations
          gtype(1) = fspac%lelem(fspac%lface(iface)%nei_elem(1))%p_geo_info%ftype
          gtype(2) =  NULL_type_id
          utype(1) = fspac%lelem(fspac%lface(iface)%nei_elem(1))%f_inf(iint)%p%ftype
          utype(2) =  NULL_type_id
          g_ord(1) = fspac%lelem(fspac%lface(iface)%nei_elem(1))%p_geo_info%order
          u_ord(1) = fspac%lelem(fspac%lface(iface)%nei_elem(1))%f_inf(iint)%p%order
          subface(1) = fspac%lface(iface)%subface(1)
          subface(2) = 0
          ! Compute_key
!!$          aux      = (max_ndime+1)*(max_FE_types+1)*(max_order+1)*(max_subfaces+2)
!!$          v_key    = ndime+(max_ndime+1)*gtype(1)+(max_ndime+1)*(max_FE_types+1)*g_ord(1)
!!$          v_key    = aux*v_key + ndime+(max_ndime+1)*gtype(2)
!!$          v_key    = aux*v_key + ndime+(max_ndime+1)*utype(1)+(max_ndime+1)*(max_FE_types+1)*u_ord(1)&
!!$               &               + (max_ndime+1)*(max_FE_types+1)*(max_order+1)*subface(1)
!!$          v_key    = aux*v_key + ndime+(max_ndime+1)*utype(2)
          aux   = (max_FE_types+1)*(max_order+1)*(max_subfaces+1)
          v_key =  utype(1) + (max_FE_types+1)*u_ord(1) + (max_FE_types+1)*(max_order+1)*subface(1)
          v_key = aux*v_key + &
               &   utype(2) + (max_FE_types+1)*u_ord(2) + (max_FE_types+1)*(max_order+1)*subface(2)
          aux   = (max_FE_types+1)*(max_order+1)
          v_key = aux*v_key + gtype(1)+(max_FE_types+1)*g_ord(1)
          v_key = aux*v_key + gtype(2)+(max_FE_types+1)*g_ord(2)
          ! Put in hash table
          call fspac%ht_pos_face_integ%put(key=v_key, val=fspac%cur_lfaci, stat = istat)
          if ( istat == now_stored ) then 
             jelem = 1    
             gfinf(jelem)%p => fspac%lelem(fspac%lface(iface)%nei_elem(jelem))%p_geo_info
             ufinf(jelem)%p => fspac%lelem(fspac%lface(iface)%nei_elem(jelem))%f_inf(iint)%p
             call integ_create(gfinf,ufinf,ndime,fspac%lfaci(fspac%cur_lfaci))
             pos_faint       = fspac%cur_lfaci
             fspac%cur_lfaci = fspac%cur_lfaci + 1
          else if ( istat == was_stored ) then
             call fspac%ht_pos_face_integ%get(key=v_key,val=pos_faint,stat=istat)
             assert ( istat == key_found )
          end if
          fspac%lface(iface)%integ(iint)%p => fspac%lfaci(pos_faint)
       end do
    end do
  end subroutine fem_space_faces_list_create


 !====================================================================================
  subroutine fill_fem_faces_list_nostruct( lface, tmesh )

    implicit none
    ! Parameters
    type(fem_mesh), intent(in)    :: tmesh
    type(fem_face), intent(inout) :: lface(tmesh%tface+tmesh%bou_tface)

    ! Local variables
    integer(ip)    :: iface, ifacb, jface, kface, jelem, kelem
    integer(ip)    :: i, nface, mface, gface, nelfa
    integer(ip)    :: faceint(tmesh%nelem,tmesh%nface)
    type(fem_mesh) :: d_msh


    call mesh_to_dual(tmesh,d_msh)
    faceint = 0
    iface = 0
    ifacb = 0
    do jelem = 1, tmesh%nelem
       nface = tmesh%pnods(jelem+1) - tmesh%p_face(jelem)
       do jface = 1, nface
          if (faceint(jelem,jface)==0) then
             gface = tmesh%lnods(tmesh%p_face(jelem)+jface-1)
             nelfa = d_msh%pnods(gface+1) - d_msh%pnods(gface)
             if (nelfa == 1) then
                !BOUNDARY FACE
                ifacb = ifacb +1    
                lface(tmesh%tface + ifacb)%nei_elem(1) = jelem
                lface(tmesh%tface + ifacb)%nei_elem(2) = 0
                lface(tmesh%tface + ifacb)%pos_elem(1) = jface
                lface(tmesh%tface + ifacb)%pos_elem(2) = 0
                lface(tmesh%tface + ifacb)%subface    = 0
                lface(tmesh%tface + ifacb)%refinement_level = 0 
             elseif (nelfa==2) then
                !INTERIOR FACE
                do i = 0,1
                   kelem = d_msh%lnods(d_msh%pnods(gface)+i)
                   if (kelem .ne. jelem) exit
                end do
                assert((i==1) .or. (d_msh%lnods(d_msh%pnods(gface)+1) == jelem))
                mface = tmesh%pnods(kelem+1) - tmesh%p_face(kelem)
                do kface = 1, mface
                   if (tmesh%lnods(tmesh%p_face(kelem)+kface-1) == gface) exit
                end do
                assert(kface .le. mface)
                iface = iface +1           
                faceint(kelem,kface) = iface
                lface(iface)%nei_elem(1) = jelem
                lface(iface)%nei_elem(2) = kelem
                lface(iface)%pos_elem(1) = jface
                lface(iface)%pos_elem(2) = kface  
                lface(iface)%subface     = 0
                lface(iface)%refinement_level = 0  
             else
                write(*,*) __FILE__,__LINE__,'ERROR! More than 2 elements per face'
             end if            
          end if
       end do
    end do
    call fem_mesh_free(d_msh)
    assert(iface == tmesh%tface)
    assert(ifacb == tmesh%bou_tface)

  end subroutine fill_fem_faces_list_nostruct

 !====================================================================================
  subroutine fill_fem_faces_list_struct( lface, tmesh, fpart)
    implicit none
    ! Parameters
    type(fem_mesh)     , intent(in)    :: tmesh
    type(fem_partition), intent(in)    :: fpart
    type(fem_face)     , intent(inout) :: lface(tmesh%tface+tmesh%bou_tface)

    ! Local variables
    integer(ip)    :: iface, ifacbe, ifacbi, jface, kface, jelem, kelem, bface
    integer(ip)    :: i, nface, mface, gface, nelfa
    integer(ip)    :: faceint(tmesh%nelem,tmesh%nface)
    type(fem_mesh) :: d_msh


    call mesh_to_dual(tmesh,d_msh)
    faceint = 0
    iface  = 0
    ifacbe = 0
    ifacbi = 0
    bface = tmesh%bou_tface - tmesh%shaface
    do jelem = 1, tmesh%nelem
       nface = tmesh%pnods(jelem+1) - tmesh%p_face(jelem)
       do jface = 1, nface
          if (faceint(jelem,jface)==0) then
             gface = tmesh%lnods(tmesh%p_face(jelem)+jface-1)
             nelfa = d_msh%pnods(gface+1) - d_msh%pnods(gface)
             if (nelfa == 1) then
                if (gface <= fpart%nmap%ni) then
                   ! EXTERNAL BOUNDARY FACE
                   ifacbe = ifacbe + 1    
                   lface(tmesh%tface + ifacbe)%nei_elem(1) = jelem
                   lface(tmesh%tface + ifacbe)%nei_elem(2) = 0
                   lface(tmesh%tface + ifacbe)%pos_elem(1) = jface
                   lface(tmesh%tface + ifacbe)%pos_elem(2) = 0
                   lface(tmesh%tface + ifacbe)%subface     = 0
                   lface(tmesh%tface + ifacbe)%refinement_level = 0 
                else
                   ! INTERFACE BOUNDARY FACE
                   ifacbi = ifacbi + 1    
                   lface(tmesh%tface + bface + ifacbi)%nei_elem(1) = jelem
                   lface(tmesh%tface + bface + ifacbi)%nei_elem(2) = 0
                   lface(tmesh%tface + bface + ifacbi)%pos_elem(1) = jface
                   lface(tmesh%tface + bface + ifacbi)%pos_elem(2) = 0
                   lface(tmesh%tface + bface + ifacbi)%subface     = 0
                   lface(tmesh%tface + bface + ifacbi)%refinement_level = 0 
                end if
             elseif (nelfa==2) then
                !INTERIOR FACE
                do i = 0,1
                   kelem = d_msh%lnods(d_msh%pnods(gface)+i)
                   if (kelem .ne. jelem) exit
                end do
                assert((i==1) .or. (d_msh%lnods(d_msh%pnods(gface)+1) == jelem))
                mface = tmesh%pnods(kelem+1) - tmesh%p_face(kelem)
                do kface = 1, mface
                   if (tmesh%lnods(tmesh%p_face(kelem)+kface-1) == gface) exit
                end do
                assert(kface .le. mface)
                iface = iface +1           
                faceint(kelem,kface) = iface
                lface(iface)%nei_elem(1) = jelem
                lface(iface)%nei_elem(2) = kelem
                lface(iface)%pos_elem(1) = jface
                lface(iface)%pos_elem(2) = kface   
                lface(iface)%subface     = 0
                lface(iface)%refinement_level = 0 
             else
                write(*,*) __FILE__,__LINE__,'ERROR! More than 2 elements per face'
             end if            
          end if
       end do
    end do
    call fem_mesh_free(d_msh)

    assert(iface == tmesh%tface)
    assert((ifacbi+ifacbe) == tmesh%bou_tface)
  end subroutine fill_fem_faces_list_struct

 !====================================================================================
  subroutine fill_fem_faces_list_nostruct_adapt( lface,lelem, tmesh )
    implicit none
    ! Parameters
    type(fem_mesh)   , intent(in)    :: tmesh
    type(fem_face)   , intent(inout) :: lface(tmesh%tface+tmesh%bou_tface)
    type(fem_element), intent(in)    :: lelem(:)

    ! Local variables
    integer(ip)    :: iface, ifacb, jface, kface, jelem, kelem
    integer(ip)    :: i, nface, mface, gface, nelfa
    integer(ip)    :: faceint(tmesh%nelem,tmesh%nface)
    type(fem_mesh) :: d_msh

    integer(ip)                   :: l_obje, nnodexface, num_constraints, p_node, l_node, g_node
    integer(ip)                   :: first_node, p_elem, neigh_elem, neigh_obje, ndime, clist_id
    logical                       :: is_the_face
    type(fem_fixed_info), pointer :: elem_info,neigh_elem_info

    call mesh_to_dual(tmesh,d_msh)
    ndime = tmesh%ndime
    faceint = 0
    iface = 0
    do jelem = 1, tmesh%nelem
       nface = tmesh%pnods(jelem+1) - tmesh%p_face(jelem)
       do jface = 1, nface
          if (faceint(jelem,jface)==0) then
             gface = tmesh%lnods(tmesh%p_face(jelem)+jface-1)
             nelfa = d_msh%pnods(gface+1) - d_msh%pnods(gface)
             if (nelfa == 1) then
                ! We are only interested in faces that are a subface of another element face
                elem_info => lelem(jelem)%p_geo_info
                l_obje = elem_info%nobje_dim(ndime) + jface - 1
                nnodexface = elem_info%crxob_i(l_obje+1) -  elem_info%crxob_i(l_obje)
                num_constraints = 0
                do p_node = elem_info%crxob_i(l_obje),elem_info%crxob_i(l_obje+1)-1
                   l_node = elem_info%crxob_j(p_node)
                   g_node = tmesh%lnods(tmesh%pnods(jelem) + l_node - 1)
                   clist_id = tmesh%constraint_list%map_node2list(g_node)
                   ! ALERT!! This is dangerous because clist_id>0 may be a BC
                   ! BUT a BC always have num_constraints = 1 and nnodexface>1
                   if (clist_id > 0 ) then
                   !if (is_constrained(tmesh%constraint_list,g_node)) then
                      num_constraints = tmesh%constraint_list%list(clist_id)%num_cing_nodes
                      if (num_constraints == nnodexface) exit
                   end if
                end do
                
                ! If there is a constrained node we must check if there is a parent face
                if (num_constraints == nnodexface) then
                   ! Check the first conrtaining node neighbouring elements and, among them if there
                   ! is one that has one face wich is defined by the set of constraining nodes
                   first_node = tmesh%constraint_list%list(clist_id)%cing_nodes(1)
                   eloop: do p_elem = d_msh%pnods(first_node), d_msh%pnods(first_node + 1) -1
                      neigh_elem = d_msh%lnods(p_elem)
                      neigh_elem_info => lelem(neigh_elem)%p_geo_info
                      kface = 0
                      floop: do neigh_obje = neigh_elem_info%nobje_dim(ndime),                      &
                           &                   neigh_elem_info%nobje_dim(ndime+1)-1
                         kface = kface + 1
                         call check_face(num_constraints,                                           &
                              &          tmesh%constraint_list%list(clist_id)%cing_nodes,           &
                              &          neigh_elem,neigh_obje, neigh_elem_info,tmesh%pnods,        &
                              &          tmesh%lnods, is_the_face)
                         if (is_the_face) exit eloop
                      end do floop
                   end do eloop
                   ! Once we found the corresponding neighbour element and face, we can assign values
                   if (is_the_face) then
                      ! TODO: if we want the first element of the pair to be the one with the 
                      ! smallest id, we should check it here.
                      iface = iface + 1
                      faceint(jelem,jface)      = iface
                      faceint(neigh_elem,kface) = iface
                      lface(iface)%nei_elem(1)  = jelem
                      lface(iface)%nei_elem(2)  = neigh_elem
                      lface(iface)%pos_elem(1)  = jface
                      lface(iface)%pos_elem(2)  = kface
                      lface(iface)%subface(1)   = 0
                      call check_subface (jelem, l_obje, neigh_elem, neigh_obje, elem_info,         &
                           &neigh_elem_info, tmesh%pnods, tmesh%lnods, lface(iface)%subface(2))
                      lface(iface)%refinement_level(1)     = 0 
                      lface(iface)%refinement_level(2)     = 1 
                   end if
                end if
             elseif (nelfa==2) then
                !INTERIOR FACE
                do i = 0,1
                   kelem = d_msh%lnods(d_msh%pnods(gface)+i)
                   if (kelem .ne. jelem) exit
                end do
                assert((i==1) .or. (d_msh%lnods(d_msh%pnods(gface)+1) == jelem))
                mface = tmesh%pnods(kelem+1) - tmesh%p_face(kelem)
                do kface = 1, mface
                   if (tmesh%lnods(tmesh%p_face(kelem)+kface-1) == gface) exit
                end do
                assert(kface .le. mface)
                iface = iface +1     
                faceint(jelem,jface) = iface      
                faceint(kelem,kface) = iface
                lface(iface)%nei_elem(1) = jelem
                lface(iface)%nei_elem(2) = kelem
                lface(iface)%pos_elem(1) = jface
                lface(iface)%pos_elem(2) = kface 
                lface(iface)%subface     = 0
                lface(iface)%refinement_level = 0 
             else
                write(*,*) __FILE__,__LINE__,'ERROR! More than 2 elements per face'
             end if            
          end if
       end do
    end do
    assert(iface == tmesh%tface)
    ifacb = 0
    do jelem = 1, tmesh%nelem
       nface = tmesh%pnods(jelem+1) - tmesh%p_face(jelem)
       do jface = 1, nface
          if (faceint(jelem,jface)==0) then
             gface = tmesh%lnods(tmesh%p_face(jelem)+jface-1)
             nelfa = d_msh%pnods(gface+1) - d_msh%pnods(gface)
             if (nelfa == 1) then
                !BOUNDARY FACE
                ifacb = ifacb +1    
                lface(tmesh%tface + ifacb)%nei_elem(1) = jelem
                lface(tmesh%tface + ifacb)%nei_elem(2) = 0
                lface(tmesh%tface + ifacb)%pos_elem(1) = jface
                lface(tmesh%tface + ifacb)%pos_elem(2) = 0
                lface(tmesh%tface + ifacb)%subface     = 0
                lface(tmesh%tface + ifacb)%refinement_level       = 0 
             elseif (nelfa==2) then
                write(*,*) __FILE__,__LINE__,'ERROR! It should be a boundary face'
             end if
          end if
       end do
    end do
    call fem_mesh_free(d_msh)
    assert(ifacb == tmesh%bou_tface)
  end subroutine fill_fem_faces_list_nostruct_adapt

  !====================================================================================
  subroutine check_face (num_nodes,list_nodes, ielem, iobje, elem_info, pnods, lnods, is_the_face)
    implicit none
    ! Parameters
    integer(ip)         , intent(in)  :: num_nodes, list_nodes(num_nodes)
    integer(ip)         , intent(in)  :: ielem, iobje
    type(fem_fixed_info), intent(in)  :: elem_info
    integer(ip)         , intent(in)  :: pnods(:),lnods(:)
    logical             , intent(out) :: is_the_face

    ! Local variables
    integer(ip) :: p_first_node, p_l_node, l_node, g_node, inode

    p_first_node = pnods(ielem) -1
    
    ! Loop over all the nodes in the iobje of ielem
    ploop: do p_l_node = elem_info%crxob_i(iobje),elem_info%crxob_i(iobje+1)-1
       ! local node id
       l_node = elem_info%crxob_j(p_l_node)
       ! Global node id
       g_node = lnods(p_first_node + l_node)
       ! Loop over all the nodes of the face; check if the node is in the list
       iloop: do inode = 1, num_nodes
          ! If in the list, check next node
          if(list_nodes(inode) == g_node) cycle ploop
       end do iloop
       ! If not in the list, this is not the face we were looking for
       is_the_face = .false.
       return
    end do ploop
    is_the_face = .true.

    
  end subroutine check_face

 !====================================================================================
  ! This subroutines checkes for neigh_elem wich subfaces of neigh_obje corresponds to the obje of 
  ! elem. It is assuming that there ist at most one level of refinement and thus, they have one 
  ! common node.
  subroutine check_subface (elem, obje, neigh_elem, neigh_obje, elem_info, neigh_elem_info, pnods,   &
       &                    lnods, subface)
    implicit none
    ! Parameters
    integer(ip)         , intent(in)  :: elem, obje
    integer(ip)         , intent(in)  :: neigh_elem, neigh_obje
    type(fem_fixed_info), intent(in)  :: elem_info, neigh_elem_info
    integer(ip)         , intent(in)  :: pnods(:),lnods(:)
    integer(ip)         , intent(out) :: subface

    ! Local variables
    integer(ip) :: p_first_node, p_l_node, l_node, g_node
    integer(ip) :: p_neigh_first_node, p_l_neigh_node, l_neigh_node, g_neigh_node

    p_first_node       = pnods(elem) - 1
    p_neigh_first_node = pnods(neigh_elem) - 1

    subface = 1
    ! Loop over all the nodes in the neigh_obje of neigh_elem 
    ploop: do p_l_neigh_node = neigh_elem_info%crxob_i(neigh_obje),                                 &
         &                     neigh_elem_info%crxob_i(neigh_obje+1)-1
       ! local node id
       l_neigh_node = elem_info%crxob_j(p_l_neigh_node)
       ! Global node id
       g_neigh_node = lnods(p_neigh_first_node + l_neigh_node)
       ! Loop over all the nodes of the obje of elem; check if the node is in the list
       do p_l_node = elem_info%crxob_i(obje),elem_info%crxob_i(obje+1)-1
          ! local node id
          l_node = elem_info%crxob_j(p_l_node)
          ! Global node id
          g_node = lnods(p_first_node + l_node)
          ! Loop over all the nodes of the face; check if the node is in the list
          if(g_neigh_node == g_node) exit ploop
       end do
       subface = subface + 1
    end do ploop

    assert(subface <= neigh_elem_info%crxob_i(neigh_obje+1)- neigh_elem_info%crxob_i(neigh_obje))
    
  end subroutine check_subface

end module fem_space_faces
