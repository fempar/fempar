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
module par_uniform_refinement_names
  ! Serial modules
  use types_names
  use memor_names
  use migratory_element_names
  use triangulation_names
  use hash_table_names
  use fe_space_types_names
  use map_names

  use stdio_names
  use mesh_io_names
  use conditions_io_names

  ! Parallel modules
  use par_context_names
  use par_environment_names
  use par_mesh_names
  use par_triangulation_names
  use par_mesh_to_triangulation_names
  use par_conditions_names
  use par_element_exchange_names

  implicit none
# include "debug.i90"
  private
  
  type, extends(migratory_element_t) :: subelems_GIDs_t
     integer(ip)               :: num_subelems
     integer(igp), allocatable :: subelems_GIDs(:)
   contains
     procedure :: size   => subelems_GIDs_size
     procedure :: pack   => subelems_GIDs_pack
     procedure :: unpack => subelems_GIDs_unpack
  end type subelems_GIDs_t

  public :: par_uniform_refinement

contains
  ! This subroutine takes an input triangulation, and generates
  ! a new (output) geometric mesh (p_mesh) obtained after one step
  ! of uniform refinement. The boundary conditions for the refined
  ! mesh are generated on the same (inout) data structure in place.
  ! In other words, it set-ups the data structures required to perform
  ! a simulation on the uniformly refined meshes. 
  ! Some current restrictions:
  !  (1) Does NOT work with mixed elements (e.g., quadrilaterals and triangles).
  !  (2) Only works with meshes of triangles (2D) and tetrahedra (3D) 
  subroutine par_uniform_refinement ( p_trian, p_mesh, p_cond )
#ifdef MPI_MOD
    use mpi
#endif
    implicit none
#ifdef MPI_H
    include 'mpif.h'
#endif
    ! Dummy arguments
    type(par_triangulation_t), target, intent(inout) :: p_trian
    type(par_mesh_t)                 , intent(out)   :: p_mesh 
    type(par_conditions_t)           , intent(inout) :: p_cond

    ! Local variables
    integer(ip), allocatable  :: old2new_vefs(:)
    integer(ip)               :: i, ivef, ielem, jelem, ivertex, isubelem, jsubelem, idime, offset, num_vertices
    integer(ip)               :: vertex1, vertex2
    integer(ip)               :: num_vertices_per_subelem
    integer(ip)               :: num_subelems
    integer(ip)               :: vef_lid, vef_dimension, vef_pos_in_neighbour
    integer(igp)              :: vef_gid
    integer(ip)               :: num_local_vertices_interface
    integer(igp)              :: num_vertices_i_am_owner, num_global_vertices
    integer(ip)               :: num_subelems_interface
    logical                   :: subelem_visited
    integer(igp), allocatable :: ptr_num_subelems_per_part(:)
    integer(ip)               :: max_mypart, elem_lid, jelem_lid
    integer(ip), allocatable  :: subelem_vertices(:,:)
    integer(ip), allocatable  :: code(:,:)
    real(rp)   , allocatable  :: valu(:,:)
    type(list_t)              :: subelems_around_vertices
    type(reference_element_t), pointer :: reference_element
    type(subelems_GIDs_t), allocatable  :: data(:)
    type(hash_table_igp_ip_t)    :: subelems_visited

    integer(ip) :: ierr, istat, lunio
    
    p_mesh%p_env => p_trian%p_env

    if ( p_trian%p_env%am_i_fine_task() ) then
       ! The triangulation must have only a single element type in triangulation
       assert ( p_trian%f_trian%pos_elem_info%last() == 1 )

       ! Element type must be triangles or tetrahedra
       assert ( p_trian%f_trian%reference_elements(1)%ftype == P_type_id )

       reference_element => p_trian%f_trian%reference_elements(1)
       
       ! call reference_element_write(reference_element)

       call generate_data_subelems ( reference_element, &
                                     num_vertices_per_subelem, &
                                     num_subelems, &
                                     subelem_vertices, &
                                     subelems_around_vertices )

       call memalloc (p_trian%f_trian%num_vefs, old2new_vefs, __FILE__, __LINE__)
       old2new_vefs = -1

       num_vertices = 0
       num_vertices_i_am_owner = 0
       do ivef=1, p_trian%f_trian%num_vefs
          if (p_trian%f_trian%vefs(ivef)%dimension < 2) then
             old2new_vefs(ivef) = num_vertices+1
             num_vertices = num_vertices + 1
             if ( p_trian%vefs(ivef)%interface /= 0 ) then ! Interface vef
                max_mypart = 0 
                do ielem=1, p_trian%f_trian%vefs(ivef)%num_elems_around
                   elem_lid = p_trian%f_trian%vefs(ivef)%elems_around(ielem)
                   if ( p_trian%elems(elem_lid)%mypart > max_mypart ) then
                      max_mypart = p_trian%elems(elem_lid)%mypart
                   end if
                end do
                if ( max_mypart == p_trian%p_env%p_context%iam+1) then
                   num_vertices_i_am_owner = num_vertices_i_am_owner + 1
                end if
             else ! Interior vef
                num_vertices_i_am_owner = num_vertices_i_am_owner + 1
             end if
          end if
       end do

       ! Determine number of global vertices in refined mesh
       call MPI_Allreduce(num_vertices_i_am_owner, num_global_vertices, 1, &
                          mpi_integer8, mpi_sum, p_trian%p_env%p_context%icontxt, ierr)
       check(ierr==0)

       ! Allocate vertices map
       call map_alloc( num_vertices, num_global_vertices, p_mesh%f_mesh_dist%nmap )

       p_mesh%f_mesh%nelem = p_trian%f_trian%num_elems * num_subelems
       p_mesh%f_mesh%nnode = num_vertices_per_subelem
       p_mesh%f_mesh%npoin = num_vertices 
       p_mesh%f_mesh%ndime = p_trian%f_trian%num_dims 

       p_mesh%f_mesh_dist%ipart  = p_trian%p_env%p_context%iam + 1 
       p_mesh%f_mesh_dist%nparts = p_trian%p_env%p_context%np

       p_cond%f_conditions%ncond = p_mesh%f_mesh%npoin

       ! Generate new pnods
       call memalloc ( p_mesh%f_mesh%nelem+1, p_mesh%f_mesh%pnods, __FILE__, __LINE__)
       p_mesh%f_mesh%pnods(1) = 1
       do ielem=2,p_mesh%f_mesh%nelem+1
          p_mesh%f_mesh%pnods(ielem)=p_mesh%f_mesh%pnods(ielem-1)+num_vertices_per_subelem
       end do

       ! Generate new lnods, coord and boundary conditions
       call memalloc ( p_mesh%f_mesh%pnods(p_mesh%f_mesh%nelem+1)-1, & 
                       p_mesh%f_mesh%lnods, __FILE__, __LINE__)

       call memalloc ( p_mesh%f_mesh%ndime, &
                       p_mesh%f_mesh%npoin, &
                       p_mesh%f_mesh%coord, __FILE__, __LINE__)

       call memalloc ( p_cond%f_conditions%ncode, &
                       p_mesh%f_mesh%npoin, &
                       code, __FILE__, __LINE__)

       call memalloc ( p_cond%f_conditions%nvalu, &
                       p_mesh%f_mesh%npoin, &
                       valu, __FILE__, __LINE__)

       offset = 1
       do ielem=1, p_trian%f_trian%num_elems
          do isubelem=1, num_subelems
             do ivertex=1, num_vertices_per_subelem
                vef_lid = p_trian%f_trian%elems(ielem)%vefs(subelem_vertices(ivertex,isubelem))
                vef_dimension = p_trian%f_trian%vefs(vef_lid)%dimension
                p_mesh%f_mesh%lnods(offset) = old2new_vefs(vef_lid)
                if ( p_trian%f_trian%vefs(vef_lid)%border /= 0 ) then
                   code(:,old2new_vefs(vef_lid)) = p_cond%f_conditions%code(:,vef_lid)
                   valu(:,old2new_vefs(vef_lid)) = p_cond%f_conditions%valu(:,vef_lid)
                else
                   code(:,old2new_vefs(vef_lid)) = 0 
                   valu(:,old2new_vefs(vef_lid)) = 0.0_rp 
                end if
                p_mesh%f_mesh_dist%nmap%l2g(old2new_vefs(vef_lid)) = p_trian%elems(ielem)%vefs_GIDs(subelem_vertices(ivertex,isubelem))
                if ( vef_dimension == 0 ) then ! Vef is a vertex
                   p_mesh%f_mesh%coord(:,old2new_vefs(vef_lid)) = & 
                        p_trian%f_trian%elems(ielem)%coordinates(:,subelem_vertices(ivertex,isubelem))
                else if ( vef_dimension == 1 ) then ! Vef is an edge
                   ! Extract local ids (within reference element) of vertices of current edge
                   vertex1 = reference_element%crxob%l(reference_element%crxob%p(subelem_vertices(ivertex,isubelem)))
                   vertex2 = reference_element%crxob%l(reference_element%crxob%p(subelem_vertices(ivertex,isubelem))+1)
                   do idime=1, p_mesh%f_mesh%ndime
                      p_mesh%f_mesh%coord(idime,old2new_vefs(vef_lid))= &
                           (p_trian%f_trian%elems(ielem)%coordinates(idime,vertex1) + &
                           p_trian%f_trian%elems(ielem)%coordinates(idime,vertex2))*0.5_rp
                   end do
                else
                   assert(.false.)
                end if
                assert ( p_mesh%f_mesh%lnods(offset) /= -1 )
                offset = offset + 1
             end do
          end do
       end do

       ! Count interface vertices in refined mesh (nnbou)
       num_local_vertices_interface = 0 
       do ivef=1, p_trian%num_itfc_vefs
          vef_lid = p_trian%lst_itfc_vefs(ivef)
          vef_dimension = p_trian%f_trian%vefs(vef_lid)%dimension
          if ( vef_dimension < 2 ) then
             num_local_vertices_interface = num_local_vertices_interface + 1
          end if
       end do
       p_mesh%f_mesh_dist%nnbou = num_local_vertices_interface

       ! List boundary vertices in refined mesh (lnbou)
       call memalloc ( p_mesh%f_mesh_dist%nnbou, p_mesh%f_mesh_dist%lnbou, __FILE__, __LINE__)
       num_local_vertices_interface = 0 
       do ivef=1, p_trian%num_itfc_vefs
          vef_lid = p_trian%lst_itfc_vefs(ivef)
          vef_dimension = p_trian%f_trian%vefs(vef_lid)%dimension
          if ( vef_dimension < 2 ) then
             num_local_vertices_interface = num_local_vertices_interface + 1
             p_mesh%f_mesh_dist%lnbou(num_local_vertices_interface) = old2new_vefs(vef_lid) 
          end if
       end do

       ! All gather number of refined elements per part
       call memalloc ( p_mesh%f_mesh_dist%nparts+1, ptr_num_subelems_per_part, __FILE__, __LINE__ )
       ptr_num_subelems_per_part(1)=1
       call mpi_allgather(int(p_mesh%f_mesh%nelem,igp), 1, MPI_INTEGER8, &
                          ptr_num_subelems_per_part(2), 1, MPI_INTEGER8, & 
                          p_trian%p_env%p_context%icontxt, ierr)
       check(ierr==0)
       do i=1, p_mesh%f_mesh_dist%nparts
          ptr_num_subelems_per_part(i+1)=ptr_num_subelems_per_part(i)+ptr_num_subelems_per_part(i+1) 
       end do
       
       ! Assign global identifiers to refined elements
       call map_alloc( p_mesh%f_mesh%nelem, & 
                       ptr_num_subelems_per_part(p_mesh%f_mesh_dist%nparts+1)-1, & 
                       p_mesh%f_mesh_dist%emap )

       do ielem = 1, p_mesh%f_mesh%nelem
          p_mesh%f_mesh_dist%emap%l2g(ielem) = int(ptr_num_subelems_per_part(p_mesh%f_mesh_dist%ipart) + ielem -1, igp)
       end do
  
       ! Communicate global identifiers of my refined elements to my neighbours
       allocate ( data(p_trian%f_el_import%nelem + p_trian%f_el_import%nghost), stat=istat )
       check(istat==0)

       ! Count interface subelems
       num_subelems_interface = 0

       do ielem = 1, p_trian%num_itfc_elems
          elem_lid = p_trian%lst_itfc_elems(ielem)
          data(elem_lid)%num_subelems = num_subelems
          call memalloc( data(elem_lid)%num_subelems, data(elem_lid)%subelems_GIDs, __FILE__, __LINE__ )
          offset = (elem_lid-1)*num_subelems
          do isubelem = 1, num_subelems
             data(elem_lid)%subelems_GIDs(isubelem) = p_mesh%f_mesh_dist%emap%l2g(offset+isubelem)
             do ivertex=1, num_vertices_per_subelem
                vef_lid = p_trian%f_trian%elems(elem_lid)%vefs(subelem_vertices(ivertex,isubelem))
                if ( p_trian%vefs(vef_lid)%interface /= 0 ) then
                   num_subelems_interface = num_subelems_interface + 1 
                   exit 
                end if
             end do
          end do
       end do
       p_mesh%f_mesh_dist%nebou = num_subelems_interface

       call memalloc ( p_mesh%f_mesh_dist%nebou, p_mesh%f_mesh_dist%lebou, __FILE__, __LINE__ )
       call memalloc ( p_mesh%f_mesh_dist%nebou+1, p_mesh%f_mesh_dist%pextn, __FILE__, __LINE__ )

       call ghost_elements_exchange( p_trian%p_env%p_context%icontxt, p_trian%f_el_import, data )

       ! List interface subelems and count its neighbours
       num_subelems_interface = 0
       p_mesh%f_mesh_dist%pextn = 0
       do ielem = 1, p_trian%num_itfc_elems
          elem_lid = p_trian%lst_itfc_elems(ielem)
          offset = (elem_lid-1)*num_subelems
          do isubelem = 1, num_subelems
             ! AFM: 20 subelements around each subelement is a reasonable estimation for the average
             !      number of subelements around each subelement?
             call subelems_visited%init(20)
             subelem_visited = .false.
             do ivertex=1, num_vertices_per_subelem
                vef_lid = p_trian%f_trian%elems(elem_lid)%vefs(subelem_vertices(ivertex,isubelem))
                if ( p_trian%vefs(vef_lid)%interface /= 0 ) then
                   if ( .not. subelem_visited ) then
                      num_subelems_interface = num_subelems_interface + 1
                      p_mesh%f_mesh_dist%lebou (num_subelems_interface) =  offset + isubelem
                      subelem_visited = .true.
                   end if
                   ! Traverse (ghost) elements around vef_lid and count with hash_table
                   do jelem=1, p_trian%f_trian%vefs(vef_lid)%num_elems_around
                      jelem_lid = p_trian%f_trian%vefs(vef_lid)%elems_around(jelem)
                      if ( p_mesh%f_mesh_dist%ipart /= p_trian%elems(jelem_lid)%mypart ) then
                         ! Identify position of vef_lid in ghost element
                         vef_gid = p_trian%vefs(vef_lid)%globalID
                         do vef_pos_in_neighbour=1, p_trian%elems(jelem_lid)%num_vefs
                            if ( vef_gid == p_trian%elems(jelem_lid)%vefs_GIDs(vef_pos_in_neighbour) ) exit
                         end do
                         assert ( vef_pos_in_neighbour <= p_trian%elems(jelem_lid)%num_vefs )
                         ! Traverse subelements around vef_lid in ghost element

                         ! Count only new subelements
                         do jsubelem=subelems_around_vertices%p(vef_pos_in_neighbour),&
                                     subelems_around_vertices%p(vef_pos_in_neighbour+1)-1  
                            call subelems_visited%put(key=data(jelem_lid)%subelems_GIDs(subelems_around_vertices%l(jsubelem)), &
                                                      val=1, &
                                                      stat=istat)
                            if ( istat == now_stored ) then
                               ! write(*,*) ielem, isubelem, subelems_around_vertices%l(jsubelem)
                               p_mesh%f_mesh_dist%pextn(num_subelems_interface+1)=&
                                      p_mesh%f_mesh_dist%pextn(num_subelems_interface+1)+1
                            end if
                         end do
                      end if
                   end do
                end if
             end do
             call subelems_visited%free()
          end do
       end do

       p_mesh%f_mesh_dist%pextn(1) = 1
       do ielem = 1, p_mesh%f_mesh_dist%nebou
          p_mesh%f_mesh_dist%pextn(ielem+1) = p_mesh%f_mesh_dist%pextn(ielem) +  p_mesh%f_mesh_dist%pextn(ielem+1)
       end do

       call memalloc ( p_mesh%f_mesh_dist%pextn(p_mesh%f_mesh_dist%nebou+1)-1, &
                       p_mesh%f_mesh_dist%lextn, __FILE__, __LINE__ )

       call memalloc ( p_mesh%f_mesh_dist%pextn(p_mesh%f_mesh_dist%nebou+1)-1, &
                       p_mesh%f_mesh_dist%lextp, __FILE__, __LINE__ )


       ! List neighbours of interface subelems
       num_subelems_interface = 0
       do ielem = 1, p_trian%num_itfc_elems
          elem_lid = p_trian%lst_itfc_elems(ielem)
          offset = (elem_lid-1)*num_subelems
          do isubelem = 1, num_subelems
             ! AFM: 20 subelements around each subelement is a reasonable estimation for the average
             !      number of subelements around each subelement?
             call subelems_visited%init(20)
             subelem_visited = .false.
             do ivertex=1, num_vertices_per_subelem
                vef_lid = p_trian%f_trian%elems(elem_lid)%vefs(subelem_vertices(ivertex,isubelem))
                if ( p_trian%vefs(vef_lid)%interface /= 0 ) then
                   if ( .not. subelem_visited ) then
                      num_subelems_interface = num_subelems_interface + 1
                      p_mesh%f_mesh_dist%lebou (num_subelems_interface) =  offset + isubelem
                      subelem_visited = .true.
                   end if
                   ! Traverse (ghost) elements around vef_lid and count with hash_table
                   do jelem=1, p_trian%f_trian%vefs(vef_lid)%num_elems_around
                      jelem_lid = p_trian%f_trian%vefs(vef_lid)%elems_around(jelem)
                      if ( p_mesh%f_mesh_dist%ipart /= p_trian%elems(jelem_lid)%mypart ) then
                         ! Identify position of vef_lid in ghost element
                         vef_gid = p_trian%vefs(vef_lid)%globalID
                         do vef_pos_in_neighbour=1, p_trian%elems(jelem_lid)%num_vefs
                            if ( vef_gid == p_trian%elems(jelem_lid)%vefs_GIDs(vef_pos_in_neighbour) ) exit
                         end do
                         assert ( vef_pos_in_neighbour <= p_trian%elems(jelem_lid)%num_vefs )
                         ! Traverse subelements around vef_lid in ghost element
                         ! Count only new subelements
                         do jsubelem=subelems_around_vertices%p(vef_pos_in_neighbour),&
                                     subelems_around_vertices%p(vef_pos_in_neighbour+1)-1  
                            call subelems_visited%put(key=data(jelem_lid)%subelems_GIDs(subelems_around_vertices%l(jsubelem)), &
                                                      val=1, &
                                                      stat=istat)
                            if ( istat == now_stored ) then
                               p_mesh%f_mesh_dist%lextn(p_mesh%f_mesh_dist%pextn(num_subelems_interface))=&
                                        data(jelem_lid)%subelems_GIDs(subelems_around_vertices%l(jsubelem))
                               p_mesh%f_mesh_dist%lextp(p_mesh%f_mesh_dist%pextn(num_subelems_interface))=&
                                    p_trian%elems(jelem_lid)%mypart
                               p_mesh%f_mesh_dist%pextn(num_subelems_interface)=&
                                    p_mesh%f_mesh_dist%pextn(num_subelems_interface)+1
                            end if
                         end do
                      end if
                   end do
                end if
             end do
             call subelems_visited%free()
          end do
       end do
       
       do ielem = p_mesh%f_mesh_dist%nebou,2,-1
          p_mesh%f_mesh_dist%pextn(ielem) = p_mesh%f_mesh_dist%pextn(ielem-1)
       end do
       p_mesh%f_mesh_dist%pextn(1) = 1

       lunio = io_open ( 'refined.' // trim(ch(p_mesh%f_mesh_dist%ipart)) // '.msh', 'write' )
       call mesh_write_file ( lunio, p_mesh%f_mesh )
       call io_close(lunio)

       ! lunio = io_open ( 'refined.' // trim(ch(p_mesh%f_mesh_dist%ipart)) // '.prt', 'write' )
       ! call mesh_distribution_write ( lunio, p_mesh%f_mesh_dist )
       ! call io_close(lunio)
       
       lunio = io_open ( 'refined.' // trim(ch(p_mesh%f_mesh_dist%ipart)) // '.cnd', 'write' )
       call conditions_write_file ( lunio, p_cond%f_conditions )
       call io_close(lunio)

       do ielem = 1, p_trian%num_itfc_elems
          elem_lid = p_trian%lst_itfc_elems(ielem)
          call memfree( data(elem_lid)%subelems_GIDs, __FILE__, __LINE__ )
       end do

       do ielem = p_trian%f_el_import%nelem+1, p_trian%f_el_import%nelem + p_trian%f_el_import%nghost
          call memfree( data(ielem)%subelems_GIDs, __FILE__, __LINE__ )
       end do

       deallocate ( data, stat=istat )
       check(istat==0)
       
       call memfree ( p_cond%f_conditions%code , __FILE__, __LINE__ )
       call memfree ( p_cond%f_conditions%valu , __FILE__, __LINE__ )
       call memmovealloc (code, p_cond%f_conditions%code, __FILE__, __LINE__)
       call memmovealloc (valu, p_cond%f_conditions%valu, __FILE__, __LINE__)

       call memfree (ptr_num_subelems_per_part, __FILE__, __LINE__ )
       call memfree (old2new_vefs, __FILE__, __LINE__)
       call memfree (subelems_around_vertices%p, __FILE__, __LINE__)
       call memfree (subelems_around_vertices%l, __FILE__, __LINE__)
       call memfree (subelem_vertices, __FILE__, __LINE__)
    end if

  end subroutine par_uniform_refinement
  
  subroutine generate_data_subelems(reference_elem, &
       num_vertices_per_subelem, &
       num_subelems, &
       subelem_vertices, & 
       subelems_around_vertices)
    implicit none
    ! Dummy arguments
    type(reference_element_t), intent(in)   :: reference_elem
    integer(ip), intent(out)                :: num_vertices_per_subelem
    integer(ip), intent(out)                :: num_subelems
    integer(ip), allocatable, intent(out)   :: subelem_vertices(:,:)
    type(list_t), intent(out)               :: subelems_around_vertices

    ! Locals
    integer(ip) :: num_vertices_per_elem 
    integer(ip) :: ielem, ivertex

    assert ( reference_elem%ftype == P_type_id )

    if(reference_elem%ndime==2) then
       num_vertices_per_elem = 6
       num_vertices_per_subelem = 3
       num_subelems = 4
       call memalloc (num_vertices_per_subelem, num_subelems, subelem_vertices, __FILE__, __LINE__)
       subelem_vertices = reshape((/1,4,5,4,2,6,4,6,5,5,6,3/),&
            (/num_vertices_per_subelem,num_subelems/))
    else if(reference_elem%ndime==3) then
       num_vertices_per_elem = 10
       num_vertices_per_subelem = 4
       num_subelems = 8
       call memalloc (num_vertices_per_subelem, num_subelems, subelem_vertices, __FILE__, __LINE__)
       !subelem_vertices = reshape((/1,5,7,8,5,2,6,9,7,6,3,10,8,9,10,4,5,8,9,7,5,9,6,7,10,7,6,9,10,8,7,9/),&
       !subelem_vertices = reshape((/1,5,6,8,5,2,7,9,6,7,3,10,8,9,10,4,5,8,9,6,5,9,7,6,10,6,7,9,10,8,6,9/),&
       subelem_vertices = reshape((/1,5,6,7,  5,2,8,9,  6,8,3,10,  7,9,10,4,  5,7,9,8,  5,6,9,8,  10,6,8,9,  10,7,6,9/),&
       ! The following are inverted:                                                        x
       !subelem_vertices = reshape((/1,5,6,7,  5,2,8,9,  6,8,3,10,  7,9,10,4,  5,7,9,8,  5,9,6,8,  10,6,8,9,  10,7,6,9/),& ! This is the original
            (/num_vertices_per_subelem,num_subelems/))
    end if

    subelems_around_vertices%n = num_vertices_per_elem
    call memalloc ( subelems_around_vertices%n+1, subelems_around_vertices%p, __FILE__, __LINE__ )
    subelems_around_vertices%p = 0
    do ielem=1, num_subelems
       do ivertex=1, num_vertices_per_subelem
          subelems_around_vertices%p(subelem_vertices(ivertex,ielem)+1) = &
               subelems_around_vertices%p(subelem_vertices(ivertex,ielem)+1) + 1
       end do
    end do

    subelems_around_vertices%p(1) = 1
    do ivertex=1, num_vertices_per_elem
       subelems_around_vertices%p(ivertex+1) = subelems_around_vertices%p(ivertex+1) + &
            subelems_around_vertices%p(ivertex)
    end do

    call memalloc ( subelems_around_vertices%p(num_vertices_per_elem+1)-1, subelems_around_vertices%l, __FILE__, __LINE__ )
    do ielem=1, num_subelems
       do ivertex=1, num_vertices_per_subelem
          subelems_around_vertices%l(subelems_around_vertices%p(subelem_vertices(ivertex,ielem))) = ielem
          subelems_around_vertices%p(subelem_vertices(ivertex,ielem)) = &
               subelems_around_vertices%p(subelem_vertices(ivertex,ielem)) + 1
       end do
    end do

    do ivertex=num_vertices_per_elem,2,-1
       subelems_around_vertices%p(ivertex) = subelems_around_vertices%p(ivertex-1)
    end do
    subelems_around_vertices%p(1) = 1

  end subroutine generate_data_subelems


  subroutine subelems_GIDs_size (my, n)
    implicit none
    class(subelems_GIDs_t), intent(in)  :: my
    integer(ip)    , intent(out) :: n
    ! Locals
    integer(ieep) :: mold(1)
    integer(ip)   :: size_of_ip, size_of_igp

    size_of_ip   = size(transfer(1_ip ,mold))
    size_of_igp  = size(transfer(1_igp,mold))

    n  = size_of_ip + & 
         size_of_igp * my%num_subelems 

  end subroutine subelems_GIDs_size
  
  subroutine subelems_GIDs_pack (my, n, buffer)
    implicit none
    class(subelems_GIDs_t), intent(in)   :: my
    integer(ip)    , intent(in)   :: n
    integer(ieep)  , intent(out)  :: buffer(n)
    
    ! Locals
    integer(ieep) :: mold(1)
    integer(ip)   :: size_of_ip, size_of_igp
    integer(ip)   :: start, end
    
    size_of_ip   = size(transfer(1_ip ,mold))
    size_of_igp  = size(transfer(1_igp,mold))
        
    start = 1
    end   = start + size_of_ip -1
    buffer(start:end) = transfer(my%num_subelems,mold)

    start = end + 1
    end   = start + size(my%subelems_GIDs)*size_of_igp - 1
    buffer(start:end) = transfer(my%subelems_GIDs,mold)
  end subroutine subelems_GIDs_pack
    
  subroutine subelems_GIDs_unpack(my, n, buffer)
    implicit none
    class(subelems_GIDs_t) , intent(inout)  :: my
    integer(ip)     , intent(in)     :: n
    integer(ieep)   , intent(in)     :: buffer(n)
        
    ! Locals
    integer(ieep) :: mold(1)
    integer(ip)   :: size_of_ip, size_of_igp
    integer(ip)   :: start, end
        
    size_of_ip   = size(transfer(1_ip ,mold))
    size_of_igp  = size(transfer(1_igp,mold))

    start = 1
    end   = start + size_of_ip -1
    my%num_subelems = transfer(buffer(start:end), my%num_subelems)

    call memalloc ( my%num_subelems, my%subelems_GIDs, __FILE__, __LINE__ )
    start = end + 1
    end   = start + size(my%subelems_GIDs)*size_of_igp - 1
    my%subelems_GIDs = transfer(buffer(start:end), my%subelems_GIDs)

  end subroutine subelems_GIDs_unpack 
end module par_uniform_refinement_names
