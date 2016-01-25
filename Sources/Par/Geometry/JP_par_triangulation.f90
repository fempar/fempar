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
module JP_par_triangulation_names
  ! Serial modules
  use types_names
  use memor_names
  use sort_names
  use migratory_element_names
  use JP_element_topology_names
  use par_element_topology_names
  use JP_triangulation_names
  use serial_triangulation_names
  use element_import_names
  use hash_table_names

  ! Parallel modules
  use par_context_names
  use par_environment_names

  implicit none
# include "debug.i90"
  private

  type, extends(vef_topology_t) :: par_vef_topology_t
     integer(ip)  :: interface  = -1 ! Interface local id of this vef
     integer(igp) :: globalID   = -1 ! Global ID of this vef
                                     ! Local ID is the position in the array of vefs
   contains
     procedure :: create  => par_vef_topology_create
     procedure :: free    => par_vef_topology_free
  end type par_vef_topology_t

  type, extends(serial_triangulation_t) :: JP_par_triangulation_t
     integer(ip)                             :: num_local   = -1    ! number of local elements
     integer(ip)                             :: num_ghost   = -1    ! number of ghost elements (remote neighbors)
     class(migratory_element_iterator_t), allocatable :: local_elements_iterator
     class(migratory_element_iterator_t), allocatable :: ghost_elements_iterator
     class(migratory_element_iterator_t), allocatable :: boundary_elements_iterator
     integer(ip)                             :: num_itfc_vefs = -1  ! Number of vefs in the interface among subdomains
     integer(ip), allocatable                :: lst_itfc_vefs(:)    ! List of vefs local IDs in the interface among subdomains 
     integer(ip)                             :: num_itfc_elems= -1  ! Number of elements in the interface
     integer(ip), allocatable                :: lst_itfc_elems(:)   ! List of elements local IDs in the interface
     type(par_environment_t), pointer        :: p_env => NULL()     ! Parallel environment describing MPI tasks among which JP_par_triangulation is distributed
     type(element_import_t)                  :: f_el_import         ! Vef describing the layout in distributed-memory of the dual graph
                                                                    ! (It is required for nearest neighbour comms on this graph)
     integer(ip)                             :: max_nparts          ! Maximum number of parts around any vef communication vef
     integer(ip)                             :: nobjs               ! Number of local vef communication vefs
     integer(ip), allocatable                :: lobjs(:,:)          ! List of local vef communication vefs
   contains
     !procedure :: create  => JP_par_triangulation_create
     procedure :: to_dual => JP_par_triangulation_to_dual
     procedure :: free    => JP_par_triangulation_free
     procedure :: create_local_elements_iterator
     procedure :: create_ghost_elements_iterator
     !procedure :: create_boundary_elements_iterator
     !procedure :: print   => JP_par_triangulation_print
     !procedure :: create_element_iterator
     !procedure :: free_element_iterator
  end type JP_par_triangulation_t

  ! Types
  public :: JP_par_triangulation_t, par_vef_topology_t


contains

  !=============================================================================
  !
  ! The following functions should overwritten in par_adaptive_triangulation...
  !
  !=============================================================================
  ! Iterator over local elements
  subroutine create_local_elements_iterator(this)
    implicit none
    class(JP_par_triangulation_t), intent(inout)  :: this
    call this%element_set%create_iterator(this%local_elements_iterator)
    select type(iterator => this%local_elements_iterator)
    class is(plain_migratory_element_iterator_t)
       call iterator%set_last(this%num_local)
    end select
  end subroutine create_local_elements_iterator
  !=============================================================================
  ! Iterator over ghost elements
  subroutine create_ghost_elements_iterator(this)
    implicit none
    class(JP_par_triangulation_t), intent(inout)  :: this
    call this%element_set%create_iterator(this%ghost_elements_iterator)
    select type(iterator => this%ghost_elements_iterator)
    class is(plain_migratory_element_iterator_t)
       call iterator%set_first(this%num_ghost+1)
    end select
  end subroutine create_ghost_elements_iterator
  !=============================================================================
  ! Iterator over boundary elements
  ! subroutine create_boundary_elements_iterator(this)
  !   implicit none
  !   class(JP_par_triangulation_t), intent(inout)  :: this
  !   call this%element_set%create_iterator(this%boundary_elements_iterator)
  !   select type(iterator => this%boundary_elements_iterator)
  !   class is(plain_migratory_element_subset_iterator_t)
  !      call this%element_set%create_subset_iterator(iterator,this%lst_itfc_elems)
  !   end select
  ! end subroutine create_boundary_elements_iterator
  !=============================================================================
  !
  ! The previous functions should overwritten in par_adaptive_triangulation...
  !
  !=============================================================================
  subroutine JP_par_triangulation_create(trian, size)
    implicit none
    class(JP_par_triangulation_t), intent(inout) :: trian

    integer(ip)                  , intent(in)    :: size
    ! Concrete types to select element and element_id in the set
    type(par_element_topology_t) :: element_mold

    ! Allocate and create element_set
    allocate(plain_migratory_element_set_t :: trian%element_set)
    call trian%element_set%create(size,element_mold)

    ! Mother class function (not type bounded by standard restriction)
    !call JP_triangulation_create(trian)
    call trian%JP_triangulation_t%create(size)

  end subroutine JP_par_triangulation_create

  !=============================================================================
  subroutine JP_par_triangulation_free(trian)
    implicit none
    class(JP_par_triangulation_t), intent(inout) :: trian

    if(trian%p_env%p_context%iam<0) return

    ! Don't change the state, that will occur in the mother class
    if ( trian%state == JP_triangulation_vefs_filled ) then
      
       trian%num_itfc_vefs = -1
       call memfree ( trian%lst_itfc_vefs, __FILE__, __LINE__ )

       trian%max_nparts = -1 
       trian%nobjs      = -1 
       call memfree ( trian%lobjs, __FILE__, __LINE__ )

    end if

    ! Don't change the state, that will occur in the mother class
    if ( trian%state == JP_triangulation_elements_filled ) then
       
       trian%num_itfc_elems = -1
       call memfree ( trian%lst_itfc_elems, __FILE__, __LINE__ )
       
       call element_import_free ( trian%f_el_import )

    end if

    ! Mother class function (not type bounded by standard restriction)
    !call JP_triangulation_free(trian)
    call trian%JP_triangulation_t%free()

    ! Deallocate the element structure array and nullify environment
    call trian%element_set%free()
    nullify(trian%p_env)

    trian%state = JP_triangulation_not_created

  end subroutine JP_par_triangulation_free

  !=============================================================================
  subroutine JP_par_triangulation_to_dual(trian)  
    implicit none
    ! Parameters
    class(JP_par_triangulation_t), intent(inout) :: trian

    ! Locals
    integer(ip)              :: iobj, jobj, istat
    class(par_element_topology_t), pointer :: elem

    assert( trian%state == JP_triangulation_elements_filled )

    ! Allocate the vef structure array 
    allocate(par_vef_topology_t :: trian%vefs(trian%num_vefs), stat=istat)
    check(istat==0)
    do iobj=1, trian%num_vefs
       !call initialize_vef_topology(trian%vefs(iobj))
       call trian%vefs(iobj)%create()
    end do

    ! Mother class function
    call trian%serial_triangulation_t%JP_triangulation_t%to_dual()

    select type( vefs => trian%vefs)
    class is(par_vef_topology_t)

       ! Assign global ID for all local vefs looping on local elements
       call trian%create_local_elements_iterator()
       associate( element_iterator =>  trian%local_elements_iterator)

         call element_iterator%begin()
         do while(.not.element_iterator%finished())
            elem => downcast_to_par_element_topology( element_iterator%current() )
            do iobj = 1, elem%num_vefs
               jobj = elem%vefs(iobj)
               if ( jobj /= -1 ) then
                  ! write(*,*) 'ZZZ', elem%vefs_GIDs(iobj)
                  vefs(jobj)%globalID = elem%vefs_GIDs(iobj)
               end if
            end do
            call element_iterator%next()
         end do

       end associate

       ! Count and list interface vefs looping over ghost elements 
       ! (all the interface vefs are there!!!)
       call trian%create_ghost_elements_iterator()
       associate( element_iterator =>  trian%ghost_elements_iterator)
         
         trian%num_itfc_vefs = 0
         call element_iterator%begin()
         do while(.not.element_iterator%finished())
            elem => downcast_to_par_element_topology( element_iterator%current() )
            do iobj = 1, elem%num_vefs
               jobj = elem%vefs(iobj)
               if ( jobj /= -1 ) then
                  if (vefs(jobj)%interface == -1) then
                     trian%num_itfc_vefs = trian%num_itfc_vefs + 1
                     vefs(jobj)%interface = trian%num_itfc_vefs
                  end if
               end if
            end do
            call element_iterator%next()
         end do
         
         call memalloc ( trian%num_itfc_vefs, trian%lst_itfc_vefs, __FILE__, __LINE__ )
         trian%num_itfc_vefs = 0 
         call element_iterator%begin()
         do while(.not.element_iterator%finished())
            elem => downcast_to_par_element_topology( element_iterator%current() )
            do iobj=1, elem%num_vefs
               jobj = elem%vefs(iobj)
               if ( jobj /= -1 ) then
                  if (vefs(jobj)%interface == (trian%num_itfc_vefs + 1)) then
                     trian%num_itfc_vefs = trian%num_itfc_vefs + 1
                     trian%lst_itfc_vefs(trian%num_itfc_vefs) = jobj 
                  end if
               end if
            end do
            call element_iterator%next()
         end do
         
       end associate
       
    class default
    end select

    ! Compute trian%max_nparts, trian%nobjs, trian%lobjs
    ! Re-order (permute) trian%lst_itfc_vefs accordingly
    call JP_par_triangulation_create_lobjs(trian)
      
  end subroutine JP_par_triangulation_to_dual

  !=============================================================================
  subroutine JP_par_triangulation_create_lobjs(trian)
    implicit none
    ! Parameters
    class(JP_par_triangulation_t), intent(inout) :: trian

    ! Locals
    integer(ip)               :: i, j, k, iobj, jelem, est_max_nparts
    integer(ip)               :: ipart, istat, count , touch
    integer(igp), allocatable :: lst_parts_per_itfc_obj (:,:)
    integer(igp), allocatable :: ws_lobjs_temp (:,:)
    integer(igp), allocatable :: sort_parts_per_itfc_obj_l1 (:)
    integer(igp), allocatable :: sort_parts_per_itfc_obj_l2 (:)
    type(hash_table_ip_ip_t)    :: ws_parts_visited
    integer(ip), parameter    :: tbl_length = 100

    class(par_element_topology_t), pointer :: elem

    ipart = trian%p_env%p_context%iam + 1

    ! Compute an estimation (upper bound) of the maximum number of parts around any local interface vef.
    ! This estimation assumes that all elements around all local interface vefs are associated to different parts.
    est_max_nparts = 0
    do i=1, trian%num_itfc_vefs
       iobj = trian%lst_itfc_vefs(i)
       est_max_nparts = max(trian%vefs(iobj)%num_elems_around, est_max_nparts)
    end do

    call memalloc ( est_max_nparts+2, trian%num_itfc_vefs, lst_parts_per_itfc_obj, __FILE__,__LINE__ )

    select type( vefs => trian%vefs)
    class is(par_vef_topology_t)

       touch = 1
       trian%max_nparts = 0
       lst_parts_per_itfc_obj = 0
       do i=1, trian%num_itfc_vefs
          call ws_parts_visited%init(tbl_length)
          !call ws_parts_visited%put(key=ipart,val=1,stat=istat)
          call ws_parts_visited%put(key=ipart,val=touch,stat=istat)
          
          lst_parts_per_itfc_obj (2,i) = ipart
          count = 1

          iobj = trian%lst_itfc_vefs(i)

          ! Count/list parts around iobj 
          do j=1, vefs(iobj)%num_elems_around 
             !call ws_parts_visited%put(key=trian%elems(jelem)%mypart,val=1,stat=istat)
             elem => downcast_to_par_element_topology( trian%element_iterator%get( vefs(iobj)%elems_around(j) ) )
             call ws_parts_visited%put(key=elem%mypart,val=touch,stat=istat)
             if ( istat == now_stored ) then
                count = count + 1
                lst_parts_per_itfc_obj (count+1,i) = elem%mypart 
             end if
          end do

          ! Finish a new column by setting up first and last entries
          lst_parts_per_itfc_obj(1,i) = count
          lst_parts_per_itfc_obj(est_max_nparts+2,i) = vefs(iobj)%globalID 
          
          ! Sort list of parts in increasing order by part identifiers
          ! This is required by the call to icomp subroutine below 
          call sort ( count, lst_parts_per_itfc_obj( 2:(count+1), i) )
          
          trian%max_nparts = max(trian%max_nparts, count)
          call ws_parts_visited%free
       end do


       ! Re-number vefs in increasing order by the number of parts that share them, 
       ! and among vefs sharing the same list of parts, in increasing order by the list 
       ! of parts shared by the vef 
       call memalloc ( est_max_nparts+2, sort_parts_per_itfc_obj_l1, __FILE__,__LINE__ )
       call memalloc ( est_max_nparts+2, sort_parts_per_itfc_obj_l2, __FILE__,__LINE__ )
       call sort_array_cols_by_row_section( est_max_nparts+2,           & ! Rows of lst_parts_per_itfc_obj 
            &                                 est_max_nparts+2,           & ! LD of lst_parts_per_itfc_obj
            &                                 trian%num_itfc_vefs,      & ! Cols of lst_parts_per_itfc_obj
            &                                 lst_parts_per_itfc_obj,     & 
            &                                 trian%lst_itfc_vefs,      &
            &                                 sort_parts_per_itfc_obj_l1, &
            &                                 sort_parts_per_itfc_obj_l2 )
       call memfree ( sort_parts_per_itfc_obj_l2, __FILE__,__LINE__ )
       call memfree ( sort_parts_per_itfc_obj_l1, __FILE__,__LINE__ )

       ! Re-compute vefs(:)%interface to reflect the new status of trian%lst_itfc_vefs
       do i=1, trian%num_itfc_vefs
          iobj = trian%lst_itfc_vefs(i)
          assert ( vefs(iobj)%interface /= -1 )
          vefs(iobj)%interface = i
       end do

       ! Identify communication vefs 
       call memalloc ( trian%max_nparts+4, trian%num_itfc_vefs, ws_lobjs_temp, __FILE__,__LINE__ )
       
       trian%nobjs=1

       ! Prepare first vef
       ws_lobjs_temp(1, trian%nobjs)= -1  ! Unused entry (maintained for historical reasons) 
       ws_lobjs_temp(2, trian%nobjs) = 1  ! Begin obj

       k  = 1 
       do i=1,trian%num_itfc_vefs-1
          if(icomp(trian%max_nparts+1,lst_parts_per_itfc_obj(:,i),lst_parts_per_itfc_obj(:,i+1)) /= 0) then
             ! Complete current vef
             ws_lobjs_temp(3,trian%nobjs)=k ! end obj
             ws_lobjs_temp(4,trian%nobjs)=lst_parts_per_itfc_obj(1,i)
             ws_lobjs_temp(5:4+lst_parts_per_itfc_obj(1,i),trian%nobjs) = &
                  & lst_parts_per_itfc_obj(2:1+lst_parts_per_itfc_obj(1,i),i)
             ws_lobjs_temp(4+lst_parts_per_itfc_obj(1,i)+1:,trian%nobjs) = 0 
             
             ! Prepare next vef
             trian%nobjs=trian%nobjs+1
             ws_lobjs_temp(1,trian%nobjs)=-1     ! Unused entry (maintained for historical reasons) 
             ws_lobjs_temp(2,trian%nobjs)=k+1    ! begin obj
          end if
          k=k+1
       end do

       ! Complete last vef
       ws_lobjs_temp(3,trian%nobjs)=k        ! end obj
       ws_lobjs_temp(4,trian%nobjs)=lst_parts_per_itfc_obj(1,i)
       ws_lobjs_temp(5:4+lst_parts_per_itfc_obj(1,i),trian%nobjs) = &
            & lst_parts_per_itfc_obj(2:1+lst_parts_per_itfc_obj(1,i),i)
       ws_lobjs_temp(4+lst_parts_per_itfc_obj(1,i)+1:,trian%nobjs) = 0 
       
       ! Reallocate lobjs and add internal vef first
       call memalloc (trian%max_nparts+4, trian%nobjs, trian%lobjs, __FILE__, __LINE__)

       ! Copy ws_lobjs_temp to lobjs ...
       do i=1,trian%nobjs
          trian%lobjs(:,i)=ws_lobjs_temp(:,i)
       end do

       call memfree ( ws_lobjs_temp, __FILE__,__LINE__ )
       call memfree ( lst_parts_per_itfc_obj, __FILE__,__LINE__ )


       !write(*,'(a,i10)') 'Number of local vefs:', &
       !   &  trian%nobjs

       !write(*,'(a)') 'List of local vefs:'
       !do i=1,trian%nobjs 
       !   write(*,'(10i10)') i, trian%lobjs(:,i)
       !end do

       
    class default
    end select

  end subroutine JP_par_triangulation_create_lobjs

  !=============================================================================
  ! vef TBP
  subroutine par_vef_topology_create (vef)
    implicit none
    class(par_vef_topology_t), intent(inout) :: vef
    vef%interface   = 0
    vef%globalID = 0
    call vef%vef_topology_t%create()
  end subroutine par_vef_topology_create
  subroutine par_vef_topology_free (vef)
    implicit none
    class(par_vef_topology_t), intent(inout) :: vef
    vef%interface   = -1
    vef%globalID = -1 
    call vef%vef_topology_t%free()
  end subroutine par_vef_topology_free

end module JP_par_triangulation_names
