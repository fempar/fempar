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
# include "debug.i90"
module fem_space_names
  ! Modules
  use types
  use memor
  use array_names
  use fem_triangulation_names
  use hash_table_names
  use problem_names
  use integration_tools_names
  !use face_integration_names
  use fem_space_types
  use dof_handler_names
  use migratory_element_names
  use fem_conditions_names
  use fem_element_names

#ifdef memcheck
  use iso_c_binding
#endif
  implicit none
  private

  integer(ip), parameter       :: max_global_interpolations  = 50   ! Maximum number of interpolations

  ! Global information of the fem_space
  type fem_space  

     integer(ip)                           :: num_continuity       ! Number of materials (maximum value)
     logical(lg)                           :: static_condensation  ! Flag for static condensation 
     logical(lg)                           :: hierarchical_basis   ! Flag for hierarchical basis
     class(migratory_element), allocatable :: mig_elems(:)         ! Migratory elements list
     type(fem_element)       , pointer     :: lelem(:)             ! List of FEs
     type(fem_face)          , allocatable :: lface(:)             ! List of active faces

     type(fem_triangulation)  , pointer :: g_trian => NULL() ! Triangulation
     type(dof_handler)        , pointer :: dof_handler

     ! Formulations used to build the space
     integer(ip) :: num_approximations
     type(discrete_problem_pointer), allocatable :: approximations(:)

     ! Array of working arrays (element matrix/vector) (to be pointed from fem_elements)
     type(position_hash_table)          :: pos_elmatvec
     type(array_rp2)                    :: lelmat(max_global_interpolations)
     type(array_rp1)                    :: lelvec(max_global_interpolations)

     ! Integrator
     type(position_hash_table)          :: pos_volume_integrator
     type(position_hash_table)          :: pos_face_integrator
     type(volume_integrator)            :: lvoli(max_global_interpolations)        
     type(face_integrator)              :: lfaci(max_global_interpolations)

     ! Array of reference elements (to be pointed from fem_elements)
     type(position_hash_table)          :: pos_elem_info
     type(fem_fixed_info)               :: lelem_info(max_elinf)

     ! Acceleration arrays
     type(list_2d), allocatable         :: object2dof(:)       ! An auxiliary array to accelerate some parts of the code
     integer(ip), allocatable           :: ndofs(:)            ! number of dofs (nblocks)
     integer(ip)                        :: time_steps_to_store ! Time steps to store in unkno

     ! List of faces where we want to integrate things
     type(fem_face), allocatable        :: interior_faces(:), boundary_faces(:), interface_faces(:)
     integer(ip) :: num_interior_faces, num_boundary_faces, num_interface_faces

     ! Much better here rather than as a module variable
     type(list) :: void_list
  end type fem_space

  ! Types
  public :: fem_space

  ! Functions
  public :: fem_space_create, fem_space_print, fem_space_free

contains

  !==================================================================================================
  ! Allocation of variables in fem_space according to the values in g_trian
  subroutine fem_space_create( g_trian, dofh, fspac, problem, approximations, bcond, continuity, order, material, &
       &                       which_approx, num_approximations, time_steps_to_store, hierarchical_basis, &
       &                       static_condensation, num_continuity, num_ghosts )
    implicit none
    type(fem_space)        , target, intent(inout) :: fspac
    type(fem_triangulation), target, intent(in)    :: g_trian   
    type(dof_handler)      , target, intent(in)    :: dofh   
    type(discrete_problem_pointer) , intent(in)    :: approximations(:)
    type(fem_conditions)           , intent(in)    :: bcond
    integer(ip), intent(in) :: material(:), order(:,:), problem(:)
    integer(ip), intent(in) :: continuity(:,:), which_approx(:)
    integer(ip), intent(in) :: num_approximations
    integer(ip), optional, intent(in) :: time_steps_to_store
    logical(lg), optional, intent(in) :: hierarchical_basis
    logical(lg), optional, intent(in) :: static_condensation
    integer(ip), optional, intent(in) :: num_continuity   
    integer(ip), optional, intent(in) :: num_ghosts         

    integer(ip) :: istat, num_ghosts_

    call fem_space_allocate_structures(  g_trian, dofh, fspac, num_approximations=num_approximations,&
         time_steps_to_store = time_steps_to_store, hierarchical_basis = hierarchical_basis, &
         static_condensation = static_condensation, num_continuity = num_continuity, &
         num_ghosts = num_ghosts )  

    assert(size(approximations)==num_approximations)
    fspac%approximations = approximations

    call fem_space_fe_list_create ( fspac, problem, which_approx, continuity, order, material, bcond )

  end subroutine fem_space_create

  !==================================================================================================
  ! Allocation of variables in fem_space according to the values in g_trian
  subroutine fem_space_allocate_structures( g_trian, dofh, fspac, num_approximations, &
       & time_steps_to_store, hierarchical_basis, static_condensation, num_continuity, &
       & num_ghosts )
    implicit none
    type(fem_space)   ,  target, intent(inout)                :: fspac
    type(fem_triangulation)   , intent(in), target   :: g_trian   
    type(dof_handler) , intent(in), target           :: dofh  
    integer(ip), optional, intent(in) :: num_approximations
    integer(ip), optional, intent(in) :: time_steps_to_store
    logical(lg), optional, intent(in) :: hierarchical_basis
    logical(lg), optional, intent(in) :: static_condensation
    integer(ip), optional, intent(in) :: num_continuity 
    integer(ip), optional, intent(in) :: num_ghosts            

    integer(ip) :: istat, num_ghosts_

    
    ! Approximations flag
    if (present(num_approximations)) then
       fspac%num_approximations = num_approximations
    else
       fspac%num_approximations = 1
    end if

    ! Hierarchical flag
    if (present(hierarchical_basis)) then
       fspac%hierarchical_basis = hierarchical_basis
    else
       fspac%hierarchical_basis = .false.
    end if

    ! Static condensation flag
    if (present(static_condensation)) then
       fspac%static_condensation = static_condensation
    else
       fspac%static_condensation = .false.
    end if

    ! Number materials flag
    if (present(num_continuity)) then
       fspac%num_continuity = num_continuity
    else
       fspac%num_continuity = 1
    end if

    ! Time steps to store flag
    if (present(time_steps_to_store)) then
       fspac%time_steps_to_store = time_steps_to_store
    else
       fspac%time_steps_to_store = 1
    end if

    ! Ghosts elements (parallel case)
    if (present(num_ghosts)) then
       num_ghosts_ = num_ghosts
    else
       num_ghosts_ = 0
    end if

    allocate( fem_element :: fspac%mig_elems(g_trian%num_elems + num_ghosts_), stat=istat)
    check ( istat == 0 )
    
    select type( this => fspac%mig_elems )
    type is(fem_element)
       fspac%lelem => this
    end select

    !  Initialization of pointer to triangulation
    fspac%g_trian => g_trian

    !  Initialization of pointer to dof_handler
    fspac%dof_handler => dofh

    ! Allocation of discretizations
    allocate(fspac%approximations(num_approximations),stat=istat)
    check(istat==0)

    ! Allocation of elemental matrix and vector parameters
    call fspac%pos_elmatvec%init(ht_length)

    ! Initialization of volume and face integrators parameters
    call fspac%pos_volume_integrator%init(ht_length)
    call fspac%pos_face_integrator%init(ht_length)

    ! Initialization of element fixed info parameters
    call fspac%pos_elem_info%init(ht_length)

  end subroutine fem_space_allocate_structures

  !==================================================================================================
  ! Fill the fem_space assuming that all elements are of type f_type but each variable has different
  ! interpolation order
  subroutine fem_space_fe_list_create( fspac, problem, which_approx, continuity, order, material, bcond )
    implicit none
    type(fem_space), intent(inout), target  :: fspac
    integer(ip)    , intent(in)       :: material(:), order(:,:), problem(:)
    integer(ip)    , intent(in)       :: continuity(:,:), which_approx(:)
    type(fem_conditions), intent(in)  :: bcond

    integer(ip) :: nunk, v_key, ltype(2), nnode, max_num_nodes, nunk_tot, dim, f_order, f_type, nvars, nvars_tot
    integer(ip) :: ielem, istat, pos_elmatvec, pos_elinf, pos_voint, ivar, lndof, iobje
    logical(lg) :: created
    integer(ip) :: aux_val

    ! nunk_tot = total amount of unknowns
    ! if ( present(prob_list_nunk) ) then
    !    nunk_tot=prob_list_nunk(1)
    !    do ielem=2, fspac%g_trian%num_elems
    !       nunk_tot=max(nunk_tot,prob_list_nunk(ielem))
    !    end do
    ! else
    !    nunk_tot = prob_nunk
    ! end if

    ! void nodes object
    fspac%void_list%n = max_nobje
    call memalloc( max_nobje+1, fspac%void_list%p, __FILE__, __LINE__)
    fspac%void_list%p = 1 
    call memalloc( 0, fspac%void_list%l, __FILE__, __LINE__ )

    ! fspac%l_nodes_object(1)%n = max_nobje
    ! call memalloc( max_nobje+1, fspac%l_nodes_object(1)%p, __FILE__, __LINE__ )
    ! fspac%l_nodes_object(1)%p = 1
    ! call memalloc( 0, fspac%l_nodes_object(1)%l, __FILE__, __LINE__ )

    ! Loop over elements
    nvars_tot = fspac%dof_handler%nvars_global
    dim = fspac%g_trian%num_dims

    ! Material
    fspac%lelem(1:fspac%g_trian%num_elems)%material = material

    ! Continuity
    do ielem = 1, fspac%g_trian%num_elems
       ! Assign type of problem and approximation to ielem
       fspac%lelem(ielem)%problem =  problem(ielem)
       fspac%lelem(ielem)%approximation =  which_approx(ielem)

       nvars = fspac%dof_handler%problems(problem(ielem))%p%nvars
       fspac%lelem(ielem)%num_vars = nvars
       f_type = fspac%g_trian%elems(ielem)%topology%ftype
       !assert(size(continuity,1)==nvars_tot)
       !assert(size(material,1)==nvars_tot)

       !SB.alert : Not OK

       ! Set continuity per unknown
       call memalloc(nvars, fspac%lelem(ielem)%continuity, __FILE__, __LINE__)

       ! Set order per unknown
       call memalloc(nvars, fspac%lelem(ielem)%order, __FILE__, __LINE__)

       ! Allocate the fixed info array
       call memalloc(nvars, fspac%lelem(ielem)%f_inf, __FILE__, __LINE__)

       ! Allocate nodes per object info
       allocate(fspac%lelem(ielem)%nodes_object(nvars), stat=istat )
       check( istat==0 )

       ! Assign pointer to interpolation fixed information and nodes per object
       do ivar=1,nvars
          !write(*,*) 'ielem,ivar',ielem,ivar,continuity(2,1)
          !write(*,*) 'POINT REACHED'
          !write(*,*) 'l2g', fspac%dof_handler%problems(problem(ielem))%p%l2g_var(ivar)
          !write(*,*) 'cont',continuity(ielem,fspac%dof_handler%problems(problem(ielem))%p%l2g_var(ivar))

          ! JP: indices of these arrays (continuity and order) should be changed to (nvars,nelem)
          fspac%lelem(ielem)%continuity(ivar) = continuity(ielem,fspac%dof_handler%problems(problem(ielem))%p%l2g_var(ivar))
          fspac%lelem(ielem)%order(ivar) = order(ielem,fspac%dof_handler%problems(problem(ielem))%p%l2g_var(ivar))

          f_order = fspac%lelem(ielem)%order(ivar)
          v_key = dim + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)*f_order
          call fspac%pos_elem_info%get(key=v_key,val=pos_elinf,stat=istat)
          if ( istat == new_index) then 
             call fem_element_fixed_info_create(fspac%lelem_info(pos_elinf),f_type,              &
                  &                             f_order,dim,created)
             assert(created)
          end if
          fspac%lelem(ielem)%f_inf(ivar)%p => fspac%lelem_info(pos_elinf)

          if ( continuity(ielem, ivar) /= 0 ) then
             fspac%lelem(ielem)%nodes_object(ivar)%p => fspac%lelem_info(pos_elinf)%ndxob
          else 
             fspac%lelem(ielem)%nodes_object(ivar)%p => fspac%void_list ! fspac%l_nodes_object(1) ! SB.alert : Think about hdG
          end if
       end do

       ! Assign pointer to geometrical fixed information (assumed to be of order 1)
       v_key = dim + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)
       call fspac%pos_elem_info%get(key=v_key,val=pos_elinf,stat=istat)
       if ( istat == new_index) then 
          call fem_element_fixed_info_create(fspac%lelem_info(pos_elinf),f_type,              &
               &                             1,dim,created)
          assert(created)
       end if
       fspac%lelem(ielem)%p_geo_info => fspac%lelem_info(pos_elinf)

       ! Compute number of DOFs in the elemental matrix associated to ielem
       lndof = 0
       max_num_nodes = 0
       do ivar=1,nvars
          if ( fspac%static_condensation ) then
             ! JP: Is this working? Because max_num_nodes is sent to integ_create to build the interpolation
             !     which is needed even with static_condensation...
             ! JP & SB: In fact we also need to think about it for the unknown. Only elem2dof can be reduced
             !          to the interface.
             nnode = fspac%lelem(ielem)%f_inf(ivar)%p%nnode -                                   &
                  &  fspac%lelem(ielem)%f_inf(ivar)%p%nodes_obj(dim+1) ! SB.alert : do not use nodes_obj
          else
             nnode = fspac%lelem(ielem)%f_inf(ivar)%p%nnode 
          end if
          lndof = lndof + nnode
          max_num_nodes = max(max_num_nodes,nnode)
       end do

       ! Assign pointer to p_mat and p_vec in ielem
       call fspac%pos_elmatvec%get(key=lndof,val=pos_elmatvec,stat=istat)
       if ( istat == new_index ) then
          call array_create ( lndof, lndof, fspac%lelmat(pos_elmatvec) )
          call array_create ( lndof, fspac%lelvec(pos_elmatvec) )
       end if
       fspac%lelem(ielem)%p_mat => fspac%lelmat(pos_elmatvec)
       fspac%lelem(ielem)%p_vec => fspac%lelvec(pos_elmatvec)

       ! Allocate elem2dof, unkno, bc_code
       write(*,*) max_num_nodes, nvars
       call memalloc( max_num_nodes, nvars, fspac%lelem(ielem)%elem2dof, __FILE__,__LINE__ )
       fspac%lelem(ielem)%elem2dof = 0
       call memalloc( max_num_nodes, nvars, fspac%time_steps_to_store, fspac%lelem(ielem)%unkno, __FILE__,__LINE__)
       fspac%lelem(ielem)%unkno = 0.0_rp
       call memalloc(nvars,fspac%lelem(ielem)%integ,__FILE__,__LINE__)
       call memalloc(nvars,fspac%lelem(ielem)%p_geo_info%nobje,fspac%lelem(ielem)%bc_code,__FILE__,__LINE__, 0)

       ! Fill bc_code
       do iobje = 1,fspac%lelem(ielem)%p_geo_info%nobje
          do ivar=1,nvars
             fspac%lelem(ielem)%bc_code(ivar,iobje) = bcond%code( fspac%dof_handler%problems(problem(ielem))%p%l2g_var(ivar),  fspac%g_trian%elems(ielem)%objects(iobje) )
          end do
       end do

       ! Assign pointers to volume integration
       ltype(2) = dim + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)
       do ivar = 1,nvars
          ! Here f_order should be taken from fspac%lelem(ielem)%order(ivar), doesn't it?
          ! f_order = fspac%lelem(ielem)%order(ivar)
          ! Otherwise all the variables have the same order.
          ltype(1) = dim + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)*f_order
          v_key    = (max_ndime+1)*(max_FE_types+1)*(max_order) * ltype(1) + ltype(2)
          call fspac%pos_volume_integrator%get(key=v_key, val=pos_voint, stat = istat)
          ! SB.alert : g_ord = 1 !!!! But only for linear geometry representation
          if ( istat == new_index ) call volume_integrator_create(f_type,f_type,dim,1,f_order,fspac%lvoli(pos_voint),     &
               &                                                  khie = fspac%hierarchical_basis, mnode=max_num_nodes)
          fspac%lelem(ielem)%integ(ivar)%p => fspac%lvoli(pos_voint)

       end do

    end do

    call integration_faces_list( fspac%g_trian, fspac )

  end subroutine fem_space_fe_list_create

  !==================================================================================================
  subroutine dof_to_elnod( el2dof, idof, ivar, nobje, pnodob, lnode )
    implicit none
    ! Parameters
    integer(ip), intent(in)  :: el2dof(:,:)
    integer(ip), intent(in)  :: idof, ivar, nobje
    integer(ip), intent(in)  :: pnodob(nobje+1)
    integer(ip), intent(out) :: lnode

    ! Locals
    integer(ip) :: iobje, jnode

    do iobje = 1,nobje
       do jnode = pnodob(iobje),pnodob(iobje+1)-1
          if ( el2dof(jnode,ivar) == idof ) then
             lnode = jnode - pnodob(iobje) + 1
             return
          end if
       end do
    end do

    ! ERROR: not found !!!!
    assert(0==1)

  end subroutine dof_to_elnod

  subroutine fem_space_print (lunou, femsp )
    implicit none
    integer(ip)      , intent(in) :: lunou
    type(fem_space), intent(in) :: femsp

    integer(ip) :: ielem

    write (lunou,*) 'Number of materials: ', femsp%num_continuity
    write (lunou,*) 'Static condensation flag: ', femsp%static_condensation

    write (lunou,*) '****PRINT ELEMENT LIST INFO****'
    do ielem = 1, femsp%g_trian%num_elems
       write (lunou,*) '****PRINT ELEMENT ',ielem,' INFO****'
       call fem_element_print ( lunou, femsp%lelem(ielem) )
       write (lunou,*) '****END PRINT ELEMENT ',ielem,' INFO****'
    end do

    write (lunou,*) 'Number boundary faces: ', femsp%num_boundary_faces
    write (lunou,*) 'Number interior faces: ', femsp%num_interior_faces
    write (lunou,*) 'Boundary faces: ', femsp%boundary_faces(:)%face_object
    write (lunou,*) 'Interior faces: ', femsp%interior_faces(:)%face_object
    

  end subroutine fem_space_print

 !==================================================================================================
  subroutine fem_space_free ( f )
    implicit none
    type(fem_space), intent(inout) :: f
    integer(ip)                    :: i,j
    integer(ip) :: istat

    ! SB.alert : To be thought for the new structure triangulation
    ! lface_free
    ! if (allocated(f%lface)) then
    !    j = 1
    !    do i = 1, f%g_trian%num_objects
    !       if ( trian%objects(i)%dimension == f%g_trian%num_dims-1 ) then
    !          f%lface(i)%nei_elem  = 0
    !          f%lface(i)%pos_elem  = 0
    !          f%lface(i)%subface   = 0
    !          f%lface(i)%refinement_level = 0
    !          nullify (f%lface(i)%p_mat)
    !          nullify (f%lface(i)%p_vec)
    !          if (allocated(f%lface(i)%integ)) call memfree(f%lface(i)%integ,__FILE__,__LINE__)
    !          if ( trian%objects(i)%border /= -1) then
    !             do j=1,f%lelem(f%lface(i)%nei_elem(1))%p%nvars
    !                call array_free(f%lface(i)%o2n(j))
    !             end do
    !             call memfree(f%lface(i)%o2n,__FILE__,__LINE__)
    !          end if
    !          j = j+1
    !       end if
    !    end do

    !    deallocate(f%lface)
    !    ! Face integrators free
    !    do i = 1,f%cur_lfaci-1 
    !       call integ_free( f%lfaci(i) )
    !    end do
    !    f%cur_lfaci = 0
    !    call f%ht_pos_face_integrator%free
    ! end if

    

    do i = 1, f%g_trian%num_elems
!       nullify ( f%lelem(i)%problem )
       call memfree( f%lelem(i)%elem2dof,__FILE__,__LINE__)
       call memfree( f%lelem(i)%unkno,__FILE__,__LINE__)
       !call memfree( f%lelem(i)%jvars,__FILE__,__LINE__)
       if (allocated(f%lelem(i)%bc_code)) then
          call memfree(f%lelem(i)%bc_code,__FILE__,__LINE__)
       end if
       nullify ( f%lelem(i)%p_geo_info )
       nullify ( f%lelem(i)%p_mat )
       nullify ( f%lelem(i)%p_vec )
       do j = 1, f%lelem(i)%num_vars
          nullify ( f%lelem(i)%nodes_object(j)%p )          
       end do
       !if(allocated(f%lelem(i)%nodes_object)) call memfree(f%lelem(i)%nodes_object,__FILE__,__LINE__)
       if(allocated(f%lelem(i)%f_inf)) call memfree(f%lelem(i)%f_inf,__FILE__,__LINE__)
       if(allocated(f%lelem(i)%integ)) call memfree(f%lelem(i)%integ,__FILE__,__LINE__)
       !if(allocated(f%lelem(i)%iv))    call memfree(f%lelem(i)%iv   ,__FILE__,__LINE__)
       if(allocated(f%lelem(i)%continuity))    call memfree(f%lelem(i)%continuity   ,__FILE__,__LINE__)
       if(allocated(f%lelem(i)%order))    call memfree(f%lelem(i)%order   ,__FILE__,__LINE__)
       !if(allocated(f%lelem(i)%material))    call memfree(f%lelem(i)%iv   ,__FILE__,__LINE__)
       !if(allocated(f%lelem(i)%p_nod)) call memfree(f%lelem(i)%p_nod,__FILE__,__LINE__)
    end do
    
    deallocate ( f%mig_elems, stat=istat )
    check(istat==0)
    nullify ( f%lelem )

    do i = 1,f%pos_elmatvec%last()
       call array_free( f%lelmat(i) )
       call array_free( f%lelvec(i) )
    end do
    !deallocate ( f%lelmat )
    !deallocate ( f%lelvec )
    call f%pos_elmatvec%free
    call f%pos_elmatvec%free

    do i = 1,f%pos_volume_integrator%last()
       call volume_integrator_free( f%lvoli(i) )
    end do
    !deallocate ( f%lvoli )
    call f%pos_volume_integrator%free

    do i = 1,f%pos_elem_info%last()
       call fem_element_fixed_info_free (f%lelem_info(i))
    end do
    call f%pos_elem_info%free

    nullify ( f%g_trian )

    call memfree( f%void_list%p, __FILE__, __LINE__ )
    call memfree( f%void_list%l, __FILE__, __LINE__ )
    f%void_list%n = 0

    do i = 1, f%dof_handler%nblocks
       call memfree (f%object2dof(i)%p , __FILE__, __LINE__ )
       call memfree (f%object2dof(i)%l , __FILE__, __LINE__ )
    end do

    call memfree (f%ndofs , __FILE__, __LINE__ )

    deallocate( f%interior_faces, stat=istat)
    check ( istat == 0 )
    deallocate( f%boundary_faces, stat=istat)
    check ( istat == 0 )

  end subroutine fem_space_free

  subroutine get_p_faces ( fspac, trian )
    implicit none
    type(fem_space), intent(in)  :: fspac
    type(fem_triangulation), intent(inout)    :: trian

    !   integer(ip) :: i


    !   num_faces

    !   call memalloc ( trian%num_elems, num_faces, __FILE__,__LINE__ )

    !   do i=1,mesh%nelem
    !         mesh%p_face(i) = mesh%pnods(i) + fspac%lelem(i)%f_inf(1)%p%nobje_dim(mesh%ndime) -1
    !   end do

  end subroutine get_p_faces

  subroutine integration_faces_list( trian, femsp ) 
    implicit none
    ! Parameters
    type(fem_triangulation), intent(in)       :: trian 
    type(fem_space), intent(inout)               :: femsp

    integer(ip) :: count_int, count_bou, mat_i, mat_j, iobje, ielem, jelem, istat
    integer(ip) :: g_var, iprob, jprob, ivars, jvars

    ! integration faces (interior / boundary)
    ! The list of boundary faces includes all faces, whereas the interior ones are only those
    ! where we expect to integrate things (based on continuity flags)
    count_int = 0
    count_bou = trian%num_boundary_faces
    do iobje = 1, trian%num_objects
       if ( trian%objects(iobje)%dimension == trian%num_dims-1 ) then
          if ( trian%objects(iobje)%border == -1 ) then
             assert( trian%objects(iobje)%num_elems_around == 2 )
             ielem = trian%objects(iobje)%elems_around(1)
             jelem = trian%objects(iobje)%elems_around(2)
             iprob = femsp%lelem(ielem)%problem
             jprob = femsp%lelem(jelem)%problem
             do ivars = 1, femsp%dof_handler%problems(iprob)%p%nvars
                g_var = femsp%dof_handler%problems(iprob)%p%l2g_var(ivars)
                jvars = femsp%dof_handler%g2l_vars(g_var,jprob)
                mat_i = femsp%lelem(ielem)%continuity(ivars)
                mat_j = femsp%lelem(jelem)%continuity(jvars)
                if ( mat_i == 0 .or. mat_i /= mat_j ) then
                   count_int = count_int + 1
                   exit
                   !femsp%interior_faces(count_int) = iobje
                end if
             end do
          else
             assert( trian%objects(iobje)%num_elems_around == 1 )
          end if
       end if
    end do

    allocate( femsp%interior_faces(count_int), stat=istat)
    check ( istat == 0 )
    allocate( femsp%boundary_faces(count_bou), stat=istat)
    check ( istat == 0 )


    femsp%num_interior_faces = count_int
    femsp%num_boundary_faces = count_bou

    count_int = 0
    count_bou = 0
    do iobje = 1, trian%num_objects
       if ( trian%objects(iobje)%dimension == trian%num_dims-1 ) then
          if ( trian%objects(iobje)%border == -1 ) then
             assert( trian%objects(iobje)%num_elems_around == 2 )
             ielem = trian%objects(iobje)%elems_around(1)
             jelem = trian%objects(iobje)%elems_around(2)
             iprob = femsp%lelem(ielem)%problem
             jprob = femsp%lelem(jelem)%problem
             do ivars = 1, femsp%dof_handler%problems(iprob)%p%nvars
                g_var = femsp%dof_handler%problems(iprob)%p%l2g_var(ivars)
                jvars = femsp%dof_handler%g2l_vars(g_var,jprob)
                mat_i = femsp%lelem(ielem)%continuity(ivars)
                mat_j = femsp%lelem(jelem)%continuity(jvars)
                if ( mat_i == 0 .or. mat_i /= mat_j ) then
                   count_int = count_int + 1
                   femsp%interior_faces(count_int)%face_object = iobje
                   exit
                end if
             end do
          else
             assert( trian%objects(iobje)%num_elems_around == 1 )
             count_bou = count_bou + 1
             femsp%boundary_faces(count_bou)%face_object = iobje
          end if
       end if
    end do


  end subroutine integration_faces_list

end module fem_space_names

