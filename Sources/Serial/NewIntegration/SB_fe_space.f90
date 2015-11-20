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
module SB_fe_space_names
  use memor_names
  use allocatable_array_names
  use matrix_names
  use array_names
  use serial_scalar_matrix_names
  use serial_scalar_array_names
  use SB_matrix_array_assembler_names
  use SB_serial_scalar_matrix_array_assembler_names
  use types_names
  use reference_fe_names
  use triangulation_names
  use reference_fe_factory_names
  use integration_tools_names
  use migratory_element_names
  use conditions_names
  use hash_table_names
  use graph_names
  use sort_names
  implicit none
# include "debug.i90"
  private

  type, extends(migratory_element_t) :: SB_finite_element_t 
  private
  ! Reference element info          
  type(p_elem_topology_t) :: cell 
  class(reference_fe_t), pointer :: geometry_reference_fe
  class(reference_fe_t), pointer :: reference_fe
  type(SB_volume_integrator_t), pointer :: volume_integrator

  integer(ip) :: field_unknowns                  ! (V_K)^field_unknowns

  !type(volume_integrator_pointer_t) :: integ(:) ! Pointer to integration parameters

  ! Local to global 
  integer(ip)     , allocatable   :: elem2dof(:)   ! Map from elem to dof   

  ! Boundary conditions
  ! In finite_element? I think it should go to the fe_space
  integer(ip), allocatable :: bc_code(:,:)   ! Boundary Condition values

contains
  procedure :: create => create_fe
  procedure :: print_fe

  procedure :: get_reference_fe
  procedure :: get_volume_integrator

  procedure :: get_elem2dof

  procedure :: size   => finite_element_size
  procedure :: pack   => finite_element_pack
  procedure :: unpack => finite_element_unpack
end type SB_finite_element_t

type :: SB_fe_space_t
  private

  type(triangulation_t), pointer :: triangulation
  type(SB_finite_element_t), allocatable :: fe_array(:)
  !type(vector_space_t)  :: fe_function_space

  type(p_reference_fe_t), allocatable :: reference_fe_phy_list(:)
  type(p_reference_fe_t), allocatable :: reference_fe_geo_list(:)

  ! Integrator
  type(SB_p_volume_integrator_t), allocatable :: volume_integrator(:)

  integer(ip) :: number_dofs     

  ! Acceleration arrays
  type(list_2d_t)       :: vef2dof       ! An auxiliary array to accelerate some parts of the code

contains

  ! TBPs
  procedure :: create => create_fe_space
  !procedure :: free
  procedure :: print => print_fe_space

  procedure :: get_fe
  procedure :: initialize_volume_integrator
  procedure :: get_max_number_nodes

  procedure :: fill_dof_info
  procedure :: create_assembler

  procedure :: symbolic_setup_assembler

  !procedure :: generate_dof_numbering

  ! procedure :: fe_iterator ! Work to be done Alberto / Javier

end type SB_fe_space_t

! Data types
public :: SB_fe_space_t, SB_finite_element_t

contains
subroutine create_fe_space( this, triangulation, topology, fe_type, number_dimensions, order, field_unknowns, boundary_conditions, continuity )
 implicit none
 class(SB_fe_space_t), intent(out) :: this
 character(*), intent(in) :: topology, fe_type
 type(triangulation_t), target, intent(in) :: triangulation
 integer(ip), intent(in)  :: number_dimensions, order, field_unknowns
 logical, optional, intent(in) :: continuity
 type(conditions_t), intent(in)  :: boundary_conditions

 integer(ip) :: i, istat,iobje
 ! Point triangulation
 ! Create reference FE (now only one)

 this%triangulation => triangulation

 ! Assuming only one reference FE for the moment  SB.alert
 allocate( this%reference_fe_phy_list(1), stat=istat )
 allocate( this%reference_fe_geo_list(1), stat=istat )
 allocate( this%volume_integrator(1), stat=istat )
 allocate( this%fe_array(triangulation%num_elems), stat=istat )

 this%reference_fe_phy_list(1)%p => &
      start_reference_fe ( topology, fe_type, number_dimensions, order , continuity )

 ! Order 1 for geometri SB.alert
 this%reference_fe_geo_list(1)%p => &
      start_reference_fe ( topology, fe_type, number_dimensions, 1 , continuity )

 allocate( this%volume_integrator(1)%p )
 call this%volume_integrator(1)%p%create( this%reference_fe_phy_list(1)%p, this%reference_fe_geo_list(1)%p )

 ! Do triangulation
 do i = 1, triangulation%num_elems

    ! Create FEs
    call this%fe_array(i)%create( this%reference_fe_phy_list(1)%p, this%reference_fe_geo_list(1)%p, &
         this%volume_integrator(1)%p, triangulation%elems(i), field_unknowns)

 end do

 ! BCs
 do i = 1, triangulation%num_elems
    call memalloc(1,this%fe_array(i)%reference_fe%get_number_vefs(),this%fe_array(i)%bc_code,__FILE__,__LINE__, 0)
    do iobje = 1,this%fe_array(i)%reference_fe%get_number_vefs()
       this%fe_array(i)%bc_code(1,iobje) = boundary_conditions%code( 1,  this%triangulation%elems(i)%vefs(iobje) )
    end do
 end do

end subroutine create_fe_space




function create_assembler(this, diagonal_blocks_symmetric_storage,&
    diagonal_blocks_symmetric, diagonal_blocks_sign)
 implicit none
 class(SB_fe_space_t)          , intent(in) :: this
 class(SB_matrix_array_assembler_t), pointer :: create_assembler
 logical                        , intent(in) :: diagonal_blocks_symmetric_storage(:)
 logical                        , intent(in) :: diagonal_blocks_symmetric(:)
 integer(ip)                    , intent(in) :: diagonal_blocks_sign(:)

 class(matrix_t), pointer :: matrix
 class(array_t) , pointer :: array

 ! 1. Select dynamically the type of class(matrix_array_assembler_t), class(matrix_t) and class(vector_t)
 ! 2. Create class(matrix_t) and class(vector_t) accordingly to their dynamic type
 !    if (this%dof_descriptor%nblocks == 1) then
 allocate ( SB_serial_scalar_matrix_array_assembler_t :: create_assembler )
 allocate ( serial_scalar_matrix_t :: matrix )
 allocate ( serial_scalar_array_t  :: array )
 select type(matrix)
    class is(serial_scalar_matrix_t)
    call matrix%create( diagonal_blocks_symmetric_storage(1),& 
         diagonal_blocks_symmetric(1),& 
         diagonal_blocks_sign(1) )
    class default
    check(.false.)
 end select
 select type(array)
    class is(serial_scalar_array_t)
    call array%create(this%number_dofs)
    class default
    check(.false.)
 end select
 !else
 !   ! allocate ( par_block_matrix_array_assembler_t :: par_fe_space_create_matrix_array_assembler )
 !   allocate ( par_block_matrix_t :: matrix )
 !   allocate ( par_block_array_t  :: array )
 !   check(.false.)
 !   ! This itinerary still to be implemented ...
 !end if
 call create_assembler%set_matrix(matrix)
 call create_assembler%set_array(array)
end function create_assembler

subroutine create_dof_info ( fe_space )
 implicit none
 type(SB_fe_space_t), intent(inout) :: fe_space 
 call create_element_to_dof_and_ndofs( fe_space )
 call create_vef2dof( fe_space )
end subroutine create_dof_info


subroutine initialize_volume_integrator( this, max_order )
 implicit none
 ! Parameters
 class(SB_fe_space_t), intent(inout) :: this  
 integer(ip), optional, intent(in)  :: max_order
 integer(ip) :: i

 do i = 1, size(this%volume_integrator)
    call this%volume_integrator(i)%p%set_integration( max_order )
 end do

end subroutine initialize_volume_integrator

!subroutine create_dof_numbering
!end subroutine create_dof_numbering

subroutine print_fe_space ( this )
 implicit none
 class(SB_fe_space_t), intent(in) :: this
 integer(ip) :: i

 write(*,*) '******************** FINITE ELEMENT SPACE ********************'
 do i = 1, this%triangulation%num_elems
    write(*,*) '******************** FINITE ELEMENT (',i,') ********************'
    call this%fe_array(i)%print_fe()
 end do
 write(*,*) ' ********************** VEF2DOF **********************'
 call print_list_2d( 6, this%vef2dof ) 

end subroutine print_fe_space

function get_fe( this, i )
 implicit none
 class(SB_fe_space_t), target, intent(in) :: this
 integer(ip) :: i
 type(SB_finite_element_t), pointer :: get_fe
 get_fe => this%fe_array(i)
end function get_fe



!*********************************************************************************
! This subroutine generates the following DOF-related data for serial runs.
! 1) It generates DOF numbering and fills the elem2dof arrays (see explanation of 
!    the subroutine below)
! 2) It generates the vef2dof array ( see explanation of the subroutine below)
! NOTE: In order to understand the following subroutines, it is useful to know that 
! currently, the *par_fe_space* (parallel) has a pointer to a *fe_space*, which 
! has not info about the distributed environment (since it is a serial object),
! but it includes the *num_elems* local elements + *num_ghosts* ghost elements.
! The ghost part is filled in *par_fe_space_create*. Thus, we know if we have a
! local or ghost element below checking whether *ielem* <=  *num_elems* (local)
! or not. We are using this trick below in some cases, to be able to use the same
! subroutines for parallel runs, and fill properly ghost info. For serial runs, 
! the only thing is that *ielem* <=  *num_elems* always, so it never goes to the
! ghost element part.
!*********************************************************************************
subroutine fill_dof_info ( fe_space )
 implicit none
 class(SB_fe_space_t), intent(inout) :: fe_space 
 !write(*,*) '********create_element_to_dof_and_ndofs********'
 call create_element_to_dof_and_ndofs( fe_space )
 call create_vef2dof( fe_space )
end subroutine fill_dof_info

function get_max_number_nodes( this)
 implicit none
 class(SB_fe_space_t), intent(in) :: this 
 integer(ip) :: get_max_number_nodes
 get_max_number_nodes = 100 ! SB.alert 
end function get_max_number_nodes



!*********************************************************************************
! This subroutine takes the triangulation and fills the element2dof structure at every 
! finite element, i.e., it labels all dofs related to local elements (not ghost), after a 
! count-list procedure, and puts the number of dofs in ndofs structure (per block).
! Note 1: The numbering is per every block independently, where the blocks are 
! defined at the dof_descriptor. A global dof numbering is not needed in the code, 
! when blocks are being used.
!*********************************************************************************
subroutine create_element_to_dof_and_ndofs( fe_space ) 
 implicit none
 ! Parameters
 class(SB_fe_space_t)     , intent(inout) :: fe_space 

 ! Local variables
 integer(ip) :: iprob, count, iobje, ielem, jelem, nvapb
 integer(ip) :: obje_l, inode, l_node, elem_ext, obje_ext, prob_ext, inode_ext, inode_l
 integer(ip) :: mater, order, nnode
 integer(ip) :: touch(1,2)

 integer(ip), allocatable :: o2n(:)

 count = 0

 mater = 1

 call memalloc ( fe_space%get_max_number_nodes(), o2n, __FILE__, __LINE__ )

 ! Part 1: Put DOFs on VEFs, taking into account that DOFs only belong to VEFs when we do not
 ! enforce continuity (continuity(ielem) /= 0). We go through all objects, elements around the
 ! object, variables of the element, and if for the value of continuity of this element no 
 ! DOFs have already been added, we add them and touch this object for this continuity value.
 ! In FEMPAR, continuity is an elemental value. If it is different from 0, the nodes/DOFs 
 ! geometrically on the interface belong to the interface objects (VEFs). Next, we only
 ! enforce continuity for elements with same continuity value (mater below), in order to 
 ! allow for situations in which we want to have continuity in patches and discontinuity among 
 ! patches based on physical arguments (every patch would involve its own value of continuity).
 ! For hp-adaptivity, we could consider the value in continuity to be p (order) and as a result
 ! not to enforce continuity among elements with different order SINCE it would lead to ERROR
 ! to enforce continuity among elements of different order.
 do iobje = 1, fe_space%triangulation%num_vefs          
    touch = 0
    !write(*,*) 'loop iobje',iobje
    do ielem = 1, fe_space%triangulation%vefs(iobje)%num_elems_around
       jelem = fe_space%triangulation%vefs(iobje)%elems_around(ielem)
       !write(*,*) 'jelem',jelem
       if ( jelem <= fe_space%triangulation%num_elems ) then ! Local elements
          if ( fe_space%fe_array(jelem)%reference_fe%get_continuity() ) then
             !mater = fe_space%fe_array(jelem)%continuity(g_var) ! SB.alert : continuity can be used as p 
             mater = 1 ! For the moment SB.alert
             do obje_l = 1, fe_space%triangulation%elems(jelem)%num_vefs
                if ( fe_space%triangulation%elems(jelem)%vefs(obje_l) == iobje ) exit
             end do
             if ( fe_space%fe_array(jelem)%bc_code(1,obje_l) == 0 ) then
                if ( touch(mater,1) == 0 ) then                            
                   touch(mater,1) = jelem
                   touch(mater,2) = obje_l
                   !write(*,*) 'obje_l new dofs',obje_l
                   call put_new_vefs_dofs_in_vef_of_element ( fe_space%triangulation, fe_space, jelem, count, obje_l )
                else
                   !write(*,*) 'obje_l old dofs',obje_l
                   call put_existing_vefs_dofs_in_vef_of_element ( fe_space%triangulation, fe_space, touch, mater, iobje, jelem, o2n, obje_l )
                end if
             end if
          end if
       else ! Ghost elements
          if ( fe_space%fe_array(jelem)%reference_fe%get_continuity() ) then
             !mater = fe_space%fe_array(jelem)%continuity(g_var) ! SB.alert : continuity can be used as p 
             mater = 1
             do obje_l = 1, fe_space%triangulation%elems(jelem)%num_vefs
                if ( fe_space%triangulation%elems(jelem)%vefs(obje_l) == iobje ) exit
             end do
             if ( touch(mater,1) /= 0) then
                call put_existing_vefs_dofs_in_vef_of_element ( fe_space%triangulation, fe_space, touch, mater, iobje, jelem, o2n, obje_l )
             end if
          end if
       end if
    end do
 end do

 ! Part 2: Put DOFs on nodes belonging to the volume object (element). For cG we only do that when 
 ! static condensation is not active. Static condensation is for all variables, elements, etc. BUT
 ! it cannot be used with dG. The following algorithm is ASSUMING that this is the case, and we are
 ! not using dG + static condensations. In any case, when creating the fe_space there is an 
 ! automatic check for satisfying that.
 ! No check about strong Dirichlet boundary conditions, because they are imposed weakly in dG, and
 ! never appear in interior nodes in cG.
 ! if ( ( .not. fe_space%static_condensation )  ) then
 !    do ielem = 1, fe_space%triangulation%num_elems
 !       iobje = fe_space%triangulation%elems(ielem)%num_vefs+1
 !       do inode = fe_space%fe_array(ielem)%nodes_per_vef%p%p(iobje), &
 !            &     fe_space%fe_array(ielem)%nodes_per_vef%p%p(iobje+1)-1 
 !          l_node = fe_space%fe_array(ielem)%nodes_per_vef%p%l(inode)
 !          count = count +1
 !          fe_space%fe_array(ielem)%elem2dof(l_node) = count
 !       end do
 !    end do
 ! end if

 ! Part 3: Assign total number of dofs created to fem space object
 fe_space%number_dofs = count


end subroutine create_element_to_dof_and_ndofs


!*********************************************************************************
! This subroutine takes the finite element space and fills the vef2dof structure. 
! The vef2dof structure puts on top of VEFs the DOFs that are meant to be continuous 
! between elements with the same continuity label. As an example, when using dG only, 
! vef2dof is void. It is more an acceleration array than a really needed structure, 
! but it is convenient when creating the dof graph. 
!*********************************************************************************
subroutine create_vef2dof ( fe_space ) 
 implicit none
 ! Parameters
 type(SB_fe_space_t), intent(inout) :: fe_space 

 ! Local variables
 integer(ip) :: iprob, count, iobje, ielem, jelem, nvapb
 integer(ip) :: obje_l, inode, l_node, mater, istat
 integer(ip) :: touch!(1,fe_space%num_continuity)

 ! Part 1: Count DOFs on VEFs, using the notion of continuity described above (in elem2dof)
 fe_space%vef2dof%n1 = fe_space%triangulation%num_vefs
 fe_space%vef2dof%n2 = 1
 call memalloc ( fe_space%triangulation%num_vefs+1, fe_space%vef2dof%p, __FILE__, __LINE__, 0 )
 do iobje = 1, fe_space%triangulation%num_vefs
    touch = 0
    do ielem = 1, fe_space%triangulation%vefs(iobje)%num_elems_around
       jelem = fe_space%triangulation%vefs(iobje)%elems_around(ielem)
       if ( jelem <= fe_space%triangulation%num_elems ) then 
          !iprob = fe_space%fe_array(jelem)%problem
          mater = 1
          if ( mater /= 0 ) then
             if (touch == 0) then
                touch = 1
                !if ( touch(g_var,mater) == 0 ) then
                !   touch(g_var,mater) = 1
                do obje_l = 1, fe_space%triangulation%elems(jelem)%num_vefs
                   if ( fe_space%triangulation%elems(jelem)%vefs(obje_l) == iobje ) exit
                end do
                if ( fe_space%fe_array(jelem)%bc_code(1,obje_l) == 0 ) then
                   fe_space%vef2dof%p(iobje+1) = fe_space%vef2dof%p(iobje+1) + &
                        fe_space%fe_array(jelem)%reference_fe%get_number_interior_nodes_vef(obje_l)
                   !+ fe_space%fe_array(jelem)%nodes_per_vef%p%p(obje_l+1) &
                   !     & - fe_space%fe_array(jelem)%nodes_per_vef%p%p(obje_l)
                end if
             end if
          end if
       end if
    end do
 end do

 fe_space%vef2dof%p(1) = 1
 do iobje = 2, fe_space%triangulation%num_vefs+1
    fe_space%vef2dof%p(iobje) = fe_space%vef2dof%p(iobje) + fe_space%vef2dof%p(iobje-1)
 end do

 call memalloc ( fe_space%vef2dof%p(fe_space%triangulation%num_vefs+1)-1, 3, fe_space%vef2dof%l, __FILE__, __LINE__ )

 ! Part 2: List DOFs on VEFs, using the notion of continuity described above (in elem2dof)
 ! We note that the vef2dof%l(:,X) is defined for X = 1,2,3
 ! vef2dof%l(:,1) : DOF LID
 ! vef2dof%l(:,2) : Variable GID associated to that DOF
 ! vef2dof%l(:,3) : Continuity value associated to that DOF (to enforce continuity)
 count = 0
 do iobje = 1, fe_space%triangulation%num_vefs
    touch = 0
    do ielem = 1, fe_space%triangulation%vefs(iobje)%num_elems_around
       jelem = fe_space%triangulation%vefs(iobje)%elems_around(ielem)
       if ( jelem <= fe_space%triangulation%num_elems ) then 
          !iprob = fe_space%fe_array(jelem)%problem
          mater = 1
          if ( mater /= 0) then
             if (touch == 0) then
                touch = 1
                !if ( touch(g_var,mater) == 0 ) then
                !   touch(g_var,mater) = 1
                do obje_l = 1, fe_space%triangulation%elems(jelem)%num_vefs
                   if ( fe_space%triangulation%elems(jelem)%vefs(obje_l) == iobje ) exit
                end do
                if ( fe_space%fe_array(jelem)%bc_code(1,obje_l) == 0 ) then
                   do inode = 1, fe_space%fe_array(jelem)%reference_fe%get_number_interior_nodes_vef(obje_l)
                      !fe_space%fe_array(jelem)%nodes_per_vef%p%p(obje_l), &
                      &     !fe_space%fe_array(jelem)%nodes_per_vef%p%p(obje_l+1)-1 
                           l_node = fe_space%fe_array(jelem)%reference_fe%get_interior_node_vef(inode,obje_l)
                      !fe_space%fe_array(jelem)%nodes_per_vef%p%l(inode)
                      count = count + 1
                      fe_space%vef2dof%l(count,1) = fe_space%fe_array(jelem)%elem2dof(l_node)
                      !fe_space%vef2dof%l(count,2) = 1
                      !fe_space%vef2dof%l(count,3) = mater
                   end do
                end if
             end if
          end if
       end if
    end do
 end do
end subroutine create_vef2dof

!*********************************************************************************
! Auxiliary function that generates new DOFs and put them in a particular VEF of a given element
!*********************************************************************************
subroutine put_new_vefs_dofs_in_vef_of_element ( trian, fe_space, jelem, &
    count, obje_l )
 implicit none
 ! Parameters
 type(triangulation_t), intent(in)         :: trian 
 type(SB_fe_space_t), intent(inout)              :: fe_space 
 integer(ip), intent(inout)                  :: count
 integer(ip), intent(in)                     :: jelem, obje_l

 ! Local variables
 integer(ip) :: inode, l_node

 do inode = 1,fe_space%fe_array(jelem)%reference_fe%get_number_interior_nodes_vef(obje_l)
    !fe_space%fe_array(jelem)%nodes_per_vef%p%p(obje_l), &
    !         &     fe_space%fe_array(jelem)%nodes_per_vef%p%p(obje_l+1)-1 
    l_node = fe_space%fe_array(jelem)%reference_fe%get_interior_node_vef(inode,obje_l)
    !l_node = fe_space%fe_array(jelem)%nodes_per_vef%p%l(inode)
    count = count + 1
    !write (*,*) '****PUT DOF**** (elem,obj_l,obj_g,node,idof) ',jelem,obje_l,l_node,count
    fe_space%fe_array(jelem)%elem2dof(l_node) = count
 end do

end subroutine put_new_vefs_dofs_in_vef_of_element

!*********************************************************************************
! Auxiliary function that puts existing DOFs in a particular VEF of a given element
!*********************************************************************************
subroutine put_existing_vefs_dofs_in_vef_of_element ( trian, fe_space, touch, mater, iobje, jelem, &
    o2n, obje_l )
 implicit none
 ! Parameters
 type(triangulation_t), intent(in)           :: trian 
 type(SB_fe_space_t), intent(inout)      :: fe_space
 integer(ip), intent(in)                     :: touch(:,:), mater, iobje, jelem, obje_l
 integer(ip), intent(out)                    :: o2n(:)

 ! Local variables
 integer(ip) :: elem_ext, obje_ext, prob_ext
 integer(ip) :: nnode, order, inode, l_node, inode_ext, inode_l

 elem_ext = touch(mater,1)
 obje_ext = touch(mater,2)
 !prob_ext = fe_space%fe_array(elem_ext)%problem
 prob_ext = 1
 !write (*,*) '****EXTRACT DOF**** (object)', iobje, ' FROM: (elem,obj_l) ',elem_ext,obje_ext, ' TO  : (elem,obj_l)', jelem,obje_l
 nnode = fe_space%fe_array(elem_ext)%reference_fe%get_number_interior_nodes_vef(obje_ext)
 !fe_space%fe_array(elem_ext)%nodes_per_vef%p%p(obje_ext+1) &
 !     &  -fe_space%fe_array(elem_ext)%nodes_per_vef%p%p(obje_ext) 
 !write(*,*) 'jelem',jelem
 !write(*,*) 'nnode',nnode
 if ( nnode > 0) then  
    order = fe_space%fe_array(elem_ext)%reference_fe%get_order()
    !write(*,*) 'order',order
    if ( trian%vefs(iobje)%dimension == trian%num_dims .and. &
         & nnode ==  (order+1)**trian%num_dims ) then
       order = order    ! hdG case
    elseif ( nnode ==  (order-1)**trian%vefs(iobje)%dimension ) then
       order = order -2 ! cG case
    else
       assert ( 0 == 1) ! SB.alert : Other situations possible when dG_continuity, cdG, hp-adaptivity ?
    end if
    call fe_space%fe_array(jelem)%reference_fe%permute_nodes_per_vef( &
         & fe_space%fe_array(elem_ext)%reference_fe,                  &
         & o2n,obje_ext,obje_l,                                                     &
         & trian%elems(elem_ext)%vefs,                                              &
         & trian%elems(jelem)%vefs,                                                 &
         & trian%vefs(iobje)%dimension,                                             &
         & order )
    do inode = 1, fe_space%fe_array(elem_ext)%reference_fe%get_number_interior_nodes_vef(obje_ext)
       !nodes_per_vef%p%p(obje_ext+1) - &
       !            fe_space%fe_array(elem_ext)%nodes_per_vef%p%p(obje_ext)
       inode_ext = fe_space%fe_array(elem_ext)%reference_fe%get_interior_node_vef(inode,obje_ext)
       !l_node = fe_space%fe_array(elem_ext)%nodes_per_vef%p%p(obje_ext) + inode - 1
       !inode_ext = fe_space%fe_array(elem_ext)%nodes_per_vef%p%l(l_node )
       inode_l = fe_space%fe_array(jelem)%reference_fe%get_interior_node_vef(o2n(inode),obje_l)
       !l_node = fe_space%fe_array(jelem)%nodes_per_vef%p%p(obje_l) + o2n(inode) - 1
       !inode_l = fe_space%fe_array(jelem)%nodes_per_vef%p%l(l_node)
       fe_space%fe_array(jelem)%elem2dof(inode_l) = fe_space%fe_array(elem_ext)%elem2dof(inode_ext)
    end do ! SB.alert : 1) face object for cG and hdG, where the face must have all their nodes
    !                   2) corner / edge only for cG
    !            * Never here for dG, continuity interface, hanging objects, etc.


 end if
end subroutine put_existing_vefs_dofs_in_vef_of_element





subroutine create_fe( this, reference_fe_unk, reference_fe_geo, volume_integrator, cell, field_unknowns )
 implicit none
 class(SB_finite_element_t), intent(inout)  :: this
 class(reference_fe_t), target, intent(in) :: reference_fe_unk, reference_fe_geo
 type(SB_volume_integrator_t), target, intent(in) :: volume_integrator
 type(elem_topology_t), target, intent(in) :: cell
 integer(ip), intent(in) :: field_unknowns

 this%reference_fe => reference_fe_unk
 this%geometry_reference_fe => reference_fe_geo
 this%volume_integrator => volume_integrator
 this%cell%p => cell  

 call memalloc( reference_fe_unk%get_number_nodes(), this%elem2dof, __FILE__, __LINE__)
 this%elem2dof = 0

end subroutine create_fe

subroutine print_fe( this )
 implicit none
 class(SB_finite_element_t), intent(in)  :: this

 write(*,*) '********************REFERENCE FE UNK********************'
 call this%reference_fe%print()
 write(*,*) '********************REFERENCE FE UNK********************'
 call this%geometry_reference_fe%print()
 write(*,*) '********************CELL********************'
 call element_print( 6, this%cell%p )
 write(*,*) '********************ELEMENT 2 DOF********************'
 write(*,*), this%elem2dof
 write(*,*) '********************BC CODE********************'
 write(*,*), this%bc_code


end subroutine print_fe

function get_reference_fe ( this )
 implicit none
 class(SB_finite_element_t), target, intent(in)  :: this
 class(reference_fe_t), pointer :: get_reference_fe
 get_reference_fe => this%reference_fe
end function get_reference_fe

function get_volume_integrator ( this )
 implicit none
 class(SB_finite_element_t), target, intent(in)  :: this
 type(SB_volume_integrator_t), pointer :: get_volume_integrator
 get_volume_integrator => this%volume_integrator
end function get_volume_integrator

function get_elem2dof ( this )
 implicit none
 class(SB_finite_element_t), target, intent(in)  :: this
 integer(ip), pointer :: get_elem2dof(:)
 get_elem2dof => this%elem2dof
end function get_elem2dof

subroutine finite_element_size (my, n)
 implicit none
 class(SB_finite_element_t), intent(in)  :: my
 integer(ip)            , intent(out) :: n

end subroutine finite_element_size

subroutine finite_element_pack (my, n, buffer)
 implicit none
 class(SB_finite_element_t), intent(in)  :: my
 integer(ip)            , intent(in)   :: n
 integer(ieep)            , intent(out)  :: buffer(n)

end subroutine finite_element_pack

subroutine finite_element_unpack(my, n, buffer)
 implicit none
 class(SB_finite_element_t), intent(inout) :: my
 integer(ip)            , intent(in)     :: n
 integer(ieep)            , intent(in)     :: buffer(n)

end subroutine finite_element_unpack

subroutine symbolic_setup_assembler(this,matrix_array_assembler)
 implicit none
 class(SB_fe_space_t)        , intent(in)    :: this
 class(SB_matrix_array_assembler_t) , intent(inout) :: matrix_array_assembler

 ! Polymorphic matrix 
 class(matrix_t), pointer :: matrix

 matrix => matrix_array_assembler%get_matrix()
 select type(matrix)
    class is(serial_scalar_matrix_t)
    call setup_dof_graph_from_block_row_col_identifiers ( this, matrix%graph )
    class default
    check(.false.)
 end select

end subroutine symbolic_setup_assembler







!*********************************************************************************
! This subroutine takes the fe_space and creates a dof_graph. 
! The dof_graph includes both the coupling by continuity like in continuous Galerkin 
! methods, and the coupling by face terms (of discontinuous Galerkin type). The algorithm 
! considers both the case with static condensation and without it. In order to call this 
! subroutine, we need to compute first element2dof and vef2dof arrays.
! 3) It generates the local graph (to be extended in *par_create_global_dof_info_names*
!    to put additional DOFs due to face integration coupling with DOFs from ghost 
!    elements) (see explanation of the subroutine below and *ghost_dofs_by_integration*
!    in *par_create_global_dof_info_names* for more insight)
!*********************************************************************************
subroutine setup_dof_graph_from_block_row_col_identifiers(  fe_space, dof_graph ) 
 implicit none
 ! Parameters
 type(SB_fe_space_t)       , intent(in)     :: fe_space 
 type(graph_t)          , intent(inout)  :: dof_graph

 ! Local variables
 integer(ip) :: iprob, count, iobje, ielem, jelem, nvapb, inter, inode, l_node
 integer(ip) :: idof, jdof, int_i, int_j, istat, jnode, job_g, jobje, touch
 integer(ip) :: l_dof, m_dof, m_node, posi, posf, l_mat, m_mat, knode

 integer(ip) :: nvapbi, nvapbj, nnode, i, iface, jprob, l_faci, l_facj, ic


 integer(ip), allocatable :: aux_ia(:)
 type(hash_table_ip_ip_t) :: visited

 touch = 1

 ! Initialize
 dof_graph%nv  = fe_space%number_dofs
 dof_graph%nv2 = fe_space%number_dofs
 call memalloc( dof_graph%nv+1, dof_graph%ia, __FILE__,__LINE__ )
 dof_graph%ia = 0

 ! COUNT PART
 call count_nnz_dofs_vefs_vs_dofs_vefs_vol_by_continuity ( fe_space%triangulation, fe_space, dof_graph )
 !call count_nnz_dofs_vol_vs_dofs_vefs_vol_by_continuity ( fe_space%triangulation, fe_space, dof_graph ) 

 dof_graph%ia(1) = 1
 do idof = 2, fe_space%number_dofs+1
    dof_graph%ia(idof) = dof_graph%ia(idof) + dof_graph%ia(idof-1)
 end do

 call memalloc ( dof_graph%ia(fe_space%number_dofs+1)-1, dof_graph%ja, __FILE__, __LINE__ )

 ! LIST PART
 call memalloc( dof_graph%nv+1, aux_ia, __FILE__,__LINE__ )
 aux_ia = dof_graph%ia
 !


 call list_nnz_dofs_vefs_vs_dofs_vefs_vol_by_continuity ( fe_space%triangulation, fe_space, dof_graph, aux_ia ) 
 !call list_nnz_dofs_vol_vs_dofs_vefs_vol_by_continuity ( fe_space%triangulation, fe_space, dof_graph, aux_ia )


 do idof = 1, fe_space%number_dofs
    ! Order increasingly column identifiers of current row 
    ! using heap sort algorithm
    posi = dof_graph%ia(idof)
    posf = dof_graph%ia(idof+1)-1
    call sort(posf-posi+1,dof_graph%ja(posi:posf))
 end do
 call dof_graph%print( 6)

 call memfree (aux_ia,__FILE__,__LINE__)

end subroutine setup_dof_graph_from_block_row_col_identifiers



!*********************************************************************************
! Count NNZ (number of nonzero entries) for DOFs on the interface (VEFs) of elements against
! both interior and interface nodes.
!*********************************************************************************
subroutine count_nnz_dofs_vefs_vs_dofs_vefs_vol_by_continuity (  trian, fe_space, dof_graph )  
 implicit none
 ! Parameters
 type(triangulation_t), intent(in)         :: trian 
 type(SB_fe_space_t), intent(in)                 :: fe_space 
 type(graph_t), intent(inout)                :: dof_graph

 ! Local variables
 type(hash_table_ip_ip_t) :: visited
 integer(ip) :: idof, ielem, inode, iobje, iprob, istat
 integer(ip) :: jdof, jelem, job_g, jobje, l_dof, l_mat
 integer(ip) :: l_node,  m_dof, m_mat,  nvapb, touch

 do iobje = 1, trian%num_vefs             
    if ( fe_space%vef2dof%p(iobje+1)-fe_space%vef2dof%p(iobje) > 0) then
       call visited%init(100) 
       do ielem = 1, trian%vefs(iobje)%num_elems_around
          jelem = trian%vefs(iobje)%elems_around(ielem)
          if ( jelem <= trian%num_elems ) then
             do jobje = 1, trian%elems(jelem)%num_vefs
                job_g = trian%elems(jelem)%vefs(jobje)
                call visited%put(key=job_g, val=touch, stat=istat)
                if ( istat == now_stored ) then   ! interface-interface
                   do idof = fe_space%vef2dof%p(iobje), fe_space%vef2dof%p(iobje+1)-1
                      l_dof = fe_space%vef2dof%l(idof,1)
                      l_mat = 1
                      do jdof = fe_space%vef2dof%p(job_g), fe_space%vef2dof%p(job_g+1)-1
                         m_dof = fe_space%vef2dof%l(jdof,1)
                         m_mat = 1
                         if ( .not. dof_graph%symmetric_storage ) then
                            dof_graph%ia(l_dof+1) = &
                                 & dof_graph%ia(l_dof+1) + 1
                         else
                            if ( m_dof >= l_dof ) then
                               dof_graph%ia(l_dof+1) = &
                                    & dof_graph%ia(l_dof+1) + 1
                            end if
                         end if
                      end do
                   end do
                end if
             end do
             !end do
             !if (.not.fe_space%static_condensation) then  ! interface-interior
             !iprob = fe_space%fe_array(jelem)%problem
             do idof = fe_space%vef2dof%p(iobje), fe_space%vef2dof%p(iobje+1)-1
                l_dof = fe_space%vef2dof%l(idof,1)
                l_mat = fe_space%vef2dof%l(idof,3)
                m_mat = 1
                if ( .not. dof_graph%symmetric_storage ) then
                   dof_graph%ia(l_dof+1) =  dof_graph%ia(l_dof+1) &
                        + fe_space%fe_array(jelem)%reference_fe%get_number_interior_nodes_vef(jobje)
                   !                         & + fe_space%fe_array(jelem)%nodes_per_vef%p%p(jobje+1) &
                   !                         & - fe_space%fe_array(jelem)%nodes_per_vef%p%p(jobje)
                else 
                   do inode = 1, fe_space%fe_array(jelem)%reference_fe%get_number_interior_nodes_vef(jobje)
                      !fe_space%fe_array(jelem)%nodes_per_vef%p%p(jobje), &
                      !  & fe_space%fe_array(jelem)%nodes_per_vef%p%p(jobje+1)-1
                      l_node = fe_space%fe_array(jelem)%reference_fe%get_interior_node_vef(inode,jobje)
                      !l_node = fe_space%fe_array(jelem)%nodes_per_vef%p%l(inode)
                      m_dof = fe_space%fe_array(jelem)%elem2dof(l_node)
                      if ( m_dof >= l_dof ) then
                         dof_graph%ia(l_dof+1) = &
                              & dof_graph%ia(l_dof+1) + 1
                      end if
                   end do
                end if
             end do
          end if
       end do
       call visited%free
    end if
 end do
end subroutine count_nnz_dofs_vefs_vs_dofs_vefs_vol_by_continuity

!*********************************************************************************
! List NNZ (number of nonzero entries) for DOFs on the interface (VEFs) of elements against
! both interior and interface nodes.
!*********************************************************************************
subroutine list_nnz_dofs_vefs_vs_dofs_vefs_vol_by_continuity ( trian, fe_space, dof_graph, aux_ia )  
  implicit none
  ! Parameters
  type(triangulation_t), intent(in)         :: trian 
  type(SB_fe_space_t), intent(in)                 :: fe_space 
  type(graph_t), intent(inout)                :: dof_graph
  integer(ip), intent(inout)                :: aux_ia(:)

  ! Local variables
  type(hash_table_ip_ip_t) :: visited
  integer(ip) :: idof, ielem, inode, iobje, iprob, istat
  integer(ip) :: jdof, jelem, job_g, jobje, l_dof, l_mat
  integer(ip) :: l_node, m_dof, m_mat, nvapb, touch, count, ic

  count = 0
  do iobje = 1, trian%num_vefs 
     if ( fe_space%vef2dof%p(iobje+1)-fe_space%vef2dof%p(iobje) > 0) then
        call visited%init(100) 
        do ielem = 1, trian%vefs(iobje)%num_elems_around
           jelem = trian%vefs(iobje)%elems_around(ielem)
           if ( jelem <= trian%num_elems ) then 
              do jobje = 1, trian%elems(jelem)%num_vefs
                 job_g = trian%elems(jelem)%vefs(jobje)
                 call visited%put(key=job_g, val=touch, stat=istat)
                 if ( istat == now_stored ) then  ! interface-interface
                    do idof = fe_space%vef2dof%p(iobje), fe_space%vef2dof%p(iobje+1)-1
                       l_dof = fe_space%vef2dof%l(idof,1)
                       do jdof = fe_space%vef2dof%p(job_g), fe_space%vef2dof%p(job_g+1)-1
                          m_dof = fe_space%vef2dof%l(jdof,1)
                          if ( .not. dof_graph%symmetric_storage ) then
                             ic = aux_ia(l_dof)
                             dof_graph%ja(ic) = m_dof
                             aux_ia(l_dof) = aux_ia(l_dof)+1
                          else 
                             if ( m_dof >= l_dof ) then
                                ic = aux_ia(l_dof)
                                dof_graph%ja(ic) = m_dof
                                aux_ia(l_dof) = aux_ia(l_dof)+1
                             end if
                          end if
                       end do
                    end do
                 end if
              end do
              !end do
              !if (.not.fe_space%static_condensation) then  ! interface-interior
              !   iprob = fe_space%fe_array(jelem)%problem
              do idof = fe_space%vef2dof%p(iobje), fe_space%vef2dof%p(iobje+1)-1
                 l_dof = fe_space%vef2dof%l(idof,1)
                 l_mat = fe_space%vef2dof%l(idof,3)
                 m_mat = 1
                 if ( .not. dof_graph%symmetric_storage ) then
                    do inode = 1, fe_space%fe_array(jelem)%reference_fe%get_number_interior_nodes_vef(jobje)
                       !fe_space%fe_array(jelem)%nodes_per_vef%p%p(jobje), &
                       ! & fe_space%fe_array(jelem)%nodes_per_vef%p%p(jobje+1)-1
                       l_node = fe_space%fe_array(jelem)%reference_fe%get_interior_node_vef(inode,jobje)
                       !l_node = fe_space%fe_array(jelem)%nodes_per_vef%p%l(inode)
                       m_dof = fe_space%fe_array(jelem)%elem2dof(l_node)
                       ic = aux_ia(l_dof)
                       dof_graph%ja(ic) = m_dof
                       aux_ia(l_dof) = aux_ia(l_dof)+1
                    end do
                 else 
                    do inode = 1, fe_space%fe_array(jelem)%reference_fe%get_number_interior_nodes_vef(jobje)
                       !fe_space%fe_array(jelem)%nodes_per_vef%p%p(jobje), &
                       !  & fe_space%fe_array(jelem)%nodes_per_vef%p%p(jobje+1)-1
                       l_node = fe_space%fe_array(jelem)%reference_fe%get_interior_node_vef(inode,jobje)
                       !l_node = fe_space%fe_array(jelem)%nodes_per_vef%p%l(inode)
                       m_dof = fe_space%fe_array(jelem)%elem2dof(l_node)
                       if ( m_dof >= l_dof ) then
                          ic = aux_ia(l_dof)
                          dof_graph%ja(ic) = m_dof
                          aux_ia(l_dof) = aux_ia(l_dof)+1
                       end if
                    end do
                 end if
              end do
           end if
        end do
        call visited%free
     end if
  end do

end subroutine list_nnz_dofs_vefs_vs_dofs_vefs_vol_by_continuity



! !*********************************************************************************
! ! Count NNZ (number of nonzero entries) for DOFs on the interface (VEFs) of elements against
! ! both interior and interface nodes.
! !*********************************************************************************
! subroutine count_nnz_dofs_vol_vs_dofs_vefs_vol_by_continuity ( iblock, jblock, trian, fe_space, dof_graph )  
! implicit none
! ! Parameters
! integer(ip), intent(in)                     :: iblock, jblock
! type(triangulation_t), intent(in)         :: trian 
! type(SB_fe_space_t), intent(in)                 :: fe_space 
! type(graph_t), intent(inout)                :: dof_graph

! ! Local variables
! integer(ip) :: ielem, inode, int_i, iobje, iprob,  jdof, jnode, job_g
! integer(ip) :: jobje, l_dof, l_mat, l_node,  m_dof, m_mat
! integer(ip) :: m_node, nvapbi, nvapbj

! ! As commented for elem2dof, static condensation is false for dG, by construction of the 
! ! fem space.
! if (.not.fe_space%static_condensation) then
! do ielem  = 1, trian%num_elems
!   iobje = trian%elems(ielem)%num_vefs+1
!   iprob = fe_space%fe_array(ielem)%problem
!   ! Interior - interior 
!   do inode = fe_space%fe_array(ielem)%nodes_per_vef%p%p(iobje), &
!        & fe_space%fe_array(ielem)%nodes_per_vef%p%p(iobje+1)-1
!      l_node = fe_space%fe_array(ielem)%nodes_per_vef%p%l(inode)
!      l_dof = fe_space%fe_array(ielem)%elem2dof(l_node)
!      if ( .not. dof_graph%symmetric_storage ) then
!         dof_graph%ia(l_dof+1) =  dof_graph%ia(l_dof+1) &
!              &  + fe_space%fe_array(ielem)%nodes_per_vef%p%p(iobje+1) &
!              & - fe_space%fe_array(ielem)%nodes_per_vef%p%p(iobje)
!      else  
!         do jnode = fe_space%fe_array(ielem)%nodes_per_vef%p%p(iobje), &
!              & fe_space%fe_array(ielem)%nodes_per_vef%p%p(iobje+1)-1
!            m_node = fe_space%fe_array(ielem)%nodes_per_vef%p%l(jnode)
!            m_dof = fe_space%fe_array(ielem)%elem2dof(m_node)
!            if ( m_dof >= l_dof ) then
!               dof_graph%ia(l_dof+1) = &
!                    & dof_graph%ia(l_dof+1) + 1
!            end if
!         end do
!      end do
!   end if
! end do
! l_mat = 1
! if ( l_mat /= 0 ) then
!    ! Interior - interface 
!   do jobje = 1, trian%elems(ielem)%num_vefs
!      job_g = trian%elems(ielem)%vefs(jobje)
!      do jdof = fe_space%vef2dof%p(job_g), fe_space%vef2dof%p(job_g+1)-1
!         m_dof = fe_space%vef2dof%l(jdof,1)
!         m_mat = fe_space%vef2dof%l(jdof,3)     
!         do inode = fe_space%fe_array(ielem)%nodes_per_vef%p%p(iobje), &
!              & fe_space%fe_array(ielem)%nodes_per_vef%p%p(iobje+1)-1
!            l_node = fe_space%fe_array(ielem)%nodes_per_vef%p%l(inode)
!            l_dof = fe_space%fe_array(ielem)%elem2dof(l_node)
!            if ( .not. dof_graph%symmetric_storage ) then
!               dof_graph%ia(l_dof+1) = &
!                    & dof_graph%ia(l_dof+1) + 1 
!            else if ( m_dof >= l_dof ) then
!               dof_graph%ia(l_dof+1) = &
!                    & dof_graph%ia(l_dof+1) + 1
!            end if
!         end do
!      end do
!   end do
! end if
! end do
! end do
! end if

! end subroutine count_nnz_dofs_vol_vs_dofs_vefs_vol_by_continuity

! !*********************************************************************************
! ! List NNZ (number of nonzero entries) for DOFs on the interface (VEFs) of elements against
! ! both interior and interface nodes.
! !*********************************************************************************
! subroutine list_nnz_dofs_vol_vs_dofs_vefs_vol_by_continuity ( iblock, jblock, trian, fe_space, dof_graph, aux_ia )  
! implicit none
! ! Parameters
! integer(ip), intent(in)                     :: iblock, jblock
! type(triangulation_t), intent(in)         :: trian 
! type(SB_fe_space_t), intent(in)                 :: fe_space 
! type(graph_t), intent(inout)              :: dof_graph
! integer(ip), intent(inout)                  :: aux_ia(:) 

! ! Local variables
! integer(ip) :: ielem, inode, iobje, iprob, jdof, jnode, job_g
! integer(ip) :: jobje, l_dof, l_mat, l_node, m_dof, m_mat
! integer(ip) :: m_node, nvapbi, nvapbj, i, ic

! if (.not.fe_space%static_condensation) then   
! do ielem  = 1, trian%num_elems
! iobje = trian%elems(ielem)%num_vefs+1
! iprob = fe_space%fe_array(ielem)%problem 
! !l_var = g2l(ivars,iprob)
! ! Interior - interior (inside element)
! do inode = fe_space%fe_array(ielem)%nodes_per_vef%p%p(iobje), &
!  & fe_space%fe_array(ielem)%nodes_per_vef%p%p(iobje+1)-1
! l_node = fe_space%fe_array(ielem)%nodes_per_vef%p%l(inode)
! l_dof = fe_space%fe_array(ielem)%elem2dof(l_node)
! if ( .not. dof_graph%symmetric_storage ) then
!   do jnode = fe_space%fe_array(ielem)%nodes_per_vef%p%p(iobje), &
!        & fe_space%fe_array(ielem)%nodes_per_vef%p%p(iobje+1)-1
!      m_node = fe_space%fe_array(ielem)%nodes_per_vef%p%l(jnode)
!      m_dof = fe_space%fe_array(ielem)%elem2dof(m_node)
!      i= aux_ia(l_dof)
!      dof_graph%ja(i) = m_dof
!      aux_ia(l_dof) = aux_ia(l_dof)+1
!   end do
! else ! ltype == csr_symm 
!   do jnode = fe_space%fe_array(ielem)%nodes_per_vef%p%p(iobje), &
!        & fe_space%fe_array(ielem)%nodes_per_vef%p%p(iobje+1)-1
!      m_node = fe_space%fe_array(ielem)%nodes_per_vef%p%l(jnode)
!      m_dof = fe_space%fe_array(ielem)%elem2dof(m_node)
!      if ( m_dof >= l_dof ) then
!         ic = aux_ia(l_dof)
!         dof_graph%ja(ic) = m_dof
!         aux_ia(l_dof) = aux_ia(l_dof)+1
!      end if
!   end do
! end if
! end do
! end do
! if ( fe_space%fe_array(ielem)%continuity /= 0 ) then
!    ! Interior - border (inside element)
! do jobje = 1, trian%elems(ielem)%num_vefs
! job_g = trian%elems(ielem)%vefs(jobje)
! do jdof = fe_space%vef2dof%p(job_g), fe_space%vef2dof%p(job_g+1)-1
!   m_dof = fe_space%vef2dof%l(jdof,1)    
!   do inode = fe_space%fe_array(ielem)%nodes_per_vef%p%p(iobje), &
!        & fe_space%fe_array(ielem)%nodes_per_vef%p%p(iobje+1)-1
!      l_node = fe_space%fe_array(ielem)%nodes_per_vef%p%l(inode)
!      l_dof = fe_space%fe_array(ielem)%elem2dof(l_node)
!      if ( .not. dof_graph%symmetric_storage ) then
!         ic = aux_ia(l_dof)
!         dof_graph%ja(ic) = m_dof
!         aux_ia(l_dof) = aux_ia(l_dof)+1
!      else if ( m_dof >= l_dof ) then
!         ic = aux_ia(l_dof)
!         dof_graph%ja(ic) = m_dof
!         aux_ia(l_dof) = aux_ia(l_dof)+1
!      end if
!   end do
! end if
! end do
! end do
! end if
! end do
! end do
! end if
! end subroutine list_nnz_dofs_vol_vs_dofs_vefs_vol_by_continuity





end module SB_fe_space_names
