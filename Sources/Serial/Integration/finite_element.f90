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
module finite_element_names
  ! Modules
  use types_names
  use memor_names
  use array_names
  use integration_tools_names
  use interpolation_tools_names
  !use face_integration_names
  use fe_space_types_names
  !use dof_descriptor_names
  use migratory_element_names
  !use conditions_names

#ifdef memcheck
  use iso_c_binding
#endif
  implicit none
  private

  ! Information of each element of the FE space
  type, extends(migratory_element_t) :: finite_element_t
     ! Reference element info          
     type(reference_element_pointer_t), allocatable :: reference_element_vars(:)    ! Topology of the reference finite element
     type(reference_element_t), pointer :: p_geo_reference_element => NULL()    ! Topology of the reference geometry ( idem fe w/ p=1)
     integer(ip),      allocatable :: order(:)                ! Order per variable
     type(volume_integrator_pointer_t), allocatable :: integ(:) ! Pointer to integration parameters
     type(interpolator_pointer_t)     , allocatable :: inter(:) ! Pointer to interpolator
     ! order in f_inf, it can be eliminated

     ! Problem and approximation
     integer(ip)                   :: problem           ! Problem to be solved
     integer(ip)                   :: num_vars          ! Number of variables of the problem
     integer(ip)                   :: approximation     ! Discretization to be used

     ! Connectivity
     integer(ip)       , allocatable :: continuity(:)     ! Continuity flag per variable
     logical           , allocatable :: enable_face_integration(:)  ! Face integration flag per variable
     type(list_pointer_t), allocatable :: nodes_per_vef(:)   ! Nodes per vefq (including interior) (nvars)
     integer(ip)                     :: material          ! Material ! SB.alert : material can be used as p   
     ! use of material still unclear

     ! Local to global 
     integer(ip)     , allocatable   :: elem2dof(:,:)   ! Map from elem to dof
     
     ! Unknown + other values
     real(rp)        , allocatable :: unkno(:,:,:)      ! Values of the solution on the nodes of the elem  (max_num_nodes, nvars, time_steps_to_store)
     real(rp)        , allocatable :: gauss_properties(:,:,:) ! Gauss point level properties with history, e.g. subscales,  rank?

     ! Boundary conditions
     integer(ip), allocatable :: bc_code(:,:)   ! Boundary Condition values
     
     ! Auxiliary working arrays (element matrix and vector)
     type(array_rp2_t), pointer :: p_mat ! Pointer to the elemental matrix
     type(array_rp1_t), pointer :: p_vec ! Pointer to the elemental vector_t

   contains
     procedure :: size   => finite_element_size
     procedure :: pack   => finite_element_pack
     procedure :: unpack => finite_element_unpack
  end type finite_element_t

  ! Information relative to the faces
  type fe_face_t
     
     ! Reference face info
     type(element_face_integrator_t) :: integ(2)  ! Pointer to face integration

     ! Face mesh info
     integer(ip)               :: face_vef
     integer(ip)               :: neighbor_element(2) ! Neighbor elements
     integer(ip)               :: local_face(2)       ! Face pos in element

     ! Auxiliary working arrays (face+element matrix and vector)
     type(array_rp2_t), pointer  :: p_mat ! Pointer to the elemental matrix
     type(array_rp1_t), pointer  :: p_vec ! Pointer to face integration vector_t

     !type(array_ip1_t), allocatable:: o2n(2)           ! permutation of the gauss points in elem2
     ! SB.alert : temporary, it is a lot of memory, and should be handled via a hash table
  end type fe_face_t

  ! Types
  public :: finite_element_t, fe_face_t

  ! Methods
  public :: finite_element_print, finite_element_free_unpacked, impose_strong_dirichlet_data

contains
  subroutine finite_element_print ( lunou, finite_element )
    implicit none
    integer(ip)      , intent(in) :: lunou
    type(finite_element_t), intent(in) :: finite_element

    integer(ip) :: ivar

    write (lunou, '(a)')     '*** begin finite element data structure ***'
    write (lunou, '(a,i10)') 'Problem: ', finite_element%problem
    write (lunou,*) 'Element dofs: ', finite_element%elem2dof

    write (lunou,*) 'Number of unknowns: ', finite_element%num_vars
    write (lunou,*) 'Order of each variable: ', finite_element%order
    write (lunou,*) 'Continuity of each variable: ', size(finite_element%continuity)
    write (lunou,*) 'Continuity of each variable: ', finite_element%continuity
    write (lunou,*) 'Enable face integration for each variable: ', finite_element%enable_face_integration
    write (lunou,*) 'Element material: ', finite_element%material
    write (lunou,*) 'Boundary conditions code: ', finite_element%bc_code

    write (lunou,*) 'Fixed info of each interpolation: '
    do ivar=1, finite_element%num_vars
       write (lunou, *) 'Type: ', finite_element%reference_element_vars(ivar)%p%ftype
       write (lunou, *) 'Order: ', finite_element%reference_element_vars(ivar)%p%order
       write (lunou, *) 'Nobje: ', finite_element%reference_element_vars(ivar)%p%nvef
       write (lunou, *) 'Nnode: ', finite_element%reference_element_vars(ivar)%p%nnode
       write (lunou, *) 'Nobje_dim: ', finite_element%reference_element_vars(ivar)%p%nvef_dim
       write (lunou, *) 'Nodes_vef: ', finite_element%reference_element_vars(ivar)%p%nodes_vef
       write (lunou, *) 'ndxob%p:  ', finite_element%reference_element_vars(ivar)%p%ndxob%p
       write (lunou, *) 'ndxob%l:  ', finite_element%reference_element_vars(ivar)%p%ndxob%l
       write (lunou, *) 'ntxob%p:  ', finite_element%reference_element_vars(ivar)%p%ntxob%p
       write (lunou, *) 'ntxob%l:  ', finite_element%reference_element_vars(ivar)%p%ntxob%l
       write (lunou, *) 'crxob%p:  ', finite_element%reference_element_vars(ivar)%p%crxob%p
       write (lunou, *) 'crxob%l:  ', finite_element%reference_element_vars(ivar)%p%crxob%l
    end do

    write (lunou,*) 'Unknown values: ', finite_element%unkno

  end subroutine finite_element_print

   ! SB.alert : to be thought now

  subroutine finite_element_size (my, n)
    implicit none
    class(finite_element_t), intent(in)  :: my
    integer(ip)            , intent(out) :: n
    
    ! Locals
    integer(ieep) :: mold(1)
    integer(ip)   :: size_of_ip
    
    size_of_ip   = size(transfer(1_ip ,mold))

    n = size_of_ip*3 + 2*size_of_ip*(my%num_vars)

  end subroutine finite_element_size

  subroutine finite_element_pack (my, n, buffer)
    implicit none
    class(finite_element_t), intent(in)  :: my
    integer(ip)            , intent(in)   :: n
    integer(ieep)            , intent(out)  :: buffer(n)
    
    ! Locals
    integer(ieep) :: mold(1)
    integer(ip) :: size_of_ip

    integer(ip) :: start, end

    size_of_ip   = size(transfer(1_ip ,mold))

    start = 1
    end   = start + size_of_ip -1
    buffer(start:end) = transfer(my%num_vars,mold)

    start = end + 1
    end   = start + size_of_ip - 1
    buffer(start:end) = transfer(my%problem,mold)

    start = end + 1
    end   = start + size_of_ip - 1
    buffer(start:end) = transfer(my%material,mold)

    start = end + 1
    end   = start + my%num_vars*size_of_ip - 1
    buffer(start:end) = transfer(my%order,mold)

    start = end + 1
    end   = start + my%num_vars*size_of_ip - 1
    buffer(start:end) = transfer(my%continuity,mold)

    start = end + 1
    end   = start + my%num_vars*size_of_ip - 1
    buffer(start:end) = transfer(my%enable_face_integration,mold)

  end subroutine finite_element_pack

  subroutine finite_element_unpack(my, n, buffer)
    implicit none
    class(finite_element_t), intent(inout) :: my
    integer(ip)            , intent(in)     :: n
    integer(ieep)            , intent(in)     :: buffer(n)

    ! Locals
    integer(ieep) :: mold(1)
    integer(ip) :: size_of_ip
    integer(ip) :: start, end
    
    size_of_ip   = size(transfer(1_ip ,mold))

    start = 1
    end   = start + size_of_ip -1
    my%num_vars  = transfer(buffer(start:end), my%num_vars)

    start = end + 1
    end   = start + size_of_ip - 1
    my%problem  = transfer(buffer(start:end), my%problem)

    start = end + 1
    end   = start + size_of_ip - 1
    my%material  = transfer(buffer(start:end), my%material)

    call memalloc( my%num_vars, my%order, __FILE__, __LINE__ )

    start = end + 1
    end   = start + my%num_vars*size_of_ip - 1
    my%order = transfer(buffer(start:end), my%order)
    
    call memalloc( my%num_vars, my%continuity, __FILE__, __LINE__ )
     
    start = end + 1
    end   = start + my%num_vars*size_of_ip - 1
    my%continuity = transfer(buffer(start:end), my%continuity)
    
    call memalloc( my%num_vars, my%enable_face_integration, __FILE__, __LINE__ )
     
    start = end + 1
    end   = start + my%num_vars*size_of_ip - 1
    my%enable_face_integration = transfer(buffer(start:end), my%enable_face_integration)
    
  end subroutine finite_element_unpack

  subroutine finite_element_free_unpacked(finite_element)
    implicit none
    type(finite_element_t), intent(inout) :: finite_element

    call memfree( finite_element%order, __FILE__, __LINE__ )
    call memfree( finite_element%continuity, __FILE__, __LINE__ )
    call memfree( finite_element%enable_face_integration, __FILE__, __LINE__ )
    
  end subroutine finite_element_free_unpacked

 !=============================================================================
  subroutine impose_strong_dirichlet_data (finite_element) 
    implicit none
    ! Parameters
    type(finite_element_t)    , intent(inout)  :: finite_element

    ! Locals
    integer(ip) :: iprob, count, ivars, inode, idof
    
    iprob = finite_element%problem
    count = 0

    !write (*,*) 'start assembly bc of matrix : ', finite_element%p_mat%a
    do ivars = 1, finite_element%num_vars
       do inode = 1,finite_element%reference_element_vars(ivars)%p%nnode
          count = count + 1
          idof = finite_element%elem2dof(inode,ivars)
          if ( idof  == 0 ) then
             finite_element%p_vec%a(:) = finite_element%p_vec%a(:) - finite_element%p_mat%a(:,count)*finite_element%unkno(inode,ivars,1)
             !write (*,*) 'add to vector', -finite_element%p_mat%a(:,count)*finite_element%unkno(inode,ivars,1)
          end if
       end do
    end do

    !write(*,*) 'finite_elementvec :', finite_element%p_vec%a

  end subroutine impose_strong_dirichlet_data

end module finite_element_names
