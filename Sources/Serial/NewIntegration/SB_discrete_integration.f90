module SB_discrete_integration_names
  use shape_values_names
  use reference_fe_names
  use types_names
  use SB_assembler_names
  use SB_fe_space_names
  use memor_names

  implicit none
# include "debug.i90"
  private

  type, abstract :: SB_discrete_integration_t
     !integer(ip)              :: domain_dimension      ! either 2 or 3
     !integer(ip), pointer     :: domain(:) => NULL()   ! List of entities over which integration is performed
   contains
     procedure (integrate_interface), deferred :: integrate
     procedure :: impose_strong_dirichlet_data
     !procedure(create_integration_interface) , deferred :: create 
     !procedure(compute_integration_interface), deferred :: compute
     !procedure(free_integration_interface)   , deferred :: free 
     !procedure :: compute_face
  end type SB_discrete_integration_t

  type SB_p_discrete_integration_t
     class(SB_discrete_integration_t), pointer :: p => NULL()  
  end type SB_p_discrete_integration_t

  public :: SB_discrete_integration_t, SB_p_discrete_integration_t

  abstract interface
     subroutine integrate_interface ( this, fe_space, assembler  )
       import :: SB_discrete_integration_t, SB_serial_fe_space_t, SB_assembler_t
       implicit none
       class(SB_discrete_integration_t) :: this
       class(SB_serial_fe_space_t), intent(inout) :: fe_space
       class(SB_assembler_t) :: assembler
     end subroutine integrate_interface
  end interface

contains

subroutine impose_strong_dirichlet_data ( this, elmat, elvec, code, value, nnode, num_fe_spaces )
 implicit none
 class(SB_discrete_integration_t) :: this
 real(rp), intent(in) :: elmat(:,:)
 real(rp), intent(inout) :: elvec(:)  
 type(i1p_t), intent(in) :: code(:)
 type(r1p_t), intent(in) :: value(:)
 integer(ip), intent(in) :: nnode(:), num_fe_spaces
 integer(ip) :: i, c, ifes

 c = 0
 do ifes = 1, num_fe_spaces
    do i = 1, nnode(ifes)
       c = c + 1
       if ( code(ifes)%p(i) /= 0 ) then
          elvec = elvec - elmat(:,c)*value(ifes)%p(i)
       end if
    end do
 end do

end subroutine impose_strong_dirichlet_data

end module SB_discrete_integration_names

module poisson_discrete_integration_names
  use shape_values_names
  use SB_assembler_names
  use SB_fe_space_names
  use SB_discrete_integration_names
  use reference_fe_names
  use types_names
  use memor_names
  
  implicit none
# include "debug.i90"
  private
  type, extends(SB_discrete_integration_t) :: poisson_discrete_integration_t
     integer(ip) :: viscosity 
     !integer(ip), parameter :: u=1, p=2
   contains
     procedure :: integrate
  end type poisson_discrete_integration_t
  
  public :: poisson_discrete_integration_t
  
contains
  
  subroutine integrate ( this, fe_space, assembler )
    implicit none
    class(poisson_discrete_integration_t) :: this
    class(SB_serial_fe_space_t), intent(inout) :: fe_space
    class(SB_assembler_t) :: assembler

    class(SB_finite_element_t), pointer :: fe
    type(SB_p_volume_integrator_t), pointer :: vol_int(:)
    real(rp), allocatable :: elmat(:,:), elvec(:)
    type(fe_map_t), pointer :: fe_map
    type(SB_quadrature_t), pointer :: quad
    integer(ip), allocatable :: number_nodes(:)

    integer(ip)  :: ndime, num_elems
    real(rp)     :: dtinv, c1, c2, mu

    integer(ip)  :: idime,igaus,inode,jnode,ngaus
    real(rp) :: factor

    !type(SB_quadrature_t), pointer :: quad

    type(shape_values_t), pointer :: shape_gradient_test_u, shape_gradient_trial_u

    real(rp), pointer :: grad_test(:,:), grad_trial(:,:)

    integer(ip) :: a, b
    integer(ip) :: i, number_blocks, number_fe_spaces

    integer(ip), pointer :: blocks(:)
    logical, pointer :: blocks_coupling(:,:)

    integer(ip) :: ielem, iapprox, nnodes
    type(i1p_t), pointer :: elem2dof(:)
    type(i1p_t), pointer :: bc_code(:)
    type(r1p_t), pointer :: bc_value(:)



    number_blocks = fe_space%get_number_blocks()
    number_fe_spaces = fe_space%get_number_fe_spaces()
    blocks => fe_space%get_field_blocks()
    blocks_coupling => fe_space%get_field_coupling()

    fe => fe_space%get_finite_element(1)
    nnodes = fe%get_number_nodes()
    call memalloc ( nnodes, nnodes, elmat, __FILE__, __LINE__ )
    call memalloc ( nnodes, elvec, __FILE__, __LINE__ )

    allocate( number_nodes(number_fe_spaces) )

    call fe_space%initialize_integration()

    !write(*,*) 'number_blocks,number_fe_spaces,blocks,blocks_coupling,nnodes',number_blocks,number_fe_spaces,blocks,blocks_coupling,nnodes

    !do iapprox = 1, size(this%approximations)
    do ielem = 1, fe_space%get_number_elements()
       elmat = 0.0_rp
       elvec = 0.0_rp

       fe => fe_space%get_finite_element(ielem)
       call fe%update_integration( number_fe_spaces )
       quad => fe%get_quadrature()
       fe_map => fe%get_fe_map()
       vol_int => fe%get_volume_integrator()

       !call fe_map%print()
       !call quad%print()
       !call vol_int(1)%p%print()
       !check(0==1)
       !gcall vol_int%print()

       call fe%get_number_nodes_field( number_nodes, number_fe_spaces )
       elem2dof => fe%get_elem2dof()
       bc_code =>  fe%get_bc_code()
       bc_value => fe%get_bc_value()
       call fe%get_number_nodes_field(number_nodes,number_fe_spaces)
       !write(*,*) 'elem2dof',elem2dof
       !write(*,*) 'bc_code',bc_code
       !write(*,*) 'bc_value',bc_value
       !write(*,*) 'number_nodes',number_nodes

       ngaus = quad%get_number_evaluation_points()
       !write(*,*) 'ngaus',ngaus

       shape_gradient_test_u => vol_int(1)%p%get_gradients()
       shape_gradient_trial_u => vol_int(1)%p%get_gradients()
       !call shape_gradient_test_u%print()
       do igaus = 1,ngaus
          factor = fe_map%get_det_jacobian(igaus) * quad%get_weight(igaus)
          do inode = 1, number_nodes(1)
             grad_trial => shape_gradient_test_u%get_value(inode,igaus)
             do jnode = 1, number_nodes(1)
                grad_test => shape_gradient_trial_u%get_value(jnode,igaus)
                ! elmat = elmat + grad_trial*grad_test
                do a = 1,size(shape_gradient_test_u%get_value(inode,igaus),1)
                   do b= 1,size(shape_gradient_trial_u%get_value(inode,igaus),2)
                      ! write(*,*) 'BLOCK 1',factor * grad_test(a,b) * grad_trial(a,b)
                      elmat(inode,jnode) = elmat(inode,jnode) &
                           + factor * grad_test(a,b) * grad_trial(a,b) 
                   end do
                end do
             end do
          end do
       end do
       !write (*,*) 'XXXXXXXXX ELMAT XXXXXXXXX'
       !write (*,*) elmat

       call this%impose_strong_dirichlet_data( elmat, elvec, bc_code, bc_value, number_nodes, number_fe_spaces )

       call assembler%assembly( number_fe_spaces, elem2dof, blocks, number_nodes, blocks_coupling, elmat, elvec )
       !end do

       !write (*,*) 'XXXXXXXXX ELMAT XXXXXXXXX'
       !write (*,*) elmat

       !do inode = 1, number_nodes(1)
       !   do jnode = 13,16! 1, number_nodes(1)
       !      aux(inode)  = aux(inode) + elmat(inode,jnode)
       !   end do
       !end do

       !write (*,*) 'XXXXXXXXX CHECK 0 XXXXXXXXX'
       !write (*,*) aux

       !write (*,*) 'XXXXXXXXX ELVEC XXXXXXXXX'
       !write (*,*) elvec

    end do



    ! Apply boundary conditions
    !call impose_strong_dirichlet_data(fe) 

    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )

  end subroutine integrate





end module poisson_discrete_integration_names











module vector_laplacian_discrete_integration_names
use shape_values_names
use SB_assembler_names
use SB_fe_space_names
use SB_discrete_integration_names
use reference_fe_names
use types_names
use memor_names

implicit none
# include "debug.i90"
private
type, extends(SB_discrete_integration_t) :: vector_laplacian_discrete_integration_t
integer(ip) :: viscosity 
!integer(ip), parameter :: u=1, p=2
contains
procedure :: integrate
end type vector_laplacian_discrete_integration_t

public :: vector_laplacian_discrete_integration_t

contains

  subroutine integrate ( this, fe_space, assembler )
    implicit none
    class(vector_laplacian_discrete_integration_t) :: this
    class(SB_serial_fe_space_t), intent(inout) :: fe_space
    class(SB_assembler_t) :: assembler

    class(SB_finite_element_t), pointer :: fe
    type(SB_p_volume_integrator_t), pointer :: vol_int(:)
    real(rp), allocatable :: elmat(:,:), elvec(:)
    type(fe_map_t), pointer :: fe_map
    type(SB_quadrature_t), pointer :: quad
    integer(ip), allocatable :: number_nodes(:)

    integer(ip)  :: ndime, num_elems
    real(rp)     :: dtinv, c1, c2, mu

    integer(ip)  :: idime,igaus,inode,jnode,ngaus
    real(rp) :: factor

    !type(SB_quadrature_t), pointer :: quad

    type(shape_values_t), pointer :: shape_gradient_test_u, shape_gradient_trial_u
    type(shape_values_t), pointer :: shape_gradient_test_v, shape_gradient_trial_v

    real(rp), pointer :: grad_test(:,:), grad_trial(:,:)

    integer(ip) :: a, b
    integer(ip) :: i, number_blocks, number_fe_spaces

    integer(ip), pointer :: field_blocks(:)
    logical, pointer :: field_coupling(:,:)

    integer(ip) :: ielem, iapprox, nnodes
    type(i1p_t), pointer :: elem2dof(:)
    type(i1p_t), pointer :: bc_code(:)
    type(r1p_t), pointer :: bc_value(:)

    number_blocks = fe_space%get_number_blocks()
    number_fe_spaces = fe_space%get_number_fe_spaces()
    field_blocks => fe_space%get_field_blocks()
    field_coupling => fe_space%get_field_coupling()

    fe => fe_space%get_finite_element(1)
    nnodes = fe%get_number_nodes()
    call memalloc ( nnodes, nnodes, elmat, __FILE__, __LINE__ )
    call memalloc ( nnodes, elvec, __FILE__, __LINE__ )

    allocate( number_nodes(number_fe_spaces) )

    call fe_space%initialize_integration()

    do ielem = 1, fe_space%get_number_elements()

       elmat = 0.0_rp
       elvec = 0.0_rp

       fe => fe_space%get_finite_element(ielem)
       call fe%update_integration( number_fe_spaces )
       fe_map => fe%get_fe_map()
       vol_int => fe%get_volume_integrator()
       call fe%get_number_nodes_field( number_nodes, number_fe_spaces )
       elem2dof => fe%get_elem2dof()
       bc_code =>  fe%get_bc_code()
       bc_value => fe%get_bc_value()
       call fe%get_number_nodes_field(number_nodes,number_fe_spaces)
       quad => fe%get_quadrature()

       ngaus = quad%get_number_evaluation_points()

       shape_gradient_test_u => vol_int(1)%p%get_gradients()
       shape_gradient_trial_u => vol_int(1)%p%get_gradients()
       shape_gradient_test_v => vol_int(2)%p%get_gradients()
       shape_gradient_trial_v => vol_int(2)%p%get_gradients()

       do igaus = 1,ngaus
          factor = fe_map%get_det_jacobian(igaus) * quad%get_weight(igaus)
          do inode = 1, number_nodes(1)
             grad_trial => shape_gradient_test_u%get_value(inode,igaus)
             do jnode = 1, number_nodes(1)
                grad_test => shape_gradient_trial_u%get_value(jnode,igaus)
                ! elmat = elmat + grad_trial*grad_test
                do a = 1,size(shape_gradient_test_u%get_value(inode,igaus),1)
                   do b= 1,size(shape_gradient_trial_u%get_value(inode,igaus),2)
                      elmat(inode,jnode) = elmat(inode,jnode) &
                           + factor * grad_test(a,b) * grad_trial(a,b) 
                   end do
                end do
             end do
          end do

          do inode = 1, number_nodes(2)
             grad_trial => shape_gradient_test_v%get_value(inode,igaus)
             do jnode = 1, number_nodes(2)
                grad_test => shape_gradient_trial_v%get_value(jnode,igaus)
                ! elmat = elmat + grad_trial*grad_test
                do a = 1,size(shape_gradient_test_v%get_value(inode,igaus),1)
                   do b= 1,size(shape_gradient_trial_v%get_value(inode,igaus),2)
                      elmat(number_nodes(1)+inode,number_nodes(1)+jnode) = elmat(number_nodes(1)+inode,number_nodes(1)+jnode) &
                           + factor * grad_test(a,b) * grad_trial(a,b) 
                   end do
                end do
             end do
          end do
       end do

       call this%impose_strong_dirichlet_data( elmat, elvec, bc_code, bc_value, number_nodes, number_fe_spaces )

       call assembler%assembly( number_fe_spaces, elem2dof, field_blocks, number_nodes, field_coupling, elmat, elvec )      
    end do

    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )

  end subroutine integrate


end module vector_laplacian_discrete_integration_names
