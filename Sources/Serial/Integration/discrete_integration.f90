module discrete_integration_names
  use field_names
  use reference_fe_names
  use types_names
  use assembler_names
  use serial_fe_space_names
  use memor_names

  implicit none
# include "debug.i90"
  private

  type, abstract :: discrete_integration_t
   contains
     procedure (integrate_interface), deferred :: integrate
     procedure                                 :: impose_strong_dirichlet_data
  end type discrete_integration_t

  type p_discrete_integration_t
     class(discrete_integration_t), pointer :: p => NULL()  
  end type p_discrete_integration_t

  public :: discrete_integration_t, p_discrete_integration_t

  abstract interface
     subroutine integrate_interface ( this, fe_space, assembler  )
       import :: discrete_integration_t, serial_fe_space_t, assembler_t
       implicit none
       class(discrete_integration_t), intent(in)    :: this
       class(serial_fe_space_t)     , intent(inout) :: fe_space
       class(assembler_t)           , intent(inout) :: assembler
     end subroutine integrate_interface
  end interface

contains

subroutine impose_strong_dirichlet_data ( this, elmat, elvec, code, value, nnode, num_fe_spaces )
 implicit none
 class(discrete_integration_t) :: this
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

end module discrete_integration_names

module poisson_discrete_integration_names
  use field_names
  use assembler_names
  use serial_fe_space_names
  use discrete_integration_names
  use reference_fe_names
  use types_names
  use memor_names
  
  implicit none
# include "debug.i90"
  private
  type, extends(discrete_integration_t) :: poisson_discrete_integration_t
     integer(ip) :: viscosity 
   contains
     procedure :: integrate
  end type poisson_discrete_integration_t
  
  public :: poisson_discrete_integration_t
  
contains
  
  subroutine integrate ( this, fe_space, assembler )
    implicit none
    class(poisson_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)          , intent(inout) :: fe_space
    class(assembler_t)                , intent(inout) :: assembler

    type(finite_element_t), pointer :: fe
    type(volume_integrator_t), pointer :: vol_int
    real(rp), allocatable :: elmat(:,:), elvec(:)
    type(fe_map_t), pointer :: fe_map
    type(quadrature_t), pointer :: quad

    integer(ip)  :: igaus,inode,jnode,ngaus
    real(rp)     :: factor

    type(vector_field_t) :: grad_test, grad_trial

    integer(ip) :: number_fe_spaces

    integer(ip), pointer :: field_blocks(:)
    logical, pointer :: field_coupling(:,:)

    integer(ip) :: ielem, iapprox, number_nodes
    type(i1p_t), pointer :: elem2dof(:)
    type(i1p_t), pointer :: bc_code(:)
    type(r1p_t), pointer :: bc_value(:)
    integer(ip), allocatable :: number_nodes_per_field(:)

    
    number_fe_spaces = fe_space%get_number_fe_spaces()
    field_blocks => fe_space%get_field_blocks()
    field_coupling => fe_space%get_field_coupling()

    fe => fe_space%get_finite_element(1)
    number_nodes = fe%get_number_nodes()
    call memalloc ( number_nodes, number_nodes, elmat, __FILE__, __LINE__ )
    call memalloc ( number_nodes, elvec, __FILE__, __LINE__ )
    call memalloc ( number_fe_spaces, number_nodes_per_field, __FILE__, __LINE__ )
    call fe%get_number_nodes_per_field( number_nodes_per_field )

    call fe_space%initialize_integration()
    
    quad => fe%get_quadrature()
    ngaus = quad%get_number_evaluation_points()
    do ielem = 1, fe_space%get_number_elements()
       elmat = 0.0_rp
       elvec = 0.0_rp

       fe => fe_space%get_finite_element(ielem)
       call fe%update_integration()
       
       fe_map   => fe%get_fe_map()
       vol_int  => fe%get_volume_integrator(1)
       elem2dof => fe%get_elem2dof()
       bc_code  => fe%get_bc_code()
       bc_value => fe%get_bc_value()

       do igaus = 1,ngaus
          factor = fe_map%get_det_jacobian(igaus) * quad%get_weight(igaus)
          do inode = 1, number_nodes
             call vol_int%get_gradient(inode,igaus,grad_trial)
             do jnode = 1, number_nodes
                call vol_int%get_gradient(jnode,igaus,grad_test)
                elmat(inode,jnode) = elmat(inode,jnode) + factor * grad_test * grad_trial
             end do
          end do
       end do
       !write (*,*) 'XXXXXXXXX ELMAT XXXXXXXXX'
       !write (*,*) elmat
       
       ! Apply boundary conditions
       call this%impose_strong_dirichlet_data( elmat, elvec, bc_code, bc_value, number_nodes_per_field, number_fe_spaces )
       call assembler%assembly( number_fe_spaces, number_nodes_per_field, elem2dof, field_blocks,  field_coupling, elmat, elvec )
    end do
    call memfree ( number_nodes_per_field, __FILE__, __LINE__ )
    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )
  end subroutine integrate
end module poisson_discrete_integration_names

module vector_laplacian_discrete_integration_names
use field_names
use assembler_names
use serial_fe_space_names
use discrete_integration_names
use reference_fe_names
use types_names
use memor_names

implicit none
# include "debug.i90"
private
type, extends(discrete_integration_t) :: vector_laplacian_discrete_integration_t
integer(ip) :: viscosity 
contains
procedure :: integrate
end type vector_laplacian_discrete_integration_t

public :: vector_laplacian_discrete_integration_t

contains
  subroutine integrate ( this, fe_space, assembler )
    implicit none
    class(vector_laplacian_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)                   , intent(inout) :: fe_space
    class(assembler_t)                         , intent(inout) :: assembler

    type(finite_element_t), pointer :: fe
    type(volume_integrator_t), pointer :: vol_int_first_fe, vol_int_second_fe
    real(rp), allocatable :: elmat(:,:), elvec(:)
    type(fe_map_t), pointer :: fe_map
    type(quadrature_t), pointer :: quad
    integer(ip), allocatable :: number_nodes_per_field(:)

    integer(ip)  :: igaus,inode,jnode,ioffset,joffset,ngaus
    real(rp) :: factor

    type(vector_field_t) :: grad_test_scalar, grad_trial_scalar
    type(tensor_field_t) :: grad_test_vector, grad_trial_vector
    
    integer(ip) :: i, number_fe_spaces

    integer(ip), pointer :: field_blocks(:)
    logical, pointer :: field_coupling(:,:)

    integer(ip) :: ielem, number_nodes
    type(i1p_t), pointer :: elem2dof(:)
    type(i1p_t), pointer :: bc_code(:)
    type(r1p_t), pointer :: bc_value(:)

    number_fe_spaces = fe_space%get_number_fe_spaces()
    field_blocks => fe_space%get_field_blocks()
    field_coupling => fe_space%get_field_coupling()

    fe => fe_space%get_finite_element(1)
    number_nodes = fe%get_number_nodes()
    call memalloc ( number_nodes, number_nodes, elmat, __FILE__, __LINE__ )
    call memalloc ( number_nodes, elvec, __FILE__, __LINE__ )
    call memalloc ( number_fe_spaces, number_nodes_per_field, __FILE__, __LINE__ )
    call fe%get_number_nodes_per_field( number_nodes_per_field )
    
    call fe_space%initialize_integration()
    
    quad  => fe%get_quadrature()
    ngaus = quad%get_number_evaluation_points()
    do ielem = 1, fe_space%get_number_elements()
       elmat = 0.0_rp
       elvec = 0.0_rp

       fe => fe_space%get_finite_element(ielem)
       call fe%update_integration()
       
       fe_map            => fe%get_fe_map()
       vol_int_first_fe  => fe%get_volume_integrator(1)
       vol_int_second_fe => fe%get_volume_integrator(2)
       elem2dof          => fe%get_elem2dof()
       bc_code           => fe%get_bc_code()
       bc_value          => fe%get_bc_value()

       do igaus = 1,ngaus
          factor = fe_map%get_det_jacobian(igaus) * quad%get_weight(igaus)
          do inode = 1, number_nodes_per_field(1)
             call vol_int_first_fe%get_gradient(inode,igaus,grad_trial_scalar)
             !call vol_int_first_fe%get_gradient(inode,igaus,grad_trial_vector)
             do jnode = 1, number_nodes_per_field(1)
                call vol_int_first_fe%get_gradient(jnode,igaus,grad_test_scalar)
                !call vol_int_first_fe%get_gradient(jnode,igaus,grad_test_vector)
                elmat(inode,jnode) = elmat(inode,jnode) + factor * grad_test_scalar * grad_trial_scalar
                !elmat(inode,jnode) = elmat(inode,jnode) + factor * (grad_test_vector .doublecontract. grad_trial_vector)
             end do
          end do

          do inode = 1, number_nodes_per_field(2)
             ioffset = number_nodes_per_field(1)+inode
             !call vol_int_second_fe%get_gradient(inode,igaus,grad_trial_scalar)
             call vol_int_second_fe%get_gradient(inode,igaus,grad_trial_vector)
             ! write(*,*) inode, grad_trial_vector%value
             do jnode = 1, number_nodes_per_field(2)
                joffset = number_nodes_per_field(1)+jnode
                call vol_int_second_fe%get_gradient(jnode,igaus,grad_test_vector)
                !call vol_int_second_fe%get_gradient(jnode,igaus,grad_test_scalar)
               elmat(ioffset,joffset) = elmat(ioffset,joffset) + factor * double_contract(grad_test_vector,grad_trial_vector)
                !elmat(ioffset,joffset) = elmat(ioffset,joffset) + factor * grad_test_scalar * grad_trial_scalar
             end do
          end do
       end do
       
       !write (*,*) 'XXXXXXXXX ELMAT field 1 XXXXXXXXX'
       !write (*,*) elmat(1:number_nodes_per_field(1),1:number_nodes_per_field(1))
       
       !write (*,*) 'XXXXXXXXX ELMAT field 2 XXXXXXXXX'
       !write (*,*) elmat(number_nodes_per_field(1)+1:,number_nodes_per_field(1)+1:)
       
       call this%impose_strong_dirichlet_data( elmat, elvec, bc_code, bc_value, number_nodes_per_field, number_fe_spaces )
       call assembler%assembly( number_fe_spaces, number_nodes_per_field, elem2dof, field_blocks,  field_coupling, elmat, elvec )      
    end do
    call memfree ( number_nodes_per_field, __FILE__, __LINE__ )
    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )
  end subroutine integrate

end module vector_laplacian_discrete_integration_names

#include "sbm_cdr_discrete_integration.i90"
#include "sbm_maxwell_discrete_integration.i90"
#include "sbm_stokes_discrete_integration.i90"
