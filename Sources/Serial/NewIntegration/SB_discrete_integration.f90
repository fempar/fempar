module SB_discrete_integration_names
  use shape_values_names
  use reference_fe_names
  use types_names

  implicit none
# include "debug.i90"
  private

  type, abstract :: SB_discrete_integration_t
     !integer(ip)              :: domain_dimension      ! either 2 or 3
     !integer(ip), pointer     :: domain(:) => NULL()   ! List of entities over which integration is performed
   contains
     procedure (integrate_interface), deferred :: integrate
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
     subroutine integrate_interface (this,vol_int,elmat,elvec)
       import :: SB_discrete_integration_t, SB_volume_integrator_t, rp
       implicit none
       class(SB_discrete_integration_t), intent(inout) :: this
       type(SB_volume_integrator_t)    , intent(inout) :: vol_int
       real(rp), intent(inout) :: elmat(:,:),elvec(:)
     end subroutine integrate_interface
  end interface

contains

end module SB_discrete_integration_names

!%%%%%%%%%%%%%%%%%%%%%%%

module poisson_discrete_integration_names
  use shape_values_names
  use SB_discrete_integration_names
  use reference_fe_names
  use types_names
  implicit none
# include "debug.i90"
  private
  type, extends(SB_discrete_integration_t) :: poisson_discrete_integration_t
  integer(ip) :: viscosity
contains
  procedure :: integrate
end type poisson_discrete_integration_t

public :: poisson_discrete_integration_t

contains

  subroutine integrate (this, vol_int, elmat, elvec)
    implicit none
    class(poisson_discrete_integration_t), intent(inout) :: this
    type(SB_volume_integrator_t)    , intent(inout) :: vol_int
    real(rp), intent(inout) :: elmat(:,:), elvec(:)

    integer(ip)  :: ndime
    real(rp)     :: dtinv, c1, c2, mu

    integer(ip)  :: idime,igaus,inode,jnode,ngaus,nnode
    real(rp) :: factor

    type(SB_quadrature_t), pointer :: quad
    type(SB_interpolation_t), pointer :: int
    class(reference_fe_t), pointer :: ref_fe
    type(fe_map_t), pointer :: fe_map

    type(shape_values_t), pointer :: shape_gradient_test, shape_gradient_trial

    real(rp), pointer :: grad_test(:,:), grad_trial(:,:)

    integer(ip) :: a, b

    quad => vol_int%get_quadrature()
    int => vol_int%get_interpolation()
    ref_fe => vol_int%get_reference_fe()
    fe_map => vol_int%get_fe_map()

    elmat = 0.0_rp
    elvec = 0.0_rp

    nnode = ref_fe%get_number_nodes()
    ngaus = quad%get_number_integration_points()
    ndime = ref_fe%get_number_dimensions()

    shape_gradient_test => vol_int%get_gradients()
    shape_gradient_trial => vol_int%get_gradients()

    do igaus = 1,ngaus
       factor = fe_map%get_det_jacobian(igaus) * quad%get_weight(igaus)
       do inode = 1, nnode
          grad_trial => shape_gradient_trial%get_value(inode,igaus)
          do jnode = 1, nnode
             grad_test => shape_gradient_test%get_value(jnode,igaus)

             do a = 1,size(shape_gradient_test%get_value(inode,igaus),1)
                do b= 1,size(shape_gradient_test%get_value(inode,igaus),2)
                   elmat(inode,jnode) = elmat(inode,jnode) &
                        + factor * grad_test(a,b) * grad_trial(a,b) 
                end do
             end do

          end do
       end do
    end do

 write (*,*) 'XXXXXXXXX ELMAT XXXXXXXXX'
 write (*,*) elmat


 ! Apply boundary conditions
 !call impose_strong_dirichlet_data(fe) 

end subroutine integrate


end module poisson_discrete_integration_names
