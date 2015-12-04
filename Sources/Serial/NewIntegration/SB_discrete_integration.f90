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
     subroutine integrate_interface (this, vol_int, elmat, elvec, quad, fe_map, nnode )
       import :: SB_discrete_integration_t, SB_p_volume_integrator_t, rp, fe_map_t, &
                 SB_quadrature_t, ip
       implicit none
       class(SB_discrete_integration_t), intent(inout) :: this
       type(SB_p_volume_integrator_t)    , intent(in) :: vol_int(:)
       real(rp), intent(inout) :: elmat(:,:), elvec(:)
       type(fe_map_t), intent(in) :: fe_map
       type(SB_quadrature_t), intent(in) :: quad
       integer(ip), intent(in) :: nnode(:)
     end subroutine integrate_interface
  end interface

contains

end module SB_discrete_integration_names

!%%%%%%%%%%%%%%%%%%%%%%%

! module poisson_discrete_integration_names
!   use shape_values_names
!   use SB_discrete_integration_names
!   use reference_fe_names
!   use types_names
!   implicit none
! # include "debug.i90"
!   private
!   type, extends(SB_discrete_integration_t) :: poisson_discrete_integration_t
!   integer(ip) :: viscosity
! contains
!   procedure :: integrate
! end type poisson_discrete_integration_t

! public :: poisson_discrete_integration_t

! contains

!   subroutine integrate (this, vol_int, elmat, elvec)
!     implicit none
!     class(poisson_discrete_integration_t), intent(inout) :: this
!     type(SB_volume_integrator_t)    , intent(inout) :: vol_int
!     real(rp), intent(inout) :: elmat(:,:), elvec(:)

!     integer(ip)  :: ndime
!     real(rp)     :: dtinv, c1, c2, mu

!     integer(ip)  :: idime,igaus,inode,jnode,ngaus,nnode
!     real(rp) :: factor

!     type(SB_quadrature_t), pointer :: quad
!     type(SB_interpolation_t), pointer :: int
!     class(reference_fe_t), pointer :: ref_fe
!     type(fe_map_t), pointer :: fe_map

!     type(shape_values_t), pointer :: shape_gradient_test, shape_gradient_trial

!     real(rp), pointer :: grad_test(:,:), grad_trial(:,:)

!     integer(ip) :: a, b

!     quad => vol_int%get_quadrature()
!     int => vol_int%get_interpolation()
!     ref_fe => vol_int%get_reference_fe()
!     fe_map => vol_int%get_fe_map()


!     nnode = ref_fe%get_number_nodes()
!     ngaus = quad%get_number_integration_points()
!     ndime = ref_fe%get_number_dimensions()

!     shape_gradient_test => vol_int%get_gradients()
!     shape_gradient_trial => vol_int%get_gradients()

!     do igaus = 1,ngaus
!        factor = fe_map%get_det_jacobian(igaus) * quad%get_weight(igaus)
!        do inode = 1, nnode
!           grad_trial => shape_gradient_trial%get_value(inode,igaus)
!           do jnode = 1, nnode
!              grad_test => shape_gradient_test%get_value(jnode,igaus)
!              ! elmat = elmat + grad_trial*grad_test
!              do a = 1,size(shape_gradient_test%get_value(inode,igaus),1)
!                 do b= 1,size(shape_gradient_test%get_value(inode,igaus),2)
!                    elmat(inode,jnode) = elmat(inode,jnode) &
!                         + factor * grad_test(a,b) * grad_trial(a,b) 
!                 end do
!              end do

!           end do
!        end do
!     end do

!  !write (*,*) 'XXXXXXXXX ELMAT XXXXXXXXX'
!  !write (*,*) elmat


!  ! Apply boundary conditions
!  !call impose_strong_dirichlet_data(fe) 

! end subroutine integrate


! end module poisson_discrete_integration_names



!%%%%%%%%%%%%%%%%%%%%%%%

module vector_laplacian_discrete_integration_names
  use shape_values_names
  use SB_discrete_integration_names
  use reference_fe_names
  use types_names
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

  subroutine integrate (this, vol_int, elmat, elvec, quad, fe_map, nnode )
    implicit none
    class(vector_laplacian_discrete_integration_t), intent(inout) :: this
    type(SB_p_volume_integrator_t)    , intent(in) :: vol_int(:)
    real(rp), intent(inout) :: elmat(:,:), elvec(:)
    type(fe_map_t), intent(in) :: fe_map
    type(SB_quadrature_t), intent(in) :: quad
    integer(ip), intent(in) :: nnode(:)

    integer(ip)  :: ndime
    real(rp)     :: dtinv, c1, c2, mu

    integer(ip)  :: idime,igaus,inode,jnode,ngaus
    real(rp) :: factor

    !type(SB_quadrature_t), pointer :: quad

    type(shape_values_t), pointer :: shape_gradient_test_u, shape_gradient_trial_u
    type(shape_values_t), pointer :: shape_gradient_test_v, shape_gradient_trial_v

    real(rp), pointer :: grad_test(:,:), grad_trial(:,:)

    integer(ip) :: a, b

    !quad => vol_int(1)%p%get_quadrature()
    !fe_map => vol_int%get_fe_map()


    !nnode = !ref_fe%get_number_nodes()
    ngaus = quad%get_number_integration_points()

    shape_gradient_test_u => vol_int(1)%p%get_gradients()
    shape_gradient_trial_u => vol_int(1)%p%get_gradients()
    shape_gradient_test_v => vol_int(2)%p%get_gradients()
    shape_gradient_trial_v => vol_int(2)%p%get_gradients()

    do igaus = 1,ngaus
       factor = fe_map%get_det_jacobian(igaus) * quad%get_weight(igaus)
       do inode = 1, nnode(1)
          grad_trial => shape_gradient_test_u%get_value(inode,igaus)
          do jnode = 1, nnode(1)
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
       do inode = 1, nnode(2)
          grad_trial => shape_gradient_test_u%get_value(inode,igaus)
          do jnode = 1, nnode(2)
             grad_test => shape_gradient_trial_u%get_value(jnode,igaus)
             ! elmat = elmat + grad_trial*grad_test
             do a = 1,size(shape_gradient_test_v%get_value(inode,igaus),1)
                do b= 1,size(shape_gradient_trial_v%get_value(inode,igaus),2)
                   ! write(*,*) 'BLOCK 2',factor * grad_test(a,b) * grad_trial(a,b)
                   elmat(nnode(1)+inode,nnode(1)+jnode) = elmat(nnode(1)+inode,nnode(1)+jnode) &
                        + factor * grad_test(a,b) * grad_trial(a,b) 
                end do
             end do
          end do
       end do
    end do

    !write (*,*) 'XXXXXXXXX ELMAT XXXXXXXXX'
    !write (*,*) elmat


    ! Apply boundary conditions
    !call impose_strong_dirichlet_data(fe) 

  end subroutine integrate


end module vector_laplacian_discrete_integration_names
