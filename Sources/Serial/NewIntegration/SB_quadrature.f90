module SB_quadrature_names
  use types_names
  use memor_names
  implicit none
# include "debug.i90"
  private

  type SB_quadrature_t
     private
     integer(ip)           :: &
          number_dimensions,    &
          number_integration_points                      ! Number of integration points
     real(rp), allocatable :: &
          coordinates(:,:),             &    ! Quadrature points position
          weight(:)                  ! Quadrature points weight
   contains
     procedure :: create
     !procedure :: free
     procedure :: print

     procedure :: get_number_dimensions
     procedure :: get_number_integration_points

     procedure :: get_pointer_coordinates
     procedure :: get_pointer_weight
  end type SB_quadrature_t

  ! Types
  public :: SB_quadrature_t

contains

  subroutine create ( this, ndime, ngaus )
    implicit none       
    integer(ip), intent(in) :: ndime, ngaus
    class(SB_quadrature_t), intent(out) :: this
    this%number_dimensions = ndime
    this%number_integration_points = ngaus
    call memalloc(ndime,ngaus,this%coordinates,__FILE__,__LINE__)
    call memalloc(ngaus,this%weight,__FILE__,__LINE__)
    this%coordinates=0.0_rp
    this%weight=0.0_rp
  end subroutine create

  function get_number_dimensions ( this )
    implicit none
    class(SB_quadrature_t), intent(in) :: this
    integer(ip) :: get_number_dimensions
    get_number_dimensions = this%number_dimensions
  end function get_number_dimensions

  function get_number_integration_points ( this )
    implicit none
    class(SB_quadrature_t), intent(in) :: this
    integer(ip) :: get_number_integration_points
    get_number_integration_points = this%number_integration_points
  end function get_number_integration_points

  function get_pointer_coordinates ( this )
    implicit none
    class(SB_quadrature_t), target, intent(in) :: this
    real(rp), pointer :: get_pointer_coordinates(:,:)
    get_pointer_coordinates => this%coordinates
  end function get_pointer_coordinates

  function get_pointer_weight ( this )
    implicit none
    class(SB_quadrature_t), target, intent(in) :: this
    real(rp), pointer :: get_pointer_weight(:)
    get_pointer_weight => this%weight
  end function get_pointer_weight

  subroutine print ( this )
    implicit none
    class(SB_quadrature_t), target, intent(in) :: this
    write(*,*) 'number_dimensions: ', this%number_dimensions
    write(*,*) 'number_integration_points: ', this%number_integration_points
    write(*,*) 'coordinates: ', this%coordinates
    write(*,*) 'weight: ', this%weight
  end subroutine print


end module SB_quadrature_names
