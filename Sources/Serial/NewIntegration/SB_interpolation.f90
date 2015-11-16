module SB_interpolation_names
  use SB_quadrature_names
  use types_names
  use memor_names
  implicit none
# include "debug.i90"
  private

  type SB_interpolation_t
     private
     integer(ip)                ::  &
          number_dimensions,        &   ! Number of dimensions
          number_shape_functions,   &   ! Number of shape functions
          number_evaluation_points      ! Number of evaluation points 
     ! (usually integration points)

     real(rp), allocatable      ::  &
          shape_functions(:,:),     &   ! Shape functions
          shape_derivatives(:,:,:), &   ! Derivatives
          hessian(:,:,:)                ! Hessian

   contains

     procedure :: create
     !procedure :: free
     procedure :: print

     !procedure :: get_shape_functions
     !procedure :: get_shape_derivatives
     !procedure :: get_hessian  

     procedure :: get_pointer_shape_functions
     procedure :: get_pointer_shape_derivatives
     procedure :: get_pointer_hessian

  end type SB_interpolation_t

  public :: SB_interpolation_t

contains

  subroutine create( this, ndime, nnode, ngaus, khes )
    implicit none
    class(SB_interpolation_t), intent(out) :: this
    integer(ip)          , intent(in)     :: nnode, ndime, ngaus
    integer(ip) :: iloc,ntens
    logical, optional :: khes

    this%number_dimensions = ndime
    this%number_shape_functions = nnode
    this%number_evaluation_points = ngaus

    call memalloc(nnode,ngaus,this%shape_functions,__FILE__,__LINE__)
    call memalloc(ndime,nnode,ngaus,this%shape_derivatives,   __FILE__,__LINE__)
    this%shape_functions = 0.0_rp
    this%shape_derivatives = 0.0_rp
    !if(present(khes).and.khes) then
       call memalloc(ntens,nnode,ngaus,this%hessian,   __FILE__,__LINE__) 
       this%hessian = 0.0_rp
    !end if

  end subroutine create

  subroutine print ( this )
    implicit none
    class(SB_interpolation_t), intent(in) :: this
    write(*,*) 'number_dimensions: ', this%number_dimensions
    write(*,*) 'number_shape_functions: ', this%number_shape_functions
    write(*,*) 'number_evaluation_points: ', this%number_evaluation_points
    write(*,*) 'shape_functions: ', this%shape_functions
    write(*,*) 'shape_derivatives: ', this%shape_derivatives
  end subroutine print

  function get_pointer_shape_functions ( this )
    implicit none
    class(SB_interpolation_t), target, intent(in) :: this
    real(rp), pointer :: get_pointer_shape_functions(:,:)
    get_pointer_shape_functions => this%shape_functions
  end function get_pointer_shape_functions

  function get_pointer_shape_derivatives ( this )
    implicit none
    class(SB_interpolation_t), target, intent(in) :: this
    real(rp), pointer :: get_pointer_shape_derivatives(:,:,:)
    get_pointer_shape_derivatives => this%shape_derivatives
  end function get_pointer_shape_derivatives

  function get_pointer_hessian ( this )
    implicit none
    class(SB_interpolation_t), target, intent(in) :: this
    real(rp), pointer :: get_pointer_hessian(:,:,:)
    get_pointer_hessian => this%hessian
  end function get_pointer_hessian

end module SB_interpolation_names
