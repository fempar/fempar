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
          number_dimensions,        &      
          number_shape_functions,   &      
          number_evaluation_points, &      
          number_entries_symmetric_tensor
     ! (usually integration points)

     real(rp), allocatable      ::  &
          shape_functions(:,:),     &   ! Shape functions
          shape_derivatives(:,:,:), &   ! Derivatives
          hessian(:,:,:)                ! Hessian

   contains

     procedure :: create
     !procedure :: free
     procedure :: print

     procedure :: get_number_dimensions
     procedure :: get_number_shape_functions
     procedure :: get_number_evaluation_points
     procedure :: get_number_entries_symmetric_tensor

     procedure :: get_shape_function
     procedure :: get_shape_derivative
     procedure :: get_hessian 

     !procedure :: get_shape_functions
     !procedure :: get_shape_derivatives
     !procedure :: get_hessian 

     procedure :: get_pointer_shape_functions
     procedure :: get_pointer_shape_derivatives
     procedure :: get_pointer_hessian

  end type SB_interpolation_t

  public :: SB_interpolation_t

contains

  function get_number_dimensions( this )
    implicit none
    class(SB_interpolation_t), intent(in) :: this
    integer(ip) :: get_number_dimensions
    get_number_dimensions = this%number_dimensions
  end function get_number_dimensions

  function get_number_shape_functions( this )
    implicit none
    class(SB_interpolation_t), intent(in) :: this
    integer(ip) :: get_number_shape_functions
    get_number_shape_functions = this%number_shape_functions
  end function get_number_shape_functions

    
  function get_number_evaluation_points( this )
    implicit none
    class(SB_interpolation_t), intent(in) :: this
    integer(ip) :: get_number_evaluation_points
    get_number_evaluation_points = this%number_evaluation_points
  end function get_number_evaluation_points


  function get_number_entries_symmetric_tensor( this )
    implicit none
    class(SB_interpolation_t), intent(in) :: this
    integer(ip) :: get_number_entries_symmetric_tensor
    get_number_entries_symmetric_tensor = this%number_entries_symmetric_tensor
  end function get_number_entries_symmetric_tensor

  subroutine create( this, ndime, nnode, ngaus, ntens, khes )
    implicit none
    class(SB_interpolation_t), intent(out) :: this
    integer(ip)          , intent(in)     :: nnode, ndime, ngaus
    integer(ip) :: iloc,ntens
    logical, optional :: khes

    this%number_dimensions = ndime
    this%number_shape_functions = nnode
    this%number_evaluation_points = ngaus
    this%number_entries_symmetric_tensor = ntens
    call memalloc(nnode,ngaus,this%shape_functions,__FILE__,__LINE__)
    call memalloc(ndime,nnode,ngaus,this%shape_derivatives,   __FILE__,__LINE__)
    this%shape_functions = 0.0_rp
    this%shape_derivatives = 0.0_rp
    if ( present(khes) ) then
       if ( khes ) then 
          call memalloc(ntens,nnode,ngaus,this%hessian,   __FILE__,__LINE__) 
          this%hessian = 0.0_rp
       end if
    end if

  end subroutine create

  subroutine print ( this )
    implicit none
    class(SB_interpolation_t), intent(in) :: this
    write(*,*) 'number_dimensions: ', this%number_dimensions
    write(*,*) 'number_shape_functions: ', this%number_shape_functions
    write(*,*) 'number_evaluation_points: ', this%number_evaluation_points
    write(*,*) 'number_entries_symmetric_tensor: ', this%number_entries_symmetric_tensor
    write(*,*) 'shape_functions: ', this%shape_functions
    write(*,*) 'shape_derivatives: ', this%shape_derivatives
    if ( allocated( this%hessian ) ) then
       write(*,*) 'hessian: ', this%hessian
    else
       write(*,*) 'hessian not computed '
    end if
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
    if ( allocated( this%hessian ) ) then
       get_pointer_hessian => this%hessian
    else
       nullify(get_pointer_hessian)
    end if
  end function get_pointer_hessian

  function get_shape_function ( this, i, j )
    implicit none
    class(SB_interpolation_t), target, intent(in) :: this
    real(rp) :: get_shape_function
    integer(ip) :: i, j
    
    get_shape_function = this%shape_functions(i,j)

  end function get_shape_function

  function get_shape_derivative ( this, i, j, k )
    implicit none
    class(SB_interpolation_t), target, intent(in) :: this
    real(rp) :: get_shape_derivative
    integer(ip) :: i, j, k

    get_shape_derivative = this%shape_derivatives(i,j,k)

  end function get_shape_derivative

  function get_hessian ( this, i, j, k )
    implicit none
    class(SB_interpolation_t), target, intent(in) :: this
    real(rp) :: get_hessian
    integer(ip) :: i, j, k
    
    get_hessian = this%hessian(i,j,k)

  end function get_hessian

end module SB_interpolation_names
