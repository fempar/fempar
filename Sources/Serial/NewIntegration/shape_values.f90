module shape_values_names
  use types_names
  use memor_names
  implicit none
# include "debug.i90"

  private

  integer(ip), parameter :: dim = 2 ! SB.alert : To be used by templatization mechanism

  type, abstract :: field_type_t
  end type field_type_t

  type, extends(field_type_t) :: scalar_field_t
     private
     real(rp) :: value(1,1)
  end type scalar_field_t
  
  type, extends(field_type_t) :: vector_field_t
     private
     real(rp) :: value(dim,1)
  end type vector_field_t
  
  type, extends(field_type_t) :: tensor_field_t
     private
     real(rp)  :: value(dim,dim)
  end type tensor_field_t

  type, extends(field_type_t) :: symmetric_tensor_field_t
     private
     real(rp)  :: value(dim,dim)
  end type symmetric_tensor_field_t

  public :: field_type_t, scalar_field_t, vector_field_t, tensor_field_t, symmetric_tensor_field_t

  type shape_values_t
     private
     character(:), allocatable        :: field_type ! scalar, vector, tensor, symmetric_tensor
     integer(ip)                      :: number_shape_functions
     integer(ip)                      :: number_evaluation_points
     class(field_type_t), allocatable :: shape_values_array(:,:)
   contains
     procedure :: create       => shape_values_create
     procedure :: free         => shape_values_free
     procedure :: print        => shape_values_print
     procedure :: get_value    => shape_values_get_value
  end type shape_values_t

  ! Types
  public :: shape_values_t, gradient_field_type

contains
  
  subroutine shape_values_create( this, &
                                 field_type, &
                                 number_shape_functions, &
                                 number_evaluation_points )
    implicit none
    class(shape_values_t), intent(inout) :: this
    character(*)         , intent(in)    :: field_type
    integer(ip)          , intent(in)    :: number_shape_functions
    integer(ip)          , intent(in)    :: number_evaluation_points
    
    this%field_type               = field_type
    this%number_shape_functions   = number_shape_functions
    this%number_evaluation_points = number_evaluation_points
    
    if ( field_type == 'scalar') then
       allocate ( scalar_field_t :: this%shape_values_array(number_shape_functions,number_evaluation_points) )
    else if( field_type == 'vector') then
       allocate ( vector_field_t :: this%shape_values_array(number_shape_functions,number_evaluation_points) )
    else if( field_type == 'tensor') then
       allocate ( tensor_field_t :: this%shape_values_array(number_shape_functions,number_evaluation_points) )
    else if( field_type == 'symmetric_tensor') then
       allocate ( symmetric_tensor_field_t :: this%shape_values_array(number_shape_functions,number_evaluation_points) )
    else
       write (*,*) 'ERROR: FIELD TYPE UNDECLARED'
       check( .false. )
    end if
    
  end subroutine shape_values_create
  
  subroutine shape_values_free( this )
    implicit none
    class(shape_values_t), intent(inout) :: this
    
    if (allocated(this%field_type)) deallocate(this%field_type)
    this%number_shape_functions   = 0
    this%number_evaluation_points = 0
    if (allocated(this%shape_values_array)) deallocate(this%shape_values_array)
  end subroutine shape_values_free
  
    subroutine shape_values_print(this)
    implicit none
    class(shape_values_t), intent(in) :: this
    integer(ip) :: i, j
    real(rp), pointer :: get_value(:,:)
    
    write(*,*) "*********SHAPE_VALUES_PRINT*********",this%number_shape_functions,this%number_evaluation_points
    do i = 1,this%number_shape_functions
       do j = 1,this%number_evaluation_points
          write(*,*) "*********SHAPE FUNCTION ",i,"*********"
          select type ( entry => this%shape_values_array(i,j) )
          class is ( scalar_field_t )
             write(*,*)  entry%value
          class is ( vector_field_t )
             write(*,*)  entry%value
          class is ( tensor_field_t )
             write(*,*)  entry%value
          class is ( symmetric_tensor_field_t )
             write(*,*)  entry%value
          class default 
             check( .false. )
          end select
       end do
    end do
  end subroutine shape_values_print
  
  function shape_values_get_value( this, i, j) 
    implicit none
    class(shape_values_t), target, intent(in) :: this
    integer(ip)                  , intent(in) :: i, j
    real(rp)                        , pointer :: shape_values_get_value(:,:)

    select type ( entry => this%shape_values_array(i,j) )
       class is ( scalar_field_t )
       shape_values_get_value => entry%value
       class is ( vector_field_t )
       shape_values_get_value => entry%value
       class is ( tensor_field_t )
       shape_values_get_value => entry%value
       class is ( symmetric_tensor_field_t )
       shape_values_get_value => entry%value
       class default 
       check( .false. )
    end select
  end function shape_values_get_value

  function gradient_field_type ( field_type)
    implicit none
    character(*), intent(in)  :: field_type
    character(:), allocatable :: gradient_field_type
    
    if ( field_type == 'scalar') then 
       gradient_field_type = "vector"
    else if( field_type == 'vector') then
       gradient_field_type = "tensor"    
    else
       write (*,*) 'ERROR: GRADIENT OF FIELD TYPE UNDECLARED'
    end if
  end function gradient_field_type
  
end module shape_values_names
