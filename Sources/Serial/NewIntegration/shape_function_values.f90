module shape_values_names
  use types_names
  use memor_names
  implicit none
# include "debug.i90"

  private

  integer(ip), parameter :: dim = 2 ! SB.alert : To be used by templatization mechanism

  type, abstract :: field_type_t
   contains
     procedure(get_value_interface), deferred :: get_value
  end type field_type_t

  abstract interface
     function get_value_interface( this )
       import :: field_type_t, rp
       class(field_type_t), target, intent(in) :: this
       real(rp), pointer :: get_value_interface(:,:)
     end function get_value_interface
  end interface
  !  type p_field_type_t
  !     class(field_type_t), pointer :: p
  !  end type p_field_type_t

  type, extends(field_type_t) :: scalar_field_t
  real(rp) :: value(1,1)
contains
  procedure :: get_value => scalar_field_get_value
end type scalar_field_t

type, extends(field_type_t) :: vector_field_t
real(rp) :: value(dim,1)
contains
procedure :: get_value => vector_field_get_value
end type vector_field_t

type, extends(field_type_t) :: tensor_field_t
real(rp)  :: value(dim,dim)
contains
procedure :: get_value => tensor_field_get_value
end type tensor_field_t

type, extends(field_type_t) :: symmetric_tensor_field_t
real(rp)  :: value(dim,dim)
contains
procedure :: get_value => symmetric_tensor_field_get_value
end type symmetric_tensor_field_t

public :: field_type_t, scalar_field_t, vector_field_t, tensor_field_t, symmetric_tensor_field_t
type shape_values_t
   private
   character(:), allocatable :: field_type ! scalar, vector, tensor, symmetric_tensor
   
   integer(ip) :: number_shape_functions
   integer(ip) :: number_evaluation_points
   integer(ip) :: number_space_dimensions
   integer(ip) :: number_field_components
   
   class(field_type_t), allocatable :: shape_values_array(:,:)
 contains
   procedure :: start => shape_values_start
   procedure :: print
   procedure :: get_value
end type shape_values_t

! Types
public :: shape_values_t, gradient_field_type

contains
  
  subroutine shape_values_start( this, field_type, number_field_components, number_evaluation_points, &
       number_shape_functions )
    implicit none
    class(shape_values_t), intent(inout) :: this
    character(*), intent(in)  :: field_type
    integer(ip), intent(in) :: number_field_components, number_evaluation_points, number_shape_functions
    integer(ip) :: istat
    
    this%number_shape_functions = number_shape_functions
    this%number_evaluation_points = number_evaluation_points
!this%number_space_dimensions = number_space_dimensions
    this%number_field_components = number_field_components
    
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
    
  end subroutine shape_values_start
  
  function get_value( this, i, j)
implicit none
class(shape_values_t), target, intent(in) :: this
integer(ip), intent(in) :: i, j
real(rp), pointer :: get_value(:,:)

select type ( entry => this%shape_values_array(i,j) )
class is ( scalar_field_t )
get_value => entry%value
class is ( vector_field_t )
get_value => entry%value
class is ( tensor_field_t )
get_value => entry%value
class is ( symmetric_tensor_field_t )
get_value => entry%value
class default 
check( .false. )
end select

!class(field_type_t), pointer :: get_entry
!get_entry => this%shape_values_array(i,j)

end function get_value

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

function scalar_field_get_value( this )
class(scalar_field_t), target, intent(in) :: this
real(rp), pointer :: scalar_field_get_value(:,:)
scalar_field_get_value => this%value
end function scalar_field_get_value

function vector_field_get_value( this )
class(vector_field_t), target, intent(in) :: this
real(rp), pointer :: vector_field_get_value(:,:)
vector_field_get_value => this%value
end function vector_field_get_value

function tensor_field_get_value( this )
class(tensor_field_t), target, intent(in) :: this
real(rp), pointer :: tensor_field_get_value(:,:)
tensor_field_get_value => this%value
end function tensor_field_get_value

function symmetric_tensor_field_get_value( this )
class(symmetric_tensor_field_t), target, intent(in) :: this
real(rp), pointer :: symmetric_tensor_field_get_value(:,:)
symmetric_tensor_field_get_value => this%value
end function symmetric_tensor_field_get_value


subroutine print( this)
implicit none
class(shape_values_t), target, intent(in) :: this
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
end subroutine print

end module shape_values_names
