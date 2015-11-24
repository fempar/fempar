module reference_fe_factory_names
  use reference_fe_names
!  use quad_lagrangian_reference_fe_names
  use types_names
  implicit none
# include "debug.i90"
  private

  public :: start_reference_fe

contains

  function start_reference_fe ( topology, fe_type, number_dimensions, order, field_type, continuity )
    implicit none 
    character(*), intent(in) :: topology, fe_type
    integer(ip), intent(in)  :: number_dimensions, order
    character(*), optional, intent(in) :: field_type
    logical, optional, intent(in) :: continuity
    
    class(reference_fe_t), pointer :: start_reference_fe
    integer(ip) :: field_components, i

    if ( topology == "quad" ) then
       if ( fe_type == "Lagrangian") then
          allocate ( quad_lagrangian_reference_fe_t :: start_reference_fe )
       else
          write (*,*) 'ERROR: ELEMENT TYPE NOT SUPPORTED'
       end if
    end if

    field_components = 1
    if ( present(field_type) ) then
       if ( field_type == "vector" ) then
          field_components = number_dimensions
       else if( field_type == "scalar" ) then
          field_components = 1
       else if ( field_type == "tensor" ) then
          field_components = number_dimensions**2
       else if ( field_type == "symmetric_tensor" ) then
          field_components = 0
          do i = 1, number_dimensions
             field_components = field_components + i
          end do
       else
          write (*,*) 'ERROR: field type not defined'
          assert ( 0 == 1)
       end if
    end if
       
    

    call start_reference_fe%create( number_dimensions, order, field_components, continuity ) 

  end function start_reference_fe

end module reference_fe_factory_names


