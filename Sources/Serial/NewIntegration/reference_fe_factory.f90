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
    
    call start_reference_fe%create( number_dimensions, order, field_type, continuity ) 

  end function start_reference_fe

end module reference_fe_factory_names


