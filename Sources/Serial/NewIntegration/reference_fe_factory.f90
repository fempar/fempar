module reference_fe_factory_names
  use reference_fe_names
!  use quad_lagrangian_reference_fe_names
  use types_names
  implicit none
# include "debug.i90"
  private

  public :: start_reference_fe

contains

  function start_reference_fe ( topology, fe_type, number_dimensions, order, continuity )
    implicit none 
    character(*), intent(in) :: topology, fe_type
    integer(ip), intent(in)  :: number_dimensions, order
    logical, optional, intent(in) :: continuity
    class(reference_fe_t), pointer :: start_reference_fe

    if ( topology == "quad" ) then
       if ( fe_type == "Lagrangian") then
          allocate ( quad_lagrangian_reference_fe_t :: start_reference_fe )
       else
          write (*,*) 'ERROR: ELEMENT TYPE NOT SUPPORTED'
       end if
    end if

    call start_reference_fe%create( number_dimensions, order, continuity ) 

  end function start_reference_fe

end module reference_fe_factory_names


