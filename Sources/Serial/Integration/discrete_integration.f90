module discrete_integration_names
  use field_names
  use reference_fe_names
  use types_names
  use assembler_names
  use serial_fe_space_names
  use memor_names

  implicit none
# include "debug.i90"
  private

  type, abstract :: discrete_integration_t
   contains
     procedure (integrate_interface), deferred :: integrate
     procedure                                 :: impose_strong_dirichlet_data
  end type discrete_integration_t

  type p_discrete_integration_t
     class(discrete_integration_t), pointer :: p => NULL()  
  end type p_discrete_integration_t

  public :: discrete_integration_t, p_discrete_integration_t

  abstract interface
     subroutine integrate_interface ( this, fe_space, assembler  )
       import :: discrete_integration_t, serial_fe_space_t, assembler_t
       implicit none
       class(discrete_integration_t), intent(in)    :: this
       class(serial_fe_space_t)     , intent(inout) :: fe_space
       class(assembler_t)           , intent(inout) :: assembler
     end subroutine integrate_interface
  end interface

contains

subroutine impose_strong_dirichlet_data ( this, elmat, elvec, code, value, nnode, num_fe_spaces )
 implicit none
 class(discrete_integration_t) :: this
 real(rp), intent(in) :: elmat(:,:)
 real(rp), intent(inout) :: elvec(:)  
 type(i1p_t), intent(in) :: code(:)
 type(r1p_t), intent(in) :: value(:)
 integer(ip), intent(in) :: nnode(:), num_fe_spaces
 integer(ip) :: i, c, ifes

 c = 0
 do ifes = 1, num_fe_spaces
    do i = 1, nnode(ifes)
       c = c + 1
       if ( code(ifes)%p(i) /= 0 ) then
          elvec = elvec - elmat(:,c)*value(ifes)%p(i)
       end if
    end do
 end do

end subroutine impose_strong_dirichlet_data

end module discrete_integration_names
