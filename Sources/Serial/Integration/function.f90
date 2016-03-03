module function_names
  use types_names
  !use memor_names
  use field_names

  implicit none
# include "debug.i90"

  private
  ! Abstract array function
  type :: array_function_t
     private
     integer(ip)  :: number_components = 1
     !character(*) :: function_type = 'space' ! 'space' or 'space_time'  
   contains
     procedure, non_overridable  :: set_number_components => array_set_number_components
     
     procedure :: get_component_value_space => array_get_component_value_space
     procedure :: get_component_value_space_time => array_get_component_value_space_time
     generic   :: get_component_value => get_component_value_space, get_component_value_space_time

     procedure :: get_component_values_set_space => array_get_component_values_set_space
     procedure :: get_component_values_set_space_time => array_get_component_values_set_space_time
     generic   :: get_component_values_set => get_component_values_set_space, get_component_values_set_space_time

     procedure :: get_value_space => array_get_value_space
     procedure :: get_value_space_time => array_get_value_space_time
     generic   :: get_value => get_value_space, get_value_space_time

     procedure :: get_values_set_space => array_get_values_set_space
     procedure :: get_values_set_space_time => array_get_values_set_space_time
     generic   :: get_values_set => get_values_set_space, get_values_set_space_time     
  end type array_function_t

  type :: scalar_function_t
     private
     !character(*) :: function_type = 'space' ! 'space' or 'space_time'  
   contains
     procedure :: get_value_space => scalar_get_value_space
     procedure :: get_value_space_time => scalar_get_value_space_time
     generic   :: get_value => get_value_space, get_value_space_time

     procedure :: get_values_set_space => scalar_get_values_set_space
     procedure :: get_values_set_space_time => scalar_get_values_set_space_time
     generic   :: get_values_set => get_values_set_space, get_values_set_space_time     
  end type scalar_function_t

  type :: vector_function_t
     private
     !character(*) :: function_type = 'space' ! 'space' or 'space_time'  
   contains
     procedure :: get_value_space => vector_get_value_space
     procedure :: get_value_space_time => vector_get_value_space_time
     generic   :: get_value => get_value_space, get_value_space_time

     procedure :: get_values_set_space => vector_get_values_set_space
     procedure :: get_values_set_space_time => vector_get_values_set_space_time
     generic   :: get_values_set => get_values_set_space, get_values_set_space_time     
  end type vector_function_t

  type :: tensor_function_t
     private
     !character(*) :: function_type = 'space' ! 'space' or 'space_time'  
   contains
     procedure :: get_value_space => tensor_get_value_space
     procedure :: get_value_space_time => tensor_get_value_space_time
     generic   :: get_value => get_value_space, get_value_space_time

     procedure :: get_values_set_space => tensor_get_values_set_space
     procedure :: get_values_set_space_time => tensor_get_values_set_space_time
     generic   :: get_values_set => get_values_set_space, get_values_set_space_time     
  end type tensor_function_t

  public :: array_function_t, scalar_function_t, vector_function_t, tensor_function_t

contains
! One only needs to fill array_get_component_value_space or array_get_component_value_space_time
! (depending on the type of function). The rest of TBPs are implemented based on that at the 
! root array_function_t. For efficiency purposes, these functions can also be filled in user
! defined functions that inherit from array_function_t
#include "sbm_array_function.i90"
! One only needs to fill scalar_get__value_space or scalar_get_value_space_time
! (depending on the type of function). The rest of TBPs are implemented based on that at the 
! root scalar_function_t. For efficiency purposes, these functions can also be filled in user
! defined functions that inherit from scalar_function_t
#include "sbm_scalar_function.i90"
! One only needs to fill vector_get__value_space or vector_get_value_space_time
! (depending on the type of function). The rest of TBPs are implemented based on that at the 
! root vector_function_t. For efficiency purposes, these functions can also be filled in user
! defined functions that inherit from vector_function_t
#include "sbm_vector_function.i90"
! One only needs to fill tensor_get__value_space or tensor_get_value_space_time
! (depending on the type of function). The rest of TBPs are implemented based on that at the 
! root tensor_function_t. For efficiency purposes, these functions can also be filled in user
! defined functions that inherit from tensor_function_t
#include "sbm_tensor_function.i90"

end module function_names
