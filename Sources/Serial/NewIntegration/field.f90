module field_names
  use types_names
  use memor_names
  implicit none
# include "debug.i90"

  private

  integer(ip), parameter :: dim = 3
  
  type :: scalar_field_t
     private
     real(rp) :: value
   contains
     procedure, non_overridable :: init  => scalar_field_init
     procedure, non_overridable :: set   => scalar_field_set
  end type scalar_field_t
  
  type :: vector_field_t
     private
     real(rp) :: value(dim)
   contains
     procedure, non_overridable :: init  => vector_field_init
     procedure, non_overridable :: set   => vector_field_set		
     procedure, non_overridable :: add   => vector_field_add
  end type vector_field_t
  
  type :: tensor_field_t
     private
     real(rp)  :: value(dim,dim)
   contains
     procedure, non_overridable :: init  => tensor_field_init
     procedure, non_overridable :: set   => tensor_field_set
  end type tensor_field_t

  type :: symmetric_tensor_field_t
     private
     real(rp)  :: value(dim,dim)
   contains			
     procedure, non_overridable :: init  => symmetric_tensor_field_init
     procedure, non_overridable :: set   => symmetric_tensor_field_set					
  end type symmetric_tensor_field_t
  
  interface operator(*)
    module procedure single_contract_vector_vector, single_contract_tensor_vector, single_contract_vector_tensor
    module procedure scal_left_vector, scal_right_vector, scal_left_tensor, scal_right_tensor
  end interface
  
  interface double_contract
    module procedure double_contract_tensor_tensor
  end interface
  
  public :: vector_field_t, tensor_field_t, symmetric_tensor_field_t 
  public :: operator(*)
  public :: double_contract
  ! public :: scalar_field_t (not actually needed, used real(rp) instead)

contains

  subroutine scalar_field_init(this,value)
    implicit none
    class(scalar_field_t), intent(inout) :: this
    real(rp)             , intent(in)    :: value
    this%value = value
  end subroutine scalar_field_init

  subroutine scalar_field_set(this,value)
    implicit none
    class(scalar_field_t), intent(inout) :: this
    real(rp)             , intent(in)    :: value
    this%value = value
  end subroutine scalar_field_set
		
  subroutine vector_field_init(this,value)
    implicit none
    class(vector_field_t), intent(inout) :: this
    real(rp)             , intent(in)    :: value
    this%value = value
  end subroutine vector_field_init

  subroutine vector_field_set(this,i,value)
    implicit none
    class(vector_field_t), intent(inout) :: this
    integer(ip)          , intent(in)    :: i
    real(rp)             , intent(in)    :: value
    this%value(i) = value
  end subroutine vector_field_set
  
  subroutine vector_field_add(this,i,value)
    implicit none
    class(vector_field_t), intent(inout) :: this
    integer(ip)          , intent(in)    :: i
    real(rp)             , intent(in)    :: value
    this%value(i) = this%value(i) + value
  end subroutine vector_field_add
		
  subroutine tensor_field_init(this,value)
    implicit none
    class(tensor_field_t), intent(inout) :: this
    real(rp)             , intent(in)    :: value
    this%value = value
  end subroutine tensor_field_init

  subroutine tensor_field_set(this,i,j,value)
    implicit none
    class(tensor_field_t), intent(inout) :: this
    integer(ip)          , intent(in)    :: i
    integer(ip)          , intent(in)    :: j
    real(rp)             , intent(in)    :: value
    this%value(i,j) = value
  end subroutine tensor_field_set
		
  subroutine symmetric_tensor_field_init(this,value)
    implicit none
    class(symmetric_tensor_field_t), intent(inout) :: this
    real(rp)                       , intent(in)    :: value
    this%value = value
  end subroutine symmetric_tensor_field_init

  subroutine symmetric_tensor_field_set(this,i,j,value)
    implicit none
    class(symmetric_tensor_field_t), intent(inout) :: this
    integer(ip)                    , intent(in)    :: i
    integer(ip)                    , intent(in)    :: j
    real(rp)                       , intent(in)    :: value
    assert(j>=i)
    this%value(i,j) = value
  end subroutine symmetric_tensor_field_set

  function single_contract_vector_vector(v1,v2) result(res)
    implicit none
    type(vector_field_t), intent(in) :: v1
    type(vector_field_t), intent(in) :: v2
    real(rp)                         :: res
    integer(ip) :: k
    res=0.0_rp
    do k=1,dim
      res = res + v1%value(k)*v2%value(k)
    end do
   end function single_contract_vector_vector
   
   function single_contract_tensor_vector(t,v) result(res)
    implicit none
    type(tensor_field_t), intent(in) :: t
    type(vector_field_t), intent(in) :: v
    type(vector_field_t)             :: res
    integer(ip) :: i, k
    res%value=0.0_rp
    do k=1,dim
      do i=1,dim
        res%value(i) = res%value(i) + t%value(i,k) * v%value(k)
      end do
    end do
   end function single_contract_tensor_vector
   
   function single_contract_vector_tensor(v,t) result(res)
    implicit none
    type(vector_field_t), intent(in) :: v
    type(tensor_field_t), intent(in) :: t
    type(vector_field_t)             :: res
    integer(ip) :: i, k
    res%value=0.0_rp
    do i=1,dim
      do k=1,dim
        res%value(i) = res%value(i) + v%value(k) * t%value(k,i)
      end do
    end do
   end function single_contract_vector_tensor
   
   function scal_left_vector(alpha,v) result(res)
    implicit none
    real(rp)            , intent(in) :: alpha
    type(vector_field_t), intent(in) :: v
    type(vector_field_t)             :: res
    res%value = alpha * v%value
   end function scal_left_vector
   
   function scal_right_vector(v,alpha) result(res)
    implicit none
    type(vector_field_t), intent(in) :: v
    real(rp)            , intent(in) :: alpha
    type(vector_field_t)             :: res
    res%value = alpha * v%value
   end function scal_right_vector
   
   function scal_left_tensor(alpha,t) result(res)
    implicit none
    real(rp)            , intent(in) :: alpha
    type(tensor_field_t), intent(in) :: t
    type(tensor_field_t)             :: res
    res%value = alpha * t%value
   end function scal_left_tensor
   
   function scal_right_tensor(t,alpha) result(res)
    implicit none
    type(tensor_field_t), intent(in) :: t
    real(rp)            , intent(in) :: alpha
    type(tensor_field_t)             :: res
    res%value = alpha * t%value
   end function scal_right_tensor
   
   function double_contract_tensor_tensor(t1,t2) result(res)
    implicit none
    type(tensor_field_t), intent(in) :: t1
    type(tensor_field_t), intent(in) :: t2
    real(rp)                         :: res
    integer(ip) :: i, j
    res = 0.0_rp
    do j=1, dim
      do i=1,dim
        res = res + t1%value(i,j)*t2%value(i,j)
      end do
    end do
   end function double_contract_tensor_tensor
  
end module field_names
