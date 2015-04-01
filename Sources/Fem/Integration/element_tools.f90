! Copyright (C) 2014 Santiago Badia, Alberto F. Martín and Javier Principe
!
! This file is part of FEMPAR (Finite Element Multiphysics PARallel library)
!
! FEMPAR is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! FEMPAR is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with FEMPAR. If not, see <http://www.gnu.org/licenses/>.
!
! Additional permission under GNU GPL version 3 section 7
!
! If you modify this Program, or any covered work, by linking or combining it 
! with the Intel Math Kernel Library and/or the Watson Sparse Matrix Package 
! and/or the HSL Mathematical Software Library (or a modified version of them), 
! containing parts covered by the terms of their respective licenses, the
! licensors of this Program grant you additional permission to convey the 
! resulting work. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# include "debug.i90"
module element_tools_names
  use types
  use fem_space_names
  use memory_guard_names
  use problem_names
  implicit none
  private

  type, extends(memory_guard) :: field
     contains
       !procedure (left_prodcut_interface), deferred :: left_prodcut
       !procedure(right_product_interface), deferred :: right_product
       !generic    :: operator(*) => left_prodcut,right_product
       procedure :: free => free_field
  end type field
  type, extends(field) :: scalar
     real(rp), allocatable :: a(:)     ! a(ngaus)
  end type scalar
  type , extends(field) :: vector
     real(rp), allocatable :: a(:,:)   ! a(:,ngaus)
  end type vector
  type , extends(field) :: tensor
     real(rp), allocatable :: a(:,:,:) ! a(:,:,ngaus)
  end type tensor

  ! We assume that all dofs in this functions are interpolated using the same basis. 
  ! It can be exteneded to the general case, complicating the machinery.
  type, extends(memory_guard) :: fem_function
     integer(ip)                :: ivar=1
     integer(ip)                :: ndof=1
     type(fem_element), pointer :: elem => NULL()
   contains
     procedure :: free => free_fem_function
  end type fem_function

  ! Shape and deriv...
  type, extends(fem_function) :: basis_function
     !real(rp), allocatable :: a(:,:)     ! shape(nnode,ngaus)
     real(rp) :: scaling=0.0_rp
     class(field), allocatable :: left_factor
     class(field), allocatable :: right_factor
   contains
     !procedure, pass(ul) :: scal_left         => scal_left_basis_function
     !procedure, pass(ur) :: scal_right        => scal_right_basis_function
     !procedure, pass(u)  :: product_by_scalar => product_scalar_x_basis_function
     !generic    :: operator(*) => scal_right, scal_left, product_by_scalar
  end type basis_function

  type, extends(fem_function) :: basis_function_gradient
     !real(rp), allocatable :: a(:,:,:)   ! deriv(ndime,nnode,ngaus)
     real(rp) :: scaling=0.0_rp
     class(field), allocatable :: left_factor
     class(field), allocatable :: right_factor
   contains
     ! procedure, pass(ul) :: scal_left => scal_left_basis_function_gradient
     ! procedure, pass(ur) :: scal_right => scal_right_basis_function_gradient
     ! procedure, pass(u)  :: product_by_scalar => product_scalar_x_basis_function
     ! procedure, pass(u)  :: product_by_vector => product_vector_x_basis_function
     ! procedure, pass(u)  :: product_by_tensor => product_tensor_x_basis_function
     ! generic    :: operator(*) => scal_right, scal_left, product_by_scalar, product_by_vector, product_by_tensor
  end type basis_function_gradient

  type, extends(fem_function) :: basis_function_divergence
     !real(rp), allocatable :: a(:,:,:)   ! deriv(ndime,nnode,ngaus)
     real(rp) :: scaling=0.0_rp
     class(field), allocatable :: left_factor
     class(field), allocatable :: right_factor
   contains
     ! procedure, pass(ul) :: scal_left => scal_left_basis_function_divergence
     ! procedure, pass(ur) :: scal_right => scal_right_basis_function_divergence
     ! procedure, pass(u)  :: product_by_scalar => product_scalar_x_basis_function
     ! procedure, pass(u)  :: product_by_vector => product_vector_x_basis_function
     ! procedure, pass(u)  :: product_by_tensor => product_tensor_x_basis_function
     ! generic    :: operator(*) => scal_right, scal_left, product_by_scalar, product_by_vector, product_by_tensor
  end type basis_function_divergence

  ! Interpolations
  type, extends(fem_function) :: given_function
     integer(ip)                :: icomp=1
  end type given_function
  type, extends(fem_function) :: given_function_gradient
     integer(ip)                :: icomp=1
  end type given_function_gradient
  type, extends(fem_function) :: given_function_divergence
     integer(ip)                :: icomp=1
  end type given_function_divergence

  ! Cosmetics...easy to implement if we can return a polymorphic 
  ! allocatable and the compiler works, so we can use the same functions
  ! for both. Currently it is a mess and we don't want to duplicate 
  ! everything.
  !
  ! type, extends(basis_function_value)    :: trial_function
  ! end type trial_function
  ! type, extends(basis_function_gradient) :: trial_function_gradient
  ! end type trial_unction_gradient
  ! type, extends(basis_function_value)    :: test_function
  ! end type test_function
  ! type, extends(basis_function_gradient) :: test_function_gradient
  ! end type test_unction_gradient
  interface basis_function
     module procedure basis_function_constructor
  end interface basis_function
  interface given_function
     module procedure given_function_constructor
  end interface given_function
  interface grad
     module procedure basis_function_gradient_constructor, given_function_gradient_constructor
  end interface grad
  interface div
     module procedure basis_function_divergence_constructor, given_function_divergence_constructor
  end interface div

  interface interpolation
     module procedure scalar_interpolation
     module procedure vector_interpolation
     module procedure scalar_gradient_interpolation
     module procedure vector_gradient_interpolation
     module procedure divergence_interpolation
  end interface interpolation

  public :: basis_function, basis_function_gradient, basis_function_divergence
  public :: given_function, given_function_gradient, given_function_divergence
  public :: scalar, vector, tensor
  public :: grad, div, interpolation

contains


  subroutine free_field ( this )
    class(field), intent(inout) :: this
  end subroutine free_field
  subroutine free_fem_function ( this )
    class(fem_function), intent(inout) :: this
  end subroutine free_fem_function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine copy_fem_function(from,to)
    implicit none
    class(fem_function), intent(in), target  :: from
    class(fem_function), intent(inout)       :: to
    to%ivar=from%ivar
    to%ndof=from%ndof
    to%elem=>from%elem
    call to%SetTemp()
  end subroutine copy_fem_function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! basis_function (some code replication to avoid functions returning polymorphic allocatables)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function basis_function_constructor(prob,ivar,elem) result(var)
    implicit none
    class(physical_problem)  , intent(in) :: prob
    integer(ip)              , intent(in) :: ivar
    type(fem_element), target, intent(in) :: elem
    type(basis_function) :: var
    integer(ip)          :: nnode,ngaus
    var%ivar = ivar
    var%elem => elem
    var%ndof = prob%vars_of_unk(ivar)
    call var%SetTemp()
  end function basis_function_constructor

 function basis_function_gradient_constructor(u) result(g)
    type(basis_function), intent(in) :: u
    type(basis_function_gradient)    :: g
    integer(ip) :: ndime, nnode, ngaus
    call u%GuardTemp()
    call copy_fem_function(u,g)
    call u%CleanTemp()
  end function basis_function_gradient_constructor

 function basis_function_divergence_constructor(u) result(g)
    type(basis_function), intent(in) :: u
    type(basis_function_divergence)    :: g
    integer(ip) :: ndime, nnode, ngaus
    call u%GuardTemp()
    call copy_fem_function(u,g)
    call u%CleanTemp()
  end function basis_function_divergence_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! given_function constructors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function given_function_constructor(prob,ivar,icomp,elem) result(var)
    implicit none
    class(physical_problem) , intent(in) :: prob
    integer(ip)             , intent(in) :: ivar
    integer(ip)             , intent(in) :: icomp
    type(fem_element),target, intent(in) :: elem
    type(given_function) :: var
    integer(ip)          :: nnode,ngaus
    var%ivar  = ivar
    var%icomp = icomp
    var%elem => elem
    var%ndof = prob%vars_of_unk(ivar)
    call var%SetTemp()
  end function given_function_constructor

 function given_function_gradient_constructor(u) result(g)
    type(given_function), intent(in) :: u
    type(given_function_gradient)    :: g
    integer(ip) :: ndime, nnode, ngaus
    call u%GuardTemp()
    call copy_fem_function(u,g)
    call u%CleanTemp()
  end function given_function_gradient_constructor

 function given_function_divergence_constructor(u) result(g)
    type(given_function), intent(in) :: u
    type(given_function_divergence)  :: g
    integer(ip) :: ndime, nnode, ngaus
    call u%GuardTemp()
    call copy_fem_function(u,g)
    call u%CleanTemp()
  end function given_function_divergence_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! given_function interpolations
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine scalar_interpolation (u,res)
    type(given_function)  , intent(in)    :: u
    type(scalar)          , intent(out) :: res
    integer(ip)  :: ndof,nnode,ngaus
    integer(ip)  :: idof,inode,igaus
    call u%GuardTemp()
    assert( u%ndof == 1)
    nnode = u%elem%integ(u%ivar)%p%uint_phy%nnode
    ngaus = u%elem%integ(u%ivar)%p%uint_phy%nlocs
    call memalloc(ngaus,res%a,__FILE__,__LINE__)
    do igaus=1,ngaus
       do inode =1,nnode
          res%a(igaus) = u%elem%integ(u%ivar)%p%uint_phy%shape(inode,igaus) * u%elem%unkno(inode,u%ivar,u%icomp)
       end do
    end do
    call res%SetTemp()
    call u%CleanTemp()
  end subroutine scalar_interpolation

  subroutine vector_interpolation (u,vec)
    type(given_function)  , intent(in)    :: u
    type(vector)          , intent(out) :: vec
    integer(ip)  :: ndof,nnode,ngaus
    integer(ip)  :: idof,inode,igaus
    call u%GuardTemp()
    nnode = u%elem%integ(u%ivar)%p%uint_phy%nnode
    ngaus = u%elem%integ(u%ivar)%p%uint_phy%nlocs
    call memalloc(u%ndof,ngaus,vec%a,__FILE__,__LINE__)
    do igaus=1,ngaus
       do idof=1,u%ndof
          do inode =1,nnode
             vec%a(idof,igaus) = u%elem%integ(u%ivar)%p%uint_phy%shape(inode,igaus) * u%elem%unkno(inode, u%ivar-1+idof, u%icomp)
          end do
       end do
    end do
    call vec%SetTemp()
    call u%CleanTemp()
  end subroutine vector_interpolation

  subroutine scalar_gradient_interpolation(g,vec)
    type(given_function_gradient), intent(in)  :: g
    type(vector)                 , intent(out) :: vec
    integer(ip) :: ndime, nnode, ngaus
    integer(ip) :: idime, inode, igaus
    call g%GuardTemp()
    assert(g%ndof==1)
    assert( associated(g%elem))
    ndime = g%elem%integ(g%ivar)%p%uint_phy%ndime
    ngaus = g%elem%integ(g%ivar)%p%uint_phy%nlocs
    call memalloc(ndime,ngaus,vec%a,__FILE__,__LINE__)
    do igaus=1,ngaus
       do inode =1,nnode
          do idime=1,ndime
             vec%a(idime,igaus) = g%elem%integ(g%ivar)%p%uint_phy%deriv(idime,inode,igaus) * g%elem%unkno(inode,g%ivar, g%icomp)
          end do
       end do
    end do
    call vec%SetTemp()
    call g%CleanTemp()
  end subroutine scalar_gradient_interpolation

  subroutine vector_gradient_interpolation(g,tens)
    type(given_function_gradient), intent(in)  :: g
    type(tensor)                 , intent(out) :: tens
    integer(ip) :: ndof, ndime, nnode, ngaus
    integer(ip) :: idof, idime, inode, igaus
    call g%GuardTemp()
    assert(associated(g%elem))
    ndime = g%elem%integ(g%ivar)%p%uint_phy%ndime
    ndof  = g%ndof
    ngaus = g%elem%integ(g%ivar)%p%uint_phy%nlocs
    call memalloc(ndime,ndof,ngaus,tens%a,__FILE__,__LINE__)
    do igaus=1,ngaus
       do idof=1,ndof
          do inode =1,nnode
             do idime=1,ndime
                tens%a(idime,idof,igaus) = g%elem%integ(g%ivar)%p%uint_phy%deriv(idime,inode,igaus) * g%elem%unkno(inode, g%ivar-1+idof, g%icomp)
             end do
          end do
       end do
    end do
    call tens%SetTemp()
    call g%CleanTemp()
  end subroutine vector_gradient_interpolation

  subroutine divergence_interpolation(u,res)
    type(given_function_divergence), intent(in)  :: u
    type(scalar)                   , intent(out) :: res
    integer(ip) :: ndof, ndime, nnode, ngaus
    integer(ip) :: idof, idime, inode, igaus
    call u%GuardTemp()
    ndime = u%elem%integ(u%ivar)%p%uint_phy%ndime
    assert( u%ndof == ndime)
    assert(associated(u%elem))
    ndof  = u%ndof
    ngaus = u%elem%integ(u%ivar)%p%uint_phy%nlocs
    call memalloc(ngaus,res%a,__FILE__,__LINE__)
    do igaus=1,ngaus
       do idof=1,ndof
          do inode =1,nnode
             do idime=1,ndime
                res%a(igaus) = u%elem%integ(u%ivar)%p%uint_phy%deriv(idime,inode,igaus) * u%elem%unkno(inode, u%ivar-1+idof, u%icomp)
             end do
          end do
       end do
    end do
    call res%SetTemp()
    call u%CleanTemp()
  end subroutine divergence_interpolation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Product by fields
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! The following abstract function could be used for
  ! basis_function, basis_function_gradient and basis_function_divergence
  ! if the last two are defined as derived from basis_functions.
  ! This, however, requires polymorphic allocatable.
  ! function product_field_x_basis_function(field_left,u) result(res)
  !   implicit none
  !   class(field)         , intent(in)  :: field
  !   class(basis_function), intent(in)  :: u
  !   class(basis_function), allocatable :: res
  !   allocate(res,mold=u)
  !   call res%SetTemp()
  !   call copy_fem_function(u,res)
  !   if(allocated(u%left_factor)) then
  !      ! Here we need * and = overloading with corresponding allocation/deallocation.
  !      res%left_factor = field_left * u%left_factor
  !   else
  !      res%left_factor = field_left
  !   end if
  ! end function product_field_x_basis_function

  ! Temporarily I replicate code
  function product_field_x_basis_function(field_left,u) result(res)
    implicit none
    class(field)        , intent(in) :: field_left
    type(basis_function), intent(in) :: u
    type(basis_function)             :: res
    call u%GuardTemp()
    call copy_fem_function(u,res)
    if(allocated(u%left_factor)) then
       ! Here we need * and = overloading with corresponding allocation/deallocation.
       !res%left_factor = field_left * u%left_factor
    else
       !res%left_factor = field_left
    end if
    call u%CleanTemp()
  end function product_field_x_basis_function
  function product_field_x_basis_function_gradient(field_left,u) result(res)
    implicit none
    class(field)        , intent(in) :: field_left
    type(basis_function_gradient), intent(in) :: u
    type(basis_function_gradient)             :: res
    call u%GuardTemp()
    call copy_fem_function(u,res)
    if(allocated(u%left_factor)) then
       ! Here we need * and = overloading with corresponding allocation/deallocation.
       !res%left_factor = field_left * u%left_factor
    else
       !res%left_factor = field_left
    end if
    call u%CleanTemp()
  end function product_field_x_basis_function_gradient
  function product_field_x_basis_function_divergence(field_left,u) result(res)
    implicit none
    class(field)        , intent(in) :: field_left
    type(basis_function_divergence), intent(in) :: u
    type(basis_function_divergence)             :: res
    call u%GuardTemp()
    call copy_fem_function(u,res)
    if(allocated(u%left_factor)) then
       ! Here we need * and = overloading with corresponding allocation/deallocation.
       !res%left_factor = field_left * u%left_factor
    else
       !res%left_factor = field_left
    end if
    call u%CleanTemp()
  end function product_field_x_basis_function_divergence

  function product_basis_function_x_field(u,field_right) result(res)
    implicit none
    class(field)        , intent(in) :: field_right
    type(basis_function), intent(in) :: u
    type(basis_function)             :: res
    call u%GuardTemp()
    call copy_fem_function(u,res)
    if(allocated(u%right_factor)) then
       ! Here we need * and = overloading with corresponding allocation/deallocation.
       !res%right_factor = u%right_factor * field_right
    else
       !res%right_factor = field_right
    end if
    call u%CleanTemp()
  end function product_basis_function_x_field
  function product_basis_function_gradient_x_field(u,field_right) result(res)
    implicit none
    class(field)        , intent(in) :: field_right
    type(basis_function_gradient), intent(in) :: u
    type(basis_function_gradient)             :: res
    call u%GuardTemp()
    call copy_fem_function(u,res)
    if(allocated(u%right_factor)) then
       ! Here we need * and = overloading with corresponding allocation/deallocation.
       !res%right_factor = u%right_factor * field_right
    else
       !res%right_factor = field_right
    end if
    call u%CleanTemp()
  end function product_basis_function_gradient_x_field
  function product_basis_function_divergence_x_field(u,field_right) result(res)
    implicit none
    class(field)        , intent(in) :: field_right
    type(basis_function_divergence), intent(in) :: u
    type(basis_function_divergence)             :: res
    call u%GuardTemp()
    call copy_fem_function(u,res)
    if(allocated(u%right_factor)) then
       ! Here we need * and = overloading with corresponding allocation/deallocation.
       !res%right_factor = u%right_factor * field_right
    else
       !res%right_factor = field_right
    end if
    call u%CleanTemp()
  end function product_basis_function_divergence_x_field

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Scaling (product by reals)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function scale_left_basis_function(alpha,ul) result(x)
    implicit none
    real(rp)            , intent(in) :: alpha
    type(basis_function), intent(in) :: ul
    type(basis_function)             :: x
    call ul%GuardTemp()
    call copy_fem_function(ul,x)
    x%scaling = alpha
    call ul%CleanTemp()
  end function scale_left_basis_function
  function scale_left_basis_function_gradient(alpha,ul) result(x)
    implicit none
    real(rp)                     , intent(in) :: alpha
    type(basis_function_gradient), intent(in) :: ul
    type(basis_function_gradient)             :: x
    call ul%GuardTemp()
    call copy_fem_function(ul,x)
    x%scaling = alpha
    call ul%CleanTemp()
  end function scale_left_basis_function_gradient
  function scale_left_basis_function_divergence(alpha,ul) result(x)
    implicit none
    real(rp)                     , intent(in) :: alpha
    type(basis_function_divergence), intent(in) :: ul
    type(basis_function_divergence)             :: x
    call ul%GuardTemp()
    call copy_fem_function(ul,x)
    x%scaling = alpha
    call ul%CleanTemp()
  end function scale_left_basis_function_divergence

  function scale_right_basis_function(ur,alpha) result(x)
    implicit none
    real(rp)            , intent(in) :: alpha
    type(basis_function), intent(in) :: ur
    type(basis_function)             :: x
    call ur%GuardTemp()
    call copy_fem_function(ur,x)
    x%scaling = alpha
    call ur%CleanTemp()
  end function scale_right_basis_function
  function scale_right_basis_function_gradient(ur,alpha) result(x)
    implicit none
    real(rp)                     , intent(in) :: alpha
    type(basis_function_gradient), intent(in) :: ur
    type(basis_function_gradient)             :: x
    call ur%GuardTemp()
    call copy_fem_function(ur,x)
    x%scaling = alpha
    call ur%CleanTemp()
  end function scale_right_basis_function_gradient
  function scale_right_basis_function_divergence(ur,alpha) result(x)
    implicit none
    real(rp)                       , intent(in) :: alpha
    type(basis_function_divergence), intent(in) :: ur
    type(basis_function_divergence)             :: x
    call ur%GuardTemp()
    call copy_fem_function(ur,x)
    x%scaling = alpha
    call ur%CleanTemp()
  end function scale_right_basis_function_divergence

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Field functions
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function product_vector_tensor(xvector,ytensor) result(zvector)
    implicit none
    type(vector) :: xvector    ! xvector%a(:,ngaus)
    type(tensor) :: ytensor    ! ytensor%a(:,:,ngaus)
    type(vector) :: zvector    ! zvector%a(:,ngaus)
    integer(ip) :: i,j,k, n1,n2,ng
    call xvector%GuardTemp()
    call ytensor%GuardTemp()
    call zvector%SetTemp()
    n1=size(ytensor%a,1)
    n2=size(ytensor%a,2)
    ng=size(ytensor%a,3)
    assert(size(xvector%a,1)==n1)
    assert(size(xvector%a,2)==ng)
    call memalloc(n2,ng,zvector%a,__FILE__,__LINE__)
    ! Change by forall, BLAS, etc...
    do k=1,ng
       do j=1,n2
          zvector%a(j,k) = 0.0_rp
          do i=1,n1
             zvector%a(j,k) = zvector%a(j,k) + xvector%a(i,k)*ytensor%a(i,j,k)
          end do
       end do
    end do
    call xvector%CleanTemp()
    call ytensor%CleanTemp()
  end function product_vector_tensor

  function product_tensor_vector(xtensor,yvector) result(zvector)
    implicit none
    type(tensor) :: xtensor    ! xtensor%a(:,:,ngaus)
    type(vector) :: yvector    ! yvector%a(:,ngaus)
    type(vector) :: zvector    ! zvector%a(:,ngaus)
    integer(ip) :: i,j,k, n1,n2,ng
    call xtensor%GuardTemp()
    call yvector%GuardTemp()
    call zvector%SetTemp()
    n1=size(xtensor%a,1)
    n2=size(xtensor%a,2)
    ng=size(xtensor%a,3)
    assert(size(yvector%a,1)==n1)
    assert(size(yvector%a,2)==ng)
    call memalloc(n1,ng,zvector%a,__FILE__,__LINE__)
    ! Change by forall, BLAS, etc...
    do k=1,ng
       do i=1,n1
          zvector%a(i,k) = 0.0_rp
          do j=1,n2
             zvector%a(i,k) = zvector%a(i,k) + xtensor%a(i,j,k)*yvector%a(j,k)
          end do
       end do
    end do
    call xtensor%CleanTemp()
    call yvector%CleanTemp()
  end function product_tensor_vector

  ! function product_scalar_vector
  ! end function product_scalar_vector
  ! function product_vector_scalar
  ! end function product_vector_scalar

  function product_scalar_scalar(xscalar,yscalar) result(zscalar)
    implicit none
    type(scalar) :: xscalar    ! xscalar(ngaus)
    type(scalar) :: yscalar    ! yscalar(ngaus)
    type(scalar) :: zscalar    ! zscalar(ngaus)
    integer(ip)  :: ng
    call xscalar%GuardTemp()
    call yscalar%GuardTemp()
    call zscalar%SetTemp()
    ng = size(xscalar%a,1)
    assert(size(yscalar%a,1)==ng)
    call memalloc(ng,zscalar%a,__FILE__,__LINE__)
    zscalar%a = xscalar%a * yscalar%a
    call xscalar%CleanTemp()
    call yscalar%CleanTemp()
  end function product_scalar_scalar

! left scaling by real
  function scale_1eft_vector(x,yscalar) result(z)
    implicit none
    real(rp)     :: x
    type(scalar) :: yscalar    ! yscalar(ngaus)
    type(scalar) :: z    ! z(ngaus)
    z%a = x*yscalar%a
  end function scale_1eft_vector

! term by term inverse
  function inv(x) result(z)
    implicit none
    type(scalar) :: x    ! x(ngaus)
    type(scalar) :: z    ! z(ngaus)
    integer(ip)  :: i,j,k
    forall(i=1:size(x%a))
       z%a(i) = 1.0_rp/x%a(i)
    end forall
  end function inv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Integration functions (sum over gauss points plus necessary contractions)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! function integral_function_function(v,u) result(mat)
  !   implicit none
  !   type(basis_function), intent(in) :: v
  !   type(basis_function), intent(in) :: u
  !   real(rp) :: mat(1,1)
  ! end function integral_function_function


!
!   function integral_scalar(K,v,u) result(z)
!     implicit none
!     type(array_rp2) :: x    ! x(ndof,ngaus)
!     type(array_rp2) :: y    ! y(ndof,ngaus)
!     type(array_rp2) :: z    ! z(ndof,ndof)
!     type(element)   :: K    ! element jacobian, etc.
!     integer(ip) :: i,j,k

!     assert(size(x,2)==size(y,2))
!     do i=1,size(x,1)
!        do j=1,size(y,1)
!           z%a(i,j) = 0.0_rp
!           do k=1,size(x,2)
!              z%a(i,j) = z%a(i,j) + K%detjm(k)*K%weight(k)*x%a(i,k)*y%a(j,k)
!           end do
!        end do
!     end do
!   end function integral_scalar
! !
!   function integral_vector(K,v,u) result(z)
!     implicit none
!     type(array_rp2) :: x    ! x(ndime,ndof,ngaus)
!     type(array_rp2) :: y    ! y(ndime,ndof,ngaus)
!     type(array_rp2) :: z    ! z(ndof,ndof)
!     type(element)   :: K    ! element jacobian, etc.
!     integer(ip) :: i,j,k,l

!     assert(size(x,2)==size(y,2))
!     do i=1,size(x,2)
!        do j=1,size(y,2)
!           z%a(i,j) = 0.0_rp
!           do k=1,size(x,3)
!              do l=1,size(x,1)
!                 z%a(i,j) = z%a(i,j) + K%detjm(k)*K%weight(k)*x%a(l,i,k)*y%a(l,j,k)
!              end do
!           end do
!        end do
!     end do
!   end function integral_vector

end module element_tools_names
