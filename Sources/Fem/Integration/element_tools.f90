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
  use memor
  use array_names
  use element_fields_names
  use fem_space_names
  use memory_guard_names
  use problem_names
  implicit none
  private

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
     real(rp) :: scaling=1.0_rp
     class(field), allocatable :: left_factor
     class(field), allocatable :: right_factor
   contains
     procedure, pass(ul) :: scale_left     => scale_left_basis_function
     procedure, pass(ur) :: scale_right    => scale_right_basis_function
     procedure, pass(u)  :: product_left  => product_field_basis_function
     procedure, pass(u)  :: product_right => product_basis_function_field
     generic :: operator(*) => scale_right, scale_left, product_left, product_right
  end type basis_function

  type, extends(basis_function) :: basis_function_gradient
     !real(rp), allocatable :: a(:,:,:)   ! deriv(ndime,nnode,ngaus)
     !real(rp) :: scaling=0.0_rp
     !class(field), allocatable :: left_factor
     !class(field), allocatable :: right_factor
   !contains
     ! procedure, pass(ul) :: scal_left => scal_left_basis_function_gradient
     ! procedure, pass(ur) :: scal_right => scal_right_basis_function_gradient
     ! procedure, pass(u)  :: product_by_scalar => product_scalar_x_basis_function
     ! procedure, pass(u)  :: product_by_vector => product_vector_x_basis_function
     ! procedure, pass(u)  :: product_by_tensor => product_tensor_x_basis_function
     ! generic    :: operator(*) => scal_right, scal_left, product_by_scalar, product_by_vector, product_by_tensor
  end type basis_function_gradient

  type, extends(fem_function) :: basis_function_divergence
     !real(rp), allocatable :: a(:,:,:)   ! deriv(ndime,nnode,ngaus)
     !real(rp) :: scaling=0.0_rp
     !class(field), allocatable :: left_factor
     !class(field), allocatable :: right_factor
   !contains
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
  type, extends(given_function) :: given_function_gradient
     !integer(ip)                :: icomp=1
  end type given_function_gradient
  type, extends(given_function) :: given_function_divergence
     !integer(ip)                :: icomp=1
  end type given_function_divergence

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

  interface integral
     module procedure integral_basis_function_basis_function
  end interface integral

  public :: basis_function, basis_function_gradient, basis_function_divergence
  public :: given_function, given_function_gradient, given_function_divergence
  public :: grad, div, integral, interpolation

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    type(scalar)          , intent(inout) :: res
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
    type(vector)          , intent(inout) :: vec
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
    type(vector)                 , intent(inout) :: vec
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
    type(tensor)                 , intent(inout) :: tens
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
    type(scalar)                   , intent(inout) :: res
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
  function product_field_basis_function(field_left,u) result(res)
    implicit none
    class(field)         , intent(in)  :: field_left
    class(basis_function), intent(in)  :: u
    class(basis_function), allocatable :: res
    allocate(res,mold=u)
    call res%SetTemp()
    call copy_fem_function(u,res)
    if(allocated(u%left_factor)) then
       ! Here we need * and = overloading with corresponding allocation/deallocation.
       res%left_factor = field_left * u%left_factor
    else
       res%left_factor = field_left
    end if
  end function product_field_basis_function

  function product_basis_function_field(u,field_right) result(res)
    implicit none
    class(field)         , intent(in)  :: field_right
    class(basis_function), intent(in)  :: u
    class(basis_function), allocatable :: res
    allocate(res,mold=u)
    call res%SetTemp()
    call copy_fem_function(u,res)
    if(allocated(u%left_factor)) then
       ! Here we need * and = overloading with corresponding allocation/deallocation.
       res%right_factor = u%right_factor * field_right 
    else
       res%right_factor = field_right
    end if
  end function product_basis_function_field

  function scale_left_basis_function(alpha,ul) result(x)
    implicit none
    real(rp)             , intent(in)  :: alpha
    class(basis_function), intent(in)  :: ul
    class(basis_function), allocatable :: x
    call ul%GuardTemp()
    allocate(x,mold=ul)
    call x%SetTemp()
    call copy_fem_function(ul,x)
    x%scaling = x%scaling * alpha
    call ul%CleanTemp()
  end function scale_left_basis_function

  function scale_right_basis_function(ur,alpha) result(x)
    implicit none
    real(rp)             , intent(in)  :: alpha
    class(basis_function), intent(in)  :: ur
    class(basis_function), allocatable :: x
    call ur%GuardTemp()
    allocate(x,mold=ur)
    call x%SetTemp()
    call copy_fem_function(ur,x)
    x%scaling = x%scaling * alpha
    call ur%CleanTemp()
  end function scale_right_basis_function


  ! ! Temporarily I replicate code
  ! function product_field_basis_function(field_left,u) result(res)
  !   implicit none
  !   class(field)        , intent(in) :: field_left
  !   type(basis_function), intent(in) :: u
  !   type(basis_function)             :: res
  !   call u%GuardTemp()
  !   call copy_fem_function(u,res)
  !   if(allocated(u%left_factor)) then
  !      ! Here we need * and = overloading with corresponding allocation/deallocation.
  !      !res%left_factor = field_left * u%left_factor
  !   else
  !      !res%left_factor = field_left
  !   end if
  !   call u%CleanTemp()
  ! end function product_field_basis_function
  ! function product_field_basis_function_gradient(field_left,u) result(res)
  !   implicit none
  !   class(field)        , intent(in) :: field_left
  !   type(basis_function_gradient), intent(in) :: u
  !   type(basis_function_gradient)             :: res
  !   call u%GuardTemp()
  !   call copy_fem_function(u,res)
  !   if(allocated(u%left_factor)) then
  !      ! Here we need * and = overloading with corresponding allocation/deallocation.
  !      !res%left_factor = field_left * u%left_factor
  !   else
  !      !res%left_factor = field_left
  !   end if
  !   call u%CleanTemp()
  ! end function product_field_basis_function_gradient
  ! function product_field_basis_function_divergence(field_left,u) result(res)
  !   implicit none
  !   class(field)        , intent(in) :: field_left
  !   type(basis_function_divergence), intent(in) :: u
  !   type(basis_function_divergence)             :: res
  !   call u%GuardTemp()
  !   call copy_fem_function(u,res)
  !   if(allocated(u%left_factor)) then
  !      ! Here we need * and = overloading with corresponding allocation/deallocation.
  !      !res%left_factor = field_left * u%left_factor
  !   else
  !      !res%left_factor = field_left
  !   end if
  !   call u%CleanTemp()
  ! end function product_field_basis_function_divergence

  ! function product_basis_function_field(u,field_right) result(res)
  !   implicit none
  !   class(field)        , intent(in) :: field_right
  !   type(basis_function), intent(in) :: u
  !   type(basis_function)             :: res
  !   call u%GuardTemp()
  !   call copy_fem_function(u,res)
  !   if(allocated(u%right_factor)) then
  !      ! Here we need * and = overloading with corresponding allocation/deallocation.
  !      !res%right_factor = u%right_factor * field_right
  !   else
  !      !res%right_factor = field_right
  !   end if
  !   call u%CleanTemp()
  ! end function product_basis_function_field
  ! function product_basis_function_gradient_field(u,field_right) result(res)
  !   implicit none
  !   class(field)        , intent(in) :: field_right
  !   type(basis_function_gradient), intent(in) :: u
  !   type(basis_function_gradient)             :: res
  !   call u%GuardTemp()
  !   call copy_fem_function(u,res)
  !   if(allocated(u%right_factor)) then
  !      ! Here we need * and = overloading with corresponding allocation/deallocation.
  !      !res%right_factor = u%right_factor * field_right
  !   else
  !      !res%right_factor = field_right
  !   end if
  !   call u%CleanTemp()
  ! end function product_basis_function_gradient_field
  ! function product_basis_function_divergence_field(u,field_right) result(res)
  !   implicit none
  !   class(field)        , intent(in) :: field_right
  !   type(basis_function_divergence), intent(in) :: u
  !   type(basis_function_divergence)             :: res
  !   call u%GuardTemp()
  !   call copy_fem_function(u,res)
  !   if(allocated(u%right_factor)) then
  !      ! Here we need * and = overloading with corresponding allocation/deallocation.
  !      !res%right_factor = u%right_factor * field_right
  !   else
  !      !res%right_factor = field_right
  !   end if
  !   call u%CleanTemp()
  ! end function product_basis_function_divergence_field

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Scaling (product by reals)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! function scale_left_basis_function(alpha,ul) result(x)
  !   implicit none
  !   real(rp)            , intent(in) :: alpha
  !   type(basis_function), intent(in) :: ul
  !   type(basis_function)             :: x
  !   call ul%GuardTemp()
  !   call copy_fem_function(ul,x)
  !   x%scaling = alpha
  !   call ul%CleanTemp()
  ! end function scale_left_basis_function
  ! function scale_left_basis_function_gradient(alpha,ul) result(x)
  !   implicit none
  !   real(rp)                     , intent(in) :: alpha
  !   type(basis_function_gradient), intent(in) :: ul
  !   type(basis_function_gradient)             :: x
  !   call ul%GuardTemp()
  !   call copy_fem_function(ul,x)
  !   x%scaling = alpha
  !   call ul%CleanTemp()
  ! end function scale_left_basis_function_gradient
  ! function scale_left_basis_function_divergence(alpha,ul) result(x)
  !   implicit none
  !   real(rp)                     , intent(in) :: alpha
  !   type(basis_function_divergence), intent(in) :: ul
  !   type(basis_function_divergence)             :: x
  !   call ul%GuardTemp()
  !   call copy_fem_function(ul,x)
  !   x%scaling = alpha
  !   call ul%CleanTemp()
  ! end function scale_left_basis_function_divergence

  ! function scale_right_basis_function(ur,alpha) result(x)
  !   implicit none
  !   real(rp)            , intent(in) :: alpha
  !   type(basis_function), intent(in) :: ur
  !   type(basis_function)             :: x
  !   call ur%GuardTemp()
  !   call copy_fem_function(ur,x)
  !   x%scaling = alpha
  !   call ur%CleanTemp()
  ! end function scale_right_basis_function
  ! function scale_right_basis_function_gradient(ur,alpha) result(x)
  !   implicit none
  !   real(rp)                     , intent(in) :: alpha
  !   type(basis_function_gradient), intent(in) :: ur
  !   type(basis_function_gradient)             :: x
  !   call ur%GuardTemp()
  !   call copy_fem_function(ur,x)
  !   x%scaling = alpha
  !   call ur%CleanTemp()
  ! end function scale_right_basis_function_gradient
  ! function scale_right_basis_function_divergence(ur,alpha) result(x)
  !   implicit none
  !   real(rp)                       , intent(in) :: alpha
  !   type(basis_function_divergence), intent(in) :: ur
  !   type(basis_function_divergence)             :: x
  !   call ur%GuardTemp()
  !   call copy_fem_function(ur,x)
  !   x%scaling = alpha
  !   call ur%CleanTemp()
  ! end function scale_right_basis_function_divergence

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Field functions
!
! Already implemented:
!
!         real    scalar  vector  tensor
! real              x       x
! scalar   x        x       x
! vector   x        x               x
! tensor                    x
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   function product_real_scalar(xreal,yscalar) result(zscalar)
!     implicit none
!     real(rp)     , intent(in) :: xreal
!     class(scalar), intent(in) :: yscalar    ! yscalar(ngaus)
!     type(scalar)  :: zscalar    ! zscalar(ngaus)
!     call real_scalar(xreal,yscalar,zscalar)
!   end function product_real_scalar
!   function product_scalar_real(xscalar,yreal) result(zscalar)
!     implicit none
!     real(rp)     , intent(in) :: yreal
!     class(scalar), intent(in) :: xscalar    ! yscalar(ngaus)
!     type(scalar)  :: zscalar    ! zscalar(ngaus)
!     call real_scalar(yreal,xscalar,zscalar)
!   end function product_scalar_real
!   subroutine real_scalar(xreal,yscalar,zscalar)
!     implicit none
!     real(rp)     , intent(in) :: xreal
!     class(scalar), intent(in) :: yscalar    ! yscalar(ngaus)
!     type(scalar) :: zscalar    ! zscalar(ngaus)
!     integer(ip)  :: ng
!     call yscalar%GuardTemp()
!     call zscalar%SetTemp()
!     ng = size(yscalar%a,1)
!     call memalloc(ng,zscalar%a,__FILE__,__LINE__)
!     zscalar%a = xreal * yscalar%a ! Intrinsic product
!     call yscalar%CleanTemp()
!   end subroutine real_scalar

!   function product_scalar_scalar(xscalar,yscalar) result(zscalar)
!     implicit none
!     class(scalar), intent(in) :: xscalar    ! xscalar(ngaus)
!     class(scalar), intent(in) :: yscalar    ! yscalar(ngaus)
!     type(scalar)  :: zscalar    ! zscalar(ngaus)
!     integer(ip)   :: ng
!     call xscalar%GuardTemp()
!     call yscalar%GuardTemp()
!     call zscalar%SetTemp()
!     ng = size(xscalar%a,1)
!     assert(size(yscalar%a,1)==ng)
!     call memalloc(ng,zscalar%a,__FILE__,__LINE__)
!     zscalar%a = xscalar%a * yscalar%a
!     call xscalar%CleanTemp()
!     call yscalar%CleanTemp()
!   end function product_scalar_scalar

!   function sum_scalar_scalar(xscalar,yscalar) result(zscalar)
!     implicit none
!     class(scalar), intent(in) :: xscalar    ! xscalar(ngaus)
!     class(scalar), intent(in) :: yscalar    ! yscalar(ngaus)
!     type(scalar)  :: zscalar    ! zscalar(ngaus)
!     integer(ip)   :: ng
!     call xscalar%GuardTemp()
!     call yscalar%GuardTemp()
!     call zscalar%SetTemp()
!     ng = size(xscalar%a,1)
!     assert(size(yscalar%a,1)==ng)
!     call memalloc(ng,zscalar%a,__FILE__,__LINE__)
!     zscalar%a = xscalar%a + yscalar%a
!     call xscalar%CleanTemp()
!     call yscalar%CleanTemp()
!   end function sum_scalar_scalar

!   function product_scalar_vector(xscalar,yvector) result(zvector)
!     implicit none
!     class(scalar), intent(in) :: xscalar    ! xscalar(ngaus)
!     class(vector), intent(in) :: yvector    ! yvector(:,ngaus)
!     type(vector)  :: zvector                ! zvector(:,ngaus)
!     integer(ip)   :: i,j,nd,ng
!     call xscalar%GuardTemp()
!     call yvector%GuardTemp()
!     call zvector%SetTemp()
!     nd = size(yvector%a,1)
!     ng = size(yvector%a,2)
!     call memalloc(nd,ng,zvector%a,__FILE__,__LINE__)
!     do i=1,ng
!        do j=1,nd
!           zvector%a(j,i) = xscalar%a(i) * yvector%a(j,i)
!        end do
!     end do

!     call xscalar%CleanTemp()
!     call yvector%CleanTemp()
!   end function product_scalar_vector

!   function product_vector_scalar(xvector,yscalar) result(zvector)
!     implicit none
!     class(vector), intent(in) :: xvector    ! xvector(:,ngaus)
!     class(scalar), intent(in) :: yscalar    ! yscalar(ngaus)
!     type(vector)  :: zvector                ! zvector(:,ngaus)
!     integer(ip)   :: i,j,nd,ng
!     call xvector%GuardTemp()
!     call yscalar%GuardTemp()
!     call zvector%SetTemp()
!     nd = size(xvector%a,1)
!     ng = size(xvector%a,2)
!     call memalloc(nd,ng,zvector%a,__FILE__,__LINE__)
!     do i=1,ng
!        do j=1,nd
!           zvector%a(j,i) = xvector%a(j,i) * yscalar%a(i)
!        end do
!     end do
!     call xvector%CleanTemp()
!     call yscalar%CleanTemp()
!   end function product_vector_scalar

!   function product_real_vector(xreal,yvector) result(zvector)
!     implicit none
!     real(rp)     , intent(in) :: xreal
!     class(vector), intent(in) :: yvector    ! yvector(ngaus)
!     type(vector)  :: zvector    ! zvector(ngaus)
!     call real_vector(xreal,yvector,zvector)
!   end function product_real_vector
!   function product_vector_real(xvector,yreal) result(zvector)
!     implicit none
!     real(rp)     , intent(in) :: yreal
!     class(vector), intent(in) :: xvector    ! yvector(ngaus)
!     type(vector)  :: zvector    ! zvector(ngaus)
!     call real_vector(yreal,xvector,zvector)
!   end function product_vector_real
!   subroutine real_vector(xreal,yvector,zvector)
!     implicit none
!     real(rp)     , intent(in) :: xreal
!     class(vector), intent(in) :: yvector    ! yvector(ngaus)
!     type(vector) :: zvector    ! zvector(ngaus)
!     integer(ip)  :: nd,ng
!     call yvector%GuardTemp()
!     call zvector%SetTemp()
!     nd = size(yvector%a,1)
!     ng = size(yvector%a,2)
!     call memalloc(nd,ng,zvector%a,__FILE__,__LINE__)
!     zvector%a = xreal * yvector%a ! Intrinsic product
!     call yvector%CleanTemp()
!   end subroutine real_vector

!   function product_vector_tensor(xvector,ytensor) result(zvector)
!     implicit none
!     class(vector), intent(in) :: xvector    ! xvector%a(:,ngaus)
!     type(tensor), intent(in) :: ytensor    ! ytensor%a(:,:,ngaus)
!     type(vector) :: zvector    ! zvector%a(:,ngaus)
!     integer(ip)  :: i,j,k, n1,n2,ng
!     call xvector%GuardTemp()
!     call ytensor%GuardTemp()
!     call zvector%SetTemp()
!     n1=size(ytensor%a,1)
!     n2=size(ytensor%a,2)
!     ng=size(ytensor%a,3)
!     assert(size(xvector%a,1)==n1)
!     assert(size(xvector%a,2)==ng)
!     call memalloc(n2,ng,zvector%a,__FILE__,__LINE__)
!     ! Change by forall, BLAS, etc...
!     do k=1,ng
!        do j=1,n2
!           zvector%a(j,k) = 0.0_rp
!           do i=1,n1
!              zvector%a(j,k) = zvector%a(j,k) + xvector%a(i,k)*ytensor%a(i,j,k)
!           end do
!        end do
!     end do
!     call xvector%CleanTemp()
!     call ytensor%CleanTemp()
!   end function product_vector_tensor

!   function product_tensor_vector(xtensor,yvector) result(zvector)
!     implicit none
!     type(tensor), intent(in) :: xtensor    ! xtensor%a(:,:,ngaus)
!     class(vector), intent(in) :: yvector    ! yvector%a(:,ngaus)
!     type(vector) :: zvector    ! zvector%a(:,ngaus)
!     integer(ip)  :: i,j,k, n1,n2,ng
!     call xtensor%GuardTemp()
!     call yvector%GuardTemp()
!     call zvector%SetTemp()
!     n1=size(xtensor%a,1)
!     n2=size(xtensor%a,2)
!     ng=size(xtensor%a,3)
!     assert(size(yvector%a,1)==n1)
!     assert(size(yvector%a,2)==ng)
!     call memalloc(n1,ng,zvector%a,__FILE__,__LINE__)
!     ! Change by forall, BLAS, etc...
!     do k=1,ng
!        do i=1,n1
!           zvector%a(i,k) = 0.0_rp
!           do j=1,n2
!              zvector%a(i,k) = zvector%a(i,k) + xtensor%a(i,j,k)*yvector%a(j,k)
!           end do
!        end do
!     end do
!     call xtensor%CleanTemp()
!     call yvector%CleanTemp()
!   end function product_tensor_vector

! ! term by term inverse
!   function inv(x) result(z)
!     implicit none
!     type(scalar) :: x    ! x(ngaus)
!     type(scalar) :: z    ! z(ngaus)
!     integer(ip)  :: i,ng
!     call x%GuardTemp()
!     call z%SetTemp()
!     ng = size(x%a,1)
!     call memalloc(ng,z%a,__FILE__,__LINE__)
!     forall(i=1:ng)
!        z%a(i) = 1.0_rp/x%a(i)
!     end forall
!     call x%CleanTemp()
!   end function inv

!   function norm(x) result(z)
!     implicit none
!     type(vector) :: x    ! x(ndime,ngaus)
!     type(scalar) :: z    ! z(ngaus)
!     integer(ip)  :: i,j,k
!     do i=1,size(x%a,2)
!        z%a(i) = 0.0_rp
!        do j=1,size(x%a,1)
!           z%a(i) = z%a(i) + x%a(j,i)
!        end do
!     end do
!   end function norm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Integration functions (sum over gauss points plus necessary contractions)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function integral_basis_function_basis_function(v,u) result(mat)
    implicit none
    type(basis_function), intent(in) :: v
    type(basis_function), intent(in) :: u
    type(array_rp2) :: mat,tmp

    integer(ip)  :: nnode,ngaus
    integer(ip)  :: istart,jstart,ipos,jpos,inode,jnode,ivar,jvar,igaus
    real(rp)     :: factor

    call u%GuardTemp()
    nnode = u%elem%integ(u%ivar)%p%uint_phy%nnode
    ngaus = u%elem%integ(u%ivar)%p%uint_phy%nlocs

    ! Allocate mat with mold=elem%p_mat
    call array_create(u%elem%p_mat%nd1,u%elem%p_mat%nd2,mat)

    ! Starting positions in the element matrix
    istart=0
    do ivar = 1, v%ivar-1
       istart = istart + v%elem%f_inf(ivar)%p%nnode 
    end do
    jstart=0
    do jvar = 1, u%ivar-1
       jstart = jstart + u%elem%f_inf(jvar)%p%nnode 
    end do

    ! Now perform operations according to left_factor and right_factor, currently only left_factor implemented
    if(allocated(v%left_factor)) then
       select type(v_left_factor => v%left_factor)
       class is(scalar)
          if(allocated(u%left_factor)) then
             select type(u_left_factor => u%left_factor)
             class is(scalar)
                ! Diagonal in ivar,ivar contraction of a basis function vector (for ndof=ndime)
                ! Assuming interpolation independent of ivar:
                call array_create(nnode,nnode,tmp)
                do igaus=1,ngaus
                   factor = v%elem%integ(v%ivar)%p%femap%detjm(igaus) * v%elem%integ(v%ivar)%p%quad%weight(igaus) * &
                        &   v_left_factor%a(igaus) * u_left_factor%a(igaus) * v%scaling * u%scaling
                   tmp%a(inode,jnode) = tmp%a(inode,jnode) + factor * &
                        & v%elem%integ(v%ivar)%p%uint_phy%shape(inode,igaus)*u%elem%integ(u%ivar)%p%uint_phy%shape(jnode,igaus)
                end do
                do inode = 1, nnode
                   do jnode = 1, nnode
                      do ivar = v%ivar, v%ivar + v%ndof - 1
                         ipos = istart + (ivar-v%ivar)*nnode + inode
                         jpos = jstart + (ivar-v%ivar)*nnode + jnode
                         mat%a(ipos,jpos) = mat%a(ipos,jpos) + tmp%a(inode,jnode) 
                      end do
                   end do
                end do
                call array_free(tmp)
                ! General case:
                ! do igaus=1,ngaus
                !    ! Quadrature is actually independent of ivar. See comments in fem_space and integration.
                !    factor = v%elem%integ(v%ivar)%p%femap%detjm(igaus)*v%elem%integ(v%ivar)%p%quad%weight(igaus)*v%left_factor(igaus)*u%left_factor(igaus)
                !    ipos = 0
                !    do ivar = v%ivar, v%ivar + v%ndof - 1
                !       do inode = 1, elem%f_inf(ivar)%p%nnode
                !          ipos = ipos + 1
                !          jpos = 0
                !          do jvar = u%ivar, u%ivar + u%ndof - 1
                !             do jnode = 1, elem%f_inf(jvar)%p%nnode
                !                jpos = jpos + 1
                !                if(jvar==ivar) then
                !                   mat%a(istart+ipos,jstart+jpos) = mat%a(istart+ipos,jstart+jpos) + factor* &
                !                     &  v%elem%integ(v%ivar)%p%uint_phy%shape(inode,igaus)*u%elem%integ(u%ivar)%p%uint_phy%shape(jnode,igaus)
                !                end if
                !             end do
                !          end do
                !       end do
                !    end do
                ! end do
             class default
                write(*,*) 'Wrong contraction'
                check(1==0)
             end select
          else
          end if
       class default
          write(*,*) 'Wrong contraction'
          check(1==0)
       end select
    else
    end if

    ! The following lines are needed to implement the non-diagonal case, e.g. a porosity tensor.
    ! do ivar = v%ivar, v%ivar + v%ndof - 1
    !    do inode = 1, elem%f_inf(ivar)%p%nnode
    !       ipos = ipos + 1
    !       do jvar = u%ivar, u%ivar + u%ndof - 1
    !          do jnode = 1, elem%f_inf(jvar)%p%nnode
    !             jpos = jpos + 1
    !             mat%a(ipos,jpos) = mat%a(ipos,jpos) + &
    !                  &  v%elem%integ(u%ivar)%p%uint_phy%shape(inode,igaus) * u%elem%integ(u%ivar)%p%uint_phy%shape(jnode,igaus)
    !          end do
    !       end do
    !    end do
    ! end do

  end function integral_basis_function_basis_function

end module element_tools_names
