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
  use types_names
  use memor_names
  use allocatable_array_names
  use volume_integration_tools_names
  use face_integration_tools_names
  use element_fields_names
  use serial_fe_space_names
  use memory_guard_names
  use problem_names
  implicit none
  private

  type, extends(memory_guard_t) :: function_t
     integer(ip)  :: ivar=1
     integer(ip)  :: nvar=1
     integer(ip)  :: idof=1 ! First dof corresponding to ivar (=sum of nnode for jvar<ivar)
     integer(ip)  :: ndof=1
     integer(ip)  :: ngaus=1
     type(volume_integrator_pointer_t), pointer :: integ(:) => NULL()
   contains
     procedure :: free => free_function
     procedure :: default_initialization => function_default_init
  end type function_t

  ! Basis functions can be left multiplied by a field.
  ! scalar -> shape(inode,igaus)
  ! vector -> shape(idime,jdime,inode,igaus) ! 2nd and 3rd index to to the matrix as one
  ! tensor -> shape(idime,jdime,kdime,ldime,inode,igaus) ! 3rd, 4th and 5th indices go to the matrix.
  type, extends(function_t) :: basis_function_t
     !real(rp), allocatable :: a(:,:)     ! shape(nnode,ngaus)
     real(rp) :: scaling=1.0_rp
     class(field_t), allocatable :: left_factor
   contains
     procedure, pass(ul) :: scale_left     => scale_left_basis_function
     procedure, pass(ur) :: scale_right    => scale_right_basis_function
     procedure, pass(u)  :: product_left  => product_field_basis_function
     !procedure, pass(u)  :: product_right => product_basis_function_field
     generic :: operator(*) => scale_right, scale_left, product_left !, product_right
     procedure :: default_initialization => basis_function_default_init
  end type basis_function_t

  ! scalar -> deriv(idime,inode,igaus)
  ! vector -> deriv(idime,jdime,kdime,inode,igaus)
  ! tensor -> deriv(idime,jdime,kdime,ldime,mdime,inode,igaus)
  type, extends(basis_function_t) :: basis_function_gradient_t
     !real(rp), allocatable :: a(:,:,:)   ! deriv(ndime,nnode,ngaus)
     !real(rp) :: scaling=0.0_rp
     !class(field_t), allocatable :: left_factor
     !class(field_t), allocatable :: right_factor
   !contains
     ! procedure, pass(ul) :: scal_left => scal_left_basis_function_gradient
     ! procedure, pass(ur) :: scal_right => scal_right_basis_function_gradient
     ! procedure, pass(u)  :: product_by_scalar => product_scalar_x_basis_function
     ! procedure, pass(u)  :: product_by_vector => product_vector_x_basis_function
     ! procedure, pass(u)  :: product_by_tensor => product_tensor_x_basis_function
     ! generic    :: operator(*) => scal_right, scal_left, product_by_scalar, product_by_vector, product_by_tensor
  end type basis_function_gradient_t

  type, extends(function_t) :: basis_function_divergence_t
     !real(rp), allocatable :: a(:,:,:)   ! deriv(ndime,nnode,ngaus)
     !real(rp) :: scaling=0.0_rp
     !class(field_t), allocatable :: left_factor
     !class(field_t), allocatable :: right_factor
   !contains
     ! procedure, pass(ul) :: scal_left => scal_left_basis_function_divergence
     ! procedure, pass(ur) :: scal_right => scal_right_basis_function_divergence
     ! procedure, pass(u)  :: product_by_scalar => product_scalar_x_basis_function
     ! procedure, pass(u)  :: product_by_vector => product_vector_x_basis_function
     ! procedure, pass(u)  :: product_by_tensor => product_tensor_x_basis_function
     ! generic    :: operator(*) => scal_right, scal_left, product_by_scalar, product_by_vector, product_by_tensor
  end type basis_function_divergence_t

  ! Interpolations
  ! type, extends(function_t) :: given_function_t
  !    integer(ip)                :: icomp=1
  ! end type given_function_t
  ! type, extends(given_function_t) :: given_function_gradient_t
  !    !integer(ip)                :: icomp=1
  ! end type given_function_gradient_t
  ! type, extends(given_function_t) :: given_function_divergence_t
  !    !integer(ip)                :: icomp=1
  ! end type given_function_divergence_t

  interface basis_function_t
     module procedure basis_function_constructor
  end interface basis_function_t
  ! interface given_function
  !    module procedure given_function_constructor
  ! end interface given_function
  interface grad
     module procedure basis_function_gradient_constructor !, given_function_gradient_constructor
  end interface grad
  interface div
     module procedure basis_function_divergence_constructor !, given_function_divergence_constructor
  end interface div

  interface interpolation
     module procedure scalar_interpolation
     module procedure vector_interpolation
     module procedure scalar_gradient_interpolation
     module procedure vector_gradient_interpolation
     ! module procedure divergence_interpolation
  end interface interpolation

  interface integral
     module procedure integral_basis_function_basis_function
     module procedure integral_basis_function_gradient_basis_function_gradient
  end interface integral

  public :: basis_function_t, basis_function_gradient_t, basis_function_divergence_t
  !public :: given_function_t, given_function_gradient_t, given_function_divergence_t
  public :: grad, div, integral, interpolation

  public :: create_scalar, create_vector, create_tensor

contains


  !=============================================================================
  subroutine function_default_init (this)
    implicit none
    class(function_t), intent(inout) :: this
    this%ivar=1
    this%nvar=1
    this%idof=1 ! First dof corresponding to ivar (=sum of nnode for jvar<ivar)
    this%ndof=1
    this%ngaus=1
    call this%NullifyTemporary()
  end subroutine function_default_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine free_function ( this )
    class(function_t), intent(inout) :: this
  end subroutine free_function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine copy_function(from,to)
    implicit none
    class(function_t), intent(in), target  :: from
    class(function_t), intent(inout)       :: to
    to%ivar=from%ivar
    to%nvar=from%nvar
    to%idof=from%idof
    to%ndof=from%ndof
    to%ngaus=from%ngaus
    to%integ=>from%integ
    call to%SetTemp()
  end subroutine copy_function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! basis_function (some code replication to avoid functions returning polymorphic allocatables)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function basis_function_constructor(prob,iunk,start,integ) result(res)
    implicit none
    class(physical_problem_t), intent(in) :: prob
    integer(ip)            , intent(in) :: iunk
    integer(ip)            , intent(in) :: start(prob%nvars+1)
    type(volume_integrator_pointer_t), target, intent(in) :: integ(:)
    type(basis_function_t) :: res
    integer(ip)          :: i,ivar
    ivar = 1
    do i=1,iunk-1
       ivar = ivar + prob%vars_of_unk(i)
    end do
    res%ivar  = ivar
    res%nvar = prob%vars_of_unk(iunk)

    res%idof  = start(ivar)                ! First dof this variable contributes to
    res%ndof  = start(prob%nvars+1) - 1    ! nvar*nnode for equal interpolation of all variables
    res%ngaus = integ(1)%p%uint_phy%nlocs  ! All variables must have same ngaus

    res%integ => integ
    call res%SetTemp()
  end function basis_function_constructor

 function basis_function_gradient_constructor(u) result(g)
    type(basis_function_t), intent(in) :: u
    type(basis_function_gradient_t)    :: g
    integer(ip) :: ndime, nnode, ngaus
    call u%GuardTemp()
    call copy_function(u,g)
    call u%CleanTemp()
  end function basis_function_gradient_constructor

 function basis_function_divergence_constructor(u) result(g)
    type(basis_function_t), intent(in) :: u
    type(basis_function_divergence_t)    :: g
    integer(ip) :: ndime, nnode, ngaus
    call u%GuardTemp()
    call copy_function(u,g)
    call u%CleanTemp()
  end function basis_function_divergence_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! given_function constructors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !  function given_function_constructor(prob,iunk,icomp,integ) result(var)
 !    implicit none
 !    class(physical_problem_t) , intent(in) :: prob
 !    integer(ip)             , intent(in) :: iunk
 !    integer(ip)             , intent(in) :: icomp
 !    type(volume_integrator_pointer_t),target, intent(in) :: integ(:)
 !    type(given_function_t) :: var
 !    integer(ip)          :: i,ivar
 !    ivar = 1
 !    do i=1,iunk-1
 !       ivar = ivar + prob%vars_of_unk(i)
 !    end do
 !    var%ivar  = ivar
 !    var%icomp = icomp
 !    var%integ => integ
 !    var%nvar = prob%vars_of_unk(iunk)
 !    call var%SetTemp()
 !  end function given_function_constructor

 ! function given_function_gradient_constructor(u) result(g)
 !    type(given_function_t), intent(in) :: u
 !    type(given_function_gradient_t)    :: g
 !    integer(ip) :: ndime, nnode, ngaus
 !    call u%GuardTemp()
 !    call copy_function(u,g)
 !    call u%CleanTemp()
 !  end function given_function_gradient_constructor

 ! function given_function_divergence_constructor(u) result(g)
 !    type(given_function_t), intent(in) :: u
 !    type(given_function_divergence_t)  :: g
 !    integer(ip) :: ndime, nnode, ngaus
 !    call u%GuardTemp()
 !    call copy_function(u,g)
 !    call u%CleanTemp()
 !  end function given_function_divergence_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Interpolations
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine create_scalar (prob, ivar, integ, res)
    implicit none
    class(physical_problem_t)         , intent(in)  :: prob
    integer(ip)                    , intent(in)  :: ivar
    type(volume_integrator_pointer_t), intent(in)  :: integ(:)
    type(scalar_t)                   , intent(out) :: res
    integer(ip)  :: nvar, ngaus
    nvar = prob%vars_of_unk(ivar)
    assert( nvar == 1)
    ngaus = integ(ivar)%p%uint_phy%nlocs
    call memalloc(ngaus,res%a,__FILE__,__LINE__)
    res%a=0.0_rp
  end subroutine create_scalar

  subroutine scalar_interpolation (unkno, ivar, icomp, integ, res)
    real(rp)     , intent(in)    :: unkno(:,:,:)
    integer(ip)  , intent(in)    :: ivar, icomp
    type(volume_integrator_pointer_t), intent(in)  :: integ(:)
    type(scalar_t) , intent(inout) :: res
    integer(ip)  :: nnode,ngaus
    integer(ip)  :: inode,igaus
    nnode = integ(ivar)%p%uint_phy%nnode
    ngaus = integ(ivar)%p%uint_phy%nlocs
    res%a = 0.0_rp
    do igaus=1,ngaus
       do inode =1,nnode
          res%a(igaus) = res%a(igaus) + integ(ivar)%p%uint_phy%shape(inode,igaus) * unkno(inode,ivar,icomp)
       end do
    end do
  end subroutine scalar_interpolation

  subroutine create_vector (prob, ivar, integ, res)
    implicit none
    class(physical_problem_t)        , intent(in)  :: prob
    integer(ip)                    , intent(in)  :: ivar
    type(volume_integrator_pointer_t), intent(in)  :: integ(:)
    type(vector_t)                   , intent(out) :: res
    integer(ip)  :: nvar, ngaus
    nvar = prob%vars_of_unk(ivar)
    ngaus = integ(ivar)%p%uint_phy%nlocs
    call memalloc(nvar,ngaus,res%a,__FILE__,__LINE__)
    res%a=0.0_rp
  end subroutine create_vector

  subroutine vector_interpolation (unkno, ivar, icomp, integ, res)
    real(rp)     , intent(in)    :: unkno(:,:,:)
    integer(ip)  , intent(in)    :: ivar, icomp
    type(volume_integrator_pointer_t), intent(in)  :: integ(:)
    type(vector_t) , intent(inout) :: res
    integer(ip)  :: nvar,nnode,ngaus
    integer(ip)  :: idof,inode,igaus
    nvar  = size(res%a,1)
    nnode = integ(ivar)%p%uint_phy%nnode
    ngaus = integ(ivar)%p%uint_phy%nlocs
    res%a = 0.0_rp
    do igaus=1,ngaus
       do idof=1,nvar
          do inode =1,nnode
             res%a(idof,igaus) = res%a(idof,igaus) + integ(ivar)%p%uint_phy%shape(inode,igaus) * &
                  &              unkno(inode, ivar-1+idof, icomp)
          end do
       end do
    end do
  end subroutine vector_interpolation

  subroutine scalar_gradient_interpolation (unkno, ivar, ndime, icomp, integ, vec)
    real(rp)     , intent(in)    :: unkno(:,:,:)
    integer(ip)  , intent(in)    :: ivar, icomp, ndime
    type(volume_integrator_pointer_t), intent(in)  :: integ(:)
    type(vector_t)                 , intent(inout) :: vec
    integer(ip) :: nnode, ngaus
    integer(ip) :: idime, inode, igaus
    assert(ndime == integ(ivar)%p%uint_phy%ndime)
    nnode = integ(ivar)%p%uint_phy%nnode
    ngaus = integ(ivar)%p%uint_phy%nlocs
    do igaus=1,ngaus
       do inode =1,nnode
          do idime=1,ndime
             vec%a(idime,igaus) = vec%a(idime,igaus)+ integ(ivar)%p%uint_phy%deriv(idime,inode,igaus) * &
                  &               unkno(inode, ivar, icomp)
          end do
       end do
    end do
  end subroutine scalar_gradient_interpolation

  subroutine create_tensor (prob, ivar, ncomp, integ, res)
    implicit none
    class(physical_problem_t)        , intent(in)  :: prob
    integer(ip)                      , intent(in)  :: ivar,ncomp
    type(volume_integrator_pointer_t), intent(in)  :: integ(:)
    type(tensor_t)                   , intent(out) :: res
    integer(ip)  :: nvar, ngaus
    nvar = prob%vars_of_unk(ivar)
    ngaus = integ(ivar)%p%uint_phy%nlocs
    call memalloc(ncomp,nvar,ngaus,res%a,__FILE__,__LINE__)
    res%a=0.0_rp
  end subroutine create_tensor

  subroutine vector_gradient_interpolation(unkno, ivar, icomp, integ, tens)
    real(rp)                         , intent(in)    :: unkno(:,:,:)
    integer(ip)                      , intent(in)    :: ivar, icomp
    type(volume_integrator_pointer_t), intent(in)    :: integ(:)
    type(tensor_t)                   , intent(inout) :: tens
    integer(ip) :: nvar, ndime, nnode, ngaus
    integer(ip) :: idof, idime, inode, igaus
    nvar  = size(tens%a,2)
    ndime = integ(ivar)%p%uint_phy%ndime
    nnode = integ(ivar)%p%uint_phy%nnode
    ngaus = integ(ivar)%p%uint_phy%nlocs
    do igaus=1,ngaus
       do idof=1,nvar
          do inode =1,nnode
             do idime=1,ndime
                tens%a(idime,idof,igaus) = tens%a(idime,idof,igaus) + &
                     integ(ivar)%p%uint_phy%deriv(idime,inode,igaus) * unkno(inode, ivar-1+idof, icomp)
             end do
          end do
       end do
    end do
  end subroutine vector_gradient_interpolation

  ! subroutine scalar_interpolation (u, res)
  !   type(given_function_t)  , intent(in)    :: u
  !   type(scalar_t)          , intent(inout) :: res
  !   integer(ip)  :: nvar,nnode,ngaus
  !   integer(ip)  :: idof,inode,igaus
  !   call u%GuardTemp()
  !   assert( u%nvar == 1)
  !   nnode = u%integ(u%ivar)%p%uint_phy%nnode
  !   ngaus = u%integ(u%ivar)%p%uint_phy%nlocs
  !   call memalloc(ngaus,res%a,__FILE__,__LINE__)
  !   do igaus=1,ngaus
  !      do inode =1,nnode
  !         res%a(igaus) = u%integ(u%ivar)%p%uint_phy%shape(inode,igaus) * u%unkno(inode,u%ivar,u%icomp)
  !      end do
  !   end do
  !   call res%SetTemp()
  !   call u%CleanTemp()
  ! end subroutine scalar_interpolation

  ! subroutine vector_interpolation (u,vec)
  !   type(given_function_t)  , intent(in)    :: u
  !   type(vector_t)          , intent(inout) :: vec
  !   integer(ip)  :: nvar,nnode,ngaus
  !   integer(ip)  :: idof,inode,igaus
  !   call u%GuardTemp()
  !   nnode = u%integ(u%ivar)%p%uint_phy%nnode
  !   ngaus = u%integ(u%ivar)%p%uint_phy%nlocs
  !   call memalloc(u%nvar,ngaus,vec%a,__FILE__,__LINE__)
  !   do igaus=1,ngaus
  !      do idof=1,u%nvar
  !         do inode =1,nnode
  !            vec%a(idof,igaus) = u%integ(u%ivar)%p%uint_phy%shape(inode,igaus) * u%unkno(inode, u%ivar-1+idof, u%icomp)
  !         end do
  !      end do
  !   end do
  !   call vec%SetTemp()
  !   call u%CleanTemp()
  ! end subroutine vector_interpolation

  ! subroutine scalar_gradient_interpolation(g,vec)
  !   type(given_function_gradient_t), intent(in)  :: g
  !   type(vector_t)                 , intent(inout) :: vec
  !   integer(ip) :: ndime, nnode, ngaus
  !   integer(ip) :: idime, inode, igaus
  !   call g%GuardTemp()
  !   assert(g%nvar==1)
  !   assert( associated(g%integ))
  !   ndime = g%integ(g%ivar)%p%uint_phy%ndime
  !   ngaus = g%integ(g%ivar)%p%uint_phy%nlocs
  !   call memalloc(ndime,ngaus,vec%a,__FILE__,__LINE__)
  !   do igaus=1,ngaus
  !      do inode =1,nnode
  !         do idime=1,ndime
  !            vec%a(idime,igaus) = g%integ(g%ivar)%p%uint_phy%deriv(idime,inode,igaus) * g%unkno(inode,g%ivar, g%icomp)
  !         end do
  !      end do
  !   end do
  !   call vec%SetTemp()
  !   call g%CleanTemp()
  ! end subroutine scalar_gradient_interpolation

  ! subroutine vector_gradient_interpolation(g,tens)
  !   type(given_function_gradient_t), intent(in)  :: g
  !   type(tensor_t)                 , intent(inout) :: tens
  !   integer(ip) :: nvar, ndime, nnode, ngaus
  !   integer(ip) :: idof, idime, inode, igaus
  !   call g%GuardTemp()
  !   assert(associated(g%integ))
  !   ndime = g%integ(g%ivar)%p%uint_phy%ndime
  !   nvar  = g%nvar
  !   ngaus = g%integ(g%ivar)%p%uint_phy%nlocs
  !   call memalloc(ndime,nvar,ngaus,tens%a,__FILE__,__LINE__)
  !   do igaus=1,ngaus
  !      do idof=1,nvar
  !         do inode =1,nnode
  !            do idime=1,ndime
  !               tens%a(idime,idof,igaus) = g%integ(g%ivar)%p%uint_phy%deriv(idime,inode,igaus) * g%unkno(inode, g%ivar-1+idof, g%icomp)
  !            end do
  !         end do
  !      end do
  !   end do
  !   call tens%SetTemp()
  !   call g%CleanTemp()
  ! end subroutine vector_gradient_interpolation

  ! subroutine divergence_interpolation(u,res)
  !   type(given_function_divergence_t), intent(in)  :: u
  !   type(scalar_t)                   , intent(inout) :: res
  !   integer(ip) :: nvar, ndime, nnode, ngaus
  !   integer(ip) :: idof, idime, inode, igaus
  !   call u%GuardTemp()
  !   ndime = u%integ(u%ivar)%p%uint_phy%ndime
  !   assert( u%nvar == ndime)
  !   assert(associated(u%integ))
  !   nvar  = u%nvar
  !   ngaus = u%integ(u%ivar)%p%uint_phy%nlocs
  !   call memalloc(ngaus,res%a,__FILE__,__LINE__)
  !   do igaus=1,ngaus
  !      do idof=1,nvar
  !         do inode =1,nnode
  !            do idime=1,ndime
  !               res%a(igaus) = u%integ(u%ivar)%p%uint_phy%deriv(idime,inode,igaus) * u%unkno(inode, u%ivar-1+idof, u%icomp)
  !            end do
  !         end do
  !      end do
  !   end do
  !   call res%SetTemp()
  !   call u%CleanTemp()
  ! end subroutine divergence_interpolation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Product by fields
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! The following abstract function could be used for
  ! basis_function, basis_function_gradient and basis_function_divergence
  ! if the last two are defined as derived from basis_functions.
  ! This, however, requires polymorphic allocatable.

  !=============================================================================
  subroutine basis_function_default_init (this)
    implicit none
    class(basis_function_t), intent(inout) :: this
    this%ivar=1
    this%nvar=1
    this%idof=1 ! First dof corresponding to ivar (=sum of nnode for jvar<ivar)
    this%ndof=1
    this%ngaus=1
    this%scaling=1.0_rp
    call this%NullifyTemporary()
  end subroutine basis_function_default_init

  function product_field_basis_function(field_left,u) result(res)
    implicit none
    class(field_t)         , intent(in)  :: field_left
    class(basis_function_t), intent(in)  :: u
    class(basis_function_t), allocatable :: res
    allocate(res,mold=u); call res%default_initialization()
    call res%SetTemp()
    call copy_function(u,res)
    if(allocated(u%left_factor)) then
       ! Here we need * and = overloading with corresponding allocation/deallocation.
       res%left_factor = field_left * u%left_factor
    else
       res%left_factor = field_left
    end if
  end function product_field_basis_function

  function scale_left_basis_function(alpha,ul) result(x)
    implicit none
    real(rp)             , intent(in)  :: alpha
    class(basis_function_t), intent(in)  :: ul
    class(basis_function_t), allocatable :: x
    call ul%GuardTemp()
    allocate(x,mold=ul); call x%default_initialization()
    call x%SetTemp()
    call copy_function(ul,x)
    x%scaling = x%scaling * alpha
    call ul%CleanTemp()
  end function scale_left_basis_function

  function scale_right_basis_function(ur,alpha) result(x)
    implicit none
    real(rp)             , intent(in)  :: alpha
    class(basis_function_t), intent(in)  :: ur
    class(basis_function_t), allocatable :: x
    call ur%GuardTemp()
    allocate(x,mold=ur); call x%default_initialization()
    call x%SetTemp()
    call copy_function(ur,x)
    x%scaling = x%scaling * alpha
    call ur%CleanTemp()
  end function scale_right_basis_function

  ! ! Temporarily I replicate code
  ! function product_field_basis_function(field_left,u) result(res)
  !   implicit none
  !   class(field_t)        , intent(in) :: field_left
  !   type(basis_function_t), intent(in) :: u
  !   type(basis_function_t)             :: res
  !   call u%GuardTemp()
  !   call res%SetTemp()
  !   call copy_function(u,res)
  !   if(allocated(u%left_factor)) then
  !      ! Here we need * and = overloading with corresponding allocation/deallocation.
  !      res%left_factor = field_left * u%left_factor
  !   else
  !      res%left_factor = field_left
  !   end if
  !   call u%CleanTemp()
  ! end function product_field_basis_function
  ! function product_field_basis_function_gradient(field_left,u) result(res)
  !   implicit none
  !   class(field_t)        , intent(in) :: field_left
  !   type(basis_function_gradient_t), intent(in) :: u
  !   type(basis_function_gradient_t)             :: res
  !   call u%GuardTemp()
  !   call res%SetTemp()
  !   call copy_function(u,res)
  !   if(allocated(u%left_factor)) then
  !      ! Here we need * and = overloading with corresponding allocation/deallocation.
  !      res%left_factor = field_left * u%left_factor
  !   else
  !      res%left_factor = field_left
  !   end if
  !   call u%CleanTemp()
  ! end function product_field_basis_function_gradient
  ! function product_field_basis_function_divergence(field_left,u) result(res)
  !   implicit none
  !   class(field_t)        , intent(in) :: field_left
  !   type(basis_function_divergence_t), intent(in) :: u
  !   type(basis_function_divergence_t)             :: res
  !   call u%GuardTemp()
  !   call res%SetTemp()
  !   call copy_function(u,res)
  !   if(allocated(u%left_factor)) then
  !      ! Here we need * and = overloading with corresponding allocation/deallocation.
  !      res%left_factor = field_left * u%left_factor
  !   else
  !      res%left_factor = field_left
  !   end if
  !   call u%CleanTemp()
  ! end function product_field_basis_function_divergence

  ! function product_basis_function_field(u,field_right) result(res)
  !   implicit none
  !   class(field_t)        , intent(in) :: field_right
  !   type(basis_function_t), intent(in) :: u
  !   type(basis_function_t)             :: res
  !   call u%GuardTemp()
  !   call copy_function(u,res)
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
  !   class(field_t)        , intent(in) :: field_right
  !   type(basis_function_gradient_t), intent(in) :: u
  !   type(basis_function_gradient_t)             :: res
  !   call u%GuardTemp()
  !   call copy_function(u,res)
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
  !   class(field_t)        , intent(in) :: field_right
  !   type(basis_function_divergence_t), intent(in) :: u
  !   type(basis_function_divergence_t)             :: res
  !   call u%GuardTemp()
  !   call copy_function(u,res)
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
  !   type(basis_function_t), intent(in) :: ul
  !   type(basis_function_t)             :: x
  !   call ul%GuardTemp()
  !   call copy_function(ul,x)
  !   x%scaling = alpha
  !   call ul%CleanTemp()
  ! end function scale_left_basis_function
  ! function scale_left_basis_function_gradient(alpha,ul) result(x)
  !   implicit none
  !   real(rp)                     , intent(in) :: alpha
  !   type(basis_function_gradient_t), intent(in) :: ul
  !   type(basis_function_gradient_t)             :: x
  !   call ul%GuardTemp()
  !   call copy_function(ul,x)
  !   x%scaling = alpha
  !   call ul%CleanTemp()
  ! end function scale_left_basis_function_gradient
  ! function scale_left_basis_function_divergence(alpha,ul) result(x)
  !   implicit none
  !   real(rp)                     , intent(in) :: alpha
  !   type(basis_function_divergence_t), intent(in) :: ul
  !   type(basis_function_divergence_t)             :: x
  !   call ul%GuardTemp()
  !   call copy_function(ul,x)
  !   x%scaling = alpha
  !   call ul%CleanTemp()
  ! end function scale_left_basis_function_divergence

  ! function scale_right_basis_function(ur,alpha) result(x)
  !   implicit none
  !   real(rp)            , intent(in) :: alpha
  !   type(basis_function_t), intent(in) :: ur
  !   type(basis_function_t)             :: x
  !   call ur%GuardTemp()
  !   call copy_function(ur,x)
  !   x%scaling = alpha
  !   call ur%CleanTemp()
  ! end function scale_right_basis_function
  ! function scale_right_basis_function_gradient(ur,alpha) result(x)
  !   implicit none
  !   real(rp)                     , intent(in) :: alpha
  !   type(basis_function_gradient_t), intent(in) :: ur
  !   type(basis_function_gradient_t)             :: x
  !   call ur%GuardTemp()
  !   call copy_function(ur,x)
  !   x%scaling = alpha
  !   call ur%CleanTemp()
  ! end function scale_right_basis_function_gradient
  ! function scale_right_basis_function_divergence(ur,alpha) result(x)
  !   implicit none
  !   real(rp)                       , intent(in) :: alpha
  !   type(basis_function_divergence_t), intent(in) :: ur
  !   type(basis_function_divergence_t)             :: x
  !   call ur%GuardTemp()
  !   call copy_function(ur,x)
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
!     class(scalar_t), intent(in) :: yscalar_t    ! yscalar_t(ngaus)
!     type(scalar_t)  :: zscalar_t    ! zscalar_t(ngaus)
!     call real_scalar(xreal,yscalar,zscalar)
!   end function product_real_scalar
!   function product_scalar_real(xscalar,yreal) result(zscalar)
!     implicit none
!     real(rp)     , intent(in) :: yreal
!     class(scalar_t), intent(in) :: xscalar_t    ! yscalar_t(ngaus)
!     type(scalar_t)  :: zscalar_t    ! zscalar_t(ngaus)
!     call real_scalar(yreal,xscalar,zscalar)
!   end function product_scalar_real
!   subroutine real_scalar(xreal,yscalar,zscalar)
!     implicit none
!     real(rp)     , intent(in) :: xreal
!     class(scalar_t), intent(in) :: yscalar_t    ! yscalar_t(ngaus)
!     type(scalar_t) :: zscalar_t    ! zscalar_t(ngaus)
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
!     class(scalar_t), intent(in) :: xscalar_t    ! xscalar_t(ngaus)
!     class(scalar_t), intent(in) :: yscalar_t    ! yscalar_t(ngaus)
!     type(scalar_t)  :: zscalar_t    ! zscalar_t(ngaus)
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
!     class(scalar_t), intent(in) :: xscalar_t    ! xscalar_t(ngaus)
!     class(scalar_t), intent(in) :: yscalar_t    ! yscalar_t(ngaus)
!     type(scalar_t)  :: zscalar_t    ! zscalar_t(ngaus)
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
!     class(scalar_t), intent(in) :: xscalar_t    ! xscalar_t(ngaus)
!     class(vector_t), intent(in) :: yvector_t    ! yvector_t(:,ngaus)
!     type(vector_t)  :: zvector_t                ! zvector_t(:,ngaus)
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
!     class(vector_t), intent(in) :: xvector_t    ! xvector_t(:,ngaus)
!     class(scalar_t), intent(in) :: yscalar_t    ! yscalar_t(ngaus)
!     type(vector_t)  :: zvector_t                ! zvector_t(:,ngaus)
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
!     class(vector_t), intent(in) :: yvector_t    ! yvector_t(ngaus)
!     type(vector_t)  :: zvector_t    ! zvector_t(ngaus)
!     call real_vector(xreal,yvector,zvector)
!   end function product_real_vector
!   function product_vector_real(xvector,yreal) result(zvector)
!     implicit none
!     real(rp)     , intent(in) :: yreal
!     class(vector_t), intent(in) :: xvector_t    ! yvector_t(ngaus)
!     type(vector_t)  :: zvector_t    ! zvector_t(ngaus)
!     call real_vector(yreal,xvector,zvector)
!   end function product_vector_real
!   subroutine real_vector(xreal,yvector,zvector)
!     implicit none
!     real(rp)     , intent(in) :: xreal
!     class(vector_t), intent(in) :: yvector_t    ! yvector_t(ngaus)
!     type(vector_t) :: zvector_t    ! zvector_t(ngaus)
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
!     class(vector_t), intent(in) :: xvector_t    ! xvector_t%a(:,ngaus)
!     type(tensor_t), intent(in) :: ytensor_t    ! ytensor_t%a(:,:,ngaus)
!     type(vector_t) :: zvector_t    ! zvector_t%a(:,ngaus)
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
!     type(tensor_t), intent(in) :: xtensor_t    ! xtensor_t%a(:,:,ngaus)
!     class(vector_t), intent(in) :: yvector_t    ! yvector_t%a(:,ngaus)
!     type(vector_t) :: zvector_t    ! zvector_t%a(:,ngaus)
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
!     type(scalar_t) :: x    ! x(ngaus)
!     type(scalar_t) :: z    ! z(ngaus)
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
!     type(vector_t) :: x    ! x(ndime,ngaus)
!     type(scalar_t) :: z    ! z(ngaus)
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
    type(basis_function_t), intent(in) :: v
    type(basis_function_t), intent(in) :: u
    type(allocatable_array_rp2_t) :: mat

    integer(ip)  :: unode,vnode,ngaus
    integer(ip)  :: ipos,jpos,inode,jnode,ivar,igaus
    real(rp)     :: factor

    call u%GuardTemp()
    call v%GuardTemp()
    ngaus = u%integ(1)%p%uint_phy%nlocs

    ! Allocate mat with mold=p_mat
    assert(v%ndof==u%ndof)
    call allocatable_array_create(v%ndof,u%ndof,mat)

    ! Now perform operations according to left_factor TODO
    if(allocated(v%left_factor)) then
       select type(v_left_factor => v%left_factor)
          class is(scalar_t)
          if(allocated(u%left_factor)) then
             select type(u_left_factor => u%left_factor)
                class is(scalar_t)
                assert(v%nvar==u%nvar)
                do ivar = 0, v%nvar -1
                   vnode = v%integ(v%ivar+ivar)%p%uint_phy%nnode
                   unode = u%integ(u%ivar+ivar)%p%uint_phy%nnode 
                   ipos = v%idof + ivar*vnode
                   jpos = u%idof + ivar*unode
                   do igaus = 1,ngaus
                      factor = v%integ(1)%p%femap%detjm(igaus) * v%integ(1)%p%quad%weight(igaus) * & ! Quadratures independent of u,v and ivar
                           &   v_left_factor%a(igaus) * u_left_factor%a(igaus) * v%scaling * u%scaling
                      do inode = 1, vnode
                         do jnode = 1, unode
                            mat%a(ipos+inode,jpos+jnode) = mat%a(ipos+inode,jpos+jnode) + factor * &
                                 & v%integ(v%ivar+ivar)%p%uint_phy%shape(inode,igaus) * &
                                 & u%integ(u%ivar+ivar)%p%uint_phy%shape(jnode,igaus)
                         end do
                      end do
                   end do
                end do
                class default
                check(.false.)
             end select
          else
          end if
          class default
          check(.false.)
       end select
    else
       assert(v%nvar==u%nvar)
       do ivar = 0, v%nvar -1
          vnode = v%integ(v%ivar+ivar)%p%uint_phy%nnode
          unode = u%integ(u%ivar+ivar)%p%uint_phy%nnode 
          ipos = v%idof + ivar*vnode
          jpos = u%idof + ivar*unode
          do igaus = 1,ngaus
             factor = v%integ(1)%p%femap%detjm(igaus)*v%integ(1)%p%quad%weight(igaus)*v%scaling*u%scaling
             do inode = 1, vnode
                do jnode = 1, unode
                   mat%a(ipos+inode,jpos+jnode) = mat%a(ipos+inode,jpos+jnode) + factor * &
                        & v%integ(v%ivar+ivar)%p%uint_phy%shape(inode,igaus) * &
                        & u%integ(u%ivar+ivar)%p%uint_phy%shape(jnode,igaus)
                end do
             end do
          end do
       end do
    end if
    call u%CleanTemp()
    call v%CleanTemp()

  end function integral_basis_function_basis_function

  function integral_basis_function_gradient_basis_function_gradient(gv,gu) result(mat)
    implicit none
    type(basis_function_gradient_t), intent(in) :: gv
    type(basis_function_gradient_t), intent(in) :: gu
    type(allocatable_array_rp2_t) :: mat

    integer(ip)  :: unode,vnode,ndime,ngaus
    integer(ip)  :: ipos,jpos,inode,jnode,idime,ivar,igaus
    real(rp)     :: factor

    call gu%GuardTemp()
    call gv%GuardTemp()
    ngaus = gu%integ(1)%p%uint_phy%nlocs
    ndime = gu%integ(1)%p%uint_phy%ndime

    ! Allocate mat with mold=p_mat
    assert(gv%ndof==gu%ndof)
    call allocatable_array_create(gv%ndof,gu%ndof,mat)

    ! Now perform operations according to left_factor TODO
    if(allocated(gv%left_factor)) then
       select type(v_left_factor => gv%left_factor)
          class is(scalar_t)
          if(allocated(gu%left_factor)) then
             select type(u_left_factor => gu%left_factor)
                class is(scalar_t)
                assert(gv%nvar==gu%nvar)
                do ivar = 0, gv%nvar -1
                   vnode = gv%integ(gv%ivar+ivar)%p%uint_phy%nnode
                   unode = gu%integ(gu%ivar+ivar)%p%uint_phy%nnode 
                   ipos = gv%idof + ivar*vnode
                   jpos = gu%idof + ivar*unode
                   do igaus = 1,ngaus
                      factor = gv%integ(1)%p%femap%detjm(igaus) * gv%integ(1)%p%quad%weight(igaus) * & ! Quadratures independent of u,v and ivar
                           &   v_left_factor%a(igaus) * u_left_factor%a(igaus) * gv%scaling * gu%scaling
                      do inode = 1, vnode
                         do jnode = 1, unode
                            do idime = 1,ndime
                               mat%a(ipos+inode,jpos+jnode) = mat%a(ipos+inode,jpos+jnode) + factor * &
                                    & gv%integ(gv%ivar+ivar)%p%uint_phy%deriv(idime,inode,igaus) * &
                                    & gu%integ(gu%ivar+ivar)%p%uint_phy%deriv(idime,jnode,igaus)
                            end do
                         end do
                      end do
                   end do
                end do
                class default
                check(.false.)
             end select
          else
          end if
          class default
          check(.false.)
       end select
    else
       assert(gv%nvar==gu%nvar)
       do ivar = 0, gv%nvar -1
          vnode = gv%integ(gv%ivar+ivar)%p%uint_phy%nnode
          unode = gu%integ(gu%ivar+ivar)%p%uint_phy%nnode 
          ipos = gv%idof + ivar*vnode
          jpos = gu%idof + ivar*unode
          do igaus = 1,ngaus
             factor = gv%integ(1)%p%femap%detjm(igaus)*gv%integ(1)%p%quad%weight(igaus)*gv%scaling*gu%scaling
             do inode = 1, vnode
                do jnode = 1, unode
                   do idime = 1,ndime
                      mat%a(ipos+inode,jpos+jnode) = mat%a(ipos+inode,jpos+jnode) + factor * &
                           & gv%integ(gv%ivar+ivar)%p%uint_phy%deriv(idime,inode,igaus) * &
                           & gu%integ(gu%ivar+ivar)%p%uint_phy%deriv(idime,jnode,igaus)
                   end do
                end do
             end do
          end do
       end do
    end if
    call gu%CleanTemp()
    call gv%CleanTemp()

  end function integral_basis_function_gradient_basis_function_gradient


end module element_tools_names
