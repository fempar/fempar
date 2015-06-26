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
module par_vector_names
  ! Serial modules
use types_names
use memor_names
use stdio_names
  use fem_vector_names
use map_apply_names

#ifdef ENABLE_BLAS       
use blas77_interfaces_names
#endif

#ifdef memcheck       
use iso_c_binding
#endif

  ! Parallel modules
  use par_context_names
  use par_environment_names
  use dof_distribution_names
use psb_penv_mod_names

  ! Abstract types
  use base_operand_names

# include "debug.i90"

  implicit none
  private

  !=============================================================
  ! TODO:
  ! 
  ! x Call to BLAS double precision or single precision 
  !   subroutines depending on the value of the rp parameter. 
  !   Currently we are always calling double precision variants 
  !   of the BLAS subroutines.
  ! 
  !=============================================================

  integer(ip), parameter :: undefined     = -1 ! Undefined. State to be assigned by the user or externally.
  integer(ip), parameter :: part_summed   = 0  ! partially summed element-based vector
  integer(ip), parameter :: full_summed   = 1  ! fully     summed element-based vector

  ! Distributed Vector
  type, extends(base_operand_t) :: par_vector_t
     ! Data structure which stores the local part 
     ! of the vector mapped to the current processor.
     ! This is required for both eb and vb data 
     ! distributions
     type( fem_vector_t ) :: f_vector

     ! Partially or fully summed
     integer(ip)  :: state 

     ! Parallel DoF distribution control info.
     type ( dof_distribution_t ), pointer  :: dof_dist => NULL()

     type ( par_environment_t ) , pointer  :: p_env => NULL()
   contains
     ! Provide type bound procedures (tbp) implementors
     procedure :: dot  => par_vector_dot_tbp
     procedure :: copy => par_vector_copy_tbp
     procedure :: init => par_vector_init_tbp
     procedure :: scal => par_vector_scal_tbp
     procedure :: axpby => par_vector_axpby_tbp
     procedure :: nrm2 => par_vector_nrm2_tbp
     procedure :: clone => par_vector_clone_tbp
     procedure :: comm  => par_vector_comm_tbp
     procedure :: free  => par_vector_free_tbp
  end type par_vector_t


  ! Types
  public :: par_vector_t

  ! Constants 
  public :: undefined, part_summed, full_summed

  ! Functions
  public :: par_vector_alloc,    par_vector_free,     par_vector_create_view, & 
       &  par_vector_clone,                                                 &
       &  par_vector_comm,     par_vector_weight,   par_vector_nrm2,        &
       &  par_vector_dot,      par_vector_copy,     par_vector_zero,        &
       &  par_vector_init,                                                  &
       &  par_vector_scale,    par_vector_mxpy,     par_vector_axpy,        &
       &  par_vector_aypx,     par_vector_pxpy,     par_vector_pxmy ,       &
       &  par_vector_print,    par_vector_print_matrix_market

!***********************************************************************
! Allocatable arrays of type(fem_vector_t)
!***********************************************************************
# define var_attr allocatable, target
# define point(a,b) call move_alloc(a,b)
# define generic_status_test             allocated
# define generic_memalloc_interface      memalloc
# define generic_memrealloc_interface    memrealloc
# define generic_memfree_interface       memfree
# define generic_memmovealloc_interface  memmovealloc

# define var_type type(par_vector_t)
# define var_size 52
# define bound_kind ip
# include "mem_header.i90"

  public :: memalloc,  memrealloc,  memfree, memmovealloc

contains

# include "mem_body.i90"

  !=============================================================================
  subroutine par_vector_alloc (dof_dist, p_env, p_vec)
    implicit none
    ! Parameters
    type(dof_distribution_t), target, intent(in)  :: dof_dist
    type(par_environment_t) , target, intent(in)  :: p_env
    type(par_vector_t)              , intent(out) :: p_vec

    ! Locals
    integer(ip)                               :: idof, id

    ! p_env%p_context is required within this subroutine
    assert ( associated(p_env%p_context) )
    assert ( p_env%p_context%created .eqv. .true.)

    p_vec%dof_dist => dof_dist
    p_vec%p_env    => p_env 
    if(p_env%p_context%iam<0) return
    call fem_vector_alloc ( dof_dist%nl, p_vec%f_vector )
    p_vec%state = undefined
  end subroutine par_vector_alloc

  !=============================================================================
  subroutine par_vector_free(p_vec)
    implicit none
    type(par_vector_t), intent(inout) :: p_vec

    ! The routine requires the partition/context info
    assert ( associated( p_vec%dof_dist ) )
    assert ( associated( p_vec%p_env%p_context ) )
    assert ( p_vec%p_env%p_context%created .eqv. .true.)
    if(p_vec%p_env%p_context%iam<0) return

    p_vec%state = undefined

    ! Free local part
    call fem_vector_free ( p_vec%f_vector )

    nullify ( p_vec%dof_dist )
    nullify ( p_vec%p_env )
  end subroutine par_vector_free

  !=============================================================================
  subroutine par_vector_create_view (s_p_vec, start, end, t_p_vec)
    implicit none
    type(par_vector_t), intent(in), target :: s_p_vec
    integer(ip)     , intent(in)         :: start
    integer(ip)     , intent(in)         :: end
    type(par_vector_t), intent(out)        :: t_p_vec

    ! The routine requires the partition/context info
    assert ( associated( s_p_vec%dof_dist ) )
    assert ( associated( s_p_vec%p_env%p_context ) )
    assert ( s_p_vec%p_env%p_context%created .eqv. .true.)

    ! Associate dof distribution and parallel environment 
    t_p_vec%dof_dist => s_p_vec%dof_dist
    t_p_vec%p_env    => s_p_vec%p_env

    assert ( s_p_vec%state /= undefined )
    t_p_vec%state = s_p_vec%state
    
    if(s_p_vec%p_env%p_context%iam<0) return

    ! Call fem_vector_create_view
    call fem_vector_create_view ( s_p_vec%f_vector, start, end, t_p_vec%f_vector ) 

  end subroutine par_vector_create_view

  !=============================================================================
  subroutine par_vector_clone (s_p_vec, t_p_vec)
    implicit none
    type(par_vector_t), intent(in), target :: s_p_vec
    type(par_vector_t), intent(out)        :: t_p_vec 

    ! p_env%p_context is required within this subroutine 
    assert ( associated(s_p_vec%dof_dist) )
    assert ( associated(s_p_vec%p_env%p_context) )
    assert ( s_p_vec%p_env%p_context%created .eqv. .true.)

    t_p_vec%dof_dist => s_p_vec%dof_dist
    t_p_vec%p_env    => s_p_vec%p_env
    t_p_vec%state = s_p_vec%state

    if(s_p_vec%p_env%p_context%iam<0) return

    call fem_vector_clone ( s_p_vec%f_vector, t_p_vec%f_vector )

  end subroutine par_vector_clone

  subroutine par_vector_comm ( p_vec )
    implicit none
    ! Parameters
    type(par_vector_t), intent(inout) :: p_vec

    ! Local variables
    integer(ip)      :: ni 
    type(par_vector_t) :: p_vec_G

    ! Pointer to part/context object is required
    assert ( associated(p_vec%dof_dist) )
    assert ( associated(p_vec%p_env%p_context) )
    assert ( p_vec%p_env%p_context%created .eqv. .true.)
    if(p_vec%p_env%p_context%iam<0) return


    ni = p_vec%f_vector%neq - p_vec%dof_dist%nb
    call par_vector_create_view ( p_vec, ni+1, p_vec%f_vector%neq, p_vec_G )
    call comm_interface ( p_vec_G )

    p_vec%state = full_summed
  end subroutine par_vector_comm

  !=============================================================================
  ! VERY IMPORTANT: comm_interface is well-defined if and only if p_vec is a 
  !                 vector on the interface
  !=============================================================================
  subroutine comm_interface (p_vec)
use par_sparse_global_collectives_names
    implicit none

    ! Parameters
    type(par_vector_t), intent(inout)         :: p_vec

    ! Pointer to part/context object is required
    assert ( associated(p_vec%dof_dist) )
    assert ( associated(p_vec%p_env%p_context) )
    assert ( p_vec%p_env%p_context%created .eqv. .true.)
    if(p_vec%p_env%p_context%iam<0) return


    ! call fem_import_print (6, p_vec%dof_dist%dof_import)
    ! write(*,*) 'SE 1', size(p_vec%f_vector%b)
    ! call fem_vector_print ( 6, p_vec%f_vector )

    ! First stage: owners receive/reduce, non-owners send
    call single_exchange ( p_vec%p_env%p_context%icontxt, &
         p_vec%dof_dist%dof_import%num_rcv,    &
         p_vec%dof_dist%dof_import%list_rcv,   &
         p_vec%dof_dist%dof_import%rcv_ptrs,   &
         p_vec%dof_dist%dof_import%unpack_idx, &
         p_vec%dof_dist%dof_import%num_snd,    &
         p_vec%dof_dist%dof_import%list_snd,   &
         p_vec%dof_dist%dof_import%snd_ptrs,   &
         p_vec%dof_dist%dof_import%pack_idx,   &
         1.0_rp,                         &
         1.0_rp,                         &
         p_vec%f_vector%b ) 

    ! write(*,*) 'SE 2'

    ! Second stage: owners send, non-owners receive/insert
    call single_exchange ( p_vec%p_env%p_context%icontxt, &
         p_vec%dof_dist%dof_import%num_snd,    &
         p_vec%dof_dist%dof_import%list_snd,   &
         p_vec%dof_dist%dof_import%snd_ptrs,   &
         p_vec%dof_dist%dof_import%pack_idx,   &
         p_vec%dof_dist%dof_import%num_rcv,    &
         p_vec%dof_dist%dof_import%list_rcv,   &
         p_vec%dof_dist%dof_import%rcv_ptrs,   &
         p_vec%dof_dist%dof_import%unpack_idx, &
         1.0_rp,                         &
         0.0_rp,                         &
         p_vec%f_vector%b )

  end subroutine comm_interface

  subroutine par_vector_weight ( p_vec, weight )
    implicit none
    ! Parameters
    type(par_vector_t), intent(inout)                :: p_vec
    real(rp)        , intent(in), target, optional :: weight(*)

    ! Local variables
    integer(ip)      :: ni 
    type(par_vector_t) :: p_vec_G

    ! Pointer to part/context object is required
    assert ( associated(p_vec%dof_dist) )
    assert ( associated(p_vec%p_env%p_context) )
    assert ( p_vec%p_env%p_context%created .eqv. .true.)
    if(p_vec%p_env%p_context%iam<0) return


    ni = p_vec%f_vector%neq - p_vec%dof_dist%nb
    call par_vector_create_view ( p_vec, ni+1, p_vec%f_vector%neq, p_vec_G )
    call weight_interface ( p_vec_G, weight )

    p_vec%state = part_summed
  end subroutine par_vector_weight

  !=============================================================================
  ! VERY IMPORTANT: weight_interface is well-defined if and only if p_vec is a 
  !                 vector on the interface
  !=============================================================================
  subroutine weight_interface ( p_vec, weight )
    implicit none
    ! Parameters
    type(par_vector_t), intent(inout)                :: p_vec
    real(rp)        , intent(in), target, optional :: weight(*)

    ! Local variables
    integer(ip) :: iobj, i, i1, i2 
    real(rp)    :: weigt

    ! Pointer to part/context object is required
    assert ( associated(p_vec%dof_dist) )
    assert ( associated(p_vec%p_env%p_context) )
    assert ( p_vec%p_env%p_context%created .eqv. .true.)
    if(p_vec%p_env%p_context%iam<0) return


    if ( present(weight) ) then
       i1 = 1
       i2 = p_vec%dof_dist%nb
       do i=i1,i2
          p_vec%f_vector%b(i) =  p_vec%f_vector%b(i) * weight(i)
       end do
    else 

       do iobj=2, p_vec%dof_dist%nobjs
          weigt=1.0_rp/real(p_vec%dof_dist%lobjs(4,iobj))
          i1 = p_vec%dof_dist%lobjs(2,iobj) - p_vec%dof_dist%ni
          i2 = p_vec%dof_dist%lobjs(3,iobj) - p_vec%dof_dist%ni
          ! write (*,*) 'yyy', i1, i2, (i2-i1+1), weigt
#ifdef ENABLE_BLAS
          call dscal ( (i2-i1+1), weigt, p_vec%f_vector%b(i1:i2), 1 )
#else
          p_vec%f_vector%b(i1:i2) = weigt * p_vec%f_vector%b(i1:i2)
#endif 
       end do

    end if
  end subroutine weight_interface

  subroutine par_vector_nrm2 (x, t)
    implicit none
    type(par_vector_t)   , intent(in), target  :: x
    real(rp)           , intent(out)         :: t
    integer(ip)                           :: ierrc

    ! Pointer to part/context object is required
    assert ( associated(x%dof_dist) )
    assert ( associated(x%p_env%p_context) )
    assert ( x%p_env%p_context%created .eqv. .true.)
    if(x%p_env%p_context%iam<0) return

    assert ( x%state /= undefined )
    ! write (*,*) 'PPP'
    ! call fem_vector_print ( 6, x%f_vector )
    call par_vector_dot (x,x,t)
    t = sqrt(t)
  end subroutine par_vector_nrm2

  !=============================================================================
  ! VERY IMPORTANT: dot_interface is well-defined if and only if x/y are 
  !                 vectors on the interface
  !============================================================================= 
  subroutine dot_interface(x,y,t)
    implicit none
    ! Parameters  
    type(par_vector_t), intent(in)  :: x
    type(par_vector_t), intent(in)  :: y
    real(rp)    , intent(out)     :: t

    ! Locals 
    integer(ip)                :: ierrc
    type(par_vector_t)              :: ws_vec

    ! Pointer to part/context object is required
    assert ( associated(x%dof_dist   ) )
    assert ( associated(x%p_env%p_context) )
    assert ( associated(y%dof_dist   ) )
    assert ( associated(y%p_env%p_context) )

    assert ( x%state /= undefined .and. y%state /= undefined  )


    if ( (x%state == part_summed .and. y%state == full_summed) .or. (x%state == full_summed .and. y%state == part_summed) ) then
       ! Perform local dot products
       call fem_vector_dot (x%f_vector, y%f_vector, t)
    else if ( (x%state == full_summed .and. y%state == full_summed) ) then
       ! Perform local weighted dot products
       call weighted_dot (x, y, t)
    else if ( (x%state == part_summed .and. y%state == part_summed) ) then

       ! Allocate space for ws_vec
       call par_vector_clone ( x, ws_vec )

       ! ws_vec <- x  
       call par_vector_copy  ( x, ws_vec )

       ! Transform ws_vec from partially summed 
       ! to fully summed
       call comm_interface ( ws_vec )
       

       call fem_vector_dot ( x%f_vector, ws_vec%f_vector, t )

       call par_vector_free ( ws_vec )
    end if


  end subroutine dot_interface

  subroutine par_vector_dot (x,y,t)
    implicit none
    ! Parameters  
    type(par_vector_t), intent(in)  :: x
    type(par_vector_t), intent(in)  :: y
    real(rp)        , intent(out) :: t

    ! Locals 
    integer(ip)                :: ierrc
    integer(ip)                   :: ni 
    type(par_vector_t)              :: x_I, x_G, y_I, y_G
    real(rp)                      :: s

    ! Pointer to part/context object is required
    assert ( associated(x%dof_dist   ) )
    assert ( associated(x%p_env%p_context) )
    assert ( associated(y%dof_dist   ) )
    assert ( associated(y%p_env%p_context) )
    assert ( x%p_env%p_context%created .eqv. .true.)
    if(x%p_env%p_context%iam<0) return

    assert ( x%state /= undefined .and. y%state /= undefined  )
    
    ni = x%f_vector%neq - x%dof_dist%nb
    if ( ni > 0 ) then
       call par_vector_create_view ( x, 1, ni, x_I )
       call par_vector_create_view ( y, 1, ni, y_I )
       call fem_vector_dot         ( x_I%f_vector, y_I%f_vector, t )
    else
       t = 0.0_rp
    end if
    
    call par_vector_create_view ( x, ni+1, x%f_vector%neq, x_G )
    call par_vector_create_view ( y, ni+1, y%f_vector%neq, y_G )
    call dot_interface          ( x_G, y_G, s )
    
    t = t + s
    
    ! Reduce-sum local dot products on all processes
    call psb_sum ( y%p_env%p_context%icontxt, t )
  end subroutine par_vector_dot

  !=============================================================================
  subroutine par_vector_copy(x,y)
    implicit none
    type(par_vector_t), intent(in)    :: x
    type(par_vector_t), intent(inout) :: y
    integer(ip)                  :: ierrc

    ! Pointer to part/context object is required
    assert ( associated(x%dof_dist   ) )
    assert ( associated(x%p_env%p_context) )
    assert ( associated(y%dof_dist   ) )
    assert ( associated(y%p_env%p_context) )
    assert ( x%p_env%p_context%created .eqv. .true.)
    if(x%p_env%p_context%iam<0) return

    assert ( x%state /= undefined  )
    ! Perform local copy
    call fem_vector_copy ( x%f_vector, y%f_vector )
    y%state = x%state
  end subroutine par_vector_copy

  !=============================================================================
  subroutine par_vector_zero(y)
    implicit none
    type(par_vector_t), intent(inout) :: y
    integer(ip)                  :: ierrc

    ! Pointer to part/context object is required
    assert ( associated(y%dof_dist   ) )
    assert ( associated(y%p_env%p_context) )
    assert ( y%p_env%p_context%created .eqv. .true.)
    if(y%p_env%p_context%iam<0) return


    ! write(*,*) 'XXX'
    assert ( y%state /= undefined  )
    ! Zero-out local copy
    call fem_vector_zero ( y%f_vector )
    
  end subroutine par_vector_zero

  !=============================================================================
  subroutine par_vector_init (t,y)
    implicit none
    real(rp)        , intent(in)     :: t
    type(par_vector_t), intent(inout)  :: y
    integer(ip)                   :: ierrc

    ! Pointer to part/context object is required
    assert ( associated(y%dof_dist   ) )
    assert ( associated(y%p_env%p_context) )
    assert ( y%p_env%p_context%created .eqv. .true.)
    if(y%p_env%p_context%iam<0) return

    if  ( t == 0.0_rp ) then
       call par_vector_zero(y)
    else 
       ! Scale local copy
       call fem_vector_init (t, y%f_vector)
       ! y%state = full_summed
    end if
  end subroutine par_vector_init


  !=============================================================================
  subroutine par_vector_scale(t,x,y)
    implicit none
    real(rp)    , intent(in)         :: t
    type(par_vector_t), intent(in)     :: x
    type(par_vector_t), intent(inout)  :: y
    integer(ip)                   :: ierrc

    ! Pointer to part/context object is required
    assert ( associated(x%dof_dist   ) )
    assert ( associated(x%p_env%p_context) )
    assert ( associated(y%dof_dist   ) )
    assert ( associated(y%p_env%p_context) )
    assert ( y%p_env%p_context%created .eqv. .true.)
    if(y%p_env%p_context%iam<0) return

    ! Check matching partition/handler
    assert ( x%state /= undefined )

    ! Scale local copy
    call fem_vector_scale (t, x%f_vector, y%f_vector )
    y%state = x%state
  end subroutine par_vector_scale

  !=============================================================================
  subroutine par_vector_mxpy(x,y)
    implicit none
    type(par_vector_t), intent(in)    :: x
    type(par_vector_t), intent(inout) :: y
    integer(ip)                  :: ierrc

    ! Pointer to part/context object is required
    assert ( associated(x%dof_dist   ) )
    assert ( associated(x%p_env%p_context) )
    assert ( associated(y%dof_dist   ) )
    assert ( associated(y%p_env%p_context) )
    assert ( y%p_env%p_context%created .eqv. .true.)
    if(y%p_env%p_context%iam<0) return

    ! Check matching partition/handler
    assert ( x%state /= undefined .and. y%state /= undefined  )
    assert ( x%state == y%state )

    call fem_vector_mxpy (x%f_vector, y%f_vector )

  end subroutine par_vector_mxpy

  !=============================================================================
  subroutine par_vector_axpy(t,x,y)
    implicit none
    real(rp)   , intent(in)         :: t
    type(par_vector_t), intent(in)    :: x
    type(par_vector_t), intent(inout) :: y
    integer(ip)                  :: ierrc

    ! Pointer to part/context object is required
    assert ( associated(x%dof_dist   ) )
    assert ( associated(x%p_env%p_context) )
    assert ( associated(y%dof_dist   ) )
    assert ( associated(y%p_env%p_context) )
    assert ( y%p_env%p_context%created .eqv. .true.)
    if(y%p_env%p_context%iam<0) return

    ! Check matching partition/handler
    assert ( x%state /= undefined .and. y%state /= undefined  )
    assert ( x%state == y%state )

    call fem_vector_axpy  ( t, x%f_vector, y%f_vector )
  end subroutine par_vector_axpy

  !=============================================================================
  subroutine par_vector_aypx(t,x,y)
    implicit none
    real(rp)    , intent(in)        :: t
    type(par_vector_t), intent(in)    :: x
    type(par_vector_t), intent(inout) :: y
    integer(ip)                  :: ierrc

    ! Pointer to part/context object is required
    assert ( associated(x%dof_dist   ) )
    assert ( associated(x%p_env%p_context) )
    assert ( associated(y%dof_dist   ) )
    assert ( associated(y%p_env%p_context) )
    assert ( y%p_env%p_context%created .eqv. .true.)
    if(y%p_env%p_context%iam<0) return

    ! Check matching partition/handler
    assert ( x%state /= undefined .and. y%state /= undefined  )
    assert ( x%state == y%state )
    call fem_vector_aypx (t, x%f_vector, y%f_vector )
  end subroutine par_vector_aypx

  !=============================================================================
  subroutine par_vector_pxpy(x,y)
    implicit none
    type(par_vector_t), intent(in)    :: x
    type(par_vector_t), intent(inout) :: y
    integer(ip)                  :: ierrc

    ! Pointer to part/context object is required
    assert ( associated(x%dof_dist   ) )
    assert ( associated(x%p_env%p_context) )
    assert ( associated(y%dof_dist   ) )
    assert ( associated(y%p_env%p_context) )
    assert ( y%p_env%p_context%created .eqv. .true.)
    if(y%p_env%p_context%iam<0) return

    ! Check matching partition/handler
    assert ( x%state /= undefined .and. y%state /= undefined  )
    assert ( x%state == y%state )
    call fem_vector_pxpy ( x%f_vector, y%f_vector )

  end subroutine par_vector_pxpy

  !=============================================================================
  subroutine par_vector_pxmy(x,y)
    implicit none
    type(par_vector_t), intent(in)     :: x
    type(par_vector_t), intent(inout)  :: y
    integer(ip)                   :: ierrc

    ! Pointer to part/context object is required
    assert ( associated(x%dof_dist   ) )
    assert ( associated(x%p_env%p_context) )
    assert ( associated(y%dof_dist   ) )
    assert ( associated(y%p_env%p_context) )
    assert ( y%p_env%p_context%created .eqv. .true.)
    if(y%p_env%p_context%iam<0) return

    ! Check matching partition/handler
    assert ( x%state /= undefined .and. y%state /= undefined  )
    assert ( x%state == y%state )
    ! write (*,*) 'XA'
    ! call fem_vector_print ( 6, x%f_vector )
    ! write (*,*) 'YA'
    ! call fem_vector_print ( 6, y%f_vector )
    call fem_vector_pxmy ( x%f_vector, y%f_vector )
    ! write (*,*) 'YD'
    ! call fem_vector_print ( 6, y%f_vector )
  end subroutine par_vector_pxmy

  subroutine par_vector_print ( luout, p_vec )
    implicit none
    integer(ip),      intent(in) :: luout
    type(par_vector_t), intent(in) :: p_vec

    ! Pointer to part/context object is required
    assert ( associated(p_vec%dof_dist   ) )
    assert ( associated(p_vec%p_env%p_context) )
    assert ( p_vec%p_env%p_context%created .eqv. .true.)
    if(p_vec%p_env%p_context%iam<0) return

    write(luout,'(a)') '*** begin par_vector data structure ***'
    call  dof_distribution_print (luout, p_vec%dof_dist)
    call  fem_vector_print    (luout, p_vec%f_vector)
    write(luout,'(a)') '*** end par_vector data structure ***'

  end subroutine par_vector_print

  subroutine par_vector_print_matrix_market ( dir_path, prefix, p_vec )
    implicit none
    ! Parameters
    character *(*)  , intent(in)  :: dir_path
    character *(*)  , intent(in)  :: prefix
    type(par_vector_t), intent(in)  :: p_vec

    ! Locals
    integer         :: iam, num_procs, lunou
    integer(ip)     :: j, ndigs_iam, ndigs_num_procs, id_map
    character(256)  :: name 
    character(256)  :: zeros
    character(256)  :: part_id

    assert ( associated(p_vec%dof_dist   ) )
    assert ( associated(p_vec%p_env%p_context) )
    assert ( p_vec%p_env%p_context%created .eqv. .true.)
    if(p_vec%p_env%p_context%iam<0) return

    name = trim(prefix) // '.par_vector' // '.mtx'

    ! Get context info
    call par_context_info ( p_vec%p_env%p_context, iam, num_procs )

    ! Form the file_path of the partition object to be read
    iam = iam + 1 ! Partition identifers start from 1 !!

    ndigs_num_procs = count_digits_par_vector (num_procs)
    zeros = ' '   
    ndigs_iam = count_digits_par_vector ( iam )

    ! write(*,*) ndgs_num_procs, ndigs_iam DBG

    do j=1,  ndigs_num_procs - ndigs_iam
       zeros (j:j) = '0'
    end do
    part_id = ch(iam)

    ! Read fem_partition data from path_file file
    lunou =  io_open (trim(dir_path) // '/' // trim(name) // '.' // trim(zeros) // trim(part_id), 'write')       
    call fem_vector_print_matrix_market ( lunou, p_vec%f_vector )
    call io_close (lunou)

  end subroutine par_vector_print_matrix_market

  function count_digits_par_vector ( i )
    implicit none
    ! Parameters
    integer(ip), intent(in) :: i 
    integer(ip)             :: count_digits_par_vector
    ! Locals   
    integer(ip)             :: x 
    x = i 
    if (x < 0) x = -x;
    count_digits_par_vector = 1;
    x = x/10;
    do while( x > 0)
       count_digits_par_vector = count_digits_par_vector + 1
       x = x/10;
    end do
  end function count_digits_par_vector


  ! Auxiliary (module private) routine
  ! for computing dot products and euclidean
  ! norms for element-based data distributions 
  subroutine weighted_dot (x,y,t)
    implicit none
    ! Parameters  
    type(par_vector_t), intent(in)  :: x
    type(par_vector_t), intent(in)  :: y
    real(rp)    , intent(out)     :: t

    ! Local variables
    integer(ip)                   :: iobj, i, id, i1, i2 
    real(rp)                      :: weigt

    ! Pointer to part/context object is required
    assert ( associated(x%dof_dist   ) )
    assert ( associated(x%p_env%p_context) )
    assert ( associated(y%dof_dist   ) )
    assert ( associated(y%p_env%p_context) )

    t = 0.0_rp 
    do iobj=2, x%dof_dist%nobjs
       weigt=1.0_rp/real(x%dof_dist%lobjs(4,iobj))
       i1=x%dof_dist%lobjs(2,iobj) - x%dof_dist%ni
       i2=x%dof_dist%lobjs(3,iobj) - x%dof_dist%ni
       call flat_weighted_dot ( i1, i2, & 
            & x%f_vector%b, y%f_vector%b, weigt, t) 
    end do
  end subroutine weighted_dot

  ! Auxiliary (module private) routine
  ! for computing dot products and euclidean
  ! norms for element-based data distributions 
  subroutine flat_weighted_dot ( i1, i2, x, y, weigt, t)
    implicit none
    ! Parameters  
    integer(ip), intent(in)    :: i1, i2
    real(rp),    intent(in)    :: x(*), y(*), weigt
    real(rp),    intent(inout) :: t

    ! Locals   
    integer(ip)                :: i, id

    do i=i1,i2
          t = t + weigt * x(i) * y(i)
    end do
  end subroutine flat_weighted_dot

  ! alpha <- op1^T * op2
 function par_vector_dot_tbp(op1,op2) result(alpha)
   implicit none
   class(par_vector_t), intent(in)    :: op1
   class(base_operand_t), intent(in)  :: op2
   real(rp) :: alpha

   ! Locals 
   integer(ip)                   :: ni 
   type(par_vector_t)              :: x_I, x_G, y_I, y_G
   real(rp)                      :: s

   ! Alpha should be defined in all tasks (not only in coarse-grid ones)
   alpha = 0.0_rp

   call op1%GuardTemp()
   call op2%GuardTemp()
   select type(op2)
   class is (par_vector_t)
      ! Pointer to part/context object is required
      assert ( associated(op1%dof_dist   ) )
      assert ( associated(op1%p_env%p_context) )
      assert ( associated(op2%dof_dist   ) )
      assert ( associated(op2%p_env%p_context) )
      assert ( op1%p_env%p_context%created .eqv. .true.)
      if(op1%p_env%p_context%iam<0) return
      
      assert ( op1%state /= undefined .and. op2%state /= undefined  )
      
      ni = op1%f_vector%neq - op1%dof_dist%nb
      if ( ni > 0 ) then
         call par_vector_create_view ( op1, 1, ni, x_I )
         call par_vector_create_view ( op2, 1, ni, y_I )
         alpha = x_I%f_vector%dot ( y_I%f_vector )
      else
         alpha = 0.0_rp
      end if
      
      call par_vector_create_view ( op1, ni+1, op1%f_vector%neq, x_G )
      call par_vector_create_view ( op2, ni+1, op2%f_vector%neq, y_G )
      call dot_interface          ( x_G, y_G, s )
      
      alpha = alpha + s
      
      ! Reduce-sum local dot products on all processes
      call psb_sum ( op2%p_env%p_context%icontxt, alpha )
   class default
      write(0,'(a)') 'par_vector_t%dot: unsupported op2 class'
      check(1==0)
   end select
   call op1%CleanTemp()
   call op2%CleanTemp()
 end function par_vector_dot_tbp

 ! op1 <- op2 
 subroutine par_vector_copy_tbp(op1,op2)
   implicit none
   class(par_vector_t), intent(inout) :: op1
   class(base_operand_t), intent(in)  :: op2
   
   call op2%GuardTemp()
   select type(op2)
   class is (par_vector_t)
      ! Pointer to part/context object is required
      assert ( associated(op2%dof_dist   ) )
      assert ( associated(op2%p_env%p_context) )
      assert ( associated(op1%dof_dist   ) )
      assert ( associated(op1%p_env%p_context) )
      assert ( op2%p_env%p_context%created .eqv. .true.)
      if(op2%p_env%p_context%iam<0) return
      
      assert ( op2%state /= undefined  )
      ! Perform local copy
      call op1%f_vector%copy ( op2%f_vector )
      op1%state = op2%state
   class default
      write(0,'(a)') 'par_vector_t%copy: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine par_vector_copy_tbp

 ! op1 <- alpha * op2
 subroutine par_vector_scal_tbp(op1,alpha,op2)
   implicit none
   class(par_vector_t), intent(inout) :: op1
   real(rp), intent(in) :: alpha
   class(base_operand_t), intent(in) :: op2

   call op2%GuardTemp()
   select type(op2)
   class is (par_vector_t)
      ! Pointer to part/context object is required
      assert ( associated(op2%dof_dist   ) )
      assert ( associated(op2%p_env%p_context) )
      assert ( associated(op1%dof_dist   ) )
      assert ( associated(op1%p_env%p_context) )
      assert ( op1%p_env%p_context%created .eqv. .true.)
      if(op1%p_env%p_context%iam<0) return
      
      ! Check matching partition/handler
      assert ( op2%state /= undefined )
      
      ! Scal local copy
      call op1%f_vector%scal ( alpha, op2%f_vector )
      op1%state = op2%state
   class default
      write(0,'(a)') 'par_vector_t%scal: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine par_vector_scal_tbp
 ! op <- alpha
 subroutine par_vector_init_tbp(op,alpha)
   implicit none
   class(par_vector_t), intent(inout) :: op
   real(rp), intent(in) :: alpha
   
   ! Pointer to part/context object is required
   assert ( associated(op%dof_dist   ) )
   assert ( associated(op%p_env%p_context) )
   assert ( op%p_env%p_context%created .eqv. .true.)
   if(op%p_env%p_context%iam<0) return
  
   ! Init local copy
   call op%f_vector%init(alpha)
   ! op%state = full_summed
 end subroutine par_vector_init_tbp

 ! op1 <- alpha*op2 + beta*op1
 subroutine par_vector_axpby_tbp(op1, alpha, op2, beta)
   implicit none
   class(par_vector_t), intent(inout) :: op1
   real(rp), intent(in) :: alpha
   class(base_operand_t), intent(in) :: op2
   real(rp), intent(in) :: beta

   call op2%GuardTemp()
   select type(op2)
   class is (par_vector_t)
      ! Pointer to part/context object is required
      assert ( associated(op2%dof_dist   ) )
      assert ( associated(op2%p_env%p_context) )
      assert ( associated(op1%dof_dist   ) )
      assert ( associated(op1%p_env%p_context) )
      assert ( op1%p_env%p_context%created )
      if(op1%p_env%p_context%iam<0) return
      
      ! Check matching partition/handler
      assert ( op1%state /= undefined .and. op2%state /= undefined  )
      assert ( op1%state == op2%state )
      
      call op1%f_vector%axpby( alpha, op2%f_vector, beta )
   class default
      write(0,'(a)') 'par_vector_t%axpby: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine par_vector_axpby_tbp

 ! alpha <- nrm2(op)
 function par_vector_nrm2_tbp(op) result(alpha)
   implicit none
   class(par_vector_t), intent(in)  :: op
   real(rp) :: alpha
   
   ! Alpha should be defined in all tasks (not only in coarse-grid ones)
   alpha = 0.0_rp

   call op%GuardTemp()
   select type(op)
   class is (par_vector_t)
      ! p_env%p_context is required within this subroutine 
      assert ( associated(op%dof_dist) )
      assert ( associated(op%p_env%p_context) )
      assert ( op%p_env%p_context%created )
      if(op%p_env%p_context%iam<0) return
      assert ( op%state /= undefined )
      alpha = op%dot(op)
      alpha = sqrt(alpha)
   class default
      write(0,'(a)') 'par_vector_t%nrm2: unsupported op2 class'
      check(1==0)
   end select   
   call op%CleanTemp()
 end function par_vector_nrm2_tbp

 ! op1 <- clone(op2) 
 subroutine par_vector_clone_tbp(op1,op2)
   implicit none
   class(par_vector_t)          , intent(inout) :: op1
   class(base_operand_t), target, intent(in)    :: op2

   call op2%GuardTemp()
   select type(op2)
   class is (par_vector_t)
      ! p_env%p_context is required within this subroutine 
      assert ( associated(op2%dof_dist) )
      assert ( associated(op2%p_env%p_context) )
      assert ( op2%p_env%p_context%created )
      
      op1%dof_dist => op2%dof_dist
      op1%p_env    => op2%p_env
      op1%state    =  op2%state

      if(op2%p_env%p_context%iam<0) return

      call op1%f_vector%clone ( op2%f_vector )

!!$      if (op1%mode == allocated) call memfreep(op1%b,__FILE__,__LINE__)
!!$      op1%neq     =  op2%neq       ! Number of equations
!!$         call memallocp(op1%neq,op1%b,__FILE__,__LINE__)
!!$      ! AFM: I think that clone should NOT init the memory just allocated.
!!$      ! The code that surrounds clone (e.g., Krylov solvers) should not
!!$      ! rely on par_vector_clone_tbp initializing the memory. I will comment
!!$      ! it out, and expect that the codes continue working.
!!$      ! op1%b = 0.0_rp 
!!$      op1%mode = allocated
   class default
      write(0,'(a)') 'par_vector_t%clone: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine par_vector_clone_tbp

 ! op <- comm(op)
 subroutine par_vector_comm_tbp(op)
   implicit none
   class(par_vector_t), intent(inout) :: op

   ! Local variables
   integer(ip)      :: ni 
   type(par_vector_t) :: op_G

   ! Pointer to part/context object is required
   assert ( associated(op%dof_dist) )
   assert ( associated(op%p_env%p_context) )
   assert ( op%p_env%p_context%created .eqv. .true.)
   if(op%p_env%p_context%iam<0) return
   
   ni = op%f_vector%neq - op%dof_dist%nb
   call par_vector_create_view ( op, ni+1, op%f_vector%neq, op_G )
   call comm_interface ( op_G )
   
   op%state = full_summed

 end subroutine par_vector_comm_tbp

 subroutine par_vector_free_tbp(this)
   implicit none
   class(par_vector_t), intent(inout) :: this

   ! The routine requires the partition/context info
   assert ( associated( this%dof_dist ) )
   assert ( associated( this%p_env%p_context ) )
   assert ( this%p_env%p_context%created .eqv. .true.)
   if(this%p_env%p_context%iam<0) return
   
   this%state = undefined
   
   ! Free local part
   call fem_vector_free ( this%f_vector )
   
   nullify ( this%dof_dist )
   nullify ( this%p_env )

 end subroutine par_vector_free_tbp


end module par_vector_names
