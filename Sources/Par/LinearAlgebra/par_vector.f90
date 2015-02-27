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
  use types
  use memor
  use fem_vector_names
  use fem_partition_names
  use fem_import_names
  use stdio
  use map_apply
#ifdef ENABLE_BLAS       
  use blas77_interfaces
#endif
#ifdef memcheck
  use iso_c_binding
#endif

  ! Module associated with the F90 interface to Trilinos.
  ! Remember: the F90 interface to Trilinos requires C
  ! interoperability (i.e., iso_c_binding module)
  !use for_trilinos_shadow_interfaces

  ! Parallel modules
  use par_partition_names
  use par_context_names
  use psb_penv_mod

# include "debug.i90"

  ! *** IMPORTANT NOTE: This cpp macro should go to a 
  ! common include file or it should be a program 
  ! subroutine otherwise
#define blk2scal(iv,idof,ndof) (((iv)-1)*(ndof)+(idof))

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

  integer(ip), parameter :: undefined     = -1 ! Undefined. State to be assigned by 
  ! the user or externally.
  integer(ip), parameter :: part_summed   = 0  ! partially summed element-based vector
  integer(ip), parameter :: full_summed   = 1  ! fully     summed element-based vector

  ! Distributed Vector
  type par_vector
     ! Local view of ONLY those components of f_vector
     ! corresponding to vertices owned by the processor  
     !type( epetra_vector ) :: epv_own

     ! Local view of BOTH: 1. those components of f_vector
     ! corresponding to vertices owned by the processor
     ! and 2. those corresponding to external nodes (owned
     ! by the neighbours of the processor)
     !type( epetra_vector ) :: epv_own_ext

     ! Data structure which stores the local part 
     ! of the vector mapped to the current processor.
     ! This is required for both eb and vb data 
     ! distributions
     type( fem_vector )    :: f_vector

     ! Partially or fully summed
     integer(ip)  :: state 

     ! Parallel partition control info.
     type ( par_partition ), pointer  :: p_part => NULL()
  end type par_vector

  interface par_vector_create_view
     module procedure par_vector_create_view_same_ndof, & 
                      par_vector_create_view_new_ndof
  end interface par_vector_create_view

  interface par_vector_l2g
     module procedure par_vector_l2g_fvec, par_vector_l2g_pvec
  end interface par_vector_l2g

  ! Types
  public :: par_vector

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
       &  par_vector_print,    par_vector_print_matrix_market,              &
       &  par_vector_l2g,      par_vector_g2l


!***********************************************************************
! Allocatable arrays of type(fem_vector)
!***********************************************************************
# define var_attr allocatable, target
# define point(a,b) call move_alloc(a,b)
# define generic_status_test             allocated
# define generic_memalloc_interface      memalloc
# define generic_memrealloc_interface    memrealloc
# define generic_memfree_interface       memfree
# define generic_memmovealloc_interface  memmovealloc

# define var_type type(par_vector)
# define var_size 52
# define bound_kind ip
# include "mem_header.i90"

  public :: memalloc,  memrealloc,  memfree, memmovealloc

contains

# include "mem_body.i90"

  !=============================================================================
  subroutine par_vector_alloc (storage, nd, p_part, p_vec)
    implicit none
    ! Parameters
    integer(ip)         ,         intent(in)  :: nd
    type(par_partition) , target, intent(in)  :: p_part
    integer(ip)         ,         intent(in)  :: storage
    type(par_vector)    ,         intent(out) :: p_vec

    ! Locals
    integer(ip)                               :: idof, id
    integer (c_int)     , allocatable         :: rc_map_values(:)

    ! p_part%p_context is required within this subroutine
    assert ( associated(p_part%p_context) )
    assert ( p_part%p_context%created .eqv. .true.)

    p_vec%p_part => p_part

    if(p_part%p_context%iam<0) return

    ! Allocate space for both own and external nodes.
    ! TO-DO: House for external nodes is not always required
    ! (e.g., for linear PDEs. right?). We should 
    ! determine when it is actually needed, and avoid allocating it 
    ! whenever it is not required 
    call fem_vector_alloc ( storage, nd, p_part%f_part%nmap%nl, p_vec%f_vector )

    ! Create epv_own and epv_own_ext as views of f_vector 
    ! if ( p_part%p_context%handler == trilinos ) then 
    !    ! epetra_vector (i.e., scalar vector)
    !    if ( storage == scal ) then
    !       assert ( nd <= max_ndofs )
    !       if ( p_vec%p_part%maps_state(nd) == map_non_created ) then
    !          call memalloc ( p_vec%p_part%f_part%nmap%nl*nd, rc_map_values,      __FILE__,__LINE__)

    !          do id=1, p_vec%p_part%f_part%nmap%nl 
    !             do idof=1, nd
    !                rc_map_values(blk2scal(id,idof,nd)) = blk2scal(p_vec%p_part%f_part%nmap%l2g(id),idof,nd)-1
    !             end do
    !          end do

    !          ! Create row-map
    !          call epetra_map_construct ( p_vec%p_part%row_map(nd), -1, & 
    !               &     (p_part%f_part%nmap%ni + p_part%f_part%nmap%nb)*nd, & 
    !               &     rc_map_values, 0, p_part%p_context%epcomm)
    !          ! Create col-map
    !          call epetra_map_construct ( p_vec%p_part%col_map(nd), -1, p_part%f_part%nmap%nl*nd, &
    !               &       rc_map_values, 0, p_part%p_context%epcomm) 
    !          ! Create importer                                          TARGET             SOURCE
    !          call epetra_import_construct ( p_vec%p_part%importer(nd), p_part%col_map(nd), p_part%row_map(nd) )

    !          call memfree ( rc_map_values,__FILE__,__LINE__)

    !          p_vec%p_part%maps_state(nd) = map_created
    !       end if
    !       call epetra_vector_construct ( p_vec%epv_own     , p_part%row_map(nd), p_vec%f_vector%b )
    !       call epetra_vector_construct ( p_vec%epv_own_ext , p_part%col_map(nd), p_vec%f_vector%b )
    !    else ! (block vector) 
    !       write (0,*) 'Error: trilinos shadow interfaces do not yet support block epetra_vector objects'
    !       stop
    !    end if
    ! else if ( p_part%p_context%handler == inhouse ) then
       p_vec%state = undefined
    ! end if

  end subroutine par_vector_alloc

  !=============================================================================
  subroutine par_vector_free(p_vec)
    implicit none
    type(par_vector), intent(inout) :: p_vec

    ! The routine requires the partition/context info
    assert ( associated( p_vec%p_part ) )
    assert ( associated( p_vec%p_part%p_context ) )
    assert ( p_vec%p_part%p_context%created .eqv. .true.)
    if(p_vec%p_part%p_context%iam<0) return

    ! if ( p_vec%p_part%p_context%handler == trilinos  ) then 
    !    ! epetra_vector (i.e., scalar vector)
    !    if ( p_vec%f_vector%storage == scal ) then 
    !       call epetra_vector_destruct ( p_vec%epv_own     )
    !       call epetra_vector_destruct ( p_vec%epv_own_ext )
    !    else ! (block vector)
    !       write (0,*) 'Error: trilinos shadow interfaces do not yet support block epetra_vector objects'
    !       stop
    !    end if
    ! elseif ( p_vec%p_part%p_context%handler == inhouse ) then
       p_vec%state = undefined
    ! end if

    ! Free local part
    call fem_vector_free ( p_vec%f_vector )

    nullify ( p_vec%p_part )
  end subroutine par_vector_free

  !=============================================================================
  subroutine par_vector_create_view_same_ndof (s_p_vec, start, end, t_p_vec)
    implicit none
    type(par_vector), intent(in), target :: s_p_vec
    integer(ip)     , intent(in)         :: start
    integer(ip)     , intent(in)         :: end
    type(par_vector), intent(out)        :: t_p_vec

    ! The routine requires the partition/context info
    assert ( associated( s_p_vec%p_part ) )
    assert ( associated( s_p_vec%p_part%p_context ) )
    assert ( s_p_vec%p_part%p_context%created .eqv. .true.)

    ! Associate parallel partition 
    t_p_vec%p_part => s_p_vec%p_part

    if ( s_p_vec%p_part%p_context%handler == inhouse ) then
       assert ( s_p_vec%state /= undefined )
       t_p_vec%state = s_p_vec%state
    end if
    
    if(s_p_vec%p_part%p_context%iam<0) return

    ! Call fem_vector_create_view
    call fem_vector_create_view ( s_p_vec%f_vector, start, end, t_p_vec%f_vector ) 

    ! Create epv_own and epv_own_ext as views of f_vector 
    ! if ( s_p_vec%p_part%p_context%handler == trilinos ) then 
    !    ! epetra_vector (i.e., scalar vector)
    !    if ( s_p_vec%f_vector%storage == scal ) then
    !       call epetra_vector_construct ( t_p_vec%epv_own,                             & 
    !                                      s_p_vec%p_part%row_map(s_p_vec%f_vector%nd), & 
    !                                      t_p_vec%f_vector%b )

    !       call epetra_vector_construct ( t_p_vec%epv_own_ext,                         & 
    !                                      s_p_vec%p_part%col_map(s_p_vec%f_vector%nd), & 
    !                                      t_p_vec%f_vector%b )
    !    else ! (block vector) 
    !       write (0,*) 'Error: trilinos shadow interfaces do not yet support block epetra_vector objects'
    !       stop
    !    end if
    ! end if
  end subroutine par_vector_create_view_same_ndof


  !=============================================================================
  subroutine par_vector_create_view_new_ndof (s_p_vec, ndstart, ndend, start, end, t_p_vec)
    implicit none
    ! Parameters
    type(par_vector), intent(in), target :: s_p_vec
    integer(ip)     , intent(in)         :: ndstart
    integer(ip)     , intent(in)         :: ndend
    integer(ip)     , intent(in)         :: start
    integer(ip)     , intent(in)         :: end
    type(par_vector), intent(out)        :: t_p_vec

    ! Locals
    integer(ip) :: nd
    integer(ip) :: idof, id
    integer (c_int) , allocatable  :: rc_map_values(:)

    ! The routine requires the partition/context info
    assert ( associated( s_p_vec%p_part ) )
    assert ( associated( s_p_vec%p_part%p_context ) )
    assert ( s_p_vec%p_part%p_context%created .eqv. .true.)

    ! Associate parallel partition 
    t_p_vec%p_part => s_p_vec%p_part

    if ( s_p_vec%p_part%p_context%handler == inhouse ) then
       assert ( s_p_vec%state /= undefined )
       t_p_vec%state = s_p_vec%state
    end if

    if(s_p_vec%p_part%p_context%iam<0) return

    ! Call fem_vector_create_view
    call fem_vector_create_view ( s_p_vec%f_vector, ndstart, ndend, start, end, t_p_vec%f_vector ) 

    ! Create epv_own and epv_own_ext as views of f_vector 
    ! if ( s_p_vec%p_part%p_context%handler == trilinos ) then 
    !    nd = ndend -ndstart + 1
    !    assert ( nd <= max_ndofs )
    !    if ( s_p_vec%p_part%maps_state(nd) == map_non_created ) then
    !       call memalloc ( s_p_vec%p_part%f_part%nmap%nl * nd, rc_map_values,   __FILE__,__LINE__)

    !       do id=1, s_p_vec%p_part%f_part%nmap%nl 
    !          do idof=1, nd
    !             rc_map_values(blk2scal(id,idof,nd)) = blk2scal(s_p_vec%p_part%f_part%nmap%l2g(id),idof,nd)-1
    !          end do
    !       end do

    !       ! Create row-map
    !       call epetra_map_construct ( t_p_vec%p_part%row_map(nd), -1, & 
    !            &                      (s_p_vec%p_part%f_part%nmap%ni + s_p_vec%p_part%f_part%nmap%nb)*nd, & 
    !            &                      rc_map_values, 0, s_p_vec%p_part%p_context%epcomm)

    !       ! Create col-map
    !       call epetra_map_construct ( t_p_vec%p_part%col_map(nd), -1, & 
    !            s_p_vec%p_part%f_part%nmap%nl*nd, &
    !            &                      rc_map_values, 0, s_p_vec%p_part%p_context%epcomm) 

    !       ! Create importer                                          
    !       call epetra_import_construct ( t_p_vec%p_part%importer(nd), & 
    !            s_p_vec%p_part%col_map(nd),   & ! TARGET
    !            s_p_vec%p_part%row_map(nd) )    ! SOURCE

    !       call memfree ( rc_map_values,__FILE__,__LINE__)

    !       t_p_vec%p_part%maps_state(nd) = map_created
    !    end if


    !    ! epetra_vector (i.e., scalar vector)
    !    if ( s_p_vec%f_vector%storage == scal ) then
    !       call epetra_vector_construct ( t_p_vec%epv_own,            & 
    !            s_p_vec%p_part%row_map(nd), & 
    !            t_p_vec%f_vector%b )

    !       call epetra_vector_construct ( t_p_vec%epv_own_ext,        & 
    !            s_p_vec%p_part%col_map(nd), & 
    !            t_p_vec%f_vector%b )
    !    else ! (block vector) 
    !       write (0,*) 'Error: trilinos shadow interfaces do not yet support block epetra_vector objects'
    !       stop
    !    end if
    ! end if
  end subroutine par_vector_create_view_new_ndof

  !=============================================================================
  subroutine par_vector_clone (s_p_vec, t_p_vec)
    implicit none
    type(par_vector), intent(in), target :: s_p_vec
    type(par_vector), intent(out)        :: t_p_vec 

    ! p_part%p_context is required within this subroutine 
    assert ( associated(s_p_vec%p_part) )
    assert ( associated(s_p_vec%p_part%p_context) )
    assert ( s_p_vec%p_part%p_context%created .eqv. .true.)

    t_p_vec%p_part => s_p_vec%p_part
    if ( t_p_vec%p_part%p_context%handler == inhouse ) then
       ! t_p_vec%state = undefined
       t_p_vec%state = s_p_vec%state
    end if

    if(s_p_vec%p_part%p_context%iam<0) return

    call fem_vector_clone ( s_p_vec%f_vector, t_p_vec%f_vector )

    ! Create epv_own and epv_own_ext as views of f_vector 
    ! if ( t_p_vec%p_part%p_context%handler == trilinos ) then 
    !    ! epetra_vector (i.e., scalar vector)
    !    if ( s_p_vec%f_vector%storage == scal ) then
    !       assert ( s_p_vec%p_part%maps_state(s_p_vec%f_vector%nd) == map_created  )
    !       call epetra_vector_construct ( t_p_vec%epv_own     , t_p_vec%p_part%row_map(s_p_vec%f_vector%nd), t_p_vec%f_vector%b )
    !       call epetra_vector_construct ( t_p_vec%epv_own_ext , t_p_vec%p_part%col_map(s_p_vec%f_vector%nd), t_p_vec%f_vector%b )
    !    else ! (block vector) 
    !       write (0,*) 'Error: trilinos shadow interfaces do not yet support block epetra_vector objects'
    !       stop
    !    end if
    ! end if

  end subroutine par_vector_clone

  subroutine par_vector_comm ( p_vec, alpha, mode )
    implicit none
    ! Parameters
    type(par_vector), intent(inout) :: p_vec
    real(rp)        , intent(in), optional  :: alpha
    integer(ip)     , intent(in), optional  :: mode

    ! Local variables
    integer(ip)      :: ni 
    type(par_vector) :: p_vec_G

    ! Pointer to part/context object is required
    assert ( associated(p_vec%p_part) )
    assert ( associated(p_vec%p_part%p_context) )
    assert ( p_vec%p_part%p_context%created .eqv. .true.)
    if(p_vec%p_part%p_context%iam<0) return

    ! assert ( p_vec%state /= undefined )
    if ( p_vec%p_part%f_part%ptype == element_based & 
         & .and. p_vec%p_part%p_context%handler == inhouse ) then
       ! Element-based partitioning/inhouse handler required
       assert ( p_vec%p_part%f_part%ptype == element_based )
       assert ( p_vec%p_part%p_context%handler == inhouse ) 

       ni = p_vec%f_vector%neq - p_vec%p_part%f_part%nmap%nb
       call par_vector_create_view ( p_vec, ni+1, p_vec%f_vector%neq, p_vec_G )
       call comm_interface ( p_vec_G, alpha, mode )

       p_vec%state = full_summed
    end if
  end subroutine par_vector_comm

  !=============================================================================
  ! VERY IMPORTANT: comm_interface is well-defined if and only if p_vec is a 
  !                 vector on the interface
  !=============================================================================
  subroutine comm_interface (p_vec, alpha, mode)
    use par_sparse_global_collectives
    implicit none

    ! Parameters
    type(par_vector), intent(inout)         :: p_vec
    real(rp)        , intent(in), optional  :: alpha
    integer(ip)     , intent(in), optional  :: mode

    ! Local variables
    real(rp)    ::  alpha_
    integer(ip) ::  mode_ 

    ! Pointer to part/context object is required
    assert ( associated(p_vec%p_part) )
    assert ( associated(p_vec%p_part%p_context) )
    assert ( p_vec%p_part%p_context%created .eqv. .true.)
    if(p_vec%p_part%p_context%iam<0) return

    ! assert ( p_vec%state /= undefined )


    ! Element-based partitioning/inhouse handler required
    assert ( p_vec%p_part%f_part%ptype == element_based )
    assert ( p_vec%p_part%p_context%handler == inhouse )  

    if (present( alpha )) then
       alpha_ = alpha
    else
       alpha_ = 1.0_rp
    end if

    if (present( mode )) then
       mode_ = mode 
    else
       mode_ = sp_all_to_all_default
    end if

    ! call fem_import_print (6, p_vec%p_part%f_import)
    ! write(*,*) 'SE 1', size(p_vec%f_vector%b)
    ! call fem_vector_print ( 6, p_vec%f_vector )

    ! First stage: owners receive/reduce, non-owners send
    call single_exchange ( p_vec%p_part%p_context%icontxt, &
         p_vec%p_part%f_import%num_rcv,    &
         p_vec%p_part%f_import%list_rcv,   &
         p_vec%p_part%f_import%rcv_ptrs,   &
         p_vec%p_part%f_import%unpack_idx, &
         p_vec%p_part%f_import%num_snd,    &
         p_vec%p_part%f_import%list_snd,   &
         p_vec%p_part%f_import%snd_ptrs,   &
         p_vec%p_part%f_import%pack_idx,   &
         alpha_,                         &
         1.0_rp,                         &
         p_vec%f_vector%nd,              &
         p_vec%f_vector%b,               &
         mode=mode_) 

    ! write(*,*) 'SE 2'

    ! Second stage: owners send, non-owners receive/insert
    call single_exchange ( p_vec%p_part%p_context%icontxt, &
         p_vec%p_part%f_import%num_snd,    &
         p_vec%p_part%f_import%list_snd,   &
         p_vec%p_part%f_import%snd_ptrs,   &
         p_vec%p_part%f_import%pack_idx,   &
         p_vec%p_part%f_import%num_rcv,    &
         p_vec%p_part%f_import%list_rcv,   &
         p_vec%p_part%f_import%rcv_ptrs,   &
         p_vec%p_part%f_import%unpack_idx, &
         1.0_rp,                         &
         0.0_rp,                         &
         p_vec%f_vector%nd,              &
         p_vec%f_vector%b,               &
         mode=mode_)

  end subroutine comm_interface

  subroutine par_vector_weight ( p_vec, weight )
    implicit none
    ! Parameters
    type(par_vector), intent(inout)                :: p_vec
    real(rp)        , intent(in), target, optional :: weight(*)

    ! Local variables
    integer(ip)      :: ni 
    type(par_vector) :: p_vec_G

    ! Pointer to part/context object is required
    assert ( associated(p_vec%p_part) )
    assert ( associated(p_vec%p_part%p_context) )
    assert ( p_vec%p_part%p_context%created .eqv. .true.)
    if(p_vec%p_part%p_context%iam<0) return

    ! assert ( p_vec%state /= undefined )


    ! Element-based partitioning/inhouse handler required
    assert ( p_vec%p_part%f_part%ptype == element_based )
    assert ( p_vec%p_part%p_context%handler == inhouse ) 

    ni = p_vec%f_vector%neq - p_vec%p_part%f_part%nmap%nb
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
    type(par_vector), intent(inout)                :: p_vec
    real(rp)        , intent(in), target, optional :: weight(*)

    ! Local variables
    integer(ip) :: iobj, i, i1, i2 
    real(rp)    :: weigt

    ! Pointer to part/context object is required
    assert ( associated(p_vec%p_part) )
    assert ( associated(p_vec%p_part%p_context) )
    assert ( p_vec%p_part%p_context%created .eqv. .true.)
    if(p_vec%p_part%p_context%iam<0) return

    ! assert ( p_vec%state /= undefined )


    ! Element-based partitioning/inhouse handler required
    assert ( p_vec%p_part%f_part%ptype == element_based )
    assert ( p_vec%p_part%p_context%handler == inhouse ) 
    
    if ( present(weight) ) then
       if ( p_vec%f_vector%storage == scal ) then
          i1 = blk2scal(                            1,                 1, p_vec%f_vector%nd)
          i2 = blk2scal(p_vec%p_part%f_part%nmap%nb  , p_vec%f_vector%nd, p_vec%f_vector%nd)
          do i=i1,i2
             p_vec%f_vector%b(:, i) =  p_vec%f_vector%b(:, i) * weight(i)
          end do
       else if ( p_vec%f_vector%storage == blk ) then
          ! Implementation pending for storage == blk !!!!
          assert (  p_vec%f_vector%storage == scal )
       end if
    else 

       if ( p_vec%f_vector%storage == scal ) then
          do iobj=2, p_vec%p_part%f_part%nobjs
             weigt=1.0_rp/real(p_vec%p_part%f_part%lobjs(4,iobj))
             i1 = blk2scal(p_vec%p_part%f_part%lobjs(2,iobj) - p_vec%p_part%f_part%nmap%ni,                 1, p_vec%f_vector%nd)
             i2 = blk2scal(p_vec%p_part%f_part%lobjs(3,iobj) - p_vec%p_part%f_part%nmap%ni, p_vec%f_vector%nd, p_vec%f_vector%nd)
             ! write (*,*) 'yyy', i1, i2, (i2-i1+1), weigt
#ifdef ENABLE_BLAS
             call dscal ( (i2-i1+1), weigt, p_vec%f_vector%b(:, i1:i2), 1 )
#else
             p_vec%f_vector%b(:, i1:i2) = weigt * p_vec%f_vector%b(:, i1:i2)
#endif 
          end do
       else if ( p_vec%f_vector%storage == blk ) then
          do iobj=2, p_vec%p_part%f_part%nobjs
             weigt=1.0_rp/real(p_vec%p_part%f_part%lobjs(4,iobj))
             i1 = p_vec%p_part%f_part%lobjs(2,iobj) - p_vec%p_part%f_part%nmap%ni
             i2 = p_vec%p_part%f_part%lobjs(3,iobj) - p_vec%p_part%f_part%nmap%ni
#ifdef ENABLE_BLAS
             call dscal ( (i2-i1+1)*p_vec%f_vector%nd, weigt, p_vec%f_vector%b(:, i1:i2), 1 )
#else
             p_vec%f_vector%b(:, i1:i2) = weigt * p_vec%f_vector%b(:, i1:i2)
#endif 
          end do
       end if
    end if
  end subroutine weight_interface

  subroutine par_vector_nrm2 (x, t)
    implicit none
    type(par_vector)   , intent(in), target  :: x
    real(rp)           , intent(out)         :: t
    integer(c_int)                           :: ierrc

    ! Pointer to part/context object is required
    assert ( associated(x%p_part) )
    assert ( associated(x%p_part%p_context) )
    assert ( x%p_part%p_context%created .eqv. .true.)
    if(x%p_part%p_context%iam<0) return

    ! if ( x%p_part%p_context%handler == inhouse ) then
       assert ( x%state /= undefined )
       ! write (*,*) 'PPP'
       ! call fem_vector_print ( 6, x%f_vector )
       call par_vector_dot (x,x,t)
       t = sqrt(t)
    ! else if ( x%p_part%p_context%handler == trilinos ) then
    !    if ( x%f_vector%storage == scal ) then 
    !       call epetra_vector_norm2 (x%epv_own, t, ierrc)
    !       assert ( ierrc == 0 )
    !    else 
    !       write (0,*) 'Error: trilinos shadow interfaces do not yet support block epetra_vector objects'
    !       stop
    !    end if
    ! end if
  end subroutine par_vector_nrm2

!!$  !=============================================================================
!!$  ! VERY IMPORTANT: nrm2_part_summed is well-defined if and only if x is a 
!!$  !                 vector on the interface
!!$  !=============================================================================
!!$  subroutine nrm2_part_summed (x, t)
!!$    implicit none
!!$    type(par_vector), intent(in), target :: x
!!$    real(rp)        , intent(out)        :: t
!!$    integer                              :: icontxt
!!$    type(par_vector)                     :: ws_vec
!!$            
!!$    ! Pointer to part/context object is required
!!$    assert ( associated(x%p_part   ) )
!!$    assert ( associated(x%p_part%p_context) )
!!$
!!$    ! Element-based partitioning/inhouse handler required
!!$    assert ( x%p_part%f_part%ptype == element_based )
!!$    assert ( x%p_part%p_context%handler == inhouse )
!!$    
!!$    ! This routine requires a partially summed interface
!!$    assert ( x%state == part_summed ) 
!!$
!!$    ! Allocate space for ws_vec
!!$    call par_vector_clone ( x, ws_vec )
!!$    
!!$    ! ws_vec <- x  
!!$    call par_vector_copy  ( x, ws_vec )
!!$
!!$    ! Transform ws_vec from partially summed 
!!$    ! to fully summed
!!$    call par_vector_comm ( ws_vec )
!!$    call par_vector_dot  ( x, ws_vec, t )
!!$    t = sqrt (t)
!!$
!!$    call par_vector_free ( ws_vec )
!!$  end subroutine nrm2_part_summed
!!$
!!$  !=============================================================================
!!$  ! VERY IMPORTANT: nrm2_full_summed is well-defined if and only if x is a 
!!$  !                 vector on the interface
!!$  !=============================================================================  
!!$  subroutine nrm2_full_summed (x, t)
!!$    implicit none
!!$    type(par_vector ), intent(in)    :: x
!!$    real(rp)         , intent(out)   :: t
!!$            
!!$    ! Pointer to part/context object is required
!!$    assert ( associated(x%p_part   ) )
!!$    assert ( associated(x%p_part%p_context) )
!!$    
!!$    ! This routine requires a partially summed interface 
!!$    assert ( x%state == full_summed ) 
!!$
!!$    ! Element-based partitioning/inhouse handler required
!!$    assert ( x%p_part%f_part%ptype == element_based )
!!$    assert ( x%p_part%p_context%handler == inhouse )
!!$
!!$    ! Perform local weighted dot products
!!$    call weighted_dot (x, x, t)
!!$    
!!$    ! Reduce-sum local dot products on all processes
!!$    call psb_sum ( x%p_part%p_context%icontxt, t )    
!!$    t = sqrt (t)
!!$  end subroutine nrm2_full_summed


  ! !=============================================================================
  ! subroutine par_vector_get_external_data (p_vec)
  !   implicit none


  !   ! Parameters
  !   type(par_vector), intent(inout)  :: p_vec

  !   ! Local variables
  !   integer(c_int)                   :: ierrc

  !   ! Pointer to part/context object is required
  !   assert ( associated(p_vec%p_part   ) )
  !   assert ( associated(p_vec%p_part%p_context) )

  !   ! Vertex-based partitioning/trilinos handler required for Epetra
  !   assert ( p_vec%p_part%f_part%ptype == vertex_based )
  !   assert ( p_vec%p_part%p_context%handler == trilinos )  

  !   if ( p_vec%f_vector%storage == scal ) then 
  !      ! Get external data !!!
  !      assert ( p_vec%p_part%maps_state(p_vec%f_vector%nd) == map_created  )
  !      call epetra_vector_import_with_importer_insert ( p_vec%epv_own_ext, p_vec%epv_own, p_vec%p_part%importer(p_vec%f_vector%nd), ierrc )
  !      assert ( ierrc == 0 )
  !   else 
  !      write (0,*) 'Error: trilinos shadow interfaces do not yet support block epetra_vector objects'
  !      stop
  !   end if

  ! end subroutine par_vector_get_external_data

  !=============================================================================
  ! VERY IMPORTANT: dot_interface is well-defined if and only if x/y are 
  !                 vectors on the interface
  !============================================================================= 
  subroutine dot_interface(x,y,t)
    implicit none
    ! Parameters  
    type(par_vector), intent(in)  :: x
    type(par_vector), intent(in)  :: y
    real(rp)    , intent(out)     :: t

    ! Locals 
    integer(c_int)                :: ierrc
    type(par_vector)              :: ws_vec

    ! Pointer to part/context object is required
    assert ( associated(x%p_part   ) )
    assert ( associated(x%p_part%p_context) )
    assert ( associated(y%p_part   ) )
    assert ( associated(y%p_part%p_context) )

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
    type(par_vector), intent(in)  :: x
    type(par_vector), intent(in)  :: y
    real(rp)        , intent(out) :: t

    ! Locals 
    integer(c_int)                :: ierrc
    integer(ip)                   :: ni 
    type(par_vector)              :: x_I, x_G, y_I, y_G
    real(rp)                      :: s

    ! Pointer to part/context object is required
    assert ( associated(x%p_part   ) )
    assert ( associated(x%p_part%p_context) )
    assert ( associated(y%p_part   ) )
    assert ( associated(y%p_part%p_context) )
    assert ( x%p_part%p_context%created .eqv. .true.)
    if(x%p_part%p_context%iam<0) return

    ! if ( y%p_part%p_context%handler == inhouse ) then
       assert ( x%state /= undefined .and. y%state /= undefined  )

       ni = x%f_vector%neq - x%p_part%f_part%nmap%nb
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
       call psb_sum ( y%p_part%p_context%icontxt, t )
    ! else if ( y%p_part%p_context%handler == trilinos ) then
    !    if ( y%f_vector%storage == scal ) then 
    !       call epetra_vector_dot (x%epv_own, y%epv_own, t, ierrc)
    !       assert ( ierrc == 0 )
    !    else 
    !       write (0,*) 'Error: trilinos shadow interfaces do not yet support block epetra_vector objects'
    !       stop
    !    end if
    ! end if

  end subroutine par_vector_dot

  !=============================================================================
  subroutine par_vector_copy(x,y)
    implicit none
    type(par_vector), intent(in)    :: x
    type(par_vector), intent(inout) :: y
    integer(c_int)                  :: ierrc

    ! Pointer to part/context object is required
    assert ( associated(x%p_part   ) )
    assert ( associated(x%p_part%p_context) )
    assert ( associated(y%p_part   ) )
    assert ( associated(y%p_part%p_context) )
    assert ( x%p_part%p_context%created .eqv. .true.)
    if(x%p_part%p_context%iam<0) return

    ! Check matching partition/handler
    assert ( x%p_part%f_part%ptype == y%p_part%f_part%ptype )
    assert ( x%p_part%p_context%handler == y%p_part%p_context%handler ) 

    ! if ( y%p_part%p_context%handler == inhouse ) then
       assert ( x%state /= undefined  )
       ! Perform local copy
       call fem_vector_copy ( x%f_vector, y%f_vector )
       y%state = x%state
    ! else if( y%p_part%p_context%handler == trilinos ) then
    !    if ( y%f_vector%storage == scal ) then 
    !       call epetra_vector_scale (y%epv_own, 1.0_c_double, x%epv_own, ierrc)
    !       assert ( ierrc == 0 )
    !    else 
    !       write (0,*) 'Error: trilinos shadow interfaces do not yet support block epetra_vector objects'
    !       stop
    !    end if
    ! end if

  end subroutine par_vector_copy

  !=============================================================================
  subroutine par_vector_zero(y)
    implicit none
    type(par_vector), intent(inout) :: y
    integer(c_int)                  :: ierrc

    ! Pointer to part/context object is required
    assert ( associated(y%p_part   ) )
    assert ( associated(y%p_part%p_context) )
    assert ( y%p_part%p_context%created .eqv. .true.)
    if(y%p_part%p_context%iam<0) return


    if ( y%p_part%p_context%handler == inhouse ) then
       ! write(*,*) 'XXX'
       assert ( y%state /= undefined  )
       ! Zero-out local copy
       call fem_vector_zero ( y%f_vector )
    else if( y%p_part%p_context%handler == trilinos ) then
       if ( y%f_vector%storage == scal ) then 
          ! call epetra_vector_scale (y%epv_own, 0.0_c_double, y%epv_own, ierrc)
          call fem_vector_zero ( y%f_vector )
          ! assert ( ierrc == 0 )
       else 
          write (0,*) 'Error: trilinos shadow interfaces do not yet support block epetra_vector objects'
          stop
       end if
    end if

  end subroutine par_vector_zero

  !=============================================================================
  subroutine par_vector_init (t,y)
    implicit none
    real(rp)        , intent(in)     :: t
    type(par_vector), intent(inout)  :: y
    integer(c_int)                   :: ierrc

    ! Pointer to part/context object is required
    assert ( associated(y%p_part   ) )
    assert ( associated(y%p_part%p_context) )
    assert ( y%p_part%p_context%created .eqv. .true.)
    if(y%p_part%p_context%iam<0) return

    if  ( t == 0.0_rp ) then
       call par_vector_zero(y)
    else 
       if (y%p_part%p_context%handler == inhouse) then
          ! Scale local copy
          call fem_vector_init (t, y%f_vector)
          y%state = full_summed
       else if( y%p_part%p_context%handler == trilinos ) then
          if ( y%f_vector%storage == scal ) then 
             call fem_vector_init (t, y%f_vector )
          else 
             write (0,*) 'Error: trilinos shadow interfaces do not yet support block epetra_vector objects'
             stop
          end if
       end if
    end if
  end subroutine par_vector_init


  !=============================================================================
  subroutine par_vector_scale(t,x,y)
    implicit none
    real(rp)    , intent(in)         :: t
    type(par_vector), intent(in)     :: x
    type(par_vector), intent(inout)  :: y
    integer(c_int)                   :: ierrc

    ! Pointer to part/context object is required
    assert ( associated(x%p_part   ) )
    assert ( associated(x%p_part%p_context) )
    assert ( associated(y%p_part   ) )
    assert ( associated(y%p_part%p_context) )
    assert ( y%p_part%p_context%created .eqv. .true.)
    if(y%p_part%p_context%iam<0) return

    ! Check matching partition/handler
    assert ( x%p_part%f_part%ptype == y%p_part%f_part%ptype )
    assert ( x%p_part%p_context%handler == y%p_part%p_context%handler ) 

    ! if (y%p_part%p_context%handler == inhouse) then
       assert ( x%state /= undefined )
       ! Scale local copy
       call fem_vector_scale (t, x%f_vector, y%f_vector )
       y%state = x%state
    ! else if( y%p_part%p_context%handler == trilinos ) then
    !    if ( y%f_vector%storage == scal ) then 
    !       call epetra_vector_scale (y%epv_own, t, x%epv_own, ierrc)
    !       assert ( ierrc == 0 )
    !    else 
    !       write (0,*) 'Error: trilinos shadow interfaces do not yet support block epetra_vector objects'
    !       stop
    !    end if
    ! end if

  end subroutine par_vector_scale

  !=============================================================================
  subroutine par_vector_mxpy(x,y)
    implicit none
    type(par_vector), intent(in)    :: x
    type(par_vector), intent(inout) :: y
    integer(c_int)                  :: ierrc

    ! Pointer to part/context object is required
    assert ( associated(x%p_part   ) )
    assert ( associated(x%p_part%p_context) )
    assert ( associated(y%p_part   ) )
    assert ( associated(y%p_part%p_context) )
    assert ( y%p_part%p_context%created .eqv. .true.)
    if(y%p_part%p_context%iam<0) return

    ! Check matching partition/handler
    assert ( x%p_part%f_part%ptype == y%p_part%f_part%ptype )
    assert ( x%p_part%p_context%handler == y%p_part%p_context%handler ) 

    ! if ( y%p_part%p_context%handler == inhouse ) then
       assert ( x%state /= undefined .and. y%state /= undefined  )
       assert ( x%state == y%state )
       call fem_vector_mxpy (x%f_vector, y%f_vector )
    ! else if( y%p_part%p_context%handler == trilinos ) then
    !    if ( y%f_vector%storage == scal ) then 
    !       call epetra_vector_update ( y%epv_own, -1.0_c_double, x%epv_own, 1.0_c_double, ierrc  ) 
    !       assert ( ierrc == 0 )
    !    else 
    !       write (0,*) 'Error: trilinos shadow interfaces do not yet support block epetra_vector objects'
    !       stop
    !    end if
    ! end if
  end subroutine par_vector_mxpy

  !=============================================================================
  subroutine par_vector_axpy(t,x,y)
    implicit none
    real(rp)   , intent(in)         :: t
    type(par_vector), intent(in)    :: x
    type(par_vector), intent(inout) :: y
    integer(c_int)                  :: ierrc

    ! Pointer to part/context object is required
    assert ( associated(x%p_part   ) )
    assert ( associated(x%p_part%p_context) )
    assert ( associated(y%p_part   ) )
    assert ( associated(y%p_part%p_context) )
    assert ( y%p_part%p_context%created .eqv. .true.)
    if(y%p_part%p_context%iam<0) return

    ! Check matching partition/handler
    assert ( x%p_part%f_part%ptype == y%p_part%f_part%ptype )
    assert ( x%p_part%p_context%handler == y%p_part%p_context%handler ) 

    ! if ( y%p_part%p_context%handler == inhouse ) then
       assert ( x%state /= undefined .and. y%state /= undefined  )
       assert ( x%state == y%state )
       call fem_vector_axpy  ( t, x%f_vector, y%f_vector )
       ! call fem_vector_print ( 6, x%f_vector )
    ! else if( y%p_part%p_context%handler == trilinos ) then
    !    if ( y%f_vector%storage == scal ) then 
    !       call epetra_vector_update ( y%epv_own, t, x%epv_own, 1.0_c_double, ierrc  )
    !       assert ( ierrc == 0 )
    !    else 
    !       write (0,*) 'Error: trilinos shadow interfaces do not yet support block epetra_vector objects'
    !       stop
    !    end if
    ! end if
  end subroutine par_vector_axpy

  !=============================================================================
  subroutine par_vector_aypx(t,x,y)
    implicit none
    real(rp)    , intent(in)        :: t
    type(par_vector), intent(in)    :: x
    type(par_vector), intent(inout) :: y
    integer(c_int)                  :: ierrc

    ! Pointer to part/context object is required
    assert ( associated(x%p_part   ) )
    assert ( associated(x%p_part%p_context) )
    assert ( associated(y%p_part   ) )
    assert ( associated(y%p_part%p_context) )
    assert ( y%p_part%p_context%created .eqv. .true.)
    if(y%p_part%p_context%iam<0) return

    ! Check matching partition/handler
    assert ( x%p_part%f_part%ptype == y%p_part%f_part%ptype )
    assert ( x%p_part%p_context%handler == y%p_part%p_context%handler ) 

    ! if ( y%p_part%p_context%handler == inhouse ) then
       assert ( x%state /= undefined .and. y%state /= undefined  )
       assert ( x%state == y%state )
       call fem_vector_aypx (t, x%f_vector, y%f_vector )
    ! else if(y%p_part%p_context%handler == trilinos) then
    !    if ( y%f_vector%storage == scal ) then 
    !       call epetra_vector_update ( y%epv_own, 1.0_c_double, x%epv_own, t, ierrc  )
    !       assert ( ierrc == 0 )
    !    else 
    !       write (0,*) 'Error: trilinos shadow interfaces do not yet support block epetra_vector objects'
    !       stop
    !    end if
    ! end if
  end subroutine par_vector_aypx

  !=============================================================================
  subroutine par_vector_pxpy(x,y)
    implicit none
    type(par_vector), intent(in)    :: x
    type(par_vector), intent(inout) :: y
    integer(c_int)                  :: ierrc

    ! Pointer to part/context object is required
    assert ( associated(x%p_part   ) )
    assert ( associated(x%p_part%p_context) )
    assert ( associated(y%p_part   ) )
    assert ( associated(y%p_part%p_context) )
    assert ( y%p_part%p_context%created .eqv. .true.)
    if(y%p_part%p_context%iam<0) return

    ! Check matching partition/handler
    assert ( x%p_part%f_part%ptype == y%p_part%f_part%ptype )
    assert ( x%p_part%p_context%handler == y%p_part%p_context%handler ) 

    ! if ( y%p_part%p_context%handler == inhouse ) then
       assert ( x%state /= undefined .and. y%state /= undefined  )
       assert ( x%state == y%state )
       call fem_vector_pxpy ( x%f_vector, y%f_vector )
    ! else if(y%p_part%p_context%handler == trilinos) then
    !    if ( y%f_vector%storage == scal ) then 
    !       call epetra_vector_update ( y%epv_own, 1.0_c_double, x%epv_own, 1.0_c_double, ierrc  ) 
    !       assert ( ierrc == 0 )
    !    else 
    !       write (0,*) 'Error: trilinos shadow interfaces do not yet support block epetra_vector objects'
    !       stop
    !    end if
    ! end if
  end subroutine par_vector_pxpy

  !=============================================================================
  subroutine par_vector_pxmy(x,y)
    implicit none
    type(par_vector), intent(in)     :: x
    type(par_vector), intent(inout)  :: y
    integer(c_int)                   :: ierrc

    ! Pointer to part/context object is required
    assert ( associated(x%p_part   ) )
    assert ( associated(x%p_part%p_context) )
    assert ( associated(y%p_part   ) )
    assert ( associated(y%p_part%p_context) )
    assert ( y%p_part%p_context%created .eqv. .true.)
    if(y%p_part%p_context%iam<0) return

    ! Check matching partition/handler
    assert ( x%p_part%f_part%ptype == y%p_part%f_part%ptype )
    assert ( x%p_part%p_context%handler == y%p_part%p_context%handler ) 

    ! if (y%p_part%p_context%handler == inhouse) then
       assert ( x%state /= undefined .and. y%state /= undefined  )
       assert ( x%state == y%state )
       ! write (*,*) 'XA'
       ! call fem_vector_print ( 6, x%f_vector )
       ! write (*,*) 'YA'
       ! call fem_vector_print ( 6, y%f_vector )
       call fem_vector_pxmy ( x%f_vector, y%f_vector )
       ! write (*,*) 'YD'
       ! call fem_vector_print ( 6, y%f_vector )
    ! else if(y%p_part%p_context%handler == trilinos) then
    !    if ( y%f_vector%storage == scal ) then 
    !       call epetra_vector_update ( y%epv_own, 1.0_c_double, x%epv_own, -1.0_c_double, ierrc  ) 
    !       assert ( ierrc == 0 )
    !    else 
    !       write (0,*) 'Error: trilinos shadow interfaces do not yet support block epetra_vector objects'
    !       stop
    !    end if
    ! end if
  end subroutine par_vector_pxmy

  subroutine par_vector_print ( luout, p_vec )
    implicit none
    integer(ip),      intent(in) :: luout
    type(par_vector), intent(in) :: p_vec

    ! Pointer to part/context object is required
    assert ( associated(p_vec%p_part   ) )
    assert ( associated(p_vec%p_part%p_context) )
    assert ( p_vec%p_part%p_context%created .eqv. .true.)
    if(p_vec%p_part%p_context%iam<0) return

    write(luout,'(a)') '*** begin par_vector data structure ***'
    if  (p_vec%p_part%p_context%handler == inhouse) then
       write (luout, '(a,i10)') 'State:', p_vec%state
    end if
    call  par_partition_print (luout, p_vec%p_part)
    call  fem_vector_print    (luout, p_vec%f_vector)
    ! if ( p_vec%p_part%p_context%handler == trilinos ) then 
    !    call epetra_vector_print ( p_vec%epv_own      )
    !    call epetra_vector_print ( p_vec%epv_own_ext  )
    ! end if
    write(luout,'(a)') '*** end par_vector data structure ***'

  end subroutine par_vector_print

  subroutine par_vector_print_matrix_market ( dir_path, prefix, p_vec )
    implicit none
    ! Parameters
    character *(*)  , intent(in)  :: dir_path
    character *(*)  , intent(in)  :: prefix
    type(par_vector), intent(in)  :: p_vec

    ! Locals
    integer         :: iam, num_procs, lunou
    integer(ip)     :: j, ndigs_iam, ndigs_num_procs, id_map
    character(256)  :: name 
    character(256)  :: zeros
    character(256)  :: part_id

    assert ( associated(p_vec%p_part   ) )
    assert ( associated(p_vec%p_part%p_context) )
    assert ( p_vec%p_part%p_context%created .eqv. .true.)
    if(p_vec%p_part%p_context%iam<0) return

    name = trim(prefix) // '.par_vector' // '.mtx'

    ! Get context info
    call par_context_info ( p_vec%p_part%p_context, iam, num_procs )

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
    type(par_vector), intent(in)  :: x
    type(par_vector), intent(in)  :: y
    real(rp)    , intent(out)     :: t

    ! Local variables
    integer(ip)                   :: iobj, i, id, i1, i2 
    real(rp)                      :: weigt

    ! Pointer to part/context object is required
    assert ( associated(x%p_part   ) )
    assert ( associated(x%p_part%p_context) )
    assert ( associated(y%p_part   ) )
    assert ( associated(y%p_part%p_context) )

    ! Element-based partitioning/inhouse handler required
    assert ( x%p_part%f_part%ptype == element_based )
    assert ( x%p_part%p_context%handler == inhouse ) 

    t = 0.0_rp 
    do iobj=2, x%p_part%f_part%nobjs
       weigt=1.0_rp/real(x%p_part%f_part%lobjs(4,iobj))
       i1=x%p_part%f_part%lobjs(2,iobj) - x%p_part%f_part%nmap%ni
       i2=x%p_part%f_part%lobjs(3,iobj) - x%p_part%f_part%nmap%ni
       call flat_weighted_dot ( x%f_vector%nd, i1, i2, & 
            & x%f_vector%b, y%f_vector%b, weigt, t) 
    end do
  end subroutine weighted_dot

  ! Auxiliary (module private) routine
  ! for computing dot products and euclidean
  ! norms for element-based data distributions 
  subroutine flat_weighted_dot (nd, i1, i2, x, y, weigt, t)
    implicit none
    ! Parameters  
    integer(ip), intent(in)    :: nd, i1, i2
    real(rp),    intent(in)    :: x(nd,*), y(nd,*), weigt
    real(rp),    intent(inout) :: t

    ! Locals   
    integer(ip)                :: i, id

    do i=i1,i2
       do id=1,nd
          t = t + weigt * x(id,i) * y(id,i)
       end do
    end do
  end subroutine flat_weighted_dot

  !================================================================================================
  subroutine par_vector_l2g_pvec(p_vec,f_vec,root,ndime)
    implicit none
    type(par_vector), intent(in)  :: p_vec
    type(fem_vector), intent(out) :: f_vec
    integer(ip),      intent(in), optional :: root
    integer(ip),      intent(in), optional :: ndime

    integer(ip) :: me,np,neq,npoin,root_pid,mpi_comm,iret
    integer(ip) :: i,idime,k

    integer(ip), allocatable :: rcvcnt_l2g(:),displ_l2g(:),l2g_gather(:)
    integer(ip), allocatable :: rcvcnt_vec(:),displ_vec(:)
    real(rp),    allocatable :: vec_gather(:)

    assert ( associated(p_vec%p_part   ) )
    assert ( associated(p_vec%p_part%p_context) )
    assert ( p_vec%p_part%p_context%created .eqv. .true.)
    if(p_vec%p_part%p_context%iam<0) return

    ! Set root
    if(.not.present(root)) then
       root_pid = 0
    else
       root_pid = root
    end if

    ! Get par_context info
    call par_context_info(p_vec%p_part%p_context,me,np)
    call psb_get_mpicomm(p_vec%p_part%p_context%icontxt,mpi_comm)

    ! Allocate recvcounts
    if(me==root_pid) then
       call memalloc(np,rcvcnt_l2g,__FILE__,__LINE__)
       call memalloc(np,rcvcnt_vec,__FILE__,__LINE__)
    else
       call memalloc(0,rcvcnt_l2g,__FILE__,__LINE__)
       call memalloc(0,rcvcnt_vec,__FILE__,__LINE__)
    end if

    ! Gather recvcounts
    call mpi_gather(p_vec%p_part%f_part%nmap%nl,1,psb_mpi_integer,rcvcnt_l2g,1,psb_mpi_integer, &
         &          root_pid,mpi_comm,iret)
    call mpi_gather(p_vec%f_vector%nd*p_vec%f_vector%neq,1,psb_mpi_integer,rcvcnt_vec,1,        &
         &          psb_mpi_integer,root_pid,mpi_comm,iret)

    ! Set displs and allocate gathered vectors l2g and vec
    if(me==root_pid) then
       call memalloc(np+1,displ_l2g,__FILE__,__LINE__)
       call memalloc(np+1,displ_vec,__FILE__,__LINE__)
       displ_l2g(1) = 0
       displ_vec(1) = 0
       do i=1,np
          displ_l2g(i+1) = displ_l2g(i) + rcvcnt_l2g(i)
          displ_vec(i+1) = displ_vec(i) + rcvcnt_vec(i)
       end do
       call memalloc(displ_l2g(np+1),l2g_gather,__FILE__,__LINE__)
       call memalloc(displ_vec(np+1),vec_gather,__FILE__,__LINE__)
    else
       call memalloc(0,displ_l2g,__FILE__,__LINE__)
       call memalloc(0,displ_vec,__FILE__,__LINE__)
       call memalloc(0,l2g_gather,__FILE__,__LINE__)
       call memalloc(0,vec_gather,__FILE__,__LINE__)
    end if

    ! Gather local to global
    call mpi_gatherv(p_vec%p_part%f_part%nmap%l2g,p_vec%p_part%f_part%nmap%nl,psb_mpi_integer, &
         l2g_gather,rcvcnt_l2g,displ_l2g,psb_mpi_integer,root_pid,mpi_comm,iret)

    ! Gather vector data
    call mpi_gatherv(p_vec%f_vector%b,p_vec%f_vector%nd*p_vec%f_vector%neq,psb_mpi_real, &
         vec_gather,rcvcnt_vec,displ_vec,psb_mpi_real,root_pid,mpi_comm,iret)

    ! Create fem vector
    if(me==root_pid) then
       ! Set number of global equations
       neq = displ_vec(np)/p_vec%f_vector%nd
       call fem_vector_alloc(p_vec%f_vector%storage,p_vec%f_vector%nd,neq,f_vec)
       if(p_vec%f_vector%storage==scal) then
          npoin = neq/ndime
          k=1
          do i=1,npoin
             do idime=1,ndime
                f_vec%b(1,(l2g_gather(i)-1)*ndime+idime) = vec_gather(k)
                k=k+1
             end do
          end do
       else if(p_vec%f_vector%storage==blk) then
          k=1
          do i=1,neq
             do idime=1,p_vec%f_vector%nd
                f_vec%b(idime,l2g_gather(i)) = vec_gather(k)
                k=k+1
             end do
          end do
       end if
    end if

    ! Deallocate objects
    call memfree(displ_l2g,__FILE__,__LINE__)
    call memfree(displ_vec,__FILE__,__LINE__)
    call memfree(rcvcnt_l2g,__FILE__,__LINE__)
    call memfree(rcvcnt_vec,__FILE__,__LINE__)
    call memfree(l2g_gather,__FILE__,__LINE__)
    call memfree(vec_gather,__FILE__,__LINE__)

  end subroutine par_vector_l2g_pvec

  !================================================================================================
  subroutine par_vector_l2g_fvec(ndime,npoin,veloc,p_part,f_vec,root,l2g,nl)
    implicit none
    integer(ip),          intent(in)  :: ndime,npoin
    real(rp),             intent(in)  :: veloc(ndime,npoin)
    type(par_partition),  intent(in)  :: p_part
    type(fem_vector),     intent(out) :: f_vec
    integer(ip), intent(in), optional :: root
    integer(ip), intent(in), optional :: l2g(:)
    integer(ip), intent(in), optional :: nl

    integer(ip) :: me,np,neq,root_pid,mpi_comm,iret,nl_
    integer(ip) :: i,idime,k

    integer(ip), allocatable :: rcvcnt_l2g(:),displ_l2g(:),l2g_gather(:)
    integer(ip), allocatable :: rcvcnt_vec(:),displ_vec(:)
    integer(ip), allocatable :: l2g_(:)
    real(rp),    allocatable :: vec_gather(:)

    assert ( associated(p_part%p_context) )
    assert ( p_part%p_context%created .eqv. .true.)
    if(p_part%p_context%iam<0) return

    ! Set root
    if(.not.present(root)) then
       root_pid = 0
    else
       root_pid = root
    end if

    ! Set l2g
    if(present(l2g)) then
       check(present(nl))
       call memalloc(nl,l2g_,__FILE__,__LINE__)
       l2g_=l2g
       nl_ =nl
    else
       call memalloc(p_part%f_part%nmap%nl,l2g_,__FILE__,__LINE__)
       l2g_=p_part%f_part%nmap%l2g
       nl_ =p_part%f_part%nmap%nl
    end if

    ! Get par_context info
    call par_context_info(p_part%p_context,me,np)
    call psb_get_mpicomm(p_part%p_context%icontxt,mpi_comm)

    ! Allocate recvcounts
    if(me==root_pid) then
       call memalloc(np,rcvcnt_l2g,__FILE__,__LINE__)
       call memalloc(np,rcvcnt_vec,__FILE__,__LINE__)
    else
       call memalloc(0,rcvcnt_l2g,__FILE__,__LINE__)
       call memalloc(0,rcvcnt_vec,__FILE__,__LINE__)
    end if

    ! Gather recvcounts
    call mpi_gather(nl_,1,psb_mpi_integer,rcvcnt_l2g,1,psb_mpi_integer,root_pid,mpi_comm,iret)
    call mpi_gather(ndime*npoin,1,psb_mpi_integer,rcvcnt_vec,1,psb_mpi_integer, &
         &          root_pid,mpi_comm,iret)

    ! Set displs and allocate gathered vectors l2g and vec
    if(me==root_pid) then
       call memalloc(np+1,displ_l2g,__FILE__,__LINE__)
       call memalloc(np+1,displ_vec,__FILE__,__LINE__)
       displ_l2g(1) = 0
       displ_vec(1) = 0
       do i=1,np
          displ_l2g(i+1) = displ_l2g(i) + rcvcnt_l2g(i)
          displ_vec(i+1) = displ_vec(i) + rcvcnt_vec(i)
       end do
       call memalloc(displ_l2g(np+1),l2g_gather,__FILE__,__LINE__)
       call memalloc(displ_vec(np+1),vec_gather,__FILE__,__LINE__)
    else
       call memalloc(0,displ_l2g,__FILE__,__LINE__)
       call memalloc(0,displ_vec,__FILE__,__LINE__)
       call memalloc(0,l2g_gather,__FILE__,__LINE__)
       call memalloc(0,vec_gather,__FILE__,__LINE__)
    end if

    ! Gather local to global
    call mpi_gatherv(l2g_,nl_,psb_mpi_integer, &
         &           l2g_gather,rcvcnt_l2g,displ_l2g,psb_mpi_integer,root_pid,mpi_comm,iret)

    ! Gather vector data
    call mpi_gatherv(veloc,ndime*npoin,psb_mpi_real,vec_gather, &
         &           rcvcnt_vec,displ_vec,psb_mpi_real,root_pid,mpi_comm,iret)

    ! Create fem vector
    if(me==root_pid) then
       ! Set number of global equations
       neq = displ_vec(np+1)/ndime
       call fem_vector_alloc(blk,ndime,neq,f_vec)
       k=1
       do i=1,neq
          do idime=1,ndime
             f_vec%b(idime,l2g_gather(i)) = vec_gather(k)
             k=k+1
          end do
       end do
    end if

    ! Deallocate objects
    call memfree(displ_l2g,__FILE__,__LINE__)
    call memfree(displ_vec,__FILE__,__LINE__)
    call memfree(rcvcnt_l2g,__FILE__,__LINE__)
    call memfree(rcvcnt_vec,__FILE__,__LINE__)
    call memfree(l2g_gather,__FILE__,__LINE__)
    call memfree(vec_gather,__FILE__,__LINE__)
    call memfree(l2g_,__FILE__,__LINE__)

  end subroutine par_vector_l2g_fvec

  !================================================================================================
  subroutine par_vector_g2l(root_pid,ndime,npoin,gvec,p_part,lvec)
    implicit none
    integer(ip),         intent(in)  :: root_pid,ndime,npoin
    type(par_partition), intent(in)  :: p_part
    real(rp),            intent(in)  :: gvec(ndime,npoin)
    real(rp),            intent(out) :: lvec(ndime,p_part%f_part%nmap%nl)
  
    integer(ip) :: me,np,mpi_comm,iret,i
    integer(ip) :: idime

    ! Get par_context info
    call par_context_info(p_part%p_context,me,np)
    call psb_get_mpicomm(p_part%p_context%icontxt,mpi_comm)

    ! Broadcast vector
    call mpi_bcast(gvec,ndime*npoin,psb_mpi_real,root_pid,mpi_comm,iret)
    !call mpi_bcast(gvec,ndime*p_part%f_part%nmap%ng,psb_mpi_real,root_pid,mpi_comm,iret)

    ! Global to local mapping
    do i=1,p_part%f_part%nmap%nl
       lvec(:,i)=gvec(:,p_part%f_part%nmap%l2g(i))
    end do
    !call map_apply_g2l(p_part%f_part%nmap,ndime,gvec,lvec)

  end subroutine par_vector_g2l

end module par_vector_names
