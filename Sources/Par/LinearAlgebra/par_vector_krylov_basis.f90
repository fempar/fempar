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
module par_vector_krylov_basis_names
  ! Serial modules
  use types
  use memor
  use fem_vector_krylov_basis_names
  
  ! Module associated with the F90 interface to Trilinos.
  ! Remember: the F90 interface to Trilinos requires C
  ! interoperability (i.e., iso_c_binding module)
  !use for_trilinos_shadow_interfaces

  ! Parallel modules
  use par_partition_names
  use par_context_names
  use psb_penv_mod
  use par_vector_names

# include "debug.i90"
  
  implicit none
  private


  ! Distributed Krylov Basis (compatible with par_vector)
  type par_vector_krylov_basis
     ! Local view of ONLY those components of b
     ! corresponding to vertices owned by the processor  
     !type( epetra_multivector ) :: epmv

     type(fem_vector_krylov_basis) :: f_basis
     
     ! Partially or fully summed
     integer(ip)  :: state = undefined

     ! Parallel partition control info.
     type ( par_partition ), pointer  :: p_part => NULL()
  end type par_vector_krylov_basis

  ! Types
  public :: par_vector_krylov_basis


  ! Functions
  public :: par_vector_krylov_basis_alloc, par_vector_krylov_basis_free,            & 
            par_vector_krylov_basis_extract_view, par_vector_krylov_basis_multidot, & 
            par_vector_krylov_basis_multiaxpy
contains

  !=============================================================================
  subroutine par_vector_krylov_basis_alloc (k, p_v, Q)
    implicit none
    integer(ip)     , intent(in)               :: k
    type(par_vector), intent(in) , target      :: p_v
    type(par_vector_krylov_basis), intent(out) :: Q 

    ! p_part and p_part%p_context is required within this subroutine
    assert ( associated(p_v%p_part)          )
    assert ( associated(p_v%p_part%p_context) )
    assert ( p_v%p_part%p_context%created .eqv. .true.)
    Q%p_part => p_v%p_part

    if(p_v%p_part%p_context%iam<0) return

    call fem_vector_krylov_basis_alloc ( k, p_v%f_vector, Q%f_basis )

    ! Create epv_mown and epv_own_ext as views of f_vector 
    ! if ( Q%p_part%p_context%handler == trilinos ) then 
    !   ! epetra_multivector (i.e., set of scalar vectors)
    !   if ( Q%f_basis%storage == scal ) then                                                
    !     call epetra_multivector_construct_from_2D_array ( Q%epmv, Q%p_part%row_map(p_v%f_vector%nd), Q%f_basis%b, & 
    !                                          ! MyLDA                     ! NumVectors
    !                                       &  Q%f_basis%nd*Q%f_basis%neq, k )
    !   else ! (block vector) 
    !     write (0,*) 'Error: trilinos shadow interfaces do not yet support block epetra_vector objects'
    !     stop
    !   end if
    ! else 
      Q%state = p_v%state
    ! end if

  end subroutine par_vector_krylov_basis_alloc

  !=============================================================================
  subroutine par_vector_krylov_basis_free (Q)
    implicit none
    type(par_vector_krylov_basis), intent(inout) :: Q

    ! The routine requires the partition/context info
    assert ( associated( Q%p_part ) )
    assert ( associated( Q%p_part%p_context ) )
    assert ( Q%p_part%p_context%created .eqv. .true.)
    if(Q%p_part%p_context%iam<0) return

    ! if ( Q%p_part%p_context%handler == trilinos  ) then 
    !    ! epetra_vector (i.e., scalar vector)
    !    if (  Q%f_basis%storage == scal ) then 
    !       call epetra_multivector_destruct ( Q%epmv )
    !    else ! (block vector)
    !       write (0,*) 'Error: trilinos shadow interfaces do not yet support block epetra_vector objects'
    !       stop
    !    end if
    ! else
       Q%state = undefined
    ! end if

    call fem_vector_krylov_basis_free ( Q%f_basis )

  end subroutine par_vector_krylov_basis_free

  !=============================================================================
  subroutine par_vector_krylov_basis_extract_view (i, Q, p_v)
    implicit none
    integer(ip)     , intent(in)                      :: i
    type(par_vector_krylov_basis), intent(in), target :: Q
    type(par_vector), intent(out)                     :: p_v

    ! The routine requires the partition/context info
    assert ( associated( Q%p_part ) )
    assert ( associated( Q%p_part%p_context ) )
    assert ( Q%p_part%p_context%created .eqv. .true.)

    ! 2. fill p_v members
    p_v%p_part => Q%p_part

    if(Q%p_part%p_context%iam<0) return

    ! 1. fill p_v%f_vector members
    call fem_vector_krylov_basis_extract_view (i, Q%f_basis, p_v%f_vector )

    ! Create epv_own and epv_own_ext as views of f_vector 
    ! if ( p_v%p_part%p_context%handler == trilinos ) then 
    !    ! epetra_vector (i.e., scalar vector)
    !    if ( Q%f_basis%storage == scal ) then
    !       assert ( Q%p_part%maps_state(p_v%f_vector%nd) == map_created  )
    !       call epetra_vector_construct ( p_v%epv_own    , p_v%p_part%row_map(p_v%f_vector%nd), p_v%f_vector%b )
    !       call epetra_vector_construct ( p_v%epv_own_ext, p_v%p_part%col_map(p_v%f_vector%nd), p_v%f_vector%b )
    !    else ! (block vector) 
    !       write (0,*) 'Error: trilinos shadow interfaces do not yet support block epetra_vector objects'
    !       stop
    !    end if
    ! else if ( p_v%p_part%p_context%handler == inhouse ) then
       p_v%state  =  Q%state
    !end if
  end subroutine par_vector_krylov_basis_extract_view

  !=============================================================================
  ! s <- Q_k^T * p_v, with Q_k = (Q(1), Q(2), .. Q(k))
  subroutine par_vector_krylov_basis_multidot (k, Q, p_v, s)
     implicit none
     ! Parameters 
     integer(ip), intent(in)                   :: k
     type(par_vector_krylov_basis), intent(in) :: Q
     type(par_vector)             , intent(in) :: p_v
     real(rp), intent(out)                     :: s(k)

     ! Locals
     integer(ip)                               :: i
     !type(epetra_map)                          :: s_map
     !type(epetra_multivector)                  :: s_epmv
     !integer (kind=c_int), allocatable         :: s_myel(:)
     !integer (kind=c_int)                      :: ierrc
     !type(epetra_multivector)                  :: p_v_m
     type(par_vector)                          :: p_v_w
     !type(epetra_multivector)                  :: Q_k


     ! The routine requires the partition/context info
     assert ( associated( Q%p_part ) )
     assert ( associated( Q%p_part%p_context ) )
     assert ( Q%p_part%p_context%created .eqv. .true.)
     if(Q%p_part%p_context%iam<0) return

     ! if ( Q%p_part%p_context%handler == trilinos  ) then 
     !    ! epetra_multivector (i.e., set of scalars vectors)
     !    if ( Q%f_basis%storage == scal ) then 
     !      call memalloc(k, s_myel, __FILE__,__LINE__)
  
     !      ! Within Trilinos, s is a local replicated vector.
     !      ! s_map and s_epvm have to be created accordingly
     !      ! to this data distribution. 
     !      do i=1,k
     !         s_myel(i) = i-1
     !      end do
          
     !      assert ( Q%p_part%maps_state(p_v%f_vector%nd) == map_created  )

     !      ! Create epetra_multivector p_v_m
     !      call epetra_multivector_construct_from_2D_array ( p_v_m, p_v%p_part%row_map(p_v%f_vector%nd),  & 
     !                                                      & p_v%f_vector%b, p_v%f_vector%nd * p_v%f_vector%neq, 1 )

     !      ! Create s_map
     !      call epetra_map_construct ( s_map, k, k, s_myel, 0, Q%p_part%p_context%epcomm )
     !      call memfree (s_myel,__FILE__,__LINE__)

     !      ! Create s_epv
     !      call epetra_multivector_construct_from_2D_array ( s_epmv, s_map, s, k, 1 )

     !      call epetra_multivector_construct_from_epetra_multivector ( Q_k, Q%epmv, 0, k ) 

     !      ! Multiply (this = ScalarThis*\e this + ScalarAB*A*B.)
     !      call epetra_multivector_multiply  ( s_epmv, c_char_'T', c_char_'N', 1.0_c_double, Q_k, p_v_m, 0.0_c_double, ierrc )
     !      assert ( ierrc == 0 )

     !      call epetra_multivector_destruct ( Q_k )

     !      call epetra_multivector_destruct ( s_epmv )
     !      call epetra_map_destruct         ( s_map )

     !      ! Create epetra_multivector p_v_m
     !      call epetra_multivector_destruct ( p_v_m )

     !    else ! (block vector)
     !       write (0,*) 'Error: trilinos shadow interfaces do not yet support block epetra_vector objects'
     !       stop
     !    end if
     ! else
        assert ( p_v%state == Q%state )
 
        if ( p_v%state == full_summed ) then
           ! ******** p_v%f_vector should be weighted in advance !!!!
           call par_vector_clone  ( p_v, p_v_w )
           call par_vector_copy   ( p_v, p_v_w )
           call par_vector_weight ( p_v_w )
        else if ( p_v%state == part_summed ) then
           ! ******** p_v%f_vector should be comm. in advance !!!!
           call par_vector_clone  ( p_v, p_v_w )
           call par_vector_copy   ( p_v, p_v_w )
           call par_vector_comm   ( p_v_w )
        end if

        ! call par_vector_print ( 6, p_v )   ! DBG:
        ! call par_vector_print ( 6, p_v_w ) ! DBG:
        
        call fem_vector_krylov_basis_multidot ( k, Q%f_basis, p_v_w%f_vector, s )
        
        call par_vector_free  ( p_v_w )

        ! Reduce-sum local dot products on all processes
        call psb_sum ( Q%p_part%p_context%icontxt, s )
     ! end if
  end subroutine par_vector_krylov_basis_multidot

  !=============================================================================
  ! p_v <- p_v + alpha * Q_k * s
  subroutine par_vector_krylov_basis_multiaxpy (k, alpha, Q, s, p_v)
     implicit none
     ! Parameters
     integer(ip), intent(in)                      :: k
     real(rp)   , intent(in)                      :: alpha
     type(par_vector_krylov_basis), intent(in)    :: Q
     real(rp), intent(in)                         :: s(k)
     type(par_vector)             , intent(inout) :: p_v

     ! Locals
     integer(ip)                                  :: i
     !type(epetra_map)                             :: s_map
     !type(epetra_multivector)                     :: s_epmv
     !integer (kind=c_int), allocatable            :: s_myel(:)
     !integer (kind=c_int)                         :: ierrc
     !type(epetra_multivector)                     :: p_v_m
     !type(epetra_multivector)                     :: Q_k


     ! The routine requires the partition/context info
     assert ( associated( Q%p_part ) )
     assert ( associated( Q%p_part%p_context ) )
     assert ( Q%p_part%p_context%created .eqv. .true.)
     if(Q%p_part%p_context%iam<0) return


     ! if ( Q%p_part%p_context%handler == trilinos ) then 
     !    ! epetra_multivector (i.e., scalar vector)
     !    if ( Q%f_basis%storage == scal ) then
     !      call memalloc (k, s_myel, __FILE__,__LINE__)
  
     !      ! Within Trilinos, s is a local replicated vector.
     !      ! s_map and s_epv have to be created accordingly
     !      ! to this data distribution. 
     !      do i=1,k
     !         s_myel(i) = i-1
     !      end do
          
     !      assert ( Q%p_part%maps_state(p_v%f_vector%nd) == map_created  )

     !      ! Create epetra_multivector p_v_m
     !      call epetra_multivector_construct_from_2D_array ( p_v_m, p_v%p_part%row_map(p_v%f_vector%nd),  & 
     !                                                        p_v%f_vector%b, p_v%f_vector%nd * p_v%f_vector%neq, 1 )

     !      ! Create s_map
     !      call epetra_map_construct ( s_map, k, k, s_myel, 0, Q%p_part%p_context%epcomm )
     !      call memfree (s_myel,__FILE__,__LINE__)

     !      ! Create s_epv
     !      call epetra_multivector_construct_from_2D_array ( s_epmv, s_map, s, k, 1 )

     !      call epetra_multivector_construct_from_epetra_multivector ( Q_k, Q%epmv, 0, k ) 

     !      ! Multiply
     !      call epetra_multivector_multiply  ( p_v_m, c_char_'N', c_char_'N', alpha, Q_k, s_epmv, 1.0_c_double, ierrc )
     !      assert ( ierrc == 0 )

     !      call epetra_multivector_destruct ( Q_k )

     !      call epetra_multivector_destruct ( s_epmv )
     !      call epetra_map_destruct         ( s_map )

     !      ! Create epetra_multivector p_v_m
     !      call epetra_multivector_destruct ( p_v_m )

     !    else ! (block vector)
     !       write (0,*) 'Error: trilinos shadow interfaces do not yet support block epetra_vector objects'
     !       stop
     !    end if
     ! else
        assert ( p_v%state == Q%state )
        call fem_vector_krylov_basis_multiaxpy (k, alpha, Q%f_basis, s, p_v%f_vector)
     !end if

  end subroutine par_vector_krylov_basis_multiaxpy
 
end module par_vector_krylov_basis_names
