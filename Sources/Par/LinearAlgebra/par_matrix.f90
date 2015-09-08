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
module par_matrix_names
  ! Serial modules
  use types_names
  use memor_names
  use matrix_names
  use array_names
  use stdio_names
#ifdef memcheck
  use iso_c_binding
#endif

  ! Parallel modules
  use par_environment_names
  use par_context_names
  use par_graph_names
  use par_vector_names
  use psb_penv_mod_names
  use dof_distribution_names

  ! Abstract types
  use abstract_vector_names
  use abstract_operator_names
  
  implicit none
# include "debug.i90"

  private

  type, extends(abstract_operator_t) :: par_matrix_t
     ! Data structure which stores the local part 
     ! of the matrix mapped to the current processor.
     ! This is required for both eb and vb data 
     ! distributions
     type( matrix_t )       :: f_matrix

     type(par_graph_t), pointer :: &
        p_graph => NULL()             ! Associated par_graph
     
     type(dof_distribution_t), pointer :: &
        dof_dist => NULL()            ! Associated (ROW) dof_distribution
     
     type(dof_distribution_t), pointer :: &
        dof_dist_cols => NULL()       ! Associated (COL) dof_distribution

     type(par_environment_t), pointer :: &
          p_env => NULL()
   contains
     procedure  :: apply     => par_matrix_apply
     procedure  :: apply_fun => par_matrix_apply_fun
     procedure  :: free      => par_matrix_free_tbp
  end type par_matrix_t

  interface par_matrix_free
     module procedure par_matrix_free_one_shot, par_matrix_free_progressively
  end interface par_matrix_free

  ! Types
  public :: par_matrix_t

  ! Functions
  public :: par_matrix_create, par_matrix_graph, par_matrix_fill_val, &
         &  par_matrix_alloc, par_matrix_free,    & 
         &  par_matrix_print, par_matrix_print_matrix_market, &
         &  par_matrix_zero, par_matvec, par_matvec_trans

!***********************************************************************
! Allocatable arrays of type(par_matrix_t)
!***********************************************************************
# define var_attr allocatable, target
# define point(a,b) call move_alloc(a,b)
# define generic_status_test             allocated
# define generic_memalloc_interface      memalloc
# define generic_memrealloc_interface    memrealloc
# define generic_memfree_interface       memfree
# define generic_memmovealloc_interface  memmovealloc

# define var_type type(par_matrix_t)
# define var_size 80
# define bound_kind ip
# include "mem_header.i90"

  public :: memalloc,  memrealloc,  memfree, memmovealloc

contains

# include "mem_body.i90"

  !=============================================================================
  subroutine par_matrix_create(is_symmetric,dof_dist,dof_dist_cols,p_env,p_matrix,def)
    implicit none
    logical                             ,intent(in)  :: is_symmetric
    type(dof_distribution_t), target    ,intent(in)  :: dof_dist
    type(dof_distribution_t), target    ,intent(in)  :: dof_dist_cols
    type(par_environment_t) , target    ,intent(in)  :: p_env
    type(par_matrix_t)                  ,intent(out) :: p_matrix
    integer(ip)             , optional  ,intent(in)  :: def
    call matrix_create(is_symmetric,p_matrix%f_matrix,def)
    p_matrix%dof_dist      => dof_dist 
    p_matrix%dof_dist_cols => dof_dist_cols 
    p_matrix%p_env => p_env
  end subroutine par_matrix_create

  subroutine par_matrix_graph(p_graph,p_matrix)
    implicit none
    type(par_graph_t) , target, intent(in)   :: p_graph
    type(par_matrix_t), intent(inout)        :: p_matrix

    ! Pointer to part/context object is required
    assert ( associated(p_graph%dof_dist) )
    assert ( associated(p_graph%p_env%p_context) )
    assert ( p_graph%p_env%p_context%created .eqv. .true.)

    ! Point to input target parallel graph
    p_matrix%p_graph => p_graph

    ! if(p_graph%p_env%p_context%iam>=0) 
    call matrix_graph ( p_graph%f_graph, p_matrix%f_matrix )
  end subroutine par_matrix_graph

  subroutine par_matrix_fill_val(p_matrix)
    implicit none
    type(par_matrix_t), intent(inout)          :: p_matrix

    if(p_matrix%p_env%p_context%iam>=0) then
       call matrix_fill_val(p_matrix%f_matrix)
    else
       return
    end if

  end subroutine par_matrix_fill_val

  subroutine par_matrix_alloc(is_symmetric,p_graph,p_matrix,def)
    implicit none
    logical                     , intent(in)   :: is_symmetric
    type(par_graph_t) , target  , intent(in)   :: p_graph
    type(par_matrix_t)          , intent(out)  :: p_matrix
    integer(ip)       , optional, intent(in)   :: def
    call par_matrix_create(is_symmetric,p_graph%dof_dist,p_graph%dof_dist_cols,p_graph%p_env,p_matrix,def)
    call par_matrix_graph(p_graph,p_matrix)
    call par_matrix_fill_val(p_matrix)
  end subroutine par_matrix_alloc

  !=============================================================================
  subroutine par_matrix_free_one_shot(p_matrix)
    implicit none

    type(par_matrix_t), intent(inout) :: p_matrix
    call par_matrix_free_progressively(p_matrix, free_values)
    call par_matrix_free_progressively(p_matrix, free_struct)
    call par_matrix_free_progressively(p_matrix, free_clean)
  end subroutine par_matrix_free_one_shot

  !=============================================================================
  subroutine par_matrix_free_progressively(p_matrix, mode)
    implicit none
    type(par_matrix_t), intent(inout) :: p_matrix
    integer(ip)     , intent(in)    :: mode

    ! The routine requires the partition/context info
    assert ( associated(p_matrix%dof_dist) )
    assert ( associated(p_matrix%p_env%p_context) )
    assert ( p_matrix%p_env%p_context%created .eqv. .true.)
    assert ( mode == free_clean .or. mode == free_struct .or. mode == free_values )

    if(p_matrix%p_env%p_context%iam<0) return

    if ( mode == free_clean ) then
       nullify ( p_matrix%dof_dist )
       nullify ( p_matrix%dof_dist_cols )
       nullify ( p_matrix%p_env )
	   nullify ( p_matrix%p_graph )
    else if ( mode == free_struct ) then
       ! AFM: This nullification cannot be here as this means that it will not be longer possible
       !      to access p_matrix%dof_dist after "free_struct"ing a par_matrix
       !      (and it is done in many parts of the code). I will move it to free_clean.
       ! AFM: The comment above NO longer applies as par_matrix now directly points to
       !      the dof_distribution instance 
       ! nullify ( p_matrix%p_graph )
    end if

    ! Free local part
    call matrix_free ( p_matrix%f_matrix, mode )

  end subroutine par_matrix_free_progressively

  !=============================================================================
  subroutine par_matrix_print(lunou, p_matrix)
    implicit none
    type(par_matrix_t)  ,  intent(in) :: p_matrix
    integer(ip)      ,  intent(in) :: lunou

    ! p_graph%dof_dist is required within this subroutine
    assert ( associated(p_matrix%dof_dist) )
    
    ! p_graph%p_env%p_context is required within this subroutine
    assert ( associated(p_matrix%p_env%p_context) )

  end subroutine par_matrix_print

  subroutine par_matrix_print_matrix_market ( dir_path, prefix, p_mat )
    implicit none
    ! Parameters
    character (*)   , intent(in) :: dir_path
    character (*)   , intent(in) :: prefix
    type(par_matrix_t), intent(in) :: p_mat

    ! Locals
    integer         :: iam, num_procs, lunou
    integer(ip)     :: j, ndigs_iam, ndigs_num_procs, id_map
    character(256)  :: name 
    character(256)  :: zeros
    character(256)  :: part_id

    assert ( associated(p_mat%p_env%p_context) )
    assert ( p_mat%p_env%p_context%created .eqv. .true.)
    if(p_mat%p_env%p_context%iam<0) return


    name = trim(prefix) // '.par_matrix' // '.mtx'

    ! Get context info
    call par_context_info ( p_mat%p_env%p_context, iam, num_procs )

    ! Form the file_path of the partition object to be read
    iam = iam + 1 ! Partition identifers start from 1 !!
     
    ndigs_num_procs = count_digits_par_matrix (num_procs)
    zeros = ' '   
    ndigs_iam = count_digits_par_matrix ( iam )
   
    ! write(*,*) ndgs_num_procs, ndigs_iam DBG
    
    do j=1,  ndigs_num_procs - ndigs_iam
       zeros (j:j) = '0'
    end do
    part_id = ch(iam)

    ! Read partition data from path_file file
    lunou =  io_open (trim(dir_path) // '/' // trim(name) // '.' // trim(zeros) // trim(part_id), 'write')


    call matrix_print_matrix_market ( lunou, p_mat%f_matrix )

    call io_close (lunou)

  end subroutine par_matrix_print_matrix_market

  function count_digits_par_matrix ( i )
    implicit none
    ! Parameters
    integer(ip), intent(in) :: i 
    integer(ip)             :: count_digits_par_matrix
    ! Locals   
    integer(ip)             :: x 
    x = i 
    if (x < 0) x = -x;
    count_digits_par_matrix = 1;
    x = x/10;
    do while( x > 0)
       count_digits_par_matrix = count_digits_par_matrix + 1
       x = x/10;
    end do
  end function count_digits_par_matrix

  subroutine par_matrix_zero (p_matrix)
    implicit none
    ! Parameters 
    type(par_matrix_t), intent(inout)    :: p_matrix

    ! p_env%p_context is required within this subroutine
    assert ( associated(p_matrix%p_env%p_context) )
    assert ( p_matrix%p_env%p_context%created .eqv. .true.)

    if(p_matrix%p_env%p_context%iam<0) return

    call matrix_zero ( p_matrix%f_matrix )
  end subroutine par_matrix_zero

  subroutine par_matvec(a,x,y)
    implicit none
    ! Parameters
    type(par_matrix_t) , intent(in)    :: a
    type(par_vector_t) , intent(in)    :: x
    type(par_vector_t) , intent(inout) :: y
    ! Locals
    !integer(c_int)                   :: ierrc
    real :: aux

    ! This routine requires the partition/context info
    assert ( associated(a%p_env) )
    assert ( associated(a%p_env%p_context) )

    assert ( a%p_env%p_context%created .eqv. .true.)
    if(a%p_env%p_context%iam<0) return

    assert ( associated(x%p_env) )
    assert ( associated(x%p_env%p_context) )

    assert ( associated(y%p_env) )
    assert ( associated(y%p_env%p_context) )
    
    assert (x%state == full_summed) 
    ! write (*,*) 'MVAX'
    ! call vector_print ( 6, x%f_vector )
    ! write (*,*) 'MVAY'
    ! call vector_print ( 6, y%f_vector )
    call matrix_matvec (a%f_matrix, x%f_vector, y%f_vector) 
    ! write (*,*) 'MVD'
    ! call vector_print ( 6, y%f_vector )
    y%state = part_summed
  end subroutine par_matvec

  subroutine par_matvec_trans(a,x,y)
    implicit none
    ! Parameters
    type(par_matrix_t) , intent(in)    :: a
    type(par_vector_t) , intent(in)    :: x
    type(par_vector_t) , intent(inout) :: y
    ! Locals
    !integer(c_int)                   :: ierrc
	
    ! This routine requires the partition/context info
    assert ( associated(a%p_env) )
    assert ( associated(a%p_env%p_context) )
    if(a%p_env%p_context%iam<0) return
	
    assert ( associated(x%p_env) )
    assert ( associated(x%p_env%p_context) )

    assert ( associated(y%p_env) )
    assert ( associated(y%p_env%p_context) )
    
    assert (x%state == full_summed) 
    call matrix_matvec_trans (a%f_matrix, x%f_vector, y%f_vector) 
    y%state = part_summed
  end subroutine par_matvec_trans

  ! op%apply(x,y) <=> y <- op*x
  ! Implicitly assumes that y is already allocated
  subroutine par_matrix_apply(op,x,y) 
    implicit none
    class(par_matrix_t), intent(in)    :: op
    class(abstract_vector_t) , intent(in)    :: x
    class(abstract_vector_t) , intent(inout) :: y 

    call x%GuardTemp()

    select type(x)
    class is (par_vector_t)
       select type(y)
       class is(par_vector_t)
          call par_matvec(op, x, y)
          ! call vector_print(6,y)
       class default
          write(0,'(a)') 'par_matrix_t%apply: unsupported y class'
          check(1==0)
       end select
    class default
       write(0,'(a)') 'par_matrix_t%apply: unsupported x class'
       check(1==0)
    end select

    call x%CleanTemp()
  end subroutine par_matrix_apply

  ! op%apply(x)
  ! Allocates room for (temporary) y
  function par_matrix_apply_fun(op,x) result(y)
    implicit none
    class(par_matrix_t), intent(in)  :: op
    class(abstract_vector_t) , intent(in)  :: x
    class(abstract_vector_t) , allocatable :: y 

    type(par_vector_t), allocatable :: local_y

    select type(x)
    class is (par_vector_t)
       allocate(local_y)
       call par_vector_alloc ( op%dof_dist, x%p_env, local_y)
       call par_matvec(op, x, local_y)
       call move_alloc(local_y, y)
       call y%SetTemp()
    class default
       write(0,'(a)') 'par_matrix_t%apply_fun: unsupported x class'
       check(1==0)
    end select
  end function par_matrix_apply_fun

  subroutine par_matrix_free_tbp(this)
    implicit none
    class(par_matrix_t), intent(inout) :: this
  end subroutine par_matrix_free_tbp


end module par_matrix_names
